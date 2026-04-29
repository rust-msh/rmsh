//! Laplacian Smoothing — iterative nodal relocation mesh optimizer.
//!
//! # Algorithm overview
//!
//! Laplacian smoothing is the simplest and most widely used mesh quality
//! improvement technique.  Each interior node is moved to the (possibly
//! weighted) centroid of its immediate neighbours:
//!
//! ```text
//! x_i ← (1 - ω) · x_i  +  ω · (1/N · Σ_j x_j)
//! ```
//!
//! where the sum is over the `N` nodes adjacent to node `i`, and `ω ∈ (0, 1]`
//! is the relaxation factor.
//!
//! The standard (unweighted) Laplacian is fast and robust but may shrink the
//! mesh or produce poor results near convex/concave regions.  Improved variants
//! include:
//!
//! * **Taubin smoothing** (λ/μ method): alternating positive and negative
//!   relaxation steps to prevent shrinkage.
//! * **Weighted Laplacian**: weight neighbours by inverse edge length, area,
//!   or cotangent to improve triangle quality.
//! * **Constrained / boundary-locked** Laplacian: skip boundary nodes entirely
//!   (the default here when `params.move_boundary_nodes = false`).
//!
//! # Reference
//!
//! G. Taubin, "A Signal Processing Approach to Fair Surface Design", *SIGGRAPH
//! '95*, 1995 (Taubin smoothing).
//! Gmsh source: `Mesh/meshSmooth.cpp`.
//!
//! # Status
//!
//! **All three variants implemented** — `Uniform`, `Cotangent`-weighted, and
//! `Taubin` λ/μ scheme with alternating relaxation passes.

use rmsh_model::{Mesh, Node};

use crate::traits::{MeshAlgoError, MeshOptimizer, OptimizeParams};

// ─── Public struct ────────────────────────────────────────────────────────────

/// Laplacian smoothing optimizer.
///
/// Relocates each interior node toward the centroid of its neighbours to
/// reduce mesh distortion and improve element quality.
#[derive(Debug, Clone)]
pub struct LaplacianSmooth {
    /// Smoothing variant to apply.
    pub variant: LaplacianVariant,

    /// Relaxation factor `ω ∈ (0, 1]`.
    ///
    /// Larger values converge faster but may be less stable.  Defaults to `0.5`.
    pub omega: f64,
}

impl Default for LaplacianSmooth {
    fn default() -> Self {
        Self {
            variant: LaplacianVariant::Uniform,
            omega: 0.5,
        }
    }
}

impl LaplacianSmooth {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_variant(mut self, variant: LaplacianVariant) -> Self {
        self.variant = variant;
        self
    }
}

// ─── Variant enum ─────────────────────────────────────────────────────────────

/// The smoothing weighting scheme.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum LaplacianVariant {
    /// Unweighted centroid average (fastest, may shrink mesh).
    #[default]
    Uniform,

    /// Cotangent-weighted (preserves mesh area better for 2-D triangle meshes).
    Cotangent,

    /// Taubin λ/μ scheme: alternating positive (`λ`) and negative (`μ`) steps
    /// to prevent volume shrinkage.
    ///
    /// Each iteration does a λ-step (positive relaxation) followed by a μ-step
    /// (negative relaxation).  With `|μ| > λ` the scheme preserves volume while
    /// smoothing high-frequency noise.
    Taubin {
        /// Positive step size (`λ > 0`).  Typical: `0.5`.
        lambda: f64,
        /// Negative step size (`μ < 0`).  Typical: `-0.53`.
        mu: f64,
    },
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl MeshOptimizer for LaplacianSmooth {
    fn name(&self) -> &'static str {
        "Laplacian Smooth"
    }

    fn optimize(&self, mesh: &mut Mesh, params: &OptimizeParams) -> Result<(), MeshAlgoError> {
        let boundary_node_ids = collect_boundary_nodes(mesh);

        match &self.variant {
            LaplacianVariant::Taubin { lambda, mu } => {
                // Taubin: alternating positive (λ) and negative (μ) passes.
                // The two-pass pair is one "iteration".
                for _iter in 0..params.iterations {
                    let d1 = laplacian_pass_uniform(
                        mesh,
                        &boundary_node_ids,
                        *lambda,
                        params.move_boundary_nodes,
                    )?;
                    let d2 = laplacian_pass_uniform(
                        mesh,
                        &boundary_node_ids,
                        *mu,
                        params.move_boundary_nodes,
                    )?;
                    if d1.max(d2) < params.tolerance {
                        break;
                    }
                }
            }
            _ => {
                for _iter in 0..params.iterations {
                    let max_displacement = laplacian_pass(
                        mesh,
                        &boundary_node_ids,
                        self.omega,
                        &self.variant,
                        params.move_boundary_nodes,
                    )?;
                    if max_displacement < params.tolerance {
                        break;
                    }
                }
            }
        }

        Ok(())
    }
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

/// Perform one Laplacian smoothing pass over all (eligible) nodes.
///
/// Returns the maximum node displacement in this pass (used for convergence
/// checking).
fn laplacian_pass(
    mesh: &mut Mesh,
    boundary_ids: &std::collections::HashSet<u64>,
    omega: f64,
    variant: &LaplacianVariant,
    move_boundary: bool,
) -> Result<f64, MeshAlgoError> {
    // Build adjacency: node_id → list of neighbour node IDs
    let adjacency = build_node_adjacency(mesh);

    let mut max_disp = 0.0_f64;
    let node_ids: Vec<u64> = mesh.nodes.keys().copied().collect();

    for id in node_ids {
        if !move_boundary && boundary_ids.contains(&id) {
            continue;
        }
        let neighbors = match adjacency.get(&id) {
            Some(nbrs) if !nbrs.is_empty() => nbrs,
            _ => continue,
        };

        let new_pos = match variant {
            LaplacianVariant::Uniform => uniform_centroid(id, neighbors, mesh)?,
            LaplacianVariant::Cotangent => {
                cotangent_centroid(id, neighbors, mesh)?
            }
            LaplacianVariant::Taubin { .. } => {
                // Taubin uses a uniform pass with the lambda/mu coefficient.
                // Each individual step is just uniform Laplacian with a specific ω.
                return laplacian_pass_uniform(mesh, boundary_ids, omega, move_boundary);
            }
        };

        // Apply relaxation
        let node = mesh.nodes.get_mut(&id).unwrap();
        let old = [node.position.x, node.position.y, node.position.z];
        let relaxed = [
            (1.0 - omega) * old[0] + omega * new_pos[0],
            (1.0 - omega) * old[1] + omega * new_pos[1],
            (1.0 - omega) * old[2] + omega * new_pos[2],
        ];
        let disp = {
            let dx = relaxed[0] - old[0];
            let dy = relaxed[1] - old[1];
            let dz = relaxed[2] - old[2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        };
        node.position.x = relaxed[0];
        node.position.y = relaxed[1];
        node.position.z = relaxed[2];
        max_disp = max_disp.max(disp);
    }

    Ok(max_disp)
}

/// Compute the unweighted centroid of a node's neighbours.
fn uniform_centroid(_id: u64, neighbors: &[u64], mesh: &Mesh) -> Result<[f64; 3], MeshAlgoError> {
    let n = neighbors.len() as f64;
    let mut sum = [0.0_f64; 3];
    for &nb_id in neighbors {
        let nb = mesh
            .nodes
            .get(&nb_id)
            .ok_or_else(|| MeshAlgoError::Generation(format!("missing neighbour node {nb_id}")))?;
        sum[0] += nb.position.x;
        sum[1] += nb.position.y;
        sum[2] += nb.position.z;
    }
    Ok([sum[0] / n, sum[1] / n, sum[2] / n])
}

/// Perform one Laplacian pass with a uniform (unweighted) centroid for all nodes.
///
/// Used by the Taubin variant where alternating positive/negative relaxation
/// factors are applied using the same centroid computation.
fn laplacian_pass_uniform(
    mesh: &mut Mesh,
    boundary_ids: &std::collections::HashSet<u64>,
    omega: f64,
    move_boundary: bool,
) -> Result<f64, MeshAlgoError> {
    let adjacency = build_node_adjacency(mesh);
    let mut max_disp = 0.0_f64;
    let node_ids: Vec<u64> = mesh.nodes.keys().copied().collect();

    for id in node_ids {
        if !move_boundary && boundary_ids.contains(&id) {
            continue;
        }
        let neighbors = match adjacency.get(&id) {
            Some(nbrs) if !nbrs.is_empty() => nbrs,
            _ => continue,
        };

        let centroid = uniform_centroid(id, neighbors, mesh)?;
        let node = mesh.nodes.get_mut(&id).unwrap();
        let old = [node.position.x, node.position.y, node.position.z];
        let relaxed = [
            (1.0 - omega) * old[0] + omega * centroid[0],
            (1.0 - omega) * old[1] + omega * centroid[1],
            (1.0 - omega) * old[2] + omega * centroid[2],
        ];
        node.position.x = relaxed[0];
        node.position.y = relaxed[1];
        node.position.z = relaxed[2];

        let dx = relaxed[0] - old[0];
        let dy = relaxed[1] - old[1];
        let dz = relaxed[2] - old[2];
        max_disp = max_disp.max((dx * dx + dy * dy + dz * dz).sqrt());
    }

    Ok(max_disp)
}

/// Compute the cotangent of the angle at vertex `opp` in triangle `(u, v, opp)`.
///
/// Returns `cot(θ) = dot(u-opp, v-opp) / |cross(u-opp, v-opp)|`.
/// Returns 0 when the triangle is degenerate (zero area).
fn cotangent(pu: [f64; 3], pv: [f64; 3], popp: [f64; 3]) -> f64 {
    let a = [pu[0] - popp[0], pu[1] - popp[1], pu[2] - popp[2]];
    let b = [pv[0] - popp[0], pv[1] - popp[1], pv[2] - popp[2]];
    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let cross_norm = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
    if cross_norm < 1e-15 {
        return 0.0;
    }
    dot / cross_norm
}

/// Build a map from each edge to its one or two opposite vertices in triangle
/// elements of the mesh.
///
/// For edge `(a, b)` the map stores `[opp1, opp2]`, where `opp1`/`opp2` are
/// the third vertices of the one or two triangles sharing that edge.  Boundary
/// edges have `opp2 == u64::MAX`.
fn build_edge_triangle_map(mesh: &Mesh) -> std::collections::HashMap<(u64, u64), [u64; 2]> {
    use std::collections::HashMap;
    let mut edge_to_opp: HashMap<(u64, u64), Vec<u64>> = HashMap::new();
    for elem in &mesh.elements {
        if elem.node_ids.len() == 3 {
            let ids = &elem.node_ids;
            let edges = [(ids[0], ids[1], ids[2]),
                         (ids[1], ids[2], ids[0]),
                         (ids[0], ids[2], ids[1])];
            for (u, v, opp) in edges {
                let key = if u < v { (u, v) } else { (v, u) };
                edge_to_opp.entry(key).or_default().push(opp);
            }
        }
    }
    edge_to_opp
        .into_iter()
        .map(|(k, v)| {
            let opp = if v.len() >= 2 {
                [v[0], v[1]]
            } else if v.len() == 1 {
                [v[0], u64::MAX]
            } else {
                [u64::MAX, u64::MAX]
            };
            (k, opp)
        })
        .collect()
}

/// Compute the cotangent-weighted centroid of a node's neighbours.
///
/// Weight `w_ij = (cot(α_ij) + cot(β_ij)) / 2` where `α_ij` and `β_ij` are
/// the angles opposite edge `(i, j)` in the two adjacent triangles.
/// Falls back to uniform centroid when no valid cotangent weights exist.
fn cotangent_centroid(id: u64, neighbors: &[u64], mesh: &Mesh) -> Result<[f64; 3], MeshAlgoError> {
    let edge_map = build_edge_triangle_map(mesh);
    let pi = {
        let n = mesh.nodes.get(&id)
            .ok_or_else(|| MeshAlgoError::Generation(format!("missing node {id}")))?;
        [n.position.x, n.position.y, n.position.z]
    };

    let mut weighted_sum = [0.0_f64; 3];
    let mut total_weight = 0.0_f64;

    for &j in neighbors {
        let key = if id < j { (id, j) } else { (j, id) };
        let opp = edge_map.get(&key).copied().unwrap_or([u64::MAX, u64::MAX]);

        let pj = {
            let n = mesh.nodes.get(&j)
                .ok_or_else(|| MeshAlgoError::Generation(format!("missing neighbour node {j}")))?;
            [n.position.x, n.position.y, n.position.z]
        };

        let mut w = 0.0_f64;
        let mut valid_count = 0u32;
        for &o in &opp {
            if o != u64::MAX {
                let po = {
                    let n = mesh.nodes.get(&o)
                        .ok_or_else(|| MeshAlgoError::Generation(format!("missing opposite node {o}")))?;
                    [n.position.x, n.position.y, n.position.z]
                };
                let c = cotangent(pi, pj, po);
                // Only accumulate positive cotangents (negative ones indicate
                // degenerate or inverted elements).
                if c > 1e-12 {
                    w += c;
                    valid_count += 1;
                }
            }
        }
        if valid_count > 0 {
            w /= valid_count as f64;
            weighted_sum[0] += w * pj[0];
            weighted_sum[1] += w * pj[1];
            weighted_sum[2] += w * pj[2];
            total_weight += w;
        }
    }

    if total_weight <= 1e-15 {
        // Fall back to uniform if cotangent weights are degenerate.
        return uniform_centroid(id, neighbors, mesh);
    }

    let inv = 1.0 / total_weight;
    Ok([weighted_sum[0] * inv, weighted_sum[1] * inv, weighted_sum[2] * inv])
}

/// Build a node-to-neighbours adjacency map from the element connectivity.
fn build_node_adjacency(mesh: &Mesh) -> std::collections::HashMap<u64, Vec<u64>> {
    use std::collections::{HashMap, HashSet};
    let mut adj: HashMap<u64, HashSet<u64>> = HashMap::new();
    for elem in &mesh.elements {
        let nodes = &elem.node_ids;
        for (i, &a) in nodes.iter().enumerate() {
            for &b in nodes.iter().skip(i + 1) {
                adj.entry(a).or_default().insert(b);
                adj.entry(b).or_default().insert(a);
            }
        }
    }
    adj.into_iter()
        .map(|(k, v)| (k, v.into_iter().collect()))
        .collect()
}

/// Collect the IDs of all nodes that lie on the mesh boundary.
///
/// A boundary node is one that belongs to at least one boundary face (a face
/// not shared by two elements).
fn collect_boundary_nodes(mesh: &Mesh) -> std::collections::HashSet<u64> {
    use std::collections::HashMap;

    // Count how many elements share each face
    let mut face_count: HashMap<[u64; 3], u32> = HashMap::new();
    for elem in &mesh.elements {
        for face in element_faces(&elem.node_ids) {
            *face_count.entry(face).or_insert(0) += 1;
        }
    }

    let mut boundary_nodes = std::collections::HashSet::new();
    for (face, count) in face_count {
        if count == 1 {
            for id in face {
                boundary_nodes.insert(id);
            }
        }
    }
    boundary_nodes
}

/// Return the triangular faces of an element given its node list.
///
/// For Tet4: 4 faces; for Tri3: 1 face (itself).  Other element types
/// are approximated by the first three nodes.
fn element_faces(nodes: &[u64]) -> Vec<[u64; 3]> {
    match nodes.len() {
        4 => {
            // Tet4: four faces
            let [a, b, c, d] = [nodes[0], nodes[1], nodes[2], nodes[3]];
            let mut faces = vec![[a, b, c], [a, b, d], [a, c, d], [b, c, d]];
            for f in &mut faces {
                f.sort_unstable();
            }
            faces
        }
        3 => {
            let mut f = [nodes[0], nodes[1], nodes[2]];
            f.sort_unstable();
            vec![f]
        }
        _ => {
            // Fallback: emit all combinations of 3 nodes from the first 4
            if nodes.len() >= 3 {
                let mut f = [nodes[0], nodes[1], nodes[2]];
                f.sort_unstable();
                vec![f]
            } else {
                vec![]
            }
        }
    }
}

// Bring Node into scope for the position field access — verify it compiles.
const _: fn() = || {
    let _: fn(&Node) -> [f64; 3] = |n| [n.position.x, n.position.y, n.position.z];
};

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rmsh_model::{Element, ElementType, Mesh, Node};

    /// A small flat triangle mesh: 4 nodes, 2 triangles.
    ///
    /// ```text
    /// 3---4
    /// |\ 2|
    /// |1\ |
    /// 1---2
    /// ```
    fn two_triangle_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 1.0, 1.0, 0.0));
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));
        mesh.add_element(Element::new(2, ElementType::Triangle3, vec![2, 4, 3]));
        mesh
    }

    /// A 3×3 uniform triangle mesh of the unit square.
    /// Interior nodes can move; corner/edge nodes are boundary-locked.
    fn uniform_grid_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        // 9 nodes in a 3x3 grid
        let mut nid = 1u64;
        for j in 0..3 {
            for i in 0..3 {
                mesh.add_node(Node::new(nid, i as f64 * 0.5, j as f64 * 0.5, 0.0));
                nid += 1;
            }
        }
        // node numbering: row-major, id = j*3 + i + 1
        // row 0: 1,2,3   row 1: 4,5,6   row 2: 7,8,9
        let tris = [
            [1u64, 2, 4],
            [2, 5, 4],
            [2, 3, 5],
            [3, 6, 5],
            [4, 5, 7],
            [5, 8, 7],
            [5, 6, 8],
            [6, 9, 8],
        ];
        for (i, tri) in tris.iter().enumerate() {
            mesh.add_element(Element::new(
                i as u64 + 1,
                ElementType::Triangle3,
                tri.to_vec(),
            ));
        }
        mesh
    }

    // ── Basic behavior ────────────────────────────────────────────────────────

    #[test]
    fn smooth_does_not_change_node_count() {
        let mut mesh = two_triangle_mesh();
        let n_before = mesh.node_count();
        LaplacianSmooth::new()
            .optimize(&mut mesh, &OptimizeParams::default())
            .expect("smooth should succeed");
        assert_eq!(mesh.node_count(), n_before);
    }

    #[test]
    fn smooth_does_not_change_element_count() {
        let mut mesh = two_triangle_mesh();
        let e_before = mesh.element_count();
        LaplacianSmooth::new()
            .optimize(&mut mesh, &OptimizeParams::default())
            .expect("smooth should succeed");
        assert_eq!(mesh.element_count(), e_before);
    }

    #[test]
    fn smooth_does_not_move_boundary_nodes_by_default() {
        let mut mesh = two_triangle_mesh();
        // Snapshot all 4 node positions before.
        let before: std::collections::HashMap<u64, [f64; 3]> = mesh
            .nodes
            .iter()
            .map(|(&id, n)| (id, [n.position.x, n.position.y, n.position.z]))
            .collect();

        LaplacianSmooth::new()
            .optimize(&mut mesh, &OptimizeParams::default())
            .expect("smooth should succeed");

        let boundary = collect_boundary_nodes(&mesh);
        for id in &boundary {
            let old = before[id];
            let new_node = &mesh.nodes[id];
            let new = [new_node.position.x, new_node.position.y, new_node.position.z];
            assert_eq!(old, new, "boundary node {id} must not move");
        }
    }

    #[test]
    fn smooth_returns_ok_on_empty_mesh() {
        let mut mesh = Mesh::new();
        let result = LaplacianSmooth::new().optimize(&mut mesh, &OptimizeParams::default());
        assert!(result.is_ok());
    }

    // ── Uniform variant ───────────────────────────────────────────────────────

    #[test]
    fn uniform_smooth_name_is_stable() {
        assert_eq!(LaplacianSmooth::new().name(), "Laplacian Smooth");
    }

    #[test]
    fn uniform_smooth_on_regular_grid_converges() {
        // A regular grid is already at its Laplacian equilibrium.
        // After smoothing (boundary locked), interior node positions should
        // change at most a small amount.
        let mut mesh = uniform_grid_mesh();
        // Snapshot the interior node (id=5, center of grid).
        let before = {
            let n = &mesh.nodes[&5];
            [n.position.x, n.position.y, n.position.z]
        };
        LaplacianSmooth { omega: 1.0, variant: LaplacianVariant::Uniform }
            .optimize(
                &mut mesh,
                &OptimizeParams {
                    iterations: 20,
                    tolerance: 1e-12,
                    move_boundary_nodes: false,
                },
            )
            .expect("smooth should succeed");
        let after = {
            let n = &mesh.nodes[&5];
            [n.position.x, n.position.y, n.position.z]
        };
        let disp = {
            let dx = after[0] - before[0];
            let dy = after[1] - before[1];
            let dz = after[2] - before[2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        };
        // Center node of a symmetric uniform grid stays near its initial position.
        assert!(disp < 0.1, "interior node moved too much: {disp}");
    }

    #[test]
    fn uniform_smooth_with_omega_zero_point_five_keeps_nodes_inside_domain() {
        let mut mesh = two_triangle_mesh();
        LaplacianSmooth { omega: 0.5, variant: LaplacianVariant::Uniform }
            .optimize(&mut mesh, &OptimizeParams { iterations: 10, ..Default::default() })
            .expect("smooth should succeed");
        // All nodes must stay in or near [0,1]x[0,1] after smoothing.
        for node in mesh.nodes.values() {
            assert!(
                node.position.x >= -0.1 && node.position.x <= 1.1,
                "x out of range: {}",
                node.position.x
            );
            assert!(
                node.position.y >= -0.1 && node.position.y <= 1.1,
                "y out of range: {}",
                node.position.y
            );
        }
    }

    // ── Unimplemented variants ────────────────────────────────────────────────

    #[test]
    fn cotangent_variant_smooths_triangle_mesh() {
        let mut mesh = uniform_grid_mesh();
        LaplacianSmooth::new()
            .with_variant(LaplacianVariant::Cotangent)
            .optimize(
                &mut mesh,
                &OptimizeParams {
                    iterations: 10,
                    move_boundary_nodes: true,
                    ..Default::default()
                },
            )
            .expect("Cotangent should succeed");
        // Node/element counts unchanged.
        assert_eq!(mesh.node_count(), 9);
        assert_eq!(mesh.elements_by_dimension(2).len(), 8);
        // Interior node (id=5) should have moved from (0.5, 0.5) due to
        // non-uniform cotangent weights on a regular grid.
        let interior = &mesh.nodes[&5];
        let moved = (interior.position.x - 0.5).abs() > 1e-6
            || (interior.position.y - 0.5).abs() > 1e-6;
        assert!(moved, "interior node should have moved after cotangent smooth");
    }

    #[test]
    fn taubin_variant_smooths_without_error() {
        let mut mesh = two_triangle_mesh();
        let result = LaplacianSmooth::new()
            .with_variant(LaplacianVariant::Taubin { lambda: 0.5, mu: -0.53 })
            .optimize(
                &mut mesh,
                &OptimizeParams {
                    move_boundary_nodes: true,
                    iterations: 5,
                    ..Default::default()
                },
            );
        assert!(result.is_ok(), "Taubin should succeed, got: {result:?}");
        // Mesh should still have the same number of nodes and elements.
        assert_eq!(mesh.node_count(), 4);
        assert_eq!(mesh.element_count(), 2);
    }

    // ── Internal helpers ──────────────────────────────────────────────────────

    #[test]
    fn build_node_adjacency_two_triangles() {
        let mesh = two_triangle_mesh();
        let adj = build_node_adjacency(&mesh);
        // Node 3 should be adjacent to 1, 2, and 4 (shared between both tris).
        let mut n3 = adj[&3].clone();
        n3.sort();
        assert_eq!(n3, vec![1, 2, 4]);
    }

    #[test]
    fn collect_boundary_nodes_two_triangles() {
        let mesh = two_triangle_mesh();
        // In this mesh ALL faces are boundary (each triangle has 3 distinct faces,
        // and the shared face [2,3] appears in both triangles so it is internal).
        // Boundary nodes = nodes on boundary faces.
        let boundary = collect_boundary_nodes(&mesh);
        // The shared edge between the two tris is [2, 3] (sorted).
        // All 4 nodes lie on at least one boundary face.
        assert_eq!(boundary.len(), 4);
    }

    #[test]
    fn element_faces_tet4_returns_four_sorted_faces() {
        let faces = element_faces(&[10, 20, 30, 40]);
        assert_eq!(faces.len(), 4);
        // Each face must have its nodes sorted.
        for f in &faces {
            let mut s = *f;
            s.sort_unstable();
            assert_eq!(*f, s);
        }
        // Check specific faces present (sorted [10,20,30], [10,20,40], [10,30,40], [20,30,40]).
        let expected: std::collections::HashSet<[u64; 3]> = [
            [10, 20, 30],
            [10, 20, 40],
            [10, 30, 40],
            [20, 30, 40],
        ]
        .into();
        let got: std::collections::HashSet<[u64; 3]> = faces.into_iter().collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn element_faces_tri3_returns_one_sorted_face() {
        let faces = element_faces(&[5, 3, 7]);
        assert_eq!(faces.len(), 1);
        assert_eq!(faces[0], [3, 5, 7]);
    }

    #[test]
    fn element_faces_degenerate_fallback() {
        // 1 node: empty
        assert!(element_faces(&[1]).is_empty());
        // 5 nodes: fallback to first 3
        let faces = element_faces(&[9, 7, 5, 3, 1]);
        assert_eq!(faces.len(), 1);
    }
}
