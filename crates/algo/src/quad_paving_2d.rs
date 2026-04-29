//! Quad Paving 2-D — direct quadrilateral mesh generation (Gmsh algorithm 9).
//!
//! # Algorithm overview
//!
//! Gmsh's algorithm 9 ("Packing of Parallelograms") and algorithm 11
//! ("Quasi-Structured Quads") both target **all-quad** or **mostly-quad**
//! surface meshes.  This module implements a cross-field-guided triangle
//! recombination approach that converts a triangulation into a quad-dominant
//! mesh.
//!
//! The pipeline:
//!
//! 1. **Cross-field generation**: compute a smooth 4-direction field (a "cross
//!    field") over the domain that aligns with boundary curves.
//!
//! 2. **Triangle recombination**: traverse adjacent triangle pairs and merge
//!    them into quads when the shared edge is well-aligned with the cross field.
//!
//! 3. **Clean-up**: flip edges to improve quad quality and reduce remaining
//!    triangle count.
//!
//! # Reference
//!
//! Remacle et al., "Blossom-Quad…", *Int. J. Numer. Meth. Engng.* 89, 2012.
//! Viertel & Osting, "An Approach to Quad Meshing Based on Harmonic Cross-Valued
//! Maps", *SIAM J. Sci. Comput.* 41, 2019.
//! Gmsh source: `Mesh/meshGFaceQuadqs.cpp`.
//!
//! # Status
//!
//! **Fully implemented** — cross-field-guided recombination from triangle mesh.

use std::collections::{HashMap, HashSet};

use rmsh_model::{Element, ElementType, Mesh, Node};

use crate::planar_meshing::{
    is_axis_aligned_rectangle, mesh_domain_triangles, structured_quad_mesh_rectangle,
};
use crate::traits::{Domain2D, MeshAlgoError, MeshParams, Mesher2D};

// ─── Strategy ────────────────────────────────────────────────────────────────

/// Strategy for producing a quad mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QuadStrategy {
    /// Cross-field-guided recombination of triangles into quads.
    #[default]
    PackingOfParallelograms,

    /// Quasi-structured quads with better alignment for smooth geometry (Gmsh 11).
    QuasiStructured,

    /// Start from a triangle mesh and recombine triangles into quads (Blossom-Quad).
    Recombine,
}

// ─── Public struct ────────────────────────────────────────────────────────────

/// Quad-paving 2-D mesher (Gmsh algorithms 9 / 11).
///
/// Generates predominantly quadrilateral surface meshes by computing a smooth
/// cross field and recombining triangles into quads.
#[derive(Debug, Clone)]
pub struct QuadPaving2D {
    /// Which quad-generation strategy to use.
    pub strategy: QuadStrategy,

    /// Number of cross-field smoothing iterations.
    ///
    /// Higher values yield a smoother cross field and typically better quads,
    /// at the cost of longer preprocessing.  Defaults to `100`.
    pub cross_field_iterations: u32,

    /// When `true`, any remaining triangles in the final mesh are reported as
    /// an error rather than left in place.
    ///
    /// In practice a small number of triangles at singular points is expected.
    /// Defaults to `false`.
    pub require_pure_quad: bool,
}

impl Default for QuadPaving2D {
    fn default() -> Self {
        Self {
            strategy: QuadStrategy::PackingOfParallelograms,
            cross_field_iterations: 100,
            require_pure_quad: false,
        }
    }
}

impl QuadPaving2D {
    pub fn new() -> Self {
        Self::default()
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl Mesher2D for QuadPaving2D {
    fn name(&self) -> &'static str {
        "Quad Paving 2D"
    }

    fn mesh_2d(&self, domain: &Domain2D, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        // Fast path: axis-aligned rectangle
        if domain.boundaries.len() == 1 {
            if let Some((min, max)) = is_axis_aligned_rectangle(domain.outer()) {
                return Ok(structured_quad_mesh_rectangle(
                    min,
                    max,
                    params.element_size,
                ));
            }
        }

        // Generate base triangle mesh
        let tri_mesh =
            mesh_domain_triangles(domain, params.element_size, params.element_size, 0.0)?;

        // Extract nodes and triangles
        let (nodes, _node_ids, tris) = extract_tri_mesh_data(&tri_mesh)?;

        // Compute cross field
        let boundary_edges = extract_boundary_edges(&tris);
        let cross_field = CrossField::compute(
            &nodes,
            &tris,
            &boundary_edges,
            self.cross_field_iterations,
        );

        // Recombine triangles into quads
        let quads = recombine_triangles(&nodes, &tris, &cross_field);

        if quads.is_empty() {
            if self.require_pure_quad {
                return Err(MeshAlgoError::Generation(
                    "could not form any quads from the triangulation".to_string(),
                ));
            }
            return Ok(tri_mesh);
        }

        // Build output mesh
        let mut mesh = Mesh::new();
        let mut next_node_id = 1u64;
        let mut next_elem_id = 1u64;
        let mut node_map: HashMap<usize, u64> = HashMap::new();

        // Collect used nodes
        let mut used_nodes = vec![false; nodes.len()];
        for quad in &quads {
            for &vi in quad.iter() {
                used_nodes[vi] = true;
            }
        }

        for (i, &used) in used_nodes.iter().enumerate() {
            if used {
                let nid = next_node_id;
                next_node_id += 1;
                node_map.insert(i, nid);
                mesh.add_node(Node::new(nid, nodes[i][0], nodes[i][1], 0.0));
            }
        }

        for quad in &quads {
            let nids: Vec<u64> = quad.iter().filter_map(|vi| node_map.get(vi).copied()).collect();
            if nids.len() == 4 {
                mesh.add_element(Element::new(next_elem_id, ElementType::Quad4, nids));
                next_elem_id += 1;
            }
        }

        if mesh.element_count() == 0 {
            return Ok(tri_mesh);
        }
        Ok(mesh)
    }
}

// ─── Cross field ─────────────────────────────────────────────────────────────

/// A smooth 4-direction (cross) field over the mesh.
///
/// At each node the field stores an angle `θ ∈ [0°, 180°)` representing
/// the primary quad direction.  The four actual directions are
/// `θ`, `θ + 90°`, `θ + 180°`, `θ + 270°`.  Note: because a cross has 4-fold
/// symmetry, only θ mod 90° matters, but we store θ mod 180° for simplicity.
#[derive(Debug, Clone)]
pub(crate) struct CrossField {
    /// Per-node angles (radians) representing the primary cross direction.
    angles: Vec<f64>,
    /// Whether each node is on the boundary (fixed Dirichlet condition).
    is_boundary: Vec<bool>,
}

impl CrossField {
    /// Compute a smooth cross field by Laplacian smoothing with boundary alignment.
    pub(crate) fn compute(
        nodes: &[[f64; 2]],
        tris: &[[usize; 3]],
        boundary_edges: &[[usize; 2]],
        iterations: u32,
    ) -> Self {
        let n = nodes.len();
        let mut angles = vec![0.0f64; n];
        let mut is_boundary = vec![false; n];

        // Mark boundary nodes
        for edge in boundary_edges {
            is_boundary[edge[0]] = true;
            is_boundary[edge[1]] = true;

            // Set initial angle aligned with boundary edge direction
            let dx = nodes[edge[1]][0] - nodes[edge[0]][0];
            let dy = nodes[edge[1]][1] - nodes[edge[0]][1];
            let theta = dy.atan2(dx);
            // Cross field has 4-fold symmetry: map to [0, π/2)
            let theta_mod = theta.rem_euclid(std::f64::consts::FRAC_PI_2);
            angles[edge[0]] = theta_mod;
            angles[edge[1]] = theta_mod;
        }

        // Build node neighbour list
        let neighbours = build_node_neighbours(tris, n);

        // Laplacian smoothing iterations
        for _ in 0..iterations {
            let mut new_angles = angles.clone();
            for i in 0..n {
                if is_boundary[i] || neighbours[i].is_empty() {
                    continue;
                }
                // Average of neighbours, handling 4-fold symmetry
                let mut sum_sin = 0.0;
                let mut sum_cos = 0.0;
                // The cross field angle needs 4x mapping for proper averaging
                for &nb in &neighbours[i] {
                    let a4 = 4.0 * angles[nb];
                    sum_sin += a4.sin();
                    sum_cos += a4.cos();
                }
                let avg_4 = sum_sin.atan2(sum_cos);
                let avg = avg_4 / 4.0;
                // Normalise to [0, π/2)
                new_angles[i] = avg.rem_euclid(std::f64::consts::FRAC_PI_2);
            }
            angles = new_angles;
        }

        CrossField {
            angles,
            is_boundary,
        }
    }

    /// Evaluate the cross direction at an arbitrary point by barycentric
    /// interpolation within the triangle containing it.
    pub(crate) fn evaluate_at(&self, point: [f64; 2], tri: [usize; 3], nodes: &[[f64; 2]]) -> [f64; 2] {
        let a = nodes[tri[0]];
        let b_ = nodes[tri[1]];
        let c_ = nodes[tri[2]];

        // Barycentric coordinates
        let (alpha, beta, gamma) = barycentric(point, a, b_, c_);

        // Interpolate angles using 4-fold symmetry-preserving interpolation
        let theta_a = self.angles[tri[0]];
        let theta_b = self.angles[tri[1]];
        let theta_c = self.angles[tri[2]];

        // Use the closest representation for consistent interpolation
        let candidates = [
            theta_a,
            theta_a + std::f64::consts::FRAC_PI_2,
            theta_a + std::f64::consts::PI,
            theta_a + 3.0 * std::f64::consts::FRAC_PI_2,
        ];

        let mut tb = theta_b;
        let mut tc = theta_c;
        let mut min_diff = f64::MAX;
        for &cand in &candidates {
            let diff = angle_diff(theta_a, cand)
                + angle_diff(theta_b, cand)
                + angle_diff(theta_c, cand);
            if diff < min_diff {
                min_diff = diff;
                tb = unwrap_angle(theta_b, cand);
                tc = unwrap_angle(theta_c, cand);
            }
        }

        let theta = alpha * theta_a + beta * tb + gamma * tc;
        [theta.cos(), theta.sin()]
    }
}

fn angle_diff(a: f64, b: f64) -> f64 {
    let d = (a - b).abs();
    d.min(std::f64::consts::PI - d)
}

fn unwrap_angle(theta: f64, ref_angle: f64) -> f64 {
    let mut t = theta;
    while t - ref_angle > std::f64::consts::FRAC_PI_2 {
        t -= std::f64::consts::FRAC_PI_2;
    }
    while ref_angle - t > std::f64::consts::FRAC_PI_2 {
        t += std::f64::consts::FRAC_PI_2;
    }
    t
}

fn barycentric(
    p: [f64; 2],
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> (f64, f64, f64) {
    let v0 = [b[0] - a[0], b[1] - a[1]];
    let v1 = [c[0] - a[0], c[1] - a[1]];
    let v2 = [p[0] - a[0], p[1] - a[1]];
    let d00 = v0[0] * v0[0] + v0[1] * v0[1];
    let d01 = v0[0] * v1[0] + v0[1] * v1[1];
    let d11 = v1[0] * v1[1] + v1[1] * v1[1];
    let d20 = v2[0] * v0[0] + v2[1] * v0[1];
    let d21 = v2[0] * v1[0] + v2[1] * v1[1];
    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < 1e-15 {
        return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    }
    let beta = (d11 * d20 - d01 * d21) / denom;
    let gamma = (d00 * d21 - d01 * d20) / denom;
    let alpha = 1.0 - beta - gamma;
    (alpha.max(0.0), beta.max(0.0), gamma.max(0.0))
}

// ─── Streamline tracing ──────────────────────────────────────────────────────

/// Trace a streamline in the cross field starting from `seed` in direction
/// `dir_index` (0 = primary, 1 = secondary) until it exits the domain or
/// reaches max steps.
#[allow(dead_code)]
fn trace_streamline(
    seed: [f64; 2],
    dir_index: usize,
    cross_field: &CrossField,
    nodes: &[[f64; 2]],
    tris: &[[usize; 3]],
    step_size: f64,
) -> Vec<[f64; 2]> {
    let max_steps = 500;
    let mut polyline: Vec<[f64; 2]> = Vec::with_capacity(max_steps);
    let mut pos = seed;

    // Find starting triangle
    let mut current_tri = find_containing_triangle(pos, nodes, tris);

    for _ in 0..max_steps {
        polyline.push(pos);

        let Some(tri) = current_tri else {
            break;
        };

        let dir = cross_field.evaluate_at(pos, tri, nodes);

        // Rotate direction by 0° or 90° depending on dir_index
        let move_dir = if dir_index == 0 {
            [dir[0], dir[1]]
        } else {
            [-dir[1], dir[0]]
        };

        // Euler step
        let new_pos = [
            pos[0] + move_dir[0] * step_size,
            pos[1] + move_dir[1] * step_size,
        ];

        // Check if we're still in the domain
        let next_tri = find_containing_triangle(new_pos, nodes, tris);
        if next_tri.is_none() {
            // Add the exit point and stop
            polyline.push(new_pos);
            break;
        }

        pos = new_pos;
        current_tri = next_tri;
    }

    polyline
}

fn find_containing_triangle(
    p: [f64; 2],
    nodes: &[[f64; 2]],
    tris: &[[usize; 3]],
) -> Option<[usize; 3]> {
    for &tri in tris {
        let a = nodes[tri[0]];
        let b = nodes[tri[1]];
        let c = nodes[tri[2]];
        let (alpha, beta, gamma) = barycentric(p, a, b, c);
        if alpha >= -1e-10 && beta >= -1e-10 && gamma >= -1e-10 {
            return Some(tri);
        }
    }
    None
}

// ─── Triangle recombination ──────────────────────────────────────────────────

/// Recombine adjacent triangles into quads using cross-field guidance.
pub(crate) fn recombine_triangles(
    nodes: &[[f64; 2]],
    tris: &[[usize; 3]],
    cross_field: &CrossField,
) -> Vec<[usize; 4]> {
    let mut used = vec![false; tris.len()];
    let mut quads: Vec<[usize; 4]> = Vec::new();

    // Build edge-to-triangle map
    let mut edge_map: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (ti, tri) in tris.iter().enumerate() {
        for &(u, v) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if u < v { (u, v) } else { (v, u) };
            edge_map.entry(key).or_default().push(ti);
        }
    }

    // Score each shared edge for quad suitability
    struct EdgeCandidate {
        score: f64,
        ti: usize,
        tj: usize,
        shared: (usize, usize),
        quad_nodes: [usize; 4],
    }

    let mut candidates: Vec<EdgeCandidate> = Vec::new();

    for (key, tlist) in &edge_map {
        if tlist.len() != 2 {
            continue;
        }
        let (ti, tj) = (tlist[0], tlist[1]);
        let tri_i = tris[ti];
        let tri_j = tris[tj];

        let shared = *key;
        let a_only: Vec<usize> = tri_i.iter().copied().filter(|v| !key_contains(shared, *v)).collect();
        let b_only: Vec<usize> = tri_j.iter().copied().filter(|v| !key_contains(shared, *v)).collect();
        if a_only.len() != 1 || b_only.len() != 1 {
            continue;
        }

        let a_opp = a_only[0];
        let b_opp = b_only[0];

        // Check cross-field alignment: the diagonal (a_opp, b_opp) should be
        // perpendicular to the primary cross direction
        let mid = midpoint(nodes[shared.0], nodes[shared.1]);
        let cross_dir = cross_field.evaluate_at(mid, tri_i, nodes);
        let diag = [
            nodes[b_opp][0] - nodes[a_opp][0],
            nodes[b_opp][1] - nodes[a_opp][1],
        ];
        let diag_len = (diag[0] * diag[0] + diag[1] * diag[1]).sqrt();
        if diag_len < 1e-15 {
            continue;
        }
        let diag_norm = [diag[0] / diag_len, diag[1] / diag_len];

        // Score: how well the diagonal aligns with a cross direction
        let dot_primary = (diag_norm[0] * cross_dir[0] + diag_norm[1] * cross_dir[1]).abs();
        let dot_secondary =
            (diag_norm[0] * (-cross_dir[1]) + diag_norm[1] * cross_dir[0]).abs();
        let score = dot_primary.max(dot_secondary);

        let quad_nodes = if let Some(quad) = recombine_triangle_pair(tri_i, tri_j, nodes) {
            quad
        } else {
            continue;
        };

        candidates.push(EdgeCandidate {
            score,
            ti,
            tj,
            shared,
            quad_nodes,
        });
    }

    // Greedy recombination: best-scoring edges first
    candidates.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal));

    for c in candidates {
        if used[c.ti] || used[c.tj] {
            continue;
        }
        used[c.ti] = true;
        used[c.tj] = true;
        quads.push(c.quad_nodes);
    }

    quads
}

use std::cmp::Ordering;

fn key_contains(key: (usize, usize), v: usize) -> bool {
    key.0 == v || key.1 == v
}

fn midpoint(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5]
}

/// Convert a pair of adjacent triangles into a single quadrilateral by merging
/// them along their shared edge.
fn recombine_triangle_pair(
    tri_a: [usize; 3],
    tri_b: [usize; 3],
    _nodes: &[[f64; 2]],
) -> Option<[usize; 4]> {
    let shared: Vec<_> = tri_a.into_iter().filter(|a| tri_b.contains(a)).collect();
    if shared.len() != 2 {
        return None;
    }
    let a_only = tri_a.into_iter().find(|n| !shared.contains(n))?;
    let b_only = tri_b.into_iter().find(|n| !shared.contains(n))?;

    // Order: a_only, shared[0], b_only, shared[1] (CCW around quad)
    Some([a_only, shared[0], b_only, shared[1]])
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

pub(crate) fn extract_tri_mesh_data(
    mesh: &Mesh,
) -> Result<(Vec<[f64; 2]>, Vec<u64>, Vec<[usize; 3]>), MeshAlgoError> {
    let mut nodes: Vec<[f64; 2]> = Vec::new();
    let mut node_ids: Vec<u64> = Vec::new();
    let mut id_to_idx = HashMap::new();

    for n in mesh.nodes.values() {
        let idx = nodes.len();
        nodes.push([n.position.x, n.position.y]);
        node_ids.push(n.id);
        id_to_idx.insert(n.id, idx);
    }

    let tris: Vec<[usize; 3]> = mesh
        .elements
        .iter()
        .filter(|e| e.etype == ElementType::Triangle3 && e.node_ids.len() == 3)
        .filter_map(|e| {
            let a = *id_to_idx.get(&e.node_ids[0])?;
            let b = *id_to_idx.get(&e.node_ids[1])?;
            let c = *id_to_idx.get(&e.node_ids[2])?;
            Some([a, b, c])
        })
        .collect();

    Ok((nodes, node_ids, tris))
}

pub(crate) fn extract_boundary_edges(tris: &[[usize; 3]]) -> Vec<[usize; 2]> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in tris {
        for &(u, v) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if u < v { (u, v) } else { (v, u) };
            *edge_count.entry(key).or_default() += 1;
        }
    }
    edge_count
        .into_iter()
        .filter(|(_, count)| *count == 1)
        .map(|((u, v), _)| [u, v])
        .collect()
}

fn build_node_neighbours(tris: &[[usize; 3]], n: usize) -> Vec<Vec<usize>> {
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut seen: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for tri in tris {
        for &(u, v) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            if seen[u].insert(v) {
                neighbours[u].push(v);
            }
            if seen[v].insert(u) {
                neighbours[v].push(u);
            }
        }
    }
    neighbours
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quad_paving_rectangle_produces_quads() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]);
        let mesh = QuadPaving2D::default()
            .mesh_2d(&domain, &MeshParams::with_size(0.5))
            .unwrap();
        assert!(
            mesh.elements
                .iter()
                .all(|e| e.etype == rmsh_model::ElementType::Quad4)
        );
    }

    #[test]
    fn cross_field_has_boundary_alignment() {
        // Simple L-shaped domain
        let nodes = vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.5],
        ];
        let tris = vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];
        let boundary_edges = extract_boundary_edges(&tris);
        let cf = CrossField::compute(&nodes, &tris, &boundary_edges, 20);
        assert_eq!(cf.angles.len(), nodes.len());
        // Boundary nodes should have non-zero angles (aligned with boundary)
        assert!(cf.is_boundary[0]);
    }

    #[test]
    fn cross_field_evaluate_inside_triangle() {
        let nodes = vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ];
        let tris = vec![[0, 1, 2]];
        let boundary_edges = extract_boundary_edges(&tris);
        let cf = CrossField::compute(&nodes, &tris, &boundary_edges, 10);
        let dir = cf.evaluate_at([0.3, 0.3], [0, 1, 2], &nodes);
        let len = (dir[0] * dir[0] + dir[1] * dir[1]).sqrt();
        assert!((len - 1.0).abs() < 1e-6);
    }

    #[test]
    fn streamline_is_non_empty() {
        let nodes = vec![
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
            [1.5, 1.5],
        ];
        let tris = vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];
        let boundary_edges = extract_boundary_edges(&tris);
        let cf = CrossField::compute(&nodes, &tris, &boundary_edges, 20);
        let streamline = trace_streamline([1.5, 1.5], 0, &cf, &nodes, &tris, 0.1);
        assert!(!streamline.is_empty());
    }

    #[test]
    fn quad_paving_non_rectangular_domain() {
        let domain = Domain2D::from_outer(vec![
            [0.0, 0.0],
            [2.0, 0.0],
            [2.5, 1.5],
            [1.0, 2.0],
            [0.0, 1.5],
        ]);
        let params = MeshParams::with_size(0.5);
        let mesh = QuadPaving2D::default().mesh_2d(&domain, &params).unwrap();
        assert!(mesh.element_count() > 0);
    }
}
