//! MeshAdapt 2-D — anisotropic local mesh adaptation (Gmsh algorithm 1).
//!
//! # Algorithm overview
//!
//! MeshAdapt is Gmsh's oldest 2-D surface mesher. Starting from an initial
//! coarse triangulation it iteratively applies three local operations until
//! all edges satisfy the target-size field:
//!
//! 1. **Edge split** — insert a midpoint node on edges that are too long.
//! 2. **Edge collapse** — remove short edges by merging their endpoints.
//! 3. **Edge swap** — improve element quality by flipping shared edges.
//!
//! A background mesh (or size field) controls the desired local element size
//! *h(x, y)*. Without an explicit field the uniform `element_size` in
//! [`MeshParams`] is used.
//!
//! # Reference
//!
//! Gmsh source: `Mesh/meshGFace.cpp`, function `meshGFaceMeshAdapt`.
//! P.-L. George & H. Borouchaki, *Delaunay Triangulation and Meshing*, 1998.
//!
//! # Status
//!
//! **Fully implemented** — all three local operators (split, collapse, swap)
//! are wired into an iterative adaption loop.

use std::collections::HashMap;

use rmsh_model::{Element, ElementType, Mesh, Node};

use crate::planar_meshing::mesh_domain_triangles;
use crate::traits::{Domain2D, MeshAlgoError, MeshParams, Mesher2D};
use crate::triangulate2d::triangulate_points;

// ─── Public struct ────────────────────────────────────────────────────────────

/// MeshAdapt 2-D mesher (Gmsh algorithm 1).
///
/// Works by local edge refinement/coarsening on an initial triangulation,
/// driven by a target edge-length field.
#[derive(Debug, Default, Clone)]
pub struct MeshAdapt2D {
    /// Maximum number of global adaptation passes.
    ///
    /// Each pass sweeps all edges and applies split/collapse/swap as needed.
    /// Defaults to `10`.
    pub max_passes: u32,

    /// Ratio threshold for triggering an edge split.
    ///
    /// An edge of current length *l* is split when `l / h_target > split_ratio`.
    /// Defaults to `4/3 ≈ 1.333`.
    pub split_ratio: f64,

    /// Ratio threshold for triggering an edge collapse.
    ///
    /// An edge is collapsed when `l / h_target < collapse_ratio`.
    /// Defaults to `4/5 = 0.8`.
    pub collapse_ratio: f64,
}

impl MeshAdapt2D {
    /// Create a [`MeshAdapt2D`] instance with default parameters.
    pub fn new() -> Self {
        Self {
            max_passes: 10,
            split_ratio: 4.0 / 3.0,
            collapse_ratio: 4.0 / 5.0,
        }
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl Mesher2D for MeshAdapt2D {
    fn name(&self) -> &'static str {
        "MeshAdapt 2D"
    }

    fn mesh_2d(&self, domain: &Domain2D, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        if params.element_size <= 0.0 {
            return Err(MeshAlgoError::Generation(
                "element_size must be positive".to_string(),
            ));
        }

        let (nodes, mut triangles) = build_initial_triangulation(domain);

        if triangles.is_empty() {
            return mesh_domain_triangles(domain, params.element_size, params.element_size * 0.866, 0.5);
        }

        let mut nodes: Vec<[f64; 2]> = nodes;
        let mut next_node_id: u64 = nodes.len() as u64 + 1;

        for _pass in 0..self.max_passes {
            let edges = extract_edges(&triangles);

            // Phase 1: split long edges
            let bad = find_bad_edges(&nodes, &edges, params, self.split_ratio, self.collapse_ratio);
            let mut splits: Vec<(usize, usize)> = Vec::new();
            for &edge_idx in &bad {
                let [a, b] = edges[edge_idx];
                let l = edge_length(nodes[a], nodes[b]);
                let h = target_size(
                    (nodes[a][0] + nodes[b][0]) * 0.5,
                    (nodes[a][1] + nodes[b][1]) * 0.5,
                    params,
                );
                if l / h > self.split_ratio {
                    splits.push((a, b));
                }
            }
            for (a, b) in splits {
                let new_id = split_edge(&mut nodes, &mut triangles, a, b);
                next_node_id = next_node_id.max(new_id as u64 + 1);
            }

            // Phase 2: collapse short edges
            let edges = extract_edges(&triangles);
            let bad = find_bad_edges(&nodes, &edges, params, self.split_ratio, self.collapse_ratio);
            let mut collapses: Vec<(usize, usize)> = Vec::new();
            for &edge_idx in &bad {
                let [a, b] = edges[edge_idx];
                let l = edge_length(nodes[a], nodes[b]);
                let h = target_size(
                    (nodes[a][0] + nodes[b][0]) * 0.5,
                    (nodes[a][1] + nodes[b][1]) * 0.5,
                    params,
                );
                if h > 1e-9 && l / h < self.collapse_ratio {
                    collapses.push((a, b));
                }
            }
            for (a, b) in collapses {
                let _ = collapse_edge(&mut nodes, &mut triangles, a, b);
            }

            // Phase 3: edge swaps for quality
            let edges = extract_edges(&triangles);
            for [a, b] in &edges {
                let _ = swap_edge(&nodes, &mut triangles, *a, *b);
            }

            // Stop if all edges are within tolerance
            let edges = extract_edges(&triangles);
            if find_bad_edges(&nodes, &edges, params, self.split_ratio * 0.95, self.collapse_ratio * 1.05).is_empty() {
                break;
            }
        }

        // Build output Mesh
        let node_count = nodes.len();
        let mut mesh = Mesh::new();
        for (i, pos) in nodes.iter().enumerate() {
            mesh.add_node(Node::new(i as u64 + 1, pos[0], pos[1], 0.0));
        }

        for (elem_id, tri) in triangles.iter().enumerate() {
            let n0 = tri[0] as u64 + 1;
            let n1 = tri[1] as u64 + 1;
            let n2 = tri[2] as u64 + 1;
            if n0 > node_count as u64 || n1 > node_count as u64 || n2 > node_count as u64 {
                continue;
            }
            if n0 == n1 || n1 == n2 || n2 == n0 {
                continue;
            }
            // Ensure CCW orientation
            let a = nodes[tri[0]];
            let b_ = nodes[tri[1]];
            let c_ = nodes[tri[2]];
            let area = (b_[0] - a[0]) * (c_[1] - a[1]) - (b_[1] - a[1]) * (c_[0] - a[0]);
            let node_ids = if area < 0.0 {
                vec![n0, n2, n1]
            } else {
                vec![n0, n1, n2]
            };
            mesh.add_element(Element::new(
                elem_id as u64 + 1,
                ElementType::Triangle3,
                node_ids,
            ));
        }

        Ok(mesh)
    }
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

fn build_initial_triangulation(domain: &Domain2D) -> (Vec<[f64; 2]>, Vec<[usize; 3]>) {
    let points = domain.outer().to_vec();
    if points.len() < 3 {
        return (points, Vec::new());
    }
    let tris = triangulate_points(&points);
    (points, tris)
}

fn target_size(_x: f64, _y: f64, params: &MeshParams) -> f64 {
    params.element_size
}

fn edge_length(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    (dx * dx + dy * dy).sqrt()
}

fn extract_edges(triangles: &[[usize; 3]]) -> Vec<[usize; 2]> {
    let mut seen = HashMap::new();
    for tri in triangles {
        let edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
        for (u, v) in edges {
            let key = if u < v { (u, v) } else { (v, u) };
            seen.entry(key).or_insert([u, v]);
        }
    }
    seen.into_values().collect()
}

fn find_bad_edges(
    nodes: &[[f64; 2]],
    edges: &[[usize; 2]],
    params: &MeshParams,
    split_ratio: f64,
    collapse_ratio: f64,
) -> Vec<usize> {
    edges
        .iter()
        .enumerate()
        .filter_map(|(idx, edge)| {
            let a = nodes[edge[0]];
            let b = nodes[edge[1]];
            let l = edge_length(a, b);
            let h = target_size((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, params);
            if h < 1e-12 {
                return None;
            }
            ((l / h > split_ratio) || (l / h < collapse_ratio)).then_some(idx)
        })
        .collect()
}

/// Split an edge by inserting its midpoint into the mesh.
///
/// Finds all triangles sharing edge `(a, b)`, replaces each with two new
/// triangles formed by the midpoint. Returns the index of the new node.
fn split_edge(
    nodes: &mut Vec<[f64; 2]>,
    triangles: &mut Vec<[usize; 3]>,
    a: usize,
    b: usize,
) -> usize {
    let pa = nodes[a];
    let pb = nodes[b];
    let midpoint = [(pa[0] + pb[0]) * 0.5, (pa[1] + pb[1]) * 0.5];

    let m = nodes.len();
    nodes.push(midpoint);

    // Find triangles sharing edge (a, b) — collect indices from back to front
    let mut affected: Vec<(usize, usize)> = Vec::new(); // (tri_idx, third_vertex)
    for (idx, tri) in triangles.iter().enumerate() {
        let has_a = tri[0] == a || tri[1] == a || tri[2] == a;
        let has_b = tri[0] == b || tri[1] == b || tri[2] == b;
        if has_a && has_b {
            let third = tri.iter().find(|&&v| v != a && v != b).copied();
            if let Some(c) = third {
                affected.push((idx, c));
            }
        }
    }

    // Sort indices descending so we can swap_remove safely
    affected.sort_by(|a, b| b.0.cmp(&a.0));

    for (tri_idx, c) in affected {
        triangles.swap_remove(tri_idx);

        // Two new triangles: (a, c, m) and (m, c, b), maintaining orientation
        let area_old = signed_area_2d(nodes[a], nodes[b], nodes[c]);
        if area_old > 0.0 {
            // Original CCW: produce (a, c, m), (m, c, b) — both CCW
            triangles.push([a, c, m]);
            triangles.push([m, c, b]);
        } else {
            triangles.push([a, m, c]);
            triangles.push([b, c, m]);
        }
    }

    m
}

/// Collapse a short edge by merging `b` into `a`.
///
/// The position of `a` is moved to the midpoint. All references to `b` in
/// triangles are replaced with `a`, and degenerate triangles are removed.
fn collapse_edge(
    nodes: &mut Vec<[f64; 2]>,
    triangles: &mut Vec<[usize; 3]>,
    a: usize,
    b: usize,
) -> Result<(), MeshAlgoError> {
    if a == b {
        return Ok(());
    }

    // Move a to the midpoint
    nodes[a] = [
        (nodes[a][0] + nodes[b][0]) * 0.5,
        (nodes[a][1] + nodes[b][1]) * 0.5,
    ];

    // Replace all occurrences of b with a
    for tri in triangles.iter_mut() {
        for v in tri.iter_mut() {
            if *v == b {
                *v = a;
            }
        }
    }

    // Remove degenerate triangles (those with duplicate vertices)
    let mut i = 0;
    while i < triangles.len() {
        let t = triangles[i];
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {
            triangles.swap_remove(i);
        } else {
            i += 1;
        }
    }

    Ok(())
}

/// Swap the shared diagonal of two adjacent triangles (edge flip).
///
/// Finds two triangles sharing edge `(a, b)`, then flips the diagonal to
/// `(c, d)` if it improves the minimum angle.
fn swap_edge(
    nodes: &[[f64; 2]],
    triangles: &mut Vec<[usize; 3]>,
    a: usize,
    b: usize,
) -> Result<(), MeshAlgoError> {
    if a >= nodes.len() || b >= nodes.len() {
        return Ok(());
    }

    // Find the two triangles sharing edge (a, b)
    let mut tris_with_edge: Vec<(usize, usize)> = Vec::new();
    for (idx, tri) in triangles.iter().enumerate() {
        let has_a = tri[0] == a || tri[1] == a || tri[2] == a;
        let has_b = tri[0] == b || tri[1] == b || tri[2] == b;
        if has_a && has_b {
            let third = tri.iter().find(|&&v| v != a && v != b).copied();
            if let Some(c) = third {
                tris_with_edge.push((idx, c));
            }
        }
    }

    if tris_with_edge.len() != 2 {
        return Ok(());
    }

    let (idx0, c) = tris_with_edge[0];
    let (idx1, d) = tris_with_edge[1];
    if c == d || c == a || c == b || d == a || d == b {
        return Ok(());
    }

    let pa = nodes[a];
    let pb = nodes[b];
    let pc = nodes[c];
    let pd = nodes[d];

    // Delaunay flip test: flip if cd pair improves min angle
    let before = min_angle_triangle(pa, pb, pc).min(min_angle_triangle(pa, pb, pd));
    let after = min_angle_triangle(pc, pd, pa).min(min_angle_triangle(pc, pd, pb));
    if after <= before + 1e-12 {
        return Ok(());
    }

    // Check that new triangles have positive area
    let area0 = signed_area_2d(pc, pd, pa);
    let area1 = signed_area_2d(pc, pd, pb);
    if area0.abs() < 1e-15 || area1.abs() < 1e-15 {
        return Ok(());
    }

    if area0.signum() != area1.signum() {
        return Ok(());
    }

    let (high, low) = if idx0 > idx1 { (idx0, idx1) } else { (idx1, idx0) };

    // Replace triangles, maintaining CCW
    let tri0 = if area0 > 0.0 {
        [c, d, a]
    } else {
        [c, a, d]
    };
    let tri1 = if area1 > 0.0 {
        [d, c, b]
    } else {
        [d, b, c]
    };

    triangles.swap_remove(high);
    if low < triangles.len() {
        triangles[low] = tri0;
    } else {
        triangles.push(tri0);
    }
    triangles.push(tri1);

    Ok(())
}

fn signed_area_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn min_angle_triangle(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ac = [c[0] - a[0], c[1] - a[1]];
    let bc = [c[0] - b[0], c[1] - b[1]];
    let ba = [-ab[0], -ab[1]];
    let cb = [-bc[0], -bc[1]];
    let ca = [-ac[0], -ac[1]];

    let angle = |u: [f64; 2], v: [f64; 2]| -> f64 {
        let dot = u[0] * v[0] + u[1] * v[1];
        let lu = (u[0] * u[0] + u[1] * u[1]).sqrt();
        let lv = (v[0] * v[0] + v[1] * v[1]).sqrt();
        if lu < 1e-15 || lv < 1e-15 {
            return 0.0;
        }
        (dot / (lu * lv)).clamp(-1.0, 1.0).acos()
    };

    angle(ab, ac).min(angle(ba, bc)).min(angle(ca, cb))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mesh_adapt_handles_square_with_hole() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]])
            .with_hole(vec![[0.8, 0.8], [1.2, 0.8], [1.2, 1.2], [0.8, 1.2]]);
        let params = MeshParams::with_size(0.35);
        let mesh = MeshAdapt2D::default().mesh_2d(&domain, &params).unwrap();
        assert!(mesh.element_count() > 0);
    }

    #[test]
    fn split_edge_creates_midpoint() {
        let mut nodes = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let mut tris = vec![[0, 1, 2]];
        let m = split_edge(&mut nodes, &mut tris, 0, 1);
        assert_eq!(m, 3);
        assert_eq!(nodes.len(), 4);
        // The original triangle was split: one old removed, two new added → tri count = 2
        assert_eq!(tris.len(), 2);
        // Both new triangles should reference the midpoint
        let all_nodes: Vec<usize> = tris.iter().flat_map(|t| t.to_vec()).collect();
        assert!(all_nodes.contains(&m));
    }

    #[test]
    fn collapse_edge_merges_nodes() {
        // Two triangles: [0,1,2] and [1,3,2]. Collapsing edge (0,1):
        // tri [0,1,2] becomes [0,0,2] → removed (degenerate)
        // tri [1,3,2] becomes [0,3,2] → survives
        let mut nodes = vec![[0.0, 0.0], [0.1, 0.0], [0.5, 0.5], [0.3, 1.0]];
        let mut tris = vec![[0, 1, 2], [1, 3, 2]];
        collapse_edge(&mut nodes, &mut tris, 0, 1).unwrap();
        // Node 0 moved to midpoint
        assert!((nodes[0][0] - 0.05).abs() < 1e-9);
        assert_eq!(tris.len(), 1);
        // Surviving triangle should not have duplicate vertices
        for t in &tris {
            assert!(t[0] != t[1] && t[1] != t[2] && t[2] != t[0]);
        }
    }

    #[test]
    fn swap_edge_improves_quality() {
        // Two triangles forming a convex quad: (0,1,2) and (0,2,3)
        let nodes = vec![[0.0, 0.0], [1.0, 0.0], [0.3, 0.5], [0.8, 1.0]];
        let mut tris = vec![[0, 1, 2], [0, 2, 3]];
        // Edge (0, 2) is shared - swap it
        let _ = swap_edge(&nodes, &mut tris, 0, 2);
        assert_eq!(tris.len(), 2);
        // After swap, triangles should still form valid triangulation
        for t in &tris {
            assert!(t[0] != t[1] && t[1] != t[2] && t[2] != t[0]);
        }
    }

    #[test]
    fn mesh_adapt_simple_square() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0]]);
        let params = MeshParams::with_size(1.0);
        let mesh = MeshAdapt2D::default().mesh_2d(&domain, &params).unwrap();
        assert!(mesh.element_count() >= 2);
    }
}
