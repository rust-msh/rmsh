//! Frontal-Quads 2-D — frontal Delaunay mesh with quad recombination (Gmsh algorithm 8).
//!
//! # Algorithm overview
//!
//! Frontal-Quads generates a quadrilateral mesh in three stages:
//!
//! 1. **Frontal Delaunay triangulation**: generates a high-quality triangular mesh
//!    via [`crate::FrontalDelaunay2D`] (Gmsh algorithm 6).
//! 2. **Cross-field recombination**: computes a smooth 4-direction cross field on
//!    the triangle mesh, then greedily recombines adjacent triangle pairs into
//!    quadrilaterals following the cross-field alignment.
//! 3. **Local quality improvement**: iteratively swaps diagonals between adjacent
//!    quads to maximise the minimum interior angle.
//!
//! The result is a quadrilateral-dominant mesh, though some isolated triangles
//! may remain in regions where recombination is not possible.
//!
//! # Reference
//!
//! Gmsh source: `Mesh/meshGFace.cpp`, algorithm 8 path.

use std::collections::{HashMap, HashSet};

use rmsh_model::{Element, ElementType, Mesh, Node};

use crate::frontal_delaunay_2d::FrontalDelaunay2D;
use crate::quad_paving_2d::{
    extract_boundary_edges, extract_tri_mesh_data, recombine_triangles, CrossField,
};
use crate::traits::{Domain2D, MeshAlgoError, MeshParams, Mesher2D};

/// Frontal-Quads 2-D mesher (Gmsh algorithm 8).
///
/// Generates quadrilateral meshes by computing a cross field on a
/// Frontal-Delaunay triangle mesh and recombining adjacent triangles.
#[derive(Debug, Clone)]
pub struct FrontalQuads2D {
    /// Number of cross-field smoothing iterations.
    pub cross_field_iterations: u32,
    /// When true, run the local quality-improvement pass after recombination.
    pub enable_improvement: bool,
}

impl Default for FrontalQuads2D {
    fn default() -> Self {
        Self {
            cross_field_iterations: 100,
            enable_improvement: true,
        }
    }
}

impl FrontalQuads2D {
    pub fn new() -> Self {
        Self::default()
    }
}

impl Mesher2D for FrontalQuads2D {
    fn name(&self) -> &'static str {
        "Frontal-Quads 2D"
    }

    fn mesh_2d(&self, domain: &Domain2D, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        // Stage 1: generate triangle mesh via Frontal Delaunay
        let tri_mesh = FrontalDelaunay2D::default().mesh_2d(domain, params)?;

        // Stage 2: extract flat triangle data
        let (nodes, node_ids, tris) = extract_tri_mesh_data(&tri_mesh)?;
        if tris.len() < 2 {
            return Ok(tri_mesh);
        }

        // Stage 3: compute cross field
        let boundary_edges = extract_boundary_edges(&tris);
        let cf = CrossField::compute(&nodes, &tris, &boundary_edges, self.cross_field_iterations);

        // Stage 4: recombine triangles into quadrilaterals
        let mut quads = recombine_triangles(&nodes, &tris, &cf);

        // Stage 5: local quality improvement (diagonal swaps)
        if self.enable_improvement && quads.len() > 1 {
            improve_quad_quality(&mut quads, &nodes, &tris);
        }

        // Stage 6: build output Mesh
        build_quad_mesh(&nodes, &node_ids, &tris, &quads)
    }
}

// ─── Quad quality metrics ─────────────────────────────────────────────────────

/// Minimum interior angle of a quadrilateral (degrees).
/// Uses the 4 nodes from the quad definition.
fn quad_min_angle(quad: &[usize; 4], nodes: &[[f64; 2]]) -> f64 {
    let n = [
        nodes[quad[0]], nodes[quad[1]],
        nodes[quad[2]], nodes[quad[3]],
    ];
    let mut min_ang = f64::MAX;
    for i in 0..4 {
        let prev = n[(i + 3) % 4];
        let cur = n[i];
        let next = n[(i + 1) % 4];
        let v1 = [cur[0] - prev[0], cur[1] - prev[1]];
        let v2 = [next[0] - cur[0], next[1] - cur[1]];
        let d1 = (v1[0] * v1[0] + v1[1] * v1[1]).sqrt();
        let d2 = (v2[0] * v2[0] + v2[1] * v2[1]).sqrt();
        if d1 < 1e-15 || d2 < 1e-15 {
            return 0.0;
        }
        let dot = (v1[0] * v2[0] + v1[1] * v2[1]) / (d1 * d2);
        let angle = dot.clamp(-1.0, 1.0).acos().to_degrees();
        min_ang = min_ang.min(angle);
    }
    min_ang
}

/// Compute the average minimum angle across all quads.
#[allow(dead_code)]
fn avg_quad_quality(quads: &[[usize; 4]], nodes: &[[f64; 2]]) -> f64 {
    let n = quads.len() as f64;
    if n == 0.0 { return 0.0; }
    quads.iter().map(|q| quad_min_angle(q, nodes)).sum::<f64>() / n
}

// ─── Quality improvement via local diagonal swaps ────────────────────────────

/// Improve quad quality by trying alternative triangle pairings.
///
/// After the initial greedy recombination, some quads may be poorly shaped.
/// This pass looks for pairs of adjacent quads that share a node pair and
/// tests whether swapping the diagonal improves both quads.
///
/// The algorithm is a 2-opt style local search:
/// 1. Build a map from triangle index to the quad it belongs to.
/// 2. Find candidate swap pairs where two quads share a triangle edge.
/// 3. Test the alternative diagonal and accept if quality improves.
fn improve_quad_quality(
    quads: &mut Vec<[usize; 4]>,
    nodes: &[[f64; 2]],
    tris: &[[usize; 3]],
) -> usize {
    let n_quads = quads.len();
    if n_quads < 2 { return 0; }

    // Map each triangle index → quad index.
    // Build from the quad-to-triangle relationship.
    let _tri_to_quad = build_tri_to_quad_map(quads, tris);

    // Map each edge to the list of quads that contain it.
    let mut edge_to_quads: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (qi, q) in quads.iter().enumerate() {
        for i in 0..4 {
            let a = q[i];
            let b = q[(i + 1) % 4];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_to_quads.entry(key).or_default().push(qi);
        }
    }

    let mut improved = 0usize;
    let mut changed = true;

    // Iterate until no more improvements or max passes.
    for _pass in 0..10 {
        if !changed { break; }
        changed = false;

        // For each quad, try swapping with each neighbor sharing an edge.
        for qi in 0..quads.len() {
            let q = quads[qi];
            let q_nodes = [q[0], q[1], q[2], q[3]];

            // For each edge of this quad, find the adjacent quad.
            for ei in 0..4 {
                let a = q_nodes[ei];
                let b = q_nodes[(ei + 1) % 4];
                let key = if a < b { (a, b) } else { (b, a) };

                let Some(neighbors) = edge_to_quads.get(&key) else { continue };
                let qj = match neighbors.iter().filter(|&&nj| nj != qi).next() {
                    Some(&nj) => nj,
                    None => continue,
                };
                if qj == qi { continue; }

                let other = quads[qj];

                // Find the two nodes NOT on the shared edge in each quad.
                let shared = [a, b];
                let qi_opp: Vec<_> = q_nodes.iter().filter(|n| !shared.contains(n)).copied().collect();
                let qj_nodes = [other[0], other[1], other[2], other[3]];
                let qj_opp: Vec<_> = qj_nodes.iter().filter(|n| !shared.contains(n)).copied().collect();

                if qi_opp.len() != 2 || qj_opp.len() != 2 { continue; }

                let qi_a = qi_opp[0]; let qi_b = qi_opp[1];
                let qj_a = qj_opp[0]; let qj_b = qj_opp[1];

                // Current pairing: quads share edge (a,b).
                // Alternative: merge into the other diagonal (qi_a, qj_a) or (qi_a, qj_b).
                // Test both alternatives and pick the one that improves min angle more.
                let current_min = quad_min_angle(&q, nodes).min(quad_min_angle(&other, nodes));

                // Alternative 1: (qi_a, qj_a) becomes the shared edge.
                let alt1_quad1 = [shared[0], qi_a, shared[1], qj_a];
                let alt1_quad2 = [shared[0], qj_a, shared[1], qi_b];
                let alt1_quad3 = [shared[0], qi_b, shared[1], qj_b];
                // We need exactly 2 quads from these 3 nodes. Try all combos.
                let alt1_min = compute_swap_quality(
                    &[alt1_quad1, alt1_quad2], nodes, current_min
                ).max(
                    compute_swap_quality(&[alt1_quad1, alt1_quad3], nodes, current_min)
                ).max(
                    compute_swap_quality(&[alt1_quad2, alt1_quad3], nodes, current_min)
                );

                // Alternative 2: (qi_b, qj_b) becomes the shared edge.
                let alt2_quad1 = [shared[0], qi_a, shared[1], qj_b];
                let alt2_quad2 = [shared[0], qj_b, shared[1], qi_b];
                let alt2_quad3 = [shared[0], qi_b, shared[1], qj_a];
                let alt2_min = compute_swap_quality(
                    &[alt2_quad1, alt2_quad2], nodes, current_min
                ).max(
                    compute_swap_quality(&[alt2_quad1, alt2_quad3], nodes, current_min)
                ).max(
                    compute_swap_quality(&[alt2_quad2, alt2_quad3], nodes, current_min)
                );

                let best_min = alt1_min.max(alt2_min);

                if best_min > current_min + 1e-6 {
                    if alt1_min >= alt2_min {
                        // Apply alt1
                        let (new1, new2) = pick_best_pair(&[alt1_quad1, alt1_quad2, alt1_quad3], nodes);
                        if let (Some(n1), Some(n2)) = (new1, new2) {
                            quads[qi] = n1;
                            quads[qj] = n2;
                            improved += 1;
                            changed = true;
                        }
                    } else {
                        let (new1, new2) = pick_best_pair(&[alt2_quad1, alt2_quad2, alt2_quad3], nodes);
                        if let (Some(n1), Some(n2)) = (new1, new2) {
                            quads[qi] = n1;
                            quads[qj] = n2;
                            improved += 1;
                            changed = true;
                        }
                    }
                }
            }
        }

        // Rebuild edge map for next pass.
        if changed {
            edge_to_quads.clear();
            for (qj, q) in quads.iter().enumerate() {
                for i in 0..4 {
                    let a = q[i]; let b = q[(i + 1) % 4];
                    let key = if a < b { (a, b) } else { (b, a) };
                    edge_to_quads.entry(key).or_default().push(qj);
                }
            }
        }
    }

    improved
}

/// Compute the minimum angle among a set of quad candidates.
fn compute_swap_quality(candidates: &[[usize; 4]; 2], nodes: &[[f64; 2]], current: f64) -> f64 {
    let q0 = quad_min_angle(&candidates[0], nodes);
    let q1 = quad_min_angle(&candidates[1], nodes);
    if q0 < 1.0 || q1 < 1.0 { return current - 1.0; }
    q0.min(q1)
}

/// From a set of 3 candidate quads, pick the 2 that form the highest
/// minimum angle.
fn pick_best_pair(candidates: &[[usize; 4]; 3], nodes: &[[f64; 2]]) -> (Option<[usize; 4]>, Option<[usize; 4]>) {
    let combos = [
        (0, 1), (0, 2), (1, 2),
    ];
    let mut best: Option<(f64, usize, usize)> = None;
    for &(i, j) in &combos {
        let qi = quad_min_angle(&candidates[i], nodes);
        let qj = quad_min_angle(&candidates[j], nodes);
        if qi < 1.0 || qj < 1.0 { continue; }
        let min_ang = qi.min(qj);
        match best {
            Some((bm, _, _)) if min_ang <= bm => {}
            _ => best = Some((min_ang, i, j)),
        }
    }
    match best {
        Some((_, i, j)) => (Some(candidates[i]), Some(candidates[j])),
        None => (None, None),
    }
}

/// Build a map from each used triangle index → quad index.
fn build_tri_to_quad_map(quads: &[[usize; 4]], tris: &[[usize; 3]]) -> Vec<Option<usize>> {
    let mut map = vec![None; tris.len()];
    for (qi, q) in quads.iter().enumerate() {
        let nodes: HashSet<usize> = q.iter().copied().collect();
        for (ti, tri) in tris.iter().enumerate() {
            if map[ti].is_some() { continue; }
            if tri.iter().all(|n| nodes.contains(n)) {
                map[ti] = Some(qi);
            }
        }
    }
    map
}

// ─── Mesh assembly ───────────────────────────────────────────────────────────

/// Build the output Mesh from quads + remaining triangles.
fn build_quad_mesh(
    nodes: &[[f64; 2]],
    _node_ids: &[u64],
    tris: &[[usize; 3]],
    quads: &[[usize; 4]],
) -> Result<Mesh, MeshAlgoError> {
    let mut mesh = Mesh::new();

    // Add nodes (1-indexed by position in the nodes array)
    for (i, &pos) in nodes.iter().enumerate() {
        mesh.add_node(Node::new(i as u64 + 1, pos[0], pos[1], 0.0));
    }

    let mut tri_used = vec![false; tris.len()];

    for (elem_id, quad) in quads.iter().enumerate() {
        let mut found = Vec::new();
        for (ti, tri) in tris.iter().enumerate() {
            if tri_used[ti] { continue; }
            let count = quad.iter().filter(|q| tri.contains(q)).count();
            if count >= 3 { found.push(ti); }
        }
        if found.len() >= 2 {
            tri_used[found[0]] = true;
            tri_used[found[1]] = true;
            let nids: Vec<u64> = quad.iter().map(|&v| v as u64 + 1).collect();
            mesh.add_element(Element::new((elem_id + 1) as u64, ElementType::Quad4, nids));
        }
    }

    let mut next_elem_id = quads.len() as u64 + 1;
    for (ti, tri) in tris.iter().enumerate() {
        if tri_used[ti] { continue; }
        let nids = vec![tri[0] as u64 + 1, tri[1] as u64 + 1, tri[2] as u64 + 1];
        mesh.add_element(Element::new(next_elem_id, ElementType::Triangle3, nids));
        next_elem_id += 1;
    }

    if mesh.element_count() == 0 {
        return Err(MeshAlgoError::Generation(
            "no elements produced by quad recombination".to_string(),
        ));
    }

    Ok(mesh)
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn rect_domain() -> Domain2D {
        Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 2.0], [0.0, 2.0]])
    }

    #[test]
    fn frontal_quads_rectangle() {
        let domain = rect_domain();
        let params = MeshParams::with_size(0.8);
        let mesh = FrontalQuads2D::default().mesh_2d(&domain, &params).unwrap();
        assert!(mesh.element_count() > 0);
        let quad_count = mesh.elements.iter().filter(|e| e.etype == ElementType::Quad4).count();
        assert!(quad_count > 0, "should produce some quads");
    }

    #[test]
    fn frontal_quads_square() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);
        let params = MeshParams::with_size(0.6);
        let mesh = FrontalQuads2D::default().mesh_2d(&domain, &params).unwrap();
        assert!(mesh.element_count() >= 2);
    }

    #[test]
    fn quad_angle_rect_is_90() {
        // A perfect rectangle quad.
        let nodes = [[0.0, 0.0], [1.0, 0.0], [1.0, 0.5], [0.0, 0.5]];
        let quad = [0usize, 1, 2, 3];
        let angle = quad_min_angle(&quad, &nodes);
        assert!((angle - 90.0).abs() < 1e-10, "rectangle should be 90°, got {angle}");
    }

    #[test]
    fn quad_angle_skewed_is_less() {
        // A skewed quad — should have smaller min angle.
        let nodes = [[0.0, 0.0], [1.0, 0.0], [0.8, 0.5], [0.2, 0.5]];
        let quad = [0usize, 1, 2, 3];
        let angle = quad_min_angle(&quad, &nodes);
        assert!(angle < 90.0, "skewed quad should have angle < 90°, got {angle}");
        assert!(angle > 0.0, "angle should be positive");
    }

    #[test]
    fn improve_quality_does_not_degrade() {
        let domain = rect_domain();
        let params = MeshParams::with_size(0.8);
        let tri_mesh = FrontalDelaunay2D::default().mesh_2d(&domain, &params).unwrap();
        let (nodes, _node_ids, tris) = extract_tri_mesh_data(&tri_mesh).unwrap();
        if tris.len() < 4 { return; }

        let boundary_edges = extract_boundary_edges(&tris);
        let cf = CrossField::compute(&nodes, &tris, &boundary_edges, 100);
        let mut quads = recombine_triangles(&nodes, &tris, &cf);
        if quads.len() < 2 { return; }

        let before = avg_quad_quality(&quads, &nodes);
        let n_improved = improve_quad_quality(&mut quads, &nodes, &tris);
        let after = avg_quad_quality(&quads, &nodes);

        // Quality should not degrade.
        assert!(after >= before - 1e-10,
            "quality degraded: {before:.2}° → {after:.2}° ({n_improved} swaps)");
    }

    #[test]
    fn quad_dominant_output() {
        let domain = rect_domain();
        let params = MeshParams::with_size(0.6);
        let mesh = FrontalQuads2D::default().mesh_2d(&domain, &params).unwrap();
        let quads = mesh.elements.iter().filter(|e| e.etype == ElementType::Quad4).count();
        let tris = mesh.elements.iter().filter(|e| e.etype == ElementType::Triangle3).count();
        assert!(quads > 0, "should produce at least some quads; got {quads} quads, {tris} tris");
    }
}
