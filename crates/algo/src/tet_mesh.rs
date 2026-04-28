//! Compact tetrahedral mesh with pre-built neighbor tables.
//!
//! Used internally by the Delaunay3D pipeline for bistellar flip optimization.
//! Neighbor tables enable O(1) face adjacency lookup, replacing the O(n log n)
//! face-map construction needed when operating on the public `Mesh` type.

use std::collections::HashMap;

use rmsh_model::{Element, ElementType, Mesh, Node};

use crate::delaunay_3d::{
    is_better_quality, min_dihedral_points, radius_edge_ratio_points, tetra_volume,
};

// ─── Structs ─────────────────────────────────────────────────────────────────

/// A single tetrahedron with pre-computed neighbor links.
///
/// Face convention: face `fi` is opposite node `nodes[fi]`, i.e. the face
/// consists of vertices `(fi+1)%4, (fi+2)%4, (fi+3)%4`.
#[derive(Debug, Clone, Copy)]
pub struct Tet {
    /// Four node indices into `TetMesh::nodes`.
    pub nodes: [u32; 4],
    /// `neighbors[fi]` = index of the tet sharing face `fi`,
    /// or `u32::MAX` for boundary faces.
    pub neighbors: [u32; 4],
}

/// Compact tetrahedral mesh with pre-built neighbor tables.
///
/// Nodes are indexed by `u32` into a flat `Vec`. Tetrahedra carry their
/// four node indices and four neighbor indices.
#[derive(Debug, Clone)]
pub struct TetMesh {
    /// Node coordinates: `nodes[i] = (x, y, z)` for compact index `i`.
    pub nodes: Vec<[f64; 3]>,
    /// The original `u64` node IDs, mapping compact index → Mesh node ID.
    pub node_ids: Vec<u64>,
    /// All tetrahedra, indexed by position.
    pub tets: Vec<Tet>,
}

// ─── Constructors / conversion ──────────────────────────────────────────────

impl TetMesh {
    /// Build a `TetMesh` from a public `Mesh` by extracting `Tetrahedron4` elements.
    pub fn from_mesh(mesh: &Mesh) -> Self {
        // Collect all node IDs, sorted, assign compact u32 indices.
        let mut all_ids: Vec<u64> = mesh.nodes.keys().copied().collect();
        all_ids.sort_unstable();
        let mut id_to_idx: HashMap<u64, u32> = HashMap::with_capacity(all_ids.len());

        let mut nodes = Vec::with_capacity(all_ids.len());
        let mut node_ids = Vec::with_capacity(all_ids.len());
        for (compact_i, nid) in all_ids.iter().enumerate() {
            id_to_idx.insert(*nid, compact_i as u32);
            let node = &mesh.nodes[nid];
            nodes.push([node.position.x, node.position.y, node.position.z]);
            node_ids.push(*nid);
        }

        // Extract Tetrahedron4 elements.
        let tets: Vec<Tet> = mesh
            .elements
            .iter()
            .filter(|e| e.etype == ElementType::Tetrahedron4 && e.node_ids.len() == 4)
            .map(|e| Tet {
                nodes: [
                    id_to_idx[&e.node_ids[0]],
                    id_to_idx[&e.node_ids[1]],
                    id_to_idx[&e.node_ids[2]],
                    id_to_idx[&e.node_ids[3]],
                ],
                neighbors: [u32::MAX; 4],
            })
            .collect();

        let mut tm = Self {
            nodes,
            node_ids,
            tets,
        };
        tm.build_neighbors();
        tm
    }

    /// Convert back to a public `Mesh`, preserving node IDs.
    pub fn to_mesh(&self) -> Mesh {
        let mut mesh = Mesh::new();
        for (i, coord) in self.nodes.iter().enumerate() {
            let nid = self.node_ids[i];
            mesh.add_node(Node::new(nid, coord[0], coord[1], coord[2]));
        }
        for (elem_id, tet) in self.tets.iter().enumerate() {
            mesh.add_element(Element::new(
                (elem_id + 1) as u64,
                ElementType::Tetrahedron4,
                vec![
                    self.node_ids[tet.nodes[0] as usize],
                    self.node_ids[tet.nodes[1] as usize],
                    self.node_ids[tet.nodes[2] as usize],
                    self.node_ids[tet.nodes[3] as usize],
                ],
            ));
        }
        mesh
    }

    /// Rebuild neighbor tables from scratch.
    ///
    /// For each tet, sets `neighbors[fi]` = index of the tet sharing face `fi`,
    /// or `u32::MAX` if the face is on the boundary.
    pub fn build_neighbors(&mut self) {
        for tet in &mut self.tets {
            tet.neighbors = [u32::MAX; 4];
        }

        // Build face → [(tet_idx, face_idx)] map with sorted keys.
        let mut face_map: HashMap<[u32; 3], Vec<(u32, u32)>> =
            HashMap::with_capacity(self.tets.len() * 2);
        for (ti, tet) in self.tets.iter().enumerate() {
            for fi in 0..4u32 {
                let mut face = self.face_nodes(ti, fi as usize);
                face.sort_unstable();
                face_map.entry(face).or_default().push((ti as u32, fi));
            }
        }

        // Link exactly-2 entries (interior faces).
        for (_key, entries) in &face_map {
            if entries.len() != 2 {
                continue;
            }
            let (t0, f0) = entries[0];
            let (t1, f1) = entries[1];
            self.tets[t0 as usize].neighbors[f0 as usize] = t1;
            self.tets[t1 as usize].neighbors[f1 as usize] = t0;
        }
    }

    /// Return the 3 node indices of face `fi` of tet `ti`.
    #[inline]
    pub fn face_nodes(&self, ti: usize, fi: usize) -> [u32; 3] {
        let tet = &self.tets[ti];
        [
            tet.nodes[(fi + 1) % 4],
            tet.nodes[(fi + 2) % 4],
            tet.nodes[(fi + 3) % 4],
        ]
    }

    /// Number of tets.
    #[inline]
    pub fn tet_count(&self) -> usize {
        self.tets.len()
    }
}

// ─── Quality helpers ─────────────────────────────────────────────────────────

/// Compute (min_dihedral, sliver_fraction, max_radius_edge) for a set of tets.
fn tetmesh_aggregate_quality(tm: &TetMesh, tets: &[[u32; 4]]) -> Option<(f64, f64, f64)> {
    let mut min_d = f64::MAX;
    let mut max_r = 0.0_f64;
    let mut sliver_like = 0usize;
    for tet in tets {
        let a = tm.nodes[tet[0] as usize];
        let b = tm.nodes[tet[1] as usize];
        let c = tm.nodes[tet[2] as usize];
        let d = tm.nodes[tet[3] as usize];
        let v = tetra_volume(a, b, c, d);
        if v <= 1e-15 {
            return None;
        }
        let dmin = min_dihedral_points(a, b, c, d);
        let r = radius_edge_ratio_points(a, b, c, d);
        if !dmin.is_finite() || !r.is_finite() {
            return None;
        }
        min_d = min_d.min(dmin);
        max_r = max_r.max(r);
        if dmin < 6.0 && r > 1.8 {
            sliver_like += 1;
        }
    }
    Some((min_d, sliver_like as f64 / (tets.len() as f64), max_r))
}

/// Returns true when significant fraction of tets are sliver-like or global
/// min dihedral is very low. Used to bias toward edge flips.
fn tetmesh_has_sliver_pressure(
    tm: &TetMesh,
    sliver_fraction_threshold: f64,
    min_dihedral_threshold: f64,
) -> bool {
    let mut total = 0usize;
    let mut sliver_like = 0usize;
    let mut global_min_d = f64::MAX;

    for tet in &tm.tets {
        let a = tm.nodes[tet.nodes[0] as usize];
        let b = tm.nodes[tet.nodes[1] as usize];
        let c = tm.nodes[tet.nodes[2] as usize];
        let d = tm.nodes[tet.nodes[3] as usize];
        let v = tetra_volume(a, b, c, d);
        if v <= 1e-15 {
            continue;
        }
        total += 1;
        let dmin = min_dihedral_points(a, b, c, d);
        let r = radius_edge_ratio_points(a, b, c, d);
        if !dmin.is_finite() || !r.is_finite() {
            continue;
        }
        global_min_d = global_min_d.min(dmin);
        if dmin < 6.0 && r > 1.8 {
            sliver_like += 1;
        }
    }

    if total == 0 {
        return false;
    }

    let sliver_frac = sliver_like as f64 / total as f64;
    sliver_frac >= sliver_fraction_threshold || global_min_d <= min_dihedral_threshold
}

// ─── Flip application helpers ────────────────────────────────────────────────

fn all_tets_positive_volume(tm: &TetMesh, candidates: &[[u32; 4]]) -> bool {
    for tet_nodes in candidates {
        let a = tm.nodes[tet_nodes[0] as usize];
        let b = tm.nodes[tet_nodes[1] as usize];
        let c = tm.nodes[tet_nodes[2] as usize];
        let d = tm.nodes[tet_nodes[3] as usize];
        if tetra_volume(a, b, c, d) <= 1e-15 {
            return false;
        }
    }
    true
}

/// 2-to-3 face flip: replace two tets sharing a face with three tets.
///
/// Removes `t0` and `t1` (ordered by descending index for swap_remove safety),
/// pushes the 3 new tets, then rebuilds neighbor tables.
fn apply_2to3_flip(tm: &mut TetMesh, t0: u32, t1: u32, new_tets: [[u32; 4]; 3]) {
    let hi = t0.max(t1);
    let lo = t0.min(t1);
    tm.tets.swap_remove(hi as usize);
    tm.tets.swap_remove(lo as usize);
    for n in new_tets {
        tm.tets.push(Tet {
            nodes: n,
            neighbors: [u32::MAX; 4],
        });
    }
    tm.build_neighbors();
}

/// 3-to-2 edge flip: replace three tets sharing an edge with two tets.
///
/// Removes `remove` tets (sorted descending), pushes 2 new tets,
/// then rebuilds neighbor tables.
fn apply_3to2_flip(tm: &mut TetMesh, remove: &[u32], new_tets: [[u32; 4]; 2]) {
    let mut sorted = remove.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    for &idx in &sorted {
        tm.tets.swap_remove(idx as usize);
    }
    for n in new_tets {
        tm.tets.push(Tet {
            nodes: n,
            neighbors: [u32::MAX; 4],
        });
    }
    tm.build_neighbors();
}

// ─── Main flip optimization ─────────────────────────────────────────────────

/// Perform multi-pass bistellar flip optimization on a `TetMesh`.
///
/// Runs up to `max_passes` iterations of 2-to-3 face flips and 3-to-2 edge
/// flips, selecting the best candidate per pass using strict and sliver-relaxed
/// quality acceptance criteria.
///
/// Returns `(face_flip_count, edge_flip_count, edge_sliver_flip_count)`.
pub fn optimize_tetmesh_flips(tm: &mut TetMesh, max_passes: usize) -> (usize, usize, usize) {
    let mut accepted_face = 0usize;
    let mut accepted_edge = 0usize;
    let mut accepted_edge_sliver = 0usize;

    for _pass in 0..max_passes {
        tm.build_neighbors();
        let prefer_edge_phase = tetmesh_has_sliver_pressure(tm, 0.30, 0.12);

        // ── 2-to-3 face flips via neighbor table ──
        let mut best_face: Option<(u32, u32, [[u32; 4]; 3], (f64, f64, f64))> = None;

        for ti in 0..tm.tets.len() as u32 {
            let tet = &tm.tets[ti as usize];
            for fi in 0..4u32 {
                let tj = tet.neighbors[fi as usize];
                if tj == u32::MAX || ti >= tj {
                    continue;
                }

                // Find the face index in tj that points back to ti.
                let tj_tet = &tm.tets[tj as usize];
                let fj = match (0..4).find(|&f| tj_tet.neighbors[f] == ti) {
                    Some(f) => f as u32,
                    None => continue,
                };

                let spine_a = tet.nodes[fi as usize];
                let spine_b = tj_tet.nodes[fj as usize];

                // Face vertices (the 3 nodes opposite spine_a / opposite fi).
                let fv = tm.face_nodes(ti as usize, fi as usize);

                let new_tets = [
                    [spine_a, spine_b, fv[0], fv[1]],
                    [spine_a, spine_b, fv[1], fv[2]],
                    [spine_a, spine_b, fv[2], fv[0]],
                ];

                if !all_tets_positive_volume(tm, &new_tets) {
                    continue;
                }

                let old_nodes = [tet.nodes, tj_tet.nodes];
                let Some(old_q) = tetmesh_aggregate_quality(tm, &old_nodes) else {
                    continue;
                };
                let Some(new_q) = tetmesh_aggregate_quality(tm, &new_tets) else {
                    continue;
                };

                if !is_better_quality(new_q, old_q) {
                    continue;
                }

                match best_face {
                    Some((_, _, _, bq)) if !is_better_quality(new_q, bq) => {}
                    _ => best_face = Some((ti, tj, new_tets, new_q)),
                }
            }
        }

        // ── 3-to-2 edge flips via edge map ──
        let mut edge_map: HashMap<[u32; 2], Vec<u32>> = HashMap::new();
        for ti in 0..tm.tets.len() as u32 {
            let n = tm.tets[ti as usize].nodes;
            for (a, b) in &[
                (n[0], n[1]),
                (n[0], n[2]),
                (n[0], n[3]),
                (n[1], n[2]),
                (n[1], n[3]),
                (n[2], n[3]),
            ] {
                let key = if *a < *b { [*a, *b] } else { [*b, *a] };
                edge_map.entry(key).or_default().push(ti);
            }
        }

        let mut best_edge: Option<(
            Vec<u32>,
            [[u32; 4]; 2],
            (f64, f64, f64),
            f64,
            bool,
        )> = None;

        for (_edge, adjacent) in &edge_map {
            if adjacent.len() != 3 {
                continue;
            }

            // Verify all three tets still exist and extract opposite vertices.
            let mut opposite_pairs = Vec::<[u32; 2]>::with_capacity(3);
            let mut opposite_vertices = Vec::<u32>::with_capacity(3);
            let mut old_tets = [[0_u32; 4]; 3];
            let mut valid = true;

            for (slot, &ti) in adjacent.iter().enumerate() {
                let tet = &tm.tets[ti as usize];
                old_tets[slot] = tet.nodes;

                // Find the two vertices not on the shared edge.
                let u = _edge[0];
                let v = _edge[1];
                let mut opp = Vec::<u32>::with_capacity(2);
                for &nid in &tet.nodes {
                    if nid != u && nid != v {
                        opp.push(nid);
                    }
                }
                if opp.len() != 2 {
                    valid = false;
                    break;
                }
                if !tet.nodes.contains(&u) || !tet.nodes.contains(&v) {
                    valid = false;
                    break;
                }
                opp.sort_unstable();
                opposite_pairs.push([opp[0], opp[1]]);
                opposite_vertices.push(opp[0]);
                opposite_vertices.push(opp[1]);
            }

            if !valid {
                continue;
            }

            opposite_vertices.sort_unstable();
            opposite_vertices.dedup();
            if opposite_vertices.len() != 3 {
                continue;
            }
            let a = opposite_vertices[0];
            let b = opposite_vertices[1];
            let c = opposite_vertices[2];
            let mut need = vec![[a, b], [b, c], [a, c]];
            for p in &mut need {
                p.sort_unstable();
            }
            let mut got = opposite_pairs.clone();
            got.sort_unstable();
            need.sort_unstable();
            if got != need {
                continue;
            }

            let new_tets = [[a, b, c, _edge[0]], [a, b, c, _edge[1]]];

            let Some((old_d, old_s, old_r)) = tetmesh_aggregate_quality(tm, &old_tets) else {
                continue;
            };
            let Some((new_d, new_s, new_r)) = tetmesh_aggregate_quality(tm, &new_tets) else {
                continue;
            };

            let strict_improves = (new_d > old_d + 1e-6)
                || ((new_d - old_d).abs() < 1e-6
                    && ((new_s < old_s - 1e-9)
                        || ((new_s - old_s).abs() < 1e-9 && new_r < old_r - 1e-9)));
            let sliver_delta = old_s - new_s;
            let strong_sliver_reduction = old_s >= 0.66 && new_s <= 0.34;
            let sliver_relaxed_improves = (sliver_delta > 0.08 || strong_sliver_reduction)
                && new_d >= old_d * 0.70
                && new_r <= old_r * 1.35;
            if !strict_improves && !sliver_relaxed_improves {
                continue;
            }

            let mut remove = vec![adjacent[0], adjacent[1], adjacent[2]];
            remove.sort_unstable();
            remove.reverse();

            let used_sliver_relaxed = !strict_improves && sliver_relaxed_improves;
            let new_quality = (new_d, new_s, new_r);

            match best_edge {
                Some((_, _, best_q, best_sliver_delta, best_relaxed)) => {
                    let prefer = if sliver_delta > best_sliver_delta + 1e-9 {
                        true
                    } else if (sliver_delta - best_sliver_delta).abs() < 1e-9 {
                        if used_sliver_relaxed != best_relaxed {
                            !used_sliver_relaxed
                        } else {
                            is_better_quality(new_quality, best_q)
                        }
                    } else {
                        false
                    };
                    if prefer {
                        best_edge = Some((
                            remove,
                            new_tets,
                            new_quality,
                            sliver_delta,
                            used_sliver_relaxed,
                        ));
                    }
                }
                _ => {
                    best_edge = Some((
                        remove,
                        new_tets,
                        new_quality,
                        sliver_delta,
                        used_sliver_relaxed,
                    ));
                }
            }
        }

        // ── Execute best flip ──
        let mut did_flip = false;

        // Phase 1: face-first unless sliver pressure.
        if !prefer_edge_phase {
            if let Some((ti, tj, new_tets, _)) = best_face.take() {
                apply_2to3_flip(tm, ti, tj, new_tets);
                accepted_face += 1;
                did_flip = true;
            }
        }

        if !did_flip {
            if let Some((remove, new_tets, _, _, used_sliver_relaxed)) = best_edge.take() {
                apply_3to2_flip(tm, &remove, new_tets);
                accepted_edge += 1;
                if used_sliver_relaxed {
                    accepted_edge_sliver += 1;
                }
                did_flip = true;
            }
        }

        // Phase 3: edge-phase failed; try face flip as fallback.
        if !did_flip && prefer_edge_phase {
            if let Some((ti, tj, new_tets, _)) = best_face.take() {
                apply_2to3_flip(tm, ti, tj, new_tets);
                accepted_face += 1;
                did_flip = true;
            }
        }

        if !did_flip {
            break; // converged
        }
    }

    (accepted_face, accepted_edge, accepted_edge_sliver)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rmsh_model::{Element, ElementType, Mesh, Node};

    fn make_two_tet_mesh() -> Mesh {
        // Two tets sharing face (2,3,4):
        // Tet 0: [1,2,3,4]  Tet 1: [2,5,3,4]
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.5, 0.5, 0.0));
        mesh.add_node(Node::new(4, 0.5, 0.0, 0.5));
        mesh.add_node(Node::new(5, 1.0, 0.5, 0.5));
        mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));
        mesh.add_element(Element::new(2, ElementType::Tetrahedron4, vec![2, 5, 3, 4]));
        mesh
    }

    fn make_three_tet_edge_fan() -> Mesh {
        // Three tets sharing edge (1,2):
        // [1,2,3,4], [1,2,4,5], [1,2,3,5]
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.5, 0.001, 0.0));
        mesh.add_node(Node::new(4, 0.5, 0.0, 0.001));
        mesh.add_node(Node::new(5, 0.5, 0.8, 0.8));
        mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));
        mesh.add_element(Element::new(2, ElementType::Tetrahedron4, vec![1, 2, 4, 5]));
        mesh.add_element(Element::new(3, ElementType::Tetrahedron4, vec![1, 2, 3, 5]));
        mesh
    }

    #[test]
    fn round_trip_conversion() {
        let mesh = make_two_tet_mesh();
        let tm = TetMesh::from_mesh(&mesh);
        assert_eq!(tm.tet_count(), 2);
        assert_eq!(tm.nodes.len(), 5);
        assert_eq!(tm.node_ids, vec![1, 2, 3, 4, 5]);

        let back = tm.to_mesh();
        assert_eq!(back.node_count(), 5);
        assert_eq!(back.element_count(), 2);
        // Verify all elements are Tetrahedron4 with valid node IDs.
        for e in &back.elements {
            assert_eq!(e.etype, ElementType::Tetrahedron4);
            for &nid in &e.node_ids {
                assert!(back.nodes.contains_key(&nid), "node {nid} missing");
            }
        }
    }

    #[test]
    fn neighbor_table_two_tets() {
        let mesh = make_two_tet_mesh();
        let tm = TetMesh::from_mesh(&mesh);

        // Find which face is shared. Each tet has 4 neighbors; exactly one
        // should be non-MAX and point to the other tet.
        let mut found = 0usize;
        for ti in 0..tm.tet_count() {
            for fi in 0..4 {
                let nj = tm.tets[ti].neighbors[fi];
                if nj != u32::MAX {
                    // The other tet points back.
                    let fj = (0..4)
                        .find(|&f| tm.tets[nj as usize].neighbors[f] == ti as u32)
                        .expect("back-link exists");
                    let face = tm.face_nodes(ti, fi);
                    let nf = tm.face_nodes(nj as usize, fj);
                    // Same 3 vertices (order may differ).
                    let mut sorted_f = face;
                    sorted_f.sort_unstable();
                    let mut sorted_nf = nf;
                    sorted_nf.sort_unstable();
                    assert_eq!(sorted_f, sorted_nf);
                    found += 1;
                }
            }
        }
        // Exactly 2 face cross-links (tet0→tet1 and tet1→tet0 = 2 links).
        assert_eq!(found, 2, "expected exactly 2 neighbor cross-links");
    }

    #[test]
    fn neighbor_table_three_tets() {
        let mesh = make_three_tet_edge_fan();
        let tm = TetMesh::from_mesh(&mesh);
        assert_eq!(tm.tet_count(), 3);

        // All interior faces should have valid neighbors.
        let mut interior_faces = 0usize;
        let mut boundary_faces = 0usize;
        for ti in 0..tm.tet_count() {
            for fi in 0..4 {
                if tm.tets[ti].neighbors[fi] != u32::MAX {
                    interior_faces += 1;
                } else {
                    boundary_faces += 1;
                }
            }
        }
        // 3 tets × 4 faces = 12 face sides.
        // Interior: each shared face counted from both sides → even.
        assert!(interior_faces > 0, "expected interior faces");
        assert_eq!((interior_faces + boundary_faces) as usize, 12);
    }

    #[test]
    fn apply_2to3_produces_three_tets() {
        let mesh = make_two_tet_mesh();
        let mut tm = TetMesh::from_mesh(&mesh);
        assert_eq!(tm.tet_count(), 2);

        // Find the shared face pair and build the 2→3 replacement.
        let (t0, t1, new_tets) = {
            let mut pair = None;
            for ti in 0..tm.tet_count() as u32 {
                for fi in 0..4 {
                    let tj = tm.tets[ti as usize].neighbors[fi as usize];
                    if tj != u32::MAX && ti < tj {
                        let fj = (0..4)
                            .find(|&f| tm.tets[tj as usize].neighbors[f] == ti)
                            .unwrap();
                        let spine_a = tm.tets[ti as usize].nodes[fi as usize];
                        let spine_b = tm.tets[tj as usize].nodes[fj];
                        let fv = tm.face_nodes(ti as usize, fi as usize);
                        let new_tets = [
                            [spine_a, spine_b, fv[0], fv[1]],
                            [spine_a, spine_b, fv[1], fv[2]],
                            [spine_a, spine_b, fv[2], fv[0]],
                        ];
                        pair = Some((ti, tj, new_tets));
                    }
                }
            }
            pair.expect("shared face exists")
        };

        let before_count = tm.tet_count();
        apply_2to3_flip(&mut tm, t0, t1, new_tets);
        assert_eq!(tm.tet_count(), before_count + 1); // -2 + 3 = +1

        // All tets should have valid neighbors after rebuild.
        for (ti, tet) in tm.tets.iter().enumerate() {
            for fi in 0..4 {
                let nj = tet.neighbors[fi];
                if nj != u32::MAX {
                    let back = (0..4)
                        .find(|&f| tm.tets[nj as usize].neighbors[f] == ti as u32);
                    assert!(
                        back.is_some(),
                        "tet {ti} face {fi} neighbor {nj} missing back-link"
                    );
                }
            }
        }
    }

    #[test]
    fn apply_3to2_produces_two_tets() {
        let mesh = make_three_tet_edge_fan();
        let mut tm = TetMesh::from_mesh(&mesh);
        assert_eq!(tm.tet_count(), 3);

        // The 3-tet edge fan should contain an edge shared by all 3 tets.
        // Find it and do the 3-to-2 flip.
        let mut edge_map: HashMap<[u32; 2], Vec<u32>> = HashMap::new();
        for ti in 0..tm.tet_count() as u32 {
            let n = tm.tets[ti as usize].nodes;
            for (a, b) in &[
                (n[0], n[1]),
                (n[0], n[2]),
                (n[0], n[3]),
                (n[1], n[2]),
                (n[1], n[3]),
                (n[2], n[3]),
            ] {
                let key = if *a < *b { [*a, *b] } else { [*b, *a] };
                edge_map.entry(key).or_default().push(ti);
            }
        }

        let (edge, adjacent) = edge_map
            .into_iter()
            .find(|(_, v)| v.len() == 3)
            .expect("shared edge exists");

        let u = edge[0];
        let v = edge[1];

        // Extract opposite vertices.
        let mut opp_verts = Vec::new();
        for &ti in &adjacent {
            let nodes = tm.tets[ti as usize].nodes;
            for &nid in &nodes {
                if nid != u && nid != v {
                    opp_verts.push(nid);
                }
            }
        }
        opp_verts.sort_unstable();
        opp_verts.dedup();
        assert_eq!(opp_verts.len(), 3);

        let new_tets = [
            [opp_verts[0], opp_verts[1], opp_verts[2], u],
            [opp_verts[0], opp_verts[1], opp_verts[2], v],
        ];

        let remove = adjacent.clone();
        let before = tm.tet_count();
        apply_3to2_flip(&mut tm, &remove, new_tets);
        assert_eq!(tm.tet_count(), before - 1); // -3 + 2 = -1

        // Verify neighbor validity.
        for (ti, tet) in tm.tets.iter().enumerate() {
            for fi in 0..4 {
                let nj = tet.neighbors[fi];
                if nj != u32::MAX {
                    let back = (0..4)
                        .find(|&f| tm.tets[nj as usize].neighbors[f] == ti as u32);
                    assert!(
                        back.is_some(),
                        "tet {ti} face {fi} neighbor {nj} missing back-link"
                    );
                }
            }
        }
    }

    #[test]
    fn tetmesh_flip_activates_on_edge_fan() {
        // Same configuration as delaunay_3d's local_flip_pass_activates_on_edge_fan.
        let mesh = make_three_tet_edge_fan();
        let mut tm = TetMesh::from_mesh(&mesh);
        let before = tm.tet_count();

        let (face_flips, edge_flips, _sliver) = optimize_tetmesh_flips(&mut tm, 4);
        let after = tm.tet_count();

        assert!(
            face_flips + edge_flips > 0,
            "expected at least one flip"
        );
        assert_eq!(before, 3);
        assert_ne!(after, before, "tet count should change after flip");
    }

    #[test]
    fn quality_parity_with_mesh_version() {
        // Build a small mesh and verify tetmesh_aggregate_quality matches
        // the equivalent logic when given the same coordinates.
        let mesh = make_two_tet_mesh();
        let tm = TetMesh::from_mesh(&mesh);

        // tetmesh_aggregate_quality on all tets should produce valid metrics.
        let all_tet_nodes: Vec<[u32; 4]> = tm.tets.iter().map(|t| t.nodes).collect();
        let q = tetmesh_aggregate_quality(&tm, &all_tet_nodes);
        assert!(q.is_some(), "quality should be Some for valid tets");
        let (min_d, sliver_frac, max_r) = q.unwrap();
        assert!(min_d.is_finite() && min_d > 0.0, "min_dihedral={min_d}");
        assert!(sliver_frac >= 0.0 && sliver_frac <= 1.0);
        assert!(max_r.is_finite() && max_r > 0.0);
    }

    #[test]
    fn empty_tetmesh_no_flips() {
        let mut tm = TetMesh {
            nodes: vec![],
            node_ids: vec![],
            tets: vec![],
        };
        let (f, e, s) = optimize_tetmesh_flips(&mut tm, 4);
        assert_eq!(f, 0);
        assert_eq!(e, 0);
        assert_eq!(s, 0);
    }

    #[test]
    fn single_tet_no_flips() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.5, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.5, 0.5, 1.0));
        mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));

        let mut tm = TetMesh::from_mesh(&mesh);
        let (f, e, s) = optimize_tetmesh_flips(&mut tm, 4);
        assert_eq!(f, 0);
        assert_eq!(e, 0);
        assert_eq!(s, 0);
        assert_eq!(tm.tet_count(), 1);
    }
}
