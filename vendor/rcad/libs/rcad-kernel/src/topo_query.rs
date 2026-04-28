//! Topology query helpers — analogous to OCCT `TopExp_Explorer` and
//! `TopExp::MapShapesAndAncestors`.
//!
//! All functions operate on `solids[0].shells[0]` for face-level queries and
//! on `brep.edges` / `brep.vertices` for edge/vertex-level queries.
//! They are safe to call on an empty BRep (return 0 or empty Vec).

use crate::BRep;

/// Returns the flat face indices (in `solids[0].shells[0].faces`) that
/// reference `edge_idx` in their outer wire.
///
/// For a well-formed manifold solid each edge is shared by exactly two faces.
/// Returns fewer if the edge is a boundary edge, more if non-manifold.
pub fn edge_adjacent_faces(brep: &BRep, edge_idx: usize) -> Vec<usize> {
    let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
        Some(sh) => sh,
        None => return Vec::new(),
    };
    shell
        .faces
        .iter()
        .enumerate()
        .filter(|(_, face)| face.outer_wire.edges.iter().any(|we| we.idx == edge_idx))
        .map(|(fi, _)| fi)
        .collect()
}

/// Returns all edge indices referenced in the outer wire of face `face_idx`
/// (in `solids[0].shells[0].faces[face_idx]`).
///
/// Duplicate edge indices (e.g. seam edges) are preserved as they appear
/// in the wire.
pub fn face_edges(brep: &BRep, face_idx: usize) -> Vec<usize> {
    brep.solids
        .first()
        .and_then(|s| s.shells.first())
        .and_then(|sh| sh.faces.get(face_idx))
        .map(|face| face.outer_wire.edges.iter().map(|we| we.idx).collect())
        .unwrap_or_default()
}

/// Returns all edge indices `ei` where `brep.edges[ei].start == vertex_idx`
/// or `brep.edges[ei].end == vertex_idx`.
pub fn vertex_adjacent_edges(brep: &BRep, vertex_idx: usize) -> Vec<usize> {
    brep.edges
        .iter()
        .enumerate()
        .filter(|(_, e)| e.start == vertex_idx || e.end == vertex_idx)
        .map(|(ei, _)| ei)
        .collect()
}

/// Number of faces in `solids[0].shells[0]`.
pub fn face_count(brep: &BRep) -> usize {
    brep.solids
        .first()
        .and_then(|s| s.shells.first())
        .map(|sh| sh.faces.len())
        .unwrap_or(0)
}

/// Number of edges in `brep.edges`.
pub fn edge_count(brep: &BRep) -> usize {
    brep.edges.len()
}

/// Number of vertices in `brep.vertices`.
pub fn vertex_count(brep: &BRep) -> usize {
    brep.vertices.len()
}

/// Returns `true` if `edge_idx` is a degenerate edge — i.e. its start and end
/// vertices are the same point (within floating-point equality), or it is
/// explicitly flagged degenerate in `brep.geom.edge_degenerated`.
///
/// Analogous to `BRep_Tool::Degenerated(edge)` in OCCT.
pub fn is_degenerate_edge(brep: &BRep, edge_idx: usize) -> bool {
    // Honour explicit flag first.
    if brep
        .geom
        .edge_degenerated
        .get(edge_idx)
        .copied()
        .unwrap_or(false)
    {
        return true;
    }
    let Some(edge) = brep.edges.get(edge_idx) else {
        return false;
    };
    let Some(v_start) = brep.vertices.get(edge.start) else {
        return false;
    };
    let Some(v_end) = brep.vertices.get(edge.end) else {
        return false;
    };
    (v_end.point - v_start.point).length_squared() < 1e-20
}

/// Returns all edge indices that are candidates for seam edges — edges whose
/// start and end vertices are the same (or nearly the same) 3D point, but that
/// are *not* flagged as degenerate.
///
/// Seam edges appear on surfaces of revolution (cylinder, sphere, cone, torus)
/// where the UV seam is a real geometric curve even though it starts and ends
/// at the same 3D location.
///
/// Analogous to walking seam edges on `BRep_Tool::IsClosed` surfaces in OCCT.
pub fn seam_edge_candidates(brep: &BRep) -> Vec<usize> {
    brep.edges
        .iter()
        .enumerate()
        .filter(|(ei, edge)| {
            // Must not be degenerate.
            if brep
                .geom
                .edge_degenerated
                .get(*ei)
                .copied()
                .unwrap_or(false)
            {
                return false;
            }
            let Some(v0) = brep.vertices.get(edge.start) else {
                return false;
            };
            let Some(v1) = brep.vertices.get(edge.end) else {
                return false;
            };
            // Seam: start ≈ end but not formally degenerate (has non-zero curve length).
            (v1.point - v0.point).length_squared() < 1e-10
        })
        .map(|(ei, _)| ei)
        .collect()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::PrimitiveSolid;

    fn box_2x2x2() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        })
    }

    #[test]
    fn box_counts() {
        let brep = box_2x2x2();
        assert_eq!(face_count(&brep), 6);
        assert_eq!(edge_count(&brep), 12);
        assert_eq!(vertex_count(&brep), 8);
    }

    #[test]
    fn box_every_edge_has_two_adjacent_faces() {
        let brep = box_2x2x2();
        for ei in 0..edge_count(&brep) {
            let adj = edge_adjacent_faces(&brep, ei);
            assert_eq!(
                adj.len(),
                2,
                "edge {ei} should have 2 adjacent faces, got {:?}",
                adj
            );
        }
    }

    #[test]
    fn box_every_vertex_has_three_adjacent_edges() {
        let brep = box_2x2x2();
        for vi in 0..vertex_count(&brep) {
            let adj = vertex_adjacent_edges(&brep, vi);
            assert_eq!(
                adj.len(),
                3,
                "vertex {vi} should have 3 adjacent edges, got {:?}",
                adj
            );
        }
    }

    #[test]
    fn box_face_has_four_edges() {
        let brep = box_2x2x2();
        // Each face of a box has 4 edges
        for fi in 0..face_count(&brep) {
            let edges = face_edges(&brep, fi);
            assert_eq!(
                edges.len(),
                4,
                "face {fi} should have 4 edges, got {:?}",
                edges
            );
        }
    }

    #[test]
    fn empty_brep_returns_zeros() {
        let brep = BRep::new();
        assert_eq!(face_count(&brep), 0);
        assert_eq!(edge_count(&brep), 0);
        assert_eq!(vertex_count(&brep), 0);
        assert!(edge_adjacent_faces(&brep, 0).is_empty());
        assert!(vertex_adjacent_edges(&brep, 0).is_empty());
        assert!(face_edges(&brep, 0).is_empty());
    }

    #[test]
    fn is_degenerate_edge_zero_length() {
        use crate::topology::Vertex;
        use glam::DVec3;

        let mut brep = BRep::new();
        // Edge 0: normal non-degenerate edge (different endpoints)
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(crate::topology::Edge { start: 0, end: 1 });
        // Edge 1: degenerate — same start and end vertex
        brep.vertices.push(Vertex { point: DVec3::Y });
        brep.vertices.push(Vertex { point: DVec3::Y }); // same point
        brep.edges.push(crate::topology::Edge { start: 2, end: 3 });
        // Edge 2: explicitly flagged degenerate
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::Z });
        brep.edges.push(crate::topology::Edge { start: 4, end: 5 });
        brep.geom.edge_degenerated = vec![false, false, true];

        assert!(!is_degenerate_edge(&brep, 0), "edge 0 should not be degenerate");
        assert!(is_degenerate_edge(&brep, 1), "edge 1 should be degenerate (zero length)");
        assert!(is_degenerate_edge(&brep, 2), "edge 2 should be degenerate (explicit flag)");
        assert!(!is_degenerate_edge(&brep, 99), "out-of-bounds should not be degenerate");
    }

    #[test]
    fn seam_edge_candidates_on_box() {
        // A plain box has no seam edges.
        let brep = box_2x2x2();
        assert!(
            seam_edge_candidates(&brep).is_empty(),
            "a box should have no seam edge candidates"
        );
    }
}
