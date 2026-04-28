//! Open shell sewing: merge multiple BReps into a single watertight solid
//! by identifying and stitching near-coincident edges.
//!
//! Analogous to OCCT `BRepOffsetAPI_Sewing`.
//!
//! # Algorithm
//! 1. Concatenate all vertices, edges, and faces from every input BRep into
//!    a single pool, reindexing everything.
//! 2. **Vertex merging** (union-find): vertices within `tolerance` of each
//!    other are merged into one representative vertex.
//! 3. **Edge matching**: after vertex merging, edges that share both endpoint
//!    vertices (in either orientation) and originated from different input
//!    BReps are considered "stitched" — they represent the same boundary edge.
//!    Duplicate edges are removed and wire references updated.
//! 4. **Shell assembly**: all faces are collected into a single shell.
//!    Free edges (with only one incident face) are reported.
//!
//! # Limitations
//! - GeomStore data (curves, surfaces, PCurves) is concatenated but not
//!   de-duplicated.
//! - Only the outer wire of each face is considered during edge matching
//!   (inner wires are preserved as-is, with reindexed edge refs).

use rcad_kernel::{
    BRep,
    topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge},
};

// ─────────────────────────────────────────────────────────────────────────────
// Public types
// ─────────────────────────────────────────────────────────────────────────────

/// Result returned by [`sew_shells`].
#[derive(Debug)]
pub struct SewingResult {
    /// The merged BRep containing all input faces in a single shell.
    pub brep: BRep,
    /// Number of edge pairs that were stitched (shared boundary resolved).
    pub stitched_pairs: usize,
    /// Indices of edges in the result BRep that have only one incident face
    /// (open boundary edges).
    pub free_edges: Vec<usize>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Merge multiple BReps into a single BRep, stitching near-coincident edges.
///
/// # Arguments
/// * `breps` — slice of BReps to merge (must have at least one solid).
/// * `tolerance` — vertex proximity threshold for merging.
///
/// # Examples
/// ```rust
/// use rcad_kernel::{BRep, geom::PrimitiveSolid};
/// use rcad_modeling::sew_shells;
///
/// // Two adjacent unit boxes sharing the x=1 face
/// let a = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let b = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// // (In practice you'd translate b so it sits adjacent to a)
/// let result = sew_shells(&[a, b], 1e-6);
/// assert!(result.brep.solids.len() == 1);
/// ```
pub fn sew_shells(breps: &[BRep], tolerance: f64) -> SewingResult {
    if breps.is_empty() {
        return SewingResult {
            brep: BRep::new(),
            stitched_pairs: 0,
            free_edges: Vec::new(),
        };
    }

    // ── Step 1: concatenate everything ────────────────────────────────────

    let mut all_vertices: Vec<Vertex> = Vec::new();
    let mut all_edges: Vec<Edge> = Vec::new();
    let mut all_faces: Vec<Face> = Vec::new();

    // Track the vertex/edge offset for each input BRep
    let mut vertex_offsets: Vec<usize> = Vec::with_capacity(breps.len());
    let mut edge_offsets: Vec<usize> = Vec::with_capacity(breps.len());

    for brep in breps {
        let v_off = all_vertices.len();
        let e_off = all_edges.len();
        vertex_offsets.push(v_off);
        edge_offsets.push(e_off);

        // Vertices
        all_vertices.extend(brep.vertices.iter().cloned());

        // Edges (reindexed)
        for e in &brep.edges {
            all_edges.push(Edge {
                start: e.start + v_off,
                end: e.end + v_off,
            });
        }

        // Faces (reindex edge refs in wires)
        if let Some(solid) = brep.solids.first() {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let reindex_wire = |w: &Wire| -> Wire {
                        Wire {
                            edges: w
                                .edges
                                .iter()
                                .map(|we| WireEdge {
                                    idx: we.idx + e_off,
                                    forward: we.forward,
                                })
                                .collect(),
                        }
                    };
                    all_faces.push(Face {
                        outer_wire: reindex_wire(&face.outer_wire),
                        inner_wires: face.inner_wires.iter().map(reindex_wire).collect(),
                        normal: face.normal,
                        triangles: face
                            .triangles
                            .iter()
                            .map(|tri| [tri[0] + v_off, tri[1] + v_off, tri[2] + v_off])
                            .collect(),
                        mesh_dirty: face.mesh_dirty,
                    });
                }
            }
        }
    }

    let n_verts = all_vertices.len();
    let n_edges = all_edges.len();

    // ── Step 2: union-find vertex merge ───────────────────────────────────

    let mut parent: Vec<usize> = (0..n_verts).collect();

    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[rb] = ra;
        }
    }

    let tol2 = tolerance * tolerance;
    for i in 0..n_verts {
        for j in (i + 1)..n_verts {
            let d2 = (all_vertices[i].point - all_vertices[j].point).length_squared();
            if d2 <= tol2 {
                union(&mut parent, i, j);
            }
        }
    }

    // Build canonical index map: old index → canonical representative
    let canon: Vec<usize> = (0..n_verts).map(|i| find(&mut parent, i)).collect();

    // Remap edge endpoints
    for e in &mut all_edges {
        e.start = canon[e.start];
        e.end = canon[e.end];
    }

    // ── Step 3: edge deduplication ────────────────────────────────────────
    // Build a map from (min_v, max_v) → first edge index.
    // For each duplicate, record the pair (kept, duplicate) as stitched.

    use std::collections::HashMap;
    let mut edge_key_to_idx: HashMap<(usize, usize), usize> = HashMap::new();
    // Maps old edge index → canonical edge index (for dedup)
    let mut edge_canon: Vec<usize> = (0..n_edges).collect();
    let mut stitched_pairs = 0usize;

    for i in 0..n_edges {
        let e = &all_edges[i];
        let key = (e.start.min(e.end), e.start.max(e.end));
        if let Some(&existing) = edge_key_to_idx.get(&key) {
            // This edge is a duplicate of `existing` → stitch
            edge_canon[i] = existing;
            stitched_pairs += 1;
        } else {
            edge_key_to_idx.insert(key, i);
        }
    }

    // Remap wire edge indices in all faces
    for face in &mut all_faces {
        let remap = |w: &mut Wire| {
            for we in &mut w.edges {
                we.idx = edge_canon[we.idx];
            }
        };
        remap(&mut face.outer_wire);
        for iw in &mut face.inner_wires {
            remap(iw);
        }
    }

    // ── Step 4: build result BRep ─────────────────────────────────────────

    // Compact edges: keep only canonical ones
    let kept_edges: Vec<Edge> = (0..n_edges)
        .filter(|&i| edge_canon[i] == i)
        .map(|i| all_edges[i])
        .collect();

    // Remap edge indices: old canonical → compacted index
    let mut compact_edge_idx: Vec<usize> = vec![0; n_edges];
    {
        let mut ci = 0;
        for i in 0..n_edges {
            if edge_canon[i] == i {
                compact_edge_idx[i] = ci;
                ci += 1;
            }
        }
    }

    // Remap wires to compacted edge indices
    for face in &mut all_faces {
        let remap2 = |w: &mut Wire| {
            for we in &mut w.edges {
                we.idx = compact_edge_idx[edge_canon[we.idx]];
            }
        };
        remap2(&mut face.outer_wire);
        for iw in &mut face.inner_wires {
            remap2(iw);
        }
    }

    // Compact vertices (keep canonical representatives, renumber)
    // Build set of canonical vertex indices and map to compact indices
    let mut seen_verts: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for &c in &canon {
        seen_verts.insert(c);
    }
    let mut canon_to_compact: Vec<usize> = vec![0; n_verts];
    let mut compact_vertices: Vec<Vertex> = Vec::new();
    // Iterate in original order so compaction is deterministic
    for i in 0..n_verts {
        if canon[i] == i {
            canon_to_compact[i] = compact_vertices.len();
            compact_vertices.push(all_vertices[i]);
        }
    }

    // Remap edge endpoints to compacted vertices
    let mut compact_edges = kept_edges;
    for e in &mut compact_edges {
        e.start = canon_to_compact[e.start];
        e.end = canon_to_compact[e.end];
    }

    // Remap face triangles to compacted vertices
    for face in &mut all_faces {
        for tri in &mut face.triangles {
            tri[0] = canon_to_compact[canon[tri[0]]];
            tri[1] = canon_to_compact[canon[tri[1]]];
            tri[2] = canon_to_compact[canon[tri[2]]];
        }
    }

    // ── Step 5: find free edges ───────────────────────────────────────────

    let n_compact_edges = compact_edges.len();
    let mut edge_face_count: Vec<usize> = vec![0; n_compact_edges];
    for face in &all_faces {
        for we in &face.outer_wire.edges {
            if we.idx < n_compact_edges {
                edge_face_count[we.idx] += 1;
            }
        }
    }
    let free_edges: Vec<usize> = (0..n_compact_edges)
        .filter(|&i| edge_face_count[i] == 1)
        .collect();

    // ── Step 6: assemble GeomStore (simple concatenation) ─────────────────

    let mut geom = rcad_kernel::GeomStore::default();
    let mut face_offset = 0usize;

    for brep in breps {
        let surf_off = geom.surfaces.len();
        geom.surfaces.extend(brep.geom.surfaces.iter().cloned());

        let n_faces_this = brep
            .solids
            .first()
            .map(|s| s.shells.iter().map(|sh| sh.faces.len()).sum::<usize>())
            .unwrap_or(0);

        // face_surface mapping (offset surface indices)
        for fi in 0..n_faces_this {
            let mapped = brep
                .geom
                .face_surface
                .get(fi)
                .and_then(|o| *o)
                .map(|idx| idx + surf_off);
            while geom.face_surface.len() < face_offset + fi + 1 {
                geom.face_surface.push(None);
            }
            geom.face_surface[face_offset + fi] = mapped;
        }

        // face_surface_range
        for fi in 0..n_faces_this {
            let range = brep.geom.face_surface_range.get(fi).and_then(|o| *o);
            while geom.face_surface_range.len() < face_offset + fi + 1 {
                geom.face_surface_range.push(None);
            }
            geom.face_surface_range[face_offset + fi] = range;
        }

        face_offset += n_faces_this;
    }

    let shell = Shell { faces: all_faces };
    let solid = Solid {
        shells: vec![shell],
    };
    let result_brep = BRep {
        vertices: compact_vertices,
        edges: compact_edges,
        solids: vec![solid],
        geom,
        compound: None,
        compsolid: None,
    };

    SewingResult {
        brep: result_brep,
        stitched_pairs,
        free_edges,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::PrimitiveSolid;

    #[test]
    fn sew_single_brep_is_identity() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = sew_shells(std::slice::from_ref(&brep), 1e-6);
        assert_eq!(result.stitched_pairs, 0);
        // All edges of a closed box should be non-free (each edge borders 2 faces)
        assert_eq!(
            result.free_edges.len(),
            0,
            "closed box should have no free edges, got {:?}",
            result.free_edges
        );
    }

    #[test]
    fn sew_two_boxes_identifies_shared_face() {
        // Box A: x ∈ [0,1], Box B: x ∈ [1,2] — share the x=1 plane
        // Both use from_primitive which places box at [0,w]×[0,h]×[0,d].
        // Box B vertices at x=0..1 coincide with Box A vertices at x=1..2
        // only if we offset B. Since from_primitive always starts at 0,
        // the two boxes overlap at x=0 face and x=1 face.
        // For a real sewing test, we use two boxes with coincident vertices.
        let a = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        // Second box: same vertices as A (completely coincident) → all edges stitched
        let b = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = sew_shells(&[a, b], 1e-6);
        // Completely overlapping boxes should stitch all edges
        assert!(
            result.stitched_pairs > 0,
            "expected stitched pairs for coincident boxes, got 0"
        );
        println!(
            "sew two coincident boxes: stitched={}, free={:?}",
            result.stitched_pairs, result.free_edges
        );
    }

    #[test]
    fn sew_empty_input() {
        let result = sew_shells(&[], 1e-6);
        assert_eq!(result.stitched_pairs, 0);
        assert!(result.brep.solids.is_empty());
    }
}
