//! BRep → surface-mesh walk: reads every face TShape's active triangulation
//! into an [`rmsh_model::Mesh`] of Triangle3 elements.
//!
//! # Node-id policy (global position dedupe — documented adjudication)
//!
//! Faces are meshed independently, so every shared BRep edge carries one
//! copy of the same node coordinates inside each adjacent face's
//! triangulation. The walk deduplicates nodes GLOBALLY by quantized
//! position (1e-7 grid — the same quantization the rcad kernel uses for its
//! vertex identity cache `VertexKey`, topods.rs), producing a single mesh
//! node per geometric point. The resulting mesh is edge-manifold for closed
//! solids and therefore directly consumable by rmsh-algo's volume meshers
//! (the tetrahedralize3d boundary machinery keys connectivity on node ids)
//! and by the rmsh-io writers. (The pro-side reference duplicated every
//! triangle corner — 3 fresh nodes per triangle — which its CAE-export
//! consumers tolerated; rmsh's mesher-facing contract needs the fused form,
//! pinned by the tests here.)
//!
//! # Locations
//!
//! Triangulation nodes live in the face-local frame (OCCT BRep_TFace:
//! "the triangulation is in the same reference system than the TFace"); the
//! walk lifts them to world coordinates through the face location
//! (`TFaceData::surface_location`). For identity locations — every case the
//! pro reference covered — this is a no-op.
//!
//! # Difference from the pro reference (indexing fix — documented)
//!
//! `TFaceData::active_triangulation` is an index INTO
//! `TFaceData::triangulations` (rcad-kernel topods.rs: "index into
//! `triangulations`; None = null active handle"), whose entries are
//! `BRep::triangulations` pool indices. The pro reference used the active
//! index directly as the pool index — correct only when the face's active
//! triangulation happens to be mounted at the same pool position. The walk
//! here resolves the full two-level indirection, which stays correct for
//! faces carrying several triangulations or later-mounted ones.
//!
//! # Per-face physical tags (adjudicated: none in v1)
//!
//! `rmsh_model::Element` now carries `physical_tag`, and v1 leaves it
//! `None`: the BRep pool carries no face-to-group metadata a tag could be
//! derived from, and inventing a pool-index-based scheme here would pin a
//! policy with no consumer contract ahead of the P1 physical-groups
//! semantics. Revisit when a tag source exists.
//!
//! # Element / node id allocation
//!
//! Mesh node ids are allocated densely from 1 in first-encounter order
//! (deterministic for a fixed BRep pool order); element ids densely from 1
//! in walk order. Ids carry no OCCT meaning — they are mesh-local handles.

use std::collections::HashMap;

use glam::DVec3;
use rcad_kernel::topods::{BRep, TShape};
use rmsh_model::{Element, ElementType, Mesh, Node};

/// Quantization grid of the node dedupe (mirrors the kernel `VertexKey`
/// quantization constant).
const NODE_QUANT_STEP: f64 = 1.0e-7;

/// Quantized world position of a triangulation node (the dedupe key).
fn position_key(p: DVec3) -> (i64, i64, i64) {
    let q = |c: f64| (c / NODE_QUANT_STEP).round() as i64;
    (q(p.x), q(p.y), q(p.z))
}

/// Convert a BRep to a surface mesh (triangles).
///
/// Walks every face TShape in the pool and reads its active triangulation
/// (the `BRep::triangulations` pool — the OCCT `BRep_Tool::Triangulation`
/// handle stand-in). See the module docs for the node-id policy, the
/// location handling, and the deliberate departures from the pro-side
/// reference implementation.
pub fn brep_surface_mesh(brep: &BRep) -> Mesh {
    let mut mesh = Mesh::new();
    let mut next_node_id: u64 = 1;
    let mut next_elem_id: u64 = 1;
    // Quantized world position -> mesh node id (the global dedupe table).
    let mut node_by_pos: HashMap<(i64, i64, i64), u64> = HashMap::new();

    for ts in &brep.tshapes {
        let TShape::Face(fd) = ts.as_ref() else {
            continue;
        };
        // Two-level indirection (see the module docs): the active index
        // selects in the face's triangulation list, whose entries are
        // BRep pool indices. Faces without a mounted triangulation are
        // skipped (nothing to read).
        let Some(active) = fd.active_triangulation else {
            continue;
        };
        let Some(pool_pos) = fd.triangulations.get(active) else {
            continue;
        };
        let Some(tri) = brep.triangulations.get(*pool_pos) else {
            continue;
        };
        // The triangulation nodes live in the face-local frame; lift to
        // world coordinates (identity for the common unlocated case).
        let location = brep.get_location(fd.surface_location);
        for a_j in 1..=tri.nb_triangles() {
            // Poly_Triangle indices are 1-based (OCCT Poly_Triangle::Get).
            let (n1, n2, n3) = tri.triangle(a_j).get();
            let mut nids = Vec::with_capacity(3);
            for n in [n1, n2, n3] {
                let world = location.transform_point3(tri.node(n as usize));
                let key = position_key(world);
                let id = match node_by_pos.get(&key) {
                    Some(&id) => id,
                    None => {
                        let id = next_node_id;
                        next_node_id += 1;
                        node_by_pos.insert(key, id);
                        mesh.add_node(Node::new(id, world.x, world.y, world.z));
                        id
                    }
                };
                nids.push(id);
            }
            mesh.add_element(Element::new(next_elem_id, ElementType::Triangle3, nids));
            next_elem_id += 1;
        }
    }
    mesh
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use glam::DVec3;
    use rmsh_model::ElementType;

    use super::brep_surface_mesh;
    use crate::test_support::meshed_unit_box;

    /// The walk over the tessellated unit box: pin what the pipeline
    /// actually produces — triangle-only, dense ids from 1, every element
    /// node resolving, and edge-manifold connectivity (every undirected
    /// triangle edge shared by exactly two triangles — the watertight pin
    /// the node-dedupe policy exists for).
    #[test]
    fn box_walk_yields_manifold_triangle_mesh() {
        let brep = meshed_unit_box();
        let brep = brep.read().unwrap();
        let mesh = brep_surface_mesh(&brep);

        assert_eq!(
            mesh.elements.len(),
            crate::test_support::EXPECTED_BOX_TRIANGLES,
            "triangle count pinned to what IncrementalMesh mounts on the box"
        );
        assert_eq!(
            mesh.nodes.len(),
            crate::test_support::EXPECTED_BOX_NODES,
            "node count pinned to the deduped corner set (NOT 3 per triangle)"
        );

        let mut edge_faces: HashMap<[u64; 2], usize> = HashMap::new();
        for (i, elem) in mesh.elements.iter().enumerate() {
            assert_eq!(elem.etype, ElementType::Triangle3, "triangle-only output");
            assert_eq!(elem.id, (i + 1) as u64, "element ids dense from 1");
            assert_eq!(elem.node_ids.len(), 3);
            for nid in &elem.node_ids {
                assert!(
                    mesh.nodes.contains_key(nid),
                    "element {} references missing node {}",
                    elem.id,
                    nid
                );
            }
            for w in 0..3 {
                let mut key = [elem.node_ids[w], elem.node_ids[(w + 1) % 3]];
                key.sort_unstable();
                *edge_faces.entry(key).or_insert(0) += 1;
            }
        }
        let bad: Vec<_> = edge_faces.iter().filter(|(_, c)| **c != 2).collect();
        assert!(
            bad.is_empty(),
            "edge-manifold (every edge exactly 2 triangles); offenders: {:?}",
            &bad[..bad.len().min(8)]
        );
        let node_ids: Vec<u64> = {
            let mut ids: Vec<u64> = mesh.nodes.keys().copied().collect();
            ids.sort_unstable();
            ids
        };
        assert_eq!(
            node_ids,
            (1..=mesh.nodes.len() as u64).collect::<Vec<_>>(),
            "node ids dense from 1 (first-encounter allocation)"
        );
    }

    /// The dedupe policy differentiator against the pro reference: shared
    /// corners are single nodes. For the unit box the 8 corners must each
    /// appear once regardless of how many triangle corners touch them.
    #[test]
    fn box_corners_are_single_nodes() {
        let brep = meshed_unit_box();
        let brep = brep.read().unwrap();
        let mesh = brep_surface_mesh(&brep);

        let corners: [DVec3; 8] = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(1.0, 0.0, 1.0),
            DVec3::new(0.0, 1.0, 1.0),
            DVec3::new(1.0, 1.0, 1.0),
        ];
        for c in corners {
            let count = mesh
                .nodes
                .values()
                .filter(|&n| (DVec3::new(n.position.x, n.position.y, n.position.z) - c).length() < 1e-9)
                .count();
            assert_eq!(count, 1, "corner {:?} must appear exactly once", c);
        }
    }

    /// An empty BRep yields an empty mesh (legal, not an error).
    #[test]
    fn empty_brep_yields_empty_mesh() {
        let brep = rcad_kernel::BRep::new();
        let mesh = brep_surface_mesh(&brep);
        assert!(mesh.nodes.is_empty());
        assert!(mesh.elements.is_empty());
    }

    /// Downstream consumability pin (the bridge-adoption behavioral
    /// contract): the walked surface mesh feeds rmsh-algo's closed-surface
    /// tetrahedralizer directly.
    #[test]
    fn walked_box_mesh_feeds_the_tetrahedralizer() {
        let brep = meshed_unit_box();
        let brep = brep.read().unwrap();
        let mesh = brep_surface_mesh(&brep);
        let volume = rmsh_algo::tetrahedralize_closed_surface(&mesh)
            .expect("the walked box mesh tetrahedralizes");
        assert!(
            volume.elements.iter().any(|e| e.dimension() == 3),
            "the volume mesh carries tetrahedra"
        );
    }

    /// Downstream consumability pin: the walked mesh round-trips through
    /// the rmsh-io MSH writer.
    #[test]
    fn walked_box_mesh_round_trips_through_the_msh_writer() {
        let brep = meshed_unit_box();
        let brep = brep.read().unwrap();
        let mesh = brep_surface_mesh(&brep);

        let mut bytes = Vec::new();
        rmsh_io::write_msh_v2(&mut bytes, &mesh).expect("the mesh writes as MSH v2");
        let parsed = rmsh_io::load_msh_from_bytes(&bytes).expect("the written file parses");
        assert_eq!(parsed.nodes.len(), mesh.nodes.len());
        assert_eq!(parsed.elements.len(), mesh.elements.len());
    }
}
