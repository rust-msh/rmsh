//! Compact tetrahedral mesh with pre-built neighbor tables.
//!
//! The node-pool form consumed by the volume optimizer passes
//! ([`crate::mesh_optimize`]: collapse / swap / 3-2 removal / relocation /
//! coloring) and the HXT parallel layer: nodes indexed by `u32` into a
//! flat `Vec`, tets carrying their four node indices and a face-neighbor
//! table. (The bistellar flip optimizer that used to live here was
//! retired when the 3-D quality card wired the delaunay_3d cluster into
//! the production refine loop — one implementation, the Mesh-side one,
//! remains.)

use std::collections::HashMap;

use rmsh_model::{ElementType, Mesh};

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

// ─── Constructors ────────────────────────────────────────────────────────────

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
        for (ti, _tet) in self.tets.iter().enumerate() {
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
        // [1,2,3,4], [1,2,4,5], [1,2,5,3]
        //
        // Corrected fixture (was [1,2,3,4],[1,2,4,5],[1,2,3,5] with nodes
        // 3=(0.5,0.001,0), 4=(0.5,0,0.001)): the old coordinates put the
        // crossing of segment (1,2) with the plane of triangle (3,4,5) at
        // (0.5, 0, 0), which lies OUTSIDE the triangle (its edge (3,4) is the
        // line y+z=0.001 in the y-z plane). The 3-2 flip on that fan is
        // bistellar-INVALID — it changed the union volume (lost 3.33e-7) —
        // and was only "passing" because the flip was unguarded. Nodes 3 and
        // 4 are now on the negative side of both y and z so the crossing
        // point is strictly inside triangle (3,4,5) and the flip is a valid
        // bistellar move (union volume exactly preserved). The invalid
        // variant lives on in the delaunay_3d flip tests
        // (local_flip_pass_skips_non_piercing_fan).
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.5, -0.002, -0.001));
        mesh.add_node(Node::new(4, 0.5, -0.001, -0.002));
        mesh.add_node(Node::new(5, 0.5, 0.8, 0.8));
        mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));
        mesh.add_element(Element::new(2, ElementType::Tetrahedron4, vec![1, 2, 4, 5]));
        mesh.add_element(Element::new(3, ElementType::Tetrahedron4, vec![1, 2, 5, 3]));
        mesh
    }

    #[test]
    fn neighbor_table_two_tets() {
        let mesh = make_two_tet_mesh();
        let tm = TetMesh::from_mesh(&mesh);

        // Find which face is shared. Each tet has 4 neighbors; exactly one
        // should be non-MAX and point to the other tet.
        let mut found = 0usize;
        for ti in 0..tm.tets.len() {
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
        assert_eq!(tm.tets.len(), 3);

        // All interior faces should have valid neighbors.
        let mut interior_faces = 0usize;
        let mut boundary_faces = 0usize;
        for ti in 0..tm.tets.len() {
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
}
