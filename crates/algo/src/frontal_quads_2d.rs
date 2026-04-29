//! Frontal-Quads 2-D — frontal Delaunay mesh with quad recombination (Gmsh algorithm 8).
//!
//! # Algorithm overview
//!
//! Frontal-Quads generates a quadrilateral mesh in two stages:
//!
//! 1. **Frontal Delaunay triangulation**: generates a high-quality triangular mesh
//!    via [`crate::FrontalDelaunay2D`] (Gmsh algorithm 6).
//! 2. **Cross-field recombination**: computes a smooth 4-direction cross field on
//!    the triangle mesh, then greedily recombines adjacent triangle pairs into
//!    quadrilaterals following the cross-field alignment.
//!
//! The result is a quadrilateral-dominant mesh, though some isolated triangles
//! may remain in regions where recombination is not possible.
//!
//! # Reference
//!
//! Gmsh source: `Mesh/meshGFace.cpp`, algorithm 8 path.

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
}

impl Default for FrontalQuads2D {
    fn default() -> Self {
        Self {
            cross_field_iterations: 100,
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
        let quads = recombine_triangles(&nodes, &tris, &cf);

        // Stage 5: build output Mesh (quads + remaining triangles)
        let mut mesh = Mesh::new();
        let id_to_idx: std::collections::HashMap<u64, usize> = node_ids
            .iter()
            .enumerate()
            .map(|(i, &id)| (id, i))
            .collect();

        // Add nodes
        for (i, &pos) in nodes.iter().enumerate() {
            mesh.add_node(Node::new(i as u64 + 1, pos[0], pos[1], 0.0));
        }

        // Track which triangles were used in quads
        let mut tri_used = vec![false; tris.len()];

        for (elem_id, quad) in quads.iter().enumerate() {
            // Find which two triangles form this quad
            let mut found = Vec::new();
            for (ti, tri) in tris.iter().enumerate() {
                if tri_used[ti] {
                    continue;
                }
                let count = quad.iter().filter(|q| tri.contains(q)).count();
                if count >= 3 {
                    found.push(ti);
                }
            }
            if found.len() >= 2 {
                tri_used[found[0]] = true;
                tri_used[found[1]] = true;
                let nids: Vec<u64> = quad.iter().map(|&v| v as u64 + 1).collect();
                mesh.add_element(Element::new(
                    (elem_id + 1) as u64,
                    ElementType::Quad4,
                    nids,
                ));
            }
        }

        // Add remaining triangles as Triangle3 elements
        let mut next_elem_id = quads.len() as u64 + 1;
        for (ti, tri) in tris.iter().enumerate() {
            if tri_used[ti] {
                continue;
            }
            let nids = vec![tri[0] as u64 + 1, tri[1] as u64 + 1, tri[2] as u64 + 1];
            mesh.add_element(Element::new(next_elem_id, ElementType::Triangle3, nids));
            next_elem_id += 1;
        }

        if mesh.element_count() == 0 {
            return Ok(tri_mesh);
        }

        let _ = id_to_idx;
        Ok(mesh)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frontal_quads_rectangle() {
        let domain =
            Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 2.0], [0.0, 2.0]]);
        let params = MeshParams::with_size(0.8);
        let mesh = FrontalQuads2D::default()
            .mesh_2d(&domain, &params)
            .unwrap();
        assert!(mesh.element_count() > 0);
        // Should have at least some quad elements
        let quad_count = mesh
            .elements
            .iter()
            .filter(|e| e.etype == ElementType::Quad4)
            .count();
        assert!(quad_count > 0, "should produce some quads");
    }

    #[test]
    fn frontal_quads_square() {
        let domain =
            Domain2D::from_outer(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);
        let params = MeshParams::with_size(0.6);
        let mesh = FrontalQuads2D::default()
            .mesh_2d(&domain, &params)
            .unwrap();
        assert!(mesh.element_count() >= 2);
    }
}
