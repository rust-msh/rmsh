//! Delaunay 2-D — Bowyer-Watson Delaunay triangulation (Gmsh algorithm 5).
//!
//! # Algorithm overview
//!
//! Pure Delaunay triangular meshing via Bowyer-Watson point insertion.
//! Boundary edges are discretised and interior points seeded on a uniform
//! staggered grid, then the Delaunay triangulation is computed and
//! elements whose centroid falls outside the domain are discarded.
//!
//! This is the simplest 2-D mesher in the Gmsh taxonomy.  For domains
//! with holes the implementation delegates to
//! [`crate::planar_meshing::mesh_domain_triangles`].

use rmsh_model::Mesh;

use crate::planar_meshing::{mesh_domain_triangles, validate_domain};
use crate::traits::{Domain2D, MeshAlgoError, MeshParams, Mesher2D};
use crate::triangulate2d::{Polygon2D, mesh_polygon};

/// Delaunay 2-D mesher (Gmsh algorithm 5).
///
/// Produces triangular meshes via Bowyer-Watson Delaunay triangulation
/// with uniform staggered-grid point sampling.
#[derive(Debug, Default, Clone)]
pub struct Delaunay2D;

impl Delaunay2D {
    pub fn new() -> Self {
        Self
    }
}

impl Mesher2D for Delaunay2D {
    fn name(&self) -> &'static str {
        "Delaunay 2D"
    }

    fn mesh_2d(&self, domain: &Domain2D, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        validate_domain(domain, params.element_size)?;

        if domain.boundaries.len() == 1 {
            let poly = Polygon2D::new(domain.outer().to_vec());
            return mesh_polygon(&poly, params.element_size)
                .map_err(|e| MeshAlgoError::Generation(e.to_string()));
        }

        mesh_domain_triangles(domain, params.element_size, params.element_size * 0.866, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delaunay_2d_meshes_square() {
        let domain =
            Domain2D::from_outer(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);
        let mesh = Delaunay2D
            .mesh_2d(&domain, &MeshParams::with_size(0.5))
            .unwrap();
        assert!(mesh.element_count() > 0);
    }

    #[test]
    fn delaunay_2d_meshes_with_hole() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0]])
            .with_hole(vec![[1.0, 1.0], [2.0, 1.0], [2.0, 2.0], [1.0, 2.0]]);
        let mesh = Delaunay2D
            .mesh_2d(&domain, &MeshParams::with_size(0.5))
            .unwrap();
        assert!(mesh.element_count() > 0);
    }

    #[test]
    fn delaunay_2d_name_is_stable() {
        assert_eq!(Delaunay2D.name(), "Delaunay 2D");
    }
}
