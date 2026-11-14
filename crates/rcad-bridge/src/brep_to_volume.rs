//! The BRep -> volume mesh bridge.
//!
//! This is the composition of the two halves the charter already owns:
//!
//! 1. `rmsh_model::brep_surface_mesh` (re-exported from `rcad-mesh`) walks a
//!    BRep's mounted face triangulations into a closed, edge-manifold
//!    triangle [`Mesh`].
//! 2. [`rmsh_algo::Delaunay3D`] (the production 3-D mesher) fills that shell
//!    with tetrahedra.
//!
//! The bridge adds the gates between them: a BRep with no triangulation, an
//! empty boundary, or a boundary that is not a closed manifold shell is
//! reported as a typed error rather than silently producing a partial volume
//! mesh. Tessellation itself stays the caller's step (the `rcad-tess`
//! `IncrementalMesh` facade mounts the triangulations); the
//! [`BridgeError::MissingTriangulation`] message names it.
//!
//! [`inspect_volume_mesh`] is the no-solver validation: total volume and the
//! tetrahedral quality statistics, so a bridge consumer can assert volume
//! conservation and quality without a solver in the loop.

use std::collections::HashMap;

use rcad_kernel::topods::{BRep, TShape};
use rmsh_algo::{Delaunay3D, MeshParams, Mesher3D, tet_quality_stats, TetQualityStats};
use rmsh_model::{brep_surface_mesh, Element, ElementType, Mesh};
use thiserror::Error;

/// A bridge failure, with the reason a consumer needs to act on.
#[derive(Debug, Error, PartialEq)]
pub enum BridgeError {
    /// The BRep has no face TShapes at all.
    #[error("the B-rep carries no faces")]
    EmptyBoundary,
    /// The BRep has faces but none carries a triangulation, so there is no
    /// boundary to fill.
    #[error("the B-rep faces carry no triangulation; mount one first with the rcad-tess IncrementalMesh facade")]
    MissingTriangulation,
    /// The walked boundary is not a closed manifold triangle shell.
    #[error("the B-rep boundary is not a closed manifold shell: {0}")]
    NotClosed(String),
    /// The volume mesher refused the surface.
    #[error("volume meshing failed: {0}")]
    Meshing(String),
}

/// The no-solver validation of a volume mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct VolumeMeshReport {
    /// The number of mesh nodes.
    pub node_count: usize,
    /// The number of tetrahedra.
    pub tet_count: usize,
    /// The total enclosed volume (the sum of the absolute tetra volumes, so
    /// an inconsistent tet winding does not cancel it away).
    pub volume: f64,
    /// The minimum dihedral angle over the tetrahedra, in degrees.
    pub min_dihedral_deg: f64,
    /// The 95th percentile of the circumradius-over-shortest-edge ratio.
    pub p95_radius_edge: f64,
    /// The fraction of sliver-like tetrahedra.
    pub sliver_fraction: f64,
}

impl VolumeMeshReport {
    /// The relative volume error against an analytically known volume.
    pub fn volume_error(&self, expected: f64) -> f64 {
        if expected.abs() <= f64::EPSILON {
            return 0.0;
        }
        (self.volume - expected).abs() / expected
    }

    /// The underlying tetrahedral quality statistics.
    pub fn quality(&self) -> TetQualityStats {
        TetQualityStats {
            min_dihedral_deg: self.min_dihedral_deg,
            p95_radius_edge: self.p95_radius_edge,
            sliver_fraction: self.sliver_fraction,
            tet_count: self.tet_count,
        }
    }
}

/// Fill a BRep's interior with a tetrahedral volume mesh.
///
/// The BRep must already carry face triangulations (see the module docs);
/// the walk reads them and the mesher fills the closed shell they form.
pub fn brep_to_volume_mesh(brep: &BRep, params: &MeshParams) -> Result<Mesh, BridgeError> {
    let face_count = brep
        .tshapes
        .iter()
        .filter(|t| matches!(t.as_ref(), TShape::Face(_)))
        .count();
    if face_count == 0 {
        return Err(BridgeError::EmptyBoundary);
    }

    let surface = brep_surface_mesh(brep);
    let triangle_count = surface
        .elements
        .iter()
        .filter(|e| e.etype == ElementType::Triangle3 && e.node_ids.len() == 3)
        .count();
    if triangle_count == 0 {
        return Err(BridgeError::MissingTriangulation);
    }
    check_closed_shell(&surface).map_err(BridgeError::NotClosed)?;

    Delaunay3D::new()
        .mesh_3d(&surface, params)
        .map_err(|e| BridgeError::Meshing(e.to_string()))
}

/// Validate a volume mesh without a solver: volume plus tet quality.
pub fn inspect_volume_mesh(mesh: &Mesh) -> VolumeMeshReport {
    let quality = tet_quality_stats(mesh);
    VolumeMeshReport {
        node_count: mesh.nodes.len(),
        tet_count: quality.tet_count,
        volume: mesh_volume(mesh),
        min_dihedral_deg: quality.min_dihedral_deg,
        p95_radius_edge: quality.p95_radius_edge,
        sliver_fraction: quality.sliver_fraction,
    }
}

/// The boundary must be a closed manifold triangle shell: every undirected
/// edge shared by exactly two triangles.
fn check_closed_shell(mesh: &Mesh) -> Result<(), String> {
    let mut edges: HashMap<(u64, u64), u32> = HashMap::new();
    let mut triangles = 0usize;
    for elem in &mesh.elements {
        if elem.etype != ElementType::Triangle3 || elem.node_ids.len() != 3 {
            continue;
        }
        triangles += 1;
        let ids = elem.node_ids.as_slice();
        for k in 0..3 {
            let a = ids[k];
            let b = ids[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            *edges.entry(key).or_insert(0) += 1;
        }
    }
    if triangles == 0 {
        return Err("no triangles".to_string());
    }
    let boundary = edges.values().filter(|c| **c == 1).count();
    let non_manifold = edges.values().filter(|c| **c > 2).count();
    if boundary > 0 || non_manifold > 0 {
        return Err(format!(
            "{boundary} boundary edges and {non_manifold} non-manifold edges over {triangles} triangles"
        ));
    }
    Ok(())
}

/// The sum of the absolute signed tetra volumes.
fn mesh_volume(mesh: &Mesh) -> f64 {
    let mut sum = 0.0_f64;
    for elem in &mesh.elements {
        if elem.etype != ElementType::Tetrahedron4 || elem.node_ids.len() != 4 {
            continue;
        }
        let p = |id: &u64| -> [f64; 3] {
            let n = &mesh.nodes[id];
            [n.position.x, n.position.y, n.position.z]
        };
        let (a, b, c, d) = (
            p(&elem.node_ids[0]),
            p(&elem.node_ids[1]),
            p(&elem.node_ids[2]),
            p(&elem.node_ids[3]),
        );
        let b = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let c = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let d = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
        let cross = [
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        ];
        sum += (d[0] * cross[0] + d[1] * cross[1] + d[2] * cross[2]).abs() / 6.0;
    }
    sum
}

/// Collect the tetrahedral elements of a volume mesh (the FEA-facing slice).
pub fn tetrahedra(mesh: &Mesh) -> Vec<[u64; 4]> {
    mesh.elements
        .iter()
        .filter(|e: &&Element| e.etype == ElementType::Tetrahedron4 && e.node_ids.len() == 4)
        .map(|e| [e.node_ids[0], e.node_ids[1], e.node_ids[2], e.node_ids[3]])
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::meshed_unit_box;
    use glam::DVec3;
    use rcad_kernel::topods::{Orientation, Shape};
    use rcad_modeling::make_box_brep;
    use std::sync::{Arc, RwLock};

    fn unit_box_params() -> MeshParams {
        MeshParams::with_size(0.4)
    }

    /// Mount a triangulation on every top-level solid / face of `brep`,
    /// exactly as a caller would before bridging.
    fn tessellate(brep: &BRep) -> Arc<RwLock<BRep>> {
        let shared = Arc::new(RwLock::new(brep.clone()));
        let roots: Vec<Shape> = {
            let guard = shared.read().unwrap();
            guard
                .tshapes
                .iter()
                .enumerate()
                .filter(|(_, t)| matches!(t.as_ref(), TShape::Solid(_) | TShape::Face(_)))
                .map(|(i, t)| Shape::from_parts(t.clone(), i, 0, Orientation::Forward))
                .collect()
        };
        for root in roots {
            let mesher = rcad_tess::b_rep_mesh::facade::IncrementalMesh::perform_with_deflection(
                Arc::clone(&shared),
                &root,
                0.1,
                false,
                0.5,
                false,
            );
            assert!(mesher.is_done(), "IncrementalMesh::Perform is done");
        }
        shared
    }

    #[test]
    fn a_closed_box_fills_with_tetrahedra() {
        let brep = meshed_unit_box();
        let brep = brep.read().unwrap();
        let mesh = brep_to_volume_mesh(&brep, &unit_box_params()).expect("the box bridges");

        let report = inspect_volume_mesh(&mesh);
        assert!(report.tet_count > 0, "the box must tetrahedralize");
        assert!(!mesh.elements.is_empty());
        assert!(tetrahedra(&mesh).len() == report.tet_count);
        // The unit box has volume 1: the mesh must conserve it.
        assert!(
            report.volume_error(1.0) < 0.02,
            "volume {} vs 1 ({}% off)",
            report.volume,
            report.volume_error(1.0) * 100.0
        );
        assert!(report.min_dihedral_deg > 0.0);
    }

    #[test]
    fn a_larger_box_scales_its_volume() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).unwrap();
        let shared = tessellate(&brep);
        let brep = shared.read().unwrap();
        let mesh =
            brep_to_volume_mesh(&brep, &MeshParams::with_size(0.6)).expect("the box bridges");
        let report = inspect_volume_mesh(&mesh);
        assert!(report.volume_error(24.0) < 0.02, "volume {}", report.volume);
    }

    #[test]
    fn an_untessellated_brep_is_reported_not_guessed() {
        let brep = crate::test_support::unit_box_brep();
        let err = brep_to_volume_mesh(&brep, &unit_box_params()).unwrap_err();
        assert_eq!(err, BridgeError::MissingTriangulation);
    }

    #[test]
    fn a_brep_without_faces_is_reported() {
        let brep = BRep::new();
        let err = brep_to_volume_mesh(&brep, &unit_box_params()).unwrap_err();
        assert_eq!(err, BridgeError::EmptyBoundary);
    }

    #[test]
    fn an_open_shell_is_reported() {
        // A box with one face left untessellated: the walk skips that face,
        // so the boundary has a rectangular hole and the closed-shell gate
        // must reject it.
        let brep = crate::test_support::unit_box_brep();
        let shared = Arc::new(RwLock::new(brep));
        let faces: Vec<Shape> = {
            let guard = shared.read().unwrap();
            guard
                .tshapes
                .iter()
                .enumerate()
                .filter(|(_, t)| matches!(t.as_ref(), TShape::Face(_)))
                .map(|(i, t)| Shape::from_parts(t.clone(), i, 0, Orientation::Forward))
                .collect()
        };
        assert_eq!(faces.len(), 6, "the box has six faces");
        for root in faces.iter().skip(1) {
            let mesher = rcad_tess::b_rep_mesh::facade::IncrementalMesh::perform_with_deflection(
                Arc::clone(&shared),
                root,
                0.1,
                false,
                0.5,
                false,
            );
            assert!(mesher.is_done());
        }
        let brep = shared.read().unwrap();
        let err = brep_to_volume_mesh(&brep, &unit_box_params()).unwrap_err();
        assert!(
            matches!(err, BridgeError::NotClosed(_)),
            "expected NotClosed, got {err:?}"
        );
    }
}
