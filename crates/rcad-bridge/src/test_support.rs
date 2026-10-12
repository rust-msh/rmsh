//! Test-only fixtures shared by the bridge modules (compiled under cfg(test)
//! only — never part of the published crate surface).

use std::sync::{Arc, RwLock};

use glam::DVec3;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, TShape};
use rcad_mesh::b_rep_mesh::facade::IncrementalMesh;

/// The triangle count the unit-box tessellation mounts (pinned after first
/// observing the pipeline output; the walk must stay stable against it).
pub(crate) const EXPECTED_BOX_TRIANGLES: usize = 12;

/// The deduped node count of the unit-box tessellation (the 8 corners, when
/// the deflection-0.1 mesh puts no interior or edge-midpoint nodes).
pub(crate) const EXPECTED_BOX_NODES: usize = 8;

/// A unit box [0,1]^3 BRep tessellated through the open-source
/// `IncrementalMesh` facade, returned shared because the mesher keeps its
/// own handle while running.
///
/// This is the documented tessellation entry point of the bridge card:
/// the walk under test only READS triangulations, the mounting happens here.
pub fn meshed_unit_box() -> Arc<RwLock<BRep>> {
    let brep = unit_box_brep();

    // MakeBox leaves the solid TShape in the pool without keeping the root
    // Shape handle; reconstruct the wrapper so the mesher's shape walk has a
    // root (the facade's own tests wrap shapes the same way).
    let solid_index = brep
        .tshapes
        .iter()
        .position(|t| matches!(t.as_ref(), TShape::Solid(_)))
        .expect("the box carries a solid TShape");
    let solid = Shape::from_parts(
        brep.tshapes[solid_index].clone(),
        solid_index,
        0,
        Orientation::Forward,
    );

    let shared: Arc<RwLock<BRep>> = Arc::new(RwLock::new(brep));
    // The `incmesh shape 0.1` DRAW mapping: linear deflection 0.1, not
    // relative, angular deflection 0.5 rad, single-threaded (the facade's
    // OCCT constructor defaults).
    let meshed =
        IncrementalMesh::perform_with_deflection(shared.clone(), &solid, 0.1, false, 0.5, false);
    assert!(meshed.is_done(), "IncrementalMesh::Perform is done");
    shared
}

/// A fresh (untessellated) unit box [0,1]^3 BRep for the geometry-source
/// queries: they project onto the face surfaces and need no triangulation.
pub fn unit_box_brep() -> BRep {
    rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("the unit box builds")
}
