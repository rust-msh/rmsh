//! Downstream consumability pins for the shared BRep -> Mesh walk
//! (`rcad_mesh::brep_surface_mesh`): the walked surface mesh feeds rmsh's
//! closed-surface tetrahedralizer directly and round-trips through the
//! rmsh-io MSH writer. These tests live on the rmsh side because they
//! consume rmsh crates; the walk itself stays rcad-side.

use rcad_mesh::brep_surface_mesh;

use crate::test_support::meshed_unit_box;

/// The walked surface mesh feeds rmsh-algo's closed-surface
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

/// The walked mesh round-trips through the rmsh-io MSH writer.
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
