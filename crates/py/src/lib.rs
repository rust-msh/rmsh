//! The `rmsh` Python extension module (`rmsh._rmsh`).
//!
//! A gmsh-compatible procedural API over the current rmsh + rcad reality:
//! CAD construction binds the real rcad kernel (primitives, the
//! OCCT-aligned boolean pipeline, STEP I/O), meshing binds the rmsh algo
//! crates through the rmsh-bridge B-Rep walker, and the msh I/O binds
//! rmsh-io with the P5 semantic element variants feeding the gmsh type-id
//! data layout natively.
//!
//! Deliberate omissions (charter: no GUI/viewer in the rmsh product):
//! `gui.*` and `plugin.*` are gone; `model.geo.*` (gmsh's built-in CAD
//! kernel namespace) is not part of v1 - build geometry through
//! `model.occ.*`.

mod io_api;
mod mesh_api;
mod model_api;
mod occ;
mod state;

use pyo3::prelude::*;

// The function name defines the exported PyInit__rmsh symbol that the
// maturin module-name `rmsh._rmsh` requires.
#[pymodule]
fn _rmsh(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    // root
    m.add_function(wrap_pyfunction!(io_api::initialize_impl, m)?)?;
    m.add_function(wrap_pyfunction!(io_api::finalize_impl, m)?)?;
    m.add_function(wrap_pyfunction!(io_api::clear_impl, m)?)?;
    m.add_function(wrap_pyfunction!(io_api::open_impl, m)?)?;
    m.add_function(wrap_pyfunction!(io_api::merge_impl, m)?)?;
    m.add_function(wrap_pyfunction!(io_api::write_impl, m)?)?;

    // model
    m.add_function(wrap_pyfunction!(model_api::model_add, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_remove, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_get_current, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_set_current, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_get_dimension, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_get_entities, m)?)?;
    m.add_function(wrap_pyfunction!(
        model_api::model_get_entities_in_bounding_box,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(model_api::model_get_bounding_box, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_get_entity_name, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_set_entity_name, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_add_physical_group, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_get_physical_groups, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_set_physical_name, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::model_get_physical_name, m)?)?;

    // model.occ
    m.add_function(wrap_pyfunction!(occ::model_occ_add_box, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_sphere, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_cylinder, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_cone, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_torus, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_rectangle, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_disk, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_point, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_line, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_circle_arc, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_spline, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_curve_loop, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_plane_surface, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_surface_loop, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_add_volume, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_cut, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_fuse, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_fragment, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_intersect, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_copy, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_remove, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_translate, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_rotate, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_dilate, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_mirror, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_import_shapes, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_synchronize, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_get_mass, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_get_properties, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_extrude, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_revolve, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_fillet, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_chamfer, m)?)?;
    m.add_function(wrap_pyfunction!(occ::model_occ_heal_shapes, m)?)?;

    // model.mesh
    m.add_function(wrap_pyfunction!(mesh_api::model_mesh_set_size, m)?)?;
    m.add_function(wrap_pyfunction!(mesh_api::model_mesh_generate, m)?)?;
    m.add_function(wrap_pyfunction!(mesh_api::model_mesh_set_order, m)?)?;
    m.add_function(wrap_pyfunction!(mesh_api::model_mesh_get_nodes, m)?)?;
    m.add_function(wrap_pyfunction!(mesh_api::model_mesh_get_elements, m)?)?;
    m.add_function(wrap_pyfunction!(mesh_api::model_mesh_clear, m)?)?;
    m.add_function(wrap_pyfunction!(mesh_api::model_mesh_optimize, m)?)?;
    m.add_function(wrap_pyfunction!(mesh_api::model_mesh_refine, m)?)?;
    m.add_function(wrap_pyfunction!(mesh_api::model_mesh_recombine, m)?)?;

    // option
    m.add_function(wrap_pyfunction!(model_api::option_set_number, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::option_get_number, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::option_set_string, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::option_get_string, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::option_set_color, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::option_get_color, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::option_restore_defaults, m)?)?;

    // logger
    m.add_function(wrap_pyfunction!(model_api::logger_start, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::logger_stop, m)?)?;
    m.add_function(wrap_pyfunction!(model_api::logger_get, m)?)?;

    Ok(())
}
