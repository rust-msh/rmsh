/// Smoke tests for rcad-render's CPU tessellation path.
/// These do NOT require a GPU — they only test the mesh-building logic.
use rcad_kernel::BRep;
use rcad_kernel::PrimitiveSolid;
use rcad_render::{EditedModelDelta, TessellationOptions, Tessellator};

fn make_box_brep() -> BRep {
    BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    })
}

/// Tessellating an empty BRep should return an empty mesh without panicking.
#[test]
fn tessellate_empty_brep_no_panic() {
    let empty = BRep::default();
    let mesh = Tessellator::tessellate(&empty);
    assert!(mesh.vertices.is_empty(), "empty BRep should yield no vertices");
    assert!(mesh.indices.is_empty(), "empty BRep should yield no indices");
}

/// Tessellating a box primitive (no geometry populated) should produce vertices.
#[test]
fn tessellate_box_has_vertices() {
    let brep = make_box_brep();
    let mesh = Tessellator::tessellate(&brep);
    // The box primitive has 8 corner vertices
    assert_eq!(mesh.vertices.len(), 8, "unit box should have 8 vertices");
}

/// Tessellating a box with triangles populated should produce triangle indices.
#[test]
fn tessellate_box_with_triangles_has_indices() {
    use rcad_algorithms::geom_populate;

    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    });
    geom_populate::populate_box_geom(&mut brep);

    let mesh = Tessellator::tessellate(&brep);
    assert!(!mesh.vertices.is_empty(), "should have vertices");
    // Triangle indices come in sets of 3
    assert!(mesh.indices.len().is_multiple_of(3), "index count must be divisible by 3");
}

/// All triangle indices must be within bounds of the vertex buffer.
#[test]
fn tessellate_box_indices_in_bounds() {
    use rcad_algorithms::geom_populate;

    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    geom_populate::populate_box_geom(&mut brep);

    let mesh = Tessellator::tessellate(&brep);
    let nv = mesh.vertices.len() as u32;
    for &idx in &mesh.indices {
        assert!(idx < nv, "triangle index {idx} out of bounds (nv={nv})");
    }
}

#[test]
fn invalidate_cache_for_edge_edit_marks_adjacent_faces_dirty() {
    use rcad_algorithms::geom_populate;

    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    geom_populate::populate_box_geom(&mut brep);

    // Prime tessellation so all faces are clean before edit invalidation.
    let opts = TessellationOptions::default();
    let _ = Tessellator::tessellate_with_options(&mut brep, &opts);
    assert!(!brep.needs_remesh(), "primed brep should be mesh-clean");

    let edits = EditedModelDelta {
        modified_edges: vec![0],
        ..EditedModelDelta::default()
    };
    let marked = Tessellator::invalidate_cache_for_edits(&mut brep, &edits);

    assert!(marked > 0, "edge edit should invalidate at least one face");
    assert!(brep.needs_remesh(), "edge edit should mark some faces dirty");
}

#[test]
fn tessellate_after_edits_clears_dirty_faces_again() {
    use rcad_algorithms::geom_populate;

    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    geom_populate::populate_box_geom(&mut brep);

    let opts = TessellationOptions::default();
    let _ = Tessellator::tessellate_with_options(&mut brep, &opts);

    let edits = EditedModelDelta {
        modified_vertices: vec![0],
        ..EditedModelDelta::default()
    };
    let _ = Tessellator::tessellate_after_edits(&mut brep, &edits, &opts);

    assert!(
        !brep.needs_remesh(),
        "post-edit incremental tessellation should clear dirty flags"
    );
}

#[test]
fn invalidate_cache_for_edits_ignores_out_of_range_indices() {
    let mut brep = make_box_brep();
    let edits = EditedModelDelta {
        modified_vertices: vec![usize::MAX],
        modified_edges: vec![usize::MAX],
        modified_faces: vec![usize::MAX],
    };

    let marked = Tessellator::invalidate_cache_for_edits(&mut brep, &edits);
    assert_eq!(marked, 0, "out-of-range edits should be ignored safely");
}
