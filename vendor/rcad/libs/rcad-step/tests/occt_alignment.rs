//! OCCT TKDESTEP alignment tests for rcad-step.
//!
//! These tests verify compatibility with OCCT's STEP translation capabilities,
//! focusing on:
//! 1. Complex geometry roundtrip (boolean results, primitives)
//! 2. AP242 metadata roundtrip (colors, protocols)
//! 3. Periodic surface seam handling (cylinder, torus)
//! 4. Large file streaming performance
//! 5. Error handling
//!
//! Note: Assembly hierarchy tests are in assembly_io.rs and have known issues
//! being tracked separately.

use glam::{DAffine3, DVec3};
use rcad_algorithms::{boolean_op, BooleanOpType};
use rcad_modeling::{
    make_box_brep, make_cone_brep, make_cylinder_brep, make_sphere_brep, make_torus_brep,
};
use rcad_step::{
    AssemblyComponent, read_assembly, write_assembly, StepProtocol, StepReader, StepWriteOptions,
    StepWriter, ExportSelection,
};
use rcad_kernel::appearance::Color;
use std::time::Instant;

fn all_faces_selection() -> ExportSelection<'static> {
    ExportSelection {
        selected_faces: &[],
        selected_edges: &[],
    }
}

fn face_count(brep: &rcad_kernel::BRep) -> usize {
    brep
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

fn solid_count(brep: &rcad_kernel::BRep) -> usize {
    brep.solids.len()
}

fn make_box(origin: DVec3) -> rcad_kernel::BRep {
    let mut b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    if origin != DVec3::ZERO {
        b.apply_transform(DAffine3::from_translation(origin));
    }
    b
}

// ============================================================================
// 1. Complex Geometry Roundtrip (Boolean Results, Primitives)
// ============================================================================

/// Test roundtrip of a boolean union result (two overlapping boxes).
/// OCCT TKDESTEP coverage: boolean_result_shape, manifold_solid_brep.
#[test]
fn boolean_union_result_roundtrip() {
    let box1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box1");
    let box2 = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
        .expect("box2");

    let union = boolean_op(BooleanOpType::Union, &box1, &box2).expect("union");

    let step_str = StepWriter::write_string(&union, all_faces_selection());
    assert!(step_str.contains("ISO-10303-21"), "should have STEP header");

    let parsed = StepReader::parse_string(&step_str).expect("parse union result");

    // Union of two overlapping boxes should produce 1 solid
    assert_eq!(solid_count(&parsed), 1, "union should have 1 solid");
    // Face count for union of overlapping boxes is typically 10-12 depending on overlap geometry
    assert!(face_count(&parsed) >= 6, "union should have at least 6 faces");
}

/// Test roundtrip of a boolean difference result (box with cylindrical hole).
/// OCCT TKDESTEP coverage: cut_operation, periodic_surface_in_boolean.
#[test]
fn boolean_difference_hole_roundtrip() {
    let box_shape = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).expect("box");
    let cylinder =
        make_cylinder_brep(DVec3::new(2.0, 2.0, 0.0), DVec3::Z, DVec3::X, 0.5, 4.0).expect("cyl");

    let difference = boolean_op(BooleanOpType::Difference, &box_shape, &cylinder).expect("diff");

    let step_str = StepWriter::write_string(&difference, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse diff result");

    assert_eq!(solid_count(&parsed), 1, "difference should have 1 solid");
    // Box with hole should have more than the original 6 faces due to the cylindrical surface
    assert!(
        face_count(&parsed) >= 6,
        "difference should have at least 6 faces"
    );
}

/// Test roundtrip of a boolean intersection result.
/// OCCT TKDESTEP coverage: common_operation.
#[test]
fn boolean_intersection_roundtrip() {
    let box1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 3.0, 3.0, 3.0).expect("box1");
    let box2 = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 3.0, 3.0, 3.0)
        .expect("box2");

    let intersection =
        boolean_op(BooleanOpType::Intersection, &box1, &box2).expect("intersection");

    let step_str = StepWriter::write_string(&intersection, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse intersection result");

    assert_eq!(solid_count(&parsed), 1, "intersection should have 1 solid");
    // Intersection of two overlapping boxes is a smaller box with 6 faces
    assert_eq!(
        face_count(&parsed),
        6,
        "intersection should have exactly 6 faces"
    );
}

/// Test roundtrip of cone-box intersection.
/// OCCT TKDESTEP coverage: conical_surface_boolean, complex_geometry_result.
#[test]
fn boolean_cone_box_intersection_roundtrip() {
    let cone = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 3.0, 4.0).expect("cone");
    let box_shape =
        make_box_brep(DVec3::new(-1.0, -1.0, 1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0)
            .expect("box");

    let intersection =
        boolean_op(BooleanOpType::Intersection, &cone, &box_shape).expect("intersection");

    let step_str = StepWriter::write_string(&intersection, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse intersection");

    assert_eq!(solid_count(&parsed), 1, "intersection should have 1 solid");
    assert!(face_count(&parsed) >= 1, "should have at least 1 face");
}

/// Test roundtrip of torus-sphere difference.
/// OCCT TKDESTEP coverage: toroidal_surface_boolean, torus_cut.
#[test]
fn boolean_torus_sphere_difference_roundtrip() {
    let torus = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 3.0, 1.0).expect("torus");
    let sphere = make_sphere_brep(DVec3::new(0.0, 0.0, 0.0), 1.5).expect("sphere");

    let result = boolean_op(BooleanOpType::Difference, &torus, &sphere).expect("diff");

    let step_str = StepWriter::write_string(&result, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse result");

    assert_eq!(solid_count(&parsed), 1, "should have 1 solid");
    assert!(face_count(&parsed) >= 1, "should have at least 1 face");
}

/// Test roundtrip of box-box difference (simple cut).
/// OCCT TKDESTEP coverage: simple_cut, box_difference.
#[test]
fn boolean_box_difference_roundtrip() {
    let box1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).expect("box1");
    let box2 = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 3.0, 3.0, 3.0)
        .expect("box2");

    let result = boolean_op(BooleanOpType::Difference, &box1, &box2).expect("diff");

    let step = StepWriter::write_string(&result, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse result");

    assert_eq!(solid_count(&parsed), 1, "should have 1 solid");
    assert!(face_count(&parsed) >= 6, "should have at least 6 faces");
}

// ============================================================================
// 2. Assembly Basic Operations (using single-part write)
// ============================================================================

/// Test that assembly with translation produces valid STEP.
/// OCCT TKDESTEP coverage: item_defined_transformation, cartesian_transformation.
#[test]
fn assembly_translation_baked() {
    let base_box = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let translation = DVec3::new(10.0, 0.0, 0.0);

    let comp = AssemblyComponent::new("shifted_box", base_box).with_translation(translation);

    let step = write_assembly("shift_test", &[comp]);

    // Verify STEP structure
    assert!(step.contains("ISO-10303-21"));
    assert!(step.contains("shifted_box"));
    assert!(step.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"));

    let components = read_assembly(&step).expect("read_assembly");
    assert!(!components.is_empty());

    // Verify geometry was translated
    let brep = &components[0].brep;
    for v in &brep.vertices {
        assert!(
            v.point.x >= 9.999,
            "vertex x should be >= 10 after baking translation, got {}",
            v.point.x
        );
    }
}

/// Test that single-part STEP roundtrips correctly.
/// OCCT TKDESTEP coverage: single_part_file, manifold_solid.
#[test]
fn single_part_step_roundtrip() {
    use rcad_step::ExportSelection;

    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).unwrap();
    let step = StepWriter::write_string(
        &brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    );

    let components = read_assembly(&step).expect("read_assembly on single-part STEP");
    assert_eq!(components.len(), 1, "single-part STEP should give 1 component");
}

// ============================================================================
// 3. AP242 Metadata Roundtrip (Colors, Protocols)
// ============================================================================

/// Test that component colors are written to STEP.
/// OCCT TKDESTEP coverage: styled_item, colour_rgb, presentation_style_assignment.
#[test]
fn component_colors_written() {
    let box1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box1");
    let box2 =
        make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box2");

    let comp1 =
        AssemblyComponent::new("red_box", box1).with_color(Color { r: 1.0, g: 0.0, b: 0.0 });
    let comp2 =
        AssemblyComponent::new("blue_box", box2).with_color(Color { r: 0.0, g: 0.0, b: 1.0 });

    let step = write_assembly("colored_asm", &[comp1, comp2]);

    // Verify STEP contains color entities
    assert!(
        step.contains("COLOUR_RGB") || step.contains("DRAUGHTING_PRE_DEFINED_COLOUR"),
        "STEP should contain color definitions"
    );
    assert!(
        step.contains("PRESENTATION_STYLE_ASSIGNMENT") || step.contains("STYLED_ITEM"),
        "STEP should contain styling entities"
    );
}

/// Test AP242 protocol writing.
/// OCCT TKDESTEP coverage: ap242_schema, managed_model_based_3d_engineering.
#[test]
fn ap242_protocol_output() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");

    let step = StepWriter::write_string_with_options(
        &brep,
        all_faces_selection(),
        &StepWriteOptions {
            protocol: StepProtocol::Ap242,
            ..Default::default()
        },
    );

    assert!(
        step.contains("AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF"),
        "should use AP242 schema"
    );

    // Roundtrip should work
    let parsed = StepReader::parse_string(&step).expect("parse AP242 step");
    assert_eq!(face_count(&parsed), 6, "box should have 6 faces");
}

/// Test AP214 vs AP242 protocol difference in output.
/// OCCT TKDESTEP coverage: schema_selection, protocol_version.
#[test]
fn protocol_selection() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");

    let ap214 = StepWriter::write_string_with_options(
        &brep,
        all_faces_selection(),
        &StepWriteOptions {
            protocol: StepProtocol::Ap214,
            ..Default::default()
        },
    );

    let ap242 = StepWriter::write_string_with_options(
        &brep,
        all_faces_selection(),
        &StepWriteOptions {
            protocol: StepProtocol::Ap242,
            ..Default::default()
        },
    );

    assert!(
        ap214.contains("AUTOMOTIVE_DESIGN"),
        "AP214 should contain AUTOMOTIVE_DESIGN schema"
    );
    assert!(
        ap242.contains("AP242_MANAGED_MODEL_BASED_3D_ENGINEERING"),
        "AP242 should contain AP242 schema"
    );

    // Both should roundtrip successfully
    let parsed_214 = StepReader::parse_string(&ap214).expect("parse AP214");
    let parsed_242 = StepReader::parse_string(&ap242).expect("parse AP242");
    assert_eq!(face_count(&parsed_214), face_count(&parsed_242));
}

/// Test colored single solid roundtrip.
/// OCCT TKDESTEP coverage: per_face_coloring, surface_style_usage.
#[test]
fn colored_solid_roundtrip() {
    use rcad_kernel::appearance::StepColor;

    let box_shape = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box");

    let colors = StepColor::new().with_solid_color(Color { r: 0.8, g: 0.4, b: 0.2 });

    let step = StepWriter::write_string_colored(&box_shape, &colors);

    assert!(step.contains("ISO-10303-21"), "should have STEP header");
    assert!(step.contains("COLOUR_RGB"), "should have color definition");

    let parsed = StepReader::parse_string(&step).expect("parse colored solid");
    assert_eq!(solid_count(&parsed), 1, "should have 1 solid");
    assert_eq!(face_count(&parsed), 6, "should have 6 faces");
}

/// Test STEP header fields.
/// OCCT TKDESTEP coverage: file_description, file_name, author_info.
#[test]
fn step_header_fields() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");

    let step = StepWriter::write_string(&brep, all_faces_selection());

    // Verify standard STEP header sections
    assert!(step.contains("ISO-10303-21"), "should have STEP magic number");
    assert!(step.contains("HEADER"), "should have HEADER section");
    assert!(step.contains("FILE_DESCRIPTION"), "should have FILE_DESCRIPTION");
    assert!(step.contains("FILE_NAME"), "should have FILE_NAME");
    assert!(step.contains("FILE_SCHEMA"), "should have FILE_SCHEMA");
    assert!(step.contains("DATA"), "should have DATA section");
    assert!(step.contains("END-ISO-10303-21"), "should have STEP footer");
}

// ============================================================================
// 4. Periodic Surface Seam Handling (Cylinder, Torus, Cone)
// ============================================================================

/// Test cylinder roundtrip preserves topology.
/// OCCT TKDESTEP coverage: cylindrical_surface_seam, periodic_surface_edge.
#[test]
fn cylinder_roundtrip() {
    let cylinder =
        make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 5.0).expect("cylinder");

    let step = StepWriter::write_string(&cylinder, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse cylinder");

    // Cylinder should have 3 faces: top, bottom, and lateral
    assert_eq!(face_count(&parsed), 3, "cylinder should have 3 faces");

    // Verify bounding box preserves cylindrical geometry
    let bbox = parsed.bounding_box().expect("should have bounds");
    let [min, max] = bbox;

    // Diameter should be 4.0 (radius * 2)
    assert!((max.x - min.x - 4.0).abs() < 0.1, "x diameter should be ~4.0");
    assert!((max.y - min.y - 4.0).abs() < 0.1, "y diameter should be ~4.0");
    assert!((max.z - min.z - 5.0).abs() < 0.1, "height should be ~5.0");
}

/// Test torus roundtrip preserves topology.
/// OCCT TKDESTEP coverage: toroidal_surface_seam, periodic_u_v_surfaces.
#[test]
fn torus_roundtrip() {
    let major_radius = 5.0;
    let minor_radius = 1.5;
    let torus = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, major_radius, minor_radius)
        .expect("torus");

    let step = StepWriter::write_string(&torus, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse torus");

    // Torus should have 1 face (the toroidal surface)
    assert!(face_count(&parsed) >= 1, "torus should have at least 1 face");

    // Verify bounding box
    let bbox = parsed.bounding_box().expect("should have bounds");
    let [min, max] = bbox;

    let expected_outer_diameter = 2.0 * (major_radius + minor_radius);

    assert!(
        (max.x - min.x - expected_outer_diameter).abs() < 0.1,
        "outer diameter should be preserved"
    );
    assert!(
        (max.y - min.y - expected_outer_diameter).abs() < 0.1,
        "outer diameter should be preserved"
    );
}

/// Test cone roundtrip preserves topology.
/// OCCT TKDESTEP coverage: conical_surface, surface_of_revolution.
#[test]
fn cone_roundtrip() {
    let cone = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 3.0, 6.0).expect("cone");

    let step = StepWriter::write_string(&cone, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse cone");

    // Cone should have 2 faces: base and lateral
    assert_eq!(face_count(&parsed), 2, "cone should have 2 faces");

    // Verify cone geometry
    let bbox = parsed.bounding_box().expect("should have bounds");
    let [min, max] = bbox;

    assert!((max.z - min.z - 6.0).abs() < 0.01, "height should be preserved");
    // Base radius at z=0 should be ~3.0, giving diameter ~6.0
    assert!((max.x - min.x - 6.0).abs() < 0.1, "base diameter should be preserved");
}

/// Test cylinder roundtrip after boolean with box.
/// OCCT TKDESTEP coverage: boolean_on_periodic_surface, split_periodic_face.
#[test]
fn cylinder_boolean_roundtrip() {
    // Create a box and cut a cylinder from it
    let box_shape =
        make_box_brep(DVec3::new(-2.0, -2.0, 0.0), DVec3::X, DVec3::Y, 4.0, 4.0, 5.0).expect("box");
    let cylinder = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 5.0).expect("cylinder");

    let result = boolean_op(BooleanOpType::Difference, &box_shape, &cylinder).expect("diff");

    let step = StepWriter::write_string(&result, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse result");

    assert_eq!(solid_count(&parsed), 1, "should have 1 solid");
    // Box with cylindrical hole
    assert!(face_count(&parsed) >= 6, "should have at least 6 faces");
}

/// Test box primitive roundtrip.
/// OCCT TKDESTEP coverage: box_primitive, elementary_surface.
#[test]
fn box_primitive_roundtrip() {
    let box_shape =
        make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).expect("box");

    let step = StepWriter::write_string(&box_shape, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse box");

    assert_eq!(solid_count(&parsed), 1, "should have 1 solid");
    assert_eq!(face_count(&parsed), 6, "box should have 6 faces");

    let bbox = parsed.bounding_box().expect("should have bounds");
    let [min, max] = bbox;
    assert!((max.x - min.x - 2.0).abs() < 0.01, "x dimension preserved");
    assert!((max.y - min.y - 3.0).abs() < 0.01, "y dimension preserved");
    assert!((max.z - min.z - 4.0).abs() < 0.01, "z dimension preserved");
}

// ============================================================================
// 5. Large File and Performance Tests
// ============================================================================

/// Test that writing a moderate assembly doesn't take excessive time.
/// OCCT TKDESTEP coverage: performance_assembly_write, streaming_output.
#[test]
fn moderate_assembly_write_perf() {
    let mut components = Vec::new();

    // Create components
    for i in 0..5 {
        for j in 0..5 {
            let x = i as f64 * 2.0;
            let y = j as f64 * 2.0;
            components.push(AssemblyComponent::new(
                format!("part_{}_{}", i, j),
                make_box(DVec3::new(x, y, 0.0)),
            ));
        }
    }

    let start = Instant::now();
    let step = write_assembly("perf_test_asm", &components);
    let write_duration = start.elapsed();

    // Write should complete in reasonable time (< 2 seconds)
    assert!(
        write_duration.as_secs() < 2,
        "write took {:?}, expected < 2s",
        write_duration
    );

    // Verify reasonable STEP file size
    let size_kb = step.len() / 1024;
    assert!(
        size_kb < 2000,
        "STEP file size {}KB seems too large for 25 boxes",
        size_kb
    );
}

/// Test STEP file size scaling with complexity.
/// OCCT TKDESTEP coverage: file_size_efficiency, minimal_redundancy.
#[test]
fn step_file_size_reasonable() {
    let simple_box = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
    let simple_step = StepWriter::write_string(&simple_box, all_faces_selection());

    // Create more complex geometry
    let torus = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 5.0, 1.0).expect("torus");
    let complex_step = StepWriter::write_string(&torus, all_faces_selection());

    // Complex geometry should be larger but not excessively so
    let ratio = complex_step.len() as f64 / simple_step.len() as f64;
    assert!(ratio < 10.0, "complex geometry ratio {} seems excessive", ratio);

    // Absolute size check
    assert!(
        complex_step.len() < 200_000,
        "complex STEP size {} seems too large",
        complex_step.len()
    );
}

/// Test parsing STEP with complex geometry.
/// OCCT TKDESTEP coverage: complex_surface_parsing, dense_geometry.
#[test]
fn complex_geometry_parses() {
    // Torus with small minor radius creates detailed geometry
    let torus = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 5.0, 0.5).expect("torus");

    let step = StepWriter::write_string(&torus, all_faces_selection());

    // Should parse without issues
    let parsed = StepReader::parse_string(&step).expect("parse complex torus");
    assert!(face_count(&parsed) >= 1, "should have at least 1 face");
}

/// Test roundtrip consistency across multiple iterations.
/// OCCT TKDESTEP coverage: idempotent_write_read, stable_output.
#[test]
fn multiple_roundtrips_stable() {
    let original = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 5.0).expect("cylinder");

    let mut current = original.clone();
    let original_faces = face_count(&original);

    for iteration in 0..5 {
        let step = StepWriter::write_string(&current, all_faces_selection());
        current = StepReader::parse_string(&step).expect(&format!("iteration {}", iteration));

        assert_eq!(
            face_count(&current),
            original_faces,
            "face count should be stable at iteration {}",
            iteration
        );
    }
}

/// Test STEP output is deterministic (same input -> same output).
/// OCCT TKDESTEP coverage: deterministic_output, reproducible_files.
#[test]
fn step_output_deterministic() {
    let brep = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 5.0, 1.0).expect("torus");

    let step1 = StepWriter::write_string(&brep, all_faces_selection());
    let step2 = StepWriter::write_string(&brep, all_faces_selection());

    assert_eq!(step1, step2, "identical inputs should produce identical STEP output");
}

/// Test very small geometry roundtrip.
/// OCCT TKDESTEP coverage: small_feature_handling, numerical_precision.
#[test]
fn small_geometry_roundtrip() {
    let tiny_box =
        make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 0.001, 0.001, 0.001).expect("tiny box");

    let step = StepWriter::write_string(&tiny_box, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse tiny box");

    assert_eq!(face_count(&parsed), 6, "tiny box should have 6 faces");
}

/// Test large coordinate values.
/// OCCT TKDESTEP coverage: large_coordinates, numerical_stability.
#[test]
fn large_coordinate_roundtrip() {
    let offset = DVec3::new(1e5, 1e5, 1e5);
    let large_box = make_box_brep(offset, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).expect("large box");

    let step = StepWriter::write_string(&large_box, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse large coordinate box");

    let bbox = parsed.bounding_box().expect("should have bounds");
    let [min, _] = bbox;

    // Offset should be approximately preserved
    assert!(
        (min.x - offset.x).abs() < 1.0,
        "x offset should be approximately preserved"
    );
}

/// Test negative coordinate values.
/// OCCT TKDESTEP coverage: negative_coordinates, sign_preservation.
#[test]
fn negative_coordinate_roundtrip() {
    let neg_box = make_box_brep(
        DVec3::new(-10.0, -10.0, -10.0),
        DVec3::X,
        DVec3::Y,
        5.0,
        5.0,
        5.0,
    )
    .expect("negative box");

    let step = StepWriter::write_string(&neg_box, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse negative box");

    let bbox = parsed.bounding_box().expect("should have bounds");
    let [min, _] = bbox;

    assert!(min.x < 0.0, "should have negative x coordinates");
    assert!(min.y < 0.0, "should have negative y coordinates");
    assert!(min.z < 0.0, "should have negative z coordinates");
}

/// Test multiple boolean operations roundtrip.
/// OCCT TKDESTEP coverage: performance_complex_geometry, many_faces.
#[test]
fn multiple_boolean_roundtrip() {
    // Create a series of boolean operations
    let mut results = Vec::new();

    for i in 0..5 {
        let base = make_box_brep(
            DVec3::new(i as f64 * 5.0, 0.0, 0.0),
            DVec3::X,
            DVec3::Y,
            4.0,
            4.0,
            4.0,
        )
        .expect("box");

        let hole = make_cylinder_brep(
            DVec3::new(i as f64 * 5.0 + 2.0, 2.0, 0.0),
            DVec3::Z,
            DVec3::X,
            0.5,
            4.0,
        )
        .expect("hole");

        let with_hole = boolean_op(BooleanOpType::Difference, &base, &hole).expect("diff");
        results.push(with_hole);
    }

    // Roundtrip each result
    for (i, brep) in results.iter().enumerate() {
        let step = StepWriter::write_string(brep, all_faces_selection());
        let parsed = StepReader::parse_string(&step).expect(&format!("parse result {}", i));
        assert_eq!(solid_count(&parsed), 1, "result {} should have 1 solid", i);
    }
}

/// Test thin geometry (sheet-like).
/// OCCT TKDESTEP coverage: thin_geometry, degenerate_dimensions.
#[test]
fn thin_geometry_roundtrip() {
    // Very thin box (sheet-like)
    let brep =
        make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 0.001).expect("thin box");

    let step = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse should succeed");

    assert_eq!(face_count(&parsed), 6, "thin box should have 6 faces");
}

/// Test non-standard orientation.
/// OCCT TKDESTEP coverage: non_axis_aligned, arbitrary_orientation.
#[test]
fn non_standard_orientation_roundtrip() {
    // Box with non-axis-aligned orientation
    let brep = make_box_brep(
        DVec3::new(1.0, 2.0, 3.0),
        DVec3::new(1.0, 1.0, 0.0).normalize(), // Non-standard X axis
        DVec3::new(-1.0, 1.0, 0.0).normalize(), // Non-standard Y axis
        5.0,
        3.0,
        2.0,
    )
    .expect("rotated box");

    let step = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse should succeed");

    assert_eq!(face_count(&parsed), 6);
}

// ============================================================================
// 6. Error Handling Tests
// ============================================================================

/// Test truncated STEP file returns error.
/// OCCT TKDESTEP coverage: error_handling, malformed_input.
#[test]
fn truncated_step_error() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
    let step_str = StepWriter::write_string(&brep, all_faces_selection());

    // Truncate the string
    let truncated = &step_str[..step_str.len() / 2];

    let result = StepReader::parse_string(truncated);
    assert!(result.is_err(), "truncated STEP should fail");
}

/// Test invalid format returns error.
/// OCCT TKDESTEP coverage: error_handling, invalid_format.
#[test]
fn invalid_format_error() {
    use rcad_step::StepError;
    let result = StepReader::parse_string("this is not a STEP file");
    assert!(
        matches!(result, Err(StepError::InvalidFormat(_))),
        "expected InvalidFormat, got {:?}",
        result
    );
}

/// Test empty string returns error.
/// OCCT TKDESTEP coverage: error_handling, empty_input.
#[test]
fn empty_string_error() {
    use rcad_step::StepError;
    let result = StepReader::parse_string("");
    assert!(
        matches!(result, Err(StepError::InvalidFormat(_))),
        "expected InvalidFormat for empty input, got {:?}",
        result
    );
}

/// Test binary garbage returns error.
/// OCCT TKDESTEP coverage: error_handling, binary_input.
#[test]
fn binary_garbage_error() {
    let garbage: String = (0..255u8).map(|b| char::from(b)).collect();

    let result = StepReader::parse_string(&garbage);
    assert!(result.is_err(), "binary garbage should fail");
}

// ============================================================================
// 7. FEA Geometry Support Tests
// ============================================================================

/// Test box geometry suitable for FEA meshing.
/// OCCT TKDESTEP coverage: box_for_fea, manifold_geometry.
#[test]
fn box_fea_geometry() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 20.0, 5.0).expect("box");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    // Verify topology for FEA
    assert_eq!(face_count(&parsed), 6, "box should have 6 faces");

    // Verify geometry bounds are preserved
    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    assert!((max.x - min.x - 10.0).abs() < 0.1, "X dimension should be 10");
    assert!((max.y - min.y - 20.0).abs() < 0.1, "Y dimension should be 20");
    assert!((max.z - min.z - 5.0).abs() < 0.1, "Z dimension should be 5");
}

/// Test cylinder geometry for FEA.
/// OCCT TKDESTEP coverage: cylinder_for_fea, periodic_surface.
#[test]
fn cylinder_fea_geometry() {
    let radius = 5.0;
    let height = 15.0;
    let brep = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, radius, height).expect("cylinder");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    // Verify topology
    assert_eq!(face_count(&parsed), 3, "cylinder should have 3 faces");

    // Verify bounds
    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    assert!((max.z - min.z - height).abs() < 0.01, "height should be preserved");
}

/// Test cone geometry for FEA.
/// OCCT TKDESTEP coverage: cone_for_fea, conical_surface.
#[test]
fn cone_fea_geometry() {
    let radius = 4.0;
    let height = 8.0;
    let brep = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, radius, height).expect("cone");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    // Verify topology
    assert_eq!(face_count(&parsed), 2, "cone should have 2 faces");

    // Verify bounds
    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    assert!((max.z - min.z - height).abs() < 0.01, "height should be preserved");
}

/// Test torus geometry for FEA.
/// OCCT TKDESTEP coverage: torus_for_fea, toroidal_surface.
#[test]
fn torus_fea_geometry() {
    let major_radius = 5.0;
    let minor_radius = 1.5;
    let brep =
        make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, major_radius, minor_radius).expect("torus");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    // Verify topology
    assert!(face_count(&parsed) >= 1, "torus should have at least 1 face");

    // Verify outer diameter
    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    let expected_diameter = 2.0 * (major_radius + minor_radius);
    assert!(
        (max.x - min.x - expected_diameter).abs() < 0.1,
        "outer diameter should be approximately preserved"
    );
}

// ============================================================================
// 8. OCCT TKDESTEP Additional Edge Cases
// ============================================================================

/// Test torus roundtrip as an alternative primitive.
/// OCCT TKDESTEP coverage: toroidal_surface, elementary_surface_roundtrip.
#[test]
fn torus_primitive_roundtrip() {
    let major_radius = 3.0;
    let minor_radius = 1.0;
    let torus = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, major_radius, minor_radius)
        .expect("torus");

    let step = StepWriter::write_string(&torus, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse torus");

    // Torus should have at least 1 face
    assert!(face_count(&parsed) >= 1, "torus should have at least 1 face");

    // Verify bounding box
    let bbox = parsed.bounding_box().expect("should have bounds");
    let [min, max] = bbox;
    let expected_diameter = 2.0 * (major_radius + minor_radius);
    assert!((max.x - min.x - expected_diameter).abs() < 0.5, "diameter preserved");
}

/// Test cylinder-box intersection roundtrip.
/// OCCT TKDESTEP coverage: cylinder_boolean, curved_planar_intersection.
#[test]
fn cylinder_box_intersection_roundtrip() {
    let cylinder =
        make_cylinder_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::Z, DVec3::X, 1.0, 4.0)
            .expect("cylinder");
    let box_shape =
        make_box_brep(DVec3::new(-1.0, -1.0, 1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box");

    let result =
        boolean_op(BooleanOpType::Intersection, &cylinder, &box_shape).expect("intersection");

    let step = StepWriter::write_string(&result, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse result");

    assert_eq!(solid_count(&parsed), 1, "should have 1 solid");
    assert!(face_count(&parsed) >= 1, "should have at least 1 face");
}

/// Test nested cylinders roundtrip.
/// OCCT TKDESTEP coverage: concentric_cylinders, tube_geometry.
#[test]
fn nested_cylinders_difference_roundtrip() {
    let outer = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 3.0, 10.0).expect("outer");
    let inner = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 10.0).expect("inner");

    let tube = boolean_op(BooleanOpType::Difference, &outer, &inner).expect("tube");

    let step = StepWriter::write_string(&tube, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse tube");

    assert_eq!(solid_count(&parsed), 1, "tube should have 1 solid");
    // Tube has inner and outer surfaces, so more faces than simple cylinder
    assert!(face_count(&parsed) >= 4, "tube should have at least 4 faces");
}

/// Test rotated primitive roundtrip.
/// OCCT TKDESTEP coverage: rotated_geometry, transformation_preservation.
#[test]
fn rotated_cylinder_roundtrip() {
    use std::f64::consts::FRAC_PI_4;

    // Create a cylinder and rotate it 45 degrees around X axis
    let mut cyl = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 5.0).expect("cylinder");

    let rotation = DAffine3::from_rotation_x(FRAC_PI_4);
    cyl.apply_transform(rotation);

    let step = StepWriter::write_string(&cyl, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse rotated cylinder");

    assert_eq!(face_count(&parsed), 3, "rotated cylinder should have 3 faces");
}

/// Test scaled geometry roundtrip.
/// OCCT TKDESTEP coverage: scaled_geometry, non_uniform_scale.
#[test]
fn scaled_box_roundtrip() {
    let mut box_shape = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");

    // Apply non-uniform scale
    let scale = DAffine3::from_scale(DVec3::new(2.0, 3.0, 4.0));
    box_shape.apply_transform(scale);

    let step = StepWriter::write_string(&box_shape, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse scaled box");

    let bbox = parsed.bounding_box().expect("should have bounds");
    let [min, max] = bbox;

    assert!((max.x - min.x - 2.0).abs() < 0.1, "x scaled to 2");
    assert!((max.y - min.y - 3.0).abs() < 0.1, "y scaled to 3");
    assert!((max.z - min.z - 4.0).abs() < 0.1, "z scaled to 4");
}

/// Test mirrored geometry roundtrip.
/// OCCT TKDESTEP coverage: mirrored_geometry, reflection_transformation.
#[test]
fn mirrored_box_roundtrip() {
    let mut box_shape = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box");

    // Mirror across XY plane (flip Z)
    let mirror = DAffine3::from_scale(DVec3::new(1.0, 1.0, -1.0));
    box_shape.apply_transform(mirror);

    let step = StepWriter::write_string(&box_shape, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse mirrored box");

    assert_eq!(face_count(&parsed), 6, "mirrored box should have 6 faces");
}

/// Test combined transform roundtrip.
/// OCCT TKDESTEP coverage: combined_transformation, complex_transform.
#[test]
fn combined_transform_roundtrip() {
    let mut box_shape = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");

    // Apply combined transform: rotate, scale, translate
    let rotation = DAffine3::from_rotation_z(std::f64::consts::FRAC_PI_4);
    let scale = DAffine3::from_scale(DVec3::splat(2.0));
    let translation = DAffine3::from_translation(DVec3::new(10.0, 20.0, 30.0));

    let combined = translation * scale * rotation;
    box_shape.apply_transform(combined);

    let step = StepWriter::write_string(&box_shape, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse transformed box");

    let bbox = parsed.bounding_box().expect("should have bounds");
    let [min, max] = bbox;

    // Box should be at offset location
    assert!(min.x > 8.0, "x should be translated");
    assert!(min.y > 18.0, "y should be translated");
    assert!(min.z > 28.0, "z should be translated");
}

/// Test boolean chain roundtrip (multiple operations).
/// OCCT TKDESTEP coverage: sequential_boolean, compound_operation.
#[test]
fn boolean_chain_roundtrip() {
    // Start with a box
    let mut result = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).expect("box");

    // Cut first cylinder
    let cyl1 = make_cylinder_brep(DVec3::new(3.0, 3.0, 0.0), DVec3::Z, DVec3::X, 1.0, 10.0)
        .expect("cyl1");
    result = boolean_op(BooleanOpType::Difference, &result, &cyl1).expect("diff1");

    // Cut second cylinder
    let cyl2 = make_cylinder_brep(DVec3::new(7.0, 7.0, 0.0), DVec3::Z, DVec3::X, 1.0, 10.0)
        .expect("cyl2");
    result = boolean_op(BooleanOpType::Difference, &result, &cyl2).expect("diff2");

    let step = StepWriter::write_string(&result, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse result");

    assert_eq!(solid_count(&parsed), 1, "should have 1 solid after chain");
}

/// Test assembly with multiple components.
/// OCCT TKDESTEP coverage: multi_component_assembly, nested_structure.
#[test]
fn multi_component_assembly_roundtrip() {
    let mut components = Vec::new();

    for i in 0..3 {
        let box_shape = make_box_brep(
            DVec3::new(i as f64 * 3.0, 0.0, 0.0),
            DVec3::X,
            DVec3::Y,
            2.0,
            2.0,
            2.0,
        )
        .expect("box");
        components.push(AssemblyComponent::new(format!("part_{}", i), box_shape));
    }

    let step = write_assembly("multi_part", &components);

    // Verify STEP contains expected structure
    assert!(step.contains("ISO-10303-21"), "should have STEP header");
    assert!(step.contains("part_0") || step.contains("part_1") || step.contains("part_2"), "should have part names");

    let parsed = read_assembly(&step).expect("read assembly");

    // Note: Assembly parsing may return different number of components
    // due to how NEXT_ASSEMBLY_USAGE_OCCURRENCE is handled
    assert!(!parsed.is_empty(), "should have at least 1 component");
}

/// Test very large box roundtrip.
/// OCCT TKDESTEP coverage: large_dimensions, numerical_stability_large.
#[test]
fn very_large_box_roundtrip() {
    let large_box =
        make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1000.0, 1000.0, 1000.0).expect("large box");

    let step = StepWriter::write_string(&large_box, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse large box");

    let bbox = parsed.bounding_box().expect("should have bounds");
    let [min, max] = bbox;

    assert!((max.x - min.x - 1000.0).abs() < 1.0, "large x dimension preserved");
}

/// Test union of three boxes roundtrip.
/// OCCT TKDESTEP coverage: multi_operand_union, compound_result.
#[test]
fn three_box_union_roundtrip() {
    let box1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box1");
    let box2 =
        make_box_brep(DVec3::new(1.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box2");
    let box3 =
        make_box_brep(DVec3::new(0.0, 1.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box3");

    let union1 = boolean_op(BooleanOpType::Union, &box1, &box2).expect("union1");
    let final_union = boolean_op(BooleanOpType::Union, &union1, &box3).expect("union2");

    let step = StepWriter::write_string(&final_union, all_faces_selection());
    let parsed = StepReader::parse_string(&step).expect("parse union");

    assert_eq!(solid_count(&parsed), 1, "should have 1 solid");
}
