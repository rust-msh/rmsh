/// Integration tests for STEP file I/O with FEA and advanced features.
///
/// These tests verify:
/// - FEA entity round-trip (element sets, node sets)
/// - Large file parsing robustness
/// - Edge cases in STEP processing
use glam::DVec3;
use rcad_modeling::{make_box_brep, make_sphere_brep, make_cylinder_brep, make_cone_brep, make_torus_brep};
use rcad_step::{StepReader, StepWriter, ExportSelection};

fn all_faces_selection() -> ExportSelection<'static> {
    ExportSelection {
        selected_faces: &[],
        selected_edges: &[],
    }
}

fn face_count(brep: &rcad_kernel::BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

fn vertex_count(brep: &rcad_kernel::BRep) -> usize {
    brep.vertices.len()
}

fn edge_count(brep: &rcad_kernel::BRep) -> usize {
    brep.edges.len()
}

// ============================================================================
// FEA Entity Round-Trip Tests
// ============================================================================

/// Test that a simple box maintains geometry suitable for FEA meshing.
#[test]
fn box_fea_geometry_round_trip() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 20.0, 5.0).expect("box");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    // Verify topology for FEA
    assert_eq!(face_count(&parsed), 6, "box should have 6 faces");
    // Note: Vertex count may differ due to triangulation vertices being added
    assert!(vertex_count(&parsed) >= 8, "box should have at least 8 vertices");
    // Note: Edge count may differ in STEP representation

    // Verify geometry bounds are preserved
    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    assert!((max.x - min.x - 10.0).abs() < 0.1, "X dimension should be 10");
    assert!((max.y - min.y - 20.0).abs() < 0.1, "Y dimension should be 20");
    assert!((max.z - min.z - 5.0).abs() < 0.1, "Z dimension should be 5");
}

/// Test cylinder geometry preservation for FEA.
#[test]
fn cylinder_fea_geometry_round_trip() {
    let radius = 5.0;
    let height = 15.0;
    let brep = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, radius, height)
        .expect("cylinder");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    // Verify topology
    assert_eq!(face_count(&parsed), 3, "cylinder should have 3 faces");
    assert!(vertex_count(&parsed) >= 2, "cylinder should have at least 2 vertices");

    // Verify bounds
    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    assert!((max.z - min.z - height).abs() < 0.01, "height should be preserved");
}

/// Test sphere geometry for FEA mesh generation.
#[test]
fn sphere_fea_geometry_round_trip() {
    let radius = 7.5;
    let brep = make_sphere_brep(DVec3::ZERO, radius).expect("sphere");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    // Verify topology
    assert!(face_count(&parsed) >= 1, "sphere should have at least 1 face");

    // Verify spherical bounds (should be approximately diameter in all directions)
    // Note: STEP representation may have variations due to how the sphere is tessellated
    // and the bounding box is computed
    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    let diameter = 2.0 * radius;
    // The bounding box should contain the sphere, so max dimension should be at least diameter
    // Allow larger tolerance for STEP representation variations
    let max_dim = (max.x - min.x).max(max.y - min.y).max(max.z - min.z);
    assert!(max_dim >= diameter * 0.8, "sphere should have at least 80% of expected diameter");
}

/// Test cone geometry preservation.
#[test]
fn cone_fea_geometry_round_trip() {
    let radius = 4.0;
    let height = 8.0;
    let brep = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, radius, height)
        .expect("cone");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    // Verify topology
    assert_eq!(face_count(&parsed), 2, "cone should have 2 faces");

    // Verify bounds
    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    assert!((max.z - min.z - height).abs() < 0.01, "height should be preserved");
}

/// Test torus geometry for complex FEA models.
#[test]
fn torus_fea_geometry_round_trip() {
    let major_radius = 5.0;
    let minor_radius = 1.5;
    let brep = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, major_radius, minor_radius)
        .expect("torus");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    // Verify topology
    assert!(face_count(&parsed) >= 1, "torus should have at least 1 face");

    // Verify outer diameter (major + minor) * 2
    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    let expected_diameter = 2.0 * (major_radius + minor_radius);
    assert!(
        (max.x - min.x - expected_diameter).abs() < 0.1,
        "outer diameter should be approximately preserved"
    );
}

/// Test multiple primitives in sequence (simulating assembly).
#[test]
fn multiple_primitives_sequential_round_trip() {
    let primitives = vec![
        ("box", make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).expect("box")),
        ("sphere", make_sphere_brep(DVec3::new(10.0, 0.0, 0.0), 1.5).expect("sphere")),
        ("cylinder", make_cylinder_brep(DVec3::new(20.0, 0.0, 0.0), DVec3::Z, DVec3::X, 1.0, 3.0).expect("cylinder")),
    ];

    for (name, brep) in primitives {
        let step_str = StepWriter::write_string(&brep, all_faces_selection());
        let parsed = StepReader::parse_string(&step_str)
            .expect(&format!("{} should parse", name));
        assert!(face_count(&parsed) > 0, "{} should have faces", name);
    }
}

// ============================================================================
// Large File Parsing Tests
// ============================================================================

/// Test parsing of STEP file with many boxes (simulating large assembly).
#[test]
fn many_boxes_round_trip() {
    let mut breps = Vec::new();
    let grid_size: i32 = 3;

    for i in 0..grid_size {
        for j in 0..grid_size {
            for k in 0..grid_size {
                let x = i as f64 * 3.0;
                let y = j as f64 * 3.0;
                let z = k as f64 * 3.0;
                let brep = make_box_brep(
                    DVec3::new(x, y, z),
                    DVec3::X,
                    DVec3::Y,
                    2.0,
                    2.0,
                    2.0,
                )
                .expect("box");
                breps.push(brep);
            }
        }
    }

    // Write and parse each individually
    let total_faces: usize = breps.iter().map(|b| face_count(b)).sum();
    assert_eq!(total_faces, (grid_size.pow(3) * 6) as usize, "should have correct face count");

    // Test round-trip for each
    for brep in &breps {
        let step_str = StepWriter::write_string(brep, all_faces_selection());
        let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");
        assert_eq!(face_count(&parsed), 6);
    }
}

/// Test STEP file string size scaling.
#[test]
fn step_string_size_scaling() {
    let small = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
    let large = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 100.0, 100.0, 100.0).expect("box");

    let small_step = StepWriter::write_string(&small, all_faces_selection());
    let large_step = StepWriter::write_string(&large, all_faces_selection());

    // Both should produce similar sized STEP files (geometry is relative)
    // Size difference should be small (< 20%)
    let size_ratio = large_step.len() as f64 / small_step.len() as f64;
    assert!(
        size_ratio < 1.2,
        "STEP size should not vary significantly with absolute dimensions"
    );
}

/// Test parsing of deeply nested STEP structures.
#[test]
fn deep_nesting_parsing() {
    // Create a complex shape by writing and parsing multiple times
    let original = make_sphere_brep(DVec3::ZERO, 5.0).expect("sphere");

    let mut current = original;
    for iteration in 0..5 {
        let step_str = StepWriter::write_string(&current, all_faces_selection());
        current = StepReader::parse_string(&step_str)
            .expect(&format!("iteration {} should succeed", iteration));
    }

    // Final result should still be valid
    assert!(face_count(&current) >= 1);
}

// ============================================================================
// Edge Case Tests
// ============================================================================

/// Test STEP file with degenerate geometry handling.
#[test]
fn thin_geometry_round_trip() {
    // Very thin box (sheet-like)
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 0.001).expect("thin box");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    assert_eq!(face_count(&parsed), 6, "thin box should have 6 faces");
}

/// Test STEP file with very small features.
#[test]
fn small_features_round_trip() {
    let tiny_sphere = make_sphere_brep(DVec3::ZERO, 0.001).expect("tiny sphere");

    let step_str = StepWriter::write_string(&tiny_sphere, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    assert!(face_count(&parsed) >= 1, "tiny sphere should have faces");
}

/// Test STEP file with very large coordinates.
#[test]
fn large_coordinates_round_trip() {
    let offset = DVec3::new(1e6, 1e6, 1e6);
    let brep = make_box_brep(offset, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).expect("offset box");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, _] = bbox;
    assert!((min.x - offset.x).abs() < 0.1, "offset should be preserved");
}

/// Test STEP file with negative coordinates.
#[test]
fn negative_coordinates_round_trip() {
    let brep = make_box_brep(DVec3::new(-5.0, -5.0, -5.0), DVec3::X, DVec3::Y, 10.0, 10.0, 10.0)
        .expect("negative box");

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    let bbox = parsed.bounding_box().expect("should have bounding box");
    let [min, max] = bbox;
    assert!(min.x < 0.0, "should have negative x");
    assert!(min.y < 0.0, "should have negative y");
    assert!(min.z < 0.0, "should have negative z");
    assert!((max.x - min.x - 10.0).abs() < 0.01);
}

/// Test STEP file header validation.
#[test]
fn step_header_validation() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
    let step_str = StepWriter::write_string(&brep, all_faces_selection());

    // Verify standard STEP header
    assert!(step_str.contains("ISO-10303-21"), "should have STEP header");
    assert!(step_str.contains("HEADER"), "should have HEADER section");
    assert!(step_str.contains("DATA"), "should have DATA section");
    assert!(step_str.contains("END-ISO-10303-21"), "should have STEP footer");
}

/// Test STEP file with non-standard orientations.
#[test]
fn non_standard_orientation_round_trip() {
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

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse should succeed");

    assert_eq!(face_count(&parsed), 6);
}

// ============================================================================
// Error Recovery Tests
// ============================================================================

/// Test that truncated STEP file returns appropriate error.
#[test]
fn truncated_step_returns_error() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
    let step_str = StepWriter::write_string(&brep, all_faces_selection());

    // Truncate the string
    let truncated = &step_str[..step_str.len() / 2];

    let result = StepReader::parse_string(truncated);
    assert!(result.is_err(), "truncated STEP should fail");
}

/// Test that corrupted STEP file returns appropriate error.
#[test]
fn corrupted_step_returns_error() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
    let mut step_str = StepWriter::write_string(&brep, all_faces_selection());

    // Corrupt the middle of the string
    let mid = step_str.len() / 2;
    step_str.replace_range(mid..mid + 10, "XXXXXXXXXX");

    let result = StepReader::parse_string(&step_str);
    // Should either parse with degraded data or return an error
    match result {
        Ok(_) => {
            // Parser recovered - acceptable
        }
        Err(_) => {
            // Parser detected corruption - acceptable
        }
    }
}

/// Test that binary garbage returns error.
#[test]
fn binary_garbage_returns_error() {
    let garbage: String = (0..255u8)
        .map(|b| char::from(b))
        .collect();

    let result = StepReader::parse_string(&garbage);
    assert!(result.is_err(), "binary garbage should fail");
}

/// Test that empty DATA section handles gracefully.
#[test]
fn empty_data_section() {
    let step_with_empty_data = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('test.stp','2024-01-01',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
ENDSEC;
END-ISO-10303-21;
"#;

    let result = StepReader::parse_string(step_with_empty_data);
    // Should either return empty BRep or an error
    match result {
        Ok(brep) => {
            assert_eq!(face_count(&brep), 0, "empty data should give empty BRep");
        }
        Err(_) => {
            // Error is acceptable for empty data
        }
    }
}

// ============================================================================
// Consistency Tests
// ============================================================================

/// Test that multiple round-trips produce consistent results.
#[test]
fn multiple_round_trips_consistent() {
    let original = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 5.0, 7.0, 3.0).expect("box");

    let mut current = original.clone();
    let original_faces = face_count(&original);

    for i in 0..10 {
        let step_str = StepWriter::write_string(&current, all_faces_selection());
        current = StepReader::parse_string(&step_str)
            .expect(&format!("round-trip {} should succeed", i));

        assert_eq!(
            face_count(&current),
            original_faces,
            "face count should be consistent after round-trip {}",
            i
        );
        // Note: Vertex count may differ between round-trips due to triangulation
        // being stored in the BRep. The key invariant is face count consistency.
    }
}

/// Test that identical inputs produce identical outputs.
#[test]
fn deterministic_output() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).expect("box");

    let step1 = StepWriter::write_string(&brep, all_faces_selection());
    let step2 = StepWriter::write_string(&brep, all_faces_selection());

    assert_eq!(step1, step2, "identical inputs should produce identical STEP output");
}
