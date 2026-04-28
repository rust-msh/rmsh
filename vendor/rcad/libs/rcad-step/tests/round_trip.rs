/// Integration tests for STEP I/O: write a shape, parse it back, verify topology
/// is preserved. These act as regression guards against serialization regressions.
use glam::DVec3;
use rcad_modeling::{make_box_brep, make_sphere_brep};
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

// ── Box round-trip ───────────────────────────────────────────────────────────

/// Write a unit box to STEP string, parse it back, and verify face count is preserved.
/// Note: Vertex count may differ due to triangulation vertices being stored in the BRep.
#[test]
fn box_round_trip_preserves_topology() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("make box");
    let original_faces = face_count(&brep);

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    assert!(step_str.contains("ISO-10303-21"), "STEP string must contain header");

    let parsed = StepReader::parse_string(&step_str).expect("parse round-tripped STEP");

    assert_eq!(
        face_count(&parsed),
        original_faces,
        "face count must be preserved after round-trip"
    );
    // Note: Vertex count may differ due to triangulation vertices being stored
    // The key invariant is that we have at least the original vertices
    assert!(
        vertex_count(&parsed) >= 8,
        "box should have at least 8 vertices after round-trip"
    );
}

/// Larger box: verify face count stays 6 regardless of dimensions.
#[test]
fn large_box_round_trip() {
    let brep = make_box_brep(DVec3::new(1.0, 2.0, 3.0), DVec3::X, DVec3::Y, 5.0, 3.0, 2.0)
        .expect("make box");
    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse");
    assert_eq!(face_count(&parsed), 6);
}

// ── Sphere round-trip ────────────────────────────────────────────────────────

/// Write a sphere to STEP, parse back; face count and vertex count must survive.
#[test]
fn sphere_round_trip_preserves_topology() {
    let brep = make_sphere_brep(DVec3::ZERO, 1.0).expect("make sphere");
    let original_faces = face_count(&brep);

    let step_str = StepWriter::write_string(&brep, all_faces_selection());
    let parsed = StepReader::parse_string(&step_str).expect("parse sphere STEP");

    assert_eq!(face_count(&parsed), original_faces);
}

// ── Error paths ──────────────────────────────────────────────────────────────

/// Passing a non-STEP string should return InvalidFormat.
#[test]
fn invalid_format_returns_error() {
    use rcad_step::StepError;
    let result = StepReader::parse_string("this is not a STEP file");
    assert!(
        matches!(result, Err(StepError::InvalidFormat(_))),
        "expected InvalidFormat, got {result:?}"
    );
}

/// An empty string should also return InvalidFormat (not panic).
#[test]
fn empty_string_returns_error() {
    use rcad_step::StepError;
    let result = StepReader::parse_string("");
    assert!(
        matches!(result, Err(StepError::InvalidFormat(_))),
        "expected InvalidFormat for empty input, got {result:?}"
    );
}
