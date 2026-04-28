//! Integration tests for `rcad-modeling` primitive builders and sweep operations.

use glam::{DVec2, DVec3};
use rcad_kernel::BRep;
use rcad_modeling::{
    chamfer_edge, extrude, fillet_edge, loft, make_box_brep, make_cone_brep, make_cylinder_brep,
    make_sphere_brep, make_torus_brep, revolve, sweep_pipe,
};

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn face_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum()
}

// ─────────────────────────────────────────────────────────────────────────────
// Primitive face counts
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn box_face_count() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).unwrap();
    assert_eq!(face_count(&brep), 6, "box must have 6 faces");
}

#[test]
fn sphere_face_count_positive() {
    let brep = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
    assert!(face_count(&brep) >= 1, "sphere must have at least 1 face");
    assert!(!brep.solids.is_empty(), "sphere must have a solid");
}

#[test]
fn cylinder_face_count() {
    let brep = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();
    // Cylinder: top cap + bottom cap + lateral face = 3
    assert_eq!(face_count(&brep), 3, "cylinder must have 3 faces");
}

#[test]
fn cone_face_count() {
    let brep = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();
    // Cone: bottom cap + lateral face = 2
    assert_eq!(face_count(&brep), 2, "cone must have 2 faces");
}

#[test]
fn torus_face_count_positive() {
    let brep = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 0.5).unwrap();
    assert!(face_count(&brep) >= 1, "torus must have at least 1 face");
}

// ─────────────────────────────────────────────────────────────────────────────
// Sweep operations
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn extrude_box_face_returns_brep() {
    // Extrude the bottom face (index 0) of a flat box (depth=0.01 ≈ sheet).
    let sheet = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 0.01).unwrap();
    let result = extrude(&sheet, 0, DVec3::Z, 1.0);
    assert!(result.is_ok(), "extrude should succeed: {:?}", result.err());
    let brep = result.unwrap();
    assert!(!brep.solids.is_empty(), "extruded brep must have solids");
    assert!(face_count(&brep) > 0, "extruded brep must have faces");
}

#[test]
fn revolve_box_face_returns_brep() {
    // Revolve a thin box face 90° around the Z axis.
    let sheet = make_box_brep(DVec3::new(1.0, 0.0, 0.0), DVec3::X, DVec3::Y, 0.5, 0.01, 1.0)
        .unwrap();
    let result = revolve(&sheet, 0, DVec3::ZERO, DVec3::Z, std::f64::consts::FRAC_PI_2);
    assert!(result.is_ok(), "revolve should succeed: {:?}", result.err());
    let brep = result.unwrap();
    assert!(!brep.solids.is_empty(), "revolved brep must have solids");
    assert!(face_count(&brep) > 0, "revolved brep must have faces");
}

#[test]
fn sweep_pipe_returns_brep() {
    // Square cross-section swept along a short straight spine.
    let profile = vec![
        DVec2::new(-0.1, -0.1),
        DVec2::new(0.1, -0.1),
        DVec2::new(0.1, 0.1),
        DVec2::new(-0.1, 0.1),
    ];
    let spine: Vec<DVec3> = (0..=8)
        .map(|i| DVec3::new(0.0, 0.0, i as f64 * 0.25))
        .collect();
    let brep = sweep_pipe(&profile, &spine).unwrap();
    assert!(!brep.solids.is_empty(), "swept brep must have solids");
    assert!(face_count(&brep) > 0, "swept brep must have faces");
}

#[test]
fn loft_two_profiles_returns_brep() {
    // Loft between a large and small regular polygon.
    let n = 8;
    let profile1: Vec<DVec3> = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            DVec3::new(a.cos(), a.sin(), 0.0)
        })
        .collect();
    let profile2: Vec<DVec3> = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            DVec3::new(0.5 * a.cos(), 0.5 * a.sin(), 2.0)
        })
        .collect();
    let brep = loft(&[profile1, profile2]).unwrap();
    assert!(!brep.solids.is_empty(), "lofted brep must have solids");
    assert!(face_count(&brep) > 0, "lofted brep must have faces");
}

// ─────────────────────────────────────────────────────────────────────────────
// Blend operations
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn fillet_box_edge_returns_brep() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let result = fillet_edge(&brep, 0, 0.1).unwrap();
    assert!(
        face_count(&result) >= 6,
        "filleted box must have at least 6 faces, got {}",
        face_count(&result)
    );
    assert!(!result.solids.is_empty(), "filleted box must have solids");
}

#[test]
fn chamfer_box_edge_returns_brep() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let result = chamfer_edge(&brep, 0, 0.1).unwrap();
    assert!(
        face_count(&result) >= 6,
        "chamfered box must have at least 6 faces, got {}",
        face_count(&result)
    );
    assert!(!result.solids.is_empty(), "chamfered box must have solids");
}

// ─────────────────────────────────────────────────────────────────────────────
// Error paths
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn invalid_radius_returns_error() {
    let result = make_sphere_brep(DVec3::ZERO, -1.0);
    assert!(result.is_err(), "negative radius should fail: {:?}", result.ok());
}

#[test]
fn zero_radius_sphere_returns_error() {
    let result = make_sphere_brep(DVec3::ZERO, 0.0);
    assert!(result.is_err(), "zero radius should fail: {:?}", result.ok());
}

#[test]
fn non_finite_dimension_returns_error() {
    let result = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, f64::NAN, 1.0, 1.0);
    assert!(result.is_err(), "NaN dimension should fail: {:?}", result.ok());
}

#[test]
fn zero_height_cylinder_returns_error() {
    let result = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 0.0);
    assert!(result.is_err(), "zero height should fail: {:?}", result.ok());
}
