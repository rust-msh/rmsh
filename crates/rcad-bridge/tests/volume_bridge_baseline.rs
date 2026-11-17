//! The tracked baseline of the general BRep -> volume-mesh half bridge.
//!
//! The corpus (geometry + mesher prescription) lives in this file; the
//! measured reference values live in `tests/baseline/volume_bridge.json`,
//! tracked in git. Three things are pinned per case:
//!
//! 1. **Analytic volume conservation** — the unit box is exactly 1.0 (its
//!    faces are planar, so the tetrahedralization is exact); the sphere and
//!    cylinder carry a tolerance because their boundary is a chordal
//!    triangulation of a curved surface.
//! 2. **Exact tetrahedron count** — the mesher is deterministic, so the
//!    connectivity is an exact regression pin.
//! 3. **Quality bands** — `inspect_volume_mesh` statistics must stay within a
//!    generous band of the measured level, so a mesher change that degrades
//!    (or grossly moves) quality owes the review an explanation without
//!    making the pin brittle to last-bit platform noise.
//!
//! This is the rmsh half only (no CAE contract): the projection into
//! `rcad_cae::VolumeMesh` is the domain consumer's step.
//!
//! Re-bless after an intentional mesher change:
//!
//! ```bash
//! RMSH_BRIDGE_BLESS=1 cargo test -p rmsh-bridge --test volume_bridge_baseline
//! ```
//!
//! # Tessellation note
//!
//! The boundary is mounted by meshing the **solid root only**. Meshing the
//! solid and then each face again (the box fixtures' habit) re-meshes faces
//! whose `active_triangulation` index has already been rewritten, which
//! panics the cylinder's pcurve reader; the solid root meshes every face
//! exactly once.

use std::sync::{Arc, RwLock};

use glam::DVec3;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, TShape};
use rcad_modeling::{make_box_brep, make_cylinder_brep, make_sphere_brep};
use rcad_tess::b_rep_mesh::facade::IncrementalMesh;
use rmsh_algo::MeshParams;
use rmsh_bridge::{VolumeMeshReport, brep_to_volume_mesh, inspect_volume_mesh};
use serde_json::{Value, json};

/// Path of the tracked baseline, relative to this test file.
const BASELINE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/baseline/volume_bridge.json");

/// One corpus case: a BRep with an analytically known volume, the mesher
/// prescription, and the tessellation deflections that mount its boundary.
struct Case {
    name: &'static str,
    brep: BRep,
    params: MeshParams,
    /// Linear deflection of the boundary triangulation.
    deflection: f64,
    /// Angular deflection of the boundary triangulation (radians). For a
    /// curved surface this dominates the chordal volume deficit, so the
    /// sphere and cylinder tighten it.
    angular_deflection: f64,
    analytic_volume: f64,
    /// Relative tolerance of the mesh volume against the analytic value.
    volume_rel_tol: f64,
}

fn corpus() -> Vec<Case> {
    vec![
        Case {
            name: "box_unit",
            brep: make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
                .expect("the unit box builds"),
            params: MeshParams::with_size(0.4),
            deflection: 0.1,
            angular_deflection: 0.5,
            analytic_volume: 1.0,
            volume_rel_tol: 1.0e-12,
        },
        Case {
            name: "sphere_r1",
            brep: make_sphere_brep(DVec3::ZERO, 1.0).expect("the sphere builds"),
            params: MeshParams::with_size(0.5),
            deflection: 0.1,
            angular_deflection: 0.3,
            analytic_volume: 4.0 * std::f64::consts::PI / 3.0,
            volume_rel_tol: 2.0e-2,
        },
        Case {
            name: "cylinder_r1_h2",
            brep: make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 2.0)
                .expect("the cylinder builds"),
            params: MeshParams::with_size(0.5),
            deflection: 0.1,
            angular_deflection: 0.8,
            analytic_volume: 2.0 * std::f64::consts::PI,
            volume_rel_tol: 3.0e-2,
        },
    ]
}

/// Mount a triangulation on every top-level solid, exactly as a caller does
/// before bridging, and return the shared BRep.
fn tessellate(brep: &BRep, deflection: f64, angular_deflection: f64) -> Arc<RwLock<BRep>> {
    let shared = Arc::new(RwLock::new(brep.clone()));
    let roots: Vec<Shape> = {
        let guard = shared.read().unwrap();
        guard
            .tshapes
            .iter()
            .enumerate()
            .filter(|(_, t)| matches!(t.as_ref(), TShape::Solid(_)))
            .map(|(i, t)| Shape::from_parts(t.clone(), i, 0, Orientation::Forward))
            .collect()
    };
    for root in roots {
        let mesher = IncrementalMesh::perform_with_deflection(
            Arc::clone(&shared),
            &root,
            deflection,
            false,
            angular_deflection,
            false,
        );
        assert!(mesher.is_done(), "IncrementalMesh::Perform is done");
    }
    shared
}

/// Bridge one case and return the mesh plus its no-solver report.
fn bridge_case(case: &Case) -> (rmsh_model::Mesh, VolumeMeshReport) {
    let shared = tessellate(&case.brep, case.deflection, case.angular_deflection);
    let brep = shared.read().unwrap();
    let mesh = brep_to_volume_mesh(&brep, &case.params)
        .unwrap_or_else(|e| panic!("{} must bridge: {e}", case.name));
    let report = inspect_volume_mesh(&mesh);
    (mesh, report)
}

/// The measured record written to / read from the tracked baseline.
fn record(report: &VolumeMeshReport) -> Value {
    json!({
        "volume": report.volume,
        "tet_count": report.tet_count,
        "min_dihedral_deg": report.min_dihedral_deg,
        "p95_radius_edge": report.p95_radius_edge,
        "sliver_fraction": report.sliver_fraction,
    })
}

/// The allowed drift of a quality field: generous, so last-bit platform
/// noise does not fail the pin, but a real mesher regression does.
fn quality_band(field: &str, want: f64) -> f64 {
    match field {
        "min_dihedral_deg" => (0.5 * want).max(0.05),
        "p95_radius_edge" => (0.5 * want).max(0.5),
        _ => 0.02,
    }
}

#[test]
fn the_half_bridge_conserves_analytic_volumes_and_holds_quality() {
    let cases = corpus();
    let mut measured = serde_json::Map::new();
    let mut checked = Vec::new();

    for case in &cases {
        let (_mesh, report) = bridge_case(case);

        // Analytic volume conservation.
        let error = report.volume_error(case.analytic_volume);
        assert!(
            error <= case.volume_rel_tol,
            "{}: volume {} vs analytic {} (rel {}, tol {})",
            case.name,
            report.volume,
            case.analytic_volume,
            error,
            case.volume_rel_tol
        );

        measured.insert(case.name.to_string(), record(&report));
        checked.push((case, report));
    }

    let bless = std::env::var("RMSH_BRIDGE_BLESS").ok().as_deref() == Some("1");
    let baseline: Value = if bless {
        let document = json!({
            "version": 1,
            "note": "Measured rmsh BRep -> volume-mesh half bridge (rmsh_algo::Delaunay3D). Volume is also checked against the analytic value in the test; tet_count is an exact pin and the quality fields are regression bands. Re-bless with RMSH_BRIDGE_BLESS=1.",
            "notes": {
                "box_unit": "Planar faces, so the tetrahedralization conserves the volume exactly.",
                "sphere_r1": "The lat-long boundary triangulation's pole fan limits min_dihedral; the mesher respects the boundary, so the pole tets stay thin. The angular deflection is the accuracy/quality knob.",
                "cylinder_r1_h2": "A coarse angular deflection (0.8 rad) trades volume accuracy for a sliver-free mesh; the caps are planar and the lateral surface chordal.",
            },
            "cases": Value::Object(measured.clone()),
        });
        std::fs::write(
            BASELINE,
            serde_json::to_string_pretty(&document).expect("the baseline serializes"),
        )
        .expect("the baseline writes");
        document
    } else {
        serde_json::from_str(&std::fs::read_to_string(BASELINE).expect("the baseline reads"))
            .expect("the baseline parses")
    };

    let recorded = baseline["cases"].as_object().expect("cases object");
    for (case, report) in &checked {
        let pin = recorded
            .get(case.name)
            .unwrap_or_else(|| panic!("{} is missing from the baseline", case.name));
        assert_eq!(
            report.tet_count,
            pin["tet_count"].as_u64().expect("tet_count") as usize,
            "{}: tet count drifted (re-bless with RMSH_BRIDGE_BLESS=1)",
            case.name
        );
        assert!(
            (report.volume - pin["volume"].as_f64().expect("volume")).abs()
                <= 1.0e-9 * pin["volume"].as_f64().expect("volume"),
            "{}: volume drifted from {} to {} (re-bless with RMSH_BRIDGE_BLESS=1)",
            case.name,
            pin["volume"],
            report.volume
        );
        for field in ["min_dihedral_deg", "p95_radius_edge", "sliver_fraction"] {
            let want = pin[field].as_f64().expect(field);
            let got = match field {
                "min_dihedral_deg" => report.min_dihedral_deg,
                "p95_radius_edge" => report.p95_radius_edge,
                _ => report.sliver_fraction,
            };
            let band = quality_band(field, want);
            assert!(
                (got - want).abs() <= band,
                "{}: {} drifted from {want} to {got} (band {band}; re-bless with RMSH_BRIDGE_BLESS=1)",
                case.name,
                field
            );
        }
    }
}
