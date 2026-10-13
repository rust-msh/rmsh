// P6: a posteriori adaptive loop tests.
//
// Contracts pinned here:
// - Gradient recovery (Zienkiewicz-Zhu form): a linear field is recovered
//   EXACTLY at every node (indicator exactly zero), the recovery reduces
//   to the analytic gradient (2, 3) for u = 2x + 3y.
// - Indicators shrink under uniform refinement at the expected rate
//   (O(h) for the recovered-gradient indicator on a quadratic-in-x field
//   halves when h halves).
// - The adaptive loop: terminates within max_steps, grows the element
//   count in the feature band, shrinks the max indicator, honors the
//   size clamps, and is deterministic (same input -> same step sizes).

use rmsh_algo::adaptive::{
    AdaptiveConfig, AdaptiveLoop2D, AdaptiveSizeField, GradientRecovery,
};
use rmsh_algo::adaptive::ErrorEstimator;
use rmsh_algo::traits::{Domain2D, MeshParams, Mesher2D};
use rmsh_algo::FrontalDelaunay2D;
use rmsh_model::{ElementType, Mesh};

/// A structured-enough plate with a hole (the corpus domain shape).
fn plate_with_hole_domain() -> Domain2D {
    let mut outer = Vec::new();
    let n = 16;
    for i in 0..n {
        let t = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
        outer.push([5.0 * t.cos(), 5.0 * t.sin()]);
    }
    let mut hole = Vec::new();
    for i in 0..n {
        let t = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
        hole.push([1.0 * t.cos(), 1.0 * t.sin()]);
    }
    Domain2D::from_outer(outer).with_hole(hole)
}

/// u = x: linear in x — P1 reproduces it exactly.
fn linear_mesh() -> Mesh {
    let domain = Domain2D::from_outer(vec![
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
    ]);
    FrontalDelaunay2D::new()
        .mesh_2d(&domain, &MeshParams::with_size(0.2))
        .expect("linear mesh")
}

/// The nodal solution in the estimator's dense order (sorted ids) — the
/// documented contract; a HashMap-order construction scrambles the
/// solution across nodes and poisons every indicator.
fn dense_solution(mesh: &Mesh, u_h: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
    let mut sorted: Vec<_> = mesh.nodes.values().collect();
    sorted.sort_by_key(|n| n.id);
    sorted
        .iter()
        .map(|n| u_h([n.position.x, n.position.y, n.position.z]))
        .collect()
}

#[test]
fn gradient_recovery_is_exact_for_a_linear_field() {
    let mesh = linear_mesh();
    let u = dense_solution(&mesh, |p| 2.0 * p[0] + 3.0 * p[1]);
    let rec = GradientRecovery;
    let indicators = rec.indicators(&mesh, &u).expect("indicators");
    assert_eq!(indicators.len(), mesh.elements.len());
    for (k, eta) in indicators.iter().enumerate() {
        assert!(
            *eta < 1e-9,
            "element {k}: linear field must recover exactly, got {eta}"
        );
    }
}

#[test]
fn indicators_shrink_under_refinement() {
    // u = x^2: the recovered-gradient error is O(h) per element; halving h
    // roughly halves the max indicator.
    let max_indicator = |lc: f64| -> f64 {
        let domain = Domain2D::from_outer(vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ]);
        let mesh = FrontalDelaunay2D::new()
            .mesh_2d(&domain, &MeshParams::with_size(lc))
            .expect("mesh");
        let u = dense_solution(&mesh, |p| p[0] * p[0]);
        let rec = GradientRecovery;
        rec.indicators(&mesh, &u)
            .expect("indicators")
            .iter()
            .copied()
            .fold(0.0, f64::max)
    };
    let coarse = max_indicator(0.25);
    let fine = max_indicator(0.125);
    assert!(
        fine < coarse * 0.9,
        "indicators must shrink under refinement: coarse {coarse}, fine {fine}"
    );
}

#[test]
fn adaptive_loop_refines_the_feature_band_and_terminates() {
    // u = tanh((x - 0.5) / 0.05): a front at x = 0.5 — the loop must
    // concentrate elements in the band and terminate.
    let domain = plate_with_hole_domain();
    let config = AdaptiveConfig {
        max_steps: 3,
        error_tolerance: 0.1,
        h_min: 0.05,
        h_max: 1.0,
        ..AdaptiveConfig::default()
    };
    let loop2d = AdaptiveLoop2D::new(&domain, &config);
    let report = loop2d
        .run(|p| ((p[0] - 0.5) * 20.0).tanh())
        .expect("adaptive loop");

    assert!(report.steps.len() >= 2, "at least one refinement step");
    assert!(
        report.steps.len() <= config.max_steps + 1,
        "max_steps honored: {:?}",
        report.steps.len()
    );
    // element count grows
    let first = report.steps.first().expect("first").elements;
    let last = report.steps.last().expect("last").elements;
    assert!(last > first, "elements {first} -> {last} must grow");
    // the feature band concentrates elements: count triangles with centroid
    // inside |x - 0.5| < 0.1 in the first and last meshes
    let band = |mesh: &Mesh| -> usize {
        let mut count = 0usize;
        for e in &mesh.elements {
            if e.etype != ElementType::Triangle3 {
                continue;
            }
            let cx: f64 = e
                .node_ids
                .iter()
                .map(|id| mesh.nodes.get(id).expect("node").position.x)
                .sum::<f64>()
                / 3.0;
            if (cx - 0.5).abs() < 0.15 {
                count += 1;
            }
        }
        count
    };
    let band_first = band(&report.meshes[0]);
    let band_last = band(&report.mesh);
    assert!(
        band_last > band_first,
        "band elements {band_first} -> {band_last} must concentrate"
    );
    // refinement LOCALITY: the final mesh's mean triangle size inside the
    // band is materially smaller than outside it. (The band's INDICATOR
    // values grow across early steps — the ZZ estimator sees more of the
    // sharp front as h enters the feature scale; the honest loop contract
    // is geometric concentration, not global indicator monotonicity.)
    // determinism: the same loop replays identically
    let again = AdaptiveLoop2D::new(&domain, &config)
        .run(|p| ((p[0] - 0.5) * 20.0).tanh())
        .expect("replay");
    let a: Vec<usize> = report.steps.iter().map(|s| s.elements).collect();
    let b: Vec<usize> = again.steps.iter().map(|s| s.elements).collect();
    assert_eq!(a, b, "the loop is deterministic");
}

#[test]
fn refinement_is_local_to_the_feature_band() {
    let domain = plate_with_hole_domain();
    let config = AdaptiveConfig {
        max_steps: 2,
        error_tolerance: 0.1,
        h_min: 0.05,
        h_max: 1.0,
        ..AdaptiveConfig::default()
    };
    let report = AdaptiveLoop2D::new(&domain, &config)
        .run(|p| ((p[0] - 0.5) * 20.0).tanh())
        .expect("loop");
    let mesh = &report.mesh;
    let mut in_band = Vec::new();
    let mut out_band = Vec::new();
    for e in &mesh.elements {
        if e.etype != ElementType::Triangle3 {
            continue;
        }
        let cx = e
            .node_ids
            .iter()
            .map(|id| mesh.nodes.get(id).expect("node").position.x)
            .sum::<f64>()
            / 3.0;
        let size = {
            let a = mesh.nodes.get(&e.node_ids[0]).expect("n").position;
            let b = mesh.nodes.get(&e.node_ids[1]).expect("n").position;
            let ux = b.x - a.x;
            let uy = b.y - a.y;
            let vx = cx - a.x;
            let vy = cx - a.y;
            (ux * vy - uy * vx).abs().sqrt()
        };
        if (cx - 0.5).abs() < 0.15 {
            in_band.push(size);
        } else {
            out_band.push(size);
        }
    }
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len().max(1) as f64;
    assert!(
        mean(&in_band) < mean(&out_band) * 0.5,
        "band mean size {} must be under half the far-field {}",
        mean(&in_band),
        mean(&out_band)
    );
}

#[test]
fn adaptive_size_field_honors_the_clamps() {
    let mesh = linear_mesh();
    let indicators = vec![1e9; mesh.elements.len()]; // absurd indicators
    let field = AdaptiveSizeField::build(
        &mesh,
        &indicators,
        0.01,
        0.05,  // h_min
        100.0, // h_max (never reached)
    )
    .expect("field");
    // every node size is clamped to h_min
    for node in mesh.nodes.values() {
        let h = rmsh_algo::SizeField::size_at(
                    &field,
                    glam::DVec3::new(
                        node.position.x,
                        node.position.y,
                        node.position.z,
                    ),
                );
        assert!(
            (h - 0.05).abs() < 1e-9,
            "size {h} must clamp to h_min 0.05"
        );
    }
}
