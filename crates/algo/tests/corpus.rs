// P6: the golden corpus pins.
//
// One test per scenario family (see rmsh_algo::corpus module docs for
// the family table). The meshers are deterministic: element counts,
// areas and quality statistics are exact and machine-independent, so
// the pins here are the corpus goldens — a mesher change that moves
// them owes the review an explanation. Time is NOT pinned here; the
// `corpus_benchmark` harness times the same families explicitly.

use rmsh_algo::{
    box_surface, format_tet_stats, format_tri_stats, l_shape_polygon, mesh_polygon,
    nested_cubes_surfaces, plate_with_hole_domain, rim_distance, tet_quality_stats,
    tri_quality_stats, AdaptiveConfig, AdaptiveLoop2D, CurvatureField, Delaunay3D, Frontal3D,
    FrontalDelaunay2D, Hxt3D, MATERIAL_TAG_CAVITY_BASE, MATERIAL_TAG_OUTER, MeshParams, Mesher2D,
    Mesher3D, Polygon2D, QualityGate, RimGeometrySource, TetQualityStats,
};
use rmsh_model::ElementType;

#[allow(unused_imports)]
use rmsh_algo::tetrahedralize_multi_material;

fn interior_edge_stats(
    mesh: &rmsh_model::Mesh,
) -> (Vec<f64>, Vec<f64>) {
    // (rim-band edge lengths |x| within 0.3 of the unit rim, far-field
    // interior edge lengths) — the graded family's locality evidence
    let pos = |id: u64| {
        let n = mesh.nodes.get(&id).expect("node");
        [n.position.x, n.position.y]
    };
    let mut count: std::collections::HashMap<([u64; 2], [u64; 2]), usize> =
        std::collections::HashMap::new();
    let mut mids: std::collections::HashMap<([u64; 2], [u64; 2]), [f64; 2]> =
        std::collections::HashMap::new();
    for elem in &mesh.elements {
        if elem.etype != ElementType::Triangle3 || elem.node_ids.len() != 3 {
            continue;
        }
        let ps = [
            pos(elem.node_ids[0]),
            pos(elem.node_ids[1]),
            pos(elem.node_ids[2]),
        ];
        for k in 0..3 {
            let a = ps[k];
            let b = ps[(k + 1) % 3];
            let ka = [a[0].to_bits(), a[1].to_bits()];
            let kb = [b[0].to_bits(), b[1].to_bits()];
            let (key, mid) = if ka <= kb {
                ((ka, kb), [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5])
            } else {
                ((kb, ka), [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5])
            };
            *count.entry(key).or_insert(0) += 1;
            mids.insert(key, mid);
        }
    }
    let mut rim = Vec::new();
    let mut far = Vec::new();
    for (key, &c) in &count {
        if c != 2 {
            continue;
        }
        let mid = mids[key];
        let len = {
            let (a, b) = key;
            let pa = [f64::from_bits(a[0]), f64::from_bits(a[1])];
            let pb = [f64::from_bits(b[0]), f64::from_bits(b[1])];
            ((pa[0] - pb[0]).powi(2) + (pa[1] - pb[1]).powi(2)).sqrt()
        };
        if rim_distance(mid) <= 0.3 {
            rim.push(len);
        } else {
            far.push(len);
        }
    }
    (rim, far)
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len().max(1) as f64
}

// ── F1: plate with a hole, uniform 0.5 ────────────────────────────────────────

#[test]
fn corpus_plate_with_hole_uniform() {
    let domain = plate_with_hole_domain();
    let mesh = FrontalDelaunay2D::new()
        .mesh_2d(&domain, &MeshParams::with_size(0.5))
        .expect("uniform plate meshes");
    let stats = tri_quality_stats(&mesh);
    println!("{}", format_tri_stats("plate_uniform", &stats));

    // analytic area: 10x10 square minus the r=1 PENTAGON hole (the hole
    // ring is a regular 5-gon: (5/2)·r²·sin(2π/5))
    let pentagon = 2.5 * (2.0 * std::f64::consts::PI / 5.0).sin();
    let expected = 100.0 - pentagon;
    assert!(
        (stats.total_area - expected).abs() < 1e-9,
        "area {} vs expected {expected}",
        stats.total_area
    );
    // golden count (deterministic mesher) — a drift owes the review an
    // explanation
    assert_eq!(
        stats.triangle_count, 1155,
        "golden triangle count drifted"
    );
    // 2-D quality band (the P3a Gmsh x1.15 oracle pins the p95 at the
    // acceptance level for the acceptance mesh; this corpus pin guards
    // the corpus mesh's own band)
    assert!(
        stats.min_angle_deg >= 30.0,
        "min angle {} below the corpus floor",
        stats.min_angle_deg
    );
    assert!(
        stats.p95_radius_edge <= 0.9,
        "p95 radius-edge {} above the corpus band",
        stats.p95_radius_edge
    );
}

// ── F2: plate with a hole, graded curvature field (the P3a matched field) ────

#[test]
fn corpus_plate_with_hole_graded() {
    let rim = RimGeometrySource {
        radius: 1.0,
        delta: 0.3,
    };
    let field = CurvatureField::with_default_factor(&rim, 0.05, 2.0)
        .expect("graded field parameters are valid");
    let params = MeshParams {
        element_size: 2.0,
        min_size: 0.05,
        max_size: 2.0,
        optimize_passes: 3,
    };
    let domain = plate_with_hole_domain();
    let mesh = FrontalDelaunay2D::new()
        .with_size_field(&field)
        .mesh_2d(&domain, &params)
        .expect("graded plate meshes");
    let stats = tri_quality_stats(&mesh);
    println!("{}", format_tri_stats("plate_graded", &stats));

    let expected = 100.0 - 2.5 * (2.0 * std::f64::consts::PI / 5.0).sin();
    assert!((stats.total_area - expected).abs() < 1e-9);
    assert_eq!(stats.triangle_count, 95, "golden triangle count drifted");
    // the P3a acceptance run measured p95 = 0.7582 on this scenario
    // against Gmsh's 0.6840 — the corpus band sits above both with margin
    assert!(
        stats.p95_radius_edge <= 1.0,
        "p95 radius-edge {} above the graded corpus band",
        stats.p95_radius_edge
    );
    assert!(
        stats.min_angle_deg >= 30.0,
        "min angle {} below the corpus floor",
        stats.min_angle_deg
    );
    // LOCALITY: the rim band's mean interior edge is materially finer
    // than the far field's. (Measured 0.644 on this prescription — the
    // graded field's rim size is 1.0 against a far field of 1.8, so the
    // honest corpus bar is 0.75; the aggressive x0.5 locality contract
    // belongs to the adaptive loop's sharp-front pin.)
    let (rim_edges, far_edges) = interior_edge_stats(&mesh);
    let (rim_mean, far_mean) = (mean(&rim_edges), mean(&far_edges));
    println!(
        "plate_graded locality: rim mean edge {rim_mean:.4} ({} edges), far mean {far_mean:.4} ({} edges)",
        rim_edges.len(),
        far_edges.len()
    );
    assert!(
        rim_mean < far_mean * 0.75,
        "rim band must be materially finer: {rim_mean} vs {far_mean}"
    );
}

// ── F3: concave L shape through the CDT path ─────────────────────────────────

#[test]
fn corpus_l_shape_concave() {
    let poly = Polygon2D::new(l_shape_polygon());
    let mesh = mesh_polygon(&poly, 0.4).expect("L shape meshes");
    let stats = tri_quality_stats(&mesh);
    println!("{}", format_tri_stats("l_shape", &stats));

    // THE P3 pin: the exact polygon area (no out-of-domain leakage)
    assert!(
        (stats.total_area - 3.0).abs() < 1e-9,
        "L-shape area {} must be exactly 3.0",
        stats.total_area
    );
    assert_eq!(stats.triangle_count, 44, "golden triangle count drifted");
    assert!(
        stats.min_angle_deg >= 20.0,
        "min angle {} below the corpus floor",
        stats.min_angle_deg
    );
    assert!(
        stats.p95_radius_edge <= 1.3,
        "p95 radius-edge {} above the corpus band",
        stats.p95_radius_edge
    );
}

// ── F4: unit cube, tetrahedral, three meshers ────────────────────────────────

fn corpus_cube_gate(name: &str, stats: &TetQualityStats) {
    println!("{}", format_tet_stats(name, &stats));
    assert!(stats.tet_count > 0, "{name} produced no tets");
    // Regression bound at the measured levels (re-pinned by the 3-D
    // quality card: wiring the flip-convergence arm into Delaunay3D
    // moved the cube family from p95 20.6 / sliver 3.3% to p95 8.9 /
    // sliver 2.6%; stretched + slender stay at their boundary-locked
    // levels — the spindle residue the local flips cannot touch). The
    // ROADMAP hard gate is the post-optimization target and is
    // REPORTED here, not asserted — the remaining gap is bounded by
    // the corner-diagonal spindle family (see ROADMAP P4).
    let hard = QualityGate::roadmap_hard_gate();
    for v in hard.violations(stats) {
        println!("  [hard-gate distance] {name}: {v}");
    }
    assert!(
        stats.min_dihedral_deg > 5.0,
        "{name} min dihedral too low: {}",
        stats.min_dihedral_deg
    );
    assert!(
        stats.p95_radius_edge <= 10.0,
        "{name} radius-edge above the measured bound: {}",
        stats.p95_radius_edge
    );
    assert!(
        stats.sliver_fraction < 0.035,
        "{name} too many slivers: {}",
        stats.sliver_fraction
    );
}

#[test]
fn corpus_cube_tet_three_meshers() {
    let surface = box_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let params = MeshParams::with_size(0.4);

    let d = Delaunay3D::default().mesh_3d(&surface, &params).unwrap();
    let f = Frontal3D::default().mesh_3d(&surface, &params).unwrap();
    let h = Hxt3D::default().mesh_3d(&surface, &params).unwrap();

    let qd = tet_quality_stats(&d);
    let qf = tet_quality_stats(&f);
    let qh = tet_quality_stats(&h);
    corpus_cube_gate("cube_delaunay3d", &qd);
    corpus_cube_gate("cube_frontal3d", &qf);
    corpus_cube_gate("cube_hxt3d", &qh);

    // Golden exact pins — all THREE meshers. The Hxt3D "cross-binary
    // drift" ledger entry resolved as an input-confusion artifact (the
    // hxt-test box fixture has an inward bottom face, the corpus box is
    // outward; 19.777 vs 21.944 were the two DIFFERENT fixtures), and
    // five-process digest equality over the whole fixture family pins
    // the pipeline cross-process bit-deterministic (2026-10-14). A
    // drift here owes the review an explanation.
    assert_eq!(qd.tet_count, 76, "delaunay3d golden count drifted");
    assert_eq!(qf.tet_count, 76, "frontal3d golden count drifted");
    assert_eq!(qh.tet_count, 73, "hxt3d golden count drifted");
    for (name, q, min_d, p95, sliver) in [
        (
            "delaunay3d",
            &qd,
            5.624_999_999_999_975,
            8.855_231_606_230_609,
            0.026_315_789_473_684_21,
        ),
        (
            "frontal3d",
            &qf,
            5.624_999_999_999_975,
            8.855_231_606_230_609,
            0.026_315_789_473_684_21,
        ),
        (
            "hxt3d",
            &qh,
            21.751_817_807_142_02,
            4.586_175_512_861_971,
            0.0,
        ),
    ] {
        assert!(
            (q.min_dihedral_deg - min_d).abs() < 1e-6,
            "{name} golden min dihedral drifted: {}",
            q.min_dihedral_deg
        );
        assert!(
            (q.p95_radius_edge - p95).abs() < 1e-6,
            "{name} golden p95 radius-edge drifted: {}",
            q.p95_radius_edge
        );
        assert!(
            (q.sliver_fraction - sliver).abs() < 1e-9,
            "{name} golden sliver fraction drifted: {}",
            q.sliver_fraction
        );
    }
}

// ── F5: nested cubes, multi-material ─────────────────────────────────────────

#[test]
fn corpus_nested_cubes_multi_material() {
    let (outer, inner) = nested_cubes_surfaces();
    let out = tetrahedralize_multi_material(&outer, &[inner]).expect("nested cubes tetrahedralize");

    let vol = |tag: i32| -> f64 {
        let mut total = 0.0;
        for e in &out.elements {
            if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
                continue;
            }
            if e.physical_tag != Some(tag) {
                continue;
            }
            let p = |i: usize| {
                let n = out.nodes.get(&e.node_ids[i]).expect("node");
                [n.position.x, n.position.y, n.position.z]
            };
            let (a, b, c, d) = (p(0), p(1), p(2), p(3));
            let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let ad = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
            let cr = [
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0],
            ];
            total += (ad[0] * cr[0] + ad[1] * cr[1] + ad[2] * cr[2]).abs() / 6.0;
        }
        total
    };
    let vol_outer = vol(MATERIAL_TAG_OUTER);
    let vol_inner = vol(MATERIAL_TAG_CAVITY_BASE);
    assert!((vol_outer - 0.578125).abs() < 1e-9, "outer {vol_outer}");
    assert!((vol_inner - 0.421875).abs() < 1e-9, "inner {vol_inner}");
    let stats = tet_quality_stats(&out);
    println!("{}", format_tet_stats("nested_cubes", &stats));
    assert!(stats.tet_count > 0);
}

// ── F6: sharp-front adaptive loop (corpus budget, not the contract pins) ─────
//
// The loop CONTRACTS (band concentration, determinism, clamp honor) are
// pinned in tests/adaptive_loop.rs; this corpus entry pins the loop's
// place in the corpus: it runs on the corpus domain and concentrates
// elements in the front band. The budget is kept light so the default
// suite stays fast.

#[test]
fn corpus_sharp_front_adaptive() {
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
        .expect("adaptive loop runs");
    assert!(report.steps.len() >= 2, "at least one refinement step");
    let first = report.steps.first().expect("first").elements;
    let last = report.steps.last().expect("last").elements;
    assert!(last > first, "elements {first} -> {last} must grow");

    // band concentration on the final mesh
    let mesh = &report.mesh;
    let mut band = 0usize;
    let mut total = 0usize;
    for e in &mesh.elements {
        if e.etype != ElementType::Triangle3 {
            continue;
        }
        total += 1;
        let cx = e
            .node_ids
            .iter()
            .map(|id| mesh.nodes.get(id).expect("node").position.x)
            .sum::<f64>()
            / 3.0;
        if (cx - 0.5).abs() < 0.15 {
            band += 1;
        }
    }
    println!("sharp_front_adaptive: steps={}, {total} tris, {band} in band", report.steps.len());
    assert!(band * 4 > total, "the band must concentrate elements");
}
