// 3-D mesher quality regression baselines (cube / stretched box / slender
// box across Delaunay3D / Frontal3D / Hxt3D). The quality statistics and
// the box geometry are the corpus module's single source (the local copies
// this file used to carry were deduplicated into `rmsh_algo::corpus`).
//
// These bounds are the raw-mesher regression baselines (measured levels
// with headroom, NOT the ROADMAP post-optimization hard gate — that gate
// is carried by `corpus::QualityGate::roadmap_hard_gate` and reported as
// a distance in the corpus tests; the quality card owns the gap).

use rmsh_algo::{
    box_surface, format_tet_stats, tet_quality_stats, Delaunay3D, Frontal3D, Hxt3D, MeshParams,
    Mesher3D, TetQualityStats,
};

fn assert_quality_baseline(name: &str, stats: TetQualityStats) {
    eprintln!("{}", format_tet_stats(name, &stats));

    assert!(stats.tet_count > 0, "{} produced no tets", name);
    assert!(stats.min_dihedral_deg > 0.10, "{} min dihedral too low", name);
    assert!(stats.p95_radius_edge < 120.0, "{} radius-edge exploded", name);
    assert!(stats.sliver_fraction < 0.40, "{} too many slivers", name);
}

#[test]
fn mesher3d_quality_baseline_cube() {
    let surface = box_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let params = MeshParams::with_size(0.4);

    let d = Delaunay3D::default().mesh_3d(&surface, &params).unwrap();
    let f = Frontal3D::default().mesh_3d(&surface, &params).unwrap();
    let h = Hxt3D::default().mesh_3d(&surface, &params).unwrap();

    assert_quality_baseline("delaunay3d", tet_quality_stats(&d));
    assert_quality_baseline("frontal3d", tet_quality_stats(&f));
    assert_quality_baseline("hxt3d", tet_quality_stats(&h));
}

#[test]
fn mesher3d_quality_baseline_stretched_box() {
    let surface = box_surface([0.0, 0.0, 0.0], [1.8, 1.1, 0.9]);
    let params = MeshParams::with_size(0.38);

    let d = Delaunay3D::default().mesh_3d(&surface, &params).unwrap();
    let f = Frontal3D::default().mesh_3d(&surface, &params).unwrap();
    let h = Hxt3D::default().mesh_3d(&surface, &params).unwrap();

    let qd = tet_quality_stats(&d);
    let qf = tet_quality_stats(&f);
    let qh = tet_quality_stats(&h);

    eprintln!("{}", format_tet_stats("delaunay3d_stretched", &qd));
    eprintln!("{}", format_tet_stats("frontal3d_stretched", &qf));
    eprintln!("{}", format_tet_stats("hxt3d_stretched", &qh));

    assert!(qd.tet_count > 0 && qf.tet_count > 0 && qh.tet_count > 0);
    assert!(qd.min_dihedral_deg.is_finite() && qf.min_dihedral_deg.is_finite() && qh.min_dihedral_deg.is_finite());
    assert!(qd.p95_radius_edge.is_finite() && qf.p95_radius_edge.is_finite() && qh.p95_radius_edge.is_finite());
    assert!((0.0..=1.0).contains(&qd.sliver_fraction));
    assert!((0.0..=1.0).contains(&qf.sliver_fraction));
    assert!((0.0..=1.0).contains(&qh.sliver_fraction));

    assert!(qf.min_dihedral_deg >= qd.min_dihedral_deg * 0.95);
    let p95_ratio = if qd.p95_radius_edge > qf.p95_radius_edge {
        qd.p95_radius_edge / qf.p95_radius_edge
    } else {
        qf.p95_radius_edge / qd.p95_radius_edge
    };
    assert!(p95_ratio <= 1.15);
    assert!(qf.sliver_fraction <= qd.sliver_fraction + 0.05);
}

#[test]
fn mesher3d_quality_slender_box_edge_pressure() {
    let surface = box_surface([0.0, 0.0, 0.0], [3.2, 0.55, 0.45]);
    let params = MeshParams::with_size(0.26);

    let d = Delaunay3D::default().mesh_3d(&surface, &params).unwrap();
    let f = Frontal3D::default().mesh_3d(&surface, &params).unwrap();
    let h = Hxt3D::default().mesh_3d(&surface, &params).unwrap();

    let qd = tet_quality_stats(&d);
    let qf = tet_quality_stats(&f);
    let qh = tet_quality_stats(&h);

    eprintln!("{}", format_tet_stats("delaunay3d_slender", &qd));
    eprintln!("{}", format_tet_stats("frontal3d_slender", &qf));
    eprintln!("{}", format_tet_stats("hxt3d_slender", &qh));

    assert!(qd.tet_count > 0 && qf.tet_count > 0 && qh.tet_count > 0);
    assert!(qd.min_dihedral_deg > 0.05 && qf.min_dihedral_deg > 0.05 && qh.min_dihedral_deg > 0.05);
    assert!(qd.p95_radius_edge.is_finite() && qf.p95_radius_edge.is_finite() && qh.p95_radius_edge.is_finite());
    assert!(qd.p95_radius_edge < 80000.0 && qf.p95_radius_edge < 80000.0 && qh.p95_radius_edge < 80000.0);
    assert!(qd.sliver_fraction < 0.60 && qf.sliver_fraction < 0.60 && qh.sliver_fraction < 0.60);
}
