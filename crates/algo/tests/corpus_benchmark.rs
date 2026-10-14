// P6: the public benchmark harness over the golden corpus families.
//
// NOT part of the default suite (`#[ignore]`): run it explicitly with
//
//     CARGO_TARGET_DIR=D:/tmp/rmsh-target cargo test -p rmsh-algo \
//         --test corpus_benchmark -- --ignored --nocapture
//
// What is reproducible: the scenario families, the meshers, the
// deterministic outputs (element counts and quality statistics — pinned
// by tests/corpus.rs). What is machine-relative: TIME. Each family runs
// 3 repetitions and reports the BEST wall-clock time (min-of-3 damps
// scheduler noise); the table is the artifact to publish alongside a
// machine spec. A regression in element counts or quality between runs
// of this harness on the same build is a mesher bug, not a benchmark
// variance.
//
// Comparisons: the P3a Gmsh x1.15 acceptance lines (uniform/graded p95
// radius-edge against a local Gmsh 4.14 CLI) live in the mesher's own
// test module and run through the Gmsh-comparison tests there; this
// harness reports rmsh-only numbers.
//
// # Archived runs
//
// 2026-10-14, release profile, 13th Gen Intel Core i7-13700 (24 logical
// cores), Windows 10.0.26200, min-of-3 — the publishable artifact:
//
// | family | elements | best time |
// |---|---|---|
// | plate_with_hole_uniform (h=0.5) | 1155 | 62.4 ms |
// | plate_with_hole_graded (rim field) | 95 | 2.70 ms |
// | l_shape_concave (h=0.4) | 44 | 624 µs |
// | cube_delaunay3d (h=0.4) | 76 | 5.35 ms |
// | cube_frontal3d (h=0.4) | 76 | 5.24 ms |
// | cube_hxt3d (h=0.4) | 73 | 83.5 ms |
// | sharp_front_adaptive (2 steps) | 2573 | 1.75 s |
//
// Same day, debug profile (the numbers the perf cards quote): uniform
// 160 ms, graded 7.2 ms, l_shape 1.9 ms, cubes 20.3/20.6/295 ms,
// sharp_front_adaptive 9.17 s (12.29 s before the AdaptiveSizeField
// grid index, same session — the Fix-C A/B).

use rmsh_algo::{
    box_surface, format_tet_stats, format_tri_stats, l_shape_polygon, mesh_polygon,
    plate_with_hole_domain, rim_distance, tet_quality_stats, tri_quality_stats, AdaptiveConfig,
    AdaptiveLoop2D, CurvatureField, Delaunay3D, Frontal3D, FrontalDelaunay2D, Hxt3D, MeshParams,
    Mesher2D, Mesher3D, Polygon2D, RimGeometrySource,
};
use rmsh_model::ElementType;
use std::time::{Duration, Instant};

fn bench<F: FnMut() -> usize>(name: &str, reps: usize, mut f: F) {
    let mut best = Duration::MAX;
    let mut last_count = 0usize;
    for _ in 0..reps {
        let t = Instant::now();
        last_count = f();
        best = best.min(t.elapsed());
    }
    println!("| {name} | {last_count} | {best:?} |");
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len().max(1) as f64
}

#[test]
#[ignore = "the public benchmark: run explicitly with --ignored --nocapture"]
fn corpus_benchmark_run() {
    println!();
    println!("# rmsh corpus benchmark (debug build)");
    println!();
    println!("| family | elements | best time (3 reps) |");
    println!("|---|---|---|");

    // F1: plate with a hole, uniform
    bench("plate_with_hole_uniform (h=0.5)", 3, || {
        let domain = plate_with_hole_domain();
        let mesh = FrontalDelaunay2D::new()
            .mesh_2d(&domain, &MeshParams::with_size(0.5))
            .expect("uniform plate");
        let st = tri_quality_stats(&mesh);
        println!("  {}", format_tri_stats("uniform", &st));
        st.triangle_count
    });

    // F2: plate with a hole, graded curvature field
    bench("plate_with_hole_graded (rim field)", 3, || {
        let rim = RimGeometrySource {
            radius: 1.0,
            delta: 0.3,
        };
        let field = CurvatureField::with_default_factor(&rim, 0.05, 2.0).expect("field");
        let domain = plate_with_hole_domain();
        let mesh = FrontalDelaunay2D::new()
            .with_size_field(&field)
            .mesh_2d(
                &domain,
                &MeshParams {
                    element_size: 2.0,
                    min_size: 0.05,
                    max_size: 2.0,
                    optimize_passes: 3,
                },
            )
            .expect("graded plate");
        let st = tri_quality_stats(&mesh);
        println!("  {}", format_tri_stats("graded", &st));
        st.triangle_count
    });

    // F3: concave L shape (CDT path)
    bench("l_shape_concave (h=0.4)", 3, || {
        let poly = Polygon2D::new(l_shape_polygon());
        let mesh = mesh_polygon(&poly, 0.4).expect("L shape");
        let st = tri_quality_stats(&mesh);
        println!("  {}", format_tri_stats("l_shape", &st));
        st.triangle_count
    });

    // F4: unit cube, three 3-D meshers
    bench("cube_delaunay3d (h=0.4)", 3, || {
        let surface = box_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let mesh = Delaunay3D::default()
            .mesh_3d(&surface, &MeshParams::with_size(0.4))
            .expect("cube");
        let st = tet_quality_stats(&mesh);
        println!("  {}", format_tet_stats("delaunay3d", &st));
        st.tet_count
    });
    bench("cube_frontal3d (h=0.4)", 3, || {
        let surface = box_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let mesh = Frontal3D::default()
            .mesh_3d(&surface, &MeshParams::with_size(0.4))
            .expect("cube");
        let st = tet_quality_stats(&mesh);
        println!("  {}", format_tet_stats("frontal3d", &st));
        st.tet_count
    });
    bench("cube_hxt3d (h=0.4)", 3, || {
        let surface = box_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let mesh = Hxt3D::default()
            .mesh_3d(&surface, &MeshParams::with_size(0.4))
            .expect("cube");
        let st = tet_quality_stats(&mesh);
        println!("  {}", format_tet_stats("hxt3d", &st));
        st.tet_count
    });

    // F6: sharp-front adaptive loop (the corpus budget)
    bench("sharp_front_adaptive (2 steps)", 3, || {
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
            .expect("adaptive loop");
        report
            .steps
            .last()
            .map(|s| s.elements)
            .unwrap_or(0)
    });

    // locality appendix for the graded family (edge-length evidence)
    {
        let rim = RimGeometrySource {
            radius: 1.0,
            delta: 0.3,
        };
        let field = CurvatureField::with_default_factor(&rim, 0.05, 2.0).expect("field");
        let domain = plate_with_hole_domain();
        let mesh = FrontalDelaunay2D::new()
            .with_size_field(&field)
            .mesh_2d(
                &domain,
                &MeshParams {
                    element_size: 2.0,
                    min_size: 0.05,
                    max_size: 2.0,
                    optimize_passes: 3,
                },
            )
            .expect("graded plate");
        let (mut rim_edges, mut far_edges) = (Vec::new(), Vec::new());
        let pos = |id: u64| {
            let n = mesh.nodes.get(&id).expect("node");
            [n.position.x, n.position.y]
        };
        let mut count: std::collections::HashMap<([u64; 2], [u64; 2]), usize> =
            std::collections::HashMap::new();
        for e in &mesh.elements {
            if e.etype != ElementType::Triangle3 {
                continue;
            }
            let ps = [
                pos(e.node_ids[0]),
                pos(e.node_ids[1]),
                pos(e.node_ids[2]),
            ];
            for k in 0..3 {
                let (a, b) = (ps[k], ps[(k + 1) % 3]);
                let ka = [a[0].to_bits(), a[1].to_bits()];
                let kb = [b[0].to_bits(), b[1].to_bits()];
                let key = if ka <= kb { (ka, kb) } else { (kb, ka) };
                *count.entry(key).or_insert(0) += 1;
            }
        }
        for ((ka, kb), &c) in &count {
            if c != 2 {
                continue;
            }
            let pa = [f64::from_bits(ka[0]), f64::from_bits(ka[1])];
            let pb = [f64::from_bits(kb[0]), f64::from_bits(kb[1])];
            let len = ((pa[0] - pb[0]).powi(2) + (pa[1] - pb[1]).powi(2)).sqrt();
            let mid = [(pa[0] + pb[0]) * 0.5, (pa[1] + pb[1]) * 0.5];
            if rim_distance(mid) <= 0.3 {
                rim_edges.push(len);
            } else {
                far_edges.push(len);
            }
        }
        println!(
            "  graded locality: rim mean {:.4} vs far mean {:.4} (ratio {:.3})",
            mean(&rim_edges),
            mean(&far_edges),
            mean(&rim_edges) / mean(&far_edges)
        );
    }

    println!();
    println!(
        "note: debug-build timings are 10-50x the release profile; publish \
         with --release for comparable numbers"
    );
}
