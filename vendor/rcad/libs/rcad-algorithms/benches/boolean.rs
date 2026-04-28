use criterion::{criterion_group, criterion_main, Criterion};
use glam::{DVec2, DVec3};
use rcad_algorithms::{boolean_op, boolean_op_par, BooleanOpType};
use rcad_modeling::{
    make_box_brep, make_cylinder_brep, make_sphere_brep,
    fillet_edge, loft, sweep_pipe,
};

// ── Boolean operations ────────────────────────────────────────────────────────

fn boolean_union_boxes(c: &mut Criterion) {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let b = make_box_brep(DVec3::new(0.5, 0.5, 0.5), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    c.bench_function("boolean_union_boxes", |bench| {
        bench.iter(|| boolean_op(BooleanOpType::Union, &a, &b).unwrap())
    });
}

fn boolean_diff_box_sphere(c: &mut Criterion) {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 1.0, 1.0), 0.8).unwrap();
    c.bench_function("boolean_diff_box_sphere", |bench| {
        bench.iter(|| boolean_op(BooleanOpType::Difference, &a, &b).unwrap())
    });
}

fn boolean_union_box_cylinder(c: &mut Criterion) {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let b = make_cylinder_brep(DVec3::new(1.0, 1.0, -0.5), DVec3::Z, DVec3::X, 0.5, 3.0).unwrap();
    c.bench_function("boolean_union_box_cylinder", |bench| {
        bench.iter(|| boolean_op(BooleanOpType::Union, &a, &b).unwrap())
    });
}

fn boolean_diff_box_cylinder(c: &mut Criterion) {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let b = make_cylinder_brep(DVec3::new(1.0, 1.0, -0.5), DVec3::Z, DVec3::X, 0.4, 3.0).unwrap();
    c.bench_function("boolean_diff_box_cylinder", |bench| {
        bench.iter(|| boolean_op(BooleanOpType::Difference, &a, &b).unwrap())
    });
}

// ── Parallel Boolean operations ───────────────────────────────────────────────

fn boolean_union_boxes_par(c: &mut Criterion) {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let b = make_box_brep(DVec3::new(0.5, 0.5, 0.5), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    c.bench_function("boolean_union_boxes_par", |bench| {
        bench.iter(|| boolean_op_par(BooleanOpType::Union, &a, &b).unwrap())
    });
}

fn boolean_diff_box_sphere_par(c: &mut Criterion) {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 1.0, 1.0), 0.8).unwrap();
    c.bench_function("boolean_diff_box_sphere_par", |bench| {
        bench.iter(|| boolean_op_par(BooleanOpType::Difference, &a, &b).unwrap())
    });
}

fn boolean_union_box_cylinder_par(c: &mut Criterion) {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let b = make_cylinder_brep(DVec3::new(1.0, 1.0, -0.5), DVec3::Z, DVec3::X, 0.5, 3.0).unwrap();
    c.bench_function("boolean_union_box_cylinder_par", |bench| {
        bench.iter(|| boolean_op_par(BooleanOpType::Union, &a, &b).unwrap())
    });
}

fn boolean_diff_box_cylinder_par(c: &mut Criterion) {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let b = make_cylinder_brep(DVec3::new(1.0, 1.0, -0.5), DVec3::Z, DVec3::X, 0.4, 3.0).unwrap();
    c.bench_function("boolean_diff_box_cylinder_par", |bench| {
        bench.iter(|| boolean_op_par(BooleanOpType::Difference, &a, &b).unwrap())
    });
}

// ── Surface intersection ──────────────────────────────────────────────────────

fn intss_plane_sphere(c: &mut Criterion) {
    use rcad_algorithms::inttools::intss::intersect_surfaces;
    use rcad_kernel::geom::{Plane, SphericalSurface, Surface3};
    let s1 = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
    let s2 = Surface3::Sphere(SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 1.0 });
    c.bench_function("intss_plane_sphere", |bench| {
        bench.iter(|| intersect_surfaces(&s1, &s2))
    });
}

// ── Fillet ────────────────────────────────────────────────────────────────────

fn fillet_box_edge(c: &mut Criterion) {
    // Fillet one edge of a box with radius 0.1
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    c.bench_function("fillet_box_edge", |bench| {
        bench.iter(|| fillet_edge(&brep, 0, 0.1).unwrap())
    });
}

// ── Loft ──────────────────────────────────────────────────────────────────────

fn loft_two_circles(c: &mut Criterion) {
    // Loft between two circular cross-sections (approximated as polygons)
    let n = 32usize;
    let profile1: Vec<DVec3> = (0..n)
        .map(|i| {
            let angle = std::f64::consts::TAU * i as f64 / n as f64;
            DVec3::new(angle.cos(), angle.sin(), 0.0)
        })
        .collect();
    let profile2: Vec<DVec3> = (0..n)
        .map(|i| {
            let angle = std::f64::consts::TAU * i as f64 / n as f64;
            DVec3::new(0.5 * angle.cos(), 0.5 * angle.sin(), 2.0)
        })
        .collect();
    c.bench_function("loft_two_circles", |bench| {
        bench.iter(|| loft(&[profile1.clone(), profile2.clone()]).unwrap())
    });
}

// ── Sweep ─────────────────────────────────────────────────────────────────────

fn sweep_square_along_arc(c: &mut Criterion) {
    // Sweep a square profile along a curved spine
    let profile: Vec<DVec2> = vec![
        DVec2::new(-0.1, -0.1),
        DVec2::new(0.1, -0.1),
        DVec2::new(0.1, 0.1),
        DVec2::new(-0.1, 0.1),
    ];
    let spine: Vec<DVec3> = (0..=32)
        .map(|i| {
            let t = i as f64 / 32.0 * std::f64::consts::PI;
            DVec3::new(t.cos(), t.sin(), 0.1 * i as f64)
        })
        .collect();
    c.bench_function("sweep_square_along_arc", |bench| {
        bench.iter(|| sweep_pipe(&profile, &spine).unwrap())
    });
}

criterion_group!(
    benches,
    boolean_union_boxes,
    boolean_diff_box_sphere,
    boolean_union_box_cylinder,
    boolean_diff_box_cylinder,
    boolean_union_boxes_par,
    boolean_diff_box_sphere_par,
    boolean_union_box_cylinder_par,
    boolean_diff_box_cylinder_par,
    intss_plane_sphere,
    fillet_box_edge,
    loft_two_circles,
    sweep_square_along_arc,
);
criterion_main!(benches);
