use criterion::{criterion_group, criterion_main, Criterion};
use glam::DVec3;
use rcad_algorithms::{Bvh, HlrCamera, hlr, mesh_brep, TessellationParams};
use rcad_modeling::{make_box_brep, make_cylinder_brep, make_sphere_brep};

// ── HLR ───────────────────────────────────────────────────────────────────────

fn hlr_box_isometric(c: &mut Criterion) {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 5.0, 5.0, 5.0).unwrap();
    let camera = HlrCamera::isometric(20.0);
    c.bench_function("hlr_box_isometric", |bench| {
        bench.iter(|| hlr(&brep, &camera, 8))
    });
}

fn hlr_cylinder_silhouette(c: &mut Criterion) {
    let brep = make_cylinder_brep(
        DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 4.0
    ).unwrap();
    let camera = HlrCamera::right(10.0);
    c.bench_function("hlr_cylinder_silhouette", |bench| {
        bench.iter(|| hlr(&brep, &camera, 8))
    });
}

fn hlr_sphere(c: &mut Criterion) {
    let brep = make_sphere_brep(DVec3::ZERO, 2.0).unwrap();
    let camera = HlrCamera::front(15.0);
    c.bench_function("hlr_sphere", |bench| {
        bench.iter(|| hlr(&brep, &camera, 8))
    });
}

// ── mesh_brep ────────────────────────────────────────────────────────────────

fn mesh_brep_box(c: &mut Criterion) {
    let mut brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).unwrap();
    let params = TessellationParams::default();
    c.bench_function("mesh_brep_box", |bench| {
        bench.iter(|| mesh_brep(&mut brep, &params))
    });
}

fn mesh_brep_sphere_coarse(c: &mut Criterion) {
    let mut brep = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
    let params = TessellationParams {
        chord_tolerance: 0.1,
        ..TessellationParams::default()
    };
    c.bench_function("mesh_brep_sphere_coarse", |bench| {
        bench.iter(|| mesh_brep(&mut brep, &params))
    });
}

fn mesh_brep_sphere_fine(c: &mut Criterion) {
    let mut brep = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
    let params = TessellationParams {
        chord_tolerance: 0.005,
        ..TessellationParams::default()
    };
    c.bench_function("mesh_brep_sphere_fine", |bench| {
        bench.iter(|| mesh_brep(&mut brep, &params))
    });
}

fn mesh_brep_cylinder(c: &mut Criterion) {
    let mut brep = make_cylinder_brep(
        DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 6.0
    ).unwrap();
    let params = TessellationParams::default();
    c.bench_function("mesh_brep_cylinder", |bench| {
        bench.iter(|| mesh_brep(&mut brep, &params))
    });
}

// ── BVH ──────────────────────────────────────────────────────────────────────

fn bvh_build_box_6_faces(c: &mut Criterion) {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).unwrap();
    c.bench_function("bvh_build_box_6_faces", |bench| {
        bench.iter(|| Bvh::build(&brep))
    });
}

fn bvh_build_sphere(c: &mut Criterion) {
    let brep = make_sphere_brep(DVec3::ZERO, 5.0).unwrap();
    c.bench_function("bvh_build_sphere", |bench| {
        bench.iter(|| Bvh::build(&brep))
    });
}

fn bvh_query_sphere(c: &mut Criterion) {
    let brep = make_sphere_brep(DVec3::ZERO, 5.0).unwrap();
    let bvh = Bvh::build(&brep);
    let aabb = rcad_algorithms::Aabb {
        min: glam::DVec3::new(-1.0, -1.0, -1.0),
        max: glam::DVec3::new(1.0, 1.0, 1.0),
    };
    c.bench_function("bvh_query_sphere", |bench| {
        bench.iter(|| bvh.query_aabb(&aabb))
    });
}

criterion_group!(
    benches,
    hlr_box_isometric,
    hlr_cylinder_silhouette,
    hlr_sphere,
    mesh_brep_box,
    mesh_brep_sphere_coarse,
    mesh_brep_sphere_fine,
    mesh_brep_cylinder,
    bvh_build_box_6_faces,
    bvh_build_sphere,
    bvh_query_sphere,
);
criterion_main!(benches);
