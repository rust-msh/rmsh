use criterion::{criterion_group, criterion_main, Criterion};
use glam::DVec3;
use rcad_kernel::properties::volume;
use rcad_kernel::{BRep, PrimitiveSolid};

fn volume_sphere(c: &mut Criterion) {
    let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
    c.bench_function("volume_sphere", |bench| bench.iter(|| volume(&brep)));
}

fn closest_point_sphere(c: &mut Criterion) {
    use rcad_kernel::geom::{SphericalSurface, Surface3};
    use rcad_kernel::projection::closest_point_on_surface;
    let surface = Surface3::Sphere(SphericalSurface {
        center: DVec3::ZERO,
        axis: DVec3::Z,
        radius: 1.0,
    });
    let points: Vec<DVec3> = (0..100)
        .map(|i| {
            let t = i as f64 * 0.063;
            DVec3::new(t.cos() * 2.0, t.sin() * 2.0, t.cos() * t.sin())
        })
        .collect();
    c.bench_function("closest_point_sphere_100", |bench| {
        bench.iter(|| {
            for p in &points {
                closest_point_on_surface(&surface, *p, 16);
            }
        })
    });
}

criterion_group!(benches, volume_sphere, closest_point_sphere);
criterion_main!(benches);
