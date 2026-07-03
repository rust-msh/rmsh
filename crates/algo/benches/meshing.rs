//! Comprehensive meshing algorithm benchmarks with multiple mesh sizes.
//!
//! Run:  `cargo bench -p rmsh-algo -- meshing`

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

use rmsh_algo::{
    Bamg2D, Delaunay2D, Delaunay3D, Domain2D, Frontal3D, FrontalDelaunay2D,
    Hxt3D, LaplacianSmooth, MeshAdapt2D, MeshOptimizer, MeshParams, Mesher2D, Mesher3D, MmgRemesh,
    QuadPaving2D, promote_to_p2,
};
use rmsh_model::{Element, ElementType, Mesh, Node};

// ─── Geometry generators ──────────────────────────────────────────────────

fn box_surface(lx: f64, ly: f64, lz: f64) -> Mesh {
    let mut mesh = Mesh::new();
    for (id, xyz) in [
        (1, [0.0, 0.0, 0.0]), (2, [lx, 0.0, 0.0]),
        (3, [lx, ly, 0.0]),   (4, [0.0, ly, 0.0]),
        (5, [0.0, 0.0, lz]),  (6, [lx, 0.0, lz]),
        (7, [lx, ly, lz]),    (8, [0.0, ly, lz]),
    ] { mesh.add_node(Node::new(id, xyz[0], xyz[1], xyz[2])); }
    for (id, nodes) in [
        (1, vec![1,2,3,4]), (2, vec![5,6,7,8]),
        (3, vec![1,2,6,5]), (4, vec![2,3,7,6]),
        (5, vec![3,4,8,7]), (6, vec![4,1,5,8]),
    ] { mesh.add_element(Element::new(id, ElementType::Quad4, nodes)); }
    mesh
}

fn cube_surface_mesh() -> Mesh { box_surface(1.0, 1.0, 1.0) }
fn big_box_surface() -> Mesh { box_surface(3.0, 2.0, 1.5) }
fn long_box_surface() -> Mesh { box_surface(5.0, 0.6, 0.6) }

fn rect_domain() -> Domain2D {
    Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 2.0], [0.0, 2.0]])
}
fn large_rect_domain() -> Domain2D {
    Domain2D::from_outer(vec![[0.0, 0.0], [10.0, 0.0], [10.0, 5.0], [0.0, 5.0]])
}
fn xlarge_rect_domain() -> Domain2D {
    Domain2D::from_outer(vec![[0.0, 0.0], [20.0, 0.0], [20.0, 10.0], [0.0, 10.0]])
}
fn huge_rect_domain() -> Domain2D {
    Domain2D::from_outer(vec![[0.0, 0.0], [50.0, 0.0], [50.0, 30.0], [0.0, 30.0]])
}

// ─── 2-D meshing benchmarks ───────────────────────────────────────────────

fn bench_2d_meshing(c: &mut Criterion) {
    let mut group = c.benchmark_group("meshing_2d");
    group.sample_size(15);
    let domains: [(&str, &Domain2D, f64); 5] = [
        ("rect_coarse", &rect_domain(), 0.4),
        ("rect_medium", &rect_domain(), 0.2),
        ("large",  &large_rect_domain(), 0.3),
        ("xlarge", &xlarge_rect_domain(), 0.3),
        ("huge",   &huge_rect_domain(), 0.4),
    ];
    let algorithms: [(&str, &dyn Fn() -> Box<dyn Mesher2D>); 5] = [
        ("adapt2d",    &(|| Box::new(MeshAdapt2D::default()))),
        ("delaunay2d", &(|| Box::new(Delaunay2D::default()))),
        ("frontal2d",  &(|| Box::new(FrontalDelaunay2D::default()))),
        ("bamg2d",     &(|| Box::new(Bamg2D::default()))),
        ("quadpaving", &(|| Box::new(QuadPaving2D::default()))),
    ];
    for (algo_name, maker) in &algorithms {
        for (dom_name, domain, mesh_size) in &domains {
            let p = MeshParams::with_size(*mesh_size);
            group.bench_with_input(BenchmarkId::new(format!("{algo_name}_{dom_name}"), ""),
                &(maker, domain, &p), |b, (mk, dm, pr)| b.iter(|| mk().mesh_2d(dm, pr).unwrap()));
        }
    }
    group.finish();
}

// ─── 3-D meshing benchmarks ───────────────────────────────────────────────

fn bench_3d_meshing(c: &mut Criterion) {
    let mut group = c.benchmark_group("meshing_3d");
    group.sample_size(10);
    let geometries: [(&str, &Mesh, f64); 5] = [
        ("cube_coarse", &cube_surface_mesh(), 0.4),
        ("cube_medium", &cube_surface_mesh(), 0.25),
        ("cube_fine",   &cube_surface_mesh(), 0.15),
        ("bigbox",      &big_box_surface(), 0.4),
        ("longbox",     &long_box_surface(), 0.3),
    ];
    let algorithms: [(&str, &dyn Fn() -> Box<dyn Mesher3D>); 4] = [
        ("delaunay3d", &(|| Box::new(Delaunay3D::default()))),
        ("frontal3d",  &(|| Box::new(Frontal3D::default()))),
        ("hxt3d",      &(|| Box::new(Hxt3D::default()))),
        ("mmg3d",      &(|| Box::new(MmgRemesh::default()))),
    ];
    for (geom_name, surface, elem_size) in &geometries {
        let p = MeshParams::with_size(*elem_size);
        for (algo_name, mk) in &algorithms {
            group.bench_with_input(BenchmarkId::new(format!("{algo_name}_{geom_name}"), ""),
                &(surface, &p, mk), |b, (s, pr, mk2)| b.iter(|| mk2().mesh_3d(s, pr).unwrap()));
        }
    }
    group.finish();
}

// ─── Mesh operation benchmarks ────────────────────────────────────────────

fn bench_mesh_ops(c: &mut Criterion) {
    let mut group = c.benchmark_group("mesh_ops");
    group.sample_size(15);
    let surface = cube_surface_mesh();
    let params = MeshParams::with_size(0.2);
    let tet_mesh = Hxt3D::default().mesh_3d(&surface, &params).unwrap();
    let tet_count = tet_mesh.elements_by_dimension(3).len();
    eprintln!("  mesh_ops mesh: {tet_count} tets");

    group.bench_with_input(BenchmarkId::new("laplacian_smooth", tet_count), &tet_mesh,
        |b, m| b.iter(|| {
            let mut mesh = m.clone();
            LaplacianSmooth::default().optimize(&mut mesh,
                &rmsh_algo::OptimizeParams { iterations: 5, ..Default::default() }).unwrap();
        }));

    group.bench_with_input(BenchmarkId::new("p2_promote", tet_count), &tet_mesh,
        |b, m| b.iter(|| { let mut mesh = m.clone(); promote_to_p2(&mut mesh); }));

    group.finish();
}

criterion_group!(benches, bench_2d_meshing, bench_3d_meshing, bench_mesh_ops);
criterion_main!(benches);
