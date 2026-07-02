//! Comprehensive meshing algorithm benchmarks.
//!
//! Measures throughput (elements/sec) for every 2-D and 3-D meshing algorithm
//! across multiple mesh sizes.
//!
//! Run:  `cargo bench -p rmsh-algo -- meshing`

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

use rmsh_algo::{
    Bamg2D, Delaunay2D, Delaunay3D, Domain2D, Frontal3D, FrontalDelaunay2D,
    Hxt3D, MeshAdapt2D, MeshOptimizer, MeshParams, Mesher2D, Mesher3D, MmgRemesh,
    QuadPaving2D, promote_to_p2,
};
use rmsh_model::{Element, ElementType, Mesh, Node};

// ─── Geometry generators ──────────────────────────────────────────────────────

fn cube_surface_mesh() -> Mesh {
    let mut mesh = Mesh::new();
    for (id, xyz) in [
        (1, [0.0, 0.0, 0.0]), (2, [1.0, 0.0, 0.0]),
        (3, [1.0, 1.0, 0.0]), (4, [0.0, 1.0, 0.0]),
        (5, [0.0, 0.0, 1.0]), (6, [1.0, 0.0, 1.0]),
        (7, [1.0, 1.0, 1.0]), (8, [0.0, 1.0, 1.0]),
    ] { mesh.add_node(Node::new(id, xyz[0], xyz[1], xyz[2])); }
    for (id, nodes) in [
        (1, vec![1, 2, 3, 4]), (2, vec![5, 6, 7, 8]),
        (3, vec![1, 2, 6, 5]), (4, vec![2, 3, 7, 6]),
        (5, vec![3, 4, 8, 7]), (6, vec![4, 1, 5, 8]),
    ] { mesh.add_element(Element::new(id, ElementType::Quad4, nodes)); }
    mesh
}

fn rect_domain() -> Domain2D {
    Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 2.0], [0.0, 2.0]])
}

fn large_rect_domain() -> Domain2D {
    Domain2D::from_outer(vec![[0.0, 0.0], [10.0, 0.0], [10.0, 5.0], [0.0, 5.0]])
}

// ─── 2-D meshing benchmarks ───────────────────────────────────────────────────

fn bench_2d_meshing(c: &mut Criterion) {
    let mut group = c.benchmark_group("meshing_2d");
    group.sample_size(20);

    let domains: [(&str, &Domain2D, f64); 3] = [
        ("rect_coarse", &rect_domain(), 0.4),
        ("rect_medium", &rect_domain(), 0.2),
        ("large",       &large_rect_domain(), 0.3),
    ];

    let algorithms: [(&str, &dyn Fn() -> Box<dyn Mesher2D>); 5] = [
        ("adapt2d",      &(|| Box::new(MeshAdapt2D::default()) as Box<dyn Mesher2D>)),
        ("delaunay2d",   &(|| Box::new(Delaunay2D::default()) as Box<dyn Mesher2D>)),
        ("frontal2d",    &(|| Box::new(FrontalDelaunay2D::default()) as Box<dyn Mesher2D>)),
        ("bamg2d",       &(|| Box::new(Bamg2D::default()) as Box<dyn Mesher2D>)),
        ("quadpaving",   &(|| Box::new(QuadPaving2D::default()) as Box<dyn Mesher2D>)),
    ];

    for (algo_name, maker) in &algorithms {
        for (dom_name, domain, mesh_size) in &domains {
            let params = MeshParams::with_size(*mesh_size);
            let label = format!("{algo_name}_{dom_name}");

            group.bench_with_input(
                BenchmarkId::new(&label, format!("h{mesh_size}")),
                &(maker, domain, &params),
                |b, (mk, dm, p)| b.iter(|| {
                    mk().mesh_2d(dm, p).unwrap()
                }),
            );
        }
    }
    group.finish();
}

// ─── 3-D meshing benchmarks ───────────────────────────────────────────────────

fn bench_3d_meshing(c: &mut Criterion) {
    let mut group = c.benchmark_group("meshing_3d");
    group.sample_size(15);

    let surface = cube_surface_mesh();
    let sizes = [0.4, 0.25, 0.15];

    let algorithms: [(&str, &dyn Fn(f64) -> Box<dyn Mesher3D>); 4] = [
        ("delaunay3d", &(|_: f64| Box::new(Delaunay3D::default()) as Box<dyn Mesher3D>)),
        ("frontal3d",  &(|_: f64| Box::new(Frontal3D::default()) as Box<dyn Mesher3D>)),
        ("hxt3d",      &(|_: f64| Box::new(Hxt3D::default()) as Box<dyn Mesher3D>)),
        ("mmg3d",      &(|_: f64| Box::new(MmgRemesh::default()) as Box<dyn Mesher3D>)),
    ];

    for &elem_size in &sizes {
        let params = MeshParams::with_size(elem_size);

        for (algo_name, make_algo) in &algorithms {
            let label = format!("{algo_name}_h{elem_size}");

            group.bench_with_input(
                BenchmarkId::new(&label, format!("{elem_size}")),
                &(&surface, &params, make_algo),
                |b, (surf, p, mk)| b.iter(|| {
                    mk(elem_size).mesh_3d(surf, p).unwrap()
                }),
            );
        }
    }
    group.finish();
}

// ─── Mesh operation benchmarks ────────────────────────────────────────────────

fn bench_mesh_ops(c: &mut Criterion) {
    let mut group = c.benchmark_group("mesh_ops");
    group.sample_size(20);

    // Generate a medium-resolution cube mesh for operations.
    let surface = cube_surface_mesh();
    let params = MeshParams::with_size(0.2);
    let tet_mesh = Hxt3D::default().mesh_3d(&surface, &params).unwrap();
    let tet_count = tet_mesh.elements_by_dimension(3).len();
    eprintln!("  mesh_ops mesh: {tet_count} tets");

    group.bench_with_input(
        BenchmarkId::new("laplacian_smooth", tet_count),
        &tet_mesh,
        |b, m| b.iter(|| {
            use rmsh_algo::LaplacianSmooth;
            let mut mesh = m.clone();
            LaplacianSmooth::default()
                .optimize(&mut mesh, &rmsh_algo::OptimizeParams {
                    iterations: 5, ..Default::default()
                }).unwrap();
        }),
    );

    group.bench_with_input(
        BenchmarkId::new("p2_promote", tet_count),
        &tet_mesh,
        |b, m| b.iter(|| {
            let mut mesh = m.clone();
            promote_to_p2(&mut mesh);
        }),
    );

    group.finish();
}

criterion_group!(benches, bench_2d_meshing, bench_3d_meshing, bench_mesh_ops);
criterion_main!(benches);
