//! Benchmark: point location in tetrahedral meshes.
//!
//! Compares the old linear-scan approach (`find_containing_tet_readonly`)
//! against the new uniform-grid spatial index (`SpatialGrid`).
//!
//! Run:  `cargo bench -p rmsh-algo -- point_location`

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};

use rmsh_algo::{
    Hxt3D, MeshParams, Mesher3D,
    SpatialGrid, cube_surface_mesh, find_containing_tet_readonly,
};
use rmsh_model::{ElementType, Mesh};

/// Generate a tet mesh and return (mesh, query_points).
fn setup_mesh(element_size: f64, num_queries: usize) -> (Mesh, Vec<[f64; 3]>) {
    let surface = cube_surface_mesh();
    let mesh = Hxt3D::default()
        .mesh_3d(&surface, &MeshParams::with_size(element_size))
        .expect("mesh generation");

    let tets: Vec<&rmsh_model::Element> = mesh
        .elements.iter()
        .filter(|e| e.etype == ElementType::Tetrahedron4 && e.node_ids.len() == 4)
        .collect();

    let mut queries = Vec::with_capacity(num_queries);
    for i in 0..num_queries.min(tets.len()) {
        let elt = tets[i % tets.len()];
        let mut sum = [0.0_f64; 3];
        for &nid in &elt.node_ids {
            if let Some(n) = mesh.nodes.get(&nid) {
                sum[0] += n.position.x; sum[1] += n.position.y; sum[2] += n.position.z;
            }
        }
        let n = elt.node_ids.len() as f64;
        queries.push([sum[0] / n, sum[1] / n, sum[2] / n]);
    }

    let vol = mesh.elements_by_dimension(3).len();
    eprintln!("  h={element_size:.2}: tets={vol}");
    (mesh, queries)
}

fn bench_point_location(c: &mut Criterion) {
    let sizes = [0.2, 0.15, 0.1];
    let mut group = c.benchmark_group("point_location");
    group.sample_size(30);

    for &elem_size in &sizes {
        let (mesh, queries) = setup_mesh(elem_size, 500);
        let tet_count = mesh.elements_by_dimension(3).len();
        if tet_count == 0 { continue; }
        let grid_res = (tet_count as f64).cbrt().ceil() as usize;
        let grid = SpatialGrid::build(&mesh, grid_res.max(4).min(32));

        group.bench_with_input(
            BenchmarkId::new(format!("linear_scan_h{elem_size}"), tet_count),
            &(mesh.clone(), queries.clone()),
            |b, (m, qs)| b.iter(|| {
                let mut f = 0usize;
                for p in qs { if find_containing_tet_readonly(m, *p).is_some() { f += 1; } }
                f
            }),
        );

        group.bench_with_input(
            BenchmarkId::new(format!("spatial_grid_h{elem_size}"), tet_count),
            &(mesh, queries, grid),
            |b, (m, qs, g)| b.iter(|| {
                let mut f = 0usize;
                for p in qs { if g.find_containing_tet(m, *p).is_some() { f += 1; } }
                f
            }),
        );
    }
    group.finish();
}

criterion_group!(benches, bench_point_location);
criterion_main!(benches);
