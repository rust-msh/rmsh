//! Comprehensive quality regression tests for all 2-D and 3-D meshing algorithms.
//!
//! Each test runs a geometry through every available algorithm, captures quality
//! metrics, and asserts they stay above a baseline.  When an algorithm improves,
//! the baseline numbers in this file should be updated.
//!
//! Run with: `cargo test -p rmsh-algo --test mesher_quality_regression`

use rmsh_algo::{
    Bamg2D, Delaunay2D, Delaunay3D, Domain2D, Frontal3D, FrontalDelaunay2D,
    FrontalQuads2D, Hxt3D, MeshAdapt2D, MeshParams, Mesher2D, Mesher3D, MmgRemesh,
    QuadPaving2D,
};
use rmsh_model::{Element, ElementType, Mesh, Node};

// ─── Quality stats ────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy)]
struct TetQualityStats {
    min_dihedral_deg: f64,
    p95_radius_edge: f64,
    sliver_fraction: f64,
    tet_count: usize,
}

#[derive(Debug, Clone, Copy)]
struct TriQualityStats {
    min_angle_deg: f64,
    p95_aspect_ratio: f64,
    avg_scaled_jacobian: f64,
    tri_count: usize,
    quad_count: usize,
}

// ─── Geometry generators ─────────────────────────────────────────────────────

fn box_surface(lx: f64, ly: f64, lz: f64) -> Mesh {
    let mut mesh = Mesh::new();
    for (id, xyz) in [
        (1, [0.0, 0.0, 0.0]), (2, [lx, 0.0, 0.0]),
        (3, [lx, ly, 0.0]),   (4, [0.0, ly, 0.0]),
        (5, [0.0, 0.0, lz]),  (6, [lx, 0.0, lz]),
        (7, [lx, ly, lz]),    (8, [0.0, ly, lz]),
    ] { mesh.add_node(Node::new(id, xyz[0], xyz[1], xyz[2])); }
    for (id, nodes) in [
        (1, vec![1, 2, 3, 4]), (2, vec![5, 6, 7, 8]),
        (3, vec![1, 2, 6, 5]), (4, vec![2, 3, 7, 6]),
        (5, vec![3, 4, 8, 7]), (6, vec![4, 1, 5, 8]),
    ] { mesh.add_element(Element::new(id, ElementType::Quad4, nodes)); }
    mesh
}

fn cube_surface() -> Mesh { box_surface(1.0, 1.0, 1.0) }

fn rect_domain() -> Domain2D {
    Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 2.0], [0.0, 2.0]])
}

fn lshape_domain() -> Domain2D {
    Domain2D::from_outer(vec![
        [0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0], [1.0, 2.0], [0.0, 2.0],
    ])
}

fn hole_domain() -> Domain2D {
    Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0]])
        .with_hole(vec![[1.0, 1.0], [2.0, 1.0], [2.0, 2.0], [1.0, 2.0]])
}

// ─── 2-D quality helpers ──────────────────────────────────────────────────────

fn p2(mesh: &Mesh, id: u64) -> [f64; 2] {
    let p = &mesh.nodes[&id].position;
    [p.x, p.y]
}

fn tri_min_angle(mesh: &Mesh, e: &Element) -> f64 {
    let n = &e.node_ids;
    let a = p2(mesh, n[0]); let b = p2(mesh, n[1]); let c = p2(mesh, n[2]);
    let ab = [b[0]-a[0], b[1]-a[1]]; let ac = [c[0]-a[0], c[1]-a[1]];
    let bc = [c[0]-b[0], c[1]-b[1]]; let ba = [-ab[0], -ab[1]];
    let cb = [-bc[0], -bc[1]]; let ca = [-ac[0], -ac[1]];
    let ang = |u: [f64;2], v: [f64;2]| -> f64 {
        let d = (u[0]*u[0]+u[1]*u[1]).sqrt().max(1e-15) * (v[0]*v[0]+v[1]*v[1]).sqrt().max(1e-15);
        ((u[0]*v[0]+u[1]*v[1]) / d).clamp(-1.0, 1.0).acos().to_degrees()
    };
    ang(ab, ac).min(ang(ba, bc)).min(ang(ca, cb))
}

fn tri_aspect_ratio(mesh: &Mesh, e: &Element) -> f64 {
    let n = &e.node_ids;
    let a = p2(mesh, n[0]); let b = p2(mesh, n[1]); let c = p2(mesh, n[2]);
    let dl = |x: [f64;2], y: [f64;2]| ((y[0]-x[0]).powi(2) + (y[1]-x[1]).powi(2)).sqrt();
    let d01 = dl(a,b); let d12 = dl(b,c); let d20 = dl(c,a);
    let mx = d01.max(d12).max(d20); let mn = d01.min(d12).min(d20).max(1e-30);
    mx / mn
}

fn tri_scaled_jacobian(mesh: &Mesh, e: &Element) -> f64 {
    let n = &e.node_ids;
    let a = p2(mesh, n[0]); let b = p2(mesh, n[1]); let c = p2(mesh, n[2]);
    let e1 = [b[0]-a[0], b[1]-a[1]]; let e2 = [c[0]-a[0], c[1]-a[1]];
    let cross = (e1[0]*e2[1] - e1[1]*e2[0]).abs();
    let l1 = (e1[0]*e1[0]+e1[1]*e1[1]).sqrt().max(1e-15);
    let l2 = (e2[0]*e2[0]+e2[1]*e2[1]).sqrt().max(1e-15);
    cross / (l1 * l2)
}

fn tri_quality_stats(mesh: &Mesh) -> TriQualityStats {
    let mut min_angles = Vec::new();
    let mut ratios = Vec::new();
    let mut jacobians = Vec::new();
    let mut tri_count = 0;
    let mut quad_count = 0;

    for e in &mesh.elements {
        if e.etype == ElementType::Triangle3 && e.node_ids.len() == 3 {
            tri_count += 1;
            min_angles.push(tri_min_angle(mesh, e));
            ratios.push(tri_aspect_ratio(mesh, e));
            jacobians.push(tri_scaled_jacobian(mesh, e));
        } else if e.etype == ElementType::Quad4 && e.node_ids.len() == 4 {
            quad_count += 1;
        }
    }

    min_angles.sort_by(|a,b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    ratios.sort_by(|a,b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p95_ratio = if ratios.is_empty() { 0.0 } else {
        let idx = ((ratios.len() - 1) as f64 * 0.95).round() as usize;
        ratios[idx]
    };
    let avg_jac = if jacobians.is_empty() { 0.0 } else {
        jacobians.iter().sum::<f64>() / jacobians.len() as f64
    };

    TriQualityStats {
        min_angle_deg: *min_angles.first().unwrap_or(&0.0),
        p95_aspect_ratio: p95_ratio,
        avg_scaled_jacobian: avg_jac,
        tri_count,
        quad_count,
    }
}

// ─── 3-D quality helpers (from mesher3d_quality_regression.rs) ────────────────

fn p3(mesh: &Mesh, id: u64) -> [f64; 3] {
    let p = &mesh.nodes[&id].position;
    [p.x, p.y, p.z]
}
fn sub3(a: [f64;3], b: [f64;3]) -> [f64;3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }
fn dot3(a: [f64;3], b: [f64;3]) -> f64 { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
fn cross3(a: [f64;3], b: [f64;3]) -> [f64;3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn len3(v: [f64;3]) -> f64 { dot3(v, v).sqrt() }

fn min_dihedral_angle_tet(a: [f64;3], b: [f64;3], c: [f64;3], d: [f64;3]) -> f64 {
    let edges = [
        (a,b,c,d), (a,c,b,d), (a,d,b,c), (b,c,a,d), (b,d,a,c), (c,d,a,b),
    ];
    edges.iter().map(|&(p,q,r,s)| {
        let n1 = cross3(sub3(q,p), sub3(r,p));
        let n2 = cross3(sub3(q,p), sub3(s,p));
        let l1 = len3(n1).max(1e-15); let l2 = len3(n2).max(1e-15);
        (dot3(n1,n2)/(l1*l2)).clamp(-1.0,1.0).acos().to_degrees()
    }).fold(f64::MAX, f64::min)
}

fn tetra_volume(a: [f64;3], b: [f64;3], c: [f64;3], d: [f64;3]) -> f64 {
    dot3(sub3(a,d), cross3(sub3(b,d), sub3(c,d))).abs() / 6.0
}

fn solve_3x3(rows: [([f64;3], f64); 3]) -> Option<[f64;3]> {
    let mut a = [[0.0;4];3];
    for i in 0..3 {
        a[i][0]=rows[i].0[0]; a[i][1]=rows[i].0[1]; a[i][2]=rows[i].0[2]; a[i][3]=rows[i].1;
    } for col in 0..3 { let mut pv=col; for r in (col+1)..3 { if a[r][col].abs()>a[pv][col].abs(){pv=r} }
    if a[pv][col].abs()<1e-15{return None} a.swap(pv,col); let inv=1.0/a[col][col];
    for j in col..4{a[col][j]*=inv} for r in 0..3{if r==col{continue} let f=a[r][col];
    for j in col..4{a[r][j]-=f*a[col][j]}} } Some([a[0][3],a[1][3],a[2][3]])
}

fn radius_edge_ratio(a: [f64;3], b: [f64;3], c: [f64;3], d: [f64;3]) -> f64 {
    let edges = [len3(sub3(a,b)), len3(sub3(a,c)), len3(sub3(a,d)),
                 len3(sub3(b,c)), len3(sub3(b,d)), len3(sub3(c,d))];
    let min_edge = edges.iter().copied().fold(f64::MAX, f64::min).max(1e-15);
    let r = match solve_3x3([(sub3(b,a),0.5*(dot3(b,b)-dot3(a,a))),
                              (sub3(c,a),0.5*(dot3(c,c)-dot3(a,a))),
                              (sub3(d,a),0.5*(dot3(d,d)-dot3(a,a)))]) {
        Some(cc) => len3(sub3(cc,a)), None => f64::INFINITY,
    }; r / min_edge
}

fn tet_quality_stats(mesh: &Mesh) -> TetQualityStats {
    let mut min_dihedral = f64::INFINITY;
    let mut ratios = Vec::new();
    let mut slivers = 0usize;
    let mut total = 0usize;
    for e in &mesh.elements {
        if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 { continue; }
        let (a,b,c,d) = (p3(mesh,e.node_ids[0]),p3(mesh,e.node_ids[1]),
                          p3(mesh,e.node_ids[2]),p3(mesh,e.node_ids[3]));
        if tetra_volume(a,b,c,d) < 1e-15 { continue; }
        total += 1;
        let d_ang = min_dihedral_angle_tet(a,b,c,d);
        let ratio = radius_edge_ratio(a,b,c,d);
        min_dihedral = min_dihedral.min(d_ang);
        ratios.push(ratio);
        if d_ang < 6.0 && ratio > 1.8 { slivers += 1; }
    } ratios.sort_by(|x,y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let p95 = if ratios.is_empty() { f64::INFINITY } else {
        ratios[((ratios.len()-1) as f64 * 0.95).round() as usize]
    }; TetQualityStats {
        min_dihedral_deg: if min_dihedral.is_finite(){min_dihedral}else{0.0},
        p95_radius_edge: p95, sliver_fraction: if total==0{1.0}else{slivers as f64/total as f64},
        tet_count: total,
    }
}

// ─── Test runners ─────────────────────────────────────────────────────────────

fn assert_tri_baseline(name: &str, s: TriQualityStats) {
    eprintln!(
        "[2D] {name}: tris={}, quads={}, min_angle={:.2}°, p95_aspect={:.2}, avg_jac={:.4}",
        s.tri_count, s.quad_count, s.min_angle_deg, s.p95_aspect_ratio, s.avg_scaled_jacobian
    );
    assert!(s.tri_count + s.quad_count > 0, "{name}: no elements");
    // BAMG is anisotropic — can produce slivers. Check existence only.
    if name.contains("bamg") { return; }
    // For all-quad meshes, skip the tri-specific assertions.
    if s.tri_count == 0 { return; }
    assert!(s.min_angle_deg > 1.0, "{name}: degenerate elements (min_angle={:.2}°)", s.min_angle_deg);
    assert!(s.avg_scaled_jacobian > 0.1, "{name}: poor avg jacobian {:.4}", s.avg_scaled_jacobian);
}

fn assert_tet_baseline(name: &str, s: TetQualityStats) {
    eprintln!(
        "[3D] {name}: tets={}, min_dihedral={:.3}°, p95_rer={:.3}, sliver={:.3}",
        s.tet_count, s.min_dihedral_deg, s.p95_radius_edge, s.sliver_fraction
    );
    assert!(s.tet_count > 0, "{name}: no tets");
    assert!(s.min_dihedral_deg > 0.10, "{name}: min dihedral too low");
    assert!(s.p95_radius_edge < 120.0, "{name}: radius-edge exploded");
    assert!(s.sliver_fraction < 0.40, "{name}: too many slivers");
}

// ─── 2-D quality tests ───────────────────────────────────────────────────────

#[test]
fn mesher2d_quality_rectangle() {
    let domain = rect_domain();
    let p = MeshParams::with_size(0.35);
    let algorithms: [(&str, Box<dyn Fn() -> Box<dyn Mesher2D>>); 5] = [
        ("adapt2d",      Box::new(|| Box::new(MeshAdapt2D::default()))),
        ("delaunay2d",   Box::new(|| Box::new(Delaunay2D::default()))),
        ("frontal2d",    Box::new(|| Box::new(FrontalDelaunay2D::default()))),
        ("bamg2d",       Box::new(|| Box::new(Bamg2D::default()))),
        ("quadpaving",   Box::new(|| Box::new(QuadPaving2D::default()))),
    ];
    for (name, maker) in &algorithms {
        let mesh = maker().mesh_2d(&domain, &p).unwrap();
        assert_tri_baseline(name, tri_quality_stats(&mesh));
    }
}

#[test]
fn mesher2d_quality_lshape() {
    let domain = lshape_domain();
    let p = MeshParams::with_size(0.3);
    let algs: [(&str, Box<dyn Fn() -> Box<dyn Mesher2D>>); 4] = [
        ("adapt2d",      Box::new(|| Box::new(MeshAdapt2D::default()))),
        ("delaunay2d",   Box::new(|| Box::new(Delaunay2D::default()))),
        ("frontal2d",    Box::new(|| Box::new(FrontalDelaunay2D::default()))),
        ("bamg2d",       Box::new(|| Box::new(Bamg2D::default()))),
    ];
    for (name, maker) in &algs {
        let mesh = maker().mesh_2d(&domain, &p).unwrap();
        let s = tri_quality_stats(&mesh);
        assert!(s.tri_count > 0, "{name}: L-shape produced no elements");
        assert_tri_baseline(name, s);
    }
}

#[test]
fn mesher2d_quality_with_hole() {
    let domain = hole_domain();
    let p = MeshParams::with_size(0.5);
    let algs: [(&str, Box<dyn Fn() -> Box<dyn Mesher2D>>); 4] = [
        ("delaunay2d",   Box::new(|| Box::new(Delaunay2D::default()))),
        ("frontal2d",    Box::new(|| Box::new(FrontalDelaunay2D::default()))),
        ("bamg2d",       Box::new(|| Box::new(Bamg2D::default()))),
        ("adapt2d",      Box::new(|| Box::new(MeshAdapt2D::default()))),
    ];
    for (name, maker) in &algs {
        let mesh = maker().mesh_2d(&domain, &p).unwrap();
        let s = tri_quality_stats(&mesh);
        assert!(s.tri_count > 0, "{name}: hole domain produced no elements");
        assert_tri_baseline(name, s);
    }
}

#[test]
fn mesher2d_quad_recombination() {
    let domain = rect_domain();
    let p = MeshParams::with_size(0.5);
    let mesh = FrontalQuads2D::default().mesh_2d(&domain, &p).unwrap();
    let s = tri_quality_stats(&mesh);
    eprintln!(
        "[2D] frontal-quads: tris={}, quads={}, min_angle={:.2}°, avg_jac={:.4}",
        s.tri_count, s.quad_count, s.min_angle_deg, s.avg_scaled_jacobian
    );
    assert!(s.tri_count + s.quad_count > 0, "frontal-quads: no elements");
    assert!(s.quad_count > 0, "frontal-quads: no quads produced");
}

// ─── 3-D quality tests (existing, extended) ───────────────────────────────────

#[test]
fn mesher3d_quality_cube() {
    let surface = cube_surface();
    let params = MeshParams::with_size(0.4);
    let algs: [(&str, Box<dyn Fn() -> Box<dyn Mesher3D>>); 4] = [
        ("delaunay3d", Box::new(|| Box::new(Delaunay3D::default()))),
        ("frontal3d",  Box::new(|| Box::new(Frontal3D::default()))),
        ("hxt3d",      Box::new(|| Box::new(Hxt3D::default()))),
        ("mmg3d",      Box::new(|| Box::new(MmgRemesh::default()))),
    ];
    for (name, maker) in &algs {
        let mesh = maker().mesh_3d(&surface, &params).unwrap();
        assert_tet_baseline(name, tet_quality_stats(&mesh));
    }
}

#[test]
fn mesher3d_quality_stretched_box() {
    let surface = box_surface(1.8, 1.1, 0.9);
    let p = MeshParams::with_size(0.38);
    let algs: [(&str, Box<dyn Fn() -> Box<dyn Mesher3D>>); 3] = [
        ("delaunay3d", Box::new(|| Box::new(Delaunay3D::default()))),
        ("frontal3d",  Box::new(|| Box::new(Frontal3D::default()))),
        ("hxt3d",      Box::new(|| Box::new(Hxt3D::default()))),
    ];
    for (name, maker) in &algs {
        let mesh = maker().mesh_3d(&surface, &p).unwrap();
        assert_tet_baseline(name, tet_quality_stats(&mesh));
    }
}

#[test]
fn mesher3d_quality_slender_box() {
    let surface = box_surface(3.2, 0.55, 0.45);
    let p = MeshParams::with_size(0.26);
    let algs: [(&str, Box<dyn Fn() -> Box<dyn Mesher3D>>); 3] = [
        ("delaunay3d", Box::new(|| Box::new(Delaunay3D::default()))),
        ("frontal3d",  Box::new(|| Box::new(Frontal3D::default()))),
        ("hxt3d",      Box::new(|| Box::new(Hxt3D::default()))),
    ];
    for (name, maker) in &algs {
        let mesh = maker().mesh_3d(&surface, &p).unwrap();
        assert_tet_baseline(name, tet_quality_stats(&mesh));
    }
}
