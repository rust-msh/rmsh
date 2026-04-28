use rmsh_algo::{Delaunay3D, Frontal3D, Hxt3D, MeshParams, Mesher3D};
use rmsh_model::{Element, ElementType, Mesh, Node};

#[derive(Debug, Clone, Copy)]
struct TetQualityStats {
    min_dihedral_deg: f64,
    p95_radius_edge: f64,
    sliver_fraction: f64,
    tet_count: usize,
}

fn box_surface(lx: f64, ly: f64, lz: f64) -> Mesh {
    let mut mesh = Mesh::new();
    for (id, xyz) in [
        (1, [0.0, 0.0, 0.0]),
        (2, [lx, 0.0, 0.0]),
        (3, [lx, ly, 0.0]),
        (4, [0.0, ly, 0.0]),
        (5, [0.0, 0.0, lz]),
        (6, [lx, 0.0, lz]),
        (7, [lx, ly, lz]),
        (8, [0.0, ly, lz]),
    ] {
        mesh.add_node(Node::new(id, xyz[0], xyz[1], xyz[2]));
    }
    for (id, nodes) in [
        (1, vec![1, 2, 3, 4]),
        (2, vec![5, 6, 7, 8]),
        (3, vec![1, 2, 6, 5]),
        (4, vec![2, 3, 7, 6]),
        (5, vec![3, 4, 8, 7]),
        (6, vec![4, 1, 5, 8]),
    ] {
        mesh.add_element(Element::new(id, ElementType::Quad4, nodes));
    }
    mesh
}

fn cube_surface() -> Mesh {
    box_surface(1.0, 1.0, 1.0)
}

fn p3(mesh: &Mesh, id: u64) -> [f64; 3] {
    let p = mesh
        .nodes
        .get(&id)
        .expect("element references existing node")
        .position;
    [p.x, p.y, p.z]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

fn min_dihedral_angle_tet(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let edges = [
        (a, b, c, d),
        (a, c, b, d),
        (a, d, b, c),
        (b, c, a, d),
        (b, d, a, c),
        (c, d, a, b),
    ];
    edges
        .iter()
        .map(|&(p, q, r, s)| {
            let pq = sub3(q, p);
            let pr = sub3(r, p);
            let ps = sub3(s, p);
            let n1 = cross3(pq, pr);
            let n2 = cross3(pq, ps);
            let l1 = len3(n1);
            let l2 = len3(n2);
            if l1 < 1e-15 || l2 < 1e-15 {
                return 0.0;
            }
            (dot3(n1, n2) / (l1 * l2)).clamp(-1.0, 1.0).acos().to_degrees()
        })
        .fold(f64::MAX, f64::min)
}

fn tetra_volume(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let ad = sub3(a, d);
    let bd = sub3(b, d);
    let cd = sub3(c, d);
    dot3(ad, cross3(bd, cd)).abs() / 6.0
}

fn solve_3x3(rows: [([f64; 3], f64); 3]) -> Option<[f64; 3]> {
    let mut a = [[0.0; 4]; 3];
    for i in 0..3 {
        a[i][0] = rows[i].0[0];
        a[i][1] = rows[i].0[1];
        a[i][2] = rows[i].0[2];
        a[i][3] = rows[i].1;
    }

    for col in 0..3 {
        let mut pivot = col;
        for r in (col + 1)..3 {
            if a[r][col].abs() > a[pivot][col].abs() {
                pivot = r;
            }
        }
        if a[pivot][col].abs() < 1e-15 {
            return None;
        }
        if pivot != col {
            a.swap(pivot, col);
        }
        let inv = 1.0 / a[col][col];
        for j in col..4 {
            a[col][j] *= inv;
        }
        for r in 0..3 {
            if r == col {
                continue;
            }
            let f = a[r][col];
            for j in col..4 {
                a[r][j] -= f * a[col][j];
            }
        }
    }

    Some([a[0][3], a[1][3], a[2][3]])
}

fn tetra_circumcenter(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Option<[f64; 3]> {
    let ba = sub3(b, a);
    let ca = sub3(c, a);
    let da = sub3(d, a);
    let rhs0 = 0.5 * (dot3(b, b) - dot3(a, a));
    let rhs1 = 0.5 * (dot3(c, c) - dot3(a, a));
    let rhs2 = 0.5 * (dot3(d, d) - dot3(a, a));
    solve_3x3([(ba, rhs0), (ca, rhs1), (da, rhs2)])
}

fn radius_edge_ratio(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let edges = [
        len3(sub3(a, b)),
        len3(sub3(a, c)),
        len3(sub3(a, d)),
        len3(sub3(b, c)),
        len3(sub3(b, d)),
        len3(sub3(c, d)),
    ];
    let min_edge = edges.iter().copied().fold(f64::MAX, f64::min).max(1e-15);
    let r = match tetra_circumcenter(a, b, c, d) {
        Some(cc) => len3(sub3(cc, a)),
        None => f64::INFINITY,
    };
    r / min_edge
}

fn tet_quality_stats(mesh: &Mesh) -> TetQualityStats {
    let mut min_dihedral = f64::INFINITY;
    let mut ratios = Vec::<f64>::new();
    let mut slivers = 0usize;
    let mut total = 0usize;

    for e in &mesh.elements {
        if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
            continue;
        }
        let a = p3(mesh, e.node_ids[0]);
        let b = p3(mesh, e.node_ids[1]);
        let c = p3(mesh, e.node_ids[2]);
        let d = p3(mesh, e.node_ids[3]);

        let volume = tetra_volume(a, b, c, d);
        if volume < 1e-15 {
            continue;
        }

        total += 1;
        let dihedral = min_dihedral_angle_tet(a, b, c, d);
        let ratio = radius_edge_ratio(a, b, c, d);
        min_dihedral = min_dihedral.min(dihedral);
        ratios.push(ratio);
        if dihedral < 6.0 && ratio > 1.8 {
            slivers += 1;
        }
    }

    ratios.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let p95 = if ratios.is_empty() {
        f64::INFINITY
    } else {
        let idx = ((ratios.len() - 1) as f64 * 0.95).round() as usize;
        ratios[idx]
    };

    TetQualityStats {
        min_dihedral_deg: if min_dihedral.is_finite() { min_dihedral } else { 0.0 },
        p95_radius_edge: p95,
        sliver_fraction: if total == 0 {
            1.0
        } else {
            slivers as f64 / total as f64
        },
        tet_count: total,
    }
}

fn assert_quality_baseline(name: &str, stats: TetQualityStats) {
    eprintln!(
        "{}: tets={}, min_dihedral={:.3} deg, p95_radius_edge={:.3}, sliver_frac={:.3}",
        name,
        stats.tet_count,
        stats.min_dihedral_deg,
        stats.p95_radius_edge,
        stats.sliver_fraction
    );

    assert!(stats.tet_count > 0, "{} produced no tets", name);
    assert!(stats.min_dihedral_deg > 0.10, "{} min dihedral too low", name);
    assert!(stats.p95_radius_edge < 120.0, "{} radius-edge exploded", name);
    assert!(stats.sliver_fraction < 0.40, "{} too many slivers", name);
}

#[test]
fn mesher3d_quality_baseline_cube() {
    let surface = cube_surface();
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
    let surface = box_surface(1.8, 1.1, 0.9);
    let params = MeshParams::with_size(0.38);

    let d = Delaunay3D::default().mesh_3d(&surface, &params).unwrap();
    let f = Frontal3D::default().mesh_3d(&surface, &params).unwrap();
    let h = Hxt3D::default().mesh_3d(&surface, &params).unwrap();

    let qd = tet_quality_stats(&d);
    let qf = tet_quality_stats(&f);
    let qh = tet_quality_stats(&h);

    eprintln!(
        "delaunay3d_stretched: tets={}, min_dihedral={:.3} deg, p95_radius_edge={:.3}, sliver_frac={:.3}",
        qd.tet_count, qd.min_dihedral_deg, qd.p95_radius_edge, qd.sliver_fraction
    );
    eprintln!(
        "frontal3d_stretched: tets={}, min_dihedral={:.3} deg, p95_radius_edge={:.3}, sliver_frac={:.3}",
        qf.tet_count, qf.min_dihedral_deg, qf.p95_radius_edge, qf.sliver_fraction
    );
    eprintln!(
        "hxt3d_stretched: tets={}, min_dihedral={:.3} deg, p95_radius_edge={:.3}, sliver_frac={:.3}",
        qh.tet_count, qh.min_dihedral_deg, qh.p95_radius_edge, qh.sliver_fraction
    );

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
    let surface = box_surface(3.2, 0.55, 0.45);
    let params = MeshParams::with_size(0.26);

    let d = Delaunay3D::default().mesh_3d(&surface, &params).unwrap();
    let f = Frontal3D::default().mesh_3d(&surface, &params).unwrap();
    let h = Hxt3D::default().mesh_3d(&surface, &params).unwrap();

    let qd = tet_quality_stats(&d);
    let qf = tet_quality_stats(&f);
    let qh = tet_quality_stats(&h);

    eprintln!(
        "delaunay3d_slender: tets={}, min_dihedral={:.3} deg, p95_radius_edge={:.3}, sliver_frac={:.3}",
        qd.tet_count, qd.min_dihedral_deg, qd.p95_radius_edge, qd.sliver_fraction
    );
    eprintln!(
        "frontal3d_slender: tets={}, min_dihedral={:.3} deg, p95_radius_edge={:.3}, sliver_frac={:.3}",
        qf.tet_count, qf.min_dihedral_deg, qf.p95_radius_edge, qf.sliver_fraction
    );
    eprintln!(
        "hxt3d_slender: tets={}, min_dihedral={:.3} deg, p95_radius_edge={:.3}, sliver_frac={:.3}",
        qh.tet_count, qh.min_dihedral_deg, qh.p95_radius_edge, qh.sliver_fraction
    );

    assert!(qd.tet_count > 0 && qf.tet_count > 0 && qh.tet_count > 0);
    assert!(qd.min_dihedral_deg > 0.05 && qf.min_dihedral_deg > 0.05 && qh.min_dihedral_deg > 0.05);
    assert!(qd.p95_radius_edge.is_finite() && qf.p95_radius_edge.is_finite() && qh.p95_radius_edge.is_finite());
    assert!(qd.p95_radius_edge < 80000.0 && qf.p95_radius_edge < 80000.0 && qh.p95_radius_edge < 80000.0);
    assert!(qd.sliver_fraction < 0.60 && qf.sliver_fraction < 0.60 && qh.sliver_fraction < 0.60);
}
