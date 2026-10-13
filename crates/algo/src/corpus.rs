//! The golden corpus (P6): the named scenario families, the quality
//! statistics, and the gates — one shared source for the corpus pins,
//! the public benchmark run, and future consumers (CLI bench verb, the
//! Python binding).
//!
//! # What a corpus entry is
//!
//! A scenario family = a geometry builder (domain / surface) + a size
//! prescription + deterministic statistics (element counts, quality
//! measures) + a gate. The meshers are deterministic, so counts and
//! quality are exactly reproducible on any machine; only TIME varies.
//! The corpus pins live in `tests/corpus.rs` (always on), the timing
//! harness in `tests/corpus_benchmark.rs` (`#[ignore]`, run explicitly).
//!
//! # Gates and honesty
//!
//! The ROADMAP's post-optimization hard gate — min dihedral ≥ 5°,
//! p95 radius-edge ≤ 2.0, sliver fraction ≤ 0.1% — is the long-term
//! target ([`QualityGate::roadmap_hard_gate`]). Where the current
//! meshers do not meet it yet (the raw 3-D output predates the quality
//! card's spindle rework), the corpus pins the MEASURED baseline as the
//! regression bound and REPORTS the distance to the hard gate — never a
//! fake green. The 2-D families already clear a strict quality band
//! (the P3a Gmsh x1.15 relative lines are the acceptance oracle for the
//! plate family; those Gmsh-CLI comparisons stay in the mesher's own
//! test module).
//!
//! # Families (v1)
//!
//! | name | geometry | prescription | mesher |
//! |---|---|---|---|
//! | `plate_with_hole_uniform` | ±5 square, r=1 5-gon hole | uniform 0.5 | Frontal-Delaunay-2D |
//! | `plate_with_hole_graded` | same | curvature rim field [0.05, 2.0] | Frontal-Delaunay-2D |
//! | `l_shape_concave` | L polygon, area 3.0 | uniform 0.4 | mesh_polygon (CDT) |
//! | `cube_tet` | unit quad-faced cube | uniform 0.4 | Delaunay3D / Frontal3D / Hxt3D |
//! | `nested_cubes_multi_material` | unit cube + [0.1, 0.85]³ cavity | — | tetrahedralize_multi_material |
//! | `sharp_front_adaptive` | ±5 square, r=1 hole | tanh front, loop budget | AdaptiveLoop2D |

use glam::DVec3;

use rmsh_model::{Element, ElementType, Mesh, Node};

use crate::geometry_source::{GeometryHit, GeometrySource};
use crate::traits::Domain2D;

// ─── geometry helpers ─────────────────────────────────────────────────────────

/// CCW square ring with `segs_per_side` segments per side.
pub fn square_ring(half: f64, segs_per_side: usize) -> Vec<[f64; 2]> {
    let mut ring = Vec::with_capacity(segs_per_side * 4);
    let side = 2.0 * half;
    for k in 0..segs_per_side {
        let t = side * k as f64 / segs_per_side as f64;
        ring.push([-half + t, -half]);
    }
    for k in 0..segs_per_side {
        let t = side * k as f64 / segs_per_side as f64;
        ring.push([half, -half + t]);
    }
    for k in 0..segs_per_side {
        let t = side * k as f64 / segs_per_side as f64;
        ring.push([half - t, half]);
    }
    for k in 0..segs_per_side {
        let t = side * k as f64 / segs_per_side as f64;
        ring.push([-half, half - t]);
    }
    ring
}

/// A regular n-gon ring approximating a circle (clockwise when `cw` —
/// the hole-ring orientation [`Domain2D::with_hole`] consumes).
pub fn circle_ring(radius: f64, n: usize, cw: bool) -> Vec<[f64; 2]> {
    (0..n)
        .map(|k| {
            let t = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
            let t = if cw { -t } else { t };
            [radius * t.cos(), radius * t.sin()]
        })
        .collect()
}

/// Distance from a point to the nearest point of the radius-1 rim circle
/// at the origin (the plate family's feature).
pub fn rim_distance(p: [f64; 2]) -> f64 {
    let rho = (p[0] * p[0] + p[1] * p[1]).sqrt();
    if rho <= 0.0 {
        return 1.0;
    }
    (rho - 1.0).abs()
}

/// The analytic rim geometry source behind the graded prescription: a
/// radius-`radius` circle in the z=0 plane answering closest hits
/// everywhere and principal curvatures 1/r within `delta` of the rim,
/// None beyond (the P3a matched-field fake, productionized).
#[derive(Debug, Clone)]
pub struct RimGeometrySource {
    pub radius: f64,
    pub delta: f64,
}

impl GeometrySource for RimGeometrySource {
    fn closest_hit(&self, query: DVec3) -> Option<GeometryHit> {
        let rho = (query.x * query.x + query.y * query.y).sqrt();
        let point = if rho <= 0.0 {
            DVec3::new(self.radius, 0.0, 0.0)
        } else {
            let scale = self.radius / rho;
            DVec3::new(query.x * scale, query.y * scale, 0.0)
        };
        let normal = DVec3::new(point.x / self.radius, point.y / self.radius, 0.0);
        Some(GeometryHit {
            point,
            normal,
            uv: [point.y.atan2(point.x), 0.0],
        })
    }

    fn principal_curvatures(&self, query: DVec3) -> Option<[f64; 2]> {
        let hit = self.closest_hit(query)?;
        let d = query.distance(hit.point);
        if d <= self.delta {
            Some([1.0 / self.radius, 1.0 / self.radius])
        } else {
            None
        }
    }
}

/// The plate-with-hole corpus domain: ±5 square outer (6 segments per
/// side), radius-1 5-gon hole (the P3a corpus domain verbatim — the
/// 5-gon is probed-and-kept, a 7-gon rim worsened the graded p95).
pub fn plate_with_hole_domain() -> Domain2D {
    Domain2D::from_outer(square_ring(5.0, 6)).with_hole(circle_ring(1.0, 5, true))
}

/// The concave L polygon (area 3.0) for the CDT family.
pub fn l_shape_polygon() -> Vec<[f64; 2]> {
    vec![
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
    ]
}

/// A quad-faced axis-aligned box surface (the 3-D families' input form).
pub fn box_surface(lo: [f64; 3], hi: [f64; 3]) -> Mesh {
    let mut mesh = Mesh::new();
    let corner_ids = [
        (1u64, [lo[0], lo[1], lo[2]]),
        (2, [hi[0], lo[1], lo[2]]),
        (3, [hi[0], hi[1], lo[2]]),
        (4, [lo[0], hi[1], lo[2]]),
        (5, [lo[0], lo[1], hi[2]]),
        (6, [hi[0], lo[1], hi[2]]),
        (7, [hi[0], hi[1], hi[2]]),
        (8, [lo[0], hi[1], hi[2]]),
    ];
    for (id, xyz) in corner_ids {
        mesh.add_node(Node::new(id, xyz[0], xyz[1], xyz[2]));
    }
    // outward-oriented quad faces (the multi-material classification and
    // the 3-D meshers consume the orientation)
    let faces: [(u64, Vec<u64>); 6] = [
        (1, vec![1, 4, 3, 2]), // z = lo (outward -z)
        (2, vec![5, 6, 7, 8]), // z = hi (outward +z)
        (3, vec![1, 2, 6, 5]), // y = lo (outward -y)
        (4, vec![2, 3, 7, 6]), // x = hi (outward +x)
        (5, vec![3, 4, 8, 7]), // y = hi (outward +y)
        (6, vec![4, 1, 5, 8]), // x = lo (outward -x)
    ];
    for (id, nodes) in faces {
        mesh.add_element(Element::new(id, ElementType::Quad4, nodes));
    }
    mesh
}

/// The nested-cubes multi-material fixture: unit outer cube, side-0.75
/// cavity at [0.1, 0.85]³ (volumes 0.578125 / 0.421875 exactly).
pub fn nested_cubes_surfaces() -> (Mesh, Mesh) {
    (
        box_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        box_surface([0.1, 0.1, 0.1], [0.85, 0.85, 0.85]),
    )
}

// ─── quality statistics ───────────────────────────────────────────────────────

/// The 95th percentile of `values` (nearest-rank).
pub fn p95(mut values: Vec<f64>) -> f64 {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if values.is_empty() {
        return f64::INFINITY;
    }
    let idx = ((values.len() - 1) as f64 * 0.95).round() as usize;
    values[idx]
}

/// 2-D triangle-mesh quality: count, worst interior angle, the p95
/// circumradius/shortest-edge ratio (equilateral = 1/sqrt(3) ~ 0.577 —
/// the P3a Gmsh-comparison measure), and the total area.
#[derive(Debug, Clone, PartialEq)]
pub struct TriQualityStats {
    pub triangle_count: usize,
    pub min_angle_deg: f64,
    pub p95_radius_edge: f64,
    pub total_area: f64,
}

/// Interior angle at vertex `p` in degrees.
fn angle_deg(p: [f64; 2], q: [f64; 2], r: [f64; 2]) -> f64 {
    let u = [q[0] - p[0], q[1] - p[1]];
    let v = [r[0] - p[0], r[1] - p[1]];
    let lu = (u[0] * u[0] + u[1] * u[1]).sqrt();
    let lv = (v[0] * v[0] + v[1] * v[1]).sqrt();
    if lu < 1e-15 || lv < 1e-15 {
        return 0.0;
    }
    ((u[0] * v[0] + u[1] * v[1]) / (lu * lv))
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}

/// Triangle-mesh quality over the `Triangle3` skeleton (degenerate
/// triangles count at angle 0 / infinite ratio — they are exactly what
/// the gate must catch).
pub fn tri_quality_stats(mesh: &Mesh) -> TriQualityStats {
    let pos = |id: u64| {
        let n = mesh.nodes.get(&id).expect("triangle node exists");
        [n.position.x, n.position.y]
    };
    let mut min_angle = f64::INFINITY;
    let mut ratios = Vec::new();
    let mut total_area = 0.0;
    let mut count = 0usize;
    for elem in &mesh.elements {
        if elem.etype != ElementType::Triangle3 || elem.node_ids.len() != 3 {
            continue;
        }
        let a = pos(elem.node_ids[0]);
        let b = pos(elem.node_ids[1]);
        let c = pos(elem.node_ids[2]);
        count += 1;
        min_angle = min_angle
            .min(angle_deg(a, b, c))
            .min(angle_deg(b, c, a))
            .min(angle_deg(c, a, b));
        let area2 = ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs();
        total_area += area2 * 0.5;
        let la = ((b[0] - c[0]).powi(2) + (b[1] - c[1]).powi(2)).sqrt();
        let lb = ((c[0] - a[0]).powi(2) + (c[1] - a[1]).powi(2)).sqrt();
        let lc = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt();
        let l_min = la.min(lb).min(lc).max(1e-15);
        let circumradius = (la * lb * lc) / (2.0 * area2.max(1e-15));
        ratios.push(circumradius / l_min);
    }
    TriQualityStats {
        triangle_count: count,
        min_angle_deg: if min_angle.is_finite() { min_angle } else { 0.0 },
        p95_radius_edge: p95(ratios),
        total_area,
    }
}

/// 3-D tetrahedron-mesh quality (the mesher3d regression definitions,
/// productionized): worst dihedral, p95 radius-edge (circumradius over
/// shortest edge), the sliver fraction (dihedral < 6° AND ratio > 1.8),
/// and the tet count.
#[derive(Debug, Clone, PartialEq)]
pub struct TetQualityStats {
    pub min_dihedral_deg: f64,
    pub p95_radius_edge: f64,
    pub sliver_fraction: f64,
    pub tet_count: usize,
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
            let n1 = cross3(sub3(q, p), sub3(r, p));
            let n2 = cross3(sub3(q, p), sub3(s, p));
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
    dot3(sub3(a, d), cross3(sub3(b, d), sub3(c, d))).abs() / 6.0
}

/// The tetrahedron circumradius over its shortest edge.
pub fn tet_radius_edge_ratio(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let mut m = [[0.0; 4]; 3];
    let rows = [
        (sub3(b, a), 0.5 * (dot3(b, b) - dot3(a, a))),
        (sub3(c, a), 0.5 * (dot3(c, c) - dot3(a, a))),
        (sub3(d, a), 0.5 * (dot3(d, d) - dot3(a, a))),
    ];
    for i in 0..3 {
        m[i][0] = rows[i].0[0];
        m[i][1] = rows[i].0[1];
        m[i][2] = rows[i].0[2];
        m[i][3] = rows[i].1;
    }
    for col in 0..3 {
        let mut pivot = col;
        for r in (col + 1)..3 {
            if m[r][col].abs() > m[pivot][col].abs() {
                pivot = r;
            }
        }
        if m[pivot][col].abs() < 1e-15 {
            return f64::INFINITY;
        }
        if pivot != col {
            m.swap(pivot, col);
        }
        let inv = 1.0 / m[col][col];
        for j in 0..4 {
            m[col][j] *= inv;
        }
        for r in 0..3 {
            if r != col {
                let f = m[r][col];
                for j in 0..4 {
                    m[r][j] -= f * m[col][j];
                }
            }
        }
    }
    let cc = [m[0][3], m[1][3], m[2][3]];
    let edges = [
        len3(sub3(a, b)),
        len3(sub3(a, c)),
        len3(sub3(a, d)),
        len3(sub3(b, c)),
        len3(sub3(b, d)),
        len3(sub3(c, d)),
    ];
    let min_edge = edges.iter().copied().fold(f64::MAX, f64::min).max(1e-15);
    len3(sub3(cc, a)) / min_edge
}

/// Tetrahedron-mesh quality over the `Tetrahedron4` skeleton (zero-
/// volume tets are skipped — they cannot be quality-sampled — every
/// other tet counts).
pub fn tet_quality_stats(mesh: &Mesh) -> TetQualityStats {
    let pos = |id: u64| {
        let n = mesh.nodes.get(&id).expect("tet node exists");
        [n.position.x, n.position.y, n.position.z]
    };
    let mut min_dihedral = f64::INFINITY;
    let mut ratios = Vec::new();
    let mut slivers = 0usize;
    let mut total = 0usize;
    for elem in &mesh.elements {
        if elem.etype != ElementType::Tetrahedron4 || elem.node_ids.len() != 4 {
            continue;
        }
        let a = pos(elem.node_ids[0]);
        let b = pos(elem.node_ids[1]);
        let c = pos(elem.node_ids[2]);
        let d = pos(elem.node_ids[3]);
        if tetra_volume(a, b, c, d) < 1e-15 {
            continue;
        }
        total += 1;
        let dihedral = min_dihedral_angle_tet(a, b, c, d);
        let ratio = tet_radius_edge_ratio(a, b, c, d);
        min_dihedral = min_dihedral.min(dihedral);
        ratios.push(ratio);
        if dihedral < 6.0 && ratio > 1.8 {
            slivers += 1;
        }
    }
    TetQualityStats {
        min_dihedral_deg: if min_dihedral.is_finite() { min_dihedral } else { 0.0 },
        p95_radius_edge: p95(ratios),
        sliver_fraction: if total == 0 {
            1.0
        } else {
            slivers as f64 / total as f64
        },
        tet_count: total,
    }
}

// ─── gates ────────────────────────────────────────────────────────────────────

/// A quality gate over [`TetQualityStats`]. The ROADMAP hard gate is the
/// post-optimization target; families that do not meet it yet pin their
/// measured baseline as the regression bound and report the distance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QualityGate {
    pub min_dihedral_deg: f64,
    pub p95_radius_edge: f64,
    pub max_sliver_fraction: f64,
}

impl QualityGate {
    /// The ROADMAP §1 post-optimization hard gate (5° / 2.0 / 0.1%).
    pub fn roadmap_hard_gate() -> Self {
        QualityGate {
            min_dihedral_deg: 5.0,
            p95_radius_edge: 2.0,
            max_sliver_fraction: 1.0e-3,
        }
    }

    /// The violated bounds, empty when the gate holds.
    pub fn violations(&self, stats: &TetQualityStats) -> Vec<String> {
        let mut out = Vec::new();
        if !(stats.min_dihedral_deg >= self.min_dihedral_deg) {
            out.push(format!(
                "min dihedral {:.3} deg < gate {:.3}",
                stats.min_dihedral_deg, self.min_dihedral_deg
            ));
        }
        if !(stats.p95_radius_edge <= self.p95_radius_edge) {
            out.push(format!(
                "p95 radius-edge {:.3} > gate {:.3}",
                stats.p95_radius_edge, self.p95_radius_edge
            ));
        }
        if !(stats.sliver_fraction <= self.max_sliver_fraction) {
            out.push(format!(
                "sliver fraction {:.5} > gate {:.5}",
                stats.sliver_fraction, self.max_sliver_fraction
            ));
        }
        out
    }
}

/// One-line quality report for a run log or benchmark table.
pub fn format_tet_stats(name: &str, stats: &TetQualityStats) -> String {
    format!(
        "{name}: tets={}, min_dihedral={:.3} deg, p95_radius_edge={:.3}, sliver_frac={:.5}",
        stats.tet_count,
        stats.min_dihedral_deg,
        stats.p95_radius_edge,
        stats.sliver_fraction
    )
}

/// One-line quality report for a 2-D family.
pub fn format_tri_stats(name: &str, stats: &TriQualityStats) -> String {
    format!(
        "{name}: tris={}, min_angle={:.3} deg, p95_radius_edge={:.3}, area={:.6}",
        stats.triangle_count,
        stats.min_angle_deg,
        stats.p95_radius_edge,
        stats.total_area
    )
}
