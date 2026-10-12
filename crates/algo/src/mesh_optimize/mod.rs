//! 3-D volume mesh quality optimization: measures + passes over the compact
//! tet mesh — the 3-D sibling of [`crate::mesh_quality_2d`], which mirrors its
//! structure (pure quality measures, then acceptance-gated passes, then a
//! round-based driver).
//!
//! # Quality measures (pure functions, per tetrahedron)
//!
//! * **Minimum dihedral angle** ([`min_dihedral_deg_tet`]): the minimum over
//!   the 6 edges of the angle between the two adjacent face normals, in
//!   degrees. Range `(0°, 70.5288°]`; the regular tetrahedron has
//!   `arccos(1/3) ≈ 70.53°` at every edge.
//! * **Radius-edge ratio** ([`radius_edge_ratio_tet`]): `R / l_min`, the
//!   circumsphere radius over the shortest edge (Chew's `q` measure;
//!   Parthasarathy-Graichen-Haber element classification). The regular
//!   tetrahedron has `R = sqrt(6)/4 · l`, so `q = sqrt(6)/4 ≈ 0.6124`
//!   (the skeleton's earlier "≈ 1.22" doc value was the *diameter*-edge
//!   ratio `2R/l`; corrected here). Degenerate (co-spherical / zero-edge)
//!   input yields `+∞`.
//! * **Sliver classification** ([`is_sliver_tet`]): the acceptance-corpus
//!   definition pinned by `tests/mesher3d_quality_regression.rs` — a tet is
//!   sliver-like iff `min_dihedral < 6° && radius-edge > 1.8` (the classic
//!   Cheng-Dey-Edelsbrunner sliver family: near-zero dihedral with bounded
//!   edge lengths, which the plain needle / cap families do not cover).
//! * **Signed volume** ([`signed_volume_tet`]): `((a-d) · ((b-d) × (c-d))) / 6`
//!   — positive for the stored orientation, used as the inversion guard in
//!   every pass (no pass may create a non-positive tetrahedron).
//! * **Scaled Jacobian** ([`scaled_jacobian_tet`]): Knupp's `min` over the 4
//!   vertices of `det(J_v) / (|e1||e2||e3|)` = `6V / (|e1||e2||e3|)` using the
//!   three edges at the vertex. Regular tet: `1/sqrt(2) ≈ 0.7071`.
//! * **Aspect ratio** ([`aspect_ratio_tet`]): longest edge / (2 · inradius),
//!   with `r_in = 3V / A_total`. Regular tet: `sqrt(6) ≈ 2.449`.
//!
//! # Passes (acceptance-gated, mirroring the 2-D module)
//!
//! * **2-3 edge swap** ([`edge_swap_pass`]): every interior face shared by
//!   exactly two tets `(a,b,c,d)` + `(a,b,e,f)` (face `(a,b,c)`-pattern, two
//!   apexes `p`, `q`) is swapped for the opposite edge `(p,q)` plus three
//!   tets, when the swap strictly improves the affected tets' lexicographic
//!   quality tuple `(worst [`QualityMetric`] goodness, worst radius-edge
//!   goodness)` — or passes the sliver-relaxed gate (>= 30% worst-ratio
//!   reduction with the primary degraded by at most 30%, the
//!   `sliver_relaxed_improves` pattern of `crate::tet_mesh`). The radius-edge
//!   tie-break keeps driving the charter `p95 radius-edge` gate after the
//!   primary metric converges (the codebase's `is_better_quality`
//!   comparison pattern). Validity gate: the
//!   new edge (the segment between the two apexes) must strictly pierce the
//!   interior of the removed shared face — i.e. it lies inside the union of
//!   the two tets. This direct segment-triangle test is the local
//!   (non-inverting) validity condition of the bistellar 2-3 flip; it is
//!   equivalent to the in-circumsphere predicate formulation (apex `q` inside
//!   the circumsphere of the opposite tet) for a convex pair, and we use the
//!   direct geometric form (George & Borouchaki, "Back to Edge Flips in 3
//!   Dimensions", IMR 2003).
//! * **3-2 edge flip** ([`edge_removal_pass_3_2`]): every interior edge
//!   shared by exactly three tets is removed, replacing the ring with the
//!   pair of tets over the opposite triangle, under the same acceptance
//!   gates. This completes the bistellar flip family — the 2-3 pass alone
//!   cannot remove edges.
//! * **Node relocation** ([`relocate_pass`]): every free (non-fixed) node is
//!   pulled toward the **volume-weighted centroid** of its star neighbours
//!   (weight of neighbour `n` = total volume of the star tets containing edge
//!   `(v,n)`) — the 3-D analog of Gmsh's weighted relocation, where element
//!   volume replaces the 2-D surface-area weighting. Adjudication:
//!   volume-weighted rather than uniform average, because volume is the
//!   natural 3-D element measure and damps refined-region nodes from being
//!   dragged across coarse regions; the plain-average 2-D behaviour is
//!   recovered on a uniform star. Mirroring the 2-D `_relocate` ladder, the
//!   move is tried at five damped fractions (`FACTOR = 1, 1/1.4, 1/1.96, …`)
//!   and a damped step is accepted only when no star tet inverts (signed
//!   volume stays positive), the total star volume does not grow (relative
//!   tolerance `1e-10`, mirroring the 2-D surface non-growth check), and the
//!   star's lexicographic quality tuple (primary metric, radius-edge
//!   tie-break) holds or improves. Accepted steps accumulate.
//! * **Edge collapse** ([`edge_collapse_pass`]): edges shorter than the
//!   collapse threshold collapse into one endpoint. Endpoint choice: when
//!   exactly one endpoint is a boundary node, the boundary endpoint survives
//!   (the domain boundary is preserved exactly); otherwise the higher-valence
//!   endpoint survives (valence = number of incident tets; tie broken by the
//!   lower compact index). Both-boundary edges are never collapsed (v1
//!   limitation, documented). Validity gates: the de Cougny-Shephard link
//!   condition in vertex-link form (the common neighbour set of the two
//!   endpoints must equal exactly the apex set of the edge ring), and a
//!   geometric simulation of the merge — every surviving affected tet must
//!   keep positive signed volume and a minimum dihedral of at least
//!   `floor = clamp(threshold/10, 0.5°, 5°)`.
//! * **Node insertion** — deliberately NOT part of the v1 volume optimizer.
//!   The skeleton's "insert into the worst element" todo is resolved to a
//!   documented total error path: interior-point insertion cannot repair the
//!   families that need it. Centroid splits worsen every probed family (cap
//!   / sliver / needle / wedge / spike / spindle: each sub-tet scores worse
//!   than the parent on both axes), and a bad tet's circumcenter lies
//!   outside it (interior circumcenters require all dihedral angles < 90° —
//!   already-decent tets), so a quality-driven insertion is a no-op by
//!   construction. Enabling `OptimizeConfig::node_insertion` raises
//!   [`MeshAlgoError::InvalidInput`]; surface-conforming refinement belongs
//!   to a future refinement card.
//!
//! # Driver
//!
//! [`MeshQualityOptimizer::optimize_with_params`] runs rounds of
//! **collapses → swaps (2-3, then 3-2) → relocations** for
//! `MeshParams::optimize_passes` rounds (trait entry: `OptimizeParams::iterations`),
//! stopping early when a round accepts no operation or the global worst
//! quality tuple (primary against `OptimizeParams::tolerance`, secondary
//! against `tolerance * 1e-3`) stops improving.
//!
//! Input contract: the mesh must be a pure tetrahedral volume mesh
//! (`Tetrahedron4` elements only — anything else is a documented
//! [`MeshAlgoError::InvalidInput`]); zero-volume tets are a documented
//! [`MeshAlgoError::Generation`] error; stored counter-clockwise-winding tets
//! are re-oriented on entry (count reported through normalization). Node IDs
//! are preserved; nodes orphaned by collapse are dropped on write-back;
//! element IDs are renumbered `1..n`.
//!
//! # Reference
//!
//! Gmsh source: `Mesh/qualityMeasures.cpp` (`qmTetrahedron::qm`),
//! `Mesh/meshGRegionDelaunayInsertion.cpp`.
//! V. Parthasarathy, C. Graichen, K. Haber, "A comparison of element quality
//! measures", *Int. J. Numer. Methods Eng.* 36, 1993.
//! H. de Cougny, M. Shephard, "Parallel refinement and coarsening of
//! tetrahedral meshes", *Int. J. Numer. Methods Eng.* 41, 1998 (link
//! condition).
//! P.-L. George, H. Borouchaki, "Back to Edge Flips in 3 Dimensions",
//! *Proc. 12th Int. Meshing Roundtable*, 2003.
//! C. Geuzaine, J.-F. Remacle, "Gmsh", *Int. J. Numer. Methods Eng.* 79,
//! 2009 (the measures + acceptance-gated passes scheme follows the Gmsh
//! optimizer architecture; reimplemented, not translated).

use std::collections::{BTreeMap, HashMap, HashSet};

use rmsh_model::{Element, ElementType, Mesh, Node};

use crate::tet_mesh::{Tet, TetMesh};
use crate::traits::{MeshAlgoError, MeshOptimizer, MeshParams, OptimizeParams};

// ─── Quality measures ─────────────────────────────────────────────────────────

/// Absolute volume epsilon for "zero-volume / degenerate" classification.
const VOL_EPS: f64 = 1e-15;

/// Signed volume of a tetrahedron: `((a-d) · ((b-d) × (c-d))) / 6`.
///
/// Positive for one orientation, negative for the reverse; the inversion
/// guard of every pass.
pub fn signed_volume_tet(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let ad = sub3(a, d);
    let bd = sub3(b, d);
    let cd = sub3(c, d);
    dot3(ad, cross3(bd, cd)) / 6.0
}

/// Dihedral angle at edge `(p, q)` between the faces `(p,q,r)` and `(p,q,s)`
/// (degrees). Zero for degenerate (collinear) input.
fn dihedral_angle(p: [f64; 3], q: [f64; 3], r: [f64; 3], s: [f64; 3]) -> f64 {
    let pq = sub3(q, p);
    let pr = sub3(r, p);
    let ps = sub3(s, p);
    let n1 = cross3(pq, pr);
    let n2 = cross3(pq, ps);
    let dot = dot3(n1, n2);
    let l1 = len3(n1);
    let l2 = len3(n2);
    if l1 < 1e-15 || l2 < 1e-15 {
        return 0.0;
    }
    (dot / (l1 * l2)).clamp(-1.0, 1.0).acos().to_degrees()
}

/// Minimum dihedral angle of a tetrahedron over its 6 edges (degrees).
pub fn min_dihedral_deg_tet(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
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
        .map(|&(p, q, r, s)| dihedral_angle(p, q, r, s))
        .fold(f64::MAX, f64::min)
}

/// Solve a 3x3 linear system by Gaussian elimination with partial pivoting.
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

/// Circumsphere center of a tetrahedron (`None` for degenerate input).
fn tetra_circumcenter(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Option<[f64; 3]> {
    let ba = sub3(b, a);
    let ca = sub3(c, a);
    let da = sub3(d, a);
    let rhs0 = 0.5 * (dot3(b, b) - dot3(a, a));
    let rhs1 = 0.5 * (dot3(c, c) - dot3(a, a));
    let rhs2 = 0.5 * (dot3(d, d) - dot3(a, a));
    solve_3x3([(ba, rhs0), (ca, rhs1), (da, rhs2)])
}

/// Radius-edge ratio `R / l_min`: circumsphere radius over the shortest edge.
///
/// `+INF` when the circumsphere does not exist (degenerate tet) or the
/// shortest edge is zero.
pub fn radius_edge_ratio_tet(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let pts = [a, b, c, d];
    let mut lmin = f64::MAX;
    for i in 0..4 {
        for j in (i + 1)..4 {
            lmin = lmin.min(len3(sub3(pts[i], pts[j])));
        }
    }
    if lmin < VOL_EPS {
        return f64::INFINITY;
    }
    match tetra_circumcenter(a, b, c, d) {
        Some(cc) => len3(sub3(cc, a)) / lmin,
        None => f64::INFINITY,
    }
}

/// Sliver classification per the acceptance-corpus definition pinned by
/// `tests/mesher3d_quality_regression.rs`: sliver-like iff
/// `min_dihedral < 6 deg && radius-edge > 1.8`.
pub fn is_sliver_tet(min_dihedral_deg: f64, radius_edge: f64) -> bool {
    min_dihedral_deg < 6.0 && radius_edge > 1.8
}

/// Scaled Jacobian (Knupp): the minimum over the 4 vertices of
/// `det(J_v) / (|e1||e2||e3|)` = `6V / (|e1||e2||e3|)` with the three edges
/// at that vertex. Zero for degenerate input.
pub fn scaled_jacobian_tet(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let vol6 = signed_volume_tet(a, b, c, d).abs() * 6.0;
    if vol6 <= VOL_EPS {
        return 0.0;
    }
    let pts = [a, b, c, d];
    let mut sj = f64::MAX;
    for v in 0..4 {
        let mut lens = [0.0; 3];
        for (k, w) in [(0, 1), (0, 2), (0, 3)].iter().enumerate() {
            let i = (v + w.0) % 4;
            let j = (v + w.1) % 4;
            lens[k] = len3(sub3(pts[i], pts[j]));
        }
        if lens.iter().any(|&l| l < VOL_EPS) {
            return 0.0;
        }
        sj = sj.min(vol6 / (lens[0] * lens[1] * lens[2]));
    }
    sj
}

/// Aspect ratio: longest edge / (2 · inradius), `r_in = 3V / A_total`.
/// `+INF` for degenerate input.
pub fn aspect_ratio_tet(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let vol = signed_volume_tet(a, b, c, d).abs();
    if vol <= VOL_EPS {
        return f64::INFINITY;
    }
    let area = triangle_area3d(a, b, c)
        + triangle_area3d(a, b, d)
        + triangle_area3d(a, c, d)
        + triangle_area3d(b, c, d);
    if area <= VOL_EPS {
        return f64::INFINITY;
    }
    let r_in = 3.0 * vol / area;
    let pts = [a, b, c, d];
    let mut lmax: f64 = 0.0;
    for i in 0..4 {
        for j in (i + 1)..4 {
            lmax = lmax.max(len3(sub3(pts[i], pts[j])));
        }
    }
    lmax / (2.0 * r_in)
}

/// Area of a 3-D triangle (half the cross-product norm).
fn triangle_area3d(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    len3(cross3(sub3(b, a), sub3(c, a))) * 0.5
}

/// The quality score of one tetrahedron, in "goodness" polarity: **higher is
/// always better**, so badness-type measures enter negated.
fn tet_goodness_points(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    metric: QualityMetric,
) -> f64 {
    match metric {
        QualityMetric::MinAngle => min_dihedral_deg_tet(a, b, c, d),
        QualityMetric::RadiusEdgeRatio => -radius_edge_ratio_tet(a, b, c, d),
        QualityMetric::ScaledJacobian => scaled_jacobian_tet(a, b, c, d),
        QualityMetric::AspectRatio => -aspect_ratio_tet(a, b, c, d),
    }
}

/// Lexicographic quality tuple of an affected tet set: `(worst primary
/// goodness, worst secondary goodness)` where the secondary is the negative
/// of the worst (largest) radius-edge ratio — unless the metric is itself
/// radius-edge, in which case the secondary is the minimum dihedral. The
/// radius-edge tie-break drives the charter `p95 radius-edge` gate once the
/// primary metric converges. `None` when any tet is non-positive-volume.
///
/// Comparison mirrors the codebase's `is_better_quality` pattern
/// (`crates/algo/src/delaunay_3d.rs`): primary with `1e-6` slack, secondary
/// with `1e-9` slack.
fn affected_quality(tm: &TetMesh, tets: &[[u32; 4]], metric: QualityMetric) -> Option<(f64, f64)> {
    let mut primary = f64::MAX;
    let mut secondary = f64::MAX;
    for nodes in tets {
        let pts = nodes.map(|n| tm.nodes[n as usize]);
        if signed_volume_tet(pts[0], pts[1], pts[2], pts[3]) <= VOL_EPS {
            return None;
        }
        primary = primary.min(tet_goodness_points(
            pts[0], pts[1], pts[2], pts[3], metric,
        ));
        secondary = secondary.min(match metric {
            QualityMetric::RadiusEdgeRatio => min_dihedral_deg_tet(pts[0], pts[1], pts[2], pts[3]),
            _ => -radius_edge_ratio_tet(pts[0], pts[1], pts[2], pts[3]),
        });
    }
    Some((primary, secondary))
}

/// Strict lexicographic improvement of a quality tuple.
fn quality_tuple_better(new: (f64, f64), old: (f64, f64)) -> bool {
    (new.0 > old.0 + 1e-6) || ((new.0 - old.0).abs() <= 1e-6 && new.1 < old.1 - 1e-9)
}

/// Hold-or-improve: the primary may not degrade (beyond slack), and when it
/// is not strictly improving the secondary may not degrade either.
fn quality_tuple_holds(new: (f64, f64), old: (f64, f64)) -> bool {
    new.0 >= old.0 - 1e-6 && (new.0 > old.0 + 1e-6 || new.1 <= old.1 + 1e-9)
}

/// Sliver-relaxed acceptance, mirroring the `sliver_relaxed_improves`
/// pattern of `crate::tet_mesh::optimize_tetmesh_flips`: the secondary
/// (radius-edge) must improve strongly (>= 30% reduction of the worst
/// ratio) while the primary degrades by at most 30% of its magnitude. This
/// unlocks flips over corner-to-corner boundary diagonals that no strictly
/// improving local move can touch.
fn quality_tuple_sliver_relaxed(new: (f64, f64), old: (f64, f64)) -> bool {
    if old.1 >= f64::MAX / 2.0 || new.1 >= f64::MAX / 2.0 {
        return false;
    }
    // secondary holds -max_ratio: improvement means new.1 > old.1
    let ratio_old = -old.1;
    let ratio_new = -new.1;
    ratio_new < ratio_old * 0.7 && new.0 >= old.0 - 0.3 * old.0.abs()
}

// ─── Vector helpers ───────────────────────────────────────────────────────────

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

// ─── Aggregate stats ──────────────────────────────────────────────────────────

/// Aggregate quality statistics of a tetrahedral mesh, computed with the same
/// formulas (and the same p95 / sliver definitions) as
/// `tests/mesher3d_quality_regression.rs`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TetMeshQualityStats {
    /// Minimum dihedral angle over all tets (degrees).
    pub min_dihedral_deg: f64,
    /// 95th percentile of the radius-edge ratio distribution.
    pub p95_radius_edge: f64,
    /// Fraction of sliver-like tets (`is_sliver_tet`).
    pub sliver_fraction: f64,
    /// Number of valid (positive-volume) tets.
    pub tet_count: usize,
}

/// Compute [`TetMeshQualityStats`] over the tetrahedral elements of `mesh`.
pub fn mesh_quality_stats(mesh: &Mesh) -> TetMeshQualityStats {
    let mut min_dihedral = f64::INFINITY;
    let mut ratios = Vec::<f64>::new();
    let mut slivers = 0usize;
    let mut total = 0usize;

    for e in &mesh.elements {
        if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
            continue;
        }
        let ids = [
            e.node_ids[0],
            e.node_ids[1],
            e.node_ids[2],
            e.node_ids[3],
        ];
        let Some(pts) = tet_points_by_id(mesh, &ids) else {
            continue;
        };
        let volume = signed_volume_tet(pts[0], pts[1], pts[2], pts[3]).abs();
        if volume < VOL_EPS {
            continue;
        }
        total += 1;
        let dihedral = min_dihedral_deg_tet(pts[0], pts[1], pts[2], pts[3]);
        let ratio = radius_edge_ratio_tet(pts[0], pts[1], pts[2], pts[3]);
        min_dihedral = min_dihedral.min(dihedral);
        ratios.push(ratio);
        if is_sliver_tet(dihedral, ratio) {
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

    TetMeshQualityStats {
        min_dihedral_deg: if min_dihedral.is_finite() {
            min_dihedral
        } else {
            0.0
        },
        p95_radius_edge: p95,
        sliver_fraction: if total == 0 {
            1.0
        } else {
            slivers as f64 / total as f64
        },
        tet_count: total,
    }
}

/// Fetch the coordinates of node IDs; `None` when any ID is missing.
fn tet_points_by_id(mesh: &Mesh, ids: &[u64; 4]) -> Option<[[f64; 3]; 4]> {
    let mut pts = [[0.0; 3]; 4];
    for (k, id) in ids.iter().enumerate() {
        let node = mesh.nodes.get(id)?;
        pts[k] = [node.position.x, node.position.y, node.position.z];
    }
    Some(pts)
}

// ─── Metric / config ──────────────────────────────────────────────────────────

/// The quality measure used to score elements and guide optimization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QualityMetric {
    /// Minimum dihedral angle over the 6 edges of the tet, in degrees.
    /// Range `(0°, 70.5288°]` for the regular tetrahedron. (For 2-D triangles
    /// the corresponding measure set lives in [`crate::mesh_quality_2d`].)
    #[default]
    MinAngle,

    /// Radius-edge ratio `R / l_min` (3-D only).
    ///
    /// Regular tet: `sqrt(6)/4 ≈ 0.6124`. Degenerate tet → `+∞`. Badness-type:
    /// lower is better.
    RadiusEdgeRatio,

    /// Scaled Jacobian (Knupp): min over vertices of `6V/(|e1||e2||e3|)`.
    /// Regular tet: `1/sqrt(2) ≈ 0.7071`. Inverted element < 0 (0 after the
    /// driver's degenerate guard).
    ScaledJacobian,

    /// Aspect ratio: longest edge / inscribed sphere diameter.
    /// Regular tet: `sqrt(6) ≈ 2.449`. Badness-type: lower is better.
    AspectRatio,
}

// ─── Enabled operations ───────────────────────────────────────────────────────

/// Controls which local mesh-modification operators are active and the
/// thresholds they use. See the module docs for the pass contracts.
#[derive(Debug, Clone)]
pub struct OptimizeConfig {
    /// Quality measure to maximise.
    pub metric: QualityMetric,

    /// Allow 2-3 edge swaps across interior faces (bistellar flips).
    pub edge_swap: bool,

    /// Allow quality-gated node relocation (volume-weighted Laplacian).
    pub laplacian_smooth: bool,

    /// Interior-point insertion into bad elements. NOT implemented in the
    /// v1 quality optimizer (see the module docs for the adjudication):
    /// centroid splits provably worsen every probed family and a bad tet's
    /// circumcenter lies outside it, so enabling this flag raises
    /// [`MeshAlgoError::InvalidInput`].
    pub node_insertion: bool,

    /// Allow edge collapse of short edges.
    pub edge_collapse: bool,

    /// Quality cutoff in the natural units of [`Self::metric`]: the collapse
    /// quality floor derives from it (`clamp(threshold/10, 0.5°, 5°)`).
    /// For `MinAngle` this is the dihedral angle in degrees (default 20°:
    /// tets with min dihedral below 20° count as poor). For `RadiusEdgeRatio`
    /// / `AspectRatio` it is the badness cutoff (above = poor); for
    /// `ScaledJacobian` it is in scaled-Jacobian units.
    pub threshold: f64,

    /// Edge-length collapse threshold. `0.0` (default) = automatic:
    /// `MeshParams::min_size` when mesh params are available via
    /// [`MeshQualityOptimizer::optimize_with_params`], otherwise `0.5 ×` the
    /// mean edge length of the mesh.
    pub collapse_threshold: f64,
}

impl Default for OptimizeConfig {
    fn default() -> Self {
        Self {
            metric: QualityMetric::MinAngle,
            edge_swap: true,
            laplacian_smooth: true,
            node_insertion: false,
            edge_collapse: true,
            threshold: 20.0, // degrees
            collapse_threshold: 0.0,
        }
    }
}

// ─── Run statistics ───────────────────────────────────────────────────────────

/// Counters of the operations a single optimize run accepted.
///
/// Every mesh change is accounted for by exactly one counter, so tests can
/// pin the identity
/// `nodes_after == nodes_before - nodes_removed` and
/// `elements_after == elements_before + face_swaps_2_3 - edge_flips_3_2
///  - tets_removed_by_collapse`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct OptimizeStats {
    /// Rounds actually executed (before early stop).
    pub rounds: usize,
    /// Edge collapses accepted (one node removed per collapse).
    pub collapses: usize,
    /// Nodes removed (== collapses: one endpoint per collapse).
    pub nodes_removed: usize,
    /// Tets deleted by collapses (the edge rings).
    pub tets_removed_by_collapse: usize,
    /// 2-3 swaps accepted (net +1 tet each).
    pub face_swaps_2_3: usize,
    /// 3-2 edge flips accepted (net -1 tet each).
    pub edge_flips_3_2: usize,
    /// Nodes whose position changed in relocation passes.
    pub relocations_accepted: usize,
}

// ─── Public struct ────────────────────────────────────────────────────────────

/// General 3-D mesh quality optimizer.
///
/// Runs collapse → swap → relocation rounds (optionally insertion) over the
/// compact tet-mesh representation, every operation acceptance-gated by the
/// configured [`QualityMetric`]. See the module docs for the pass contracts.
#[derive(Debug, Clone)]
pub struct MeshQualityOptimizer {
    /// Configuration controlling which operations are active.
    pub config: OptimizeConfig,
}

impl Default for MeshQualityOptimizer {
    fn default() -> Self {
        Self {
            config: OptimizeConfig::default(),
        }
    }
}

impl MeshQualityOptimizer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_config(mut self, config: OptimizeConfig) -> Self {
        self.config = config;
        self
    }

    /// Optimize `mesh`, taking the round count from
    /// [`MeshParams::optimize_passes`] and the collapse threshold from
    /// [`MeshParams::min_size`] (unless [`OptimizeConfig::collapse_threshold`]
    /// overrides it). Returns the accepted-operation statistics.
    pub fn optimize_with_params(
        &self,
        mesh: &mut Mesh,
        params: &OptimizeParams,
        mesh_params: &MeshParams,
    ) -> Result<OptimizeStats, MeshAlgoError> {
        self.run(mesh, params, Some(mesh_params))
    }

    /// Shared driver. `mesh_params`, when given, supplies the round count
    /// (`optimize_passes`) and the collapse threshold (`min_size`).
    fn run(
        &self,
        mesh: &mut Mesh,
        params: &OptimizeParams,
        mesh_params: Option<&MeshParams>,
    ) -> Result<OptimizeStats, MeshAlgoError> {
        let rounds = match mesh_params {
            Some(mp) => mp.optimize_passes as usize,
            None => params.iterations as usize,
        };
        if rounds == 0 {
            // Documented semantics: a zero round count skips optimization.
            return Ok(OptimizeStats::default());
        }

        // Documented total error path (see module docs): interior-point
        // insertion cannot repair the families that need it — centroid
        // splits provably worsen every probed family (cap / sliver / needle
        // / wedge / spike / spindle) and a bad tet's circumcenter lies
        // outside it, so `node_insertion` refuses instead of pretending.
        if self.config.node_insertion {
            return Err(MeshAlgoError::InvalidInput(
                "node_insertion is not part of the v1 volume optimizer \
                 (collapse / swap / relocate only)"
                    .to_string(),
            ));
        }
        // Input contract: pure tetrahedral volume mesh.
        if mesh.elements.is_empty() {
            return Err(MeshAlgoError::InvalidInput(
                "volume optimizer requires tetrahedral elements; mesh has none".to_string(),
            ));
        }
        for e in &mesh.elements {
            if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
                return Err(MeshAlgoError::InvalidInput(format!(
                    "volume optimizer requires a pure Tetrahedron4 mesh; element {} is {:?}",
                    e.id, e.etype
                )));
            }
        }

        let physical_names = mesh.physical_names.clone();
        let periodic = mesh.periodic.clone();

        let mut tm = TetMesh::from_mesh(mesh);
        normalize_orientations(&mut tm)?;

        let mut stats = OptimizeStats::default();
        let metric = self.config.metric;
        let floor_deg = (self.config.threshold * 0.1).clamp(0.5, 5.0);

        let mut collapse_threshold = self.config.collapse_threshold;
        if collapse_threshold <= 0.0 {
            if let Some(mp) = mesh_params {
                collapse_threshold = mp.min_size;
            }
        }
        if collapse_threshold <= 0.0 {
            collapse_threshold = 0.5 * mean_edge_length(&tm);
        }

        let mut prev = mesh_worst_quality_tuple(&tm, metric);

        for round in 0..rounds {
            // Cumulative counters at round start, for the no-op early stop.
            let before_round = (
                stats.face_swaps_2_3,
                stats.edge_flips_3_2,
                stats.relocations_accepted,
                stats.collapses,
            );

            if self.config.edge_collapse {
                let out = edge_collapse_pass(&mut tm, collapse_threshold, floor_deg);
                stats.collapses += out.collapses;
                stats.tets_removed_by_collapse += out.tets_removed;
            }
            if self.config.edge_swap {
                stats.face_swaps_2_3 += edge_swap_pass(&mut tm, metric);
                stats.edge_flips_3_2 += edge_removal_pass_3_2(&mut tm, metric);
            }
            if self.config.laplacian_smooth {
                // Relocate to per-round stability (the 2-D module's
                // laplaceSmoothing loop: repeat relocate_pass until no node
                // moves, capped at 8 sweeps).
                let free = free_nodes(&tm, params.move_boundary_nodes);
                for _ in 0..8 {
                    let moved = relocate_pass(&mut tm, &free, metric);
                    stats.relocations_accepted += moved;
                    if moved == 0 {
                        break;
                    }
                }
            }

            stats.rounds = round + 1;
            stats.nodes_removed = stats.collapses;

            // Early stop: no accepted operation this round, or the global
            // worst quality tuple stopped improving (primary against
            // `params.tolerance`, secondary against `tolerance * 1e-3` —
            // the same 1e-6 / 1e-9 slack pair as the tuple comparisons).
            let after_round = (
                stats.face_swaps_2_3,
                stats.edge_flips_3_2,
                stats.relocations_accepted,
                stats.collapses,
            );
            let cur = mesh_worst_quality_tuple(&tm, metric);
            let gained = (cur.0 - prev.0, cur.1 - prev.1);
            prev = cur;
            if after_round == before_round
                || (gained.0 < params.tolerance && gained.1 < params.tolerance * 1e-3)
            {
                break;
            }
        }

        let mut rebuilt = tetmesh_to_mesh_referenced(&tm);
        rebuilt.physical_names = physical_names;
        rebuilt.periodic = periodic;
        *mesh = rebuilt;
        Ok(stats)
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl MeshOptimizer for MeshQualityOptimizer {
    fn name(&self) -> &'static str {
        "Mesh Quality Optimizer"
    }

    fn optimize(&self, mesh: &mut Mesh, params: &OptimizeParams) -> Result<(), MeshAlgoError> {
        self.run(mesh, params, None).map(|_| ())
    }
}

// ─── TetMesh plumbing ─────────────────────────────────────────────────────────

/// Re-orient every tet to positive signed volume; error on zero-volume tets.
///
/// Returns the number of tets whose node order was flipped.
pub(crate) fn normalize_orientations(tm: &mut TetMesh) -> Result<usize, MeshAlgoError> {
    let mut flipped = 0usize;
    for tet in &mut tm.tets {
        let a = tm.nodes[tet.nodes[0] as usize];
        let b = tm.nodes[tet.nodes[1] as usize];
        let c = tm.nodes[tet.nodes[2] as usize];
        let d = tm.nodes[tet.nodes[3] as usize];
        let v = signed_volume_tet(a, b, c, d);
        if v < -VOL_EPS {
            tet.nodes.swap(2, 3);
            flipped += 1;
        } else if v <= VOL_EPS {
            return Err(MeshAlgoError::Generation(format!(
                "degenerate (zero-volume) tetrahedron with nodes {:?}",
                tet.nodes
            )));
        }
    }
    Ok(flipped)
}

/// Convert back to a public `Mesh`, dropping orphan nodes (no tet references
/// them — collapse leftovers) and renumbering elements `1..n`. Node IDs are
/// preserved.
fn tetmesh_to_mesh_referenced(tm: &TetMesh) -> Mesh {
    let mut referenced = vec![false; tm.nodes.len()];
    for tet in &tm.tets {
        for &n in &tet.nodes {
            referenced[n as usize] = true;
        }
    }
    let mut mesh = Mesh::new();
    for (i, coord) in tm.nodes.iter().enumerate() {
        if !referenced[i] {
            continue;
        }
        mesh.add_node(Node::new(tm.node_ids[i], coord[0], coord[1], coord[2]));
    }
    for (elem_id, tet) in tm.tets.iter().enumerate() {
        mesh.add_element(Element::new(
            (elem_id + 1) as u64,
            ElementType::Tetrahedron4,
            vec![
                tm.node_ids[tet.nodes[0] as usize],
                tm.node_ids[tet.nodes[1] as usize],
                tm.node_ids[tet.nodes[2] as usize],
                tm.node_ids[tet.nodes[3] as usize],
            ],
        ));
    }
    mesh
}

/// Node → incident tet indices (built once per pass, deterministic order).
fn star_tets(tm: &TetMesh) -> HashMap<u32, Vec<u32>> {
    let mut star: HashMap<u32, Vec<u32>> = HashMap::new();
    for (ti, tet) in tm.tets.iter().enumerate() {
        for &n in &tet.nodes {
            star.entry(n).or_default().push(ti as u32);
        }
    }
    star
}

/// Free nodes for relocation: all referenced nodes, minus the boundary set
/// unless `move_boundary` allows moving boundary nodes.
fn free_nodes(tm: &TetMesh, move_boundary: bool) -> Vec<u32> {
    let mut free: Vec<u32> = Vec::new();
    if move_boundary {
        let mut seen = vec![false; tm.nodes.len()];
        for tet in &tm.tets {
            for &n in &tet.nodes {
                let i = n as usize;
                if !seen[i] {
                    seen[i] = true;
                    free.push(n);
                }
            }
        }
    } else {
        let bnd = boundary_node_set(tm);
        let mut seen = vec![false; tm.nodes.len()];
        for tet in &tm.tets {
            for &n in &tet.nodes {
                let i = n as usize;
                if !seen[i] && !bnd.contains(&n) {
                    seen[i] = true;
                    free.push(n);
                }
            }
        }
    }
    free.sort_unstable();
    free
}

/// Boundary node set: nodes appearing on any face with no adjacent tet.
fn boundary_node_set(tm: &TetMesh) -> HashSet<u32> {
    let mut set = HashSet::new();
    for (ti, tet) in tm.tets.iter().enumerate() {
        for fi in 0..4 {
            if tet.neighbors[fi] == u32::MAX {
                for &n in tm.face_nodes(ti, fi).iter() {
                    set.insert(n);
                }
            }
        }
    }
    set
}

/// Unique edges as sorted `[min, max]` keys in ascending order.
fn tet_edge_map(tm: &TetMesh) -> BTreeMap<[u32; 2], Vec<u32>> {
    let mut map: BTreeMap<[u32; 2], Vec<u32>> = BTreeMap::new();
    for (ti, tet) in tm.tets.iter().enumerate() {
        let n = tet.nodes;
        for (a, b) in [
            (n[0], n[1]),
            (n[0], n[2]),
            (n[0], n[3]),
            (n[1], n[2]),
            (n[1], n[3]),
            (n[2], n[3]),
        ] {
            map.entry([a.min(b), a.max(b)]).or_default().push(ti as u32);
        }
    }
    map
}

/// Mean edge length over the unique edges.
fn mean_edge_length(tm: &TetMesh) -> f64 {
    let map = tet_edge_map(tm);
    if map.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0;
    for (key, _) in &map {
        sum += len3(sub3(tm.nodes[key[0] as usize], tm.nodes[key[1] as usize]));
    }
    sum / map.len() as f64
}

/// Global worst lexicographic quality tuple over all tets — the driver's
/// convergence signal.
fn mesh_worst_quality_tuple(tm: &TetMesh, metric: QualityMetric) -> (f64, f64) {
    let nodes: Vec<[u32; 4]> = tm.tets.iter().map(|t| t.nodes).collect();
    affected_quality(tm, &nodes, metric).unwrap_or((f64::MAX, f64::MAX))
}

/// Nodes sharing a tet with `x`, excluding `x` itself.
fn vertex_neighbors(tm: &TetMesh, x: u32) -> HashSet<u32> {
    let mut set = HashSet::new();
    for tet in &tm.tets {
        if tet.nodes.contains(&x) {
            for &n in &tet.nodes {
                if n != x {
                    set.insert(n);
                }
            }
        }
    }
    set
}

/// Number of tets containing node `x`.
fn count_tets_with_node(tm: &TetMesh, x: u32) -> usize {
    tm.tets.iter().filter(|t| t.nodes.contains(&x)).count()
}

// ─── Pass 1: edge collapse ────────────────────────────────────────────────────

/// Outcome of one collapse pass.
#[derive(Debug, Clone, Copy, Default)]
struct CollapseOutcome {
    collapses: usize,
    tets_removed: usize,
}

/// Collapse every edge shorter than `threshold` whose merge passes the link
/// condition and the geometric gates. Sweeps repeat until no collapse is
/// accepted (cap 16).
fn edge_collapse_pass(tm: &mut TetMesh, threshold: f64, floor_deg: f64) -> CollapseOutcome {
    let mut out = CollapseOutcome::default();
    for _sweep in 0..16 {
        let bnd = boundary_node_set(tm);
        let edge_keys: Vec<[u32; 2]> = tet_edge_map(tm).into_keys().collect();
        let mut accepted = 0usize;

        for [u, v] in edge_keys {
            let len = len3(sub3(tm.nodes[u as usize], tm.nodes[v as usize]));
            if len >= threshold {
                continue;
            }
            // Current edge ring: all tets containing both endpoints. A
            // manifold interior edge has a ring of >= 3 tets.
            let ring: Vec<u32> = tm
                .tets
                .iter()
                .enumerate()
                .filter(|(_, t)| t.nodes.contains(&u) && t.nodes.contains(&v))
                .map(|(i, _)| i as u32)
                .collect();
            if ring.len() < 3 {
                continue;
            }
            let ub = bnd.contains(&u);
            let vb = bnd.contains(&v);
            if ub && vb {
                // Both endpoints on the boundary: never collapsed (v1
                // documented limitation — keeps the domain boundary intact).
                continue;
            }
            // Survivor choice: the boundary endpoint when there is exactly
            // one; otherwise the higher-valence endpoint (tie: lower index).
            let (keep, remove) = if ub {
                (u, v)
            } else if vb {
                (v, u)
            } else {
                let cu = count_tets_with_node(tm, u);
                let cv = count_tets_with_node(tm, v);
                if cu > cv || (cu == cv && u < v) {
                    (u, v)
                } else {
                    (v, u)
                }
            };

            // Link condition (vertex-link form): the common neighbours of the
            // endpoints must be exactly the apex set of the ring.
            let n_keep = vertex_neighbors(tm, keep);
            let n_remove = vertex_neighbors(tm, remove);
            let common: HashSet<u32> = n_keep.intersection(&n_remove).copied().collect();
            let mut apexes: HashSet<u32> = HashSet::new();
            for &ti in &ring {
                for &n in tm.tets[ti as usize].nodes.iter() {
                    if n != keep && n != remove {
                        apexes.insert(n);
                    }
                }
            }
            if common != apexes {
                continue;
            }

            // Geometric simulation: every tet touched by the merge (contains
            // keep or remove, outside the ring) must keep positive volume and
            // at least `floor_deg` minimum dihedral.
            let ring_set: HashSet<u32> = ring.iter().copied().collect();
            let mut ok = true;
            for (ti, t) in tm.tets.iter().enumerate() {
                if ring_set.contains(&(ti as u32)) {
                    continue;
                }
                if t.nodes.contains(&keep) || t.nodes.contains(&remove) {
                    let nodes = t.nodes.map(|n| if n == remove { keep } else { n });
                    let a = tm.nodes[nodes[0] as usize];
                    let b = tm.nodes[nodes[1] as usize];
                    let c = tm.nodes[nodes[2] as usize];
                    let d = tm.nodes[nodes[3] as usize];
                    if signed_volume_tet(a, b, c, d) <= VOL_EPS
                        || min_dihedral_deg_tet(a, b, c, d) < floor_deg
                    {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok {
                continue;
            }

            // Apply: delete the ring (descending for swap_remove safety),
            // retarget `remove` → `keep`, rebuild.
            let mut ring_sorted = ring;
            ring_sorted.sort_unstable();
            ring_sorted.reverse();
            for &ti in &ring_sorted {
                tm.tets.swap_remove(ti as usize);
            }
            for t in tm.tets.iter_mut() {
                for n in t.nodes.iter_mut() {
                    if *n == remove {
                        *n = keep;
                    }
                }
            }
            tm.build_neighbors();
            accepted += 1;
            out.tets_removed += ring_set.len();
        }

        out.collapses += accepted;
        if accepted == 0 {
            break;
        }
    }
    out
}

// ─── Pass 2: 2-3 edge swap ────────────────────────────────────────────────────

/// Validity of the 2-3 flip: the segment between the two apexes must strictly
/// pierce the interior of the removed shared face (the new edge lies inside
/// the union of the two tets). Direct segment-triangle form of the
/// in-circumsphere predicate; `false` for coplanar / grazing configurations.
/// Shared with `crate::tet_mesh::optimize_tetmesh_flips`, whose 2-3 and 3-2
/// candidate loops must apply the same bistellar validity condition.
pub(crate) fn segment_pierces_triangle(p: [f64; 3], q: [f64; 3], tri: [[f64; 3]; 3]) -> bool {
    let n = cross3(sub3(tri[1], tri[0]), sub3(tri[2], tri[0]));
    let ln = len3(n);
    if ln < 1e-15 {
        return false;
    }
    let d0 = dot3(sub3(p, tri[0]), n);
    let d1 = dot3(sub3(q, tri[0]), n);
    if d0 * d1 >= 0.0 {
        return false; // same side (or touching the plane): no strict piercing
    }
    let t = d0 / (d0 - d1);
    let hit = [
        p[0] + t * (q[0] - p[0]),
        p[1] + t * (q[1] - p[1]),
        p[2] + t * (q[2] - p[2]),
    ];
    // Barycentric coordinates of `hit` in the triangle (Cramer form).
    let v0 = sub3(tri[1], tri[0]);
    let v1 = sub3(tri[2], tri[0]);
    let v2 = sub3(hit, tri[0]);
    let d00 = dot3(v0, v0);
    let d01 = dot3(v0, v1);
    let d11 = dot3(v1, v1);
    let d20 = dot3(v2, v0);
    let d21 = dot3(v2, v1);
    let den = d00 * d11 - d01 * d01;
    if den.abs() < 1e-30 {
        return false;
    }
    let bc_u = (d11 * d20 - d01 * d21) / den;
    let bc_v = (d00 * d21 - d01 * d20) / den;
    const INSIDE_EPS: f64 = 1e-12;
    bc_u > INSIDE_EPS && bc_v > INSIDE_EPS && bc_u + bc_v < 1.0 - INSIDE_EPS
}

/// Locate the current tets sharing `face` (as a node set), one containing
/// `apex_a`, the other `apex_b`. Returns their indices or `None`.
fn find_face_pair(
    tm: &TetMesh,
    face: [u32; 3],
    apex_a: u32,
    apex_b: u32,
) -> Option<(usize, usize)> {
    let mut ti = None;
    let mut tj = None;
    for (i, t) in tm.tets.iter().enumerate() {
        let has_face = face.iter().all(|n| t.nodes.contains(n));
        if !has_face {
            continue;
        }
        if t.nodes.contains(&apex_a) {
            ti = Some(i);
        } else if t.nodes.contains(&apex_b) {
            tj = Some(i);
        }
        if ti.is_some() && tj.is_some() {
            break;
        }
    }
    match (ti, tj) {
        (Some(a), Some(b)) => Some((a, b)),
        _ => None,
    }
}

/// 2-3 swap pass: up to 8 sweeps over the interior faces; each face shared by
/// exactly two tets is swapped for its opposite edge when that strictly
/// improves the worst affected-tet goodness.
fn edge_swap_pass(tm: &mut TetMesh, metric: QualityMetric) -> usize {
    let mut swaps = 0usize;
    for _sweep in 0..8 {
        tm.build_neighbors();
        // Candidate snapshot (deterministic: tet index, then face index).
        let mut candidates: Vec<([u32; 3], u32, u32)> = Vec::new();
        for ti in 0..tm.tets.len() {
            for fi in 0..4 {
                let tj = tm.tets[ti].neighbors[fi];
                if tj == u32::MAX || ti as u32 >= tj {
                    continue;
                }
                let Some(fj) = (0..4).find(|&f| tm.tets[tj as usize].neighbors[f] == ti as u32)
                else {
                    continue;
                };
                let apex_a = tm.tets[ti].nodes[fi];
                let apex_b = tm.tets[tj as usize].nodes[fj];
                let face = tm.face_nodes(ti, fi);
                candidates.push((face, apex_a, apex_b));
            }
        }

        let mut any = false;
        for (face, apex_a, apex_b) in candidates {
            if apex_a == apex_b {
                continue;
            }
            // Re-locate the pair in the current (churned) mesh.
            let Some((ti, tj)) = find_face_pair(tm, face, apex_a, apex_b) else {
                continue;
            };
            let pa = tm.nodes[apex_a as usize];
            let pb = tm.nodes[apex_b as usize];
            let tri = [
                tm.nodes[face[0] as usize],
                tm.nodes[face[1] as usize],
                tm.nodes[face[2] as usize],
            ];
            if !segment_pierces_triangle(pa, pb, tri) {
                continue;
            }
            // Build the three replacement tets around edge (apex_a, apex_b),
            // each normalized to positive orientation.
            let mut new_tets = [[0u32; 4]; 3];
            let mut valid = true;
            for k in 0..3 {
                let mut t = [apex_a, apex_b, face[k], face[(k + 1) % 3]];
                let a = tm.nodes[t[0] as usize];
                let b = tm.nodes[t[1] as usize];
                let c = tm.nodes[t[2] as usize];
                let d = tm.nodes[t[3] as usize];
                let vol = signed_volume_tet(a, b, c, d);
                if vol < 0.0 {
                    t.swap(2, 3);
                } else if vol <= VOL_EPS {
                    valid = false;
                    break;
                }
                new_tets[k] = t;
            }
            if !valid {
                continue;
            }
            // Acceptance: strict lexicographic improvement of the affected
            // set, or the sliver-relaxed gate (strong radius-edge reduction
            // with bounded primary regression) for corner-diagonal slivers.
            let old_set = [tm.tets[ti].nodes, tm.tets[tj].nodes];
            let Some(old_q) = affected_quality(tm, &old_set, metric) else {
                continue;
            };
            let Some(new_q) = affected_quality(tm, &new_tets, metric) else {
                continue;
            };
            if quality_tuple_better(new_q, old_q)
                || quality_tuple_sliver_relaxed(new_q, old_q)
            {
                let (hi, lo) = if ti > tj { (ti, tj) } else { (tj, ti) };
                tm.tets.swap_remove(hi);
                tm.tets.swap_remove(lo);
                for t in new_tets {
                    tm.tets.push(Tet {
                        nodes: t,
                        neighbors: [u32::MAX; 4],
                    });
                }
                tm.build_neighbors();
                swaps += 1;
                any = true;
            }
        }
        if !any {
            break;
        }
    }
    swaps
}

// ─── Pass 2b: 3-2 edge flip (edge removal) ────────────────────────────────────

/// 3-2 edge flip pass: every interior edge shared by exactly three tets
/// `(u,v,a)`-ring is removed, replacing the three tets with the pair
/// `(a,b,c,u)` + `(a,b,c,v)` over the opposite face, when the affected-set
/// quality tuple strictly improves (or passes the sliver-relaxed gate). This
/// completes the bistellar flip family named in the module docs — the 2-3
/// pass alone cannot remove edges. Returns the number of flips.
fn edge_removal_pass_3_2(tm: &mut TetMesh, metric: QualityMetric) -> usize {
    let mut flips = 0usize;
    for _sweep in 0..8 {
        // Candidate edges with exactly three incident tets (deterministic
        // order via the sorted edge map).
        let mut candidates: Vec<[u32; 2]> = Vec::new();
        for (key, ring) in tet_edge_map(tm) {
            if ring.len() == 3 {
                candidates.push(key);
            }
        }
        let mut any = false;
        for [u, v] in candidates {
            // Re-derive the ring in the current (churned) mesh; the three
            // opposite vertices must be distinct and pairwise connected by
            // the shared faces (i.e. they form the opposite triangle).
            let mut ring: Vec<[u32; 4]> = Vec::new();
            for t in &tm.tets {
                if t.nodes.contains(&u) && t.nodes.contains(&v) {
                    ring.push(t.nodes);
                }
            }
            if ring.len() != 3 {
                continue;
            }
            let mut opposite: Vec<u32> = Vec::new();
            for nodes in &ring {
                let opp: Vec<u32> = nodes
                    .iter()
                    .copied()
                    .filter(|&n| n != u && n != v)
                    .collect();
                if opp.len() != 2 {
                    opposite.clear();
                    break;
                }
                opposite.push(opp[0]);
                opposite.push(opp[1]);
            }
            if opposite.len() != 6 {
                continue;
            }
            opposite.sort_unstable();
            opposite.dedup();
            if opposite.len() != 3 {
                continue;
            }
            let [a, b, c] = [opposite[0], opposite[1], opposite[2]];
            if a == u || a == v || b == u || b == v || c == u || c == v {
                continue;
            }
            // Bistellar 3-2 validity: the removed edge must pierce the
            // interior of the opposite triangle, so that the bipyramid
            // (a,b,c,u,v) is decomposed exactly — the replacement preserves
            // the union of the three ring tets.
            let pu = tm.nodes[u as usize];
            let pv = tm.nodes[v as usize];
            let tri = [
                tm.nodes[a as usize],
                tm.nodes[b as usize],
                tm.nodes[c as usize],
            ];
            if !segment_pierces_triangle(pu, pv, tri) {
                continue;
            }
            let mut new_tets = [[a, b, c, u], [a, b, c, v]];
            let mut valid = true;
            for t in new_tets.iter_mut() {
                let pts = t.map(|n| tm.nodes[n as usize]);
                let vol = signed_volume_tet(pts[0], pts[1], pts[2], pts[3]);
                if vol < 0.0 {
                    t.swap(2, 3);
                } else if vol <= VOL_EPS {
                    valid = false;
                    break;
                }
            }
            if !valid {
                continue;
            }
            let old_set = ring;
            let Some(old_q) = affected_quality(tm, &old_set, metric) else {
                continue;
            };
            let Some(new_q) = affected_quality(tm, &new_tets, metric) else {
                continue;
            };
            if quality_tuple_better(new_q, old_q)
                || quality_tuple_sliver_relaxed(new_q, old_q)
            {
                // Remove the ring tets (descending for swap_remove safety),
                // push the pair, rebuild.
                let mut remove: Vec<u32> = Vec::new();
                for (ti, t) in tm.tets.iter().enumerate() {
                    if t.nodes.contains(&u) && t.nodes.contains(&v) {
                        remove.push(ti as u32);
                    }
                }
                remove.sort_unstable();
                remove.reverse();
                for ti in remove {
                    tm.tets.swap_remove(ti as usize);
                }
                for t in new_tets {
                    tm.tets.push(Tet {
                        nodes: t,
                        neighbors: [u32::MAX; 4],
                    });
                }
                tm.build_neighbors();
                flips += 1;
                any = true;
            }
        }
        if !any {
            break;
        }
    }
    flips
}

// ─── Pass 3: node relocation ──────────────────────────────────────────────────

/// One golden-section maximization of `f` over `[lo, hi]` (unimodal
/// assumption; 24 iterations, ratio 0.618 — the `optimizeMeshGRegion`
/// search scheme).
fn golden_max<F: Fn(f64) -> f64>(f: F, lo: f64, hi: f64) -> f64 {
    const PHI: f64 = 0.618_033_988_749_894_9;
    let (mut a, mut b) = (lo, hi);
    let mut c = b - PHI * (b - a);
    let mut d = a + PHI * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    for _ in 0..24 {
        if fc > fd {
            b = d;
            d = c;
            fd = fc;
            c = b - PHI * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + PHI * (b - a);
            fd = f(d);
        }
    }
    (a + b) * 0.5
}

/// Relocation pass: Gmsh `optimizeMeshGRegion`-style per-axis golden-section
/// line search over every free node's star. The objective is the star's
/// worst primary goodness; a found move is accepted only through the same
/// gates as every other pass — no star tet may invert, and the star's
/// lexicographic tuple (primary metric, radius-edge tie-break) must hold or
/// improve. Returns the number of nodes whose position changed.
fn relocate_pass(tm: &mut TetMesh, free: &[u32], metric: QualityMetric) -> usize {
    let star = star_tets(tm);
    let mut moved = 0usize;
    for &v in free {
        let tets = match star.get(&v) {
            Some(t) if !t.is_empty() => t,
            _ => continue,
        };
        let star_nodes: Vec<[u32; 4]> = tets.iter().map(|&ti| tm.tets[ti as usize].nodes).collect();
        let Some(base_q) = affected_quality(tm, &star_nodes, metric) else {
            continue;
        };
        // Local search bracket: half the mean star edge length (the star's
        // own scale, so the search cannot jump across neighbouring nodes).
        let pv = tm.nodes[v as usize];
        let mut h = 0.0f64;
        let mut hcount = 0usize;
        for &ti in tets {
            let nodes = tm.tets[ti as usize].nodes;
            for &n in nodes.iter() {
                if n != v {
                    h += len3(sub3(tm.nodes[n as usize], pv));
                    hcount += 1;
                }
            }
        }
        if hcount == 0 {
            continue;
        }
        h *= 0.5 / hcount as f64;
        if h <= 0.0 {
            continue;
        }

        // Star quality at a trial position: `None` when any star tet would
        // invert (signalled as the worst possible objective).
        let quality_at = |p: [f64; 3]| -> (f64, f64) {
            let mut primary = f64::MAX;
            let mut secondary = f64::MAX;
            for &ti in tets {
                let nodes = tm.tets[ti as usize].nodes;
                let pts: [[f64; 3]; 4] =
                    nodes.map(|n| if n == v { p } else { tm.nodes[n as usize] });
                if signed_volume_tet(pts[0], pts[1], pts[2], pts[3]) <= VOL_EPS {
                    return (f64::MIN, f64::MIN);
                }
                primary = primary.min(tet_goodness_points(
                    pts[0], pts[1], pts[2], pts[3], metric,
                ));
                secondary = secondary.min(match metric {
                    QualityMetric::RadiusEdgeRatio => {
                        min_dihedral_deg_tet(pts[0], pts[1], pts[2], pts[3])
                    }
                    _ => -radius_edge_ratio_tet(pts[0], pts[1], pts[2], pts[3]),
                });
            }
            (primary, secondary)
        };

        // Per-axis line search, starting from the current position. The
        // search objective is a bounded scalarization of the lexicographic
        // tuple — squashed primary + 1e-3 · squashed secondary — so the
        // search actively improves the radius-edge leg wherever the primary
        // is flat, while the primary always dominates (both squashes map to
        // (-1, 1); the weight gap cannot be bridged by the secondary leg).
        let squash = |x: f64| x / (1.0 + x.abs());
        let objective = |p: [f64; 3]| -> f64 {
            let (q0, q1) = quality_at(p);
            squash(q0) + 1.0e-3 * squash(q1)
        };
        let mut best = pv;
        let mut best_q = quality_at(best);
        for axis in 0..3 {
            let lo = best[axis] - h;
            let hi = best[axis] + h;
            let t = golden_max(
                |x| {
                    let mut p = best;
                    p[axis] = x;
                    objective(p)
                },
                lo,
                hi,
            );
            let mut p = best;
            p[axis] = t;
            let q = quality_at(p);
            if q.0 > best_q.0 || ((q.0 - best_q.0).abs() <= 1e-6 && q.1 < best_q.1) {
                best = p;
                best_q = q;
            }
        }

        if best != pv && quality_tuple_holds(best_q, base_q) {
            tm.nodes[v as usize] = best;
            moved += 1;
        }
    }
    moved
}

#[cfg(test)]
mod tests;
