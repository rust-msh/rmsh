//! Curve-curve extrema: find parameter pairs (s, t) minimising |C1(s) − C2(t)|.
//!
//! Analogous to OCCT `GeomAPI_ExtremaCurveCurve`.
//!
//! Algorithm:
//!   1. Coarse `n_samples × n_samples` grid scan over both curve domains.
//!   2. Newton-Raphson refinement from each grid minimum using finite-difference
//!      first and second derivatives.
//!   3. Deduplication within a parameter tolerance; sort by distance ascending.

use crate::geom::{Curve3, CurveEval};
use glam::DVec3;

/// A single local-minimum pair returned by [`extrema_curve_curve`].
#[derive(Debug, Clone)]
pub struct ExtremaPair {
    /// Parameter on the first curve at the closest approach.
    pub param1: f64,
    /// Parameter on the second curve at the closest approach.
    pub param2: f64,
    /// Point on the first curve.
    pub point1: DVec3,
    /// Point on the second curve.
    pub point2: DVec3,
    /// Euclidean distance at closest approach.
    pub distance: f64,
}

/// Collection of all local minima found between two curves.
#[derive(Debug, Clone)]
pub struct CurveCurveExtrema {
    /// All local minima, sorted by distance ascending.
    pub pairs: Vec<ExtremaPair>,
}

impl CurveCurveExtrema {
    /// Convenience: distance of the global (closest) minimum.
    /// Returns `f64::INFINITY` if no pairs were found.
    pub fn min_distance(&self) -> f64 {
        self.pairs
            .first()
            .map(|p| p.distance)
            .unwrap_or(f64::INFINITY)
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────────────────────

fn curve_domain(c: &Curve3) -> [f64; 2] {
    match c {
        Curve3::Line(l) => {
            let _ = l;
            [-1e6, 1e6] // infinite line: clamp to large range for sampling
        }
        other => other.default_domain(),
    }
}

/// Squared distance at parameter pair (s, t).
#[inline]
fn sq_dist(c1: &Curve3, c2: &Curve3, s: f64, t: f64) -> f64 {
    (c1.point_at(s) - c2.point_at(t)).length_squared()
}

const H: f64 = 1e-6; // finite-difference step

/// Gradient of f(s,t) = |C1(s) − C2(t)|²
///   df/ds = 2 · (C1(s)−C2(t)) · C1'(s)
///   df/dt = 2 · (C1(s)−C2(t)) · (−C2'(t))
fn gradient(c1: &Curve3, c2: &Curve3, s: f64, t: f64) -> [f64; 2] {
    let p1 = c1.point_at(s);
    let p2 = c2.point_at(t);
    let diff = p1 - p2;
    // finite-difference tangents (works for all Curve3 variants)
    let d1 = (c1.point_at(s + H) - c1.point_at(s - H)) / (2.0 * H);
    let d2 = (c2.point_at(t + H) - c2.point_at(t - H)) / (2.0 * H);
    [2.0 * diff.dot(d1), -2.0 * diff.dot(d2)]
}

/// Approximate Gauss-Newton Hessian diagonal (for robust step-size limiting).
fn hessian_diag(c1: &Curve3, c2: &Curve3, s: f64, t: f64) -> [f64; 2] {
    let d1 = (c1.point_at(s + H) - c1.point_at(s - H)) / (2.0 * H);
    let d2 = (c2.point_at(t + H) - c2.point_at(t - H)) / (2.0 * H);
    // Gauss-Newton diagonal: 2 * ||d_i||^2
    [
        2.0 * d1.length_squared().max(1e-30),
        2.0 * d2.length_squared().max(1e-30),
    ]
}

// ──────────────────────────────────────────────────────────────────────────────
// Newton refinement
// ──────────────────────────────────────────────────────────────────────────────

const MAX_ITER: usize = 50;
const GRAD_TOL: f64 = 1e-12;
const PARAM_TOL: f64 = 1e-9;

fn newton_refine(
    c1: &Curve3,
    c2: &Curve3,
    dom1: [f64; 2],
    dom2: [f64; 2],
    s0: f64,
    t0: f64,
) -> (f64, f64) {
    let mut s = s0;
    let mut t = t0;

    for _ in 0..MAX_ITER {
        let g = gradient(c1, c2, s, t);
        if g[0].abs() < GRAD_TOL && g[1].abs() < GRAD_TOL {
            break;
        }
        let h = hessian_diag(c1, c2, s, t);
        let ds = -g[0] / h[0];
        let dt = -g[1] / h[1];

        // Line-search: halve step until objective decreases
        let f0 = sq_dist(c1, c2, s, t);
        let mut alpha = 1.0;
        for _ in 0..8 {
            let ns = (s + alpha * ds).clamp(dom1[0], dom1[1]);
            let nt = (t + alpha * dt).clamp(dom2[0], dom2[1]);
            if sq_dist(c1, c2, ns, nt) < f0 {
                s = ns;
                t = nt;
                break;
            }
            alpha *= 0.5;
        }

        if (ds * alpha).abs() < PARAM_TOL && (dt * alpha).abs() < PARAM_TOL {
            break;
        }
    }
    (s, t)
}

// ──────────────────────────────────────────────────────────────────────────────
// Public API
// ──────────────────────────────────────────────────────────────────────────────

/// Find all local minima of the curve-curve distance function.
///
/// `n_samples` controls the coarse grid density used to seed the Newton
/// refinement. A value of 16–32 is sufficient for most analytic curves.
/// For complex B-splines use 64+.
pub fn extrema_curve_curve(c1: &Curve3, c2: &Curve3, n_samples: usize) -> CurveCurveExtrema {
    let dom1 = curve_domain(c1);
    let dom2 = curve_domain(c2);
    let n = n_samples.max(2);

    // ── 1. Coarse grid ────────────────────────────────────────────────────────
    // Evaluate sq_dist on an n×n grid; collect indices of local grid-minima.
    let ss: Vec<f64> = (0..n)
        .map(|i| dom1[0] + (dom1[1] - dom1[0]) * i as f64 / (n - 1) as f64)
        .collect();
    let tt: Vec<f64> = (0..n)
        .map(|j| dom2[0] + (dom2[1] - dom2[0]) * j as f64 / (n - 1) as f64)
        .collect();

    let mut grid = vec![vec![0.0f64; n]; n];
    for (i, &s) in ss.iter().enumerate() {
        for (j, &t) in tt.iter().enumerate() {
            grid[i][j] = sq_dist(c1, c2, s, t);
        }
    }

    // Collect local grid minima (cell smaller than all 8-connected neighbours).
    let mut seeds: Vec<(f64, f64)> = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let v = grid[i][j];
            let is_min = (0usize..2)
                .flat_map(|di| (0usize..2).map(move |dj| (di, dj)))
                .all(|(di, dj)| {
                    let ni = i.wrapping_add(di).wrapping_sub(1);
                    let nj = j.wrapping_add(dj).wrapping_sub(1);
                    if ni == i && nj == j {
                        return true;
                    }
                    if ni >= n || nj >= n {
                        return true;
                    }
                    grid[ni][nj] >= v
                });
            if is_min {
                seeds.push((ss[i], tt[j]));
            }
        }
    }
    // Always include corners and edges to catch boundary minima.
    for &s in &[dom1[0], dom1[1]] {
        for &t in &[dom2[0], dom2[1]] {
            seeds.push((s, t));
        }
    }

    // ── 2. Newton refinement ──────────────────────────────────────────────────
    let mut pairs: Vec<ExtremaPair> = seeds
        .iter()
        .map(|&(s0, t0)| {
            let (s, t) = newton_refine(c1, c2, dom1, dom2, s0, t0);
            let p1 = c1.point_at(s);
            let p2 = c2.point_at(t);
            ExtremaPair {
                param1: s,
                param2: t,
                point1: p1,
                point2: p2,
                distance: (p1 - p2).length(),
            }
        })
        .collect();

    // ── 3. Deduplicate within parameter tolerance ─────────────────────────────
    const DEDUP_TOL: f64 = 1e-4;
    pairs.sort_by(|a, b| {
        a.distance
            .partial_cmp(&b.distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut kept: Vec<ExtremaPair> = Vec::new();
    'outer: for p in pairs {
        for k in &kept {
            if (p.param1 - k.param1).abs() < DEDUP_TOL && (p.param2 - k.param2).abs() < DEDUP_TOL {
                continue 'outer;
            }
        }
        kept.push(p);
    }

    CurveCurveExtrema { pairs: kept }
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Circle3, Line3};
    use glam::DVec3;

    fn line(origin: DVec3, dir: DVec3) -> Curve3 {
        Curve3::Line(Line3 {
            origin,
            direction: dir.normalize(),
        })
    }

    fn circle(center: DVec3, normal: DVec3, radius: f64) -> Curve3 {
        Curve3::Circle(Circle3 {
            center,
            normal: normal.normalize(),
            radius,
        })
    }

    #[test]
    fn parallel_lines() {
        // Two parallel lines, d apart → distance = d
        let c1 = line(DVec3::ZERO, DVec3::X);
        let c2 = line(DVec3::new(0.0, 3.0, 0.0), DVec3::X);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        assert!(!ex.pairs.is_empty());
        let d = ex.min_distance();
        assert!((d - 3.0).abs() < 0.01, "expected 3.0, got {d}");
    }

    #[test]
    fn skew_lines() {
        // c1 along X, c2 along Y at height 5 → distance = 5
        let c1 = line(DVec3::ZERO, DVec3::X);
        let c2 = line(DVec3::new(0.0, 0.0, 5.0), DVec3::Y);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        let d = ex.min_distance();
        assert!((d - 5.0).abs() < 0.01, "expected 5.0, got {d}");
    }

    #[test]
    fn line_and_circle() {
        // Line along X at y=0, z=0; circle in XY plane r=2 centred at (10,0,0)
        // Closest: line point (8,0,0), circle point (8,0,0) — line passes through circle
        // Actually distance should be 0 when line intersects the circle region...
        // Use a line that doesn't pass through: along Z axis at (5,0,0) → circle in XY, r=2 centred origin
        let c1 = line(DVec3::new(5.0, 0.0, 0.0), DVec3::Z);
        let c2 = circle(DVec3::ZERO, DVec3::Z, 2.0);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        let d = ex.min_distance();
        // Line at x=5, circle r=2 at origin → closest circle point (2,0,0), line point (5,0,0) → d=3
        assert!((d - 3.0).abs() < 0.05, "expected ~3.0, got {d}");
    }

    #[test]
    fn concentric_circles() {
        // Two circles same center same normal, r=2 and r=5 → distance = 3
        let c1 = circle(DVec3::ZERO, DVec3::Z, 2.0);
        let c2 = circle(DVec3::ZERO, DVec3::Z, 5.0);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        let d = ex.min_distance();
        assert!((d - 3.0).abs() < 0.01, "expected 3.0, got {d}");
    }

    #[test]
    fn intersecting_lines_have_zero_distance() {
        // Two lines that cross at the origin → distance = 0
        let c1 = line(DVec3::ZERO, DVec3::X);
        let c2 = line(DVec3::ZERO, DVec3::Y);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        let d = ex.min_distance();
        assert!(d < 0.01, "crossing lines should have distance ≈ 0, got {d}");
    }

    #[test]
    fn same_circle_has_zero_min_distance() {
        // The same circle compared to itself → min distance = 0
        let c = circle(DVec3::ZERO, DVec3::Z, 3.0);
        let ex = extrema_curve_curve(&c, &c, 32);
        let d = ex.min_distance();
        assert!(d < 0.01, "same circle min distance should be 0, got {d}");
    }

    #[test]
    fn extrema_pairs_are_sorted_by_distance() {
        let c1 = line(DVec3::ZERO, DVec3::X);
        let c2 = circle(DVec3::ZERO, DVec3::Z, 3.0);
        let ex = extrema_curve_curve(&c1, &c2, 32);
        let distances: Vec<f64> = ex.pairs.iter().map(|p| p.distance).collect();
        for w in distances.windows(2) {
            assert!(w[0] <= w[1] + 1e-10, "pairs should be sorted ascending by distance");
        }
    }
}
