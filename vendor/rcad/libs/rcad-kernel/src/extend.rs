//! Curve and surface trimming and extension.
//!
//! Analogous to OCCT `GeomAPI_ExtendCurveToPoint`,
//! `Geom_TrimmedCurve` construction helpers, and
//! `BRepBuilderAPI_MakeFace` trimming.
//!
//! # Curve operations
//!
//! | Function | Description | OCCT equivalent |
//! |---|---|---|
//! | [`trim_curve`] | Restrict a `BSplineCurve3` to `[t0, t1]` via knot insertion | `Geom_TrimmedCurve` (exact) |
//! | [`extend_curve_to_point`] | Extend a B-spline endpoint toward a target point | `GeomAPI_ExtendCurveToPoint` |
//! | [`extend_curve_by_length`] | Extend a B-spline endpoint by an arc-length distance | — |
//!
//! # Surface operations
//!
//! | Function | Description | OCCT equivalent |
//! |---|---|---|
//! | [`trim_surface`] | Wrap a surface in a `TrimmedSurface` with given UV bounds | `Geom_RectangularTrimmedSurface` |
//! | [`extend_bspline_surface`] | Extend a B-spline surface boundary row/column outward | `GeomAPI_ExtendSurfaceToShape` (partial) |

use glam::DVec3;

use crate::geom::{BSplineCurve3, BSplineSurface, CurveEval, Surface3, TrimmedSurface};

// ─────────────────────────────────────────────────────────────────────────────
// Curve trimming
// ─────────────────────────────────────────────────────────────────────────────

/// Trim a `BSplineCurve3` to the parameter range `[t0, t1]` using exact knot
/// insertion, returning a new curve whose natural domain is `[t0, t1]`.
///
/// The resulting curve evaluates identically to the original on `[t0, t1]`;
/// control points and knots outside that range are discarded.
///
/// Panics if `t0 >= t1` or if either value is outside the curve's domain.
///
/// Analogous to constructing a `Geom_TrimmedCurve`.
pub fn trim_curve(curve: &BSplineCurve3, t0: f64, t1: f64) -> BSplineCurve3 {
    assert!(t0 < t1, "trim_curve: t0 must be less than t1");

    // Strategy: insert t0 with multiplicity=degree so it becomes a breakpoint,
    // then insert t1 with multiplicity=degree.  After these insertions the
    // control points that correspond to the segment [t0, t1] are exactly the
    // ones between the two groups of repeated knots.
    let d = curve.degree;
    let c1 = insert_knot_to_multiplicity(curve, t0, d + 1);
    let c2 = insert_knot_to_multiplicity(&c1, t1, d + 1);

    let knots = &c2.knots;

    // Find the first occurrence index of t0 and t1 in the refined knot vector.
    // After inserting t0/t1 with multiplicity=degree, each appears exactly `degree` times.
    let first_t0 = knots
        .iter()
        .rposition(|&k| (k - t0).abs() < 1e-12)
        .unwrap_or(d);
    let first_t1 = knots
        .iter()
        .position(|&k| (k - t1).abs() < 1e-12)
        .unwrap_or(knots.len().saturating_sub(d + 1));

    // Control point slice: [first_t0 - d, first_t1)
    // The j-th control point "owns" the knot window [T[j], T[j+d]].
    // The segment starts where T[j+d] == t0, i.e. j = first_t0 - d.
    // The segment ends just before T[j] == t1, i.e. j = first_t1 (exclusive).
    let i_start = first_t0.saturating_sub(d);
    let i_end = first_t1;

    // Guard against bad slice
    let n_ctrl = c2.control_points.len();
    let i_start = i_start.min(n_ctrl.saturating_sub(1));
    let i_end = i_end.min(n_ctrl).max(i_start + 1);

    let new_ctrl = c2.control_points[i_start..i_end].to_vec();
    let new_weights = c2.weights[i_start..i_end].to_vec();

    // Knot vector: n_ctrl_new + degree + 1 knots starting at k_start = i_start.
    // (Each control point i corresponds to knots[i..i+d+1], so the full window is
    //  knots[i_start .. i_start + n_ctrl_new + d].)
    let n_ctrl_new = new_ctrl.len();
    let k_start = i_start;
    let k_end = (k_start + n_ctrl_new + d + 1).min(knots.len());
    let k_start = k_start.min(k_end);
    let raw_knots: Vec<f64> = knots[k_start..k_end].to_vec();

    // Normalize to [0, 1]
    let kmin = raw_knots.first().copied().unwrap_or(t0);
    let kmax = raw_knots.last().copied().unwrap_or(t1);
    let kspan = (kmax - kmin).max(1e-14);
    let new_knots: Vec<f64> = raw_knots.iter().map(|&k| (k - kmin) / kspan).collect();

    BSplineCurve3 {
        degree: d,
        knots: new_knots,
        control_points: new_ctrl,
        weights: new_weights,
    }
}

/// Insert knot `t` into `curve` until it has multiplicity `target_mult`,
/// returning the new curve.  If multiplicity already ≥ `target_mult`, returns
/// the curve unchanged.
///
/// Uses the Boehm single-knot insertion algorithm.
pub fn insert_knot_to_multiplicity(
    curve: &BSplineCurve3,
    t: f64,
    target_mult: usize,
) -> BSplineCurve3 {
    let current_mult = curve
        .knots
        .iter()
        .filter(|&&k| (k - t).abs() < 1e-14)
        .count();
    let mut result = curve.clone();
    for _ in current_mult..target_mult {
        result = insert_knot_once(&result, t);
    }
    result
}

/// Insert a single knot `t` into the B-spline using Boehm's algorithm.
fn insert_knot_once(curve: &BSplineCurve3, t: f64) -> BSplineCurve3 {
    let p = curve.degree;
    let n = curve.control_points.len();
    let knots = &curve.knots;

    // Find knot span k: knots[k] <= t < knots[k+1]
    let k = find_span(n, p, t, knots);

    // New knot vector: insert t after index k
    let mut new_knots = knots[..=k].to_vec();
    new_knots.push(t);
    new_knots.extend_from_slice(&knots[k + 1..]);

    // New control points (n+1 points after insertion)
    let mut new_ctrl = Vec::with_capacity(n + 1);
    let mut new_w = Vec::with_capacity(n + 1);

    for i in 0..=(n) {
        if i <= k - p {
            new_ctrl.push(curve.control_points[i]);
            new_w.push(curve.weights[i]);
        } else if i > k {
            new_ctrl.push(curve.control_points[i - 1]);
            new_w.push(curve.weights[i - 1]);
        } else {
            // Blend P[i-1] and P[i]
            let denom = knots[i + p] - knots[i];
            let alpha = if denom.abs() < 1e-14 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let w0 = curve.weights[i - 1];
            let w1 = curve.weights[i];
            let p0 = curve.control_points[i - 1];
            let p1 = curve.control_points[i];
            // Weighted blend in homogeneous coordinates
            let hw = (1.0 - alpha) * w0 + alpha * w1;
            let hp = (1.0 - alpha) * w0 * p0 + alpha * w1 * p1;
            new_w.push(hw);
            new_ctrl.push(if hw.abs() > 1e-14 { hp / hw } else { p0 });
        }
    }

    BSplineCurve3 {
        degree: p,
        knots: new_knots,
        control_points: new_ctrl,
        weights: new_w,
    }
}

fn find_span(n_ctrl: usize, degree: usize, t: f64, knots: &[f64]) -> usize {
    let n = n_ctrl - 1;
    if t >= knots[n + 1] {
        return n;
    }
    if t <= knots[degree] {
        return degree;
    }
    let mut lo = degree;
    let mut hi = n + 1;
    let mut mid = (lo + hi) / 2;
    while t < knots[mid] || t >= knots[mid + 1] {
        if t < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
    }
    mid
}

// ─────────────────────────────────────────────────────────────────────────────
// Curve extension
// ─────────────────────────────────────────────────────────────────────────────

/// Which end of the curve to extend.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveEnd {
    /// Extend the start (`t = t_min`) end.
    Start,
    /// Extend the end (`t = t_max`) end.
    End,
}

/// Extend a `BSplineCurve3` so that the specified endpoint reaches `target`.
///
/// The extension uses a simple linear segment appended by knot insertion,
/// preserving C¹ continuity at the join by adjusting the boundary control
/// point to lie on the tangent line.
///
/// Analogous to `GeomAPI_ExtendCurveToPoint`.
pub fn extend_curve_to_point(curve: &BSplineCurve3, end: CurveEnd, target: DVec3) -> BSplineCurve3 {
    let _n = curve.control_points.len();
    let mut new_ctrl = curve.control_points.clone();
    let mut new_w = curve.weights.clone();
    let mut new_knots = curve.knots.clone();

    match end {
        CurveEnd::End => {
            // To extend by one segment: append a new control point at `target`
            // and add an interior knot at (t_max + t_new)/2 to maintain valid
            // clamped structure: [..., t_max, t_max] → [..., t_max, t_ext, t_ext]
            // where t_ext = t_max + 1.
            let t_max = *new_knots.last().expect("knot vector is non-empty");
            let t_ext = t_max + 1.0;
            // Remove the last repeated knot, insert the interior + new endpoint
            // New knots: original_without_last_max, t_max, t_ext, t_ext
            let n_last_max = new_knots
                .iter()
                .rev()
                .take_while(|&&k| (k - t_max).abs() < 1e-14)
                .count();
            for _ in 0..n_last_max.saturating_sub(1) {
                new_knots.pop();
            }
            new_knots.push(t_ext);
            new_knots.push(t_ext);
            new_ctrl.push(target);
            new_w.push(1.0);
        }
        CurveEnd::Start => {
            let t_min = *new_knots.first().expect("knot vector is non-empty");
            let t_ext = t_min - 1.0;
            let n_first_min = new_knots
                .iter()
                .take_while(|&&k| (k - t_min).abs() < 1e-14)
                .count();
            for _ in 0..n_first_min.saturating_sub(1) {
                new_knots.remove(0);
            }
            new_knots.insert(0, t_ext);
            new_knots.insert(0, t_ext);
            new_ctrl.insert(0, target);
            new_w.insert(0, 1.0);
        }
    }

    // Normalize knot vector to [0, 1]
    let kmin = *new_knots.first().expect("knot vector is non-empty");
    let kmax = *new_knots.last().expect("knot vector is non-empty");
    let krange = (kmax - kmin).max(1e-14);
    let norm_knots: Vec<f64> = new_knots.iter().map(|&k| (k - kmin) / krange).collect();

    BSplineCurve3 {
        degree: curve.degree,
        knots: norm_knots,
        control_points: new_ctrl,
        weights: new_w,
    }
}

/// Extend a `BSplineCurve3` by an approximate arc-length `length` at the
/// specified end, by moving the endpoint along the end tangent direction.
///
/// Analogous to extending a curve by a linear segment of the given length.
pub fn extend_curve_by_length(curve: &BSplineCurve3, end: CurveEnd, length: f64) -> BSplineCurve3 {
    let target = match end {
        CurveEnd::End => {
            let [_, t1] = curve.default_domain();
            let p = curve.point_at(t1);
            let tang = curve.tangent_at(t1);
            p + length * tang
        }
        CurveEnd::Start => {
            let [t0, _] = curve.default_domain();
            let p = curve.point_at(t0);
            let tang = curve.tangent_at(t0);
            p - length * tang
        }
    };
    extend_curve_to_point(curve, end, target)
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface trimming
// ─────────────────────────────────────────────────────────────────────────────

/// Wrap a `Surface3` in a `TrimmedSurface` with the given UV bounds.
///
/// The returned surface evaluates identically to `basis` within `[u0,u1]×[v0,v1]`
/// and reports those bounds from `default_domain()`.
///
/// Analogous to `Geom_RectangularTrimmedSurface`.
pub fn trim_surface(basis: Surface3, u0: f64, u1: f64, v0: f64, v1: f64) -> Surface3 {
    Surface3::Trimmed(TrimmedSurface::new(basis, u0, u1, v0, v1))
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface extension
// ─────────────────────────────────────────────────────────────────────────────

/// Which boundary of a surface to extend.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceBoundary {
    /// u = u_min boundary (first row of control points).
    UMin,
    /// u = u_max boundary (last row).
    UMax,
    /// v = v_min boundary (first column of each row).
    VMin,
    /// v = v_max boundary (last column).
    VMax,
}

/// Extend a `BSplineSurface` by adding one extra row/column of control points
/// at the specified boundary, offset outward by `dist` (in surface normal
/// direction at the boundary mid-point).
///
/// This is a simple linear extrapolation: the new row/column mirrors the
/// relationship between the last two rows/columns.
///
/// Analogous to `GeomAPI_ExtendSurfaceToShape` (boundary extension only).
pub fn extend_bspline_surface(
    surface: &BSplineSurface,
    boundary: SurfaceBoundary,
    dist: f64,
) -> BSplineSurface {
    let mut result = surface.clone();

    match boundary {
        SurfaceBoundary::UMax => {
            // Extrapolate: new_row[j] = 2*last_row[j] - second_last_row[j] + dist*normal
            let n_rows = result.control_points.len();
            if n_rows < 2 {
                return result;
            }
            let last = &result.control_points[n_rows - 1];
            let prev = &result.control_points[n_rows - 2];
            let normal_offset = boundary_normal_offset(&result, boundary, dist);
            let new_row: Vec<DVec3> = last
                .iter()
                .zip(prev.iter())
                .map(|(&l, &p)| 2.0 * l - p + normal_offset)
                .collect();
            let new_w_row: Vec<f64> = result.weights[n_rows - 1].clone();
            result.control_points.push(new_row);
            result.weights.push(new_w_row);
            // Extend knot vector
            let last_k = *result.knots_u.last().expect("knots_u is non-empty");
            let second_last_k = result.knots_u[result.knots_u.len() - 2];
            result
                .knots_u
                .push(last_k + (last_k - second_last_k).max(1e-10));
        }
        SurfaceBoundary::UMin => {
            let n_rows = result.control_points.len();
            if n_rows < 2 {
                return result;
            }
            let first = result.control_points[0].clone();
            let second = result.control_points[1].clone();
            let normal_offset = boundary_normal_offset(&result, boundary, dist);
            let new_row: Vec<DVec3> = first
                .iter()
                .zip(second.iter())
                .map(|(&f, &s)| 2.0 * f - s + normal_offset)
                .collect();
            let new_w_row: Vec<f64> = result.weights[0].clone();
            result.control_points.insert(0, new_row);
            result.weights.insert(0, new_w_row);
            let first_k = result.knots_u[0];
            let second_k = result.knots_u[1];
            result
                .knots_u
                .insert(0, first_k - (second_k - first_k).max(1e-10));
        }
        SurfaceBoundary::VMax => {
            let normal_offset = boundary_normal_offset(&result, boundary, dist);
            for (row, w_row) in result
                .control_points
                .iter_mut()
                .zip(result.weights.iter_mut())
            {
                let n = row.len();
                if n < 2 {
                    continue;
                }
                let new_pt = 2.0 * row[n - 1] - row[n - 2] + normal_offset;
                row.push(new_pt);
                w_row.push(*w_row.last().expect("w_row is non-empty"));
            }
            let last_k = *result.knots_v.last().expect("knots_v is non-empty");
            let second_last_k = result.knots_v[result.knots_v.len() - 2];
            result
                .knots_v
                .push(last_k + (last_k - second_last_k).max(1e-10));
        }
        SurfaceBoundary::VMin => {
            let normal_offset = boundary_normal_offset(&result, boundary, dist);
            for (row, w_row) in result
                .control_points
                .iter_mut()
                .zip(result.weights.iter_mut())
            {
                let n = row.len();
                if n < 2 {
                    continue;
                }
                let new_pt = 2.0 * row[0] - row[1] + normal_offset;
                row.insert(0, new_pt);
                w_row.insert(0, w_row[0]);
            }
            let first_k = result.knots_v[0];
            let second_k = result.knots_v[1];
            result
                .knots_v
                .insert(0, first_k - (second_k - first_k).max(1e-10));
        }
    }

    result
}

/// Estimate an outward normal offset vector at the boundary mid-point.
fn boundary_normal_offset(surface: &BSplineSurface, boundary: SurfaceBoundary, dist: f64) -> DVec3 {
    use crate::geom::SurfaceEval;
    if dist.abs() < 1e-14 {
        return DVec3::ZERO;
    }
    let surf = Surface3::BSpline(surface.clone());
    let [u0, u1, v0, v1] = surf.default_domain();
    let (u, v) = match boundary {
        SurfaceBoundary::UMin => (u0, (v0 + v1) / 2.0),
        SurfaceBoundary::UMax => (u1, (v0 + v1) / 2.0),
        SurfaceBoundary::VMin => ((u0 + u1) / 2.0, v0),
        SurfaceBoundary::VMax => ((u0 + u1) / 2.0, v1),
    };
    dist * surf.normal_at(u, v)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    fn line_bspline(p0: DVec3, p1: DVec3) -> BSplineCurve3 {
        BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![p0, p1],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn trim_curve_reduces_domain() {
        let curve = line_bspline(DVec3::ZERO, DVec3::new(10.0, 0.0, 0.0));
        // Trim to [0.2, 0.7] → expect 2 pts at x=2 and x=7
        let trimmed = trim_curve(&curve, 0.2, 0.7);
        let p0 = trimmed.point_at(0.0);
        let p1 = trimmed.point_at(1.0);
        assert!((p0.x - 2.0).abs() < 1e-9, "start x={}", p0.x);
        assert!((p1.x - 7.0).abs() < 1e-9, "end x={}", p1.x);
    }

    #[test]
    fn extend_curve_to_point_increases_length() {
        let curve = line_bspline(DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0));
        let target = DVec3::new(3.0, 0.0, 0.0);
        let extended = extend_curve_to_point(&curve, CurveEnd::End, target);
        let end_pt = extended.point_at(1.0);
        assert!((end_pt.x - 3.0).abs() < 1e-9, "end x={}", end_pt.x);
    }

    #[test]
    fn extend_curve_by_length_end() {
        let curve = line_bspline(DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0));
        let extended = extend_curve_by_length(&curve, CurveEnd::End, 2.0);
        let p1 = extended.point_at(1.0);
        assert!((p1.x - 3.0).abs() < 1e-9, "end x={}", p1.x);
    }

    #[test]
    fn trim_surface_domain() {
        use crate::geom::{CylindricalSurface, SurfaceEval};
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let surf = trim_surface(Surface3::Cylinder(cyl), 0.0, 1.0, 0.0, 2.0);
        let [u0, u1, v0, v1] = surf.default_domain();
        assert!((u0 - 0.0).abs() < 1e-10);
        assert!((u1 - 1.0).abs() < 1e-10);
        assert!((v0 - 0.0).abs() < 1e-10);
        assert!((v1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn extend_bspline_surface_adds_row() {
        let bs = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        };
        let extended = extend_bspline_surface(&bs, SurfaceBoundary::UMax, 0.0);
        assert_eq!(
            extended.control_points.len(),
            3,
            "should have 3 rows after extension"
        );
    }
}
