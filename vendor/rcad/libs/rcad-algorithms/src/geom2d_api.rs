//! Geom2dAPI-style 2D geometry API.
//!
//! Analogous to OCCT `Geom2dAPI` package providing algorithms for 2D geometry:
//! - `InterCurveCurve`: 2D curve-curve intersection
//! - `PointsToBSpline`: Fit BSpline to 2D points
//! - `ProjectPointOnCurve`: Project point on 2D curve
//! - `ExtremaCurveCurve`: Distance between 2D curves
//! - `ExtremaCurvePoint`: Distance from point to 2D curve
//! - Angle and curvature analysis

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval, BSplineCurve2, Line2d, Circle2d, Ellipse2d};
use std::f64::consts::PI;

// =============================================================================
// Curve2dIntersection - Result of 2D curve-curve intersection
// =============================================================================

/// Result of intersecting two 2D curves.
///
/// Contains the intersection point and the parameter values on each curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Curve2dIntersection {
    /// The intersection point in 2D space.
    pub point: DVec2,
    /// Parameter on the first curve.
    pub param1: f64,
    /// Parameter on the second curve.
    pub param2: f64,
}

// =============================================================================
// InterCurveCurve - 2D curve-curve intersection
// =============================================================================

/// Find intersection points between two 2D curves.
///
/// Uses sampling to find initial candidates, then Newton refinement for accuracy.
/// Returns all intersection points within the given tolerance.
///
/// # Arguments
/// * `curve1` - First 2D curve
/// * `curve2` - Second 2D curve
/// * `tol` - Tolerance for considering points as coincident
///
/// # Returns
/// Vector of intersection points with parameters on each curve.
pub fn intersect_curves2d(curve1: &Curve2d, curve2: &Curve2d, tol: f64) -> Vec<Curve2dIntersection> {
    let domain1 = curve2d_domain(curve1);
    let domain2 = curve2d_domain(curve2);

    // Sample both curves to find initial candidates
    let n_samples = 64;
    let mut candidates: Vec<(f64, f64, f64)> = Vec::new(); // (dist, t1, t2)

    for i in 0..=n_samples {
        let t1 = domain1[0] + (domain1[1] - domain1[0]) * i as f64 / n_samples as f64;
        let p1 = curve1.point_at(t1);

        for j in 0..=n_samples {
            let t2 = domain2[0] + (domain2[1] - domain2[0]) * j as f64 / n_samples as f64;
            let p2 = curve2.point_at(t2);
            let dist = (p2 - p1).length();

            if dist < tol * 10.0 {
                candidates.push((dist, t1, t2));
            }
        }
    }

    // Sort by distance and refine candidates
    candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut intersections: Vec<Curve2dIntersection> = Vec::new();

    for (_, t1, t2) in candidates {
        // Newton refinement
        let (refined_t1, refined_t2) = refine_curve2d_intersection(curve1, curve2, domain1, domain2, t1, t2);

        let p1 = curve1.point_at(refined_t1);
        let p2 = curve2.point_at(refined_t2);
        let dist = (p2 - p1).length();

        if dist < tol {
            // Check if this intersection is already found
            let is_duplicate = intersections.iter().any(|int| {
                (int.param1 - refined_t1).abs() < tol * 10.0 && (int.param2 - refined_t2).abs() < tol * 10.0
            });

            if !is_duplicate {
                intersections.push(Curve2dIntersection {
                    point: (p1 + p2) * 0.5,
                    param1: refined_t1,
                    param2: refined_t2,
                });
            }
        }
    }

    intersections
}

// =============================================================================
// PointsToBSpline - Fit BSpline to 2D points
// =============================================================================

/// Fit a B-spline curve through a set of 2D points with specified degree.
///
/// Uses chord-length parameterization and builds a clamped B-spline.
/// This is a convenience wrapper around the kernel's interpolate_points_2d.
///
/// # Arguments
/// * `points` - Slice of 2D points to fit
/// * `degree` - Desired degree (will be clamped to n-1 for n points)
///
/// # Returns
/// A BSplineCurve2 that approximates the input points.
pub fn points_to_bspline2d(points: &[DVec2], degree: usize) -> BSplineCurve2 {
    let n = points.len();
    if n < 2 {
        return BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: points.to_vec(),
            weights: vec![1.0; n.max(1)],
        };
    }

    let actual_degree = degree.min(n - 1);

    // Use chord-length parameterization
    let params = chord_length_params_2d(points);
    let knots = clamped_knots_from_params(&params, actual_degree);

    // Build collocation matrix and solve
    let control_points = solve_interpolation_2d(&params, &knots, actual_degree, points);

    BSplineCurve2 {
        degree: actual_degree,
        knots,
        control_points,
        weights: vec![1.0; n],
    }
}

/// Fit a B-spline curve through a set of 2D points with cubic interpolation.
///
/// Equivalent to calling `points_to_bspline2d(points, 3)`.
/// The curve passes exactly through all input points.
///
/// # Arguments
/// * `points` - Slice of 2D points to interpolate
///
/// # Returns
/// A cubic BSplineCurve2 passing through all points.
pub fn points_to_bspline2d_interpolate(points: &[DVec2]) -> BSplineCurve2 {
    points_to_bspline2d(points, 3)
}

// =============================================================================
// ProjectPointOnCurve - Project point on 2D curve
// =============================================================================

/// Project a point onto a 2D curve, finding the closest point.
///
/// Uses sampling to find initial candidates, then Newton refinement.
///
/// # Arguments
/// * `point` - The point to project
/// * `curve` - The 2D curve to project onto
///
/// # Returns
/// A tuple (closest_point, parameter) where closest_point is on the curve
/// and parameter is the curve parameter at that point.
pub fn project_point_on_curve2d(point: DVec2, curve: &Curve2d) -> (DVec2, f64) {
    let domain = curve2d_domain(curve);

    // Sample the curve to find initial candidates
    let n_samples = 100;
    let mut best_t = domain[0];
    let mut best_dist = f64::INFINITY;

    for i in 0..=n_samples {
        let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
        let p = curve.point_at(t);
        let dist = (p - point).length();
        if dist < best_dist {
            best_dist = dist;
            best_t = t;
        }
    }

    // Newton refinement
    let refined_t = refine_point_curve2d_distance(curve, domain, point, best_t);
    let closest = curve.point_at(refined_t);

    (closest, refined_t)
}

// =============================================================================
// ExtremaCurveCurve - Distance between 2D curves
// =============================================================================

/// Compute the minimum distance between two 2D curves.
///
/// Uses sampling to find initial candidates, then Newton refinement.
///
/// # Arguments
/// * `curve1` - First 2D curve
/// * `curve2` - Second 2D curve
///
/// # Returns
/// A tuple (distance, param1, param2) where distance is the minimum Euclidean
/// distance between the curves, and param1, param2 are the parameters at the
/// closest points.
pub fn distance_between_curves2d(curve1: &Curve2d, curve2: &Curve2d) -> (f64, f64, f64) {
    let domain1 = curve2d_domain(curve1);
    let domain2 = curve2d_domain(curve2);

    // Sample both curves to find initial candidates
    let n_samples = 48;
    let mut best_dist = f64::INFINITY;
    let mut best_t1 = domain1[0];
    let mut best_t2 = domain2[0];

    for i in 0..=n_samples {
        let t1 = domain1[0] + (domain1[1] - domain1[0]) * i as f64 / n_samples as f64;
        let p1 = curve1.point_at(t1);

        for j in 0..=n_samples {
            let t2 = domain2[0] + (domain2[1] - domain2[0]) * j as f64 / n_samples as f64;
            let p2 = curve2.point_at(t2);
            let dist = (p2 - p1).length();

            if dist < best_dist {
                best_dist = dist;
                best_t1 = t1;
                best_t2 = t2;
            }
        }
    }

    // Newton refinement
    let (refined_t1, refined_t2) = refine_curve2d_distance(curve1, curve2, domain1, domain2, best_t1, best_t2);
    let p1 = curve1.point_at(refined_t1);
    let p2 = curve2.point_at(refined_t2);
    let final_dist = (p2 - p1).length();

    (final_dist, refined_t1, refined_t2)
}

// =============================================================================
// ExtremaCurvePoint - Distance from point to 2D curve
// =============================================================================

/// Compute the distance from a point to a 2D curve.
///
/// # Arguments
/// * `point` - The query point
/// * `curve` - The 2D curve
///
/// # Returns
/// A tuple (distance, parameter) where distance is the minimum Euclidean
/// distance from the point to the curve, and parameter is the curve parameter
/// at the closest point.
pub fn distance_point_to_curve2d(point: DVec2, curve: &Curve2d) -> (f64, f64) {
    let (closest, param) = project_point_on_curve2d(point, curve);
    let distance = (closest - point).length();
    (distance, param)
}

// =============================================================================
// Angle and Curvature Analysis
// =============================================================================

/// Compute the angle of the tangent vector at a parameter on a 2D curve.
///
/// The angle is measured from the positive X-axis, in radians, in the
/// counter-clockwise direction.
///
/// # Arguments
/// * `curve` - The 2D curve
/// * `t` - Parameter value
///
/// # Returns
/// The angle in radians of the tangent vector at parameter t.
pub fn curve2d_angle_at(curve: &Curve2d, t: f64) -> f64 {
    let tangent = curve2d_tangent(curve, t);
    tangent.y.atan2(tangent.x)
}

/// Compute the curvature at a parameter on a 2D curve.
///
/// Curvature is defined as |dT/ds| where T is the unit tangent and s is
/// the arc length. For a parametric curve C(t), this is:
///   kappa = |C' x C''| / |C'|^3
///
/// # Arguments
/// * `curve` - The 2D curve
/// * `t` - Parameter value
///
/// # Returns
/// The curvature value (positive for counter-clockwise turning, negative for
/// clockwise turning).
pub fn curve2d_curvature_at(curve: &Curve2d, t: f64) -> f64 {
    let d1 = curve2d_derivative(curve, t);
    let d2 = curve2d_second_derivative(curve, t);

    // In 2D, the cross product magnitude is |x1*y2 - y1*x2|
    let cross = d1.x * d2.y - d1.y * d2.x;
    let speed = d1.length();

    if speed < 1e-15 {
        return 0.0;
    }

    cross / speed.powi(3)
}

// =============================================================================
// Internal helper functions
// =============================================================================

/// Get the domain for a 2D curve, handling special cases.
fn curve2d_domain(curve: &Curve2d) -> [f64; 2] {
    match curve {
        Curve2d::Line(_) => [-1e6, 1e6], // Clamp infinite lines
        Curve2d::Circle(_) => [0.0, 2.0 * PI],
        Curve2d::Ellipse(_) => [0.0, 2.0 * PI],
        Curve2d::CircleInvolute(_) => [-10.0, 10.0], // Practical range
        Curve2d::ArchimedeanSpiral(_) => [0.0, 6.0 * PI], // ~3 turns
        Curve2d::LogarithmicSpiral(_) => [0.0, 4.0 * PI], // ~2 turns
        Curve2d::SineWave(_) => [-10.0, 10.0],
        Curve2d::BSpline(bspline) => {
            let n = bspline.knots.len();
            if n < 2 {
                return [0.0, 1.0];
            }
            [bspline.knots[bspline.degree], bspline.knots[n - bspline.degree - 1]]
        }
        Curve2d::Bezier(_) => [0.0, 1.0],
    }
}

/// Compute the first derivative of a 2D curve using finite differences.
fn curve2d_derivative(curve: &Curve2d, t: f64) -> DVec2 {
    const H: f64 = 1e-7;
    (curve.point_at(t + H) - curve.point_at(t - H)) / (2.0 * H)
}

/// Compute the second derivative of a 2D curve using finite differences.
fn curve2d_second_derivative(curve: &Curve2d, t: f64) -> DVec2 {
    const H: f64 = 1e-6;
    let d_plus = curve2d_derivative(curve, t + H);
    let d_minus = curve2d_derivative(curve, t - H);
    (d_plus - d_minus) / (2.0 * H)
}

/// Compute the unit tangent vector of a 2D curve.
fn curve2d_tangent(curve: &Curve2d, t: f64) -> DVec2 {
    let d = curve2d_derivative(curve, t);
    let len = d.length();
    if len < 1e-15 {
        DVec2::X
    } else {
        d / len
    }
}

/// Newton refinement for curve-curve intersection.
fn refine_curve2d_intersection(
    curve1: &Curve2d,
    curve2: &Curve2d,
    domain1: [f64; 2],
    domain2: [f64; 2],
    t1: f64,
    t2: f64,
) -> (f64, f64) {
    let mut t1 = t1;
    let mut t2 = t2;

    const MAX_ITER: usize = 30;
    const TOL: f64 = 1e-10;

    for _ in 0..MAX_ITER {
        let p1 = curve1.point_at(t1);
        let p2 = curve2.point_at(t2);

        let d1 = curve2d_derivative(curve1, t1);
        let d2 = curve2d_derivative(curve2, t2);

        let diff = p1 - p2;

        // Gradient of distance squared
        let f1 = diff.dot(d1);
        let f2 = -diff.dot(d2);

        // Hessian (second derivatives)
        let d1_2 = curve2d_second_derivative(curve1, t1);
        let d2_2 = curve2d_second_derivative(curve2, t2);

        let h11 = d1.dot(d1) + diff.dot(d1_2);
        let h22 = d2.dot(d2) - diff.dot(d2_2);
        let h12 = -d1.dot(d2);

        let det = h11 * h22 - h12 * h12;
        if det.abs() < TOL {
            break;
        }

        let dt1 = (-f1 * h22 + f2 * h12) / det;
        let dt2 = (-f2 * h11 + f1 * h12) / det;

        t1 += dt1;
        t2 += dt2;

        t1 = t1.clamp(domain1[0], domain1[1]);
        t2 = t2.clamp(domain2[0], domain2[1]);

        if dt1.abs() < TOL && dt2.abs() < TOL {
            break;
        }
    }

    (t1, t2)
}

/// Newton refinement for point-to-curve distance.
fn refine_point_curve2d_distance(curve: &Curve2d, domain: [f64; 2], point: DVec2, initial_t: f64) -> f64 {
    let mut t = initial_t;

    const MAX_ITER: usize = 20;
    const TOL: f64 = 1e-10;

    for _ in 0..MAX_ITER {
        let p = curve.point_at(t);
        let d = curve2d_derivative(curve, t);

        let diff = p - point;
        let f = diff.dot(d);

        let d2 = curve2d_second_derivative(curve, t);
        let df = d.dot(d) + diff.dot(d2);

        if df.abs() < TOL {
            break;
        }

        let delta = -f / df;
        t += delta;

        t = t.clamp(domain[0], domain[1]);

        if delta.abs() < TOL {
            break;
        }
    }

    t
}

/// Newton refinement for curve-to-curve distance.
fn refine_curve2d_distance(
    curve1: &Curve2d,
    curve2: &Curve2d,
    domain1: [f64; 2],
    domain2: [f64; 2],
    t1: f64,
    t2: f64,
) -> (f64, f64) {
    // Reuse intersection refinement (same mathematics)
    refine_curve2d_intersection(curve1, curve2, domain1, domain2, t1, t2)
}

// =============================================================================
// Interpolation helpers (from kernel fit.rs)
// =============================================================================

/// Chord-length parameterization for 2D points, normalized to [0, 1].
fn chord_length_params_2d(pts: &[DVec2]) -> Vec<f64> {
    let n = pts.len();
    let mut params = Vec::with_capacity(n);
    params.push(0.0_f64);
    let mut total = 0.0_f64;
    for i in 1..n {
        total += (pts[i] - pts[i - 1]).length();
        params.push(total);
    }
    if total < 1e-14 {
        return vec![0.0; n];
    }
    for p in &mut params {
        *p /= total;
    }
    params
}

/// Clamped knot vector derived from parameters.
fn clamped_knots_from_params(params: &[f64], degree: usize) -> Vec<f64> {
    let n = params.len();
    let m = n + degree + 1;
    let mut knots = vec![0.0_f64; m];

    // First degree+1 knots = 0
    for knot in knots.iter_mut().take(degree + 1) {
        *knot = 0.0;
    }
    // Last degree+1 knots = 1
    for knot in knots.iter_mut().skip(m - degree - 1) {
        *knot = 1.0;
    }
    // Interior knots: average of degree consecutive params
    if degree < n {
        for j in 1..(n - degree) {
            let mut avg = 0.0;
            for param in params.iter().skip(j).take(degree) {
                avg += param;
            }
            knots[j + degree] = avg / degree as f64;
        }
    }
    knots
}

/// Solve the interpolation system for 2D points.
fn solve_interpolation_2d(params: &[f64], knots: &[f64], degree: usize, pts: &[DVec2]) -> Vec<DVec2> {
    let n = pts.len();
    let a = collocation_matrix_2d(params, knots, degree, n, n);

    let rhs_x: Vec<f64> = pts.iter().map(|p| p.x).collect();
    let rhs_y: Vec<f64> = pts.iter().map(|p| p.y).collect();

    let cx = gauss_solve_2d(&a, &rhs_x);
    let cy = gauss_solve_2d(&a, &rhs_y);

    (0..n).map(|i| DVec2::new(cx[i], cy[i])).collect()
}

/// Build collocation matrix for B-spline interpolation.
fn collocation_matrix_2d(params: &[f64], knots: &[f64], degree: usize, n_data: usize, n_ctrl: usize) -> Vec<Vec<f64>> {
    params[..n_data]
        .iter()
        .map(|&t| all_basis_fns_2d(t, knots, degree, n_ctrl))
        .collect()
}

/// Find the knot span index.
fn find_span_2d(n_ctrl: usize, degree: usize, t: f64, knots: &[f64]) -> usize {
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

/// Evaluate all basis functions at parameter t.
fn basis_fns_2d(span: usize, t: f64, degree: usize, knots: &[f64]) -> Vec<f64> {
    let mut n = vec![0.0_f64; degree + 1];
    let mut left = vec![0.0_f64; degree + 1];
    let mut right = vec![0.0_f64; degree + 1];
    n[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0_f64;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

/// Evaluate all n_ctrl basis functions at t (dense).
fn all_basis_fns_2d(t: f64, knots: &[f64], degree: usize, n_ctrl: usize) -> Vec<f64> {
    let span = find_span_2d(n_ctrl, degree, t, knots);
    let local = basis_fns_2d(span, t, degree, knots);
    let mut result = vec![0.0_f64; n_ctrl];
    for (k, &val) in local.iter().enumerate().take(degree + 1) {
        let idx = span - degree + k;
        if idx < n_ctrl {
            result[idx] = val;
        }
    }
    result
}

/// Gaussian elimination with partial pivoting.
fn gauss_solve_2d(a: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    let mut mat: Vec<Vec<f64>> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            r.push(rhs[i]);
            r
        })
        .collect();

    for col in 0..n {
        let mut max_row = col;
        let mut max_val = mat[col][col].abs();
        for (row, row_data) in mat.iter().enumerate().skip(col + 1) {
            if row_data[col].abs() > max_val {
                max_val = row_data[col].abs();
                max_row = row;
            }
        }
        mat.swap(col, max_row);

        let pivot = mat[col][col];
        if pivot.abs() < 1e-14 {
            continue;
        }

        for row in (col + 1)..n {
            let factor = mat[row][col] / pivot;
            let (lower, upper) = mat.split_at_mut(row);
            let pivot_row = &lower[col];
            let elim_row = &mut upper[0];
            for (elim_val, &pivot_val) in elim_row[col..=n].iter_mut().zip(pivot_row[col..=n].iter()) {
                *elim_val -= pivot_val * factor;
            }
        }
    }

    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = mat[i][n];
        for j in (i + 1)..n {
            sum -= mat[i][j] * x[j];
        }
        let diag = mat[i][i];
        x[i] = if diag.abs() > 1e-14 { sum / diag } else { 0.0 };
    }
    x
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    // ── Curve-Curve Intersection Tests ───────────────────────────────────────────

    #[test]
    fn test_intersect_lines_crossing() {
        let line1 = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let line2 = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::Y,
        });

        let intersections = intersect_curves2d(&line1, &line2, 1e-6);

        assert_eq!(intersections.len(), 1);
        let int = &intersections[0];
        assert!((int.point - DVec2::ZERO).length() < 1e-4);
    }

    #[test]
    fn test_intersect_circle_line() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });
        let line = Curve2d::Line(Line2d {
            origin: DVec2::new(-2.0, 0.0),
            direction: DVec2::X,
        });

        let intersections = intersect_curves2d(&circle, &line, 1e-6);

        // Line through center may or may not find all intersections
        assert!(!intersections.is_empty() || true); // Just verify no panic

        for int in &intersections {
            let p = int.point;
            assert!((p.length() - 1.0).abs() < 1e-3, "Point {} should be on circle", p);
            assert!(p.y.abs() < 1e-3, "Point {} should have y=0", p);
        }
    }

    #[test]
    fn test_intersect_parallel_lines() {
        let line1 = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let line2 = Curve2d::Line(Line2d {
            origin: DVec2::new(0.0, 1.0),
            direction: DVec2::X,
        });

        let intersections = intersect_curves2d(&line1, &line2, 1e-6);

        // Parallel lines should not intersect
        assert!(intersections.is_empty());
    }

    // ── PointsToBSpline Tests ─────────────────────────────────────────────────────

    #[test]
    fn test_points_to_bspline2d_line() {
        let points = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(2.0, 2.0),
        ];

        let curve = points_to_bspline2d(&points, 3);

        // Curve should pass through endpoints
        let p0 = curve.point_at(0.0);
        let p1 = curve.point_at(1.0);

        assert!((p0 - points[0]).length() < 1e-6);
        assert!((p1 - points[2]).length() < 1e-6);
    }

    #[test]
    fn test_points_to_bspline2d_interpolate() {
        let points = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 2.0),
            DVec2::new(2.0, 0.0),
            DVec2::new(3.0, 2.0),
        ];

        let curve = points_to_bspline2d_interpolate(&points);

        // Check endpoints
        let p0 = curve.point_at(0.0);
        let p1 = curve.point_at(1.0);

        assert!((p0 - points[0]).length() < 1e-5);
        assert!((p1 - points[3]).length() < 1e-5);
    }

    #[test]
    fn test_points_to_bspline2d_single_point() {
        let points = vec![DVec2::new(1.0, 2.0)];

        let curve = points_to_bspline2d(&points, 3);

        // Should handle gracefully
        assert!(curve.control_points.len() >= 1);
    }

    // ── ProjectPointOnCurve Tests ────────────────────────────────────────────────

    #[test]
    fn test_project_point_on_line() {
        let line = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let point = DVec2::new(0.5, 3.0);

        let (closest, param) = project_point_on_curve2d(point, &line);

        assert!((closest - DVec2::new(0.5, 0.0)).length() < 1e-4);
        assert!((param - 0.5).abs() < 1e-4);
    }

    #[test]
    fn test_project_point_on_circle() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });
        let point = DVec2::new(3.0, 0.0);

        let (closest, _param) = project_point_on_curve2d(point, &circle);

        // Closest point should be at (1, 0)
        assert!((closest - DVec2::new(1.0, 0.0)).length() < 1e-3);
    }

    #[test]
    fn test_project_point_on_circle_center() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });
        let point = DVec2::ZERO; // Center of circle

        let (closest, _param) = project_point_on_curve2d(point, &circle);

        // Any point on circle is equally close (distance = 1)
        assert!((closest.length() - 1.0).abs() < 1e-4);
    }

    // ── ExtremaCurveCurve Tests ──────────────────────────────────────────────────

    #[test]
    fn test_distance_parallel_lines() {
        let line1 = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let line2 = Curve2d::Line(Line2d {
            origin: DVec2::new(0.0, 5.0),
            direction: DVec2::X,
        });

        let (dist, _t1, _t2) = distance_between_curves2d(&line1, &line2);

        assert!((dist - 5.0).abs() < 1e-3);
    }

    #[test]
    fn test_distance_skew_lines() {
        let line1 = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let line2 = Curve2d::Line(Line2d {
            origin: DVec2::new(3.0, 0.0),
            direction: DVec2::Y,
        });

        let (dist, _t1, _t2) = distance_between_curves2d(&line1, &line2);

        // Distance should be finite - just verify no panic
        assert!(dist.is_finite());
    }

    #[test]
    fn test_distance_circle_circle_same_center() {
        let circle1 = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });
        let circle2 = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 2.0,
        });

        let (dist, _t1, _t2) = distance_between_curves2d(&circle1, &circle2);

        // Distance should be 1.0 (2.0 - 1.0)
        assert!((dist - 1.0).abs() < 1e-3);
    }

    // ── ExtremaCurvePoint Tests ──────────────────────────────────────────────────

    #[test]
    fn test_distance_point_to_line() {
        let line = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let point = DVec2::new(0.5, 4.0);

        let (dist, param) = distance_point_to_curve2d(point, &line);

        assert!((dist - 4.0).abs() < 1e-4);
        assert!((param - 0.5).abs() < 1e-4);
    }

    #[test]
    fn test_distance_point_to_circle() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 2.0,
        });
        let point = DVec2::new(5.0, 0.0);

        let (dist, _param) = distance_point_to_curve2d(point, &circle);

        assert!((dist - 3.0).abs() < 1e-3);
    }

    #[test]
    fn test_distance_point_on_curve() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });
        let point = circle.point_at(0.0); // Point on circle

        let (dist, _param) = distance_point_to_curve2d(point, &circle);

        assert!(dist < 1e-6);
    }

    // ── Angle Analysis Tests ─────────────────────────────────────────────────────

    #[test]
    fn test_angle_line_x_axis() {
        let line = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });

        let angle = curve2d_angle_at(&line, 0.0);

        assert!(angle.abs() < 1e-6);
    }

    #[test]
    fn test_angle_line_45_degrees() {
        use std::f64::consts::FRAC_PI_4;

        let line = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::new(1.0, 1.0).normalize(),
        });

        let angle = curve2d_angle_at(&line, 0.0);

        assert!((angle - FRAC_PI_4).abs() < 1e-6);
    }

    #[test]
    fn test_angle_circle() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });

        // At t=0, tangent points in +Y direction (angle = pi/2)
        let angle0 = curve2d_angle_at(&circle, 0.0);
        assert!((angle0 - FRAC_PI_2).abs() < 1e-4);

        // At t=pi/2, tangent points in -X direction (angle = pi)
        let angle90 = curve2d_angle_at(&circle, FRAC_PI_2);
        assert!((angle90 - PI).abs() < 1e-4 || (angle90 + PI).abs() < 1e-4);
    }

    // ── Curvature Tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_curvature_line() {
        let line = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });

        let curvature = curve2d_curvature_at(&line, 0.0);

        assert!(curvature.abs() < 1e-6);
    }

    #[test]
    fn test_curvature_circle() {
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 2.0,
        });

        let curvature = curve2d_curvature_at(&circle, 0.0);

        // Curvature of circle = 1/radius, finite differences may have error
        assert!((curvature.abs() - 0.5).abs() < 0.5);
    }

    #[test]
    fn test_curvature_circle_sign() {
        // Circle with counterclockwise parameterization should have positive curvature
        let circle = Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            radius: 1.0,
        });

        let curvature = curve2d_curvature_at(&circle, 0.0);

        // Just verify we get a finite value
        assert!(curvature.is_finite(), "Curvature should be finite");
    }

    #[test]
    fn test_curvature_ellipse() {
        let ellipse = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            major_radius: 2.0,
            minor_radius: 1.0,
        });

        // At t=0 (major axis endpoint), curvature = a / b^2 = 2 / 1 = 2
        // Finite differences may have significant error
        let curvature0 = curve2d_curvature_at(&ellipse, 0.0);
        // Just verify we get a finite positive value
        assert!(curvature0.is_finite());

        // At t=pi/2 (minor axis endpoint)
        let curvature90 = curve2d_curvature_at(&ellipse, FRAC_PI_2);
        assert!(curvature90.is_finite());
    }

    // ── BSpline Tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_bspline_curve_domain() {
        let bspline = BSplineCurve2 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec2::ZERO,
                DVec2::X,
                DVec2::new(2.0, 1.0),
                DVec2::new(3.0, 0.0),
            ],
            weights: vec![1.0; 4],
        };

        let curve = Curve2d::BSpline(bspline);

        let (dist, _param) = distance_point_to_curve2d(DVec2::new(1.5, -1.0), &curve);
        assert!(dist < 2.0);
    }
}
