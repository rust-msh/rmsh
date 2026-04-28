//! Curve fitting: construct a B-spline curve through a set of 3D points.
//!
//! Analogous to OCCT `GeomAPI_Interpolate` and `GeomAPI_PointsToBSpline`.
//!
//! Public functions:
//! - [`interpolate_points`] — exact interpolation (passes through all points, 3D)
//! - [`interpolate_points_2d`] — exact interpolation for 2D points
//! - [`approximate_points`] — least-squares B-spline approximation
//!
//! Both use **chord-length parameterization** and a **cubic (degree-3)** B-spline
//! with clamped knots.  The interpolation builds and solves a tridiagonal linear
//! system (Thomas algorithm) for the control points.

use glam::{DVec2, DVec3};

use crate::geom::{BSplineCurve2, BSplineCurve3};

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Error type for fitting operations.
#[derive(Debug, Clone, PartialEq)]
pub enum FitError {
    /// Fewer than 2 points were given.
    TooFewPoints,
    /// All points are coincident (zero total chord length).
    DegeneratePoints,
}

impl std::fmt::Display for FitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TooFewPoints => write!(f, "at least 2 points are required"),
            Self::DegeneratePoints => write!(f, "all points are coincident"),
        }
    }
}

impl std::error::Error for FitError {}

/// Build a cubic B-spline that passes **exactly** through each point in `pts`.
///
/// # Parameterization
/// Uses chord-length parameterization: `t[i] = sum of chord lengths up to i`,
/// normalized to `[0, 1]`.
///
/// # Knot vector
/// Clamped cubic knots are constructed from the chord-length parameters using
/// the averaging formula so that the system matrix is well-conditioned.
///
/// # Returns
/// A [`BSplineCurve3`] with degree 3, non-rational (all weights = 1).
///
/// # Errors
/// Returns [`FitError::TooFewPoints`] when fewer than 2 points are given.
/// Returns [`FitError::DegeneratePoints`] when all points coincide.
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::fit::interpolate_points;
/// let pts = vec![
///     DVec3::new(0.0, 0.0, 0.0),
///     DVec3::new(1.0, 1.0, 0.0),
///     DVec3::new(2.0, 0.0, 0.0),
/// ];
/// let curve = interpolate_points(&pts).unwrap();
/// ```
pub fn interpolate_points(pts: &[DVec3]) -> Result<BSplineCurve3, FitError> {
    let n = pts.len();
    if n < 2 {
        return Err(FitError::TooFewPoints);
    }

    // Chord-length parameters: t[0]=0, t[n-1]=1
    let params = chord_length_params(pts)?;

    let degree = 3_usize.min(n - 1); // degree ≤ n-1

    // Build clamped knot vector of size n + degree + 1
    let knots = clamped_knots_from_params(&params, degree);

    // Build the collocation matrix B (n×n banded):
    // B[i][j] = N_{j,degree}(t[i])
    // and solve B * ctrl = pts component-by-component using the Thomas algorithm
    // (valid because the matrix is strictly diagonally dominant for cubic splines).
    let ctrl = solve_interpolation(&params, &knots, degree, pts);

    Ok(BSplineCurve3 {
        degree,
        knots,
        control_points: ctrl,
        weights: vec![1.0; n],
    })
}

/// Build a cubic B-spline that passes **exactly** through each 2D point in `pts`.
///
/// This is the 2D analogue of [`interpolate_points`], producing a [`BSplineCurve2`].
///
/// # Parameterization
/// Uses chord-length parameterization normalized to `[0, 1]`.
///
/// # Knot vector
/// Clamped cubic knots constructed from the chord-length parameters.
///
/// # Returns
/// A [`BSplineCurve2`] with degree 3 (or `n-1` for small inputs), non-rational
/// (all weights = 1).
///
/// # Errors
/// Returns [`FitError::TooFewPoints`] when fewer than 2 points are given.
/// Returns [`FitError::DegeneratePoints`] when all points coincide.
///
/// # Examples
/// ```rust
/// use glam::DVec2;
/// use rcad_kernel::fit::interpolate_points_2d;
/// let pts = vec![
///     DVec2::new(0.0, 0.0),
///     DVec2::new(1.0, 1.0),
///     DVec2::new(2.0, 0.0),
/// ];
/// let curve = interpolate_points_2d(&pts).unwrap();
/// ```
pub fn interpolate_points_2d(pts: &[DVec2]) -> Result<BSplineCurve2, FitError> {
    let n = pts.len();
    if n < 2 {
        return Err(FitError::TooFewPoints);
    }

    // Chord-length parameters: t[0]=0, t[n-1]=1
    let params = chord_length_params_2d(pts)?;

    let degree = 3_usize.min(n - 1);

    // Build clamped knot vector
    let knots = clamped_knots_from_params(&params, degree);

    // Solve interpolation system via Gaussian elimination with partial pivoting
    let ctrl = solve_interpolation_2d(&params, &knots, degree, pts);

    Ok(BSplineCurve2 {
        degree,
        knots,
        control_points: ctrl,
        weights: vec![1.0; n],
    })
}

/// Build a cubic B-spline that **approximates** `pts` with `n_ctrl` control points
/// (must satisfy `2 ≤ n_ctrl < pts.len()`).
///
/// Uses least-squares fitting: minimises `Σ |B(t_i) - pts[i]|²` with a
/// uniformly-spaced knot vector.
///
/// The two end control points are pinned to the first and last data points,
/// so the curve always passes through the endpoints.
///
/// # Errors
/// Returns [`FitError::TooFewPoints`] when fewer than 2 points are given.
/// Returns [`FitError::DegeneratePoints`] when all points coincide.
pub fn approximate_points(pts: &[DVec3], n_ctrl: usize) -> Result<BSplineCurve3, FitError> {
    let n = pts.len();
    if n < 2 {
        return Err(FitError::TooFewPoints);
    }
    // Clamp n_ctrl to [2, n] — if n_ctrl >= n we fall back to exact interpolation
    let n_ctrl = n_ctrl.clamp(2, n);
    if n_ctrl == n {
        return interpolate_points(pts);
    }

    let params = chord_length_params(pts)?;
    let degree = 3_usize.min(n_ctrl - 1);
    let knots = uniform_clamped_knots(n_ctrl, degree);

    // Normal equations: AᵀA * ctrl = Aᵀ * pts  (per component)
    // Compute collocation matrix A (n × n_ctrl)
    let a = collocation_matrix(&params, &knots, degree, n, n_ctrl);

    // Pin endpoints
    // ctrl[0] = pts[0], ctrl[n_ctrl-1] = pts[n-1]
    // The inner n_ctrl-2 control points are the unknowns.
    let m = n_ctrl - 2; // number of free control points
    let ctrl_x = solve_least_squares(&a, pts, 0, pts[0].x, pts[n - 1].x, m, n_ctrl);
    let ctrl_y = solve_least_squares(&a, pts, 1, pts[0].y, pts[n - 1].y, m, n_ctrl);
    let ctrl_z = solve_least_squares(&a, pts, 2, pts[0].z, pts[n - 1].z, m, n_ctrl);

    let control_points: Vec<DVec3> = (0..n_ctrl)
        .map(|i| DVec3::new(ctrl_x[i], ctrl_y[i], ctrl_z[i]))
        .collect();

    Ok(BSplineCurve3 {
        degree,
        knots,
        control_points,
        weights: vec![1.0; n_ctrl],
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Parameter computation
// ─────────────────────────────────────────────────────────────────────────────

/// Chord-length parameterization normalized to [0, 1].
fn chord_length_params(pts: &[DVec3]) -> Result<Vec<f64>, FitError> {
    let n = pts.len();
    let mut params = Vec::with_capacity(n);
    params.push(0.0_f64);
    let mut total = 0.0_f64;
    for i in 1..n {
        total += (pts[i] - pts[i - 1]).length();
        params.push(total);
    }
    if total < 1e-14 {
        return Err(FitError::DegeneratePoints);
    }
    for p in &mut params {
        *p /= total;
    }
    Ok(params)
}

/// Chord-length parameterization for 2D points, normalized to [0, 1].
fn chord_length_params_2d(pts: &[DVec2]) -> Result<Vec<f64>, FitError> {
    let n = pts.len();
    let mut params = Vec::with_capacity(n);
    params.push(0.0_f64);
    let mut total = 0.0_f64;
    for i in 1..n {
        total += (pts[i] - pts[i - 1]).length();
        params.push(total);
    }
    if total < 1e-14 {
        return Err(FitError::DegeneratePoints);
    }
    for p in &mut params {
        *p /= total;
    }
    Ok(params)
}

// ─────────────────────────────────────────────────────────────────────────────
// Knot vector construction
// ─────────────────────────────────────────────────────────────────────────────

/// Clamped cubic knot vector derived from chord-length parameters.
///
/// For n data points (n control points in interpolation), the knot vector has
/// size `n + degree + 1`. Interior knots are averages of consecutive params
/// as per Piegl & Tiller §9.3.
fn clamped_knots_from_params(params: &[f64], degree: usize) -> Vec<f64> {
    let n = params.len(); // number of control points = number of data points
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

/// Uniform clamped knot vector for approximation with `n_ctrl` control points.
fn uniform_clamped_knots(n_ctrl: usize, degree: usize) -> Vec<f64> {
    let m = n_ctrl + degree + 1;
    let mut knots = vec![0.0_f64; m];
    let n_interior = n_ctrl - degree - 1;
    for knot in knots.iter_mut().take(degree + 1) {
        *knot = 0.0;
    }
    for knot in knots.iter_mut().skip(m - degree - 1) {
        *knot = 1.0;
    }
    if n_interior > 0 {
        for j in 1..=n_interior {
            knots[j + degree] = j as f64 / (n_interior + 1) as f64;
        }
    }
    knots
}

// ─────────────────────────────────────────────────────────────────────────────
// B-spline basis evaluation (Cox-de Boor recursion)
// ─────────────────────────────────────────────────────────────────────────────

/// Find the knot span index `i` such that `knots[i] <= t < knots[i+1]`.
/// Clamps to valid range.
fn find_span(n_ctrl: usize, degree: usize, t: f64, knots: &[f64]) -> usize {
    let n = n_ctrl - 1; // last control point index
    if t >= knots[n + 1] {
        return n;
    }
    if t <= knots[degree] {
        return degree;
    }
    // Binary search
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

/// Evaluate all non-zero basis functions N_{span-degree..=span, degree}(t).
/// Returns a slice of degree+1 values.
fn basis_fns(span: usize, t: f64, degree: usize, knots: &[f64]) -> Vec<f64> {
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

/// Evaluate all n_ctrl basis functions N_{i,degree}(t) (dense, n_ctrl values).
fn all_basis_fns(t: f64, knots: &[f64], degree: usize, n_ctrl: usize) -> Vec<f64> {
    let span = find_span(n_ctrl, degree, t, knots);
    let local = basis_fns(span, t, degree, knots);
    let mut result = vec![0.0_f64; n_ctrl];
    for (k, &val) in local.iter().enumerate().take(degree + 1) {
        let idx = span - degree + k;
        if idx < n_ctrl {
            result[idx] = val;
        }
    }
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Interpolation via Thomas algorithm (tridiagonal)
// ─────────────────────────────────────────────────────────────────────────────

/// Build full collocation matrix (n_data × n_ctrl).
fn collocation_matrix(
    params: &[f64],
    knots: &[f64],
    degree: usize,
    n_data: usize,
    n_ctrl: usize,
) -> Vec<Vec<f64>> {
    params[..n_data]
        .iter()
        .map(|&t| all_basis_fns(t, knots, degree, n_ctrl))
        .collect()
}

/// Solve the exact-interpolation system B * ctrl = data using the Thomas algorithm
/// (works because the cubic interpolation matrix is tridiagonal/banded with
/// strict diagonal dominance for chord-length parameterization).
///
/// Falls back to Gauss elimination for small n.
fn solve_interpolation(params: &[f64], knots: &[f64], degree: usize, pts: &[DVec3]) -> Vec<DVec3> {
    let n = pts.len();
    let a = collocation_matrix(params, knots, degree, n, n);

    // Solve three scalar tridiagonal systems (x, y, z) simultaneously
    // using forward-elimination + back-substitution.
    let rhs_x: Vec<f64> = pts.iter().map(|p| p.x).collect();
    let rhs_y: Vec<f64> = pts.iter().map(|p| p.y).collect();
    let rhs_z: Vec<f64> = pts.iter().map(|p| p.z).collect();

    let cx = gauss_solve(&a, &rhs_x);
    let cy = gauss_solve(&a, &rhs_y);
    let cz = gauss_solve(&a, &rhs_z);

    (0..n).map(|i| DVec3::new(cx[i], cy[i], cz[i])).collect()
}

/// Solve the 2D exact-interpolation system B * ctrl = data via Gaussian elimination.
fn solve_interpolation_2d(
    params: &[f64],
    knots: &[f64],
    degree: usize,
    pts: &[DVec2],
) -> Vec<DVec2> {
    let n = pts.len();
    let a = collocation_matrix(params, knots, degree, n, n);

    let rhs_x: Vec<f64> = pts.iter().map(|p| p.x).collect();
    let rhs_y: Vec<f64> = pts.iter().map(|p| p.y).collect();

    let cx = gauss_solve(&a, &rhs_x);
    let cy = gauss_solve(&a, &rhs_y);

    (0..n).map(|i| DVec2::new(cx[i], cy[i])).collect()
}

/// Gaussian elimination with partial pivoting for a dense n×n system.
/// Returns the solution vector.
fn gauss_solve(a: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    // Augmented matrix [A | b]
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
        // Partial pivot
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
            continue; // degenerate
        }

        for row in (col + 1)..n {
            let factor = mat[row][col] / pivot;
            let (lower, upper) = mat.split_at_mut(row);
            let pivot_row = &lower[col];
            let elim_row = &mut upper[0];
            for (elim_val, &pivot_val) in
                elim_row[col..=n].iter_mut().zip(pivot_row[col..=n].iter())
            {
                *elim_val -= pivot_val * factor;
            }
        }
    }

    // Back substitution
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

// ─────────────────────────────────────────────────────────────────────────────
// Least-squares approximation
// ─────────────────────────────────────────────────────────────────────────────

/// Solve the least-squares normal equations for one coordinate component.
/// Endpoints are pinned to `p0` (first) and `p1` (last).
/// Returns a Vec of length `n_ctrl`.
fn solve_least_squares(
    a: &[Vec<f64>], // collocation matrix n_data × n_ctrl
    pts: &[DVec3],
    coord: usize, // 0=x, 1=y, 2=z
    p0: f64,
    p1: f64,
    m: usize, // number of free unknowns = n_ctrl - 2
    n_ctrl: usize,
) -> Vec<f64> {
    let n_data = pts.len();

    // Rhs after pinning endpoints: r[i] = data[i] - A[i][0]*p0 - A[i][n_ctrl-1]*p1
    let rhs: Vec<f64> = (0..n_data)
        .map(|i| {
            let val = match coord {
                0 => pts[i].x,
                1 => pts[i].y,
                _ => pts[i].z,
            };
            val - a[i][0] * p0 - a[i][n_ctrl - 1] * p1
        })
        .collect();

    // Sub-matrix A_inner = A[:, 1..n_ctrl-1]  (n_data × m)
    let a_inner: Vec<Vec<f64>> = (0..n_data).map(|i| a[i][1..n_ctrl - 1].to_vec()).collect();

    // Normal equations: (A_inner^T A_inner) x = A_inner^T rhs
    let mut ata = vec![vec![0.0_f64; m]; m];
    let mut atr = vec![0.0_f64; m];
    for i in 0..n_data {
        for j in 0..m {
            atr[j] += a_inner[i][j] * rhs[i];
            for k in 0..m {
                ata[j][k] += a_inner[i][j] * a_inner[i][k];
            }
        }
    }

    let inner = gauss_solve(&ata, &atr);

    // Reconstruct full ctrl vector
    let mut ctrl = vec![0.0_f64; n_ctrl];
    ctrl[0] = p0;
    ctrl[1..(m + 1)].copy_from_slice(&inner[..m]);
    ctrl[n_ctrl - 1] = p1;
    ctrl
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::CurveEval;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn interpolate_two_points() {
        let pts = vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)];
        let curve = interpolate_points(&pts).unwrap();
        let p0 = curve.point_at(0.0);
        let p1 = curve.point_at(1.0);
        assert!(p0.distance(DVec3::ZERO) < 1e-6, "p0={p0}");
        assert!(p1.distance(DVec3::new(1.0, 0.0, 0.0)) < 1e-6, "p1={p1}");
    }

    #[test]
    fn interpolate_three_collinear() {
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
        ];
        let curve = interpolate_points(&pts).unwrap();
        for (i, &_pt) in pts.iter().enumerate() {
            let t = i as f64 / (pts.len() - 1) as f64;
            // The interpolated t isn't exactly equidistant in parameter space
            // but the curve passes through each point at its chord-length param.
            let _p = curve.point_at(t); // just verify no panic
        }
        // Endpoints must be exact
        let p0 = curve.point_at(0.0);
        let p1 = curve.point_at(1.0);
        assert!(p0.distance(pts[0]) < 1e-6);
        assert!(p1.distance(pts[2]) < 1e-6);
    }

    #[test]
    fn interpolate_arc_endpoints() {
        // Four points on a circular arc
        let pts: Vec<DVec3> = (0..4)
            .map(|i| {
                let a = std::f64::consts::FRAC_PI_2 * i as f64 / 3.0;
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect();
        let curve = interpolate_points(&pts).unwrap();
        // Endpoints
        assert!(curve.point_at(0.0).distance(pts[0]) < 1e-5);
        assert!(curve.point_at(1.0).distance(pts[3]) < 1e-5);
        // Interior points at their parameter values
        let params = chord_length_params(&pts).unwrap();
        for i in 1..3 {
            let p = curve.point_at(params[i]);
            assert!(
                p.distance(pts[i]) < 1e-4,
                "Point {i}: expected {:?} got {p} (err {})",
                pts[i],
                p.distance(pts[i])
            );
        }
    }

    #[test]
    fn interpolate_rejects_single_point() {
        assert!(matches!(
            interpolate_points(&[DVec3::ZERO]),
            Err(FitError::TooFewPoints)
        ));
    }

    #[test]
    fn interpolate_rejects_coincident_points() {
        let pts = vec![DVec3::ZERO; 3];
        assert!(matches!(
            interpolate_points(&pts),
            Err(FitError::DegeneratePoints)
        ));
    }

    #[test]
    fn approximate_falls_back_to_interpolate_when_n_ctrl_eq_n() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
        ];
        let c_approx = approximate_points(&pts, 3).unwrap();
        let c_interp = interpolate_points(&pts).unwrap();
        // Both should give same endpoints
        assert!(c_approx.point_at(0.0).distance(c_interp.point_at(0.0)) < 1e-6);
        assert!(c_approx.point_at(1.0).distance(c_interp.point_at(1.0)) < 1e-6);
    }

    #[test]
    fn approximate_points_basic() {
        // 10 points on a sine wave, approximated with 5 control points
        let pts: Vec<DVec3> = (0..10)
            .map(|i| {
                let x = i as f64 / 9.0;
                DVec3::new(x, (x * std::f64::consts::PI).sin(), 0.0)
            })
            .collect();
        let curve = approximate_points(&pts, 5).unwrap();
        // Endpoints pinned
        assert!(curve.point_at(0.0).distance(pts[0]) < 1e-6);
        assert!(curve.point_at(1.0).distance(pts[9]) < 1e-6);
        // Curve spans x ∈ [0,1]
        let mid = curve.point_at(0.5);
        assert!(mid.x > 0.1 && mid.x < 0.9, "midpoint x={}", mid.x);
    }

    #[test]
    fn chord_params_normalized() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
        ];
        let p = chord_length_params(&pts).unwrap();
        assert_eq!(p[0], 0.0);
        assert_eq!(p[2], 1.0);
        assert!(approx_eq(p[1], 1.0 / 3.0, 1e-12));
    }

    // ── 2D interpolation tests ──────────────────────────────────────────────

    #[test]
    fn interpolate_2d_line() {
        // 3 collinear points along y = x
        let pts = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(2.0, 2.0),
        ];
        use crate::geom::Curve2dEval;
        let curve = interpolate_points_2d(&pts).unwrap();

        // Endpoints must be exact
        let p0 = curve.point_at(0.0);
        let p1 = curve.point_at(1.0);
        assert!(p0.distance(DVec2::new(0.0, 0.0)) < 1e-6, "start={p0}");
        assert!(p1.distance(DVec2::new(2.0, 2.0)) < 1e-6, "end={p1}");

        // Midpoint should lie at (1, 1) since points are equally spaced
        let params = chord_length_params_2d(&pts).unwrap();
        let mid = curve.point_at(params[1]);
        assert!(mid.distance(DVec2::new(1.0, 1.0)) < 1e-5, "midpoint={mid}");
    }

    #[test]
    fn interpolate_2d_circle_arc() {
        // 9 points on a quarter circle of radius 1
        use crate::geom::Curve2dEval;
        let n = 9_usize;
        let pts: Vec<DVec2> = (0..n)
            .map(|i| {
                let a = std::f64::consts::FRAC_PI_2 * i as f64 / (n - 1) as f64;
                DVec2::new(a.cos(), a.sin())
            })
            .collect();

        let curve = interpolate_points_2d(&pts).unwrap();

        // Endpoints
        assert!(
            curve.point_at(0.0).distance(pts[0]) < 1e-5,
            "start err={}",
            curve.point_at(0.0).distance(pts[0])
        );
        assert!(
            curve.point_at(1.0).distance(*pts.last().unwrap()) < 1e-5,
            "end err={}",
            curve.point_at(1.0).distance(*pts.last().unwrap())
        );

        // Midpoint radius should be ≈ 1.0
        let mid_pt = curve.point_at(0.5);
        let radius = mid_pt.length();
        assert!((radius - 1.0).abs() < 0.01, "midpoint radius={radius}");
    }
}
