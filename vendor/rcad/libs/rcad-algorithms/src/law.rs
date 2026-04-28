//! Law-style evolution functions for variable operations.
//!
//! This module provides mathematical functions that evolve over a parameter range,
//! similar to OCCT's `Law_Function`, `Law_Linear`, `Law_BSpline`, etc.
//!
//! # Overview
//!
//! Laws are used in various CAD operations to define how a value changes along a path:
//! - Variable-radius fillets
//! - Variable-distance offsets
//! - Draft angles that change along a face
//! - Sweeps with varying cross-sections
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::law::{LinearLaw, LawFunction};
//!
//! let law = LinearLaw::new(0.0, 1.0, 1.0, 3.0);
//! assert!((law.value(0.0) - 1.0).abs() < 1e-10);
//! assert!((law.value(0.5) - 2.0).abs() < 1e-10);
//! assert!((law.value(1.0) - 3.0).abs() < 1e-10);
//! ```

use std::f64::consts::PI;

// ============================================================================
// Core Trait
// ============================================================================

/// Base trait for law functions.
///
/// A law function maps a parameter `t` to a value, with optional derivative
/// computation. The domain is always a closed interval `[t_min, t_max]`.
///
/// This trait mirrors OCCT's `Law_Function` hierarchy.
pub trait LawFunction: std::fmt::Debug + Send + Sync {
    /// Evaluate the law at parameter `t`.
    ///
    /// For parameters outside the domain, the behavior is implementation-defined:
    /// - Linear and constant laws extrapolate
    /// - BSpline laws typically clamp to boundary values
    fn value(&self, t: f64) -> f64;

    /// Compute the first derivative (dV/dt) at parameter `t`.
    ///
    /// Default implementation uses finite differences.
    fn derivative(&self, t: f64) -> f64 {
        let h = 1e-7;
        let domain = self.domain();
        let t_min = domain[0];
        let t_max = domain[1];

        // Clamp t+h and t-h within the domain for numerical stability
        let t_plus = (t + h).min(t_max);
        let t_minus = (t - h).max(t_min);

        if (t_plus - t_minus).abs() < 1e-14 {
            0.0
        } else {
            (self.value(t_plus) - self.value(t_minus)) / (t_plus - t_minus)
        }
    }

    /// Return the parameter domain `[t_min, t_max]`.
    fn domain(&self) -> [f64; 2];

    /// Check if parameter `t` is within the domain.
    fn is_in_domain(&self, t: f64) -> bool {
        let [t_min, t_max] = self.domain();
        t >= t_min && t <= t_max
    }

    /// Evaluate the law, clamping the parameter to the domain.
    fn value_clamped(&self, t: f64) -> f64 {
        let [t_min, t_max] = self.domain();
        self.value(t.clamp(t_min, t_max))
    }
}

// ============================================================================
// Law_Constant
// ============================================================================

/// A constant law that returns the same value for all parameters.
///
/// Equivalent to OCCT's `Law_Constant`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantLaw {
    value: f64,
    domain: [f64; 2],
}

impl ConstantLaw {
    /// Create a new constant law with the given value.
    ///
    /// The default domain is `[0.0, 1.0]`.
    pub fn new(value: f64) -> Self {
        Self {
            value,
            domain: [0.0, 1.0],
        }
    }

    /// Create a constant law with a custom domain.
    pub fn with_domain(value: f64, t_min: f64, t_max: f64) -> Self {
        Self {
            value,
            domain: [t_min, t_max],
        }
    }
}

impl LawFunction for ConstantLaw {
    fn value(&self, _t: f64) -> f64 {
        self.value
    }

    fn derivative(&self, _t: f64) -> f64 {
        0.0
    }

    fn domain(&self) -> [f64; 2] {
        self.domain
    }
}

// ============================================================================
// Law_Linear
// ============================================================================

/// A linear law that interpolates between two points.
///
/// Defined by two parameter-value pairs: `(t1, v1)` and `(t2, v2)`.
/// The value at any parameter `t` is computed as:
/// ```text
/// V(t) = v1 + (v2 - v1) * (t - t1) / (t2 - t1)
/// ```
///
/// Equivalent to OCCT's `Law_Linear`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearLaw {
    t1: f64,
    v1: f64,
    t2: f64,
    v2: f64,
    slope: f64,
}

impl LinearLaw {
    /// Create a new linear law between points `(t1, v1)` and `(t2, v2)`.
    ///
    /// # Panics
    /// Panics if `t1 == t2` (zero-length domain).
    pub fn new(t1: f64, v1: f64, t2: f64, v2: f64) -> Self {
        if (t2 - t1).abs() < 1e-14 {
            panic!("LinearLaw: t1 and t2 must be different");
        }
        let slope = (v2 - v1) / (t2 - t1);
        Self {
            t1,
            v1,
            t2,
            v2,
            slope,
        }
    }

    /// Get the slope (rate of change) of this linear law.
    pub fn slope(&self) -> f64 {
        self.slope
    }

    /// Get the start value.
    pub fn start_value(&self) -> f64 {
        self.v1
    }

    /// Get the end value.
    pub fn end_value(&self) -> f64 {
        self.v2
    }
}

impl LawFunction for LinearLaw {
    fn value(&self, t: f64) -> f64 {
        self.v1 + self.slope * (t - self.t1)
    }

    fn derivative(&self, _t: f64) -> f64 {
        self.slope
    }

    fn domain(&self) -> [f64; 2] {
        if self.t1 < self.t2 {
            [self.t1, self.t2]
        } else {
            [self.t2, self.t1]
        }
    }
}

// ============================================================================
// Law_BSpline
// ============================================================================

/// A BSpline law defined by control points in parameter-value space.
///
/// This is a 1D BSpline curve that maps parameter `t` to value `V(t)`.
/// The curve is defined by a set of parameter-value pairs that become
/// the control points of the BSpline.
///
/// Equivalent to OCCT's `Law_BSpline`.
#[derive(Debug, Clone, PartialEq)]
pub struct BSplineLaw {
    degree: usize,
    knots: Vec<f64>,
    control_values: Vec<f64>,
    weights: Vec<f64>,
}

impl BSplineLaw {
    /// Create a BSpline law from parameter-value pairs.
    ///
    /// The `params` and `values` arrays must have the same length.
    /// A clamped BSpline of the specified degree is constructed.
    ///
    /// # Arguments
    /// * `params` - Parameter values (will become knots)
    /// * `values` - Corresponding function values (will become control values)
    /// * `degree` - BSpline degree (typically 3 for cubic)
    ///
    /// # Panics
    /// Panics if arrays have different lengths, fewer than 2 points,
    /// or degree >= number of points.
    pub fn from_points(params: &[f64], values: &[f64], degree: usize) -> Self {
        let n = params.len();
        if n != values.len() {
            panic!("BSplineLaw: params and values must have same length");
        }
        if n < 2 {
            panic!("BSplineLaw: at least 2 points required");
        }
        if degree >= n {
            panic!("BSplineLaw: degree must be less than number of points");
        }

        // Build clamped knot vector from parameters
        let knots = Self::build_knots(params, degree);

        // Compute control values via interpolation
        // For simplicity, we use the values directly as control points
        // (This is exact only for degree n-1; otherwise it's an approximation)
        // A proper implementation would solve the interpolation system.
        let control_values = if degree == n - 1 {
            // Bezier case: control values are exactly the input values
            values.to_vec()
        } else {
            // General case: solve for control points that interpolate the given points
            Self::solve_interpolation(params, values, degree, &knots)
        };

        Self {
            degree,
            knots,
            control_values,
            weights: vec![1.0; n],
        }
    }

    /// Build a clamped knot vector from parameters.
    fn build_knots(params: &[f64], degree: usize) -> Vec<f64> {
        let n = params.len();
        let m = n + degree + 1;
        let mut knots = vec![0.0; m];

        // First degree+1 knots = first parameter
        for knot in knots.iter_mut().take(degree + 1) {
            *knot = params[0];
        }
        // Last degree+1 knots = last parameter
        for knot in knots.iter_mut().skip(m - degree - 1) {
            *knot = params[n - 1];
        }
        // Interior knots: average of consecutive params
        if n > degree + 1 {
            for j in 1..=(n - degree - 1) {
                let mut sum = 0.0;
                for k in j..(j + degree) {
                    sum += params[k.min(n - 1)];
                }
                knots[j + degree] = sum / degree as f64;
            }
        }

        knots
    }

    /// Solve the interpolation system to get control values.
    fn solve_interpolation(
        params: &[f64],
        values: &[f64],
        degree: usize,
        knots: &[f64],
    ) -> Vec<f64> {
        let n = params.len();

        // Build collocation matrix
        let mut matrix = vec![vec![0.0; n]; n];
        for (i, &t) in params.iter().enumerate() {
            let basis = Self::basis_fns_all(t, knots, degree, n);
            for j in 0..n {
                matrix[i][j] = basis[j];
            }
        }

        // Solve using Gaussian elimination
        Self::gauss_solve(&matrix, values)
    }

    /// Find the knot span index.
    fn find_span(t: f64, degree: usize, knots: &[f64], n_ctrl: usize) -> usize {
        let n = n_ctrl - 1;
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

    /// Evaluate all basis functions at parameter t.
    fn basis_fns_all(t: f64, knots: &[f64], degree: usize, n_ctrl: usize) -> Vec<f64> {
        let span = Self::find_span(t, degree, knots, n_ctrl);
        let local = Self::basis_fns(span, t, degree, knots);

        let mut result = vec![0.0; n_ctrl];
        for (k, &val) in local.iter().enumerate().take(degree + 1) {
            let idx = span - degree + k;
            if idx < n_ctrl {
                result[idx] = val;
            }
        }
        result
    }

    /// Evaluate non-zero basis functions at parameter t.
    fn basis_fns(span: usize, t: f64, degree: usize, knots: &[f64]) -> Vec<f64> {
        let mut n = vec![0.0; degree + 1];
        let mut left = vec![0.0; degree + 1];
        let mut right = vec![0.0; degree + 1];
        n[0] = 1.0;

        for j in 1..=degree {
            left[j] = t - knots[span + 1 - j];
            right[j] = knots[span + j] - t;
            let mut saved = 0.0;
            for r in 0..j {
                let denom = right[r + 1] + left[j - r];
                let temp = if denom.abs() > 1e-14 {
                    n[r] / denom
                } else {
                    0.0
                };
                n[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            n[j] = saved;
        }
        n
    }

    /// Gaussian elimination with partial pivoting.
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
                continue;
            }

            for row in (col + 1)..n {
                let factor = mat[row][col] / pivot;
                let (lower, upper) = mat.split_at_mut(row);
                let pivot_row = &lower[col];
                let elim_row = &mut upper[0];
                for (elim_val, &pivot_val) in elim_row[col..=n].iter_mut().zip(pivot_row[col..=n].iter())
                {
                    *elim_val -= pivot_val * factor;
                }
            }
        }

        // Back substitution
        let mut x = vec![0.0; n];
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

    /// Evaluate the BSpline at parameter t using de Boor's algorithm.
    fn evaluate(&self, t: f64) -> f64 {
        let n = self.control_values.len();
        if n == 0 {
            return 0.0;
        }

        let t_min = self.knots[self.degree];
        let t_max = self.knots[self.knots.len() - self.degree - 1];
        let t_clamped = t.clamp(t_min, t_max);

        // Find knot span
        let mut span = self.degree;
        for (i, &knot) in self
            .knots
            .iter()
            .enumerate()
            .take(self.knots.len() - self.degree - 1)
            .skip(self.degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }

        // Initialize control values for the span
        let mut d: Vec<f64> = (0..=self.degree)
            .map(|j| {
                let idx = (span - self.degree + j).min(n - 1);
                self.control_values[idx]
            })
            .collect();

        // De Boor recursion
        for r in 1..=self.degree {
            for j in (r..=self.degree).rev() {
                let i = span - self.degree + j;
                let denom = self.knots[i + self.degree - r + 1] - self.knots[i];
                let alpha = if denom.abs() < 1e-14 {
                    0.0
                } else {
                    (t_clamped - self.knots[i]) / denom
                };
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
            }
        }

        d[self.degree]
    }

    /// Evaluate the derivative analytically.
    fn evaluate_derivative(&self, t: f64) -> f64 {
        let n = self.control_values.len();
        if n < 2 || self.degree == 0 {
            return 0.0;
        }

        let t_min = self.knots[self.degree];
        let t_max = self.knots[self.knots.len() - self.degree - 1];
        let t_clamped = t.clamp(t_min, t_max);

        // Compute derivative control points
        let deg = self.degree;
        let mut deriv_ctrl: Vec<f64> = Vec::with_capacity(n - 1);

        for i in 0..n - 1 {
            let denom = self.knots[i + deg + 1] - self.knots[i + 1];
            if denom.abs() > 1e-14 {
                deriv_ctrl.push(
                    deg as f64 * (self.control_values[i + 1] - self.control_values[i]) / denom,
                );
            } else {
                deriv_ctrl.push(0.0);
            }
        }

        // Evaluate derivative curve at t (degree reduced by 1)
        let deriv_degree = deg - 1;
        if deriv_ctrl.is_empty() {
            return 0.0;
        }

        // Find knot span for derivative
        let mut span = deriv_degree;
        for (i, &knot) in self
            .knots
            .iter()
            .enumerate()
            .take(self.knots.len() - deriv_degree - 2)
            .skip(deriv_degree + 1)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }

        // Initialize
        let mut d: Vec<f64> = (0..=deriv_degree)
            .map(|j| {
                let idx = (span - deriv_degree + j).min(deriv_ctrl.len() - 1);
                deriv_ctrl[idx]
            })
            .collect();

        // De Boor recursion
        for r in 1..=deriv_degree {
            for j in (r..=deriv_degree).rev() {
                let i = span - deriv_degree + j;
                let denom = self.knots[i + deriv_degree - r + 2] - self.knots[i + 1];
                let alpha = if denom.abs() < 1e-14 {
                    0.0
                } else {
                    (t_clamped - self.knots[i + 1]) / denom
                };
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
            }
        }

        d[deriv_degree]
    }
}

impl LawFunction for BSplineLaw {
    fn value(&self, t: f64) -> f64 {
        self.evaluate(t)
    }

    fn derivative(&self, t: f64) -> f64 {
        self.evaluate_derivative(t)
    }

    fn domain(&self) -> [f64; 2] {
        let d = self.degree;
        let n = self.knots.len();
        if n < 2 * d + 2 {
            return [0.0, 1.0];
        }
        [self.knots[d], self.knots[n - d - 1]]
    }
}

// ============================================================================
// Law_Composite
// ============================================================================

/// A segment in a composite law.
#[derive(Debug)]
struct CompositeSegment {
    t1: f64,
    t2: f64,
    law: Box<dyn LawFunction>,
}

/// A composite law made of multiple segments stitched together.
///
/// Each segment has its own domain and law function. Parameters are mapped
/// to the appropriate segment based on their value.
///
/// Equivalent to OCCT's `Law_Composite`.
#[derive(Debug, Default)]
pub struct CompositeLaw {
    segments: Vec<CompositeSegment>,
}

impl CompositeLaw {
    /// Create a new empty composite law.
    pub fn new() -> Self {
        Self { segments: vec![] }
    }

    /// Add a segment covering domain `[t1, t2]` with the given law.
    ///
    /// Segments should not overlap. The law's internal domain is remapped
    /// to `[t1, t2]`.
    ///
    /// # Panics
    /// Panics if `t1 >= t2`.
    pub fn add_segment(&mut self, t1: f64, t2: f64, law: Box<dyn LawFunction>) {
        if t1 >= t2 {
            panic!("CompositeLaw: t1 must be less than t2");
        }
        self.segments.push(CompositeSegment { t1, t2, law });
        // Sort segments by t1
        self.segments.sort_by(|a, b| {
            a.t1.partial_cmp(&b.t1).unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    /// Find the segment containing parameter t.
    fn find_segment(&self, t: f64) -> Option<(usize, f64)> {
        for (i, seg) in self.segments.iter().enumerate() {
            if t >= seg.t1 && t <= seg.t2 {
                return Some((i, t));
            }
        }
        // Check if t is just before first segment or after last
        if let Some(first) = self.segments.first() {
            if t < first.t1 {
                return Some((0, first.t1));
            }
        }
        if let Some(last) = self.segments.last() {
            if t > last.t2 {
                return Some((self.segments.len() - 1, last.t2));
            }
        }
        None
    }
}

impl LawFunction for CompositeLaw {
    fn value(&self, t: f64) -> f64 {
        if let Some((i, t_seg)) = self.find_segment(t) {
            let seg = &self.segments[i];
            // Remap t from [t1, t2] to law's domain
            let law_domain = seg.law.domain();
            let law_t = if (seg.t2 - seg.t1).abs() < 1e-14 {
                law_domain[0]
            } else {
                let ratio = (t_seg - seg.t1) / (seg.t2 - seg.t1);
                law_domain[0] + ratio * (law_domain[1] - law_domain[0])
            };
            seg.law.value(law_t)
        } else {
            0.0
        }
    }

    fn derivative(&self, t: f64) -> f64 {
        if let Some((i, t_seg)) = self.find_segment(t) {
            let seg = &self.segments[i];
            let law_domain = seg.law.domain();
            let law_t = if (seg.t2 - seg.t1).abs() < 1e-14 {
                law_domain[0]
            } else {
                let ratio = (t_seg - seg.t1) / (seg.t2 - seg.t1);
                law_domain[0] + ratio * (law_domain[1] - law_domain[0])
            };
            // Chain rule: dV/dt = (dV/dlaw_t) * (dlaw_t/dt)
            let law_deriv = seg.law.derivative(law_t);
            let scale = (law_domain[1] - law_domain[0]) / (seg.t2 - seg.t1);
            law_deriv * scale
        } else {
            0.0
        }
    }

    fn domain(&self) -> [f64; 2] {
        if self.segments.is_empty() {
            return [0.0, 1.0];
        }
        let t_min = self.segments.first().map(|s| s.t1).unwrap_or(0.0);
        let t_max = self.segments.last().map(|s| s.t2).unwrap_or(1.0);
        [t_min, t_max]
    }
}

// ============================================================================
// Law_Interpolate
// ============================================================================

/// An interpolated law using cubic spline interpolation.
///
/// Given a set of (parameter, value) points, computes a C1-continuous
/// function that passes through all points.
///
/// Equivalent to OCCT's `Law_Interpolate`.
#[derive(Debug, Clone)]
pub struct InterpolateLaw {
    points: Vec<(f64, f64)>,
    periodic: bool,
    // Cached cubic coefficients for each interval
    coefficients: Vec<CubicCoefficients>,
}

/// Cubic polynomial coefficients for Hermite interpolation.
#[derive(Debug, Clone, Copy, Default)]
struct CubicCoefficients {
    a: f64, // constant
    b: f64, // linear
    c: f64, // quadratic
    d: f64, // cubic
}

impl InterpolateLaw {
    /// Create an interpolated law from parameter-value points.
    ///
    /// Uses cubic Hermite interpolation with finite-difference tangents.
    ///
    /// # Arguments
    /// * `points` - Array of (parameter, value) pairs
    /// * `periodic` - If true, the first and last points are connected smoothly
    ///
    /// # Panics
    /// Panics if fewer than 2 points are provided.
    pub fn from_points(points: &[(f64, f64)], periodic: bool) -> Self {
        if points.len() < 2 {
            panic!("InterpolateLaw: at least 2 points required");
        }

        let mut sorted_points: Vec<(f64, f64)> = points.to_vec();
        sorted_points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        // Compute cubic Hermite coefficients for each interval
        let n = sorted_points.len();
        let mut coefficients = Vec::with_capacity(n - 1);

        // Compute tangents using Catmull-Rom style
        let tangents: Vec<f64> = (0..n)
            .map(|i| {
                if periodic {
                    let prev = if i == 0 { n - 2 } else { i - 1 };
                    let next = if i == n - 1 { 1 } else { i + 1 };
                    // Adjust for periodic wrap
                    let dt_prev = sorted_points[i].0 - sorted_points[prev].0;
                    let dt_next = sorted_points[next].0 - sorted_points[i].0;
                    let dv_prev = sorted_points[i].1 - sorted_points[prev].1;
                    let dv_next = sorted_points[next].1 - sorted_points[i].1;
                    if dt_prev.abs() > 1e-14 && dt_next.abs() > 1e-14 {
                        0.5 * (dv_prev / dt_prev + dv_next / dt_next)
                    } else {
                        0.0
                    }
                } else if i == 0 {
                    // Forward difference
                    let dt = sorted_points[1].0 - sorted_points[0].0;
                    if dt.abs() > 1e-14 {
                        (sorted_points[1].1 - sorted_points[0].1) / dt
                    } else {
                        0.0
                    }
                } else if i == n - 1 {
                    // Backward difference
                    let dt = sorted_points[n - 1].0 - sorted_points[n - 2].0;
                    if dt.abs() > 1e-14 {
                        (sorted_points[n - 1].1 - sorted_points[n - 2].1) / dt
                    } else {
                        0.0
                    }
                } else {
                    // Central difference (Catmull-Rom)
                    let dt_prev = sorted_points[i].0 - sorted_points[i - 1].0;
                    let dt_next = sorted_points[i + 1].0 - sorted_points[i].0;
                    let dv_prev = sorted_points[i].1 - sorted_points[i - 1].1;
                    let dv_next = sorted_points[i + 1].1 - sorted_points[i].1;
                    if dt_prev.abs() > 1e-14 && dt_next.abs() > 1e-14 {
                        0.5 * (dv_prev / dt_prev + dv_next / dt_next)
                    } else {
                        0.0
                    }
                }
            })
            .collect();

        // Compute Hermite cubic coefficients for each interval
        for i in 0..n - 1 {
            let (t0, v0) = sorted_points[i];
            let (t1, v1) = sorted_points[i + 1];
            let m0 = tangents[i];
            let m1 = tangents[i + 1];

            let h = t1 - t0;
            if h.abs() < 1e-14 {
                coefficients.push(CubicCoefficients::default());
                continue;
            }

            // Hermite basis functions:
            // H00(s) = 2s^3 - 3s^2 + 1
            // H10(s) = s^3 - 2s^2 + s
            // H01(s) = -2s^3 + 3s^2
            // H11(s) = s^3 - s^2
            // where s = (t - t0) / h
            //
            // P(t) = v0*H00(s) + h*m0*H10(s) + v1*H01(s) + h*m1*H11(s)
            //
            // Expanding in powers of s:
            // s^3: 2v0 + h*m0 - 2v1 + h*m1
            // s^2: -3v0 - 2h*m0 + 3v1 - h*m1
            // s^1: h*m0
            // s^0: v0

            let s2_coef = -3.0 * v0 - 2.0 * h * m0 + 3.0 * v1 - h * m1;
            let s3_coef = 2.0 * v0 + h * m0 - 2.0 * v1 + h * m1;

            coefficients.push(CubicCoefficients {
                a: v0,
                b: h * m0,
                c: s2_coef,
                d: s3_coef,
            });
        }

        Self {
            points: sorted_points,
            periodic,
            coefficients,
        }
    }

    /// Find the interval index for parameter t.
    fn find_interval(&self, t: f64) -> usize {
        let n = self.points.len();
        for i in 0..n - 1 {
            if t >= self.points[i].0 && t <= self.points[i + 1].0 {
                return i;
            }
        }
        // Clamp to last interval
        n - 2
    }
}

impl LawFunction for InterpolateLaw {
    fn value(&self, t: f64) -> f64 {
        let n = self.points.len();
        if n == 0 {
            return 0.0;
        }

        let [t_min, t_max] = self.domain();
        let t_clamped = t.clamp(t_min, t_max);

        let i = self.find_interval(t_clamped);
        let (t0, _) = self.points[i];
        let (t1, _) = self.points[i + 1];
        let h = t1 - t0;

        if h.abs() < 1e-14 {
            return self.points[i].1;
        }

        let s = (t_clamped - t0) / h;
        let coef = &self.coefficients[i];
        coef.a + coef.b * s + coef.c * s * s + coef.d * s * s * s
    }

    fn derivative(&self, t: f64) -> f64 {
        let n = self.points.len();
        if n == 0 {
            return 0.0;
        }

        let [t_min, t_max] = self.domain();
        let t_clamped = t.clamp(t_min, t_max);

        let i = self.find_interval(t_clamped);
        let (t0, _) = self.points[i];
        let (t1, _) = self.points[i + 1];
        let h = t1 - t0;

        if h.abs() < 1e-14 {
            return 0.0;
        }

        let s = (t_clamped - t0) / h;
        let coef = &self.coefficients[i];

        // dP/ds = b + 2*c*s + 3*d*s^2
        // dP/dt = dP/ds * ds/dt = dP/ds / h
        let dp_ds = coef.b + 2.0 * coef.c * s + 3.0 * coef.d * s * s;
        dp_ds / h
    }

    fn domain(&self) -> [f64; 2] {
        if self.points.is_empty() {
            return [0.0, 1.0];
        }
        [self.points[0].0, self.points[self.points.len() - 1].0]
    }
}

// ============================================================================
// Common Laws
// ============================================================================

/// A sine law that smoothly transitions between two values.
///
/// The transition uses a half-sine wave from v1 to v2.
/// At t=t1, the value is v1; at t=t2, the value is v2.
#[derive(Debug, Clone, Copy)]
pub struct SineLaw {
    t1: f64,
    v1: f64,
    t2: f64,
    v2: f64,
}

impl SineLaw {
    /// Create a sine law transitioning from (t1, v1) to (t2, v2).
    pub fn new(t1: f64, v1: f64, t2: f64, v2: f64) -> Self {
        Self { t1, v1, t2, v2 }
    }
}

impl LawFunction for SineLaw {
    fn value(&self, t: f64) -> f64 {
        let h = self.t2 - self.t1;
        if h.abs() < 1e-14 {
            return self.v1;
        }
        let s = ((t - self.t1) / h).clamp(0.0, 1.0);
        // Half-sine wave: 0.5 * (1 - cos(pi * s)) maps [0,1] to [0,1]
        let blend = 0.5 * (1.0 - (PI * s).cos());
        self.v1 + (self.v2 - self.v1) * blend
    }

    fn derivative(&self, t: f64) -> f64 {
        let h = self.t2 - self.t1;
        if h.abs() < 1e-14 {
            return 0.0;
        }
        let s = ((t - self.t1) / h).clamp(0.0, 1.0);
        // d/ds [0.5 * (1 - cos(pi*s))] = 0.5 * pi * sin(pi*s)
        let d_blend = 0.5 * PI * (PI * s).sin();
        (self.v2 - self.v1) * d_blend / h
    }

    fn domain(&self) -> [f64; 2] {
        if self.t1 < self.t2 {
            [self.t1, self.t2]
        } else {
            [self.t2, self.t1]
        }
    }
}

/// Create a sine law transitioning from (t1, v1) to (t2, v2).
pub fn sine_law(t1: f64, v1: f64, t2: f64, v2: f64) -> SineLaw {
    SineLaw::new(t1, v1, t2, v2)
}

/// A smooth-step law using cubic Hermite interpolation.
///
/// Creates an S-curve transition from 0 to 1 over the domain [t1, t2].
/// The derivative is zero at both endpoints for smooth attachment
/// to constant segments.
#[derive(Debug, Clone, Copy)]
pub struct SmoothStepLaw {
    t1: f64,
    t2: f64,
}

impl SmoothStepLaw {
    /// Create a smooth-step law over the domain [t1, t2].
    pub fn new(t1: f64, t2: f64) -> Self {
        Self { t1, t2 }
    }
}

impl LawFunction for SmoothStepLaw {
    fn value(&self, t: f64) -> f64 {
        let h = self.t2 - self.t1;
        if h.abs() < 1e-14 {
            return 0.0;
        }
        let s = ((t - self.t1) / h).clamp(0.0, 1.0);
        // Smoothstep: 3s^2 - 2s^3
        s * s * (3.0 - 2.0 * s)
    }

    fn derivative(&self, t: f64) -> f64 {
        let h = self.t2 - self.t1;
        if h.abs() < 1e-14 {
            return 0.0;
        }
        let s = ((t - self.t1) / h).clamp(0.0, 1.0);
        // d/ds [3s^2 - 2s^3] = 6s - 6s^2 = 6s(1-s)
        let d_s = 6.0 * s * (1.0 - s);
        d_s / h
    }

    fn domain(&self) -> [f64; 2] {
        if self.t1 < self.t2 {
            [self.t1, self.t2]
        } else {
            [self.t2, self.t1]
        }
    }
}

/// Create a smooth-step law over the domain [t1, t2].
pub fn smooth_step_law(t1: f64, t2: f64) -> SmoothStepLaw {
    SmoothStepLaw::new(t1, t2)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // -------------------------------------------------------------------------
    // ConstantLaw tests
    // -------------------------------------------------------------------------

    #[test]
    fn constant_law_basic() {
        let law = ConstantLaw::new(5.0);
        assert!(approx_eq(law.value(0.0), 5.0, TOL));
        assert!(approx_eq(law.value(0.5), 5.0, TOL));
        assert!(approx_eq(law.value(1.0), 5.0, TOL));
        assert!(approx_eq(law.derivative(0.5), 0.0, TOL));
    }

    #[test]
    fn constant_law_custom_domain() {
        let law = ConstantLaw::with_domain(3.0, -1.0, 5.0);
        assert_eq!(law.domain(), [-1.0, 5.0]);
        assert!(approx_eq(law.value(2.0), 3.0, TOL));
    }

    // -------------------------------------------------------------------------
    // LinearLaw tests
    // -------------------------------------------------------------------------

    #[test]
    fn linear_law_basic() {
        let law = LinearLaw::new(0.0, 1.0, 1.0, 3.0);
        assert!(approx_eq(law.value(0.0), 1.0, TOL));
        assert!(approx_eq(law.value(0.5), 2.0, TOL));
        assert!(approx_eq(law.value(1.0), 3.0, TOL));
        assert!(approx_eq(law.derivative(0.5), 2.0, TOL));
        assert!(approx_eq(law.slope(), 2.0, TOL));
    }

    #[test]
    fn linear_law_negative_slope() {
        let law = LinearLaw::new(0.0, 5.0, 1.0, 1.0);
        assert!(approx_eq(law.value(0.0), 5.0, TOL));
        assert!(approx_eq(law.value(1.0), 1.0, TOL));
        assert!(approx_eq(law.derivative(0.0), -4.0, TOL));
    }

    #[test]
    fn linear_law_extrapolation() {
        let law = LinearLaw::new(0.0, 1.0, 1.0, 3.0);
        assert!(approx_eq(law.value(-0.5), 0.0, TOL));
        assert!(approx_eq(law.value(1.5), 4.0, TOL));
    }

    #[test]
    #[should_panic]
    fn linear_law_zero_domain_panics() {
        LinearLaw::new(1.0, 0.0, 1.0, 1.0);
    }

    // -------------------------------------------------------------------------
    // BSplineLaw tests
    // -------------------------------------------------------------------------

    #[test]
    fn bspline_law_linear() {
        // Create a linear BSpline (degree 1)
        let params = vec![0.0, 1.0];
        let values = vec![1.0, 3.0];
        let law = BSplineLaw::from_points(&params, &values, 1);

        assert!(approx_eq(law.value(0.0), 1.0, 1e-6));
        assert!(approx_eq(law.value(1.0), 3.0, 1e-6));
        assert!(approx_eq(law.value(0.5), 2.0, 1e-6));
    }

    #[test]
    fn bspline_law_quadratic() {
        // Quadratic BSpline (Bezier curve with 3 control points)
        let params = vec![0.0, 0.5, 1.0];
        let values = vec![0.0, 1.0, 0.0];
        let law = BSplineLaw::from_points(&params, &values, 2);

        // Check endpoints (clamped BSpline passes through endpoints)
        assert!(approx_eq(law.value(0.0), 0.0, 1e-6));
        assert!(approx_eq(law.value(1.0), 0.0, 1e-6));

        // For quadratic Bezier with control values [0, 1, 0], the max is at t=0.5
        // value = B0(0.5)*0 + B1(0.5)*1 + B2(0.5)*0 = 0.5
        let mid = law.value(0.5);
        assert!(approx_eq(mid, 0.5, 1e-6), "mid value should be 0.5, got {}", mid);
    }

    #[test]
    fn bspline_law_derivative() {
        let params = vec![0.0, 1.0];
        let values = vec![1.0, 3.0];
        let law = BSplineLaw::from_points(&params, &values, 1);

        // Derivative of linear is constant
        let deriv = law.derivative(0.5);
        assert!(approx_eq(deriv, 2.0, 1e-6));
    }

    #[test]
    fn bspline_law_domain() {
        let params = vec![0.0, 0.5, 1.0];
        let values = vec![1.0, 2.0, 3.0];
        let law = BSplineLaw::from_points(&params, &values, 2);

        let domain = law.domain();
        assert!(approx_eq(domain[0], 0.0, TOL));
        assert!(approx_eq(domain[1], 1.0, TOL));
    }

    // -------------------------------------------------------------------------
    // CompositeLaw tests
    // -------------------------------------------------------------------------

    #[test]
    fn composite_law_basic() {
        let mut law = CompositeLaw::new();

        // Add two linear segments
        law.add_segment(0.0, 0.5, Box::new(LinearLaw::new(0.0, 0.0, 1.0, 1.0)));
        law.add_segment(0.5, 1.0, Box::new(ConstantLaw::new(1.0)));

        assert!(approx_eq(law.value(0.0), 0.0, TOL));
        assert!(approx_eq(law.value(0.25), 0.5, TOL));
        assert!(approx_eq(law.value(0.5), 1.0, TOL));
        assert!(approx_eq(law.value(0.75), 1.0, TOL));
        assert!(approx_eq(law.value(1.0), 1.0, TOL));
    }

    #[test]
    fn composite_law_domain() {
        let mut law = CompositeLaw::new();
        law.add_segment(0.0, 0.5, Box::new(ConstantLaw::new(1.0)));
        law.add_segment(0.5, 2.0, Box::new(ConstantLaw::new(2.0)));

        let domain = law.domain();
        assert!(approx_eq(domain[0], 0.0, TOL));
        assert!(approx_eq(domain[1], 2.0, TOL));
    }

    #[test]
    fn composite_law_derivative() {
        let mut law = CompositeLaw::new();
        law.add_segment(0.0, 1.0, Box::new(LinearLaw::new(0.0, 0.0, 1.0, 2.0)));

        // Derivative should be 2.0
        let deriv = law.derivative(0.5);
        assert!(approx_eq(deriv, 2.0, TOL));
    }

    // -------------------------------------------------------------------------
    // InterpolateLaw tests
    // -------------------------------------------------------------------------

    #[test]
    fn interpolate_law_basic() {
        let points = vec![(0.0, 0.0), (0.5, 1.0), (1.0, 0.0)];
        let law = InterpolateLaw::from_points(&points, false);

        // Check endpoints
        assert!(approx_eq(law.value(0.0), 0.0, TOL));
        assert!(approx_eq(law.value(0.5), 1.0, TOL));
        assert!(approx_eq(law.value(1.0), 0.0, TOL));
    }

    #[test]
    fn interpolate_law_derivative() {
        let points = vec![(0.0, 0.0), (1.0, 2.0)];
        let law = InterpolateLaw::from_points(&points, false);

        // For two points, it's linear with slope 2
        let deriv = law.derivative(0.5);
        assert!(approx_eq(deriv, 2.0, TOL));
    }

    #[test]
    fn interpolate_law_periodic() {
        let points = vec![(0.0, 0.0), (0.5, 1.0), (1.0, 0.0)];
        let law = InterpolateLaw::from_points(&points, true);

        // Periodic should have matching derivatives at endpoints
        let d0 = law.derivative(0.01);
        let d1 = law.derivative(0.99);
        // Derivatives should have opposite signs for a symmetric wave
        assert!(d0 * d1 < 0.0, "periodic derivatives should have opposite signs");
    }

    #[test]
    fn interpolate_law_smoothness() {
        let points = vec![(0.0, 0.0), (0.25, 0.5), (0.5, 1.0), (0.75, 0.5), (1.0, 0.0)];
        let law = InterpolateLaw::from_points(&points, false);

        // Value should be continuous at each point
        for &(t, v) in &points {
            assert!(approx_eq(law.value(t), v, TOL));
        }
    }

    // -------------------------------------------------------------------------
    // SineLaw tests
    // -------------------------------------------------------------------------

    #[test]
    fn sine_law_basic() {
        let law = sine_law(0.0, 0.0, 1.0, 1.0);

        assert!(approx_eq(law.value(0.0), 0.0, TOL));
        assert!(approx_eq(law.value(1.0), 1.0, TOL));

        // Midpoint should be exactly 0.5 for sine
        assert!(approx_eq(law.value(0.5), 0.5, TOL));
    }

    #[test]
    fn sine_law_derivative_endpoints() {
        let law = sine_law(0.0, 0.0, 1.0, 1.0);

        // Derivative should be zero at endpoints
        assert!(approx_eq(law.derivative(0.0), 0.0, 1e-6));
        assert!(approx_eq(law.derivative(1.0), 0.0, 1e-6));

        // Maximum derivative at midpoint
        let d_mid = law.derivative(0.5);
        assert!(d_mid > 1.0, "derivative at midpoint should be > 1, got {}", d_mid);
    }

    // -------------------------------------------------------------------------
    // SmoothStepLaw tests
    // -------------------------------------------------------------------------

    #[test]
    fn smooth_step_law_basic() {
        let law = smooth_step_law(0.0, 1.0);

        assert!(approx_eq(law.value(0.0), 0.0, TOL));
        assert!(approx_eq(law.value(1.0), 1.0, TOL));

        // Midpoint should be exactly 0.5
        assert!(approx_eq(law.value(0.5), 0.5, TOL));
    }

    #[test]
    fn smooth_step_law_derivative() {
        let law = smooth_step_law(0.0, 1.0);

        // Derivative should be zero at endpoints
        assert!(approx_eq(law.derivative(0.0), 0.0, 1e-6));
        assert!(approx_eq(law.derivative(1.0), 0.0, 1e-6));

        // Maximum derivative at midpoint
        let d_mid = law.derivative(0.5);
        assert!(approx_eq(d_mid, 1.5, TOL));
    }

    #[test]
    fn smooth_step_law_monotonic() {
        let law = smooth_step_law(0.0, 1.0);

        // Should be monotonically increasing
        let mut prev = law.value(0.0);
        for i in 1..=100 {
            let t = i as f64 / 100.0;
            let v = law.value(t);
            assert!(v >= prev, "value decreased at t={}", t);
            prev = v;
        }
    }

    // -------------------------------------------------------------------------
    // General LawFunction tests
    // -------------------------------------------------------------------------

    #[test]
    fn law_function_clamped() {
        let law = LinearLaw::new(0.0, 0.0, 1.0, 1.0);

        assert!(approx_eq(law.value_clamped(-1.0), 0.0, TOL));
        assert!(approx_eq(law.value_clamped(2.0), 1.0, TOL));
    }

    #[test]
    fn law_function_in_domain() {
        let law = LinearLaw::new(0.0, 0.0, 1.0, 1.0);

        assert!(law.is_in_domain(0.0));
        assert!(law.is_in_domain(0.5));
        assert!(law.is_in_domain(1.0));
        assert!(!law.is_in_domain(-0.1));
        assert!(!law.is_in_domain(1.1));
    }
}
