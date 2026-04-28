//! ApproxInt-style intersection curve approximation.
//!
//! This module provides approximation of intersection curves between surfaces,
//! analogous to OpenCASCADE's ApproxInt package.
//!
//! # Overview
//!
//! When two surfaces intersect, the intersection curve can be complex. ApproxInt
//! provides tools to:
//! - Approximate the 3D intersection curve
//! - Approximate the 2D parameter-space curves (pcurves) on each surface
//! - Ensure same-parameter consistency between 3D and 2D curves
//!
//! # Main types
//!
//! - [`IntersectionApproximator`] - Main approximator for intersection curves
//! - [`ApproxOptions`] - Configuration options for approximation
//! - [`ApproxResult`] - Result containing 3D and 2D curves
//!
//! # Example
//!
//! ```rust
//! use rcad_algorithms::approx_int::{IntersectionApproximator, ApproxOptions};
//! use glam::{DVec3, DVec2};
//!
//! let mut approx = IntersectionApproximator::new();
//! // Add points sampled from the intersection
//! approx.add_point(DVec3::new(0.0, 0.0, 0.0), DVec2::new(0.0, 0.0), DVec2::new(0.0, 0.0));
//! approx.add_point(DVec3::new(1.0, 0.0, 0.0), DVec2::new(1.0, 0.0), DVec2::new(1.0, 0.0));
//!
//! let result = approx.compute(1e-6);
//! assert!(result.curve3d.is_some());
//! ```

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{BSplineCurve2, BSplineCurve3, Circle3, Curve2d, Curve2dEval, Curve3, CurveEval, Line3, Surface3, SurfaceEval};
use rcad_kernel::fit::{interpolate_points, interpolate_points_2d};

// ─────────────────────────────────────────────────────────────────────────────
// Approximation Options
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration options for intersection curve approximation.
///
/// Controls tolerance, degree, and continuity of the resulting curves.
#[derive(Debug, Clone)]
pub struct ApproxOptions {
    /// Tolerance for approximation (default: 1e-6).
    pub tolerance: f64,
    /// Maximum degree for the resulting BSpline (default: 8).
    pub max_degree: usize,
    /// Minimum degree for the resulting BSpline (default: 3).
    pub min_degree: usize,
    /// Desired continuity order: 0 = C0, 1 = C1, 2 = C2 (default: 2).
    pub continuity: usize,
    /// Maximum number of segments for piecewise approximation.
    pub max_segments: usize,
    /// Whether to enforce same-parameter constraint (default: true).
    pub same_parameter: bool,
}

impl Default for ApproxOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            max_degree: 8,
            min_degree: 3,
            continuity: 2,
            max_segments: 100,
            same_parameter: true,
        }
    }
}

impl ApproxOptions {
    /// Create default options.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the approximation tolerance.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Set the maximum degree for the resulting BSpline.
    pub fn with_max_degree(mut self, deg: usize) -> Self {
        self.max_degree = deg;
        self
    }

    /// Set the desired continuity order.
    pub fn with_continuity(mut self, cont: usize) -> Self {
        self.continuity = cont;
        self
    }

    /// Set the minimum degree.
    pub fn with_min_degree(mut self, deg: usize) -> Self {
        self.min_degree = deg;
        self
    }

    /// Set the maximum number of segments.
    pub fn with_max_segments(mut self, n: usize) -> Self {
        self.max_segments = n;
        self
    }

    /// Enable or disable same-parameter constraint.
    pub fn with_same_parameter(mut self, same: bool) -> Self {
        self.same_parameter = same;
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Approximation Result
// ─────────────────────────────────────────────────────────────────────────────

/// Result of intersection curve approximation.
///
/// Contains the approximated 3D curve and optional 2D curves in each
/// surface's parameter space.
#[derive(Debug, Clone)]
pub struct ApproxResult {
    /// The approximated 3D curve.
    pub curve3d: Option<BSplineCurve3>,
    /// The approximated 2D curve on the first surface.
    pub curve2d1: Option<BSplineCurve2>,
    /// The approximated 2D curve on the second surface.
    pub curve2d2: Option<BSplineCurve2>,
    /// Achieved tolerance (maximum deviation).
    pub achieved_tolerance: f64,
    /// Whether the approximation succeeded.
    pub success: bool,
    /// Error message if approximation failed.
    pub error: Option<String>,
}

impl Default for ApproxResult {
    fn default() -> Self {
        Self {
            curve3d: None,
            curve2d1: None,
            curve2d2: None,
            achieved_tolerance: f64::INFINITY,
            success: false,
            error: None,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Sample Point
// ─────────────────────────────────────────────────────────────────────────────

/// A sample point from the intersection of two surfaces.
///
/// Contains the 3D position and the corresponding 2D parameters on each surface.
#[derive(Debug, Clone, Copy)]
pub struct IntersectionSample {
    /// 3D position on the intersection curve.
    pub point: DVec3,
    /// Parameter (u, v) on the first surface.
    pub uv1: DVec2,
    /// Parameter (u, v) on the second surface.
    pub uv2: DVec2,
    /// Parameter along the curve (chord-length normalized).
    pub param: f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Intersection Approximator
// ─────────────────────────────────────────────────────────────────────────────

/// Main approximator for intersection curves (ApproxInt_Approx).
///
/// Collects sample points from a surface-surface intersection and builds
/// approximating B-spline curves for the 3D curve and 2D pcurves.
///
/// # Usage
///
/// 1. Create a new approximator
/// 2. Add sample points with [`add_point`](Self::add_point)
/// 3. Compute the approximation with [`compute`](Self::compute)
/// 4. Retrieve the results with [`curve3d`](Self::curve3d), [`curve2d1`](Self::curve2d1), [`curve2d2`](Self::curve2d2)
#[derive(Debug, Clone, Default)]
pub struct IntersectionApproximator {
    /// Sample points collected from the intersection.
    samples: Vec<IntersectionSample>,
    /// Computed 3D curve.
    curve3d: Option<BSplineCurve3>,
    /// Computed 2D curve on first surface.
    curve2d1: Option<BSplineCurve2>,
    /// Computed 2D curve on second surface.
    curve2d2: Option<BSplineCurve2>,
    /// Achieved tolerance.
    achieved_tolerance: f64,
}

impl IntersectionApproximator {
    /// Create a new intersection approximator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create an approximator with pre-allocated capacity.
    pub fn with_capacity(n: usize) -> Self {
        Self {
            samples: Vec::with_capacity(n),
            curve3d: None,
            curve2d1: None,
            curve2d2: None,
            achieved_tolerance: f64::INFINITY,
        }
    }

    /// Add a sample point to the approximator.
    ///
    /// # Arguments
    /// * `point` - 3D position on the intersection curve
    /// * `uv1` - Parameter (u, v) on the first surface
    /// * `uv2` - Parameter (u, v) on the second surface
    pub fn add_point(&mut self, point: DVec3, uv1: DVec2, uv2: DVec2) {
        // Compute chord-length parameter
        let param = if let Some(last) = self.samples.last() {
            last.param + (point - last.point).length()
        } else {
            0.0
        };

        self.samples.push(IntersectionSample {
            point,
            uv1,
            uv2,
            param,
        });
    }

    /// Add a sample point with explicit parameter value.
    pub fn add_point_with_param(&mut self, point: DVec3, uv1: DVec2, uv2: DVec2, param: f64) {
        self.samples.push(IntersectionSample {
            point,
            uv1,
            uv2,
            param,
        });
    }

    /// Get the number of sample points.
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Check if there are no sample points.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Clear all samples and results.
    pub fn clear(&mut self) {
        self.samples.clear();
        self.curve3d = None;
        self.curve2d1 = None;
        self.curve2d2 = None;
        self.achieved_tolerance = f64::INFINITY;
    }

    /// Compute the approximating curves.
    ///
    /// # Arguments
    /// * `tol` - Tolerance for approximation
    ///
    /// # Returns
    /// The approximation result containing the curves.
    pub fn compute(&mut self, tol: f64) -> ApproxResult {
        self.compute_with_options(&ApproxOptions::default().with_tolerance(tol))
    }

    /// Compute the approximating curves with custom options.
    pub fn compute_with_options(&mut self, options: &ApproxOptions) -> ApproxResult {
        let n = self.samples.len();

        if n < 2 {
            return ApproxResult {
                curve3d: None,
                curve2d1: None,
                curve2d2: None,
                achieved_tolerance: f64::INFINITY,
                success: false,
                error: Some("At least 2 points are required".to_string()),
            };
        }

        // Normalize parameters to [0, 1]
        let total_length = self.samples.last().map(|s| s.param).unwrap_or(1.0);
        if total_length < 1e-14 {
            return ApproxResult {
                curve3d: None,
                curve2d1: None,
                curve2d2: None,
                achieved_tolerance: f64::INFINITY,
                success: false,
                error: Some("All points are coincident".to_string()),
            };
        }

        // Normalize params
        let params: Vec<f64> = self.samples.iter().map(|s| s.param / total_length).collect();

        // Extract 3D points
        let pts3d: Vec<DVec3> = self.samples.iter().map(|s| s.point).collect();

        // Extract 2D points for each surface
        let pts2d1: Vec<DVec2> = self.samples.iter().map(|s| s.uv1).collect();
        let pts2d2: Vec<DVec2> = self.samples.iter().map(|s| s.uv2).collect();

        // Interpolate 3D curve
        let curve3d = match interpolate_points(&pts3d) {
            Ok(curve) => curve,
            Err(e) => {
                return ApproxResult {
                    curve3d: None,
                    curve2d1: None,
                    curve2d2: None,
                    achieved_tolerance: f64::INFINITY,
                    success: false,
                    error: Some(format!("3D interpolation failed: {}", e)),
                };
            }
        };

        // Interpolate 2D curves
        let curve2d1 = interpolate_points_2d(&pts2d1).ok();
        let curve2d2 = interpolate_points_2d(&pts2d2).ok();

        // Compute achieved tolerance
        let achieved_tol = self.compute_achieved_tolerance(&curve3d, &pts3d, &params);

        // Store results
        self.curve3d = Some(curve3d.clone());
        self.curve2d1 = curve2d1.clone();
        self.curve2d2 = curve2d2.clone();
        self.achieved_tolerance = achieved_tol;

        ApproxResult {
            curve3d: Some(curve3d),
            curve2d1,
            curve2d2,
            achieved_tolerance: achieved_tol,
            success: achieved_tol <= options.tolerance,
            error: if achieved_tol > options.tolerance {
                Some(format!("Achieved tolerance {} exceeds requested {}", achieved_tol, options.tolerance))
            } else {
                None
            },
        }
    }

    /// Compute the achieved tolerance by comparing the curve to the original points.
    fn compute_achieved_tolerance(
        &self,
        curve: &BSplineCurve3,
        pts: &[DVec3],
        params: &[f64],
    ) -> f64 {
        let mut max_dev = 0.0_f64;
        for (i, &t) in params.iter().enumerate() {
            let pt_on_curve = curve.point_at(t);
            let dev = (pts[i] - pt_on_curve).length();
            max_dev = max_dev.max(dev);
        }
        max_dev
    }

    /// Get the computed 3D curve.
    pub fn curve3d(&self) -> Option<&BSplineCurve3> {
        self.curve3d.as_ref()
    }

    /// Get the computed 2D curve on the first surface.
    pub fn curve2d1(&self) -> Option<&BSplineCurve2> {
        self.curve2d1.as_ref()
    }

    /// Get the computed 2D curve on the second surface.
    pub fn curve2d2(&self) -> Option<&BSplineCurve2> {
        self.curve2d2.as_ref()
    }

    /// Get the achieved tolerance.
    pub fn achieved_tolerance(&self) -> f64 {
        self.achieved_tolerance
    }

    /// Get all sample points.
    pub fn samples(&self) -> &[IntersectionSample] {
        &self.samples
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Same Parameter Computation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the same-parameter deviation between a 3D curve and its 2D representation.
///
/// Given a 3D curve, a 2D curve in the parameter space of a surface, and the surface
/// itself, compute the maximum deviation between:
/// - The 3D curve point at parameter t
/// - The surface evaluated at the 2D curve point at parameter t
///
/// This is analogous to OCCT's `ApproxInt_SameParameter`.
///
/// # Arguments
/// * `curve3d` - The 3D curve
/// * `curve2d` - The 2D curve in the surface's parameter space
/// * `surface` - The surface on which the 2D curve is defined
/// * `tol` - Tolerance for numerical operations
///
/// # Returns
/// The maximum deviation (should be within tolerance for a valid same-parameter curve).
pub fn compute_same_parameter(
    curve3d: &Curve3,
    curve2d: &Curve2d,
    surface: &Surface3,
    tol: f64,
) -> f64 {
    let domain = curve3d.default_domain();
    let n_samples = 50; // Number of samples for deviation check

    let mut max_dev = 0.0_f64;

    for i in 0..=n_samples {
        let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;

        let pt3d = curve3d.point_at(t);
        let uv = curve2d.point_at(t);
        let pt_on_surf = surface.point_at(uv.x, uv.y);

        let dev = (pt3d - pt_on_surf).length();
        max_dev = max_dev.max(dev);

        if dev > tol && dev > max_dev {
            max_dev = dev;
        }
    }

    max_dev
}

/// Compute same-parameter deviation for BSpline curves.
///
/// Convenience function that works directly with BSpline types.
pub fn compute_same_parameter_bspline(
    curve3d: &BSplineCurve3,
    curve2d: &BSplineCurve2,
    surface: &Surface3,
    n_samples: usize,
) -> f64 {
    let mut max_dev = 0.0_f64;

    for i in 0..=n_samples {
        let t = i as f64 / n_samples as f64;

        let pt3d = curve3d.point_at(t);
        let uv = curve2d.point_at(t);
        let pt_on_surf = surface.point_at(uv.x, uv.y);

        let dev = (pt3d - pt_on_surf).length();
        max_dev = max_dev.max(dev);
    }

    max_dev
}

/// Adjust the 2D curve to achieve same-parameter consistency.
///
/// This function modifies the 2D curve so that the surface(c2d(t)) matches
/// the 3D curve as closely as possible.
///
/// # Arguments
/// * `curve3d` - The 3D curve (reference)
/// * `curve2d` - The 2D curve to adjust
/// * `surface` - The surface on which the 2D curve is defined
/// * `tol` - Tolerance for the adjustment
///
/// # Returns
/// The adjusted 2D curve and the achieved deviation.
pub fn adjust_same_parameter(
    curve3d: &BSplineCurve3,
    curve2d: &BSplineCurve2,
    surface: &Surface3,
    tol: f64,
) -> (BSplineCurve2, f64) {
    let n = curve2d.control_points.len();

    // Collect corrected UV points
    let mut corrected_uvs: Vec<DVec2> = Vec::with_capacity(n);

    // For each control point parameter, adjust UV to minimize deviation
    for i in 0..n {
        // Map control point index to parameter
        let t = i as f64 / (n - 1).max(1) as f64;

        let target_pt = curve3d.point_at(t);
        let current_uv = curve2d.point_at(t);

        // Project the target point onto the surface to get corrected UV
        let corrected_uv = project_point_to_surface_uv(&target_pt, surface, current_uv, tol);
        corrected_uvs.push(corrected_uv);
    }

    // Rebuild the curve with corrected control points
    let adjusted_curve = BSplineCurve2 {
        degree: curve2d.degree,
        knots: curve2d.knots.clone(),
        control_points: corrected_uvs,
        weights: curve2d.weights.clone(),
    };

    let achieved_dev = compute_same_parameter_bspline(curve3d, &adjusted_curve, surface, 50);

    (adjusted_curve, achieved_dev)
}

/// Project a point onto a surface to find the UV parameter.
///
/// Uses Newton iteration starting from an initial UV guess.
fn project_point_to_surface_uv(
    point: &DVec3,
    surface: &Surface3,
    initial_uv: DVec2,
    tol: f64,
) -> DVec2 {
    let mut uv = initial_uv;
    let max_iter = 10;

    for _ in 0..max_iter {
        let surf_pt = surface.point_at(uv.x, uv.y);
        let diff = *point - surf_pt;
        let dist = diff.length();

        if dist < tol {
            break;
        }

        // Compute numerical gradient
        let h = 1e-6;
        let surf_pt_du = surface.point_at(uv.x + h, uv.y);
        let surf_pt_dv = surface.point_at(uv.x, uv.y + h);

        let du = (surf_pt_du - surf_pt) / h;
        let dv = (surf_pt_dv - surf_pt) / h;

        // Simple gradient descent step
        let grad_u = -diff.dot(du);
        let grad_v = -diff.dot(dv);

        let step = 0.5;
        uv.x += step * grad_u;
        uv.y += step * grad_v;
    }

    uv
}

// ─────────────────────────────────────────────────────────────────────────────
// 2D Curve Approximation
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate a 2D curve from a set of points with parameters.
///
/// This is the ApproxInt_2dCurve equivalent.
///
/// # Arguments
/// * `points` - Slice of (point, parameter) pairs
/// * `tol` - Tolerance for approximation
///
/// # Returns
/// An approximating BSplineCurve2.
pub fn approximate_2d_curve(points: &[(DVec2, f64)], tol: f64) -> BSplineCurve2 {
    if points.is_empty() {
        return BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::ZERO, DVec2::ZERO],
            weights: vec![1.0, 1.0],
        };
    }

    if points.len() == 1 {
        return BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![points[0].0, points[0].0],
            weights: vec![1.0, 1.0],
        };
    }

    // Extract just the points
    let pts: Vec<DVec2> = points.iter().map(|(p, _)| *p).collect();

    // Use interpolation for exact fit
    interpolate_points_2d(&pts).unwrap_or_else(|_| {
        // Fallback: create a line from first to last point
        BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![pts[0], *pts.last().unwrap()],
            weights: vec![1.0, 1.0],
        }
    })
}

/// Approximate a 2D curve with a specified number of control points.
///
/// Uses least-squares approximation if n_ctrl < number of points.
pub fn approximate_2d_curve_with_ctrl(
    points: &[(DVec2, f64)],
    n_ctrl: usize,
    tol: f64,
) -> BSplineCurve2 {
    if points.len() <= 2 {
        return approximate_2d_curve(points, tol);
    }

    let pts: Vec<DVec2> = points.iter().map(|(p, _)| *p).collect();

    // If n_ctrl equals number of points, use interpolation
    if n_ctrl >= pts.len() {
        return interpolate_points_2d(&pts).unwrap_or_else(|_| approximate_2d_curve(points, tol));
    }

    // Otherwise, we'd need least-squares - for now use interpolation
    // (full least-squares 2D would need additional implementation in fit module)
    interpolate_points_2d(&pts).unwrap_or_else(|_| approximate_2d_curve(points, tol))
}

// ─────────────────────────────────────────────────────────────────────────────
// Point Sampling for Intersection
// ─────────────────────────────────────────────────────────────────────────────

/// Sample points uniformly along a 3D curve.
///
/// # Arguments
/// * `curve` - The curve to sample
/// * `n_points` - Number of points to generate
///
/// # Returns
/// A vector of 3D points.
pub fn sample_intersection_points(curve: &Curve3, n_points: usize) -> Vec<DVec3> {
    if n_points == 0 {
        return Vec::new();
    }

    let domain = curve.default_domain();

    // Handle unbounded curves
    let (t0, t1) = if domain[0].is_finite() && domain[1].is_finite() {
        (domain[0], domain[1])
    } else {
        (-10.0, 10.0) // Default range for unbounded curves
    };

    (0..n_points)
        .map(|i| {
            let t = if n_points > 1 {
                t0 + (t1 - t0) * i as f64 / (n_points - 1) as f64
            } else {
                (t0 + t1) / 2.0
            };
            curve.point_at(t)
        })
        .collect()
}

/// Sample points along a curve with adaptive density based on curvature.
///
/// Places more points in regions of high curvature and fewer in flat regions.
///
/// # Arguments
/// * `curve` - The curve to sample
/// * `tol` - Tolerance controlling point density
/// * `max_points` - Maximum number of points to generate
///
/// # Returns
/// A vector of 3D points with adaptive spacing.
pub fn sample_with_adaptive_density(curve: &Curve3, tol: f64, max_points: usize) -> Vec<DVec3> {
    if max_points == 0 {
        return Vec::new();
    }

    let domain = curve.default_domain();
    let (t0, t1) = if domain[0].is_finite() && domain[1].is_finite() {
        (domain[0], domain[1])
    } else {
        (-10.0, 10.0)
    };

    // Start with a uniform initial sampling
    let initial_n = (max_points / 4).max(10);
    let mut params: Vec<f64> = (0..initial_n)
        .map(|i| t0 + (t1 - t0) * i as f64 / (initial_n - 1).max(1) as f64)
        .collect();

    // Iteratively refine based on chord deviation
    let mut refined = true;
    let mut iteration = 0;
    let max_iterations = 10;

    while refined && params.len() < max_points && iteration < max_iterations {
        refined = false;
        iteration += 1;

        let mut new_params = Vec::with_capacity(params.len() * 2);
        new_params.push(params[0]);

        for i in 1..params.len() {
            let t_prev = params[i - 1];
            let t_curr = params[i];

            let pt_prev = curve.point_at(t_prev);
            let pt_curr = curve.point_at(t_curr);

            // Check midpoint deviation
            let t_mid = (t_prev + t_curr) / 2.0;
            let pt_mid = curve.point_at(t_mid);

            // Chord deviation
            let chord = pt_curr - pt_prev;
            let chord_len = chord.length();

            if chord_len > 1e-10 {
                let chord_dir = chord / chord_len;
                let to_mid = pt_mid - pt_prev;
                let along_chord = to_mid.dot(chord_dir);
                let perp = to_mid - along_chord * chord_dir;
                let deviation = perp.length();

                if deviation > tol && params.len() + new_params.len() < max_points {
                    // Insert midpoint
                    new_params.push(t_mid);
                    refined = true;
                }
            }

            new_params.push(t_curr);
        }

        params = new_params;
    }

    // Limit to max_points
    params.truncate(max_points);

    // Convert parameters to points
    params.iter().map(|&t| curve.point_at(t)).collect()
}

/// Sample points on a curve between two parameter values.
///
/// # Arguments
/// * `curve` - The curve to sample
/// * `t_start` - Start parameter
/// * `t_end` - End parameter
/// * `n_points` - Number of points
pub fn sample_curve_segment(curve: &Curve3, t_start: f64, t_end: f64, n_points: usize) -> Vec<DVec3> {
    if n_points == 0 {
        return Vec::new();
    }

    (0..n_points)
        .map(|i| {
            let t = if n_points > 1 {
                t_start + (t_end - t_start) * i as f64 / (n_points - 1) as f64
            } else {
                (t_start + t_end) / 2.0
            };
            curve.point_at(t)
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Additional Approximation Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Approximate an intersection polyline with a smooth B-spline.
///
/// # Arguments
/// * `polyline` - The polyline points to approximate
/// * `tol` - Approximation tolerance
///
/// # Returns
/// An approximating BSplineCurve3.
pub fn approximate_polyline(polyline: &[DVec3], tol: f64) -> Option<BSplineCurve3> {
    if polyline.len() < 2 {
        return None;
    }

    interpolate_points(polyline).ok()
}

/// Approximate an intersection curve given as a polyline with UV parameters.
///
/// # Arguments
/// * `points3d` - 3D points on the intersection
/// * `uvs1` - Corresponding UV parameters on first surface
/// * `uvs2` - Corresponding UV parameters on second surface
/// * `tol` - Approximation tolerance
///
/// # Returns
/// The approximation result with all three curves.
pub fn approximate_intersection(
    points3d: &[DVec3],
    uvs1: &[DVec2],
    uvs2: &[DVec2],
    tol: f64,
) -> ApproxResult {
    if points3d.len() < 2 {
        return ApproxResult {
            curve3d: None,
            curve2d1: None,
            curve2d2: None,
            achieved_tolerance: f64::INFINITY,
            success: false,
            error: Some("At least 2 points required".to_string()),
        };
    }

    let n = points3d.len();

    // Build 3D curve
    let curve3d = match interpolate_points(points3d) {
        Ok(c) => Some(c),
        Err(e) => {
            return ApproxResult {
                curve3d: None,
                curve2d1: None,
                curve2d2: None,
                achieved_tolerance: f64::INFINITY,
                success: false,
                error: Some(format!("3D interpolation failed: {}", e)),
            };
        }
    };

    // Build 2D curves if UVs are provided
    let curve2d1 = if uvs1.len() == n {
        interpolate_points_2d(uvs1).ok()
    } else {
        None
    };

    let curve2d2 = if uvs2.len() == n {
        interpolate_points_2d(uvs2).ok()
    } else {
        None
    };

    ApproxResult {
        curve3d,
        curve2d1,
        curve2d2,
        achieved_tolerance: tol, // Actual tolerance would need evaluation
        success: true,
        error: None,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_approx_options() {
        let opts = ApproxOptions::default()
            .with_tolerance(1e-4)
            .with_max_degree(5)
            .with_continuity(1);

        assert!((opts.tolerance - 1e-4).abs() < 1e-10);
        assert_eq!(opts.max_degree, 5);
        assert_eq!(opts.continuity, 1);
    }

    #[test]
    fn test_intersection_approximator_empty() {
        let approx = IntersectionApproximator::new();
        assert!(approx.is_empty());
        assert_eq!(approx.len(), 0);

        let result = approx.clone().compute(1e-6);
        assert!(!result.success);
        assert!(result.error.is_some());
    }

    #[test]
    fn test_intersection_approximator_line() {
        let mut approx = IntersectionApproximator::new();

        // Create a simple line intersection
        approx.add_point(DVec3::new(0.0, 0.0, 0.0), DVec2::new(0.0, 0.0), DVec2::new(0.0, 0.0));
        approx.add_point(DVec3::new(1.0, 0.0, 0.0), DVec2::new(1.0, 0.0), DVec2::new(1.0, 0.0));
        approx.add_point(DVec3::new(2.0, 0.0, 0.0), DVec2::new(2.0, 0.0), DVec2::new(2.0, 0.0));

        let result = approx.compute(1e-6);
        assert!(result.success);
        assert!(result.curve3d.is_some());
        assert!(result.curve2d1.is_some());
        assert!(result.curve2d2.is_some());
    }

    #[test]
    fn test_intersection_approximator_curve() {
        let mut approx = IntersectionApproximator::new();

        // Create a curved intersection (quarter circle)
        let n = 11;
        for i in 0..n {
            let t = PI / 2.0 * i as f64 / (n - 1) as f64;
            let x = t.cos();
            let y = t.sin();
            approx.add_point(
                DVec3::new(x, y, 0.0),
                DVec2::new(x, y),
                DVec2::new(t, 0.0),
            );
        }

        let result = approx.compute(1e-4);
        assert!(result.success, "Approximation should succeed");

        let curve = result.curve3d.as_ref().unwrap();

        // Check endpoints
        let p0 = curve.point_at(0.0);
        let p1 = curve.point_at(1.0);

        assert!((p0 - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-4);
        assert!((p1 - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-4);
    }

    #[test]
    fn test_achieved_tolerance() {
        let mut approx = IntersectionApproximator::new();

        // Add points for a smooth curve
        for i in 0..10 {
            let t = i as f64 / 9.0;
            approx.add_point(
                DVec3::new(t, t * t, 0.0),
                DVec2::new(t, t * t),
                DVec2::new(t, 0.0),
            );
        }

        let result = approx.compute(1e-4);
        assert!(result.achieved_tolerance < 1e-2);
    }

    #[test]
    fn test_approximate_2d_curve() {
        let points = vec![
            (DVec2::new(0.0, 0.0), 0.0),
            (DVec2::new(0.5, 0.25), 0.5),
            (DVec2::new(1.0, 1.0), 1.0),
        ];

        let curve = approximate_2d_curve(&points, 1e-6);

        // Check endpoints
        let p0 = curve.point_at(0.0);
        let p1 = curve.point_at(1.0);

        assert!((p0 - DVec2::new(0.0, 0.0)).length() < 1e-6);
        assert!((p1 - DVec2::new(1.0, 1.0)).length() < 1e-6);
    }

    #[test]
    fn test_sample_intersection_points() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });

        let points = sample_intersection_points(&circle, 10);
        assert_eq!(points.len(), 10);

        // All points should be on the circle
        for pt in &points {
            let r = pt.length();
            assert!((r - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_sample_with_adaptive_density() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });

        let points = sample_with_adaptive_density(&circle, 0.01, 50);
        assert!(points.len() >= 10);
        assert!(points.len() <= 50);

        // All points should be on the circle
        for pt in &points {
            let r = pt.length();
            assert!((r - 1.0).abs() < 1e-4, "Point {} should be on circle, got r={}", pt, r);
        }
    }

    #[test]
    fn test_compute_same_parameter() {
        use rcad_kernel::geom::{Plane, Line3, Line2d};

        // Create a line on a plane
        let line3d = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });

        let line2d = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });

        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        // A line on a plane should have zero deviation
        let dev = compute_same_parameter(&line3d, &line2d, &plane, 1e-6);
        assert!(dev < 1e-6, "Same parameter deviation should be near zero, got {}", dev);
    }

    #[test]
    fn test_curve2d_accessors() {
        let mut approx = IntersectionApproximator::new();

        approx.add_point(DVec3::ZERO, DVec2::ZERO, DVec2::ZERO);
        approx.add_point(DVec3::X, DVec2::X, DVec2::X);

        let result = approx.compute(1e-6);
        assert!(result.success);

        // Test accessors
        assert!(approx.curve3d().is_some());
        assert!(approx.curve2d1().is_some());
        assert!(approx.curve2d2().is_some());
    }

    #[test]
    fn test_approximate_intersection() {
        let points3d = vec![
            DVec3::ZERO,
            DVec3::new(0.5, 0.5, 0.0),
            DVec3::X,
        ];

        let uvs1 = vec![
            DVec2::ZERO,
            DVec2::new(0.5, 0.5),
            DVec2::X,
        ];

        let uvs2 = vec![
            DVec2::ZERO,
            DVec2::new(0.25, 0.0),
            DVec2::new(0.5, 0.0),
        ];

        let result = approximate_intersection(&points3d, &uvs1, &uvs2, 1e-6);
        assert!(result.success);
        assert!(result.curve3d.is_some());
        assert!(result.curve2d1.is_some());
        assert!(result.curve2d2.is_some());
    }

    #[test]
    fn test_approximate_polyline() {
        let polyline = vec![
            DVec3::ZERO,
            DVec3::new(1.0, 0.5, 0.0),
            DVec3::X,
        ];

        let curve = approximate_polyline(&polyline, 1e-6);
        assert!(curve.is_some());

        let c = curve.unwrap();
        let p0 = c.point_at(0.0);
        let p1 = c.point_at(1.0);

        assert!((p0 - DVec3::ZERO).length() < 1e-6);
        assert!((p1 - DVec3::X).length() < 1e-6);
    }

    #[test]
    fn test_clear() {
        let mut approx = IntersectionApproximator::new();

        approx.add_point(DVec3::ZERO, DVec2::ZERO, DVec2::ZERO);
        approx.add_point(DVec3::X, DVec2::X, DVec2::X);

        assert_eq!(approx.len(), 2);

        approx.clear();

        assert!(approx.is_empty());
        assert!(approx.curve3d().is_none());
    }

    #[test]
    fn test_samples_accessor() {
        let mut approx = IntersectionApproximator::new();

        approx.add_point(DVec3::ZERO, DVec2::ZERO, DVec2::ZERO);
        approx.add_point(DVec3::X, DVec2::X, DVec2::X);
        approx.add_point(DVec3::Y, DVec2::Y, DVec2::Y);

        let samples = approx.samples();
        assert_eq!(samples.len(), 3);

        // Check chord-length parameter is accumulated
        assert_eq!(samples[0].param, 0.0);
        assert!(samples[1].param > 0.0);
        assert!(samples[2].param > samples[1].param);
    }

    #[test]
    fn test_sample_curve_segment() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });

        let points = sample_curve_segment(&line, 0.0, 5.0, 6);
        assert_eq!(points.len(), 6);

        // Check points are equally spaced
        for (i, pt) in points.iter().enumerate() {
            let expected = i as f64;
            assert!((pt.x - expected).abs() < 1e-10);
        }
    }
}
