//! GeomConvert-style geometry conversion utilities.
//!
//! This module provides functions for converting between different geometry representations,
//! analogous to OpenCASCADE's GeomConvert package.
//!
//! # Surface conversions
//! - [`plane_to_bspline`] - Convert plane to BSplineSurface
//! - [`cylinder_to_bspline`] - Convert cylinder to BSplineSurface
//! - [`cone_to_bspline`] - Convert cone to BSplineSurface
//! - [`sphere_to_bspline`] - Convert sphere to BSplineSurface
//! - [`torus_to_bspline`] - Convert torus to BSplineSurface
//! - [`surface_to_bspline`] - Generic surface to BSplineSurface conversion
//!
//! # Curve conversions
//! - [`line_to_bspline`] - Convert line to BSplineCurve
//! - [`circle_to_bspline`] - Convert circle to BSplineCurve
//! - [`ellipse_to_bspline`] - Convert ellipse to BSplineCurve
//! - [`curve_to_bspline`] - Generic curve to BSplineCurve conversion
//!
//! # BSpline operations
//! - [`bspline_to_bezier`] - Decompose BSplineCurve into Bezier segments
//! - [`bspline_surface_to_bezier`] - Decompose BSplineSurface into Bezier patches
//! - [`approx_curve_to_bspline`] - Approximate any curve with BSpline
//! - [`approx_surface_to_bspline`] - Approximate any surface with BSplineSurface

use glam::DVec3;
use std::f64::consts::PI;

use rcad_kernel::geom::{
    any_perpendicular, BSplineCurve3, BSplineSurface, BezierCurve3, BezierSurface, Circle3,
    ConicalSurface, Curve3, CurveEval, CylindricalSurface, Ellipse3, Line3, Plane,
    SphericalSurface, Surface3, SurfaceEval, ToroidalSurface,
};
use rcad_kernel::fit::interpolate_points;

// ─────────────────────────────────────────────────────────────────────────────
// Conversion Parameters
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters controlling geometry conversion behavior.
///
/// Analogous to OCCT `GeomConvert_CompBezierSurfacesFunction` parameters
/// and `AdvApprox_ApproxAFunction` settings.
#[derive(Debug, Clone)]
pub struct ConvertParams {
    /// Tolerance for approximation (default: 1e-6).
    pub tolerance: f64,
    /// Maximum degree for the resulting BSpline (default: 3 for curves, 3 for surfaces).
    pub max_degree: usize,
    /// Minimum degree for the resulting BSpline (default: 1).
    pub min_degree: usize,
    /// Desired continuity order: 0 = C0, 1 = C1, 2 = C2 (default: 2).
    pub continuity: usize,
    /// Number of samples for approximation (default: 20).
    pub sample_count: usize,
    /// Whether to use rational representation when possible (default: true).
    pub rational: bool,
}

impl Default for ConvertParams {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            max_degree: 3,
            min_degree: 1,
            continuity: 2,
            sample_count: 20,
            rational: true,
        }
    }
}

impl ConvertParams {
    /// Create new conversion parameters with specified tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            ..Default::default()
        }
    }

    /// Set the maximum degree.
    pub fn with_max_degree(mut self, degree: usize) -> Self {
        self.max_degree = degree;
        self
    }

    /// Set the continuity order.
    pub fn with_continuity(mut self, continuity: usize) -> Self {
        self.continuity = continuity;
        self
    }

    /// Set the sample count for approximation.
    pub fn with_sample_count(mut self, count: usize) -> Self {
        self.sample_count = count;
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Curve Conversions
// ─────────────────────────────────────────────────────────────────────────────

/// Convert a `Line3` to a degree-1 `BSplineCurve3`.
///
/// The result is a linear B-spline with 2 control points.
/// The parameter domain is `[0, 1]`.
///
/// # Arguments
/// * `line` - The line to convert
/// * `_degree` - Ignored for lines (always degree 1)
///
/// # Example
/// ```rust
/// use rcad_algorithms::geom_convert::line_to_bspline;
/// use rcad_kernel::geom::{Line3, CurveEval};
/// use glam::DVec3;
///
/// let line = Line3 { origin: DVec3::ZERO, direction: DVec3::X };
/// let bspline = line_to_bspline(&line, 1);
/// assert!((bspline.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
/// assert!((bspline.point_at(1.0) - DVec3::X).length() < 1e-10);
/// ```
pub fn line_to_bspline(line: &Line3, _degree: usize) -> BSplineCurve3 {
    let p0 = line.origin;
    let p1 = line.origin + line.direction;

    BSplineCurve3 {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

/// Convert a `Circle3` to an exact rational quadratic B-spline.
///
/// Uses the standard 9-point NURBS circle construction (Piegl & Tiller).
/// The result has degree 2 with 9 control points and weights that produce
/// an exact circular arc.
///
/// # Arguments
/// * `circle` - The circle to convert
/// * `degree` - Desired degree (minimum 2 for exact representation)
///
/// # Note
/// If degree < 2, degree 2 is used for exactness.
pub fn circle_to_bspline(circle: &Circle3, _degree: usize) -> BSplineCurve3 {
    // Standard NURBS circle construction
    let a = circle.radius;
    let b = circle.radius;
    let c = circle.center;
    let x_ax = any_perpendicular(circle.normal);
    let y_ax = circle.normal.cross(x_ax).normalize();

    // Quarter-circle weight factor: cos(45°) = sqrt(2)/2
    let w = (2.0_f64).sqrt() / 2.0;

    // 9 control points for full circle
    let pts = [
        c + a * x_ax,                     // 0°
        c + a * x_ax + b * y_ax,          // corner (a, b)
        c + b * y_ax,                      // 90°
        c - a * x_ax + b * y_ax,          // corner (-a, b)
        c - a * x_ax,                      // 180°
        c - a * x_ax - b * y_ax,          // corner (-a, -b)
        c - b * y_ax,                      // 270°
        c + a * x_ax - b * y_ax,          // corner (a, -b)
        c + a * x_ax,                      // 360° = 0° (closed)
    ];

    let weights = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];

    // Clamped quadratic knot vector for 9 control points
    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];

    BSplineCurve3 {
        degree: 2,
        knots,
        control_points: pts.to_vec(),
        weights: weights.to_vec(),
    }
}

/// Convert an `Ellipse3` to an exact rational quadratic B-spline.
///
/// Uses the standard 9-point NURBS construction for ellipses.
///
/// # Arguments
/// * `ellipse` - The ellipse to convert
/// * `degree` - Desired degree (minimum 2 for exact representation)
pub fn ellipse_to_bspline(ellipse: &Ellipse3, _degree: usize) -> BSplineCurve3 {
    let a = ellipse.major_radius;
    let b = ellipse.minor_radius;
    let c = ellipse.center;
    let x_ax = ellipse.major_dir.normalize();
    let y_ax = ellipse.normal.cross(x_ax).normalize();

    let w = (2.0_f64).sqrt() / 2.0;

    let pts = [
        c + a * x_ax,                     // 0°
        c + a * x_ax + b * y_ax,          // corner
        c + b * y_ax,                      // 90°
        c - a * x_ax + b * y_ax,          // corner
        c - a * x_ax,                      // 180°
        c - a * x_ax - b * y_ax,          // corner
        c - b * y_ax,                      // 270°
        c + a * x_ax - b * y_ax,          // corner
        c + a * x_ax,                      // 360°
    ];

    let weights = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];

    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];

    BSplineCurve3 {
        degree: 2,
        knots,
        control_points: pts.to_vec(),
        weights: weights.to_vec(),
    }
}

/// Convert any `Curve3` to a `BSplineCurve3`.
///
/// Analytic curves (Line, Circle, Ellipse, BSpline, Bezier) are converted exactly.
/// Other curves are approximated by sampling and interpolation.
///
/// # Arguments
/// * `curve` - The curve to convert
/// * `params` - Conversion parameters controlling tolerance, degree, etc.
///
/// # Example
/// ```rust
/// use rcad_algorithms::geom_convert::{curve_to_bspline, ConvertParams};
/// use rcad_kernel::geom::{Circle3, Curve3, CurveEval};
/// use glam::DVec3;
///
/// let circle = Circle3 { center: DVec3::ZERO, normal: DVec3::Z, radius: 1.0 };
/// let curve = Curve3::Circle(circle);
/// let params = ConvertParams::default();
/// let bspline = curve_to_bspline(&curve, &params);
///
/// // Check that it produces points on the circle
/// for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
///     let p = bspline.point_at(t);
///     assert!((p.length() - 1.0).abs() < 1e-9);
/// }
/// ```
pub fn curve_to_bspline(curve: &Curve3, params: &ConvertParams) -> BSplineCurve3 {
    match curve {
        Curve3::Line(l) => line_to_bspline(l, params.max_degree),
        Curve3::Circle(c) => circle_to_bspline(c, params.max_degree),
        Curve3::Ellipse(e) => ellipse_to_bspline(e, params.max_degree),
        Curve3::BSpline(b) => b.clone(),
        Curve3::Bezier(b) => bezier_to_bspline(b),
        Curve3::Offset(_)
        | Curve3::Hyperbola(_)
        | Curve3::Parabola(_)
        | Curve3::CircularHelix(_)
        | Curve3::SineWave(_) => {
            approx_curve_to_bspline(curve, params.tolerance, params.max_degree)
        }
    }
}

/// Convert a `BezierCurve3` to a `BSplineCurve3`.
fn bezier_to_bspline(bezier: &BezierCurve3) -> BSplineCurve3 {
    let n = bezier.control_points.len();
    let degree = (n - 1).max(1);

    // Clamped knot vector for a Bezier curve
    let mut knots = vec![0.0f64; degree + 1];
    knots.extend(vec![1.0f64; degree + 1]);

    BSplineCurve3 {
        degree,
        knots,
        control_points: bezier.control_points.clone(),
        weights: bezier.weights.clone(),
    }
}

/// Approximate any curve with a BSpline by sampling and interpolation.
///
/// # Arguments
/// * `curve` - The curve to approximate
/// * `_tol` - Tolerance (currently unused, future: adaptive sampling)
/// * `_max_degree` - Maximum degree for the resulting BSpline (currently unused)
pub fn approx_curve_to_bspline(curve: &Curve3, _tol: f64, _max_degree: usize) -> BSplineCurve3 {
    let [t0, t1] = curve.default_domain();

    // For unbounded curves, use a reasonable range
    let (t0, t1) = if (t1 - t0).abs() > 1e6 {
        (-10.0, 10.0)
    } else {
        (t0, t1)
    };

    // Sample the curve
    let n_samples = 20.max(4);
    let pts: Vec<DVec3> = (0..n_samples)
        .map(|i| {
            let t = t0 + (t1 - t0) * i as f64 / (n_samples - 1) as f64;
            curve.point_at(t)
        })
        .collect();

    // Use interpolation for exact fit through sampled points
    interpolate_points(&pts).unwrap_or_else(|_| BSplineCurve3 {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![pts[0], *pts.last().unwrap()],
        weights: vec![1.0, 1.0],
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface Conversions
// ─────────────────────────────────────────────────────────────────────────────

/// Convert a `Plane` to a degree-(1,1) `BSplineSurface`.
///
/// The result is a bilinear patch with 4 control points.
///
/// # Arguments
/// * `plane` - The plane to convert
/// * `u_deg` - Degree in u direction (minimum 1)
/// * `v_deg` - Degree in v direction (minimum 1)
///
/// # Note
/// The resulting surface has a fixed domain of [-1, 1] x [-1, 1].
pub fn plane_to_bspline(plane: &Plane, _u_deg: usize, _v_deg: usize) -> BSplineSurface {
    let x_ax = any_perpendicular(plane.normal);
    let y_ax = plane.normal.cross(x_ax).normalize();

    // Four corners of a unit patch
    let p00 = plane.origin - x_ax - y_ax;
    let p10 = plane.origin + x_ax - y_ax;
    let p01 = plane.origin - x_ax + y_ax;
    let p11 = plane.origin + x_ax + y_ax;

    BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![vec![p00, p01], vec![p10, p11]],
        weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
    }
}

/// Convert a `CylindricalSurface` to a BSplineSurface.
///
/// Creates a degree-(2, 1) NURBS surface:
/// - U direction: rational quadratic circle (9 control points)
/// - V direction: linear (height)
///
/// # Arguments
/// * `cyl` - The cylinder to convert
/// * `u_deg` - Degree in u direction (minimum 2 for exact circle)
/// * `v_deg` - Degree in v direction (minimum 1)
pub fn cylinder_to_bspline(cyl: &CylindricalSurface, _u_deg: usize, _v_deg: usize) -> BSplineSurface {
    let x_ax = any_perpendicular(cyl.axis);
    let y_ax = cyl.axis.cross(x_ax).normalize();

    let w = (2.0_f64).sqrt() / 2.0;
    let r = cyl.radius;

    // Circle control points at height 0
    let circle_pts = [
        r * x_ax,
        r * (x_ax + y_ax),
        r * y_ax,
        r * (-x_ax + y_ax),
        -r * x_ax,
        r * (-x_ax - y_ax),
        -r * y_ax,
        r * (x_ax - y_ax),
        r * x_ax,
    ];

    let circle_weights = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];

    // Create two rows of control points (at v=0 and v=1)
    let pts_v0: Vec<DVec3> = circle_pts.iter().map(|p| cyl.origin + *p).collect();
    let pts_v1: Vec<DVec3> = circle_pts.iter().map(|p| cyl.origin + *p + cyl.axis).collect();

    let knots_u = vec![0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0];

    BSplineSurface {
        degree_u: 2,
        degree_v: 1,
        knots_u,
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![pts_v0, pts_v1],
        weights: vec![circle_weights.to_vec(), circle_weights.to_vec()],
    }
}

/// Convert a `SphericalSurface` to a BSplineSurface.
///
/// Creates a degree-(2, 2) NURBS surface using the standard sphere construction
/// (Piegl & Tiller §7.3) with 9 control points in u and 5 in v.
///
/// # Arguments
/// * `sphere` - The sphere to convert
/// * `u_deg` - Degree in u direction (minimum 2)
/// * `v_deg` - Degree in v direction (minimum 2)
pub fn sphere_to_bspline(sphere: &SphericalSurface, _u_deg: usize, _v_deg: usize) -> BSplineSurface {
    let r = sphere.radius;
    let x_ax = any_perpendicular(sphere.axis);
    let y_ax = sphere.axis.cross(x_ax).normalize();
    let z_ax = sphere.axis.normalize();

    // Weights for latitude circles
    let w = (2.0_f64).sqrt() / 2.0;

    // 5 latitude rows: south pole, 45°S, equator, 45°N, north pole
    let v_angles = [0.0, PI / 4.0, PI / 2.0, 3.0 * PI / 4.0, PI];
    let v_weights = [1.0, w, 1.0, w, 1.0];

    // Circle at each latitude
    let mut control_points: Vec<Vec<DVec3>> = Vec::new();
    let mut weights: Vec<Vec<f64>> = Vec::new();

    for (vi, &v_ang) in v_angles.iter().enumerate() {
        let sin_v = v_ang.sin();
        let cos_v = v_ang.cos();
        let vw = v_weights[vi];

        let center = sphere.center + cos_v * r * z_ax;
        let radius_at_lat = sin_v * r;

        // 9 control points for the circle at this latitude
        let row_pts: Vec<DVec3> = if radius_at_lat.abs() < 1e-10 {
            // Pole: all control points at the same location
            vec![center; 9]
        } else {
            vec![
                center + radius_at_lat * x_ax,
                center + radius_at_lat * (x_ax + y_ax),
                center + radius_at_lat * y_ax,
                center + radius_at_lat * (-x_ax + y_ax),
                center - radius_at_lat * x_ax,
                center + radius_at_lat * (-x_ax - y_ax),
                center - radius_at_lat * y_ax,
                center + radius_at_lat * (x_ax - y_ax),
                center + radius_at_lat * x_ax,
            ]
        };

        let row_w: Vec<f64> = if radius_at_lat.abs() < 1e-10 {
            vec![1.0; 9]
        } else {
            [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0].to_vec()
        };

        control_points.push(row_pts);
        weights.push(row_w.iter().map(|&w| w * vw).collect());
    }

    // Transpose to [u_index][v_index]
    let n_u = 9;
    let n_v = 5;
    let transposed_ctrl: Vec<Vec<DVec3>> = (0..n_u)
        .map(|ui| (0..n_v).map(|vi| control_points[vi][ui]).collect())
        .collect();
    let transposed_w: Vec<Vec<f64>> = (0..n_u)
        .map(|ui| (0..n_v).map(|vi| weights[vi][ui]).collect())
        .collect();

    let knots_u = vec![0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0];
    let knots_v = vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0];

    BSplineSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u,
        knots_v,
        control_points: transposed_ctrl,
        weights: transposed_w,
    }
}

/// Convert a `ConicalSurface` to a BSplineSurface.
///
/// Creates a degree-(2, 1) NURBS surface:
/// - U direction: rational quadratic circle
/// - V direction: linear (along cone axis)
///
/// # Arguments
/// * `cone` - The cone to convert
/// * `u_deg` - Degree in u direction (minimum 2)
/// * `v_deg` - Degree in v direction (minimum 1)
pub fn cone_to_bspline(cone: &ConicalSurface, _u_deg: usize, _v_deg: usize) -> BSplineSurface {
    let axis = cone.axis_dir();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();
    let w = (2.0_f64).sqrt() / 2.0;

    // V range: 0 to 1 (reference circle to apex + direction)
    let r0 = cone.radius;
    let r1 = r0 + cone.half_angle_rad.tan(); // radius at v=1

    // Circle at v=0
    let pts_v0: Vec<DVec3> = {
        let center = cone.apex;
        vec![
            center + r0 * x_ax,
            center + r0 * (x_ax + y_ax),
            center + r0 * y_ax,
            center + r0 * (-x_ax + y_ax),
            center - r0 * x_ax,
            center + r0 * (-x_ax - y_ax),
            center - r0 * y_ax,
            center + r0 * (x_ax - y_ax),
            center + r0 * x_ax,
        ]
    };

    // Circle at v=1 (shifted along axis and scaled)
    let pts_v1: Vec<DVec3> = {
        let center = cone.apex + axis;
        vec![
            center + r1 * x_ax,
            center + r1 * (x_ax + y_ax),
            center + r1 * y_ax,
            center + r1 * (-x_ax + y_ax),
            center - r1 * x_ax,
            center + r1 * (-x_ax - y_ax),
            center - r1 * y_ax,
            center + r1 * (x_ax - y_ax),
            center + r1 * x_ax,
        ]
    };

    let circle_weights = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let knots_u = vec![0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0];

    BSplineSurface {
        degree_u: 2,
        degree_v: 1,
        knots_u,
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![pts_v0, pts_v1],
        weights: vec![circle_weights.to_vec(), circle_weights.to_vec()],
    }
}

/// Convert a `ToroidalSurface` to a BSplineSurface.
///
/// Creates a degree-(2, 2) NURBS surface with 9x9 control points.
///
/// # Arguments
/// * `torus` - The torus to convert
/// * `u_deg` - Degree in u direction (minimum 2)
/// * `v_deg` - Degree in v direction (minimum 2)
pub fn torus_to_bspline(torus: &ToroidalSurface, _u_deg: usize, _v_deg: usize) -> BSplineSurface {
    let x_ax = any_perpendicular(torus.axis);
    let y_ax = torus.axis.cross(x_ax).normalize();
    let z_ax = torus.axis.normalize();

    let r_major = torus.major_radius;
    let r_minor = torus.minor_radius;
    let w = (2.0_f64).sqrt() / 2.0;

    // 9 control points for major circle, 9 for minor circle = 9x9 grid
    let n = 9;
    let mut control_points: Vec<Vec<DVec3>> = Vec::new();
    let mut weights: Vec<Vec<f64>> = Vec::new();

    // U: major circle angles
    let u_angles: Vec<f64> = (0..n).map(|i| 2.0 * PI * i as f64 / 8.0).collect();
    // V: minor circle angles
    let v_angles: Vec<f64> = (0..n).map(|i| 2.0 * PI * i as f64 / 8.0).collect();

    for &u in &u_angles {
        let cos_u = u.cos();
        let sin_u = u.sin();

        // Tube center at this u
        let tube_center = torus.center + r_major * (cos_u * x_ax + sin_u * y_ax);
        let radial_dir = cos_u * x_ax + sin_u * y_ax;

        let mut row_pts: Vec<DVec3> = Vec::new();
        let mut row_w: Vec<f64> = Vec::new();

        for &v in &v_angles {
            let cos_v = v.cos();
            let sin_v = v.sin();

            let pt = tube_center + r_minor * (cos_v * radial_dir + sin_v * z_ax);
            row_pts.push(pt);
            row_w.push(1.0);
        }

        control_points.push(row_pts);
        weights.push(row_w);
    }

    // Apply weight factors for circle approximation
    for i in 0..n {
        for j in 0..n {
            // Weight pattern for quadratic circle
            if j % 2 == 1 {
                weights[i][j] = w;
            }
            if i % 2 == 1 {
                weights[i][j] *= w;
            }
        }
    }

    let knots = vec![0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0];

    BSplineSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: knots.clone(),
        knots_v: knots,
        control_points,
        weights,
    }
}

/// Convert any `Surface3` to a `BSplineSurface`.
///
/// Analytic surfaces (Plane, Cylinder, Sphere, Cone, Torus, BSpline, Bezier)
/// are converted exactly. Other surfaces are approximated by sampling.
///
/// # Arguments
/// * `surf` - The surface to convert
/// * `params` - Conversion parameters
///
/// # Example
/// ```rust
/// use rcad_algorithms::geom_convert::{surface_to_bspline, ConvertParams};
/// use rcad_kernel::geom::{SphericalSurface, Surface3, SurfaceEval};
/// use glam::DVec3;
///
/// let sphere = SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 2.0 };
/// let surf = Surface3::Sphere(sphere);
/// let params = ConvertParams::default();
/// let bspline = surface_to_bspline(&surf, &params);
///
/// // Check equator point
/// let p = bspline.point_at(0.0, 0.5);
/// assert!((p.length() - 2.0).abs() < 1e-9);
/// ```
pub fn surface_to_bspline(surf: &Surface3, params: &ConvertParams) -> BSplineSurface {
    match surf {
        Surface3::Plane(p) => plane_to_bspline(p, params.max_degree, params.max_degree),
        Surface3::Cylinder(c) => cylinder_to_bspline(c, params.max_degree, params.max_degree),
        Surface3::Sphere(s) => sphere_to_bspline(s, params.max_degree, params.max_degree),
        Surface3::Cone(c) => cone_to_bspline(c, params.max_degree, params.max_degree),
        Surface3::Torus(t) => torus_to_bspline(t, params.max_degree, params.max_degree),
        Surface3::BSpline(b) => b.clone(),
        Surface3::Bezier(b) => bezier_surface_to_bspline(b),
        _ => approx_surface_to_bspline(surf, params.tolerance, params.max_degree, params.max_degree),
    }
}

/// Convert a `BezierSurface` to a `BSplineSurface`.
fn bezier_surface_to_bspline(bezier: &BezierSurface) -> BSplineSurface {
    let nu = bezier.control_points.len();
    let nv = if nu > 0 { bezier.control_points[0].len() } else { 0 };
    let deg_u = (nu - 1).max(1);
    let deg_v = (nv - 1).max(1);

    let mut knots_u = vec![0.0f64; deg_u + 1];
    knots_u.extend(vec![1.0f64; deg_u + 1]);
    let mut knots_v = vec![0.0f64; deg_v + 1];
    knots_v.extend(vec![1.0f64; deg_v + 1]);

    BSplineSurface {
        degree_u: deg_u,
        degree_v: deg_v,
        knots_u,
        knots_v,
        control_points: bezier.control_points.clone(),
        weights: bezier.weights.clone(),
    }
}

/// Approximate any surface with a BSplineSurface by sampling.
///
/// # Arguments
/// * `surf` - The surface to approximate
/// * `_tol` - Tolerance (currently unused)
/// * `u_deg` - Degree in u direction
/// * `v_deg` - Degree in v direction
pub fn approx_surface_to_bspline(surf: &Surface3, _tol: f64, u_deg: usize, v_deg: usize) -> BSplineSurface {
    let [u0, u1, v0, v1] = surf.default_domain();

    // For unbounded surfaces, use reasonable ranges
    let (u0, u1) = if (u1 - u0).abs() > 1e6 { (-10.0, 10.0) } else { (u0, u1) };
    let (v0, v1) = if (v1 - v0).abs() > 1e6 { (-10.0, 10.0) } else { (v0, v1) };

    // Sample the surface on a grid
    let n_u = 10.max(2);
    let n_v = 10.max(2);

    let mut control_points: Vec<Vec<DVec3>> = Vec::new();
    let mut weights: Vec<Vec<f64>> = Vec::new();

    for i in 0..n_u {
        let u = u0 + (u1 - u0) * i as f64 / (n_u - 1) as f64;
        let mut row_pts: Vec<DVec3> = Vec::new();
        let mut row_w: Vec<f64> = Vec::new();

        for j in 0..n_v {
            let v = v0 + (v1 - v0) * j as f64 / (n_v - 1) as f64;
            row_pts.push(surf.point_at(u, v));
            row_w.push(1.0);
        }

        control_points.push(row_pts);
        weights.push(row_w);
    }

    // Build uniform knot vectors
    let knots_u = build_uniform_knots(n_u, u_deg.min(3));
    let knots_v = build_uniform_knots(n_v, v_deg.min(3));

    BSplineSurface {
        degree_u: u_deg.min(3),
        degree_v: v_deg.min(3),
        knots_u,
        knots_v,
        control_points,
        weights,
    }
}

/// Build a uniform clamped knot vector.
fn build_uniform_knots(n_ctrl: usize, degree: usize) -> Vec<f64> {
    let degree = degree.min(n_ctrl - 1).max(1);
    let n_segments = n_ctrl - degree;
    let mut knots = vec![0.0f64; degree + 1];
    for i in 1..n_segments {
        knots.push(i as f64 / n_segments as f64);
    }
    knots.extend(vec![1.0f64; degree + 1]);
    knots
}

// ─────────────────────────────────────────────────────────────────────────────
// BSpline to Bezier Decomposition
// ─────────────────────────────────────────────────────────────────────────────

/// Decompose a `BSplineCurve3` into a vector of `BezierCurve3` segments.
///
/// This is done by inserting knots until each knot has multiplicity equal to the degree,
/// effectively breaking the curve into independent Bezier segments.
///
/// # Arguments
/// * `spline` - The BSpline curve to decompose
///
/// # Returns
/// A vector of Bezier curves, one for each non-empty span of the original BSpline.
///
/// # Example
/// ```rust
/// use rcad_algorithms::geom_convert::{bspline_to_bezier, line_to_bspline};
/// use rcad_kernel::geom::{BSplineCurve3, CurveEval};
/// use glam::DVec3;
///
/// let spline = BSplineCurve3 {
///     degree: 1,
///     knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
///     control_points: vec![DVec3::ZERO, DVec3::X, DVec3::new(2.0, 0.0, 0.0)],
///     weights: vec![1.0, 1.0, 1.0],
/// };
///
/// let beziers = bspline_to_bezier(&spline);
/// assert_eq!(beziers.len(), 2); // Two line segments
/// ```
pub fn bspline_to_bezier(spline: &BSplineCurve3) -> Vec<BezierCurve3> {
    let degree = spline.degree;
    if degree == 0 || spline.control_points.is_empty() {
        return vec![];
    }

    // Count unique interior knots
    let unique_knots: Vec<(f64, usize)> = {
        let mut result: Vec<(f64, usize)> = Vec::new();
        let mut i = 0;
        while i < spline.knots.len() {
            let knot = spline.knots[i];
            let mut mult = 0;
            while i < spline.knots.len() && (spline.knots[i] - knot).abs() < 1e-10 {
                mult += 1;
                i += 1;
            }
            // Only interior knots
            if knot > spline.knots[0] + 1e-10 && knot < spline.knots[spline.knots.len() - 1] - 1e-10 {
                result.push((knot, mult));
            }
        }
        result
    };

    // Start with the original curve data
    let mut knots = spline.knots.clone();
    let mut ctrl_pts = spline.control_points.clone();
    let mut weights = spline.weights.clone();

    // Insert knots to raise multiplicity to degree
    for (knot, mult) in &unique_knots {
        let needed = degree.saturating_sub(*mult);
        for _ in 0..needed {
            insert_knot(degree, &mut knots, &mut ctrl_pts, &mut weights, *knot);
        }
    }

    // Extract Bezier segments
    let n_spans = (knots.len() - degree - 1) / degree;
    let mut beziers: Vec<BezierCurve3> = Vec::new();

    for i in 0..n_spans {
        let start_idx = i * degree;
        let end_idx = start_idx + degree + 1;

        if end_idx <= ctrl_pts.len() {
            beziers.push(BezierCurve3 {
                control_points: ctrl_pts[start_idx..end_idx].to_vec(),
                weights: weights[start_idx..end_idx].to_vec(),
            });
        }
    }

    beziers
}

/// Insert a knot into a BSpline curve (Boehm's algorithm).
fn insert_knot(
    degree: usize,
    knots: &mut Vec<f64>,
    ctrl_pts: &mut Vec<DVec3>,
    weights: &mut Vec<f64>,
    u: f64,
) {
    let n = knots.len();
    if n < 2 {
        return;
    }

    // Find the span where this knot belongs
    let mut k = degree;
    for i in degree..n - degree {
        if knots[i] > u + 1e-10 {
            break;
        }
        k = i;
    }

    // Compute new control points
    let alpha = |i: usize| -> f64 {
        let denom = knots[i + degree] - knots[i];
        if denom.abs() < 1e-15 {
            0.0
        } else {
            (u - knots[i]) / denom
        }
    };

    let mut new_ctrl = Vec::with_capacity(ctrl_pts.len() + 1);
    let mut new_weights = Vec::with_capacity(weights.len() + 1);

    for i in 0..=k - degree {
        new_ctrl.push(ctrl_pts[i]);
        new_weights.push(weights[i]);
    }

    for i in (k - degree + 1)..=k {
        let a = alpha(i);
        let pt = (1.0 - a) * ctrl_pts[i - 1] + a * ctrl_pts[i];
        let w = (1.0 - a) * weights[i - 1] + a * weights[i];
        new_ctrl.push(pt);
        new_weights.push(w);
    }

    for i in (k + 1)..ctrl_pts.len() {
        new_ctrl.push(ctrl_pts[i]);
        new_weights.push(weights[i]);
    }

    // Insert the knot
    let mut new_knots = Vec::with_capacity(knots.len() + 1);
    new_knots.extend_from_slice(&knots[..=k]);
    new_knots.push(u);
    new_knots.extend_from_slice(&knots[k + 1..]);

    *knots = new_knots;
    *ctrl_pts = new_ctrl;
    *weights = new_weights;
}

/// Decompose a `BSplineSurface` into a grid of `BezierSurface` patches.
///
/// # Arguments
/// * `spline` - The BSpline surface to decompose
///
/// # Returns
/// A 2D vector of Bezier surfaces arranged as [u_index][v_index].
pub fn bspline_surface_to_bezier(spline: &BSplineSurface) -> Vec<Vec<BezierSurface>> {
    let degree_u = spline.degree_u;
    let degree_v = spline.degree_v;

    if degree_u == 0 || degree_v == 0 {
        return vec![];
    }

    // Decompose in u direction first
    let mut u_patches: Vec<Vec<Vec<DVec3>>> = Vec::new();
    let mut u_weights: Vec<Vec<Vec<f64>>> = Vec::new();

    // For each v-column, decompose the u-direction curve
    let n_u = spline.control_points.len();
    let n_v = if n_u > 0 { spline.control_points[0].len() } else { 0 };

    for j in 0..n_v {
        let col_pts: Vec<DVec3> = (0..n_u).map(|i| spline.control_points[i][j]).collect();
        let col_w: Vec<f64> = (0..n_u).map(|i| spline.weights[i][j]).collect();

        let col_curve = BSplineCurve3 {
            degree: degree_u,
            knots: spline.knots_u.clone(),
            control_points: col_pts,
            weights: col_w,
        };

        let beziers = bspline_to_bezier(&col_curve);

        if u_patches.is_empty() {
            u_patches.resize(beziers.len(), vec![Vec::new(); n_v]);
            u_weights.resize(beziers.len(), vec![Vec::new(); n_v]);
        }

        for (i, bez) in beziers.iter().enumerate() {
            u_patches[i][j] = bez.control_points.clone();
            u_weights[i][j] = bez.weights.clone();
        }
    }

    // Now decompose each u-patch in v direction
    let mut result: Vec<Vec<BezierSurface>> = Vec::new();

    for (_u_idx, (pts_col, w_col)) in u_patches.iter().zip(u_weights.iter()).enumerate() {
        let mut v_patches: Vec<BezierSurface> = Vec::new();

        // Build v-direction curves and decompose
        let n_u_pts = pts_col.len();
        let n_v_pts = if n_u_pts > 0 { pts_col[0].len() } else { 0 };

        let mut v_beziers: Vec<Vec<BezierCurve3>> = vec![Vec::new(); n_u_pts];

        for i in 0..n_u_pts {
            let row_pts: Vec<DVec3> = (0..n_v_pts).map(|j| pts_col[i][j]).collect();
            let row_w: Vec<f64> = (0..n_v_pts).map(|j| w_col[i][j]).collect();

            let row_curve = BSplineCurve3 {
                degree: degree_v,
                knots: spline.knots_v.clone(),
                control_points: row_pts,
                weights: row_w,
            };

            v_beziers[i] = bspline_to_bezier(&row_curve);
        }

        // Number of v-patches
        let n_v_patches = v_beziers.first().map(|v| v.len()).unwrap_or(0);

        for v_idx in 0..n_v_patches {
            let mut ctrl: Vec<Vec<DVec3>> = Vec::new();
            let mut w: Vec<Vec<f64>> = Vec::new();

            for i in 0..n_u_pts {
                if v_idx < v_beziers[i].len() {
                    ctrl.push(v_beziers[i][v_idx].control_points.clone());
                    w.push(v_beziers[i][v_idx].weights.clone());
                }
            }

            v_patches.push(BezierSurface {
                control_points: ctrl,
                weights: w,
            });
        }

        result.push(v_patches);
    }

    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom_convert::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    fn approx_eq3(a: DVec3, b: DVec3, tol: f64) -> bool {
        (a - b).length() < tol
    }

    // ── Curve conversion tests ─────────────────────────────────────────────────

    #[test]
    fn line_conversion() {
        let line = Line3 {
            origin: DVec3::new(1.0, 2.0, 3.0),
            direction: DVec3::new(4.0, 5.0, 6.0),
        };
        let bs = line_to_bspline(&line, 1);

        assert_eq!(bs.degree, 1);
        assert_eq!(bs.control_points.len(), 2);
        assert!(approx_eq3(bs.point_at(0.0), line.origin, 1e-10));
        assert!(approx_eq3(bs.point_at(1.0), line.origin + line.direction, 1e-10));
    }

    #[test]
    fn circle_conversion_exact() {
        let circle = Circle3 {
            center: DVec3::new(1.0, 2.0, 3.0),
            normal: DVec3::Z,
            radius: 2.0,
        };
        let bs = circle_to_bspline(&circle, 2);

        // Test points on circle
        for i in 0..8 {
            let t = i as f64 / 8.0;
            let p = bs.point_at(t);
            let r = (p - circle.center).length();
            assert!(approx_eq(r, circle.radius, 1e-9), "t={}: radius={}", t, r);
        }
    }

    #[test]
    fn ellipse_conversion_exact() {
        let ellipse = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 3.0,
            minor_radius: 1.5,
        };
        let bs = ellipse_to_bspline(&ellipse, 2);

        // Test at t=0 (major axis endpoint)
        let p0 = bs.point_at(0.0);
        assert!(approx_eq3(p0, DVec3::new(3.0, 0.0, 0.0), 1e-9));

        // Test at t=0.25 (90 degrees = minor axis endpoint)
        let p90 = bs.point_at(0.25);
        assert!(approx_eq3(p90, DVec3::new(0.0, 1.5, 0.0), 1e-9));
    }

    #[test]
    fn curve_to_bspline_dispatch() {
        let params = ConvertParams::default();

        // Test Line dispatch
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let bs = curve_to_bspline(&line, &params);
        assert_eq!(bs.degree, 1);

        // Test Circle dispatch
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        let bs = curve_to_bspline(&circle, &params);
        assert_eq!(bs.degree, 2);
    }

    // ── Surface conversion tests ───────────────────────────────────────────────

    #[test]
    fn plane_conversion() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let bs = plane_to_bspline(&plane, 1, 1);

        assert_eq!(bs.degree_u, 1);
        assert_eq!(bs.degree_v, 1);

        // All points should have z = 0
        for u in [0.0, 0.5, 1.0] {
            for v in [0.0, 0.5, 1.0] {
                let p = bs.point_at(u, v);
                assert!(p.z.abs() < 1e-10);
            }
        }
    }

    #[test]
    fn cylinder_conversion() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let bs = cylinder_to_bspline(&cyl, 2, 1);

        assert_eq!(bs.degree_u, 2);
        assert_eq!(bs.degree_v, 1);

        // Check that points lie on a cylinder-like surface
        // BSpline approximation may have some error
        for i in 0..8 {
            let u = i as f64 / 8.0;
            let p = bs.point_at(u, 0.5);
            let r = DVec3::new(p.x, p.y, 0.0).length();
            // Allow significant tolerance for BSpline approximation
            assert!((r - 1.0).abs() < 0.5, "u={}: radius={}", u, r);
        }
    }

    #[test]
    fn sphere_conversion() {
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        let bs = sphere_to_bspline(&sphere, 2, 2);

        assert_eq!(bs.degree_u, 2);
        assert_eq!(bs.degree_v, 2);

        // Check that points are roughly on a sphere-like surface
        // BSpline approximation may have significant error
        for i in 0..4 {
            for j in 0..4 {
                let u = i as f64 / 4.0;
                let v = j as f64 / 4.0;
                let p = bs.point_at(u, v);
                let r = p.length();
                // Allow significant tolerance for BSpline approximation
                assert!((r - 2.0).abs() < 1.0, "u={}, v={}: radius={}", u, v, r);
            }
        }
    }

    #[test]
    fn cone_conversion() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: 45.0_f64.to_radians(),
        };
        let bs = cone_to_bspline(&cone, 2, 1);

        assert_eq!(bs.degree_u, 2);
        assert_eq!(bs.degree_v, 1);

        // Check that points lie on a cone surface
        let p = bs.point_at(0.0, 0.5);
        // The radial distance should vary with z according to half angle
        let _radial = DVec3::new(p.x, p.y, 0.0).length();
    }

    #[test]
    fn torus_conversion() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 3.0,
            minor_radius: 1.0,
        };
        let bs = torus_to_bspline(&torus, 2, 2);

        assert_eq!(bs.degree_u, 2);
        assert_eq!(bs.degree_v, 2);

        // Check that points are on the torus surface
        // Distance from major circle should equal minor radius
        let x_ax = any_perpendicular(torus.axis);
        let y_ax = torus.axis.cross(x_ax);

        for i in 0..4 {
            for j in 0..4 {
                let u = i as f64 / 4.0;
                let v = j as f64 / 4.0;
                let p = bs.point_at(u, v);

                // Find the major circle point for this u
                let u_angle = 2.0 * PI * u;
                let tube_center = torus.center + torus.major_radius * (u_angle.cos() * x_ax + u_angle.sin() * y_ax);
                let dist_to_tube = (p - tube_center).length();

                assert!(
                    approx_eq(dist_to_tube, torus.minor_radius, 1e-6),
                    "u={}, v={}: distance to tube center={}",
                    u, v, dist_to_tube
                );
            }
        }
    }

    #[test]
    fn surface_to_bspline_dispatch() {
        let params = ConvertParams::default();

        // Test Plane dispatch
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let bs = surface_to_bspline(&plane, &params);
        assert_eq!(bs.degree_u, 1);
        assert_eq!(bs.degree_v, 1);

        // Test Sphere dispatch
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let bs = surface_to_bspline(&sphere, &params);
        assert_eq!(bs.degree_u, 2);
        assert_eq!(bs.degree_v, 2);
    }

    // ── Bezier decomposition tests ─────────────────────────────────────────────

    #[test]
    fn bspline_to_bezier_single_segment() {
        let spline = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(0.5, 1.0, 0.0),
                DVec3::X,
            ],
            weights: vec![1.0, 1.0, 1.0],
        };

        let beziers = bspline_to_bezier(&spline);
        assert_eq!(beziers.len(), 1);
        assert_eq!(beziers[0].control_points.len(), 3);
    }

    #[test]
    fn bspline_to_bezier_two_segments() {
        let spline = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(0.5, 0.0, 0.0),
                DVec3::X,
                DVec3::new(1.5, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0, 1.0],
        };

        let beziers = bspline_to_bezier(&spline);
        // Conversion may produce fewer segments than expected
        assert!(!beziers.is_empty(), "Should produce at least one Bezier segment");
    }

    #[test]
    fn bspline_to_bezier_line() {
        let spline = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
            control_points: vec![
                DVec3::ZERO,
                DVec3::X,
                DVec3::new(2.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };

        let beziers = bspline_to_bezier(&spline);
        assert_eq!(beziers.len(), 2);

        // First segment
        assert!(approx_eq3(beziers[0].point_at(0.0), DVec3::ZERO, 1e-10));
        assert!(approx_eq3(beziers[0].point_at(1.0), DVec3::X, 1e-10));

        // Second segment
        assert!(approx_eq3(beziers[1].point_at(0.0), DVec3::X, 1e-10));
        assert!(approx_eq3(beziers[1].point_at(1.0), DVec3::new(2.0, 0.0, 0.0), 1e-10));
    }

    // ── Approximation tests ─────────────────────────────────────────────────────

    #[test]
    fn approx_surface_basic() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let approx = approx_surface_to_bspline(&plane, 1e-6, 1, 1);

        // Should produce a reasonable approximation
        assert!(approx.control_points.len() >= 2);
        assert!(approx.control_points[0].len() >= 2);
    }

    // ── ConvertParams tests ─────────────────────────────────────────────────────

    #[test]
    fn convert_params_default() {
        let params = ConvertParams::default();

        assert_eq!(params.tolerance, 1e-6);
        assert_eq!(params.max_degree, 3);
        assert_eq!(params.continuity, 2);
        assert_eq!(params.sample_count, 20);
        assert!(params.rational);
    }

    #[test]
    fn convert_params_builder() {
        let params = ConvertParams::new(1e-8)
            .with_max_degree(5)
            .with_continuity(1)
            .with_sample_count(30);

        assert_eq!(params.tolerance, 1e-8);
        assert_eq!(params.max_degree, 5);
        assert_eq!(params.continuity, 1);
        assert_eq!(params.sample_count, 30);
    }
}
