//! GCPnts-style point sampling algorithms.
//!
//! Analogous to OCCT `GCPnts` package, provides algorithms for sampling points
//! on curves and surfaces:
//! - `AbscissaPoint`: Point at arc length distance
//! - `UniformAbscissa`: Uniform arc length sampling
//! - `UniformDeflection`: Sample by deviation tolerance
//! - `TangentialDeflection`: Sample by tangent deviation
//! - `QuasiUniformAbscissa`: Quasi-uniform sampling
//! - `SurfaceSampling`: Sample on surfaces
//!
//! Uses adaptive subdivision and numerical integration for robust point placement.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};

use crate::tolerance::TOLERANCE_ABS;

/// Finite difference step size for derivative computation.
const H: f64 = 1e-6;

// =============================================================================
// Derivative Computation (Finite Differences)
// =============================================================================

/// Compute curve derivative via finite differences.
fn curve_derivative(curve: &Curve3, t: f64) -> DVec3 {
    (curve.point_at(t + H) - curve.point_at(t - H)) / (2.0 * H)
}

/// Compute curve second derivative via finite differences.
fn curve_second_derivative(curve: &Curve3, t: f64) -> DVec3 {
    let d_plus = curve_derivative(curve, t + H);
    let d_minus = curve_derivative(curve, t - H);
    (d_plus - d_minus) / (2.0 * H)
}

/// Compute surface partial derivatives via finite differences.
fn surface_derivatives(surface: &Surface3, u: f64, v: f64) -> (DVec3, DVec3) {
    let du = (surface.point_at(u + H, v) - surface.point_at(u - H, v)) / (2.0 * H);
    let dv = (surface.point_at(u, v + H) - surface.point_at(u, v - H)) / (2.0 * H);
    (du, dv)
}

// =============================================================================
// Arc Length Computation
// =============================================================================

/// Compute the arc length of a curve between two parameters.
///
/// Uses adaptive Simpson's rule for numerical integration.
pub fn arc_length(curve: &Curve3, t0: f64, t1: f64) -> f64 {
    let domain = curve_domain(curve);
    let t0_clamped = t0.clamp(domain[0], domain[1]);
    let t1_clamped = t1.clamp(domain[0], domain[1]);

    if (t1_clamped - t0_clamped).abs() < TOLERANCE_ABS {
        return 0.0;
    }

    // Adaptive Simpson's rule integration
    adaptive_simpson_arc_length(curve, t0_clamped, t1_clamped, 1e-8)
}

/// Compute the total arc length of a curve.
pub fn total_arc_length(curve: &Curve3) -> f64 {
    let domain = curve_domain(curve);
    arc_length(curve, domain[0], domain[1])
}

/// Adaptive Simpson's rule for arc length integration.
fn adaptive_simpson_arc_length(curve: &Curve3, a: f64, b: f64, tol: f64) -> f64 {
    fn speed(curve: &Curve3, t: f64) -> f64 {
        curve_derivative(curve, t).length()
    }

    let mid = (a + b) / 2.0;
    let fa = speed(curve, a);
    let fb = speed(curve, b);
    let fm = speed(curve, mid);

    let whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);

    // Left half
    let lm = (a + mid) / 2.0;
    let flm = speed(curve, lm);
    let left = (mid - a) / 6.0 * (fa + 4.0 * flm + fm);

    // Right half
    let rm = (mid + b) / 2.0;
    let frm = speed(curve, rm);
    let right = (b - mid) / 6.0 * (fm + 4.0 * frm + fb);

    let sum = left + right;

    if (sum - whole).abs() < 15.0 * tol {
        sum + (sum - whole) / 15.0
    } else {
        adaptive_simpson_arc_length(curve, a, mid, tol / 2.0)
            + adaptive_simpson_arc_length(curve, mid, b, tol / 2.0)
    }
}

// =============================================================================
// AbscissaPoint - Point at arc length
// =============================================================================

/// Find the point at a given arc length distance from the start of the curve.
///
/// Returns the point and its parameter value.
/// Uses binary search to find the parameter that achieves the desired arc length.
pub fn point_at_arc_length(curve: &Curve3, distance: f64) -> (DVec3, f64) {
    let domain = curve_domain(curve);
    let total_len = arc_length(curve, domain[0], domain[1]);

    // Clamp distance to valid range
    let clamped_dist = distance.clamp(0.0, total_len);

    if clamped_dist <= TOLERANCE_ABS {
        return (curve.point_at(domain[0]), domain[0]);
    }
    if (clamped_dist - total_len).abs() <= TOLERANCE_ABS {
        return (curve.point_at(domain[1]), domain[1]);
    }

    // Binary search for parameter
    let param = find_parameter_at_arc_length(curve, domain[0], domain[1], clamped_dist);
    (curve.point_at(param), param)
}

/// Find points at equal arc length intervals.
///
/// Returns a vector of (point, parameter) pairs.
pub fn points_at_equal_arc_length(curve: &Curve3, n_points: usize) -> Vec<(DVec3, f64)> {
    if n_points == 0 {
        return vec![];
    }

    let domain = curve_domain(curve);
    let total_len = arc_length(curve, domain[0], domain[1]);

    if total_len < TOLERANCE_ABS {
        // Degenerate curve, return start point only
        return vec![(curve.point_at(domain[0]), domain[0])];
    }

    let mut result = Vec::with_capacity(n_points);

    if n_points == 1 {
        let (pt, param) = point_at_arc_length(curve, total_len / 2.0);
        result.push((pt, param));
    } else {
        for i in 0..n_points {
            let dist = total_len * i as f64 / (n_points - 1) as f64;
            let (pt, param) = point_at_arc_length(curve, dist);
            result.push((pt, param));
        }
    }

    result
}

/// Binary search to find parameter at a given arc length from the start.
fn find_parameter_at_arc_length(curve: &Curve3, t_min: f64, t_max: f64, target_len: f64) -> f64 {
    let mut lo = t_min;
    let mut hi = t_max;
    let tol = TOLERANCE_ABS;

    for _ in 0..50 {
        // Max iterations for convergence
        let mid = (lo + hi) / 2.0;
        let len = arc_length(curve, t_min, mid);

        if (len - target_len).abs() < tol {
            return mid;
        }

        if len < target_len {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    (lo + hi) / 2.0
}

// =============================================================================
// UniformAbscissa - Uniform arc length sampling
// =============================================================================

/// Compute uniform arc length parameter distribution.
///
/// Returns parameter values for n_points uniformly distributed along arc length.
pub fn uniform_abscissa(curve: &Curve3, n_points: usize) -> Vec<f64> {
    if n_points == 0 {
        return vec![];
    }

    let domain = curve_domain(curve);
    let total_len = arc_length(curve, domain[0], domain[1]);

    if total_len < TOLERANCE_ABS || n_points == 1 {
        let mid = (domain[0] + domain[1]) / 2.0;
        return vec![mid];
    }

    let mut params = Vec::with_capacity(n_points);

    for i in 0..n_points {
        let dist = total_len * i as f64 / (n_points - 1) as f64;
        let param = find_parameter_at_arc_length(curve, domain[0], domain[1], dist);
        params.push(param);
    }

    params
}

/// Compute uniform arc length points on a curve.
///
/// Returns n_points uniformly distributed along arc length.
pub fn uniform_abscissa_points(curve: &Curve3, n_points: usize) -> Vec<DVec3> {
    uniform_abscissa(curve, n_points)
        .into_iter()
        .map(|t| curve.point_at(t))
        .collect()
}

// =============================================================================
// UniformDeflection - Sample by deviation
// =============================================================================

/// Sample a curve based on maximum deviation from chord.
///
/// Returns parameters where the curve deviates from straight lines by at most max_deviation.
/// Uses recursive subdivision.
pub fn uniform_deflection(curve: &Curve3, max_deviation: f64) -> Vec<f64> {
    let domain = curve_domain(curve);
    let mut params = vec![domain[0]];

    uniform_deflection_recursive(curve, domain[0], domain[1], max_deviation, &mut params);

    params.push(domain[1]); // Add endpoint
    params.sort_by(|a, b| a.partial_cmp(b).unwrap());
    params.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE_ABS);
    params
}

/// Recursive helper for uniform_deflection.
fn uniform_deflection_recursive(
    curve: &Curve3,
    t0: f64,
    t1: f64,
    max_dev: f64,
    params: &mut Vec<f64>,
) {
    let p0 = curve.point_at(t0);
    let p1 = curve.point_at(t1);

    // Find the point of maximum deviation from the chord
    let mid = (t0 + t1) / 2.0;
    let pmid = curve.point_at(mid);

    let deviation = point_to_segment_distance(pmid, p0, p1);

    if deviation > max_dev {
        // Subdivide
        uniform_deflection_recursive(curve, t0, mid, max_dev, params);
        params.push(mid);
        uniform_deflection_recursive(curve, mid, t1, max_dev, params);
    }
}

/// Compute the distance from a point to a line segment.
fn point_to_segment_distance(point: DVec3, seg_start: DVec3, seg_end: DVec3) -> f64 {
    let seg = seg_end - seg_start;
    let seg_len_sq = seg.length_squared();

    if seg_len_sq < TOLERANCE_ABS * TOLERANCE_ABS {
        return (point - seg_start).length();
    }

    // Project point onto line
    let t = ((point - seg_start).dot(seg) / seg_len_sq).clamp(0.0, 1.0);
    let projection = seg_start + t * seg;
    (point - projection).length()
}

/// Adaptive curve sampling with tolerance.
///
/// Samples the curve adaptively to capture all features within tolerance.
pub fn adaptive_sample_curve(curve: &Curve3, tol: f64) -> Vec<f64> {
    let domain = curve_domain(curve);

    // Start with endpoints
    let mut params = vec![domain[0]];

    // Use both deviation and angular criteria
    adaptive_sample_recursive(curve, domain[0], domain[1], tol, &mut params);

    params.push(domain[1]);
    params.sort_by(|a, b| a.partial_cmp(b).unwrap());
    params.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE_ABS);

    params
}

/// Recursive helper for adaptive sampling.
fn adaptive_sample_recursive(
    curve: &Curve3,
    t0: f64,
    t1: f64,
    tol: f64,
    params: &mut Vec<f64>,
) {
    let p0 = curve.point_at(t0);
    let p1 = curve.point_at(t1);

    // Check deviation at midpoint
    let mid = (t0 + t1) / 2.0;
    let pmid = curve.point_at(mid);

    let deviation = point_to_segment_distance(pmid, p0, p1);

    // Check angular deviation (tangent change)
    let d0 = curve_derivative(curve, t0).normalize_or_zero();
    let d1 = curve_derivative(curve, t1).normalize_or_zero();
    let angle_dev = (1.0 - d0.dot(d1)).max(0.0); // 0 for same direction, 2 for opposite

    if deviation > tol || angle_dev > 0.1 {
        // Subdivide
        adaptive_sample_recursive(curve, t0, mid, tol, params);
        params.push(mid);
        adaptive_sample_recursive(curve, mid, t1, tol, params);
    }
}

// =============================================================================
// TangentialDeflection - Sample by tangent deviation
// =============================================================================

/// Sample a curve based on tangent deviation.
///
/// Returns parameters where the tangent direction changes by at most angle_tol
/// and curvature is within curvature_tol.
pub fn tangential_deflection(curve: &Curve3, angle_tol: f64, curvature_tol: f64) -> Vec<f64> {
    let domain = curve_domain(curve);
    let mut params = vec![domain[0]];

    tangential_deflection_recursive(curve, domain[0], domain[1], angle_tol, curvature_tol, &mut params);

    params.push(domain[1]);
    params.sort_by(|a, b| a.partial_cmp(b).unwrap());
    params.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE_ABS);

    params
}

/// Recursive helper for tangential deflection sampling.
fn tangential_deflection_recursive(
    curve: &Curve3,
    t0: f64,
    t1: f64,
    angle_tol: f64,
    curvature_tol: f64,
    params: &mut Vec<f64>,
) {
    // Compute tangents at endpoints
    let tan0 = curve_derivative(curve, t0).normalize_or_zero();
    let tan1 = curve_derivative(curve, t1).normalize_or_zero();

    // Angular deviation
    let dot = tan0.dot(tan1);
    let angle = dot.acos().max(0.0);

    // Compute curvature at midpoint
    let mid = (t0 + t1) / 2.0;
    let d1 = curve_derivative(curve, mid);
    let d2 = curve_second_derivative(curve, mid);

    let d1_len = d1.length();
    let curvature = if d1_len > TOLERANCE_ABS {
        (d1.cross(d2)).length() / (d1_len * d1_len * d1_len)
    } else {
        0.0
    };

    // Check if subdivision is needed
    if angle > angle_tol || curvature > curvature_tol {
        let mid = (t0 + t1) / 2.0;
        tangential_deflection_recursive(curve, t0, mid, angle_tol, curvature_tol, params);
        params.push(mid);
        tangential_deflection_recursive(curve, mid, t1, angle_tol, curvature_tol, params);
    }
}

// =============================================================================
// QuasiUniformAbscissa - Quasi-uniform sampling
// =============================================================================

/// Quasi-uniform parameter distribution.
///
/// Distributes points quasi-uniformly, accounting for curvature to ensure
/// better spacing in high-curvature regions.
pub fn quasi_uniform(curve: &Curve3, n_points: usize) -> Vec<f64> {
    if n_points == 0 {
        return vec![];
    }

    let domain = curve_domain(curve);

    if n_points == 1 {
        return vec![(domain[0] + domain[1]) / 2.0];
    }

    // Compute curvature-weighted distribution
    let n_samples = n_points * 10;
    let mut weights = Vec::with_capacity(n_samples + 1);

    for i in 0..=n_samples {
        let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
        let d1 = curve_derivative(curve, t);
        let d2 = curve_second_derivative(curve, t);
        let d1_len = d1.length();

        let curvature = if d1_len > TOLERANCE_ABS {
            (d1.cross(d2)).length() / (d1_len * d1_len * d1_len)
        } else {
            0.0
        };

        // Weight is 1 + curvature factor for better distribution in curved regions
        weights.push(1.0 + curvature * 10.0);
    }

    // Compute cumulative weight for inverse transform sampling
    let mut cum_weight = vec![0.0_f64; n_samples + 1];
    let mut total_weight = 0.0;

    for i in 0..=n_samples {
        total_weight += weights[i];
        cum_weight[i] = total_weight;
    }

    // Sample at regular weight intervals
    let mut params = Vec::with_capacity(n_points);
    for i in 0..n_points {
        let target_weight = total_weight * i as f64 / (n_points - 1) as f64;

        // Find corresponding parameter
        let idx = cum_weight
            .binary_search_by(|probe| probe.partial_cmp(&target_weight).unwrap())
            .unwrap_or_else(|x| x.min(n_samples));

        let t = domain[0] + (domain[1] - domain[0]) * idx as f64 / n_samples as f64;
        params.push(t);
    }

    // Ensure endpoints are included
    if !params.is_empty() {
        params[0] = domain[0];
        params[n_points - 1] = domain[1];
    }

    params
}

// =============================================================================
// SurfaceSampling - Sample on surfaces
// =============================================================================

/// Sample a surface uniformly in UV space.
///
/// Returns a grid of n_u x n_v points on the surface.
pub fn sample_surface_uniform(surface: &Surface3, n_u: usize, n_v: usize) -> Vec<DVec3> {
    if n_u == 0 || n_v == 0 {
        return vec![];
    }

    let domain = surface.default_domain();

    // Handle infinite domains by clamping to reasonable range
    let (u0, u1) = if domain[0].is_infinite() || domain[1].is_infinite() {
        (-10.0, 10.0)
    } else {
        (domain[0], domain[1])
    };
    let (v0, v1) = if domain[2].is_infinite() || domain[3].is_infinite() {
        (-10.0, 10.0)
    } else {
        (domain[2], domain[3])
    };

    let mut points = Vec::with_capacity(n_u * n_v);

    for i in 0..n_u {
        let u = u0 + (u1 - u0) * i as f64 / (n_u - 1).max(1) as f64;
        for j in 0..n_v {
            let v = v0 + (v1 - v0) * j as f64 / (n_v - 1).max(1) as f64;
            points.push(surface.point_at(u, v));
        }
    }

    points
}

/// Sample a surface at regular UV steps.
///
/// Returns points sampled at u_step and v_step intervals in parameter space.
pub fn sample_surface_grid(surface: &Surface3, u_step: f64, v_step: f64) -> Vec<DVec3> {
    if u_step <= 0.0 || v_step <= 0.0 {
        return vec![];
    }

    let domain = surface.default_domain();

    // Handle infinite domains by returning empty (grid sampling doesn't make sense)
    if domain[0].is_infinite() || domain[1].is_infinite()
        || domain[2].is_infinite() || domain[3].is_infinite()
    {
        return vec![];
    }

    let u0 = domain[0];
    let u1 = domain[1];
    let v0 = domain[2];
    let v1 = domain[3];

    let n_u = ((u1 - u0) / u_step).ceil() as usize + 1;
    let n_v = ((v1 - v0) / v_step).ceil() as usize + 1;

    let mut points = Vec::with_capacity(n_u * n_v);

    let mut u = u0;
    while u <= u1 + TOLERANCE_ABS {
        let mut v = v0;
        while v <= v1 + TOLERANCE_ABS {
            let u_clamped = u.min(u1);
            let v_clamped = v.min(v1);
            points.push(surface.point_at(u_clamped, v_clamped));
            v += v_step;
        }
        u += u_step;
    }

    points
}

/// Sample a surface adaptively based on curvature.
///
/// Uses Gaussian curvature to refine sampling in high-curvature regions.
pub fn sample_surface_adaptive(surface: &Surface3, tol: f64, max_points: usize) -> Vec<DVec3> {
    let domain = surface.default_domain();

    // Start with a coarse grid
    let n_init = 5;
    let mut points = Vec::new();

    for i in 0..=n_init {
        for j in 0..=n_init {
            let u = domain[0] + (domain[1] - domain[0]) * i as f64 / n_init as f64;
            let v = domain[2] + (domain[3] - domain[2]) * j as f64 / n_init as f64;
            points.push(surface.point_at(u, v));
        }
    }

    // Subdivide based on surface deviation
    let mut count = (n_init + 1) * (n_init + 1);
    sample_surface_adaptive_recursive(
        surface,
        domain[0],
        domain[1],
        domain[2],
        domain[3],
        tol,
        &mut points,
        &mut count,
        max_points,
    );

    points
}

/// Recursive helper for adaptive surface sampling.
fn sample_surface_adaptive_recursive(
    surface: &Surface3,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    tol: f64,
    points: &mut Vec<DVec3>,
    count: &mut usize,
    max_points: usize,
) {
    if *count >= max_points {
        return;
    }

    // Check deviation at patch center
    let um = (u0 + u1) / 2.0;
    let vm = (v0 + v1) / 2.0;

    let p00 = surface.point_at(u0, v0);
    let p01 = surface.point_at(u0, v1);
    let p10 = surface.point_at(u1, v0);
    let p11 = surface.point_at(u1, v1);
    let pm = surface.point_at(um, vm);

    // Bilinear interpolation at center
    let bilinear_center = (p00 + p01 + p10 + p11) / 4.0;
    let deviation = (pm - bilinear_center).length();

    if deviation > tol && *count < max_points {
        // Add center point
        points.push(pm);
        *count += 1;

        // Recurse into 4 quadrants
        sample_surface_adaptive_recursive(surface, u0, um, v0, vm, tol, points, count, max_points);
        sample_surface_adaptive_recursive(surface, um, u1, v0, vm, tol, points, count, max_points);
        sample_surface_adaptive_recursive(surface, u0, um, vm, v1, tol, points, count, max_points);
        sample_surface_adaptive_recursive(surface, um, u1, vm, v1, tol, points, count, max_points);
    }
}

/// Sample isoparametric curves on a surface.
///
/// Returns points along u-isolines at given v parameters.
pub fn sample_u_isolines(surface: &Surface3, n_iso: usize, n_points_per_iso: usize) -> Vec<DVec3> {
    if n_iso == 0 || n_points_per_iso == 0 {
        return vec![];
    }

    let domain = surface.default_domain();

    // Handle infinite domains by clamping to reasonable range
    let (u0, u1) = if domain[0].is_infinite() || domain[1].is_infinite() {
        (-10.0, 10.0)
    } else {
        (domain[0], domain[1])
    };
    let (v0, v1) = if domain[2].is_infinite() || domain[3].is_infinite() {
        (-10.0, 10.0)
    } else {
        (domain[2], domain[3])
    };

    let mut points = Vec::with_capacity(n_iso * n_points_per_iso);

    for i in 0..n_iso {
        let v = v0 + (v1 - v0) * i as f64 / (n_iso - 1).max(1) as f64;
        for j in 0..n_points_per_iso {
            let u = u0 + (u1 - u0) * j as f64 / (n_points_per_iso - 1).max(1) as f64;
            points.push(surface.point_at(u, v));
        }
    }

    points
}

/// Sample isoparametric curves on a surface.
///
/// Returns points along v-isolines at given u parameters.
pub fn sample_v_isolines(surface: &Surface3, n_iso: usize, n_points_per_iso: usize) -> Vec<DVec3> {
    if n_iso == 0 || n_points_per_iso == 0 {
        return vec![];
    }

    let domain = surface.default_domain();

    // Handle infinite domains by clamping to reasonable range
    let (u0, u1) = if domain[0].is_infinite() || domain[1].is_infinite() {
        (-10.0, 10.0)
    } else {
        (domain[0], domain[1])
    };
    let (v0, v1) = if domain[2].is_infinite() || domain[3].is_infinite() {
        (-10.0, 10.0)
    } else {
        (domain[2], domain[3])
    };

    let mut points = Vec::with_capacity(n_iso * n_points_per_iso);

    for i in 0..n_iso {
        let u = u0 + (u1 - u0) * i as f64 / (n_iso - 1).max(1) as f64;
        for j in 0..n_points_per_iso {
            let v = v0 + (v1 - v0) * j as f64 / (n_points_per_iso - 1).max(1) as f64;
            points.push(surface.point_at(u, v));
        }
    }

    points
}

// =============================================================================
// Utility Functions
// =============================================================================

/// Get the domain for a curve, handling infinite lines specially.
fn curve_domain(curve: &Curve3) -> [f64; 2] {
    match curve {
        Curve3::Line(_) => [-1e6, 1e6], // Clamp infinite lines to large range
        other => other.default_domain(),
    }
}

/// Compute the bounding box of sampled points.
pub fn sampled_points_bounds(points: &[DVec3]) -> (DVec3, DVec3) {
    if points.is_empty() {
        return (DVec3::ZERO, DVec3::ZERO);
    }

    let mut min_pt = points[0];
    let mut max_pt = points[0];

    for pt in points.iter().skip(1) {
        min_pt = min_pt.min(*pt);
        max_pt = max_pt.max(*pt);
    }

    (min_pt, max_pt)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle3, Line3, Plane, SphericalSurface};

    fn create_test_line() -> Curve3 {
        let start = DVec3::new(0.0, 0.0, 0.0);
        let end = DVec3::new(10.0, 0.0, 0.0);
        Curve3::Line(Line3 {
            origin: start,
            direction: (end - start).normalize(),
        })
    }

    fn create_test_circle() -> Curve3 {
        Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, 0.0, 0.0),
            normal: DVec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        })
    }

    fn create_test_plane() -> Surface3 {
        Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 0.0, 0.0),
            normal: DVec3::new(0.0, 0.0, 1.0),
        })
    }

    fn create_test_sphere() -> Surface3 {
        Surface3::Sphere(SphericalSurface {
            center: DVec3::new(0.0, 0.0, 0.0),
            axis: DVec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        })
    }

    #[test]
    fn test_arc_length_line() {
        let line = create_test_line();
        let domain = curve_domain(&line);
        // Domain is [-1e6, 1e6], so arc length = 2e6
        let len = arc_length(&line, domain[0], domain[1]);
        assert!((len - 2e6).abs() < 100.0, "Expected ~2e6, got {}", len);
    }

    #[test]
    fn test_arc_length_circle() {
        let circle = create_test_circle();
        let len = total_arc_length(&circle);
        let expected = 2.0 * std::f64::consts::PI * 5.0;
        assert!((len - expected).abs() < 0.001, "Expected {}, got {}", expected, len);
    }

    #[test]
    fn test_point_at_arc_length_circle() {
        let circle = create_test_circle();
        let quarter_circ = 0.25 * 2.0 * std::f64::consts::PI * 5.0;

        let (pt, _param) = point_at_arc_length(&circle, quarter_circ);

        // Based on the circle parameterization:
        // At t=0: point is at (0, radius, 0) = (0, 5, 0)
        // At t=pi/2: point is at (-radius, 0, 0) = (-5, 0, 0)
        // So quarter arc length from t=0 should give approximately (-5, 0, 0)
        assert!((pt.x - (-5.0)).abs() < 0.5, "Expected x~-5, got {}", pt.x);
        assert!(pt.y.abs() < 0.5, "Expected y~0, got {}", pt.y);
    }

    #[test]
    fn test_points_at_equal_arc_length() {
        let circle = create_test_circle();
        let points = points_at_equal_arc_length(&circle, 5);

        assert_eq!(points.len(), 5);

        // For 5 points equally spaced by arc length on a circle,
        // the Euclidean distance (chord length) between consecutive points is:
        // 2 * r * sin(angle/2) where angle = 2π / 4 = π/2
        // So chord length = 2 * 5 * sin(π/4) = 10 * sqrt(2)/2 = 5 * sqrt(2)
        let expected_chord = 5.0 * std::f64::consts::SQRT_2;
        for i in 1..points.len() {
            let dist = (points[i].0 - points[i - 1].0).length();
            assert!((dist - expected_chord).abs() < 0.2, "Expected chord {}, got {}", expected_chord, dist);
        }
    }

    #[test]
    fn test_uniform_abscissa() {
        let circle = create_test_circle();
        let params = uniform_abscissa(&circle, 10);

        assert_eq!(params.len(), 10);

        // Check monotonicity
        for i in 1..params.len() {
            assert!(params[i] > params[i - 1], "Parameters should be increasing");
        }
    }

    #[test]
    fn test_uniform_abscissa_points() {
        let circle = create_test_circle();
        let points = uniform_abscissa_points(&circle, 8);

        assert_eq!(points.len(), 8);

        // All points should be at radius 5 from center
        for pt in &points {
            let r = pt.length();
            assert!((r - 5.0).abs() < 0.01, "Expected radius 5, got {}", r);
        }
    }

    #[test]
    fn test_uniform_deflection() {
        let circle = create_test_circle();
        let params = uniform_deflection(&circle, 0.1);

        // Circle should require many points for small deviation
        assert!(params.len() > 10, "Expected more than 10 points, got {}", params.len());

        // First and last should be close to domain bounds
        let domain = curve_domain(&circle);
        assert!((params[0] - domain[0]).abs() < 0.01);
        assert!((params[params.len() - 1] - domain[1]).abs() < 0.01);
    }

    #[test]
    fn test_adaptive_sample_curve() {
        let circle = create_test_circle();
        let params = adaptive_sample_curve(&circle, 0.01);

        // Should have endpoints
        assert!(params.len() >= 2, "Expected at least 2 points, got {}", params.len());

        // Check that we have a reasonable number of sample points
        assert!(params.len() > 5, "Expected more than 5 points for circle");
    }

    #[test]
    fn test_tangential_deflection() {
        let circle = create_test_circle();
        let params = tangential_deflection(&circle, 0.1, 100.0);

        // Should have some points
        assert!(!params.is_empty());

        // Should include endpoints
        let domain = curve_domain(&circle);
        assert!((params[0] - domain[0]).abs() < 0.01);
        assert!((params[params.len() - 1] - domain[1]).abs() < 0.01);
    }

    #[test]
    fn test_quasi_uniform() {
        let circle = create_test_circle();
        let params = quasi_uniform(&circle, 10);

        assert_eq!(params.len(), 10);

        // Check endpoints
        let domain = curve_domain(&circle);
        assert!((params[0] - domain[0]).abs() < 0.01);
        assert!((params[9] - domain[1]).abs() < 0.01);
    }

    #[test]
    fn test_sample_surface_uniform() {
        let plane = create_test_plane();
        let points = sample_surface_uniform(&plane, 5, 5);

        assert_eq!(points.len(), 25);

        // All points should be at z=0 for the plane
        for pt in &points {
            assert!(pt.z.abs() < 0.01, "Expected z~0, got {}", pt.z);
        }
    }

    #[test]
    fn test_sample_surface_grid() {
        // Use sphere instead of plane since plane has infinite domain
        let sphere = create_test_sphere();
        let points = sample_surface_grid(&sphere, 0.5, 0.5);

        // Should have multiple points
        assert!(!points.is_empty());

        // All points should be at radius ~5
        for pt in &points {
            let r = pt.length();
            assert!((r - 5.0).abs() < 0.1, "Expected radius 5, got {}", r);
        }
    }

    #[test]
    fn test_sample_surface_uniform_sphere() {
        let sphere = create_test_sphere();
        let points = sample_surface_uniform(&sphere, 11, 11);

        assert_eq!(points.len(), 121);

        // All points should be at radius ~5 from center
        for pt in &points {
            let r = pt.length();
            assert!((r - 5.0).abs() < 0.01, "Expected radius 5, got {}", r);
        }
    }

    #[test]
    fn test_sample_u_isolines() {
        let plane = create_test_plane();
        let points = sample_u_isolines(&plane, 3, 5);

        assert_eq!(points.len(), 15);
    }

    #[test]
    fn test_sample_v_isolines() {
        let plane = create_test_plane();
        let points = sample_v_isolines(&plane, 4, 6);

        assert_eq!(points.len(), 24);
    }

    #[test]
    fn test_sampled_points_bounds() {
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(-1.0, -2.0, -3.0),
        ];

        let (min_pt, max_pt) = sampled_points_bounds(&points);

        assert_eq!(min_pt, DVec3::new(-1.0, -2.0, -3.0));
        assert_eq!(max_pt, DVec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_sample_surface_adaptive() {
        let sphere = create_test_sphere();
        let points = sample_surface_adaptive(&sphere, 0.1, 1000);

        // Should have some points
        assert!(!points.is_empty());

        // All points should be at radius ~5
        for pt in &points {
            let r = pt.length();
            assert!((r - 5.0).abs() < 0.1, "Expected radius 5, got {}", r);
        }
    }

    #[test]
    fn test_point_to_segment_distance() {
        let seg_start = DVec3::new(0.0, 0.0, 0.0);
        let seg_end = DVec3::new(10.0, 0.0, 0.0);

        // Point on the segment
        let dist = point_to_segment_distance(DVec3::new(5.0, 0.0, 0.0), seg_start, seg_end);
        assert!(dist.abs() < 0.001, "Expected 0, got {}", dist);

        // Point perpendicular to segment
        let dist = point_to_segment_distance(DVec3::new(5.0, 5.0, 0.0), seg_start, seg_end);
        assert!((dist - 5.0).abs() < 0.001, "Expected 5, got {}", dist);

        // Point beyond segment end
        let dist = point_to_segment_distance(DVec3::new(15.0, 5.0, 0.0), seg_start, seg_end);
        assert!((dist - 7.071).abs() < 0.01, "Expected ~7.071, got {}", dist);
    }

    #[test]
    fn test_arc_length_partial() {
        let circle = create_test_circle();

        // Test half circle
        let half_len = arc_length(&circle, 0.0, std::f64::consts::PI);
        let expected = std::f64::consts::PI * 5.0;
        assert!((half_len - expected).abs() < 0.001, "Expected {}, got {}", expected, half_len);
    }

    #[test]
    fn test_zero_points() {
        let circle = create_test_circle();

        assert!(uniform_abscissa(&circle, 0).is_empty());
        assert!(points_at_equal_arc_length(&circle, 0).is_empty());
        assert!(quasi_uniform(&circle, 0).is_empty());
    }

    #[test]
    fn test_single_point() {
        let circle = create_test_circle();

        let params = uniform_abscissa(&circle, 1);
        assert_eq!(params.len(), 1);

        let params = quasi_uniform(&circle, 1);
        assert_eq!(params.len(), 1);
    }

    #[test]
    fn test_empty_surface_samples() {
        let plane = create_test_plane();

        assert!(sample_surface_uniform(&plane, 0, 5).is_empty());
        assert!(sample_surface_uniform(&plane, 5, 0).is_empty());
        assert!(sample_surface_grid(&plane, 0.0, 1.0).is_empty());
        assert!(sample_surface_grid(&plane, 1.0, 0.0).is_empty());
    }
}
