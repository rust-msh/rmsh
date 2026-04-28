use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

/// A numerically sampled intersection curve.
#[derive(Debug, Clone)]
pub struct SampledCurve {
    pub points: Vec<DVec3>,
    pub is_closed: bool,
    /// Diagnostic: number of oscillation events detected during marching.
    pub oscillation_count: usize,
    /// Diagnostic: whether step size was reduced during marching.
    pub step_reduced: bool,
}

impl Default for SampledCurve {
    fn default() -> Self {
        Self {
            points: Vec::new(),
            is_closed: false,
            oscillation_count: 0,
            step_reduced: false,
        }
    }
}

/// Configuration for adaptive marching behavior.
#[derive(Debug, Clone, Copy)]
pub struct MarchingConfig {
    /// Initial step size for curve tracing.
    pub step_size: f64,
    /// Minimum allowed step size (for convergence failure fallback).
    pub min_step_size: f64,
    /// Maximum number of steps per direction.
    pub max_steps: usize,
    /// Maximum allowed oscillations before step reduction.
    pub max_oscillations: usize,
    /// Factor to reduce step size when oscillation is detected.
    pub step_reduction_factor: f64,
    /// Enable multi-scale seed detection.
    pub multiscale_seeds: bool,
}

impl Default for MarchingConfig {
    fn default() -> Self {
        Self {
            step_size: 0.1,
            min_step_size: 1e-10,
            max_steps: 500,
            max_oscillations: 3,
            step_reduction_factor: 0.5,
            multiscale_seeds: false,
        }
    }
}

/// Result of adaptive sampling density calculation.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveSampling {
    /// Number of samples in u-direction.
    pub n_u: usize,
    /// Number of samples in v-direction.
    pub n_v: usize,
    /// Estimated characteristic length for step size.
    pub characteristic_length: f64,
}

/// Compute adaptive sampling density based on surface geometry.
pub fn adaptive_sampling_density(surface: &Surface3, base_density: usize) -> AdaptiveSampling {
    let default = AdaptiveSampling {
        n_u: base_density,
        n_v: base_density,
        characteristic_length: 1.0,
    };

    match surface {
        Surface3::Cylinder(c) => {
            // Cylinder: u = azimuth, v = height
            // Azimuth should be proportional to circumference (2πR)
            // Height should cover the expected range
            let circumference = std::f64::consts::TAU * c.radius;
            let n_u = (base_density as f64 * (circumference / c.radius).sqrt()).ceil() as usize;
            let n_v = base_density;
            AdaptiveSampling {
                n_u: n_u.max(base_density / 2).min(base_density * 2),
                n_v: n_v.max(base_density / 2),
                characteristic_length: c.radius * 0.1,
            }
        }
        Surface3::Sphere(s) => {
            // Sphere: uniform sampling based on radius
            let n = (base_density as f64 * (s.radius / 1.0).sqrt()).ceil() as usize;
            AdaptiveSampling {
                n_u: n.max(base_density),
                n_v: n.max(base_density),
                characteristic_length: s.radius * 0.1,
            }
        }
        Surface3::Torus(t) => {
            // Torus: major radius for u, minor radius for v
            let ratio = t.major_radius / t.minor_radius.max(1e-10);
            let n_u = (base_density as f64 * ratio.sqrt()).ceil() as usize;
            let n_v = base_density;
            AdaptiveSampling {
                n_u: n_u.max(base_density).min(base_density * 3),
                n_v: n_v.max(base_density / 2),
                characteristic_length: t.minor_radius * 0.1,
            }
        }
        Surface3::Cone(c) => {
            // Cone: similar to cylinder but with varying radius
            let avg_radius = c.radius * 0.5; // approximate average
            let n_u = (base_density as f64 * (avg_radius / 1.0).sqrt()).ceil() as usize;
            let n_v = base_density;
            AdaptiveSampling {
                n_u: n_u.max(base_density / 2),
                n_v: n_v.max(base_density / 2),
                characteristic_length: avg_radius * 0.1,
            }
        }
        Surface3::BSpline(bs) => {
            // BSpline: estimate from control point bounding box
            let bbox = estimate_bspline_bbox(bs);
            let max_extent = (bbox.1 - bbox.0).max_element();
            let n = (base_density as f64 * (max_extent / 1.0).sqrt()).ceil() as usize;
            AdaptiveSampling {
                n_u: n.max(base_density),
                n_v: n.max(base_density),
                characteristic_length: max_extent * 0.05,
            }
        }
        _ => default,
    }
}

/// Estimate bounding box of a BSpline surface from control points.
fn estimate_bspline_bbox(bs: &BSplineSurface) -> (DVec3, DVec3) {
    let mut min_pt = DVec3::splat(f64::INFINITY);
    let mut max_pt = DVec3::splat(f64::NEG_INFINITY);

    for row in &bs.control_points {
        for pt in row {
            min_pt = min_pt.min(*pt);
            max_pt = max_pt.max(*pt);
        }
    }

    if !min_pt.is_finite() {
        min_pt = DVec3::splat(-1.0);
    }
    if !max_pt.is_finite() {
        max_pt = DVec3::splat(1.0);
    }

    (min_pt, max_pt)
}

/// Coarse UV grid search (5×5 over [0,1]²) to find the closest surface sample.
/// Returns (u, v) of the closest grid point.
fn closest_uv_coarse(surface: &Surface3, point: DVec3) -> (f64, f64) {
    const N: usize = 5;
    let mut best_u = 0.5_f64;
    let mut best_v = 0.5_f64;
    let mut best_dist_sq = f64::MAX;
    for i in 0..N {
        for j in 0..N {
            let u = i as f64 / (N - 1) as f64;
            let v = j as f64 / (N - 1) as f64;
            let p = surface.point_at(u, v);
            let d = (p - point).length_squared();
            if d < best_dist_sq {
                best_dist_sq = d;
                best_u = u;
                best_v = v;
            }
        }
    }
    (best_u, best_v)
}

/// Evaluate the implicit function F(P) for a surface: F=0 on surface.
pub fn surface_implicit(surface: &Surface3, point: DVec3) -> f64 {
    match surface {
        Surface3::Plane(p) => (point - p.origin).dot(p.normal),
        Surface3::Cylinder(c) => {
            let v = point - c.origin;
            let along = v.dot(c.axis);
            let perp = v - c.axis * along;
            perp.length() - c.radius
        }
        Surface3::Sphere(s) => (point - s.center).length() - s.radius,
        Surface3::Cone(c) => {
            let axis = c.axis_dir();
            let v = point - c.apex;
            let along = v.dot(axis);
            let perp_len = (v - axis * along).length();
            perp_len - c.radius_at_axial(along)
        }
        Surface3::Torus(t) => {
            let v = point - t.center;
            let along = v.dot(t.axis);
            let perp = v - t.axis * along;
            let perp_len = perp.length();
            let d = perp_len - t.major_radius;
            (d * d + along * along).sqrt() - t.minor_radius
        }
        _ => {
            // Closest-point signed distance: project onto the nearest surface sample.
            let (u, v) = closest_uv_coarse(surface, point);
            let closest = surface.point_at(u, v);
            let normal = surface.normal_at(u, v);
            let n_len = normal.length();
            if n_len < 1e-14 {
                return (point - closest).length();
            }
            (point - closest).dot(normal / n_len)
        }
    }
}

/// Compute the gradient ∇F at a point for a surface.
fn surface_gradient(surface: &Surface3, point: DVec3) -> DVec3 {
    match surface {
        Surface3::Plane(p) => p.normal,
        Surface3::Cylinder(c) => {
            let v = point - c.origin;
            let along = v.dot(c.axis);
            let perp = v - c.axis * along;
            let perp_len = perp.length();
            if perp_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            perp / perp_len
        }
        Surface3::Sphere(s) => {
            let v = point - s.center;
            let len = v.length();
            if len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            v / len
        }
        Surface3::Cone(c) => {
            let axis = c.axis_dir();
            let v = point - c.apex;
            let along = v.dot(axis);
            let perp = v - axis * along;
            let perp_len = perp.length();
            if perp_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            let tan_a = c.half_angle_rad.tan();
            perp / perp_len - axis * tan_a
        }
        Surface3::Torus(t) => {
            let v = point - t.center;
            let along = v.dot(t.axis);
            let perp = v - t.axis * along;
            let perp_len = perp.length();
            if perp_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            let tube_center = t.center + perp / perp_len * t.major_radius;
            let tv = point - tube_center;
            let tv_len = tv.length();
            if tv_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            tv / tv_len
        }
        _ => {
            let (u, v) = closest_uv_coarse(surface, point);
            let normal = surface.normal_at(u, v);
            let n_len = normal.length();
            if n_len < 1e-14 {
                return DVec3::ZERO;
            }
            normal / n_len
        }
    }
}

/// Project a point onto a surface using Newton iteration.
pub fn project_onto_surface(surface: &Surface3, point: DVec3, max_iter: usize) -> DVec3 {
    let mut p = point;
    for _ in 0..max_iter {
        let f = surface_implicit(surface, p);
        if f.abs() < TOLERANCE_ABS {
            break;
        }
        let g = surface_gradient(surface, p);
        let g_len_sq = g.length_squared();
        if g_len_sq < TOLERANCE_ABS * TOLERANCE_ABS {
            break;
        }
        p -= g * (f / g_len_sq);
    }
    p
}

/// Project a point onto the intersection of two surfaces.
fn project_onto_intersection(s1: &Surface3, s2: &Surface3, point: DVec3) -> DVec3 {
    let mut p = point;
    for _ in 0..50 {
        let f1 = surface_implicit(s1, p);
        let f2 = surface_implicit(s2, p);
        if f1.abs() < TOLERANCE_ABS && f2.abs() < TOLERANCE_ABS {
            break;
        }
        let g1 = surface_gradient(s1, p);
        let g2 = surface_gradient(s2, p);

        // Solve 2x2 system: move by λ1*g1 + λ2*g2 to satisfy both constraints
        let a11 = g1.dot(g1);
        let a12 = g1.dot(g2);
        let a22 = g2.dot(g2);
        let det = a11 * a22 - a12 * a12;
        if det.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
            // Degenerate — just project onto each surface alternately
            p = project_onto_surface(s1, p, 5);
            p = project_onto_surface(s2, p, 5);
            continue;
        }
        let lambda1 = (a22 * f1 - a12 * f2) / det;
        let lambda2 = (a11 * f2 - a12 * f1) / det;
        p -= g1 * lambda1 + g2 * lambda2;
    }
    p
}

/// Find seed points for intersection curve marching by sampling one surface.
pub fn find_seed_points(s1: &Surface3, s2: &Surface3, sample_points: &[DVec3]) -> Vec<DVec3> {
    let mut seeds = Vec::new();

    // Look for sign changes of F2 along the sample points
    let values: Vec<f64> = sample_points
        .iter()
        .map(|&p| surface_implicit(s2, p))
        .collect();

    for i in 0..values.len().saturating_sub(1) {
        if values[i] * values[i + 1] < 0.0 {
            // Sign change — interpolate
            let t = values[i] / (values[i] - values[i + 1]);
            let p = sample_points[i] + (sample_points[i + 1] - sample_points[i]) * t;
            let seed = project_onto_intersection(s1, s2, p);
            seeds.push(seed);
        }
    }

    seeds
}

/// Like `find_seed_points` but treats `sample_points` as a `n_u × n_v` grid
/// (row-major: index = iu * n_v + iv) and checks sign changes along BOTH the
/// u-direction and v-direction edges. This avoids missing seeds when the
/// intersection curve runs along one of the grid directions.
pub fn find_seed_points_grid(
    s1: &Surface3,
    s2: &Surface3,
    sample_points: &[DVec3],
    n_u: usize,
    n_v: usize,
) -> Vec<DVec3> {
    assert_eq!(sample_points.len(), n_u * n_v, "grid size mismatch");
    let mut seeds = Vec::new();

    let values: Vec<f64> = sample_points
        .iter()
        .map(|&p| surface_implicit(s2, p))
        .collect();

    let idx = |iu: usize, iv: usize| iu * n_v + iv;

    // Check u-direction edges (vary iv, fixed iu)
    for iu in 0..n_u {
        for iv in 0..n_v.saturating_sub(1) {
            let a = idx(iu, iv);
            let b = idx(iu, iv + 1);
            if values[a] * values[b] < 0.0 {
                let t = values[a] / (values[a] - values[b]);
                let p = sample_points[a].lerp(sample_points[b], t);
                seeds.push(project_onto_intersection(s1, s2, p));
            }
        }
    }

    // Check v-direction edges (vary iu, fixed iv)
    for iu in 0..n_u.saturating_sub(1) {
        for iv in 0..n_v {
            let a = idx(iu, iv);
            let b = idx(iu + 1, iv);
            if values[a] * values[b] < 0.0 {
                let t = values[a] / (values[a] - values[b]);
                let p = sample_points[a].lerp(sample_points[b], t);
                seeds.push(project_onto_intersection(s1, s2, p));
            }
        }
    }

    seeds
}

/// Multi-scale seed point detection with deduplication.
/// Runs seed detection at multiple grid resolutions and merges results.
pub fn find_seed_points_multiscale(
    s1: &Surface3,
    s2: &Surface3,
    sampler: impl Fn(usize, usize) -> Vec<DVec3>,
    scales: &[usize],
    dedup_tolerance: f64,
) -> Vec<DVec3> {
    let mut all_seeds = Vec::new();
    let dedup_tol_sq = dedup_tolerance * dedup_tolerance;

    for &n in scales {
        let n_u = n;
        let n_v = n;
        let samples = sampler(n_u, n_v);
        let seeds = find_seed_points_grid(s1, s2, &samples, n_u, n_v);

        for seed in seeds {
            // Deduplicate: skip if too close to an existing seed
            let is_dup = all_seeds.iter().any(|s: &DVec3| (*s - seed).length_squared() < dedup_tol_sq);
            if !is_dup {
                all_seeds.push(seed);
            }
        }
    }

    all_seeds
}

/// March an intersection curve starting from a seed point.
/// Traces in both directions along the curve until it returns to start
/// (closed) or exits bounds.
pub fn march_intersection(
    s1: &Surface3,
    s2: &Surface3,
    seed: DVec3,
    step_size: f64,
    max_steps: usize,
    bounds_check: impl Fn(DVec3) -> bool,
) -> SampledCurve {
    let config = MarchingConfig {
        step_size,
        max_steps,
        ..Default::default()
    };
    // Use the simpler implementation without Clone bound
    let mut result = SampledCurve::default();

    // Try forward direction
    let forward = march_one_direction_monitored_simple(
        s1, s2, seed, config.step_size, config.max_steps,
        &bounds_check, config.max_oscillations, config.min_step_size,
        config.step_reduction_factor,
    );

    // Try backward direction
    let backward = march_one_direction_monitored_simple(
        s1, s2, seed, -config.step_size, config.max_steps,
        &bounds_check, config.max_oscillations, config.min_step_size,
        config.step_reduction_factor,
    );

    // Combine: reverse backward, then append forward (excluding duplicate seed)
    result.points = backward.points.into_iter().rev().collect();
    if !forward.points.is_empty() {
        result.points.extend(forward.points.into_iter().skip(1));
    }

    // Check closure
    result.is_closed = result.points.len() > 2
        && points_coincide(result.points[0], *result.points.last().unwrap());

    if result.is_closed {
        result.points.pop();
    }

    result
}

/// March an intersection curve with full configuration and convergence monitoring.
pub fn march_intersection_with_config(
    s1: &Surface3,
    s2: &Surface3,
    seed: DVec3,
    config: &MarchingConfig,
    bounds_check: impl Fn(DVec3) -> bool,
) -> SampledCurve {
    let mut result = SampledCurve::default();

    // Try forward direction with convergence monitoring
    let forward = march_one_direction_monitored_simple(
        s1, s2, seed, config.step_size, config.max_steps,
        &bounds_check, config.max_oscillations, config.min_step_size,
        config.step_reduction_factor,
    );
    result.oscillation_count += forward.oscillation_count;
    result.step_reduced = result.step_reduced || forward.step_reduced;

    // Try backward direction
    let backward = march_one_direction_monitored_simple(
        s1, s2, seed, -config.step_size, config.max_steps,
        &bounds_check, config.max_oscillations, config.min_step_size,
        config.step_reduction_factor,
    );
    result.oscillation_count += backward.oscillation_count;
    result.step_reduced = result.step_reduced || backward.step_reduced;

    // Combine: reverse backward, then append forward (excluding duplicate seed)
    result.points = backward.points.into_iter().rev().collect();
    if !forward.points.is_empty() {
        // Skip the seed which is duplicated
        result.points.extend(forward.points.into_iter().skip(1));
    }

    // Check closure
    result.is_closed = result.points.len() > 2
        && points_coincide(result.points[0], *result.points.last().unwrap());

    if result.is_closed {
        result.points.pop();
    }

    result
}

/// Result of monitored single-direction marching.
struct MonitoredMarchResult {
    points: Vec<DVec3>,
    oscillation_count: usize,
    step_reduced: bool,
}

/// March in one direction with oscillation detection and step reduction.
/// Simple version without Clone bound on bounds_check.
fn march_one_direction_monitored_simple(
    s1: &Surface3,
    s2: &Surface3,
    seed: DVec3,
    step_size: f64,
    max_steps: usize,
    bounds_check: &impl Fn(DVec3) -> bool,
    max_oscillations: usize,
    min_step_size: f64,
    step_reduction_factor: f64,
) -> MonitoredMarchResult {
    let mut points = vec![seed];
    let mut current = seed;
    let mut oscillation_count = 0usize;
    let mut step_reduced = false;
    let mut current_step = step_size;
    let mut consecutive_oscillations = 0usize;
    let mut last_dir = DVec3::ZERO;

    // Closure tolerance: within 2× step_size of start is considered closed.
    let closure_tol_sq = (step_size * 2.0) * (step_size * 2.0);
    // Track arc length to avoid infinite loops.
    let mut arc_len = 0.0_f64;

    for _step_idx in 0..max_steps {
        let g1 = surface_gradient(s1, current);
        let g2 = surface_gradient(s2, current);
        let tangent = g1.cross(g2);
        let t_len = tangent.length();
        if t_len < TOLERANCE_ABS {
            break; // tangent surfaces — can't march
        }
        let dir = tangent / t_len * step_size.signum();

        // Oscillation detection: direction reversal
        if last_dir.length_squared() > 0.5 {
            let alignment = dir.dot(last_dir);
            if alignment < -0.9 {
                oscillation_count += 1;
                consecutive_oscillations += 1;

                // If too many consecutive oscillations, reduce step size
                if consecutive_oscillations >= max_oscillations && current_step > min_step_size {
                    current_step *= step_reduction_factor;
                    step_reduced = true;
                    consecutive_oscillations = 0;

                    // Reset current position to last good point and continue
                    if points.len() > 1 {
                        current = points[points.len() - 1];
                    }
                    continue;
                }
            } else {
                consecutive_oscillations = 0;
            }
        }

        let next_raw = current + dir * current_step.abs();
        let next = project_onto_intersection(s1, s2, next_raw);

        if !bounds_check(next) {
            break;
        }

        let step_dist = (next - current).length();
        arc_len += step_dist;

        // Check if we've returned to start (closed curve).
        if points.len() > 10 && (next - points[0]).length_squared() < closure_tol_sq {
            points.push(points[0]); // seal the loop
            break;
        }

        // Cap arc length to prevent runaway on infinite/very long open curves.
        let arc_cap = step_size * (max_steps as f64).min(400.0);
        if arc_len >= arc_cap {
            break;
        }

        points.push(next);
        current = next;
        last_dir = dir;
    }

    MonitoredMarchResult {
        points,
        oscillation_count,
        step_reduced,
    }
}

/// Generate sample points on a cylinder surface for seed finding.
pub fn sample_cylinder(
    cyl: &CylindricalSurface,
    height_range: [f64; 2],
    n_theta: usize,
    n_h: usize,
) -> Vec<DVec3> {
    let u = if cyl.axis.x.abs() < 0.9 {
        cyl.axis.cross(DVec3::X).normalize()
    } else {
        cyl.axis.cross(DVec3::Y).normalize()
    };
    let v = cyl.axis.cross(u);

    let mut points = Vec::with_capacity(n_theta * n_h);
    for ih in 0..n_h {
        let h = height_range[0]
            + (height_range[1] - height_range[0]) * ih as f64 / (n_h - 1).max(1) as f64;
        for it in 0..n_theta {
            let theta = 2.0 * std::f64::consts::PI * it as f64 / n_theta as f64;
            let p = cyl.origin + cyl.axis * h + (u * theta.cos() + v * theta.sin()) * cyl.radius;
            points.push(p);
        }
    }
    points
}

/// Generate sample points on a sphere surface for seed finding.
pub fn sample_sphere(sphere: &SphericalSurface, n_theta: usize, n_phi: usize) -> Vec<DVec3> {
    let u = if sphere.axis.x.abs() < 0.9 {
        sphere.axis.cross(DVec3::X).normalize()
    } else {
        sphere.axis.cross(DVec3::Y).normalize()
    };
    let v = sphere.axis.cross(u);

    let mut points = Vec::with_capacity(n_theta * n_phi);
    for ip in 0..n_phi {
        let phi = std::f64::consts::PI * ip as f64 / (n_phi - 1).max(1) as f64;
        for it in 0..n_theta {
            let theta = 2.0 * std::f64::consts::PI * it as f64 / n_theta as f64;
            let p = sphere.center
                + sphere.radius
                    * (sphere.axis * phi.cos() + (u * theta.cos() + v * theta.sin()) * phi.sin());
            points.push(p);
        }
    }
    points
}

/// Generate sample points on a torus surface for seed finding.
pub fn sample_torus(torus: &ToroidalSurface, n_u: usize, n_v: usize) -> Vec<DVec3> {
    let u_dir = if torus.axis.x.abs() < 0.9 {
        torus.axis.cross(DVec3::X).normalize()
    } else {
        torus.axis.cross(DVec3::Y).normalize()
    };
    let v_dir = torus.axis.cross(u_dir);

    let mut points = Vec::with_capacity(n_u * n_v);
    for iu in 0..n_u {
        let u = 2.0 * std::f64::consts::PI * iu as f64 / n_u as f64;
        let cu = u.cos();
        let su = u.sin();
        let ring_center = torus.center + (u_dir * cu + v_dir * su) * torus.major_radius;
        let ring_outward = u_dir * cu + v_dir * su;

        for iv in 0..n_v {
            let v = 2.0 * std::f64::consts::PI * iv as f64 / n_v as f64;
            let p =
                ring_center + (ring_outward * v.cos() + torus.axis * v.sin()) * torus.minor_radius;
            points.push(p);
        }
    }
    points
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn implicit_plane() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        });
        assert!((surface_implicit(&plane, DVec3::ZERO)).abs() < TOLERANCE_ABS);
        assert!((surface_implicit(&plane, DVec3::new(0.0, 1.0, 0.0)) - 1.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn implicit_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        });
        assert!((surface_implicit(&sphere, DVec3::new(2.0, 0.0, 0.0))).abs() < TOLERANCE_ABS);
        assert!((surface_implicit(&sphere, DVec3::new(1.0, 0.0, 0.0)) + 1.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn implicit_cylinder() {
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        });
        assert!((surface_implicit(&cyl, DVec3::new(3.0, 5.0, 0.0))).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn implicit_torus() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        });
        // Point on the outer equator: (6, 0, 0)
        assert!((surface_implicit(&torus, DVec3::new(6.0, 0.0, 0.0))).abs() < TOLERANCE_ABS);
        // Point on the inner equator: (4, 0, 0)
        assert!((surface_implicit(&torus, DVec3::new(4.0, 0.0, 0.0))).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn project_onto_sphere_test() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        });
        let p = project_onto_surface(&sphere, DVec3::new(3.0, 0.0, 0.0), 20);
        assert!((p.length() - 2.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn march_sphere_cylinder_intersection() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        });
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(1.0, 0.0, 0.0),
            axis: DVec3::Y,
            radius: 1.5,
        });

        // Find a seed by sampling the cylinder and looking for sign changes on sphere
        let cyl_surf = match &cyl {
            Surface3::Cylinder(c) => c,
            _ => unreachable!(),
        };
        let samples = sample_cylinder(cyl_surf, [-1.0, 1.0], 32, 4);
        let seeds = find_seed_points(&cyl, &sphere, &samples);

        assert!(!seeds.is_empty(), "Should find at least one seed point");
        let seed = seeds[0];

        // Verify seed is approximately on both surfaces
        assert!(
            surface_implicit(&sphere, seed).abs() < 0.1,
            "seed not near sphere: F={}",
            surface_implicit(&sphere, seed)
        );
        assert!(
            surface_implicit(&cyl, seed).abs() < 0.1,
            "seed not near cylinder: F={}",
            surface_implicit(&cyl, seed)
        );

        let curve = march_intersection(&sphere, &cyl, seed, 0.1, 200, |_| true);
        assert!(
            curve.points.len() > 5,
            "Expected marched curve with several points, got {}",
            curve.points.len()
        );

        // All points should be approximately on both surfaces
        for p in &curve.points {
            assert!(
                surface_implicit(&sphere, *p).abs() < 0.05,
                "point not on sphere: F={}",
                surface_implicit(&sphere, *p)
            );
            assert!(
                surface_implicit(&cyl, *p).abs() < 0.05,
                "point not on cylinder: F={}",
                surface_implicit(&cyl, *p)
            );
        }
    }

    #[test]
    fn sample_cylinder_test() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        };
        let pts = sample_cylinder(&cyl, [0.0, 2.0], 8, 4);
        assert_eq!(pts.len(), 32);
        for p in &pts {
            let r = (p.x * p.x + p.z * p.z).sqrt();
            assert!((r - 1.0).abs() < TOLERANCE_ABS);
        }
    }
}
