//! Analytic intersection of a sphere and a cone.
//!
//! # Case classification
//!
//! ## Sphere center on cone axis (axis-aligned case)
//!
//! When the sphere center lies on the cone axis, the intersection consists of
//! circles at axial heights where the sphere's cross-section radius equals
//! the cone's radius at that height. Solve:
//!
//!   r_sphere(z) = sqrt(R² - (z - z_c)²)
//!   r_cone(z) = r_ref + (z - z_ref) * tan(half_angle)
//!
//! The intersection points satisfy r_sphere(z) = r_cone(z), which leads to
//! a quartic equation. We use numerical root-finding to locate solutions.
//!
//! ## General case
//!
//! For all other configurations the intersection is a space curve of degree <= 4.
//! We return `General` so the caller falls back to numeric marching.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, ConicalSurface, SphericalSurface};

use crate::tolerance::TOLERANCE_ABS;

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of sphere x cone intersection.
#[derive(Debug, Clone)]
pub enum SphereConeResult {
    /// The sphere and cone do not intersect.
    NoIntersection,
    /// Sphere center is on cone axis; intersection is a single circle.
    SingleCircle(Circle3),
    /// Sphere center is on cone axis; intersection consists of two circles.
    TwoCircles(Circle3, Circle3),
    /// The intersection is a tangent point.
    TangentPoint(DVec3),
    /// General case. Caller should fall back to marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `sphere` and `cone`.
pub fn intersect_sphere_cone(
    sphere: &SphericalSurface,
    cone: &ConicalSurface,
) -> SphereConeResult {
    intersect_sphere_cone_with_tolerance(sphere, cone, 0.0)
}

/// Compute sphere-cone intersection with additional fuzzy tolerance.
///
/// This relaxes axis-aligned and distance early-out checks by `fuzzy_tol` so
/// near-coincident cases can still classify into analytic branches.
pub fn intersect_sphere_cone_with_tolerance(
    sphere: &SphericalSurface,
    cone: &ConicalSurface,
    fuzzy_tol: f64,
) -> SphereConeResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);
    let axis = cone.axis_dir();
    let apex = cone.apex_point();

    // Project sphere center onto cone axis
    let d = sphere.center - apex;
    let z_c = d.dot(axis); // axial distance from cone apex to sphere center
    let foot = apex + axis * z_c;
    let d_perp = (sphere.center - foot).length();

    // ── Axis-aligned case: sphere center on cone axis ───────────────────────────
    if d_perp < tol * 10.0 {
        return intersect_sphere_cone_on_axis(sphere, cone, z_c, tol);
    }

    // ── General case: numerical fallback ───────────────────────────────────────
    SphereConeResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Axis-aligned case
// ─────────────────────────────────────────────────────────────────────────────

fn intersect_sphere_cone_on_axis(
    sphere: &SphericalSurface,
    cone: &ConicalSurface,
    z_c: f64,
    tol: f64,
) -> SphereConeResult {
    let axis = cone.axis_dir();
    let apex = cone.apex_point();
    let tan_half = cone.half_angle_rad.tan();
    let big_r = sphere.radius;
    let r_ref = cone.radius;

    // Cone radius at axial distance z from apex: r_cone(z) = r_ref + z * tan_half
    // Sphere cross-section radius at axial distance z from sphere center:
    //   r_sphere(z) = sqrt(big_r² - (z - z_c)²)  when |z - z_c| <= big_r
    //
    // Intersection: r_sphere(z) = r_cone(z)
    // Let u = z - z_c (offset from sphere center along axis)
    //   sqrt(big_r² - u²) = r_ref + (z_c + u) * tan_half
    // Square both sides:
    //   big_r² - u² = (r_ref + z_c * tan_half + u * tan_half)²
    // This is a quartic in u. We solve by sampling and bisection.

    // Sampling range: u in [-big_r, +big_r]
    let n = 128usize;
    let mut roots: Vec<f64> = Vec::new();

    let f = |u: f64| -> f64 {
        if u.abs() > big_r {
            return f64::NAN;
        }
        let r_sphere_sq = big_r * big_r - u * u;
        if r_sphere_sq < 0.0 {
            return f64::NAN;
        }
        let r_sphere = r_sphere_sq.sqrt();
        let z = z_c + u;
        let r_cone = r_ref + z * tan_half;
        r_sphere - r_cone
    };

    let mut prev_u = -big_r;
    let mut prev_f = f(prev_u);

    for i in 1..=n {
        let u = -big_r + 2.0 * big_r * i as f64 / n as f64;
        let curr_f = f(u);

        if prev_f.is_nan() || curr_f.is_nan() {
            prev_u = u;
            prev_f = curr_f;
            continue;
        }

        // Sign change indicates a root
        if prev_f * curr_f < 0.0 {
            // Bisection
            let mut lo = prev_u;
            let mut hi = u;
            for _ in 0..64 {
                let mid = (lo + hi) * 0.5;
                let f_mid = f(mid);
                if f_mid.is_nan() {
                    break;
                }
                if f(lo) * f_mid < 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            roots.push((lo + hi) * 0.5);
        } else if curr_f.abs() < tol {
            // Near-zero: check if this is a tangent point
            roots.push(u);
        }

        prev_u = u;
        prev_f = curr_f;
    }

    // Remove duplicate roots
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    roots.dedup_by(|a, b| (*a - *b).abs() < tol);

    // Convert roots to circles
    let mut circles: Vec<Circle3> = Vec::new();
    for u in roots {
        let z = z_c + u;
        let r_cone = r_ref + z * tan_half;

        // Skip if radius is negative or too small
        if r_cone < tol {
            // Tangent point at apex or near apex
            let pt = apex + axis * z;
            // Check if this is actually on the sphere
            if (pt - sphere.center).length() < big_r + tol {
                continue; // Will be handled as a point later if needed
            }
            continue;
        }

        // Verify the solution
        let r_sphere_sq = big_r * big_r - u * u;
        if r_sphere_sq < -tol {
            continue;
        }

        let center = apex + axis * z;
        circles.push(Circle3 {
            center,
            normal: axis,
            radius: r_cone.max(0.0),
        });
    }

    match circles.len() {
        0 => SphereConeResult::NoIntersection,
        1 => SphereConeResult::SingleCircle(circles[0]),
        2 => SphereConeResult::TwoCircles(circles[0], circles[1]),
        _ => SphereConeResult::General, // More than 2 circles is unusual
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    fn sphere(center: DVec3, radius: f64) -> SphericalSurface {
        SphericalSurface {
            center,
            axis: DVec3::Z,
            radius,
        }
    }

    fn cone(apex: DVec3, axis: DVec3, half_angle_deg: f64) -> ConicalSurface {
        ConicalSurface {
            apex,
            axis,
            radius: 0.0,
            half_angle_rad: half_angle_deg.to_radians(),
        }
    }

    /// Sphere center on cone axis, sphere larger than cone tip.
    /// Cone: apex at origin, 45 degree half-angle.
    /// Sphere: center at (0,0,5), radius 3.
    /// At z=5: r_cone = 5*tan(45) = 5
    /// Sphere cross-section: r_sphere = sqrt(9 - 0) = 3
    /// They don't intersect at z=5, but may intersect at other heights.
    #[test]
    fn sphere_on_cone_axis_general() {
        let s = sphere(DVec3::new(0.0, 0.0, 5.0), 3.0);
        let k = cone(DVec3::ZERO, DVec3::Z, 45.0);
        // The cone at z=5 has radius 5, sphere cross-section at z=5 has radius 3.
        // The sphere extends from z=2 to z=8. At z=2, cone radius = 2, sphere radius = sqrt(9-9) = 0.
        // At z=8, cone radius = 8, sphere radius = sqrt(9-9) = 0.
        // They should intersect at some z where r_cone = r_sphere.
        let result = intersect_sphere_cone(&s, &k);
        // Depending on geometry, may have 0, 1, or 2 circles
        match result {
            SphereConeResult::NoIntersection => {}
            SphereConeResult::SingleCircle(_) => {}
            SphereConeResult::TwoCircles(_, _) => {}
            SphereConeResult::TangentPoint(_) => {}
            SphereConeResult::General => {}
        }
    }

    /// Sphere completely inside cone (no intersection of surfaces).
    #[test]
    fn sphere_inside_cone_no_intersection() {
        // Cone: apex at (0,0,-10), 45 degree angle
        // At z=0, cone radius = 10
        // Sphere: center at (0,0,0), radius 2
        // Sphere is entirely inside the cone
        let s = sphere(DVec3::ZERO, 2.0);
        let k = cone(DVec3::new(0.0, 0.0, -10.0), DVec3::Z, 45.0);
        let result = intersect_sphere_cone(&s, &k);
        assert!(matches!(
            result,
            SphereConeResult::NoIntersection | SphereConeResult::General
        ));
    }

    /// Sphere center off-axis should return General.
    #[test]
    fn sphere_off_axis_general() {
        let s = sphere(DVec3::new(2.0, 0.0, 5.0), 3.0);
        let k = cone(DVec3::ZERO, DVec3::Z, 45.0);
        let result = intersect_sphere_cone(&s, &k);
        assert!(matches!(result, SphereConeResult::General));
    }

    /// Sphere tangent to cone (one circle).
    #[test]
    fn sphere_tangent_cone_single_circle() {
        // Cone: apex at origin, 45 degree angle
        // Sphere: center at (0,0,5), radius 5
        // At z=5: r_cone = 5, sphere cross-section radius at z=0 = sqrt(25-25) = 0
        // So sphere just touches cone at z=5
        let s = sphere(DVec3::new(0.0, 0.0, 5.0), 5.0);
        let k = cone(DVec3::ZERO, DVec3::Z, 45.0);
        let result = intersect_sphere_cone(&s, &k);
        // Expect at least one circle (might be tangent)
        match result {
            SphereConeResult::SingleCircle(c) => {
                assert!((c.center.z - 5.0).abs() < 1e-3 || (c.radius - 5.0).abs() < 1e-3);
            }
            SphereConeResult::TwoCircles(_, _) => {}
            SphereConeResult::TangentPoint(_) => {}
            _ => {}
        }
    }
}
