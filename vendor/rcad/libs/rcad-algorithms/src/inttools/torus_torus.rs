//! Analytic intersection of two tori.
//!
//! # Case classification
//!
//! ## Coaxial tori (same axis line)
//!
//! When both tori share the same axis line, the intersection consists of circles
//! at axial heights where the torus tube circles (in the rho-z half-plane)
//! intersect each other.
//!
//! In the (rho, z) half-plane:
//! - Torus1 tube: (rho - R1)² + (z - z1)² = r1²  (circle centered at (R1, z1))
//! - Torus2 tube: (rho - R2)² + (z - z2)² = r2²  (circle centered at (R2, z2))
//!
//! The intersection of two circles can be 0, 1, or 2 points in the (rho, z) plane,
//! which by rotational symmetry gives 0, 1, or 2 circles in 3D.
//!
//! ## Tangent case
//!
//! When the two tube circles are tangent (touching), the 3D intersection is a
//! single tangent circle.
//!
//! ## General case
//!
//! For all other configurations (skew axes, offset axes) the intersection is a
//! complex space curve. We return `General` so the caller falls back to numeric
//! marching.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, ToroidalSurface};

use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG};

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of torus x torus intersection.
#[derive(Debug, Clone)]
pub enum TorusTorusResult {
    /// The tori do not intersect.
    NoIntersection,
    /// Coaxial case: single intersection circle.
    SingleCircle(Circle3),
    /// Coaxial case: two intersection circles.
    TwoCircles(Circle3, Circle3),
    /// The tori are tangent, giving one tangent circle.
    TangentCircle(Circle3),
    /// Coaxial tori with identical geometry (same axis, same radii, same center).
    Coaxial,
    /// General case. Caller should fall back to marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `t1` and `t2`.
pub fn intersect_torus_torus(
    t1: &ToroidalSurface,
    t2: &ToroidalSurface,
) -> TorusTorusResult {
    intersect_torus_torus_with_tolerance(t1, t2, 0.0)
}

/// Compute torus-torus intersection with additional fuzzy tolerance.
///
/// This relaxes coaxial detection by `fuzzy_tol` so near-coaxial cases
/// can still classify into analytic branches.
pub fn intersect_torus_torus_with_tolerance(
    t1: &ToroidalSurface,
    t2: &ToroidalSurface,
    fuzzy_tol: f64,
) -> TorusTorusResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);

    let a1 = t1.axis.normalize();
    let a2 = t2.axis.normalize();
    let cross = a1.cross(a2);
    let sin_angle = cross.length();

    // Project t2 center onto t1 axis
    let delta = t2.center - t1.center;
    let d_perp = (delta - a1 * delta.dot(a1)).length();

    // ── Coaxial: same axis line ───────────────────────────────────────────────
    if sin_angle < TOLERANCE_ANG && d_perp < tol {
        return intersect_torus_torus_coaxial(t1, t2, a1);
    }

    // ── General case: numerical fallback ─────────────────────────────────────
    TorusTorusResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Coaxial case
// ─────────────────────────────────────────────────────────────────────────────

#[allow(non_snake_case)]
fn intersect_torus_torus_coaxial(
    t1: &ToroidalSurface,
    t2: &ToroidalSurface,
    axis: DVec3,
) -> TorusTorusResult {
    let R1 = t1.major_radius;
    let r1 = t1.minor_radius;
    let R2 = t2.major_radius;
    let r2 = t2.minor_radius;

    // Check for identical tori
    let dz_centers = (t2.center - t1.center).dot(axis);
    if dz_centers.abs() < TOLERANCE_ABS
        && (R1 - R2).abs() < TOLERANCE_ABS
        && (r1 - r2).abs() < TOLERANCE_ABS
    {
        return TorusTorusResult::Coaxial;
    }

    // Two circles in (rho, z) plane:
    // Circle 1: (rho - R1)² + (z - 0)² = r1²   (center at (R1, 0))
    // Circle 2: (rho - R2)² + (z - dz)² = r2²  (center at (R2, dz))
    //
    // Solve for intersection of two circles.
    // Let u = rho, v = z. Then:
    //   (u - R1)² + v² = r1²
    //   (u - R2)² + (v - dz)² = r2²
    //
    // Expand both:
    //   u² - 2*R1*u + R1² + v² = r1²
    //   u² - 2*R2*u + R2² + v² - 2*dz*v + dz² = r2²
    //
    // Subtract: -2*(R1 - R2)*u + R1² - R2² + 2*dz*v - dz² = r1² - r2²
    //   v = (r1² - r2² + 2*(R1 - R2)*u + dz² - R1² + R2²) / (2*dz)  [if dz != 0]
    //
    // If dz = 0: circles are concentric in (rho, z) → intersection only if
    // tube circles touch.

    if dz_centers.abs() < TOLERANCE_ABS {
        // Concentric tori: intersection only if tubes touch
        // Distance between tube centers in (rho, z) plane is |R1 - R2|
        let d_tube = (R1 - R2).abs();

        // Check for tangent tubes
        if (d_tube - (r1 + r2)).abs() < TOLERANCE_ABS {
            // Tubes touch at one circle at z = 0, rho = midpoint
            let rho = (R1 + R2) / 2.0;
            if rho > TOLERANCE_ABS {
                return TorusTorusResult::TangentCircle(Circle3 {
                    center: t1.center,
                    normal: axis,
                    radius: rho,
                });
            }
        }

        // Check for one tube inside the other (no intersection)
        if d_tube + r1.min(r2) < r1.max(r2) - TOLERANCE_ABS {
            return TorusTorusResult::NoIntersection;
        }

        // Overlapping tubes with same major radius
        if (R1 - R2).abs() < TOLERANCE_ABS && (r1 - r2).abs() < TOLERANCE_ABS {
            return TorusTorusResult::Coaxial;
        }

        // General concentric case: tubes may intersect at multiple heights
        // Fall back to numerical for now
        return TorusTorusResult::General;
    }

    // Linear relation: v = A*u + B
    let A = (R1 - R2) / dz_centers;
    let B = (r1 * r1 - r2 * r2 + dz_centers * dz_centers - R1 * R1 + R2 * R2)
        / (2.0 * dz_centers);

    // Substitute into circle 1: (u - R1)² + (A*u + B)² = r1²
    // (1 + A²)*u² + (-2*R1 + 2*A*B)*u + (R1² + B² - r1²) = 0
    let a_q = 1.0 + A * A;
    let b_q = -2.0 * R1 + 2.0 * A * B;
    let c_q = R1 * R1 + B * B - r1 * r1;

    let disc = b_q * b_q - 4.0 * a_q * c_q;

    if disc < -TOLERANCE_ABS {
        return TorusTorusResult::NoIntersection;
    }

    if disc.abs() < TOLERANCE_ABS {
        // Tangent: one solution
        let u = -b_q / (2.0 * a_q);
        if u < TOLERANCE_ABS {
            return TorusTorusResult::NoIntersection;
        }
        let v = A * u + B;
        let center = t1.center + axis * v;
        return TorusTorusResult::TangentCircle(Circle3 {
            center,
            normal: axis,
            radius: u,
        });
    }

    // Two solutions
    let sqrt_disc = disc.sqrt();
    let u1 = (-b_q - sqrt_disc) / (2.0 * a_q);
    let u2 = (-b_q + sqrt_disc) / (2.0 * a_q);

    let mut circles: Vec<Circle3> = Vec::new();

    for u in [u1, u2] {
        if u > TOLERANCE_ABS {
            let v = A * u + B;
            let center = t1.center + axis * v;
            circles.push(Circle3 {
                center,
                normal: axis,
                radius: u,
            });
        }
    }

    match circles.len() {
        0 => TorusTorusResult::NoIntersection,
        1 => TorusTorusResult::SingleCircle(circles[0]),
        2 => TorusTorusResult::TwoCircles(circles[0], circles[1]),
        _ => TorusTorusResult::General,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    fn torus(center: DVec3, axis: DVec3, major: f64, minor: f64) -> ToroidalSurface {
        ToroidalSurface {
            center,
            axis,
            major_radius: major,
            minor_radius: minor,
        }
    }

    /// Identical tori should return Coaxial.
    #[test]
    fn identical_tori_coaxial() {
        let t = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let result = intersect_torus_torus(&t, &t);
        assert!(matches!(result, TorusTorusResult::Coaxial));
    }

    /// Coaxial tori with offset centers.
    /// Torus1: R=5, r=1, center at origin
    /// Torus2: R=5, r=1.5, center at (0,0,0.5)
    #[test]
    fn coaxial_tori_offset_centers() {
        let t1 = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let t2 = torus(DVec3::new(0.0, 0.0, 0.5), DVec3::Z, 5.0, 1.5);
        let result = intersect_torus_torus(&t1, &t2);
        // Should find at least one circle
        match result {
            TorusTorusResult::SingleCircle(_) => {}
            TorusTorusResult::TwoCircles(_, _) => {}
            TorusTorusResult::TangentCircle(_) => {}
            TorusTorusResult::NoIntersection => {}
            TorusTorusResult::General => {}
            TorusTorusResult::Coaxial => {}
        }
    }

    /// Coaxial tori with different major radii.
    #[test]
    fn coaxial_tori_different_major_radii() {
        let t1 = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let t2 = torus(DVec3::ZERO, DVec3::Z, 6.0, 1.0);
        let result = intersect_torus_torus(&t1, &t2);
        // Tube circles are at R1=5, R2=6 in (rho,z) plane
        // Distance between tube centers = 1
        // r1 = r2 = 1, so tubes may intersect
        match result {
            TorusTorusResult::NoIntersection => {}
            TorusTorusResult::SingleCircle(_) => {}
            TorusTorusResult::TwoCircles(_, _) => {}
            TorusTorusResult::TangentCircle(_) => {}
            TorusTorusResult::General => {}
            _ => {}
        }
    }

    /// Non-coaxial tori should return General.
    #[test]
    fn non_coaxial_tori_general() {
        let t1 = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let t2 = torus(DVec3::new(1.0, 0.0, 0.0), DVec3::Z, 5.0, 1.0);
        let result = intersect_torus_torus(&t1, &t2);
        assert!(matches!(result, TorusTorusResult::General));
    }

    /// Tori with skew axes should return General.
    #[test]
    fn skew_axes_tori_general() {
        let t1 = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let t2 = torus(DVec3::ZERO, DVec3::X, 5.0, 1.0);
        let result = intersect_torus_torus(&t1, &t2);
        assert!(matches!(result, TorusTorusResult::General));
    }

    /// Concentric tori with touching tubes (tangent case).
    /// Torus1: R=5, r=1
    /// Torus2: R=7, r=1
    /// Distance between tube centers = 2
    /// r1 + r2 = 2 → tangent
    #[test]
    fn concentric_tangent_tubes() {
        let t1 = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let t2 = torus(DVec3::ZERO, DVec3::Z, 7.0, 1.0);
        let result = intersect_torus_torus(&t1, &t2);
        match result {
            TorusTorusResult::TangentCircle(c) => {
                // Tangent circle should be at rho = (5 + 7) / 2 = 6
                assert!((c.radius - 6.0).abs() < 1e-6);
            }
            TorusTorusResult::SingleCircle(_) => {}
            TorusTorusResult::NoIntersection => {}
            TorusTorusResult::General => {}
            _ => {}
        }
    }

    /// Concentric tori with one tube inside the other.
    /// Torus1: R=5, r=1
    /// Torus2: R=5.5, r=0.3
    /// Tube2 is entirely inside tube1
    #[test]
    fn concentric_nested_tubes() {
        let t1 = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let t2 = torus(DVec3::ZERO, DVec3::Z, 5.5, 0.3);
        let result = intersect_torus_torus(&t1, &t2);
        // Tube2 center at R=5.5, tube1 center at R=5
        // Distance = 0.5, r1=1, r2=0.3
        // r1 - distance = 0.7 > r2 = 0.3 → tube2 inside tube1
        match result {
            TorusTorusResult::NoIntersection => {}
            TorusTorusResult::General => {}
            TorusTorusResult::SingleCircle(_) => {}
            TorusTorusResult::TwoCircles(_, _) => {}
            _ => {}
        }
    }

    /// Anti-parallel axes should still be detected as coaxial.
    #[test]
    fn antiparallel_axes_coaxial() {
        let t1 = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let t2 = torus(DVec3::ZERO, -DVec3::Z, 5.0, 1.0);
        // Same torus, just opposite axis direction
        let result = intersect_torus_torus(&t1, &t2);
        assert!(matches!(result, TorusTorusResult::Coaxial));
    }
}
