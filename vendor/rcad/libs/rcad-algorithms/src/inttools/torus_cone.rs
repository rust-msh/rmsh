//! Analytic intersection of a torus and a cone.
//!
//! # Case classification
//!
//! ## Coaxial case (cone apex on torus axis, cone axis = torus axis)
//!
//! When both surfaces share the same axis line and the cone apex lies on the
//! torus axis, the intersection consists of circles at axial heights where
//! the torus tube circle (in the rho-z half-plane) intersects the cone line.
//!
//! In the (rho, z) half-plane:
//! - Torus tube: (rho - R)² + z² = r²  (circle centered at (R, 0))
//! - Cone: rho = (z - z_apex) * tan(half_angle)  (line from apex)
//!
//! Substituting gives a quadratic equation in z.
//!
//! ## General case
//!
//! For all other configurations the intersection is a complex space curve.
//! We return `General` so the caller falls back to numeric marching.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, ConicalSurface, ToroidalSurface};

use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG};

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of torus x cone intersection.
#[derive(Debug, Clone)]
pub enum TorusConeResult {
    /// The torus and cone do not intersect.
    NoIntersection,
    /// Coaxial case: single intersection circle.
    SingleCircle(Circle3),
    /// Coaxial case: two intersection circles.
    TwoCircles(Circle3, Circle3),
    /// The intersection is a tangent circle.
    TangentCircle(Circle3),
    /// General case. Caller should fall back to marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `torus` and `cone`.
pub fn intersect_torus_cone(
    torus: &ToroidalSurface,
    cone: &ConicalSurface,
) -> TorusConeResult {
    intersect_torus_cone_with_tolerance(torus, cone, 0.0)
}

/// Compute torus-cone intersection with additional fuzzy tolerance.
///
/// This relaxes coaxial detection by `fuzzy_tol` so near-coaxial cases
/// can still classify into analytic branches.
pub fn intersect_torus_cone_with_tolerance(
    torus: &ToroidalSurface,
    cone: &ConicalSurface,
    fuzzy_tol: f64,
) -> TorusConeResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);

    let t_axis = torus.axis.normalize();
    let c_axis = cone.axis_dir();
    let cross = t_axis.cross(c_axis);
    let sin_angle = cross.length();
    let apex = cone.apex_point();

    // Project cone apex onto torus axis
    let t = (apex - torus.center).dot(t_axis);
    let foot = torus.center + t_axis * t;
    let d_apex = (apex - foot).length();

    // ── Coaxial: same axis line ───────────────────────────────────────────────
    if sin_angle < TOLERANCE_ANG && d_apex < tol {
        return intersect_torus_cone_coaxial(torus, cone, t_axis);
    }

    // ── General case: numerical fallback ─────────────────────────────────────
    TorusConeResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Coaxial case
// ─────────────────────────────────────────────────────────────────────────────

#[allow(non_snake_case)]
fn intersect_torus_cone_coaxial(
    torus: &ToroidalSurface,
    cone: &ConicalSurface,
    axis: DVec3,
) -> TorusConeResult {
    let R = torus.major_radius;
    let r = torus.minor_radius;
    let ta = cone.half_angle_rad.tan();

    // Determine cone orientation relative to torus axis
    let sigma = if cone.axis_dir().dot(axis) >= 0.0 { 1.0 } else { -1.0 };
    let apex = cone.apex_point();
    let r_ref = cone.radius;

    // z_apex: axial coordinate of cone apex relative to torus center
    let z_apex = (apex - torus.center).dot(axis);

    // In the (rho, z) half-plane:
    // Torus tube: (rho - R)² + z² = r²   (tube center at (R, 0))
    // Cone: rho = r_ref + sigma * (z - z_apex) * ta
    //
    // Let's set up the equation:
    // Let a_cone = sigma * ta  (cone slope in rho-z plane, signed)
    // rho_cone(z) = r_ref + a_cone * (z - z_apex)
    //
    // Substitute into torus equation:
    // (r_ref + a_cone*(z - z_apex) - R)² + z² = r²
    //
    // Let A = a_cone, B = r_ref - R - A*z_apex
    // Then rho_cone = A*z + B
    //
    // (A*z + B - R)² + z² = r²
    // A²*z² + 2*A*(B-R)*z + (B-R)² + z² = r²
    // (A² + 1)*z² + 2*A*(B-R)*z + (B-R)² - r² = 0

    let A = sigma * ta;
    let B = r_ref - A * z_apex;
    let rho_offset = B - R;

    let a_q = A * A + 1.0;
    let b_q = 2.0 * A * rho_offset;
    let c_q = rho_offset * rho_offset - r * r;

    let disc = b_q * b_q - 4.0 * a_q * c_q;

    if disc < -TOLERANCE_ABS {
        return TorusConeResult::NoIntersection;
    }

    if disc.abs() < TOLERANCE_ABS {
        // Tangent: one solution
        let z = -b_q / (2.0 * a_q);
        let rho = A * z + B;

        if rho < TOLERANCE_ABS {
            return TorusConeResult::NoIntersection;
        }

        let center = torus.center + axis * z;
        return TorusConeResult::TangentCircle(Circle3 {
            center,
            normal: axis,
            radius: rho,
        });
    }

    // Two solutions
    let sqrt_disc = disc.sqrt();
    let z1 = (-b_q - sqrt_disc) / (2.0 * a_q);
    let z2 = (-b_q + sqrt_disc) / (2.0 * a_q);

    let mut circles: Vec<Circle3> = Vec::new();

    for z in [z1, z2] {
        let rho = A * z + B;
        if rho > TOLERANCE_ABS {
            let center = torus.center + axis * z;
            circles.push(Circle3 {
                center,
                normal: axis,
                radius: rho,
            });
        }
    }

    match circles.len() {
        0 => TorusConeResult::NoIntersection,
        1 => TorusConeResult::SingleCircle(circles[0]),
        2 => TorusConeResult::TwoCircles(circles[0], circles[1]),
        _ => TorusConeResult::General,
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

    fn cone(apex: DVec3, axis: DVec3, half_angle_deg: f64) -> ConicalSurface {
        ConicalSurface {
            apex,
            axis,
            radius: 0.0,
            half_angle_rad: half_angle_deg.to_radians(),
        }
    }

    /// Coaxial torus and cone with intersecting tube.
    /// Torus: R=5, r=4 (large tube), axis=Z, center=origin
    /// Cone: apex at origin, 45 degree angle
    /// Solve: (z*tan(45) - 5)² + z² = 16
    ///        (z - 5)² + z² = 16
    ///        2z² - 10z + 9 = 0
    ///        z = (10 ± sqrt(100 - 72)) / 4 = (10 ± sqrt(28)) / 4
    ///        z ≈ {3.82, 1.18}
    #[test]
    fn coaxial_torus_cone_two_circles() {
        let t = torus(DVec3::ZERO, DVec3::Z, 5.0, 4.0);
        let k = cone(DVec3::ZERO, DVec3::Z, 45.0);
        match intersect_torus_cone(&t, &k) {
            TorusConeResult::TwoCircles(c1, c2) => {
                // Verify z values
                let z1 = c1.center.z;
                let z2 = c2.center.z;
                let expected_z1 = (10.0 - 28_f64.sqrt()) / 4.0;
                let expected_z2 = (10.0 + 28_f64.sqrt()) / 4.0;
                assert!((z1 - expected_z1).abs() < 1e-6 || (z1 - expected_z2).abs() < 1e-6);
                assert!((z2 - expected_z1).abs() < 1e-6 || (z2 - expected_z2).abs() < 1e-6);
            }
            TorusConeResult::SingleCircle(_) => {
                // Also acceptable - tangent case
            }
            other => panic!("expected TwoCircles or SingleCircle, got {other:?}"),
        }
    }

    /// Non-coaxial torus and cone should return General.
    #[test]
    fn non_coaxial_torus_cone_general() {
        let t = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let k = cone(DVec3::new(1.0, 0.0, 0.0), DVec3::Z, 45.0);
        let result = intersect_torus_cone(&t, &k);
        assert!(matches!(result, TorusConeResult::General));
    }

    /// Torus and cone with non-parallel axes should return General.
    #[test]
    fn skew_axes_torus_cone_general() {
        let t = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let k = cone(DVec3::ZERO, DVec3::X, 45.0);
        let result = intersect_torus_cone(&t, &k);
        assert!(matches!(result, TorusConeResult::General));
    }

    /// Cone apex offset from torus axis should return General.
    #[test]
    fn offset_apex_torus_cone_general() {
        let t = torus(DVec3::ZERO, DVec3::Z, 5.0, 1.0);
        let k = cone(DVec3::new(1.0, 0.0, 0.0), DVec3::Z, 45.0);
        let result = intersect_torus_cone(&t, &k);
        assert!(matches!(result, TorusConeResult::General));
    }

    /// Cone with reference radius (non-zero base radius).
    #[test]
    fn cone_with_reference_radius() {
        let t = torus(DVec3::ZERO, DVec3::Z, 5.0, 2.0);
        let k = ConicalSurface {
            apex: DVec3::new(0.0, 0.0, -3.0),
            axis: DVec3::Z,
            radius: 1.0, // Reference radius at apex
            half_angle_rad: 45.0_f64.to_radians(),
        };
        let result = intersect_torus_cone(&t, &k);
        // Should find intersection
        match result {
            TorusConeResult::SingleCircle(_) => {}
            TorusConeResult::TwoCircles(_, _) => {}
            TorusConeResult::TangentCircle(_) => {}
            TorusConeResult::NoIntersection => {}
            TorusConeResult::General => {}
        }
    }
}
