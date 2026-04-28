//! Analytic intersection of a cylinder with a torus.
//!
//! # Cases
//!
//! - **Coaxial**: Two circles when cylinder radius intersects torus tube
//! - **Perpendicular axes**: Complex quartic curve, numerical fallback
//! - **General**: Numerical marching

use rcad_kernel::geom::{Circle3, CylindricalSurface, ToroidalSurface};

use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG};

/// Result of cylinder x torus intersection.
#[derive(Debug, Clone)]
pub enum CylinderTorusResult {
    /// No intersection.
    NoIntersection,
    /// Single tangent circle.
    TangentCircle(Circle3),
    /// Two circles (coaxial case).
    TwoCircles(Circle3, Circle3),
    /// Complex intersection, fall back to numerical marching.
    General,
}

/// Compute the analytic intersection of `cylinder` and `torus`.
pub fn intersect_cylinder_torus(
    cylinder: &CylindricalSurface,
    torus: &ToroidalSurface,
) -> CylinderTorusResult {
    intersect_cylinder_torus_with_tolerance(cylinder, torus, 0.0)
}

/// Cylinder x torus intersection with fuzzy tolerance.
pub fn intersect_cylinder_torus_with_tolerance(
    cylinder: &CylindricalSurface,
    torus: &ToroidalSurface,
    fuzzy_tol: f64,
) -> CylinderTorusResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);

    let a_cyl = cylinder.axis.normalize();
    let a_tor = torus.axis.normalize();

    // Check for coaxial case
    let cross = a_cyl.cross(a_tor);
    let sin_angle = cross.length();

    // Project cylinder origin onto torus axis
    let delta = cylinder.origin - torus.center;
    let d_perp = (delta - a_tor * delta.dot(a_tor)).length();

    // Coaxial: same axis line
    if sin_angle < TOLERANCE_ANG && d_perp < tol {
        return intersect_cylinder_torus_coaxial(cylinder, torus, tol);
    }

    // General case: numerical fallback
    CylinderTorusResult::General
}

fn intersect_cylinder_torus_coaxial(
    cylinder: &CylindricalSurface,
    torus: &ToroidalSurface,
    tol: f64,
) -> CylinderTorusResult {
    // In the (rho, z) half-plane:
    // Torus tube: (rho - major_r)^2 + z^2 = minor_r^2
    // Cylinder: rho = r_cyl

    let major_r = torus.major_radius;
    let minor_r = torus.minor_radius;
    let r_cyl = cylinder.radius;

    // Cylinder radius must intersect the tube circle
    // Tube circle center at (major_r, 0) in (rho, z) plane
    // Distance from tube center to cylinder: |major_r - r_cyl|
    let d = (major_r - r_cyl).abs();

    if d > minor_r + tol {
        return CylinderTorusResult::NoIntersection;
    }

    if (d - minor_r).abs() < tol {
        // Tangent: one circle
        let z = 0.0;
        let center = torus.center + torus.axis * z;
        return CylinderTorusResult::TangentCircle(Circle3 {
            center,
            normal: torus.axis,
            radius: r_cyl,
        });
    }

    // Two intersection points in (rho, z): z = +/-sqrt(minor_r^2 - d^2)
    let z_offset = (minor_r * minor_r - d * d).sqrt();

    let center1 = torus.center + torus.axis * z_offset;
    let center2 = torus.center - torus.axis * z_offset;

    CylinderTorusResult::TwoCircles(
        Circle3 { center: center1, normal: torus.axis, radius: r_cyl },
        Circle3 { center: center2, normal: torus.axis, radius: r_cyl },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{CylindricalSurface, ToroidalSurface};

    #[test]
    fn coaxial_cylinder_torus_produces_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Cylinder coaxial with torus, radius between inner and outer
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 5.0, // Same as major radius
        };

        let result = intersect_cylinder_torus(&cylinder, &torus);

        // Should produce two circles at z = +/-minor_radius
        match result {
            CylinderTorusResult::TwoCircles(c1, c2) => {
                assert!((c1.radius - 5.0).abs() < 1e-6, "c1 radius expected 5.0, got {}", c1.radius);
                assert!((c2.radius - 5.0).abs() < 1e-6, "c2 radius expected 5.0, got {}", c2.radius);
                // Circle centers should be at +/-minor_radius along the axis
                assert!((c1.center - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-6, "c1 center expected (0,1,0)");
                assert!((c2.center - DVec3::new(0.0, -1.0, 0.0)).length() < 1e-6, "c2 center expected (0,-1,0)");
            }
            other => panic!("Expected TwoCircles for coaxial case, got {:?}", other),
        }
    }

    #[test]
    fn coaxial_cylinder_torus_no_intersection() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Cylinder coaxial but radius outside torus tube range
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 7.0, // Outside R + r = 6.0
        };

        let result = intersect_cylinder_torus(&cylinder, &torus);
        assert!(matches!(result, CylinderTorusResult::NoIntersection));
    }

    #[test]
    fn coaxial_cylinder_torus_tangent() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Cylinder coaxial, tangent to outer edge of tube
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 6.0, // Exactly R + r = 6.0
        };

        let result = intersect_cylinder_torus(&cylinder, &torus);

        match result {
            CylinderTorusResult::TangentCircle(c) => {
                assert!((c.radius - 6.0).abs() < 1e-6);
                // Center should be at torus center (z=0)
                assert!((c.center - DVec3::ZERO).length() < 1e-6);
            }
            other => panic!("Expected TangentCircle, got {:?}", other),
        }
    }

    #[test]
    fn coaxial_cylinder_torus_inner_tangent() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Cylinder coaxial, tangent to inner edge of tube
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 4.0, // Exactly R - r = 4.0
        };

        let result = intersect_cylinder_torus(&cylinder, &torus);

        match result {
            CylinderTorusResult::TangentCircle(c) => {
                assert!((c.radius - 4.0).abs() < 1e-6);
                assert!((c.center - DVec3::ZERO).length() < 1e-6);
            }
            other => panic!("Expected TangentCircle, got {:?}", other),
        }
    }

    #[test]
    fn non_coaxial_returns_general() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Cylinder with perpendicular axis
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::X,
            radius: 5.0,
        };

        let result = intersect_cylinder_torus(&cylinder, &torus);
        assert!(matches!(result, CylinderTorusResult::General));
    }

    #[test]
    fn coaxial_with_offset_origin() {
        // Torus centered at origin
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Cylinder with origin offset along axis (but still coaxial)
        let cylinder = CylindricalSurface {
            origin: DVec3::new(0.0, 2.0, 0.0),
            axis: DVec3::Y,
            radius: 5.0,
        };

        let result = intersect_cylinder_torus(&cylinder, &torus);
        // Should still produce two circles (coaxial geometry is the same)
        assert!(matches!(result, CylinderTorusResult::TwoCircles(_, _)));
    }
}
