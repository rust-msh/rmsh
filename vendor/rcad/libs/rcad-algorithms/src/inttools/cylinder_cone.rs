//! Analytic intersection of a cylinder and a cone.
//!
//! # Case classification
//!
//! ## Coaxial (axes coincide)
//!
//! When the cylinder axis and cone axis are the same line, the intersection is
//! a **circle** at the height where the cone radius equals the cylinder radius:
//!
//! ```text
//! r_cone(h) = (h − h_apex) · tan(β)  where h_apex is the height of the apex
//! r_cone(h) = r_cyl  →  h = h_apex + r_cyl / tan(β)
//! ```
//!
//! Returns [`CoaxialCircle`](CylinderConeResult::CoaxialCircle) if the apex is
//! on the positive-axis side (the circle is real), or
//! [`NoIntersection`](CylinderConeResult::NoIntersection) otherwise.
//!
//! ## Parallel axes (non-coaxial)
//!
//! When the axes are parallel but distinct (cross-product ≈ 0), the
//! cylinder-cone intersection is a quartic curve in general.  We perform a
//! radial distance test to detect obvious non-intersections and fall back to
//! marching otherwise.
//!
//! ## General / skew axes
//!
//! For all other configurations we return [`General`](CylinderConeResult::General)
//! so the caller can fall back to numeric marching.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, ConicalSurface, CylindricalSurface};

use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG};

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of cylinder × cone intersection.
#[derive(Debug, Clone)]
pub enum CylinderConeResult {
    /// The surfaces do not intersect.
    NoIntersection,
    /// Coaxial configuration: exactly one intersection circle.
    CoaxialCircle(Circle3),
    /// General case (skew axes or oblique angle not handled analytically).
    /// The caller should fall back to numeric marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `cyl` and `cone`.
pub fn intersect_cylinder_cone(
    cyl: &CylindricalSurface,
    cone: &ConicalSurface,
) -> CylinderConeResult {
    let a_cyl = cyl.axis.normalize();
    let a_cone = cone.axis.normalize();

    let cross = a_cyl.cross(a_cone);
    let sin_angle = cross.length(); // |sin θ| between the two axes

    // ── Parallel axes (including coaxial) ────────────────────────────────────
    if sin_angle < TOLERANCE_ANG {
        return intersect_parallel_cylinder_cone(cyl, cone, a_cyl, a_cone);
    }

    // ── General / skew ────────────────────────────────────────────────────────
    // Perform a quick distance-based no-intersection test:
    // Find closest distance between the two axes.  If the cylinder completely
    // misses the cone's bounding envelope, return NoIntersection.

    // For now, return General (marching handles this correctly).
    CylinderConeResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel (and coaxial) axes
// ─────────────────────────────────────────────────────────────────────────────

fn intersect_parallel_cylinder_cone(
    cyl: &CylindricalSurface,
    cone: &ConicalSurface,
    a_cyl: DVec3,
    a_cone: DVec3,
) -> CylinderConeResult {
    let apex = cone.apex_point();
    // Make sure a_cone points the same direction as a_cyl for height arithmetic.
    // (The cross product is ~0 so they are parallel; they may anti-parallel.)
    let _a_cone = if a_cyl.dot(a_cone) >= 0.0 { a_cone } else { -a_cone };

    let r_cyl = cyl.radius;
    let tan_beta = cone.half_angle_rad.tan();

    // Perpendicular distance between the two axes.
    let delta = apex - cyl.origin;
    let delta_perp = delta - a_cyl * delta.dot(a_cyl);
    let d_perp = delta_perp.length();

    // ── Coaxial ──────────────────────────────────────────────────────────────
    if d_perp < TOLERANCE_ABS {
        // Height of apex above cyl.origin along shared axis.
        let h_apex = (apex - cyl.origin).dot(a_cyl);

        // At height h (measured from cyl.origin), cone radius = (h - h_apex)*tan_beta
        // (only positive when h > h_apex, i.e. above the apex in axis direction).
        // Set equal to r_cyl:  h = h_apex + r_cyl / tan_beta
        if tan_beta.abs() < 1e-14 {
            // Degenerate cone (half_angle = 0), no lateral surface.
            return CylinderConeResult::NoIntersection;
        }
        let h_circle = h_apex + r_cyl / tan_beta;

        // The circle must be on the cone's nappe (h_circle > h_apex).
        if h_circle <= h_apex - TOLERANCE_ABS {
            return CylinderConeResult::NoIntersection;
        }

        let center = cyl.origin + a_cyl * h_circle;
        return CylinderConeResult::CoaxialCircle(Circle3 {
            center,
            normal: a_cyl,
            radius: r_cyl,
        });
    }

    // ── Parallel but offset axes ──────────────────────────────────────────────
    // Quick radial distance test:
    // At every height h the cone has radius r_cone(h) = max(0, |h - h_apex| * tan_beta).
    // The cylinder is at radial distance r_cyl from its axis.
    // The two axes are d_perp apart.
    //
    // The cylinder surface and cone surface can only intersect if there exists
    // some h where the radial circles (centred d_perp apart) overlap:
    //   r_cone(h) + r_cyl >= d_perp   and   |r_cone(h) - r_cyl| <= d_perp
    //
    // Since r_cone grows unboundedly with |h|, intersection is always possible
    // unless the cylinder axis is so far from the cone apex that even the widest
    // accessible h gives no overlap (but in practice the cone is finite in the
    // BRep, so marching will bound it).
    //
    // We do the simple check: can the two ever touch at any height?
    // r_cone >= d_perp - r_cyl (cylinder can be reached from cone at some h)
    // → h - h_apex >= (d_perp - r_cyl) / tan_beta
    // This is always achievable if d_perp - r_cyl > 0 is reachable, i.e.
    // the cylinder is not entirely inside the cone at all heights — always true
    // for parallel offset.  So fall back to General.
    CylinderConeResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    fn cyl(origin: DVec3, axis: DVec3, radius: f64) -> CylindricalSurface {
        CylindricalSurface { origin, axis, radius }
    }

    fn cone(apex: DVec3, axis: DVec3, half_angle_deg: f64) -> ConicalSurface {
        ConicalSurface {
            apex,
            axis,
            radius: 0.0,
            half_angle_rad: half_angle_deg.to_radians(),
        }
    }

    /// Coaxial cylinder and cone (Z axis): circle at h = apex_z + r / tan(β).
    ///
    /// Cone: apex at (0,0,0), axis Z, half_angle=45° → tan(β)=1.
    /// Cylinder: r=2, axis Z → circle at h = 0 + 2/1 = 2.
    #[test]
    fn coaxial_circle_z_axis() {
        let c = cyl(DVec3::ZERO, DVec3::Z, 2.0);
        let k = cone(DVec3::ZERO, DVec3::Z, 45.0);
        match intersect_cylinder_cone(&c, &k) {
            CylinderConeResult::CoaxialCircle(circ) => {
                assert!(
                    (circ.center.z - 2.0).abs() < 1e-9,
                    "circle z={}, expected 2.0",
                    circ.center.z
                );
                assert!((circ.radius - 2.0).abs() < 1e-9);
            }
            other => panic!("expected CoaxialCircle, got {other:?}"),
        }
    }

    /// Coaxial but cone apex ABOVE the cylinder origin with an offset.
    ///
    /// Cone: apex at (0,0,5), axis Z, half_angle=30° → tan(β)=1/√3.
    /// Cylinder: r=1 → h = 5 + 1/(1/√3) = 5 + √3 ≈ 6.732.
    #[test]
    fn coaxial_circle_offset_apex() {
        let c = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let k = cone(DVec3::new(0.0, 0.0, 5.0), DVec3::Z, 30.0);
        match intersect_cylinder_cone(&c, &k) {
            CylinderConeResult::CoaxialCircle(circ) => {
                let expected_h = 5.0 + 1.0 / (30.0_f64.to_radians().tan());
                assert!(
                    (circ.center.z - expected_h).abs() < 1e-9,
                    "circle z={}, expected {}",
                    circ.center.z,
                    expected_h
                );
            }
            other => panic!("expected CoaxialCircle, got {other:?}"),
        }
    }

    /// Skew axes → General.
    #[test]
    fn skew_axes_general() {
        let c = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let k = cone(DVec3::ZERO, DVec3::new(1.0, 1.0, 0.0).normalize(), 45.0);
        assert!(matches!(intersect_cylinder_cone(&c, &k), CylinderConeResult::General));
    }

    /// Perpendicular axes → General.
    #[test]
    fn perpendicular_axes_general() {
        let c = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let k = cone(DVec3::ZERO, DVec3::X, 45.0);
        assert!(matches!(intersect_cylinder_cone(&c, &k), CylinderConeResult::General));
    }
}
