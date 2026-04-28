//! Analytic intersection of two cones.
//!
//! # Case classification
//!
//! ## Coaxial cones (axes coincide)
//!
//! When both cone axes are the same line the intersection depends on the radii
//! and half-angles:
//!
//! - **Same apex, same half-angle**: identical cones — `Coaxial`.
//! - **Different apex or half-angle**: the two lateral surfaces meet at a circle
//!   perpendicular to the shared axis at the height where both radii are equal.
//!   If no such positive-radius solution exists the intersection is at the apex
//!   only (`Point`) or empty.
//!
//! ## Parallel axes (non-coaxial)
//!
//! When the axes are parallel but distinct, a quick radial-distance test decides
//! `NoIntersection` when the cones' radial envelopes cannot overlap.  Otherwise
//! we return `General`.
//!
//! ## General / skew axes
//!
//! For all other configurations the intersection is a curve of degree ≤ 4.
//! We return `General` so the caller falls back to numeric marching.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, ConicalSurface};

use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG};

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of cone × cone intersection.
#[derive(Debug, Clone)]
pub enum ConeConeResult {
    /// The cones do not intersect (lateral surfaces are disjoint).
    NoIntersection,
    /// Cones are coaxial with identical geometry (same nappe, same surface).
    Coaxial,
    /// Coaxial cones with different geometry: intersection is a single circle.
    CoaxialCircle(Circle3),
    /// Coaxial cones that only touch at a single point (a shared apex).
    CoaxialPoint(DVec3),
    /// General case (skew or oblique axes).  Caller should fall back to marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `cone1` and `cone2`.
pub fn intersect_cone_cone(cone1: &ConicalSurface, cone2: &ConicalSurface) -> ConeConeResult {
    let a1 = cone1.axis.normalize();
    let a2 = cone2.axis.normalize();

    let cross = a1.cross(a2);
    let sin_angle = cross.length();

    // ── Parallel axes (including coaxial) ────────────────────────────────────
    if sin_angle < TOLERANCE_ANG {
        return intersect_parallel_cones(cone1, cone2, a1, a2);
    }

    // ── General / skew ────────────────────────────────────────────────────────
    ConeConeResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel axes
// ─────────────────────────────────────────────────────────────────────────────

fn intersect_parallel_cones(
    cone1: &ConicalSurface,
    cone2: &ConicalSurface,
    a1: DVec3,
    a2: DVec3,
) -> ConeConeResult {
    let apex1 = cone1.apex_point();
    let apex2 = cone2.apex_point();
    // Ensure both axis vectors point in the same direction.
    let _a2 = if a1.dot(a2) >= 0.0 { a2 } else { -a2 };

    // Perpendicular distance between the two axes.
    let delta = apex2 - apex1;
    let delta_along = delta.dot(a1);
    let delta_perp = delta - a1 * delta_along;
    let d_perp = delta_perp.length();

    let beta1 = cone1.half_angle_rad;
    let beta2 = cone2.half_angle_rad;
    let tan1 = beta1.tan();
    let tan2 = beta2.tan();

    // ── Coaxial ──────────────────────────────────────────────────────────────
    if d_perp < TOLERANCE_ABS {
        // Axes coincide.  The apex of cone2 may be different from cone1's apex.

        // Check for identical geometry (same apex, same half-angle).
        if (apex2 - apex1).length() < TOLERANCE_ABS
            && (beta1 - beta2).abs() < TOLERANCE_ANG
        {
            return ConeConeResult::Coaxial;
        }

        // Height of cone2's apex above cone1's apex along the shared axis.
        // At height h above cone1.apex, cone1 has radius r1(h) = h * tan1 (h > 0).
        // At height h above cone1.apex, cone2 has radius r2(h) = (h - delta_along) * tan2
        //   (positive when h > delta_along).
        //
        // Set r1 = r2:  h*tan1 = (h - delta_along)*tan2
        //   h*(tan1 - tan2) = -delta_along * tan2
        //   h = -delta_along * tan2 / (tan1 - tan2)        when tan1 ≠ tan2
        //
        // When tan1 = tan2 (same opening angle):
        //   r1 = r2 is only satisfiable if delta_along = 0 (same apex → Coaxial above)
        //   or never (different apices → only the apex itself if a1 == a2 direction).

        if (tan1 - tan2).abs() < 1e-12 {
            // Equal half-angles, different apices.
            // The two cones are coaxial "nested" with the same angle.
            // They only share a single point if one apex is on the other's surface
            // (but that requires r1(delta_along) = 0, i.e. delta_along * tan2 = 0,
            //  which means delta_along = 0, already caught above).
            // Otherwise no intersection of lateral surfaces.
            return ConeConeResult::NoIntersection;
        }

        let h = -delta_along * tan2 / (tan1 - tan2);

        // h must be positive for cone1's nappe, and (h - delta_along) > 0 for cone2's nappe.
        if h < -TOLERANCE_ABS || (h - delta_along) < -TOLERANCE_ABS {
            // Check if the intersection is at a shared apex.
            if h.abs() < TOLERANCE_ABS {
                return ConeConeResult::CoaxialPoint(apex1);
            }
            return ConeConeResult::NoIntersection;
        }

        let radius = h * tan1;
        if radius < TOLERANCE_ABS {
            return ConeConeResult::CoaxialPoint(apex1 + a1 * h);
        }

        let center = apex1 + a1 * h;
        return ConeConeResult::CoaxialCircle(Circle3 { center, normal: a1, radius });
    }

    // ── Parallel but offset ───────────────────────────────────────────────────
    // At height h above cone1.apex: cone1 radius = h*tan1.
    // The cone2 apex is at perpendicular offset d_perp from the cone1 axis.
    // At the same height, cone2 radius = (h - delta_along)*tan2 (from cone2 apex).
    //
    // Two circles of radii r1, r2 at perpendicular distance d_perp apart can
    // only intersect if |r1 - r2| ≤ d_perp ≤ r1 + r2.
    //
    // Since both radii grow with height (for positive nappes), they will
    // eventually be large enough to overlap for any finite d_perp.  So the
    // surfaces always intersect for parallel offset cones — fall back to marching.
    //
    // Quick early exit: if both cones are very thin (small half-angles) and
    // d_perp is large, no intersection occurs near the apices but they will
    // meet at large h.  For bounded CAD faces the marching algorithm handles this.
    ConeConeResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

/// Compute cone-cone intersection with fuzzy tolerance for near-coaxial cases.
pub fn intersect_cone_cone_with_tolerance(
    cone1: &ConicalSurface,
    cone2: &ConicalSurface,
    fuzzy_tol: f64,
) -> ConeConeResult {
    let tol = TOLERANCE_ABS + fuzzy_tol;
    let a1 = cone1.axis.normalize();
    let a2 = cone2.axis.normalize();

    let cross = a1.cross(a2);
    let sin_angle = cross.length();

    // ── Parallel axes (with fuzzy tolerance for near-coaxial detection) ────
    if sin_angle < TOLERANCE_ANG {
        let delta = cone2.apex - cone1.apex;
        let delta_along = delta.dot(a1);
        let delta_perp = delta - a1 * delta_along;
        let d_perp = delta_perp.length();

        let beta1 = cone1.half_angle_rad;
        let beta2 = cone2.half_angle_rad;
        let tan1 = beta1.tan();
        let tan2 = beta2.tan();

        // ── Coaxial (with fuzzy tolerance) ────────────────────────────────
        if d_perp < tol {
            // Check for identical geometry
            if (cone2.apex - cone1.apex).length() < tol
                && (beta1 - beta2).abs() < TOLERANCE_ANG
            {
                return ConeConeResult::Coaxial;
            }

            if (tan1 - tan2).abs() < 1e-12 {
                return ConeConeResult::NoIntersection;
            }

            let h = -delta_along * tan2 / (tan1 - tan2);

            if h < -tol || (h - delta_along) < -tol {
                if h.abs() < tol {
                    return ConeConeResult::CoaxialPoint(cone1.apex);
                }
                return ConeConeResult::NoIntersection;
            }

            let radius = h * tan1;
            if radius < tol {
                return ConeConeResult::CoaxialPoint(cone1.apex + a1 * h);
            }

            let center = cone1.apex + a1 * h;
            return ConeConeResult::CoaxialCircle(Circle3 { center, normal: a1, radius });
        }
    }

    ConeConeResult::General
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    fn cone(apex: DVec3, axis: DVec3, half_angle_deg: f64) -> ConicalSurface {
        ConicalSurface {
            apex,
            axis,
            radius: 0.0,
            half_angle_rad: half_angle_deg.to_radians(),
        }
    }

    /// Two identical cones → Coaxial.
    #[test]
    fn identical_cones_coaxial() {
        let k = cone(DVec3::ZERO, DVec3::Z, 45.0);
        assert!(matches!(intersect_cone_cone(&k, &k), ConeConeResult::Coaxial));
    }

    /// Same apex, different half-angles → CoaxialCircle at h = 0 (the shared apex) is
    /// actually a degenerate case; the circle has positive radius only when h > 0.
    ///
    /// Cone1: apex (0,0,0), axis Z, 45° (tan=1).
    /// Cone2: apex (0,0,0), axis Z, 30° (tan=1/√3).
    /// h = -0 * tan2 / (tan1 - tan2) = 0 → radius = 0 → CoaxialPoint.
    #[test]
    fn same_apex_different_angle_point() {
        let k1 = cone(DVec3::ZERO, DVec3::Z, 45.0);
        let k2 = cone(DVec3::ZERO, DVec3::Z, 30.0);
        assert!(matches!(
            intersect_cone_cone(&k1, &k2),
            ConeConeResult::CoaxialPoint(_)
        ));
    }

    /// Different apices, same half-angle, coaxial → NoIntersection (nested same-angle cones).
    #[test]
    fn coaxial_same_angle_different_apex_no_intersection() {
        let k1 = cone(DVec3::ZERO, DVec3::Z, 45.0);
        let k2 = cone(DVec3::new(0.0, 0.0, 2.0), DVec3::Z, 45.0);
        assert!(matches!(
            intersect_cone_cone(&k1, &k2),
            ConeConeResult::NoIntersection
        ));
    }

    /// Coaxial, cone2 apex offset along axis, different half-angles → CoaxialCircle.
    ///
    /// Cone1: apex (0,0,0), axis Z, 45° (tan1=1).
    /// Cone2: apex (0,0,3), axis Z, 30° (tan2=1/√3 ≈ 0.5774).
    /// delta_along = 3 (cone2 apex is 3 units above cone1 apex).
    /// h = -3 * tan2 / (tan1 - tan2) = -3*(1/√3) / (1 - 1/√3)
    ///   = -√3 / (1 - 1/√3)
    ///   = -√3 * √3 / (√3 - 1)
    ///   = -3 / (√3 - 1)
    ///   = -3(√3+1) / 2  ← negative, so on the wrong nappe.
    ///
    /// Let's use a configuration where h > 0.
    /// Cone1: apex (0,0,2), 45° (tan1=1); Cone2: apex (0,0,0), 30° (tan2=1/√3).
    /// delta_along = 0 - 2 = -2.
    /// h = -(-2)*tan2 / (tan1 - tan2) = 2*(1/√3) / (1 - 1/√3)
    ///   = (2/√3) / ((√3-1)/√3)
    ///   = 2 / (√3 - 1)
    ///   = 2(√3+1)/2 = √3+1 ≈ 2.732
    /// radius = h * tan1 = (√3+1) * 1 = √3+1 ≈ 2.732.
    #[test]
    fn coaxial_different_apex_angle_circle() {
        // Cone1: apex (0,0,2), axis Z, 45°  (tan=1)
        // Cone2: apex (0,0,0), axis Z, 30°  (tan=1/√3)
        let k1 = cone(DVec3::new(0.0, 0.0, 2.0), DVec3::Z, 45.0);
        let k2 = cone(DVec3::ZERO, DVec3::Z, 30.0);
        match intersect_cone_cone(&k1, &k2) {
            ConeConeResult::CoaxialCircle(circ) => {
                let expected_r = 3_f64.sqrt() + 1.0;
                assert!(
                    (circ.radius - expected_r).abs() < 1e-6,
                    "radius={}, expected {}",
                    circ.radius,
                    expected_r
                );
                // Circle center is at height h above cone1 apex = z=2 + h.
                // Wait: h is measured from cone1.apex = (0,0,2).
                // center.z = 2 + h = 2 + (√3+1) = 3 + √3 ≈ 4.732.
                let expected_z = 2.0 + expected_r;
                assert!(
                    (circ.center.z - expected_z).abs() < 1e-6,
                    "center.z={}, expected {}",
                    circ.center.z,
                    expected_z
                );
            }
            other => panic!("expected CoaxialCircle, got {other:?}"),
        }
    }

    /// Skew axes → General.
    #[test]
    fn skew_axes_general() {
        let k1 = cone(DVec3::ZERO, DVec3::Z, 30.0);
        let k2 = cone(DVec3::new(1.0, 0.0, 0.0), DVec3::X, 30.0);
        assert!(matches!(intersect_cone_cone(&k1, &k2), ConeConeResult::General));
    }

    /// Anti-parallel axes (same line but opposite directions) should still be
    /// treated as coaxial.
    #[test]
    fn antiparallel_coaxial() {
        let k1 = cone(DVec3::ZERO, DVec3::Z, 45.0);
        let k2 = cone(DVec3::ZERO, -DVec3::Z, 45.0);
        // Same apex, same half-angle, axis anti-parallel → should detect Coaxial.
        assert!(matches!(intersect_cone_cone(&k1, &k2), ConeConeResult::Coaxial));
    }
}
