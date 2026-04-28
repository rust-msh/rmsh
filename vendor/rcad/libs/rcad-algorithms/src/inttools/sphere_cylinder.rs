//! Analytic intersection of a sphere and a cylinder.
//!
//! # Case classification
//!
//! ## Axis-aligned case (sphere centre on cylinder axis)
//!
//! When the sphere centre **C** lies on the cylinder axis (`d_perp ≈ 0`), the
//! intersection degenerates to one or two circles perpendicular to the axis:
//!
//! ```text
//! h = h_c ± √(R² − r²)
//! ```
//!
//! where
//! - `h_c = (C − O) · â`  (height of sphere centre above cylinder origin),
//! - `R` = sphere radius,
//! - `r` = cylinder radius,
//! - `â` = unit cylinder axis,
//! - `h` = height of the intersection circle on the cylinder axis.
//!
//! Each such `h` yields a circle of radius `r` centred at `O + h · â` with
//! normal `â`.
//!
//! ## Parallel-axis offset case
//!
//! When the sphere centre is off-axis but the sphere and cylinder axes are
//! parallel (or the cylinder has no preferred axis direction relative to the
//! sphere), we can still decide:
//!
//! - Let `d` = perpendicular distance from sphere centre to cylinder axis.
//! - The sphere surface is at radial distances `[d − R, d + R]` from the axis.
//! - The cylinder surface is at radial distance `r` from the axis.
//!
//! Therefore:
//! - If `d − R > r` or `d + R < r` (and `d > r` for the latter): **no intersection**.
//! - If `|d − r| ≤ R`: the sphere surface intersects the cylinder surface;
//!   the exact intersection is a quartic (Viviani-type) curve — return `General`.
//!
//! ## General case
//!
//! For all other configurations (arbitrary relative orientation of sphere centre
//! and cylinder axis) the intersection is a quartic space curve.  We return
//! `General` so the caller can fall back to numeric marching.

use rcad_kernel::geom::{Circle3, CylindricalSurface, SphericalSurface};

use crate::tolerance::TOLERANCE_ABS;

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of sphere × cylinder intersection.
#[derive(Debug, Clone)]
pub enum SphereCylinderResult {
    /// Sphere and cylinder do not intersect (disjoint or one fully inside the
    /// other when the axis-aligned condition also holds).
    NoIntersection,
    /// Exactly one tangent circle (R = r and sphere centre on axis, or the two
    /// roots coincide).
    TangentCircle(Circle3),
    /// Two distinct intersection circles (axis-aligned, `R > r`).
    TwoCircles(Circle3, Circle3),
    /// The intersection is a quartic space curve.  The caller should fall back
    /// to numeric marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `sphere` and `cyl`.
///
/// Returns one of [`SphereCylinderResult`]'s variants:
///
/// - [`NoIntersection`](SphereCylinderResult::NoIntersection) — disjoint.
/// - [`TangentCircle`](SphereCylinderResult::TangentCircle) — one tangent circle
///   (axis-aligned case, discriminant = 0).
/// - [`TwoCircles`](SphereCylinderResult::TwoCircles) — two circles (axis-aligned).
/// - [`General`](SphereCylinderResult::General) — quartic; fall back to marching.
///
/// The axis-aligned tolerance is ten times the absolute position tolerance.
pub fn intersect_sphere_cylinder(
    sphere: &SphericalSurface,
    cyl: &CylindricalSurface,
) -> SphereCylinderResult {
    intersect_sphere_cylinder_with_tolerance(sphere, cyl, 0.0)
}

/// Compute sphere-cylinder intersection with additional fuzzy tolerance.
///
/// This relaxes axis-aligned and distance early-out checks by `fuzzy_tol` so
/// near-coincident cases can still classify into analytic branches.
pub fn intersect_sphere_cylinder_with_tolerance(
    sphere: &SphericalSurface,
    cyl: &CylindricalSurface,
    fuzzy_tol: f64,
) -> SphereCylinderResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);
    let axis = cyl.axis.normalize();
    let d = sphere.center - cyl.origin;
    let d_along = d.dot(axis);
    let d_perp_vec = d - axis * d_along;
    let d_perp = d_perp_vec.length(); // perpendicular distance: sphere centre → cyl axis

    let r = cyl.radius;
    let big_r = sphere.radius;

    // ── Axis-aligned case ─────────────────────────────────────────────────────
    if d_perp < tol * 10.0 {
        // Sphere centre is on (or extremely close to) the cylinder axis.
        let disc = big_r * big_r - r * r;

        if disc < -tol {
            return SphereCylinderResult::NoIntersection;
        }

        let h_c = d_along;

        if disc.abs() < tol {
            let center = cyl.origin + axis * h_c;
            return SphereCylinderResult::TangentCircle(Circle3 {
                center,
                normal: axis,
                radius: r,
            });
        }

        let delta_h = disc.sqrt();
        let c1 = Circle3 {
            center: cyl.origin + axis * (h_c - delta_h),
            normal: axis,
            radius: r,
        };
        let c2 = Circle3 {
            center: cyl.origin + axis * (h_c + delta_h),
            normal: axis,
            radius: r,
        };
        return SphereCylinderResult::TwoCircles(c1, c2);
    }

    // ── Off-axis: early-out distance test ─────────────────────────────────────
    //
    // The cylinder lateral surface is everywhere at distance `r` from the axis.
    // The sphere surface spans radial distances [d_perp − R, d_perp + R] from
    // the axis (considering all points on the sphere surface).
    //
    // No intersection when:
    //   (a) d_perp - R > r  →  sphere is entirely outside the cylinder (far side)
    //   (b) d_perp + R < r  →  sphere is entirely inside the cylinder (near side)
    //       but only when the sphere is smaller than the cylinder radius + offset
    //
    // Case (a): closest radial approach of sphere to axis exceeds cylinder radius.
    if d_perp - big_r > r + tol {
        return SphereCylinderResult::NoIntersection;
    }
    // Case (b): sphere is fully enclosed inside the cylinder laterally.
    if d_perp + big_r < r - tol {
        return SphereCylinderResult::NoIntersection;
    }

    // ── Quartic (Viviani-type) intersection ───────────────────────────────────
    SphereCylinderResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    fn sphere(center: DVec3, radius: f64) -> SphericalSurface {
        SphericalSurface { center, axis: DVec3::Z, radius }
    }

    fn cyl(origin: DVec3, radius: f64) -> CylindricalSurface {
        CylindricalSurface { origin, axis: DVec3::Z, radius }
    }

    // ── Axis-aligned ──────────────────────────────────────────────────────────

    /// R > r, centre on axis → two circles
    #[test]
    fn two_circles_axis_aligned() {
        let sph = sphere(DVec3::new(0.0, 0.0, 3.0), 5.0);
        let c = cyl(DVec3::ZERO, 4.0);
        match intersect_sphere_cylinder(&sph, &c) {
            SphereCylinderResult::TwoCircles(c1, c2) => {
                // delta_h = sqrt(25 - 16) = 3
                // h1 = 3 - 3 = 0, h2 = 3 + 3 = 6
                assert!((c1.center.z - 0.0).abs() < 1e-9, "c1.z={}", c1.center.z);
                assert!((c2.center.z - 6.0).abs() < 1e-9, "c2.z={}", c2.center.z);
                assert!((c1.radius - 4.0).abs() < 1e-9);
                assert!((c2.radius - 4.0).abs() < 1e-9);
            }
            other => panic!("expected TwoCircles, got {other:?}"),
        }
    }

    /// R = r, centre on axis → tangent circle
    #[test]
    fn tangent_circle_equal_radii() {
        let sph = sphere(DVec3::new(0.0, 0.0, 5.0), 3.0);
        let c = cyl(DVec3::ZERO, 3.0);
        match intersect_sphere_cylinder(&sph, &c) {
            SphereCylinderResult::TangentCircle(tc) => {
                assert!((tc.center.z - 5.0).abs() < 1e-9, "tc.z={}", tc.center.z);
                assert!((tc.radius - 3.0).abs() < 1e-9);
            }
            other => panic!("expected TangentCircle, got {other:?}"),
        }
    }

    /// R < r, centre on axis → no intersection
    #[test]
    fn no_intersection_sphere_smaller_axis_aligned() {
        let sph = sphere(DVec3::ZERO, 1.0);
        let c = cyl(DVec3::ZERO, 2.0);
        assert!(matches!(
            intersect_sphere_cylinder(&sph, &c),
            SphereCylinderResult::NoIntersection
        ));
    }

    // ── Off-axis distance tests ───────────────────────────────────────────────

    /// Sphere centre far off-axis, sphere entirely outside cylinder.
    /// d_perp=10, R=2, r=1 → d_perp - R = 8 > r = 1 → NoIntersection.
    #[test]
    fn no_intersection_off_axis_too_far() {
        let sph = sphere(DVec3::new(10.0, 0.0, 0.0), 2.0);
        let c = cyl(DVec3::ZERO, 1.0);
        assert!(matches!(
            intersect_sphere_cylinder(&sph, &c),
            SphereCylinderResult::NoIntersection
        ));
    }

    /// Sphere centre off-axis but sphere is small and entirely inside cylinder.
    /// d_perp=0.5, R=0.1, r=2 → d_perp + R = 0.6 < r = 2 → NoIntersection.
    #[test]
    fn no_intersection_off_axis_inside_cylinder() {
        let sph = sphere(DVec3::new(0.5, 0.0, 0.0), 0.1);
        let c = cyl(DVec3::ZERO, 2.0);
        assert!(matches!(
            intersect_sphere_cylinder(&sph, &c),
            SphereCylinderResult::NoIntersection
        ));
    }

    /// Sphere centre off-axis, sphere large enough to reach the cylinder surface.
    /// d_perp=1, R=5, r=2 → d_perp - R = -4 < r; d_perp + R = 6 > r → General.
    #[test]
    fn general_off_axis_intersecting() {
        let sph = sphere(DVec3::new(1.0, 0.0, 0.0), 5.0);
        let c = cyl(DVec3::ZERO, 2.0);
        assert!(matches!(intersect_sphere_cylinder(&sph, &c), SphereCylinderResult::General));
    }

    /// Classic test: sphere centre at (1,0,0), far from cylinder → General (was already).
    #[test]
    fn general_off_axis_original() {
        let sph = sphere(DVec3::new(1.0, 0.0, 0.0), 5.0);
        let c = cyl(DVec3::ZERO, 2.0);
        assert!(matches!(intersect_sphere_cylinder(&sph, &c), SphereCylinderResult::General));
    }
}
