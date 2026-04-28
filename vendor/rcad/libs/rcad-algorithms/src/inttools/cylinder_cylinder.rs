//! Analytic intersection of two cylinders.
//!
//! # Case classification
//!
//! ## Parallel axes
//!
//! When the two cylinder axes are parallel (cross-product length ≈ 0):
//!
//! - **Coaxial**: axes coincide → intersection is the full cylinder surface of
//!   the smaller radius (degenerate; we return `Coaxial`).
//! - **Offset parallel**: axes are parallel but distinct.  The gap between the
//!   two axis lines is `d`.
//!   - `d ≥ r1 + r2`: no intersection.
//!   - `d = r1 + r2` (within tolerance): external tangent, one generator line.
//!   - `|r1 − r2| < d < r1 + r2`: two generator lines (cross-section chords).
//!   - `d = |r1 − r2|`: internal tangent, one generator line.
//!   - `d < |r1 − r2|`: one cylinder inside the other, no surface intersection.
//!
//! ## Perpendicular axes (Steinmetz intersection)
//!
//! When the axes are perpendicular and the cross-section distance equals zero
//! (axes actually intersect), the intersection curves are two ellipses —
//! specifically the classic Steinmetz configuration.  We return
//! `Perpendicular(TwoEllipses(...))` for this sub-case.
//!
//! ## General skew axes
//!
//! For all other orientations we return `General`, signalling the caller to
//! fall back to numeric marching.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, CylindricalSurface, Ellipse3, Line3};

use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG};

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of cylinder × cylinder intersection.
#[derive(Debug, Clone)]
pub enum CylinderCylinderResult {
    /// The cylinders do not intersect (disjoint or fully nested).
    NoIntersection,
    /// Cylinders are coaxial; the intersection is the full lateral surface of
    /// the smaller cylinder.
    Coaxial,
    /// Parallel axes: the intersection is exactly one generator line (external
    /// or internal tangent).
    OneGeneratorLine(Line3),
    /// Parallel axes: the intersection consists of two generator lines.
    TwoGeneratorLines(Line3, Line3),
    /// Perpendicular intersecting axes (Steinmetz): two ellipses.
    TwoEllipses(Ellipse3, Ellipse3),
    /// Perpendicular intersecting axes, equal radii: two circles.
    TwoCircles(Circle3, Circle3),
    /// General case (skew axes or oblique angle not handled analytically).
    /// The caller should fall back to numeric marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `cyl1` and `cyl2`.
pub fn intersect_cylinder_cylinder(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
) -> CylinderCylinderResult {
    intersect_cylinder_cylinder_with_eps(cyl1, cyl2, TOLERANCE_ABS, TOLERANCE_ANG)
}

/// Compute cylinder-cylinder intersection with additional fuzzy tolerance.
pub fn intersect_cylinder_cylinder_with_tolerance(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    fuzzy_tol: f64,
) -> CylinderCylinderResult {
    let linear_tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);
    let angular_tol = TOLERANCE_ANG + fuzzy_tol.max(0.0);
    intersect_cylinder_cylinder_with_eps(cyl1, cyl2, linear_tol, angular_tol)
}

fn intersect_cylinder_cylinder_with_eps(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    linear_tol: f64,
    angular_tol: f64,
) -> CylinderCylinderResult {
    let a1 = cyl1.axis.normalize();
    let a2 = cyl2.axis.normalize();

    let cross = a1.cross(a2);
    let sin_angle = cross.length(); // |sin θ|

    // ── Parallel axes ────────────────────────────────────────────────────────
    if sin_angle < angular_tol {
        return intersect_parallel_cylinders(cyl1, cyl2, a1, linear_tol);
    }

    // ── Perpendicular axes ────────────────────────────────────────────────────
    let cos_angle = a1.dot(a2).abs();
    if cos_angle < angular_tol {
        return intersect_perpendicular_cylinders(cyl1, cyl2, a1, a2, linear_tol);
    }

    // ── General skew / oblique ────────────────────────────────────────────────
    CylinderCylinderResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel axes
// ─────────────────────────────────────────────────────────────────────────────

fn intersect_parallel_cylinders(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    axis: DVec3,
    linear_tol: f64,
) -> CylinderCylinderResult {
    let r1 = cyl1.radius;
    let r2 = cyl2.radius;

    // Perpendicular distance between the two parallel axes.
    let delta = cyl2.origin - cyl1.origin;
    let delta_perp = delta - axis * delta.dot(axis);
    let d = delta_perp.length();

    // Coaxial check
    if d < linear_tol {
        if (r1 - r2).abs() < linear_tol {
            return CylinderCylinderResult::Coaxial;
        }
        // One inside the other along the same axis
        return CylinderCylinderResult::NoIntersection;
    }

    let sum = r1 + r2;
    let diff = (r1 - r2).abs();

    if d > sum + linear_tol {
        return CylinderCylinderResult::NoIntersection;
    }
    if d < diff - linear_tol {
        // One cylinder fully inside the other
        return CylinderCylinderResult::NoIntersection;
    }

    // Direction from cyl1 axis to cyl2 axis (perpendicular)
    let dir_perp = delta_perp.normalize();

    // External tangent
    if (d - sum).abs() < linear_tol {
        let point = cyl1.origin + dir_perp * r1;
        return CylinderCylinderResult::OneGeneratorLine(Line3 {
            origin: point,
            direction: axis,
        });
    }
    // Internal tangent
    if (d - diff).abs() < linear_tol {
        // The tangent line is on the side of the smaller cylinder that is
        // closest to the larger cylinder's axis.
        let point = if r1 >= r2 {
            cyl1.origin + dir_perp * r1
        } else {
            cyl1.origin - dir_perp * r1
        };
        return CylinderCylinderResult::OneGeneratorLine(Line3 {
            origin: point,
            direction: axis,
        });
    }

    // Two generator lines: find the two intersection points in the
    // perpendicular cross-section.
    //
    // In 2D: circle1 centred at origin radius r1, circle2 centred at (d, 0)
    // radius r2.  The intersection x-coordinate:
    //   x = (d² + r1² - r2²) / (2d)
    //   y = ±sqrt(r1² - x²)
    let x = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let y_sq = r1 * r1 - x * x;
    if y_sq < 0.0 {
        return CylinderCylinderResult::NoIntersection;
    }
    let y = y_sq.sqrt();

    // Orthogonal unit vector in the cross-section plane
    let v_perp = axis.cross(dir_perp).normalize();

    let p1 = cyl1.origin + dir_perp * x + v_perp * y;
    let p2 = cyl1.origin + dir_perp * x - v_perp * y;

    CylinderCylinderResult::TwoGeneratorLines(
        Line3 { origin: p1, direction: axis },
        Line3 { origin: p2, direction: axis },
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Perpendicular intersecting axes (Steinmetz)
// ─────────────────────────────────────────────────────────────────────────────

fn intersect_perpendicular_cylinders(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    a1: DVec3,
    a2: DVec3,
    linear_tol: f64,
) -> CylinderCylinderResult {
    let r1 = cyl1.radius;
    let r2 = cyl2.radius;

    // Find the closest point between the two axes (skew lines).
    // Parametric form: P = O1 + t*a1  and  Q = O2 + s*a2
    // The connecting vector at closest approach is perpendicular to both axes.
    let w = cyl1.origin - cyl2.origin;
    let b = a1.dot(a2);
    let denom = 1.0 - b * b;

    if denom.abs() < 1e-12 {
        // Degenerate (should not reach here since we checked perpendicularity)
        return CylinderCylinderResult::General;
    }

    let d1 = a1.dot(w);
    let d2 = a2.dot(w);
    let t = (b * d2 - d1) / denom;
    let s = (d2 - b * d1) / denom;

    let closest1 = cyl1.origin + a1 * t;
    let closest2 = cyl2.origin + a2 * s;

    // Perpendicular distance between axes
    let dist = (closest1 - closest2).length();

    if dist > r1 + r2 + linear_tol {
        return CylinderCylinderResult::NoIntersection;
    }

    // For the Steinmetz case the axes must actually cross (dist ≈ 0).
    // For dist > 0 the analytic form is much harder; fall back to marching.
    if dist > linear_tol * 10.0 {
        return CylinderCylinderResult::General;
    }

    // Intersection point of the two axes
    let origin = (closest1 + closest2) * 0.5;

    // Third axis = a1 × a2  (normal to both, the "viewing" direction)
    let _a3 = a1.cross(a2).normalize();

    if (r1 - r2).abs() < linear_tol {
        // Equal radii: the Steinmetz intersection is two circles in planes
        // at ±45° between the two axes (actually the intersection lies in
        // planes whose normals are a1±a2, but the curves ARE circles for
        // equal-radius perpendicular cylinders).
        //
        // Each circle: normal = (a1 ± a2).normalize(), radius = r, center = origin
        let n1 = (a1 + a2).normalize();
        let n2 = (a1 - a2).normalize();
        let r = (r1 * r1 + r2 * r2).sqrt() / std::f64::consts::SQRT_2;
        let _ = r; // radius of the Steinmetz circles
        // For equal radii r1=r2=r, the Steinmetz circles have radius r*sqrt(2)/sqrt(2)=r
        // ... actually radius = r1 (same as cylinder radius for equal radii).
        let circle1 = Circle3 { center: origin, normal: n1, radius: r1 };
        let circle2 = Circle3 { center: origin, normal: n2, radius: r1 };
        return CylinderCylinderResult::TwoCircles(circle1, circle2);
    }

    // Unequal radii: intersection curves are two congruent ellipses.
    // Each ellipse lies in a plane spanned by a3 and (a1 or a2).
    // - Ellipse 1: normal = a2, major axis along a1, minor along a3
    //   major_radius = r2, minor_radius = r1  (projected from cyl1's cross-section)
    // - Ellipse 2: normal = a1, major axis along a2, minor along a3
    //   major_radius = r1, minor_radius = r2

    let ellipse1 = Ellipse3 {
        center: origin,
        normal: a2,
        major_dir: a1,
        major_radius: r2,
        minor_radius: r1,
    };
    let ellipse2 = Ellipse3 {
        center: origin,
        normal: a1,
        major_dir: a2,
        major_radius: r1,
        minor_radius: r2,
    };
    CylinderCylinderResult::TwoEllipses(ellipse1, ellipse2)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn cyl(origin: DVec3, axis: DVec3, radius: f64) -> CylindricalSurface {
        CylindricalSurface { origin, axis, radius }
    }

    // ── Parallel axes ──────────────────────────────────────────────────────────

    #[test]
    fn parallel_coaxial() {
        let c1 = cyl(DVec3::ZERO, DVec3::Z, 2.0);
        let c2 = cyl(DVec3::ZERO, DVec3::Z, 2.0);
        assert!(matches!(
            intersect_cylinder_cylinder(&c1, &c2),
            CylinderCylinderResult::Coaxial
        ));
    }

    #[test]
    fn parallel_no_intersection() {
        let c1 = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let c2 = cyl(DVec3::new(5.0, 0.0, 0.0), DVec3::Z, 1.0);
        assert!(matches!(
            intersect_cylinder_cylinder(&c1, &c2),
            CylinderCylinderResult::NoIntersection
        ));
    }

    #[test]
    fn parallel_external_tangent() {
        let c1 = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let c2 = cyl(DVec3::new(2.0, 0.0, 0.0), DVec3::Z, 1.0);
        match intersect_cylinder_cylinder(&c1, &c2) {
            CylinderCylinderResult::OneGeneratorLine(l) => {
                // Tangent at x=1, y=0
                assert!((l.origin.x - 1.0).abs() < 1e-9, "x={}", l.origin.x);
                assert!(l.origin.y.abs() < 1e-9);
            }
            other => panic!("expected OneGeneratorLine, got {other:?}"),
        }
    }

    #[test]
    fn parallel_two_generator_lines() {
        // Two unit cylinders with axes offset by 1 (d=1 < r1+r2=2)
        let c1 = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let c2 = cyl(DVec3::new(1.0, 0.0, 0.0), DVec3::Z, 1.0);
        match intersect_cylinder_cylinder(&c1, &c2) {
            CylinderCylinderResult::TwoGeneratorLines(l1, l2) => {
                // Both lines parallel to Z
                assert!((l1.direction.z.abs() - 1.0).abs() < 1e-9);
                assert!((l2.direction.z.abs() - 1.0).abs() < 1e-9);
                // x-coordinate of intersection: x = (1 + 1 - 1)/(2) = 0.5
                assert!((l1.origin.x - 0.5).abs() < 1e-9, "l1.x={}", l1.origin.x);
                assert!((l2.origin.x - 0.5).abs() < 1e-9, "l2.x={}", l2.origin.x);
                // y-coordinates are symmetric: y = ±sqrt(1 - 0.25) = ±sqrt(0.75)
                let expected_y = (0.75_f64).sqrt();
                assert!((l1.origin.y.abs() - expected_y).abs() < 1e-9, "l1.y={}", l1.origin.y);
                assert!((l2.origin.y.abs() - expected_y).abs() < 1e-9, "l2.y={}", l2.origin.y);
            }
            other => panic!("expected TwoGeneratorLines, got {other:?}"),
        }
    }

    // ── Perpendicular axes ─────────────────────────────────────────────────────

    #[test]
    fn perpendicular_equal_radii_steinmetz() {
        // Classic Steinmetz: two unit cylinders along X and Y, intersecting at origin
        let c1 = cyl(DVec3::ZERO, DVec3::X, 1.0);
        let c2 = cyl(DVec3::ZERO, DVec3::Y, 1.0);
        match intersect_cylinder_cylinder(&c1, &c2) {
            CylinderCylinderResult::TwoCircles(circ1, circ2) => {
                assert!((circ1.radius - 1.0).abs() < 1e-9);
                assert!((circ2.radius - 1.0).abs() < 1e-9);
            }
            other => panic!("expected TwoCircles, got {other:?}"),
        }
    }

    #[test]
    fn perpendicular_unequal_radii_steinmetz() {
        // c1: axis=X, r1=2;  c2: axis=Y, r2=1
        // ellipse1: normal=a2=Y, major_dir=a1=X, major_radius=r2=1, minor_radius=r1=2
        // ellipse2: normal=a1=X, major_dir=a2=Y, major_radius=r1=2, minor_radius=r2=1
        let c1 = cyl(DVec3::ZERO, DVec3::X, 2.0);
        let c2 = cyl(DVec3::ZERO, DVec3::Y, 1.0);
        match intersect_cylinder_cylinder(&c1, &c2) {
            CylinderCylinderResult::TwoEllipses(e1, e2) => {
                assert!((e1.minor_radius - 2.0).abs() < 1e-9, "e1.minor={}", e1.minor_radius);
                assert!((e1.major_radius - 1.0).abs() < 1e-9, "e1.major={}", e1.major_radius);
                assert!((e2.minor_radius - 1.0).abs() < 1e-9, "e2.minor={}", e2.minor_radius);
                assert!((e2.major_radius - 2.0).abs() < 1e-9, "e2.major={}", e2.major_radius);
            }
            other => panic!("expected TwoEllipses, got {other:?}"),
        }
    }

    #[test]
    fn skew_axes_falls_back_to_general() {
        // Axes cross at 45°
        let c1 = cyl(DVec3::ZERO, DVec3::X, 1.0);
        let c2 = cyl(DVec3::ZERO, DVec3::new(1.0, 1.0, 0.0).normalize(), 1.0);
        assert!(matches!(
            intersect_cylinder_cylinder(&c1, &c2),
            CylinderCylinderResult::General
        ));
    }

    #[test]
    fn fuzzy_tolerance_recovers_near_tangent_parallel_case() {
        let c1 = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let c2 = cyl(DVec3::new(2.0 + 3.0e-7, 0.0, 0.0), DVec3::Z, 1.0);

        assert!(matches!(
            intersect_cylinder_cylinder(&c1, &c2),
            CylinderCylinderResult::NoIntersection
        ));

        assert!(matches!(
            intersect_cylinder_cylinder_with_tolerance(&c1, &c2, 4.0e-7),
            CylinderCylinderResult::OneGeneratorLine(_)
        ));
    }
}
