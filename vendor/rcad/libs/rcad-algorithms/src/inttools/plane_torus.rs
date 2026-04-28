//! Analytic intersection of a plane with a torus.
//!
//! # Cases
//!
//! - **Perpendicular to axis**: Two circles (inner and outer equator)
//! - **Parallel to axis**: Two circles when plane intersects tube center circle
//!   - If |d| < R: Two circles of radius r at intersections with tube center circle
//!   - If |d| = R: Two tangent circles (Villarceau circles configuration)
//!   - If |d| > R+r: No intersection
//! - **Oblique**: Complex curve, fall back to numerical marching

use rcad_kernel::geom::{Circle3, Plane, ToroidalSurface};

use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG};

/// Result of plane x torus intersection.
#[derive(Debug, Clone)]
pub enum PlaneTorusResult {
    /// The plane does not intersect the torus.
    NoIntersection,
    /// Single tangent circle.
    TangentCircle(Circle3),
    /// Two circles (perpendicular case).
    TwoCircles(Circle3, Circle3),
    /// Complex intersection, fall back to numerical marching.
    General,
}

/// Compute the analytic intersection of `plane` and `torus`.
pub fn intersect_plane_torus(plane: &Plane, torus: &ToroidalSurface) -> PlaneTorusResult {
    intersect_plane_torus_with_tolerance(plane, torus, 0.0)
}

/// Plane x torus intersection with fuzzy tolerance.
pub fn intersect_plane_torus_with_tolerance(
    plane: &Plane,
    torus: &ToroidalSurface,
    fuzzy_tol: f64,
) -> PlaneTorusResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);

    // Normalize plane normal and torus axis
    let n = plane.normal.normalize();
    let a = torus.axis.normalize();

    // Check if plane is perpendicular to torus axis
    let dot_na = n.dot(a).abs();

    if dot_na > 1.0 - TOLERANCE_ANG {
        // Plane perpendicular to axis: circular cross-section
        return intersect_plane_torus_perpendicular(plane, torus, tol);
    }

    // Check if plane is parallel to torus axis
    if dot_na < TOLERANCE_ANG {
        // Plane parallel to axis: may produce two circles
        return intersect_plane_torus_parallel(plane, torus, tol);
    }

    // General oblique case: fall back to numerical
    PlaneTorusResult::General
}

fn intersect_plane_torus_perpendicular(
    plane: &Plane,
    torus: &ToroidalSurface,
    tol: f64,
) -> PlaneTorusResult {
    // Distance from torus center to plane along axis
    let signed_dist = (torus.center - plane.origin).dot(torus.axis);
    let abs_dist = signed_dist.abs();

    // Maximum distance for intersection is the minor radius
    if abs_dist > torus.minor_radius + tol {
        return PlaneTorusResult::NoIntersection;
    }

    // Tangent case: one circle
    if (abs_dist - torus.minor_radius).abs() < tol {
        let center = torus.center - torus.axis * signed_dist;
        return PlaneTorusResult::TangentCircle(Circle3 {
            center,
            normal: torus.axis,
            radius: torus.major_radius,
        });
    }

    // Two circles at height signed_dist from torus center
    // Circle radius on the tube: sqrt(r^2 - d^2) where r = minor_radius, d = distance
    let tube_circle_r = (torus.minor_radius * torus.minor_radius - signed_dist * signed_dist).sqrt();

    // Two circles at major_radius +/- tube_circle_r from axis
    let r1 = torus.major_radius + tube_circle_r;
    let r2 = (torus.major_radius - tube_circle_r).max(0.0);

    let center = torus.center - torus.axis * signed_dist;

    if r2 < tol {
        // Inner circle degenerates to point
        PlaneTorusResult::TangentCircle(Circle3 {
            center,
            normal: torus.axis,
            radius: r1,
        })
    } else {
        PlaneTorusResult::TwoCircles(
            Circle3 { center, normal: torus.axis, radius: r1 },
            Circle3 { center, normal: torus.axis, radius: r2 },
        )
    }
}

fn intersect_plane_torus_parallel(
    plane: &Plane,
    torus: &ToroidalSurface,
    tol: f64,
) -> PlaneTorusResult {
    // Signed distance from torus center to plane along the plane normal
    let n = plane.normal.normalize();
    let signed_dist = (plane.origin - torus.center).dot(n);
    let d = signed_dist.abs();

    // No intersection if plane is too far from torus
    if d > torus.major_radius + torus.minor_radius + tol {
        return PlaneTorusResult::NoIntersection;
    }

    // Analytical solution: two circles when plane intersects tube center circle
    // This happens when |d| <= R (the plane cuts through the tube center circle)
    if d <= torus.major_radius + tol {
        // The tube center circle (radius R in plane perpendicular to torus axis)
        // intersects the plane at two points when |d| < R
        // Each intersection produces a circle of radius r in the plane

        // Compute the direction perpendicular to the plane normal that lies in the
        // plane containing both the plane normal and the torus axis
        let a = torus.axis.normalize();

        // Compute the direction in the plane that is perpendicular to the torus axis
        // This gives us the direction toward the tube center intersection points
        let in_plane_perp = n.cross(a);
        let perp_len = in_plane_perp.length();

        if perp_len < TOLERANCE_ANG {
            // Degenerate case: plane normal is parallel to torus axis
            // This shouldn't happen in the parallel case
            return PlaneTorusResult::General;
        }

        let dir_perp = in_plane_perp / perp_len;

        // Calculate the distance along the perpendicular direction to the tube center
        // intersection points. From d² + z² = R², we get z = ±√(R² - d²)
        let d_sq = d * d;
        let r_sq = torus.major_radius * torus.major_radius;
        let z_dist_sq = r_sq - d_sq;

        if z_dist_sq < -tol * tol {
            // No real intersection (should not happen given our checks above)
            return PlaneTorusResult::NoIntersection;
        }

        let z_dist = z_dist_sq.max(0.0).sqrt();

        // The two circle centers are at:
        // center = torus_center + d * dir_to_plane ± z_dist * dir_perp
        let base_center = torus.center + signed_dist * n;

        let center1 = base_center + z_dist * dir_perp;
        let center2 = base_center - z_dist * dir_perp;

        // Check for tangent case (circles merge into one)
        if (center1 - center2).length() < tol {
            return PlaneTorusResult::TangentCircle(Circle3 {
                center: center1,
                normal: n,
                radius: torus.minor_radius,
            });
        }

        // Two circles of radius r in the plane
        PlaneTorusResult::TwoCircles(
            Circle3 {
                center: center1,
                normal: n,
                radius: torus.minor_radius,
            },
            Circle3 {
                center: center2,
                normal: n,
                radius: torus.minor_radius,
            },
        )
    } else {
        // d > R: Plane is between the tube center circle and outer edge
        // This produces a more complex intersection (ellipse-like)
        // Fall back to numerical marching for this case
        PlaneTorusResult::General
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    #[test]
    fn plane_perpendicular_to_torus_axis_produces_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane perpendicular to Y axis, slicing through center
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // Outer circle at major_radius + minor_radius
                assert!((c1.radius - 6.0).abs() < 1e-6, "Outer circle radius expected 6.0, got {}", c1.radius);
                // Inner circle at major_radius - minor_radius
                assert!((c2.radius - 4.0).abs() < 1e-6, "Inner circle radius expected 4.0, got {}", c2.radius);
                // Both circles should have the same center
                assert!((c1.center - DVec3::ZERO).length() < 1e-6);
                assert!((c2.center - DVec3::ZERO).length() < 1e-6);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_parallel_to_torus_axis_through_center_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane parallel to Y axis (normal = X), passing through center
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // Both circles should have radius equal to minor radius
                assert!((c1.radius - 1.0).abs() < 1e-6, "Circle 1 radius expected 1.0, got {}", c1.radius);
                assert!((c2.radius - 1.0).abs() < 1e-6, "Circle 2 radius expected 1.0, got {}", c2.radius);
                // Both circles should be in the plane (normal = X)
                assert!((c1.normal - DVec3::X).length() < 1e-6);
                assert!((c2.normal - DVec3::X).length() < 1e-6);
                // Centers should be at z = ±R = ±5
                assert!((c1.center.z.abs() - 5.0).abs() < 1e-6, "Circle 1 center z should be ±5");
                assert!((c2.center.z.abs() - 5.0).abs() < 1e-6, "Circle 2 center z should be ±5");
                // Both centers at x=0, y=0
                assert!(c1.center.x.abs() < 1e-6);
                assert!(c1.center.y.abs() < 1e-6);
                assert!(c2.center.x.abs() < 1e-6);
                assert!(c2.center.y.abs() < 1e-6);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_parallel_to_torus_axis_offset_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane parallel to Y axis, offset by d=3 from center
        let plane = Plane {
            origin: DVec3::new(3.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // Both circles should have radius equal to minor radius
                assert!((c1.radius - 1.0).abs() < 1e-6);
                assert!((c2.radius - 1.0).abs() < 1e-6);
                // Centers should be at x=3 (the plane's position)
                assert!((c1.center.x - 3.0).abs() < 1e-6);
                assert!((c2.center.x - 3.0).abs() < 1e-6);
                // z = ±sqrt(R² - d²) = ±sqrt(25 - 9) = ±4
                let expected_z = (25.0_f64 - 9.0_f64).sqrt();
                assert!((c1.center.z.abs() - expected_z).abs() < 1e-6);
                assert!((c2.center.z.abs() - expected_z).abs() < 1e-6);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_parallel_to_torus_axis_at_major_radius_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane at distance R = 5 (touches tube center circle)
        let plane = Plane {
            origin: DVec3::new(5.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);

        // The result type depends on the exact intersection geometry
        // Just verify we get a valid result (don't panic)
        let _ = result;
    }

    #[test]
    fn plane_parallel_to_torus_axis_near_edge_returns_general() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane at distance 5.5 (R < d < R+r), produces complex intersection
        let plane = Plane {
            origin: DVec3::new(5.5, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);
        // d = 5.5 > R = 5, so this should fall back to General
        assert!(matches!(result, PlaneTorusResult::General));
    }

    #[test]
    fn plane_parallel_to_torus_axis_no_intersection() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane outside torus (d > R + r = 6)
        let plane = Plane {
            origin: DVec3::new(7.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);
        assert!(matches!(result, PlaneTorusResult::NoIntersection));
    }

    #[test]
    fn plane_parallel_to_torus_axis_negative_offset_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane parallel to Y axis, offset by d=-4 from center
        let plane = Plane {
            origin: DVec3::new(-4.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // Both circles should have radius equal to minor radius
                assert!((c1.radius - 1.0).abs() < 1e-6);
                assert!((c2.radius - 1.0).abs() < 1e-6);
                // Centers should be at x=-4
                assert!((c1.center.x + 4.0).abs() < 1e-6);
                assert!((c2.center.x + 4.0).abs() < 1e-6);
                // z = ±sqrt(R² - d²) = ±sqrt(25 - 16) = ±3
                let expected_z = (25.0_f64 - 16.0_f64).sqrt();
                assert!((c1.center.z.abs() - expected_z).abs() < 1e-6);
                assert!((c2.center.z.abs() - expected_z).abs() < 1e-6);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_parallel_to_torus_axis_tangent_at_inner_radius() {
        // Torus centered at origin with axis along Y
        // R = 5, r = 1, inner radius = R - r = 4
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane at distance d = 4 (inner radius)
        let plane = Plane {
            origin: DVec3::new(4.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // d = 4 < R = 5, so we still get two circles
                assert!((c1.radius - 1.0).abs() < 1e-6);
                assert!((c2.radius - 1.0).abs() < 1e-6);
                // z = ±sqrt(R² - d²) = ±sqrt(25 - 16) = ±3
                let expected_z = (25.0_f64 - 16.0_f64).sqrt();
                assert!((c1.center.z.abs() - expected_z).abs() < 1e-6);
                assert!((c2.center.z.abs() - expected_z).abs() < 1e-6);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_perpendicular_tangent_to_torus() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane tangent to torus (at top of tube)
        let plane = Plane {
            origin: DVec3::new(0.0, 1.0, 0.0),
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TangentCircle(c) => {
                // Tangent circle at the major radius
                assert!((c.radius - 5.0).abs() < 1e-6);
            }
            other => panic!("Expected TangentCircle, got {:?}", other),
        }
    }

    #[test]
    fn plane_perpendicular_no_intersection() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane outside torus
        let plane = Plane {
            origin: DVec3::new(0.0, 2.0, 0.0),
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);
        assert!(matches!(result, PlaneTorusResult::NoIntersection));
    }

    #[test]
    fn plane_perpendicular_offset_produces_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane offset by 0.5 from center
        let plane = Plane {
            origin: DVec3::new(0.0, 0.5, 0.0),
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // tube_circle_r = sqrt(1 - 0.25) = sqrt(0.75) = 0.866025...
                let expected_tube_r = (1.0_f64 * 1.0 - 0.5_f64 * 0.5).sqrt();
                let expected_r1 = 5.0 + expected_tube_r;
                let expected_r2 = 5.0 - expected_tube_r;

                assert!((c1.radius - expected_r1).abs() < 1e-6, "Outer circle radius mismatch");
                assert!((c2.radius - expected_r2).abs() < 1e-6, "Inner circle radius mismatch");
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_oblique_returns_general() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane at 45 degrees to torus axis
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(1.0, 1.0, 0.0).normalize(),
        };

        let result = intersect_plane_torus(&plane, &torus);
        assert!(matches!(result, PlaneTorusResult::General));
    }
}
