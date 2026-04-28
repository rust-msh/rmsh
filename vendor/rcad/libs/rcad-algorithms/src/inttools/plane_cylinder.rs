#[allow(unused_imports)]
use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

#[derive(Debug, Clone)]
pub enum PlaneCylinderResult {
    NoIntersection,
    TangentLine(Line3),
    TwoLines(Line3, Line3),
    Circle(Circle3),
    Ellipse(Ellipse3),
}

pub fn intersect_plane_cylinder(plane: &Plane, cyl: &CylindricalSurface) -> PlaneCylinderResult {
    let cos_angle = plane.normal.dot(cyl.axis).abs();

    if cos_angle < TOLERANCE_ANG {
        // Plane parallel to cylinder axis
        let axis_to_plane = (plane.origin - cyl.origin).dot(plane.normal);
        let dist = axis_to_plane.abs();

        if dist > cyl.radius + TOLERANCE_ABS {
            return PlaneCylinderResult::NoIntersection;
        }
        if (dist - cyl.radius).abs() < TOLERANCE_ABS {
            let tang_point = cyl.origin + plane.normal * (-axis_to_plane);
            return PlaneCylinderResult::TangentLine(Line3 {
                origin: tang_point,
                direction: cyl.axis,
            });
        }
        let offset_dir = plane.normal.cross(cyl.axis).normalize();
        let half_chord = (cyl.radius * cyl.radius - dist * dist).sqrt();
        let center_on_plane = cyl.origin - plane.normal * axis_to_plane;

        let l1_origin = center_on_plane + offset_dir * half_chord;
        let l2_origin = center_on_plane - offset_dir * half_chord;

        return PlaneCylinderResult::TwoLines(
            Line3 {
                origin: l1_origin,
                direction: cyl.axis,
            },
            Line3 {
                origin: l2_origin,
                direction: cyl.axis,
            },
        );
    }

    if (cos_angle - 1.0).abs() < TOLERANCE_ANG {
        // Plane perpendicular to cylinder axis → circle
        let t = (plane.origin - cyl.origin).dot(cyl.axis);
        let center = cyl.origin + cyl.axis * t;
        return PlaneCylinderResult::Circle(Circle3 {
            center,
            normal: cyl.axis,
            radius: cyl.radius,
        });
    }

    // General oblique case → ellipse
    let major_radius = cyl.radius / cos_angle;
    let minor_radius = cyl.radius;

    let t = (plane.origin - cyl.origin).dot(plane.normal) / cyl.axis.dot(plane.normal);
    let center = cyl.origin + cyl.axis * t;

    let major_dir = (cyl.axis - plane.normal * cyl.axis.dot(plane.normal)).normalize();

    PlaneCylinderResult::Ellipse(Ellipse3 {
        center,
        normal: plane.normal,
        major_dir,
        major_radius,
        minor_radius,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perpendicular_plane_gives_circle() {
        let plane = Plane {
            origin: DVec3::new(0.0, 5.0, 0.0),
            normal: DVec3::Y,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        };
        match intersect_plane_cylinder(&plane, &cyl) {
            PlaneCylinderResult::Circle(c) => {
                assert!((c.radius - 2.0).abs() < TOLERANCE_ABS);
                assert!((c.center.y - 5.0).abs() < TOLERANCE_ABS);
            }
            other => panic!("Expected Circle, got {other:?}"),
        }
    }

    #[test]
    fn parallel_plane_two_lines() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        };
        match intersect_plane_cylinder(&plane, &cyl) {
            PlaneCylinderResult::TwoLines(l1, l2) => {
                // Both lines should be along Y at z = ±2
                assert!(vectors_parallel(l1.direction, DVec3::Y));
                assert!(vectors_parallel(l2.direction, DVec3::Y));
                assert!((l1.origin.x).abs() < TOLERANCE_ABS);
                assert!((l2.origin.x).abs() < TOLERANCE_ABS);
            }
            other => panic!("Expected TwoLines, got {other:?}"),
        }
    }

    #[test]
    fn parallel_plane_no_intersection() {
        let plane = Plane {
            origin: DVec3::new(5.0, 0.0, 0.0),
            normal: DVec3::X,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        };
        assert!(matches!(
            intersect_plane_cylinder(&plane, &cyl),
            PlaneCylinderResult::NoIntersection
        ));
    }

    #[test]
    fn oblique_plane_gives_ellipse() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(0.0, 1.0, 1.0).normalize(),
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        };
        match intersect_plane_cylinder(&plane, &cyl) {
            PlaneCylinderResult::Ellipse(e) => {
                assert!((e.minor_radius - 1.0).abs() < TOLERANCE_ABS);
                assert!(e.major_radius > 1.0);
            }
            other => panic!("Expected Ellipse, got {other:?}"),
        }
    }

    #[test]
    fn tangent_plane() {
        let plane = Plane {
            origin: DVec3::new(2.0, 0.0, 0.0),
            normal: DVec3::X,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        };
        assert!(matches!(
            intersect_plane_cylinder(&plane, &cyl),
            PlaneCylinderResult::TangentLine(_)
        ));
    }
}
