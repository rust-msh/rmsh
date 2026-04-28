use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

#[derive(Debug, Clone)]
pub enum PlaneSphereResult {
    NoIntersection,
    TangentPoint(DVec3),
    Circle(Circle3),
}

pub fn intersect_plane_sphere(plane: &Plane, sphere: &SphericalSurface) -> PlaneSphereResult {
    let signed_dist = (sphere.center - plane.origin).dot(plane.normal);
    let abs_dist = signed_dist.abs();

    if abs_dist > sphere.radius + TOLERANCE_ABS {
        return PlaneSphereResult::NoIntersection;
    }
    if (abs_dist - sphere.radius).abs() < TOLERANCE_ABS {
        let point = sphere.center - plane.normal * signed_dist;
        return PlaneSphereResult::TangentPoint(point);
    }

    let circle_radius = (sphere.radius * sphere.radius - signed_dist * signed_dist).sqrt();
    let center = sphere.center - plane.normal * signed_dist;

    PlaneSphereResult::Circle(Circle3 {
        center,
        normal: plane.normal,
        radius: circle_radius,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plane_through_center() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        };
        match intersect_plane_sphere(&plane, &sphere) {
            PlaneSphereResult::Circle(c) => {
                assert!((c.radius - 3.0).abs() < TOLERANCE_ABS);
                assert!(points_coincide(c.center, DVec3::ZERO));
            }
            other => panic!("Expected Circle, got {other:?}"),
        }
    }

    #[test]
    fn plane_offset() {
        let plane = Plane {
            origin: DVec3::new(0.0, 2.0, 0.0),
            normal: DVec3::Y,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        };
        match intersect_plane_sphere(&plane, &sphere) {
            PlaneSphereResult::Circle(c) => {
                let expected_r = (9.0_f64 - 4.0).sqrt();
                assert!((c.radius - expected_r).abs() < TOLERANCE_ABS);
            }
            other => panic!("Expected Circle, got {other:?}"),
        }
    }

    #[test]
    fn tangent() {
        let plane = Plane {
            origin: DVec3::new(0.0, 3.0, 0.0),
            normal: DVec3::Y,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        };
        assert!(matches!(
            intersect_plane_sphere(&plane, &sphere),
            PlaneSphereResult::TangentPoint(_)
        ));
    }

    #[test]
    fn no_intersection() {
        let plane = Plane {
            origin: DVec3::new(0.0, 10.0, 0.0),
            normal: DVec3::Y,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        };
        assert!(matches!(
            intersect_plane_sphere(&plane, &sphere),
            PlaneSphereResult::NoIntersection
        ));
    }
}
