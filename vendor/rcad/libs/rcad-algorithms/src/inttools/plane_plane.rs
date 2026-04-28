use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

#[derive(Debug, Clone)]
pub enum PlanePlaneResult {
    Parallel,
    Coincident,
    Line(Line3),
}

pub fn intersect_plane_plane(p1: &Plane, p2: &Plane) -> PlanePlaneResult {
    let n1 = p1.normal;
    let n2 = p2.normal;
    let cross = n1.cross(n2);

    if is_zero_vec(cross) {
        let d = (p2.origin - p1.origin).dot(n1);
        if d.abs() < TOLERANCE_ABS {
            return PlanePlaneResult::Coincident;
        }
        return PlanePlaneResult::Parallel;
    }

    let direction = cross.normalize();

    let d1 = n1.dot(p1.origin);
    let d2 = n2.dot(p2.origin);
    let origin = solve_two_plane_point(n1, d1, n2, d2, direction);

    PlanePlaneResult::Line(Line3 { origin, direction })
}

/// Find a point on the intersection line of two planes by zeroing the largest
/// component of the line direction and solving the resulting 2x2 system.
fn solve_two_plane_point(n1: DVec3, d1: f64, n2: DVec3, d2: f64, dir: DVec3) -> DVec3 {
    let abs_dir = DVec3::new(dir.x.abs(), dir.y.abs(), dir.z.abs());

    if abs_dir.x >= abs_dir.y && abs_dir.x >= abs_dir.z {
        // Set x = 0
        let det = n1.y * n2.z - n1.z * n2.y;
        let y = (d1 * n2.z - d2 * n1.z) / det;
        let z = (n1.y * d2 - n2.y * d1) / det;
        DVec3::new(0.0, y, z)
    } else if abs_dir.y >= abs_dir.z {
        // Set y = 0
        let det = n1.x * n2.z - n1.z * n2.x;
        let x = (d1 * n2.z - d2 * n1.z) / det;
        let z = (n1.x * d2 - n2.x * d1) / det;
        DVec3::new(x, 0.0, z)
    } else {
        // Set z = 0
        let det = n1.x * n2.y - n1.y * n2.x;
        let x = (d1 * n2.y - d2 * n1.y) / det;
        let y = (n1.x * d2 - n2.x * d1) / det;
        DVec3::new(x, y, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parallel_planes() {
        let p1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let p2 = Plane {
            origin: DVec3::new(0.0, 0.0, 1.0),
            normal: DVec3::Z,
        };
        assert!(matches!(
            intersect_plane_plane(&p1, &p2),
            PlanePlaneResult::Parallel
        ));
    }

    #[test]
    fn coincident_planes() {
        let p1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let p2 = Plane {
            origin: DVec3::new(1.0, 2.0, 0.0),
            normal: DVec3::Z,
        };
        assert!(matches!(
            intersect_plane_plane(&p1, &p2),
            PlanePlaneResult::Coincident
        ));
    }

    #[test]
    fn anti_parallel_coincident() {
        let p1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let p2 = Plane {
            origin: DVec3::ZERO,
            normal: -DVec3::Z,
        };
        assert!(matches!(
            intersect_plane_plane(&p1, &p2),
            PlanePlaneResult::Coincident
        ));
    }

    #[test]
    fn xy_vs_xz_planes() {
        let p1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let p2 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        };
        match intersect_plane_plane(&p1, &p2) {
            PlanePlaneResult::Line(line) => {
                // Intersection should be along X axis
                assert!(vectors_parallel(line.direction, DVec3::X));
                // Point should lie on both planes (z=0 and y=0)
                assert!(line.origin.y.abs() < TOLERANCE_ABS);
                assert!(line.origin.z.abs() < TOLERANCE_ABS);
            }
            _ => panic!("Expected Line"),
        }
    }

    #[test]
    fn arbitrary_planes() {
        let p1 = Plane {
            origin: DVec3::new(0.0, 0.0, 1.0),
            normal: DVec3::Z,
        };
        let p2 = Plane {
            origin: DVec3::new(0.0, 2.0, 0.0),
            normal: DVec3::Y,
        };
        match intersect_plane_plane(&p1, &p2) {
            PlanePlaneResult::Line(line) => {
                // Point must lie on both planes
                let d1 = (line.origin - p1.origin).dot(p1.normal);
                let d2 = (line.origin - p2.origin).dot(p2.normal);
                assert!(d1.abs() < TOLERANCE_ABS, "point not on plane 1: d={d1}");
                assert!(d2.abs() < TOLERANCE_ABS, "point not on plane 2: d={d2}");
                // Direction perpendicular to both normals
                assert!(line.direction.dot(p1.normal).abs() < TOLERANCE_ABS);
                assert!(line.direction.dot(p2.normal).abs() < TOLERANCE_ABS);
            }
            _ => panic!("Expected Line"),
        }
    }
}
