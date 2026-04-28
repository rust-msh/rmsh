use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

/// V-V: test if two vertices coincide.
pub fn vertex_vertex_coincide(p1: DVec3, p2: DVec3) -> bool {
    points_coincide(p1, p2)
}

/// V-E: test if a vertex lies on a line edge within its parametric range.
/// Returns the parameter if it does.
pub fn vertex_on_line(point: DVec3, line: &Line3, t_range: [f64; 2]) -> Option<f64> {
    let v = point - line.origin;
    let t = v.dot(line.direction);
    let closest = line.origin + line.direction * t;
    if !points_coincide(closest, point) {
        return None;
    }
    if t < t_range[0] - TOLERANCE_ABS || t > t_range[1] + TOLERANCE_ABS {
        return None;
    }
    Some(t)
}

/// V-F: test if a vertex lies on a plane (within tolerance).
pub fn vertex_on_plane(point: DVec3, plane: &Plane) -> bool {
    let d = (point - plane.origin).dot(plane.normal);
    d.abs() < TOLERANCE_ABS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vertex_on_line_segment() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let t = vertex_on_line(DVec3::new(0.5, 0.0, 0.0), &line, [0.0, 1.0]);
        assert!(t.is_some());
        assert!((t.unwrap() - 0.5).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn vertex_off_line() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        assert!(vertex_on_line(DVec3::new(0.5, 1.0, 0.0), &line, [0.0, 1.0]).is_none());
    }

    #[test]
    fn vertex_on_plane_test() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        assert!(vertex_on_plane(DVec3::new(5.0, 3.0, 0.0), &plane));
        assert!(!vertex_on_plane(DVec3::new(5.0, 3.0, 1.0), &plane));
    }
}
