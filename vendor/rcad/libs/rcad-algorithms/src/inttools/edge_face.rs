use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

pub struct EdgeFaceHit {
    pub point: DVec3,
    pub edge_param: f64,
}

/// Intersect a line segment (bounded by t_range) with a plane.
/// Does NOT check face boundary containment — caller must do that.
pub fn intersect_line_plane(line: &Line3, t_range: [f64; 2], plane: &Plane) -> Option<EdgeFaceHit> {
    let denom = line.direction.dot(plane.normal);
    if denom.abs() < TOLERANCE_ABS {
        return None;
    }
    let t = (plane.origin - line.origin).dot(plane.normal) / denom;
    if t < t_range[0] - TOLERANCE_ABS || t > t_range[1] + TOLERANCE_ABS {
        return None;
    }
    let point = line.origin + line.direction * t;
    Some(EdgeFaceHit {
        point,
        edge_param: t,
    })
}

/// Build a local 2D basis on the plane (u, v axes).
pub fn plane_local_basis(plane: &Plane) -> (DVec3, DVec3) {
    let n = plane.normal;
    let ref_dir = if n.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let u = n.cross(ref_dir).normalize();
    let v = n.cross(u).normalize();
    (u, v)
}

/// Check if `point` lies inside a planar face whose boundary vertices are given
/// in order. Uses 2D projection + ray-casting.
pub fn point_in_planar_face(point: DVec3, plane: &Plane, face_verts: &[DVec3]) -> bool {
    if face_verts.len() < 3 {
        return false;
    }
    let (u_axis, v_axis) = plane_local_basis(plane);

    let project = |p: DVec3| -> (f64, f64) {
        let d = p - plane.origin;
        (d.dot(u_axis), d.dot(v_axis))
    };

    let (px, py) = project(point);
    let poly: Vec<(f64, f64)> = face_verts.iter().map(|v| project(*v)).collect();

    ray_cast_contains(px, py, &poly)
}

fn ray_cast_contains(px: f64, py: f64, poly: &[(f64, f64)]) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = poly[i];
        let (xj, yj) = poly[j];
        if ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Clip an infinite line to a convex polygon on a plane.
/// Returns the parametric interval `(t_min, t_max)` of the line inside the polygon,
/// or None if the line doesn't cross the polygon.
pub fn clip_line_to_convex_polygon(
    line: &Line3,
    plane: &Plane,
    face_verts: &[DVec3],
) -> Option<(f64, f64)> {
    if face_verts.len() < 3 {
        return None;
    }
    let (u_axis, v_axis) = plane_local_basis(plane);

    // Project line direction and origin onto 2D
    let line_u = line.direction.dot(u_axis);
    let line_v = line.direction.dot(v_axis);
    let d = line.origin - plane.origin;
    let origin_u = d.dot(u_axis);
    let origin_v = d.dot(v_axis);

    // Determine polygon winding (signed area)
    let n = face_verts.len();
    let pts_2d: Vec<(f64, f64)> = face_verts
        .iter()
        .map(|v| {
            let d = *v - plane.origin;
            (d.dot(u_axis), d.dot(v_axis))
        })
        .collect();

    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += pts_2d[i].0 * pts_2d[j].1 - pts_2d[j].0 * pts_2d[i].1;
    }
    // sign: +1 for CCW, -1 for CW. Inward normal = sign * (ey, -ex)
    let sign = if area >= 0.0 { 1.0 } else { -1.0 };

    // Cyrus-Beck clipping against each polygon edge
    let mut t_enter = f64::NEG_INFINITY;
    let mut t_exit = f64::INFINITY;

    for i in 0..n {
        let j = (i + 1) % n;
        let (ax, ay) = pts_2d[i];
        let (bx, by) = pts_2d[j];

        let ex = bx - ax;
        let ey = by - ay;
        // Inward-facing edge normal: sign * (ey, -ex)
        let nx = sign * ey;
        let ny = sign * (-ex);

        let denom = nx * line_u + ny * line_v;
        let num = nx * (origin_u - ax) + ny * (origin_v - ay);

        if denom.abs() < TOLERANCE_ABS {
            // Line parallel to edge
            if num > TOLERANCE_ABS {
                // Line is outside this edge
                return None;
            }
            // Line is inside or on the edge, continue
        } else {
            let t = -num / denom;
            if denom < 0.0 {
                // Entering
                t_enter = t_enter.max(t);
            } else {
                // Exiting
                t_exit = t_exit.min(t);
            }
        }
    }

    if t_enter > t_exit + TOLERANCE_ABS {
        return None;
    }
    Some((t_enter, t_exit))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_through_plane() {
        let line = Line3 {
            origin: DVec3::new(0.5, 0.5, -1.0),
            direction: DVec3::Z,
        };
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let hit = intersect_line_plane(&line, [-10.0, 10.0], &plane).unwrap();
        assert!((hit.edge_param - 1.0).abs() < TOLERANCE_ABS);
        assert!(points_coincide(hit.point, DVec3::new(0.5, 0.5, 0.0)));
    }

    #[test]
    fn line_parallel_to_plane() {
        let line = Line3 {
            origin: DVec3::new(0.0, 0.0, 1.0),
            direction: DVec3::X,
        };
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        assert!(intersect_line_plane(&line, [-10.0, 10.0], &plane).is_none());
    }

    #[test]
    fn point_inside_square() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let verts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        assert!(point_in_planar_face(
            DVec3::new(0.5, 0.5, 0.0),
            &plane,
            &verts
        ));
        assert!(!point_in_planar_face(
            DVec3::new(1.5, 0.5, 0.0),
            &plane,
            &verts
        ));
    }

    #[test]
    fn clip_line_to_square() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let verts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let line = Line3 {
            origin: DVec3::new(0.5, -1.0, 0.0),
            direction: DVec3::Y,
        };
        let (t_min, t_max) = clip_line_to_convex_polygon(&line, &plane, &verts).unwrap();
        assert!((t_min - 1.0).abs() < 1e-6, "t_min={t_min}");
        assert!((t_max - 2.0).abs() < 1e-6, "t_max={t_max}");
    }
}
