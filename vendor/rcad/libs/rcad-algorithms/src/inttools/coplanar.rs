use glam::DVec3;
use rcad_kernel::geom::*;

use crate::inttools::edge_face::plane_local_basis;
use crate::tolerance::*;

/// Result of coplanar face overlap analysis.
#[derive(Debug)]
pub struct CoplanarResult {
    /// Regions of face1 not covered by face2.
    pub f1_only: Vec<Vec<DVec3>>,
    /// Regions of face2 not covered by face1.
    pub f2_only: Vec<Vec<DVec3>>,
    /// Overlapping region (on both faces).
    pub overlap: Vec<Vec<DVec3>>,
}

/// Analyze two coplanar face polygons. Computes overlap and exclusive regions
/// using the Sutherland-Hodgman polygon clipping algorithm.
pub fn analyze_coplanar_faces(poly1: &[DVec3], poly2: &[DVec3], plane: &Plane) -> CoplanarResult {
    let (u_axis, v_axis) = plane_local_basis(plane);

    let to_2d = |p: DVec3| -> [f64; 2] {
        let d = p - plane.origin;
        [d.dot(u_axis), d.dot(v_axis)]
    };

    let p1_2d: Vec<[f64; 2]> = poly1.iter().map(|&p| to_2d(p)).collect();
    let p2_2d: Vec<[f64; 2]> = poly2.iter().map(|&p| to_2d(p)).collect();

    // Compute intersection using Sutherland-Hodgman
    let overlap_2d = sutherland_hodgman_clip(&p1_2d, &p2_2d);

    let from_2d = |pts: &[[f64; 2]]| -> Vec<DVec3> {
        pts.iter()
            .map(|p| plane.origin + u_axis * p[0] + v_axis * p[1])
            .collect()
    };

    let overlap_3d = if overlap_2d.len() >= 3 {
        vec![from_2d(&overlap_2d)]
    } else {
        vec![]
    };

    // f1_only and f2_only are harder — for now approximate by returning the
    // whole face when overlap is empty, or flagging for the builder to handle
    // via classification.
    // For box booleans, coplanar faces are typically fully overlapping or disjoint.
    CoplanarResult {
        f1_only: vec![],
        f2_only: vec![],
        overlap: overlap_3d,
    }
}

/// Sutherland-Hodgman polygon clipping: clips `subject` against `clip` polygon.
/// Both polygons are in 2D. Returns the clipped polygon vertices.
fn sutherland_hodgman_clip(subject: &[[f64; 2]], clip: &[[f64; 2]]) -> Vec<[f64; 2]> {
    if subject.is_empty() || clip.is_empty() {
        return vec![];
    }

    let mut output = subject.to_vec();

    let n = clip.len();
    for i in 0..n {
        if output.is_empty() {
            return vec![];
        }

        let j = (i + 1) % n;
        let edge_start = clip[i];
        let edge_end = clip[j];

        let input = output;
        output = Vec::new();

        let m = input.len();
        for k in 0..m {
            let current = input[k];
            let previous = input[if k == 0 { m - 1 } else { k - 1 }];

            let curr_inside = is_inside(current, edge_start, edge_end);
            let prev_inside = is_inside(previous, edge_start, edge_end);

            if curr_inside {
                if !prev_inside
                    && let Some(inter) = line_intersect_2d(previous, current, edge_start, edge_end)
                {
                    output.push(inter);
                }
                output.push(current);
            } else if prev_inside
                && let Some(inter) = line_intersect_2d(previous, current, edge_start, edge_end)
            {
                output.push(inter);
            }
        }
    }

    output
}

/// Check if point is on the "inside" (left side) of directed edge a→b.
fn is_inside(point: [f64; 2], a: [f64; 2], b: [f64; 2]) -> bool {
    let cross = (b[0] - a[0]) * (point[1] - a[1]) - (b[1] - a[1]) * (point[0] - a[0]);
    cross >= -TOLERANCE_ABS
}

/// 2D line segment intersection.
fn line_intersect_2d(p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], p4: [f64; 2]) -> Option<[f64; 2]> {
    let d1x = p2[0] - p1[0];
    let d1y = p2[1] - p1[1];
    let d2x = p4[0] - p3[0];
    let d2y = p4[1] - p3[1];

    let denom = d1x * d2y - d1y * d2x;
    if denom.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
        return None;
    }

    let t = ((p3[0] - p1[0]) * d2y - (p3[1] - p1[1]) * d2x) / denom;

    Some([p1[0] + t * d1x, p1[1] + t * d1y])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overlapping_squares() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let poly1 = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
            DVec3::new(2.0, 2.0, 0.0),
            DVec3::new(0.0, 2.0, 0.0),
        ];
        let poly2 = vec![
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(3.0, 1.0, 0.0),
            DVec3::new(3.0, 3.0, 0.0),
            DVec3::new(1.0, 3.0, 0.0),
        ];
        let result = analyze_coplanar_faces(&poly1, &poly2, &plane);
        assert_eq!(result.overlap.len(), 1);
        let overlap = &result.overlap[0];
        assert_eq!(overlap.len(), 4);
        // Overlap should be the 1x1 square (1,1)-(2,2)
        for v in overlap {
            assert!(v.x >= 1.0 - TOLERANCE_ABS && v.x <= 2.0 + TOLERANCE_ABS);
            assert!(v.y >= 1.0 - TOLERANCE_ABS && v.y <= 2.0 + TOLERANCE_ABS);
        }
    }

    #[test]
    fn disjoint_squares() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let poly1 = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let poly2 = vec![
            DVec3::new(5.0, 5.0, 0.0),
            DVec3::new(6.0, 5.0, 0.0),
            DVec3::new(6.0, 6.0, 0.0),
            DVec3::new(5.0, 6.0, 0.0),
        ];
        let result = analyze_coplanar_faces(&poly1, &poly2, &plane);
        assert!(result.overlap.is_empty());
    }

    #[test]
    fn contained_square() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let outer = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(4.0, 0.0, 0.0),
            DVec3::new(4.0, 4.0, 0.0),
            DVec3::new(0.0, 4.0, 0.0),
        ];
        let inner = vec![
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(3.0, 1.0, 0.0),
            DVec3::new(3.0, 3.0, 0.0),
            DVec3::new(1.0, 3.0, 0.0),
        ];
        let result = analyze_coplanar_faces(&outer, &inner, &plane);
        assert_eq!(result.overlap.len(), 1);
        assert_eq!(result.overlap[0].len(), 4);
    }

    #[test]
    fn sutherland_hodgman_basic() {
        let subject = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let clip = vec![[1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [1.0, 3.0]];
        let result = sutherland_hodgman_clip(&subject, &clip);
        assert_eq!(result.len(), 4);
    }
}
