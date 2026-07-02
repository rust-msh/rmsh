//! Bowyer-Watson Delaunay triangulation kernel.
//!
//! [`triangulate_points`] is the core Delaunay triangulation routine.
//! Polygon meshing (`mesh_polygon`, `Polygon2D`) lives in
//! [`crate::delaunay_2d`].

use rmsh_model::{ElementType, Mesh, Node};

/// Bowyer-Watson incremental Delaunay triangulation.
///
/// Returns a list of triangles as `[i, j, k]` index triples into `pts`.
/// All points must be distinct (within floating-point tolerance).
pub fn triangulate_points(pts: &[[f64; 2]]) -> Vec<[usize; 3]> {
    let n = pts.len();
    if n < 3 {
        return Vec::new();
    }

    let (mut min_x, mut min_y) = (f64::MAX, f64::MAX);
    let (mut max_x, mut max_y) = (f64::MIN, f64::MIN);
    for p in pts {
        min_x = min_x.min(p[0]);
        min_y = min_y.min(p[1]);
        max_x = max_x.max(p[0]);
        max_y = max_y.max(p[1]);
    }
    let dx = (max_x - min_x).max(1e-9);
    let dy = (max_y - min_y).max(1e-9);
    let d = dx.max(dy);
    let mx = (min_x + max_x) / 2.0;
    let my = (min_y + max_y) / 2.0;

    let st0 = [mx - 20.0 * d, my - d];
    let st1 = [mx, my + 20.0 * d];
    let st2 = [mx + 20.0 * d, my - d];

    let mut all: Vec<[f64; 2]> = pts.to_vec();
    let st_start = all.len();
    all.push(st0);
    all.push(st1);
    all.push(st2);

    let mut triangles: Vec<[usize; 3]> = vec![[st_start, st_start + 1, st_start + 2]];

    for i in 0..n {
        let p = pts[i];
        let bad: Vec<[usize; 3]> = triangles
            .iter()
            .filter(|&&tri| circumcircle_contains(all[tri[0]], all[tri[1]], all[tri[2]], p))
            .copied()
            .collect();

        let mut boundary: Vec<[usize; 2]> = Vec::new();
        for &tri in &bad {
            let edges = [[tri[0], tri[1]], [tri[1], tri[2]], [tri[2], tri[0]]];
            for edge in edges {
                let shared = bad
                    .iter()
                    .any(|&other| other != tri && tri_has_edge(other, edge[0], edge[1]));
                if !shared {
                    boundary.push(edge);
                }
            }
        }

        triangles.retain(|t| !bad.contains(t));

        for edge in boundary {
            triangles.push([edge[0], edge[1], i]);
        }
    }

    triangles.retain(|t| t[0] < st_start && t[1] < st_start && t[2] < st_start);
    triangles
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

/// Returns `true` if the circumcircle of triangle `(a, b, c)` contains `p`.
fn circumcircle_contains(a: [f64; 2], b: [f64; 2], c: [f64; 2], p: [f64; 2]) -> bool {
    let ax = a[0] - p[0];
    let ay = a[1] - p[1];
    let bx = b[0] - p[0];
    let by = b[1] - p[1];
    let cx = c[0] - p[0];
    let cy = c[1] - p[1];

    let a2 = ax * ax + ay * ay;
    let b2 = bx * bx + by * by;
    let c2 = cx * cx + cy * cy;

    let det = ax * (by * c2 - cy * b2) - ay * (bx * c2 - cx * b2) + a2 * (bx * cy - by * cx);
    let orient = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    det * orient > 0.0
}

fn tri_has_edge(tri: [usize; 3], a: usize, b: usize) -> bool {
    [[tri[0], tri[1]], [tri[1], tri[2]], [tri[2], tri[0]]]
        .iter()
        .any(|e| (e[0] == a && e[1] == b) || (e[0] == b && e[1] == a))
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_square_delaunay() {
        let pts = [[0.0f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let tris = triangulate_points(&pts);
        assert!(!tris.is_empty());
        for t in &tris { assert!(t[0] < pts.len() && t[1] < pts.len() && t[2] < pts.len()); }
    }

    #[test]
    fn fewer_than_three_points_returns_empty() {
        assert!(triangulate_points(&[]).is_empty());
        assert!(triangulate_points(&[[0.0, 0.0]]).is_empty());
        assert!(triangulate_points(&[[0.0, 0.0], [1.0, 0.0]]).is_empty());
    }

    #[test]
    fn exactly_three_points_gives_one_triangle() {
        let pts = [[0.0f64, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let tris = triangulate_points(&pts);
        assert_eq!(tris.len(), 1);
        let mut sorted = tris[0];
        sorted.sort();
        assert_eq!(sorted, [0, 1, 2]);
    }

    #[test]
    fn all_output_indices_are_valid() {
        let pts: Vec<[f64; 2]> = (0..20)
            .map(|i| { let a = i as f64 * std::f64::consts::TAU / 20.0; [a.cos(), a.sin()] })
            .collect();
        let tris = triangulate_points(&pts);
        for t in &tris { for &idx in t { assert!(idx < pts.len()); } }
    }

    #[test]
    fn no_degenerate_triangles_in_output() {
        let pts: Vec<[f64; 2]> = (0..15)
            .map(|i| { let a = i as f64 * std::f64::consts::TAU / 15.0; [a.cos(), a.sin()] })
            .collect();
        let tris = triangulate_points(&pts);
        for t in &tris {
            let a = pts[t[0]]; let b = pts[t[1]]; let c = pts[t[2]];
            let area2 = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            assert!(area2.abs() > 1e-10, "degenerate triangle");
        }
    }

    #[test]
    fn no_duplicate_triangles_in_output() {
        let pts: Vec<[f64; 2]> = (0..12)
            .map(|i| { let a = i as f64 * std::f64::consts::TAU / 12.0; [a.cos(), a.sin()] })
            .collect();
        let tris = triangulate_points(&pts);
        let mut sorted_tris: Vec<[usize; 3]> =
            tris.iter().map(|t| { let mut s = *t; s.sort(); s }).collect();
        sorted_tris.sort();
        let orig_len = sorted_tris.len();
        sorted_tris.dedup();
        assert_eq!(sorted_tris.len(), orig_len, "duplicate triangles");
    }

    #[test]
    fn all_input_points_appear_in_output_triangulation() {
        let pts: Vec<[f64; 2]> = (0..10)
            .map(|i| { let a = i as f64 * std::f64::consts::TAU / 10.0; [a.cos(), a.sin()] })
            .collect();
        let tris = triangulate_points(&pts);
        let used: std::collections::HashSet<usize> = tris.iter().flatten().copied().collect();
        for i in 0..pts.len() { assert!(used.contains(&i), "point {i} not used"); }
    }
}
