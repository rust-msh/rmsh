//! 2-D mesh quality optimization: the Gmsh surface-mesh smoothing scheme.
//!
//! # Algorithm overview
//!
//! Two quality measures and two passes, following Gmsh's 2-D optimization
//! machinery (`Mesh/meshGFaceOptimize.cpp`, `Mesh/qualityMeasures.cpp`):
//!
//! * **eta shape measure** (`qmTriangle::eta`): `1 - |60° - amin| / 60°`
//!   over the minimum interior angle `amin`. 1 for equilateral, 0 degenerate.
//! * **gamma shape measure** (`qmTriangle::gamma`): `2·r/R` (twice the
//!   inradius over circumradius), computed as
//!   `2·(2·sina·sinb·sinc)/(sina+sinb+sinc)`.
//!
//! * **Relocation pass** (`_relocate` / `_isItAGoodIdeaToMoveThatVertex`):
//!   every free (interior) vertex is pulled toward the centroid of its
//!   neighbours; the move is tried at five damped fractions
//!   (`FACTOR = 1, 1/1.4, 1/1.96, ...`) and accepted only when the total
//!   surface area of the vertex star does not grow (relative tolerance
//!   `1e-10`) and the star's minimum eta does not decrease. Accepted
//!   damped steps accumulate within the same vertex visit.
//! * **Delaunay-restore swap pass** (`PolyMeshDelaunayize`): flip every
//!   interior edge that fails the in-circumcircle test (the Delaunay
//!   triangulation is the max-min-angle triangulation of a planar point
//!   set). Constraint edges are never swapped.
//!
//! # Reference
//!
//! Gmsh source: `Mesh/meshGFaceOptimize.cpp` (`_relocate`, `vertexStar`,
//! `_isItAGoodIdeaToMoveThatVertex`, `laplaceSmoothing`),
//! `Mesh/qualityMeasures.cpp` (`qmTriangle::eta`, `qmTriangle::gamma`).
//! C. Geuzaine, J.-F. Remacle, "Gmsh", *Int. J. Numer. Methods Eng.* 79,
//! 2009 (the two-quality-threshold Laplacian + swap scheme follows Gmsh's
//! meshGFaceOptimize approach; reimplemented, not translated).

use crate::cdt_2d::TriMesh2d;

/// Minimum interior angle of a triangle in degrees.
fn min_angle_deg(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    angle_deg(a, b, c).min(angle_deg(b, c, a)).min(angle_deg(c, a, b))
}

/// Angle at vertex `p` between the directions to `q` and `r` (radians).
fn angle_rad(p: [f64; 2], q: [f64; 2], r: [f64; 2]) -> f64 {
    let u = [q[0] - p[0], q[1] - p[1]];
    let v = [r[0] - p[0], r[1] - p[1]];
    let lu = (u[0] * u[0] + u[1] * u[1]).sqrt();
    let lv = (v[0] * v[0] + v[1] * v[1]).sqrt();
    if lu < 1e-15 || lv < 1e-15 {
        return 0.0;
    }
    ((u[0] * v[0] + u[1] * v[1]) / (lu * lv)).clamp(-1.0, 1.0).acos()
}

fn angle_deg(p: [f64; 2], q: [f64; 2], r: [f64; 2]) -> f64 {
    angle_rad(p, q, r).to_degrees()
}

/// Gmsh `qmTriangle::eta`: `1 - |60° - amin| / 60` over the minimum interior
/// angle. 1 for the equilateral triangle, 0 for a degenerate one.
pub(crate) fn eta_shape_measure(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let amin = min_angle_deg(a, b, c);
    1.0 - (60.0 - amin).abs() / 60.0
}

/// Gmsh `qmTriangle::gamma`: `2 r / R` computed through the sines of the
/// corner angles: `2 · (2·sina·sinb·sinc) / (sina + sinb + sinc)` — the
/// swap criterion of the gamma-driven pass.
pub(crate) fn gamma_shape_measure(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let sina = angle_rad(b, a, c).sin();
    let sinb = angle_rad(c, b, a).sin();
    let sinc = angle_rad(a, c, b).sin();
    if sina == 0.0 && sinb == 0.0 && sinc == 0.0 {
        return 0.0;
    }
    2.0 * (2.0 * sina * sinb * sinc / (sina + sinb + sinc))
}

/// Twice the signed area (positive = CCW).
fn orient(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// One Gmsh `laplaceSmoothing` pass over the live mesh: relocate every free
/// vertex toward its neighbour centroid with the damped, quality-gain-checked
/// acceptance of `_relocate` / `_isItAGoodIdeaToMoveThatVertex`.
///
/// `free` marks the vertices allowed to move (interior, non-pinned); all
/// other vertices stay fixed. Returns the number of accepted moves.
pub(crate) fn relocate_pass(mesh: &mut TriMesh2d, free: &std::collections::HashSet<usize>) -> usize {
    let mut accepted_total = 0usize;
    // vertices in increasing index order (Gmsh accumulates in vertex-num order)
    let mut vs: Vec<usize> = mesh.live_vertices();
    vs.sort_unstable();
    // the incident-triangle adjacency, built once for the whole pass
    let mut star_of: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    for t in mesh.live() {
        for &v in &mesh.tri[t] {
            star_of.entry(v).or_default().push(t);
        }
    }
    for v in vs {
        if !free.contains(&v) {
            continue;
        }
        // neighbours from the incident triangles (two occurrences per edge,
        // deduplicated in place)
        let star = match star_of.get(&v) {
            Some(s) => s,
            None => continue,
        };
        let mut neighbours: Vec<usize> = Vec::with_capacity(star.len() * 2);
        for &t in star {
            for &w in &mesh.tri[t] {
                if w != v {
                    neighbours.push(w);
                }
            }
        }
        neighbours.sort_unstable();
        neighbours.dedup();
        if neighbours.len() < 3 {
            continue;
        }
        let before = mesh.pts[v];
        // plain centroid target (the planar Euclidean metric reduces Gmsh's
        // metric-weighted accumulation to the uniform centroid)
        let (mut sx, mut sy) = (0.0f64, 0.0f64);
        for &n in &neighbours {
            sx += mesh.pts[n][0];
            sy += mesh.pts[n][1];
        }
        let count = neighbours.len() as f64;
        let after = [sx / count, sy / count];
        // star of the vertex: incident triangles
        let star = match star_of.get(&v) {
            Some(s) => s.clone(),
            None => continue,
        };
        // current star surface area and minimum eta
        let mut surface = 0.0f64;
        let mut minq = f64::INFINITY;
        for &t in &star {
            let [a, b, c] = mesh.tri[t];
            surface += orient(mesh.pts[a], mesh.pts[b], mesh.pts[c]).abs() * 0.5;
            minq = minq.min(eta_shape_measure(mesh.pts[a], mesh.pts[b], mesh.pts[c]));
        }
        // damped trial ladder (MAXITER 5, FACTOR /= 1.4); accepted steps
        // accumulate in `actual`, and the running (surface, minq) move with
        // them
        let mut actual = before;
        let mut factor = 1.0f64;
        for _ in 0..5 {
            let trial = [
                after[0] * factor + actual[0] * (1.0 - factor),
                after[1] * factor + actual[1] * (1.0 - factor),
            ];
            // evaluate the star at the trial position
            let mut new_surface = 0.0f64;
            let mut new_minq = f64::INFINITY;
            let mut inverted = false;
            for &t in &star {
                let corners: Vec<[f64; 2]> = mesh.tri[t]
                    .iter()
                    .map(|&w| if w == v { trial } else { mesh.pts[w] })
                    .collect();
                let o = orient(corners[0], corners[1], corners[2]);
                if o <= 0.0 {
                    inverted = true;
                    break;
                }
                new_surface += o * 0.5;
                new_minq = new_minq.min(eta_shape_measure(corners[0], corners[1], corners[2]));
            }
            if !inverted && (1.0 + 1e-10) * surface >= new_surface && new_minq >= minq {
                actual = trial;
                surface = new_surface;
                minq = new_minq;
                accepted_total += 1;
            }
            factor /= 1.4;
        }
        if actual != before {
            mesh.pts[v] = actual;
        }
    }
    accepted_total
}

/// One quality-driven swap pass over the live mesh (Gmsh `edgeSwapTest` /
/// `swapEdgePass` in `Mesh/meshGFaceBDS.cpp`): every interior edge whose
/// adjacent quad is strictly convex is swapped when the swap strictly
/// improves the minimum gamma of the two triangles (the swap is not limited
/// to Delaunay-restoring flips — quality-improving swaps also reshape
/// boundary-anchored flat spans that the circumcircle criterion leaves
/// alone). Constraint edges are never swapped. Returns the number of swaps.
pub(crate) fn gamma_swap_pass(mesh: &mut TriMesh2d) -> usize {
    let mut swaps = 0usize;
    for _pass in 0..8 {
        let mut any = false;
        // interior edges of the live mesh with their adjacent triangles,
        // built once per pass (deterministic order)
        let mut edges: Vec<((usize, usize), (usize, usize), (usize, usize))> = Vec::new();
        let mut seen: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();
        for t in mesh.live() {
            for k in 0..3 {
                let (u, v) = (mesh.tri[t][k], mesh.tri[t][(k + 1) % 3]);
                let key = if u < v { (u, v) } else { (v, u) };
                if mesh.constraint.contains(&crate::cdt_2d::edge_key(u, v))
                    || mesh.neigh[t][k] < 0
                    || !seen.insert(key)
                {
                    continue;
                }
                let apex = mesh
                    .tri[t]
                    .iter()
                    .copied()
                    .find(|&x| x != u && x != v)
                    .expect("edge endpoints are triangle vertices");
                let n = mesh.neigh[t][k] as usize;
                let apex1 = mesh
                    .tri[n]
                    .iter()
                    .copied()
                    .find(|&x| x != u && x != v)
                    .expect("edge endpoints are triangle vertices");
                edges.push((key, (t, apex), (n, apex1)));
            }
        }
        for ((u, v), (t0, w), (t1, x)) in edges {
            if mesh.dead[t0] || mesh.dead[t1] || t0 == t1 {
                continue;
            }
            if !mesh_is_quad_convex(mesh.pts[u], mesh.pts[w], mesh.pts[v], mesh.pts[x]) {
                continue;
            }
            let gamma = |tri: [usize; 3]| {
                gamma_shape_measure(mesh.pts[tri[0]], mesh.pts[tri[1]], mesh.pts[tri[2]])
            };
            let qa = gamma(mesh.tri[t0]).min(gamma(mesh.tri[t1]));
            let mut n0 = [u, w, x];
            if mesh_orient(mesh.pts[n0[0]], mesh.pts[n0[1]], mesh.pts[n0[2]]) < 0.0 {
                n0 = [u, x, w];
            }
            let mut n1 = [w, v, x];
            if mesh_orient(mesh.pts[n1[0]], mesh.pts[n1[1]], mesh.pts[n1[2]]) < 0.0 {
                n1 = [w, x, v];
            }
            let qb = gamma(n0).min(gamma(n1));
            if qb > qa + 1e-12 {
                mesh.flip(t0, t1, u, v, w, x);
                swaps += 1;
                any = true;
            }
        }
        if !any {
            break;
        }
    }
    swaps
}

fn mesh_is_quad_convex(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let o = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| {
        (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    };
    (o(a, b, c) > 0.0 && o(b, c, d) > 0.0 && o(c, d, a) > 0.0 && o(d, a, b) > 0.0)
        || (o(a, b, c) < 0.0 && o(b, c, d) < 0.0 && o(c, d, a) < 0.0 && o(d, a, b) < 0.0)
}

fn mesh_orient(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eta_measures_are_one_for_equilateral() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        // equilateral triangle of side 1: height sqrt(3)/2
        let c = [0.5, 0.5 * 3.0f64.sqrt()];
        assert!((eta_shape_measure(a, b, c) - 1.0).abs() < 1e-12);
        assert!((gamma_shape_measure(a, b, c) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn eta_is_low_for_a_sliver() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.5, 0.01];
        let e = eta_shape_measure(a, b, c);
        assert!(e < 0.1, "sliver eta {e} too high");
        let g = gamma_shape_measure(a, b, c);
        assert!(g < 0.1, "sliver gamma {g} too high");
    }

    #[test]
    fn min_angle_handles_degenerate_input() {
        assert_eq!(min_angle_deg([0.0, 0.0], [0.0, 0.0], [1.0, 1.0]), 0.0);
        assert_eq!(eta_shape_measure([0.0, 0.0], [0.0, 0.0], [1.0, 1.0]), 0.0);
    }
}
