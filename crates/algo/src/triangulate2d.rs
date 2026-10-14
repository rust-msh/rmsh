//! Simple 2-D mesh generation via Bowyer-Watson Delaunay triangulation.
//!
//! # Example
//! ```rust
//! use rmsh_algo::triangulate2d::{Polygon2D, mesh_polygon};
//!
//! // A unit-square polygon
//! let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
//! let mesh = mesh_polygon(&poly, 0.25).unwrap();
//! assert!(mesh.element_count() > 0);
//! ```

use std::collections::HashMap;

use crate::cdt_2d::edge_key;
use rmsh_model::{Element, ElementType, Mesh, Node};
use thiserror::Error;

/// Error type for 2-D mesh generation.
#[derive(Error, Debug)]
pub enum MeshError {
    #[error("Mesh generation failed: {0}")]
    Generation(String),
}

/// A 2-D polygon defined by an ordered list of boundary vertices.
pub struct Polygon2D {
    pub vertices: Vec<[f64; 2]>,
}

impl Polygon2D {
    /// Create a new polygon from a list of 2-D vertices (CCW or CW order).
    pub fn new(vertices: Vec<[f64; 2]>) -> Self {
        Self { vertices }
    }

    /// Point-in-polygon test using the ray casting algorithm.
    pub fn contains(&self, p: [f64; 2]) -> bool {
        let n = self.vertices.len();
        let (px, py) = (p[0], p[1]);
        let mut inside = false;
        let mut j = n - 1;
        for i in 0..n {
            let (xi, yi) = (self.vertices[i][0], self.vertices[i][1]);
            let (xj, yj) = (self.vertices[j][0], self.vertices[j][1]);
            if ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi) {
                inside = !inside;
            }
            j = i;
        }
        inside
    }

    /// Axis-aligned bounding box: returns `(min, max)`.
    pub fn bounding_box(&self) -> ([f64; 2], [f64; 2]) {
        let mut min = [f64::MAX; 2];
        let mut max = [f64::MIN; 2];
        for v in &self.vertices {
            min[0] = min[0].min(v[0]);
            min[1] = min[1].min(v[1]);
            max[0] = max[0].max(v[0]);
            max[1] = max[1].max(v[1]);
        }
        (min, max)
    }
}

// ─── Delaunay triangulation (Lawson incremental) ─────────────────────────────

/// Incremental Delaunay triangulation (Lawson construction with exact
/// predicates).
///
/// Returns a list of triangles as `[i, j, k]` index triples into `pts`.
/// All points must be distinct (within floating-point tolerance).
///
/// The computation runs on the point set translated and scaled into the
/// unit box: the super-triangle is sized relative to the bbox, and large /
/// offset coordinates make the geometric determinants lose precision
/// through cancellation. Index triples are unaffected. Points inserted
/// exactly onto an existing edge split it; exactly cocircular ties are left
/// in place (weakly Delaunay). The output is a valid triangulation of the
/// point set — no overlapping triangles — with the Delaunay property where
/// the input is non-degenerate.
pub fn triangulate_points(pts: &[[f64; 2]]) -> Vec<[usize; 3]> {
    let n = pts.len();
    if n < 3 {
        return Vec::new();
    }

    // Translate + scale into the unit box (indices unchanged)
    let (min_x, min_y) = (f64::MAX, f64::MAX);
    let (max_x, max_y) = (f64::MIN, f64::MIN);
    let (min_x, min_y) = pts.iter().fold((min_x, min_y), |(a, b), p| {
        (a.min(p[0]), b.min(p[1]))
    });
    let (max_x, max_y) = pts.iter().fold((max_x, max_y), |(a, b), p| {
        (a.max(p[0]), b.max(p[1]))
    });
    let scale = (max_x - min_x).max(max_y - min_y).max(1e-12);
    let scaled: Vec<[f64; 2]> = pts
        .iter()
        .map(|p| [(p[0] - min_x) / scale, (p[1] - min_y) / scale])
        .collect();

    // Compute bounding box for super-triangle sizing
    let (mut min_x, mut min_y) = (f64::MAX, f64::MAX);
    let (mut max_x, mut max_y) = (f64::MIN, f64::MIN);
    for p in &scaled {
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

    // Super-triangle covering all input points
    let st0 = [mx - 20.0 * d, my - d];
    let st1 = [mx, my + 20.0 * d];
    let st2 = [mx + 20.0 * d, my - d];

    let mut all: Vec<[f64; 2]> = scaled;
    let st_start = all.len();
    all.push(st0);
    all.push(st1);
    all.push(st2);

    /// The unordered edge key of (u, v).
    /// The CCW edge slot carrying the unordered edge (u, v).
    fn edge_slot_of(tri: &[usize; 3], u: usize, v: usize) -> Option<usize> {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            if (a == u && b == v) || (a == v && b == u) {
                return Some(k);
            }
        }
        None
    }

    // ── Lawson incremental construction ──────────────────────────────────
    // Every point is inserted into its containing triangle (a fan of 3, or
    // an exact edge split when the point lies EXACTLY on an edge) and the
    // result is restored to the Delaunay property by stack-driven
    // legalization flips (Lawson 1977, in the robust form of Shewchuk 1997
    // with the exact orientation / in-circumcircle predicates of
    // `crate::predicates`). Unlike the Bowyer-Watson cavity scheme this
    // replaces, validity never depends on a cavity being star-shaped: a
    // flip is taken only on a strictly convex quad, so collinear chains and
    // cocircular ties leave a valid weakly-Delaunay triangulation instead
    // of overlapping triangles (the concave-domain leak the P3 quality card
    // traced to cavity ties on the super-triangle circles).
    //
    // Per-triangle neighbour links (`nbr[ti][k]` = the live carrier across
    // the CCW edge (tri[k], tri[k+1]), -1 = none) are maintained through
    // every mutation: the point location walks from the previous host
    // (expected sub-linear), and the legalization carrier lookup reads the
    // link in O(1) (full-scan fallbacks preserve the old behaviour on any
    // stale hint). Both were full O(T) scans and dominated the build time
    // (measured 94% at n = 5000).
    let mut triangles: Vec<[usize; 3]> = vec![tri_ccw([st_start, st_start + 1, st_start + 2], &all)];
    let mut nbr: Vec<[i32; 3]> = vec![[-1; 3]];
    // vertex -> incident live triangles (Gmsh's adjv): the legalize
    // carrier lookup scans one endpoint's star instead of the whole
    // triangle vector
    let mut adjv: Vec<Vec<usize>> = vec![Vec::new(); all.len()];
    for v in &triangles[0] {
        adjv[*v].push(0);
    }

    // Rewire the neighbourhoods after replacing `dead` triangles (index +
    // old value + old links) with `new` ones (index + value): interior
    // edges pair the new triangles mutually; the surviving external edges
    // relink to their old neighbours and fix the neighbours' back links.
    // In-place index reuse keeps most back links valid; the back-link fix
    // handles the slots that moved.
    fn patch_neighbours(
        nbr: &mut Vec<[i32; 3]>,
        adjv: &mut [Vec<usize>],
        dead: &[(usize, [usize; 3], [i32; 3])],
        new: &[(usize, [usize; 3])],
    ) {
        // star maintenance: the replaced triangles leave their corners'
        // stars, the new ones enter them
        for (di, dtri, _) in dead {
            for &v in dtri {
                adjv[v].retain(|t| *t != *di);
            }
        }
        for (ni, ntri) in new {
            for &v in ntri {
                adjv[v].push(*ni);
            }
        }
        // grow alongside the triangle vector (new slots start unlinked)
        let need = new.iter().map(|(i, _)| i + 1).max().unwrap_or(0);
        if nbr.len() < need {
            nbr.resize(need, [-1; 3]);
        }
        let mut external: std::collections::HashMap<(usize, usize), (usize, usize)> =
            std::collections::HashMap::new();
        for (di, dtri, _dn) in dead {
            for k in 0..3 {
                external.insert(edge_key(dtri[k], dtri[(k + 1) % 3]), (*di, k));
            }
        }
        for (_ni, ntri) in new {
            for k in 0..3 {
                let key = edge_key(ntri[k], ntri[(k + 1) % 3]);
                external.remove(&key);
            }
        }
        // interior pairing among the new triangles
        for (ni, ntri) in new {
            for k in 0..3 {
                for (nj, ntri2) in new {
                    if nj == ni {
                        continue;
                    }
                    if let Some(k2) = edge_slot_of(ntri2, ntri[k], ntri[(k + 1) % 3]) {
                        nbr[*ni][k] = *nj as i32;
                        nbr[*nj][k2] = *ni as i32;
                        break;
                    }
                }
            }
        }
        // external relinks
        for (ni, ntri) in new {
            for k in 0..3 {
                let key = edge_key(ntri[k], ntri[(k + 1) % 3]);
                if let Some(&(di, dk)) = external.get(&key) {
                    let Some((_, _, dn)) = dead.iter().find(|(di2, _, _)| *di2 == di) else {
                        continue;
                    };
                    let old_nb = dn[dk];
                    nbr[*ni][k] = old_nb;
                    if old_nb >= 0 {
                        let on = old_nb as usize;
                        if let Some(kk) = (0..3).find(|&kk| nbr[on][kk] == di as i32) {
                            nbr[on][kk] = *ni as i32;
                        }
                    }
                }
            }
        }
    }

    let mut hint = 0usize;
    for i in 0..n {
        // locate: walk from the hint triangle across the first violated
        // edge (deterministic CCW-first rule; the mesh is weakly Delaunay
        // between insertions, so the walk terminates), falling back to the
        // exhaustive scan on a tripped guard
        let mut host: Option<(usize, [i32; 3])> = None;
        {
            let mut t = hint;
            let mut guard = 0usize;
            loop {
                guard += 1;
                if guard > 4 * triangles.len() + 64 {
                    break;
                }
                let tri = triangles[t];
                let mut signs = [0i32; 3];
                let mut violated = None;
                for k in 0..3 {
                    let s = crate::predicates::orient2d_sign(
                        all[tri[k]],
                        all[tri[(k + 1) % 3]],
                        all[i],
                    );
                    signs[k] = s;
                    if s < 0 {
                        violated = Some(k);
                        break;
                    }
                }
                match violated {
                    None => {
                        host = Some((t, signs));
                        break;
                    }
                    Some(k) => {
                        let nx = nbr[t][k];
                        if nx < 0 {
                            break; // walked off the hull: exhaustive fallback
                        }
                        t = nx as usize;
                    }
                }
            }
        }
        if host.is_none() {
            for (ti, tri) in triangles.iter().enumerate() {
                let s = [
                    crate::predicates::orient2d_sign(all[tri[0]], all[tri[1]], all[i]),
                    crate::predicates::orient2d_sign(all[tri[1]], all[tri[2]], all[i]),
                    crate::predicates::orient2d_sign(all[tri[2]], all[tri[0]], all[i]),
                ];
                if s[0] >= 0 && s[1] >= 0 && s[2] >= 0 {
                    host = Some((ti, s));
                    break;
                }
            }
        }
        let Some((ti, s)) = host else {
            // unreachable for distinct inputs: the scaled point set lies in
            // the unit box, which the super-triangle covers with margin
            continue;
        };
        hint = ti;
        let tri = triangles[ti];
        let mut new_tris: Vec<(usize, [usize; 3])> = Vec::with_capacity(4);
        if s[0] == 0 || s[1] == 0 || s[2] == 0 {
            // exactly on edge k (exact collinearity, e.g. an axis-aligned
            // boundary chain): split the edge in both carrying triangles
            let k = if s[0] == 0 {
                0
            } else if s[1] == 0 {
                1
            } else {
                2
            };
            let (u, v) = (tri[k], tri[(k + 1) % 3]);
            let x = tri[(k + 2) % 3];
            let dot = (all[i][0] - all[u][0]) * (all[i][0] - all[v][0])
                + (all[i][1] - all[u][1]) * (all[i][1] - all[v][1]);
            if dot >= 0.0 {
                // i coincides with u or v: a duplicate input point, already
                // a vertex
                continue;
            }
            let t0 = tri_ccw([u, i, x], &all);
            let t1 = tri_ccw([i, v, x], &all);
            // the neighbour across (u, v): the link, validated
            let oj_find = {
                let o = nbr[ti][k];
                if o >= 0 && tri_has_edge(triangles[o as usize], u, v) {
                    Some(o as usize)
                } else {
                    (0..triangles.len())
                        .find(|&oj| oj != ti && tri_has_edge(triangles[oj], u, v))
                }
            };
            triangles[ti] = t0;
            new_tris.push((ti, t0));
            triangles.push(t1);
            new_tris.push((triangles.len() - 1, t1));
            if let Some(oj) = oj_find {
                let ntri = triangles[oj];
                let y = ntri
                    .iter()
                    .copied()
                    .find(|&q| q != u && q != v)
                    .expect("edge endpoints are triangle vertices");
                let t2 = tri_ccw([u, y, i], &all);
                let t3 = tri_ccw([i, y, v], &all);
                let dead_nbr = nbr[oj];
                triangles[oj] = t2;
                new_tris.push((oj, t2));
                triangles.push(t3);
                new_tris.push((triangles.len() - 1, t3));
                let dn0 = nbr[ti];
                patch_neighbours(
                    &mut nbr,
                    &mut adjv,
                    &[
                        (ti, tri, dn0),
                        (oj, ntri, dead_nbr),
                    ],
                    &new_tris,
                );
            } else {
                let dn0 = nbr[ti];
                patch_neighbours(&mut nbr, &mut adjv, &[(ti, tri, dn0)], &new_tris);
            }
        } else {
            // strictly inside: replace the host with the fan of 3
            let t0 = tri_ccw([tri[0], tri[1], i], &all);
            let t1 = tri_ccw([tri[1], tri[2], i], &all);
            let t2 = tri_ccw([tri[2], tri[0], i], &all);
            let dead_nbr = nbr[ti];
            triangles[ti] = t0;
            new_tris.push((ti, t0));
            triangles.push(t1);
            new_tris.push((triangles.len() - 1, t1));
            triangles.push(t2);
            new_tris.push((triangles.len() - 1, t2));
            patch_neighbours(&mut nbr, &mut adjv, &[(ti, tri, dead_nbr)], &new_tris);
        }
        // legalize: every edge of every new triangle is checked against its
        // opposite triangle; the stack grows with each flip's new edges and
        // drains when no edge violates the in-circumcircle criterion
        let mut stack: Vec<[usize; 3]> = Vec::with_capacity(16);
        for (_tidx, t) in &new_tris {
            for k in 0..3 {
                stack.push([t[k], t[(k + 1) % 3], t[(k + 2) % 3]]);
            }
        }
        let mut guard = 0usize;
        while let Some([a, b, w]) = stack.pop() {
            guard += 1;
            if guard > 64 * triangles.len() + 4096 {
                break;
            }
            // the carriers of (a, b) live in endpoint a's star: the same
            // census the old full scan took, at star size instead of O(T)
            let mut found: Vec<(usize, [usize; 3])> = adjv[a]
                .iter()
                .copied()
                .filter(|t| tri_has_edge(triangles[*t], a, b))
                .map(|ti| (ti, triangles[ti]))
                .collect();
            if found.len() != 2 {
                // stale star (should not happen with maintained adjacency):
                // the old exhaustive census, unchanged
                found = triangles
                    .iter()
                    .enumerate()
                    .filter(|(_, t)| tri_has_edge(**t, a, b))
                    .map(|(ti, t)| (ti, *t))
                    .collect();
            }
            let mut it = found.into_iter();
            let (Some((ti0, tri0)), Some((ti1, tri1))) = (it.next(), it.next()) else {
                continue;
            };
            if it.next().is_some() {
                continue;
            }
            let apex = |t: [usize; 3]| t.iter().copied().find(|&q| q != a && q != b);
            // the carrier holding w as its apex owns the circle; the other
            // holds the opposite vertex x. (Stack entries may be stale
            // after earlier flips — the carriers are re-identified here.)
            let (w_tri, oti) = if apex(tri0) == Some(w) {
                (tri0, ti1)
            } else if apex(tri1) == Some(w) {
                (tri1, ti0)
            } else {
                continue;
            };
            let Some(x) = apex(triangles[oti]) else {
                continue;
            };
            // ghost-vertex rule (Shewchuk's triangle.c dummy semantics;
            // cf. Gmsh, which avoids super-triangles entirely by inserting
            // into the boundary-only initial mesh): a flip is never taken
            // on a quad touching a super-triangle vertex. Circles through
            // two real points and one super-triangle vertex are enormous
            // and would pull the interior triangles toward the ghost —
            // leaving, after the ghost discard, holes instead of the real
            // hull triangulation. Skipping keeps the mesh valid (flips are
            // the only topology changes besides splits) and the hull
            // interior tiled by real-only triangles.
            if w >= st_start || x >= st_start {
                continue;
            }
            // the directed edge (da, db) in the w-carrier's CCW order, so
            // (da, db, w) is CCW and the in-circumcircle sign reads direct
            let mut da = a;
            let mut db = b;
            for k in 0..3 {
                if w_tri[k] == a && w_tri[(k + 1) % 3] == b {
                    break;
                }
                if w_tri[k] == b && w_tri[(k + 1) % 3] == a {
                    da = b;
                    db = a;
                    break;
                }
            }
            // flip only when x lies STRICTLY inside the circumcircle of the
            // CCW triangle (da, db, w) (exact zero = cocircular tie: leave)
            if crate::predicates::incircle_sign(all[da], all[db], all[w], all[x]) <= 0 {
                continue;
            }
            // and only on a strictly convex quad. The CCW cyclic order
            // around the directed edge (da, db) is (da, x, db, w) — w sits
            // on the LEFT of da->db, x on the right — but the test accepts
            // either uniform polarity so the geometry, not the labeling,
            // decides.
            let o1 = crate::predicates::orient2d_sign(all[da], all[w], all[db]);
            let o2 = crate::predicates::orient2d_sign(all[w], all[db], all[x]);
            let o3 = crate::predicates::orient2d_sign(all[db], all[x], all[da]);
            let o4 = crate::predicates::orient2d_sign(all[x], all[da], all[w]);
            if !((o1 > 0 && o2 > 0 && o3 > 0 && o4 > 0)
                || (o1 < 0 && o2 < 0 && o3 < 0 && o4 < 0))
            {
                continue;
            }
            let n0 = tri_ccw([da, w, x], &all);
            let n1 = tri_ccw([w, db, x], &all);
            let d0 = nbr[ti0];
            let d1 = nbr[ti1];
            triangles[ti0] = n0;
            triangles[ti1] = n1;
            patch_neighbours(
                &mut nbr,
                &mut adjv,
                &[(ti0, tri0, d0), (ti1, tri1, d1)],
                &[(ti0, n0), (ti1, n1)],
            );
            for t in [n0, n1] {
                for k in 0..3 {
                    stack.push([t[k], t[(k + 1) % 3], t[(k + 2) % 3]]);
                }
            }
        }
    }

    // Drop any triangle that uses super-triangle vertices
    triangles.retain(|t| t[0] < st_start && t[1] < st_start && t[2] < st_start);
    triangles
}

/// Normalize a triangle to counter-clockwise storage.
fn tri_ccw(tri: [usize; 3], all: &[[f64; 2]]) -> [usize; 3] {
    if orient2f(all[tri[0]], all[tri[1]], all[tri[2]]) < 0.0 {
        [tri[0], tri[2], tri[1]]
    } else {
        tri
    }
}

/// Signed area cross product (twice the signed area).
fn orient2f(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

// ─── Polygon meshing ─────────────────────────────────────────────────────────

/// Generate a 2-D triangular mesh inside `polygon` with approximate target
/// edge length `mesh_size`.
///
/// The algorithm:
/// 1. Discretises the boundary edges (respecting `mesh_size`); the boundary
///    points form the constraint chain, verbatim.
/// 2. Scatters interior seed points on a staggered uniform grid (seeds
///    closer than `0.1 x mesh_size` to an accepted point are dropped).
/// 3. Runs the constrained Delaunay build ([`crate::cdt_2d::TriMesh2d`]:
///    exact-predicate Lawson initial triangulation, Sloan constraint
///    recovery, flood classification) — the live triangles cover exactly
///    the domain, concave included.
///
/// Returns a [`Mesh`] with `Triangle3` elements at z = 0.
pub fn mesh_polygon(polygon: &Polygon2D, mesh_size: f64) -> Result<Mesh, MeshError> {
    if polygon.vertices.len() < 3 {
        return Err(MeshError::Generation(
            "Polygon must have at least 3 vertices".to_string(),
        ));
    }
    if mesh_size <= 0.0 {
        return Err(MeshError::Generation(
            "mesh_size must be positive".to_string(),
        ));
    }

    let (bb_min, bb_max) = polygon.bounding_box();

    // ── Boundary points (the constraint chain, no dedup) ──────────────────
    let mut points: Vec<[f64; 2]> = Vec::new();
    let mut ring: Vec<usize> = Vec::new();
    let nv = polygon.vertices.len();
    for i in 0..nv {
        let a = polygon.vertices[i];
        let b = polygon.vertices[(i + 1) % nv];
        let len = dist2(a, b);
        let nseg = ((len / mesh_size).ceil() as usize).max(1);
        for k in 0..nseg {
            let t = k as f64 / nseg as f64;
            ring.push(points.len());
            points.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
        }
    }

    // ── Interior seed points on staggered grid ────────────────────────────
    // A seed too close to an accepted point (boundary point or earlier
    // seed) is dropped; chain vertices are never merged, so every input
    // boundary vertex survives into the mesh.
    let tol = mesh_size * 0.1;
    let mut iy = 1usize;
    let mut py = bb_min[1] + mesh_size;
    while py < bb_max[1] {
        // Offset every other row by half a cell for better triangle quality
        let offset_x = if iy % 2 == 0 { 0.0 } else { mesh_size * 0.5 };
        let mut px = bb_min[0] + offset_x + mesh_size;
        while px < bb_max[0] {
            if polygon.contains([px, py]) {
                let mut clear = true;
                for q in &points {
                    let dx = px - q[0];
                    let dy = py - q[1];
                    if dx * dx + dy * dy < tol * tol {
                        clear = false;
                        break;
                    }
                }
                if clear {
                    points.push([px, py]);
                }
            }
            px += mesh_size;
        }
        py += mesh_size * 0.866; // sqrt(3)/2 for equilateral row spacing
        iy += 1;
    }

    if points.len() < 3 {
        return Err(MeshError::Generation(
            "Too few distinct points after sampling".to_string(),
        ));
    }

    // ── Constrained Delaunay build with exact domain classification ───────
    let mut constraints: Vec<(usize, usize)> = Vec::with_capacity(ring.len());
    let n = ring.len();
    for k in 0..n {
        constraints.push((ring[k], ring[(k + 1) % n]));
    }
    let tmesh = crate::cdt_2d::TriMesh2d::build(points, &constraints, &[ring.clone()]);

    // ── Build Mesh from the live (inside-domain) triangles ─────────────────
    let mut mesh = Mesh::new();
    let mut node_id: u64 = 1;
    let mut elem_id: u64 = 1;
    let mut pt_to_node: HashMap<usize, u64> = HashMap::new();

    for t in tmesh.live() {
        let tri = tmesh.tri[t];
        let mut nids: Vec<u64> = Vec::with_capacity(3);
        for &pi in &tri {
            let nid = *pt_to_node.entry(pi).or_insert_with(|| {
                let id = node_id;
                node_id += 1;
                mesh.add_node(Node::new(id, tmesh.pts[pi][0], tmesh.pts[pi][1], 0.0));
                id
            });
            nids.push(nid);
        }
        mesh.add_element(Element::new(elem_id, ElementType::Triangle3, nids));
        elem_id += 1;
    }

    if mesh.element_count() == 0 {
        return Err(MeshError::Generation(
            "No interior triangles were generated".to_string(),
        ));
    }

    Ok(mesh)
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

pub(crate) fn tri_has_edge(tri: [usize; 3], a: usize, b: usize) -> bool {
    [[tri[0], tri[1]], [tri[1], tri[2]], [tri[2], tri[0]]]
        .iter()
        .any(|e| (e[0] == a && e[1] == b) || (e[0] == b && e[1] == a))
}

fn dist2(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    (dx * dx + dy * dy).sqrt()
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_square_delaunay() {
        // Four corners → 2 triangles (minimum)
        let pts = [[0.0f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let tris = triangulate_points(&pts);
        // trace: replicate scaled set + st to see the pre-discard mesh
        eprintln!("n={}", pts.len());
        assert!(!tris.is_empty(), "should produce at least one triangle");
        // All index values must be valid
        for t in &tris {
            assert!(t[0] < pts.len() && t[1] < pts.len() && t[2] < pts.len());
        }
    }

    #[test]
    fn mesh_unit_square() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let mesh = mesh_polygon(&poly, 0.25).expect("meshing should succeed");
        assert!(mesh.node_count() > 0);
        assert!(mesh.element_count() > 0);
        // All elements are Triangle3
        for elem in &mesh.elements {
            assert_eq!(elem.etype, ElementType::Triangle3);
            assert_eq!(elem.node_ids.len(), 3);
        }
        // All node IDs in elements must exist in the mesh
        for elem in &mesh.elements {
            for &nid in &elem.node_ids {
                assert!(mesh.nodes.contains_key(&nid));
            }
        }
    }

    #[test]
    fn mesh_l_shape() {
        // L-shaped domain (non-convex)
        let poly = Polygon2D::new(vec![
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ]);
        let mesh = mesh_polygon(&poly, 0.3).expect("L-shape meshing should succeed");
        assert!(mesh.element_count() > 0);
    }

    #[test]
    fn mesh_polygon_produces_planar_2d_mesh() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]);
        let mesh = mesh_polygon(&poly, 0.4).expect("meshing should succeed");

        assert!(mesh.node_count() > 0);
        assert!(mesh.element_count() > 0);

        for node in mesh.nodes.values() {
            assert!(node.position.z.abs() < 1e-12, "2D meshing should keep z=0");
        }

        for elem in &mesh.elements {
            assert_eq!(elem.dimension(), 2, "all elements should be 2D");
            assert_eq!(elem.etype, ElementType::Triangle3);
        }
    }

    #[test]
    fn mesh_rejects_bad_inputs() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0]]);
        assert!(mesh_polygon(&poly, 0.5).is_err(), "< 3 vertices must fail");

        let poly2 = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]]);
        assert!(
            mesh_polygon(&poly2, -1.0).is_err(),
            "negative size must fail"
        );
    }

    #[test]
    fn point_in_polygon() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        assert!(poly.contains([0.5, 0.5]));
        assert!(!poly.contains([1.5, 0.5]));
        assert!(!poly.contains([-0.1, 0.5]));
    }

    // ── P1: triangulate_points edge cases & validity ──────────────────────────

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
        // trace: replicate scaled set + st to see the pre-discard mesh
        eprintln!("n={}", pts.len());
        assert_eq!(tris.len(), 1);
        let t = tris[0];
        // All indices must be 0, 1, 2 (any order).
        let mut sorted = t;
        sorted.sort();
        assert_eq!(sorted, [0, 1, 2]);
    }

    #[test]
    fn all_output_indices_are_valid() {
        let pts: Vec<[f64; 2]> = (0..20)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 20.0;
                [angle.cos(), angle.sin()]
            })
            .collect();
        let tris = triangulate_points(&pts);
        // trace: replicate scaled set + st to see the pre-discard mesh
        eprintln!("n={}", pts.len());
        assert!(!tris.is_empty());
        for t in &tris {
            for &idx in t {
                assert!(idx < pts.len(), "index {idx} out of bounds");
            }
        }
    }

    #[test]
    fn no_degenerate_triangles_in_output() {
        let pts: Vec<[f64; 2]> = (0..15)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 15.0;
                [angle.cos(), angle.sin()]
            })
            .collect();
        let tris = triangulate_points(&pts);
        // trace: replicate scaled set + st to see the pre-discard mesh
        eprintln!("n={}", pts.len());
        for t in &tris {
            let a = pts[t[0]];
            let b = pts[t[1]];
            let c = pts[t[2]];
            // Signed area via cross product; degenerate if ~0.
            let area2 = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            assert!(area2.abs() > 1e-10, "degenerate triangle: area2={area2}");
        }
    }

    #[test]
    fn no_duplicate_triangles_in_output() {
        let pts: Vec<[f64; 2]> = (0..12)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 12.0;
                [angle.cos(), angle.sin()]
            })
            .collect();
        let tris = triangulate_points(&pts);
        // trace: replicate scaled set + st to see the pre-discard mesh
        eprintln!("n={}", pts.len());
        let mut sorted_tris: Vec<[usize; 3]> =
            tris.iter().map(|t| { let mut s = *t; s.sort(); s }).collect();
        sorted_tris.sort();
        let orig_len = sorted_tris.len();
        sorted_tris.dedup();
        assert_eq!(sorted_tris.len(), orig_len, "duplicate triangles found");
    }

    #[test]
    fn all_input_points_appear_in_output_triangulation() {
        let pts: Vec<[f64; 2]> = (0..10)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 10.0;
                [angle.cos(), angle.sin()]
            })
            .collect();
        let tris = triangulate_points(&pts);
        // trace: replicate scaled set + st to see the pre-discard mesh
        eprintln!("n={}", pts.len());
        let used: std::collections::HashSet<usize> = tris.iter().flatten().copied().collect();
        for i in 0..pts.len() {
            assert!(used.contains(&i), "point {i} not used in any triangle");
        }
    }

    // ── P1: Polygon2D helpers ─────────────────────────────────────────────────

    #[test]
    fn bounding_box_unit_square() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let (min, max) = poly.bounding_box();
        assert_eq!(min, [0.0, 0.0]);
        assert_eq!(max, [1.0, 1.0]);
    }

    #[test]
    fn bounding_box_non_axis_aligned() {
        let poly = Polygon2D::new(vec![[-2.0, 1.0], [3.0, -1.0], [1.0, 4.0]]);
        let (min, max) = poly.bounding_box();
        assert!((min[0] - (-2.0)).abs() < 1e-12);
        assert!((min[1] - (-1.0)).abs() < 1e-12);
        assert!((max[0] - 3.0).abs() < 1e-12);
        assert!((max[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn contains_on_boundary_is_consistent() {
        // Ray casting: points exactly on edges may or may not be "inside" — just
        // ensure the function doesn't panic and returns a bool.
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let _: bool = poly.contains([0.0, 0.5]); // left edge
        let _: bool = poly.contains([0.5, 0.0]); // bottom edge
    }

    // ── P1: mesh_polygon quality ──────────────────────────────────────────────

    #[test]
    fn finer_mesh_size_produces_more_elements() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let coarse = mesh_polygon(&poly, 0.5).expect("coarse should succeed");
        let fine = mesh_polygon(&poly, 0.15).expect("fine should succeed");
        assert!(
            fine.element_count() > coarse.element_count(),
            "fine={} coarse={}",
            fine.element_count(),
            coarse.element_count()
        );
    }

    #[test]
    fn mesh_polygon_all_centroids_inside_polygon() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]);
        let mesh = mesh_polygon(&poly, 0.3).expect("meshing should succeed");
        for elem in &mesh.elements {
            let pts: Vec<_> = elem
                .node_ids
                .iter()
                .map(|id| &mesh.nodes[id])
                .collect();
            let cx = pts.iter().map(|n| n.position.x).sum::<f64>() / pts.len() as f64;
            let cy = pts.iter().map(|n| n.position.y).sum::<f64>() / pts.len() as f64;
            assert!(
                poly.contains([cx, cy]),
                "element centroid ({cx},{cy}) is outside polygon"
            );
        }
    }

    #[test]
    fn mesh_polygon_zero_size_fails() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        assert!(mesh_polygon(&poly, 0.0).is_err());
    }

    /// REGRESSION GUARD (P3 quality tail card, 2026-10-13): on CONCAVE
    /// domains the old Bowyer-Watson triangulator leaked out-of-domain
    /// triangles at most mesh sizes (census on this L shape, area 3.0:
    /// 0.30 → +0.087, 0.35 → +0.213, 0.38 → +0.913, 0.40 → +0.371,
    /// 0.48 → +0.414, 0.50 → +1.350). Root cause (per-insertion probe,
    /// 2026-10-13): the super-triangle cavity scheme went non-Delaunay
    /// from the first collinear-chain insertion on (fp-level cavity ties
    /// on the huge super-triangle circles) and structurally corrupted the
    /// mesh at the notch chains — 12 overlapping triangle pairs, raw area
    /// 4.153 against a hull of 4.0. The fix is structural, not a filter
    /// patch: the exact-predicate Lawson incremental construction
    /// (splits + strictly-convex legalization flips) plus the exact
    /// triangle-in-polygon acceptance. This test pins the exact-domain
    /// area at the worst census size; the quad_recombine irregular
    /// fixture keeps its census-clean 0.45 pin.
    #[test]
    fn concave_polygon_triangulates_to_the_exact_polygon_area() {
        let poly = Polygon2D::new(vec![
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ]);
        let mesh = mesh_polygon(&poly, 0.4).expect("triangulation succeeds");
        let total: f64 = mesh
            .elements
            .iter()
            .filter_map(|e| {
                if e.node_ids.len() != 3 {
                    return None;
                }
                let p: Vec<[f64; 2]> = e
                    .node_ids
                    .iter()
                    .map(|id| {
                        let n = mesh.nodes.get(id).expect("node");
                        [n.position.x, n.position.y]
                    })
                    .collect();
                Some(
                    ((p[1][0] - p[0][0]) * (p[2][1] - p[0][1])
                        - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0]))
                        .abs()
                        / 2.0,
                )
            })
            .sum();
        assert!(
            (total - 3.0).abs() < 1e-9,
            "the L shape is 3.0 mm^2: got {total} (out-of-domain triangles)"
        );
    }
}
