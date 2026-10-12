//! Constrained Delaunay triangulation in 2-D with incremental point
//! insertion — the live-connectivity mesh structure behind the Frontal-
//! Delaunay mesher.
//!
//! # Algorithm overview
//!
//! The structure mirrors the BDS/MTri3 machinery Gmsh uses for its 2-D
//! frontal-Delaunay algorithms (`Mesh/meshGFaceDelaunay.cpp`,
//! `Mesh/meshGFaceTri3.h`): triangles carry vertex indices plus per-edge
//! neighbour links, points are inserted by Bowyer-Watson cavity retriangulation
//! (`recurFindCavity` / `insertVertexB`), and the domain boundary is a set of
//! constraint edges that cavity growth and swaps may never cross.
//!
//! Build pipeline:
//! 1. Unconstrained Bowyer-Watson triangulation of the point set (the
//!    super-triangle Delaunay in [`crate::triangulate2d`]).
//! 2. Constraint recovery (Sloan's flip-based edge insertion).
//! 3. Inside/outside classification by flood fill from the convex-hull
//!    boundary, with constraint edges as walls — the exact-domain analogue of
//!    Gmsh's `colorExteriorTriangles`.
//!
//! # Reference
//!
//! S. Sloan, "A fast algorithm for constructing constrained Delaunay
//! triangulations", *Int. J. Numer. Methods Eng.* 37, 1993.
//! Gmsh source: `Mesh/meshGFaceDelaunay.cpp` (`insertVertexB`,
//! `recurFindCavity`), `Mesh/meshGFace.cpp` (`colorExteriorTriangles`).

use std::collections::{HashMap, HashSet};

use crate::triangulate2d::triangulate_points;

/// A triangle-based 2-D mesh with neighbour links and constraint edges.
///
/// Triangles are kept counter-clockwise. `neigh[t][k]` is the index of the
/// neighbour across the edge `(tri[t][k], tri[t][(k + 1) % 3])`, or a negative
/// value when the edge is on the mesh boundary. Dead triangles are deleted
/// slots retained for index stability.
#[derive(Debug, Clone, Default)]
pub(crate) struct TriMesh2d {
    /// Vertices of the triangulation (the caller's point set).
    pub pts: Vec<[f64; 2]>,
    /// Triangle corner indices, counter-clockwise.
    pub tri: Vec<[usize; 3]>,
    /// Per-edge neighbour triangle indices (negative = none).
    pub neigh: Vec<[i32; 3]>,
    /// Deleted slots (compacted away by [`TriMesh2d::live`]).
    pub dead: Vec<bool>,
    /// Constraint (domain boundary) edges as normalized `(min, max)` pairs.
    pub constraint: HashSet<(usize, usize)>,
    /// Per-vertex mesh-implied size, FROZEN at the initial mesh (Gmsh
    /// `vSizes`: `setLcs` over the initial 1D discretization; new vertices
    /// get the barycentric interpolation of their host triangle — the
    /// min-incident-edge field is never re-shrunk during insertion).
    pub vsize: Vec<f64>,
    /// Per-vertex background (size-field) size (Gmsh `vSizesBGM`); NaN when
    /// no field is in play.
    pub bgm: Vec<f64>,
}

/// Normalize an undirected edge to its `(min, max)` key.
pub(crate) fn edge_key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// Twice the signed area of `(a, b, c)` (positive = CCW).
fn orient(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// Whether the circumcircle of `(a, b, c)` contains `p` (orientation-safe).
fn in_circumcircle(a: [f64; 2], b: [f64; 2], c: [f64; 2], p: [f64; 2]) -> bool {
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
    let ori = orient(a, b, c);
    det * ori > 0.0
}

/// Euclidean distance between two points.
fn dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}

/// Proper crossing test between open segments `(p, q)` and `(a, b)`:
/// strict sign change on both sides, no touching endpoints.
fn segments_properly_cross(p: [f64; 2], q: [f64; 2], a: [f64; 2], b: [f64; 2]) -> bool {
    let o1 = orient(p, q, a);
    let o2 = orient(p, q, b);
    let o3 = orient(a, b, p);
    let o4 = orient(a, b, q);
    ((o1 > 0.0 && o2 < 0.0) || (o1 < 0.0 && o2 > 0.0))
        && ((o3 > 0.0 && o4 < 0.0) || (o3 < 0.0 && o4 > 0.0))
}

/// Closed-segment intersection: the segments share at least one point
/// (proper crossing or an endpoint lying on the other segment). The
/// classification flood uses this so a centroid jump cannot slip through a
/// boundary corner.
fn segments_intersect_closed(p: [f64; 2], q: [f64; 2], a: [f64; 2], b: [f64; 2]) -> bool {
    let o1 = orient(p, q, a);
    let o2 = orient(p, q, b);
    let o3 = orient(a, b, p);
    let o4 = orient(a, b, q);
    if ((o1 > 0.0 && o2 < 0.0) || (o1 < 0.0 && o2 > 0.0))
        && ((o3 > 0.0 && o4 < 0.0) || (o3 < 0.0 && o4 > 0.0))
    {
        return true;
    }
    let on = |s: [f64; 2], e: [f64; 2], x: [f64; 2]| {
        x[0] >= s[0].min(e[0]) && x[0] <= s[0].max(e[0])
            && x[1] >= s[1].min(e[1]) && x[1] <= s[1].max(e[1])
    };
    if o1 == 0.0 && on(p, q, a) {
        return true;
    }
    if o2 == 0.0 && on(p, q, b) {
        return true;
    }
    if o3 == 0.0 && on(a, b, p) {
        return true;
    }
    if o4 == 0.0 && on(a, b, q) {
        return true;
    }
    false
}

/// Circumcenter of `(a, b, c)`, `None` for (near-)degenerate triangles.
pub(crate) fn circumcenter(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Option<[f64; 2]> {
    let d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    if d.abs() < 1e-30 {
        return None;
    }
    let a2 = a[0] * a[0] + a[1] * a[1];
    let b2 = b[0] * b[0] + b[1] * b[1];
    let c2 = c[0] * c[0] + c[1] * c[1];
    Some([
        (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d,
        (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d,
    ])
}

impl TriMesh2d {
    /// Build the constrained triangulation of `pts` with the given constraint
    /// edges (undirected index pairs), then classify triangles inside/outside
    /// the constraint loops and drop the outside ones.
    ///
    /// `rings` are the boundary vertex rings in input order (`rings[0]` the
    /// outer boundary, the rest holes) — they carry the side information for
    /// the classification seeds.
    ///
    /// Returns the mesh with only interior triangles live (the boundary of the
    /// live set is exactly the constraint set).
    pub fn build(
        pts: Vec<[f64; 2]>,
        constraints: &[(usize, usize)],
        rings: &[Vec<usize>],
    ) -> TriMesh2d {
        let mut mesh = TriMesh2d {
            vsize: vec![f64::NAN; pts.len()],
            bgm: vec![f64::NAN; pts.len()],
            pts,
            ..TriMesh2d::default()
        };
        // 1. Unconstrained Bowyer-Watson triangulation, canonicalized CCW.
        let raw = triangulate_points(&mesh.pts);
        for t in raw {
            let t = if orient(mesh.pts[t[0]], mesh.pts[t[1]], mesh.pts[t[2]]) < 0.0 {
                [t[0], t[2], t[1]]
            } else {
                t
            };
            mesh.tri.push(t);
            mesh.neigh.push([-1; 3]);
            mesh.dead.push(false);
        }
        mesh.rebuild_neighbours();
        mesh.constraint = constraints.iter().map(|&(a, b)| edge_key(a, b)).collect();
        // 2. Constraint recovery (Sloan's flip-based insertion).
        for &(a, b) in constraints {
            mesh.recover_constraint(a, b);
        }
        // 3. Classification: flood fill from the ring seeds, constraint edges
        //    as walls (the exact-domain analogue of colorExteriorTriangles).
        let outside = mesh.classify_outside(rings);
        for ti in outside {
            mesh.dead[ti] = true;
        }
        mesh.init_vsizes();
        mesh
    }

    /// Rebuild all neighbour links from the triangle set (edge-key matching).
    pub(crate) fn rebuild_neighbours(&mut self) {
        let mut edge_tri: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for t in 0..self.tri.len() {
            for k in 0..3 {
                edge_tri
                    .entry(edge_key(self.tri[t][k], self.tri[t][(k + 1) % 3]))
                    .or_default()
                    .push(t);
            }
        }
        self.neigh = vec![[-1; 3]; self.tri.len()];
        for (&key, ts) in &edge_tri {
            if ts.len() != 2 {
                continue;
            }
            let (t0, t1) = (ts[0], ts[1]);
            for t in [t0, t1] {
                for k in 0..3 {
                    if edge_key(self.tri[t][k], self.tri[t][(k + 1) % 3]) == key {
                        self.neigh[t][k] = if t == t0 { t1 as i32 } else { t0 as i32 };
                    }
                }
            }
        }
    }

    /// The live (not deleted) triangle indices, in index order.
    pub fn live(&self) -> Vec<usize> {
        (0..self.tri.len()).filter(|&t| !self.dead[t]).collect()
    }

    /// Oriented area of a live triangle (positive, CCW).
    pub fn area2(&self, t: usize) -> f64 {
        let [a, b, c] = self.tri[t];
        orient(self.pts[a], self.pts[b], self.pts[c])
    }

    /// Circumradius of a triangle.
    pub fn circumradius(&self, t: usize) -> f64 {
        let [a, b, c] = self.tri[t];
        match circumcenter(self.pts[a], self.pts[b], self.pts[c]) {
            Some(cc) => dist(cc, self.pts[a]),
            None => f64::INFINITY,
        }
    }

    /// Local size of a triangle (Gmsh `insertVertexB`): the smaller of the
    /// mean `vSizes` (frozen mesh-implied) and the mean `vSizesBGM` (field):
    /// `LL = min(lc, lcBGM)` with `lc = ONE_THIRD * sum(sizes)`.
    pub fn lc(&self, t: usize) -> f64 {
        let [a, b, c] = self.tri[t];
        let vs = |arr: &Vec<f64>, v: usize| {
            let sv = arr[v];
            if sv.is_finite() && sv > 0.0 {
                sv
            } else {
                f64::INFINITY
            }
        };
        let mean_vs = (vs(&self.vsize, a) + vs(&self.vsize, b) + vs(&self.vsize, c)) / 3.0;
        let mean_bgm = (vs(&self.bgm, a) + vs(&self.bgm, b) + vs(&self.bgm, c)) / 3.0;
        mean_vs.min(mean_bgm)
    }

    /// Normalized triangle radius (Gmsh `MTri3::getRadius` with
    /// `radiusNorm == 2`): circumradius divided by the local size.
    pub fn radius(&self, t: usize) -> f64 {
        let l = self.lc(t);
        if !l.is_finite() || l <= 0.0 {
            return f64::INFINITY;
        }
        self.circumradius(t) / l
    }

    /// Whether `p` lies strictly inside the open segment `(a, b)`
    /// (collinear and between the endpoints).
    fn point_on_open_segment(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> bool {
        let o = orient(a, b, p);
        let ab2 = (b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2);
        if o.abs() > 1e-12 * ab2.sqrt().max(1e-30) {
            return false;
        }
        let dot = (p[0] - a[0]) * (p[0] - b[0]) + (p[1] - a[1]) * (p[1] - b[1]);
        dot < 0.0
    }

    /// Recover the constraint edge `(a, b)`: split it at any vertex lying on
    /// it, then flip intersecting edges until the edge exists (Sloan 1993).
    fn recover_constraint(&mut self, a: usize, b: usize) {
        // a vertex strictly inside the segment forces a split (collinear
        // boundary chains)
        let split = (0..self.pts.len())
            .find(|&i| i != a && i != b && Self::point_on_open_segment(self.pts[i], self.pts[a], self.pts[b]));
        if let Some(m) = split {
            self.recover_constraint(a, m);
            self.recover_constraint(m, b);
            return;
        }
        if self.has_edge(a, b) {
            return;
        }
        // collect edges properly crossing (a, b)
        let mut crossing: Vec<(usize, usize)> = Vec::new();
        for t in self.live() {
            for k in 0..3 {
                let u = self.tri[t][k];
                let v = self.tri[t][(k + 1) % 3];
                if u == a || u == b || v == a || v == b {
                    continue;
                }
                if segments_properly_cross(self.pts[a], self.pts[b], self.pts[u], self.pts[v]) {
                    let key = edge_key(u, v);
                    if !crossing.contains(&key) {
                        crossing.push(key);
                    }
                }
            }
        }
        // flip crossing edges until (a, b) appears. Each flip removes one
        // crossing for the flipped edge but may create new crossings further
        // from (a, b); the loop terminates because every flip strictly
        // increases the total area of the two triangles on the (a, b) side
        // already reached (standard Sloan argument).
        let mut guard = 0usize;
        while !self.has_edge(a, b) {
            guard += 1;
            if guard > 10000 || crossing.is_empty() {
                break;
            }
            let Some((u, v)) = crossing.pop() else {
                break;
            };
            if !self.has_edge(u, v) {
                continue;
            }
            let Some((t0, k0)) = self.find_edge(u, v) else {
                continue;
            };
            if self.constraint.contains(&edge_key(u, v)) {
                // never flip another constraint edge
                continue;
            }
            let t1 = self.neigh[t0][k0];
            if t1 < 0 {
                continue;
            }
            let t1 = t1 as usize;
            let w = self.opposite(t0, u, v);
            let x = self.opposite(t1, u, v);
            let (Some(w), Some(x)) = (w, x) else {
                continue;
            };
            // the quad (u, w, v, x) must be strictly convex for the flip to
            // preserve the triangulation
            if !Self::quad_strictly_convex(
                self.pts[u], self.pts[w], self.pts[v], self.pts[x],
            ) {
                // not flippable yet: retry later (after other flips reshape
                // the pocket)
                crossing.insert(0, (u, v));
                continue;
            }
            self.flip(t0, t1, u, v, w, x);
            // after the flip, the new edges (w, x) may cross (a, b)
            if segments_properly_cross(self.pts[a], self.pts[b], self.pts[w], self.pts[x]) {
                crossing.push(edge_key(w, x));
            }
        }
        // Recovery stall: the flip pocket cannot be reshaped (non-convex
        // quads against dense boundary chains). Gmsh's remedy
        // (remeshUnrecoveredEdges) refines the 1D mesh — here: force a
        // vertex at the midpoint of the unrecovered segment and recover the
        // two halves (recursion halves the segment, so it terminates).
        if !self.has_edge(a, b) {
            if let Some(m) = self.force_constraint_vertex(a, b) {
                self.recover_constraint(a, m);
                self.recover_constraint(m, b);
            }
        }
    }

    /// Force a vertex at the midpoint of the segment `(a, b)`: when the
    /// midpoint already lies on a mesh edge, that edge is split; otherwise
    /// the midpoint is inserted into its host triangle with a forced
    /// insertion (the insertVertexB `force` path). Returns the new index.
    fn force_constraint_vertex(&mut self, a: usize, b: usize) -> Option<usize> {
        let mid = [
            (self.pts[a][0] + self.pts[b][0]) * 0.5,
            (self.pts[a][1] + self.pts[b][1]) * 0.5,
        ];
        // on an existing live edge? split it
        for t in self.live() {
            for k in 0..3 {
                let (u, v) = (self.tri[t][k], self.tri[t][(k + 1) % 3]);
                if point_segment_distance(mid, self.pts[u], self.pts[v])
                    < 1e-9 * dist(self.pts[u], self.pts[v]).max(1e-30)
                {
                    let p = self.split_edge_stale(u, v);
                    self.rebuild_neighbours();
                    let _ = p;
                    return Some(p);
                }
            }
        }
        // insert into the host triangle (forced: no too-close gate)
        let from = self.locate(mid, a % self.tri.len().max(1))?;
        self.pts.push(mid);
        let p = self.pts.len() - 1;
        if self.vsize.len() <= p {
            self.vsize.push(f64::NAN);
            self.bgm.push(f64::NAN);
        }
        let inserted = self.insert(p, from, &HashSet::new(), f64::NAN, true);
        if !inserted {
            self.pts.pop();
            return None;
        }
        Some(p)
    }

    /// Whether the quad `a, b, c, d` (in cyclic order) is strictly convex.
    fn quad_strictly_convex(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
        let o1 = orient(a, b, c);
        let o2 = orient(b, c, d);
        let o3 = orient(c, d, a);
        let o4 = orient(d, a, b);
        (o1 > 0.0 && o2 > 0.0 && o3 > 0.0 && o4 > 0.0)
            || (o1 < 0.0 && o2 < 0.0 && o3 < 0.0 && o4 < 0.0)
    }

    /// Whether the undirected edge `(a, b)` exists in the live triangulation.
    pub(crate) fn has_edge(&self, a: usize, b: usize) -> bool {
        self.find_edge(a, b).is_some()
    }

    /// Find a live triangle carrying the edge `(a, b)`: returns
    /// `(triangle, edge index k)` with `tri[t][k] == a, tri[t][k+1] == b`
    /// in some rotation.
    pub(crate) fn find_edge(&self, a: usize, b: usize) -> Option<(usize, usize)> {
        for t in self.live() {
            for k in 0..3 {
                let (u, v) = (self.tri[t][k], self.tri[t][(k + 1) % 3]);
                if (u == a && v == b) || (u == b && v == a) {
                    return Some((t, k));
                }
            }
        }
        None
    }

    /// The vertex of triangle `t` opposite the edge `(u, v)`.
    fn opposite(&self, t: usize, u: usize, v: usize) -> Option<usize> {
        self.tri[t]
            .iter()
            .find(|&&k| k != u && k != v)
            .copied()
    }

    /// Flip the shared edge `(u, v)` of triangles `t0`/`t1` (with apices
    /// `w`/`x`) to the edge `(w, x)`. The quad must be strictly convex.
    pub(crate) fn flip(&mut self, t0: usize, t1: usize, u: usize, v: usize, w: usize, x: usize) {
        // outer neighbours before rewiring
        let nb = |t: usize, a: usize, b: usize| -> i32 {
            for k in 0..3 {
                let (p, q) = (self.tri[t][k], self.tri[t][(k + 1) % 3]);
                if (p == a && q == b) || (p == b && q == a) {
                    return self.neigh[t][k];
                }
            }
            -1
        };
        let n_uw = nb(t0, u, w);
        let n_wv = nb(t0, w, v);
        let n_vx = nb(t1, v, x);
        let n_xu = nb(t1, x, u);
        // new triangles (u, w, x) and (w, v, x), CCW by quad convexity
        let mut n0 = [u, w, x];
        if orient(self.pts[n0[0]], self.pts[n0[1]], self.pts[n0[2]]) < 0.0 {
            n0 = [u, x, w];
        }
        let mut n1 = [w, v, x];
        if orient(self.pts[n1[0]], self.pts[n1[1]], self.pts[n1[2]]) < 0.0 {
            n1 = [w, x, v];
        }
        self.tri[t0] = n0;
        self.tri[t1] = n1;
        self.neigh[t0] = [-1; 3];
        self.neigh[t1] = [-1; 3];
        // outer links (back links fixed on the neighbour side too)
        self.link_edge(t0, u, w, n_uw);
        self.link_edge(t0, x, u, n_xu);
        self.link_edge(t1, w, v, n_wv);
        self.link_edge(t1, v, x, n_vx);
        // the two internal edges face each other
        self.link_edge(t0, w, x, t1 as i32);
        self.link_edge(t1, w, x, t0 as i32);
    }

    /// Set the neighbour of the edge of triangle `t` that contains vertices
    /// `a`/`b`, and fix the corresponding back link on that neighbour.
    fn link_edge(&mut self, t: usize, a: usize, b: usize, n: i32) {
        let mut k_found = None;
        for k in 0..3 {
            let (p, q) = (self.tri[t][k], self.tri[t][(k + 1) % 3]);
            if (p == a && q == b) || (p == b && q == a) {
                k_found = Some(k);
                break;
            }
        }
        let Some(k) = k_found else {
            return;
        };
        self.neigh[t][k] = n;
        if n >= 0 {
            let nt = n as usize;
            for kk in 0..3 {
                let (p, q) = (self.tri[nt][kk], self.tri[nt][(kk + 1) % 3]);
                if (p == a && q == b) || (p == b && q == a) {
                    self.neigh[nt][kk] = t as i32;
                    return;
                }
            }
        }
    }

    /// Walk from triangle `from` to the live triangle containing point `p`
    /// (visibility walk; `from` should be near `p`).
    fn locate(&self, p: [f64; 2], from: usize) -> Option<usize> {
        let mut t = from;
        let mut steps = 0usize;
        loop {
            steps += 1;
            if steps > 4 * self.tri.len() + 16 {
                return None;
            }
            if self.dead[t] {
                // walk to any live neighbour
                let Some(nt) = self.neigh[t].iter().copied().find(|&n| n >= 0 && !self.dead[n as usize])
                else {
                    return None;
                };
                t = nt as usize;
                continue;
            }
            let mut moved = false;
            for k in 0..3 {
                let (a, b) = (self.tri[t][k], self.tri[t][(k + 1) % 3]);
                let o = orient(self.pts[a], self.pts[b], p);
                if o < -1e-12 {
                    let n = self.neigh[t][k];
                    if n >= 0 && !self.dead[n as usize] {
                        t = n as usize;
                        moved = true;
                        break;
                    }
                }
            }
            if !moved {
                return Some(t);
            }
        }
    }

    /// Insert point index `p` (already present in `pts`) into the live
    /// triangulation via Bowyer-Watson cavity retriangulation.
    ///
    /// `boundary_pt` marks vertices classified on the domain boundary: the
    /// too-close edge-distance check of Gmsh's `insertVertexB` applies only
    /// to shell edges whose BOTH endpoints are boundary vertices.
    /// `ll` is the local size scale for the rejection thresholds.
    /// Returns `true` when the insertion was accepted.
    pub(crate) fn insert(
        &mut self,
        p: usize,
        from: usize,
        boundary_pt: &HashSet<usize>,
        bnew: f64,
        force: bool,
    ) -> bool {
        let Some(host) = self.locate(self.pts[p], from) else {
            return false;
        };
        // the new vertex's mesh-implied size: barycentric interpolation of
        // the host triangle's FROZEN vSizes (Gmsh insertAPoint:
        // lc1 = (1-u-v)*vSizes[0] + u*vSizes[1] + v*vSizes[2])
        let [h0, h1, h2] = self.tri[host];
        let (p0, p1, p2) = (self.pts[h0], self.pts[h1], self.pts[h2]);
        let pp = self.pts[p];
        let area = orient(p0, p1, p2);
        let vnew = if area.abs() > 1e-30 {
            let l0 = orient(p1, p2, pp) / area;
            let l1 = orient(p2, p0, pp) / area;
            let l2 = (1.0 - l0 - l1).clamp(0.0, 1.0);
            let l0 = l0.clamp(0.0, 1.0);
            let l1 = l1.clamp(0.0, 1.0);
            let s = |v: usize| {
                let sv = self.vsize[v];
                if sv.is_finite() && sv > 0.0 {
                    sv
                } else {
                    bnew
                }
            };
            (l0 * s(h0) + l1 * s(h1) + l2 * s(h2)).max(1e-6)
        } else {
            bnew
        };
        // 1. Cavity (Gmsh recurFindCavity): the located host triangle enters
        //    the cavity UNCONDITIONALLY (a point inside an obtuse triangle
        //    can lie outside its circumcircle), and growth crosses every
        //    live neighbour whose circumcircle contains p; constraint edges
        //    block growth (recurFindCavityAniso internalEdges).
        let mut in_cavity: HashSet<usize> = std::collections::HashSet::new();
        in_cavity.insert(host);
        let mut stack = vec![host];
        while let Some(t) = stack.pop() {
            for k in 0..3 {
                if self
                    .constraint
                    .contains(&edge_key(self.tri[t][k], self.tri[t][(k + 1) % 3]))
                {
                    continue;
                }
                let n = self.neigh[t][k];
                if n < 0 || self.dead[n as usize] || in_cavity.contains(&(n as usize)) {
                    continue;
                }
                let nt = n as usize;
                let [a, b, c] = self.tri[nt];
                if in_circumcircle(self.pts[a], self.pts[b], self.pts[c], self.pts[p]) {
                    in_cavity.insert(nt);
                    stack.push(nt);
                }
            }
        }
        if in_cavity.is_empty() {
            return false;
        }
        // 2-4. Shell, cycle order, and the star-shapeness test — with the
        //      classic cavity shrink: a cavity whose fan does not close
        //      (broken cycle, Euler mismatch, or area not conserved) loses
        //      its farthest triangle and the shell is rebuilt, until the fan
        //      closes or only the host remains (the host's own fan always
        //      conserves area, so the insertion terminates). Gmsh rejects
        //      such cavities outright (insertVertexB -2/-3); the shrink is
        //      the standard Bowyer-Watson remedy for non-star-shaped cavities
        //      and replaces that rejection because fp-level cavity differences
        //      left unrefined triangles near domain corners.
        let sz = |arr: &Vec<f64>, v: usize| {
            let sv = arr[v];
            if sv.is_finite() && sv > 0.0 {
                sv
            } else {
                f64::INFINITY
            }
        };
        let pp = self.pts[p];
        let mut live_cavity: Vec<usize> = in_cavity.iter().copied().collect();
        live_cavity.sort_unstable();
        let ordered: Vec<(usize, usize)> = loop {
            let cavity_set: HashSet<usize> = live_cavity.iter().copied().collect();
            let mut shell: Vec<(usize, usize)> = Vec::new();
            for &t in &live_cavity {
                for k in 0..3 {
                    let (u, v) = (self.tri[t][k], self.tri[t][(k + 1) % 3]);
                    if self.constraint.contains(&edge_key(u, v)) {
                        shell.push((u, v));
                        continue;
                    }
                    let n = self.neigh[t][k];
                    if n < 0 || self.dead[n as usize] || !cavity_set.contains(&(n as usize)) {
                        shell.push((u, v));
                    }
                }
            }
            // Euler sanity (insertVertexB: shell.size() != cavity.size()+2)
            let mut candidate = None;
            if shell.len() == live_cavity.len() + 2 {
                shell.sort_unstable();
                // head-to-tail cycle from the smallest edge
                let mut by_tail: HashMap<usize, Vec<usize>> = HashMap::new();
                for (i, &(u, _v)) in shell.iter().enumerate() {
                    by_tail.entry(u).or_default().push(i);
                }
                let mut cur = shell[0];
                let mut ord = vec![cur];
                let mut closed = true;
                for _ in 1..shell.len() {
                    match by_tail.get(&cur.1).and_then(|v| v.first()) {
                        Some(&j) => {
                            cur = shell[j];
                            ord.push(cur);
                        }
                        None => {
                            closed = false;
                            break;
                        }
                    }
                }
                if closed {
                    // area conservation + too-close rejections (the too-close
                    // gate is lifted by `force`, Gmsh insertVertexB's forced
                    // insertion used by edge-recovery refinement)
                    let old_volume: f64 =
                        live_cavity.iter().map(|&t| self.area2(t).abs() * 0.5).sum();
                    let mut new_volume = 0.0f64;
                    let mut too_close = false;
                    for &(u, v) in &ord {
                        let (pu, pv) = (self.pts[u], self.pts[v]);
                        let a = 0.5 * orient(pu, pv, pp).abs();
                        if !force && a < 1e-25 {
                            too_close = true;
                        }
                        new_volume += a;
                        let ll = ((sz(&self.vsize, u) + sz(&self.vsize, v) + vnew) / 3.0)
                            .min((sz(&self.bgm, u) + sz(&self.bgm, v) + bnew) / 3.0);
                        let ll = if ll.is_finite() { ll } else { f64::INFINITY };
                        let d1 = dist(pp, pu);
                        let d2 = dist(pp, pv);
                        let d3 = dist(pu, pv);
                        let cosv = if d1 > 0.0 && d2 > 0.0 {
                            (d1 * d1 + d2 * d2 - d3 * d3) / (2.0 * d1 * d2)
                        } else {
                            -1.0
                        };
                        if !force && (d1 < 0.5 * ll || d2 < 0.5 * ll || cosv < -0.9999) {
                            too_close = true;
                        }
                        if !force
                            && boundary_pt.contains(&u)
                            && boundary_pt.contains(&v)
                        {
                            let d4 = point_segment_distance(pp, pu, pv);
                            if d4 < 0.4 * ll {
                                too_close = true;
                            }
                        }
                    }
                    if !too_close
                        && (old_volume - new_volume).abs() <= 1e-12 * old_volume
                    {
                        candidate = Some(ord);
                    }
                }
            }
            if let Some(ord) = candidate {
                break ord;
            }
            // too-close rejections do not heal by shrinking (the offending
            // edge stays near p): reject outright — unless `force` lifts the
            // gate
            let mut too_close_any = false;
            {
                // recompute the too-close flag for the current shell
                let cavity_set: HashSet<usize> = live_cavity.iter().copied().collect();
                let mut shell: Vec<(usize, usize)> = Vec::new();
                for &t in &live_cavity {
                    for k in 0..3 {
                        let (u, v) = (self.tri[t][k], self.tri[t][(k + 1) % 3]);
                        if self.constraint.contains(&edge_key(u, v)) {
                            shell.push((u, v));
                            continue;
                        }
                        let n = self.neigh[t][k];
                        if n < 0 || self.dead[n as usize] || !cavity_set.contains(&(n as usize)) {
                            shell.push((u, v));
                        }
                    }
                }
                for &(u, v) in &shell {
                    if force {
                        break;
                    }
                    let (pu, pv) = (self.pts[u], self.pts[v]);
                    let ll = ((sz(&self.vsize, u) + sz(&self.vsize, v) + vnew) / 3.0)
                        .min((sz(&self.bgm, u) + sz(&self.bgm, v) + bnew) / 3.0);
                    let ll = if ll.is_finite() { ll } else { f64::INFINITY };
                    let d1 = dist(pp, pu);
                    let d2 = dist(pp, pv);
                    let d3 = dist(pu, pv);
                    let cosv = if d1 > 0.0 && d2 > 0.0 {
                        (d1 * d1 + d2 * d2 - d3 * d3) / (2.0 * d1 * d2)
                    } else {
                        -1.0
                    };
                    if d1 < 0.5 * ll || d2 < 0.5 * ll || cosv < -0.9999 {
                        too_close_any = true;
                    }
                    if boundary_pt.contains(&u) && boundary_pt.contains(&v)
                        && point_segment_distance(pp, pu, pv) < 0.4 * ll
                    {
                        too_close_any = true;
                    }
                }
            }
            if too_close_any && !force {
                return false;
            }
            if live_cavity.is_empty() {
                return false;
            }
            // shrink: drop the cavity triangle whose centroid is farthest
            // from the insertion point
            let mut far = live_cavity[0];
            let mut far_d = -1.0f64;
            for &t in &live_cavity {
                let c = [
                    (self.pts[self.tri[t][0]][0] + self.pts[self.tri[t][1]][0]
                        + self.pts[self.tri[t][2]][0])
                        / 3.0,
                    (self.pts[self.tri[t][0]][1] + self.pts[self.tri[t][1]][1]
                        + self.pts[self.tri[t][2]][1])
                        / 3.0,
                ];
                let d = (c[0] - pp[0]).powi(2) + (c[1] - pp[1]).powi(2);
                if d > far_d {
                    far_d = d;
                    far = t;
                }
            }
            live_cavity.retain(|&t| t != far);
        };
        // 5. Commit: delete the (possibly shrunken) cavity, add the fan
        //    (u_i, v_i, p) in cycle order; edge 0 = shell edge, edge 1 =
        //    (v, p), edge 2 = (p, u).
        if self.vsize.len() <= p {
            let grow = p + 1 - self.vsize.len();
            self.vsize.extend(std::iter::repeat(f64::NAN).take(grow));
            self.bgm.extend(std::iter::repeat(f64::NAN).take(grow));
        }
        self.vsize[p] = vnew;
        self.bgm[p] = bnew;
        for &t in &live_cavity {
            self.dead[t] = true;
        }
        let new_first = self.tri.len();
        for &(u, v) in &ordered {
            self.tri.push([u, v, p]);
            self.neigh.push([-1; 3]);
            self.dead.push(false);
        }
        let count = ordered.len() as i32;
        for i in 0..count {
            let t = new_first + i as usize;
            // edge 1 (v, p) faces the next fan triangle; edge 2 (p, u) the
            // previous one (the cycle closes head-to-tail)
            self.neigh[t][1] = new_first as i32 + (i + 1) % count;
            self.neigh[t][2] = new_first as i32 + (i + count - 1) % count;
        }
        // edge 0 of fan triangle i faces the live triangle beyond the shell
        // edge: that neighbour still points at the deleted cavity triangle —
        // fix its back link, and record ours
        for (i, &(u, v)) in ordered.iter().enumerate() {
            let t = new_first + i;
            let key = edge_key(u, v);
            let Some(cavity_t) = live_cavity.iter().copied().find(|&ct| {
                (0..3).any(|k| edge_key(self.tri[ct][k], self.tri[ct][(k + 1) % 3]) == key)
            }) else {
                continue;
            };
            let Some(k) = (0..3).find(|&k| {
                edge_key(self.tri[cavity_t][k], self.tri[cavity_t][(k + 1) % 3]) == key
            }) else {
                continue;
            };
            let n = self.neigh[cavity_t][k];
            self.neigh[t][0] = n;
            if n >= 0 && !self.dead[n as usize] {
                let nt = n as usize;
                if let Some(kk) = (0..3).find(|&kk| {
                    edge_key(self.tri[nt][kk], self.tri[nt][(kk + 1) % 3]) == key
                }) {
                    self.neigh[nt][kk] = t as i32;
                }
            }
        }
        // 6. vSizes stay frozen (no setLcs during insertion — Gmsh keeps the
        //    initial 1D discretization sizes and interpolates for new
        //    vertices); the new vertex's sizes were recorded at commit.
        true
    }

    /// Gmsh `setLcs`: keep the shortest-incident-edge value per vertex
    /// (initialization only).
    fn shrink_vsize(&mut self, u: usize, v: usize) {
        let l = dist(self.pts[u], self.pts[v]);
        for &x in &[u, v] {
            let sv = self.vsize[x];
            if !sv.is_finite() || sv < 0.0 || sv > l {
                self.vsize[x] = l;
            }
        }
    }

    /// Initialize per-vertex shortest incident edge lengths from the live
    /// triangulation (Gmsh `setLcsInit` + `setLcs` over the initial mesh).
    fn init_vsizes(&mut self) {
        self.vsize = vec![f64::NAN; self.pts.len()];
        for t in self.live() {
            for k in 0..3 {
                self.shrink_vsize(self.tri[t][k], self.tri[t][(k + 1) % 3]);
            }
        }
    }

    /// Swap the shared edge of two live triangles when it restores the
    /// Delaunay property (in-circumcircle test) — the max-min-angle flip.
    /// Constraint edges are never swapped.
    pub(crate) fn delaunay_swap_at(&mut self, t0: usize, k: usize) -> bool {
        if self.dead[t0] || self.constraint.contains(&edge_key(
            self.tri[t0][k],
            self.tri[t0][(k + 1) % 3],
        )) {
            return false;
        }
        let n = self.neigh[t0][k];
        if n < 0 {
            return false;
        }
        let t1 = n as usize;
        if self.dead[t1] {
            return false;
        }
        let (u, v) = (self.tri[t0][k], self.tri[t0][(k + 1) % 3]);
        if self.constraint.contains(&edge_key(u, v)) {
            return false;
        }
        let (Some(w), Some(x)) = (self.opposite(t0, u, v), self.opposite(t1, u, v)) else {
            return false;
        };
        if !in_circumcircle(self.pts[u], self.pts[v], self.pts[w], self.pts[x]) {
            return false;
        }
        if !Self::quad_strictly_convex(
            self.pts[u], self.pts[w], self.pts[v], self.pts[x],
        ) {
            return false;
        }
        self.flip(t0, t1, u, v, w, x);
        true
    }

    /// One full Delaunay-restore swap pass over the live mesh. Returns the
    /// number of flips.
    pub(crate) fn delaunay_swap_pass(&mut self) -> usize {
        let mut flips = 0usize;
        let mut pass_guard = 0usize;
        loop {
            let mut changed = false;
            for t in self.live() {
                for k in 0..3 {
                    if self.dead[t] {
                        break;
                    }
                    if self.delaunay_swap_at(t, k) {
                        changed = true;
                        flips += 1;
                    }
                }
            }
            pass_guard += 1;
            if !changed || pass_guard > 64 {
                break;
            }
        }
        flips
    }

    /// Classify the outside triangles: one flood seed per boundary ring (the
    /// triangle across the ring's first carried edge, on the non-domain side),
    /// then flood fill across non-constraint edges. The outer ring seeds the
    /// exterior region; each hole ring seeds its own interior region.
    fn classify_outside(&self, rings: &[Vec<usize>]) -> Vec<usize> {
        let mut seeds: Vec<usize> = Vec::new();
        for (ri, ring) in rings.iter().enumerate() {
            let n = ring.len();
            if n < 3 {
                continue;
            }
            // ring orientation by signed area
            let mut area2 = 0.0;
            for i in 0..n {
                let a = self.pts[ring[i]];
                let b = self.pts[ring[(i + 1) % n]];
                area2 += a[0] * b[1] - b[0] * a[1];
            }
            let ccw = area2 > 0.0;
            let outer = ri == 0;
            // domain side relative to the ring direction u -> v: left for a
            // CCW outer ring or a CW hole ring, right otherwise
            let domain_left = if outer { ccw } else { !ccw };
            for i in 0..n {
                let (u, v) = (ring[i], ring[(i + 1) % n]);
                let Some((t, k)) = self.find_edge(u, v) else {
                    continue;
                };
                let (eu, ev) = (self.tri[t][k], self.tri[t][(k + 1) % 3]);
                let forward = eu == u && ev == v;
                let w = self.opposite(t, eu, ev);
                let Some(w) = w else { continue };
                let t_left_of_carried = orient(self.pts[eu], self.pts[ev], self.pts[w]) > 0.0;
                let t_left_of_ring = if forward { t_left_of_carried } else { !t_left_of_carried };
                let seeded = if t_left_of_ring == domain_left {
                    // t is the domain-side triangle: seed the one across
                    let across = self.neigh[t][k];
                    if across >= 0 && !self.dead[across as usize] {
                        seeds.push(across as usize);
                        true
                    } else {
                        false
                    }
                } else {
                    seeds.push(t);
                    true
                };
                if seeded {
                    break;
                }
            }
        }
        let mut out = self.flood_outside(&seeds);
        out.sort_unstable();
        out
    }

    /// Flood fill the outside region from `seeds`. A neighbour is reachable
    /// only when the segment between the two triangle centroids does not
    /// intersect any constraint chain segment (closed-segment test).
    ///
    /// Geometric walls instead of edge-key walls: the unconstrained
    /// Bowyer-Watson stage can leave "cut-corner" edges — non-constraint
    /// edges sharing BOTH endpoints with the boundary chain (skipping over a
    /// corner) — which no edge-key wall can block, while the centroid
    /// segment between the two sides geometrically crosses the chain.
    fn flood_outside(&self, seeds: &[usize]) -> Vec<usize> {
        let chain: Vec<([f64; 2], [f64; 2])> = self
            .constraint
            .iter()
            .map(|&(a, b)| (self.pts[a], self.pts[b]))
            .collect();
        let centroid = |t: usize| {
            [
                (self.pts[self.tri[t][0]][0] + self.pts[self.tri[t][1]][0] + self.pts[self.tri[t][2]][0]) / 3.0,
                (self.pts[self.tri[t][0]][1] + self.pts[self.tri[t][1]][1] + self.pts[self.tri[t][2]][1]) / 3.0,
            ]
        };
        let mut outside: HashSet<usize> = HashSet::new();
        let mut stack: Vec<usize> = Vec::new();
        for &s in seeds {
            if !self.dead[s] {
                outside.insert(s);
                stack.push(s);
            }
        }
        while let Some(t) = stack.pop() {
            let ct = centroid(t);
            for k in 0..3 {
                let n = self.neigh[t][k];
                if n < 0 || self.dead[n as usize] || outside.contains(&(n as usize)) {
                    continue;
                }
                let nt = n as usize;
                let cn = centroid(nt);
                // wall iff the centroid-to-centroid segment hits the chain
                let blocked = chain
                    .iter()
                    .any(|&(a, b)| segments_intersect_closed(ct, cn, a, b));
                if !blocked {
                    outside.insert(nt);
                    stack.push(nt);
                }
            }
        }
        outside.into_iter().collect()
    }

    /// Index of the live triangle with the largest normalized radius (ties
    /// broken by smaller index via the ascending scan), with the Gmsh
    /// `isActive` condition on at least one of its edges, skipping the
    /// already-processed triangles.
    pub(crate) fn worst_active(&self, limit: f64, skip: &HashSet<usize>) -> Option<usize> {
        let mut best: Option<usize> = None;
        let mut best_r = limit;
        for t in 0..self.tri.len() {
            if self.dead[t] || skip.contains(&t) {
                continue;
            }
            let r = self.radius(t);
            if r <= best_r {
                continue;
            }
            if !self.is_active(t, limit) {
                continue;
            }
            best = Some(t);
            best_r = r;
        }
        best
    }

    /// Gmsh `isActive`: the triangle has an edge whose neighbour is missing
    /// or already finished (radius below the limit). Deleted neighbours count
    /// as missing (Gmsh's initial mesh has no triangles where ours carries
    /// classified-out dead slots).
    pub(crate) fn is_active(&self, t: usize, limit: f64) -> bool {
        if self.dead[t] {
            return false;
        }
        for k in 0..3 {
            let n = self.neigh[t][k];
            if n < 0 || self.dead[n as usize] {
                return true;
            }
            let nr = self.radius(n as usize);
            if nr < limit && nr > 0.0 {
                return true;
            }
        }
        false
    }

    /// The first active edge index of a triangle (Gmsh `active_edge`).
    pub(crate) fn active_edge(&self, t: usize, limit: f64) -> Option<usize> {
        for k in 0..3 {
            let n = self.neigh[t][k];
            if n < 0 || self.dead[n as usize] {
                return Some(k);
            }
            let nr = self.radius(n as usize);
            if nr < limit && nr > 0.0 {
                return Some(k);
            }
        }
        None
    }

    /// Split the edge `(u, v)` at its midpoint: the (one or two) live
    /// triangles carrying the edge are replaced by two each, and a constraint
    /// edge is replaced by its two halves. Neighbour links are left stale —
    /// callers batch splits and call [`TriMesh2d::rebuild_neighbours`] once
    /// afterwards. Returns the midpoint index.
    pub(crate) fn split_edge_stale(&mut self, u: usize, v: usize) -> usize {
        let mid = [
            (self.pts[u][0] + self.pts[v][0]) * 0.5,
            (self.pts[u][1] + self.pts[v][1]) * 0.5,
        ];
        self.pts.push(mid);
        let p = self.pts.len() - 1;
        // interpolate the frozen sizes (barycentric on a segment = the mean)
        let interp = |arr: &Vec<f64>, a: usize, b: usize| -> f64 {
            let (sa, sb) = (arr[a], arr[b]);
            match (sa.is_finite() && sa > 0.0, sb.is_finite() && sb > 0.0) {
                (true, true) => 0.5 * (sa + sb),
                (true, false) => sa,
                (false, true) => sb,
                (false, false) => f64::NAN,
            }
        };
        let vs_p = interp(&self.vsize, u, v);
        let bg_p = interp(&self.bgm, u, v);
        self.vsize.push(vs_p);
        self.bgm.push(bg_p);
        // the carriers of (u, v) among the live triangles
        let key = edge_key(u, v);
        let mut carriers: Vec<(usize, usize)> = Vec::new(); // (triangle, apex)
        for t in self.live() {
            if (0..3).any(|k| edge_key(self.tri[t][k], self.tri[t][(k + 1) % 3]) == key) {
                let apex = self
                    .tri[t]
                    .iter()
                    .copied()
                    .find(|&x| x != u && x != v)
                    .expect("edge endpoints are triangle vertices");
                carriers.push((t, apex));
            }
        }
        // kill the carriers, add the split pairs
        for &(carrier, apex) in &carriers {
            self.dead[carrier] = true;
            let mut n0 = [u, p, apex];
            if orient(self.pts[n0[0]], self.pts[n0[1]], self.pts[n0[2]]) < 0.0 {
                n0 = [u, apex, p];
            }
            self.tri.push(n0);
            self.neigh.push([-1; 3]);
            self.dead.push(false);
            let mut n1 = [p, v, apex];
            if orient(self.pts[n1[0]], self.pts[n1[1]], self.pts[n1[2]]) < 0.0 {
                n1 = [p, apex, v];
            }
            self.tri.push(n1);
            self.neigh.push([-1; 3]);
            self.dead.push(false);
        }
        // constraints follow the split
        if self.constraint.remove(&key) {
            self.constraint.insert(edge_key(u, p));
            self.constraint.insert(edge_key(p, v));
        }
        p
    }

    /// The live neighbour vertices of a vertex (unique, sorted).
    #[cfg(test)]
    pub(crate) fn vertex_neighbours(&self, v: usize) -> Vec<usize> {
        use std::collections::BTreeSet;
        let mut s: BTreeSet<usize> = BTreeSet::new();
        for t in self.live() {
            if self.tri[t].contains(&v) {
                for k in 0..3 {
                    if self.tri[t][k] != v {
                        s.insert(self.tri[t][k]);
                    }
                }
            }
        }
        s.into_iter().collect()
    }

    /// The live triangles incident to a vertex.
    #[cfg(test)]
    pub(crate) fn vertex_tris(&self, v: usize) -> Vec<usize> {
        self.live()
            .into_iter()
            .filter(|&t| self.tri[t].contains(&v))
            .collect()
    }

    /// The distinct vertices referenced by live triangles.
    pub(crate) fn live_vertices(&self) -> Vec<usize> {
        let mut s: Vec<usize> = self
            .live()
            .iter()
            .flat_map(|&t| self.tri[t])
            .collect();
        s.sort_unstable();
        s.dedup();
        s
    }
}

/// Distance from point `p` to segment `(a, b)`.
pub(crate) fn point_segment_distance(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ap = [p[0] - a[0], p[1] - a[1]];
    let ab2 = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if ab2 <= 0.0 {
        0.0
    } else {
        ((ap[0] * ab[0] + ap[1] * ab[1]) / ab2).clamp(0.0, 1.0)
    };
    let dx = p[0] - (a[0] + t * ab[0]);
    let dy = p[1] - (a[1] + t * ab[1]);
    (dx * dx + dy * dy).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring_edges(ring: &[usize]) -> Vec<(usize, usize)> {
        (0..ring.len())
            .map(|i| (ring[i], ring[(i + 1) % ring.len()]))
            .collect()
    }

    fn ring_area2(mesh: &TriMesh2d, ring: &[usize]) -> f64 {
        let n = ring.len();
        (0..n)
            .map(|i| {
                let a = mesh.pts[ring[i]];
                let b = mesh.pts[ring[(i + 1) % n]];
                a[0] * b[1] - b[0] * a[1]
            })
            .sum::<f64>()
            / 2.0
    }

    #[test]
    fn concave_l_shape_triangulates_to_the_exact_area() {
        let ring = vec![
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ];
        let rings = vec![(0..ring.len()).collect::<Vec<_>>()];
        let mesh = TriMesh2d::build(ring.clone(), &ring_edges(&rings[0]), &rings);
        // every constraint edge survives recovery
        for &(a, b) in ring_edges(&rings[0]).iter() {
            assert!(mesh.has_edge(a, b), "constraint ({a},{b}) missing");
        }
        // exact area: sum of live triangle areas matches the ring area
        let total: f64 = mesh.live().iter().map(|&t| mesh.area2(t).abs() / 2.0).sum();
        assert!(
            (total - 3.0).abs() < 1e-9,
            "L shape is 3.0: got {total} ({:?} tris)",
            mesh.live().len()
        );
        // orientation: ring CCW, live triangles CCW
        assert!(ring_area2(&mesh, &rings[0]) > 0.0);
        for &t in mesh.live().iter() {
            assert!(mesh.area2(t) > 0.0, "triangle {t} not CCW");
        }
    }

    #[test]
    fn plate_with_hole_classifies_the_hole_outside() {
        let outer = vec![[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]];
        let hole = vec![[1.5, 1.5], [1.5, 2.5], [2.5, 2.5], [2.5, 1.5]];
        let mut pts = outer.clone();
        let hole_start = pts.len();
        pts.extend(hole.iter().copied());
        let outer_ring: Vec<usize> = (0..outer.len()).collect();
        let hole_ring: Vec<usize> = (hole_start..pts.len()).collect();
        let mut constraints = ring_edges(&outer_ring);
        constraints.extend(ring_edges(&hole_ring));
        let rings = vec![outer_ring, hole_ring];
        let mesh = TriMesh2d::build(pts, &constraints, &rings);
        let total: f64 = mesh.live().iter().map(|&t| mesh.area2(t).abs() / 2.0).sum();
        assert!(
            (total - (16.0 - 1.0)).abs() < 1e-9,
            "area should be 15.0: got {total}"
        );
        // no live triangle crosses the hole
        for &t in mesh.live().iter() {
            let c = [
                (mesh.pts[mesh.tri[t][0]][0] + mesh.pts[mesh.tri[t][1]][0] + mesh.pts[mesh.tri[t][2]][0]) / 3.0,
                (mesh.pts[mesh.tri[t][0]][1] + mesh.pts[mesh.tri[t][1]][1] + mesh.pts[mesh.tri[t][2]][1]) / 3.0,
            ];
            assert!(
                !(c[0] > 1.5 && c[0] < 2.5 && c[1] > 1.5 && c[1] < 2.5),
                "triangle {t} centroid inside the hole"
            );
        }
    }

    #[test]
    fn insert_adds_a_point_with_area_conservation() {
        let outer = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let rings = vec![(0..outer.len()).collect::<Vec<_>>()];
        let mut pts = outer.clone();
        let constraints = ring_edges(&rings[0]);
        let mut mesh = TriMesh2d::build(pts.clone(), &constraints, &rings);
        let before: f64 = mesh.live().iter().map(|&t| mesh.area2(t).abs() / 2.0).sum();
        pts.push([1.0, 1.0]);
        mesh.pts = pts.clone();
        let center = pts.len() - 1;
        assert!(mesh.insert(center, 0, &rings[0].iter().copied().collect(), 1.0, false));
        let after: f64 = mesh.live().iter().map(|&t| mesh.area2(t).abs() / 2.0).sum();
        assert!(
            (before - after).abs() < 1e-9,
            "area changed on insertion: {before} -> {after}"
        );
        assert_eq!(mesh.live().len(), 4);
        for &t in mesh.live().iter() {
            assert!(mesh.area2(t) > 0.0);
        }
    }

    #[test]
    fn manifold_neighbour_links_are_symmetric() {
        let outer = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let rings = vec![(0..outer.len()).collect::<Vec<_>>()];
        let mesh = TriMesh2d::build(outer, &ring_edges(&rings[0]), &rings);
        for t in mesh.live() {
            for k in 0..3 {
                let n = mesh.neigh[t][k];
                if n < 0 {
                    continue;
                }
                let nt = n as usize;
                let key = edge_key(mesh.tri[t][k], mesh.tri[t][(k + 1) % 3]);
                let back = (0..3).any(|kk| {
                    edge_key(mesh.tri[nt][kk], mesh.tri[nt][(kk + 1) % 3]) == key
                });
                assert!(back, "asymmetric neighbour link {t}.{k} -> {nt}");
            }
        }
    }
}

