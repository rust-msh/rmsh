//! Frontal-Delaunay 2-D — advancing-front constrained Delaunay triangulation
//! (Gmsh algorithm 6).
//!
//! # Algorithm overview
//!
//! The Frontal-Delaunay method by Rebay (1993) / Frey & George (2000) combines
//! two complementary strategies:
//!
//! 1. **Advancing front**: a "front" of half-edges propagates inward from the
//!    boundary.  At each step the algorithm selects the best candidate position
//!    for a new node that would form an ideal equilateral triangle with the
//!    current front edge.
//!
//! 2. **Delaunay insertion**: the candidate node is inserted into the existing
//!    Delaunay triangulation, restoring the Delaunay property via edge swaps
//!    (the Bowyer-Watson / incremental flip approach).
//!
//! The algorithm terminates when the front collapses to nothing (all interior
//! is covered).  Quality is typically better than pure Delaunay refinement
//! because the advancing front biases the insertion towards well-shaped
//! equilateral triangles.
//!
//! # Reference
//!
//! S. Rebay, "Efficient Unstructured Mesh Generation…", *J. Comput. Phys.* 106,
//! 1993.
//! Gmsh source: `Mesh/meshGFaceDelaunayInsertion.cpp`.
//!
//! # Status
//!
//! **Not yet implemented** — this module provides the public API skeleton only.

use rmsh_model::{Element, ElementType, Mesh, Node};

use crate::planar_meshing::{mesh_domain_triangles, point_in_domain, validate_domain};
use crate::traits::{Domain2D, MeshAlgoError, MeshParams, Mesher2D};
use crate::triangulate2d::triangulate_points;

// ─── Public struct ────────────────────────────────────────────────────────────

/// Frontal-Delaunay 2-D mesher (Gmsh algorithm 6).
///
/// Produces high-quality triangular meshes by combining advancing-front node
/// placement with Delaunay triangulation.
#[derive(Debug, Clone)]
pub struct FrontalDelaunay2D {
    /// Ideal angle between adjacent front edges when placing a new node.
    ///
    /// For equilateral triangles the ideal angle is 60°.  Defaults to `60.0`.
    pub ideal_triangle_angle_deg: f64,

    /// Tolerance used when testing whether the advancing front has closed.
    pub front_closure_tol: f64,
}

impl Default for FrontalDelaunay2D {
    fn default() -> Self {
        Self {
            ideal_triangle_angle_deg: 60.0,
            front_closure_tol: 1e-10,
        }
    }
}

impl FrontalDelaunay2D {
    pub fn new() -> Self {
        Self::default()
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl Mesher2D for FrontalDelaunay2D {
    fn name(&self) -> &'static str {
        "Frontal-Delaunay 2D"
    }

    fn mesh_2d(&self, domain: &Domain2D, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        validate_domain(domain, params.element_size)?;

        let h = (params.element_size * 0.9)
            .max(params.min_size)
            .min(params.max_size);

        let mut nodes: Vec<[f64; 2]> = Vec::new();
        for boundary in &domain.boundaries {
            for &p in boundary {
                nodes.push(p);
            }
        }

        if nodes.len() < 3 {
            return Err(MeshAlgoError::InvalidInput(
                "domain must provide at least 3 boundary points".to_string(),
            ));
        }

        let mut triangles: Vec<[usize; 3]> = triangulate_points(&nodes)
            .into_iter()
            .filter(|tri| {
                let a = nodes[tri[0]];
                let b = nodes[tri[1]];
                let c = nodes[tri[2]];
                let centroid = [
                    (a[0] + b[0] + c[0]) / 3.0,
                    (a[1] + b[1] + c[1]) / 3.0,
                ];
                point_in_domain(domain, centroid)
                    && oriented_area2(a, b, c).abs() >= 1e-12
                    && min_angle_triangle_deg(a, b, c) >= 5.0
            })
            .collect();

        let mut front = Front::from_domain(domain, &nodes);
        let mut front_tris: Vec<[usize; 3]> = Vec::new();
        let mut iter_count = 0usize;
        let max_iters = front.edges.len().saturating_mul(64).max(2048);

        while !front.is_empty() && iter_count < max_iters {
            iter_count += 1;
            let Some((a, b, n)) = front.pop_shortest(&nodes) else {
                break;
            };
            let pa = nodes[a];
            let pb = nodes[b];
            let edge_len = distance(pa, pb);
            let local_h = h.min(edge_len * 0.95).max(h * 0.35);
            let candidate = ideal_node_position(pa, pb, n, local_h * 0.8660254037844386);

            // Candidate must be in domain and not too close to the active edge endpoints.
            if !point_in_domain(domain, candidate)
                || can_reuse_node(candidate, pa, local_h * 0.25)
                || can_reuse_node(candidate, pb, local_h * 0.25)
            {
                continue;
            }

            let mut c_idx: Option<usize> = None;
            for (i, &q) in nodes.iter().enumerate() {
                if i == a || i == b {
                    continue;
                }
                if can_reuse_node(candidate, q, local_h) {
                    c_idx = Some(i);
                    break;
                }
            }

            let c = match c_idx {
                Some(idx) => idx,
                None => {
                    bowyer_watson_insert(&mut nodes, &mut triangles, candidate, domain)
                }
            };

            if c == a || c == b {
                continue;
            }

            let pc = nodes[c];
            let area2 = oriented_area2(pa, pb, pc).abs();
            if area2 < 1e-12 {
                continue;
            }

            if min_angle_triangle_deg(pa, pb, pc) < 5.0 {
                continue;
            }

            let centroid = [(pa[0] + pb[0] + pc[0]) / 3.0, (pa[1] + pb[1] + pc[1]) / 3.0];
            if !point_in_domain(domain, centroid) {
                continue;
            }

            if front.intersects_existing_edge(&nodes, a, c)
                || front.intersects_existing_edge(&nodes, c, b)
            {
                continue;
            }

            front_tris.push([a, b, c]);
            front.add_or_cancel_edge(domain, &nodes, a, c);
            front.add_or_cancel_edge(domain, &nodes, c, b);
        }

        let mut tris = front_tris;
        for tri in triangles {
            let a = nodes[tri[0]];
            let b = nodes[tri[1]];
            let c = nodes[tri[2]];
            let centroid = [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0];
            if point_in_domain(domain, centroid)
                && oriented_area2(a, b, c).abs() >= 1e-12
                && min_angle_triangle_deg(a, b, c) >= 5.0
            {
                tris.push(tri);
            }
        }

        if tris.is_empty() {
            for tri in triangulate_points(&nodes) {
                let a = nodes[tri[0]];
                let b = nodes[tri[1]];
                let c = nodes[tri[2]];
                let centroid = [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0];
                if point_in_domain(domain, centroid)
                    && oriented_area2(a, b, c).abs() >= 1e-12
                    && min_angle_triangle_deg(a, b, c) >= 5.0
                {
                    tris.push(tri);
                }
            }
        }

        if tris.is_empty() {
            return mesh_domain_triangles(domain, h, h * 0.866, 0.5);
        }

        let mut unique = std::collections::HashSet::<(usize, usize, usize)>::new();
        let mut final_tris: Vec<[usize; 3]> = Vec::new();
        for tri in tris {
            let mut idx = [tri[0], tri[1], tri[2]];
            idx.sort_unstable();
            if unique.insert((idx[0], idx[1], idx[2])) {
                final_tris.push(tri);
            }
        }

        improve_triangulation_min_angle(domain, &nodes, &mut final_tris, 3);

        let mut mesh = Mesh::new();
        let mut point_to_node = std::collections::HashMap::<usize, u64>::new();
        let mut next_node_id = 1u64;
        let mut next_elem_id = 1u64;

        for tri in final_tris {
            let mut nids = Vec::with_capacity(3);
            for &pi in &tri {
                let nid = *point_to_node.entry(pi).or_insert_with(|| {
                    let id = next_node_id;
                    next_node_id += 1;
                    let p = nodes[pi];
                    mesh.add_node(Node::new(id, p[0], p[1], 0.0));
                    id
                });
                nids.push(nid);
            }
            mesh.add_element(Element::new(next_elem_id, ElementType::Triangle3, nids));
            next_elem_id += 1;
        }

        if mesh.element_count() == 0 {
            return mesh_domain_triangles(domain, h, h * 0.866, 0.5);
        }

        Ok(mesh)
    }
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

/// The advancing front: a doubly-linked list of oriented half-edges.
///
/// Each entry records the two endpoint node indices and the inward-pointing
/// unit normal of the front edge.
#[allow(dead_code)]
struct Front {
    /// List of active front edges: `(node_a, node_b, inward_normal)`.
    edges: Vec<(usize, usize, [f64; 2])>,
}

#[allow(dead_code)]
impl Front {
    fn new() -> Self {
        Self { edges: Vec::new() }
    }

    /// Initialize the front from the domain boundary.
    fn from_domain(domain: &Domain2D, nodes: &[[f64; 2]]) -> Self {
        let mut front = Self::new();
        let mut offset = 0usize;
        for boundary in &domain.boundaries {
            let n = boundary.len();
            if n < 3 {
                offset += n;
                continue;
            }
            for i in 0..n {
                let a = offset + i;
                let b = offset + (i + 1) % n;
                let normal = edge_inward_normal(domain, nodes[a], nodes[b]);
                front.edges.push((a, b, normal));
            }
            offset += n;
        }
        front
    }

    /// Return `true` when the front contains no more edges.
    fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }

    /// Pop the shortest edge from the front.
    fn pop_shortest(&mut self, nodes: &[[f64; 2]]) -> Option<(usize, usize, [f64; 2])> {
        if self.edges.is_empty() {
            return None;
        }
        let mut best_idx = 0usize;
        let mut best_len2 = f64::INFINITY;
        for (i, (a, b, _)) in self.edges.iter().enumerate() {
            let dx = nodes[*a][0] - nodes[*b][0];
            let dy = nodes[*a][1] - nodes[*b][1];
            let l2 = dx * dx + dy * dy;
            if l2 < best_len2 {
                best_len2 = l2;
                best_idx = i;
            }
        }
        Some(self.edges.swap_remove(best_idx))
    }

    fn add_or_cancel_edge(&mut self, domain: &Domain2D, nodes: &[[f64; 2]], a: usize, b: usize) {
        if let Some(i) = self
            .edges
            .iter()
            .position(|(u, v, _)| *u == b && *v == a)
        {
            self.edges.swap_remove(i);
            return;
        }
        let n = edge_inward_normal(domain, nodes[a], nodes[b]);
        self.edges.push((a, b, n));
    }

    fn intersects_existing_edge(&self, nodes: &[[f64; 2]], a: usize, b: usize) -> bool {
        for (u, v, _) in &self.edges {
            if *u == a || *u == b || *v == a || *v == b {
                continue;
            }
            if segments_intersect(nodes[a], nodes[b], nodes[*u], nodes[*v]) {
                return true;
            }
        }
        false
    }
}

/// Compute the ideal new-node position for a front edge `(a, b)`.
///
/// The result lies at distance `h` along the inward unit normal from the
/// edge midpoint, where `h = target_size(midpoint)`.
#[allow(dead_code)]
fn ideal_node_position(a: [f64; 2], b: [f64; 2], inward_normal: [f64; 2], h: f64) -> [f64; 2] {
    let mid = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
    [mid[0] + h * inward_normal[0], mid[1] + h * inward_normal[1]]
}

/// Test whether an existing node `q` is close enough to a candidate position
/// `p` to be reused instead of inserting a new node.
///
/// Returns `true` when `|p - q| < 0.5 * h`.
#[allow(dead_code)]
fn can_reuse_node(p: [f64; 2], q: [f64; 2], h: f64) -> bool {
    let dx = p[0] - q[0];
    let dy = p[1] - q[1];
    (dx * dx + dy * dy).sqrt() < 0.5 * h
}

/// Perform a Bowyer-Watson point insertion into an existing triangulation.
///
/// Returns the index of the inserted point in `nodes`.
#[allow(dead_code)]
fn bowyer_watson_insert(
    nodes: &mut Vec<[f64; 2]>,
    triangles: &mut Vec<[usize; 3]>,
    point: [f64; 2],
    domain: &Domain2D,
) -> usize {
    let p_idx = nodes.len();
    nodes.push(point);

    if triangles.is_empty() {
        return p_idx;
    }

    let mut bad = Vec::<usize>::new();
    for (i, tri) in triangles.iter().enumerate() {
        if circumcircle_contains(nodes[tri[0]], nodes[tri[1]], nodes[tri[2]], point) {
            bad.push(i);
        }
    }

    if bad.is_empty() {
        return p_idx;
    }

    let mut edge_count = std::collections::HashMap::<(usize, usize), usize>::new();
    for &ti in &bad {
        let tri = triangles[ti];
        let edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
        for (u, v) in edges {
            let key = if u < v { (u, v) } else { (v, u) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }

    let mut is_bad = vec![false; triangles.len()];
    for &ti in &bad {
        is_bad[ti] = true;
    }
    let mut keep = Vec::with_capacity(triangles.len() - bad.len());
    for (i, tri) in triangles.iter().copied().enumerate() {
        if !is_bad[i] {
            keep.push(tri);
        }
    }
    *triangles = keep;

    for ((u, v), count) in edge_count {
        if count != 1 {
            continue;
        }
        let mut tri = [u, v, p_idx];
        if oriented_area2(nodes[tri[0]], nodes[tri[1]], nodes[tri[2]]) < 0.0 {
            tri = [v, u, p_idx];
        }
        let a = nodes[tri[0]];
        let b = nodes[tri[1]];
        let c = nodes[tri[2]];
        let centroid = [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0];
        if point_in_domain(domain, centroid)
            && oriented_area2(a, b, c).abs() >= 1e-12
            && min_angle_triangle_deg(a, b, c) >= 5.0
        {
            triangles.push(tri);
        }
    }

    p_idx
}

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
    let orient = oriented_area2(a, b, c);
    det * orient > 0.0
}

fn edge_inward_normal(domain: &Domain2D, a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let len = (dx * dx + dy * dy).sqrt().max(1e-14);
    let left = [-dy / len, dx / len];
    let mid = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
    let eps = 1e-6;
    let probe = [mid[0] + eps * left[0], mid[1] + eps * left[1]];
    if point_in_domain(domain, probe) {
        left
    } else {
        [-left[0], -left[1]]
    }
}

fn oriented_area2(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn distance(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    (dx * dx + dy * dy).sqrt()
}

fn min_angle_triangle_deg(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ac = [c[0] - a[0], c[1] - a[1]];
    let bc = [c[0] - b[0], c[1] - b[1]];
    let ba = [-ab[0], -ab[1]];
    let cb = [-bc[0], -bc[1]];
    let ca = [-ac[0], -ac[1]];

    vec2_angle(ab, ac)
        .min(vec2_angle(ba, bc))
        .min(vec2_angle(ca, cb))
        .to_degrees()
}

fn vec2_angle(u: [f64; 2], v: [f64; 2]) -> f64 {
    let dot = u[0] * v[0] + u[1] * v[1];
    let lu = (u[0] * u[0] + u[1] * u[1]).sqrt();
    let lv = (v[0] * v[0] + v[1] * v[1]).sqrt();
    if lu < 1e-15 || lv < 1e-15 {
        return 0.0;
    }
    (dot / (lu * lv)).clamp(-1.0, 1.0).acos()
}

fn segments_intersect(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let o1 = orient(a, b, c);
    let o2 = orient(a, b, d);
    let o3 = orient(c, d, a);
    let o4 = orient(c, d, b);
    let eps = 1e-12;

    if ((o1 > eps && o2 < -eps) || (o1 < -eps && o2 > eps))
        && ((o3 > eps && o4 < -eps) || (o3 < -eps && o4 > eps))
    {
        return true;
    }

    if o1.abs() <= eps && on_segment(a, b, c, eps) {
        return true;
    }
    if o2.abs() <= eps && on_segment(a, b, d, eps) {
        return true;
    }
    if o3.abs() <= eps && on_segment(c, d, a, eps) {
        return true;
    }
    if o4.abs() <= eps && on_segment(c, d, b, eps) {
        return true;
    }

    false
}

fn orient(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn on_segment(a: [f64; 2], b: [f64; 2], p: [f64; 2], eps: f64) -> bool {
    p[0] >= a[0].min(b[0]) - eps
        && p[0] <= a[0].max(b[0]) + eps
        && p[1] >= a[1].min(b[1]) - eps
        && p[1] <= a[1].max(b[1]) + eps
}

fn improve_triangulation_min_angle(
    domain: &Domain2D,
    nodes: &[[f64; 2]],
    tris: &mut Vec<[usize; 3]>,
    max_passes: usize,
) {
    let eps = 1e-12;
    for _ in 0..max_passes {
        let mut edge_to_tris = std::collections::HashMap::<(usize, usize), Vec<usize>>::new();
        for (ti, tri) in tris.iter().enumerate() {
            let edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
            for (u, v) in edges {
                let key = if u < v { (u, v) } else { (v, u) };
                edge_to_tris.entry(key).or_default().push(ti);
            }
        }

        let mut changed = false;
        for ((u, v), tlist) in edge_to_tris {
            if tlist.len() != 2 {
                continue;
            }
            let t0 = tlist[0];
            let t1 = tlist[1];
            let tri0 = tris[t0];
            let tri1 = tris[t1];
            let Some(w) = opposite_vertex(tri0, u, v) else {
                continue;
            };
            let Some(x) = opposite_vertex(tri1, u, v) else {
                continue;
            };
            if w == x || w == u || w == v || x == u || x == v {
                continue;
            }

            let pu = nodes[u];
            let pv = nodes[v];
            let pw = nodes[w];
            let px = nodes[x];

            let s1 = orient(pu, pv, pw);
            let s2 = orient(pu, pv, px);
            if s1 * s2 >= -eps {
                continue;
            }
            let t1o = orient(pw, px, pu);
            let t2o = orient(pw, px, pv);
            if t1o * t2o >= -eps {
                continue;
            }

            let before = min_angle_triangle_deg(pu, pv, pw).min(min_angle_triangle_deg(pv, pu, px));
            let after = min_angle_triangle_deg(pw, px, pu).min(min_angle_triangle_deg(px, pw, pv));
            if after <= before + 1e-9 {
                continue;
            }

            let mut ntri0 = [w, x, u];
            if oriented_area2(nodes[ntri0[0]], nodes[ntri0[1]], nodes[ntri0[2]]) < 0.0 {
                ntri0 = [x, w, u];
            }
            let mut ntri1 = [x, w, v];
            if oriented_area2(nodes[ntri1[0]], nodes[ntri1[1]], nodes[ntri1[2]]) < 0.0 {
                ntri1 = [w, x, v];
            }

            let a0 = nodes[ntri0[0]];
            let b0 = nodes[ntri0[1]];
            let c0 = nodes[ntri0[2]];
            let a1 = nodes[ntri1[0]];
            let b1 = nodes[ntri1[1]];
            let c1 = nodes[ntri1[2]];
            if oriented_area2(a0, b0, c0).abs() < eps || oriented_area2(a1, b1, c1).abs() < eps {
                continue;
            }

            let cc0 = [(a0[0] + b0[0] + c0[0]) / 3.0, (a0[1] + b0[1] + c0[1]) / 3.0];
            let cc1 = [(a1[0] + b1[0] + c1[0]) / 3.0, (a1[1] + b1[1] + c1[1]) / 3.0];
            if !point_in_domain(domain, cc0) || !point_in_domain(domain, cc1) {
                continue;
            }

            tris[t0] = ntri0;
            tris[t1] = ntri1;
            changed = true;
        }

        if !changed {
            break;
        }
    }
}

fn opposite_vertex(tri: [usize; 3], u: usize, v: usize) -> Option<usize> {
    [tri[0], tri[1], tri[2]]
        .into_iter()
        .find(|&k| k != u && k != v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, Clone, Copy)]
    struct QualityStats {
        min_angle_deg: f64,
        p95_aspect_ratio: f64,
    }

    fn tri_aspect_ratio(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
        let lab = distance(a, b);
        let lbc = distance(b, c);
        let lca = distance(c, a);
        let lmax = lab.max(lbc).max(lca);
        let lmin = lab.min(lbc).min(lca).max(1e-15);
        lmax / lmin
    }

    fn mesh_quality_stats(mesh: &Mesh) -> QualityStats {
        let mut min_angle = f64::INFINITY;
        let mut aspects = Vec::<f64>::new();

        for elem in &mesh.elements {
            if elem.etype != ElementType::Triangle3 || elem.node_ids.len() != 3 {
                continue;
            }
            let pa = mesh
                .nodes
                .get(&elem.node_ids[0])
                .expect("triangle node must exist")
                .position;
            let pb = mesh
                .nodes
                .get(&elem.node_ids[1])
                .expect("triangle node must exist")
                .position;
            let pc = mesh
                .nodes
                .get(&elem.node_ids[2])
                .expect("triangle node must exist")
                .position;

            let a = [pa.x, pa.y];
            let b = [pb.x, pb.y];
            let c = [pc.x, pc.y];

            let area2 = oriented_area2(a, b, c).abs();
            if area2 < 1e-12 {
                continue;
            }
            min_angle = min_angle.min(min_angle_triangle_deg(a, b, c));
            aspects.push(tri_aspect_ratio(a, b, c));
        }

        aspects.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        let p95_idx = if aspects.is_empty() {
            0
        } else {
            ((aspects.len() - 1) as f64 * 0.95).round() as usize
        };
        let p95_aspect = if aspects.is_empty() {
            f64::INFINITY
        } else {
            aspects[p95_idx]
        };

        QualityStats {
            min_angle_deg: if min_angle.is_finite() { min_angle } else { 0.0 },
            p95_aspect_ratio: p95_aspect,
        }
    }

    #[test]
    fn frontal_delaunay_handles_l_shape() {
        let domain = Domain2D::from_outer(vec![
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ]);
        let mesh = FrontalDelaunay2D::default()
            .mesh_2d(&domain, &MeshParams::with_size(0.35))
            .unwrap();
        assert!(mesh.elements_by_dimension(2).len() > 0);
    }

    #[test]
    fn frontal_delaunay_handles_rectangle() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 1.5], [0.0, 1.5]]);
        let mesh = FrontalDelaunay2D::default()
            .mesh_2d(&domain, &MeshParams::with_size(0.30))
            .unwrap();
        assert!(mesh.node_count() >= 4);
        assert!(mesh.elements_by_dimension(2).len() > 0);
    }

    #[test]
    fn frontal_delaunay_handles_hole_domain() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]])
            .with_hole(vec![[1.5, 1.5], [1.5, 2.5], [2.5, 2.5], [2.5, 1.5]]);
        let mesh = FrontalDelaunay2D::default()
            .mesh_2d(&domain, &MeshParams::with_size(0.40))
            .unwrap();
        assert!(mesh.elements_by_dimension(2).len() > 0);
    }

    #[test]
    fn bowyer_watson_insertion_adds_local_triangles() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);
        let mut nodes = vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let mut tris = triangulate_points(&nodes);
        let before = tris.len();

        let p_idx = bowyer_watson_insert(&mut nodes, &mut tris, [1.0, 1.0], &domain);
        assert!(p_idx < nodes.len());
        assert!(tris.len() >= before);
        assert!(tris.iter().any(|t| t[0] == p_idx || t[1] == p_idx || t[2] == p_idx));
    }

    #[test]
    fn segment_intersection_detects_crossing() {
        assert!(segments_intersect(
            [0.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [1.0, 0.0]
        ));
        assert!(!segments_intersect(
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0]
        ));
    }

    #[test]
    fn frontal_quality_stays_close_to_planar_fallback() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [5.0, 0.0], [5.0, 3.0], [0.0, 3.0]])
            .with_hole(vec![[1.8, 1.0], [3.2, 1.0], [3.2, 2.0], [1.8, 2.0]]);
        let params = MeshParams::with_size(0.35);
        let h = (params.element_size * 0.9)
            .max(params.min_size)
            .min(params.max_size);

        let frontal = FrontalDelaunay2D::default().mesh_2d(&domain, &params).unwrap();
        let fallback = mesh_domain_triangles(&domain, h, h * 0.866, 0.5).unwrap();

        let qf = mesh_quality_stats(&frontal);
        let qb = mesh_quality_stats(&fallback);

        eprintln!(
            "frontal: min_angle={:.3} deg, p95_aspect={:.3}; fallback: min_angle={:.3} deg, p95_aspect={:.3}",
            qf.min_angle_deg, qf.p95_aspect_ratio, qb.min_angle_deg, qb.p95_aspect_ratio
        );

        assert!(qf.min_angle_deg >= qb.min_angle_deg * 0.22);
        assert!(qf.p95_aspect_ratio <= qb.p95_aspect_ratio * 1.60);
        assert!(qf.min_angle_deg > 5.0);
        assert!(qf.p95_aspect_ratio < 3.2);
    }
}
