//! Watertight pre-check and repair hints for triangle surface meshes.
//!
//! This is the P4 roadmap card "watertight pre-check + repair hints": the
//! census/audit a tetrahedralization front-end needs before it accepts a
//! surface mesh, reported as structured defect records plus one actionable
//! hint string per defect. It is deliberately report-only: the input
//! [`Mesh`] is never mutated; repair actions belong to the consumer.
//!
//! Scope: the audit inspects the triangle soup formed by the mesh's 2D
//! elements that carry exactly three node ids (Triangle3). A volume mesh
//! must have its boundary extracted into surface elements before auditing.
//! Triangles referencing node ids that are missing from the mesh are skipped
//! (they cannot be positioned in space).
//!
//! Detection concepts follow mesh-repair tooling (AdMesh/JMesh style) and
//! Gmsh's mesh classification needs (GModel::classifyFaces walks the same
//! edge-incidence structure to build boundary loops):
//! - boundary edges (1 incidence) grouped into loops = holes;
//! - non-manifold edges (>= 3 incidence) listed with their incidence;
//! - orientation inconsistencies = an undirected edge whose directed
//!   traversals do not balance (a consistent closed mesh traverses every
//!   edge once in each direction);
//! - duplicate triangles (same node set);
//! - degenerate triangles (repeated node id or area below 1e-12);
//! - self-intersections = triangle pairs sharing NO node whose geometry
//!   overlaps (Moeller's triangle-triangle test; see `triangles_intersect`).
//!
//! Scaling note: the self-intersection sweep is O(n^2) over triangle pairs
//! with an axis-aligned bounding-box pre-filter. That is fine at the
//! precheck/fixture scale (up to a few thousand triangles). For larger
//! meshes the pair enumeration can be replaced by a uniform spatial hash or
//! a BVH over the same per-triangle AABBs without touching the test itself.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use rmsh_model::Mesh;

/// Area below this threshold counts as a degenerate (zero-area) triangle.
const DEGENERATE_AREA_EPS: f64 = 1e-12;
/// Planar signed-distance magnitude below this counts as "on the plane" in
/// the triangle-triangle intersection test.
const PLANE_DIST_EPS: f64 = 1e-12;

/// Why a triangle was flagged degenerate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegenerateReason {
    /// A node id appears more than once in the triangle connectivity.
    RepeatedNode,
    /// The triangle's 3D area is below [`DEGENERATE_AREA_EPS`].
    ZeroArea,
}

/// One boundary loop of the triangle soup: a hole to cap.
#[derive(Debug, Clone)]
pub struct BoundaryLoop {
    /// Whether the walk closed back onto its start node. An open walk means
    /// the boundary is not a manifold curve (odd incidence at some node).
    pub closed: bool,
    /// Boundary edges directed along the traversal.
    pub edges: Vec<(u64, u64)>,
    /// Nodes along the traversal (when `closed`, the cycle closes back to
    /// `node_ids[0]` without repeating it).
    pub node_ids: Vec<u64>,
    /// Sum of the 3D lengths of the loop edges.
    pub total_length: f64,
}

impl BoundaryLoop {
    /// Number of boundary edges on this loop.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Actionable repair hint for this loop.
    pub fn hint(&self) -> String {
        if self.closed {
            format!(
                "hole: cap the boundary loop of {} edges {:?} (total length {:.6})",
                self.edge_count(),
                self.edges,
                self.total_length
            )
        } else {
            format!(
                "open boundary chain: close or cap the {} edges {:?} (total length {:.6})",
                self.edge_count(),
                self.edges,
                self.total_length
            )
        }
    }
}

/// An edge shared by 3 or more triangles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NonManifoldEdge {
    /// The edge as an ordered (min, max) node id pair.
    pub nodes: (u64, u64),
    /// Number of incident triangles.
    pub incidence: usize,
}

impl NonManifoldEdge {
    /// Actionable repair hint for this edge.
    pub fn hint(&self) -> String {
        format!(
            "non-manifold edge ({}, {}) has {} incident triangles: merge or explode the incidence",
            self.nodes.0, self.nodes.1, self.incidence
        )
    }
}

/// An undirected edge whose directed traversals do not balance.
///
/// A consistently oriented closed mesh sends every undirected edge once in
/// each direction; `forward` counts `(min, max)` traversals and `backward`
/// counts `(max, min)` traversals. Edges with a single traversal are
/// boundary (reported as holes), not orientation mismatches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OrientationMismatch {
    /// The edge as an ordered (min, max) node id pair.
    pub nodes: (u64, u64),
    /// Traversals directed `(min, max)`.
    pub forward: usize,
    /// Traversals directed `(max, min)`.
    pub backward: usize,
}

impl OrientationMismatch {
    /// Actionable repair hint for this mismatch.
    pub fn hint(&self) -> String {
        format!(
            "edge ({}, {}) traversed {} times forward / {} times backward: flip one triangle using this edge",
            self.nodes.0, self.nodes.1, self.forward, self.backward
        )
    }
}

/// A group of two or more triangles built from the same node set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DuplicateTriangle {
    /// Element ids of every triangle in the group (encounter order).
    pub element_ids: Vec<u64>,
    /// The shared node set, sorted.
    pub nodes: [u64; 3],
}

impl DuplicateTriangle {
    /// Actionable repair hint for this group.
    pub fn hint(&self) -> String {
        format!(
            "remove duplicate triangle element(s) {:?} with nodes {:?} (keep exactly one)",
            self.element_ids, self.nodes
        )
    }
}

/// A degenerate triangle flagged by the audit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DegenerateTriangle {
    /// The offending element's id.
    pub element_id: u64,
    /// Why it was flagged.
    pub reason: DegenerateReason,
    /// The triangle's 3D area.
    pub area: f64,
}

impl DegenerateTriangle {
    /// Actionable repair hint for this triangle.
    pub fn hint(&self) -> String {
        let why = match self.reason {
            DegenerateReason::RepeatedNode => "repeated node id",
            DegenerateReason::ZeroArea => "zero area",
        };
        format!(
            "remove or repair degenerate triangle element {} ({}: area {:.3e})",
            self.element_id, why, self.area
        )
    }
}

/// A pair of triangles that geometrically intersect while sharing no node.
///
/// Pairs that share a node are exempt: adjacent triangles legitimately touch
/// along edges and at nodes; the defect class is the non-adjacent overlap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SelfIntersection {
    /// Element id of the first triangle.
    pub element_a: u64,
    /// Element id of the second triangle.
    pub element_b: u64,
}

impl SelfIntersection {
    /// Actionable repair hint for this pair.
    pub fn hint(&self) -> String {
        format!(
            "triangles element {} and element {} intersect while sharing no node: split or remesh the overlap region",
            self.element_a, self.element_b
        )
    }
}

/// The full watertight census of a triangle mesh.
#[derive(Debug, Clone, Default)]
pub struct WatertightReport {
    /// Boundary loops (holes), in deterministic order.
    pub boundary_loops: Vec<BoundaryLoop>,
    /// Over-shared edges, ordered by node pair.
    pub non_manifold_edges: Vec<NonManifoldEdge>,
    /// Unbalanced directed edges, ordered by node pair.
    pub orientation_mismatches: Vec<OrientationMismatch>,
    /// Same-node-set triangle groups.
    pub duplicate_triangles: Vec<DuplicateTriangle>,
    /// Degenerate triangles, in element order.
    pub degenerate_triangles: Vec<DegenerateTriangle>,
    /// Non-adjacent intersecting triangle pairs.
    pub self_intersections: Vec<SelfIntersection>,
    /// True only when every defect class above is empty.
    pub is_watertight: bool,
    /// One actionable hint per detected defect, in class order.
    pub hints: Vec<String>,
}

/// One audited triangle: connectivity plus 3D positions.
#[derive(Debug, Clone)]
struct Tri {
    element_id: u64,
    nodes: [u64; 3],
    pos: [[f64; 3]; 3],
}

/// Audit a triangle mesh for watertightness defects (see module docs).
///
/// The signature is stable: it never mutates the mesh and never fails; every
/// finding is structured data in the returned [`WatertightReport`].
pub fn audit_watertight(mesh: &Mesh) -> WatertightReport {
    let tris = collect_triangles(mesh);

    let mut non_manifold_edges = Vec::new();
    let mut orientation_mismatches = Vec::new();
    let mut duplicate_triangles = Vec::new();
    let mut degenerate_triangles = Vec::new();
    let mut self_intersections = Vec::new();

    // Degenerate triangles: repeated node ids, else area below threshold.
    for t in &tris {
        let area = triangle_area(&t.pos);
        if t.nodes[0] == t.nodes[1] || t.nodes[1] == t.nodes[2] || t.nodes[0] == t.nodes[2] {
            degenerate_triangles.push(DegenerateTriangle {
                element_id: t.element_id,
                reason: DegenerateReason::RepeatedNode,
                area,
            });
        } else if area < DEGENERATE_AREA_EPS {
            degenerate_triangles.push(DegenerateTriangle {
                element_id: t.element_id,
                reason: DegenerateReason::ZeroArea,
                area,
            });
        }
    }

    // Duplicates: group triangles by their sorted node set.
    let mut by_nodes: BTreeMap<[u64; 3], Vec<u64>> = BTreeMap::new();
    for t in &tris {
        let mut key = t.nodes;
        key.sort_unstable();
        by_nodes.entry(key).or_default().push(t.element_id);
    }
    for (nodes, ids) in &by_nodes {
        if ids.len() >= 2 {
            duplicate_triangles.push(DuplicateTriangle {
                element_ids: ids.clone(),
                nodes: *nodes,
            });
        }
    }

    // Edge incidence with directed traversal counts: key -> (forward, backward)
    // where forward counts (min, max) traversals and backward (max, min).
    let mut edge_dirs: HashMap<(u64, u64), (usize, usize)> = HashMap::new();
    for t in &tris {
        for i in 0..3 {
            let a = t.nodes[i];
            let b = t.nodes[(i + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            let e = edge_dirs.entry(key).or_insert((0, 0));
            if a < b {
                e.0 += 1;
            } else {
                e.1 += 1;
            }
        }
    }

    // Boundary edges (incidence 1), sorted for deterministic loop tracing.
    let mut boundary_edges: Vec<(u64, u64)> = edge_dirs
        .iter()
        .filter(|(_, c)| c.0 + c.1 == 1)
        .map(|(k, _)| *k)
        .collect();
    boundary_edges.sort_unstable();
    let boundary_loops = trace_boundary_loops(&boundary_edges, mesh);

    // Non-manifold edges (incidence >= 3).
    for (&nodes, &(forward, backward)) in &edge_dirs {
        let incidence = forward + backward;
        if incidence >= 3 {
            non_manifold_edges.push(NonManifoldEdge { nodes, incidence });
        }
    }
    non_manifold_edges.sort_by_key(|e| e.nodes);

    // Orientation mismatches: incidence >= 2 with unbalanced directions.
    // (Incidence-1 edges are boundary, already reported as holes.)
    for (&nodes, &(forward, backward)) in &edge_dirs {
        if forward + backward >= 2 && forward != backward {
            orientation_mismatches.push(OrientationMismatch {
                nodes,
                forward,
                backward,
            });
        }
    }
    orientation_mismatches.sort_by_key(|m| m.nodes);

    // Self-intersections: O(n^2) sweep over triangle pairs (see the module
    // scaling note for the spatial-hash / BVH upgrade path). Pairs sharing a
    // node are exempt: adjacent triangles legitimately touch along edges and
    // at nodes; the defect class is the non-adjacent overlap.
    let aabbs: Vec<([f64; 3], [f64; 3])> = tris.iter().map(triangle_aabb).collect();
    for i in 0..tris.len() {
        for j in (i + 1)..tris.len() {
            if tris[i].nodes.iter().any(|n| tris[j].nodes.contains(n)) {
                continue;
            }
            if !aabb_overlap(&aabbs[i], &aabbs[j]) {
                continue;
            }
            if triangles_intersect(&tris[i], &tris[j]) {
                self_intersections.push(SelfIntersection {
                    element_a: tris[i].element_id,
                    element_b: tris[j].element_id,
                });
            }
        }
    }

    let is_watertight = boundary_loops.is_empty()
        && non_manifold_edges.is_empty()
        && orientation_mismatches.is_empty()
        && duplicate_triangles.is_empty()
        && degenerate_triangles.is_empty()
        && self_intersections.is_empty();

    // One actionable hint per defect, in class order.
    let mut hints = Vec::new();
    hints.extend(boundary_loops.iter().map(|d| d.hint()));
    hints.extend(non_manifold_edges.iter().map(|d| d.hint()));
    hints.extend(orientation_mismatches.iter().map(|d| d.hint()));
    hints.extend(duplicate_triangles.iter().map(|d| d.hint()));
    hints.extend(degenerate_triangles.iter().map(|d| d.hint()));
    hints.extend(self_intersections.iter().map(|d| d.hint()));

    WatertightReport {
        boundary_loops,
        non_manifold_edges,
        orientation_mismatches,
        duplicate_triangles,
        degenerate_triangles,
        self_intersections,
        is_watertight,
        hints,
    }
}

/// Collect the auditable triangles: 2D elements with exactly three node ids
/// whose nodes all exist in the mesh (unpositionable triangles are skipped).
fn collect_triangles(mesh: &Mesh) -> Vec<Tri> {
    let mut tris = Vec::new();
    for elem in &mesh.elements {
        if elem.dimension() != 2 || elem.node_ids.len() != 3 {
            continue;
        }
        let mut pos = [[0.0f64; 3]; 3];
        let mut positioned = true;
        for (k, &nid) in elem.node_ids.iter().enumerate() {
            match mesh.nodes.get(&nid) {
                Some(n) => pos[k] = [n.position.x, n.position.y, n.position.z],
                None => {
                    positioned = false;
                    break;
                }
            }
        }
        if !positioned {
            continue;
        }
        tris.push(Tri {
            element_id: elem.id,
            nodes: [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]],
            pos,
        });
    }
    tris
}

/// Trace boundary edges into loops: consume edges node-by-node starting from
/// the smallest unused edge until the walk closes onto its start or dead-ends
/// (odd boundary incidence => an open chain, flagged via `closed: false`).
fn trace_boundary_loops(boundary_edges: &[(u64, u64)], mesh: &Mesh) -> Vec<BoundaryLoop> {
    let mut adj: HashMap<u64, BTreeSet<(u64, u64)>> = HashMap::new();
    for &key in boundary_edges {
        adj.entry(key.0).or_default().insert(key);
        adj.entry(key.1).or_default().insert(key);
    }

    let mut consumed: HashSet<(u64, u64)> = HashSet::new();
    let mut loops = Vec::new();
    for &start in boundary_edges {
        if consumed.contains(&start) {
            continue;
        }
        consumed.insert(start);
        let (s, first) = start;
        let mut edges = vec![start];
        let mut nodes = vec![s, first];
        let mut closed = false;
        let mut node = first;
        loop {
            let mut next_key = None;
            if let Some(set) = adj.get(&node) {
                for &key in set {
                    if !consumed.contains(&key) {
                        next_key = Some(key);
                        break;
                    }
                }
            }
            let Some(key) = next_key else {
                break; // dead end: open boundary chain
            };
            consumed.insert(key);
            let other = if key.0 == node { key.1 } else { key.0 };
            edges.push((node, other));
            if other == s {
                closed = true;
                break;
            }
            nodes.push(other);
            node = other;
        }
        let total_length: f64 = edges
            .iter()
            .map(|(a, b)| match (mesh.nodes.get(a), mesh.nodes.get(b)) {
                (Some(pa), Some(pb)) => dist3(
                    &[pa.position.x, pa.position.y, pa.position.z],
                    &[pb.position.x, pb.position.y, pb.position.z],
                ),
                _ => 0.0,
            })
            .sum();
        loops.push(BoundaryLoop {
            closed,
            edges,
            node_ids: nodes,
            total_length,
        });
    }
    loops
}

/// Triangle-triangle intersection test after Moeller (1997), "A fast
/// triangle-triangle intersection test", Journal of Graphics Tools 2(2):
/// reject on coherent plane-sidedness, then compare the two crossing
/// intervals along the line of intersection; coplanar pairs fall back to a
/// 2D projection test. Reimplemented from the published method (no code
/// copied); shared-node pairs never reach this function.
fn triangles_intersect(a: &Tri, b: &Tri) -> bool {
    let Some((n2, d2)) = plane_of(b) else {
        return false; // degenerate b: no plane, flagged separately
    };
    let da: [f64; 3] = a.pos.map(|p| dot3(n2, p) + d2);
    if da.iter().all(|d| d.abs() <= PLANE_DIST_EPS) {
        return coplanar_intersect(a, b, &n2);
    }
    if da.iter().all(|d| *d > PLANE_DIST_EPS) || da.iter().all(|d| *d < -PLANE_DIST_EPS) {
        return false;
    }

    let Some((n1, d1)) = plane_of(a) else {
        return false;
    };
    let db: [f64; 3] = b.pos.map(|p| dot3(n1, p) + d1);
    if db.iter().all(|d| d.abs() <= PLANE_DIST_EPS) {
        return coplanar_intersect(a, b, &n1);
    }
    if db.iter().all(|d| *d > PLANE_DIST_EPS) || db.iter().all(|d| *d < -PLANE_DIST_EPS) {
        return false;
    }

    // Both triangles straddle the other's plane: compare crossing intervals
    // along the intersection line N1 x N2, projected on its dominant axis.
    let dir = cross3(n1, n2);
    let axis = dominant_axis(&dir);
    let (t_lo, t_hi) = crossing_interval(a, &da, axis);
    let (s_lo, s_hi) = crossing_interval(b, &db, axis);
    t_lo <= s_hi && s_lo <= t_hi
}

/// The plane of a triangle as (normal, d) with dot(n, x) + d = 0; `None` for
/// a degenerate (zero-normal) triangle.
fn plane_of(t: &Tri) -> Option<([f64; 3], f64)> {
    let n = cross3(sub3(t.pos[1], t.pos[0]), sub3(t.pos[2], t.pos[0]));
    if norm3(n) <= PLANE_DIST_EPS {
        return None;
    }
    Some((n, -dot3(n, t.pos[0])))
}

/// Interval of `t` projected on `axis` where the triangle crosses the plane
/// of the other triangle (signed distances `d`). Vertices lying on the plane
/// contribute themselves; strict sign changes contribute the edge crossings.
fn crossing_interval(t: &Tri, d: &[f64; 3], axis: usize) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
        let di = d[i];
        let dj = d[j];
        if di.abs() <= PLANE_DIST_EPS {
            // Vertex i lies on the plane: the crossing passes through it.
            lo = lo.min(t.pos[i][axis]);
            hi = hi.max(t.pos[i][axis]);
        } else if dj.abs() > PLANE_DIST_EPS && (di > 0.0) != (dj > 0.0) {
            // Edge (i, j) strictly crosses the plane.
            let param = di / (di - dj);
            let v = t.pos[i][axis] + param * (t.pos[j][axis] - t.pos[i][axis]);
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    if lo > hi {
        (0.0, 0.0) // no crossing found; harmless default
    } else {
        (lo, hi)
    }
}

/// Coplanar triangle-triangle overlap: project on the plane perpendicular to
/// the dominant component of the shared normal, then test edge pairs and
/// vertex containment (boundary contact counts as intersecting).
fn coplanar_intersect(a: &Tri, b: &Tri, n: &[f64; 3]) -> bool {
    let axis = dominant_axis(n);
    let (ax, ay) = ((axis + 1) % 3, (axis + 2) % 3);
    let a2: [(f64, f64); 3] = a.pos.map(|p| (p[ax], p[ay]));
    let b2: [(f64, f64); 3] = b.pos.map(|p| (p[ax], p[ay]));
    for i in 0..3 {
        let p1 = a2[i];
        let p2 = a2[(i + 1) % 3];
        for j in 0..3 {
            if seg2_intersect(p1, p2, b2[j], b2[(j + 1) % 3]) {
                return true;
            }
        }
    }
    point_in_tri2(a2, b2[0]) || point_in_tri2(b2, a2[0])
}

/// Index of the largest-magnitude component.
fn dominant_axis(v: &[f64; 3]) -> usize {
    let mut axis = 0;
    for k in 1..3 {
        if v[k].abs() > v[axis].abs() {
            axis = k;
        }
    }
    axis
}

/// Inclusive 2D segment-segment intersection (touching counts).
fn seg2_intersect(p1: (f64, f64), p2: (f64, f64), q1: (f64, f64), q2: (f64, f64)) -> bool {
    let d1 = orient2(q1, q2, p1);
    let d2 = orient2(q1, q2, p2);
    let d3 = orient2(p1, p2, q1);
    let d4 = orient2(p1, p2, q2);
    if ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
        && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
    {
        return true;
    }
    (d1 == 0.0 && on_seg2(q1, q2, p1))
        || (d2 == 0.0 && on_seg2(q1, q2, p2))
        || (d3 == 0.0 && on_seg2(p1, p2, q1))
        || (d4 == 0.0 && on_seg2(p1, p2, q2))
}

/// Inclusive point-in-triangle test (boundary counts as inside).
fn point_in_tri2(tri: [(f64, f64); 3], p: (f64, f64)) -> bool {
    let d1 = orient2(tri[0], tri[1], p);
    let d2 = orient2(tri[1], tri[2], p);
    let d3 = orient2(tri[2], tri[0], p);
    !((d1 < 0.0 || d2 < 0.0 || d3 < 0.0) && (d1 > 0.0 || d2 > 0.0 || d3 > 0.0))
}

/// 2D orientation determinant.
fn orient2(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}

/// Is `p` within the axis-aligned bounds of segment `ab` (with eps slack)?
fn on_seg2(a: (f64, f64), b: (f64, f64), p: (f64, f64)) -> bool {
    p.0 >= a.0.min(b.0) - PLANE_DIST_EPS
        && p.0 <= a.0.max(b.0) + PLANE_DIST_EPS
        && p.1 >= a.1.min(b.1) - PLANE_DIST_EPS
        && p.1 <= a.1.max(b.1) + PLANE_DIST_EPS
}

/// Twice the area magnitude folded into the actual area.
fn triangle_area(pos: &[[f64; 3]; 3]) -> f64 {
    norm3(cross3(sub3(pos[1], pos[0]), sub3(pos[2], pos[0]))) * 0.5
}

/// Per-triangle axis-aligned bounding box (min, max).
fn triangle_aabb(t: &Tri) -> ([f64; 3], [f64; 3]) {
    let mut min = t.pos[0];
    let mut max = t.pos[0];
    for p in &t.pos[1..] {
        for k in 0..3 {
            min[k] = min[k].min(p[k]);
            max[k] = max[k].max(p[k]);
        }
    }
    (min, max)
}

/// Inclusive AABB overlap on all three axes.
fn aabb_overlap(a: &([f64; 3], [f64; 3]), b: &([f64; 3], [f64; 3])) -> bool {
    for k in 0..3 {
        if a.1[k] < b.0[k] || b.1[k] < a.0[k] {
            return false;
        }
    }
    true
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

fn dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    norm3([b[0] - a[0], b[1] - a[1], b[2] - a[2]])
}

#[cfg(test)]
mod tests {
    use super::audit_watertight;
    use super::DegenerateReason;
    use rmsh_model::{Element, ElementType, Mesh, Node};

    const EPS: f64 = 1e-9;

    fn add_node(mesh: &mut Mesh, id: u64, x: f64, y: f64, z: f64) {
        mesh.add_node(Node::new(id, x, y, z));
    }

    fn add_tri(mesh: &mut Mesh, id: u64, a: u64, b: u64, c: u64) {
        mesh.add_element(Element::new(id, ElementType::Triangle3, vec![a, b, c]));
    }

    /// Unit cube [0,1]^3 as 12 consistently outward-oriented triangles.
    /// Node layout: 1..4 bottom (z=0) counterclockwise from above,
    /// 5..8 top (z=1) in the same order.
    fn cube_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        add_node(&mut mesh, 1, 0.0, 0.0, 0.0);
        add_node(&mut mesh, 2, 1.0, 0.0, 0.0);
        add_node(&mut mesh, 3, 1.0, 1.0, 0.0);
        add_node(&mut mesh, 4, 0.0, 1.0, 0.0);
        add_node(&mut mesh, 5, 0.0, 0.0, 1.0);
        add_node(&mut mesh, 6, 1.0, 0.0, 1.0);
        add_node(&mut mesh, 7, 1.0, 1.0, 1.0);
        add_node(&mut mesh, 8, 0.0, 1.0, 1.0);

        // bottom (normal -z)
        add_tri(&mut mesh, 1, 1, 3, 2);
        add_tri(&mut mesh, 2, 1, 4, 3);
        // top (normal +z)
        add_tri(&mut mesh, 3, 5, 6, 7);
        add_tri(&mut mesh, 4, 5, 7, 8);
        // front y=0 (normal -y)
        add_tri(&mut mesh, 5, 1, 2, 6);
        add_tri(&mut mesh, 6, 1, 6, 5);
        // back y=1 (normal +y)
        add_tri(&mut mesh, 7, 3, 4, 8);
        add_tri(&mut mesh, 8, 3, 8, 7);
        // left x=0 (normal -x)
        add_tri(&mut mesh, 9, 1, 5, 8);
        add_tri(&mut mesh, 10, 1, 8, 4);
        // right x=1 (normal +x)
        add_tri(&mut mesh, 11, 2, 3, 7);
        add_tri(&mut mesh, 12, 2, 7, 6);
        mesh
    }

    #[test]
    fn clean_unit_cube_audits_watertight() {
        let report = audit_watertight(&cube_mesh());
        assert!(
            report.is_watertight,
            "the closed consistent cube must audit watertight: {report:?}"
        );
        assert!(report.boundary_loops.is_empty());
        assert!(report.non_manifold_edges.is_empty());
        assert!(report.orientation_mismatches.is_empty());
        assert!(report.duplicate_triangles.is_empty());
        assert!(report.degenerate_triangles.is_empty());
        assert!(report.self_intersections.is_empty());
        assert!(report.hints.is_empty());
    }

    #[test]
    fn punched_cube_reports_one_triangular_hole() {
        let mut mesh = cube_mesh();
        // Drop the last triangle (element id 12, nodes 2,7,6); its three
        // edges become boundary.
        mesh.elements.pop();
        let report = audit_watertight(&mesh);
        assert!(!report.is_watertight);
        assert_eq!(report.boundary_loops.len(), 1, "{report:?}");
        let loop0 = &report.boundary_loops[0];
        assert!(loop0.closed, "the punched hole must close into a loop");
        assert_eq!(loop0.edge_count(), 3);
        // Directed along the traversal from the smallest boundary edge.
        assert_eq!(loop0.edges, vec![(2, 6), (6, 7), (7, 2)]);
        assert_eq!(loop0.node_ids, vec![2, 6, 7]);
        assert!(
            (loop0.total_length - (1.0 + 2f64.sqrt() + 1.0)).abs() < EPS,
            "loop length = 1 + sqrt(2) + 1, got {}",
            loop0.total_length
        );
        let hint = &report.hints[0];
        assert!(hint.contains("cap the boundary loop of 3 edges"), "{hint}");
        // Nothing else is wrong with the punched cube.
        assert!(report.non_manifold_edges.is_empty());
        assert!(report.orientation_mismatches.is_empty());
        assert!(report.duplicate_triangles.is_empty());
        assert!(report.self_intersections.is_empty());
    }

    #[test]
    fn flipped_triangle_cube_reports_three_orientation_mismatches() {
        let mut mesh = cube_mesh();
        // Reverse the winding of the last triangle (element id 12,
        // nodes 2,7,6 -> 6,7,2); connectivity is unchanged.
        let last = mesh.elements.last_mut().unwrap();
        last.node_ids.reverse();
        let report = audit_watertight(&mesh);
        assert!(!report.is_watertight);
        assert_eq!(
            report.orientation_mismatches.len(),
            3,
            "each of the flipped triangle's 3 edges is traversed twice in one direction: {report:?}"
        );
        let mut keys: Vec<(u64, u64)> = report
            .orientation_mismatches
            .iter()
            .map(|m| m.nodes)
            .collect();
        keys.sort_unstable();
        assert_eq!(keys, vec![(2, 6), (2, 7), (6, 7)]);
        for m in &report.orientation_mismatches {
            assert!(
                (m.forward == 2 && m.backward == 0) || (m.forward == 0 && m.backward == 2),
                "unbalanced edge {m:?}"
            );
        }
        // Connectivity is untouched: still closed and manifold.
        assert!(report.boundary_loops.is_empty());
        assert!(report.non_manifold_edges.is_empty());
        assert!(report.self_intersections.is_empty());
        assert_eq!(report.hints.len(), 3);
    }

    #[test]
    fn duplicated_triangle_cube_reports_duplicate_and_nonmanifold() {
        let mut mesh = cube_mesh();
        // Duplicate element 1 (nodes 1,3,2) as element 13 with same winding.
        add_tri(&mut mesh, 13, 1, 3, 2);
        let report = audit_watertight(&mesh);
        assert!(!report.is_watertight);

        // 1 duplicate group: elements [1, 13], nodes {1,2,3}.
        assert_eq!(report.duplicate_triangles.len(), 1, "{report:?}");
        let dup = &report.duplicate_triangles[0];
        assert_eq!(dup.nodes, [1, 2, 3]);
        assert_eq!(dup.element_ids, vec![1, 13]);

        // The duplicated triangle has 3 edges; each goes from 2 to 3
        // incidences. (The card draft expected 2; a triangle has three
        // distinct edges and every one of them becomes over-shared.)
        assert_eq!(report.non_manifold_edges.len(), 3, "{report:?}");
        let mut nm: Vec<(u64, u64)> = report
            .non_manifold_edges
            .iter()
            .map(|e| e.nodes)
            .collect();
        nm.sort_unstable();
        assert_eq!(nm, vec![(1, 2), (1, 3), (2, 3)]);
        for e in &report.non_manifold_edges {
            assert_eq!(e.incidence, 3);
        }

        // The duplicate also unbalances the 3 edges it touches.
        assert_eq!(report.orientation_mismatches.len(), 3, "{report:?}");

        // The duplicate shares ALL nodes with its twin: the adjacency
        // exemption must not flag the pair as a self-intersection.
        assert!(report.self_intersections.is_empty());
        assert!(report.boundary_loops.is_empty());
    }

    #[test]
    fn crossing_triangles_are_detected_as_self_intersection() {
        // Two triangles crossing through each other in 3D, sharing no node:
        // T1 lies in z=0, T2 lies in the vertical plane x=1.
        let mut mesh = Mesh::new();
        add_node(&mut mesh, 1, 0.0, 0.0, 0.0);
        add_node(&mut mesh, 2, 4.0, 0.0, 0.0);
        add_node(&mut mesh, 3, 0.0, 4.0, 0.0);
        add_node(&mut mesh, 4, 1.0, -1.0, -1.0);
        add_node(&mut mesh, 5, 1.0, 3.0, -1.0);
        add_node(&mut mesh, 6, 1.0, -1.0, 3.0);
        add_tri(&mut mesh, 1, 1, 2, 3);
        add_tri(&mut mesh, 2, 4, 5, 6);

        let report = audit_watertight(&mesh);
        assert_eq!(report.self_intersections.len(), 1, "{report:?}");
        let ix = report.self_intersections[0];
        assert_eq!(ix.element_a, 1);
        assert_eq!(ix.element_b, 2);
        assert!(report.hints.iter().any(|h| h.contains("element 1 and element 2")));
        assert!(!report.is_watertight);
    }

    #[test]
    fn coplanar_overlapping_triangles_are_detected() {
        // KNOWN-BAD defect class of the P0-era L-shape fixture: the concave
        // domain leak at mesh size 0.4 produces EXTRA IN-PLANE triangles that
        // overlap valid ones (drifting area). rmsh-geo cannot build that
        // polygon mesh (mesh_polygon lives in rmsh-algo and algo does not
        // depend on geo), so the equivalent overlapping coplanar soup is
        // built by hand: two same-plane triangles overlapping over a region.
        let mut mesh = Mesh::new();
        add_node(&mut mesh, 1, 0.0, 0.0, 0.0);
        add_node(&mut mesh, 2, 2.0, 0.0, 0.0);
        add_node(&mut mesh, 3, 0.0, 2.0, 0.0);
        add_node(&mut mesh, 4, 1.0, 0.0, 0.0);
        add_node(&mut mesh, 5, 3.0, 0.0, 0.0);
        add_node(&mut mesh, 6, 1.0, 2.0, 0.0);
        add_tri(&mut mesh, 1, 1, 2, 3);
        add_tri(&mut mesh, 2, 4, 5, 6);

        let report = audit_watertight(&mesh);
        assert_eq!(report.self_intersections.len(), 1, "{report:?}");
        assert_eq!(report.self_intersections[0].element_a, 1);
        assert_eq!(report.self_intersections[0].element_b, 2);
    }

    #[test]
    fn node_touching_triangles_are_exempt_from_intersection() {
        // Two triangles touching at exactly one shared node: adjacent
        // geometry, never a self-intersection.
        let mut mesh = Mesh::new();
        add_node(&mut mesh, 1, 0.0, 0.0, 0.0);
        add_node(&mut mesh, 2, 1.0, 0.0, 0.0);
        add_node(&mut mesh, 3, 0.0, 1.0, 0.0);
        add_node(&mut mesh, 4, -1.0, 0.0, 0.0);
        add_node(&mut mesh, 5, 0.0, -1.0, 0.0);
        add_tri(&mut mesh, 1, 1, 2, 3);
        add_tri(&mut mesh, 2, 1, 4, 5);

        let report = audit_watertight(&mesh);
        assert!(
            report.self_intersections.is_empty(),
            "the shared-node adjacency exemption must suppress the pair: {report:?}"
        );
        // Two lone triangles are of course not watertight (boundary edges).
        assert!(!report.is_watertight);
    }

    #[test]
    fn degenerate_triangles_are_flagged() {
        // Zero-area: three distinct but collinear nodes.
        let mut mesh = Mesh::new();
        add_node(&mut mesh, 1, 0.0, 0.0, 0.0);
        add_node(&mut mesh, 2, 1.0, 0.0, 0.0);
        add_node(&mut mesh, 3, 2.0, 0.0, 0.0);
        add_node(&mut mesh, 4, 0.0, 1.0, 0.0);
        add_tri(&mut mesh, 1, 1, 2, 3);
        let report = audit_watertight(&mesh);
        assert_eq!(report.degenerate_triangles.len(), 1, "{report:?}");
        assert_eq!(report.degenerate_triangles[0].element_id, 1);
        assert_eq!(
            report.degenerate_triangles[0].reason,
            DegenerateReason::ZeroArea
        );
        assert!(report.degenerate_triangles[0].area < 1e-12);

        // Repeated node id inside one triangle.
        let mut mesh2 = Mesh::new();
        add_node(&mut mesh2, 1, 0.0, 0.0, 0.0);
        add_node(&mut mesh2, 2, 1.0, 0.0, 0.0);
        add_node(&mut mesh2, 3, 0.0, 1.0, 0.0);
        add_tri(&mut mesh2, 7, 1, 2, 2);
        let report2 = audit_watertight(&mesh2);
        assert_eq!(report2.degenerate_triangles.len(), 1, "{report2:?}");
        assert_eq!(report2.degenerate_triangles[0].element_id, 7);
        assert_eq!(
            report2.degenerate_triangles[0].reason,
            DegenerateReason::RepeatedNode
        );
        assert!(!report2.is_watertight);
    }
}
