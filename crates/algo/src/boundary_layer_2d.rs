//! 2-D boundary layers (ROADMAP P3c): structured quad/tri element rows
//! inserted along a domain boundary, growing inward, with collision
//! handling for thin features.
//!
//! # Algorithm provenance
//!
//! The algorithm follows Gmsh's 2-D boundary-layer machinery
//! (`src/Mesh/meshGFaceBoundaryLayer.cpp` and
//! `src/geo/boundaryLayersData.cpp`, GPL-2+). It was understood from those
//! sources and re-implemented from scratch in Rust for rmsh (MIT); no code
//! text was copied. What is carried over is the *architecture*:
//!
//! - the existing face mesh's boundary discretization becomes ring 0 of the
//!   layer structure (Gmsh: the meshed `GEdge` lines);
//! - inward normals come from the boundary geometry, one per boundary edge,
//!   and per-vertex directions average the two adjacent edge normals
//!   (Gmsh `interiorNormal` / angle-bisector columns);
//! - layer elements are quads between consecutive rings where both rings
//!   move, plus triangles where the front degenerates (Gmsh L93 quads and
//!   the L106-120 transition triangle when the two columns differ in
//!   height);
//! - the old interior triangles are discarded and the region inside the
//!   final ring is retriangulated (Gmsh: `deMeshGFace` kill + remesh with
//!   the "hop" boundary, L468-470).
//!
//! # v1 scope (honest)
//!
//! - Single outer boundary loop, no holes. A mesh whose boundary does not
//!   form exactly one closed loop on the polygon is rejected.
//! - Collision handling is a *preventive local clamp*, not Gmsh's full
//!   quality-layer engine: each boundary vertex's total advance is capped
//!   at `clamp_safety * clearance(vertex)` where `clearance` is the vertex's
//!   ray-cast distance along its own offset direction to the nearest
//!   *non-adjacent* boundary segment (the first thing the advancing column
//!   would hit). With `clamp_safety <= 0.5` two opposing walls cannot
//!   cross. A ring self-intersection re-check (segment-segment tests)
//!   catches the remaining cases — notably averaged-normal corner notches —
//!   and clamps ONLY the columns involved in a crossing, with a global
//!   depth-halving as the last resort. Layer quality optimization, Gmsh's
//!   fan/cross/tail corner treatments, and per-node `hWall` profiles along
//!   edges are out of v1 scope.
//! - The input mesh contributes its *boundary discretization only*; the
//!   interior is re-seeded and retriangulated (deterministically).

use std::collections::HashMap;
use std::sync::Arc;

use glam::DVec3;
use rmsh_model::{Element, ElementType, Mesh, Node};
use thiserror::Error;

use crate::size_field::SizeField;
use crate::triangulate2d::{Polygon2D, triangulate_points};

/// Everything that can go wrong while inserting boundary layers.
#[derive(Error, Debug, Clone, PartialEq)]
pub enum BoundaryLayerError {
    /// The layer specification itself is invalid.
    #[error("invalid boundary layer spec: {0}")]
    InvalidSpec(String),
    /// The input mesh / boundary pair is unusable.
    #[error("invalid boundary layer input: {0}")]
    InvalidInput(String),
    /// Layer construction failed (collision fallback exhausted, degenerate
    /// rings, ...).
    #[error("boundary layer generation failed: {0}")]
    Generation(String),
}

/// Specification of a 2-D boundary layer.
///
/// Layer *k* (1-based) has thickness `first_height * growth^(k-1)`: a
/// geometric progression. The cumulative advance after `k` layers is
/// [`cumulative_height`].
#[derive(Debug, Clone)]
pub struct BoundaryLayer2D {
    /// Number of layers to insert (0 = no-op passthrough).
    pub layers: usize,
    /// Height of the first (innermost-contacting) layer. Ignored when
    /// [`BoundaryLayer2D::thickness_field`] is set.
    pub first_height: f64,
    /// Growth ratio between consecutive layer thicknesses. Must be finite
    /// and positive; 1.0 gives uniform layers.
    pub growth: f64,
    /// Target element size for the retriangulated interior remainder.
    /// `None` derives it from the input mesh's mean triangle edge length.
    pub interior_mesh_size: Option<f64>,
    /// Optional size field sampled at each boundary node to give the LOCAL
    /// first-layer height (Gmsh `hWall` semantics); later layers grow from
    /// that local anchor by `growth`. When set, `first_height` is unused.
    pub thickness_field: Option<Arc<dyn SizeField>>,
    /// Fraction of the local clearance a vertex may advance at most.
    /// Must lie in `(0, 0.5]`: with 0.5, two opposing walls exactly touch;
    /// the default 0.45 keeps a 10% margin.
    pub clamp_safety: f64,
}

impl BoundaryLayer2D {
    /// A spec with default collision parameters; `growth = 1.0` (uniform
    /// layers), no thickness field, interior size derived from the input
    /// mesh.
    pub fn new(layers: usize, first_height: f64, growth: f64) -> Self {
        Self {
            layers,
            first_height,
            growth,
            interior_mesh_size: None,
            thickness_field: None,
            clamp_safety: 0.45,
        }
    }

    fn validate(&self) -> Result<(), BoundaryLayerError> {
        if !self.growth.is_finite() || self.growth <= 0.0 {
            return Err(BoundaryLayerError::InvalidSpec(format!(
                "growth must be finite and positive, got {}",
                self.growth
            )));
        }
        if self.thickness_field.is_none()
            && (!self.first_height.is_finite() || self.first_height <= 0.0)
        {
            return Err(BoundaryLayerError::InvalidSpec(format!(
                "first_height must be finite and positive, got {}",
                self.first_height
            )));
        }
        if !self.clamp_safety.is_finite() || self.clamp_safety <= 0.0 || self.clamp_safety > 0.5 {
            return Err(BoundaryLayerError::InvalidSpec(format!(
                "clamp_safety must be finite in (0, 0.5], got {}",
                self.clamp_safety
            )));
        }
        if let Some(s) = self.interior_mesh_size {
            if !s.is_finite() || s <= 0.0 {
                return Err(BoundaryLayerError::InvalidSpec(format!(
                    "interior_mesh_size must be finite and positive, got {s}"
                )));
            }
        }
        Ok(())
    }
}

/// Cumulative inward advance after `k` layers for first height `h` and
/// growth ratio `g`: `h * (g^k - 1) / (g - 1)`, or `h * k` when `g == 1`.
pub fn cumulative_height(h: f64, g: f64, k: usize) -> f64 {
    if k == 0 {
        return 0.0;
    }
    if g == 1.0 {
        h * k as f64
    } else {
        h * (g.powi(k as i32) - 1.0) / (g - 1.0)
    }
}

// ─── Small 2-D geometry helpers ──────────────────────────────────────────────

/// Twice the signed area of triangle (a, b, c); positive for CCW.
fn orient(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// Shoelace sum over a closed polyline (twice its signed area).
fn signed_area2(pts: &[[f64; 2]]) -> f64 {
    let n = pts.len();
    let mut s = 0.0;
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        s += a[0] * b[1] - b[0] * a[1];
    }
    s
}

fn norm2(v: [f64; 2]) -> f64 {
    (v[0] * v[0] + v[1] * v[1]).sqrt()
}

/// Squared distance from point `p` to segment `(a, b)`.
fn point_segment_distance_sq(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ap = [p[0] - a[0], p[1] - a[1]];
    let len_sq = ab[0] * ab[0] + ab[1] * ab[1];
    if len_sq <= 0.0 {
        return ap[0] * ap[0] + ap[1] * ap[1];
    }
    let t = ((ap[0] * ab[0] + ap[1] * ab[1]) / len_sq).clamp(0.0, 1.0);
    let dx = ap[0] - t * ab[0];
    let dy = ap[1] - t * ab[1];
    dx * dx + dy * dy
}

fn on_segment(a: [f64; 2], p: [f64; 2], b: [f64; 2]) -> bool {
    let cross = (p[0] - a[0]) * (b[1] - a[1]) - (p[1] - a[1]) * (b[0] - a[0]);
    let dot = (p[0] - a[0]) * (b[0] - a[0]) + (p[1] - a[1]) * (b[1] - a[1]);
    cross == 0.0 && dot >= 0.0 && dot <= (b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)
}

/// Inclusive segment-segment intersection (proper crossings and
/// endpoint-on-segment touches). Used for the offset-ring self-intersection
/// check; adjacent segments sharing an endpoint are excluded by the caller.
fn segments_intersect(p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], p4: [f64; 2]) -> bool {
    let o1 = orient(p1, p2, p3);
    let o2 = orient(p1, p2, p4);
    let o3 = orient(p3, p4, p1);
    let o4 = orient(p3, p4, p2);
    if o1 * o2 < 0.0 && o3 * o4 < 0.0 {
        return true;
    }
    (o1 == 0.0 && on_segment(p1, p3, p2))
        || (o2 == 0.0 && on_segment(p1, p4, p2))
        || (o3 == 0.0 && on_segment(p3, p1, p4))
        || (o4 == 0.0 && on_segment(p3, p2, p4))
}

/// Position of a mesh node as an XY pair.
fn pos2(mesh: &Mesh, id: u64) -> Result<[f64; 2], BoundaryLayerError> {
    mesh.nodes
        .get(&id)
        .map(|n| [n.position.x, n.position.y])
        .ok_or_else(|| BoundaryLayerError::InvalidInput(format!("node {id} missing from mesh")))
}

// ─── Boundary loop extraction ────────────────────────────────────────────────

/// Extract the single closed boundary loop of a triangle mesh as an ordered
/// node-id list. The mesh must be a pure triangle mesh whose 1-incidence
/// boundary forms one closed manifold loop.
fn extract_boundary_loop(mesh: &Mesh) -> Result<Vec<u64>, BoundaryLayerError> {
    let mut edge_count: HashMap<[u64; 2], usize> = HashMap::new();
    for e in &mesh.elements {
        match e.etype {
            ElementType::Triangle3 => {
                for k in 0..3 {
                    let u = e.node_ids[k];
                    let v = e.node_ids[(k + 1) % 3];
                    *edge_count.entry([u.min(v), u.max(v)]).or_insert(0) += 1;
                }
            }
            ElementType::Line2 | ElementType::Point1 => {}
            other => {
                return Err(BoundaryLayerError::InvalidInput(format!(
                    "input mesh must be a pure triangle mesh, found {other:?}"
                )));
            }
        }
    }
    if edge_count.is_empty() {
        return Err(BoundaryLayerError::InvalidInput(
            "input mesh has no triangles".to_string(),
        ));
    }

    // Boundary adjacency (edges used by exactly one triangle).
    let mut adj: HashMap<u64, Vec<u64>> = HashMap::new();
    for (key, count) in &edge_count {
        if *count == 1 {
            adj.entry(key[0]).or_default().push(key[1]);
            adj.entry(key[1]).or_default().push(key[0]);
        }
    }
    if adj.is_empty() {
        return Err(BoundaryLayerError::InvalidInput(
            "mesh boundary is empty (closed surface?)".to_string(),
        ));
    }
    for (id, nbrs) in &adj {
        if nbrs.len() != 2 {
            return Err(BoundaryLayerError::InvalidInput(format!(
                "boundary node {id} has {} boundary edges (not a single manifold loop)",
                nbrs.len()
            )));
        }
    }

    // Walk the loop from the smallest node id (deterministic).
    let start = *adj.keys().min().expect("non-empty");
    let mut loop_ids = vec![start];
    let mut prev: Option<u64> = None;
    let mut cur = start;
    loop {
        let mut next: Option<u64> = None;
        for &nbr in &adj[&cur] {
            if Some(nbr) != prev {
                next = Some(nbr);
                break;
            }
        }
        let next = next.ok_or_else(|| {
            BoundaryLayerError::InvalidInput("boundary loop walk dead-ended".to_string())
        })?;
        if next == start {
            break;
        }
        loop_ids.push(next);
        prev = Some(cur);
        cur = next;
        if loop_ids.len() > adj.len() {
            return Err(BoundaryLayerError::InvalidInput(
                "boundary loop walk did not close".to_string(),
            ));
        }
    }
    if loop_ids.len() < 3 {
        return Err(BoundaryLayerError::InvalidInput(
            "boundary loop has fewer than 3 nodes".to_string(),
        ));
    }
    Ok(loop_ids)
}

/// Every loop node must lie on the polygon outline (within a small fraction
/// of the polygon's bbox diagonal).
fn validate_loop_on_polygon(
    pts: &[[f64; 2]],
    boundary: &Polygon2D,
) -> Result<(), BoundaryLayerError> {
    if boundary.vertices.len() < 3 {
        return Err(BoundaryLayerError::InvalidInput(
            "boundary polygon has fewer than 3 vertices".to_string(),
        ));
    }
    let (bmin, bmax) = boundary.bounding_box();
    let diag = ((bmax[0] - bmin[0]).powi(2) + (bmax[1] - bmin[1]).powi(2)).sqrt();
    let tol = 1e-6 * diag.max(1e-30);
    let nv = boundary.vertices.len();
    for p in pts {
        let mut best = f64::INFINITY;
        for i in 0..nv {
            let a = boundary.vertices[i];
            let b = boundary.vertices[(i + 1) % nv];
            best = best.min(point_segment_distance_sq(*p, a, b));
        }
        if best.sqrt() > tol {
            return Err(BoundaryLayerError::InvalidInput(format!(
                "mesh boundary node ({}, {}) does not lie on the boundary polygon",
                p[0], p[1]
            )));
        }
    }
    Ok(())
}

// ─── Ring construction ───────────────────────────────────────────────────────

/// One candidate layer front: per-vertex advance and position.
struct Ring {
    /// Advance of each vertex (min of scaled cumulative height and cap).
    advance: Vec<f64>,
    /// Ring node positions (loop order).
    positions: Vec<[f64; 2]>,
}

/// Build rings 0..=layers for a global depth `factor` (the collision
/// fallback bisects this). Ring k advance = min(factor * cumulative, cap),
/// which is monotone non-decreasing in k for every vertex.
fn build_rings(
    pts: &[[f64; 2]],
    dirs: &[[f64; 2]],
    h0: &[f64],
    caps: &[f64],
    growth: f64,
    layers: usize,
    factor: f64,
) -> Vec<Ring> {
    let n = pts.len();
    let mut rings = Vec::with_capacity(layers + 1);
    rings.push(Ring {
        advance: vec![0.0; n],
        positions: pts.to_vec(),
    });
    for k in 1..=layers {
        let advance: Vec<f64> = (0..n)
            .map(|i| (factor * cumulative_height(h0[i], growth, k)).min(caps[i]))
            .collect();
        let positions: Vec<[f64; 2]> = (0..n)
            .map(|i| [
                pts[i][0] + dirs[i][0] * advance[i],
                pts[i][1] + dirs[i][1] * advance[i],
            ])
            .collect();
        rings.push(Ring { advance, positions });
    }
    rings
}

/// First self-crossing segment pair of a ring polyline (O(n^2),
/// non-adjacent pairs only; zero-length segments from stalled columns are
/// skipped), in deterministic scan order.
fn find_ring_crossing(positions: &[[f64; 2]], scale: f64) -> Option<(usize, usize)> {
    let n = positions.len();
    let zero = 1e-14 * scale;
    let seg = |i: usize| (positions[i], positions[(i + 1) % n]);
    let nondegenerate: Vec<usize> = (0..n)
        .filter(|&i| {
            let (a, b) = seg(i);
            norm2([b[0] - a[0], b[1] - a[1]]) > zero
        })
        .collect();
    for (ai, &i) in nondegenerate.iter().enumerate() {
        for &j in nondegenerate.iter().skip(ai + 1) {
            // skip cyclic neighbours (they share an endpoint by construction)
            let cyclic_adjacent = (i + 1) % n == j || (j + 1) % n == i;
            if cyclic_adjacent || (i == 0 && j == n - 1) {
                continue;
            }
            let (a1, a2) = seg(i);
            let (b1, b2) = seg(j);
            if segments_intersect(a1, a2, b1, b2) {
                return Some((i, j));
            }
        }
    }
    None
}

/// Area measure for the positive-area contract: |2A| for a triangle,
/// signed shoelace for a quad (CCW expected).
fn element_area2(etype: ElementType, p: &[[f64; 2]]) -> f64 {
    match etype {
        ElementType::Quad4 => {
            (p[0][0] * p[1][1] - p[1][0] * p[0][1])
                + (p[1][0] * p[2][1] - p[2][0] * p[1][1])
                + (p[2][0] * p[3][1] - p[3][0] * p[2][1])
                + (p[3][0] * p[0][1] - p[0][0] * p[3][1])
        }
        _ => orient(p[0], p[1], p[2]).abs(),
    }
}

/// First intersection parameter t > 0 of the ray `p + t*d` with segment
/// `(a, b)`, or infinity when the ray misses. Grazing a segment endpoint
/// counts as a hit (both neighbouring segments report it).
fn ray_segment_hit(p: [f64; 2], d: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let e = [b[0] - a[0], b[1] - a[1]];
    let denom = d[0] * e[1] - d[1] * e[0];
    if denom.abs() <= 1e-15 {
        return f64::INFINITY;
    }
    let ap = [a[0] - p[0], a[1] - p[1]];
    let t = (ap[0] * e[1] - ap[1] * e[0]) / denom;
    let s = (ap[0] * d[1] - ap[1] * d[0]) / denom;
    if t > 0.0 && (0.0..=1.0).contains(&s) {
        t
    } else {
        f64::INFINITY
    }
}

/// Insert boundary layers into a triangulated domain (see module docs).
pub fn insert_boundary_layers(
    mesh: &Mesh,
    boundary: &Polygon2D,
    spec: &BoundaryLayer2D,
) -> Result<Mesh, BoundaryLayerError> {
    spec.validate()?;
    if spec.layers == 0 {
        return Ok(mesh.clone());
    }

    // ── Ring 0: the mesh's own boundary discretization, canonicalized CCW ──
    let loop_ids = extract_boundary_loop(mesh)?;
    let mut pts: Vec<[f64; 2]> = loop_ids
        .iter()
        .map(|&id| pos2(mesh, id))
        .collect::<Result<Vec<_>, _>>()?;
    if signed_area2(&pts) < 0.0 {
        pts.reverse();
    }
    validate_loop_on_polygon(&pts, boundary)?;
    let n = pts.len();
    let scale: f64 = {
        let (mut mn, mut mx) = ([f64::MAX; 2], [f64::MIN; 2]);
        for p in &pts {
            mn[0] = mn[0].min(p[0]);
            mn[1] = mn[1].min(p[1]);
            mx[0] = mx[0].max(p[0]);
            mx[1] = mx[1].max(p[1]);
        }
        norm2([mx[0] - mn[0], mx[1] - mn[1]]).max(1e-30)
    };

    // ── Inward normals: per-edge LEFT normals of the CCW loop ─────────────
    let mut edge_normals = Vec::with_capacity(n);
    let mut edge_lengths = Vec::with_capacity(n);
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        let d = [b[0] - a[0], b[1] - a[1]];
        let len = norm2(d);
        if len <= 1e-14 * scale {
            return Err(BoundaryLayerError::InvalidInput(format!(
                "degenerate boundary edge {i} (length {len})"
            )));
        }
        edge_normals.push([-d[1] / len, d[0] / len]);
        edge_lengths.push(len);
    }
    // Per-vertex offset directions: normalized sum of adjacent edge normals
    // (angle bisector); fall back to the longer adjacent edge's normal when
    // the sum degenerates (needle vertex).
    let mut dirs = Vec::with_capacity(n);
    for i in 0..n {
        let prev = edge_normals[(i + n - 1) % n];
        let next = edge_normals[i];
        let s = [prev[0] + next[0], prev[1] + next[1]];
        let len = norm2(s);
        if len > 1e-12 {
            dirs.push([s[0] / len, s[1] / len]);
        } else if edge_lengths[(i + n - 1) % n] >= edge_lengths[i] {
            dirs.push(prev);
        } else {
            dirs.push(next);
        }
    }

    // ── Local first-layer heights (SizeField overrides the constant) ──────
    let h0: Vec<f64> = (0..n)
        .map(|i| match &spec.thickness_field {
            Some(field) => {
                let v = field.size_at(DVec3::new(pts[i][0], pts[i][1], 0.0));
                if !v.is_finite() || v <= 0.0 {
                    return Err(BoundaryLayerError::Generation(format!(
                        "thickness field answered {v} at boundary node {i}"
                    )));
                }
                Ok(v)
            }
            None => Ok(spec.first_height),
        })
        .collect::<Result<Vec<_>, _>>()?;

    // ── Preventive collision caps (ray clearance along the offset direction)
    // Clearance of vertex i = distance along its offset ray to the nearest
    // NON-ADJACENT boundary segment (the first thing the advancing column
    // would hit). Cap = clamp_safety * clearance; with clamp_safety <= 0.5
    // two opposing walls cannot cross. Rays that miss everything give an
    // infinite cap; the ring self-intersection re-check below remains the
    // safety net for the estimate's blind spots.
    let caps: Vec<f64> = (0..n)
        .map(|i| {
            let mut clearance = f64::INFINITY;
            for s in 0..n {
                // segment (s, s+1) is adjacent to vertex i iff it touches it
                if s == i || s == (i + n - 1) % n {
                    continue;
                }
                let t = ray_segment_hit(pts[i], dirs[i], pts[s], pts[(s + 1) % n]);
                if t < clearance {
                    clearance = t;
                }
            }
            clearance * spec.clamp_safety
        })
        .collect();

    // ── Rings with the self-intersection re-check: local clamp first, then
    // a global depth-halving last resort ──────────────────────────────────
    // A detected ring crossing clamps the advance caps of ONLY the columns
    // involved in the crossing ("just below the collision"); every ring is
    // rebuilt so per-column advances stay monotone. Straight columns keep
    // their full spec depth. If local clamping does not settle, the global
    // depth factor halves.
    let mut caps = caps;
    let mut rounds = 0usize;
    let mut factor = 1.0;
    let rings = loop {
        let candidate = build_rings(&pts, &dirs, &h0, &caps, spec.growth, spec.layers, factor);
        let crossings: Vec<(usize, usize)> = candidate
            .iter()
            .skip(1)
            .filter_map(|r| find_ring_crossing(&r.positions, scale))
            .collect();
        if crossings.is_empty() {
            break candidate;
        }
        for &(s1, s2) in &crossings {
            for &v in &[s1, (s1 + 1) % n, s2, (s2 + 1) % n] {
                caps[v] *= 0.5;
            }
        }
        rounds += 1;
        if rounds > 40 {
            factor *= 0.5;
            if factor < 1e-3 {
                return Err(BoundaryLayerError::Generation(
                    "offset rings self-intersect at every probed depth".to_string(),
                ));
            }
        }
    };

    // ── Node allocation (stalled columns reuse the previous ring's node) ──
    // Node ids are (index into node_positions) + 1 throughout.
    let mut node_positions: Vec<[f64; 2]> = pts.clone();
    let mut ring_node_ids: Vec<Vec<u64>> = vec![Vec::with_capacity(n); rings.len()];
    ring_node_ids[0] = (1..=n as u64).collect();
    for k in 1..rings.len() {
        let reuse = ring_node_ids[k - 1].clone();
        for i in 0..n {
            if rings[k].advance[i] == rings[k - 1].advance[i] {
                ring_node_ids[k].push(reuse[i]);
            } else {
                ring_node_ids[k].push((node_positions.len() + 1) as u64);
                node_positions.push(rings[k].positions[i]);
            }
        }
    }

    // ── Output mesh: nodes first (allocation order), then elements ────────
    let mut out = Mesh::new();
    for (idx, p) in node_positions.iter().enumerate() {
        out.add_node(Node::new(idx as u64 + 1, p[0], p[1], 0.0));
    }
    let mut elem_id: u64 = 1;
    let mut add_element = |out: &mut Mesh, etype: ElementType, ids: Vec<u64>| -> Result<(), BoundaryLayerError> {
        let p: Vec<[f64; 2]> = ids
            .iter()
            .map(|&id| node_positions[(id - 1) as usize])
            .collect();
        let area2 = element_area2(etype, &p);
        if area2 <= 1e-12 * scale * scale {
            return Err(BoundaryLayerError::Generation(format!(
                "degenerate layer element {ids:?} (area {area2})"
            )));
        }
        out.add_element(Element::new(elem_id, etype, ids));
        elem_id += 1;
        Ok(())
    };
    for k in 1..rings.len() {
        for i in 0..n {
            let j = (i + 1) % n;
            let o1 = ring_node_ids[k - 1][i];
            let o2 = ring_node_ids[k - 1][j];
            let i1 = ring_node_ids[k][i];
            let i2 = ring_node_ids[k][j];
            let moved_a = i1 != o1;
            let moved_b = i2 != o2;
            match (moved_a, moved_b) {
                (true, true) => add_element(&mut out, ElementType::Quad4, vec![o1, o2, i2, i1])?,
                (true, false) => {
                    add_element(&mut out, ElementType::Triangle3, vec![o1, o2, i1])?
                }
                (false, true) => {
                    add_element(&mut out, ElementType::Triangle3, vec![o1, o2, i2])?
                }
                (false, false) => {}
            }
        }
    }

    // ── Interior remainder: retriangulate inside the final ring ───────────
    let inner = &rings[rings.len() - 1];
    let inner_ids = &ring_node_ids[rings.len() - 1];
    if signed_area2(&inner.positions) <= 1e-12 * scale * scale {
        return Err(BoundaryLayerError::Generation(
            "final offset ring is degenerate".to_string(),
        ));
    }
    let remainder = Polygon2D::new(inner.positions.clone());

    let interior_size = match spec.interior_mesh_size {
        Some(s) => s,
        None => {
            // mean triangle edge length of the input mesh
            let mut sum = 0.0;
            let mut count = 0usize;
            for e in &mesh.elements {
                if e.etype != ElementType::Triangle3 {
                    continue;
                }
                for k in 0..3 {
                    let a = pos2(mesh, e.node_ids[k])?;
                    let b = pos2(mesh, e.node_ids[(k + 1) % 3])?;
                    sum += norm2([b[0] - a[0], b[1] - a[1]]);
                    count += 1;
                }
            }
            if count == 0 {
                return Err(BoundaryLayerError::InvalidInput(
                    "input mesh has no triangles".to_string(),
                ));
            }
            (sum / count as f64).max(1e-12 * scale)
        }
    };

    // Staggered interior seed grid (same pattern as mesh_polygon), dropping
    // seeds that land too close to a ring node.
    let (bb_min, bb_max) = remainder.bounding_box();
    let mut seeds: Vec<[f64; 2]> = Vec::new();
    let min_seed_dist_sq = (0.35 * interior_size).powi(2);
    let mut iy = 1usize;
    let mut py = bb_min[1] + interior_size;
    while py < bb_max[1] {
        let offset_x = if iy % 2 == 0 { 0.0 } else { interior_size * 0.5 };
        let mut px = bb_min[0] + offset_x + interior_size;
        while px < bb_max[0] {
            if remainder.contains([px, py]) {
                let p = [px, py];
                let too_close = inner
                    .positions
                    .iter()
                    .any(|q| (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) < min_seed_dist_sq);
                if !too_close {
                    seeds.push(p);
                }
            }
            px += interior_size;
        }
        py += interior_size * 0.866_025_403_784_438_6; // sqrt(3)/2 row spacing
        iy += 1;
    }

    // Triangulate ring points + seeds; keep inside triangles, CCW-oriented.
    let mut all_pts = inner.positions.clone();
    all_pts.extend_from_slice(&seeds);
    let tris = triangulate_points(&all_pts);
    let mut pt_to_node: HashMap<usize, u64> = HashMap::new();
    for (i, &id) in inner_ids.iter().enumerate() {
        pt_to_node.insert(i, id);
    }
    for tri in tris {
        let a = all_pts[tri[0]];
        let b = all_pts[tri[1]];
        let c = all_pts[tri[2]];
        let centroid = [
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
        ];
        if !remainder.contains(centroid) {
            continue;
        }
        let area2 = orient(a, b, c);
        // degenerate triangles contribute zero area; dropping them keeps the
        // area-coverage contract exact
        if area2.abs() <= 1e-12 * scale * scale {
            continue;
        }
        let (t0, t1, t2) = if area2 < 0.0 { (tri[0], tri[2], tri[1]) } else { (tri[0], tri[1], tri[2]) };
        let mut ids = Vec::with_capacity(3);
        for &pi in &[t0, t1, t2] {
            let id = match pt_to_node.get(&pi) {
                Some(&id) => id,
                None => {
                    let id = (node_positions.len() + 1) as u64;
                    node_positions.push(all_pts[pi]);
                    out.add_node(Node::new(id, all_pts[pi][0], all_pts[pi][1], 0.0));
                    pt_to_node.insert(pi, id);
                    id
                }
            };
            ids.push(id);
        }
        let p: Vec<[f64; 2]> = ids
            .iter()
            .map(|&id| node_positions[(id - 1) as usize])
            .collect();
        if element_area2(ElementType::Triangle3, &p) <= 1e-12 * scale * scale {
            continue;
        }
        out.add_element(Element::new(elem_id, ElementType::Triangle3, ids));
        elem_id += 1;
    }

    Ok(out)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::size_field::FieldFn;
    use crate::triangulate2d::mesh_polygon;

    const EPS: f64 = 1e-9;

    fn rect_mesh(w: f64, h: f64, size: f64) -> (Polygon2D, Mesh) {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [w, 0.0], [w, h], [0.0, h]]);
        let mesh = mesh_polygon(&poly, size).expect("rect mesh");
        (poly, mesh)
    }

    /// Every output element must have strictly positive signed area (CCW).
    fn assert_all_positive_area(mesh: &Mesh) {
        let area2 = |ids: &[u64]| -> f64 {
            let p: Vec<[f64; 2]> = ids
                .iter()
                .map(|id| {
                    let n = mesh.nodes.get(id).expect("node exists");
                    [n.position.x, n.position.y]
                })
                .collect();
            (p[1][0] - p[0][0]) * (p[2][1] - p[0][1])
                - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0])
        };
        for e in &mesh.elements {
            match e.etype {
                ElementType::Triangle3 => {
                    assert!(
                        area2(&e.node_ids) > 0.0,
                        "inverted/degenerate triangle {:?}",
                        e.node_ids
                    );
                }
                ElementType::Quad4 => {
                    let p: Vec<[f64; 2]> = e
                        .node_ids
                        .iter()
                        .map(|id| {
                            let n = mesh.nodes.get(id).expect("node exists");
                            [n.position.x, n.position.y]
                        })
                        .collect();
                    let shoelace = (p[0][0] * p[1][1] - p[1][0] * p[0][1])
                        + (p[1][0] * p[2][1] - p[2][0] * p[1][1])
                        + (p[2][0] * p[3][1] - p[3][0] * p[2][1])
                        + (p[3][0] * p[0][1] - p[0][0] * p[3][1]);
                    assert!(shoelace > 0.0, "inverted/degenerate quad {:?}", e.node_ids);
                }
                _ => panic!("unexpected element type {:?}", e.etype),
            }
        }
    }

    /// Total covered area of all 2-D elements.
    fn total_area(mesh: &Mesh) -> f64 {
        mesh.elements
            .iter()
            .map(|e| {
                let p: Vec<[f64; 2]> = e
                    .node_ids
                    .iter()
                    .map(|id| {
                        let n = mesh.nodes.get(id).expect("node exists");
                        [n.position.x, n.position.y]
                    })
                    .collect();
                let area2 = match e.etype {
                    ElementType::Triangle3 => {
                        (p[1][0] - p[0][0]) * (p[2][1] - p[0][1])
                            - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0])
                    }
                    ElementType::Quad4 => {
                        (p[0][0] * p[1][1] - p[1][0] * p[0][1])
                            + (p[1][0] * p[2][1] - p[2][0] * p[1][1])
                            + (p[2][0] * p[3][1] - p[3][0] * p[2][1])
                            + (p[3][0] * p[0][1] - p[0][0] * p[3][1])
                    }
                    _ => panic!("unexpected element type"),
                };
                area2.abs() / 2.0
            })
            .sum()
    }

    /// Is there a node at position `p` within `eps`?
    fn has_node_at(mesh: &Mesh, x: f64, y: f64, eps: f64) -> bool {
        mesh.nodes.values().any(|n| {
            (n.position.x - x).abs() <= eps && (n.position.y - y).abs() <= eps
        })
    }

    /// Is `p` within `tol` of the polygon outline?
    fn near_boundary(poly: &Polygon2D, p: [f64; 2], tol: f64) -> bool {
        let nv = poly.vertices.len();
        (0..nv).any(|i| {
            let a = poly.vertices[i];
            let b = poly.vertices[(i + 1) % nv];
            point_segment_distance_sq(p, a, b) <= tol * tol
        })
    }

    fn count_at_y(mesh: &Mesh, y: f64, eps: f64) -> usize {
        mesh.nodes
            .values()
            .filter(|n| (n.position.y - y).abs() <= eps)
            .count()
    }

    fn count_type(mesh: &Mesh, t: ElementType) -> usize {
        mesh.elements.iter().filter(|e| e.etype == t).count()
    }

    // ── Pin: geometric progression helper ────────────────────────────────────

    #[test]
    fn cumulative_height_geometric_progression_pinned() {
        // thickness of layer k = first_height * growth^(k-1); cumulative:
        (0.05 * (1.0)).assert_eq(cumulative_height(0.05, 2.0, 1));
        (0.05 * (1.0 + 2.0)).assert_eq(cumulative_height(0.05, 2.0, 2));
        (0.05 * (1.0 + 2.0 + 4.0)).assert_eq(cumulative_height(0.05, 2.0, 3));
        // uniform growth
        (0.1 * 4.0).assert_eq(cumulative_height(0.1, 1.0, 4));
        (0.0).assert_eq(cumulative_height(0.1, 1.0, 0));
    }

    trait AssertF64 {
        fn assert_eq(self, other: f64);
    }
    impl AssertF64 for f64 {
        fn assert_eq(self, other: f64) {
            assert!((self - other).abs() < 1e-12, "{self} != {other}");
        }
    }

    // ── Pin: orientation convention ──────────────────────────────────────────

    /// CCW boundary: inward normal of the bottom edge (a->b = +x) is the
    /// LEFT normal (+y). The layer grows to +y off the bottom edge, and a
    /// CW `boundary` argument changes nothing (the loop is canonicalized).
    #[test]
    fn orientation_convention_left_normal_is_inward() {
        let (poly, mesh) = rect_mesh(4.0, 1.0, 0.4);
        let spec = BoundaryLayer2D {
            interior_mesh_size: Some(0.3),
            ..BoundaryLayer2D::new(1, 0.05, 1.0)
        };
        let out = insert_boundary_layers(&mesh, &poly, &spec).expect("layers insert");

        // 9 interior bottom-edge boundary nodes (x = 0.4 .. 3.6): their
        // first-layer nodes sit at y = +first_height exactly.
        let bottom_boundary = mesh
            .nodes
            .values()
            .filter(|n| n.position.y.abs() < 1e-12 && n.position.x > 1e-12 && 4.0 - n.position.x > 1e-12)
            .count();
        assert_eq!(bottom_boundary, 9, "expected 9 interior bottom boundary nodes");
        assert_eq!(count_at_y(&out, 0.05, EPS), bottom_boundary);
        for k in 1..=9usize {
            let x = 0.4 * k as f64;
            assert!(has_node_at(&out, x, 0.05, EPS), "no ring node at ({x}, 0.05)");
        }
        // nothing grew outward (inward = +y for the CCW bottom edge)
        assert!(
            out.nodes.values().all(|n| n.position.y >= -1e-12),
            "a node lies below the boundary: outward growth"
        );

        // CW boundary argument: identical structure, same growth side.
        let poly_cw = Polygon2D::new(vec![
            [0.0, 1.0],
            [4.0, 1.0],
            [4.0, 0.0],
            [0.0, 0.0],
        ]);
        let out_cw = insert_boundary_layers(&mesh, &poly_cw, &spec).expect("cw layers insert");
        assert_eq!(out_cw.node_count(), out.node_count());
        assert_eq!(count_type(&out_cw, ElementType::Quad4), count_type(&out, ElementType::Quad4));
        assert_eq!(count_at_y(&out_cw, 0.05, EPS), bottom_boundary);
    }

    // ── Pin: sharp corner, averaged-normal offset ────────────────────────────

    #[test]
    fn sharp_corner_averaged_normal_offset_pinned() {
        // CCW triangle with a 45 degree corner at the origin (and a sharp
        // spike at (1,0)). Depth is kept below every column's ray clearance
        // (the tightest is 0.00526 near the spike, cap 0.00237) so the
        // offset ring stays strictly inside and CONVEX — making the exact
        // area pin meaningful. The collision path itself is exercised by
        // the sliver test.
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [0.05, 0.05]]);
        let mesh = mesh_polygon(&poly, 0.1).expect("triangle mesh");
        let spec = BoundaryLayer2D {
            interior_mesh_size: Some(0.1),
            ..BoundaryLayer2D::new(1, 0.002, 1.4)
        };
        let out = insert_boundary_layers(&mesh, &poly, &spec).expect("layers insert");

        // Corner (0,0): n_prev = left normal of closing edge
        // (0.05,0.05)->(0,0) = (0.70710678, -0.70710678);
        // n_next = left normal of (0,0)->(1,0) = (0, 1);
        // averaged + normalized = (0.92387953, 0.38268343).
        let d = (0.923_879_532_511_286_7, 0.382_683_432_365_089_8);
        assert!(
            has_node_at(&out, 0.002 * d.0, 0.002 * d.1, EPS),
            "corner ring node not at the averaged-normal offset"
        );
        // every ring node stays strictly inside the polygon
        for n in out.nodes.values() {
            assert!(
                poly.contains([n.position.x, n.position.y])
                    || near_boundary(&poly, [n.position.x, n.position.y], 1e-9),
                "ring/interior node ({}, {}) left the domain",
                n.position.x,
                n.position.y
            );
        }
        assert_all_positive_area(&out);
        assert!((total_area(&out) - 0.025).abs() < EPS, "area {}", total_area(&out));
    }

    // ── Pin: collision clamp on a thin sliver ────────────────────────────────

    #[test]
    fn sliver_collision_clamped_no_overlap_no_inversion() {
        // 1 x 0.1 sliver: 4 layers of 0.04 growing 1.2 would reach 0.3373
        // total per side; the domain only allows 0.045 (0.45 * 0.1 gap).
        let (poly, mesh) = rect_mesh(1.0, 0.1, 0.2);
        let spec = BoundaryLayer2D {
            interior_mesh_size: Some(0.05),
            ..BoundaryLayer2D::new(4, 0.04, 1.2)
        };
        let out = insert_boundary_layers(&mesh, &poly, &spec).expect("layers insert");

        // Boundary loop: 5 + 1 + 5 + 1 = 12 nodes. Ray clearance is 0.1 for
        // the long-edge vertices (0.1414 along the corner bisectors):
        // cap 0.045 / 0.0636. Ring1 advances 0.04, ring2 clamps (0.045 on
        // the edges, 0.0636 at the corners = the miter points), rings 3-4
        // stall and reuse the ring2 nodes.
        assert_eq!(count_type(&out, ElementType::Quad4), 12 * 2);
        // Interior remainder: thin convex 12-gon ring -> exactly 12 - 2 = 10
        // fan triangles, no seeds fit. No transition triangles (uniform
        // stall across every column).
        assert_eq!(count_type(&out, ElementType::Triangle3), 10);
        assert_eq!(out.node_count(), 12 * 3);

        // No overlap: every node stays strictly inside the sliver, and the
        // clamped fronts stop short of the midline y = 0.05 from both sides.
        for n in out.nodes.values() {
            let (x, y) = (n.position.x, n.position.y);
            assert!(x >= -EPS && x <= 1.0 + EPS && y >= -EPS && y <= 0.1 + EPS);
            let in_forbidden_band = y > 0.046 && y < 0.054;
            assert!(
                !in_forbidden_band,
                "node at ({x},{y}) crossed the clamped mid-band"
            );
        }
        assert_all_positive_area(&out);
        assert!((total_area(&out) - 0.1).abs() < EPS, "area {}", total_area(&out));
    }

    // ── Pin: SizeField-driven local thickness ────────────────────────────────

    #[test]
    fn size_field_driven_local_thickness_pinned() {
        let (poly, mesh) = rect_mesh(4.0, 1.0, 0.5);
        // h = 0.01 on the left half, 0.03 on the right half.
        let field = FieldFn::new(|p: DVec3| if p.x <= 2.0 { 0.01 } else { 0.03 });
        let spec = BoundaryLayer2D {
            thickness_field: Some(Arc::new(field)),
            interior_mesh_size: Some(0.3),
            ..BoundaryLayer2D::new(2, 999.0, 1.0) // first_height must be ignored
        };
        let out = insert_boundary_layers(&mesh, &poly, &spec).expect("layers insert");

        // Boundary nodes at x = 1.0 and x = 3.0 (exact grid points of the
        // 0.5 discretization): local first heights 0.01 / 0.03, then grow
        // by growth = 1.0 -> ring2 at 0.02 / 0.06.
        assert!(has_node_at(&out, 1.0, 0.01, EPS), "left-edge ring1");
        assert!(has_node_at(&out, 1.0, 0.02, EPS), "left-edge ring2");
        assert!(has_node_at(&out, 3.0, 0.03, EPS), "right-edge ring1");
        assert!(has_node_at(&out, 3.0, 0.06, EPS), "right-edge ring2");
        assert_all_positive_area(&out);
    }

    // ── Pin: CW boundary loop from the mesh is canonicalized ─────────────────

    #[test]
    fn cw_boundary_loop_canonicalized_to_ccw() {
        // Hand-built square mesh (5 triangles, fan around an interior
        // point). The boundary loop walk direction is an implementation
        // detail; the result must always grow inward (+y off the bottom).
        let pts: Vec<[f64; 2]> = vec![
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.5],
        ];
        let tris: [[usize; 3]; 5] = [[0, 1, 5], [1, 2, 5], [2, 3, 5], [3, 4, 5], [4, 0, 5]];
        let mut mesh = Mesh::new();
        for (i, p) in pts.iter().enumerate() {
            mesh.add_node(Node::new(i as u64 + 1, p[0], p[1], 0.0));
        }
        for (i, t) in tris.iter().enumerate() {
            mesh.add_element(Element::new(
                i as u64 + 1,
                ElementType::Triangle3,
                t.iter().map(|&k| k as u64 + 1).collect(),
            ));
        }
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let spec = BoundaryLayer2D {
            interior_mesh_size: Some(1.0),
            ..BoundaryLayer2D::new(1, 0.1, 1.0)
        };
        let out = insert_boundary_layers(&mesh, &poly, &spec).expect("layers insert");

        // Mid-bottom boundary node (0.5, 0) grows to (0.5, +0.1).
        assert!(has_node_at(&out, 0.5, 0.1, EPS));
        // 5 boundary segments x 1 layer, no collision.
        assert_eq!(count_type(&out, ElementType::Quad4), 5);
        assert_all_positive_area(&out);
        assert!((total_area(&out) - 1.0).abs() < EPS, "area {}", total_area(&out));
    }

    #[test]
    fn geometric_thickness_progression_and_quad_rows_pinned() {
        let (poly, mesh) = rect_mesh(4.0, 1.0, 0.4);
        // Total depth 0.21 stays below the corner-notch crossing threshold,
        // so no column clamps and every ring row sits at the exact spec
        // depth (the corner-clamp path is exercised by the dedicated test
        // below).
        let spec = BoundaryLayer2D {
            interior_mesh_size: Some(0.3),
            ..BoundaryLayer2D::new(3, 0.03, 2.0)
        };
        let out = insert_boundary_layers(&mesh, &poly, &spec).expect("layers insert");

        // Bottom-edge columns advance 0.03, 0.09, 0.21 (cumulative).
        for k in 1..=9usize {
            let x = 0.4 * k as f64;
            assert!(has_node_at(&out, x, 0.03, EPS), "missing ({x}, 0.03)");
            assert!(has_node_at(&out, x, 0.09, EPS), "missing ({x}, 0.09)");
            assert!(has_node_at(&out, x, 0.21, EPS), "missing ({x}, 0.21)");
        }

        // One quad per boundary segment per layer, no clamping: the segment
        // count is measured from the input mesh's 1-incidence edges.
        let mut census: HashMap<[u64; 2], usize> = HashMap::new();
        for e in &mesh.elements {
            for k in 0..3 {
                let u = e.node_ids[k];
                let v = e.node_ids[(k + 1) % 3];
                *census.entry([u.min(v), u.max(v)]).or_insert(0) += 1;
            }
        }
        let m = census.values().filter(|&&c| c == 1).count();
        assert_eq!(count_type(&out, ElementType::Quad4), m * 3);
        assert_all_positive_area(&out);
        // Area coverage: layers + interior == polygon area.
        assert!(
            (total_area(&out) - 4.0).abs() < EPS,
            "total area {} != 4.0",
            total_area(&out)
        );
        // Interior is still triangulated (no holes).
        assert!(count_type(&out, ElementType::Triangle3) >= 1);
    }

    #[test]
    fn corner_notch_self_crossing_clamps_locally() {
        // 3 layers of 0.05/0.10/0.20 (total 0.35) on the 4x1 rectangle: the
        // averaged-normal corner columns notch outward enough that ring 3
        // self-crosses near the corners. The reactive clamp must stay LOCAL:
        // straight columns keep the full spec depth (rows at 0.05/0.15/0.35
        // pinned), corner columns clamp below the collision, and the mesh
        // stays valid with exact area coverage.
        let (poly, mesh) = rect_mesh(4.0, 1.0, 0.4);
        let spec = BoundaryLayer2D {
            interior_mesh_size: Some(0.3),
            ..BoundaryLayer2D::new(3, 0.05, 2.0)
        };
        let out = insert_boundary_layers(&mesh, &poly, &spec).expect("layers insert");

        assert_all_positive_area(&out);
        assert!((total_area(&out) - 4.0).abs() < EPS, "area {}", total_area(&out));
        assert!(has_node_at(&out, 0.4, 0.05, EPS), "row 1 not at spec depth");
        assert!(has_node_at(&out, 0.4, 0.15, EPS), "row 2 not at spec depth");
        assert!(has_node_at(&out, 0.4, 0.35, EPS), "row 3 not at spec depth");
        // no node is pushed outside the domain by the clamping
        for nd in out.nodes.values() {
            let p = [nd.position.x, nd.position.y];
            assert!(
                poly.contains(p) || near_boundary(&poly, p, 1e-9),
                "node ({}, {}) left the domain",
                p[0],
                p[1]
            );
        }
    }

    // ── Spec validation / no-op ──────────────────────────────────────────────

    #[test]
    fn zero_layers_is_a_passthrough() {
        let (poly, mesh) = rect_mesh(1.0, 1.0, 0.5);
        let out = insert_boundary_layers(&mesh, &poly, &BoundaryLayer2D::new(0, 0.1, 1.0))
            .expect("zero layers is ok");
        assert_eq!(out.node_count(), mesh.node_count());
        assert_eq!(out.element_count(), mesh.element_count());
    }

    #[test]
    fn invalid_specs_are_rejected() {
        let (poly, mesh) = rect_mesh(1.0, 1.0, 0.5);
        let bad = BoundaryLayer2D::new(2, 0.0, 1.0);
        assert!(matches!(
            insert_boundary_layers(&mesh, &poly, &bad),
            Err(BoundaryLayerError::InvalidSpec(_))
        ));
        let bad = BoundaryLayer2D::new(2, 0.1, 0.0);
        assert!(matches!(
            insert_boundary_layers(&mesh, &poly, &bad),
            Err(BoundaryLayerError::InvalidSpec(_))
        ));
        let mut bad = BoundaryLayer2D::new(2, 0.1, 1.0);
        bad.clamp_safety = 0.9; // must be <= 0.5
        assert!(matches!(
            insert_boundary_layers(&mesh, &poly, &bad),
            Err(BoundaryLayerError::InvalidSpec(_))
        ));
    }

    #[test]
    fn mesh_not_on_polygon_is_rejected() {
        let (_, mesh) = rect_mesh(1.0, 1.0, 0.5);
        let elsewhere = Polygon2D::new(vec![[10.0, 10.0], [11.0, 10.0], [11.0, 11.0]]);
        assert!(matches!(
            insert_boundary_layers(&mesh, &elsewhere, &BoundaryLayer2D::new(1, 0.1, 1.0)),
            Err(BoundaryLayerError::InvalidInput(_))
        ));
    }

    // ── Round trip through rmsh-io (write MSH v2 + parse) ────────────────────

    #[test]
    fn layered_mesh_round_trips_through_msh_v2() {
        let (poly, mesh) = rect_mesh(4.0, 1.0, 0.4);
        let spec = BoundaryLayer2D {
            interior_mesh_size: Some(0.3),
            ..BoundaryLayer2D::new(3, 0.05, 2.0)
        };
        let out = insert_boundary_layers(&mesh, &poly, &spec).expect("layers insert");

        let mut bytes = Vec::new();
        rmsh_io::write_msh_v2(&mut bytes, &out).expect("msh v2 write");
        let back = rmsh_io::parse_msh(bytes.as_slice()).expect("msh v2 parse");

        assert_eq!(back.node_count(), out.node_count());
        assert_eq!(back.element_count(), out.element_count());
        assert_eq!(count_type(&back, ElementType::Quad4), count_type(&out, ElementType::Quad4));
        assert_eq!(
            count_type(&back, ElementType::Triangle3),
            count_type(&out, ElementType::Triangle3)
        );
    }
}
