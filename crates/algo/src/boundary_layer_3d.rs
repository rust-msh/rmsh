//! 3-D boundary layers (ROADMAP P4): structured prism layers inserted from a
//! closed triangle surface mesh, growing inward, with per-vertex collision
//! clamping for thin features.
//!
//! # Algorithm provenance
//!
//! The algorithm follows Gmsh's 3-D boundary-layer machinery
//! (`src/Mesh/meshGRegionBoundaryLayer.cpp`, building on the 2-D concepts of
//! `meshGFaceBoundaryLayer.cpp`, GPL-2+). It was understood from those sources
//! and re-implemented from scratch in Rust for rmsh (MIT); no code text was
//! copied. What is carried over is the *architecture*:
//!
//! - the closed surface mesh's triangles become ring 0 of the layer structure
//!   (Gmsh: the boundary faces of the region);
//! - each surface vertex carries an offset direction derived from the
//!   incident triangle normals (Gmsh `blyr_mvertex::average_normal`), and the
//!   layer elements are prisms between consecutive rings;
//! - the remaining interior cavity after the final ring is filled by a
//!   separate volume mesher (Gmsh: the region interior is remeshed with the
//!   layer front as its boundary).
//!
//! # v1 scope (honest)
//!
//! - One averaged normal per vertex: Gmsh's richer treatment — per-face
//!   normals, ridge classification (INTERNAL / EXTERNAL / FLAT), fan rows of
//!   sub-normals across external ridges, and per-corner retriangulation — is
//!   out of v1 scope. The single area-weighted vertex normal reproduces Gmsh's
//!   "one normal per vertex" baseline (`_threshold_angle = M_PI / 6` with
//!   `N_SUBNORMALS = 0`).
//! - The input must be a CLOSED, manifold, consistently oriented triangle
//!   surface (watertight). v1 checks the cheap structural consequences (every
//!   edge shared by exactly two triangles, non-zero enclosed volume) but the
//!   full watertight precheck / repair-hint card is a separate P4 item.
//! - Collision handling is the *preventive ray-cap* only (the reactive
//!   ring-self-intersection clamp of the 2-D module has no cheap 3-D analog
//!   and is deferred): each vertex's TOTAL advance is capped at
//!   `clamp_safety * clearance(vertex)` where `clearance` is the distance
//!   along the vertex's inward ray to the nearest non-incident surface
//!   triangle. This inherits the 2-D module's ray-vs-point-distance lesson:
//!   the ray-cast is the honest metric — a straight-line clearance, not a
//!   nearest-point distance, because the column advances along the ray. With
//!   `clamp_safety <= 0.5` two opposing walls cannot cross (the cap is at
//!   most half of each column's own distance to the opposing wall).
//! - The clamp is UNIFORM per vertex across its whole column: ring `k`'s
//!   advance is `min(cumulative_height(first_height, growth, k), cap)`, so a
//!   clamped column keeps its shape and simply stops earlier; there is no
//!   per-layer twist / anisotropic shortening in v1. Stalled columns (advance
//!   unchanged between consecutive rings) reuse the previous ring's node id,
//!   and a triangle whose three columns all stall produces no prism for that
//!   transition (the exact 3-D analog of the 2-D module's skipped flat quad).
//! - The tetrahedralizer used for the interior cavity
//!   ([`tetrahedralize_closed_surface`]) is a centroid-star mesher: it works
//!   best for convex / star-shaped cavities, which is the fixture class of
//!   this card. Non-star cavities are caught by the volume-conservation
//!   guard and reported, not silently mis-meshed.
//! - No thickness field: unlike [`crate::boundary_layer_2d::BoundaryLayer2D`],
//!   the v1 3-D spec has a constant `first_height` only. A per-node `hWall` /
//!   SizeField profile is deferred until the 2-D field semantics settle.
//!
//! # Orientation convention (pinned)
//!
//! The input surface's orientation is CANONICALIZED to outward: the enclosed
//! signed volume (divergence theorem) is computed; if it is negative every
//! triangle's winding is flipped for the rest of the pipeline, so the output
//! is identical for an inward- or outward-oriented input. A near-zero
//! enclosed volume is rejected (not a closed solid).
//!
//! With outward orientation, the inward offset of a vertex is `-normal`, and
//! the prism for surface triangle `(a, b, c)` between ring `k` (surface side)
//! and ring `k + 1` (deeper) is emitted as
//! `[a_{k+1}, b_{k+1}, c_{k+1}, a_k, b_k, c_k]`
//! (Gmsh Prism6 node order: nodes 0-2 the deeper triangle, nodes 3-5 the
//! shallower one, node 3 above node 0 and so on). Its signed volume
//! (3-tetrahedron decomposition, see [`prism_volume6`]) is positive exactly
//! when the offset advances against the bottom triangle's normal — i.e. all
//! layer prisms of an outward-oriented closed surface come out with POSITIVE
//! signed volume. The final offset surface is emitted with the same
//! (canonical, outward) winding as the input triangles, so it is the
//! outward-oriented boundary of the remaining cavity and its divergence-theorem
//! volume is positive.

use std::collections::{HashMap, HashSet};

use rmsh_model::{Element, ElementType, Mesh, Node};
use thiserror::Error;

use crate::boundary_layer_2d::cumulative_height;
use crate::tetrahedralize3d::tetrahedralize_closed_surface;

/// Everything that can go wrong while inserting 3-D boundary layers.
#[derive(Error, Debug, Clone, PartialEq)]
pub enum BoundaryLayer3DError {
    /// The layer specification itself is invalid.
    #[error("invalid boundary layer 3D spec: {0}")]
    InvalidSpec(String),
    /// The input surface mesh is unusable.
    #[error("invalid boundary layer 3D input: {0}")]
    InvalidInput(String),
    /// Layer construction failed (degenerate prisms, volume conservation
    /// violated, interior tetrahedralization failed, ...).
    #[error("boundary layer 3D generation failed: {0}")]
    Generation(String),
}

/// Specification of a 3-D boundary layer.
///
/// Layer *k* (1-based) has thickness `first_height * growth^(k-1)`: a
/// geometric progression (the same convention as the 2-D module). The
/// cumulative advance after `k` layers is
/// [`cumulative_height`](crate::boundary_layer_2d::cumulative_height).
#[derive(Debug, Clone)]
pub struct BoundaryLayer3D {
    /// Number of layers to insert (0 = no-op passthrough of the surface).
    pub layers: usize,
    /// Height of the first (surface-contacting) layer.
    pub first_height: f64,
    /// Growth ratio between consecutive layer thicknesses. Must be finite
    /// and positive; 1.0 gives uniform layers.
    pub growth: f64,
    /// Fraction of the local ray clearance a vertex may advance at most.
    /// Must lie in `(0, 0.5]`: with 0.5, two opposing walls exactly touch;
    /// the default 0.45 keeps a 10% margin.
    pub clamp_safety: f64,
}

impl BoundaryLayer3D {
    /// A spec with default collision parameters (`clamp_safety = 0.45`).
    pub fn new(layers: usize, first_height: f64, growth: f64) -> Self {
        Self {
            layers,
            first_height,
            growth,
            clamp_safety: 0.45,
        }
    }

    fn validate(&self) -> Result<(), BoundaryLayer3DError> {
        if !self.growth.is_finite() || self.growth <= 0.0 {
            return Err(BoundaryLayer3DError::InvalidSpec(format!(
                "growth must be finite and positive, got {}",
                self.growth
            )));
        }
        if !self.first_height.is_finite() || self.first_height <= 0.0 {
            return Err(BoundaryLayer3DError::InvalidSpec(format!(
                "first_height must be finite and positive, got {}",
                self.first_height
            )));
        }
        if !self.clamp_safety.is_finite()
            || self.clamp_safety <= 0.0
            || self.clamp_safety > 0.5
        {
            return Err(BoundaryLayer3DError::InvalidSpec(format!(
                "clamp_safety must be finite in (0, 0.5], got {}",
                self.clamp_safety
            )));
        }
        Ok(())
    }
}

// ─── Small 3-D geometry helpers ──────────────────────────────────────────────

type Vec3 = [f64; 3];

fn vsub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vadd(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn vscale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn vdot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vcross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn vnorm(a: Vec3) -> f64 {
    vdot(a, a).sqrt()
}

/// Six times the signed volume of tetrahedron `(a, b, c, d)`; positive when
/// `d` lies on the side of the right-hand normal of `(a, b, c)`.
fn tet_volume6(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> f64 {
    vdot(vcross(vsub(b, a), vsub(c, a)), vsub(d, a))
}

/// Signed volume of the prism `(p0..p5)` (Gmsh Prism6 order: `p0p1p2` the
/// bottom triangle, `p3p4p5` the top one, `p3` above `p0` and so on), via the
/// standard 3-tetrahedron decomposition. Positive when the bottom triangle's
/// right-hand normal points from the bottom toward the top.
fn prism_volume6(p: &[Vec3; 6]) -> f64 {
    (tet_volume6(p[0], p[1], p[2], p[5])
        + tet_volume6(p[0], p[1], p[5], p[4])
        + tet_volume6(p[0], p[4], p[5], p[3]))
        / 6.0
}

/// Enclosed volume of a closed triangle surface (divergence theorem):
/// `V = (1/6) * sum dot(a, cross(b, c))` over all triangles. Positive for an
/// outward-oriented (counter-clockwise seen from outside) surface.
fn closed_surface_volume(tris: &[[u64; 3]], pos: &HashMap<u64, Vec3>) -> f64 {
    let mut v = 0.0;
    for t in tris {
        let a = pos[&t[0]];
        let b = pos[&t[1]];
        let c = pos[&t[2]];
        v += vdot(a, vcross(b, c));
    }
    v / 6.0
}

/// First hit distance `t > 0` of the ray `orig + t * dir` (`dir` unit) against
/// the triangle `(a, b, c)` (two-sided), or infinity when the ray misses.
/// Moeller-Trumbore with a small barycentric slack so that vertex/edge grazes
/// count as hits (deterministic for the fixture class; a missed graze only
/// loosens the cap, it cannot invert an element).
fn ray_triangle_hit(orig: Vec3, dir: Vec3, a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let e1 = vsub(b, a);
    let e2 = vsub(c, a);
    let p = vcross(dir, e2);
    let det = vdot(e1, p);
    // parallel: threshold scales with the triangle area scale
    if det.abs() <= 1e-14 * (vnorm(e1) * vnorm(e2)) {
        return f64::INFINITY;
    }
    let inv = 1.0 / det;
    let tvec = vsub(orig, a);
    let u = vdot(tvec, p) * inv;
    if u < -1e-9 || u > 1.0 + 1e-9 {
        return f64::INFINITY;
    }
    let q = vcross(tvec, e1);
    let v = vdot(dir, q) * inv;
    if v < -1e-9 || u + v > 1.0 + 1e-9 {
        return f64::INFINITY;
    }
    let t = vdot(e2, q) * inv;
    if t > 1e-12 {
        t
    } else {
        f64::INFINITY
    }
}

// ─── Input surface plumbing ──────────────────────────────────────────────────

/// Collect the triangle soup of a closed surface mesh: triangle connectivity
/// (canonicalized to outward winding), node positions, the sorted surface
/// vertex ids and per-vertex area-weighted outward normals.
struct SurfaceData {
    /// Triangles in canonical (outward) winding.
    tris: Vec<[u64; 3]>,
    /// Positions of all surface nodes.
    pos: HashMap<u64, Vec3>,
    /// Distinct surface vertex ids, ascending (deterministic order).
    verts: Vec<u64>,
    /// Unit outward normal per surface vertex (area-weighted average of the
    /// incident triangle normals: the sum of the unnormalized cross products,
    /// which is exactly area weighting since `|cross(b-a, c-a)| = 2 * area`).
    normals: HashMap<u64, Vec3>,
    /// Bounding-box diagonal (epsilon scale).
    scale: f64,
    /// Enclosed volume of the canonical (outward) surface.
    volume: f64,
}

fn node_pos(mesh: &Mesh, id: u64) -> Result<Vec3, BoundaryLayer3DError> {
    mesh.nodes
        .get(&id)
        .map(|n| [n.position.x, n.position.y, n.position.z])
        .ok_or_else(|| BoundaryLayer3DError::InvalidInput(format!("node {id} missing from mesh")))
}

fn collect_surface(mesh: &Mesh) -> Result<SurfaceData, BoundaryLayer3DError> {
    // ── Triangle soup: pure triangle surface, Line/Point decorations ignored ──
    let mut tris: Vec<[u64; 3]> = Vec::new();
    for e in &mesh.elements {
        match e.etype {
            ElementType::Triangle3 => tris.push([e.node_ids[0], e.node_ids[1], e.node_ids[2]]),
            ElementType::Line2 | ElementType::Point1 => {}
            other => {
                return Err(BoundaryLayer3DError::InvalidInput(format!(
                    "input mesh must be a pure triangle surface, found {other:?}"
                )));
            }
        }
    }
    if tris.is_empty() {
        return Err(BoundaryLayer3DError::InvalidInput(
            "input mesh has no triangles".to_string(),
        ));
    }

    // ── Closed-manifold structural check: every undirected edge exactly twice ──
    let mut edge_count: HashMap<[u64; 2], usize> = HashMap::new();
    let mut pos_set: HashSet<u64> = HashSet::new();
    for t in &tris {
        for k in 0..3 {
            let u = t[k];
            let v = t[(k + 1) % 3];
            if u == v {
                return Err(BoundaryLayer3DError::InvalidInput(format!(
                    "degenerate triangle edge ({u}, {v})"
                )));
            }
            *edge_count.entry([u.min(v), u.max(v)]).or_insert(0) += 1;
            pos_set.insert(u);
        }
    }
    for (key, count) in &edge_count {
        if *count != 2 {
            return Err(BoundaryLayer3DError::InvalidInput(format!(
                "surface edge ({}, {}) is used {} times (not a closed manifold)",
                key[0], key[1], count
            )));
        }
    }

    // ── Node positions + epsilon scale ────────────────────────────────────────
    let mut pos: HashMap<u64, Vec3> = HashMap::with_capacity(pos_set.len());
    for &id in &pos_set {
        pos.insert(id, node_pos(mesh, id)?);
    }
    let mut mn = [f64::MAX; 3];
    let mut mx = [f64::MIN; 3];
    for p in pos.values() {
        for k in 0..3 {
            mn[k] = mn[k].min(p[k]);
            mx[k] = mx[k].max(p[k]);
        }
    }
    let scale = vnorm([mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]]).max(1e-30);

    // ── Orientation canonicalization: flip everything if the enclosed volume
    // is negative; reject a near-zero enclosed volume ─────────────────────────
    let volume = closed_surface_volume(&tris, &pos);
    let vol_eps = 1e-9 * scale * scale * scale;
    if volume.abs() <= vol_eps {
        return Err(BoundaryLayer3DError::InvalidInput(format!(
            "surface encloses no volume ({volume}; not a closed solid)"
        )));
    }
    if volume < 0.0 {
        for t in &mut tris {
            t.swap(1, 2);
        }
    }
    let volume = volume.abs();

    // ── Sorted vertex list ────────────────────────────────────────────────────
    let mut verts: Vec<u64> = pos_set.into_iter().collect();
    verts.sort_unstable();

    // ── Area-weighted vertex normals ──────────────────────────────────────────
    let mut sums: HashMap<u64, Vec3> = HashMap::with_capacity(verts.len());
    for v in &verts {
        sums.insert(*v, [0.0; 3]);
    }
    for t in &tris {
        let a = pos[&t[0]];
        let b = pos[&t[1]];
        let c = pos[&t[2]];
        let n = vcross(vsub(b, a), vsub(c, a)); // 2 * area * unit normal
        for &id in t {
            let s = sums.get_mut(&id).expect("vertex registered above");
            *s = vadd(*s, n);
        }
    }
    let mut normals: HashMap<u64, Vec3> = HashMap::with_capacity(verts.len());
    for v in &verts {
        let s = sums[v];
        let len = vnorm(s);
        if !len.is_finite() || len <= 1e-12 * scale {
            return Err(BoundaryLayer3DError::InvalidInput(format!(
                "vertex {v} has a degenerate area-weighted normal (norm {len})"
            )));
        }
        normals.insert(*v, vscale(s, 1.0 / len));
    }

    Ok(SurfaceData {
        tris,
        pos,
        verts,
        normals,
        scale,
        volume,
    })
}

/// Insert boundary layers into a closed triangle surface mesh (see module
/// docs): prism columns growing inward from every surface triangle with
/// per-vertex collision caps, then the remaining interior cavity filled by
/// [`tetrahedralize_closed_surface`] over the final offset surface.
pub fn insert_boundary_layers_3d(
    surface: &Mesh,
    spec: &BoundaryLayer3D,
) -> Result<Mesh, BoundaryLayer3DError> {
    spec.validate()?;
    if spec.layers == 0 {
        return Ok(surface.clone());
    }

    // ── Ring 0: the closed surface, canonicalized to outward orientation ──
    let sd = collect_surface(surface)?;
    let n = sd.verts.len();
    let vidx: HashMap<u64, usize> = sd
        .verts
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i))
        .collect();

    // ── Preventive collision caps: per-vertex ray-cast along the inward ray ──
    // Clearance of vertex v = distance along its inward ray (-normal) to the
    // nearest NON-INCIDENT surface triangle (the first thing the advancing
    // column would hit; incident triangles are exempt — the same adjacency
    // rule as the 2-D module). Cap = clamp_safety * clearance: an honest
    // along-the-ray metric, so with clamp_safety <= 0.5 two opposing walls
    // cannot cross. Rays that miss everything give an infinite cap.
    let mut caps: Vec<f64> = Vec::with_capacity(n);
    for v in &sd.verts {
        let orig = sd.pos[v];
        let dir = vscale(sd.normals[v], -1.0); // inward
        let mut clearance = f64::INFINITY;
        for t in &sd.tris {
            if t[0] == *v || t[1] == *v || t[2] == *v {
                continue; // incident-triangle exemption
            }
            let hit = ray_triangle_hit(orig, dir, sd.pos[&t[0]], sd.pos[&t[1]], sd.pos[&t[2]]);
            if hit < clearance {
                clearance = hit;
            }
        }
        caps.push(clearance * spec.clamp_safety);
    }

    // ── Ring advances: min(cumulative, cap) — a UNIFORM clamp per column, so
    // ring k's advance stays monotone in k and clamped columns keep shape ────
    let mut advances: Vec<Vec<f64>> = Vec::with_capacity(spec.layers + 1);
    advances.push(vec![0.0; n]);
    for k in 1..=spec.layers {
        let row: Vec<f64> = (0..n)
            .map(|i| cumulative_height(spec.first_height, spec.growth, k).min(caps[i]))
            .collect();
        advances.push(row);
    }

    // ── Ring node allocation: stalled columns (advance unchanged between
    // consecutive rings) reuse the previous ring's node id ───────────────────
    let mut all_pos: HashMap<u64, Vec3> = sd.pos.clone();
    let mut new_id = sd.verts[n - 1] + 1; // surface ids are the largest input ids here
    let mut ring_ids: Vec<Vec<u64>> = Vec::with_capacity(spec.layers + 1);
    ring_ids.push(sd.verts.clone());
    for k in 1..=spec.layers {
        let mut row = Vec::with_capacity(n);
        for i in 0..n {
            if advances[k][i] == advances[k - 1][i] {
                row.push(ring_ids[k - 1][i]);
            } else {
                let p = vadd(
                    sd.pos[&sd.verts[i]],
                    vscale(sd.normals[&sd.verts[i]], -advances[k][i]),
                );
                all_pos.insert(new_id, p);
                row.push(new_id);
                new_id += 1;
            }
        }
        ring_ids.push(row);
    }

    // ── Output mesh: nodes first (surface ids, then ring allocation order) ──
    let mut out = Mesh::new();
    for v in &sd.verts {
        let p = sd.pos[v];
        out.add_node(Node::new(*v, p[0], p[1], p[2]));
    }
    for k in 1..=spec.layers {
        for i in 0..n {
            let id = ring_ids[k][i];
            if ring_ids[k - 1][i] != id {
                let p = all_pos[&id];
                out.add_node(Node::new(id, p[0], p[1], p[2]));
            }
        }
    }

    let vol_eps = 1e-12 * sd.scale * sd.scale * sd.scale;
    let mut eid: u64 = 1;

    // ── Prism columns: one prism per surface triangle per layer transition ──
    // Node order [a_{k}, b_{k}, c_{k}, a_{k-1}, b_{k-1}, c_{k-1}]: the DEEPER
    // ring first, so the bottom triangle's outward normal points from bottom
    // to top and the signed volume is positive (see module docs). A triangle
    // whose three columns all stall would duplicate ids — skipped, the exact
    // analog of the 2-D module's skipped flat quad.
    let mut prisms_volume = 0.0;
    for k in 1..=spec.layers {
        for t in &sd.tris {
            let (ia, ib, ic) = (vidx[&t[0]], vidx[&t[1]], vidx[&t[2]]);
            let stalled = advances[k][ia] == advances[k - 1][ia]
                && advances[k][ib] == advances[k - 1][ib]
                && advances[k][ic] == advances[k - 1][ic];
            if stalled {
                continue;
            }
            let ids = [
                ring_ids[k][ia],
                ring_ids[k][ib],
                ring_ids[k][ic],
                ring_ids[k - 1][ia],
                ring_ids[k - 1][ib],
                ring_ids[k - 1][ic],
            ];
            let p: [Vec3; 6] = ids.map(|id| all_pos[&id]);
            let vol = prism_volume6(&p);
            if !(vol > vol_eps) {
                return Err(BoundaryLayer3DError::Generation(format!(
                    "degenerate or inverted prism {ids:?} (volume {vol})"
                )));
            }
            prisms_volume += vol;
            out.add_element(Element::new(eid, ElementType::Prism6, ids.to_vec()));
            eid += 1;
        }
    }

    // ── Final offset surface: the cavity boundary, canonical outward winding,
    // outward-oriented w.r.t. the remaining interior ─────────────────────────
    let last = spec.layers;
    for t in &sd.tris {
        let ids = [
            ring_ids[last][vidx[&t[0]]],
            ring_ids[last][vidx[&t[1]]],
            ring_ids[last][vidx[&t[2]]],
        ];
        out.add_element(Element::new(eid, ElementType::Triangle3, ids.to_vec()));
        eid += 1;
    }

    // ── Interior cavity: centroid-star tetrahedralization of the final ring ──
    let mut core = Mesh::new();
    let mut core_nodes: HashSet<u64> = HashSet::new();
    let mut core_tris: Vec<[u64; 3]> = Vec::with_capacity(sd.tris.len());
    for t in &sd.tris {
        let ids = [
            ring_ids[last][vidx[&t[0]]],
            ring_ids[last][vidx[&t[1]]],
            ring_ids[last][vidx[&t[2]]],
        ];
        core_nodes.extend(ids);
        core_tris.push(ids);
    }
    let mut core_node_ids: Vec<u64> = core_nodes.into_iter().collect();
    core_node_ids.sort_unstable();
    for id in &core_node_ids {
        let p = all_pos[id];
        core.add_node(Node::new(*id, p[0], p[1], p[2]));
    }
    for (i, ids) in core_tris.iter().enumerate() {
        core.add_element(Element::new(i as u64 + 1, ElementType::Triangle3, ids.to_vec()));
    }
    let tets = tetrahedralize_closed_surface(&core).map_err(|e| {
        BoundaryLayer3DError::Generation(format!("interior tetrahedralization failed: {e}"))
    })?;

    let mut tets_volume = 0.0;
    for node in tets.nodes.values() {
        if !all_pos.contains_key(&node.id) {
            all_pos.insert(
                node.id,
                [node.position.x, node.position.y, node.position.z],
            );
            out.add_node(node.clone());
        }
    }
    for e in &tets.elements {
        let p: Vec<Vec3> = e
            .node_ids
            .iter()
            .map(|id| {
                all_pos.get(id).copied().ok_or_else(|| {
                    BoundaryLayer3DError::Generation(format!(
                        "tetrahedralizer referenced unknown node {id}"
                    ))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        tets_volume += tet_volume6(p[0], p[1], p[2], p[3]).abs() / 6.0;
        out.add_element(Element::new(eid, ElementType::Tetrahedron4, e.node_ids.clone()));
        eid += 1;
    }

    // ── Volume conservation guards: prisms + tets must reproduce the input's
    // enclosed volume, and the tets must exactly fill the cavity the final
    // surface encloses (a non-star cavity mis-meshes and trips this) ────────
    let vol_tol = 1e-6 * sd.scale * sd.scale * sd.scale;
    if (prisms_volume + tets_volume - sd.volume).abs() > vol_tol {
        return Err(BoundaryLayer3DError::Generation(format!(
            "volume conservation violated: prisms {prisms_volume} + tets {tets_volume} != enclosed {}",
            sd.volume
        )));
    }
    let cavity_volume = closed_surface_volume(&core_tris, &all_pos);
    if (cavity_volume - tets_volume).abs() > vol_tol {
        return Err(BoundaryLayer3DError::Generation(format!(
            "interior tets ({tets_volume}) do not fill the final-surface cavity ({cavity_volume}); \
             cavity is likely not star-shaped"
        )));
    }

    Ok(out)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rmsh_io::{parse_msh, write_msh_v2};

    const EPS: f64 = 1e-9;

    trait AssertF64 {
        fn assert_close(self, other: f64);
    }
    impl AssertF64 for f64 {
        fn assert_close(self, other: f64) {
            assert!(
                (self - other).abs() < 1e-12,
                "{self} != {other}"
            );
        }
    }

    /// Closed box surface `[0, dx] x [0, dy] x [0, dz]`: 8 vertices, 12
    /// outward-oriented triangles (each quad face split by one diagonal).
    /// The diagonals are chosen to join each face's two even-parity corners
    /// (parity of x+Y+z), so every corner sees an EQUAL weight from each of
    /// its three incident faces and the area-weighted normal is the exact
    /// body diagonal (+-1,+-1,+-1)/sqrt(3) for all 8 corners.
    fn box_surface(dx: f64, dy: f64, dz: f64) -> Mesh {
        let mut m = Mesh::new();
        let pts: [[f64; 3]; 8] = [
            [0.0, 0.0, 0.0], // 1 (even)
            [dx, 0.0, 0.0],  // 2
            [dx, dy, 0.0],   // 3 (even)
            [0.0, dy, 0.0],  // 4
            [0.0, 0.0, dz],  // 5
            [dx, 0.0, dz],   // 6 (even)
            [dx, dy, dz],    // 7
            [0.0, dy, dz],   // 8 (even)
        ];
        for (i, p) in pts.iter().enumerate() {
            m.add_node(Node::new(i as u64 + 1, p[0], p[1], p[2]));
        }
        // outward windings: bottom (0,0,-1), top (0,0,+1), y=0 (0,-1,0),
        // x=dx (+1,0,0), y=dy (0,+1,0), x=0 (-1,0,0); even-parity diagonals
        let tris: [[u64; 3]; 12] = [
            [1, 3, 2],
            [1, 4, 3], // bottom, diagonal (1,3)
            [5, 6, 8],
            [6, 7, 8], // top, diagonal (6,8)
            [1, 2, 6],
            [1, 6, 5], // front y=0, diagonal (1,6)
            [2, 3, 6],
            [3, 7, 6], // right x=dx, diagonal (3,6)
            [3, 4, 8],
            [3, 8, 7], // back y=dy, diagonal (3,8)
            [4, 1, 8],
            [1, 5, 8], // left x=0, diagonal (1,8)
        ];
        for (i, t) in tris.iter().enumerate() {
            m.add_element(Element::new(i as u64 + 1, ElementType::Triangle3, t.to_vec()));
        }
        m
    }

    fn count_type(mesh: &Mesh, t: ElementType) -> usize {
        mesh.elements.iter().filter(|e| e.etype == t).count()
    }

    /// Is there a node at position `p` within `eps`?
    fn has_node_at(mesh: &Mesh, x: f64, y: f64, z: f64, eps: f64) -> bool {
        mesh.nodes.values().any(|n| {
            (n.position.x - x).abs() <= eps
                && (n.position.y - y).abs() <= eps
                && (n.position.z - z).abs() <= eps
        })
    }

    fn pos_of(mesh: &Mesh, id: u64) -> Vec3 {
        let n = mesh.nodes.get(&id).expect("node exists");
        [n.position.x, n.position.y, n.position.z]
    }

    /// Signed volume of one output prism element.
    fn prism_volume_of(mesh: &Mesh, e: &Element) -> f64 {
        let p: Vec<Vec3> = e.node_ids.iter().map(|&id| pos_of(mesh, id)).collect();
        prism_volume6(p.as_slice().try_into().expect("6 nodes"))
    }

    /// Signed volume of one output tet element (absolute: the centroid-star
    /// tetrahedralizer does not guarantee connectivity sign).
    fn tet_volume_of(mesh: &Mesh, e: &Element) -> f64 {
        let p: Vec<Vec3> = e.node_ids.iter().map(|&id| pos_of(mesh, id)).collect();
        tet_volume6(p[0], p[1], p[2], p[3]).abs() / 6.0
    }

    /// Unsigned area2 of one output triangle element.
    fn triangle_area2_of(mesh: &Mesh, e: &Element) -> f64 {
        let p: Vec<Vec3> = e.node_ids.iter().map(|&id| pos_of(mesh, id)).collect();
        vnorm(vcross(vsub(p[1], p[0]), vsub(p[2], p[0])))
    }

    /// Volumes summed by element family over an output mesh.
    fn volume_split(mesh: &Mesh) -> (f64, f64) {
        let mut prisms = 0.0;
        let mut tets = 0.0;
        for e in &mesh.elements {
            match e.etype {
                ElementType::Prism6 => prisms += prism_volume_of(mesh, e),
                ElementType::Tetrahedron4 => tets += tet_volume_of(mesh, e),
                _ => {}
            }
        }
        (prisms, tets)
    }

    /// Signed volume of the final offset surface: the sum over the mesh's
    /// Triangle3 elements (they form the closed cavity boundary).
    fn triangle_surface_volume(mesh: &Mesh) -> f64 {
        let mut v = 0.0;
        for e in &mesh.elements {
            if e.etype != ElementType::Triangle3 {
                continue;
            }
            let p: Vec<Vec3> = e.node_ids.iter().map(|&id| pos_of(mesh, id)).collect();
            v += vdot(p[0], vcross(p[1], p[2]));
        }
        v / 6.0
    }

    /// All volume elements positive (prisms: signed per the documented
    /// inward convention; tets: non-degenerate), all triangles non-degenerate.
    fn assert_all_positive(mesh: &Mesh) {
        for e in &mesh.elements {
            match e.etype {
                ElementType::Prism6 => {
                    let v = prism_volume_of(mesh, e);
                    assert!(v > 0.0, "non-positive prism {} (vol {v})", e.id);
                }
                ElementType::Tetrahedron4 => {
                    let v = tet_volume_of(mesh, e);
                    assert!(v > 0.0, "degenerate tet {} (vol {v})", e.id);
                }
                ElementType::Triangle3 => {
                    let a2 = triangle_area2_of(mesh, e);
                    assert!(a2 > 0.0, "degenerate triangle {} (area2 {a2})", e.id);
                }
                other => panic!("unexpected element type {other:?}"),
            }
        }
    }

    // ── Pin: the geometric progression helper is reused from the 2-D module ──

    #[test]
    fn cumulative_height_import_pinned() {
        cumulative_height(0.05, 1.2, 3).assert_close(0.05 * (1.0 + 1.2 + 1.44));
        cumulative_height(0.05, 1.2, 0).assert_close(0.0);
    }

    // ── Pin: prism signed-volume convention (module docs "Orientation") ──────

    #[test]
    fn prism_signed_volume_unit_example_pinned() {
        // bottom (0,0,0),(1,0,0),(0,1,0), top translated by +z: positive.
        let p: [Vec3; 6] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        prism_volume6(&p).assert_close(0.5); // base area 1/2, height 1
        // the reversed prism flips BOTH triangle windings (swapping one
        // vertex pair would build a twisted, not reversed, hexahedron)
        let mut q = p;
        q.swap(1, 2);
        q.swap(4, 5);
        prism_volume6(&q).assert_close(-0.5);
    }

    // ── Pin: corner vertex normal is the exact diagonal (unit cube) ──────────

    /// Ring-1 node of corner (0,0,0) sits at first_height * (1,1,1)/sqrt(3)
    /// (the area-weighted normal of the three incident faces), and the whole
    /// progression 0.05 / 0.11 / 0.182 (cumulative) is pinned at 1e-9: the
    /// unit cube's ray clearance sqrt(3) caps at 0.45*sqrt(3) > 0.182, so no
    /// column clamps.
    #[test]
    fn unit_cube_corner_ring_progression_pinned() {
        let cube = box_surface(1.0, 1.0, 1.0);
        let spec = BoundaryLayer3D::new(3, 0.05, 1.2);
        let out = insert_boundary_layers_3d(&cube, &spec).expect("layers insert");

        let s3 = 3.0f64.sqrt();
        for (k, cum) in [0.05, 0.11, 0.182].iter().enumerate() {
            let s = cum / s3;
            // ring ids are allocated after the 8 surface ids, in vertex order:
            // corner (0,0,0) is surface vertex 1 -> its ring-k node id is
            // 8 + 3*(k) + 1 (ring-major, vertex-ascending allocation).
            let id = 8 + (k as u64) * 8 + 1;
            let p = pos_of(&out, id);
            assert!((p[0] - s).abs() < EPS, "ring{} x {} != {}", k + 1, p[0], s);
            assert!((p[1] - s).abs() < EPS, "ring{} y {} != {}", k + 1, p[1], s);
            assert!((p[2] - s).abs() < EPS, "ring{} z {} != {}", k + 1, p[2], s);
        }
        // every ring node strictly inside the unit cube (inward growth)
        for n in out.nodes.values() {
            for c in [n.position.x, n.position.y, n.position.z] {
                assert!(
                    c > -EPS && c < 1.0 + EPS,
                    "node {} at {} left the cube",
                    n.id,
                    c
                );
            }
        }
        assert_eq!(count_type(&out, ElementType::Prism6), 36, "12 tris x 3 layers");
    }

    // ── Pin: counts, positivity, volume conservation (unit cube) ─────────────

    /// 36 prisms + 12 final-surface triangles + 12 interior tets; every
    /// element positive; prisms + tets = 1.0 (the enclosed volume) within
    /// 1e-9; the final offset surface is the shrunk box
    /// (1 - 2*0.182/sqrt(3))^3 and its volume is exactly covered by the tets
    /// (the void check).
    #[test]
    fn unit_cube_counts_positivity_volume_conservation() {
        let cube = box_surface(1.0, 1.0, 1.0);
        let spec = BoundaryLayer3D::new(3, 0.05, 1.2);
        let out = insert_boundary_layers_3d(&cube, &spec).expect("layers insert");

        assert_eq!(count_type(&out, ElementType::Prism6), 36);
        assert_eq!(count_type(&out, ElementType::Triangle3), 12);
        assert_eq!(count_type(&out, ElementType::Tetrahedron4), 12);
        assert_eq!(out.node_count(), 8 + 8 * 3 + 1, "surface + 3 rings + centroid");

        assert_all_positive(&out);

        let (prisms, tets) = volume_split(&out);
        assert!(
            (prisms + tets - 1.0).abs() < EPS,
            "prisms {prisms} + tets {tets} != 1.0"
        );

        let s = 0.182 / 3.0f64.sqrt();
        let cavity = (1.0 - 2.0 * s) * (1.0 - 2.0 * s) * (1.0 - 2.0 * s);
        let surf_vol = triangle_surface_volume(&out);
        assert!(
            (surf_vol - cavity).abs() < EPS,
            "final surface volume {surf_vol} != shrunk box {cavity}"
        );
        assert!(
            (tets - surf_vol).abs() < EPS,
            "tets {tets} do not exactly fill the cavity {surf_vol}"
        );
    }

    // ── Pin: AREA weighting (not unit-normal averaging) of vertex normals ────

    /// Irregular tet (0,0,0),(1,0,0),(0,1,0),(0,0,2): vertex 1's incident
    /// faces have area vectors (0,0,-1), (0,-2,0), (-2,0,0), so the
    /// area-weighted inward direction is (2,2,1)/3 — unit-normal averaging
    /// would give (1,1,1)/sqrt(3). The ring-1 node position discriminates
    /// the two schemes.
    #[test]
    fn area_weighted_normal_pinned_on_irregular_tet() {
        let mut m = Mesh::new();
        let pts: [[f64; 3]; 4] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]];
        for (i, p) in pts.iter().enumerate() {
            m.add_node(Node::new(i as u64 + 1, p[0], p[1], p[2]));
        }
        // outward windings (see the derivation in the test body comment)
        let tris: [[u64; 3]; 4] = [[1, 3, 2], [1, 2, 4], [1, 4, 3], [2, 3, 4]];
        for (i, t) in tris.iter().enumerate() {
            m.add_element(Element::new(i as u64 + 1, ElementType::Triangle3, t.to_vec()));
        }
        let spec = BoundaryLayer3D::new(1, 0.01, 1.0);
        let out = insert_boundary_layers_3d(&m, &spec).expect("layers insert");

        let s = 0.01 / 3.0; // per-axis offsets of (2,2,1)/3 scaled by first_height
        // vertex 1 is the smallest id; its ring-1 node is allocated first
        // (id = 4 surface nodes + 1 = 5)
        let p = pos_of(&out, 5);
        assert!((p[0] - 2.0 * s).abs() < EPS, "x {} != {}", p[0], 2.0 * s);
        assert!((p[1] - 2.0 * s).abs() < EPS, "y {} != {}", p[1], 2.0 * s);
        assert!((p[2] - s).abs() < EPS, "z {} != {}", p[2], s);
        assert_all_positive(&out);
    }

    // ── Pin: collision clamp on a thin slab ──────────────────────────────────

    /// 1 x 1 x 0.1 slab: one closed surface (12 triangles); 3 layers of 0.03
    /// growing 1.2 want a 0.1092 total advance per corner column. The slab
    /// corner normals are dominated by the two 1x1 faces (area-weighted raw
    /// sum (-0.2,-0.2,-2): about 5.7 degrees off -z, NOT a body diagonal —
    /// that pin belongs to the cube only), so the inward ray is nearly
    /// vertical and hits the opposite face after 0.1/n_z ~ 0.101. The cap =
    /// clamp_safety * clearance then clamps ring 2 onward, and the axial
    /// progress at the cap is EXACTLY clamp_safety * gap = 0.045 (the
    /// honest-ray property: t_hit * n_z = gap), so no prism crosses the
    /// mid-plane z = 0.05 from either side. The cavity (and its centroid)
    /// legitimately straddles the mid-plane — only layer elements are checked.
    #[test]
    fn thin_slab_collision_clamped_no_midplane_crossing() {
        let slab = box_surface(1.0, 1.0, 0.1);
        let spec = BoundaryLayer3D::new(3, 0.03, 1.2);
        let out = insert_boundary_layers_3d(&slab, &spec).expect("layers insert");

        // no layer node crosses the mid-plane: every RING node stays on its
        // origin side of z = 0.05 (bottom-grown below 0.046, top-grown above
        // 0.054). Side-wall prisms legitimately span the full slab height
        // (their SURFACE nodes sit at both z = 0 and z = 0.1), and the cavity
        // centroid legitimately sits at the mid-plane — so the check is over
        // ring nodes only: node ids referenced by prisms minus the input
        // surface ids.
        let surface_ids: HashSet<u64> = slab.nodes.keys().copied().collect();
        let mut layer_node_ids: HashSet<u64> = HashSet::new();
        for e in &out.elements {
            if e.etype == ElementType::Prism6 {
                layer_node_ids.extend(e.node_ids.iter().copied());
            }
        }
        for id in layer_node_ids.difference(&surface_ids) {
            let z = pos_of(&out, *id)[2];
            assert!(
                z < 0.046 || z > 0.054,
                "layer node {id} at z = {z} crossed the clamped mid-band"
            );
        }

        // clamped ring-2 corner positions: the cap is reached axially at
        // exactly 0.45 * 0.1 = 0.045, with a small lateral drift
        // 0.0045 = 0.045 * |n_xy / n_z| = 0.045 * 0.1
        assert!(
            has_node_at(&out, 0.0045, 0.0045, 0.045, EPS),
            "bottom clamped corner missing"
        );
        assert!(
            has_node_at(&out, 0.0045, 0.0045, 0.055, EPS),
            "top clamped corner missing"
        );

        // ring 1 is unclamped (0.03 < cap), ring 2 clamps, ring 3 stalls
        // (advance equals ring 2's): 12 tris x 2 moving transitions = 24
        assert_eq!(count_type(&out, ElementType::Prism6), 24);
        assert_eq!(count_type(&out, ElementType::Triangle3), 12);
        assert_eq!(out.node_count(), 8 + 8 + 8 + 1, "surface + 2 rings + centroid");

        assert_all_positive(&out);
        let (prisms, tets) = volume_split(&out);
        assert!(
            (prisms + tets - 0.1).abs() < EPS,
            "prisms {prisms} + tets {tets} != 0.1"
        );
        // cavity: the shrunk box [0.0045, 0.9955]^2 x [0.045, 0.055]
        let cavity = 0.991 * 0.991 * 0.01;
        let surf_vol = triangle_surface_volume(&out);
        assert!(
            (surf_vol - cavity).abs() < EPS,
            "cavity {surf_vol} != {cavity}"
        );
    }

    // ── Pin: stalled columns reuse ring nodes and emit no prisms ─────────────

    /// 4 layers of 0.05 growing 1.5 on the slab: the cap (0.45 * 0.1/n_z ~
    /// 0.0454) is already below ring 1's 0.05, so ring 1 clamps and rings
    /// 2-4 stall: no new nodes, no prisms for the stalled transitions
    /// (12 prisms only), and the final surface sits at the ring-1 positions.
    #[test]
    fn thin_slab_deep_spec_stalls_and_reuses_nodes() {
        let slab = box_surface(1.0, 1.0, 0.1);
        let spec = BoundaryLayer3D::new(4, 0.05, 1.5);
        let out = insert_boundary_layers_3d(&slab, &spec).expect("layers insert");

        assert_eq!(count_type(&out, ElementType::Prism6), 12);
        assert_eq!(count_type(&out, ElementType::Triangle3), 12);
        assert_eq!(out.node_count(), 8 + 8 + 1, "surface + ring 1 (3 stalled) + centroid");

        // the clamped position exists EXACTLY ONCE per corner (node reuse):
        // corner (0,0,0) clamps to (0.0045, 0.0045, 0.045)
        let clamp_hits = out
            .nodes
            .values()
            .filter(|n| {
                (n.position.x - 0.0045).abs() < EPS
                    && (n.position.y - 0.0045).abs() < EPS
                    && (n.position.z - 0.045).abs() < EPS
            })
            .count();
        assert_eq!(clamp_hits, 1, "stalled column must reuse one node");

        assert_all_positive(&out);
        let (prisms, tets) = volume_split(&out);
        assert!(
            (prisms + tets - 0.1).abs() < EPS,
            "prisms {prisms} + tets {tets} != 0.1"
        );
    }

    // ── Pin: inward-oriented input is canonicalized ──────────────────────────

    /// The same cube with every triangle winding flipped (negative enclosed
    /// volume) must produce the identical layered result.
    #[test]
    fn inward_oriented_input_canonicalized() {
        let cube = box_surface(1.0, 1.0, 1.0);
        let mut flipped = Mesh::new();
        for n in cube.nodes.values() {
            flipped.add_node(n.clone());
        }
        for (i, e) in cube.elements.iter().enumerate() {
            flipped.add_element(Element::new(
                i as u64 + 1,
                e.etype,
                vec![e.node_ids[0], e.node_ids[2], e.node_ids[1]],
            ));
        }
        let spec = BoundaryLayer3D::new(3, 0.05, 1.2);
        let out = insert_boundary_layers_3d(&flipped, &spec).expect("layers insert");

        assert_eq!(count_type(&out, ElementType::Prism6), 36);
        assert_eq!(count_type(&out, ElementType::Tetrahedron4), 12);
        // corner (0,0,0) ring-1 node still at +0.05/sqrt(3) per axis
        let s = 0.05 / 3.0f64.sqrt();
        let p = pos_of(&out, 9);
        assert!((p[0] - s).abs() < EPS && (p[1] - s).abs() < EPS && (p[2] - s).abs() < EPS);
        let (prisms, tets) = volume_split(&out);
        assert!((prisms + tets - 1.0).abs() < EPS);
        assert_all_positive(&out);
    }

    // ── Spec validation / no-op / invalid input ──────────────────────────────

    #[test]
    fn zero_layers_passthrough_and_invalid_specs_rejected() {
        let cube = box_surface(1.0, 1.0, 1.0);
        let out = insert_boundary_layers_3d(&cube, &BoundaryLayer3D::new(0, 0.1, 1.0))
            .expect("zero layers is ok");
        assert_eq!(out.node_count(), cube.node_count());
        assert_eq!(out.element_count(), cube.element_count());

        let bad = BoundaryLayer3D::new(2, 0.0, 1.0);
        assert!(matches!(
            insert_boundary_layers_3d(&cube, &bad),
            Err(BoundaryLayer3DError::InvalidSpec(_))
        ));
        let bad = BoundaryLayer3D::new(2, 0.1, 0.0);
        assert!(matches!(
            insert_boundary_layers_3d(&cube, &bad),
            Err(BoundaryLayer3DError::InvalidSpec(_))
        ));
        let mut bad = BoundaryLayer3D::new(2, 0.1, 1.0);
        bad.clamp_safety = 0.9; // must be <= 0.5
        assert!(matches!(
            insert_boundary_layers_3d(&cube, &bad),
            Err(BoundaryLayer3DError::InvalidSpec(_))
        ));
    }

    #[test]
    fn open_surface_rejected() {
        let mut m = Mesh::new();
        m.add_node(Node::new(1, 0.0, 0.0, 0.0));
        m.add_node(Node::new(2, 1.0, 0.0, 0.0));
        m.add_node(Node::new(3, 0.0, 1.0, 0.0));
        m.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));
        assert!(matches!(
            insert_boundary_layers_3d(&m, &BoundaryLayer3D::new(1, 0.1, 1.0)),
            Err(BoundaryLayer3DError::InvalidInput(_))
        ));
    }

    // ── Round trip through rmsh-io (write MSH v2 + parse) ────────────────────

    #[test]
    fn layered_mesh_round_trips_through_msh_v2() {
        let cube = box_surface(1.0, 1.0, 1.0);
        let spec = BoundaryLayer3D::new(3, 0.05, 1.2);
        let out = insert_boundary_layers_3d(&cube, &spec).expect("layers insert");

        let mut bytes = Vec::new();
        write_msh_v2(&mut bytes, &out).expect("msh v2 write");
        let back = parse_msh(bytes.as_slice()).expect("msh v2 parse");

        assert_eq!(back.node_count(), out.node_count());
        assert_eq!(back.element_count(), out.element_count());
        assert_eq!(count_type(&back, ElementType::Prism6), 36);
        assert_eq!(count_type(&back, ElementType::Triangle3), 12);
        assert_eq!(count_type(&back, ElementType::Tetrahedron4), 12);
    }
}
