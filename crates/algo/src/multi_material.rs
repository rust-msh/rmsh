//! Multi-material / partitioned volume meshing (ROADMAP P4): material
//! classification and cavity subtraction over closed surface meshes.
//!
//! # v1 scope limitation (read this first)
//!
//! TRUE interface conformity — adjacent regions' tet meshes sharing the
//! interface triangulation, the way Gmsh's `classifyTetrahedraInRegions`
//! achieves it by flood-filling tet adjacency between recovered constraint
//! faces — requires constrained insertion of the inner surfaces into the
//! outer tetrahedralization, i.e. a 3-D constrained Delaunay triangulation
//! (3-D CDT). The 3-D CDT is a future card. This module is the v1
//! deliverable on top of the existing centroid-star tetrahedralizer:
//!
//! 1. tetrahedralize the outer closed surface (unconstrained, existing path);
//! 2. classify every tet by its CENTROID against each inner closed surface
//!    by ray parity (policy pinned below);
//! 3. remove the tets inside an inner surface — the cavity void. The cavity
//!    wall is conceptually the inner surface's own triangulation; v1 does
//!    NOT emit it as elements and does NOT mesh the void region (a future
//!    card meshes the void region; the 3-D CDT card upgrades classification
//!    to constrained insertion for true conformity);
//! 4. surviving tets carry `physical_tag = Some(MATERIAL_TAG_OUTER)` (the
//!    outer material). Per-material meshes of the inner regions are produced
//!    by calling [`tetrahedralize_region`] on each inner surface separately;
//!    those meshes share no nodes or faces with the outer mesh (v1 is
//!    non-conforming at interfaces by construction).
//!
//! # Gmsh provenance (algorithm concepts only)
//!
//! The multi-region concepts follow Gmsh's volume pipeline (C. Geuzaine,
//! J.-F. Remacle, Gmsh, GPL): `classifyTetrahedraInRegions`
//! (meshGRegionDelaunay.cpp) classifies a whole-group recovered tet mesh
//! into per-region volumes by adjacency flood fill bounded by constraint
//! faces, deleting the tets that fall into "the void" (outside every named
//! region); meshGRegion.cpp meshes a GROUP of regions through one boundary
//! recovery over the union boundary. This module re-implements the
//! ALGORITHM CONCEPTS from that description — no source text was copied
//! (GPL -> MIT clean-room port). The v1 delta versus Gmsh is explicit:
//! Gmsh's flood fill is bounded by constraint faces (which is what makes
//! interfaces conforming); v1 has no constrained faces yet, so per-tet
//! centroid ray parity stands in for the flood fill and only Gmsh's
//! void-deletion branch is ported (cavity subtraction).
//!
//! # Classification policy (pinned)
//!
//! A query point (the tet centroid) is classified against each inner
//! surface's fan triangulation by CROSSING PARITY:
//!
//! - Rays are cast along +X, then +Y, then +Z, in that fixed order.
//! - Each ray-triangle test is Moeller-Trumbore. A triangle hit counts one
//!   crossing iff it lies strictly beyond the origin (t above a
//!   scene-relative tolerance) AND strictly inside the triangle
//!   (barycentric coordinates above a 1e-9 margin from every edge).
//! - A hit ON a triangle edge or vertex (within the margin) marks the
//!   parity for that direction UNRELIABLE (the classic ray-through-edge
//!   degeneracy); the cast is retried along the next axis.
//! - If all three axis directions are unreliable (adversarial geometry),
//!   the final arbiter is the winding number via signed solid angle
//!   (van Oosterom-Strackee): inside iff winding > 0.5.
//! - A ray parallel to a triangle's plane contributes no crossing (the
//!   plane never crosses the ray transversally). Points lying exactly ON a
//!   surface are tie-broken arbitrarily but deterministically by the
//!   cascade; on-surface INPUT nodes are the validation layer's concern
//!   (see the enclosure policy), not the classifier's.
//! - Parity is orientation-INDEPENDENT (crossing counts do not care about
//!   winding direction); only the winding fallback and the divergence
//!   volumes consume the consistent-orientation requirement.
//!
//! # Enclosure / disjointness policy (pinned, all refusals name the violation)
//!
//! v1 validates at NODE granularity and refuses conservatively:
//!
//! 1. every input (outer and each inner) must be a PURE surface mesh (2-D
//!    elements only) that is closed, manifold and CONSISTENTLY oriented
//!    (every undirected edge traversed exactly once in each direction);
//! 2. each inner's axis-aligned bounding box must be strictly inside the
//!    outer's (by a scene-relative tolerance). An inner surface TOUCHING
//!    the outer boundary necessarily touches at a node, hence touches the
//!    bounding boxes, hence is refused here naming the bounding-box
//!    violation;
//! 3. pairwise inner bounding boxes must not strictly overlap (cheap
//!    disjointness check; boxes merely sharing a face plane are allowed —
//!    v1 does not prove true surface disjointness, documented limitation);
//! 4. every inner node classifies INSIDE the outer surface by ray parity.
//!    For convex outers this layer is unreachable (a node outside a convex
//!    outer forces its box outside too); it is defense in depth for
//!    non-convex outers. v1 cannot detect an inner FACE piercing a
//!    non-convex outer with all its nodes inside — full surface-intersection
//!    validation is future work.
//!
//! # Volume bookkeeping
//!
//! [`region_volumes`] answers from the closed SURFACES by the divergence
//! theorem (signed tetrahedral sums over each fan triangulation):
//! `vec[0]` is the outer material volume (outer enclosed volume minus all
//! cavity volumes) and `vec[1..]` are the cavity volumes. On surfaces the
//! identity `material + sum(cavities) == outer enclosed volume` holds to
//! floating-point exactness (pinned 1e-9). The TET-level identity is
//! coarser by construction: cavity subtraction removes WHOLE tets, so
//! `sum(surviving tets) + sum(cavity surfaces)` deviates from the outer
//! volume by the whole-tet straddle error; that deviation is measured and
//! pinned per case in the tests, alongside the exact invariant
//! `surviving + removed == outer tet sum`.
//!
//! # Determinism
//!
//! The pipeline is deterministic within a process: fixed ray order, fixed
//! tie-break cascade, element generation in the deterministic boundary
//! order of [`crate::tetrahedralize3d::tetrahedralize_closed_surface`].
//! Cross-process BIT identity is not pinned (the node-hash iteration inside
//! the centroid computation uses std's per-process random seed); tags and
//! volumes are pinned at 1e-9 there.

use std::collections::HashMap;

use rmsh_model::{Element, Mesh, Point3};
use thiserror::Error;

use crate::traits::MeshAlgoError;
use crate::tetrahedralize3d::{Mesh3DError, tetrahedralize_closed_surface};

/// Physical tag carried by every surviving tet of
/// [`tetrahedralize_multi_material`]: the outer material.
pub const MATERIAL_TAG_OUTER: i32 = 0;

/// Errors of the multi-material pipeline. Every validation variant names
/// the violating input and the reason.
#[derive(Error, Debug)]
pub enum MultiMaterialError {
    /// The outer surface input is invalid.
    #[error("outer surface is invalid: {0}")]
    OuterSurface(String),

    /// An inner surface input is invalid.
    #[error("inner surface {index} is invalid: {reason}")]
    InnerSurface { index: usize, reason: String },

    /// Two inner surfaces violate the disjointness precondition.
    #[error("inner surface {index} violates disjointness: {reason}")]
    Disjointness { index: usize, reason: String },

    /// An inner surface violates the enclosure precondition.
    #[error("inner surface {index} is not enclosed: {reason}")]
    Enclosure { index: usize, reason: String },

    /// The volume meshing step failed.
    #[error("3D meshing failed: {0}")]
    Generation(String),
}

impl From<MultiMaterialError> for MeshAlgoError {
    fn from(value: MultiMaterialError) -> Self {
        match value {
            MultiMaterialError::Generation(msg) => MeshAlgoError::Generation(msg),
            other => MeshAlgoError::InvalidInput(other.to_string()),
        }
    }
}

impl From<Mesh3DError> for MultiMaterialError {
    fn from(value: Mesh3DError) -> Self {
        MultiMaterialError::Generation(value.to_string())
    }
}

/// Tetrahedralize the outer material of a multi-material layout.
///
/// `outer` is the enclosing closed surface; `inners` are zero or more
/// DISJOINT closed surfaces fully inside it (the cavity boundaries). The
/// outer volume is tetrahedralized on the existing centroid-star path and
/// every tet whose centroid classifies inside an inner surface is REMOVED
/// (the cavity void). Surviving tets carry
/// [`MATERIAL_TAG_OUTER`]; the cavity wall is not emitted as elements and
/// the void region is not meshed (v1 — see the module docs for the
/// conformity limitation and upgrade path).
///
/// An empty `inners` slice degenerates to
/// [`tetrahedralize_region`]'s result with the material tag applied.
/// A layout whose cavities consume every tet yields a legal empty mesh
/// (documented v1 policy: the caller validates material presence).
pub fn tetrahedralize_multi_material(
    outer: &Mesh,
    inners: &[Mesh],
) -> Result<Mesh, MultiMaterialError> {
    let outer_box = outer.bounding_box().ok_or_else(|| {
        MultiMaterialError::OuterSurface("outer surface has no nodes".to_string())
    })?;
    let tol = 1e-9 * (outer_box.1 - outer_box.0).norm();

    // Validation layer 1: the outer must be a pure, closed, manifold,
    // consistently oriented surface (the classification consumes its
    // triangulation; the volume identity consumes its orientation).
    let outer_polys = collect_surface_polygons(outer, None)?;
    check_oriented_closed(&outer_polys).map_err(MultiMaterialError::OuterSurface)?;
    let outer_tris = resolve_points(outer, &fan_triangulate(&outer_polys));

    // Validation layers 2-4 per inner: oriented closure, strict bounding-box
    // containment, pairwise disjointness, per-node parity against the outer
    // surface (module docs, "Enclosure / disjointness policy").
    let mut inner_boxes: Vec<(Point3, Point3)> = Vec::with_capacity(inners.len());
    let mut inner_tris: Vec<Vec<[Point3; 3]>> = Vec::with_capacity(inners.len());
    for (k, inner) in inners.iter().enumerate() {
        let polys = collect_surface_polygons(inner, Some(k))?;
        check_oriented_closed(&polys)
            .map_err(|reason| MultiMaterialError::InnerSurface { index: k, reason })?;

        let bbox = inner.bounding_box().ok_or_else(|| MultiMaterialError::InnerSurface {
            index: k,
            reason: "inner surface has no nodes".to_string(),
        })?;
        for (j, &other) in inner_boxes.iter().enumerate() {
            if aabb_overlap_strict(bbox, other, tol) {
                return Err(MultiMaterialError::Disjointness {
                    index: k,
                    reason: format!("bounding box overlaps inner surface {j}"),
                });
            }
        }
        if !aabb_strictly_inside(bbox, outer_box, tol) {
            return Err(MultiMaterialError::Enclosure {
                index: k,
                reason: "bounding box is not strictly inside the outer bounding box \
                         (the inner surface touches or crosses the outer boundary)"
                    .to_string(),
            });
        }

        // Layer 4: every inner node must classify inside the outer surface
        // (first-appearance order keeps the refusal deterministic).
        let tris = resolve_points(inner, &fan_triangulate(&polys));
        let mut node_order: Vec<u64> = Vec::new();
        for poly in &polys {
            for nid in poly {
                if !node_order.contains(nid) {
                    node_order.push(*nid);
                }
            }
        }
        for nid in node_order {
            let pos = inner.nodes.get(&nid).expect("validated node").position;
            if classify_point(pos, &outer_tris, tol) == Classification::Outside {
                return Err(MultiMaterialError::Enclosure {
                    index: k,
                    reason: format!("node {nid} classifies outside the outer surface"),
                });
            }
        }

        inner_boxes.push(bbox);
        inner_tris.push(tris);
    }

    // Step 1 of the pipeline: tetrahedralize the outer surface on the
    // existing centroid-star path.
    let tet_mesh = tetrahedralize_closed_surface(outer)?;

    // Steps 2-3: classify every tet by centroid parity against each inner
    // surface and drop the tets inside a cavity. Step 4: tag survivors as
    // the outer material and renumber deterministically in generation order.
    let mut out = Mesh::new();
    let mut used_nodes: Vec<u64> = Vec::new();
    let mut next_id = 1u64;
    for elem in &tet_mesh.elements {
        if elem.dimension() != 3 {
            continue;
        }
        let coords = |id: u64| tet_mesh.nodes.get(&id).expect("tet node present").position;
        let c = Point3::from(
            (coords(elem.node_ids[0]).coords
                + coords(elem.node_ids[1]).coords
                + coords(elem.node_ids[2]).coords
                + coords(elem.node_ids[3]).coords)
                / 4.0,
        );
        let mut removed = false;
        for tris in &inner_tris {
            if classify_point(c, tris, tol) == Classification::Inside {
                removed = true;
                break;
            }
        }
        if removed {
            continue;
        }
        let mut kept = Element::new(next_id, elem.etype.clone(), elem.node_ids.clone());
        kept.physical_tag = Some(MATERIAL_TAG_OUTER);
        for nid in &kept.node_ids {
            if !used_nodes.contains(nid) {
                used_nodes.push(*nid);
            }
        }
        out.add_element(kept);
        next_id += 1;
    }
    used_nodes.sort_unstable();
    for nid in used_nodes {
        let node = tet_mesh.nodes.get(&nid).expect("tet node present");
        out.add_node(node.clone());
    }
    Ok(out)
}

/// The per-material meshing path: mesh ONE closed region surface on its own.
///
/// This is the composition contract of the multi-material v1 entry: to also
/// mesh the inner regions, call this on each inner surface separately. The
/// resulting meshes share no nodes or faces with each other or with the
/// [`tetrahedralize_multi_material`] output — v1 interfaces are
/// non-conforming by construction (module docs, first paragraph).
pub fn tetrahedralize_region(surface: &Mesh) -> Result<Mesh, Mesh3DError> {
    tetrahedralize_closed_surface(surface)
}

/// Divergence-theorem volume bookkeeping for a multi-material layout.
///
/// Returns `[material_volume, cavity_volume_1, ..., cavity_volume_k]` where
/// `material_volume = enclosed(outer) - sum(enclosed(inner_k))` and each
/// cavity volume is the signed divergence sum over the inner surface's own
/// fan triangulation. The surfaces must be closed, manifold and
/// consistently oriented (validated); semantic enclosure is NOT validated
/// here — [`tetrahedralize_multi_material`] is the validating entry.
///
/// Adjudication note: the roadmap sketch named `region_volumes(outer_tets,
/// inners)`; the entry deliberately takes the closed SURFACES instead — tet
/// sums would make the bookkeeping identity inexact by the whole-tet
/// straddle, while surface divergence keeps
/// `sum(answer) == enclosed(outer)` floating-point exact (pinned 1e-9).
/// The tet-level identity is pinned separately in the tests at the measured
/// straddle tolerance.
pub fn region_volumes(outer: &Mesh, inners: &[Mesh]) -> Result<Vec<f64>, MultiMaterialError> {
    let outer_polys = collect_surface_polygons(outer, None)?;
    check_oriented_closed(&outer_polys).map_err(MultiMaterialError::OuterSurface)?;
    let outer_vol =
        enclosed_signed_volume(&resolve_points(outer, &fan_triangulate(&outer_polys)));

    let mut vols = Vec::with_capacity(inners.len() + 1);
    let mut cavity_sum = 0.0;
    for (k, inner) in inners.iter().enumerate() {
        let polys = collect_surface_polygons(inner, Some(k))?;
        check_oriented_closed(&polys)
            .map_err(|reason| MultiMaterialError::InnerSurface { index: k, reason })?;
        let vol = enclosed_signed_volume(&resolve_points(inner, &fan_triangulate(&polys)));
        vols.push(vol);
        cavity_sum += vol;
    }
    vols.insert(0, outer_vol - cavity_sum);
    Ok(vols)
}

// ─── Surface extraction and validation ────────────────────────────────────────

fn surface_error(index: Option<usize>, reason: String) -> MultiMaterialError {
    match index {
        None => MultiMaterialError::OuterSurface(reason),
        Some(index) => MultiMaterialError::InnerSurface { index, reason },
    }
}

/// Extract the 2-D boundary polygons (node ids, fan order) of a surface
/// mesh. Refuses volume elements (v1 requires pure surface input for every
/// multi-material participant) and missing nodes.
fn collect_surface_polygons(mesh: &Mesh, index: Option<usize>) -> Result<Vec<Vec<u64>>, MultiMaterialError> {
    for elem in &mesh.elements {
        if elem.dimension() == 3 {
            return Err(surface_error(
                index,
                format!(
                    "must be a pure surface mesh (found a volume element {:?})",
                    elem.etype
                ),
            ));
        }
    }

    let mut polys = Vec::new();
    for elem in &mesh.elements {
        if elem.dimension() != 2 {
            continue;
        }
        if elem.node_ids.len() < 3 {
            return Err(surface_error(
                index,
                format!("2D element {} has fewer than 3 nodes", elem.id),
            ));
        }
        for nid in &elem.node_ids {
            if !mesh.nodes.contains_key(nid) {
                return Err(surface_error(
                    index,
                    format!("node {} missing in source mesh", nid),
                ));
            }
        }
        polys.push(elem.node_ids.clone());
    }
    Ok(polys)
}

/// Check that a fan-ordered polygon set forms a closed, manifold,
/// consistently oriented surface: every undirected edge is traversed
/// exactly once in each direction.
fn check_oriented_closed(polys: &[Vec<u64>]) -> Result<(), String> {
    if polys.is_empty() {
        return Err("surface has no 2D elements".to_string());
    }

    // (min_id, max_id) -> (forward_count, backward_count) by traversal.
    let mut edges: HashMap<(u64, u64), (usize, usize)> = HashMap::new();
    for poly in polys {
        for i in 1..(poly.len() - 1) {
            for (a, b) in [(poly[0], poly[i]), (poly[i], poly[i + 1]), (poly[i + 1], poly[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                let entry = edges.entry(key).or_insert((0, 0));
                if a < b {
                    entry.0 += 1;
                } else {
                    entry.1 += 1;
                }
            }
        }
    }

    for ((a, b), &(forward, backward)) in &edges {
        if forward + backward != 2 {
            return Err(format!(
                "edge ({}, {}) is used {} times (surface is open or non-manifold)",
                a,
                b,
                forward + backward
            ));
        }
        if forward != 1 {
            return Err(format!(
                "edge ({}, {}) is traversed twice in the same direction (surface orientation is inconsistent)",
                a,
                b
            ));
        }
    }
    Ok(())
}

/// Fan-triangulate polygons (same fan convention as the centroid-star
/// tetrahedralizer: apex at `poly[0]`). Preserves orientation.
fn fan_triangulate(polys: &[Vec<u64>]) -> Vec<[u64; 3]> {
    let mut tris = Vec::new();
    for poly in polys {
        for i in 1..(poly.len() - 1) {
            tris.push([poly[0], poly[i], poly[i + 1]]);
        }
    }
    tris
}

/// Resolve triangle node ids to positions. The ids were validated by
/// [`collect_surface_polygons`], so a miss is a programming error.
fn resolve_points(mesh: &Mesh, tris: &[[u64; 3]]) -> Vec<[Point3; 3]> {
    tris.iter()
        .map(|t| {
            [
                mesh.nodes.get(&t[0]).expect("validated node").position,
                mesh.nodes.get(&t[1]).expect("validated node").position,
                mesh.nodes.get(&t[2]).expect("validated node").position,
            ]
        })
        .collect()
}

/// Signed enclosed volume of a fan triangulation by the divergence theorem:
/// sum of `dot(v0, cross(v1, v2)) / 6`. Positive for consistently
/// outward-oriented surfaces, negative for inward; meaningless garbage for
/// inconsistently oriented input (the closure check refuses those).
fn enclosed_signed_volume(tris: &[[Point3; 3]]) -> f64 {
    let mut sum6 = 0.0;
    for tri in tris {
        let v1 = tri[1] - tri[0];
        let v2 = tri[2] - tri[0];
        sum6 += tri[0].coords.dot(&v1.cross(&v2));
    }
    sum6 / 6.0
}

// ─── Point classification: ray parity with a documented degeneracy cascade ────

/// Outcome of one ray cast against a triangle soup.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RayCast {
    /// No grazing hit: the crossing count is reliable parity evidence.
    Parity(usize),
    /// The ray hit a triangle edge or vertex — parity is unreliable.
    Grazing,
}

/// The classification of a point against a closed surface triangulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Classification {
    /// Odd crossing parity (or winding > 0.5 in the fallback).
    Inside,
    /// Even crossing parity (or winding <= 0.5 in the fallback).
    Outside,
}

/// Barycentric margin separating an interior hit from an edge/vertex hit
/// (scale-free — barycentric coordinates).
const GRAZE_EPS: f64 = 1e-9;

/// The fixed tie-break ray order of the classification policy.
fn ray_directions() -> [Point3; 3] {
    [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
    ]
}

/// Moeller-Trumbore ray cast against a triangle soup. Hits at or behind the
/// origin (within `t_tol`) never count — the surface face the origin may lie
/// on is excluded by construction. A transversal hit within `GRAZE_EPS` of a
/// triangle edge marks the whole cast as grazing (short-circuit).
fn cast_ray(origin: Point3, dir: Point3, tris: &[[Point3; 3]], t_tol: f64) -> RayCast {
    let mut crossings = 0usize;
    for tri in tris {
        let e1 = tri[1] - tri[0];
        let e2 = tri[2] - tri[0];
        let p = dir.coords.cross(&e2);
        let det = e1.dot(&p);
        // Parallel (including degenerate zero-area triangles): no
        // transversal crossing, no graze — the plane never crosses the ray.
        if det.abs() <= 1e-14 * e1.norm() * e2.norm() {
            continue;
        }
        let inv_det = 1.0 / det;
        let s = origin - tri[0];
        let u = s.dot(&p) * inv_det;
        if u < -GRAZE_EPS || u > 1.0 + GRAZE_EPS {
            continue;
        }
        let q = s.cross(&e1);
        let v = dir.coords.dot(&q) * inv_det;
        if v < -GRAZE_EPS || u + v > 1.0 + GRAZE_EPS {
            continue;
        }
        let t = e2.dot(&q) * inv_det;
        if t <= t_tol {
            continue;
        }
        // A hit this close to an edge or vertex is a graze: parity along
        // this direction is unreliable, the caller retries the next axis.
        if u <= GRAZE_EPS || v <= GRAZE_EPS || u + v >= 1.0 - GRAZE_EPS {
            return RayCast::Grazing;
        }
        crossings += 1;
    }
    RayCast::Parity(crossings)
}

/// Signed solid angle winding number (van Oosterom-Strackee): the sum of
/// per-triangle signed solid angles divided by 4*pi. Inside a closed,
/// consistently oriented surface the value is 1; outside, 0. The final
/// arbiter of the classification cascade.
fn winding_number(origin: Point3, tris: &[[Point3; 3]]) -> f64 {
    let mut solid_angle_sum = 0.0;
    for tri in tris {
        let a = tri[0] - origin;
        let b = tri[1] - origin;
        let c = tri[2] - origin;
        let (la, lb, lc) = (a.norm(), b.norm(), c.norm());
        if la == 0.0 || lb == 0.0 || lc == 0.0 {
            // The origin coincides with a vertex: undefined, contributes 0
            // (on-surface points are the validation layer's concern).
            continue;
        }
        let triple = a.dot(&b.cross(&c));
        let denom = la * lb * lc + a.dot(&b) * lc + a.dot(&c) * lb + b.dot(&c) * la;
        if denom == 0.0 {
            continue;
        }
        let omega = 2.0 * triple.abs().atan2(denom);
        solid_angle_sum += triple.signum() * omega;
    }
    solid_angle_sum / (4.0 * std::f64::consts::PI)
}

/// Classify a point against a closed surface triangulation (module docs,
/// "Classification policy"): ray parity along +X/+Y/+Z with the grazing
/// cascade, winding-number fallback as the final arbiter.
fn classify_point(origin: Point3, tris: &[[Point3; 3]], t_tol: f64) -> Classification {
    for dir in ray_directions() {
        match cast_ray(origin, dir, tris, t_tol) {
            RayCast::Parity(n) => {
                return if n % 2 == 1 {
                    Classification::Inside
                } else {
                    Classification::Outside
                };
            }
            RayCast::Grazing => continue,
        }
    }
    if winding_number(origin, tris) > 0.5 {
        Classification::Inside
    } else {
        Classification::Outside
    }
}

// ─── Bounding-box predicates (scene-relative tolerances) ──────────────────────

/// Strict containment of `inner` in `outer` with a tolerance margin.
fn aabb_strictly_inside(
    inner: (Point3, Point3),
    outer: (Point3, Point3),
    tol: f64,
) -> bool {
    inner.0.x > outer.0.x + tol
        && inner.0.y > outer.0.y + tol
        && inner.0.z > outer.0.z + tol
        && inner.1.x < outer.1.x - tol
        && inner.1.y < outer.1.y - tol
        && inner.1.z < outer.1.z - tol
}

/// Strict overlap of two boxes (touching faces do not overlap).
fn aabb_overlap_strict(a: (Point3, Point3), b: (Point3, Point3), tol: f64) -> bool {
    a.0.x < b.1.x - tol
        && b.0.x < a.1.x - tol
        && a.0.y < b.1.y - tol
        && b.0.y < a.1.y - tol
        && a.0.z < b.1.z - tol
        && b.0.z < a.1.z - tol
}

#[cfg(test)]
mod tests {
    use super::*;
    use rmsh_model::{Element, ElementType, Node};

    /// A closed unit-ish cube surface `[lo, hi]` with 8 nodes and 6
    /// consistently OUTWARD-oriented quads.
    fn cube_surface(lo: [f64; 3], hi: [f64; 3]) -> Mesh {
        let [lx, ly, lz] = lo;
        let [hx, hy, hz] = hi;
        let mut m = Mesh::new();
        m.add_node(Node::new(1, lx, ly, lz));
        m.add_node(Node::new(2, hx, ly, lz));
        m.add_node(Node::new(3, hx, hy, lz));
        m.add_node(Node::new(4, lx, hy, lz));
        m.add_node(Node::new(5, lx, ly, hz));
        m.add_node(Node::new(6, hx, ly, hz));
        m.add_node(Node::new(7, hx, hy, hz));
        m.add_node(Node::new(8, lx, hy, hz));

        m.add_element(Element::new(1, ElementType::Quad4, vec![1, 4, 3, 2])); // z=lo, outward -z
        m.add_element(Element::new(2, ElementType::Quad4, vec![5, 6, 7, 8])); // z=hi, outward +z
        m.add_element(Element::new(3, ElementType::Quad4, vec![1, 2, 6, 5])); // y=lo, outward -y
        m.add_element(Element::new(4, ElementType::Quad4, vec![2, 3, 7, 6])); // x=hi, outward +x
        m.add_element(Element::new(5, ElementType::Quad4, vec![3, 4, 8, 7])); // y=hi, outward +y
        m.add_element(Element::new(6, ElementType::Quad4, vec![4, 1, 5, 8])); // x=lo, outward -x
        m
    }

    /// A hand-built outward-oriented tetrahedron surface:
    /// (0,0,0), (1,0,0), (0,1,0), (0,0,1).
    fn tet_tris() -> Vec<[Point3; 3]> {
        let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
        vec![
            // z=0 face, outward -z
            [p(0.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(1.0, 0.0, 0.0)],
            // y=0 face, outward -y
            [p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 0.0, 1.0)],
            // x=0 face, outward -x
            [p(0.0, 0.0, 0.0), p(0.0, 0.0, 1.0), p(0.0, 1.0, 0.0)],
            // slant x+y+z=1, outward +(+1,+1,+1)
            [p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(0.0, 0.0, 1.0)],
        ]
    }

    /// The unit octahedron |x|+|y|+|z| <= 1 as 8 outward triangles. A ray
    /// along any axis from the center exits through a VERTEX, so all three
    /// axis directions graze and only the winding fallback can classify.
    fn octahedron_tris() -> Vec<[Point3; 3]> {
        let mut tris = Vec::new();
        for sx in [-1.0, 1.0] {
            for sy in [-1.0, 1.0] {
                for sz in [-1.0, 1.0] {
                    let a = Point3::new(sx, 0.0, 0.0);
                    let b = Point3::new(0.0, sy, 0.0);
                    let c = Point3::new(0.0, 0.0, sz);
                    // The cross of (b-a, c-a) is (sy*sz, sx*sz, sx*sy);
                    // swap b/c for negative octants to keep it outward.
                    if sx * sy * sz > 0.0 {
                        tris.push([a, b, c]);
                    } else {
                        tris.push([a, c, b]);
                    }
                }
            }
        }
        tris
    }

    /// Sum of tet volumes of a volume mesh (absolute).
    fn tet_volume_sum(mesh: &Mesh) -> f64 {
        mesh.elements
            .iter()
            .filter(|e| e.dimension() == 3)
            .map(|e| {
                let get = |id: u64| mesh.nodes[&id].position;
                let a = get(e.node_ids[0]);
                let b = get(e.node_ids[1]);
                let c = get(e.node_ids[2]);
                let d = get(e.node_ids[3]);
                (b - a).cross(&(c - a)).dot(&(d - a)).abs() / 6.0
            })
            .sum()
    }

    // ── Ray-parity unit tests (hand-built tet with known answers) ─────────────

    #[test]
    fn ray_parity_tet_inside_outside_and_on_boundary_tie_break() {
        let tris = tet_tris();

        // Centroid region: strictly inside on every axis ray.
        assert_eq!(
            classify_point(Point3::new(0.25, 0.25, 0.25), &tris, 1e-12),
            Classification::Inside
        );
        // Far field: clearly outside.
        assert_eq!(
            classify_point(Point3::new(2.0, 2.0, 2.0), &tris, 1e-12),
            Classification::Outside
        );
        // A point ON the slant edge (1,0,0)-(0,0,1): the +X cast hits the
        // slant plane at t = 0 (excluded by the t tolerance), the y=0 and
        // z=0 planes are parallel to the ray, and the x=0 plane is behind.
        // Parity(0) -> Outside: the documented deterministic tie-break for
        // on-surface points (arbitrary but fixed).
        assert_eq!(
            classify_point(Point3::new(0.5, 0.0, 0.5), &tris, 1e-12),
            Classification::Outside
        );
    }

    #[test]
    fn ray_parity_cascade_octahedron_center_via_winding_fallback() {
        let tris = octahedron_tris();
        // Every axis ray from the origin exits through a vertex of the
        // octahedron: all three directions report Grazing and the winding
        // fallback must classify the center as inside.
        assert_eq!(
            classify_point(Point3::new(0.0, 0.0, 0.0), &tris, 1e-12),
            Classification::Inside
        );
        // Sanity: a point outside beyond a face classifies by plain parity.
        assert_eq!(
            classify_point(Point3::new(0.9, 0.9, 0.9), &tris, 1e-12),
            Classification::Outside
        );
    }

    #[test]
    fn winding_number_tet_center_and_far_field() {
        let tris = tet_tris();
        let w_in = winding_number(Point3::new(0.25, 0.25, 0.25), &tris);
        assert!(
            (w_in - 1.0).abs() < 1e-9,
            "winding inside a closed surface should be 1, got {w_in}"
        );
        let w_out = winding_number(Point3::new(3.0, 3.0, 3.0), &tris);
        assert!(
            w_out.abs() < 1e-9,
            "winding outside a closed surface should be 0, got {w_out}"
        );
    }

    #[test]
    fn ray_cast_grazing_is_detected_on_vertex_hits() {
        let tris = octahedron_tris();
        // The +X ray from the origin exits through the vertex (1,0,0):
        // every face hit lands on an edge/vertex, so the cast must report
        // Grazing rather than a (bogus) parity count.
        assert_eq!(
            cast_ray(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0), &tris, 1e-12),
            RayCast::Grazing
        );
    }

    // ── Surface validation unit tests ─────────────────────────────────────────

    #[test]
    fn oriented_closure_accepts_good_cube_and_detects_flip_and_hole() {
        let cube = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let polys = collect_surface_polygons(&cube, None).expect("cube is a pure surface");
        check_oriented_closed(&polys).expect("the outward cube is closed and oriented");

        // Flipping one face traversal must be caught as an orientation fault.
        let mut flipped = polys.clone();
        flipped[2].reverse();
        let err = check_oriented_closed(&flipped).expect_err("flipped face must be refused");
        assert!(
            err.contains("direction"),
            "flip should name the direction violation, got: {err}"
        );

        // Dropping one face must be caught as an open surface.
        let mut holed = polys.clone();
        holed.remove(1);
        let err = check_oriented_closed(&holed).expect_err("open surface must be refused");
        assert!(
            err.contains("open"),
            "hole should name the open/non-manifold violation, got: {err}"
        );

        // An empty polygon set has no 2D elements at all.
        let err = check_oriented_closed(&[]).expect_err("empty surface must be refused");
        assert!(err.contains("no 2D elements"), "got: {err}");
    }

    #[test]
    fn collect_surface_polygons_refuses_volume_elements_and_missing_nodes() {
        let mut mesh = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        mesh.add_element(Element::new(7, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));
        let err = collect_surface_polygons(&mesh, Some(0))
            .expect_err("volume elements must be refused for a surface input");
        match err {
            MultiMaterialError::InnerSurface { index, reason } => {
                assert_eq!(index, 0);
                assert!(reason.contains("pure surface"), "got: {reason}");
            }
            other => panic!("unexpected error: {other:?}"),
        }

        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 99]));
        let err = collect_surface_polygons(&mesh, None)
            .expect_err("missing nodes must be refused");
        match err {
            MultiMaterialError::OuterSurface(reason) => assert!(reason.contains("Node") || reason.contains("node 99") || reason.contains("99")),
            other => panic!("unexpected error: {other:?}"),
        }
    }

    // ── Volume and box predicate unit tests ───────────────────────────────────

    #[test]
    fn enclosed_signed_volume_cube_inner_and_reversed() {
        let cube = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let tris = {
            let polys = collect_surface_polygons(&cube, None).expect("pure surface");
            let pts = resolve_points(&cube, &fan_triangulate(&polys));
            pts
        };
        let vol = enclosed_signed_volume(&tris);
        assert!((vol - 1.0).abs() < 1e-9, "unit cube volume should be 1, got {vol}");

        let inner = cube_surface([0.3, 0.3, 0.3], [0.7, 0.7, 0.7]);
        let polys = collect_surface_polygons(&inner, None).expect("pure surface");
        let vol = enclosed_signed_volume(&resolve_points(&inner, &fan_triangulate(&polys)));
        assert!(
            (vol - 0.064).abs() < 1e-9,
            "0.4-cube volume should be 0.064, got {vol}"
        );

        // Reversed orientation flips the sign (documented semantics).
        let reversed: Vec<[Point3; 3]> = tris
            .iter()
            .map(|t| [t[0], t[2], t[1]])
            .collect();
        let vol = enclosed_signed_volume(&reversed);
        assert!((vol + 1.0).abs() < 1e-9, "reversed cube should be -1, got {vol}");
    }

    #[test]
    fn aabb_predicates_strict_containment_and_strict_overlap() {
        let tol = 1e-9;
        let outer = (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
        let inside = (Point3::new(0.3, 0.3, 0.3), Point3::new(0.7, 0.7, 0.7));
        let flush = (Point3::new(0.0, 0.0, 0.0), Point3::new(0.4, 0.4, 0.4));
        let poking = (Point3::new(0.5, 0.5, 0.5), Point3::new(1.5, 1.5, 1.5));

        assert!(aabb_strictly_inside(inside, outer, tol));
        // A flush box touches the outer boundary: not strictly inside.
        assert!(!aabb_strictly_inside(flush, outer, tol));
        assert!(!aabb_strictly_inside(poking, outer, tol));

        let left = (Point3::new(0.05, 0.05, 0.05), Point3::new(0.45, 0.45, 0.45));
        let right_touching = (Point3::new(0.45, 0.05, 0.05), Point3::new(0.85, 0.45, 0.45));
        let right_overlapping = (Point3::new(0.3, 0.05, 0.05), Point3::new(0.7, 0.45, 0.45));
        // Face-touching boxes do not strictly overlap (documented allowance).
        assert!(!aabb_overlap_strict(left, right_touching, tol));
        assert!(aabb_overlap_strict(left, right_overlapping, tol));
    }

    // ── Multi-material pipeline (red until the entries are implemented) ──────

    /// Nested cubes: unit-cube outer, 0.75-cube inner at the center. Six of
    /// the twelve centroid-star tets have all centroid coordinates in
    /// (0.1, 0.85) and must be REMOVED; the remaining six sum to exactly
    /// 0.5. Whole-tet straddle: |6/12 - 0.421875| = 0.078125.
    #[test]
    fn nested_cubes_cavity_removed_and_volume_identity_pinned() {
        let outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let inner = cube_surface([0.1, 0.1, 0.1], [0.85, 0.85, 0.85]);

        let out = tetrahedralize_multi_material(&outer, &[inner.clone()])
            .expect("nested cubes are a valid layout");

        let tets = out.elements_by_dimension(3);
        assert_eq!(tets.len(), 6, "6 of 12 tets must be removed, got {}", tets.len());

        let vol = tet_volume_sum(&out);
        assert!(
            (vol - 0.5).abs() < 1e-9,
            "surviving tets must sum to 6/12 = 0.5, got {vol}"
        );
        assert!(
            (vol + 0.421875 - 1.0).abs() <= 0.08,
            "surviving + cavity theory must match outer within the measured \
             whole-tet straddle 0.078125, got {}",
            vol + 0.421875 - 1.0
        );

        // Every surviving tet is tagged the outer material.
        for e in &out.elements {
            assert_eq!(e.physical_tag, Some(MATERIAL_TAG_OUTER));
        }

        // Consistency: every surviving centroid classifies outside the cavity.
        let polys = collect_surface_polygons(&inner, Some(0)).expect("pure surface");
        let tris = resolve_points(&inner, &fan_triangulate(&polys));
        for e in &tets {
            let get = |id: u64| out.nodes[&id].position;
            let c = Point3::from(
                (get(e.node_ids[0]).coords
                    + get(e.node_ids[1]).coords
                    + get(e.node_ids[2]).coords
                    + get(e.node_ids[3]).coords)
                    / 4.0,
            );
            assert_eq!(
                classify_point(c, &tris, 1e-9),
                Classification::Outside,
                "a surviving tet centroid classifies inside the cavity"
            );
        }
    }

    /// Two disjoint inner cavities: both are removed (3 tets total: T11 by
    /// cavity 1, T3+T4 by cavity 2 — exact centroid table), the cavity
    /// volumes are additive, and the identity is pinned at the measured
    /// straddle |3/12 - 0.088| = 0.162.
    #[test]
    fn two_disjoint_cavities_both_removed_volumes_additive() {
        let outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let inner1 = cube_surface([0.05, 0.05, 0.05], [0.45, 0.45, 0.45]);
        let inner2 = cube_surface([0.3, 0.3, 0.8], [0.7, 0.7, 0.95]);

        let out = tetrahedralize_multi_material(&outer, &[inner1, inner2])
            .expect("two disjoint cavities are a valid layout");

        let tets = out.elements_by_dimension(3);
        assert_eq!(
            tets.len(),
            9,
            "3 of 12 tets must be removed (1 by cavity 1, 2 by cavity 2)"
        );

        let vol = tet_volume_sum(&out);
        assert!((vol - 9.0 / 12.0).abs() < 1e-9, "9 surviving tets sum to 9/12");

        let cavity_theory: f64 = 0.4 * 0.4 * 0.4 + 0.4 * 0.4 * 0.15;
        assert!(
            (cavity_theory - 0.088).abs() < 1e-12,
            "cavity volumes are additive: 0.064 + 0.024 = 0.088"
        );
        assert!(
            (vol + cavity_theory - 1.0).abs() <= 0.17,
            "identity pinned at the measured straddle 0.162, got {}",
            vol + cavity_theory - 1.0
        );
    }

    /// An inner surface touching the outer boundary is refused, naming the
    /// violation (the flush face forces bounding-box contact, which the
    /// strict-containment layer rejects).
    #[test]
    fn inner_touching_outer_boundary_is_refused() {
        let outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        // Flush with the outer faces x=0, y=0, z=0.
        let inner = cube_surface([0.0, 0.0, 0.0], [0.4, 0.4, 0.4]);

        let err = tetrahedralize_multi_material(&outer, &[inner])
            .expect_err("a touching inner surface must be refused");
        match &err {
            MultiMaterialError::Enclosure { index, reason } => {
                assert_eq!(*index, 0);
                assert!(
                    reason.contains("bounding box"),
                    "refusal must name the violation, got: {reason}"
                );
            }
            other => panic!("unexpected error: {other:?}"),
        }
        let msg = err.to_string();
        assert!(msg.contains("inner surface 0"), "got: {msg}");
    }

    /// An inner surface whose bounding box is not strictly inside the
    /// outer's is refused (partial overlap with the outer boundary).
    #[test]
    fn inner_aabb_not_strictly_inside_is_refused() {
        let outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let inner = cube_surface([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);

        let err = tetrahedralize_multi_material(&outer, &[inner]).expect_err("must refuse");
        match &err {
            MultiMaterialError::Enclosure { index, reason } => {
                assert_eq!(*index, 0);
                assert!(reason.contains("bounding box"), "got: {reason}");
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    /// Pairwise inner bounding boxes that strictly overlap are refused,
    /// naming both indices; boxes merely sharing a face plane are allowed
    /// (documented v1 limitation: true surface disjointness is not proven).
    #[test]
    fn pairwise_inner_disjointness_is_enforced_and_touching_allowed() {
        let outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let a = cube_surface([0.05, 0.05, 0.05], [0.45, 0.45, 0.45]);
        let b = cube_surface([0.3, 0.05, 0.05], [0.7, 0.45, 0.45]);

        let err = tetrahedralize_multi_material(&outer, &[a, b])
            .expect_err("overlapping cavities must be refused");
        match &err {
            MultiMaterialError::Disjointness { index, reason } => {
                assert_eq!(*index, 1, "the second inner is the offender");
                assert!(reason.contains("overlaps"), "got: {reason}");
            }
            other => panic!("unexpected error: {other:?}"),
        }

        // Face-touching boxes are allowed through (limitation pin). The
        // exact centroid table gives 3 removals: T11 by box A and T2+T5 by
        // box B ([0.45,0.85]x[0.05,0.45]x[0.05,0.45] catches the centroids
        // (0.625,0.375,0.125) and (0.625,0.125,0.375)).
        let a = cube_surface([0.05, 0.05, 0.05], [0.45, 0.45, 0.45]);
        let b = cube_surface([0.45, 0.05, 0.05], [0.85, 0.45, 0.45]);
        let out = tetrahedralize_multi_material(&outer, &[a, b])
            .expect("face-touching cavity boxes are allowed in v1");
        assert_eq!(out.elements_by_dimension(3).len(), 9);
    }

    /// An empty inners slice must reproduce the single-surface call
    /// exactly (connectivity, nodes, volumes) — with the documented delta
    /// that the multi-material entry tags every tet as the outer material.
    #[test]
    fn empty_inners_slice_equals_single_surface_call() {
        let outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);

        let plain = tetrahedralize_region(&outer).expect("single-surface call must succeed");
        let multi = tetrahedralize_multi_material(&outer, &[])
            .expect("empty inners must succeed");

        let plain_tets = plain.elements_by_dimension(3);
        let multi_tets = multi.elements_by_dimension(3);
        assert_eq!(plain_tets.len(), multi_tets.len());
        for (p, m) in plain_tets.iter().zip(multi_tets.iter()) {
            assert_eq!(p.node_ids, m.node_ids, "connectivity must be identical");
        }
        assert_eq!(plain.node_count(), multi.node_count());
        for (id, node) in &plain.nodes {
            let other = multi.nodes.get(id).expect("same node set");
            assert_eq!(node.position, other.position);
        }
        assert!(
            (tet_volume_sum(&plain) - tet_volume_sum(&multi)).abs() < 1e-12,
            "volumes must be identical"
        );
        // The documented delta: the multi-material entry tags the material.
        assert!(plain_tets.iter().all(|e| e.physical_tag.is_none()));
        assert!(
            multi_tets
                .iter()
                .all(|e| e.physical_tag == Some(MATERIAL_TAG_OUTER))
        );
    }

    /// Same input, two calls: identical elements (ids, tags, connectivity)
    /// and identical node positions — single-process determinism.
    #[test]
    fn classification_is_deterministic_across_calls() {
        let outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let inner = cube_surface([0.1, 0.1, 0.1], [0.85, 0.85, 0.85]);

        let a = tetrahedralize_multi_material(&outer, &[inner.clone()])
            .expect("call A must succeed");
        let b = tetrahedralize_multi_material(&outer, &[inner])
            .expect("call B must succeed");

        assert_eq!(a.elements.len(), b.elements.len());
        for (ea, eb) in a.elements.iter().zip(b.elements.iter()) {
            assert_eq!(ea.id, eb.id);
            assert_eq!(ea.physical_tag, eb.physical_tag);
            assert_eq!(ea.node_ids, eb.node_ids);
        }
        let mut na: Vec<_> = a.nodes.iter().map(|(id, n)| (*id, n.position)).collect();
        let mut nb: Vec<_> = b.nodes.iter().map(|(id, n)| (*id, n.position)).collect();
        na.sort_unstable_by_key(|(id, _)| *id);
        nb.sort_unstable_by_key(|(id, _)| *id);
        assert_eq!(na, nb, "node positions must be bit-identical across calls");
    }

    /// A cavity consuming every tet is a legal empty mesh (documented v1
    /// policy: the caller validates material presence).
    #[test]
    fn fully_consumed_material_yields_empty_mesh() {
        let outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let inner = cube_surface([0.1, 0.1, 0.1], [0.9, 0.9, 0.9]);

        let out = tetrahedralize_multi_material(&outer, &[inner])
            .expect("a fully consumed material is legal");
        assert_eq!(out.elements_by_dimension(3).len(), 0);
        assert_eq!(out.node_count(), 0);
    }

    /// Volume-element inputs are refused for outer and inner alike (v1
    /// requires pure surface input; the boundary-extraction path belongs to
    /// the single-surface entry).
    #[test]
    fn volume_element_inputs_are_refused() {
        let mut outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        outer.add_element(Element::new(7, ElementType::Tetrahedron4, vec![1, 2, 3, 5]));
        let err = tetrahedralize_multi_material(&outer, &[])
            .expect_err("volume-element outer must be refused");
        assert!(
            matches!(&err, MultiMaterialError::OuterSurface(r) if r.contains("pure surface")),
            "got: {err:?}"
        );

        let outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let mut inner = cube_surface([0.3, 0.3, 0.3], [0.7, 0.7, 0.7]);
        inner.add_element(Element::new(7, ElementType::Tetrahedron4, vec![1, 2, 3, 5]));
        let err = tetrahedralize_multi_material(&outer, &[inner])
            .expect_err("volume-element inner must be refused");
        assert!(
            matches!(&err, MultiMaterialError::InnerSurface { index: 0, reason }
                if reason.contains("pure surface")),
            "got: {err:?}"
        );
    }

    /// Divergence-theorem bookkeeping: material + cavities == outer, pinned
    /// at 1e-9, with each surface's volume matching its analytic value.
    #[test]
    fn region_volumes_identity_pinned_1e9() {
        let outer = cube_surface([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let inner = cube_surface([0.3, 0.3, 0.3], [0.7, 0.7, 0.7]);

        let vols = region_volumes(&outer, &[inner]).expect("bookkeeping must succeed");
        assert_eq!(vols.len(), 2, "one material volume + one cavity volume");
        assert!((vols[0] - 0.936).abs() < 1e-9, "material = 1 - 0.064, got {}", vols[0]);
        assert!((vols[1] - 0.064).abs() < 1e-9, "cavity = 0.064, got {}", vols[1]);
        let total: f64 = vols.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-9,
            "material + cavities must equal the outer enclosed volume, got {total}"
        );

        // Empty inners: a single entry equal to the outer enclosed volume.
        let vols = region_volumes(&outer, &[]).expect("bookkeeping must succeed");
        assert_eq!(vols.len(), 1);
        assert!((vols[0] - 1.0).abs() < 1e-9);
    }
}
