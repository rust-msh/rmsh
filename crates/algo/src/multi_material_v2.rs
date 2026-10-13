//! The conforming multi-material carve (v2, ROADMAP P4 upgrade): a tagged
//! tetrahedral mesh of ALL material regions sharing EXACT interface
//! triangulations.
//!
//! # What this module adds over v1
//!
//! v1 ([`crate::multi_material`]) classifies whole tets by centroid parity
//! and deletes the cavity tets — the cavity wall is a jagged
//! centroid-carve surface, the whole-tet straddle corrupts the volume
//! bookkeeping, and the void is not meshed. v2 upgrades the pipeline for
//! CONVEX cavity surfaces:
//!
//! 1. tetrahedralize the outer surface (the existing centroid-star kernel,
//!    unchanged);
//! 2. split every kernel tet recursively by the cavity face planes — the
//!    classic tetrahedron-by-plane cut with CLOSED-FORM decompositions: a
//!    1+3 corner split yields one corner tet plus a triangular prism
//!    (3 tets); a 2+2 split yields two prisms (3 tets each). Every piece
//!    stays a LIST OF TETS — no polyhedron face bookkeeping, no winding
//!    drift, no cracked loops: a cut corner (the crossing of one tet edge
//!    with the plane) is computed ONCE and shared by both sides verbatim;
//! 3. classify each final tet by its corners: all four on the inner side
//!    of every plane of cavity k -> cavity k's material; otherwise the
//!    outer material. Tets are TAGGED, nothing is deleted.
//!
//! # Why the interfaces come out conforming
//!
//! The tets on the two sides of a cavity plane share the cut corners
//! bit-identically (same edge, same interpolation, one computation), so
//! their faces on the plane weld into shared triangles. Over all cut tets,
//! the cuts tile each input cavity triangle (every point of the triangle
//! is on the boundary between cavity-side and outer-side material): the
//! input cavity triangulation is preserved as REGIONS (each triangle is
//! exactly tiled by interface faces — a mutual refinement, not the
//! verbatim fan). A per-point symbolic tie-break (a coordinate-hash
//! micro-offset, the same for the same geometric point in every tet)
//! resolves the knife-edge degeneracies of cascading cuts consistently.
//!
//! The whole-tet straddle vanishes: `sum(outer tets) + sum(cavity k tets)
//! == enclosed(outer)` holds to floating-point exactness (pinned 1e-9) —
//! every cut is volume-conserving by construction (the prism split is
//! verified against the tet volume in the tests, including a 200-plane
//! pseudo-random sweep).
//!
//! # The degeneracy band and the measured-coverage gate
//!
//! Cavity planes grazing the kernel structure (passing within the sliver
//! band of kernel vertices/edges) generate knife-edge cut cascades whose
//! near-degenerate shards cannot all be emitted (the near-degenerate
//! volume floor removes them, and with them their corner patches). The
//! routing measures the achieved interface coverage per cavity (shared
//! cross-material face area against the cavity surface area) and falls
//! back to the v1 path when ANY cavity falls below 98% — an honest
//! straddled wall instead of a leaky conforming one. Clean layouts (the
//! nested-cube corpus) pass the gate and get the exact carve.
//!
//! # Scope (honest)
//!
//! Convex cavities only. A non-convex cavity cannot be carved by its own
//! face half-spaces (the complement pieces stop being convex); the entry
//! routes such layouts to the v1 path with its documented straddle. The
//! GENERAL cavity case is the 3-D constrained-Delaunay face-recovery
//! problem; the reference route for that card is TetGen's constrained
//! edge/face recovery by flip sequences (H. Si, "TetGen, a Delaunay-Based
//! Quality Tetrahedral Mesh Generator", ACM TOMS 41(2), 2015 — algorithm
//! concepts only, clean-room).
//!
//! # Node welding
//!
//! Corner coordinates reach the output through many cut paths, so output
//! nodes are welded through a quantized key grid at `1e-12 x scene
//! diagonal`; the first occurrence's exact coordinates win. Distinct nodes
//! closer than the grid merge (documented tolerance); the volume identity
//! is pinned far above the grid.
//!
//! # Determinism
//!
//! Tets are produced in kernel-tet order x plane order x decomposition
//! order; node ids in first-touch order. No HashMap iteration reaches the
//! output.

use std::collections::HashMap;

use rmsh_model::{Element, ElementType, Mesh, Node};
use thiserror::Error;

/// Quantized node-welding grid: 1e-12 of the scene diagonal.
const WELD_GRID_REL: f64 = 1e-12;

/// The deduplicated face planes of one convex cavity.
#[derive(Debug)]
pub(crate) struct CavityPlan {
    pub planes: Vec<([f64; 3], f64)>,
}

pub(crate) type CarveResult<T> = Result<T, CarveError>;

#[derive(Error, Debug)]
pub(crate) enum CarveError {
    #[error("cavity carve failed: {0}")]
    Generation(String),
}

/// Normalized plane of an outward triangle `(a, b, c)`.
fn tri_plane(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Option<([f64; 3], f64)> {
    let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let n = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    let ln = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if ln <= 1e-30 {
        return None;
    }
    let n = [n[0] / ln, n[1] / ln, n[2] / ln];
    Some((n, -(n[0] * a[0] + n[1] * a[1] + n[2] * a[2])))
}

/// Quantized plane key for deduplication (1e-9 in normal components,
/// 1e-9 x scene scale in the offset).
fn plane_key(n: [f64; 3], d: f64, scale: f64) -> (i64, i64, i64, i64) {
    const Q: f64 = 1e9;
    (
        (n[0] * Q).round() as i64,
        (n[1] * Q).round() as i64,
        (n[2] * Q).round() as i64,
        (d / scale.max(1e-30) * Q).round() as i64,
    )
}

/// Convexity of a fan triangulation: every triangle plane must support the
/// whole surface (all vertices on the non-positive side, tolerance `tol`).
pub(crate) fn check_convex(
    sigmas: &[[[f64; 3]; 3]],
    points: &[[f64; 3]],
    tol: f64,
) -> Result<(), String> {
    for s in sigmas {
        let Some((n, d)) = tri_plane(s[0], s[1], s[2]) else {
            continue;
        };
        for p in points {
            if n[0] * p[0] + n[1] * p[1] + n[2] * p[2] + d > tol {
                return Err(format!(
                    "cavity surface is not convex (vertex ({}, {}, {}) lies \
                     outside a face plane; the conforming carve requires \
                     convex cavities)",
                    p[0], p[1], p[2]
                ));
            }
        }
    }
    Ok(())
}

/// Build the carve plan of one validated convex cavity from its fan
/// triangles (the deduplicated face planes only).
pub(crate) fn build_plan(sigmas: &[[[f64; 3]; 3]], scale: f64) -> CavityPlan {
    let mut planes: Vec<([f64; 3], f64)> = Vec::new();
    let mut seen: HashMap<(i64, i64, i64, i64), ()> = HashMap::new();
    for s in sigmas {
        let Some((n, d)) = tri_plane(s[0], s[1], s[2]) else {
            continue;
        };
        let key = plane_key(n, d, scale);
        if seen.insert(key, ()).is_some() {
            continue;
        }
        planes.push((n, d));
    }
    CavityPlan { planes }
}

/// A welded node pool: quantized keys, first-occurrence coordinates.
pub(crate) struct NodePool {
    scale: f64,
    map: HashMap<(i64, i64, i64), usize>,
    pub coords: Vec<[f64; 3]>,
}

impl NodePool {
    pub fn new(scale: f64) -> Self {
        NodePool {
            scale: scale.max(1e-30),
            map: HashMap::new(),
            coords: Vec::new(),
        }
    }

    fn key(&self, p: [f64; 3]) -> (i64, i64, i64) {
        let g = WELD_GRID_REL * self.scale;
        (
            (p[0] / g).round() as i64,
            (p[1] / g).round() as i64,
            (p[2] / g).round() as i64,
        )
    }

    pub fn node(&mut self, p: [f64; 3]) -> usize {
        let k = self.key(p);
        if let Some(&i) = self.map.get(&k) {
            return i;
        }
        self.coords.push(p);
        self.map.insert(k, self.coords.len() - 1);
        self.coords.len() - 1
    }
}

/// One edge crossing of the plane (shared by both cut pieces verbatim).
struct Cut {
    p: [f64; 3],
}

fn plane_value(plane: &([f64; 3], f64), p: [f64; 3]) -> f64 {
    plane.0[0] * p[0] + plane.0[1] * p[1] + plane.0[2] * p[2] + plane.1
}

/// Split one tet by one plane. Returns `(neg_tets, pos_tets)` — the tets on
/// the negative side (the cavity interior side) and on the positive side.
/// Conserves volume exactly: the 1+3 case emits the corner tet plus the
/// canonical 3-tet prism split; the 2+2 case emits two canonical 3-tet
/// prism splits; on-plane corners snap to the positive side (the negative
/// side stays strictly interior), and crossings onto a snapped corner
/// degenerate to the corner itself (volumeless shards are dropped).
fn split_tet(
    tet: &[[f64; 3]; 4],
    plane: &([f64; 3], f64),
    eps: f64,
) -> (Vec<[[f64; 3]; 4]>, Vec<[[f64; 3]; 4]>) {
    // Symbolic tie-break: a corner lying (numerically) ON the plane gets a
    // deterministic micro-offset derived from its QUANTIZED COORDINATES —
    // the same geometric point receives the same offset in every tet that
    // touches it, so cascading cuts stay consistent (no cracks) and no
    // sliver-generating snaps are needed beyond true degeneracies.
    let scale = eps * 1e9;
    let bias = |p: &[f64; 3]| -> f64 {
        let g = WELD_GRID_REL * scale.max(1e-30);
        let kx = (p[0] / g).round() as i64;
        let ky = (p[1] / g).round() as i64;
        let kz = (p[2] / g).round() as i64;
        let h = (kx.wrapping_mul(73856093) ^ ky.wrapping_mul(19349663) ^ kz.wrapping_mul(83492791))
            .rem_euclid(7);
        // the magnitude must stay small enough that edge-crossings computed
        // from different edges of the same geometric corner land within the
        // node-welding grid (they weld into ONE shared node), yet large
        // enough to dominate the fp-level cascade noise (~1e-16 relative)
        (h as f64 - 3.0) * 0.34 * 1e-14 * scale
    };
    let vals = [
        plane_value(plane, tet[0]) + bias(&tet[0]),
        plane_value(plane, tet[1]) + bias(&tet[1]),
        plane_value(plane, tet[2]) + bias(&tet[2]),
        plane_value(plane, tet[3]) + bias(&tet[3]),
    ];
    // true-degeneracy snapping threshold (far below the classification eps)
    let eps_snap = 1e-15 * scale.max(1e-30);
    let mut neg_idx: Vec<usize> = Vec::new();
    let mut pos_idx: Vec<usize> = Vec::new();
    for (i, &v) in vals.iter().enumerate() {
        if v < -eps {
            neg_idx.push(i);
        } else {
            pos_idx.push(i);
        }
    }
    if neg_idx.is_empty() {
        return (Vec::new(), vec![*tet]);
    }
    if pos_idx.is_empty() {
        return (vec![*tet], Vec::new());
    }

    // crossing point of the plane with the edge (a, b): a
    // true-degenerate endpoint snaps to itself
    let crossing = |a: usize, b: usize| -> [f64; 3] {
        let (va, vb) = (vals[a], vals[b]);
        let (pa, pb) = (tet[a], tet[b]);
        if va.abs() <= eps_snap {
            return pa;
        }
        if vb.abs() <= eps_snap {
            return pb;
        }
        let t = va / (va - vb);
        [
            pa[0] + t * (pb[0] - pa[0]),
            pa[1] + t * (pb[1] - pa[1]),
            pa[2] + t * (pb[2] - pa[2]),
        ]
    };
    // canonical 3-tet decomposition of the triangular prism
    // (P1,P2,P3 ; Q1,Q2,Q3) with the corner correspondence Pk <-> Qk
    let prism = |p1: [f64; 3],
                 p2: [f64; 3],
                 p3: [f64; 3],
                 q1: [f64; 3],
                 q2: [f64; 3],
                 q3: [f64; 3]|
     -> Vec<[[f64; 3]; 4]> {
        vec![[p1, p2, p3, q3], [p1, p2, q3, q2], [p1, q3, q2, q1]]
    };

    let mut neg_out: Vec<[[f64; 3]; 4]> = Vec::new();
    let mut pos_out: Vec<[[f64; 3]; 4]> = Vec::new();
    if neg_idx.len() == 1 {
        // one negative corner A; positive B, C, D; crossings on AB, AC, AD
        let a = neg_idx[0];
        let (b, c, d) = (pos_idx[0], pos_idx[1], pos_idx[2]);
        let u = crossing(a, b);
        let v = crossing(a, c);
        let w = crossing(a, d);
        neg_out.push([tet[a], u, v, w]);
        // the truncated rest: prism (B, C, D ; u, v, w)
        pos_out.extend(prism(tet[b], tet[c], tet[d], u, v, w));
    } else if pos_idx.len() == 1 {
        // mirror: one positive corner D; negative A, B, C; crossings on
        // AD, BD, CD
        let d = pos_idx[0];
        let (a, b, c) = (neg_idx[0], neg_idx[1], neg_idx[2]);
        let u = crossing(a, d);
        let v = crossing(b, d);
        let w = crossing(c, d);
        pos_out.push([tet[d], u, v, w]);
        // the truncated rest: prism (A, B, C ; u, v, w)
        neg_out.extend(prism(tet[a], tet[b], tet[c], u, v, w));
    } else {
        // two negative corners A, B; positive C, D; crossings on AC, AD,
        // BC, BD
        let (a, b) = (neg_idx[0], neg_idx[1]);
        let (c, d) = (pos_idx[0], pos_idx[1]);
        let pac = crossing(a, c);
        let pad = crossing(a, d);
        let pbc = crossing(b, c);
        let pbd = crossing(b, d);
        // negative wedge: prism (A, pac, pad ; B, pbc, pbd)
        neg_out.extend(prism(tet[a], pac, pad, tet[b], pbc, pbd));
        // positive wedge: prism (C, pac, pbc ; D, pad, pbd)
        pos_out.extend(prism(tet[c], pac, pbc, tet[d], pad, pbd));
    }
    // NOTE: no volume filter here — a flattened (near-zero-volume) piece
    // can still carry REAL interface area along cavity edges/corners;
    // dropping it would punch pinholes into the shared wall. Degenerate
    // pieces are filtered at emission by the welded-node distinctness
    // check, and their volume contribution is below the pinned identity
    // tolerance either way.
    (neg_out, pos_out)
}

/// Run the full conforming carve.
///
/// `outer_tets` is the kernel tetrahedralization of the outer surface
/// (Tetrahedron4 only, validated by the caller); `plans` carry the
/// validated convex cavities in input order; `scale` is the scene
/// diagonal. Cavity k's tets carry the tag returned by
/// `material_of_cavity`; tets outside every cavity carry
/// [`crate::multi_material::MATERIAL_TAG_OUTER`].
pub(crate) fn carve_conforming(
    outer_tets: &Mesh,
    plans: &[CavityPlan],
    scale: f64,
    material_of_cavity: impl Fn(usize) -> i32,
) -> CarveResult<CarvedMesh> {
    let eps = 1e-9 * scale;
    let mut pool = NodePool::new(scale);
    let mut tets: Vec<([usize; 4], i32)> = Vec::new();

    for elem in &outer_tets.elements {
        if elem.dimension() != 3 || elem.node_ids.len() != 4 {
            continue;
        }
        let get = |id: u64| -> [f64; 3] {
            let n = outer_tets.nodes.get(&id).expect("kernel tet node");
            [n.position.x, n.position.y, n.position.z]
        };
        let seed = [
            get(elem.node_ids[0]),
            get(elem.node_ids[1]),
            get(elem.node_ids[2]),
            get(elem.node_ids[3]),
        ];
        // normalize the winding once (positive signed volume)
        let (a, b, c, d) = (seed[0], seed[1], seed[2], seed[3]);
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let ad = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
        let cx = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        let signed = cx[0] * ad[0] + cx[1] * ad[1] + cx[2] * ad[2];
        let mut current: Vec<[[f64; 3]; 4]> =
            vec![if signed < 0.0 { [a, c, b, d] } else { [a, b, c, d] }];
        for plan in plans {
            for plane in &plan.planes {
                let mut next: Vec<[[f64; 3]; 4]> = Vec::new();
                for t in &current {
                    let (neg, pos) = split_tet(t, plane, eps);
                    next.extend(neg);
                    next.extend(pos);
                }
                current = next;
            }
        }
        // classify by corners: ALL corners on the inner side of every plane
        // of cavity k -> cavity k's material
        for t in current {
            let mut owner = None;
            for (k, plan) in plans.iter().enumerate() {
                if plan
                    .planes
                    .iter()
                    .all(|p| t.iter().all(|&v| plane_value(p, v) <= eps))
                {
                    owner = Some(k);
                    break;
                }
            }
            let tag = match owner {
                None => crate::multi_material::MATERIAL_TAG_OUTER,
                Some(k) => material_of_cavity(k),
            };
            let mut ids = [pool.node(t[0]), pool.node(t[1]), pool.node(t[2]), pool.node(t[3])];
            // welding can collapse near-degenerate shards into repeated
            // node ids: such a tet carries no volume and no real face
            // area — drop it
            if ids[0] == ids[1]
                || ids[0] == ids[2]
                || ids[0] == ids[3]
                || ids[1] == ids[2]
                || ids[1] == ids[3]
                || ids[2] == ids[3]
            {
                continue;
            }
            // normalize the winding to positive signed volume
            let (p0, p1, p2, p3) = (
                pool.coords[ids[0]],
                pool.coords[ids[1]],
                pool.coords[ids[2]],
                pool.coords[ids[3]],
            );
            let ab = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            let ac = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
            let ad = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
            let cx = [
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0],
            ];
            if cx[0] * ad[0] + cx[1] * ad[1] + cx[2] * ad[2] < 0.0 {
                ids.swap(2, 3);
            }
            // near-degenerate shards (the cascade's knife-edge cuts) carry
            // no volume and stack duplicate faces around cavity edges and
            // corners: drop them. The floor sits three orders below the
            // pinned volume-identity tolerance.
            let vol_floor = 1e-13 * scale * scale * scale;
            let vol = (cx[0] * ad[0] + cx[1] * ad[1] + cx[2] * ad[2]).abs();
            if vol <= vol_floor {
                continue;
            }
            tets.push((ids, tag));
        }
    }

    // assemble
    let mut mesh = Mesh::new();
    for (i, c) in pool.coords.iter().enumerate() {
        mesh.add_node(Node::new(i as u64 + 1, c[0], c[1], c[2]));
    }
    for (i, (nodes, tag)) in tets.iter().enumerate() {
        let mut e = Element::new(
            i as u64 + 1,
            ElementType::Tetrahedron4,
            nodes.iter().map(|&n| n as u64 + 1).collect(),
        );
        e.physical_tag = Some(*tag);
        mesh.add_element(e);
    }
    Ok(CarvedMesh { mesh })
}

/// The conforming carve output.
pub(crate) struct CarvedMesh {
    pub mesh: Mesh,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vol6(t: [[f64; 3]; 4]) -> f64 {
        let (a, b, c, d) = (t[0], t[1], t[2], t[3]);
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let ad = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
        let cx = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        (cx[0] * ad[0] + cx[1] * ad[1] + cx[2] * ad[2]).abs()
    }

    #[test]
    fn split_tet_conserves_volume_one_vs_three() {
        let tet = [[0.5, 0.5, 0.5], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let plane = ([1.0, 0.0, 0.0], -0.25);
        let (neg, pos) = split_tet(&tet, &plane, 1e-9);
        let seed = vol6(tet);
        let sum: f64 = neg.iter().chain(pos.iter()).map(|t| vol6(*t)).sum();
        assert!(
            (sum - seed).abs() < 1e-12,
            "split must conserve: seed {seed} vs sum {sum}"
        );
        let pos_sum: f64 = pos.iter().map(|t| vol6(*t)).sum();
        assert!(pos_sum > seed * 0.5, "the truncated rest dominates");
    }

    #[test]
    fn split_tet_conserves_volume_two_plus_two() {
        let tet = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let plane = ([0.0, 1.0, 1.0], -0.6);
        let (neg, pos) = split_tet(&tet, &plane, 1e-9);
        let sum: f64 = neg.iter().chain(pos.iter()).map(|t| vol6(*t)).sum();
        assert!((sum - vol6(tet)).abs() < 1e-12, "2+2 split must conserve");
        assert_eq!(neg.len(), 3, "the wedge is a 3-tet prism split");
        assert_eq!(pos.len(), 3);
    }

    #[test]
    fn split_tet_random_planes_conserve() {
        let tet = [[0.1, 0.1, 0.1], [0.9, 0.2, 0.3], [0.3, 0.8, 0.2], [0.25, 0.3, 0.85]];
        let seed = vol6(tet);
        let mut state = 0x2545F4914F6CDD1Du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..200 {
            let n = [
                (next() % 1000) as f64 / 1000.0 - 0.5,
                (next() % 1000) as f64 / 1000.0 - 0.5,
                (next() % 1000) as f64 / 1000.0 - 0.5,
            ];
            let d = (next() % 1000) as f64 / 1000.0 - 0.5;
            let ln = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            if ln < 1e-3 {
                continue;
            }
            let plane = ([n[0] / ln, n[1] / ln, n[2] / ln], d);
            let (neg, pos) = split_tet(&tet, &plane, 1e-9);
            let sum: f64 = neg.iter().chain(pos.iter()).map(|t| vol6(*t)).sum();
            assert!(
                (sum - seed).abs() < 1e-12,
                "random-plane split must conserve: seed {seed} vs sum {sum}"
            );
        }
    }
}

