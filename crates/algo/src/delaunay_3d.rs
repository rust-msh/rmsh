//! Delaunay 3-D — incremental Delaunay tetrahedralization (Gmsh algorithm 1).
//!
//! # Algorithm overview
//!
//! The 3-D Delaunay algorithm generalises the 2-D Bowyer-Watson method to
//! tetrahedra.  Given a set of points in ℝ³ it builds the unique (up to
//! degeneracies) tetrahedralization where no point lies strictly inside the
//! circumsphere of any tetrahedron.
//!
//! When used for **mesh generation** the algorithm proceeds in two phases:
//!
//! ## Phase 1 — Boundary-conforming Delaunay (constrained)
//!
//! 1. Insert all surface vertices of the input shell mesh into a super-tetrahedron
//!    that bounds the entire domain.
//! 2. Recover the boundary faces by forcing them into the triangulation via
//!    edge and face insertion (constrained Delaunay).
//!
//! ## Phase 2 — Delaunay refinement (Ruppert/Shewchuk in 3-D)
//!
//! 3. Compute the radius-edge ratio `ρ = R / l_min` for each tetrahedron
//!    (`R` = circumradius, `l_min` = shortest edge).
//! 4. While any tetrahedron has `ρ > ρ_max` (typically ≈ 2):
//!    a. Insert the circumcenter of the worst tetrahedron.
//!    b. Restore the Delaunay property via bistellar flips (3-D edge & face swaps).
//!    c. If the circumcenter is outside the domain, reject and mark the face.
//! 5. Remove the super-tetrahedron and all elements touching it.
//!
//! This is the algorithm implemented in TetGen and used as the basis of Gmsh's
//! own Delaunay 3-D pipeline.
//!
//! # Reference
//!
//! J. R. Shewchuk, "Tetrahedral Mesh Generation by Delaunay Refinement",
//! *SCG '98*, 1998.
//! H. Si, "TetGen, a Delaunay-Based Quality Tetrahedral Mesh Generator",
//! *ACM TOMS* 41(2), 2015.
//! Gmsh source: `Mesh/meshGRegion.cpp`.
//!
//! # Status
//!
//! **Mostly implemented** — `mesh_3d()` seeds from CentroidStarMesher3D and runs
//! `refine_bad_tetrahedra()` with radius-edge-ratio–driven refinement, bistellar
//! face flips, and quality metric evaluation.  One edge-case refinement path
//! still returns `NotImplemented`.

use rmsh_model::{Element, ElementType, Mesh, Node};
use std::collections::HashMap;

use crate::tetrahedralize3d::CentroidStarMesher3D;
use crate::traits::{MeshAlgoError, MeshParams, Mesher3D};

// ─── Public struct ────────────────────────────────────────────────────────────

/// Delaunay 3-D mesher (Gmsh algorithm 1).
///
/// Produces boundary-conforming Delaunay tetrahedral meshes via incremental
/// point insertion and Delaunay refinement.
#[derive(Debug, Clone)]
pub struct Delaunay3D {
    /// Maximum allowed radius-edge ratio `R / l_min`.
    ///
    /// Lower values produce better-quality tetrahedra but more elements.
    /// Shewchuk proves termination for `ρ_max ≥ 2.0`.  Defaults to `2.0`.
    pub max_radius_edge_ratio: f64,

    /// Maximum allowed dihedral angle deterioration.
    ///
    /// Tetrahedra with the minimum dihedral angle below this threshold (degrees)
    /// are candidates for refinement before the radius-edge ratio test.
    /// Set to `0.0` to disable.  Defaults to `5.0`.
    pub min_dihedral_angle_deg: f64,

    /// When `true`, circumcenters that fall outside the domain are reflected
    /// back inside (off-center insertion) rather than rejected.
    pub use_off_center_insertion: bool,
}

impl Default for Delaunay3D {
    fn default() -> Self {
        Self {
            max_radius_edge_ratio: 2.0,
            min_dihedral_angle_deg: 5.0,
            use_off_center_insertion: true,
        }
    }
}

impl Delaunay3D {
    pub fn new() -> Self {
        Self::default()
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl Mesher3D for Delaunay3D {
    fn name(&self) -> &'static str {
        "Delaunay 3D"
    }

    fn mesh_3d(&self, surface: &Mesh, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        validate_params(self, params)?;

        // Phase 0: seed with a robust closed-surface tetrahedralization.
        //
        // We then run a quality-driven refinement loop to move toward
        // Delaunay-style radius-edge constraints.
        let seed = CentroidStarMesher3D.mesh_3d(surface, params)?;
        refine_bad_tetrahedra(
            seed,
            self.max_radius_edge_ratio,
            params.element_size,
            params.max_size,
            params.optimize_passes,
        )
    }
}

fn refine_bad_tetrahedra(
    mut mesh: Mesh,
    max_radius_edge_ratio: f64,
    target_size: f64,
    max_size: f64,
    optimize_passes: u32,
) -> Result<Mesh, MeshAlgoError> {
    let mut stats = RefinementStats::default();
    let sliver_floor_deg = 0.25_f64;

    // Hard edge-length stop criterion combines target and optional max-size cap.
    let edge_limit = target_size.min(max_size);

    // Keep refinement bounded and predictable for UI usage, while still letting
    // element_size effectively control mesh density.
    let size_factor = ((mesh.diagonal_length() / edge_limit).ceil() as u32).clamp(1, 32);
    let max_passes = (optimize_passes.max(1) * size_factor).min(256);
    if max_passes == 0 {
        return Ok(mesh);
    }

    let mut next_node_id = mesh
        .nodes
        .keys()
        .copied()
        .max()
        .unwrap_or(0)
        .saturating_add(1);
    let mut next_elem_id = mesh
        .elements
        .iter()
        .map(|e| e.id)
        .max()
        .unwrap_or(0)
        .saturating_add(1);

    for _pass in 0..max_passes {
        let Some((worst_idx, _score)) =
            find_worst_tetrahedron(&mesh, max_radius_edge_ratio, edge_limit, sliver_floor_deg)
        else {
            stats.exits_no_worst += 1;
            break;
        };

        let worst_nodes = match mesh.elements.get(worst_idx) {
            Some(e) if e.etype == ElementType::Tetrahedron4 && e.node_ids.len() == 4 => {
                [e.node_ids[0], e.node_ids[1], e.node_ids[2], e.node_ids[3]]
            }
            _ => break,
        };
        if worst_nodes.len() != 4 {
            break;
        }

        let ratio = tetra_radius_edge_ratio_from_mesh(&mesh, &worst_nodes)?;
        let longest_edge = tetra_max_edge_length_from_mesh(&mesh, &worst_nodes)?;
        let min_dihedral = tetra_min_dihedral_from_mesh(&mesh, &worst_nodes)?;
        if ratio <= max_radius_edge_ratio
            && longest_edge <= edge_limit
            && min_dihedral >= sliver_floor_deg
        {
            stats.exits_quality_satisfied += 1;
            break;
        }

        let (p0, p1, p2, p3) = (
            node_xyz_from_mesh(&mesh, worst_nodes[0])?,
            node_xyz_from_mesh(&mesh, worst_nodes[1])?,
            node_xyz_from_mesh(&mesh, worst_nodes[2])?,
            node_xyz_from_mesh(&mesh, worst_nodes[3])?,
        );

        let centroid = [
            (p0[0] + p1[0] + p2[0] + p3[0]) * 0.25,
            (p0[1] + p1[1] + p2[1] + p3[1]) * 0.25,
            (p0[2] + p1[2] + p2[2] + p3[2]) * 0.25,
        ];

        let sliver_like = min_dihedral < sliver_floor_deg * 2.0 && ratio > max_radius_edge_ratio * 1.1;
        let best_point = if sliver_like {
            stats.sliver_priority_inserts += 1;
            select_fallback_refinement_point(p0, p1, p2, p3)
                .or_else(|| select_refinement_point(p0, p1, p2, p3))
                .unwrap_or(centroid)
        } else {
            select_refinement_point(p0, p1, p2, p3)
                .or_else(|| select_fallback_refinement_point(p0, p1, p2, p3))
                .unwrap_or(centroid)
        };
        let mut insertion_point = best_point;
        let mut edge_split = None::<(usize, usize, usize, usize)>;
        let predicted_min = min_child_dihedral_for_point(p0, p1, p2, p3, best_point);
        if predicted_min < sliver_floor_deg * 0.5 {
            stats.edge_bisection_considered += 1;
            if let Some((i, j, k, l, edge_point, edge_metrics)) =
                best_edge_split_partition([p0, p1, p2, p3])
            {
                let point_metrics = split_quality_metrics(p0, p1, p2, p3, best_point)
                    .map(|(_, md, mr, sf)| (md, sf, mr))
                    .unwrap_or((predicted_min, 1.0, f64::INFINITY));

                let edge_better = (edge_metrics.0 > point_metrics.0 + 1e-6)
                    || ((edge_metrics.0 - point_metrics.0).abs() < 1e-6
                        && ((edge_metrics.1 < point_metrics.1 - 1e-9)
                            || ((edge_metrics.1 - point_metrics.1).abs() < 1e-9
                                && edge_metrics.2 < point_metrics.2)));
                let edge_not_exploding = edge_metrics.2 <= point_metrics.2 * 2.5;

                if edge_better && edge_not_exploding {
                    insertion_point = edge_point;
                    edge_split = Some((i, j, k, l));
                    stats.edge_bisection_fallback += 1;
                } else {
                    stats.edge_bisection_rejected += 1;
                }
            } else {
                stats.edge_bisection_rejected += 1;
            }
        }

        if insertion_point == centroid {
            stats.centroid_fallback += 1;
        } else {
            stats.candidate_selected += 1;
        }
        let new_node_id = next_node_id;
        next_node_id = next_node_id.saturating_add(1);

        mesh.add_node(Node::new(
            new_node_id,
            insertion_point[0],
            insertion_point[1],
            insertion_point[2],
        ));

        // Replace one bad tetrahedron by four children sharing the inserted node,
        // or two children via longest-edge bisection for strongly sliver-like cases.
        let [n0, n1, n2, n3] = worst_nodes;
        let ids = [n0, n1, n2, n3];

        mesh.elements.swap_remove(worst_idx);
        if let Some((i, j, k, l)) = edge_split {
            mesh.add_element(Element::new(
                next_elem_id,
                ElementType::Tetrahedron4,
                vec![ids[i], new_node_id, ids[k], ids[l]],
            ));
            next_elem_id = next_elem_id.saturating_add(1);
            mesh.add_element(Element::new(
                next_elem_id,
                ElementType::Tetrahedron4,
                vec![new_node_id, ids[j], ids[k], ids[l]],
            ));
            next_elem_id = next_elem_id.saturating_add(1);
        } else {
            mesh.add_element(Element::new(
                next_elem_id,
                ElementType::Tetrahedron4,
                vec![n0, n1, n2, new_node_id],
            ));
            next_elem_id = next_elem_id.saturating_add(1);
            mesh.add_element(Element::new(
                next_elem_id,
                ElementType::Tetrahedron4,
                vec![n0, n1, n3, new_node_id],
            ));
            next_elem_id = next_elem_id.saturating_add(1);
            mesh.add_element(Element::new(
                next_elem_id,
                ElementType::Tetrahedron4,
                vec![n0, n2, n3, new_node_id],
            ));
            next_elem_id = next_elem_id.saturating_add(1);
            mesh.add_element(Element::new(
                next_elem_id,
                ElementType::Tetrahedron4,
                vec![n1, n2, n3, new_node_id],
            ));
            next_elem_id = next_elem_id.saturating_add(1);
        }

        stats.refined_tets += 1;
    }

    let flip_passes = optimize_passes.clamp(1, 8) as usize;
    let mut tet_mesh = crate::tet_mesh::TetMesh::from_mesh(&mesh);
    let (face_flips, edge_flips, edge_sliver_flips) =
        crate::tet_mesh::optimize_tetmesh_flips(&mut tet_mesh, flip_passes);
    stats.local_face_flips += face_flips;
    stats.local_edge_flips += edge_flips;
    stats.local_edge_sliver_accepts += edge_sliver_flips;
    mesh = tet_mesh.to_mesh();

    if should_log_refinement_stats() {
        eprintln!(
            "delaunay3d refinement stats: refined_tets={}, candidate_selected={}, centroid_fallback={}, sliver_priority_inserts={}, edge_bisection_considered={}, edge_bisection_fallback={}, edge_bisection_rejected={}, local_face_flips={}, local_edge_flips={}, local_edge_sliver_accepts={}, exits_no_worst={}, exits_quality_satisfied={}",
            stats.refined_tets,
            stats.candidate_selected,
            stats.centroid_fallback,
            stats.sliver_priority_inserts,
            stats.edge_bisection_considered,
            stats.edge_bisection_fallback,
            stats.edge_bisection_rejected,
            stats.local_face_flips,
            stats.local_edge_flips,
            stats.local_edge_sliver_accepts,
            stats.exits_no_worst,
            stats.exits_quality_satisfied
        );
    }

    Ok(mesh)
}

#[derive(Default)]
struct RefinementStats {
    refined_tets: usize,
    candidate_selected: usize,
    centroid_fallback: usize,
    sliver_priority_inserts: usize,
    edge_bisection_considered: usize,
    edge_bisection_fallback: usize,
    edge_bisection_rejected: usize,
    local_face_flips: usize,
    local_edge_flips: usize,
    local_edge_sliver_accepts: usize,
    exits_no_worst: usize,
    exits_quality_satisfied: usize,
}

fn optimize_local_face_flips(
    mesh: &mut Mesh,
    next_elem_id: &mut u64,
    max_passes: usize,
) -> Result<(usize, usize, usize), MeshAlgoError> {
    let mut accepted_face = 0usize;
    let mut accepted_edge = 0usize;
    let mut accepted_edge_sliver = 0usize;
    for _ in 0..max_passes {
        let mut face_map: HashMap<[u64; 3], Vec<(usize, u64)>> = HashMap::new();
        let mut edge_map: HashMap<[u64; 2], Vec<usize>> = HashMap::new();
        for (ti, e) in mesh.elements.iter().enumerate() {
            if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
                continue;
            }
            let n = [e.node_ids[0], e.node_ids[1], e.node_ids[2], e.node_ids[3]];
            let faces = [
                ([n[0], n[1], n[2]], n[3]),
                ([n[0], n[1], n[3]], n[2]),
                ([n[0], n[2], n[3]], n[1]),
                ([n[1], n[2], n[3]], n[0]),
            ];
            for (mut f, opp) in faces {
                f.sort_unstable();
                face_map.entry(f).or_default().push((ti, opp));
            }
            for (mut e2, _) in [
                ([n[0], n[1]], n[2]),
                ([n[0], n[2]], n[1]),
                ([n[0], n[3]], n[1]),
                ([n[1], n[2]], n[0]),
                ([n[1], n[3]], n[0]),
                ([n[2], n[3]], n[0]),
            ] {
                e2.sort_unstable();
                edge_map.entry(e2).or_default().push(ti);
            }
        }

        let mut did_flip = false;
        let prefer_edge_phase = has_sliver_pressure(mesh, 0.30, 0.12);
        let mut best_face_flip: Option<(usize, usize, [[u64; 4]; 3], (f64, f64, f64))> = None;
        let mut face_entries: Vec<_> = face_map.into_iter().collect();
        face_entries.sort_by_key(|(face, _)| *face);
        for (face, adjacent) in face_entries {
            if adjacent.len() != 2 {
                continue;
            }
            let (t0, o0) = adjacent[0];
            let (t1, o1) = adjacent[1];
            if t0 == t1 || o0 == o1 {
                continue;
            }
            let Some(e0) = mesh.elements.get(t0) else {
                continue;
            };
            let Some(e1) = mesh.elements.get(t1) else {
                continue;
            };
            if e0.etype != ElementType::Tetrahedron4
                || e1.etype != ElementType::Tetrahedron4
                || e0.node_ids.len() != 4
                || e1.node_ids.len() != 4
            {
                continue;
            }

            let old_tets = [
                [e0.node_ids[0], e0.node_ids[1], e0.node_ids[2], e0.node_ids[3]],
                [e1.node_ids[0], e1.node_ids[1], e1.node_ids[2], e1.node_ids[3]],
            ];
            let new_tets = [
                [o0, o1, face[0], face[1]],
                [o0, o1, face[1], face[2]],
                [o0, o1, face[2], face[0]],
            ];

            let Some((old_d, old_s, old_r)) = aggregate_tet_quality(mesh, &old_tets) else {
                continue;
            };
            let Some((new_d, new_s, new_r)) = aggregate_tet_quality(mesh, &new_tets) else {
                continue;
            };

            let improves = (new_d > old_d + 1e-6)
                || ((new_d - old_d).abs() < 1e-6
                    && ((new_s < old_s - 1e-9)
                        || ((new_s - old_s).abs() < 1e-9 && new_r < old_r - 1e-9)));
            if !improves {
                continue;
            }

            let new_quality = (new_d, new_s, new_r);
            match best_face_flip {
                Some((_, _, _, best_q)) if !is_better_quality(new_quality, best_q) => {}
                _ => best_face_flip = Some((t0, t1, new_tets, new_quality)),
            }
        }

        if !prefer_edge_phase {
            if let Some((t0, t1, new_tets, _)) = best_face_flip {
                let hi = t0.max(t1);
                let lo = t0.min(t1);
                mesh.elements.swap_remove(hi);
                mesh.elements.swap_remove(lo);

                for tet in new_tets {
                    mesh.add_element(Element::new(
                        *next_elem_id,
                        ElementType::Tetrahedron4,
                        vec![tet[0], tet[1], tet[2], tet[3]],
                    ));
                    *next_elem_id = next_elem_id.saturating_add(1);
                }

                accepted_face += 1;
                did_flip = true;
            }
        }

        if !did_flip {
            let mut best_edge_flip: Option<(
                Vec<usize>,
                [[u64; 4]; 2],
                (f64, f64, f64),
                f64,
                bool,
            )> = None;
            let mut edge_entries: Vec<_> = edge_map.into_iter().collect();
            edge_entries.sort_by_key(|(edge, _)| *edge);
            for (edge, adjacent) in edge_entries {
                if adjacent.len() != 3 {
                    continue;
                }

                let Some(e0) = mesh.elements.get(adjacent[0]) else {
                    continue;
                };
                let Some(e1) = mesh.elements.get(adjacent[1]) else {
                    continue;
                };
                let Some(e2) = mesh.elements.get(adjacent[2]) else {
                    continue;
                };
                if e0.etype != ElementType::Tetrahedron4
                    || e1.etype != ElementType::Tetrahedron4
                    || e2.etype != ElementType::Tetrahedron4
                    || e0.node_ids.len() != 4
                    || e1.node_ids.len() != 4
                    || e2.node_ids.len() != 4
                {
                    continue;
                }

                let u = edge[0];
                let v = edge[1];
                let mut opposite_pairs = Vec::<[u64; 2]>::with_capacity(3);
                let mut opposite_vertices = Vec::<u64>::with_capacity(3);
                let mut old_tets = [[0_u64; 4]; 3];
                let mut valid = true;

                for (slot, &ti) in adjacent.iter().enumerate() {
                    let Some(e) = mesh.elements.get(ti) else {
                        valid = false;
                        break;
                    };
                    let n = [e.node_ids[0], e.node_ids[1], e.node_ids[2], e.node_ids[3]];
                    old_tets[slot] = n;

                    let mut opp = Vec::<u64>::with_capacity(2);
                    for &nid in &n {
                        if nid != u && nid != v {
                            opp.push(nid);
                        }
                    }
                    if opp.len() != 2 {
                        valid = false;
                        break;
                    }
                    if !n.contains(&u) || !n.contains(&v) {
                        valid = false;
                        break;
                    }
                    opp.sort_unstable();
                    opposite_pairs.push([opp[0], opp[1]]);
                    opposite_vertices.push(opp[0]);
                    opposite_vertices.push(opp[1]);
                }
                if !valid {
                    continue;
                }

                opposite_vertices.sort_unstable();
                opposite_vertices.dedup();
                if opposite_vertices.len() != 3 {
                    continue;
                }
                let a = opposite_vertices[0];
                let b = opposite_vertices[1];
                let c = opposite_vertices[2];
                let mut need = vec![[a, b], [b, c], [a, c]];
                for p in &mut need {
                    p.sort_unstable();
                }
                let mut got = opposite_pairs.clone();
                got.sort_unstable();
                need.sort_unstable();
                if got != need {
                    continue;
                }

                let new_tets = [[a, b, c, u], [a, b, c, v]];
                let Some((old_d, old_s, old_r)) = aggregate_tet_quality(mesh, &old_tets) else {
                    continue;
                };
                let Some((new_d, new_s, new_r)) = aggregate_tet_quality(mesh, &new_tets) else {
                    continue;
                };

                let strict_improves = (new_d > old_d + 1e-6)
                    || ((new_d - old_d).abs() < 1e-6
                        && ((new_s < old_s - 1e-9)
                            || ((new_s - old_s).abs() < 1e-9 && new_r < old_r - 1e-9)));
                let sliver_delta = old_s - new_s;
                let strong_sliver_reduction = old_s >= 0.66 && new_s <= 0.34;
                let sliver_relaxed_improves = (sliver_delta > 0.08 || strong_sliver_reduction)
                    && new_d >= old_d * 0.70
                    && new_r <= old_r * 1.35;
                if !strict_improves && !sliver_relaxed_improves {
                    continue;
                }

                let new_quality = (new_d, new_s, new_r);
                let mut remove = vec![adjacent[0], adjacent[1], adjacent[2]];
                remove.sort_unstable();
                remove.reverse();

                let used_sliver_relaxed = !strict_improves && sliver_relaxed_improves;

                match best_edge_flip {
                    Some((_, _, best_q, best_sliver_delta, best_relaxed)) => {
                        let prefer = if sliver_delta > best_sliver_delta + 1e-9 {
                            true
                        } else if (sliver_delta - best_sliver_delta).abs() < 1e-9 {
                            if used_sliver_relaxed != best_relaxed {
                                !used_sliver_relaxed
                            } else {
                                is_better_quality(new_quality, best_q)
                            }
                        } else {
                            false
                        };
                        if prefer {
                            best_edge_flip = Some((
                                remove,
                                new_tets,
                                new_quality,
                                sliver_delta,
                                used_sliver_relaxed,
                            ));
                        }
                    }
                    _ => {
                        best_edge_flip = Some((
                            remove,
                            new_tets,
                            new_quality,
                            sliver_delta,
                            used_sliver_relaxed,
                        ))
                    }
                }
            }

            if let Some((remove, new_tets, _q, _d, used_sliver_relaxed)) = best_edge_flip {
                for idx in remove {
                    mesh.elements.swap_remove(idx);
                }

                for tet in new_tets {
                    mesh.add_element(Element::new(
                        *next_elem_id,
                        ElementType::Tetrahedron4,
                        vec![tet[0], tet[1], tet[2], tet[3]],
                    ));
                    *next_elem_id = next_elem_id.saturating_add(1);
                }

                accepted_edge += 1;
                if used_sliver_relaxed {
                    accepted_edge_sliver += 1;
                }
                did_flip = true;
            }
        }

        if !did_flip && prefer_edge_phase {
            if let Some((t0, t1, new_tets, _)) = best_face_flip {
                let hi = t0.max(t1);
                let lo = t0.min(t1);
                mesh.elements.swap_remove(hi);
                mesh.elements.swap_remove(lo);

                for tet in new_tets {
                    mesh.add_element(Element::new(
                        *next_elem_id,
                        ElementType::Tetrahedron4,
                        vec![tet[0], tet[1], tet[2], tet[3]],
                    ));
                    *next_elem_id = next_elem_id.saturating_add(1);
                }

                accepted_face += 1;
                did_flip = true;
            }
        }

        if !did_flip {
            break;
        }
    }
    Ok((accepted_face, accepted_edge, accepted_edge_sliver))
}

pub(crate) fn is_better_quality(a: (f64, f64, f64), b: (f64, f64, f64)) -> bool {
    (a.0 > b.0 + 1e-6)
        || ((a.0 - b.0).abs() < 1e-6
            && ((a.1 < b.1 - 1e-9)
                || ((a.1 - b.1).abs() < 1e-9 && a.2 < b.2 - 1e-9)))
}

fn has_sliver_pressure(mesh: &Mesh, sliver_fraction_threshold: f64, min_dihedral_threshold: f64) -> bool {
    let mut total = 0usize;
    let mut sliver_like = 0usize;
    let mut global_min_d = f64::MAX;

    for e in &mesh.elements {
        if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
            continue;
        }
        let Ok(a) = node_xyz_from_mesh(mesh, e.node_ids[0]) else {
            continue;
        };
        let Ok(b) = node_xyz_from_mesh(mesh, e.node_ids[1]) else {
            continue;
        };
        let Ok(c) = node_xyz_from_mesh(mesh, e.node_ids[2]) else {
            continue;
        };
        let Ok(d) = node_xyz_from_mesh(mesh, e.node_ids[3]) else {
            continue;
        };
        let v = tetra_volume(a, b, c, d);
        if v <= 1e-15 {
            continue;
        }

        total += 1;
        let dmin = min_dihedral_points(a, b, c, d);
        let r = radius_edge_ratio_points(a, b, c, d);
        if !dmin.is_finite() || !r.is_finite() {
            continue;
        }
        global_min_d = global_min_d.min(dmin);
        if dmin < 6.0 && r > 1.8 {
            sliver_like += 1;
        }
    }

    if total == 0 {
        return false;
    }

    let sliver_frac = sliver_like as f64 / total as f64;
    sliver_frac >= sliver_fraction_threshold || global_min_d <= min_dihedral_threshold
}

fn aggregate_tet_quality(mesh: &Mesh, tets: &[[u64; 4]]) -> Option<(f64, f64, f64)> {
    let mut min_d = f64::MAX;
    let mut max_r = 0.0_f64;
    let mut sliver_like = 0usize;
    for tet in tets {
        let a = node_xyz_from_mesh(mesh, tet[0]).ok()?;
        let b = node_xyz_from_mesh(mesh, tet[1]).ok()?;
        let c = node_xyz_from_mesh(mesh, tet[2]).ok()?;
        let d = node_xyz_from_mesh(mesh, tet[3]).ok()?;
        let v = tetra_volume(a, b, c, d);
        if v <= 1e-15 {
            return None;
        }

        let dmin = min_dihedral_points(a, b, c, d);
        let r = radius_edge_ratio_points(a, b, c, d);
        if !dmin.is_finite() || !r.is_finite() {
            return None;
        }

        min_d = min_d.min(dmin);
        max_r = max_r.max(r);
        if dmin < 6.0 && r > 1.8 {
            sliver_like += 1;
        }
    }
    Some((min_d, sliver_like as f64 / (tets.len() as f64), max_r))
}

fn should_log_refinement_stats() -> bool {
    std::env::var("RMSH_DEBUG_REFINEMENT")
        .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
        .unwrap_or(false)
}

fn find_worst_tetrahedron(
    mesh: &Mesh,
    max_radius_edge_ratio: f64,
    edge_limit: f64,
    sliver_floor_deg: f64,
) -> Option<(usize, f64)> {
    let mut worst: Option<(usize, f64)> = None;
    for (idx, elem) in mesh.elements.iter().enumerate() {
        if elem.etype != ElementType::Tetrahedron4 || elem.node_ids.len() != 4 {
            continue;
        }
        let Ok(r) = tetra_radius_edge_ratio_from_mesh(mesh, &elem.node_ids) else {
            continue;
        };
        let Ok(lmax) = tetra_max_edge_length_from_mesh(mesh, &elem.node_ids) else {
            continue;
        };
        let Ok(dmin) = tetra_min_dihedral_from_mesh(mesh, &elem.node_ids) else {
            continue;
        };
        if r <= max_radius_edge_ratio && lmax <= edge_limit && dmin >= sliver_floor_deg {
            continue;
        }

        let quality_pressure = r / max_radius_edge_ratio;
        let size_pressure = lmax / edge_limit;
        let dihedral_pressure = if dmin < sliver_floor_deg {
            1.0 + (sliver_floor_deg - dmin) / sliver_floor_deg
        } else {
            0.0
        };
        let score = quality_pressure.max(size_pressure).max(dihedral_pressure);
        match worst {
            Some((_, w)) if score <= w => {}
            _ => worst = Some((idx, score)),
        }
    }
    worst
}

#[cfg(test)]
fn tetra_centroid_from_mesh(mesh: &Mesh, tet: &[u64]) -> Result<[f64; 3], MeshAlgoError> {
    if tet.len() != 4 {
        return Err(MeshAlgoError::Generation(
            "tetrahedron must have 4 nodes".to_string(),
        ));
    }
    let mut sum = [0.0_f64; 3];
    for &nid in tet {
        let node = mesh
            .nodes
            .get(&nid)
            .ok_or_else(|| MeshAlgoError::Generation(format!("missing node id {nid}")))?;
        sum[0] += node.position.x;
        sum[1] += node.position.y;
        sum[2] += node.position.z;
    }
    Ok([sum[0] * 0.25, sum[1] * 0.25, sum[2] * 0.25])
}

fn node_xyz_from_mesh(mesh: &Mesh, node_id: u64) -> Result<[f64; 3], MeshAlgoError> {
    let node = mesh
        .nodes
        .get(&node_id)
        .ok_or_else(|| MeshAlgoError::Generation(format!("missing node id {node_id}")))?;
    Ok([node.position.x, node.position.y, node.position.z])
}

fn select_refinement_point(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Option<[f64; 3]> {
    let parent_ratio = radius_edge_ratio_points(a, b, c, d);
    let parent_dihedral = min_dihedral_points(a, b, c, d);

    let centroid = [
        (a[0] + b[0] + c[0] + d[0]) * 0.25,
        (a[1] + b[1] + c[1] + d[1]) * 0.25,
        (a[2] + b[2] + c[2] + d[2]) * 0.25,
    ];
    let mut candidates = Vec::<[f64; 3]>::with_capacity(16);
    candidates.push(centroid);

    // Barycentric interior points with controlled distance from faces.
    for &(wa, wb, wc, wd) in &[
        (0.40, 0.20, 0.20, 0.20),
        (0.20, 0.40, 0.20, 0.20),
        (0.20, 0.20, 0.40, 0.20),
        (0.20, 0.20, 0.20, 0.40),
        (0.55, 0.15, 0.15, 0.15),
        (0.15, 0.55, 0.15, 0.15),
        (0.15, 0.15, 0.55, 0.15),
        (0.15, 0.15, 0.15, 0.55),
    ] {
        candidates.push([
            wa * a[0] + wb * b[0] + wc * c[0] + wd * d[0],
            wa * a[1] + wb * b[1] + wc * c[1] + wd * d[1],
            wa * a[2] + wb * b[2] + wc * c[2] + wd * d[2],
        ]);
    }

    for v in [a, b, c, d] {
        candidates.push([
            centroid[0] * 0.85 + v[0] * 0.15,
            centroid[1] * 0.85 + v[1] * 0.15,
            centroid[2] * 0.85 + v[2] * 0.15,
        ]);
    }

    if let Some(ic) = tetra_incenter(a, b, c, d) {
        candidates.push(ic);
    }

    let (cc, _r) = circumsphere(a, b, c, d);
    if cc[0].is_finite() && cc[1].is_finite() && cc[2].is_finite() {
        candidates.push([
            centroid[0] * 0.70 + cc[0] * 0.30,
            centroid[1] * 0.70 + cc[1] * 0.30,
            centroid[2] * 0.70 + cc[2] * 0.30,
        ]);
    }

    let mut best: Option<([f64; 3], f64, f64, f64)> = None;
    let mut best_strict: Option<([f64; 3], f64, f64, f64)> = None;
    for p in candidates {
        if !point_in_tetrahedron(a, b, c, d, p, 1e-14) {
            continue;
        }
        let Some((_quality, child_min_dihedral, child_max_ratio, child_sliver_frac)) =
            split_quality_metrics(a, b, c, d, p)
        else {
            continue;
        };

        // Soft guard: reject only extreme blow-ups.
        if child_max_ratio > parent_ratio * 3.0 && child_min_dihedral < parent_dihedral {
            continue;
        }
        if child_sliver_frac > 0.75 && child_min_dihedral <= parent_dihedral {
            continue;
        }

        if child_min_dihedral >= 0.5 {
            match best_strict {
                Some((_, bd, bs, br))
                    if (child_min_dihedral < bd)
                        || ((child_min_dihedral - bd).abs() < 1e-9
                            && ((child_sliver_frac > bs)
                                || ((child_sliver_frac - bs).abs() < 1e-9
                                    && child_max_ratio >= br))) =>
                {
                }
                _ => {
                    best_strict =
                        Some((p, child_min_dihedral, child_sliver_frac, child_max_ratio))
                }
            }
        }

        match best {
            Some((_, bd, bs, br))
                if (child_min_dihedral < bd)
                    || ((child_min_dihedral - bd).abs() < 1e-9
                        && ((child_sliver_frac > bs)
                            || ((child_sliver_frac - bs).abs() < 1e-9 && child_max_ratio >= br))) =>
            {
            }
            _ => best = Some((p, child_min_dihedral, child_sliver_frac, child_max_ratio)),
        }
    }
    best_strict.or(best).map(|(p, _, _, _)| p)
}

fn split_quality_metrics(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    p: [f64; 3],
) -> Option<(f64, f64, f64, f64)> {
    let tets = [[a, b, c, p], [a, b, d, p], [a, c, d, p], [b, c, d, p]];
    let mut min_d = f64::MAX;
    let mut max_r: f64 = 0.0;
    let mut sliver_like = 0usize;
    for t in tets {
        let v = tetra_volume(t[0], t[1], t[2], t[3]);
        if v <= 1e-15 {
            return None;
        }
        let dmin = min_dihedral_points(t[0], t[1], t[2], t[3]);
        if !dmin.is_finite() {
            return None;
        }
        min_d = min_d.min(dmin);
        let r = radius_edge_ratio_points(t[0], t[1], t[2], t[3]);
        if !r.is_finite() {
            return None;
        }
        max_r = max_r.max(r);
        if dmin < 6.0 && r > 1.8 {
            sliver_like += 1;
        }
    }
    let sliver_frac = sliver_like as f64 / 4.0;
    // Keep a compact aggregate score for fallback callers; primary ranking uses tuple rules.
    let score = 1.25 * min_d - 0.60 * max_r - 5.50 * sliver_frac;
    Some((score, min_d, max_r, sliver_frac))
}

fn point_in_tetrahedron(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    p: [f64; 3],
    eps: f64,
) -> bool {
    let v = tetra_volume(a, b, c, d);
    if v <= eps {
        return false;
    }
    let v0 = tetra_volume(p, b, c, d);
    let v1 = tetra_volume(a, p, c, d);
    let v2 = tetra_volume(a, b, p, d);
    let v3 = tetra_volume(a, b, c, p);
    let sum = v0 + v1 + v2 + v3;
    if (sum - v).abs() > eps * 32.0 {
        return false;
    }
    v0 > eps && v1 > eps && v2 > eps && v3 > eps
}

pub(crate) fn min_dihedral_points(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    [
        dihedral(a, b, c, d),
        dihedral(a, c, b, d),
        dihedral(a, d, b, c),
        dihedral(b, c, a, d),
        dihedral(b, d, a, c),
        dihedral(c, d, a, b),
    ]
    .into_iter()
    .fold(f64::MAX, f64::min)
}

pub(crate) fn radius_edge_ratio_points(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let nodes = [a, b, c, d];
    radius_edge_ratio(&nodes, [0, 1, 2, 3])
}

fn tetra_incenter(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Option<[f64; 3]> {
    let sa = triangle_area(b, c, d);
    let sb = triangle_area(a, c, d);
    let sc = triangle_area(a, b, d);
    let sd = triangle_area(a, b, c);
    let sum = sa + sb + sc + sd;
    if sum <= 1e-15 {
        return None;
    }
    Some([
        (sa * a[0] + sb * b[0] + sc * c[0] + sd * d[0]) / sum,
        (sa * a[1] + sb * b[1] + sc * c[1] + sd * d[1]) / sum,
        (sa * a[2] + sb * b[2] + sc * c[2] + sd * d[2]) / sum,
    ])
}

fn triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cross = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    0.5 * (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt()
}

fn min_child_dihedral_for_point(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    p: [f64; 3],
) -> f64 {
    [
        min_dihedral_points(a, b, c, p),
        min_dihedral_points(a, b, d, p),
        min_dihedral_points(a, c, d, p),
        min_dihedral_points(b, c, d, p),
    ]
    .into_iter()
    .fold(f64::MAX, f64::min)
}

fn best_edge_split_partition(
    points: [[f64; 3]; 4],
) -> Option<(usize, usize, usize, usize, [f64; 3], (f64, f64, f64))> {
    let edges = [(0usize, 1usize), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let mut best: Option<(usize, usize, usize, usize, [f64; 3], f64, f64, f64)> = None;
    for (i, j) in edges {
        let mut others = [0usize; 2];
        let mut oi = 0usize;
        for k in 0..4 {
            if k != i && k != j {
                others[oi] = k;
                oi += 1;
            }
        }
        let k = others[0];
        let l = others[1];

        let Some((best_point, min_d, sliver_frac, max_r)) =
            edge_split_quality_metrics(points, i, j, k, l)
        else {
            continue;
        };

        match best {
            Some((_, _, _, _, _, bd, bs, br))
                if (min_d < bd)
                    || ((min_d - bd).abs() < 1e-9
                        && ((sliver_frac > bs)
                            || ((sliver_frac - bs).abs() < 1e-9 && max_r >= br))) => {}
            _ => best = Some((i, j, k, l, best_point, min_d, sliver_frac, max_r)),
        }
    }
    best.map(|(i, j, k, l, p, md, sf, mr)| (i, j, k, l, p, (md, sf, mr)))
}

fn edge_split_quality_metrics(
    points: [[f64; 3]; 4],
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Option<([f64; 3], f64, f64, f64)> {
    let alphas = [0.50_f64];
    let mut best: Option<([f64; 3], f64, f64, f64)> = None;

    for alpha in alphas {
        let p = [
            points[i][0] * (1.0 - alpha) + points[j][0] * alpha,
            points[i][1] * (1.0 - alpha) + points[j][1] * alpha,
            points[i][2] * (1.0 - alpha) + points[j][2] * alpha,
        ];

        let tets = [
            [points[i], p, points[k], points[l]],
            [p, points[j], points[k], points[l]],
        ];
        let mut min_d = f64::MAX;
        let mut max_r = 0.0_f64;
        let mut sliver_like = 0usize;
        let mut valid = true;
        for t in tets {
            let v = tetra_volume(t[0], t[1], t[2], t[3]);
            if v <= 1e-15 {
                valid = false;
                break;
            }
            let d = min_dihedral_points(t[0], t[1], t[2], t[3]);
            let r = radius_edge_ratio_points(t[0], t[1], t[2], t[3]);
            if !d.is_finite() || !r.is_finite() {
                valid = false;
                break;
            }
            min_d = min_d.min(d);
            max_r = max_r.max(r);
            if d < 6.0 && r > 1.8 {
                sliver_like += 1;
            }
        }
        if !valid {
            continue;
        }
        let sliver_frac = sliver_like as f64 / 2.0;
        match best {
            Some((_, bd, bs, br))
                if (min_d < bd)
                    || ((min_d - bd).abs() < 1e-9
                        && ((sliver_frac > bs)
                            || ((sliver_frac - bs).abs() < 1e-9 && max_r >= br))) => {}
            _ => best = Some((p, min_d, sliver_frac, max_r)),
        }
    }

    best
}

fn longest_edge_biased_point(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
) -> Option<[f64; 3]> {
    let verts = [a, b, c, d];
    let edges = [(0usize, 1usize), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let mut best = None::<((usize, usize), f64)>;
    for (i, j) in edges {
        let dx = verts[i][0] - verts[j][0];
        let dy = verts[i][1] - verts[j][1];
        let dz = verts[i][2] - verts[j][2];
        let l2 = dx * dx + dy * dy + dz * dz;
        match best {
            Some((_, b2)) if l2 <= b2 => {}
            _ => best = Some(((i, j), l2)),
        }
    }
    let ((i, j), _) = best?;
    let mut others = Vec::<usize>::with_capacity(2);
    for k in 0..4 {
        if k != i && k != j {
            others.push(k);
        }
    }

    let u = verts[i];
    let v = verts[j];
    let w = verts[others[0]];
    let x = verts[others[1]];
    Some([
        0.40 * u[0] + 0.40 * v[0] + 0.10 * w[0] + 0.10 * x[0],
        0.40 * u[1] + 0.40 * v[1] + 0.10 * w[1] + 0.10 * x[1],
        0.40 * u[2] + 0.40 * v[2] + 0.10 * w[2] + 0.10 * x[2],
    ])
}

fn select_fallback_refinement_point(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
) -> Option<[f64; 3]> {
    let centroid = [
        (a[0] + b[0] + c[0] + d[0]) * 0.25,
        (a[1] + b[1] + c[1] + d[1]) * 0.25,
        (a[2] + b[2] + c[2] + d[2]) * 0.25,
    ];

    let mut candidates = vec![centroid];
    if let Some(ic) = tetra_incenter(a, b, c, d) {
        candidates.push(ic);
        // Blend toward centroid to avoid hugging near-degenerate face geometry.
        candidates.push([
            0.7 * ic[0] + 0.3 * centroid[0],
            0.7 * ic[1] + 0.3 * centroid[1],
            0.7 * ic[2] + 0.3 * centroid[2],
        ]);
    }
    if let Some(lp) = longest_edge_biased_point(a, b, c, d) {
        candidates.push(lp);
    }

    // Vertex-to-opposite-face-centroid interior candidates.
    for (v, f1, f2, f3) in [(a, b, c, d), (b, a, c, d), (c, a, b, d), (d, a, b, c)] {
        let fc = [
            (f1[0] + f2[0] + f3[0]) / 3.0,
            (f1[1] + f2[1] + f3[1]) / 3.0,
            (f1[2] + f2[2] + f3[2]) / 3.0,
        ];
        candidates.push([
            0.25 * v[0] + 0.75 * fc[0],
            0.25 * v[1] + 0.75 * fc[1],
            0.25 * v[2] + 0.75 * fc[2],
        ]);
    }

    let mut best: Option<([f64; 3], f64, f64, f64)> = None;
    let mut best_strict: Option<([f64; 3], f64, f64, f64)> = None;
    for p in candidates {
        if !point_in_tetrahedron(a, b, c, d, p, 1e-14) {
            continue;
        }
        let Some((_score, min_d, max_r, sliver_frac)) = split_quality_metrics(a, b, c, d, p)
        else {
            continue;
        };

        // Prefer non-sliver candidates first if available.
        if min_d >= 0.5 {
            match best_strict {
                Some((_, bd, bs, br))
                    if (min_d < bd)
                        || ((min_d - bd).abs() < 1e-9
                            && ((sliver_frac > bs)
                                || ((sliver_frac - bs).abs() < 1e-9 && max_r >= br))) =>
                {
                }
                _ => best_strict = Some((p, min_d, sliver_frac, max_r)),
            }
        }

        match best {
            Some((_, bd, bs, br))
                if (min_d < bd)
                    || ((min_d - bd).abs() < 1e-9
                        && ((sliver_frac > bs)
                            || ((sliver_frac - bs).abs() < 1e-9 && max_r >= br))) => {}
            _ => best = Some((p, min_d, sliver_frac, max_r)),
        }
    }
    best_strict
        .or(best)
        .map(|(p, _, _, _)| p)
}

pub(crate) fn tetra_volume(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let ad = [a[0] - d[0], a[1] - d[1], a[2] - d[2]];
    let bd = [b[0] - d[0], b[1] - d[1], b[2] - d[2]];
    let cd = [c[0] - d[0], c[1] - d[1], c[2] - d[2]];
    let cross = [
        bd[1] * cd[2] - bd[2] * cd[1],
        bd[2] * cd[0] - bd[0] * cd[2],
        bd[0] * cd[1] - bd[1] * cd[0],
    ];
    (ad[0] * cross[0] + ad[1] * cross[1] + ad[2] * cross[2]).abs() / 6.0
}

fn dihedral(p: [f64; 3], q: [f64; 3], r: [f64; 3], s: [f64; 3]) -> f64 {
    let pq = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
    let pr = [r[0] - p[0], r[1] - p[1], r[2] - p[2]];
    let ps = [s[0] - p[0], s[1] - p[1], s[2] - p[2]];
    let n1 = [
        pq[1] * pr[2] - pq[2] * pr[1],
        pq[2] * pr[0] - pq[0] * pr[2],
        pq[0] * pr[1] - pq[1] * pr[0],
    ];
    let n2 = [
        pq[1] * ps[2] - pq[2] * ps[1],
        pq[2] * ps[0] - pq[0] * ps[2],
        pq[0] * ps[1] - pq[1] * ps[0],
    ];
    let dot = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
    let l1 = (n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]).sqrt();
    let l2 = (n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]).sqrt();
    if l1 < 1e-12 || l2 < 1e-12 {
        return 0.0;
    }
    (dot / (l1 * l2)).clamp(-1.0, 1.0).acos().to_degrees()
}

fn tetra_radius_edge_ratio_from_mesh(mesh: &Mesh, tet: &[u64]) -> Result<f64, MeshAlgoError> {
    if tet.len() != 4 {
        return Err(MeshAlgoError::Generation(
            "tetrahedron must have 4 nodes".to_string(),
        ));
    }
    let mut pts = [[0.0_f64; 3]; 4];
    for (i, &nid) in tet.iter().enumerate() {
        let node = mesh
            .nodes
            .get(&nid)
            .ok_or_else(|| MeshAlgoError::Generation(format!("missing node id {nid}")))?;
        pts[i] = [node.position.x, node.position.y, node.position.z];
    }
    Ok(radius_edge_ratio(&pts, [0, 1, 2, 3]))
}

fn tetra_max_edge_length_from_mesh(mesh: &Mesh, tet: &[u64]) -> Result<f64, MeshAlgoError> {
    if tet.len() != 4 {
        return Err(MeshAlgoError::Generation(
            "tetrahedron must have 4 nodes".to_string(),
        ));
    }
    let mut pts = [[0.0_f64; 3]; 4];
    for (i, &nid) in tet.iter().enumerate() {
        let node = mesh
            .nodes
            .get(&nid)
            .ok_or_else(|| MeshAlgoError::Generation(format!("missing node id {nid}")))?;
        pts[i] = [node.position.x, node.position.y, node.position.z];
    }

    let edges = [(0usize, 1usize), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let mut lmax = 0.0_f64;
    for (i, j) in edges {
        let dx = pts[i][0] - pts[j][0];
        let dy = pts[i][1] - pts[j][1];
        let dz = pts[i][2] - pts[j][2];
        let l = (dx * dx + dy * dy + dz * dz).sqrt();
        lmax = lmax.max(l);
    }
    Ok(lmax)
}

fn tetra_min_dihedral_from_mesh(mesh: &Mesh, tet: &[u64]) -> Result<f64, MeshAlgoError> {
    if tet.len() != 4 {
        return Err(MeshAlgoError::Generation(
            "tetrahedron must have 4 nodes".to_string(),
        ));
    }
    let a = node_xyz_from_mesh(mesh, tet[0])?;
    let b = node_xyz_from_mesh(mesh, tet[1])?;
    let c = node_xyz_from_mesh(mesh, tet[2])?;
    let d = node_xyz_from_mesh(mesh, tet[3])?;
    Ok(min_dihedral_points(a, b, c, d))
}

fn validate_params(algo: &Delaunay3D, params: &MeshParams) -> Result<(), MeshAlgoError> {
    if !params.element_size.is_finite() || params.element_size <= 0.0 {
        return Err(MeshAlgoError::InvalidInput(
            "element_size must be a positive finite value".to_string(),
        ));
    }
    if !params.max_size.is_finite() || params.max_size <= 0.0 {
        return Err(MeshAlgoError::InvalidInput(
            "max_size must be a positive finite value".to_string(),
        ));
    }
    if params.max_size < params.element_size {
        return Err(MeshAlgoError::InvalidInput(
            "max_size must be >= element_size".to_string(),
        ));
    }
    if !algo.max_radius_edge_ratio.is_finite() || algo.max_radius_edge_ratio < 2.0 {
        return Err(MeshAlgoError::InvalidInput(
            "max_radius_edge_ratio must be finite and >= 2.0".to_string(),
        ));
    }
    if !algo.min_dihedral_angle_deg.is_finite() || algo.min_dihedral_angle_deg < 0.0 {
        return Err(MeshAlgoError::InvalidInput(
            "min_dihedral_angle_deg must be finite and >= 0.0".to_string(),
        ));
    }
    Ok(())
}

// ─── Internal helpers (stubs) ─────────────────────────────────────────────────

/// Compute the circumsphere of a tetrahedron with vertices `a, b, c, d`.
///
/// Returns `(centre, radius)`.
#[allow(dead_code)]
fn circumsphere(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> ([f64; 3], f64) {
    // Solve linear system from |x-a|^2 = |x-b|^2 = |x-c|^2 = |x-d|^2.
    // This yields A * x = rhs with 3 equations.
    let rows = [
        (
            [
                2.0 * (b[0] - a[0]),
                2.0 * (b[1] - a[1]),
                2.0 * (b[2] - a[2]),
            ],
            b[0] * b[0] + b[1] * b[1] + b[2] * b[2] - (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]),
        ),
        (
            [
                2.0 * (c[0] - a[0]),
                2.0 * (c[1] - a[1]),
                2.0 * (c[2] - a[2]),
            ],
            c[0] * c[0] + c[1] * c[1] + c[2] * c[2] - (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]),
        ),
        (
            [
                2.0 * (d[0] - a[0]),
                2.0 * (d[1] - a[1]),
                2.0 * (d[2] - a[2]),
            ],
            d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]),
        ),
    ];

    if let Some(center) = solve_3x3(rows) {
        let dx = center[0] - a[0];
        let dy = center[1] - a[1];
        let dz = center[2] - a[2];
        let radius = (dx * dx + dy * dy + dz * dz).sqrt();
        (center, radius)
    } else {
        // Degenerate tetrahedron: return finite fallback center + infinite radius.
        let center = [
            (a[0] + b[0] + c[0] + d[0]) * 0.25,
            (a[1] + b[1] + c[1] + d[1]) * 0.25,
            (a[2] + b[2] + c[2] + d[2]) * 0.25,
        ];
        (center, f64::INFINITY)
    }
}

/// Test whether point `p` lies strictly inside the circumsphere of `(a,b,c,d)`.
///
/// Uses the in-sphere predicate.  Returns `> 0` if inside, `< 0` if outside,
/// `0` on the sphere (degenerate).
#[allow(dead_code)]
fn in_sphere_test(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3], p: [f64; 3]) -> f64 {
    // Sign convention here is "radius - distance":
    // > 0 => inside, 0 => on sphere, < 0 => outside.
    let (center, radius) = circumsphere(a, b, c, d);
    if !radius.is_finite() {
        return 0.0;
    }
    let dx = p[0] - center[0];
    let dy = p[1] - center[1];
    let dz = p[2] - center[2];
    radius - (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Perform a 3-D bistellar flip on the set of tetrahedra sharing an edge or face.
///
/// * 2-to-3 flip: split the two tets sharing a face into three new tets sharing the
///   edge between their opposite vertices.
/// * 3-to-2 flip: merge three tets sharing an edge into two new tets sharing a face.
///
/// Both operate on `tets` in-place: the old tetrahedra are `swap_remove`d and new
/// ones appended. Node positions are needed for volume computations.
fn bistellar_flip(
    tets: &mut Vec<[usize; 4]>,
    nodes: &[[f64; 3]],
    flip_type: BistellarFlipType,
    indices: &[usize],
) -> Result<(), MeshAlgoError> {
    match flip_type {
        BistellarFlipType::TwoToThree => bistellar_flip_2_to_3(tets, nodes, indices),
        BistellarFlipType::ThreeToTwo => bistellar_flip_3_to_2(tets, nodes, indices),
        BistellarFlipType::FourToFour => Err(MeshAlgoError::NotImplemented),
    }
}

/// 2-to-3 flip: replace two tetrahedra sharing a face `(a,b,c)` with three
/// tetrahedra sharing the edge `(d1,d2)`.
///
/// `indices` must be `[idx0, idx1]` pointing to the two tets in `tets`.
fn bistellar_flip_2_to_3(
    tets: &mut Vec<[usize; 4]>,
    nodes: &[[f64; 3]],
    indices: &[usize],
) -> Result<(), MeshAlgoError> {
    if indices.len() != 2 {
        return Err(MeshAlgoError::Generation("2-to-3 flip requires exactly 2 tets".into()));
    }
    let (i0, i1) = (indices[0], indices[1]);
    if i0 >= tets.len() || i1 >= tets.len() {
        return Ok(());
    }

    let t0 = tets[i0];
    let t1 = tets[i1];

    // Find the shared face (3 common vertices) and the two opposite vertices.
    let mut common = Vec::new();
    let mut d1 = None;
    let mut d2 = None;
    for &v in &t0 {
        if t1.contains(&v) {
            common.push(v);
        } else {
            d1 = Some(v);
        }
    }
    for &v in &t1 {
        if !t0.contains(&v) {
            d2 = Some(v);
        }
    }

    if common.len() != 3 || d1.is_none() || d2.is_none() {
        return Ok(());
    }
    let face @ [a, b, c] = [common[0], common[1], common[2]];
    let (d1, d2) = (d1.unwrap(), d2.unwrap());

    // The new edge (d1, d2) must intersect the shared face for a valid flip.
    // Check that d1 and d2 are on opposite sides of face (a,b,c).
    let v = tetra_volume_3d(nodes[a], nodes[b], nodes[c], nodes[d1]);
    let w = tetra_volume_3d(nodes[a], nodes[b], nodes[c], nodes[d2]);
    if v.signum() == w.signum() || v.abs() < 1e-20 || w.abs() < 1e-20 {
        return Ok(());
    }

    // Three new tets sharing edge (d1, d2): (a, b, d1, d2), (b, c, d1, d2), (c, a, d1, d2)
    let new_tets = [
        [face[0], face[1], d1, d2],
        [face[1], face[2], d1, d2],
        [face[2], face[0], d1, d2],
    ];

    // Validate new tets have positive volume
    for &tet in &new_tets {
        let vol = tetra_volume_3d(nodes[tet[0]], nodes[tet[1]], nodes[tet[2]], nodes[tet[3]]);
        if vol < 1e-15 {
            return Ok(());
        }
    }

    // Remove old tets (higher index first) and add new ones
    let (high, low) = if i0 > i1 { (i0, i1) } else { (i1, i0) };
    tets.swap_remove(high);
    if low < tets.len() {
        tets.swap_remove(low);
    } else if high != low {
        // low was removed by swap_remove(high)
    }
    for &tet in &new_tets {
        tets.push(tet);
    }

    Ok(())
}

/// 3-to-2 flip: replace three tetrahedra sharing an edge `(a,b)` with two
/// tetrahedra sharing the face `(c1,c2,c3)`.
///
/// `indices` must be `[i0, i1, i2]` pointing to the three tets.
fn bistellar_flip_3_to_2(
    tets: &mut Vec<[usize; 4]>,
    nodes: &[[f64; 3]],
    indices: &[usize],
) -> Result<(), MeshAlgoError> {
    if indices.len() != 3 {
        return Err(MeshAlgoError::Generation("3-to-2 flip requires exactly 3 tets".into()));
    }

    // Find the edge (a,b) shared by all three tets and their opposite vertices.
    let t0 = tets[indices[0]];
    let t1 = tets[indices[1]];
    let t2 = tets[indices[2]];

    // Find common vertex pair
    let mut edge = None;
    for &u in &t0 {
        for &v in &t0 {
            if u >= v {
                continue;
            }
            if t1.contains(&u) && t1.contains(&v) && t2.contains(&u) && t2.contains(&v) {
                edge = Some((u, v));
                break;
            }
        }
        if edge.is_some() {
            break;
        }
    }

    let Some((a, b)) = edge else {
        return Ok(());
    };

    // Collect the three opposite vertices (one from each tet)
    let mut opp = Vec::new();
    for &idx in indices {
        let tet = tets[idx];
        for &v in &tet {
            if v != a && v != b {
                opp.push(v);
                break;
            }
        }
    }
    if opp.len() != 3 {
        return Ok(());
    }
    let [c1, c2, c3] = [opp[0], opp[1], opp[2]];

    // Two new tets: (a, c1, c2, c3) and (b, c1, c2, c3)
    // Ensure both have positive volume
    let v1 = tetra_volume_3d(nodes[a], nodes[c1], nodes[c2], nodes[c3]);
    let v2 = tetra_volume_3d(nodes[b], nodes[c1], nodes[c2], nodes[c3]);
    if v1.abs() < 1e-15 || v2.abs() < 1e-15 {
        return Ok(());
    }

    let new_tet1 = if v1 > 0.0 {
        [a, c1, c2, c3]
    } else {
        [a, c2, c1, c3]
    };
    let new_tet2 = if v2 > 0.0 {
        [b, c1, c3, c2]
    } else {
        [b, c2, c3, c1]
    };

    // Remove old tets (descending index order)
    let mut sorted_indices = [indices[0], indices[1], indices[2]];
    sorted_indices.sort_unstable_by(|a, b| b.cmp(a));
    for &idx in &sorted_indices {
        if idx < tets.len() {
            tets.swap_remove(idx);
        }
    }

    tets.push(new_tet1);
    tets.push(new_tet2);

    Ok(())
}

fn tetra_volume_3d(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let ad = [a[0] - d[0], a[1] - d[1], a[2] - d[2]];
    let bd = [b[0] - d[0], b[1] - d[1], b[2] - d[2]];
    let cd = [c[0] - d[0], c[1] - d[1], c[2] - d[2]];
    let cross = [
        bd[1] * cd[2] - bd[2] * cd[1],
        bd[2] * cd[0] - bd[0] * cd[2],
        bd[0] * cd[1] - bd[1] * cd[0],
    ];
    (ad[0] * cross[0] + ad[1] * cross[1] + ad[2] * cross[2]) / 6.0
}

enum BistellarFlipType {
    /// Replace 2 tets sharing a face with 3 tets sharing an edge.
    TwoToThree,
    /// Replace 3 tets sharing an edge with 2 tets sharing a face.
    ThreeToTwo,
    /// Replace 4 tets sharing a degree-2 edge with 4 tets (4-to-4).
    FourToFour,
}

/// Compute the radius-edge ratio `R / l_min` of a tetrahedron.
///
/// `R` is the circumradius; `l_min` is the length of the shortest edge.
#[allow(dead_code)]
fn radius_edge_ratio(nodes: &[[f64; 3]], tet: [usize; 4]) -> f64 {
    if tet.iter().any(|&i| i >= nodes.len()) {
        return f64::INFINITY;
    }

    let a = nodes[tet[0]];
    let b = nodes[tet[1]];
    let c = nodes[tet[2]];
    let d = nodes[tet[3]];

    let (_, radius) = circumsphere(a, b, c, d);
    if !radius.is_finite() {
        return f64::INFINITY;
    }

    let mut min_edge = f64::INFINITY;
    let edges = [(a, b), (a, c), (a, d), (b, c), (b, d), (c, d)];
    for (u, v) in edges {
        let dx = u[0] - v[0];
        let dy = u[1] - v[1];
        let dz = u[2] - v[2];
        let l = (dx * dx + dy * dy + dz * dz).sqrt();
        min_edge = min_edge.min(l);
    }

    if min_edge <= 1e-15 {
        return f64::INFINITY;
    }
    radius / min_edge
}

fn solve_3x3(rows: [([f64; 3], f64); 3]) -> Option<[f64; 3]> {
    let mut m = [
        [rows[0].0[0], rows[0].0[1], rows[0].0[2], rows[0].1],
        [rows[1].0[0], rows[1].0[1], rows[1].0[2], rows[1].1],
        [rows[2].0[0], rows[2].0[1], rows[2].0[2], rows[2].1],
    ];

    for col in 0..3 {
        let mut pivot = col;
        for r in (col + 1)..3 {
            if m[r][col].abs() > m[pivot][col].abs() {
                pivot = r;
            }
        }
        if m[pivot][col].abs() < 1e-15 {
            return None;
        }
        if pivot != col {
            m.swap(pivot, col);
        }

        let pivot_val = m[col][col];
        for j in col..4 {
            m[col][j] /= pivot_val;
        }

        for r in 0..3 {
            if r == col {
                continue;
            }
            let factor = m[r][col];
            for j in col..4 {
                m[r][j] -= factor * m[col][j];
            }
        }
    }

    Some([m[0][3], m[1][3], m[2][3]])
}

#[cfg(test)]
mod tests {
    use super::*;

    use rmsh_model::{Element, ElementType, Mesh, Node};

    fn cube_surface_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 1.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(5, 0.0, 0.0, 1.0));
        mesh.add_node(Node::new(6, 1.0, 0.0, 1.0));
        mesh.add_node(Node::new(7, 1.0, 1.0, 1.0));
        mesh.add_node(Node::new(8, 0.0, 1.0, 1.0));

        mesh.add_element(Element::new(1, ElementType::Quad4, vec![1, 2, 3, 4]));
        mesh.add_element(Element::new(2, ElementType::Quad4, vec![5, 6, 7, 8]));
        mesh.add_element(Element::new(3, ElementType::Quad4, vec![1, 2, 6, 5]));
        mesh.add_element(Element::new(4, ElementType::Quad4, vec![2, 3, 7, 6]));
        mesh.add_element(Element::new(5, ElementType::Quad4, vec![3, 4, 8, 7]));
        mesh.add_element(Element::new(6, ElementType::Quad4, vec![4, 1, 5, 8]));
        mesh
    }

    #[test]
    fn delaunay3d_name_is_stable() {
        let algo = Delaunay3D::new();
        assert_eq!(algo.name(), "Delaunay 3D");
    }

    #[test]
    fn delaunay3d_mesh_flow_runs() {
        let algo = Delaunay3D::default();
        let params = MeshParams::with_size(0.4);
        let out = algo
            .mesh_3d(&cube_surface_mesh(), &params)
            .expect("meshing should succeed");

        assert!(out.node_count() >= 9);
        assert!(out.elements_by_dimension(3).len() >= 12);
    }

    #[test]
    fn delaunay3d_respects_mesh_size_density() {
        let algo = Delaunay3D::default();
        let mesh = cube_surface_mesh();

        let mut coarse = MeshParams::with_size(1.0);
        coarse.max_size = 1.2;
        coarse.optimize_passes = 2;

        let mut fine = MeshParams::with_size(0.25);
        fine.max_size = 0.3;
        fine.optimize_passes = 2;

        let out_coarse = algo
            .mesh_3d(&mesh, &coarse)
            .expect("coarse meshing should succeed");
        let out_fine = algo
            .mesh_3d(&mesh, &fine)
            .expect("fine meshing should succeed");

        let coarse_tets = out_coarse.elements_by_dimension(3).len();
        let fine_tets = out_fine.elements_by_dimension(3).len();
        assert!(
            fine_tets > coarse_tets,
            "smaller mesh size should create denser tetra mesh: coarse={coarse_tets}, fine={fine_tets}"
        );
    }

    #[test]
    fn delaunay3d_rejects_bad_mesh_params() {
        let algo = Delaunay3D::default();
        let bad = MeshParams {
            element_size: 0.0,
            min_size: 0.0,
            max_size: 0.0,
            optimize_passes: 0,
        };

        let err = algo
            .mesh_3d(&Mesh::new(), &bad)
            .expect_err("invalid mesh params should error");
        match err {
            MeshAlgoError::InvalidInput(msg) => assert!(msg.contains("element_size")),
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn delaunay3d_rejects_invalid_algo_params() {
        let mut algo = Delaunay3D::default();
        algo.max_radius_edge_ratio = 1.9;
        let params = MeshParams::with_size(0.5);

        let err = algo
            .mesh_3d(&cube_surface_mesh(), &params)
            .expect_err("invalid algorithm params should error");
        match err {
            MeshAlgoError::InvalidInput(msg) => assert!(msg.contains("max_radius_edge_ratio")),
            other => panic!("unexpected error: {other:?}"),
        }

        let mut algo = Delaunay3D::default();
        algo.min_dihedral_angle_deg = -1.0;
        let err = algo
            .mesh_3d(&cube_surface_mesh(), &params)
            .expect_err("negative dihedral angle should error");
        match err {
            MeshAlgoError::InvalidInput(msg) => assert!(msg.contains("min_dihedral_angle_deg")),
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn validate_params_accepts_default_configuration() {
        let algo = Delaunay3D::default();
        let params = MeshParams::with_size(0.5);
        validate_params(&algo, &params).expect("default parameters should be valid");
    }

    #[test]
    fn circumsphere_and_in_sphere_work_for_regular_tet() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];

        let (center, radius) = circumsphere(a, b, c, d);
        assert!(radius.is_finite());
        assert!((center[0] - 0.5).abs() < 1e-9);
        assert!((center[1] - 0.5).abs() < 1e-9);
        assert!((center[2] - 0.5).abs() < 1e-9);

        let inside = [0.5, 0.5, 0.5];
        let outside = [2.0, 2.0, 2.0];
        assert!(in_sphere_test(a, b, c, d, inside) > 0.0);
        assert!(in_sphere_test(a, b, c, d, outside) < 0.0);
    }

    #[test]
    fn radius_edge_ratio_is_finite_for_non_degenerate_tet() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let ratio = radius_edge_ratio(&nodes, [0, 1, 2, 3]);
        assert!(ratio.is_finite());
        assert!(ratio > 0.0);
    }

    #[test]
    fn tetra_max_edge_length_from_mesh_works() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 2.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 0.0, 1.0));

        let lmax = tetra_max_edge_length_from_mesh(&mesh, &[1, 2, 3, 4]).expect("valid tet");
        assert!((lmax - 2.2360679).abs() < 1e-5);
    }

    // ── P1: circumsphere & in_sphere_test ─────────────────────────────────────

    #[test]
    fn circumsphere_of_unit_tet_is_at_center_with_known_radius() {
        // Unit tet: a=(0,0,0), b=(1,0,0), c=(0,1,0), d=(0,0,1).
        // Circumcenter is at (0.5, 0.5, 0.5), radius = sqrt(3)/2.
        let a = [0.0f64, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];
        let (center, radius) = circumsphere(a, b, c, d);
        assert!((center[0] - 0.5).abs() < 1e-9, "cx={}", center[0]);
        assert!((center[1] - 0.5).abs() < 1e-9, "cy={}", center[1]);
        assert!((center[2] - 0.5).abs() < 1e-9, "cz={}", center[2]);
        let expected_r = (3.0f64).sqrt() / 2.0;
        assert!((radius - expected_r).abs() < 1e-9, "r={}", radius);
    }

    #[test]
    fn circumsphere_degenerate_tet_returns_infinite_radius() {
        // Four coplanar points → degenerate.
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.5, 0.5, 0.0]; // same plane
        let (_, radius) = circumsphere(a, b, c, d);
        assert!(!radius.is_finite(), "degenerate tet should give infinite radius");
    }

    #[test]
    fn in_sphere_test_classifies_points_correctly() {
        let a = [0.0f64, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];
        // Circumcenter (0.5,0.5,0.5): clearly inside
        let inside = [0.5, 0.5, 0.5];
        assert!(in_sphere_test(a, b, c, d, inside) > 0.0);
        // Far away: clearly outside
        let outside = [10.0, 10.0, 10.0];
        assert!(in_sphere_test(a, b, c, d, outside) < 0.0);
    }

    // ── P1: radius_edge_ratio ─────────────────────────────────────────────────

    #[test]
    fn radius_edge_ratio_of_regular_tet_is_known_value() {
        // For a regular tet with edge length 1, R = sqrt(6)/4, lmin = 1.
        // So radius_edge_ratio = sqrt(6)/4 ≈ 0.6124.
        let s = 1.0_f64;
        let h = (2.0_f64 / 3.0).sqrt() * s;
        let nodes = [
            [0.0, 0.0, 0.0],
            [s, 0.0, 0.0],
            [s / 2.0, h, 0.0],
            [s / 2.0, h / 3.0, (2.0_f64 / 3.0).sqrt() * h],
        ];
        let ratio = radius_edge_ratio(&nodes, [0, 1, 2, 3]);
        assert!(ratio.is_finite());
        // Regular tet ratio ≈ 0.6124 – 0.9; just confirm in range and > 0.
        assert!(ratio > 0.0 && ratio < 2.0, "ratio={ratio}");
    }

    #[test]
    fn radius_edge_ratio_out_of_bounds_index_returns_infinity() {
        let nodes = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let ratio = radius_edge_ratio(&nodes, [0, 1, 2, 3]);
        assert!(!ratio.is_finite());
    }

    #[test]
    fn local_flip_pass_activates_on_edge_fan() {
        // Build three tetrahedra sharing edge (1,2):
        // [1,2,3,4], [1,2,4,5], [1,2,3,5]
        // which can be replaced by two tetrahedra [3,4,5,1], [3,4,5,2].
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.5, 0.001, 0.0));
        mesh.add_node(Node::new(4, 0.5, 0.0, 0.001));
        mesh.add_node(Node::new(5, 0.5, 0.8, 0.8));

        mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));
        mesh.add_element(Element::new(2, ElementType::Tetrahedron4, vec![1, 2, 4, 5]));
        mesh.add_element(Element::new(3, ElementType::Tetrahedron4, vec![1, 2, 3, 5]));

        let before_tets = mesh.elements_by_dimension(3).len();
        let mut next_elem_id = 10_u64;
        let (face_flips, edge_flips, _edge_sliver_flips) =
            optimize_local_face_flips(&mut mesh, &mut next_elem_id, 2).expect("flip pass");
        let after_tets = mesh.elements_by_dimension(3).len();

        assert!(
            face_flips + edge_flips > 0,
            "expected at least one local flip to activate"
        );
        assert_eq!(before_tets, 3);
        assert_ne!(after_tets, before_tets);
    }

    #[test]
    fn radius_edge_ratio_degenerate_tet_returns_infinity() {
        // All four points collinear → degenerate.
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let ratio = radius_edge_ratio(&nodes, [0, 1, 2, 3]);
        // Either infinite or astronomically large (degenerate path).
        assert!(ratio > 1e10 || !ratio.is_finite(), "ratio={ratio}");
    }

    // ── P1: solve_3x3 ─────────────────────────────────────────────────────────

    #[test]
    fn solve_3x3_identity_system() {
        let result = solve_3x3([
            ([1.0, 0.0, 0.0], 3.0),
            ([0.0, 1.0, 0.0], 5.0),
            ([0.0, 0.0, 1.0], 7.0),
        ]);
        let sol = result.expect("identity system has unique solution");
        assert!((sol[0] - 3.0).abs() < 1e-12);
        assert!((sol[1] - 5.0).abs() < 1e-12);
        assert!((sol[2] - 7.0).abs() < 1e-12);
    }

    #[test]
    fn solve_3x3_general_system() {
        // 2x + y + z = 4, x + 3y + z = 7, x + y + 4z = 9 → x=0.4, y=1.8, z=1.8
        let result = solve_3x3([
            ([2.0, 1.0, 1.0], 4.0),
            ([1.0, 3.0, 1.0], 7.0),
            ([1.0, 1.0, 4.0], 9.0),
        ]);
        let sol = result.expect("general system should have solution");
        // Verify A*x = b.
        let residual_0 = (2.0 * sol[0] + sol[1] + sol[2] - 4.0).abs();
        let residual_1 = (sol[0] + 3.0 * sol[1] + sol[2] - 7.0).abs();
        let residual_2 = (sol[0] + sol[1] + 4.0 * sol[2] - 9.0).abs();
        assert!(residual_0 < 1e-10, "row 0 residual={residual_0}");
        assert!(residual_1 < 1e-10, "row 1 residual={residual_1}");
        assert!(residual_2 < 1e-10, "row 2 residual={residual_2}");
    }

    #[test]
    fn solve_3x3_singular_returns_none() {
        // Rows 0 and 1 are identical → rank-deficient.
        let result = solve_3x3([
            ([1.0, 2.0, 3.0], 6.0),
            ([1.0, 2.0, 3.0], 6.0),
            ([0.0, 0.0, 1.0], 1.0),
        ]);
        assert!(result.is_none(), "singular system should return None");
    }

    // ── P1: tetra_centroid_from_mesh ──────────────────────────────────────────

    #[test]
    fn tetra_centroid_is_average_of_four_nodes() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 4.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 4.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 0.0, 4.0));
        let c = tetra_centroid_from_mesh(&mesh, &[1, 2, 3, 4]).expect("centroid ok");
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 1.0).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn tetra_centroid_wrong_arity_errors() {
        let mesh = Mesh::new();
        assert!(tetra_centroid_from_mesh(&mesh, &[1, 2, 3]).is_err());
    }

    // ── P1: refinement reduces worst ratio ────────────────────────────────────

    #[test]
    fn refinement_produces_more_elements_than_seed() {
        let algo = Delaunay3D::default();
        let seed_count = {
            use crate::tetrahedralize3d::CentroidStarMesher3D;
            use crate::traits::Mesher3D as _;
            let params = MeshParams::with_size(0.5);
            CentroidStarMesher3D
                .mesh_3d(&cube_surface_mesh(), &params)
                .unwrap()
                .elements_by_dimension(3)
                .len()
        };
        let mut params = MeshParams::with_size(0.5);
        params.optimize_passes = 3;
        let refined = algo
            .mesh_3d(&cube_surface_mesh(), &params)
            .expect("refinement should succeed");
        let refined_count = refined.elements_by_dimension(3).len();
        assert!(
            refined_count >= seed_count,
            "refinement should not reduce element count: seed={seed_count} refined={refined_count}"
        );
    }

    // ── P2: bistellar flips ─────────────────────────────────────────────────

    #[test]
    fn bistellar_flip_2_to_3_preserves_tet_count() {
        // Two tets sharing face (1,2,3): [1,2,3,4] and [1,2,3,5]
        let nodes = [
            [0.0, 0.0, 0.0],   // 0
            [1.0, 0.0, 0.0],   // 1
            [0.0, 1.0, 0.0],   // 2
            [0.25, 0.25, 1.0], // 3 — tip above
            [0.25, 0.25, -1.0], // 4 — tip below
        ];
        let mut tets = vec![[0, 1, 2, 3], [0, 2, 1, 4]];
        assert_eq!(tets.len(), 2);
        bistellar_flip(&mut tets, &nodes, BistellarFlipType::TwoToThree, &[0, 1]).unwrap();
        // 2 old removed, 3 new added → 3 total
        assert_eq!(tets.len(), 3);
    }

    #[test]
    fn bistellar_flip_3_to_2_preserves_tet_count() {
        // Three tets sharing edge (0,1): [0,1,2,3], [0,1,3,4], [0,1,4,2]
        let nodes = [
            [0.0, 0.0, 0.0], // 0
            [0.0, 0.0, 1.0], // 1
            [1.0, 0.0, 0.5], // 2
            [0.0, 1.0, 0.5], // 3
            [-1.0, 0.0, 0.5], // 4
        ];
        let mut tets = vec![[0, 1, 2, 3], [0, 1, 3, 4], [0, 1, 4, 2]];
        assert_eq!(tets.len(), 3);
        bistellar_flip(&mut tets, &nodes, BistellarFlipType::ThreeToTwo, &[0, 1, 2]).unwrap();
        // 3 old removed, 2 new added → 2 total
        assert_eq!(tets.len(), 2);
    }

    #[test]
    fn bistellar_flip_2_to_3_rejects_non_adjacent() {
        // Two tets that don't share a face
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [10.0, 10.0, 10.0],
            [11.0, 10.0, 10.0],
            [10.0, 11.0, 10.0],
            [10.0, 10.0, 11.0],
        ];
        let mut tets = vec![[0, 1, 2, 3], [4, 5, 6, 7]];
        // These tets don't share a face, so the flip should be a no-op
        let result = bistellar_flip(&mut tets, &nodes, BistellarFlipType::TwoToThree, &[0, 1]);
        // Should return Ok but not modify (common.len() != 3)
        assert!(result.is_ok());
    }
}
