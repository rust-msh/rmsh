//! Phase R.B — Face imprinting and gap/overlap detection.
//!
//! **Face imprinting** (`imprint_brep`): splits each face of `target` wherever
//! the boundary of `tool` crosses it, without performing a boolean classification.
//! The result is a new BRep whose faces share edges with the tool boundary — a
//! prerequisite for conformal meshing (FEM/FDTD).
//!
//! **Gap/overlap detection** (`detect_gaps_overlaps`): reports pairs of faces
//! from two BReps that are either too close (gap) or interpenetrating (overlap),
//! using face bounding-box pre-filtering and `closest_point_on_surface`.

use glam::{DVec2, DVec3};
use rcad_kernel::BRep;
use rcad_kernel::geom::*;
use rcad_kernel::projection::closest_point_on_surface;
use rcad_kernel::topology::*;

use crate::bopds::ds::{DS, ShapeOrigin};
use crate::builder::SubFace;
use crate::bvh::{Aabb, Bvh};
use crate::pave_filler::PaveFiller;
use crate::triangulate::triangulate_polygon;

// ─────────────────────────────────────────────────────────────────────────────
// Public types
// ─────────────────────────────────────────────────────────────────────────────

/// Result of imprinting `tool` geometry onto `target`.
#[derive(Debug)]
pub struct ImprintResult {
    /// Modified target BRep whose faces are split wherever the tool boundary crosses.
    pub brep: BRep,
    /// Pairs of (target face index in result, source tool face index) that share
    /// a seam (imprinted edge).
    pub seam_edges: Vec<(usize, usize)>,
}

/// Detected gap between two faces.
#[derive(Debug, Clone)]
pub struct Gap {
    /// Face index in BRep A.
    pub face_a: usize,
    /// Face index in BRep B.
    pub face_b: usize,
    /// Maximum gap distance found between the two faces.
    pub max_gap: f64,
    /// A world-space point on face A that is closest to face B.
    pub sample_point: DVec3,
}

/// Detected overlap (interpenetration) between two faces.
#[derive(Debug, Clone)]
pub struct Overlap {
    /// Face index in BRep A.
    pub face_a: usize,
    /// Face index in BRep B.
    pub face_b: usize,
    /// Estimated penetration depth (positive = overlapping).
    pub penetration_depth: f64,
}

/// Report from gap/overlap detection.
#[derive(Debug, Default)]
pub struct GapOverlapReport {
    pub gaps: Vec<Gap>,
    pub overlaps: Vec<Overlap>,
    /// Pairs of faces (a, b) that are perfectly coincident and coplanar.
    pub shared_faces: Vec<(usize, usize)>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Face imprinting
// ─────────────────────────────────────────────────────────────────────────────

/// Imprint the boundary of `tool` onto the faces of `target`.
///
/// This runs the PaveFiller intersection pass between the two BReps, then splits
/// each target face by the intersection curves recorded in its `FaceInfo`.
/// No boolean classification is performed — all faces of `target` are preserved,
/// but split where the tool boundary crosses them.
///
/// Analogy: OCCT `BRepAlgoAPI_Splitter` (lightweight variant — keeps all target faces).
pub fn imprint_brep(target: &BRep, tool: &BRep) -> ImprintResult {
    // Run PaveFiller to compute intersections
    let mut ds = DS::new(target, tool);
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();

    // Identify which DS faces came from target (ShapeA)
    let target_face_indices: Vec<usize> = ds
        .faces
        .iter()
        .enumerate()
        .filter(|(_, f)| f.origin == ShapeOrigin::ShapeA)
        .map(|(i, _)| i)
        .collect();

    // Identify tool face indices per DS face (for seam tracking)
    let tool_face_indices: Vec<usize> = ds
        .faces
        .iter()
        .enumerate()
        .filter(|(_, f)| f.origin == ShapeOrigin::ShapeB)
        .map(|(i, _)| i)
        .collect();

    let mut result_faces: Vec<Face> = Vec::new();
    let mut seam_edges: Vec<(usize, usize)> = Vec::new();

    for &dfi in &target_face_indices {
        let sub_faces = split_face_by_curves(&ds, dfi);

        let has_intersection = !ds.faces[dfi].face_info.curves_in.is_empty();
        let result_face_start = result_faces.len();

        for sf in sub_faces {
            let triangles = triangulate_polygon(&sf.boundary, sf.normal);
            result_faces.push(Face {
                outer_wire: Wire {
                    edges: sf
                        .boundary
                        .windows(2)
                        .enumerate()
                        .map(|(i, _)| WireEdge {
                            idx: i,
                            forward: true,
                        })
                        .collect(),
                },
                inner_wires: vec![],
                normal: sf.normal,
                triangles,
                    mesh_dirty: false,
            });
        }

        // Record seam edges for every tool face that has curves on this target face
        if has_intersection {
            for &tfi in &tool_face_indices {
                // Check if any curve on this target face came from a FF interference with this tool face
                let shares_curve = ds.interferences.iter().any(|iv| {
                    if let crate::bopds::ds::Interference::FaceFace { f1, f2, curves, .. } = iv {
                        let (ta, tb) = if *f1 == dfi { (*f1, *f2) } else { (*f2, *f1) };
                        ta == dfi && tb == tfi && !curves.is_empty()
                    } else {
                        false
                    }
                });
                if shares_curve {
                    for ri in result_face_start..result_faces.len() {
                        seam_edges.push((ri, ds.faces[tfi].source_face_idx));
                    }
                }
            }
        }
    }

    // Assemble result BRep from split faces
    let brep = BRep {
        vertices: target.vertices.clone(),
        edges: target.edges.clone(),
        solids: vec![Solid {
            shells: vec![Shell {
                faces: result_faces,
            }],
        }],
        geom: target.geom.clone(),
        compound: None,
        compsolid: None,
    };

    ImprintResult { brep, seam_edges }
}

/// Split a single DS face by its intersection curves.
/// Shared with builder logic — produces a list of SubFace.
fn split_face_by_curves(ds: &DS, face_idx: usize) -> Vec<SubFace> {
    let face = &ds.faces[face_idx];
    let fi = &face.face_info;

    if fi.curves_in.is_empty() {
        let boundary = face
            .boundary_verts
            .iter()
            .map(|&vi| ds.vertices[vi].point)
            .collect();
        return vec![SubFace {
            boundary,
            surface: face.surface.clone(),
            normal: face.normal,
            uv_centroid: None,
            sample_override: None,
            uv_domain: None,
            inner_wires: vec![],
            analytic_outer_circle: None,
        }];
    }

    match &face.surface.clone() {
        Surface3::Plane(plane) => split_planar_face_simple(ds, face_idx, plane),
        Surface3::Cylinder(_)
        | Surface3::Sphere(_)
        | Surface3::Cone(_)
        | Surface3::Torus(_) => split_curved_face(ds, face_idx),
        _ => {
            // For other curved surfaces: return whole face for now
            let boundary = face
                .boundary_verts
                .iter()
                .map(|&vi| ds.vertices[vi].point)
                .collect();
            vec![SubFace {
                boundary,
                surface: face.surface.clone(),
                normal: face.normal,
                uv_centroid: None,
                sample_override: None,
            uv_domain: None,
            inner_wires: vec![],
            analytic_outer_circle: None,
            }]
        }
    }
}

fn split_planar_face_simple(ds: &DS, face_idx: usize, plane: &Plane) -> Vec<SubFace> {
    use crate::inttools::edge_face::plane_local_basis;

    let face = &ds.faces[face_idx];
    let boundary_3d: Vec<DVec3> = face
        .boundary_verts
        .iter()
        .map(|&vi| ds.vertices[vi].point)
        .collect();

    let mut segments: Vec<(DVec3, DVec3)> = Vec::new();
    for &ci in &face.face_info.curves_in {
        let ic = &ds.intersection_curves[ci];
        let p0 = ds.vertices[ic.start_vertex].point;
        let p1 = ds.vertices[ic.end_vertex].point;
        if (p1 - p0).length_squared()
            > crate::tolerance::TOLERANCE_ABS * crate::tolerance::TOLERANCE_ABS
        {
            segments.push((p0, p1));
        }
    }

    if segments.is_empty() {
        return vec![SubFace {
            boundary: boundary_3d,
            surface: face.surface.clone(),
            normal: face.normal,
            uv_centroid: None,
            sample_override: None,
            uv_domain: None,
            inner_wires: vec![],
            analytic_outer_circle: None,
        }];
    }

    let (u_axis, v_axis) = plane_local_basis(plane);

    let project = |p: DVec3| -> [f64; 2] {
        let d = p - plane.origin;
        [d.dot(u_axis), d.dot(v_axis)]
    };
    let unproject = |uv: [f64; 2]| -> DVec3 { plane.origin + u_axis * uv[0] + v_axis * uv[1] };

    let mut polygons_2d: Vec<Vec<[f64; 2]>> =
        vec![boundary_3d.iter().map(|&p| project(p)).collect()];

    for (seg_a, seg_b) in &segments {
        let sa = project(*seg_a);
        let sb = project(*seg_b);
        let mut next: Vec<Vec<[f64; 2]>> = Vec::new();
        for poly in polygons_2d.drain(..) {
            let split = split_poly_2d(&poly, sa, sb);
            next.extend(split);
        }
        polygons_2d = next;
    }

    polygons_2d
        .into_iter()
        .filter(|p| p.len() >= 3)
        .map(|poly_2d| {
            let boundary: Vec<DVec3> = poly_2d.iter().map(|&uv| unproject(uv)).collect();
            SubFace {
                boundary,
                surface: face.surface.clone(),
                normal: face.normal,
                uv_centroid: None,
            sample_override: None,
            uv_domain: None,
            inner_wires: vec![],
            analytic_outer_circle: None,
            }
        })
        .collect()
}

/// Split a 2D polygon by a directed segment. Returns 1 polygon if no split, 2 if split.
fn split_poly_2d(poly: &[[f64; 2]], sa: [f64; 2], sb: [f64; 2]) -> Vec<Vec<[f64; 2]>> {
    let n = poly.len();
    let seg_dir = [sb[0] - sa[0], sb[1] - sa[1]];

    let signed_dist =
        |p: [f64; 2]| -> f64 { seg_dir[0] * (p[1] - sa[1]) - seg_dir[1] * (p[0] - sa[0]) };

    let sides: Vec<f64> = poly.iter().map(|&p| signed_dist(p)).collect();

    // Find crossings
    let mut crossings: Vec<(usize, [f64; 2])> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let di = sides[i];
        let dj = sides[j];
        if di * dj < 0.0 {
            let t = di / (di - dj);
            let cx = poly[i][0] + t * (poly[j][0] - poly[i][0]);
            let cy = poly[i][1] + t * (poly[j][1] - poly[i][1]);
            crossings.push((i, [cx, cy]));
        }
    }

    if crossings.len() < 2 {
        return vec![poly.to_vec()];
    }

    // Use first two crossings
    let (i0, c0) = crossings[0];
    let (i1, c1) = crossings[1];

    let (ia, ib, ca, cb) = if i0 <= i1 {
        (i0, i1, c0, c1)
    } else {
        (i1, i0, c1, c0)
    };

    // Sub-poly A: [0..=ia] + ca + cb + [ib+1..]
    let mut sub_a = poly[..=ia].to_vec();
    sub_a.push(ca);
    sub_a.push(cb);
    sub_a.extend_from_slice(&poly[ib + 1..]);

    // Sub-poly B: [ia+1..=ib] + cb + ca
    let mut sub_b = poly[ia + 1..=ib].to_vec();
    sub_b.push(cb);
    sub_b.push(ca);

    let mut result = Vec::new();
    if sub_a.len() >= 3 {
        result.push(sub_a);
    }
    if sub_b.len() >= 3 {
        result.push(sub_b);
    }
    if result.is_empty() {
        result.push(poly.to_vec());
    }
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Gap / overlap detection
// ─────────────────────────────────────────────────────────────────────────────

/// Detect gaps and overlaps between two BReps.
///
/// For each pair of faces (one from each BRep) that are within `tolerance` of
/// each other, samples points on face A and measures the distance to surface B.
///
/// - Distance ∈ (0, tolerance]: **Gap**
/// - Distance ≈ 0 and normals anti-parallel: **SharedFace**
/// - Distance < 0 (interpenetration, estimated): **Overlap**
pub fn detect_gaps_overlaps(a: &BRep, b: &BRep, tolerance: f64) -> GapOverlapReport {
    let mut report = GapOverlapReport::default();

    // Flatten faces with their surface indices
    let faces_a = collect_faces_with_surfaces(a);
    let faces_b = collect_faces_with_surfaces(b);

    for (fa_idx, fa_pts, _fa_surf, fa_normal) in &faces_a {
        // Bounding box of face A
        let (a_min, a_max) = aabb(fa_pts);

        for (fb_idx, _fb_pts, fb_surf, fb_normal) in &faces_b {
            let (b_min, b_max) = {
                let fb_pts2 = collect_face_points(b, *fb_idx);
                aabb(&fb_pts2)
            };

            // AABB pre-filter: skip if clearly too far
            let gap_max = (b_min - a_max).max(a_min - b_max).max_element();
            if gap_max > tolerance * 2.0 + 1.0 {
                continue;
            }

            // Sample up to 5 points on face A and measure distance to surface B
            let samples = sample_face_points(fa_pts, 5);
            let mut max_dist: f64 = f64::NEG_INFINITY;
            let mut min_dist: f64 = f64::INFINITY;
            let mut closest_sample = fa_pts[0];

            for &sp in &samples {
                let proj = closest_point_on_surface(fb_surf, sp, 8);
                let d = proj.distance;
                if d < min_dist {
                    min_dist = d;
                    closest_sample = sp;
                }
                if d > max_dist {
                    max_dist = d;
                }
            }

            // Classify
            let normals_antiparallel = fa_normal.dot(*fb_normal) < -0.9;

            if min_dist.abs() < tolerance * 0.1 && normals_antiparallel {
                report.shared_faces.push((*fa_idx, *fb_idx));
            } else if min_dist > 0.0 && min_dist <= tolerance {
                report.gaps.push(Gap {
                    face_a: *fa_idx,
                    face_b: *fb_idx,
                    max_gap: max_dist,
                    sample_point: closest_sample,
                });
            } else if min_dist < -tolerance * 0.1 {
                report.overlaps.push(Overlap {
                    face_a: *fa_idx,
                    face_b: *fb_idx,
                    penetration_depth: -min_dist,
                });
            }
        }
    }

    report
}

// ─────────────────────────────────────────────────────────────────────────────
// Minimum distance
// ─────────────────────────────────────────────────────────────────────────────

/// 计算两个 BRep 之间的最小距离。
///
/// 使用 BVH 加速面对候选筛选，然后对候选面对采样计算距离。
/// 返回两个 BRep 表面之间的最短距离（0.0 表示相交或接触）。
pub fn min_distance(a: &BRep, b: &BRep) -> f64 {
    use crate::bvh::{Aabb, Bvh};

    // 空形状视为无穷远
    if a.solids.is_empty() || b.solids.is_empty() {
        return f64::INFINITY;
    }
    if a.solids[0].shells.is_empty() || b.solids[0].shells.is_empty() {
        return f64::INFINITY;
    }

    let bvh_b = Bvh::build(b);

    let faces_a = collect_faces_with_surfaces(a);
    let faces_b = collect_faces_with_surfaces(b);

    let mut global_min = f64::INFINITY;

    for (fa_idx, fa_pts, _fa_surf, _fa_normal) in &faces_a {
        if fa_pts.is_empty() {
            continue;
        }
        // 用面 A 的 AABB 在 BVH B 中查询候选面
        let (a_min, a_max) = aabb(fa_pts);
        let query = Aabb { min: a_min, max: a_max };
        let candidate_b_faces = bvh_b.query_aabb(&query);

        // 若没有候选，也检查 B 最近面（用面中心查）
        let centroid_a = fa_pts.iter().copied().sum::<DVec3>() / fa_pts.len() as f64;
        let nearest = bvh_b.nearest_faces(centroid_a, f64::INFINITY, 3);
        let mut check_faces: Vec<usize> = candidate_b_faces.clone();
        for (fi, _) in &nearest {
            if !check_faces.contains(fi) {
                check_faces.push(*fi);
            }
        }

        for fb_idx in check_faces {
            if fb_idx >= faces_b.len() {
                continue;
            }
            let (_, _fb_pts, fb_surf, _) = &faces_b[fb_idx];

            // 对面 A 采样点，计算到面 B 曲面的距离
            let samples = sample_face_points(fa_pts, 8);
            for &sp in &samples {
                let proj = closest_point_on_surface(fb_surf, sp, 8);
                let d = proj.distance.abs();
                if d < global_min {
                    global_min = d;
                }
            }
        }
        let _ = fa_idx;
    }

    global_min
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Collect (flat_face_idx, vertex_points, surface, normal) for every face in brep.
fn collect_faces_with_surfaces(brep: &BRep) -> Vec<(usize, Vec<DVec3>, Surface3, DVec3)> {
    let mut result = Vec::new();
    let mut flat_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let pts: Vec<DVec3> = face
                    .outer_wire
                    .edges
                    .iter()
                    .filter_map(|we| brep.edges.get(we.idx))
                    .map(|e| brep.vertices[e.start].point)
                    .collect();

                let surface = brep
                    .geom
                    .face_surface
                    .get(flat_idx)
                    .and_then(|&si| si)
                    .map(|si| brep.geom.surfaces[si].clone())
                    .unwrap_or_else(|| {
                        Surface3::Plane(Plane {
                            origin: DVec3::ZERO,
                            normal: face.normal,
                        })
                    });

                result.push((flat_idx, pts, surface, face.normal));
                flat_idx += 1;
            }
        }
    }
    result
}

fn collect_face_points(brep: &BRep, flat_idx: usize) -> Vec<DVec3> {
    let mut idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if idx == flat_idx {
                    return face
                        .outer_wire
                        .edges
                        .iter()
                        .filter_map(|we| brep.edges.get(we.idx))
                        .map(|e| brep.vertices[e.start].point)
                        .collect();
                }
                idx += 1;
            }
        }
    }
    vec![]
}

fn aabb(pts: &[DVec3]) -> (DVec3, DVec3) {
    if pts.is_empty() {
        return (DVec3::ZERO, DVec3::ZERO);
    }
    let mut mn = pts[0];
    let mut mx = pts[0];
    for &p in pts.iter().skip(1) {
        mn = mn.min(p);
        mx = mx.max(p);
    }
    (mn, mx)
}

/// Pick up to `n` evenly spaced sample points from a face boundary.
fn sample_face_points(pts: &[DVec3], n: usize) -> Vec<DVec3> {
    if pts.is_empty() {
        return vec![];
    }
    let step = (pts.len() as f64 / n as f64).ceil() as usize;
    let step = step.max(1);
    pts.iter().step_by(step).copied().collect()
}

// Split a curved face (Cylinder, Sphere, Cone, Torus) using parameter-space (UV) 2D clipping.
// Falls back to legacy method when UV data or PCurves are missing.
fn split_curved_face(ds: &DS, face_idx: usize) -> Vec<SubFace> {
    let face = &ds.faces[face_idx];

    // Need UV boundary to operate in parameter space
    let uv_boundary = match &face.uv_boundary {
        Some(b) if b.len() >= 3 => b.clone(),
        _ => return split_curved_face_legacy(ds, face_idx),
    };

    let surface = face.surface.clone();
    let normal = face.normal;
    let is_periodic_u = matches!(&surface, Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Torus(_));
    let is_sphere = matches!(&surface, Surface3::Sphere(_));
    let u_period = if is_periodic_u || is_sphere {
        std::f64::consts::TAU
    } else {
        0.0
    };

    // Collect 2D trim polylines from PCurves for each intersection curve
    let mut trim_polylines: Vec<Vec<DVec2>> = Vec::new();
    for &ci in &face.face_info.curves_in {
        if let Some(pcurve) = find_pcurve_for_face(ds, ci, face_idx) {
            let ic = &ds.intersection_curves[ci];
            let [t0, t1] = ic.t_range;
            const N: usize = 32;
            let raw_pts: Vec<DVec2> = match &pcurve {
                // BSpline PCurves from polyline_pcurve_by_projection use
                // chord-length parameterization normalized to [0,1].
                rcad_kernel::geom::Curve2d::BSpline(_) => (0..=N)
                    .map(|i| {
                        let t = i as f64 / N as f64;
                        pcurve.point_at(t)
                    })
                    .collect(),
                // Analytic curves use the same t parameterization as the 3D intersection curve.
                _ => (0..=N)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * i as f64 / N as f64;
                        pcurve.point_at(t)
                    })
                    .collect(),
            };
            let pts = if u_period > 0.0 {
                let pts = unwrap_u_polyline(raw_pts, u_period);
                if pts.len() >= 2 {
                    let u_span = pts.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max)
                        - pts.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                    if u_span > std::f64::consts::PI {
                        let u_mid = pts.iter().map(|p| p.x).sum::<f64>() / pts.len() as f64;
                        let offset = (u_mid / u_period).floor() * u_period;
                        pts.into_iter()
                            .map(|p| DVec2::new(p.x - offset, p.y))
                            .collect::<Vec<_>>()
                    } else {
                        pts
                    }
                } else {
                    pts
                }
            } else {
                raw_pts
            };
            if pts.len() >= 2 {
                trim_polylines.push(pts);
            }
        }
    }

    // If no PCurves available, fall back to legacy method
    if trim_polylines.is_empty() {
        return split_curved_face_legacy(ds, face_idx);
    }

    // Split UV polygon by each trim polyline
    let mut uv_polygons: Vec<Vec<DVec2>> = vec![uv_boundary];

    for trim in &trim_polylines {
        let mut next: Vec<Vec<DVec2>> = Vec::new();
        for poly in uv_polygons.drain(..) {
            let effective_trim = if u_period > 0.0 {
                periodic_trim_to_open_isoline(&poly, trim, u_period)
                    .unwrap_or_else(|| trim.clone())
            } else {
                trim.clone()
            };
            let halves = split_uv_polygon_by_trim(&poly, &effective_trim);
            next.extend(halves);
        }
        uv_polygons = next;
    }

    // Map each UV sub-polygon back to 3D
    uv_polygons
        .into_iter()
        .filter(|p| p.len() >= 3)
        .map(|uv_poly| {
            let n = uv_poly.len() as f64;
            let centroid_uv = uv_poly.iter().copied().sum::<DVec2>() / n;

            let boundary: Vec<DVec3> = match &surface {
                Surface3::Sphere(_) | Surface3::Cone(_) => {
                    curved_subface_boundary_3d(&uv_poly, &trim_polylines, &surface)
                }
                _ => uv_poly.iter().map(|uv| surface.point_at(uv.x, uv.y)).collect(),
            };
            // For curved surfaces, compute the actual surface normal at the centroid UV
            let sub_normal = {
                let computed = surface.normal_at(centroid_uv.x, centroid_uv.y);
                if computed.length_squared() > 0.5 {
                    computed
                } else {
                    normal
                }
            };
            SubFace {
                boundary,
                surface: surface.clone(),
                normal: sub_normal,
                uv_centroid: Some(centroid_uv),
                sample_override: None,
            uv_domain: None,
            inner_wires: vec![],
            analytic_outer_circle: None,
            }
        })
        .collect()
}

// Legacy approximate method for curved face splitting.
fn split_curved_face_legacy(ds: &DS, face_idx: usize) -> Vec<SubFace> {
    let face = &ds.faces[face_idx];
    let surface = face.surface.clone();
    let normal = face.normal;

    // Collect all intersection polylines for this face
    let mut all_polylines: Vec<Vec<DVec3>> = Vec::new();
    for &ci in &face.face_info.curves_in {
        let ic = &ds.intersection_curves[ci];
        if ic.polyline.len() >= 2 {
            all_polylines.push(ic.polyline.clone());
        } else {
            // Analytic curve - sample it into a polyline
            let pts: Vec<DVec3> = (0..=16)
                .map(|i| {
                    let t = ic.t_range[0] + (ic.t_range[1] - ic.t_range[0]) * i as f64 / 16.0;
                    use rcad_kernel::CurveEval;
                    ic.curve.point_at(t)
                })
                .collect();
            all_polylines.push(pts);
        }
    }

    if all_polylines.is_empty() {
        let boundary = face
            .boundary_verts
            .iter()
            .map(|&vi| ds.vertices[vi].point)
            .collect();
        return vec![SubFace {
            boundary,
            surface,
            normal,
            uv_centroid: None,
            sample_override: None,
            uv_domain: None,
            inner_wires: vec![],
            analytic_outer_circle: None,
        }];
    }

    // Collect boundary vertices
    let boundary_pts: Vec<DVec3> = face
        .boundary_verts
        .iter()
        .map(|&vi| ds.vertices[vi].point)
        .collect();

    if boundary_pts.len() < 3 {
        return vec![SubFace {
            boundary: boundary_pts,
            surface,
            normal,
            uv_centroid: None,
            sample_override: None,
            uv_domain: None,
            inner_wires: vec![],
            analytic_outer_circle: None,
        }];
    }

    // For each intersection polyline, split the boundary
    let mut result_boundaries: Vec<Vec<DVec3>> = vec![boundary_pts];

    for polyline in &all_polylines {
        let (Some(&seg_start), Some(&seg_end)) = (polyline.first(), polyline.last()) else {
            continue;
        };

        let mut next_result: Vec<Vec<DVec3>> = Vec::new();
        for bnd in result_boundaries.drain(..) {
            let n = bnd.len();
            if n < 3 {
                next_result.push(bnd);
                continue;
            }

            let Some((i_start, _)) = bnd
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    a.distance_squared(seg_start)
                        .partial_cmp(&b.distance_squared(seg_start))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
            else {
                next_result.push(bnd);
                continue;
            };
            let Some((i_end, _)) = bnd
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    a.distance_squared(seg_end)
                        .partial_cmp(&b.distance_squared(seg_end))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
            else {
                next_result.push(bnd);
                continue;
            };

            let (ia, ib, p_a, p_b) = if i_start <= i_end {
                (i_start, i_end, seg_start, seg_end)
            } else {
                (i_end, i_start, seg_end, seg_start)
            };

            if ia == ib {
                next_result.push(bnd);
                continue;
            }

            // Sub-face A
            let mut sub_a: Vec<DVec3> = bnd[..=ia].to_vec();
            sub_a.push(p_a);
            for &p in polyline.iter().skip(1).rev().skip(1) {
                sub_a.push(p);
            }
            sub_a.push(p_b);
            sub_a.extend_from_slice(&bnd[ib..]);

            // Sub-face B
            let mut sub_b: Vec<DVec3> = bnd[ia..=ib].to_vec();
            sub_b.push(p_b);
            for &p in polyline.iter().skip(1).rev().skip(1) {
                sub_b.push(p);
            }
            sub_b.push(p_a);

            if sub_a.len() >= 3 {
                next_result.push(sub_a);
            }
            if sub_b.len() >= 3 {
                next_result.push(sub_b);
            }
        }
        result_boundaries = next_result;
    }

    result_boundaries
        .into_iter()
        .filter(|b| b.len() >= 3)
        .map(|boundary| SubFace {
            boundary,
            surface: surface.clone(),
            normal,
            uv_centroid: None,
            sample_override: None,
            uv_domain: None,
            inner_wires: vec![],
            analytic_outer_circle: None,
        })
        .collect()
}

// Find the PCurve for the given intersection curve as it lies on the given face.
fn find_pcurve_for_face(
    ds: &DS,
    curve_idx: usize,
    face_idx: usize,
) -> Option<rcad_kernel::geom::Curve2d> {
    use crate::bopds::ds::Interference;
    for interference in &ds.interferences {
        if let Interference::FaceFace { f1, f2, curves, .. } = interference {
            if curves.contains(&curve_idx) {
                let ic = &ds.intersection_curves[curve_idx];
                if *f1 == face_idx {
                    return ic.pcurve_on_a.clone();
                } else if *f2 == face_idx {
                    return ic.pcurve_on_b.clone();
                }
            }
        }
    }
    None
}

/// Compute a robust 3D boundary for a curved sub-face given its UV polygon
/// and trim polylines. Samples UV edges to avoid degenerate polygons at
/// surface singularities (sphere poles, cone apex).
fn curved_subface_boundary_3d(
    uv_poly: &[DVec2],
    trim_polylines_uv: &[Vec<DVec2>],
    surface: &Surface3,
) -> Vec<DVec3> {
    use crate::tolerance::TOLERANCE_ABS;
    const EDGE_SAMPLES: usize = 8;

    let mut pts: Vec<DVec3> = Vec::new();

    // 1. Sample each UV edge and evaluate 3D positions
    let n = uv_poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = uv_poly[i];
        let b = uv_poly[j];
        for k in 0..EDGE_SAMPLES {
            let t = k as f64 / EDGE_SAMPLES as f64;
            let uv = DVec2::new(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
            pts.push(surface.point_at(uv.x, uv.y));
        }
    }

    // 2. Consecutive deduplication — collapse runs of pole/apex samples
    let mut deduped: Vec<DVec3> = Vec::new();
    for p in &pts {
        if deduped.is_empty() || (*p - deduped[deduped.len() - 1]).length_squared() > TOLERANCE_ABS * TOLERANCE_ABS {
            deduped.push(*p);
        }
    }
    // Close the loop: remove last point if it equals the first
    if deduped.len() > 1 && (deduped[0] - deduped[deduped.len() - 1]).length_squared() < TOLERANCE_ABS * TOLERANCE_ABS {
        deduped.pop();
    }

    // 3. If still degenerate, supplement with trim polyline 3D points
    if deduped.len() < 3 {
        for trim_uv in trim_polylines_uv {
            if trim_uv.len() < 2 {
                continue;
            }
            for uv in trim_uv {
                let p3 = surface.point_at(uv.x, uv.y);
                if point_in_polygon_2d(uv_poly, *uv) || point_near_polygon_2d(uv_poly, *uv, 0.1) {
                    if deduped.iter().all(|q| (p3 - *q).length_squared() > TOLERANCE_ABS * TOLERANCE_ABS) {
                        deduped.push(p3);
                    }
                }
            }
        }
    }

    // 4. Final global dedup
    let mut result: Vec<DVec3> = Vec::new();
    for p in &deduped {
        if result.iter().all(|q| (*p - *q).length_squared() > TOLERANCE_ABS * TOLERANCE_ABS) {
            result.push(*p);
        }
    }

    result
}

// Check if a 2D point is within margin of any edge of a polygon.
fn point_near_polygon_2d(poly: &[DVec2], pt: DVec2, margin: f64) -> bool {
    let n = poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = poly[i];
        let b = poly[j];
        let ab = b - a;
        let len_sq = ab.length_squared();
        let t = if len_sq < 1e-14 { 0.0 } else { ((pt - a).dot(ab) / len_sq).clamp(0.0, 1.0) };
        let closest = a + t * ab;
        if (pt - closest).length() < margin {
            return true;
        }
    }
    false
}

fn unwrap_u_polyline(pts: Vec<DVec2>, period: f64) -> Vec<DVec2> {
    if pts.len() < 2 {
        return pts;
    }
    let mut result = Vec::with_capacity(pts.len());
    result.push(pts[0]);
    let mut offset = 0.0_f64;
    for i in 1..pts.len() {
        let prev_u = result[i - 1].x;
        let curr_u = pts[i].x + offset;
        let diff = curr_u - prev_u;
        if diff > period * 0.5 {
            offset -= period;
        } else if diff < -period * 0.5 {
            offset += period;
        }
        result.push(DVec2::new(pts[i].x + offset, pts[i].y));
    }
    result
}

fn periodic_trim_to_open_isoline(poly: &[DVec2], trim: &[DVec2], u_period: f64) -> Option<Vec<DVec2>> {
    if poly.len() < 3 || trim.len() < 3 || u_period <= 0.0 {
        return None;
    }

    let trim_start = trim[0];
    let trim_end = trim[trim.len() - 1];
    let is_closed = (trim_start - trim_end).length_squared() < 1e-6;
    if !is_closed {
        return None;
    }

    let u_min_trim = trim.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let u_max_trim = trim.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let v_min_trim = trim.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let v_max_trim = trim.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    let u_span = u_max_trim - u_min_trim;
    let v_span = v_max_trim - v_min_trim;

    let poly_u_min = poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let poly_u_max = poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let poly_v_min = poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let poly_v_max = poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    let poly_v_span = poly_v_max - poly_v_min;

    if u_span < 0.9 * u_period {
        return None;
    }
    if poly_v_span <= 1e-12 || v_span > 0.1 * poly_v_span {
        return None;
    }

    let v_level = trim.iter().map(|p| p.y).sum::<f64>() / trim.len() as f64;
    if v_level <= poly_v_min + 1e-9 || v_level >= poly_v_max - 1e-9 {
        return None;
    }

    Some(vec![
        DVec2::new(poly_u_min, v_level),
        DVec2::new(poly_u_max, v_level),
    ])
}

// Split a 2D UV polygon by a 2D trim polyline.
fn split_uv_polygon_by_trim(poly: &[DVec2], trim: &[DVec2]) -> Vec<Vec<DVec2>> {
    let n = poly.len();
    if n < 3 || trim.len() < 2 {
        return vec![poly.to_vec()];
    }

    let trim_start = trim[0];
    let trim_end = trim[trim.len() - 1];

    // Detect closed trim
    let is_closed_trim = (trim_start - trim_end).length_squared() < 1e-6;
    if is_closed_trim {
        let trim_centroid = trim.iter().copied().sum::<DVec2>() / trim.len() as f64;
        let is_inside = point_in_polygon_2d(poly, trim_centroid);
        if is_inside {
            let trim_dedup: Vec<DVec2> = {
                let mut v = trim.to_vec();
                if v.len() > 1 && (v[0] - v[v.len() - 1]).length_squared() < 1e-12 {
                    v.pop();
                }
                v
            };
            if trim_dedup.len() >= 3 {
                return vec![trim_dedup, poly.to_vec()];
            }
        }
        return vec![poly.to_vec()];
    }

    // Find closest point on polygon boundary for trim endpoints
    let closest_on_boundary = |q: DVec2| -> (usize, f64, DVec2) {
        let mut best_edge = 0usize;
        let mut best_t = 0.0f64;
        let mut best_pt = poly[0];
        let mut best_dist = f64::INFINITY;

        for i in 0..n {
            let j = (i + 1) % n;
            let a = poly[i];
            let b = poly[j];
            let ab = b - a;
            let len_sq = ab.dot(ab);
            let t = if len_sq < 1e-14 {
                0.0
            } else {
                ((q - a).dot(ab) / len_sq).clamp(0.0, 1.0)
            };
            let proj = a + t * ab;
            let dist = (q - proj).length_squared();
            if dist < best_dist {
                best_dist = dist;
                best_edge = i;
                best_t = t;
                best_pt = proj;
            }
        }
        (best_edge, best_t, best_pt)
    };

    let (edge_s, _t_s, pt_s) = closest_on_boundary(trim_start);
    let (edge_e, _t_e, pt_e) = closest_on_boundary(trim_end);

    let (ia, ib, p_a, p_b, trim_forward) = if edge_s <= edge_e {
        (edge_s, edge_e, pt_s, pt_e, true)
    } else {
        (edge_e, edge_s, pt_e, pt_s, false)
    };

    if ia == ib {
        return vec![poly.to_vec()];
    }

    let trim_pts: Vec<DVec2> = if trim_forward {
        trim.to_vec()
    } else {
        trim.iter().copied().rev().collect()
    };

    // Sub-polygon A
    let mut sub_a: Vec<DVec2> = poly[..=ia].to_vec();
    sub_a.push(p_a);
    for &p in trim_pts.iter().skip(1).rev().skip(1) {
        sub_a.push(p);
    }
    sub_a.push(p_b);
    sub_a.extend_from_slice(&poly[ib + 1..]);

    // Sub-polygon B
    let mut sub_b: Vec<DVec2> = vec![p_a];
    sub_b.extend_from_slice(&poly[ia + 1..=ib]);
    sub_b.push(p_b);
    for &p in trim_pts.iter().skip(1).rev().skip(1).rev() {
        sub_b.push(p);
    }

    let dedup_2d = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut result: Vec<DVec2> = Vec::new();
        for p in v {
            if result.is_empty() || (p - result[result.len() - 1]).length_squared() > 1e-18 {
                result.push(p);
            }
        }
        if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < 1e-18 {
            result.pop();
        }
        result
    };

    let sub_a = dedup_2d(sub_a);
    let sub_b = dedup_2d(sub_b);

    let mut out = Vec::new();
    if sub_a.len() >= 3 {
        out.push(sub_a);
    }
    if sub_b.len() >= 3 {
        out.push(sub_b);
    }

    if out.is_empty() {
        vec![poly.to_vec()]
    } else {
        out
    }
}

// Check if a 2D point is inside a 2D polygon using ray casting.
fn point_in_polygon_2d(poly: &[DVec2], pt: DVec2) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let vi = poly[i];
        let vj = poly[j];
        if ((vi.y > pt.y) != (vj.y > pt.y))
            && (pt.x < (vj.x - vi.x) * (pt.y - vi.y) / (vj.y - vi.y) + vi.x)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_modeling::builder::{box_brep, sphere_brep};
    use crate::geom_populate::populate_box_geom;

    fn make_box(origin: DVec3, x: DVec3, y: DVec3, w: f64, h: f64, d: f64) -> BRep {
        box_brep(origin, x, y, w, h, d).expect("box creation should succeed")
    }

    fn make_sphere(center: DVec3, radius: f64) -> BRep {
        sphere_brep(center, radius).expect("sphere creation should succeed")
    }

    #[test]
    fn test_imprint_box_onto_box() {
        // Two overlapping boxes
        let mut target = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0);
        let mut tool = make_box(DVec3::new(1.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0);
        
        populate_box_geom(&mut target);
        populate_box_geom(&mut tool);
        
        let result = imprint_brep(&target, &tool);
        
        // Result should have faces split where tool boundary crosses target
        assert!(!result.brep.solids[0].shells[0].faces.is_empty());
        // Should have seam edges where target and tool faces meet
        assert!(!result.seam_edges.is_empty());
    }

    #[test]
    fn test_imprint_no_intersection() {
        // Two non-overlapping boxes
        let mut target = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0);
        let mut tool = make_box(DVec3::new(5.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0);
        
        populate_box_geom(&mut target);
        populate_box_geom(&mut tool);
        
        let result = imprint_brep(&target, &tool);
        
        // No intersection means no seam edges
        assert!(result.seam_edges.is_empty());
        // Target faces should remain unchanged (6 faces)
        assert_eq!(result.brep.solids[0].shells[0].faces.len(), 6);
    }

    #[test]
    fn test_imprint_sphere_onto_box() {
        // Sphere intersecting a box - tests curved surface handling
        let mut target = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0);
        let tool = make_sphere(DVec3::new(1.0, 1.0, 1.0), 1.5);
        
        populate_box_geom(&mut target);
        // Sphere already has geometry populated
        
        let result = imprint_brep(&target, &tool);
        
        // Should have faces (even if just original faces)
        assert!(!result.brep.solids[0].shells[0].faces.is_empty());
    }

    #[test]
    fn test_gap_detection_overlapping_boxes() {
        let a = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0);
        let b = make_box(DVec3::new(2.1, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0);
        
        let report = detect_gaps_overlaps(&a, &b, 0.5);
        
        // There should be a small gap between the boxes
        assert!(!report.gaps.is_empty() || !report.overlaps.is_empty() || !report.shared_faces.is_empty());
    }

    #[test]
    fn test_gap_detection_touching_boxes() {
        let a = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0);
        let b = make_box(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0);
        
        let report = detect_gaps_overlaps(&a, &b, 0.01);
        
        // Boxes are touching, might have shared faces
        // Just verify it doesn't panic
        let _ = report.gaps.len() + report.overlaps.len() + report.shared_faces.len();
    }

    #[test]
    fn test_split_face_by_curves_empty() {
        // Test that split_face_by_curves handles empty curves correctly
        // This is tested indirectly through imprint_brep with non-overlapping shapes
        let target = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0);
        let tool = make_box(DVec3::new(10.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0);
        
        let result = imprint_brep(&target, &tool);
        
        // No intersection means faces should not be split
        assert_eq!(result.brep.solids[0].shells[0].faces.len(), 6);
    }

    #[test]
    fn test_detect_gaps_overlaps_empty_brep() {
        // Empty BRep should not panic
        let empty = BRep {
            vertices: vec![],
            edges: vec![],
            solids: vec![],
            geom: Default::default(),
            compound: None,
            compsolid: None,
        };
        let box_brep = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0);
        
        let report = detect_gaps_overlaps(&empty, &box_brep, 0.1);
        
        // Should return empty report
        assert!(report.gaps.is_empty());
        assert!(report.overlaps.is_empty());
        assert!(report.shared_faces.is_empty());
    }
}
