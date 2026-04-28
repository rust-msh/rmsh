//! BRepExtrema-style distance/extrema calculations.
//!
//! Analogous to OCCT `BRepExtrema` package:
//! - `DistShapeShape`: distance computation between shapes
//! - `Extrema`: find extremum points (closest/furthest)
//! - `ClosestPoint`: closest point queries on geometry
//! - `SupportShapes`: find supporting geometry for a point
//!
//! Uses Newton iteration and sampling-based approaches for robust convergence.
//! Derivatives are computed via finite differences for generality.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, CurveEval, Line3, Surface3, SurfaceEval};
use rcad_kernel::BRep;

use crate::bvh::{Aabb, Bvh};
use crate::tolerance::TOLERANCE_ABS;

/// Finite difference step size for derivative computation.
const H: f64 = 1e-6;

// =============================================================================
// DistShapeShape - Distance between shapes
// =============================================================================

/// Compute the Euclidean distance between two points.
pub fn distance_point_point(p1: DVec3, p2: DVec3) -> f64 {
    (p2 - p1).length()
}

/// Compute the distance from a point to a curve.
///
/// Returns the distance and the parameter value on the curve where the closest point lies.
/// Uses Newton iteration for accurate parameter refinement.
pub fn distance_point_curve(point: DVec3, curve: &Curve3) -> (f64, f64) {
    let (param, closest_pt) = closest_point_on_curve(curve, point);
    let distance = distance_point_point(point, closest_pt);
    (distance, param)
}

/// Compute the distance from a point to a surface.
///
/// Returns the distance and the UV parameters where the closest point lies.
/// Uses Newton iteration for accurate parameter refinement.
pub fn distance_point_surface(point: DVec3, surface: &Surface3) -> (f64, f64, f64) {
    let (uv, closest_pt) = closest_point_on_surface(surface, point);
    let distance = distance_point_point(point, closest_pt);
    (distance, uv.x, uv.y)
}

/// Compute the minimum distance between two curves.
///
/// Returns the distance and the two parameter values on each curve.
/// Uses sampling to find initial candidates, then Newton refinement.
pub fn distance_curve_curve(curve1: &Curve3, curve2: &Curve3) -> (f64, f64, f64) {
    let domain1 = curve_domain(curve1);
    let domain2 = curve_domain(curve2);

    // Sample both curves to find initial candidates
    let n_samples = 32;
    let mut best_dist = f64::INFINITY;
    let mut best_t1 = domain1[0];
    let mut best_t2 = domain2[0];

    for i in 0..=n_samples {
        let t1 = domain1[0] + (domain1[1] - domain1[0]) * i as f64 / n_samples as f64;
        let p1 = curve1.point_at(t1);

        for j in 0..=n_samples {
            let t2 = domain2[0] + (domain2[1] - domain2[0]) * j as f64 / n_samples as f64;
            let p2 = curve2.point_at(t2);
            let dist = (p2 - p1).length();

            if dist < best_dist {
                best_dist = dist;
                best_t1 = t1;
                best_t2 = t2;
            }
        }
    }

    // Newton refinement
    let (refined_t1, refined_t2) = refine_curve_curve_distance(curve1, curve2, domain1, domain2, best_t1, best_t2);
    let p1 = curve1.point_at(refined_t1);
    let p2 = curve2.point_at(refined_t2);
    let final_dist = (p2 - p1).length();

    (final_dist, refined_t1, refined_t2)
}

/// Compute the minimum distance between a curve and a surface.
///
/// Returns the distance and the parameter on the curve plus UV on the surface.
/// Uses sampling to find initial candidates, then Newton refinement.
pub fn distance_curve_surface(curve: &Curve3, surface: &Surface3) -> (f64, f64, f64, f64) {
    let curve_domain = curve_domain(curve);
    let surf_domain = surface_domain(surface);

    // Sample curve and surface to find initial candidates
    let n_curve = 24;
    let n_surf = 12;
    let mut best_dist = f64::INFINITY;
    let mut best_t = curve_domain[0];
    let mut best_u = surf_domain[0];
    let mut best_v = surf_domain[2];

    for i in 0..=n_curve {
        let t = curve_domain[0] + (curve_domain[1] - curve_domain[0]) * i as f64 / n_curve as f64;
        let p_curve = curve.point_at(t);

        for j in 0..=n_surf {
            let u = surf_domain[0] + (surf_domain[1] - surf_domain[0]) * j as f64 / n_surf as f64;
            for k in 0..=n_surf {
                let v = surf_domain[2] + (surf_domain[3] - surf_domain[2]) * k as f64 / n_surf as f64;
                let p_surf = surface.point_at(u, v);
                let dist = (p_surf - p_curve).length();

                if dist < best_dist {
                    best_dist = dist;
                    best_t = t;
                    best_u = u;
                    best_v = v;
                }
            }
        }
    }

    // Newton refinement
    let (refined_t, refined_u, refined_v) = refine_curve_surface_distance(
        curve, surface, curve_domain, surf_domain, best_t, best_u, best_v,
    );
    let p_curve = curve.point_at(refined_t);
    let p_surf = surface.point_at(refined_u, refined_v);
    let final_dist = (p_surf - p_curve).length();

    (final_dist, refined_t, refined_u, refined_v)
}

/// Compute the minimum distance between two surfaces.
///
/// Returns the distance and UV parameters on both surfaces.
/// Uses sampling to find initial candidates, then Newton refinement.
pub fn distance_surface_surface(surf1: &Surface3, surf2: &Surface3) -> (f64, f64, f64, f64, f64) {
    let domain1 = surface_domain(surf1);
    let domain2 = surface_domain(surf2);

    // Sample both surfaces to find initial candidates
    let n_samples = 10;
    let mut best_dist = f64::INFINITY;
    let mut best_u1 = domain1[0];
    let mut best_v1 = domain1[2];
    let mut best_u2 = domain2[0];
    let mut best_v2 = domain2[2];

    for i1 in 0..=n_samples {
        let u1 = domain1[0] + (domain1[1] - domain1[0]) * i1 as f64 / n_samples as f64;
        for j1 in 0..=n_samples {
            let v1 = domain1[2] + (domain1[3] - domain1[2]) * j1 as f64 / n_samples as f64;
            let p1 = surf1.point_at(u1, v1);

            for i2 in 0..=n_samples {
                let u2 = domain2[0] + (domain2[1] - domain2[0]) * i2 as f64 / n_samples as f64;
                for j2 in 0..=n_samples {
                    let v2 = domain2[2] + (domain2[3] - domain2[2]) * j2 as f64 / n_samples as f64;
                    let p2 = surf2.point_at(u2, v2);
                    let dist = (p2 - p1).length();

                    if dist < best_dist {
                        best_dist = dist;
                        best_u1 = u1;
                        best_v1 = v1;
                        best_u2 = u2;
                        best_v2 = v2;
                    }
                }
            }
        }
    }

    // Newton refinement
    let (refined_u1, refined_v1, refined_u2, refined_v2) = refine_surface_surface_distance(
        surf1, surf2, domain1, domain2, best_u1, best_v1, best_u2, best_v2,
    );

    let p1 = surf1.point_at(refined_u1, refined_v1);
    let p2 = surf2.point_at(refined_u2, refined_v2);
    let final_dist = (p2 - p1).length();

    (final_dist, refined_u1, refined_v1, refined_u2, refined_v2)
}

/// Compute the minimum distance between two BRep shapes.
///
/// Returns the distance and the two closest points.
/// Uses BVH acceleration for efficiency.
pub fn distance_brep_brep(brep1: &BRep, brep2: &BRep) -> (f64, DVec3, DVec3) {
    let bvh1 = Bvh::build(brep1);
    let bvh2 = Bvh::build(brep2);

    // Get candidate face pairs
    let candidate_pairs = Bvh::candidate_pairs(&bvh1, &bvh2);

    if candidate_pairs.is_empty() {
        // Fallback: compute bounding box centers distance
        let bb1 = compute_brep_aabb(brep1);
        let bb2 = compute_brep_aabb(brep2);
        let center1 = bb1.center();
        let center2 = bb2.center();
        return ((center2 - center1).length(), center1, center2);
    }

    let mut best_dist = f64::INFINITY;
    let mut best_pt1 = DVec3::ZERO;
    let mut best_pt2 = DVec3::ZERO;

    // Check each candidate pair
    for (fi1, fi2) in candidate_pairs {
        let surf1 = get_brep_surface(brep1, fi1);
        let surf2 = get_brep_surface(brep2, fi2);

        if let (Some(s1), Some(s2)) = (surf1, surf2) {
            let (dist, u1, v1, u2, v2) = distance_surface_surface(&s1, &s2);
            if dist < best_dist {
                best_dist = dist;
                best_pt1 = s1.point_at(u1, v1);
                best_pt2 = s2.point_at(u2, v2);
            }
        }
    }

    // Also check vertex-to-vertex distances
    for v1 in &brep1.vertices {
        for v2 in &brep2.vertices {
            let dist = (v2.point - v1.point).length();
            if dist < best_dist {
                best_dist = dist;
                best_pt1 = v1.point;
                best_pt2 = v2.point;
            }
        }
    }

    (best_dist, best_pt1, best_pt2)
}

// =============================================================================
// Extrema - Find extremum points
// =============================================================================

/// Find the n closest points on a curve to a given point.
///
/// Returns a vector of (parameter, distance) pairs sorted by distance.
pub fn find_closest_points(curve: &Curve3, point: DVec3, n_points: usize) -> Vec<(f64, f64)> {
    let domain = curve_domain(curve);

    // Sample the curve to find local minima
    let n_samples = 100;
    let mut candidates: Vec<(f64, f64)> = Vec::new();

    for i in 0..=n_samples {
        let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
        let p = curve.point_at(t);
        let dist = (p - point).length();
        candidates.push((t, dist));
    }

    // Find local minima
    let mut local_minima: Vec<(f64, f64)> = Vec::new();
    for i in 1..candidates.len() - 1 {
        if candidates[i].1 < candidates[i - 1].1 && candidates[i].1 < candidates[i + 1].1 {
            // Refine using Newton
            let refined_t = refine_point_curve_distance(curve, domain, point, candidates[i].0);
            let refined_dist = (curve.point_at(refined_t) - point).length();
            local_minima.push((refined_t, refined_dist));
        }
    }

    // Also include endpoints
    let (t0, d0) = candidates[0];
    let (tn, dn) = candidates[candidates.len() - 1];
    local_minima.push((t0, d0));
    local_minima.push((tn, dn));

    // Sort by distance and take top n
    local_minima.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    local_minima.truncate(n_points);

    local_minima
}

/// Find the furthest points on a BRep in a given direction.
///
/// Returns the two points that are furthest apart when projected onto the direction.
pub fn find_furthest_points(brep: &BRep, direction: DVec3) -> (DVec3, DVec3) {
    let dir = direction.normalize();
    let mut min_proj = f64::INFINITY;
    let mut max_proj = f64::NEG_INFINITY;
    let mut min_point = DVec3::ZERO;
    let mut max_point = DVec3::ZERO;

    // Check all vertices
    for v in &brep.vertices {
        let proj = v.point.dot(dir);
        if proj < min_proj {
            min_proj = proj;
            min_point = v.point;
        }
        if proj > max_proj {
            max_proj = proj;
            max_point = v.point;
        }
    }

    // Also sample face interiors for more accuracy
    let face_indices = get_all_face_indices(brep);
    for face_idx in face_indices {
        if let Some(surf) = get_brep_surface(brep, face_idx) {
            let domain = surface_domain(&surf);
            let n_samples = 5;
            for i in 0..=n_samples {
                for j in 0..=n_samples {
                    let u = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
                    let v = domain[2] + (domain[3] - domain[2]) * j as f64 / n_samples as f64;
                    let p = surf.point_at(u, v);
                    let proj = p.dot(dir);
                    if proj < min_proj {
                        min_proj = proj;
                        min_point = p;
                    }
                    if proj > max_proj {
                        max_proj = proj;
                        max_point = p;
                    }
                }
            }
        }
    }

    (min_point, max_point)
}

// =============================================================================
// ClosestPoint - Closest point queries
// =============================================================================

/// Find the closest point on a curve to a given point.
///
/// Returns the parameter value and the closest point on the curve.
/// Uses Newton iteration for accuracy.
pub fn closest_point_on_curve(curve: &Curve3, point: DVec3) -> (f64, DVec3) {
    let domain = curve_domain(curve);

    // Initial guess by sampling
    let n_samples = 50;
    let mut best_t = domain[0];
    let mut best_dist = f64::INFINITY;

    for i in 0..=n_samples {
        let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
        let p = curve.point_at(t);
        let dist = (p - point).length();
        if dist < best_dist {
            best_dist = dist;
            best_t = t;
        }
    }

    // Newton refinement
    let refined_t = refine_point_curve_distance(curve, domain, point, best_t);
    let closest = curve.point_at(refined_t);

    (refined_t, closest)
}

/// Find the closest point on a surface to a given point.
///
/// Returns the UV parameters and the closest point on the surface.
/// Uses Newton iteration for accuracy.
pub fn closest_point_on_surface(surface: &Surface3, point: DVec3) -> (DVec2, DVec3) {
    let domain = surface_domain(surface);

    // Initial guess by sampling
    let n_samples = 20;
    let mut best_u = domain[0];
    let mut best_v = domain[2];
    let mut best_dist = f64::INFINITY;

    for i in 0..=n_samples {
        let u = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
        for j in 0..=n_samples {
            let v = domain[2] + (domain[3] - domain[2]) * j as f64 / n_samples as f64;
            let p = surface.point_at(u, v);
            let dist = (p - point).length();
            if dist < best_dist {
                best_dist = dist;
                best_u = u;
                best_v = v;
            }
        }
    }

    // Newton refinement
    let (refined_u, refined_v) = refine_point_surface_distance(surface, domain, point, best_u, best_v);
    let closest = surface.point_at(refined_u, refined_v);

    (DVec2::new(refined_u, refined_v), closest)
}

// =============================================================================
// SupportShapes - Find supporting geometry
// =============================================================================

/// Find the face that supports a point (closest face).
///
/// Returns the face index if found.
pub fn find_supporting_face(brep: &BRep, point: DVec3) -> Option<usize> {
    let mut best_face = None;
    let mut best_dist = f64::INFINITY;
    let tolerance = TOLERANCE_ABS * 100.0;

    let face_indices = get_all_face_indices(brep);
    for face_idx in face_indices {
        if let Some(surf) = get_brep_surface(brep, face_idx) {
            let (uv, closest) = closest_point_on_surface(&surf, point);
            let dist = (closest - point).length();

            // Check if point is within the face boundary (UV domain)
            let domain = surface_domain(&surf);
            if uv.x >= domain[0] - tolerance && uv.x <= domain[1] + tolerance
                && uv.y >= domain[2] - tolerance && uv.y <= domain[3] + tolerance
            {
                if dist < best_dist {
                    best_dist = dist;
                    best_face = Some(face_idx);
                }
            }
        }
    }

    best_face
}

/// Find the edge that supports a point (closest edge).
///
/// Returns the edge index if found.
pub fn find_supporting_edge(brep: &BRep, point: DVec3) -> Option<usize> {
    let mut best_edge = None;
    let mut best_dist = f64::INFINITY;
    let tolerance = TOLERANCE_ABS * 100.0;

    for (edge_idx, edge) in brep.edges.iter().enumerate() {
        // Try to get the curve from geometry store, or create a line from vertices
        let curve = if let Some(c) = get_brep_curve(brep, edge_idx) {
            c
        } else {
            // Create an implicit line from edge vertices
            let start_pt = brep.vertices.get(edge.start).map(|v| v.point)?;
            let end_pt = brep.vertices.get(edge.end).map(|v| v.point)?;
            let dir = (end_pt - start_pt).normalize();
            Curve3::Line(Line3 {
                origin: start_pt,
                direction: dir,
            })
        };

        let (t, closest) = closest_point_on_curve(&curve, point);
        let dist = (closest - point).length();

        // Check if parameter is within edge range
        let edge_range = brep.geom.edge_curve_range.get(edge_idx)
            .and_then(|r| *r)
            .unwrap_or_else(|| curve.default_domain());

        if t >= edge_range[0] - tolerance && t <= edge_range[1] + tolerance {
            if dist < best_dist {
                best_dist = dist;
                best_edge = Some(edge_idx);
            }
        }
    }

    best_edge
}

// =============================================================================
// Internal helper functions
// =============================================================================

/// Get the domain for a curve, handling infinite lines specially.
fn curve_domain(curve: &Curve3) -> [f64; 2] {
    match curve {
        Curve3::Line(_) => [-1e6, 1e6], // Clamp infinite lines to large range
        other => other.default_domain(),
    }
}

/// Get the domain for a surface, handling infinite domains (planes) specially.
fn surface_domain(surface: &Surface3) -> [f64; 4] {
    let domain = surface.default_domain();
    let u0 = if domain[0].is_infinite() { -10.0 } else { domain[0] };
    let u1 = if domain[1].is_infinite() { 10.0 } else { domain[1] };
    let v0 = if domain[2].is_infinite() { -10.0 } else { domain[2] };
    let v1 = if domain[3].is_infinite() { 10.0 } else { domain[3] };
    [u0, u1, v0, v1]
}

/// Compute the AABB of a BRep.
fn compute_brep_aabb(brep: &BRep) -> Aabb {
    let mut aabb = Aabb::empty();
    for v in &brep.vertices {
        aabb.expand_point(v.point);
    }
    aabb
}

/// Get a surface from a BRep by face index.
fn get_brep_surface(brep: &BRep, face_idx: usize) -> Option<Surface3> {
    brep.geom.face_surface.get(face_idx)
        .and_then(|s| *s)
        .and_then(|idx| brep.geom.surfaces.get(idx).cloned())
}

/// Get all face indices from a BRep.
fn get_all_face_indices(brep: &BRep) -> Vec<usize> {
    let mut count = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            count += shell.faces.len();
        }
    }
    (0..count).collect()
}

/// Get a curve from a BRep by edge index.
fn get_brep_curve(brep: &BRep, edge_idx: usize) -> Option<Curve3> {
    brep.geom.edge_curve.get(edge_idx)
        .and_then(|c| *c)
        .and_then(|idx| brep.geom.curves.get(idx).cloned())
}

/// Compute curve derivative via finite differences.
fn curve_derivative(curve: &Curve3, t: f64) -> DVec3 {
    (curve.point_at(t + H) - curve.point_at(t - H)) / (2.0 * H)
}

/// Compute curve second derivative via finite differences.
fn curve_second_derivative(curve: &Curve3, t: f64) -> DVec3 {
    let d_plus = curve_derivative(curve, t + H);
    let d_minus = curve_derivative(curve, t - H);
    (d_plus - d_minus) / (2.0 * H)
}

/// Compute surface partial derivatives via finite differences.
fn surface_derivatives(surface: &Surface3, u: f64, v: f64) -> (DVec3, DVec3) {
    let du = (surface.point_at(u + H, v) - surface.point_at(u - H, v)) / (2.0 * H);
    let dv = (surface.point_at(u, v + H) - surface.point_at(u, v - H)) / (2.0 * H);
    (du, dv)
}

/// Newton refinement for point-to-curve distance.
fn refine_point_curve_distance(curve: &Curve3, domain: [f64; 2], point: DVec3, initial_t: f64) -> f64 {
    let mut t = initial_t;

    const MAX_ITER: usize = 20;
    const TOL: f64 = 1e-10;

    for _ in 0..MAX_ITER {
        let p = curve.point_at(t);
        let d = curve_derivative(curve, t);

        let diff = p - point;
        let f = diff.dot(d);

        let d2 = curve_second_derivative(curve, t);
        let df = d.dot(d) + diff.dot(d2);

        if df.abs() < TOL {
            break;
        }

        let delta = -f / df;
        t += delta;

        // Clamp to domain
        t = t.clamp(domain[0], domain[1]);

        if delta.abs() < TOL {
            break;
        }
    }

    t
}

/// Newton refinement for point-to-surface distance.
fn refine_point_surface_distance(surface: &Surface3, domain: [f64; 4], point: DVec3, initial_u: f64, initial_v: f64) -> (f64, f64) {
    let mut u = initial_u;
    let mut v = initial_v;

    const MAX_ITER: usize = 20;
    const TOL: f64 = 1e-10;

    for _ in 0..MAX_ITER {
        let p = surface.point_at(u, v);
        let (du, dv) = surface_derivatives(surface, u, v);

        let diff = p - point;

        // Gradient of distance squared
        let fu = diff.dot(du);
        let fv = diff.dot(dv);

        // Hessian approximation using finite differences for second derivatives
        let (du_du, du_dv) = surface_derivatives(surface, u + H, v);
        let (_dv_du, dv_dv) = surface_derivatives(surface, u, v + H);

        let d2uu = (du_du - du) / H;
        let d2vv = (dv_dv - dv) / H;
        let d2uv = (du_dv - du) / H;

        let fuu = du.dot(du) + diff.dot(d2uu);
        let fvv = dv.dot(dv) + diff.dot(d2vv);
        let fuv = du.dot(dv) + diff.dot(d2uv);

        // Solve 2x2 system
        let det = fuu * fvv - fuv * fuv;
        if det.abs() < TOL {
            break;
        }

        let du_param = (-fu * fvv + fv * fuv) / det;
        let dv_param = (-fv * fuu + fu * fuv) / det;

        u += du_param;
        v += dv_param;

        // Clamp to domain
        u = u.clamp(domain[0], domain[1]);
        v = v.clamp(domain[2], domain[3]);

        if du_param.abs() < TOL && dv_param.abs() < TOL {
            break;
        }
    }

    (u, v)
}

/// Newton refinement for curve-to-curve distance.
fn refine_curve_curve_distance(
    curve1: &Curve3, curve2: &Curve3,
    domain1: [f64; 2], domain2: [f64; 2],
    t1: f64, t2: f64,
) -> (f64, f64) {
    let mut t1 = t1;
    let mut t2 = t2;

    const MAX_ITER: usize = 30;
    const TOL: f64 = 1e-10;

    for _ in 0..MAX_ITER {
        let p1 = curve1.point_at(t1);
        let p2 = curve2.point_at(t2);

        let d1 = curve_derivative(curve1, t1);
        let d2 = curve_derivative(curve2, t2);

        let diff = p1 - p2;

        // Gradient
        let f1 = diff.dot(d1);
        let f2 = -diff.dot(d2);

        // Hessian
        let d1_2 = curve_second_derivative(curve1, t1);
        let d2_2 = curve_second_derivative(curve2, t2);

        let h11 = d1.dot(d1) + diff.dot(d1_2);
        let h22 = d2.dot(d2) - diff.dot(d2_2);
        let h12 = -d1.dot(d2);

        let det = h11 * h22 - h12 * h12;
        if det.abs() < TOL {
            break;
        }

        let dt1 = (-f1 * h22 + f2 * h12) / det;
        let dt2 = (-f2 * h11 + f1 * h12) / det;

        t1 += dt1;
        t2 += dt2;

        t1 = t1.clamp(domain1[0], domain1[1]);
        t2 = t2.clamp(domain2[0], domain2[1]);

        if dt1.abs() < TOL && dt2.abs() < TOL {
            break;
        }
    }

    (t1, t2)
}

/// Newton refinement for curve-to-surface distance.
fn refine_curve_surface_distance(
    curve: &Curve3, surface: &Surface3,
    curve_domain: [f64; 2], surf_domain: [f64; 4],
    t: f64, u: f64, v: f64,
) -> (f64, f64, f64) {
    let mut t = t;
    let mut u = u;
    let mut v = v;

    const MAX_ITER: usize = 30;
    const TOL: f64 = 1e-10;

    for _ in 0..MAX_ITER {
        let pc = curve.point_at(t);
        let ps = surface.point_at(u, v);

        let dc = curve_derivative(curve, t);
        let (ds_u, ds_v) = surface_derivatives(surface, u, v);

        let diff = pc - ps;

        // Gradient
        let ft = diff.dot(dc);
        let fu = -diff.dot(ds_u);
        let fv = -diff.dot(ds_v);

        // Simple gradient descent step (more robust than full Newton for 3D problems)
        let step = 0.5;
        let htt = dc.dot(dc).max(TOL);
        let huu = ds_u.dot(ds_u).max(TOL);
        let hvv = ds_v.dot(ds_v).max(TOL);

        t -= step * ft / htt;
        u -= step * fu / huu;
        v -= step * fv / hvv;

        t = t.clamp(curve_domain[0], curve_domain[1]);
        u = u.clamp(surf_domain[0], surf_domain[1]);
        v = v.clamp(surf_domain[2], surf_domain[3]);

        if ft.abs() < TOL && fu.abs() < TOL && fv.abs() < TOL {
            break;
        }
    }

    (t, u, v)
}

/// Newton refinement for surface-to-surface distance.
fn refine_surface_surface_distance(
    surf1: &Surface3, surf2: &Surface3,
    domain1: [f64; 4], domain2: [f64; 4],
    u1: f64, v1: f64, u2: f64, v2: f64,
) -> (f64, f64, f64, f64) {
    let mut u1 = u1;
    let mut v1 = v1;
    let mut u2 = u2;
    let mut v2 = v2;

    const MAX_ITER: usize = 30;
    const TOL: f64 = 1e-10;
    const STEP: f64 = 0.3;

    for _ in 0..MAX_ITER {
        let p1 = surf1.point_at(u1, v1);
        let p2 = surf2.point_at(u2, v2);

        let (du1, dv1) = surface_derivatives(surf1, u1, v1);
        let (du2, dv2) = surface_derivatives(surf2, u2, v2);

        let diff = p1 - p2;

        // Gradient
        let fu1 = diff.dot(du1);
        let fv1 = diff.dot(dv1);
        let fu2 = -diff.dot(du2);
        let fv2 = -diff.dot(dv2);

        // Simple gradient descent
        u1 -= STEP * fu1 / (du1.dot(du1) + TOL);
        v1 -= STEP * fv1 / (dv1.dot(dv1) + TOL);
        u2 -= STEP * fu2 / (du2.dot(du2) + TOL);
        v2 -= STEP * fv2 / (dv2.dot(dv2) + TOL);

        u1 = u1.clamp(domain1[0], domain1[1]);
        v1 = v1.clamp(domain1[2], domain1[3]);
        u2 = u2.clamp(domain2[0], domain2[1]);
        v2 = v2.clamp(domain2[2], domain2[3]);

        if fu1.abs() < TOL && fv1.abs() < TOL && fu2.abs() < TOL && fv2.abs() < TOL {
            break;
        }
    }

    (u1, v1, u2, v2)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle3, CylindricalSurface, Line3, Plane, SphericalSurface};
    use rcad_kernel::{BRep, PrimitiveSolid};
    use std::f64::consts::PI;

    // ── Point-Point Tests ───────────────────────────────────────────────────────

    #[test]
    fn test_distance_point_point() {
        let p1 = DVec3::ZERO;
        let p2 = DVec3::new(3.0, 4.0, 0.0);
        let dist = distance_point_point(p1, p2);
        assert!((dist - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_distance_point_point_zero() {
        let p = DVec3::new(1.0, 2.0, 3.0);
        let dist = distance_point_point(p, p);
        assert!(dist.abs() < 1e-10);
    }

    // ── Point-Curve Tests ───────────────────────────────────────────────────────

    #[test]
    fn test_distance_point_curve_line() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let point = DVec3::new(0.5, 3.0, 0.0);

        let (dist, param) = distance_point_curve(point, &line);

        assert!((dist - 3.0).abs() < 1e-4);
        assert!((param - 0.5).abs() < 1e-4);
    }

    #[test]
    fn test_distance_point_curve_circle() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        let point = DVec3::new(3.0, 0.0, 0.0);

        let (dist, _param) = distance_point_curve(point, &circle);

        // Distance should be 2.0 (3.0 - 1.0)
        assert!((dist - 2.0).abs() < 1e-4);
    }

    #[test]
    fn test_closest_point_on_curve() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X.normalize(),
        });
        let point = DVec3::new(2.5, 5.0, 0.0);

        let (param, closest) = closest_point_on_curve(&line, point);

        assert!((param - 2.5).abs() < 1e-4);
        assert!((closest - DVec3::new(2.5, 0.0, 0.0)).length() < 1e-4);
    }

    // ── Point-Surface Tests ──────────────────────────────────────────────────────

    #[test]
    fn test_distance_point_surface_plane() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let point = DVec3::new(1.0, 2.0, 5.0);

        let (dist, _u, _v) = distance_point_surface(point, &plane);

        assert!((dist - 5.0).abs() < 1e-4);
        // Note: UV coordinates depend on the plane's internal parameterization
        // which uses any_perpendicular for the x-axis direction.
    }

    #[test]
    fn test_distance_point_surface_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        });
        let point = DVec3::new(5.0, 0.0, 0.0);

        let (dist, _u, _v) = distance_point_surface(point, &sphere);

        // Distance should be 3.0 (5.0 - 2.0)
        assert!((dist - 3.0).abs() < 1e-3);
    }

    #[test]
    fn test_closest_point_on_surface() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let point = DVec3::new(3.0, 4.0, 10.0);

        let (_uv, closest) = closest_point_on_surface(&plane, point);

        // The closest point should be the projection onto the plane
        // Note: UV coordinates depend on the plane's internal parameterization
        // which uses any_perpendicular for the x-axis direction.
        assert!((closest - DVec3::new(3.0, 4.0, 0.0)).length() < 1e-4);
    }

    // ── Curve-Curve Tests ────────────────────────────────────────────────────────

    #[test]
    fn test_distance_curve_curve_parallel_lines() {
        let line1 = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let line2 = Curve3::Line(Line3 {
            origin: DVec3::new(0.0, 5.0, 0.0),
            direction: DVec3::X,
        });

        let (dist, _t1, _t2) = distance_curve_curve(&line1, &line2);

        assert!((dist - 5.0).abs() < 1e-3);
    }

    #[test]
    fn test_distance_curve_curve_skew_lines() {
        let line1 = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let line2 = Curve3::Line(Line3 {
            origin: DVec3::new(0.0, 0.0, 1.0),
            direction: DVec3::Y,
        });

        let (dist, _t1, _t2) = distance_curve_curve(&line1, &line2);

        // Skew lines: minimum distance is 1.0
        assert!((dist - 1.0).abs() < 1e-3);
    }

    // ── Curve-Surface Tests ──────────────────────────────────────────────────────

    #[test]
    fn test_distance_curve_surface_line_plane() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::new(0.0, 0.0, 5.0),
            direction: DVec3::X,
        });
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let (dist, _t, _u, _v) = distance_curve_surface(&line, &plane);

        assert!((dist - 5.0).abs() < 1e-3);
    }

    // ── Surface-Surface Tests ────────────────────────────────────────────────────

    #[test]
    fn test_distance_surface_surface_parallel_planes() {
        let plane1 = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let plane2 = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 0.0, 3.0),
            normal: DVec3::Z,
        });

        let (dist, _u1, _v1, _u2, _v2) = distance_surface_surface(&plane1, &plane2);

        assert!((dist - 3.0).abs() < 1e-3);
    }

    #[test]
    fn test_distance_surface_surface_sphere_plane() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(0.0, 0.0, 5.0),
            axis: DVec3::Z,
            radius: 2.0,
        });
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let (dist, _u1, _v1, _u2, _v2) = distance_surface_surface(&sphere, &plane);

        // Distance should be 3.0 (5.0 - 2.0)
        assert!((dist - 3.0).abs() < 0.02);
    }

    // ── BRep-BRep Tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_distance_brep_brep_boxes() {
        let box1 = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let box2 = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let (dist, _p1, _p2) = distance_brep_brep(&box1, &box2);

        // Both boxes at origin, should be 0 distance (they overlap)
        assert!(dist < 0.1);
    }

    // ── Find Closest Points Tests ────────────────────────────────────────────────

    #[test]
    fn test_find_closest_points_circle() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        let point = DVec3::new(0.0, 0.0, 0.0); // Center of circle

        let closest = find_closest_points(&circle, point, 2);

        // For a point at the center, all points are equally distant
        // Should return valid parameters
        assert!(!closest.is_empty());
        for (_, dist) in &closest {
            assert!((dist - 1.0).abs() < 1e-3);
        }
    }

    // ── Find Furthest Points Tests ───────────────────────────────────────────────

    #[test]
    fn test_find_furthest_points_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        let (p1, p2) = find_furthest_points(&brep, DVec3::X);

        // Should find opposite corners along X direction
        let proj1 = p1.dot(DVec3::X);
        let proj2 = p2.dot(DVec3::X);

        // One should be minimum, one should be maximum
        assert!((proj2 - proj1).abs() > 1.0);
    }

    // ── Support Shapes Tests ─────────────────────────────────────────────────────

    #[test]
    fn test_find_supporting_face() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Point on the +Z face
        let point = DVec3::new(0.5, 0.5, 0.5);

        let face = find_supporting_face(&brep, point);
        assert!(face.is_some());
    }

    #[test]
    fn test_find_supporting_edge() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Point near an edge
        let point = DVec3::new(0.5, 0.0, 0.5);

        let edge = find_supporting_edge(&brep, point);
        assert!(edge.is_some());
    }

    // ── Newton Refinement Tests ─────────────────────────────────────────────────

    #[test]
    fn test_refinement_convergence() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::new(5.0, 5.0, 0.0),
            normal: DVec3::Z,
            radius: 3.0,
        });
        let point = DVec3::new(10.0, 5.0, 0.0);

        let (dist, _param) = distance_point_curve(point, &circle);

        // Distance should be 2.0 (10 - 5 - 3)
        assert!((dist - 2.0).abs() < 1e-3);
        // Note: The exact parameter value depends on the circle's internal parameterization
        // which uses any_perpendicular for the x-axis direction.
        // For a circle with Z normal, the parameterization is:
        // point(t) = center + radius * (cos(t) * Y + sin(t) * (-X))
        // So the "rightmost" point corresponds to t = 3*pi/2
    }
}
