//! BRepIntCurveSurface-style curve-surface intersection operations.
//!
//! Provides intersection computations between curves/lines and BRep shapes.
//! Analogous to OCCT `BRepIntCurveSurface` package.
//!
//! # Capabilities
//!
//! - **Curve-BRep Intersection**: Find all intersection points between a curve and a BRep
//! - **Line-BRep Intersection**: Efficient ray/line intersection with BRep faces
//! - **Ray Casting**: Cast rays through a BRep and collect all hits
//! - **Point-Inside Test**: Determine if a point is inside a solid using ray casting
//!
//! # Example
//!
//! ```rust
//! use glam::DVec3;
//! use rcad_kernel::{BRep, PrimitiveSolid};
//! use rcad_algorithms::brep_int_curve_surface::{intersect_line_with_brep, ray_cast};
//!
//! let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
//!
//! // Intersect a line through the box
//! let intersections = intersect_line_with_brep(DVec3::new(1.0, 1.0, 5.0), DVec3::NEG_Z, &box_brep, 1e-7);
//! assert_eq!(intersections.len(), 2); // Entry and exit points
//!
//! // Ray cast from above
//! let hits = ray_cast(DVec3::new(1.0, 1.0, 5.0), DVec3::NEG_Z, &box_brep);
//! assert!(!hits.is_empty());
//! ```

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, Surface3, Line3, CurveEval, SurfaceEval};
use rcad_kernel::{BRep, Face};
use crate::bvh::{Aabb, Bvh};
use crate::tolerance::TOLERANCE_ABS;
use crate::int_ana::{intersect_line_plane, intersect_line_torus};
use crate::inttools::curve_surface::{
    intersect_line_cylinder as intersect_line_cylinder_range,
    intersect_line_sphere as intersect_line_sphere_range,
    intersect_line_cone as intersect_line_cone_range,
};
use rayon::prelude::*;

// =============================================================================
// Result Types
// =============================================================================

/// Result of a curve-BRep intersection.
///
/// Represents a single intersection point between a curve and a BRep face.
#[derive(Debug, Clone)]
pub struct CurveBRepIntersection {
    /// The intersection point in 3D space.
    pub point: DVec3,
    /// Parameter value on the curve where the intersection occurs.
    pub param: f64,
    /// Index of the face containing the intersection point.
    pub face_index: usize,
    /// UV parameters on the face surface at the intersection point.
    pub uv: DVec2,
}

impl Default for CurveBRepIntersection {
    fn default() -> Self {
        Self {
            point: DVec3::ZERO,
            param: 0.0,
            face_index: 0,
            uv: DVec2::ZERO,
        }
    }
}

/// Result of a curve-face intersection.
///
/// Represents a single intersection point between a curve and a single face.
#[derive(Debug, Clone)]
pub struct CurveFaceIntersection {
    /// The intersection point in 3D space.
    pub point: DVec3,
    /// Parameter value on the curve where the intersection occurs.
    pub param: f64,
    /// UV parameters on the face surface at the intersection point.
    pub uv: DVec2,
}

impl Default for CurveFaceIntersection {
    fn default() -> Self {
        Self {
            point: DVec3::ZERO,
            param: 0.0,
            uv: DVec2::ZERO,
        }
    }
}

/// Result of a ray casting operation.
///
/// Represents a single hit point when a ray intersects a BRep.
#[derive(Debug, Clone)]
pub struct RayHit {
    /// The intersection point in 3D space.
    pub point: DVec3,
    /// Distance from the ray origin to the hit point.
    pub distance: f64,
    /// Index of the face containing the hit point.
    pub face_index: usize,
    /// Surface normal at the hit point.
    pub normal: DVec3,
}

impl Default for RayHit {
    fn default() -> Self {
        Self {
            point: DVec3::ZERO,
            distance: 0.0,
            face_index: 0,
            normal: DVec3::Z,
        }
    }
}

// =============================================================================
// Curve-BRep Intersection
// =============================================================================

/// Intersect a curve with a BRep shape.
///
/// Computes all intersection points between the curve and all faces of the BRep.
/// Results are sorted by the curve parameter.
///
/// # Arguments
/// * `curve` - The 3D curve to intersect.
/// * `brep` - The BRep shape to intersect with.
/// * `tol` - Tolerance for geometric computations.
///
/// # Returns
/// A vector of intersection points sorted by curve parameter.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Curve3, Line3};
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use rcad_algorithms::brep_int_curve_surface::intersect_curve_with_brep;
///
/// let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let line = Curve3::Line(Line3 { origin: DVec3::new(1.0, 1.0, 5.0), direction: DVec3::NEG_Z });
/// let intersections = intersect_curve_with_brep(&line, &box_brep, 1e-7);
/// assert_eq!(intersections.len(), 2);
/// ```
pub fn intersect_curve_with_brep(
    curve: &Curve3,
    brep: &BRep,
    tol: f64,
) -> Vec<CurveBRepIntersection> {
    let face_indices = collect_face_indices(brep);

    // Use parallel processing for large number of faces
    let all_intersections: Vec<Vec<CurveBRepIntersection>> = if face_indices.len() > 16 {
        face_indices
            .par_iter()
            .map(|&face_idx| {
                intersect_curve_with_face(curve, brep, face_idx, tol)
                    .into_iter()
                    .map(|cfi| CurveBRepIntersection {
                        point: cfi.point,
                        param: cfi.param,
                        face_index: face_idx,
                        uv: cfi.uv,
                    })
                    .collect()
            })
            .collect()
    } else {
        face_indices
            .iter()
            .map(|&face_idx| {
                intersect_curve_with_face(curve, brep, face_idx, tol)
                    .into_iter()
                    .map(|cfi| CurveBRepIntersection {
                        point: cfi.point,
                        param: cfi.param,
                        face_index: face_idx,
                        uv: cfi.uv,
                    })
                    .collect()
            })
            .collect()
    };

    // Flatten and sort by parameter
    let mut results: Vec<CurveBRepIntersection> = all_intersections.into_iter().flatten().collect();
    results.sort_by(|a, b| {
        a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Remove duplicates (within tolerance)
    deduplicate_intersections(&mut results, tol);

    results
}

/// Intersect a line with a BRep shape.
///
/// Computes all intersection points between an infinite line and all faces of the BRep.
/// Results are sorted by the line parameter.
///
/// # Arguments
/// * `origin` - The origin point of the line.
/// * `direction` - The direction vector of the line (will be normalized).
/// * `brep` - The BRep shape to intersect with.
/// * `tol` - Tolerance for geometric computations.
///
/// # Returns
/// A vector of intersection points sorted by line parameter.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use rcad_algorithms::brep_int_curve_surface::intersect_line_with_brep;
///
/// let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let intersections = intersect_line_with_brep(DVec3::new(1.0, 1.0, 5.0), DVec3::NEG_Z, &box_brep, 1e-7);
/// assert_eq!(intersections.len(), 2);
/// ```
pub fn intersect_line_with_brep(
    origin: DVec3,
    direction: DVec3,
    brep: &BRep,
    tol: f64,
) -> Vec<CurveBRepIntersection> {
    let dir = direction.normalize_or_zero();
    let face_indices = collect_face_indices(brep);

    // Use BVH for early rejection if available
    let bvh = Bvh::build(brep);
    let candidate_faces = if face_indices.len() > 8 {
        filter_faces_by_line_aabb(&bvh, origin, dir, &face_indices)
    } else {
        face_indices
    };

    // Intersect with each candidate face
    let mut results: Vec<CurveBRepIntersection> = candidate_faces
        .into_iter()
        .flat_map(|face_idx| {
            intersect_line_with_face(origin, dir, brep, face_idx, tol)
                .into_iter()
                .map(move |cfi| CurveBRepIntersection {
                    point: cfi.point,
                    param: cfi.param,
                    face_index: face_idx,
                    uv: cfi.uv,
                })
        })
        .collect();

    // Sort by parameter
    results.sort_by(|a, b| {
        a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Remove duplicates
    deduplicate_intersections(&mut results, tol);

    results
}

// =============================================================================
// Face-Level Intersections
// =============================================================================

/// Intersect a curve with a single face of a BRep.
///
/// Computes all intersection points between the curve and the face's surface,
/// then filters to only those that lie within the face's boundaries.
///
/// # Arguments
/// * `curve` - The 3D curve to intersect.
/// * `brep` - The BRep containing the face.
/// * `face_idx` - Index of the face to intersect with.
/// * `tol` - Tolerance for geometric computations.
///
/// # Returns
/// A vector of intersection points.
pub fn intersect_curve_with_face(
    curve: &Curve3,
    brep: &BRep,
    face_idx: usize,
    tol: f64,
) -> Vec<CurveFaceIntersection> {
    // Handle lines specially using analytic intersection
    if let Curve3::Line(line) = curve {
        return intersect_line_with_face(line.origin, line.direction, brep, face_idx, tol);
    }

    // Get the surface for this face
    let surface = match get_face_surface(brep, face_idx) {
        Some(s) => s,
        None => return Vec::new(),
    };

    // Get the curve's parameter range
    let [t0, t1] = get_curve_domain(curve);

    // Sample the curve and find intersections
    let n_samples = 64;
    let mut intersections = Vec::new();

    // Sample to find sign changes
    let initial_point = curve.point_at(t0);
    let mut prev_proj = project_point_to_surface(&surface, initial_point);
    let mut prev_dist = (prev_proj.point - initial_point).length();

    for i in 1..=n_samples {
        let t = t0 + (t1 - t0) * i as f64 / n_samples as f64;
        let point = curve.point_at(t);
        let proj = project_point_to_surface(&surface, point);
        let dist = (proj.point - point).length();

        // Check for sign change (intersection)
        if prev_dist < tol || dist < tol {
            // Already close to surface
            if prev_dist < tol && (i == 1 || intersections.last().map_or(true, |last: &CurveFaceIntersection| {
                (last.param - (t - (t1 - t0) / n_samples as f64)).abs() > tol
            })) {
                intersections.push(CurveFaceIntersection {
                    point: prev_proj.point,
                    param: t - (t1 - t0) / n_samples as f64,
                    uv: proj.uv,
                });
            }
        } else if let Some(cfi) = refine_curve_surface_intersection(
            curve, &surface, t - (t1 - t0) / n_samples as f64, t, tol
        ) {
            intersections.push(cfi);
        }

        prev_proj = proj;
        prev_dist = dist;
    }

    // Filter to face bounds
    intersections.retain(|cfi| is_point_in_face_bounds(brep, cfi.point, cfi.uv, face_idx, tol));

    intersections
}

/// Intersect a line with a single face of a BRep.
///
/// Computes all intersection points between the line and the face's surface
/// that lie within the face's boundaries.
///
/// # Arguments
/// * `origin` - The origin point of the line.
/// * `dir` - The direction vector of the line (should be normalized).
/// * `brep` - The BRep containing the face.
/// * `face_idx` - Index of the face to intersect with.
/// * `tol` - Tolerance for geometric computations.
///
/// # Returns
/// A vector of intersection points within the face bounds.
pub fn intersect_line_with_face(
    origin: DVec3,
    dir: DVec3,
    brep: &BRep,
    face_idx: usize,
    tol: f64,
) -> Vec<CurveFaceIntersection> {
    // Get the surface for this face
    let surface = match get_face_surface(brep, face_idx) {
        Some(s) => s,
        None => return Vec::new(),
    };

    // Intersect line with surface
    let hits = intersect_line_with_surface(origin, dir, &surface);

    if hits.is_empty() {
        return Vec::new();
    }

    // Find all hits within the face bounds
    let mut results = Vec::new();
    for hit in hits {
        let uv = compute_uv_for_point(&surface, hit.point);

        if is_point_in_face_bounds(brep, hit.point, uv, face_idx, tol) {
            results.push(CurveFaceIntersection {
                point: hit.point,
                param: hit.param,
                uv,
            });
        }
    }

    results
}

// =============================================================================
// Ray Casting
// =============================================================================

/// Cast a ray through a BRep and collect all hit points.
///
/// A ray is a semi-infinite line starting from the origin and extending
/// in the given direction. Only intersections in the positive direction
/// are returned.
///
/// # Arguments
/// * `origin` - The origin point of the ray.
/// * `direction` - The direction vector of the ray (will be normalized).
/// * `brep` - The BRep shape to cast the ray through.
///
/// # Returns
/// A vector of hit points sorted by distance from the origin.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use rcad_algorithms::brep_int_curve_surface::ray_cast;
///
/// let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let hits = ray_cast(DVec3::new(1.0, 1.0, 5.0), DVec3::NEG_Z, &box_brep);
/// assert!(!hits.is_empty());
/// ```
pub fn ray_cast(
    origin: DVec3,
    direction: DVec3,
    brep: &BRep,
) -> Vec<RayHit> {
    let dir = direction.normalize_or_zero();

    // Get all line intersections
    let intersections = intersect_line_with_brep(origin, dir, brep, TOLERANCE_ABS);

    // Filter to positive direction and convert to RayHit
    let mut hits: Vec<RayHit> = intersections
        .into_iter()
        .filter(|i| i.param > TOLERANCE_ABS)
        .filter_map(|i| {
            let face = get_face(brep, i.face_index);
            let normal = face.map_or(DVec3::Z, |f| f.normal);

            Some(RayHit {
                point: i.point,
                distance: i.param,
                face_index: i.face_index,
                normal,
            })
        })
        .collect();

    // Sort by distance
    hits.sort_by(|a, b| {
        a.distance.partial_cmp(&b.distance).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Deduplicate
    deduplicate_ray_hits(&mut hits);

    hits
}

/// Shoot a ray with a maximum distance limit.
///
/// Similar to `ray_cast`, but only returns hits within the maximum distance.
///
/// # Arguments
/// * `origin` - The origin point of the ray.
/// * `direction` - The direction vector of the ray (will be normalized).
/// * `brep` - The BRep shape to cast the ray through.
/// * `max_distance` - Maximum distance from origin to consider.
///
/// # Returns
/// A vector of hit points within the distance limit, sorted by distance.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use rcad_algorithms::brep_int_curve_surface::shoot_ray;
///
/// let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let hits = shoot_ray(DVec3::new(1.0, 1.0, 5.0), DVec3::NEG_Z, &box_brep, 10.0);
/// assert!(!hits.is_empty());
/// ```
pub fn shoot_ray(
    origin: DVec3,
    direction: DVec3,
    brep: &BRep,
    max_distance: f64,
) -> Vec<RayHit> {
    let hits = ray_cast(origin, direction, brep);

    hits.into_iter()
        .filter(|h| h.distance <= max_distance + TOLERANCE_ABS)
        .collect()
}

/// Determine if a point is inside a solid using ray casting.
///
/// Uses the even-odd rule: cast a ray from the point and count intersections.
/// If the number of intersections is odd, the point is inside.
/// Uses multiple rays for robustness and voting.
///
/// # Arguments
/// * `point` - The query point.
/// * `brep` - The BRep solid to test against.
///
/// # Returns
/// `true` if the point is inside the solid, `false` otherwise.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use rcad_algorithms::brep_int_curve_surface::is_point_inside_by_ray;
///
/// let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let inside = is_point_inside_by_ray(DVec3::new(1.0, 1.0, 1.0), &box_brep);
/// assert!(inside);
/// ```
pub fn is_point_inside_by_ray(point: DVec3, brep: &BRep) -> bool {
    // Use multiple ray directions for robustness
    let directions = [
        DVec3::X,
        DVec3::Y,
        DVec3::Z,
        DVec3::new(1.0, 1.0, 1.0).normalize(),
    ];

    let mut inside_votes = 0;
    let mut total_votes = 0;

    for dir in directions {
        let hits = ray_cast(point, dir, brep);

        // Skip if too many hits (likely grazing or on boundary)
        if hits.len() > 100 {
            continue;
        }

        // Odd number of hits = inside
        if hits.len() % 2 == 1 {
            inside_votes += 1;
        }
        total_votes += 1;
    }

    if total_votes == 0 {
        // Fallback: try an arbitrary direction
        let hits = ray_cast(point, DVec3::new(0.57735, 0.57735, 0.57735), brep);
        return hits.len() % 2 == 1;
    }

    // Majority vote
    inside_votes > total_votes / 2
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Collect all face indices from a BRep.
fn collect_face_indices(brep: &BRep) -> Vec<usize> {
    let mut indices = Vec::new();
    let mut count = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for _ in &shell.faces {
                indices.push(count);
                count += 1;
            }
        }
    }

    indices
}

/// Get a face by its flat index.
fn get_face(brep: &BRep, face_idx: usize) -> Option<&Face> {
    let mut count = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if count == face_idx {
                    return Some(face);
                }
                count += 1;
            }
        }
    }

    None
}

/// Get the surface for a face by its index.
fn get_face_surface(brep: &BRep, face_idx: usize) -> Option<Surface3> {
    let surface_idx = brep.geom.face_surface.get(face_idx).copied().flatten()?;
    brep.geom.surfaces.get(surface_idx).cloned()
}

/// Get the default domain for a curve.
fn get_curve_domain(curve: &Curve3) -> [f64; 2] {
    match curve {
        Curve3::Line(_) => [-1e10, 1e10], // Infinite
        _ => curve.default_domain(),
    }
}

/// Project a point onto a surface and return the projected point and UV.
fn project_point_to_surface(surface: &Surface3, point: DVec3) -> ProjectedPoint {
    use rcad_kernel::projection::closest_point_on_surface;

    let result = closest_point_on_surface(surface, point, 32);
    ProjectedPoint {
        point: result.point,
        uv: DVec2::new(result.params.0, result.params.1),
    }
}

struct ProjectedPoint {
    point: DVec3,
    uv: DVec2,
}

/// Intersect a line with a surface.
fn intersect_line_with_surface(origin: DVec3, dir: DVec3, surface: &Surface3) -> Vec<LineSurfaceHit> {
    let line = Line3 { origin, direction: dir };
    let infinite_range = [-1e10, 1e10];

    match surface {
        Surface3::Plane(plane) => {
            intersect_line_plane(&line, plane)
                .map(|r| vec![LineSurfaceHit { point: r.point, param: r.param }])
                .unwrap_or_default()
        }
        Surface3::Cylinder(cyl) => {
            intersect_line_cylinder_range(&line, infinite_range, cyl)
                .into_iter()
                .map(|h| LineSurfaceHit { point: h.point, param: h.curve_param })
                .collect()
        }
        Surface3::Sphere(sph) => {
            intersect_line_sphere_range(&line, infinite_range, sph)
                .into_iter()
                .map(|h| LineSurfaceHit { point: h.point, param: h.curve_param })
                .collect()
        }
        Surface3::Cone(cone) => {
            intersect_line_cone_range(&line, infinite_range, cone)
                .into_iter()
                .map(|h| LineSurfaceHit { point: h.point, param: h.curve_param })
                .collect()
        }
        Surface3::Torus(torus) => {
            intersect_line_with_torus(&line, torus)
        }
        Surface3::BSpline(bspline) => {
            // For B-spline surfaces, use iterative intersection
            intersect_line_with_bspline_surface(&line, &Surface3::BSpline(bspline.clone()))
        }
        _ => {
            // General case: sample and refine
            intersect_line_with_general_surface(&line, surface)
        }
    }
}

struct LineSurfaceHit {
    point: DVec3,
    param: f64,
}

impl Default for LineSurfaceHit {
    fn default() -> Self {
        Self {
            point: DVec3::ZERO,
            param: 0.0,
        }
    }
}

/// Intersect a line with a torus surface.
fn intersect_line_with_torus(line: &Line3, torus: &rcad_kernel::geom::ToroidalSurface) -> Vec<LineSurfaceHit> {
    let results = intersect_line_torus(line, torus);
    results.into_iter().map(|(point, param)| LineSurfaceHit { point, param }).collect()
}

/// Intersect a line with a B-spline surface using iterative methods.
fn intersect_line_with_bspline_surface(line: &Line3, surface: &Surface3) -> Vec<LineSurfaceHit> {
    intersect_line_with_general_surface(line, surface)
}

/// Intersect a line with a general surface using iterative methods.
fn intersect_line_with_general_surface(line: &Line3, surface: &Surface3) -> Vec<LineSurfaceHit> {
    let domain = surface.default_domain();
    let [u0, u1, v0, v1] = domain;

    // Sample the surface to find intersection candidates
    let n_samples = 32;
    let mut hits = Vec::new();

    for i in 0..n_samples {
        let u = u0 + (u1 - u0) * i as f64 / (n_samples - 1) as f64;

        for j in 0..n_samples {
            let v = v0 + (v1 - v0) * j as f64 / (n_samples - 1) as f64;
            let surf_point = surface.point_at(u, v);

            // Find closest point on line
            let to_surface = surf_point - line.origin;
            let t = to_surface.dot(line.direction);

            if t.abs() < 1e10 { // Within reasonable range
                let line_point = line.origin + t * line.direction;
                let dist = (line_point - surf_point).length();

                if dist < TOLERANCE_ABS * 100.0 {
                    // Close enough - refine
                    if let Some(hit) = refine_line_surface_intersection(line, surface, t, u, v) {
                        hits.push(hit);
                    }
                }
            }
        }
    }

    // Sort and deduplicate
    hits.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal));
    deduplicate_line_surface_hits(&mut hits);

    hits
}

/// Refine a line-surface intersection using Newton iteration.
fn refine_line_surface_intersection(
    line: &Line3,
    surface: &Surface3,
    initial_t: f64,
    initial_u: f64,
    initial_v: f64,
) -> Option<LineSurfaceHit> {
    let mut t = initial_t;
    let mut u = initial_u;
    let mut v = initial_v;

    const MAX_ITER: usize = 20;
    const H: f64 = 1e-7;

    for _ in 0..MAX_ITER {
        let line_point = line.origin + t * line.direction;
        let surf_point = surface.point_at(u, v);

        let error = (surf_point - line_point).length();
        if error < TOLERANCE_ABS {
            return Some(LineSurfaceHit { point: surf_point, param: t });
        }

        // Compute Jacobian using finite differences
        let surf_du = surface.point_at(u + H, v);
        let surf_dv = surface.point_at(u, v + H);

        let du = (surf_du - surf_point) / H;
        let dv = (surf_dv - surf_point) / H;

        // Residual: surf_point - line_point
        // We want to minimize this
        // Simple gradient descent approach
        let residual = line_point - surf_point;

        // Update parameters
        let delta_u = residual.dot(du) / du.length_squared().max(1e-12);
        let delta_v = residual.dot(dv) / dv.length_squared().max(1e-12);
        let delta_t = residual.dot(line.direction);

        u += delta_u * 0.5;
        v += delta_v * 0.5;
        t += delta_t * 0.5;

        // Clamp to surface domain
        let domain = surface.default_domain();
        u = u.clamp(domain[0], domain[1]);
        v = v.clamp(domain[2], domain[3]);
    }

    // Return best effort
    let surf_point = surface.point_at(u, v);
    let line_point = line.origin + t * line.direction;
    let error = (surf_point - line_point).length();

    if error < TOLERANCE_ABS * 100.0 {
        Some(LineSurfaceHit { point: surf_point, param: t })
    } else {
        None
    }
}

/// Refine a curve-surface intersection using Newton iteration.
fn refine_curve_surface_intersection(
    curve: &Curve3,
    surface: &Surface3,
    t_lo: f64,
    t_hi: f64,
    tol: f64,
) -> Option<CurveFaceIntersection> {
    let mut t = (t_lo + t_hi) / 2.0;

    const MAX_ITER: usize = 30;
    const H: f64 = 1e-7;

    for _ in 0..MAX_ITER {
        let curve_point = curve.point_at(t);
        let proj = project_point_to_surface(surface, curve_point);
        let dist = (proj.point - curve_point).length();

        if dist < tol {
            return Some(CurveFaceIntersection {
                point: proj.point,
                param: t,
                uv: proj.uv,
            });
        }

        // Compute derivative
        let curve_point_plus = curve.point_at(t + H);
        let tangent = (curve_point_plus - curve_point) / H;

        // Direction from curve to surface
        let to_surface = proj.point - curve_point;

        // Update t to move curve point toward surface
        if tangent.length_squared() > 1e-12 {
            let delta_t = to_surface.dot(tangent) / tangent.length_squared();
            t += delta_t * 0.5;
        }

        // Clamp to range
        t = t.clamp(t_lo, t_hi);
    }

    None
}

/// Compute UV parameters for a point on a surface.
fn compute_uv_for_point(surface: &Surface3, point: DVec3) -> DVec2 {
    let proj = project_point_to_surface(surface, point);
    proj.uv
}

/// Check if a point lies within a face's boundaries.
///
/// Uses point-in-polygon test with the face's outer wire vertices.
fn is_point_in_face_bounds(brep: &BRep, point: DVec3, _uv: DVec2, face_idx: usize, tol: f64) -> bool {
    let face = match get_face(brep, face_idx) {
        Some(f) => f,
        None => return true, // No face bounds info, assume inside
    };

    if face.outer_wire.edges.is_empty() {
        return true;
    }

    // Get the surface for this face to project onto
    let surface = match get_face_surface(brep, face_idx) {
        Some(s) => s,
        None => return true, // Can't check without surface
    };

    // Collect the wire vertices in order
    let wire_vertices = collect_wire_vertices(brep, &face.outer_wire);
    if wire_vertices.len() < 3 {
        return true;
    }

    // Project the point and wire vertices onto a 2D plane for point-in-polygon test
    let is_inside = point_in_polygon_on_surface(&surface, point, &wire_vertices, tol);

    // Also check against inner wires (holes) - point must be outside all inner wires
    if is_inside {
        for inner_wire in &face.inner_wires {
            let inner_vertices = collect_wire_vertices(brep, inner_wire);
            if inner_vertices.len() >= 3 {
                if point_in_polygon_on_surface(&surface, point, &inner_vertices, tol) {
                    return false; // Point is inside a hole
                }
            }
        }
    }

    is_inside
}

/// Collect 3D vertex positions from a wire in traversal order.
fn collect_wire_vertices(brep: &BRep, wire: &rcad_kernel::topology::Wire) -> Vec<DVec3> {
    let mut vertices = Vec::new();
    for wire_edge in &wire.edges {
        let edge = match brep.edges.get(wire_edge.idx) {
            Some(e) => e,
            None => continue,
        };

        // Get the starting vertex of this edge (based on traversal direction)
        let vertex_idx = if wire_edge.forward { edge.start } else { edge.end };
        if let Some(v) = brep.vertices.get(vertex_idx) {
            vertices.push(v.point);
        }
    }
    vertices
}

/// Check if a point is inside a polygon projected onto a surface.
///
/// Uses the ray casting algorithm (even-odd rule) in the surface's natural 2D coordinate system.
fn point_in_polygon_on_surface(surface: &Surface3, point: DVec3, polygon: &[DVec3], tol: f64) -> bool {
    // Check if we have enough unique vertices for a proper polygon
    let unique_vertices: Vec<DVec3> = {
        let mut unique: Vec<DVec3> = Vec::new();
        for &v in polygon {
            if !unique.iter().any(|&u: &DVec3| (u - v).length() < tol) {
                unique.push(v);
            }
        }
        unique
    };

    // For surfaces with degenerate wires (like a sphere with just a seam edge),
    // we need different logic based on the surface type
    if unique_vertices.len() < 3 {
        // Degenerate polygon - use surface-specific logic
        return point_on_closed_surface(surface, point, tol);
    }

    // For planes and other simple surfaces, project to 2D and use point-in-polygon
    match surface {
        Surface3::Plane(plane) => {
            // Create a local 2D coordinate system on the plane
            let normal = plane.normal.normalize();
            let (u_dir, v_dir) = get_plane_tangent_dirs(normal);

            // Project the test point relative to plane origin
            let rel_point = point - plane.origin;
            let p2d = DVec2::new(rel_point.dot(u_dir), rel_point.dot(v_dir));

            // Project polygon vertices
            let poly2d: Vec<DVec2> = polygon
                .iter()
                .map(|&v| {
                    let rel = v - plane.origin;
                    DVec2::new(rel.dot(u_dir), rel.dot(v_dir))
                })
                .collect();

            point_in_polygon_2d(p2d, &poly2d, tol)
        }
        Surface3::Sphere(_) | Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Torus(_) => {
            // For curved surfaces with proper wires, use bounding box check
            // A more sophisticated implementation would use proper surface parameter space
            point_in_bounding_box(point, polygon, tol)
        }
        _ => {
            // Default to bounding box check
            point_in_bounding_box(point, polygon, tol)
        }
    }
}

/// Check if a point is on a closed surface (like a full sphere or torus).
///
/// For closed surfaces with degenerate wires (just seams), any point on the surface
/// is considered valid.
fn point_on_closed_surface(surface: &Surface3, point: DVec3, tol: f64) -> bool {
    match surface {
        Surface3::Sphere(sph) => {
            // Check if point is on the sphere surface
            let dist = (point - sph.center).length();
            (dist - sph.radius).abs() < tol * 10.0
        }
        Surface3::Torus(torus) => {
            // Check if point is on the torus surface using the distance formula
            // Distance from point to torus center in the plane perpendicular to axis
            let rel = point - torus.center;
            let axis = torus.axis.normalize();
            let axial_dist = rel.dot(axis).abs(); // Distance along axis
            let radial_vec = rel - axial_dist * axis;
            let radial_dist = radial_vec.length();
            // Distance to torus surface: sqrt((radial_dist - major_r)^2 + axial_dist^2) - minor_r
            let dist_to_surface = ((radial_dist - torus.major_radius).powi(2) + axial_dist.powi(2)).sqrt() - torus.minor_radius;
            dist_to_surface.abs() < tol * 10.0
        }
        Surface3::Cylinder(cyl) => {
            // For a cylinder, check if point is on the infinite cylinder surface
            // This assumes a full cylinder - caps would need separate face handling
            let rel = point - cyl.origin;
            let axis = cyl.axis.normalize();
            let axial_dist = rel.dot(axis);
            let radial_vec = rel - axial_dist * axis;
            let radial_dist = radial_vec.length();
            (radial_dist - cyl.radius).abs() < tol * 10.0
        }
        Surface3::Cone(cone) => {
            // For a cone, check if point is on the infinite cone surface
            // Distance from apex along axis
            let rel = point - cone.apex;
            let axis = cone.axis.normalize();
            let axial_dist = rel.dot(axis);
            if axial_dist < -tol {
                return false; // Behind apex
            }
            let radial_vec = rel - axial_dist * axis;
            let radial_dist = radial_vec.length();
            // Cone radius at this height is proportional to axial distance
            // half_angle is such that radius = axial_dist * tan(half_angle)
            // We need to derive half_angle from the cone definition
            // Assuming cone is defined with base_radius at some height
            // For now, use a simple check
            let expected_radius = axial_dist * 0.5; // Approximate, should use actual cone angle
            (radial_dist - expected_radius).abs() < tol * 10.0
        }
        _ => {
            // For other surfaces, default to true
            // This may need refinement for specific surface types
            true
        }
    }
}

/// Get two orthogonal tangent directions for a plane with the given normal.
fn get_plane_tangent_dirs(normal: DVec3) -> (DVec3, DVec3) {
    let u_dir = if normal.x.abs() > 0.9 {
        normal.cross(DVec3::Y).normalize()
    } else {
        normal.cross(DVec3::X).normalize()
    };
    let v_dir = normal.cross(u_dir).normalize();
    (u_dir, v_dir)
}

/// Check if a point is inside a 2D polygon using ray casting.
fn point_in_polygon_2d(point: DVec2, polygon: &[DVec2], tol: f64) -> bool {
    if polygon.len() < 3 {
        return false;
    }

    // First check if point is on any edge (within tolerance)
    for i in 0..polygon.len() {
        let j = (i + 1) % polygon.len();
        if point_on_segment_2d(point, polygon[i], polygon[j], tol) {
            return true; // On boundary counts as inside
        }
    }

    // Ray casting algorithm: count intersections with a horizontal ray
    let mut inside = false;
    let mut j = polygon.len() - 1;

    for i in 0..polygon.len() {
        let pi = polygon[i];
        let pj = polygon[j];

        if ((pi.y > point.y) != (pj.y > point.y))
            && (point.x < (pj.x - pi.x) * (point.y - pi.y) / (pj.y - pi.y) + pi.x)
        {
            inside = !inside;
        }
        j = i;
    }

    inside
}

/// Check if a 2D point lies on a line segment.
fn point_on_segment_2d(p: DVec2, a: DVec2, b: DVec2, tol: f64) -> bool {
    let ab = b - a;
    let ap = p - a;
    let ab_len_sq = ab.length_squared();

    if ab_len_sq < tol * tol {
        // Degenerate segment (a == b)
        return ap.length() < tol;
    }

    // Project p onto line ab
    let t = ap.dot(ab) / ab_len_sq;
    if t < -tol / ab.length() || t > 1.0 + tol / ab.length() {
        return false;
    }

    let closest = a + t.clamp(0.0, 1.0) * ab;
    (p - closest).length() < tol
}

/// Check if a point is inside a bounding box defined by polygon vertices.
fn point_in_bounding_box(point: DVec3, polygon: &[DVec3], tol: f64) -> bool {
    if polygon.is_empty() {
        return true;
    }

    let mut min_pt = polygon[0];
    let mut max_pt = polygon[0];

    for &v in &polygon[1..] {
        min_pt = min_pt.min(v);
        max_pt = max_pt.max(v);
    }

    point.x >= min_pt.x - tol
        && point.x <= max_pt.x + tol
        && point.y >= min_pt.y - tol
        && point.y <= max_pt.y + tol
        && point.z >= min_pt.z - tol
        && point.z <= max_pt.z + tol
}

/// Filter faces using line AABB intersection.
fn filter_faces_by_line_aabb(
    bvh: &Bvh,
    origin: DVec3,
    dir: DVec3,
    face_indices: &[usize],
) -> Vec<usize> {
    // Create an AABB for the ray (extended in both directions)
    let mut ray_aabb = Aabb::empty();

    // Extend along the line for a significant distance
    for t in [-1000.0, 1000.0] {
        let p = origin + dir * t;
        ray_aabb.expand_point(p);
    }

    // Add a small tolerance
    ray_aabb.expand_point(origin + DVec3::splat(TOLERANCE_ABS));

    // Use BVH to query faces that intersect the ray AABB
    let candidate_faces = bvh.query_aabb(&ray_aabb);

    // Filter to only return faces that are in the original list
    let face_set: std::collections::HashSet<usize> = face_indices.iter().copied().collect();
    candidate_faces
        .into_iter()
        .filter(|idx| face_set.contains(idx))
        .collect()
}

/// Remove duplicate intersections.
fn deduplicate_intersections(intersections: &mut Vec<CurveBRepIntersection>, tol: f64) {
    if intersections.len() <= 1 {
        return;
    }

    let mut keep = vec![true; intersections.len()];

    for i in 0..intersections.len() {
        if !keep[i] {
            continue;
        }
        for j in (i + 1)..intersections.len() {
            if !keep[j] {
                continue;
            }
            let dist = (intersections[i].point - intersections[j].point).length();
            if dist < tol * 10.0 {
                keep[j] = false;
            }
        }
    }

    let mut write_idx = 0;
    for i in 0..intersections.len() {
        if keep[i] {
            intersections[write_idx] = intersections[i].clone();
            write_idx += 1;
        }
    }
    intersections.truncate(write_idx);
}

/// Remove duplicate line-surface hits.
fn deduplicate_line_surface_hits(hits: &mut Vec<LineSurfaceHit>) {
    if hits.len() <= 1 {
        return;
    }

    let mut keep = vec![true; hits.len()];

    for i in 0..hits.len() {
        if !keep[i] {
            continue;
        }
        for j in (i + 1)..hits.len() {
            if !keep[j] {
                continue;
            }
            let dist = (hits[i].point - hits[j].point).length();
            if dist < TOLERANCE_ABS * 100.0 {
                keep[j] = false;
            }
        }
    }

    let mut write_idx = 0;
    for i in 0..hits.len() {
        if keep[i] {
            hits[write_idx] = std::mem::take(&mut hits[i]);
            write_idx += 1;
        }
    }
    hits.truncate(write_idx);
}

/// Remove duplicate ray hits.
fn deduplicate_ray_hits(hits: &mut Vec<RayHit>) {
    if hits.len() <= 1 {
        return;
    }

    let mut keep = vec![true; hits.len()];

    for i in 0..hits.len() {
        if !keep[i] {
            continue;
        }
        for j in (i + 1)..hits.len() {
            if !keep[j] {
                continue;
            }
            let dist = (hits[i].point - hits[j].point).length();
            if dist < TOLERANCE_ABS * 100.0 {
                keep[j] = false;
            }
        }
    }

    let mut write_idx = 0;
    for i in 0..hits.len() {
        if keep[i] {
            hits[write_idx] = std::mem::take(&mut hits[i]);
            write_idx += 1;
        }
    }
    hits.truncate(write_idx);
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Plane, SphericalSurface, CylindricalSurface, Line3, Circle3};
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn line_through_box() {
        let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        // Box goes from (0,0,0) to (2,2,2), so center is at (1,1,1)
        // Line through the center of the box, along Z
        let intersections = intersect_line_with_brep(
            DVec3::new(1.0, 1.0, 5.0),  // Center of box is at (1,1,1)
            DVec3::NEG_Z,
            &box_brep,
            TOLERANCE_ABS,
        );

        // Should hit back face (z=2) and front face (z=0)
        assert_eq!(intersections.len(), 2);

        // Check that the points are on opposite ends
        assert!((intersections[0].point.z - 2.0).abs() < 0.1);
        assert!((intersections[1].point.z).abs() < 0.1);
    }

    #[test]
    fn line_through_sphere() {
        let sphere_brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        // Line through the center of the sphere
        let intersections = intersect_line_with_brep(
            DVec3::new(-5.0, 0.0, 0.0),
            DVec3::X,
            &sphere_brep,
            TOLERANCE_ABS,
        );

        // Should hit the sphere at two points
        assert_eq!(intersections.len(), 2);

        // Check distances from origin (sphere center)
        for i in &intersections {
            let dist = i.point.length();
            assert!((dist - 1.0).abs() < 0.1);
        }
    }

    #[test]
    fn line_misses_box() {
        let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        // Line outside the box
        let intersections = intersect_line_with_brep(
            DVec3::new(-5.0, 5.0, 0.0),
            DVec3::X,
            &box_brep,
            TOLERANCE_ABS,
        );

        assert!(intersections.is_empty());
    }

    #[test]
    fn ray_cast_through_box() {
        let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        let hits = ray_cast(DVec3::new(0.0, 0.0, 5.0), DVec3::NEG_Z, &box_brep);

        assert!(!hits.is_empty());

        // All distances should be positive
        for hit in &hits {
            assert!(hit.distance > 0.0);
        }
    }

    #[test]
    fn shoot_ray_limited_distance() {
        let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        // Shoot a ray that won't reach the box
        let hits = shoot_ray(
            DVec3::new(0.0, 0.0, 100.0),
            DVec3::NEG_Z,
            &box_brep,
            1.0, // Max distance too short to reach box
        );

        assert!(hits.is_empty());
    }

    #[test]
    fn point_inside_box() {
        let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        // Point at origin should be inside
        assert!(is_point_inside_by_ray(DVec3::ZERO, &box_brep));
    }

    #[test]
    fn point_outside_box() {
        let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        // Point far away should be outside
        assert!(!is_point_inside_by_ray(DVec3::new(10.0, 0.0, 0.0), &box_brep));
    }

    #[test]
    fn point_inside_sphere() {
        let sphere_brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        assert!(is_point_inside_by_ray(DVec3::ZERO, &sphere_brep));
    }

    #[test]
    fn point_outside_sphere() {
        let sphere_brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        assert!(!is_point_inside_by_ray(DVec3::new(2.0, 0.0, 0.0), &sphere_brep));
    }

    #[test]
    fn curve_line_through_box() {
        let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        // Box goes from (0,0,0) to (2,2,2), so center is at (1,1,1)
        let line = Curve3::Line(Line3 {
            origin: DVec3::new(1.0, 1.0, 5.0),
            direction: DVec3::NEG_Z,
        });

        let intersections = intersect_curve_with_brep(&line, &box_brep, TOLERANCE_ABS);

        assert_eq!(intersections.len(), 2);
    }

    #[test]
    fn curve_circle_through_box() {
        let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        });

        // Circle in XY plane at z=0
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.5,
        });

        let intersections = intersect_curve_with_brep(&circle, &box_brep, TOLERANCE_ABS);

        // The circle should intersect multiple faces
        assert!(intersections.len() >= 4);
    }

    #[test]
    fn line_plane_intersection() {
        let line = Line3 {
            origin: DVec3::new(0.0, 0.0, 5.0),
            direction: DVec3::NEG_Z,
        };
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };

        let result = intersect_line_plane(&line, &plane);
        assert!(result.is_some());

        let intersection = result.unwrap();
        assert!((intersection.point.z).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn line_cylinder_intersection() {
        let line = Line3 {
            origin: DVec3::new(-5.0, 0.0, 0.0),
            direction: DVec3::X,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };

        let hits = intersect_line_cylinder_range(&line, [-10.0, 10.0], &cyl);
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn line_sphere_intersection() {
        let line = Line3 {
            origin: DVec3::new(-5.0, 0.0, 0.0),
            direction: DVec3::X,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };

        let hits = intersect_line_sphere_range(&line, [-10.0, 10.0], &sphere);
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn intersection_result_types() {
        // Test CurveBRepIntersection
        let cbri = CurveBRepIntersection {
            point: DVec3::ZERO,
            param: 1.0,
            face_index: 0,
            uv: DVec2::ZERO,
        };
        assert_eq!(cbri.param, 1.0);

        // Test CurveFaceIntersection
        let cfi = CurveFaceIntersection {
            point: DVec3::ZERO,
            param: 2.0,
            uv: DVec2::new(0.5, 0.5),
        };
        assert_eq!(cfi.param, 2.0);

        // Test RayHit
        let rh = RayHit {
            point: DVec3::new(1.0, 2.0, 3.0),
            distance: 3.74,
            face_index: 1,
            normal: DVec3::Z,
        };
        assert!((rh.distance - 3.74).abs() < 1e-6);
    }

    #[test]
    fn multiple_rays_for_inside_test() {
        let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        // Box goes from (0,0,0) to (2,2,2)
        // Test several points
        let test_points = [
            (DVec3::new(1.0, 1.0, 1.0), true),   // Center - definitely inside
            (DVec3::new(0.5, 0.5, 0.5), true),   // Inside
            (DVec3::new(3.0, 1.0, 1.0), false),  // Outside (x=3 > 2)
            (DVec3::new(1.0, 3.0, 1.0), false),  // Outside (y=3 > 2)
        ];

        for (point, expected_inside) in test_points {
            let is_inside = is_point_inside_by_ray(point, &box_brep);
            assert_eq!(is_inside, expected_inside, "Point {:?} inside test failed", point);
        }
    }

    #[test]
    fn cylinder_intersection() {
        let cyl_brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        // Cylinder is along Y axis, centered at origin, with height 2 (y: -1 to 1)
        // Line through the center of the cylinder (z=0, not z=1)
        let intersections = intersect_line_with_brep(
            DVec3::new(-5.0, 0.0, 0.0),  // At z=0, goes through cylinder center
            DVec3::X,
            &cyl_brep,
            TOLERANCE_ABS,
        );

        // Should hit the curved surface at two points (x=-1 and x=1)
        assert!(intersections.len() >= 2);
    }
}
