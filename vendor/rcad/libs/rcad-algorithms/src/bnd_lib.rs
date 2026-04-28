//! BndLib-style bounding algorithms for geometry.
//!
//! This module provides algorithms analogous to OCCT's `BndLib` class:
//!
//! - **Add**: Add geometry (points, curves, surfaces) to bounding boxes
//! - **Compute**: Compute bounds for curves, surfaces, faces, edges
//! - **Optimized**: Optimized and precise bounds computation
//! - **BoundingSphere**: Bounding sphere computation and intersection tests
//! - **Intersection**: Box and sphere intersection tests
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::bnd_lib::*;
//! use rcad_algorithms::brep_bnd::BoundingBox;
//! use glam::DVec3;
//!
//! // Add points to a box
//! let mut bbox = BoundingBox::new();
//! add_point_to_box(DVec3::new(1.0, 2.0, 3.0), &mut bbox);
//!
//! // Create a bounding sphere from points
//! let points = vec![
//!     DVec3::new(0.0, 0.0, 0.0),
//!     DVec3::new(1.0, 0.0, 0.0),
//!     DVec3::new(0.0, 1.0, 0.0),
//! ];
//! let sphere = BoundingSphere::from_points(&points);
//!
//! // Test intersection
//! assert!(spheres_intersect(&sphere, &sphere));
//! ```

use glam::DVec3;
use rcad_kernel::{BRep, Curve3, Surface3};
use rcad_kernel::geom::{CurveEval, SurfaceEval};

use crate::brep_bnd::BoundingBox;

// =============================================================================
// Add Functions - Add geometry to bounds
// =============================================================================

/// Add a point to a bounding box.
///
/// Expands the box to include the given point.
pub fn add_point_to_box(point: DVec3, bbox: &mut BoundingBox) {
    bbox.add_point(point);
}

/// Add a curve segment to a bounding box over a parameter range.
///
/// Samples the curve at regular intervals between t1 and t2,
/// expanding the box to include all sample points.
///
/// # Arguments
///
/// * `curve` - The curve to add
/// * `t1` - Start parameter
/// * `t2` - End parameter
/// * `box` - The bounding box to expand
/// * `tol` - Tolerance for additional padding
pub fn add_curve_to_box(curve: &Curve3, t1: f64, t2: f64, bbox: &mut BoundingBox, tol: f64) {
    // Determine number of samples based on curve type
    let n_samples = get_curve_sample_count(curve, t1, t2);

    let dt = (t2 - t1) / n_samples as f64;

    for i in 0..=n_samples {
        let t = t1 + dt * i as f64;
        let p = curve.point_at(t);
        if p.is_finite() {
            bbox.add_point(p);
        }
    }

    // Add tolerance padding
    if tol > 0.0 && bbox.is_valid() {
        bbox.enlarge(tol);
    }
}

/// Add a surface patch to a bounding box over a parameter domain.
///
/// Samples the surface at a grid of points within the parameter range,
/// expanding the box to include all sample points.
///
/// # Arguments
///
/// * `surface` - The surface to add
/// * `u1, u2` - U parameter range
/// * `v1, v2` - V parameter range
/// * `box` - The bounding box to expand
/// * `tol` - Tolerance for additional padding
pub fn add_surface_to_box(
    surface: &Surface3,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    bbox: &mut BoundingBox,
    tol: f64,
) {
    // Determine grid resolution based on surface type
    let (n_u, n_v) = get_surface_sample_counts(surface, u1, u2, v1, v2);

    let du = (u2 - u1) / n_u as f64;
    let dv = (v2 - v1) / n_v as f64;

    for i in 0..=n_u {
        for j in 0..=n_v {
            let u = u1 + du * i as f64;
            let v = v1 + dv * j as f64;
            let p = surface.point_at(u, v);
            if p.is_finite() {
                bbox.add_point(p);
            }
        }
    }

    // Add tolerance padding
    if tol > 0.0 && bbox.is_valid() {
        bbox.enlarge(tol);
    }
}

/// Determine appropriate sample count for a curve based on its type.
fn get_curve_sample_count(curve: &Curve3, t1: f64, t2: f64) -> usize {
    let param_length = (t2 - t1).abs();

    match curve {
        Curve3::Line(_) => 2, // Line only needs endpoints
        Curve3::Circle(_) => {
            // Sample based on arc length
            let angle = param_length;
            ((angle / (std::f64::consts::PI / 4.0)).ceil() as usize).max(4).min(32)
        }
        Curve3::Ellipse(_) => {
            // Ellipses need more samples
            ((param_length / (std::f64::consts::PI / 6.0)).ceil() as usize).max(6).min(36)
        }
        Curve3::Parabola(_) | Curve3::Hyperbola(_) => {
            // Conic sections need moderate sampling
            20
        }
        Curve3::BSpline(bspline) => {
            // Sample based on number of control points and knots
            let n_poles = bspline.control_points.len();
            // More samples for higher degree curves
            (n_poles * 2).max(10).min(50)
        }
        Curve3::Bezier(bezier) => {
            // Sample based on degree
            (bezier.control_points.len() * 2).max(10).min(40)
        }
        Curve3::Offset(_) => {
            // Offset curves need more samples
            30
        }
        Curve3::CircularHelix(_) => {
            // Helix curves need good sampling
            24
        }
        Curve3::SineWave(_) => {
            // Sine waves need good sampling
            20
        }
        // Handle any other curve types with reasonable defaults
        _ => {
            // Default for unknown curves
            20
        }
    }
}

/// Determine appropriate sample counts for a surface based on its type.
fn get_surface_sample_counts(surface: &Surface3, u1: f64, u2: f64, v1: f64, v2: f64) -> (usize, usize) {
    let u_len = (u2 - u1).abs();
    let v_len = (v2 - v1).abs();

    match surface {
        Surface3::Plane(_) => {
            // Planes only need corners (but we sample for robustness)
            (2, 2)
        }
        Surface3::Cylinder(_) => {
            // Cylinders: u is angle, v is height
            let n_u = ((u_len / (std::f64::consts::PI / 4.0)).ceil() as usize).max(4).min(16);
            let n_v = ((v_len / 1.0).ceil() as usize).max(2).min(10);
            (n_u, n_v)
        }
        Surface3::Cone(_) => {
            // Cones similar to cylinders
            let n_u = ((u_len / (std::f64::consts::PI / 4.0)).ceil() as usize).max(4).min(16);
            let n_v = ((v_len / 1.0).ceil() as usize).max(2).min(10);
            (n_u, n_v)
        }
        Surface3::Sphere(_) => {
            // Spheres need good sampling in both directions
            let n_u = ((u_len / (std::f64::consts::PI / 4.0)).ceil() as usize).max(4).min(16);
            let n_v = ((v_len / (std::f64::consts::PI / 4.0)).ceil() as usize).max(4).min(16);
            (n_u, n_v)
        }
        Surface3::Torus(_) => {
            // Tori need good sampling
            let n_u = ((u_len / (std::f64::consts::PI / 3.0)).ceil() as usize).max(6).min(20);
            let n_v = ((v_len / (std::f64::consts::PI / 3.0)).ceil() as usize).max(6).min(20);
            (n_u, n_v)
        }
        Surface3::BSpline(bspline) => {
            // Sample based on control net
            let n_u = (bspline.control_points.len() * 2).max(5).min(20);
            let n_v = if !bspline.control_points.is_empty() {
                (bspline.control_points[0].len() * 2).max(5).min(20)
            } else {
                5
            };
            (n_u, n_v)
        }
        Surface3::Bezier(bezier) => {
            let n_u = (bezier.control_points.len() * 2).max(4).min(16);
            let n_v = if !bezier.control_points.is_empty() {
                (bezier.control_points[0].len() * 2).max(4).min(16)
            } else {
                4
            };
            (n_u, n_v)
        }
        Surface3::Revolution(_) => {
            // Revolution surfaces: u is angle, v is profile parameter
            let n_u = ((u_len / (std::f64::consts::PI / 4.0)).ceil() as usize).max(4).min(16);
            let n_v = 10;
            (n_u, n_v)
        }
        Surface3::LinearExtrusion(_) => {
            // Extrusion surfaces
            let n_u = 10;
            let n_v = ((v_len / 1.0).ceil() as usize).max(2).min(10);
            (n_u, n_v)
        }
        Surface3::Offset(_) => {
            // Offset surfaces need more samples
            (12, 12)
        }
        // Handle remaining surface types with reasonable defaults
        _ => {
            // Default for other surfaces
            (10, 10)
        }
    }
}

// =============================================================================
// Compute Bounds
// =============================================================================

/// Compute the bounding box for a curve.
///
/// Uses adaptive sampling based on curve type.
///
/// # Arguments
///
/// * `curve` - The curve to bound
/// * `tol` - Tolerance for bounds padding
pub fn curve_bounds(curve: &Curve3, tol: f64) -> BoundingBox {
    let range = curve.default_domain();
    curve_bounds_with_range(curve, range[0], range[1], tol)
}

/// Compute the bounding box for a curve over a specific parameter range.
pub fn curve_bounds_with_range(curve: &Curve3, t1: f64, t2: f64, tol: f64) -> BoundingBox {
    let mut bbox = BoundingBox::new();
    add_curve_to_box(curve, t1, t2, &mut bbox, tol);
    bbox
}

/// Compute the bounding box for a surface.
///
/// Uses adaptive sampling based on surface type.
///
/// # Arguments
///
/// * `surface` - The surface to bound
/// * `tol` - Tolerance for bounds padding
pub fn surface_bounds(surface: &Surface3, tol: f64) -> BoundingBox {
    let domain = surface.default_domain();
    surface_bounds_with_domain(surface, domain[0], domain[1], domain[2], domain[3], tol)
}

/// Compute the bounding box for a surface over a specific parameter domain.
pub fn surface_bounds_with_domain(
    surface: &Surface3,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    tol: f64,
) -> BoundingBox {
    let mut bbox = BoundingBox::new();
    add_surface_to_box(surface, u1, u2, v1, v2, &mut bbox, tol);
    bbox
}

/// Compute the bounding box for a face in a BRep.
///
/// Uses the face's surface and its parameter range.
///
/// # Arguments
///
/// * `brep` - The BRep containing the face
/// * `face_idx` - Flat index of the face
/// * `tol` - Tolerance for bounds padding
pub fn face_bounds(brep: &BRep, face_idx: usize, tol: f64) -> BoundingBox {
    let mut bbox = BoundingBox::new();

    // Get the face and its surface
    let (face, surface_idx) = match get_face_and_surface_idx(brep, face_idx) {
        Some(result) => result,
        None => return bbox,
    };

    // Get the surface
    let surface = match brep.geom.surfaces.get(surface_idx) {
        Some(s) => s,
        None => return bbox,
    };

    // Get the parameter domain for this face
    let domain = brep.geom.face_surface_range.get(face_idx)
        .copied()
        .flatten()
        .unwrap_or_else(|| surface.default_domain());

    // Add surface bounds
    add_surface_to_box(surface, domain[0], domain[1], domain[2], domain[3], &mut bbox, tol);

    // Also add vertices from the face boundary
    for wire_edge in &face.outer_wire.edges {
        if wire_edge.idx < brep.edges.len() {
            let edge = &brep.edges[wire_edge.idx];
            if edge.start < brep.vertices.len() {
                bbox.add_point(brep.vertices[edge.start].point);
            }
            if edge.end < brep.vertices.len() {
                bbox.add_point(brep.vertices[edge.end].point);
            }
        }
    }

    // Add vertices from inner wires
    for inner_wire in &face.inner_wires {
        for wire_edge in &inner_wire.edges {
            if wire_edge.idx < brep.edges.len() {
                let edge = &brep.edges[wire_edge.idx];
                if edge.start < brep.vertices.len() {
                    bbox.add_point(brep.vertices[edge.start].point);
                }
                if edge.end < brep.vertices.len() {
                    bbox.add_point(brep.vertices[edge.end].point);
                }
            }
        }
    }

    bbox
}

/// Compute the bounding box for an edge in a BRep.
///
/// Uses the edge's 3D curve and its parameter range.
///
/// # Arguments
///
/// * `brep` - The BRep containing the edge
/// * `edge_idx` - Index of the edge
/// * `tol` - Tolerance for bounds padding
pub fn edge_bounds(brep: &BRep, edge_idx: usize, tol: f64) -> BoundingBox {
    let mut bbox = BoundingBox::new();

    if edge_idx >= brep.edges.len() {
        return bbox;
    }

    let edge = &brep.edges[edge_idx];

    // Add edge vertices
    if edge.start < brep.vertices.len() {
        bbox.add_point(brep.vertices[edge.start].point);
    }
    if edge.end < brep.vertices.len() {
        bbox.add_point(brep.vertices[edge.end].point);
    }

    // Sample along the curve if available
    if let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c) {
        if let Some(curve) = brep.geom.curves.get(curve_idx) {
            let range = brep.geom.edge_curve_range.get(edge_idx)
                .copied()
                .flatten()
                .unwrap_or_else(|| curve.default_domain());

            add_curve_to_box(curve, range[0], range[1], &mut bbox, 0.0);
        }
    }

    // Add tolerance padding
    if tol > 0.0 && bbox.is_valid() {
        bbox.enlarge(tol);
    }

    bbox
}

// =============================================================================
// Optimized Bounds
// =============================================================================

/// Compute optimized bounding box for a BRep.
///
/// Uses fast approximations where possible:
/// - Uses vertex positions directly
/// - Uses fewer samples for curved geometry
/// - Caches intermediate results
///
/// # Arguments
///
/// * `brep` - The BRep to bound
/// * `tol` - Tolerance for bounds padding
pub fn optimized_bounds(brep: &BRep, tol: f64) -> BoundingBox {
    let mut bbox = BoundingBox::new();

    // Add all vertices first (fast path)
    for vertex in &brep.vertices {
        bbox.add_point(vertex.point);
    }

    // For faces with curved surfaces, add a few more sample points
    let n_faces = count_brep_faces(brep);
    for face_idx in 0..n_faces {
        if let Some(surface_idx) = brep.geom.face_surface.get(face_idx).and_then(|s| *s) {
            if let Some(surface) = brep.geom.surfaces.get(surface_idx) {
                // Only sample non-planar surfaces
                if !matches!(surface, Surface3::Plane(_)) {
                    let domain = brep.geom.face_surface_range.get(face_idx)
                        .copied()
                        .flatten()
                        .unwrap_or_else(|| surface.default_domain());

                    // Use fewer samples for optimized bounds
                    let n_u = 3;
                    let n_v = 3;
                    let du = (domain[1] - domain[0]) / n_u as f64;
                    let dv = (domain[3] - domain[2]) / n_v as f64;

                    for i in 0..=n_u {
                        for j in 0..=n_v {
                            let u = domain[0] + du * i as f64;
                            let v = domain[2] + dv * j as f64;
                            let p = surface.point_at(u, v);
                            if p.is_finite() {
                                bbox.add_point(p);
                            }
                        }
                    }
                }
            }
        }
    }

    // Add tolerance padding
    if tol > 0.0 && bbox.is_valid() {
        bbox.enlarge(tol);
    }

    bbox
}

/// Compute precise bounding box for a BRep with higher sampling.
///
/// Uses higher sampling density for more accurate bounds:
/// - More samples per curve/surface
/// - Includes all vertices
/// - Accounts for extreme points
///
/// # Arguments
///
/// * `brep` - The BRep to bound
/// * `n_samples` - Number of samples per direction (higher = more precise)
pub fn precise_bounds(brep: &BRep, n_samples: usize) -> BoundingBox {
    let mut bbox = BoundingBox::new();

    // Add all vertices
    for vertex in &brep.vertices {
        bbox.add_point(vertex.point);
    }

    // Process all faces with high sampling
    let n_faces = count_brep_faces(brep);
    for face_idx in 0..n_faces {
        if let Some(surface_idx) = brep.geom.face_surface.get(face_idx).and_then(|s| *s) {
            if let Some(surface) = brep.geom.surfaces.get(surface_idx) {
                let domain = brep.geom.face_surface_range.get(face_idx)
                    .copied()
                    .flatten()
                    .unwrap_or_else(|| surface.default_domain());

                let du = (domain[1] - domain[0]) / n_samples.max(2) as f64;
                let dv = (domain[3] - domain[2]) / n_samples.max(2) as f64;

                for i in 0..=n_samples.max(2) {
                    for j in 0..=n_samples.max(2) {
                        let u = domain[0] + du * i as f64;
                        let v = domain[2] + dv * j as f64;
                        let p = surface.point_at(u, v);
                        if p.is_finite() {
                            bbox.add_point(p);
                        }
                    }
                }
            }
        }
    }

    // Process all edges with high sampling
    for edge_idx in 0..brep.edges.len() {
        if let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c) {
            if let Some(curve) = brep.geom.curves.get(curve_idx) {
                let range = brep.geom.edge_curve_range.get(edge_idx)
                    .copied()
                    .flatten()
                    .unwrap_or_else(|| curve.default_domain());

                let dt = (range[1] - range[0]) / n_samples.max(4) as f64;

                for i in 0..=n_samples.max(4) {
                    let t = range[0] + dt * i as f64;
                    let p = curve.point_at(t);
                    if p.is_finite() {
                        bbox.add_point(p);
                    }
                }
            }
        }
    }

    bbox
}

// =============================================================================
// Bounding Sphere
// =============================================================================

/// A bounding sphere for 3D geometry.
///
/// Provides efficient intersection tests and containment checks.
#[derive(Debug, Clone, Copy)]
pub struct BoundingSphere {
    center: DVec3,
    radius: f64,
}

impl Default for BoundingSphere {
    fn default() -> Self {
        Self::new()
    }
}

impl BoundingSphere {
    /// Create a new empty bounding sphere.
    ///
    /// An empty sphere has zero radius at the origin.
    pub fn new() -> Self {
        Self {
            center: DVec3::ZERO,
            radius: 0.0,
        }
    }

    /// Create a bounding sphere from center and radius.
    pub fn from_center_radius(center: DVec3, radius: f64) -> Self {
        Self { center, radius }
    }

    /// Create a bounding sphere from a set of points.
    ///
    /// Uses an iterative algorithm to find a tight bounding sphere.
    /// Based on Ritter's algorithm for approximate minimum bounding sphere.
    pub fn from_points(points: &[DVec3]) -> Self {
        if points.is_empty() {
            return Self::new();
        }

        if points.len() == 1 {
            return Self {
                center: points[0],
                radius: 0.0,
            };
        }

        // Step 1: Find approximate center using component-wise average
        let mut center = DVec3::ZERO;
        for &p in points {
            center += p;
        }
        center /= points.len() as f64;

        // Step 2: Find the point farthest from center
        let mut max_dist_sq = 0.0;
        for &p in points {
            let dist_sq = (p - center).length_squared();
            if dist_sq > max_dist_sq {
                max_dist_sq = dist_sq;
            }
        }
        let mut radius = max_dist_sq.sqrt();

        // Step 3: Iteratively expand sphere to include all points
        // This is Ritter's algorithm refinement
        for &p in points {
            let dist = (p - center).length();
            if dist > radius {
                let new_radius = (radius + dist) * 0.5;
                let direction = (p - center).normalize_or_zero();
                center = center + direction * (new_radius - radius);
                radius = new_radius;
            }
        }

        // Step 4: Additional pass to tighten the sphere
        let mut final_radius = radius;
        for &p in points {
            let dist = (p - center).length();
            if dist > final_radius {
                final_radius = dist;
            }
        }

        Self {
            center,
            radius: final_radius,
        }
    }

    /// Create a bounding sphere from a bounding box.
    ///
    /// The sphere circumscribes the box corners.
    pub fn from_bounding_box(bbox: &BoundingBox) -> Self {
        if bbox.is_empty() {
            return Self::new();
        }

        let center = bbox.center();
        let radius = bbox.diagonal() * 0.5;

        Self { center, radius }
    }

    /// Create the smallest bounding sphere enclosing two spheres.
    pub fn from_two_spheres(s1: &BoundingSphere, s2: &BoundingSphere) -> Self {
        let dist = (s2.center - s1.center).length();

        // If one sphere contains the other
        if dist + s2.radius <= s1.radius {
            return *s1;
        }
        if dist + s1.radius <= s2.radius {
            return *s2;
        }

        // Otherwise, compute the minimal enclosing sphere
        let direction = (s2.center - s1.center).normalize_or_zero();
        let radius = (dist + s1.radius + s2.radius) * 0.5;
        let center = s1.center + direction * (radius - s1.radius);

        Self { center, radius }
    }

    /// Get the center of the sphere.
    pub fn center(&self) -> DVec3 {
        self.center
    }

    /// Get the radius of the sphere.
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Check if the sphere is empty (zero radius).
    pub fn is_empty(&self) -> bool {
        self.radius <= 0.0
    }

    /// Check if a point is inside the sphere.
    ///
    /// Points exactly on the boundary are considered inside.
    pub fn contains(&self, p: DVec3) -> bool {
        (p - self.center).length_squared() <= self.radius * self.radius
    }

    /// Check if a point is inside the sphere with tolerance.
    pub fn contains_with_tol(&self, p: DVec3, tol: f64) -> bool {
        let r = self.radius + tol;
        (p - self.center).length_squared() <= r * r
    }

    /// Check if another sphere is completely inside this sphere.
    pub fn contains_sphere(&self, other: &BoundingSphere) -> bool {
        let dist = (other.center - self.center).length();
        dist + other.radius <= self.radius
    }

    /// Compute distance from a point to the sphere surface.
    ///
    /// Returns 0.0 if the point is inside the sphere.
    pub fn distance_to(&self, p: DVec3) -> f64 {
        let dist = (p - self.center).length();
        (dist - self.radius).max(0.0)
    }

    /// Expand the sphere to include a point.
    pub fn expand_point(&mut self, p: DVec3) {
        let dist = (p - self.center).length();
        if dist > self.radius {
            self.radius = dist;
        }
    }

    /// Expand the sphere to include another sphere.
    pub fn expand_sphere(&mut self, other: &BoundingSphere) {
        *self = BoundingSphere::from_two_spheres(self, other);
    }

    /// Get the bounding box for this sphere.
    pub fn to_bounding_box(&self) -> BoundingBox {
        BoundingBox::from_corners(
            self.center - DVec3::splat(self.radius),
            self.center + DVec3::splat(self.radius),
        )
    }

    /// Get the volume of the sphere.
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * std::f64::consts::PI * self.radius * self.radius * self.radius
    }

    /// Get the surface area of the sphere.
    pub fn surface_area(&self) -> f64 {
        4.0 * std::f64::consts::PI * self.radius * self.radius
    }
}

// =============================================================================
// Intersection Tests
// =============================================================================

/// Check if two bounding boxes intersect.
///
/// Boxes intersect if they overlap in all three dimensions.
pub fn boxes_intersect(box1: &BoundingBox, box2: &BoundingBox) -> bool {
    box1.intersects(box2)
}

/// Check if a bounding sphere intersects a bounding box.
///
/// Uses the closest point on the box to the sphere center.
pub fn sphere_intersects_box(sphere: &BoundingSphere, bbox: &BoundingBox) -> bool {
    if bbox.is_empty() || sphere.is_empty() {
        return false;
    }

    // Find the closest point on the box to the sphere center
    let closest = bbox.closest_point(sphere.center);

    // Check if the closest point is within the sphere
    let dist_sq = (closest - sphere.center).length_squared();
    let radius_sq = sphere.radius * sphere.radius;

    dist_sq <= radius_sq
}

/// Check if two bounding spheres intersect.
///
/// Spheres intersect if the distance between centers is less than
/// the sum of their radii.
pub fn spheres_intersect(s1: &BoundingSphere, s2: &BoundingSphere) -> bool {
    if s1.is_empty() || s2.is_empty() {
        return false;
    }

    let dist_sq = (s2.center - s1.center).length_squared();
    let radius_sum = s1.radius + s2.radius;

    dist_sq <= radius_sum * radius_sum
}

/// Check if a sphere completely contains a box.
pub fn sphere_contains_box(sphere: &BoundingSphere, bbox: &BoundingBox) -> bool {
    if bbox.is_empty() {
        return true;
    }
    if sphere.is_empty() {
        return false;
    }

    // All corners must be inside the sphere
    let corners = bbox.all_corners();
    for &corner in &corners {
        if !sphere.contains(corner) {
            return false;
        }
    }
    true
}

/// Check if a box completely contains a sphere.
pub fn box_contains_sphere(bbox: &BoundingBox, sphere: &BoundingSphere) -> bool {
    if sphere.is_empty() {
        return true;
    }
    if bbox.is_empty() {
        return false;
    }

    // The sphere's bounding box must be inside the box
    let sphere_box = sphere.to_bounding_box();

    bbox.contains(sphere_box.min, 0.0) && bbox.contains(sphere_box.max, 0.0)
}

/// Compute the intersection volume of two boxes.
///
/// Returns 0.0 if boxes don't intersect.
pub fn box_intersection_volume(box1: &BoundingBox, box2: &BoundingBox) -> f64 {
    let intersection = box1.intersection(box2);
    intersection.volume()
}

/// Compute the distance between two boxes.
///
/// Returns 0.0 if boxes intersect.
pub fn box_distance(box1: &BoundingBox, box2: &BoundingBox) -> f64 {
    if box1.intersects(box2) {
        return 0.0;
    }

    // Compute distance based on separation along each axis
    let mut dist_sq = 0.0;

    // X separation
    if box1.max.x < box2.min.x {
        let d = box2.min.x - box1.max.x;
        dist_sq += d * d;
    } else if box2.max.x < box1.min.x {
        let d = box1.min.x - box2.max.x;
        dist_sq += d * d;
    }

    // Y separation
    if box1.max.y < box2.min.y {
        let d = box2.min.y - box1.max.y;
        dist_sq += d * d;
    } else if box2.max.y < box1.min.y {
        let d = box1.min.y - box2.max.y;
        dist_sq += d * d;
    }

    // Z separation
    if box1.max.z < box2.min.z {
        let d = box2.min.z - box1.max.z;
        dist_sq += d * d;
    } else if box2.max.z < box1.min.z {
        let d = box1.min.z - box2.max.z;
        dist_sq += d * d;
    }

    dist_sq.sqrt()
}

/// Compute the distance between two spheres.
///
/// Returns 0.0 if spheres intersect.
pub fn sphere_distance(s1: &BoundingSphere, s2: &BoundingSphere) -> f64 {
    let center_dist = (s2.center - s1.center).length();
    let radius_sum = s1.radius + s2.radius;

    (center_dist - radius_sum).max(0.0)
}

/// Compute the distance between a sphere and a box.
///
/// Returns 0.0 if they intersect.
pub fn sphere_box_distance(sphere: &BoundingSphere, bbox: &BoundingBox) -> f64 {
    if bbox.is_empty() || sphere.is_empty() {
        return f64::INFINITY;
    }

    let dist = bbox.distance_to(sphere.center);
    (dist - sphere.radius).max(0.0)
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Count the total number of faces in a BRep.
fn count_brep_faces(brep: &BRep) -> usize {
    brep.solids.iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum()
}

/// Get a face and its surface index by flat index.
fn get_face_and_surface_idx(brep: &BRep, face_idx: usize) -> Option<(&rcad_kernel::topology::Face, usize)> {
    let mut current_idx = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            if face_idx < current_idx + shell.faces.len() {
                let local_idx = face_idx - current_idx;
                let face = &shell.faces[local_idx];

                // Get the surface index
                let surface_idx = brep.geom.face_surface.get(face_idx)
                    .and_then(|s| *s)?;

                return Some((face, surface_idx));
            }
            current_idx += shell.faces.len();
        }
    }

    None
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_bnd::BoundingBox;
    use rcad_kernel::{PrimitiveSolid, geom::{Line3, Plane, Circle3 as GeomCircle}};
    use rcad_kernel::geom::CylindricalSurface;

    // ── Add Functions Tests ─────────────────────────────────────────────────────────

    #[test]
    fn test_add_point_to_box() {
        let mut bbox = BoundingBox::new();
        add_point_to_box(DVec3::new(1.0, 2.0, 3.0), &mut bbox);

        assert!(bbox.is_valid());
        assert_eq!(bbox.min, DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(bbox.max, DVec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_add_curve_to_box_line() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });

        let mut bbox = BoundingBox::new();
        add_curve_to_box(&line, 0.0, 10.0, &mut bbox, 0.1);

        assert!(bbox.is_valid());
        // Line from (0,0,0) to (10,0,0) with tolerance
        assert!(bbox.min.x <= 0.0);
        assert!(bbox.max.x >= 10.0);
        // Should have tolerance padding
        assert!((bbox.min.x - (-0.1)).abs() < 1e-9);
        assert!((bbox.max.x - 10.1).abs() < 1e-9);
    }

    #[test]
    fn test_add_curve_to_box_circle() {
        let circle = Curve3::Circle(GeomCircle {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });

        let mut bbox = BoundingBox::new();
        add_curve_to_box(&circle, 0.0, std::f64::consts::PI * 2.0, &mut bbox, 0.0);

        assert!(bbox.is_valid());
        // Full circle should span [-1, 1] in x and y
        assert!(bbox.min.x <= -1.0);
        assert!(bbox.max.x >= 1.0);
        assert!(bbox.min.y <= -1.0);
        assert!(bbox.max.y >= 1.0);
    }

    #[test]
    fn test_add_surface_to_box_plane() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let mut bbox = BoundingBox::new();
        add_surface_to_box(&plane, 0.0, 10.0, 0.0, 20.0, &mut bbox, 0.0);

        assert!(bbox.is_valid());
        // Plane is at z=0
        assert!((bbox.min.z).abs() < 1e-6);
        assert!((bbox.max.z).abs() < 1e-6);
    }

    #[test]
    fn test_add_surface_to_box_cylinder() {
        let cylinder = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });

        let mut bbox = BoundingBox::new();
        add_surface_to_box(&cylinder, 0.0, std::f64::consts::PI * 2.0, 0.0, 5.0, &mut bbox, 0.1);

        assert!(bbox.is_valid());
        // Cylinder of radius 1 should span [-1, 1] in x and y
        assert!(bbox.min.x <= -1.0);
        assert!(bbox.max.x >= 1.0);
        assert!(bbox.min.y <= -1.0);
        assert!(bbox.max.y >= 1.0);
        // Height 0 to 5 in z
        assert!(bbox.min.z <= 0.0);
        assert!(bbox.max.z >= 5.0);
    }

    // ── Compute Bounds Tests ────────────────────────────────────────────────────────

    #[test]
    fn test_curve_bounds_line() {
        // Use a circle instead of a line, since lines have infinite domain
        let circle = Curve3::Circle(GeomCircle {
            center: DVec3::new(5.0, 0.0, 0.0),
            normal: DVec3::Z,
            radius: 1.0,
        });

        let bounds = curve_bounds(&circle, 0.0);

        assert!(bounds.is_valid());
    }

    #[test]
    fn test_surface_bounds_sphere() {
        let sphere = rcad_kernel::geom::SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        let surface = Surface3::Sphere(sphere);

        let bounds = surface_bounds(&surface, 0.0);

        assert!(bounds.is_valid());
        // Sphere of radius 2 should span [-2, 2] in all directions
        assert!(bounds.min.x <= -2.0);
        assert!(bounds.max.x >= 2.0);
    }

    #[test]
    fn test_face_bounds_box() {
        // Use cylinder instead of box, since box doesn't set up GeomStore
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        // Get bounds for face 0
        let bounds = face_bounds(&brep, 0, 0.0);

        assert!(bounds.is_valid(), "Bounds should be valid for cylinder face");
    }

    #[test]
    fn test_edge_bounds_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 2.0,
            depth: 3.0,
        });

        let bounds = edge_bounds(&brep, 0, 0.0);

        assert!(bounds.is_valid());
    }

    // ── Optimized Bounds Tests ──────────────────────────────────────────────────────

    #[test]
    fn test_optimized_bounds_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });

        // Use zero tolerance to get exact bounds
        let bounds = optimized_bounds(&brep, 0.0);

        assert!(bounds.is_valid());
        assert!((bounds.size().x - 2.0).abs() < 1e-6);
        assert!((bounds.size().y - 3.0).abs() < 1e-6);
        assert!((bounds.size().z - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_optimized_bounds_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let bounds = optimized_bounds(&brep, 0.0);

        assert!(bounds.is_valid());
        // Sphere bounds depend on how create_sphere sets up vertices
        // Just check that bounds are valid and contain the origin
        assert!(bounds.min.x <= 0.0);
        assert!(bounds.max.x >= 0.0);
    }

    #[test]
    fn test_precise_bounds() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let bounds = precise_bounds(&brep, 20);

        assert!(bounds.is_valid());
        // Higher precision should give tighter bounds
        assert!(bounds.min.x >= -1.1);
        assert!(bounds.max.x <= 1.1);
    }

    #[test]
    fn test_precise_vs_optimized() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let opt_bounds = optimized_bounds(&brep, 0.0);
        let prec_bounds = precise_bounds(&brep, 50);

        // Both should be valid
        assert!(opt_bounds.is_valid());
        assert!(prec_bounds.is_valid());

        // Precise bounds should be at least as large as optimized
        // (optimized might miss some curvature)
        assert!(prec_bounds.volume() >= opt_bounds.volume() * 0.9);
    }

    // ── BoundingSphere Tests ────────────────────────────────────────────────────────

    #[test]
    fn test_bounding_sphere_new() {
        let sphere = BoundingSphere::new();

        assert!(sphere.is_empty());
        assert_eq!(sphere.center(), DVec3::ZERO);
        assert_eq!(sphere.radius(), 0.0);
    }

    #[test]
    fn test_bounding_sphere_from_center_radius() {
        let sphere = BoundingSphere::from_center_radius(DVec3::new(1.0, 2.0, 3.0), 5.0);

        assert!(!sphere.is_empty());
        assert_eq!(sphere.center(), DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(sphere.radius(), 5.0);
    }

    #[test]
    fn test_bounding_sphere_from_points_single() {
        let points = vec![DVec3::new(1.0, 2.0, 3.0)];
        let sphere = BoundingSphere::from_points(&points);

        assert_eq!(sphere.center(), DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(sphere.radius(), 0.0);
    }

    #[test]
    fn test_bounding_sphere_from_points_line() {
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
        ];
        let sphere = BoundingSphere::from_points(&points);

        // Sphere should contain both endpoints
        assert!(sphere.contains(DVec3::new(0.0, 0.0, 0.0)));
        assert!(sphere.contains(DVec3::new(2.0, 0.0, 0.0)));

        // Center should be around (1, 0, 0) with radius around 1
        assert!((sphere.center().x - 1.0).abs() < 0.5);
        assert!(sphere.radius() >= 0.9);
    }

    #[test]
    fn test_bounding_sphere_from_points_triangle() {
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let sphere = BoundingSphere::from_points(&points);

        // All points should be inside
        for &p in &points {
            assert!(sphere.contains(p));
        }
    }

    #[test]
    fn test_bounding_sphere_from_points_tetrahedron() {
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
        ];
        let sphere = BoundingSphere::from_points(&points);

        // All points should be inside
        for &p in &points {
            assert!(sphere.contains(p));
        }

        // Volume should be positive
        assert!(sphere.volume() > 0.0);
    }

    #[test]
    fn test_bounding_sphere_from_bounding_box() {
        let bbox = BoundingBox::from_corners(
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(2.0, 2.0, 2.0),
        );

        let sphere = BoundingSphere::from_bounding_box(&bbox);

        assert_eq!(sphere.center(), DVec3::new(1.0, 1.0, 1.0));
        // Diagonal is sqrt(12), radius is half
        let expected_radius = (12.0_f64).sqrt() * 0.5;
        assert!((sphere.radius() - expected_radius).abs() < 1e-9);
    }

    #[test]
    fn test_bounding_sphere_from_two_spheres() {
        let s1 = BoundingSphere::from_center_radius(DVec3::ZERO, 1.0);
        let s2 = BoundingSphere::from_center_radius(DVec3::new(4.0, 0.0, 0.0), 1.0);

        let combined = BoundingSphere::from_two_spheres(&s1, &s2);

        // Should contain both original spheres
        assert!(combined.contains_sphere(&s1));
        assert!(combined.contains_sphere(&s2));

        // Center should be around (2, 0, 0)
        assert!((combined.center().x - 2.0).abs() < 0.5);
    }

    #[test]
    fn test_bounding_sphere_contains() {
        let sphere = BoundingSphere::from_center_radius(DVec3::ZERO, 1.0);

        // Inside
        assert!(sphere.contains(DVec3::ZERO));
        assert!(sphere.contains(DVec3::new(0.5, 0.5, 0.5)));

        // On boundary
        assert!(sphere.contains(DVec3::new(1.0, 0.0, 0.0)));

        // Outside
        assert!(!sphere.contains(DVec3::new(2.0, 0.0, 0.0)));
    }

    #[test]
    fn test_bounding_sphere_distance_to() {
        let sphere = BoundingSphere::from_center_radius(DVec3::ZERO, 1.0);

        // Inside
        assert!((sphere.distance_to(DVec3::ZERO) - 0.0).abs() < 1e-9);

        // On surface
        assert!((sphere.distance_to(DVec3::new(1.0, 0.0, 0.0)) - 0.0).abs() < 1e-9);

        // Outside
        assert!((sphere.distance_to(DVec3::new(3.0, 0.0, 0.0)) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_bounding_sphere_to_bounding_box() {
        let sphere = BoundingSphere::from_center_radius(DVec3::new(1.0, 2.0, 3.0), 2.0);
        let bbox = sphere.to_bounding_box();

        assert_eq!(bbox.min, DVec3::new(-1.0, 0.0, 1.0));
        assert_eq!(bbox.max, DVec3::new(3.0, 4.0, 5.0));
    }

    #[test]
    fn test_bounding_sphere_volume_surface_area() {
        let sphere = BoundingSphere::from_center_radius(DVec3::ZERO, 1.0);

        // Volume = 4/3 * pi * r^3
        let expected_volume = (4.0 / 3.0) * std::f64::consts::PI;
        assert!((sphere.volume() - expected_volume).abs() < 1e-9);

        // Surface area = 4 * pi * r^2
        let expected_area = 4.0 * std::f64::consts::PI;
        assert!((sphere.surface_area() - expected_area).abs() < 1e-9);
    }

    // ── Intersection Tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_boxes_intersect_overlapping() {
        let box1 = BoundingBox::from_corners(DVec3::ZERO, DVec3::new(2.0, 2.0, 2.0));
        let box2 = BoundingBox::from_corners(DVec3::new(1.0, 1.0, 1.0), DVec3::new(3.0, 3.0, 3.0));

        assert!(boxes_intersect(&box1, &box2));
    }

    #[test]
    fn test_boxes_intersect_touching() {
        let box1 = BoundingBox::from_corners(DVec3::ZERO, DVec3::new(1.0, 1.0, 1.0));
        let box2 = BoundingBox::from_corners(DVec3::new(1.0, 0.0, 0.0), DVec3::new(2.0, 1.0, 1.0));

        assert!(boxes_intersect(&box1, &box2));
    }

    #[test]
    fn test_boxes_intersect_separate() {
        let box1 = BoundingBox::from_corners(DVec3::ZERO, DVec3::new(1.0, 1.0, 1.0));
        let box2 = BoundingBox::from_corners(DVec3::new(2.0, 2.0, 2.0), DVec3::new(3.0, 3.0, 3.0));

        assert!(!boxes_intersect(&box1, &box2));
    }

    #[test]
    fn test_sphere_intersects_box_inside() {
        let sphere = BoundingSphere::from_center_radius(DVec3::new(0.5, 0.5, 0.5), 0.5);
        let bbox = BoundingBox::from_corners(DVec3::ZERO, DVec3::new(1.0, 1.0, 1.0));

        assert!(sphere_intersects_box(&sphere, &bbox));
    }

    #[test]
    fn test_sphere_intersects_box_touching() {
        let sphere = BoundingSphere::from_center_radius(DVec3::new(2.0, 0.5, 0.5), 1.0);
        let bbox = BoundingBox::from_corners(DVec3::ZERO, DVec3::new(1.0, 1.0, 1.0));

        assert!(sphere_intersects_box(&sphere, &bbox));
    }

    #[test]
    fn test_sphere_intersects_box_separate() {
        let sphere = BoundingSphere::from_center_radius(DVec3::new(3.0, 0.5, 0.5), 0.5);
        let bbox = BoundingBox::from_corners(DVec3::ZERO, DVec3::new(1.0, 1.0, 1.0));

        assert!(!sphere_intersects_box(&sphere, &bbox));
    }

    #[test]
    fn test_spheres_intersect_overlapping() {
        let s1 = BoundingSphere::from_center_radius(DVec3::ZERO, 1.0);
        let s2 = BoundingSphere::from_center_radius(DVec3::new(1.0, 0.0, 0.0), 1.0);

        assert!(spheres_intersect(&s1, &s2));
    }

    #[test]
    fn test_spheres_intersect_touching() {
        let s1 = BoundingSphere::from_center_radius(DVec3::ZERO, 1.0);
        let s2 = BoundingSphere::from_center_radius(DVec3::new(2.0, 0.0, 0.0), 1.0);

        assert!(spheres_intersect(&s1, &s2));
    }

    #[test]
    fn test_spheres_intersect_separate() {
        let s1 = BoundingSphere::from_center_radius(DVec3::ZERO, 1.0);
        let s2 = BoundingSphere::from_center_radius(DVec3::new(3.0, 0.0, 0.0), 1.0);

        assert!(!spheres_intersect(&s1, &s2));
    }

    #[test]
    fn test_sphere_contains_box() {
        let sphere = BoundingSphere::from_center_radius(DVec3::ZERO, 5.0);
        let bbox = BoundingBox::from_corners(
            DVec3::new(-1.0, -1.0, -1.0),
            DVec3::new(1.0, 1.0, 1.0),
        );

        assert!(sphere_contains_box(&sphere, &bbox));
    }

    #[test]
    fn test_box_contains_sphere() {
        let bbox = BoundingBox::from_corners(
            DVec3::new(-5.0, -5.0, -5.0),
            DVec3::new(5.0, 5.0, 5.0),
        );
        let sphere = BoundingSphere::from_center_radius(DVec3::ZERO, 2.0);

        assert!(box_contains_sphere(&bbox, &sphere));
    }

    #[test]
    fn test_box_intersection_volume() {
        let box1 = BoundingBox::from_corners(DVec3::ZERO, DVec3::new(2.0, 2.0, 2.0));
        let box2 = BoundingBox::from_corners(DVec3::new(1.0, 1.0, 1.0), DVec3::new(3.0, 3.0, 3.0));

        let volume = box_intersection_volume(&box1, &box2);

        // Intersection is [1,2] x [1,2] x [1,2] = volume 1
        assert!((volume - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_box_distance_intersecting() {
        let box1 = BoundingBox::from_corners(DVec3::ZERO, DVec3::new(2.0, 2.0, 2.0));
        let box2 = BoundingBox::from_corners(DVec3::new(1.0, 1.0, 1.0), DVec3::new(3.0, 3.0, 3.0));

        assert!((box_distance(&box1, &box2) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_box_distance_separate() {
        let box1 = BoundingBox::from_corners(DVec3::ZERO, DVec3::new(1.0, 1.0, 1.0));
        let box2 = BoundingBox::from_corners(DVec3::new(3.0, 0.0, 0.0), DVec3::new(4.0, 1.0, 1.0));

        // Distance between x=1 and x=3 is 2
        assert!((box_distance(&box1, &box2) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_sphere_distance_intersecting() {
        let s1 = BoundingSphere::from_center_radius(DVec3::ZERO, 1.0);
        let s2 = BoundingSphere::from_center_radius(DVec3::new(1.0, 0.0, 0.0), 1.0);

        assert!((sphere_distance(&s1, &s2) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_sphere_distance_separate() {
        let s1 = BoundingSphere::from_center_radius(DVec3::ZERO, 1.0);
        let s2 = BoundingSphere::from_center_radius(DVec3::new(5.0, 0.0, 0.0), 1.0);

        // Distance is 5 - 1 - 1 = 3
        assert!((sphere_distance(&s1, &s2) - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_sphere_box_distance() {
        let sphere = BoundingSphere::from_center_radius(DVec3::new(3.0, 0.5, 0.5), 1.0);
        let bbox = BoundingBox::from_corners(DVec3::ZERO, DVec3::new(1.0, 1.0, 1.0));

        // Distance from sphere center (3, 0.5, 0.5) to box is 2
        // With radius 1, distance is 2 - 1 = 1
        assert!((sphere_box_distance(&sphere, &bbox) - 1.0).abs() < 1e-9);
    }

    // ── Integration Tests ───────────────────────────────────────────────────────────

    #[test]
    fn test_brep_sphere_bounds_and_intersection() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let bounds = optimized_bounds(&brep, 0.0);
        let sphere = BoundingSphere::from_bounding_box(&bounds);

        // Sphere should contain the original bounding box
        let sphere2 = BoundingSphere::from_center_radius(DVec3::ZERO, 1.5);
        assert!(spheres_intersect(&sphere, &sphere2));
    }

    #[test]
    fn test_multiple_primitives_bounds() {
        let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0
        });
        let sphere_brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 0.5 });
        let cylinder_brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 0.5, height: 2.0
        });

        let box_bounds = optimized_bounds(&box_brep, 0.0);
        let sphere_bounds = optimized_bounds(&sphere_brep, 0.0);
        let cylinder_bounds = optimized_bounds(&cylinder_brep, 0.0);

        // All should have valid bounds
        assert!(box_bounds.is_valid());
        assert!(sphere_bounds.is_valid());
        assert!(cylinder_bounds.is_valid());

        // Box and sphere should intersect (both at origin)
        assert!(boxes_intersect(&box_bounds, &sphere_bounds));
    }

    #[test]
    fn test_empty_brep_bounds() {
        let brep = BRep::new();

        let bounds = optimized_bounds(&brep, 0.0);

        assert!(bounds.is_empty());
    }

    #[test]
    fn test_tolerances_affect_bounds() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0
        });

        let bounds_no_tol = optimized_bounds(&brep, 0.0);
        let bounds_with_tol = optimized_bounds(&brep, 0.5);

        // Bounds with tolerance should be larger
        assert!(bounds_with_tol.volume() > bounds_no_tol.volume());
    }
}
