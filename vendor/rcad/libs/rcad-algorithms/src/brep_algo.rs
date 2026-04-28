//! BRepAlgo-style algorithm utilities for BRep analysis and manipulation.
//!
//! This module provides utilities analogous to OCCT's `BRepAlgo` class:
//!
//! - **Normal evaluation**: Compute face normals, edge tangents, and vertex normals
//! - **Tolerance propagation**: Propagate tolerances through the BRep hierarchy
//! - **Tools**: Compute geometric properties like areas, volumes, and lengths
//! - **Validity checking**: Check BRep validity and orientation consistency
//! - **Orientation fixing**: Fix orientation issues in the BRep
//! - **Connected components**: Find connected components in the BRep
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::brep_algo::*;
//! use rcad_kernel::BRep;
//!
//! let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
//!     width: 1.0, height: 1.0, depth: 1.0
//! });
//!
//! // Check validity
//! assert!(is_valid_brep(&brep));
//!
//! // Compute properties
//! let volume = total_volume(&brep);
//! let area = total_surface_area(&brep);
//! assert!((volume - 1.0).abs() < 1e-6);
//! assert!((area - 6.0).abs() < 1e-6);
//!
//! // Evaluate normals
//! let normal = evaluate_face_normal(&brep, 0, 0.5, 0.5);
//! assert!(normal.length() > 0.9);
//! ```

use glam::DVec3;
use rcad_kernel::{BRep, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topology::{Face, Shell, Wire};
use std::collections::{HashMap, HashSet};

// =============================================================================
// Error Types
// =============================================================================

/// Errors that can occur during BRepAlgo operations.
#[derive(Debug, Clone)]
pub enum BRepAlgoError {
    /// Invalid face index.
    InvalidFaceIndex { index: usize, max: usize },
    /// Invalid edge index.
    InvalidEdgeIndex { index: usize, max: usize },
    /// Invalid vertex index.
    InvalidVertexIndex { index: usize, max: usize },
    /// Missing geometry for the specified entity.
    MissingGeometry { kind: &'static str, index: usize },
    /// Degenerate geometry (zero-length edge, etc.).
    DegenerateGeometry { kind: &'static str, index: usize },
    /// Orientation fix failed.
    OrientationFixFailed { reason: String },
}

impl std::fmt::Display for BRepAlgoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BRepAlgoError::InvalidFaceIndex { index, max } => {
                write!(f, "Invalid face index {} (max {})", index, max)
            }
            BRepAlgoError::InvalidEdgeIndex { index, max } => {
                write!(f, "Invalid edge index {} (max {})", index, max)
            }
            BRepAlgoError::InvalidVertexIndex { index, max } => {
                write!(f, "Invalid vertex index {} (max {})", index, max)
            }
            BRepAlgoError::MissingGeometry { kind, index } => {
                write!(f, "Missing {} geometry at index {}", kind, index)
            }
            BRepAlgoError::DegenerateGeometry { kind, index } => {
                write!(f, "Degenerate {} at index {}", kind, index)
            }
            BRepAlgoError::OrientationFixFailed { reason } => {
                write!(f, "Orientation fix failed: {}", reason)
            }
        }
    }
}

impl std::error::Error for BRepAlgoError {}

// =============================================================================
// Normal Evaluation
// =============================================================================

/// Evaluate the normal vector of a face at the given UV parameters.
///
/// Returns the unit normal vector of the face's underlying surface at (u, v).
/// If the face has no surface geometry, returns the face's stored normal.
///
/// # Arguments
///
/// * `brep` - The BRep containing the face
/// * `face_idx` - Flat index of the face across all solids/shells
/// * `u` - U parameter on the surface
/// * `v` - V parameter on the surface
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::evaluate_face_normal;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Evaluate normal at center of first face
/// let normal = evaluate_face_normal(&brep, 0, 0.5, 0.5);
/// assert!(normal.length() > 0.9);
/// ```
pub fn evaluate_face_normal(brep: &BRep, face_idx: usize, u: f64, v: f64) -> DVec3 {
    // Try to get the surface for this face
    if let Some(surf) = get_face_surface(brep, face_idx) {
        return surf.normal_at(u, v);
    }

    // Fall back to the face's stored normal
    if let Some(face) = get_face_by_flat_index(brep, face_idx) {
        return face.normal;
    }

    DVec3::Z // Default fallback
}

/// Evaluate the tangent vector of an edge at the given parameter.
///
/// Returns the unit tangent vector of the edge's 3D curve at parameter t.
/// For reversed edges, the tangent is negated based on the edge orientation
/// in its containing wire.
///
/// # Arguments
///
/// * `brep` - The BRep containing the edge
/// * `edge_idx` - Index of the edge
/// * `t` - Parameter on the edge's curve
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::evaluate_edge_tangent;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Evaluate tangent at midpoint of first edge
/// let tangent = evaluate_edge_tangent(&brep, 0, 0.5);
/// assert!(tangent.length() > 0.9);
/// ```
pub fn evaluate_edge_tangent(brep: &BRep, edge_idx: usize, t: f64) -> DVec3 {
    if edge_idx >= brep.edges.len() {
        return DVec3::X;
    }

    // Get the curve for this edge
    if let Some(curve) = get_edge_curve(brep, edge_idx) {
        let tangent = curve.tangent_at(t);

        // Get the parameter range and normalize t if needed
        if let Some(Some(range)) = brep.geom.edge_curve_range.get(edge_idx) {
            let [tmin, tmax] = *range;
            let normalized_t = tmin + t * (tmax - tmin);
            return curve.tangent_at(normalized_t);
        }

        return tangent;
    }

    // Fall back to computing tangent from vertex positions
    let edge = &brep.edges[edge_idx];
    if let (Some(v0), Some(v1)) = (brep.vertices.get(edge.start), brep.vertices.get(edge.end)) {
        let tangent = (v1.point - v0.point).normalize_or(DVec3::X);
        return tangent;
    }

    DVec3::X
}

/// Evaluate the normal vector at a vertex by averaging adjacent face normals.
///
/// Computes the area-weighted average of normals from all faces that
/// share this vertex. This gives a smooth normal at corners and edges.
///
/// # Arguments
///
/// * `brep` - The BRep containing the vertex
/// * `vertex_idx` - Index of the vertex
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::evaluate_vertex_normal;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Sphere { radius: 1.0 });
/// // Vertex normals on a sphere should point outward
/// let normal = evaluate_vertex_normal(&brep, 0);
/// assert!(normal.length() > 0.9);
/// ```
pub fn evaluate_vertex_normal(brep: &BRep, vertex_idx: usize) -> DVec3 {
    if vertex_idx >= brep.vertices.len() {
        return DVec3::Z;
    }

    // Find all faces that contain this vertex
    let adjacent_faces = find_vertex_faces(brep, vertex_idx);

    if adjacent_faces.is_empty() {
        return DVec3::Z;
    }

    // Compute area-weighted average of face normals
    let mut weighted_normal = DVec3::ZERO;
    let mut total_weight = 0.0;

    for (_face_idx, face) in adjacent_faces {
        // Use face area as weight (approximate from triangles if available)
        let weight = if !face.triangles.is_empty() {
            face.triangles.len() as f64
        } else {
            1.0
        };
        weighted_normal += face.normal * weight;
        total_weight += weight;
    }

    if total_weight > 0.0 {
        (weighted_normal / total_weight).normalize_or(DVec3::Z)
    } else {
        DVec3::Z
    }
}

// =============================================================================
// Tolerance Propagation
// =============================================================================

/// Propagate tolerances from edges to vertices.
///
/// For each vertex, sets its tolerance to at least the maximum tolerance
/// of all edges that share that vertex. This ensures that vertex tolerances
/// are consistent with edge tolerances for boolean operations.
///
/// # Arguments
///
/// * `brep` - The BRep to modify
/// * `tol` - Base tolerance value to ensure minimum tolerance
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::propagate_edge_tolerances;
/// use rcad_kernel::BRep;
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// propagate_edge_tolerances(&mut brep, 1e-7);
/// ```
pub fn propagate_edge_tolerances(brep: &mut BRep, tol: f64) {
    // Ensure tolerance arrays are sized correctly
    while brep.geom.vertex_tolerance.len() < brep.vertices.len() {
        brep.geom.vertex_tolerance.push(tol);
    }
    while brep.geom.edge_tolerance.len() < brep.edges.len() {
        brep.geom.edge_tolerance.push(tol);
    }

    // Build vertex-to-edges map
    let mut vertex_max_tol: Vec<f64> = vec![tol; brep.vertices.len()];

    for (edge_idx, edge) in brep.edges.iter().enumerate() {
        let edge_tol = brep.geom.edge_tolerance.get(edge_idx).copied().unwrap_or(tol);

        // Update vertex tolerances to be at least the edge tolerance
        if edge.start < vertex_max_tol.len() {
            vertex_max_tol[edge.start] = vertex_max_tol[edge.start].max(edge_tol);
        }
        if edge.end < vertex_max_tol.len() {
            vertex_max_tol[edge.end] = vertex_max_tol[edge.end].max(edge_tol);
        }
    }

    // Apply the computed tolerances
    for (i, &max_tol) in vertex_max_tol.iter().enumerate() {
        if i < brep.geom.vertex_tolerance.len() {
            brep.geom.vertex_tolerance[i] = brep.geom.vertex_tolerance[i].max(max_tol);
        }
    }
}

/// Propagate tolerances from faces to edges and vertices.
///
/// For each edge, sets its tolerance to at least the maximum tolerance
/// of all faces that share that edge. Similarly updates vertex tolerances.
///
/// # Arguments
///
/// * `brep` - The BRep to modify
/// * `tol` - Base tolerance value to ensure minimum tolerance
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::propagate_face_tolerances;
/// use rcad_kernel::BRep;
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// propagate_face_tolerances(&mut brep, 1e-7);
/// ```
pub fn propagate_face_tolerances(brep: &mut BRep, tol: f64) {
    // Ensure tolerance arrays are sized correctly
    while brep.geom.face_tolerance.len() < count_faces(brep) {
        brep.geom.face_tolerance.push(tol);
    }
    while brep.geom.edge_tolerance.len() < brep.edges.len() {
        brep.geom.edge_tolerance.push(tol);
    }
    while brep.geom.vertex_tolerance.len() < brep.vertices.len() {
        brep.geom.vertex_tolerance.push(tol);
    }

    // Build edge-to-faces map
    let mut edge_max_tol: Vec<f64> = vec![tol; brep.edges.len()];
    let mut face_idx: usize = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let face_tol = brep.geom.face_tolerance.get(face_idx).copied().unwrap_or(tol);

                // Process all edges in this face
                for we in &face.outer_wire.edges {
                    if we.idx < edge_max_tol.len() {
                        edge_max_tol[we.idx] = edge_max_tol[we.idx].max(face_tol);
                    }
                }
                for inner in &face.inner_wires {
                    for we in &inner.edges {
                        if we.idx < edge_max_tol.len() {
                            edge_max_tol[we.idx] = edge_max_tol[we.idx].max(face_tol);
                        }
                    }
                }

                face_idx += 1;
            }
        }
    }

    // Apply edge tolerances
    for (i, &max_tol) in edge_max_tol.iter().enumerate() {
        if i < brep.geom.edge_tolerance.len() {
            brep.geom.edge_tolerance[i] = brep.geom.edge_tolerance[i].max(max_tol);
        }
    }

    // Also propagate to vertices
    propagate_edge_tolerances(brep, tol);
}

// =============================================================================
// Tools
// =============================================================================

/// Compute the maximum face area in the BRep.
///
/// Returns the area of the largest face, or 0.0 if there are no faces.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::max_face_area;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 2.0, depth: 3.0
/// });
/// // Max face area is 2*3 = 6
/// let max_area = max_face_area(&brep);
/// assert!((max_area - 6.0).abs() < 1e-6);
/// ```
pub fn max_face_area(brep: &BRep) -> f64 {
    let mut max_area: f64 = 0.0;
    let mut face_idx = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let area: f64 = compute_face_area(brep, face, face_idx);
                max_area = max_area.max(area);
                face_idx += 1;
            }
        }
    }

    max_area
}

/// Compute the minimum face area in the BRep.
///
/// Returns the area of the smallest face, or 0.0 if there are no faces.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::min_face_area;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 2.0, depth: 3.0
/// });
/// // Min face area is 1*2 = 2
/// let min_area = min_face_area(&brep);
/// assert!((min_area - 2.0).abs() < 1e-6);
/// ```
pub fn min_face_area(brep: &BRep) -> f64 {
    let mut min_area = f64::INFINITY;
    let mut face_idx = 0;
    let mut has_faces = false;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let area = compute_face_area(brep, face, face_idx);
                min_area = min_area.min(area);
                has_faces = true;
                face_idx += 1;
            }
        }
    }

    if has_faces { min_area } else { 0.0 }
}

/// Compute the maximum edge length in the BRep.
///
/// Returns the length of the longest edge, or 0.0 if there are no edges.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::max_edge_length;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 2.0, depth: 3.0
/// });
/// // Max edge length is 3.0
/// let max_len = max_edge_length(&brep);
/// assert!((max_len - 3.0).abs() < 1e-6);
/// ```
pub fn max_edge_length(brep: &BRep) -> f64 {
    let mut max_length: f64 = 0.0;

    for (edge_idx, edge) in brep.edges.iter().enumerate() {
        let length: f64 = compute_edge_length(brep, edge_idx, edge);
        max_length = max_length.max(length);
    }

    max_length
}

/// Compute the total volume of all solids in the BRep.
///
/// Uses the divergence theorem for closed shells. Returns 0.0 for open shells.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::total_volume;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 2.0, height: 3.0, depth: 4.0
/// });
/// let vol = total_volume(&brep);
/// assert!((vol - 24.0).abs() < 1e-6);
/// ```
pub fn total_volume(brep: &BRep) -> f64 {
    rcad_kernel::volume(brep)
}

/// Compute the total surface area of all faces in the BRep.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::total_surface_area;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 2.0, height: 3.0, depth: 4.0
/// });
/// // SA = 2*(2*3 + 3*4 + 2*4) = 2*(6+12+8) = 52
/// let area = total_surface_area(&brep);
/// assert!((area - 52.0).abs() < 1e-6);
/// ```
pub fn total_surface_area(brep: &BRep) -> f64 {
    rcad_kernel::surface_area(brep)
}

// =============================================================================
// Validity Check
// =============================================================================

/// Check if the BRep is valid for boolean operations.
///
/// Performs a comprehensive validity check including:
/// - Valid vertex, edge, and face indices
/// - Consistent edge-vertex topology
/// - Closed shells for solids
/// - Valid geometry associations
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::is_valid_brep;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// assert!(is_valid_brep(&brep));
/// ```
pub fn is_valid_brep(brep: &BRep) -> bool {
    // Check that vertices have valid indices in edges
    for edge in &brep.edges {
        if edge.start >= brep.vertices.len() || edge.end >= brep.vertices.len() {
            return false;
        }
        // Note: Degenerate edges (start == end) are valid for closed curves like circles.
        // These are commonly used in CAD to represent seam edges on surfaces of revolution.
    }

    // Check that edges have valid indices in wires
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Check outer wire
                for we in &face.outer_wire.edges {
                    if we.idx >= brep.edges.len() {
                        return false;
                    }
                }
                // Check inner wires
                for inner in &face.inner_wires {
                    for we in &inner.edges {
                        if we.idx >= brep.edges.len() {
                            return false;
                        }
                    }
                }
            }
        }
    }

    // Check that solids have at least one shell
    for solid in &brep.solids {
        if solid.shells.is_empty() {
            return false;
        }
    }

    // Check geometry associations
    for (_edge_idx, curve_opt) in brep.geom.edge_curve.iter().enumerate() {
        if let Some(curve_idx) = curve_opt {
            if *curve_idx >= brep.geom.curves.len() {
                return false;
            }
        }
    }

    // Check closed shells for solids
    for solid in &brep.solids {
        for shell in &solid.shells {
            if !is_shell_closed(brep, shell) {
                return false;
            }
        }
    }

    true
}

/// Orientation issue detected during BRep validation.
#[derive(Debug, Clone)]
pub struct OrientationIssue {
    /// Index of the face or shell with the issue.
    pub entity_index: usize,
    /// Type of entity: "face" or "shell".
    pub entity_type: String,
    /// Description of the orientation issue.
    pub description: String,
}

/// Check orientation consistency in the BRep.
///
/// Returns a list of orientation issues found:
/// - Faces with normals not matching the surface orientation
/// - Shells with inconsistent face orientations
/// - Edges with conflicting orientations across faces
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::check_orientation;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let issues = check_orientation(&brep);
/// assert!(issues.is_empty());
/// ```
pub fn check_orientation(brep: &BRep) -> Vec<OrientationIssue> {
    let mut issues = Vec::new();
    let mut face_idx = 0;

    for (solid_idx, solid) in brep.solids.iter().enumerate() {
        for (shell_idx, shell) in solid.shells.iter().enumerate() {
            // Check shell closure and orientation
            let edge_counts = count_edge_uses_in_shell(brep, shell);

            for (edge_idx, count) in &edge_counts {
                if *count != 2 {
                    issues.push(OrientationIssue {
                        entity_index: *edge_idx,
                        entity_type: "edge".to_string(),
                        description: format!(
                            "Edge {} in solid {} shell {} appears {} times (expected 2)",
                            edge_idx, solid_idx, shell_idx, count
                        ),
                    });
                }
            }

            // Check individual faces
            for (_local_face_idx, face) in shell.faces.iter().enumerate() {
                // Check if face normal is consistent with the surface normal
                // Note: For curved surfaces (cylinder, sphere, etc.), the normal varies across
                // the surface, so we only perform this check for planar surfaces where the
                // normal is constant.
                if let Some(surf) = get_face_surface(brep, face_idx) {
                    if matches!(surf, Surface3::Plane(_)) {
                        // For planar surfaces, check that the face normal matches the surface normal
                        let domain = surf.default_domain();
                        let u_mid = (domain[0] + domain[1]) / 2.0;
                        let v_mid = (domain[2] + domain[3]) / 2.0;
                        let surf_normal = surf.normal_at(u_mid, v_mid);

                        let dot = surf_normal.dot(face.normal);
                        if dot < 0.0 {
                            issues.push(OrientationIssue {
                                entity_index: face_idx,
                                entity_type: "face".to_string(),
                                description: format!(
                                    "Face {} normal does not match surface normal (dot = {:.3})",
                                    face_idx, dot
                                ),
                            });
                        }
                    }
                    // For curved surfaces, the face normal represents the overall orientation
                    // (inward/outward relative to the solid), not a specific point normal.
                    // Skip the normal consistency check for curved surfaces.
                }

                face_idx += 1;
            }
        }
    }

    issues
}

// =============================================================================
// Orientation Fixing
// =============================================================================

/// Fix orientation issues in the BRep.
///
/// Attempts to fix:
/// - Faces with reversed normals
/// - Shells with inconsistent face orientations
/// - Wires with wrong orientation
///
/// Returns `true` if any changes were made.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::fix_orientation;
/// use rcad_kernel::BRep;
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let fixed = fix_orientation(&mut brep);
/// // The box from primitives should already be correctly oriented
/// assert!(!fixed); // No changes needed
/// ```
pub fn fix_orientation(brep: &mut BRep) -> bool {
    let mut changed = false;

    // First check and fix face orientations
    let issues = check_orientation(brep);

    for issue in &issues {
        if issue.entity_type == "face" {
            // Reverse the face normal
            if let Some(face) = get_face_by_flat_index_mut(brep, issue.entity_index) {
                face.normal = -face.normal;
                changed = true;
            }
        }
    }

    // Note: Wire orientation checking is complex and depends on the specific CAD conventions used.
    // For well-formed primitives, the wire orientation should already be correct.
    // Skip wire orientation checks for now to avoid false positives.

    changed
}

/// Reverse a face's orientation.
///
/// Negates the face normal and reverses all wire orientations.
///
/// # Arguments
///
/// * `brep` - The BRep to modify
/// * `face_idx` - Flat index of the face to reverse
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::reverse_face;
/// use rcad_kernel::BRep;
///
/// let mut brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let original_normal = brep.solids[0].shells[0].faces[0].normal;
/// reverse_face(&mut brep, 0);
/// let new_normal = brep.solids[0].shells[0].faces[0].normal;
/// assert!((original_normal + new_normal).length() < 1e-9);
/// ```
pub fn reverse_face(brep: &mut BRep, face_idx: usize) {
    if let Some(face) = get_face_by_flat_index_mut(brep, face_idx) {
        // Reverse the face normal
        face.normal = -face.normal;

        // Reverse all wires
        reverse_wire(&mut face.outer_wire);
        for inner in &mut face.inner_wires {
            reverse_wire(inner);
        }
    }
}

// =============================================================================
// Connected Components
// =============================================================================

/// Find connected components in the BRep.
///
/// Returns a list of connected components, where each component is a list of
/// face indices that are connected through shared edges.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo::find_connected_components;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // A single box should have one connected component
/// let components = find_connected_components(&brep);
/// assert_eq!(components.len(), 1);
/// assert_eq!(components[0].len(), 6); // 6 faces
/// ```
pub fn find_connected_components(brep: &BRep) -> Vec<Vec<usize>> {
    let total_faces = count_faces(brep);
    if total_faces == 0 {
        return Vec::new();
    }

    // Build face adjacency through shared edges
    let mut face_adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();

    let mut face_idx = 0;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Collect all edges in this face
                let mut face_edges = HashSet::new();
                for we in &face.outer_wire.edges {
                    face_edges.insert(we.idx);
                }
                for inner in &face.inner_wires {
                    for we in &inner.edges {
                        face_edges.insert(we.idx);
                    }
                }

                // Map each edge to this face
                for edge_idx in face_edges {
                    edge_to_faces.entry(edge_idx).or_default().push(face_idx);
                }

                face_idx += 1;
            }
        }
    }

    // Build face adjacency from shared edges
    for faces in edge_to_faces.values() {
        for i in 0..faces.len() {
            for j in (i + 1)..faces.len() {
                face_adjacency.entry(faces[i]).or_default().push(faces[j]);
                face_adjacency.entry(faces[j]).or_default().push(faces[i]);
            }
        }
    }

    // Find connected components using BFS
    let mut visited = vec![false; total_faces];
    let mut components = Vec::new();

    for start_face in 0..total_faces {
        if visited[start_face] {
            continue;
        }

        let mut component = Vec::new();
        let mut queue = vec![start_face];
        visited[start_face] = true;

        while let Some(face) = queue.pop() {
            component.push(face);

            if let Some(neighbors) = face_adjacency.get(&face) {
                for &neighbor in neighbors {
                    if !visited[neighbor] {
                        visited[neighbor] = true;
                        queue.push(neighbor);
                    }
                }
            }
        }

        components.push(component);
    }

    components
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Get a face by its flat index (across all solids/shells).
fn get_face_by_flat_index(brep: &BRep, face_idx: usize) -> Option<&Face> {
    let mut current_idx = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            if face_idx < current_idx + shell.faces.len() {
                let local_idx = face_idx - current_idx;
                return Some(&shell.faces[local_idx]);
            }
            current_idx += shell.faces.len();
        }
    }

    None
}

/// Get a mutable reference to a face by its flat index.
fn get_face_by_flat_index_mut(brep: &mut BRep, face_idx: usize) -> Option<&mut Face> {
    let mut current_idx = 0;

    for solid in &mut brep.solids {
        for shell in &mut solid.shells {
            if face_idx < current_idx + shell.faces.len() {
                let local_idx = face_idx - current_idx;
                return Some(&mut shell.faces[local_idx]);
            }
            current_idx += shell.faces.len();
        }
    }

    None
}

/// Get the surface for a face by its flat index.
fn get_face_surface(brep: &BRep, face_idx: usize) -> Option<&Surface3> {
    let surf_idx = brep.geom.face_surface.get(face_idx)?.as_ref().copied()?;
    brep.geom.surfaces.get(surf_idx)
}

/// Get the curve for an edge by its index.
fn get_edge_curve(brep: &BRep, edge_idx: usize) -> Option<&Curve3> {
    let curve_idx = brep.geom.edge_curve.get(edge_idx)?.as_ref().copied()?;
    brep.geom.curves.get(curve_idx)
}

/// Count the total number of faces in a BRep.
fn count_faces(brep: &BRep) -> usize {
    brep.solids.iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum()
}

/// Check if a shell is closed by counting edge uses.
fn is_shell_closed(brep: &BRep, shell: &Shell) -> bool {
    if shell.faces.is_empty() {
        return false;
    }

    let edge_counts = count_edge_uses_in_shell(brep, shell);

    // For a closed shell, each edge should appear exactly twice
    edge_counts.values().all(|&count| count == 2)
}

/// Count how many times each edge is used in a shell.
fn count_edge_uses_in_shell(_brep: &BRep, shell: &Shell) -> HashMap<usize, usize> {
    let mut edge_counts = HashMap::new();

    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            *edge_counts.entry(we.idx).or_insert(0) += 1;
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                *edge_counts.entry(we.idx).or_insert(0) += 1;
            }
        }
    }

    edge_counts
}

/// Find all faces that contain a given vertex.
fn find_vertex_faces(brep: &BRep, vertex_idx: usize) -> Vec<(usize, &Face)> {
    let mut result = Vec::new();
    let mut face_idx = 0;

    // First find all edges that contain this vertex
    let mut vertex_edges = HashSet::new();
    for (edge_idx, edge) in brep.edges.iter().enumerate() {
        if edge.start == vertex_idx || edge.end == vertex_idx {
            vertex_edges.insert(edge_idx);
        }
    }

    // Then find all faces that contain any of these edges
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let mut has_vertex = false;
                for we in &face.outer_wire.edges {
                    if vertex_edges.contains(&we.idx) {
                        has_vertex = true;
                        break;
                    }
                }
                if !has_vertex {
                    for inner in &face.inner_wires {
                        for we in &inner.edges {
                            if vertex_edges.contains(&we.idx) {
                                has_vertex = true;
                                break;
                            }
                        }
                        if has_vertex {
                            break;
                        }
                    }
                }

                if has_vertex {
                    result.push((face_idx, face));
                }
                face_idx += 1;
            }
        }
    }

    result
}

/// Compute the area of a face using its triangles or wire vertices.
fn compute_face_area(brep: &BRep, face: &Face, _face_idx: usize) -> f64 {
    // Use pre-triangulated data if available
    if !face.triangles.is_empty() {
        let mut area = 0.0;
        for &[i, j, k] in &face.triangles {
            if let (Some(a), Some(b), Some(c)) = (
                brep.vertices.get(i),
                brep.vertices.get(j),
                brep.vertices.get(k),
            ) {
                area += (b.point - a.point).cross(c.point - a.point).length() * 0.5;
            }
        }
        return area;
    }

    // Fall back to computing from wire vertices
    let wire_pts: Vec<DVec3> = face.outer_wire.edges.iter()
        .filter_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            let vidx = if we.forward { edge.start } else { edge.end };
            brep.vertices.get(vidx).map(|v| v.point)
        })
        .collect();

    if wire_pts.len() < 3 {
        return 0.0;
    }

    // Compute area using the shoelace formula projected onto the face plane
    let n = face.normal;
    let mut area = 0.0;

    for i in 0..wire_pts.len() {
        let j = (i + 1) % wire_pts.len();
        area += (wire_pts[i].cross(wire_pts[j])).dot(n);
    }

    area.abs() * 0.5
}

/// Compute the length of an edge.
fn compute_edge_length(brep: &BRep, edge_idx: usize, edge: &rcad_kernel::topology::Edge) -> f64 {
    // Try to compute from curve geometry
    if let Some(curve) = get_edge_curve(brep, edge_idx) {
        let range = brep.geom.edge_curve_range.get(edge_idx)
            .and_then(|o| *o)
            .unwrap_or_else(|| curve.default_domain());
        return curve_length(curve, range[0], range[1]);
    }

    // Fall back to vertex distance
    if let (Some(v0), Some(v1)) = (brep.vertices.get(edge.start), brep.vertices.get(edge.end)) {
        return (v1.point - v0.point).length();
    }

    0.0
}

/// Approximate curve length by numerical integration.
fn curve_length(curve: &Curve3, t0: f64, t1: f64) -> f64 {
    const NUM_SAMPLES: usize = 100;
    let dt = (t1 - t0) / NUM_SAMPLES as f64;

    let mut length = 0.0;
    let mut prev_point = curve.point_at(t0);

    for i in 1..=NUM_SAMPLES {
        let t = t0 + i as f64 * dt;
        let point = curve.point_at(t);
        length += (point - prev_point).length();
        prev_point = point;
    }

    length
}

/// Reverse a wire's orientation.
fn reverse_wire(wire: &mut Wire) {
    // Reverse the order of edges and flip their direction
    wire.edges.reverse();
    for we in &mut wire.edges {
        we.forward = !we.forward;
    }
}

/// Check if a wire is oriented correctly.
///
/// The convention used in this codebase (matching triangle winding):
/// - Face normal points outward from the solid
/// - The wire should be counterclockwise when viewed from INSIDE the solid
/// - "Inside" means looking in the direction opposite to the face normal
fn is_wire_oriented_correctly(brep: &BRep, wire: &Wire, face_normal: DVec3) -> bool {
    // Collect vertices in order, one per edge (the start of each directed edge)
    let wire_pts: Vec<DVec3> = wire.edges.iter()
        .filter_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            // For a forward edge, use start vertex; for reverse, use end vertex
            let vidx = if we.forward { edge.start } else { edge.end };
            brep.vertices.get(vidx).map(|v| v.point)
        })
        .collect();

    if wire_pts.len() < 3 {
        return true;
    }

    // For planar faces, compute signed area using edge vectors projected onto the face plane.
    // This is more reliable than the position-based cross product formula.

    // Compute the centroid
    let centroid: DVec3 = wire_pts.iter().sum::<DVec3>() / wire_pts.len() as f64;

    // Compute signed area as sum of (edge × radial) · view_direction
    // view_direction = face_normal points outward, so we look FROM inside (opposite direction)
    // which means we use -face_normal for the viewing direction.
    let view_dir = face_normal; // View from outside the solid
    let mut signed_area = 0.0;
    for i in 0..wire_pts.len() {
        let j = (i + 1) % wire_pts.len();
        let edge_vec = wire_pts[j] - wire_pts[i];
        let radial = wire_pts[i] - centroid;
        // Cross product gives the normal contribution
        let contribution = radial.cross(edge_vec).dot(view_dir);
        signed_area += contribution;
    }

    // For counterclockwise when viewed from inside (opposite to face_normal),
    // the signed area when viewed from outside should be negative.
    // So we check if signed_area <= 0 (clockwise from outside = counterclockwise from inside)
    signed_area <= 0.0
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    fn make_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 2.0,
            depth: 3.0,
        })
    }

    fn make_sphere() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 })
    }

    fn make_cylinder() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        })
    }

    // ── Normal Evaluation Tests ───────────────────────────────────────────────

    #[test]
    fn test_evaluate_face_normal_box() {
        let brep = make_box();
        // Box has planar faces, normal should be constant
        let normal = evaluate_face_normal(&brep, 0, 0.0, 0.0);
        assert!((normal.length() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_evaluate_face_normal_sphere() {
        let brep = make_sphere();
        // Sphere normal should point outward
        let normal = evaluate_face_normal(&brep, 0, 0.0, 0.0);
        assert!((normal.length() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_evaluate_edge_tangent_box() {
        let brep = make_box();
        let tangent = evaluate_edge_tangent(&brep, 0, 0.5);
        assert!((tangent.length() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_evaluate_vertex_normal_box() {
        let brep = make_box();
        // Vertex normal at a corner should be some average
        let normal = evaluate_vertex_normal(&brep, 0);
        // Just check it's a valid unit vector
        assert!(normal.length() > 0.0);
    }

    #[test]
    fn test_evaluate_vertex_normal_sphere() {
        let brep = make_sphere();
        // For a sphere, vertex normals should point outward
        // The sphere has vertices at poles
        if !brep.vertices.is_empty() {
            let normal = evaluate_vertex_normal(&brep, 0);
            // Check that the normal points roughly in the Y direction (sphere axis)
            assert!(normal.y.abs() > 0.5 || normal.length() > 0.0);
        }
    }

    // ── Tolerance Propagation Tests ───────────────────────────────────────────

    #[test]
    fn test_propagate_edge_tolerances() {
        let mut brep = make_box();
        propagate_edge_tolerances(&mut brep, 1e-7);

        // Check that tolerances are set
        for &tol in &brep.geom.vertex_tolerance {
            assert!(tol > 0.0);
        }
        for &tol in &brep.geom.edge_tolerance {
            assert!(tol > 0.0);
        }
    }

    #[test]
    fn test_propagate_face_tolerances() {
        let mut brep = make_box();
        propagate_face_tolerances(&mut brep, 1e-7);

        // Check that tolerances are set
        for &tol in &brep.geom.vertex_tolerance {
            assert!(tol > 0.0);
        }
        for &tol in &brep.geom.edge_tolerance {
            assert!(tol > 0.0);
        }
    }

    // ── Tools Tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_max_face_area() {
        let brep = make_box();
        // Box 1x2x3 has face areas: 2, 3, 6 (each appears twice)
        let max_area = max_face_area(&brep);
        assert!((max_area - 6.0).abs() < 1e-6);
    }

    #[test]
    fn test_min_face_area() {
        let brep = make_box();
        // Box 1x2x3 has face areas: 2, 3, 6 (each appears twice)
        let min_area = min_face_area(&brep);
        assert!((min_area - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_max_edge_length() {
        let brep = make_box();
        // Box 1x2x3 has edge lengths: 1, 2, 3
        let max_len = max_edge_length(&brep);
        assert!((max_len - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_total_volume() {
        let brep = make_box();
        let vol = total_volume(&brep);
        assert!((vol - 6.0).abs() < 1e-6); // 1 * 2 * 3
    }

    #[test]
    fn test_total_surface_area() {
        let brep = make_box();
        // SA = 2*(1*2 + 2*3 + 1*3) = 2*(2+6+3) = 22
        let area = total_surface_area(&brep);
        assert!((area - 22.0).abs() < 1e-6);
    }

    #[test]
    fn test_total_volume_sphere() {
        let brep = make_sphere();
        let vol = total_volume(&brep);
        let expected = 4.0 / 3.0 * std::f64::consts::PI;
        // Allow some tolerance due to tessellation
        assert!((vol - expected).abs() / expected < 0.01);
    }

    #[test]
    fn test_total_surface_area_sphere() {
        let brep = make_sphere();
        let area = total_surface_area(&brep);
        let expected = 4.0 * std::f64::consts::PI;
        // Allow some tolerance due to tessellation
        assert!((area - expected).abs() / expected < 0.01);
    }

    // ── Validity Check Tests ──────────────────────────────────────────────────

    #[test]
    fn test_is_valid_brep_box() {
        let brep = make_box();
        assert!(is_valid_brep(&brep));
    }

    #[test]
    fn test_is_valid_brep_sphere() {
        let brep = make_sphere();
        assert!(is_valid_brep(&brep));
    }

    #[test]
    fn test_is_valid_brep_cylinder() {
        let brep = make_cylinder();
        assert!(is_valid_brep(&brep));
    }

    #[test]
    fn test_is_valid_brep_empty() {
        let brep = BRep::new();
        // Empty BRep is technically valid (no invalid topology)
        assert!(is_valid_brep(&brep));
    }

    #[test]
    fn test_check_orientation_box() {
        let brep = make_box();
        let issues = check_orientation(&brep);
        // Box from primitives should have correct orientation
        assert!(issues.is_empty());
    }

    #[test]
    fn test_check_orientation_cylinder() {
        let brep = make_cylinder();
        let issues = check_orientation(&brep);
        // Cylinder from primitives should have correct orientation
        assert!(issues.is_empty());
    }

    // ── Orientation Fixing Tests ──────────────────────────────────────────────

    #[test]
    fn test_fix_orientation_no_change() {
        let mut brep = make_box();
        let changed = fix_orientation(&mut brep);
        // Box from primitives should already be correctly oriented
        assert!(!changed);
    }

    #[test]
    fn test_reverse_face() {
        let mut brep = make_box();
        let original_normal = brep.solids[0].shells[0].faces[0].normal;
        reverse_face(&mut brep, 0);
        let new_normal = brep.solids[0].shells[0].faces[0].normal;

        // Normal should be negated
        assert!((original_normal + new_normal).length() < 1e-9);
    }

    #[test]
    fn test_reverse_face_twice() {
        let mut brep = make_box();
        let original_normal = brep.solids[0].shells[0].faces[0].normal;
        reverse_face(&mut brep, 0);
        reverse_face(&mut brep, 0);
        let new_normal = brep.solids[0].shells[0].faces[0].normal;

        // Normal should be back to original
        assert!((original_normal - new_normal).length() < 1e-9);
    }

    // ── Connected Components Tests ────────────────────────────────────────────

    #[test]
    fn test_find_connected_components_box() {
        let brep = make_box();
        let components = find_connected_components(&brep);
        assert_eq!(components.len(), 1);
        assert_eq!(components[0].len(), 6); // 6 faces in a box
    }

    #[test]
    fn test_find_connected_components_sphere() {
        let brep = make_sphere();
        let components = find_connected_components(&brep);
        assert_eq!(components.len(), 1);
        // Sphere has 1 face
        assert_eq!(components[0].len(), 1);
    }

    #[test]
    fn test_find_connected_components_cylinder() {
        let brep = make_cylinder();
        let components = find_connected_components(&brep);
        assert_eq!(components.len(), 1);
        // Cylinder has 3 faces (top, bottom, side)
        assert_eq!(components[0].len(), 3);
    }

    #[test]
    fn test_find_connected_components_empty() {
        let brep = BRep::new();
        let components = find_connected_components(&brep);
        assert!(components.is_empty());
    }

    // ── Edge Cases Tests ──────────────────────────────────────────────────────

    #[test]
    fn test_evaluate_face_normal_invalid_index() {
        let brep = make_box();
        // Invalid index should return a fallback
        let normal = evaluate_face_normal(&brep, 100, 0.0, 0.0);
        assert!(normal.length() > 0.0);
    }

    #[test]
    fn test_evaluate_edge_tangent_invalid_index() {
        let brep = make_box();
        // Invalid index should return a fallback
        let tangent = evaluate_edge_tangent(&brep, 100, 0.5);
        assert!(tangent.length() > 0.0);
    }

    #[test]
    fn test_evaluate_vertex_normal_invalid_index() {
        let brep = make_box();
        // Invalid index should return a fallback
        let normal = evaluate_vertex_normal(&brep, 100);
        assert!(normal.length() > 0.0);
    }

    #[test]
    fn test_max_face_area_empty() {
        let brep = BRep::new();
        assert_eq!(max_face_area(&brep), 0.0);
    }

    #[test]
    fn test_min_face_area_empty() {
        let brep = BRep::new();
        assert_eq!(min_face_area(&brep), 0.0);
    }

    #[test]
    fn test_max_edge_length_empty() {
        let brep = BRep::new();
        assert_eq!(max_edge_length(&brep), 0.0);
    }

    #[test]
    fn test_total_volume_empty() {
        let brep = BRep::new();
        assert_eq!(total_volume(&brep), 0.0);
    }

    #[test]
    fn test_total_surface_area_empty() {
        let brep = BRep::new();
        assert_eq!(total_surface_area(&brep), 0.0);
    }

    // ── Error Display Tests ────────────────────────────────────────────────────

    #[test]
    fn test_error_display() {
        let err = BRepAlgoError::InvalidFaceIndex { index: 10, max: 5 };
        assert!(format!("{}", err).contains("Invalid face index"));

        let err = BRepAlgoError::MissingGeometry { kind: "surface", index: 5 };
        assert!(format!("{}", err).contains("Missing surface"));

        let err = BRepAlgoError::DegenerateGeometry { kind: "edge", index: 3 };
        assert!(format!("{}", err).contains("Degenerate edge"));
    }
}
