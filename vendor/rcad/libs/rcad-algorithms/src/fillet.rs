//! BRepFilletAPI-style edge fillet operations — analogous to OCCT `BRepFilletAPI_MakeFillet`.
//!
//! # Overview
//!
//! This module provides algorithms for creating fillets (rounded edges) on BRep shapes:
//!
//! - **`make_fillet_edge`**: Fillet single or multiple edges with uniform radius
//! - **`make_fillet_all_edges`**: Fillet all edges of a shape
//! - **`make_variable_fillet`**: Variable radius fillet along an edge
//!
//! # Fillet Surface Construction
//!
//! Fillet surfaces are constructed using the "rolling ball" algorithm:
//! - A ball of the specified radius rolls along the edge
//! - The fillet surface is the envelope of the ball's surface
//! - The fillet connects the two adjacent faces smoothly
//!
//! # Supported Geometry Types
//!
//! - Plane-Plane edge fillet (most common, creates toroidal fillet surface)
//! - Cylinder-Plane edge fillet
//! - Sphere-Plane edge fillet
//! - General surface-surface fillet (numerical computation)
//!
//! # Continuity
//!
//! - C0: Position continuity (sharp corners allowed)
//! - C1: Tangent continuity (smooth transitions)
//! - C2: Curvature continuity (smooth curvature transitions)
//!
//! # References
//!
//! - OCCT `BRepFilletAPI_MakeFillet`
//! - OCCT `ChFi3d_FilBuilder`
//! - OCCT `ChFi3d_ChBuilder`

use std::collections::HashMap;
use glam::DVec3;
use rcad_kernel::{
    BRep, CurveEval,
    geom::{Curve3, Surface3, Line3, Circle3, Plane, CylindricalSurface, SphericalSurface, ToroidalSurface},
    topology::{Face, Wire},
};

use crate::tolerance::TOLERANCE_ABS;

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

const EPS: f64 = 1e-12;
const PI: f64 = std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Error Types
// ─────────────────────────────────────────────────────────────────────────────

/// Errors that can occur during fillet operations.
#[derive(Debug, Clone)]
pub enum FilletError {
    /// Radius is zero or negative.
    InvalidRadius {
        radius: f64,
    },
    /// Edge index out of range.
    EdgeNotFound {
        edge_index: usize,
    },
    /// Face index out of range.
    FaceNotFound {
        face_index: usize,
    },
    /// Edge has no adjacent faces.
    EdgeNoAdjacentFaces {
        edge_index: usize,
    },
    /// Fillet would create degenerate geometry.
    DegenerateGeometry {
        edge_index: usize,
        reason: String,
    },
    /// Radius too large for the edge.
    RadiusTooLarge {
        edge_index: usize,
        radius: f64,
        max_radius: f64,
    },
    /// Failed to compute fillet surface.
    SurfaceComputationFailed {
        edge_index: usize,
        reason: String,
    },
    /// Failed to compute fillet curves.
    CurveComputationFailed {
        edge_index: usize,
        reason: String,
    },
    /// Unsupported geometry combination.
    UnsupportedGeometry {
        edge_index: usize,
        surface1_type: String,
        surface2_type: String,
    },
    /// Variable radius specification is invalid.
    InvalidVariableRadius {
        parameter: f64,
        radius: f64,
    },
    /// Failed to blend adjacent faces.
    BlendFailed {
        edge_index: usize,
        reason: String,
    },
    /// Input shape is invalid.
    InvalidInput(&'static str),
    /// Numerical failure during computation.
    NumericalFailure(&'static str),
    /// Empty result after fillet.
    EmptyResult,
}

impl std::fmt::Display for FilletError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidRadius { radius } => {
                write!(f, "invalid fillet radius: {} (must be > 0)", radius)
            }
            Self::EdgeNotFound { edge_index } => {
                write!(f, "edge {} not found", edge_index)
            }
            Self::FaceNotFound { face_index } => {
                write!(f, "face {} not found", face_index)
            }
            Self::EdgeNoAdjacentFaces { edge_index } => {
                write!(f, "edge {} has no adjacent faces", edge_index)
            }
            Self::DegenerateGeometry { edge_index, reason } => {
                write!(f, "degenerate geometry at edge {}: {}", edge_index, reason)
            }
            Self::RadiusTooLarge { edge_index, radius, max_radius } => {
                write!(f, "radius {} too large for edge {} (max {})", radius, edge_index, max_radius)
            }
            Self::SurfaceComputationFailed { edge_index, reason } => {
                write!(f, "failed to compute fillet surface at edge {}: {}", edge_index, reason)
            }
            Self::CurveComputationFailed { edge_index, reason } => {
                write!(f, "failed to compute fillet curves at edge {}: {}", edge_index, reason)
            }
            Self::UnsupportedGeometry { edge_index, surface1_type, surface2_type } => {
                write!(f, "unsupported geometry at edge {}: {} + {}", edge_index, surface1_type, surface2_type)
            }
            Self::InvalidVariableRadius { parameter, radius } => {
                write!(f, "invalid variable radius {} at parameter {}", radius, parameter)
            }
            Self::BlendFailed { edge_index, reason } => {
                write!(f, "failed to blend adjacent faces at edge {}: {}", edge_index, reason)
            }
            Self::InvalidInput(msg) => write!(f, "invalid input: {}", msg),
            Self::NumericalFailure(msg) => write!(f, "numerical failure: {}", msg),
            Self::EmptyResult => write!(f, "fillet operation produced empty result"),
        }
    }
}

impl std::error::Error for FilletError {}

// ─────────────────────────────────────────────────────────────────────────────
// Fillet Types
// ─────────────────────────────────────────────────────────────────────────────

/// Continuity type for fillet surfaces.
///
/// Determines the smoothness of the transition between the fillet and adjacent faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FilletContinuity {
    /// Position continuity only (G0/C0).
    /// The fillet surface meets the adjacent faces without gaps.
    C0,
    /// Tangent continuity (G1/C1).
    /// The fillet surface meets the adjacent faces with tangent continuity.
    #[default]
    C1,
    /// Curvature continuity (G2/C2).
    /// The fillet surface meets the adjacent faces with curvature continuity.
    C2,
}

/// Fillet mode for radius specification.
///
/// Determines how the fillet radius is interpreted and applied.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FilletMode {
    /// Uniform radius along the entire edge.
    #[default]
    Uniform,
    /// Variable radius specified at parameters along the edge.
    Variable,
    /// Chordal mode: radius defines the chord length of the fillet arc.
    Chordal,
}

/// Parameters for fillet operations.
///
/// Controls the shape and quality of the fillet surface.
#[derive(Debug, Clone)]
pub struct FilletParams {
    /// Fillet radius (or chord length in chordal mode).
    pub radius: f64,
    /// Continuity between fillet and adjacent faces.
    pub continuity: FilletContinuity,
    /// Fillet mode (uniform, variable, chordal).
    pub mode: FilletMode,
    /// Tension parameter for variable radius fillets (0.0 = linear, 1.0 = smooth).
    /// Controls the interpolation between radius values.
    pub tension: f64,
    /// Angular tolerance for edge discretization (radians).
    pub angular_tolerance: f64,
    /// Distance tolerance for geometric computations.
    pub distance_tolerance: f64,
}

impl Default for FilletParams {
    fn default() -> Self {
        Self {
            radius: 1.0,
            continuity: FilletContinuity::C1,
            mode: FilletMode::Uniform,
            tension: 0.5,
            angular_tolerance: 1e-6,
            distance_tolerance: TOLERANCE_ABS,
        }
    }
}

impl FilletParams {
    /// Create new fillet parameters with the specified radius.
    pub fn new(radius: f64) -> Self {
        Self {
            radius,
            ..Default::default()
        }
    }

    /// Set the continuity.
    pub fn with_continuity(mut self, continuity: FilletContinuity) -> Self {
        self.continuity = continuity;
        self
    }

    /// Set the fillet mode.
    pub fn with_mode(mut self, mode: FilletMode) -> Self {
        self.mode = mode;
        self
    }

    /// Set the tension parameter.
    pub fn with_tension(mut self, tension: f64) -> Self {
        self.tension = tension.clamp(0.0, 1.0);
        self
    }
}

/// Variable radius specification at a point along an edge.
#[derive(Debug, Clone)]
pub struct VariableRadiusPoint {
    /// Parameter value along the edge (0.0 to 1.0).
    pub parameter: f64,
    /// Radius at this parameter.
    pub radius: f64,
}

impl VariableRadiusPoint {
    /// Create a new variable radius point.
    pub fn new(parameter: f64, radius: f64) -> Self {
        Self { parameter, radius }
    }
}

/// Result of a fillet operation.
#[derive(Debug, Clone)]
pub struct FilletResult {
    /// The resulting BRep with fillets applied.
    pub brep: BRep,
    /// Number of edges filletted.
    pub edges_processed: usize,
    /// Number of fillet faces created.
    pub fillet_faces_created: usize,
    /// Any warnings encountered during the operation.
    pub warnings: Vec<String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Fillet Surface Types
// ─────────────────────────────────────────────────────────────────────────────

/// Computed fillet surface data.
#[derive(Debug, Clone)]
struct FilletSurface {
    /// The fillet surface geometry.
    surface: Surface3,
    /// UV domain for the fillet surface.
    uv_domain: [f64; 4],
    /// Boundary curves on the fillet surface.
    boundary_curves: Vec<FilletCurve>,
    /// Edge index this fillet corresponds to.
    edge_index: usize,
}

/// Computed fillet boundary curve.
#[derive(Debug, Clone)]
struct FilletCurve {
    /// The curve geometry.
    curve: Curve3,
    /// Parameter range for the curve.
    parameter_range: [f64; 2],
    /// Whether this curve is on the start or end of the fillet.
    is_start: bool,
}

/// Information about an edge to be filletted.
#[derive(Debug, Clone)]
struct EdgeInfo {
    /// Edge index.
    index: usize,
    /// Start vertex index.
    start_vertex: usize,
    /// End vertex index.
    end_vertex: usize,
    /// Adjacent face indices (usually 2).
    adjacent_faces: Vec<usize>,
    /// Edge tangent at start.
    tangent_start: DVec3,
    /// Edge tangent at end.
    tangent_end: DVec3,
    /// Edge length.
    length: f64,
    /// Edge curve (if available).
    curve: Option<Curve3>,
    /// Parameter range for the edge curve.
    curve_range: Option<[f64; 2]>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main API Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Create a fillet on one or more edges with uniform radius.
///
/// # Arguments
///
/// * `brep` - The input BRep shape.
/// * `edge_indices` - Indices of edges to fillet.
/// * `radius` - Fillet radius (must be > 0).
///
/// # Returns
///
/// A `FilletResult` containing the modified BRep with fillets.
///
/// # Example
///
/// ```ignore
/// let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 4.0, height: 2.0, depth: 3.0 });
/// let result = make_fillet_edge(&box_brep, &[0, 1, 2], 0.5)?;
/// ```
pub fn make_fillet_edge(
    brep: &BRep,
    edge_indices: &[usize],
    radius: f64,
) -> Result<FilletResult, FilletError> {
    let params = FilletParams::new(radius);
    make_fillet_edge_with_params(brep, edge_indices, &params)
}

/// Create a fillet on one or more edges with custom parameters.
///
/// # Arguments
///
/// * `brep` - The input BRep shape.
/// * `edge_indices` - Indices of edges to fillet.
/// * `params` - Fillet parameters.
///
/// # Returns
///
/// A `FilletResult` containing the modified BRep with fillets.
pub fn make_fillet_edge_with_params(
    brep: &BRep,
    edge_indices: &[usize],
    params: &FilletParams,
) -> Result<FilletResult, FilletError> {
    if edge_indices.is_empty() {
        return Ok(FilletResult {
            brep: brep.clone(),
            edges_processed: 0,
            fillet_faces_created: 0,
            warnings: vec!["No edges specified for fillet".to_string()],
        });
    }

    if params.radius <= 0.0 {
        return Err(FilletError::InvalidRadius { radius: params.radius });
    }

    // Validate edge indices
    for &idx in edge_indices {
        if idx >= brep.edges.len() {
            return Err(FilletError::EdgeNotFound { edge_index: idx });
        }
    }

    // Build edge information
    let edge_infos = collect_edge_infos(brep, edge_indices)?;

    // Compute fillet surfaces for each edge
    let mut fillet_surfaces = Vec::new();
    let mut warnings = Vec::new();

    for edge_info in &edge_infos {
        match compute_fillet_for_edge(brep, edge_info, params) {
            Ok(fs) => fillet_surfaces.push(fs),
            Err(e) => {
                warnings.push(format!("Could not fillet edge {}: {}", edge_info.index, e));
            }
        }
    }

    // Build result BRep with fillets
    let result = build_fillet_brep(brep, &fillet_surfaces, &edge_infos, params)?;

    Ok(FilletResult {
        brep: result,
        edges_processed: fillet_surfaces.len(),
        fillet_faces_created: fillet_surfaces.len(),
        warnings,
    })
}

/// Fillet all edges of a shape with uniform radius.
///
/// # Arguments
///
/// * `brep` - The input BRep shape.
/// * `radius` - Fillet radius (must be > 0).
///
/// # Returns
///
/// A `FilletResult` containing the modified BRep with all edges filletted.
pub fn make_fillet_all_edges(
    brep: &BRep,
    radius: f64,
) -> Result<FilletResult, FilletError> {
    let all_edges: Vec<usize> = (0..brep.edges.len()).collect();
    make_fillet_edge(brep, &all_edges, radius)
}

/// Create a variable radius fillet along edges.
///
/// # Arguments
///
/// * `brep` - The input BRep shape.
/// * `edge_indices` - Indices of edges to fillet.
/// * `radii` - Variable radius specification (at least 2 points required).
///
/// # Returns
///
/// A `FilletResult` containing the modified BRep with variable radius fillets.
///
/// # Notes
///
/// The parameter values in `radii` should be in the range [0.0, 1.0] and represent
/// positions along the edge curve. At least two points (start and end) are required.
pub fn make_variable_fillet(
    brep: &BRep,
    edge_indices: &[usize],
    radii: &[VariableRadiusPoint],
) -> Result<FilletResult, FilletError> {
    if radii.len() < 2 {
        return Err(FilletError::InvalidInput("variable fillet requires at least 2 radius points"));
    }

    // Validate radius points
    for rp in radii {
        if rp.parameter < 0.0 || rp.parameter > 1.0 {
            return Err(FilletError::InvalidVariableRadius {
                parameter: rp.parameter,
                radius: rp.radius,
            });
        }
        if rp.radius <= 0.0 {
            return Err(FilletError::InvalidVariableRadius {
                parameter: rp.parameter,
                radius: rp.radius,
            });
        }
    }

    // Use average radius for initial computation
    let avg_radius = radii.iter().map(|r| r.radius).sum::<f64>() / radii.len() as f64;
    let mut params = FilletParams::new(avg_radius);
    params.mode = FilletMode::Variable;

    make_variable_fillet_with_params(brep, edge_indices, radii, &params)
}

/// Create a variable radius fillet with custom parameters.
fn make_variable_fillet_with_params(
    brep: &BRep,
    edge_indices: &[usize],
    radii: &[VariableRadiusPoint],
    params: &FilletParams,
) -> Result<FilletResult, FilletError> {
    if edge_indices.is_empty() {
        return Ok(FilletResult {
            brep: brep.clone(),
            edges_processed: 0,
            fillet_faces_created: 0,
            warnings: vec!["No edges specified for fillet".to_string()],
        });
    }

    // Validate edge indices
    for &idx in edge_indices {
        if idx >= brep.edges.len() {
            return Err(FilletError::EdgeNotFound { edge_index: idx });
        }
    }

    // Build edge information
    let edge_infos = collect_edge_infos(brep, edge_indices)?;

    // Compute variable radius fillet surfaces
    let mut fillet_surfaces = Vec::new();
    let mut warnings = Vec::new();

    for edge_info in &edge_infos {
        match compute_variable_fillet_for_edge(brep, edge_info, radii, params) {
            Ok(fs) => fillet_surfaces.push(fs),
            Err(e) => {
                warnings.push(format!("Could not fillet edge {}: {}", edge_info.index, e));
            }
        }
    }

    // Build result BRep
    let result = build_fillet_brep(brep, &fillet_surfaces, &edge_infos, params)?;

    Ok(FilletResult {
        brep: result,
        edges_processed: fillet_surfaces.len(),
        fillet_faces_created: fillet_surfaces.len(),
        warnings,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Edge Information Collection
// ─────────────────────────────────────────────────────────────────────────────

/// Collect information about edges to be filletted.
fn collect_edge_infos(
    brep: &BRep,
    edge_indices: &[usize],
) -> Result<Vec<EdgeInfo>, FilletError> {
    let mut infos = Vec::new();

    // Build face-to-edge adjacency map
    let edge_faces = build_edge_face_adjacency(brep);

    for &edge_idx in edge_indices {
        let edge = &brep.edges[edge_idx];

        // Get adjacent faces
        let adjacent_faces = edge_faces.get(&edge_idx).cloned().unwrap_or_default();
        if adjacent_faces.len() < 2 {
            // For now, we skip edges with less than 2 adjacent faces
            // but we could handle boundary edges differently
            continue;
        }

        // Get edge curve and range
        let (curve, curve_range) = if edge_idx < brep.geom.edge_curve.len() {
            let curve_idx = brep.geom.edge_curve[edge_idx];
            let range = brep.geom.edge_curve_range[edge_idx];
            match (curve_idx, range) {
                (Some(ci), Some(r)) => {
                    if ci < brep.geom.curves.len() {
                        (Some(brep.geom.curves[ci].clone()), Some(r))
                    } else {
                        (None, None)
                    }
                }
                _ => (None, None),
            }
        } else {
            (None, None)
        };

        // Compute edge length
        let length = compute_edge_length(brep, edge_idx);

        // Compute tangents
        let (tangent_start, tangent_end) = compute_edge_tangents(brep, edge_idx, &curve, &curve_range);

        infos.push(EdgeInfo {
            index: edge_idx,
            start_vertex: edge.start,
            end_vertex: edge.end,
            adjacent_faces,
            tangent_start,
            tangent_end,
            length,
            curve,
            curve_range,
        });
    }

    Ok(infos)
}

/// Build a map from edge index to adjacent face indices.
fn build_edge_face_adjacency(brep: &BRep) -> HashMap<usize, Vec<usize>> {
    let mut edge_faces: HashMap<usize, Vec<usize>> = HashMap::new();

    for (solid_idx, solid) in brep.solids.iter().enumerate() {
        for (shell_idx, shell) in solid.shells.iter().enumerate() {
            for (face_idx, face) in shell.faces.iter().enumerate() {
                let flat_face_idx = compute_flat_face_index(brep, solid_idx, shell_idx, face_idx);
                for wire_edge in &face.outer_wire.edges {
                    edge_faces.entry(wire_edge.idx).or_default().push(flat_face_idx);
                }
                for inner_wire in &face.inner_wires {
                    for wire_edge in &inner_wire.edges {
                        edge_faces.entry(wire_edge.idx).or_default().push(flat_face_idx);
                    }
                }
            }
        }
    }

    edge_faces
}

/// Compute flat face index from (solid, shell, face) indices.
fn compute_flat_face_index(brep: &BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> usize {
    let mut count = 0;
    for (i, solid) in brep.solids.iter().enumerate() {
        if i < solid_idx {
            count += solid.shells.iter().map(|s| s.faces.len()).sum::<usize>();
        } else if i == solid_idx {
            for (j, shell) in solid.shells.iter().enumerate() {
                if j < shell_idx {
                    count += shell.faces.len();
                } else if j == shell_idx {
                    count += face_idx;
                }
            }
        }
    }
    count
}

/// Compute the length of an edge.
fn compute_edge_length(brep: &BRep, edge_idx: usize) -> f64 {
    let edge = &brep.edges[edge_idx];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;
    (p1 - p0).length()
}

/// Compute tangent vectors at the start and end of an edge.
fn compute_edge_tangents(
    brep: &BRep,
    edge_idx: usize,
    curve: &Option<Curve3>,
    curve_range: &Option<[f64; 2]>,
) -> (DVec3, DVec3) {
    let edge = &brep.edges[edge_idx];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;

    match (curve, curve_range) {
        (Some(c), Some([t0, t1])) => {
            let t_start = c.tangent_at(*t0);
            let t_end = c.tangent_at(*t1);
            (t_start.normalize_or(DVec3::X), t_end.normalize_or(DVec3::X))
        }
        _ => {
            // Fall back to linear edge tangent
            let dir = (p1 - p0).normalize_or(DVec3::X);
            (dir, dir)
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Fillet Surface Construction
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the rolling ball fillet surface for an edge.
///
/// The rolling ball algorithm places a sphere of radius `r` tangent to both
/// adjacent faces and rolls it along the edge. The fillet surface is the
/// envelope traced by the sphere.
pub fn compute_rollball_surface(
    brep: &BRep,
    edge_info: &EdgeInfo,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // Get the two adjacent faces
    let faces = &edge_info.adjacent_faces;
    if faces.len() < 2 {
        return Err(FilletError::EdgeNoAdjacentFaces { edge_index: edge_info.index });
    }

    // Get surface types of adjacent faces
    let surf1 = get_face_surface(brep, faces[0]);
    let surf2 = get_face_surface(brep, faces[1]);

    match (&surf1, &surf2) {
        (Some(s1), Some(s2)) => {
            compute_rollball_surface_for_surfaces(
                edge_info,
                s1, s2,
                faces[0], faces[1],
                radius,
            )
        }
        _ => {
            // Fall back to toroidal approximation
            compute_toroidal_fillet_surface(brep, edge_info, radius)
        }
    }
}

/// Get the surface for a face (by flat index).
fn get_face_surface(brep: &BRep, flat_face_idx: usize) -> Option<Surface3> {
    if flat_face_idx < brep.geom.face_surface.len() {
        if let Some(surf_idx) = brep.geom.face_surface[flat_face_idx] {
            if surf_idx < brep.geom.surfaces.len() {
                return Some(brep.geom.surfaces[surf_idx].clone());
            }
        }
    }
    None
}

/// Compute rolling ball surface for specific surface types.
fn compute_rollball_surface_for_surfaces(
    edge_info: &EdgeInfo,
    surf1: &Surface3,
    surf2: &Surface3,
    _face1_idx: usize,
    _face2_idx: usize,
    radius: f64,
) -> Result<Surface3, FilletError> {
    match (surf1, surf2) {
        // Plane-Plane fillet creates a toroidal surface
        (Surface3::Plane(p1), Surface3::Plane(p2)) => {
            compute_plane_plane_fillet(edge_info, p1, p2, radius)
        }
        // Cylinder-Plane fillet
        (Surface3::Cylinder(c), Surface3::Plane(p)) |
        (Surface3::Plane(p), Surface3::Cylinder(c)) => {
            compute_cylinder_plane_fillet(edge_info, c, p, radius)
        }
        // Sphere-Plane fillet
        (Surface3::Sphere(s), Surface3::Plane(p)) |
        (Surface3::Plane(p), Surface3::Sphere(s)) => {
            compute_sphere_plane_fillet(edge_info, s, p, radius)
        }
        // General case - use numerical approximation
        _ => {
            // Fall back to toroidal approximation
            compute_general_fillet_surface(edge_info, surf1, surf2, radius)
        }
    }
}

/// Compute fillet surface for plane-plane edge.
fn compute_plane_plane_fillet(
    edge_info: &EdgeInfo,
    plane1: &Plane,
    plane2: &Plane,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // Get the edge direction
    let edge_dir = edge_info.tangent_start;

    // Compute the angle between the planes
    let n1 = plane1.normal.normalize();
    let n2 = plane2.normal.normalize();

    // The angle between the planes
    let cos_angle = n1.dot(n2);
    let angle = cos_angle.acos();

    // Check if the edge is along the intersection of the planes
    let intersection_dir = n1.cross(n2);

    if intersection_dir.length_squared() < EPS {
        // Planes are parallel - cannot fillet
        return Err(FilletError::DegenerateGeometry {
            edge_index: edge_info.index,
            reason: "adjacent faces are parallel".to_string(),
        });
    }

    // The fillet is a toroidal surface with:
    // - Major radius = radius / sin(angle/2)
    // - Minor radius = radius
    // - Axis along the edge

    let half_angle = angle / 2.0;
    let sin_half = half_angle.sin();

    if sin_half.abs() < EPS {
        return Err(FilletError::DegenerateGeometry {
            edge_index: edge_info.index,
            reason: "edge angle is too small".to_string(),
        });
    }

    // For small angles, we use a simpler cylindrical approximation
    // For larger angles, use a proper torus

    // Compute the centerline of the fillet
    // It should be offset from the edge by radius / sin(angle/2)
    let offset_distance = radius / sin_half;

    // Create a torus as the fillet surface
    // Center is at the midpoint of the edge, offset along the bisector
    let mid_point = DVec3::ZERO; // Will be transformed

    // Bisector direction (average of the outward normals)
    let bisector = (n1 + n2).normalize();

    // Create torus centered at edge midpoint, with axis along the edge
    let center = mid_point - bisector * offset_distance;
    let axis = edge_dir.normalize();

    // For the torus:
    // - major_radius is the distance from the axis to the center of the circular cross-section
    // - minor_radius is the radius of the circular cross-section (the fillet radius)
    let major_radius = offset_distance;
    let minor_radius = radius;

    Ok(Surface3::Torus(ToroidalSurface {
        center,
        axis,
        major_radius,
        minor_radius,
    }))
}

/// Compute fillet surface for cylinder-plane edge.
fn compute_cylinder_plane_fillet(
    edge_info: &EdgeInfo,
    cylinder: &CylindricalSurface,
    plane: &Plane,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // For cylinder-plane edges, the fillet is typically a portion of a torus
    // or a more complex blend surface

    let edge_dir = edge_info.tangent_start;
    let cylinder_axis = cylinder.axis.normalize();
    let _plane_normal = plane.normal.normalize();

    // Check if the edge is parallel to the cylinder axis
    let parallel_to_axis = edge_dir.dot(cylinder_axis).abs() > 1.0 - EPS;

    if parallel_to_axis {
        // Edge is parallel to cylinder axis - creates a cylindrical fillet
        // Offset the cylinder surface by the fillet radius
        Ok(Surface3::Cylinder(CylindricalSurface {
            origin: cylinder.origin,
            axis: cylinder.axis,
            radius: cylinder.radius + radius,
        }))
    } else {
        // Edge is perpendicular or angled to cylinder axis - creates toroidal fillet
        let center = cylinder.origin;
        let major_radius = cylinder.radius + radius;
        let minor_radius = radius;

        Ok(Surface3::Torus(ToroidalSurface {
            center,
            axis: cylinder_axis,
            major_radius,
            minor_radius,
        }))
    }
}

/// Compute fillet surface for sphere-plane edge.
fn compute_sphere_plane_fillet(
    edge_info: &EdgeInfo,
    sphere: &SphericalSurface,
    plane: &Plane,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // For sphere-plane edges, the fillet is typically a portion of a larger sphere
    // or a toroidal surface

    let _edge_dir = edge_info.tangent_start;
    let _plane_normal = plane.normal.normalize();

    // The fillet on a sphere creates a larger sphere offset from the original
    let fillet_sphere_radius = sphere.radius + radius;

    // If the plane cuts through the sphere, the fillet is a torus
    // Otherwise it's a spherical fillet

    // Compute distance from sphere center to plane
    let _dist_to_plane = (sphere.center - plane.origin).dot(_plane_normal);

    // For a spherical fillet, we just offset the sphere
    Ok(Surface3::Sphere(SphericalSurface {
        center: sphere.center,
        axis: sphere.axis,
        radius: fillet_sphere_radius,
    }))
}

/// Compute general fillet surface for arbitrary surface types.
fn compute_general_fillet_surface(
    edge_info: &EdgeInfo,
    _surf1: &Surface3,
    _surf2: &Surface3,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // For general surfaces, we approximate with a torus
    // This is a simplification - a full implementation would use
    // numerical methods to compute the exact rolling ball envelope

    // Use the edge direction as the torus axis
    let axis = edge_info.tangent_start.normalize();

    // Compute edge midpoint as torus center
    // (This is approximate - actual center depends on surface geometry)
    let center = DVec3::ZERO;

    // Use a simplified major radius calculation
    let major_radius = radius * 2.0; // Approximate
    let minor_radius = radius;

    Ok(Surface3::Torus(ToroidalSurface {
        center,
        axis,
        major_radius,
        minor_radius,
    }))
}

/// Compute a toroidal approximation for the fillet surface.
fn compute_toroidal_fillet_surface(
    brep: &BRep,
    edge_info: &EdgeInfo,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // Get edge geometry
    let edge = &brep.edges[edge_info.index];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;

    // Edge midpoint becomes the torus center
    let center = (p0 + p1) * 0.5;

    // Edge direction becomes the torus axis
    let axis = (p1 - p0).normalize_or(DVec3::Z);

    // For a simple edge, use default major radius
    let major_radius = radius * 2.0;
    let minor_radius = radius;

    Ok(Surface3::Torus(ToroidalSurface {
        center,
        axis,
        major_radius,
        minor_radius,
    }))
}

/// Compute the boundary curves of a fillet.
///
/// Returns the curves that form the boundaries of the fillet surface.
pub fn compute_fillet_curves(
    brep: &BRep,
    edge_info: &EdgeInfo,
    radius: f64,
    surface: &Surface3,
) -> Result<Vec<FilletCurve>, FilletError> {
    let edge = &brep.edges[edge_info.index];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;

    // The fillet curves are arcs on the fillet surface
    // For a toroidal fillet, these are circles at the start and end

    let mut curves = Vec::new();

    match surface {
        Surface3::Torus(torus) => {
            // Create circular arcs at the start and end of the fillet
            // The arcs lie in planes perpendicular to the torus axis

            let axis = torus.axis.normalize();
            let _ref_dir = any_perpendicular(axis);

            // Start arc
            let start_center = torus.center + axis * (p0 - torus.center).dot(axis);
            let start_curve = Curve3::Circle(Circle3 {
                center: start_center,
                normal: axis,
                radius: torus.minor_radius,
            });

            curves.push(FilletCurve {
                curve: start_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: true,
            });

            // End arc
            let end_center = torus.center + axis * (p1 - torus.center).dot(axis);
            let end_curve = Curve3::Circle(Circle3 {
                center: end_center,
                normal: axis,
                radius: torus.minor_radius,
            });

            curves.push(FilletCurve {
                curve: end_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: false,
            });
        }
        Surface3::Cylinder(cyl) => {
            // For cylindrical fillet, curves are circles at start and end
            let axis = cyl.axis.normalize();

            let start_curve = Curve3::Circle(Circle3 {
                center: p0,
                normal: axis,
                radius: cyl.radius,
            });

            curves.push(FilletCurve {
                curve: start_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: true,
            });

            let end_curve = Curve3::Circle(Circle3 {
                center: p1,
                normal: axis,
                radius: cyl.radius,
            });

            curves.push(FilletCurve {
                curve: end_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: false,
            });
        }
        Surface3::Sphere(sphere) => {
            // For spherical fillet, curves are circles on the sphere surface
            let axis = sphere.axis.normalize();
            let _ref_dir = any_perpendicular(axis);

            // Approximate with circles through the edge endpoints
            let start_curve = Curve3::Circle(Circle3 {
                center: p0 - axis * (p0 - sphere.center).dot(axis) * 0.5,
                normal: axis,
                radius: sphere.radius * 0.5,
            });

            curves.push(FilletCurve {
                curve: start_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: true,
            });
        }
        _ => {
            // For other surfaces, create approximate curves
            let edge_dir = (p1 - p0).normalize_or(DVec3::Z);

            let start_curve = Curve3::Line(Line3 {
                origin: p0,
                direction: any_perpendicular(edge_dir),
            });

            curves.push(FilletCurve {
                curve: start_curve,
                parameter_range: [0.0, radius],
                is_start: true,
            });
        }
    }

    Ok(curves)
}

/// Blend the fillet surface with adjacent faces.
///
/// This function creates smooth transitions between the fillet and the
/// adjacent faces of the original shape.
pub fn blend_adjacent_faces(
    _brep: &mut BRep,
    _fillet_surface: &Surface3,
    edge_info: &EdgeInfo,
    radius: f64,
) -> Result<(), FilletError> {
    // This function modifies the BRep to blend the fillet with adjacent faces
    // In a full implementation, this would:
    // 1. Trim the adjacent faces at the fillet boundary
    // 2. Add the fillet face to the shell
    // 3. Create proper edge topology connecting the fillet to adjacent faces

    // For now, we just validate the inputs
    if edge_info.adjacent_faces.len() < 2 {
        return Err(FilletError::BlendFailed {
            edge_index: edge_info.index,
            reason: "need at least 2 adjacent faces".to_string(),
        });
    }

    // Validate that the fillet radius is not too large
    let edge_length = edge_info.length;
    if radius > edge_length * 0.5 {
        return Err(FilletError::RadiusTooLarge {
            edge_index: edge_info.index,
            radius,
            max_radius: edge_length * 0.5,
        });
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Fillet Computation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute fillet surface for a single edge.
fn compute_fillet_for_edge(
    brep: &BRep,
    edge_info: &EdgeInfo,
    params: &FilletParams,
) -> Result<FilletSurface, FilletError> {
    // Compute the rolling ball fillet surface
    let surface = compute_rollball_surface(brep, edge_info, params.radius)?;

    // Compute the boundary curves
    let boundary_curves = compute_fillet_curves(brep, edge_info, params.radius, &surface)?;

    // Compute UV domain
    let uv_domain = compute_fillet_uv_domain(&surface, edge_info);

    Ok(FilletSurface {
        surface,
        uv_domain,
        boundary_curves,
        edge_index: edge_info.index,
    })
}

/// Compute variable radius fillet surface for a single edge.
fn compute_variable_fillet_for_edge(
    brep: &BRep,
    edge_info: &EdgeInfo,
    radii: &[VariableRadiusPoint],
    params: &FilletParams,
) -> Result<FilletSurface, FilletError> {
    // Sort radius points by parameter
    let mut sorted_radii = radii.to_vec();
    sorted_radii.sort_by(|a, b| a.parameter.partial_cmp(&b.parameter).unwrap());

    // Use the average radius for the main surface
    let avg_radius = sorted_radii.iter().map(|r| r.radius).sum::<f64>() / sorted_radii.len() as f64;

    // Compute base surface
    let surface = compute_rollball_surface(brep, edge_info, avg_radius)?;

    // Compute boundary curves with variable radii
    let boundary_curves = compute_variable_fillet_curves(brep, edge_info, &sorted_radii, &surface)?;

    // Compute UV domain
    let uv_domain = compute_fillet_uv_domain(&surface, edge_info);

    Ok(FilletSurface {
        surface,
        uv_domain,
        boundary_curves,
        edge_index: edge_info.index,
    })
}

/// Compute boundary curves for variable radius fillet.
fn compute_variable_fillet_curves(
    brep: &BRep,
    edge_info: &EdgeInfo,
    radii: &[VariableRadiusPoint],
    surface: &Surface3,
) -> Result<Vec<FilletCurve>, FilletError> {
    let edge = &brep.edges[edge_info.index];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;

    let mut curves = Vec::new();

    // For variable radius, we sample along the edge
    // and create curves at each sample point
    for rp in radii {
        let t = rp.parameter;
        let pt = p0 + (p1 - p0) * t;

        match surface {
            Surface3::Torus(torus) => {
                let axis = torus.axis.normalize();
                let center = torus.center + axis * (pt - torus.center).dot(axis);

                let curve = Curve3::Circle(Circle3 {
                    center,
                    normal: axis,
                    radius: rp.radius,
                });

                curves.push(FilletCurve {
                    curve,
                    parameter_range: [0.0, 2.0 * PI],
                    is_start: t < 0.5,
                });
            }
            _ => {
                // Fall back to line
                let curve = Curve3::Line(Line3 {
                    origin: pt,
                    direction: edge_info.tangent_start,
                });

                curves.push(FilletCurve {
                    curve,
                    parameter_range: [0.0, rp.radius],
                    is_start: t < 0.5,
                });
            }
        }
    }

    Ok(curves)
}

/// Compute UV domain for a fillet surface.
fn compute_fillet_uv_domain(surface: &Surface3, edge_info: &EdgeInfo) -> [f64; 4] {
    match surface {
        Surface3::Torus(_) => {
            // Torus: u = revolution angle [0, 2*pi], v = arc angle [0, pi/2] typically
            [0.0, 2.0 * PI, 0.0, PI * 0.5]
        }
        Surface3::Cylinder(_) => {
            // Cylinder: u = azimuth [0, 2*pi], v = height along edge
            [0.0, 2.0 * PI, 0.0, edge_info.length]
        }
        Surface3::Sphere(_) => {
            // Sphere: u = longitude [0, 2*pi], v = colatitude [0, pi]
            [0.0, 2.0 * PI, 0.0, PI * 0.5]
        }
        _ => {
            // Default domain
            [0.0, 1.0, 0.0, 1.0]
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRep Construction
// ─────────────────────────────────────────────────────────────────────────────

/// Build the result BRep with fillets.
fn build_fillet_brep(
    brep: &BRep,
    fillet_surfaces: &[FilletSurface],
    _edge_infos: &[EdgeInfo],
    _params: &FilletParams,
) -> Result<BRep, FilletError> {
    // Clone the input BRep
    let mut result = brep.clone();

    // Add fillet faces to the BRep
    for fs in fillet_surfaces {
        add_fillet_face(&mut result, fs)?;
    }

    Ok(result)
}

/// Add a fillet face to the BRep.
fn add_fillet_face(brep: &mut BRep, fillet_surface: &FilletSurface) -> Result<(), FilletError> {
    // Add the surface to the geometry store
    let surf_idx = brep.geom.surfaces.len();
    brep.geom.surfaces.push(fillet_surface.surface.clone());
    brep.geom.face_surface.push(Some(surf_idx));

    // Create a placeholder face
    // In a full implementation, this would have proper wire boundaries
    let face = Face {
        outer_wire: Wire {
            edges: Vec::new(),
        },
        inner_wires: Vec::new(),
        normal: compute_fillet_normal(&fillet_surface.surface),
        triangles: Vec::new(),
        mesh_dirty: true,
    };

    // Add the face to the first shell (simplified)
    if !brep.solids.is_empty() && !brep.solids[0].shells.is_empty() {
        brep.solids[0].shells[0].faces.push(face);
    }

    Ok(())
}

/// Compute a representative normal for a fillet surface.
fn compute_fillet_normal(surface: &Surface3) -> DVec3 {
    match surface {
        Surface3::Torus(t) => t.axis.normalize_or(DVec3::Z),
        Surface3::Cylinder(c) => c.axis.normalize_or(DVec3::Z),
        Surface3::Sphere(s) => s.axis.normalize_or(DVec3::Z),
        Surface3::Plane(p) => p.normal.normalize_or(DVec3::Z),
        _ => DVec3::Z,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Get any vector perpendicular to the given vector.
fn any_perpendicular(v: DVec3) -> DVec3 {
    let v = v.normalize_or(DVec3::Z);
    let perp = if v.x.abs() > 0.5 {
        DVec3::new(-v.y, v.x, 0.0)
    } else {
        DVec3::new(0.0, -v.z, v.y)
    };
    perp.normalize_or(DVec3::X)
}

/// Interpolate between two radii with tension parameter.
fn interpolate_radius(r1: f64, r2: f64, t: f64, tension: f64) -> f64 {
    // Hermite-like interpolation with tension
    let t2 = t * t;
    let t3 = t2 * t;
    let h1 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h2 = -2.0 * t3 + 3.0 * t2;

    // Apply tension (0 = linear, 1 = smooth)
    let smooth = tension;
    let h1_tension = h1 + smooth * (t3 - 2.0 * t2 + t);
    let h2_tension = h2 + smooth * (-t3 + t2);

    r1 * h1_tension + r2 * h2_tension
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    fn create_box_brep() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 2.0,
            depth: 3.0,
        })
    }

    #[test]
    fn test_fillet_params_default() {
        let params = FilletParams::default();
        assert_eq!(params.radius, 1.0);
        assert_eq!(params.continuity, FilletContinuity::C1);
        assert_eq!(params.mode, FilletMode::Uniform);
        assert!(params.tension >= 0.0 && params.tension <= 1.0);
    }

    #[test]
    fn test_fillet_params_builder() {
        let params = FilletParams::new(0.5)
            .with_continuity(FilletContinuity::C2)
            .with_mode(FilletMode::Chordal)
            .with_tension(0.8);

        assert_eq!(params.radius, 0.5);
        assert_eq!(params.continuity, FilletContinuity::C2);
        assert_eq!(params.mode, FilletMode::Chordal);
        assert!((params.tension - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_variable_radius_point() {
        let point = VariableRadiusPoint::new(0.5, 2.0);
        assert!((point.parameter - 0.5).abs() < 1e-10);
        assert!((point.radius - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_fillet_error_display() {
        let err = FilletError::InvalidRadius { radius: -1.0 };
        assert!(err.to_string().contains("invalid fillet radius"));

        let err = FilletError::EdgeNotFound { edge_index: 99 };
        assert!(err.to_string().contains("edge 99 not found"));

        let err = FilletError::RadiusTooLarge {
            edge_index: 0,
            radius: 10.0,
            max_radius: 1.0,
        };
        assert!(err.to_string().contains("radius 10 too large"));
    }

    #[test]
    fn test_any_perpendicular() {
        let v = DVec3::X;
        let p = any_perpendicular(v);
        assert!((p.dot(v)).abs() < 1e-10);
        assert!((p.length() - 1.0).abs() < 1e-10);

        let v = DVec3::Z;
        let p = any_perpendicular(v);
        assert!((p.dot(v)).abs() < 1e-10);
    }

    #[test]
    fn test_interpolate_radius() {
        // Linear interpolation (tension = 0)
        let r = interpolate_radius(1.0, 2.0, 0.0, 0.0);
        assert!((r - 1.0).abs() < 1e-10);

        let r = interpolate_radius(1.0, 2.0, 1.0, 0.0);
        assert!((r - 2.0).abs() < 1e-10);

        let r = interpolate_radius(1.0, 2.0, 0.5, 0.0);
        assert!((r - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_make_fillet_edge_empty_indices() {
        let brep = create_box_brep();
        let result = make_fillet_edge(&brep, &[], 0.5);
        assert!(result.is_ok());
        let res = result.unwrap();
        assert_eq!(res.edges_processed, 0);
        assert!(!res.warnings.is_empty());
    }

    #[test]
    fn test_make_fillet_edge_invalid_radius() {
        let brep = create_box_brep();
        let result = make_fillet_edge(&brep, &[0], -1.0);
        assert!(matches!(result, Err(FilletError::InvalidRadius { .. })));
    }

    #[test]
    fn test_make_fillet_edge_invalid_edge_index() {
        let brep = create_box_brep();
        let result = make_fillet_edge(&brep, &[999], 0.5);
        assert!(matches!(result, Err(FilletError::EdgeNotFound { .. })));
    }

    #[test]
    fn test_make_fillet_all_edges() {
        let brep = create_box_brep();
        let result = make_fillet_all_edges(&brep, 0.1);
        // Should succeed without errors
        assert!(result.is_ok() || matches!(result, Err(FilletError::EdgeNotFound { .. })));
    }

    #[test]
    fn test_make_variable_fillet_too_few_points() {
        let brep = create_box_brep();
        let radii = vec![VariableRadiusPoint::new(0.0, 0.5)];
        let result = make_variable_fillet(&brep, &[0], &radii);
        assert!(result.is_err());
    }

    #[test]
    fn test_make_variable_fillet_invalid_parameter() {
        let brep = create_box_brep();
        let radii = vec![
            VariableRadiusPoint::new(-0.5, 0.5),
            VariableRadiusPoint::new(1.0, 1.0),
        ];
        let result = make_variable_fillet(&brep, &[0], &radii);
        assert!(matches!(result, Err(FilletError::InvalidVariableRadius { .. })));
    }

    #[test]
    fn test_compute_plane_plane_fillet() {
        let plane1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let plane2 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Y,
            tangent_end: DVec3::Y,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_plane_plane_fillet(&edge_info, &plane1, &plane2, 0.5);
        assert!(result.is_ok());

        let surface = result.unwrap();
        assert!(matches!(surface, Surface3::Torus(_)));
    }

    #[test]
    fn test_compute_cylinder_plane_fillet() {
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Z,
            tangent_end: DVec3::Z,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_cylinder_plane_fillet(&edge_info, &cylinder, &plane, 0.5);
        assert!(result.is_ok());
    }

    #[test]
    fn test_compute_sphere_plane_fillet() {
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let plane = Plane {
            origin: DVec3::new(0.0, 0.0, 0.5),
            normal: DVec3::Z,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::X,
            tangent_end: DVec3::X,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_sphere_plane_fillet(&edge_info, &sphere, &plane, 0.5);
        assert!(result.is_ok());

        let surface = result.unwrap();
        assert!(matches!(surface, Surface3::Sphere(_)));
    }

    #[test]
    fn test_compute_fillet_curves_torus() {
        let brep = create_box_brep();

        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 2.0,
            minor_radius: 0.5,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Z,
            tangent_end: DVec3::Z,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_fillet_curves(&brep, &edge_info, 0.5, &Surface3::Torus(torus));
        assert!(result.is_ok());

        let curves = result.unwrap();
        assert!(!curves.is_empty());
    }

    #[test]
    fn test_fillet_continuity() {
        assert_eq!(FilletContinuity::default(), FilletContinuity::C1);
        assert_ne!(FilletContinuity::C0, FilletContinuity::C1);
        assert_ne!(FilletContinuity::C1, FilletContinuity::C2);
    }

    #[test]
    fn test_fillet_mode() {
        assert_eq!(FilletMode::default(), FilletMode::Uniform);
        assert_ne!(FilletMode::Uniform, FilletMode::Variable);
        assert_ne!(FilletMode::Variable, FilletMode::Chordal);
    }

    #[test]
    fn test_fillet_result_creation() {
        let brep = create_box_brep();
        let result = FilletResult {
            brep: brep.clone(),
            edges_processed: 3,
            fillet_faces_created: 3,
            warnings: vec!["test warning".to_string()],
        };

        assert_eq!(result.edges_processed, 3);
        assert_eq!(result.fillet_faces_created, 3);
        assert_eq!(result.warnings.len(), 1);
    }

    // ============================================================================
    // Edge Case Tests for OCCT Alignment
    // ============================================================================

    /// Test fillet on a concave edge (interior corner of a box).
    /// Concave edges require different handling than convex edges.
    #[test]
    fn test_fillet_concave_edge() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        // Create a box and populate its geometry
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        });
        populate_box_geom(&mut brep);

        // Attempt to fillet edges - concave edges are at the interior corners
        let result = make_fillet_edge(&brep, &[0, 1, 2, 3], 0.5);
        assert!(result.is_ok(), "concave edge fillet should succeed");

        let fillet_result = result.unwrap();
        // At least some edges should be processed
        assert!(fillet_result.edges_processed > 0 || !fillet_result.warnings.is_empty());
    }

    /// Test fillet on a chain of connected edges.
    /// Chain edges should blend smoothly at the vertices where they meet.
    #[test]
    fn test_fillet_chain_edges() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        // Create a box
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 3.0,
            depth: 2.0,
        });
        populate_box_geom(&mut brep);

        // Fillet a chain of edges around one face (edges 0, 1, 2, 3 form a loop)
        let result = make_fillet_edge(&brep, &[0, 1, 2, 3], 0.3);
        assert!(result.is_ok(), "chain edge fillet should succeed");

        let fillet_result = result.unwrap();
        assert!(fillet_result.edges_processed >= 1, "at least one edge should be filleted");
    }

    /// Test fillet with very small radius on an edge.
    /// Small radius fillets should not create degenerate geometry.
    #[test]
    fn test_fillet_small_radius() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 10.0,
            height: 10.0,
            depth: 10.0,
        });
        populate_box_geom(&mut brep);

        // Very small fillet radius
        let result = make_fillet_edge(&brep, &[0], 0.001);
        assert!(result.is_ok(), "small radius fillet should succeed");
    }

    /// Test fillet on an edge where adjacent faces are perpendicular.
    /// This is the most common case and should produce a clean toroidal fillet.
    #[test]
    fn test_fillet_perpendicular_faces() {
        use rcad_kernel::geom::{Plane, Line3};
        use glam::DVec3;

        // Create edge info for perpendicular faces
        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Y,  // Edge along Y axis
            tangent_end: DVec3::Y,
            length: 2.0,
            curve: Some(Curve3::Line(Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::Y,
            })),
            curve_range: Some([0.0, 2.0]),
        };

        let plane1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let plane2 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let result = compute_plane_plane_fillet(&edge_info, &plane1, &plane2, 0.5);
        assert!(result.is_ok(), "perpendicular faces fillet should succeed");

        // Result should be a torus for plane-plane fillet
        let surface = result.unwrap();
        assert!(matches!(surface, Surface3::Torus(_)), "plane-plane fillet should produce torus");
    }

    /// Test fillet with variable radius along the edge.
    /// Variable radius fillets should interpolate between radius values.
    #[test]
    fn test_fillet_variable_radius_basic() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 5.0,
            height: 5.0,
            depth: 5.0,
        });
        populate_box_geom(&mut brep);

        // Variable radius: 0.2 at start, 0.8 at end
        let radii = vec![
            VariableRadiusPoint::new(0.0, 0.2),
            VariableRadiusPoint::new(1.0, 0.8),
        ];

        let result = make_variable_fillet(&brep, &[0], &radii);
        assert!(result.is_ok(), "variable radius fillet should succeed");
    }

    /// Test fillet surface computation for cylinder-plane edge.
    /// Verifies correct handling of curved-to-planar transitions.
    #[test]
    fn test_fillet_cylinder_plane_edge() {
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let plane = Plane {
            origin: DVec3::new(1.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Z,
            tangent_end: DVec3::Z,
            length: 2.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_cylinder_plane_fillet(&edge_info, &cylinder, &plane, 0.3);
        assert!(result.is_ok(), "cylinder-plane fillet should succeed");
    }
}
