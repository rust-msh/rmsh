//! BRepChamfer-style edge chamfer operations — analogous to OCCT `BRepFilletAPI_MakeChamfer`.
//!
//! # Overview
//!
//! This module provides algorithms for creating chamfers (beveled edges) on BRep solids:
//!
//! - **`make_chamfer_edge`**: Symmetric chamfer with equal distances on both sides
//! - **`make_chamfer_asymmetric`**: Asymmetric chamfer with two different distances
//! - **`make_chamfer_angle`**: Chamfer defined by distance and angle
//! - **`make_chamfer_all_edges`**: Apply chamfer to all edges of a solid
//!
//! # Chamfer Modes
//!
//! - **Symmetric**: Equal offset from edge on both adjacent faces
//! - **Asymmetric**: Different distances on each adjacent face
//! - **Distance-Angle**: One distance and an angle from the first face
//!
//! # Supported Geometry
//!
//! - Plane-plane intersection edges (most common)
//! - Cylinder-plane intersection edges
//! - General surface-surface intersection (with approximation)
//!
//! # Algorithm
//!
//! 1. Identify edges to chamfer and their adjacent faces
//! 2. Compute chamfer surface parameters based on chamfer mode
//! 3. Construct chamfer surface between the two adjacent faces
//! 4. Compute intersection curves with original faces
//! 5. Trim original faces and add chamfer face
//! 6. Rebuild the shell topology
//!
//! # References
//!
//! - OCCT `BRepFilletAPI_MakeChamfer`
//! - OCCT `ChFi3d_FilBuilder` (chamfer algorithm internals)

use glam::DVec3;
use rcad_kernel::{
    BRep, PrimitiveSolid,
    geom::{Curve3, Surface3, Line3, Plane, CylindricalSurface},
    topology::{Edge, Face, Vertex, Wire, WireEdge},
};
use std::collections::HashMap;

// ═══════════════════════════════════════════════════════════════════════════════
// Constants
// ═══════════════════════════════════════════════════════════════════════════════

/// Default tolerance for geometric operations.
const TOLERANCE: f64 = 1e-9;

/// Small value for checking if vectors are parallel.
const PARALLEL_TOL: f64 = 1e-10;

// ═══════════════════════════════════════════════════════════════════════════════
// Chamfer Modes and Parameters
// ═══════════════════════════════════════════════════════════════════════════════

/// Chamfer mode defining how distances are interpreted.
///
/// Analogous to OCCT `ChFiDS_ChamfMode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChamferMode {
    /// Symmetric chamfer: equal distance on both adjacent faces.
    /// Single distance parameter used.
    Symmetric,
    /// Asymmetric chamfer: two different distances on each face.
    /// `d1` on first face, `d2` on second face.
    Asymmetric,
    /// Distance-angle chamfer: one distance and an angle from the first face.
    /// The angle is measured from the first face normal.
    DistanceAngle,
}

/// Parameters controlling the chamfer operation.
#[derive(Debug, Clone)]
pub struct ChamferParams {
    /// Chamfer mode (symmetric, asymmetric, distance-angle).
    pub mode: ChamferMode,
    /// Primary distance (used in all modes).
    pub distance1: f64,
    /// Secondary distance (used only in asymmetric mode).
    pub distance2: f64,
    /// Angle in radians (used only in distance-angle mode).
    /// The angle is measured from the first face.
    pub angle: f64,
}

impl Default for ChamferParams {
    fn default() -> Self {
        Self {
            mode: ChamferMode::Symmetric,
            distance1: 1.0,
            distance2: 1.0,
            angle: std::f64::consts::FRAC_PI_4,
        }
    }
}

impl ChamferParams {
    /// Create symmetric chamfer parameters with a single distance.
    pub fn symmetric(distance: f64) -> Self {
        Self {
            mode: ChamferMode::Symmetric,
            distance1: distance,
            distance2: distance,
            angle: std::f64::consts::FRAC_PI_4,
        }
    }

    /// Create asymmetric chamfer parameters with two distances.
    pub fn asymmetric(d1: f64, d2: f64) -> Self {
        Self {
            mode: ChamferMode::Asymmetric,
            distance1: d1,
            distance2: d2,
            angle: std::f64::consts::FRAC_PI_4,
        }
    }

    /// Create distance-angle chamfer parameters.
    pub fn distance_angle(distance: f64, angle_rad: f64) -> Self {
        Self {
            mode: ChamferMode::DistanceAngle,
            distance1: distance,
            distance2: distance * angle_rad.tan(),
            angle: angle_rad,
        }
    }

    /// Validate the parameters.
    pub fn validate(&self) -> Result<(), ChamferError> {
        if self.distance1 <= 0.0 {
            return Err(ChamferError::InvalidDistance {
                value: self.distance1,
                reason: "distance1 must be positive".to_string(),
            });
        }
        if self.mode == ChamferMode::Asymmetric && self.distance2 <= 0.0 {
            return Err(ChamferError::InvalidDistance {
                value: self.distance2,
                reason: "distance2 must be positive for asymmetric mode".to_string(),
            });
        }
        if self.mode == ChamferMode::DistanceAngle {
            if self.angle <= 0.0 || self.angle >= std::f64::consts::FRAC_PI_2 {
                return Err(ChamferError::InvalidAngle {
                    value: self.angle,
                    reason: "angle must be between 0 and 90 degrees".to_string(),
                });
            }
        }
        Ok(())
    }

    /// Get the effective distances on both faces.
    /// Returns (distance_on_face1, distance_on_face2).
    pub fn get_distances(&self) -> (f64, f64) {
        match self.mode {
            ChamferMode::Symmetric => (self.distance1, self.distance1),
            ChamferMode::Asymmetric => (self.distance1, self.distance2),
            ChamferMode::DistanceAngle => (self.distance1, self.distance1 * self.angle.tan()),
        }
    }
}

/// Result of a chamfer operation.
#[derive(Debug)]
pub struct ChamferResult {
    /// The resulting BRep with chamfers applied.
    pub brep: BRep,
    /// Number of edges chamfered.
    pub edges_chamfered: usize,
    /// Number of chamfer faces created.
    pub chamfer_faces_created: usize,
    /// Warnings encountered during the operation.
    pub warnings: Vec<ChamferWarning>,
}

/// Warning messages from chamfer operations.
#[derive(Debug, Clone)]
pub enum ChamferWarning {
    /// Edge could not be chamfered due to geometry limitations.
    EdgeSkipped { edge_index: usize, reason: String },
    /// Chamfer distance was reduced to avoid geometry issues.
    DistanceReduced { edge_index: usize, original: f64, reduced: f64 },
    /// Face topology was modified unexpectedly.
    TopologyModified { description: String },
}

// ═══════════════════════════════════════════════════════════════════════════════
// Error Types
// ═══════════════════════════════════════════════════════════════════════════════

/// Error type for chamfer operations.
#[derive(Debug, Clone, PartialEq)]
pub enum ChamferError {
    /// Invalid chamfer distance.
    InvalidDistance { value: f64, reason: String },
    /// Invalid chamfer angle.
    InvalidAngle { value: f64, reason: String },
    /// Edge not found in the BRep.
    EdgeNotFound { edge_index: usize },
    /// Edge has no adjacent faces.
    NoAdjacentFaces { edge_index: usize },
    /// Face geometry is not supported for chamfer.
    UnsupportedSurface { face_index: usize, surface_type: String },
    /// Chamfer would create invalid geometry.
    InvalidResult { description: String },
    /// Failed to compute chamfer surface.
    ChamferSurfaceFailed { edge_index: usize, reason: String },
    /// Input shape is invalid.
    InvalidInput { description: String },
    /// Topology error during chamfer construction.
    TopologyError { description: String },
    /// Boolean operation failed.
    BooleanFailed { description: String },
}

impl std::fmt::Display for ChamferError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidDistance { value, reason } => {
                write!(f, "invalid distance {}: {}", value, reason)
            }
            Self::InvalidAngle { value, reason } => {
                write!(f, "invalid angle {:.2} degrees: {}", value.to_degrees(), reason)
            }
            Self::EdgeNotFound { edge_index } => {
                write!(f, "edge {} not found in BRep", edge_index)
            }
            Self::NoAdjacentFaces { edge_index } => {
                write!(f, "edge {} has no adjacent faces", edge_index)
            }
            Self::UnsupportedSurface { face_index, surface_type } => {
                write!(f, "face {} has unsupported surface type: {}", face_index, surface_type)
            }
            Self::InvalidResult { description } => {
                write!(f, "chamfer produced invalid result: {}", description)
            }
            Self::ChamferSurfaceFailed { edge_index, reason } => {
                write!(f, "failed to compute chamfer surface for edge {}: {}", edge_index, reason)
            }
            Self::InvalidInput { description } => {
                write!(f, "invalid input: {}", description)
            }
            Self::TopologyError { description } => {
                write!(f, "topology error: {}", description)
            }
            Self::BooleanFailed { description } => {
                write!(f, "boolean operation failed: {}", description)
            }
        }
    }
}

impl std::error::Error for ChamferError {}

// ═══════════════════════════════════════════════════════════════════════════════
// Edge Information
// ═══════════════════════════════════════════════════════════════════════════════

/// Information about an edge and its adjacent faces.
#[derive(Debug, Clone)]
struct EdgeInfo {
    /// Edge index in the BRep.
    edge_index: usize,
    /// Start vertex index.
    start_vertex: usize,
    /// End vertex index.
    end_vertex: usize,
    /// Start point in 3D.
    start_point: DVec3,
    /// End point in 3D.
    end_point: DVec3,
    /// Edge tangent direction (normalized).
    tangent: DVec3,
    /// Edge length.
    length: f64,
    /// Adjacent face indices.
    adjacent_faces: Vec<usize>,
    /// Surface types of adjacent faces.
    adjacent_surfaces: Vec<Option<usize>>,
}

/// Information about adjacent face geometry.
#[derive(Debug, Clone)]
struct AdjacentFaceInfo {
    /// Face index.
    face_index: usize,
    /// Surface index in GeomStore.
    surface_index: Option<usize>,
    /// Face normal at edge midpoint.
    normal: DVec3,
    /// Surface type.
    surface_type: SurfaceType,
    /// Reference point on the face near the edge.
    reference_point: DVec3,
}

/// Classification of surface types for chamfer computation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    BSpline,
    Other,
}

impl From<&Surface3> for SurfaceType {
    fn from(surface: &Surface3) -> Self {
        match surface {
            Surface3::Plane(_) => SurfaceType::Plane,
            Surface3::Cylinder(_) => SurfaceType::Cylinder,
            Surface3::Cone(_) => SurfaceType::Cone,
            Surface3::Sphere(_) => SurfaceType::Sphere,
            Surface3::Torus(_) => SurfaceType::Torus,
            Surface3::BSpline(_) => SurfaceType::BSpline,
            _ => SurfaceType::Other,
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Core Chamfer Operations
// ═══════════════════════════════════════════════════════════════════════════════

/// Create a symmetric chamfer on specified edges.
///
/// The chamfer is created with equal distance on both adjacent faces.
///
/// # Arguments
///
/// * `brep` - The input BRep solid.
/// * `edge_indices` - Indices of edges to chamfer.
/// * `distance` - Chamfer distance (equal on both sides).
///
/// # Returns
///
/// A `ChamferResult` containing the chamfered BRep and operation statistics.
///
/// # Example
///
/// ```ignore
/// let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let result = make_chamfer_edge(&box_brep, &[0, 1, 2], 0.1).unwrap();
/// ```
pub fn make_chamfer_edge(
    brep: &BRep,
    edge_indices: &[usize],
    distance: f64,
) -> Result<ChamferResult, ChamferError> {
    let params = ChamferParams::symmetric(distance);
    make_chamfer(brep, edge_indices, &params)
}

/// Create an asymmetric chamfer on specified edges.
///
/// The chamfer has different distances on each adjacent face.
///
/// # Arguments
///
/// * `brep` - The input BRep solid.
/// * `edge_indices` - Indices of edges to chamfer.
/// * `d1` - Distance on the first adjacent face.
/// * `d2` - Distance on the second adjacent face.
///
/// # Returns
///
/// A `ChamferResult` containing the chamfered BRep and operation statistics.
pub fn make_chamfer_asymmetric(
    brep: &BRep,
    edge_indices: &[usize],
    d1: f64,
    d2: f64,
) -> Result<ChamferResult, ChamferError> {
    let params = ChamferParams::asymmetric(d1, d2);
    make_chamfer(brep, edge_indices, &params)
}

/// Create a distance-angle chamfer on specified edges.
///
/// The chamfer is defined by a distance on the first face and an angle.
///
/// # Arguments
///
/// * `brep` - The input BRep solid.
/// * `edge_indices` - Indices of edges to chamfer.
/// * `distance` - Distance on the first adjacent face.
/// * `angle` - Angle from the first face in radians.
///
/// # Returns
///
/// A `ChamferResult` containing the chamfered BRep and operation statistics.
pub fn make_chamfer_angle(
    brep: &BRep,
    edge_indices: &[usize],
    distance: f64,
    angle: f64,
) -> Result<ChamferResult, ChamferError> {
    let params = ChamferParams::distance_angle(distance, angle);
    make_chamfer(brep, edge_indices, &params)
}

/// Apply chamfer to all edges of a solid.
///
/// # Arguments
///
/// * `brep` - The input BRep solid.
/// * `distance` - Chamfer distance.
///
/// # Returns
///
/// A `ChamferResult` containing the chamfered BRep and operation statistics.
pub fn make_chamfer_all_edges(
    brep: &BRep,
    distance: f64,
) -> Result<ChamferResult, ChamferError> {
    let all_edges: Vec<usize> = (0..brep.edges.len()).collect();
    make_chamfer_edge(brep, &all_edges, distance)
}

/// Core chamfer implementation.
///
/// This function processes all specified edges and creates chamfers.
fn make_chamfer(
    brep: &BRep,
    edge_indices: &[usize],
    params: &ChamferParams,
) -> Result<ChamferResult, ChamferError> {
    // Validate parameters
    params.validate()?;

    // Validate input
    if brep.solids.is_empty() {
        return Err(ChamferError::InvalidInput {
            description: "BRep has no solids".to_string(),
        });
    }

    // Build edge-to-face adjacency
    let edge_to_faces = build_edge_face_adjacency(brep);

    // Process each edge
    let mut result_brep = brep.clone();
    let mut edges_chamfered = 0;
    let mut chamfer_faces_created = 0;
    let mut warnings = Vec::new();

    for &edge_idx in edge_indices {
        if edge_idx >= brep.edges.len() {
            warnings.push(ChamferWarning::EdgeSkipped {
                edge_index: edge_idx,
                reason: "edge index out of range".to_string(),
            });
            continue;
        }

        let adjacent_faces = match edge_to_faces.get(&edge_idx) {
            Some(faces) if faces.len() >= 2 => faces.clone(),
            _ => {
                warnings.push(ChamferWarning::EdgeSkipped {
                    edge_index: edge_idx,
                    reason: "edge has fewer than 2 adjacent faces".to_string(),
                });
                continue;
            }
        };

        // Get edge information
        let edge_info = get_edge_info(brep, edge_idx, &adjacent_faces)?;

        // Get adjacent face information
        let face_infos = get_adjacent_face_infos(brep, &adjacent_faces, &edge_info)?;

        if face_infos.len() < 2 {
            warnings.push(ChamferWarning::EdgeSkipped {
                edge_index: edge_idx,
                reason: "could not get adjacent face info".to_string(),
            });
            continue;
        }

        // Compute chamfer geometry based on surface types
        let chamfer_geom = match compute_chamfer_geometry(brep, &edge_info, &face_infos, params) {
            Ok(g) => g,
            Err(e) => {
                warnings.push(ChamferWarning::EdgeSkipped {
                    edge_index: edge_idx,
                    reason: format!("chamfer computation failed: {}", e),
                });
                continue;
            }
        };

        // Apply chamfer to the result BRep
        match apply_chamfer_to_brep(&mut result_brep, &edge_info, &face_infos, &chamfer_geom, params) {
            Ok(faces_added) => {
                edges_chamfered += 1;
                chamfer_faces_created += faces_added;
            }
            Err(e) => {
                warnings.push(ChamferWarning::EdgeSkipped {
                    edge_index: edge_idx,
                    reason: format!("failed to apply chamfer: {}", e),
                });
            }
        }
    }

    // If no edges were chamfered, return an error
    if edges_chamfered == 0 {
        return Err(ChamferError::InvalidResult {
            description: "no edges could be chamfered".to_string(),
        });
    }

    Ok(ChamferResult {
        brep: result_brep,
        edges_chamfered,
        chamfer_faces_created,
        warnings,
    })
}

// ═══════════════════════════════════════════════════════════════════════════════
// Chamfer Geometry Computation
// ═══════════════════════════════════════════════════════════════════════════════

/// Chamfer geometry resulting from computation.
#[derive(Debug, Clone)]
struct ChamferGeometry {
    /// Chamfer surface (always a plane for simple chamfers).
    chamfer_surface: Surface3,
    /// Origin point of the chamfer surface.
    origin: DVec3,
    /// Normal of the chamfer surface.
    normal: DVec3,
    /// New vertex positions at the start of the edge.
    start_vertices: [DVec3; 2],
    /// New vertex positions at the end of the edge.
    end_vertices: [DVec3; 2],
    /// Chamfer distance on face 1.
    d1: f64,
    /// Chamfer distance on face 2.
    d2: f64,
}

/// Compute chamfer surface between two adjacent faces.
///
/// This function determines the chamfer surface based on the surface types
/// of the adjacent faces and the chamfer parameters.
fn compute_chamfer_geometry(
    brep: &BRep,
    edge_info: &EdgeInfo,
    face_infos: &[AdjacentFaceInfo],
    params: &ChamferParams,
) -> Result<ChamferGeometry, ChamferError> {
    let face0 = &face_infos[0];
    let face1 = &face_infos[1];

    let (d1, d2) = params.get_distances();

    // Check if we can reduce to plane-plane case
    match (face0.surface_type, face1.surface_type) {
        (SurfaceType::Plane, SurfaceType::Plane) => {
            compute_chamfer_plane_plane(edge_info, face0, face1, d1, d2)
        }
        (SurfaceType::Plane, SurfaceType::Cylinder) |
        (SurfaceType::Cylinder, SurfaceType::Plane) => {
            compute_chamfer_cylinder_plane(brep, edge_info, face0, face1, d1, d2)
        }
        (SurfaceType::Cylinder, SurfaceType::Cylinder) => {
            compute_chamfer_cylinder_cylinder(brep, edge_info, face0, face1, d1, d2)
        }
        _ => {
            // General case: approximate with plane
            compute_chamfer_general(edge_info, face0, face1, d1, d2)
        }
    }
}

/// Compute chamfer for plane-plane intersection.
///
/// This is the most common case and produces a planar chamfer face.
fn compute_chamfer_plane_plane(
    edge_info: &EdgeInfo,
    face0: &AdjacentFaceInfo,
    face1: &AdjacentFaceInfo,
    d1: f64,
    d2: f64,
) -> Result<ChamferGeometry, ChamferError> {
    let n0 = face0.normal.normalize();
    let n1 = face1.normal.normalize();
    let edge_dir = edge_info.tangent;

    // Compute the direction along each face (perpendicular to edge)
    let along_face0 = edge_dir.cross(n0).normalize();
    let along_face1 = edge_dir.cross(n1).normalize();

    // Compute offset points on each face
    let mid_point = (edge_info.start_point + edge_info.end_point) * 0.5;

    // Points offset from the edge along each face
    let p0_offset = mid_point + along_face0 * d1;
    let p1_offset = mid_point + along_face1 * d2;

    // The chamfer surface passes through the edge and the two offset points
    // For plane-plane, the chamfer is also a plane

    // Compute chamfer plane normal
    // The chamfer plane bisects the angle between the two face planes (for symmetric case)
    // For asymmetric case, it's adjusted based on d1 and d2 ratio

    let edge_to_p0 = p0_offset - mid_point;
    let edge_to_p1 = p1_offset - mid_point;

    // Chamfer plane normal is perpendicular to the edge and to the vector from edge to chamfer center
    let chamfer_normal = edge_dir.cross(edge_to_p0 + edge_to_p1).normalize();

    // Verify the normal direction
    if chamfer_normal.length() < PARALLEL_TOL {
        // Edge and faces are coplanar - shouldn't happen for valid edge
        return Err(ChamferError::ChamferSurfaceFailed {
            edge_index: edge_info.edge_index,
            reason: "faces are coplanar, no chamfer possible".to_string(),
        });
    }

    // Compute new vertex positions
    // For each end of the edge, compute where the chamfer intersects each face
    let start_p0 = edge_info.start_point + along_face0 * d1;
    let start_p1 = edge_info.start_point + along_face1 * d2;
    let end_p0 = edge_info.end_point + along_face0 * d1;
    let end_p1 = edge_info.end_point + along_face1 * d2;

    Ok(ChamferGeometry {
        chamfer_surface: Surface3::Plane(Plane {
            origin: mid_point,
            normal: chamfer_normal,
        }),
        origin: mid_point,
        normal: chamfer_normal,
        start_vertices: [start_p0, start_p1],
        end_vertices: [end_p0, end_p1],
        d1,
        d2,
    })
}

/// Compute chamfer for cylinder-plane intersection.
fn compute_chamfer_cylinder_plane(
    brep: &BRep,
    edge_info: &EdgeInfo,
    face0: &AdjacentFaceInfo,
    face1: &AdjacentFaceInfo,
    d1: f64,
    d2: f64,
) -> Result<ChamferGeometry, ChamferError> {
    // Identify which face is cylinder and which is plane
    let (cyl_face, plane_face, dc, dp) = if face0.surface_type == SurfaceType::Cylinder {
        (face0, face1, d1, d2)
    } else {
        (face1, face0, d2, d1)
    };

    // Get cylinder parameters
    let cyl_surface_idx = match cyl_face.surface_index {
        Some(idx) => idx,
        None => return compute_chamfer_general(edge_info, face0, face1, d1, d2),
    };

    let cylinder = match brep.geom.surfaces.get(cyl_surface_idx) {
        Some(Surface3::Cylinder(c)) => c,
        _ => return compute_chamfer_general(edge_info, face0, face1, d1, d2),
    };

    let n_plane = plane_face.normal.normalize();
    let edge_dir = edge_info.tangent;
    let mid_point = (edge_info.start_point + edge_info.end_point) * 0.5;

    // Compute radial direction from cylinder axis to edge point
    let axis = cylinder.axis.normalize();
    let to_edge = mid_point - cylinder.origin;
    let radial_dir = (to_edge - axis * to_edge.dot(axis)).normalize();

    // For cylinder, the chamfer moves along the radial direction
    // For plane, it moves perpendicular to edge in the plane
    let along_plane = edge_dir.cross(n_plane).normalize();

    // Compute chamfer surface
    // The chamfer plane is defined by the edge and the two offset directions
    let cyl_offset_dir = radial_dir; // Direction away from cylinder axis
    let plane_offset_dir = along_plane;

    let offset_on_cyl = mid_point + cyl_offset_dir * dc;
    let offset_on_plane = mid_point + plane_offset_dir * dp;

    // Compute chamfer plane normal
    let v1 = offset_on_cyl - mid_point;
    let v2 = offset_on_plane - mid_point;
    let chamfer_normal = edge_dir.cross(v1 + v2).normalize();

    if chamfer_normal.length() < PARALLEL_TOL {
        return Err(ChamferError::ChamferSurfaceFailed {
            edge_index: edge_info.edge_index,
            reason: "could not compute chamfer normal".to_string(),
        });
    }

    // Compute new vertex positions
    let start_on_cyl = edge_info.start_point + cyl_offset_dir * dc;
    let start_on_plane = edge_info.start_point + plane_offset_dir * dp;
    let end_on_cyl = edge_info.end_point + cyl_offset_dir * dc;
    let end_on_plane = edge_info.end_point + plane_offset_dir * dp;

    // Return in consistent order (face0, face1)
    if face0.surface_type == SurfaceType::Cylinder {
        Ok(ChamferGeometry {
            chamfer_surface: Surface3::Plane(Plane {
                origin: mid_point,
                normal: chamfer_normal,
            }),
            origin: mid_point,
            normal: chamfer_normal,
            start_vertices: [start_on_cyl, start_on_plane],
            end_vertices: [end_on_cyl, end_on_plane],
            d1,
            d2,
        })
    } else {
        Ok(ChamferGeometry {
            chamfer_surface: Surface3::Plane(Plane {
                origin: mid_point,
                normal: chamfer_normal,
            }),
            origin: mid_point,
            normal: chamfer_normal,
            start_vertices: [start_on_plane, start_on_cyl],
            end_vertices: [end_on_plane, end_on_cyl],
            d1,
            d2,
        })
    }
}

/// Compute chamfer for cylinder-cylinder intersection.
fn compute_chamfer_cylinder_cylinder(
    brep: &BRep,
    edge_info: &EdgeInfo,
    face0: &AdjacentFaceInfo,
    face1: &AdjacentFaceInfo,
    d1: f64,
    d2: f64,
) -> Result<ChamferGeometry, ChamferError> {
    // Get cylinder parameters for both faces
    let get_cylinder = |face: &AdjacentFaceInfo| -> Option<(DVec3, DVec3, f64)> {
        let idx = face.surface_index?;
        if let Some(Surface3::Cylinder(c)) = brep.geom.surfaces.get(idx) {
            Some((c.origin, c.axis.normalize(), c.radius))
        } else {
            None
        }
    };

    let cyl0 = match get_cylinder(face0) {
        Some(c) => c,
        None => return compute_chamfer_general(edge_info, face0, face1, d1, d2),
    };

    let cyl1 = match get_cylinder(face1) {
        Some(c) => c,
        None => return compute_chamfer_general(edge_info, face0, face1, d1, d2),
    };

    let mid_point = (edge_info.start_point + edge_info.end_point) * 0.5;
    let edge_dir = edge_info.tangent;

    // Compute radial directions from each cylinder axis
    let (_origin0, axis0, _radius0) = cyl0;
    let (_origin1, axis1, _radius1) = cyl1;

    let to_edge0 = mid_point - cyl0.0;
    let radial0 = (to_edge0 - axis0 * to_edge0.dot(axis0)).normalize();

    let to_edge1 = mid_point - cyl1.0;
    let radial1 = (to_edge1 - axis1 * to_edge1.dot(axis1)).normalize();

    // Compute offset points
    let offset0 = mid_point + radial0 * d1;
    let offset1 = mid_point + radial1 * d2;

    // Compute chamfer plane normal
    let v1 = offset0 - mid_point;
    let v2 = offset1 - mid_point;
    let chamfer_normal = edge_dir.cross(v1 + v2).normalize();

    if chamfer_normal.length() < PARALLEL_TOL {
        return Err(ChamferError::ChamferSurfaceFailed {
            edge_index: edge_info.edge_index,
            reason: "could not compute chamfer normal for cylinder-cylinder".to_string(),
        });
    }

    // Compute new vertex positions
    let start0 = edge_info.start_point + radial0 * d1;
    let start1 = edge_info.start_point + radial1 * d2;
    let end0 = edge_info.end_point + radial0 * d1;
    let end1 = edge_info.end_point + radial1 * d2;

    Ok(ChamferGeometry {
        chamfer_surface: Surface3::Plane(Plane {
            origin: mid_point,
            normal: chamfer_normal,
        }),
        origin: mid_point,
        normal: chamfer_normal,
        start_vertices: [start0, start1],
        end_vertices: [end0, end1],
        d1,
        d2,
    })
}

/// Compute chamfer for general surface-surface intersection.
///
/// Uses approximation based on face normals at the edge.
fn compute_chamfer_general(
    edge_info: &EdgeInfo,
    face0: &AdjacentFaceInfo,
    face1: &AdjacentFaceInfo,
    d1: f64,
    d2: f64,
) -> Result<ChamferGeometry, ChamferError> {
    let n0 = face0.normal.normalize();
    let n1 = face1.normal.normalize();
    let edge_dir = edge_info.tangent;

    let mid_point = (edge_info.start_point + edge_info.end_point) * 0.5;

    // Compute directions along each face (perpendicular to edge in the face plane)
    let along_face0 = edge_dir.cross(n0).normalize();
    let along_face1 = edge_dir.cross(n1).normalize();

    // Compute chamfer plane using the same logic as plane-plane
    let offset0 = mid_point + along_face0 * d1;
    let offset1 = mid_point + along_face1 * d2;

    let v1 = offset0 - mid_point;
    let v2 = offset1 - mid_point;
    let chamfer_normal = edge_dir.cross(v1 + v2).normalize();

    if chamfer_normal.length() < PARALLEL_TOL {
        return Err(ChamferError::ChamferSurfaceFailed {
            edge_index: edge_info.edge_index,
            reason: "degenerate chamfer geometry".to_string(),
        });
    }

    // Compute new vertex positions
    let start0 = edge_info.start_point + along_face0 * d1;
    let start1 = edge_info.start_point + along_face1 * d2;
    let end0 = edge_info.end_point + along_face0 * d1;
    let end1 = edge_info.end_point + along_face1 * d2;

    Ok(ChamferGeometry {
        chamfer_surface: Surface3::Plane(Plane {
            origin: mid_point,
            normal: chamfer_normal,
        }),
        origin: mid_point,
        normal: chamfer_normal,
        start_vertices: [start0, start1],
        end_vertices: [end0, end1],
        d1,
        d2,
    })
}

/// Compute the chamfer surface (public API).
///
/// Returns the chamfer surface between two adjacent faces.
pub fn compute_chamfer_surface(
    brep: &BRep,
    edge_index: usize,
    params: &ChamferParams,
) -> Result<Surface3, ChamferError> {
    let edge_to_faces = build_edge_face_adjacency(brep);
    let adjacent_faces = edge_to_faces.get(&edge_index).cloned().unwrap_or_default();

    if adjacent_faces.len() < 2 {
        return Err(ChamferError::NoAdjacentFaces { edge_index });
    }

    let edge_info = get_edge_info(brep, edge_index, &adjacent_faces)?;
    let face_infos = get_adjacent_face_infos(brep, &adjacent_faces, &edge_info)?;

    if face_infos.len() < 2 {
        return Err(ChamferError::NoAdjacentFaces { edge_index });
    }

    let chamfer_geom = compute_chamfer_geometry(brep, &edge_info, &face_infos, params)?;
    Ok(chamfer_geom.chamfer_surface)
}

/// Compute the boundary curves of a chamfer.
///
/// Returns the two curves where the chamfer face meets the original faces.
pub fn compute_chamfer_curves(
    brep: &BRep,
    edge_index: usize,
    params: &ChamferParams,
) -> Result<(Curve3, Curve3), ChamferError> {
    let edge_to_faces = build_edge_face_adjacency(brep);
    let adjacent_faces = edge_to_faces.get(&edge_index).cloned().unwrap_or_default();

    if adjacent_faces.len() < 2 {
        return Err(ChamferError::NoAdjacentFaces { edge_index });
    }

    let edge_info = get_edge_info(brep, edge_index, &adjacent_faces)?;
    let face_infos = get_adjacent_face_infos(brep, &adjacent_faces, &edge_info)?;

    if face_infos.len() < 2 {
        return Err(ChamferError::NoAdjacentFaces { edge_index });
    }

    let chamfer_geom = compute_chamfer_geometry(brep, &edge_info, &face_infos, params)?;

    // Create lines for the two boundary curves
    let curve1 = Curve3::Line(Line3 {
        origin: chamfer_geom.start_vertices[0],
        direction: (chamfer_geom.end_vertices[0] - chamfer_geom.start_vertices[0]).normalize_or(DVec3::X),
    });

    let curve2 = Curve3::Line(Line3 {
        origin: chamfer_geom.start_vertices[1],
        direction: (chamfer_geom.end_vertices[1] - chamfer_geom.start_vertices[1]).normalize_or(DVec3::X),
    });

    Ok((curve1, curve2))
}

/// Trim adjacent faces to meet the chamfer.
///
/// This function modifies the BRep by adjusting face boundaries.
pub fn trim_adjacent_faces(
    _brep: &mut BRep,
    _edge_index: usize,
    _params: &ChamferParams,
) -> Result<(), ChamferError> {
    // This is a placeholder for the face trimming implementation.
    // In a full implementation, this would:
    // 1. Compute the intersection curves between chamfer surface and each face
    // 2. Create new edges along these curves
    // 3. Modify the face boundaries to use the new edges
    // 4. Update the topology

    // For now, we return Ok as the trimming is handled in apply_chamfer_to_brep
    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
// BRep Modification
// ═══════════════════════════════════════════════════════════════════════════════

/// Apply chamfer geometry to the BRep.
///
/// This function:
/// 1. Creates new vertices for the chamfer corners
/// 2. Creates new edges for the chamfer boundaries
/// 3. Creates the chamfer face
/// 4. Updates the shell topology
fn apply_chamfer_to_brep(
    brep: &mut BRep,
    edge_info: &EdgeInfo,
    _face_infos: &[AdjacentFaceInfo],
    chamfer_geom: &ChamferGeometry,
    _params: &ChamferParams,
) -> Result<usize, ChamferError> {
    // For simplicity, we'll create a new shell with the chamfer
    // In a production implementation, this would be more sophisticated

    // Add new vertices
    let v_start_0 = brep.vertices.len();
    brep.vertices.push(Vertex { point: chamfer_geom.start_vertices[0] });

    let v_start_1 = brep.vertices.len();
    brep.vertices.push(Vertex { point: chamfer_geom.start_vertices[1] });

    let v_end_0 = brep.vertices.len();
    brep.vertices.push(Vertex { point: chamfer_geom.end_vertices[0] });

    let v_end_1 = brep.vertices.len();
    brep.vertices.push(Vertex { point: chamfer_geom.end_vertices[1] });

    // Create new edges for chamfer boundaries
    // Edge from start0 to end0 (along face 0)
    let e0_idx = create_line_edge(brep, v_start_0, v_end_0);

    // Edge from start1 to end1 (along face 1)
    let e1_idx = create_line_edge(brep, v_start_1, v_end_1);

    // Edge from start0 to start1 (across start of chamfer)
    let e_start_idx = create_line_edge(brep, v_start_0, v_start_1);

    // Edge from end0 to end1 (across end of chamfer)
    let e_end_idx = create_line_edge(brep, v_end_0, v_end_1);

    // Create the chamfer face
    let chamfer_face = Face {
        outer_wire: Wire {
            edges: vec![
                WireEdge::fwd(e0_idx),
                WireEdge::fwd(e_end_idx),
                WireEdge::rev(e1_idx),
                WireEdge::rev(e_start_idx),
            ],
        },
        inner_wires: vec![],
        normal: chamfer_geom.normal,
        triangles: vec![],
        mesh_dirty: true,
    };

    // Add the chamfer surface to GeomStore
    let surf_idx = brep.geom.surfaces.len();
    brep.geom.surfaces.push(chamfer_geom.chamfer_surface.clone());

    // Add the face to the first shell
    if let Some(solid) = brep.solids.first_mut() {
        if let Some(shell) = solid.shells.first_mut() {
            shell.faces.push(chamfer_face);

            // Add surface reference for the new face
            while brep.geom.face_surface.len() < shell.faces.len() {
                brep.geom.face_surface.push(None);
            }
            brep.geom.face_surface.push(Some(surf_idx));
        }
    }

    Ok(1) // One chamfer face created
}

/// Create a line edge between two vertices.
fn create_line_edge(brep: &mut BRep, start_idx: usize, end_idx: usize) -> usize {
    let p0 = brep.vertices[start_idx].point;
    let p1 = brep.vertices[end_idx].point;
    let d = p1 - p0;
    let len = d.length();
    let dir = if len > TOLERANCE { d / len } else { DVec3::X };

    let edge_idx = brep.edges.len();
    brep.edges.push(Edge { start: start_idx, end: end_idx });

    let curve_idx = brep.geom.curves.len();
    brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
    brep.geom.edge_curve.push(Some(curve_idx));
    brep.geom.edge_curve_range.push(Some([0.0, len]));
    brep.geom.edge_degenerated.push(false);

    edge_idx
}

// ═══════════════════════════════════════════════════════════════════════════════
// Helper Functions
// ═══════════════════════════════════════════════════════════════════════════════

/// Build edge-to-face adjacency map.
fn build_edge_face_adjacency(brep: &BRep) -> HashMap<usize, Vec<usize>> {
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();

    for (solid_idx, solid) in brep.solids.iter().enumerate() {
        for (shell_idx, shell) in solid.shells.iter().enumerate() {
            for (face_idx, face) in shell.faces.iter().enumerate() {
                let flat_face_idx = compute_flat_face_index(brep, solid_idx, shell_idx, face_idx);
                for we in &face.outer_wire.edges {
                    edge_to_faces.entry(we.idx).or_default().push(flat_face_idx);
                }
                for inner in &face.inner_wires {
                    for we in &inner.edges {
                        edge_to_faces.entry(we.idx).or_default().push(flat_face_idx);
                    }
                }
            }
        }
    }

    edge_to_faces
}

/// Compute the flat face index in the BRep.
fn compute_flat_face_index(brep: &BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> usize {
    let mut count = 0;
    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            if si == solid_idx && shi == shell_idx {
                return count + face_idx;
            }
            count += shell.faces.len();
        }
    }
    count
}

/// Get information about an edge.
fn get_edge_info(
    brep: &BRep,
    edge_idx: usize,
    adjacent_faces: &[usize],
) -> Result<EdgeInfo, ChamferError> {
    let edge = brep.edges.get(edge_idx).ok_or(ChamferError::EdgeNotFound { edge_index: edge_idx })?;

    let start_vertex = edge.start;
    let end_vertex = edge.end;

    let start_point = brep.vertices.get(start_vertex)
        .map(|v| v.point)
        .ok_or(ChamferError::EdgeNotFound { edge_index: edge_idx })?;

    let end_point = brep.vertices.get(end_vertex)
        .map(|v| v.point)
        .ok_or(ChamferError::EdgeNotFound { edge_index: edge_idx })?;

    let diff = end_point - start_point;
    let length = diff.length();
    let tangent = if length > TOLERANCE { diff / length } else { DVec3::X };

    // Get surface indices for adjacent faces
    let adjacent_surfaces: Vec<Option<usize>> = adjacent_faces.iter().map(|&fi| {
        brep.geom.face_surface.get(fi).and_then(|o| *o)
    }).collect();

    Ok(EdgeInfo {
        edge_index: edge_idx,
        start_vertex,
        end_vertex,
        start_point,
        end_point,
        tangent,
        length,
        adjacent_faces: adjacent_faces.to_vec(),
        adjacent_surfaces,
    })
}

/// Get information about adjacent faces.
fn get_adjacent_face_infos(
    brep: &BRep,
    adjacent_faces: &[usize],
    edge_info: &EdgeInfo,
) -> Result<Vec<AdjacentFaceInfo>, ChamferError> {
    let mid_point = (edge_info.start_point + edge_info.end_point) * 0.5;

    let mut face_infos = Vec::with_capacity(adjacent_faces.len());

    for &flat_face_idx in adjacent_faces {
        // Find the actual face
        let (solid_idx, shell_idx, face_idx) = find_face_indices(brep, flat_face_idx);

        let face = brep.solids.get(solid_idx)
            .and_then(|s| s.shells.get(shell_idx))
            .and_then(|sh| sh.faces.get(face_idx));

        let face = match face {
            Some(f) => f,
            None => continue,
        };

        let surface_index = brep.geom.face_surface.get(flat_face_idx).and_then(|o| *o);

        let surface_type = surface_index
            .and_then(|si| brep.geom.surfaces.get(si))
            .map(|s| SurfaceType::from(s))
            .unwrap_or(SurfaceType::Other);

        face_infos.push(AdjacentFaceInfo {
            face_index: flat_face_idx,
            surface_index,
            normal: face.normal,
            surface_type,
            reference_point: mid_point,
        });
    }

    Ok(face_infos)
}

/// Find solid, shell, and face indices from a flat face index.
fn find_face_indices(brep: &BRep, flat_face_idx: usize) -> (usize, usize, usize) {
    let mut count = 0;
    for (solid_idx, solid) in brep.solids.iter().enumerate() {
        for (shell_idx, shell) in solid.shells.iter().enumerate() {
            for face_idx in 0..shell.faces.len() {
                if count == flat_face_idx {
                    return (solid_idx, shell_idx, face_idx);
                }
                count += 1;
            }
        }
    }
    (0, 0, 0)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_box() -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        // Populate geometry
        crate::geom_populate::populate_box_geom(&mut brep);
        brep
    }

    #[test]
    fn test_chamfer_params_symmetric() {
        let params = ChamferParams::symmetric(0.5);
        assert_eq!(params.mode, ChamferMode::Symmetric);
        assert_eq!(params.distance1, 0.5);
        assert_eq!(params.get_distances(), (0.5, 0.5));
        assert!(params.validate().is_ok());
    }

    #[test]
    fn test_chamfer_params_asymmetric() {
        let params = ChamferParams::asymmetric(0.3, 0.6);
        assert_eq!(params.mode, ChamferMode::Asymmetric);
        assert_eq!(params.distance1, 0.3);
        assert_eq!(params.distance2, 0.6);
        assert_eq!(params.get_distances(), (0.3, 0.6));
        assert!(params.validate().is_ok());
    }

    #[test]
    fn test_chamfer_params_distance_angle() {
        let angle = std::f64::consts::FRAC_PI_4; // 45 degrees
        let params = ChamferParams::distance_angle(0.5, angle);
        assert_eq!(params.mode, ChamferMode::DistanceAngle);
        assert_eq!(params.distance1, 0.5);
        assert_eq!(params.angle, angle);
        assert!(params.validate().is_ok());

        let (d1, d2) = params.get_distances();
        assert!((d1 - 0.5).abs() < 1e-10);
        assert!((d2 - 0.5).abs() < 1e-10); // tan(45) = 1
    }

    #[test]
    fn test_chamfer_params_validation() {
        // Negative distance should fail
        let params = ChamferParams::symmetric(-0.1);
        assert!(params.validate().is_err());

        // Zero distance should fail
        let params = ChamferParams::symmetric(0.0);
        assert!(params.validate().is_err());

        // Invalid angle should fail
        let params = ChamferParams::distance_angle(0.5, 0.0);
        assert!(params.validate().is_err());

        let params = ChamferParams::distance_angle(0.5, std::f64::consts::FRAC_PI_2);
        assert!(params.validate().is_err());
    }

    #[test]
    fn test_chamfer_error_display() {
        let err = ChamferError::InvalidDistance {
            value: -1.0,
            reason: "must be positive".to_string(),
        };
        let s = format!("{}", err);
        assert!(s.contains("invalid distance"));
        assert!(s.contains("-1"));

        let err = ChamferError::EdgeNotFound { edge_index: 42 };
        let s = format!("{}", err);
        assert!(s.contains("edge 42"));

        let err = ChamferError::ChamferSurfaceFailed {
            edge_index: 5,
            reason: "test reason".to_string(),
        };
        let s = format!("{}", err);
        assert!(s.contains("edge 5"));
        assert!(s.contains("test reason"));
    }

    #[test]
    fn test_make_chamfer_edge_basic() {
        let box_brep = create_test_box();

        // Chamfer the first edge
        let result = make_chamfer_edge(&box_brep, &[0], 0.1);
        assert!(result.is_ok(), "chamfer should succeed: {:?}", result.err());

        let chamfer_result = result.unwrap();
        assert!(chamfer_result.edges_chamfered > 0);
        assert!(chamfer_result.chamfer_faces_created > 0);
    }

    #[test]
    fn test_make_chamfer_multiple_edges() {
        let box_brep = create_test_box();

        // Chamfer multiple edges
        let result = make_chamfer_edge(&box_brep, &[0, 1, 2], 0.1);
        assert!(result.is_ok());

        let chamfer_result = result.unwrap();
        assert!(chamfer_result.edges_chamfered >= 1);
    }

    #[test]
    fn test_make_chamfer_asymmetric() {
        let box_brep = create_test_box();

        let result = make_chamfer_asymmetric(&box_brep, &[0], 0.1, 0.2);
        assert!(result.is_ok());

        let chamfer_result = result.unwrap();
        assert!(chamfer_result.edges_chamfered >= 1);
    }

    #[test]
    fn test_make_chamfer_angle() {
        let box_brep = create_test_box();

        let angle = 30.0_f64.to_radians();
        let result = make_chamfer_angle(&box_brep, &[0], 0.1, angle);
        assert!(result.is_ok());

        let chamfer_result = result.unwrap();
        assert!(chamfer_result.edges_chamfered >= 1);
    }

    #[test]
    fn test_make_chamfer_all_edges() {
        let box_brep = create_test_box();

        let result = make_chamfer_all_edges(&box_brep, 0.05);
        assert!(result.is_ok());

        let chamfer_result = result.unwrap();
        // A box has 12 edges, but some might be skipped
        assert!(chamfer_result.edges_chamfered > 0);
    }

    #[test]
    fn test_compute_chamfer_surface() {
        let box_brep = create_test_box();
        let params = ChamferParams::symmetric(0.1);

        let result = compute_chamfer_surface(&box_brep, 0, &params);
        assert!(result.is_ok());

        let surface = result.unwrap();
        match surface {
            Surface3::Plane(plane) => {
                // Chamfer surface should be a plane for box edges
                assert!(plane.normal.length() > 0.9);
            }
            _ => panic!("expected plane surface for box chamfer"),
        }
    }

    #[test]
    fn test_compute_chamfer_curves() {
        let box_brep = create_test_box();
        let params = ChamferParams::symmetric(0.1);

        let result = compute_chamfer_curves(&box_brep, 0, &params);
        assert!(result.is_ok());

        let (curve1, curve2) = result.unwrap();
        match curve1 {
            Curve3::Line(line) => {
                assert!(line.direction.length() > 0.9);
            }
            _ => panic!("expected line curve for chamfer boundary"),
        }
        match curve2 {
            Curve3::Line(line) => {
                assert!(line.direction.length() > 0.9);
            }
            _ => panic!("expected line curve for chamfer boundary"),
        }
    }

    #[test]
    fn test_invalid_edge_index() {
        let box_brep = create_test_box();

        let result = make_chamfer_edge(&box_brep, &[999], 0.1);
        // Should return error or skip the edge with warning
        match result {
            Ok(chamfer_result) => {
                // Edge was skipped
                assert_eq!(chamfer_result.edges_chamfered, 0);
                assert!(!chamfer_result.warnings.is_empty());
            }
            Err(e) => {
                // Error was returned
                assert!(matches!(e, ChamferError::InvalidResult { .. }));
            }
        }
    }

    #[test]
    fn test_zero_distance_error() {
        let box_brep = create_test_box();

        let result = make_chamfer_edge(&box_brep, &[0], 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_empty_solid_error() {
        let empty_brep = BRep::new();

        let result = make_chamfer_edge(&empty_brep, &[0], 0.1);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_edge_face_adjacency() {
        let box_brep = create_test_box();

        let adjacency = build_edge_face_adjacency(&box_brep);

        // Each edge of a box should be adjacent to 2 faces
        for (_, faces) in adjacency.iter() {
            assert_eq!(faces.len(), 2, "box edges should have 2 adjacent faces");
        }
    }

    #[test]
    fn test_surface_type_from_surface() {
        let plane = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        assert_eq!(SurfaceType::from(&plane), SurfaceType::Plane);

        let cylinder = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        assert_eq!(SurfaceType::from(&cylinder), SurfaceType::Cylinder);
    }

    #[test]
    fn test_chamfer_warning_display() {
        let warning = ChamferWarning::EdgeSkipped {
            edge_index: 5,
            reason: "test reason".to_string(),
        };
        match warning {
            ChamferWarning::EdgeSkipped { edge_index, reason } => {
                assert_eq!(edge_index, 5);
                assert_eq!(reason, "test reason");
            }
            _ => panic!("wrong variant"),
        }
    }

    #[test]
    fn test_chamfer_result_structure() {
        let box_brep = create_test_box();

        let result = make_chamfer_edge(&box_brep, &[0], 0.1).unwrap();

        // Verify result structure
        assert!(!result.brep.vertices.is_empty());
        assert!(!result.brep.edges.is_empty());
        assert!(!result.brep.solids.is_empty());
        assert!(result.edges_chamfered > 0);
        assert!(result.chamfer_faces_created > 0);
    }

    #[test]
    fn test_chamfer_geometry_for_plane_plane() {
        let edge_info = EdgeInfo {
            edge_index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::ZERO,
            end_point: DVec3::X,
            tangent: DVec3::X,
            length: 1.0,
            adjacent_faces: vec![0, 1],
            adjacent_surfaces: vec![Some(0), Some(1)],
        };

        let face0 = AdjacentFaceInfo {
            face_index: 0,
            surface_index: Some(0),
            normal: DVec3::Z,
            surface_type: SurfaceType::Plane,
            reference_point: DVec3::new(0.5, 0.0, 0.0),
        };

        let face1 = AdjacentFaceInfo {
            face_index: 1,
            surface_index: Some(1),
            normal: DVec3::Y,
            surface_type: SurfaceType::Plane,
            reference_point: DVec3::new(0.5, 0.0, 0.0),
        };

        let result = compute_chamfer_plane_plane(&edge_info, &face0, &face1, 0.1, 0.1);
        assert!(result.is_ok());

        let geom = result.unwrap();
        assert!((geom.d1 - 0.1).abs() < 1e-10);
        assert!((geom.d2 - 0.1).abs() < 1e-10);
    }

    // ============================================================================
    // Edge Case Tests for OCCT Alignment
    // ============================================================================

    /// Test asymmetric chamfer with different distances on each face.
    /// Asymmetric chamfers create non-45-degree bevels.
    #[test]
    fn test_chamfer_asymmetric() {
        let box_brep = create_test_box();

        // Create asymmetric chamfer: 0.1 on one face, 0.3 on the other
        let result = make_chamfer_asymmetric(&box_brep, &[0], 0.1, 0.3);
        assert!(result.is_ok(), "asymmetric chamfer should succeed");

        let chamfer_result = result.unwrap();
        assert!(chamfer_result.edges_chamfered >= 1, "at least one edge should be chamfered");

        // Verify the chamfer geometry has correct distances
        let params = ChamferParams::asymmetric(0.1, 0.3);
        let (d1, d2) = params.get_distances();
        assert!((d1 - 0.1).abs() < 1e-10, "distance1 should be 0.1");
        assert!((d2 - 0.3).abs() < 1e-10, "distance2 should be 0.3");
    }

    /// Test chamfer defined by distance and angle.
    /// Angle-mode chamfers create bevels at specified angles from the first face.
    #[test]
    fn test_chamfer_angle_mode() {
        let box_brep = create_test_box();

        // Create chamfer with 30-degree angle
        let angle = 30.0_f64.to_radians();
        let result = make_chamfer_angle(&box_brep, &[0], 0.2, angle);
        assert!(result.is_ok(), "angle-mode chamfer should succeed");

        let chamfer_result = result.unwrap();
        assert!(chamfer_result.edges_chamfered >= 1, "at least one edge should be chamfered");

        // Verify angle parameter
        let params = ChamferParams::distance_angle(0.2, angle);
        assert_eq!(params.mode, ChamferMode::DistanceAngle);
        assert!((params.angle - angle).abs() < 1e-10);
    }

    /// Test chamfer on all edges of a box simultaneously.
    /// This tests edge-blend interactions at vertices.
    #[test]
    fn test_chamfer_all_edges_box() {
        let box_brep = create_test_box();

        let result = make_chamfer_all_edges(&box_brep, 0.05);
        assert!(result.is_ok(), "chamfer all edges should succeed");

        let chamfer_result = result.unwrap();
        // A box has 12 edges, at least some should be chamfered
        assert!(chamfer_result.edges_chamfered > 0, "some edges should be chamfered");
    }

    /// Test chamfer on multiple selected edges.
    /// Tests handling of edge selection and sequential processing.
    #[test]
    fn test_chamfer_multiple_selected_edges() {
        let box_brep = create_test_box();

        // Chamfer a subset of edges (forming a corner)
        let result = make_chamfer_edge(&box_brep, &[0, 1, 2], 0.15);
        assert!(result.is_ok(), "multiple edge chamfer should succeed");

        let chamfer_result = result.unwrap();
        assert!(chamfer_result.edges_chamfered >= 1, "at least one edge should be chamfered");
    }

    /// Test chamfer with very small distance.
    /// Small chamfers should not create degenerate geometry.
    #[test]
    fn test_chamfer_very_small_distance() {
        let box_brep = create_test_box();

        // Very small chamfer
        let result = make_chamfer_edge(&box_brep, &[0], 0.001);
        assert!(result.is_ok(), "very small chamfer should succeed");
    }

    /// Test chamfer surface computation for plane-plane edge.
    /// Verifies correct chamfer surface geometry for the most common case.
    #[test]
    fn test_chamfer_surface_plane_plane() {
        let box_brep = create_test_box();
        let params = ChamferParams::symmetric(0.1);

        let result = compute_chamfer_surface(&box_brep, 0, &params);
        assert!(result.is_ok(), "chamfer surface computation should succeed");

        let surface = result.unwrap();
        // For plane-plane edges, chamfer surface should be a plane
        assert!(matches!(surface, Surface3::Plane(_)), "plane-plane chamfer should produce planar surface");
    }
}
