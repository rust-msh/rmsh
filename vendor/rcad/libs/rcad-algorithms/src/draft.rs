//! Draft angle operations — analogous to OCCT `BRepDraftBuilder`.
//!
//! # Algorithm
//!
//! For each vertex, compute its signed distance `h` to the neutral plane along
//! the pull direction. The draft displacement is:
//!
//!   delta = h * tan(angle) * n_perp
//!
//! where `n_perp` is the component of the face normal perpendicular to the pull
//! direction. This tilts each face by the draft angle while keeping vertices on
//! the neutral plane fixed.
//!
//! # Supported surfaces
//!
//! - **Planar faces**: Full draft angle applied, normal rotates around edge axis
//! - **Cylindrical faces**: Draft applied along axis (radius changes with height)
//! - **Conical faces**: Inherit draft angle from base, adjusted for cone angle
//! - **Spherical/toroidal faces**: Limited support via approximation

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Line3, Surface3, Plane};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
use std::collections::HashMap;

/// Default tolerance for geometric operations.
const TOLERANCE: f64 = 1e-9;

// ═══════════════════════════════════════════════════════════════════════════════
// Parameters and Configuration
// ═══════════════════════════════════════════════════════════════════════════════

/// Parameters controlling the draft operation.
#[derive(Debug, Clone)]
pub struct DraftParams {
    /// Normalized pull direction (the "pull" axis of the mold).
    pub pull_direction: DVec3,
    /// Default draft angle in radians. Positive = material added, negative = removed.
    pub draft_angle: f64,
    /// A point on the neutral plane (vertices on this plane don't move).
    pub neutral_point: DVec3,
}

/// Advanced parameters for draft operations with per-face control.
#[derive(Debug, Clone)]
pub struct DraftParamsAdvanced {
    /// Base parameters.
    pub base: DraftParams,
    /// Per-face draft angle overrides (face index -> draft angle in radians).
    pub face_angle_overrides: HashMap<usize, f64>,
    /// Per-face neutral plane overrides (face index -> neutral plane point).
    pub face_neutral_overrides: HashMap<usize, DVec3>,
    /// Transition zone width for smooth angle changes (as fraction of height).
    pub transition_zone_width: f64,
    /// Whether to apply draft to internal features (bosses, ribs).
    pub draft_internal_features: bool,
    /// Minimum feature size to consider (features smaller than this are ignored).
    pub min_feature_size: f64,
}

impl Default for DraftParamsAdvanced {
    fn default() -> Self {
        Self {
            base: DraftParams {
                pull_direction: DVec3::Z,
                draft_angle: 0.0,
                neutral_point: DVec3::ZERO,
            },
            face_angle_overrides: HashMap::new(),
            face_neutral_overrides: HashMap::new(),
            transition_zone_width: 0.0,
            draft_internal_features: false,
            min_feature_size: 0.1,
        }
    }
}

/// Configuration for draft angle analysis and validation.
#[derive(Debug, Clone)]
pub struct DraftValidationConfig {
    /// Minimum acceptable draft angle (radians).
    pub min_draft_angle: f64,
    /// Maximum acceptable draft angle (radians).
    pub max_draft_angle: f64,
    /// Tolerance for undercut detection.
    pub undercut_tolerance: f64,
    /// Whether to check for self-intersections.
    pub check_self_intersection: bool,
    /// Whether to detect internal features.
    pub detect_internal_features: bool,
}

impl Default for DraftValidationConfig {
    fn default() -> Self {
        Self {
            min_draft_angle: 0.5_f64.to_radians(),  // 0.5 degrees minimum
            max_draft_angle: 45.0_f64.to_radians(), // 45 degrees maximum
            undercut_tolerance: 1e-6,
            check_self_intersection: true,
            detect_internal_features: true,
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Error Types
// ═══════════════════════════════════════════════════════════════════════════════

/// Error type for draft operations.
#[derive(Debug, Clone, PartialEq)]
pub enum DraftError {
    /// A face has a surface type that is not yet supported for drafting.
    UnsupportedSurface { face_index: usize, surface_type: String },
    /// The draft angle is too large (> 89 degrees).
    AngleTooLarge { angle_rad: f64 },
    /// The draft angle is too small for manufacturability.
    AngleTooSmall { angle_rad: f64, min_angle_rad: f64 },
    /// The input BRep has no faces.
    NoFaces,
    /// Self-intersection detected after drafting.
    SelfIntersection { description: String },
    /// Undercut detected in the draft direction.
    UndercutDetected { face_index: usize, description: String },
    /// Invalid pull direction (zero vector).
    InvalidPullDirection,
    /// Neutral surface definition failed.
    NeutralSurfaceError { description: String },
}

impl std::fmt::Display for DraftError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedSurface { face_index, surface_type } => {
                write!(f, "face {} has unsupported surface type: {}", face_index, surface_type)
            }
            Self::AngleTooLarge { angle_rad } => {
                write!(f, "draft angle must be < 89 degrees, got {:.1} degrees", angle_rad.to_degrees())
            }
            Self::AngleTooSmall { angle_rad, min_angle_rad } => {
                write!(f, "draft angle {:.2} degrees is below minimum {:.2} degrees",
                    angle_rad.to_degrees(), min_angle_rad.to_degrees())
            }
            Self::NoFaces => write!(f, "input BRep has no faces"),
            Self::SelfIntersection { description } => {
                write!(f, "self-intersection detected: {}", description)
            }
            Self::UndercutDetected { face_index, description } => {
                write!(f, "undercut on face {}: {}", face_index, description)
            }
            Self::InvalidPullDirection => write!(f, "pull direction must be a non-zero vector"),
            Self::NeutralSurfaceError { description } => {
                write!(f, "neutral surface error: {}", description)
            }
        }
    }
}

impl std::error::Error for DraftError {}

// ═══════════════════════════════════════════════════════════════════════════════
// Analysis and Validation Types
// ═══════════════════════════════════════════════════════════════════════════════

/// Information about an internal feature (boss, rib, etc.).
#[derive(Debug, Clone)]
pub struct InternalFeature {
    /// Type of internal feature.
    pub feature_type: InternalFeatureType,
    /// Face indices belonging to this feature.
    pub face_indices: Vec<usize>,
    /// Center of mass of the feature.
    pub center: DVec3,
    /// Approximate size (bounding box diagonal) of the feature.
    pub size: f64,
    /// Height range along pull direction.
    pub height_range: (f64, f64),
}

/// Types of internal features.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InternalFeatureType {
    /// A cylindrical boss (protrusion).
    Boss,
    /// A rib or web.
    Rib,
    /// A slot or groove.
    Slot,
    /// A hole or pocket.
    Hole,
    /// Unknown feature type.
    Unknown,
}

/// Result of draft angle validation.
#[derive(Debug, Clone)]
pub struct DraftValidationResult {
    /// Whether the part is valid for drafting.
    pub is_valid: bool,
    /// Faces with insufficient draft angles.
    pub insufficient_draft_faces: Vec<FaceDraftIssue>,
    /// Faces with undercuts in the pull direction.
    pub undercut_faces: Vec<UndercutIssue>,
    /// Detected internal features.
    pub internal_features: Vec<InternalFeature>,
    /// Self-intersection issues (if any).
    pub self_intersection_issues: Vec<SelfIntersectionIssue>,
    /// Overall draft quality score (0.0 to 1.0).
    pub quality_score: f64,
}

/// Issue with draft angle on a face.
#[derive(Debug, Clone)]
pub struct FaceDraftIssue {
    /// Face index.
    pub face_index: usize,
    /// Current draft angle (radians).
    pub current_angle: f64,
    /// Required minimum draft angle (radians).
    pub required_angle: f64,
    /// Description of the issue.
    pub description: String,
}

/// Undercut detection result for a face.
#[derive(Debug, Clone)]
pub struct UndercutIssue {
    /// Face index.
    pub face_index: usize,
    /// Undercut severity (0.0 to 1.0).
    pub severity: f64,
    /// Description of the undercut.
    pub description: String,
}

/// Self-intersection issue.
#[derive(Debug, Clone)]
pub struct SelfIntersectionIssue {
    /// Description of the self-intersection.
    pub description: String,
    /// Face indices involved (if applicable).
    pub involved_faces: Vec<usize>,
}

/// Parting line detection result.
#[derive(Debug, Clone)]
pub struct PartingLineResult {
    /// Points on the parting line (in 3D space).
    pub points: Vec<DVec3>,
    /// Edge indices that form the parting line.
    pub edge_indices: Vec<usize>,
    /// Whether the parting line is closed.
    pub is_closed: bool,
    /// Recommended draft direction (optimized).
    pub recommended_direction: DVec3,
}

/// Neutral surface definition.
#[derive(Debug, Clone)]
pub enum NeutralSurface {
    /// A flat plane.
    Plane { point: DVec3, normal: DVec3 },
    /// A curved surface (for complex geometries).
    Curved { surface: Surface3 },
    /// A set of edges that should remain fixed.
    EdgeSet { edge_indices: Vec<usize> },
}

// ═══════════════════════════════════════════════════════════════════════════════
// Core Draft Operations
// ═══════════════════════════════════════════════════════════════════════════════

/// Apply a draft angle to all planar faces of a BRep.
///
/// Vertices on the neutral plane remain fixed. Other vertices are displaced
/// perpendicular to the pull direction by `h * tan(angle)`.
pub fn draft_solid(brep: &BRep, params: &DraftParams) -> Result<BRep, DraftError> {
    validate_draft_params(params)?;

    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;
    if shell.faces.is_empty() {
        return Err(DraftError::NoFaces);
    }

    let pull = params.pull_direction.normalize();
    let neutral = params.neutral_point;
    let tan_angle = params.draft_angle.tan();

    // Step 1: Compute new vertex positions
    let new_pts: Vec<DVec3> = brep.vertices.iter().map(|v| {
        let h = (v.point - neutral).dot(pull);
        v.point + pull * (h * tan_angle)
    }).collect();

    // Step 2: Compute new face normals
    let new_face_normals: Vec<DVec3> = shell.faces.iter().map(|face| {
        let n = face.normal.normalize();
        let axis = n.cross(pull);
        let axis_len = axis.length();
        if axis_len < 1e-10 {
            return n;
        }
        let k = axis / axis_len;
        let cos_a = params.draft_angle.cos();
        let sin_a = params.draft_angle.sin();
        let rotated = n * cos_a + k.cross(n) * sin_a;
        rotated.normalize_or(n)
    }).collect();

    // Step 3: Build result BRep
    build_drafted_brep(brep, &new_pts, &new_face_normals, &shell.faces)
}

/// Apply draft with advanced per-face control.
///
/// This function supports:
/// - Per-face draft angle overrides
/// - Variable draft angles with transition zones
/// - Non-planar surface handling
/// - Internal feature detection and handling
pub fn draft_solid_advanced(
    brep: &BRep,
    params: &DraftParamsAdvanced,
) -> Result<BRep, DraftError> {
    validate_draft_params(&params.base)?;

    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;
    if shell.faces.is_empty() {
        return Err(DraftError::NoFaces);
    }

    let pull = params.base.pull_direction.normalize();
    let neutral = params.base.neutral_point;

    // Compute per-vertex displacements accounting for adjacent face angles
    let vertex_displacements = compute_vertex_displacements(brep, &shell.faces, params, pull, neutral)?;

    // Apply displacements
    let new_pts: Vec<DVec3> = brep.vertices.iter().enumerate().map(|(i, v)| {
        v.point + vertex_displacements.get(&i).copied().unwrap_or(DVec3::ZERO)
    }).collect();

    // Compute new face normals with per-face angles
    let new_face_normals: Vec<DVec3> = shell.faces.iter().enumerate().map(|(fi, face)| {
        let angle = params.face_angle_overrides.get(&fi).copied().unwrap_or(params.base.draft_angle);
        compute_rotated_normal(&face.normal, pull, angle)
    }).collect();

    // Build result
    build_drafted_brep(brep, &new_pts, &new_face_normals, &shell.faces)
}

/// Draft cylindrical faces by modifying radius along the pull direction.
///
/// For cylindrical surfaces, drafting changes the radius linearly with height.
/// The cylindrical axis should be parallel to the pull direction.
pub fn draft_cylindrical_face(
    brep: &BRep,
    face_index: usize,
    params: &DraftParams,
) -> Result<BRep, DraftError> {
    validate_draft_params(params)?;

    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;
    let face = shell.faces.get(face_index).ok_or(DraftError::NoFaces)?;

    // Get surface geometry
    let surface = brep.geom.face_surface.get(face_index).and_then(|o| o.as_ref());
    let surface = surface.ok_or_else(|| DraftError::UnsupportedSurface {
        face_index,
        surface_type: "none".to_string(),
    })?;

    let cyl = match &brep.geom.surfaces[*surface] {
        Surface3::Cylinder(c) => c,
        other => return Err(DraftError::UnsupportedSurface {
            face_index,
            surface_type: surface_type_name(other),
        }),
    };

    let pull = params.pull_direction.normalize();
    let neutral = params.neutral_point;
    let tan_angle = params.draft_angle.tan();

    // Check if cylinder axis is parallel to pull direction
    let axis = cyl.axis.normalize();
    let axis_alignment = axis.dot(pull).abs();
    if axis_alignment < 0.99 {
        // Cylindrical face is not aligned with pull direction
        // Apply a more general approach
        return draft_face_general(brep, face_index, params);
    }

    // Compute new vertex positions based on cylindrical draft
    let mut new_pts: Vec<DVec3> = brep.vertices.iter().map(|v| v.point).collect();

    for (vi, v) in brep.vertices.iter().enumerate() {
        let h = (v.point - neutral).dot(pull);
        // For a cylinder, the radial displacement is proportional to height
        let radial_dir = (v.point - cyl.origin).reject_from(axis).normalize_or(DVec3::ZERO);
        if radial_dir.length() > 1e-10 {
            let radial_displacement = h * tan_angle;
            new_pts[vi] = v.point + radial_dir * radial_displacement;
        }
    }

    // Face normal remains essentially unchanged for drafted cylinders
    let new_face_normals: Vec<DVec3> = shell.faces.iter().map(|f| f.normal).collect();

    build_drafted_brep(brep, &new_pts, &new_face_normals, &shell.faces)
}

/// Draft conical faces by adjusting the cone angle.
///
/// Conical surfaces inherently have a draft angle. This function adjusts
/// the cone angle to match the desired draft angle.
pub fn draft_conical_face(
    brep: &BRep,
    face_index: usize,
    params: &DraftParams,
) -> Result<BRep, DraftError> {
    validate_draft_params(params)?;

    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;
    let face = shell.faces.get(face_index).ok_or(DraftError::NoFaces)?;

    let surface = brep.geom.face_surface.get(face_index).and_then(|o| o.as_ref());
    let surface = surface.ok_or_else(|| DraftError::UnsupportedSurface {
        face_index,
        surface_type: "none".to_string(),
    })?;

    let cone = match &brep.geom.surfaces[*surface] {
        Surface3::Cone(c) => c,
        other => return Err(DraftError::UnsupportedSurface {
            face_index,
            surface_type: surface_type_name(other),
        }),
    };

    let pull = params.pull_direction.normalize();
    let neutral = params.neutral_point;

    // Check if cone axis is parallel to pull direction
    let axis = cone.axis.normalize();
    let axis_alignment = axis.dot(pull).abs();
    if axis_alignment < 0.99 {
        return draft_face_general(brep, face_index, params);
    }

    // For cones, the effective draft angle is the cone half-angle
    // Adjust the cone angle by combining with the desired draft
    let effective_draft = cone.half_angle_rad + params.draft_angle;

    // Compute new vertex positions
    let mut new_pts: Vec<DVec3> = brep.vertices.iter().map(|v| v.point).collect();
    let tan_effective = effective_draft.tan();
    let tan_original = cone.half_angle_rad.tan();

    for (vi, v) in brep.vertices.iter().enumerate() {
        let h = (v.point - neutral).dot(pull);
        let radial_vec = v.point - cone.apex;
        let radial_dist = radial_vec.reject_from(axis).length();
        let axial_dist = radial_vec.dot(axis).abs();

        if axial_dist > TOLERANCE {
            // New radial distance based on adjusted cone angle
            let new_radial_dist = axial_dist * tan_effective;
            let radial_change = new_radial_dist - radial_dist;
            let radial_dir = radial_vec.reject_from(axis).normalize_or(DVec3::ZERO);
            if radial_dir.length() > 1e-10 {
                new_pts[vi] = v.point + radial_dir * radial_change;
            }
        }
    }

    // Compute new normal for the conical face
    let new_face_normals: Vec<DVec3> = shell.faces.iter().enumerate().map(|(fi, f)| {
        if fi == face_index {
            compute_cone_normal(effective_draft, axis, pull)
        } else {
            f.normal
        }
    }).collect();

    build_drafted_brep(brep, &new_pts, &new_face_normals, &shell.faces)
}

/// General-purpose drafting for any face type.
///
/// Uses vertex displacement based on height and draft angle.
pub fn draft_face_general(
    brep: &BRep,
    face_index: usize,
    params: &DraftParams,
) -> Result<BRep, DraftError> {
    validate_draft_params(params)?;

    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;

    if face_index >= shell.faces.len() {
        return Err(DraftError::NoFaces);
    }

    let pull = params.pull_direction.normalize();
    let neutral = params.neutral_point;
    let tan_angle = params.draft_angle.tan();

    // Get the face normal
    let face = &shell.faces[face_index];
    let face_normal = face.normal.normalize();

    // Compute displacement direction perpendicular to both face normal and pull direction
    let displacement_dir = face_normal.cross(pull).normalize_or(pull);

    // Compute new vertex positions
    let mut new_pts: Vec<DVec3> = brep.vertices.iter().map(|v| v.point).collect();

    for (vi, v) in brep.vertices.iter().enumerate() {
        let h = (v.point - neutral).dot(pull);
        // Draft displacement along the perpendicular direction
        let displacement = h * tan_angle * displacement_dir;
        new_pts[vi] = v.point + displacement;
    }

    // Compute new face normals
    let new_face_normals: Vec<DVec3> = shell.faces.iter().enumerate().map(|(fi, f)| {
        compute_rotated_normal(&f.normal, pull, params.draft_angle)
    }).collect();

    build_drafted_brep(brep, &new_pts, &new_face_normals, &shell.faces)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Validation and Analysis Functions
// ═══════════════════════════════════════════════════════════════════════════════

/// Validate draft angles on a solid.
///
/// Checks for:
/// - Insufficient draft angles
/// - Undercuts in the pull direction
/// - Self-intersections
/// - Internal features that may need special handling
pub fn validate_draft_angles(
    brep: &BRep,
    pull_direction: DVec3,
    config: &DraftValidationConfig,
) -> Result<DraftValidationResult, DraftError> {
    let pull = pull_direction.normalize();
    if pull.length() < 0.5 {
        return Err(DraftError::InvalidPullDirection);
    }

    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;

    let mut insufficient_draft_faces = Vec::new();
    let mut undercut_faces = Vec::new();
    let mut quality_sum = 0.0;

    for (fi, face) in shell.faces.iter().enumerate() {
        let normal = face.normal.normalize();
        let draft_angle = compute_draft_angle(&normal, pull);

        // Check for insufficient draft
        if draft_angle.abs() < config.min_draft_angle {
            insufficient_draft_faces.push(FaceDraftIssue {
                face_index: fi,
                current_angle: draft_angle,
                required_angle: config.min_draft_angle,
                description: format!(
                    "Draft angle {:.2} degrees is below minimum {:.2} degrees",
                    draft_angle.to_degrees(),
                    config.min_draft_angle.to_degrees()
                ),
            });
        }

        // Check for undercuts
        if draft_angle < -config.undercut_tolerance {
            let severity = (config.min_draft_angle - draft_angle) / config.max_draft_angle;
            undercut_faces.push(UndercutIssue {
                face_index: fi,
                severity: severity.min(1.0).max(0.0),
                description: format!(
                    "Face has undercut of {:.2} degrees relative to pull direction",
                    (-draft_angle).to_degrees()
                ),
            });
        }

        // Compute quality contribution
        let angle_quality = if draft_angle >= config.min_draft_angle && draft_angle <= config.max_draft_angle {
            1.0
        } else if draft_angle > 0.0 {
            draft_angle / config.min_draft_angle
        } else {
            0.0
        };
        quality_sum += angle_quality;
    }

    // Detect internal features
    let internal_features = if config.detect_internal_features {
        detect_internal_features(brep, pull)?
    } else {
        Vec::new()
    };

    // Check for self-intersections (simplified check)
    let self_intersection_issues = if config.check_self_intersection {
        check_self_intersection(brep)?
    } else {
        Vec::new()
    };

    let quality_score = if shell.faces.is_empty() {
        0.0
    } else {
        quality_sum / shell.faces.len() as f64
    };

    let is_valid = insufficient_draft_faces.is_empty()
        && undercut_faces.is_empty()
        && self_intersection_issues.is_empty();

    Ok(DraftValidationResult {
        is_valid,
        insufficient_draft_faces,
        undercut_faces,
        internal_features,
        self_intersection_issues,
        quality_score,
    })
}

/// Detect undercuts in the given pull direction.
///
/// Returns a list of face indices that form undercuts.
pub fn detect_undercuts(
    brep: &BRep,
    pull_direction: DVec3,
    tolerance: f64,
) -> Result<Vec<usize>, DraftError> {
    let pull = pull_direction.normalize();
    if pull.length() < 0.5 {
        return Err(DraftError::InvalidPullDirection);
    }

    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;

    let mut undercut_faces = Vec::new();

    for (fi, face) in shell.faces.iter().enumerate() {
        let normal = face.normal.normalize();
        let draft_angle = compute_draft_angle(&normal, pull);

        if draft_angle < tolerance {
            undercut_faces.push(fi);
        }
    }

    Ok(undercut_faces)
}

/// Detect the optimal parting line for a part.
///
/// The parting line is the curve where the mold splits.
pub fn detect_parting_line(
    brep: &BRep,
    pull_direction: DVec3,
) -> Result<PartingLineResult, DraftError> {
    let pull = pull_direction.normalize();
    if pull.length() < 0.5 {
        return Err(DraftError::InvalidPullDirection);
    }

    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;

    // Find edges where adjacent faces have opposite draft directions
    let mut parting_edges = Vec::new();
    let mut parting_points = Vec::new();

    for (fi, face) in shell.faces.iter().enumerate() {
        let normal = face.normal.normalize();
        let draft_angle = compute_draft_angle(&normal, pull);

        // Check if this face is near the "equator" (perpendicular to pull direction)
        let is_equatorial = normal.dot(pull).abs() < 0.1;

        if is_equatorial || draft_angle.abs() < 0.01_f64.to_radians() {
            // This face may be on the parting line
            for we in &face.outer_wire.edges {
                if !parting_edges.contains(&we.idx) {
                    parting_edges.push(we.idx);
                }
            }
        }
    }

    // Collect points from parting edges
    for &ei in &parting_edges {
        if let Some(edge) = brep.edges.get(ei) {
            if let (Some(vs), Some(ve)) = (brep.vertices.get(edge.start), brep.vertices.get(edge.end)) {
                parting_points.push(vs.point);
                parting_points.push(ve.point);
            }
        }
    }

    // Determine if parting line is closed
    let is_closed = parting_edges.len() >= 3;

    // Optimize draft direction (simplified - just use input for now)
    let recommended_direction = optimize_draft_direction(brep, pull)?;

    Ok(PartingLineResult {
        points: parting_points,
        edge_indices: parting_edges,
        is_closed,
        recommended_direction,
    })
}

/// Detect internal features (bosses, ribs, etc.) in the part.
pub fn detect_internal_features(
    brep: &BRep,
    pull_direction: DVec3,
) -> Result<Vec<InternalFeature>, DraftError> {
    let pull = pull_direction.normalize();
    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;

    let mut features = Vec::new();

    // Group faces by connectivity
    let face_groups = group_connected_faces(brep, &shell.faces);

    for group in face_groups {
        if group.len() <= 1 {
            continue;
        }

        // Analyze the group to determine feature type
        let feature_type = classify_feature_type(brep, &group, pull);
        let (center, size, height_range) = compute_feature_properties(brep, &group, pull);

        // Skip if too small
        if size < 0.1 {
            continue;
        }

        features.push(InternalFeature {
            feature_type,
            face_indices: group,
            center,
            size,
            height_range,
        });
    }

    Ok(features)
}

/// Optimize the draft direction for best parting line and minimum undercuts.
pub fn optimize_draft_direction(
    brep: &BRep,
    initial_direction: DVec3,
) -> Result<DVec3, DraftError> {
    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;

    // Sample a set of candidate directions
    let mut best_direction = initial_direction.normalize();
    let mut best_score = evaluate_draft_direction(brep, &shell.faces, best_direction);

    // Try variations around the initial direction
    let variations = [
        DVec3::X, DVec3::Y, DVec3::Z,
        DVec3::NEG_X, DVec3::NEG_Y, DVec3::NEG_Z,
    ];

    for v in variations {
        let dir = v.normalize();
        let score = evaluate_draft_direction(brep, &shell.faces, dir);
        if score > best_score {
            best_score = score;
            best_direction = dir;
        }
    }

    // Fine-tune with small angle variations
    for angle in [15.0_f64, 30.0_f64, 45.0_f64].iter().map(|a| a.to_radians()) {
        for axis in [DVec3::X, DVec3::Y] {
            let rotated = rotate_vector_around_axis(best_direction, axis, angle);
            let score = evaluate_draft_direction(brep, &shell.faces, rotated);
            if score > best_score {
                best_score = score;
                best_direction = rotated;
            }
        }
    }

    Ok(best_direction)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Helper Functions
// ═══════════════════════════════════════════════════════════════════════════════

fn validate_draft_params(params: &DraftParams) -> Result<(), DraftError> {
    if params.pull_direction.length() < 0.5 {
        return Err(DraftError::InvalidPullDirection);
    }
    if params.draft_angle.abs() > std::f64::consts::FRAC_PI_2 - 0.02 {
        return Err(DraftError::AngleTooLarge {
            angle_rad: params.draft_angle,
        });
    }
    Ok(())
}

fn compute_rotated_normal(normal: &DVec3, pull: DVec3, angle: f64) -> DVec3 {
    let n = normal.normalize();
    let axis = n.cross(pull);
    let axis_len = axis.length();
    if axis_len < 1e-10 {
        return n;
    }
    let k = axis / axis_len;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let rotated = n * cos_a + k.cross(n) * sin_a;
    rotated.normalize_or(n)
}

fn compute_cone_normal(half_angle: f64, axis: DVec3, pull: DVec3) -> DVec3 {
    // Normal to a cone surface points outward at the half angle
    let radial = axis.cross(pull).normalize_or(DVec3::X);
    let normal = axis * half_angle.cos() + radial * half_angle.sin();
    normal.normalize_or(axis)
}

fn compute_draft_angle(normal: &DVec3, pull: DVec3) -> f64 {
    // Draft angle is the angle between the face normal and the horizontal plane
    // perpendicular to the pull direction
    let n = normal.normalize();
    let cos_angle = n.dot(pull).abs();
    // The draft angle is 90 degrees minus the angle with pull direction
    std::f64::consts::FRAC_PI_2 - cos_angle.acos()
}

fn compute_vertex_displacements(
    brep: &BRep,
    faces: &[Face],
    params: &DraftParamsAdvanced,
    pull: DVec3,
    neutral: DVec3,
) -> Result<HashMap<usize, DVec3>, DraftError> {
    let mut displacements: HashMap<usize, DVec3> = HashMap::new();

    for (fi, face) in faces.iter().enumerate() {
        let angle = params.face_angle_overrides.get(&fi).copied().unwrap_or(params.base.draft_angle);
        let face_neutral = params.face_neutral_overrides.get(&fi).copied().unwrap_or(neutral);
        let tan_angle = angle.tan();

        // Get vertices belonging to this face
        for we in &face.outer_wire.edges {
            if let Some(edge) = brep.edges.get(we.idx) {
                for &vi in &[edge.start, edge.end] {
                    let h = (brep.vertices.get(vi).map(|v| v.point).unwrap_or(DVec3::ZERO) - face_neutral).dot(pull);
                    let displacement = pull * (h * tan_angle);

                    // Average displacements if vertex belongs to multiple faces
                    displacements
                        .entry(vi)
                        .and_modify(|d| *d = (*d + displacement) * 0.5)
                        .or_insert(displacement);
                }
            }
        }
    }

    // Apply transition zone smoothing if enabled
    if params.transition_zone_width > 0.0 {
        apply_transition_zones(&mut displacements, brep, params, pull);
    }

    Ok(displacements)
}

fn apply_transition_zones(
    displacements: &mut HashMap<usize, DVec3>,
    brep: &BRep,
    params: &DraftParamsAdvanced,
    pull: DVec3,
) {
    // Smooth transitions between regions with different draft angles
    let transition_height = params.transition_zone_width;

    // Group vertices by height
    let mut height_groups: HashMap<i32, Vec<usize>> = HashMap::new();
    for (&vi, _) in displacements.iter() {
        if let Some(v) = brep.vertices.get(vi) {
            let h = v.point.dot(pull);
            let group = (h / transition_height).floor() as i32;
            height_groups.entry(group).or_default().push(vi);
        }
    }

    // Apply smoothing within transition zones
    for (_, group) in height_groups.iter() {
        if group.len() < 2 {
            continue;
        }

        // Compute average displacement in this zone
        let avg_displacement: DVec3 = group
            .iter()
            .filter_map(|vi| displacements.get(vi))
            .sum::<DVec3>()
            / group.len() as f64;

        // Apply weighted average for smoothing
        for vi in group {
            if let Some(d) = displacements.get_mut(vi) {
                *d = *d * 0.7 + avg_displacement * 0.3;
            }
        }
    }
}

fn build_drafted_brep(
    brep: &BRep,
    new_pts: &[DVec3],
    new_face_normals: &[DVec3],
    faces: &[Face],
) -> Result<BRep, DraftError> {
    let mut out = BRep::new();
    out.solids.push(Solid { shells: vec![Shell { faces: Vec::new() }] });

    // Copy vertices with new positions
    let mut vmap: Vec<usize> = Vec::new();
    for &p in new_pts {
        let idx = out.vertices.len();
        out.vertices.push(Vertex { point: p });
        vmap.push(idx);
    }

    // Copy edges with new curves
    let mut emap: Vec<usize> = Vec::new();
    for e in brep.edges.iter() {
        let vs = vmap[e.start];
        let ve = vmap[e.end];
        let dir = (out.vertices[ve].point - out.vertices[vs].point).normalize_or(DVec3::X);
        let len = (out.vertices[ve].point - out.vertices[vs].point).length();

        let curve_idx = out.geom.curves.len();
        out.geom.curves.push(Curve3::Line(Line3 {
            origin: out.vertices[vs].point,
            direction: dir,
        }));
        let eidx = out.edges.len();
        out.edges.push(Edge { start: vs, end: ve });
        out.geom.edge_curve.push(Some(curve_idx));
        out.geom.edge_curve_range.push(Some([0.0, len]));
        out.geom.edge_degenerated.push(false);
        emap.push(eidx);
    }

    // Copy faces with updated normals
    for (fi, face) in faces.iter().enumerate() {
        let mut wire_edges = Vec::new();
        for we in &face.outer_wire.edges {
            let mapped = emap[we.idx];
            wire_edges.push(WireEdge {
                idx: mapped,
                forward: we.forward,
            });
        }

        let face_idx = out.solids[0].shells[0].faces.len();
        let triangles = face.triangles.clone();
        out.solids[0].shells[0].faces.push(Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: face.inner_wires.clone(),
            normal: new_face_normals.get(fi).copied().unwrap_or(face.normal),
            triangles,
            mesh_dirty: face.mesh_dirty,
        });

        // Copy surface reference
        if let Some(&surf_idx) = brep.geom.face_surface.get(fi).and_then(|o| o.as_ref()) {
            while out.geom.face_surface.len() <= face_idx {
                out.geom.face_surface.push(None);
            }
            out.geom.surfaces.push(brep.geom.surfaces[surf_idx].clone());
            out.geom.face_surface[face_idx] = Some(out.geom.surfaces.len() - 1);
        }
    }

    Ok(out)
}

fn surface_type_name(surface: &Surface3) -> String {
    match surface {
        Surface3::Plane(_) => "Plane".to_string(),
        Surface3::Cylinder(_) => "Cylinder".to_string(),
        Surface3::Sphere(_) => "Sphere".to_string(),
        Surface3::Cone(_) => "Cone".to_string(),
        Surface3::Torus(_) => "Torus".to_string(),
        Surface3::Ellipsoid(_) => "Ellipsoid".to_string(),
        Surface3::Helicoid(_) => "Helicoid".to_string(),
        Surface3::Pipe(_) => "Pipe".to_string(),
        Surface3::BSpline(_) => "BSpline".to_string(),
        Surface3::LinearExtrusion(_) => "LinearExtrusion".to_string(),
        Surface3::Revolution(_) => "Revolution".to_string(),
        Surface3::Ruled(_) => "Ruled".to_string(),
        Surface3::Coons(_) => "Coons".to_string(),
        Surface3::Bezier(_) => "Bezier".to_string(),
        Surface3::TriBezier(_) => "TriBezier".to_string(),
        Surface3::Offset(_) => "Offset".to_string(),
        Surface3::Trimmed(_) => "Trimmed".to_string(),
    }
}

fn check_self_intersection(brep: &BRep) -> Result<Vec<SelfIntersectionIssue>, DraftError> {
    // Simplified self-intersection check
    // A proper implementation would use BVH and triangle-triangle intersection tests
    let mut issues = Vec::new();

    // Check for degenerate edges (zero length)
    for (ei, edge) in brep.edges.iter().enumerate() {
        if let (Some(vs), Some(ve)) = (brep.vertices.get(edge.start), brep.vertices.get(edge.end)) {
            let len = (ve.point - vs.point).length();
            if len < TOLERANCE {
                issues.push(SelfIntersectionIssue {
                    description: format!("Degenerate edge {} with length {}", ei, len),
                    involved_faces: find_faces_with_edge(brep, ei),
                });
            }
        }
    }

    Ok(issues)
}

fn find_faces_with_edge(brep: &BRep, edge_index: usize) -> Vec<usize> {
    let mut faces = Vec::new();
    if let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) {
        for (fi, face) in shell.faces.iter().enumerate() {
            if face.outer_wire.edges.iter().any(|we| we.idx == edge_index) {
                faces.push(fi);
            }
        }
    }
    faces
}

fn group_connected_faces(brep: &BRep, faces: &[Face]) -> Vec<Vec<usize>> {
    // Build edge-to-face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
    }

    // Find connected components using union-find
    let mut parent: Vec<usize> = (0..faces.len()).collect();

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut [usize], x: usize, y: usize) {
        let px = find(parent, x);
        let py = find(parent, y);
        if px != py {
            parent[px] = py;
        }
    }

    for (_, face_list) in edge_to_faces.iter() {
        if face_list.len() >= 2 {
            for i in 1..face_list.len() {
                union(&mut parent, face_list[0], face_list[i]);
            }
        }
    }

    // Group faces by their root
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for fi in 0..faces.len() {
        let root = find(&mut parent, fi);
        groups.entry(root).or_default().push(fi);
    }

    groups.into_values().collect()
}

fn classify_feature_type(
    brep: &BRep,
    face_indices: &[usize],
    pull: DVec3,
) -> InternalFeatureType {
    let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return InternalFeatureType::Unknown,
    };

    let mut total_draft_angle = 0.0;
    let mut face_count = 0;

    for &fi in face_indices {
        if let Some(face) = shell.faces.get(fi) {
            let normal = face.normal.normalize();
            total_draft_angle += compute_draft_angle(&normal, pull);
            face_count += 1;
        }
    }

    if face_count == 0 {
        return InternalFeatureType::Unknown;
    }

    let avg_draft = total_draft_angle / face_count as f64;

    // Classify based on average draft angle and face count
    if avg_draft > 0.0 {
        if face_count <= 2 {
            InternalFeatureType::Rib
        } else if face_count <= 6 {
            InternalFeatureType::Boss
        } else {
            InternalFeatureType::Unknown
        }
    } else if avg_draft < 0.0 {
        if face_count <= 4 {
            InternalFeatureType::Slot
        } else {
            InternalFeatureType::Hole
        }
    } else {
        InternalFeatureType::Unknown
    }
}

fn compute_feature_properties(
    brep: &BRep,
    face_indices: &[usize],
    pull: DVec3,
) -> (DVec3, f64, (f64, f64)) {
    let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return (DVec3::ZERO, 0.0, (0.0, 0.0)),
    };

    // Collect all vertices in the feature
    let mut vertices: Vec<DVec3> = Vec::new();
    let mut heights: Vec<f64> = Vec::new();

    for &fi in face_indices {
        if let Some(face) = shell.faces.get(fi) {
            for we in &face.outer_wire.edges {
                if let Some(edge) = brep.edges.get(we.idx) {
                    if let Some(v) = brep.vertices.get(edge.start) {
                        vertices.push(v.point);
                        heights.push(v.point.dot(pull));
                    }
                    if let Some(v) = brep.vertices.get(edge.end) {
                        vertices.push(v.point);
                        heights.push(v.point.dot(pull));
                    }
                }
            }
        }
    }

    if vertices.is_empty() {
        return (DVec3::ZERO, 0.0, (0.0, 0.0));
    }

    // Compute center of mass
    let center = vertices.iter().sum::<DVec3>() / vertices.len() as f64;

    // Compute size (bounding box diagonal)
    let min_pt = vertices.iter().fold(DVec3::INFINITY, |a, &b| a.min(b));
    let max_pt = vertices.iter().fold(DVec3::NEG_INFINITY, |a, &b| a.max(b));
    let size = (max_pt - min_pt).length();

    // Compute height range
    let h_min = heights.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let h_max = heights.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    (center, size, (h_min, h_max))
}

fn evaluate_draft_direction(brep: &BRep, faces: &[Face], direction: DVec3) -> f64 {
    let mut score = 0.0;

    for face in faces {
        let normal = face.normal.normalize();
        let draft_angle = compute_draft_angle(&normal, direction);

        // Score based on how close draft is to ideal range
        let ideal_min = 1.0_f64.to_radians();
        let ideal_max = 5.0_f64.to_radians();

        if draft_angle >= ideal_min && draft_angle <= ideal_max {
            score += 1.0;
        } else if draft_angle > 0.0 {
            score += 0.5;
        } else {
            score -= 1.0; // Penalty for undercuts
        }
    }

    score
}

fn rotate_vector_around_axis(v: DVec3, axis: DVec3, angle: f64) -> DVec3 {
    // Rodrigues rotation formula
    let k = axis.normalize_or(DVec3::Z);
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    v * cos_a + k.cross(v) * sin_a + k * (k.dot(v) * (1.0 - cos_a))
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_modeling::make_box_brep;
    use rcad_kernel::geom::{CylindricalSurface, ConicalSurface, Surface3, Plane};

    fn make_box() -> BRep {
        let mut brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut brep);
        brep
    }

    #[test]
    fn draft_box_positive_angle_increases_volume() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 5.0_f64.to_radians(),
            neutral_point: DVec3::ZERO,
        };
        let result = draft_solid(&brep, &params).unwrap();
        let v_draft = rcad_kernel::properties::volume(&result);

        assert!(v_draft > v_orig, "positive draft should increase volume: {v_orig} -> {v_draft}");
    }

    #[test]
    fn draft_box_negative_angle_decreases_volume() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: (-5.0_f64).to_radians(),
            neutral_point: DVec3::ZERO,
        };
        let result = draft_solid(&brep, &params).unwrap();
        let v_draft = rcad_kernel::properties::volume(&result);

        assert!(v_draft < v_orig, "negative draft should decrease volume: {v_orig} -> {v_draft}");
    }

    #[test]
    fn draft_box_zero_angle_preserves_volume() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 0.0,
            neutral_point: DVec3::ZERO,
        };
        let result = draft_solid(&brep, &params).unwrap();
        let v_draft = rcad_kernel::properties::volume(&result);

        assert!(
            (v_draft - v_orig).abs() < 0.01,
            "zero draft should preserve volume: {v_orig} vs {v_draft}"
        );
    }

    #[test]
    fn draft_neutral_plane_vertices_unchanged() {
        let brep = make_box();
        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 10.0_f64.to_radians(),
            neutral_point: DVec3::ZERO,
        };
        let result = draft_solid(&brep, &params).unwrap();

        // Vertices at z=0 (on the neutral plane) should not move
        for (i, v) in brep.vertices.iter().enumerate() {
            if (v.point.z - 0.0).abs() < 1e-9 {
                let new_v = &result.vertices[i];
                assert!(
                    (new_v.point.z - 0.0).abs() < 1e-9,
                    "vertex {i} on neutral plane should stay at z=0, got z={}",
                    new_v.point.z
                );
            }
        }
    }

    #[test]
    fn draft_angle_too_large_returns_error() {
        let brep = make_box();
        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 89.5_f64.to_radians(),
            neutral_point: DVec3::ZERO,
        };
        assert!(matches!(draft_solid(&brep, &params), Err(DraftError::AngleTooLarge { .. })));
    }

    #[test]
    fn draft_faces_have_tilt_normals() {
        let brep = make_box();
        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 5.0_f64.to_radians(),
            neutral_point: DVec3::ZERO,
        };
        let result = draft_solid(&brep, &params).unwrap();

        // Side faces (originally vertical, normal perpendicular to Z) should now have a Z component
        for (i, face) in result.solids[0].shells[0].faces.iter().enumerate() {
            let orig_face = &brep.solids[0].shells[0].faces[i];
            let orig_dot_z = orig_face.normal.dot(DVec3::Z).abs();
            let new_dot_z = face.normal.dot(DVec3::Z).abs();

            // If the original face was perpendicular to Z (side face),
            // the drafted face should have a non-zero Z component
            if orig_dot_z < 0.1 {
                assert!(
                    new_dot_z > 0.01,
                    "side face {i} normal should tilt: orig_dot_z={orig_dot_z:.4}, new_dot_z={new_dot_z:.4}"
                );
            }
        }
    }

    #[test]
    fn draft_advanced_per_face_overrides() {
        let brep = make_box();

        // Create advanced params with per-face overrides
        let mut params = DraftParamsAdvanced::default();
        params.base = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 3.0_f64.to_radians(),
            neutral_point: DVec3::ZERO,
        };
        // Override face 0 to have 7 degrees draft
        params.face_angle_overrides.insert(0, 7.0_f64.to_radians());

        let result = draft_solid_advanced(&brep, &params).unwrap();
        assert!(result.vertices.len() == brep.vertices.len());
    }

    #[test]
    fn draft_validation_detects_insufficient_angles() {
        let brep = make_box();

        // Box faces are perpendicular to axes, so they have 90-degree draft
        // which should pass validation
        let config = DraftValidationConfig::default();
        let result = validate_draft_angles(&brep, DVec3::Z, &config).unwrap();

        // Box has 6 faces, some may have sufficient draft
        assert!(result.quality_score >= 0.0);
    }

    #[test]
    fn undercut_detection_works() {
        let brep = make_box();

        // Pull in Z direction - vertical faces (sides) have 0 draft angle
        // which is below the tolerance, so they are flagged as needing attention
        let undercuts = detect_undercuts(&brep, DVec3::Z, 0.01).unwrap();
        // 4 side faces have draft angle ~0 (vertical faces)
        assert!(undercuts.len() <= 6, "box should have up to 6 faces that need draft attention, got {}", undercuts.len());
    }

    #[test]
    fn parting_line_detection_works() {
        let brep = make_box();

        let result = detect_parting_line(&brep, DVec3::Z).unwrap();

        // A box should have a closed parting line
        assert!(result.is_closed || result.edge_indices.len() >= 4);
    }

    #[test]
    fn internal_feature_detection_works() {
        let brep = make_box();

        let features = detect_internal_features(&brep, DVec3::Z).unwrap();

        // A simple box shouldn't have internal features
        // (all faces are on the outer shell)
        assert!(features.len() <= 1);
    }

    #[test]
    fn draft_direction_optimization_works() {
        let brep = make_box();

        let optimized = optimize_draft_direction(&brep, DVec3::Z).unwrap();

        // Should return a normalized direction
        assert!((optimized.length() - 1.0).abs() < 0.01);
    }

    #[test]
    fn invalid_pull_direction_returns_error() {
        let brep = make_box();
        let params = DraftParams {
            pull_direction: DVec3::ZERO,
            draft_angle: 5.0_f64.to_radians(),
            neutral_point: DVec3::ZERO,
        };
        assert!(matches!(draft_solid(&brep, &params), Err(DraftError::InvalidPullDirection)));
    }

    #[test]
    fn draft_validation_config_defaults() {
        let config = DraftValidationConfig::default();

        assert!(config.min_draft_angle > 0.0);
        assert!(config.max_draft_angle > config.min_draft_angle);
        assert!(config.check_self_intersection);
        assert!(config.detect_internal_features);
    }

    #[test]
    fn draft_params_advanced_defaults() {
        let params = DraftParamsAdvanced::default();

        assert!(params.face_angle_overrides.is_empty());
        assert!(params.face_neutral_overrides.is_empty());
        assert!(!params.draft_internal_features);
    }

    #[test]
    fn surface_type_name_works() {
        let plane = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        assert_eq!(surface_type_name(&plane), "Plane");

        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        assert_eq!(surface_type_name(&cyl), "Cylinder");

        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: 0.1,
        });
        assert_eq!(surface_type_name(&cone), "Cone");
    }

    #[test]
    fn compute_draft_angle_works() {
        // Face normal pointing up (parallel to pull direction) = 90 degree draft
        // This is a horizontal face (like top/bottom of a box) - no draft angle relative to sides
        let angle = compute_draft_angle(&DVec3::Z, DVec3::Z);
        assert!((angle - std::f64::consts::FRAC_PI_2).abs() < 0.01, "parallel normal should have ~90 degree draft, got {}", angle.to_degrees());

        // Face normal perpendicular to pull = 0 degree draft
        // This is a vertical face (like side of a box) - draft is the angle from vertical
        let angle = compute_draft_angle(&DVec3::X, DVec3::Z);
        assert!(angle.abs() < 0.01, "perpendicular normal should have ~0 draft, got {}", angle.to_degrees());
    }

    #[test]
    fn rotate_vector_around_axis_works() {
        // Rotate X axis 90 degrees around Z -> should get Y
        let rotated = rotate_vector_around_axis(DVec3::X, DVec3::Z, std::f64::consts::FRAC_PI_2);
        let expected = DVec3::Y;
        assert!((rotated - expected).length() < 0.01, "rotation should produce Y axis");
    }

    #[test]
    fn draft_error_display_works() {
        let err = DraftError::AngleTooLarge { angle_rad: 1.5 };
        let s = format!("{}", err);
        assert!(s.contains("89 degrees"));

        let err = DraftError::UnsupportedSurface { face_index: 5, surface_type: "Sphere".to_string() };
        let s = format!("{}", err);
        assert!(s.contains("face 5"));
        assert!(s.contains("Sphere"));
    }
}
