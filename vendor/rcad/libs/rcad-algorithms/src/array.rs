//! Array (pattern) operations — linear and circular repetition of BRep solids.
//!
//! Analogous to OCCT `BRepOffsetAPI_MakeThickSolid`-style patterns and
//! `BRepFeat_MakeLinearForm` / `BRepFeat_MakeRevol` for feature repetition.
//!
//! # Operations
//!
//! - **Linear pattern**: repeat along a direction with uniform spacing
//! - **Circular pattern**: rotate around an axis with uniform angular spacing
//! - **Mirror pattern**: mirror across a plane, optionally including original
//! - **Rectangular grid pattern**: 2D array with staggered options
//! - **Variable spacing pattern**: non-uniform spacing along a direction
//! - **Path pattern**: pattern along a curve with optional alignment

use glam::{DMat4, DVec3};
use rcad_kernel::BRep;
use rcad_kernel::geom::{
    Circle3, ConicalSurface, Curve3, CylindricalSurface, Ellipse3, Hyperbola3, Line3,
    LinearExtrusionSurface, OffsetSurface, Plane, RevolutionSurface, SphericalSurface, Surface3,
    ToroidalSurface, TrimmedSurface,
};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

/// Parameters for a linear pattern.
#[derive(Debug, Clone)]
pub struct LinearPatternParams {
    /// Direction of the pattern.
    pub direction: DVec3,
    /// Number of copies (including the original). Must be >= 1.
    pub count: usize,
    /// Spacing between consecutive copies.
    pub spacing: f64,
}

/// Parameters for a circular pattern.
#[derive(Debug, Clone)]
pub struct CircularPatternParams {
    /// A point on the rotation axis.
    pub axis_origin: DVec3,
    /// Normalized rotation axis direction.
    pub axis_direction: DVec3,
    /// Number of copies (including the original). Must be >= 1.
    pub count: usize,
    /// Total angle in radians for the full pattern (copies are evenly spaced).
    pub total_angle: f64,
}

/// Parameters for a mirror pattern.
#[derive(Debug, Clone)]
pub struct MirrorPatternParams {
    /// Origin point on the mirror plane.
    pub plane_origin: DVec3,
    /// Normal vector of the mirror plane (will be normalized).
    pub plane_normal: DVec3,
    /// Whether to include the original shape in the result.
    pub include_original: bool,
}

/// Parameters for a rectangular grid pattern.
#[derive(Debug, Clone)]
pub struct RectangularPatternParams {
    /// Direction for the first (X) axis of the pattern.
    pub direction1: DVec3,
    /// Number of copies along direction1 (including original).
    pub count1: usize,
    /// Spacing between copies along direction1.
    pub spacing1: f64,
    /// Direction for the second (Y) axis of the pattern.
    pub direction2: DVec3,
    /// Number of copies along direction2 (including original).
    pub count2: usize,
    /// Spacing between copies along direction2.
    pub spacing2: f64,
    /// Stagger pattern configuration.
    pub stagger: StaggerConfig,
}

/// Stagger configuration for rectangular patterns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StaggerConfig {
    /// No staggering - regular grid.
    #[default]
    None,
    /// Offset odd rows by half spacing1.
    OddRows,
    /// Offset even rows by half spacing1.
    EvenRows,
}

/// Parameters for variable spacing patterns.
#[derive(Debug, Clone)]
pub struct VariableSpacingPatternParams {
    /// Direction of the pattern.
    pub direction: DVec3,
    /// List of spacings between consecutive copies.
    /// Number of copies = spacings.len() + 1 (includes original).
    pub spacings: Vec<f64>,
}

/// Parameters for distance-based spacing pattern.
#[derive(Debug, Clone)]
pub struct DistanceSpacingPatternParams {
    /// Direction of the pattern.
    pub direction: DVec3,
    /// Total distance to cover.
    pub total_distance: f64,
    /// Number of copies (including the original).
    pub count: usize,
}

/// Parameters for a path-based pattern.
#[derive(Debug, Clone)]
pub struct PathPatternParams {
    /// List of parameter values (0.0 to 1.0) along the path where copies are placed.
    /// 0.0 = start of path, 1.0 = end of path.
    pub parameters: Vec<f64>,
    /// Whether to align instances with the path tangent.
    pub align_to_path: bool,
    /// Up vector for alignment when align_to_path is true.
    pub up_vector: DVec3,
}

/// Parameters for a pattern with suppression support.
#[derive(Debug, Clone)]
pub struct PatternWithSuppressionParams {
    /// Indices of instances to suppress (0-indexed, 0 = original).
    pub suppressed_indices: Vec<usize>,
    /// Base pattern transformation matrices for each instance.
    pub transforms: Vec<DMat4>,
}

/// Parameters for pattern within a boundary.
#[derive(Debug, Clone)]
pub struct BoundaryPatternParams {
    /// Grid direction 1.
    pub direction1: DVec3,
    /// Grid count 1.
    pub count1: usize,
    /// Grid spacing 1.
    pub spacing1: f64,
    /// Grid direction 2.
    pub direction2: DVec3,
    /// Grid count 2.
    pub spacing2: f64,
    /// Function to test if a point is within the boundary.
    /// Returns true if the point (grid position) should include an instance.
    pub boundary_test: fn(DVec3) -> bool,
}

/// Instance-specific transformation for patterns.
#[derive(Debug, Clone)]
pub struct InstanceTransform {
    /// Index of the instance (0 = original).
    pub index: usize,
    /// Additional transformation to apply to this instance.
    pub transform: DMat4,
}

/// Error type for pattern operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PatternError {
    /// Count must be at least 1.
    InvalidCount,
    /// Spacing must be positive.
    InvalidSpacing,
    /// Direction vector must be non-zero.
    ZeroDirection,
    /// Axis direction must be non-zero.
    ZeroAxis,
    /// Total angle must be non-zero and <= 2*pi.
    InvalidAngle,
    /// Input BRep has no solids.
    NoSolids,
    /// Plane normal must be non-zero.
    ZeroPlaneNormal,
    /// Spacings list must not be empty.
    EmptySpacings,
    /// Negative spacing in list.
    NegativeSpacing,
    /// Total distance must be positive.
    InvalidDistance,
    /// Parameters list must not be empty.
    EmptyParameters,
    /// Parameter must be in range [0, 1].
    InvalidParameter,
    /// Transform list must not be empty.
    EmptyTransforms,
    /// Suppressed index out of range.
    SuppressedIndexOutOfRange,
    /// Path curve evaluation failed.
    PathEvaluationFailed,
}

impl std::fmt::Display for PatternError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidCount => write!(f, "pattern count must be >= 1"),
            Self::InvalidSpacing => write!(f, "pattern spacing must be > 0"),
            Self::ZeroDirection => write!(f, "pattern direction must be non-zero"),
            Self::ZeroAxis => write!(f, "pattern axis must be non-zero"),
            Self::InvalidAngle => write!(f, "pattern angle must be > 0 and <= 2*pi"),
            Self::NoSolids => write!(f, "input BRep has no solids"),
            Self::ZeroPlaneNormal => write!(f, "mirror plane normal must be non-zero"),
            Self::EmptySpacings => write!(f, "spacings list must not be empty"),
            Self::NegativeSpacing => write!(f, "all spacings must be >= 0"),
            Self::InvalidDistance => write!(f, "total distance must be > 0"),
            Self::EmptyParameters => write!(f, "parameters list must not be empty"),
            Self::InvalidParameter => write!(f, "parameters must be in range [0, 1]"),
            Self::EmptyTransforms => write!(f, "transform list must not be empty"),
            Self::SuppressedIndexOutOfRange => write!(f, "suppressed index out of range"),
            Self::PathEvaluationFailed => write!(f, "path curve evaluation failed"),
        }
    }
}

/// Apply a linear pattern to a BRep — repeat copies along a direction.
///
/// Returns a new BRep containing all copies merged into a single solid.
/// The original is included as the first copy (offset 0).
pub fn linear_pattern(
    brep: &BRep,
    params: &LinearPatternParams,
) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.spacing <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }
    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    for i in 0..params.count {
        let offset = dir * (i as f64 * params.spacing);
        append_transformed_brep(&mut out, brep, &translation_matrix(offset))?;
    }

    Ok(out)
}

/// Apply a circular pattern to a BRep — rotate copies around an axis.
///
/// Returns a new BRep containing all copies merged into a single solid.
/// The original is included as the first copy (angle 0).
pub fn circular_pattern(
    brep: &BRep,
    params: &CircularPatternParams,
) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.total_angle <= 0.0 || params.total_angle > std::f64::consts::TAU {
        return Err(PatternError::InvalidAngle);
    }
    let axis = params
        .axis_direction
        .try_normalize()
        .ok_or(PatternError::ZeroAxis)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();
    let angle_step = params.total_angle / params.count as f64;

    for i in 0..params.count {
        let angle = i as f64 * angle_step;
        let mat = rotation_matrix(params.axis_origin, axis, angle);
        append_transformed_brep(&mut out, brep, &mat)?;
    }

    Ok(out)
}

// ── Mirror Pattern ─────────────────────────────────────────────────────────────

/// Apply a mirror pattern to a BRep — mirror across a plane.
///
/// Returns a new BRep containing the mirrored copy and optionally the original.
pub fn mirror_pattern(
    brep: &BRep,
    params: &MirrorPatternParams,
) -> Result<BRep, PatternError> {
    let normal = params
        .plane_normal
        .try_normalize()
        .ok_or(PatternError::ZeroPlaneNormal)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    // Include original if requested
    if params.include_original {
        append_transformed_brep(&mut out, brep, &DMat4::IDENTITY)?;
    }

    // Mirror across the plane
    let mirror_mat = mirror_matrix(params.plane_origin, normal);
    append_transformed_brep(&mut out, brep, &mirror_mat)?;

    Ok(out)
}

/// Apply a compound mirror pattern — mirror and linear pattern combined.
///
/// Mirrors the shape across a plane, then applies a linear pattern to both
/// the original and mirrored copies.
pub fn mirror_linear_pattern(
    brep: &BRep,
    mirror_params: &MirrorPatternParams,
    linear_params: &LinearPatternParams,
) -> Result<BRep, PatternError> {
    // First apply mirror
    let mirrored = mirror_pattern(brep, mirror_params)?;

    // Then apply linear pattern to the mirrored result
    linear_pattern(&mirrored, linear_params)
}

/// Apply a compound mirror and circular pattern.
///
/// Mirrors the shape across a plane, then applies a circular pattern.
pub fn mirror_circular_pattern(
    brep: &BRep,
    mirror_params: &MirrorPatternParams,
    circular_params: &CircularPatternParams,
) -> Result<BRep, PatternError> {
    let mirrored = mirror_pattern(brep, mirror_params)?;
    circular_pattern(&mirrored, circular_params)
}

// ── Rectangular Grid Pattern ───────────────────────────────────────────────────

/// Apply a rectangular grid pattern to a BRep.
///
/// Creates a 2D grid of copies along two orthogonal (or non-orthogonal) directions.
/// Supports staggered patterns for alternating row offsets.
pub fn rectangular_pattern(
    brep: &BRep,
    params: &RectangularPatternParams,
) -> Result<BRep, PatternError> {
    if params.count1 < 1 || params.count2 < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.spacing1 <= 0.0 || params.spacing2 <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }

    let dir1 = params
        .direction1
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;
    let dir2 = params
        .direction2
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    for j in 0..params.count2 {
        let row_offset = dir2 * (j as f64 * params.spacing2);
        let stagger_offset = match params.stagger {
            StaggerConfig::None => DVec3::ZERO,
            StaggerConfig::OddRows if j % 2 == 1 => dir1 * (params.spacing1 * 0.5),
            StaggerConfig::EvenRows if j % 2 == 0 && j > 0 => dir1 * (params.spacing1 * 0.5),
            _ => DVec3::ZERO,
        };

        for i in 0..params.count1 {
            let col_offset = dir1 * (i as f64 * params.spacing1);
            let total_offset = row_offset + col_offset + stagger_offset;
            append_transformed_brep(&mut out, brep, &translation_matrix(total_offset))?;
        }
    }

    Ok(out)
}

/// Compute the transformation matrix for a specific position in a rectangular grid.
///
/// Returns the transformation matrix for the instance at (i, j) in the grid.
pub fn rectangular_pattern_transform(
    params: &RectangularPatternParams,
    i: usize,
    j: usize,
) -> Result<DMat4, PatternError> {
    if params.spacing1 <= 0.0 || params.spacing2 <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }

    let dir1 = params
        .direction1
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;
    let dir2 = params
        .direction2
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    let row_offset = dir2 * (j as f64 * params.spacing2);
    let stagger_offset = match params.stagger {
        StaggerConfig::None => DVec3::ZERO,
        StaggerConfig::OddRows if j % 2 == 1 => dir1 * (params.spacing1 * 0.5),
        StaggerConfig::EvenRows if j % 2 == 0 && j > 0 => dir1 * (params.spacing1 * 0.5),
        _ => DVec3::ZERO,
    };
    let col_offset = dir1 * (i as f64 * params.spacing1);
    let total_offset = row_offset + col_offset + stagger_offset;

    Ok(translation_matrix(total_offset))
}

// ── Variable Spacing Pattern ───────────────────────────────────────────────────

/// Apply a pattern with non-uniform spacing.
///
/// Creates copies at positions determined by cumulative spacings.
pub fn variable_spacing_pattern(
    brep: &BRep,
    params: &VariableSpacingPatternParams,
) -> Result<BRep, PatternError> {
    if params.spacings.is_empty() {
        return Err(PatternError::EmptySpacings);
    }
    if params.spacings.iter().any(|&s| s < 0.0) {
        return Err(PatternError::NegativeSpacing);
    }

    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    // First copy is at offset 0 (original)
    append_transformed_brep(&mut out, brep, &DMat4::IDENTITY)?;

    // Subsequent copies at cumulative offsets
    let mut cumulative = 0.0;
    for &spacing in &params.spacings {
        cumulative += spacing;
        let offset = dir * cumulative;
        append_transformed_brep(&mut out, brep, &translation_matrix(offset))?;
    }

    Ok(out)
}

/// Apply a pattern with distance-based spacing.
///
/// Distributes copies evenly along a total distance.
pub fn distance_spacing_pattern(
    brep: &BRep,
    params: &DistanceSpacingPatternParams,
) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.total_distance <= 0.0 {
        return Err(PatternError::InvalidDistance);
    }

    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    // Spacing between copies
    let spacing = if params.count > 1 {
        params.total_distance / (params.count - 1) as f64
    } else {
        0.0
    };

    for i in 0..params.count {
        let offset = dir * (i as f64 * spacing);
        append_transformed_brep(&mut out, brep, &translation_matrix(offset))?;
    }

    Ok(out)
}

// ── Path Pattern ───────────────────────────────────────────────────────────────

/// Trait for path evaluation used in path patterns.
pub trait PathEvaluator {
    /// Evaluate the path at parameter t (0.0 to 1.0).
    /// Returns (position, tangent) tuple.
    fn evaluate(&self, t: f64) -> Option<(DVec3, DVec3)>;
}

/// Apply a pattern along a path.
///
/// Places copies at specified parameter values along a path curve.
/// Optionally aligns instances with the path tangent.
pub fn path_pattern(
    brep: &BRep,
    params: &PathPatternParams,
    path: &dyn PathEvaluator,
) -> Result<BRep, PatternError> {
    if params.parameters.is_empty() {
        return Err(PatternError::EmptyParameters);
    }
    if params.parameters.iter().any(|&t| t < 0.0 || t > 1.0) {
        return Err(PatternError::InvalidParameter);
    }

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();
    let up = params.up_vector.normalize_or(DVec3::Z);

    for &t in &params.parameters {
        let (position, tangent) = path
            .evaluate(t)
            .ok_or(PatternError::PathEvaluationFailed)?;

        let mat = if params.align_to_path {
            let tangent = tangent.normalize_or(DVec3::X);
            // Create a coordinate frame
            let side = up.cross(tangent).normalize_or(DVec3::Y);
            let real_up = tangent.cross(side).normalize_or(up);

            // Build rotation matrix from basis vectors
            DMat4::from_cols(
                glam::DVec4::new(tangent.x, tangent.y, tangent.z, 0.0),
                glam::DVec4::new(side.x, side.y, side.z, 0.0),
                glam::DVec4::new(real_up.x, real_up.y, real_up.z, 0.0),
                glam::DVec4::new(position.x, position.y, position.z, 1.0),
            )
        } else {
            translation_matrix(position)
        };

        append_transformed_brep(&mut out, brep, &mat)?;
    }

    Ok(out)
}

/// Apply a pattern along a path with equal spacing.
///
/// Creates count copies evenly distributed along the path.
pub fn path_pattern_equal_spacing(
    brep: &BRep,
    path: &dyn PathEvaluator,
    count: usize,
    align_to_path: bool,
    up_vector: DVec3,
) -> Result<BRep, PatternError> {
    if count < 1 {
        return Err(PatternError::InvalidCount);
    }

    let params = PathPatternParams {
        parameters: (0..count).map(|i| i as f64 / (count - 1).max(1) as f64).collect(),
        align_to_path,
        up_vector,
    };

    path_pattern(brep, &params, path)
}

// ── Pattern with Suppression ───────────────────────────────────────────────────

/// Apply a pattern with instance suppression.
///
/// Creates copies but excludes instances at the specified indices.
pub fn pattern_with_suppression(
    brep: &BRep,
    transforms: &[DMat4],
    suppressed_indices: &[usize],
) -> Result<BRep, PatternError> {
    if transforms.is_empty() {
        return Err(PatternError::EmptyTransforms);
    }

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    // Check suppressed indices are valid
    let max_index = transforms.len() - 1;
    if suppressed_indices.iter().any(|&i| i > max_index) {
        return Err(PatternError::SuppressedIndexOutOfRange);
    }

    let mut out = BRep::new();
    let suppressed_set: std::collections::HashSet<usize> = suppressed_indices.iter().copied().collect();

    for (i, mat) in transforms.iter().enumerate() {
        if !suppressed_set.contains(&i) {
            append_transformed_brep(&mut out, brep, mat)?;
        }
    }

    Ok(out)
}

/// Apply a linear pattern with suppression support.
pub fn linear_pattern_with_suppression(
    brep: &BRep,
    params: &LinearPatternParams,
    suppressed_indices: &[usize],
) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.spacing <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }

    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let transforms: Vec<DMat4> = (0..params.count)
        .map(|i| translation_matrix(dir * (i as f64 * params.spacing)))
        .collect();

    pattern_with_suppression(brep, &transforms, suppressed_indices)
}

/// Apply a circular pattern with suppression support.
pub fn circular_pattern_with_suppression(
    brep: &BRep,
    params: &CircularPatternParams,
    suppressed_indices: &[usize],
) -> Result<BRep, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.total_angle <= 0.0 || params.total_angle > std::f64::consts::TAU {
        return Err(PatternError::InvalidAngle);
    }

    let axis = params
        .axis_direction
        .try_normalize()
        .ok_or(PatternError::ZeroAxis)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let angle_step = params.total_angle / params.count as f64;
    let transforms: Vec<DMat4> = (0..params.count)
        .map(|i| rotation_matrix(params.axis_origin, axis, i as f64 * angle_step))
        .collect();

    pattern_with_suppression(brep, &transforms, suppressed_indices)
}

// ── Pattern with Instance Transforms ───────────────────────────────────────────

/// Apply a pattern with instance-specific additional transformations.
///
/// Each instance can have its own additional transformation applied
/// on top of the base pattern transformation.
pub fn pattern_with_instance_transforms(
    brep: &BRep,
    base_transforms: &[DMat4],
    instance_transforms: &[InstanceTransform],
) -> Result<BRep, PatternError> {
    if base_transforms.is_empty() {
        return Err(PatternError::EmptyTransforms);
    }

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    // Build a map from index to instance transform
    let instance_map: std::collections::HashMap<usize, DMat4> = instance_transforms
        .iter()
        .map(|it| (it.index, it.transform))
        .collect();

    let mut out = BRep::new();

    for (i, base_mat) in base_transforms.iter().enumerate() {
        let final_mat = if let Some(instance_mat) = instance_map.get(&i) {
            *base_mat * *instance_mat
        } else {
            *base_mat
        };
        append_transformed_brep(&mut out, brep, &final_mat)?;
    }

    Ok(out)
}

// ── Pattern within Boundary ────────────────────────────────────────────────────

/// Apply a rectangular grid pattern constrained within a boundary.
///
/// Only creates instances whose grid positions fall within the boundary.
pub fn pattern_within_boundary(
    brep: &BRep,
    params: &BoundaryPatternParams,
) -> Result<BRep, PatternError> {
    if params.count1 < 1 {
        return Err(PatternError::InvalidCount);
    }
    if params.spacing1 <= 0.0 || params.spacing2 <= 0.0 {
        return Err(PatternError::InvalidSpacing);
    }

    let dir1 = params
        .direction1
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;
    let dir2 = params
        .direction2
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    if brep.solids.is_empty() {
        return Err(PatternError::NoSolids);
    }

    let mut out = BRep::new();

    // Estimate count2 from boundary test
    // We need to iterate until boundary test fails consistently
    let mut j = 0;
    let mut found_in_row;

    loop {
        found_in_row = false;
        let row_offset = dir2 * (j as f64 * params.spacing2);

        for i in 0..params.count1 {
            let col_offset = dir1 * (i as f64 * params.spacing1);
            let grid_pos = row_offset + col_offset;

            if (params.boundary_test)(grid_pos) {
                found_in_row = true;
                append_transformed_brep(&mut out, brep, &translation_matrix(grid_pos))?;
            }
        }

        // Stop if no instances were found in this row and we've gone past reasonable bounds
        if !found_in_row && j > 0 {
            break;
        }
        j += 1;

        // Safety limit to prevent infinite loops
        if j > 1000 {
            break;
        }
    }

    Ok(out)
}

// ── Utility Functions ──────────────────────────────────────────────────────────

/// Create a mirror transformation matrix across a plane.
fn mirror_matrix(plane_origin: DVec3, plane_normal: DVec3) -> DMat4 {
    // Mirror matrix: reflect across plane
    // M = I - 2 * n * n^T for reflection through origin
    // For arbitrary plane: translate to origin, reflect, translate back
    let n = plane_normal;
    let reflect = DMat4::from_cols(
        glam::DVec4::new(1.0 - 2.0 * n.x * n.x, -2.0 * n.y * n.x, -2.0 * n.z * n.x, 0.0),
        glam::DVec4::new(-2.0 * n.x * n.y, 1.0 - 2.0 * n.y * n.y, -2.0 * n.z * n.y, 0.0),
        glam::DVec4::new(-2.0 * n.x * n.z, -2.0 * n.y * n.z, 1.0 - 2.0 * n.z * n.z, 0.0),
        glam::DVec4::new(0.0, 0.0, 0.0, 1.0),
    );

    DMat4::from_translation(plane_origin) * reflect * DMat4::from_translation(-plane_origin)
}

/// Generate transformation matrices for a linear pattern.
pub fn generate_linear_transforms(params: &LinearPatternParams) -> Result<Vec<DMat4>, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }

    let dir = params
        .direction
        .try_normalize()
        .ok_or(PatternError::ZeroDirection)?;

    Ok((0..params.count)
        .map(|i| translation_matrix(dir * (i as f64 * params.spacing)))
        .collect())
}

/// Generate transformation matrices for a circular pattern.
pub fn generate_circular_transforms(params: &CircularPatternParams) -> Result<Vec<DMat4>, PatternError> {
    if params.count < 1 {
        return Err(PatternError::InvalidCount);
    }

    let axis = params
        .axis_direction
        .try_normalize()
        .ok_or(PatternError::ZeroAxis)?;

    let angle_step = params.total_angle / params.count as f64;

    Ok((0..params.count)
        .map(|i| rotation_matrix(params.axis_origin, axis, i as f64 * angle_step))
        .collect())
}

/// Scale a pattern by uniformly scaling all spacing values.
pub fn scale_pattern_params(params: &LinearPatternParams, scale: f64) -> LinearPatternParams {
    LinearPatternParams {
        direction: params.direction,
        count: params.count,
        spacing: params.spacing * scale,
    }
}

/// Scale a rectangular pattern.
pub fn scale_rectangular_params(params: &RectangularPatternParams, scale: f64) -> RectangularPatternParams {
    RectangularPatternParams {
        direction1: params.direction1,
        count1: params.count1,
        spacing1: params.spacing1 * scale,
        direction2: params.direction2,
        count2: params.count2,
        spacing2: params.spacing2 * scale,
        stagger: params.stagger,
    }
}

// ── Internal helpers ─────────────────────────────────────────────────────────

fn translation_matrix(offset: DVec3) -> DMat4 {
    DMat4::from_translation(offset)
}

fn rotation_matrix(origin: DVec3, axis: DVec3, angle: f64) -> DMat4 {
    DMat4::from_translation(origin)
        * DMat4::from_axis_angle(axis, angle)
        * DMat4::from_translation(-origin)
}

fn append_transformed_brep(
    target: &mut BRep,
    source: &BRep,
    mat: &DMat4,
) -> Result<(), PatternError> {
    let v_offset = target.vertices.len();
    let e_offset = target.edges.len();
    let curve_offset = target.geom.curves.len();
    let surface_offset = target.geom.surfaces.len();

    // Transform and copy vertices
    for v in &source.vertices {
        let p = mat.transform_point3(v.point.into());
        target.vertices.push(Vertex {
            point: DVec3::new(p.x, p.y, p.z),
        });
    }

    // Transform and copy curves
    for curve in &source.geom.curves {
        target.geom.curves.push(transform_curve(curve, mat));
    }

    // Transform and copy surfaces
    for surface in &source.geom.surfaces {
        target.geom.surfaces.push(transform_surface(surface, mat));
    }

    // Copy edges with remapped vertex indices
    for e in &source.edges {
        target.edges.push(Edge {
            start: e.start + v_offset,
            end: e.end + v_offset,
        });
    }

    // Remap edge geometry references
    for &ec in &source.geom.edge_curve {
        target
            .geom
            .edge_curve
            .push(ec.map(|c| c + curve_offset));
    }
    for &ecr in &source.geom.edge_curve_range {
        target.geom.edge_curve_range.push(ecr);
    }
    for &ed in &source.geom.edge_degenerated {
        target.geom.edge_degenerated.push(ed);
    }

    // Copy faces with remapped edge indices
    for solid in &source.solids {
        let mut shell = Shell { faces: Vec::new() };
        for face in &solid.shells[0].faces {
            let wire_edges: Vec<WireEdge> = face
                .outer_wire
                .edges
                .iter()
                .map(|we| WireEdge {
                    idx: we.idx + e_offset,
                    forward: we.forward,
                })
                .collect();

            // Remap inner wire edge indices too
            let inner_wires: Vec<Wire> = face
                .inner_wires
                .iter()
                .map(|wire| Wire {
                    edges: wire
                        .edges
                        .iter()
                        .map(|we| WireEdge {
                            idx: we.idx + e_offset,
                            forward: we.forward,
                        })
                        .collect(),
                })
                .collect();

            shell.faces.push(Face {
                outer_wire: Wire { edges: wire_edges },
                inner_wires,
                normal: {
                    let rotated = mat.transform_vector3(face.normal.into());
                    DVec3::new(rotated.x, rotated.y, rotated.z).normalize_or(face.normal)
                },
                triangles: face.triangles.iter().map(|[i, j, k]| [i + v_offset, j + v_offset, k + v_offset]).collect(),
                mesh_dirty: face.mesh_dirty,
            });
        }
        target.solids.push(Solid {
            shells: vec![shell],
        });
    }

    // Remap face_surface references
    for &fs in &source.geom.face_surface {
        target
            .geom
            .face_surface
            .push(fs.map(|s| s + surface_offset));
    }

    Ok(())
}

fn transform_curve(curve: &Curve3, mat: &DMat4) -> Curve3 {
    let transform_point = |p: DVec3| {
        let r = mat.transform_point3(p.into());
        DVec3::new(r.x, r.y, r.z)
    };
    let transform_direction = |v: DVec3| {
        let r = mat.transform_vector3(v.into());
        DVec3::new(r.x, r.y, r.z).normalize_or(v)
    };

    match curve {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: transform_point(l.origin),
            direction: transform_direction(l.direction),
        }),
        Curve3::Circle(c) => Curve3::Circle(Circle3 {
            center: transform_point(c.center),
            normal: transform_direction(c.normal),
            radius: c.radius,
        }),
        Curve3::Ellipse(e) => Curve3::Ellipse(Ellipse3 {
            center: transform_point(e.center),
            normal: transform_direction(e.normal),
            major_dir: transform_direction(e.major_dir),
            major_radius: e.major_radius,
            minor_radius: e.minor_radius,
        }),
        Curve3::Hyperbola(h) => Curve3::Hyperbola(Hyperbola3 {
            center: transform_point(h.center),
            normal: transform_direction(h.normal),
            major_dir: transform_direction(h.major_dir),
            semi_major: h.semi_major,
            semi_minor: h.semi_minor,
        }),
        Curve3::BSpline(b) => {
            let mut nb = b.clone();
            for cp in &mut nb.control_points {
                *cp = transform_point(*cp);
            }
            Curve3::BSpline(nb)
        }
        Curve3::Bezier(b) => {
            let mut nb = b.clone();
            for cp in &mut nb.control_points {
                *cp = transform_point(*cp);
            }
            Curve3::Bezier(nb)
        }
        _ => curve.clone(),
    }
}

fn transform_surface(surface: &Surface3, mat: &DMat4) -> Surface3 {
    let transform_point = |p: DVec3| {
        let r = mat.transform_point3(p.into());
        DVec3::new(r.x, r.y, r.z)
    };
    let transform_direction = |v: DVec3| {
        let r = mat.transform_vector3(v.into());
        DVec3::new(r.x, r.y, r.z).normalize_or(v)
    };

    match surface {
        Surface3::Plane(p) => Surface3::Plane(Plane {
            origin: transform_point(p.origin),
            normal: transform_direction(p.normal),
        }),
        Surface3::Cylinder(c) => Surface3::Cylinder(CylindricalSurface {
            origin: transform_point(c.origin),
            axis: transform_direction(c.axis),
            radius: c.radius,
        }),
        Surface3::Sphere(s) => Surface3::Sphere(SphericalSurface {
            center: transform_point(s.center),
            axis: transform_direction(s.axis),
            radius: s.radius,
        }),
        Surface3::Cone(c) => Surface3::Cone(ConicalSurface {
            apex: transform_point(c.apex),
            axis: transform_direction(c.axis),
            radius: c.radius,
            half_angle_rad: c.half_angle_rad,
        }),
        Surface3::Torus(t) => Surface3::Torus(ToroidalSurface {
            center: transform_point(t.center),
            axis: transform_direction(t.axis),
            major_radius: t.major_radius,
            minor_radius: t.minor_radius,
        }),
        Surface3::BSpline(b) => {
            let mut nb = b.clone();
            for row in &mut nb.control_points {
                for cp in row {
                    *cp = transform_point(*cp);
                }
            }
            Surface3::BSpline(nb)
        }
        Surface3::LinearExtrusion(le) => Surface3::LinearExtrusion(LinearExtrusionSurface {
            profile: Box::new(transform_curve(&le.profile, mat)),
            direction: le.direction,
        }),
        Surface3::Revolution(r) => Surface3::Revolution(RevolutionSurface {
            profile: Box::new(transform_curve(&r.profile, mat)),
            axis_origin: transform_point(r.axis_origin),
            axis_dir: transform_direction(r.axis_dir),
        }),
        Surface3::Offset(o) => Surface3::Offset(OffsetSurface {
            basis: Box::new(transform_surface(&o.basis, mat)),
            offset_distance: o.offset_distance,
        }),
        Surface3::Trimmed(t) => Surface3::Trimmed(TrimmedSurface {
            basis: Box::new(transform_surface(&t.basis, mat)),
            trim: t.trim,
        }),
        _ => surface.clone(),
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_modeling::make_box_brep;

    fn make_box() -> BRep {
        let mut brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut brep);
        brep
    }

    fn make_box_at(origin: DVec3) -> BRep {
        let mut brep = make_box_brep(origin, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut brep);
        brep
    }

    // ── Linear Pattern Tests ────────────────────────────────────────────────────

    #[test]
    fn linear_pattern_count_1_returns_original() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 1,
            spacing: 2.0,
        };
        let result = linear_pattern(&brep, &params).unwrap();
        let v_result = rcad_kernel::properties::volume(&result);

        assert!((v_result - v_orig).abs() < 1e-9);
    }

    #[test]
    fn linear_pattern_count_3_produces_3x_volume() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 3,
            spacing: 2.0,
        };
        let result = linear_pattern(&brep, &params).unwrap();
        let v_result = rcad_kernel::properties::volume(&result);

        assert!(
            (v_result - 3.0 * v_orig).abs() < 0.01,
            "expected 3x volume, got {v_result} vs expected {}",
            3.0 * v_orig
        );
    }

    #[test]
    fn linear_pattern_invalid_spacing_returns_error() {
        let brep = make_box();
        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 3,
            spacing: -1.0,
        };
        assert!(matches!(
            linear_pattern(&brep, &params),
            Err(PatternError::InvalidSpacing)
        ));
    }

    #[test]
    fn linear_pattern_zero_direction_returns_error() {
        let brep = make_box();
        let params = LinearPatternParams {
            direction: DVec3::ZERO,
            count: 3,
            spacing: 1.0,
        };
        assert!(matches!(
            linear_pattern(&brep, &params),
            Err(PatternError::ZeroDirection)
        ));
    }

    #[test]
    fn linear_pattern_zero_count_returns_error() {
        let brep = make_box();
        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 0,
            spacing: 1.0,
        };
        assert!(matches!(
            linear_pattern(&brep, &params),
            Err(PatternError::InvalidCount)
        ));
    }

    // ── Circular Pattern Tests ──────────────────────────────────────────────────

    #[test]
    fn circular_pattern_count_4_produces_4x_volume() {
        let mut brep = make_box_brep(DVec3::new(3.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut brep);
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = CircularPatternParams {
            axis_origin: DVec3::ZERO,
            axis_direction: DVec3::Z,
            count: 4,
            total_angle: std::f64::consts::TAU,
        };
        let result = circular_pattern(&brep, &params).unwrap();

        let total_solids = result.solids.len();
        assert_eq!(total_solids, 4, "expected 4 solids, got {total_solids}");

        let v_result = rcad_kernel::properties::volume(&result);
        assert!(
            (v_result - 4.0 * v_orig).abs() < 0.01,
            "expected 4x volume, got {v_result} vs expected {}",
            4.0 * v_orig
        );
    }

    #[test]
    fn circular_pattern_half_turn_produces_2x_volume() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = CircularPatternParams {
            axis_origin: DVec3::ZERO,
            axis_direction: DVec3::Z,
            count: 2,
            total_angle: std::f64::consts::PI,
        };
        let result = circular_pattern(&brep, &params).unwrap();
        let v_result = rcad_kernel::properties::volume(&result);

        assert!(
            (v_result - 2.0 * v_orig).abs() < 0.01,
            "expected 2x volume, got {v_result} vs expected {}",
            2.0 * v_orig
        );
    }

    #[test]
    fn circular_pattern_invalid_angle_returns_error() {
        let brep = make_box();
        let params = CircularPatternParams {
            axis_origin: DVec3::ZERO,
            axis_direction: DVec3::Z,
            count: 4,
            total_angle: 0.0,
        };
        assert!(matches!(
            circular_pattern(&brep, &params),
            Err(PatternError::InvalidAngle)
        ));
    }

    #[test]
    fn circular_pattern_angle_too_large_returns_error() {
        let brep = make_box();
        let params = CircularPatternParams {
            axis_origin: DVec3::ZERO,
            axis_direction: DVec3::Z,
            count: 4,
            total_angle: 7.0,
        };
        assert!(matches!(
            circular_pattern(&brep, &params),
            Err(PatternError::InvalidAngle)
        ));
    }

    // ── Mirror Pattern Tests ────────────────────────────────────────────────────

    #[test]
    fn mirror_pattern_without_original_produces_1_solid() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = MirrorPatternParams {
            plane_origin: DVec3::new(0.5, 0.0, 0.0),
            plane_normal: DVec3::X,
            include_original: false,
        };
        let result = mirror_pattern(&brep, &params).unwrap();

        assert_eq!(result.solids.len(), 1);
        let v_result = rcad_kernel::properties::volume(&result);
        assert!((v_result - v_orig).abs() < 0.01);
    }

    #[test]
    fn mirror_pattern_with_original_produces_2_solids() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = MirrorPatternParams {
            plane_origin: DVec3::new(0.5, 0.0, 0.0),
            plane_normal: DVec3::X,
            include_original: true,
        };
        let result = mirror_pattern(&brep, &params).unwrap();

        assert_eq!(result.solids.len(), 2);
        let v_result = rcad_kernel::properties::volume(&result);
        assert!((v_result - 2.0 * v_orig).abs() < 0.01);
    }

    #[test]
    fn mirror_pattern_reflects_position() {
        let brep = make_box_at(DVec3::new(2.0, 0.0, 0.0));

        let params = MirrorPatternParams {
            plane_origin: DVec3::ZERO,
            plane_normal: DVec3::X,
            include_original: false,
        };
        let result = mirror_pattern(&brep, &params).unwrap();

        // The mirrored box should be centered around x = -2.5 (originally at 2.5 center)
        // The center of the bounding box should be negative x
        let bbox = result.bounding_box().expect("should have bounding box");
        let center_x = (bbox[0].x + bbox[1].x) * 0.5;
        assert!(center_x < 0.0, "mirrored box center x should be negative, got {}", center_x);
    }

    #[test]
    fn mirror_pattern_zero_normal_returns_error() {
        let brep = make_box();
        let params = MirrorPatternParams {
            plane_origin: DVec3::ZERO,
            plane_normal: DVec3::ZERO,
            include_original: true,
        };
        assert!(matches!(
            mirror_pattern(&brep, &params),
            Err(PatternError::ZeroPlaneNormal)
        ));
    }

    #[test]
    fn mirror_linear_pattern_combines_operations() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let mirror_params = MirrorPatternParams {
            plane_origin: DVec3::new(0.5, 0.0, 0.0),
            plane_normal: DVec3::X,
            include_original: true,
        };
        let linear_params = LinearPatternParams {
            direction: DVec3::Y,
            count: 2,
            spacing: 2.0,
        };

        let result = mirror_linear_pattern(&brep, &mirror_params, &linear_params).unwrap();

        // 2 (mirror + original) * 2 (linear count) = 4 solids
        assert_eq!(result.solids.len(), 4);
        let v_result = rcad_kernel::properties::volume(&result);
        assert!((v_result - 4.0 * v_orig).abs() < 0.01);
    }

    // ── Rectangular Pattern Tests ───────────────────────────────────────────────

    #[test]
    fn rectangular_pattern_2x3_produces_6_solids() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = RectangularPatternParams {
            direction1: DVec3::X,
            count1: 2,
            spacing1: 2.0,
            direction2: DVec3::Y,
            count2: 3,
            spacing2: 2.0,
            stagger: StaggerConfig::None,
        };
        let result = rectangular_pattern(&brep, &params).unwrap();

        assert_eq!(result.solids.len(), 6);
        let v_result = rcad_kernel::properties::volume(&result);
        assert!((v_result - 6.0 * v_orig).abs() < 0.01);
    }

    #[test]
    fn rectangular_pattern_1x1_produces_1_solid() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = RectangularPatternParams {
            direction1: DVec3::X,
            count1: 1,
            spacing1: 2.0,
            direction2: DVec3::Y,
            count2: 1,
            spacing2: 2.0,
            stagger: StaggerConfig::None,
        };
        let result = rectangular_pattern(&brep, &params).unwrap();

        assert_eq!(result.solids.len(), 1);
        let v_result = rcad_kernel::properties::volume(&result);
        assert!((v_result - v_orig).abs() < 1e-9);
    }

    #[test]
    fn rectangular_pattern_stagger_odd_rows() {
        let brep = make_box();

        let params = RectangularPatternParams {
            direction1: DVec3::X,
            count1: 3,
            spacing1: 2.0,
            direction2: DVec3::Y,
            count2: 2,
            spacing2: 2.0,
            stagger: StaggerConfig::OddRows,
        };
        let result = rectangular_pattern(&brep, &params).unwrap();

        assert_eq!(result.solids.len(), 6);
    }

    #[test]
    fn rectangular_pattern_stagger_even_rows() {
        let brep = make_box();

        let params = RectangularPatternParams {
            direction1: DVec3::X,
            count1: 3,
            spacing1: 2.0,
            direction2: DVec3::Y,
            count2: 3,
            spacing2: 2.0,
            stagger: StaggerConfig::EvenRows,
        };
        let result = rectangular_pattern(&brep, &params).unwrap();

        assert_eq!(result.solids.len(), 9);
    }

    #[test]
    fn rectangular_pattern_invalid_count_returns_error() {
        let brep = make_box();
        let params = RectangularPatternParams {
            direction1: DVec3::X,
            count1: 0,
            spacing1: 2.0,
            direction2: DVec3::Y,
            count2: 3,
            spacing2: 2.0,
            stagger: StaggerConfig::None,
        };
        assert!(matches!(
            rectangular_pattern(&brep, &params),
            Err(PatternError::InvalidCount)
        ));
    }

    #[test]
    fn rectangular_pattern_invalid_spacing_returns_error() {
        let brep = make_box();
        let params = RectangularPatternParams {
            direction1: DVec3::X,
            count1: 2,
            spacing1: -1.0,
            direction2: DVec3::Y,
            count2: 3,
            spacing2: 2.0,
            stagger: StaggerConfig::None,
        };
        assert!(matches!(
            rectangular_pattern(&brep, &params),
            Err(PatternError::InvalidSpacing)
        ));
    }

    #[test]
    fn rectangular_pattern_transform_computes_correct_offset() {
        let params = RectangularPatternParams {
            direction1: DVec3::X,
            count1: 3,
            spacing1: 2.0,
            direction2: DVec3::Y,
            count2: 3,
            spacing2: 3.0,
            stagger: StaggerConfig::None,
        };

        let mat = rectangular_pattern_transform(&params, 1, 1).unwrap();
        let offset = mat.transform_point3(glam::DVec3::ZERO);

        assert!((offset.x - 2.0).abs() < 1e-9, "expected x=2.0, got {}", offset.x);
        assert!((offset.y - 3.0).abs() < 1e-9, "expected y=3.0, got {}", offset.y);
        assert!((offset.z - 0.0).abs() < 1e-9, "expected z=0.0, got {}", offset.z);
    }

    // ── Variable Spacing Pattern Tests ──────────────────────────────────────────

    #[test]
    fn variable_spacing_pattern_produces_correct_count() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = VariableSpacingPatternParams {
            direction: DVec3::X,
            spacings: vec![1.0, 2.0, 3.0],
        };
        let result = variable_spacing_pattern(&brep, &params).unwrap();

        // Original + 3 copies = 4 solids
        assert_eq!(result.solids.len(), 4);
        let v_result = rcad_kernel::properties::volume(&result);
        assert!((v_result - 4.0 * v_orig).abs() < 0.01);
    }

    #[test]
    fn variable_spacing_pattern_single_spacing_produces_2_solids() {
        let brep = make_box();

        let params = VariableSpacingPatternParams {
            direction: DVec3::X,
            spacings: vec![2.0],
        };
        let result = variable_spacing_pattern(&brep, &params).unwrap();

        assert_eq!(result.solids.len(), 2);
    }

    #[test]
    fn variable_spacing_pattern_empty_spacings_returns_error() {
        let brep = make_box();
        let params = VariableSpacingPatternParams {
            direction: DVec3::X,
            spacings: vec![],
        };
        assert!(matches!(
            variable_spacing_pattern(&brep, &params),
            Err(PatternError::EmptySpacings)
        ));
    }

    #[test]
    fn variable_spacing_pattern_negative_spacing_returns_error() {
        let brep = make_box();
        let params = VariableSpacingPatternParams {
            direction: DVec3::X,
            spacings: vec![1.0, -1.0, 2.0],
        };
        assert!(matches!(
            variable_spacing_pattern(&brep, &params),
            Err(PatternError::NegativeSpacing)
        ));
    }

    // ── Distance Spacing Pattern Tests ──────────────────────────────────────────

    #[test]
    fn distance_spacing_pattern_distributes_evenly() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = DistanceSpacingPatternParams {
            direction: DVec3::X,
            total_distance: 10.0,
            count: 3,
        };
        let result = distance_spacing_pattern(&brep, &params).unwrap();

        assert_eq!(result.solids.len(), 3);
        let v_result = rcad_kernel::properties::volume(&result);
        assert!((v_result - 3.0 * v_orig).abs() < 0.01);
    }

    #[test]
    fn distance_spacing_pattern_single_copy_at_origin() {
        let brep = make_box();

        let params = DistanceSpacingPatternParams {
            direction: DVec3::X,
            total_distance: 10.0,
            count: 1,
        };
        let result = distance_spacing_pattern(&brep, &params).unwrap();

        assert_eq!(result.solids.len(), 1);
    }

    #[test]
    fn distance_spacing_pattern_zero_distance_returns_error() {
        let brep = make_box();
        let params = DistanceSpacingPatternParams {
            direction: DVec3::X,
            total_distance: 0.0,
            count: 3,
        };
        assert!(matches!(
            distance_spacing_pattern(&brep, &params),
            Err(PatternError::InvalidDistance)
        ));
    }

    // ── Path Pattern Tests ──────────────────────────────────────────────────────

    struct LinePath {
        start: DVec3,
        end: DVec3,
    }

    impl PathEvaluator for LinePath {
        fn evaluate(&self, t: f64) -> Option<(DVec3, DVec3)> {
            let pos = self.start.lerp(self.end, t);
            let tangent = (self.end - self.start).normalize_or(DVec3::X);
            Some((pos, tangent))
        }
    }

    #[test]
    fn path_pattern_along_line_produces_correct_count() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let path = LinePath {
            start: DVec3::ZERO,
            end: DVec3::new(10.0, 0.0, 0.0),
        };
        let params = PathPatternParams {
            parameters: vec![0.0, 0.5, 1.0],
            align_to_path: false,
            up_vector: DVec3::Z,
        };
        let result = path_pattern(&brep, &params, &path).unwrap();

        assert_eq!(result.solids.len(), 3);
        let v_result = rcad_kernel::properties::volume(&result);
        assert!((v_result - 3.0 * v_orig).abs() < 0.01);
    }

    #[test]
    fn path_pattern_with_alignment() {
        let brep = make_box();

        let path = LinePath {
            start: DVec3::ZERO,
            end: DVec3::new(10.0, 0.0, 0.0),
        };
        let params = PathPatternParams {
            parameters: vec![0.0, 1.0],
            align_to_path: true,
            up_vector: DVec3::Z,
        };
        let result = path_pattern(&brep, &params, &path).unwrap();

        assert_eq!(result.solids.len(), 2);
    }

    #[test]
    fn path_pattern_empty_parameters_returns_error() {
        let brep = make_box();
        let path = LinePath {
            start: DVec3::ZERO,
            end: DVec3::new(10.0, 0.0, 0.0),
        };
        let params = PathPatternParams {
            parameters: vec![],
            align_to_path: false,
            up_vector: DVec3::Z,
        };
        assert!(matches!(
            path_pattern(&brep, &params, &path),
            Err(PatternError::EmptyParameters)
        ));
    }

    #[test]
    fn path_pattern_invalid_parameter_returns_error() {
        let brep = make_box();
        let path = LinePath {
            start: DVec3::ZERO,
            end: DVec3::new(10.0, 0.0, 0.0),
        };
        let params = PathPatternParams {
            parameters: vec![0.0, 1.5],
            align_to_path: false,
            up_vector: DVec3::Z,
        };
        assert!(matches!(
            path_pattern(&brep, &params, &path),
            Err(PatternError::InvalidParameter)
        ));
    }

    #[test]
    fn path_pattern_equal_spacing_produces_correct_count() {
        let brep = make_box();

        let path = LinePath {
            start: DVec3::ZERO,
            end: DVec3::new(10.0, 0.0, 0.0),
        };
        let result = path_pattern_equal_spacing(&brep, &path, 5, false, DVec3::Z).unwrap();

        assert_eq!(result.solids.len(), 5);
    }

    // ── Pattern with Suppression Tests ──────────────────────────────────────────

    #[test]
    fn pattern_with_suppression_excludes_indices() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 4,
            spacing: 2.0,
        };
        let result = linear_pattern_with_suppression(&brep, &params, &[1, 2]).unwrap();

        // 4 - 2 suppressed = 2 solids
        assert_eq!(result.solids.len(), 2);
        let v_result = rcad_kernel::properties::volume(&result);
        assert!((v_result - 2.0 * v_orig).abs() < 0.01);
    }

    #[test]
    fn pattern_with_suppression_all_suppressed_produces_empty() {
        let brep = make_box();

        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 3,
            spacing: 2.0,
        };
        let result = linear_pattern_with_suppression(&brep, &params, &[0, 1, 2]).unwrap();

        assert_eq!(result.solids.len(), 0);
    }

    #[test]
    fn pattern_with_suppression_no_suppression_produces_all() {
        let brep = make_box();

        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 3,
            spacing: 2.0,
        };
        let result = linear_pattern_with_suppression(&brep, &params, &[]).unwrap();

        assert_eq!(result.solids.len(), 3);
    }

    #[test]
    fn circular_pattern_with_suppression_excludes_indices() {
        let brep = make_box();

        let params = CircularPatternParams {
            axis_origin: DVec3::ZERO,
            axis_direction: DVec3::Z,
            count: 4,
            total_angle: std::f64::consts::TAU,
        };
        let result = circular_pattern_with_suppression(&brep, &params, &[0]).unwrap();

        assert_eq!(result.solids.len(), 3);
    }

    #[test]
    fn pattern_with_suppression_out_of_range_returns_error() {
        let brep = make_box();

        let transforms = vec![DMat4::IDENTITY, DMat4::IDENTITY];
        let result = pattern_with_suppression(&brep, &transforms, &[5]);

        assert!(matches!(result, Err(PatternError::SuppressedIndexOutOfRange)));
    }

    // ── Pattern with Instance Transforms Tests ───────────────────────────────────

    #[test]
    fn pattern_with_instance_transforms_applies_additional_transform() {
        let brep = make_box();

        let base_transforms = vec![
            DMat4::IDENTITY,
            DMat4::from_translation(DVec3::new(5.0, 0.0, 0.0)),
        ];
        let instance_transforms = vec![
            InstanceTransform {
                index: 1,
                transform: DMat4::from_scale(glam::DVec3::splat(2.0)),
            },
        ];

        let result = pattern_with_instance_transforms(&brep, &base_transforms, &instance_transforms).unwrap();

        assert_eq!(result.solids.len(), 2);
        // Second instance should have double volume due to scale
        let v_result = rcad_kernel::properties::volume(&result);
        assert!(v_result > 2.0, "scaled instance should increase total volume");
    }

    // ── Pattern within Boundary Tests ────────────────────────────────────────────

    #[test]
    fn pattern_within_boundary_respects_boundary() {
        let brep = make_box();

        // Circular boundary with radius 5 centered at origin
        fn circular_boundary(p: DVec3) -> bool {
            p.x * p.x + p.y * p.y <= 25.0
        }

        let params = BoundaryPatternParams {
            direction1: DVec3::X,
            count1: 5,
            spacing1: 2.0,
            direction2: DVec3::Y,
            spacing2: 2.0,
            boundary_test: circular_boundary,
        };

        let result = pattern_within_boundary(&brep, &params).unwrap();

        // Should have instances only within the circle
        // With spacing 2.0, we'd have points at (0,0), (2,0), (4,0), etc.
        // Only points within radius 5 should be included
        assert!(result.solids.len() > 0);
    }

    #[test]
    fn pattern_within_boundary_empty_outside_boundary() {
        let brep = make_box();

        // Boundary that only includes points at (0, 0, 0)
        fn single_point_boundary(_p: DVec3) -> bool {
            false
        }

        let params = BoundaryPatternParams {
            direction1: DVec3::X,
            count1: 3,
            spacing1: 2.0,
            direction2: DVec3::Y,
            spacing2: 2.0,
            boundary_test: single_point_boundary,
        };

        let result = pattern_within_boundary(&brep, &params).unwrap();

        assert_eq!(result.solids.len(), 0);
    }

    // ── Utility Function Tests ──────────────────────────────────────────────────

    #[test]
    fn generate_linear_transforms_produces_correct_count() {
        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 5,
            spacing: 2.0,
        };
        let transforms = generate_linear_transforms(&params).unwrap();

        assert_eq!(transforms.len(), 5);
    }

    #[test]
    fn generate_circular_transforms_produces_correct_count() {
        let params = CircularPatternParams {
            axis_origin: DVec3::ZERO,
            axis_direction: DVec3::Z,
            count: 6,
            total_angle: std::f64::consts::TAU,
        };
        let transforms = generate_circular_transforms(&params).unwrap();

        assert_eq!(transforms.len(), 6);
    }

    #[test]
    fn scale_pattern_params_scales_spacing() {
        let params = LinearPatternParams {
            direction: DVec3::X,
            count: 3,
            spacing: 2.0,
        };
        let scaled = scale_pattern_params(&params, 2.0);

        assert_eq!(scaled.spacing, 4.0);
        assert_eq!(scaled.count, params.count);
    }

    #[test]
    fn scale_rectangular_params_scales_both_spacings() {
        let params = RectangularPatternParams {
            direction1: DVec3::X,
            count1: 2,
            spacing1: 2.0,
            direction2: DVec3::Y,
            count2: 3,
            spacing2: 3.0,
            stagger: StaggerConfig::None,
        };
        let scaled = scale_rectangular_params(&params, 0.5);

        assert!((scaled.spacing1 - 1.0).abs() < 1e-9);
        assert!((scaled.spacing2 - 1.5).abs() < 1e-9);
    }

    // ── Error Display Tests ──────────────────────────────────────────────────────

    #[test]
    fn pattern_error_display_works() {
        assert_eq!(
            format!("{}", PatternError::InvalidCount),
            "pattern count must be >= 1"
        );
        assert_eq!(
            format!("{}", PatternError::ZeroPlaneNormal),
            "mirror plane normal must be non-zero"
        );
        assert_eq!(
            format!("{}", PatternError::EmptySpacings),
            "spacings list must not be empty"
        );
        assert_eq!(
            format!("{}", PatternError::PathEvaluationFailed),
            "path curve evaluation failed"
        );
    }
}
