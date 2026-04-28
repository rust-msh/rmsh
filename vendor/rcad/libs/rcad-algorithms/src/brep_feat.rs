//! BRepFeat-style feature-based modeling operations.
//!
//! This module provides feature-based modeling operations analogous to OCCT's TKFeat
//! (BRepFeat package). Features are operations that add or remove material from a
//! base shape while maintaining design intent.
//!
//! # Feature Types
//!
//! - **Rib**: Add reinforcing rib features from a wire profile
//! - **Groove**: Create slot/groove features by removing material
//! - **Prism**: Feature-based prism (extrusion) with fuse modes
//! - **Revol**: Feature-based revolution with fuse modes
//! - **Pipe**: Feature-based pipe along a spine curve
//! - **Draft**: Apply draft angle to faces for moldability

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

use crate::{BooleanError, BooleanOpType, boolean_op};

// ═══════════════════════════════════════════════════════════════════════════════
// Error Types
// ═══════════════════════════════════════════════════════════════════════════════

/// Errors returned by BRepFeat operations.
#[derive(Debug)]
pub enum BRepFeatError {
    /// Input contains non-finite values.
    NonFiniteInput(&'static str),
    /// Input value must be positive.
    NonPositiveInput(&'static str),
    /// Invalid input geometry or parameters.
    InvalidInput(String),
    /// Zero-length vector where non-zero is required.
    ZeroVector(&'static str),
    /// Vectors are parallel when they should not be.
    ParallelVectors(&'static str, &'static str),
    /// Boolean operation failed.
    BooleanFailed(BooleanError),
    /// Modeling operation failed.
    ModelingFailed(String),
    /// Profile wire is invalid (not closed, too few edges, etc.).
    InvalidProfile(String),
    /// Feature does not intersect with target shape.
    NoIntersection,
    /// Draft angle is out of valid range.
    InvalidDraftAngle { angle_rad: f64 },
    /// Face index out of range.
    FaceNotFound { face_index: usize },
    /// Neutral plane is invalid.
    InvalidNeutralPlane,
}

impl std::fmt::Display for BRepFeatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFiniteInput(name) => write!(f, "{name} must be finite"),
            Self::NonPositiveInput(name) => write!(f, "{name} must be > 0"),
            Self::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
            Self::ZeroVector(name) => write!(f, "{name} must be non-zero"),
            Self::ParallelVectors(a, b) => write!(f, "{a} must not be parallel to {b}"),
            Self::BooleanFailed(err) => write!(f, "boolean operation failed: {err}"),
            Self::ModelingFailed(msg) => write!(f, "modeling operation failed: {msg}"),
            Self::InvalidProfile(msg) => write!(f, "invalid profile: {msg}"),
            Self::NoIntersection => write!(f, "feature does not intersect with target shape"),
            Self::InvalidDraftAngle { angle_rad } => {
                write!(f, "draft angle {:.1} degrees is out of valid range (-89, 89)", angle_rad.to_degrees())
            }
            Self::FaceNotFound { face_index } => write!(f, "face index {face_index} not found"),
            Self::InvalidNeutralPlane => write!(f, "neutral plane definition is invalid"),
        }
    }
}

impl std::error::Error for BRepFeatError {}

impl From<BooleanError> for BRepFeatError {
    fn from(value: BooleanError) -> Self {
        Self::BooleanFailed(value)
    }
}

impl From<rcad_modeling::BuildError> for BRepFeatError {
    fn from(value: rcad_modeling::BuildError) -> Self {
        Self::ModelingFailed(value.to_string())
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Fuse Mode and Parameters
// ═══════════════════════════════════════════════════════════════════════════════

/// Defines how a feature interacts with the base shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FuseMode {
    /// Add material to the base shape (boolean union).
    Add,
    /// Remove material from the base shape (boolean difference).
    Remove,
    /// Compute the intersection with the base shape.
    Common,
}

impl From<FuseMode> for BooleanOpType {
    fn from(mode: FuseMode) -> Self {
        match mode {
            FuseMode::Add => BooleanOpType::Union,
            FuseMode::Remove => BooleanOpType::Difference,
            FuseMode::Common => BooleanOpType::Intersection,
        }
    }
}

/// Parameters for feature operations.
#[derive(Debug, Clone)]
pub struct FeatureParams {
    /// Tolerance for merging vertices and edges after the operation.
    pub merge_tolerance: f64,
    /// Whether to perform validity checks after the operation.
    pub validate_result: bool,
    /// Whether to simplify the result (merge coplanar faces, etc.).
    pub simplify_result: bool,
}

impl Default for FeatureParams {
    fn default() -> Self {
        Self {
            merge_tolerance: 1e-6,
            validate_result: true,
            simplify_result: true,
        }
    }
}

/// Parameters specific to rib features.
#[derive(Debug, Clone)]
pub struct RibParams {
    /// Thickness of the rib (perpendicular to the profile).
    pub thickness: f64,
    /// Height of the rib (along the extrusion direction).
    pub height: f64,
    /// Draft angle for the rib sides (radians).
    pub draft_angle: f64,
    /// Whether to merge the rib with the target shape.
    pub fuse: bool,
}

impl Default for RibParams {
    fn default() -> Self {
        Self {
            thickness: 1.0,
            height: 1.0,
            draft_angle: 0.0,
            fuse: true,
        }
    }
}

/// Parameters specific to groove features.
#[derive(Debug, Clone)]
pub struct GrooveParams {
    /// Depth of the groove.
    pub depth: f64,
    /// Width of the groove (if applicable).
    pub width: Option<f64>,
    /// Whether the groove should go through the entire shape.
    pub through_all: bool,
    /// Draft angle for the groove sides (radians).
    pub draft_angle: f64,
}

impl Default for GrooveParams {
    fn default() -> Self {
        Self {
            depth: 1.0,
            width: None,
            through_all: false,
            draft_angle: 0.0,
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Helper Functions
// ═══════════════════════════════════════════════════════════════════════════════

const EPS: f64 = 1e-12;

fn validate_finite(name: &'static str, v: f64) -> Result<f64, BRepFeatError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(BRepFeatError::NonFiniteInput(name))
    }
}

fn validate_positive(name: &'static str, v: f64) -> Result<f64, BRepFeatError> {
    let v = validate_finite(name, v)?;
    if v > 0.0 {
        Ok(v)
    } else {
        Err(BRepFeatError::NonPositiveInput(name))
    }
}

fn normalize(name: &'static str, v: DVec3) -> Result<DVec3, BRepFeatError> {
    if !v.is_finite() {
        return Err(BRepFeatError::NonFiniteInput(name));
    }
    if v.length_squared() <= EPS {
        return Err(BRepFeatError::ZeroVector(name));
    }
    Ok(v.normalize())
}

fn axis_ref_basis(axis: DVec3, ref_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), BRepFeatError> {
    let y_axis = normalize("axis", axis)?;
    let ref_dir = normalize("ref_dir", ref_dir)?;
    let x_reject = ref_dir - y_axis * ref_dir.dot(y_axis);
    if x_reject.length_squared() <= EPS {
        return Err(BRepFeatError::ParallelVectors("ref_dir", "axis"));
    }
    let x_axis = x_reject.normalize();
    let z_axis = x_axis.cross(y_axis).normalize();
    Ok((x_axis, y_axis, z_axis))
}

/// Build a prism solid from bottom and top polygon sections.
fn build_prism_from_sections(bot: &[DVec3], top: &[DVec3], dir: DVec3) -> Result<BRep, BRepFeatError> {
    let n = bot.len();
    if n < 3 || top.len() != n {
        return Err(BRepFeatError::InvalidInput("section vertex count mismatch or too few vertices".to_string()));
    }

    let mut brep = BRep::new();

    // Add vertices: bot[0..n] then top[0..n]
    for &p in bot {
        brep.vertices.push(Vertex { point: p });
    }
    for &p in top {
        brep.vertices.push(Vertex { point: p });
    }

    fn add_line_edge(brep: &mut BRep, start: usize, end: usize) -> usize {
        let p0 = brep.vertices[start].point;
        let p1 = brep.vertices[end].point;
        let d = p1 - p0;
        let len = d.length();
        let dir_vec = if len > 0.0 { d / len } else { DVec3::X };
        let ei = brep.edges.len();
        brep.edges.push(Edge { start, end });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir_vec }));
        brep.geom.edge_curve.push(Some(ci));
        brep.geom.edge_curve_range.push(Some([0.0, len]));
        brep.geom.edge_degenerated.push(false);
        ei
    }

    // Bottom-cap edges
    let bot_edges: Vec<usize> = (0..n).map(|i| add_line_edge(&mut brep, i, (i + 1) % n)).collect();
    // Top-cap edges
    let top_edges: Vec<usize> = (0..n).map(|i| add_line_edge(&mut brep, n + i, n + (i + 1) % n)).collect();
    // Vertical edges
    let vert_edges: Vec<usize> = (0..n).map(|i| add_line_edge(&mut brep, i, n + i)).collect();

    let mut faces = Vec::with_capacity(n + 2);

    // Bottom cap (outward normal = -dir)
    {
        let wire_edges: Vec<WireEdge> = (0..n)
            .map(|i| WireEdge { idx: bot_edges[n - 1 - i], forward: false })
            .collect();
        faces.push(Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: vec![],
            normal: -dir,
            triangles: vec![],
            mesh_dirty: true,
        });
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane { origin: bot[0], normal: -dir }));
        brep.geom.face_surface.push(Some(si));
    }

    // Top cap (outward normal = +dir)
    {
        let wire_edges: Vec<WireEdge> = (0..n)
            .map(|i| WireEdge { idx: top_edges[i], forward: true })
            .collect();
        faces.push(Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: vec![],
            normal: dir,
            triangles: vec![],
            mesh_dirty: true,
        });
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane { origin: top[0], normal: dir }));
        brep.geom.face_surface.push(Some(si));
    }

    // Lateral quad faces
    for i in 0..n {
        let j = (i + 1) % n;
        let a = bot[i];
        let b = bot[j];
        let c = top[j];
        let face_normal = {
            let ab = b - a;
            let ac = c - a;
            let nv = ab.cross(ac);
            if nv.length_squared() > 1e-24 { nv.normalize() } else { -dir.cross(ab).normalize_or(DVec3::X) }
        };
        let wire_edges = vec![
            WireEdge { idx: bot_edges[i], forward: true },
            WireEdge { idx: vert_edges[j], forward: true },
            WireEdge { idx: top_edges[i], forward: false },
            WireEdge { idx: vert_edges[i], forward: false },
        ];
        faces.push(Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: vec![],
            normal: face_normal,
            triangles: vec![],
            mesh_dirty: true,
        });
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane { origin: a, normal: face_normal }));
        brep.geom.face_surface.push(Some(si));
    }

    brep.solids.push(Solid { shells: vec![Shell { faces }] });
    Ok(brep)
}

/// Build a polygon face BRep from vertices.
fn build_polygon_face_brep(profile_verts: &[DVec3]) -> Result<BRep, BRepFeatError> {
    if profile_verts.len() < 3 {
        return Err(BRepFeatError::InvalidInput("profile needs >= 3 vertices".to_string()));
    }

    let n = profile_verts.len();
    let mut brep = BRep::new();

    for &p in profile_verts {
        brep.vertices.push(Vertex { point: p });
    }

    for i in 0..n {
        let j = (i + 1) % n;
        brep.edges.push(Edge { start: i, end: j });
    }

    let normal = {
        let a = profile_verts[0];
        let b = profile_verts[1];
        let c = profile_verts[2];
        let n_vec = (b - a).cross(c - a);
        if n_vec.length_squared() <= EPS {
            return Err(BRepFeatError::InvalidInput("profile vertices are degenerate".to_string()));
        }
        n_vec.normalize()
    };

    let face = Face {
        outer_wire: Wire {
            edges: (0..n).map(WireEdge::fwd).collect(),
        },
        inner_wires: vec![],
        normal,
        triangles: vec![],
        mesh_dirty: true,
    };

    brep.solids.push(Solid {
        shells: vec![Shell { faces: vec![face] }],
    });

    Ok(brep)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Rib Operations
// ═══════════════════════════════════════════════════════════════════════════════

/// Create a rib feature from a wire profile.
///
/// A rib is a thin-wall feature that reinforces a part. The profile wire defines
/// the cross-section of the rib, and it is extruded in the given direction with
/// the specified thickness.
///
/// # Arguments
///
/// * `target` - The base shape to add the rib to.
/// * `profile_wire` - Vertices defining the rib profile (closed polygon).
/// * `direction` - Direction of rib extrusion.
/// * `thickness` - Thickness of the rib perpendicular to the profile.
///
/// # Returns
///
/// The resulting shape with the rib added.
///
/// # Example
///
/// ```ignore
/// let result = make_rib(&base_shape, &profile, DVec3::Y, 2.0)?;
/// ```
pub fn make_rib(
    target: &BRep,
    profile_wire: &[DVec3],
    direction: DVec3,
    thickness: f64,
) -> Result<BRep, BRepFeatError> {
    if profile_wire.len() < 3 {
        return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
    }

    let dir = normalize("direction", direction)?;
    let thickness = validate_positive("thickness", thickness)?;

    // Extrude the profile in both directions for thickness
    let half_thickness = thickness / 2.0;

    // Find the profile normal
    let a = profile_wire[0];
    let b = profile_wire[1];
    let c = profile_wire[2];
    let profile_normal = (b - a).cross(c - a).normalize();

    // Create two offset copies of the profile
    let profile_offset_neg: Vec<DVec3> = profile_wire.iter()
        .map(|&p| p - profile_normal * half_thickness)
        .collect();
    let profile_offset_pos: Vec<DVec3> = profile_wire.iter()
        .map(|&p| p + profile_normal * half_thickness)
        .collect();

    // Build a thick prism for the rib
    let rib_solid = build_rib_solid(&profile_offset_neg, &profile_offset_pos, dir, thickness)?;

    // Fuse with target
    Ok(boolean_op(BooleanOpType::Union, target, &rib_solid)?)
}

/// Create a linear rib from a profile with specified height.
///
/// Similar to `make_rib` but with explicit control over the rib height.
///
/// # Arguments
///
/// * `target` - The base shape to add the rib to.
/// * `profile` - Vertices defining the rib profile.
/// * `direction` - Direction of rib extrusion.
///
/// # Returns
///
/// The resulting shape with the linear rib added.
pub fn make_linear_rib(
    target: &BRep,
    profile: &[DVec3],
    direction: DVec3,
) -> Result<BRep, BRepFeatError> {
    if profile.len() < 3 {
        return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
    }

    let dir = normalize("direction", direction)?;

    // Compute profile centroid and height
    let centroid: DVec3 = profile.iter().sum::<DVec3>() / profile.len() as f64;

    // Find the profile extents along the direction
    let heights: Vec<f64> = profile.iter().map(|&p| (p - centroid).dot(dir)).collect();
    let min_h = heights.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max_h = heights.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let profile_height = max_h - min_h;

    // Create a prism from the profile
    let bottom: Vec<DVec3> = profile.iter().map(|&p| p - dir * min_h).collect();
    let top: Vec<DVec3> = bottom.iter().map(|&p| p + dir * profile_height).collect();

    let rib_solid = build_prism_from_sections(&bottom, &top, dir)?;

    // Fuse with target
    Ok(boolean_op(BooleanOpType::Union, target, &rib_solid)?)
}

/// Build a solid rib from two offset profile sections.
fn build_rib_solid(
    bot: &[DVec3],
    top: &[DVec3],
    dir: DVec3,
    _thickness: f64,
) -> Result<BRep, BRepFeatError> {
    build_prism_from_sections(bot, top, dir)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Groove Operations
// ═══════════════════════════════════════════════════════════════════════════════

/// Create a groove (slot) feature from a profile wire.
///
/// A groove is a depression cut into a part. The profile wire defines the
/// cross-section of the groove.
///
/// # Arguments
///
/// * `target` - The base shape to cut the groove into.
/// * `profile_wire` - Vertices defining the groove profile.
/// * `direction` - Direction of groove extrusion.
/// * `depth` - Depth of the groove.
///
/// # Returns
///
/// The resulting shape with the groove cut.
pub fn make_groove(
    target: &BRep,
    profile_wire: &[DVec3],
    direction: DVec3,
    depth: f64,
) -> Result<BRep, BRepFeatError> {
    if profile_wire.len() < 3 {
        return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
    }

    let dir = normalize("direction", direction)?;
    let depth = validate_positive("depth", depth)?;

    // Build the groove tool
    let bottom: Vec<DVec3> = profile_wire.to_vec();
    let top: Vec<DVec3> = bottom.iter().map(|&p| p + dir * depth).collect();

    let groove_tool = build_prism_from_sections(&bottom, &top, dir)?;

    // Subtract from target
    Ok(boolean_op(BooleanOpType::Difference, target, &groove_tool)?)
}

/// Create a through groove (slot) that goes through the entire shape.
///
/// Similar to `make_groove` but the groove extends through the entire target shape.
///
/// # Arguments
///
/// * `target` - The base shape to cut the groove into.
/// * `profile` - Vertices defining the groove profile.
/// * `direction` - Direction of groove extrusion.
///
/// # Returns
///
/// The resulting shape with the through groove cut.
pub fn make_through_groove(
    target: &BRep,
    profile: &[DVec3],
    direction: DVec3,
) -> Result<BRep, BRepFeatError> {
    if profile.len() < 3 {
        return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
    }

    let dir = normalize("direction", direction)?;

    // Compute the bounding box of the target to determine the through distance
    let (min_pt, max_pt) = compute_bounding_box(target);
    let extent = (max_pt - min_pt).length() + 1.0;

    // Build the groove tool that extends beyond the target
    let bottom: Vec<DVec3> = profile.iter().map(|&p| p - dir * extent).collect();
    let top: Vec<DVec3> = profile.iter().map(|&p| p + dir * extent).collect();

    let groove_tool = build_prism_from_sections(&bottom, &top, dir)?;

    // Subtract from target
    Ok(boolean_op(BooleanOpType::Difference, target, &groove_tool)?)
}

/// Compute the axis-aligned bounding box of a BRep.
fn compute_bounding_box(brep: &BRep) -> (DVec3, DVec3) {
    if brep.vertices.is_empty() {
        return (DVec3::ZERO, DVec3::ZERO);
    }

    let min_pt = brep.vertices.iter().fold(DVec3::INFINITY, |acc, v| acc.min(v.point));
    let max_pt = brep.vertices.iter().fold(DVec3::NEG_INFINITY, |acc, v| acc.max(v.point));

    (min_pt, max_pt)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Prism Feature
// ═══════════════════════════════════════════════════════════════════════════════

/// Create a feature-based prism from a profile.
///
/// Creates a prism by extruding a profile in the given direction and combining
/// it with the target shape using the specified fuse mode.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `profile` - Vertices defining the prism profile.
/// * `direction` - Direction of extrusion.
/// * `fuse_mode` - How to combine with the target (Add, Remove, Common).
///
/// # Returns
///
/// The resulting shape after the prism operation.
pub fn make_prism_feature(
    target: &BRep,
    profile: &[DVec3],
    direction: DVec3,
    fuse_mode: FuseMode,
) -> Result<BRep, BRepFeatError> {
    if profile.len() < 3 {
        return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
    }

    let dir = normalize("direction", direction)?;

    // Compute extrusion depth based on profile extents
    let centroid: DVec3 = profile.iter().sum::<DVec3>() / profile.len() as f64;
    let heights: Vec<f64> = profile.iter().map(|&p| (p - centroid).dot(dir)).collect();
    let max_h = heights.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    // Extrude the profile
    let bottom: Vec<DVec3> = profile.to_vec();
    let top: Vec<DVec3> = bottom.iter().map(|&p| p + dir * max_h.abs().max(1.0)).collect();

    let prism_tool = build_prism_from_sections(&bottom, &top, dir)?;

    // Apply boolean operation based on fuse mode
    let op = BooleanOpType::from(fuse_mode);
    Ok(boolean_op(op, target, &prism_tool)?)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Revolution Feature
// ═══════════════════════════════════════════════════════════════════════════════

/// Create a feature-based revolution from a profile.
///
/// Creates a revolution by rotating a profile around an axis and combining
/// it with the target shape using the specified fuse mode.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `profile` - Vertices defining the revolution profile.
/// * `axis` - Origin point of the revolution axis.
/// * `axis_dir` - Direction of the revolution axis.
/// * `angle` - Revolution angle in radians (full circle = 2*PI).
/// * `fuse_mode` - How to combine with the target (Add, Remove, Common).
///
/// # Returns
///
/// The resulting shape after the revolution operation.
pub fn make_revol_feature(
    target: &BRep,
    profile: &[DVec3],
    axis: DVec3,
    axis_dir: DVec3,
    angle: f64,
    fuse_mode: FuseMode,
) -> Result<BRep, BRepFeatError> {
    if profile.len() < 3 {
        return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
    }

    if !axis.is_finite() {
        return Err(BRepFeatError::NonFiniteInput("axis"));
    }

    let axis_dir = normalize("axis_dir", axis_dir)?;
    let angle = validate_finite("angle", angle)?;

    // Build the profile face
    let profile_brep = build_polygon_face_brep(profile)?;

    // Create the revolution
    let revol_tool = rcad_modeling::revolve(&profile_brep, 0, axis, axis_dir, angle)?;

    // Apply boolean operation based on fuse mode
    let op = BooleanOpType::from(fuse_mode);
    Ok(boolean_op(op, target, &revol_tool)?)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Pipe Feature
// ═══════════════════════════════════════════════════════════════════════════════

/// Create a feature-based pipe along a spine curve.
///
/// Creates a pipe by sweeping a profile along a spine path and combining
/// it with the target shape using the specified fuse mode.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `profile` - Vertices defining the pipe cross-section.
/// * `spine` - Vertices defining the spine path.
/// * `fuse_mode` - How to combine with the target (Add, Remove, Common).
///
/// # Returns
///
/// The resulting shape after the pipe operation.
pub fn make_pipe_feature(
    target: &BRep,
    profile: &[DVec3],
    spine: &[DVec3],
    fuse_mode: FuseMode,
) -> Result<BRep, BRepFeatError> {
    if profile.len() < 3 {
        return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
    }
    if spine.len() < 2 {
        return Err(BRepFeatError::InvalidInput("spine needs >= 2 vertices".to_string()));
    }

    // Build the pipe by sweeping the profile along the spine
    let pipe_tool = build_pipe_solid(profile, spine)?;

    // Apply boolean operation based on fuse mode
    let op = BooleanOpType::from(fuse_mode);
    Ok(boolean_op(op, target, &pipe_tool)?)
}

/// Build a pipe solid by sweeping a profile along a spine.
fn build_pipe_solid(profile: &[DVec3], spine: &[DVec3]) -> Result<BRep, BRepFeatError> {
    if spine.len() < 2 {
        return Err(BRepFeatError::InvalidInput("spine needs at least 2 points".to_string()));
    }

    // Compute the spine direction at each point
    let mut frames: Vec<(DVec3, DVec3, DVec3)> = Vec::with_capacity(spine.len());

    for i in 0..spine.len() {
        let tangent = if i == 0 {
            (spine[1] - spine[0]).normalize_or(DVec3::Z)
        } else if i == spine.len() - 1 {
            (spine[i] - spine[i - 1]).normalize_or(DVec3::Z)
        } else {
            (spine[i + 1] - spine[i - 1]).normalize_or(DVec3::Z)
        };

        // Build a local coordinate frame
        let up = if tangent.cross(DVec3::Y).length() > 0.1 {
            tangent.cross(DVec3::Y).normalize()
        } else {
            tangent.cross(DVec3::X).normalize()
        };
        let right = tangent.cross(up).normalize();

        frames.push((spine[i], right, up));
    }

    // Build cross-sections at each spine point
    let sections: Vec<Vec<DVec3>> = frames.iter().map(|(origin, right, up)| {
        profile.iter().map(|&p| {
            origin + right * p.x + up * p.y
        }).collect()
    }).collect();

    // Build the pipe by lofting through sections
    build_loft_solid(&sections)
}

/// Build a loft solid through multiple sections.
fn build_loft_solid(sections: &[Vec<DVec3>]) -> Result<BRep, BRepFeatError> {
    if sections.len() < 2 {
        return Err(BRepFeatError::InvalidInput("need at least 2 sections for loft".to_string()));
    }

    let n = sections[0].len();
    if n < 3 {
        return Err(BRepFeatError::InvalidInput("each section needs at least 3 vertices".to_string()));
    }

    // Check all sections have the same number of vertices
    for (i, s) in sections.iter().enumerate() {
        if s.len() != n {
            return Err(BRepFeatError::InvalidInput(format!(
                "section {} has {} vertices, expected {}", i, s.len(), n
            )));
        }
    }

    let mut brep = BRep::new();

    let num_sections = sections.len();

    // Add all vertices
    for si in 0..num_sections {
        for &p in &sections[si] {
            brep.vertices.push(Vertex { point: p });
        }
    }

    // Add edges and faces
    // Create edges between consecutive vertices in each section (caps)
    // and between sections (lateral faces)

    // Build edge tables
    // Cap edges for each section
    let mut cap_edges: Vec<Vec<usize>> = Vec::with_capacity(num_sections);
    for si in 0..num_sections {
        let base = si * n;
        let mut edges: Vec<usize> = Vec::with_capacity(n);
        for i in 0..n {
            let start = base + i;
            let end = base + (i + 1) % n;
            let p0 = brep.vertices[start].point;
            let p1 = brep.vertices[end].point;
            let d = p1 - p0;
            let len = d.length();
            let dir = if len > 0.0 { d / len } else { DVec3::X };
            let ei = brep.edges.len();
            brep.edges.push(Edge { start, end });
            let ci = brep.geom.curves.len();
            brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
            brep.geom.edge_curve.push(Some(ci));
            brep.geom.edge_curve_range.push(Some([0.0, len]));
            brep.geom.edge_degenerated.push(false);
            edges.push(ei);
        }
        cap_edges.push(edges);
    }

    // Lateral edges
    let mut lateral_edges: Vec<Vec<usize>> = Vec::with_capacity(num_sections - 1);
    for si in 0..num_sections - 1 {
        let base0 = si * n;
        let base1 = (si + 1) * n;
        let mut edges: Vec<usize> = Vec::with_capacity(n);
        for i in 0..n {
            let start = base0 + i;
            let end = base1 + i;
            let p0 = brep.vertices[start].point;
            let p1 = brep.vertices[end].point;
            let d = p1 - p0;
            let len = d.length();
            let dir = if len > 0.0 { d / len } else { DVec3::X };
            let ei = brep.edges.len();
            brep.edges.push(Edge { start, end });
            let ci = brep.geom.curves.len();
            brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
            brep.geom.edge_curve.push(Some(ci));
            brep.geom.edge_curve_range.push(Some([0.0, len]));
            brep.geom.edge_degenerated.push(false);
            edges.push(ei);
        }
        lateral_edges.push(edges);
    }

    // Build faces
    let mut faces = Vec::new();

    // Bottom cap
    let bottom_normal = {
        let a = sections[0][0];
        let b = sections[0][1];
        let c = sections[0][2];
        (b - a).cross(c - a).normalize_or(DVec3::Z)
    };
    let bottom_wire: Vec<WireEdge> = (0..n)
        .map(|i| WireEdge { idx: cap_edges[0][n - 1 - i], forward: false })
        .collect();
    faces.push(Face {
        outer_wire: Wire { edges: bottom_wire },
        inner_wires: vec![],
        normal: bottom_normal,
        triangles: vec![],
        mesh_dirty: true,
    });

    // Top cap
    let top_normal = {
        let a = sections[num_sections - 1][0];
        let b = sections[num_sections - 1][1];
        let c = sections[num_sections - 1][2];
        (b - a).cross(c - a).normalize_or(DVec3::Z)
    };
    let top_wire: Vec<WireEdge> = (0..n)
        .map(|i| WireEdge { idx: cap_edges[num_sections - 1][i], forward: true })
        .collect();
    faces.push(Face {
        outer_wire: Wire { edges: top_wire },
        inner_wires: vec![],
        normal: top_normal,
        triangles: vec![],
        mesh_dirty: true,
    });

    // Lateral faces (quads between sections)
    for si in 0..num_sections - 1 {
        for i in 0..n {
            let j = (i + 1) % n;
            let base0 = si * n;
            let base1 = (si + 1) * n;

            // Compute face normal
            let p0 = sections[si][i];
            let p1 = sections[si][j];
            let p2 = sections[si + 1][j];
            let normal = (p1 - p0).cross(p2 - p0).normalize_or(top_normal);

            let wire_edges = vec![
                WireEdge { idx: cap_edges[si][i], forward: true },
                WireEdge { idx: lateral_edges[si][j], forward: true },
                WireEdge { idx: cap_edges[si + 1][i], forward: false },
                WireEdge { idx: lateral_edges[si][i], forward: false },
            ];

            faces.push(Face {
                outer_wire: Wire { edges: wire_edges },
                inner_wires: vec![],
                normal,
                triangles: vec![],
                mesh_dirty: true,
            });
        }
    }

    brep.solids.push(Solid { shells: vec![Shell { faces }] });
    Ok(brep)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Draft Feature
// ═══════════════════════════════════════════════════════════════════════════════

/// Parameters for draft application.
#[derive(Debug, Clone)]
pub struct DraftFeatureParams {
    /// Draft angle in radians.
    pub angle: f64,
    /// Neutral plane origin point.
    pub neutral_point: DVec3,
    /// Neutral plane normal (pull direction).
    pub pull_direction: DVec3,
}

impl Default for DraftFeatureParams {
    fn default() -> Self {
        Self {
            angle: 2.0_f64.to_radians(),
            neutral_point: DVec3::ZERO,
            pull_direction: DVec3::Z,
        }
    }
}

/// Apply draft angle to specified faces.
///
/// Draft angle is applied to allow parts to be removed from molds. The draft
/// tilts the faces so that the part can be pulled out in the pull direction.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `face_indices` - Indices of faces to apply draft to.
/// * `angle` - Draft angle in radians. Positive values add draft for easy removal.
/// * `neutral_plane` - A point on the neutral plane (vertices here stay fixed).
///
/// # Returns
///
/// The shape with draft applied.
pub fn apply_draft_feature(
    target: &BRep,
    face_indices: &[usize],
    angle: f64,
    neutral_plane: DVec3,
) -> Result<BRep, BRepFeatError> {
    // Validate angle
    if angle.abs() >= std::f64::consts::FRAC_PI_2 - 1e-6 {
        return Err(BRepFeatError::InvalidDraftAngle { angle_rad: angle });
    }

    let shell = target.solids.first()
        .and_then(|s| s.shells.first())
        .ok_or_else(|| BRepFeatError::InvalidInput("target has no solids".to_string()))?;

    // Validate face indices
    for &fi in face_indices {
        if fi >= shell.faces.len() {
            return Err(BRepFeatError::FaceNotFound { face_index: fi });
        }
    }

    // If no faces specified, apply draft to all vertical faces
    let faces_to_draft: Vec<usize> = if face_indices.is_empty() {
        // Find all faces that are not horizontal (have a component perpendicular to Z)
        shell.faces.iter().enumerate()
            .filter(|(_, face)| {
                let normal = face.normal.normalize_or(DVec3::Z);
                normal.dot(DVec3::Z).abs() < 0.99
            })
            .map(|(i, _)| i)
            .collect()
    } else {
        face_indices.to_vec()
    };

    if faces_to_draft.is_empty() {
        // No faces to draft, return a clone
        return Ok(target.clone());
    }

    // Apply draft transformation
    let pull_dir = DVec3::Z; // Default pull direction
    let tan_angle = angle.tan();

    // Compute new vertex positions
    let mut new_positions: Vec<DVec3> = target.vertices.iter().map(|v| v.point).collect();

    for &fi in &faces_to_draft {
        let face = &shell.faces[fi];

        // Get vertices belonging to this face
        for we in &face.outer_wire.edges {
            if let Some(edge) = target.edges.get(we.idx) {
                for &vi in &[edge.start, edge.end] {
                    if vi < target.vertices.len() {
                        let v = target.vertices[vi].point;
                        let h = (v - neutral_plane).dot(pull_dir);
                        // Apply draft displacement
                        let radial_dir = (v - neutral_plane).reject_from(pull_dir).normalize_or(DVec3::ZERO);
                        if radial_dir.length() > 1e-10 {
                            new_positions[vi] = v + radial_dir * (h * tan_angle);
                        }
                    }
                }
            }
        }
    }

    // Build the new BRep with modified vertex positions
    build_drafted_brep(target, &new_positions, &shell.faces)
}

/// Build a new BRep with drafted vertex positions.
fn build_drafted_brep(
    original: &BRep,
    new_positions: &[DVec3],
    faces: &[Face],
) -> Result<BRep, BRepFeatError> {
    let mut brep = BRep::new();
    brep.solids.push(Solid { shells: vec![Shell { faces: Vec::new() }] });

    // Copy vertices with new positions
    for &p in new_positions {
        brep.vertices.push(Vertex { point: p });
    }

    // Copy edges with updated curves
    for e in original.edges.iter() {
        let p0 = new_positions.get(e.start).copied().unwrap_or(DVec3::ZERO);
        let p1 = new_positions.get(e.end).copied().unwrap_or(DVec3::ZERO);
        let d = p1 - p0;
        let len = d.length();
        let dir = if len > 0.0 { d / len } else { DVec3::X };

        brep.edges.push(Edge { start: e.start, end: e.end });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
        brep.geom.edge_curve.push(Some(ci));
        brep.geom.edge_curve_range.push(Some([0.0, len]));
        brep.geom.edge_degenerated.push(false);

        // Preserve edge mapping
        while brep.edges.len() > brep.geom.edge_curve.len() {
            brep.geom.edge_curve.push(None);
            brep.geom.edge_curve_range.push(None);
            brep.geom.edge_degenerated.push(false);
        }
    }

    // Copy faces
    for face in faces {
        let wire_edges: Vec<WireEdge> = face.outer_wire.edges.iter().map(|we| {
            WireEdge { idx: we.idx, forward: we.forward }
        }).collect();

        brep.solids[0].shells[0].faces.push(Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: face.inner_wires.clone(),
            normal: face.normal,
            triangles: vec![],
            mesh_dirty: true,
        });
    }

    Ok(brep)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Advanced Feature Operations
// ═══════════════════════════════════════════════════════════════════════════════

/// Create a drafted prism feature (tapered extrusion).
///
/// Creates a prism with draft angle, useful for molded parts.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `profile` - Vertices defining the prism profile.
/// * `direction` - Direction of extrusion.
/// * `depth` - Extrusion depth.
/// * `draft_angle` - Draft angle in radians.
/// * `fuse_mode` - How to combine with the target.
///
/// # Returns
///
/// The resulting shape with the drafted prism.
pub fn make_drafted_prism(
    target: &BRep,
    profile: &[DVec3],
    direction: DVec3,
    depth: f64,
    draft_angle: f64,
    fuse_mode: FuseMode,
) -> Result<BRep, BRepFeatError> {
    if profile.len() < 3 {
        return Err(BRepFeatError::InvalidProfile("profile needs >= 3 vertices".to_string()));
    }

    let dir = normalize("direction", direction)?;
    let depth = validate_positive("depth", depth)?;

    if draft_angle.abs() >= std::f64::consts::FRAC_PI_2 - 1e-6 {
        return Err(BRepFeatError::InvalidDraftAngle { angle_rad: draft_angle });
    }

    // Compute centroid
    let centroid: DVec3 = profile.iter().sum::<DVec3>() / profile.len() as f64;

    // Apply draft by scaling the top profile
    let taper = depth * draft_angle.tan();

    let bottom: Vec<DVec3> = profile.to_vec();
    let top: Vec<DVec3> = profile.iter().map(|&p| {
        let radial = p - centroid;
        let radial_2d = radial - dir * radial.dot(dir);
        let radial_dir = if radial_2d.length() > EPS {
            radial_2d.normalize()
        } else {
            DVec3::ZERO
        };
        p + dir * depth + radial_dir * taper
    }).collect();

    let prism_tool = build_prism_from_sections(&bottom, &top, dir)?;

    let op = BooleanOpType::from(fuse_mode);
    Ok(boolean_op(op, target, &prism_tool)?)
}

/// Create a multi-profile pipe (loft) feature.
///
/// Creates a solid by lofting through multiple profiles.
///
/// # Arguments
///
/// * `target` - The base shape.
/// * `profiles` - Vector of profiles, each being a vector of vertices.
/// * `fuse_mode` - How to combine with the target.
///
/// # Returns
///
/// The resulting shape with the loft feature.
pub fn make_loft_feature(
    target: &BRep,
    profiles: &[Vec<DVec3>],
    fuse_mode: FuseMode,
) -> Result<BRep, BRepFeatError> {
    if profiles.len() < 2 {
        return Err(BRepFeatError::InvalidInput("need at least 2 profiles for loft".to_string()));
    }

    for (i, profile) in profiles.iter().enumerate() {
        if profile.len() < 3 {
            return Err(BRepFeatError::InvalidProfile(format!(
                "profile {} needs >= 3 vertices", i
            )));
        }
    }

    let loft_tool = build_loft_solid(profiles)?;

    let op = BooleanOpType::from(fuse_mode);
    Ok(boolean_op(op, target, &loft_tool)?)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    fn make_test_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        })
    }

    #[test]
    fn test_make_rib_adds_material() {
        let target = make_test_box();

        let profile = vec![
            DVec3::new(-0.5, 2.0, -0.5),
            DVec3::new(0.5, 2.0, -0.5),
            DVec3::new(0.5, 2.0, 0.5),
            DVec3::new(-0.5, 2.0, 0.5),
        ];

        let result = make_rib(&target, &profile, DVec3::Y, 1.0);

        assert!(result.is_ok(), "make_rib should succeed: {:?}", result);
        let result = result.unwrap();
        assert!(!result.solids.is_empty(), "result should have solids");
    }

    #[test]
    fn test_make_linear_rib_adds_material() {
        let target = make_test_box();

        let profile = vec![
            DVec3::new(-0.5, 2.0, -0.5),
            DVec3::new(0.5, 2.0, -0.5),
            DVec3::new(0.5, 2.0, 0.5),
            DVec3::new(-0.5, 2.0, 0.5),
        ];

        let result = make_linear_rib(&target, &profile, DVec3::Y);

        assert!(result.is_ok(), "make_linear_rib should succeed: {:?}", result);
        let result = result.unwrap();
        assert!(!result.solids.is_empty(), "result should have solids");
    }

    #[test]
    fn test_make_groove_removes_material() {
        let target = make_test_box();

        let profile = vec![
            DVec3::new(-0.5, 0.0, -0.5),
            DVec3::new(0.5, 0.0, -0.5),
            DVec3::new(0.5, 0.0, 0.5),
            DVec3::new(-0.5, 0.0, 0.5),
        ];

        let result = make_groove(&target, &profile, DVec3::Y, 2.0);

        assert!(result.is_ok(), "make_groove should succeed: {:?}", result);
        let result = result.unwrap();
        assert!(!result.solids.is_empty(), "result should have solids");
    }

    #[test]
    fn test_make_through_groove() {
        let target = make_test_box();

        let profile = vec![
            DVec3::new(-0.3, 0.0, -0.3),
            DVec3::new(0.3, 0.0, -0.3),
            DVec3::new(0.3, 0.0, 0.3),
            DVec3::new(-0.3, 0.0, 0.3),
        ];

        let result = make_through_groove(&target, &profile, DVec3::Y);

        assert!(result.is_ok(), "make_through_groove should succeed: {:?}", result);
    }

    #[test]
    fn test_make_prism_feature_add() {
        let target = make_test_box();

        let profile = vec![
            DVec3::new(-0.5, 2.0, -0.5),
            DVec3::new(0.5, 2.0, -0.5),
            DVec3::new(0.5, 2.0, 0.5),
            DVec3::new(-0.5, 2.0, 0.5),
        ];

        let result = make_prism_feature(&target, &profile, DVec3::Y, FuseMode::Add);

        assert!(result.is_ok(), "make_prism_feature Add should succeed: {:?}", result);
    }

    #[test]
    fn test_make_prism_feature_remove() {
        let target = make_test_box();

        let profile = vec![
            DVec3::new(-0.5, 0.0, -0.5),
            DVec3::new(0.5, 0.0, -0.5),
            DVec3::new(0.5, 0.0, 0.5),
            DVec3::new(-0.5, 0.0, 0.5),
        ];

        let result = make_prism_feature(&target, &profile, DVec3::Y, FuseMode::Remove);

        assert!(result.is_ok(), "make_prism_feature Remove should succeed: {:?}", result);
    }

    #[test]
    fn test_make_revol_feature() {
        let target = make_test_box();

        let profile = vec![
            DVec3::new(1.5, -0.3, 0.0),
            DVec3::new(2.0, -0.3, 0.0),
            DVec3::new(2.0, 0.3, 0.0),
            DVec3::new(1.5, 0.3, 0.0),
        ];

        let result = make_revol_feature(
            &target,
            &profile,
            DVec3::ZERO,
            DVec3::Z,
            std::f64::consts::TAU,
            FuseMode::Add,
        );

        assert!(result.is_ok(), "make_revol_feature should succeed: {:?}", result);
    }

    #[test]
    fn test_make_pipe_feature() {
        let target = make_test_box();

        let profile = vec![
            DVec3::new(0.0, -0.3, -0.3),
            DVec3::new(0.0, 0.3, -0.3),
            DVec3::new(0.0, 0.3, 0.3),
            DVec3::new(0.0, -0.3, 0.3),
        ];

        let spine = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(0.0, 2.0, 0.0),
            DVec3::new(0.0, 4.0, 0.0),
        ];

        let result = make_pipe_feature(&target, &profile, &spine, FuseMode::Add);

        assert!(result.is_ok(), "make_pipe_feature should succeed: {:?}", result);
    }

    #[test]
    fn test_apply_draft_feature() {
        let target = make_test_box();

        // Apply draft to all vertical faces (indices 0-5 are the 6 faces of a box)
        let result = apply_draft_feature(&target, &[0, 1, 2, 3], 5.0_f64.to_radians(), DVec3::ZERO);

        assert!(result.is_ok(), "apply_draft_feature should succeed: {:?}", result);
    }

    #[test]
    fn test_make_drafted_prism() {
        let target = make_test_box();

        let profile = vec![
            DVec3::new(-0.5, 2.0, -0.5),
            DVec3::new(0.5, 2.0, -0.5),
            DVec3::new(0.5, 2.0, 0.5),
            DVec3::new(-0.5, 2.0, 0.5),
        ];

        let result = make_drafted_prism(
            &target,
            &profile,
            DVec3::Y,
            1.0,
            5.0_f64.to_radians(),
            FuseMode::Add,
        );

        assert!(result.is_ok(), "make_drafted_prism should succeed: {:?}", result);
    }

    #[test]
    fn test_make_loft_feature() {
        let target = make_test_box();

        let profile1 = vec![
            DVec3::new(-0.5, 0.0, -0.5),
            DVec3::new(0.5, 0.0, -0.5),
            DVec3::new(0.5, 0.0, 0.5),
            DVec3::new(-0.5, 0.0, 0.5),
        ];

        let profile2 = vec![
            DVec3::new(-1.0, 2.0, -1.0),
            DVec3::new(1.0, 2.0, -1.0),
            DVec3::new(1.0, 2.0, 1.0),
            DVec3::new(-1.0, 2.0, 1.0),
        ];

        let profiles = vec![profile1, profile2];

        let result = make_loft_feature(&target, &profiles, FuseMode::Add);

        assert!(result.is_ok(), "make_loft_feature should succeed: {:?}", result);
    }

    #[test]
    fn test_invalid_profile_rejected() {
        let target = make_test_box();

        let profile = vec![DVec3::ZERO, DVec3::X];

        let result = make_rib(&target, &profile, DVec3::Y, 1.0);

        assert!(result.is_err(), "make_rib should reject profile with < 3 vertices");
    }

    #[test]
    fn test_zero_direction_rejected() {
        let target = make_test_box();

        let profile = vec![
            DVec3::new(-0.5, 0.0, -0.5),
            DVec3::new(0.5, 0.0, -0.5),
            DVec3::new(0.5, 0.0, 0.5),
        ];

        let result = make_groove(&target, &profile, DVec3::ZERO, 1.0);

        assert!(result.is_err(), "make_groove should reject zero direction");
    }

    #[test]
    fn test_invalid_draft_angle_rejected() {
        let target = make_test_box();

        let result = apply_draft_feature(&target, &[], std::f64::consts::FRAC_PI_2, DVec3::ZERO);

        assert!(result.is_err(), "apply_draft_feature should reject 90 degree angle");
    }

    #[test]
    fn test_fuse_mode_conversion() {
        assert_eq!(BooleanOpType::from(FuseMode::Add), BooleanOpType::Union);
        assert_eq!(BooleanOpType::from(FuseMode::Remove), BooleanOpType::Difference);
        assert_eq!(BooleanOpType::from(FuseMode::Common), BooleanOpType::Intersection);
    }

    #[test]
    fn test_feature_params_default() {
        let params = FeatureParams::default();

        assert!(params.merge_tolerance > 0.0);
        assert!(params.validate_result);
        assert!(params.simplify_result);
    }

    #[test]
    fn test_rib_params_default() {
        let params = RibParams::default();

        assert!(params.thickness > 0.0);
        assert!(params.height > 0.0);
        assert_eq!(params.draft_angle, 0.0);
        assert!(params.fuse);
    }

    #[test]
    fn test_groove_params_default() {
        let params = GrooveParams::default();

        assert!(params.depth > 0.0);
        assert!(params.width.is_none());
        assert!(!params.through_all);
    }

    #[test]
    fn test_error_display() {
        let err = BRepFeatError::InvalidDraftAngle { angle_rad: 1.5 };
        let s = format!("{}", err);
        assert!(s.contains("degrees"));

        let err = BRepFeatError::InvalidProfile("test".to_string());
        let s = format!("{}", err);
        assert!(s.contains("test"));

        let err = BRepFeatError::FaceNotFound { face_index: 99 };
        let s = format!("{}", err);
        assert!(s.contains("99"));
    }
}
