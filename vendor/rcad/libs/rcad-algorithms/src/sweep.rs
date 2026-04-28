//! BRepSweep-style sweep operations for creating shapes by sweeping profiles.
//!
//! This module provides comprehensive sweep operations analogous to OCCT's BRepSweep:
//! - **LinearSweep**: Extrude profiles along a direction
//! - **RotationalSweep**: Revolve profiles around an axis
//! - **PipeSweep**: Sweep profiles along a spine path with Frenet frame support
//!
//! # Example
//!
//! ```
//! use glam::DVec3;
//! use rcad_algorithms::sweep::{linear_sweep, SweepOptions, SweepMode};
//!
//! // Create a simple box by extruding a square profile
//! let profile_pts = vec![
//!     DVec3::new(0.0, 0.0, 0.0),
//!     DVec3::new(1.0, 0.0, 0.0),
//!     DVec3::new(1.0, 1.0, 0.0),
//!     DVec3::new(0.0, 1.0, 0.0),
//! ];
//! let result = linear_sweep(&profile_pts, DVec3::Z, 2.0);
//! ```

use std::collections::HashMap;

use glam::{DVec2, DVec3};
use rcad_kernel::BRep;
use rcad_kernel::geom::{
    Circle3, Curve3, Line3, Plane, Surface3,
};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

// ─────────────────────────────────────────────────────────────────────────────
// Error Types
// ─────────────────────────────────────────────────────────────────────────────

/// Errors returned by sweep operations.
#[derive(Debug, Clone)]
pub enum SweepError {
    /// Input vector is zero-length.
    ZeroVector(&'static str),
    /// Input value is not finite (NaN or infinity).
    NonFiniteInput(&'static str),
    /// Input value must be positive.
    NonPositiveInput(&'static str),
    /// Profile has insufficient vertices.
    InsufficientVertices { minimum: usize, actual: usize },
    /// Spine path has insufficient points.
    InsufficientSpinePoints { minimum: usize, actual: usize },
    /// Profile vertex counts don't match for multi-section sweep.
    VertexCountMismatch { expected: usize, actual: usize },
    /// Degenerate geometry encountered.
    DegenerateGeometry(&'static str),
    /// Invalid parameter value.
    InvalidParameter(&'static str),
    /// Corner handling failed.
    CornerHandlingFailed(String),
    /// Modeling operation failed.
    ModelingError(String),
}

impl std::fmt::Display for SweepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ZeroVector(name) => write!(f, "{} must be a non-zero vector", name),
            Self::NonFiniteInput(name) => write!(f, "{} must be finite", name),
            Self::NonPositiveInput(name) => write!(f, "{} must be positive", name),
            Self::InsufficientVertices { minimum, actual } => {
                write!(f, "profile needs at least {} vertices, got {}", minimum, actual)
            }
            Self::InsufficientSpinePoints { minimum, actual } => {
                write!(f, "spine needs at least {} points, got {}", minimum, actual)
            }
            Self::VertexCountMismatch { expected, actual } => {
                write!(f, "profile vertex count mismatch: expected {}, got {}", expected, actual)
            }
            Self::DegenerateGeometry(msg) => write!(f, "degenerate geometry: {}", msg),
            Self::InvalidParameter(msg) => write!(f, "invalid parameter: {}", msg),
            Self::CornerHandlingFailed(msg) => write!(f, "corner handling failed: {}", msg),
            Self::ModelingError(msg) => write!(f, "modeling error: {}", msg),
        }
    }
}

impl std::error::Error for SweepError {}

impl From<rcad_modeling::BuildError> for SweepError {
    fn from(value: rcad_modeling::BuildError) -> Self {
        Self::ModelingError(value.to_string())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Sweep Configuration Types
// ─────────────────────────────────────────────────────────────────────────────

/// Sweep mode determining how the profile is transformed along the path.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SweepMode {
    /// Simple translation along a direction (linear extrusion).
    #[default]
    Translation,
    /// Rotation around an axis (revolution).
    Rotation,
    /// Pipe sweep along a spine curve.
    Pipe,
}

/// Corner handling mode for pipe sweeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CornerMode {
    /// Sharp corners with no rounding.
    #[default]
    Sharp,
    /// Round corners with specified radius.
    Rounded,
    /// Extend and intersect at corners.
    Extended,
}

/// Configuration options for sweep operations.
#[derive(Debug, Clone)]
pub struct SweepOptions {
    /// Sweep mode.
    pub mode: SweepMode,
    /// Use Frenet frame for pipe sweeps (tangent, normal, binormal).
    /// When false, uses a fixed reference direction.
    pub is_frenet: bool,
    /// How to handle corners in pipe sweeps.
    pub corner_mode: CornerMode,
    /// Corner radius for rounded corners (only used when corner_mode is Rounded).
    pub corner_radius: f64,
    /// Maintain continuous normal across segments.
    pub continuous_normal: bool,
    /// Close the sweep (connect last section back to first).
    pub closed: bool,
    /// Tolerance for geometric operations.
    pub tolerance: f64,
}

impl Default for SweepOptions {
    fn default() -> Self {
        Self {
            mode: SweepMode::Translation,
            is_frenet: true,
            corner_mode: CornerMode::Sharp,
            corner_radius: 0.0,
            continuous_normal: false,
            closed: false,
            tolerance: 1e-6,
        }
    }
}

impl SweepOptions {
    /// Create options for a linear sweep.
    pub fn linear() -> Self {
        Self {
            mode: SweepMode::Translation,
            ..Default::default()
        }
    }

    /// Create options for a rotational sweep.
    pub fn rotational() -> Self {
        Self {
            mode: SweepMode::Rotation,
            ..Default::default()
        }
    }

    /// Create options for a pipe sweep with Frenet frame.
    pub fn pipe_frenet() -> Self {
        Self {
            mode: SweepMode::Pipe,
            is_frenet: true,
            ..Default::default()
        }
    }

    /// Create options for a pipe sweep with fixed reference direction.
    pub fn pipe_fixed() -> Self {
        Self {
            mode: SweepMode::Pipe,
            is_frenet: false,
            ..Default::default()
        }
    }

    /// Set the corner mode.
    pub fn with_corner_mode(mut self, mode: CornerMode) -> Self {
        self.corner_mode = mode;
        self
    }

    /// Set the corner radius.
    pub fn with_corner_radius(mut self, radius: f64) -> Self {
        self.corner_radius = radius;
        self
    }

    /// Enable/disable continuous normal.
    pub fn with_continuous_normal(mut self, enabled: bool) -> Self {
        self.continuous_normal = enabled;
        self
    }

    /// Enable/disable closed sweep.
    pub fn with_closed(mut self, closed: bool) -> Self {
        self.closed = closed;
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// History Tracking
// ─────────────────────────────────────────────────────────────────────────────

/// Tracks which result faces came from which parts of the sweep operation.
#[derive(Debug, Clone, Default)]
pub struct SweepHistory {
    /// Result face indices that are the start cap.
    pub start_cap: Vec<usize>,
    /// Result face indices that are the end cap.
    pub end_cap: Vec<usize>,
    /// Result face indices that are lateral swept faces.
    pub lateral_faces: Vec<usize>,
    /// Mapping from profile edge index to generated lateral face index.
    pub profile_edge_to_lateral: HashMap<usize, usize>,
    /// Corner face indices (for rounded corners).
    pub corner_faces: Vec<usize>,
}

impl SweepHistory {
    /// All result face indices tracked by this history.
    pub fn all_faces(&self) -> impl Iterator<Item = usize> + '_ {
        self.start_cap
            .iter()
            .chain(&self.end_cap)
            .chain(&self.lateral_faces)
            .chain(&self.corner_faces)
            .copied()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal Helpers
// ─────────────────────────────────────────────────────────────────────────────

const EPS: f64 = 1e-12;

fn validate_finite(name: &'static str, v: f64) -> Result<f64, SweepError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(SweepError::NonFiniteInput(name))
    }
}

fn validate_positive(name: &'static str, v: f64) -> Result<f64, SweepError> {
    let v = validate_finite(name, v)?;
    if v > 0.0 {
        Ok(v)
    } else {
        Err(SweepError::NonPositiveInput(name))
    }
}

fn normalize_vector(name: &'static str, v: DVec3) -> Result<DVec3, SweepError> {
    if !v.is_finite() {
        return Err(SweepError::NonFiniteInput(name));
    }
    let len_sq = v.length_squared();
    if len_sq <= EPS {
        return Err(SweepError::ZeroVector(name));
    }
    Ok(v / len_sq.sqrt())
}

/// Rodrigues' rotation formula: rotate point `p` around axis through `axis_origin`.
fn rotate_point(p: DVec3, axis_origin: DVec3, axis_dir: DVec3, angle: f64) -> DVec3 {
    let v = p - axis_origin;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let rotated = v * cos_a + axis_dir.cross(v) * sin_a + axis_dir * axis_dir.dot(v) * (1.0 - cos_a);
    rotated + axis_origin
}

/// Compute the outward normal of a planar polygon from first 3 vertices.
fn polygon_normal(pts: &[DVec3]) -> DVec3 {
    if pts.len() < 3 {
        return DVec3::Z;
    }
    let n = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
    if n.length_squared() > EPS {
        n.normalize()
    } else {
        DVec3::Z
    }
}

/// Compute the centroid of a set of points.
fn centroid(pts: &[DVec3]) -> DVec3 {
    if pts.is_empty() {
        return DVec3::ZERO;
    }
    pts.iter().sum::<DVec3>() / pts.len() as f64
}

// ─────────────────────────────────────────────────────────────────────────────
// BRep Builder Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Internal helper to add a vertex and return its index.
fn add_vertex(brep: &mut BRep, point: DVec3) -> usize {
    let idx = brep.vertices.len();
    brep.vertices.push(Vertex { point });
    idx
}

/// Internal helper to add a linear edge and return its index.
fn add_line_edge(brep: &mut BRep, start: usize, end: usize) -> usize {
    let p0 = brep.vertices[start].point;
    let p1 = brep.vertices[end].point;
    let d = p1 - p0;
    let len = d.length();
    let dir = if len > EPS { d / len } else { DVec3::X };

    let ei = brep.edges.len();
    brep.edges.push(Edge { start, end });

    let ci = brep.geom.curves.len();
    brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
    brep.geom.edge_curve.push(Some(ci));
    brep.geom.edge_curve_range.push(Some([0.0, len]));
    brep.geom.edge_degenerated.push(false);
    ei
}

/// Internal helper to add a circular arc edge.
fn add_arc_edge(brep: &mut BRep, circle: Circle3, start_angle: f64, end_angle: f64,
                start_v: usize, end_v: usize) -> usize {
    let ei = brep.edges.len();
    brep.edges.push(Edge { start: start_v, end: end_v });

    let ci = brep.geom.curves.len();
    brep.geom.curves.push(Curve3::Circle(circle));
    brep.geom.edge_curve.push(Some(ci));
    brep.geom.edge_curve_range.push(Some([start_angle, end_angle]));
    brep.geom.edge_degenerated.push(false);
    ei
}

// ─────────────────────────────────────────────────────────────────────────────
// Linear Sweep
// ─────────────────────────────────────────────────────────────────────────────

/// Extrude a profile polygon along a direction by a distance.
///
/// Creates a solid by sweeping a closed polygon profile along a linear path.
/// The result is a closed shell with:
/// - Start cap (original profile)
/// - End cap (translated copy)
/// - Lateral faces (one per profile edge)
///
/// # Arguments
/// * `profile_pts` - Polygon vertices in CCW order when viewed from +normal direction
/// * `direction` - Extrusion direction (will be normalized)
/// * `distance` - Extrusion distance (must be positive)
///
/// # Returns
/// A closed BRep solid.
pub fn linear_sweep(
    profile_pts: &[DVec3],
    direction: DVec3,
    distance: f64,
) -> Result<BRep, SweepError> {
    let (brep, _) = linear_sweep_with_history(profile_pts, direction, distance)?;
    Ok(brep)
}

/// Extrude a profile polygon along a direction with history tracking.
pub fn linear_sweep_with_history(
    profile_pts: &[DVec3],
    direction: DVec3,
    distance: f64,
) -> Result<(BRep, SweepHistory), SweepError> {
    linear_sweep_with_options(profile_pts, direction, distance, SweepOptions::default())
}

/// Extrude a profile polygon along a direction with custom options.
pub fn linear_sweep_with_options(
    profile_pts: &[DVec3],
    direction: DVec3,
    distance: f64,
    _options: SweepOptions,
) -> Result<(BRep, SweepHistory), SweepError> {
    if profile_pts.len() < 3 {
        return Err(SweepError::InsufficientVertices {
            minimum: 3,
            actual: profile_pts.len(),
        });
    }

    let dir = normalize_vector("direction", direction)?;
    let distance = validate_positive("distance", distance)?;

    let n = profile_pts.len();
    let offset = dir * distance;

    // Create the end profile points
    let end_pts: Vec<DVec3> = profile_pts.iter().map(|&p| p + offset).collect();

    // Compute profile normal
    let profile_normal = polygon_normal(profile_pts);

    let mut brep = BRep {
        vertices: Vec::with_capacity(2 * n),
        edges: Vec::new(),
        solids: Vec::new(),
        geom: rcad_kernel::GeomStore::default(),
        compound: None,
        compsolid: None,
    };

    // Add vertices: start[0..n], then end[0..n]
    let start_vi: Vec<usize> = profile_pts.iter()
        .map(|&p| add_vertex(&mut brep, p))
        .collect();
    let end_vi: Vec<usize> = end_pts.iter()
        .map(|&p| add_vertex(&mut brep, p))
        .collect();

    let mut faces = Vec::with_capacity(n + 2);
    let mut history = SweepHistory::default();

    // Start cap (normal pointing outward = -dir for solid interior)
    {
        let mut wire_edges = Vec::with_capacity(n);
        for i in 0..n {
            let j = (i + 1) % n;
            let ei = add_line_edge(&mut brep, start_vi[i], start_vi[j]);
            wire_edges.push(WireEdge::fwd(ei));
        }
        let cap_normal = -dir;
        let face = Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: vec![],
            normal: cap_normal,
            triangles: vec![],
            mesh_dirty: true,
        };
        let fi = faces.len();
        faces.push(face);
        history.start_cap.push(fi);

        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: profile_pts[0],
            normal: cap_normal,
        }));
        brep.geom.face_surface.push(Some(si));
    }

    // End cap (normal pointing outward = +dir)
    {
        let mut wire_edges = Vec::with_capacity(n);
        for i in 0..n {
            let j = (i + 1) % n;
            let ei = add_line_edge(&mut brep, end_vi[i], end_vi[j]);
            // Reverse for correct winding (CCW from +dir view)
            wire_edges.push(WireEdge::rev(ei));
        }
        let face = Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: vec![],
            normal: dir,
            triangles: vec![],
            mesh_dirty: true,
        };
        let fi = faces.len();
        faces.push(face);
        history.end_cap.push(fi);

        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: end_pts[0],
            normal: dir,
        }));
        brep.geom.face_surface.push(Some(si));
    }

    // Lateral faces: one quad per profile edge
    for i in 0..n {
        let j = (i + 1) % n;

        let a = profile_pts[i];
        let b = profile_pts[j];
        let c = end_pts[j];

        // Compute quad normal (outward)
        let lat_normal = (b - a).cross(c - a).normalize_or(profile_normal);

        // Create 4 edges for the quad
        let e_bot = add_line_edge(&mut brep, start_vi[i], start_vi[j]);
        let e_right = add_line_edge(&mut brep, start_vi[j], end_vi[j]);
        let e_top = add_line_edge(&mut brep, end_vi[j], end_vi[i]);
        let e_left = add_line_edge(&mut brep, end_vi[i], start_vi[i]);

        let wire_edges = vec![
            WireEdge::fwd(e_bot),
            WireEdge::fwd(e_right),
            WireEdge::fwd(e_top),
            WireEdge::fwd(e_left),
        ];

        let face = Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: vec![],
            normal: lat_normal,
            triangles: vec![],
            mesh_dirty: true,
        };
        let fi = faces.len();
        faces.push(face);
        history.lateral_faces.push(fi);
        history.profile_edge_to_lateral.insert(i, fi);

        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: a,
            normal: lat_normal,
        }));
        brep.geom.face_surface.push(Some(si));
    }

    brep.solids.push(Solid {
        shells: vec![Shell { faces }],
    });

    Ok((brep, history))
}

/// Extrude an existing BRep face along a direction.
///
/// # Arguments
/// * `brep` - Source BRep containing the face to extrude
/// * `face_idx` - Index of the face in the BRep
/// * `direction` - Extrusion direction
/// * `distance` - Extrusion distance
pub fn linear_sweep_face(
    brep: &BRep,
    face_idx: usize,
    direction: DVec3,
    distance: f64,
) -> Result<BRep, SweepError> {
    // Extract the face boundary points
    let face = brep.solids.first()
        .and_then(|s| s.shells.first())
        .and_then(|sh| sh.faces.get(face_idx))
        .ok_or_else(|| SweepError::InvalidParameter("face_idx out of range"))?;

    let pts: Vec<DVec3> = face.outer_wire.edges.iter()
        .filter_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            let vi = if we.forward { edge.start } else { edge.end };
            brep.vertices.get(vi).map(|v| v.point)
        })
        .collect();

    linear_sweep(&pts, direction, distance)
}

/// Extrude an existing BRep wire (closed wire becomes a shell).
///
/// # Arguments
/// * `brep` - Source BRep containing the wire
/// * `wire_idx` - Index identifying which wire (0 = outer wire of face 0)
/// * `direction` - Extrusion direction
/// * `distance` - Extrusion distance
pub fn linear_sweep_wire(
    brep: &BRep,
    wire_idx: usize,
    direction: DVec3,
    distance: f64,
) -> Result<BRep, SweepError> {
    // For simplicity, extract from face's outer wire
    let face = brep.solids.first()
        .and_then(|s| s.shells.first())
        .and_then(|sh| sh.faces.get(wire_idx))
        .ok_or_else(|| SweepError::InvalidParameter("wire_idx out of range"))?;

    let pts: Vec<DVec3> = face.outer_wire.edges.iter()
        .filter_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            let vi = if we.forward { edge.start } else { edge.end };
            brep.vertices.get(vi).map(|v| v.point)
        })
        .collect();

    linear_sweep(&pts, direction, distance)
}

// ─────────────────────────────────────────────────────────────────────────────
// Rotational Sweep
// ─────────────────────────────────────────────────────────────────────────────

/// Revolve a profile polygon around an axis.
///
/// Creates a solid by rotating a closed polygon profile around an axis.
/// For a full revolution (angle = 2*pi), the start and end are identified.
///
/// # Arguments
/// * `profile_pts` - Polygon vertices defining the profile
/// * `axis_origin` - A point on the rotation axis
/// * `axis_dir` - Direction of the rotation axis
/// * `angle` - Rotation angle in radians (positive = CCW when viewed from +axis)
///
/// # Returns
/// A closed BRep solid.
pub fn rotational_sweep(
    profile_pts: &[DVec3],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle: f64,
) -> Result<BRep, SweepError> {
    let (brep, _) = rotational_sweep_with_history(profile_pts, axis_origin, axis_dir, angle)?;
    Ok(brep)
}

/// Revolve a profile polygon around an axis with history tracking.
pub fn rotational_sweep_with_history(
    profile_pts: &[DVec3],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle: f64,
) -> Result<(BRep, SweepHistory), SweepError> {
    rotational_sweep_with_options(profile_pts, axis_origin, axis_dir, angle, SweepOptions::default())
}

/// Revolve a profile polygon around an axis with custom options.
pub fn rotational_sweep_with_options(
    profile_pts: &[DVec3],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle: f64,
    _options: SweepOptions,
) -> Result<(BRep, SweepHistory), SweepError> {
    if profile_pts.len() < 2 {
        return Err(SweepError::InsufficientVertices {
            minimum: 2,
            actual: profile_pts.len(),
        });
    }

    let axis_dir = normalize_vector("axis_dir", axis_dir)?;
    let angle = validate_positive("angle", angle)?;

    if !axis_origin.is_finite() {
        return Err(SweepError::NonFiniteInput("axis_origin"));
    }

    let n = profile_pts.len();
    let full_revolution = (angle - std::f64::consts::TAU).abs() < 1e-6;

    // Compute rotated positions
    let rot_pts: Vec<DVec3> = profile_pts.iter()
        .map(|&p| rotate_point(p, axis_origin, axis_dir, angle))
        .collect();

    let mut brep = BRep {
        vertices: Vec::new(),
        edges: Vec::new(),
        solids: Vec::new(),
        geom: rcad_kernel::GeomStore::default(),
        compound: None,
        compsolid: None,
    };

    // Add start vertices
    let start_vi: Vec<usize> = profile_pts.iter()
        .map(|&p| add_vertex(&mut brep, p))
        .collect();

    // Add end vertices (same as start for full revolution)
    let end_vi: Vec<usize> = if full_revolution {
        start_vi.clone()
    } else {
        rot_pts.iter()
            .map(|&p| add_vertex(&mut brep, p))
            .collect()
    };

    let mut faces = Vec::new();
    let mut history = SweepHistory::default();

    // Profile normal
    let profile_normal = polygon_normal(profile_pts);

    // Create start and end caps for partial revolution
    if !full_revolution {
        // Start cap
        {
            let mut wire_edges = Vec::with_capacity(n);
            for i in 0..n {
                let j = (i + 1) % n;
                let ei = add_line_edge(&mut brep, start_vi[i], start_vi[j]);
                wire_edges.push(WireEdge::fwd(ei));
            }
            let face = Face {
                outer_wire: Wire { edges: wire_edges },
                inner_wires: vec![],
                normal: profile_normal,
                triangles: vec![],
                mesh_dirty: true,
            };
            let fi = faces.len();
            faces.push(face);
            history.start_cap.push(fi);

            let si = brep.geom.surfaces.len();
            brep.geom.surfaces.push(Surface3::Plane(Plane {
                origin: profile_pts[0],
                normal: profile_normal,
            }));
            brep.geom.face_surface.push(Some(si));
        }

        // End cap
        {
            let mut wire_edges = Vec::with_capacity(n);
            for i in 0..n {
                let j = (i + 1) % n;
                let ei = add_line_edge(&mut brep, end_vi[i], end_vi[j]);
                wire_edges.push(WireEdge::rev(ei));
            }
            let rot_normal = rotate_point(profile_normal, axis_origin, axis_dir, angle)
                - axis_origin;
            let rot_normal = rot_normal.normalize_or(profile_normal);
            let face = Face {
                outer_wire: Wire { edges: wire_edges },
                inner_wires: vec![],
                normal: rot_normal,
                triangles: vec![],
                mesh_dirty: true,
            };
            let fi = faces.len();
            faces.push(face);
            history.end_cap.push(fi);

            let si = brep.geom.surfaces.len();
            brep.geom.surfaces.push(Surface3::Plane(Plane {
                origin: rot_pts[0],
                normal: rot_normal,
            }));
            brep.geom.face_surface.push(Some(si));
        }
    }

    // Create lateral faces
    for i in 0..n {
        let j = (i + 1) % n;

        let p0 = profile_pts[i];
        let p1 = profile_pts[j];
        let p1_rot = rot_pts[j];

        let lat_normal = (p1 - p0).cross(p1_rot - p0).normalize_or(profile_normal);

        // Create edges
        let e_bot = add_line_edge(&mut brep, start_vi[i], start_vi[j]);
        let e_right = add_line_edge(&mut brep, start_vi[j], end_vi[j]);
        let e_top = add_line_edge(&mut brep, end_vi[j], end_vi[i]);
        let e_left = add_line_edge(&mut brep, end_vi[i], start_vi[i]);

        let wire_edges = vec![
            WireEdge::fwd(e_bot),
            WireEdge::fwd(e_right),
            WireEdge::fwd(e_top),
            WireEdge::fwd(e_left),
        ];

        let face = Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: vec![],
            normal: lat_normal,
            triangles: vec![],
            mesh_dirty: true,
        };
        let fi = faces.len();
        faces.push(face);
        history.lateral_faces.push(fi);
        history.profile_edge_to_lateral.insert(i, fi);

        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: p0,
            normal: lat_normal,
        }));
        brep.geom.face_surface.push(Some(si));
    }

    brep.solids.push(Solid {
        shells: vec![Shell { faces }],
    });

    Ok((brep, history))
}

/// Revolve a BRep face around an axis.
pub fn rotational_sweep_face(
    brep: &BRep,
    face_idx: usize,
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle: f64,
) -> Result<BRep, SweepError> {
    let face = brep.solids.first()
        .and_then(|s| s.shells.first())
        .and_then(|sh| sh.faces.get(face_idx))
        .ok_or_else(|| SweepError::InvalidParameter("face_idx out of range"))?;

    let pts: Vec<DVec3> = face.outer_wire.edges.iter()
        .filter_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            let vi = if we.forward { edge.start } else { edge.end };
            brep.vertices.get(vi).map(|v| v.point)
        })
        .collect();

    rotational_sweep(&pts, axis_origin, axis_dir, angle)
}

/// Revolve a BRep wire around an axis.
pub fn rotational_sweep_wire(
    brep: &BRep,
    wire_idx: usize,
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle: f64,
) -> Result<BRep, SweepError> {
    rotational_sweep_face(brep, wire_idx, axis_origin, axis_dir, angle)
}

// ─────────────────────────────────────────────────────────────────────────────
// Pipe Sweep
// ─────────────────────────────────────────────────────────────────────────────

/// Sweep a 2D profile along a 3D spine path.
///
/// Creates a solid by sweeping a profile along a path. The profile orientation
/// is computed using Frenet frames (tangent, normal, binormal) at each spine point.
///
/// # Arguments
/// * `profile_2d` - Profile polygon in XY plane (local coordinate system)
/// * `spine` - Path points in 3D (must have at least 2 points)
/// * `mode` - Sweep mode (currently only Frenet supported)
///
/// # Returns
/// A closed BRep solid.
pub fn pipe_sweep(
    profile_2d: &[DVec2],
    spine: &[DVec3],
    _mode: SweepMode,
) -> Result<BRep, SweepError> {
    let (brep, _) = pipe_sweep_with_history(profile_2d, spine)?;
    Ok(brep)
}

/// Sweep a 2D profile along a 3D spine with history tracking.
pub fn pipe_sweep_with_history(
    profile_2d: &[DVec2],
    spine: &[DVec3],
) -> Result<(BRep, SweepHistory), SweepError> {
    pipe_sweep_with_options(profile_2d, spine, SweepOptions::pipe_frenet())
}

/// Sweep a 2D profile along a 3D spine with custom options.
pub fn pipe_sweep_with_options(
    profile_2d: &[DVec2],
    spine: &[DVec3],
    options: SweepOptions,
) -> Result<(BRep, SweepHistory), SweepError> {
    if profile_2d.len() < 3 {
        return Err(SweepError::InsufficientVertices {
            minimum: 3,
            actual: profile_2d.len(),
        });
    }
    if spine.len() < 2 {
        return Err(SweepError::InsufficientSpinePoints {
            minimum: 2,
            actual: spine.len(),
        });
    }

    // Compute tangents at each spine station
    let tangents = compute_spine_tangents(spine)?;

    // Compute Frenet frames at each station
    let frames = compute_frenet_frames(&tangents, options.is_frenet)?;

    // Transform 2D profile to 3D at each station
    let cross_sections: Vec<Vec<DVec3>> = spine.iter()
        .enumerate()
        .map(|(i, &origin)| {
            let frame = &frames[i];
            profile_2d.iter()
                .map(|p2| {
                    origin + p2.x * frame.right + p2.y * frame.up
                })
                .collect()
        })
        .collect();

    // Build the BRep by lofting between sections
    build_lofted_solid(&cross_sections, options.closed)
}

/// Compute tangent vectors at each spine point.
fn compute_spine_tangents(spine: &[DVec3]) -> Result<Vec<DVec3>, SweepError> {
    let n = spine.len();
    let mut tangents = Vec::with_capacity(n);

    for i in 0..n {
        let tan = if i == 0 {
            (spine[1] - spine[0]).normalize_or_zero()
        } else if i == n - 1 {
            (spine[n - 1] - spine[n - 2]).normalize_or_zero()
        } else {
            (spine[i + 1] - spine[i - 1]).normalize_or_zero()
        };

        if tan.length_squared() < EPS {
            return Err(SweepError::DegenerateGeometry("spine has zero-length segment"));
        }
        tangents.push(tan);
    }

    Ok(tangents)
}

/// Frenet frame (tangent, right, up) at a point on the spine.
#[derive(Debug, Clone, Copy)]
struct FrenetFrame {
    tangent: DVec3,
    right: DVec3,
    up: DVec3,
}

/// Compute Frenet frames along the spine.
fn compute_frenet_frames(tangents: &[DVec3], is_frenet: bool) -> Result<Vec<FrenetFrame>, SweepError> {
    let n = tangents.len();
    let mut frames = Vec::with_capacity(n);

    let world_up_primary = DVec3::Y;
    let world_up_fallback = DVec3::Z;

    let mut prev_up = DVec3::ZERO;

    for (i, &tan) in tangents.iter().enumerate() {
        let (right, up) = if is_frenet && i > 0 && prev_up.length_squared() > EPS {
            // Use previous up for smooth frame transition
            let right = tan.cross(prev_up).normalize_or_zero();
            if right.length_squared() < EPS {
                // Fallback for nearly parallel case
                let right_raw = tan.cross(world_up_primary);
                if right_raw.length_squared() > EPS {
                    (right_raw.normalize(), right_raw.normalize().cross(tan).normalize_or_zero())
                } else {
                    let r = tan.cross(world_up_fallback).normalize_or_zero();
                    (r, r.cross(tan).normalize_or_zero())
                }
            } else {
                (right, right.cross(tan).normalize_or_zero())
            }
        } else {
            // Initial frame or fixed mode
            let right_raw = tan.cross(world_up_primary);
            if right_raw.length_squared() > EPS {
                let r = right_raw.normalize();
                (r, r.cross(tan).normalize_or_zero())
            } else {
                let r = tan.cross(world_up_fallback).normalize_or_zero();
                (r, r.cross(tan).normalize_or_zero())
            }
        };

        prev_up = up;
        frames.push(FrenetFrame { tangent: tan, right, up });
    }

    Ok(frames)
}

/// Build a solid by lofting between cross-sections.
fn build_lofted_solid(sections: &[Vec<DVec3>], _closed: bool) -> Result<(BRep, SweepHistory), SweepError> {
    if sections.len() < 2 {
        return Err(SweepError::InsufficientSpinePoints {
            minimum: 2,
            actual: sections.len(),
        });
    }

    let n_verts = sections[0].len();
    let n_sections = sections.len();

    // Validate all sections have same vertex count
    for s in sections.iter() {
        if s.len() != n_verts {
            return Err(SweepError::VertexCountMismatch {
                expected: n_verts,
                actual: s.len(),
            });
        }
    }

    let mut brep = BRep {
        vertices: Vec::new(),
        edges: Vec::new(),
        solids: Vec::new(),
        geom: rcad_kernel::GeomStore::default(),
        compound: None,
        compsolid: None,
    };

    // Add all vertices
    let vi: Vec<Vec<usize>> = sections.iter()
        .map(|sec| {
            sec.iter().map(|&p| add_vertex(&mut brep, p)).collect()
        })
        .collect();

    let mut faces = Vec::new();
    let mut history = SweepHistory::default();

    // Start cap
    {
        let pts = &sections[0];
        let normal = -polygon_normal(pts);
        let mut wire_edges = Vec::with_capacity(n_verts);
        for i in 0..n_verts {
            let j = (i + 1) % n_verts;
            let ei = add_line_edge(&mut brep, vi[0][i], vi[0][j]);
            wire_edges.push(WireEdge::fwd(ei));
        }
        let face = Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: vec![],
            normal,
            triangles: vec![],
            mesh_dirty: true,
        };
        let fi = faces.len();
        faces.push(face);
        history.start_cap.push(fi);

        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: pts[0],
            normal,
        }));
        brep.geom.face_surface.push(Some(si));
    }

    // End cap
    {
        let pts = &sections[n_sections - 1];
        let normal = polygon_normal(pts);
        let mut wire_edges = Vec::with_capacity(n_verts);
        for i in 0..n_verts {
            let j = (i + 1) % n_verts;
            let ei = add_line_edge(&mut brep, vi[n_sections - 1][i], vi[n_sections - 1][j]);
            wire_edges.push(WireEdge::rev(ei));
        }
        let face = Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: vec![],
            normal,
            triangles: vec![],
            mesh_dirty: true,
        };
        let fi = faces.len();
        faces.push(face);
        history.end_cap.push(fi);

        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: pts[0],
            normal,
        }));
        brep.geom.face_surface.push(Some(si));
    }

    // Lateral faces between consecutive sections
    for sec in 0..n_sections - 1 {
        let pts_bot = &sections[sec];
        let pts_top = &sections[sec + 1];

        for i in 0..n_verts {
            let j = (i + 1) % n_verts;

            let a = pts_bot[i];
            let b = pts_bot[j];
            let c = pts_top[j];

            let lat_normal = (b - a).cross(c - a).normalize_or(DVec3::Z);

            let e_bot = add_line_edge(&mut brep, vi[sec][i], vi[sec][j]);
            let e_right = add_line_edge(&mut brep, vi[sec][j], vi[sec + 1][j]);
            let e_top = add_line_edge(&mut brep, vi[sec + 1][j], vi[sec + 1][i]);
            let e_left = add_line_edge(&mut brep, vi[sec + 1][i], vi[sec][i]);

            let wire_edges = vec![
                WireEdge::fwd(e_bot),
                WireEdge::fwd(e_right),
                WireEdge::fwd(e_top),
                WireEdge::fwd(e_left),
            ];

            let face = Face {
                outer_wire: Wire { edges: wire_edges },
                inner_wires: vec![],
                normal: lat_normal,
                triangles: vec![],
                mesh_dirty: true,
            };
            let fi = faces.len();
            faces.push(face);
            history.lateral_faces.push(fi);

            let si = brep.geom.surfaces.len();
            brep.geom.surfaces.push(Surface3::Plane(Plane {
                origin: a,
                normal: lat_normal,
            }));
            brep.geom.face_surface.push(Some(si));
        }
    }

    brep.solids.push(Solid {
        shells: vec![Shell { faces }],
    });

    Ok((brep, history))
}

/// Sweep a wire along a spine path.
pub fn pipe_sweep_wire(
    profile_2d: &[DVec2],
    spine: &[DVec3],
) -> Result<BRep, SweepError> {
    pipe_sweep(profile_2d, spine, SweepMode::Pipe)
}

/// Sweep with Frenet frame rotation (smooth rotation along spine).
pub fn pipe_with_rotation(
    profile_2d: &[DVec2],
    spine: &[DVec3],
) -> Result<BRep, SweepError> {
    let options = SweepOptions {
        mode: SweepMode::Pipe,
        is_frenet: true,
        continuous_normal: true,
        ..Default::default()
    };
    let (brep, _) = pipe_sweep_with_options(profile_2d, spine, options)?;
    Ok(brep)
}

// ─────────────────────────────────────────────────────────────────────────────
// Corner Handling
// ─────────────────────────────────────────────────────────────────────────────

/// Handle corners in a pipe sweep.
///
/// When the spine has sharp corners, this function creates smooth transitions
/// by inserting fillet-like geometry.
///
/// # Arguments
/// * `spine` - The path curve points
/// * `profile_2d` - The profile to sweep
/// * `corner_radius` - Radius for corner rounding
pub fn handle_pipe_corners(
    spine: &[DVec3],
    profile_2d: &[DVec2],
    corner_radius: f64,
) -> Result<BRep, SweepError> {
    if spine.len() < 3 {
        return Err(SweepError::InsufficientSpinePoints {
            minimum: 3,
            actual: spine.len(),
        });
    }
    if profile_2d.len() < 3 {
        return Err(SweepError::InsufficientVertices {
            minimum: 3,
            actual: profile_2d.len(),
        });
    }

    let corner_radius = validate_positive("corner_radius", corner_radius)?;

    // Find corners (points where direction changes significantly)
    let corners = find_corners(spine, 5.0_f64.to_radians());

    if corners.is_empty() {
        // No corners, do regular sweep
        return pipe_sweep(profile_2d, spine, SweepMode::Pipe);
    }

    // Build a modified spine with rounded corners
    let rounded_spine = round_spine_corners(spine, &corners, corner_radius)?;

    pipe_sweep(profile_2d, &rounded_spine, SweepMode::Pipe)
}

/// Find corner points where the spine direction changes significantly.
fn find_corners(spine: &[DVec3], angle_threshold: f64) -> Vec<usize> {
    let n = spine.len();
    if n < 3 {
        return vec![];
    }

    let mut corners = Vec::new();

    for i in 1..n - 1 {
        let dir_in = (spine[i] - spine[i - 1]).normalize_or_zero();
        let dir_out = (spine[i + 1] - spine[i]).normalize_or_zero();

        let dot = dir_in.dot(dir_out);
        let angle = dot.acos();

        if angle > angle_threshold {
            corners.push(i);
        }
    }

    corners
}

/// Round corners by inserting additional points.
fn round_spine_corners(
    spine: &[DVec3],
    corners: &[usize],
    radius: f64,
) -> Result<Vec<DVec3>, SweepError> {
    let mut result = Vec::new();
    let n = spine.len();

    let mut last_added = 0;

    for &corner_idx in corners {
        // Add points before the corner
        result.extend_from_slice(&spine[last_added..corner_idx]);

        // Compute the corner geometry
        let p_prev = spine[corner_idx - 1];
        let p_corner = spine[corner_idx];
        let p_next = spine[corner_idx + 1];

        let dir_in = (p_corner - p_prev).normalize_or_zero();
        let dir_out = (p_next - p_corner).normalize_or_zero();

        // Calculate arc points
        let arc_points = compute_corner_arc(p_prev, p_corner, p_next, dir_in, dir_out, radius);
        result.extend(arc_points);

        last_added = corner_idx + 1;
    }

    // Add remaining points
    if last_added < n {
        result.extend_from_slice(&spine[last_added..]);
    }

    Ok(result)
}

/// Compute arc points for a rounded corner.
#[allow(clippy::too_many_arguments)]
fn compute_corner_arc(
    _p_prev: DVec3,
    p_corner: DVec3,
    _p_next: DVec3,
    dir_in: DVec3,
    dir_out: DVec3,
    radius: f64,
) -> Vec<DVec3> {
    // Calculate the angle at the corner
    let dot = dir_in.dot(dir_out);
    let angle = dot.acos();

    if angle < 1e-6 || angle > std::f64::consts::PI - 1e-6 {
        return vec![p_corner];
    }

    // Calculate arc center and angles
    let half_angle = angle / 2.0;
    let dist_to_arc_start = radius / half_angle.tan();

    let p_start = p_corner - dir_in * dist_to_arc_start;
    let p_end = p_corner + dir_out * dist_to_arc_start;

    // Compute the bisector direction (center is along this)
    let bisector = (dir_in + dir_out).normalize_or_zero();
    let center_dist = radius / (angle / 2.0).sin();
    let center = p_corner + bisector * center_dist;

    // Generate arc points
    let n_arc_points = ((angle * radius) / 0.1).ceil() as usize;
    let n_arc_points = n_arc_points.max(3);

    let mut arc_points = Vec::with_capacity(n_arc_points);

    // Normal to the plane of the arc
    let normal = dir_in.cross(dir_out).normalize_or_zero();

    for i in 0..=n_arc_points {
        let t = i as f64 / n_arc_points as f64;
        let pt = p_start.lerp(p_end, t);

        // Project onto arc
        let to_pt = pt - center;
        let radial = to_pt - normal * to_pt.dot(normal);
        if radial.length_squared() > EPS {
            arc_points.push(center + radial.normalize() * radius);
        } else {
            arc_points.push(pt);
        }
    }

    arc_points
}

// ─────────────────────────────────────────────────────────────────────────────
// Law-Based Variable Section Sweeps
// ─────────────────────────────────────────────────────────────────────────────

/// A law function that defines how a parameter varies along the spine.
pub trait Law: std::fmt::Debug {
    /// Evaluate the law at parameter t (0 to 1).
    fn evaluate(&self, t: f64) -> f64;
}

/// A linear law: scales linearly from start_value to end_value.
#[derive(Debug, Clone)]
pub struct LinearLaw {
    pub start_value: f64,
    pub end_value: f64,
}

impl Law for LinearLaw {
    fn evaluate(&self, t: f64) -> f64 {
        self.start_value + t * (self.end_value - self.start_value)
    }
}

/// A constant law: returns a fixed value.
#[derive(Debug, Clone)]
pub struct ConstantLaw {
    pub value: f64,
}

impl Law for ConstantLaw {
    fn evaluate(&self, _t: f64) -> f64 {
        self.value
    }
}

/// A sinusoidal law.
#[derive(Debug, Clone)]
pub struct SineLaw {
    pub amplitude: f64,
    pub frequency: f64,
    pub phase: f64,
    pub offset: f64,
}

impl Law for SineLaw {
    fn evaluate(&self, t: f64) -> f64 {
        self.offset + self.amplitude * (self.frequency * t * std::f64::consts::TAU + self.phase).sin()
    }
}

/// A piecewise linear law defined by key points.
#[derive(Debug, Clone)]
pub struct PiecewiseLinearLaw {
    pub points: Vec<(f64, f64)>, // (t, value) pairs
}

impl Law for PiecewiseLinearLaw {
    fn evaluate(&self, t: f64) -> f64 {
        if self.points.is_empty() {
            return 0.0;
        }
        if self.points.len() == 1 {
            return self.points[0].1;
        }

        // Find the segment
        let mut i = 0;
        while i < self.points.len() - 1 && self.points[i + 1].0 < t {
            i += 1;
        }

        if i >= self.points.len() - 1 {
            return self.points.last().map(|p| p.1).unwrap_or(0.0);
        }

        let (t0, v0) = self.points[i];
        let (t1, v1) = self.points[i + 1];

        if (t1 - t0).abs() < EPS {
            return v0;
        }

        let local_t = (t - t0) / (t1 - t0);
        v0 + local_t * (v1 - v0)
    }
}

/// Perform a linear law sweep where the profile scales along the spine.
///
/// # Arguments
/// * `profile_2d` - Base profile in XY plane
/// * `spine` - Path points
/// * `law` - Law defining scale factor along the spine
pub fn linear_law_sweep<L: Law>(
    profile_2d: &[DVec2],
    spine: &[DVec3],
    law: L,
) -> Result<BRep, SweepError> {
    if profile_2d.len() < 3 {
        return Err(SweepError::InsufficientVertices {
            minimum: 3,
            actual: profile_2d.len(),
        });
    }
    if spine.len() < 2 {
        return Err(SweepError::InsufficientSpinePoints {
            minimum: 2,
            actual: spine.len(),
        });
    }

    let n_spine = spine.len();

    // Compute Frenet frames
    let tangents = compute_spine_tangents(spine)?;
    let frames = compute_frenet_frames(&tangents, true)?;

    // Transform and scale profile at each station
    let cross_sections: Vec<Vec<DVec3>> = spine.iter()
        .enumerate()
        .map(|(i, &origin)| {
            let frame = &frames[i];
            let t = i as f64 / (n_spine - 1) as f64;
            let scale = law.evaluate(t);

            profile_2d.iter()
                .map(|p2| {
                    origin + scale * p2.x * frame.right + scale * p2.y * frame.up
                })
                .collect()
        })
        .collect();

    build_lofted_solid(&cross_sections, false).map(|(brep, _)| brep)
}

/// Variable section sweep with different profiles at each station.
///
/// # Arguments
/// * `profiles` - Profiles at each station (in XY plane)
/// * `spine` - Path points (must have same length as profiles)
pub fn variable_section_sweep(
    profiles: &[Vec<DVec2>],
    spine: &[DVec3],
) -> Result<BRep, SweepError> {
    if profiles.len() != spine.len() {
        return Err(SweepError::InvalidParameter(
            "profiles and spine must have the same length"
        ));
    }
    if profiles.len() < 2 {
        return Err(SweepError::InsufficientSpinePoints {
            minimum: 2,
            actual: profiles.len(),
        });
    }

    let n_verts = profiles.first().map(|p| p.len()).unwrap_or(0);
    if n_verts < 3 {
        return Err(SweepError::InsufficientVertices {
            minimum: 3,
            actual: n_verts,
        });
    }

    // Validate all profiles have same vertex count
    for p in profiles.iter() {
        if p.len() != n_verts {
            return Err(SweepError::VertexCountMismatch {
                expected: n_verts,
                actual: p.len(),
            });
        }
    }

    // Compute Frenet frames
    let tangents = compute_spine_tangents(spine)?;
    let frames = compute_frenet_frames(&tangents, true)?;

    // Transform each profile to 3D
    let cross_sections: Vec<Vec<DVec3>> = spine.iter()
        .enumerate()
        .map(|(i, &origin)| {
            let frame = &frames[i];
            profiles[i].iter()
                .map(|p2| {
                    origin + p2.x * frame.right + p2.y * frame.up
                })
                .collect()
        })
        .collect();

    build_lofted_solid(&cross_sections, false).map(|(brep, _)| brep)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn square_profile() -> Vec<DVec3> {
        vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ]
    }

    fn triangle_profile_2d() -> Vec<DVec2> {
        vec![
            DVec2::new(-0.5, 0.0),
            DVec2::new(0.5, 0.0),
            DVec2::new(0.0, 1.0),
        ]
    }

    fn square_profile_2d() -> Vec<DVec2> {
        vec![
            DVec2::new(-0.5, -0.5),
            DVec2::new(0.5, -0.5),
            DVec2::new(0.5, 0.5),
            DVec2::new(-0.5, 0.5),
        ]
    }

    fn linear_spine() -> Vec<DVec3> {
        vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(0.0, 0.0, 2.0),
        ]
    }

    fn arc_spine() -> Vec<DVec3> {
        let n = 10;
        (0..n).map(|i| {
            let t = i as f64 / (n - 1) as f64;
            let angle = t * std::f64::consts::FRAC_PI_2;
            DVec3::new(angle.cos(), angle.sin(), 0.0)
        }).collect()
    }

    // Linear sweep tests
    #[test]
    fn linear_sweep_creates_box() {
        let profile = square_profile();
        let result = linear_sweep(&profile, DVec3::Z, 2.0).unwrap();

        assert!(!result.solids.is_empty(), "should have solids");
        let shell = &result.solids[0].shells[0];
        assert_eq!(shell.faces.len(), 6, "box should have 6 faces");
    }

    #[test]
    fn linear_sweep_rejects_zero_direction() {
        let profile = square_profile();
        let err = linear_sweep(&profile, DVec3::ZERO, 1.0).unwrap_err();
        assert!(matches!(err, SweepError::ZeroVector(_)));
    }

    #[test]
    fn linear_sweep_rejects_negative_distance() {
        let profile = square_profile();
        let err = linear_sweep(&profile, DVec3::Z, -1.0).unwrap_err();
        assert!(matches!(err, SweepError::NonPositiveInput(_)));
    }

    #[test]
    fn linear_sweep_rejects_insufficient_vertices() {
        let profile = vec![DVec3::ZERO, DVec3::X];
        let err = linear_sweep(&profile, DVec3::Z, 1.0).unwrap_err();
        assert!(matches!(err, SweepError::InsufficientVertices { .. }));
    }

    // Rotational sweep tests
    #[test]
    fn rotational_sweep_creates_solid() {
        // Create a simple rectangle profile
        let profile = vec![
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 1.0),
            DVec3::new(1.0, 0.0, 1.0),
        ];

        let result = rotational_sweep(&profile, DVec3::ZERO, DVec3::Z,
                                      std::f64::consts::FRAC_PI_2).unwrap();

        assert!(!result.solids.is_empty(), "should have solids");
        assert!(!result.solids[0].shells[0].faces.is_empty(), "should have faces");
    }

    #[test]
    fn rotational_sweep_full_revolution() {
        let profile = vec![
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 1.0),
            DVec3::new(1.0, 0.0, 1.0),
        ];

        let result = rotational_sweep(&profile, DVec3::ZERO, DVec3::Z,
                                      std::f64::consts::TAU).unwrap();

        assert!(!result.solids.is_empty());
    }

    #[test]
    fn rotational_sweep_rejects_zero_axis() {
        let profile = vec![DVec3::X, DVec3::new(2.0, 0.0, 0.0)];
        let err = rotational_sweep(&profile, DVec3::ZERO, DVec3::ZERO, 1.0).unwrap_err();
        assert!(matches!(err, SweepError::ZeroVector(_)));
    }

    // Pipe sweep tests
    #[test]
    fn pipe_sweep_linear_spine() {
        let profile = triangle_profile_2d();
        let spine = linear_spine();

        let result = pipe_sweep(&profile, &spine, SweepMode::Pipe).unwrap();

        assert!(!result.solids.is_empty(), "should have solids");
        // Should have start cap + end cap + lateral faces
        assert!(result.solids[0].shells[0].faces.len() >= 5);
    }

    #[test]
    fn pipe_sweep_arc_spine() {
        let profile = square_profile_2d();
        let spine = arc_spine();

        let result = pipe_sweep(&profile, &spine, SweepMode::Pipe).unwrap();

        assert!(!result.solids.is_empty());
    }

    #[test]
    fn pipe_sweep_rejects_short_spine() {
        let profile = triangle_profile_2d();
        let spine = vec![DVec3::ZERO];

        let err = pipe_sweep(&profile, &spine, SweepMode::Pipe).unwrap_err();
        assert!(matches!(err, SweepError::InsufficientSpinePoints { .. }));
    }

    #[test]
    fn pipe_with_rotation_creates_solid() {
        let profile = square_profile_2d();
        let spine = arc_spine();

        let result = pipe_with_rotation(&profile, &spine).unwrap();
        assert!(!result.solids.is_empty());
    }

    // Corner handling tests
    #[test]
    fn handle_corners_no_corners() {
        let profile = square_profile_2d();
        let spine = linear_spine();

        let result = handle_pipe_corners(&spine, &profile, 0.1).unwrap();
        assert!(!result.solids.is_empty());
    }

    #[test]
    fn handle_corners_with_sharp_corner() {
        let profile = square_profile_2d();
        let spine = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(1.0, 2.0, 0.0),
        ];

        let result = handle_pipe_corners(&spine, &profile, 0.1).unwrap();
        assert!(!result.solids.is_empty());
    }

    // Law tests
    #[test]
    fn linear_law_evaluates_correctly() {
        let law = LinearLaw { start_value: 1.0, end_value: 2.0 };
        assert!((law.evaluate(0.0) - 1.0).abs() < 1e-10);
        assert!((law.evaluate(0.5) - 1.5).abs() < 1e-10);
        assert!((law.evaluate(1.0) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn constant_law_evaluates_correctly() {
        let law = ConstantLaw { value: 5.0 };
        assert!((law.evaluate(0.0) - 5.0).abs() < 1e-10);
        assert!((law.evaluate(0.5) - 5.0).abs() < 1e-10);
        assert!((law.evaluate(1.0) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn sine_law_evaluates_correctly() {
        let law = SineLaw {
            amplitude: 1.0,
            frequency: 1.0,
            phase: 0.0,
            offset: 0.0,
        };
        assert!((law.evaluate(0.0) - 0.0).abs() < 1e-10);
        assert!((law.evaluate(0.25) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn linear_law_sweep_creates_tapered_shape() {
        let profile = square_profile_2d();
        let spine = linear_spine();
        let law = LinearLaw { start_value: 1.0, end_value: 0.5 };

        let result = linear_law_sweep(&profile, &spine, law).unwrap();
        assert!(!result.solids.is_empty());
    }

    #[test]
    fn variable_section_sweep_creates_solid() {
        let profile1: Vec<DVec2> = vec![
            DVec2::new(-1.0, -1.0),
            DVec2::new(1.0, -1.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(-1.0, 1.0),
        ];
        let profile2: Vec<DVec2> = vec![
            DVec2::new(-0.5, -0.5),
            DVec2::new(0.5, -0.5),
            DVec2::new(0.5, 0.5),
            DVec2::new(-0.5, 0.5),
        ];

        let profiles = vec![profile1, profile2];
        let spine = vec![DVec3::ZERO, DVec3::Z];

        let result = variable_section_sweep(&profiles, &spine).unwrap();
        assert!(!result.solids.is_empty());
    }

    // Options tests
    #[test]
    fn sweep_options_defaults() {
        let opts = SweepOptions::default();
        assert_eq!(opts.mode, SweepMode::Translation);
        assert!(opts.is_frenet);
        assert_eq!(opts.corner_mode, CornerMode::Sharp);
        assert!(!opts.closed);
    }

    #[test]
    fn sweep_options_builders() {
        let opts = SweepOptions::linear();
        assert_eq!(opts.mode, SweepMode::Translation);

        let opts = SweepOptions::rotational();
        assert_eq!(opts.mode, SweepMode::Rotation);

        let opts = SweepOptions::pipe_frenet();
        assert_eq!(opts.mode, SweepMode::Pipe);
        assert!(opts.is_frenet);

        let opts = SweepOptions::pipe_fixed();
        assert_eq!(opts.mode, SweepMode::Pipe);
        assert!(!opts.is_frenet);
    }

    #[test]
    fn sweep_options_chain() {
        let opts = SweepOptions::pipe_frenet()
            .with_corner_mode(CornerMode::Rounded)
            .with_corner_radius(0.5)
            .with_continuous_normal(true)
            .with_closed(true);

        assert_eq!(opts.corner_mode, CornerMode::Rounded);
        assert!((opts.corner_radius - 0.5).abs() < 1e-10);
        assert!(opts.continuous_normal);
        assert!(opts.closed);
    }

    // History tests
    #[test]
    fn linear_sweep_history() {
        let profile = square_profile();
        let (_, history) = linear_sweep_with_history(&profile, DVec3::Z, 1.0).unwrap();

        assert_eq!(history.start_cap.len(), 1);
        assert_eq!(history.end_cap.len(), 1);
        assert_eq!(history.lateral_faces.len(), 4);
        assert_eq!(history.all_faces().count(), 6);
    }

    // Edge cases
    #[test]
    fn find_corners_detects_sharp_turn() {
        let spine = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(2.0, 1.0, 0.0),
        ];

        let corners = find_corners(&spine, 5.0_f64.to_radians());
        // Should detect at least one corner at the 90-degree turn
        assert!(!corners.is_empty());
    }

    #[test]
    fn find_corners_no_sharp_turn() {
        let spine = linear_spine();
        let corners = find_corners(&spine, 5.0_f64.to_radians());
        assert!(corners.is_empty());
    }
}
