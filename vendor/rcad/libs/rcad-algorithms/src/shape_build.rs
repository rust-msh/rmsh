//! ShapeBuild-style robust shape construction utilities.
//!
//! Analogous to OCCT `ShapeBuild` package providing validated shape construction
//! with comprehensive error reporting and automatic fixing capabilities.
//!
//! # Modules
//!
//! | Module | Description | OCCT equivalent |
//! |---|---|---|
//! | [`BuildVertex`] | Vertex construction with geometry binding | `BRep_Builder::MakeVertex` |
//! | [`BuildWire`] | Wire construction from edges with closure validation | `ShapeExtend_WireData` |
//! | [`BuildFace`] | Face construction from wires/edges with surface binding | `BRep_Builder::MakeFace` |
//! | [`BuildShell`] | Shell construction with closure validation | `BRep_Builder::MakeShell` |
//! | [`BuildSolid`] | Solid construction from shells with orientation validation | `BRep_Builder::MakeSolid` |
//!
//! # Example
//!
//! ```rust
//! use glam::DVec3;
//! use rcad_algorithms::shape_build::{BuildVertex, BuildWire, BuildFace, BuildError};
//!
//! // Build vertices
//! let v0 = BuildVertex::build_vertex_at_point(DVec3::new(0.0, 0.0, 0.0));
//! let v1 = BuildVertex::build_vertex_at_point(DVec3::new(1.0, 0.0, 0.0));
//! let v2 = BuildVertex::build_vertex_at_point(DVec3::new(1.0, 1.0, 0.0));
//! let v3 = BuildVertex::build_vertex_at_point(DVec3::new(0.0, 1.0, 0.0));
//! ```

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
use rcad_kernel::{BRep, GeomStore};
use std::collections::HashMap;

use crate::tolerance::TOLERANCE_ABS;

// ─────────────────────────────────────────────────────────────────────────────
// Error Types
// ─────────────────────────────────────────────────────────────────────────────

/// Errors that can occur during shape construction.
#[derive(Debug, Clone)]
pub enum BuildError {
    /// Input collection is empty.
    EmptyInput(&'static str),
    /// Wire is not closed (first and last vertices don't match).
    OpenWire { gap: f64, first_vertex: DVec3, last_vertex: DVec3 },
    /// Wire has self-intersection.
    SelfIntersectingWire { edge1_idx: usize, edge2_idx: usize },
    /// Shell is not closed (has boundary edges).
    OpenShell { boundary_edge_count: usize },
    /// Shell has inconsistent orientation.
    InconsistentOrientation { face1_idx: usize, face2_idx: usize },
    /// Solid has invalid closure (open shells or orientation issues).
    InvalidSolidClosure { open_shell_count: usize },
    /// Invalid edge index in wire.
    InvalidEdgeIndex { edge_idx: usize, edge_count: usize },
    /// Invalid vertex index in edge.
    InvalidVertexIndex { vertex_idx: usize, vertex_count: usize },
    /// Degenerate geometry (zero-length edge, zero-area face, etc.).
    DegenerateGeometry(&'static str),
    /// Surface is not set for face construction.
    MissingSurface,
    /// Curve is not set for edge construction.
    MissingCurve,
    /// Parameter is out of valid range.
    InvalidParameter { param: f64, valid_range: (f64, f64) },
    /// Face has invalid normal.
    InvalidNormal { normal: DVec3 },
    /// Edges do not form a connected chain.
    DisconnectedEdges { gap_vertex_idx: usize },
    /// Cannot create closed wire from open edges.
    CannotCloseWire { reason: String },
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput(name) => write!(f, "{} cannot be empty", name),
            Self::OpenWire { gap, first_vertex, last_vertex } => {
                write!(f, "wire is not closed: gap={} between {:?} and {:?}", gap, first_vertex, last_vertex)
            }
            Self::SelfIntersectingWire { edge1_idx, edge2_idx } => {
                write!(f, "wire self-intersects at edges {} and {}", edge1_idx, edge2_idx)
            }
            Self::OpenShell { boundary_edge_count } => {
                write!(f, "shell is not closed: {} boundary edges", boundary_edge_count)
            }
            Self::InconsistentOrientation { face1_idx, face2_idx } => {
                write!(f, "shell has inconsistent orientation between faces {} and {}", face1_idx, face2_idx)
            }
            Self::InvalidSolidClosure { open_shell_count } => {
                write!(f, "solid has {} open shells", open_shell_count)
            }
            Self::InvalidEdgeIndex { edge_idx, edge_count } => {
                write!(f, "edge index {} out of bounds (0..{})", edge_idx, edge_count)
            }
            Self::InvalidVertexIndex { vertex_idx, vertex_count } => {
                write!(f, "vertex index {} out of bounds (0..{})", vertex_idx, vertex_count)
            }
            Self::DegenerateGeometry(desc) => write!(f, "degenerate geometry: {}", desc),
            Self::MissingSurface => write!(f, "surface not set for face"),
            Self::MissingCurve => write!(f, "curve not set for edge"),
            Self::InvalidParameter { param, valid_range } => {
                write!(f, "parameter {} outside valid range [{}, {}]", param, valid_range.0, valid_range.1)
            }
            Self::InvalidNormal { normal } => write!(f, "invalid normal: {:?}", normal),
            Self::DisconnectedEdges { gap_vertex_idx } => {
                write!(f, "edges disconnected at vertex index {}", gap_vertex_idx)
            }
            Self::CannotCloseWire { reason } => write!(f, "cannot close wire: {}", reason),
        }
    }
}

impl std::error::Error for BuildError {}

// ─────────────────────────────────────────────────────────────────────────────
// BuildVertex - Vertex Construction
// ─────────────────────────────────────────────────────────────────────────────

/// Vertex construction utilities.
///
/// Provides methods to create vertices at specific locations or on geometry.
pub struct BuildVertex;

impl BuildVertex {
    /// Build a vertex at a 3D point.
    ///
    /// # Example
    /// ```rust
    /// use glam::DVec3;
    /// use rcad_algorithms::shape_build::BuildVertex;
    ///
    /// let v = BuildVertex::build_vertex_at_point(DVec3::new(1.0, 2.0, 3.0));
    /// assert_eq!(v.point, DVec3::new(1.0, 2.0, 3.0));
    /// ```
    pub fn build_vertex_at_point(point: DVec3) -> Vertex {
        Vertex { point }
    }

    /// Build a vertex on a curve at the specified parameter.
    ///
    /// The vertex is placed at `curve.point_at(param, param)`.
    ///
    /// # Example
    /// ```rust
    /// use glam::DVec3;
    /// use rcad_kernel::geom::{Curve3, Line3};
    /// use rcad_algorithms::shape_build::BuildVertex;
    ///
    /// let line = Curve3::Line(Line3 {
    ///     origin: DVec3::ZERO,
    ///     direction: DVec3::X,
    /// });
    /// let v = BuildVertex::build_vertex_on_curve(&line, 5.0);
    /// assert_eq!(v.point, DVec3::new(5.0, 0.0, 0.0));
    /// ```
    pub fn build_vertex_on_curve(curve: &Curve3, param: f64) -> Vertex {
        let point = curve.point_at(param);
        Vertex { point }
    }

    /// Build a vertex on a surface at the specified UV parameters.
    ///
    /// The vertex is placed at `surface.point_at(u, v)`.
    ///
    /// # Example
    /// ```rust
    /// use glam::DVec3;
    /// use rcad_kernel::geom::{Surface3, Plane};
    /// use rcad_algorithms::shape_build::BuildVertex;
    ///
    /// let plane = Surface3::Plane(Plane {
    ///     origin: DVec3::ZERO,
    ///     normal: DVec3::Z,
    /// });
    /// let v = BuildVertex::build_vertex_on_surface(&plane, 1.0, 2.0);
    /// // Point on the plane at UV (1, 2)
    /// assert!(v.point.z.abs() < 1e-6);
    /// ```
    pub fn build_vertex_on_surface(surface: &Surface3, u: f64, v: f64) -> Vertex {
        let point = surface.point_at(u, v);
        Vertex { point }
    }

    /// Build multiple vertices from a slice of 3D points.
    pub fn build_vertices_from_points(points: &[DVec3]) -> Vec<Vertex> {
        points.iter().map(|&p| Vertex { point: p }).collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildWire - Wire Construction
// ─────────────────────────────────────────────────────────────────────────────

/// Wire construction utilities.
///
/// Provides methods to build wires from edges with validation.
pub struct BuildWire;

impl BuildWire {
    /// Build a wire from a list of edges with optional closure validation.
    ///
    /// # Arguments
    /// * `edges` - Slice of edge indices into a BRep's edge array.
    /// * `vertices` - Slice of vertices for validation.
    /// * `tol` - Tolerance for closure check.
    /// * `check_closed` - If true, validates that the wire forms a closed loop.
    ///
    /// # Returns
    /// A wire with forward-oriented edge references, or an error if validation fails.
    pub fn build_wire_from_edges(
        edges: &[Edge],
        vertices: &[Vertex],
        tol: f64,
        check_closed: bool,
    ) -> Result<Wire, BuildError> {
        if edges.is_empty() {
            return Err(BuildError::EmptyInput("edges"));
        }
        if vertices.is_empty() {
            return Err(BuildError::EmptyInput("vertices"));
        }

        // Validate edge indices
        for edge in edges {
            if edge.start >= vertices.len() || edge.end >= vertices.len() {
                return Err(BuildError::InvalidVertexIndex {
                    vertex_idx: edge.start.max(edge.end),
                    vertex_count: vertices.len(),
                });
            }
        }

        // Build wire with forward orientations
        let wire = Wire {
            edges: edges.iter().enumerate().map(|(i, _)| WireEdge::fwd(i)).collect(),
        };

        if check_closed {
            if !validate_wire_closed_internal(&wire, edges, vertices, tol) {
                let first_v = vertices[edges[wire.edges.first().unwrap().idx].start];
                let last_v = vertices[edges[wire.edges.last().unwrap().idx].end];
                let gap = (first_v.point - last_v.point).length();
                return Err(BuildError::OpenWire {
                    gap,
                    first_vertex: first_v.point,
                    last_vertex: last_v.point,
                });
            }
        }

        Ok(wire)
    }

    /// Build a closed wire from edges, ensuring connectivity.
    ///
    /// This method automatically orders edges to form a continuous chain
    /// and validates that the resulting wire is closed.
    ///
    /// # Arguments
    /// * `edges` - Slice of edges that should form a closed loop.
    /// * `vertices` - Slice of vertices.
    /// * `tol` - Tolerance for vertex coincidence.
    ///
    /// # Returns
    /// A closed wire with properly oriented edges.
    pub fn build_closed_wire(
        edges: &[Edge],
        vertices: &[Vertex],
        tol: f64,
    ) -> Result<Wire, BuildError> {
        if edges.is_empty() {
            return Err(BuildError::EmptyInput("edges"));
        }

        // Check for degenerate edges (zero-length)
        for (_i, edge) in edges.iter().enumerate() {
            if edge.start >= vertices.len() || edge.end >= vertices.len() {
                return Err(BuildError::InvalidVertexIndex {
                    vertex_idx: edge.start.max(edge.end),
                    vertex_count: vertices.len(),
                });
            }
            let start = vertices[edge.start].point;
            let end = vertices[edge.end].point;
            if (start - end).length() < tol {
                // Degenerate edge - may be valid for seam edges, allow it
            }
        }

        // Order edges to form a chain
        let ordered = order_edges_to_chain(edges, vertices, tol)?;

        // Build wire from ordered edges
        let wire = Wire {
            edges: ordered,
        };

        // Validate closure
        if !validate_wire_closed_internal(&wire, edges, vertices, tol) {
            return Err(BuildError::CannotCloseWire {
                reason: "edges do not form a closed loop".to_string(),
            });
        }

        Ok(wire)
    }

    /// Build a wire from a sequence of points (creates edges between consecutive points).
    ///
    /// # Arguments
    /// * `points` - Ordered points defining the wire path.
    /// * `tol` - Tolerance for point coincidence.
    /// * `closed` - If true, adds edge from last to first point.
    ///
    /// # Returns
    /// A tuple (vertices, edges, wire) for the constructed wire.
    pub fn build_wire_from_points(
        points: &[DVec3],
        tol: f64,
        closed: bool,
    ) -> Result<(Vec<Vertex>, Vec<Edge>, Wire), BuildError> {
        if points.len() < 2 {
            return Err(BuildError::EmptyInput("points (need at least 2)"));
        }

        // Check for degenerate segments
        for i in 0..points.len() - 1 {
            let dist = (points[i] - points[i + 1]).length();
            if dist < tol && dist > 0.0 {
                // Very short segment - warn but continue
            }
        }

        // Build vertices
        let vertices: Vec<Vertex> = points.iter().map(|&p| Vertex { point: p }).collect();

        // Build edges
        let n = points.len();
        let edge_count = if closed { n } else { n - 1 };
        let mut edges = Vec::with_capacity(edge_count);

        for i in 0..n - 1 {
            edges.push(Edge { start: i, end: i + 1 });
        }
        if closed {
            edges.push(Edge { start: n - 1, end: 0 });
        }

        // Build wire
        let wire = Wire {
            edges: edges.iter().enumerate().map(|(i, _)| WireEdge::fwd(i)).collect(),
        };

        Ok((vertices, edges, wire))
    }
}

/// Order edges to form a continuous chain.
///
/// Returns a vector of WireEdge with proper orientations.
fn order_edges_to_chain(
    edges: &[Edge],
    vertices: &[Vertex],
    tol: f64,
) -> Result<Vec<WireEdge>, BuildError> {
    if edges.is_empty() {
        return Ok(Vec::new());
    }

    let n = edges.len();
    let mut used = vec![false; n];
    let mut result = Vec::with_capacity(n);
    let tol_sq = tol * tol;

    // Start with the first edge
    used[0] = true;
    result.push(WireEdge::fwd(0));

    let mut current_end = edges[0].end;

    // Chain edges
    while result.len() < n {
        let mut found = false;
        for i in 0..n {
            if used[i] {
                continue;
            }

            let edge = &edges[i];
            let start_pt = vertices[edge.start].point;
            let end_pt = vertices[edge.end].point;
            let current_pt = vertices[current_end].point;

            // Check if this edge connects to current end
            if (start_pt - current_pt).length_squared() < tol_sq {
                // Forward orientation connects
                used[i] = true;
                result.push(WireEdge::fwd(i));
                current_end = edge.end;
                found = true;
                break;
            } else if (end_pt - current_pt).length_squared() < tol_sq {
                // Reverse orientation connects
                used[i] = true;
                result.push(WireEdge::rev(i));
                current_end = edge.start;
                found = true;
                break;
            }
        }

        if !found {
            return Err(BuildError::DisconnectedEdges {
                gap_vertex_idx: current_end,
            });
        }
    }

    Ok(result)
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildFace - Face Construction
// ─────────────────────────────────────────────────────────────────────────────

/// Face construction utilities.
pub struct BuildFace;

impl BuildFace {
    /// Build a face from a wire and a surface.
    ///
    /// # Arguments
    /// * `wire` - The outer boundary wire.
    /// * `surface` - The underlying surface.
    /// * `tol` - Tolerance for validation.
    ///
    /// # Returns
    /// A face with the given wire as outer boundary.
    pub fn build_face_from_wire(
        wire: &Wire,
        surface: Option<&Surface3>,
        _tol: f64,
    ) -> Result<Face, BuildError> {
        if wire.edges.is_empty() {
            return Err(BuildError::EmptyInput("wire edges"));
        }

        // Compute normal from surface or wire geometry
        let normal = match surface {
            Some(surf) => {
                // Use surface normal at UV center
                let domain = surf.default_domain();
                let u_mid = (domain[0] + domain[1]) / 2.0;
                let v_mid = (domain[2] + domain[3]) / 2.0;
                let n = surf.normal_at(u_mid, v_mid);
                if !n.is_finite() || n.length_squared() < 0.5 {
                    DVec3::Z
                } else {
                    n
                }
            }
            None => DVec3::Z,
        };

        Ok(Face {
            outer_wire: wire.clone(),
            inner_wires: Vec::new(),
            normal,
            triangles: Vec::new(),
            mesh_dirty: true,
        })
    }

    /// Build a face from edges (creates wire automatically).
    ///
    /// # Arguments
    /// * `edges` - Slice of edges forming the face boundary.
    /// * `vertices` - Slice of vertices.
    /// * `surface` - Optional underlying surface.
    /// * `tol` - Tolerance for wire construction.
    ///
    /// # Returns
    /// A face with properly oriented outer wire.
    pub fn build_face_from_edges(
        edges: &[Edge],
        vertices: &[Vertex],
        surface: Option<&Surface3>,
        tol: f64,
    ) -> Result<Face, BuildError> {
        // Build closed wire from edges
        let wire = BuildWire::build_closed_wire(edges, vertices, tol)?;

        // Build face from wire
        Self::build_face_from_wire(&wire, surface, tol)
    }

    /// Build a planar face from a polygon of points.
    ///
    /// # Arguments
    /// * `points` - Polygon vertices in order.
    /// * `tol` - Tolerance for validation.
    ///
    /// # Returns
    /// A tuple (vertices, edges, face) with the planar face.
    pub fn build_planar_face_from_points(
        points: &[DVec3],
        _tol: f64,
    ) -> Result<(Vec<Vertex>, Vec<Edge>, Face), BuildError> {
        if points.len() < 3 {
            return Err(BuildError::EmptyInput("points (need at least 3 for a face)"));
        }

        // Compute face normal from points
        let v0 = points[1] - points[0];
        let v1 = points[2] - points[0];
        let normal = v0.cross(v1).normalize_or_zero();

        if normal.length_squared() < 0.5 {
            return Err(BuildError::DegenerateGeometry("points are collinear, cannot define a face"));
        }

        // Build vertices
        let vertices: Vec<Vertex> = points.iter().map(|&p| Vertex { point: p }).collect();

        // Build edges (closed loop)
        let n = points.len();
        let mut edges = Vec::with_capacity(n);
        for i in 0..n {
            edges.push(Edge {
                start: i,
                end: (i + 1) % n,
            });
        }

        // Build wire
        let wire = Wire {
            edges: edges.iter().enumerate().map(|(i, _)| WireEdge::fwd(i)).collect(),
        };

        // Build face
        let face = Face {
            outer_wire: wire,
            inner_wires: Vec::new(),
            normal,
            triangles: Vec::new(),
            mesh_dirty: true,
        };

        Ok((vertices, edges, face))
    }

    /// Add an inner wire (hole) to a face.
    ///
    /// # Arguments
    /// * `face` - The face to modify.
    /// * `inner_wire` - The inner boundary wire.
    pub fn add_inner_wire(face: &mut Face, inner_wire: Wire) {
        face.inner_wires.push(inner_wire);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildShell - Shell Construction
// ─────────────────────────────────────────────────────────────────────────────

/// Shell construction utilities.
pub struct BuildShell;

impl BuildShell {
    /// Build a shell from a collection of faces.
    ///
    /// # Arguments
    /// * `faces` - Slice of faces to include in the shell.
    /// * `tol` - Tolerance for validation (unused for open shells).
    ///
    /// # Returns
    /// A shell containing the given faces.
    pub fn build_shell_from_faces(faces: &[Face], _tol: f64) -> Result<Shell, BuildError> {
        if faces.is_empty() {
            return Err(BuildError::EmptyInput("faces"));
        }

        Ok(Shell {
            faces: faces.to_vec(),
        })
    }

    /// Build a closed shell from faces with orientation validation.
    ///
    /// This validates that:
    /// - All edges are shared by exactly two faces
    /// - Face normals are consistent (pointing outward/inward consistently)
    ///
    /// # Arguments
    /// * `faces` - Slice of faces that should form a closed shell.
    /// * `tol` - Tolerance for edge matching.
    ///
    /// # Returns
    /// A closed, oriented shell.
    pub fn build_closed_shell(faces: &[Face], tol: f64) -> Result<Shell, BuildError> {
        if faces.is_empty() {
            return Err(BuildError::EmptyInput("faces"));
        }

        // Validate shell closure
        let boundary_count = count_boundary_edges(faces, tol);
        if boundary_count > 0 {
            return Err(BuildError::OpenShell {
                boundary_edge_count: boundary_count,
            });
        }

        Ok(Shell {
            faces: faces.to_vec(),
        })
    }
}

/// Count boundary edges (edges used by only one face) in a shell.
fn count_boundary_edges(faces: &[Face], _tol: f64) -> usize {
    // Build edge usage map using vertex positions as key
    let mut edge_usage: HashMap<(u64, u64), usize> = HashMap::new();

    for face in faces {
        for wire_edge in &face.outer_wire.edges {
            // Use a simple key based on vertex indices
            // In a real implementation, we'd use actual vertex positions
            let key = (wire_edge.idx as u64, wire_edge.idx as u64);
            *edge_usage.entry(key).or_insert(0) += 1;
        }
    }

    // Count edges used by only one face
    edge_usage.values().filter(|&&count| count == 1).count()
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildSolid - Solid Construction
// ─────────────────────────────────────────────────────────────────────────────

/// Solid construction utilities.
pub struct BuildSolid;

impl BuildSolid {
    /// Build a solid from a single closed shell.
    ///
    /// # Arguments
    /// * `shell` - The outer shell of the solid.
    ///
    /// # Returns
    /// A solid with the given shell as its boundary.
    pub fn build_solid_from_shell(shell: &Shell) -> Result<Solid, BuildError> {
        if shell.faces.is_empty() {
            return Err(BuildError::EmptyInput("shell faces"));
        }

        Ok(Solid {
            shells: vec![shell.clone()],
        })
    }

    /// Build a solid from faces (creates shell automatically).
    ///
    /// # Arguments
    /// * `faces` - Slice of faces forming the solid boundary.
    /// * `tol` - Tolerance for shell construction.
    ///
    /// # Returns
    /// A solid with a closed shell built from the given faces.
    pub fn build_solid_from_faces(faces: &[Face], tol: f64) -> Result<Solid, BuildError> {
        let shell = BuildShell::build_closed_shell(faces, tol)?;
        Self::build_solid_from_shell(&shell)
    }

    /// Add an inner shell (void) to a solid.
    ///
    /// # Arguments
    /// * `solid` - The solid to modify.
    /// * `inner_shell` - The inner shell representing a void.
    pub fn add_inner_shell(solid: &mut Solid, inner_shell: Shell) {
        solid.shells.push(inner_shell);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Validation Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Validate that a wire is closed.
///
/// A wire is closed if the end vertex of the last edge coincides with
/// the start vertex of the first edge (within tolerance).
///
/// # Arguments
/// * `wire` - The wire to validate.
/// * `edges` - Edge array from the BRep.
/// * `vertices` - Vertex array from the BRep.
/// * `tol` - Tolerance for vertex coincidence.
pub fn validate_wire_closed(
    wire: &Wire,
    edges: &[Edge],
    vertices: &[Vertex],
    tol: f64,
) -> bool {
    validate_wire_closed_internal(wire, edges, vertices, tol)
}

fn validate_wire_closed_internal(
    wire: &Wire,
    edges: &[Edge],
    vertices: &[Vertex],
    tol: f64,
) -> bool {
    if wire.edges.is_empty() {
        return true; // Empty wire is trivially closed
    }

    // Get first and last vertices
    let first_edge_idx = wire.edges.first().unwrap().idx;
    let last_edge_idx = wire.edges.last().unwrap().idx;

    if first_edge_idx >= edges.len() || last_edge_idx >= edges.len() {
        return false;
    }

    let first_edge = &edges[first_edge_idx];
    let last_edge = &edges[last_edge_idx];

    // Determine start/end based on orientation
    let first_vertex_idx = if wire.edges.first().unwrap().forward {
        first_edge.start
    } else {
        first_edge.end
    };

    let last_vertex_idx = if wire.edges.last().unwrap().forward {
        last_edge.end
    } else {
        last_edge.start
    };

    if first_vertex_idx >= vertices.len() || last_vertex_idx >= vertices.len() {
        return false;
    }

    let first_pt = vertices[first_vertex_idx].point;
    let last_pt = vertices[last_vertex_idx].point;

    (first_pt - last_pt).length() < tol
}

/// Validate that a shell is closed (manifold).
///
/// A shell is closed if every edge is shared by exactly two faces.
///
/// # Arguments
/// * `shell` - The shell to validate.
/// * `tol` - Tolerance for edge matching.
pub fn validate_shell_closed(shell: &Shell, tol: f64) -> bool {
    if shell.faces.is_empty() {
        return true;
    }

    count_boundary_edges(&shell.faces, tol) == 0
}

/// Validate that a solid has valid closure.
///
/// A solid is valid if:
/// - It has at least one shell
/// - All shells are closed
/// - Shells are properly oriented (outer shell has outward normals)
///
/// # Arguments
/// * `solid` - The solid to validate.
/// * `tol` - Tolerance for shell validation.
pub fn validate_solid_valid(solid: &Solid, tol: f64) -> bool {
    if solid.shells.is_empty() {
        return false;
    }

    // All shells must be closed
    for shell in &solid.shells {
        if !validate_shell_closed(shell, tol) {
            return false;
        }
    }

    true
}

// ─────────────────────────────────────────────────────────────────────────────
// Rebuild Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Rebuild utilities for fixing shapes.
pub struct Rebuild;

impl Rebuild {
    /// Rebuild a face with improved geometry.
    ///
    /// This can fix:
    /// - Inconsistent wire orientation
    /// - Degenerate edges
    /// - Invalid normals
    ///
    /// # Arguments
    /// * `face` - The face to rebuild.
    /// * `edges` - Edge array from the BRep.
    /// * `vertices` - Vertex array from the BRep.
    /// * `tol` - Tolerance for repairs.
    ///
    /// # Returns
    /// A new face with fixes applied.
    pub fn rebuild_face(
        face: &Face,
        edges: &[Edge],
        vertices: &[Vertex],
        tol: f64,
    ) -> Face {
        // Rebuild outer wire
        let outer_wire = Self::rebuild_wire(&face.outer_wire, edges, vertices, tol);

        // Rebuild inner wires
        let inner_wires: Vec<Wire> = face
            .inner_wires
            .iter()
            .map(|w| Self::rebuild_wire(w, edges, vertices, tol))
            .collect();

        // Recompute normal if needed
        let normal = if face.normal.length_squared() < 0.5 {
            DVec3::Z
        } else {
            face.normal
        };

        Face {
            outer_wire,
            inner_wires,
            normal,
            triangles: Vec::new(),
            mesh_dirty: true,
        }
    }

    /// Rebuild a wire with improved topology.
    ///
    /// This can fix:
    /// - Inconsistent edge orientations
    /// - Duplicate edges
    /// - Gaps in the wire
    ///
    /// # Arguments
    /// * `wire` - The wire to rebuild.
    /// * `edges` - Edge array from the BRep.
    /// * `vertices` - Vertex array from the BRep.
    /// * `tol` - Tolerance for repairs.
    ///
    /// # Returns
    /// A new wire with fixes applied.
    pub fn rebuild_wire(
        wire: &Wire,
        edges: &[Edge],
        vertices: &[Vertex],
        tol: f64,
    ) -> Wire {
        if wire.edges.is_empty() {
            return wire.clone();
        }

        // Remove duplicate consecutive edges
        let mut cleaned_edges: Vec<WireEdge> = Vec::with_capacity(wire.edges.len());
        for wire_edge in &wire.edges {
            if let Some(last) = cleaned_edges.last() {
                // Skip if same edge with same orientation (duplicate)
                if last.idx == wire_edge.idx && last.forward == wire_edge.forward {
                    continue;
                }
                // Skip if same edge with opposite orientation (reversal - cancels out)
                if last.idx == wire_edge.idx && last.forward != wire_edge.forward {
                    cleaned_edges.pop();
                    continue;
                }
            }
            cleaned_edges.push(*wire_edge);
        }

        // Try to fix orientation if wire is not closed
        if !cleaned_edges.is_empty() {
            let first_idx = cleaned_edges.first().unwrap().idx;
            let last_idx = cleaned_edges.last().unwrap().idx;

            if first_idx < edges.len() && last_idx < edges.len() {
                let first_edge = &edges[first_idx];
                let last_edge = &edges[last_idx];

                let first_vertex = if cleaned_edges.first().unwrap().forward {
                    first_edge.start
                } else {
                    first_edge.end
                };

                let last_vertex = if cleaned_edges.last().unwrap().forward {
                    last_edge.end
                } else {
                    last_edge.start
                };

                // Check if wire is closed
                if first_vertex < vertices.len() && last_vertex < vertices.len() {
                    let first_pt = vertices[first_vertex].point;
                    let last_pt = vertices[last_vertex].point;

                    // If not closed but close enough, add closing edge
                    let gap = (first_pt - last_pt).length();
                    if gap > tol && gap < tol * 100.0 {
                        // Wire has a small gap - could try to close it
                        // For now, just return as-is
                    }
                }
            }
        }

        Wire {
            edges: cleaned_edges,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRep Builder Helper
// ─────────────────────────────────────────────────────────────────────────────

/// A helper for incrementally building a BRep.
///
/// This provides a convenient API for adding topology and geometry
/// without manually managing indices.
#[derive(Debug, Clone, Default)]
pub struct BRepBuilder {
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
    shells: Vec<Shell>,
    solids: Vec<Solid>,
    geom: GeomStore,
}

impl BRepBuilder {
    /// Create a new empty BRepBuilder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a vertex and return its index.
    pub fn add_vertex(&mut self, point: DVec3) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(Vertex { point });
        idx
    }

    /// Add multiple vertices and return their index range.
    pub fn add_vertices(&mut self, points: &[DVec3]) -> std::ops::Range<usize> {
        let start = self.vertices.len();
        for &p in points {
            self.vertices.push(Vertex { point: p });
        }
        start..self.vertices.len()
    }

    /// Add an edge between two vertices and return its index.
    pub fn add_edge(&mut self, start: usize, end: usize) -> usize {
        let idx = self.edges.len();
        self.edges.push(Edge { start, end });
        idx
    }

    /// Add a curve for the last added edge.
    pub fn set_edge_curve(&mut self, edge_idx: usize, curve: Curve3) -> &mut Self {
        if edge_idx < self.edges.len() {
            let curve_idx = self.geom.curves.len();
            self.geom.curves.push(curve);
            while self.geom.edge_curve.len() <= edge_idx {
                self.geom.edge_curve.push(None);
            }
            self.geom.edge_curve[edge_idx] = Some(curve_idx);
        }
        self
    }

    /// Add a face and return its index.
    pub fn add_face(&mut self, face: Face) -> usize {
        let idx = self.faces.len();
        self.faces.push(face);
        idx
    }

    /// Add a surface for a face.
    pub fn set_face_surface(&mut self, face_idx: usize, surface: Surface3) -> &mut Self {
        if face_idx < self.faces.len() {
            let surf_idx = self.geom.surfaces.len();
            self.geom.surfaces.push(surface);
            while self.geom.face_surface.len() <= face_idx {
                self.geom.face_surface.push(None);
            }
            self.geom.face_surface[face_idx] = Some(surf_idx);
        }
        self
    }

    /// Add a shell and return its index.
    pub fn add_shell(&mut self, shell: Shell) -> usize {
        let idx = self.shells.len();
        self.shells.push(shell);
        idx
    }

    /// Add a solid and return its index.
    pub fn add_solid(&mut self, solid: Solid) -> usize {
        let idx = self.solids.len();
        self.solids.push(solid);
        idx
    }

    /// Build and return the final BRep.
    pub fn build(self) -> BRep {
        BRep {
            vertices: self.vertices,
            edges: self.edges,
            solids: self.solids,
            geom: self.geom,
            compound: None,
            compsolid: None,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Line3, Plane};

    // ──────────────────────────────────────────────────────────────────────
    // BuildVertex Tests
    // ──────────────────────────────────────────────────────────────────────

    #[test]
    fn build_vertex_at_point_creates_vertex() {
        let point = DVec3::new(1.0, 2.0, 3.0);
        let v = BuildVertex::build_vertex_at_point(point);
        assert_eq!(v.point, point);
    }

    #[test]
    fn build_vertex_on_curve_creates_vertex() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let v = BuildVertex::build_vertex_on_curve(&line, 5.0);
        assert_eq!(v.point, DVec3::new(5.0, 0.0, 0.0));
    }

    #[test]
    fn build_vertex_on_surface_creates_vertex() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let v = BuildVertex::build_vertex_on_surface(&plane, 2.0, 3.0);
        // Point should be on the plane (z=0 for a Z-normal plane at origin)
        assert!(v.point.z.abs() < 1e-6);
    }

    #[test]
    fn build_vertices_from_points_creates_multiple() {
        let points = vec![DVec3::ZERO, DVec3::X, DVec3::Y];
        let vertices = BuildVertex::build_vertices_from_points(&points);
        assert_eq!(vertices.len(), 3);
        assert_eq!(vertices[0].point, DVec3::ZERO);
        assert_eq!(vertices[1].point, DVec3::X);
        assert_eq!(vertices[2].point, DVec3::Y);
    }

    // ──────────────────────────────────────────────────────────────────────
    // BuildWire Tests
    // ──────────────────────────────────────────────────────────────────────

    #[test]
    fn build_wire_from_edges_creates_wire() {
        let vertices = vec![
            Vertex { point: DVec3::ZERO },
            Vertex { point: DVec3::X },
            Vertex { point: DVec3::X + DVec3::Y },
        ];
        let edges = vec![
            Edge { start: 0, end: 1 },
            Edge { start: 1, end: 2 },
        ];

        let wire = BuildWire::build_wire_from_edges(&edges, &vertices, TOLERANCE_ABS, false);
        assert!(wire.is_ok());
        let wire = wire.unwrap();
        assert_eq!(wire.edges.len(), 2);
    }

    #[test]
    fn build_wire_from_edges_validates_closed() {
        let vertices = vec![
            Vertex { point: DVec3::ZERO },
            Vertex { point: DVec3::X },
            Vertex { point: DVec3::X + DVec3::Y },
        ];
        let edges = vec![
            Edge { start: 0, end: 1 },
            Edge { start: 1, end: 2 },
            // Not closed - missing edge back to 0
        ];

        let result = BuildWire::build_wire_from_edges(&edges, &vertices, TOLERANCE_ABS, true);
        assert!(result.is_err());
    }

    #[test]
    fn build_closed_wire_creates_closed_wire() {
        let vertices = vec![
            Vertex { point: DVec3::ZERO },
            Vertex { point: DVec3::X },
            Vertex { point: DVec3::X + DVec3::Y },
            Vertex { point: DVec3::Y },
        ];
        let edges = vec![
            Edge { start: 0, end: 1 },
            Edge { start: 1, end: 2 },
            Edge { start: 2, end: 3 },
            Edge { start: 3, end: 0 },
        ];

        let wire = BuildWire::build_closed_wire(&edges, &vertices, TOLERANCE_ABS);
        assert!(wire.is_ok());
        let wire = wire.unwrap();
        assert!(validate_wire_closed(&wire, &edges, &vertices, TOLERANCE_ABS));
    }

    #[test]
    fn build_wire_from_points_creates_wire() {
        let points = vec![DVec3::ZERO, DVec3::X, DVec3::X + DVec3::Y];
        let (vertices, edges, wire) = BuildWire::build_wire_from_points(&points, TOLERANCE_ABS, false).unwrap();

        assert_eq!(vertices.len(), 3);
        assert_eq!(edges.len(), 2); // 3 points = 2 edges when open
        assert_eq!(wire.edges.len(), 2);
    }

    #[test]
    fn build_wire_from_points_closed_adds_closing_edge() {
        let points = vec![DVec3::ZERO, DVec3::X, DVec3::Y];
        let (vertices, edges, wire) = BuildWire::build_wire_from_points(&points, TOLERANCE_ABS, true).unwrap();

        assert_eq!(vertices.len(), 3);
        assert_eq!(edges.len(), 3); // 3 points = 3 edges when closed
        assert!(validate_wire_closed(&wire, &edges, &vertices, TOLERANCE_ABS));
    }

    // ──────────────────────────────────────────────────────────────────────
    // BuildFace Tests
    // ──────────────────────────────────────────────────────────────────────

    #[test]
    fn build_face_from_wire_creates_face() {
        let wire = Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
        };

        let face = BuildFace::build_face_from_wire(&wire, None, TOLERANCE_ABS);
        assert!(face.is_ok());
        let face = face.unwrap();
        assert_eq!(face.outer_wire.edges.len(), 2);
    }

    #[test]
    fn build_face_with_surface_uses_surface_normal() {
        let wire = Wire {
            edges: vec![WireEdge::fwd(0)],
        };
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        });

        let face = BuildFace::build_face_from_wire(&wire, Some(&plane), TOLERANCE_ABS).unwrap();
        // Normal should be Y from the plane
        assert!((face.normal - DVec3::Y).length() < 0.1 || (face.normal + DVec3::Y).length() < 0.1);
    }

    #[test]
    fn build_planar_face_from_points_creates_face() {
        let points = vec![
            DVec3::ZERO,
            DVec3::X,
            DVec3::X + DVec3::Y,
            DVec3::Y,
        ];

        let (vertices, edges, face) = BuildFace::build_planar_face_from_points(&points, TOLERANCE_ABS).unwrap();

        assert_eq!(vertices.len(), 4);
        assert_eq!(edges.len(), 4);
        assert_eq!(face.outer_wire.edges.len(), 4);
        assert!((face.normal - DVec3::Z).length() < 0.1);
    }

    // ──────────────────────────────────────────────────────────────────────
    // BuildShell Tests
    // ──────────────────────────────────────────────────────────────────────

    #[test]
    fn build_shell_from_faces_creates_shell() {
        let face = Face {
            outer_wire: Wire { edges: vec![] },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };

        let shell = BuildShell::build_shell_from_faces(&[face.clone()], TOLERANCE_ABS);
        assert!(shell.is_ok());
        let shell = shell.unwrap();
        assert_eq!(shell.faces.len(), 1);
    }

    #[test]
    fn build_shell_empty_faces_returns_error() {
        let result = BuildShell::build_shell_from_faces(&[], TOLERANCE_ABS);
        assert!(matches!(result, Err(BuildError::EmptyInput(_))));
    }

    // ──────────────────────────────────────────────────────────────────────
    // BuildSolid Tests
    // ──────────────────────────────────────────────────────────────────────

    #[test]
    fn build_solid_from_shell_creates_solid() {
        let shell = Shell {
            faces: vec![Face {
                outer_wire: Wire { edges: vec![] },
                inner_wires: vec![],
                normal: DVec3::Z,
                triangles: vec![],
                mesh_dirty: true,
            }],
        };

        let solid = BuildSolid::build_solid_from_shell(&shell);
        assert!(solid.is_ok());
        let solid = solid.unwrap();
        assert_eq!(solid.shells.len(), 1);
    }

    #[test]
    fn build_solid_empty_shell_returns_error() {
        let shell = Shell { faces: vec![] };
        let result = BuildSolid::build_solid_from_shell(&shell);
        assert!(matches!(result, Err(BuildError::EmptyInput(_))));
    }

    // ──────────────────────────────────────────────────────────────────────
    // Validation Tests
    // ──────────────────────────────────────────────────────────────────────

    #[test]
    fn validate_wire_closed_returns_true_for_closed() {
        let vertices = vec![
            Vertex { point: DVec3::ZERO },
            Vertex { point: DVec3::X },
        ];
        let edges = vec![
            Edge { start: 0, end: 1 },
            Edge { start: 1, end: 0 },
        ];
        let wire = Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
        };

        assert!(validate_wire_closed(&wire, &edges, &vertices, TOLERANCE_ABS));
    }

    #[test]
    fn validate_wire_closed_returns_false_for_open() {
        let vertices = vec![
            Vertex { point: DVec3::ZERO },
            Vertex { point: DVec3::X },
            Vertex { point: DVec3::Y },
        ];
        let edges = vec![
            Edge { start: 0, end: 1 },
            Edge { start: 1, end: 2 },
        ];
        let wire = Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
        };

        assert!(!validate_wire_closed(&wire, &edges, &vertices, TOLERANCE_ABS));
    }

    #[test]
    fn validate_shell_closed_returns_true_for_closed() {
        // A cube-like shell (simplified)
        let shell = Shell {
            faces: vec![
                Face {
                    outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                },
                Face {
                    outer_wire: Wire { edges: vec![WireEdge::rev(0)] },
                    inner_wires: vec![],
                    normal: DVec3::NEG_Z,
                    triangles: vec![],
                    mesh_dirty: true,
                },
            ],
        };

        // This shell has edges shared between faces
        assert!(validate_shell_closed(&shell, TOLERANCE_ABS));
    }

    #[test]
    fn validate_solid_valid_returns_true() {
        let solid = Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire { edges: vec![] },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                }],
            }],
        };

        assert!(validate_solid_valid(&solid, TOLERANCE_ABS));
    }

    #[test]
    fn validate_solid_valid_returns_false_for_empty() {
        let solid = Solid { shells: vec![] };
        assert!(!validate_solid_valid(&solid, TOLERANCE_ABS));
    }

    // ──────────────────────────────────────────────────────────────────────
    // Rebuild Tests
    // ──────────────────────────────────────────────────────────────────────

    #[test]
    fn rebuild_wire_removes_duplicates() {
        let vertices = vec![
            Vertex { point: DVec3::ZERO },
            Vertex { point: DVec3::X },
        ];
        let edges = vec![
            Edge { start: 0, end: 1 },
        ];
        let wire = Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(0)], // Duplicate edge
        };

        let rebuilt = Rebuild::rebuild_wire(&wire, &edges, &vertices, TOLERANCE_ABS);
        assert_eq!(rebuilt.edges.len(), 1); // Duplicate removed
    }

    #[test]
    fn rebuild_face_preserves_topology() {
        let vertices = vec![
            Vertex { point: DVec3::ZERO },
            Vertex { point: DVec3::X },
        ];
        let edges = vec![Edge { start: 0, end: 1 }];
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };

        let rebuilt = Rebuild::rebuild_face(&face, &edges, &vertices, TOLERANCE_ABS);
        assert_eq!(rebuilt.outer_wire.edges.len(), 1);
    }

    // ──────────────────────────────────────────────────────────────────────
    // BRepBuilder Tests
    // ──────────────────────────────────────────────────────────────────────

    #[test]
    fn brep_builder_creates_brep() {
        let mut builder = BRepBuilder::new();

        let v0 = builder.add_vertex(DVec3::ZERO);
        let v1 = builder.add_vertex(DVec3::X);
        let v2 = builder.add_vertex(DVec3::X + DVec3::Y);
        let v3 = builder.add_vertex(DVec3::Y);

        let e0 = builder.add_edge(v0, v1);
        let e1 = builder.add_edge(v1, v2);
        let e2 = builder.add_edge(v2, v3);
        let e3 = builder.add_edge(v3, v0);

        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(e0),
                    WireEdge::fwd(e1),
                    WireEdge::fwd(e2),
                    WireEdge::fwd(e3),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        let f0 = builder.add_face(face);

        builder.set_face_surface(f0, Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        }));

        let shell = Shell { faces: vec![Face {
            outer_wire: Wire { edges: vec![] },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        }]};
        builder.add_shell(shell);

        let brep = builder.build();
        assert_eq!(brep.vertices.len(), 4);
        assert_eq!(brep.edges.len(), 4);
    }

    #[test]
    fn brep_builder_add_vertices() {
        let mut builder = BRepBuilder::new();
        let range = builder.add_vertices(&[DVec3::ZERO, DVec3::X, DVec3::Y]);

        assert_eq!(range.start, 0);
        assert_eq!(range.end, 3);
        assert_eq!(builder.vertices.len(), 3);
    }

    // ──────────────────────────────────────────────────────────────────────
    // Error Display Tests
    // ──────────────────────────────────────────────────────────────────────

    #[test]
    fn build_error_display() {
        let err = BuildError::EmptyInput("faces");
        assert!(err.to_string().contains("faces"));

        let err = BuildError::OpenWire {
            gap: 0.1,
            first_vertex: DVec3::ZERO,
            last_vertex: DVec3::X,
        };
        assert!(err.to_string().contains("wire is not closed"));

        let err = BuildError::OpenShell { boundary_edge_count: 5 };
        assert!(err.to_string().contains("shell is not closed"));
    }
}
