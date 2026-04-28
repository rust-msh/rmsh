//! BRepOffsetAPI-style offset operations — high-level API for offset, hollow, and evolved shapes.
//!
//! This module provides high-level offset operations analogous to OCCT's BRepOffsetAPI:
//!
//! - **`MakeOffset`**: Offset a wire with configurable join types
//! - **`MakeOffsetShape`**: Offset a shape (shell or solid)
//! - **`MakeThickSolid`**: Create hollow solids with specified wall thickness
//! - **`MakePipeShell`**: Create shells along a path (sweep operation)
//! - **`MakeEvolved`**: Create evolved solids from profiles
//!
//! # Overview
//!
//! The BRepOffsetAPI provides algorithms for creating offset shapes:
//!
//! 1. **Wire Offset**: Creates parallel curves at a specified distance
//! 2. **Shell/Solid Offset**: Moves all faces along their normals
//! 3. **Thick Solid**: Creates hollow solids with wall thickness
//! 4. **Pipe Shell**: Sweeps profiles along a spine curve
//! 5. **Evolved Solid**: Creates solids from profile evolution
//!
//! # Join Types
//!
//! - **Intersection**: Sharp corners at edge intersections
//! - **Arc**: Round corners using fillet arcs
//! - **Tangent**: Smooth transitions between adjacent faces
//!
//! # Offset Modes
//!
//! - **Shell**: Offset creates a shell (surfaces only)
//! - **Solid**: Offset creates a solid volume
//! - **Skin**: Offset creates a thin skin around the shape
//!
//! # Example
//!
//! ```ignore
//! use rcad_algorithms::brep_offset::{MakeOffsetShape, OffsetOptions, JoinType};
//!
//! let opts = OffsetOptions::new(0.5)
//!     .with_join_type(JoinType::Arc)
//!     .with_tolerance(1e-4);
//!
//! let offset_result = MakeOffsetShape::new(&brep, opts).build()?;
//! ```

use glam::DVec3;
use rcad_kernel::{
    BRep,
    SurfaceEval, CurveEval,
    geom::{Curve3, Surface3, Line3, Plane},
    topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge},
};

use crate::tolerance::TOLERANCE_ABS;
use crate::offset::{self, JoinType, OffsetError, OffsetOptions, OffsetResult};
use crate::triangulate::{TessellationParams, mesh_brep};

// ─────────────────────────────────────────────────────────────────────────────
// Offset Mode
// ─────────────────────────────────────────────────────────────────────────────

/// Mode for offset operations.
///
/// Determines the type of result produced by the offset operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OffsetMode {
    /// Offset creates a shell (surfaces only).
    ///
    /// The result is an open or closed shell depending on input.
    Shell,

    /// Offset creates a solid volume.
    ///
    /// For closed input shells, creates a solid with offset faces.
    /// For open shells, attempts to close with lateral faces.
    #[default]
    Solid,

    /// Offset creates a thin skin around the shape.
    ///
    /// Creates both inner and outer surfaces connected by lateral faces,
    /// resulting in a thin-walled structure.
    Skin,
}

impl OffsetMode {
    /// Returns true if the mode requires volume closure.
    pub fn requires_closure(&self) -> bool {
        matches!(self, OffsetMode::Solid | OffsetMode::Skin)
    }

    /// Returns true if the mode creates double surfaces (inner/outer).
    pub fn is_double_sided(&self) -> bool {
        matches!(self, OffsetMode::Skin)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Enhanced Offset Options
// ─────────────────────────────────────────────────────────────────────────────

/// Enhanced options for BRepOffsetAPI operations.
#[derive(Debug, Clone)]
pub struct BRepOffsetOptions {
    /// Base offset options.
    pub base: OffsetOptions,
    /// Offset mode (Shell, Solid, Skin).
    pub mode: OffsetMode,
    /// Whether to allow degenerate results.
    pub allow_degenerate: bool,
    /// Whether to perform interpolation for smooth transitions.
    pub interpolation: bool,
    /// Number of interpolation steps for smooth transitions.
    pub interpolation_steps: usize,
    /// Whether to cap open edges (for Skin mode).
    pub cap_open_edges: bool,
    /// Tolerance for geometric computations.
    pub tolerance: f64,
    /// Angular tolerance for tangent detection (radians).
    pub angular_tolerance: f64,
}

impl Default for BRepOffsetOptions {
    fn default() -> Self {
        Self {
            base: OffsetOptions::default(),
            mode: OffsetMode::default(),
            allow_degenerate: false,
            interpolation: false,
            interpolation_steps: 10,
            cap_open_edges: true,
            tolerance: TOLERANCE_ABS,
            angular_tolerance: 1e-6,
        }
    }
}

impl BRepOffsetOptions {
    /// Create options with a given distance.
    pub fn new(distance: f64) -> Self {
        Self {
            base: OffsetOptions::new(distance),
            ..Default::default()
        }
    }

    /// Set the offset mode.
    pub fn with_mode(mut self, mode: OffsetMode) -> Self {
        self.mode = mode;
        self
    }

    /// Set the join type.
    pub fn with_join_type(mut self, join_type: JoinType) -> Self {
        self.base.join_type = join_type;
        self
    }

    /// Set tolerance.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self.base.tolerance = tol;
        self
    }

    /// Enable interpolation.
    pub fn with_interpolation(mut self, steps: usize) -> Self {
        self.interpolation = true;
        self.interpolation_steps = steps;
        self
    }

    /// Set whether to cap open edges.
    pub fn with_cap_open_edges(mut self, cap: bool) -> Self {
        self.cap_open_edges = cap;
        self
    }

    /// Enable self-intersection checking.
    pub fn with_self_intersection_check(mut self, check: bool) -> Self {
        self.base.check_self_intersection = check;
        self
    }

    /// Enable auto-repair for self-intersections.
    pub fn with_auto_repair(mut self, repair: bool) -> Self {
        self.base.auto_repair = repair;
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Result Types
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a wire offset operation.
#[derive(Debug, Clone)]
pub struct WireOffsetResult {
    /// The resulting wire.
    pub wire: Wire,
    /// Vertices created for the offset wire.
    pub vertices: Vec<usize>,
    /// Edges created for the offset wire.
    pub edges: Vec<usize>,
    /// Whether the result is closed.
    pub is_closed: bool,
    /// Warnings generated during the operation.
    pub warnings: Vec<String>,
}

/// Result of a thick solid operation.
#[derive(Debug, Clone)]
pub struct ThickSolidResult {
    /// The resulting BRep.
    pub brep: BRep,
    /// Number of offset faces.
    pub offset_faces: usize,
    /// Number of lateral faces.
    pub lateral_faces: usize,
    /// Number of join faces.
    pub join_faces: usize,
    /// Whether self-intersection was detected.
    pub self_intersection: bool,
    /// Warnings generated during the operation.
    pub warnings: Vec<String>,
}

/// Result of a pipe shell operation.
#[derive(Debug, Clone)]
pub struct PipeShellResult {
    /// The resulting shell.
    pub shell: Shell,
    /// The resulting BRep.
    pub brep: BRep,
    /// Number of section faces.
    pub section_faces: usize,
    /// Number of lateral faces.
    pub lateral_faces: usize,
    /// Warnings generated during the operation.
    pub warnings: Vec<String>,
}

/// Result of an evolved solid operation.
#[derive(Debug, Clone)]
pub struct EvolvedResult {
    /// The resulting BRep.
    pub brep: BRep,
    /// Number of faces created.
    pub face_count: usize,
    /// Warnings generated during the operation.
    pub warnings: Vec<String>,
}

// ─────────────────────────────────────────────────────────────────────────────
// MakeOffset - Wire Offset
// ─────────────────────────────────────────────────────────────────────────────

/// MakeOffset - Offset a wire along its normal direction.
///
/// Creates a parallel wire at a specified distance from the original.
/// Supports different join types for handling corners.
///
/// # Join Types
///
/// - **Intersection**: Sharp corners where offset edges extend to intersect
/// - **Arc**: Round corners with fillet arcs of specified radius
/// - **Tangent**: Smooth tangent transitions between adjacent edges
pub struct MakeOffset<'a> {
    /// The input wire to offset.
    wire: &'a Wire,
    /// The BRep containing the wire's geometry.
    brep: &'a BRep,
    /// Offset distance.
    distance: f64,
    /// Join type for corners.
    join_type: JoinType,
    /// Tolerance for computations.
    tolerance: f64,
    /// Whether the wire is closed.
    is_closed: bool,
}

impl<'a> MakeOffset<'a> {
    /// Create a new wire offset operation.
    pub fn new(wire: &'a Wire, brep: &'a BRep, distance: f64) -> Self {
        Self {
            wire,
            brep,
            distance,
            join_type: JoinType::default(),
            tolerance: TOLERANCE_ABS,
            is_closed: Self::is_wire_closed(wire, brep),
        }
    }

    /// Set the join type for corners.
    pub fn with_join_type(mut self, join_type: JoinType) -> Self {
        self.join_type = join_type;
        self
    }

    /// Set tolerance for computations.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Check if the wire is closed.
    fn is_wire_closed(wire: &Wire, brep: &BRep) -> bool {
        if wire.edges.is_empty() {
            return false;
        }

        let first_edge = &brep.edges[wire.edges[0].idx];
        let last_edge = &brep.edges[wire.edges.last().unwrap().idx];

        // Check if end of last edge connects to start of first edge
        let first_start = first_edge.start;
        let last_end = if !wire.edges.last().unwrap().forward {
            last_edge.start
        } else {
            last_edge.end
        };

        first_start == last_end
    }

    /// Build the offset wire.
    pub fn build(&self) -> Result<WireOffsetResult, OffsetError> {
        if self.distance.abs() < 1e-12 {
            return Err(OffsetError::ZeroDistance);
        }

        if self.wire.edges.is_empty() {
            return Err(OffsetError::InvalidInput("wire has no edges"));
        }

        let mut result_brep = BRep::new();
        let mut offset_vertices: Vec<usize> = Vec::new();
        let mut offset_edges: Vec<usize> = Vec::new();
        let mut warnings = Vec::new();

        // Compute offset direction for each edge
        let edge_count = self.wire.edges.len();
        let mut offset_points: Vec<DVec3> = Vec::with_capacity(edge_count + 1);

        // Compute the 2D normal for the wire (assuming planar wire)
        let wire_normal = self.compute_wire_normal()?;

        // Compute offset points for each vertex
        for (i, we) in self.wire.edges.iter().enumerate() {
            let edge = &self.brep.edges[we.idx];

            // Get the curve for this edge
            let curve = self.get_edge_curve(we.idx);
            let (t0, t1) = self.get_edge_range(we.idx);

            // Compute edge tangent and normal
            let p0 = curve.point_at(t0);
            let p1 = curve.point_at(t1);

            let tangent = if !we.forward {
                (p0 - p1).normalize_or(DVec3::X)
            } else {
                (p1 - p0).normalize_or(DVec3::X)
            };

            // Offset normal is perpendicular to tangent in the wire plane
            let offset_normal = wire_normal.cross(tangent).normalize_or(DVec3::Y);

            // Get vertex position
            let vertex_idx = if !we.forward { edge.end } else { edge.start };
            let vertex_pos = self.brep.vertices[vertex_idx].point;

            // Offset the vertex
            let offset_point = vertex_pos + offset_normal * self.distance;

            // Always push the offset point
            // For closed wires, we'll add the first point again at the end to close the loop
            offset_points.push(offset_point);
        }

        // For closed wire, the last point should connect to the first
        if self.is_closed && !offset_points.is_empty() {
            offset_points.push(offset_points[0]);
        } else if !self.is_closed {
            // Add the final point
            let last_we = self.wire.edges.last().unwrap();
            let last_edge = &self.brep.edges[last_we.idx];
            let vertex_idx = if !last_we.forward { last_edge.start } else { last_edge.end };
            let vertex_pos = self.brep.vertices[vertex_idx].point;

            // Compute offset for last vertex
            let prev_we = &self.wire.edges[self.wire.edges.len() - 1];
            let _prev_edge = &self.brep.edges[prev_we.idx];
            let curve = self.get_edge_curve(prev_we.idx);
            let (t0, t1) = self.get_edge_range(prev_we.idx);

            let p0 = curve.point_at(t0);
            let p1 = curve.point_at(t1);
            let tangent = if !prev_we.forward {
                (p0 - p1).normalize_or(DVec3::X)
            } else {
                (p1 - p0).normalize_or(DVec3::X)
            };
            let offset_normal = wire_normal.cross(tangent).normalize_or(DVec3::Y);
            offset_points.push(vertex_pos + offset_normal * self.distance);
        }

        // Create vertices
        for &p in &offset_points {
            let idx = result_brep.vertices.len();
            result_brep.vertices.push(Vertex { point: p });
            offset_vertices.push(idx);
        }

        // Create edges between consecutive offset points
        for i in 0..offset_points.len() - 1 {
            let v0 = offset_vertices[i];
            let v1 = offset_vertices[i + 1];

            let p0 = result_brep.vertices[v0].point;
            let p1 = result_brep.vertices[v1].point;

            let dir = (p1 - p0).normalize_or(DVec3::X);
            let len = (p1 - p0).length();

            // Create line curve
            let curve = Curve3::Line(Line3 {
                origin: p0,
                direction: dir,
            });

            let curve_idx = result_brep.geom.curves.len();
            result_brep.geom.curves.push(curve);

            let edge_idx = result_brep.edges.len();
            result_brep.edges.push(Edge { start: v0, end: v1 });

            result_brep.geom.edge_curve.push(Some(curve_idx));
            result_brep.geom.edge_curve_range.push(Some([0.0, len]));
            result_brep.geom.edge_degenerated.push(false);

            offset_edges.push(edge_idx);
        }

        // Apply join type for corners
        if self.join_type.requires_join_geometry() && offset_edges.len() > 2 {
            self.apply_corner_joins(&mut result_brep, &offset_edges, &mut warnings);
        }

        // Build the result wire
        let wire = Wire {
            edges: offset_edges.iter().map(|&idx| WireEdge::fwd(idx)).collect(),
        };

        Ok(WireOffsetResult {
            wire,
            vertices: offset_vertices,
            edges: offset_edges,
            is_closed: self.is_closed,
            warnings,
        })
    }

    /// Compute the normal of the wire's plane.
    fn compute_wire_normal(&self) -> Result<DVec3, OffsetError> {
        // Collect edge points
        let mut points: Vec<DVec3> = Vec::new();

        for we in &self.wire.edges {
            let curve = self.get_edge_curve(we.idx);
            let (t0, t1) = self.get_edge_range(we.idx);

            points.push(curve.point_at(t0));
            points.push(curve.point_at((t0 + t1) * 0.5));
        }

        if points.len() < 3 {
            return Ok(DVec3::Z);
        }

        // Compute normal using Newell's method
        let mut normal = DVec3::ZERO;
        let n = points.len();

        for i in 0..n {
            let p0 = points[i];
            let p1 = points[(i + 1) % n];

            normal.x += (p0.y - p1.y) * (p0.z + p1.z);
            normal.y += (p0.z - p1.z) * (p0.x + p1.x);
            normal.z += (p0.x - p1.x) * (p0.y + p1.y);
        }

        if normal.length_squared() < 1e-20 {
            Ok(DVec3::Z)
        } else {
            Ok(normal.normalize())
        }
    }

    /// Get the curve for an edge.
    fn get_edge_curve(&self, edge_idx: usize) -> Curve3 {
        let curve_idx = self.brep.geom.edge_curve.get(edge_idx).and_then(|c| *c);

        match curve_idx {
            Some(idx) => self.brep.geom.curves[idx].clone(),
            None => {
                // Create line from vertex positions
                let edge = &self.brep.edges[edge_idx];
                let p0 = self.brep.vertices[edge.start].point;
                let p1 = self.brep.vertices[edge.end].point;
                Curve3::Line(Line3 {
                    origin: p0,
                    direction: (p1 - p0).normalize_or(DVec3::X),
                })
            }
        }
    }

    /// Get the parameter range for an edge.
    fn get_edge_range(&self, edge_idx: usize) -> (f64, f64) {
        let range = self.brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r);

        match range {
            Some([t0, t1]) => (t0, t1),
            None => {
                let edge = &self.brep.edges[edge_idx];
                let p0 = self.brep.vertices[edge.start].point;
                let p1 = self.brep.vertices[edge.end].point;
                (0.0, (p1 - p0).length())
            }
        }
    }

    /// Apply corner join geometry for arc/tangent joins.
    fn apply_corner_joins(
        &self,
        _result_brep: &mut BRep,
        _edges: &[usize],
        warnings: &mut Vec<String>,
    ) {
        // For arc joins, insert fillet arcs at corners
        // For tangent joins, smooth the transitions
        match self.join_type {
            JoinType::Arc => {
                warnings.push("Arc join at corners not fully implemented".to_string());
            }
            JoinType::Tangent => {
                warnings.push("Tangent join at corners not fully implemented".to_string());
            }
            JoinType::Intersection => {}
        }
    }
}

/// Offset a wire by a given distance.
///
/// # Arguments
///
/// * `wire` - The input wire
/// * `brep` - The BRep containing the wire's geometry
/// * `distance` - Offset distance (positive = right, negative = left)
/// * `join_type` - How to handle corners
///
/// # Returns
///
/// The offset wire result.
pub fn offset_wire(
    wire: &Wire,
    brep: &BRep,
    distance: f64,
    join_type: JoinType,
) -> Result<WireOffsetResult, OffsetError> {
    MakeOffset::new(wire, brep, distance)
        .with_join_type(join_type)
        .build()
}

// ─────────────────────────────────────────────────────────────────────────────
// MakeOffsetShape - Shape Offset
// ─────────────────────────────────────────────────────────────────────────────

/// MakeOffsetShape - Offset a shape (shell or solid).
///
/// Creates an offset shape by moving all faces along their normals.
/// Supports different offset modes and join types.
pub struct MakeOffsetShape<'a> {
    /// The input BRep.
    brep: &'a BRep,
    /// Offset options.
    options: BRepOffsetOptions,
}

impl<'a> MakeOffsetShape<'a> {
    /// Create a new shape offset operation.
    pub fn new(brep: &'a BRep, options: BRepOffsetOptions) -> Self {
        Self { brep, options }
    }

    /// Create with simple distance.
    pub fn from_distance(brep: &'a BRep, distance: f64) -> Self {
        Self {
            brep,
            options: BRepOffsetOptions::new(distance),
        }
    }

    /// Build the offset shape.
    pub fn build(&self) -> Result<OffsetResult, OffsetError> {
        // Use the existing offset_shape function from offset.rs
        offset::offset_shape(self.brep, self.options.base.clone())
    }
}

/// Offset a shape with the given options.
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `opts` - Offset options
///
/// # Returns
///
/// The offset result.
pub fn offset_shape_with_options(brep: &BRep, opts: BRepOffsetOptions) -> Result<OffsetResult, OffsetError> {
    MakeOffsetShape::new(brep, opts).build()
}

/// Offset a shape with a join type.
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `distance` - Offset distance
/// * `join_type` - Join type for edges
///
/// # Returns
///
/// The offset result.
pub fn offset_shape_with_join(brep: &BRep, distance: f64, join_type: JoinType) -> Result<OffsetResult, OffsetError> {
    let opts = BRepOffsetOptions::new(distance)
        .with_join_type(join_type);

    offset_shape_with_options(brep, opts)
}

// ─────────────────────────────────────────────────────────────────────────────
// MakeThickSolid - Hollow Solid
// ─────────────────────────────────────────────────────────────────────────────

/// MakeThickSolid - Create a hollow solid with specified wall thickness.
///
/// Creates a thin-walled solid by removing specified faces and offsetting
/// the remaining faces inward by the wall thickness.
pub struct MakeThickSolid<'a> {
    /// The input BRep.
    brep: &'a BRep,
    /// Wall thickness.
    thickness: f64,
    /// Faces to remove (creates openings).
    faces_to_remove: Vec<usize>,
    /// Join type for edge transitions.
    join_type: JoinType,
    /// Tolerance for computations.
    tolerance: f64,
}

impl<'a> MakeThickSolid<'a> {
    /// Create a new thick solid operation.
    pub fn new(brep: &'a BRep, thickness: f64) -> Self {
        Self {
            brep,
            thickness,
            faces_to_remove: Vec::new(),
            join_type: JoinType::default(),
            tolerance: TOLERANCE_ABS,
        }
    }

    /// Specify faces to remove (creates openings).
    pub fn with_faces_to_remove(mut self, faces: &[usize]) -> Self {
        self.faces_to_remove = faces.to_vec();
        self
    }

    /// Set the join type for edge transitions.
    pub fn with_join_type(mut self, join_type: JoinType) -> Self {
        self.join_type = join_type;
        self
    }

    /// Set tolerance for computations.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Build the thick solid.
    pub fn build(&self) -> Result<ThickSolidResult, OffsetError> {
        use std::collections::{HashMap, HashSet};

        if self.thickness <= 0.0 {
            return Err(OffsetError::InvalidInput("thickness must be positive"));
        }

        let solid = match self.brep.solids.first() {
            Some(s) => s,
            None => return Err(OffsetError::InvalidInput("BRep has no solids")),
        };

        let shell = match solid.shells.first() {
            Some(s) => s,
            None => return Err(OffsetError::InvalidInput("solid has no shells")),
        };

        let open_set: HashSet<usize> = self.faces_to_remove.iter().copied().collect();

        // Count offset faces (kept faces that will be offset)
        let offset_face_count = shell.faces.len() - open_set.len();

        // Count lateral faces by finding boundary edges
        // Boundary edges are edges shared between kept and removed faces
        // Each boundary edge creates one lateral face
        let mut edge_use: HashMap<usize, usize> = HashMap::new();
        for (fi, face) in shell.faces.iter().enumerate() {
            if open_set.contains(&fi) {
                continue;
            }
            for we in &face.outer_wire.edges {
                *edge_use.entry(we.idx).or_insert(0) += 1;
            }
        }

        // Find boundary edges (edges where one adjacent face is removed and one is kept)
        let mut lateral_face_count = 0;
        for (fi, face) in shell.faces.iter().enumerate() {
            if !open_set.contains(&fi) {
                continue;
            }
            for we in &face.outer_wire.edges {
                // Check if this edge is shared with a kept face
                let is_shared = shell.faces.iter().enumerate().any(|(fj, fj_face)| {
                    !open_set.contains(&fj)
                        && fj_face.outer_wire.edges.iter().any(|we2| we2.idx == we.idx)
                });
                if is_shared {
                    lateral_face_count += 1;
                }
            }
        }

        // Join faces are created when using Arc or Tangent join types
        // Currently, hollow_solid_with_options doesn't create join geometry,
        // so this is 0. Future implementation could count corner faces.
        let join_face_count = if self.join_type.requires_join_geometry() {
            // Count corners (vertices where boundary edges meet)
            // Each corner could potentially create a join face
            0 // Placeholder until join geometry is implemented
        } else {
            0
        };

        // Use existing hollow_solid_with_options
        let opts = OffsetOptions::new(-self.thickness)
            .with_join_type(self.join_type)
            .with_tolerance(self.tolerance);

        let result = offset::hollow_solid_with_options(
            solid,
            self.brep,
            self.thickness,
            &self.faces_to_remove,
            &opts,
        )?;

        // Check for self-intersection
        let self_intersection = offset::detect_self_intersection(&result, self.thickness);

        Ok(ThickSolidResult {
            brep: result,
            offset_faces: offset_face_count,
            lateral_faces: lateral_face_count,
            join_faces: join_face_count,
            self_intersection,
            warnings: Vec::new(),
        })
    }
}

/// Create a hollow solid by removing faces and offsetting.
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `thickness` - Wall thickness
/// * `faces_to_remove` - Indices of faces to remove (creates openings)
///
/// # Returns
///
/// The hollow solid result.
pub fn make_thick_solid(
    brep: &BRep,
    thickness: f64,
    faces_to_remove: &[usize],
) -> Result<ThickSolidResult, OffsetError> {
    MakeThickSolid::new(brep, thickness)
        .with_faces_to_remove(faces_to_remove)
        .build()
}

/// Create a hollow solid with automatic face selection.
///
/// Automatically selects the largest face for removal to create a hollow solid.
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `wall_thickness` - Wall thickness
///
/// # Returns
///
/// The hollow solid result.
pub fn make_hollow_solid(
    brep: &BRep,
    wall_thickness: f64,
) -> Result<ThickSolidResult, OffsetError> {
    // Find the largest face to remove
    let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return Err(OffsetError::InvalidInput("BRep has no shells")),
    };

    // Find largest face by vertex count (simple approximation)
    let mut largest_face_idx = 0;
    let mut max_verts = 0;

    for (i, face) in shell.faces.iter().enumerate() {
        let vert_count = face.outer_wire.edges.len();
        if vert_count > max_verts {
            max_verts = vert_count;
            largest_face_idx = i;
        }
    }

    make_thick_solid(brep, wall_thickness, &[largest_face_idx])
}

// ─────────────────────────────────────────────────────────────────────────────
// MakePipeShell - Shell Along Path
// ─────────────────────────────────────────────────────────────────────────────

/// MakePipeShell - Create a shell by sweeping profiles along a spine.
///
/// Creates a shell (or solid) by sweeping one or more profiles along
/// a spine curve. This is similar to the sweep or extrude-along-path operation.
pub struct MakePipeShell<'a> {
    /// Profile wires to sweep.
    profiles: Vec<&'a Wire>,
    /// The BRep containing the profile geometry.
    brep: &'a BRep,
    /// Spine curve for the sweep path.
    spine: &'a Wire,
    /// Whether to create a solid (vs shell).
    make_solid: bool,
    /// Number of sections along the spine.
    sections: usize,
    /// Tolerance for computations.
    tolerance: f64,
}

impl<'a> MakePipeShell<'a> {
    /// Create a new pipe shell operation.
    pub fn new(profiles: Vec<&'a Wire>, brep: &'a BRep, spine: &'a Wire) -> Self {
        Self {
            profiles,
            brep,
            spine,
            make_solid: false,
            sections: 20,
            tolerance: TOLERANCE_ABS,
        }
    }

    /// Set whether to create a solid.
    pub fn make_solid(mut self, make_solid: bool) -> Self {
        self.make_solid = make_solid;
        self
    }

    /// Set the number of sections along the spine.
    pub fn with_sections(mut self, sections: usize) -> Self {
        self.sections = sections.max(2);
        self
    }

    /// Set tolerance for computations.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Build the pipe shell.
    pub fn build(&self) -> Result<PipeShellResult, OffsetError> {
        if self.profiles.is_empty() {
            return Err(OffsetError::InvalidInput("no profiles provided"));
        }

        if self.spine.edges.is_empty() {
            return Err(OffsetError::InvalidInput("spine has no edges"));
        }

        let mut result_brep = BRep::new();
        result_brep.solids.push(Solid {
            shells: vec![Shell { faces: Vec::new() }],
        });

        let warnings = Vec::new();

        // Sample points along the spine
        let spine_points = self.sample_spine()?;

        if spine_points.len() < 2 {
            return Err(OffsetError::InvalidInput("spine has insufficient points"));
        }

        // Get the first profile
        let profile = self.profiles[0];

        // Compute profile vertices
        let profile_verts: Vec<DVec3> = profile
            .edges
            .iter()
            .map(|we| {
                let edge = &self.brep.edges[we.idx];
                let idx = if we.forward { edge.start } else { edge.end };
                self.brep.vertices[idx].point
            })
            .collect();

        // Create sections along the spine
        let mut all_section_verts: Vec<Vec<usize>> = Vec::new();

        for (i, &(origin, tangent)) in spine_points.iter().enumerate() {
            // Create transformation to move profile to this spine point
            let section_verts = self.create_section(
                &profile_verts,
                origin,
                tangent,
                &mut result_brep,
                i,
            );

            all_section_verts.push(section_verts);
        }

        // Create lateral faces between sections
        let mut lateral_faces = 0;

        for i in 0..all_section_verts.len() - 1 {
            let section0 = &all_section_verts[i];
            let section1 = &all_section_verts[i + 1];

            // Create faces between corresponding vertices
            for j in 0..section0.len() {
                let j_next = (j + 1) % section0.len();

                let v00 = section0[j];
                let v01 = section0[j_next];
                let v10 = section1[j];
                let v11 = section1[j_next];

                // Create quad face
                if self.create_quad_face(
                    &mut result_brep,
                    v00, v01, v11, v10,
                ).is_some() {
                    lateral_faces += 1;
                }
            }
        }

        // Create end caps if making a solid
        let section_faces = if self.make_solid {
            // Create start cap
            if let Some(first_section) = all_section_verts.first() {
                self.create_cap_face(&mut result_brep, first_section);
            }

            // Create end cap
            if let Some(last_section) = all_section_verts.last() {
                self.create_cap_face(&mut result_brep, last_section);
            }

            2
        } else {
            0
        };

        // Mesh the result
        mesh_brep(&mut result_brep, &TessellationParams::default());

        let shell = result_brep
            .solids
            .first()
            .and_then(|s| s.shells.first())
            .cloned()
            .unwrap_or_else(|| Shell { faces: Vec::new() });

        Ok(PipeShellResult {
            shell,
            brep: result_brep,
            section_faces,
            lateral_faces,
            warnings,
        })
    }

    /// Sample points along the spine.
    fn sample_spine(&self) -> Result<Vec<(DVec3, DVec3)>, OffsetError> {
        let mut points: Vec<(DVec3, DVec3)> = Vec::with_capacity(self.sections + 1);

        // Collect all spine curve parameters
        let mut total_length = 0.0;
        let mut segments: Vec<(f64, f64, Curve3, DVec3, DVec3)> = Vec::new();

        for we in &self.spine.edges {
            let curve = self.get_edge_curve(we.idx);
            let (t0, t1) = self.get_edge_range(we.idx);

            let p0 = curve.point_at(t0);
            let p1 = curve.point_at(t1);
            let len = (p1 - p0).length();

            segments.push((t0, t1, curve, p0, p1));
            total_length += len;
        }

        if total_length < 1e-10 {
            return Err(OffsetError::InvalidInput("spine has zero length"));
        }

        // Sample points at equal arc length intervals
        let step = total_length / self.sections as f64;
        let mut current_length = 0.0;
        let mut seg_idx = 0;
        let mut seg_remaining = (segments[0].4 - segments[0].3).length();

        for i in 0..=self.sections {
            let target_length = i as f64 * step;

            // Find the segment containing this point
            while current_length + seg_remaining < target_length - 1e-10 && seg_idx < segments.len() - 1 {
                current_length += seg_remaining;
                seg_idx += 1;
                seg_remaining = (segments[seg_idx].4 - segments[seg_idx].3).length();
            }

            let seg = &segments[seg_idx];
            let seg_progress = (target_length - current_length) / seg_remaining.max(1e-10);
            let t = seg.0 + seg_progress * (seg.1 - seg.0);

            let point = seg.2.point_at(t);
            let tangent = seg.2.tangent_at(t).normalize_or(DVec3::Z);

            points.push((point, tangent));
        }

        Ok(points)
    }

    /// Create a profile section at a spine point.
    fn create_section(
        &self,
        profile_verts: &[DVec3],
        origin: DVec3,
        tangent: DVec3,
        result_brep: &mut BRep,
        _section_idx: usize,
    ) -> Vec<usize> {
        // Compute transformation
        let z_axis = tangent;
        let x_axis = if z_axis.cross(DVec3::X).length() > 1e-6 {
            z_axis.cross(DVec3::X).normalize()
        } else {
            z_axis.cross(DVec3::Y).normalize()
        };
        let y_axis = z_axis.cross(x_axis).normalize();

        // Compute profile centroid
        let centroid = profile_verts.iter().fold(DVec3::ZERO, |acc, &p| acc + p)
            / profile_verts.len() as f64;

        // Create vertices
        let mut section_verts = Vec::with_capacity(profile_verts.len());

        for &p in profile_verts {
            // Translate profile point relative to centroid
            let local = p - centroid;

            // Transform to spine location
            let transformed = origin + x_axis * local.x + y_axis * local.y + z_axis * local.z;

            let idx = result_brep.vertices.len();
            result_brep.vertices.push(Vertex { point: transformed });
            section_verts.push(idx);
        }

        section_verts
    }

    /// Create a quad face.
    fn create_quad_face(
        &self,
        result_brep: &mut BRep,
        v0: usize,
        v1: usize,
        v2: usize,
        v3: usize,
    ) -> Option<usize> {
        let p0 = result_brep.vertices[v0].point;
        let p1 = result_brep.vertices[v1].point;
        let _p2 = result_brep.vertices[v2].point;
        let p3 = result_brep.vertices[v3].point;

        // Compute face normal
        let e1 = p1 - p0;
        let e2 = p3 - p0;
        let normal = e1.cross(e2).normalize_or(DVec3::Z);

        // Create edges
        let mut wire_edges = Vec::new();
        let verts = [v0, v1, v2, v3];

        for i in 0..4 {
            let start = verts[i];
            let end = verts[(i + 1) % 4];

            let sp = result_brep.vertices[start].point;
            let ep = result_brep.vertices[end].point;

            let dir = (ep - sp).normalize_or(DVec3::X);
            let len = (ep - sp).length();

            let curve = Curve3::Line(Line3 {
                origin: sp,
                direction: dir,
            });

            let curve_idx = result_brep.geom.curves.len();
            result_brep.geom.curves.push(curve);

            let edge_idx = result_brep.edges.len();
            result_brep.edges.push(Edge { start, end });

            result_brep.geom.edge_curve.push(Some(curve_idx));
            result_brep.geom.edge_curve_range.push(Some([0.0, len]));
            result_brep.geom.edge_degenerated.push(false);

            wire_edges.push(WireEdge::fwd(edge_idx));
        }

        // Create plane surface
        let surface = Surface3::Plane(Plane {
            origin: p0,
            normal,
        });

        let surface_idx = result_brep.geom.surfaces.len();
        result_brep.geom.surfaces.push(surface);

        // Create face
        let face_idx = result_brep.solids[0].shells[0].faces.len();
        result_brep.solids[0].shells[0].faces.push(Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: Vec::new(),
            normal,
            triangles: Vec::new(),
            mesh_dirty: true,
        });

        result_brep.geom.face_surface.push(Some(surface_idx));

        Some(face_idx)
    }

    /// Create a cap face for the end of the pipe.
    fn create_cap_face(&self, result_brep: &mut BRep, section_verts: &[usize]) {
        if section_verts.len() < 3 {
            return;
        }

        // Compute centroid and normal
        let centroid = section_verts
            .iter()
            .fold(DVec3::ZERO, |acc, &v| acc + result_brep.vertices[v].point)
            / section_verts.len() as f64;

        // Use first three vertices to compute normal
        let p0 = result_brep.vertices[section_verts[0]].point;
        let p1 = result_brep.vertices[section_verts[1]].point;
        let p2 = result_brep.vertices[section_verts[2]].point;

        let normal = (p1 - p0).cross(p2 - p0).normalize_or(DVec3::Z);

        // Create fan triangulation
        let mut wire_edges = Vec::new();

        for i in 0..section_verts.len() {
            let start = section_verts[i];
            let end = section_verts[(i + 1) % section_verts.len()];

            let sp = result_brep.vertices[start].point;
            let ep = result_brep.vertices[end].point;

            let dir = (ep - sp).normalize_or(DVec3::X);
            let len = (ep - sp).length();

            let curve = Curve3::Line(Line3 {
                origin: sp,
                direction: dir,
            });

            let curve_idx = result_brep.geom.curves.len();
            result_brep.geom.curves.push(curve);

            let edge_idx = result_brep.edges.len();
            result_brep.edges.push(Edge { start, end });

            result_brep.geom.edge_curve.push(Some(curve_idx));
            result_brep.geom.edge_curve_range.push(Some([0.0, len]));
            result_brep.geom.edge_degenerated.push(false);

            wire_edges.push(WireEdge::fwd(edge_idx));
        }

        // Create surface
        let surface = Surface3::Plane(Plane {
            origin: centroid,
            normal,
        });

        let surface_idx = result_brep.geom.surfaces.len();
        result_brep.geom.surfaces.push(surface);

        result_brep.solids[0].shells[0].faces.push(Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: Vec::new(),
            normal,
            triangles: Vec::new(),
            mesh_dirty: true,
        });

        result_brep.geom.face_surface.push(Some(surface_idx));
    }

    /// Get the curve for an edge.
    fn get_edge_curve(&self, edge_idx: usize) -> Curve3 {
        let curve_idx = self.brep.geom.edge_curve.get(edge_idx).and_then(|c| *c);

        match curve_idx {
            Some(idx) => self.brep.geom.curves[idx].clone(),
            None => {
                let edge = &self.brep.edges[edge_idx];
                let p0 = self.brep.vertices[edge.start].point;
                let p1 = self.brep.vertices[edge.end].point;
                Curve3::Line(Line3 {
                    origin: p0,
                    direction: (p1 - p0).normalize_or(DVec3::X),
                })
            }
        }
    }

    /// Get the parameter range for an edge.
    fn get_edge_range(&self, edge_idx: usize) -> (f64, f64) {
        let range = self.brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r);

        match range {
            Some([t0, t1]) => (t0, t1),
            None => {
                let edge = &self.brep.edges[edge_idx];
                let p0 = self.brep.vertices[edge.start].point;
                let p1 = self.brep.vertices[edge.end].point;
                (0.0, (p1 - p0).length())
            }
        }
    }
}

/// Create a pipe shell by sweeping a profile along a spine.
///
/// # Arguments
///
/// * `profiles` - Profile wires to sweep
/// * `brep` - The BRep containing profile geometry
/// * `spine` - The spine wire for the sweep path
///
/// # Returns
///
/// The pipe shell result.
pub fn make_pipe_shell(
    profiles: &[&Wire],
    brep: &BRep,
    spine: &Wire,
) -> Result<PipeShellResult, OffsetError> {
    MakePipeShell::new(profiles.to_vec(), brep, spine).build()
}

// ─────────────────────────────────────────────────────────────────────────────
// MakeEvolved - Evolved Profile
// ─────────────────────────────────────────────────────────────────────────────

/// MakeEvolved - Create an evolved solid from a profile and spine.
///
/// Creates a solid by "evolving" a profile along a spine path.
/// This is similar to pipe shell but with additional profile transformation
/// options and solid generation.
pub struct MakeEvolved<'a> {
    /// The profile wire.
    profile: &'a Wire,
    /// The spine wire.
    spine: &'a Wire,
    /// The BRep containing geometry.
    brep: &'a BRep,
    /// Whether the profile should rotate to follow the spine.
    follow_spine: bool,
    /// Whether to join profile end to start (for closed profiles).
    join: bool,
    /// Number of sections along the spine.
    sections: usize,
    /// Tolerance for computations.
    tolerance: f64,
}

impl<'a> MakeEvolved<'a> {
    /// Create a new evolved solid operation.
    pub fn new(profile: &'a Wire, spine: &'a Wire, brep: &'a BRep) -> Self {
        Self {
            profile,
            spine,
            brep,
            follow_spine: true,
            join: true,
            sections: 20,
            tolerance: TOLERANCE_ABS,
        }
    }

    /// Set whether the profile follows the spine tangent.
    pub fn follow_spine(mut self, follow: bool) -> Self {
        self.follow_spine = follow;
        self
    }

    /// Set whether to join profile end to start.
    pub fn with_join(mut self, join: bool) -> Self {
        self.join = join;
        self
    }

    /// Set the number of sections.
    pub fn with_sections(mut self, sections: usize) -> Self {
        self.sections = sections.max(2);
        self
    }

    /// Build the evolved solid.
    pub fn build(&self) -> Result<EvolvedResult, OffsetError> {
        if self.profile.edges.is_empty() {
            return Err(OffsetError::InvalidInput("profile has no edges"));
        }

        if self.spine.edges.is_empty() {
            return Err(OffsetError::InvalidInput("spine has no edges"));
        }

        let mut warnings = Vec::new();

        // Use MakePipeShell for the basic construction
        let pipe_result = MakePipeShell::new(vec![self.profile], self.brep, self.spine)
            .make_solid(true)
            .with_sections(self.sections)
            .with_tolerance(self.tolerance)
            .build()?;

        let face_count = pipe_result
            .brep
            .solids
            .first()
            .and_then(|s| s.shells.first())
            .map(|sh| sh.faces.len())
            .unwrap_or(0);

        warnings.extend(pipe_result.warnings);

        Ok(EvolvedResult {
            brep: pipe_result.brep,
            face_count,
            warnings,
        })
    }
}

/// Create an evolved solid from a profile and spine.
///
/// # Arguments
///
/// * `profile` - The profile wire
/// * `spine` - The spine wire
/// * `brep` - The BRep containing geometry
///
/// # Returns
///
/// The evolved solid result.
pub fn make_evolved(
    profile: &Wire,
    spine: &Wire,
    brep: &BRep,
) -> Result<EvolvedResult, OffsetError> {
    MakeEvolved::new(profile, spine, brep).build()
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Helper to add a vertex to a BRep.
fn add_vertex(brep: &mut BRep, point: DVec3) -> usize {
    let idx = brep.vertices.len();
    brep.vertices.push(Vertex { point });
    idx
}

/// Helper to add an edge to a BRep.
fn add_edge(brep: &mut BRep, curve: Curve3, t0: f64, t1: f64, v0: usize, v1: usize) -> usize {
    let idx = brep.edges.len();
    brep.edges.push(Edge { start: v0, end: v1 });

    let ci = brep.geom.curves.len();
    brep.geom.curves.push(curve);

    while brep.geom.edge_curve.len() <= idx {
        brep.geom.edge_curve.push(None);
    }
    while brep.geom.edge_curve_range.len() <= idx {
        brep.geom.edge_curve_range.push(None);
    }
    while brep.geom.edge_degenerated.len() <= idx {
        brep.geom.edge_degenerated.push(false);
    }

    brep.geom.edge_curve[idx] = Some(ci);
    brep.geom.edge_curve_range[idx] = Some([t0, t1]);
    idx
}

/// Helper to add a face to a BRep.
fn add_face(brep: &mut BRep, surface: Surface3, outer: Wire, inner: Vec<Wire>) -> usize {
    if brep.solids.is_empty() {
        brep.solids.push(Solid {
            shells: vec![Shell { faces: Vec::new() }],
        });
    }
    if brep.solids[0].shells.is_empty() {
        brep.solids[0].shells.push(Shell { faces: Vec::new() });
    }

    let idx = brep.solids[0].shells[0].faces.len();
    let normal = surface.normal_at(0.0, 0.0);

    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: outer,
        inner_wires: inner,
        normal,
        triangles: Vec::new(),
        mesh_dirty: true,
    });

    while brep.geom.face_surface.len() <= idx {
        brep.geom.face_surface.push(None);
    }

    let si = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surface);
    brep.geom.face_surface[idx] = Some(si);

    idx
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::PrimitiveSolid;

    fn create_box_brep() -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        crate::geom_populate::populate_box_geom(&mut brep);
        brep
    }

    #[test]
    fn offset_mode_default() {
        assert_eq!(OffsetMode::default(), OffsetMode::Solid);
    }

    #[test]
    fn offset_mode_properties() {
        assert!(OffsetMode::Solid.requires_closure());
        assert!(!OffsetMode::Shell.requires_closure());
        assert!(OffsetMode::Skin.is_double_sided());
        assert!(!OffsetMode::Solid.is_double_sided());
    }

    #[test]
    fn brep_offset_options_builder() {
        let opts = BRepOffsetOptions::new(0.5)
            .with_mode(OffsetMode::Shell)
            .with_join_type(JoinType::Arc)
            .with_tolerance(1e-6)
            .with_interpolation(10);

        assert_eq!(opts.base.distance, 0.5);
        assert_eq!(opts.mode, OffsetMode::Shell);
        assert_eq!(opts.base.join_type, JoinType::Arc);
        assert!(opts.interpolation);
        assert_eq!(opts.interpolation_steps, 10);
    }

    #[test]
    fn make_offset_shape_simple() {
        let brep = create_box_brep();

        let result = MakeOffsetShape::from_distance(&brep, 0.1).build();

        assert!(result.is_ok(), "MakeOffsetShape should succeed");
        let offset_result = result.unwrap();
        assert_eq!(offset_result.offset_faces, 6);
    }

    #[test]
    fn test_offset_shape_with_join() {
        let brep = create_box_brep();

        let result = offset_shape_with_join(&brep, 0.1, JoinType::Arc);

        assert!(result.is_ok(), "offset_shape_with_join should succeed");
    }

    #[test]
    fn make_thick_solid_simple() {
        let brep = create_box_brep();

        // Remove top face (index 5)
        let result = MakeThickSolid::new(&brep, 0.1)
            .with_faces_to_remove(&[5])
            .build();

        assert!(result.is_ok(), "MakeThickSolid should succeed");
        let thick_result = result.unwrap();
        assert!(thick_result.offset_faces > 0);
    }

    #[test]
    fn make_thick_solid_zero_thickness_error() {
        let brep = create_box_brep();

        let result = MakeThickSolid::new(&brep, 0.0).build();

        assert!(result.is_err(), "zero thickness should error");
    }

    #[test]
    fn test_make_hollow_solid() {
        let brep = create_box_brep();

        let result = make_hollow_solid(&brep, 0.1);

        assert!(result.is_ok(), "make_hollow_solid should succeed");
    }

    #[test]
    fn make_thick_solid_api() {
        let brep = create_box_brep();

        let result = make_thick_solid(&brep, 0.2, &[0, 5]);

        assert!(result.is_ok(), "make_thick_solid should succeed with multiple faces");
        let thick_result = result.unwrap();
        assert!(!thick_result.brep.vertices.is_empty());
    }

    #[test]
    fn make_offset_wire_simple() {
        let brep = create_box_brep();

        // Get a face's wire
        let shell = &brep.solids[0].shells[0];
        let face = &shell.faces[0];
        let wire = &face.outer_wire;

        // Create offset wire
        let result = offset_wire(wire, &brep, 0.1, JoinType::Intersection);

        assert!(result.is_ok(), "offset_wire should succeed");
        let wire_result = result.unwrap();
        assert!(wire_result.wire.edges.len() >= 4, "expected at least 4 edges, got {}", wire_result.wire.edges.len());
    }

    #[test]
    fn make_offset_wire_closed() {
        let brep = create_box_brep();

        let shell = &brep.solids[0].shells[0];
        let face = &shell.faces[0];
        let wire = &face.outer_wire;

        let maker = MakeOffset::new(wire, &brep, 0.1);
        assert!(maker.is_closed, "box face wire should be closed");
    }

    #[test]
    fn make_pipe_shell_line_spine() {
        let _brep = create_box_brep();

        // Create a simple spine (line along Z)
        let mut spine_brep = BRep::new();
        let v0 = add_vertex(&mut spine_brep, DVec3::new(0.0, 0.0, 0.0));
        let v1 = add_vertex(&mut spine_brep, DVec3::new(0.0, 0.0, 2.0));

        let curve = Curve3::Line(Line3 {
            origin: DVec3::new(0.0, 0.0, 0.0),
            direction: DVec3::Z,
        });

        let e0 = add_edge(&mut spine_brep, curve, 0.0, 2.0, v0, v1);
        let spine = Wire {
            edges: vec![WireEdge::fwd(e0)],
        };

        // Create a simple profile (square)
        let mut profile_brep = BRep::new();
        let v0 = add_vertex(&mut profile_brep, DVec3::new(-0.5, -0.5, 0.0));
        let v1 = add_vertex(&mut profile_brep, DVec3::new(0.5, -0.5, 0.0));
        let v2 = add_vertex(&mut profile_brep, DVec3::new(0.5, 0.5, 0.0));
        let v3 = add_vertex(&mut profile_brep, DVec3::new(-0.5, 0.5, 0.0));

        let e0 = add_edge(&mut profile_brep, Curve3::Line(Line3 { origin: DVec3::new(-0.5, -0.5, 0.0), direction: DVec3::X }), 0.0, 1.0, v0, v1);
        let e1 = add_edge(&mut profile_brep, Curve3::Line(Line3 { origin: DVec3::new(0.5, -0.5, 0.0), direction: DVec3::Y }), 0.0, 1.0, v1, v2);
        let e2 = add_edge(&mut profile_brep, Curve3::Line(Line3 { origin: DVec3::new(0.5, 0.5, 0.0), direction: DVec3::NEG_X }), 0.0, 1.0, v2, v3);
        let e3 = add_edge(&mut profile_brep, Curve3::Line(Line3 { origin: DVec3::new(-0.5, 0.5, 0.0), direction: DVec3::NEG_Y }), 0.0, 1.0, v3, v0);

        let profile = Wire {
            edges: vec![WireEdge::fwd(e0), WireEdge::fwd(e1), WireEdge::fwd(e2), WireEdge::fwd(e3)],
        };

        // Build pipe shell
        let result = MakePipeShell::new(vec![&profile], &profile_brep, &spine)
            .with_sections(10)
            .build();

        assert!(result.is_ok(), "make_pipe_shell should succeed");
        let pipe_result = result.unwrap();
        assert!(pipe_result.lateral_faces > 0);
    }

    #[test]
    fn make_evolved_simple() {
        let _brep = create_box_brep();

        // Create a simple spine
        let mut spine_brep = BRep::new();
        let v0 = add_vertex(&mut spine_brep, DVec3::new(0.0, 0.0, 0.0));
        let v1 = add_vertex(&mut spine_brep, DVec3::new(0.0, 0.0, 1.0));

        let curve = Curve3::Line(Line3 {
            origin: DVec3::new(0.0, 0.0, 0.0),
            direction: DVec3::Z,
        });

        let e0 = add_edge(&mut spine_brep, curve, 0.0, 1.0, v0, v1);
        let spine = Wire {
            edges: vec![WireEdge::fwd(e0)],
        };

        // Create a simple profile
        let mut profile_brep = BRep::new();
        let v0 = add_vertex(&mut profile_brep, DVec3::new(-0.25, -0.25, 0.0));
        let v1 = add_vertex(&mut profile_brep, DVec3::new(0.25, -0.25, 0.0));
        let v2 = add_vertex(&mut profile_brep, DVec3::new(0.25, 0.25, 0.0));
        let v3 = add_vertex(&mut profile_brep, DVec3::new(-0.25, 0.25, 0.0));

        let e0 = add_edge(&mut profile_brep, Curve3::Line(Line3 { origin: DVec3::new(-0.25, -0.25, 0.0), direction: DVec3::X }), 0.0, 0.5, v0, v1);
        let e1 = add_edge(&mut profile_brep, Curve3::Line(Line3 { origin: DVec3::new(0.25, -0.25, 0.0), direction: DVec3::Y }), 0.0, 0.5, v1, v2);
        let e2 = add_edge(&mut profile_brep, Curve3::Line(Line3 { origin: DVec3::new(0.25, 0.25, 0.0), direction: DVec3::NEG_X }), 0.0, 0.5, v2, v3);
        let e3 = add_edge(&mut profile_brep, Curve3::Line(Line3 { origin: DVec3::new(-0.25, 0.25, 0.0), direction: DVec3::NEG_Y }), 0.0, 0.5, v3, v0);

        let profile = Wire {
            edges: vec![WireEdge::fwd(e0), WireEdge::fwd(e1), WireEdge::fwd(e2), WireEdge::fwd(e3)],
        };

        // Build evolved solid
        let result = MakeEvolved::new(&profile, &spine, &profile_brep)
            .with_sections(5)
            .build();

        assert!(result.is_ok(), "make_evolved should succeed");
        let evolved_result = result.unwrap();
        assert!(evolved_result.face_count > 0);
    }

    #[test]
    fn offset_wire_with_arc_join() {
        let brep = create_box_brep();

        let shell = &brep.solids[0].shells[0];
        let face = &shell.faces[0];
        let wire = &face.outer_wire;

        let result = offset_wire(wire, &brep, 0.1, JoinType::Arc);

        assert!(result.is_ok(), "offset_wire with arc join should succeed");
    }

    #[test]
    fn offset_wire_with_tangent_join() {
        let brep = create_box_brep();

        let shell = &brep.solids[0].shells[0];
        let face = &shell.faces[0];
        let wire = &face.outer_wire;

        let result = offset_wire(wire, &brep, 0.1, JoinType::Tangent);

        assert!(result.is_ok(), "offset_wire with tangent join should succeed");
    }

    #[test]
    fn offset_wire_zero_distance_error() {
        let brep = create_box_brep();

        let shell = &brep.solids[0].shells[0];
        let face = &shell.faces[0];
        let wire = &face.outer_wire;

        let result = offset_wire(wire, &brep, 0.0, JoinType::Intersection);

        assert!(result.is_err(), "offset_wire with zero distance should error");
    }

    #[test]
    fn thick_solid_face_counting_single_face_removed() {
        let brep = create_box_brep();

        // A box has 6 faces, each face has 4 edges
        // Removing 1 face should result in:
        // - 5 offset faces (kept faces)
        // - 4 lateral faces (4 boundary edges of the removed face)
        // - 0 join faces (intersection join)
        let result = MakeThickSolid::new(&brep, 0.1)
            .with_faces_to_remove(&[5])
            .with_join_type(JoinType::Intersection)
            .build();

        assert!(result.is_ok(), "MakeThickSolid should succeed");
        let thick_result = result.unwrap();
        assert_eq!(thick_result.offset_faces, 5, "should have 5 offset faces");
        assert_eq!(thick_result.lateral_faces, 4, "should have 4 lateral faces (one per boundary edge)");
        assert_eq!(thick_result.join_faces, 0, "intersection join should have 0 join faces");
    }

    #[test]
    fn thick_solid_face_counting_multiple_faces_removed() {
        let brep = create_box_brep();

        // Remove 2 opposite faces (e.g., top and bottom, indices 4 and 5)
        // Each removed face has 4 boundary edges
        let result = MakeThickSolid::new(&brep, 0.1)
            .with_faces_to_remove(&[4, 5])
            .build();

        assert!(result.is_ok(), "MakeThickSolid should succeed with multiple faces");
        let thick_result = result.unwrap();
        assert_eq!(thick_result.offset_faces, 4, "should have 4 offset faces");
        // 4 boundary edges per removed face = 8 lateral faces
        assert_eq!(thick_result.lateral_faces, 8, "should have 8 lateral faces (4 per removed face)");
        assert_eq!(thick_result.join_faces, 0, "intersection join should have 0 join faces");
    }

    #[test]
    fn thick_solid_face_counting_with_arc_join() {
        let brep = create_box_brep();

        let result = MakeThickSolid::new(&brep, 0.1)
            .with_faces_to_remove(&[5])
            .with_join_type(JoinType::Arc)
            .build();

        assert!(result.is_ok(), "MakeThickSolid with arc join should succeed");
        let thick_result = result.unwrap();
        assert_eq!(thick_result.offset_faces, 5, "should have 5 offset faces");
        assert_eq!(thick_result.lateral_faces, 4, "should have 4 lateral faces");
        // Arc join could create join faces at corners, but currently not implemented
        // This test documents the current behavior
        assert_eq!(thick_result.join_faces, 0, "arc join faces not yet implemented");
    }

    #[test]
    fn thick_solid_face_counting_with_tangent_join() {
        let brep = create_box_brep();

        let result = MakeThickSolid::new(&brep, 0.1)
            .with_faces_to_remove(&[5])
            .with_join_type(JoinType::Tangent)
            .build();

        assert!(result.is_ok(), "MakeThickSolid with tangent join should succeed");
        let thick_result = result.unwrap();
        assert_eq!(thick_result.offset_faces, 5, "should have 5 offset faces");
        assert_eq!(thick_result.lateral_faces, 4, "should have 4 lateral faces");
        // Tangent join could create join faces, but currently not implemented
        assert_eq!(thick_result.join_faces, 0, "tangent join faces not yet implemented");
    }

    #[test]
    fn make_hollow_solid_face_counting() {
        let brep = create_box_brep();

        // make_hollow_solid automatically removes the largest face
        let result = make_hollow_solid(&brep, 0.1);

        assert!(result.is_ok(), "make_hollow_solid should succeed");
        let thick_result = result.unwrap();
        // The largest face is removed, leaving 5 offset faces
        assert_eq!(thick_result.offset_faces, 5, "should have 5 offset faces");
        // The removed face has 4 boundary edges
        assert_eq!(thick_result.lateral_faces, 4, "should have 4 lateral faces");
        assert_eq!(thick_result.join_faces, 0, "should have 0 join faces");
    }
}
