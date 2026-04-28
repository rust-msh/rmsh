//! ShapeExtend-style shape extension utilities.
//!
//! Analogous to OCCT `ShapeExtend` package providing extended data structures
//! for shape analysis and repair:
//! - `ShapeExtend_WireData`: Extended wire data with edge management
//! - `ShapeExtend_CompositeSurface`: Composite surface made of patches
//! - `ShapeExtend_BasicMsgRegistrator`: Basic message registration
//! - `ShapeExtend_MsgRegistrator`: Message registration with shape context
//! - `ShapeExtend_Explorer`: Extended shape exploration utilities

use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topology::{Face, Shell, Wire};
use rcad_kernel::BRep;

use crate::brep_tools::ShapeType;

// ─────────────────────────────────────────────────────────────────────────────
// ShapeExtend_WireData - Extended Wire Data
// ─────────────────────────────────────────────────────────────────────────────

/// Extended wire data structure for wire manipulation.
///
/// Analogous to OCCT `ShapeExtend_WireData`, this provides a mutable
/// wire representation that supports edge insertion, removal, and
/// ordering operations.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::shape_extend::WireData;
///
/// let mut wire = WireData::new();
/// wire.add_edge(0, true);  // Add edge 0 with forward orientation
/// wire.add_edge(1, false); // Add edge 1 with reversed orientation
/// assert_eq!(wire.edge_count(), 2);
/// ```
#[derive(Debug, Clone)]
pub struct WireData {
    /// List of (edge_idx, orientation) pairs.
    /// orientation: true = forward, false = reversed.
    edges: Vec<(usize, bool)>,
    /// Cached wire length (sum of edge lengths).
    cached_length: Option<f64>,
    /// Flag indicating if the wire is closed.
    closed: Option<bool>,
}

impl Default for WireData {
    fn default() -> Self {
        Self::new()
    }
}

impl WireData {
    /// Create a new empty wire data.
    pub fn new() -> Self {
        WireData {
            edges: Vec::new(),
            cached_length: None,
            closed: None,
        }
    }

    /// Create a wire data with pre-allocated capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        WireData {
            edges: Vec::with_capacity(capacity),
            cached_length: None,
            closed: None,
        }
    }

    /// Add an edge to the end of the wire.
    ///
    /// # Arguments
    /// * `edge_idx` - Index of the edge in the BRep.
    /// * `orientation` - true for forward, false for reversed.
    pub fn add_edge(&mut self, edge_idx: usize, orientation: bool) {
        self.edges.push((edge_idx, orientation));
        self.cached_length = None;
        self.closed = None;
    }

    /// Add an edge at a specific position in the wire.
    ///
    /// # Arguments
    /// * `pos` - Position to insert at (0-indexed).
    /// * `edge_idx` - Index of the edge in the BRep.
    /// * `orientation` - true for forward, false for reversed.
    ///
    /// # Panics
    /// Panics if `pos` is greater than the current edge count.
    pub fn add_edge_at(&mut self, pos: usize, edge_idx: usize, orientation: bool) {
        if pos > self.edges.len() {
            panic!(
                "Position {} out of bounds (wire has {} edges)",
                pos,
                self.edges.len()
            );
        }
        self.edges.insert(pos, (edge_idx, orientation));
        self.cached_length = None;
        self.closed = None;
    }

    /// Remove an edge at the specified position.
    ///
    /// # Arguments
    /// * `pos` - Position of the edge to remove (0-indexed).
    ///
    /// # Returns
    /// The removed (edge_idx, orientation) pair.
    ///
    /// # Panics
    /// Panics if `pos` is out of bounds.
    pub fn remove_edge(&mut self, pos: usize) -> (usize, bool) {
        if pos >= self.edges.len() {
            panic!(
                "Position {} out of bounds (wire has {} edges)",
                pos,
                self.edges.len()
            );
        }
        let removed = self.edges.remove(pos);
        self.cached_length = None;
        self.closed = None;
        removed
    }

    /// Get the edge at the specified position.
    ///
    /// # Arguments
    /// * `pos` - Position of the edge (0-indexed).
    ///
    /// # Returns
    /// The (edge_idx, orientation) pair, or None if out of bounds.
    pub fn edge_at(&self, pos: usize) -> Option<(usize, bool)> {
        self.edges.get(pos).copied()
    }

    /// Get a slice of all edges.
    ///
    /// Returns `&[(edge_idx, orientation)]`.
    pub fn edges(&self) -> &[(usize, bool)] {
        &self.edges
    }

    /// Get the number of edges in the wire.
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Check if the wire has no edges.
    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }

    /// Set the orientation of an edge at the specified position.
    ///
    /// # Arguments
    /// * `pos` - Position of the edge.
    /// * `orientation` - New orientation (true = forward, false = reversed).
    pub fn set_orientation(&mut self, pos: usize, orientation: bool) {
        if pos < self.edges.len() {
            self.edges[pos].1 = orientation;
            self.closed = None;
        }
    }

    /// Reverse the orientation of all edges.
    ///
    /// This also reverses the order of edges to maintain wire continuity.
    pub fn reverse(&mut self) {
        self.edges.reverse();
        for edge in &mut self.edges {
            edge.1 = !edge.1;
        }
        self.closed = None;
    }

    /// Check if the wire forms a closed loop.
    ///
    /// For a wire to be closed, it must have at least one edge
    /// and the edges must form a continuous chain.
    ///
    /// Note: This method returns the cached value if available.
    /// Use `is_closed_with_brep` for accurate checking with topology.
    pub fn is_closed(&self) -> bool {
        self.closed.unwrap_or_else(|| self.edges.len() > 0)
    }

    /// Check if the wire is closed using BRep topology.
    ///
    /// # Arguments
    /// * `brep` - The BRep containing the edge topology.
    ///
    /// # Returns
    /// true if the wire forms a closed loop.
    pub fn is_closed_with_brep(&self, brep: &BRep) -> bool {
        if self.edges.is_empty() {
            return false;
        }

        // Get the first and last vertices
        let first_edge_idx = self.edges[0].0;
        let last_edge_idx = self.edges[self.edges.len() - 1].0;

        if first_edge_idx >= brep.edges.len() || last_edge_idx >= brep.edges.len() {
            return false;
        }

        let first_edge = &brep.edges[first_edge_idx];
        let last_edge = &brep.edges[last_edge_idx];

        let first_orient = self.edges[0].1;
        let last_orient = self.edges[self.edges.len() - 1].1;

        // First vertex of wire
        let first_vertex = if first_orient {
            first_edge.start
        } else {
            first_edge.end
        };

        // Last vertex of wire
        let last_vertex = if last_orient {
            last_edge.end
        } else {
            last_edge.start
        };

        first_vertex == last_vertex
    }

    /// Compute the total length of the wire.
    ///
    /// For a wire with no edges, returns 0.0.
    /// This is a placeholder that returns edge count * 1.0 for now,
    /// as actual length computation requires curve geometry.
    pub fn length(&self) -> f64 {
        self.cached_length.unwrap_or(0.0)
    }

    /// Compute the wire length using BRep geometry.
    ///
    /// # Arguments
    /// * `brep` - The BRep containing edge geometry.
    ///
    /// # Returns
    /// The total length of all edges in the wire.
    pub fn length_with_brep(&self, brep: &BRep) -> f64 {
        let mut total = 0.0;

        for &(edge_idx, _orientation) in &self.edges {
            if edge_idx >= brep.edges.len() {
                continue;
            }

            // Try to get curve bounds from geometry store
            if let Some(range) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) {
                if let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c) {
                    if let Some(curve) = brep.geom.curves.get(curve_idx) {
                        // Approximate length by sampling
                        let steps = 10;
                        let du = (range[1] - range[0]) / steps as f64;
                        for i in 0..steps {
                            let u = range[0] + du * i as f64;
                            let u_next = range[0] + du * (i + 1) as f64;
                            let p1 = curve.point_at(u);
                            let p2 = curve.point_at(u_next);
                            total += (p2 - p1).length();
                        }
                    }
                }
            }
        }

        total
    }

    /// Clear all edges from the wire.
    pub fn clear(&mut self) {
        self.edges.clear();
        self.cached_length = None;
        self.closed = None;
    }

    /// Check if the wire contains a specific edge.
    ///
    /// # Arguments
    /// * `edge_idx` - Edge index to search for.
    ///
    /// # Returns
    /// true if the edge is in the wire.
    pub fn contains_edge(&self, edge_idx: usize) -> bool {
        self.edges.iter().any(|(idx, _)| *idx == edge_idx)
    }

    /// Find the position of an edge in the wire.
    ///
    /// # Arguments
    /// * `edge_idx` - Edge index to search for.
    ///
    /// # Returns
    /// The position of the edge, or None if not found.
    pub fn find_edge(&self, edge_idx: usize) -> Option<usize> {
        self.edges.iter().position(|(idx, _)| *idx == edge_idx)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ShapeExtend_CompositeSurface - Composite Surface
// ─────────────────────────────────────────────────────────────────────────────

/// A patch in a composite surface.
#[derive(Debug, Clone)]
struct SurfacePatch {
    /// The surface geometry.
    surface: Surface3,
    /// U parameter range for this patch.
    u_range: [f64; 2],
    /// V parameter range for this patch.
    v_range: [f64; 2],
}

/// Composite surface made of multiple surface patches.
///
/// Analogous to OCCT `ShapeExtend_CompositeSurface`, this provides
/// a unified interface for a collection of surface patches with
/// defined parameter ranges.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::shape_extend::CompositeSurface;
/// use rcad_kernel::geom::{Surface3, Plane};
/// use glam::DVec3;
///
/// let mut composite = CompositeSurface::new();
/// let plane = Surface3::Plane(rcad_kernel::geom::Plane {
///     origin: DVec3::ZERO,
///     normal: DVec3::Z,
/// });
/// composite.add_surface(plane, [0.0, 1.0], [0.0, 1.0]);
/// assert_eq!(composite.patch_count(), 1);
/// ```
#[derive(Debug, Clone)]
pub struct CompositeSurface {
    /// List of surface patches.
    patches: Vec<SurfacePatch>,
    /// Global U range.
    u_range: [f64; 2],
    /// Global V range.
    v_range: [f64; 2],
}

impl Default for CompositeSurface {
    fn default() -> Self {
        Self::new()
    }
}

impl CompositeSurface {
    /// Create a new empty composite surface.
    pub fn new() -> Self {
        CompositeSurface {
            patches: Vec::new(),
            u_range: [0.0, 0.0],
            v_range: [0.0, 0.0],
        }
    }

    /// Create a composite surface with pre-allocated capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        CompositeSurface {
            patches: Vec::with_capacity(capacity),
            u_range: [0.0, 0.0],
            v_range: [0.0, 0.0],
        }
    }

    /// Add a surface patch to the composite.
    ///
    /// # Arguments
    /// * `surface` - The surface geometry.
    /// * `u_range` - U parameter range [u_min, u_max].
    /// * `v_range` - V parameter range [v_min, v_max].
    pub fn add_surface(&mut self, surface: Surface3, u_range: [f64; 2], v_range: [f64; 2]) {
        // Update global range
        if self.patches.is_empty() {
            self.u_range = u_range;
            self.v_range = v_range;
        } else {
            self.u_range[0] = self.u_range[0].min(u_range[0]);
            self.u_range[1] = self.u_range[1].max(u_range[1]);
            self.v_range[0] = self.v_range[0].min(v_range[0]);
            self.v_range[1] = self.v_range[1].max(v_range[1]);
        }

        self.patches.push(SurfacePatch {
            surface,
            u_range,
            v_range,
        });
    }

    /// Get the number of surface patches.
    pub fn patch_count(&self) -> usize {
        self.patches.len()
    }

    /// Check if the composite has no patches.
    pub fn is_empty(&self) -> bool {
        self.patches.is_empty()
    }

    /// Get the global U parameter range.
    pub fn global_u_range(&self) -> [f64; 2] {
        self.u_range
    }

    /// Get the global V parameter range.
    pub fn global_v_range(&self) -> [f64; 2] {
        self.v_range
    }

    /// Get the surface patch at the given global parameters.
    ///
    /// # Arguments
    /// * `u` - Global U parameter.
    /// * `v` - Global V parameter.
    ///
    /// # Returns
    /// Reference to the surface containing the given parameters, or None.
    pub fn surface_at(&self, u: f64, v: f64) -> Option<&Surface3> {
        for patch in &self.patches {
            if u >= patch.u_range[0] && u <= patch.u_range[1]
                && v >= patch.v_range[0] && v <= patch.v_range[1]
            {
                return Some(&patch.surface);
            }
        }
        None
    }

    /// Convert global parameters to local patch parameters.
    ///
    /// # Arguments
    /// * `u` - Global U parameter.
    /// * `v` - Global V parameter.
    ///
    /// # Returns
    /// A tuple (patch_index, local_u, local_v), or (0, u, v) if no patch found.
    pub fn local_params(&self, u: f64, v: f64) -> (usize, f64, f64) {
        for (idx, patch) in self.patches.iter().enumerate() {
            if u >= patch.u_range[0] && u <= patch.u_range[1]
                && v >= patch.v_range[0] && v <= patch.v_range[1]
            {
                // Map global params to local params (identity for now)
                // In a real implementation, this would apply a transformation
                let local_u = (u - patch.u_range[0]) / (patch.u_range[1] - patch.u_range[0]);
                let local_v = (v - patch.v_range[0]) / (patch.v_range[1] - patch.v_range[0]);
                return (idx, local_u, local_v);
            }
        }
        (0, u, v)
    }

    /// Get a patch by index.
    ///
    /// # Arguments
    /// * `idx` - Patch index.
    ///
    /// # Returns
    /// The patch surface and ranges, or None if out of bounds.
    pub fn patch(&self, idx: usize) -> Option<(&Surface3, [f64; 2], [f64; 2])> {
        self.patches.get(idx).map(|p| (&p.surface, p.u_range, p.v_range))
    }

    /// Evaluate the composite surface at the given parameters.
    ///
    /// # Arguments
    /// * `u` - Global U parameter.
    /// * `v` - Global V parameter.
    ///
    /// # Returns
    /// The point on the surface, or None if parameters are outside all patches.
    pub fn point_at(&self, u: f64, v: f64) -> Option<glam::DVec3> {
        self.surface_at(u, v).map(|surf| surf.point_at(u, v))
    }

    /// Evaluate the normal at the given parameters.
    ///
    /// # Arguments
    /// * `u` - Global U parameter.
    /// * `v` - Global V parameter.
    ///
    /// # Returns
    /// The normal vector, or None if parameters are outside all patches.
    pub fn normal_at(&self, u: f64, v: f64) -> Option<glam::DVec3> {
        self.surface_at(u, v).map(|surf| surf.normal_at(u, v))
    }

    /// Clear all patches from the composite surface.
    pub fn clear(&mut self) {
        self.patches.clear();
        self.u_range = [0.0, 0.0];
        self.v_range = [0.0, 0.0];
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ShapeExtend_BasicMsgRegistrator - Basic Message Registration
// ─────────────────────────────────────────────────────────────────────────────

/// Message severity level.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MessageSeverity {
    /// Informational message.
    Info,
    /// Warning message.
    Warning,
    /// Error message.
    Error,
    /// Failure message (operation failed).
    Fail,
}

impl std::fmt::Display for MessageSeverity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MessageSeverity::Info => write!(f, "Info"),
            MessageSeverity::Warning => write!(f, "Warning"),
            MessageSeverity::Error => write!(f, "Error"),
            MessageSeverity::Fail => write!(f, "Fail"),
        }
    }
}

/// A registered message with severity.
#[derive(Debug, Clone)]
pub struct ShapeMessage {
    /// The message text.
    pub message: String,
    /// The message severity.
    pub severity: MessageSeverity,
}

impl ShapeMessage {
    /// Create a new shape message.
    pub fn new(message: impl Into<String>, severity: MessageSeverity) -> Self {
        ShapeMessage {
            message: message.into(),
            severity,
        }
    }

    /// Create an informational message.
    pub fn info(message: impl Into<String>) -> Self {
        Self::new(message, MessageSeverity::Info)
    }

    /// Create a warning message.
    pub fn warning(message: impl Into<String>) -> Self {
        Self::new(message, MessageSeverity::Warning)
    }

    /// Create an error message.
    pub fn error(message: impl Into<String>) -> Self {
        Self::new(message, MessageSeverity::Error)
    }

    /// Create a failure message.
    pub fn fail(message: impl Into<String>) -> Self {
        Self::new(message, MessageSeverity::Fail)
    }
}

/// Basic message registrator for collecting messages.
///
/// Analogous to OCCT `ShapeExtend_BasicMsgRegistrator`, this provides
/// a simple message collection mechanism.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::shape_extend::{MessageRegistrator, MessageSeverity};
///
/// let mut reg = MessageRegistrator::new();
/// reg.add_message("Processing started", MessageSeverity::Info);
/// reg.add_message("Found small edge", MessageSeverity::Warning);
/// assert_eq!(reg.message_count(), 2);
/// ```
#[derive(Debug, Clone, Default)]
pub struct MessageRegistrator {
    /// List of registered messages.
    messages: Vec<ShapeMessage>,
}

impl MessageRegistrator {
    /// Create a new empty message registrator.
    pub fn new() -> Self {
        MessageRegistrator {
            messages: Vec::new(),
        }
    }

    /// Create a registrator with pre-allocated capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        MessageRegistrator {
            messages: Vec::with_capacity(capacity),
        }
    }

    /// Add a message with the given severity.
    ///
    /// # Arguments
    /// * `msg` - The message text.
    /// * `severity` - The message severity.
    pub fn add_message(&mut self, msg: impl Into<String>, severity: MessageSeverity) {
        self.messages.push(ShapeMessage::new(msg, severity));
    }

    /// Add an informational message.
    pub fn add_info(&mut self, msg: impl Into<String>) {
        self.add_message(msg, MessageSeverity::Info);
    }

    /// Add a warning message.
    pub fn add_warning(&mut self, msg: impl Into<String>) {
        self.add_message(msg, MessageSeverity::Warning);
    }

    /// Add an error message.
    pub fn add_error(&mut self, msg: impl Into<String>) {
        self.add_message(msg, MessageSeverity::Error);
    }

    /// Add a failure message.
    pub fn add_fail(&mut self, msg: impl Into<String>) {
        self.add_message(msg, MessageSeverity::Fail);
    }

    /// Get all registered messages.
    pub fn messages(&self) -> &[ShapeMessage] {
        &self.messages
    }

    /// Get the number of messages.
    pub fn message_count(&self) -> usize {
        self.messages.len()
    }

    /// Check if there are no messages.
    pub fn is_empty(&self) -> bool {
        self.messages.is_empty()
    }

    /// Get messages of a specific severity.
    ///
    /// # Arguments
    /// * `severity` - The severity to filter by.
    ///
    /// # Returns
    /// Vector of messages with the given severity.
    pub fn messages_by_severity(&self, severity: MessageSeverity) -> Vec<&ShapeMessage> {
        self.messages.iter().filter(|m| m.severity == severity).collect()
    }

    /// Count messages of a specific severity.
    ///
    /// # Arguments
    /// * `severity` - The severity to count.
    ///
    /// # Returns
    /// Number of messages with the given severity.
    pub fn count_by_severity(&self, severity: MessageSeverity) -> usize {
        self.messages.iter().filter(|m| m.severity == severity).count()
    }

    /// Check if there are any errors or failures.
    pub fn has_errors(&self) -> bool {
        self.messages.iter().any(|m| {
            m.severity == MessageSeverity::Error || m.severity == MessageSeverity::Fail
        })
    }

    /// Check if there are any warnings.
    pub fn has_warnings(&self) -> bool {
        self.messages.iter().any(|m| m.severity == MessageSeverity::Warning)
    }

    /// Clear all messages.
    pub fn clear(&mut self) {
        self.messages.clear();
    }

    /// Merge messages from another registrator.
    pub fn merge(&mut self, other: &MessageRegistrator) {
        self.messages.extend(other.messages.clone());
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ShapeExtend_MsgRegistrator - Message Registration with Shape Context
// ─────────────────────────────────────────────────────────────────────────────

/// A message associated with a shape index.
#[derive(Debug, Clone)]
pub struct ShapeContextMessage {
    /// The shape index the message is associated with.
    pub shape_idx: usize,
    /// The message.
    pub message: ShapeMessage,
}

/// Message registrator with shape context.
///
/// Analogous to OCCT `ShapeExtend_MsgRegistrator`, this extends
/// the basic registrator to associate messages with specific shapes.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::shape_extend::{ShapeMessageRegistrator, MessageSeverity};
///
/// let mut reg = ShapeMessageRegistrator::new();
/// reg.add_shape_message(0, "Edge is too short", MessageSeverity::Warning);
/// reg.add_shape_message(1, "Face has bad normal", MessageSeverity::Error);
/// assert_eq!(reg.message_count(), 2);
/// ```
#[derive(Debug, Clone, Default)]
pub struct ShapeMessageRegistrator {
    /// List of messages with shape context.
    messages: Vec<ShapeContextMessage>,
}

impl ShapeMessageRegistrator {
    /// Create a new empty shape message registrator.
    pub fn new() -> Self {
        ShapeMessageRegistrator {
            messages: Vec::new(),
        }
    }

    /// Create a registrator with pre-allocated capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        ShapeMessageRegistrator {
            messages: Vec::with_capacity(capacity),
        }
    }

    /// Add a message associated with a shape.
    ///
    /// # Arguments
    /// * `shape_idx` - Index of the shape the message is about.
    /// * `msg` - The message text.
    /// * `severity` - The message severity.
    pub fn add_shape_message(&mut self, shape_idx: usize, msg: impl Into<String>, severity: MessageSeverity) {
        self.messages.push(ShapeContextMessage {
            shape_idx,
            message: ShapeMessage::new(msg, severity),
        });
    }

    /// Add an informational message for a shape.
    pub fn add_info(&mut self, shape_idx: usize, msg: impl Into<String>) {
        self.add_shape_message(shape_idx, msg, MessageSeverity::Info);
    }

    /// Add a warning message for a shape.
    pub fn add_warning(&mut self, shape_idx: usize, msg: impl Into<String>) {
        self.add_shape_message(shape_idx, msg, MessageSeverity::Warning);
    }

    /// Add an error message for a shape.
    pub fn add_error(&mut self, shape_idx: usize, msg: impl Into<String>) {
        self.add_shape_message(shape_idx, msg, MessageSeverity::Error);
    }

    /// Add a failure message for a shape.
    pub fn add_fail(&mut self, shape_idx: usize, msg: impl Into<String>) {
        self.add_shape_message(shape_idx, msg, MessageSeverity::Fail);
    }

    /// Get all messages with shape context.
    pub fn messages(&self) -> &[ShapeContextMessage] {
        &self.messages
    }

    /// Get the number of messages.
    pub fn message_count(&self) -> usize {
        self.messages.len()
    }

    /// Check if there are no messages.
    pub fn is_empty(&self) -> bool {
        self.messages.is_empty()
    }

    /// Get all messages for a specific shape.
    ///
    /// # Arguments
    /// * `shape_idx` - The shape index to filter by.
    ///
    /// # Returns
    /// Vector of messages for the given shape.
    pub fn messages_for_shape(&self, shape_idx: usize) -> Vec<&ShapeContextMessage> {
        self.messages.iter().filter(|m| m.shape_idx == shape_idx).collect()
    }

    /// Get all shapes that have messages of a specific severity.
    ///
    /// # Arguments
    /// * `severity` - The severity to filter by.
    ///
    /// # Returns
    /// Vector of shape indices with messages of the given severity.
    pub fn shapes_with_severity(&self, severity: MessageSeverity) -> Vec<usize> {
        let mut shapes: Vec<usize> = self.messages
            .iter()
            .filter(|m| m.message.severity == severity)
            .map(|m| m.shape_idx)
            .collect();
        shapes.sort();
        shapes.dedup();
        shapes
    }

    /// Check if there are any errors or failures.
    pub fn has_errors(&self) -> bool {
        self.messages.iter().any(|m| {
            m.message.severity == MessageSeverity::Error || m.message.severity == MessageSeverity::Fail
        })
    }

    /// Check if there are any warnings.
    pub fn has_warnings(&self) -> bool {
        self.messages.iter().any(|m| m.message.severity == MessageSeverity::Warning)
    }

    /// Clear all messages.
    pub fn clear(&mut self) {
        self.messages.clear();
    }

    /// Merge messages from another registrator.
    pub fn merge(&mut self, other: &ShapeMessageRegistrator) {
        self.messages.extend(other.messages.clone());
    }

    /// Convert to a basic message registrator (without shape context).
    pub fn to_basic(&self) -> MessageRegistrator {
        let mut basic = MessageRegistrator::with_capacity(self.messages.len());
        for msg in &self.messages {
            basic.add_message(msg.message.message.clone(), msg.message.severity);
        }
        basic
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ShapeExtend_Explorer - Extended Shape Exploration
// ─────────────────────────────────────────────────────────────────────────────

/// Extended shape exploration utilities.
///
/// Analogous to OCCT `ShapeExtend_Explorer`, this provides
/// enhanced shape traversal and query capabilities.
pub struct ShapeExplorer;

impl ShapeExplorer {
    /// Helper to iterate over all faces in a BRep.
    fn iter_faces(brep: &BRep) -> impl Iterator<Item = (usize, &Face)> {
        brep.solids.iter()
            .flat_map(|solid| solid.shells.iter())
            .flat_map(|shell| shell.faces.iter())
            .enumerate()
    }

    /// Helper to iterate over all wires in a face (outer and inner).
    fn iter_wires(face: &Face) -> impl Iterator<Item = &Wire> {
        std::iter::once(&face.outer_wire).chain(face.inner_wires.iter())
    }

    /// Explore all subshapes of a given shape.
    ///
    /// # Arguments
    /// * `brep` - The BRep to explore.
    /// * `shape_idx` - The shape index (solid index in this simplified implementation).
    ///
    /// # Returns
    /// Vector of shell indices for the given solid.
    ///
    /// # Note
    /// In this simplified implementation, shape_idx is interpreted as:
    /// - For a solid: returns the shell indices (0-based, consecutive)
    pub fn explore_shape(brep: &BRep, shape_idx: usize) -> Vec<usize> {
        if shape_idx < brep.solids.len() {
            // Return shell count for this solid
            (0..brep.solids[shape_idx].shells.len()).collect()
        } else {
            Vec::new()
        }
    }

    /// Count subshapes of a given type in a BRep.
    ///
    /// # Arguments
    /// * `brep` - The BRep to analyze.
    /// * `shape_type` - The type of shape to count.
    ///
    /// # Returns
    /// Number of shapes of the given type.
    pub fn count_subshapes(brep: &BRep, shape_type: ShapeType) -> usize {
        match shape_type {
            ShapeType::Vertex => brep.vertices.len(),
            ShapeType::Edge => brep.edges.len(),
            ShapeType::Wire => {
                // Count all wires (outer and inner) across all faces
                brep.solids.iter()
                    .flat_map(|solid| solid.shells.iter())
                    .flat_map(|shell| shell.faces.iter())
                    .map(|face| 1 + face.inner_wires.len())
                    .sum()
            }
            ShapeType::Face => {
                brep.solids.iter()
                    .flat_map(|solid| solid.shells.iter())
                    .map(|shell| shell.faces.len())
                    .sum()
            }
            ShapeType::Shell => {
                brep.solids.iter()
                    .map(|solid| solid.shells.len())
                    .sum()
            }
            ShapeType::Solid => brep.solids.len(),
            ShapeType::CompSolid => brep.compsolid.as_ref().map_or(1, |_| 1),
            ShapeType::Compound => brep.compound.as_ref().map_or(1, |_| 1),
            ShapeType::Empty => 0,
        }
    }

    /// Get all unique edge indices in the BRep.
    ///
    /// # Arguments
    /// * `brep` - The BRep to explore.
    ///
    /// # Returns
    /// Vector of all edge indices referenced in wires.
    pub fn all_edges(brep: &BRep) -> Vec<usize> {
        let mut edges: Vec<usize> = Vec::new();

        // Collect edges from all faces through the topology hierarchy
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for wire in Self::iter_wires(face) {
                        for wire_edge in &wire.edges {
                            edges.push(wire_edge.idx);
                        }
                    }
                }
            }
        }

        // Sort and deduplicate
        edges.sort();
        edges.dedup();
        edges
    }

    /// Get all unique vertex indices in the BRep.
    ///
    /// # Arguments
    /// * `brep` - The BRep to explore.
    ///
    /// # Returns
    /// Vector of all vertex indices used by edges.
    pub fn all_vertices(brep: &BRep) -> Vec<usize> {
        let mut vertices: Vec<usize> = Vec::new();

        for edge in &brep.edges {
            vertices.push(edge.start);
            vertices.push(edge.end);
        }

        vertices.sort();
        vertices.dedup();
        vertices
    }

    /// Find faces that share a given edge.
    ///
    /// # Arguments
    /// * `brep` - The BRep to search.
    /// * `edge_idx` - The edge index to find sharing faces for.
    ///
    /// # Returns
    /// Vector of face indices (0-based global face index) that reference the given edge.
    pub fn faces_sharing_edge(brep: &BRep, edge_idx: usize) -> Vec<usize> {
        let mut faces = Vec::new();
        let mut face_counter = 0;

        for solid in &brep.solids {
            for shell in &solid.shells {
                for _face in &shell.faces {
                    // Check if this face contains the edge
                    let has_edge = Self::iter_wires(_face)
                        .any(|wire| wire.edges.iter().any(|we| we.idx == edge_idx));
                    if has_edge {
                        faces.push(face_counter);
                    }
                    face_counter += 1;
                }
            }
        }

        faces
    }

    /// Find edges that share a given vertex.
    ///
    /// # Arguments
    /// * `brep` - The BRep to search.
    /// * `vertex_idx` - The vertex index to find sharing edges for.
    ///
    /// # Returns
    /// Vector of edge indices that reference the given vertex.
    pub fn edges_sharing_vertex(brep: &BRep, vertex_idx: usize) -> Vec<usize> {
        let mut edges = Vec::new();

        for (edge_idx, edge) in brep.edges.iter().enumerate() {
            if edge.start == vertex_idx || edge.end == vertex_idx {
                edges.push(edge_idx);
            }
        }

        edges
    }

    /// Find the boundary edges of a shell or solid.
    ///
    /// # Arguments
    /// * `brep` - The BRep to analyze.
    ///
    /// # Returns
    /// Vector of edge indices that appear exactly once (boundary edges).
    pub fn boundary_edges(brep: &BRep) -> Vec<usize> {
        let mut edge_counts: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();

        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for wire in Self::iter_wires(face) {
                        for wire_edge in &wire.edges {
                            *edge_counts.entry(wire_edge.idx).or_insert(0) += 1;
                        }
                    }
                }
            }
        }

        edge_counts
            .into_iter()
            .filter(|(_, count)| *count == 1)
            .map(|(idx, _)| idx)
            .collect()
    }

    /// Find non-manifold edges in the BRep.
    ///
    /// A non-manifold edge is shared by more than two faces.
    ///
    /// # Arguments
    /// * `brep` - The BRep to analyze.
    ///
    /// # Returns
    /// Vector of non-manifold edge indices.
    pub fn non_manifold_edges(brep: &BRep) -> Vec<usize> {
        let mut edge_counts: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();

        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for wire in Self::iter_wires(face) {
                        for wire_edge in &wire.edges {
                            *edge_counts.entry(wire_edge.idx).or_insert(0) += 1;
                        }
                    }
                }
            }
        }

        edge_counts
            .into_iter()
            .filter(|(_, count)| *count > 2)
            .map(|(idx, _)| idx)
            .collect()
    }

    /// Get a summary of the BRep topology.
    ///
    /// # Arguments
    /// * `brep` - The BRep to summarize.
    ///
    /// # Returns
    /// A string describing the topology.
    pub fn topology_summary(brep: &BRep) -> String {
        let face_count = Self::count_subshapes(brep, ShapeType::Face);
        let shell_count = Self::count_subshapes(brep, ShapeType::Shell);
        let wire_count = Self::count_subshapes(brep, ShapeType::Wire);

        format!(
            "BRep topology: {} vertices, {} edges, {} wires, {} faces, {} shells, {} solids",
            brep.vertices.len(),
            brep.edges.len(),
            wire_count,
            face_count,
            shell_count,
            brep.solids.len()
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{Plane, Surface3};
    use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

    // Helper to create a simple test BRep with a triangular face
    fn create_triangular_brep() -> BRep {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.vertices.push(Vertex { point: DVec3::Y });

        // Closed triangle edges
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let wire = Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
        };

        let face = Face {
            outer_wire: wire,
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };

        let shell = Shell { faces: vec![face] };
        let solid = Solid { shells: vec![shell] };
        brep.solids.push(solid);

        brep
    }

    // ───────────────────────────────────────────────────────────────────────
    // WireData Tests
    // ───────────────────────────────────────────────────────────────────────

    #[test]
    fn test_wire_data_new() {
        let wire = WireData::new();
        assert!(wire.is_empty());
        assert_eq!(wire.edge_count(), 0);
    }

    #[test]
    fn test_wire_data_add_edge() {
        let mut wire = WireData::new();
        wire.add_edge(0, true);
        wire.add_edge(1, false);
        wire.add_edge(2, true);

        assert_eq!(wire.edge_count(), 3);
        assert_eq!(wire.edges(), &[(0, true), (1, false), (2, true)]);
    }

    #[test]
    fn test_wire_data_add_edge_at() {
        let mut wire = WireData::new();
        wire.add_edge(0, true);
        wire.add_edge(2, true);
        wire.add_edge_at(1, 1, false);

        assert_eq!(wire.edges(), &[(0, true), (1, false), (2, true)]);
    }

    #[test]
    fn test_wire_data_remove_edge() {
        let mut wire = WireData::new();
        wire.add_edge(0, true);
        wire.add_edge(1, false);
        wire.add_edge(2, true);

        let removed = wire.remove_edge(1);
        assert_eq!(removed, (1, false));
        assert_eq!(wire.edges(), &[(0, true), (2, true)]);
    }

    #[test]
    fn test_wire_data_edge_at() {
        let mut wire = WireData::new();
        wire.add_edge(0, true);
        wire.add_edge(1, false);

        assert_eq!(wire.edge_at(0), Some((0, true)));
        assert_eq!(wire.edge_at(1), Some((1, false)));
        assert_eq!(wire.edge_at(2), None);
    }

    #[test]
    fn test_wire_data_contains_edge() {
        let mut wire = WireData::new();
        wire.add_edge(0, true);
        wire.add_edge(2, false);

        assert!(wire.contains_edge(0));
        assert!(!wire.contains_edge(1));
        assert!(wire.contains_edge(2));
    }

    #[test]
    fn test_wire_data_find_edge() {
        let mut wire = WireData::new();
        wire.add_edge(0, true);
        wire.add_edge(1, false);
        wire.add_edge(2, true);

        assert_eq!(wire.find_edge(0), Some(0));
        assert_eq!(wire.find_edge(1), Some(1));
        assert_eq!(wire.find_edge(3), None);
    }

    #[test]
    fn test_wire_data_reverse() {
        let mut wire = WireData::new();
        wire.add_edge(0, true);
        wire.add_edge(1, false);
        wire.add_edge(2, true);
        wire.reverse();

        assert_eq!(wire.edges(), &[(2, false), (1, true), (0, false)]);
    }

    #[test]
    fn test_wire_data_set_orientation() {
        let mut wire = WireData::new();
        wire.add_edge(0, true);
        wire.add_edge(1, false);
        wire.set_orientation(0, false);

        assert_eq!(wire.edges(), &[(0, false), (1, false)]);
    }

    #[test]
    fn test_wire_data_clear() {
        let mut wire = WireData::new();
        wire.add_edge(0, true);
        wire.add_edge(1, false);
        wire.clear();

        assert!(wire.is_empty());
        assert_eq!(wire.edge_count(), 0);
    }

    #[test]
    fn test_wire_data_is_closed_with_brep() {
        let brep = create_triangular_brep();

        let mut wire = WireData::new();
        wire.add_edge(0, true);
        wire.add_edge(1, true);
        wire.add_edge(2, true);

        assert!(wire.is_closed_with_brep(&brep));

        // Open wire
        let mut open_wire = WireData::new();
        open_wire.add_edge(0, true);
        open_wire.add_edge(1, true);

        assert!(!open_wire.is_closed_with_brep(&brep));
    }

    // ───────────────────────────────────────────────────────────────────────
    // CompositeSurface Tests
    // ───────────────────────────────────────────────────────────────────────

    #[test]
    fn test_composite_surface_new() {
        let composite = CompositeSurface::new();
        assert!(composite.is_empty());
        assert_eq!(composite.patch_count(), 0);
    }

    #[test]
    fn test_composite_surface_add_surface() {
        let mut composite = CompositeSurface::new();
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        composite.add_surface(plane.clone(), [0.0, 1.0], [0.0, 1.0]);
        assert_eq!(composite.patch_count(), 1);
        assert_eq!(composite.global_u_range(), [0.0, 1.0]);
        assert_eq!(composite.global_v_range(), [0.0, 1.0]);

        composite.add_surface(plane, [1.0, 2.0], [0.0, 1.0]);
        assert_eq!(composite.patch_count(), 2);
        assert_eq!(composite.global_u_range(), [0.0, 2.0]);
    }

    #[test]
    fn test_composite_surface_at() {
        let mut composite = CompositeSurface::new();
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        composite.add_surface(plane, [0.0, 1.0], [0.0, 1.0]);

        // Inside patch
        assert!(composite.surface_at(0.5, 0.5).is_some());

        // Outside patch
        assert!(composite.surface_at(1.5, 0.5).is_none());
    }

    #[test]
    fn test_composite_surface_local_params() {
        let mut composite = CompositeSurface::new();
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        composite.add_surface(plane, [0.0, 1.0], [0.0, 1.0]);

        let (idx, u, v) = composite.local_params(0.5, 0.5);
        assert_eq!(idx, 0);
        assert!((u - 0.5).abs() < 1e-10);
        assert!((v - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_composite_surface_clear() {
        let mut composite = CompositeSurface::new();
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        composite.add_surface(plane, [0.0, 1.0], [0.0, 1.0]);
        composite.clear();

        assert!(composite.is_empty());
        assert_eq!(composite.patch_count(), 0);
    }

    // ───────────────────────────────────────────────────────────────────────
    // MessageRegistrator Tests
    // ───────────────────────────────────────────────────────────────────────

    #[test]
    fn test_message_registrator_new() {
        let reg = MessageRegistrator::new();
        assert!(reg.is_empty());
        assert_eq!(reg.message_count(), 0);
    }

    #[test]
    fn test_message_registrator_add_message() {
        let mut reg = MessageRegistrator::new();
        reg.add_message("Test message", MessageSeverity::Info);
        reg.add_message("Warning message", MessageSeverity::Warning);
        reg.add_message("Error message", MessageSeverity::Error);

        assert_eq!(reg.message_count(), 3);
    }

    #[test]
    fn test_message_registrator_convenience_methods() {
        let mut reg = MessageRegistrator::new();
        reg.add_info("Info");
        reg.add_warning("Warning");
        reg.add_error("Error");
        reg.add_fail("Fail");

        assert_eq!(reg.message_count(), 4);
        assert!(reg.has_errors());
        assert!(reg.has_warnings());
    }

    #[test]
    fn test_message_registrator_filter_by_severity() {
        let mut reg = MessageRegistrator::new();
        reg.add_info("Info 1");
        reg.add_warning("Warning 1");
        reg.add_error("Error 1");
        reg.add_warning("Warning 2");

        assert_eq!(reg.count_by_severity(MessageSeverity::Info), 1);
        assert_eq!(reg.count_by_severity(MessageSeverity::Warning), 2);
        assert_eq!(reg.count_by_severity(MessageSeverity::Error), 1);
    }

    #[test]
    fn test_message_registrator_clear() {
        let mut reg = MessageRegistrator::new();
        reg.add_message("Test", MessageSeverity::Info);
        reg.clear();

        assert!(reg.is_empty());
    }

    #[test]
    fn test_message_registrator_merge() {
        let mut reg1 = MessageRegistrator::new();
        reg1.add_info("Message 1");

        let mut reg2 = MessageRegistrator::new();
        reg2.add_warning("Message 2");

        reg1.merge(&reg2);
        assert_eq!(reg1.message_count(), 2);
    }

    // ───────────────────────────────────────────────────────────────────────
    // ShapeMessageRegistrator Tests
    // ───────────────────────────────────────────────────────────────────────

    #[test]
    fn test_shape_message_registrator_new() {
        let reg = ShapeMessageRegistrator::new();
        assert!(reg.is_empty());
        assert_eq!(reg.message_count(), 0);
    }

    #[test]
    fn test_shape_message_registrator_add_shape_message() {
        let mut reg = ShapeMessageRegistrator::new();
        reg.add_shape_message(0, "Edge too short", MessageSeverity::Warning);
        reg.add_shape_message(1, "Invalid face", MessageSeverity::Error);

        assert_eq!(reg.message_count(), 2);
    }

    #[test]
    fn test_shape_message_registrator_messages_for_shape() {
        let mut reg = ShapeMessageRegistrator::new();
        reg.add_shape_message(0, "Message 1", MessageSeverity::Info);
        reg.add_shape_message(1, "Message 2", MessageSeverity::Warning);
        reg.add_shape_message(0, "Message 3", MessageSeverity::Error);

        let msgs = reg.messages_for_shape(0);
        assert_eq!(msgs.len(), 2);
    }

    #[test]
    fn test_shape_message_registrator_shapes_with_severity() {
        let mut reg = ShapeMessageRegistrator::new();
        reg.add_shape_message(0, "Info", MessageSeverity::Info);
        reg.add_shape_message(1, "Warning", MessageSeverity::Warning);
        reg.add_shape_message(2, "Error", MessageSeverity::Error);
        reg.add_shape_message(3, "Warning 2", MessageSeverity::Warning);

        let warning_shapes = reg.shapes_with_severity(MessageSeverity::Warning);
        assert_eq!(warning_shapes, vec![1, 3]);
    }

    #[test]
    fn test_shape_message_registrator_to_basic() {
        let mut reg = ShapeMessageRegistrator::new();
        reg.add_shape_message(0, "Message 1", MessageSeverity::Info);
        reg.add_shape_message(1, "Message 2", MessageSeverity::Warning);

        let basic = reg.to_basic();
        assert_eq!(basic.message_count(), 2);
    }

    // ───────────────────────────────────────────────────────────────────────
    // ShapeExplorer Tests
    // ───────────────────────────────────────────────────────────────────────

    #[test]
    fn test_shape_explorer_count_subshapes() {
        let brep = create_triangular_brep();

        assert_eq!(ShapeExplorer::count_subshapes(&brep, ShapeType::Vertex), 3);
        assert_eq!(ShapeExplorer::count_subshapes(&brep, ShapeType::Edge), 3);
        assert_eq!(ShapeExplorer::count_subshapes(&brep, ShapeType::Face), 1);
        assert_eq!(ShapeExplorer::count_subshapes(&brep, ShapeType::Shell), 1);
        assert_eq!(ShapeExplorer::count_subshapes(&brep, ShapeType::Solid), 1);
    }

    #[test]
    fn test_shape_explorer_all_edges() {
        let brep = create_triangular_brep();
        let edges = ShapeExplorer::all_edges(&brep);
        assert_eq!(edges, vec![0, 1, 2]);
    }

    #[test]
    fn test_shape_explorer_all_vertices() {
        let brep = create_triangular_brep();
        let vertices = ShapeExplorer::all_vertices(&brep);
        assert_eq!(vertices, vec![0, 1, 2]);
    }

    #[test]
    fn test_shape_explorer_edges_sharing_vertex() {
        let brep = create_triangular_brep();

        let edges = ShapeExplorer::edges_sharing_vertex(&brep, 0);
        assert_eq!(edges, vec![0, 2]); // edges 0 (0->1) and 2 (2->0)

        let edges = ShapeExplorer::edges_sharing_vertex(&brep, 1);
        assert_eq!(edges, vec![0, 1]); // edges 0 (0->1) and 1 (1->2)
    }

    #[test]
    fn test_shape_explorer_boundary_edges() {
        let brep = create_triangular_brep();
        let boundary = ShapeExplorer::boundary_edges(&brep);
        // All edges appear exactly once, so all are boundary edges
        assert_eq!(boundary.len(), 3);
    }

    #[test]
    fn test_shape_explorer_non_manifold_edges() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(Edge { start: 0, end: 1 });

        // Create multiple faces sharing the same edge
        let wire = Wire {
            edges: vec![WireEdge::fwd(0)],
        };

        // 3 faces sharing the same edge -> non-manifold
        let mut faces = Vec::new();
        for _ in 0..3 {
            faces.push(Face {
                outer_wire: wire.clone(),
                inner_wires: vec![],
                normal: DVec3::Z,
                triangles: vec![],
                mesh_dirty: true,
            });
        }

        let shell = Shell { faces };
        let solid = Solid { shells: vec![shell] };
        brep.solids.push(solid);

        let nm_edges = ShapeExplorer::non_manifold_edges(&brep);
        assert_eq!(nm_edges, vec![0]);
    }

    #[test]
    fn test_shape_explorer_topology_summary() {
        let brep = create_triangular_brep();
        let summary = ShapeExplorer::topology_summary(&brep);
        assert!(summary.contains("3 vertices"));
        assert!(summary.contains("3 edges"));
        assert!(summary.contains("1 faces"));
    }

    #[test]
    fn test_shape_explorer_explore_shape() {
        let brep = create_triangular_brep();
        // Explore the solid (index 0) - should return shell indices
        let subshapes = ShapeExplorer::explore_shape(&brep, 0);
        assert_eq!(subshapes, vec![0]); // 1 shell at index 0
    }

    // ───────────────────────────────────────────────────────────────────────
    // ShapeMessage Tests
    // ───────────────────────────────────────────────────────────────────────

    #[test]
    fn test_shape_message_convenience_constructors() {
        let info = ShapeMessage::info("Info message");
        assert_eq!(info.severity, MessageSeverity::Info);

        let warning = ShapeMessage::warning("Warning message");
        assert_eq!(warning.severity, MessageSeverity::Warning);

        let error = ShapeMessage::error("Error message");
        assert_eq!(error.severity, MessageSeverity::Error);

        let fail = ShapeMessage::fail("Fail message");
        assert_eq!(fail.severity, MessageSeverity::Fail);
    }

    #[test]
    fn test_message_severity_display() {
        assert_eq!(format!("{}", MessageSeverity::Info), "Info");
        assert_eq!(format!("{}", MessageSeverity::Warning), "Warning");
        assert_eq!(format!("{}", MessageSeverity::Error), "Error");
        assert_eq!(format!("{}", MessageSeverity::Fail), "Fail");
    }
}
