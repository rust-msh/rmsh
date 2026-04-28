//! BRepTopAdaptor-style topology adapters for high-level topology traversal.
//!
//! This module provides high-level adapters for exploring BRep topology,
//! analogous to OCCT's `TopExp_Explorer` and `BRepAdaptor` classes.
//!
//! # Overview
//!
//! - **Explorers**: `FaceExplorer`, `EdgeExplorer`, `VertexExplorer`, `WireExplorer`
//!   provide forward-only iteration over topology elements.
//! - **ShapeIterator**: Generic iterator implementing `std::iter::Iterator` for all shape types.
//! - **Topology queries**: Helper functions for adjacency queries.
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::brep_top_adaptor::*;
//! use rcad_kernel::BRep;
//!
//! let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
//!     width: 1.0, height: 1.0, depth: 1.0
//! });
//!
//! // Count faces using FaceExplorer
//! let mut explorer = FaceExplorer::new(&brep);
//! let mut face_count = 0;
//! while explorer.next().is_some() {
//!     face_count += 1;
//! }
//! assert_eq!(face_count, 6);
//! ```

use rcad_kernel::BRep;
use rcad_kernel::topology::{Face, WireEdge};
pub use crate::brep_tools::ShapeType;

// =============================================================================
// Face Adaptor
// =============================================================================

/// Adapter providing convenient access to face data.
///
/// Analogous to OCCT `BRepAdaptor_Face`.
#[derive(Debug, Clone)]
pub struct FaceAdaptor<'a> {
    brep: &'a BRep,
    face_idx: usize,
    shell_idx: usize,
    solid_idx: usize,
}

impl<'a> FaceAdaptor<'a> {
    /// Create a new face adaptor.
    pub fn new(brep: &'a BRep, face_idx: usize, shell_idx: usize, solid_idx: usize) -> Self {
        Self {
            brep,
            face_idx,
            shell_idx,
            solid_idx,
        }
    }

    /// Returns the face index in the flattened face list.
    pub fn index(&self) -> usize {
        self.face_idx
    }

    /// Returns the shell index containing this face.
    pub fn shell_index(&self) -> usize {
        self.shell_idx
    }

    /// Returns the solid index containing this face.
    pub fn solid_index(&self) -> usize {
        self.solid_idx
    }

    /// Returns a reference to the face topology.
    pub fn face(&self) -> Option<&'a Face> {
        self.brep.solids.get(self.solid_idx)
            .and_then(|s| s.shells.get(self.shell_idx))
            .and_then(|sh| sh.faces.get(self.face_idx))
    }

    /// Returns the surface index for this face, if available.
    pub fn surface_index(&self) -> Option<usize> {
        // Compute the flattened face index
        let mut flat_idx = 0usize;
        for (si, solid) in self.brep.solids.iter().enumerate() {
            for (shi, shell) in solid.shells.iter().enumerate() {
                for fi in 0..shell.faces.len() {
                    if si == self.solid_idx && shi == self.shell_idx && fi == self.face_idx {
                        return self.brep.geom.face_surface.get(flat_idx).copied().flatten();
                    }
                    flat_idx += 1;
                }
            }
        }
        None
    }

    /// Returns the number of edges in the outer wire.
    pub fn edge_count(&self) -> usize {
        self.face()
            .map(|f| f.outer_wire.edges.len())
            .unwrap_or(0)
    }

    /// Returns the number of inner wires (holes).
    pub fn inner_wire_count(&self) -> usize {
        self.face()
            .map(|f| f.inner_wires.len())
            .unwrap_or(0)
    }

    /// Returns the face tolerance.
    pub fn tolerance(&self) -> f64 {
        // Default tolerance, could be extended to read from GeomStore
        1e-6
    }
}

// =============================================================================
// Edge Adaptor
// =============================================================================

/// Adapter providing convenient access to edge data.
///
/// Analogous to OCCT `BRepAdaptor_Curve` / `BRepAdaptor_Edge`.
#[derive(Debug, Clone)]
pub struct EdgeAdaptor<'a> {
    brep: &'a BRep,
    edge_idx: usize,
}

impl<'a> EdgeAdaptor<'a> {
    /// Create a new edge adaptor.
    pub fn new(brep: &'a BRep, edge_idx: usize) -> Self {
        Self { brep, edge_idx }
    }

    /// Returns the edge index.
    pub fn index(&self) -> usize {
        self.edge_idx
    }

    /// Returns the edge topology (start and end vertex indices).
    pub fn edge(&self) -> Option<rcad_kernel::topology::Edge> {
        self.brep.edges.get(self.edge_idx).copied()
    }

    /// Returns the start vertex index.
    pub fn start_vertex(&self) -> Option<usize> {
        self.brep.edges.get(self.edge_idx).map(|e| e.start)
    }

    /// Returns the end vertex index.
    pub fn end_vertex(&self) -> Option<usize> {
        self.brep.edges.get(self.edge_idx).map(|e| e.end)
    }

    /// Returns the 3D curve index for this edge, if available.
    pub fn curve_index(&self) -> Option<usize> {
        self.brep.geom.edge_curve.get(self.edge_idx).copied().flatten()
    }

    /// Returns the parameter range for this edge.
    pub fn parameter_range(&self) -> Option<[f64; 2]> {
        self.brep.geom.edge_curve_range.get(self.edge_idx).copied().flatten()
    }

    /// Returns true if this edge is degenerate.
    pub fn is_degenerate(&self) -> bool {
        self.brep.geom.edge_degenerated.get(self.edge_idx).copied().unwrap_or(false)
    }

    /// Returns true if this edge is closed (start == end vertex).
    pub fn is_closed(&self) -> bool {
        self.brep.edges.get(self.edge_idx)
            .map(|e| e.start == e.end)
            .unwrap_or(false)
    }

    /// Returns the edge tolerance.
    pub fn tolerance(&self) -> f64 {
        self.brep.geom.edge_tolerance.get(self.edge_idx).copied().unwrap_or(1e-6)
    }
}

// =============================================================================
// Vertex Adaptor
// =============================================================================

/// Adapter providing convenient access to vertex data.
///
/// Analogous to OCCT `BRepAdaptor_Point` / `BRep_Tool` for vertices.
#[derive(Debug, Clone)]
pub struct VertexAdaptor<'a> {
    brep: &'a BRep,
    vertex_idx: usize,
}

impl<'a> VertexAdaptor<'a> {
    /// Create a new vertex adaptor.
    pub fn new(brep: &'a BRep, vertex_idx: usize) -> Self {
        Self { brep, vertex_idx }
    }

    /// Returns the vertex index.
    pub fn index(&self) -> usize {
        self.vertex_idx
    }

    /// Returns the 3D point location of the vertex.
    pub fn point(&self) -> Option<glam::DVec3> {
        self.brep.vertices.get(self.vertex_idx).map(|v| v.point)
    }

    /// Returns the vertex tolerance.
    pub fn tolerance(&self) -> f64 {
        self.brep.geom.vertex_tolerance.get(self.vertex_idx).copied().unwrap_or(1e-6)
    }
}

// =============================================================================
// Face Explorer
// =============================================================================

/// Forward-only explorer for faces in a BRep.
///
/// Analogous to OCCT `TopExp_Explorer(shape, TopAbs_FACE)`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::FaceExplorer;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// let mut explorer = FaceExplorer::new(&brep);
/// let mut faces = Vec::new();
/// while let Some(idx) = explorer.next() {
///     faces.push(idx);
/// }
/// assert_eq!(faces.len(), 6);
/// ```
#[derive(Debug, Clone)]
pub struct FaceExplorer<'a> {
    brep: &'a BRep,
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    flat_idx: usize,
    current: Option<FaceAdaptor<'a>>,
}

impl<'a> FaceExplorer<'a> {
    /// Create a new face explorer.
    pub fn new(brep: &'a BRep) -> Self {
        Self {
            brep,
            solid_idx: 0,
            shell_idx: 0,
            face_idx: 0,
            flat_idx: 0,
            current: None,
        }
    }

    /// Advance to the next face and return its flattened index.
    ///
    /// Returns `None` when all faces have been visited.
    pub fn next(&mut self) -> Option<usize> {
        loop {
            // Try to get current solid
            let solid = self.brep.solids.get(self.solid_idx)?;

            // Try to get current shell
            let shell = solid.shells.get(self.shell_idx);

            match shell {
                Some(sh) => {
                    // Try to get current face
                    if self.face_idx < sh.faces.len() {
                        let flat_idx = self.flat_idx;
                        self.current = Some(FaceAdaptor::new(
                            self.brep,
                            self.face_idx,
                            self.shell_idx,
                            self.solid_idx,
                        ));
                        self.face_idx += 1;
                        self.flat_idx += 1;
                        return Some(flat_idx);
                    } else {
                        // Move to next shell
                        self.shell_idx += 1;
                        self.face_idx = 0;
                    }
                }
                None => {
                    // Move to next solid
                    self.solid_idx += 1;
                    self.shell_idx = 0;
                    self.face_idx = 0;
                }
            }
        }
    }

    /// Returns the adaptor for the current face.
    ///
    /// Only valid after a successful call to `next()`.
    pub fn current_adaptor(&self) -> Option<&FaceAdaptor<'a>> {
        self.current.as_ref()
    }

    /// Reset the explorer to start from the beginning.
    pub fn reset(&mut self) {
        self.solid_idx = 0;
        self.shell_idx = 0;
        self.face_idx = 0;
        self.flat_idx = 0;
        self.current = None;
    }
}

// =============================================================================
// Edge Explorer
// =============================================================================

/// Forward-only explorer for edges in a BRep.
///
/// Analogous to OCCT `TopExp_Explorer(shape, TopAbs_EDGE)`.
#[derive(Debug, Clone)]
pub struct EdgeExplorer<'a> {
    brep: &'a BRep,
    edge_idx: usize,
    current: Option<EdgeAdaptor<'a>>,
}

impl<'a> EdgeExplorer<'a> {
    /// Create a new edge explorer.
    pub fn new(brep: &'a BRep) -> Self {
        Self {
            brep,
            edge_idx: 0,
            current: None,
        }
    }

    /// Advance to the next edge and return its index.
    ///
    /// Returns `None` when all edges have been visited.
    pub fn next(&mut self) -> Option<usize> {
        if self.edge_idx < self.brep.edges.len() {
            let idx = self.edge_idx;
            self.current = Some(EdgeAdaptor::new(self.brep, idx));
            self.edge_idx += 1;
            Some(idx)
        } else {
            self.current = None;
            None
        }
    }

    /// Returns the adaptor for the current edge.
    ///
    /// Only valid after a successful call to `next()`.
    pub fn current_adaptor(&self) -> Option<&EdgeAdaptor<'a>> {
        self.current.as_ref()
    }

    /// Reset the explorer to start from the beginning.
    pub fn reset(&mut self) {
        self.edge_idx = 0;
        self.current = None;
    }
}

// =============================================================================
// Vertex Explorer
// =============================================================================

/// Forward-only explorer for vertices in a BRep.
///
/// Analogous to OCCT `TopExp_Explorer(shape, TopAbs_VERTEX)`.
#[derive(Debug, Clone)]
pub struct VertexExplorer<'a> {
    brep: &'a BRep,
    vertex_idx: usize,
}

impl<'a> VertexExplorer<'a> {
    /// Create a new vertex explorer.
    pub fn new(brep: &'a BRep) -> Self {
        Self {
            brep,
            vertex_idx: 0,
        }
    }

    /// Advance to the next vertex and return its index.
    ///
    /// Returns `None` when all vertices have been visited.
    pub fn next(&mut self) -> Option<usize> {
        if self.vertex_idx < self.brep.vertices.len() {
            let idx = self.vertex_idx;
            self.vertex_idx += 1;
            Some(idx)
        } else {
            None
        }
    }

    /// Reset the explorer to start from the beginning.
    pub fn reset(&mut self) {
        self.vertex_idx = 0;
    }
}

// =============================================================================
// Wire Explorer
// =============================================================================

/// An edge reference with orientation from wire traversal.
#[derive(Debug, Clone, Copy)]
pub struct OrientedEdge {
    /// Edge index in `BRep.edges`.
    pub edge_idx: usize,
    /// True if traversed in forward direction (start to end).
    pub forward: bool,
}

impl OrientedEdge {
    /// Create a new oriented edge.
    pub fn new(edge_idx: usize, forward: bool) -> Self {
        Self { edge_idx, forward }
    }
}

/// Forward-only explorer for edges in a wire (face boundary).
///
/// Analogous to OCCT `TopExp_Explorer(face, TopAbs_EDGE)` or `BRepTools_WireExplorer`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::WireExplorer;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// // Explore the wire of the first face
/// let mut explorer = WireExplorer::new(&brep, 0);
/// let mut edge_count = 0;
/// while explorer.next().is_some() {
///     edge_count += 1;
/// }
/// assert_eq!(edge_count, 4); // Each box face has 4 edges
/// ```
#[derive(Debug, Clone)]
pub struct WireExplorer<'a> {
    brep: &'a BRep,
    face_idx: usize,
    wire_idx: usize,  // 0 = outer wire, 1+ = inner wires
    edge_idx: usize,
    current: Option<OrientedEdge>,
}

impl<'a> WireExplorer<'a> {
    /// Create a new wire explorer for a specific face.
    ///
    /// `face_idx` is the flattened face index across all solids/shells.
    pub fn new(brep: &'a BRep, face_idx: usize) -> Self {
        Self {
            brep,
            face_idx,
            wire_idx: 0,
            edge_idx: 0,
            current: None,
        }
    }

    /// Get the face topology for the current face index.
    fn get_face(&self) -> Option<&'a Face> {
        let mut flat_idx = 0;
        for solid in &self.brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    if flat_idx == self.face_idx {
                        return Some(face);
                    }
                    flat_idx += 1;
                }
            }
        }
        None
    }

    /// Advance to the next edge and return its oriented reference.
    ///
    /// Returns `None` when all edges (outer and inner wires) have been visited.
    pub fn next(&mut self) -> Option<OrientedEdge> {
        let face = self.get_face()?;

        loop {
            // Determine which wire we're on
            let wire = if self.wire_idx == 0 {
                &face.outer_wire
            } else {
                face.inner_wires.get(self.wire_idx - 1)?
            };

            // Try to get current edge
            if self.edge_idx < wire.edges.len() {
                let we = &wire.edges[self.edge_idx];
                self.current = Some(OrientedEdge::new(we.idx, we.forward));
                self.edge_idx += 1;
                return self.current;
            } else {
                // Move to next wire
                self.wire_idx += 1;
                self.edge_idx = 0;
            }
        }
    }

    /// Returns the current oriented edge.
    pub fn current(&self) -> Option<OrientedEdge> {
        self.current
    }

    /// Returns true if currently iterating over the outer wire.
    pub fn is_outer_wire(&self) -> bool {
        self.wire_idx == 0
    }

    /// Returns the current wire index (0 = outer, 1+ = inner).
    pub fn wire_index(&self) -> usize {
        self.wire_idx
    }

    /// Reset the explorer to start from the beginning.
    pub fn reset(&mut self) {
        self.wire_idx = 0;
        self.edge_idx = 0;
        self.current = None;
    }
}

// =============================================================================
// Shape Iterator
// =============================================================================

/// Internal state for shape iteration.
#[derive(Debug, Clone)]
struct ShapeIterState {
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    wire_idx: usize,
    edge_idx: usize,
    vertex_idx: usize,
}

impl ShapeIterState {
    fn new() -> Self {
        Self {
            solid_idx: 0,
            shell_idx: 0,
            face_idx: 0,
            wire_idx: 0,
            edge_idx: 0,
            vertex_idx: 0,
        }
    }
}

/// Generic iterator over shapes of a specific type in a BRep.
///
/// Implements `std::iter::Iterator` for idiomatic iteration.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::{ShapeIterator, ShapeType};
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// // Iterate over all faces
/// let faces: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Face).collect();
/// assert_eq!(faces.len(), 6);
///
/// // Iterate over all edges
/// let edges: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Edge).collect();
/// assert_eq!(edges.len(), 12);
///
/// // Iterate over all vertices
/// let vertices: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Vertex).collect();
/// assert_eq!(vertices.len(), 8);
/// ```
#[derive(Debug, Clone)]
pub struct ShapeIterator<'a> {
    brep: &'a BRep,
    shape_type: ShapeType,
    state: ShapeIterState,
    done: bool,
}

impl<'a> ShapeIterator<'a> {
    /// Create a new shape iterator for the given shape type.
    pub fn new(brep: &'a BRep, shape_type: ShapeType) -> Self {
        Self {
            brep,
            shape_type,
            state: ShapeIterState::new(),
            done: false,
        }
    }
}

impl<'a> Iterator for ShapeIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        match self.shape_type {
            ShapeType::Vertex => {
                if self.state.vertex_idx < self.brep.vertices.len() {
                    let idx = self.state.vertex_idx;
                    self.state.vertex_idx += 1;
                    Some(idx)
                } else {
                    self.done = true;
                    None
                }
            }
            ShapeType::Edge => {
                if self.state.edge_idx < self.brep.edges.len() {
                    let idx = self.state.edge_idx;
                    self.state.edge_idx += 1;
                    Some(idx)
                } else {
                    self.done = true;
                    None
                }
            }
            ShapeType::Face => {
                // Iterate over faces in order
                loop {
                    let solid = self.brep.solids.get(self.state.solid_idx)?;
                    let shell = solid.shells.get(self.state.shell_idx);

                    match shell {
                        Some(sh) => {
                            if self.state.face_idx < sh.faces.len() {
                                // Compute flattened index
                                let flat_idx = self.compute_flat_face_index();
                                self.state.face_idx += 1;
                                return Some(flat_idx);
                            } else {
                                self.state.shell_idx += 1;
                                self.state.face_idx = 0;
                            }
                        }
                        None => {
                            self.state.solid_idx += 1;
                            self.state.shell_idx = 0;
                            self.state.face_idx = 0;
                        }
                    }
                }
            }
            ShapeType::Shell => {
                // Iterate over shells
                loop {
                    let solid = self.brep.solids.get(self.state.solid_idx)?;
                    if self.state.shell_idx < solid.shells.len() {
                        // Compute flattened shell index
                        let flat_idx = self.compute_flat_shell_index();
                        self.state.shell_idx += 1;
                        return Some(flat_idx);
                    } else {
                        self.state.solid_idx += 1;
                        self.state.shell_idx = 0;
                    }
                }
            }
            ShapeType::Solid => {
                if self.state.solid_idx < self.brep.solids.len() {
                    let idx = self.state.solid_idx;
                    self.state.solid_idx += 1;
                    Some(idx)
                } else {
                    self.done = true;
                    None
                }
            }
            ShapeType::Wire => {
                // Iterate over all wires (outer + inner) in face order
                loop {
                    // Get current face
                    let face = self.get_current_face()?;

                    // Determine which wire
                    if self.state.wire_idx == 0 {
                        // Outer wire
                        let wire_idx = self.compute_flat_wire_index();
                        self.state.wire_idx += 1;
                        return Some(wire_idx);
                    } else if self.state.wire_idx - 1 < face.inner_wires.len() {
                        // Inner wire
                        let wire_idx = self.compute_flat_wire_index();
                        self.state.wire_idx += 1;
                        return Some(wire_idx);
                    } else {
                        // Move to next face
                        self.advance_face();
                        self.state.wire_idx = 0;
                    }
                }
            }
            ShapeType::Compound | ShapeType::CompSolid | ShapeType::Empty => {
                self.done = true;
                None
            }
        }
    }
}

impl<'a> ShapeIterator<'a> {
    /// Compute the flattened face index for the current state.
    fn compute_flat_face_index(&self) -> usize {
        let mut flat_idx = 0;
        for (si, solid) in self.brep.solids.iter().enumerate() {
            for (shi, shell) in solid.shells.iter().enumerate() {
                if si < self.state.solid_idx
                    || (si == self.state.solid_idx && shi < self.state.shell_idx)
                {
                    flat_idx += shell.faces.len();
                } else if si == self.state.solid_idx && shi == self.state.shell_idx {
                    flat_idx += self.state.face_idx;
                    break;
                }
            }
        }
        flat_idx
    }

    /// Compute the flattened shell index for the current state.
    fn compute_flat_shell_index(&self) -> usize {
        let mut flat_idx = 0;
        for (si, solid) in self.brep.solids.iter().enumerate() {
            if si < self.state.solid_idx {
                flat_idx += solid.shells.len();
            } else if si == self.state.solid_idx {
                flat_idx += self.state.shell_idx;
                break;
            }
        }
        flat_idx
    }

    /// Compute the flattened wire index for the current state.
    fn compute_flat_wire_index(&self) -> usize {
        let mut flat_idx = 0;
        for (si, solid) in self.brep.solids.iter().enumerate() {
            for (shi, shell) in solid.shells.iter().enumerate() {
                for (fi, face) in shell.faces.iter().enumerate() {
                    if si < self.state.solid_idx
                        || (si == self.state.solid_idx && shi < self.state.shell_idx)
                        || (si == self.state.solid_idx
                            && shi == self.state.shell_idx
                            && fi < self.state.face_idx)
                    {
                        flat_idx += 1 + face.inner_wires.len(); // outer + inner wires
                    } else if si == self.state.solid_idx
                        && shi == self.state.shell_idx
                        && fi == self.state.face_idx
                    {
                        flat_idx += self.state.wire_idx;
                        break;
                    }
                }
            }
        }
        flat_idx
    }

    /// Get the current face based on iterator state.
    fn get_current_face(&self) -> Option<&'a Face> {
        self.brep.solids.get(self.state.solid_idx)
            .and_then(|s| s.shells.get(self.state.shell_idx))
            .and_then(|sh| sh.faces.get(self.state.face_idx))
    }

    /// Advance to the next face.
    fn advance_face(&mut self) {
        let solid = match self.brep.solids.get(self.state.solid_idx) {
            Some(s) => s,
            None => return,
        };

        let shell = match solid.shells.get(self.state.shell_idx) {
            Some(sh) => sh,
            None => return,
        };

        if self.state.face_idx + 1 < shell.faces.len() {
            self.state.face_idx += 1;
        } else if self.state.shell_idx + 1 < solid.shells.len() {
            self.state.shell_idx += 1;
            self.state.face_idx = 0;
        } else if self.state.solid_idx + 1 < self.brep.solids.len() {
            self.state.solid_idx += 1;
            self.state.shell_idx = 0;
            self.state.face_idx = 0;
        } else {
            // End of iteration - advance past the last face so get_current_face returns None
            self.state.face_idx += 1;
        }
    }
}

// =============================================================================
// Topology Queries
// =============================================================================

/// Returns all edge indices referenced by a face (including inner wires).
///
/// Duplicate edge indices are preserved as they appear in the wire.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::edges_of_face;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// let edges = edges_of_face(&brep, 0);
/// assert_eq!(edges.len(), 4);
/// ```
pub fn edges_of_face(brep: &BRep, face_idx: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let mut flat_idx = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if flat_idx == face_idx {
                    // Add edges from outer wire
                    for we in &face.outer_wire.edges {
                        result.push(we.idx);
                    }
                    // Add edges from inner wires
                    for inner in &face.inner_wires {
                        for we in &inner.edges {
                            result.push(we.idx);
                        }
                    }
                    return result;
                }
                flat_idx += 1;
            }
        }
    }

    result
}

/// Returns all face indices that reference the given edge.
///
/// For a manifold solid, each edge is typically shared by 2 faces.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::faces_of_edge;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// // Each edge of a box is shared by exactly 2 faces
/// let faces = faces_of_edge(&brep, 0);
/// assert_eq!(faces.len(), 2);
/// ```
pub fn faces_of_edge(brep: &BRep, edge_idx: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let mut flat_idx = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Check outer wire
                for we in &face.outer_wire.edges {
                    if we.idx == edge_idx {
                        result.push(flat_idx);
                        break;
                    }
                }
                // Check inner wires (if not already found in outer)
                if !result.contains(&flat_idx) {
                    for inner in &face.inner_wires {
                        for we in &inner.edges {
                            if we.idx == edge_idx {
                                result.push(flat_idx);
                                break;
                            }
                        }
                        if result.contains(&flat_idx) {
                            break;
                        }
                    }
                }
                flat_idx += 1;
            }
        }
    }

    result
}

/// Returns the (start, end) vertex indices of an edge.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::vertices_of_edge;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// let (start, end) = vertices_of_edge(&brep, 0);
/// assert!(start < 8); // Box has 8 vertices (0-7)
/// assert!(end < 8);
/// ```
pub fn vertices_of_edge(brep: &BRep, edge_idx: usize) -> (usize, usize) {
    brep.edges
        .get(edge_idx)
        .map(|e| (e.start, e.end))
        .unwrap_or((usize::MAX, usize::MAX))
}

/// Returns all edge indices that reference the given vertex.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::edges_of_vertex;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// // Each vertex of a box has 3 incident edges
/// let edges = edges_of_vertex(&brep, 0);
/// assert_eq!(edges.len(), 3);
/// ```
pub fn edges_of_vertex(brep: &BRep, vertex_idx: usize) -> Vec<usize> {
    brep.edges
        .iter()
        .enumerate()
        .filter(|(_, e)| e.start == vertex_idx || e.end == vertex_idx)
        .map(|(ei, _)| ei)
        .collect()
}

/// Returns all faces that share the given vertex (through their edges).
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::faces_of_vertex;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// // Each vertex of a box is shared by 3 faces
/// let faces = faces_of_vertex(&brep, 0);
/// assert_eq!(faces.len(), 3);
/// ```
pub fn faces_of_vertex(brep: &BRep, vertex_idx: usize) -> Vec<usize> {
    let mut result = Vec::new();

    // Get all edges that reference this vertex
    let vertex_edges = edges_of_vertex(brep, vertex_idx);

    // Get all faces that reference these edges
    for edge_idx in vertex_edges {
        for face_idx in faces_of_edge(brep, edge_idx) {
            if !result.contains(&face_idx) {
                result.push(face_idx);
            }
        }
    }

    result
}

/// Returns the number of faces in a BRep.
pub fn face_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum()
}

/// Returns the number of shells in a BRep.
pub fn shell_count(brep: &BRep) -> usize {
    brep.solids.iter().map(|s| s.shells.len()).sum()
}

/// Returns the number of wires in a BRep (including inner wires).
pub fn wire_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .map(|f| 1 + f.inner_wires.len())
        .sum()
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::PrimitiveSolid;

    fn create_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        })
    }

    fn create_cylinder() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        })
    }

    // -------------------------------------------------------------------------
    // FaceExplorer tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_face_explorer_box() {
        let brep = create_box();
        let mut explorer = FaceExplorer::new(&brep);

        let mut faces = Vec::new();
        while let Some(idx) = explorer.next() {
            faces.push(idx);
            assert!(explorer.current_adaptor().is_some());
        }

        assert_eq!(faces.len(), 6);
        assert!(explorer.next().is_none());
    }

    #[test]
    fn test_face_explorer_cylinder() {
        let brep = create_cylinder();
        let mut explorer = FaceExplorer::new(&brep);

        let mut faces = Vec::new();
        while let Some(idx) = explorer.next() {
            faces.push(idx);
        }

        assert_eq!(faces.len(), 3); // Cylinder: lateral + top + bottom
    }

    #[test]
    fn test_face_explorer_reset() {
        let brep = create_box();
        let mut explorer = FaceExplorer::new(&brep);

        let mut count1 = 0;
        while explorer.next().is_some() {
            count1 += 1;
        }

        explorer.reset();

        let mut count2 = 0;
        while explorer.next().is_some() {
            count2 += 1;
        }

        assert_eq!(count1, count2);
        assert_eq!(count1, 6);
    }

    #[test]
    fn test_face_adaptor() {
        let brep = create_box();
        let mut explorer = FaceExplorer::new(&brep);

        while let Some(_) = explorer.next() {
            let adaptor = explorer.current_adaptor().unwrap();
            assert!(adaptor.edge_count() > 0);
            assert!(adaptor.face().is_some());
        }
    }

    // -------------------------------------------------------------------------
    // EdgeExplorer tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_edge_explorer_box() {
        let brep = create_box();
        let mut explorer = EdgeExplorer::new(&brep);

        let mut edges = Vec::new();
        while let Some(idx) = explorer.next() {
            edges.push(idx);
            assert!(explorer.current_adaptor().is_some());
        }

        assert_eq!(edges.len(), 12);
        assert!(explorer.next().is_none());
    }

    #[test]
    fn test_edge_adaptor() {
        let brep = create_box();
        let mut explorer = EdgeExplorer::new(&brep);

        while let Some(_) = explorer.next() {
            let adaptor = explorer.current_adaptor().unwrap();
            assert!(adaptor.edge().is_some());
            assert!(adaptor.start_vertex().is_some());
            assert!(adaptor.end_vertex().is_some());
        }
    }

    // -------------------------------------------------------------------------
    // VertexExplorer tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_vertex_explorer_box() {
        let brep = create_box();
        let mut explorer = VertexExplorer::new(&brep);

        let mut vertices = Vec::new();
        while let Some(idx) = explorer.next() {
            vertices.push(idx);
        }

        assert_eq!(vertices.len(), 8);
        assert!(explorer.next().is_none());
    }

    #[test]
    fn test_vertex_explorer_reset() {
        let brep = create_box();
        let mut explorer = VertexExplorer::new(&brep);

        let mut count1 = 0;
        while explorer.next().is_some() {
            count1 += 1;
        }

        explorer.reset();

        let mut count2 = 0;
        while explorer.next().is_some() {
            count2 += 1;
        }

        assert_eq!(count1, count2);
    }

    // -------------------------------------------------------------------------
    // WireExplorer tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_wire_explorer_box_face() {
        let brep = create_box();

        // Each face of a box has 4 edges in its outer wire
        for face_idx in 0..6 {
            let mut explorer = WireExplorer::new(&brep, face_idx);
            let mut edge_count = 0;

            while explorer.next().is_some() {
                edge_count += 1;
                assert!(explorer.is_outer_wire());
            }

            assert_eq!(edge_count, 4);
        }
    }

    #[test]
    fn test_wire_explorer_oriented_edge() {
        let brep = create_box();
        let mut explorer = WireExplorer::new(&brep, 0);

        let mut edges = Vec::new();
        while let Some(oriented) = explorer.next() {
            edges.push(oriented);
        }

        assert_eq!(edges.len(), 4);
        // All edges should have valid indices
        for e in &edges {
            assert!(e.edge_idx < 12);
        }
    }

    // -------------------------------------------------------------------------
    // ShapeIterator tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_shape_iterator_faces() {
        let brep = create_box();
        let faces: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Face).collect();
        assert_eq!(faces.len(), 6);
    }

    #[test]
    fn test_shape_iterator_edges() {
        let brep = create_box();
        let edges: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Edge).collect();
        assert_eq!(edges.len(), 12);
    }

    #[test]
    fn test_shape_iterator_vertices() {
        let brep = create_box();
        let vertices: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Vertex).collect();
        assert_eq!(vertices.len(), 8);
    }

    #[test]
    fn test_shape_iterator_solids() {
        let brep = create_box();
        let solids: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Solid).collect();
        assert_eq!(solids.len(), 1);
    }

    #[test]
    fn test_shape_iterator_shells() {
        let brep = create_box();
        let shells: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Shell).collect();
        assert_eq!(shells.len(), 1);
    }

    #[test]
    fn test_shape_iterator_wires() {
        let brep = create_box();
        let wires: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Wire).collect();
        assert_eq!(wires.len(), 6); // 6 faces, each with 1 outer wire
    }

    // -------------------------------------------------------------------------
    // Topology query tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_edges_of_face() {
        let brep = create_box();

        for face_idx in 0..6 {
            let edges = edges_of_face(&brep, face_idx);
            assert_eq!(edges.len(), 4);
        }
    }

    #[test]
    fn test_faces_of_edge() {
        let brep = create_box();

        // Each edge of a box is shared by exactly 2 faces
        for edge_idx in 0..12 {
            let faces = faces_of_edge(&brep, edge_idx);
            assert_eq!(faces.len(), 2);
        }
    }

    #[test]
    fn test_vertices_of_edge() {
        let brep = create_box();

        let (start, end) = vertices_of_edge(&brep, 0);
        assert!(start < 8);
        assert!(end < 8);
        assert_ne!(start, end); // Non-degenerate edge
    }

    #[test]
    fn test_edges_of_vertex() {
        let brep = create_box();

        // Each vertex of a box has 3 incident edges
        for vertex_idx in 0..8 {
            let edges = edges_of_vertex(&brep, vertex_idx);
            assert_eq!(edges.len(), 3);
        }
    }

    #[test]
    fn test_faces_of_vertex() {
        let brep = create_box();

        // Each vertex of a box is shared by 3 faces
        for vertex_idx in 0..8 {
            let faces = faces_of_vertex(&brep, vertex_idx);
            assert_eq!(faces.len(), 3);
        }
    }

    // -------------------------------------------------------------------------
    // Count tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_face_count() {
        let brep = create_box();
        assert_eq!(face_count(&brep), 6);

        let cylinder = create_cylinder();
        assert_eq!(face_count(&cylinder), 3);
    }

    #[test]
    fn test_shell_count() {
        let brep = create_box();
        assert_eq!(shell_count(&brep), 1);
    }

    #[test]
    fn test_wire_count() {
        let brep = create_box();
        assert_eq!(wire_count(&brep), 6);
    }

    // -------------------------------------------------------------------------
    // Empty BRep tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_empty_brep() {
        let brep = BRep::new();

        let mut face_explorer = FaceExplorer::new(&brep);
        assert!(face_explorer.next().is_none());

        let mut edge_explorer = EdgeExplorer::new(&brep);
        assert!(edge_explorer.next().is_none());

        let mut vertex_explorer = VertexExplorer::new(&brep);
        assert!(vertex_explorer.next().is_none());

        let faces: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Face).collect();
        assert!(faces.is_empty());
    }

    // -------------------------------------------------------------------------
    // Cylinder seam edge test
    // -------------------------------------------------------------------------

    #[test]
    fn test_cylinder_seam_edge() {
        let brep = create_cylinder();

        // Count faces of each edge
        let mut edges_with_two_faces = 0;
        let mut edges_with_one_face = 0;

        for edge_idx in 0..brep.edges.len() {
            let faces = faces_of_edge(&brep, edge_idx);
            if faces.len() == 2 {
                edges_with_two_faces += 1;
            } else if faces.len() == 1 {
                edges_with_one_face += 1;
            }
        }

        // Cylinder has seam edges that appear once in the face's wire
        // and circular edges at top/bottom
        assert!(edges_with_two_faces > 0);
    }

    // -------------------------------------------------------------------------
    // VertexAdaptor tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_vertex_adaptor() {
        let brep = create_box();

        for idx in 0..8 {
            let adaptor = VertexAdaptor::new(&brep, idx);
            assert_eq!(adaptor.index(), idx);
            assert!(adaptor.point().is_some());
        }
    }
}
