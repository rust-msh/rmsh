//! Graph-topology wrapper for `BRep`, analogous to OCCT's `BRepGraph` module.
//!
//! `BRepGraph` pre-computes and caches all adjacency relations so that repeated
//! topology queries run in O(1) instead of the O(n) scans performed by the
//! free functions in [`crate::topo_query`].
//!
//! # What is cached
//!
//! | Query | OCCT equivalent |
//! |---|---|
//! | edge → adjacent face indices | `TopExp::MapShapesAndAncestors` |
//! | face → edge indices | `TopExp::MapShapes` |
//! | vertex → adjacent edge indices | `TopExp::MapShapesAndAncestors` |
//! | vertex → adjacent face indices | derived from above |
//! | edge → (start_vertex, end_vertex) | `BRep_Tool::Vertices` |
//!
//! # Mutation tracking
//!
//! Each vertex, edge, and face has a corresponding dirty bit.  Algorithms that
//! mutate a BRep entity can set these bits via the `mark_*_modified` methods,
//! and downstream steps can query them to decide what needs re-evaluation.
//! This is a lightweight analogue to OCCT's `BRepCheck_Analyzer` invalidation
//! mechanism and the history-event bus.
//!
//! # Traversal
//!
//! `BRepGraph` provides DFS and BFS iterators over face/edge/vertex connectivity:
//! - [`BRepGraph::dfs_faces`] — visits all faces reachable (via shared edges) from a seed face.
//! - [`BRepGraph::bfs_faces`] — same graph, breadth-first order.
//! - [`BRepGraph::dfs_edges_from_vertex`] — walks edge-adjacency from a seed vertex DFS.
//!
//! # Examples
//!
//! ```rust
//! use rcad_kernel::{BRep, BRepGraph, PrimitiveSolid};
//!
//! let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
//! let graph = BRepGraph::from_brep(&brep);
//!
//! // O(1) adjacency
//! let adj_faces = graph.edge_adjacent_faces(0);
//! assert_eq!(adj_faces.len(), 2);
//!
//! // Traversal: all 6 faces of a box are connected
//! let visited = graph.bfs_faces(0).collect::<Vec<_>>();
//! assert_eq!(visited.len(), 6);
//!
//! // Manifold / closed inspection
//! assert!(graph.is_manifold());
//! assert!(graph.is_closed());
//! ```

use crate::{BRep, Edge, Face, PCurve, Vertex};
use glam::DVec3;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet, VecDeque};

// ─────────────────────────────────────────────────────────────────────────────
// BRepGraph
// ─────────────────────────────────────────────────────────────────────────────

/// Cached graph-topology view of a [`BRep`].
///
/// Build with [`BRepGraph::from_brep`].  The graph is a **snapshot** — it does
/// not update automatically when the source `BRep` is mutated.  Call
/// `from_brep` again to refresh.
#[derive(Debug, Clone)]
pub struct BRepGraph {
    // ── Counts ────────────────────────────────────────────────────────────────
    pub vertex_count: usize,
    pub edge_count: usize,
    pub face_count: usize,

    // ── Adjacency tables (indexed by the corresponding entity index) ──────────
    /// `edge_to_faces[edge_idx]` = list of flat face indices that include this edge
    /// in their outer wire.
    edge_to_faces: Vec<Vec<usize>>,

    /// `face_to_edges[face_idx]` = list of edge indices in the outer wire of that face.
    face_to_edges: Vec<Vec<usize>>,

    /// `vertex_to_edges[vertex_idx]` = list of edge indices where start == vertex_idx
    /// or end == vertex_idx.
    vertex_to_edges: Vec<Vec<usize>>,

    /// `vertex_to_faces[vertex_idx]` = list of flat face indices whose outer wire
    /// references a wire-edge that touches this vertex.
    vertex_to_faces: Vec<Vec<usize>>,

    /// `edge_endpoints[edge_idx]` = (start_vertex_idx, end_vertex_idx).
    edge_endpoints: Vec<(usize, usize)>,

    // ── Dirty / modified bits ─────────────────────────────────────────────────
    /// Dirty flags for vertices (by vertex index).
    vertex_dirty: Vec<bool>,
    /// Dirty flags for edges (by edge index).
    edge_dirty: Vec<bool>,
    /// Dirty flags for faces (by flat face index).
    face_dirty: Vec<bool>,
}

/// Programmatic builder for a [`BRepGraph`] without first materializing a full
/// [`BRep`].
///
/// Analogous to OCCT `BRepGraph_Builder`.
#[derive(Debug, Clone)]
pub struct BRepGraphBuilder {
    vertex_count: usize,
    edge_count: usize,
    face_count: usize,
    edge_to_faces: Vec<Vec<usize>>,
    face_to_edges: Vec<Vec<usize>>,
    edge_endpoints: Vec<(usize, usize)>,
    vertex_dirty: Vec<bool>,
    edge_dirty: Vec<bool>,
    face_dirty: Vec<bool>,
}

/// Geometry/topology accessor over graph node indices.
///
/// This bridges flat `BRepGraph` indices back to the source [`BRep`] entities
/// and geometry pools.
///
/// Analogous to OCCT `BRepGraph_Tool` / `BRep_Tool` accessors.
#[derive(Debug, Clone, Copy)]
pub struct BRepGraphTool<'a> {
    graph: &'a BRepGraph,
    brep: &'a BRep,
}

/// Non-manifold topology summary derived from edge-face adjacency.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct NonManifoldSummary {
    /// Edge indices with exactly one adjacent face.
    pub boundary_edges: Vec<usize>,
    /// Edge indices with zero adjacent faces.
    pub orphan_edges: Vec<usize>,
    /// Edge indices with more than two adjacent faces.
    pub multi_face_edges: Vec<usize>,
    /// Vertex indices touched by at least one multi-face edge.
    pub non_manifold_vertices: Vec<usize>,
}

impl NonManifoldSummary {
    pub fn is_clean(&self) -> bool {
        self.boundary_edges.is_empty() && self.orphan_edges.is_empty() && self.multi_face_edges.is_empty()
    }
}

/// Categorised repair hints for a non-manifold BRep.
///
/// Each hint describes one problem and suggests a concrete remedy action.
/// Analogous to OCCT `BRepCheck_Analyzer` diagnostic entries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RepairHint {
    /// Two boundary edges share both endpoints — they can be stitched together.
    StitchablePair {
        edge_a: usize,
        edge_b: usize,
        face_a: usize,
        face_b: usize,
    },
    /// A boundary edge has no stitch candidate — the hole must be capped.
    UnmatchedBoundaryEdge { edge_idx: usize, face_idx: usize },
    /// An orphan edge is attached to no face — it should be removed.
    OrphanEdge { edge_idx: usize },
    /// An edge is shared by more than two faces — it must be split into
    /// separate copies, one per face pair.
    MultiManifoldEdge { edge_idx: usize, face_count: usize },
    /// A vertex lies on a multi-face edge — it may require duplication when
    /// the surrounding edge is split.
    NonManifoldVertex { vertex_idx: usize, connected_multi_edges: Vec<usize> },
}

/// Detailed actionable repair hints derived from a `NonManifoldSummary`.
///
/// Build with [`BRepGraph::repair_hints`].
#[derive(Debug, Clone, Default)]
pub struct ManifoldRepairHints {
    pub hints: Vec<RepairHint>,
}

impl ManifoldRepairHints {
    /// Returns `true` when there are no repair items.
    pub fn is_empty(&self) -> bool {
        self.hints.is_empty()
    }

    /// All hints of the `StitchablePair` variant.
    pub fn stitchable_pairs(&self) -> impl Iterator<Item = &RepairHint> {
        self.hints.iter().filter(|h| matches!(h, RepairHint::StitchablePair { .. }))
    }

    /// All hints of the `OrphanEdge` variant.
    pub fn orphan_edges(&self) -> impl Iterator<Item = &RepairHint> {
        self.hints.iter().filter(|h| matches!(h, RepairHint::OrphanEdge { .. }))
    }
}

impl BRepGraph {
    // ── Construction ──────────────────────────────────────────────────────────

    /// Build a `BRepGraph` from the current state of `brep`.
    ///
    /// Iterates all topology entities once (O(V + E + F + wire_edges)) to
    /// populate the adjacency tables.  Subsequent adjacency queries are O(1).
    pub fn from_brep(brep: &BRep) -> Self {
        let vc = brep.vertices.len();
        let ec = brep.edges.len();

        // ── Pre-fill edge endpoints ───────────────────────────────────────────
        let edge_endpoints: Vec<(usize, usize)> = brep
            .edges
            .iter()
            .map(|e| (e.start, e.end))
            .collect();

        // ── vertex → edges ────────────────────────────────────────────────────
        let mut vertex_to_edges: Vec<Vec<usize>> = vec![Vec::new(); vc];
        for (ei, ep) in edge_endpoints.iter().enumerate() {
            if ep.0 < vc {
                vertex_to_edges[ep.0].push(ei);
            }
            if ep.1 < vc && ep.1 != ep.0 {
                vertex_to_edges[ep.1].push(ei);
            }
        }

        // ── Build flat face list and face/edge/vertex adjacency ───────────────
        // Count faces first.
        let fc: usize = brep
            .solids
            .iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();

        let mut edge_to_faces: Vec<Vec<usize>> = vec![Vec::new(); ec];
        let mut face_to_edges: Vec<Vec<usize>> = vec![Vec::new(); fc];
        let mut vertex_to_faces: Vec<Vec<usize>> = vec![Vec::new(); vc];

        let mut flat_fi = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let mut face_edge_set = Vec::new();
                    let mut face_vert_set = HashSet::new();
                    for we in &face.outer_wire.edges {
                        if we.idx < ec {
                            // edge → face
                            if !edge_to_faces[we.idx].contains(&flat_fi) {
                                edge_to_faces[we.idx].push(flat_fi);
                            }
                            // face → edge
                            if !face_edge_set.contains(&we.idx) {
                                face_edge_set.push(we.idx);
                            }
                            // vertices touched by this wire edge
                            let (vs, ve) = edge_endpoints[we.idx];
                            if vs < vc {
                                face_vert_set.insert(vs);
                            }
                            if ve < vc {
                                face_vert_set.insert(ve);
                            }
                        }
                    }
                    // Also walk inner wires.
                    for inner in &face.inner_wires {
                        for we in &inner.edges {
                            if we.idx < ec {
                                if !edge_to_faces[we.idx].contains(&flat_fi) {
                                    edge_to_faces[we.idx].push(flat_fi);
                                }
                                if !face_edge_set.contains(&we.idx) {
                                    face_edge_set.push(we.idx);
                                }
                                let (vs, ve) = edge_endpoints[we.idx];
                                if vs < vc {
                                    face_vert_set.insert(vs);
                                }
                                if ve < vc {
                                    face_vert_set.insert(ve);
                                }
                            }
                        }
                    }
                    face_to_edges[flat_fi] = face_edge_set;
                    // vertex → faces
                    for vi in face_vert_set {
                        if !vertex_to_faces[vi].contains(&flat_fi) {
                            vertex_to_faces[vi].push(flat_fi);
                        }
                    }
                    flat_fi += 1;
                }
            }
        }

        BRepGraph {
            vertex_count: vc,
            edge_count: ec,
            face_count: fc,
            edge_to_faces,
            face_to_edges,
            vertex_to_edges,
            vertex_to_faces,
            edge_endpoints,
            vertex_dirty: vec![false; vc],
            edge_dirty: vec![false; ec],
            face_dirty: vec![false; fc],
        }
    }

    /// Start building a graph programmatically.
    ///
    /// The caller provides the intended vertex/edge/face counts up front, then
    /// fills edge endpoints and edge/face incidence through the returned
    /// [`BRepGraphBuilder`].
    pub fn builder(vertex_count: usize, edge_count: usize, face_count: usize) -> BRepGraphBuilder {
        BRepGraphBuilder::new(vertex_count, edge_count, face_count)
    }

    /// Create a geometry/topology access wrapper over this graph and `brep`.
    pub fn tool<'a>(&'a self, brep: &'a BRep) -> BRepGraphTool<'a> {
        BRepGraphTool::new(self, brep)
    }

    // ── O(1) adjacency queries ────────────────────────────────────────────────

    /// Faces that share `edge_idx` in their wire.  Expected length 2 for a
    /// manifold edge, 1 for a boundary edge, 0 if the index is out of range.
    #[inline]
    pub fn edge_adjacent_faces(&self, edge_idx: usize) -> &[usize] {
        self.edge_to_faces
            .get(edge_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Edge indices in the outer (+inner) wires of `face_idx`.
    #[inline]
    pub fn face_edges(&self, face_idx: usize) -> &[usize] {
        self.face_to_edges
            .get(face_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Edge indices incident to `vertex_idx`.
    #[inline]
    pub fn vertex_adjacent_edges(&self, vertex_idx: usize) -> &[usize] {
        self.vertex_to_edges
            .get(vertex_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Face indices adjacent to `vertex_idx` (faces whose wire touches it).
    #[inline]
    pub fn vertex_adjacent_faces(&self, vertex_idx: usize) -> &[usize] {
        self.vertex_to_faces
            .get(vertex_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// (start_vertex, end_vertex) for `edge_idx`.
    #[inline]
    pub fn edge_endpoints(&self, edge_idx: usize) -> Option<(usize, usize)> {
        self.edge_endpoints.get(edge_idx).copied()
    }

    // ── Manifold / closed inspection ──────────────────────────────────────────

    /// Returns `true` if every edge is shared by exactly 2 faces (no boundary
    /// or non-manifold edges).
    ///
    /// Analogous to `BRepCheck_Shell::Closed()` manifold test in OCCT.
    pub fn is_manifold(&self) -> bool {
        self.edge_to_faces.iter().all(|adj| adj.len() == 2)
    }

    /// Returns `true` if no edge has fewer than 2 adjacent faces (i.e., no
    /// boundary / free edges exist).
    pub fn is_closed(&self) -> bool {
        self.edge_to_faces.iter().all(|adj| adj.len() >= 2)
    }

    /// Returns edge indices where the number of adjacent faces is not exactly 2.
    ///
    /// - `adj.len() == 1` → free (boundary) edge
    /// - `adj.len() > 2`  → non-manifold edge
    /// - `adj.len() == 0` → orphan edge (not used in any face wire)
    pub fn non_manifold_edges(&self) -> Vec<usize> {
        self.edge_to_faces
            .iter()
            .enumerate()
            .filter(|(_, adj)| adj.len() != 2)
            .map(|(ei, _)| ei)
            .collect()
    }

    /// Edge indices with exactly one adjacent face.
    pub fn boundary_edges(&self) -> Vec<usize> {
        self.edge_to_faces
            .iter()
            .enumerate()
            .filter(|(_, adj)| adj.len() == 1)
            .map(|(ei, _)| ei)
            .collect()
    }

    /// Edge indices with zero adjacent faces.
    pub fn orphan_edges(&self) -> Vec<usize> {
        self.edge_to_faces
            .iter()
            .enumerate()
            .filter(|(_, adj)| adj.is_empty())
            .map(|(ei, _)| ei)
            .collect()
    }

    /// Edge indices with more than two adjacent faces.
    pub fn multi_face_edges(&self) -> Vec<usize> {
        self.edge_to_faces
            .iter()
            .enumerate()
            .filter(|(_, adj)| adj.len() > 2)
            .map(|(ei, _)| ei)
            .collect()
    }

    /// Vertex indices touched by multi-face edges (>2 adjacent faces).
    pub fn non_manifold_vertices(&self) -> Vec<usize> {
        let mut verts = HashSet::new();
        for ei in self.multi_face_edges() {
            if let Some((vs, ve)) = self.edge_endpoints(ei) {
                verts.insert(vs);
                verts.insert(ve);
            }
        }
        let mut out: Vec<usize> = verts.into_iter().collect();
        out.sort_unstable();
        out
    }

    /// Combined non-manifold summary report.
    pub fn non_manifold_summary(&self) -> NonManifoldSummary {
        NonManifoldSummary {
            boundary_edges: self.boundary_edges(),
            orphan_edges: self.orphan_edges(),
            multi_face_edges: self.multi_face_edges(),
            non_manifold_vertices: self.non_manifold_vertices(),
        }
    }

    /// The number of faces sharing `edge_idx` (the edge "valence").
    ///
    /// - 0 → orphan edge
    /// - 1 → boundary / free edge
    /// - 2 → manifold edge (expected)
    /// - >2 → non-manifold / T-junction edge
    #[inline]
    pub fn edge_valence(&self, edge_idx: usize) -> usize {
        self.edge_to_faces
            .get(edge_idx)
            .map(|v| v.len())
            .unwrap_or(0)
    }

    /// The number of edges incident to `vertex_idx` (vertex degree).
    #[inline]
    pub fn vertex_degree(&self, vertex_idx: usize) -> usize {
        self.vertex_to_edges
            .get(vertex_idx)
            .map(|v| v.len())
            .unwrap_or(0)
    }

    /// Generate actionable [`ManifoldRepairHints`] for this graph.
    ///
    /// The algorithm classifies each problem edge and groups compatible
    /// boundary edges as stitchable pairs (same orientation·length within 1e-6).
    ///
    /// Analogous to OCCT `BRepCheck_Analyzer` hint emission.
    pub fn repair_hints(&self, brep: &crate::BRep) -> ManifoldRepairHints {
        let mut hints = Vec::new();

        // ── Multi-face (non-manifold) edges ───────────────────────────────────
        let multi = self.multi_face_edges();
        for ei in &multi {
            hints.push(RepairHint::MultiManifoldEdge {
                edge_idx: *ei,
                face_count: self.edge_valence(*ei),
            });
        }

        // ── Non-manifold vertices ─────────────────────────────────────────────
        for vi in self.non_manifold_vertices() {
            let connected: Vec<usize> = self
                .vertex_adjacent_edges(vi)
                .iter()
                .copied()
                .filter(|&ei| self.edge_valence(ei) > 2)
                .collect();
            if !connected.is_empty() {
                hints.push(RepairHint::NonManifoldVertex {
                    vertex_idx: vi,
                    connected_multi_edges: connected,
                });
            }
        }

        // ── Orphan edges ──────────────────────────────────────────────────────
        for ei in self.orphan_edges() {
            hints.push(RepairHint::OrphanEdge { edge_idx: ei });
        }

        // ── Boundary edges: find stitchable pairs ─────────────────────────────
        let boundary = self.boundary_edges();
        let mut paired: HashSet<usize> = HashSet::new();

        for i in 0..boundary.len() {
            let ei = boundary[i];
            if paired.contains(&ei) {
                continue;
            }
            let face_i = self.edge_adjacent_faces(ei).first().copied().unwrap_or(0);
            let (va, vb) = match self.edge_endpoints(ei) {
                Some(ep) => ep,
                None => continue,
            };
            let pa = brep.vertices.get(va).map(|v| v.point).unwrap_or_default();
            let pb = brep.vertices.get(vb).map(|v| v.point).unwrap_or_default();
            let len_i = (pb - pa).length();

            let mut matched = false;
            for j in (i + 1)..boundary.len() {
                let ej = boundary[j];
                if paired.contains(&ej) {
                    continue;
                }
                let (vc2, vd) = match self.edge_endpoints(ej) {
                    Some(ep) => ep,
                    None => continue,
                };
                let pc = brep.vertices.get(vc2).map(|v| v.point).unwrap_or_default();
                let pd = brep.vertices.get(vd).map(|v| v.point).unwrap_or_default();
                let len_j = (pd - pc).length();

                // Stitchable if both endpoints are coincident (within 1e-6),
                // possibly with reversed orientation.
                let direct  = (pa - pc).length() < 1e-6 && (pb - pd).length() < 1e-6;
                let reverse = (pa - pd).length() < 1e-6 && (pb - pc).length() < 1e-6;
                let same_len = (len_i - len_j).abs() < 1e-6;

                if same_len && (direct || reverse) {
                    let face_j = self.edge_adjacent_faces(ej).first().copied().unwrap_or(0);
                    hints.push(RepairHint::StitchablePair {
                        edge_a: ei,
                        edge_b: ej,
                        face_a: face_i,
                        face_b: face_j,
                    });
                    paired.insert(ei);
                    paired.insert(ej);
                    matched = true;
                    break;
                }
            }
            if !matched {
                hints.push(RepairHint::UnmatchedBoundaryEdge {
                    edge_idx: ei,
                    face_idx: face_i,
                });
            }
        }

        ManifoldRepairHints { hints }
    }

    /// Returns vertex indices that are referenced by a number of edges other
    /// than the expected valence.  For a closed manifold each vertex has
    /// valence ≥ 2 (connected by at least 2 edges).
    pub fn low_valence_vertices(&self) -> Vec<usize> {
        self.vertex_to_edges
            .iter()
            .enumerate()
            .filter(|(_, adj)| adj.len() < 2)
            .map(|(vi, _)| vi)
            .collect()
    }

    // ── Dirty / modification tracking ─────────────────────────────────────────

    /// Mark vertex `vertex_idx` as having been modified.
    pub fn mark_vertex_modified(&mut self, vertex_idx: usize) {
        if let Some(d) = self.vertex_dirty.get_mut(vertex_idx) {
            *d = true;
        }
    }

    /// Mark edge `edge_idx` as having been modified.
    pub fn mark_edge_modified(&mut self, edge_idx: usize) {
        if let Some(d) = self.edge_dirty.get_mut(edge_idx) {
            *d = true;
        }
    }

    /// Mark face `face_idx` (flat index) as having been modified.
    pub fn mark_face_modified(&mut self, face_idx: usize) {
        if let Some(d) = self.face_dirty.get_mut(face_idx) {
            *d = true;
        }
    }

    /// Clear all dirty flags.
    pub fn clear_dirty(&mut self) {
        self.vertex_dirty.iter_mut().for_each(|d| *d = false);
        self.edge_dirty.iter_mut().for_each(|d| *d = false);
        self.face_dirty.iter_mut().for_each(|d| *d = false);
    }

    /// All vertex indices currently marked dirty.
    pub fn modified_vertices(&self) -> Vec<usize> {
        self.vertex_dirty
            .iter()
            .enumerate()
            .filter(|(_, d)| **d)
            .map(|(vi, _)| vi)
            .collect()
    }

    /// All edge indices currently marked dirty.
    pub fn modified_edges(&self) -> Vec<usize> {
        self.edge_dirty
            .iter()
            .enumerate()
            .filter(|(_, d)| **d)
            .map(|(ei, _)| ei)
            .collect()
    }

    /// All face indices (flat) currently marked dirty.
    pub fn modified_faces(&self) -> Vec<usize> {
        self.face_dirty
            .iter()
            .enumerate()
            .filter(|(_, d)| **d)
            .map(|(fi, _)| fi)
            .collect()
    }

    /// Returns `true` if any entity is dirty.
    pub fn has_modifications(&self) -> bool {
        self.vertex_dirty.iter().any(|&d| d)
            || self.edge_dirty.iter().any(|&d| d)
            || self.face_dirty.iter().any(|&d| d)
    }

    // ── Accessor methods for checkpoint/rollback (used by rcad-advanced-modeling) ─

    /// Get a reference to the edge-to-faces adjacency table.
    pub fn edge_to_faces_table(&self) -> &[Vec<usize>] {
        &self.edge_to_faces
    }

    /// Get a reference to the face-to-edges adjacency table.
    pub fn face_to_edges_table(&self) -> &[Vec<usize>] {
        &self.face_to_edges
    }

    /// Get a reference to the vertex-to-edges adjacency table.
    pub fn vertex_to_edges_table(&self) -> &[Vec<usize>] {
        &self.vertex_to_edges
    }

    /// Get a reference to the vertex-to-faces adjacency table.
    pub fn vertex_to_faces_table(&self) -> &[Vec<usize>] {
        &self.vertex_to_faces
    }

    /// Get a reference to the edge endpoints table.
    pub fn edge_endpoints_table(&self) -> &[(usize, usize)] {
        &self.edge_endpoints
    }

    /// Get a reference to the vertex dirty flags.
    pub fn vertex_dirty_flags(&self) -> &[bool] {
        &self.vertex_dirty
    }

    /// Get a reference to the edge dirty flags.
    pub fn edge_dirty_flags(&self) -> &[bool] {
        &self.edge_dirty
    }

    /// Get a reference to the face dirty flags.
    pub fn face_dirty_flags(&self) -> &[bool] {
        &self.face_dirty
    }

    // ── Graph traversal ───────────────────────────────────────────────────────

    /// Depth-first traversal of faces, starting from `seed_face_idx`.
    ///
    /// Faces are considered adjacent when they share at least one edge.
    /// Returns a `DfsFaces` iterator that yields flat face indices in DFS
    /// discovery order.
    ///
    /// Analogous to `BRepTools_WireExplorer` / `TopExp_Explorer` over `TopAbs_FACE`.
    pub fn dfs_faces(&self, seed_face_idx: usize) -> DfsFaces<'_> {
        DfsFaces::new(self, seed_face_idx)
    }

    /// Breadth-first traversal of faces, starting from `seed_face_idx`.
    ///
    /// Returns a `BfsFaces` iterator that yields flat face indices in BFS
    /// discovery order.
    pub fn bfs_faces(&self, seed_face_idx: usize) -> BfsFaces<'_> {
        BfsFaces::new(self, seed_face_idx)
    }

    /// Depth-first traversal of edges, starting from `seed_vertex_idx`.
    ///
    /// Visits each edge reachable by following edge-adjacency from the seed
    /// vertex.  Yields edge indices in DFS discovery order.
    pub fn dfs_edges_from_vertex(&self, seed_vertex_idx: usize) -> DfsEdgesFromVertex<'_> {
        DfsEdgesFromVertex::new(self, seed_vertex_idx)
    }
}

impl BRepGraphBuilder {
    /// Create an empty builder with pre-sized incidence tables.
    pub fn new(vertex_count: usize, edge_count: usize, face_count: usize) -> Self {
        Self {
            vertex_count,
            edge_count,
            face_count,
            edge_to_faces: vec![Vec::new(); edge_count],
            face_to_edges: vec![Vec::new(); face_count],
            edge_endpoints: vec![(0, 0); edge_count],
            vertex_dirty: vec![false; vertex_count],
            edge_dirty: vec![false; edge_count],
            face_dirty: vec![false; face_count],
        }
    }

    /// Set the start/end vertices of an edge.
    pub fn set_edge_endpoints(
        &mut self,
        edge_idx: usize,
        start_vertex: usize,
        end_vertex: usize,
    ) -> &mut Self {
        if let Some(slot) = self.edge_endpoints.get_mut(edge_idx) {
            *slot = (start_vertex, end_vertex);
        }
        self
    }

    /// Add an edge→face incidence entry.
    pub fn add_edge_face(&mut self, edge_idx: usize, face_idx: usize) -> &mut Self {
        if let Some(adj) = self.edge_to_faces.get_mut(edge_idx) {
            if !adj.contains(&face_idx) {
                adj.push(face_idx);
            }
        }
        self
    }

    /// Add a face→edge incidence entry.
    pub fn add_face_edge(&mut self, face_idx: usize, edge_idx: usize) -> &mut Self {
        if let Some(adj) = self.face_to_edges.get_mut(face_idx) {
            if !adj.contains(&edge_idx) {
                adj.push(edge_idx);
            }
        }
        self
    }

    /// Mark a vertex dirty in the initial graph state.
    pub fn mark_vertex_modified(&mut self, vertex_idx: usize) -> &mut Self {
        if let Some(dirty) = self.vertex_dirty.get_mut(vertex_idx) {
            *dirty = true;
        }
        self
    }

    /// Mark an edge dirty in the initial graph state.
    pub fn mark_edge_modified(&mut self, edge_idx: usize) -> &mut Self {
        if let Some(dirty) = self.edge_dirty.get_mut(edge_idx) {
            *dirty = true;
        }
        self
    }

    /// Mark a face dirty in the initial graph state.
    pub fn mark_face_modified(&mut self, face_idx: usize) -> &mut Self {
        if let Some(dirty) = self.face_dirty.get_mut(face_idx) {
            *dirty = true;
        }
        self
    }

    /// Build the graph, deriving vertex incidence from edge endpoints and
    /// face-edge incidence, and validating internal consistency.
    pub fn build(self) -> Result<BRepGraph, Vec<String>> {
        let graph = self.build_unchecked();
        let errors = graph.validate_invariants();
        if errors.is_empty() {
            Ok(graph)
        } else {
            Err(errors)
        }
    }

    /// Build the graph without validation.
    pub fn build_unchecked(self) -> BRepGraph {
        let mut vertex_to_edges = vec![Vec::new(); self.vertex_count];
        for (edge_idx, &(start_vertex, end_vertex)) in self.edge_endpoints.iter().enumerate() {
            if start_vertex < self.vertex_count {
                vertex_to_edges[start_vertex].push(edge_idx);
            }
            if end_vertex < self.vertex_count && end_vertex != start_vertex {
                vertex_to_edges[end_vertex].push(edge_idx);
            }
        }

        let mut vertex_to_faces = vec![Vec::new(); self.vertex_count];
        for (face_idx, edges) in self.face_to_edges.iter().enumerate() {
            let mut vertices_on_face = HashSet::new();
            for &edge_idx in edges {
                if let Some(&(start_vertex, end_vertex)) = self.edge_endpoints.get(edge_idx) {
                    if start_vertex < self.vertex_count {
                        vertices_on_face.insert(start_vertex);
                    }
                    if end_vertex < self.vertex_count {
                        vertices_on_face.insert(end_vertex);
                    }
                }
            }
            for vertex_idx in vertices_on_face {
                vertex_to_faces[vertex_idx].push(face_idx);
            }
        }

        BRepGraph {
            vertex_count: self.vertex_count,
            edge_count: self.edge_count,
            face_count: self.face_count,
            edge_to_faces: self.edge_to_faces,
            face_to_edges: self.face_to_edges,
            vertex_to_edges,
            vertex_to_faces,
            edge_endpoints: self.edge_endpoints,
            vertex_dirty: self.vertex_dirty,
            edge_dirty: self.edge_dirty,
            face_dirty: self.face_dirty,
        }
    }
}

impl<'a> BRepGraphTool<'a> {
    pub fn new(graph: &'a BRepGraph, brep: &'a BRep) -> Self {
        Self { graph, brep }
    }

    /// Raw vertex by graph vertex index.
    pub fn vertex(&self, vertex_idx: usize) -> Option<&'a Vertex> {
        self.brep.vertices.get(vertex_idx)
    }

    /// Raw edge by graph edge index.
    pub fn edge(&self, edge_idx: usize) -> Option<&'a Edge> {
        self.brep.edges.get(edge_idx)
    }

    /// Raw face by graph flat-face index.
    pub fn face(&self, face_idx: usize) -> Option<&'a Face> {
        self.face_location(face_idx)
            .and_then(|(solid_idx, shell_idx, local_face_idx)| {
                self.brep
                    .solids
                    .get(solid_idx)?
                    .shells
                    .get(shell_idx)?
                    .faces
                    .get(local_face_idx)
            })
    }

    /// Map a flat face index back to `(solid_idx, shell_idx, local_face_idx)`.
    pub fn face_location(&self, face_idx: usize) -> Option<(usize, usize, usize)> {
        let mut flat_face_idx = 0usize;
        for (solid_idx, solid) in self.brep.solids.iter().enumerate() {
            for (shell_idx, shell) in solid.shells.iter().enumerate() {
                for local_face_idx in 0..shell.faces.len() {
                    if flat_face_idx == face_idx {
                        return Some((solid_idx, shell_idx, local_face_idx));
                    }
                    flat_face_idx += 1;
                }
            }
        }
        None
    }

    /// Vertex position by graph vertex index.
    pub fn vertex_point(&self, vertex_idx: usize) -> Option<DVec3> {
        self.vertex(vertex_idx).map(|vertex| vertex.point)
    }

    /// Edge endpoint positions by graph edge index.
    pub fn edge_points(&self, edge_idx: usize) -> Option<(DVec3, DVec3)> {
        let (start_vertex, end_vertex) = self.graph.edge_endpoints(edge_idx)?;
        Some((self.vertex_point(start_vertex)?, self.vertex_point(end_vertex)?))
    }

    /// Face normal by graph flat-face index.
    pub fn face_normal(&self, face_idx: usize) -> Option<DVec3> {
        self.face(face_idx).map(|face| face.normal)
    }

    /// Associated 3D curve for an edge, if the source `BRep` carries one.
    pub fn edge_curve(&self, edge_idx: usize) -> Option<&'a crate::geom::Curve3> {
        let curve_idx = self.brep.geom.edge_curve.get(edge_idx).copied().flatten()?;
        self.brep.geom.curves.get(curve_idx)
    }

    /// Associated surface for a face, if the source `BRep` carries one.
    pub fn face_surface(&self, face_idx: usize) -> Option<&'a crate::geom::Surface3> {
        let surface_idx = self.brep.geom.face_surface.get(face_idx).copied().flatten()?;
        self.brep.geom.surfaces.get(surface_idx)
    }

    /// PCurves attached to an edge.
    pub fn edge_pcurves(&self, edge_idx: usize) -> Option<&'a [PCurve]> {
        self.brep.geom.edge_pcurves.get(edge_idx).map(|pcurves| pcurves.as_slice())
    }

    /// Edge parameter range on its 3D curve, if present.
    pub fn edge_curve_range(&self, edge_idx: usize) -> Option<[f64; 2]> {
        self.brep.geom.edge_curve_range.get(edge_idx).copied().flatten()
    }

    /// Face surface parameter domain override, if present.
    pub fn face_surface_range(&self, face_idx: usize) -> Option<[f64; 4]> {
        self.brep.geom.face_surface_range.get(face_idx).copied().flatten()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Iterators
// ─────────────────────────────────────────────────────────────────────────────

/// DFS iterator over faces reachable from a seed face via shared edges.
pub struct DfsFaces<'g> {
    graph: &'g BRepGraph,
    stack: Vec<usize>,
    visited: HashSet<usize>,
}

impl<'g> DfsFaces<'g> {
    fn new(graph: &'g BRepGraph, seed: usize) -> Self {
        let mut visited = HashSet::new();
        let mut stack = Vec::new();
        if seed < graph.face_count {
            stack.push(seed);
            visited.insert(seed);
        }
        DfsFaces { graph, stack, visited }
    }
}

impl Iterator for DfsFaces<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        let fi = self.stack.pop()?;
        // Push unvisited neighbors: faces reachable via edges of `fi`.
        for &ei in self.graph.face_edges(fi) {
            for &adj_fi in self.graph.edge_adjacent_faces(ei) {
                if self.visited.insert(adj_fi) {
                    self.stack.push(adj_fi);
                }
            }
        }
        Some(fi)
    }
}

/// BFS iterator over faces reachable from a seed face via shared edges.
pub struct BfsFaces<'g> {
    graph: &'g BRepGraph,
    queue: VecDeque<usize>,
    visited: HashSet<usize>,
}

impl<'g> BfsFaces<'g> {
    fn new(graph: &'g BRepGraph, seed: usize) -> Self {
        let mut visited = HashSet::new();
        let mut queue = VecDeque::new();
        if seed < graph.face_count {
            queue.push_back(seed);
            visited.insert(seed);
        }
        BfsFaces { graph, queue, visited }
    }
}

impl Iterator for BfsFaces<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        let fi = self.queue.pop_front()?;
        for &ei in self.graph.face_edges(fi) {
            for &adj_fi in self.graph.edge_adjacent_faces(ei) {
                if self.visited.insert(adj_fi) {
                    self.queue.push_back(adj_fi);
                }
            }
        }
        Some(fi)
    }
}

/// DFS iterator over edges reachable from a seed vertex via edge-vertex adjacency.
pub struct DfsEdgesFromVertex<'g> {
    graph: &'g BRepGraph,
    edge_stack: Vec<usize>,
    vertex_stack: Vec<usize>,
    visited_edges: HashSet<usize>,
    visited_verts: HashSet<usize>,
}

impl<'g> DfsEdgesFromVertex<'g> {
    fn new(graph: &'g BRepGraph, seed_vertex: usize) -> Self {
        let mut visited_verts = HashSet::new();
        let mut vertex_stack = Vec::new();
        if seed_vertex < graph.vertex_count {
            vertex_stack.push(seed_vertex);
            visited_verts.insert(seed_vertex);
        }
        DfsEdgesFromVertex {
            graph,
            edge_stack: Vec::new(),
            vertex_stack,
            visited_edges: HashSet::new(),
            visited_verts,
        }
    }
}

impl Iterator for DfsEdgesFromVertex<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        // Expand vertices on the vertex stack into unvisited edges.
        while let Some(vi) = self.vertex_stack.pop() {
            for &ei in self.graph.vertex_adjacent_edges(vi) {
                if self.visited_edges.insert(ei) {
                    self.edge_stack.push(ei);
                }
            }
        }
        let ei = self.edge_stack.pop()?;
        // Push far-end vertex.
        if let Some((vs, ve)) = self.graph.edge_endpoints(ei) {
            for &vn in &[vs, ve] {
                if self.visited_verts.insert(vn) {
                    self.vertex_stack.push(vn);
                }
            }
        }
        // Drain new vertices immediately.
        while let Some(vi) = self.vertex_stack.pop() {
            for &nei in self.graph.vertex_adjacent_edges(vi) {
                if self.visited_edges.insert(nei) {
                    self.edge_stack.push(nei);
                }
            }
        }
        Some(ei)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Checkpoint Data (for external checkpoint/rollback implementations)
// ─────────────────────────────────────────────────────────────────────────────

/// Public checkpoint data for external checkpoint/rollback implementations.
///
/// This is a public data structure that can be used by external crates
/// (like `rcad-advanced-modeling`) to implement checkpoint and rollback functionality.
///
/// Use [`BRepGraph::checkpoint_data`] to create and [`BRepGraph::restore_from_data`]
/// to restore.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRepGraphCheckpointData {
    pub vertex_count: usize,
    pub edge_count: usize,
    pub face_count: usize,
    pub edge_to_faces: Vec<Vec<usize>>,
    pub face_to_edges: Vec<Vec<usize>>,
    pub vertex_to_edges: Vec<Vec<usize>>,
    pub vertex_to_faces: Vec<Vec<usize>>,
    pub edge_endpoints: Vec<(usize, usize)>,
    pub vertex_dirty: Vec<bool>,
    pub edge_dirty: Vec<bool>,
    pub face_dirty: Vec<bool>,
}

impl BRepGraph {
    // ── Checkpoint data methods (for external checkpoint/rollback implementations) ──

    /// Create a public checkpoint data snapshot.
    ///
    /// This returns a [`BRepGraphCheckpointData`] that can be used by external crates
    /// to implement checkpoint and rollback functionality.
    pub fn checkpoint_data(&self) -> BRepGraphCheckpointData {
        BRepGraphCheckpointData {
            vertex_count: self.vertex_count,
            edge_count: self.edge_count,
            face_count: self.face_count,
            edge_to_faces: self.edge_to_faces.clone(),
            face_to_edges: self.face_to_edges.clone(),
            vertex_to_edges: self.vertex_to_edges.clone(),
            vertex_to_faces: self.vertex_to_faces.clone(),
            edge_endpoints: self.edge_endpoints.clone(),
            vertex_dirty: self.vertex_dirty.clone(),
            edge_dirty: self.edge_dirty.clone(),
            face_dirty: self.face_dirty.clone(),
        }
    }

    /// Restore the graph from public checkpoint data.
    ///
    /// After this call the graph's adjacency tables, counts, and dirty flags
    /// are exactly as they were when the checkpoint data was created.
    pub fn restore_from_data(&mut self, data: &BRepGraphCheckpointData) {
        self.vertex_count = data.vertex_count;
        self.edge_count = data.edge_count;
        self.face_count = data.face_count;
        self.edge_to_faces = data.edge_to_faces.clone();
        self.face_to_edges = data.face_to_edges.clone();
        self.vertex_to_edges = data.vertex_to_edges.clone();
        self.vertex_to_faces = data.vertex_to_faces.clone();
        self.edge_endpoints = data.edge_endpoints.clone();
        self.vertex_dirty = data.vertex_dirty.clone();
        self.edge_dirty = data.edge_dirty.clone();
        self.face_dirty = data.face_dirty.clone();
    }

    // ── Invariant validation ──────────────────────────────────────────────────

    /// Validate internal topology invariants and return a list of error strings.
    ///
    /// An empty list means the graph is consistent.  Checks include:
    /// - Table lengths match the declared counts.
    /// - Every (start, end) vertex index in `edge_endpoints` is in bounds.
    /// - Every face index referenced by `edge_to_faces` is in bounds.
    /// - Every edge index referenced by `face_to_edges` is in bounds.
    ///
    /// Analogous to `BRepGraph_Validate` in OCCT 8.0.
    pub fn validate_invariants(&self) -> Vec<String> {
        let mut errors: Vec<String> = Vec::new();

        if self.edge_to_faces.len() != self.edge_count {
            errors.push(format!(
                "edge_to_faces length {} ≠ edge_count {}",
                self.edge_to_faces.len(),
                self.edge_count
            ));
        }
        if self.face_to_edges.len() != self.face_count {
            errors.push(format!(
                "face_to_edges length {} ≠ face_count {}",
                self.face_to_edges.len(),
                self.face_count
            ));
        }
        if self.vertex_to_edges.len() != self.vertex_count {
            errors.push(format!(
                "vertex_to_edges length {} ≠ vertex_count {}",
                self.vertex_to_edges.len(),
                self.vertex_count
            ));
        }
        if self.vertex_to_faces.len() != self.vertex_count {
            errors.push(format!(
                "vertex_to_faces length {} ≠ vertex_count {}",
                self.vertex_to_faces.len(),
                self.vertex_count
            ));
        }
        if self.edge_endpoints.len() != self.edge_count {
            errors.push(format!(
                "edge_endpoints length {} ≠ edge_count {}",
                self.edge_endpoints.len(),
                self.edge_count
            ));
        }
        if self.vertex_dirty.len() != self.vertex_count {
            errors.push(format!(
                "vertex_dirty length {} ≠ vertex_count {}",
                self.vertex_dirty.len(),
                self.vertex_count
            ));
        }
        if self.edge_dirty.len() != self.edge_count {
            errors.push(format!(
                "edge_dirty length {} ≠ edge_count {}",
                self.edge_dirty.len(),
                self.edge_count
            ));
        }
        if self.face_dirty.len() != self.face_count {
            errors.push(format!(
                "face_dirty length {} ≠ face_count {}",
                self.face_dirty.len(),
                self.face_count
            ));
        }

        // Check edge_endpoints vertex indices are in bounds.
        for (ei, &(vs, ve)) in self.edge_endpoints.iter().enumerate() {
            if vs >= self.vertex_count {
                errors.push(format!(
                    "edge {} start vertex {} out of range (vertex_count={})",
                    ei, vs, self.vertex_count
                ));
            }
            if ve >= self.vertex_count {
                errors.push(format!(
                    "edge {} end vertex {} out of range (vertex_count={})",
                    ei, ve, self.vertex_count
                ));
            }
        }

        // Check face indices referenced by edge_to_faces are in bounds.
        for (ei, faces) in self.edge_to_faces.iter().enumerate() {
            for &fi in faces {
                if fi >= self.face_count {
                    errors.push(format!(
                        "edge {} references face {} out of range (face_count={})",
                        ei, fi, self.face_count
                    ));
                }
            }
        }

        // Check edge indices referenced by face_to_edges are in bounds.
        for (fi, edges) in self.face_to_edges.iter().enumerate() {
            for &ei in edges {
                if ei >= self.edge_count {
                    errors.push(format!(
                        "face {} references edge {} out of range (edge_count={})",
                        fi, ei, self.edge_count
                    ));
                }
            }
        }

        errors
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        BRep, Edge, Face, Shell, Solid, Surface3, Vertex, Wire, WireEdge,
        geom::PrimitiveSolid,
    };
    use glam::DVec3;

    fn unit_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        })
    }

    /// Build a minimal non-manifold BRep where edge 0 is shared by 3 faces.
    fn non_manifold_tripod() -> BRep {
        let vertices = vec![
            Vertex { point: DVec3::new(0.0, 0.0, 0.0) }, // 0
            Vertex { point: DVec3::new(1.0, 0.0, 0.0) }, // 1
            Vertex { point: DVec3::new(0.0, 1.0, 0.0) }, // 2
            Vertex { point: DVec3::new(0.0, 0.0, 1.0) }, // 3
            Vertex { point: DVec3::new(0.0, -1.0, 0.0) }, // 4
        ];
        // Edge 0 is the shared spine. Other edges are unique per triangle face.
        let edges = vec![
            Edge { start: 0, end: 1 }, // shared by 3 faces
            Edge { start: 1, end: 2 },
            Edge { start: 2, end: 0 },
            Edge { start: 1, end: 3 },
            Edge { start: 3, end: 0 },
            Edge { start: 1, end: 4 },
            Edge { start: 4, end: 0 },
        ];

        let f0 = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::new(0, true),
                    WireEdge::new(1, true),
                    WireEdge::new(2, true),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        let f1 = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::new(0, true),
                    WireEdge::new(3, true),
                    WireEdge::new(4, true),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Y,
            triangles: vec![],
            mesh_dirty: true,
        };
        let f2 = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::new(0, true),
                    WireEdge::new(5, true),
                    WireEdge::new(6, true),
                ],
            },
            inner_wires: vec![],
            normal: -DVec3::Y,
            triangles: vec![],
            mesh_dirty: true,
        };

        BRep {
            vertices,
            edges,
            solids: vec![Solid {
                shells: vec![Shell {
                    faces: vec![f0, f1, f2],
                }],
            }],
            geom: Default::default(),
            compound: None,
            compsolid: None,
        }
    }

    // ── Construction & counts ─────────────────────────────────────────────────

    #[test]
    fn box_entity_counts_correct() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        assert_eq!(g.vertex_count, 8, "box: 8 vertices");
        assert_eq!(g.edge_count, 12, "box: 12 edges");
        assert_eq!(g.face_count, 6, "box: 6 faces");
    }

    #[test]
    fn empty_brep_yields_zero_counts() {
        let brep = BRep::new();
        let g = BRepGraph::from_brep(&brep);
        assert_eq!(g.vertex_count, 0);
        assert_eq!(g.edge_count, 0);
        assert_eq!(g.face_count, 0);
        assert!(g.edge_adjacent_faces(0).is_empty());
        assert!(!g.has_modifications());
    }

    // ── O(1) adjacency ────────────────────────────────────────────────────────

    #[test]
    fn box_every_edge_has_two_adjacent_faces() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for ei in 0..g.edge_count {
            let adj = g.edge_adjacent_faces(ei);
            assert_eq!(adj.len(), 2, "edge {ei}: expected 2 adjacent faces, got {adj:?}");
        }
    }

    #[test]
    fn box_every_face_has_four_edges() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for fi in 0..g.face_count {
            let edges = g.face_edges(fi);
            assert_eq!(edges.len(), 4, "face {fi}: expected 4 edges, got {edges:?}");
        }
    }

    #[test]
    fn box_every_vertex_has_three_adjacent_edges() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for vi in 0..g.vertex_count {
            let adj = g.vertex_adjacent_edges(vi);
            assert_eq!(adj.len(), 3, "vertex {vi}: expected 3 adjacent edges, got {adj:?}");
        }
    }

    #[test]
    fn box_every_vertex_has_three_adjacent_faces() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for vi in 0..g.vertex_count {
            let adj = g.vertex_adjacent_faces(vi);
            assert_eq!(adj.len(), 3, "vertex {vi}: expected 3 adjacent faces, got {adj:?}");
        }
    }

    #[test]
    fn edge_endpoints_match_brep_edges() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for (ei, edge) in brep.edges.iter().enumerate() {
            let ep = g.edge_endpoints(ei).expect("endpoints should exist");
            assert_eq!(ep.0, edge.start);
            assert_eq!(ep.1, edge.end);
        }
    }

    // ── Manifold / closed inspection ──────────────────────────────────────────

    #[test]
    fn box_is_manifold_and_closed() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        assert!(g.is_manifold(), "unit box should be manifold");
        assert!(g.is_closed(), "unit box should be closed");
        assert!(g.non_manifold_edges().is_empty(), "no non-manifold edges");
        assert!(g.low_valence_vertices().is_empty(), "no low-valence vertices");

        let summary = g.non_manifold_summary();
        assert!(summary.is_clean());
        assert!(summary.boundary_edges.is_empty());
        assert!(summary.orphan_edges.is_empty());
        assert!(summary.multi_face_edges.is_empty());
        assert!(summary.non_manifold_vertices.is_empty());
    }

    #[test]
    fn detects_multi_face_non_manifold_edge() {
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);

        assert!(!g.is_manifold(), "tripod should be non-manifold");
        assert!(!g.is_closed(), "tripod has boundary edges on each side face");

        let multi = g.multi_face_edges();
        assert_eq!(multi, vec![0], "edge 0 should be the only multi-face edge");
        assert_eq!(g.boundary_edges(), vec![1, 2, 3, 4, 5, 6]);
        assert!(g.orphan_edges().is_empty(), "no orphan edges expected");

        let verts = g.non_manifold_vertices();
        assert_eq!(verts, vec![0, 1], "shared edge endpoints should be non-manifold vertices");

        let summary = g.non_manifold_summary();
        assert!(!summary.is_clean());
        assert_eq!(summary.multi_face_edges, vec![0]);
        assert_eq!(summary.non_manifold_vertices, vec![0, 1]);
    }

    // ── Dirty / modification tracking ─────────────────────────────────────────

    #[test]
    fn dirty_marking_and_query() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);
        assert!(!g.has_modifications());

        g.mark_vertex_modified(0);
        g.mark_edge_modified(3);
        g.mark_face_modified(5);

        assert!(g.has_modifications());
        assert_eq!(g.modified_vertices(), vec![0]);
        assert_eq!(g.modified_edges(), vec![3]);
        assert_eq!(g.modified_faces(), vec![5]);

        g.clear_dirty();
        assert!(!g.has_modifications());
        assert!(g.modified_vertices().is_empty());
        assert!(g.modified_edges().is_empty());
        assert!(g.modified_faces().is_empty());
    }

    #[test]
    fn out_of_range_mark_is_safe() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);
        // Should not panic.
        g.mark_vertex_modified(9999);
        g.mark_edge_modified(9999);
        g.mark_face_modified(9999);
        assert!(!g.has_modifications());
    }

    // ── DFS / BFS face traversal ──────────────────────────────────────────────

    #[test]
    fn bfs_faces_visits_all_box_faces() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let mut visited: Vec<usize> = g.bfs_faces(0).collect();
        visited.sort_unstable();
        assert_eq!(visited, vec![0, 1, 2, 3, 4, 5], "BFS should visit all 6 faces");
    }

    #[test]
    fn dfs_faces_visits_all_box_faces() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let mut visited: Vec<usize> = g.dfs_faces(0).collect();
        visited.sort_unstable();
        assert_eq!(visited, vec![0, 1, 2, 3, 4, 5], "DFS should visit all 6 faces");
    }

    #[test]
    fn dfs_bfs_face_traversal_no_duplicates() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let bfs: Vec<usize> = g.bfs_faces(2).collect();
        let dfs: Vec<usize> = g.dfs_faces(2).collect();
        // No duplicates.
        let bfs_set: HashSet<usize> = bfs.iter().copied().collect();
        let dfs_set: HashSet<usize> = dfs.iter().copied().collect();
        assert_eq!(bfs.len(), bfs_set.len(), "BFS duplicates");
        assert_eq!(dfs.len(), dfs_set.len(), "DFS duplicates");
    }

    #[test]
    fn out_of_range_seed_returns_empty_traversal() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        assert!(g.bfs_faces(999).next().is_none());
        assert!(g.dfs_faces(999).next().is_none());
    }

    // ── DFS edge traversal ────────────────────────────────────────────────────

    #[test]
    fn dfs_edges_from_vertex_visits_all_box_edges() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let mut visited: Vec<usize> = g.dfs_edges_from_vertex(0).collect();
        visited.sort_unstable();
        assert_eq!(visited.len(), 12, "DFS from any vertex should reach all 12 edges");
        assert_eq!(visited, (0..12).collect::<Vec<_>>());
    }

    // ── Edge valence and vertex degree ────────────────────────────────────────

    #[test]
    fn edge_valence_manifold_box() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for ei in 0..brep.edges.len() {
            assert_eq!(g.edge_valence(ei), 2, "every box edge should have valence 2");
        }
        assert_eq!(g.edge_valence(9999), 0, "out-of-range returns 0");
    }

    #[test]
    fn edge_valence_tripod_multi_face_edge() {
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);
        assert_eq!(g.edge_valence(0), 3, "shared spine should have valence 3");
        assert_eq!(g.edge_valence(1), 1, "side edge should have valence 1 (boundary)");
    }

    #[test]
    fn vertex_degree_box() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        // Every vertex of a box is incident to exactly 3 edges.
        for vi in 0..brep.vertices.len() {
            assert_eq!(g.vertex_degree(vi), 3, "box vertex {vi} should have degree 3");
        }
    }

    // ── repair_hints ──────────────────────────────────────────────────────────

    #[test]
    fn repair_hints_empty_for_manifold_box() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let hints = g.repair_hints(&brep);
        assert!(
            hints.is_empty(),
            "a closed manifold box should have no repair hints"
        );
    }

    #[test]
    fn repair_hints_detect_multi_manifold_edge_in_tripod() {
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);
        let hints = g.repair_hints(&brep);
        // Must include at least one MultiManifoldEdge hint for edge 0.
        let has_multi = hints.hints.iter().any(|h| {
            matches!(h, RepairHint::MultiManifoldEdge { edge_idx: 0, .. })
        });
        assert!(has_multi, "tripod must generate a MultiManifoldEdge hint for edge 0");
    }

    #[test]
    fn repair_hints_detect_non_manifold_vertex_in_tripod() {
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);
        let hints = g.repair_hints(&brep);
        let has_nm_vert = hints.hints.iter().any(|h| {
            matches!(h, RepairHint::NonManifoldVertex { .. })
        });
        assert!(has_nm_vert, "tripod must generate NonManifoldVertex hints");
    }

    #[test]
    fn repair_hints_detect_unmatched_boundary_in_tripod() {
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);
        let hints = g.repair_hints(&brep);
        // Tripod side edges (1-6) are all boundary edges;
        // they cannot be stitched to each other because there are no matching
        // partner edges with coincident endpoints.
        let has_unmatched = hints.hints.iter().any(|h| {
            matches!(h, RepairHint::UnmatchedBoundaryEdge { .. })
        });
        assert!(has_unmatched, "tripod must have unmatched boundary edges");
    }

    // ── Checkpoint data (public API for external checkpoint/rollback) ───────────

    #[test]
    fn checkpoint_data_captures_clean_state() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let cp = g.checkpoint_data();
        // Checkpoint counts should match the live graph.
        assert_eq!(cp.vertex_count, g.vertex_count);
        assert_eq!(cp.edge_count, g.edge_count);
        assert_eq!(cp.face_count, g.face_count);
        // Dirty flags are all false in a freshly built graph.
        assert!(cp.vertex_dirty.iter().all(|&d| !d));
        assert!(cp.edge_dirty.iter().all(|&d| !d));
        assert!(cp.face_dirty.iter().all(|&d| !d));
    }

    #[test]
    fn restore_from_data_undoes_dirty_marks() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);
        // Nothing dirty initially.
        assert!(!g.has_modifications());

        let cp = g.checkpoint_data();

        // Mark several entities dirty.
        g.mark_vertex_modified(2);
        g.mark_edge_modified(7);
        g.mark_face_modified(4);
        assert!(g.has_modifications());

        // Restore.
        g.restore_from_data(&cp);
        assert!(!g.has_modifications(), "all dirty bits should be cleared after restore");
        assert!(g.modified_vertices().is_empty());
        assert!(g.modified_edges().is_empty());
        assert!(g.modified_faces().is_empty());
    }

    // ── validate_invariants ───────────────────────────────────────────────────

    #[test]
    fn validate_invariants_passes_for_unit_box() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let errors = g.validate_invariants();
        assert!(errors.is_empty(), "unit box graph should have no invariant violations: {errors:?}");
    }

    #[test]
    fn validate_invariants_passes_for_non_manifold_tripod() {
        // The tripod is non-manifold but structurally consistent (table lengths match,
        // indices are in bounds).
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);
        let errors = g.validate_invariants();
        assert!(errors.is_empty(), "tripod graph should have no index invariant violations: {errors:?}");
    }

    // ── BRepGraphBuilder / BRepGraphTool ────────────────────────────────────

    #[test]
    fn builder_constructs_same_adjacency_as_box() {
        let brep = unit_box();
        let graph_from_brep = BRepGraph::from_brep(&brep);

        let mut builder = BRepGraph::builder(
            graph_from_brep.vertex_count,
            graph_from_brep.edge_count,
            graph_from_brep.face_count,
        );
        for edge_idx in 0..graph_from_brep.edge_count {
            let (start_vertex, end_vertex) = graph_from_brep.edge_endpoints(edge_idx).unwrap();
            builder.set_edge_endpoints(edge_idx, start_vertex, end_vertex);
            for &face_idx in graph_from_brep.edge_adjacent_faces(edge_idx) {
                builder.add_edge_face(edge_idx, face_idx);
            }
        }
        for face_idx in 0..graph_from_brep.face_count {
            for &edge_idx in graph_from_brep.face_edges(face_idx) {
                builder.add_face_edge(face_idx, edge_idx);
            }
        }
        builder.mark_face_modified(2).mark_edge_modified(1);

        let built = builder.build().expect("builder output should validate");
        assert_eq!(built.vertex_count, graph_from_brep.vertex_count);
        assert_eq!(built.edge_count, graph_from_brep.edge_count);
        assert_eq!(built.face_count, graph_from_brep.face_count);
        for edge_idx in 0..built.edge_count {
            assert_eq!(built.edge_adjacent_faces(edge_idx), graph_from_brep.edge_adjacent_faces(edge_idx));
            assert_eq!(built.edge_endpoints(edge_idx), graph_from_brep.edge_endpoints(edge_idx));
        }
        for face_idx in 0..built.face_count {
            assert_eq!(built.face_edges(face_idx), graph_from_brep.face_edges(face_idx));
        }
        assert_eq!(built.vertex_adjacent_edges(0), graph_from_brep.vertex_adjacent_edges(0));
        assert!(built.modified_faces().contains(&2));
        assert!(built.modified_edges().contains(&1));
    }

    #[test]
    fn builder_reports_invalid_indices() {
        let mut builder = BRepGraph::builder(2, 1, 1);
        builder
            .set_edge_endpoints(0, 0, 7)
            .add_edge_face(0, 0)
            .add_face_edge(0, 0);

        let errors = builder.build().expect_err("out-of-range endpoint should fail validation");
        assert!(
            errors.iter().any(|error| error.contains("end vertex 7 out of range")),
            "expected endpoint validation error, got {errors:?}"
        );
    }

    #[test]
    fn tool_exposes_flat_face_and_geometry_access() {
        let mut brep = unit_box();
        brep.geom.surfaces.push(Surface3::Plane(crate::geom::Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        }));
        brep.geom.face_surface = vec![Some(0); 6];
        brep.geom.face_surface_range = vec![Some([0.0, 1.0, 0.0, 1.0]); 6];

        let graph = BRepGraph::from_brep(&brep);
        let tool = graph.tool(&brep);

        let face = tool.face(0).expect("flat face 0 should resolve");
        assert_eq!(face.outer_wire.edges.len(), 4);
        assert_eq!(tool.face_location(0), Some((0, 0, 0)));
        assert!(tool.face_normal(0).is_some());
        assert!(matches!(tool.face_surface(0), Some(Surface3::Plane(_))));
        assert_eq!(tool.face_surface_range(0), Some([0.0, 1.0, 0.0, 1.0]));

        let (start_point, end_point) = tool.edge_points(0).expect("edge 0 should have endpoints");
        assert_eq!(start_point, brep.vertices[brep.edges[0].start].point);
        assert_eq!(end_point, brep.vertices[brep.edges[0].end].point);
    }
}
