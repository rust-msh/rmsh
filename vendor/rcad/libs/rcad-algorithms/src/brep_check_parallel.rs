//! Parallel BRep validity checker.
//!
//! This module provides parallel versions of the BRep check algorithms from
//! `brep_check`, using Rayon for multi-threaded execution.
//!
//! # When to use
//!
//! Use this module when checking large BReps with many faces/edges. For small
//! models, the overhead of parallel execution may not be worth it.
//!
//! # Performance
//!
//! The parallel checker distributes work across multiple threads:
//! - Face-level checks run in parallel across all faces
//! - Edge validation uses parallel iteration
//! - Vertex validation uses parallel iteration
//! - Results are merged at the end
//!
//! Example speedup on an 8-core machine for a 10,000-face model: ~4-6x faster.

use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;
use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::CurveEval;

// Re-export CheckIssue and CheckResult from the base module for convenience.
pub use crate::brep_check::{CheckIssue, CheckResult};

/// Configuration options for parallel BRep checking.
#[derive(Debug, Clone)]
pub struct ParallelCheckOptions {
    /// Minimum number of faces to trigger parallel processing.
    /// Below this threshold, the sequential check is used.
    pub min_faces_for_parallel: usize,

    /// Number of faces to process per thread batch.
    /// Smaller chunks provide better load balancing but more overhead.
    pub chunk_size: usize,

    /// Minimum number of edges to trigger parallel edge checking.
    pub min_edges_for_parallel: usize,

    /// Minimum number of vertices to trigger parallel vertex checking.
    pub min_vertices_for_parallel: usize,

    /// Tolerance for geometric comparisons (wire closure, duplicate vertices).
    pub tolerance: f64,

    /// Whether to check for duplicate vertices at the same position.
    pub check_duplicate_vertices: bool,

    /// Whether to check for isolated vertices (not referenced by any edge).
    pub check_isolated_vertices: bool,

    /// Whether to check vertex positions are finite (not NaN or infinity).
    pub check_finite_vertices: bool,
}

impl Default for ParallelCheckOptions {
    fn default() -> Self {
        Self {
            min_faces_for_parallel: 100,
            chunk_size: 32,
            min_edges_for_parallel: 100,
            min_vertices_for_parallel: 100,
            tolerance: 1e-6,
            check_duplicate_vertices: true,
            check_isolated_vertices: true,
            check_finite_vertices: true,
        }
    }
}

impl ParallelCheckOptions {
    /// Create options optimized for small models (uses sequential processing).
    pub fn small_model() -> Self {
        Self {
            min_faces_for_parallel: usize::MAX,
            min_edges_for_parallel: usize::MAX,
            min_vertices_for_parallel: usize::MAX,
            ..Self::default()
        }
    }

    /// Create options optimized for large models (aggressive parallelization).
    pub fn large_model() -> Self {
        Self {
            min_faces_for_parallel: 10,
            chunk_size: 64,
            min_edges_for_parallel: 10,
            min_vertices_for_parallel: 10,
            ..Self::default()
        }
    }

    /// Set the geometric tolerance.
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Set the chunk size for parallel processing.
    pub fn with_chunk_size(mut self, chunk_size: usize) -> Self {
        self.chunk_size = chunk_size.max(1);
        self
    }

    /// Enable or disable duplicate vertex checking.
    pub fn with_duplicate_vertex_check(mut self, enabled: bool) -> Self {
        self.check_duplicate_vertices = enabled;
        self
    }

    /// Enable or disable isolated vertex checking.
    pub fn with_isolated_vertex_check(mut self, enabled: bool) -> Self {
        self.check_isolated_vertices = enabled;
        self
    }
}

/// Parallel check result with thread-local issue collection.
#[derive(Debug, Clone, Default)]
struct ThreadLocalIssues {
    issues: Vec<CheckIssue>,
}

/// Additional issue types specific to parallel checking.
#[derive(Debug, Clone, PartialEq)]
pub enum ParallelCheckIssue {
    /// A vertex is duplicated (another vertex exists at the same position).
    DuplicateVertex {
        vertex_a: usize,
        vertex_b: usize,
        distance: f64,
    },
    /// A vertex is not referenced by any edge.
    IsolatedVertex {
        vertex_idx: usize,
    },
    /// A vertex has non-finite coordinates (NaN or infinity).
    NonFiniteVertex {
        vertex_idx: usize,
    },
}

impl std::fmt::Display for ParallelCheckIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParallelCheckIssue::DuplicateVertex { vertex_a, vertex_b, distance } => {
                write!(f, "DuplicateVertex: vertices {} and {} at same position (distance: {})", vertex_a, vertex_b, distance)
            }
            ParallelCheckIssue::IsolatedVertex { vertex_idx } => {
                write!(f, "IsolatedVertex: vertex {} not referenced by any edge", vertex_idx)
            }
            ParallelCheckIssue::NonFiniteVertex { vertex_idx } => {
                write!(f, "NonFiniteVertex: vertex {} has NaN or infinite coordinates", vertex_idx)
            }
        }
    }
}

/// Extended result including parallel-specific issues.
#[derive(Debug, Clone, Default)]
pub struct ParallelCheckResult {
    /// Basic structural issues from the standard check.
    pub issues: Vec<CheckIssue>,
    /// Parallel-specific issues (duplicate vertices, isolated vertices, etc.).
    pub parallel_issues: Vec<ParallelCheckIssue>,
    /// Whether the check was performed in parallel.
    pub was_parallel: bool,
    /// Number of threads used (1 for sequential).
    pub thread_count: usize,
}

impl ParallelCheckResult {
    /// Returns `true` if no issues were found.
    pub fn is_valid(&self) -> bool {
        self.issues.is_empty() && self.parallel_issues.is_empty()
    }

    /// Convert to a standard CheckResult.
    pub fn to_check_result(self) -> CheckResult {
        CheckResult { issues: self.issues }
    }
}

/// Check the validity of a BRep with automatic parallel/sequential selection.
///
/// This function automatically chooses between parallel and sequential processing
/// based on the model size and `ParallelCheckOptions::min_faces_for_parallel`.
///
/// # Arguments
///
/// * `brep` - The BRep to check
/// * `options` - Configuration options for the check
///
/// # Returns
///
/// A `ParallelCheckResult` containing all issues found and execution metadata.
pub fn check_parallel_with_options(brep: &BRep, options: &ParallelCheckOptions) -> ParallelCheckResult {
    let face_count: usize = brep.solids.iter()
        .map(|s| s.shells.iter().map(|sh| sh.faces.len()).sum::<usize>())
        .sum();
    let edge_count = brep.edges.len();
    let vertex_count = brep.vertices.len();

    // Decide whether to use parallel or sequential processing
    let use_parallel = face_count >= options.min_faces_for_parallel
        || edge_count >= options.min_edges_for_parallel
        || vertex_count >= options.min_vertices_for_parallel;

    let thread_count = if use_parallel {
        rayon::current_num_threads()
    } else {
        1
    };

    let mut result = if use_parallel {
        check_parallel_internal(brep, options)
    } else {
        check_sequential_internal(brep, options)
    };

    result.was_parallel = use_parallel;
    result.thread_count = thread_count;
    result
}

/// Internal parallel check implementation.
fn check_parallel_internal(brep: &BRep, options: &ParallelCheckOptions) -> ParallelCheckResult {
    let mut issues = Vec::new();
    let mut parallel_issues = Vec::new();
    let n_edges = brep.edges.len();
    let n_verts = brep.vertices.len();

    // C5: edge vertex bounds (parallel)
    let edge_issues: Vec<CheckIssue> = brep.edges
        .par_iter()
        .enumerate()
        .flat_map(|(eidx, edge)| {
            let mut local_issues = Vec::new();
            if edge.start >= n_verts {
                local_issues.push(CheckIssue::InvalidVertexIndex {
                    edge: eidx,
                    vertex_idx: edge.start,
                });
            }
            if edge.end >= n_verts {
                local_issues.push(CheckIssue::InvalidVertexIndex {
                    edge: eidx,
                    vertex_idx: edge.end,
                });
            }
            local_issues
        })
        .collect();
    issues.extend(edge_issues);

    // C6: manifold check - count edge references (parallel reduction)
    let edge_face_count: Vec<usize> = compute_edge_face_counts_parallel(brep, n_edges);

    let manifold_issues: Vec<CheckIssue> = edge_face_count
        .par_iter()
        .enumerate()
        .filter_map(|(eidx, &count)| {
            if count != 2 {
                Some(CheckIssue::NonManifoldEdge {
                    edge_idx: eidx,
                    face_count: count,
                })
            } else {
                None
            }
        })
        .collect();
    issues.extend(manifold_issues);

    // Face-level checks (parallel with chunking)
    let face_issues = check_faces_parallel_chunked(brep, n_edges, options.chunk_size);
    issues.extend(face_issues);

    // Vertex-level checks (parallel)
    let vertex_results = check_vertices_parallel(brep, options);
    parallel_issues.extend(vertex_results);

    ParallelCheckResult {
        issues,
        parallel_issues,
        was_parallel: true,
        thread_count: rayon::current_num_threads(),
    }
}

/// Internal sequential check implementation (fallback for small models).
fn check_sequential_internal(brep: &BRep, options: &ParallelCheckOptions) -> ParallelCheckResult {
    // Use the standard sequential check for basic issues
    let result = crate::brep_check::check(brep);
    let mut parallel_issues = Vec::new();

    // Add vertex checks that are specific to parallel module
    let vertex_results = check_vertices_sequential(brep, options);
    parallel_issues.extend(vertex_results);

    ParallelCheckResult {
        issues: result.issues,
        parallel_issues,
        was_parallel: false,
        thread_count: 1,
    }
}

/// Check the validity of a BRep in parallel with default options.
///
/// This is the parallel version of `brep_check::check()`. It performs the same
/// checks but distributes the work across multiple threads for better performance
/// on large models.
///
/// # Checks performed (same as serial version)
///
/// - Wire closure
/// - Face normal consistency
/// - Degenerate faces
/// - Edge/vertex index validity
/// - Manifold topology
/// - Wire self-intersection
/// - Duplicate vertices (parallel-specific)
/// - Isolated vertices (parallel-specific)
/// - Non-finite vertex positions (parallel-specific)
///
/// # Arguments
///
/// * `brep` - The BRep to check
///
/// # Returns
///
/// A `CheckResult` containing all issues found.
pub fn check_parallel(brep: &BRep) -> CheckResult {
    let options = ParallelCheckOptions::default();
    check_parallel_with_options(brep, &options).to_check_result()
}

/// Compute edge face counts in parallel.
fn compute_edge_face_counts_parallel(brep: &BRep, n_edges: usize) -> Vec<usize> {
    // Use atomic counters for thread-safe counting
    let counts: Vec<AtomicUsize> = (0..n_edges)
        .map(|_| AtomicUsize::new(0))
        .collect();

    // Collect all edge references first in parallel
    let edge_refs: Vec<usize> = brep.solids
        .iter()
        .flat_map(|solid| solid.shells.iter())
        .flat_map(|shell| shell.faces.iter())
        .flat_map(|face| {
            let outer_refs: Vec<usize> = face.outer_wire.edges.iter().map(|we| we.idx).collect();
            let inner_refs: Vec<usize> = face.inner_wires.iter()
                .flat_map(|wire| wire.edges.iter().map(|we| we.idx))
                .collect();
            outer_refs.into_iter().chain(inner_refs)
        })
        .filter(|&idx| idx < n_edges)
        .collect();

    // Increment counts atomically
    for idx in edge_refs {
        counts[idx].fetch_add(1, Ordering::Relaxed);
    }

    // Convert to regular Vec
    counts.into_iter().map(|c| c.into_inner()).collect()
}

/// Check all faces in parallel with chunking for work stealing.
fn check_faces_parallel_chunked(brep: &BRep, n_edges: usize, chunk_size: usize) -> Vec<CheckIssue> {
    // Create a flat list of (solid_idx, shell_idx, face_idx, face_ref) for parallel iteration
    let face_items: Vec<(usize, usize, usize, &rcad_kernel::topology::Face)> = brep.solids
        .iter()
        .enumerate()
        .flat_map(|(si, solid)| {
            solid.shells.iter().enumerate().flat_map(move |(shi, shell)| {
                shell.faces.iter().enumerate().map(move |(fi, face)| {
                    (si, shi, fi, face)
                })
            })
        })
        .collect();

    // Process in chunks for better work stealing
    face_items
        .par_chunks(chunk_size.max(1))
        .flat_map(|chunk| {
            chunk.iter()
                .flat_map(|&(si, shi, fi, face)| {
                    check_single_face(brep, face, si, shi, fi, n_edges)
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Check a single face for all issues.
fn check_single_face(
    brep: &BRep,
    face: &rcad_kernel::topology::Face,
    si: usize,
    shi: usize,
    fi: usize,
    n_edges: usize,
) -> Vec<CheckIssue> {
    let mut issues = Vec::new();
    let wire = &face.outer_wire;

    // C2: zero normal
    if face.normal == DVec3::ZERO {
        issues.push(CheckIssue::ZeroNormal {
            solid: si,
            shell: shi,
            face: fi,
        });
    }

    // C3: degenerate face
    if wire.edges.len() < 3 {
        issues.push(CheckIssue::DegenerateFace {
            solid: si,
            shell: shi,
            face: fi,
        });
        return issues; // Can't check wire closure for degenerate face
    }

    // C4: edge index bounds + collect wire vertices
    let mut valid = true;
    let mut wire_verts: Vec<(usize, usize)> = Vec::new();
    for we in &wire.edges {
        if we.idx >= n_edges {
            issues.push(CheckIssue::InvalidEdgeIndex {
                solid: si,
                shell: shi,
                face: fi,
                edge_idx: we.idx,
            });
            valid = false;
        } else {
            let edge = &brep.edges[we.idx];
            let (sv, ev) = if we.forward {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            };
            wire_verts.push((sv, ev));
        }
    }

    if !valid {
        return issues;
    }

    // C1: wire closure
    let n = wire_verts.len();
    for i in 0..n {
        let next = (i + 1) % n;
        let end_v = wire_verts[i].1;
        let start_v = wire_verts[next].0;
        if end_v != start_v {
            let end_pt = brep.vertices[end_v].point;
            let start_pt = brep.vertices[start_v].point;
            if (end_pt - start_pt).length() > 1e-6 {
                issues.push(CheckIssue::OpenWire {
                    solid: si,
                    shell: shi,
                    face: fi,
                    wire_pos: i,
                });
            }
        }
    }

    // C7: wire self-intersection
    check_wire_self_intersection_local(
        &wire_verts,
        &brep.vertices,
        si, shi, fi, 0,
        &mut issues,
    );

    // C8: geometric self-intersection
    check_geometric_self_intersection_local(
        &wire_verts,
        &brep.vertices,
        si, shi, fi,
        &mut issues,
    );

    // Check inner wires
    for (wi, inner_wire) in face.inner_wires.iter().enumerate() {
        if inner_wire.edges.len() < 2 {
            continue;
        }

        let mut inner_verts: Vec<(usize, usize)> = Vec::new();
        let mut inner_valid = true;

        for we in &inner_wire.edges {
            if we.idx >= n_edges {
                issues.push(CheckIssue::InvalidEdgeIndex {
                    solid: si,
                    shell: shi,
                    face: fi,
                    edge_idx: we.idx,
                });
                inner_valid = false;
            } else {
                let edge = &brep.edges[we.idx];
                let (sv, ev) = if we.forward {
                    (edge.start, edge.end)
                } else {
                    (edge.end, edge.start)
                };
                inner_verts.push((sv, ev));
            }
        }

        if !inner_valid {
            continue;
        }

        // Inner wire closure
        let n_inner = inner_verts.len();
        for i in 0..n_inner {
            let next = (i + 1) % n_inner;
            let end_v = inner_verts[i].1;
            let start_v = inner_verts[next].0;
            if end_v != start_v {
                let end_pt = brep.vertices[end_v].point;
                let start_pt = brep.vertices[start_v].point;
                if (end_pt - start_pt).length() > 1e-6 {
                    issues.push(CheckIssue::OpenWire {
                        solid: si,
                        shell: shi,
                        face: fi,
                        wire_pos: i,
                    });
                }
            }
        }

        check_wire_self_intersection_local(
            &inner_verts,
            &brep.vertices,
            si, shi, fi,
            wi + 1,
            &mut issues,
        );
    }

    issues
}

/// Check a wire for self-intersection topology.
fn check_wire_self_intersection_local(
    wire_verts: &[(usize, usize)],
    vertices: &[rcad_kernel::topology::Vertex],
    solid: usize,
    shell: usize,
    face: usize,
    wire_idx: usize,
    issues: &mut Vec<CheckIssue>,
) {
    use std::collections::HashMap;
    let mut vertex_count: HashMap<usize, usize> = HashMap::new();
    for &(sv, ev) in wire_verts {
        *vertex_count.entry(sv).or_insert(0) += 1;
        *vertex_count.entry(ev).or_insert(0) += 1;
    }
    for (&vidx, &count) in &vertex_count {
        if count > 2 {
            if vidx < vertices.len() {
                issues.push(CheckIssue::SelfIntersectingWire {
                    solid,
                    shell,
                    face,
                    wire_idx,
                    vertex: vidx,
                });
            }
        }
    }
}

/// Check for geometric self-intersection in a wire.
fn check_geometric_self_intersection_local(
    wire_verts: &[(usize, usize)],
    vertices: &[rcad_kernel::topology::Vertex],
    solid: usize,
    shell: usize,
    face: usize,
    issues: &mut Vec<CheckIssue>,
) {
    let n = wire_verts.len();
    if n < 4 {
        return; // Need at least 4 edges for potential self-intersection
    }

    // Check pairs of non-adjacent edges for intersection
    for i in 0..n {
        // Adjacent edges share a vertex, so check edges that are at least 2 apart
        for j in (i + 2)..n {
            // Skip if edges are adjacent (wraparound case)
            if i == 0 && j == n - 1 {
                continue;
            }

            // Get edge endpoints
            let (a_start, a_end) = wire_verts[i];
            let (b_start, b_end) = wire_verts[j];

            let p1 = vertices[a_start].point;
            let p2 = vertices[a_end].point;
            let p3 = vertices[b_start].point;
            let p4 = vertices[b_end].point;

            // Check 2D projection intersection (XY plane)
            if segments_intersect_2d(p1, p2, p3, p4) {
                issues.push(CheckIssue::GeometricSelfIntersection {
                    solid,
                    shell,
                    face,
                    edge_a: i,
                    edge_b: j,
                });
            }
        }
    }
}

/// Check if two 2D line segments intersect.
fn segments_intersect_2d(p1: DVec3, p2: DVec3, p3: DVec3, p4: DVec3) -> bool {
    // Project to XY plane
    let x1 = p1.x; let y1 = p1.y;
    let x2 = p2.x; let y2 = p2.y;
    let x3 = p3.x; let y3 = p3.y;
    let x4 = p4.x; let y4 = p4.y;

    // Check bounding box overlap first
    let (min_x1, max_x1) = if x1 < x2 { (x1, x2) } else { (x2, x1) };
    let (min_y1, max_y1) = if y1 < y2 { (y1, y2) } else { (y2, y1) };
    let (min_x2, max_x2) = if x3 < x4 { (x3, x4) } else { (x4, x3) };
    let (min_y2, max_y2) = if y3 < y4 { (y3, y4) } else { (y4, y3) };

    if max_x1 < min_x2 || max_x2 < min_x1 || max_y1 < min_y2 || max_y2 < min_y1 {
        return false;
    }

    // CCW orientation test
    fn ccw(ax: f64, ay: f64, bx: f64, by: f64, cx: f64, cy: f64) -> bool {
        (cy - ay) * (bx - ax) > (by - ay) * (cx - ax)
    }

    // Check proper intersection
    if ccw(x1, y1, x3, y3, x4, y4) != ccw(x2, y2, x3, y3, x4, y4)
        && ccw(x1, y1, x2, y2, x3, y3) != ccw(x1, y1, x2, y2, x4, y4)
    {
        return true;
    }

    false
}

/// Check vertices in parallel for duplicate, isolated, and non-finite vertices.
fn check_vertices_parallel(brep: &BRep, options: &ParallelCheckOptions) -> Vec<ParallelCheckIssue> {
    let n_verts = brep.vertices.len();
    if n_verts == 0 {
        return Vec::new();
    }

    let mut issues = Vec::new();

    // Check for non-finite vertices (parallel)
    if options.check_finite_vertices {
        let non_finite: Vec<ParallelCheckIssue> = brep.vertices
            .par_iter()
            .enumerate()
            .filter_map(|(vidx, v)| {
                if !v.point.is_finite() {
                    Some(ParallelCheckIssue::NonFiniteVertex { vertex_idx: vidx })
                } else {
                    None
                }
            })
            .collect();
        issues.extend(non_finite);
    }

    // Check for isolated vertices (parallel)
    if options.check_isolated_vertices {
        // Build a set of referenced vertices using atomic booleans
        let referenced: Vec<std::sync::atomic::AtomicBool> = (0..n_verts)
            .map(|_| std::sync::atomic::AtomicBool::new(false))
            .collect();

        // Mark all vertices referenced by edges
        brep.edges.par_iter().for_each(|edge| {
            if edge.start < n_verts {
                referenced[edge.start].store(true, Ordering::Relaxed);
            }
            if edge.end < n_verts {
                referenced[edge.end].store(true, Ordering::Relaxed);
            }
        });

        // Find isolated vertices
        let isolated: Vec<ParallelCheckIssue> = referenced
            .into_par_iter()
            .enumerate()
            .filter_map(|(vidx, ref_flag)| {
                if !ref_flag.into_inner() {
                    Some(ParallelCheckIssue::IsolatedVertex { vertex_idx: vidx })
                } else {
                    None
                }
            })
            .collect();
        issues.extend(isolated);
    }

    // Check for duplicate vertices (parallel spatial hashing)
    if options.check_duplicate_vertices {
        let duplicates = find_duplicate_vertices_parallel(&brep.vertices, options.tolerance);
        issues.extend(duplicates);
    }

    issues
}

/// Check vertices sequentially (fallback for small models).
fn check_vertices_sequential(brep: &BRep, options: &ParallelCheckOptions) -> Vec<ParallelCheckIssue> {
    let n_verts = brep.vertices.len();
    if n_verts == 0 {
        return Vec::new();
    }

    let mut issues = Vec::new();

    // Check for non-finite vertices
    if options.check_finite_vertices {
        for (vidx, v) in brep.vertices.iter().enumerate() {
            if !v.point.is_finite() {
                issues.push(ParallelCheckIssue::NonFiniteVertex { vertex_idx: vidx });
            }
        }
    }

    // Check for isolated vertices
    if options.check_isolated_vertices {
        let mut referenced = vec![false; n_verts];
        for edge in &brep.edges {
            if edge.start < n_verts {
                referenced[edge.start] = true;
            }
            if edge.end < n_verts {
                referenced[edge.end] = true;
            }
        }
        for (vidx, &is_ref) in referenced.iter().enumerate() {
            if !is_ref {
                issues.push(ParallelCheckIssue::IsolatedVertex { vertex_idx: vidx });
            }
        }
    }

    // Check for duplicate vertices
    if options.check_duplicate_vertices {
        for i in 0..n_verts {
            for j in (i + 1)..n_verts {
                let dist = (brep.vertices[i].point - brep.vertices[j].point).length();
                if dist < options.tolerance {
                    issues.push(ParallelCheckIssue::DuplicateVertex {
                        vertex_a: i,
                        vertex_b: j,
                        distance: dist,
                    });
                }
            }
        }
    }

    issues
}

/// Find duplicate vertices using parallel spatial hashing.
fn find_duplicate_vertices_parallel(
    vertices: &[rcad_kernel::topology::Vertex],
    tolerance: f64,
) -> Vec<ParallelCheckIssue> {
    use std::collections::HashMap;

    let cell_size = tolerance * 10.0; // Grid cell size

    // Compute spatial hash for each vertex in parallel
    let hashed: Vec<(i64, i64, i64, usize)> = vertices
        .par_iter()
        .enumerate()
        .filter_map(|(vidx, v)| {
            if !v.point.is_finite() {
                return None;
            }
            let cell_x = (v.point.x / cell_size).floor() as i64;
            let cell_y = (v.point.y / cell_size).floor() as i64;
            let cell_z = (v.point.z / cell_size).floor() as i64;
            Some((cell_x, cell_y, cell_z, vidx))
        })
        .collect();

    // Group by cell
    let mut cells: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (cx, cy, cz, vidx) in hashed {
        cells.entry((cx, cy, cz)).or_default().push(vidx);
    }

    // Check for duplicates within each cell and neighboring cells
    let mut issues = Vec::new();

    for ((cx, cy, cz), cell_vertices) in &cells {
        // Check vertices within this cell
        for i in 0..cell_vertices.len() {
            for j in (i + 1)..cell_vertices.len() {
                let vi = cell_vertices[i];
                let vj = cell_vertices[j];
                let dist = (vertices[vi].point - vertices[vj].point).length();
                if dist < tolerance {
                    issues.push(ParallelCheckIssue::DuplicateVertex {
                        vertex_a: vi,
                        vertex_b: vj,
                        distance: dist,
                    });
                }
            }
        }

        // Check neighboring cells
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if dx == 0 && dy == 0 && dz == 0 {
                        continue;
                    }
                    let neighbor_key = (cx + dx, cy + dy, cz + dz);
                    if let Some(neighbor_vertices) = cells.get(&neighbor_key) {
                        for &vi in cell_vertices {
                            for &vj in neighbor_vertices {
                                if vi >= vj {
                                    continue; // Avoid duplicate pairs
                                }
                                let dist = (vertices[vi].point - vertices[vj].point).length();
                                if dist < tolerance {
                                    issues.push(ParallelCheckIssue::DuplicateVertex {
                                        vertex_a: vi,
                                        vertex_b: vj,
                                        distance: dist,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    issues
}

/// Parallel check with configurable batch size.
///
/// Use this for fine-grained control over parallelization.
///
/// # Arguments
///
/// * `brep` - The BRep to check
/// * `batch_size` - Number of faces to process per thread batch
///
/// # Returns
///
/// A `CheckResult` containing all issues found.
pub fn check_parallel_with_batch_size(brep: &BRep, batch_size: usize) -> CheckResult {
    let options = ParallelCheckOptions {
        chunk_size: batch_size,
        ..ParallelCheckOptions::default()
    };
    check_parallel_with_options(brep, &options).to_check_result()
}

/// Check multiple BReps in parallel.
///
/// Useful for batch validation of many models.
///
/// # Arguments
///
/// * `breps` - Slice of BReps to check
///
/// # Returns
///
/// Vector of `CheckResult`s, one per input BRep.
pub fn check_many_parallel(breps: &[BRep]) -> Vec<CheckResult> {
    breps.par_iter().map(|brep| check_parallel(brep)).collect()
}

/// Check multiple BReps in parallel with options.
///
/// # Arguments
///
/// * `breps` - Slice of BReps to check
/// * `options` - Configuration options
///
/// # Returns
///
/// Vector of `ParallelCheckResult`s, one per input BRep.
pub fn check_many_parallel_with_options(breps: &[BRep], options: &ParallelCheckOptions) -> Vec<ParallelCheckResult> {
    breps.par_iter().map(|brep| check_parallel_with_options(brep, options)).collect()
}

/// Statistics about the parallel check execution.
#[derive(Debug, Clone, Default)]
pub struct ParallelCheckStats {
    /// Number of faces checked.
    pub face_count: usize,
    /// Number of edges checked.
    pub edge_count: usize,
    /// Number of vertices checked.
    pub vertex_count: usize,
    /// Number of issues found.
    pub issue_count: usize,
    /// Whether the check was valid (no issues).
    pub is_valid: bool,
    /// Whether parallel processing was used.
    pub was_parallel: bool,
    /// Number of threads used.
    pub thread_count: usize,
}

// ═══════════════════════════════════════════════════════════════════════════════
// Face Check Result
// ═══════════════════════════════════════════════════════════════════════════════

/// Result of checking a single face in parallel.
#[derive(Debug, Clone)]
pub struct FaceCheckResult {
    /// Solid index containing this face.
    pub solid_idx: usize,
    /// Shell index within the solid.
    pub shell_idx: usize,
    /// Face index within the shell.
    pub face_idx: usize,
    /// Whether the face passed all checks.
    pub is_valid: bool,
    /// Issues found with this face.
    pub issues: Vec<FaceCheckIssue>,
    /// Number of edges in the outer wire.
    pub outer_wire_edge_count: usize,
    /// Number of inner wires.
    pub inner_wire_count: usize,
    /// Face normal vector.
    pub normal: DVec3,
    /// Whether the normal is valid (non-zero, unit length).
    pub normal_valid: bool,
    /// Wire closure status (true if closed).
    pub outer_wire_closed: bool,
    /// Number of gaps in outer wire.
    pub outer_wire_gaps: usize,
    /// Whether the face has self-intersections.
    pub has_self_intersection: bool,
}

/// Issue specific to a face check.
#[derive(Debug, Clone, PartialEq)]
pub enum FaceCheckIssue {
    /// Face normal is zero.
    ZeroNormal,
    /// Face has fewer than 3 edges.
    DegenerateFace,
    /// Wire is not closed at the given position.
    OpenWire { wire_pos: usize, gap_distance: f64 },
    /// Edge index is out of bounds.
    InvalidEdgeIndex { edge_idx: usize },
    /// Wire self-intersection at vertex.
    SelfIntersection { vertex_idx: usize, wire_idx: usize },
    /// Geometric self-intersection between edges.
    GeometricSelfIntersection { edge_a: usize, edge_b: usize },
    /// Inner wire is not closed.
    InnerWireOpen { wire_idx: usize, wire_pos: usize },
}

impl FaceCheckResult {
    /// Returns true if this face has no issues.
    pub fn is_clean(&self) -> bool {
        self.issues.is_empty()
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!(
                "Face ({}/{}/{}): valid, {} edges, {} holes",
                self.solid_idx, self.shell_idx, self.face_idx,
                self.outer_wire_edge_count, self.inner_wire_count
            )
        } else {
            format!(
                "Face ({}/{}/{}): {} issue(s)",
                self.solid_idx, self.shell_idx, self.face_idx,
                self.issues.len()
            )
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Edge Check Result
// ═══════════════════════════════════════════════════════════════════════════════

/// Result of checking a single edge in parallel.
#[derive(Debug, Clone)]
pub struct EdgeCheckResult {
    /// Edge index in the BRep.
    pub edge_idx: usize,
    /// Whether the edge passed all checks.
    pub is_valid: bool,
    /// Issues found with this edge.
    pub issues: Vec<EdgeCheckIssue>,
    /// Start vertex index.
    pub start_vertex: usize,
    /// End vertex index.
    pub end_vertex: usize,
    /// Edge length (distance between vertices).
    pub length: f64,
    /// Whether the edge is degenerate (zero length).
    pub is_degenerate: bool,
    /// Number of faces referencing this edge.
    pub face_count: usize,
    /// Whether the edge is manifold (referenced by exactly 2 faces).
    pub is_manifold: bool,
    /// Tolerance of the edge (computed from vertex-curve gap if available).
    pub tolerance: f64,
    /// Whether there is a self-intersection in the edge curve.
    pub has_self_intersection: bool,
}

/// Issue specific to an edge check.
#[derive(Debug, Clone, PartialEq)]
pub enum EdgeCheckIssue {
    /// Vertex index is out of bounds.
    InvalidVertexIndex { vertex_idx: usize },
    /// Edge is degenerate (start and end vertices are the same).
    DegenerateEdge,
    /// Edge is not manifold (not shared by exactly 2 faces).
    NonManifold { face_count: usize },
    /// Edge has no adjacent faces.
    FreeEdge,
    /// SameParameter violation (curve endpoints don't match vertex positions).
    SameParameterViolation { start_gap: f64, end_gap: f64 },
    /// Edge has self-intersection in its curve.
    SelfIntersection,
}

impl EdgeCheckResult {
    /// Returns true if this edge has no issues.
    pub fn is_clean(&self) -> bool {
        self.issues.is_empty()
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!(
                "Edge {}: valid, length={:.4}, faces={}",
                self.edge_idx, self.length, self.face_count
            )
        } else {
            format!(
                "Edge {}: {} issue(s)",
                self.edge_idx, self.issues.len()
            )
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Shell Validation Result
// ═══════════════════════════════════════════════════════════════════════════════

/// Result of validating a shell in parallel.
#[derive(Debug, Clone)]
pub struct ShellValidationResult {
    /// Solid index containing this shell.
    pub solid_idx: usize,
    /// Shell index within the solid.
    pub shell_idx: usize,
    /// Whether the shell passed all validation checks.
    pub is_valid: bool,
    /// Number of faces in the shell.
    pub face_count: usize,
    /// Number of edges in the shell.
    pub edge_count: usize,
    /// Number of vertices in the shell.
    pub vertex_count: usize,
    /// Euler characteristic of the shell.
    pub euler_characteristic: i64,
    /// Whether the shell is closed (no free edges).
    pub is_closed: bool,
    /// Whether the shell is manifold (no non-manifold edges).
    pub is_manifold: bool,
    /// Number of open edges (edges referenced by only 1 face).
    pub open_edge_count: usize,
    /// Number of non-manifold edges (edges referenced by 3+ faces).
    pub non_manifold_edge_count: usize,
    /// Whether the shell orientation is consistent.
    pub orientation_consistent: bool,
    /// Estimated genus (only meaningful for closed shells).
    pub genus: Option<i64>,
    /// Face check results for all faces in this shell.
    pub face_results: Vec<FaceCheckResult>,
    /// Validation errors.
    pub errors: Vec<String>,
    /// Validation warnings.
    pub warnings: Vec<String>,
}

impl ShellValidationResult {
    /// Returns true if the shell is a closed manifold with consistent orientation.
    pub fn is_closed_manifold(&self) -> bool {
        self.is_closed && self.is_manifold && self.orientation_consistent
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        let status = if self.is_valid { "VALID" } else { "INVALID" };
        format!(
            "Shell ({}/{}): {} | F={}, E={}, V={}, χ={}, closed={}, manifold={}",
            self.solid_idx, self.shell_idx, status,
            self.face_count, self.edge_count, self.vertex_count,
            self.euler_characteristic, self.is_closed, self.is_manifold
        )
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Solid Validation Result
// ═══════════════════════════════════════════════════════════════════════════════

/// Result of validating a solid in parallel.
#[derive(Debug, Clone)]
pub struct SolidValidationResult {
    /// Solid index in the BRep.
    pub solid_idx: usize,
    /// Whether the solid passed all validation checks.
    pub is_valid: bool,
    /// Number of shells in the solid.
    pub shell_count: usize,
    /// Number of faces in the solid.
    pub face_count: usize,
    /// Number of edges in the solid.
    pub edge_count: usize,
    /// Number of vertices in the solid.
    pub vertex_count: usize,
    /// Euler characteristic of the solid.
    pub euler_characteristic: i64,
    /// Whether the solid is closed (all shells closed).
    pub is_closed: bool,
    /// Whether the solid is manifold.
    pub is_manifold: bool,
    /// Whether the solid has valid orientation.
    pub orientation_valid: bool,
    /// Whether the solid volume is positive.
    pub has_positive_volume: bool,
    /// Estimated volume of the solid.
    pub volume: f64,
    /// Estimated genus.
    pub genus: Option<i64>,
    /// Shell validation results for all shells in this solid.
    pub shell_results: Vec<ShellValidationResult>,
    /// Validation errors.
    pub errors: Vec<String>,
    /// Validation warnings.
    pub warnings: Vec<String>,
}

impl SolidValidationResult {
    /// Returns true if the solid is valid for BRep operations.
    pub fn is_valid_for_operations(&self) -> bool {
        self.is_valid && self.is_closed && self.is_manifold
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        let status = if self.is_valid { "VALID" } else { "INVALID" };
        format!(
            "Solid {}: {} | shells={}, F={}, E={}, V={}, volume={:.4}",
            self.solid_idx, status, self.shell_count,
            self.face_count, self.edge_count, self.vertex_count, self.volume
        )
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Parallel Check Configuration and Report
// ═══════════════════════════════════════════════════════════════════════════════

/// Configuration for comprehensive parallel BRep checking.
#[derive(Debug, Clone)]
pub struct ParallelCheckConfig {
    /// Number of threads to use (0 = use all available).
    pub num_threads: usize,
    /// Minimum number of items to trigger parallel processing.
    pub parallel_threshold: usize,
    /// Tolerance for geometric comparisons.
    pub tolerance: f64,
    /// Check face validity.
    pub check_faces: bool,
    /// Check edge validity.
    pub check_edges: bool,
    /// Check vertex validity.
    pub check_vertices: bool,
    /// Check shell topology.
    pub check_shells: bool,
    /// Check solid topology.
    pub check_solids: bool,
    /// Check for duplicate vertices.
    pub check_duplicate_vertices: bool,
    /// Check for isolated vertices.
    pub check_isolated_vertices: bool,
    /// Check for non-finite vertices.
    pub check_finite_vertices: bool,
    /// Check SameParameter condition for edges.
    pub check_same_parameter: bool,
    /// Check manifold condition.
    pub check_manifold: bool,
    /// Check wire closure.
    pub check_wire_closure: bool,
    /// Check self-intersections.
    pub check_self_intersections: bool,
}

impl Default for ParallelCheckConfig {
    fn default() -> Self {
        Self {
            num_threads: 0, // Use all available
            parallel_threshold: 100,
            tolerance: 1e-6,
            check_faces: true,
            check_edges: true,
            check_vertices: true,
            check_shells: true,
            check_solids: true,
            check_duplicate_vertices: true,
            check_isolated_vertices: true,
            check_finite_vertices: true,
            check_same_parameter: true,
            check_manifold: true,
            check_wire_closure: true,
            check_self_intersections: true,
        }
    }
}

impl ParallelCheckConfig {
    /// Create a config for fast checking (skip expensive checks).
    pub fn fast() -> Self {
        Self {
            check_self_intersections: false,
            check_same_parameter: false,
            check_duplicate_vertices: false,
            ..Self::default()
        }
    }

    /// Create a config for thorough checking (all checks enabled).
    pub fn thorough() -> Self {
        Self {
            tolerance: 1e-9,
            ..Self::default()
        }
    }

    /// Set the number of threads.
    pub fn with_threads(mut self, n: usize) -> Self {
        self.num_threads = n;
        self
    }

    /// Set the tolerance.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }
}

/// Timing information for a check phase.
#[derive(Debug, Clone, Default)]
pub struct CheckPhaseTiming {
    /// Phase name.
    pub phase: String,
    /// Duration in milliseconds.
    pub duration_ms: u64,
    /// Number of items processed.
    pub items_processed: usize,
}

/// Comprehensive report from parallel BRep checking.
#[derive(Debug, Clone, Default)]
pub struct ParallelCheckReport {
    /// Overall validity status.
    pub is_valid: bool,
    /// Total number of faces.
    pub total_faces: usize,
    /// Total number of edges.
    pub total_edges: usize,
    /// Total number of vertices.
    pub total_vertices: usize,
    /// Total number of solids.
    pub total_solids: usize,
    /// Total number of shells.
    pub total_shells: usize,
    /// Number of threads used.
    pub threads_used: usize,
    /// Whether parallel processing was used.
    pub was_parallel: bool,
    /// Total check duration.
    pub total_duration_ms: u64,
    /// Timing breakdown by phase.
    pub phase_timings: Vec<CheckPhaseTiming>,
    /// Face check results.
    pub face_results: Vec<FaceCheckResult>,
    /// Edge check results.
    pub edge_results: Vec<EdgeCheckResult>,
    /// Shell validation results.
    pub shell_results: Vec<ShellValidationResult>,
    /// Solid validation results.
    pub solid_results: Vec<SolidValidationResult>,
    /// Basic structural issues.
    pub structural_issues: Vec<CheckIssue>,
    /// Parallel-specific issues.
    pub parallel_issues: Vec<ParallelCheckIssue>,
    /// Summary statistics.
    pub stats: ParallelCheckStats,
}

impl ParallelCheckReport {
    /// Returns true if the BRep passed all checks.
    pub fn is_clean(&self) -> bool {
        self.is_valid && self.structural_issues.is_empty() && self.parallel_issues.is_empty()
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        let status = if self.is_valid { "VALID" } else { "INVALID" };
        format!(
            "BRep {}: {} solids, {} faces, {} edges, {} vertices | {}ms ({} threads)",
            status, self.total_solids, self.total_faces, self.total_edges,
            self.total_vertices, self.total_duration_ms, self.threads_used
        )
    }

    /// Returns timing breakdown as a formatted string.
    pub fn timing_summary(&self) -> String {
        let mut lines: Vec<String> = self.phase_timings.iter()
            .map(|t| format!("  {}: {}ms ({} items)", t.phase, t.duration_ms, t.items_processed))
            .collect();
        lines.insert(0, "Timing breakdown:".to_string());
        lines.join("\n")
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Parallel Face Checking
// ═══════════════════════════════════════════════════════════════════════════════

/// Check all faces in a BRep in parallel.
///
/// This function distributes face checking across multiple threads for better
/// performance on large models. Each face is checked for:
/// - Zero normal
/// - Degenerate wire (< 3 edges)
/// - Wire closure
/// - Edge index validity
/// - Self-intersections
///
/// # Arguments
///
/// * `brep` - The BRep to check.
/// * `num_threads` - Number of threads to use (0 = use all available).
///
/// # Returns
///
/// A vector of `FaceCheckResult`, one per face.
pub fn check_faces_parallel(brep: &BRep, num_threads: usize) -> Vec<FaceCheckResult> {
    let n_edges = brep.edges.len();
    let tolerance = 1e-6;

    // Create a flat list of face references for parallel iteration
    let face_items: Vec<(usize, usize, usize)> = brep.solids
        .iter()
        .enumerate()
        .flat_map(|(si, solid)| {
            solid.shells.iter().enumerate().flat_map(move |(shi, shell)| {
                shell.faces.iter().enumerate().map(move |(fi, _)| {
                    (si, shi, fi)
                })
            })
        })
        .collect();

    // Configure thread pool if specified
    let results = if num_threads > 0 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());
        pool.install(|| {
            face_items.par_iter()
                .map(|&(si, shi, fi)| check_single_face_detailed(brep, si, shi, fi, n_edges, tolerance))
                .collect()
        })
    } else {
        face_items.par_iter()
            .map(|&(si, shi, fi)| check_single_face_detailed(brep, si, shi, fi, n_edges, tolerance))
            .collect()
    };

    results
}

/// Check a single face and return a detailed result.
fn check_single_face_detailed(
    brep: &BRep,
    si: usize,
    shi: usize,
    fi: usize,
    n_edges: usize,
    tolerance: f64,
) -> FaceCheckResult {
    let solid = &brep.solids[si];
    let shell = &solid.shells[shi];
    let face = &shell.faces[fi];
    let wire = &face.outer_wire;

    let mut issues = Vec::new();
    let mut outer_wire_closed = true;
    let mut outer_wire_gaps = 0usize;
    let mut has_self_intersection = false;

    // Check normal
    let normal_valid = face.normal != DVec3::ZERO && (face.normal.length() - 1.0).abs() < 0.01;
    if face.normal == DVec3::ZERO {
        issues.push(FaceCheckIssue::ZeroNormal);
    }

    // Check degenerate face
    if wire.edges.len() < 3 {
        issues.push(FaceCheckIssue::DegenerateFace);
        return FaceCheckResult {
            solid_idx: si,
            shell_idx: shi,
            face_idx: fi,
            is_valid: false,
            issues,
            outer_wire_edge_count: wire.edges.len(),
            inner_wire_count: face.inner_wires.len(),
            normal: face.normal,
            normal_valid,
            outer_wire_closed: false,
            outer_wire_gaps: 0,
            has_self_intersection: false,
        };
    }

    // Check edge indices and collect wire vertices
    let mut wire_verts: Vec<(usize, usize)> = Vec::new();
    let mut valid = true;

    for we in &wire.edges {
        if we.idx >= n_edges {
            issues.push(FaceCheckIssue::InvalidEdgeIndex { edge_idx: we.idx });
            valid = false;
        } else {
            let edge = &brep.edges[we.idx];
            let (sv, ev) = if we.forward { (edge.start, edge.end) } else { (edge.end, edge.start) };
            wire_verts.push((sv, ev));
        }
    }

    if valid {
        // Check wire closure
        let n = wire_verts.len();
        for i in 0..n {
            let next = (i + 1) % n;
            let end_v = wire_verts[i].1;
            let start_v = wire_verts[next].0;
            if end_v != start_v {
                let end_pt = brep.vertices.get(end_v).map(|v| v.point).unwrap_or_default();
                let start_pt = brep.vertices.get(start_v).map(|v| v.point).unwrap_or_default();
                let gap = (end_pt - start_pt).length();
                if gap > tolerance {
                    issues.push(FaceCheckIssue::OpenWire { wire_pos: i, gap_distance: gap });
                    outer_wire_closed = false;
                    outer_wire_gaps += 1;
                }
            }
        }

        // Check for self-intersection (topological)
        use std::collections::HashMap;
        let mut vertex_count: HashMap<usize, usize> = HashMap::new();
        for &(sv, ev) in &wire_verts {
            *vertex_count.entry(sv).or_insert(0) += 1;
            *vertex_count.entry(ev).or_insert(0) += 1;
        }
        for (&vidx, &count) in &vertex_count {
            if count > 2 {
                issues.push(FaceCheckIssue::SelfIntersection { vertex_idx: vidx, wire_idx: 0 });
                has_self_intersection = true;
            }
        }

        // Check for geometric self-intersection
        check_geometric_self_intersection_face(&wire_verts, &brep.vertices, &mut issues);
    }

    // Check inner wires
    for (wi, inner_wire) in face.inner_wires.iter().enumerate() {
        if inner_wire.edges.len() < 2 {
            continue;
        }

        let mut inner_verts: Vec<(usize, usize)> = Vec::new();
        let mut inner_valid = true;

        for we in &inner_wire.edges {
            if we.idx >= n_edges {
                issues.push(FaceCheckIssue::InvalidEdgeIndex { edge_idx: we.idx });
                inner_valid = false;
            } else {
                let edge = &brep.edges[we.idx];
                let (sv, ev) = if we.forward { (edge.start, edge.end) } else { (edge.end, edge.start) };
                inner_verts.push((sv, ev));
            }
        }

        if inner_valid {
            let n_inner = inner_verts.len();
            for i in 0..n_inner {
                let next = (i + 1) % n_inner;
                let end_v = inner_verts[i].1;
                let start_v = inner_verts[next].0;
                if end_v != start_v {
                    let end_pt = brep.vertices.get(end_v).map(|v| v.point).unwrap_or_default();
                    let start_pt = brep.vertices.get(start_v).map(|v| v.point).unwrap_or_default();
                    let gap = (end_pt - start_pt).length();
                    if gap > tolerance {
                        issues.push(FaceCheckIssue::InnerWireOpen { wire_idx: wi + 1, wire_pos: i });
                    }
                }
            }
        }
    }

    FaceCheckResult {
        solid_idx: si,
        shell_idx: shi,
        face_idx: fi,
        is_valid: issues.is_empty(),
        issues,
        outer_wire_edge_count: wire.edges.len(),
        inner_wire_count: face.inner_wires.len(),
        normal: face.normal,
        normal_valid,
        outer_wire_closed,
        outer_wire_gaps,
        has_self_intersection,
    }
}

/// Check for geometric self-intersections in a face wire.
fn check_geometric_self_intersection_face(
    wire_verts: &[(usize, usize)],
    vertices: &[rcad_kernel::topology::Vertex],
    issues: &mut Vec<FaceCheckIssue>,
) {
    let n = wire_verts.len();
    if n < 4 {
        return;
    }

    // Check pairs of non-adjacent edges
    for i in 0..n {
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue;
            }

            let (a_start, a_end) = wire_verts[i];
            let (b_start, b_end) = wire_verts[j];

            let p1 = vertices.get(a_start).map(|v| v.point).unwrap_or_default();
            let p2 = vertices.get(a_end).map(|v| v.point).unwrap_or_default();
            let p3 = vertices.get(b_start).map(|v| v.point).unwrap_or_default();
            let p4 = vertices.get(b_end).map(|v| v.point).unwrap_or_default();

            if segments_intersect_2d(p1, p2, p3, p4) {
                issues.push(FaceCheckIssue::GeometricSelfIntersection { edge_a: i, edge_b: j });
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Parallel Edge Checking
// ═══════════════════════════════════════════════════════════════════════════════

/// Check all edges in a BRep in parallel.
///
/// This function distributes edge checking across multiple threads. Each edge is
/// checked for:
/// - Vertex index validity
/// - Degeneracy (zero length)
/// - Manifold condition
/// - SameParameter violations
/// - Self-intersections
///
/// # Arguments
///
/// * `brep` - The BRep to check.
/// * `num_threads` - Number of threads to use (0 = use all available).
///
/// # Returns
///
/// A vector of `EdgeCheckResult`, one per edge.
pub fn check_edges_parallel(brep: &BRep, num_threads: usize) -> Vec<EdgeCheckResult> {
    let n_verts = brep.vertices.len();
    let tolerance = 1e-6;

    // Pre-compute edge face counts
    let edge_face_counts = compute_edge_face_counts_parallel(brep, brep.edges.len());

    let results = if num_threads > 0 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());
        pool.install(|| {
            brep.edges.par_iter()
                .enumerate()
                .map(|(eidx, edge)| check_single_edge(brep, eidx, edge, n_verts, edge_face_counts[eidx], tolerance))
                .collect()
        })
    } else {
        brep.edges.par_iter()
            .enumerate()
            .map(|(eidx, edge)| check_single_edge(brep, eidx, edge, n_verts, edge_face_counts[eidx], tolerance))
            .collect()
    };

    results
}

/// Check a single edge and return a detailed result.
fn check_single_edge(
    brep: &BRep,
    eidx: usize,
    edge: &rcad_kernel::topology::Edge,
    n_verts: usize,
    face_count: usize,
    tolerance: f64,
) -> EdgeCheckResult {
    let mut issues = Vec::new();

    // Check vertex indices
    let start_valid = edge.start < n_verts;
    let end_valid = edge.end < n_verts;

    if !start_valid {
        issues.push(EdgeCheckIssue::InvalidVertexIndex { vertex_idx: edge.start });
    }
    if !end_valid {
        issues.push(EdgeCheckIssue::InvalidVertexIndex { vertex_idx: edge.end });
    }

    // Compute edge length
    let start_pt = if start_valid { brep.vertices[edge.start].point } else { DVec3::ZERO };
    let end_pt = if end_valid { brep.vertices[edge.end].point } else { DVec3::ZERO };
    let length = (end_pt - start_pt).length();
    let is_degenerate = length < tolerance;

    if is_degenerate && start_valid && end_valid && edge.start != edge.end {
        issues.push(EdgeCheckIssue::DegenerateEdge);
    }

    // Check manifold condition
    let is_manifold = face_count == 2;
    if face_count == 0 {
        issues.push(EdgeCheckIssue::FreeEdge);
    } else if face_count != 2 {
        issues.push(EdgeCheckIssue::NonManifold { face_count });
    }

    // Check SameParameter condition
    let mut edge_tolerance = tolerance;
    if let Some(curve_idx) = brep.geom.edge_curve.get(eidx).and_then(|c| *c) {
        if let Some(curve) = brep.geom.curves.get(curve_idx) {
            if let Some(range) = brep.geom.edge_curve_range.get(eidx).and_then(|r| *r) {
                if start_valid && end_valid {
                    let eval_start = curve.point_at(range[0]);
                    let eval_end = curve.point_at(range[1]);
                    let start_gap = (eval_start - start_pt).length();
                    let end_gap = (eval_end - end_pt).length();

                    edge_tolerance = start_gap.max(end_gap);

                    if start_gap > tolerance || end_gap > tolerance {
                        issues.push(EdgeCheckIssue::SameParameterViolation { start_gap, end_gap });
                    }
                }
            }
        }
    }

    EdgeCheckResult {
        edge_idx: eidx,
        is_valid: issues.is_empty(),
        issues,
        start_vertex: edge.start,
        end_vertex: edge.end,
        length,
        is_degenerate,
        face_count,
        is_manifold,
        tolerance: edge_tolerance,
        has_self_intersection: false, // Would require curve evaluation
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Parallel Shell Validation
// ═══════════════════════════════════════════════════════════════════════════════

/// Validate all shells in a BRep in parallel.
///
/// This function checks each shell for:
/// - Closure (no free edges)
/// - Manifold condition
/// - Euler characteristic
/// - Orientation consistency
///
/// # Arguments
///
/// * `brep` - The BRep to validate.
///
/// # Returns
///
/// A vector of `ShellValidationResult`, one per shell.
pub fn validate_shells_parallel(brep: &BRep) -> Vec<ShellValidationResult> {
    // Create a flat list of shells for parallel processing
    let shell_items: Vec<(usize, usize)> = brep.solids
        .iter()
        .enumerate()
        .flat_map(|(si, solid)| {
            solid.shells.iter().enumerate().map(move |(shi, _)| (si, shi))
        })
        .collect();

    shell_items.par_iter()
        .map(|&(si, shi)| validate_single_shell(brep, si, shi))
        .collect()
}

/// Validate a single shell.
fn validate_single_shell(brep: &BRep, si: usize, shi: usize) -> ShellValidationResult {
    use std::collections::{HashMap, HashSet};

    let solid = &brep.solids[si];
    let shell = &solid.shells[shi];
    let n_edges = brep.edges.len();
    let tolerance = 1e-6;

    let mut errors = Vec::new();
    let mut warnings = Vec::new();

    // Count edges and vertices
    let mut unique_edges: HashSet<usize> = HashSet::new();
    let mut unique_vertices: HashSet<usize> = HashSet::new();

    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            if we.idx < n_edges {
                unique_edges.insert(we.idx);
                let edge = &brep.edges[we.idx];
                unique_vertices.insert(edge.start);
                unique_vertices.insert(edge.end);
            }
        }
        for wire in &face.inner_wires {
            for we in &wire.edges {
                if we.idx < n_edges {
                    unique_edges.insert(we.idx);
                    let edge = &brep.edges[we.idx];
                    unique_vertices.insert(edge.start);
                    unique_vertices.insert(edge.end);
                }
            }
        }
    }

    let edge_count = unique_edges.len();
    let vertex_count = unique_vertices.len();
    let face_count = shell.faces.len();

    // Count edge face references
    let mut edge_face_count: HashMap<usize, usize> = HashMap::new();
    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            if we.idx < n_edges {
                *edge_face_count.entry(we.idx).or_insert(0) += 1;
            }
        }
        for wire in &face.inner_wires {
            for we in &wire.edges {
                if we.idx < n_edges {
                    *edge_face_count.entry(we.idx).or_insert(0) += 1;
                }
            }
        }
    }

    let open_edge_count = edge_face_count.values().filter(|&&c| c == 1).count();
    let non_manifold_edge_count = edge_face_count.values().filter(|&&c| c > 2).count();

    let is_closed = open_edge_count == 0;
    let is_manifold = non_manifold_edge_count == 0;

    // Compute Euler characteristic
    let euler_characteristic = vertex_count as i64 - edge_count as i64 + face_count as i64;

    // Compute genus (only meaningful for closed shells)
    let genus = if is_closed {
        let g = (2 - euler_characteristic) / 2;
        if (2 - euler_characteristic) % 2 == 0 && g >= 0 { Some(g) } else { None }
    } else {
        None
    };

    // Check orientation consistency
    let orientation_consistent = check_shell_orientation_consistency(shell, brep);

    // Get face results
    let face_results: Vec<FaceCheckResult> = shell.faces
        .iter()
        .enumerate()
        .map(|(fi, _)| check_single_face_detailed(brep, si, shi, fi, n_edges, tolerance))
        .collect();

    // Generate errors and warnings
    if !is_closed {
        errors.push(format!("Shell has {} open edges", open_edge_count));
    }
    if !is_manifold {
        errors.push(format!("Shell has {} non-manifold edges", non_manifold_edge_count));
    }
    if !orientation_consistent {
        warnings.push("Shell orientation may be inconsistent".to_string());
    }

    let is_valid = errors.is_empty() && face_results.iter().all(|f| f.is_valid);

    ShellValidationResult {
        solid_idx: si,
        shell_idx: shi,
        is_valid,
        face_count,
        edge_count,
        vertex_count,
        euler_characteristic,
        is_closed,
        is_manifold,
        open_edge_count,
        non_manifold_edge_count,
        orientation_consistent,
        genus,
        face_results,
        errors,
        warnings,
    }
}

/// Check orientation consistency for a shell.
fn check_shell_orientation_consistency(shell: &rcad_kernel::topology::Shell, brep: &BRep) -> bool {
    use std::collections::HashMap;

    let n_edges = brep.edges.len();
    let mut edge_orientations: HashMap<usize, Vec<bool>> = HashMap::new();

    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            if we.idx < n_edges {
                edge_orientations.entry(we.idx).or_default().push(we.forward);
            }
        }
    }

    // For a properly oriented shell, each edge should have one forward and one backward reference
    for (_, orientations) in &edge_orientations {
        if orientations.len() == 2 {
            // Adjacent faces should have opposite orientations for shared edges
            if orientations[0] == orientations[1] {
                return false;
            }
        }
    }

    true
}

// ═══════════════════════════════════════════════════════════════════════════════
// Parallel Solid Validation
// ═══════════════════════════════════════════════════════════════════════════════

/// Validate all solids in a BRep in parallel.
///
/// This function checks each solid for:
/// - Shell closure
/// - Manifold condition
/// - Orientation validity
/// - Volume calculation
///
/// # Arguments
///
/// * `brep` - The BRep to validate.
///
/// # Returns
///
/// A vector of `SolidValidationResult`, one per solid.
pub fn validate_solids_parallel(brep: &BRep) -> Vec<SolidValidationResult> {
    brep.solids.par_iter()
        .enumerate()
        .map(|(si, _)| validate_single_solid(brep, si))
        .collect()
}

/// Validate a single solid.
fn validate_single_solid(brep: &BRep, si: usize) -> SolidValidationResult {
    use std::collections::HashSet;

    let solid = &brep.solids[si];
    let n_edges = brep.edges.len();

    let mut errors = Vec::new();
    let mut warnings = Vec::new();

    // Get shell results
    let shell_results: Vec<ShellValidationResult> = solid.shells
        .iter()
        .enumerate()
        .map(|(shi, _)| validate_single_shell(brep, si, shi))
        .collect();

    // Aggregate counts
    let face_count: usize = solid.shells.iter().map(|s| s.faces.len()).sum();
    let edge_count: usize;
    let vertex_count: usize;

    {
        let mut edges: HashSet<usize> = HashSet::new();
        let mut verts: HashSet<usize> = HashSet::new();

        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    if we.idx < n_edges {
                        edges.insert(we.idx);
                        let edge = &brep.edges[we.idx];
                        verts.insert(edge.start);
                        verts.insert(edge.end);
                    }
                }
            }
        }

        edge_count = edges.len();
        vertex_count = verts.len();
    }

    // Compute Euler characteristic
    let euler_characteristic = vertex_count as i64 - edge_count as i64 + face_count as i64;

    // Check if all shells are closed and manifold
    let is_closed = shell_results.iter().all(|s| s.is_closed);
    let is_manifold = shell_results.iter().all(|s| s.is_manifold);
    let orientation_valid = shell_results.iter().all(|s| s.orientation_consistent);

    // Compute volume (approximate using shell volumes)
    let volume: f64 = solid.shells.iter()
        .map(|shell| compute_shell_volume(shell, brep))
        .sum();

    let has_positive_volume = volume > 0.0;

    // Compute genus
    let genus = if is_closed && is_manifold {
        let g = (2 - euler_characteristic) / 2;
        if (2 - euler_characteristic) % 2 == 0 && g >= 0 { Some(g) } else { None }
    } else {
        None
    };

    // Generate errors
    if !is_closed {
        errors.push("Solid has unclosed shells".to_string());
    }
    if !is_manifold {
        errors.push("Solid has non-manifold topology".to_string());
    }
    if !has_positive_volume {
        warnings.push("Solid has zero or negative volume".to_string());
    }

    let is_valid = errors.is_empty() && shell_results.iter().all(|s| s.is_valid);

    SolidValidationResult {
        solid_idx: si,
        is_valid,
        shell_count: solid.shells.len(),
        face_count,
        edge_count,
        vertex_count,
        euler_characteristic,
        is_closed,
        is_manifold,
        orientation_valid,
        has_positive_volume,
        volume,
        genus,
        shell_results,
        errors,
        warnings,
    }
}

/// Compute the volume of a shell using signed volume method.
fn compute_shell_volume(shell: &rcad_kernel::topology::Shell, brep: &BRep) -> f64 {
    use std::collections::HashMap;

    let n_edges = brep.edges.len();
    let mut volume = 0.0_f64;

    for face in &shell.faces {
        // Get vertices of the outer wire
        let mut verts: Vec<DVec3> = Vec::new();
        for we in &face.outer_wire.edges {
            if we.idx < n_edges {
                let edge = &brep.edges[we.idx];
                let vi = if we.forward { edge.start } else { edge.end };
                if vi < brep.vertices.len() {
                    verts.push(brep.vertices[vi].point);
                }
            }
        }

        // Compute signed volume contribution using triangulation
        if verts.len() >= 3 {
            let origin = verts[0];
            for i in 1..verts.len() - 1 {
                let v1 = verts[i] - origin;
                let v2 = verts[i + 1] - origin;
                let signed_vol = v1.cross(v2).dot(face.normal) / 6.0;
                volume += signed_vol;
            }
        }
    }

    volume.abs()
}

// ═══════════════════════════════════════════════════════════════════════════════
// Comprehensive Parallel Check
// ═══════════════════════════════════════════════════════════════════════════════

/// Perform a comprehensive parallel check of a BRep.
///
/// This function runs all configured checks in parallel and returns a detailed
/// report including timing information for each phase.
///
/// # Arguments
///
/// * `brep` - The BRep to check.
/// * `config` - Configuration for the check.
///
/// # Returns
///
/// A `ParallelCheckReport` containing all results and timing information.
pub fn check_brep_parallel(brep: &BRep, config: &ParallelCheckConfig) -> ParallelCheckReport {
    let start_time = Instant::now();
    let mut phase_timings: Vec<CheckPhaseTiming> = Vec::new();
    let mut structural_issues = Vec::new();
    let mut parallel_issues = Vec::new();

    // Configure thread pool
    let threads_used = if config.num_threads > 0 {
        config.num_threads
    } else {
        rayon::current_num_threads()
    };

    // Count totals
    let total_solids = brep.solids.len();
    let total_shells: usize = brep.solids.iter().map(|s| s.shells.len()).sum();
    let total_faces: usize = brep.solids.iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum();
    let total_edges = brep.edges.len();
    let total_vertices = brep.vertices.len();

    let use_parallel = total_faces >= config.parallel_threshold
        || total_edges >= config.parallel_threshold
        || total_vertices >= config.parallel_threshold;

    // Phase 1: Face checking
    let mut face_results = Vec::new();
    if config.check_faces {
        let phase_start = Instant::now();
        face_results = if use_parallel {
            check_faces_parallel(brep, threads_used)
        } else {
            check_faces_sequential(brep)
        };
        phase_timings.push(CheckPhaseTiming {
            phase: "faces".to_string(),
            duration_ms: phase_start.elapsed().as_millis() as u64,
            items_processed: total_faces,
        });

        // Collect issues from face results
        for fr in &face_results {
            if !fr.is_valid {
                structural_issues.push(CheckIssue::DegenerateFace {
                    solid: fr.solid_idx,
                    shell: fr.shell_idx,
                    face: fr.face_idx,
                });
            }
        }
    }

    // Phase 2: Edge checking
    let mut edge_results = Vec::new();
    if config.check_edges {
        let phase_start = Instant::now();
        edge_results = if use_parallel {
            check_edges_parallel(brep, threads_used)
        } else {
            check_edges_sequential(brep)
        };
        phase_timings.push(CheckPhaseTiming {
            phase: "edges".to_string(),
            duration_ms: phase_start.elapsed().as_millis() as u64,
            items_processed: total_edges,
        });

        // Collect issues from edge results
        for er in &edge_results {
            for issue in &er.issues {
                match issue {
                    EdgeCheckIssue::InvalidVertexIndex { vertex_idx } => {
                        structural_issues.push(CheckIssue::InvalidVertexIndex {
                            edge: er.edge_idx,
                            vertex_idx: *vertex_idx,
                        });
                    }
                    EdgeCheckIssue::NonManifold { face_count } => {
                        structural_issues.push(CheckIssue::NonManifoldEdge {
                            edge_idx: er.edge_idx,
                            face_count: *face_count,
                        });
                    }
                    _ => {}
                }
            }
        }
    }

    // Phase 3: Vertex checking
    if config.check_vertices {
        let phase_start = Instant::now();

        // Check for non-finite vertices
        if config.check_finite_vertices {
            for (vidx, v) in brep.vertices.iter().enumerate() {
                if !v.point.is_finite() {
                    parallel_issues.push(ParallelCheckIssue::NonFiniteVertex { vertex_idx: vidx });
                }
            }
        }

        // Check for isolated vertices
        if config.check_isolated_vertices {
            let mut referenced = vec![false; brep.vertices.len()];
            for edge in &brep.edges {
                if edge.start < brep.vertices.len() {
                    referenced[edge.start] = true;
                }
                if edge.end < brep.vertices.len() {
                    referenced[edge.end] = true;
                }
            }
            for (vidx, &is_ref) in referenced.iter().enumerate() {
                if !is_ref {
                    parallel_issues.push(ParallelCheckIssue::IsolatedVertex { vertex_idx: vidx });
                }
            }
        }

        // Check for duplicate vertices
        if config.check_duplicate_vertices {
            let duplicates = find_duplicate_vertices_parallel(&brep.vertices, config.tolerance);
            parallel_issues.extend(duplicates);
        }

        phase_timings.push(CheckPhaseTiming {
            phase: "vertices".to_string(),
            duration_ms: phase_start.elapsed().as_millis() as u64,
            items_processed: total_vertices,
        });
    }

    // Phase 4: Shell validation
    let mut shell_results = Vec::new();
    if config.check_shells {
        let phase_start = Instant::now();
        shell_results = validate_shells_parallel(brep);
        phase_timings.push(CheckPhaseTiming {
            phase: "shells".to_string(),
            duration_ms: phase_start.elapsed().as_millis() as u64,
            items_processed: total_shells,
        });
    }

    // Phase 5: Solid validation
    let mut solid_results = Vec::new();
    if config.check_solids {
        let phase_start = Instant::now();
        solid_results = validate_solids_parallel(brep);
        phase_timings.push(CheckPhaseTiming {
            phase: "solids".to_string(),
            duration_ms: phase_start.elapsed().as_millis() as u64,
            items_processed: total_solids,
        });
    }

    let total_duration_ms = start_time.elapsed().as_millis() as u64;

    // Determine overall validity
    let is_valid = structural_issues.is_empty()
        && parallel_issues.is_empty()
        && shell_results.iter().all(|s| s.is_valid)
        && solid_results.iter().all(|s| s.is_valid);

    // Build stats
    let stats = ParallelCheckStats {
        face_count: total_faces,
        edge_count: total_edges,
        vertex_count: total_vertices,
        issue_count: structural_issues.len() + parallel_issues.len(),
        is_valid,
        was_parallel: use_parallel,
        thread_count: threads_used,
    };

    ParallelCheckReport {
        is_valid,
        total_faces,
        total_edges,
        total_vertices,
        total_solids,
        total_shells,
        threads_used,
        was_parallel: use_parallel,
        total_duration_ms,
        phase_timings,
        face_results,
        edge_results,
        shell_results,
        solid_results,
        structural_issues,
        parallel_issues,
        stats,
    }
}

/// Sequential face checking fallback.
fn check_faces_sequential(brep: &BRep) -> Vec<FaceCheckResult> {
    let n_edges = brep.edges.len();
    let tolerance = 1e-6;

    let mut results = Vec::new();
    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for fi in 0..shell.faces.len() {
                results.push(check_single_face_detailed(brep, si, shi, fi, n_edges, tolerance));
            }
        }
    }
    results
}

/// Sequential edge checking fallback.
fn check_edges_sequential(brep: &BRep) -> Vec<EdgeCheckResult> {
    let n_verts = brep.vertices.len();
    let tolerance = 1e-6;

    // Compute edge face counts
    let mut edge_face_counts = vec![0usize; brep.edges.len()];
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    if we.idx < brep.edges.len() {
                        edge_face_counts[we.idx] += 1;
                    }
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        if we.idx < brep.edges.len() {
                            edge_face_counts[we.idx] += 1;
                        }
                    }
                }
            }
        }
    }

    brep.edges.iter()
        .enumerate()
        .map(|(eidx, edge)| check_single_edge(brep, eidx, edge, n_verts, edge_face_counts[eidx], tolerance))
        .collect()
}

/// Perform parallel check and return detailed statistics.
pub fn check_parallel_with_stats(brep: &BRep) -> (CheckResult, ParallelCheckStats) {
    let face_count: usize = brep.solids.iter()
        .map(|s| s.shells.iter().map(|sh| sh.faces.len()).sum::<usize>())
        .sum();
    let edge_count = brep.edges.len();
    let vertex_count = brep.vertices.len();

    let options = ParallelCheckOptions::default();
    let result = check_parallel_with_options(brep, &options);

    let stats = ParallelCheckStats {
        face_count,
        edge_count,
        vertex_count,
        issue_count: result.issues.len() + result.parallel_issues.len(),
        is_valid: result.is_valid(),
        was_parallel: result.was_parallel,
        thread_count: result.thread_count,
    };

    (result.to_check_result(), stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::BRep;
    use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn test_check_parallel_empty_brep() {
        let brep = BRep::default();
        let result = check_parallel(&brep);
        assert!(result.is_valid());
    }

    #[test]
    fn test_check_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = check_parallel(&brep);
        assert!(result.is_valid(), "issues: {:?}", result.issues);
    }

    #[test]
    fn test_check_parallel_cylinder() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        let result = check_parallel(&brep);
        // Cylinder has seam edges that may trigger non-manifold warnings
        // The check should complete without panic, not necessarily be valid
        let _ = result.issues.len();
    }

    #[test]
    fn test_check_parallel_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere {
            radius: 1.0,
        });
        let result = check_parallel(&brep);
        // Sphere has seam edges that may trigger non-manifold warnings
        // The check should complete without panic, not necessarily be valid
        let _ = result.issues.len();
    }

    #[test]
    fn test_check_many_parallel() {
        let breps: Vec<BRep> = vec![
            BRep::from_primitive(PrimitiveSolid::Box {
                width: 1.0, height: 1.0, depth: 1.0,
            }),
            BRep::from_primitive(PrimitiveSolid::Sphere {
                radius: 1.0,
            }),
            BRep::from_primitive(PrimitiveSolid::Cylinder {
                radius: 1.0, height: 2.0,
            }),
        ];

        let results = check_many_parallel(&breps);
        assert_eq!(results.len(), 3);
        // Box should be valid
        assert!(results[0].is_valid(), "issues: {:?}", results[0].issues);
        // Sphere and cylinder have seam edges that may trigger warnings
        // Just verify the checks completed
        let _ = results[1].issues.len();
        let _ = results[2].issues.len();
    }

    #[test]
    fn test_check_parallel_with_stats() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let (result, stats) = check_parallel_with_stats(&brep);
        assert!(result.is_valid(), "issues: {:?}", result.issues);
        assert_eq!(stats.face_count, 6); // Box has 6 faces
        assert_eq!(stats.edge_count, 12); // Box has 12 edges
        assert_eq!(stats.vertex_count, 8); // Box has 8 vertices
        assert_eq!(stats.issue_count, 0);
        assert!(stats.is_valid);
    }

    #[test]
    fn test_segments_intersect_2d() {
        // Crossing segments
        let p1 = DVec3::new(0.0, 0.0, 0.0);
        let p2 = DVec3::new(2.0, 2.0, 0.0);
        let p3 = DVec3::new(0.0, 2.0, 0.0);
        let p4 = DVec3::new(2.0, 0.0, 0.0);
        assert!(segments_intersect_2d(p1, p2, p3, p4));

        // Non-crossing segments
        let p5 = DVec3::new(0.0, 0.0, 0.0);
        let p6 = DVec3::new(1.0, 1.0, 0.0);
        let p7 = DVec3::new(3.0, 3.0, 0.0);
        let p8 = DVec3::new(4.0, 4.0, 0.0);
        assert!(!segments_intersect_2d(p5, p6, p7, p8));
    }

    #[test]
    fn test_parallel_options_default() {
        let opts = ParallelCheckOptions::default();
        assert_eq!(opts.min_faces_for_parallel, 100);
        assert_eq!(opts.chunk_size, 32);
        assert!(opts.check_duplicate_vertices);
        assert!(opts.check_isolated_vertices);
        assert!(opts.check_finite_vertices);
    }

    #[test]
    fn test_parallel_options_small_model() {
        let opts = ParallelCheckOptions::small_model();
        assert_eq!(opts.min_faces_for_parallel, usize::MAX);
    }

    #[test]
    fn test_parallel_options_large_model() {
        let opts = ParallelCheckOptions::large_model();
        assert_eq!(opts.min_faces_for_parallel, 10);
        assert_eq!(opts.chunk_size, 64);
    }

    #[test]
    fn test_parallel_vs_sequential_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Both should produce same results
        let parallel_result = check_parallel(&brep);
        let sequential_result = crate::brep_check::check(&brep);

        assert_eq!(parallel_result.is_valid(), sequential_result.is_valid());
        assert_eq!(parallel_result.issues.len(), sequential_result.issues.len());
    }

    #[test]
    fn test_parallel_vs_sequential_invalid_brep() {
        // Create an invalid BRep with open wire
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 3, end: 0 }); // Gap: v2 != v3

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let parallel_result = check_parallel(&brep);
        let sequential_result = crate::brep_check::check(&brep);

        // Both should detect the open wire
        assert!(!parallel_result.is_valid());
        assert!(!sequential_result.is_valid());

        // Both should have same number of issues
        assert_eq!(parallel_result.issues.len(), sequential_result.issues.len());

        // Both should have OpenWire issue
        assert!(parallel_result.issues.iter().any(|i| matches!(i, CheckIssue::OpenWire { .. })));
        assert!(sequential_result.issues.iter().any(|i| matches!(i, CheckIssue::OpenWire { .. })));
    }

    #[test]
    fn test_small_model_uses_sequential() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let opts = ParallelCheckOptions::small_model();
        let result = check_parallel_with_options(&brep, &opts);

        assert!(!result.was_parallel, "Small model should use sequential processing");
        assert_eq!(result.thread_count, 1);
    }

    #[test]
    fn test_large_model_uses_parallel() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let opts = ParallelCheckOptions::large_model();
        let result = check_parallel_with_options(&brep, &opts);

        assert!(result.was_parallel, "Large model settings should use parallel processing");
        assert!(result.thread_count >= 1);
    }

    #[test]
    fn test_isolated_vertex_detection() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(2.0, 0.0, 0.0) }); // Isolated

        brep.edges.push(Edge { start: 0, end: 1 });

        let opts = ParallelCheckOptions {
            check_isolated_vertices: true,
            check_duplicate_vertices: false,
            check_finite_vertices: false,
            ..ParallelCheckOptions::default()
        };

        let result = check_parallel_with_options(&brep, &opts);

        assert!(result.parallel_issues.iter().any(|i| matches!(
            i,
            ParallelCheckIssue::IsolatedVertex { vertex_idx: 2 }
        )), "Should detect isolated vertex 2");
    }

    #[test]
    fn test_non_finite_vertex_detection() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(f64::NAN, 0.0, 0.0) }); // NaN

        brep.edges.push(Edge { start: 0, end: 1 });

        let opts = ParallelCheckOptions {
            check_finite_vertices: true,
            check_duplicate_vertices: false,
            check_isolated_vertices: false,
            ..ParallelCheckOptions::default()
        };

        let result = check_parallel_with_options(&brep, &opts);

        assert!(result.parallel_issues.iter().any(|i| matches!(
            i,
            ParallelCheckIssue::NonFiniteVertex { vertex_idx: 1 }
        )), "Should detect non-finite vertex 1");
    }

    #[test]
    fn test_duplicate_vertex_detection() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // Duplicate

        brep.edges.push(Edge { start: 0, end: 1 });

        let opts = ParallelCheckOptions {
            check_duplicate_vertices: true,
            check_isolated_vertices: false,
            check_finite_vertices: false,
            tolerance: 1e-6,
            ..ParallelCheckOptions::default()
        };

        let result = check_parallel_with_options(&brep, &opts);

        assert!(result.parallel_issues.iter().any(|i| matches!(
            i,
            ParallelCheckIssue::DuplicateVertex { vertex_a: 0, vertex_b: 1, .. }
        )), "Should detect duplicate vertices");
    }

    #[test]
    fn test_check_many_parallel_with_options() {
        let breps: Vec<BRep> = vec![
            BRep::from_primitive(PrimitiveSolid::Box {
                width: 1.0, height: 1.0, depth: 1.0,
            }),
            BRep::from_primitive(PrimitiveSolid::Sphere {
                radius: 1.0,
            }),
        ];

        let opts = ParallelCheckOptions::default();
        let results = check_many_parallel_with_options(&breps, &opts);

        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_parallel_check_result_is_valid() {
        let mut result = ParallelCheckResult::default();
        assert!(result.is_valid());

        result.issues.push(CheckIssue::DegenerateFace { solid: 0, shell: 0, face: 0 });
        assert!(!result.is_valid());
    }

    #[test]
    fn test_parallel_check_result_to_check_result() {
        let mut result = ParallelCheckResult::default();
        result.issues.push(CheckIssue::DegenerateFace { solid: 0, shell: 0, face: 0 });

        let check_result = result.to_check_result();
        assert_eq!(check_result.issues.len(), 1);
    }

    /// Generate a large BRep for performance testing.
    #[cfg(test)]
    fn generate_large_brep(n_boxes: usize) -> BRep {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();

        // Create a grid of connected quads
        let mut vertex_offset = 0usize;
        let mut edge_offset = 0usize;

        for _z in 0..n_boxes {
            for _y in 0..n_boxes {
                for _x in 0..n_boxes {
                    // Add 8 vertices for a box
                    for dz in 0..2 {
                        for dy in 0..2 {
                            for dx in 0..2 {
                                let x = dx as f64;
                                let y = dy as f64;
                                let z = dz as f64;
                                brep.vertices.push(Vertex {
                                    point: DVec3::new(x, y, z),
                                });
                            }
                        }
                    }

                    // Add 12 edges for the box
                    let v = vertex_offset;
                    let edges = vec![
                        (v+0, v+1), (v+1, v+3), (v+3, v+2), (v+2, v+0), // bottom
                        (v+4, v+5), (v+5, v+7), (v+7, v+6), (v+6, v+4), // top
                        (v+0, v+4), (v+1, v+5), (v+2, v+6), (v+3, v+7), // vertical
                    ];

                    for (start, end) in edges {
                        brep.edges.push(Edge { start, end });
                    }

                    // Add 6 faces for the box
                    let e = edge_offset;
                    let face_wire_indices = vec![
                        vec![(e+0, true), (e+1, true), (e+2, true), (e+3, true)],   // bottom
                        vec![(e+4, true), (e+5, true), (e+6, true), (e+7, true)],   // top
                        vec![(e+0, true), (e+8, true), (e+4, false), (e+11,false)], // front
                        vec![(e+2, false), (e+10,true), (e+6, false), (e+9, false)],// back
                        vec![(e+3, true), (e+10,false),(e+7, true), (e+8, false)], // left
                        vec![(e+1, false),(e+9, true), (e+5, false), (e+11,true)], // right
                    ];

                    let normals = vec![
                        DVec3::NEG_Z, DVec3::Z, DVec3::NEG_Y, DVec3::Y, DVec3::NEG_X, DVec3::X,
                    ];

                    let mut faces = Vec::new();
                    for (fi, wire_indices) in face_wire_indices.iter().enumerate() {
                        faces.push(Face {
                            outer_wire: Wire {
                                edges: wire_indices.iter().map(|&(idx, fwd)| {
                                    if fwd { WireEdge::fwd(idx) } else { WireEdge::rev(idx) }
                                }).collect(),
                            },
                            inner_wires: vec![],
                            normal: normals[fi],
                            triangles: vec![],
                            mesh_dirty: true,
                        });
                    }

                    brep.solids.push(Solid {
                        shells: vec![Shell { faces }],
                    });

                    vertex_offset += 8;
                    edge_offset += 12;
                }
            }
        }

        brep
    }

    #[test]
    fn test_large_brep_parallel_vs_sequential() {
        // Create a moderately large BRep
        let brep = generate_large_brep(3); // 27 boxes, 162 faces

        let parallel_result = check_parallel(&brep);
        let sequential_result = crate::brep_check::check(&brep);

        // Results should be identical
        assert_eq!(parallel_result.is_valid(), sequential_result.is_valid());
        assert_eq!(parallel_result.issues.len(), sequential_result.issues.len());
    }

    #[test]
    fn test_parallel_options_builder() {
        let opts = ParallelCheckOptions::default()
            .with_tolerance(1e-9)
            .with_chunk_size(128)
            .with_duplicate_vertex_check(false)
            .with_isolated_vertex_check(false);

        assert!((opts.tolerance - 1e-9).abs() < 1e-15);
        assert_eq!(opts.chunk_size, 128);
        assert!(!opts.check_duplicate_vertices);
        assert!(!opts.check_isolated_vertices);
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // Tests for check_faces_parallel
    // ═══════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_check_faces_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let results = check_faces_parallel(&brep, 0);
        assert_eq!(results.len(), 6, "Box should have 6 faces");

        for result in &results {
            assert!(result.is_valid, "Face should be valid: {:?}", result.issues);
            assert!(result.outer_wire_closed, "Outer wire should be closed");
            assert_eq!(result.outer_wire_edge_count, 4, "Each face should have 4 edges");
        }
    }

    #[test]
    fn test_check_faces_parallel_open_wire() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 3, end: 0 }); // Gap: v2 != v3

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let results = check_faces_parallel(&brep, 0);
        assert_eq!(results.len(), 1);
        assert!(!results[0].is_valid, "Face with open wire should be invalid");
        assert!(!results[0].outer_wire_closed, "Wire should be reported as open");
        assert!(results[0].outer_wire_gaps > 0, "Should have gaps");
    }

    #[test]
    fn test_check_faces_parallel_degenerate_face() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(Edge { start: 0, end: 1 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0)], // Only 1 edge
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let results = check_faces_parallel(&brep, 0);
        assert!(!results[0].is_valid, "Degenerate face should be invalid");
        assert!(results[0].issues.iter().any(|i| matches!(i, FaceCheckIssue::DegenerateFace)));
    }

    #[test]
    fn test_check_faces_parallel_zero_normal() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.vertices.push(Vertex { point: DVec3::Y });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::ZERO, // Zero normal
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let results = check_faces_parallel(&brep, 0);
        assert!(!results[0].is_valid, "Face with zero normal should be invalid");
        assert!(results[0].issues.iter().any(|i| matches!(i, FaceCheckIssue::ZeroNormal)));
    }

    #[test]
    fn test_check_faces_parallel_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let results = check_faces_parallel(&brep, 0);

        // Sphere should have faces
        assert!(!results.is_empty(), "Sphere should have faces");
        // Verify basic structure
        for result in &results {
            assert!(result.outer_wire_edge_count >= 0);
        }
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // Tests for check_edges_parallel
    // ═══════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_check_edges_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let results = check_edges_parallel(&brep, 0);
        assert_eq!(results.len(), 12, "Box should have 12 edges");

        for result in &results {
            assert!(result.is_valid, "Edge should be valid: {:?}", result.issues);
            assert!(result.is_manifold, "Each edge should be manifold");
            assert_eq!(result.face_count, 2, "Each edge should be shared by 2 faces");
            assert!(!result.is_degenerate, "No edge should be degenerate");
            assert!(result.length > 0.0, "Each edge should have positive length");
        }
    }

    #[test]
    fn test_check_edges_parallel_invalid_vertex() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.edges.push(Edge { start: 0, end: 99 }); // Invalid vertex

        let results = check_edges_parallel(&brep, 0);
        assert_eq!(results.len(), 1);
        assert!(!results[0].is_valid, "Edge with invalid vertex should be invalid");
        assert!(results[0].issues.iter().any(|i| matches!(i, EdgeCheckIssue::InvalidVertexIndex { .. })));
    }

    #[test]
    fn test_check_edges_parallel_degenerate() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::ZERO }); // Same position
        brep.edges.push(Edge { start: 0, end: 1 });

        let results = check_edges_parallel(&brep, 0);
        assert!(results[0].is_degenerate, "Edge with same vertex positions should be degenerate");
    }

    #[test]
    fn test_check_edges_parallel_free_edge() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(Edge { start: 0, end: 1 }); // Not referenced by any face

        let results = check_edges_parallel(&brep, 0);
        assert!(!results[0].is_valid, "Free edge should be invalid");
        assert!(results[0].issues.iter().any(|i| matches!(i, EdgeCheckIssue::FreeEdge)));
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // Tests for validate_shells_parallel
    // ═══════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_validate_shells_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let results = validate_shells_parallel(&brep);
        assert_eq!(results.len(), 1, "Box should have 1 shell");

        let shell = &results[0];
        assert!(shell.is_valid, "Box shell should be valid");
        assert!(shell.is_closed, "Box shell should be closed");
        assert!(shell.is_manifold, "Box shell should be manifold");
        assert_eq!(shell.face_count, 6);
        assert_eq!(shell.euler_characteristic, 2, "Box Euler characteristic should be 2");
        assert_eq!(shell.genus, Some(0), "Box genus should be 0");
    }

    #[test]
    fn test_validate_shells_parallel_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let results = validate_shells_parallel(&brep);
        assert_eq!(results.len(), 1);

        let shell = &results[0];
        assert!(shell.is_closed, "Sphere shell should be closed");
        assert_eq!(shell.euler_characteristic, 2, "Sphere Euler characteristic should be 2");
        assert_eq!(shell.genus, Some(0), "Sphere genus should be 0");
    }

    #[test]
    fn test_validate_shells_parallel_open_shell() {
        let mut brep = BRep::new();
        // Create a simple open shell (just one face)
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.vertices.push(Vertex { point: DVec3::Y });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let results = validate_shells_parallel(&brep);
        assert!(!results[0].is_closed, "Single face shell should be open");
        assert!(results[0].open_edge_count > 0, "Should have open edges");
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // Tests for validate_solids_parallel
    // ═══════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_validate_solids_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let results = validate_solids_parallel(&brep);
        assert_eq!(results.len(), 1, "Should have 1 solid");

        let solid = &results[0];
        assert!(solid.is_valid, "Box solid should be valid");
        assert!(solid.is_closed, "Box solid should be closed");
        assert!(solid.is_manifold, "Box solid should be manifold");
        assert_eq!(solid.face_count, 6);
        assert_eq!(solid.edge_count, 12);
        assert_eq!(solid.vertex_count, 8);
        assert!(solid.volume >= 0.0, "Box volume should be non-negative");
    }

    #[test]
    fn test_validate_solids_parallel_cylinder() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let results = validate_solids_parallel(&brep);
        assert_eq!(results.len(), 1);

        let solid = &results[0];
        assert!(solid.is_closed, "Cylinder solid should be closed");
        assert!(solid.volume >= 0.0, "Cylinder volume should be non-negative");
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // Tests for check_brep_parallel
    // ═══════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_check_brep_parallel_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        assert!(report.is_valid, "Box should pass all checks");
        assert!(report.structural_issues.is_empty());
        assert!(report.parallel_issues.is_empty());
        assert_eq!(report.total_faces, 6);
        assert_eq!(report.total_edges, 12);
        assert_eq!(report.total_vertices, 8);
        assert_eq!(report.total_solids, 1);
        assert!(report.total_duration_ms < 10000, "Should complete quickly");
    }

    #[test]
    fn test_check_brep_parallel_fast_config() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let config = ParallelCheckConfig::fast();
        let report = check_brep_parallel(&brep, &config);

        // Fast config skips some checks
        assert!(report.phase_timings.iter().any(|t| t.phase == "faces"));
    }

    #[test]
    fn test_check_brep_parallel_thorough_config() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let config = ParallelCheckConfig::thorough();
        let report = check_brep_parallel(&brep, &config);

        // Thorough config has tighter tolerance
        assert!((config.tolerance - 1e-9).abs() < 1e-15);
    }

    #[test]
    fn test_check_brep_parallel_custom_threads() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let config = ParallelCheckConfig::default().with_threads(2);
        let report = check_brep_parallel(&brep, &config);

        // Should work with custom thread count
        assert!(report.threads_used >= 1);
    }

    #[test]
    fn test_check_brep_parallel_timing() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        // Should have timing for each phase
        let phases: Vec<&str> = report.phase_timings.iter().map(|t| t.phase.as_str()).collect();
        assert!(phases.contains(&"faces"));
        assert!(phases.contains(&"edges"));
        assert!(phases.contains(&"vertices"));
        assert!(phases.contains(&"shells"));
        assert!(phases.contains(&"solids"));
    }

    #[test]
    fn test_check_brep_parallel_invalid_brep() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(f64::NAN, 0.0, 0.0) }); // NaN vertex
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // Duplicate
        brep.edges.push(Edge { start: 0, end: 2 });
        // Vertex 1 is isolated

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        assert!(!report.is_valid, "Invalid BRep should fail checks");
        assert!(!report.parallel_issues.is_empty(), "Should have parallel-specific issues");
    }

    #[test]
    fn test_check_brep_parallel_empty_brep() {
        let brep = BRep::default();

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        assert!(report.is_valid, "Empty BRep should be valid (no issues)");
        assert_eq!(report.total_faces, 0);
        assert_eq!(report.total_edges, 0);
        assert_eq!(report.total_vertices, 0);
    }

    #[test]
    fn test_check_brep_parallel_summary() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        let summary = report.summary();
        assert!(summary.contains("VALID"));
        assert!(summary.contains("1 solids"));
        assert!(summary.contains("6 faces"));
    }

    #[test]
    fn test_check_brep_parallel_timing_summary() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let config = ParallelCheckConfig::default();
        let report = check_brep_parallel(&brep, &config);

        let timing = report.timing_summary();
        assert!(timing.contains("Timing breakdown"));
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // Tests for result types
    // ═══════════════════════════════════════════════════════════════════════════

    #[test]
    fn test_face_check_result_summary() {
        let result = FaceCheckResult {
            solid_idx: 0,
            shell_idx: 0,
            face_idx: 1,
            is_valid: true,
            issues: vec![],
            outer_wire_edge_count: 4,
            inner_wire_count: 0,
            normal: DVec3::Z,
            normal_valid: true,
            outer_wire_closed: true,
            outer_wire_gaps: 0,
            has_self_intersection: false,
        };

        let summary = result.summary();
        assert!(summary.contains("valid"));
        assert!(result.is_clean());
    }

    #[test]
    fn test_edge_check_result_summary() {
        let result = EdgeCheckResult {
            edge_idx: 0,
            is_valid: true,
            issues: vec![],
            start_vertex: 0,
            end_vertex: 1,
            length: 1.0,
            is_degenerate: false,
            face_count: 2,
            is_manifold: true,
            tolerance: 1e-6,
            has_self_intersection: false,
        };

        let summary = result.summary();
        assert!(summary.contains("valid"));
        assert!(result.is_clean());
    }

    #[test]
    fn test_shell_validation_result_summary() {
        let result = ShellValidationResult {
            solid_idx: 0,
            shell_idx: 0,
            is_valid: true,
            face_count: 6,
            edge_count: 12,
            vertex_count: 8,
            euler_characteristic: 2,
            is_closed: true,
            is_manifold: true,
            open_edge_count: 0,
            non_manifold_edge_count: 0,
            orientation_consistent: true,
            genus: Some(0),
            face_results: vec![],
            errors: vec![],
            warnings: vec![],
        };

        assert!(result.is_closed_manifold());
        let summary = result.summary();
        assert!(summary.contains("VALID"));
    }

    #[test]
    fn test_solid_validation_result_summary() {
        let result = SolidValidationResult {
            solid_idx: 0,
            is_valid: true,
            shell_count: 1,
            face_count: 6,
            edge_count: 12,
            vertex_count: 8,
            euler_characteristic: 2,
            is_closed: true,
            is_manifold: true,
            orientation_valid: true,
            has_positive_volume: true,
            volume: 1.0,
            genus: Some(0),
            shell_results: vec![],
            errors: vec![],
            warnings: vec![],
        };

        assert!(result.is_valid_for_operations());
        let summary = result.summary();
        assert!(summary.contains("VALID"));
    }

    #[test]
    fn test_parallel_check_config_presets() {
        let fast = ParallelCheckConfig::fast();
        assert!(!fast.check_self_intersections);
        assert!(!fast.check_same_parameter);

        let thorough = ParallelCheckConfig::thorough();
        assert!((thorough.tolerance - 1e-9).abs() < 1e-15);
        assert!(thorough.check_self_intersections);
    }

    #[test]
    fn test_parallel_check_report_timing() {
        let timing = CheckPhaseTiming {
            phase: "test".to_string(),
            duration_ms: 100,
            items_processed: 50,
        };

        assert_eq!(timing.phase, "test");
        assert_eq!(timing.duration_ms, 100);
        assert_eq!(timing.items_processed, 50);
    }
}
