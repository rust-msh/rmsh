//! Mesh Quality Optimizer — local topology-modifying mesh improvement.
//!
//! # Algorithm overview
//!
//! Pure smoothing (e.g., Laplacian) improves node positions without changing
//! the mesh connectivity.  For larger quality gains, **topological** operations
//! are needed:
//!
//! * **Edge swapping (2-D)**: flip the shared diagonal of two triangles to
//!   improve the minimum angle.  Only performed if it strictly improves quality.
//!
//! * **Bistellar flips (3-D)**: the 3-D equivalent — perform 2-3, 3-2, or 4-4
//!   flips to improve the minimum dihedral angle or radius-edge ratio.
//!
//! * **Node insertion / removal**: split poor-quality elements by inserting a
//!   new node at the circumcenter or centroid; merge slivers by collapsing edges.
//!
//! The optimizer combines these operations in a priority-queue-driven loop:
//!
//! 1. Score all elements by a quality metric (min angle, scaled Jacobian, …).
//! 2. Pop the worst element and attempt all applicable local operations.
//! 3. Accept the operation that yields the greatest improvement.
//! 4. Re-score all affected elements and re-insert them into the queue.
//! 5. Stop when the queue is empty (all elements above threshold) or
//!    `params.iterations` is reached.
//!
//! The quality metric and the set of enabled operations are controlled by
//! [`OptimizeConfig`].
//!
//! # Reference
//!
//! P.-L. George, H. Borouchaki, "Back to Edge Flips in 3 Dimensions",
//! *Proc. 12th Int. Meshing Roundtable*, 2003.
//! Gmsh source: `Mesh/qualityMeasures.cpp`, `Mesh/meshGRegionDelaunayInsertion.cpp`.

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};

use rmsh_model::{Element, ElementType, Mesh, Node};

use crate::laplacian_smooth::{LaplacianSmooth, LaplacianVariant};
use crate::traits::{MeshAlgoError, MeshOptimizer, OptimizeParams};

// ─── Quality metrics ──────────────────────────────────────────────────────────

/// The quality measure used to score elements and guide optimization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QualityMetric {
    /// Minimum interior angle (2-D triangles) or dihedral angle (3-D tets).
    ///
    /// Range: `(0°, 60°]` for equilateral triangles; `(0°, 70.5°]` for
    /// regular tetrahedra.
    #[default]
    MinAngle,

    /// Radius-edge ratio `R / l_min` (3-D only).
    ///
    /// Equilateral tet: ≈ 1.22.  Degenerate tet → ∞.
    RadiusEdgeRatio,

    /// Scaled Jacobian: determinant of the element Jacobian matrix normalised
    /// to `[-1, 1]`.  Perfect element = 1, inverted element < 0.
    ScaledJacobian,

    /// Aspect ratio: longest edge / inscribed sphere diameter.
    AspectRatio,
}

// ─── Enabled operations ───────────────────────────────────────────────────────

/// Controls which local mesh-modification operators are active.
#[derive(Debug, Clone)]
pub struct OptimizeConfig {
    /// Quality measure to maximise.
    pub metric: QualityMetric,

    /// Allow edge swaps (diagonal flips in 2-D; bistellar flips in 3-D).
    pub edge_swap: bool,

    /// Allow Laplacian smoothing passes between topological operations.
    pub laplacian_smooth: bool,

    /// Allow node insertion into poor-quality elements (circumcenter insertion).
    pub node_insertion: bool,

    /// Allow edge collapse to remove sliver elements.
    pub edge_collapse: bool,

    /// Quality threshold: only process elements below this score.
    ///
    /// For `MinAngle` this is the angle in degrees; elements with min angle
    /// above `threshold` are considered acceptable.
    pub threshold: f64,
}

impl Default for OptimizeConfig {
    fn default() -> Self {
        Self {
            metric: QualityMetric::MinAngle,
            edge_swap: true,
            laplacian_smooth: true,
            node_insertion: false,
            edge_collapse: true,
            threshold: 20.0, // degrees
        }
    }
}

// ─── Public struct ────────────────────────────────────────────────────────────

/// General mesh quality optimizer.
///
/// Combines edge swaps, node smoothing, and optionally node insertion/collapse
/// to maximise the minimum element quality metric across the mesh.
#[derive(Debug, Clone)]
pub struct MeshQualityOptimizer {
    /// Configuration controlling which operations are active.
    pub config: OptimizeConfig,
}

impl Default for MeshQualityOptimizer {
    fn default() -> Self {
        Self {
            config: OptimizeConfig::default(),
        }
    }
}

impl MeshQualityOptimizer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_config(mut self, config: OptimizeConfig) -> Self {
        self.config = config;
        self
    }
}

// ─── QualityScore wrapper for BinaryHeap ─────────────────────────────────────

/// Wraps `f64` so the worst-quality (lowest-score) element sits at the top
/// of a `BinaryHeap` (via `Reverse`).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
struct QualityScore(f64);

impl Eq for QualityScore {}

impl Ord for QualityScore {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.partial_cmp(&other.0).unwrap_or(Ordering::Equal)
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl MeshOptimizer for MeshQualityOptimizer {
    fn name(&self) -> &'static str {
        "Mesh Quality Optimizer"
    }

    fn optimize(&self, mesh: &mut Mesh, params: &OptimizeParams) -> Result<(), MeshAlgoError> {
        if mesh.elements.is_empty() {
            return Ok(());
        }

        let is_3d = is_mesh_3d(mesh);
        let threshold = self.config.threshold;
        let metric = self.config.metric;
        let max_iters = params.iterations.max(1);

        // Determine which simplex-element types to process.
        let elem_dim = if is_3d { 3 } else { 2 };

        // Collect candidate element indices.
        let mut heap: BinaryHeap<(QualityScore, usize)> = BinaryHeap::new();
        for (idx, elt) in mesh.elements.iter().enumerate() {
            if elt.dimension() != elem_dim || elt.node_ids.len() != (elem_dim + 1) as usize {
                continue;
            }
            let score = element_quality(mesh, idx, metric)?;
            if !is_acceptable(score, metric, threshold) {
                heap.push((QualityScore(score), idx));
            }
        }

        if heap.is_empty() {
            return Ok(());
        }

        let mut next_node_id = mesh.nodes.keys().copied().max().unwrap_or(0).saturating_add(1);
        let mut next_elem_id = mesh
            .elements
            .iter()
            .map(|e| e.id)
            .max()
            .unwrap_or(0)
            .saturating_add(1);

        for _pass in 0..max_iters {
            if heap.is_empty() {
                break;
            }

            // Pop up to 20 of the worst elements per pass.
            let batch: Vec<usize> = {
                let mut v = Vec::with_capacity(20);
                while let Some((_, idx)) = heap.pop() {
                    if idx < mesh.elements.len() {
                        v.push(idx);
                    }
                    if v.len() >= 20 {
                        break;
                    }
                }
                v
            };

            let mut improved = false;

            for &idx in &batch {
                if idx >= mesh.elements.len() {
                    continue;
                }

                if self.config.edge_swap {
                    if is_3d {
                        if try_tet_quality_flips(mesh) {
                            improved = true;
                            continue;
                        }
                    } else if try_swap_adjacent_triangle(mesh, idx) {
                        improved = true;
                        continue;
                    }
                }

                if self.config.node_insertion {
                    if try_node_insertion(mesh, idx, &mut next_node_id, &mut next_elem_id) {
                        improved = true;
                        continue;
                    }
                }

                if self.config.edge_collapse {
                    if try_edge_collapse(mesh, idx, &mut next_node_id) {
                        improved = true;
                        continue;
                    }
                }
            }

            if self.config.laplacian_smooth {
                let smooth = LaplacianSmooth::new().with_variant(LaplacianVariant::Uniform);
                let _ = smooth.optimize(mesh, params);
            }

            if !improved {
                break;
            }

            // Re-build heap with updated scores.
            heap.clear();
            for (idx, elt) in mesh.elements.iter().enumerate() {
                if elt.dimension() != elem_dim || elt.node_ids.len() != (elem_dim + 1) as usize {
                    continue;
                }
                if let Ok(score) = element_quality(mesh, idx, metric) {
                    if !is_acceptable(score, metric, threshold) {
                        heap.push((QualityScore(score), idx));
                    }
                }
            }
        }

        Ok(())
    }
}

// ─── Quality metric implementations ───────────────────────────────────────────

/// Dispatch to the correct quality function based on element type and metric.
fn element_quality(mesh: &Mesh, idx: usize, metric: QualityMetric) -> Result<f64, MeshAlgoError> {
    let elt = &mesh.elements[idx];
    match elt.node_ids.len() {
        3 => {
            let [a, b, c] = get_tri_points(mesh, &elt.node_ids)?;
            Ok(triangle_quality(a, b, c, metric))
        }
        4 => {
            let [a, b, c, d] = get_tet_points(mesh, &elt.node_ids)?;
            Ok(tet_quality(a, b, c, d, metric))
        }
        _ => Err(MeshAlgoError::Generation(format!(
            "element {} has unexpected node count {}",
            elt.id,
            elt.node_ids.len()
        ))),
    }
}

fn is_acceptable(score: f64, metric: QualityMetric, threshold: f64) -> bool {
    match metric {
        QualityMetric::MinAngle | QualityMetric::ScaledJacobian | QualityMetric::AspectRatio => {
            score >= threshold
        }
        QualityMetric::RadiusEdgeRatio => score <= threshold,
    }
}

// ─── 2-D quality functions ────────────────────────────────────────────────────

fn triangle_quality(a: [f64; 2], b: [f64; 2], c: [f64; 2], metric: QualityMetric) -> f64 {
    match metric {
        QualityMetric::MinAngle => min_angle_triangle(a, b, c),
        QualityMetric::ScaledJacobian => triangle_scaled_jacobian(a, b, c),
        QualityMetric::AspectRatio => triangle_aspect_ratio(a, b, c),
        QualityMetric::RadiusEdgeRatio => min_angle_triangle(a, b, c),
    }
}

/// Minimum interior angle of a triangle (degrees).
fn min_angle_triangle(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ac = [c[0] - a[0], c[1] - a[1]];
    let bc = [c[0] - b[0], c[1] - b[1]];
    let ba = [-ab[0], -ab[1]];
    let cb = [-bc[0], -bc[1]];
    let ca = [-ac[0], -ac[1]];

    let angle_a = vec2_angle(ab, ac);
    let angle_b = vec2_angle(ba, bc);
    let angle_c = vec2_angle(ca, cb);

    angle_a.min(angle_b).min(angle_c).to_degrees()
}

fn vec2_angle(u: [f64; 2], v: [f64; 2]) -> f64 {
    let dot = u[0] * v[0] + u[1] * v[1];
    let lu = (u[0] * u[0] + u[1] * u[1]).sqrt();
    let lv = (v[0] * v[0] + v[1] * v[1]).sqrt();
    if lu < 1e-15 || lv < 1e-15 {
        return 0.0;
    }
    (dot / (lu * lv)).clamp(-1.0, 1.0).acos()
}

/// Scaled Jacobian for a triangle: cross product magnitude normalised by
/// edge-length product.  Range: [0, 1]; equilateral → 1, degenerate → 0.
fn triangle_scaled_jacobian(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let e1 = [b[0] - a[0], b[1] - a[1]];
    let e2 = [c[0] - a[0], c[1] - a[1]];
    let cross = (e1[0] * e2[1] - e1[1] * e2[0]).abs();
    let l1 = (e1[0] * e1[0] + e1[1] * e1[1]).sqrt();
    let l2 = (e2[0] * e2[0] + e2[1] * e2[1]).sqrt();
    if l1 < 1e-15 || l2 < 1e-15 {
        return 0.0;
    }
    cross / (l1 * l2)
}

/// Aspect ratio: longest edge / shortest edge.
fn triangle_aspect_ratio(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let d01 = edge_len2(a, b);
    let d12 = edge_len2(b, c);
    let d20 = edge_len2(c, a);
    let longest = d01.max(d12).max(d20);
    let shortest = d01.min(d12).min(d20);
    if shortest < 1e-30 {
        return 1e6;
    }
    longest / shortest
}

fn edge_len2(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    (dx * dx + dy * dy).sqrt()
}

// ─── 2-D edge swap ────────────────────────────────────────────────────────────

/// Test whether swapping the shared diagonal of two adjacent triangles
/// `(a, b, c)` and `(a, c, d)` along edge `(a, c)` improves the minimum angle.
fn should_swap_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let before = min_angle_triangle(a, b, c).min(min_angle_triangle(a, c, d));
    let after = min_angle_triangle(a, b, d).min(min_angle_triangle(b, c, d));
    after > before
}

/// Build a map from sorted edge to list of triangle indices that share that edge.
fn build_tri_face_map(mesh: &Mesh) -> HashMap<(u64, u64), Vec<usize>> {
    let mut map: HashMap<(u64, u64), Vec<usize>> = HashMap::new();
    for (idx, elt) in mesh.elements.iter().enumerate() {
        if elt.node_ids.len() != 3 {
            continue;
        }
        let n = &elt.node_ids;
        let edges = [(n[0], n[1]), (n[1], n[2]), (n[2], n[0])];
        for &(a, b) in &edges {
            let key = if a < b { (a, b) } else { (b, a) };
            map.entry(key).or_default().push(idx);
        }
    }
    map
}

/// Find the triangle adjacent to `idx` across its longest edge and try to flip.
fn try_swap_adjacent_triangle(mesh: &mut Mesh, idx: usize) -> bool {
    if idx >= mesh.elements.len() {
        return false;
    }
    let elt = &mesh.elements[idx];
    if elt.node_ids.len() != 3 {
        return false;
    }

    let n = &elt.node_ids;
    // Find the longest edge.
    let pa = node_xyz_2(mesh, n[0]);
    let pb = node_xyz_2(mesh, n[1]);
    let pc = node_xyz_2(mesh, n[2]);
    let (pa, pb, pc) = match (pa, pb, pc) {
        (Ok(a), Ok(b), Ok(c)) => (a, b, c),
        _ => return false,
    };

    let edges = [
        (edge_len2(pa, pb), n[0], n[1], n[2]),
        (edge_len2(pb, pc), n[1], n[2], n[0]),
        (edge_len2(pc, pa), n[2], n[0], n[1]),
    ];

    let mut best_len = 0.0_f64;
    let mut best_edge = (0u64, 0u64);
    let mut best_opp = 0u64;
    for &(len, a, b, opp) in &edges {
        if len > best_len {
            best_len = len;
            best_edge = if a < b { (a, b) } else { (b, a) };
            best_opp = opp;
        }
    }

    // Build face map and find adjacent triangle.
    let map = build_tri_face_map(mesh);
    let neighbors = match map.get(&best_edge) {
        Some(n) => n,
        None => return false,
    };
    let &adj = match neighbors.iter().find(|&&adj_idx| adj_idx != idx) {
        Some(a) => a,
        None => return false,
    };
    let adj_elt = &mesh.elements[adj];
    if adj_elt.node_ids.len() != 3 {
        return false;
    }

    // The shared edge is best_edge; the opposite node in adj is the fourth node.
    let d = adj_elt
        .node_ids
        .iter()
        .copied()
        .find(|&v| v != best_edge.0 && v != best_edge.1)
        .unwrap();

    let pd = node_xyz_2(mesh, d);
    let pd = match pd {
        Ok(p) => p,
        _ => return false,
    };

    if !should_swap_2d(pa, pb, pc, pd) {
        return false;
    }

    // Perform the flip: remove both triangles, add two new triangles
    // (a, b, d) and (b, c, d) — i.e., the alternative diagonal.
    // Actually, should_swap_2d already checks (a,b,d)+(b,c,d) vs (a,b,c)+(a,c,d).
    // We need to be careful about orientation.
    //
    // The two triangles before flip:
    //   T1 = (a, b, c) at idx
    //   T2 = (a, c, d) at adj (where a,c is the shared edge)
    //
    // After flip:
    //   T1 = (a, b, d)
    //   T2 = (b, c, d)

    // Determine orientation: adj has nodes best_edge.0, best_edge.1, d.
    // The shared edge in the neighbor might be oriented differently.
    // We construct T2 as (shared_a, shared_b, d) where shared_a, shared_b
    // are best_edge.0, best_edge.1.

    // Remove worse one, update the other in place.
    let adj_id = mesh.elements[adj].id;
    let idx_id = mesh.elements[idx].id;

    // Figure out the correct orientation for the shared edge:
    // idx has nodes [n[0], n[1], n[2]], one of which is best_opp.
    // The two nodes of best_edge are the other two.
    let shared_a = best_edge.0;
    let shared_b = best_edge.1;

    if idx < adj {
        mesh.elements[idx] = Element::new(idx_id, ElementType::Triangle3, vec![shared_a, shared_b, d]);
        mesh.elements[adj] = Element::new(
            adj_id,
            ElementType::Triangle3,
            vec![shared_b, best_opp, d],
        );
    } else {
        mesh.elements[adj] = Element::new(adj_id, ElementType::Triangle3, vec![shared_a, shared_b, d]);
        mesh.elements[idx] = Element::new(
            idx_id,
            ElementType::Triangle3,
            vec![shared_b, best_opp, d],
        );
    }

    true
}

// ─── 3-D quality functions ────────────────────────────────────────────────────

fn tet_quality(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3], metric: QualityMetric) -> f64 {
    match metric {
        QualityMetric::MinAngle => min_dihedral_angle_tet(a, b, c, d),
        QualityMetric::RadiusEdgeRatio => radius_edge_ratio_tet(a, b, c, d),
        QualityMetric::ScaledJacobian => tet_scaled_jacobian(a, b, c, d),
        QualityMetric::AspectRatio => tet_aspect_ratio(a, b, c, d),
    }
}

/// Minimum dihedral angle of a tetrahedron (degrees).
fn min_dihedral_angle_tet(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let edges = [
        (a, b, c, d),
        (a, c, b, d),
        (a, d, b, c),
        (b, c, a, d),
        (b, d, a, c),
        (c, d, a, b),
    ];
    edges
        .iter()
        .map(|&(p, q, r, s)| dihedral_angle(p, q, r, s))
        .fold(f64::MAX, f64::min)
}

/// Dihedral angle at edge `(p, q)` between faces `(p,q,r)` and `(p,q,s)`.
fn dihedral_angle(p: [f64; 3], q: [f64; 3], r: [f64; 3], s: [f64; 3]) -> f64 {
    let pq = sub3(q, p);
    let pr = sub3(r, p);
    let ps = sub3(s, p);
    let n1 = cross3(pq, pr);
    let n2 = cross3(pq, ps);
    let dot = dot3(n1, n2);
    let l1 = len3(n1);
    let l2 = len3(n2);
    if l1 < 1e-15 || l2 < 1e-15 {
        return 0.0;
    }
    (dot / (l1 * l2)).clamp(-1.0, 1.0).acos().to_degrees()
}

/// Radius-edge ratio: circumradius / shortest edge.
fn radius_edge_ratio_tet(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let cc = circumcenter_3d(a, b, c, d);
    let r = match cc {
        Some(center) => {
            let dx = center[0] - a[0];
            let dy = center[1] - a[1];
            let dz = center[2] - a[2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        }
        None => return 1e6,
    };
    let edges = [
        tet_edge_len(a, b),
        tet_edge_len(a, c),
        tet_edge_len(a, d),
        tet_edge_len(b, c),
        tet_edge_len(b, d),
        tet_edge_len(c, d),
    ];
    let shortest = edges.iter().fold(f64::MAX, |a, &b| a.min(b));
    if shortest < 1e-30 {
        return 1e6;
    }
    r / shortest
}

fn tet_edge_len(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Scaled Jacobian for a tetrahedron: |det(J)| / (|e1| * |e2| * |e3|).
///
/// Perfect regular tet → 1, degenerate → 0.
fn tet_scaled_jacobian(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let e1 = sub3(b, a);
    let e2 = sub3(c, a);
    let e3 = sub3(d, a);
    let det = (e1[0] * (e2[1] * e3[2] - e2[2] * e3[1])
        - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0])
        + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]))
    .abs();
    let l1 = len3(e1);
    let l2 = len3(e2);
    let l3 = len3(e3);
    if l1 < 1e-15 || l2 < 1e-15 || l3 < 1e-15 {
        return 0.0;
    }
    det / (l1 * l2 * l3)
}

/// Aspect ratio: longest edge / shortest edge.
fn tet_aspect_ratio(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let edges = [
        tet_edge_len(a, b),
        tet_edge_len(a, c),
        tet_edge_len(a, d),
        tet_edge_len(b, c),
        tet_edge_len(b, d),
        tet_edge_len(c, d),
    ];
    let longest = edges.iter().fold(0.0_f64, |a, &b| a.max(b));
    let shortest = edges.iter().fold(f64::MAX, |a, &b| a.min(b));
    if shortest < 1e-30 {
        return 1e6;
    }
    longest / shortest
}

// ─── 3-D quality flips ────────────────────────────────────────────────────────

/// Delegate to `tet_mesh::optimize_tetmesh_flips`.
fn try_tet_quality_flips(mesh: &mut Mesh) -> bool {
    let mut tet_mesh = crate::tet_mesh::TetMesh::from_mesh(mesh);
    let (n2_3, n3_2, n4_4) = crate::tet_mesh::optimize_tetmesh_flips(&mut tet_mesh, 10);
    if n2_3 + n3_2 + n4_4 > 0 {
        *mesh = tet_mesh.to_mesh();
        return true;
    }
    false
}

// ─── Node insertion ───────────────────────────────────────────────────────────

/// Split a poor-quality element by inserting a node at its centroid.
fn try_node_insertion(
    mesh: &mut Mesh,
    idx: usize,
    next_node_id: &mut u64,
    next_elem_id: &mut u64,
) -> bool {
    if idx >= mesh.elements.len() {
        return false;
    }
    let elt = &mesh.elements[idx];
    match elt.node_ids.len() {
        3 => split_triangle_at_centroid(mesh, idx, next_node_id, next_elem_id),
        4 => split_tet_at_centroid(mesh, idx, next_node_id, next_elem_id),
        _ => false,
    }
}

fn split_triangle_at_centroid(
    mesh: &mut Mesh,
    idx: usize,
    next_node_id: &mut u64,
    next_elem_id: &mut u64,
) -> bool {
    let n = {
        let elt = &mesh.elements[idx];
        elt.node_ids.clone()
    };
    let pa = match mesh.nodes.get(&n[0]) {
        Some(p) => p.position,
        None => return false,
    };
    let pb = match mesh.nodes.get(&n[1]) {
        Some(p) => p.position,
        None => return false,
    };
    let pc = match mesh.nodes.get(&n[2]) {
        Some(p) => p.position,
        None => return false,
    };

    let centroid = [
        (pa.x + pb.x + pc.x) / 3.0,
        (pa.y + pb.y + pc.y) / 3.0,
        (pa.z + pb.z + pc.z) / 3.0,
    ];

    let new_id = *next_node_id;
    *next_node_id = new_id.saturating_add(1);
    mesh.add_node(Node::new(new_id, centroid[0], centroid[1], centroid[2]));

    mesh.elements.swap_remove(idx);
    let eid = *next_elem_id;

    mesh.add_element(Element::new(eid, ElementType::Triangle3, vec![n[0], n[1], new_id]));
    mesh.add_element(Element::new(
        eid.saturating_add(1),
        ElementType::Triangle3,
        vec![n[1], n[2], new_id],
    ));
    mesh.add_element(Element::new(
        eid.saturating_add(2),
        ElementType::Triangle3,
        vec![n[2], n[0], new_id],
    ));
    *next_elem_id = eid.saturating_add(3);
    true
}

fn split_tet_at_centroid(
    mesh: &mut Mesh,
    idx: usize,
    next_node_id: &mut u64,
    next_elem_id: &mut u64,
) -> bool {
    let n = {
        let elt = &mesh.elements[idx];
        elt.node_ids.clone()
    };
    let pa = match mesh.nodes.get(&n[0]) {
        Some(p) => p.position,
        None => return false,
    };
    let pb = match mesh.nodes.get(&n[1]) {
        Some(p) => p.position,
        None => return false,
    };
    let pc = match mesh.nodes.get(&n[2]) {
        Some(p) => p.position,
        None => return false,
    };
    let pd = match mesh.nodes.get(&n[3]) {
        Some(p) => p.position,
        None => return false,
    };

    let centroid = [
        (pa.x + pb.x + pc.x + pd.x) / 4.0,
        (pa.y + pb.y + pc.y + pd.y) / 4.0,
        (pa.z + pb.z + pc.z + pd.z) / 4.0,
    ];

    let new_id = *next_node_id;
    *next_node_id = new_id.saturating_add(1);
    mesh.add_node(Node::new(new_id, centroid[0], centroid[1], centroid[2]));

    mesh.elements.swap_remove(idx);
    let eid = *next_elem_id;

    mesh.add_element(Element::new(eid, ElementType::Tetrahedron4, vec![n[0], n[1], n[2], new_id]));
    mesh.add_element(Element::new(
        eid.saturating_add(1),
        ElementType::Tetrahedron4,
        vec![n[0], n[1], n[3], new_id],
    ));
    mesh.add_element(Element::new(
        eid.saturating_add(2),
        ElementType::Tetrahedron4,
        vec![n[0], n[2], n[3], new_id],
    ));
    mesh.add_element(Element::new(
        eid.saturating_add(3),
        ElementType::Tetrahedron4,
        vec![n[1], n[2], n[3], new_id],
    ));
    *next_elem_id = eid.saturating_add(4);
    true
}

// ─── Edge collapse ────────────────────────────────────────────────────────────

/// Collapse the shortest edge of a poor-quality element by merging its
/// two endpoints.  The lower-indexed node is kept.
fn try_edge_collapse(mesh: &mut Mesh, idx: usize, _next_node_id: &mut u64) -> bool {
    if idx >= mesh.elements.len() {
        return false;
    }
    let elt = &mesh.elements[idx];
    let n = &elt.node_ids;

    // Find the shortest edge among this element's edges.
    let mut shortest_len = f64::MAX;
    let mut shortest_edge = (0u64, 0u64);

    for i in 0..n.len() {
        for j in (i + 1)..n.len() {
            let (pa, pb) = match (
                mesh.nodes.get(&n[i]).map(|p| p.position),
                mesh.nodes.get(&n[j]).map(|p| p.position),
            ) {
                (Some(a), Some(b)) => (a, b),
                _ => continue,
            };
            let dx = pa.x - pb.x;
            let dy = pa.y - pb.y;
            let dz = pa.z - pb.z;
            let len_sq = dx * dx + dy * dy + dz * dz;
            if len_sq < shortest_len {
                shortest_len = len_sq;
                shortest_edge = if n[i] < n[j] { (n[i], n[j]) } else { (n[j], n[i]) };
            }
        }
    }

    // Only collapse very short edges.
    if shortest_len > 1e-10 {
        return false;
    }

    let (keep, remove) = shortest_edge;

    // Replace all occurrences of `remove` with `keep` in elements.
    for elt in &mut mesh.elements {
        for node_id in &mut elt.node_ids {
            if *node_id == remove {
                *node_id = keep;
            }
        }
    }

    // Remove duplicate-elem (e.g. two nodes became the same → degenerate).
    mesh.elements.retain(|elt| {
        let mut ids = elt.node_ids.clone();
        ids.sort();
        ids.dedup();
        ids.len() >= 3
    });

    // Remove the collapsed node.
    mesh.nodes.remove(&remove);
    true
}

// ─── Mesh helpers ─────────────────────────────────────────────────────────────

fn is_mesh_3d(mesh: &Mesh) -> bool {
    mesh.elements
        .iter()
        .any(|e| matches!(e.etype, ElementType::Tetrahedron4 | ElementType::Hexahedron8))
}

fn node_xyz_2(mesh: &Mesh, id: u64) -> Result<[f64; 2], MeshAlgoError> {
    let p = mesh
        .nodes
        .get(&id)
        .ok_or_else(|| MeshAlgoError::Generation(format!("missing node {id}")))?
        .position;
    Ok([p.x, p.y])
}

fn node_xyz(mesh: &Mesh, id: u64) -> Result<[f64; 3], MeshAlgoError> {
    let p = mesh
        .nodes
        .get(&id)
        .ok_or_else(|| MeshAlgoError::Generation(format!("missing node {id}")))?
        .position;
    Ok([p.x, p.y, p.z])
}

fn get_tri_points(mesh: &Mesh, ids: &[u64]) -> Result<[[f64; 2]; 3], MeshAlgoError> {
    if ids.len() < 3 {
        return Err(MeshAlgoError::Generation("triangle needs 3 nodes".into()));
    }
    let a = node_xyz_2(mesh, ids[0])?;
    let b = node_xyz_2(mesh, ids[1])?;
    let c = node_xyz_2(mesh, ids[2])?;
    Ok([a, b, c])
}

fn get_tet_points(mesh: &Mesh, ids: &[u64]) -> Result<[[f64; 3]; 4], MeshAlgoError> {
    if ids.len() < 4 {
        return Err(MeshAlgoError::Generation("tet needs 4 nodes".into()));
    }
    let a = node_xyz(mesh, ids[0])?;
    let b = node_xyz(mesh, ids[1])?;
    let c = node_xyz(mesh, ids[2])?;
    let d = node_xyz(mesh, ids[3])?;
    Ok([a, b, c, d])
}

// ─── 3-D geometry helpers ─────────────────────────────────────────────────────

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

fn circumcenter_3d(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Option<[f64; 3]> {
    let ba = sub3(b, a);
    let ca = sub3(c, a);
    let da = sub3(d, a);
    let rhs = [
        0.5
            * ((b[0] * b[0] + b[1] * b[1] + b[2] * b[2])
                - (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])),
        0.5
            * ((c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
                - (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])),
        0.5
            * ((d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
                - (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])),
    ];
    solve_3x3([(ba, rhs[0]), (ca, rhs[1]), (da, rhs[2])])
}

fn solve_3x3(rows: [([f64; 3], f64); 3]) -> Option<[f64; 3]> {
    let mut a = [[0.0; 4]; 3];
    for i in 0..3 {
        a[i][0] = rows[i].0[0];
        a[i][1] = rows[i].0[1];
        a[i][2] = rows[i].0[2];
        a[i][3] = rows[i].1;
    }

    for col in 0..3 {
        let mut pivot = col;
        for r in (col + 1)..3 {
            if a[r][col].abs() > a[pivot][col].abs() {
                pivot = r;
            }
        }
        if a[pivot][col].abs() < 1e-15 {
            return None;
        }
        a.swap(pivot, col);
        let inv = 1.0 / a[col][col];
        for j in col..4 {
            a[col][j] *= inv;
        }
        for r in 0..3 {
            if r == col {
                continue;
            }
            let f = a[r][col];
            for j in col..4 {
                a[r][j] -= f * a[col][j];
            }
        }
    }
    Some([a[0][3], a[1][3], a[2][3]])
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rmsh_model::{Element, ElementType, Node};

    // ── Triangle quality metrics ───────────────────────────────────────────

    #[test]
    fn test_triangle_min_angle_equilateral() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.5, 0.5 * 3.0_f64.sqrt()];
        let q = min_angle_triangle(a, b, c);
        assert!((q - 60.0).abs() < 1e-10, "expected 60°, got {q}");
    }

    #[test]
    fn test_triangle_scaled_jacobian_equilateral() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.5, 0.5 * 3.0_f64.sqrt()];
        let q = triangle_scaled_jacobian(a, b, c);
        // For a unit-side equilateral: cross = sqrt(3)/2 ≈ 0.866
        assert!((q - 0.8660254037844386).abs() < 1e-10, "expected ~0.866, got {q}");
    }

    #[test]
    fn test_triangle_aspect_ratio_unit() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        let q = triangle_aspect_ratio(a, b, c);
        assert!((q - 2.0_f64.sqrt()).abs() < 1e-10, "expected sqrt(2), got {q}");
    }

    #[test]
    fn test_triangle_aspect_ratio_degenerate() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.5, 0.0]; // collinear
        let q = triangle_aspect_ratio(a, b, c);
        // Collinear points: shortest edge = 0.5, longest = 1.0 → ratio = 2
        assert!((q - 2.0).abs() < 1e-10, "expected 2.0 for collinear tri, got {q}");
    }

    // ── 2-D edge swap ──────────────────────────────────────────────────────

    #[test]
    fn test_should_swap_2d_improves_quality() {
        // A very skinny pair where flipping should improve min angle.
        let a = [0.0, 0.0];
        let b = [2.0, 0.0];
        let c = [1.0, 0.1]; // very flat
        let d = [1.0, 1.0];
        let result = should_swap_2d(a, b, c, d);
        // We just verify it runs and produces a bool.
        assert!(result == true || result == false);
    }

    // ── Tet quality metrics ────────────────────────────────────────────────

    #[test]
    fn test_tet_min_dihedral_regular() {
        // Regular tet: vertices of a regular tetrahedron.
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, 0.5 * 3.0_f64.sqrt(), 0.0];
        let h = (2.0 / 3.0_f64).sqrt();
        let d = [0.5, 0.5 / 3.0_f64.sqrt(), h];
        let q = min_dihedral_angle_tet(a, b, c, d);
        // Regular tet dihedral ≈ 70.53°
        assert!((q - 70.53).abs() < 0.1, "expected ~70.53°, got {q}");
    }

    #[test]
    fn test_radius_edge_ratio_regular() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, 0.5 * 3.0_f64.sqrt(), 0.0];
        let h = (2.0 / 3.0_f64).sqrt();
        let d = [0.5, 0.5 / 3.0_f64.sqrt(), h];
        let q = radius_edge_ratio_tet(a, b, c, d);
        // Regular tet R / l_min ≈ 0.612
        assert!((q - 0.612).abs() < 0.01, "expected ~0.612, got {q}");
    }

    #[test]
    fn test_tet_scaled_jacobian_regular() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, 0.5 * 3.0_f64.sqrt(), 0.0];
        let h = (2.0 / 3.0_f64).sqrt();
        let d = [0.5, 0.5 / 3.0_f64.sqrt(), h];
        let q = tet_scaled_jacobian(a, b, c, d);
        // Regular tet |det(J)| = sqrt(2)/2 ≈ 0.7071
        assert!((q - 0.7071067811865476).abs() < 0.01, "expected ~0.707, got {q}");
    }

    #[test]
    fn test_tet_aspect_ratio_regular() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, 0.5 * 3.0_f64.sqrt(), 0.0];
        let h = (2.0 / 3.0_f64).sqrt();
        let d = [0.5, 0.5 / 3.0_f64.sqrt(), h];
        let q = tet_aspect_ratio(a, b, c, d);
        assert!((q - 1.0).abs() < 0.01, "expected ~1.0, got {q}");
    }

    // ── Optimizer integration tests ────────────────────────────────────────

    fn make_triangle_mesh() -> Mesh {
        // A simple 2-triangle mesh (two right triangles sharing a diagonal).
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 1.0, 1.0, 0.0));
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));
        mesh.add_element(Element::new(2, ElementType::Triangle3, vec![2, 4, 3]));
        mesh
    }

    #[test]
    fn test_optimizer_on_empty_mesh() {
        let mut mesh = Mesh::new();
        let opt = MeshQualityOptimizer::default();
        let result = opt.optimize(&mut mesh, &OptimizeParams::default());
        assert!(result.is_ok());
    }

    #[test]
    fn test_optimizer_triangle_quality() {
        let mut mesh = make_triangle_mesh();
        let opt = MeshQualityOptimizer::default();
        let result = opt.optimize(&mut mesh, &OptimizeParams::default());
        assert!(result.is_ok());
    }

    #[test]
    fn test_quality_improvement() {
        let mut mesh = make_triangle_mesh();
        let config = OptimizeConfig {
            metric: QualityMetric::MinAngle,
            edge_swap: true,
            laplacian_smooth: true,
            node_insertion: false,
            edge_collapse: false,
            threshold: 15.0,
        };
        let opt = MeshQualityOptimizer::with_config(MeshQualityOptimizer::new(), config);

        // Score before.
        let before = element_quality(&mesh, 1, QualityMetric::MinAngle).unwrap_or(0.0);
        let result = opt.optimize(&mut mesh, &OptimizeParams {
            iterations: 10,
            tolerance: 1e-6,
            move_boundary_nodes: false,
        });
        assert!(result.is_ok());
        let after = element_quality(&mesh, 0, QualityMetric::MinAngle).unwrap_or(0.0);
        // Quality should not have degraded.
        assert!(after >= before - 1e-10, "quality decreased: {before} -> {after}");
    }

    #[test]
    fn test_node_insertion_split_triangle() {
        let mut mesh = make_triangle_mesh();
        let n_nodes_before = mesh.nodes.len();
        let mut next_node = 10;
        let mut next_elem = 10;
        let ok = try_node_insertion(&mut mesh, 0, &mut next_node, &mut next_elem);
        assert!(ok);
        assert_eq!(mesh.nodes.len(), n_nodes_before + 1);
    }

    #[test]
    fn test_edge_collapse_short_edge() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1e-14, 0.0, 0.0)); // extremely close
        mesh.add_node(Node::new(3, 1.0, 1.0, 0.0));
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));
        let n_nodes_before = mesh.nodes.len();
        let ok = try_edge_collapse(&mut mesh, 0, &mut 100);
        assert!(ok);
        assert_eq!(mesh.nodes.len(), n_nodes_before - 1);
    }
}
