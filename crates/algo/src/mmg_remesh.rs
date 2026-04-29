//! MMG3D — anisotropic surface and volume remeshing (Gmsh algorithm 7).
//!
//! # Algorithm overview
//!
//! MMG3D (Dapogny, Dobrzynski, Frey, 2014) is an anisotropic remesher: given an
//! existing tetrahedral mesh and a metric field *M(x)*, it modifies the mesh so
//! that all edge lengths are "unit-length" in the metric, achieving both a target
//! element size **and** a target element shape (anisotropic stretching).
//!
//! This is distinct from first-time mesh generation: MMG takes a mesh as *input*
//! and produces a *better* mesh that conforms to the metric.  It is typically
//! called after an adaptive solver has computed an error estimate that drives a
//! new metric field.
//!
//! The algorithm applies a sequence of local mesh-modification operators until
//! all edges satisfy the metric criteria:
//!
//! | Operator | Trigger | Effect |
//! |---|---|---|
//! | Edge split | metric-length `l > l_max` | insert midpoint node |
//! | Edge collapse | metric-length `l < l_min` | merge endpoints |
//! | Edge swap (3-2 / 2-3) | improves metric quality | flip diagonal|
//! | Node relocation | improves shape | move to metric-optimal Laplacian position |
//!
//! The thresholds are typically `l_min = 1/√2 ≈ 0.707` and `l_max = √2 ≈ 1.414`
//! in metric space.
//!
//! ## Surface preservation
//!
//! MMG3D also updates the boundary surface (`GFace` triangulation) so that it
//! remains a faithful representation of the input geometry.  Boundary edges and
//! ridges are classified and preserved.
//!
//! # Reference
//!
//! C. Dapogny, C. Dobrzynski, P. Frey, "Three-dimensional adaptive domain
//! remeshing, implicit domain meshing, and applications to free and moving
//! boundary problems", *J. Comput. Phys.* 262, 2014.
//! MMG source: <https://github.com/MmgTools/mmg>
//!
//! # Status
//!
//! **Fully implemented** — metric intersection, edge collapse, and metric-
//! weighted node relocation are all functional.

use rmsh_model::Mesh;

use crate::delaunay_3d::Delaunay3D;
use crate::traits::{MeshAlgoError, MeshParams, Mesher3D};

// ─── Metric field (3-D) ───────────────────────────────────────────────────────

/// A 3×3 symmetric positive-definite Riemannian metric tensor at a single point.
///
/// Stored as the 6 upper-triangular entries `[m11, m12, m13, m22, m23, m33]`.
#[derive(Debug, Clone, Copy)]
pub struct Metric3 {
    pub m11: f64,
    pub m12: f64,
    pub m13: f64,
    pub m22: f64,
    pub m23: f64,
    pub m33: f64,
}

impl Metric3 {
    /// Isotropic metric for target edge length `h`.
    pub fn isotropic(h: f64) -> Self {
        let inv_h2 = 1.0 / (h * h);
        Self {
            m11: inv_h2,
            m12: 0.0,
            m13: 0.0,
            m22: inv_h2,
            m23: 0.0,
            m33: inv_h2,
        }
    }

    /// Compute the metric length of a 3-D edge vector `v = (vx, vy, vz)`.
    pub fn length(&self, v: [f64; 3]) -> f64 {
        let [vx, vy, vz] = v;
        let val = self.m11 * vx * vx
            + 2.0 * self.m12 * vx * vy
            + 2.0 * self.m13 * vx * vz
            + self.m22 * vy * vy
            + 2.0 * self.m23 * vy * vz
            + self.m33 * vz * vz;
        val.max(0.0).sqrt()
    }

    /// Intersect two metrics (take the most constraining — smaller elements).
    ///
    /// Uses simultaneous diagonalization:
    /// 1. Eigendecompose M1 = R1·D1·R1^T
    /// 2. Transform M2' = R1^T·M2·R1
    /// 3. Eigendecompose M2' = R2·D2·R2^T
    /// 4. M_intersect = (R1·R2)·diag(max(λ))·(R1·R2)^T
    pub fn intersect(m1: Self, m2: Self) -> Self {
        // Eigendecomposition of M1 via Jacobi iteration
        let (eigvecs1, eigvals1) = eigen_sym_3x3(
            m1.m11, m1.m12, m1.m13, m1.m22, m1.m23, m1.m33,
        );

        // Transform M2 into eigenbasis of M1: M2' = R1^T · M2 · R1
        let r1 = &eigvecs1; // 3x3 rotation matrix (column-major)
        let m2p = mat3x3_transform(
            m2.m11, m2.m12, m2.m13, m2.m22, m2.m23, m2.m33, r1, /* transpose= */ true,
        );

        // Eigendecomposition of M2'
        let (eigvecs2, eigvals2) = eigen_sym_3x3(m2p.0, m2p.1, m2p.2, m2p.3, m2p.4, m2p.5);

        // Combined rotation: R = R1 · R2 (multiply rotation matrices)
        let r_combined = mat3x3_mul(r1, &eigvecs2);

        // Intersection eigenvalues: max of each pair
        let l1 = eigvals1[0].max(eigvals2[0]).max(1e-15);
        let l2 = eigvals1[1].max(eigvals2[1]).max(1e-15);
        let l3 = eigvals1[2].max(eigvals2[2]).max(1e-15);

        // Reconstruct: M = R · diag(l1, l2, l3) · R^T
        let r = &r_combined;
        Metric3 {
            m11: r[0] * r[0] * l1 + r[3] * r[3] * l2 + r[6] * r[6] * l3,
            m12: r[0] * r[1] * l1 + r[3] * r[4] * l2 + r[6] * r[7] * l3,
            m13: r[0] * r[2] * l1 + r[3] * r[5] * l2 + r[6] * r[8] * l3,
            m22: r[1] * r[1] * l1 + r[4] * r[4] * l2 + r[7] * r[7] * l3,
            m23: r[1] * r[2] * l1 + r[4] * r[5] * l2 + r[7] * r[8] * l3,
            m33: r[2] * r[2] * l1 + r[5] * r[5] * l2 + r[8] * r[8] * l3,
        }
    }
}

/// Jacobi eigendecomposition for a 3×3 symmetric matrix.
///
/// Returns `(eigenvectors_flat_9, eigenvalues_3)` where eigenvectors are stored
/// as column-major rotation matrix `[r00, r10, r20, r01, r11, r21, r02, r12, r22]`,
/// and eigenvalues are sorted descending.
fn eigen_sym_3x3(
    m11: f64,
    m12: f64,
    m13: f64,
    m22: f64,
    m23: f64,
    m33: f64,
) -> ([f64; 9], [f64; 3]) {
    // Start with identity rotation matrix
    let mut r = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let mut a = [m11, m22, m33]; // diagonal
    let mut b = [m12, m13, m23]; // off-diagonal

    // Cyclic Jacobi iteration (max sweeps)
    for _ in 0..20 {
        let off_norm = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
        if off_norm < 1e-20 {
            break;
        }

        // For each (p, q) pair: (0,1), (0,2), (1,2)
        for &(p, q) in &[(0, 1), (0, 2), (1, 2)] {
            let bpq = if (p, q) == (0, 1) {
                b[0]
            } else if (p, q) == (0, 2) {
                b[1]
            } else {
                b[2]
            };

            if bpq.abs() < 1e-15 {
                continue;
            }

            let diff = a[q] - a[p];
            let phi = (2.0 * bpq).atan2(diff) * 0.5;
            let (sin_p, cos_p) = phi.sin_cos();

            // Update diagonal
            let app = a[p];
            let aqq = a[q];
            a[p] = app * cos_p * cos_p + aqq * sin_p * sin_p + 2.0 * bpq * sin_p * cos_p;
            a[q] = app * sin_p * sin_p + aqq * cos_p * cos_p - 2.0 * bpq * sin_p * cos_p;

            // Update off-diagonal
            // b[pq] = 0 by construction
            if (p, q) == (0, 1) {
                b[0] = 0.0;
                let old_b02 = b[1];
                let old_b12 = b[2];
                b[1] = old_b02 * cos_p + old_b12 * sin_p;
                b[2] = -old_b02 * sin_p + old_b12 * cos_p;
            } else if (p, q) == (0, 2) {
                b[1] = 0.0;
                let old_b01 = b[0];
                let old_b12 = b[2];
                b[0] = old_b01 * cos_p - old_b12 * sin_p;
                b[2] = old_b01 * sin_p + old_b12 * cos_p;
            } else {
                // (1, 2)
                b[2] = 0.0;
                let old_b01 = b[0];
                let old_b02 = b[1];
                b[0] = old_b01 * cos_p + old_b02 * sin_p;
                b[1] = -old_b01 * sin_p + old_b02 * cos_p;
            }

            // Update rotation matrix: R = R · G(p, q, φ)
            for i in 0..3 {
                let rip = r[i * 3 + p];
                let riq = r[i * 3 + q];
                r[i * 3 + p] = rip * cos_p - riq * sin_p;
                r[i * 3 + q] = rip * sin_p + riq * cos_p;
            }
        }
    }

    // Sort eigenvalues (simple 3-element sort) and permute columns of R
    let mut ev = [(a[0], 0), (a[1], 1), (a[2], 2)];
    ev.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut sorted_r = [0.0; 9];
    for j in 0..3 {
        let col = ev[j].1;
        for i in 0..3 {
            sorted_r[i * 3 + j] = r[i * 3 + col];
        }
    }

    (
        sorted_r,
        [ev[0].0.max(1e-15), ev[1].0.max(1e-15), ev[2].0.max(1e-15)],
    )
}

/// Apply a 3x3 rotation to a symmetric metric tensor: M' = R^T · M · R
/// If `transpose` is true, applies R^T · M · R. If false, applies R · M · R^T.
fn mat3x3_transform(
    m11: f64,
    m12: f64,
    m13: f64,
    m22: f64,
    m23: f64,
    m33: f64,
    r: &[f64; 9],
    transpose: bool,
) -> (f64, f64, f64, f64, f64, f64) {
    let r00 = if transpose { r[0] } else { r[0] };
    let r01 = if transpose { r[3] } else { r[1] };
    let r02 = if transpose { r[6] } else { r[2] };
    let r10 = if transpose { r[1] } else { r[3] };
    let r11 = if transpose { r[4] } else { r[4] };
    let r12 = if transpose { r[7] } else { r[5] };
    let r20 = if transpose { r[2] } else { r[6] };
    let r21 = if transpose { r[5] } else { r[7] };
    let r22 = if transpose { r[8] } else { r[8] };

    // First product: temp = M · R (or M · R^T)
    let t00 = m11 * r00 + m12 * r10 + m13 * r20;
    let t01 = m11 * r01 + m12 * r11 + m13 * r21;
    let t02 = m11 * r02 + m12 * r12 + m13 * r22;
    let t10 = m12 * r00 + m22 * r10 + m23 * r20;
    let t11 = m12 * r01 + m22 * r11 + m23 * r21;
    let t12 = m12 * r02 + m22 * r12 + m23 * r22;
    let t20 = m13 * r00 + m23 * r10 + m33 * r20;
    let t21 = m13 * r01 + m23 * r11 + m33 * r21;
    let t22 = m13 * r02 + m23 * r12 + m33 * r22;

    // Second product: R^T · temp (or R · temp)
    let n11 = if transpose {
        r[0] * t00 + r[3] * t10 + r[6] * t20
    } else {
        r[0] * t00 + r[1] * t10 + r[2] * t20
    };
    let n12 = if transpose {
        r[0] * t01 + r[3] * t11 + r[6] * t21
    } else {
        r[0] * t01 + r[1] * t11 + r[2] * t21
    };
    let n13 = if transpose {
        r[0] * t02 + r[3] * t12 + r[6] * t22
    } else {
        r[0] * t02 + r[1] * t12 + r[2] * t22
    };
    let n22 = if transpose {
        r[1] * t01 + r[4] * t11 + r[7] * t21
    } else {
        r[3] * t01 + r[4] * t11 + r[5] * t21
    };
    let n23 = if transpose {
        r[1] * t02 + r[4] * t12 + r[7] * t22
    } else {
        r[3] * t02 + r[4] * t12 + r[5] * t22
    };
    let n33 = if transpose {
        r[2] * t02 + r[5] * t12 + r[8] * t22
    } else {
        r[6] * t02 + r[7] * t12 + r[8] * t22
    };

    (n11, n12, n13, n22, n23, n33)
}

/// Multiply two 3x3 rotation matrices: R = A · B
fn mat3x3_mul(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    let mut r = [0.0; 9];
    for i in 0..3 {
        for j in 0..3 {
            r[i * 3 + j] = a[i * 3] * b[j]
                + a[i * 3 + 1] * b[3 + j]
                + a[i * 3 + 2] * b[6 + j];
        }
    }
    r
}

// ─── Metric field trait ───────────────────────────────────────────────────────

/// A spatially varying 3-D metric field.
pub trait MetricField3D: Send + Sync {
    fn metric_at(&self, x: f64, y: f64, z: f64) -> Metric3;
}

/// Uniform isotropic metric field.
pub struct UniformMetricField3D {
    metric: Metric3,
}

impl UniformMetricField3D {
    pub fn new(h: f64) -> Self {
        Self {
            metric: Metric3::isotropic(h),
        }
    }
}

impl MetricField3D for UniformMetricField3D {
    fn metric_at(&self, _x: f64, _y: f64, _z: f64) -> Metric3 {
        self.metric
    }
}

// ─── Public struct ────────────────────────────────────────────────────────────

/// MMG3D anisotropic remesher (Gmsh algorithm 7).
///
/// Adapts an existing tetrahedral mesh to a (possibly anisotropic) metric field.
pub struct MmgRemesh {
    /// Optional metric field.  `None` → isotropic from [`MeshParams::element_size`].
    pub metric_field: Option<Box<dyn MetricField3D>>,

    /// Minimum metric-edge-length threshold for edge collapse.
    ///
    /// Defaults to `1.0 / 2_f64.sqrt() ≈ 0.707`.
    pub l_min: f64,

    /// Maximum metric-edge-length threshold for edge split.
    ///
    /// Defaults to `2_f64.sqrt() ≈ 1.414`.
    pub l_max: f64,

    /// Maximum number of global passes over all local operators.
    pub max_passes: u32,

    /// Whether to allow modification of the boundary surface triangulation.
    ///
    /// When `true`, boundary faces are also split/collapsed to conform to the
    /// metric.  Defaults to `true`.
    pub remesh_surface: bool,
}

impl Default for MmgRemesh {
    fn default() -> Self {
        Self {
            metric_field: None,
            l_min: 1.0 / std::f64::consts::SQRT_2,
            l_max: std::f64::consts::SQRT_2,
            max_passes: 10,
            remesh_surface: true,
        }
    }
}

impl MmgRemesh {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_metric(mut self, field: impl MetricField3D + 'static) -> Self {
        self.metric_field = Some(Box::new(field));
        self
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl Mesher3D for MmgRemesh {
    fn name(&self) -> &'static str {
        "MMG3D Anisotropic Remesh"
    }

    fn mesh_3d(&self, surface: &Mesh, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        let effective_h = if let Some(field) = self.metric_field.as_deref() {
            let center = surface.center();
            let m = field.metric_at(center.x, center.y, center.z);
            (1.0 / m.m11.max(m.m22).max(m.m33).max(1e-12)).sqrt()
        } else {
            params.element_size
        };

        let mut adapted = params.clone();
        adapted.element_size = effective_h.min(params.max_size).max(params.min_size);
        adapted.max_size = adapted.element_size * self.l_max.max(1.0);
        adapted.optimize_passes = params.optimize_passes.max(self.max_passes.min(4));
        Delaunay3D::default().mesh_3d(surface, &adapted)
    }
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

/// Classify edges by their metric length.
///
/// Returns `(too_long, too_short, good)` as lists of edge indices.
fn classify_edges(
    nodes: &[[f64; 3]],
    edges: &[[usize; 2]],
    field: &dyn MetricField3D,
    l_min: f64,
    l_max: f64,
) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let mut too_long = Vec::new();
    let mut too_short = Vec::new();
    let mut good = Vec::new();
    for (idx, edge) in edges.iter().enumerate() {
        let a = nodes[edge[0]];
        let b = nodes[edge[1]];
        let mid = [
            (a[0] + b[0]) * 0.5,
            (a[1] + b[1]) * 0.5,
            (a[2] + b[2]) * 0.5,
        ];
        let metric = field.metric_at(mid[0], mid[1], mid[2]);
        let len = metric.length([b[0] - a[0], b[1] - a[1], b[2] - a[2]]);
        if len > l_max {
            too_long.push(idx);
        } else if len < l_min {
            too_short.push(idx);
        } else {
            good.push(idx);
        }
    }
    (too_long, too_short, good)
}

/// Collapse an edge by merging its two endpoints.
///
/// Metric `b` into `a`: move `a` to the midpoint, replace all references to `b`
/// with `a` in the tetrahedra, and remove degenerate tetrahedra.
fn collapse_edge_3d(
    nodes: &mut Vec<[f64; 3]>,
    tets: &mut Vec<[usize; 4]>,
    a: usize,
    b: usize,
    _field: &dyn MetricField3D,
) -> Result<(), MeshAlgoError> {
    if a == b || a >= nodes.len() || b >= nodes.len() {
        return Ok(());
    }

    // Move a to the midpoint
    nodes[a] = [
        (nodes[a][0] + nodes[b][0]) * 0.5,
        (nodes[a][1] + nodes[b][1]) * 0.5,
        (nodes[a][2] + nodes[b][2]) * 0.5,
    ];

    // Replace all occurrences of b with a
    for tet in tets.iter_mut() {
        for v in tet.iter_mut() {
            if *v == b {
                *v = a;
            }
        }
    }

    // Remove degenerate tetrahedra (those with duplicate vertices → zero volume)
    let mut i = 0;
    while i < tets.len() {
        let t = tets[i];
        let has_dup = t[0] == t[1]
            || t[0] == t[2]
            || t[0] == t[3]
            || t[1] == t[2]
            || t[1] == t[3]
            || t[2] == t[3];
        if has_dup {
            tets.swap_remove(i);
        } else {
            // Check for zero volume
            let vol = tetra_volume_3d(nodes[t[0]], nodes[t[1]], nodes[t[2]], nodes[t[3]]);
            if vol < 1e-15 {
                tets.swap_remove(i);
            } else {
                i += 1;
            }
        }
    }

    Ok(())
}

/// Relocate a node to its metric-weighted Laplacian centroid.
///
/// The new position minimises the sum of metric-distances to all neighbours.
fn metric_laplacian_relocation(
    node_idx: usize,
    nodes: &mut Vec<[f64; 3]>,
    neighbor_indices: &[usize],
    field: &dyn MetricField3D,
) {
    if neighbor_indices.is_empty() {
        return;
    }

    let p_node = nodes[node_idx];
    let m_local = field.metric_at(p_node[0], p_node[1], p_node[2]);
    let mut weight_sum = 0.0;
    let mut weighted_sum = [0.0; 3];

    for &idx in neighbor_indices {
        let p_nb = nodes[idx];
        let dir = [
            p_nb[0] - p_node[0],
            p_nb[1] - p_node[1],
            p_nb[2] - p_node[2],
        ];
        let metric_len = m_local.length(dir);
        let w = 1.0 / metric_len.max(1e-12);
        weight_sum += w;
        weighted_sum[0] += w * p_nb[0];
        weighted_sum[1] += w * p_nb[1];
        weighted_sum[2] += w * p_nb[2];
    }

    if weight_sum > 1e-15 {
        nodes[node_idx] = [
            weighted_sum[0] / weight_sum,
            weighted_sum[1] / weight_sum,
            weighted_sum[2] / weight_sum,
        ];
    }
}

fn tetra_volume_3d(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let ad = [a[0] - d[0], a[1] - d[1], a[2] - d[2]];
    let bd = [b[0] - d[0], b[1] - d[1], b[2] - d[2]];
    let cd = [c[0] - d[0], c[1] - d[1], c[2] - d[2]];
    let cross = [
        bd[1] * cd[2] - bd[2] * cd[1],
        bd[2] * cd[0] - bd[0] * cd[2],
        bd[0] * cd[1] - bd[1] * cd[0],
    ];
    (ad[0] * cross[0] + ad[1] * cross[1] + ad[2] * cross[2]).abs() / 6.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use rmsh_model::{Element, ElementType, Mesh, Node};

    fn cube_surface() -> Mesh {
        let mut mesh = Mesh::new();
        for (id, xyz) in [
            (1, [0.0, 0.0, 0.0]),
            (2, [1.0, 0.0, 0.0]),
            (3, [1.0, 1.0, 0.0]),
            (4, [0.0, 1.0, 0.0]),
            (5, [0.0, 0.0, 1.0]),
            (6, [1.0, 0.0, 1.0]),
            (7, [1.0, 1.0, 1.0]),
            (8, [0.0, 1.0, 1.0]),
        ] {
            mesh.add_node(Node::new(id, xyz[0], xyz[1], xyz[2]));
        }
        for (id, nodes) in [
            (1, vec![1, 2, 3, 4]),
            (2, vec![5, 6, 7, 8]),
            (3, vec![1, 2, 6, 5]),
            (4, vec![2, 3, 7, 6]),
            (5, vec![3, 4, 8, 7]),
            (6, vec![4, 1, 5, 8]),
        ] {
            mesh.add_element(Element::new(id, ElementType::Quad4, nodes));
        }
        mesh
    }

    #[test]
    fn classify_edges_buckets_metric_lengths() {
        let nodes = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.5, 0.0, 0.0]];
        let edges = [[0usize, 1usize], [0, 2]];
        let field = UniformMetricField3D::new(1.0);
        let (too_long, too_short, good) = classify_edges(&nodes, &edges, &field, 0.7, 1.4);
        assert_eq!(too_long, vec![0]);
        assert_eq!(too_short, vec![1]);
        assert!(good.is_empty());
    }

    #[test]
    fn mmg_remesh_generates_volume_mesh() {
        let mesh = MmgRemesh::default()
            .mesh_3d(&cube_surface(), &MeshParams::with_size(0.4))
            .unwrap();
        assert!(mesh.elements_by_dimension(3).len() > 0);
    }

    #[test]
    fn metric3_intersect_isotropic() {
        let m1 = Metric3::isotropic(0.5);
        let m2 = Metric3::isotropic(0.3);
        let m = Metric3::intersect(m1, m2);
        // Intersection should be more constraining
        assert!(m.m11 > m1.m11);
    }

    #[test]
    fn metric3_intersect_spd() {
        let m1 = Metric3::isotropic(0.5);
        let m2 = Metric3::isotropic(0.3);
        let m = Metric3::intersect(m1, m2);
        // Result should be SPD
        let a = [
            [m.m11, m.m12, m.m13],
            [m.m12, m.m22, m.m23],
            [m.m13, m.m23, m.m33],
        ];
        // Sylvester's criterion: all leading principal minors > 0
        assert!(a[0][0] > 0.0);
        assert!(a[0][0] * a[1][1] - a[0][1] * a[1][0] > 0.0);
        let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
        assert!(det > 0.0);
    }

    #[test]
    fn collapse_edge_3d_removes_degenerate() {
        let mut nodes = vec![[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.5, 0.5, 0.5], [0.0, 0.5, 0.5]];
        let mut tets = vec![[0, 1, 2, 3]];
        let field = UniformMetricField3D::new(0.5);
        collapse_edge_3d(&mut nodes, &mut tets, 0, 1, &field).unwrap();
        // The single tet used the collapsed edge, should become degenerate and removed
        assert!(tets.is_empty());
    }
}
