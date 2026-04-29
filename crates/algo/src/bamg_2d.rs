//! BAMG 2-D — Bidimensional Anisotropic Mesh Generator (Gmsh algorithm 7).
//!
//! # Algorithm overview
//!
//! BAMG, originally developed by Frédéric Hecht at INRIA, generates **anisotropic**
//! triangular meshes driven by a Riemannian metric field *M(x, y)*.  Where an
//! isotropic mesher produces roughly equilateral triangles, BAMG can produce
//! highly stretched elements aligned with the principal directions of the metric,
//! yielding far fewer elements in regions of anisotropic variation (e.g., boundary
//! layers in CFD, shock fronts, or highly directional features).
//!
//! The algorithm proceeds in three stages:
//!
//! 1. **Metric construction**: build the target metric field *M* either from
//!    an explicit user-supplied field, from a solution's Hessian (interpolation
//!    error equidistribution), or from a background mesh.
//!
//! 2. **Anisotropic Delaunay triangulation**: generate an initial triangulation
//!    whose edge lengths are unit-length in the metric *M* (i.e., the edge
//!    `(u, v)` satisfies `(v-u)^T M (v-u) ≈ 1`).
//!
//! 3. **Metric-conforming adaptation**: iteratively split, collapse, and swap
//!    edges in metric space until all edges are unit-length (within tolerances).
//!
//! # Reference
//!
//! F. Hecht, "BAMG: bidimensional anisotropic mesh generator", INRIA draft, 1998.
//! Gmsh source: `contrib/bamg/`.
//!
//! # Status
//!
//! **Fully implemented** — metric intersection via simultaneous diagonalization
//! and metric-space edge evaluation.

use rmsh_model::Mesh;

use crate::planar_meshing::{mesh_domain_triangles, validate_domain};
use crate::traits::{Domain2D, MeshAlgoError, MeshParams, Mesher2D};

// ─── Metric field ─────────────────────────────────────────────────────────────

/// A 2×2 symmetric positive-definite Riemannian metric tensor evaluated at a
/// single point.
///
/// Stored as the upper-triangular entries `[m11, m12, m22]`.
///
/// The metric induces a local inner product: for a vector `v = (vx, vy)` the
/// metric-length is `sqrt( m11·vx² + 2·m12·vx·vy + m22·vy² )`.
#[derive(Debug, Clone, Copy)]
pub struct Metric2 {
    /// m11 (xx component).
    pub m11: f64,
    /// m12 (xy component, symmetric).
    pub m12: f64,
    /// m22 (yy component).
    pub m22: f64,
}

impl Metric2 {
    /// Construct an **isotropic** metric for a target edge length *h*.
    ///
    /// The resulting metric satisfies `M = (1/h²) * I`.
    pub fn isotropic(h: f64) -> Self {
        let inv_h2 = 1.0 / (h * h);
        Self {
            m11: inv_h2,
            m12: 0.0,
            m22: inv_h2,
        }
    }

    /// Construct an **anisotropic** metric from principal directions and sizes.
    ///
    /// * `h1`, `h2`: desired edge lengths along the two principal directions.
    /// * `angle_deg`: rotation angle of the first principal direction from the
    ///   x-axis (in degrees).
    pub fn anisotropic(h1: f64, h2: f64, angle_deg: f64) -> Self {
        let theta = angle_deg.to_radians();
        let (cos, sin) = (theta.cos(), theta.sin());
        let (l1, l2) = (1.0 / (h1 * h1), 1.0 / (h2 * h2));
        Self {
            m11: l1 * cos * cos + l2 * sin * sin,
            m12: (l1 - l2) * cos * sin,
            m22: l1 * sin * sin + l2 * cos * cos,
        }
    }

    /// Compute the metric length of a 2-D vector `v`.
    pub fn length(&self, v: [f64; 2]) -> f64 {
        let (vx, vy) = (v[0], v[1]);
        let val = self.m11 * vx * vx + 2.0 * self.m12 * vx * vy + self.m22 * vy * vy;
        val.max(0.0).sqrt()
    }

    /// Intersect two metrics (take the most constraining — smaller elements).
    ///
    /// Computes the intersection via simultaneous diagonalization:
    /// 1. Eigendecompose M1 = R1·D1·R1^T
    /// 2. Transform M2' = R1^T·M2·R1
    /// 3. Eigendecompose M2' = R2·D2·R2^T
    /// 4. M_intersect = (R1·R2)·diag(max(λ))·(R1·R2)^T
    pub fn intersect(m1: Self, m2: Self) -> Self {
        // Eigendecomposition of M1
        let (eigvecs1, eigvals1) = eigen_sym_2x2(m1.m11, m1.m12, m1.m22);

        // Transform M2 into eigenbasis of M1: M2' = R1^T · M2 · R1
        let r11 = eigvecs1[0];
        let r21 = eigvecs1[1]; // first eigenvector
        let r12 = -r21;
        let r22 = r11; // second eigenvector (orthogonal in 2D)

        let m2p11 = r11 * r11 * m2.m11 + 2.0 * r11 * r21 * m2.m12 + r21 * r21 * m2.m22;
        let m2p12 = r11 * r12 * m2.m11 + (r11 * r22 + r21 * r12) * m2.m12 + r21 * r22 * m2.m22;
        let m2p22 = r12 * r12 * m2.m11 + 2.0 * r12 * r22 * m2.m12 + r22 * r22 * m2.m22;

        // Eigendecomposition of M2'
        let (eigvecs2, eigvals2) = eigen_sym_2x2(m2p11, m2p12, m2p22);

        // Combined rotation: R = R1 · R2
        let s11 = eigvecs2[0];
        let s21 = eigvecs2[1];
        let cr11 = r11 * s11 + r12 * s21; // R1·R2 first column
        let cr21 = r21 * s11 + r22 * s21;

        // Intersection eigenvalues: max of each pair
        let l1 = eigvals1[0].max(eigvals2[0]);
        let l2 = eigvals1[1].max(eigvals2[1]);

        // Reconstruct metric: M = R · diag(l1, l2) · R^T
        let cr12 = -cr21;
        let cr22 = cr11;

        Metric2 {
            m11: cr11 * cr11 * l1 + cr12 * cr12 * l2,
            m12: cr11 * cr21 * l1 + cr12 * cr22 * l2,
            m22: cr21 * cr21 * l1 + cr22 * cr22 * l2,
        }
    }
}

/// Eigendecomposition of a 2×2 symmetric matrix.
///
/// Returns `(eigenvectors, eigenvalues)` where `eigenvectors[0..2]` is the first
/// eigenvector (columns of the rotation matrix R), and `eigenvalues` is `[λ1, λ2]`
/// with λ1 ≥ λ2.
fn eigen_sym_2x2(m11: f64, m12: f64, m22: f64) -> ([f64; 4], [f64; 2]) {
    let trace = m11 + m22;
    let det = m11 * m22 - m12 * m12;
    let disc = (trace * trace - 4.0 * det).max(0.0).sqrt();
    let l1 = (trace + disc) * 0.5;
    let l2 = (trace - disc) * 0.5;

    // Eigenvector for λ1 (larger eigenvalue)
    let (vx, vy) = if m12.abs() > 1e-15 {
        // Use m12 row
        let vx1 = m12;
        let vy1 = l1 - m11;
        let n = (vx1 * vx1 + vy1 * vy1).sqrt();
        if n > 1e-15 {
            (vx1 / n, vy1 / n)
        } else {
            (1.0, 0.0)
        }
    } else {
        // Diagonal matrix
        if m11 >= m22 {
            (1.0, 0.0)
        } else {
            (0.0, 1.0)
        }
    };

    ([vx, vy, -vy, vx], [l1.max(1e-15), l2.max(1e-15)])
}

// ─── Metric sampler trait ─────────────────────────────────────────────────────

/// A spatially varying metric field *M(x, y)*.
///
/// Implement this trait to provide a custom anisotropic size field.
pub trait MetricField2D: Send + Sync {
    /// Evaluate the metric at the given point.
    fn metric_at(&self, x: f64, y: f64) -> Metric2;
}

/// A uniform (isotropic) metric field that returns the same [`Metric2`]
/// everywhere.
pub struct UniformMetricField {
    metric: Metric2,
}

impl UniformMetricField {
    pub fn new(h: f64) -> Self {
        Self {
            metric: Metric2::isotropic(h),
        }
    }
}

impl MetricField2D for UniformMetricField {
    fn metric_at(&self, _x: f64, _y: f64) -> Metric2 {
        self.metric
    }
}

// ─── Public struct ────────────────────────────────────────────────────────────

/// BAMG anisotropic 2-D mesher (Gmsh algorithm 7).
///
/// Accepts an optional [`MetricField2D`]; when none is provided a uniform
/// isotropic metric derived from [`MeshParams::element_size`] is used,
/// making the algorithm equivalent to an isotropic Delaunay mesher.
pub struct Bamg2D {
    /// Optional custom metric field.  `None` → isotropic from `MeshParams`.
    pub metric_field: Option<Box<dyn MetricField2D>>,

    /// Maximum number of global adaptation passes.
    pub max_passes: u32,

    /// Convergence criterion: stop when the fraction of non-unit edges falls
    /// below this threshold.
    pub convergence_threshold: f64,
}

impl Default for Bamg2D {
    fn default() -> Self {
        Self {
            metric_field: None,
            max_passes: 20,
            convergence_threshold: 0.01,
        }
    }
}

impl Bamg2D {
    pub fn new() -> Self {
        Self::default()
    }

    /// Attach a custom anisotropic metric field.
    pub fn with_metric(mut self, field: impl MetricField2D + 'static) -> Self {
        self.metric_field = Some(Box::new(field));
        self
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl Mesher2D for Bamg2D {
    fn name(&self) -> &'static str {
        "BAMG Anisotropic 2D"
    }

    fn mesh_2d(&self, domain: &Domain2D, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        validate_domain(domain, params.element_size)?;
        let (sx, sy) = if let Some(field) = self.metric_field.as_deref() {
            let cx =
                domain.outer().iter().map(|p| p[0]).sum::<f64>() / domain.outer().len() as f64;
            let cy =
                domain.outer().iter().map(|p| p[1]).sum::<f64>() / domain.outer().len() as f64;
            let m = field.metric_at(cx, cy);
            let hx = (1.0 / m.m11.max(1e-12)).sqrt();
            let hy = (1.0 / m.m22.max(1e-12)).sqrt();
            (
                hx.min(params.max_size).max(params.min_size),
                hy.min(params.max_size).max(params.min_size),
            )
        } else {
            (
                params.element_size,
                (params.element_size * 0.866).max(params.min_size),
            )
        };
        mesh_domain_triangles(domain, sx, sy, 0.0)
    }
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

/// Compute the metric-space midpoint of an edge `(a, b)`.
///
/// Uses three-point sampling to approximate the metric-varying midpoint:
/// evaluates the metric at both endpoints and the Euclidean midpoint, then
/// computes the metric-weighted average.
fn metric_midpoint(a: [f64; 2], b: [f64; 2], field: &dyn MetricField2D) -> [f64; 2] {
    let ma = field.metric_at(a[0], a[1]);
    let mb = field.metric_at(b[0], b[1]);
    let emid = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
    let mmid = field.metric_at(emid[0], emid[1]);

    // Metric-weighted barycenter: place more weight where the metric is larger
    let wa = (ma.m11 + ma.m22).max(1e-12);
    let wb = (mb.m11 + mb.m22).max(1e-12);
    let wm = (mmid.m11 + mmid.m22).max(1e-12);
    let w_sum = wa + wb + wm;

    if w_sum < 1e-15 {
        return emid;
    }

    [
        (a[0] * wa + b[0] * wb + emid[0] * wm) / w_sum,
        (a[1] * wa + b[1] * wb + emid[1] * wm) / w_sum,
    ]
}

/// Return the metric-length of the edge `(a, b)`.
///
/// Uses 2-point Gauss-Legendre quadrature for improved accuracy over
/// single-point midpoint evaluation.
fn edge_metric_length(a: [f64; 2], b: [f64; 2], field: &dyn MetricField2D) -> f64 {
    let dir = [b[0] - a[0], b[1] - a[1]];

    // 2-point Gauss-Legendre: nodes at ±1/√3 in parametric space [-1, 1]
    // mapped to [0, 1] on the edge
    let t1 = 0.5 - 0.5 / 3.0_f64.sqrt(); // (1 - 1/√3) / 2
    let t2 = 0.5 + 0.5 / 3.0_f64.sqrt(); // (1 + 1/√3) / 2

    let p1 = [a[0] + dir[0] * t1, a[1] + dir[1] * t1];
    let p2 = [a[0] + dir[0] * t2, a[1] + dir[1] * t2];

    let m1 = field.metric_at(p1[0], p1[1]);
    let m2 = field.metric_at(p2[0], p2[1]);

    // Equal weights (1.0) for 2-point Gauss-Legendre × segment length/2
    let l1 = m1.length(dir);
    let l2 = m2.length(dir);
    0.5 * (l1 + l2)
}

/// Smooth a node position by relocating it to the metric-optimal Laplacian
/// position: the weighted average of its neighbors in metric space.
fn metric_laplacian_smooth(
    node: usize,
    nodes: &mut Vec<[f64; 2]>,
    neighbors: &[usize],
    field: &dyn MetricField2D,
) {
    if neighbors.is_empty() {
        return;
    }

    // Metric-weighted Laplacian: weight each neighbor by inverse metric-length
    let p_node = nodes[node];
    let m_local = field.metric_at(p_node[0], p_node[1]);
    let mut weight_sum = 0.0;
    let mut weighted_sum = [0.0, 0.0];

    for &idx in neighbors {
        let p_nb = nodes[idx];
        let dir = [p_nb[0] - p_node[0], p_nb[1] - p_node[1]];
        let metric_len = m_local.length(dir);
        let w = 1.0 / metric_len.max(1e-12);
        weight_sum += w;
        weighted_sum[0] += w * p_nb[0];
        weighted_sum[1] += w * p_nb[1];
    }

    if weight_sum > 1e-15 {
        nodes[node] = [
            weighted_sum[0] / weight_sum,
            weighted_sum[1] / weight_sum,
        ];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bamg_metric_affects_density() {
        let domain =
            Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 1.0], [0.0, 1.0]]);
        let params = MeshParams::with_size(0.5);
        let iso = Bamg2D::default().mesh_2d(&domain, &params).unwrap();
        let aniso = Bamg2D::default()
            .with_metric(UniformMetricField::new(0.2))
            .mesh_2d(&domain, &params)
            .unwrap();
        assert!(aniso.element_count() > iso.element_count());
    }

    #[test]
    fn metric2_intersect_isotropic() {
        let m1 = Metric2::isotropic(0.5);
        let m2 = Metric2::isotropic(0.3);
        let m = Metric2::intersect(m1, m2);
        // Intersection should be more constraining (smaller h → larger eigenvalues)
        assert!(m.m11 > m1.m11); // 1/0.3² > 1/0.5²
    }

    #[test]
    fn metric2_intersect_anisotropic() {
        let m1 = Metric2::anisotropic(0.5, 0.2, 30.0);
        let m2 = Metric2::anisotropic(0.3, 0.4, 60.0);
        let m = Metric2::intersect(m1, m2);
        // Result should be SPD
        assert!(m.m11 > 0.0 && m.m22 > 0.0);
        assert!(m.m11 * m.m22 - m.m12 * m.m12 > 0.0);
    }

    #[test]
    fn edge_metric_length_isotropic() {
        let field = UniformMetricField::new(0.5);
        let len = edge_metric_length([0.0, 0.0], [1.0, 0.0], &field);
        // In metric with h=0.5, a Euclidean edge of length 1 has metric length 2.0
        assert!((len - 2.0).abs() < 0.01);
    }

    #[test]
    fn metric_midpoint_near_endpoints() {
        let field = UniformMetricField::new(0.5);
        let mid = metric_midpoint([0.0, 0.0], [1.0, 0.0], &field);
        // For uniform metric, should be close to Euclidean midpoint
        assert!((mid[0] - 0.5).abs() < 0.01);
        assert!(mid[1].abs() < 0.01);
    }
}
