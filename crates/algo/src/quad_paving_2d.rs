//! Quad Paving 2-D — direct quadrilateral mesh generation (Gmsh algorithm 9).
//!
//! # Algorithm overview
//!
//! Gmsh's algorithm 9 ("Packing of Parallelograms") and algorithm 11
//! ("Quasi-Structured Quads") both target **all-quad** or **mostly-quad**
//! surface meshes.  This module implements the packing-of-parallelograms
//! approach as the primary strategy, with quad-dominant output.
//!
//! The family of algorithms proceeds as follows:
//!
//! 1. **Cross-field generation**: compute a smooth 4-direction field (a "cross
//!    field") over the domain that aligns with boundary curves and smoothly
//!    interpolates inward.  Each cross specifies the preferred quad orientation
//!    at that point.
//!
//! 2. **Streamline tracing**: trace streamlines (integral curves) in the
//!    cross-field directions to produce two families of curves that form the
//!    "skeleton" of the quad mesh.
//!
//! 3. **Quad patch construction**: identify closed loops in the streamline
//!    network that form quadrilateral patch boundaries, then fill each patch
//!    with a structured quad grid.
//!
//! 4. **Clean-up**: fix remaining irregular nodes, insert triangles at
//!    singularities of the cross field (unavoidable for non-topological disks).
//!
//! An alternative, simpler approach is **Q-Morph** (advancing-front quads):
//! start from a triangle mesh and convert triangles to quads via local
//! recombination operations guided by the same cross-field.  This is enabled
//! by the [`QuadStrategy::Recombine`] variant.
//!
//! # Reference
//!
//! Remacle et al., "Blossom-Quad…", *Int. J. Numer. Meth. Engng.* 89, 2012.
//! Viertel & Osting, "An Approach to Quad Meshing Based on Harmonic Cross-Valued
//! Maps", *SIAM J. Sci. Comput.* 41, 2019.
//! Gmsh source: `Mesh/meshGFaceQuadqs.cpp`.
//!
//! # Status
//!
//! The cross-field machinery is implemented: boundary-aligned field
//! computation (double-angle Laplacian smoothing, [`CrossField::compute`]),
//! barycentric evaluation ([`CrossField::evaluate_at`]), streamline tracing
//! ([`trace_streamline`]), and single-pair recombination
//! ([`recombine_triangle_pair`]). What is not wired yet is the full
//! packing-of-parallelograms pipeline (cross-field streamline network ->
//! quad patch identification -> structured patch filling -> clean-up);
//! [`QuadPaving2D::mesh_2d`] currently produces structured quads on
//! axis-aligned rectangles and falls back to a triangle mesh elsewhere
//! (or errors under `require_pure_quad`).

use std::collections::{BTreeSet, HashMap};

use rmsh_model::Mesh;

use crate::planar_meshing::{
    is_axis_aligned_rectangle, mesh_domain_triangles, structured_quad_mesh_rectangle,
};
use crate::traits::{Domain2D, MeshAlgoError, MeshParams, Mesher2D};

// ─── Strategy ────────────────────────────────────────────────────────────────

/// Strategy for producing a quad mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QuadStrategy {
    /// Packing of parallelograms: cross-field + streamline + patch (Gmsh 9).
    #[default]
    PackingOfParallelograms,

    /// Quasi-structured quads with better alignment for smooth geometry (Gmsh 11).
    QuasiStructured,

    /// Start from a triangle mesh and recombine triangles into quads (Blossom-Quad).
    Recombine,
}

// ─── Public struct ────────────────────────────────────────────────────────────

/// Quad-paving 2-D mesher (Gmsh algorithms 9 / 11).
///
/// Generates predominantly quadrilateral surface meshes by tracing a smooth
/// cross field and packing parallelogram-shaped quads along it.
#[derive(Debug, Clone)]
pub struct QuadPaving2D {
    /// Which quad-generation strategy to use.
    pub strategy: QuadStrategy,

    /// Number of cross-field smoothing iterations.
    ///
    /// Higher values yield a smoother cross field and typically better quads,
    /// at the cost of longer preprocessing.  Defaults to `100`.
    pub cross_field_iterations: u32,

    /// When `true`, any remaining triangles in the final mesh are reported as
    /// an error rather than left in place.
    ///
    /// In practice a small number of triangles at singular points is expected.
    /// Defaults to `false`.
    pub require_pure_quad: bool,
}

impl Default for QuadPaving2D {
    fn default() -> Self {
        Self {
            strategy: QuadStrategy::PackingOfParallelograms,
            cross_field_iterations: 100,
            require_pure_quad: false,
        }
    }
}

impl QuadPaving2D {
    pub fn new() -> Self {
        Self::default()
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl Mesher2D for QuadPaving2D {
    fn name(&self) -> &'static str {
        "Quad Paving 2D"
    }

    fn mesh_2d(&self, domain: &Domain2D, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        if domain.boundaries.len() == 1 {
            if let Some((min, max)) = is_axis_aligned_rectangle(domain.outer()) {
                return Ok(structured_quad_mesh_rectangle(
                    min,
                    max,
                    params.element_size,
                ));
            }
        }

        let tri_mesh =
            mesh_domain_triangles(domain, params.element_size, params.element_size, 0.0)?;
        if self.require_pure_quad {
            return Err(MeshAlgoError::Generation(
                "pure quad output currently requires an axis-aligned rectangular domain"
                    .to_string(),
            ));
        }
        Ok(tri_mesh)
    }
}

// ─── Internal helpers (stubs) ─────────────────────────────────────────────────

/// A smooth 4-direction (cross) field over the mesh.
///
/// At each triangle the field stores an angle `θ ∈ [0°, 90°)` representing
/// the primary quad direction.  The four actual directions are
/// `θ`, `θ + 90°`, `θ + 180°`, `θ + 270°`.
#[allow(dead_code)]
struct CrossField {
    /// Per-node angles (radians) representing the cross orientation.
    angles: Vec<f64>,
}

#[allow(dead_code)]
impl CrossField {
    /// Compute a smooth cross field on the interior of the given triangle mesh,
    /// with boundary alignment as a Dirichlet condition.
    ///
    /// Seeding: every boundary node takes the double-angle circular mean of
    /// its incident boundary-edge tangent directions (the mean of directions
    /// `t_i` lives at `atan2(Σ sin 2t_i, Σ cos 2t_i) / 2`), normalized into
    /// the canonical quadrant `[0, π/2)` — a cross is 4-fold symmetric, so
    /// opposite edge orientations are equivalent. Interior nodes start at 0.
    ///
    /// Smoothing: `iterations` Jacobi sweeps of Laplacian smoothing over the
    /// interior nodes in double-angle space `(cos 2θ, sin 2θ)` — the correct
    /// linear space for π/2-periodic angle data, where a plain average of raw
    /// angles would be wrong across the 0/90-degree wrap. Boundary nodes are
    /// never updated (Dirichlet).
    fn compute(
        nodes: &[[f64; 2]],
        triangles: &[[usize; 3]],
        boundary_edges: &[[usize; 2]],
        iterations: u32,
    ) -> Self {
        let n = nodes.len();

        // Node adjacency from the triangle edges (sorted sets: deterministic
        // sweeps).
        let mut neighbors: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
        for tri in triangles {
            for k in 0..3 {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                neighbors[a].insert(b);
                neighbors[b].insert(a);
            }
        }

        // Boundary node -> incident boundary-edge tangent directions.
        let mut boundary: HashMap<usize, Vec<[f64; 2]>> = HashMap::new();
        for edge in boundary_edges {
            let (a, b) = (edge[0], edge[1]);
            let t = [nodes[b][0] - nodes[a][0], nodes[b][1] - nodes[a][1]];
            boundary.entry(a).or_default().push(t);
            boundary.entry(b).or_default().push(t);
        }

        let mut angles = vec![0.0f64; n];
        for (&node, tangents) in &boundary {
            let (mut c, mut s) = (0.0f64, 0.0f64);
            for t in tangents {
                let len = (t[0] * t[0] + t[1] * t[1]).sqrt();
                if len <= f64::EPSILON {
                    continue;
                }
                let theta = (t[1] / len).atan2(t[0] / len);
                c += (2.0 * theta).cos();
                s += (2.0 * theta).sin();
            }
            if c != 0.0 || s != 0.0 {
                angles[node] = normalize_cross_angle(s.atan2(c) / 2.0);
            }
        }

        // Jacobi sweeps over interior nodes in double-angle space.
        for _ in 0..iterations {
            let mut next = angles.clone();
            for i in 0..n {
                if boundary.contains_key(&i) {
                    continue; // Dirichlet boundary condition
                }
                let (mut c, mut s) = (0.0f64, 0.0f64);
                for &j in &neighbors[i] {
                    c += (2.0 * angles[j]).cos();
                    s += (2.0 * angles[j]).sin();
                }
                if c != 0.0 || s != 0.0 {
                    next[i] = normalize_cross_angle(s.atan2(c) / 2.0);
                }
            }
            angles = next;
        }

        CrossField { angles }
    }

    /// Evaluate the cross direction at an arbitrary point by interpolating
    /// within the triangle containing it.
    ///
    /// The angle is interpolated as a circular mean in double-angle space
    /// `(cos 2θ, sin 2θ)` with standard barycentric weights, then returned as
    /// a unit direction vector. For a degenerate (zero-area) triangle the
    /// node-0 direction is returned.
    ///
    /// (Signature note: the skeleton passed no connectivity, which made the
    /// documented contract impossible — the triangle list is required to find
    /// the interpolated vertices. It is part of the signature now.)
    fn evaluate_at(
        &self,
        point: [f64; 2],
        tri_idx: usize,
        nodes: &[[f64; 2]],
        triangles: &[[usize; 3]],
    ) -> [f64; 2] {
        let tri = &triangles[tri_idx];
        let a = nodes[tri[0]];
        let b = nodes[tri[1]];
        let c = nodes[tri[2]];

        // Standard barycentric coordinates via the dot-product solve.
        let v0 = [b[0] - a[0], b[1] - a[1]];
        let v1 = [c[0] - a[0], c[1] - a[1]];
        let v2 = [point[0] - a[0], point[1] - a[1]];
        let d00 = v0[0] * v0[0] + v0[1] * v0[1];
        let d01 = v0[0] * v1[0] + v0[1] * v1[1];
        let d11 = v1[0] * v1[0] + v1[1] * v1[1];
        let d20 = v2[0] * v0[0] + v2[1] * v0[1];
        let d21 = v2[0] * v1[0] + v2[1] * v1[1];
        let denom = d00 * d11 - d01 * d01;
        if denom.abs() <= 1e-30 {
            // Degenerate triangle: fall back to node 0's direction.
            let theta = self.angles.get(tri[0]).copied().unwrap_or(0.0);
            return [theta.cos(), theta.sin()];
        }
        let w_b = (d11 * d20 - d01 * d21) / denom;
        let w_c = (d00 * d21 - d01 * d20) / denom;
        let w_a = 1.0 - w_b - w_c;

        // Circular mean in double-angle space, then back to a unit vector.
        let mut cs = 0.0f64;
        let mut sn = 0.0f64;
        for (w, node) in [(w_a, tri[0]), (w_b, tri[1]), (w_c, tri[2])] {
            let theta = 2.0 * self.angles[node];
            cs += w * theta.cos();
            sn += w * theta.sin();
        }
        let phi = sn.atan2(cs) / 2.0;
        [phi.cos(), phi.sin()]
    }
}

/// Trace a streamline in the cross field starting from `seed` in direction
/// `dir_index` (0 or 1) until it exits the domain or loops back.
///
/// Returns a polyline as a list of 2-D points (the seed is the first vertex).
///
/// * Integration: midpoint method (2-stage Runge-Kutta — plain Euler spirals
///   out of circular streamlines); `step_size` is the step length.
/// * Branch continuity: the cross field is an *axis* (unoriented), so the raw
///   interpolated direction flips sign across the ±180-degree wrap of the
///   double-angle representation. Each step's direction is therefore flipped
///   against the previous heading to keep the streamline continuous.
/// * Triangle tracking: the containing triangle is carried between steps and
///   refound through the edge-adjacency neighborhood, with a brute-force scan
///   as numerical fallback.
/// * Termination: the traced point leaves every triangle (domain exit — the
///   exit point is pushed as the last vertex), returns within one
///   `step_size` of the seed after at least four steps (loop closure), or a
///   10_000-step safety cap fires.
/// * A seed outside every triangle yields the single-point polyline
///   `[seed]`.
#[allow(dead_code)]
fn trace_streamline(
    seed: [f64; 2],
    dir_index: usize,
    cross_field: &CrossField,
    nodes: &[[f64; 2]],
    triangles: &[[usize; 3]],
    step_size: f64,
) -> Vec<[f64; 2]> {
    const MAX_STEPS: usize = 10_000;

    // Edge -> incident triangles, for neighbor-tracking refinds.
    let mut edge_tris: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (ti, tri) in triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            edge_tris.entry((a.min(b), a.max(b))).or_default().push(ti);
        }
    }

    // Evaluate the cross axis at `p`, take the `dir_index` branch (0 = primary
    // direction, 1 = perpendicular), and orient it to continue `reference`.
    let axis = |p: [f64; 2], tri: usize, reference: [f64; 2]| -> [f64; 2] {
        let mut d = cross_field.evaluate_at(p, tri, nodes, triangles);
        if dir_index == 1 {
            d = [-d[1], d[0]];
        }
        if d[0] * reference[0] + d[1] * reference[1] < 0.0 {
            d = [-d[0], -d[1]];
        }
        d
    };

    let mut polyline = vec![seed];
    let mut current = match locate_triangle(seed, nodes, triangles) {
        Some(t) => t,
        None => return polyline, // seed outside the domain
    };
    let mut p = seed;
    let mut heading = axis(p, current, [1.0, 0.0]);
    for step in 0..MAX_STEPS {
        // Midpoint step (2-stage Runge-Kutta), branch-continuous.
        let mid = [
            p[0] + 0.5 * step_size * heading[0],
            p[1] + 0.5 * step_size * heading[1],
        ];
        let d2 = axis(mid, current, heading);
        p = [
            p[0] + step_size * d2[0],
            p[1] + step_size * d2[1],
        ];
        heading = d2;

        // Triangle tracking: neighbors first, brute-force fallback.
        current = match refind_triangle(p, current, nodes, triangles, &edge_tris) {
            Some(t) => t,
            None => {
                polyline.push(p); // exit point on (or just past) the boundary
                return polyline;
            }
        };
        polyline.push(p);

        // Loop closure: back within one step of the seed.
        if step >= 4 {
            let dx = p[0] - seed[0];
            let dy = p[1] - seed[1];
            if dx * dx + dy * dy <= step_size * step_size {
                return polyline;
            }
        }
    }
    polyline
}

/// Barycentric coordinates of `p` in triangle `(a, b, c)`, or `None` for a
/// degenerate (zero-area) triangle.
fn barycentric(p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Option<[f64; 3]> {
    let d = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    if d.abs() <= 1e-30 {
        return None;
    }
    let la = ((b[0] - p[0]) * (c[1] - p[1]) - (b[1] - p[1]) * (c[0] - p[0])) / d;
    let lb = ((c[0] - p[0]) * (a[1] - p[1]) - (c[1] - p[1]) * (a[0] - p[0])) / d;
    let lc = 1.0 - la - lb;
    Some([la, lb, lc])
}

/// Containment test with a small boundary tolerance.
fn point_in_triangle(p: [f64; 2], tri: &[usize; 3], nodes: &[[f64; 2]]) -> bool {
    const TOL: f64 = -1e-9;
    match barycentric(p, nodes[tri[0]], nodes[tri[1]], nodes[tri[2]]) {
        Some([la, lb, lc]) => la >= TOL && lb >= TOL && lc >= TOL,
        None => false,
    }
}

/// Brute-force containing-triangle search.
fn locate_triangle(p: [f64; 2], nodes: &[[f64; 2]], triangles: &[[usize; 3]]) -> Option<usize> {
    triangles.iter().position(|tri| point_in_triangle(p, tri, nodes))
}

/// Refind the triangle containing `p`: the tracked triangle first, then its
/// edge-neighbors, then a full scan.
fn refind_triangle(
    p: [f64; 2],
    current: usize,
    nodes: &[[f64; 2]],
    triangles: &[[usize; 3]],
    edge_tris: &HashMap<(usize, usize), Vec<usize>>,
) -> Option<usize> {
    if point_in_triangle(p, &triangles[current], nodes) {
        return Some(current);
    }
    let tri = &triangles[current];
    for k in 0..3 {
        let (a, b) = (tri[k], tri[(k + 1) % 3]);
        if let Some(candidates) = edge_tris.get(&(a.min(b), a.max(b))) {
            for &t in candidates {
                if t != current && point_in_triangle(p, &triangles[t], nodes) {
                    return Some(t);
                }
            }
        }
    }
    locate_triangle(p, nodes, triangles)
}

/// Map an angle into the cross field's canonical half-open quadrant `[0, π/2)`.
fn normalize_cross_angle(theta: f64) -> f64 {
    let half = std::f64::consts::FRAC_PI_2;
    let mut t = theta % half;
    if t < 0.0 {
        t += half;
    }
    if t == 0.0 {
        return 0.0; // canonicalize -0.0
    }
    t
}

/// Convert a pair of adjacent triangles into a single quadrilateral by merging
/// them along their shared edge (the Blossom-Quad recombination step).
///
/// The quad lists its nodes in perimeter order: for triangles `ABC` and `DCB`
/// sharing edge `BC` the result is `A-B-D-C` (apex of the first triangle,
/// first shared vertex in the first triangle's vertex order, apex of the
/// second triangle, second shared vertex).
///
/// Returns `None` when the triangles do not share exactly one edge, or when
/// the union is not a strictly convex quadrilateral (same-side apexes,
/// collinear corners, zero area).
#[allow(dead_code)]
fn recombine_triangle_pair(
    tri_a: [usize; 3],
    tri_b: [usize; 3],
    nodes: &[[f64; 2]],
) -> Option<[usize; 4]> {
    let shared: Vec<_> = tri_a.into_iter().filter(|a| tri_b.contains(a)).collect();
    if shared.len() != 2 {
        return None;
    }
    let a_only = tri_a.into_iter().find(|n| !shared.contains(n))?;
    let b_only = tri_b.into_iter().find(|n| !shared.contains(n))?;
    let quad = [a_only, shared[0], b_only, shared[1]];

    // Contract gate: the union must be a strictly convex quadrilateral.
    let p = [
        nodes[quad[0]],
        nodes[quad[1]],
        nodes[quad[2]],
        nodes[quad[3]],
    ];
    let mut orientation = 0.0f64;
    for i in 0..4 {
        let (a, b, c) = (p[i], p[(i + 1) % 4], p[(i + 2) % 4]);
        let cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
        if cross.abs() <= 1e-12 {
            return None; // collinear corner (zero-area union)
        }
        if orientation == 0.0 {
            orientation = cross.signum();
        } else if cross.signum() != orientation {
            return None; // reflex corner or self-intersecting union
        }
    }
    Some(quad)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quad_paving_rectangle_produces_quads() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]);
        let mesh = QuadPaving2D::default()
            .mesh_2d(&domain, &MeshParams::with_size(0.5))
            .unwrap();
        assert!(
            mesh.elements
                .iter()
                .all(|e| e.etype == rmsh_model::ElementType::Quad4)
        );
    }

    // ── P3b: cross-field acceptance (CrossField::compute / evaluate_at) ──────

    /// Equilateral triangle with its three boundary edges, for pinning the
    /// boundary-alignment seeds exactly.
    fn equilateral_triangle() -> (
        Vec<[f64; 2]>,
        Vec<[usize; 3]>,
        Vec<[usize; 2]>,
    ) {
        (
            vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [0.5, 0.8660254037844386],
            ],
            vec![[0, 1, 2]],
            vec![[0, 1], [1, 2], [2, 0]],
        )
    }

    /// 3x3 node grid over the unit square (8 triangles, consistent split),
    /// with the perimeter as boundary edges.
    fn square_grid_3x3() -> (Vec<[f64; 2]>, Vec<[usize; 3]>, Vec<[usize; 2]>) {
        let mut nodes = Vec::new();
        for iy in 0..3 {
            for ix in 0..3 {
                nodes.push([0.5 * ix as f64, 0.5 * iy as f64]);
            }
        }
        let idx = |ix: usize, iy: usize| iy * 3 + ix;
        let mut tris = Vec::new();
        for iy in 0..2 {
            for ix in 0..2 {
                let (a, b, c, d) = (idx(ix, iy), idx(ix + 1, iy), idx(ix + 1, iy + 1), idx(ix, iy + 1));
                tris.push([a, b, c]);
                tris.push([a, c, d]);
            }
        }
        let mut boundary = Vec::new();
        for i in 0..2 {
            boundary.push([idx(i, 0), idx(i + 1, 0)]);
            boundary.push([idx(2, i), idx(2, i + 1)]);
            boundary.push([idx(2 - i, 2), idx(1 - i, 2)]);
            boundary.push([idx(0, 2 - i), idx(0, 1 - i)]);
        }
        (nodes, tris, boundary)
    }

    #[test]
    fn cross_field_compute_seeds_boundary_alignment_and_keeps_it_dirichlet() {
        let (nodes, tris, boundary) = equilateral_triangle();
        // With zero smoothing iterations the field is exactly the boundary
        // seed. Pinned by the double-angle mean of each node's two incident
        // boundary tangents: node 0 mixes 0 and 240 degrees -> 30 degrees;
        // node 1 mixes 0 and 120 -> 60 degrees; node 2 mixes 120 and 240
        // -> 180/2 -> normalizes to 0 in [0, 90).
        let f0 = CrossField::compute(&nodes, &tris, &boundary, 0);
        assert_eq!(f0.angles.len(), 3);
        let expect = [
            std::f64::consts::FRAC_PI_6,
            std::f64::consts::FRAC_PI_3,
            0.0,
        ];
        for i in 0..3 {
            assert!(
                (f0.angles[i] - expect[i]).abs() < 1e-9,
                "seed angle {i}: {} != {}",
                f0.angles[i],
                expect[i]
            );
        }
        // All nodes are boundary here: smoothing must not move them
        // (Dirichlet condition).
        let f10 = CrossField::compute(&nodes, &tris, &boundary, 10);
        for i in 0..3 {
            assert!((f10.angles[i] - expect[i]).abs() < 1e-9);
        }
    }

    #[test]
    fn cross_field_compute_stays_in_the_half_open_quadrant() {
        let (nodes, tris, boundary) = square_grid_3x3();
        let f = CrossField::compute(&nodes, &tris, &boundary, 50);
        assert_eq!(f.angles.len(), nodes.len());
        for (i, a) in f.angles.iter().enumerate() {
            assert!(
                *a >= 0.0 && *a < std::f64::consts::FRAC_PI_2,
                "angle at node {i} out of [0, 90): {a}"
            );
        }
    }

    #[test]
    fn cross_field_evaluate_at_interpolates_in_double_angle_space() {
        let nodes = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let tris = vec![[0usize, 1, 2]];
        let f = CrossField {
            angles: vec![0.0, std::f64::consts::PI / 12.0, std::f64::consts::PI / 6.0],
        };
        // Exactly at a vertex: that node's direction, unit norm.
        let d0 = f.evaluate_at([0.0, 0.0], 0, &nodes, &tris);
        assert!((d0[0] - 1.0).abs() < 1e-12 && d0[1].abs() < 1e-12);
        // At the centroid: circular (double-angle) mean of {0, 15, 30}deg
        // is 15 degrees.
        let c = f.evaluate_at([1.0 / 3.0, 1.0 / 3.0], 0, &nodes, &tris);
        let want = std::f64::consts::PI / 12.0;
        assert!((c[0] - want.cos()).abs() < 1e-9, "cx={} vs {}", c[0], want.cos());
        assert!((c[1] - want.sin()).abs() < 1e-9, "cy={} vs {}", c[1], want.sin());
        assert!((c[0] * c[0] + c[1] * c[1] - 1.0).abs() < 1e-12, "unit norm");
        // Degenerate triangle: falls back to node 0's direction, no panic.
        let degenerate_nodes = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]];
        let d2 = f.evaluate_at([0.5, 0.0], 0, &degenerate_nodes, &tris);
        assert!((d2[0] - 1.0).abs() < 1e-12 && d2[1].abs() < 1e-12);
    }

    // ── P3b: streamline acceptance ───────────────────────────────────────────

    /// A 2x1 rectangle split into two triangles.
    fn two_tri_rectangle() -> (Vec<[f64; 2]>, Vec<[usize; 3]>) {
        (
            vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]],
            vec![[0, 1, 2], [0, 2, 3]],
        )
    }

    #[test]
    fn trace_streamline_follows_the_cross_and_exits_the_domain() {
        let (nodes, tris) = two_tri_rectangle();
        let f = CrossField { angles: vec![0.0; 4] };

        // dir 0 of a zero cross is +x: the streamline runs to the right
        // boundary at constant y.
        let line = trace_streamline([0.2, 0.5], 0, &f, &nodes, &tris, 0.1);
        assert_eq!(line[0], [0.2, 0.5]);
        let last = line.last().unwrap();
        assert!(last[0] > 1.9, "must reach the right boundary, last={last:?}");
        for p in &line {
            assert!((p[1] - 0.5).abs() < 1e-12, "constant y, got {p:?}");
        }
        assert!(line.len() >= 10, "sensible step count, len={}", line.len());

        // dir 1 is the perpendicular branch: straight up at constant x.
        let up = trace_streamline([0.2, 0.5], 1, &f, &nodes, &tris, 0.1);
        let lastu = up.last().unwrap();
        assert!(lastu[1] > 0.95, "must reach the top boundary, last={lastu:?}");
        for p in &up {
            assert!((p[0] - 0.2).abs() < 1e-12, "constant x, got {p:?}");
        }
    }

    #[test]
    fn trace_streamline_closes_loops_in_a_swirling_field() {
        // Annulus with 24 spokes: inner ring r=0.7, outer ring r=0.8, every
        // node's angle = its polar angle + 90 degrees (exact CCW tangent).
        // Both rings agree at each spoke, so the field is exactly tangential
        // on the spokes and the residual wedge drift cancels by symmetry:
        // the orbit circulates at mid radius and terminates by loop closure
        // near the seed instead of exiting the domain.
        const SPOKES: usize = 24;
        let (r_in, r_out) = (0.7f64, 0.8f64);
        let mut nodes: Vec<[f64; 2]> = Vec::with_capacity(2 * SPOKES);
        let mut angles: Vec<f64> = Vec::with_capacity(2 * SPOKES);
        for k in 0..SPOKES {
            let phi = (k as f64) * std::f64::consts::TAU / SPOKES as f64;
            nodes.push([r_in * phi.cos(), r_in * phi.sin()]);
            angles.push(phi + std::f64::consts::FRAC_PI_2);
        }
        for k in 0..SPOKES {
            let phi = (k as f64) * std::f64::consts::TAU / SPOKES as f64;
            nodes.push([r_out * phi.cos(), r_out * phi.sin()]);
            angles.push(phi + std::f64::consts::FRAC_PI_2);
        }
        let mut tris: Vec<[usize; 3]> = Vec::with_capacity(2 * SPOKES);
        for k in 0..SPOKES {
            let k1 = (k + 1) % SPOKES;
            tris.push([k, SPOKES + k, SPOKES + k1]);
            tris.push([k, SPOKES + k1, k1]);
        }
        let f = CrossField { angles };

        let step = 0.05;
        let line = trace_streamline([0.75, 0.0], 0, &f, &nodes, &tris, step);
        let first = line[0];
        let last = line.last().unwrap();
        let disp = ((last[0] - first[0]).powi(2) + (last[1] - first[1]).powi(2)).sqrt();
        assert!(
            disp <= step + 1e-9,
            "loop must close within one step of the seed, disp={disp}, last={last:?}"
        );
        assert!(
            line.len() >= 20 && line.len() <= 1000,
            "a full circulation takes a sane number of steps, len={}",
            line.len()
        );
        for p in &line {
            let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!(
                r >= r_in - 1e-9 && r <= r_out + 1e-9,
                "streamline stayed inside the annulus, got {p:?} (r={r})"
            );
        }
    }

    #[test]
    fn trace_streamline_from_a_seed_outside_the_domain_returns_the_seed() {
        let (nodes, tris) = two_tri_rectangle();
        let f = CrossField { angles: vec![0.0; 4] };
        let line = trace_streamline([5.0, 0.5], 0, &f, &nodes, &tris, 0.1);
        assert_eq!(line.len(), 1);
        assert_eq!(line[0], [5.0, 0.5]);
    }

    // ── P3b: recombine_triangle_pair contract (convexity gate) ───────────────

    #[test]
    fn recombine_triangle_pair_orders_the_perimeter_and_gates_convexity() {
        // ABC=(0,1,2), DCB=(3,2,1) sharing BC: strictly convex union
        // (apexes on opposite sides of the shared edge).
        let nodes = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [2.0, 0.5]];
        let quad = recombine_triangle_pair([0, 1, 2], [3, 2, 1], &nodes)
            .expect("convex pair merges");
        assert_eq!(quad, [0, 1, 3, 2], "perimeter walk A-B-D-C");

        // Same-side apexes: the union is non-convex -> None.
        let same_side = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.5, 0.5]];
        assert!(recombine_triangle_pair([0, 1, 2], [3, 2, 1], &same_side).is_none());

        // Collinear corners: zero-area union -> None.
        let collinear = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]];
        assert!(recombine_triangle_pair([0, 1, 2], [1, 3, 2], &collinear).is_none());

        // Non-adjacent triangles -> None.
        let far = vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [10.0, 10.0],
            [11.0, 10.0],
            [12.0, 12.0],
        ];
        assert!(recombine_triangle_pair([0, 1, 2], [3, 4, 5], &far).is_none());
    }
}
