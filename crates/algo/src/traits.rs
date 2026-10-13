//! Abstract traits for finite element mesh generation and optimization algorithms.
//!
//! Mirrors the Gmsh algorithm taxonomy:
//!
//! | Trait | Purpose |
//! |---|---|
//! | [`Mesher2D`] | 2-D surface mesh generation (triangles, quads) |
//! | [`Mesher3D`] | 3-D volume mesh generation (tetrahedra, hexahedra) |
//! | [`MeshOptimizer`] | Post-generation mesh quality improvement |
//!
//! Each concrete algorithm lives in its own module and implements one of these
//! traits, making them independently replaceable and testable.
//!
//! # Usage pattern
//!
//! ```
//! use rmsh_algo::traits::{Domain2D, MeshParams, Mesher2D};
//! use rmsh_algo::FrontalDelaunay2D;
//!
//! let mesher = FrontalDelaunay2D::default();
//! let domain = Domain2D::from_outer(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
//! let params = MeshParams::with_size(0.2);
//! let mesh = mesher.mesh_2d(&domain, &params).expect("planar domain meshes");
//! assert!(!mesh.elements.is_empty());
//! ```

use rmsh_model::Mesh;
use thiserror::Error;

// ─── Error Type ───────────────────────────────────────────────────────────────

/// Errors produced by mesh generation and optimization algorithms.
#[derive(Error, Debug, Clone)]
pub enum MeshAlgoError {
    /// The algorithm implementation is a placeholder stub (not yet written).
    #[error("algorithm not yet implemented")]
    NotImplemented,

    /// The supplied input geometry or parameters are invalid.
    #[error("invalid input: {0}")]
    InvalidInput(String),

    /// An internal mesh generation failure occurred.
    #[error("mesh generation failed: {0}")]
    Generation(String),

    /// Mesh optimization could not converge or encountered an error.
    #[error("mesh optimization failed: {0}")]
    Optimization(String),
}

// ─── Common Parameter Types ────────────────────────────────────────────────────

/// Parameters shared by all mesh generation algorithms.
///
/// All size constraints are expressed in the same length unit as the input geometry.
#[derive(Debug, Clone)]
pub struct MeshParams {
    /// Target element edge length (characteristic length *h*).
    ///
    /// Controls the coarseness / density of the output mesh.
    pub element_size: f64,

    /// Minimum allowed element edge length.
    ///
    /// Elements smaller than this are not normally produced. Defaults to `element_size / 10`.
    pub min_size: f64,

    /// Maximum allowed element edge length.
    ///
    /// Edges longer than this are subdivided. Defaults to `element_size * 2`.
    pub max_size: f64,

    /// Number of mesh-quality optimization passes executed *after* the initial
    /// mesh is generated. Set to `0` to skip optimization entirely.
    pub optimize_passes: u32,
}

impl MeshParams {
    /// Construct a [`MeshParams`] with sensible defaults for the given target size.
    pub fn with_size(element_size: f64) -> Self {
        Self {
            element_size,
            min_size: element_size / 10.0,
            max_size: element_size * 2.0,
            optimize_passes: 3,
        }
    }
}

/// Parameters for mesh-quality optimization algorithms.
#[derive(Debug, Clone)]
pub struct OptimizeParams {
    /// Maximum number of smoothing or optimization iterations.
    pub iterations: u32,

    /// Convergence tolerance: stop early when the average node displacement or
    /// quality improvement falls below this value.
    pub tolerance: f64,

    /// When `true`, boundary nodes may be relocated along the boundary curve.
    /// When `false` (the default), boundary nodes are kept fixed.
    pub move_boundary_nodes: bool,
}

impl Default for OptimizeParams {
    fn default() -> Self {
        Self {
            iterations: 10,
            tolerance: 1e-6,
            move_boundary_nodes: false,
        }
    }
}

// ─── Domain Types ──────────────────────────────────────────────────────────────

/// A 2-D meshing domain defined by one or more closed boundary polylines.
///
/// * The **first** polyline is the outer boundary — vertices listed
///   counter-clockwise.
/// * Any additional polylines represent **holes** — vertices listed clockwise.
///
/// All coordinates are in 2-D (XY plane, z = 0).
pub struct Domain2D {
    /// Outer boundary followed by zero or more hole boundaries.
    pub boundaries: Vec<Vec<[f64; 2]>>,
}

impl Domain2D {
    /// Create a domain with just an outer boundary (no holes).
    pub fn from_outer(outer: Vec<[f64; 2]>) -> Self {
        Self {
            boundaries: vec![outer],
        }
    }

    /// Return the outer boundary polyline.
    pub fn outer(&self) -> &[[f64; 2]] {
        &self.boundaries[0]
    }

    /// Return hole polylines (may be empty).
    pub fn holes(&self) -> &[Vec<[f64; 2]>] {
        if self.boundaries.len() > 1 {
            &self.boundaries[1..]
        } else {
            &[]
        }
    }

    /// Append a hole to the domain and return `self` (builder pattern).
    pub fn with_hole(mut self, hole: Vec<[f64; 2]>) -> Self {
        self.boundaries.push(hole);
        self
    }
}

// ─── Traits ───────────────────────────────────────────────────────────────────

/// A 2-D surface meshing algorithm.
///
/// Implementors produce a triangular or quadrilateral surface mesh that fills
/// the given planar [`Domain2D`].
///
/// # Gmsh counterparts
///
/// | Gmsh algo # | Module |
/// |---|---|
/// | 1 — MeshAdapt | [`crate::mesh_adapt_2d`] |
/// | 6 — Frontal-Delaunay | [`crate::frontal_delaunay_2d`] |
/// | 7 — BAMG | [`crate::bamg_2d`] |
/// | 9 — Packing of Parallelograms | [`crate::quad_paving_2d`] |
pub trait Mesher2D {
    /// Human-readable algorithm name, e.g. `"Frontal-Delaunay 2D"`.
    fn name(&self) -> &'static str;

    /// Generate a 2-D mesh for `domain` using the specified `params`.
    ///
    /// # Returns
    ///
    /// A [`Mesh`] whose nodes live in the XY plane (`z = 0`).
    fn mesh_2d(&self, domain: &Domain2D, params: &MeshParams) -> Result<Mesh, MeshAlgoError>;
}

/// A 3-D volume meshing algorithm.
///
/// Implementors produce a tetrahedral or hexahedral volume mesh that fills
/// the interior of the closed surface provided as input.
///
/// # Gmsh counterparts
///
/// | Gmsh algo # | Module |
/// |---|---|
/// | 1 — Delaunay 3D | [`crate::delaunay_3d`] |
/// | 4 — Frontal | [`crate::frontal_3d`] |
/// | 7 — MMG3D | [`crate::mmg_remesh`] |
/// | 10 — HXT | [`crate::hxt_3d`] |
pub trait Mesher3D {
    /// Human-readable algorithm name, e.g. `"HXT Parallel Delaunay"`.
    fn name(&self) -> &'static str;

    /// Generate a volume mesh inside the given closed `surface`.
    ///
    /// # Parameters
    ///
    /// The `surface` mesh must be a closed, manifold triangular shell.
    /// All vertices of the output mesh that coincide with `surface` vertices
    /// retain the same node IDs.
    ///
    /// # Returns
    ///
    /// A [`Mesh`] containing the generated volume elements (Tet4, Hex8, …).
    fn mesh_3d(&self, surface: &Mesh, params: &MeshParams) -> Result<Mesh, MeshAlgoError>;
}

/// A mesh quality optimizer or smoother.
///
/// Implementors improve mesh quality **in-place**, for example via:
/// - Node relocation (smoothing)
/// - Local topological operations (edge / face swapping)
///
/// # Gmsh counterparts
///
/// | Optimizer | Module |
/// |---|---|
/// | Laplacian Smooth | [`crate::laplacian_smooth`] |
/// | Quality Optimizer | [`crate::mesh_optimize`] |
pub trait MeshOptimizer {
    /// Human-readable optimizer name, e.g. `"Laplacian Smooth"`.
    fn name(&self) -> &'static str;

    /// Improve the quality of `mesh` in-place.
    fn optimize(&self, mesh: &mut Mesh, params: &OptimizeParams) -> Result<(), MeshAlgoError>;
}
