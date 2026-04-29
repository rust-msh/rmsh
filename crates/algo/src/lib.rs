pub use rmsh_io::{MshError, parse_msh};

// ─── Existing algorithms ──────────────────────────────────────────────────────

pub mod triangulate2d;
pub use triangulate2d::{MeshError, Polygon2D, mesh_polygon, triangulate_points};

pub(crate) mod planar_meshing;

pub mod tetrahedralize3d;
pub use tetrahedralize3d::{CentroidStarMesher3D, Mesh3DError, tetrahedralize_closed_surface};

// ─── Abstract traits ──────────────────────────────────────────────────────────

/// Abstract traits for mesh generation and optimization algorithms.
///
/// All concrete algorithms implement one of [`traits::Mesher2D`],
/// [`traits::Mesher3D`], or [`traits::MeshOptimizer`].
pub mod traits;
pub use traits::{
    Domain2D, MeshAlgoError, MeshOptimizer, MeshParams, Mesher2D, Mesher3D, OptimizeParams,
};

// ─── 2-D surface meshing algorithms ──────────────────────────────────────────

/// MeshAdapt 2-D: local edge-split/collapse/swap refinement (Gmsh algo 1).
pub mod mesh_adapt_2d;
pub use mesh_adapt_2d::MeshAdapt2D;

/// Frontal-Delaunay 2-D: advancing front + Delaunay insertion (Gmsh algo 6).
pub mod frontal_delaunay_2d;
pub use frontal_delaunay_2d::FrontalDelaunay2D;

/// BAMG: bidimensional anisotropic mesh generator (Gmsh algo 7).
pub mod bamg_2d;
pub use bamg_2d::{Bamg2D, Metric2, MetricField2D, UniformMetricField};

/// Delaunay 2-D: Bowyer-Watson triangulation (Gmsh algo 5).
pub mod delaunay_2d;
pub use delaunay_2d::Delaunay2D;

/// Quad Paving 2-D: packing of parallelograms / cross-field quads (Gmsh algo 9/11).
pub mod quad_paving_2d;
pub use quad_paving_2d::{QuadPaving2D, QuadStrategy};

/// Frontal-Quads 2-D: frontal Delaunay + cross-field quad recombination (Gmsh algo 8).
pub mod frontal_quads_2d;
pub use frontal_quads_2d::FrontalQuads2D;

// ─── 3-D volume meshing algorithms ───────────────────────────────────────────

/// Delaunay 3-D: incremental Bowyer-Watson + Delaunay refinement (Gmsh algo 1).
pub mod delaunay_3d;
pub use delaunay_3d::Delaunay3D;

/// Compact tetrahedral mesh with neighbor tables (used internally by Delaunay3D).
pub mod tet_mesh;

/// Frontal-Delaunay 3-D: advancing-front tetrahedralization (Gmsh algo 4).
pub mod frontal_3d;
pub use frontal_3d::Frontal3D;

/// HXT: high-performance parallel Delaunay tetrahedralization (Gmsh algo 10).
pub mod hxt_3d;
pub use hxt_3d::Hxt3D;

/// MMG3D: anisotropic surface and volume remeshing (Gmsh algo 7).
pub mod mmg_remesh;
pub use mmg_remesh::{Metric3, MetricField3D, MmgRemesh, UniformMetricField3D};

// ─── Mesh optimization ────────────────────────────────────────────────────────

/// Laplacian smoothing: iterative node relocation toward neighbour centroid.
pub mod laplacian_smooth;
pub use laplacian_smooth::{LaplacianSmooth, LaplacianVariant};

/// Mesh quality optimizer: edge swaps, node insertion/collapse, smoothing.
pub mod mesh_optimize;
pub use mesh_optimize::{MeshQualityOptimizer, OptimizeConfig, QualityMetric};
