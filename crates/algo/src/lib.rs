pub use rmsh_io::{MshError, parse_msh};

// ─── Existing algorithms ──────────────────────────────────────────────────────

pub mod triangulate2d;
pub use triangulate2d::{MeshError, Polygon2D, mesh_polygon, triangulate_points};

/// Robust adaptive-exact geometric predicates (internal machinery of the 2-D
/// Delaunay codes).
pub(crate) mod predicates;

pub(crate) mod planar_meshing;

/// Constrained Delaunay triangulation with incremental insertion (internal
/// machinery of the 2-D meshers).
pub(crate) mod cdt_2d;

/// 2-D quality optimization: Gmsh-style relocation + swap passes (internal).
pub(crate) mod mesh_quality_2d;

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

// ─── Geometry-source abstraction (dependency inversion, P2) ───────────────────

/// The meshing-oriented geometry abstraction: closest-point / normal /
/// curvature queries, defined here (zero-rcad) and implemented by rcad
/// adapters on the consuming side (see the `rmsh-bridge` crate).
pub mod geometry_source;
pub use geometry_source::{GeometryHit, GeometrySource};

// ─── Size-field framework (adaptive sizing, P2) ───────────────────────────────

/// The size-field framework: uniform / analytic-closure / curvature-derived /
/// distance-derived element-size fields, min-max composition, and gradient
/// limiting over a background lattice (the sizing ground for adaptive
/// meshing, P3+).
pub mod size_field;
pub use size_field::{
    Composition, ComposedField, CurvatureField, DistanceField, FieldFn, GradientLimitedField,
    SizeField, SizeFieldError, UniformField,
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

/// Triangle-to-quadrilateral recombination (P3b): Blossom-style pairing of
/// adjacent triangle pairs scored by quality.
pub mod quad_recombine;

/// Quad Paving 2-D: packing of parallelograms / cross-field quads (Gmsh algo 9/11).
pub mod quad_paving_2d;
pub use quad_paving_2d::{QuadPaving2D, QuadStrategy};

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

// ─── 2-D boundary layers (P3c) ────────────────────────────────────────────────

/// 2-D boundary layers: structured quad/tri rows grown inward from a domain
/// boundary with local collision clamping (Gmsh boundary-layer architecture,
/// re-implemented; see the module docs for provenance and v1 scope).
pub mod boundary_layer_2d;
pub use boundary_layer_2d::{BoundaryLayer2D, BoundaryLayerError, cumulative_height, insert_boundary_layers};

/// 3-D boundary layers (P4): structured prism columns grown inward from a
/// closed triangle surface with per-vertex ray-cap collision clamping, the
/// remaining cavity filled by the centroid-star tetrahedralizer (Gmsh
/// boundary-layer architecture, re-implemented; see the module docs for
/// provenance and v1 scope).
pub mod boundary_layer_3d;
pub use boundary_layer_3d::{BoundaryLayer3D, BoundaryLayer3DError, insert_boundary_layers_3d};

/// Multi-material / partitioned volumes (P4): material classification and
/// cavity subtraction over closed surface meshes (v1 — see the module docs
/// for the interface-conformity limitation and the 3-D CDT upgrade path).
pub mod multi_material;
pub use multi_material::{
    MATERIAL_TAG_OUTER, MultiMaterialError, region_volumes, tetrahedralize_multi_material,
    tetrahedralize_region,
};
