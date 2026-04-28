use glam::DVec3;
use serde::{Deserialize, Serialize};

/// Geometric (analytic) model types: position, curve, surface, primitive descriptors.
///
/// This module describes *what shape is*.
pub mod geom;

/// 装配体：层级结构、组件实例化、世界变换展平。
///
/// 类比 OCCT `XCAFDoc_ShapeTool` 管理的 shape 层级。
pub mod assembly;

/// Topology model types: vertex/edge/face/shell/solid incidence relationships.
///
/// This module describes *how things are connected*.
pub mod topology;

/// Shape properties: surface area, volume, centroid.
///
/// Analogous to OCCT `GProp_GProps` + `BRepGProp`.
pub mod properties;

/// Topology query helpers: edge adjacency, vertex adjacency, shape counts.
///
/// Analogous to OCCT `TopExp_Explorer` and `TopExp::MapShapesAndAncestors`.
pub mod topo_query;

/// Cached graph-topology wrapper with O(1) adjacency, DFS/BFS traversal,
/// and mutation-dirty tracking.
///
/// Analogous to OCCT `BRepGraph` module (new in OCCT 7.7+).
pub mod brep_graph;

/// Persistent naming hooks for stable user-level topology labels.
///
/// Analogous to OCCT OCAF/TopoNaming-style name tables.
pub mod naming;

/// Persistent naming semantics for BRepGraph topology entities.
///
/// Provides stable, operation-surviving identifiers for topology entities.
pub mod persistent_naming;

/// Differential geometry: principal curvatures, Gaussian curvature, mean curvature.
///
/// Analogous to OCCT `GeomLProp_SLProps`.
pub mod curvature;

/// Curve arc-length computation.
///
/// Analogous to OCCT `GCPnts_AbscissaPoint` / `CPnts_AbscissaPoint::Length`.
pub mod arc_length;

/// Visual appearance: per-face/solid RGB color and basic material.
///
/// Analogous to OCCT `XCAFDoc_ColorTool`.
pub mod appearance;

/// Dimension and tolerance object model (GD&T).
///
/// Analogous to OCCT `XCAFDimTolObjects`.
pub mod dim_tol;

/// Annotation object model for CAD annotations (PMI).
///
/// Analogous to OCCT `XCAFNoteObjects`.
pub mod annotation;

/// Precision constants and per-entity tolerance query helpers.
///
/// Analogous to OCCT `Precision` class and `BRep_Tool::Tolerance`.
pub mod tolerance;

/// Curve fitting: B-spline interpolation and approximation through point sets.
///
/// Analogous to OCCT `GeomAPI_Interpolate` and `GeomAPI_PointsToBSpline`.
pub mod fit;

/// Closest-point projection from a 3D point onto a curve or surface.
///
/// Analogous to OCCT `GeomAPI_ProjectPointOnCurve` and
/// `GeomAPI_ProjectPointOnSurf`.
pub mod projection;

/// Shape-to-shape and point-to-shape minimum distance.
///
/// Analogous to OCCT `BRepExtrema_DistShapeShape`.
pub mod distance;

/// Curve-curve extrema: find (s,t) minimising |C1(s) − C2(t)|.
///
/// Analogous to OCCT `GeomAPI_ExtremaCurveCurve`.
pub mod extrema;

/// NURBS interoperability: convert analytic curves/surfaces to BSpline.
///
/// Analogous to OCCT `GeomConvert::CurveToBSplineCurve` /
/// `GeomConvert::SurfaceToBSplineSurface`.
pub mod nurbs_convert;

/// Curve and surface trimming and extension.
///
/// Analogous to OCCT `Geom_TrimmedCurve` construction helpers,
/// `GeomAPI_ExtendCurveToPoint`, and `Geom_RectangularTrimmedSurface`.
pub mod extend;

pub use distance::{ShapeDistance, min_distance, point_to_shape_distance};
pub use extend::{
    CurveEnd, SurfaceBoundary, extend_bspline_surface, extend_curve_by_length,
    extend_curve_to_point, insert_knot_to_multiplicity, trim_curve, trim_surface,
};
pub use extrema::{CurveCurveExtrema, ExtremaPair, extrema_curve_curve};
pub use fit::{FitError, approximate_points, interpolate_points, interpolate_points_2d};
pub use nurbs_convert::{
    bezier_curve_to_bspline, bezier_surface_to_bspline, circle_to_bspline, curve_to_bspline,
    cylinder_to_bspline, ellipse_to_bspline, line_to_bspline, line_to_bspline_range,
    plane_to_bspline, sphere_to_bspline, surface_to_bspline,
};
pub use projection::{
    CurveProjection, SurfaceProjection, closest_point_on_curve, closest_point_on_surface,
};

pub use appearance::{Color, FaceColor, StepColor};
pub use dim_tol::{
    DimensionType, GeometricToleranceType, ToleranceModifier,
    DimensionalTolerance, GeometricToleranceObject, DatumReference, DatumSystem, DimTolStore,
};
pub use annotation::{
    NoteType, ArrowType, WeldType,
    AnnotationNote, TextAnnotation, LeaderLine, SurfaceTextureSymbol,
    WeldSymbol, BalloonAnnotation, AnnotationStore,
    Annotation, AnnotationKind, Note, NoteCategory, NoteTarget, View, ViewProjection,
};
pub use arc_length::arc_length;
pub use curvature::{gaussian_curvature, mean_curvature, principal_curvatures};
pub use geom::{Point3, Vec3, Point2, Vec2};
pub use geom::PrimitiveSolid;
pub use geom::TrimmedSurface;
pub use geom::{
    ArchimedeanSpiral2d, BSplineCurve2, CircleInvolute2d, Ellipse2d, LogarithmicSpiral2d,
    SineWave2d,
};
pub use geom::{BSplineSurface, CoonsSurface, LinearExtrusionSurface, RevolutionSurface, RuledSurface};
pub use geom::{BezierCurve2, BezierCurve3, BezierSurface, TriBezierSurface};
pub use geom::{Curve2d, Curve3, Surface3};
pub use geom::{Curve2dEval, CurveEval, SurfaceEval, any_perpendicular};
pub use geom::{EllipsoidalSurface, HelicoidSurface, PipeSurface};
pub use geom::{CircularHelix3, Hyperbola3, Parabola3, SineWave3};
pub use geom::{OffsetCurve3, OffsetSurface};
pub use properties::{InertiaTensor, centroid, inertia_tensor, surface_area, volume};
pub use tolerance::{
    ANGULAR, APPROXIMATION, CONFUSION, edge_same_parameter, edge_same_range, edge_tolerance,
    face_domain, face_tolerance, model_tolerance, vertex_tolerance,
    resize_tolerance_arrays,
    set_vertex_tolerance, update_vertex_tolerance,
    set_edge_tolerance, update_edge_tolerance,
    set_face_tolerance, update_face_tolerance,
    finalize_tolerance_hierarchy,
};
pub use topo_query::{
    edge_adjacent_faces, edge_count, face_count, face_edges, is_degenerate_edge,
    seam_edge_candidates, vertex_adjacent_edges, vertex_count,
};
pub use brep_graph::{
    BfsFaces, BRepGraph, BRepGraphBuilder, BRepGraphCheckpointData,
    BRepGraphTool, DfsEdgesFromVertex, DfsFaces,
    ManifoldRepairHints, NonManifoldSummary, RepairHint,
};
pub use naming::{PersistentNamingHooks, TopoEntityRef};
pub use persistent_naming::{
    ConflictResolution, CrossOperationHistory, CrossOperationStabilityReport, EntityGenealogy,
    EntityType, EntityTypeStability, IssueSeverity, NamingConflictResolution, NamingContext,
    NamingEvent, NamingHistory, NamingIssue, NamingRule, NamingStabilityReport,
    NamePropagationPolicy, OperationId, OperationRecord, OperationStats, OperationType,
    PersistentId, PersistentNamingEngine, PersistentNamingHooksExt,
};
pub use topology::{Compound, CompSolid, Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

/// A parameter-space curve binding that ties a 3D edge to an adjacent face's
/// surface parameter domain (u, v).  Analogous to OCCT `BRep_CurveOnSurface`.
///
/// `surface_idx` indexes into `GeomStore.surfaces`.
/// `curve2d_idx` indexes into `GeomStore.curve2ds`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PCurve {
    pub surface_idx: usize,
    pub curve2d_idx: usize,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GeomStore {
    /// Pool of 3D analytic curves.
    pub curves: Vec<Curve3>,
    /// Pool of analytic surfaces.
    pub surfaces: Vec<Surface3>,
    /// Pool of 2D parameter-space curves used by PCurves.
    pub curve2ds: Vec<Curve2d>,
    /// Indexed by `BRep.edges` index; value is index into `curves`.
    pub edge_curve: Vec<Option<usize>>,
    /// Flattened face order across solids/shells; value is index into `surfaces`.
    pub face_surface: Vec<Option<usize>>,
    /// Indexed by `BRep.edges` index; each entry is the list of PCurves for
    /// that edge on its adjacent faces (usually 1, seam edges have 2).
    pub edge_pcurves: Vec<Vec<PCurve>>,
    /// Parallel to `edge_curve`: the parameter range [t1, t2] of the edge on its
    /// 3D curve. `None` = unknown (algorithms fall back to `CurveEval::default_domain`).
    /// Analogous to `BRep_Edge::Range()` in OCCT.
    #[serde(default)]
    pub edge_curve_range: Vec<Option<[f64; 2]>>,
    /// Parallel to `BRep.edges`: `true` if this is a degenerate edge (zero-length,
    /// e.g. a polar singularity). Analogous to `BRep_Edge::Degenerated()` in OCCT.
    #[serde(default)]
    pub edge_degenerated: Vec<bool>,
    /// Per-vertex tolerance (falls back to `tolerance::CONFUSION` when absent or zero).
    /// Parallel to `BRep.vertices`. Analogous to `BRep_Tool::Tolerance(vertex)` in OCCT.
    #[serde(default)]
    pub vertex_tolerance: Vec<f64>,
    /// Per-edge tolerance (falls back to `tolerance::CONFUSION` when absent or zero).
    /// Parallel to `BRep.edges`. Analogous to `BRep_Tool::Tolerance(edge)` in OCCT.
    #[serde(default)]
    pub edge_tolerance: Vec<f64>,
    /// Per-face tolerance (falls back to `tolerance::CONFUSION` when absent or zero).
    /// Parallel to the flattened face order (same indexing as `face_surface`).
    /// Analogous to `BRep_Tool::Tolerance(face)` in OCCT.
    #[serde(default)]
    pub face_tolerance: Vec<f64>,
    /// Per-curve2d parameter range [t1, t2].
    ///
    /// Used when the PCurve originates from a STEP `TRIMMED_CURVE` entity in
    /// 2D parameter space. `None` means the natural domain of the curve is used.
    /// Parallel to `GeomStore.curve2ds`. Analogous to `edge_curve_range` for 3D.
    #[serde(default)]
    pub curve2d_range: Vec<Option<[f64; 2]>>,
    /// Per-face surface parameter domain override [u1, u2, v1, v2].
    ///
    /// When populated (e.g. from a STEP `RECTANGULAR_TRIMMED_SURFACE`), the face
    /// is restricted to this subdomain of its underlying surface. `None` means
    /// `SurfaceEval::default_domain()` is used. Parallel to `face_surface`.
    /// Analogous to `edge_curve_range` for 3D curves.
    #[serde(default)]
    pub face_surface_range: Vec<Option<[f64; 4]>>,
    /// Per-edge SameParameter flag.
    ///
    /// `true` if the 3D curve and all PCurves share the same parameterization
    /// (i.e. the parameter `t` on the 3D curve maps directly to the same `t`
    /// on every PCurve). Analogous to `BRep_Edge::SameParameter()` in OCCT.
    /// When absent or empty, assumed `true` for analytic primitives we generate.
    #[serde(default)]
    pub edge_same_parameter: Vec<bool>,
    /// Per-edge SameRange flag.
    ///
    /// `true` if all PCurves on this edge share the same `[t1, t2]` parameter
    /// range as the 3D curve's `edge_curve_range`. Analogous to
    /// `BRep_Edge::SameRange()` in OCCT.
    /// When absent or empty, assumed `true` for analytic primitives we generate.
    #[serde(default)]
    pub edge_same_range: Vec<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub solids: Vec<Solid>,
    #[serde(default)]
    pub geom: GeomStore,
    /// Optional compound container for multi-shape assemblies.
    ///
    /// When set, this BRep represents a compound shape. The `solids` field
    /// contains flattened solids for backward compatibility.
    #[serde(default)]
    pub compound: Option<topology::Compound>,
    /// Optional CompSolid container for connected multi-region solids.
    #[serde(default)]
    pub compsolid: Option<topology::CompSolid>,
}

impl Default for BRep {
    fn default() -> Self {
        Self::new()
    }
}

impl BRep {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            solids: Vec::new(),
            geom: GeomStore::default(),
            compound: None,
            compsolid: None,
        }
    }

    /// Create a BRep representing a compound of shapes.
    ///
    /// The compound's solids are also flattened into `self.solids` for
    /// backward compatibility with code that expects `brep.solids`.
    pub fn from_compound(compound: topology::Compound) -> Self {
        let solids = compound.flatten_solids().into_iter().cloned().collect();
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            solids,
            geom: GeomStore::default(),
            compound: Some(compound),
            compsolid: None,
        }
    }

    /// Create a BRep representing a CompSolid (connected multi-region solid).
    pub fn from_compsolid(compsolid: topology::CompSolid) -> Self {
        let solids = compsolid.solids.clone();
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            solids,
            geom: GeomStore::default(),
            compound: None,
            compsolid: Some(compsolid),
        }
    }

    /// Create a compound from multiple BReps.
    ///
    /// Each input BRep's solids are extracted and added to the compound.
    pub fn compound_from_shapes(shapes: &[BRep]) -> BRep {
        let mut compound = topology::Compound::new();
        for shape in shapes {
            for solid in &shape.solids {
                compound.add_solid(None, solid.clone());
            }
        }
        Self::from_compound(compound)
    }

    /// Explode this BRep into constituent shapes.
    ///
    /// If this is a compound, returns one BRep per top-level shape.
    /// If not a compound, returns a single-element Vec containing self.
    pub fn explode(&self) -> Vec<BRep> {
        if let Some(ref compound) = self.compound {
            let mut result = Vec::new();

            // Explode solids
            for (_, solid) in &compound.solids {
                let mut brep = BRep::new();
                brep.solids.push(solid.clone());
                result.push(brep);
            }

            // Explode comp_solids
            for (_, cs) in &compound.comp_solids {
                result.push(BRep::from_compsolid(cs.clone()));
            }

            // Explode shells
            for (_, shell) in &compound.shells {
                let mut brep = BRep::new();
                brep.solids.push(topology::Solid {
                    shells: vec![shell.clone()],
                });
                result.push(brep);
            }

            // Explode nested compounds
            for (_, nested) in &compound.compounds {
                result.push(BRep::from_compound(nested.clone()));
            }

            result
        } else if let Some(ref cs) = self.compsolid {
            // Explode CompSolid into individual solids
            cs.solids
                .iter()
                .map(|solid| {
                    let mut brep = BRep::new();
                    brep.solids.push(solid.clone());
                    brep
                })
                .collect()
        } else {
            vec![self.clone()]
        }
    }

    /// Add a shape to this BRep's compound.
    ///
    /// If this BRep is not already a compound, it will be converted to one.
    pub fn add_shape(&mut self, shape: BRep) {
        // Ensure we have a compound
        if self.compound.is_none() {
            let mut compound = topology::Compound::new();
            // Move existing solids into compound
            for solid in std::mem::take(&mut self.solids) {
                compound.add_solid(None, solid);
            }
            self.compound = Some(compound);
        }

        // Add the shape
        if let Some(ref mut compound) = self.compound {
            for solid in &shape.solids {
                compound.add_solid(None, solid.clone());
                self.solids.push(solid.clone());
            }
        }
    }

    /// Remove a shape from this BRep's compound by index.
    ///
    /// Returns `true` if a shape was removed.
    pub fn remove_shape(&mut self, index: usize) -> bool {
        if let Some(ref mut compound) = self.compound {
            if index < compound.solids.len() {
                compound.remove_solid(index);
                // Rebuild flattened solids
                self.solids = compound.flatten_solids().into_iter().cloned().collect();
                return true;
            }
        }
        false
    }

    /// Returns `true` if this BRep represents a compound.
    pub fn is_compound(&self) -> bool {
        self.compound.is_some()
    }

    /// Returns `true` if this BRep represents a CompSolid.
    pub fn is_compsolid(&self) -> bool {
        self.compsolid.is_some()
    }

    /// Get the compound if this BRep is one.
    pub fn as_compound(&self) -> Option<&topology::Compound> {
        self.compound.as_ref()
    }

    /// Get the CompSolid if this BRep is one.
    pub fn as_compsolid(&self) -> Option<&topology::CompSolid> {
        self.compsolid.as_ref()
    }

    /// Iterate over all solids, including those in compounds and compsolids.
    pub fn iter_solids(&self) -> impl Iterator<Item = &topology::Solid> {
        self.solids.iter()
    }

    /// Flatten all solids from this BRep.
    ///
    /// If this is a compound, returns all nested solids.
    /// Otherwise returns the direct solids.
    pub fn flatten_to_solids(&self) -> Vec<&topology::Solid> {
        if let Some(ref compound) = self.compound {
            compound.flatten_solids()
        } else {
            self.solids.iter().collect()
        }
    }

    /// Creates a unit box B-Rep.
    ///
    /// Vertex layout:
    ///   0:(0,0,0)  1:(w,0,0)  2:(w,h,0)  3:(0,h,0)   <- front face (z=0)
    ///   4:(0,0,d)  5:(w,0,d)  6:(w,h,d)  7:(0,h,d)   <- back face  (z=d)
    fn create_box(width: f64, height: f64, depth: f64) -> Self {
        let (w, h, d) = (width, height, depth);

        let vertices = vec![
            Vertex {
                point: DVec3::new(0.0, 0.0, 0.0),
            }, // 0
            Vertex {
                point: DVec3::new(w, 0.0, 0.0),
            }, // 1
            Vertex {
                point: DVec3::new(w, h, 0.0),
            }, // 2
            Vertex {
                point: DVec3::new(0.0, h, 0.0),
            }, // 3
            Vertex {
                point: DVec3::new(0.0, 0.0, d),
            }, // 4
            Vertex {
                point: DVec3::new(w, 0.0, d),
            }, // 5
            Vertex {
                point: DVec3::new(w, h, d),
            }, // 6
            Vertex {
                point: DVec3::new(0.0, h, d),
            }, // 7
        ];

        // 12 edges: 4 front + 4 back + 4 lateral
        let edges = vec![
            Edge { start: 0, end: 1 }, // 0  front-bottom
            Edge { start: 1, end: 2 }, // 1  front-right
            Edge { start: 2, end: 3 }, // 2  front-top
            Edge { start: 3, end: 0 }, // 3  front-left
            Edge { start: 4, end: 5 }, // 4  back-bottom
            Edge { start: 5, end: 6 }, // 5  back-right
            Edge { start: 6, end: 7 }, // 6  back-top
            Edge { start: 7, end: 4 }, // 7  back-left
            Edge { start: 0, end: 4 }, // 8  lateral-bl
            Edge { start: 1, end: 5 }, // 9  lateral-br
            Edge { start: 2, end: 6 }, // 10 lateral-tr
            Edge { start: 3, end: 7 }, // 11 lateral-tl
        ];

        let faces = vec![
            // Front  (z=0, normal -Z)
            Face {
                outer_wire: Wire {
                    edges: vec![
                        WireEdge::fwd(0),
                        WireEdge::fwd(1),
                        WireEdge::fwd(2),
                        WireEdge::fwd(3),
                    ],
                },
                inner_wires: vec![],
                normal: DVec3::new(0.0, 0.0, -1.0),
                triangles: vec![[0, 1, 2], [0, 2, 3]],
                mesh_dirty: false,
            },
            // Back   (z=d, normal +Z)
            Face {
                outer_wire: Wire {
                    edges: vec![
                        WireEdge::fwd(4),
                        WireEdge::fwd(5),
                        WireEdge::fwd(6),
                        WireEdge::fwd(7),
                    ],
                },
                inner_wires: vec![],
                normal: DVec3::new(0.0, 0.0, 1.0),
                triangles: vec![[5, 4, 7], [5, 7, 6]],
                mesh_dirty: false,
            },
            // Bottom (y=0, normal -Y)
            Face {
                outer_wire: Wire {
                    edges: vec![
                        WireEdge::fwd(0),
                        WireEdge::fwd(9),
                        WireEdge::rev(4),
                        WireEdge::rev(8),
                    ],
                },
                inner_wires: vec![],
                normal: DVec3::new(0.0, -1.0, 0.0),
                triangles: vec![[0, 1, 5], [0, 5, 4]],
                mesh_dirty: false,
            },
            // Top    (y=h, normal +Y)
            Face {
                outer_wire: Wire {
                    edges: vec![
                        WireEdge::rev(2),
                        WireEdge::fwd(10),
                        WireEdge::fwd(6),
                        WireEdge::rev(11),
                    ],
                },
                inner_wires: vec![],
                normal: DVec3::new(0.0, 1.0, 0.0),
                triangles: vec![[3, 2, 6], [3, 6, 7]],
                mesh_dirty: false,
            },
            // Left   (x=0, normal -X)
            Face {
                outer_wire: Wire {
                    edges: vec![
                        WireEdge::rev(3),
                        WireEdge::fwd(11),
                        WireEdge::fwd(7),
                        WireEdge::rev(8),
                    ],
                },
                inner_wires: vec![],
                normal: DVec3::new(-1.0, 0.0, 0.0),
                triangles: vec![[0, 3, 7], [0, 7, 4]],
                mesh_dirty: false,
            },
            // Right  (x=w, normal +X)
            Face {
                outer_wire: Wire {
                    edges: vec![
                        WireEdge::fwd(1),
                        WireEdge::fwd(10),
                        WireEdge::rev(5),
                        WireEdge::rev(9),
                    ],
                },
                inner_wires: vec![],
                normal: DVec3::new(1.0, 0.0, 0.0),
                triangles: vec![[1, 2, 6], [1, 6, 5]],
                mesh_dirty: false,
            },
        ];

        // Populate GeomStore with planes (surfaces) for face geometry
        // Note: We populate surfaces but NOT curves/pcurves for the box.
        // The box's edges are implicit line segments defined by vertex positions.
        // This allows primitive detection (is_box) to work without requiring
        // geometrically accurate pcurves.
        use geom::{Plane, Surface3};

        // Create 6 plane surfaces for faces
        // F0: Front (z=0, normal -Z), F1: Back (z=d, normal +Z)
        // F2: Bottom (y=0, normal -Y), F3: Top (y=h, normal +Y)
        // F4: Left (x=0, normal -X), F5: Right (x=w, normal +X)
        let surfaces: Vec<Surface3> = vec![
            Surface3::Plane(Plane { origin: DVec3::new(0.0, 0.0, 0.0), normal: -DVec3::Z }), // Front
            Surface3::Plane(Plane { origin: DVec3::new(0.0, 0.0, d), normal: DVec3::Z }),   // Back
            Surface3::Plane(Plane { origin: DVec3::new(0.0, 0.0, 0.0), normal: -DVec3::Y }), // Bottom
            Surface3::Plane(Plane { origin: DVec3::new(0.0, h, 0.0), normal: DVec3::Y }),    // Top
            Surface3::Plane(Plane { origin: DVec3::new(0.0, 0.0, 0.0), normal: -DVec3::X }), // Left
            Surface3::Plane(Plane { origin: DVec3::new(w, 0.0, 0.0), normal: DVec3::X }),    // Right
        ];

        // Face surface indices
        let face_surface: Vec<Option<usize>> = (0..6).map(|i| Some(i)).collect();

        let geom = GeomStore {
            curves: Vec::new(),
            surfaces,
            curve2ds: Vec::new(),
            edge_curve: Vec::new(),
            face_surface,
            edge_pcurves: Vec::new(),
            edge_curve_range: Vec::new(),
            edge_degenerated: Vec::new(),
            vertex_tolerance: Vec::new(),
            edge_tolerance: Vec::new(),
            face_tolerance: Vec::new(),
            curve2d_range: Vec::new(),
            face_surface_range: Vec::new(),
            edge_same_parameter: Vec::new(),
            edge_same_range: Vec::new(),
        };

        BRep {
            vertices,
            edges,
            solids: vec![Solid {
                shells: vec![Shell { faces }],
            }],
            geom,
            compound: None,
            compsolid: None,
        }
    }
    ///
    /// Topology (OCCT-compatible single-seam representation):
    ///   Vertices: north (0, r, 0), south (0, -r, 0)
    ///   Edge E0:  seam meridian (north → south), Circle3 in XZ plane (normal = +Z)
    ///   Face F0:  SphericalSurface, outer_wire = [E0_fwd, E0_rev] (seam edge repeated)
    ///   PCurves:  E0 forward  → Line2d u=0,  v: 0 → π
    ///             E0 reversed → Line2d u=2π, v: π → 0
    fn create_sphere(radius: f64) -> Self {
        use geom::*;
        use std::f64::consts::PI;

        let r = radius;
        // Vertices
        let north = DVec3::new(0.0, r, 0.0);
        let south = DVec3::new(0.0, -r, 0.0);
        let vertices = vec![Vertex { point: north }, Vertex { point: south }];

        // Edge E0: seam meridian (north→south) — Circle3 in XZ plane
        // The seam lies at theta=0, i.e. x>0, z=0 plane.
        let edges = vec![Edge { start: 0, end: 1 }]; // E0

        // Face F0: outer_wire uses E0 twice (forward then reversed = seam)
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::rev(0)],
            },
            inner_wires: vec![],
            normal: DVec3::X, // outward, approximate
            triangles: vec![],
            mesh_dirty: true,
        };
        let shell = Shell { faces: vec![face] };
        let solid = Solid {
            shells: vec![shell],
        };

        // GeomStore
        let seam_curve = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: r,
        });
        let sphere_surf = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: r,
        });
        // PCurves for the seam edge on the sphere:
        //   Sphere param: u = longitude [0, 2π], v = colatitude [0, π] (phi from north pole)
        //   Forward half (north→south at u=0): Line2d origin=(0,0) dir=(0,1) extent π
        //   Reversed half (south→north at u=2π): Line2d origin=(2π,π) dir=(0,-1) extent π
        let pc_fwd = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(0.0, 1.0),
        });
        let pc_rev = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(2.0 * PI, PI),
            direction: glam::DVec2::new(0.0, -1.0),
        });
        let geom = GeomStore {
            curves: vec![seam_curve],
            surfaces: vec![sphere_surf],
            curve2ds: vec![pc_fwd, pc_rev],
            edge_curve: vec![Some(0)],   // E0 → seam_curve
            face_surface: vec![Some(0)], // F0 → sphere_surf
            edge_pcurves: vec![vec![
                // E0: two pcurves (fwd and rev sides of seam)
                PCurve {
                    surface_idx: 0,
                    curve2d_idx: 0,
                },
                PCurve {
                    surface_idx: 0,
                    curve2d_idx: 1,
                },
            ]],
            // E0 is the half-meridian: t ∈ [0, π] on Circle3 (north→south)
            edge_curve_range: vec![Some([0.0, PI])],
            edge_degenerated: vec![false],
            vertex_tolerance: Vec::new(),
            edge_tolerance: Vec::new(),
            face_tolerance: Vec::new(),
            curve2d_range: Vec::new(),
            face_surface_range: Vec::new(),
            edge_same_parameter: Vec::new(),
            edge_same_range: Vec::new(),
        };

        Self {
            vertices,
            edges,
            solids: vec![solid],
            geom,
            compound: None,
            compsolid: None,
        }
    }

    /// Creates an analytic cylinder BRep along +Y axis, centered at origin.
    ///
    /// Topology:
    ///   Vertices: top_p (r, h/2, 0), bot_p (r, -h/2, 0)
    ///   Edges:
    ///     E0: top circle (Circle3, top_p → top_p seam, center=(0,h/2,0), normal=+Y)
    ///     E1: bot circle (Circle3, bot_p → bot_p seam, center=(0,-h/2,0), normal=-Y)
    ///     E2: seam line  (Line3,   top_p → bot_p)
    ///   Faces:
    ///     F0: CylindricalSurface, wire=[E2, E1_rev, E2_rev, E0]
    ///     F1: Plane +Y cap,       wire=[E0]
    ///     F2: Plane -Y cap,       wire=[E1_rev]  (stored as E1 with wire handling orientation)
    fn create_cylinder(radius: f64, height: f64) -> Self {
        use geom::*;
        use std::f64::consts::PI;

        let r = radius;
        let h = height;
        let half_h = h * 0.5;

        let top_p = DVec3::new(r, half_h, 0.0);
        let bot_p = DVec3::new(r, -half_h, 0.0);
        let vertices = vec![Vertex { point: top_p }, Vertex { point: bot_p }];

        let edges = vec![
            Edge { start: 0, end: 0 },
            Edge { start: 1, end: 1 },
            Edge { start: 0, end: 1 },
        ];

        let f0 = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(2),
                    WireEdge::rev(1),
                    WireEdge::rev(2),
                    WireEdge::fwd(0),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::X,
            triangles: vec![],
            mesh_dirty: true,
        };
        let f1 = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0)],
            },
            inner_wires: vec![],
            normal: DVec3::Y,
            triangles: vec![],
            mesh_dirty: true,
        };
        let f2 = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::rev(1)],
            },
            inner_wires: vec![],
            normal: -DVec3::Y,
            triangles: vec![],
            mesh_dirty: true,
        };
        let shell = Shell {
            faces: vec![f0, f1, f2],
        };
        let solid = Solid {
            shells: vec![shell],
        };

        let top_circle = Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, half_h, 0.0),
            normal: DVec3::Y,
            radius: r,
        });
        let bot_circle = Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, -half_h, 0.0),
            normal: -DVec3::Y,
            radius: r,
        });
        let seam_line = Curve3::Line(Line3 {
            origin: top_p,
            direction: (bot_p - top_p).normalize(),
        });

        let cyl_surf = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(0.0, -half_h, 0.0),
            axis: DVec3::Y,
            radius: r,
        });
        let top_plane = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, half_h, 0.0),
            normal: DVec3::Y,
        });
        let bot_plane = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, -half_h, 0.0),
            normal: -DVec3::Y,
        });

        let e0_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(2.0 * PI, h),
            direction: glam::DVec2::new(-1.0, 0.0),
        });
        let e0_on_f1 = Curve2d::Circle(Circle2d {
            center: glam::DVec2::ZERO,
            radius: r,
        });
        let e1_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(1.0, 0.0),
        });
        let e1_on_f2 = Curve2d::Circle(Circle2d {
            center: glam::DVec2::ZERO,
            radius: r,
        });
        let e2_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, h),
            direction: glam::DVec2::new(0.0, -1.0),
        });

        let geom = GeomStore {
            curves: vec![top_circle, bot_circle, seam_line],
            surfaces: vec![cyl_surf, top_plane, bot_plane],
            curve2ds: vec![e0_on_f0, e0_on_f1, e1_on_f0, e1_on_f2, e2_on_f0],
            edge_curve: vec![Some(0), Some(1), Some(2)],
            face_surface: vec![Some(0), Some(1), Some(2)],
            edge_pcurves: vec![
                vec![
                    PCurve {
                        surface_idx: 0,
                        curve2d_idx: 0,
                    },
                    PCurve {
                        surface_idx: 1,
                        curve2d_idx: 1,
                    },
                ],
                vec![
                    PCurve {
                        surface_idx: 0,
                        curve2d_idx: 2,
                    },
                    PCurve {
                        surface_idx: 2,
                        curve2d_idx: 3,
                    },
                ],
                vec![PCurve {
                    surface_idx: 0,
                    curve2d_idx: 4,
                }],
            ],
            edge_curve_range: vec![
                Some([0.0, 2.0 * PI]),
                Some([0.0, 2.0 * PI]),
                Some([0.0, h]),
            ],
            edge_degenerated: vec![false, false, false],
            vertex_tolerance: Vec::new(),
            edge_tolerance: Vec::new(),
            face_tolerance: Vec::new(),
            curve2d_range: Vec::new(),
            face_surface_range: Vec::new(),
            edge_same_parameter: Vec::new(),
            edge_same_range: Vec::new(),
        };

        Self {
            vertices,
            edges,
            solids: vec![solid],
            geom,
            compound: None,
            compsolid: None,
        }
    }

    /// Creates an analytic cone BRep along +Y axis, apex at +Y, centered at origin.
    ///
    /// Topology:
    ///   Vertices: apex (0, h/2, 0), base_p (R, -h/2, 0)
    ///   Edges:
    ///     E0: base circle (Circle3, base_p → base_p seam, normal=-Y)
    ///     E1: slant line  (Line3,   apex → base_p)
    ///   Faces:
    ///     F0: ConicalSurface, wire=[E1, E0, E1_rev]  (seam)
    ///     F1: Plane -Y cap,  wire=[E0]
    fn create_cone(base_radius: f64, height: f64) -> Self {
        use geom::*;
        use std::f64::consts::PI;
        let r = base_radius;
        let h = height;
        let half_h = h * 0.5;

        let apex_pt = DVec3::new(0.0, half_h, 0.0);
        let base_pt = DVec3::new(r, -half_h, 0.0);
        let vertices = vec![Vertex { point: apex_pt }, Vertex { point: base_pt }];

        let edges = vec![
            Edge { start: 1, end: 1 }, // E0 base circle seam
            Edge { start: 0, end: 1 }, // E1 slant line apex→base_p
        ];

        // F0 conical lateral: E1 fwd (apex→base seam), E0 fwd (base circle), E1 rev (base→apex seam)
        let f0 = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(1), WireEdge::fwd(0), WireEdge::rev(1)],
            },
            inner_wires: vec![],
            normal: DVec3::X,
            triangles: vec![],
            mesh_dirty: true,
        };
        // F1 base cap: E0 reversed (base circle CW from -Y view)
        let f1 = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::rev(0)],
            },
            inner_wires: vec![],
            normal: -DVec3::Y,
            triangles: vec![],
            mesh_dirty: true,
        };
        let shell = Shell {
            faces: vec![f0, f1],
        };
        let solid = Solid {
            shells: vec![shell],
        };

        // Curves
        let base_circle = Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, -half_h, 0.0),
            normal: -DVec3::Y,
            radius: r,
        });
        let slant = Curve3::Line(Line3 {
            origin: apex_pt,
            direction: (base_pt - apex_pt).normalize(),
        });

        // half-angle = atan(R / h)
        let half_angle = (r / h).atan();
        // ConicalSurface: apex at top, axis pointing down (-Y), radius at apex = 0
        let cone_surf = Surface3::Cone(geom::ConicalSurface {
            apex: apex_pt,
            axis: -DVec3::Y,
            radius: 0.0, // radius at apex
            half_angle_rad: half_angle,
        });
        let base_plane = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, -half_h, 0.0),
            normal: -DVec3::Y,
        });

        // PCurves
        // Cone param: u=azimuth [0,2π], v=slant distance from apex
        let slant_len = (r * r + h * h).sqrt();
        // E0 (base circle) on F0 (cone): iso-line at v=slant_len, u from 0 to 2π
        let e0_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, slant_len),
            direction: glam::DVec2::new(1.0, 0.0),
        });
        let e0_on_f1 = Curve2d::Circle(Circle2d {
            center: glam::DVec2::ZERO,
            radius: r,
        });
        // E1 (slant seam) on F0: iso-line u=0, v from 0 to slant_len
        let e1_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(0.0, 1.0),
        });

        let geom = GeomStore {
            curves: vec![base_circle, slant],
            surfaces: vec![cone_surf, base_plane],
            curve2ds: vec![e0_on_f0, e0_on_f1, e1_on_f0],
            edge_curve: vec![Some(0), Some(1)],
            face_surface: vec![Some(0), Some(1)],
            edge_pcurves: vec![
                // E0: base circle → on cone and base plane
                vec![
                    PCurve {
                        surface_idx: 0,
                        curve2d_idx: 0,
                    },
                    PCurve {
                        surface_idx: 1,
                        curve2d_idx: 1,
                    },
                ],
                // E1: slant line → on cone only
                vec![PCurve {
                    surface_idx: 0,
                    curve2d_idx: 2,
                }],
            ],
            // E0: full base circle [0, 2π]; E1: slant from apex to base [0, slant_len]
            edge_curve_range: vec![
                Some([0.0, 2.0 * PI]),  // E0 base circle
                Some([0.0, slant_len]), // E1 slant line
            ],
            edge_degenerated: vec![false, false],
            vertex_tolerance: Vec::new(),
            edge_tolerance: Vec::new(),
            face_tolerance: Vec::new(),
            curve2d_range: Vec::new(),
            face_surface_range: Vec::new(),
            edge_same_parameter: Vec::new(),
            edge_same_range: Vec::new(),
        };

        Self {
            vertices,
            edges,
            solids: vec![solid],
            geom,
            compound: None,
            compsolid: None,
        }
    }

    /// Creates an analytic torus BRep around +Y axis, centered at origin.
    ///
    /// Topology (double-seam representation):
    ///   Vertex: seam_pt = (R+r, 0, 0) — intersection of major and minor seams
    ///   Edges:
    ///     E0: major seam circle (major radius R, center=origin, normal=+Y)
    ///     E1: minor seam circle (minor radius r, in XY plane at (R,0,0))
    ///   Face:
    ///     F0: ToroidalSurface, outer_wire=[E0, E1, E0_rev, E1_rev]
    ///   PCurves:
    ///     E0 on F0: Line2d (0,0)→(2π,0)   [major seam: u from 0..2π at v=0]
    ///     E1 on F0: Line2d (0,0)→(0,2π)   [minor seam: v from 0..2π at u=0]
    fn create_torus(major_radius: f64, minor_radius: f64) -> Self {
        use geom::*;
        use std::f64::consts::PI;

        let big_r = major_radius;
        let small_r = minor_radius;

        // Single vertex at the seam intersection
        let seam_pt = DVec3::new(big_r + small_r, 0.0, 0.0);
        let vertices = vec![Vertex { point: seam_pt }];

        // E0: major seam (full circle of radius R in XZ plane)
        // E1: minor seam (full circle of radius r in a plane through the tube)
        let edges = vec![
            Edge { start: 0, end: 0 }, // E0 major circle seam
            Edge { start: 0, end: 0 }, // E1 minor circle seam
        ];

        // F0: outer_wire — E0 fwd (major seam), E1 fwd (minor seam),
        //     E0 rev (major seam reversed), E1 rev (minor seam reversed)
        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(1),
                    WireEdge::rev(0),
                    WireEdge::rev(1),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::X,
            triangles: vec![],
            mesh_dirty: true,
        };
        let shell = Shell { faces: vec![face] };
        let solid = Solid {
            shells: vec![shell],
        };

        // Curves
        let major_circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Y,
            radius: big_r,
        });
        // Minor circle: centered at (R,0,0), in the YZ plane (normal = +X)
        let minor_circle = Curve3::Circle(Circle3 {
            center: DVec3::new(big_r, 0.0, 0.0),
            normal: DVec3::X,
            radius: small_r,
        });

        let torus_surf = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: big_r,
            minor_radius: small_r,
        });

        // PCurves — torus param: u=major angle [0,2π], v=minor angle [0,2π]
        let e0_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(1.0, 0.0),
        });
        let e0_on_f0_rev = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(2.0 * PI, 0.0),
            direction: glam::DVec2::new(-1.0, 0.0),
        });
        let e1_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(0.0, 1.0),
        });
        let e1_on_f0_rev = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 2.0 * PI),
            direction: glam::DVec2::new(0.0, -1.0),
        });

        let geom = GeomStore {
            curves: vec![major_circle, minor_circle],
            surfaces: vec![torus_surf],
            curve2ds: vec![e0_on_f0, e0_on_f0_rev, e1_on_f0, e1_on_f0_rev],
            edge_curve: vec![Some(0), Some(1)],
            face_surface: vec![Some(0)],
            edge_pcurves: vec![
                // E0: major seam — two pcurves (forward and reverse passes)
                vec![
                    PCurve {
                        surface_idx: 0,
                        curve2d_idx: 0,
                    },
                    PCurve {
                        surface_idx: 0,
                        curve2d_idx: 1,
                    },
                ],
                // E1: minor seam — two pcurves
                vec![
                    PCurve {
                        surface_idx: 0,
                        curve2d_idx: 2,
                    },
                    PCurve {
                        surface_idx: 0,
                        curve2d_idx: 3,
                    },
                ],
            ],
            // Both seams are full circles [0, 2π]
            edge_curve_range: vec![
                Some([0.0, 2.0 * PI]), // E0 major seam circle
                Some([0.0, 2.0 * PI]), // E1 minor seam circle
            ],
            edge_degenerated: vec![false, false],
            vertex_tolerance: Vec::new(),
            edge_tolerance: Vec::new(),
            face_tolerance: Vec::new(),
            curve2d_range: Vec::new(),
            face_surface_range: Vec::new(),
            edge_same_parameter: Vec::new(),
            edge_same_range: Vec::new(),
        };

        Self {
            vertices,
            edges,
            solids: vec![solid],
            geom,
            compound: None,
            compsolid: None,
        }
    }

    /// Materializes a primitive solid descriptor into an analytic B-Rep.
    ///
    /// The resulting BRep has fully populated GeomStore entries (edge_curve,
    /// face_surface, edge_pcurves). Triangles are NOT pre-populated; the render
    /// layer tessellates on demand.
    ///
    /// User-facing code should prefer `rcad-modeling` construction helpers.
    pub fn from_primitive(primitive: PrimitiveSolid) -> Self {
        match primitive {
            PrimitiveSolid::Box {
                width,
                height,
                depth,
            } => Self::create_box(width, height, depth),
            PrimitiveSolid::Sphere { radius } => Self::create_sphere(radius),
            PrimitiveSolid::Cylinder { radius, height } => Self::create_cylinder(radius, height),
            PrimitiveSolid::Cone {
                base_radius,
                height,
            } => Self::create_cone(base_radius, height),
            PrimitiveSolid::Torus {
                major_radius,
                minor_radius,
            } => Self::create_torus(major_radius, minor_radius),
        }
    }

    pub fn center(&self) -> DVec3 {
        if self.vertices.is_empty() {
            return DVec3::ZERO;
        }
        let mut sum = DVec3::ZERO;
        for v in &self.vertices {
            sum += v.point;
        }
        sum / self.vertices.len() as f64
    }

    /// Returns the axis-aligned bounding box of all vertices as `[min, max]`,
    /// or `None` if the BRep has no vertices.
    pub fn bounding_box(&self) -> Option<[DVec3; 2]> {
        if self.vertices.is_empty() {
            return None;
        }
        let mut mn = DVec3::splat(f64::INFINITY);
        let mut mx = DVec3::splat(f64::NEG_INFINITY);
        for v in &self.vertices {
            mn = mn.min(v.point);
            mx = mx.max(v.point);
        }
        Some([mn, mx])
    }

    /// Apply a rigid-body (or general affine) transform to this BRep **in-place**.
    ///
    /// Transforms:
    /// - All vertex positions via `mat.transform_point3(p)`.
    /// - All analytic curve origins/directions and surface origins/axes.
    /// - B-spline / Bezier control points.
    /// - Recursive offset curve bases.
    ///
    /// # Notes
    /// - Direction vectors (normals, axes) are transformed with the matrix's
    ///   linear part and then renormalized — this handles rotations correctly.
    /// - Radii and offset distances are **not** scaled; callers must ensure the
    ///   transform is isometric (rotation + translation) for analytic shapes to
    ///   remain exact. Uniform scaling is supported for BSpline/Bezier (control
    ///   points are scale-transformed directly).
    /// - `glam::DAffine3` covers any combination of rotation, translation, and
    ///   uniform/non-uniform scale.
    pub fn apply_transform(&mut self, mat: glam::DAffine3) {
        use geom::{Curve3, Surface3};

        // ── Vertices ─────────────────────────────────────────────────────
        for v in &mut self.vertices {
            v.point = mat.transform_point3(v.point);
        }

        // ── Analytic curves ───────────────────────────────────────────────
        fn xf_curve(c: &mut Curve3, mat: glam::DAffine3) {
            match c {
                Curve3::Line(l) => {
                    l.origin = mat.transform_point3(l.origin);
                    l.direction = mat.transform_vector3(l.direction).normalize_or_zero();
                }
                Curve3::Circle(c3) => {
                    c3.center = mat.transform_point3(c3.center);
                    c3.normal = mat.transform_vector3(c3.normal).normalize_or_zero();
                }
                Curve3::Ellipse(e) => {
                    e.center = mat.transform_point3(e.center);
                    e.normal = mat.transform_vector3(e.normal).normalize_or_zero();
                    e.major_dir = mat.transform_vector3(e.major_dir).normalize_or_zero();
                }
                Curve3::BSpline(b) => {
                    for p in &mut b.control_points {
                        *p = mat.transform_point3(*p);
                    }
                }
                Curve3::Bezier(b) => {
                    for p in &mut b.control_points {
                        *p = mat.transform_point3(*p);
                    }
                }
                Curve3::Offset(o) => {
                    xf_curve(&mut o.basis, mat);
                    o.offset_dir = mat.transform_vector3(o.offset_dir).normalize_or_zero();
                }
                Curve3::Hyperbola(h) => {
                    h.center = mat.transform_point3(h.center);
                    h.normal = mat.transform_vector3(h.normal).normalize_or_zero();
                    h.major_dir = mat.transform_vector3(h.major_dir).normalize_or_zero();
                }
                Curve3::Parabola(p) => {
                    p.vertex = mat.transform_point3(p.vertex);
                    p.normal = mat.transform_vector3(p.normal).normalize_or_zero();
                    p.axis_dir = mat.transform_vector3(p.axis_dir).normalize_or_zero();
                }
                Curve3::CircularHelix(h) => {
                    h.origin = mat.transform_point3(h.origin);
                    h.axis = mat.transform_vector3(h.axis).normalize_or_zero();
                    h.ref_dir = mat.transform_vector3(h.ref_dir).normalize_or_zero();
                }
                Curve3::SineWave(s) => {
                    s.origin = mat.transform_point3(s.origin);
                    s.baseline_dir = mat.transform_vector3(s.baseline_dir).normalize_or_zero();
                    s.amplitude_dir = mat.transform_vector3(s.amplitude_dir).normalize_or_zero();
                }
            }
        }
        for c in &mut self.geom.curves {
            xf_curve(c, mat);
        }

        // ── Analytic surfaces ─────────────────────────────────────────────
        fn xf_surface(s: &mut Surface3, mat: glam::DAffine3) {
            match s {
                Surface3::Plane(p) => {
                    p.origin = mat.transform_point3(p.origin);
                    p.normal = mat.transform_vector3(p.normal).normalize_or_zero();
                }
                Surface3::Cylinder(c) => {
                    c.origin = mat.transform_point3(c.origin);
                    c.axis = mat.transform_vector3(c.axis).normalize_or_zero();
                }
                Surface3::Sphere(s) => {
                    s.center = mat.transform_point3(s.center);
                    s.axis = mat.transform_vector3(s.axis).normalize_or_zero();
                }
                Surface3::Cone(c) => {
                    c.apex = mat.transform_point3(c.apex);
                    c.axis = mat.transform_vector3(c.axis).normalize_or_zero();
                }
                Surface3::Torus(t) => {
                    t.center = mat.transform_point3(t.center);
                    t.axis = mat.transform_vector3(t.axis).normalize_or_zero();
                }
                Surface3::Ellipsoid(e) => {
                    e.center = mat.transform_point3(e.center);
                    e.axis = mat.transform_vector3(e.axis).normalize_or_zero();
                    e.ref_dir = mat.transform_vector3(e.ref_dir).normalize_or_zero();
                }
                Surface3::Helicoid(h) => {
                    h.origin = mat.transform_point3(h.origin);
                    h.axis = mat.transform_vector3(h.axis).normalize_or_zero();
                    h.ref_dir = mat.transform_vector3(h.ref_dir).normalize_or_zero();
                }
                Surface3::Pipe(p) => {
                    xf_surface_curve(&mut p.spine, mat);
                    p.ref_dir = mat.transform_vector3(p.ref_dir).normalize_or_zero();
                }
                Surface3::BSpline(b) => {
                    for row in &mut b.control_points {
                        for p in row.iter_mut() {
                            *p = mat.transform_point3(*p);
                        }
                    }
                }
                Surface3::Bezier(b) => {
                    for row in &mut b.control_points {
                        for p in row.iter_mut() {
                            *p = mat.transform_point3(*p);
                        }
                    }
                }
                Surface3::TriBezier(b) => {
                    for row in &mut b.control_points {
                        for p in row.iter_mut() {
                            *p = mat.transform_point3(*p);
                        }
                    }
                }
                Surface3::LinearExtrusion(le) => {
                    le.direction = mat.transform_vector3(le.direction).normalize_or_zero();
                    xf_surface_curve(&mut le.profile, mat);
                }
                Surface3::Revolution(r) => {
                    r.axis_origin = mat.transform_point3(r.axis_origin);
                    r.axis_dir = mat.transform_vector3(r.axis_dir).normalize_or_zero();
                    xf_surface_curve(&mut r.profile, mat);
                }
                Surface3::Ruled(r) => {
                    xf_surface_curve(&mut r.start, mat);
                    xf_surface_curve(&mut r.end, mat);
                }
                Surface3::Coons(c) => {
                    xf_surface_curve(&mut c.south, mat);
                    xf_surface_curve(&mut c.north, mat);
                    xf_surface_curve(&mut c.west, mat);
                    xf_surface_curve(&mut c.east, mat);
                }
                Surface3::Offset(o) => {
                    xf_surface(&mut o.basis, mat);
                }
                Surface3::Trimmed(t) => {
                    xf_surface(&mut t.basis, mat);
                    // trim domain is in parameter space — unchanged by transform
                }
            }
        }
        fn xf_surface_curve(c: &mut Box<geom::Curve3>, mat: glam::DAffine3) {
            // reuse the standalone curve transformer
            match c.as_mut() {
                geom::Curve3::Line(l) => {
                    l.origin = mat.transform_point3(l.origin);
                    l.direction = mat.transform_vector3(l.direction).normalize_or_zero();
                }
                geom::Curve3::Circle(c3) => {
                    c3.center = mat.transform_point3(c3.center);
                    c3.normal = mat.transform_vector3(c3.normal).normalize_or_zero();
                }
                geom::Curve3::Ellipse(e) => {
                    e.center = mat.transform_point3(e.center);
                    e.normal = mat.transform_vector3(e.normal).normalize_or_zero();
                    e.major_dir = mat.transform_vector3(e.major_dir).normalize_or_zero();
                }
                geom::Curve3::BSpline(b) => {
                    for p in &mut b.control_points {
                        *p = mat.transform_point3(*p);
                    }
                }
                geom::Curve3::Bezier(b) => {
                    for p in &mut b.control_points {
                        *p = mat.transform_point3(*p);
                    }
                }
                geom::Curve3::Offset(o) => {
                    o.offset_dir = mat.transform_vector3(o.offset_dir).normalize_or_zero();
                }
                geom::Curve3::Hyperbola(h) => {
                    h.center = mat.transform_point3(h.center);
                    h.normal = mat.transform_vector3(h.normal).normalize_or_zero();
                    h.major_dir = mat.transform_vector3(h.major_dir).normalize_or_zero();
                }
                geom::Curve3::Parabola(p) => {
                    p.vertex = mat.transform_point3(p.vertex);
                    p.normal = mat.transform_vector3(p.normal).normalize_or_zero();
                    p.axis_dir = mat.transform_vector3(p.axis_dir).normalize_or_zero();
                }
                geom::Curve3::CircularHelix(h) => {
                    h.origin = mat.transform_point3(h.origin);
                    h.axis = mat.transform_vector3(h.axis).normalize_or_zero();
                    h.ref_dir = mat.transform_vector3(h.ref_dir).normalize_or_zero();
                }
                geom::Curve3::SineWave(s) => {
                    s.origin = mat.transform_point3(s.origin);
                    s.baseline_dir = mat.transform_vector3(s.baseline_dir).normalize_or_zero();
                    s.amplitude_dir = mat.transform_vector3(s.amplitude_dir).normalize_or_zero();
                }
            }
        }
        for s in &mut self.geom.surfaces {
            xf_surface(s, mat);
        }

        // ── Face normals ──────────────────────────────────────────────────
        for solid in &mut self.solids {
            for shell in &mut solid.shells {
                for face in &mut shell.faces {
                    face.normal = mat.transform_vector3(face.normal).normalize_or_zero();
                }
            }
        }
    }

    /// Return a new BRep with all geometry transformed by `mat`.
    ///
    /// The original BRep is unchanged.
    pub fn transformed(&self, mat: glam::DAffine3) -> BRep {
        let mut result = self.clone();
        result.apply_transform(mat);
        result
    }

    /// Mark every face's cached triangulation as stale.
    ///
    /// After calling this, the next [`mesh_brep`] invocation will
    /// re-tessellate all faces unconditionally.
    ///
    /// Call this whenever geometry (vertices, edges, surfaces) has changed
    /// but the topology remains valid — for example after a deformation or
    /// parameter update.
    ///
    /// [`mesh_brep`]: https://docs.rs/rcad-algorithms/latest/rcad_algorithms/fn.mesh_brep.html
    pub fn invalidate_mesh(&mut self) {
        for solid in &mut self.solids {
            for shell in &mut solid.shells {
                for face in &mut shell.faces {
                    face.mesh_dirty = true;
                }
            }
        }
    }

    /// Returns `true` if any face has a stale (dirty) cached triangulation.
    ///
    /// Use this to decide whether to call [`mesh_brep`] before rendering.
    pub fn needs_remesh(&self) -> bool {
        self.solids
            .iter()
            .flat_map(|s| s.shells.iter())
            .flat_map(|sh| sh.faces.iter())
            .any(|f| f.mesh_dirty)
    }

    /// Build a deterministic baseline persistent naming table for this BRep.
    ///
    /// Labels are generated as `v{idx}`, `e{idx}`, `f{idx}`, `s{idx}`.
    pub fn persistent_naming_hooks(&self) -> PersistentNamingHooks {
        PersistentNamingHooks::with_default_labels_for_brep(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn geom_populated(brep: &BRep) -> bool {
        !brep.geom.face_surface.is_empty()
            && brep.geom.face_surface.iter().all(|s| s.is_some())
            && !brep.geom.edge_pcurves.is_empty()
    }

    #[test]
    fn creates_sphere_with_analytic_geom() {
        let brep = BRep::create_sphere(1.0);
        assert!(!brep.vertices.is_empty());
        assert!(geom_populated(&brep));
        assert!(matches!(
            brep.geom.surfaces.first(),
            Some(Surface3::Sphere(_))
        ));
        // triangles are empty — render layer tessellates on demand
        assert!(
            brep.solids
                .iter()
                .flat_map(|s| &s.shells)
                .flat_map(|sh| &sh.faces)
                .all(|f| f.triangles.is_empty())
        );
    }

    #[test]
    fn creates_cylinder_with_analytic_geom() {
        let brep = BRep::create_cylinder(1.0, 2.0);
        assert!(!brep.vertices.is_empty());
        assert!(geom_populated(&brep));
        assert!(
            brep.geom
                .surfaces
                .iter()
                .any(|s| matches!(s, Surface3::Cylinder(_)))
        );
    }

    #[test]
    fn creates_cone_with_analytic_geom() {
        let brep = BRep::create_cone(1.0, 2.0);
        assert!(!brep.vertices.is_empty());
        assert!(geom_populated(&brep));
        assert!(
            brep.geom
                .surfaces
                .iter()
                .any(|s| matches!(s, Surface3::Cone(_)))
        );
    }

    #[test]
    fn creates_torus_with_analytic_geom() {
        let brep = BRep::create_torus(1.0, 0.3);
        assert!(!brep.vertices.is_empty());
        assert!(geom_populated(&brep));
        assert!(
            brep.geom
                .surfaces
                .iter()
                .any(|s| matches!(s, Surface3::Torus(_)))
        );
    }
}
