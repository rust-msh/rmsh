//! Generic CAD kernel abstraction + rcad implementation.
//!
//! All rcad-specific types live here.  The rest of `rmsh-py` goes through
//! [`CadKernel`] and never imports rcad directly.

use std::path::Path;

use glam::{DAffine3, DVec2, DVec3};
use rmsh_model::Mesh;

// ---------------------------------------------------------------------------
// Option types �?simple enums / structs with no rcad dependency
// ---------------------------------------------------------------------------

/// STEP export protocol.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepProtocol {
    Ap214,
    Ap242,
}

impl StepProtocol {
    fn to_rcad(self) -> rcad_step::StepProtocol {
        match self {
            Self::Ap214 => rcad_step::StepProtocol::Ap214,
            Self::Ap242 => rcad_step::StepProtocol::Ap242,
        }
    }
}

/// Options for STEP export.
#[derive(Debug, Clone)]
pub struct StepExportOptions {
    pub protocol: StepProtocol,
    pub solid_color: Option<(u8, u8, u8)>,
    pub gmsh_strict: bool,
}

/// A Boolean operation type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOp {
    Cut,
    Fuse,
    Intersect,
}

// ---------------------------------------------------------------------------
// Boundary-representation types returned by the kernel.
//
// These mirror rcad's topology types but use simple Rust tuples/lists so
// that callers do not depend on rcad.
// ---------------------------------------------------------------------------

/// A 3-D curve kind returned by [`CadKernel::inspect_shape`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CurveKind {
    Line,
    Circle,
    Ellipse,
    BSpline,
    Other(String),
}

/// A surface kind returned by [`CadKernel::inspect_shape`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SurfaceKind {
    Plane,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    BSpline,
    Other(String),
}

/// Summary geometry statistics of a CAD shape.
#[derive(Debug, Clone)]
pub struct ShapeStats {
    pub vertices: usize,
    pub edges: usize,
    pub faces: usize,
    pub solids: usize,
    pub closed_edges: usize,
    pub curves_3d: usize,
    pub surfaces_3d: usize,
    pub curves_2d: usize,
    pub edge_curve_some: usize,
    pub edge_curve_none: usize,
    pub edge_pcurve_slots: usize,
    pub edges_with_pcurves: usize,
    pub total_pcurves: usize,
    pub outer_wire_edges: usize,
    pub outer_wire_edge_refs_total: usize,
    pub outer_wire_edges_per_face_min: usize,
    pub outer_wire_edges_per_face_max: usize,
    pub outer_wire_face_max_surface_kind: SurfaceKind,
    pub faces_with_outer_wire_over_100: usize,
    pub outer_wire_edges_with_curve: usize,
    pub outer_wire_edges_without_curve: usize,
    pub outer_wire_edges_with_pcurves: usize,
    pub curve3_kinds: Vec<(CurveKind, usize)>,
    pub surface3_kinds: Vec<(SurfaceKind, usize)>,
    pub pcurve_surface_kinds: Vec<(SurfaceKind, usize)>,
    pub outer_wire_curve3_kinds: Vec<(CurveKind, usize)>,
}

// ---------------------------------------------------------------------------
// The trait
// ---------------------------------------------------------------------------

/// Abstract CAD kernel �?encapsulates geometry creation, Boolean operations,
/// STEP import/export, tessellation, and shape queries.
///
/// The associated [`Shape`](Self::Shape) type is the kernel's native
/// boundary-representation.  Callers store shapes by tag and pass them back
/// into the kernel for manipulation.
pub trait CadKernel: Send + Sync {
    type Shape: Clone + Send + 'static;

    // -- Primitives ---------------------------------------------------------

    fn create_box(
        &self,
        origin: DVec3,
        dx: f64,
        dy: f64,
        dz: f64,
    ) -> Result<Self::Shape, String>;

    fn create_sphere(&self, center: DVec3, radius: f64) -> Result<Self::Shape, String>;

    fn create_cylinder(
        &self,
        origin: DVec3,
        axis: DVec3,
        ref_dir: DVec3,
        radius: f64,
        height: f64,
    ) -> Result<Self::Shape, String>;

    fn create_cone(
        &self,
        origin: DVec3,
        axis: DVec3,
        ref_dir: DVec3,
        r1: f64,
        r2: f64,
        height: f64,
    ) -> Result<Self::Shape, String>;

    fn create_torus(
        &self,
        origin: DVec3,
        axis: DVec3,
        ref_dir: DVec3,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Self::Shape, String>;

    // -- Low-level BRep building (gmsh-style API) --------------------------

    /// Create a BRep with a single vertex at `point`.
    fn make_point_shape(&self, point: DVec3) -> Self::Shape;

    /// Create a BRep with a single line edge from `p0` to `p1`.
    fn make_line_shape(&self, p0: DVec3, p1: DVec3) -> Result<Self::Shape, String>;

    /// Create a BRep with a single circular-arc edge from `p0` via `center` to `p1`.
    fn make_circle_arc_shape(
        &self,
        p0: DVec3,
        center: DVec3,
        p1: DVec3,
    ) -> Result<Self::Shape, String>;

    /// Create a BRep with a single BSpline edge interpolating `points`.
    fn make_spline_shape(&self, points: &[DVec3]) -> Result<Self::Shape, String>;

    /// Create a rectangular planar face.
    fn make_rectangle_shape(
        &self,
        x: f64,
        y: f64,
        z: f64,
        dx: f64,
        dy: f64,
    ) -> Result<Self::Shape, String>;

    /// Create a disk / elliptical planar face.
    fn make_disk_shape(
        &self,
        center: DVec3,
        rx: f64,
        ry: f64,
    ) -> Result<Self::Shape, String>;

    /// Build a planar face from one outer and zero-or-more inner curve loops.
    ///
    /// `outer_curves` and `inner_curves` are tags whose shapes are wire edges
    /// previously created via [`make_line_shape`] / [`make_circle_arc_shape`] /
    /// [`make_spline_shape`].
    fn make_plane_surface_from_curves(
        &self,
        outer_curves: &[Self::Shape],
        inner_curves: &[Vec<Self::Shape>],
    ) -> Result<Self::Shape, String>;

    // -- Boolean operations -------------------------------------------------

    fn boolean_op(
        &self,
        op: BooleanOp,
        a: &Self::Shape,
        b: &Self::Shape,
    ) -> Result<Self::Shape, String>;

    /// Fragment (split) `objects` with `tools`; returns (object_parts, tool_parts).
    fn fragment(
        &self,
        objects: &[Self::Shape],
        tools: &[Self::Shape],
    ) -> Result<(Vec<Self::Shape>, Vec<Self::Shape>), String>;

    // -- Extrude / Revolve --------------------------------------------------

    fn extrude_face(
        &self,
        shape: &Self::Shape,
        face_idx: usize,
        direction: DVec3,
        distance: f64,
    ) -> Result<Self::Shape, String>;

    fn revolve_face(
        &self,
        shape: &Self::Shape,
        face_idx: usize,
        axis_origin: DVec3,
        axis_dir: DVec3,
        angle: f64,
    ) -> Result<Self::Shape, String>;

    // -- Fillet / Chamfer ---------------------------------------------------

    fn fillet_edges(
        &self,
        shape: &Self::Shape,
        edges: &[(usize, f64)],
    ) -> Result<Self::Shape, String>;

    fn chamfer_edges(
        &self,
        shape: &Self::Shape,
        edges: &[(usize, f64)],
    ) -> Result<Self::Shape, String>;

    // -- Heal / Repair ------------------------------------------------------

    /// Repair a shape and return (repaired_shape, report_string).
    fn heal_shape(
        &self,
        shape: &Self::Shape,
        tolerance: f64,
    ) -> Result<(Self::Shape, String), String>;

    // -- STEP I/O -----------------------------------------------------------

    fn read_step_file(&self, path: &Path) -> Result<Self::Shape, String>;

    fn write_step_string(
        &self,
        shape: &Self::Shape,
        options: &StepExportOptions,
    ) -> Result<String, String>;

    // -- Tessellation -------------------------------------------------------

    /// Convert a shape into a tri-mesh (`Mesh`).
    fn tessellate(&self, shape: &Self::Shape) -> Mesh;

    // -- Transforms ---------------------------------------------------------

    fn apply_transform(&self, shape: &mut Self::Shape, xf: DAffine3);

    // -- Properties ---------------------------------------------------------

    fn volume(&self, shape: &Self::Shape) -> f64;
    fn surface_area(&self, shape: &Self::Shape) -> f64;
    fn centroid(&self, shape: &Self::Shape) -> DVec3;

    // -- Queries ------------------------------------------------------------

    /// Topological dimension (0=point, 1=curve, 2=surface, 3=solid).
    fn dimension(&self, shape: &Self::Shape) -> i32;

    fn bounding_box(&self, shape: &Self::Shape) -> Option<([f64; 3], [f64; 3])>;

    fn vertex_count(&self, shape: &Self::Shape) -> usize;
    fn edge_count(&self, shape: &Self::Shape) -> usize;
    fn has_spherical_surface(&self, shape: &Self::Shape) -> bool;

    /// Detailed shape-stats for debugging / introspection.
    fn inspect_shape(&self, shape: &Self::Shape) -> ShapeStats;

    /// Merge multiple shapes into one for export (re-indexes all entities).
    fn merge_for_export(&self, shapes: &[Self::Shape]) -> Self::Shape;

    /// Explode a compound shape into individual solid components.
    fn explode_solids(&self, shape: &Self::Shape) -> Vec<Self::Shape>;
}

// ---------------------------------------------------------------------------
// rcad implementation
// ---------------------------------------------------------------------------

/// The default backend backed by the `rcad` kernel (OCCT-based).
#[derive(Default)]
pub struct RcadKernel;

impl CadKernel for RcadKernel {
    type Shape = rcad_kernel::BRep;

    // -- Primitives ---------------------------------------------------------

    fn create_box(
        &self,
        origin: DVec3,
        dx: f64,
        dy: f64,
        dz: f64,
    ) -> Result<Self::Shape, String> {
        rcad_modeling::box_brep(origin, DVec3::X, DVec3::Y, dx, dy, dz)
            .map_err(|e| e.to_string())
    }

    fn create_sphere(&self, center: DVec3, radius: f64) -> Result<Self::Shape, String> {
        rcad_modeling::sphere_brep(center, radius)
            .map_err(|e| e.to_string())
    }

    fn create_cylinder(
        &self,
        origin: DVec3,
        axis: DVec3,
        ref_dir: DVec3,
        radius: f64,
        height: f64,
    ) -> Result<Self::Shape, String> {
        rcad_modeling::cylinder_brep(origin, axis, ref_dir, radius, height)
            .map_err(|e| e.to_string())
    }

    fn create_cone(
        &self,
        origin: DVec3,
        axis: DVec3,
        ref_dir: DVec3,
        r1: f64,
        r2: f64,
        height: f64,
    ) -> Result<Self::Shape, String> {
        if (r1 - r2).abs() <= 1e-12 {
            rcad_modeling::cylinder_brep(origin, axis, ref_dir, r1, height)
                .map_err(|e| e.to_string())
        } else if r2 <= 1e-12 {
            rcad_modeling::cone_brep(origin, axis, ref_dir, r1, height)
                .map_err(|e| e.to_string())
        } else if r1 <= 1e-12 {
            rcad_modeling::cone_brep(origin, -axis, ref_dir, r2, height)
                .map_err(|e| e.to_string())
        } else {
            Ok(frustum_brep_impl(origin, axis, ref_dir, r1, r2, height))
        }
    }

    fn create_torus(
        &self,
        origin: DVec3,
        axis: DVec3,
        ref_dir: DVec3,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Self::Shape, String> {
        rcad_modeling::torus_brep(origin, axis, ref_dir, major_radius, minor_radius)
            .map_err(|e| e.to_string())
    }

    // -- Low-level BRep building -------------------------------------------

    fn make_point_shape(&self, point: DVec3) -> Self::Shape {
        let mut brep = rcad_kernel::BRep::new();
        make_vertex_internal(&mut brep, point);
        brep
    }

    fn make_line_shape(&self, p0: DVec3, p1: DVec3) -> Result<Self::Shape, String> {
        let seg = p1 - p0;
        let len = seg.length();
        if len < 1e-12 {
            return Err("degenerate line: start and end points are coincident".into());
        }
        let dir = seg / len;
        let mut brep = rcad_kernel::BRep::new();
        let v0 = make_vertex_internal(&mut brep, p0);
        let v1 = make_vertex_internal(&mut brep, p1);
        rcad_modeling::builder::make_edge(
            &mut brep,
            rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: p0,
                direction: dir,
            }),
            0.0,
            len,
            v0,
            v1,
        )
        .map_err(|e| e.to_string())?;
        Ok(brep)
    }

    fn make_circle_arc_shape(
        &self,
        p0: DVec3,
        center: DVec3,
        p1: DVec3,
    ) -> Result<Self::Shape, String> {
        let r0 = p0 - center;
        let r1 = p1 - center;
        let radius = r0.length();
        if radius < 1e-12 {
            return Err("circle arc radius is zero".into());
        }
        if (r1.length() - radius).abs() > 1e-6 * radius.max(1.0) {
            return Err("start/end points are not equidistant from center".into());
        }

        let normal = r0.cross(r1);
        if normal.length_squared() < 1e-20 {
            return Err("circle arc points are collinear".into());
        }
        let normal = normal.normalize();
        let x_axis = rcad_kernel::geom::any_perpendicular(normal);
        let y_axis = normal.cross(x_axis);

        use std::f64::consts::PI;
        let t0 = r0.dot(y_axis).atan2(r0.dot(x_axis));
        let mut t1 = r1.dot(y_axis).atan2(r1.dot(x_axis));
        while t1 <= t0 {
            t1 += 2.0 * PI;
        }

        let mut brep = rcad_kernel::BRep::new();
        let v0 = make_vertex_internal(&mut brep, p0);
        let v1 = make_vertex_internal(&mut brep, p1);
        rcad_modeling::builder::make_edge(
            &mut brep,
            rcad_kernel::geom::Curve3::Circle(rcad_kernel::geom::Circle3::new(center, normal, radius)),
            t0,
            t1,
            v0,
            v1,
        )
        .map_err(|e| e.to_string())?;
        Ok(brep)
    }

    fn make_spline_shape(&self, points: &[DVec3]) -> Result<Self::Shape, String> {
        if points.len() < 2 {
            return Err("spline requires at least 2 points".into());
        }
        let bspline = rcad_kernel::fit::interpolate_points(points)
            .map_err(|e| format!("spline interpolation failed: {e}"))?;
        let t0 = bspline
            .knots
            .get(bspline.degree)
            .copied()
            .ok_or("spline: invalid knot vector")?;
        let t1 = bspline
            .knots
            .get(bspline.knots.len().saturating_sub(bspline.degree + 1))
            .copied()
            .ok_or("spline: invalid knot vector")?;

        let mut brep = rcad_kernel::BRep::new();
        let v0 = make_vertex_internal(&mut brep, points[0]);
        let v1 = make_vertex_internal(&mut brep, *points.last().expect("len >= 2"));
        rcad_modeling::builder::make_edge(
            &mut brep,
            rcad_kernel::geom::Curve3::BSpline(bspline),
            t0,
            t1,
            v0,
            v1,
        )
        .map_err(|e| e.to_string())?;
        Ok(brep)
    }

    fn make_rectangle_shape(
        &self,
        x: f64,
        y: f64,
        z: f64,
        dx: f64,
        dy: f64,
    ) -> Result<Self::Shape, String> {
        if dx.abs() < 1e-12 || dy.abs() < 1e-12 {
            return Err("dx and dy must be non-zero".into());
        }

        let p0 = DVec3::new(x, y, z);
        let p1 = DVec3::new(x + dx, y, z);
        let p2 = DVec3::new(x + dx, y + dy, z);
        let p3 = DVec3::new(x, y + dy, z);

        let mut shape = rcad_kernel::BRep::new();

        let v0 = make_vertex_internal(&mut shape, p0);
        let v1 = make_vertex_internal(&mut shape, p1);
        let v2 = make_vertex_internal(&mut shape, p2);
        let v3 = make_vertex_internal(&mut shape, p3);

        let mk_edge = |brep: &mut rcad_kernel::BRep, a: DVec3, b: DVec3, va: usize, vb: usize| {
            let seg = b - a;
            let len = seg.length();
            let dir = seg / len;
            rcad_modeling::builder::make_edge(
                brep,
                rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: a,
                    direction: dir,
                }),
                0.0,
                len,
                va,
                vb,
            )
        };

        let e0 = mk_edge(&mut shape, p0, p1, v0, v1).map_err(|e| e.to_string())?;
        let e1 = mk_edge(&mut shape, p1, p2, v1, v2).map_err(|e| e.to_string())?;
        let e2 = mk_edge(&mut shape, p2, p3, v2, v3).map_err(|e| e.to_string())?;
        let e3 = mk_edge(&mut shape, p3, p0, v3, v0).map_err(|e| e.to_string())?;

        let outer = rcad_kernel::topology::Wire {
            edges: vec![
                rcad_kernel::topology::WireEdge::fwd(e0),
                rcad_kernel::topology::WireEdge::fwd(e1),
                rcad_kernel::topology::WireEdge::fwd(e2),
                rcad_kernel::topology::WireEdge::fwd(e3),
            ],
        };

        let normal = (p1 - p0).cross(p3 - p0).normalize_or_zero();
        rcad_modeling::builder::make_face(
            &mut shape,
            rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane {
                origin: p0,
                normal,
            }),
            outer,
            Vec::new(),
        )
        .map_err(|e| e.to_string())?;

        Ok(shape)
    }

    fn make_disk_shape(
        &self,
        center: DVec3,
        rx: f64,
        ry: f64,
    ) -> Result<Self::Shape, String> {
        use rcad_kernel::geom::*;
        use rcad_kernel::topology::*;

        let major_is_x = rx.abs() >= ry.abs();
        let major_radius = if major_is_x { rx.abs() } else { ry.abs() };
        let minor_radius = if major_is_x { ry.abs() } else { rx.abs() };
        let major_dir = if major_is_x { DVec3::X } else { DVec3::Y };
        let mut shape = rcad_kernel::BRep::new();

        let edge_curve: Curve3 = if (rx.abs() - ry.abs()).abs() < 1e-12 {
            Curve3::Circle(Circle3::new(center, DVec3::Z, rx.abs()))
        } else {
            Curve3::Ellipse(Ellipse3 {
                center,
                normal: DVec3::Z,
                major_dir,
                major_radius,
                minor_radius,
            })
        };

        let v0 = shape.vertices.len();
        let v1 = shape.vertices.len();
        // Single edge: make_edge returns edge index.
        let edge = rcad_modeling::builder::make_edge(
            &mut shape,
            edge_curve,
            0.0,
            2.0 * std::f64::consts::PI,
            v0,
            v1,
        )
        .map_err(|e| e.to_string())?;

        let wire = rcad_kernel::topology::Wire {
            edges: vec![rcad_kernel::topology::WireEdge::fwd(edge)],
        };

        rcad_modeling::builder::make_face(
            &mut shape,
            rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane {
                origin: center,
                normal: DVec3::Z,
            }),
            wire,
            Vec::new(),
        )
        .map_err(|e| e.to_string())?;

        // Add pcurve for the edge.
        let curve2d_idx = if (rx.abs() - ry.abs()).abs() < 1e-12 {
            shape.geom.curve2ds.push(Curve2d::Circle(Circle2d::new(glam::DVec2::new(center.x, center.y), rx.abs())));
            shape.geom.curve2ds.len() - 1
        } else {
            shape.geom.curve2ds.push(Curve2d::Ellipse(Ellipse2d {
                center: glam::DVec2::new(center.x, center.y),
                major_dir: if major_is_x {
                    glam::DVec2::X
                } else {
                    glam::DVec2::Y
                },
                major_radius,
                minor_radius,
            }));
            shape.geom.curve2ds.len() - 1
        };
        while shape.geom.edge_pcurves.len() <= edge {
            shape.geom.edge_pcurves.push(Vec::new());
        }
        shape.geom.edge_pcurves[edge].push(rcad_kernel::PCurve {
            surface_idx: 0,
            curve2d_idx,
        });

        Ok(shape)
    }

    fn make_plane_surface_from_curves(
        &self,
        outer_curves: &[Self::Shape],
        inner_curves: &[Vec<Self::Shape>],
    ) -> Result<Self::Shape, String> {
        if outer_curves.is_empty() {
            return Err("plane surface requires at least one outer curve".into());
        }

        use rcad_kernel::geom::*;
        use rcad_kernel::topology::*;

        let mut shape = rcad_kernel::BRep::new();
        let mut vertex_map: Vec<(DVec3, usize)> = Vec::new();

        let get_or_add_vertex = |brep: &mut rcad_kernel::BRep,
                                 map: &mut Vec<(DVec3, usize)>,
                                 p: DVec3|
         -> usize {
            if let Some((_, idx)) = map.iter().find(|(q, _)| (*q - p).length_squared() <= 1e-18) {
                *idx
            } else {
                let idx = make_vertex_internal(brep, p);
                map.push((p, idx));
                idx
            }
        };

        let copy_edge =
            |brep: &mut rcad_kernel::BRep,
             src: &rcad_kernel::BRep,
             src_edge_idx: usize,
             reverse: bool,
             map: &mut Vec<(DVec3, usize)>,
             get_or_add: &dyn Fn(&mut rcad_kernel::BRep, &mut Vec<(DVec3, usize)>, DVec3) -> usize|
             -> Result<usize, String> {
                let curve_idx = src
                    .geom
                    .edge_curve
                    .get(src_edge_idx)
                    .copied()
                    .flatten()
                    .ok_or("edge has no 3D curve")?;
                let curve = src
                    .geom
                    .curves
                    .get(curve_idx)
                    .cloned()
                    .ok_or("invalid curve index")?;

                let range = src
                    .geom
                    .edge_curve_range
                    .get(src_edge_idx)
                    .copied()
                    .flatten()
                    .unwrap_or_else(|| {
                        let [t0, t1] = curve.default_domain();
                        [t0, t1]
                    });

                let mut p_start = src
                    .vertices
                    .get(src.edges[src_edge_idx].start)
                    .map(|v| v.point)
                    .ok_or("invalid edge start")?;
                let mut p_end = src
                    .vertices
                    .get(src.edges[src_edge_idx].end)
                    .map(|v| v.point)
                    .ok_or("invalid edge end")?;

                let (mut t0, mut t1) = (range[0], range[1]);
                if reverse {
                    std::mem::swap(&mut p_start, &mut p_end);
                    std::mem::swap(&mut t0, &mut t1);
                }

                let v0 = get_or_add(brep, map, p_start);
                let v1 = get_or_add(brep, map, p_end);
                rcad_modeling::builder::make_edge(brep, curve, t0, t1, v0, v1)
                    .map_err(|e| e.to_string())
            };

        let mut loop_wires: Vec<Wire> = Vec::new();
        let mut first_loop_points: Vec<DVec3> = Vec::new();

        // Outer loop
        let mut wire_edges = Vec::with_capacity(outer_curves.len());
        for src in outer_curves {
            let e = copy_edge(
                &mut shape,
                src,
                0,
                false,
                &mut vertex_map,
                &get_or_add_vertex,
            )?;
            if first_loop_points.is_empty() {
                if let Some(v) = src.vertices.get(src.edges[0].start) {
                    first_loop_points.push(v.point);
                }
            }
            wire_edges.push(WireEdge::fwd(e));
        }
        loop_wires.push(Wire { edges: wire_edges });

        // Inner loops
        for inner in inner_curves {
            let mut wire_edges = Vec::with_capacity(inner.len());
            for src in inner {
                let e = copy_edge(
                    &mut shape,
                    src,
                    0,
                    false,
                    &mut vertex_map,
                    &get_or_add_vertex,
                )?;
                wire_edges.push(WireEdge::fwd(e));
            }
            loop_wires.push(Wire { edges: wire_edges });
        }

        let outer = loop_wires.remove(0);
        let mut normal = DVec3::Z;
        let mut origin = first_loop_points.first().copied().unwrap_or(DVec3::ZERO);
        if first_loop_points.len() >= 3 {
            origin = first_loop_points[0];
            for i in 1..first_loop_points.len().saturating_sub(1) {
                let a = first_loop_points[i] - origin;
                let b = first_loop_points[i + 1] - origin;
                let n = a.cross(b);
                if n.length_squared() > 1e-20 {
                    normal = n.normalize();
                    break;
                }
            }
        }

        rcad_modeling::builder::make_face(
            &mut shape,
            Surface3::Plane(Plane { origin, normal }),
            outer,
            loop_wires,
        )
        .map_err(|e| e.to_string())?;

        Ok(shape)
    }

    // -- Boolean operations -------------------------------------------------

    fn boolean_op(
        &self,
        op: BooleanOp,
        a: &Self::Shape,
        b: &Self::Shape,
    ) -> Result<Self::Shape, String> {
        use rcad_algorithms::{BooleanOpType, SimplifyOptions, boolean_op_simplified};
        use rcad_algorithms::healing::{ShapeProcessConfig, run_shape_process};

        let rcad_op = match op {
            BooleanOp::Cut => BooleanOpType::Difference,
            BooleanOp::Fuse => BooleanOpType::Union,
            BooleanOp::Intersect => BooleanOpType::Intersection,
        };

        let options = SimplifyOptions::default();
        boolean_op_simplified(rcad_op, a, b, options)
            .map(|(brep, _)| {
                use rcad_algorithms::brep_repair::{
                    make_connected_iterative_with_growth_cap,
                    remove_internal_faces_post_boolean,
                };

                let (cleaned, _) =
                    run_shape_process(&brep, &ShapeProcessConfig::boolean_cleanup_preset());

                if matches!(op, BooleanOp::Intersect) {
                    let (no_internal, _) = remove_internal_faces_post_boolean(&cleaned);
                    let (topo_reduced, _) = make_connected_iterative_with_growth_cap(
                        &no_internal, 1.0e-7, 4, 2.0, 1.0e-4,
                    );

                    if has_exportable_shell_topology(&topo_reduced) {
                        return topo_reduced;
                    }
                    if has_exportable_shell_topology(&cleaned) {
                        return cleaned;
                    }
                    return brep;
                }

                let (no_internal, _) = remove_internal_faces_post_boolean(&cleaned);
                let (topo_reduced, _) = make_connected_iterative_with_growth_cap(
                    &no_internal, 1.0e-7, 4, 2.0, 1.0e-4,
                );
                topo_reduced
            })
            .map_err(|e| e.to_string())
    }

    fn fragment(
        &self,
        objects: &[Self::Shape],
        tools: &[Self::Shape],
    ) -> Result<(Vec<Self::Shape>, Vec<Self::Shape>), String> {
        use rcad_algorithms::split_objects_with_tools;

        if objects.len() == 1 && tools.len() == 1 {
            let obj = &objects[0];
            let tool = &tools[0];

            let object_parts = self.explode_solids(obj);
            let tool_parts = self.explode_solids(tool);

            let split_objects = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let (splits, _) = split_objects_with_tools(&[obj.clone()], &[tool.clone()]);
                splits
                    .into_iter()
                    .next()
                    .map(|b| self.explode_solids(&b))
                    .unwrap_or_else(|| object_parts.clone())
            }))
            .unwrap_or(object_parts);

            let split_tools = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let (splits, _) = split_objects_with_tools(&[tool.clone()], &[obj.clone()]);
                splits
                    .into_iter()
                    .next()
                    .map(|b| self.explode_solids(&b))
                    .unwrap_or_else(|| tool_parts.clone())
            }))
            .unwrap_or(tool_parts);

            Ok((split_objects, split_tools))
        } else {
            let (split_objects, _) = split_objects_with_tools(objects, tools);
            let (split_tools, _) = split_objects_with_tools(tools, objects);
            Ok((split_objects, split_tools))
        }
    }

    // -- Extrude / Revolve --------------------------------------------------

    fn extrude_face(
        &self,
        shape: &Self::Shape,
        face_idx: usize,
        direction: DVec3,
        distance: f64,
    ) -> Result<Self::Shape, String> {
        rcad_modeling::builder::ops::extrude(shape, face_idx, direction, distance)
            .map_err(|e| e.to_string())
    }

    fn revolve_face(
        &self,
        shape: &Self::Shape,
        face_idx: usize,
        axis_origin: DVec3,
        axis_dir: DVec3,
        angle: f64,
    ) -> Result<Self::Shape, String> {
        rcad_modeling::builder::ops::revolve(shape, face_idx, axis_origin, axis_dir, angle)
            .map_err(|e| e.to_string())
    }

    // -- Fillet / Chamfer ---------------------------------------------------

    fn fillet_edges(
        &self,
        shape: &Self::Shape,
        edges: &[(usize, f64)],
    ) -> Result<Self::Shape, String> {
        rcad_modeling::builder::fillet::fillet_edges(shape, edges)
            .map_err(|e| e.to_string())
    }

    fn chamfer_edges(
        &self,
        shape: &Self::Shape,
        edges: &[(usize, f64)],
    ) -> Result<Self::Shape, String> {
        let mut edges_sorted: Vec<(usize, f64)> = edges.to_vec();
        edges_sorted.sort_by(|a, b| b.0.cmp(&a.0));
        let mut current = shape.clone();
        for (edge_idx, dist) in &edges_sorted {
            current = rcad_modeling::builder::fillet::chamfer_edge(&current, *edge_idx, *dist)
                .map_err(|e| e.to_string())?;
        }
        Ok(current)
    }

    // -- Heal / Repair ------------------------------------------------------

    fn heal_shape(
        &self,
        shape: &Self::Shape,
        tolerance: f64,
    ) -> Result<(Self::Shape, String), String> {
        let (repaired, report) = rcad_algorithms::brep_repair::repair(shape, tolerance);
        let report_str = format!(
            "vertices_merged={} degenerate_faces_removed={} normals_recomputed={} wires_fixed={}",
            report.vertices_merged,
            report.degenerate_faces_removed,
            report.normals_recomputed,
            report.wires_fixed,
        );
        Ok((repaired, report_str))
    }

    // -- STEP I/O -----------------------------------------------------------

    fn read_step_file(&self, path: &Path) -> Result<Self::Shape, String> {
        rcad_step::StepReader::read_file(path).map_err(|e| e.to_string())
    }

    fn write_step_string(
        &self,
        shape: &Self::Shape,
        options: &StepExportOptions,
    ) -> Result<String, String> {
        use rcad_kernel::appearance::Color;
        use rcad_step::writer::{StepHeader, StepWriteOptions, StepWriter};
        use rcad_step::ExportSelection;

        let colors = options.solid_color.map(|(r, g, b)| {
            rcad_kernel::appearance::StepColor {
                solid_color: Some(Color::from_rgb8(r, g, b)),
                face_colors: Vec::new(),
            }
        });

        if options.gmsh_strict {
            let fixed = normalize_for_strict_step_export(shape);
            let step_options = StepWriteOptions {
                protocol: options.protocol.to_rcad(),
                colors,
                properties: Vec::new(),
                ap242_metadata: None,
                header: StepHeader::default(),
                export_standalone_wire_overlay: true,

            };
            Ok(StepWriter::write_string_with_options(
                &fixed,
                ExportSelection {
                    selected_faces: &[],
                    selected_edges: &[],
                },
                &step_options,
            ))
        } else {
            let step_options = StepWriteOptions {
                protocol: options.protocol.to_rcad(),
                colors,
                properties: Vec::new(),
                ap242_metadata: None,
                header: StepHeader::default(),
                export_standalone_wire_overlay: true,

            };
            Ok(StepWriter::write_string_with_options(
                shape,
                ExportSelection {
                    selected_faces: &[],
                    selected_edges: &[],
                },
                &step_options,
            ))
        }
    }

    // -- Tessellation -------------------------------------------------------

    fn tessellate(&self, shape: &Self::Shape) -> Mesh {
        let mut mesh = Mesh::new();
        let mut node_id: u64 = 1;
        let mut elem_id: u64 = 1;
        let mut vertex_to_node: std::collections::HashMap<usize, u64> = std::collections::HashMap::new();

        for solid in &shape.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for &[i0, i1, i2] in &face.triangles {
                        let mut nids = Vec::with_capacity(3);
                        for vi in [i0, i1, i2] {
                            let nid = *vertex_to_node.entry(vi).or_insert_with(|| {
                                let id = node_id;
                                node_id += 1;
                                if let Some(v) = shape.vertices.get(vi) {
                                    mesh.add_node(rmsh_model::Node::new(id, v.point.x, v.point.y, v.point.z));
                                }
                                id
                            });
                            nids.push(nid);
                        }
                        mesh.add_element(rmsh_model::Element::new(
                            elem_id,
                            rmsh_model::ElementType::Triangle3,
                            nids,
                        ));
                        elem_id += 1;
                    }
                }
            }
        }
        mesh
    }

    // -- Transforms ---------------------------------------------------------

    fn apply_transform(&self, shape: &mut Self::Shape, xf: DAffine3) {
        shape.apply_transform(xf);
    }

    // -- Properties ---------------------------------------------------------

    fn volume(&self, shape: &Self::Shape) -> f64 {
        rcad_kernel::properties::volume(shape).abs()
    }

    fn surface_area(&self, shape: &Self::Shape) -> f64 {
        rcad_kernel::properties::surface_area(shape)
    }

    fn centroid(&self, shape: &Self::Shape) -> DVec3 {
        rcad_kernel::properties::centroid(shape)
    }

    // -- Queries ------------------------------------------------------------

    fn dimension(&self, shape: &Self::Shape) -> i32 {
        cad_shape_dimension_impl(shape)
    }

    fn bounding_box(&self, shape: &Self::Shape) -> Option<([f64; 3], [f64; 3])> {
        if shape.vertices.is_empty() {
            return None;
        }
        let mut min = [f64::MAX; 3];
        let mut max = [f64::MIN; 3];
        for v in &shape.vertices {
            min[0] = min[0].min(v.point.x);
            min[1] = min[1].min(v.point.y);
            min[2] = min[2].min(v.point.z);
            max[0] = max[0].max(v.point.x);
            max[1] = max[1].max(v.point.y);
            max[2] = max[2].max(v.point.z);
        }
        Some((min, max))
    }

    fn vertex_count(&self, shape: &Self::Shape) -> usize {
        shape.vertices.len()
    }

    fn edge_count(&self, shape: &Self::Shape) -> usize {
        shape.edges.len()
    }

    fn has_spherical_surface(&self, shape: &Self::Shape) -> bool {
        shape
            .geom
            .surfaces
            .iter()
            .any(|s| matches!(s, rcad_kernel::geom::Surface3::Sphere(_)))
    }

    fn inspect_shape(&self, shape: &Self::Shape) -> ShapeStats {
        use std::collections::{BTreeMap, HashSet};

        let face_count: usize = shape
            .solids
            .iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();

        let mut total_outer_wire_edge_refs = 0usize;
        let mut min_outer_wire_edges_per_face = usize::MAX;
        let mut max_outer_wire_edges_per_face = 0usize;
        let mut max_outer_wire_face_surface_kind_str = "UnknownSurface";
        let mut faces_with_outer_wire_over_100 = 0usize;
        let mut face_cursor = 0usize;

        for solid in &shape.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let c = face.outer_wire.edges.len();
                    total_outer_wire_edge_refs += c;
                    min_outer_wire_edges_per_face = min_outer_wire_edges_per_face.min(c);
                    if c > max_outer_wire_edges_per_face {
                        max_outer_wire_edges_per_face = c;
                        max_outer_wire_face_surface_kind_str = shape
                            .geom
                            .face_surface
                            .get(face_cursor)
                            .copied()
                            .flatten()
                            .and_then(|si| shape.geom.surfaces.get(si))
                            .map(surface3_kind_name_impl)
                            .unwrap_or("UnknownSurface");
                    }
                    if c >= 100 {
                        faces_with_outer_wire_over_100 += 1;
                    }
                    face_cursor += 1;
                }
            }
        }
        let min_outer_wire_edges_per_face = if face_count == 0 {
            0
        } else {
            min_outer_wire_edges_per_face
        };

        let edge_curve_some = shape.geom.edge_curve.iter().filter(|c| c.is_some()).count();
        let edge_curve_none = shape.edges.len().saturating_sub(edge_curve_some);
        let edge_closed = shape.edges.iter().filter(|e| e.start == e.end).count();

        let mut outer_wire_edge_indices: HashSet<usize> = HashSet::new();
        for solid in &shape.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for we in &face.outer_wire.edges {
                        outer_wire_edge_indices.insert(we.idx);
                    }
                }
            }
        }
        let outer_wire_edge_count = outer_wire_edge_indices.len();
        let outer_wire_edges_with_curve = outer_wire_edge_indices
            .iter()
            .filter(|&&ei| shape.geom.edge_curve.get(ei).copied().flatten().is_some())
            .count();
        let outer_wire_edges_without_curve =
            outer_wire_edge_count.saturating_sub(outer_wire_edges_with_curve);
        let outer_wire_edges_with_pcurves = outer_wire_edge_indices
            .iter()
            .filter(|&&ei| {
                shape
                    .geom
                    .edge_pcurves
                    .get(ei)
                    .map(|pcs| !pcs.is_empty())
                    .unwrap_or(false)
            })
            .count();

        let edges_with_pcurves = shape
            .geom
            .edge_pcurves
            .iter()
            .filter(|pcs| !pcs.is_empty())
            .count();
        let total_pcurves: usize = shape.geom.edge_pcurves.iter().map(|pcs| pcs.len()).sum();

        let mut curve_kinds: BTreeMap<&'static str, usize> = BTreeMap::new();
        for c in &shape.geom.curves {
            *curve_kinds.entry(curve3_kind_name_impl(c)).or_insert(0) += 1;
        }

        let mut surface_kinds: BTreeMap<&'static str, usize> = BTreeMap::new();
        for s in &shape.geom.surfaces {
            *surface_kinds.entry(surface3_kind_name_impl(s)).or_insert(0) += 1;
        }

        let mut pcurve_surface_kinds: BTreeMap<&'static str, usize> = BTreeMap::new();
        for pcs in &shape.geom.edge_pcurves {
            for pc in pcs {
                let key = shape
                    .geom
                    .surfaces
                    .get(pc.surface_idx)
                    .map(surface3_kind_name_impl)
                    .unwrap_or("UnknownSurface");
                *pcurve_surface_kinds.entry(key).or_insert(0) += 1;
            }
        }

        let mut outer_wire_curve_kinds: BTreeMap<&'static str, usize> = BTreeMap::new();
        for &ei in &outer_wire_edge_indices {
            if let Some(curve_idx) = shape.geom.edge_curve.get(ei).copied().flatten()
                && let Some(curve) = shape.geom.curves.get(curve_idx)
            {
                *outer_wire_curve_kinds
                    .entry(curve3_kind_name_impl(curve))
                    .or_insert(0) += 1;
            }
        }

        let to_kind_pairs = |map: BTreeMap<&'static str, usize>| -> Vec<(CurveKind, usize)> {
            map.into_iter()
                .map(|(k, v)| (str_to_curve_kind(k), v))
                .collect()
        };
        let to_surface_pairs =
            |map: BTreeMap<&'static str, usize>| -> Vec<(SurfaceKind, usize)> {
                map.into_iter()
                    .map(|(k, v)| (str_to_surface_kind(k), v))
                    .collect()
            };

        ShapeStats {
            vertices: shape.vertices.len(),
            edges: shape.edges.len(),
            faces: face_count,
            solids: shape.solids.len(),
            closed_edges: edge_closed,
            curves_3d: shape.geom.curves.len(),
            surfaces_3d: shape.geom.surfaces.len(),
            curves_2d: shape.geom.curve2ds.len(),
            edge_curve_some,
            edge_curve_none,
            edge_pcurve_slots: shape.geom.edge_pcurves.len(),
            edges_with_pcurves,
            total_pcurves,
            outer_wire_edges: outer_wire_edge_count,
            outer_wire_edge_refs_total: total_outer_wire_edge_refs,
            outer_wire_edges_per_face_min: min_outer_wire_edges_per_face,
            outer_wire_edges_per_face_max: max_outer_wire_edges_per_face,
            outer_wire_face_max_surface_kind: str_to_surface_kind(
                max_outer_wire_face_surface_kind_str,
            ),
            faces_with_outer_wire_over_100,
            outer_wire_edges_with_curve,
            outer_wire_edges_without_curve,
            outer_wire_edges_with_pcurves,
            curve3_kinds: to_kind_pairs(curve_kinds),
            surface3_kinds: to_surface_pairs(surface_kinds),
            pcurve_surface_kinds: to_surface_pairs(pcurve_surface_kinds),
            outer_wire_curve3_kinds: to_kind_pairs(outer_wire_curve_kinds),
        }
    }

    fn merge_for_export(&self, shapes: &[Self::Shape]) -> Self::Shape {
        merge_breps_for_export_impl(shapes)
    }

    fn explode_solids(&self, shape: &Self::Shape) -> Vec<Self::Shape> {
        let connected = explode_by_connected_components_impl(shape);
        if connected.len() > 1 {
            return connected;
        }

        let mut out: Vec<rcad_kernel::BRep> = Vec::new();
        for solid in &shape.solids {
            if solid.shells.is_empty() {
                let mut part = shape.clone();
                part.solids = vec![solid.clone()];
                out.push(part);
                continue;
            }
            for shell in &solid.shells {
                let mut part = shape.clone();
                part.solids = vec![rcad_kernel::topology::Solid {
                    shells: vec![shell.clone()],
                }];
                out.push(part);
            }
        }
        if out.is_empty() {
            out.push(shape.clone());
        }
        out
    }
}

// ---------------------------------------------------------------------------
// Internal helpers (not exposed through the trait)
// ---------------------------------------------------------------------------

fn make_vertex_internal(brep: &mut rcad_kernel::BRep, point: DVec3) -> usize {
    let idx = brep.vertices.len();
    brep.vertices.push(rcad_kernel::topology::Vertex { point });
    idx
}

fn cad_shape_dimension_impl(shape: &rcad_kernel::BRep) -> i32 {
    if shape.solids.is_empty() {
        if !shape.edges.is_empty() {
            return 1;
        }
        if !shape.vertices.is_empty() {
            return 0;
        }
        return 0;
    }

    let face_count: usize = shape
        .solids
        .iter()
        .flat_map(|solid| solid.shells.iter())
        .map(|shell| shell.faces.len())
        .sum();
    if face_count == 1 {
        2
    } else {
        3
    }
}

fn curve3_kind_name_impl(curve: &rcad_kernel::geom::Curve3) -> &'static str {
    use rcad_kernel::geom::Curve3;
    match curve {
        Curve3::Line(_) => "Line",
        Curve3::Circle(_) => "Circle",
        Curve3::Ellipse(_) => "Ellipse",
        Curve3::BSpline(_) => "BSpline",
        Curve3::Hyperbola(_) => "Hyperbola",
        Curve3::Parabola(_) => "Parabola",
        Curve3::CircularHelix(_) => "CircularHelix",
        Curve3::SineWave(_) => "SineWave",
        Curve3::Offset(_) => "Offset",
        Curve3::Bezier(_) => "Bezier",
    }
}

fn surface3_kind_name_impl(surface: &rcad_kernel::geom::Surface3) -> &'static str {
    use rcad_kernel::geom::Surface3;
    match surface {
        Surface3::Plane(_) => "Plane",
        Surface3::Cylinder(_) => "Cylinder",
        Surface3::Sphere(_) => "Sphere",
        Surface3::Cone(_) => "Cone",
        Surface3::Torus(_) => "Torus",
        Surface3::BSpline(_) => "BSpline",
        Surface3::Ellipsoid(_) => "Ellipsoid",
        Surface3::Helicoid(_) => "Helicoid",
        Surface3::Pipe(_) => "Pipe",
        Surface3::LinearExtrusion(_) => "LinearExtrusion",
        Surface3::Revolution(_) => "Revolution",
        Surface3::Ruled(_) => "Ruled",
        Surface3::Coons(_) => "Coons",
        Surface3::TriBezier(_) => "TriBezier",
        Surface3::Bezier(_) => "Bezier",
        Surface3::Offset(_) => "Offset",
        Surface3::Trimmed(_) => "Trimmed",
    }
}

fn str_to_curve_kind(s: &str) -> CurveKind {
    match s {
        "Line" => CurveKind::Line,
        "Circle" => CurveKind::Circle,
        "Ellipse" => CurveKind::Ellipse,
        "BSpline" => CurveKind::BSpline,
        other => CurveKind::Other(other.to_string()),
    }
}

fn str_to_surface_kind(s: &str) -> SurfaceKind {
    match s {
        "Plane" => SurfaceKind::Plane,
        "Cylinder" => SurfaceKind::Cylinder,
        "Sphere" => SurfaceKind::Sphere,
        "Cone" => SurfaceKind::Cone,
        "Torus" => SurfaceKind::Torus,
        "BSpline" => SurfaceKind::BSpline,
        other => SurfaceKind::Other(other.to_string()),
    }
}

fn has_exportable_shell_topology(brep: &rcad_kernel::BRep) -> bool {
    if brep.edges.is_empty() {
        return false;
    }
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if face.outer_wire.edges.is_empty() {
                    return false;
                }
                if face
                    .outer_wire
                    .edges
                    .iter()
                    .any(|we| we.idx >= brep.edges.len())
                {
                    return false;
                }
                for wire in &face.inner_wires {
                    if wire.edges.iter().any(|we| we.idx >= brep.edges.len()) {
                        return false;
                    }
                }
            }
        }
    }
    true
}

fn explode_by_connected_components_impl(brep: &rcad_kernel::BRep) -> Vec<rcad_kernel::BRep> {
    use rcad_algorithms::find_connected_components;

    let components = find_connected_components(brep);
    if components.len() <= 1 {
        return vec![brep.clone()];
    }

    let flat_faces: Vec<rcad_kernel::topology::Face> = brep
        .solids
        .iter()
        .flat_map(|solid| solid.shells.iter())
        .flat_map(|shell| shell.faces.iter().cloned())
        .collect();

    let mut out: Vec<rcad_kernel::BRep> = Vec::new();
    for component in components {
        if component.is_empty() {
            continue;
        }
        let mut faces: Vec<rcad_kernel::topology::Face> = Vec::with_capacity(component.len());
        let mut complete = true;
        for face_idx in component {
            if let Some(face) = flat_faces.get(face_idx) {
                faces.push(face.clone());
            } else {
                complete = false;
                break;
            }
        }
        if !complete || faces.is_empty() {
            continue;
        }
        if let Ok(solid) =
            rcad_algorithms::shape_build::BuildSolid::build_solid_from_faces(&faces, 1e-7)
        {
            let mut part = brep.clone();
            part.solids = vec![solid];
            out.push(part);
        }
    }
    if out.is_empty() {
        vec![brep.clone()]
    } else {
        out
    }
}

/// Frustum (truncated cone) construction with full analytic geometry.
fn frustum_brep_impl(
    center: DVec3,
    axis_norm: DVec3,
    ref_dir_hint: DVec3,
    r1: f64,
    r2: f64,
    height: f64,
) -> rcad_kernel::BRep {
    use rcad_kernel::geom::*;
    use rcad_kernel::topology::*;
    use rcad_kernel::{GeomStore, PCurve};
    use std::f64::consts::PI;

    let (axis_eff, rb, rt) = if r1 >= r2 {
        (axis_norm, r1, r2)
    } else {
        (-axis_norm, r2, r1)
    };

    let mut x_axis = ref_dir_hint - axis_eff * ref_dir_hint.dot(axis_eff);
    if x_axis.length_squared() < 1e-18 {
        let fallback = if axis_eff.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
        x_axis = fallback - axis_eff * fallback.dot(axis_eff);
    }
    x_axis = x_axis.normalize_or_zero();
    let z_axis = axis_eff.cross(x_axis).normalize_or_zero();

    let map_point = |p: DVec3| center + x_axis * p.x + axis_eff * p.y + z_axis * p.z;
    let half_h = height * 0.5;
    let tan_half = (rb - rt) / height;
    let half_angle = tan_half.atan();
    let apex_dist = rt / tan_half;
    let slant_len = ((rb - rt) * (rb - rt) + height * height).sqrt();
    let v_top = apex_dist / half_angle.cos();
    let v_base = v_top + slant_len;

    let top_pt = map_point(DVec3::new(rt, half_h, 0.0));
    let base_pt = map_point(DVec3::new(rb, -half_h, 0.0));

    let vertices = vec![Vertex { point: top_pt }, Vertex { point: base_pt }];
    let edges = vec![
        Edge { start: 0, end: 0 },
        Edge { start: 1, end: 1 },
        Edge { start: 0, end: 1 },
    ];

    let side_face = Face {
        outer_wire: Wire {
            edges: vec![WireEdge::fwd(2), WireEdge::rev(1), WireEdge::rev(2), WireEdge::fwd(0)],
        },
        inner_wires: vec![],
        normal: x_axis,
        triangles: vec![],
        mesh_dirty: true,
                sample_point: None,
                surface_idx: None,
    };
    let top_face = Face {
        outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
        inner_wires: vec![],
        normal: axis_eff,
        triangles: vec![],
        mesh_dirty: true,
                sample_point: None,
                surface_idx: None,
    };
    let base_face = Face {
        outer_wire: Wire { edges: vec![WireEdge::rev(1)] },
        inner_wires: vec![],
        normal: -axis_eff,
        triangles: vec![],
        mesh_dirty: true,
                sample_point: None,
                surface_idx: None,
    };

    let shell = Shell { faces: vec![side_face, top_face, base_face] };
    let solid = Solid { shells: vec![shell] };

    let top_center = map_point(DVec3::new(0.0, half_h, 0.0));
    let base_center = map_point(DVec3::new(0.0, -half_h, 0.0));
    let apex = map_point(DVec3::new(0.0, half_h + apex_dist, 0.0));

    let top_circle = Curve3::Circle(Circle3::new(top_center, axis_eff, rt));
    let base_circle = Curve3::Circle(Circle3::new(base_center, -axis_eff, rb));
    let seam_line = Curve3::Line(Line3 {
        origin: top_pt,
        direction: (base_pt - top_pt).normalize_or_zero(),
    });

    let side_surface = Surface3::Cone(ConicalSurface {
        apex,
        axis: -axis_eff,
        radius: 0.0,
        half_angle_rad: half_angle,
    });
    let top_plane = Surface3::Plane(Plane {
        origin: top_center, normal: axis_eff,
    });
    let base_plane = Surface3::Plane(Plane {
        origin: base_center, normal: -axis_eff,
    });

    let e0_on_side = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(2.0 * PI, v_top),
        direction: glam::DVec2::new(-1.0, 0.0),
    });
    let e0_on_top = Curve2d::Circle(Circle2d::new(glam::DVec2::ZERO, rt));
    let e1_on_side = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, v_base),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    let e1_on_base = Curve2d::Circle(Circle2d::new(glam::DVec2::ZERO, rb));
    let e2_on_side = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, v_top),
        direction: glam::DVec2::new(0.0, 1.0),
    });

    let geom = GeomStore {
        curves: vec![top_circle, base_circle, seam_line],
        surfaces: vec![side_surface, top_plane, base_plane],
        curve2ds: vec![e0_on_side, e0_on_top, e1_on_side, e1_on_base, e2_on_side],
        edge_curve: vec![Some(0), Some(1), Some(2)],
        face_surface: vec![Some(0), Some(1), Some(2)],
        edge_pcurves: vec![
            vec![PCurve { surface_idx: 0, curve2d_idx: 0 }, PCurve { surface_idx: 1, curve2d_idx: 1 }],
            vec![PCurve { surface_idx: 0, curve2d_idx: 2 }, PCurve { surface_idx: 2, curve2d_idx: 3 }],
            vec![PCurve { surface_idx: 0, curve2d_idx: 4 }],
        ],
        edge_curve_range: vec![Some([0.0, 2.0 * PI]), Some([0.0, 2.0 * PI]), Some([0.0, slant_len])],
        edge_degenerated: vec![false, false, false],
        vertex_tolerance: Vec::new(),
        edge_tolerance: Vec::new(),
        face_tolerance: Vec::new(),
        curve2d_range: Vec::new(),
        face_surface_range: Vec::new(),
        edge_same_parameter: Vec::new(),
        edge_same_range: Vec::new(),
        face_internal_vertices: Vec::new(),
        edge_vertex_params: Vec::new(),
    };

    rcad_kernel::BRep {
        vertices,
        edges,
        solids: vec![solid],
        geom,
        compound: None,
        compsolid: None,
    }
}

fn merge_breps_for_export_impl(shapes: &[rcad_kernel::BRep]) -> rcad_kernel::BRep {
    let mut out = rcad_kernel::BRep::new();

    for shape in shapes {
        let vertex_count = shape.vertices.len();
        let edge_count = shape.edges.len();
        let face_count: usize = shape
            .solids
            .iter()
            .flat_map(|solid| solid.shells.iter())
            .map(|shell| shell.faces.len())
            .sum();
        let curve2d_count = shape.geom.curve2ds.len();

        let vertex_offset = out.vertices.len();
        let edge_offset = out.edges.len();
        let curve_offset = out.geom.curves.len();
        let surface_offset = out.geom.surfaces.len();
        let curve2d_offset = out.geom.curve2ds.len();

        out.vertices.extend(shape.vertices.iter().cloned());
        out.edges.extend(shape.edges.iter().map(|e| rcad_kernel::topology::Edge {
            start: e.start + vertex_offset,
            end: e.end + vertex_offset,
        }));

        for mut solid in shape.solids.clone() {
            for shell in &mut solid.shells {
                for face in &mut shell.faces {
                    for we in &mut face.outer_wire.edges {
                        we.idx += edge_offset;
                    }
                    for iw in &mut face.inner_wires {
                        for we in &mut iw.edges {
                            we.idx += edge_offset;
                        }
                    }
                    for tri in &mut face.triangles {
                        tri[0] += vertex_offset;
                        tri[1] += vertex_offset;
                        tri[2] += vertex_offset;
                    }
                }
            }
            out.solids.push(solid);
        }

        out.geom.curves.extend(shape.geom.curves.iter().cloned());
        out.geom.surfaces.extend(shape.geom.surfaces.iter().cloned());
        out.geom.curve2ds.extend(shape.geom.curve2ds.iter().cloned());

        for ei in 0..edge_count {
            out.geom.edge_curve.push(
                shape.geom.edge_curve.get(ei).copied().flatten().map(|i| i + curve_offset),
            );
            out.geom.edge_curve_range.push(
                shape.geom.edge_curve_range.get(ei).copied().flatten(),
            );
            out.geom.edge_degenerated.push(
                shape.geom.edge_degenerated.get(ei).copied().unwrap_or(false),
            );
            out.geom.edge_tolerance.push(
                shape.geom.edge_tolerance.get(ei).copied().unwrap_or(0.0),
            );
            out.geom.edge_same_parameter.push(
                shape.geom.edge_same_parameter.get(ei).copied().unwrap_or(false),
            );
            out.geom.edge_same_range.push(
                shape.geom.edge_same_range.get(ei).copied().unwrap_or(false),
            );

            let pcs = shape
                .geom
                .edge_pcurves
                .get(ei)
                .map(|pcs| {
                    pcs.iter()
                        .map(|pc| rcad_kernel::PCurve {
                            surface_idx: pc.surface_idx + surface_offset,
                            curve2d_idx: pc.curve2d_idx + curve2d_offset,
                        })
                        .collect()
                })
                .unwrap_or_default();
            out.geom.edge_pcurves.push(pcs);
        }

        for fi in 0..face_count {
            out.geom.face_surface.push(
                shape.geom.face_surface.get(fi).copied().flatten().map(|i| i + surface_offset),
            );
            out.geom.face_tolerance.push(
                shape.geom.face_tolerance.get(fi).copied().unwrap_or(0.0),
            );
            out.geom.face_surface_range.push(
                shape.geom.face_surface_range.get(fi).copied().unwrap_or(None),
            );
        }

        for vi in 0..vertex_count {
            out.geom.vertex_tolerance.push(
                shape.geom.vertex_tolerance.get(vi).copied().unwrap_or(0.0),
            );
        }

        for ci in 0..curve2d_count {
            out.geom.curve2d_range.push(
                shape.geom.curve2d_range.get(ci).copied().unwrap_or(None),
            );
        }
    }

    out
}

fn normalize_for_strict_step_export(brep: &rcad_kernel::BRep) -> rcad_kernel::BRep {
    use rcad_kernel::geom::{Curve3, Line3};

    let mut out = brep.clone();

    if out.geom.edge_curve.len() < out.edges.len() {
        out.geom.edge_curve.resize(out.edges.len(), None);
    }
    if out.geom.edge_curve_range.len() < out.edges.len() {
        out.geom.edge_curve_range.resize(out.edges.len(), None);
    }

    let mut referenced = vec![false; out.edges.len()];
    for solid in &out.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    if we.idx < referenced.len() {
                        referenced[we.idx] = true;
                    }
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        if we.idx < referenced.len() {
                            referenced[we.idx] = true;
                        }
                    }
                }
            }
        }
    }

    for (edge_idx, edge) in out.edges.iter().enumerate() {
        if !referenced[edge_idx] || out.geom.edge_curve[edge_idx].is_some() {
            continue;
        }
        let Some(ps) = out.vertices.get(edge.start).map(|v| v.point) else { continue; };
        let Some(pe) = out.vertices.get(edge.end).map(|v| v.point) else { continue; };
        let d = pe - ps;
        let len = d.length();
        if !len.is_finite() || len <= 1e-12 {
            continue;
        }
        let dir = d / len;
        let cid = out.geom.curves.len();
        out.geom.curves.push(Curve3::Line(Line3 {
            origin: ps,
            direction: dir,
        }));
        out.geom.edge_curve[edge_idx] = Some(cid);
        out.geom.edge_curve_range[edge_idx] = Some([0.0, len]);
    }

    let edge_curve = out.geom.edge_curve.clone();
    for solid in &mut out.solids {
        for shell in &mut solid.shells {
            for face in &mut shell.faces {
                sanitize_wire_for_strict_export(&mut face.outer_wire, &out.edges, &edge_curve);
                for wire in &mut face.inner_wires {
                    sanitize_wire_for_strict_export(wire, &out.edges, &edge_curve);
                }
            }
        }
    }

    out
}

fn sanitize_wire_for_strict_export(
    wire: &mut rcad_kernel::topology::Wire,
    edges: &[rcad_kernel::topology::Edge],
    edge_curve: &[Option<usize>],
) {
    let mut filtered = Vec::with_capacity(wire.edges.len());
    for we in &wire.edges {
        if we.idx >= edges.len() {
            continue;
        }
        let has_curve = edge_curve.get(we.idx).and_then(|v| *v).is_some();
        if !has_curve {
            let e = &edges[we.idx];
            if e.start == e.end {
                continue;
            }
        }
        if filtered
            .last()
            .is_some_and(|prev: &rcad_kernel::topology::WireEdge| prev.idx == we.idx && prev.forward == we.forward)
        {
            continue;
        }
        filtered.push(*we);
    }
    wire.edges = filtered;
}

