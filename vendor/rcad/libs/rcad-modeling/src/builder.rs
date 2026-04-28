//! User-facing geometry construction helpers.
//!
//! This module is the public modeling entry layer for RCAD.
//! The API intentionally prefers OCCT-style direct constructor functions
//! over fluent builder structs.

pub mod brep_builder;
mod curve;
pub mod fillet;
pub mod ops;
mod solid;
mod surface;
pub mod wire_ops;

pub use brep_builder::*;
pub use curve::*;
pub use fillet::{chamfer_edge, chamfer_edge_angle, chamfer_edge_safe, corner_blend, fillet_edge, fillet_edge_safe, fillet_edges, fillet_edge_variable_radius};
pub use fillet::{chamfer_edge_with_history, chamfer_edge_angle_with_history, fillet_edge_with_history, fillet_edge_variable_radius_with_history, fillet_edges_with_history, corner_blend_with_history};
pub use fillet::{FilletHistory, MultiFilletHistory, CornerBlendHistory, SafeFilletResult};
pub use ops::*;
pub use solid::*;
pub use surface::*;
pub use wire_ops::{chamfer_wire_2d, fillet_wire_2d, project_wire_onto_surface};

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Surface3};
use std::error::Error;
use std::fmt;

const EPS: f64 = 1e-12;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BuildError {
    NonFiniteValue(&'static str),
    NonPositiveValue(&'static str),
    ZeroVector(&'static str),
    ParallelVectors(&'static str, &'static str),
    DegenerateGeometry(&'static str),
    InvalidIndex(usize),
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteValue(name) => write!(f, "{name} must be finite"),
            Self::NonPositiveValue(name) => write!(f, "{name} must be > 0"),
            Self::ZeroVector(name) => write!(f, "{name} must be non-zero"),
            Self::ParallelVectors(a, b) => write!(f, "{a} must not be parallel to {b}"),
            Self::DegenerateGeometry(msg) => write!(f, "degenerate geometry: {msg}"),
            Self::InvalidIndex(idx) => write!(f, "invalid index: {idx}"),
        }
    }
}

impl Error for BuildError {}

fn validate_point(name: &'static str, point: DVec3) -> Result<DVec3, BuildError> {
    if point.is_finite() {
        Ok(point)
    } else {
        Err(BuildError::NonFiniteValue(name))
    }
}

fn validate_positive(name: &'static str, value: f64) -> Result<f64, BuildError> {
    if !value.is_finite() {
        Err(BuildError::NonFiniteValue(name))
    } else if value <= 0.0 {
        Err(BuildError::NonPositiveValue(name))
    } else {
        Ok(value)
    }
}

fn normalize_vector(name: &'static str, vector: DVec3) -> Result<DVec3, BuildError> {
    validate_point(name, vector)?;
    if vector.length_squared() <= EPS {
        Err(BuildError::ZeroVector(name))
    } else {
        Ok(vector.normalize())
    }
}

fn normalize_rejection(
    name: &'static str,
    vector: DVec3,
    reference_name: &'static str,
    reference: DVec3,
) -> Result<DVec3, BuildError> {
    let vector = normalize_vector(name, vector)?;
    let rejected = vector - reference * vector.dot(reference);
    if rejected.length_squared() <= EPS {
        Err(BuildError::ParallelVectors(name, reference_name))
    } else {
        Ok(rejected.normalize())
    }
}

fn basis_from_x_y(x_dir: DVec3, y_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), BuildError> {
    let x_axis = normalize_vector("x_dir", x_dir)?;
    let y_axis = normalize_rejection("y_dir", y_dir, "x_dir", x_axis)?;
    let z_axis = x_axis.cross(y_axis).normalize();
    Ok((x_axis, y_axis, z_axis))
}

fn basis_from_axis_ref(axis: DVec3, ref_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), BuildError> {
    let y_axis = normalize_vector("axis", axis)?;
    let x_axis = normalize_rejection("ref_dir", ref_dir, "axis", y_axis)?;
    let z_axis = x_axis.cross(y_axis).normalize();
    Ok((x_axis, y_axis, z_axis))
}

fn translate_brep(brep: &mut BRep, offset: DVec3) {
    // Translate vertices
    for vertex in &mut brep.vertices {
        vertex.point += offset;
    }
    // Translate analytic geometry in GeomStore
    for curve in &mut brep.geom.curves {
        match curve {
            Curve3::Line(l) => l.origin += offset,
            Curve3::Circle(c) => c.center += offset,
            Curve3::Ellipse(e) => e.center += offset,
            Curve3::Hyperbola(h) => h.center += offset,
            Curve3::BSpline(b) => {
                for cp in &mut b.control_points {
                    *cp += offset;
                }
            }
            Curve3::Bezier(b) => {
                for cp in &mut b.control_points {
                    *cp += offset;
                }
            }
            _ => {}
        }
    }
    for surface in &mut brep.geom.surfaces {
        match surface {
            Surface3::Plane(p) => p.origin += offset,
            Surface3::Cylinder(c) => c.origin += offset,
            Surface3::Sphere(s) => s.center += offset,
            Surface3::Cone(c) => c.apex += offset,
            Surface3::Torus(t) => t.center += offset,
            Surface3::BSpline(b) => {
                for row in &mut b.control_points {
                    for cp in row {
                        *cp += offset;
                    }
                }
            }
            _ => {}
        }
    }
    for solid in &mut brep.solids {
        for shell in &mut solid.shells {
            for face in &mut shell.faces {
                face.normal = face.normal.normalize_or_zero();
            }
        }
    }
}

fn do_mirror_brep(brep: &BRep, plane_origin: DVec3, plane_normal: DVec3) -> BRep {
    let n = plane_normal.normalize();
    let mirror_point = |p: DVec3| -> DVec3 {
        let d = (p - plane_origin).dot(n);
        p - n * (2.0 * d)
    };
    let mirror_vec = |v: DVec3| -> DVec3 {
        // Reflect direction: v - 2*(v·n)*n
        v - n * (2.0 * v.dot(n))
    };

    let mut out = BRep::new();

    // Mirror vertices
    for vertex in &brep.vertices {
        out.vertices.push(rcad_kernel::topology::Vertex {
            point: mirror_point(vertex.point),
        });
    }

    // Mirror curves
    for curve in &brep.geom.curves {
        out.geom.curves.push(match curve {
            Curve3::Line(l) => Curve3::Line(rcad_kernel::geom::Line3 {
                origin: mirror_point(l.origin),
                direction: mirror_vec(l.direction),
            }),
            Curve3::Circle(c) => Curve3::Circle(rcad_kernel::geom::Circle3 {
                center: mirror_point(c.center),
                normal: mirror_vec(c.normal),
                radius: c.radius,
            }),
            Curve3::Ellipse(e) => Curve3::Ellipse(rcad_kernel::geom::Ellipse3 {
                center: mirror_point(e.center),
                normal: mirror_vec(e.normal),
                major_dir: mirror_vec(e.major_dir),
                major_radius: e.major_radius,
                minor_radius: e.minor_radius,
            }),
            Curve3::Hyperbola(h) => Curve3::Hyperbola(rcad_kernel::geom::Hyperbola3 {
                center: mirror_point(h.center),
                normal: mirror_vec(h.normal),
                major_dir: mirror_vec(h.major_dir),
                semi_major: h.semi_major,
                semi_minor: h.semi_minor,
            }),
            Curve3::BSpline(b) => {
                let mut nb = b.clone();
                for cp in &mut nb.control_points {
                    *cp = mirror_point(*cp);
                }
                Curve3::BSpline(nb)
            }
            Curve3::Bezier(b) => {
                let mut nb = b.clone();
                for cp in &mut nb.control_points {
                    *cp = mirror_point(*cp);
                }
                Curve3::Bezier(nb)
            }
            _ => curve.clone(),
        });
    }

    // Mirror surfaces
    for surface in &brep.geom.surfaces {
        out.geom.surfaces.push(match surface {
            Surface3::Plane(p) => Surface3::Plane(rcad_kernel::geom::Plane {
                origin: mirror_point(p.origin),
                normal: mirror_vec(p.normal),
            }),
            Surface3::Cylinder(c) => Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
                origin: mirror_point(c.origin),
                axis: mirror_vec(c.axis),
                radius: c.radius,
            }),
            Surface3::Sphere(s) => Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
                center: mirror_point(s.center),
                axis: mirror_vec(s.axis),
                radius: s.radius,
            }),
            Surface3::Cone(c) => Surface3::Cone(rcad_kernel::geom::ConicalSurface {
                apex: mirror_point(c.apex),
                axis: mirror_vec(c.axis),
                radius: c.radius,
                half_angle_rad: c.half_angle_rad,
            }),
            Surface3::Torus(t) => Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
                center: mirror_point(t.center),
                axis: mirror_vec(t.axis),
                major_radius: t.major_radius,
                minor_radius: t.minor_radius,
            }),
            Surface3::BSpline(b) => {
                let mut nb = b.clone();
                for row in &mut nb.control_points {
                    for cp in row {
                        *cp = mirror_point(*cp);
                    }
                }
                Surface3::BSpline(nb)
            }
            _ => surface.clone(),
        });
    }

    // Mirror edges (vertex indices are 1:1 mapped)
    for e in &brep.edges {
        out.edges.push(rcad_kernel::topology::Edge {
            start: e.start,
            end: e.end,
        });
    }
    out.geom.edge_curve = brep.geom.edge_curve.clone();
    out.geom.edge_curve_range = brep.geom.edge_curve_range.clone();
    out.geom.edge_degenerated = brep.geom.edge_degenerated.clone();

    // Mirror faces — flip normal (mirror reverses orientation)
    for solid in &brep.solids {
        let mut shell = rcad_kernel::topology::Shell { faces: Vec::new() };
        for face in &solid.shells[0].faces {
            let wire_edges: Vec<rcad_kernel::topology::WireEdge> = face
                .outer_wire
                .edges
                .iter()
                .map(|we| rcad_kernel::topology::WireEdge {
                    idx: we.idx,
                    forward: !we.forward, // flip orientation for mirror
                })
                .collect();
            shell.faces.push(rcad_kernel::topology::Face {
                outer_wire: rcad_kernel::topology::Wire { edges: wire_edges },
                inner_wires: face.inner_wires.clone(),
                normal: mirror_vec(face.normal),
                // Flip triangle winding order to maintain outward normals after mirror
                triangles: face.triangles.iter().map(|[i, j, k]| [*i, *k, *j]).collect(),
                mesh_dirty: face.mesh_dirty,
            });
        }
        out.solids.push(rcad_kernel::topology::Solid {
            shells: vec![shell],
        });
    }

    // Mirror face_surface references
    out.geom.face_surface = brep.geom.face_surface.clone();

    out
}

fn transform_brep(brep: &mut BRep, origin: DVec3, x_axis: DVec3, y_axis: DVec3, z_axis: DVec3) {
    let xform_point = |p: DVec3| -> DVec3 { origin + x_axis * p.x + y_axis * p.y + z_axis * p.z };
    let xform_vec =
        |v: DVec3| -> DVec3 { (x_axis * v.x + y_axis * v.y + z_axis * v.z).normalize_or_zero() };

    // Transform vertices
    for vertex in &mut brep.vertices {
        vertex.point = xform_point(vertex.point);
    }

    // Transform analytic geometry in GeomStore
    for curve in &mut brep.geom.curves {
        match curve {
            Curve3::Line(l) => {
                l.origin = xform_point(l.origin);
                l.direction = xform_vec(l.direction);
            }
            Curve3::Circle(c) => {
                c.center = xform_point(c.center);
                c.normal = xform_vec(c.normal);
            }
            Curve3::Ellipse(e) => {
                e.center = xform_point(e.center);
                e.normal = xform_vec(e.normal);
                e.major_dir = xform_vec(e.major_dir);
            }
            Curve3::Hyperbola(h) => {
                h.center = xform_point(h.center);
                h.normal = xform_vec(h.normal);
                h.major_dir = xform_vec(h.major_dir);
            }
            Curve3::BSpline(b) => {
                for cp in &mut b.control_points {
                    *cp = xform_point(*cp);
                }
            }
            Curve3::Bezier(b) => {
                for cp in &mut b.control_points {
                    *cp = xform_point(*cp);
                }
            }
            _ => {}
        }
    }
    for surface in &mut brep.geom.surfaces {
        match surface {
            Surface3::Plane(p) => {
                p.origin = xform_point(p.origin);
                p.normal = xform_vec(p.normal);
            }
            Surface3::Cylinder(c) => {
                c.origin = xform_point(c.origin);
                c.axis = xform_vec(c.axis);
            }
            Surface3::Sphere(s) => {
                s.center = xform_point(s.center);
                s.axis = xform_vec(s.axis);
            }
            Surface3::Cone(c) => {
                c.apex = xform_point(c.apex);
                c.axis = xform_vec(c.axis);
            }
            Surface3::Torus(t) => {
                t.center = xform_point(t.center);
                t.axis = xform_vec(t.axis);
            }
            Surface3::BSpline(b) => {
                for row in &mut b.control_points {
                    for cp in row {
                        *cp = xform_point(*cp);
                    }
                }
            }
            _ => {}
        }
    }

    for solid in &mut brep.solids {
        for shell in &mut solid.shells {
            for face in &mut shell.faces {
                let transformed =
                    x_axis * face.normal.x + y_axis * face.normal.y + z_axis * face.normal.z;
                face.normal = transformed.normalize_or_zero();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::{PrimitiveSolid, Surface3};

    #[test]
    fn line_rejects_zero_direction() {
        let err = line(DVec3::ZERO, DVec3::ZERO).unwrap_err();
        assert_eq!(err, BuildError::ZeroVector("direction"));
    }

    #[test]
    fn ellipse_rejects_parallel_major_direction() {
        let err = ellipse(DVec3::ZERO, DVec3::Z, DVec3::Z, 2.0, 1.0).unwrap_err();
        assert_eq!(err, BuildError::ParallelVectors("major_dir", "normal"));
    }

    #[test]
    fn box_brep_builds_transformed_vertices() {
        let brep = box_brep(DVec3::new(1.0, 2.0, 3.0), DVec3::Y, DVec3::Z, 2.0, 3.0, 4.0).unwrap();

        assert_eq!(brep.vertices.len(), 8);
        assert!(
            brep.vertices
                .iter()
                .any(|v| v.point == DVec3::new(1.0, 2.0, 3.0))
        );
        assert!(
            brep.vertices
                .iter()
                .any(|v| v.point == DVec3::new(5.0, 4.0, 6.0))
        );
    }

    #[test]
    fn sphere_brep_translates_bounds() {
        let brep = sphere_brep(DVec3::new(10.0, -2.0, 4.0), 2.0).unwrap();

        let min_y = brep
            .vertices
            .iter()
            .map(|v| v.point.y)
            .fold(f64::INFINITY, f64::min);
        let max_y = brep
            .vertices
            .iter()
            .map(|v| v.point.y)
            .fold(f64::NEG_INFINITY, f64::max);

        assert!((min_y - (-4.0)).abs() < 1e-6);
        assert!((max_y - 0.0).abs() < 1e-6);
    }

    #[test]
    fn cylinder_primitive_returns_expected_shape() {
        let primitive = cylinder_primitive(3.0, 5.0).unwrap();

        match primitive {
            PrimitiveSolid::Cylinder { radius, height } => {
                assert_eq!(radius, 3.0);
                assert_eq!(height, 5.0);
            }
            other => panic!("expected cylinder primitive, got {other:?}"),
        }
    }

    #[test]
    fn make_plane_alias_matches_plane_constructor() {
        let surface = make_plane(DVec3::new(1.0, 2.0, 3.0), DVec3::Z).unwrap();

        match surface {
            Surface3::Plane(plane) => {
                assert_eq!(plane.origin, DVec3::new(1.0, 2.0, 3.0));
                assert_eq!(plane.normal, DVec3::Z);
            }
            other => panic!("expected plane surface, got {other:?}"),
        }
    }

    #[test]
    fn mirror_box_across_xy_plane() {
        let mut brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        // Populate geometry manually (geom_populate lives in rcad-algorithms)
        use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
        for edge in &brep.edges {
            let p0 = brep.vertices[edge.start].point;
            let p1 = brep.vertices[edge.end].point;
            let delta = p1 - p0;
            let len = delta.length();
            let dir = if len > 1e-12 { delta / len } else { DVec3::X };
            let curve_idx = brep.geom.curves.len();
            brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
            brep.geom.edge_curve.push(Some(curve_idx));
            brep.geom.edge_curve_range.push(Some([0.0, len]));
            brep.geom.edge_degenerated.push(false);
        }
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let origin = face.outer_wire.edges.first()
                        .and_then(|we| brep.edges.get(we.idx))
                        .map(|e| brep.vertices[e.start].point)
                        .unwrap_or(DVec3::ZERO);
                    let surf_idx = brep.geom.surfaces.len();
                    brep.geom.surfaces.push(Surface3::Plane(Plane { origin, normal: face.normal }));
                    let _face_idx = brep.geom.face_surface.len();
                    brep.geom.face_surface.push(Some(surf_idx));
                }
            }
        }
        let v_orig = rcad_kernel::properties::volume(&brep);

        let mirrored = do_mirror_brep(&brep, DVec3::ZERO, DVec3::Z);
        let v_mirrored = rcad_kernel::properties::volume(&mirrored);

        // Debug: check face triangles
        for (fi, face) in mirrored.solids[0].shells[0].faces.iter().enumerate() {
            eprintln!("mirrored face {fi}: normal={:?} tris={}", face.normal, face.triangles.len());
            for tri in &face.triangles {
                let a = mirrored.vertices[tri[0]].point;
                let b = mirrored.vertices[tri[1]].point;
                let c = mirrored.vertices[tri[2]].point;
                let gn = (b - a).cross(c - a);
                eprintln!("  tri: gn={gn:?}");
            }
        }

        assert!(
            (v_mirrored - v_orig).abs() < 0.01,
            "mirror should preserve volume: {v_orig} vs {v_mirrored}"
        );

        for (i, v) in brep.vertices.iter().enumerate() {
            let mv = &mirrored.vertices[i];
            assert!((mv.point.z - (-v.point.z)).abs() < 1e-9, "vertex {i} z should be negated");
        }
    }
}
