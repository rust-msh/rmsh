//! `model.occ.*` bindings: OCC-style CAD construction over the real rcad
//! kernel — primitives from rcad-modeling, booleans from the OCCT-aligned
//! rcad-algo BOP pipeline, transforms via `BRep::apply_transform`, STEP
//! import via rcad-step and surface synchronization through the
//! rmsh-bridge B-Rep walker.

use std::collections::BTreeMap;
use std::sync::{Arc, RwLock};

use glam::{DAffine3, DMat3, DVec3, DQuat};
use pyo3::exceptions::{PyNotImplementedError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;

use crate::state::{with_current_model, OccShape};

fn shape_brep(model: &mut crate::state::Model, tag: i32) -> PyResult<&mut rcad_kernel::BRep> {
    match model.shapes.get_mut(&tag) {
        Some(OccShape::BRep { brep, .. }) => Ok(brep),
        other => Err(PyValueError::new_err(format!(
            "entity {tag} is not a CAD shape (found dimension {:?})",
            other.map(|s| s.dimension()).unwrap_or(-1)
        ))),
    }
}

fn insert_brep(model: &mut crate::state::Model, dim: i32, brep: rcad_kernel::BRep) -> i32 {
    let tag = model.alloc_tag();
    model.shapes.insert(tag, OccShape::BRep { dim, brep });
    tag
}

fn build_err(e: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(format!("rcad primitive build failed: {e}"))
}

fn ref_dir_for(axis: DVec3) -> DVec3 {
    let a = axis.normalize_or_zero();
    let candidate = if a.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    candidate.cross(a).normalize_or_zero()
}

#[pyfunction]
#[pyo3(signature = (x, y, z, dx, dy, dz, tag=-1))]
pub fn model_occ_add_box(
    x: f64, y: f64, z: f64, dx: f64, dy: f64, dz: f64, tag: i32,
) -> PyResult<i32> {
    with_current_model(|model| {
        let mut brep = rcad_modeling::prim::primapi::box_brep(dx, dy, dz).map_err(build_err)?;
        brep.apply_transform(DAffine3::from_translation(DVec3::new(x, y, z)));
        Ok(match tag {
            -1 => insert_brep(model, 3, brep),
            t => {
                model.shapes.insert(t, OccShape::BRep { dim: 3, brep });
                t
            }
        })
    })
}

#[pyfunction]
#[pyo3(signature = (x, y, z, r, tag=-1))]
pub fn model_occ_add_sphere(x: f64, y: f64, z: f64, r: f64, tag: i32) -> PyResult<i32> {
    with_current_model(|model| {
        let brep =
            rcad_modeling::prim::primapi::sphere_brep(DVec3::new(x, y, z), r).map_err(build_err)?;
        Ok(match tag {
            -1 => insert_brep(model, 3, brep),
            t => {
                model.shapes.insert(t, OccShape::BRep { dim: 3, brep });
                t
            }
        })
    })
}

/// gmsh: (x, y, z) is the base center, (dx, dy, dz) the axis vector whose
/// LENGTH is the cylinder height.
#[pyfunction]
#[pyo3(signature = (x, y, z, dx, dy, dz, r, tag=-1))]
pub fn model_occ_add_cylinder(
    x: f64, y: f64, z: f64, dx: f64, dy: f64, dz: f64, r: f64, tag: i32,
) -> PyResult<i32> {
    with_current_model(|model| {
        let axis = DVec3::new(dx, dy, dz);
        let height = axis.length();
        if height <= 0.0 {
            return Err(PyValueError::new_err("cylinder axis vector is zero"));
        }
        let brep = rcad_modeling::prim::primapi::cylinder_brep(
            DVec3::new(x, y, z),
            axis,
            ref_dir_for(axis),
            r,
            height,
        )
        .map_err(build_err)?;
        Ok(match tag {
            -1 => insert_brep(model, 3, brep),
            t => {
                model.shapes.insert(t, OccShape::BRep { dim: 3, brep });
                t
            }
        })
    })
}

#[pyfunction]
#[pyo3(signature = (x, y, z, dx, dy, dz, r1, r2, tag=-1))]
pub fn model_occ_add_cone(
    x: f64, y: f64, z: f64, dx: f64, dy: f64, dz: f64, r1: f64, r2: f64, tag: i32,
) -> PyResult<i32> {
    with_current_model(|model| {
        let axis = DVec3::new(dx, dy, dz);
        let height = axis.length();
        if height <= 0.0 {
            return Err(PyValueError::new_err("cone axis vector is zero"));
        }
        let brep = rcad_modeling::prim::primapi::cone_brep(
            DVec3::new(x, y, z),
            axis,
            ref_dir_for(axis),
            r1,
            r2,
            height,
        )
        .map_err(build_err)?;
        Ok(match tag {
            -1 => insert_brep(model, 3, brep),
            t => {
                model.shapes.insert(t, OccShape::BRep { dim: 3, brep });
                t
            }
        })
    })
}

/// Legacy rmsh signature `(x, y, z, dx, dy, dz, r1, r2, tag=-1)`: the axis
/// vector form.
#[pyfunction]
#[pyo3(signature = (x, y, z, dx, dy, dz, r1, r2, tag=-1))]
pub fn model_occ_add_torus(
    x: f64, y: f64, z: f64, dx: f64, dy: f64, dz: f64, r1: f64, r2: f64, tag: i32,
) -> PyResult<i32> {
    with_current_model(|model| {
        let axis = DVec3::new(dx, dy, dz);
        if axis.length() <= 0.0 {
            return Err(PyValueError::new_err("torus axis vector is zero"));
        }
        let brep = rcad_modeling::prim::primapi::torus_brep(
            DVec3::new(x, y, z),
            axis,
            ref_dir_for(axis),
            r1,
            r2,
        )
        .map_err(build_err)?;
        Ok(match tag {
            -1 => insert_brep(model, 3, brep),
            t => {
                model.shapes.insert(t, OccShape::BRep { dim: 3, brep });
                t
            }
        })
    })
}

#[pyfunction]
#[pyo3(signature = (x, y, z, dx, dy, tag=-1))]
pub fn model_occ_add_rectangle(
    x: f64, y: f64, z: f64, dx: f64, dy: f64, tag: i32,
) -> PyResult<i32> {
    with_current_model(|model| {
        let brep = rcad_modeling::prim::primapi::make_planar_rect_brep(
            DVec3::new(x, y, z),
            DVec3::X,
            DVec3::Y,
            0.0,
            dx,
            0.0,
            dy,
        )
        .map_err(build_err)?;
        Ok(match tag {
            -1 => insert_brep(model, 2, brep),
            t => {
                model.shapes.insert(t, OccShape::BRep { dim: 2, brep });
                t
            }
        })
    })
}

/// v1 fidelity note: the disk is a 64-gon planar face (rcad has no analytic
/// ellipse face primitive yet) - the same documented-deviation pattern as
/// the P3 corpus field discretization.
#[pyfunction]
#[pyo3(signature = (xc, yc, zc, rx, ry, tag=-1))]
pub fn model_occ_add_disk(xc: f64, yc: f64, zc: f64, rx: f64, ry: f64, tag: i32) -> PyResult<i32> {
    with_current_model(|model| {
        const N: usize = 64;
        let outer: Vec<DVec3> = (0..N)
            .map(|i| {
                let t = 2.0 * std::f64::consts::PI * (i as f64) / (N as f64);
                DVec3::new(xc + rx * t.cos(), yc + ry * t.sin(), zc)
            })
            .collect();
        let brep = rcad_modeling::prim::primapi::make_planar_polygon_brep(
            DVec3::new(xc, yc, zc),
            DVec3::Z,
            &outer,
            &[],
        )
        .map_err(build_err)?;
        Ok(match tag {
            -1 => insert_brep(model, 2, brep),
            t => {
                model.shapes.insert(t, OccShape::BRep { dim: 2, brep });
                t
            }
        })
    })
}

// ── construction geometry (tagged, reported by getEntities, ignored by
//    meshing in v1) ────────────────────────────────────────────────────

#[pyfunction]
#[pyo3(signature = (x, y, z, mesh_size=0.0, tag=-1))]
pub fn model_occ_add_point(x: f64, y: f64, z: f64, mesh_size: f64, tag: i32) -> PyResult<i32> {
    let _ = mesh_size; // v1: point mesh sizes do not feed the field
    with_current_model(|model| {
        let t = match tag {
            -1 => model.alloc_tag(),
            t => t,
        };
        model.shapes.insert(t, OccShape::Point(DVec3::new(x, y, z)));
        Ok(t)
    })
}

#[pyfunction]
#[pyo3(signature = (start_tag, end_tag, tag=-1))]
pub fn model_occ_add_line(start_tag: i32, end_tag: i32, tag: i32) -> PyResult<i32> {
    with_current_model(|model| {
        let (a, b) = match (model.shapes.get(&start_tag), model.shapes.get(&end_tag)) {
            (Some(OccShape::Point(a)), Some(OccShape::Point(b))) => (*a, *b),
            _ => {
                return Err(PyValueError::new_err(
                    "addLine endpoints must be point entities",
                ))
            }
        };
        let t = match tag {
            -1 => model.alloc_tag(),
            t => t,
        };
        model.shapes.insert(t, OccShape::Line(a, b));
        Ok(t)
    })
}

#[pyfunction]
#[pyo3(signature = (start_tag, center_tag, end_tag, tag=-1))]
pub fn model_occ_add_circle_arc(
    start_tag: i32, center_tag: i32, end_tag: i32, tag: i32,
) -> PyResult<i32> {
    with_current_model(|model| {
        let pts = match (
            model.shapes.get(&start_tag),
            model.shapes.get(&center_tag),
            model.shapes.get(&end_tag),
        ) {
            (
                Some(OccShape::Point(a)),
                Some(OccShape::Point(c)),
                Some(OccShape::Point(b)),
            ) => (*a, *c, *b),
            _ => return Err(PyValueError::new_err("addCircleArc requires point entities")),
        };
        let t = match tag {
            -1 => model.alloc_tag(),
            t => t,
        };
        model.shapes.insert(t, OccShape::Arc(pts.0, pts.1, pts.2));
        Ok(t)
    })
}

#[pyfunction]
#[pyo3(signature = (point_tags, tag=-1))]
pub fn model_occ_add_spline(point_tags: Vec<i32>, tag: i32) -> PyResult<i32> {
    with_current_model(|model| {
        let mut pts = Vec::new();
        for t in &point_tags {
            match model.shapes.get(t) {
                Some(OccShape::Point(p)) => pts.push(*p),
                _ => return Err(PyValueError::new_err("addSpline requires point entities")),
            }
        }
        let t = match tag {
            -1 => model.alloc_tag(),
            t => t,
        };
        model.shapes.insert(t, OccShape::Spline(pts));
        Ok(t)
    })
}

#[pyfunction]
#[pyo3(signature = (curve_tags, tag=-1))]
pub fn model_occ_add_curve_loop(curve_tags: Vec<i32>, tag: i32) -> PyResult<i32> {
    with_current_model(|model| {
        let t = match tag {
            -1 => model.alloc_tag(),
            t => t,
        };
        model
            .shapes
            .insert(t, OccShape::Loop { dim: 1, members: curve_tags });
        Ok(t)
    })
}

#[pyfunction]
#[pyo3(signature = (wire_tags, tag=-1))]
pub fn model_occ_add_plane_surface(wire_tags: Vec<i32>, tag: i32) -> PyResult<i32> {
    with_current_model(|model| {
        // v1: a plane surface from exactly one curve loop of addLine
        // segments builds a planar polygon face; anything richer is not
        // wired yet.
        let surface_err = || {
            PyNotImplementedError::new_err(
                "occ.addPlaneSurface v1 supports a single curve loop of addLine segments",
            )
        };
        if wire_tags.len() != 1 {
            return Err(surface_err());
        }
        let loop_members = match model.shapes.get(&wire_tags[0]) {
            Some(OccShape::Loop { members, .. }) => members.clone(),
            _ => return Err(surface_err()),
        };
        if loop_members.len() < 3 {
            return Err(surface_err());
        }
        let mut pts = Vec::new();
        for c in &loop_members {
            match model.shapes.get(c) {
                Some(OccShape::Line(a, _)) => pts.push(*a),
                _ => return Err(surface_err()),
            }
        }
        let origin = pts[0];
        let normal = (pts[1] - origin).cross(pts[2] - origin).normalize_or_zero();
        let brep = rcad_modeling::prim::primapi::make_planar_polygon_brep(origin, normal, &pts, &[])
            .map_err(build_err)?;
        Ok(match tag {
            -1 => insert_brep(model, 2, brep),
            t => {
                model.shapes.insert(t, OccShape::BRep { dim: 2, brep });
                t
            }
        })
    })
}

#[pyfunction]
#[pyo3(signature = (surface_tags, tag=-1))]
pub fn model_occ_add_surface_loop(surface_tags: Vec<i32>, tag: i32) -> PyResult<i32> {
    with_current_model(|model| {
        let t = match tag {
            -1 => model.alloc_tag(),
            t => t,
        };
        model
            .shapes
            .insert(t, OccShape::Loop { dim: 2, members: surface_tags });
        Ok(t)
    })
}

#[pyfunction]
#[pyo3(signature = (shell_tags, tag=-1))]
pub fn model_occ_add_volume(shell_tags: Vec<i32>, tag: i32) -> PyResult<i32> {
    let _ = (shell_tags, tag);
    Err(PyNotImplementedError::new_err(
        "occ.addVolume from surface loops requires shell sewing - not wired in v1; \
         use the primitive constructors or cut/fuse/intersect",
    ))
}

// ── booleans ─────────────────────────────────────────────────────────
//
// v1 folds lists sequentially (object list fused, then tools applied one
// by one) instead of a general n-ary BOP - a documented fidelity note.

#[derive(Clone, Copy)]
enum BooleanFlavor {
    Cut,
    Fuse,
    Intersect,
}

fn fold_boolean(
    objects: Vec<(i32, i32)>,
    tools: Vec<(i32, i32)>,
    flavor: BooleanFlavor,
) -> PyResult<(Vec<(i32, i32)>, Vec<Vec<(i32, i32)>>)> {
    with_current_model(move |model| {
        if objects.is_empty() {
            return Err(PyValueError::new_err("boolean object list is empty"));
        }
        let mut objects_iter = objects.iter();
        let (_, first_tag) = objects_iter.next().unwrap();
        let mut current = match model.shapes.get(first_tag) {
            Some(OccShape::BRep { brep, .. }) => brep.clone(),
            _ => {
                return Err(PyValueError::new_err(format!(
                    "boolean object entity {first_tag} is not a CAD shape"
                )))
            }
        };
        for (_, t) in objects_iter {
            let next = match model.shapes.get(t) {
                Some(OccShape::BRep { brep, .. }) => brep.clone(),
                _ => {
                    return Err(PyValueError::new_err(format!(
                        "boolean object entity {t} is not a CAD shape"
                    )))
                }
            };
            current = match flavor {
                // Within the object list, cut/intersect group with fuse.
                BooleanFlavor::Cut | BooleanFlavor::Fuse => rcad_algo::fuse(&current, &next),
                BooleanFlavor::Intersect => rcad_algo::common(&current, &next),
            }
            .map_err(|e| PyRuntimeError::new_err(format!("boolean failed: {e:?}")))?;
        }
        for (_, t) in &tools {
            let tool = match model.shapes.get(t) {
                Some(OccShape::BRep { brep, .. }) => brep.clone(),
                _ => {
                    return Err(PyValueError::new_err(format!(
                        "boolean tool entity {t} is not a CAD shape"
                    )))
                }
            };
            current = match flavor {
                BooleanFlavor::Cut => rcad_algo::cut(&current, &tool),
                BooleanFlavor::Fuse => rcad_algo::fuse(&current, &tool),
                BooleanFlavor::Intersect => rcad_algo::common(&current, &tool),
            }
            .map_err(|e| PyRuntimeError::new_err(format!("boolean failed: {e:?}")))?;
        }
        let tag = insert_brep(model, 3, current);
        let out = vec![(3, tag)];
        // v1: the per-input-output map is a single pass-through entry.
        let map = objects.iter().map(|&dt| vec![dt]).collect();
        Ok((out, map))
    })
}

#[pyfunction]
#[pyo3(signature = (object_dim_tags, tool_dim_tags, tag=-1, remove_object=true, remove_tool=true))]
pub fn model_occ_cut(
    object_dim_tags: Vec<(i32, i32)>,
    tool_dim_tags: Vec<(i32, i32)>,
    tag: i32,
    remove_object: bool,
    remove_tool: bool,
) -> PyResult<(Vec<(i32, i32)>, Vec<Vec<(i32, i32)>>)> {
    let _ = tag;
    let result = fold_boolean(
        object_dim_tags.clone(),
        tool_dim_tags.clone(),
        BooleanFlavor::Cut,
    );
    with_current_model(move |model| {
        let r = result?;
        if remove_object {
            for (_, t) in &object_dim_tags {
                model.shapes.remove(t);
            }
        }
        if remove_tool {
            for (_, t) in &tool_dim_tags {
                model.shapes.remove(t);
            }
        }
        Ok(r)
    })
}

#[pyfunction]
#[pyo3(signature = (object_dim_tags, tool_dim_tags, tag=-1, remove_object=true, remove_tool=true))]
pub fn model_occ_fuse(
    object_dim_tags: Vec<(i32, i32)>,
    tool_dim_tags: Vec<(i32, i32)>,
    tag: i32,
    remove_object: bool,
    remove_tool: bool,
) -> PyResult<(Vec<(i32, i32)>, Vec<Vec<(i32, i32)>>)> {
    let _ = tag;
    let result = fold_boolean(
        object_dim_tags.clone(),
        tool_dim_tags.clone(),
        BooleanFlavor::Fuse,
    );
    with_current_model(move |model| {
        let r = result?;
        if remove_object {
            for (_, t) in &object_dim_tags {
                model.shapes.remove(t);
            }
        }
        if remove_tool {
            for (_, t) in &tool_dim_tags {
                model.shapes.remove(t);
            }
        }
        Ok(r)
    })
}

#[pyfunction]
#[pyo3(signature = (object_dim_tags, tool_dim_tags, tag=-1, remove_object=true, remove_tool=true))]
pub fn model_occ_intersect(
    object_dim_tags: Vec<(i32, i32)>,
    tool_dim_tags: Vec<(i32, i32)>,
    tag: i32,
    remove_object: bool,
    remove_tool: bool,
) -> PyResult<(Vec<(i32, i32)>, Vec<Vec<(i32, i32)>>)> {
    let _ = tag;
    let result = fold_boolean(
        object_dim_tags.clone(),
        tool_dim_tags.clone(),
        BooleanFlavor::Intersect,
    );
    with_current_model(move |model| {
        let r = result?;
        if remove_object {
            for (_, t) in &object_dim_tags {
                model.shapes.remove(t);
            }
        }
        if remove_tool {
            for (_, t) in &tool_dim_tags {
                model.shapes.remove(t);
            }
        }
        Ok(r)
    })
}

#[pyfunction]
#[pyo3(signature = (object_dim_tags, tool_dim_tags, tag=-1, remove_object=true, remove_tool=true))]
pub fn model_occ_fragment(
    object_dim_tags: Vec<(i32, i32)>,
    tool_dim_tags: Vec<(i32, i32)>,
    tag: i32,
    remove_object: bool,
    remove_tool: bool,
) -> PyResult<(Vec<(i32, i32)>, Vec<Vec<(i32, i32)>>)> {
    let _ = (object_dim_tags, tool_dim_tags, tag, remove_object, remove_tool);
    Err(PyNotImplementedError::new_err(
        "occ.fragment needs the general n-ary boolean with output splitting - not wired in v1",
    ))
}

// ── transforms / copies ──────────────────────────────────────────────

#[pyfunction]
pub fn model_occ_translate(dim_tags: Vec<(i32, i32)>, dx: f64, dy: f64, dz: f64) -> PyResult<()> {
    with_current_model(|model| {
        for (_, t) in &dim_tags {
            match model.shapes.get_mut(t) {
                Some(OccShape::BRep { brep, .. }) => {
                    brep.apply_transform(DAffine3::from_translation(DVec3::new(dx, dy, dz)));
                }
                Some(_) => {}
                None => return Err(PyValueError::new_err(format!("entity {t} not found"))),
            }
        }
        Ok(())
    })
}

#[pyfunction]
#[pyo3(signature = (dim_tags, x, y, z, ax, ay, az, angle))]
pub fn model_occ_rotate(
    dim_tags: Vec<(i32, i32)>, x: f64, y: f64, z: f64, ax: f64, ay: f64, az: f64, angle: f64,
) -> PyResult<()> {
    with_current_model(|model| {
        let axis = DVec3::new(ax, ay, az);
        if axis.length() == 0.0 {
            return Err(PyValueError::new_err("rotation axis is zero"));
        }
        let pivot = DVec3::new(x, y, z);
        let rot = DAffine3::from_rotation_translation(
            DQuat::from_axis_angle(axis.normalize(), angle),
            pivot,
        ) * DAffine3::from_translation(-pivot);
        for (_, t) in &dim_tags {
            match model.shapes.get_mut(t) {
                Some(OccShape::BRep { brep, .. }) => brep.apply_transform(rot),
                Some(_) => {}
                None => return Err(PyValueError::new_err(format!("entity {t} not found"))),
            }
        }
        Ok(())
    })
}

#[pyfunction]
#[pyo3(signature = (dim_tags, x, y, z, a, b, c))]
pub fn model_occ_dilate(
    dim_tags: Vec<(i32, i32)>, x: f64, y: f64, z: f64, a: f64, b: f64, c: f64,
) -> PyResult<()> {
    with_current_model(|model| {
        let pivot = DVec3::new(x, y, z);
        let scale = DMat3::from_cols(
            DVec3::new(a, 0.0, 0.0),
            DVec3::new(0.0, b, 0.0),
            DVec3::new(0.0, 0.0, c),
        );
        let m = DAffine3::from_translation(pivot)
            * DAffine3::from_mat3(scale)
            * DAffine3::from_translation(-pivot);
        for (_, t) in &dim_tags {
            match model.shapes.get_mut(t) {
                Some(OccShape::BRep { brep, .. }) => brep.apply_transform(m),
                Some(_) => {}
                None => return Err(PyValueError::new_err(format!("entity {t} not found"))),
            }
        }
        Ok(())
    })
}

#[pyfunction]
#[pyo3(signature = (dim_tags, a, b, c, d))]
pub fn model_occ_mirror(dim_tags: Vec<(i32, i32)>, a: f64, b: f64, c: f64, d: f64) -> PyResult<()> {
    with_current_model(|model| {
        let n = DVec3::new(a, b, c);
        let len = n.length();
        if len == 0.0 {
            return Err(PyValueError::new_err("mirror plane normal is zero"));
        }
        let n = n / len;
        let d = d / len;
        // Reflection across the plane n.x + d = 0: x' = x - 2 (n.x + d) n.
        let nn = DMat3::from_cols(n * n.x, n * n.y, n * n.z);
        let m = DAffine3::from_mat3_translation(DMat3::IDENTITY - 2.0 * nn, -2.0 * d * n);
        for (_, t) in &dim_tags {
            match model.shapes.get_mut(t) {
                Some(OccShape::BRep { brep, .. }) => brep.apply_transform(m),
                Some(_) => {}
                None => return Err(PyValueError::new_err(format!("entity {t} not found"))),
            }
        }
        Ok(())
    })
}

#[pyfunction]
pub fn model_occ_copy(dim_tags: Vec<(i32, i32)>) -> PyResult<Vec<(i32, i32)>> {
    with_current_model(|model| {
        let mut out = Vec::new();
        for (d, t) in &dim_tags {
            let cloned = match model.shapes.get(t) {
                Some(OccShape::BRep { dim, brep }) => OccShape::BRep {
                    dim: *dim,
                    brep: brep.clone(),
                },
                Some(OccShape::Point(p)) => OccShape::Point(*p),
                Some(OccShape::Line(a, b)) => OccShape::Line(*a, *b),
                Some(OccShape::Arc(a, c, b)) => OccShape::Arc(*a, *c, *b),
                Some(OccShape::Spline(pts)) => OccShape::Spline(pts.clone()),
                Some(OccShape::Loop { dim, members }) => OccShape::Loop {
                    dim: *dim,
                    members: members.clone(),
                },
                None => return Err(PyValueError::new_err(format!("entity {t} not found"))),
            };
            let tag = model.alloc_tag();
            model.shapes.insert(tag, cloned);
            out.push((*d, tag));
        }
        Ok(out)
    })
}

#[pyfunction]
#[pyo3(signature = (dim_tags, recursive=false))]
pub fn model_occ_remove(dim_tags: Vec<(i32, i32)>, recursive: bool) -> PyResult<()> {
    let _ = recursive;
    with_current_model(|model| {
        for (_, t) in &dim_tags {
            model.shapes.remove(t);
        }
        Ok(())
    })
}

#[pyfunction]
#[pyo3(signature = (file_name, highest_dim_only=true, format=""))]
pub fn model_occ_import_shapes(
    file_name: String,
    highest_dim_only: bool,
    format: &str,
) -> PyResult<Vec<(i32, i32)>> {
    let _ = (highest_dim_only, format);
    with_current_model(|model| {
        let brep = rcad_step::StepReader::read_file(&file_name)
            .map_err(|e| PyValueError::new_err(format!("STEP import failed: {e}")))?;
        let tag = insert_brep(model, 3, brep);
        Ok(vec![(3, tag)])
    })
}

// ── properties ───────────────────────────────────────────────────────

#[pyfunction]
pub fn model_occ_get_mass(tag: i32) -> PyResult<f64> {
    with_current_model(|model| match model.shapes.get(&tag) {
        Some(OccShape::BRep { brep, .. }) => Ok(rcad_kernel::volume(brep)),
        Some(_) => Ok(0.0),
        None => Err(PyValueError::new_err(format!("entity {tag} not found"))),
    })
}

#[pyfunction]
pub fn model_occ_get_properties(tag: i32) -> PyResult<(f64, f64, f64, f64, f64)> {
    with_current_model(|model| match model.shapes.get(&tag) {
        Some(OccShape::BRep { brep, .. }) => {
            let vol = rcad_kernel::volume(brep);
            let area = rcad_kernel::surface_area(brep);
            let c = rcad_kernel::centroid(brep);
            Ok((vol, area, c.x, c.y, c.z))
        }
        Some(_) => Ok((0.0, 0.0, 0.0, 0.0, 0.0)),
        None => Err(PyValueError::new_err(format!("entity {tag} not found"))),
    })
}

// ── v1 honest omissions ──────────────────────────────────────────────

#[pyfunction]
#[pyo3(signature = (tag, curve_tags, radii))]
pub fn model_occ_fillet(tag: i32, curve_tags: Vec<i32>, radii: Vec<f64>) -> PyResult<i32> {
    let _ = (tag, curve_tags, radii);
    Err(PyNotImplementedError::new_err(
        "occ.fillet needs the fillet kernel (TKFillet port) - not available yet",
    ))
}

#[pyfunction]
#[pyo3(signature = (tag, curve_tags, distances))]
pub fn model_occ_chamfer(tag: i32, curve_tags: Vec<i32>, distances: Vec<f64>) -> PyResult<i32> {
    let _ = (tag, curve_tags, distances);
    Err(PyNotImplementedError::new_err(
        "occ.chamfer needs the fillet kernel (TKFillet port) - not available yet",
    ))
}

#[pyfunction]
#[pyo3(signature = (tag, face_idx, dx, dy, dz, distance))]
pub fn model_occ_extrude(
    tag: i32, face_idx: i32, dx: f64, dy: f64, dz: f64, distance: f64,
) -> PyResult<i32> {
    let _ = (tag, face_idx, dx, dy, dz, distance);
    Err(PyNotImplementedError::new_err(
        "occ.extrude from a solid's face requires face extraction + prism - not wired in v1; \
         build prisms from 2-D profiles instead",
    ))
}

#[pyfunction]
#[pyo3(signature = (tag, face_idx, ax, ay, az, dx, dy, dz, angle))]
pub fn model_occ_revolve(tag: i32, face_idx: i32, ax: f64, ay: f64, az: f64, dx: f64, dy: f64, dz: f64, angle: f64) -> PyResult<i32> {
    let _ = (tag, face_idx, ax, ay, az, dx, dy, dz, angle);
    Err(PyNotImplementedError::new_err(
        "occ.revolve from a solid's face is not wired in v1",
    ))
}

#[pyfunction]
#[pyo3(signature = (tag, tolerance=1e-8))]
pub fn model_occ_heal_shapes(tag: i32, tolerance: f64) -> PyResult<()> {
    let _ = (tag, tolerance);
    Err(PyNotImplementedError::new_err(
        "occ.healShapes is not wired in v1 (rcad shape healing runs inside STEP import)",
    ))
}

// ── synchronize ──────────────────────────────────────────────────────

/// Triangulate a B-Rep through the kernel's OCCT-port incremental mesher,
/// so the rmsh-bridge walker can read the face triangulations. Pattern
/// adapted from rcad-render's tessellate route (the pool is moved into the
/// shared lock and taken back — face TShapes are Arc-shared, so the mounts
/// stay consistent).
fn triangulate_brep(brep: &mut rcad_kernel::BRep, deflection: f64) -> bool {
    let Some(root) = top_level_root(brep) else {
        return false;
    };
    let mut params = rcad_mesh::i_mesh_tools::IMeshToolsParameters::default();
    params.deflection = deflection;
    params.deflection_interior = deflection;
    let shared: Arc<RwLock<rcad_kernel::BRep>> =
        Arc::new(RwLock::new(std::mem::take(brep)));
    // The mesher's worker threads can panic on retired/degenerate face
    // data (observed: the boolean result's retired sphere-half faces). A
    // Rust panic crossing the FFI boundary aborts the Python process, so
    // the perform is contained here - the rcad-render R23 pattern.
    let perform = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mesher = rcad_mesh::b_rep_mesh::facade::IncrementalMesh::perform(
            Arc::clone(&shared),
            &root,
            params,
        );
        let done = mesher.is_done();
        drop(mesher);
        done
    }));
    let done = perform.unwrap_or(false);
    // best-effort pool recovery even on panic (the lock may be poisoned)
    let mut meshed = match Arc::try_unwrap(shared) {
        Ok(lock) => lock
            .into_inner()
            .unwrap_or_else(|poisoned| poisoned.into_inner()),
        Err(shared) => shared
            .read()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .clone(),
    };
    std::mem::swap(brep, &mut meshed);
    done
}

/// The unreferenced top-level shapes of the pool (a single one, or a
/// compound of several) — the mesher's traversal root. Adapted from
/// rcad-render tessellate's top_level_root.
pub(crate) fn top_level_root(brep: &mut rcad_kernel::BRep) -> Option<rcad_kernel::topods::Shape> {
    use rcad_kernel::topods::Shape;
    if brep.tshapes.is_empty() {
        return None;
    }
    let mut referenced = vec![false; brep.tshapes.len()];
    for ts in &brep.tshapes {
        for_each_child_shape(ts, &mut |child: &Shape| {
            if child.index < referenced.len() {
                referenced[child.index] = true;
            }
        });
    }
    let tops: Vec<Shape> = (0..brep.tshapes.len())
        .filter(|&i| !referenced[i])
        .map(|i| brep.shape_at(i))
        .collect();
    match tops.len() {
        0 => None,
        1 => Some(tops.into_iter().next().expect("one top shape")),
        _ => Some(brep.add_tcompound(tops)),
    }
}

/// Visit the topological containment children of a TShape. Adapted from
/// rcad-render tessellate (the pool layout lives in the kernel; render
/// keeps its own walker local, and so does this binding).
fn for_each_child_shape(the_ts: &rcad_kernel::topods::TShape, the_f: &mut dyn FnMut(&rcad_kernel::topods::Shape)) {
    use rcad_kernel::topods::TShape;
    match the_ts {
        TShape::Vertex(_) => {}
        TShape::Edge(a_ed) => {
            the_f(&a_ed.first);
            the_f(&a_ed.last);
        }
        TShape::Wire(a_wd) => {
            for a_e in &a_wd.edges {
                the_f(a_e);
            }
        }
        TShape::Face(a_fd) => {
            the_f(&a_fd.outer_wire);
            for a_w in &a_fd.inner_wires {
                the_f(a_w);
            }
            for a_v in &a_fd.internal_vertices {
                the_f(a_v);
            }
        }
        TShape::Shell(a_shd) => {
            for a_face in &a_shd.faces {
                the_f(a_face);
            }
        }
        TShape::Solid(a_sd) => {
            for a_s in &a_sd.shells {
                the_f(a_s);
            }
            for a_v in &a_sd.internal_vertices {
                the_f(a_v);
            }
            for a_e in &a_sd.internal_edges {
                the_f(a_e);
            }
        }
        TShape::Compound(a_cd) => {
            for a_s in a_cd.iter() {
                the_f(a_s);
            }
        }
        TShape::CompSolid(a_cd) => {
            for a_s in a_cd.iter() {
                the_f(a_s);
            }
        }
    }
}

/// Synchronize: triangulate every CAD shape through the kernel mesher,
/// walk it through the rmsh-bridge B-Rep reader into the model's mesh
/// (welded Triangle3 per entity), and record the element -> entity map so
/// physical groups apply at write time.
#[pyfunction]
pub fn model_occ_synchronize() -> PyResult<()> {
    with_current_model(|model| {
        let mut mesh = rmsh_model::Mesh::new();
        let mut element_entity: BTreeMap<u64, i32> = BTreeMap::new();
        let deflection = model
            .mesh_size
            .map(|s| (s * 0.1).max(1.0e-4))
            .unwrap_or(0.01);
        for (&tag, shape) in &model.shapes {
            if let OccShape::BRep { brep, .. } = shape {
                let mut working = brep.clone();
                triangulate_brep(&mut working, deflection);
                let part = rmsh_bridge::brep_surface_mesh(&working);
                let offset = mesh.nodes.keys().copied().max().unwrap_or(0);
                let mut map: BTreeMap<u64, u64> = BTreeMap::new();
                for node in part.nodes.values() {
                    let new_id = node.id + offset;
                    map.insert(node.id, new_id);
                    mesh.add_node(rmsh_model::Node::new(
                        new_id,
                        node.position.x,
                        node.position.y,
                        node.position.z,
                    ));
                }
                let element_offset = mesh.elements.iter().map(|e| e.id).max().unwrap_or(0);
                for element in &part.elements {
                    let new_id = element.id + element_offset;
                    let node_ids: Vec<u64> = element.node_ids.iter().map(|n| map[n]).collect();
                    let el = rmsh_model::Element::new(new_id, element.etype, node_ids);
                    mesh.add_element(el);
                    element_entity.insert(new_id, tag);
                }
            }
        }
        model.mesh = mesh;
        model.element_entity = element_entity;
        Ok(())
    })
}

// keep shape_brep reachable for future cards (mesh entity queries);
// currently unused because all v1 paths clone the brep.
#[allow(dead_code)]
fn _keep_shape_brep(model: &mut crate::state::Model, tag: i32) -> PyResult<()> {
    shape_brep(model, tag).map(|_| ())
}
