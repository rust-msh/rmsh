//! `model.*`, `option.*` and `logger.*` bindings.

use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;

use crate::state::{
    self, lock, with_current_model, with_current_model_ref, OptionValue,
};

// ── model lifecycle ──────────────────────────────────────────────────

#[pyfunction]
pub fn model_add(name: &str) -> PyResult<()> {
    let mut state = lock();
    state::require_initialized(&state)?;
    if state.models.contains_key(name) {
        return Err(PyValueError::new_err(format!(
            "model {name:?} already exists"
        )));
    }
    state.models.insert(name.to_string(), state::Model::new(name.to_string()));
    state.current = Some(name.to_string());
    Ok(())
}

#[pyfunction]
pub fn model_remove() -> PyResult<()> {
    let mut state = lock();
    state::require_initialized(&state)?;
    if let Some(name) = state.current.take() {
        state.models.remove(&name);
    }
    Ok(())
}

#[pyfunction]
pub fn model_get_current() -> PyResult<String> {
    let state = lock();
    state::require_initialized(&state)?;
    Ok(state.current.clone().unwrap_or_default())
}

#[pyfunction]
pub fn model_set_current(name: &str) -> PyResult<()> {
    let mut state = lock();
    state::require_initialized(&state)?;
    if !state.models.contains_key(name) {
        return Err(PyValueError::new_err(format!(
            "model {name:?} does not exist"
        )));
    }
    state.current = Some(name.to_string());
    Ok(())
}

#[pyfunction]
pub fn model_get_dimension() -> PyResult<i32> {
    with_current_model_ref(|model| {
        Ok(model
            .shapes
            .values()
            .map(|s| s.dimension())
            .chain(model.mesh.elements.iter().map(|e| e.dimension() as i32))
            .max()
            .unwrap_or(0))
    })
}

#[pyfunction]
#[pyo3(signature = (dim=-1))]
pub fn model_get_entities(dim: i32) -> PyResult<Vec<(i32, i32)>> {
    with_current_model_ref(|model| {
        let mut out: Vec<(i32, i32)> = model
            .shapes
            .iter()
            .map(|(tag, shape)| (shape.dimension(), *tag))
            .filter(|(d, _)| dim < 0 || *d == dim)
            .collect();
        out.sort();
        Ok(out)
    })
}

#[pyfunction]
#[pyo3(signature = (xmin, ymin, zmin, xmax, ymax, zmax, dim=-1))]
pub fn model_get_entities_in_bounding_box(
    xmin: f64, ymin: f64, zmin: f64, xmax: f64, ymax: f64, zmax: f64, dim: i32,
) -> PyResult<Vec<(i32, i32)>> {
    let _ = (xmin, ymin, zmin, xmax, ymax, zmax, dim);
    Err(PyNotImplementedError::new_err(
        "model.getEntitiesInBoundingBox is not wired in v1 (use getEntities)",
    ))
}

#[pyfunction]
#[pyo3(signature = (dim=-1, tag=-1))]
pub fn model_get_bounding_box(dim: i32, tag: i32) -> PyResult<(f64, f64, f64, f64, f64, f64)> {
    with_current_model_ref(|model| {
        if dim < 0 {
            let (min, max) = model.bbox();
            return Ok((min[0], min[1], min[2], max[0], max[1], max[2]));
        }
        match model.shapes.get(&tag) {
            Some(shape) => {
                let (min, max) = match shape {
                    state::OccShape::BRep { brep, .. } => {
                        let bb = brep
                            .bounding_box()
                            .ok_or_else(|| {
                                PyValueError::new_err("empty shape has no bounding box")
                            })?;
                        (
                            [bb[0].x, bb[0].y, bb[0].z],
                            [bb[1].x, bb[1].y, bb[1].z],
                        )
                    }
                    state::OccShape::Point(p) => ([p.x, p.y, p.z], [p.x, p.y, p.z]),
                    _ => {
                        return Err(PyValueError::new_err(
                            "bounding box of wire-frame entities is not available",
                        ))
                    }
                };
                Ok((min[0], min[1], min[2], max[0], max[1], max[2]))
            }
            None => Err(PyValueError::new_err(format!(
                "entity (dim={dim}, tag={tag}) not found"
            ))),
        }
    })
}

#[pyfunction]
#[pyo3(signature = (dim, tag))]
pub fn model_get_entity_name(dim: i32, tag: i32) -> PyResult<String> {
    with_current_model_ref(|model| {
        for group in &model.physicals {
            if group.dim == dim && group.entities.contains(&tag) {
                if let Some(name) = &group.name {
                    return Ok(name.clone());
                }
            }
        }
        Ok(String::new())
    })
}

#[pyfunction]
pub fn model_set_entity_name(_dim: i32, _tag: i32, _name: &str) -> PyResult<()> {
    Err(PyNotImplementedError::new_err(
        "model.setEntityName is not wired in v1 (use setPhysicalName)",
    ))
}

// ── physical groups ──────────────────────────────────────────────────

#[pyfunction]
#[pyo3(signature = (dim, tags, tag=-1, name=""))]
pub fn model_add_physical_group(dim: i32, tags: Vec<i32>, tag: i32, name: &str) -> PyResult<i32> {
    with_current_model(|model| {
        let tag = match tag {
            -1 => {
                let t = model.next_physical_tag;
                model.next_physical_tag += 1;
                t
            }
            t => {
                if t >= model.next_physical_tag {
                    model.next_physical_tag = t + 1;
                }
                t
            }
        };
        model.physicals.push(state::PhysicalGroup {
            dim,
            tag,
            entities: tags,
            name: if name.is_empty() { None } else { Some(name.to_string()) },
        });
        Ok(tag)
    })
}

#[pyfunction]
#[pyo3(signature = (dim=-1))]
pub fn model_get_physical_groups(dim: i32) -> PyResult<Vec<(i32, i32)>> {
    with_current_model_ref(|model| {
        Ok(model
            .physicals
            .iter()
            .filter(|g| dim < 0 || g.dim == dim)
            .map(|g| (g.dim, g.tag))
            .collect())
    })
}

#[pyfunction]
#[pyo3(signature = (dim, tag, name))]
pub fn model_set_physical_name(dim: i32, tag: i32, name: &str) -> PyResult<()> {
    with_current_model(|model| {
        let group = model
            .physicals
            .iter_mut()
            .find(|g| g.dim == dim && g.tag == tag)
            .ok_or_else(|| {
                PyValueError::new_err(format!("physical group ({dim}, {tag}) not found"))
            })?;
        group.name = Some(name.to_string());
        model.mesh.physical_names.insert(
            (dim, tag),
            name.to_string(),
        );
        Ok(())
    })
}

#[pyfunction]
#[pyo3(signature = (dim, tag))]
pub fn model_get_physical_name(dim: i32, tag: i32) -> PyResult<String> {
    with_current_model_ref(|model| {
        Ok(model
            .physicals
            .iter()
            .find(|g| g.dim == dim && g.tag == tag)
            .and_then(|g| g.name.clone())
            .unwrap_or_default())
    })
}

// ── options ──────────────────────────────────────────────────────────

#[pyfunction]
pub fn option_set_number(name: &str, value: f64) -> PyResult<()> {
    state::check_option_name(name)?;
    let mut state = lock();
    state
        .options
        .insert(name.to_string(), OptionValue::Number(value));
    Ok(())
}

#[pyfunction]
pub fn option_get_number(name: &str) -> PyResult<f64> {
    state::check_option_name(name)?;
    let state = lock();
    state::option_number(&state.options, name).ok_or_else(|| {
        PyValueError::new_err(format!("option {name:?} is not a number option"))
    })
}

#[pyfunction]
pub fn option_set_string(name: &str, value: &str) -> PyResult<()> {
    state::check_option_name(name)?;
    let mut state = lock();
    state
        .options
        .insert(name.to_string(), OptionValue::String(value.to_string()));
    Ok(())
}

#[pyfunction]
pub fn option_get_string(name: &str) -> PyResult<String> {
    state::check_option_name(name)?;
    let state = lock();
    match state.options.get(name) {
        Some(OptionValue::String(s)) => Ok(s.clone()),
        _ => match state::default_option(name) {
            Some(OptionValue::String(s)) => Ok(s),
            _ => Ok(String::new()),
        },
    }
}

#[pyfunction]
pub fn option_set_color(name: &str, r: i32, g: i32, b: i32, a: i32) -> PyResult<()> {
    state::check_option_name(name)?;
    let mut state = lock();
    state
        .options
        .insert(name.to_string(), OptionValue::Color((r, g, b, a)));
    Ok(())
}

#[pyfunction]
pub fn option_get_color(name: &str) -> PyResult<(i32, i32, i32, i32)> {
    state::check_option_name(name)?;
    let state = lock();
    match state.options.get(name) {
        Some(OptionValue::Color(c)) => Ok(*c),
        _ => Ok((255, 255, 255, 255)),
    }
}

#[pyfunction]
pub fn option_restore_defaults() -> PyResult<()> {
    let mut state = lock();
    state.options.clear();
    Ok(())
}

// ── logger ───────────────────────────────────────────────────────────

#[pyfunction]
pub fn logger_start() -> PyResult<()> {
    let mut state = lock();
    state::require_initialized(&state)?;
    state.logger_on = true;
    state.logs.clear();
    Ok(())
}

#[pyfunction]
pub fn logger_stop() -> PyResult<()> {
    let mut state = lock();
    state.logger_on = false;
    Ok(())
}

#[pyfunction]
pub fn logger_get() -> PyResult<Vec<String>> {
    let state = lock();
    Ok(state.logs.clone())
}
