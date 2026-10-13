//! `model.mesh.*` bindings: mesh generation over the synchronized CAD
//! surface (rmsh-bridge output) and the gmsh-shaped data access surface.
//! Element type numbers are the gmsh MSH type ids natively - the P5
//! semantic variants map through `ElementType::to_gmsh_type_id`.

use pyo3::exceptions::{PyNotImplementedError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;

use rmsh_algo::traits::{Mesher3D, MeshParams};

use crate::state::option_number;

fn mesh_params(model: &crate::state::Model, options: &std::collections::BTreeMap<String, crate::state::OptionValue>) -> MeshParams {
    let size = model
        .mesh_size
        .or_else(|| option_number(options, "Mesh.MeshSizeMax"))
        .unwrap_or(0.25);
    let mut params = MeshParams::with_size(size);
    if let Some(min) = option_number(options, "Mesh.MeshSizeMin") {
        if min > 0.0 {
            params.min_size = min;
        }
    }
    params
}

#[pyfunction]
#[pyo3(signature = (dim_tags, size))]
pub fn model_mesh_set_size(dim_tags: Vec<(i32, i32)>, size: f64) -> PyResult<()> {
    if size <= 0.0 {
        return Err(PyValueError::new_err("mesh size must be positive"));
    }
    let _ = &dim_tags; // v1: the size applies model-wide (documented)
    crate::state::with_current_model(|model| {
        model.mesh_size = Some(size);
        Ok(())
    })
}

#[pyfunction]
#[pyo3(signature = (dim))]
pub fn model_mesh_generate(dim: i32) -> PyResult<()> {
    crate::state::with_current_model_options(|options, model| {
        match dim {
            2 => {
                // The synchronized surface mesh IS the 2-D mesh; regenerate
                // it from CAD when the caller changed nothing (idempotent).
                if model.mesh.nodes.is_empty() {
                    return Err(pyo3::exceptions::PyRuntimeError::new_err(
                        "empty mesh; call model.occ.synchronize() after building CAD entities",
                    ));
                }
                Ok(())
            }
            3 => {
                let params = mesh_params(model, options);
                let volume = rmsh_algo::CentroidStarMesher3D
                    .mesh_3d(&model.mesh, &params)
                    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
                // Keep the surface elements: re-add the 2-D elements of the
                // synchronized mesh on top of the volume mesh (boundary node
                // ids are preserved by the centroid-star mesher).
                let mut combined = volume;
                for element in &model.mesh.elements {
                    if element.dimension() == 2 {
                        combined.add_element(element.clone());
                    }
                }
                model.mesh = combined;
                Ok(())
            }
            other => Err(PyValueError::new_err(format!(
                "mesh.generate supports dim 2 or 3, got {other}"
            ))),
        }
    })
}

#[pyfunction]
#[pyo3(signature = (dim=-1, tag=-1, include_boundary=false, return_parametric_coord=false))]
pub fn model_mesh_get_nodes(
    dim: i32,
    tag: i32,
    include_boundary: bool,
    return_parametric_coord: bool,
) -> PyResult<(Vec<u64>, Vec<f64>, Vec<f64>)> {
    let _ = (tag, include_boundary, return_parametric_coord);
    let _ = dim; // v1: nodes are reported model-wide
    crate::state::with_current_model_ref(|model| {
        let mut tags: Vec<u64> = model.mesh.nodes.keys().copied().collect();
        tags.sort_unstable();
        let mut coords = Vec::with_capacity(tags.len() * 3);
        for t in &tags {
            let node = &model.mesh.nodes[t];
            coords.push(node.position.x);
            coords.push(node.position.y);
            coords.push(node.position.z);
        }
        let parametric = vec![0.0; coords.len()];
        Ok((tags, coords, parametric))
    })
}

#[pyfunction]
#[pyo3(signature = (dim=-1, tag=-1))]
pub fn model_mesh_get_elements(
    dim: i32,
    tag: i32,
) -> PyResult<(Vec<i32>, Vec<Vec<u64>>, Vec<Vec<u64>>)> {
    let _ = tag;
    crate::state::with_current_model_ref(|model| {
        // group by gmsh type id, preserving ascending element id order
        let mut by_type: std::collections::BTreeMap<i32, Vec<&rmsh_model::Element>> =
            std::collections::BTreeMap::new();
        for element in model.mesh.elements.iter() {
            if dim >= 0 && (element.dimension() as i32) != dim {
                continue;
            }
            let type_id = element
                .etype
                .to_gmsh_type_id()
                .unwrap_or_else(|| match element.etype {
                    rmsh_model::ElementType::Unknown(id) => id,
                    other => panic!("element {other:?} lacks a gmsh type id"),
                });
            by_type.entry(type_id).or_default().push(element);
        }
        let mut types = Vec::new();
        let mut tags = Vec::new();
        let mut node_tags = Vec::new();
        for (type_id, elements) in by_type {
            types.push(type_id);
            tags.push(elements.iter().map(|e| e.id).collect());
            node_tags.push(
                elements
                    .iter()
                    .flat_map(|e| e.node_ids.iter().copied())
                    .collect(),
            );
        }
        Ok((types, tags, node_tags))
    })
}

#[pyfunction]
#[pyo3(signature = (order))]
pub fn model_mesh_set_order(order: i32) -> PyResult<()> {
    match order {
        1 => Ok(()),
        n => Err(PyNotImplementedError::new_err(format!(
            "mesh.setOrder({n}) needs the geometry-projection lift of midside nodes - \
             the P5 curved-element machinery is in place, the converter is a follow-up"
        ))),
    }
}

#[pyfunction]
#[pyo3(signature = (method="", force=false, niter=1, dim_tags=Vec::new(), dim_tags_n=0))]
pub fn model_mesh_optimize(
    method: &str,
    force: bool,
    niter: i32,
    dim_tags: Vec<(i32, i32)>,
    dim_tags_n: i32,
) -> PyResult<()> {
    let _ = (method, force, niter, dim_tags, dim_tags_n);
    let _ = niter;
    Err(PyNotImplementedError::new_err(
        "mesh.optimize is not wired in v1 (the rmsh quality optimizer is exposed through \
         the CLI quality verb; the binding lands with the quality-stat surface)",
    ))
}

#[pyfunction]
#[pyo3(signature = (dim=2, tag=-1, angle=45.0))]
pub fn model_mesh_recombine(dim: i32, tag: i32, angle: f64) -> PyResult<()> {
    let _ = (dim, tag, angle);
    Err(PyNotImplementedError::new_err(
        "mesh.recombine over the synchronized surface is a follow-up (the P3b quad \
         recombination kernel exists; the binding is not wired)",
    ))
}

#[pyfunction]
pub fn model_mesh_refine() -> PyResult<()> {
    Err(PyNotImplementedError::new_err(
        "mesh.refine (uniform h-refinement) is not wired in v1",
    ))
}

#[pyfunction]
#[pyo3(signature = (dim_tags=Vec::new(), dim_tags_n=0))]
pub fn model_mesh_clear(dim_tags: Vec<(i32, i32)>, dim_tags_n: i32) -> PyResult<()> {
    let _ = (dim_tags, dim_tags_n);
    crate::state::with_current_model(|model| {
        model.mesh = rmsh_model::Mesh::new();
        model.element_entity.clear();
        Ok(())
    })
}
