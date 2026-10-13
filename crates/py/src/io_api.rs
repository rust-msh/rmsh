//! Root-level `initialize / finalize / clear / open / merge / write`.

use std::path::Path;

use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;

use crate::state::{self, lock, Model, OccShape};

fn load_msh(path: &Path) -> PyResult<rmsh_model::Mesh> {
    rmsh_io::load_msh_from_path(path).map_err(|e| PyIOError::new_err(e.to_string()))
}

#[pyfunction]
#[pyo3(signature = (argv0="", read_config_files=true, run_in_background=false))]
pub fn initialize_impl(argv0: &str, read_config_files: bool, run_in_background: bool) -> PyResult<()> {
    let _ = (argv0, read_config_files, run_in_background);
    let mut state = lock();
    if state.initialized {
        return Err(PyRuntimeError::new_err(
            "rmsh is already initialized; call finalize() first",
        ));
    }
    state.initialized = true;
    Ok(())
}

#[pyfunction]
pub fn finalize_impl() -> PyResult<()> {
    let mut state = lock();
    if !state.initialized {
        return Err(PyRuntimeError::new_err(
            "rmsh is not initialized; call initialize() first",
        ));
    }
    *state = state::RuntimeState::default();
    Ok(())
}

#[pyfunction]
pub fn clear_impl() -> PyResult<()> {
    let mut state = lock();
    state::require_initialized(&state)?;
    state.models.clear();
    state.current = None;
    Ok(())
}

/// `open` replaces the workspace with a new model loaded from the file
/// (gmsh semantics: a fresh model named after the file stem).
#[pyfunction]
pub fn open_impl(file_name: &str) -> PyResult<()> {
    let path = Path::new(file_name);
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();
    let stem = path
        .file_stem()
        .and_then(|s| Some(s.to_string_lossy().into_owned()))
        .unwrap_or_else(|| "model".to_string());

    let mut state = lock();
    state::require_initialized(&state)?;

    match ext.as_str() {
        "msh" => {
            let mesh = load_msh(path)?;
            let mut model = Model::new(stem.clone());
            model.mesh = mesh;
            state.models.insert(stem.clone(), model);
            state.current = Some(stem);
            Ok(())
        }
        "step" | "stp" => {
            let brep = rcad_step::StepReader::read_file(path)
                .map_err(|e| PyValueError::new_err(format!("STEP import failed: {e}")))?;
            let mut model = Model::new(stem.clone());
            let tag = model.alloc_tag();
            model.shapes.insert(tag, OccShape::BRep { dim: 3, brep });
            state.models.insert(stem.clone(), model);
            state.current = Some(stem);
            Ok(())
        }
        other => Err(PyValueError::new_err(format!(
            "unsupported file extension .{other}; expected .msh, .step, or .stp"
        ))),
    }
}

/// `merge` appends the file's content into the current model.
#[pyfunction]
pub fn merge_impl(file_name: &str) -> PyResult<()> {
    let path = Path::new(file_name);
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();

    let mut state = lock();
    state::require_initialized(&state)?;
    let name = state
        .current
        .clone()
        .ok_or_else(|| PyRuntimeError::new_err("no current model; call model.add() first"))?;

    match ext.as_str() {
        "msh" => {
            let mesh = load_msh(path)?;
            let model = state.models.get_mut(&name).unwrap();
            merge_mesh(model, mesh)?;
            Ok(())
        }
        "step" | "stp" => {
            let brep = rcad_step::StepReader::read_file(path)
                .map_err(|e| PyValueError::new_err(format!("STEP import failed: {e}")))?;
            let model = state.models.get_mut(&name).unwrap();
            let tag = model.alloc_tag();
            model.shapes.insert(tag, OccShape::BRep { dim: 3, brep });
            Ok(())
        }
        other => Err(PyValueError::new_err(format!(
            "unsupported file extension .{other}; expected .msh, .step, or .stp"
        ))),
    }
}

/// Append `patch`'s nodes/elements into `model`, renumbering to avoid id
/// collisions (nodes keep their coordinates; element ids shift past the
/// existing max).
fn merge_mesh(model: &mut Model, patch: rmsh_model::Mesh) -> PyResult<()> {
    let node_offset = model.mesh.nodes.keys().copied().max().unwrap_or(0);
    let mut map = std::collections::BTreeMap::new();
    for node in patch.nodes.values() {
        let new_id = node.id + node_offset;
        map.insert(node.id, new_id);
        model.mesh.add_node(rmsh_model::Node::new(
            new_id,
            node.position.x,
            node.position.y,
            node.position.z,
        ));
    }
    let element_offset = model.mesh.elements.iter().map(|e| e.id).max().unwrap_or(0);
    for element in &patch.elements {
        let node_ids: Vec<u64> = element.node_ids.iter().map(|n| map[n]).collect();
        let el = rmsh_model::Element::new(element.id + element_offset, element.etype, node_ids);
        model.mesh.add_element(el);
    }
    for ((dim, tag), name) in &patch.physical_names {
        model
            .mesh
            .physical_names
            .entry((*dim, *tag))
            .or_insert_with(|| name.clone());
    }
    Ok(())
}

#[pyfunction]
pub fn write_impl(file_name: &str) -> PyResult<()> {
    let path = Path::new(file_name);
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();

    let mut state = lock();
    state::require_initialized(&state)?;
    let name = state
        .current
        .clone()
        .ok_or_else(|| PyRuntimeError::new_err("no current model; call model.add() first"))?;
    let version = state::option_number(&state.options, "Mesh.MshFileVersion").unwrap_or(2.2);
    let binary = state::option_number(&state.options, "Mesh.Binary").unwrap_or(0.0) >= 0.5;
    let model = state.models.get_mut(&name).unwrap();

    match ext.as_str() {
        "msh" => {
            apply_physicals_to_mesh(model);
            let mut file = std::fs::File::create(path)
                .map_err(|e| PyIOError::new_err(e.to_string()))?;
            let result = match (version, binary) {
                (v, _) if v < 4.0 => {
                    if binary {
                        rmsh_io::write_msh_v2_binary(&mut file, &model.mesh)
                    } else {
                        rmsh_io::write_msh_v2(&mut file, &model.mesh)
                    }
                }
                _ => {
                    if binary {
                        rmsh_io::write_msh_v4_binary(&mut file, &model.mesh)
                    } else {
                        rmsh_io::write_msh_v4(&mut file, &model.mesh)
                    }
                }
            };
            result.map_err(|e| PyIOError::new_err(e.to_string()))
        }
        "step" | "stp" => {
            // v1: exactly one CAD solid is exportable (a multi-solid model
            // needs an assembly writer).
            let solids: Vec<&rcad_kernel::BRep> = model
                .shapes
                .values()
                .filter_map(|s| match s {
                    OccShape::BRep { dim: 3, brep } => Some(brep),
                    _ => None,
                })
                .collect();
            if solids.len() != 1 {
                return Err(PyValueError::new_err(format!(
                    "STEP export v1 writes exactly one solid; the model has {}",
                    solids.len()
                )));
            }
            let selection = rcad_step::ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            };
            let content = rcad_step::StepWriter::write_string(solids[0], selection);
            std::fs::write(path, content).map_err(|e| PyIOError::new_err(e.to_string()))
        }
        other => Err(PyValueError::new_err(format!(
            "unsupported file extension .{other}; expected .msh, .step, or .stp"
        ))),
    }
}

/// Physical groups: elements tagged with their source-entity tag at
/// synchronize time are rewritten to the FIRST matching group's tag (and
/// the group's name lands in the msh PhysicalNames table); unclaimed
/// entities strip their tag.
fn apply_physicals_to_mesh(model: &mut Model) {
    let mut claimed: std::collections::BTreeMap<u64, i32> = Default::default();
    for group in &model.physicals {
        for element in &model.mesh.elements {
            if element.dimension() as i32 != group.dim {
                continue;
            }
            if let Some(&entity) = model.element_entity.get(&element.id) {
                if group.entities.contains(&entity) && !claimed.contains_key(&element.id) {
                    claimed.insert(element.id, group.tag);
                }
            }
        }
    }
    for element in &mut model.mesh.elements {
        element.physical_tag = claimed.get(&element.id).copied();
    }
    // Wire-frame entities that carry no elements still contribute their
    // group through the PhysicalNames table via setPhysicalName.
}
