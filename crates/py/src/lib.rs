use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::PathBuf;
use std::sync::{LazyLock, Mutex};
use std::time::Duration;

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_algorithms::{
    BooleanOpType, boolean_op, geom_populate,
    brep_repair::repair,
};
use rcad_modeling::builder::{
    cone_brep, torus_brep,
    fillet::{chamfer_edge, fillet_edges},
    ops::{extrude, revolve},
};
use rmsh_algo::{
    CentroidStarMesher3D, Delaunay3D, FrontalDelaunay2D, Frontal3D, Hxt3D,
    Bamg2D, QuadPaving2D,
    LaplacianSmooth,
    MeshAlgoError, MeshOptimizer, MeshParams, Mesher2D, Mesher3D,
    OptimizeParams, Polygon2D, mesh_polygon,
    Domain2D,
};
use rmsh_model::{Element, ElementType, Mesh, Node};

macro_rules! stub_pyfunction {
    ($rust_name:ident, $py_name:literal, $prototype:literal) => {
        #[pyfunction]
        #[pyo3(name = $py_name, signature = (*args, **kwargs))]
        fn $rust_name(
            args: &pyo3::Bound<'_, PyTuple>,
            kwargs: Option<&pyo3::Bound<'_, PyDict>>,
        ) -> pyo3::PyResult<()> {
            let _ = (args, kwargs);
            Err(PyNotImplementedError::new_err(format!(
                "{} is not implemented yet",
                $prototype
            )))
        }
    };
}

#[derive(Default)]
struct RuntimeState {
    initialized: bool,
    current_mesh: Option<Mesh>,
    current_path: Option<PathBuf>,
    option_numbers: HashMap<String, f64>,
    option_strings: HashMap<String, String>,
    option_colors: HashMap<String, (i32, i32, i32, i32)>,
    /// CAD shapes created via model.occ.add* functions, keyed by tag.
    cad_shapes: HashMap<i32, BRep>,
    /// Next auto-assigned tag for CAD shapes.
    next_cad_tag: i32,
}

static STATE: LazyLock<Mutex<RuntimeState>> = LazyLock::new(|| Mutex::new(RuntimeState::default()));

fn ensure_initialized(state: &RuntimeState) -> PyResult<()> {
    if state.initialized {
        Ok(())
    } else {
        Err(pyo3::exceptions::PyRuntimeError::new_err(
            "rmsh is not initialized; call initialize() first",
        ))
    }
}

fn load_mesh_from_path(path: &PathBuf) -> PyResult<Mesh> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase());

    match ext.as_deref() {
        Some("msh") => rmsh_io::load_msh_from_path(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string())),
        Some("step") | Some("stp") => rmsh_io::load_step_from_path(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string())),
        _ => Err(pyo3::exceptions::PyValueError::new_err(
            "unsupported file extension; expected .msh, .step, or .stp",
        )),
    }
}

fn boundary_loop_from_surface_mesh(mesh: &Mesh) -> PyResult<Vec<[f64; 2]>> {
    let mut edge_count: HashMap<(u64, u64), usize> = HashMap::new();
    for e in &mesh.elements {
        if e.dimension() != 2 || e.node_ids.len() < 3 {
            continue;
        }
        for i in 0..e.node_ids.len() {
            let a = e.node_ids[i];
            let b = e.node_ids[(i + 1) % e.node_ids.len()];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }

    let boundary_edges: Vec<(u64, u64)> = edge_count
        .into_iter()
        .filter_map(|(edge, count)| if count == 1 { Some(edge) } else { None })
        .collect();

    if boundary_edges.len() < 3 {
        return Err(pyo3::exceptions::PyRuntimeError::new_err(
            "cannot extract boundary loop from current 2D mesh",
        ));
    }

    let mut adjacency: HashMap<u64, Vec<u64>> = HashMap::new();
    for (a, b) in &boundary_edges {
        adjacency.entry(*a).or_default().push(*b);
        adjacency.entry(*b).or_default().push(*a);
    }

    let start = *adjacency.keys().min().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("cannot find boundary start node")
    })?;

    let mut loop_ids = vec![start];
    let mut visited_edges: HashSet<(u64, u64)> = HashSet::new();
    let mut prev: Option<u64> = None;
    let mut current = start;

    for _ in 0..(boundary_edges.len() + 2) {
        let neighbors = adjacency.get(&current).ok_or_else(|| {
            pyo3::exceptions::PyRuntimeError::new_err("broken boundary adjacency")
        })?;
        let next = neighbors
            .iter()
            .copied()
            .find(|n| Some(*n) != prev)
            .or_else(|| neighbors.first().copied())
            .ok_or_else(|| {
                pyo3::exceptions::PyRuntimeError::new_err("boundary traversal failed")
            })?;

        let edge_key = if current < next {
            (current, next)
        } else {
            (next, current)
        };
        if visited_edges.contains(&edge_key) {
            break;
        }
        visited_edges.insert(edge_key);

        current = next;
        if current == start {
            break;
        }
        loop_ids.push(current);
        prev = loop_ids.iter().rev().nth(1).copied();
    }

    if loop_ids.len() < 3 {
        return Err(pyo3::exceptions::PyRuntimeError::new_err(
            "extracted boundary loop is degenerate",
        ));
    }

    let mut polygon = Vec::with_capacity(loop_ids.len());
    for nid in loop_ids {
        let node = mesh.nodes.get(&nid).ok_or_else(|| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("missing node id {}", nid))
        })?;
        polygon.push([node.position.x, node.position.y]);
    }
    Ok(polygon)
}

fn merge_meshes(base: &mut Mesh, incoming: &Mesh) {
    let mut next_node_id = base.nodes.keys().max().copied().unwrap_or(0) + 1;
    let mut next_elem_id = base.elements.iter().map(|e| e.id).max().unwrap_or(0) + 1;

    let mut node_remap: HashMap<u64, u64> = HashMap::new();
    for node in incoming.nodes.values() {
        let new_id = next_node_id;
        next_node_id += 1;
        node_remap.insert(node.id, new_id);
        base.add_node(Node {
            id: new_id,
            position: node.position,
        });
    }

    for elem in &incoming.elements {
        let new_nodes: Vec<u64> = elem
            .node_ids
            .iter()
            .filter_map(|nid| node_remap.get(nid).copied())
            .collect();
        if new_nodes.len() != elem.node_ids.len() {
            continue;
        }
        let mut new_elem = Element::new(next_elem_id, elem.etype, new_nodes);
        new_elem.physical_tag = elem.physical_tag;
        base.add_element(new_elem);
        next_elem_id += 1;
    }
}

/// Convert a `BRep` into a `Mesh` by extracting its triangles.
fn tessellate_brep(brep: &BRep) -> Mesh {
    let mut mesh = Mesh::new();
    let mut node_id: u64 = 1;
    let mut elem_id: u64 = 1;

    // Map BRep vertex index → mesh node id
    let mut vertex_to_node: HashMap<usize, u64> = HashMap::new();

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for &[i0, i1, i2] in &face.triangles {
                    let mut nids = Vec::with_capacity(3);
                    for vi in [i0, i1, i2] {
                        let nid = *vertex_to_node.entry(vi).or_insert_with(|| {
                            let id = node_id;
                            node_id += 1;
                            if let Some(v) = brep.vertices.get(vi) {
                                mesh.add_node(Node::new(id, v.point.x, v.point.y, v.point.z));
                            }
                            id
                        });
                        nids.push(nid);
                    }
                    mesh.add_element(Element::new(elem_id, ElementType::Triangle3, nids));
                    elem_id += 1;
                }
            }
        }
    }
    mesh
}

fn extract_required<T>(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
    index: usize,
    kw_names: &[&str],
    expected: &str,
) -> PyResult<T>
where
    T: for<'a> FromPyObject<'a>,
{
    if let Some(kwargs) = kwargs {
        for name in kw_names {
            if let Some(value) = kwargs.get_item(name)? {
                return value.extract::<T>();
            }
        }
    }
    if index < args.len() {
        return args.get_item(index)?.extract::<T>();
    }
    Err(pyo3::exceptions::PyTypeError::new_err(format!(
        "missing required argument '{}': expected {}",
        kw_names.first().copied().unwrap_or("arg"),
        expected
    )))
}

fn option_number(state: &RuntimeState, name: &str) -> Option<f64> {
    state.option_numbers.get(name).copied()
}

fn gmsh_type_id(etype: ElementType) -> i32 {
    match etype {
        ElementType::Point1 => 15,
        ElementType::Line2 => 1,
        ElementType::Triangle3 => 2,
        ElementType::Quad4 => 3,
        ElementType::Tetrahedron4 => 4,
        ElementType::Hexahedron8 => 5,
        ElementType::Prism6 => 6,
        ElementType::Pyramid5 => 7,
        ElementType::Unknown(id) => id,
    }
}

fn mesh_max_dimension(mesh: &Mesh) -> i32 {
    mesh.elements
        .iter()
        .map(|e| i32::from(e.dimension()))
        .max()
        .unwrap_or(0)
}

#[pyfunction]
#[pyo3(name = "initialize", signature = (*args, **kwargs))]
fn initialize_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let _ = (args, kwargs);
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    state.initialized = true;
    state.current_mesh = None;
    state.current_path = None;
    Ok(())
}

#[pyfunction]
#[pyo3(name = "finalize", signature = (*args, **kwargs))]
fn finalize_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let _ = (args, kwargs);
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    state.initialized = false;
    state.current_mesh = None;
    state.current_path = None;
    state.option_numbers.clear();
    state.option_strings.clear();
    state.option_colors.clear();
    state.cad_shapes.clear();
    state.next_cad_tag = 0;
    Ok(())
}

#[pyfunction]
#[pyo3(name = "clear", signature = (*args, **kwargs))]
fn clear_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let _ = (args, kwargs);
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.current_mesh = None;
    state.current_path = None;
    state.cad_shapes.clear();
    state.next_cad_tag = 0;
    Ok(())
}

#[pyfunction]
#[pyo3(name = "open", signature = (*args, **kwargs))]
fn open_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let file_name: String = extract_required(args, kwargs, 0, &["fileName", "file_name"], "str")?;
    let path = PathBuf::from(&file_name);
    let mesh = load_mesh_from_path(&path)?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.current_mesh = Some(mesh);
    state.current_path = Some(path);
    Ok(())
}

#[pyfunction]
#[pyo3(name = "merge", signature = (*args, **kwargs))]
fn merge_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let file_name: String = extract_required(args, kwargs, 0, &["fileName", "file_name"], "str")?;
    let path = PathBuf::from(&file_name);
    let incoming = load_mesh_from_path(&path)?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    match state.current_mesh.as_mut() {
        Some(current) => merge_meshes(current, &incoming),
        None => state.current_mesh = Some(incoming),
    }
    Ok(())
}

#[pyfunction]
#[pyo3(name = "write", signature = (*args, **kwargs))]
fn write_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let file_name: String = extract_required(args, kwargs, 0, &["fileName", "file_name"], "str")?;
    let path = PathBuf::from(&file_name);
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase());

    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let mesh = state.current_mesh.as_ref().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded; call open() or generate() first")
    })?;

    match ext.as_deref() {
        Some("msh") => rmsh_io::save_msh_v4_to_path(&path, mesh)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?,
        Some("step") | Some("stp") => rmsh_io::save_step_to_path(&path, mesh)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?,
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "unsupported write format; only .msh and .step/.stp are currently supported",
            ));
        }
    }

    Ok(())
}

#[pyfunction]
#[pyo3(name = "option_set_number", signature = (*args, **kwargs))]
fn option_set_number_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let value: f64 = extract_required(args, kwargs, 1, &["value"], "float")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    state.option_numbers.insert(name, value);
    Ok(())
}

#[pyfunction]
#[pyo3(name = "option_get_number", signature = (*args, **kwargs))]
fn option_get_number_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<f64> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    state.option_numbers.get(&name).copied().ok_or_else(|| {
        pyo3::exceptions::PyKeyError::new_err(format!("number option not set: {}", name))
    })
}

#[pyfunction]
#[pyo3(name = "option_set_string", signature = (*args, **kwargs))]
fn option_set_string_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let value: String = extract_required(args, kwargs, 1, &["value"], "str")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    state.option_strings.insert(name, value);
    Ok(())
}

#[pyfunction]
#[pyo3(name = "option_get_string", signature = (*args, **kwargs))]
fn option_get_string_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    state.option_strings.get(&name).cloned().ok_or_else(|| {
        pyo3::exceptions::PyKeyError::new_err(format!("string option not set: {}", name))
    })
}

#[pyfunction]
#[pyo3(name = "option_set_color", signature = (*args, **kwargs))]
fn option_set_color_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let r: i32 = extract_required(args, kwargs, 1, &["r"], "int")?;
    let g: i32 = extract_required(args, kwargs, 2, &["g"], "int")?;
    let b: i32 = extract_required(args, kwargs, 3, &["b"], "int")?;
    let a: i32 = extract_required(args, kwargs, 4, &["a"], "int")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    state.option_colors.insert(name, (r, g, b, a));
    Ok(())
}

#[pyfunction]
#[pyo3(name = "option_get_color", signature = (*args, **kwargs))]
fn option_get_color_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(i32, i32, i32, i32)> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    state.option_colors.get(&name).copied().ok_or_else(|| {
        pyo3::exceptions::PyKeyError::new_err(format!("color option not set: {}", name))
    })
}

#[pyfunction]
#[pyo3(name = "option_restore_defaults", signature = (*_args, **_kwargs))]
fn option_restore_defaults_impl(
    _args: &Bound<'_, PyTuple>,
    _kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    state.option_numbers.clear();
    state.option_strings.clear();
    state.option_colors.clear();
    Ok(())
}

stub_pyfunction!(logger_start_impl, "logger_start", "rmshLoggerStart(int *ierr)");
stub_pyfunction!(logger_stop_impl, "logger_stop", "rmshLoggerStop(int *ierr)");
stub_pyfunction!(logger_get_impl, "logger_get", "rmshLoggerGet(char ***log, size_t *log_n, int *ierr)");

stub_pyfunction!(model_add_impl, "model_add", "rmshModelAdd(const char *name, int *ierr)");
stub_pyfunction!(model_remove_impl, "model_remove", "rmshModelRemove(int *ierr)");
stub_pyfunction!(model_get_current_impl, "model_get_current", "rmshModelGetCurrent(char **name, int *ierr)");
stub_pyfunction!(model_set_current_impl, "model_set_current", "rmshModelSetCurrent(const char *name, int *ierr)");

#[pyfunction]
#[pyo3(name = "model_get_dimension", signature = (*args, **kwargs))]
fn model_get_dimension_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let _ = (args, kwargs);
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    if !state.cad_shapes.is_empty() {
        return Ok(3);
    }
    if let Some(mesh) = state.current_mesh.as_ref() {
        return Ok(mesh_max_dimension(mesh));
    }
    Ok(0)
}

#[pyfunction]
#[pyo3(name = "model_get_entities", signature = (*args, **kwargs))]
fn model_get_entities_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<(i32, i32)>> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int").unwrap_or(-1);

    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut out = Vec::new();
    for tag in state.cad_shapes.keys().copied() {
        if dim < 0 || dim == 3 {
            out.push((3, tag));
        }
    }

    if out.is_empty() {
        if let Some(mesh) = state.current_mesh.as_ref() {
            let mut dims: HashSet<i32> = HashSet::new();
            for e in &mesh.elements {
                dims.insert(i32::from(e.dimension()));
            }
            let mut dims_sorted: Vec<i32> = dims.into_iter().collect();
            dims_sorted.sort_unstable();
            for d in dims_sorted {
                if dim < 0 || dim == d {
                    out.push((d, 1));
                }
            }
        }
    }

    Ok(out)
}

stub_pyfunction!(model_get_entity_name_impl, "model_get_entity_name", "rmshModelGetEntityName(int dim, int tag, char **name, int *ierr)");
stub_pyfunction!(model_set_entity_name_impl, "model_set_entity_name", "rmshModelSetEntityName(int dim, int tag, const char *name, int *ierr)");

#[pyfunction]
#[pyo3(name = "model_get_bounding_box", signature = (*args, **kwargs))]
fn model_get_bounding_box_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(f64, f64, f64, f64, f64, f64)> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int").unwrap_or(-1);
    let tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int").unwrap_or(-1);

    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    if (dim < 0 || dim == 3) && tag > 0 {
        if let Some(shape) = state.cad_shapes.get(&tag) {
            if shape.vertices.is_empty() {
                return Err(pyo3::exceptions::PyRuntimeError::new_err("shape has no vertices"));
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
            return Ok((min[0], min[1], min[2], max[0], max[1], max[2]));
        }
    }

    let mesh = state.current_mesh.as_ref().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded")
    })?;
    let (min, max) = mesh.bounding_box().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("mesh is empty")
    })?;
    Ok((min.x, min.y, min.z, max.x, max.y, max.z))
}

stub_pyfunction!(model_add_physical_group_impl, "model_add_physical_group", "rmshModelAddPhysicalGroup(int dim, const int *tags, size_t tags_n, int tag, const char *name, int *ierr)");
stub_pyfunction!(model_get_physical_groups_impl, "model_get_physical_groups", "rmshModelGetPhysicalGroups(int **dimTags, size_t *dimTags_n, int dim, int *ierr)");
stub_pyfunction!(model_set_physical_name_impl, "model_set_physical_name", "rmshModelSetPhysicalName(int dim, int tag, const char *name, int *ierr)");
stub_pyfunction!(model_get_physical_name_impl, "model_get_physical_name", "rmshModelGetPhysicalName(int dim, int tag, char **name, int *ierr)");

#[pyfunction]
#[pyo3(name = "model_occ_add_box", signature = (*args, **kwargs))]
fn model_occ_add_box_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let x: f64 = extract_required(args, kwargs, 0, &["x"], "float")?;
    let y: f64 = extract_required(args, kwargs, 1, &["y"], "float")?;
    let z: f64 = extract_required(args, kwargs, 2, &["z"], "float")?;
    let dx: f64 = extract_required(args, kwargs, 3, &["dx"], "float")?;
    let dy: f64 = extract_required(args, kwargs, 4, &["dy"], "float")?;
    let dz: f64 = extract_required(args, kwargs, 5, &["dz"], "float")?;
    let tag: i32 = extract_required(args, kwargs, 6, &["tag"], "int").unwrap_or(-1);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut shape = rcad_modeling::box_brep(
        DVec3::new(x, y, z),
        DVec3::X,
        DVec3::Y,
        dx,
        dy,
        dz,
    )
    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    geom_populate::populate_box_geom(&mut shape);

    let assigned_tag = if tag > 0 { tag } else { state.next_cad_tag + 1 };
    state.next_cad_tag = assigned_tag.max(state.next_cad_tag);
    state.cad_shapes.insert(assigned_tag, shape);
    Ok(assigned_tag)
}

#[pyfunction]
#[pyo3(name = "model_occ_add_sphere", signature = (*args, **kwargs))]
fn model_occ_add_sphere_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let x: f64 = extract_required(args, kwargs, 0, &["xc", "x"], "float")?;
    let y: f64 = extract_required(args, kwargs, 1, &["yc", "y"], "float")?;
    let z: f64 = extract_required(args, kwargs, 2, &["zc", "z"], "float")?;
    let r: f64 = extract_required(args, kwargs, 3, &["radius", "r"], "float")?;
    let tag: i32 = extract_required(args, kwargs, 4, &["tag"], "int").unwrap_or(-1);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let shape = rcad_modeling::sphere_brep(DVec3::new(x, y, z), r)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let assigned_tag = if tag > 0 { tag } else { state.next_cad_tag + 1 };
    state.next_cad_tag = assigned_tag.max(state.next_cad_tag);
    state.cad_shapes.insert(assigned_tag, shape);
    Ok(assigned_tag)
}

#[pyfunction]
#[pyo3(name = "model_occ_add_cylinder", signature = (*args, **kwargs))]
fn model_occ_add_cylinder_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let x: f64 = extract_required(args, kwargs, 0, &["x"], "float")?;
    let y: f64 = extract_required(args, kwargs, 1, &["y"], "float")?;
    let z: f64 = extract_required(args, kwargs, 2, &["z"], "float")?;
    let dx: f64 = extract_required(args, kwargs, 3, &["dx"], "float")?;
    let dy: f64 = extract_required(args, kwargs, 4, &["dy"], "float")?;
    let dz: f64 = extract_required(args, kwargs, 5, &["dz"], "float")?;
    let r: f64 = extract_required(args, kwargs, 6, &["r"], "float")?;
    let tag: i32 = extract_required(args, kwargs, 7, &["tag"], "int").unwrap_or(-1);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let axis_vec = DVec3::new(dx, dy, dz);
    let height = axis_vec.length();
    if height < 1e-15 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "cylinder axis direction (dx, dy, dz) must be non-zero",
        ));
    }
    let axis_norm = axis_vec.normalize();
    // Pick a reference direction perpendicular to the axis
    let ref_dir = if axis_norm.x.abs() < 0.9 {
        DVec3::X
    } else {
        DVec3::Y
    };
    let shape = rcad_modeling::cylinder_brep(
        DVec3::new(x, y, z),
        axis_norm,
        ref_dir,
        r,
        height,
    )
    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let assigned_tag = if tag > 0 { tag } else { state.next_cad_tag + 1 };
    state.next_cad_tag = assigned_tag.max(state.next_cad_tag);
    state.cad_shapes.insert(assigned_tag, shape);
    Ok(assigned_tag)
}

#[pyfunction]
#[pyo3(name = "model_occ_cut", signature = (*args, **kwargs))]
fn model_occ_cut_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<(i32, i32)>> {
    let obj_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["objectDimTags"], "list of (dim, tag)")?;
    let tool_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 1, &["toolDimTags"], "list of (dim, tag)")?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut result_shape: Option<BRep> = None;
    for &(_, tag) in &obj_dim_tags {
        if let Some(s) = state.cad_shapes.get(&tag) {
            result_shape = Some(s.clone());
            break;
        }
    }
    let mut base = result_shape.ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err("no valid object shape found for cut")
    })?;

    for &(_, tag) in &tool_dim_tags {
        if let Some(tool) = state.cad_shapes.get(&tag) {
            base = boolean_op(BooleanOpType::Difference, &base, tool)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("boolean cut failed: {e}")))?;
        }
    }

    // Store the result back and tessellate into current mesh
    let first_obj_tag = obj_dim_tags.first().map(|t| t.1).unwrap_or(1);
    let mesh = tessellate_brep(&base);
    state.current_mesh = Some(mesh);
    state.cad_shapes.insert(first_obj_tag, base);

    // Remove consumed tools
    for &(_, tag) in &tool_dim_tags {
        state.cad_shapes.remove(&tag);
    }

    let result_tags: Vec<(i32, i32)> = obj_dim_tags.clone();
    Ok(result_tags)
}

#[pyfunction]
#[pyo3(name = "model_occ_fuse", signature = (*args, **kwargs))]
fn model_occ_fuse_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<(i32, i32)>> {
    let obj_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["objectDimTags"], "list of (dim, tag)")?;
    let tool_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 1, &["toolDimTags"], "list of (dim, tag)")?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut result_shape: Option<BRep> = None;
    for &(_, tag) in &obj_dim_tags {
        if let Some(s) = state.cad_shapes.get(&tag) {
            result_shape = Some(s.clone());
            break;
        }
    }
    let mut base = result_shape.ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err("no valid object shape found for fuse")
    })?;

    for &(_, tag) in &tool_dim_tags {
        if let Some(tool) = state.cad_shapes.get(&tag) {
            base = boolean_op(BooleanOpType::Union, &base, tool)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("boolean fuse failed: {e}")))?;
        }
    }

    let first_obj_tag = obj_dim_tags.first().map(|t| t.1).unwrap_or(1);
    let mesh = tessellate_brep(&base);
    state.current_mesh = Some(mesh);
    state.cad_shapes.insert(first_obj_tag, base);

    for &(_, tag) in &tool_dim_tags {
        state.cad_shapes.remove(&tag);
    }

    let result_tags: Vec<(i32, i32)> = obj_dim_tags.clone();
    Ok(result_tags)
}

#[pyfunction]
#[pyo3(name = "model_occ_fragment", signature = (*args, **kwargs))]
fn model_occ_fragment_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<(i32, i32)>> {
    let obj_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["objectDimTags"], "list of (dim, tag)")?;
    let tool_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 1, &["toolDimTags"], "list of (dim, tag)")?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut result_shape: Option<BRep> = None;
    for &(_, tag) in &obj_dim_tags {
        if let Some(s) = state.cad_shapes.get(&tag) {
            result_shape = Some(s.clone());
            break;
        }
    }
    let mut base = result_shape.ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err("no valid object shape found for fragment")
    })?;

    for &(_, tag) in &tool_dim_tags {
        if let Some(tool) = state.cad_shapes.get(&tag) {
            base = boolean_op(BooleanOpType::Intersection, &base, tool)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("boolean fragment failed: {e}")))?;
        }
    }

    let first_obj_tag = obj_dim_tags.first().map(|t| t.1).unwrap_or(1);
    let mesh = tessellate_brep(&base);
    state.current_mesh = Some(mesh);
    state.cad_shapes.insert(first_obj_tag, base);

    let mut result_tags: Vec<(i32, i32)> = obj_dim_tags.clone();
    result_tags.extend(tool_dim_tags.iter().copied());
    Ok(result_tags)
}

#[pyfunction]
#[pyo3(name = "model_occ_synchronize", signature = (*args, **kwargs))]
fn model_occ_synchronize_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let _ = (args, kwargs);
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    // Tessellate all CAD shapes and merge into current_mesh
    let tags: Vec<i32> = state.cad_shapes.keys().copied().collect();
    for tag in tags {
        if let Some(shape) = state.cad_shapes.get(&tag) {
            let mesh = tessellate_brep(shape);
            match state.current_mesh.as_mut() {
                Some(current) => merge_meshes(current, &mesh),
                None => state.current_mesh = Some(mesh),
            }
        }
    }
    state.cad_shapes.clear();
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_mesh_set_size", signature = (*args, **kwargs))]
fn model_mesh_set_size_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let size: f64 = extract_required(args, kwargs, 1, &["size"], "float")?;
    if size <= 0.0 {
        return Err(pyo3::exceptions::PyValueError::new_err("size must be > 0"));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.option_numbers.insert("Mesh.MeshSizeMin".to_string(), size);
    state.option_numbers.insert("Mesh.MeshSizeMax".to_string(), size);
    state.option_numbers.insert("Mesh.CharacteristicLengthMin".to_string(), size);
    state.option_numbers.insert("Mesh.CharacteristicLengthMax".to_string(), size);
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_mesh_generate", signature = (*args, **kwargs))]
fn model_mesh_generate_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int")?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let surface = state.current_mesh.clone().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded; call open() first")
    })?;

    // Resolve element size from options (new names take priority over deprecated ones)
    let size_max = option_number(&state, "Mesh.MeshSizeMax")
        .or_else(|| option_number(&state, "Mesh.CharacteristicLengthMax"))
        .filter(|v| *v > 0.0 && *v < 1e20);
    let size_min = option_number(&state, "Mesh.MeshSizeMin")
        .or_else(|| option_number(&state, "Mesh.CharacteristicLengthMin"))
        .filter(|v| *v > 0.0);
    let size_factor = option_number(&state, "Mesh.MeshSizeFactor")
        .or_else(|| option_number(&state, "Mesh.CharacteristicLengthFactor"))
        .filter(|v| *v > 0.0)
        .unwrap_or(1.0);

    let base_size = size_max.unwrap_or(1.0) * size_factor;
    let mut params = MeshParams::with_size(base_size);
    if let Some(v) = size_min { params.min_size = v; }
    if let Some(v) = size_max { params.max_size = v * size_factor; }
    if params.max_size < params.min_size {
        std::mem::swap(&mut params.max_size, &mut params.min_size);
    }

    // Read algorithm selectors
    let algo_2d = option_number(&state, "Mesh.Algorithm").map(|v| v as i32).unwrap_or(6);
    let algo_3d = option_number(&state, "Mesh.Algorithm3D").map(|v| v as i32).unwrap_or(1);

    let convert_err = |e: MeshAlgoError| pyo3::exceptions::PyRuntimeError::new_err(e.to_string());

    let generated = if dim == 3 {
        // 3D algorithms: 1=Delaunay, 4=Frontal, 10=HXT, else=CentroidStar(default)
        match algo_3d {
            1 => Delaunay3D::default().mesh_3d(&surface, &params).map_err(convert_err)?,
            4 => Frontal3D::default().mesh_3d(&surface, &params).map_err(convert_err)?,
            10 => Hxt3D::default().mesh_3d(&surface, &params).map_err(convert_err)?,
            _ => CentroidStarMesher3D.mesh_3d(&surface, &params).map_err(convert_err)?,
        }
    } else if dim == 2 {
        let polygon = boundary_loop_from_surface_mesh(&surface)?;
        let domain = Domain2D::from_outer(polygon);
        // 2D algorithms: 5/6=Frontal-Delaunay, 7=BAMG, 8/9=Quad, else=basic triangulate
        match algo_2d {
            5 | 6 => FrontalDelaunay2D::default().mesh_2d(&domain, &params).map_err(convert_err)?,
            7 => Bamg2D::default().mesh_2d(&domain, &params).map_err(convert_err)?,
            8 | 9 => QuadPaving2D::default().mesh_2d(&domain, &params).map_err(convert_err)?,
            _ => mesh_polygon(&Polygon2D::new(domain.outer().to_vec()), params.element_size)
                    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?,
        }
    } else {
        return Err(PyNotImplementedError::new_err(
            "only dim=2 and dim=3 are currently implemented",
        ));
    };
    state.current_mesh = Some(generated);
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_mesh_set_order", signature = (*args, **kwargs))]
fn model_mesh_set_order_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let order: i32 = extract_required(args, kwargs, 0, &["order"], "int")?;
    if order != 1 {
        return Err(PyNotImplementedError::new_err(
            "only first-order elements are currently supported",
        ));
    }
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_mesh_get_nodes", signature = (*args, **kwargs))]
fn model_mesh_get_nodes_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(Vec<u64>, Vec<f64>, Vec<f64>)> {
    let dim: i32 = extract_required(args, kwargs, 3, &["dim"], "int").unwrap_or(-1);
    let return_parametric: i32 =
        extract_required(args, kwargs, 6, &["returnParametricCoord"], "int").unwrap_or(0);

    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mesh = state.current_mesh.as_ref().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded")
    })?;

    if dim >= 0
        && !mesh
            .elements
            .iter()
            .any(|e| i32::from(e.dimension()) == dim)
    {
        return Ok((Vec::new(), Vec::new(), Vec::new()));
    }

    let mut node_ids: Vec<u64> = mesh.nodes.keys().copied().collect();
    node_ids.sort_unstable();

    let mut coords = Vec::with_capacity(node_ids.len() * 3);
    for id in &node_ids {
        if let Some(n) = mesh.nodes.get(id) {
            coords.push(n.position.x);
            coords.push(n.position.y);
            coords.push(n.position.z);
        }
    }

    let parametric = if return_parametric != 0 {
        vec![0.0; node_ids.len()]
    } else {
        Vec::new()
    };

    Ok((node_ids, coords, parametric))
}

#[pyfunction]
#[pyo3(name = "model_mesh_get_elements", signature = (*args, **kwargs))]
fn model_mesh_get_elements_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(Vec<i32>, Vec<Vec<u64>>, Vec<Vec<u64>>)> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int").unwrap_or(-1);

    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mesh = state.current_mesh.as_ref().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded")
    })?;

    let mut grouped: BTreeMap<i32, (Vec<u64>, Vec<u64>)> = BTreeMap::new();
    for elem in &mesh.elements {
        if dim >= 0 && i32::from(elem.dimension()) != dim {
            continue;
        }
        let type_id = gmsh_type_id(elem.etype);
        let entry = grouped.entry(type_id).or_insert_with(|| (Vec::new(), Vec::new()));
        entry.0.push(elem.id);
        entry.1.extend(elem.node_ids.iter().copied());
    }

    let mut element_types = Vec::with_capacity(grouped.len());
    let mut element_tags = Vec::with_capacity(grouped.len());
    let mut node_tags = Vec::with_capacity(grouped.len());
    for (etype, (tags, nodes)) in grouped {
        element_types.push(etype);
        element_tags.push(tags);
        node_tags.push(nodes);
    }

    Ok((element_types, element_tags, node_tags))
}

#[pyfunction]
#[pyo3(name = "model_mesh_clear", signature = (*args, **kwargs))]
fn model_mesh_clear_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let _ = (args, kwargs);
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.current_mesh = None;
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_mesh_optimize", signature = (*args, **kwargs))]
fn model_mesh_optimize_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let method: String = if args.len() > 0 || kwargs.map(|k| k.contains("method").unwrap_or(false)).unwrap_or(false) {
        extract_required(args, kwargs, 0, &["method"], "str").unwrap_or_else(|_| "Laplace".to_string())
    } else {
        "Laplace".to_string()
    };
    let niter: u32 = extract_required(args, kwargs, 2, &["niter"], "int").unwrap_or(10);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let mesh = state.current_mesh.as_mut().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded")
    })?;

    let params = OptimizeParams {
        iterations: niter,
        ..Default::default()
    };

    match method.as_str() {
        "Laplace" | "Laplacian" | "" => {
            LaplacianSmooth::default()
                .optimize(mesh, &params)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        }
        other => {
            return Err(PyNotImplementedError::new_err(format!(
                "optimizer '{}' not yet implemented; available: Laplace",
                other
            )));
        }
    }
    Ok(())
}

stub_pyfunction!(model_mesh_refine_impl, "model_mesh_refine", "rmshModelMeshRefine(int *ierr)");
stub_pyfunction!(model_mesh_recombine_impl, "model_mesh_recombine", "rmshModelMeshRecombine(int dim, int tag, double angle, int *ierr)");

stub_pyfunction!(plugin_set_number_impl, "plugin_set_number", "rmshPluginSetNumber(const char *name, const char *option, double value, int *ierr)");
stub_pyfunction!(plugin_set_string_impl, "plugin_set_string", "rmshPluginSetString(const char *name, const char *option, const char *value, int *ierr)");
stub_pyfunction!(plugin_run_impl, "plugin_run", "rmshPluginRun(const char *name, int *ierr)");

#[pyfunction]
#[pyo3(name = "gui_initialize", signature = (*args, **kwargs))]
fn gui_initialize_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let _ = (args, kwargs);
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)
}

#[pyfunction]
#[pyo3(name = "gui_run", signature = (*args, **kwargs))]
fn gui_run_impl(
    py: Python<'_>,
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let _ = (args, kwargs);

    let (mesh, mesh_name) = {
        let state = STATE
            .lock()
            .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
        ensure_initialized(&state)?;
        let mesh = state.current_mesh.clone().ok_or_else(|| {
            pyo3::exceptions::PyRuntimeError::new_err("no current model/mesh; call open() or generate() first")
        })?;
        let mesh_name = state
            .current_path
            .as_ref()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .map(ToOwned::to_owned)
            .unwrap_or_else(|| "rmsh-python-session.msh".to_string());
        (mesh, mesh_name)
    };

    py.allow_threads(move || rmsh_viewer::run_native_viewer(None, Some((mesh, mesh_name))))
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

#[pyfunction]
#[pyo3(name = "gui_wait", signature = (*args, **kwargs))]
fn gui_wait_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let timeout: f64 = if args.len() > 0 || kwargs.is_some() {
        extract_required(args, kwargs, 0, &["time", "timeout"], "float")?
    } else {
        0.0
    };

    if timeout > 0.0 {
        std::thread::sleep(Duration::from_secs_f64(timeout));
    }
    Ok(())
}

// ── Shape properties ──────────────────────────────────────────────────────────

/// Return the surface area of the current CAD shape or mesh.
#[pyfunction]
#[pyo3(name = "model_occ_get_mass", signature = (*args, **kwargs))]
fn model_occ_get_mass_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<f64> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let brep = state.cad_shapes.get(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?;
    Ok(rcad_kernel::properties::volume(brep).abs())
}

/// Return (volume, surface_area, cx, cy, cz) for a CAD shape.
#[pyfunction]
#[pyo3(name = "model_occ_get_properties", signature = (*args, **kwargs))]
fn model_occ_get_properties_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(f64, f64, f64, f64, f64)> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let brep = state.cad_shapes.get(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?;
    let vol = rcad_kernel::properties::volume(brep).abs();
    let area = rcad_kernel::properties::surface_area(brep);
    let c = rcad_kernel::properties::centroid(brep);
    Ok((vol, area, c.x, c.y, c.z))
}

// ── Extrude / Revolve ─────────────────────────────────────────────────────────

/// Extrude face `face_idx` of shape `tag` along `(dx,dy,dz)` by `distance`.
/// Returns a new tag for the resulting solid.
#[pyfunction]
#[pyo3(name = "model_occ_extrude", signature = (*args, **kwargs))]
fn model_occ_extrude_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;
    let face_idx: usize = extract_required(args, kwargs, 1, &["face_idx"], "int")?;
    let dx: f64 = extract_required(args, kwargs, 2, &["dx"], "float")?;
    let dy: f64 = extract_required(args, kwargs, 3, &["dy"], "float")?;
    let dz: f64 = extract_required(args, kwargs, 4, &["dz"], "float")?;
    let distance: f64 = extract_required(args, kwargs, 5, &["distance"], "float")?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let base = state.cad_shapes.get(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?.clone();

    let result = extrude(&base, face_idx, DVec3::new(dx, dy, dz), distance)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let new_tag = state.next_cad_tag + 1;
    state.next_cad_tag = new_tag;
    state.cad_shapes.insert(new_tag, result);
    Ok(new_tag)
}

/// Revolve face `face_idx` of shape `tag` around axis through `(ax,ay,az)`
/// in direction `(dx,dy,dz)` by `angle` radians. Returns a new tag.
#[pyfunction]
#[pyo3(name = "model_occ_revolve", signature = (*args, **kwargs))]
fn model_occ_revolve_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;
    let face_idx: usize = extract_required(args, kwargs, 1, &["face_idx"], "int")?;
    let ax: f64 = extract_required(args, kwargs, 2, &["ax"], "float")?;
    let ay: f64 = extract_required(args, kwargs, 3, &["ay"], "float")?;
    let az: f64 = extract_required(args, kwargs, 4, &["az"], "float")?;
    let dx: f64 = extract_required(args, kwargs, 5, &["dx"], "float")?;
    let dy: f64 = extract_required(args, kwargs, 6, &["dy"], "float")?;
    let dz: f64 = extract_required(args, kwargs, 7, &["dz"], "float")?;
    let angle: f64 = extract_required(args, kwargs, 8, &["angle"], "float")?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let base = state.cad_shapes.get(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?.clone();

    let result = revolve(&base, face_idx, DVec3::new(ax, ay, az), DVec3::new(dx, dy, dz), angle)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let new_tag = state.next_cad_tag + 1;
    state.next_cad_tag = new_tag;
    state.cad_shapes.insert(new_tag, result);
    Ok(new_tag)
}

// ── Cone / Torus ──────────────────────────────────────────────────────────────

/// Add a cone with base center at (x,y,z), axis direction (dx,dy,dz), base radius r, height h.
/// Returns an integer tag for the new shape.
#[pyfunction]
#[pyo3(name = "model_occ_add_cone", signature = (*args, **kwargs))]
fn model_occ_add_cone_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let x: f64 = extract_required(args, kwargs, 0, &["x"], "float")?;
    let y: f64 = extract_required(args, kwargs, 1, &["y"], "float")?;
    let z: f64 = extract_required(args, kwargs, 2, &["z"], "float")?;
    let dx: f64 = extract_required(args, kwargs, 3, &["dx"], "float")?;
    let dy: f64 = extract_required(args, kwargs, 4, &["dy"], "float")?;
    let dz: f64 = extract_required(args, kwargs, 5, &["dz"], "float")?;
    let r: f64 = extract_required(args, kwargs, 6, &["r"], "float")?;
    let tag: i32 = extract_required(args, kwargs, 7, &["tag"], "int").unwrap_or(-1);

    let axis_vec = DVec3::new(dx, dy, dz);
    let height = axis_vec.length();
    if height < 1e-15 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "cone axis direction (dx, dy, dz) must be non-zero",
        ));
    }
    let axis_norm = axis_vec.normalize();
    let ref_dir = if axis_norm.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut shape = cone_brep(DVec3::new(x, y, z), axis_norm, ref_dir, r, height)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    geom_populate::populate_box_geom(&mut shape);

    let assigned_tag = if tag > 0 { tag } else { state.next_cad_tag + 1 };
    state.next_cad_tag = assigned_tag.max(state.next_cad_tag);
    state.cad_shapes.insert(assigned_tag, shape);
    Ok(assigned_tag)
}

/// Add a torus with center at (x,y,z), axis direction (dx,dy,dz),
/// major radius R and tube radius r.
/// Returns an integer tag for the new shape.
#[pyfunction]
#[pyo3(name = "model_occ_add_torus", signature = (*args, **kwargs))]
fn model_occ_add_torus_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let x: f64 = extract_required(args, kwargs, 0, &["x"], "float")?;
    let y: f64 = extract_required(args, kwargs, 1, &["y"], "float")?;
    let z: f64 = extract_required(args, kwargs, 2, &["z"], "float")?;
    let dx: f64 = extract_required(args, kwargs, 3, &["dx"], "float").unwrap_or(0.0);
    let dy: f64 = extract_required(args, kwargs, 4, &["dy"], "float").unwrap_or(0.0);
    let dz: f64 = extract_required(args, kwargs, 5, &["dz"], "float").unwrap_or(1.0);
    let r1: f64 = extract_required(args, kwargs, 6, &["r1"], "float")?;
    let r2: f64 = extract_required(args, kwargs, 7, &["r2"], "float")?;
    let tag: i32 = extract_required(args, kwargs, 8, &["tag"], "int").unwrap_or(-1);

    let axis_vec = DVec3::new(dx, dy, dz);
    let axis_norm = if axis_vec.length_squared() < 1e-20 {
        DVec3::Z
    } else {
        axis_vec.normalize()
    };
    let ref_dir = if axis_norm.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let shape = torus_brep(DVec3::new(x, y, z), axis_norm, ref_dir, r1, r2)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let assigned_tag = if tag > 0 { tag } else { state.next_cad_tag + 1 };
    state.next_cad_tag = assigned_tag.max(state.next_cad_tag);
    state.cad_shapes.insert(assigned_tag, shape);
    Ok(assigned_tag)
}

// ── Rectangle (2D domain for surface meshing) ─────────────────────────────────

/// Add a planar rectangle in the XY plane.
/// Signature matches gmsh.model.occ.addRectangle:
///   addRectangle(x, y, z, dx, dy, tag=-1) -> tag
///
/// Unlike addBox this creates a *surface* entity (2D boundary mesh in current_mesh)
/// suitable for 2D meshing with model.mesh.generate(2).
#[pyfunction]
#[pyo3(name = "model_occ_add_rectangle", signature = (*args, **kwargs))]
fn model_occ_add_rectangle_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let x: f64  = extract_required(args, kwargs, 0, &["x"],  "float")?;
    let y: f64  = extract_required(args, kwargs, 1, &["y"],  "float")?;
    let _z: f64 = extract_required(args, kwargs, 2, &["z"],  "float").unwrap_or(0.0);
    let dx: f64 = extract_required(args, kwargs, 3, &["dx"], "float")?;
    let dy: f64 = extract_required(args, kwargs, 4, &["dy"], "float")?;
    let tag: i32 = extract_required(args, kwargs, 5, &["tag"], "int").unwrap_or(-1);

    // Build a single quad boundary mesh of 2 triangles representing the rectangle
    let mut mesh = Mesh::new();
    // Nodes: four corners
    mesh.add_node(Node::new(1, x,      y,      0.0));
    mesh.add_node(Node::new(2, x + dx, y,      0.0));
    mesh.add_node(Node::new(3, x + dx, y + dy, 0.0));
    mesh.add_node(Node::new(4, x,      y + dy, 0.0));
    // Two triangles to fill the rectangle
    mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));
    mesh.add_element(Element::new(2, ElementType::Triangle3, vec![1, 3, 4]));

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let assigned_tag = if tag > 0 { tag } else { state.next_cad_tag + 1 };
    state.next_cad_tag = assigned_tag.max(state.next_cad_tag);
    // Store surface mesh as current_mesh so generate(2) can use it
    state.current_mesh = Some(mesh);
    Ok(assigned_tag)
}

/// Round the edges of a volume.
/// Signature matches gmsh.model.occ.fillet:
///   fillet(tag, curveTags, radii) -> new_tag
/// `curveTags`: list of edge indices (0-based) to fillet.
/// `radii`: list of radii, one per edge, or a single value applied to all.
/// Returns a new tag for the modified shape.
#[pyfunction]
#[pyo3(name = "model_occ_fillet", signature = (*args, **kwargs))]
fn model_occ_fillet_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;

    let curve_tags: Vec<usize> = args.get_item(1).ok()
        .or_else(|| kwargs.and_then(|kw| kw.get_item("curveTags").ok().flatten()))
        .ok_or_else(|| pyo3::exceptions::PyTypeError::new_err("curveTags required"))?
        .extract::<Vec<usize>>()
        .map_err(|_| pyo3::exceptions::PyTypeError::new_err("curveTags must be a list of ints"))?;

    let radii_raw: Vec<f64> = args.get_item(2).ok()
        .or_else(|| kwargs.and_then(|kw| kw.get_item("radii").ok().flatten()))
        .ok_or_else(|| pyo3::exceptions::PyTypeError::new_err("radii required"))?
        .extract::<Vec<f64>>()
        .map_err(|_| pyo3::exceptions::PyTypeError::new_err("radii must be a list of floats"))?;

    // Expand scalar radius to per-edge list
    let radii: Vec<f64> = if radii_raw.len() == 1 {
        vec![radii_raw[0]; curve_tags.len()]
    } else {
        radii_raw
    };
    if radii.len() != curve_tags.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "radii length must match curveTags length (or be a single value)"
        ));
    }

    let edges: Vec<(usize, f64)> = curve_tags.into_iter().zip(radii).collect();

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let base = state.cad_shapes.get(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?.clone();

    let result = fillet_edges(&base, &edges)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let new_tag = state.next_cad_tag + 1;
    state.next_cad_tag = new_tag;
    state.cad_shapes.insert(new_tag, result);
    Ok(new_tag)
}

/// Chamfer the edges of a volume.
/// Signature matches gmsh.model.occ.chamfer:
///   chamfer(tag, curveTags, distances) -> new_tag
/// `curveTags`: list of edge indices (0-based) to chamfer.
/// `distances`: list of chamfer distances, one per edge, or a single value.
/// Returns a new tag for the modified shape.
#[pyfunction]
#[pyo3(name = "model_occ_chamfer", signature = (*args, **kwargs))]
fn model_occ_chamfer_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;

    let curve_tags: Vec<usize> = args.get_item(1).ok()
        .or_else(|| kwargs.and_then(|kw| kw.get_item("curveTags").ok().flatten()))
        .ok_or_else(|| pyo3::exceptions::PyTypeError::new_err("curveTags required"))?
        .extract::<Vec<usize>>()
        .map_err(|_| pyo3::exceptions::PyTypeError::new_err("curveTags must be a list of ints"))?;

    let distances_raw: Vec<f64> = args.get_item(2).ok()
        .or_else(|| kwargs.and_then(|kw| kw.get_item("distances").ok().flatten()))
        .ok_or_else(|| pyo3::exceptions::PyTypeError::new_err("distances required"))?
        .extract::<Vec<f64>>()
        .map_err(|_| pyo3::exceptions::PyTypeError::new_err("distances must be a list of floats"))?;

    let distances: Vec<f64> = if distances_raw.len() == 1 {
        vec![distances_raw[0]; curve_tags.len()]
    } else {
        distances_raw
    };
    if distances.len() != curve_tags.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "distances length must match curveTags length (or be a single value)"
        ));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let base = state.cad_shapes.get(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?.clone();

    // Apply chamfers sequentially (descending edge index to preserve indices)
    let mut edges: Vec<(usize, f64)> = curve_tags.into_iter().zip(distances).collect();
    edges.sort_by(|a, b| b.0.cmp(&a.0));
    let mut current = base;
    for (edge_idx, dist) in edges {
        current = chamfer_edge(&current, edge_idx, dist)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
    }

    let new_tag = state.next_cad_tag + 1;
    state.next_cad_tag = new_tag;
    state.cad_shapes.insert(new_tag, current);
    Ok(new_tag)
}

/// Heal (repair) a shape: merge close vertices, recompute normals, fix wire orientations.
/// Signature matches gmsh.model.occ.healShapes:
///   heal_shapes(tag, tolerance=1e-8) -> report_dict
/// Updates the shape in-place. Returns a dict with repair counts.
#[pyfunction]
#[pyo3(name = "model_occ_heal_shapes", signature = (*args, **kwargs))]
fn model_occ_heal_shapes_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<pyo3::PyObject> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;
    let tolerance: f64 = extract_required(args, kwargs, 1, &["tolerance"], "float")
        .unwrap_or(1e-8);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let brep = state.cad_shapes.get_mut(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?;

    let (repaired, report) = repair(brep, tolerance);
    *brep = repaired;

    Python::with_gil(|py| {
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("vertices_merged", report.vertices_merged)?;
        dict.set_item("degenerate_faces_removed", report.degenerate_faces_removed)?;
        dict.set_item("normals_recomputed", report.normals_recomputed)?;
        dict.set_item("wires_fixed", report.wires_fixed)?;
        Ok(dict.into())
    })
}

#[pyo3::prelude::pymodule]
fn _rmsh(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(pyo3::wrap_pyfunction!(initialize_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(finalize_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(clear_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(open_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(merge_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(write_impl, m)?)?;

    m.add_function(pyo3::wrap_pyfunction!(option_set_number_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(option_get_number_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(option_set_string_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(option_get_string_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(option_set_color_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(option_get_color_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(option_restore_defaults_impl, m)?)?;

    m.add_function(pyo3::wrap_pyfunction!(logger_start_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(logger_stop_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(logger_get_impl, m)?)?;

    m.add_function(pyo3::wrap_pyfunction!(model_add_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_remove_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_get_current_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_set_current_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_get_dimension_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_get_entities_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_get_entity_name_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_set_entity_name_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_get_bounding_box_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_add_physical_group_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_get_physical_groups_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_set_physical_name_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_get_physical_name_impl, m)?)?;

    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_box_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_sphere_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_cylinder_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_cone_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_torus_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_rectangle_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_cut_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_fuse_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_fragment_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_synchronize_impl, m)?)?;

    m.add_function(pyo3::wrap_pyfunction!(model_mesh_set_size_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_generate_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_set_order_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_get_nodes_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_get_elements_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_clear_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_optimize_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_refine_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_recombine_impl, m)?)?;

    m.add_function(pyo3::wrap_pyfunction!(plugin_set_number_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(plugin_set_string_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(plugin_run_impl, m)?)?;

    m.add_function(pyo3::wrap_pyfunction!(gui_initialize_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(gui_run_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(gui_wait_impl, m)?)?;

    m.add_function(pyo3::wrap_pyfunction!(model_occ_get_mass_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_get_properties_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_extrude_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_revolve_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_fillet_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_chamfer_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_heal_shapes_impl, m)?)?;

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
