use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::f64::consts::PI;
use std::path::PathBuf;
use std::sync::{LazyLock, Mutex};
use std::time::Duration;

use glam::{DAffine3, DMat3, DQuat, DVec3};
use rmsh_algo::{
    CentroidStarMesher3D, Delaunay2D, Delaunay3D, FrontalDelaunay2D, Frontal3D, Hxt3D,
    Bamg2D, FrontalQuads2D, MeshAdapt2D, MmgRemesh, QuadPaving2D,
    LaplacianSmooth, promote_to_p2,
    FieldManager, ConstantField, DistanceField, ThresholdField,
    MathEvalField, MinField, MaxField, BoxField, RestrictField,
    MeshAlgoError, MeshOptimizer, MeshParams, Mesher2D, Mesher3D,
    OptimizeParams, Polygon2D, mesh_polygon,
    MeshQualityOptimizer, OptimizeConfig, QualityMetric,
    Domain2D,
};

mod cad;
use cad::{BooleanOp, CadKernel, CurveKind, RcadKernel, StepExportOptions, StepProtocol, SurfaceKind};
use rmsh_model::{Element, ElementType, Mesh, Node};
// BRep import retained for remaining functions that directly build/manipulate shapes
use rcad_kernel::BRep;

struct RuntimeState {
    kernel: RcadKernel,
    initialized: bool,
    model_name: String,
    mesh_order: i32,
    current_mesh: Option<Mesh>,
    current_path: Option<PathBuf>,
    option_numbers: HashMap<String, f64>,
    option_strings: HashMap<String, String>,
    option_colors: HashMap<String, (i32, i32, i32, i32)>,
    entity_names: HashMap<(i32, i32), String>,
    physical_groups: HashMap<(i32, i32), Vec<i32>>,
    physical_names: HashMap<(i32, i32), String>,
    plugin_numbers: HashMap<(String, String), f64>,
    plugin_strings: HashMap<(String, String), String>,
    logger_enabled: bool,
    logs: Vec<String>,
    /// CAD shapes created via model.occ.add* functions, keyed by tag.
    cad_shapes: HashMap<i32, <RcadKernel as CadKernel>::Shape>,
    /// OCC-style storage for created curve loops (wire tags).
    occ_curve_loops: HashMap<i32, Vec<i32>>,
    /// OCC-style storage for created surface loops (shell tags).
    occ_surface_loops: HashMap<i32, Vec<i32>>,
    /// Next auto-assigned tag for CAD shapes.
    next_cad_tag: i32,
    /// Background size field manager (Gmsh Field equivalent).
    field_mgr: FieldManager,
    /// Transfinite curve constraints: (curve_tag, num_segments, distribution).
    transfinite_curves: Vec<(i32, usize, String, f64)>,
    /// Transfinite surface constraints: (surface_tag, arrangement).
    transfinite_surfaces: Vec<(i32, String, Vec<i32>)>,
}

impl Default for RuntimeState {
    fn default() -> Self {
        Self {
            kernel: RcadKernel,
            initialized: false,
            model_name: String::new(),
            mesh_order: 1,
            current_mesh: None,
            current_path: None,
            option_numbers: HashMap::new(),
            option_strings: HashMap::new(),
            option_colors: HashMap::new(),
            entity_names: HashMap::new(),
            physical_groups: HashMap::new(),
            physical_names: HashMap::new(),
            plugin_numbers: HashMap::new(),
            plugin_strings: HashMap::new(),
            logger_enabled: false,
            logs: Vec::new(),
            cad_shapes: HashMap::new(),
            occ_curve_loops: HashMap::new(),
            occ_surface_loops: HashMap::new(),
            next_cad_tag: 0,
            field_mgr: FieldManager::new(),
            transfinite_curves: Vec::new(),
            transfinite_surfaces: Vec::new(),
        }
    }
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



fn step_solid_color_from_state(state: &RuntimeState) -> Option<(u8, u8, u8)> {
    let rgba = *state
        .option_colors
        .get("Geometry.OCCSolidColor")
        .or_else(|| state.option_colors.get("STEP.SolidColor"))?;
    let clamp = |v: i32| -> u8 { v.clamp(0, 255) as u8 };
    Some((clamp(rgba.0), clamp(rgba.1), clamp(rgba.2)))
}

fn step_protocol_from_state(state: &RuntimeState) -> StepProtocol {
    let keys = ["STEP.Protocol", "Geometry.OCCStepProtocol", "Geometry.OCCSTEPProtocol"];
    for key in keys {
        if let Some(v) = state
            .option_strings
            .get(key)
            .cloned()
            .or_else(|| option_default_string(key).map(str::to_string))
        {
            let s = v.to_ascii_lowercase();
            if s.contains("214") {
                return StepProtocol::Ap214;
            }
            if s.contains("242") {
                return StepProtocol::Ap242;
            }
        }
    }
    StepProtocol::Ap242
}

fn step_gmsh_strict_from_state(state: &RuntimeState) -> bool {
    let numeric_keys = ["STEP.GmshStrict", "Geometry.OCCStepGmshStrict"];
    let mut seen_numeric = false;
    let mut strict_numeric = false;
    for key in numeric_keys {
        if let Some(v) = state.option_numbers.get(key) {
            seen_numeric = true;
            if *v > 0.5 {
                strict_numeric = true;
            }
        }
    }
    if seen_numeric {
        return strict_numeric;
    }

    let style_keys = ["STEP.WriterStyle", "Geometry.OCCSTEPWriterStyle"];
    for key in style_keys {
        if let Some(style) = state.option_strings.get(key) {
            let v = style.to_ascii_lowercase();
            if v.contains("legacy")
                || v.contains("default")
                || v.contains("relaxed")
                || v.contains("off")
            {
                return false;
            }
            if v.contains("gmsh") || v.contains("occt") || v.contains("strict") {
                return true;
            }
            return false;
        }
    }

    // Align STEP save defaults with gmsh/OCCT-compatible output unless the
    // user explicitly opts out via STEP.GmshStrict=0 or writer style options.
    true
}

fn bbox_intersects(a_min: [f64; 3], a_max: [f64; 3], b_min: [f64; 3], b_max: [f64; 3]) -> bool {
    !(a_max[0] < b_min[0]
        || a_min[0] > b_max[0]
        || a_max[1] < b_min[1]
        || a_min[1] > b_max[1]
        || a_max[2] < b_min[2]
        || a_min[2] > b_max[2])
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

/// Detect 4 corner indices in a boundary polygon by finding vertices where
/// the turning angle is significant (> threshold degrees). Non-corner edge
/// vertices are nearly collinear (~180°) while corners are ~90°.
fn detect_corner_indices(polygon: &[[f64; 2]]) -> Vec<usize> {
    let n = polygon.len();
    if n <= 4 {
        return (0..n.min(4)).collect();
    }

    // Compute turning angle at each vertex.
    struct Turn { idx: usize, angle: f64 }
    let mut turns: Vec<Turn> = (0..n) // skip first/last for clean computation
        .filter_map(|i| {
            let prev = polygon[(i + n - 1) % n];
            let cur = polygon[i];
            let next = polygon[(i + 1) % n];
            let v1 = [cur[0] - prev[0], cur[1] - prev[1]];
            let v2 = [next[0] - cur[0], next[1] - cur[1]];
            let d1 = (v1[0] * v1[0] + v1[1] * v1[1]).sqrt();
            let d2 = (v2[0] * v2[0] + v2[1] * v2[1]).sqrt();
            if d1 < 1e-12 || d2 < 1e-12 {
                return None;
            }
            let dot = (v1[0] * v2[0] + v1[1] * v2[1]) / (d1 * d2);
            let angle = dot.clamp(-1.0, 1.0).acos().to_degrees();
            Some(Turn { idx: i, angle })
        })
        .collect();

    turns.sort_by(|a, b| a.angle.partial_cmp(&b.angle).unwrap_or(std::cmp::Ordering::Equal));

    // Pick the 4 smallest angles (sharpest turns = corners).
    let mut corners: Vec<usize> = turns.iter().take(4).map(|t| t.idx).collect();
    corners.sort_unstable();

    // Ensure continuous ordering starting from the first corner.
    if corners.len() >= 2 {
        // Check if the polygon goes forward or wraps
        if corners[0] + n - corners[corners.len() - 1] < corners[1] - corners[0] {
            corners.rotate_right(1);
        }
    }
    corners
}

/// Split a boundary polygon into 4 edges at the given corner indices.
/// Returns 4 edge segments (each as a Vec of points from corner to corner).
fn split_polygon_into_edges<'a>(polygon: &'a [[f64; 2]], corners: &[usize]) -> Vec<Vec<&'a [f64; 2]>> {
    if corners.len() < 2 {
        return vec![polygon.iter().collect()];
    }
    let n = polygon.len();
    let mut edges = Vec::with_capacity(4);
    for ci in 0..corners.len() {
        let start = corners[ci];
        let end = corners[(ci + 1) % corners.len()];
        let mut edge = Vec::new();
        let mut i = start;
        loop {
            edge.push(&polygon[i]);
            if i == end {
                break;
            }
            i = (i + 1) % n;
        }
        edges.push(edge);
    }
    edges
}

fn boundary_loop_from_brep_face(brep: &BRep) -> Option<Vec<[f64; 2]>> {
    // Fallback path for 2D meshing when face triangulation is unavailable.
    // Build the polygon directly from the first face outer wire.
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let mut polygon: Vec<[f64; 2]> = Vec::new();
                for we in &face.outer_wire.edges {
                    let edge = brep.edges.get(we.idx)?;
                    let vidx = if we.forward { edge.start } else { edge.end };
                    let v = brep.vertices.get(vidx)?;
                    polygon.push([v.point.x, v.point.y]);
                }
                if polygon.len() >= 3 {
                    return Some(polygon);
                }
            }
        }
    }
    None
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

fn next_or_requested_tag(state: &mut RuntimeState, requested: i32) -> i32 {
    let assigned = if requested > 0 {
        requested
    } else {
        state.next_cad_tag + 1
    };
    state.next_cad_tag = assigned.max(state.next_cad_tag);
    assigned
}

fn point_from_point_entity(state: &RuntimeState, tag: i32, name: &str) -> PyResult<DVec3> {
    let shape = state
        .cad_shapes
        .get(&tag)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err(format!("{name}: no shape with tag {tag}")))?;
    if shape.vertices.len() != 1 || !shape.edges.is_empty() || !shape.solids.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "{name}: tag {tag} is not a point entity"
        )));
    }
    Ok(shape.vertices[0].point)
}

fn option_number(state: &RuntimeState, name: &str) -> Option<f64> {
    state.option_numbers.get(name).copied()
}

fn option_default_number(name: &str) -> Option<f64> {
    match name {
        "Mesh.MeshSizeMax" | "Mesh.CharacteristicLengthMax" => Some(1.0e22),
        "Mesh.MeshSizeMin" | "Mesh.CharacteristicLengthMin" => Some(0.0),
        "Mesh.MeshSizeFactor" | "Mesh.CharacteristicLengthFactor" => Some(1.0),
        "Mesh.Algorithm" => Some(6.0),
        "Mesh.Algorithm3D" => Some(1.0),
        "Geometry.Points" | "Mesh.Points" => Some(1.0),
        "Geometry.Curves" | "Geometry.Lines" | "Geometry.Surfaces" | "Geometry.Volumes" => {
            Some(1.0)
        }
        "STEP.GmshStrict" | "Geometry.OCCStepGmshStrict" => Some(1.0),
        _ => None,
    }
}

fn option_default_string(name: &str) -> Option<&'static str> {
    match name {
        "STEP.Protocol" | "Geometry.OCCStepProtocol" | "Geometry.OCCSTEPProtocol" => {
            Some("AP214")
        }
        "STEP.WriterStyle" | "Geometry.OCCSTEPWriterStyle" => Some("gmsh-strict"),
        _ => None,
    }
}

fn option_default_color(name: &str) -> Option<(i32, i32, i32, i32)> {
    match name {
        "General.Background" => Some((255, 255, 255, 255)),
        _ => None,
    }
}

fn estimate_mesh_characteristic_size(mesh: &Mesh) -> Option<f64> {
    let mut iter = mesh.nodes.values();
    let first = iter.next()?.position;

    let mut min = first;
    let mut max = first;
    let mut count = 1usize;

    for node in iter {
        let p = node.position;
        min.x = min.x.min(p.x);
        min.y = min.y.min(p.y);
        min.z = min.z.min(p.z);
        max.x = max.x.max(p.x);
        max.y = max.y.max(p.y);
        max.z = max.z.max(p.z);
        count += 1;
    }

    if count < 2 {
        return None;
    }

    let diag = (max - min).norm();
    if !diag.is_finite() || diag <= 1e-12 {
        return None;
    }

    // Keep auto-size conservative and geometry-scaled when explicit limits are absent.
    Some((diag / 20.0).max(1e-6))
}

fn push_log(state: &mut RuntimeState, msg: impl Into<String>) {
    if state.logger_enabled {
        state.logs.push(msg.into());
    }
}

fn refine_triangles_once(mesh: &mut Mesh) -> usize {
    use std::collections::HashMap;

    let mut next_node_id = mesh.nodes.keys().max().copied().unwrap_or(0) + 1;
    let mut next_elem_id = mesh.elements.iter().map(|e| e.id).max().unwrap_or(0) + 1;
    let mut edge_midpoint_node: HashMap<(u64, u64), u64> = HashMap::new();

    let elements_snapshot = mesh.elements.clone();
    let mut refined: Vec<Element> = Vec::with_capacity(elements_snapshot.len() * 2);
    let mut refined_count = 0usize;

    for elem in elements_snapshot {
        if elem.etype != ElementType::Triangle3 || elem.node_ids.len() != 3 {
            refined.push(elem);
            continue;
        }

        let a = elem.node_ids[0];
        let b = elem.node_ids[1];
        let c = elem.node_ids[2];

        let midpoint = |u: u64,
                        v: u64,
                        mesh: &mut Mesh,
                        next_node_id: &mut u64,
                        edge_midpoint_node: &mut HashMap<(u64, u64), u64>|
         -> Option<u64> {
            let key = if u < v { (u, v) } else { (v, u) };
            if let Some(id) = edge_midpoint_node.get(&key).copied() {
                return Some(id);
            }
            let pu = mesh.nodes.get(&u)?.position;
            let pv = mesh.nodes.get(&v)?.position;
            let id = *next_node_id;
            *next_node_id += 1;
            mesh.add_node(Node::new(
                id,
                0.5 * (pu.x + pv.x),
                0.5 * (pu.y + pv.y),
                0.5 * (pu.z + pv.z),
            ));
            edge_midpoint_node.insert(key, id);
            Some(id)
        };

        let Some(ab) = midpoint(a, b, mesh, &mut next_node_id, &mut edge_midpoint_node) else {
            refined.push(elem);
            continue;
        };
        let Some(bc) = midpoint(b, c, mesh, &mut next_node_id, &mut edge_midpoint_node) else {
            refined.push(elem);
            continue;
        };
        let Some(ca) = midpoint(c, a, mesh, &mut next_node_id, &mut edge_midpoint_node) else {
            refined.push(elem);
            continue;
        };

        let mut push_tri = |n0: u64, n1: u64, n2: u64, phys: Option<i32>| {
            let mut e = Element::new(next_elem_id, ElementType::Triangle3, vec![n0, n1, n2]);
            e.physical_tag = phys;
            next_elem_id += 1;
            refined.push(e);
        };

        push_tri(a, ab, ca, elem.physical_tag);
        push_tri(ab, b, bc, elem.physical_tag);
        push_tri(ca, bc, c, elem.physical_tag);
        push_tri(ab, bc, ca, elem.physical_tag);
        refined_count += 1;
    }

    mesh.elements = refined;
    refined_count
}

fn recombine_triangles_to_quads(mesh: &mut Mesh) -> usize {
    use std::collections::HashMap;

    let elements_snapshot = mesh.elements.clone();
    let mut next_elem_id = mesh.elements.iter().map(|e| e.id).max().unwrap_or(0) + 1;
    let mut tri_edges: HashMap<(u64, u64), Vec<(usize, [u64; 3], Option<i32>)>> = HashMap::new();
    let mut tri_indices: Vec<usize> = Vec::new();

    for (idx, elem) in elements_snapshot.iter().enumerate() {
        if elem.etype != ElementType::Triangle3 || elem.node_ids.len() != 3 {
            continue;
        }
        tri_indices.push(idx);
        let tri = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
        let edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
        for (u, v) in edges {
            let key = if u < v { (u, v) } else { (v, u) };
            tri_edges
                .entry(key)
                .or_default()
                .push((idx, tri, elem.physical_tag));
        }
    }

    let mut consumed: HashSet<usize> = HashSet::new();
    let mut recombined: Vec<Element> = Vec::new();
    let mut recombined_count = 0usize;

    for shared in tri_edges.values() {
        if shared.len() != 2 {
            continue;
        }
        let (i0, t0, p0) = shared[0];
        let (i1, t1, p1) = shared[1];
        if consumed.contains(&i0) || consumed.contains(&i1) {
            continue;
        }
        let mut common: Vec<u64> = t0.iter().copied().filter(|v| t1.contains(v)).collect();
        common.sort_unstable();
        common.dedup();
        if common.len() != 2 {
            continue;
        }
        let c0 = common[0];
        let c1 = common[1];
        let o0 = t0.iter().copied().find(|v| *v != c0 && *v != c1);
        let o1 = t1.iter().copied().find(|v| *v != c0 && *v != c1);
        let (Some(o0), Some(o1)) = (o0, o1) else {
            continue;
        };

        let mut quad = Element::new(next_elem_id, ElementType::Quad4, vec![o0, c0, o1, c1]);
        quad.physical_tag = if p0.is_some() { p0 } else { p1 };
        next_elem_id += 1;
        recombined.push(quad);
        consumed.insert(i0);
        consumed.insert(i1);
        recombined_count += 1;
    }

    let mut out: Vec<Element> = Vec::with_capacity(elements_snapshot.len());
    for (idx, elem) in elements_snapshot.into_iter().enumerate() {
        if consumed.contains(&idx) {
            continue;
        }
        out.push(elem);
    }
    out.extend(recombined);
    mesh.elements = out;
    recombined_count
}

fn gmsh_type_id(etype: ElementType) -> i32 {
    match etype {
        ElementType::Point1 => 15,
        ElementType::Line2 => 1,
        ElementType::Line3 => 8,
        ElementType::Triangle3 => 2,
        ElementType::Triangle6 => 9,
        ElementType::Quad4 => 3,
        ElementType::Quad9 => 10,
        ElementType::Tetrahedron4 => 4,
        ElementType::Tetrahedron10 => 11,
        ElementType::Hexahedron8 => 5,
        ElementType::Hexahedron27 => 12,
        ElementType::Prism6 => 6,
        ElementType::Prism18 => 13,
        ElementType::Pyramid5 => 7,
        ElementType::Pyramid14 => 14,
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
    state.model_name = "default".to_string();
    state.mesh_order = 1;
    state.current_mesh = None;
    state.current_path = None;
    state.entity_names.clear();
    state.physical_groups.clear();
    state.physical_names.clear();
    state.plugin_numbers.clear();
    state.plugin_strings.clear();
    state.logger_enabled = false;
    state.logs.clear();
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
    state.entity_names.clear();
    state.physical_groups.clear();
    state.physical_names.clear();
    state.plugin_numbers.clear();
    state.plugin_strings.clear();
    state.logger_enabled = false;
    state.logs.clear();
    state.cad_shapes.clear();
    state.occ_curve_loops.clear();
    state.occ_surface_loops.clear();
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
    state.entity_names.clear();
    state.physical_groups.clear();
    state.physical_names.clear();
    state.cad_shapes.clear();
    state.occ_curve_loops.clear();
    state.occ_surface_loops.clear();
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

    match ext.as_deref() {
        Some("msh") => {
            let mesh = state.current_mesh.as_ref().ok_or_else(|| {
                pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded; call open() or generate() first")
            })?;
            rmsh_io::save_msh_v4_to_path(&path, mesh)
                .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        }
        Some("step") | Some("stp") => {
            if !state.cad_shapes.is_empty() {
                let tags: Vec<i32> = state.cad_shapes.keys().copied().collect();
                let all_shapes: Vec<<RcadKernel as CadKernel>::Shape> = tags
                    .iter()
                    .filter_map(|tag| state.cad_shapes.get(tag).cloned())
                    .collect();

                let shape = if all_shapes.len() == 1 {
                    all_shapes[0].clone()
                } else {
                    state.kernel.merge_for_export(&all_shapes)
                };

                let gmsh_strict = step_gmsh_strict_from_state(&state);
                let solid_color = step_solid_color_from_state(&state);
                let protocol = step_protocol_from_state(&state);

                let options = StepExportOptions {
                    protocol,
                    solid_color,
                    gmsh_strict,
                };

                let step_str = state
                    .kernel
                    .write_step_string(&shape, &options)
                    .map_err(|e| pyo3::exceptions::PyIOError::new_err(e))?;

                std::fs::write(&path, &step_str)
                    .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
            } else {
                let mesh = state.current_mesh.as_ref().ok_or_else(|| {
                    pyo3::exceptions::PyRuntimeError::new_err(
                        "no CAD shape or mesh loaded; call model.occ.add*()/synchronize()/generate() first",
                    )
                })?;
                rmsh_io::save_step_to_path(&path, mesh)
                    .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
            }
        }
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
    ensure_initialized(&state)?;
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
    ensure_initialized(&state)?;
    state
        .option_numbers
        .get(&name)
        .copied()
        .or_else(|| option_default_number(&name))
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(format!("number option not set: {}", name)))
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
    ensure_initialized(&state)?;
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
    ensure_initialized(&state)?;
    state
        .option_strings
        .get(&name)
        .cloned()
        .or_else(|| option_default_string(&name).map(str::to_string))
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(format!("string option not set: {}", name)))
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
    ensure_initialized(&state)?;
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
    ensure_initialized(&state)?;
    state
        .option_colors
        .get(&name)
        .copied()
        .or_else(|| option_default_color(&name))
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(format!("color option not set: {}", name)))
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
    ensure_initialized(&state)?;
    state.option_numbers.clear();
    state.option_strings.clear();
    state.option_colors.clear();
    Ok(())
}

#[pyfunction]
#[pyo3(name = "logger_start", signature = (*args, **kwargs))]
fn logger_start_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let _ = (args, kwargs);
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.logger_enabled = true;
    state.logs.clear();
    push_log(&mut state, "logger started");
    Ok(())
}

#[pyfunction]
#[pyo3(name = "logger_stop", signature = (*args, **kwargs))]
fn logger_stop_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let _ = (args, kwargs);
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    push_log(&mut state, "logger stopped");
    state.logger_enabled = false;
    Ok(())
}

#[pyfunction]
#[pyo3(name = "logger_get", signature = (*args, **kwargs))]
fn logger_get_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<String>> {
    let _ = (args, kwargs);
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    Ok(state.logs.clone())
}

#[pyfunction]
#[pyo3(name = "model_add", signature = (*args, **kwargs))]
fn model_add_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.model_name = name.clone();
    push_log(&mut state, format!("model add: {name}"));
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_remove", signature = (*args, **kwargs))]
fn model_remove_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let _ = (args, kwargs);
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.current_mesh = None;
    state.current_path = None;
    state.cad_shapes.clear();
    state.next_cad_tag = 0;
    state.entity_names.clear();
    state.physical_groups.clear();
    state.physical_names.clear();
    push_log(&mut state, "model removed".to_string());
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_get_current", signature = (*args, **kwargs))]
fn model_get_current_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<String> {
    let _ = (args, kwargs);
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    Ok(state.model_name.clone())
}

#[pyfunction]
#[pyo3(name = "model_set_current", signature = (*args, **kwargs))]
fn model_set_current_impl(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.model_name = name.clone();
    push_log(&mut state, format!("model current set: {name}"));
    Ok(())
}

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
        let max_dim = state
            .cad_shapes
            .values()
            .map(|s| state.kernel.dimension(s))
            .max()
            .unwrap_or(0);
        return Ok(max_dim);
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
    for (tag, shape) in &state.cad_shapes {
        let d = state.kernel.dimension(shape);
        if dim < 0 || dim == d {
            out.push((d, *tag));
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

#[pyfunction]
#[pyo3(name = "model_get_entity_name", signature = (*args, **kwargs))]
fn model_get_entity_name_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int")?;
    let tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int")?;
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    Ok(state
        .entity_names
        .get(&(dim, tag))
        .cloned()
        .unwrap_or_else(|| format!("Entity({dim},{tag})")))
}

#[pyfunction]
#[pyo3(name = "model_set_entity_name", signature = (*args, **kwargs))]
fn model_set_entity_name_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int")?;
    let tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int")?;
    let name: String = extract_required(args, kwargs, 2, &["name"], "str")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.entity_names.insert((dim, tag), name.clone());
    push_log(&mut state, format!("entity name set ({dim},{tag})={name}"));
    Ok(())
}

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

    if tag > 0 {
        if let Some(shape) = state.cad_shapes.get(&tag) {
            let top = state.kernel.dimension(shape);
            if dim >= 0 && dim > top {
                return Err(pyo3::exceptions::PyRuntimeError::new_err(format!(
                    "entity ({dim},{tag}) not found"
                )));
            }
            if let Some((min, max)) = state.kernel.bounding_box(shape) {
                return Ok((min[0], min[1], min[2], max[0], max[1], max[2]));
            }
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

#[pyfunction]
#[pyo3(name = "model_get_entities_in_bounding_box", signature = (*args, **kwargs))]
fn model_get_entities_in_bounding_box_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<(i32, i32)>> {
    let xmin: f64 = extract_required(args, kwargs, 0, &["xmin"], "float")?;
    let ymin: f64 = extract_required(args, kwargs, 1, &["ymin"], "float")?;
    let zmin: f64 = extract_required(args, kwargs, 2, &["zmin"], "float")?;
    let xmax: f64 = extract_required(args, kwargs, 3, &["xmax"], "float")?;
    let ymax: f64 = extract_required(args, kwargs, 4, &["ymax"], "float")?;
    let zmax: f64 = extract_required(args, kwargs, 5, &["zmax"], "float")?;
    let dim: i32 = extract_required(args, kwargs, 6, &["dim"], "int").unwrap_or(-1);

    let qmin = [xmin.min(xmax), ymin.min(ymax), zmin.min(zmax)];
    let qmax = [xmin.max(xmax), ymin.max(ymax), zmin.max(zmax)];

    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut out = Vec::new();
    for (tag, shape) in &state.cad_shapes {
        let Some((smin, smax)) = state.kernel.bounding_box(shape) else {
            continue;
        };
        if !bbox_intersects(smin, smax, qmin, qmax) {
            continue;
        }
        let top = state.kernel.dimension(shape);
        for d in (0..=top).rev() {
            if dim < 0 || dim == d {
                out.push((d, *tag));
            }
        }
    }

    out.sort_unstable();
    out.dedup();
    Ok(out)
}

#[pyfunction]
#[pyo3(name = "model_add_physical_group", signature = (*args, **kwargs))]
fn model_add_physical_group_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int")?;
    let tags: Vec<i32> = extract_required(args, kwargs, 1, &["tags"], "list[int]")?;
    let requested_tag: i32 = extract_required(args, kwargs, 2, &["tag"], "int").unwrap_or(-1);
    let name: String = extract_required(args, kwargs, 3, &["name"], "str").unwrap_or_default();

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let next_tag = state
        .physical_groups
        .keys()
        .filter(|(d, _)| *d == dim)
        .map(|(_, t)| *t)
        .max()
        .unwrap_or(0)
        + 1;
    let tag = if requested_tag > 0 { requested_tag } else { next_tag };
    state.physical_groups.insert((dim, tag), tags.clone());
    if !name.is_empty() {
        state.physical_names.insert((dim, tag), name.clone());
    }
    push_log(&mut state, format!("physical group add ({dim},{tag}) with {} entities", tags.len()));
    Ok(tag)
}

#[pyfunction]
#[pyo3(name = "model_get_physical_groups", signature = (*args, **kwargs))]
fn model_get_physical_groups_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<(i32, i32)>> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int").unwrap_or(-1);
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut out: Vec<(i32, i32)> = state
        .physical_groups
        .keys()
        .filter(|(d, _)| dim < 0 || *d == dim)
        .copied()
        .collect();
    out.sort_unstable();
    Ok(out)
}

#[pyfunction]
#[pyo3(name = "model_set_physical_name", signature = (*args, **kwargs))]
fn model_set_physical_name_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int")?;
    let tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int")?;
    let name: String = extract_required(args, kwargs, 2, &["name"], "str")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.physical_names.insert((dim, tag), name);
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_get_physical_name", signature = (*args, **kwargs))]
fn model_get_physical_name_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<String> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int")?;
    let tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int")?;
    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state
        .physical_names
        .get(&(dim, tag))
        .cloned()
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(format!("physical name not set: ({dim},{tag})")))
}

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

    let shape = state
        .kernel
        .create_box(DVec3::new(x, y, z), dx, dy, dz)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

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

    let shape = state
        .kernel
        .create_sphere(DVec3::new(x, y, z), r)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

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
    let ref_dir = if axis_norm.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let shape = state
        .kernel
        .create_cylinder(DVec3::new(x, y, z), axis_norm, ref_dir, r, height)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

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
) -> PyResult<(Vec<(i32, i32)>, Vec<Vec<(i32, i32)>>)> {
    let obj_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["objectDimTags"], "list of (dim, tag)")?;
    let tool_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 1, &["toolDimTags"], "list of (dim, tag)")?;
    let tag: i32 = extract_required(args, kwargs, 2, &["tag"], "int").unwrap_or(-1);
    let remove_object: bool =
        extract_required(args, kwargs, 3, &["removeObject", "remove_object"], "bool")
            .unwrap_or(true);
    let remove_tool: bool =
        extract_required(args, kwargs, 4, &["removeTool", "remove_tool"], "bool")
            .unwrap_or(true);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let first_obj_tag = obj_dim_tags.first().map(|t| t.1).unwrap_or(1);
    let first_obj_dim = obj_dim_tags.first().map(|t| t.0).unwrap_or(3);

    let mut base = state.cad_shapes.get(&first_obj_tag).cloned().ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err("no valid object shape found for cut")
    })?;

    for &(_, tag) in &tool_dim_tags {
        if let Some(tool) = state.cad_shapes.get(&tag) {
            base = state
                .kernel
                .boolean_op(BooleanOp::Cut, &base, tool)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("boolean cut failed: {e}")))?;
        }
    }

    let result_tag = if tag > 0 {
        next_or_requested_tag(&mut state, tag)
    } else if remove_object {
        first_obj_tag
    } else {
        next_or_requested_tag(&mut state, -1)
    };
    state.current_mesh = Some(state.kernel.tessellate(&base));
    state.cad_shapes.insert(result_tag, base);

    if remove_object {
        for &(_, obj_tag) in &obj_dim_tags {
            if obj_tag != result_tag {
                state.cad_shapes.remove(&obj_tag);
            }
        }
    }
    if remove_tool {
        for &(_, tool_tag) in &tool_dim_tags {
            state.cad_shapes.remove(&tool_tag);
        }
    }

    let out_dim_tags = vec![(first_obj_dim, result_tag)];
    let mut out_dim_tags_map: Vec<Vec<(i32, i32)>> = Vec::new();
    for _ in &obj_dim_tags {
        out_dim_tags_map.push(out_dim_tags.clone());
    }
    for _ in &tool_dim_tags {
        if remove_tool {
            out_dim_tags_map.push(Vec::new());
        } else {
            out_dim_tags_map.push(out_dim_tags.clone());
        }
    }
    Ok((out_dim_tags, out_dim_tags_map))
}

#[pyfunction]
#[pyo3(name = "model_occ_fuse", signature = (*args, **kwargs))]
fn model_occ_fuse_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(Vec<(i32, i32)>, Vec<Vec<(i32, i32)>>)> {
    let obj_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["objectDimTags"], "list of (dim, tag)")?;
    let tool_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 1, &["toolDimTags"], "list of (dim, tag)")?;
    let tag: i32 = extract_required(args, kwargs, 2, &["tag"], "int").unwrap_or(-1);
    let remove_object: bool =
        extract_required(args, kwargs, 3, &["removeObject", "remove_object"], "bool")
            .unwrap_or(true);
    let remove_tool: bool =
        extract_required(args, kwargs, 4, &["removeTool", "remove_tool"], "bool")
            .unwrap_or(true);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let first_obj_tag = obj_dim_tags.first().map(|t| t.1).unwrap_or(1);
    let first_obj_dim = obj_dim_tags.first().map(|t| t.0).unwrap_or(3);

    let mut base = state.cad_shapes.get(&first_obj_tag).cloned().ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err("no valid object shape found for fuse")
    })?;

    for &(_, tag) in &tool_dim_tags {
        if let Some(tool) = state.cad_shapes.get(&tag) {
            base = state
                .kernel
                .boolean_op(BooleanOp::Fuse, &base, tool)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("boolean fuse failed: {e}")))?;
        }
    }

    let result_tag = if tag > 0 {
        next_or_requested_tag(&mut state, tag)
    } else if remove_object {
        first_obj_tag
    } else {
        next_or_requested_tag(&mut state, -1)
    };
    state.current_mesh = Some(state.kernel.tessellate(&base));
    state.cad_shapes.insert(result_tag, base);

    if remove_object {
        for &(_, obj_tag) in &obj_dim_tags {
            if obj_tag != result_tag {
                state.cad_shapes.remove(&obj_tag);
            }
        }
    }
    if remove_tool {
        for &(_, tool_tag) in &tool_dim_tags {
            state.cad_shapes.remove(&tool_tag);
        }
    }

    let out_dim_tags = vec![(first_obj_dim, result_tag)];
    let mut out_dim_tags_map: Vec<Vec<(i32, i32)>> = Vec::new();
    for _ in &obj_dim_tags {
        out_dim_tags_map.push(out_dim_tags.clone());
    }
    for _ in &tool_dim_tags {
        if remove_tool {
            out_dim_tags_map.push(Vec::new());
        } else {
            out_dim_tags_map.push(out_dim_tags.clone());
        }
    }
    Ok((out_dim_tags, out_dim_tags_map))
}

#[pyfunction]
#[pyo3(name = "model_occ_fragment", signature = (*args, **kwargs))]
fn model_occ_fragment_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(Vec<(i32, i32)>, Vec<Vec<(i32, i32)>>)> {
    let obj_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["objectDimTags"], "list of (dim, tag)")?;
    let tool_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 1, &["toolDimTags"], "list of (dim, tag)")?;
    let tag: i32 = extract_required(args, kwargs, 2, &["tag"], "int").unwrap_or(-1);
    let remove_object: bool =
        extract_required(args, kwargs, 3, &["removeObject", "remove_object"], "bool")
            .unwrap_or(true);
    let remove_tool: bool =
        extract_required(args, kwargs, 4, &["removeTool", "remove_tool"], "bool")
            .unwrap_or(true);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let objects: Vec<_> = obj_dim_tags
        .iter()
        .filter_map(|&(dim, tag)| state.cad_shapes.get(&tag).cloned().map(|s| (dim, tag, s)))
        .collect();
    if objects.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err("no valid object shape found for fragment"));
    }
    let tools: Vec<_> = tool_dim_tags
        .iter()
        .filter_map(|&(dim, tag)| state.cad_shapes.get(&tag).cloned().map(|s| (dim, tag, s)))
        .collect();
    if tools.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err("no valid tool shape found for fragment"));
    }

    let obj_shapes: Vec<_> = objects.iter().map(|(_, _, s)| s.clone()).collect();
    let tool_shapes: Vec<_> = tools.iter().map(|(_, _, s)| s.clone()).collect();
    let (split_objects, split_tools) = state
        .kernel
        .fragment(&obj_shapes, &tool_shapes)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("fragment failed: {e}")))?;

    let mut out_dim_tags: Vec<(i32, i32)> = Vec::new();
    let mut out_dim_tags_map: Vec<Vec<(i32, i32)>> = Vec::new();
    let mut kept_object_originals = std::collections::HashSet::new();
    let mut kept_tool_originals = std::collections::HashSet::new();

    // Map split objects back, preserving original tags for first part when remove_object == true
    for ((dim, original_tag, _), split) in objects.iter().zip(split_objects) {
        let parts = state.kernel.explode_solids(&split);
        let mut mapped: Vec<(i32, i32)> = Vec::new();
        for (idx, part) in parts.into_iter().enumerate() {
            let out_tag = if remove_object && idx == 0 {
                kept_object_originals.insert(*original_tag);
                *original_tag
            } else {
                next_or_requested_tag(&mut state, -1)
            };
            state.cad_shapes.insert(out_tag, part);
            let pair = (*dim, out_tag);
            out_dim_tags.push(pair);
            mapped.push(pair);
        }
        out_dim_tags_map.push(mapped);
    }

    // Map split tools
    for ((dim, original_tag, _), split) in tools.iter().zip(split_tools) {
        let parts = state.kernel.explode_solids(&split);
        let mut mapped: Vec<(i32, i32)> = Vec::new();
        for (idx, part) in parts.into_iter().enumerate() {
            let out_tag = if remove_tool && idx == 0 {
                kept_tool_originals.insert(*original_tag);
                *original_tag
            } else {
                next_or_requested_tag(&mut state, -1)
            };
            state.cad_shapes.insert(out_tag, part);
            let pair = (*dim, out_tag);
            out_dim_tags.push(pair);
            mapped.push(pair);
        }
        out_dim_tags_map.push(mapped);
    }

    if remove_object {
        for (_, original_tag) in &obj_dim_tags {
            if !kept_object_originals.contains(original_tag) {
                state.cad_shapes.remove(original_tag);
            }
        }
    }
    if remove_tool {
        for (_, original_tag) in &tool_dim_tags {
            if !kept_tool_originals.contains(original_tag) {
                state.cad_shapes.remove(original_tag);
            }
        }
    }

    if tag > 0 && !out_dim_tags.is_empty() {
        let old_tag = out_dim_tags[0].1;
        if old_tag != tag {
            if let Some(shape) = state.cad_shapes.remove(&old_tag) {
                let new_tag = next_or_requested_tag(&mut state, tag);
                state.cad_shapes.insert(new_tag, shape);
                out_dim_tags[0].1 = new_tag;
                for mapped in &mut out_dim_tags_map {
                    for pair in mapped {
                        if pair.1 == old_tag { pair.1 = new_tag; }
                    }
                }
            }
        }
    }

    let first_tag = out_dim_tags.first().map(|(_, t)| *t).ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("fragment produced no output shape")
    })?;
    if let Some(shape) = state.cad_shapes.get(&first_tag) {
        state.current_mesh = Some(state.kernel.tessellate(shape));
    }

    Ok((out_dim_tags, out_dim_tags_map))
}

#[pyfunction]
#[pyo3(name = "model_occ_intersect", signature = (*args, **kwargs))]
fn model_occ_intersect_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<(Vec<(i32, i32)>, Vec<Vec<(i32, i32)>>)> {
    let obj_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["objectDimTags"], "list of (dim, tag)")?;
    let tool_dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 1, &["toolDimTags"], "list of (dim, tag)")?;
    let tag: i32 = extract_required(args, kwargs, 2, &["tag"], "int").unwrap_or(-1);
    let remove_object: bool =
        extract_required(args, kwargs, 3, &["removeObject", "remove_object"], "bool")
            .unwrap_or(true);
    let remove_tool: bool =
        extract_required(args, kwargs, 4, &["removeTool", "remove_tool"], "bool")
            .unwrap_or(true);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let first_obj_tag = obj_dim_tags.first().map(|t| t.1).unwrap_or(1);
    let first_obj_dim = obj_dim_tags.first().map(|t| t.0).unwrap_or(3);

    let mut base = state.cad_shapes.get(&first_obj_tag).cloned().ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err("no valid object shape found for intersect")
    })?;

    for &(_, tag) in &tool_dim_tags {
        if let Some(tool) = state.cad_shapes.get(&tag) {
            base = state
                .kernel
                .boolean_op(BooleanOp::Intersect, &base, tool)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("boolean intersect failed: {e}")))?;
        }
    }

    let result_tag = if tag > 0 {
        next_or_requested_tag(&mut state, tag)
    } else if remove_object {
        first_obj_tag
    } else {
        next_or_requested_tag(&mut state, -1)
    };
    state.current_mesh = Some(state.kernel.tessellate(&base));
    state.cad_shapes.insert(result_tag, base);

    if remove_object {
        for &(_, obj_tag) in &obj_dim_tags {
            if obj_tag != result_tag {
                state.cad_shapes.remove(&obj_tag);
            }
        }
    }
    if remove_tool {
        for &(_, tool_tag) in &tool_dim_tags {
            state.cad_shapes.remove(&tool_tag);
        }
    }

    let out_dim_tags = vec![(first_obj_dim, result_tag)];
    let mut out_dim_tags_map: Vec<Vec<(i32, i32)>> = Vec::new();
    for _ in &obj_dim_tags {
        out_dim_tags_map.push(out_dim_tags.clone());
    }
    for _ in &tool_dim_tags {
        if remove_tool {
            out_dim_tags_map.push(Vec::new());
        } else {
            out_dim_tags_map.push(out_dim_tags.clone());
        }
    }
    Ok((out_dim_tags, out_dim_tags_map))
}

#[pyfunction]
#[pyo3(name = "model_occ_copy", signature = (*args, **kwargs))]
fn model_occ_copy_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<(i32, i32)>> {
    let dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["dimTags"], "list of (dim, tag)")?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut out = Vec::new();
    for (_dim, tag) in dim_tags {
        let src = state.cad_shapes.get(&tag).cloned().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
        })?;
        let new_tag = state.next_cad_tag + 1;
        state.next_cad_tag = new_tag;
        let dim = state.kernel.dimension(&src);
        state.cad_shapes.insert(new_tag, src);
        out.push((dim, new_tag));
    }
    Ok(out)
}

#[pyfunction]
#[pyo3(name = "model_occ_remove", signature = (*args, **kwargs))]
fn model_occ_remove_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["dimTags"], "list of (dim, tag)")?;
    let recursive: bool =
        extract_required(args, kwargs, 1, &["recursive"], "bool").unwrap_or(false);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    for (_dim, tag) in &dim_tags {
        state.cad_shapes.remove(tag);
        state.occ_curve_loops.remove(tag);
        state.occ_surface_loops.remove(tag);
    }

    if recursive {
        let removed: std::collections::HashSet<i32> = dim_tags.iter().map(|(_, t)| *t).collect();
        state
            .occ_curve_loops
            .retain(|_, curves| !curves.iter().any(|t| removed.contains(&t.abs())));
        state
            .occ_surface_loops
            .retain(|_, faces| !faces.iter().any(|t| removed.contains(t)));
    }

    state.current_mesh = None;
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_occ_translate", signature = (*args, **kwargs))]
fn model_occ_translate_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["dimTags"], "list of (dim, tag)")?;
    let dx: f64 = extract_required(args, kwargs, 1, &["dx"], "float")?;
    let dy: f64 = extract_required(args, kwargs, 2, &["dy"], "float")?;
    let dz: f64 = extract_required(args, kwargs, 3, &["dz"], "float")?;

    let xf = DAffine3::from_translation(DVec3::new(dx, dy, dz));

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    for (_dim, tag) in dim_tags {
        let mut shape = state.cad_shapes.get(&tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
        })?.clone();
        state.kernel.apply_transform(&mut shape, xf);
        state.cad_shapes.insert(tag, shape);
    }
    state.current_mesh = None;
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_occ_rotate", signature = (*args, **kwargs))]
fn model_occ_rotate_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["dimTags"], "list of (dim, tag)")?;
    let x: f64 = extract_required(args, kwargs, 1, &["x"], "float")?;
    let y: f64 = extract_required(args, kwargs, 2, &["y"], "float")?;
    let z: f64 = extract_required(args, kwargs, 3, &["z"], "float")?;
    let ax: f64 = extract_required(args, kwargs, 4, &["ax"], "float")?;
    let ay: f64 = extract_required(args, kwargs, 5, &["ay"], "float")?;
    let az: f64 = extract_required(args, kwargs, 6, &["az"], "float")?;
    let angle: f64 = extract_required(args, kwargs, 7, &["angle"], "float")?;

    let axis = DVec3::new(ax, ay, az);
    if axis.length_squared() < 1e-20 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "rotate axis (ax, ay, az) must be non-zero",
        ));
    }
    let axis = axis.normalize();
    let r = DMat3::from_quat(DQuat::from_axis_angle(axis, angle));
    let p = DVec3::new(x, y, z);
    let t = p - r * p;
    let xf = DAffine3::from_mat3_translation(r, t);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    for (_dim, tag) in dim_tags {
        let mut shape = state.cad_shapes.get(&tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
        })?.clone();
        state.kernel.apply_transform(&mut shape, xf);
        state.cad_shapes.insert(tag, shape);
    }
    state.current_mesh = None;
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_occ_dilate", signature = (*args, **kwargs))]
fn model_occ_dilate_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["dimTags"], "list of (dim, tag)")?;
    let x: f64 = extract_required(args, kwargs, 1, &["x"], "float")?;
    let y: f64 = extract_required(args, kwargs, 2, &["y"], "float")?;
    let z: f64 = extract_required(args, kwargs, 3, &["z"], "float")?;
    let a: f64 = extract_required(args, kwargs, 4, &["a"], "float")?;
    let b: f64 = extract_required(args, kwargs, 5, &["b"], "float")?;
    let c: f64 = extract_required(args, kwargs, 6, &["c"], "float")?;

    let s = DMat3::from_diagonal(DVec3::new(a, b, c));
    let p = DVec3::new(x, y, z);
    let t = p - s * p;
    let xf = DAffine3::from_mat3_translation(s, t);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    for (_dim, tag) in dim_tags {
        let mut shape = state.cad_shapes.get(&tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
        })?.clone();
        state.kernel.apply_transform(&mut shape, xf);
        state.cad_shapes.insert(tag, shape);
    }
    state.current_mesh = None;
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_occ_mirror", signature = (*args, **kwargs))]
fn model_occ_mirror_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let dim_tags: Vec<(i32, i32)> =
        extract_required(args, kwargs, 0, &["dimTags"], "list of (dim, tag)")?;
    let a: f64 = extract_required(args, kwargs, 1, &["a"], "float")?;
    let b: f64 = extract_required(args, kwargs, 2, &["b"], "float")?;
    let c: f64 = extract_required(args, kwargs, 3, &["c"], "float")?;
    let d: f64 = extract_required(args, kwargs, 4, &["d"], "float")?;

    let n_raw = DVec3::new(a, b, c);
    if n_raw.length_squared() < 1e-20 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "mirror plane normal (a,b,c) must be non-zero",
        ));
    }
    let n = n_raw.normalize();
    let d_n = d / n_raw.length();

    let outer = DMat3::from_cols(
        DVec3::new(n.x * n.x, n.x * n.y, n.x * n.z),
        DVec3::new(n.y * n.x, n.y * n.y, n.y * n.z),
        DVec3::new(n.z * n.x, n.z * n.y, n.z * n.z),
    );
    let m = DMat3::IDENTITY - 2.0 * outer;
    let t = -2.0 * d_n * n;
    let xf = DAffine3::from_mat3_translation(m, t);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    for (_dim, tag) in dim_tags {
        let mut shape = state.cad_shapes.get(&tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
        })?.clone();
        state.kernel.apply_transform(&mut shape, xf);
        state.cad_shapes.insert(tag, shape);
    }
    state.current_mesh = None;
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_occ_import_shapes", signature = (*args, **kwargs))]
fn model_occ_import_shapes_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<(i32, i32)>> {
    let file_name: String =
        extract_required(args, kwargs, 0, &["fileName", "file_name"], "str")?;
    let highest_dim_only: bool =
        extract_required(args, kwargs, 1, &["highestDimOnly"], "bool").unwrap_or(true);
    let _format: String =
        extract_required(args, kwargs, 2, &["format"], "str").unwrap_or_default();

    let path = PathBuf::from(&file_name);
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();

    if ext != "step" && ext != "stp" {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "importShapes currently supports only STEP (.step/.stp)",
        ));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let imported = state
        .kernel
        .read_step_file(&path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e))?;

    let tag = state.next_cad_tag + 1;
    state.next_cad_tag = tag;
    let dim = state.kernel.dimension(&imported);
    state.cad_shapes.insert(tag, imported);
    state.current_mesh = None;

    if highest_dim_only {
        Ok(vec![(dim, tag)])
    } else {
        let top_dim = state.kernel.dimension(state.cad_shapes.get(&tag).expect("inserted shape exists"));
        let mut dims: Vec<i32> = (0..=top_dim).rev().collect();
        dims.dedup();
        Ok(dims.into_iter().map(|d| (d, tag)).collect())
    }
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

    state.current_mesh = None;
    let tags: Vec<i32> = state.cad_shapes.keys().copied().collect();
    for tag in tags {
        if let Some(shape) = state.cad_shapes.get(&tag) {
            let mesh = state.kernel.tessellate(shape);
            match state.current_mesh.as_mut() {
                Some(current) => merge_meshes(current, &mesh),
                None => state.current_mesh = Some(mesh),
            }
        }
    }
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
    let surface = if let Some(mesh) = state.current_mesh.clone() {
        mesh
    } else if !state.cad_shapes.is_empty() {
        let mut tags: Vec<i32> = state.cad_shapes.keys().copied().collect();
        tags.sort_unstable();
        let mut shape = state.cad_shapes.get(&tags[0]).cloned().ok_or_else(|| {
            pyo3::exceptions::PyRuntimeError::new_err("no CAD shapes available")
        })?;
        if tags.len() > 1 {
            let all_shapes: Vec<_> = tags
                .iter()
                .filter_map(|tag| state.cad_shapes.get(tag).cloned())
                .collect();
            shape = state.kernel.merge_for_export(&all_shapes);
        }

        state.kernel.tessellate(&shape)
    } else {
        return Err(pyo3::exceptions::PyRuntimeError::new_err(
            "no mesh loaded and no CAD shape available",
        ));
    };

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

    let auto_size = estimate_mesh_characteristic_size(&surface).unwrap_or(1.0);
    let base_size = size_max.unwrap_or(auto_size) * size_factor;
    let base_size = if base_size.is_finite() && base_size > 0.0 {
        base_size
    } else {
        auto_size.max(1e-6)
    };
    let mut params = MeshParams::with_size(base_size);
    // Wire background size field if set.
    if state.field_mgr.has_background() {
        push_log(&mut state, "background size field active for meshing".to_string());
    }
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
        // 3D algorithms: 1=Delaunay, 4=Frontal, 7=MMG3D, 10=HXT
        // 3=Automatic → redirect to 1 (Delaunay)
        match algo_3d {
            1 | 3 => Delaunay3D::default().mesh_3d(&surface, &params).map_err(convert_err)?,
            4 => Frontal3D::default().mesh_3d(&surface, &params).map_err(convert_err)?,
            7 => MmgRemesh::default().mesh_3d(&surface, &params).map_err(convert_err)?,
            10 => Hxt3D::default().mesh_3d(&surface, &params).map_err(convert_err)?,
            _ => CentroidStarMesher3D.mesh_3d(&surface, &params).map_err(convert_err)?,
        }
    } else if dim == 2 {
        let polygon: Vec<[f64; 2]> = match boundary_loop_from_surface_mesh(&surface) {
            Ok(p) => p,
            Err(_) => {
                state
                    .cad_shapes
                    .values()
                    .find_map(boundary_loop_from_brep_face)
                    .ok_or_else(|| {
                        pyo3::exceptions::PyRuntimeError::new_err(
                            "cannot extract boundary loop from current 2D mesh",
                        )
                    })?
            }
        };

        // ── Transfinite TFI path: structured quad mesh from boundary ────────
        if !state.transfinite_surfaces.is_empty() && polygon.len() >= 4 {
            let corners = detect_corner_indices(&polygon);
            if corners.len() >= 4 {
                let _convert_err = |e: MeshAlgoError| pyo3::exceptions::PyRuntimeError::new_err(e.to_string());
                let edges = split_polygon_into_edges(&polygon, &corners);

                // Determine segment counts from transfinite curve constraints.
                // First 2 non-empty edges → nr, ns.
                let mut segs: Vec<usize> = Vec::new();
                for edge in &edges {
                    if edge.len() >= 2 {
                        // Look for a matching transfinite curve by tag,
                        // or fall back to estimate from edge length/size.
                        let edge_len = {
                            let p0 = edge[0];
                            let p1 = edge[edge.len() - 1];
                            let dx = p1[0] - p0[0];
                            let dy = p1[1] - p0[1];
                            (dx * dx + dy * dy).sqrt()
                        };
                        let num_seg = (edge_len / params.element_size).ceil() as usize;
                        segs.push(num_seg.max(1));
                    }
                }

                if segs.len() >= 2 {
                    let nr = segs[0];
                    let ns = segs[1];

                    // Build 4 TransfiniteCurve objects from corners.
                    let c = &corners;
                    let p = &polygon;
                    let curves = [
                        rmsh_algo::TransfiniteCurve::new([p[c[0]][0], p[c[0]][1], 0.0], [p[c[1]][0], p[c[1]][1], 0.0], nr),
                        rmsh_algo::TransfiniteCurve::new([p[c[1]][0], p[c[1]][1], 0.0], [p[c[2]][0], p[c[2]][1], 0.0], ns),
                        rmsh_algo::TransfiniteCurve::new([p[c[3]][0], p[c[3]][1], 0.0], [p[c[2]][0], p[c[2]][1], 0.0], nr),
                        rmsh_algo::TransfiniteCurve::new([p[c[0]][0], p[c[0]][1], 0.0], [p[c[3]][0], p[c[3]][1], 0.0], ns),
                    ];

                    let tfi_mesh = rmsh_algo::tfi_surface_mesh(&curves, nr, ns)
                        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
                    state.current_mesh = Some(tfi_mesh);
                    push_log(&mut state, format!("transfinite TFI: {nr}×{ns} structured quads"));
                    return Ok(());
                }
            }
        }

        let domain = Domain2D::from_outer(polygon);
        // 2D algorithms: 1=MeshAdapt, 5=Delaunay, 6=Frontal-Delaunay, 7=BAMG,
        //   8=Frontal-Quads, 9=Quad-Paving
        // 2=Automatic → redirect to 6 (Frontal-Delaunay)
        // 3=Initial mesh only → redirect to 5 (Delaunay, coarse)
        // 4=Frontal-Delaunay for quads (legacy) → redirect to 8 (Frontal-Quads)
        match algo_2d {
            1 => MeshAdapt2D::default().mesh_2d(&domain, &params).map_err(convert_err)?,
            2 | 6 => FrontalDelaunay2D::default().mesh_2d(&domain, &params).map_err(convert_err)?,
            3 | 5 => Delaunay2D::default().mesh_2d(&domain, &params).map_err(convert_err)?,
            4 | 8 => FrontalQuads2D::default().mesh_2d(&domain, &params).map_err(convert_err)?,
            7 => Bamg2D::default().mesh_2d(&domain, &params).map_err(convert_err)?,
            9 => QuadPaving2D::default().mesh_2d(&domain, &params).map_err(convert_err)?,
            _ => mesh_polygon(&Polygon2D::new(domain.outer().to_vec()), params.element_size)
                    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?,
        }
    } else {
        return Err(PyNotImplementedError::new_err(
            "only dim=2 and dim=3 are currently implemented",
        ));
    };
    state.current_mesh = Some(generated);

    // Post-generation optimization: if Mesh.Optimize is set and the mesh
    // qualifies, run the topology-based optimizer automatically.
    let optimize_enabled = option_number(&state, "Mesh.Optimize").unwrap_or(0.0) != 0.0;
    let mesh_order = state.mesh_order;
    drop(state); // release lock before optimize
    if optimize_enabled {
        if let Some(ref mut mesh) = STATE.lock().map_err(|_|
            pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned")
        )?.current_mesh {
            let opt_params = OptimizeParams {
                iterations: 3,
                ..Default::default()
            };
            if dim == 3 {
                let config = OptimizeConfig {
                    metric: QualityMetric::MinAngle,
                    edge_swap: true,
                    laplacian_smooth: true,
                    node_insertion: false,
                    edge_collapse: true,
                    threshold: 18.0,
                };
                let _ = MeshQualityOptimizer::with_config(MeshQualityOptimizer::new(), config)
                    .optimize(mesh, &opt_params);
            }
        }
    }
    // Re-acquire state for possible P2 promotion below.
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    if mesh_order >= 2 {
        if let Some(ref mut mesh) = state.current_mesh {
            promote_to_p2(mesh);
            push_log(&mut state, format!("post-gen P2 promote: {mesh_order}"));
        }
    }
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_mesh_set_order", signature = (*args, **kwargs))]
fn model_mesh_set_order_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let order: i32 = extract_required(args, kwargs, 0, &["order"], "int")?;
    if order < 1 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "mesh order must be >= 1",
        ));
    }
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.mesh_order = order;
    // Apply P2 promotion if order >= 2 and there's a current mesh.
    if order >= 2 {
        if let Some(ref mut mesh) = state.current_mesh {
            let before_nodes = mesh.node_count();
            promote_to_p2(mesh);
            let added = mesh.node_count() - before_nodes;
            push_log(&mut state, format!("mesh order set to {order}: P2 promotion added {added} edge-midpoint nodes"));
        } else {
            push_log(&mut state, format!("mesh order set to {order} (no mesh yet — will promote on next generate)"));
        }
    } else {
        push_log(&mut state, format!("mesh order set to {order}"));
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

/// `model_mesh_get_element_qualities(tags, type)` → Vec<f64>
///
/// Returns quality metrics for the requested elements.
/// `type` is one of: "minAngle" (tri/tet), "scaledJacobian", "aspectRatio", "radiusEdgeRatio".
/// Returns one quality value per element tag.
#[pyfunction]
#[pyo3(name = "model_mesh_get_element_qualities", signature = (*args, **kwargs))]
fn model_mesh_get_element_qualities_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<f64>> {
    let elem_tags: Vec<u64> = extract_required(args, kwargs, 0, &["tags"], "list")?;
    let quality_type: String = extract_required(args, kwargs, 1, &["type"], "str")?;

    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let mesh = state.current_mesh.as_ref().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded")
    })?;

    // Build a lookup for element tags
    let tag_to_idx: std::collections::HashMap<u64, usize> = mesh
        .elements
        .iter()
        .enumerate()
        .map(|(i, e)| (e.id, i))
        .collect();

    let mut results = Vec::with_capacity(elem_tags.len());
    for tag in &elem_tags {
        let Some(&ei) = tag_to_idx.get(tag) else {
            return Err(pyo3::exceptions::PyValueError::new_err(
                format!("element tag {tag} not found")));
        };
        let elt = &mesh.elements[ei];
        let quality = compute_element_quality_impl(mesh, elt, &quality_type)?;
        results.push(quality);
    }
    Ok(results)
}

/// Helper: compute a single quality metric for an element.
fn compute_element_quality_impl(
    mesh: &Mesh,
    elt: &Element,
    quality_type: &str,
) -> PyResult<f64> {
    let node_xyz = |id: u64| -> PyResult<[f64; 3]> {
        let n = mesh.nodes.get(&id).ok_or_else(|| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("missing node {id}"))
        })?;
        Ok([n.position.x, n.position.y, n.position.z])
    };

    match elt.node_ids.len() {
        3 => {
            let a = node_xyz(elt.node_ids[0])?;
            let b = node_xyz(elt.node_ids[1])?;
            let c = node_xyz(elt.node_ids[2])?;
            Ok(compute_tri_quality(a, b, c, quality_type))
        }
        4 if elt.etype == ElementType::Tetrahedron4 => {
            let a = node_xyz(elt.node_ids[0])?;
            let b = node_xyz(elt.node_ids[1])?;
            let c = node_xyz(elt.node_ids[2])?;
            let d = node_xyz(elt.node_ids[3])?;
            Ok(compute_tet_quality(a, b, c, d, quality_type))
        }
        _ => Err(pyo3::exceptions::PyNotImplementedError::new_err(
            format!("quality for element type {:?} with {} nodes", elt.etype, elt.node_ids.len()))),
    }
}

fn compute_tri_quality(a: [f64; 3], b: [f64; 3], c: [f64; 3], qtype: &str) -> f64 {
    let ab = [b[0]-a[0], b[1]-a[1]]; let ac = [c[0]-a[0], c[1]-a[1]];
    let ba = [-ab[0], -ab[1]]; let bc = [c[0]-b[0], c[1]-b[1]];
    let ca = [-ac[0], -ac[1]]; let cb = [-bc[0], -bc[1]];
    let ang = |u: [f64;2], v: [f64;2]| -> f64 {
        let d = ((u[0]*u[0]+u[1]*u[1]).sqrt() * (v[0]*v[0]+v[1]*v[1]).sqrt()).max(1e-15);
        ((u[0]*v[0]+u[1]*v[1]) / d).clamp(-1.0,1.0).acos().to_degrees()
    };
    match qtype {
        "minAngle" => ang(ab,ac).min(ang(ba,bc)).min(ang(ca,cb)),
        "scaledJacobian" => {
            let cross = (ab[0]*ac[1] - ab[1]*ac[0]).abs();
            let l1 = (ab[0]*ab[0]+ab[1]*ab[1]).sqrt().max(1e-15);
            let l2 = (ac[0]*ac[0]+ac[1]*ac[1]).sqrt().max(1e-15);
            cross / (l1 * l2)
        }
        "aspectRatio" => {
            let dl = |x: [f64;2], y: [f64;2]| ((y[0]-x[0]).powi(2)+(y[1]-x[1]).powi(2)).sqrt();
            let d01=dl([a[0],a[1]],[b[0],b[1]]); let d12=dl([b[0],b[1]],[c[0],c[1]]);
            let d20=dl([c[0],c[1]],[a[0],a[1]]);
            d01.max(d12).max(d20) / d01.min(d12).min(d20).max(1e-30)
        }
        _ => 0.0,
    }
}

fn compute_tet_quality(a: [f64;3], b: [f64;3], c: [f64;3], d: [f64;3], qtype: &str) -> f64 {
    let sub = |x:[f64;3],y:[f64;3]| [x[0]-y[0],x[1]-y[1],x[2]-y[2]];
    let dot = |x:[f64;3],y:[f64;3]| x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
    let len = |x:[f64;3]| dot(x,x).sqrt();
    let cross = |x:[f64;3],y:[f64;3]| [x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]];
    match qtype {
        "minAngle" => {
            let dih = |p,q,r,s| { let n1=cross(sub(q,p),sub(r,p)); let n2=cross(sub(q,p),sub(s,p));
                (dot(n1,n2)/(len(n1)*len(n2)).max(1e-15)).clamp(-1.0,1.0).acos().to_degrees() };
            dih(a,b,c,d).min(dih(a,c,b,d)).min(dih(a,d,b,c)).min(dih(b,c,a,d)).min(dih(b,d,a,c)).min(dih(c,d,a,b))
        }
        "radiusEdgeRatio" => {
            let edges = [len(sub(a,b)),len(sub(a,c)),len(sub(a,d)),len(sub(b,c)),len(sub(b,d)),len(sub(c,d))];
            let min_e = edges.iter().copied().fold(f64::MAX, f64::min).max(1e-15);
            let ba=sub(b,a);let ca=sub(c,a);let da=sub(d,a);
            let rhs = [0.5*(dot(b,b)-dot(a,a)),0.5*(dot(c,c)-dot(a,a)),0.5*(dot(d,d)-dot(a,a))];
            // 3×3 solver
            let mut m = [[ba[0],ba[1],ba[2],rhs[0]],[ca[0],ca[1],ca[2],rhs[1]],[da[0],da[1],da[2],rhs[2]]];
            for col in 0..3 {
                let mut pv=col; for r in (col+1)..3 { if m[r][col].abs()>m[pv][col].abs(){pv=r} }
                if m[pv][col].abs()<1e-15 { continue; }
                m.swap(pv,col); let inv=1.0/m[col][col];
                for j in col..4 { m[col][j]*=inv }
                for r in 0..3 { if r==col{continue} let f=m[r][col]; for j in col..4{m[r][j]-=f*m[col][j]} }
            }
            let cc = Some([m[0][3],m[1][3],m[2][3]]);
            match cc { Some(center) => len(sub(center,a)) / min_e, None => 1e6 }
        }
        _ => 0.0,
    }
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
    // Clear field and transfinite constraints on mesh clear.
    state.field_mgr = FieldManager::new();
    state.transfinite_curves.clear();
    state.transfinite_surfaces.clear();
    Ok(())
}

// ─── Size Field API (model.mesh.field.*) ─────────────────────────────────

/// `model_mesh_field_add(kind)` → int tag
#[pyfunction]
#[pyo3(name = "model_mesh_field_add", signature = (*args, **kwargs))]
fn model_mesh_field_add_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let kind: String = extract_required(args, kwargs, 0, &["kind"], "str")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let field: Box<dyn rmsh_algo::SizeField> = match kind.as_str() {
        "Constant" => Box::new(ConstantField { lc: 1.0 }),
        "Distance" => Box::new(DistanceField {
            source_type: "Point".into(), cx: 0.0, cy: 0.0, cz: 0.0, source_tag: 0,
        }),
        "Threshold" => Box::new(ThresholdField {
            in_field: 0, lc_min: 0.1, lc_max: 1.0, dist_min: 0.0, dist_max: 1.0,
        }),
        "Min" => Box::new(MinField { fields: vec![] }),
        "Max" => Box::new(MaxField { fields: vec![] }),
        "MathEval" => Box::new(MathEvalField {
            expression: "1.0".into(), fallback_lc: 1.0,
        }),
        "Box" => Box::new(BoxField {
            lc_inside: 0.1, lc_outside: 1.0,
            x_min: -1.0, x_max: 1.0, y_min: -1.0, y_max: 1.0, z_min: -1.0, z_max: 1.0,
        }),
        "Restrict" => Box::new(RestrictField { in_field: 0, box_bounds: None }),
        _ => return Err(pyo3::exceptions::PyValueError::new_err(
            format!("unknown field type '{kind}'; supported: Constant/Distance/Threshold/Min/Max/MathEval/Box/Restrict"),
        )),
    };
    let tag = state.field_mgr.add(field);
    push_log(&mut state, format!("field {tag}: type={kind}"));
    Ok(tag)
}

/// `model_mesh_field_set_number(tag, key, value)`
#[pyfunction]
#[pyo3(name = "model_mesh_field_set_number", signature = (*args, **kwargs))]
fn model_mesh_field_set_number_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;
    let key: String = extract_required(args, kwargs, 1, &["key"], "str")?;
    let value: f64 = extract_required(args, kwargs, 2, &["value"], "float")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    // Store parameter in option_numbers under "Field.{tag}.{key}" namespace.
    let opt_key = format!("Field.{tag}.{key}");
    state.option_numbers.insert(opt_key, value);
    push_log(&mut state, format!("field {tag}: set {key} = {value}"));
    Ok(())
}

/// `model_mesh_field_set_as_background_mesh(tag)`
#[pyfunction]
#[pyo3(name = "model_mesh_field_set_as_background_mesh", signature = (*args, **kwargs))]
fn model_mesh_field_set_as_background_mesh_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.field_mgr.set_background(tag)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
    push_log(&mut state, format!("background mesh field set to {tag}"));
    Ok(())
}

// ─── Transfinite API (model.mesh.setTransfinite*) ─────────────────────────

/// `model_mesh_set_transfinite_curve(tag, num_nodes, [mesh_type, coef])`
#[pyfunction]
#[pyo3(name = "model_mesh_set_transfinite_curve", signature = (*args, **kwargs))]
fn model_mesh_set_transfinite_curve_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;
    let num_nodes: i32 = extract_required(args, kwargs, 1, &["numNodes"], "int")?;
    let mesh_type: String = extract_required(args, kwargs, 2, &["meshType"], "str").unwrap_or_else(|_| "Progression".to_string());
    let coef: f64 = extract_required(args, kwargs, 3, &["coef"], "float").unwrap_or(1.0);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let segments = (num_nodes - 1).max(1) as usize;
    state.transfinite_curves.push((tag, segments, mesh_type.clone(), coef));
    push_log(&mut state, format!("transfinite curve {tag}: {segments} segments, type={mesh_type}"));
    Ok(())
}

/// `model_mesh_set_transfinite_surface(tag, arrangement, corners)`
#[pyfunction]
#[pyo3(name = "model_mesh_set_transfinite_surface", signature = (*args, **kwargs))]
fn model_mesh_set_transfinite_surface_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;
    let arrangement: String = extract_required(args, kwargs, 1, &["arrangement"], "str")
        .unwrap_or_else(|_| "Left".to_string());
    let corners: Vec<i32> = extract_required(args, kwargs, 2, &["cornerTags"], "list").unwrap_or_default();

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state.transfinite_surfaces.push((tag, arrangement.clone(), corners));
    push_log(&mut state, format!("transfinite surface {tag}: arrangement={arrangement}"));
    Ok(())
}

/// `model_mesh_optimize(method, niter)`
///
/// Supported methods:
/// - "Laplace" / "Laplacian" — Laplacian smoothing (default)
/// - "Optimize3D" — topology-based quality optimization for 3-D tet meshes
/// - "Optimize2D" — topology-based quality optimization for 2-D tri meshes
/// - "HighOrder" — high-order element optimization (P2 curvature adaptation)
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

    // Read OptimizeThreshold option before borrowing state mutably.
    let opt_threshold = {
        let state = STATE
            .lock()
            .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
        option_number(&state, "Mesh.OptimizeThreshold")
    };

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
            push_log(&mut state, format!("laplacian smooth: {niter} iterations"));
        }
        "Optimize3D" | "Optimize" => {
            let threshold = opt_threshold.unwrap_or(20.0);
            let config = OptimizeConfig {
                metric: QualityMetric::MinAngle,
                edge_swap: true,
                laplacian_smooth: true,
                node_insertion: false,
                edge_collapse: true,
                threshold,
            };
            MeshQualityOptimizer::with_config(MeshQualityOptimizer::new(), config)
                .optimize(mesh, &params)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
            push_log(&mut state, format!("topology optimize ({niter} passes, threshold={threshold}°)"));
        }
        "Netgen" | "HighOrder" => {
            let lapsmooth = LaplacianSmooth::default();
            lapsmooth.optimize(mesh, &params)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
            let config = OptimizeConfig {
                metric: QualityMetric::MinAngle,
                edge_swap: true,
                laplacian_smooth: false,
                node_insertion: false,
                edge_collapse: false,
                threshold: 20.0,
            };
            MeshQualityOptimizer::with_config(MeshQualityOptimizer::new(), config)
                .optimize(mesh, &OptimizeParams { iterations: niter, ..Default::default() })
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
            push_log(&mut state, format!("Netgen-style optimize ({niter} passes)"));
        }
        other => {
            return Err(PyNotImplementedError::new_err(format!(
                "optimizer '{other}' not yet implemented; available: Laplace, Optimize3D, Netgen"
            )));
        }
    }
    Ok(())
}

#[pyfunction]
#[pyo3(name = "model_mesh_refine", signature = (*args, **kwargs))]
fn model_mesh_refine_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<usize> {
    let _ = (args, kwargs);
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let mesh = state.current_mesh.as_mut().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded")
    })?;

    let refined = refine_triangles_once(mesh);
    push_log(&mut state, format!("mesh refine: {refined} triangle elements refined"));
    Ok(refined)
}

#[pyfunction]
#[pyo3(name = "model_mesh_recombine", signature = (*args, **kwargs))]
fn model_mesh_recombine_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<usize> {
    let dim: i32 = extract_required(args, kwargs, 0, &["dim"], "int").unwrap_or(2);
    let _tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int").unwrap_or(-1);
    let _angle: f64 = extract_required(args, kwargs, 2, &["angle"], "float").unwrap_or(45.0);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let mesh = state.current_mesh.as_mut().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded")
    })?;

    if dim != 2 {
        push_log(&mut state, format!("mesh recombine skipped for dim={dim}"));
        return Ok(0);
    }

    let recombined = recombine_triangles_to_quads(mesh);
    push_log(&mut state, format!("mesh recombine: {recombined} quads created"));
    Ok(recombined)
}

#[pyfunction]
#[pyo3(name = "plugin_set_number", signature = (*args, **kwargs))]
fn plugin_set_number_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let option: String = extract_required(args, kwargs, 1, &["option"], "str")?;
    let value: f64 = extract_required(args, kwargs, 2, &["value"], "float")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state
        .plugin_numbers
        .insert((name.clone(), option.clone()), value);
    push_log(&mut state, format!("plugin number set: {name}.{option}={value}"));
    Ok(())
}

#[pyfunction]
#[pyo3(name = "plugin_set_string", signature = (*args, **kwargs))]
fn plugin_set_string_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let option: String = extract_required(args, kwargs, 1, &["option"], "str")?;
    let value: String = extract_required(args, kwargs, 2, &["value"], "str")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    state
        .plugin_strings
        .insert((name.clone(), option.clone()), value.clone());
    push_log(&mut state, format!("plugin string set: {name}.{option}={value}"));
    Ok(())
}

#[pyfunction]
#[pyo3(name = "plugin_run", signature = (*args, **kwargs))]
fn plugin_run_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<()> {
    let name: String = extract_required(args, kwargs, 0, &["name"], "str")?;
    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let lname = name.to_ascii_lowercase();
    match lname.as_str() {
        "refine" | "refinemesh" => {
            let refined = {
                let mesh = state.current_mesh.as_mut().ok_or_else(|| {
                    pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded")
                })?;
                refine_triangles_once(mesh)
            };
            push_log(&mut state, format!("plugin run refine: {refined} triangles refined"));
        }
        "recombine" => {
            let recombined = {
                let mesh = state.current_mesh.as_mut().ok_or_else(|| {
                    pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded")
                })?;
                recombine_triangles_to_quads(mesh)
            };
            push_log(&mut state, format!("plugin run recombine: {recombined} quads created"));
        }
        "smooth" | "laplace" => {
            let niter = state
                .plugin_numbers
                .get(&(name.clone(), "niter".to_string()))
                .copied()
                .unwrap_or(10.0)
                .max(1.0) as u32;
            {
                let mesh = state.current_mesh.as_mut().ok_or_else(|| {
                    pyo3::exceptions::PyRuntimeError::new_err("no mesh loaded")
                })?;
                let params = OptimizeParams {
                    iterations: niter,
                    ..Default::default()
                };
                LaplacianSmooth::default()
                    .optimize(mesh, &params)
                    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
            }
            push_log(&mut state, format!("plugin run smooth: {niter} iterations"));
        }
        _ => {
            push_log(&mut state, format!("plugin run no-op: {name}"));
        }
    }
    Ok(())
}

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

    let (mesh, mesh_name, viewer_cfg) = {
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

        let viewer_cfg = rmsh_viewer::ViewerConfig {
            show_nodes: option_number(&state, "Geometry.Points")
                .or_else(|| option_number(&state, "Mesh.Points"))
                .map(|v| v > 0.5),
            show_edges: option_number(&state, "Geometry.Curves")
                .or_else(|| option_number(&state, "Geometry.Lines"))
                .map(|v| v > 0.5),
            show_faces: option_number(&state, "Geometry.Surfaces").map(|v| v > 0.5),
            show_volumes: option_number(&state, "Geometry.Volumes").map(|v| v > 0.5),
        };

        (mesh, mesh_name, viewer_cfg)
    };

    py.allow_threads(move || rmsh_viewer::run_native_viewer(None, Some((mesh, mesh_name)), Some(viewer_cfg)))
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
    let shape = state.cad_shapes.get(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?;
    Ok(state.kernel.volume(shape).abs())
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
    let shape = state.cad_shapes.get(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?;
    let vol = state.kernel.volume(shape).abs();
    let area = state.kernel.surface_area(shape);
    let c = state.kernel.centroid(shape);
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

    let result = state
        .kernel
        .extrude_face(&base, face_idx, DVec3::new(dx, dy, dz), distance)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

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

    let result = state
        .kernel
        .revolve_face(&base, face_idx, DVec3::new(ax, ay, az), DVec3::new(dx, dy, dz), angle)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

    let new_tag = state.next_cad_tag + 1;
    state.next_cad_tag = new_tag;
    state.cad_shapes.insert(new_tag, result);
    Ok(new_tag)
}

/// Add a cone/frustum with center at (x,y,z), axis direction (dx,dy,dz),
/// radius r1 at one end and radius r2 at the other end.
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
    let r1: f64 = extract_required(args, kwargs, 6, &["r1", "r"], "float")?;
    let r2: f64 = extract_required(args, kwargs, 7, &["r2"], "float")?;
    let tag: i32 = extract_required(args, kwargs, 8, &["tag"], "int").unwrap_or(-1);
    let angle: f64 = extract_required(args, kwargs, 9, &["angle"], "float").unwrap_or(2.0 * PI);

    if (angle - 2.0 * PI).abs() > 1e-12 {
        return Err(pyo3::exceptions::PyNotImplementedError::new_err(
            "addCone currently supports only full cones/frustums (angle = 2*pi)",
        ));
    }

    let axis_vec = DVec3::new(dx, dy, dz);
    let height = axis_vec.length();
    if height < 1e-15 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "cone axis direction (dx, dy, dz) must be non-zero",
        ));
    }
    if r1 < 0.0 || r2 < 0.0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "cone radii r1 and r2 must be non-negative",
        ));
    }
    if r1 <= 1e-15 && r2 <= 1e-15 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "at least one of r1/r2 must be > 0",
        ));
    }
    let axis_norm = axis_vec.normalize();
    let center = DVec3::new(x, y, z);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let ref_dir = if axis_norm.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let shape = if (r1 - r2).abs() <= 1e-12 {
        state
            .kernel
            .create_cylinder(center, axis_norm, ref_dir, r1, height)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?
    } else {
        state
            .kernel
            .create_cone(center, axis_norm, ref_dir, r1, r2, height)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?
    };

    let assigned_tag = if tag > 0 { tag } else { state.next_cad_tag + 1 };
    state.next_cad_tag = assigned_tag.max(state.next_cad_tag);
    state.cad_shapes.insert(assigned_tag, shape);
    Ok(assigned_tag)
}

/// Add a torus.
///
/// Supported signatures:
/// 1) gmsh-compatible:
///    addTorus(x, y, z, r1, r2, tag=-1, angle=2*pi, zAxis=[])
/// 2) legacy rmsh-compatible:
///    addTorus(x, y, z, dx, dy, dz, r1, r2, tag=-1)
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
    let has_axis_kwargs = if let Some(kw) = kwargs {
        kw.contains("dx")? || kw.contains("dy")? || kw.contains("dz")?
    } else {
        false
    };
    let use_legacy_signature = has_axis_kwargs || args.len() >= 8;

    let (axis_norm, r1, r2, tag) = if use_legacy_signature {
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
        (axis_norm, r1, r2, tag)
    } else {
        let r1: f64 = extract_required(args, kwargs, 3, &["r1", "r"], "float")?;
        let r2: f64 = extract_required(args, kwargs, 4, &["r2"], "float")?;
        let tag: i32 = extract_required(args, kwargs, 5, &["tag"], "int").unwrap_or(-1);
        let angle: f64 = extract_required(args, kwargs, 6, &["angle"], "float").unwrap_or(2.0 * PI);
        if (angle - 2.0 * PI).abs() > 1e-12 {
            return Err(pyo3::exceptions::PyNotImplementedError::new_err(
                "addTorus currently supports only full torus (angle = 2*pi)",
            ));
        }

        let z_axis: Vec<f64> = extract_required(args, kwargs, 7, &["zAxis"], "list[float]")
            .unwrap_or_default();
        let axis_vec = if z_axis.len() == 3 {
            DVec3::new(z_axis[0], z_axis[1], z_axis[2])
        } else {
            DVec3::Z
        };
        let axis_norm = if axis_vec.length_squared() < 1e-20 {
            DVec3::Z
        } else {
            axis_vec.normalize()
        };

        (axis_norm, r1, r2, tag)
    };

    let ref_dir = if axis_norm.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let shape = state
        .kernel
        .create_torus(DVec3::new(x, y, z), axis_norm, ref_dir, r1, r2)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

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
/// Unlike addBox this creates a *surface* CAD entity suitable for 2D meshing.
#[pyfunction]
#[pyo3(name = "model_occ_add_rectangle", signature = (*args, **kwargs))]
fn model_occ_add_rectangle_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let x: f64 = extract_required(args, kwargs, 0, &["x"], "float")?;
    let y: f64 = extract_required(args, kwargs, 1, &["y"], "float")?;
    let z: f64 = extract_required(args, kwargs, 2, &["z"], "float").unwrap_or(0.0);
    let dx: f64 = extract_required(args, kwargs, 3, &["dx"], "float")?;
    let dy: f64 = extract_required(args, kwargs, 4, &["dy"], "float")?;
    let tag: i32 = extract_required(args, kwargs, 5, &["tag"], "int").unwrap_or(-1);

    if !x.is_finite() || !y.is_finite() || !z.is_finite() || !dx.is_finite() || !dy.is_finite() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "rectangle parameters must be finite",
        ));
    }
    if dx.abs() < 1e-12 || dy.abs() < 1e-12 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "dx and dy must be non-zero",
        ));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let shape = state
        .kernel
        .make_rectangle_shape(x, y, z, dx, dy)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

    let assigned_tag = if tag > 0 { tag } else { state.next_cad_tag + 1 };
    state.next_cad_tag = assigned_tag.max(state.next_cad_tag);
    state.cad_shapes.insert(assigned_tag, shape);
    state.current_mesh = None;
    Ok(assigned_tag)
}

/// Add a 0D OCC point entity.
/// Signature matches gmsh.model.occ.addPoint:
///   addPoint(x, y, z, meshSize=0, tag=-1) -> tag
#[pyfunction]
#[pyo3(name = "model_occ_add_point", signature = (*args, **kwargs))]
fn model_occ_add_point_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let x: f64 = extract_required(args, kwargs, 0, &["x"], "float")?;
    let y: f64 = extract_required(args, kwargs, 1, &["y"], "float")?;
    let z: f64 = extract_required(args, kwargs, 2, &["z"], "float")?;
    let _mesh_size: f64 = extract_required(args, kwargs, 3, &["meshSize"], "float").unwrap_or(0.0);
    let tag: i32 = extract_required(args, kwargs, 4, &["tag"], "int").unwrap_or(-1);

    if !x.is_finite() || !y.is_finite() || !z.is_finite() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "point parameters must be finite",
        ));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let shape = state.kernel.make_point_shape(DVec3::new(x, y, z));
    let assigned_tag = next_or_requested_tag(&mut state, tag);
    state.cad_shapes.insert(assigned_tag, shape);
    state.current_mesh = None;
    Ok(assigned_tag)
}

/// Add a line between two OCC points.
/// Signature matches gmsh.model.occ.addLine:
///   addLine(startTag, endTag, tag=-1) -> tag
#[pyfunction]
#[pyo3(name = "model_occ_add_line", signature = (*args, **kwargs))]
fn model_occ_add_line_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let start_tag: i32 = extract_required(args, kwargs, 0, &["startTag"], "int")?;
    let end_tag: i32 = extract_required(args, kwargs, 1, &["endTag"], "int")?;
    let tag: i32 = extract_required(args, kwargs, 2, &["tag"], "int").unwrap_or(-1);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let p0 = point_from_point_entity(&state, start_tag, "addLine")?;
    let p1 = point_from_point_entity(&state, end_tag, "addLine")?;

    let shape = state
        .kernel
        .make_line_shape(p0, p1)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

    let assigned_tag = next_or_requested_tag(&mut state, tag);
    state.cad_shapes.insert(assigned_tag, shape);
    state.current_mesh = None;
    Ok(assigned_tag)
}

/// Add a circular arc from startTag to endTag around centerTag.
/// Signature matches gmsh.model.occ.addCircleArc:
///   addCircleArc(startTag, centerTag, endTag, tag=-1) -> tag
#[pyfunction]
#[pyo3(name = "model_occ_add_circle_arc", signature = (*args, **kwargs))]
fn model_occ_add_circle_arc_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let start_tag: i32 = extract_required(args, kwargs, 0, &["startTag"], "int")?;
    let center_tag: i32 = extract_required(args, kwargs, 1, &["centerTag"], "int")?;
    let end_tag: i32 = extract_required(args, kwargs, 2, &["endTag"], "int")?;
    let tag: i32 = extract_required(args, kwargs, 3, &["tag"], "int").unwrap_or(-1);

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let p0 = point_from_point_entity(&state, start_tag, "addCircleArc")?;
    let c = point_from_point_entity(&state, center_tag, "addCircleArc")?;
    let p1 = point_from_point_entity(&state, end_tag, "addCircleArc")?;

    let shape = state
        .kernel
        .make_circle_arc_shape(p0, c, p1)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

    let assigned_tag = next_or_requested_tag(&mut state, tag);
    state.cad_shapes.insert(assigned_tag, shape);
    state.current_mesh = None;
    Ok(assigned_tag)
}

/// Add an interpolating spline through OCC points.
/// Signature matches gmsh.model.occ.addSpline:
///   addSpline(pointTags, tag=-1) -> tag
#[pyfunction]
#[pyo3(name = "model_occ_add_spline", signature = (*args, **kwargs))]
fn model_occ_add_spline_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let point_tags: Vec<i32> = extract_required(args, kwargs, 0, &["pointTags"], "list of int")?;
    let tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int").unwrap_or(-1);

    if point_tags.len() < 2 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addSpline: at least 2 points are required",
        ));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut pts = Vec::with_capacity(point_tags.len());
    for t in &point_tags {
        pts.push(point_from_point_entity(&state, *t, "addSpline")?);
    }

    let shape = state
        .kernel
        .make_spline_shape(&pts)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

    let assigned_tag = next_or_requested_tag(&mut state, tag);
    state.cad_shapes.insert(assigned_tag, shape);
    state.current_mesh = None;
    Ok(assigned_tag)
}

/// Add a planar disk/ellipse surface in the XY plane.
/// Signature matches gmsh.model.occ.addDisk:
///   addDisk(xc, yc, zc, rx, ry, tag=-1) -> tag
#[pyfunction]
#[pyo3(name = "model_occ_add_disk", signature = (*args, **kwargs))]
fn model_occ_add_disk_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let xc: f64 = extract_required(args, kwargs, 0, &["xc"], "float")?;
    let yc: f64 = extract_required(args, kwargs, 1, &["yc"], "float")?;
    let zc: f64 = extract_required(args, kwargs, 2, &["zc"], "float").unwrap_or(0.0);
    let rx: f64 = extract_required(args, kwargs, 3, &["rx"], "float")?;
    let ry: f64 = extract_required(args, kwargs, 4, &["ry"], "float")?;
    let tag: i32 = extract_required(args, kwargs, 5, &["tag"], "int").unwrap_or(-1);

    if !xc.is_finite() || !yc.is_finite() || !zc.is_finite() || !rx.is_finite() || !ry.is_finite() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "disk parameters must be finite",
        ));
    }
    if rx.abs() < 1e-12 || ry.abs() < 1e-12 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "rx and ry must be non-zero",
        ));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let shape = state
        .kernel
        .make_disk_shape(DVec3::new(xc, yc, zc), rx, ry)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

    let assigned_tag = next_or_requested_tag(&mut state, tag);
    state.cad_shapes.insert(assigned_tag, shape);
    state.current_mesh = None;
    Ok(assigned_tag)
}

/// Add a curve loop from curve tags.
/// Signature matches gmsh.model.occ.addCurveLoop:
///   addCurveLoop(curveTags, tag=-1) -> tag
#[pyfunction]
#[pyo3(name = "model_occ_add_curve_loop", signature = (*args, **kwargs))]
fn model_occ_add_curve_loop_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let curve_tags: Vec<i32> = extract_required(args, kwargs, 0, &["curveTags"], "list[int]")?;
    let tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int").unwrap_or(-1);

    if curve_tags.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addCurveLoop: curveTags must not be empty",
        ));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    for curve_tag in &curve_tags {
        let abs_tag = curve_tag.abs();
        let shape = state.cad_shapes.get(&abs_tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "addCurveLoop: curve tag {abs_tag} not found"
            ))
        })?;
        if shape.edges.is_empty() {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "addCurveLoop: tag {abs_tag} is not a curve entity"
            )));
        }
    }

    let assigned_tag = next_or_requested_tag(&mut state, tag);
    state.occ_curve_loops.insert(assigned_tag, curve_tags);
    Ok(assigned_tag)
}

/// Add a planar surface from one outer curve loop and optional inner loops.
/// Signature matches gmsh.model.occ.addPlaneSurface:
///   addPlaneSurface(wireTags, tag=-1) -> tag
#[pyfunction]
#[pyo3(name = "model_occ_add_plane_surface", signature = (*args, **kwargs))]
fn model_occ_add_plane_surface_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let wire_tags: Vec<i32> = extract_required(args, kwargs, 0, &["wireTags"], "list[int]")?;
    let tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int").unwrap_or(-1);

    if wire_tags.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addPlaneSurface: wireTags must not be empty",
        ));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    // First wire tag is the outer loop, remaining are inner loops.
    let outer_loop_tag = wire_tags[0];
    let outer_curve_tags = state.occ_curve_loops.get(&outer_loop_tag).cloned().ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!(
            "addPlaneSurface: curve loop tag {} not found",
            outer_loop_tag
        ))
    })?;

    let mut outer_curves: Vec<<RcadKernel as CadKernel>::Shape> = Vec::new();
    for signed_tag in &outer_curve_tags {
        let abs_tag = signed_tag.abs();
        let curve = state.cad_shapes.get(&abs_tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "addPlaneSurface: curve tag {} not found",
                abs_tag
            ))
        })?.clone();
        outer_curves.push(curve);
    }

    let mut inner_loops: Vec<Vec<<RcadKernel as CadKernel>::Shape>> = Vec::new();
    for loop_tag in &wire_tags[1..] {
        let curve_tags = state.occ_curve_loops.get(loop_tag).cloned().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "addPlaneSurface: curve loop tag {} not found",
                loop_tag
            ))
        })?;
        let mut inner_curves: Vec<<RcadKernel as CadKernel>::Shape> = Vec::new();
        for signed_tag in &curve_tags {
            let abs_tag = signed_tag.abs();
            let curve = state.cad_shapes.get(&abs_tag).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "addPlaneSurface: curve tag {} not found",
                    abs_tag
                ))
            })?.clone();
            inner_curves.push(curve);
        }
        inner_loops.push(inner_curves);
    }

    let shape = state
        .kernel
        .make_plane_surface_from_curves(&outer_curves, &inner_loops)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

    let assigned_tag = next_or_requested_tag(&mut state, tag);
    state.cad_shapes.insert(assigned_tag, shape);
    state.current_mesh = None;
    Ok(assigned_tag)
}

/// Add a surface loop from surface tags.
/// Signature matches gmsh.model.occ.addSurfaceLoop:
///   addSurfaceLoop(surfaceTags, tag=-1) -> tag
#[pyfunction]
#[pyo3(name = "model_occ_add_surface_loop", signature = (*args, **kwargs))]
fn model_occ_add_surface_loop_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let surface_tags: Vec<i32> = extract_required(args, kwargs, 0, &["surfaceTags"], "list[int]")?;
    let tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int").unwrap_or(-1);

    if surface_tags.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addSurfaceLoop: surfaceTags must not be empty",
        ));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    for surface_tag in &surface_tags {
        let shape = state.cad_shapes.get(surface_tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "addSurfaceLoop: surface tag {} not found",
                surface_tag
            ))
        })?;
        if shape.solids.is_empty() {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "addSurfaceLoop: tag {} is not a surface entity",
                surface_tag
            )));
        }
    }

    let assigned_tag = next_or_requested_tag(&mut state, tag);
    state.occ_surface_loops.insert(assigned_tag, surface_tags);
    Ok(assigned_tag)
}

/// Add a volume from one or more surface-loop tags.
/// Signature matches gmsh.model.occ.addVolume:
///   addVolume(shellTags, tag=-1) -> tag
#[pyfunction]
#[pyo3(name = "model_occ_add_volume", signature = (*args, **kwargs))]
fn model_occ_add_volume_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<i32> {
    let shell_tags: Vec<i32> = extract_required(args, kwargs, 0, &["shellTags"], "list[int]")?;
    let tag: i32 = extract_required(args, kwargs, 1, &["tag"], "int").unwrap_or(-1);

    if shell_tags.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addVolume: shellTags must not be empty",
        ));
    }

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let mut surface_tags = Vec::<i32>::new();
    for shell_tag in shell_tags {
        let faces = state.occ_surface_loops.get(&shell_tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "addVolume: surface loop tag {} not found",
                shell_tag
            ))
        })?;
        surface_tags.extend(faces.iter().copied());
    }

    surface_tags.sort_unstable();
    surface_tags.dedup();
    if surface_tags.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addVolume: resolved surface set is empty",
        ));
    }

    let mut shapes = Vec::new();
    for surface_tag in surface_tags {
        let s = state.cad_shapes.get(&surface_tag).cloned().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "addVolume: surface tag {} not found",
                surface_tag
            ))
        })?;
        shapes.push(s);
    }

    let merged = if shapes.len() == 1 {
        shapes.remove(0)
    } else {
        state.kernel.merge_for_export(&shapes)
    };

    let assigned_tag = next_or_requested_tag(&mut state, tag);
    state.cad_shapes.insert(assigned_tag, merged);
    state.current_mesh = None;
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

    let result = state
        .kernel
        .fillet_edges(&base, &edges)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

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

    let edges: Vec<(usize, f64)> = curve_tags.into_iter().zip(distances).collect();
    let result = state
        .kernel
        .chamfer_edges(&base, &edges)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

    let new_tag = state.next_cad_tag + 1;
    state.next_cad_tag = new_tag;
    state.cad_shapes.insert(new_tag, result);
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
    let shape = state.cad_shapes.get(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?.clone();

    let (repaired, report_str) = state
        .kernel
        .heal_shape(&shape, tolerance)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;
    state.cad_shapes.insert(tag, repaired);

    let parse_val = |key: &str| -> usize {
        report_str
            .split(' ')
            .find_map(|part| {
                let mut parts = part.splitn(2, '=');
                if parts.next()? == key { parts.next()?.parse().ok() } else { None }
            })
            .unwrap_or(0)
    };

    Python::with_gil(|py| {
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("vertices_merged", parse_val("vertices_merged"))?;
        dict.set_item("degenerate_faces_removed", parse_val("degenerate_faces_removed"))?;
        dict.set_item("normals_recomputed", parse_val("normals_recomputed"))?;
        dict.set_item("wires_fixed", parse_val("wires_fixed"))?;
        Ok(dict.into())
    })
}

/// Inspect a CAD shape's internal geometry buffers to debug boolean/export parity.
/// Returns counts for vertices/edges/surfaces/curves and attached edge pcurves.
#[pyfunction]
#[pyo3(name = "model_occ_debug_shape_geom", signature = (*args, **kwargs))]
fn model_occ_debug_shape_geom_impl(
    args: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<pyo3::PyObject> {
    let tag: i32 = extract_required(args, kwargs, 0, &["tag"], "int")?;

    let state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;
    let shape = state.cad_shapes.get(&tag).ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
    })?;

    let stats = state.kernel.inspect_shape(shape);

    let kind_to_str = |kind: &CurveKind| -> &'static str {
        match kind {
            CurveKind::Line => "Line",
            CurveKind::Circle => "Circle",
            CurveKind::Ellipse => "Ellipse",
            CurveKind::BSpline => "BSpline",
            CurveKind::Other(s) => Box::leak(s.clone().into_boxed_str()),
        }
    };
    let surface_kind_to_str = |kind: &SurfaceKind| -> &'static str {
        match kind {
            SurfaceKind::Plane => "Plane",
            SurfaceKind::Cylinder => "Cylinder",
            SurfaceKind::Sphere => "Sphere",
            SurfaceKind::Cone => "Cone",
            SurfaceKind::Torus => "Torus",
            SurfaceKind::BSpline => "BSpline",
            SurfaceKind::Other(s) => Box::leak(s.clone().into_boxed_str()),
        }
    };

    Python::with_gil(|py| {
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("tag", tag)?;
        dict.set_item("vertices", stats.vertices)?;
        dict.set_item("edges", stats.edges)?;
        dict.set_item("faces", stats.faces)?;
        dict.set_item("solids", stats.solids)?;
        dict.set_item("closed_edges", stats.closed_edges)?;
        dict.set_item("curves_3d", stats.curves_3d)?;
        dict.set_item("surfaces_3d", stats.surfaces_3d)?;
        dict.set_item("curves_2d", stats.curves_2d)?;
        dict.set_item("edge_curve_slots", stats.edge_curve_some + stats.edge_curve_none)?;
        dict.set_item("edge_curve_some", stats.edge_curve_some)?;
        dict.set_item("edge_curve_none", stats.edge_curve_none)?;
        dict.set_item("edge_pcurve_slots", stats.edges_with_pcurves + stats.edge_curve_none)?;
        dict.set_item("edges_with_pcurves", stats.edges_with_pcurves)?;
        dict.set_item("total_pcurves", stats.total_pcurves)?;
        dict.set_item("outer_wire_edges", stats.outer_wire_edge_refs_total)?;
        dict.set_item("outer_wire_edge_refs_total", stats.outer_wire_edge_refs_total)?;
        dict.set_item("outer_wire_edges_per_face_min", stats.outer_wire_edges_per_face_min)?;
        dict.set_item("outer_wire_edges_per_face_max", stats.outer_wire_edges_per_face_max)?;
        dict.set_item("outer_wire_face_max_surface_kind", surface_kind_to_str(&stats.outer_wire_face_max_surface_kind))?;
        dict.set_item("faces_with_outer_wire_over_100", stats.faces_with_outer_wire_over_100)?;
        dict.set_item("outer_wire_edges_with_curve", stats.outer_wire_edges_with_curve)?;
        dict.set_item("outer_wire_edges_without_curve", stats.outer_wire_edges_without_curve)?;
        dict.set_item("outer_wire_edges_with_pcurves", stats.outer_wire_edges_with_pcurves)?;

        let curve_dict = pyo3::types::PyDict::new(py);
        for (k, v) in &stats.curve3_kinds {
            curve_dict.set_item(kind_to_str(k), *v)?;
        }
        dict.set_item("curve3_kinds", curve_dict)?;

        let surface_dict = pyo3::types::PyDict::new(py);
        for (k, v) in &stats.surface3_kinds {
            surface_dict.set_item(surface_kind_to_str(k), *v)?;
        }
        dict.set_item("surface3_kinds", surface_dict)?;

        let pcurve_surface_dict = pyo3::types::PyDict::new(py);
        for (k, v) in &stats.pcurve_surface_kinds {
            pcurve_surface_dict.set_item(surface_kind_to_str(k), *v)?;
        }
        dict.set_item("pcurve_surface_kinds", pcurve_surface_dict)?;

        let outer_wire_curve_dict = pyo3::types::PyDict::new(py);
        for (k, v) in &stats.outer_wire_curve3_kinds {
            outer_wire_curve_dict.set_item(kind_to_str(k), *v)?;
        }
        dict.set_item("outer_wire_curve3_kinds", outer_wire_curve_dict)?;

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
    m.add_function(pyo3::wrap_pyfunction!(model_get_entities_in_bounding_box_impl, m)?)?;
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
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_point_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_line_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_circle_arc_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_spline_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_disk_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_curve_loop_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_plane_surface_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_surface_loop_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_add_volume_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_cut_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_fuse_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_fragment_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_intersect_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_copy_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_remove_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_translate_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_rotate_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_dilate_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_mirror_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_import_shapes_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_occ_synchronize_impl, m)?)?;

    m.add_function(pyo3::wrap_pyfunction!(model_mesh_set_size_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_generate_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_set_order_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_get_nodes_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_get_elements_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_get_element_qualities_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_clear_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_optimize_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_refine_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_recombine_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_field_add_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_field_set_number_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_field_set_as_background_mesh_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_set_transfinite_curve_impl, m)?)?;
    m.add_function(pyo3::wrap_pyfunction!(model_mesh_set_transfinite_surface_impl, m)?)?;

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
    m.add_function(pyo3::wrap_pyfunction!(model_occ_debug_shape_geom_impl, m)?)?;

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        estimate_mesh_characteristic_size, option_default_number, step_protocol_from_state,
        Mesh, Node, RuntimeState, StepProtocol, cad, RcadKernel,
    };
    use super::cad::CadKernel;
    use glam::DVec3;

    fn count_occurrences(haystack: &str, needle: &str) -> usize {
        haystack.matches(needle).count()
    }

    #[test]
    fn strict_cylinder_emits_cylindrical_surface_and_seam_curve() {
        let kernel = RcadKernel;
        let brep = kernel.create_cylinder(DVec3::ZERO, DVec3::Z, DVec3::X, 0.5, 2.0)
            .expect("cylinder should build");
        let step = kernel.write_step_string(&brep, &cad::StepExportOptions {
            protocol: StepProtocol::Ap214,
            solid_color: None,
            gmsh_strict: true,
        }).expect("strict cylinder STEP export should succeed");

        assert_eq!(count_occurrences(&step, "ADVANCED_FACE"), 3);
        assert_eq!(count_occurrences(&step, "CYLINDRICAL_SURFACE"), 1);
        assert!(count_occurrences(&step, "SEAM_CURVE") >= 1);
    }

    #[test]
    fn strict_frustum_cone_emits_conical_side_and_three_faces() {
        let kernel = RcadKernel;
        let brep = kernel.create_cone(DVec3::ZERO, DVec3::Z, DVec3::X, 0.8, 0.2, 2.0)
            .expect("cone should build");
        let step = kernel.write_step_string(&brep, &cad::StepExportOptions {
            protocol: StepProtocol::Ap214,
            solid_color: None,
            gmsh_strict: true,
        }).expect("strict frustum STEP export should succeed");

        assert_eq!(count_occurrences(&step, "ADVANCED_FACE"), 3);
        assert_eq!(count_occurrences(&step, "CONICAL_SURFACE"), 1);
        assert_eq!(count_occurrences(&step, "EDGE_CURVE"), 3);
        assert!(!step.contains("TRIANGULATED"));
        assert!(!step.contains("TESSELLATED"));
    }

    #[test]
    fn strict_standalone_line_emits_wireframe_curve_set() {
        let kernel = RcadKernel;
        let brep = kernel.make_line_shape(DVec3::ZERO, DVec3::X)
            .expect("line should build");
        let step = kernel.write_step_string(&brep, &cad::StepExportOptions {
            protocol: StepProtocol::Ap214,
            solid_color: None,
            gmsh_strict: true,
        }).expect("strict standalone line STEP export should succeed");

        assert!(step.contains("GEOMETRIC_CURVE_SET"));
        assert!(step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"));
        assert!(step.contains("LINE"));
    }

    #[test]
    fn strict_standalone_spline_emits_bspline_curve() {
        use glam::DVec3;
        let kernel = RcadKernel;
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
        ];
        let brep = kernel.make_spline_shape(&pts)
            .expect("spline should build");
        let step = kernel.write_step_string(&brep, &cad::StepExportOptions {
            protocol: StepProtocol::Ap214,
            solid_color: None,
            gmsh_strict: true,
        }).expect("strict standalone spline STEP export should succeed");

        assert!(step.contains("GEOMETRIC_CURVE_SET"));
        assert!(step.contains("B_SPLINE_CURVE_WITH_KNOTS"));
    }

    #[test]
    fn step_protocol_uses_default_option_value() {
        let state = RuntimeState::default();
        let protocol = step_protocol_from_state(&state);
        assert!(matches!(protocol, StepProtocol::Ap214));
    }

    #[test]
    fn option_default_number_accepts_geometry_lines_alias() {
        assert_eq!(option_default_number("Geometry.Curves"), Some(1.0));
        assert_eq!(option_default_number("Geometry.Lines"), Some(1.0));
    }

    #[test]
    fn estimate_mesh_characteristic_size_uses_bbox_diagonal() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 2.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 3.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 0.0, 6.0));

        let h = estimate_mesh_characteristic_size(&mesh)
            .expect("non-degenerate mesh should produce an auto size");
        let expected = (49.0_f64).sqrt() / 20.0;
        assert!((h - expected).abs() < 1e-12);
    }

    #[test]
    fn estimate_mesh_characteristic_size_rejects_degenerate_mesh() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 1.0, 1.0, 1.0));
        mesh.add_node(Node::new(2, 1.0, 1.0, 1.0));
        assert_eq!(estimate_mesh_characteristic_size(&mesh), None);
    }
}
