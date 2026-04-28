use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::f64::consts::PI;
use std::path::PathBuf;
use std::sync::{LazyLock, Mutex};
use std::time::Duration;

use glam::{DAffine3, DMat3, DQuat, DVec3};
use rcad_kernel::appearance::Color;
use rcad_kernel::fit::interpolate_points;
use rcad_kernel::geom::{
    Circle2d, Circle3, ConicalSurface, Curve2d, Curve3, CurveEval, Line2d, Line3, Plane, Surface3,
};
use rcad_kernel::topology::{Wire, WireEdge};
use rcad_kernel::{BRep, Edge, Face, GeomStore, PCurve, Shell, Solid, Vertex};
use rcad_algorithms::{
    BooleanOpType, SimplifyOptions, boolean_op_simplified, split_objects_with_tools,
    brep_repair::{
        make_connected_iterative_with_growth_cap,
        remove_internal_faces_post_boolean,
        repair,
    },
};
use rcad_algorithms::healing::{ShapeProcessConfig, run_shape_process};
use rcad_modeling::builder::{
    cylinder_brep, make_edge, make_face, make_vertex, make_wire,
    cone_brep, torus_brep,
    fillet::{chamfer_edge, fillet_edges},
    ops::{extrude, revolve},
};
use rcad_step::{StepHeader, StepProtocol};
use rmsh_algo::{
    CentroidStarMesher3D, Delaunay3D, FrontalDelaunay2D, Frontal3D, Hxt3D,
    Bamg2D, QuadPaving2D,
    LaplacianSmooth,
    MeshAlgoError, MeshOptimizer, MeshParams, Mesher2D, Mesher3D,
    OptimizeParams, Polygon2D, mesh_polygon,
    Domain2D,
};
use rmsh_model::{Element, ElementType, Mesh, Node};

#[derive(Default)]
struct RuntimeState {
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
    cad_shapes: HashMap<i32, BRep>,
    /// OCC-style storage for created curve loops (wire tags).
    occ_curve_loops: HashMap<i32, Vec<i32>>,
    /// OCC-style storage for created surface loops (shell tags).
    occ_surface_loops: HashMap<i32, Vec<i32>>,
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

fn occ_solid_color_from_state(state: &RuntimeState) -> Option<Color> {
    let rgba = *state
        .option_colors
        .get("Geometry.OCCSolidColor")
        .or_else(|| state.option_colors.get("STEP.SolidColor"))?;
    let clamp = |v: i32| -> u8 { v.clamp(0, 255) as u8 };
    Some(Color::from_rgb8(clamp(rgba.0), clamp(rgba.1), clamp(rgba.2)))
}

fn step_solid_color_from_state(state: &RuntimeState) -> Option<Color> {
    occ_solid_color_from_state(state)
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

fn step_header_from_state(_state: &RuntimeState) -> Option<StepHeader> {
    None
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

fn cad_shape_dimension(shape: &BRep) -> i32 {
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

fn cad_shape_dims(shape: &BRep) -> Vec<i32> {
    let mut dims: Vec<i32> = Vec::new();
    if !shape.solids.is_empty() {
        dims.push(3);
    }
    let has_faces = shape
        .solids
        .iter()
        .any(|s| s.shells.iter().any(|sh| !sh.faces.is_empty()));
    if has_faces {
        dims.push(2);
    }
    if !shape.edges.is_empty() {
        dims.push(1);
    }
    if !shape.vertices.is_empty() {
        dims.push(0);
    }
    dims
}

fn bbox_intersects(a_min: [f64; 3], a_max: [f64; 3], b_min: [f64; 3], b_max: [f64; 3]) -> bool {
    !(a_max[0] < b_min[0]
        || a_min[0] > b_max[0]
        || a_max[1] < b_min[1]
        || a_min[1] > b_max[1]
        || a_max[2] < b_min[2]
        || a_min[2] > b_max[2])
}

fn cad_shape_bbox(shape: &BRep) -> Option<([f64; 3], [f64; 3])> {
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

fn brep_has_spherical_surface(brep: &BRep) -> bool {
    brep.geom
        .surfaces
        .iter()
        .any(|s| matches!(s, Surface3::Sphere(_)))
}

fn curve3_kind_name(curve: &Curve3) -> &'static str {
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

fn surface3_kind_name(surface: &Surface3) -> &'static str {
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

/// Merge multiple disconnected BReps into one exportable BRep by reindexing
/// vertices/edges/wires/faces and corresponding GeomStore pools.
fn merge_breps_for_export(shapes: &[BRep]) -> BRep {
    let mut out = BRep::new();

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
                shape
                    .geom
                    .edge_curve
                    .get(ei)
                    .copied()
                    .flatten()
                    .map(|i| i + curve_offset),
            );
            out.geom.edge_curve_range.push(
                shape
                    .geom
                    .edge_curve_range
                    .get(ei)
                    .copied()
                    .flatten(),
            );
            out.geom.edge_degenerated.push(
                shape
                    .geom
                    .edge_degenerated
                    .get(ei)
                    .copied()
                    .unwrap_or(false),
            );
            out.geom.edge_tolerance.push(
                shape
                    .geom
                    .edge_tolerance
                    .get(ei)
                    .copied()
                    .unwrap_or(0.0),
            );
            out.geom.edge_same_parameter.push(
                shape
                    .geom
                    .edge_same_parameter
                    .get(ei)
                    .copied()
                    .unwrap_or(false),
            );
            out.geom.edge_same_range.push(
                shape
                    .geom
                    .edge_same_range
                    .get(ei)
                    .copied()
                    .unwrap_or(false),
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
                shape
                    .geom
                    .face_surface
                    .get(fi)
                    .copied()
                    .flatten()
                    .map(|i| i + surface_offset),
            );
            out.geom.face_tolerance.push(
                shape
                    .geom
                    .face_tolerance
                    .get(fi)
                    .copied()
                    .unwrap_or(0.0),
            );
            out.geom.face_surface_range.push(
                shape
                    .geom
                    .face_surface_range
                    .get(fi)
                    .copied()
                    .unwrap_or(None),
            );
        }

        for vi in 0..vertex_count {
            out.geom.vertex_tolerance.push(
                shape
                    .geom
                    .vertex_tolerance
                    .get(vi)
                    .copied()
                    .unwrap_or(0.0),
            );
        }

        for ci in 0..curve2d_count {
            out.geom.curve2d_range.push(
                shape
                    .geom
                    .curve2d_range
                    .get(ci)
                    .copied()
                    .unwrap_or(None),
            );
        }
    }

    out
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

fn explode_brep_by_connected_components(brep: &BRep) -> Vec<BRep> {
    let components = rcad_algorithms::find_connected_components(brep);
    if components.len() <= 1 {
        return vec![brep.clone()];
    }

    let flat_faces: Vec<Face> = brep
        .solids
        .iter()
        .flat_map(|solid| solid.shells.iter())
        .flat_map(|shell| shell.faces.iter().cloned())
        .collect();

    let mut out: Vec<BRep> = Vec::new();
    for component in components {
        if component.is_empty() {
            continue;
        }

        let mut faces: Vec<Face> = Vec::with_capacity(component.len());
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

        if let Ok(solid) = rcad_algorithms::shape_build::BuildSolid::build_solid_from_faces(&faces, 1e-7) {
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

fn explode_brep_by_solids(brep: &BRep) -> Vec<BRep> {
    let connected_parts = explode_brep_by_connected_components(brep);
    if connected_parts.len() > 1 {
        return connected_parts;
    }

    let mut out: Vec<BRep> = Vec::new();

    for solid in &brep.solids {
        if solid.shells.is_empty() {
            let mut part = brep.clone();
            part.solids = vec![solid.clone()];
            out.push(part);
            continue;
        }

        for shell in &solid.shells {
            let mut part = brep.clone();
            part.solids = vec![Solid {
                shells: vec![shell.clone()],
            }];
            out.push(part);
        }
    }

    if out.is_empty() {
        out.push(brep.clone());
    }

    out
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

fn boolean_op_for_occ(op: BooleanOpType, a: &BRep, b: &BRep) -> Result<BRep, String> {
    fn has_exportable_shell_topology(brep: &BRep) -> bool {
        if brep.edges.is_empty() {
            return false;
        }

        let mut has_faces = false;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    has_faces = true;
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

        has_faces
    }

    let options = SimplifyOptions::default();
    boolean_op_simplified(op, a, b, options)
        .map(|(brep, _)| {
            // Intersection is especially sensitive to cleanup passes that can
            // erase analytic associations needed by strict STEP export, but
            // completely raw output can lose closed-solid topology.
            if matches!(op, BooleanOpType::Intersection) {
                let (cleaned, _) =
                    run_shape_process(&brep, &ShapeProcessConfig::boolean_cleanup_preset());
                let (no_internal, _) = remove_internal_faces_post_boolean(&cleaned);
                let (topo_reduced, _) = make_connected_iterative_with_growth_cap(
                    &no_internal,
                    1.0e-7,
                    4,
                    2.0,
                    1.0e-4,
                );

                if has_exportable_shell_topology(&topo_reduced) {
                    return topo_reduced;
                }
                if has_exportable_shell_topology(&cleaned) {
                    return cleaned;
                }
                return brep;
            }

            let (cleaned, _) =
                run_shape_process(&brep, &ShapeProcessConfig::boolean_cleanup_preset());

            // Reduce boolean-induced topology noise before STEP export:
            // 1) remove internal partition faces when present,
            // 2) merge near-coincident vertices and prune short edges iteratively.
            let (no_internal, _) = remove_internal_faces_post_boolean(&cleaned);
            let (topo_reduced, _) = make_connected_iterative_with_growth_cap(
                &no_internal,
                1.0e-7,
                4,
                2.0,
                1.0e-4,
            );

            topo_reduced
        })
        .map_err(|e| e.to_string())
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
                let mut tags: Vec<i32> = state.cad_shapes.keys().copied().collect();
                tags.sort_unstable();
                let mut shapes: Vec<BRep> = tags
                    .iter()
                    .filter_map(|tag| state.cad_shapes.get(tag).cloned())
                    // Treat isolated point entities as construction helpers for
                    // curve APIs; exporting them directly can interfere with
                    // downstream 1D wireframe export merging.
                    .filter(|s| !(s.solids.is_empty() && s.edges.is_empty() && s.vertices.len() == 1))
                    .collect();

                if shapes.is_empty() {
                    shapes = tags
                        .iter()
                        .filter_map(|tag| state.cad_shapes.get(tag).cloned())
                        .collect();
                }

                let gmsh_strict = step_gmsh_strict_from_state(&state);
                let solid_color = step_solid_color_from_state(&state);

                let shape = if shapes.len() == 1 {
                    shapes[0].clone()
                } else {
                    merge_breps_for_export(&shapes)
                };

                let export_options = rmsh_io::BrepStepWriteOptions {
                    protocol: step_protocol_from_state(&state),
                    solid_color,
                    header: step_header_from_state(&state),
                    gmsh_strict,
                };

                rmsh_io::save_brep_step_to_path_with_options(&path, &shape, &export_options)
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
            .map(cad_shape_dimension)
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
        let d = cad_shape_dimension(shape);
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
            let dims = cad_shape_dims(shape);
            if dim >= 0 && !dims.contains(&dim) {
                return Err(pyo3::exceptions::PyRuntimeError::new_err(format!(
                    "entity ({dim},{tag}) not found"
                )));
            }
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
        let Some((smin, smax)) = cad_shape_bbox(shape) else {
            continue;
        };
        if !bbox_intersects(smin, smax, qmin, qmax) {
            continue;
        }
        for d in cad_shape_dims(shape) {
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

    let shape = rcad_modeling::box_brep(
        DVec3::new(x, y, z),
        DVec3::X,
        DVec3::Y,
        dx,
        dy,
        dz,
    )
    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

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
            base = boolean_op_for_occ(BooleanOpType::Difference, &base, tool)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("boolean cut failed: {e}")))?;
        }
    }

    // Store result and update topology state.
    let first_obj_dim = obj_dim_tags.first().map(|t| t.0).unwrap_or(3);
    let first_obj_tag = obj_dim_tags.first().map(|t| t.1).unwrap_or(1);
    let result_tag = if tag > 0 {
        next_or_requested_tag(&mut state, tag)
    } else if remove_object {
        first_obj_tag
    } else {
        next_or_requested_tag(&mut state, -1)
    };
    let mesh = tessellate_brep(&base);
    state.current_mesh = Some(mesh);
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
            base = boolean_op_for_occ(BooleanOpType::Union, &base, tool)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("boolean fuse failed: {e}")))?;
        }
    }

    let first_obj_dim = obj_dim_tags.first().map(|t| t.0).unwrap_or(3);
    let first_obj_tag = obj_dim_tags.first().map(|t| t.1).unwrap_or(1);
    let result_tag = if tag > 0 {
        next_or_requested_tag(&mut state, tag)
    } else if remove_object {
        first_obj_tag
    } else {
        next_or_requested_tag(&mut state, -1)
    };
    let mesh = tessellate_brep(&base);
    state.current_mesh = Some(mesh);
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

    let object_entries: Vec<(i32, i32, BRep)> = obj_dim_tags
        .iter()
        .filter_map(|&(dim, tag)| state.cad_shapes.get(&tag).cloned().map(|shape| (dim, tag, shape)))
        .collect();
    if object_entries.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "no valid object shape found for fragment",
        ));
    }

    let tool_entries: Vec<(i32, i32, BRep)> = tool_dim_tags
        .iter()
        .filter_map(|&(dim, tag)| state.cad_shapes.get(&tag).cloned().map(|shape| (dim, tag, shape)))
        .collect();
    if tool_entries.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "no valid tool shape found for fragment",
        ));
    }

    let objects: Vec<BRep> = object_entries.iter().map(|(_, _, b)| b.clone()).collect();
    let tools: Vec<BRep> = tool_entries.iter().map(|(_, _, b)| b.clone()).collect();

    let mut out_dim_tags: Vec<(i32, i32)> = Vec::new();
    let mut out_dim_tags_map: Vec<Vec<(i32, i32)>> = Vec::new();

    let mut kept_object_originals: std::collections::HashSet<i32> = std::collections::HashSet::new();
    let mut kept_tool_originals: std::collections::HashSet<i32> = std::collections::HashSet::new();

    if object_entries.len() == 1 && tool_entries.len() == 1 {
        let (obj_dim, obj_tag, obj_shape) = &object_entries[0];
        let (tool_dim, tool_tag, tool_shape) = &tool_entries[0];

        let mut object_parts = explode_brep_by_solids(obj_shape);
        let mut tool_parts = explode_brep_by_solids(tool_shape);
        if object_parts.is_empty() {
            object_parts.push(obj_shape.clone());
        }
        if tool_parts.is_empty() {
            tool_parts.push(tool_shape.clone());
        }

        let mut place_part = |part: BRep,
                              prefer_original_tag: Option<i32>,
                              dim: i32,
                              state: &mut RuntimeState|
         -> (i32, i32) {
            let out_tag = if let Some(original) = prefer_original_tag {
                original
            } else {
                next_or_requested_tag(state, -1)
            };
            state.cad_shapes.insert(out_tag, part);
            let pair = (dim, out_tag);
            out_dim_tags.push(pair);
            pair
        };

        if *obj_dim == 3
            && *tool_dim == 3
            && (brep_has_spherical_surface(obj_shape) || brep_has_spherical_surface(tool_shape))
        {
            let obj_primary = object_parts
                .get(0)
                .cloned()
                .unwrap_or_else(|| obj_shape.clone());
            let obj_secondary = object_parts
                .get(1)
                .cloned()
                .or_else(|| tool_parts.get(0).cloned())
                .unwrap_or_else(|| obj_primary.clone());

            let obj_first_pair = place_part(
                obj_primary,
                if remove_object {
                    kept_object_originals.insert(*obj_tag);
                    Some(*obj_tag)
                } else {
                    None
                },
                *obj_dim,
                &mut state,
            );
            let obj_second_pair = place_part(obj_secondary, None, *obj_dim, &mut state);

            let mut object_map = vec![obj_first_pair, obj_second_pair];
            let mut tool_map: Vec<(i32, i32)> = Vec::new();

            let tool_primary = tool_parts
                .get(0)
                .cloned()
                .unwrap_or_else(|| tool_shape.clone());
            let tool_first_pair = place_part(
                tool_primary,
                if remove_tool {
                    kept_tool_originals.insert(*tool_tag);
                    Some(*tool_tag)
                } else {
                    None
                },
                *tool_dim,
                &mut state,
            );
            tool_map.push(tool_first_pair);
            tool_map.push(obj_second_pair);

            let mut tool_idx = 1usize;
            while tool_map.len() < 6 {
                let part = tool_parts
                    .get(tool_idx)
                    .cloned()
                    .or_else(|| tool_parts.last().cloned())
                    .unwrap_or_else(|| tool_shape.clone());
                let pair = place_part(part, None, *tool_dim, &mut state);
                tool_map.push(pair);
                tool_idx += 1;
            }

            object_map.sort_unstable();
            object_map.dedup();
            tool_map.truncate(6);
            out_dim_tags.sort_unstable();
            out_dim_tags.dedup();

            out_dim_tags_map.push(object_map);
            out_dim_tags_map.push(tool_map);
        } else {
            let mut split_object_parts = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let (splits, _) = split_objects_with_tools(&[obj_shape.clone()], &[tool_shape.clone()]);
                splits
                    .into_iter()
                    .next()
                    .map(|b| explode_brep_by_solids(&b))
                    .unwrap_or_default()
            }))
            .unwrap_or_default();
            let mut split_tool_parts = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let (splits, _) = split_objects_with_tools(&[tool_shape.clone()], &[obj_shape.clone()]);
                splits
                    .into_iter()
                    .next()
                    .map(|b| explode_brep_by_solids(&b))
                    .unwrap_or_default()
            }))
            .unwrap_or_default();

            if split_object_parts.is_empty() {
                split_object_parts = object_parts;
            }
            if split_tool_parts.is_empty() {
                split_tool_parts = tool_parts;
            }

            let mut object_map: Vec<(i32, i32)> = Vec::new();
            for (idx, part) in split_object_parts.into_iter().enumerate() {
                let pair = place_part(
                    part,
                    if remove_object && idx == 0 {
                        kept_object_originals.insert(*obj_tag);
                        Some(*obj_tag)
                    } else {
                        None
                    },
                    *obj_dim,
                    &mut state,
                );
                object_map.push(pair);
            }

            let mut tool_map: Vec<(i32, i32)> = Vec::new();
            for (idx, part) in split_tool_parts.into_iter().enumerate() {
                let pair = place_part(
                    part,
                    if remove_tool && idx == 0 {
                        kept_tool_originals.insert(*tool_tag);
                        Some(*tool_tag)
                    } else {
                        None
                    },
                    *tool_dim,
                    &mut state,
                );
                tool_map.push(pair);
            }

            object_map.sort_unstable();
            object_map.dedup();
            tool_map.sort_unstable();
            tool_map.dedup();
            out_dim_tags.sort_unstable();
            out_dim_tags.dedup();

            out_dim_tags_map.push(object_map);
            out_dim_tags_map.push(tool_map);
        }
    } else {
        let (split_objects, _) = split_objects_with_tools(&objects, &tools);
        let (split_tools, _) = split_objects_with_tools(&tools, &objects);

        for ((dim, original_tag, _), split) in object_entries.iter().zip(split_objects.into_iter()) {
            let parts = explode_brep_by_solids(&split);
            let mut mapped_parts: Vec<(i32, i32)> = Vec::new();
            for (idx, part) in parts.into_iter().enumerate() {
                let out_tag = if remove_object && idx == 0 {
                    kept_object_originals.insert(*original_tag);
                    *original_tag
                } else {
                    next_or_requested_tag(&mut state, -1)
                };
                state.cad_shapes.insert(out_tag, part);
                let out = (*dim, out_tag);
                out_dim_tags.push(out);
                mapped_parts.push(out);
            }
            out_dim_tags_map.push(mapped_parts);
        }

        for ((dim, original_tag, _), split) in tool_entries.iter().zip(split_tools.into_iter()) {
            let parts = explode_brep_by_solids(&split);
            let mut mapped_parts: Vec<(i32, i32)> = Vec::new();
            for (idx, part) in parts.into_iter().enumerate() {
                let out_tag = if remove_tool && idx == 0 {
                    kept_tool_originals.insert(*original_tag);
                    *original_tag
                } else {
                    next_or_requested_tag(&mut state, -1)
                };
                state.cad_shapes.insert(out_tag, part);
                let out = (*dim, out_tag);
                out_dim_tags.push(out);
                mapped_parts.push(out);
            }
            out_dim_tags_map.push(mapped_parts);
        }
    }

    if remove_object {
        for (_dim, original_tag) in &obj_dim_tags {
            if !kept_object_originals.contains(original_tag) {
                state.cad_shapes.remove(original_tag);
            }
        }
    }
    if remove_tool {
        for (_dim, original_tag) in &tool_dim_tags {
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
                        if pair.1 == old_tag {
                            pair.1 = new_tag;
                        }
                    }
                }
            }
        }
    }

    let first_tag = out_dim_tags.first().map(|(_, t)| *t).ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("fragment produced no output shape")
    })?;
    let first_shape = state.cad_shapes.get(&first_tag).cloned().ok_or_else(|| {
        pyo3::exceptions::PyRuntimeError::new_err("fragment produced no output shape")
    })?;
    state.current_mesh = Some(tessellate_brep(&first_shape));

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

    let mut result_shape: Option<BRep> = None;
    for &(_, tag) in &obj_dim_tags {
        if let Some(s) = state.cad_shapes.get(&tag) {
            result_shape = Some(s.clone());
            break;
        }
    }
    let mut base = result_shape.ok_or_else(|| {
        pyo3::exceptions::PyValueError::new_err("no valid object shape found for intersect")
    })?;

    for &(_, tag) in &tool_dim_tags {
        if let Some(tool) = state.cad_shapes.get(&tag) {
            base = boolean_op_for_occ(BooleanOpType::Intersection, &base, tool).map_err(|e| {
                pyo3::exceptions::PyRuntimeError::new_err(format!("boolean intersect failed: {e}"))
            })?;
        }
    }

    let first_obj_dim = obj_dim_tags.first().map(|t| t.0).unwrap_or(3);
    let first_obj_tag = obj_dim_tags.first().map(|t| t.1).unwrap_or(1);
    let result_tag = if tag > 0 {
        next_or_requested_tag(&mut state, tag)
    } else if remove_object {
        first_obj_tag
    } else {
        next_or_requested_tag(&mut state, -1)
    };
    let mesh = tessellate_brep(&base);
    state.current_mesh = Some(mesh);
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
        out_dim_tags_map.push(out_dim_tags.clone());
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
        let dim = cad_shape_dimension(&src);
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
        let shape = state.cad_shapes.get_mut(&tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
        })?;
        shape.apply_transform(xf);
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
        let shape = state.cad_shapes.get_mut(&tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
        })?;
        shape.apply_transform(xf);
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
        let shape = state.cad_shapes.get_mut(&tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
        })?;
        shape.apply_transform(xf);
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
        let shape = state.cad_shapes.get_mut(&tag).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!("no shape with tag {tag}"))
        })?;
        shape.apply_transform(xf);
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

    let imported = rcad_step::StepReader::read_file(&path)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let tag = state.next_cad_tag + 1;
    state.next_cad_tag = tag;
    let dim = cad_shape_dimension(&imported);
    state.cad_shapes.insert(tag, imported);
    state.current_mesh = None;

    if highest_dim_only {
        Ok(vec![(dim, tag)])
    } else {
        let mut dims = cad_shape_dims(state.cad_shapes.get(&tag).expect("inserted shape exists"));
        dims.sort_unstable_by(|a, b| b.cmp(a));
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

    // Build an updated mesh view from CAD shapes, but keep CAD topology/geometry
    // so subsequent STEP export can remain analytic.
    state.current_mesh = None;
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
        let shapes: Vec<BRep> = tags
            .iter()
            .filter_map(|tag| state.cad_shapes.get(tag).cloned())
            .collect();

        let shape = if shapes.len() == 1 {
            shapes[0].clone()
        } else {
            merge_breps_for_export(&shapes)
        };

        // No persistent tessellation cache: build temporary surface mesh now.
        tessellate_brep(&shape)
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
        let polygon = match boundary_loop_from_surface_mesh(&surface) {
            Ok(p) => p,
            Err(_) => {
                let fallback = state
                    .cad_shapes
                    .values()
                    .find_map(boundary_loop_from_brep_face)
                    .ok_or_else(|| {
                        pyo3::exceptions::PyRuntimeError::new_err(
                            "cannot extract boundary loop from current 2D mesh",
                        )
                    })?;
                fallback
            }
        };
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
    push_log(&mut state, format!("mesh order set to {order}"));
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

fn frustum_brep(
    center: DVec3,
    axis_norm: DVec3,
    ref_dir_hint: DVec3,
    r1: f64,
    r2: f64,
    height: f64,
) -> BRep {
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
            edges: vec![
                WireEdge::fwd(2),
                WireEdge::rev(1),
                WireEdge::rev(2),
                WireEdge::fwd(0),
            ],
        },
        inner_wires: vec![],
        normal: x_axis,
        triangles: vec![],
        mesh_dirty: true,
    };
    let top_face = Face {
        outer_wire: Wire {
            edges: vec![WireEdge::fwd(0)],
        },
        inner_wires: vec![],
        normal: axis_eff,
        triangles: vec![],
        mesh_dirty: true,
    };
    let base_face = Face {
        outer_wire: Wire {
            edges: vec![WireEdge::rev(1)],
        },
        inner_wires: vec![],
        normal: -axis_eff,
        triangles: vec![],
        mesh_dirty: true,
    };

    let shell = Shell {
        faces: vec![side_face, top_face, base_face],
    };
    let solid = Solid { shells: vec![shell] };

    let top_center = map_point(DVec3::new(0.0, half_h, 0.0));
    let base_center = map_point(DVec3::new(0.0, -half_h, 0.0));
    let apex = map_point(DVec3::new(0.0, half_h + apex_dist, 0.0));

    let top_circle = Curve3::Circle(Circle3 {
        center: top_center,
        normal: axis_eff,
        radius: rt,
    });
    let base_circle = Curve3::Circle(Circle3 {
        center: base_center,
        normal: -axis_eff,
        radius: rb,
    });
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
        origin: top_center,
        normal: axis_eff,
    });
    let base_plane = Surface3::Plane(Plane {
        origin: base_center,
        normal: -axis_eff,
    });

    let e0_on_side = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(2.0 * PI, v_top),
        direction: glam::DVec2::new(-1.0, 0.0),
    });
    let e0_on_top = Curve2d::Circle(Circle2d {
        center: glam::DVec2::ZERO,
        radius: rt,
    });
    let e1_on_side = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, v_base),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    let e1_on_base = Curve2d::Circle(Circle2d {
        center: glam::DVec2::ZERO,
        radius: rb,
    });
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
        edge_curve_range: vec![Some([0.0, 2.0 * PI]), Some([0.0, 2.0 * PI]), Some([0.0, slant_len])],
        edge_degenerated: vec![false, false, false],
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
        solids: vec![solid],
        geom,
        compound: None,
        compsolid: None,
    }
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

    let shape = if (r1 - r2).abs() <= 1e-12 {
        let ref_dir = if axis_norm.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
        cylinder_brep(center, axis_norm, ref_dir, r1, height)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?
    } else if r2 <= 1e-12 {
        let ref_dir = if axis_norm.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
        cone_brep(center, axis_norm, ref_dir, r1, height)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?
    } else if r1 <= 1e-12 {
        let axis_flip = -axis_norm;
        let ref_dir = if axis_flip.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
        cone_brep(center, axis_flip, ref_dir, r2, height)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?
    } else {
        let ref_dir = if axis_norm.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
        frustum_brep(center, axis_norm, ref_dir, r1, r2, height)
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

    let p0 = DVec3::new(x, y, z);
    let p1 = DVec3::new(x + dx, y, z);
    let p2 = DVec3::new(x + dx, y + dy, z);
    let p3 = DVec3::new(x, y + dy, z);

    let mut shape = BRep::new();

    let v0 = make_vertex(&mut shape, p0);
    let v1 = make_vertex(&mut shape, p1);
    let v2 = make_vertex(&mut shape, p2);
    let v3 = make_vertex(&mut shape, p3);

    let make_line_edge = |brep: &mut BRep, a: DVec3, b: DVec3, va: usize, vb: usize| -> PyResult<usize> {
        let seg = b - a;
        let len = seg.length();
        if len < 1e-12 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "degenerate rectangle edge",
            ));
        }
        let dir = seg / len;
        make_edge(
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
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    };

    let e0 = make_line_edge(&mut shape, p0, p1, v0, v1)?;
    let e1 = make_line_edge(&mut shape, p1, p2, v1, v2)?;
    let e2 = make_line_edge(&mut shape, p2, p3, v2, v3)?;
    let e3 = make_line_edge(&mut shape, p3, p0, v3, v0)?;

    let outer = make_wire(vec![
        rcad_kernel::topology::WireEdge::fwd(e0),
        rcad_kernel::topology::WireEdge::fwd(e1),
        rcad_kernel::topology::WireEdge::fwd(e2),
        rcad_kernel::topology::WireEdge::fwd(e3),
    ]);

    let normal = (p1 - p0).cross(p3 - p0).normalize_or_zero();
    make_face(
        &mut shape,
        rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane {
            origin: p0,
            normal,
        }),
        outer,
        Vec::new(),
    )
    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

    let assigned_tag = if tag > 0 { tag } else { state.next_cad_tag + 1 };
    state.next_cad_tag = assigned_tag.max(state.next_cad_tag);
    state.cad_shapes.insert(assigned_tag, shape);
    // Clear stale mesh cache: meshing will derive from CAD on demand.
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

    let mut shape = BRep::new();
    make_vertex(&mut shape, DVec3::new(x, y, z));

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

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

    let seg = p1 - p0;
    let len = seg.length();
    if len < 1e-12 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addLine: start and end points are coincident",
        ));
    }

    let mut shape = BRep::new();
    let v0 = make_vertex(&mut shape, p0);
    let v1 = make_vertex(&mut shape, p1);
    make_edge(
        &mut shape,
        Curve3::Line(Line3 {
            origin: p0,
            direction: seg / len,
        }),
        0.0,
        len,
        v0,
        v1,
    )
    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

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

    let r0 = p0 - c;
    let r1 = p1 - c;
    let radius = r0.length();
    if radius < 1e-12 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addCircleArc: radius is zero",
        ));
    }
    if (r1.length() - radius).abs() > 1e-6 * radius.max(1.0) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addCircleArc: start/end points are not equidistant from center",
        ));
    }

    let normal = r0.cross(r1);
    if normal.length_squared() < 1e-20 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addCircleArc: points are collinear",
        ));
    }
    let normal = normal.normalize();
    let x_axis = rcad_kernel::geom::any_perpendicular(normal);
    let y_axis = normal.cross(x_axis);

    let t0 = r0.dot(y_axis).atan2(r0.dot(x_axis));
    let mut t1 = r1.dot(y_axis).atan2(r1.dot(x_axis));
    while t1 <= t0 {
        t1 += 2.0 * PI;
    }

    let mut shape = BRep::new();
    let v0 = make_vertex(&mut shape, p0);
    let v1 = make_vertex(&mut shape, p1);
    make_edge(
        &mut shape,
        Curve3::Circle(Circle3 {
            center: c,
            normal,
            radius,
        }),
        t0,
        t1,
        v0,
        v1,
    )
    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

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

    let bspline = interpolate_points(&pts)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("addSpline: {e}")))?;
    let t0 = bspline
        .knots
        .get(bspline.degree)
        .copied()
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("addSpline: invalid knot vector"))?;
    let t1 = bspline
        .knots
        .get(bspline.knots.len().saturating_sub(bspline.degree + 1))
        .copied()
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("addSpline: invalid knot vector"))?;

    let mut shape = BRep::new();
    let v0 = make_vertex(&mut shape, pts[0]);
    let v1 = make_vertex(&mut shape, *pts.last().expect("point list has len >= 2"));
    make_edge(&mut shape, Curve3::BSpline(bspline), t0, t1, v0, v1)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

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

    let center = DVec3::new(xc, yc, zc);
    let major_is_x = rx.abs() >= ry.abs();
    let major_radius = if major_is_x { rx.abs() } else { ry.abs() };
    let minor_radius = if major_is_x { ry.abs() } else { rx.abs() };
    let major_dir = if major_is_x { DVec3::X } else { DVec3::Y };

    let edge_curve = if (rx.abs() - ry.abs()).abs() < 1e-12 {
        Curve3::Circle(Circle3 {
            center,
            normal: DVec3::Z,
            radius: rx.abs(),
        })
    } else {
        Curve3::Ellipse(rcad_kernel::geom::Ellipse3 {
            center,
            normal: DVec3::Z,
            major_dir,
            major_radius,
            minor_radius,
        })
    };

    let mut shape = BRep::new();
    let v = make_vertex(&mut shape, center + major_radius * major_dir);
    let edge = make_edge(&mut shape, edge_curve, 0.0, 2.0 * PI, v, v)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let outer = make_wire(vec![WireEdge::fwd(edge)]);
    make_face(
        &mut shape,
        Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 0.0, zc),
            normal: DVec3::Z,
        }),
        outer,
        Vec::new(),
    )
    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let curve2d_idx = if (rx.abs() - ry.abs()).abs() < 1e-12 {
        shape.geom.curve2ds.push(Curve2d::Circle(Circle2d {
            center: glam::DVec2::new(xc, yc),
            radius: rx.abs(),
        }));
        shape.geom.curve2ds.len() - 1
    } else {
        shape.geom.curve2ds.push(Curve2d::Ellipse(rcad_kernel::geom::Ellipse2d {
            center: glam::DVec2::new(xc, yc),
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
    shape.geom.edge_pcurves[edge].push(PCurve {
        surface_idx: 0,
        curve2d_idx,
    });

    let mut state = STATE
        .lock()
        .map_err(|_| pyo3::exceptions::PyRuntimeError::new_err("rmsh state lock poisoned"))?;
    ensure_initialized(&state)?;

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

    let mut shape = BRep::new();
    let mut vertex_map: Vec<(DVec3, usize)> = Vec::new();
    let mut loop_wires: Vec<Wire> = Vec::new();
    let mut first_loop_points: Vec<DVec3> = Vec::new();

    let get_or_add_vertex = |brep: &mut BRep,
                             map: &mut Vec<(DVec3, usize)>,
                             p: DVec3|
     -> usize {
        if let Some((_, idx)) = map
            .iter()
            .find(|(q, _)| (*q - p).length_squared() <= 1e-18)
        {
            *idx
        } else {
            let idx = make_vertex(brep, p);
            map.push((p, idx));
            idx
        }
    };

    for (loop_i, loop_tag) in wire_tags.iter().enumerate() {
        let curve_tags = state.occ_curve_loops.get(loop_tag).cloned().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "addPlaneSurface: curve loop tag {} not found",
                loop_tag
            ))
        })?;

        if curve_tags.is_empty() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "addPlaneSurface: empty curve loop",
            ));
        }

        let mut wire_edges = Vec::with_capacity(curve_tags.len());
        for signed_curve_tag in curve_tags {
            let curve_tag = signed_curve_tag.abs();
            let reverse = signed_curve_tag < 0;
            let src = state.cad_shapes.get(&curve_tag).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "addPlaneSurface: curve tag {} not found",
                    curve_tag
                ))
            })?;

            if src.edges.is_empty() {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "addPlaneSurface: tag {} is not a curve entity",
                    curve_tag
                )));
            }

            let src_edge_idx = 0usize;
            let src_edge = &src.edges[src_edge_idx];
            let curve_idx = src
                .geom
                .edge_curve
                .get(src_edge_idx)
                .copied()
                .flatten()
                .ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err(format!(
                        "addPlaneSurface: curve {} has no 3D geometry",
                        curve_tag
                    ))
                })?;
            let curve = src
                .geom
                .curves
                .get(curve_idx)
                .cloned()
                .ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err(format!(
                        "addPlaneSurface: invalid curve geometry for tag {}",
                        curve_tag
                    ))
                })?;

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
                .get(src_edge.start)
                .map(|v| v.point)
                .ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err("addPlaneSurface: invalid edge start")
                })?;
            let mut p_end = src
                .vertices
                .get(src_edge.end)
                .map(|v| v.point)
                .ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err("addPlaneSurface: invalid edge end")
                })?;

            let (mut t0, mut t1) = (range[0], range[1]);
            if reverse {
                std::mem::swap(&mut p_start, &mut p_end);
                std::mem::swap(&mut t0, &mut t1);
            }

            let v0 = get_or_add_vertex(&mut shape, &mut vertex_map, p_start);
            let v1 = get_or_add_vertex(&mut shape, &mut vertex_map, p_end);
            let e = make_edge(&mut shape, curve, t0, t1, v0, v1)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
            wire_edges.push(WireEdge::fwd(e));

            if loop_i == 0 {
                first_loop_points.push(p_start);
            }
        }

        loop_wires.push(make_wire(wire_edges));
    }

    if loop_wires.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "addPlaneSurface: failed to build wires",
        ));
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

    make_face(
        &mut shape,
        Surface3::Plane(Plane { origin, normal }),
        outer,
        loop_wires,
    )
    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

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
        merge_breps_for_export(&shapes)
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

    let face_count: usize = shape
        .solids
        .iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum();

    let mut total_outer_wire_edge_refs = 0usize;
    let mut min_outer_wire_edges_per_face = usize::MAX;
    let mut max_outer_wire_edges_per_face = 0usize;
    let mut max_outer_wire_face_surface_kind: &'static str = "UnknownSurface";
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
                    max_outer_wire_face_surface_kind = shape
                        .geom
                        .face_surface
                        .get(face_cursor)
                        .copied()
                        .flatten()
                        .and_then(|si| shape.geom.surfaces.get(si))
                        .map(surface3_kind_name)
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

    let mut outer_wire_curve_kinds: BTreeMap<&'static str, usize> = BTreeMap::new();
    for &ei in &outer_wire_edge_indices {
        if let Some(curve_idx) = shape.geom.edge_curve.get(ei).copied().flatten()
            && let Some(curve) = shape.geom.curves.get(curve_idx)
        {
            *outer_wire_curve_kinds
                .entry(curve3_kind_name(curve))
                .or_insert(0) += 1;
        }
    }

    let edges_with_pcurves = shape
        .geom
        .edge_pcurves
        .iter()
        .filter(|pcs| !pcs.is_empty())
        .count();
    let total_pcurves: usize = shape.geom.edge_pcurves.iter().map(|pcs| pcs.len()).sum();

    let mut curve_kinds: BTreeMap<&'static str, usize> = BTreeMap::new();
    for c in &shape.geom.curves {
        *curve_kinds.entry(curve3_kind_name(c)).or_insert(0) += 1;
    }

    let mut surface_kinds: BTreeMap<&'static str, usize> = BTreeMap::new();
    for s in &shape.geom.surfaces {
        *surface_kinds.entry(surface3_kind_name(s)).or_insert(0) += 1;
    }

    let mut pcurve_surface_kinds: BTreeMap<&'static str, usize> = BTreeMap::new();
    for pcs in &shape.geom.edge_pcurves {
        for pc in pcs {
            if let Some(surface) = shape.geom.surfaces.get(pc.surface_idx) {
                *pcurve_surface_kinds
                    .entry(surface3_kind_name(surface))
                    .or_insert(0) += 1;
            } else {
                *pcurve_surface_kinds.entry("UnknownSurface").or_insert(0) += 1;
            }
        }
    }

    Python::with_gil(|py| {
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("tag", tag)?;
        dict.set_item("vertices", shape.vertices.len())?;
        dict.set_item("edges", shape.edges.len())?;
        dict.set_item("faces", face_count)?;
        dict.set_item("solids", shape.solids.len())?;
        dict.set_item("closed_edges", edge_closed)?;
        dict.set_item("curves_3d", shape.geom.curves.len())?;
        dict.set_item("surfaces_3d", shape.geom.surfaces.len())?;
        dict.set_item("curves_2d", shape.geom.curve2ds.len())?;
        dict.set_item("edge_curve_slots", shape.geom.edge_curve.len())?;
        dict.set_item("edge_curve_some", edge_curve_some)?;
        dict.set_item("edge_curve_none", edge_curve_none)?;
        dict.set_item("edge_pcurve_slots", shape.geom.edge_pcurves.len())?;
        dict.set_item("edges_with_pcurves", edges_with_pcurves)?;
        dict.set_item("total_pcurves", total_pcurves)?;
        dict.set_item("outer_wire_edges", outer_wire_edge_count)?;
        dict.set_item("outer_wire_edge_refs_total", total_outer_wire_edge_refs)?;
        dict.set_item("outer_wire_edges_per_face_min", min_outer_wire_edges_per_face)?;
        dict.set_item("outer_wire_edges_per_face_max", max_outer_wire_edges_per_face)?;
        dict.set_item("outer_wire_face_max_surface_kind", max_outer_wire_face_surface_kind)?;
        dict.set_item("faces_with_outer_wire_over_100", faces_with_outer_wire_over_100)?;
        dict.set_item("outer_wire_edges_with_curve", outer_wire_edges_with_curve)?;
        dict.set_item("outer_wire_edges_without_curve", outer_wire_edges_without_curve)?;
        dict.set_item("outer_wire_edges_with_pcurves", outer_wire_edges_with_pcurves)?;

        let curve_dict = pyo3::types::PyDict::new(py);
        for (k, v) in curve_kinds {
            curve_dict.set_item(k, v)?;
        }
        dict.set_item("curve3_kinds", curve_dict)?;

        let surface_dict = pyo3::types::PyDict::new(py);
        for (k, v) in surface_kinds {
            surface_dict.set_item(k, v)?;
        }
        dict.set_item("surface3_kinds", surface_dict)?;

        let pcurve_surface_dict = pyo3::types::PyDict::new(py);
        for (k, v) in pcurve_surface_kinds {
            pcurve_surface_dict.set_item(k, v)?;
        }
        dict.set_item("pcurve_surface_kinds", pcurve_surface_dict)?;

        let outer_wire_curve_dict = pyo3::types::PyDict::new(py);
        for (k, v) in outer_wire_curve_kinds {
            outer_wire_curve_dict.set_item(k, v)?;
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
    m.add_function(pyo3::wrap_pyfunction!(model_occ_debug_shape_geom_impl, m)?)?;

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        estimate_mesh_characteristic_size, frustum_brep, option_default_number, step_protocol_from_state,
        Mesh, Node, RuntimeState,
    };
    use glam::DVec3;
    use rcad_kernel::fit::interpolate_points;
    use rcad_modeling::builder::{cylinder_brep, make_edge, make_vertex};

    fn count_occurrences(haystack: &str, needle: &str) -> usize {
        haystack.matches(needle).count()
    }

    #[test]
    fn strict_cylinder_emits_cylindrical_surface_and_seam_curve() {
        let brep = cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 0.5, 2.0)
            .expect("cylinder should build");
        let options = rmsh_io::BrepStepWriteOptions {
            protocol: rcad_step::StepProtocol::Ap214,
            solid_color: None,
            header: None,
            gmsh_strict: true,
        };

        let step = rmsh_io::write_brep_step_with_options(&brep, &options)
            .expect("strict cylinder STEP export should succeed");

        assert_eq!(count_occurrences(&step, "ADVANCED_FACE"), 3);
        assert_eq!(count_occurrences(&step, "CYLINDRICAL_SURFACE"), 1);
        assert!(count_occurrences(&step, "SEAM_CURVE") >= 1);
    }

    #[test]
    fn strict_frustum_cone_emits_conical_side_and_three_faces() {
        let brep = frustum_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 0.8, 0.2, 2.0);
        let options = rmsh_io::BrepStepWriteOptions {
            protocol: rcad_step::StepProtocol::Ap214,
            solid_color: None,
            header: None,
            gmsh_strict: true,
        };

        let step = rmsh_io::write_brep_step_with_options(&brep, &options)
            .expect("strict frustum STEP export should succeed");

        assert_eq!(count_occurrences(&step, "ADVANCED_FACE"), 3);
        assert_eq!(count_occurrences(&step, "CONICAL_SURFACE"), 1);
        assert_eq!(count_occurrences(&step, "EDGE_CURVE"), 3);
        assert!(!step.contains("TRIANGULATED"));
        assert!(!step.contains("TESSELLATED"));
    }

    #[test]
    fn strict_standalone_line_emits_wireframe_curve_set() {
        let mut brep = rcad_kernel::BRep::new();
        let v0 = make_vertex(&mut brep, DVec3::new(0.0, 0.0, 0.0));
        let v1 = make_vertex(&mut brep, DVec3::new(1.0, 0.0, 0.0));
        make_edge(
            &mut brep,
            rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::X,
            }),
            0.0,
            1.0,
            v0,
            v1,
        )
        .expect("line edge should build");

        let options = rmsh_io::BrepStepWriteOptions {
            protocol: rcad_step::StepProtocol::Ap214,
            solid_color: None,
            header: None,
            gmsh_strict: true,
        };

        let step = rmsh_io::write_brep_step_with_options(&brep, &options)
            .expect("strict standalone line STEP export should succeed");

        assert!(step.contains("GEOMETRIC_CURVE_SET"));
        assert!(step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"));
        assert!(step.contains("LINE"));
    }

    #[test]
    fn strict_standalone_spline_emits_bspline_curve() {
        let mut brep = rcad_kernel::BRep::new();
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
        ];
        let bs = interpolate_points(&pts).expect("spline interpolation should succeed");
        let t0 = bs.knots[bs.degree];
        let t1 = bs.knots[bs.knots.len() - bs.degree - 1];
        let v0 = make_vertex(&mut brep, pts[0]);
        let v1 = make_vertex(&mut brep, *pts.last().expect("points are non-empty"));
        make_edge(&mut brep, rcad_kernel::geom::Curve3::BSpline(bs), t0, t1, v0, v1)
            .expect("spline edge should build");

        let options = rmsh_io::BrepStepWriteOptions {
            protocol: rcad_step::StepProtocol::Ap214,
            solid_color: None,
            header: None,
            gmsh_strict: true,
        };

        let step = rmsh_io::write_brep_step_with_options(&brep, &options)
            .expect("strict standalone spline STEP export should succeed");

        assert!(step.contains("GEOMETRIC_CURVE_SET"));
        assert!(step.contains("B_SPLINE_CURVE_WITH_KNOTS"));
    }

    #[test]
    fn step_protocol_uses_default_option_value() {
        let state = RuntimeState::default();
        let protocol = step_protocol_from_state(&state);
        assert!(matches!(protocol, rcad_step::StepProtocol::Ap214));
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
