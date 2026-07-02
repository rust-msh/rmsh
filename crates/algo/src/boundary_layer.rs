//! Boundary layer mesh generation — prism/hex extrusion from surface.
//!
//! Given a closed surface mesh, extrude a layered prismatic boundary layer
//! along the surface normal, then fill the remaining interior with a tetrahedral
//! mesh.  This is the standard "mesh + BL + fill" workflow for CFD.
//!
//! # Algorithm
//!
//! 1. Compute smoothed vertex normals (area-weighted average of adjacent
//!    face normals).
//! 2. For each layer `k = 0..n-1`, offset vertices by `h_k * normal`,
//!    where `h_k` follows a geometric progression.
//! 3. Generate prism elements connecting layer `k` to layer `k+1`.
//! 4. The offset surface at layer `n-1` becomes input for volume meshing.
//!
//! # Gmsh counterpart
//!
//! `Mesh.BoundaryLayerField` + `Field.BoundaryLayer`.
//! Also `model.mesh.set_extrude_layers`.

use std::collections::HashMap;

use rmsh_model::{Element, ElementType, Mesh, Node};

// ─── Error ────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub enum BlError {
    InvalidInput(String),
    Generation(String),
}

impl std::fmt::Display for BlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BlError::InvalidInput(msg) => write!(f, "boundary layer: {msg}"),
            BlError::Generation(msg) => write!(f, "boundary layer generation: {msg}"),
        }
    }
}

// ─── Parameters ───────────────────────────────────────────────────────────────

/// Parameters controlling boundary layer generation.
#[derive(Debug, Clone)]
pub struct BoundaryLayerParams {
    /// Number of layers.
    pub num_layers: usize,
    /// Height of the first layer (adjacent to the wall).
    pub first_height: f64,
    /// Geometric progression ratio: h_{k+1} = h_k * growth_rate.
    /// Typical values: 1.1 – 1.3.
    pub growth_rate: f64,
    /// When true, the boundary layer is extruded outward (external flow).
    /// When false, inward (internal flow / domain interior).
    pub outward: bool,
}

impl Default for BoundaryLayerParams {
    fn default() -> Self {
        Self {
            num_layers: 5,
            first_height: 0.01,
            growth_rate: 1.2,
            outward: true,
        }
    }
}

// ─── Core function ────────────────────────────────────────────────────────────

/// Generate a boundary layer around `surface` and return the combined mesh
/// (prism layers + original surface).
///
/// `surface` must contain Triangle3 elements forming a closed (or oriented)
/// surface.  The function works on any triangle mesh — the normals are
/// computed from the mesh geometry, not from CAD topology.
///
/// Returns a mesh containing:
/// - All original surface nodes.
/// - `num_layers × n_triangles` prism elements (Triangle3 → Prism6).
/// - New offset nodes at each layer.
pub fn generate_boundary_layer(
    surface: &Mesh,
    params: &BoundaryLayerParams,
) -> Result<Mesh, BlError> {
    if params.num_layers == 0 {
        return Err(BlError::InvalidInput("num_layers must be >= 1".into()));
    }

    // Collect all triangles (Tri3 elements).
    let tri_indices: Vec<usize> = surface
        .elements
        .iter()
        .enumerate()
        .filter(|(_, e)| e.etype == ElementType::Triangle3 && e.node_ids.len() == 3)
        .map(|(i, _)| i)
        .collect();

    if tri_indices.is_empty() {
        // Try getting triangles from Quad4 elements by splitting.
        let quad_indices: Vec<usize> = surface
            .elements
            .iter()
            .enumerate()
            .filter(|(_, e)| e.etype == ElementType::Quad4 && e.node_ids.len() == 4)
            .map(|(i, _)| i)
            .collect();
        if quad_indices.is_empty() {
            return Err(BlError::InvalidInput(
                "surface mesh must contain Triangle3 or Quad4 elements".into(),
            ));
        }
        // Return error for now — Quad4 splitting will be handled in a future update.
        return Err(BlError::InvalidInput(
            "Quad4 surfaces not yet supported; triangulate first".into(),
        ));
    }

    // Step 1: compute area-weighted vertex normals.
    let vertex_normals = compute_vertex_normals(surface, &tri_indices)?;

    // Step 2: for each layer, compute node offsets.
    // Compute total thickness per layer.
    let mut layer_heights: Vec<f64> = Vec::with_capacity(params.num_layers);
    let mut h = params.first_height;
    for _ in 0..params.num_layers {
        layer_heights.push(h);
        h *= params.growth_rate;
    }

    // Cumulative offset distance for the top of each layer.
    let mut cumulative: Vec<f64> = Vec::with_capacity(params.num_layers + 1);
    cumulative.push(0.0);
    for hk in &layer_heights {
        cumulative.push(cumulative.last().unwrap() + hk);
    }

    // Step 3: build the output mesh.
    // Strategy: create a layer mesh for each layer, then combine.
    let mut out = Mesh::new();

    // Copy all nodes from the original surface.
    let mut node_map: HashMap<u64, Vec<u64>> = HashMap::new();
    // node_map[original_id] = [layer_0_id, layer_1_id, ..., layer_n_id]
    // layer_0_id = original node (not offset)

    for (nid, node) in &surface.nodes {
        let mut layer_ids = Vec::with_capacity(params.num_layers + 1);
        // Layer 0: original node.
        layer_ids.push(*nid);
        out.add_node(Node::new(*nid, node.position.x, node.position.y, node.position.z));
        node_map.insert(*nid, layer_ids);
    }

    let mut next_id = surface.nodes.keys().copied().max().unwrap_or(0).saturating_add(1);

    // For each layer 1..num_layers, create offset nodes.
    let sign = if params.outward { 1.0 } else { -1.0 };
    for layer in 1..=params.num_layers {
        let offset_dist = cumulative[layer]; // total offset at this layer
        for &nid in surface.nodes.keys() {
            let normal = vertex_normals.get(&nid).copied().unwrap_or([0.0, 0.0, 0.0]);
            let node = &surface.nodes[&nid];
            let xn = node.position.x + sign * offset_dist * normal[0];
            let yn = node.position.y + sign * offset_dist * normal[1];
            let zn = node.position.z + sign * offset_dist * normal[2];
            let new_id = next_id;
            next_id += 1;
            out.add_node(Node::new(new_id, xn, yn, zn));
            node_map.get_mut(&nid).unwrap().push(new_id);
        }
    }

    // Step 4: generate prism elements between consecutive layers.
    let mut eid = surface.elements.iter().map(|e| e.id).max().unwrap_or(0).saturating_add(1);

    for layer in 0..params.num_layers {
        for &ti in &tri_indices {
            let elt = &surface.elements[ti];
            let n = &elt.node_ids;
            // Original nodes: n[0], n[1], n[2]
            // Layer-lower nodes: node_map[&n[0]][layer], ...
            // Layer-upper nodes: node_map[&n[0]][layer+1], ...
            let l0 = |idx: usize| node_map[&n[idx]][layer];
            let l1 = |idx: usize| node_map[&n[idx]][layer + 1];

            let lower = [l0(0), l0(1), l0(2)];
            let upper = [l1(0), l1(1), l1(2)];

            // Prism6 node ordering: bottom triangle ccw, top triangle ccw.
            out.add_element(Element::new(eid, ElementType::Prism6, vec![
                lower[0], lower[1], lower[2],
                upper[0], upper[1], upper[2],
            ]));
            eid += 1;
        }
    }

    Ok(out)
}

// ─── Normal computation ───────────────────────────────────────────────────────

/// Compute area-weighted vertex normals for all nodes referenced by triangles.
///
/// Returns a map from node ID to unit normal vector.
fn compute_vertex_normals(
    mesh: &Mesh,
    tri_indices: &[usize],
) -> Result<HashMap<u64, [f64; 3]>, BlError> {
    // First pass: accumulate face normals × area for each vertex.
    let mut normals: HashMap<u64, [f64; 3]> = HashMap::new();
    let mut areas: HashMap<u64, f64> = HashMap::new();

    for &ti in tri_indices {
        let elt = &mesh.elements[ti];
        let n = &elt.node_ids;
        let get = |i: usize| -> Result<[f64; 3], BlError> {
            let p = mesh.nodes.get(&n[i])
                .ok_or_else(|| BlError::Generation(format!("missing node {}", n[i])))?;
            Ok([p.position.x, p.position.y, p.position.z])
        };
        let a = get(0)?;
        let b = get(1)?;
        let c = get(2)?;

        // Face normal = (b-a) × (c-a)
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let mut nrm = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        let area = 0.5 * (nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]).sqrt();
        if area < 1e-20 {
            continue;
        }
        // Weight = area (larger faces contribute more to vertex normal).
        for &vid in &[n[0], n[1], n[2]] {
            let entry = normals.entry(vid).or_insert([0.0, 0.0, 0.0]);
            entry[0] += nrm[0];
            entry[1] += nrm[1];
            entry[2] += nrm[2];
            *areas.entry(vid).or_insert(0.0) += area;
        }
    }

    // Second pass: normalise each vertex normal.
    for (_, nrm) in normals.iter_mut() {
        let len = (nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]).sqrt();
        if len > 1e-15 {
            nrm[0] /= len;
            nrm[1] /= len;
            nrm[2] /= len;
        }
    }

    Ok(normals)
}

// ─── Interior intersection helper ────────────────────────────────────────────

/// Given a surface mesh with a boundary layer and an inner point known to be
/// inside the volume, compute the signed distance from the boundary layer
/// to the inner point.  This is used to determine the `outward` direction.
///
/// Returns `true` if the average normal points toward the inner point (meaning
/// outward = false for internal flow).
pub fn detect_inward_normals(surface: &Mesh, inner_point: [f64; 3]) -> Result<bool, BlError> {
    let n_tris = surface.elements.len();
    if n_tris == 0 {
        return Err(BlError::InvalidInput("empty surface mesh".into()));
    }

    let mut avg_normal = [0.0_f64; 3];
    let mut count = 0usize;

    for elt in &surface.elements {
        if elt.etype != ElementType::Triangle3 || elt.node_ids.len() != 3 {
            continue;
        }
        let n = &elt.node_ids;
        let get = |i: usize| -> Option<[f64; 3]> {
            let p = surface.nodes.get(&n[i])?;
            Some([p.position.x, p.position.y, p.position.z])
        };
        let (a, b, c) = match (get(0), get(1), get(2)) {
            (Some(a), Some(b), Some(c)) => (a, b, c),
            _ => continue,
        };
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let nrm = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        let len = (nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]).sqrt();
        if len < 1e-20 {
            continue;
        }
        avg_normal[0] += nrm[0] / len;
        avg_normal[1] += nrm[1] / len;
        avg_normal[2] += nrm[2] / len;
        count += 1;
    }

    if count == 0 {
        return Err(BlError::Generation("no valid triangles in surface".into()));
    }

    let len = (avg_normal[0] * avg_normal[0] + avg_normal[1] * avg_normal[1] + avg_normal[2] * avg_normal[2]).sqrt();
    if len < 1e-15 {
        return Err(BlError::Generation("surface has zero average normal".into()));
    }
    avg_normal[0] /= len;
    avg_normal[1] /= len;
    avg_normal[2] /= len;

    // Compute centroid of the surface.
    let mut cx = 0.0_f64; let mut cy = 0.0_f64; let mut cz = 0.0_f64;
    for node in surface.nodes.values() {
        cx += node.position.x;
        cy += node.position.y;
        cz += node.position.z;
    }
    let n_nodes = surface.nodes.len() as f64;
    cx /= n_nodes; cy /= n_nodes; cz /= n_nodes;

    // Vector from centroid to inner_point.
    let to_inner = [inner_point[0] - cx, inner_point[1] - cy, inner_point[2] - cz];
    let dot = to_inner[0] * avg_normal[0] + to_inner[1] * avg_normal[1] + to_inner[2] * avg_normal[2];

    // If dot > 0, the average normal points toward the inner point (normals point inward).
    Ok(dot > 0.0)
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rmsh_model::{Element, ElementType, Mesh, Node};

    /// A flat square surface (2 triangles) at z=0, pointing upward.
    fn flat_surface() -> Mesh {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 1.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 1.0, 0.0));
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));
        mesh.add_element(Element::new(2, ElementType::Triangle3, vec![1, 3, 4]));
        mesh
    }

    #[test]
    fn vertex_normals_flat_surface() {
        let mesh = flat_surface();
        let tris: Vec<usize> = (0..mesh.elements.len()).collect();
        let normals = compute_vertex_normals(&mesh, &tris).unwrap();
        for (_, n) in &normals {
            assert!((n[0]).abs() < 1e-10, "x component should be near 0: {}", n[0]);
            assert!((n[1]).abs() < 1e-10, "y component should be near 0: {}", n[1]);
            assert!((n[2]).abs() > 0.99, "z component should be ~±1: {}", n[2]);
        }
    }

    #[test]
    fn generate_bl_simple() {
        let mesh = flat_surface();
        let params = BoundaryLayerParams {
            num_layers: 3,
            first_height: 0.1,
            growth_rate: 1.5,
            outward: true,
        };
        let bl = generate_boundary_layer(&mesh, &params).unwrap();
        // Each original node gets `num_layers` offset copies (plus original).
        // Total nodes: 4 + 4*3 = 16
        assert_eq!(bl.nodes.len(), 4 + 4 * 3);
        // Each triangle × num_layers = 2 × 3 = 6 prism elements
        assert_eq!(bl.elements.len(), 6);
        for elt in &bl.elements {
            assert_eq!(elt.etype, ElementType::Prism6);
            assert_eq!(elt.node_ids.len(), 6);
        }
    }

    #[test]
    fn generate_bl_zero_layers_errors() {
        let mesh = flat_surface();
        let params = BoundaryLayerParams {
            num_layers: 0,
            ..Default::default()
        };
        let result = generate_boundary_layer(&mesh, &params);
        assert!(result.is_err());
    }

    #[test]
    fn generate_bl_no_triangles_errors() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        let result = generate_boundary_layer(&mesh, &BoundaryLayerParams::default());
        assert!(result.is_err());
    }

    #[test]
    fn detect_inward_normal_flat() {
        let mesh = flat_surface();
        // Inner point above the surface → normals point toward it → inward detected.
        let inward = detect_inward_normals(&mesh, [0.5, 0.5, 1.0]).unwrap();
        assert!(inward, "normals should point toward inner point (above z=0)");
    }

    #[test]
    fn detect_outward_normal_flat() {
        let mesh = flat_surface();
        // Inner point below the surface → normals point away from it → not inward.
        let inward = detect_inward_normals(&mesh, [0.5, 0.5, -1.0]).unwrap();
        assert!(!inward, "normals should point away from inner point (below z=0)");
    }

    #[test]
    fn bl_outward_offset_z() {
        let mesh = flat_surface();
        let params = BoundaryLayerParams {
            num_layers: 2,
            first_height: 0.5,
            growth_rate: 1.0,
            outward: true,
        };
        let bl = generate_boundary_layer(&mesh, &params).unwrap();
        // Total nodes = 4 surface + 4×2 offset = 12
        assert_eq!(bl.nodes.len(), 12);
        // Node IDs 1-4: original (z=0). Check by scanning.
        let mut count_z0 = 0usize;
        let mut count_z05 = 0usize;
        let mut count_z1 = 0usize;
        for node in bl.nodes.values() {
            let z = node.position.z;
            if (z - 0.0).abs() < 1e-12 { count_z0 += 1; }
            if (z - 0.5).abs() < 1e-12 { count_z05 += 1; }
            if (z - 1.0).abs() < 1e-12 { count_z1 += 1; }
        }
        assert_eq!(count_z0, 4, "expected 4 nodes at z=0");
        assert_eq!(count_z05, 4, "expected 4 nodes at z=0.5");
        assert_eq!(count_z1, 4, "expected 4 nodes at z=1.0");
    }
}
