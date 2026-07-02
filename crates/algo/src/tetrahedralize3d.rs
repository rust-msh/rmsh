use std::collections::{BTreeMap, HashSet};

use rmsh_model::{Element, ElementType, Mesh, Node, Point3};
use thiserror::Error;

use crate::traits::{MeshAlgoError, MeshParams, Mesher3D};

#[derive(Error, Debug)]
pub enum Mesh3DError {
    #[error("3D meshing failed: {0}")]
    Generation(String),
}

impl From<Mesh3DError> for MeshAlgoError {
    fn from(value: Mesh3DError) -> Self {
        match value {
            Mesh3DError::Generation(msg) => MeshAlgoError::Generation(msg),
        }
    }
}

/// Centroid-star tetrahedralization algorithm.
///
/// This is a lightweight volume mesher that creates one interior point at the
/// centroid of all boundary nodes and then connects each boundary triangle to
/// that point. It is robust and simple, and works well for star-shaped domains.
#[derive(Debug, Default, Clone, Copy)]
pub struct CentroidStarMesher3D;

impl Mesher3D for CentroidStarMesher3D {
    fn name(&self) -> &'static str {
        "Centroid Star 3D"
    }

    fn mesh_3d(&self, surface: &Mesh, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        if !params.element_size.is_finite() || params.element_size <= 0.0 {
            return Err(MeshAlgoError::InvalidInput(
                "element_size must be a positive finite value".to_string(),
            ));
        }
        tetrahedralize_closed_surface(surface).map_err(Into::into)
    }
}

/// Build a simple tetrahedral volume mesh from a closed boundary surface.
///
/// The method is intentionally simple:
/// - Build boundary polygons from either 3D element boundary faces or existing 2D elements.
/// - Triangulate polygons with a fan split.
/// - Add one interior node at the boundary centroid.
/// - Create one tetrahedron per boundary triangle: `(a, b, c, centroid)`.
///
/// This works best for star-shaped or convex-like closed surfaces.
pub fn tetrahedralize_closed_surface(mesh: &Mesh) -> Result<Mesh, Mesh3DError> {
    let boundary_polys = collect_boundary_polygons(mesh)?;
    if boundary_polys.is_empty() {
        return Err(Mesh3DError::Generation(
            "No boundary polygons found for 3D meshing".to_string(),
        ));
    }

    let mut boundary_tris: Vec<[u64; 3]> = Vec::new();
    for poly in &boundary_polys {
        if poly.len() < 3 {
            continue;
        }
        for i in 1..(poly.len() - 1) {
            boundary_tris.push([poly[0], poly[i], poly[i + 1]]);
        }
    }
    if boundary_tris.is_empty() {
        return Err(Mesh3DError::Generation(
            "Boundary triangulation produced no triangles".to_string(),
        ));
    }

    let mut boundary_nodes: HashSet<u64> = HashSet::new();
    for tri in &boundary_tris {
        boundary_nodes.insert(tri[0]);
        boundary_nodes.insert(tri[1]);
        boundary_nodes.insert(tri[2]);
    }

    let centroid = centroid_of_nodes(mesh, &boundary_nodes)?;
    let centroid_id = mesh
        .nodes
        .keys()
        .copied()
        .max()
        .unwrap_or(0)
        .saturating_add(1);

    let mut out = Mesh::new();
    for nid in &boundary_nodes {
        let node = mesh.nodes.get(nid).ok_or_else(|| {
            Mesh3DError::Generation(format!("Node {} missing in source mesh", nid))
        })?;
        out.add_node(node.clone());
    }
    out.add_node(Node::new(centroid_id, centroid.x, centroid.y, centroid.z));

    let mut eid: u64 = 1;
    let mut accepted = 0usize;
    for tri in &boundary_tris {
        let pa = out
            .nodes
            .get(&tri[0])
            .ok_or_else(|| Mesh3DError::Generation("Missing tetra node A".to_string()))?
            .position;
        let pb = out
            .nodes
            .get(&tri[1])
            .ok_or_else(|| Mesh3DError::Generation("Missing tetra node B".to_string()))?
            .position;
        let pc = out
            .nodes
            .get(&tri[2])
            .ok_or_else(|| Mesh3DError::Generation("Missing tetra node C".to_string()))?
            .position;

        let vol6 = tetra_signed_volume6(pa, pb, pc, centroid).abs();
        if vol6 < 1e-14 {
            continue;
        }

        out.add_element(Element::new(
            eid,
            ElementType::Tetrahedron4,
            vec![tri[0], tri[1], tri[2], centroid_id],
        ));
        eid += 1;
        accepted += 1;
    }

    if accepted == 0 {
        return Err(Mesh3DError::Generation(
            "Generated tetrahedra are degenerate".to_string(),
        ));
    }

    Ok(out)
}

fn collect_boundary_polygons(mesh: &Mesh) -> Result<Vec<Vec<u64>>, Mesh3DError> {
    let mut has_volume = false;
    let mut face_counts: BTreeMap<Vec<u64>, (usize, Vec<u64>)> = BTreeMap::new();

    for elem in &mesh.elements {
        if elem.dimension() != 3 {
            continue;
        }
        has_volume = true;

        let faces = elem.etype.faces();
        if faces.is_empty() {
            continue;
        }

        for local_face in faces {
            let mut face_nodes = Vec::with_capacity(local_face.len());
            for li in *local_face {
                let Some(nid) = elem.node_ids.get(*li) else {
                    return Err(Mesh3DError::Generation(format!(
                        "Element {} face connectivity is inconsistent",
                        elem.id
                    )));
                };
                face_nodes.push(*nid);
            }

            let mut key = face_nodes.clone();
            key.sort_unstable();
            let entry = face_counts.entry(key).or_insert((0usize, face_nodes));
            entry.0 += 1;
        }
    }

    if has_volume {
        let polys = face_counts
            .into_values()
            .filter_map(|(count, face)| (count == 1).then_some(face))
            .collect::<Vec<_>>();
        return Ok(polys);
    }

    let polys = mesh
        .elements
        .iter()
        .filter(|e| e.dimension() == 2 && e.node_ids.len() >= 3)
        .map(|e| e.node_ids.clone())
        .collect::<Vec<_>>();
    Ok(polys)
}

fn centroid_of_nodes(mesh: &Mesh, ids: &HashSet<u64>) -> Result<Point3, Mesh3DError> {
    if ids.is_empty() {
        return Err(Mesh3DError::Generation(
            "Cannot compute centroid from an empty node set".to_string(),
        ));
    }

    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sz = 0.0;
    let mut n = 0usize;
    for nid in ids {
        let node = mesh
            .nodes
            .get(nid)
            .ok_or_else(|| Mesh3DError::Generation(format!("Node {} not found", nid)))?;
        sx += node.position.x;
        sy += node.position.y;
        sz += node.position.z;
        n += 1;
    }

    let inv = 1.0 / n as f64;
    Ok(Point3::new(sx * inv, sy * inv, sz * inv))
}

fn tetra_signed_volume6(a: Point3, b: Point3, c: Point3, d: Point3) -> f64 {
    let ab = b - a;
    let ac = c - a;
    let ad = d - a;
    ab.cross(&ac).dot(&ad)
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;
    use crate::traits::{MeshParams, Mesher3D};
    use rmsh_io::{load_msh_from_bytes, load_step_from_bytes, write_msh_v2, write_msh_v4};

    fn step_file_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("..")
            .join("testdata")
            .join(name)
    }

    #[test]
    fn tetrahedralize_cube_surface() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 1.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(5, 0.0, 0.0, 1.0));
        mesh.add_node(Node::new(6, 1.0, 0.0, 1.0));
        mesh.add_node(Node::new(7, 1.0, 1.0, 1.0));
        mesh.add_node(Node::new(8, 0.0, 1.0, 1.0));

        // 6 quad faces
        mesh.add_element(Element::new(1, ElementType::Quad4, vec![1, 2, 3, 4]));
        mesh.add_element(Element::new(2, ElementType::Quad4, vec![5, 6, 7, 8]));
        mesh.add_element(Element::new(3, ElementType::Quad4, vec![1, 2, 6, 5]));
        mesh.add_element(Element::new(4, ElementType::Quad4, vec![2, 3, 7, 6]));
        mesh.add_element(Element::new(5, ElementType::Quad4, vec![3, 4, 8, 7]));
        mesh.add_element(Element::new(6, ElementType::Quad4, vec![4, 1, 5, 8]));

        let out = tetrahedralize_closed_surface(&mesh).expect("tetrahedralization should succeed");
        assert_eq!(out.elements_by_dimension(3).len(), 12);
        assert!(out.node_count() >= 9);
    }

    #[test]
    fn centroid_star_mesher_trait_path_works() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 1.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(5, 0.0, 0.0, 1.0));
        mesh.add_node(Node::new(6, 1.0, 0.0, 1.0));
        mesh.add_node(Node::new(7, 1.0, 1.0, 1.0));
        mesh.add_node(Node::new(8, 0.0, 1.0, 1.0));

        mesh.add_element(Element::new(1, ElementType::Quad4, vec![1, 2, 3, 4]));
        mesh.add_element(Element::new(2, ElementType::Quad4, vec![5, 6, 7, 8]));
        mesh.add_element(Element::new(3, ElementType::Quad4, vec![1, 2, 6, 5]));
        mesh.add_element(Element::new(4, ElementType::Quad4, vec![2, 3, 7, 6]));
        mesh.add_element(Element::new(5, ElementType::Quad4, vec![3, 4, 8, 7]));
        mesh.add_element(Element::new(6, ElementType::Quad4, vec![4, 1, 5, 8]));

        let params = MeshParams::with_size(0.5);
        let mesher = CentroidStarMesher3D;
        let out = mesher
            .mesh_3d(&mesh, &params)
            .expect("Mesher3D trait pipeline should succeed");

        assert_eq!(mesher.name(), "Centroid Star 3D");
        assert_eq!(out.elements_by_dimension(3).len(), 12);
        assert_eq!(out.node_count(), 9);
    }

    #[test]
    fn centroid_star_mesher_rejects_invalid_params() {
        let mesher = CentroidStarMesher3D;
        let mesh = Mesh::new();
        let bad = MeshParams {
            element_size: 0.0,
            min_size: 0.0,
            max_size: 0.0,
            optimize_passes: 0,
            size_field: None,
        };

        let err = mesher
            .mesh_3d(&mesh, &bad)
            .expect_err("invalid params should return error");

        match err {
            MeshAlgoError::InvalidInput(msg) => {
                assert!(msg.contains("element_size"));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn tetrahedralize_rejects_empty_input() {
        let mesh = Mesh::new();
        let err = tetrahedralize_closed_surface(&mesh).expect_err("empty mesh should fail");
        match err {
            Mesh3DError::Generation(msg) => {
                assert!(msg.contains("No boundary polygons"));
            }
        }
    }

    #[test]
    fn tetrahedralize_reports_missing_nodes_in_surface_elements() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        // Node 99 is missing on purpose.
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 99]));

        let err = tetrahedralize_closed_surface(&mesh).expect_err("should fail on missing nodes");
        match err {
            Mesh3DError::Generation(msg) => {
                assert!(msg.contains("Node 99"));
            }
        }
    }

    #[test]
    fn tetrahedralize_rejects_degenerate_generated_tets() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 1.0, 0.0));
        // Single planar triangle surface => centroid is coplanar with triangle.
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));

        let err = tetrahedralize_closed_surface(&mesh)
            .expect_err("planar open surface should produce only degenerate tets");
        match err {
            Mesh3DError::Generation(msg) => {
                assert!(msg.contains("degenerate"));
            }
        }
    }

    #[test]
    fn collect_boundary_polygons_from_volume_mesh_keeps_only_boundary_faces() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 0.0, 1.0));
        mesh.add_node(Node::new(5, 1.0, 1.0, 1.0));

        // Two tets sharing one face [1,2,3].
        mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));
        mesh.add_element(Element::new(2, ElementType::Tetrahedron4, vec![1, 2, 3, 5]));

        let polys = collect_boundary_polygons(&mesh).expect("boundary extraction should succeed");
        // 2 tetrahedra * 4 faces - 2 shared copies of one face = 6 boundary faces
        assert_eq!(polys.len(), 6);
        // The shared face should not be part of the boundary.
        let mut shared = vec![1, 2, 3];
        shared.sort_unstable();
        assert!(!polys.iter().any(|f| {
            let mut k = f.clone();
            k.sort_unstable();
            k == shared
        }));
    }

    #[test]
    fn collect_boundary_polygons_reports_inconsistent_volume_face_connectivity() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 1.0, 0.0));
        // Invalid Tet4: only 3 node ids -> face local indexing should fail.
        mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3]));

        let err =
            collect_boundary_polygons(&mesh).expect_err("invalid face connectivity should fail");
        match err {
            Mesh3DError::Generation(msg) => {
                assert!(msg.contains("connectivity is inconsistent"));
            }
        }
    }

    #[test]
    fn centroid_of_nodes_reports_empty_and_missing_node_errors() {
        let mesh = Mesh::new();
        let empty = HashSet::<u64>::new();
        let err = centroid_of_nodes(&mesh, &empty).expect_err("empty centroid set should fail");
        match err {
            Mesh3DError::Generation(msg) => assert!(msg.contains("empty node set")),
        }

        let mut ids = HashSet::<u64>::new();
        ids.insert(42);
        let err = centroid_of_nodes(&mesh, &ids).expect_err("missing node should fail");
        match err {
            Mesh3DError::Generation(msg) => assert!(msg.contains("Node 42")),
        }
    }

    #[test]
    fn tetra_signed_volume6_has_expected_sign() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(0.0, 1.0, 0.0);
        let d_pos = Point3::new(0.0, 0.0, 1.0);
        let d_neg = Point3::new(0.0, 0.0, -1.0);

        let v_pos = tetra_signed_volume6(a, b, c, d_pos);
        let v_neg = tetra_signed_volume6(a, b, c, d_neg);

        assert!(v_pos > 0.0);
        assert!(v_neg < 0.0);
        assert!((v_pos + v_neg).abs() < 1e-12);
    }

    #[test]
    fn mesh3d_error_converts_to_framework_error() {
        let err = Mesh3DError::Generation("boom".to_string());
        let framework: MeshAlgoError = err.into();
        match framework {
            MeshAlgoError::Generation(msg) => assert_eq!(msg, "boom"),
            other => panic!("unexpected mapping: {other:?}"),
        }
    }

    // ── P0: Volume conservation ───────────────────────────────────────────────

    /// Sum of absolute tet volumes must equal the bounding volume of the input
    /// surface.  For a unit cube the enclosed volume is exactly 1.0.
    #[test]
    fn output_volume_sum_equals_cube_volume() {
        let mut surface = Mesh::new();
        surface.add_node(Node::new(1, 0.0, 0.0, 0.0));
        surface.add_node(Node::new(2, 1.0, 0.0, 0.0));
        surface.add_node(Node::new(3, 1.0, 1.0, 0.0));
        surface.add_node(Node::new(4, 0.0, 1.0, 0.0));
        surface.add_node(Node::new(5, 0.0, 0.0, 1.0));
        surface.add_node(Node::new(6, 1.0, 0.0, 1.0));
        surface.add_node(Node::new(7, 1.0, 1.0, 1.0));
        surface.add_node(Node::new(8, 0.0, 1.0, 1.0));

        surface.add_element(Element::new(1, ElementType::Quad4, vec![1, 2, 3, 4]));
        surface.add_element(Element::new(2, ElementType::Quad4, vec![5, 6, 7, 8]));
        surface.add_element(Element::new(3, ElementType::Quad4, vec![1, 2, 6, 5]));
        surface.add_element(Element::new(4, ElementType::Quad4, vec![2, 3, 7, 6]));
        surface.add_element(Element::new(5, ElementType::Quad4, vec![3, 4, 8, 7]));
        surface.add_element(Element::new(6, ElementType::Quad4, vec![4, 1, 5, 8]));

        let vol_mesh = tetrahedralize_closed_surface(&surface).expect("tetrahedralization failed");

        let total_vol: f64 = vol_mesh
            .elements
            .iter()
            .filter(|e| e.etype == ElementType::Tetrahedron4 && e.node_ids.len() == 4)
            .map(|e| {
                let get = |id: u64| -> Point3 {
                    let n = &vol_mesh.nodes[&id];
                    Point3::new(n.position.x, n.position.y, n.position.z)
                };
                let a = get(e.node_ids[0]);
                let b = get(e.node_ids[1]);
                let c = get(e.node_ids[2]);
                let d = get(e.node_ids[3]);
                // |V| = |det[b-a, c-a, d-a]| / 6
                (tetra_signed_volume6(a, b, c, d).abs()) / 6.0
            })
            .sum();

        assert!(
            (total_vol - 1.0).abs() < 1e-9,
            "total volume should be 1.0, got {total_vol}"
        );
    }

    // ── P0: Topology correctness ──────────────────────────────────────────────

    /// Every node in the output mesh must belong to at least one tetrahedron.
    #[test]
    fn every_node_used_in_at_least_one_tet() {
        let mut surface = Mesh::new();
        surface.add_node(Node::new(1, 0.0, 0.0, 0.0));
        surface.add_node(Node::new(2, 1.0, 0.0, 0.0));
        surface.add_node(Node::new(3, 1.0, 1.0, 0.0));
        surface.add_node(Node::new(4, 0.0, 1.0, 0.0));
        surface.add_node(Node::new(5, 0.0, 0.0, 1.0));
        surface.add_node(Node::new(6, 1.0, 0.0, 1.0));
        surface.add_node(Node::new(7, 1.0, 1.0, 1.0));
        surface.add_node(Node::new(8, 0.0, 1.0, 1.0));
        surface.add_element(Element::new(1, ElementType::Quad4, vec![1, 2, 3, 4]));
        surface.add_element(Element::new(2, ElementType::Quad4, vec![5, 6, 7, 8]));
        surface.add_element(Element::new(3, ElementType::Quad4, vec![1, 2, 6, 5]));
        surface.add_element(Element::new(4, ElementType::Quad4, vec![2, 3, 7, 6]));
        surface.add_element(Element::new(5, ElementType::Quad4, vec![3, 4, 8, 7]));
        surface.add_element(Element::new(6, ElementType::Quad4, vec![4, 1, 5, 8]));

        let vol_mesh = tetrahedralize_closed_surface(&surface).expect("tetrahedralization failed");

        let used: std::collections::HashSet<u64> = vol_mesh
            .elements
            .iter()
            .flat_map(|e| e.node_ids.iter().copied())
            .collect();

        for id in vol_mesh.nodes.keys() {
            assert!(used.contains(id), "node {id} is unreferenced in any element");
        }
    }

    /// All tetrahedra must have non-zero volume (no degenerate elements).
    #[test]
    fn no_degenerate_tets_in_cube_output() {
        let mut surface = Mesh::new();
        surface.add_node(Node::new(1, 0.0, 0.0, 0.0));
        surface.add_node(Node::new(2, 1.0, 0.0, 0.0));
        surface.add_node(Node::new(3, 1.0, 1.0, 0.0));
        surface.add_node(Node::new(4, 0.0, 1.0, 0.0));
        surface.add_node(Node::new(5, 0.0, 0.0, 1.0));
        surface.add_node(Node::new(6, 1.0, 0.0, 1.0));
        surface.add_node(Node::new(7, 1.0, 1.0, 1.0));
        surface.add_node(Node::new(8, 0.0, 1.0, 1.0));
        surface.add_element(Element::new(1, ElementType::Quad4, vec![1, 2, 3, 4]));
        surface.add_element(Element::new(2, ElementType::Quad4, vec![5, 6, 7, 8]));
        surface.add_element(Element::new(3, ElementType::Quad4, vec![1, 2, 6, 5]));
        surface.add_element(Element::new(4, ElementType::Quad4, vec![2, 3, 7, 6]));
        surface.add_element(Element::new(5, ElementType::Quad4, vec![3, 4, 8, 7]));
        surface.add_element(Element::new(6, ElementType::Quad4, vec![4, 1, 5, 8]));

        let vol_mesh = tetrahedralize_closed_surface(&surface).expect("tetrahedralization failed");

        for e in vol_mesh
            .elements
            .iter()
            .filter(|e| e.etype == ElementType::Tetrahedron4 && e.node_ids.len() == 4)
        {
            let get = |id: u64| -> Point3 {
                let n = &vol_mesh.nodes[&id];
                Point3::new(n.position.x, n.position.y, n.position.z)
            };
            let vol6 = tetra_signed_volume6(
                get(e.node_ids[0]),
                get(e.node_ids[1]),
                get(e.node_ids[2]),
                get(e.node_ids[3]),
            )
            .abs();
            assert!(
                vol6 > 1e-14,
                "element {} is degenerate (vol6 = {vol6})",
                e.id
            );
        }
    }

    /// The number of output 3D elements must equal 12 for a cube surface
    /// (6 quad faces × 2 triangles each = 12 boundary triangles).
    #[test]
    fn element_count_is_correct_for_cube() {
        let mut surface = Mesh::new();
        surface.add_node(Node::new(1, 0.0, 0.0, 0.0));
        surface.add_node(Node::new(2, 1.0, 0.0, 0.0));
        surface.add_node(Node::new(3, 1.0, 1.0, 0.0));
        surface.add_node(Node::new(4, 0.0, 1.0, 0.0));
        surface.add_node(Node::new(5, 0.0, 0.0, 1.0));
        surface.add_node(Node::new(6, 1.0, 0.0, 1.0));
        surface.add_node(Node::new(7, 1.0, 1.0, 1.0));
        surface.add_node(Node::new(8, 0.0, 1.0, 1.0));
        surface.add_element(Element::new(1, ElementType::Quad4, vec![1, 2, 3, 4]));
        surface.add_element(Element::new(2, ElementType::Quad4, vec![5, 6, 7, 8]));
        surface.add_element(Element::new(3, ElementType::Quad4, vec![1, 2, 6, 5]));
        surface.add_element(Element::new(4, ElementType::Quad4, vec![2, 3, 7, 6]));
        surface.add_element(Element::new(5, ElementType::Quad4, vec![3, 4, 8, 7]));
        surface.add_element(Element::new(6, ElementType::Quad4, vec![4, 1, 5, 8]));

        let vol_mesh = tetrahedralize_closed_surface(&surface).expect("tetrahedralization failed");
        assert_eq!(
            vol_mesh.elements_by_dimension(3).len(),
            12,
            "cube should produce exactly 12 tetrahedra"
        );
    }

    // ── P0: centroid_of_nodes direct tests ────────────────────────────────────

    #[test]
    fn centroid_of_nodes_computes_correct_value() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 2.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 2.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 0.0, 2.0));

        let ids: std::collections::HashSet<u64> = [1, 2, 3, 4].into();
        let c = centroid_of_nodes(&mesh, &ids).expect("centroid should succeed");
        assert!((c.x - 0.5).abs() < 1e-12);
        assert!((c.y - 0.5).abs() < 1e-12);
        assert!((c.z - 0.5).abs() < 1e-12);
    }

    // ── P0: tetra_signed_volume6 additional tests ─────────────────────────────

    #[test]
    fn tetra_signed_volume6_unit_tet_value() {
        // Unit tet: a=(0,0,0), b=(1,0,0), c=(0,1,0), d=(0,0,1)
        // Volume = 1/6, so vol6 = 1.0.
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(0.0, 1.0, 0.0);
        let d = Point3::new(0.0, 0.0, 1.0);
        let v = tetra_signed_volume6(a, b, c, d);
        assert!((v.abs() - 1.0).abs() < 1e-12, "unit tet vol6 should be 1.0, got {v}");
    }

    #[test]
    fn step_mesh_can_be_saved_as_gmsh_v2_and_v4_after_3d_meshing() {
        let step_path = step_file_path("my_cube.step");
        let step_bytes = std::fs::read(&step_path)
            .unwrap_or_else(|e| panic!("failed to read {}: {}", step_path.display(), e));

        let surface_mesh = load_step_from_bytes(&step_bytes).expect("STEP should parse");
        let volume_mesh = tetrahedralize_closed_surface(&surface_mesh)
            .expect("3D meshing should succeed for cube surface");

        assert!(volume_mesh.node_count() > 0);
        assert!(volume_mesh.elements_by_dimension(3).len() > 0);

        // Save to Gmsh v4 and load back.
        let mut v4_bytes = Vec::new();
        write_msh_v4(&mut v4_bytes, &volume_mesh).expect("MSH v4 write should succeed");
        let v4_loaded = load_msh_from_bytes(&v4_bytes).expect("MSH v4 readback should succeed");
        assert_eq!(v4_loaded.node_count(), volume_mesh.node_count());
        assert_eq!(v4_loaded.element_count(), volume_mesh.element_count());

        // Save to Gmsh v2 and load back.
        let mut v2_bytes = Vec::new();
        write_msh_v2(&mut v2_bytes, &volume_mesh).expect("MSH v2 write should succeed");
        let v2_loaded = load_msh_from_bytes(&v2_bytes).expect("MSH v2 readback should succeed");
        assert_eq!(v2_loaded.node_count(), volume_mesh.node_count());
        assert_eq!(v2_loaded.element_count(), volume_mesh.element_count());
    }
}
