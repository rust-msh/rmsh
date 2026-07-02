//! P1 → P2 element promotion.
//!
//! Takes a first-order (P1) mesh and produces a second-order (P2) mesh by
//! inserting edge-midpoint nodes.  Shared edges between adjacent elements
//! produce a single shared midpoint node.
//!
//! # Algorithm
//!
//! 1. Collect all unique edges across all elements, keyed by `(min_node, max_node)`.
//! 2. For each unique edge, create a new node at the geometric midpoint.
//! 3. For each P1 element, create the corresponding P2 element by inserting
//!    edge-midpoint IDs in the Gmsh canonical order for the target type.
//!
//! # Supported promotions
//!
//! | P1 type | P2 type | New nodes per element |
//! |---|---|---|
//! | Line2 | Line3 | 1 (edge midpoint) |
//! | Triangle3 | Triangle6 | 3 (edge midpoints) |
//! | Quad4 | Quad9 | 5 (4 edges + centre) |
//! | Tetrahedron4 | Tetrahedron10 | 6 (edge midpoints) |
//! | Hexahedron8 | Hexahedron27 | 19 (12 edges + 6 faces + centre) |
//! | Prism6 | Prism18 | 12 (9 edges + 3 faces) |
//! | Pyramid5 | Pyramid14 | 9 (8 edges + 1 face centre) |

use std::collections::HashMap;

use rmsh_model::{Element, ElementType, Mesh, Node};

/// Map each P1 [`ElementType`] to its P2 counterpart.
pub fn p2_type(p1: ElementType) -> Option<ElementType> {
    match p1 {
        ElementType::Line2 => Some(ElementType::Line3),
        ElementType::Triangle3 => Some(ElementType::Triangle6),
        ElementType::Quad4 => Some(ElementType::Quad9),
        ElementType::Tetrahedron4 => Some(ElementType::Tetrahedron10),
        ElementType::Hexahedron8 => Some(ElementType::Hexahedron27),
        ElementType::Prism6 => Some(ElementType::Prism18),
        ElementType::Pyramid5 => Some(ElementType::Pyramid14),
        _ => None,
    }
}

/// Promote all mesh elements from P1 to P2.
///
/// Line, triangle, quadrilateral, tetrahedron, hexahedron, prism, and pyramid
/// elements are promoted.  Unknown element types are left unchanged.
///
/// Newly created nodes receive IDs starting from `next_node_id` and are added
/// to the mesh.  Existing element IDs are preserved (the element count does
/// not change — only node counts increase).
pub fn promote_to_p2(mesh: &mut Mesh) {
    // Phase 1: collect all unique edges + their midpoints (immutable read).
    let mut edge_midpoints: HashMap<(u64, u64), [f64; 3]> = HashMap::new();
    for element in &mesh.elements {
        let Some(_target_type) = p2_type(element.etype) else {
            continue;
        };
        let edges = element.etype.edges();
        for &[i, j] in edges {
            let a = element.node_ids[i];
            let b = element.node_ids[j];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_midpoints.entry(key).or_insert_with(|| {
                let na = &mesh.nodes[&a].position;
                let nb = &mesh.nodes[&b].position;
                [
                    (na.x + nb.x) * 0.5,
                    (na.y + nb.y) * 0.5,
                    (na.z + nb.z) * 0.5,
                ]
            });
        }
    }

    // Phase 2: insert all new midpoint nodes (mutable write).
    let mut next_node_id: u64 = mesh.nodes.keys().copied().max().unwrap_or(0) + 1;
    let mut edge_map: HashMap<(u64, u64), u64> = HashMap::new();
    for (key, pos) in &edge_midpoints {
        let nid = next_node_id;
        next_node_id += 1;
        mesh.add_node(Node::new(nid, pos[0], pos[1], pos[2]));
        edge_map.insert(*key, nid);
    }

    // Phase 3: rebuild each element (mutable write, no borrows on nodes).
    // For Hexahedron27: face centres + interior need additional node creation.
    let mut extra_nodes: Vec<(u64, u64, [f64; 3])> = Vec::new(); // (elem_idx, nid, pos)

    for element in &mut mesh.elements {
        let Some(target_type) = p2_type(element.etype) else {
            continue;
        };
        let edges = element.etype.edges();
        let mut new_nodes = element.node_ids.clone(); // corners first

        // Append edge-midpoint IDs.
        for &[i, j] in edges {
            let a = element.node_ids[i];
            let b = element.node_ids[j];
            let key = if a < b { (a, b) } else { (b, a) };
            if let Some(&mid) = edge_map.get(&key) {
                new_nodes.push(mid);
            }
        }

        // Hexahedron27: additional face-centre + interior nodes.
        if target_type == ElementType::Hexahedron27 {
            let faces = element.etype.faces();
            for face in faces {
                let nid = next_node_id;
                next_node_id += 1;
                let cx = face.iter().map(|&idx| element.node_ids[idx]).map(|id| mesh.nodes[&id].position)
                    .map(|p| p.x).sum::<f64>() / face.len() as f64;
                let cy = face.iter().map(|&idx| element.node_ids[idx]).map(|id| mesh.nodes[&id].position)
                    .map(|p| p.y).sum::<f64>() / face.len() as f64;
                let cz = face.iter().map(|&idx| element.node_ids[idx]).map(|id| mesh.nodes[&id].position)
                    .map(|p| p.z).sum::<f64>() / face.len() as f64;
                extra_nodes.push((element.id, nid, [cx, cy, cz]));
                new_nodes.push(nid);
            }
            // Interior centroid.
            let nid = next_node_id;
            next_node_id += 1;
            let (cx, cy, cz) = {
                let mut sx = 0.0_f64; let mut sy = 0.0_f64; let mut sz = 0.0_f64;
                let n = element.node_ids.len() as f64;
                for &id in &element.node_ids {
                    let p = &mesh.nodes[&id].position;
                    sx += p.x; sy += p.y; sz += p.z;
                }
                (sx / n, sy / n, sz / n)
            };
            extra_nodes.push((element.id, nid, [cx, cy, cz]));
            new_nodes.push(nid);
        }

        element.etype = target_type;
        element.node_ids = new_nodes;
    }

    // Phase 4: insert all extra (face-centre, interior) nodes.
    // We skip deduplication (no two Hex27 elements share a face centre).
    for (_elem_id, nid, pos) in extra_nodes {
        mesh.add_node(Node::new(nid, pos[0], pos[1], pos[2]));
    }
}

/// Compute the target P2 node count for a given P1 type.
/// Returns `None` for non-promotable types.
pub fn p2_node_count(p1: ElementType) -> Option<usize> {
    p2_type(p1).map(|t| t.node_count())
}

#[cfg(test)]
mod tests {
    use super::*;
use rmsh_model::{ElementType, Mesh, Node};

    fn make_line_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 2.0, 0.0, 0.0));
        mesh.add_element(Element::new(1, ElementType::Line2, vec![1, 2]));
        mesh
    }

    fn make_tri_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 2.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 2.0, 0.0));
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));
        mesh
    }

    fn make_two_tri_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 2.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 2.0, 2.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 2.0, 0.0));
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));
        mesh.add_element(Element::new(2, ElementType::Triangle3, vec![1, 3, 4]));
        mesh
    }

    fn make_tet_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 0.0, 1.0));
        mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));
        mesh
    }

    fn make_cube_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        for i in 0..8 {
            let x = if i & 1 != 0 { 1.0 } else { 0.0 };
            let y = if i & 2 != 0 { 1.0 } else { 0.0 };
            let z = if i & 4 != 0 { 1.0 } else { 0.0 };
            mesh.add_node(Node::new((i + 1) as u64, x, y, z));
        }
        mesh.add_element(Element::new(1, ElementType::Hexahedron8, vec![1, 2, 4, 3, 5, 6, 8, 7]));
        mesh
    }

    #[test]
    fn p2_type_mapping_correct() {
        assert_eq!(p2_type(ElementType::Line2), Some(ElementType::Line3));
        assert_eq!(p2_type(ElementType::Triangle3), Some(ElementType::Triangle6));
        assert_eq!(p2_type(ElementType::Tetrahedron4), Some(ElementType::Tetrahedron10));
        assert_eq!(p2_type(ElementType::Quad4), Some(ElementType::Quad9));
        assert_eq!(p2_type(ElementType::Hexahedron8), Some(ElementType::Hexahedron27));
        assert_eq!(p2_type(ElementType::Prism6), Some(ElementType::Prism18));
        assert_eq!(p2_type(ElementType::Pyramid5), Some(ElementType::Pyramid14));
        assert!(p2_type(ElementType::Point1).is_none());
    }

    #[test]
    fn promote_line_adds_midpoint() {
        let mut mesh = make_line_mesh();
        promote_to_p2(&mut mesh);
        assert_eq!(mesh.nodes.len(), 3);
        assert_eq!(mesh.elements.len(), 1);
        assert_eq!(mesh.elements[0].etype, ElementType::Line3);
        assert_eq!(mesh.elements[0].node_ids.len(), 3);
        // Midpoint should be at (1, 0, 0)
        let mid_id = mesh.elements[0].node_ids[2];
        let mid = &mesh.nodes[&mid_id];
        assert!((mid.position.x - 1.0).abs() < 1e-12);
        assert!((mid.position.y - 0.0).abs() < 1e-12);
    }

    #[test]
    fn promote_triangle_adds_three_midpoints() {
        let mut mesh = make_tri_mesh();
        promote_to_p2(&mut mesh);
        assert_eq!(mesh.nodes.len(), 6);
        assert_eq!(mesh.elements[0].etype, ElementType::Triangle6);
        assert_eq!(mesh.elements[0].node_ids.len(), 6);
    }

    #[test]
    fn promote_shared_edge_creates_one_midpoint() {
        let mut mesh = make_two_tri_mesh();
        promote_to_p2(&mut mesh);
        // shared edge (1,3) should have only one midpoint
        // total unique edges: 5 → 5 new nodes + 4 original = 9
        assert_eq!(mesh.nodes.len(), 9);
        assert_eq!(mesh.elements.len(), 2);
        for elt in &mesh.elements {
            assert_eq!(elt.etype, ElementType::Triangle6);
            assert_eq!(elt.node_ids.len(), 6);
        }
    }

    #[test]
    fn promote_tet_adds_six_midpoints() {
        let mut mesh = make_tet_mesh();
        promote_to_p2(&mut mesh);
        assert_eq!(mesh.nodes.len(), 10);
        assert_eq!(mesh.elements[0].etype, ElementType::Tetrahedron10);
        assert_eq!(mesh.elements[0].node_ids.len(), 10);
    }

    #[test]
    fn promote_hex_adds_19_nodes() {
        let mut mesh = make_cube_mesh();
        promote_to_p2(&mut mesh);
        // 8 corners + 12 edge midpoints + 6 face centres + 1 interior = 27
        assert_eq!(mesh.nodes.len(), 27);
        assert_eq!(mesh.elements[0].etype, ElementType::Hexahedron27);
        assert_eq!(mesh.elements[0].node_ids.len(), 27);
    }

    #[test]
    fn p2_node_count_returns_correct_values() {
        assert_eq!(p2_node_count(ElementType::Triangle3), Some(6));
        assert_eq!(p2_node_count(ElementType::Tetrahedron4), Some(10));
        assert_eq!(p2_node_count(ElementType::Hexahedron8), Some(27));
        assert_eq!(p2_node_count(ElementType::Point1), None);
    }

    #[test]
    fn promote_empty_mesh_does_not_panic() {
        let mut mesh = Mesh::new();
        promote_to_p2(&mut mesh);
        assert!(mesh.nodes.is_empty());
        assert!(mesh.elements.is_empty());
    }

    #[test]
    fn promote_tet_midpoint_coordinates() {
        let mut mesh = make_tet_mesh();
        promote_to_p2(&mut mesh);
        let elt = &mesh.elements[0];
        // Edges of Tet4 in order: (0,1), (1,2), (2,0), (0,3), (1,3), (2,3)
        // Nodes 1-4 are corners, 5-10 are midpoints
        // Edge 0-1 (nodes 1,2) → midpoint at (0.5, 0.0, 0.0) should be node 5
        let mid = &mesh.nodes[&elt.node_ids[4]];
        assert!((mid.position.x - 0.5).abs() < 1e-12);
        assert!((mid.position.y - 0.0).abs() < 1e-12);
        assert!((mid.position.z - 0.0).abs() < 1e-12);
    }
}
