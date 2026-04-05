// ---------------------------------------------------------------------------
// Surface Extraction — Tet mesh → boundary triangles → FieldMesh
// ---------------------------------------------------------------------------
//
// Extracts the outer surface of a tetrahedral mesh by finding faces that
// appear exactly once (boundary faces). Computes vertex normals and produces
// a FieldMesh suitable for the wgpu rendering pipeline.

use std::collections::HashMap;

use emstudio_domain::msh_loader::{self, MshMesh};

use crate::mesh_data::FieldMesh;
use crate::mesh_data::FieldVertex;

/// A face key: 3 sorted node tags for unique identification.
type FaceKey = (u64, u64, u64);

fn make_face_key(a: u64, b: u64, c: u64) -> FaceKey {
    let mut arr = [a, b, c];
    arr.sort();
    (arr[0], arr[1], arr[2])
}

/// Extract boundary surface from a tetrahedral mesh and convert to FieldMesh.
///
/// For tet4 elements (4 nodes), each tet has 4 triangular faces.
/// Boundary faces are those that appear in exactly one tetrahedron.
pub fn extract_surface(mesh: &MshMesh) -> FieldMesh {
    // Count face occurrences across all tetrahedra
    let mut face_count: HashMap<FaceKey, Vec<[u64; 3]>> = HashMap::new();

    for elem in &mesh.elements {
        if elem.element_type != msh_loader::element_types::TET4 {
            continue;
        }
        if elem.node_tags.len() < 4 {
            continue;
        }

        let n = &elem.node_tags;
        // 4 faces of a tetrahedron (outward-facing convention)
        let faces: [[u64; 3]; 4] = [
            [n[0], n[2], n[1]], // face 0: nodes 0,2,1
            [n[0], n[1], n[3]], // face 1: nodes 0,1,3
            [n[1], n[2], n[3]], // face 2: nodes 1,2,3
            [n[0], n[3], n[2]], // face 3: nodes 0,3,2
        ];

        for face in &faces {
            let key = make_face_key(face[0], face[1], face[2]);
            face_count.entry(key).or_default().push(*face);
        }
    }

    // Collect boundary faces (appear exactly once)
    let boundary_faces: Vec<[u64; 3]> = face_count
        .into_iter()
        .filter(|(_, faces)| faces.len() == 1)
        .map(|(_, faces)| faces[0])
        .collect();

    // Build vertex list and index mapping
    let mut node_to_vertex: HashMap<u64, u32> = HashMap::new();
    let mut vertices: Vec<FieldVertex> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();

    for face in &boundary_faces {
        let mut tri_indices = [0u32; 3];
        for (i, &node_tag) in face.iter().enumerate() {
            let vertex_idx = if let Some(&idx) = node_to_vertex.get(&node_tag) {
                idx
            } else {
                let idx = vertices.len() as u32;
                let pos = mesh
                    .node_position(node_tag)
                    .unwrap_or([0.0, 0.0, 0.0]);
                vertices.push(FieldVertex {
                    position: [pos[0] as f32, pos[1] as f32, pos[2] as f32],
                    normal: [0.0, 0.0, 0.0], // computed below
                    field_value: 0.0,
                });
                node_to_vertex.insert(node_tag, idx);
                idx
            };
            tri_indices[i] = vertex_idx;
        }
        indices.extend_from_slice(&tri_indices);
    }

    // Compute face normals and accumulate for smooth vertex normals
    for tri in indices.chunks(3) {
        if tri.len() < 3 {
            continue;
        }
        let (i0, i1, i2) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
        let p0 = vertices[i0].position;
        let p1 = vertices[i1].position;
        let p2 = vertices[i2].position;

        let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];

        let nx = e1[1] * e2[2] - e1[2] * e2[1];
        let ny = e1[2] * e2[0] - e1[0] * e2[2];
        let nz = e1[0] * e2[1] - e1[1] * e2[0];

        for &idx in tri {
            let v = &mut vertices[idx as usize];
            v.normal[0] += nx;
            v.normal[1] += ny;
            v.normal[2] += nz;
        }
    }

    // Normalize vertex normals
    for v in &mut vertices {
        let len = (v.normal[0] * v.normal[0]
            + v.normal[1] * v.normal[1]
            + v.normal[2] * v.normal[2])
            .sqrt();
        if len > 1e-8 {
            v.normal[0] /= len;
            v.normal[1] /= len;
            v.normal[2] /= len;
        }
    }

    // Build wireframe indices (unique edges)
    let mut edge_set: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    let mut wire_indices = Vec::new();
    for tri in indices.chunks(3) {
        if tri.len() < 3 {
            continue;
        }
        for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            if edge_set.insert(key) {
                wire_indices.push(a);
                wire_indices.push(b);
            }
        }
    }

    FieldMesh {
        vertices,
        indices,
        wire_indices,
        field_range: [0.0, 0.0],
        field_imag: None,
        vector_field: None,
    }
}

/// Extract only triangle surface elements (for Q3D MoM meshes) to FieldMesh.
pub fn extract_triangles(mesh: &MshMesh) -> FieldMesh {
    let mut node_to_vertex: HashMap<u64, u32> = HashMap::new();
    let mut vertices: Vec<FieldVertex> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();

    for elem in mesh.triangles() {
        if elem.node_tags.len() < 3 {
            continue;
        }
        let mut tri_indices = [0u32; 3];
        for (i, &node_tag) in elem.node_tags.iter().take(3).enumerate() {
            let vertex_idx = if let Some(&idx) = node_to_vertex.get(&node_tag) {
                idx
            } else {
                let idx = vertices.len() as u32;
                let pos = mesh.node_position(node_tag).unwrap_or([0.0, 0.0, 0.0]);
                vertices.push(FieldVertex {
                    position: [pos[0] as f32, pos[1] as f32, pos[2] as f32],
                    normal: [0.0, 0.0, 0.0],
                    field_value: 0.0,
                });
                node_to_vertex.insert(node_tag, idx);
                idx
            };
            tri_indices[i] = vertex_idx;
        }
        indices.extend_from_slice(&tri_indices);
    }

    // Compute normals (same as above)
    for tri in indices.chunks(3) {
        if tri.len() < 3 {
            continue;
        }
        let (i0, i1, i2) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
        let p0 = vertices[i0].position;
        let p1 = vertices[i1].position;
        let p2 = vertices[i2].position;
        let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let nx = e1[1] * e2[2] - e1[2] * e2[1];
        let ny = e1[2] * e2[0] - e1[0] * e2[2];
        let nz = e1[0] * e2[1] - e1[1] * e2[0];
        for &idx in tri {
            let v = &mut vertices[idx as usize];
            v.normal[0] += nx;
            v.normal[1] += ny;
            v.normal[2] += nz;
        }
    }
    for v in &mut vertices {
        let len = (v.normal[0] * v.normal[0]
            + v.normal[1] * v.normal[1]
            + v.normal[2] * v.normal[2])
            .sqrt();
        if len > 1e-8 {
            v.normal[0] /= len;
            v.normal[1] /= len;
            v.normal[2] /= len;
        }
    }

    let mut edge_set: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    let mut wire_indices = Vec::new();
    for tri in indices.chunks(3) {
        if tri.len() < 3 {
            continue;
        }
        for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            if edge_set.insert(key) {
                wire_indices.push(a);
                wire_indices.push(b);
            }
        }
    }

    FieldMesh {
        vertices,
        indices,
        wire_indices,
        field_range: [0.0, 0.0],
        field_imag: None,
        vector_field: None,
    }
}

/// Get the node tag to vertex index mapping from an extraction.
/// Useful for mapping field data that is node-indexed to vertex-indexed.
pub fn build_node_to_vertex_map(mesh: &MshMesh, field_mesh: &FieldMesh) -> HashMap<u64, u32> {
    let mut map = HashMap::new();
    // Re-extract to build the mapping (same order as extract_surface)
    let mut vertex_idx = 0u32;
    for node in &mesh.nodes {
        // Check if this node's position matches any vertex
        for (vi, v) in field_mesh.vertices.iter().enumerate() {
            let pos = [node.x as f32, node.y as f32, node.z as f32];
            if (v.position[0] - pos[0]).abs() < 1e-6
                && (v.position[1] - pos[1]).abs() < 1e-6
                && (v.position[2] - pos[2]).abs() < 1e-6
            {
                map.insert(node.tag, vi as u32);
                break;
            }
        }
        vertex_idx += 1;
    }
    let _ = vertex_idx; // suppress unused warning
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn make_test_mesh() -> MshMesh {
        // A single tetrahedron
        let msh_data = "\
$MeshFormat\n\
4.1 0 8\n\
$EndMeshFormat\n\
$Nodes\n\
1 4 1 4\n\
3 1 0 4\n\
1\n\
2\n\
3\n\
4\n\
0.0 0.0 0.0\n\
1.0 0.0 0.0\n\
0.0 1.0 0.0\n\
0.0 0.0 1.0\n\
$EndNodes\n\
$Elements\n\
1 1 1 1\n\
3 1 4 1\n\
1 1 2 3 4\n\
$EndElements\n";
        let cursor = Cursor::new(msh_data.as_bytes());
        let mut reader = BufReader::new(cursor);
        MshMesh::read_from(&mut reader).unwrap()
    }

    #[test]
    fn extract_single_tet_surface() {
        let mesh = make_test_mesh();
        let surface = extract_surface(&mesh);

        // A single tet has 4 boundary faces = 4 triangles
        assert_eq!(surface.indices.len(), 12); // 4 triangles × 3 indices
        assert_eq!(surface.vertices.len(), 4); // 4 unique nodes

        // All normals should be non-zero
        for v in &surface.vertices {
            let len = (v.normal[0] * v.normal[0]
                + v.normal[1] * v.normal[1]
                + v.normal[2] * v.normal[2])
                .sqrt();
            assert!(len > 0.5, "normal should be normalized, got len={}", len);
        }
    }

    #[test]
    fn extract_empty_mesh() {
        // Mesh with no elements
        let msh_data = "\
$MeshFormat\n\
4.1 0 8\n\
$EndMeshFormat\n\
$Nodes\n\
1 1 1 1\n\
3 1 0 1\n\
1\n\
0.0 0.0 0.0\n\
$EndNodes\n\
$Elements\n\
0 0 1 0\n\
$EndElements\n";
        let cursor = Cursor::new(msh_data.as_bytes());
        let mut reader = BufReader::new(cursor);
        let mesh = MshMesh::read_from(&mut reader).unwrap();
        let surface = extract_surface(&mesh);

        assert!(surface.vertices.is_empty());
        assert!(surface.indices.is_empty());
    }
}
