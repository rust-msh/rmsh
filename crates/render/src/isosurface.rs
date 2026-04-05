// ---------------------------------------------------------------------------
// Isosurface Extraction — Marching Tetrahedra algorithm
// ---------------------------------------------------------------------------
//
// Extracts an isosurface from a tetrahedral mesh at a given field threshold.
// Each tetrahedron has 4 vertices → 16 sign patterns → at most 2 triangles
// per tet.

use std::collections::HashMap;

use emstudio_domain::emsfld_loader::FieldBlock;
use emstudio_domain::msh_loader::{self, MshMesh};

use crate::field_mapping::{self, FieldComponent};
use crate::mesh_data::{FieldMesh, FieldVertex};

/// Edge table for Marching Tetrahedra.
/// For each of the 16 sign patterns (4 bits, one per vertex),
/// lists the edges that are intersected. Each edge is a pair of vertex
/// local indices (0-3). Up to 4 edges → 2 triangles.
/// Format: [num_triangles, e0a, e0b, e1a, e1b, e2a, e2b, ...]
/// -1 = unused.
const EDGE_TABLE: [[i8; 7]; 16] = [
    [0, -1, -1, -1, -1, -1, -1], // 0000: all outside
    [1, 0, 1, 0, 2, 0, 3],       // 0001: v0 inside
    [1, 1, 0, 1, 3, 1, 2],       // 0010: v1 inside
    [2, 0, 2, 1, 3, 0, 3],       // 0011: v0,v1 inside (quad → 2 tris)
    [1, 2, 0, 2, 1, 2, 3],       // 0100: v2 inside
    [2, 0, 1, 2, 1, 0, 3],       // 0101: v0,v2 inside
    [2, 1, 0, 2, 0, 1, 3],       // 0110: v1,v2 inside
    [1, 0, 3, 1, 3, 2, 3],       // 0111: v0,v1,v2 inside → v3 outside
    [1, 3, 0, 3, 2, 3, 1],       // 1000: v3 inside
    [2, 0, 1, 3, 2, 0, 2],       // 1001: v0,v3 inside
    [2, 1, 0, 3, 0, 1, 2],       // 1010: v1,v3 inside
    [1, 2, 0, 2, 3, 2, 1],       // 1011: v0,v1,v3 inside → v2 outside
    [2, 2, 0, 3, 0, 2, 1],       // 1100: v2,v3 inside
    [1, 1, 0, 1, 2, 1, 3],       // 1101: v0,v2,v3 inside → v1 outside
    [1, 0, 1, 0, 3, 0, 2],       // 1110: v1,v2,v3 inside → v0 outside
    [0, -1, -1, -1, -1, -1, -1], // 1111: all inside
];

/// Extract an isosurface at the given threshold from a tet mesh with field data.
pub fn extract_isosurface(
    mesh: &MshMesh,
    field: &FieldBlock,
    component: FieldComponent,
    threshold: f64,
) -> FieldMesh {
    // Pre-compute field values at each node
    let mut node_values: HashMap<u64, f64> = HashMap::new();
    for node in &mesh.nodes {
        let node_idx = (node.tag - 1) as usize;
        if node_idx < field.num_nodes {
            let val = extract_isosurface_scalar(field, node_idx, component);
            node_values.insert(node.tag, val);
        }
    }

    let mut vertices: Vec<FieldVertex> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    let mut field_min = f32::MAX;
    let mut field_max = f32::MIN;

    // Deduplicate edge intersection points
    let mut edge_verts: HashMap<(u64, u64), u32> = HashMap::new();

    for elem in &mesh.elements {
        if elem.element_type != msh_loader::element_types::TET4 {
            continue;
        }
        if elem.node_tags.len() < 4 {
            continue;
        }

        let tags = &elem.node_tags;
        let vals: [f64; 4] = [
            node_values.get(&tags[0]).copied().unwrap_or(0.0),
            node_values.get(&tags[1]).copied().unwrap_or(0.0),
            node_values.get(&tags[2]).copied().unwrap_or(0.0),
            node_values.get(&tags[3]).copied().unwrap_or(0.0),
        ];

        let positions: [[f64; 3]; 4] = [
            mesh.node_position(tags[0]).unwrap_or([0.0; 3]),
            mesh.node_position(tags[1]).unwrap_or([0.0; 3]),
            mesh.node_position(tags[2]).unwrap_or([0.0; 3]),
            mesh.node_position(tags[3]).unwrap_or([0.0; 3]),
        ];

        // Classify vertices: inside (>= threshold) = 1, outside = 0
        let mut sign = 0u8;
        for i in 0..4 {
            if vals[i] >= threshold {
                sign |= 1 << i;
            }
        }

        let entry = &EDGE_TABLE[sign as usize];
        let num_tris = entry[0] as usize;
        if num_tris == 0 {
            continue;
        }

        // Generate triangles from edge intersections
        for tri_idx in 0..num_tris {
            let base = 1 + tri_idx * 6;
            let mut tri_indices = [0u32; 3];

            for vi in 0..3 {
                let va = entry[base + vi * 2] as usize;
                let vb = entry[base + vi * 2 + 1] as usize;
                let tag_a = tags[va];
                let tag_b = tags[vb];

                let edge_key = if tag_a < tag_b {
                    (tag_a, tag_b)
                } else {
                    (tag_b, tag_a)
                };

                let idx = if let Some(&existing) = edge_verts.get(&edge_key) {
                    existing
                } else {
                    // Interpolate position and field value along edge
                    let t = if (vals[vb] - vals[va]).abs() > 1e-12 {
                        (threshold - vals[va]) / (vals[vb] - vals[va])
                    } else {
                        0.5
                    };
                    let t = t.clamp(0.0, 1.0);

                    let pos = [
                        positions[va][0] + t * (positions[vb][0] - positions[va][0]),
                        positions[va][1] + t * (positions[vb][1] - positions[va][1]),
                        positions[va][2] + t * (positions[vb][2] - positions[va][2]),
                    ];
                    let fv = (vals[va] + t * (vals[vb] - vals[va])) as f32;
                    field_min = field_min.min(fv);
                    field_max = field_max.max(fv);

                    let new_idx = vertices.len() as u32;
                    vertices.push(FieldVertex {
                        position: [pos[0] as f32, pos[1] as f32, pos[2] as f32],
                        normal: [0.0, 0.0, 0.0], // computed below
                        field_value: fv,
                    });
                    edge_verts.insert(edge_key, new_idx);
                    new_idx
                };

                tri_indices[vi] = idx;
            }

            indices.extend_from_slice(&tri_indices);
        }
    }

    // Compute vertex normals from face normals
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
        let len = (v.normal[0] * v.normal[0] + v.normal[1] * v.normal[1] + v.normal[2] * v.normal[2]).sqrt();
        if len > 1e-8 {
            v.normal[0] /= len;
            v.normal[1] /= len;
            v.normal[2] /= len;
        }
    }

    // Build wireframe
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

    if field_min > field_max {
        field_min = 0.0;
        field_max = 1.0;
    }

    FieldMesh {
        vertices,
        indices,
        wire_indices,
        field_range: [field_min, field_max],
        field_imag: None,
        vector_field: None,
    }
}

/// Extract a scalar field value for isosurface/slice operations.
pub fn extract_isosurface_scalar(field: &FieldBlock, node_idx: usize, component: FieldComponent) -> f64 {
    match component {
        FieldComponent::Magnitude => {
            if field.num_components >= 3 {
                field.vector_magnitude(node_idx).unwrap_or(0.0)
            } else {
                field.data.get(node_idx).map(|c| c.magnitude()).unwrap_or(0.0)
            }
        }
        FieldComponent::RealPart => {
            if field.num_components >= 3 {
                let base = node_idx * field.num_components;
                let rx = field.data.get(base).map(|c| c.real).unwrap_or(0.0);
                let ry = field.data.get(base + 1).map(|c| c.real).unwrap_or(0.0);
                let rz = field.data.get(base + 2).map(|c| c.real).unwrap_or(0.0);
                (rx * rx + ry * ry + rz * rz).sqrt()
            } else {
                field.data.get(node_idx).map(|c| c.real).unwrap_or(0.0)
            }
        }
        _ => {
            field.data.get(node_idx * field.num_components).map(|c| c.magnitude()).unwrap_or(0.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use emstudio_domain::emsfld_loader::ComplexValue;
    use std::io::{BufReader, Cursor};

    #[test]
    fn isosurface_single_tet() {
        // Create a single tet with field values 0, 0, 0, 1
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
        let mesh = MshMesh::read_from(&mut reader).unwrap();

        // Field: node 4 has value 1.0, others 0.0
        let field = FieldBlock {
            frequency: 1.0,
            data: vec![
                ComplexValue { real: 0.0, imag: 0.0 },
                ComplexValue { real: 0.0, imag: 0.0 },
                ComplexValue { real: 0.0, imag: 0.0 },
                ComplexValue { real: 1.0, imag: 0.0 },
            ],
            num_nodes: 4,
            num_components: 1,
        };

        let iso = extract_isosurface(&mesh, &field, FieldComponent::Magnitude, 0.5);
        // One vertex inside (node 4) → 1 triangle
        assert_eq!(iso.indices.len(), 3);
        assert_eq!(iso.vertices.len(), 3);
    }
}
