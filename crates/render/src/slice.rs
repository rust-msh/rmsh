use emstudio_domain::emsfld_loader::FieldBlock;
use emstudio_domain::msh_loader::{self, MshMesh};

use crate::field_mapping::FieldComponent;
use crate::mesh_data::{FieldMesh, FieldVertex};

/// Axis-aligned slice plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SliceAxis {
    /// Slice at Z = value (XY plane)
    Z,
    /// Slice at Y = value (XZ plane)
    Y,
    /// Slice at X = value (YZ plane)
    X,
}

/// Generate a planar mesh at the given slice position, evaluating a field function
/// at each grid point.
///
/// The plane spans `[-extent, extent]` in both tangent directions, subdivided
/// into `resolution x resolution` quads.
pub fn generate_slice_mesh(
    axis: SliceAxis,
    value: f32,
    extent: f32,
    resolution: u32,
    field_fn: &dyn Fn(f32, f32, f32) -> f32,
) -> FieldMesh {
    let n = resolution + 1;
    let step = 2.0 * extent / resolution as f32;
    let mut vertices = Vec::with_capacity((n * n) as usize);
    let mut field_min = f32::MAX;
    let mut field_max = f32::MIN;

    for iv in 0..n {
        for iu in 0..n {
            let u = -extent + iu as f32 * step;
            let v = -extent + iv as f32 * step;

            let (pos, normal) = match axis {
                SliceAxis::Z => ([u, v, value], [0.0, 0.0, 1.0]),
                SliceAxis::Y => ([u, value, v], [0.0, 1.0, 0.0]),
                SliceAxis::X => ([value, u, v], [1.0, 0.0, 0.0]),
            };

            let fv = field_fn(pos[0], pos[1], pos[2]);
            field_min = field_min.min(fv);
            field_max = field_max.max(fv);

            vertices.push(FieldVertex {
                position: pos,
                normal,
                field_value: fv,
            });
        }
    }

    let mut indices = Vec::new();
    for iv in 0..(n - 1) {
        for iu in 0..(n - 1) {
            let i00 = iv * n + iu;
            let i01 = i00 + 1;
            let i10 = (iv + 1) * n + iu;
            let i11 = i10 + 1;
            indices.extend_from_slice(&[i00, i10, i01]);
            indices.extend_from_slice(&[i01, i10, i11]);
        }
    }

    // Border wireframe only
    let mut wire_indices = Vec::new();
    for iu in 0..resolution {
        // bottom edge
        wire_indices.push(iu);
        wire_indices.push(iu + 1);
        // top edge
        wire_indices.push((n - 1) * n + iu);
        wire_indices.push((n - 1) * n + iu + 1);
        // left edge
        wire_indices.push(iu * n);
        wire_indices.push((iu + 1) * n);
        // right edge
        wire_indices.push(iu * n + (n - 1));
        wire_indices.push((iu + 1) * n + (n - 1));
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

/// Synthetic field function: a wave-like pattern in 3D space.
pub fn synthetic_volume_field(x: f32, y: f32, z: f32) -> f32 {
    let r = (x * x + y * y + z * z).sqrt();
    if r < 0.01 {
        return 1.0;
    }
    (4.0 * r).sin() / r
}

/// Generate a slice mesh by intersecting a plane with tetrahedra from real data.
///
/// For each tet that the plane intersects, compute the intersection polygon
/// (3 or 4 vertices), triangulate it, and interpolate field values.
pub fn generate_slice_mesh_from_tets(
    mesh: &MshMesh,
    field: &FieldBlock,
    component: FieldComponent,
    axis: SliceAxis,
    value: f32,
) -> FieldMesh {
    let mut vertices: Vec<FieldVertex> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    let mut field_min = f32::MAX;
    let mut field_max = f32::MIN;

    // Pre-compute field values per node
    let mut node_field: std::collections::HashMap<u64, f64> = std::collections::HashMap::new();
    for node in &mesh.nodes {
        let node_idx = (node.tag - 1) as usize;
        if node_idx < field.num_nodes {
            let val = crate::isosurface::extract_isosurface_scalar(field, node_idx, component);
            node_field.insert(node.tag, val);
        }
    }

    let axis_idx = match axis {
        SliceAxis::X => 0,
        SliceAxis::Y => 1,
        SliceAxis::Z => 2,
    };
    let normal: [f32; 3] = match axis {
        SliceAxis::X => [1.0, 0.0, 0.0],
        SliceAxis::Y => [0.0, 1.0, 0.0],
        SliceAxis::Z => [0.0, 0.0, 1.0],
    };

    for elem in &mesh.elements {
        if elem.element_type != msh_loader::element_types::TET4 {
            continue;
        }
        if elem.node_tags.len() < 4 {
            continue;
        }

        let tags = &elem.node_tags;
        let positions: Vec<[f64; 3]> = tags
            .iter()
            .take(4)
            .map(|&t| mesh.node_position(t).unwrap_or([0.0; 3]))
            .collect();
        let field_vals: Vec<f64> = tags
            .iter()
            .take(4)
            .map(|&t| node_field.get(&t).copied().unwrap_or(0.0))
            .collect();

        // Classify vertices: above or below plane
        let plane_val = value as f64;
        let signs: Vec<bool> = positions.iter().map(|p| p[axis_idx] >= plane_val).collect();

        let above_count = signs.iter().filter(|&&s| s).count();
        if above_count == 0 || above_count == 4 {
            continue; // No intersection
        }

        // Find intersection points along edges
        let mut iso_points: Vec<([f32; 3], f32)> = Vec::new(); // (position, field_value)

        let tet_edges: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

        for &(va, vb) in &tet_edges {
            if signs[va] == signs[vb] {
                continue; // Both on same side
            }
            let pa = positions[va][axis_idx];
            let pb = positions[vb][axis_idx];
            let t = if (pb - pa).abs() > 1e-12 {
                (plane_val - pa) / (pb - pa)
            } else {
                0.5
            };
            let t = t.clamp(0.0, 1.0);

            let pos = [
                (positions[va][0] + t * (positions[vb][0] - positions[va][0])) as f32,
                (positions[va][1] + t * (positions[vb][1] - positions[va][1])) as f32,
                (positions[va][2] + t * (positions[vb][2] - positions[va][2])) as f32,
            ];
            let fv = (field_vals[va] + t * (field_vals[vb] - field_vals[va])) as f32;
            iso_points.push((pos, fv));
        }

        if iso_points.len() < 3 {
            continue;
        }

        // Triangulate the intersection polygon (3 or 4 points)
        let base_idx = vertices.len() as u32;
        for &(pos, fv) in &iso_points {
            field_min = field_min.min(fv);
            field_max = field_max.max(fv);
            vertices.push(FieldVertex {
                position: pos,
                normal,
                field_value: fv,
            });
        }

        // Fan triangulation from first vertex
        for i in 1..iso_points.len() - 1 {
            indices.push(base_idx);
            indices.push(base_idx + i as u32);
            indices.push(base_idx + i as u32 + 1);
        }
    }

    if field_min > field_max {
        field_min = 0.0;
        field_max = 1.0;
    }

    // Build border wireframe
    let mut wire_indices = Vec::new();
    let mut edge_set: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
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
        field_range: [field_min, field_max],
        field_imag: None,
        vector_field: None,
    }
}
