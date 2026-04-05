// ---------------------------------------------------------------------------
// Picking — CPU ray-triangle intersection for field value probing
// ---------------------------------------------------------------------------

use crate::mesh_data::FieldMesh;

/// Result of a field probe pick operation.
#[derive(Debug, Clone)]
pub struct PickResult {
    pub position: [f32; 3],
    pub field_value: f32,
    pub triangle_idx: usize,
}

/// Pick the nearest triangle in a FieldMesh by ray-triangle intersection.
/// Returns the hit position and interpolated field value.
pub fn pick_field(
    mesh: &FieldMesh,
    ray_origin: [f32; 3],
    ray_dir: [f32; 3],
) -> Option<PickResult> {
    let mut nearest: Option<(f32, PickResult)> = None;

    for (tri_idx, tri) in mesh.indices.chunks(3).enumerate() {
        if tri.len() < 3 {
            continue;
        }
        let v0 = &mesh.vertices[tri[0] as usize];
        let v1 = &mesh.vertices[tri[1] as usize];
        let v2 = &mesh.vertices[tri[2] as usize];

        if let Some((t, u, v)) = ray_triangle_intersect(
            ray_origin,
            ray_dir,
            v0.position,
            v1.position,
            v2.position,
        ) {
            if t > 0.0 {
                let should_update = match &nearest {
                    Some((best_t, _)) => t < *best_t,
                    None => true,
                };
                if should_update {
                    let w = 1.0 - u - v;
                    let pos = [
                        w * v0.position[0] + u * v1.position[0] + v * v2.position[0],
                        w * v0.position[1] + u * v1.position[1] + v * v2.position[1],
                        w * v0.position[2] + u * v1.position[2] + v * v2.position[2],
                    ];
                    let field_value =
                        w * v0.field_value + u * v1.field_value + v * v2.field_value;

                    nearest = Some((
                        t,
                        PickResult {
                            position: pos,
                            field_value,
                            triangle_idx: tri_idx,
                        },
                    ));
                }
            }
        }
    }

    nearest.map(|(_, result)| result)
}

/// Möller–Trumbore ray-triangle intersection.
/// Returns (t, u, v) where t is the distance along the ray and (u, v) are
/// barycentric coordinates.
fn ray_triangle_intersect(
    origin: [f32; 3],
    dir: [f32; 3],
    v0: [f32; 3],
    v1: [f32; 3],
    v2: [f32; 3],
) -> Option<(f32, f32, f32)> {
    let edge1 = sub3(v1, v0);
    let edge2 = sub3(v2, v0);
    let h = cross3(dir, edge2);
    let a = dot3(edge1, h);

    if a.abs() < 1e-7 {
        return None; // Parallel
    }

    let f = 1.0 / a;
    let s = sub3(origin, v0);
    let u = f * dot3(s, h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let q = cross3(s, edge1);
    let v = f * dot3(dir, q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    let t = f * dot3(edge2, q);
    Some((t, u, v))
}

fn sub3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_data::FieldVertex;

    #[test]
    fn ray_hits_triangle() {
        let (t, u, v) = ray_triangle_intersect(
            [0.25, 0.25, -1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        )
        .unwrap();

        assert!((t - 1.0).abs() < 1e-4);
        assert!((u - 0.25).abs() < 1e-4);
        assert!((v - 0.25).abs() < 1e-4);
    }

    #[test]
    fn ray_misses_triangle() {
        let result = ray_triangle_intersect(
            [2.0, 2.0, -1.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        assert!(result.is_none());
    }

    #[test]
    fn pick_interpolates_field() {
        let mesh = FieldMesh {
            vertices: vec![
                FieldVertex { position: [0.0, 0.0, 0.0], normal: [0.0, 0.0, 1.0], field_value: 0.0 },
                FieldVertex { position: [1.0, 0.0, 0.0], normal: [0.0, 0.0, 1.0], field_value: 1.0 },
                FieldVertex { position: [0.0, 1.0, 0.0], normal: [0.0, 0.0, 1.0], field_value: 2.0 },
            ],
            indices: vec![0, 1, 2],
            wire_indices: vec![],
            field_range: [0.0, 2.0],
            field_imag: None,
            vector_field: None,
        };

        let result = pick_field(&mesh, [0.25, 0.25, -1.0], [0.0, 0.0, 1.0]).unwrap();
        // Barycentric: u=0.25 (v1 weight), v=0.25 (v2 weight), w=0.5 (v0 weight)
        // field = 0.5*0.0 + 0.25*1.0 + 0.25*2.0 = 0.75
        assert!((result.field_value - 0.75).abs() < 0.05);
    }
}
