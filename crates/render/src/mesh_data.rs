use std::collections::HashSet;

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct FieldVertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
    pub field_value: f32,
}

impl FieldVertex {
    pub fn buffer_layout() -> wgpu::VertexBufferLayout<'static> {
        static ATTRS: &[wgpu::VertexAttribute] = &[
            wgpu::VertexAttribute {
                format: wgpu::VertexFormat::Float32x3,
                offset: 0,
                shader_location: 0,
            },
            wgpu::VertexAttribute {
                format: wgpu::VertexFormat::Float32x3,
                offset: 12,
                shader_location: 1,
            },
            wgpu::VertexAttribute {
                format: wgpu::VertexFormat::Float32,
                offset: 24,
                shader_location: 2,
            },
        ];
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<FieldVertex>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: ATTRS,
        }
    }
}

/// Arrow instance data for instanced rendering.
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct ArrowInstance {
    pub position: [f32; 3],
    pub direction: [f32; 3],
    pub magnitude: f32,
    pub _pad: f32,
}

impl ArrowInstance {
    pub fn buffer_layout() -> wgpu::VertexBufferLayout<'static> {
        static ATTRS: &[wgpu::VertexAttribute] = &[
            wgpu::VertexAttribute {
                format: wgpu::VertexFormat::Float32x3,
                offset: 0,
                shader_location: 3,
            },
            wgpu::VertexAttribute {
                format: wgpu::VertexFormat::Float32x3,
                offset: 12,
                shader_location: 4,
            },
            wgpu::VertexAttribute {
                format: wgpu::VertexFormat::Float32,
                offset: 24,
                shader_location: 5,
            },
        ];
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<ArrowInstance>() as u64,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: ATTRS,
        }
    }
}

#[derive(Clone)]
pub struct FieldMesh {
    pub vertices: Vec<FieldVertex>,
    pub indices: Vec<u32>,
    pub wire_indices: Vec<u32>,
    pub field_range: [f32; 2],
    /// Imaginary part of field per-vertex (for phase animation).
    pub field_imag: Option<Vec<f32>>,
    /// Per-vertex 3D vector field (for arrow visualization).
    pub vector_field: Option<Vec<[f32; 3]>>,
}

impl FieldMesh {
    /// Generate a UV sphere with synthetic complex field values.
    pub fn uv_sphere(n_lat: u32, n_lon: u32, radius: f32) -> Self {
        let mut vertices = Vec::with_capacity(((n_lat + 1) * (n_lon + 1)) as usize);
        let mut field_min = f32::MAX;
        let mut field_max = f32::MIN;
        let mut field_imag_data = Vec::new();

        for lat in 0..=n_lat {
            let theta = std::f32::consts::PI * lat as f32 / n_lat as f32;
            let sin_theta = theta.sin();
            let cos_theta = theta.cos();

            for lon in 0..=n_lon {
                let phi = 2.0 * std::f32::consts::PI * lon as f32 / n_lon as f32;
                let sin_phi = phi.sin();
                let cos_phi = phi.cos();

                let x = sin_theta * cos_phi;
                let y = cos_theta;
                let z = sin_theta * sin_phi;

                let field_real = (3.0 * phi).sin() * (2.0 * theta).cos();
                let field_im = (2.0 * phi).cos() * (3.0 * theta).sin();

                field_min = field_min.min(field_real);
                field_max = field_max.max(field_real);

                vertices.push(FieldVertex {
                    position: [x * radius, y * radius, z * radius],
                    normal: [x, y, z],
                    field_value: field_real,
                });
                field_imag_data.push(field_im);
            }
        }

        let (indices, wire_indices) = build_sphere_indices(n_lat, n_lon);

        Self {
            vertices,
            indices,
            wire_indices,
            field_range: [field_min, field_max],
            field_imag: Some(field_imag_data),
            vector_field: None,
        }
    }

    /// Generate a subdivided cube with 3D vector field for arrow/slice demos.
    pub fn cube(subdivisions: u32, half_size: f32) -> Self {
        let n = subdivisions + 1;
        let step = 2.0 * half_size / subdivisions as f32;
        let mut vertices = Vec::new();
        let mut vector_field = Vec::new();
        let mut field_min = f32::MAX;
        let mut field_max = f32::MIN;

        // Generate the 6 faces of a cube
        let faces: &[([f32; 3], [f32; 3], [f32; 3])] = &[
            // (normal, u_axis, v_axis) — u and v span the face
            ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),  // +Z
            ([0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]), // -Z
            ([1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]),  // +X
            ([-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]),  // -X
            ([0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]),  // +Y
            ([0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),  // -Y
        ];

        let mut indices = Vec::new();

        for &(normal, u_axis, v_axis) in faces {
            let base = vertices.len() as u32;
            for iv in 0..n {
                for iu in 0..n {
                    let u = -half_size + iu as f32 * step;
                    let v = -half_size + iv as f32 * step;
                    let pos = [
                        normal[0] * half_size + u_axis[0] * u + v_axis[0] * v,
                        normal[1] * half_size + u_axis[1] * u + v_axis[1] * v,
                        normal[2] * half_size + u_axis[2] * u + v_axis[2] * v,
                    ];

                    // Synthetic vector field: a rotating vortex
                    let (vx, vy, vz) = synthetic_vector_field(pos);
                    let mag = (vx * vx + vy * vy + vz * vz).sqrt();

                    field_min = field_min.min(mag);
                    field_max = field_max.max(mag);

                    vertices.push(FieldVertex {
                        position: pos,
                        normal,
                        field_value: mag,
                    });
                    vector_field.push([vx, vy, vz]);
                }
            }

            for iv in 0..(n - 1) {
                for iu in 0..(n - 1) {
                    let i00 = base + iv * n + iu;
                    let i01 = i00 + 1;
                    let i10 = i00 + n;
                    let i11 = i10 + 1;
                    indices.extend_from_slice(&[i00, i10, i01]);
                    indices.extend_from_slice(&[i01, i10, i11]);
                }
            }
        }

        let wire_indices = build_wire_indices(&indices);

        Self {
            vertices,
            indices,
            wire_indices,
            field_range: [field_min, field_max],
            field_imag: None,
            vector_field: Some(vector_field),
        }
    }

    /// Generate arrow instance data from the vector field at a subsampled rate.
    pub fn generate_arrows(&self, every_n: u32) -> Vec<ArrowInstance> {
        let vf = match &self.vector_field {
            Some(v) => v,
            None => return Vec::new(),
        };

        let mut arrows = Vec::new();
        let max_mag = self.field_range[1].max(0.001);

        for (i, (vtx, dir)) in self.vertices.iter().zip(vf.iter()).enumerate() {
            if i as u32 % every_n != 0 {
                continue;
            }
            let mag = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
            if mag < 1e-6 {
                continue;
            }
            let norm_dir = [dir[0] / mag, dir[1] / mag, dir[2] / mag];
            arrows.push(ArrowInstance {
                position: vtx.position,
                direction: norm_dir,
                magnitude: mag / max_mag,
                _pad: 0.0,
            });
        }
        arrows
    }
}

/// Synthetic rotating vortex vector field.
fn synthetic_vector_field(pos: [f32; 3]) -> (f32, f32, f32) {
    let [x, y, z] = pos;
    // Curl-like pattern: field rotates around Y axis with some Z component
    let vx = -z * 0.8 + y * 0.3;
    let vy = (x * x + z * z).sqrt().sin() * 0.5;
    let vz = x * 0.8 - y * 0.3;
    (vx, vy, vz)
}

fn build_sphere_indices(n_lat: u32, n_lon: u32) -> (Vec<u32>, Vec<u32>) {
    let row_len = n_lon + 1;
    let mut indices = Vec::new();
    for lat in 0..n_lat {
        for lon in 0..n_lon {
            let i00 = lat * row_len + lon;
            let i01 = i00 + 1;
            let i10 = (lat + 1) * row_len + lon;
            let i11 = i10 + 1;
            indices.extend_from_slice(&[i00, i10, i01]);
            indices.extend_from_slice(&[i01, i10, i11]);
        }
    }
    let wire_indices = build_wire_indices(&indices);
    (indices, wire_indices)
}

fn build_wire_indices(indices: &[u32]) -> Vec<u32> {
    let mut edge_set: HashSet<(u32, u32)> = HashSet::new();
    let mut wire = Vec::new();
    for tri in indices.chunks(3) {
        for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            if edge_set.insert(key) {
                wire.push(a);
                wire.push(b);
            }
        }
    }
    wire
}

/// Generate a unit arrow mesh (shaft + cone head) oriented along +Y.
/// Returns (vertices, indices) for a triangle list.
pub fn generate_arrow_base_mesh() -> (Vec<[f32; 3]>, Vec<u32>) {
    let segments = 8u32;
    let shaft_radius = 0.02;
    let shaft_length = 0.7;
    let head_radius = 0.06;
    let head_length = 0.3;

    let mut verts = Vec::new();
    let mut idxs = Vec::new();

    // Shaft cylinder
    for i in 0..=segments {
        let angle = 2.0 * std::f32::consts::PI * i as f32 / segments as f32;
        let cx = angle.cos() * shaft_radius;
        let cz = angle.sin() * shaft_radius;
        verts.push([cx, 0.0, cz]);              // bottom ring
        verts.push([cx, shaft_length, cz]);      // top ring
    }
    for i in 0..segments {
        let b = i * 2;
        idxs.extend_from_slice(&[b, b + 2, b + 1]);
        idxs.extend_from_slice(&[b + 1, b + 2, b + 3]);
    }

    // Cone head
    let tip_idx = verts.len() as u32;
    verts.push([0.0, shaft_length + head_length, 0.0]); // tip
    let cone_base_start = verts.len() as u32;
    for i in 0..=segments {
        let angle = 2.0 * std::f32::consts::PI * i as f32 / segments as f32;
        let cx = angle.cos() * head_radius;
        let cz = angle.sin() * head_radius;
        verts.push([cx, shaft_length, cz]);
    }
    for i in 0..segments {
        idxs.extend_from_slice(&[cone_base_start + i, cone_base_start + i + 1, tip_idx]);
    }

    (verts, idxs)
}
