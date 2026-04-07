use rcad_render::Camera;
use wgpu::util::DeviceExt;

#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
struct GizmoVertex {
    position: [f32; 3],
    color: [f32; 3],
}

/// Axis gizmo -- renders XYZ coordinate axes with labels in a corner of the viewport.
pub struct GizmoRenderer {
    pipeline: wgpu::RenderPipeline,
    vertex_buffer: wgpu::Buffer,
    uniform_buffer: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    vertex_count: u32,
}

impl GizmoRenderer {
    pub fn new(device: &wgpu::Device, surface_format: wgpu::TextureFormat) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Gizmo Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("gizmo.wgsl").into()),
        });

        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Gizmo Uniform Buffer"),
            contents: bytemuck::cast_slice(&glam::Mat4::IDENTITY.to_cols_array_2d()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("Gizmo Bind Group Layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Gizmo Bind Group"),
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Gizmo Pipeline Layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Gizmo Pipeline"),
            layout: Some(&pipeline_layout),
            cache: None,
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_main"),
                compilation_options: Default::default(),
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<GizmoVertex>() as wgpu::BufferAddress,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3],
                }],
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_main"),
                compilation_options: Default::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format: surface_format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::LineList,
                ..Default::default()
            },
            depth_stencil: None,
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
        });

        let vertices = Self::build_vertices();
        let vertex_count = vertices.len() as u32;
        let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Gizmo Vertex Buffer"),
            contents: bytemuck::cast_slice(&vertices),
            usage: wgpu::BufferUsages::VERTEX,
        });

        Self {
            pipeline,
            vertex_buffer,
            uniform_buffer,
            bind_group,
            vertex_count,
        }
    }

    fn build_vertices() -> Vec<GizmoVertex> {
        let axis_len = 1.0f32;
        let gap = 0.08f32;
        let ls = 0.16f32;
        let hs = ls * 0.5;

        let xc = [1.0f32, 0.15, 0.15];
        let yc = [0.2f32, 1.0, 0.2];
        let zc = [0.2f32, 0.65, 1.0];

        let mut v = Vec::new();
        let mut seg = |a: [f32; 3], b: [f32; 3], c: [f32; 3]| {
            v.push(GizmoVertex { position: a, color: c });
            v.push(GizmoVertex { position: b, color: c });
        };

        // Axes
        seg([0.0, 0.0, 0.0], [axis_len, 0.0, 0.0], xc);
        seg([0.0, 0.0, 0.0], [0.0, axis_len, 0.0], yc);
        seg([0.0, 0.0, 0.0], [0.0, 0.0, axis_len], zc);

        // X label at +X
        let x0 = axis_len + gap;
        seg([x0, -hs, -hs], [x0, hs, hs], xc);
        seg([x0, -hs, hs], [x0, hs, -hs], xc);

        // Y label at +Y
        let y0 = axis_len + gap;
        seg([-hs, y0 + hs * 0.6, 0.0], [0.0, y0, 0.0], yc);
        seg([hs, y0 + hs * 0.6, 0.0], [0.0, y0, 0.0], yc);
        seg([0.0, y0, 0.0], [0.0, y0 - hs, 0.0], yc);

        // Z label at +Z
        let z0 = axis_len + gap;
        seg([-hs, hs, z0], [hs, hs, z0], zc);
        seg([hs, hs, z0], [-hs, -hs, z0], zc);
        seg([-hs, -hs, z0], [hs, -hs, z0], zc);

        v
    }

    /// Update the gizmo's view-projection matrix from the camera.
    pub fn update(
        &self,
        queue: &wgpu::Queue,
        camera: &Camera,
        viewport_width: u32,
        viewport_height: u32,
    ) {
        let gizmo_size_px = 80.0f32;
        let w = viewport_width.max(1) as f32;
        let h = viewport_height.max(1) as f32;

        let eye_dir = glam::Vec3::new(
            camera.rot_y.cos() * camera.rot_x.cos(),
            camera.rot_x.sin(),
            camera.rot_y.sin() * camera.rot_x.cos(),
        )
        .normalize_or_zero();
        let eye = eye_dir * 3.0;
        let view = glam::Mat4::look_at_rh(eye, glam::Vec3::ZERO, glam::Vec3::Y);

        let half = 1.5f32;
        let scale_x = gizmo_size_px / w;
        let scale_y = gizmo_size_px / h;
        let offset_x = 1.0 - scale_x;
        let offset_y = -1.0 + scale_y;

        let proj = glam::Mat4::orthographic_rh(-half, half, -half, half, 0.1, 100.0);
        let scale_mat = glam::Mat4::from_scale(glam::Vec3::new(scale_x, scale_y, 1.0));
        let translate_mat =
            glam::Mat4::from_translation(glam::Vec3::new(offset_x, offset_y, 0.0));

        let vp = translate_mat * scale_mat * proj * view;
        queue.write_buffer(
            &self.uniform_buffer,
            0,
            bytemuck::cast_slice(&vp.to_cols_array_2d()),
        );
    }

    /// Draw the gizmo into an active render pass.
    pub fn draw(&self, render_pass: &mut wgpu::RenderPass<'_>) {
        render_pass.set_pipeline(&self.pipeline);
        render_pass.set_bind_group(0, &self.bind_group, &[]);
        render_pass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
        render_pass.draw(0..self.vertex_count, 0..1);
    }
}
