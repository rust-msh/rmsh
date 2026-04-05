use wgpu;
use wgpu::util::DeviceExt;

use crate::mesh_data::{ArrowInstance, FieldVertex, generate_arrow_base_mesh};

const SHADER_SRC: &str = include_str!("field_shader.wgsl");

/// Pipeline for instanced arrow rendering.
/// Shares the same bind group layout and uniform buffer as `FieldPipeline`.
pub struct ArrowPipeline {
    pipeline: wgpu::RenderPipeline,
    base_vertex_buf: wgpu::Buffer,
    base_index_buf: wgpu::Buffer,
    instance_buf: wgpu::Buffer,
    num_base_indices: u32,
    num_instances: u32,
    max_instances: u32,
}

impl ArrowPipeline {
    pub fn new(
        device: &wgpu::Device,
        scene_bind_group_layout: &wgpu::BindGroupLayout,
        color_format: wgpu::TextureFormat,
        max_instances: u32,
    ) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("arrow-shader"),
            source: wgpu::ShaderSource::Wgsl(SHADER_SRC.into()),
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("arrow-pl"),
            bind_group_layouts: &[scene_bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("arrow-rp"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_arrow"),
                buffers: &[FieldVertex::buffer_layout(), ArrowInstance::buffer_layout()],
                compilation_options: Default::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: color_format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: Default::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                cull_mode: None,
                ..Default::default()
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: Default::default(),
                bias: Default::default(),
            }),
            multisample: Default::default(),
            multiview: None,
            cache: None,
        });

        // Generate base arrow geometry
        let (base_verts, base_indices) = generate_arrow_base_mesh();
        // Convert to FieldVertex (normal = rough approximation, field_value = 0)
        let base_field_verts: Vec<FieldVertex> = base_verts
            .iter()
            .map(|&p| {
                let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt().max(0.001);
                FieldVertex {
                    position: p,
                    normal: [p[0] / len, 0.5, p[2] / len], // approximate
                    field_value: 0.0,
                }
            })
            .collect();

        let base_vertex_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("arrow-base-vb"),
            contents: bytemuck::cast_slice(&base_field_verts),
            usage: wgpu::BufferUsages::VERTEX,
        });
        let base_index_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("arrow-base-ib"),
            contents: bytemuck::cast_slice(&base_indices),
            usage: wgpu::BufferUsages::INDEX,
        });

        // Pre-allocate instance buffer
        let instance_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("arrow-instance-buf"),
            size: (max_instances as u64) * std::mem::size_of::<ArrowInstance>() as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        Self {
            pipeline,
            base_vertex_buf,
            base_index_buf,
            instance_buf,
            num_base_indices: base_indices.len() as u32,
            num_instances: 0,
            max_instances,
        }
    }

    /// Upload arrow instances. Returns the number actually uploaded (capped at max).
    pub fn upload_instances(&mut self, queue: &wgpu::Queue, instances: &[ArrowInstance]) -> u32 {
        let count = (instances.len() as u32).min(self.max_instances);
        if count > 0 {
            queue.write_buffer(
                &self.instance_buf,
                0,
                bytemuck::cast_slice(&instances[..count as usize]),
            );
        }
        self.num_instances = count;
        count
    }

    /// Draw the arrows into an existing render pass.
    pub fn draw<'rp>(&'rp self, render_pass: &mut wgpu::RenderPass<'rp>, bind_group: &'rp wgpu::BindGroup) {
        if self.num_instances == 0 {
            return;
        }
        render_pass.set_pipeline(&self.pipeline);
        render_pass.set_bind_group(0, bind_group, &[]);
        render_pass.set_vertex_buffer(0, self.base_vertex_buf.slice(..));
        render_pass.set_vertex_buffer(1, self.instance_buf.slice(..));
        render_pass.set_index_buffer(self.base_index_buf.slice(..), wgpu::IndexFormat::Uint32);
        render_pass.draw_indexed(0..self.num_base_indices, 0, 0..self.num_instances);
    }
}
