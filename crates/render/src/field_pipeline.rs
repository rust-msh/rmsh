use wgpu;

use crate::colormap::ColormapType;
use crate::mesh_data::{FieldMesh, FieldVertex};

const SHADER_SRC: &str = include_str!("field_shader.wgsl");

/// Uniform data uploaded to the GPU each frame.
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct FieldUniforms {
    pub mvp: [f32; 16],
    pub eye_pos: [f32; 3],
    pub _pad0: f32,
    pub light_dir: [f32; 3],
    pub _pad1: f32,
    pub field_min: f32,
    pub field_max: f32,
    pub opacity: f32,
    pub _pad2: f32,
}

/// Owns all wgpu resources for rendering a field-colored mesh.
///
/// Architecture: renders the 3D scene to an offscreen framebuffer (with depth),
/// then blits the result to egui's render pass (which has no depth attachment).
pub struct FieldPipeline {
    // Scene rendering
    scene_pipeline: wgpu::RenderPipeline,
    wire_pipeline: wgpu::RenderPipeline,
    vertex_buf: wgpu::Buffer,
    index_buf: wgpu::Buffer,
    wire_index_buf: wgpu::Buffer,
    num_indices: u32,
    num_wire_indices: u32,
    uniform_buf: wgpu::Buffer,
    scene_bind_group: wgpu::BindGroup,

    // Colormap
    colormap_texture: wgpu::Texture,
    colormap_view: wgpu::TextureView,
    colormap_sampler: wgpu::Sampler,
    scene_bind_group_layout: wgpu::BindGroupLayout,

    // Offscreen framebuffer
    color_texture: wgpu::Texture,
    color_view: wgpu::TextureView,
    depth_view: wgpu::TextureView,
    fb_size: [u32; 2],

    // Blit pass
    blit_pipeline: wgpu::RenderPipeline,
    blit_bind_group_layout: wgpu::BindGroupLayout,
    blit_bind_group: wgpu::BindGroup,
    blit_sampler: wgpu::Sampler,

    // Target format for blit output
    #[allow(dead_code)]
    target_format: wgpu::TextureFormat,
}

impl FieldPipeline {
    pub fn new(
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        target_format: wgpu::TextureFormat,
        mesh: &FieldMesh,
        colormap: ColormapType,
    ) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("field-shader"),
            source: wgpu::ShaderSource::Wgsl(SHADER_SRC.into()),
        });

        // --- Buffers ---
        let vertex_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("field-vertex-buf"),
            contents: bytemuck::cast_slice(&mesh.vertices),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        });
        let index_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("field-index-buf"),
            contents: bytemuck::cast_slice(&mesh.indices),
            usage: wgpu::BufferUsages::INDEX,
        });
        let wire_index_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("field-wire-index-buf"),
            contents: bytemuck::cast_slice(&mesh.wire_indices),
            usage: wgpu::BufferUsages::INDEX,
        });

        let uniform_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("field-uniform-buf"),
            size: std::mem::size_of::<FieldUniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        // --- Colormap texture (256 x 1, 2D) ---
        let lut = colormap.generate_lut(256);
        let (colormap_texture, colormap_view, colormap_sampler) =
            create_colormap_texture(device, queue, &lut);

        // --- Scene bind group layout ---
        let scene_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("field-scene-bgl"),
                entries: &[
                    // uniforms
                    wgpu::BindGroupLayoutEntry {
                        binding: 0,
                        visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Uniform,
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    },
                    // colormap texture
                    wgpu::BindGroupLayoutEntry {
                        binding: 1,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Texture {
                            sample_type: wgpu::TextureSampleType::Float { filterable: true },
                            view_dimension: wgpu::TextureViewDimension::D2,
                            multisampled: false,
                        },
                        count: None,
                    },
                    // colormap sampler
                    wgpu::BindGroupLayoutEntry {
                        binding: 2,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                        count: None,
                    },
                ],
            });

        let scene_bind_group = create_scene_bind_group(
            device,
            &scene_bind_group_layout,
            &uniform_buf,
            &colormap_view,
            &colormap_sampler,
        );

        // --- Scene pipeline layout ---
        let scene_pipeline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("field-scene-pl"),
                bind_group_layouts: &[&scene_bind_group_layout],
                push_constant_ranges: &[],
            });

        // Offscreen format
        let offscreen_format = wgpu::TextureFormat::Rgba8UnormSrgb;

        // --- Scene pipeline (filled triangles) ---
        let scene_pipeline =
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("field-scene-rp"),
                layout: Some(&scene_pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &shader,
                    entry_point: Some("vs_main"),
                    buffers: &[FieldVertex::buffer_layout()],
                    compilation_options: Default::default(),
                },
                fragment: Some(wgpu::FragmentState {
                    module: &shader,
                    entry_point: Some("fs_main"),
                    targets: &[Some(wgpu::ColorTargetState {
                        format: offscreen_format,
                        blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                    compilation_options: Default::default(),
                }),
                primitive: wgpu::PrimitiveState {
                    topology: wgpu::PrimitiveTopology::TriangleList,
                    cull_mode: None, // show both sides
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

        // --- Wireframe pipeline (lines) ---
        let wire_pipeline =
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("field-wire-rp"),
                layout: Some(&scene_pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &shader,
                    entry_point: Some("vs_main"),
                    buffers: &[FieldVertex::buffer_layout()],
                    compilation_options: Default::default(),
                },
                fragment: Some(wgpu::FragmentState {
                    module: &shader,
                    entry_point: Some("fs_wire"),
                    targets: &[Some(wgpu::ColorTargetState {
                        format: offscreen_format,
                        blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                    compilation_options: Default::default(),
                }),
                primitive: wgpu::PrimitiveState {
                    topology: wgpu::PrimitiveTopology::LineList,
                    ..Default::default()
                },
                depth_stencil: Some(wgpu::DepthStencilState {
                    format: wgpu::TextureFormat::Depth32Float,
                    depth_write_enabled: false,
                    depth_compare: wgpu::CompareFunction::LessEqual,
                    stencil: Default::default(),
                    bias: wgpu::DepthBiasState {
                        constant: -2,
                        slope_scale: -1.0,
                        clamp: 0.0,
                    },
                }),
                multisample: Default::default(),
                multiview: None,
                cache: None,
            });

        // --- Offscreen framebuffer (initial 64x64, resized on first frame) ---
        let fb_size = [64, 64];
        let (color_texture, color_view, depth_view) =
            create_framebuffer(device, fb_size, offscreen_format);

        // --- Blit pipeline ---
        let blit_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("field-blit-bgl"),
                entries: &[
                    wgpu::BindGroupLayoutEntry {
                        binding: 0,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Texture {
                            sample_type: wgpu::TextureSampleType::Float { filterable: true },
                            view_dimension: wgpu::TextureViewDimension::D2,
                            multisampled: false,
                        },
                        count: None,
                    },
                    wgpu::BindGroupLayoutEntry {
                        binding: 1,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                        count: None,
                    },
                ],
            });

        let blit_pipeline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("field-blit-pl"),
                bind_group_layouts: &[&blit_bind_group_layout],
                push_constant_ranges: &[],
            });

        let blit_pipeline =
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("field-blit-rp"),
                layout: Some(&blit_pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &shader,
                    entry_point: Some("vs_blit"),
                    buffers: &[],
                    compilation_options: Default::default(),
                },
                fragment: Some(wgpu::FragmentState {
                    module: &shader,
                    entry_point: Some("fs_blit"),
                    targets: &[Some(wgpu::ColorTargetState {
                        format: target_format,
                        blend: Some(wgpu::BlendState::REPLACE),
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                    compilation_options: Default::default(),
                }),
                primitive: wgpu::PrimitiveState {
                    topology: wgpu::PrimitiveTopology::TriangleList,
                    ..Default::default()
                },
                depth_stencil: None,
                multisample: Default::default(),
                multiview: None,
                cache: None,
            });

        let blit_sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("field-blit-sampler"),
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            ..Default::default()
        });

        let blit_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("field-blit-bg"),
            layout: &blit_bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(&color_view),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::Sampler(&blit_sampler),
                },
            ],
        });

        Self {
            scene_pipeline,
            wire_pipeline,
            vertex_buf,
            index_buf,
            wire_index_buf,
            num_indices: mesh.indices.len() as u32,
            num_wire_indices: mesh.wire_indices.len() as u32,
            uniform_buf,
            scene_bind_group,
            colormap_texture,
            colormap_view,
            colormap_sampler,
            scene_bind_group_layout,
            color_texture,
            color_view,
            depth_view,
            fb_size,
            blit_pipeline,
            blit_bind_group_layout,
            blit_bind_group,
            blit_sampler,
            target_format,
        }
    }

    pub fn update_uniforms(&self, queue: &wgpu::Queue, uniforms: &FieldUniforms) {
        queue.write_buffer(&self.uniform_buf, 0, bytemuck::bytes_of(uniforms));
    }

    pub fn update_colormap(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        colormap: ColormapType,
    ) {
        let lut = colormap.generate_lut(256);
        let (tex, view, sampler) = create_colormap_texture(device, queue, &lut);
        self.colormap_texture = tex;
        self.colormap_view = view;
        self.colormap_sampler = sampler;
        // Recreate scene bind group with new texture
        self.scene_bind_group = create_scene_bind_group(
            device,
            &self.scene_bind_group_layout,
            &self.uniform_buf,
            &self.colormap_view,
            &self.colormap_sampler,
        );
    }

    /// Resize the offscreen framebuffer if the viewport size changed.
    pub fn resize_if_needed(&mut self, device: &wgpu::Device, size: [u32; 2]) {
        // Round up to multiples of 16 to reduce thrashing
        let w = ((size[0].max(1) + 15) / 16) * 16;
        let h = ((size[1].max(1) + 15) / 16) * 16;
        if [w, h] == self.fb_size {
            return;
        }
        self.fb_size = [w, h];
        let (ct, cv, dv) =
            create_framebuffer(device, [w, h], wgpu::TextureFormat::Rgba8UnormSrgb);
        self.color_texture = ct;
        self.color_view = cv;
        self.depth_view = dv;

        // Recreate blit bind group with new color view
        self.blit_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("field-blit-bg"),
            layout: &self.blit_bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(&self.color_view),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::Sampler(&self.blit_sampler),
                },
            ],
        });
    }

    /// Render the 3D scene to the offscreen framebuffer.
    /// If an `ArrowPipeline` is provided, arrows are drawn after the mesh.
    pub fn render_scene(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        show_wireframe: bool,
        show_solid: bool,
        arrow_pipeline: Option<&crate::arrow_pipeline::ArrowPipeline>,
    ) {
        let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("field-scene-pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: &self.color_view,
                depth_slice: None,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color {
                        r: 0.12,
                        g: 0.12,
                        b: 0.15,
                        a: 1.0,
                    }),
                    store: wgpu::StoreOp::Store,
                },
            })],
            depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                view: &self.depth_view,
                depth_ops: Some(wgpu::Operations {
                    load: wgpu::LoadOp::Clear(1.0),
                    store: wgpu::StoreOp::Store,
                }),
                stencil_ops: None,
            }),
            occlusion_query_set: None,
            timestamp_writes: None,
        });

        // Draw filled mesh
        if show_solid {
            rpass.set_pipeline(&self.scene_pipeline);
            rpass.set_bind_group(0, &self.scene_bind_group, &[]);
            rpass.set_vertex_buffer(0, self.vertex_buf.slice(..));
            rpass.set_index_buffer(self.index_buf.slice(..), wgpu::IndexFormat::Uint32);
            rpass.draw_indexed(0..self.num_indices, 0, 0..1);
        }

        // Draw wireframe
        if show_wireframe {
            rpass.set_pipeline(&self.wire_pipeline);
            rpass.set_bind_group(0, &self.scene_bind_group, &[]);
            rpass.set_vertex_buffer(0, self.vertex_buf.slice(..));
            rpass.set_index_buffer(self.wire_index_buf.slice(..), wgpu::IndexFormat::Uint32);
            rpass.draw_indexed(0..self.num_wire_indices, 0, 0..1);
        }

        // Draw arrows (instanced)
        if let Some(arrows) = arrow_pipeline {
            arrows.draw(&mut rpass, &self.scene_bind_group);
        }
    }

    /// Blit the offscreen color texture onto egui's render pass.
    pub fn blit(&self, render_pass: &mut wgpu::RenderPass<'_>) {
        render_pass.set_pipeline(&self.blit_pipeline);
        render_pass.set_bind_group(0, &self.blit_bind_group, &[]);
        render_pass.draw(0..3, 0..1);
    }

    /// Update vertex buffer data (e.g., for phase animation).
    pub fn update_vertices(&self, queue: &wgpu::Queue, vertices: &[FieldVertex]) {
        queue.write_buffer(&self.vertex_buf, 0, bytemuck::cast_slice(vertices));
    }

    /// Replace the mesh geometry entirely (for mode switching).
    pub fn swap_mesh(&mut self, device: &wgpu::Device, mesh: &FieldMesh) {
        self.vertex_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("field-vertex-buf"),
            contents: bytemuck::cast_slice(&mesh.vertices),
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
        });
        self.index_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("field-index-buf"),
            contents: bytemuck::cast_slice(&mesh.indices),
            usage: wgpu::BufferUsages::INDEX,
        });
        self.wire_index_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("field-wire-index-buf"),
            contents: bytemuck::cast_slice(&mesh.wire_indices),
            usage: wgpu::BufferUsages::INDEX,
        });
        self.num_indices = mesh.indices.len() as u32;
        self.num_wire_indices = mesh.wire_indices.len() as u32;
    }

    /// Access the bind group layout (for creating compatible pipelines like ArrowPipeline).
    pub fn bind_group_layout(&self) -> &wgpu::BindGroupLayout {
        &self.scene_bind_group_layout
    }

    /// Access the offscreen color texture (for screenshot readback).
    pub fn color_texture(&self) -> &wgpu::Texture {
        &self.color_texture
    }

    /// Get the current framebuffer size.
    pub fn framebuffer_size(&self) -> [u32; 2] {
        self.fb_size
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn create_colormap_texture(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    lut: &[[u8; 4]],
) -> (wgpu::Texture, wgpu::TextureView, wgpu::Sampler) {
    let width = lut.len() as u32;
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("field-colormap-tex"),
        size: wgpu::Extent3d {
            width,
            height: 1,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Rgba8UnormSrgb,
        usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
        view_formats: &[],
    });

    queue.write_texture(
        wgpu::TexelCopyTextureInfo {
            texture: &texture,
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        bytemuck::cast_slice(lut),
        wgpu::TexelCopyBufferLayout {
            offset: 0,
            bytes_per_row: Some(width * 4),
            rows_per_image: None,
        },
        wgpu::Extent3d {
            width,
            height: 1,
            depth_or_array_layers: 1,
        },
    );

    let view = texture.create_view(&Default::default());
    let sampler = device.create_sampler(&wgpu::SamplerDescriptor {
        label: Some("field-colormap-sampler"),
        address_mode_u: wgpu::AddressMode::ClampToEdge,
        mag_filter: wgpu::FilterMode::Linear,
        min_filter: wgpu::FilterMode::Linear,
        ..Default::default()
    });

    (texture, view, sampler)
}

fn create_scene_bind_group(
    device: &wgpu::Device,
    layout: &wgpu::BindGroupLayout,
    uniform_buf: &wgpu::Buffer,
    colormap_view: &wgpu::TextureView,
    colormap_sampler: &wgpu::Sampler,
) -> wgpu::BindGroup {
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("field-scene-bg"),
        layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buf.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: wgpu::BindingResource::TextureView(colormap_view),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: wgpu::BindingResource::Sampler(colormap_sampler),
            },
        ],
    })
}

fn create_framebuffer(
    device: &wgpu::Device,
    size: [u32; 2],
    color_format: wgpu::TextureFormat,
) -> (wgpu::Texture, wgpu::TextureView, wgpu::TextureView) {
    let color_tex = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("field-fb-color"),
        size: wgpu::Extent3d {
            width: size[0],
            height: size[1],
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: color_format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let color_view = color_tex.create_view(&Default::default());

    let depth_tex = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("field-fb-depth"),
        size: wgpu::Extent3d {
            width: size[0],
            height: size[1],
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Depth32Float,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
        view_formats: &[],
    });
    let depth_view = depth_tex.create_view(&Default::default());

    (color_tex, color_view, depth_view)
}

use wgpu::util::DeviceExt;
