use egui::{Color32, Sense, Ui, Vec2};

use emstudio_domain::Project;

#[cfg(not(target_arch = "wasm32"))]
const TRIANGLE_WGSL: &str = r#"
@vertex
fn vs_main(@builtin(vertex_index) vertex_index: u32) -> @builtin(position) vec4<f32> {
    var positions = array<vec2<f32>, 3>(
        vec2<f32>(0.0, 0.7),
        vec2<f32>(-0.7, -0.6),
        vec2<f32>(0.7, -0.6)
    );

    let pos = positions[vertex_index];
    return vec4<f32>(pos, 0.0, 1.0);
}

@fragment
fn fs_main() -> @location(0) vec4<f32> {
    return vec4<f32>(0.13, 0.66, 0.83, 1.0);
}
"#;

#[derive(Debug, Clone)]
pub struct WgpuRenderConfig {
    pub use_webgpu: bool,
    pub msaa_samples: u32,
}

impl Default for WgpuRenderConfig {
    fn default() -> Self {
        Self {
            use_webgpu: cfg!(target_arch = "wasm32"),
            msaa_samples: 1,
        }
    }
}

#[derive(Debug, Clone)]
enum RuntimeStatus {
    PendingInit,
    Ready,
    Unsupported(&'static str),
    Failed(String),
}

pub struct SceneViewport {
    pub title: String,
    pub config: WgpuRenderConfig,
    runtime_status: RuntimeStatus,
    render_state: Option<OffscreenRenderer>,
    frame_counter: u64,
}

impl Default for SceneViewport {
    fn default() -> Self {
        Self {
            title: "3D View".to_string(),
            config: WgpuRenderConfig::default(),
            runtime_status: RuntimeStatus::PendingInit,
            render_state: None,
            frame_counter: 0,
        }
    }
}

impl SceneViewport {
    pub fn ui(&mut self, ui: &mut Ui, project: &Project) {
        self.ensure_ready();

        if let Some(renderer) = &mut self.render_state {
            if let Err(err) = renderer.render_frame() {
                self.runtime_status = RuntimeStatus::Failed(err);
                self.render_state = None;
            } else {
                self.frame_counter = self.frame_counter.saturating_add(1);
            }
        }

        let desired_size = Vec2::new(ui.available_width(), ui.available_height().max(180.0));
        let (rect, response) = ui.allocate_exact_size(desired_size, Sense::drag());
        let painter = ui.painter_at(rect);

        painter.rect_filled(rect, 6.0, Color32::from_rgb(16, 24, 34));
        painter.text(
            rect.left_top() + egui::vec2(14.0, 14.0),
            egui::Align2::LEFT_TOP,
            format!("{}", self.title),
            egui::TextStyle::Heading.resolve(ui.style()),
            Color32::LIGHT_BLUE,
        );

        let status_line = match &self.runtime_status {
            RuntimeStatus::PendingInit => "wgpu: pending init".to_string(),
            RuntimeStatus::Ready => format!("wgpu: ready (offscreen), frames={}", self.frame_counter),
            RuntimeStatus::Unsupported(msg) => format!("wgpu: unsupported ({msg})"),
            RuntimeStatus::Failed(msg) => format!("wgpu: failed ({msg})"),
        };

        painter.text(
            rect.left_top() + egui::vec2(14.0, 44.0),
            egui::Align2::LEFT_TOP,
            status_line,
            egui::TextStyle::Body.resolve(ui.style()),
            Color32::WHITE,
        );

        painter.text(
            rect.left_top() + egui::vec2(14.0, 64.0),
            egui::Align2::LEFT_TOP,
            format!("Project: {}", project.title),
            egui::TextStyle::Body.resolve(ui.style()),
            Color32::WHITE,
        );

        painter.text(
            rect.left_top() + egui::vec2(14.0, 84.0),
            egui::Align2::LEFT_TOP,
            format!("Objects: {}", project.model.objects.len()),
            egui::TextStyle::Body.resolve(ui.style()),
            Color32::WHITE,
        );

        ui.ctx().request_repaint();
        if response.dragged() {
            ui.ctx().request_repaint();
        }
    }

    fn ensure_ready(&mut self) {
        if self.render_state.is_some() {
            return;
        }

        match try_create_renderer() {
            Ok(renderer) => {
                self.render_state = Some(renderer);
                self.runtime_status = RuntimeStatus::Ready;
            }
            Err(CreateRendererError::Unsupported(message)) => {
                self.runtime_status = RuntimeStatus::Unsupported(message);
            }
            Err(CreateRendererError::Failed(message)) => {
                self.runtime_status = RuntimeStatus::Failed(message);
            }
        }
    }
}

struct OffscreenRenderer {
    device: wgpu::Device,
    queue: wgpu::Queue,
    render_pipeline: wgpu::RenderPipeline,
    target_texture: wgpu::Texture,
}

impl OffscreenRenderer {
    fn render_frame(&mut self) -> Result<(), String> {
        let view = self
            .target_texture
            .create_view(&wgpu::TextureViewDescriptor::default());

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("emstudio-offscreen-encoder"),
            });

        {
            let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("emstudio-offscreen-pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.06,
                            g: 0.10,
                            b: 0.16,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                occlusion_query_set: None,
                timestamp_writes: None,
            });
            rpass.set_pipeline(&self.render_pipeline);
            rpass.draw(0..3, 0..1);
        }

        self.queue.submit([encoder.finish()]);
        Ok(())
    }
}

#[allow(dead_code)]
enum CreateRendererError {
    Unsupported(&'static str),
    Failed(String),
}

#[cfg(not(target_arch = "wasm32"))]
fn try_create_renderer() -> Result<OffscreenRenderer, CreateRendererError> {
    pollster::block_on(create_renderer_async())
}

#[cfg(target_arch = "wasm32")]
fn try_create_renderer() -> Result<OffscreenRenderer, CreateRendererError> {
    Err(CreateRendererError::Unsupported(
        "sync init disabled on wasm; use host for now",
    ))
}

#[cfg(not(target_arch = "wasm32"))]
async fn create_renderer_async() -> Result<OffscreenRenderer, CreateRendererError> {
    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());

    let adapter = instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            compatible_surface: None,
            ..Default::default()
        })
        .await
        .map_err(|err| CreateRendererError::Failed(format!("adapter request failed: {err}")))?;

    let (device, queue) = adapter
        .request_device(&wgpu::DeviceDescriptor {
            label: Some("emstudio-offscreen-device"),
            required_features: wgpu::Features::empty(),
            required_limits: wgpu::Limits::default(),
            memory_hints: wgpu::MemoryHints::default(),
            trace: wgpu::Trace::default(),
        })
        .await
        .map_err(|err| CreateRendererError::Failed(format!("device request failed: {err}")))?;

    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("emstudio-triangle-shader"),
        source: wgpu::ShaderSource::Wgsl(TRIANGLE_WGSL.into()),
    });

    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("emstudio-offscreen-pipeline-layout"),
        bind_group_layouts: &[],
        push_constant_ranges: &[],
    });

    let render_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("emstudio-offscreen-pipeline"),
        layout: Some(&pipeline_layout),
        vertex: wgpu::VertexState {
            module: &shader,
            entry_point: Some("vs_main"),
            buffers: &[],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: &shader,
            entry_point: Some("fs_main"),
            targets: &[Some(wgpu::ColorTargetState {
                format: wgpu::TextureFormat::Rgba8UnormSrgb,
                blend: Some(wgpu::BlendState::REPLACE),
                write_mask: wgpu::ColorWrites::ALL,
            })],
            compilation_options: wgpu::PipelineCompilationOptions::default(),
        }),
        primitive: wgpu::PrimitiveState::default(),
        depth_stencil: None,
        multisample: wgpu::MultisampleState::default(),
        multiview: None,
        cache: None,
    });

    let target_texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("emstudio-offscreen-target"),
        size: wgpu::Extent3d {
            width: 640,
            height: 360,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Rgba8UnormSrgb,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });

    Ok(OffscreenRenderer {
        device,
        queue,
        render_pipeline,
        target_texture,
    })
}
