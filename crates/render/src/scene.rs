use egui;
use egui_wgpu;
use wgpu;

use crate::animation::PhaseAnimator;
use crate::arrow_pipeline::ArrowPipeline;
use crate::camera::{OrbitCamera, ViewPreset};
use crate::colormap::ColormapType;
use crate::far_field;
use crate::field_mapping::{self, FieldComponent};
use crate::field_pipeline::{FieldPipeline, FieldUniforms};
use crate::isosurface;
use crate::mesh_data::FieldMesh;
use crate::mesh_quality::{self, QualityMetric};
use crate::picking;
use crate::slice;
use crate::surface_extraction;

// ---------------------------------------------------------------------------
// Visualization modes
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VisMode {
    /// Surface colormap on UV sphere
    Surface,
    /// Vector arrows on cube surface
    Arrows,
    /// Slice plane through volume
    Slice,
    /// 3D far-field radiation pattern
    FarField,
    /// Phase-animated complex field
    Animation,
    /// Wireframe only (no filled mesh)
    Wireframe,
    /// Isosurface extraction (Marching Tetrahedra)
    Isosurface,
    /// Mesh quality visualization
    MeshQuality,
}

impl VisMode {
    pub const ALL: &[VisMode] = &[
        VisMode::Surface,
        VisMode::Arrows,
        VisMode::Slice,
        VisMode::FarField,
        VisMode::Animation,
        VisMode::Wireframe,
        VisMode::Isosurface,
        VisMode::MeshQuality,
    ];

    pub fn name(&self) -> &'static str {
        match self {
            Self::Surface => "Surface",
            Self::Arrows => "Arrows",
            Self::Slice => "Slice",
            Self::FarField => "Far Field",
            Self::Animation => "Animation",
            Self::Wireframe => "Wireframe",
            Self::Isosurface => "Isosurface",
            Self::MeshQuality => "Mesh Quality",
        }
    }
}

// ---------------------------------------------------------------------------
// Callback resources
// ---------------------------------------------------------------------------

pub struct FieldSceneResources {
    pub pipeline: FieldPipeline,
    pub arrow_pipeline: Option<ArrowPipeline>,
}

// ---------------------------------------------------------------------------
// Paint callback
// ---------------------------------------------------------------------------

struct FieldSceneCallback {
    uniforms: FieldUniforms,
    show_wireframe: bool,
    show_solid: bool,
    show_arrows: bool,
    viewport_size: [u32; 2],
}

impl egui_wgpu::CallbackTrait for FieldSceneCallback {
    fn prepare(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        _screen_descriptor: &egui_wgpu::ScreenDescriptor,
        encoder: &mut wgpu::CommandEncoder,
        callback_resources: &mut egui_wgpu::CallbackResources,
    ) -> Vec<wgpu::CommandBuffer> {
        let res = callback_resources.get_mut::<FieldSceneResources>().unwrap();
        res.pipeline.resize_if_needed(device, self.viewport_size);
        res.pipeline.update_uniforms(queue, &self.uniforms);
        let arrow_ref = if self.show_arrows {
            res.arrow_pipeline.as_ref()
        } else {
            None
        };
        res.pipeline
            .render_scene(encoder, self.show_wireframe, self.show_solid, arrow_ref);
        Vec::new()
    }

    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        render_pass: &mut wgpu::RenderPass<'static>,
        callback_resources: &egui_wgpu::CallbackResources,
    ) {
        let res = callback_resources.get::<FieldSceneResources>().unwrap();
        res.pipeline.blit(render_pass);
    }
}

// ---------------------------------------------------------------------------
// High-level UI state
// ---------------------------------------------------------------------------

pub struct FieldSceneState {
    pub camera: OrbitCamera,
    pub colormap: ColormapType,
    pub opacity: f32,
    pub show_wireframe: bool,
    pub field_range: [f32; 2],
    pub vis_mode: VisMode,

    // Phase animation
    pub animator: PhaseAnimator,
    last_frame_time: Option<std::time::Instant>,

    // Slice
    pub slice_z: f32,
    pub slice_axis: slice::SliceAxis,

    // Isosurface
    pub iso_threshold: f32,
    iso_dirty: bool,

    // Arrow subsample rate
    pub arrow_subsample: u32,

    // Mesh quality
    pub quality_metric: QualityMetric,

    // Probe result (from picking)
    pub probe_result: Option<picking::PickResult>,

    // Pre-generated meshes
    sphere_mesh: Option<FieldMesh>,
    cube_mesh: Option<FieldMesh>,
    far_field_mesh: Option<FieldMesh>,

    // Real data mesh (loaded from files)
    loaded_mesh: Option<FieldMesh>,
    /// Original MSH mesh data (for tet-based operations: slicing, isosurface, quality, picking).
    loaded_msh: Option<emstudio_domain::msh_loader::MshMesh>,
    /// Field data file handle for lazy block loading.
    loaded_field: Option<emstudio_domain::emsfld_loader::EmsFldFile>,
    /// Node tag to vertex index mapping for the loaded mesh.
    node_to_vertex: Option<std::collections::HashMap<u64, u32>>,
    /// Currently selected frequency index.
    pub selected_frequency: usize,
    /// Field component to visualize.
    pub field_component: FieldComponent,
    /// Real far-field data (loaded from result store).
    loaded_far_field: Option<emstudio_domain::result_store::FarFieldData>,

    // State
    render_state: Option<egui_wgpu::RenderState>,
    colormap_dirty: bool,
    mode_dirty: bool,
    slice_dirty: bool,
}

impl Default for FieldSceneState {
    fn default() -> Self {
        Self {
            camera: OrbitCamera::default(),
            colormap: ColormapType::Rainbow,
            opacity: 1.0,
            show_wireframe: false,
            field_range: [0.0, 1.0],
            vis_mode: VisMode::Surface,
            animator: PhaseAnimator::default(),
            last_frame_time: None,
            slice_z: 0.0,
            slice_axis: slice::SliceAxis::Z,
            iso_threshold: 0.5,
            iso_dirty: false,
            arrow_subsample: 3,
            quality_metric: QualityMetric::AspectRatio,
            probe_result: None,
            sphere_mesh: None,
            cube_mesh: None,
            far_field_mesh: None,
            loaded_mesh: None,
            loaded_msh: None,
            loaded_field: None,
            node_to_vertex: None,
            selected_frequency: 0,
            field_component: FieldComponent::Magnitude,
            loaded_far_field: None,
            render_state: None,
            colormap_dirty: false,
            mode_dirty: false,
            slice_dirty: false,
        }
    }
}

impl FieldSceneState {
    pub fn new() -> Self {
        Self::default()
    }

    /// Initialize GPU resources. Call once from `App::new()`.
    pub fn init_gpu(&mut self, render_state: &egui_wgpu::RenderState, mesh: &FieldMesh) {
        self.field_range = mesh.field_range;
        self.render_state = Some(render_state.clone());

        // Pre-generate all meshes
        let cube = FieldMesh::cube(10, 1.0);
        let far_field = far_field::generate_pattern_mesh(60, 120, &far_field::patch_gain);

        let pipeline = FieldPipeline::new(
            &render_state.device,
            &render_state.queue,
            render_state.target_format,
            mesh,
            self.colormap,
        );

        // Create arrow pipeline
        let arrow_pipeline = ArrowPipeline::new(
            &render_state.device,
            pipeline.bind_group_layout(),
            wgpu::TextureFormat::Rgba8UnormSrgb,
            4096,
        );

        // Upload initial arrow data from cube mesh
        let arrows = cube.generate_arrows(3);
        let mut arrow_pipeline = arrow_pipeline;
        arrow_pipeline.upload_instances(&render_state.queue, &arrows);

        render_state
            .renderer
            .write()
            .callback_resources
            .insert(FieldSceneResources {
                pipeline,
                arrow_pipeline: Some(arrow_pipeline),
            });

        // Store meshes for mode switching
        self.sphere_mesh = Some(mesh.clone());
        self.cube_mesh = Some(cube);
        self.far_field_mesh = Some(far_field);
    }

    /// Show the 3D viewport.
    pub fn show_viewport(&mut self, ui: &mut egui::Ui) {
        let desired = egui::vec2(ui.available_width(), ui.available_height().max(100.0));
        let (rect, response) =
            ui.allocate_exact_size(desired, egui::Sense::click_and_drag());

        // Mouse interaction
        if response.dragged_by(egui::PointerButton::Primary) {
            let delta = response.drag_delta();
            self.camera.rotate(delta.x, -delta.y);
        }
        if response.dragged_by(egui::PointerButton::Middle) {
            let delta = response.drag_delta();
            self.camera.pan(delta.x, delta.y);
        }
        if response.dragged_by(egui::PointerButton::Secondary) {
            let delta = response.drag_delta();
            self.camera.pan(delta.x, delta.y);
        }
        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll.abs() > 0.01 {
            self.camera.zoom(scroll);
        }

        // Handle mode switch
        if self.mode_dirty {
            self.mode_dirty = false;
            self.apply_mode_switch();
        }

        // Handle colormap change
        if self.colormap_dirty {
            self.colormap_dirty = false;
            if let Some(rs) = &self.render_state {
                let mut renderer = rs.renderer.write();
                if let Some(res) = renderer.callback_resources.get_mut::<FieldSceneResources>() {
                    res.pipeline.update_colormap(&rs.device, &rs.queue, self.colormap);
                }
            }
        }

        // Handle slice position change
        if self.slice_dirty && self.vis_mode == VisMode::Slice {
            self.slice_dirty = false;
            self.update_slice_mesh();
        }

        // Handle isosurface threshold change
        if self.iso_dirty && self.vis_mode == VisMode::Isosurface {
            self.iso_dirty = false;
            self.update_isosurface_mesh();
        }

        // Phase animation tick
        if self.vis_mode == VisMode::Animation {
            let now = std::time::Instant::now();
            if let Some(last) = self.last_frame_time {
                let dt = now.duration_since(last).as_secs_f32();
                self.animator.tick(dt);
                self.apply_phase_animation();
            }
            self.last_frame_time = Some(now);
        } else {
            self.last_frame_time = None;
        }

        // Compute uniforms
        let pixels_per_point = ui.ctx().pixels_per_point();
        let vp_w = (rect.width() * pixels_per_point) as u32;
        let vp_h = (rect.height() * pixels_per_point) as u32;
        let aspect = vp_w as f32 / vp_h.max(1) as f32;

        let mvp = self.camera.view_projection(aspect);
        let eye = self.camera.eye_position();
        let light_dir = eye.normalize();

        let uniforms = FieldUniforms {
            mvp: mvp.to_cols_array(),
            eye_pos: eye.into(),
            _pad0: 0.0,
            light_dir: light_dir.into(),
            _pad1: 0.0,
            field_min: self.field_range[0],
            field_max: self.field_range[1],
            opacity: self.opacity,
            _pad2: 0.0,
        };

        let callback = FieldSceneCallback {
            uniforms,
            show_wireframe: self.show_wireframe || self.vis_mode == VisMode::Wireframe,
            show_solid: self.vis_mode != VisMode::Wireframe,
            show_arrows: self.vis_mode == VisMode::Arrows,
            viewport_size: [vp_w.max(1), vp_h.max(1)],
        };

        ui.painter()
            .add(egui_wgpu::Callback::new_paint_callback(rect, callback));

        ui.ctx().request_repaint();
    }

    /// Show control panel.
    pub fn show_controls(&mut self, ui: &mut egui::Ui) {
        ui.heading("Field Visualization");
        ui.separator();

        // Mode selector
        ui.horizontal(|ui| {
            ui.label("Mode:");
            let prev = self.vis_mode;
            egui::ComboBox::from_id_salt("vis_mode")
                .selected_text(self.vis_mode.name())
                .show_ui(ui, |ui| {
                    for &mode in VisMode::ALL {
                        ui.selectable_value(&mut self.vis_mode, mode, mode.name());
                    }
                });
            if self.vis_mode != prev {
                self.mode_dirty = true;
            }
        });

        ui.separator();

        // Colormap
        ui.horizontal(|ui| {
            ui.label("Colormap:");
            let prev = self.colormap;
            egui::ComboBox::from_id_salt("colormap_select")
                .selected_text(self.colormap.name())
                .show_ui(ui, |ui| {
                    for &cmap in ColormapType::ALL {
                        ui.selectable_value(&mut self.colormap, cmap, cmap.name());
                    }
                });
            if self.colormap != prev {
                self.colormap_dirty = true;
            }
        });

        ui.add(egui::Slider::new(&mut self.opacity, 0.0..=1.0).text("Opacity"));
        ui.checkbox(&mut self.show_wireframe, "Show wireframe");

        // Mode-specific controls
        ui.separator();
        match self.vis_mode {
            VisMode::Animation => {
                ui.label("Phase Animation");
                if ui
                    .button(if self.animator.playing { "Pause" } else { "Play" })
                    .clicked()
                {
                    self.animator.playing = !self.animator.playing;
                }
                ui.add(
                    egui::Slider::new(&mut self.animator.phase_deg, 0.0..=360.0)
                        .text("Phase"),
                );
                ui.add(
                    egui::Slider::new(&mut self.animator.speed_deg_per_sec, 10.0..=720.0)
                        .text("Speed (deg/s)"),
                );
            }
            VisMode::Slice => {
                ui.label("Slice Plane");
                let prev_z = self.slice_z;
                ui.add(
                    egui::Slider::new(&mut self.slice_z, -1.0..=1.0)
                        .text("Position"),
                );
                let prev_axis = self.slice_axis;
                ui.horizontal(|ui| {
                    ui.label("Axis:");
                    ui.selectable_value(&mut self.slice_axis, slice::SliceAxis::X, "X");
                    ui.selectable_value(&mut self.slice_axis, slice::SliceAxis::Y, "Y");
                    ui.selectable_value(&mut self.slice_axis, slice::SliceAxis::Z, "Z");
                });
                if (self.slice_z - prev_z).abs() > 0.001 || self.slice_axis != prev_axis {
                    self.slice_dirty = true;
                }
            }
            VisMode::Isosurface => {
                ui.label("Isosurface");
                let prev = self.iso_threshold;
                ui.add(
                    egui::Slider::new(&mut self.iso_threshold, self.field_range[0]..=self.field_range[1])
                        .text("Threshold"),
                );
                if (self.iso_threshold - prev).abs() > 0.001 {
                    self.iso_dirty = true;
                }
            }
            VisMode::MeshQuality => {
                ui.label("Mesh Quality");
                let prev = self.quality_metric;
                egui::ComboBox::from_id_salt("quality_metric")
                    .selected_text(self.quality_metric.name())
                    .show_ui(ui, |ui| {
                        for &metric in QualityMetric::ALL {
                            ui.selectable_value(&mut self.quality_metric, metric, metric.name());
                        }
                    });
                if self.quality_metric != prev {
                    self.mode_dirty = true;
                }
            }
            VisMode::FarField => {
                ui.label("3D Radiation Pattern");
                if self.loaded_far_field.is_some() {
                    ui.label("Showing loaded far-field data");
                } else {
                    ui.label("Patch antenna gain pattern (demo)");
                }
            }
            VisMode::Arrows => {
                ui.label("Vector Field Arrows");
                let prev = self.arrow_subsample;
                ui.add(
                    egui::Slider::new(&mut self.arrow_subsample, 1..=20)
                        .text("Subsample"),
                );
                if self.arrow_subsample != prev {
                    self.mode_dirty = true;
                }
            }
            VisMode::Surface => {
                if self.has_real_data() {
                    ui.label("Real field data on mesh surface");
                } else {
                    ui.label("Spherical harmonic field on sphere");
                }
            }
            VisMode::Wireframe => {
                ui.label("Wireframe mesh view");
            }
        }

        // Real data controls (frequency/component selectors)
        if self.has_real_data() && matches!(self.vis_mode, VisMode::Surface | VisMode::Arrows | VisMode::Slice | VisMode::Animation | VisMode::Isosurface) {
            ui.separator();
            ui.label("Field Data");

            // Frequency selector
            if self.num_frequencies() > 0 {
                let freqs = self.frequencies_hz();
                ui.horizontal(|ui| {
                    ui.label("Frequency:");
                    let prev = self.selected_frequency;
                    egui::ComboBox::from_id_salt("freq_select")
                        .selected_text(if self.selected_frequency < freqs.len() {
                            format!("{:.3} GHz", freqs[self.selected_frequency] / 1e9)
                        } else {
                            "N/A".to_string()
                        })
                        .show_ui(ui, |ui| {
                            for (i, &f) in freqs.iter().enumerate() {
                                ui.selectable_value(
                                    &mut self.selected_frequency,
                                    i,
                                    format!("{:.3} GHz", f / 1e9),
                                );
                            }
                        });
                    if self.selected_frequency != prev {
                        self.apply_field_to_loaded_mesh();
                    }
                });
            }

            // Component selector
            ui.horizontal(|ui| {
                ui.label("Component:");
                let prev = self.field_component;
                egui::ComboBox::from_id_salt("field_comp")
                    .selected_text(format!("{:?}", self.field_component))
                    .show_ui(ui, |ui| {
                        for &comp in &[
                            FieldComponent::Magnitude,
                            FieldComponent::RealPart,
                            FieldComponent::ImagPart,
                            FieldComponent::Phase,
                            FieldComponent::ComponentX,
                            FieldComponent::ComponentY,
                            FieldComponent::ComponentZ,
                        ] {
                            ui.selectable_value(&mut self.field_component, comp, format!("{:?}", comp));
                        }
                    });
                if self.field_component != prev {
                    self.apply_field_to_loaded_mesh();
                }
            });
        }

        // Probe result display
        if let Some(ref probe) = self.probe_result {
            ui.separator();
            ui.label(format!(
                "Probe: ({:.3}, {:.3}, {:.3})",
                probe.position[0], probe.position[1], probe.position[2]
            ));
            ui.label(format!("Value: {:.6}", probe.field_value));
        }

        // View presets
        ui.separator();
        ui.label("View Presets");
        ui.horizontal_wrapped(|ui| {
            for &preset in &[
                ViewPreset::Front,
                ViewPreset::Back,
                ViewPreset::Left,
                ViewPreset::Right,
                ViewPreset::Top,
                ViewPreset::Iso,
            ] {
                if ui.button(format!("{:?}", preset)).clicked() {
                    self.camera.set_preset(preset);
                }
            }
        });

        ui.separator();
        ui.label(format!(
            "Distance: {:.2}  Az: {:.1}  El: {:.1}",
            self.camera.distance,
            self.camera.azimuth.to_degrees(),
            self.camera.elevation.to_degrees(),
        ));
    }

    /// Show colorbar legend.
    pub fn show_colorbar(&self, ui: &mut egui::Ui) {
        ui.separator();
        ui.label("Field Range");

        let lut = self.colormap.generate_lut(64);
        let available_height = ui.available_height().min(200.0).max(60.0);
        let bar_width = 20.0;

        let (rect, _) = ui.allocate_exact_size(
            egui::vec2(bar_width + 60.0, available_height),
            egui::Sense::hover(),
        );

        let painter = ui.painter_at(rect);
        let bar_rect = egui::Rect::from_min_size(rect.min, egui::vec2(bar_width, available_height));

        let n = lut.len();
        for i in 0..n {
            let frac_top = i as f32 / n as f32;
            let frac_bot = (i + 1) as f32 / n as f32;
            let color_idx = n - 1 - i;
            let c = lut[color_idx];
            let color = egui::Color32::from_rgb(c[0], c[1], c[2]);

            let y0 = bar_rect.top() + frac_top * bar_rect.height();
            let y1 = bar_rect.top() + frac_bot * bar_rect.height();
            painter.rect_filled(
                egui::Rect::from_min_max(
                    egui::pos2(bar_rect.left(), y0),
                    egui::pos2(bar_rect.right(), y1),
                ),
                0.0,
                color,
            );
        }

        painter.rect_stroke(
            bar_rect,
            0.0,
            egui::Stroke::new(1.0, egui::Color32::GRAY),
            egui::StrokeKind::Outside,
        );

        let text_x = bar_rect.right() + 4.0;
        let style = egui::TextStyle::Small.resolve(ui.style());
        painter.text(
            egui::pos2(text_x, bar_rect.top()),
            egui::Align2::LEFT_TOP,
            format!("{:.3}", self.field_range[1]),
            style.clone(),
            egui::Color32::WHITE,
        );
        painter.text(
            egui::pos2(text_x, bar_rect.center().y),
            egui::Align2::LEFT_CENTER,
            format!("{:.3}", (self.field_range[0] + self.field_range[1]) / 2.0),
            style.clone(),
            egui::Color32::WHITE,
        );
        painter.text(
            egui::pos2(text_x, bar_rect.bottom()),
            egui::Align2::LEFT_BOTTOM,
            format!("{:.3}", self.field_range[0]),
            style,
            egui::Color32::WHITE,
        );
    }

    // -----------------------------------------------------------------------
    // Real data loading (Milestone 7)
    // -----------------------------------------------------------------------

    /// Load a mesh from a .msh file and extract its surface for rendering.
    pub fn load_mesh_file(&mut self, msh_mesh: &emstudio_domain::msh_loader::MshMesh) {
        let mut field_mesh = surface_extraction::extract_surface(msh_mesh);
        let node_map = surface_extraction::build_node_to_vertex_map(msh_mesh, &field_mesh);

        // If we already have field data loaded, apply it
        if let Some(ref fld) = self.loaded_field {
            if let Ok(block) = fld.load_block(self.selected_frequency) {
                field_mapping::map_field_to_mesh(
                    &mut field_mesh,
                    &block,
                    self.field_component,
                    &node_map,
                );
            }
        }

        self.node_to_vertex = Some(node_map);
        self.loaded_mesh = Some(field_mesh);
        self.loaded_msh = Some(msh_mesh.clone());
        self.mode_dirty = true;
    }

    /// Load field data from an .emsfld file.
    pub fn load_field_file(&mut self, fld: emstudio_domain::emsfld_loader::EmsFldFile) {
        self.loaded_field = Some(fld);
        self.selected_frequency = 0;
        // Apply to loaded mesh if available
        self.apply_field_to_loaded_mesh();
    }

    /// Change the selected frequency index and update field visualization.
    pub fn set_frequency(&mut self, freq_idx: usize) {
        self.selected_frequency = freq_idx;
        self.apply_field_to_loaded_mesh();
    }

    /// Change the field component being visualized.
    pub fn set_field_component(&mut self, component: FieldComponent) {
        self.field_component = component;
        self.apply_field_to_loaded_mesh();
    }

    fn apply_field_to_loaded_mesh(&mut self) {
        let fld = match &self.loaded_field {
            Some(f) => f,
            None => return,
        };
        let mesh = match &mut self.loaded_mesh {
            Some(m) => m,
            None => return,
        };
        let node_map = match &self.node_to_vertex {
            Some(m) => m,
            None => return,
        };

        if let Ok(block) = fld.load_block(self.selected_frequency) {
            let range = field_mapping::map_field_to_mesh(
                mesh,
                &block,
                self.field_component,
                node_map,
            );
            self.field_range = range;

            // Upload to GPU
            if let Some(rs) = &self.render_state {
                let mut renderer = rs.renderer.write();
                if let Some(res) = renderer
                    .callback_resources
                    .get_mut::<FieldSceneResources>()
                {
                    res.pipeline.swap_mesh(&rs.device, mesh);
                }
            }
        }
    }

    /// Get the number of available frequencies from loaded field data.
    pub fn num_frequencies(&self) -> usize {
        self.loaded_field
            .as_ref()
            .map_or(0, |f| f.num_frequencies())
    }

    /// Get frequency values in Hz from loaded field data.
    pub fn frequencies_hz(&self) -> Vec<f64> {
        self.loaded_field
            .as_ref()
            .map_or(Vec::new(), |f| f.frequencies_hz())
    }

    /// Whether real data is loaded (vs synthetic demo data).
    pub fn has_real_data(&self) -> bool {
        self.loaded_mesh.is_some()
    }

    // -----------------------------------------------------------------------
    // Internal helpers
    // -----------------------------------------------------------------------

    fn apply_mode_switch(&mut self) {
        let rs = match &self.render_state {
            Some(r) => r.clone(),
            None => return,
        };

        // Use real data when available
        let has_real = self.has_real_data();

        match self.vis_mode {
            VisMode::Slice => {
                self.slice_dirty = true;
                self.update_slice_mesh();
                return;
            }
            VisMode::Isosurface => {
                self.iso_dirty = true;
                self.update_isosurface_mesh();
                return;
            }
            VisMode::MeshQuality => {
                self.update_mesh_quality();
                return;
            }
            _ => {}
        }

        let mesh = match self.vis_mode {
            VisMode::Surface | VisMode::Animation | VisMode::Wireframe => {
                if has_real {
                    self.loaded_mesh.as_ref()
                } else {
                    self.sphere_mesh.as_ref()
                }
            }
            VisMode::Arrows => {
                if has_real {
                    // For real data, also generate arrows from vector field
                    if let (Some(fld), Some(loaded_mesh), Some(node_map)) =
                        (&self.loaded_field, &self.loaded_mesh, &self.node_to_vertex)
                    {
                        if let Ok(block) = fld.load_block(self.selected_frequency) {
                            let mut mesh_copy = loaded_mesh.clone();
                            field_mapping::map_vector_field(&mut mesh_copy, &block, node_map);
                            let arrows = mesh_copy.generate_arrows(self.arrow_subsample);
                            let mut renderer = rs.renderer.write();
                            if let Some(res) = renderer.callback_resources.get_mut::<FieldSceneResources>() {
                                if let Some(ref mut ap) = res.arrow_pipeline {
                                    ap.upload_instances(&rs.queue, &arrows);
                                }
                            }
                        }
                    }
                    self.loaded_mesh.as_ref()
                } else {
                    // Upload synthetic arrows from cube mesh
                    if let Some(ref cube) = self.cube_mesh {
                        let arrows = cube.generate_arrows(self.arrow_subsample);
                        let mut renderer = rs.renderer.write();
                        if let Some(res) = renderer.callback_resources.get_mut::<FieldSceneResources>() {
                            if let Some(ref mut ap) = res.arrow_pipeline {
                                ap.upload_instances(&rs.queue, &arrows);
                            }
                        }
                    }
                    self.cube_mesh.as_ref()
                }
            }
            VisMode::FarField => {
                if self.loaded_far_field.is_some() {
                    // Generate mesh from real far-field data
                    let ff_data = self.loaded_far_field.as_ref().unwrap();
                    if let Some(mesh) = far_field::generate_pattern_mesh_from_data(ff_data, "GainTotal") {
                        self.far_field_mesh = Some(mesh);
                    }
                }
                self.far_field_mesh.as_ref()
            }
            VisMode::Slice | VisMode::Isosurface | VisMode::MeshQuality => unreachable!(),
        };

        if let Some(m) = mesh {
            self.field_range = m.field_range;
            if self.vis_mode == VisMode::Animation {
                let source = if has_real { &self.loaded_mesh } else { &self.sphere_mesh };
                if let Some(src) = source {
                    if let Some(ref imag) = src.field_imag {
                        self.field_range = PhaseAnimator::envelope_range(
                            &src.vertices.iter().map(|v| v.field_value).collect::<Vec<_>>(),
                            imag,
                        );
                    }
                }
            }
            let mut renderer = rs.renderer.write();
            if let Some(res) = renderer.callback_resources.get_mut::<FieldSceneResources>() {
                res.pipeline.swap_mesh(&rs.device, m);
            }
        }
    }

    fn update_slice_mesh(&mut self) {
        let rs = match &self.render_state {
            Some(r) => r.clone(),
            None => return,
        };

        let mesh = if let (Some(msh), Some(fld)) = (&self.loaded_msh, &self.loaded_field) {
            // Use real data: intersect slice plane with tetrahedra
            if let Ok(block) = fld.load_block(self.selected_frequency) {
                slice::generate_slice_mesh_from_tets(
                    msh,
                    &block,
                    self.field_component,
                    self.slice_axis,
                    self.slice_z,
                )
            } else {
                slice::generate_slice_mesh(
                    self.slice_axis,
                    self.slice_z,
                    1.0,
                    40,
                    &slice::synthetic_volume_field,
                )
            }
        } else {
            slice::generate_slice_mesh(
                self.slice_axis,
                self.slice_z,
                1.0,
                40,
                &slice::synthetic_volume_field,
            )
        };

        self.field_range = mesh.field_range;
        let mut renderer = rs.renderer.write();
        if let Some(res) = renderer.callback_resources.get_mut::<FieldSceneResources>() {
            res.pipeline.swap_mesh(&rs.device, &mesh);
        }
    }

    fn update_isosurface_mesh(&mut self) {
        let rs = match &self.render_state {
            Some(r) => r.clone(),
            None => return,
        };

        let mesh = if let (Some(msh), Some(fld)) = (&self.loaded_msh, &self.loaded_field) {
            if let Ok(block) = fld.load_block(self.selected_frequency) {
                isosurface::extract_isosurface(
                    msh,
                    &block,
                    self.field_component,
                    self.iso_threshold as f64,
                )
            } else {
                return;
            }
        } else {
            // No real data available for isosurface
            return;
        };

        self.field_range = mesh.field_range;
        let mut renderer = rs.renderer.write();
        if let Some(res) = renderer.callback_resources.get_mut::<FieldSceneResources>() {
            res.pipeline.swap_mesh(&rs.device, &mesh);
        }
    }

    fn update_mesh_quality(&mut self) {
        let rs = match &self.render_state {
            Some(r) => r.clone(),
            None => return,
        };

        let mesh = if let Some(ref msh) = self.loaded_msh {
            mesh_quality::compute_mesh_quality(msh, self.quality_metric)
        } else {
            return;
        };

        self.field_range = mesh.field_range;
        let mut renderer = rs.renderer.write();
        if let Some(res) = renderer.callback_resources.get_mut::<FieldSceneResources>() {
            res.pipeline.swap_mesh(&rs.device, &mesh);
        }
    }

    fn apply_phase_animation(&mut self) {
        let rs = match &self.render_state {
            Some(r) => r.clone(),
            None => return,
        };

        let source = if self.has_real_data() {
            self.loaded_mesh.as_ref()
        } else {
            self.sphere_mesh.as_ref()
        };
        let source = match source {
            Some(m) => m,
            None => return,
        };
        let imag = match &source.field_imag {
            Some(im) => im,
            None => return,
        };

        let real_values: Vec<f32> = source.vertices.iter().map(|v| v.field_value).collect();
        let animated = self.animator.apply(&real_values, imag);

        let mut verts = source.vertices.clone();
        for (v, &new_val) in verts.iter_mut().zip(animated.iter()) {
            v.field_value = new_val;
        }

        let renderer = rs.renderer.read();
        if let Some(res) = renderer.callback_resources.get::<FieldSceneResources>() {
            res.pipeline.update_vertices(&rs.queue, &verts);
        }
    }

    // -----------------------------------------------------------------------
    // Far-field data loading
    // -----------------------------------------------------------------------

    /// Load real far-field data for 3D radiation pattern visualization.
    pub fn load_far_field_data(&mut self, data: emstudio_domain::result_store::FarFieldData) {
        self.loaded_far_field = Some(data);
        if self.vis_mode == VisMode::FarField {
            self.mode_dirty = true;
        }
    }

    // -----------------------------------------------------------------------
    // Q3D field overlay loading
    // -----------------------------------------------------------------------

    /// Load a Q3D field overlay (current density, charge distribution, ohmic loss).
    /// Routes to the correct surface extraction method based on the field quantity.
    pub fn load_q3d_overlay(
        &mut self,
        quantity: emstudio_domain::field_overlay::FieldQuantity,
        fld: emstudio_domain::emsfld_loader::EmsFldFile,
        msh: &emstudio_domain::msh_loader::MshMesh,
    ) {
        use emstudio_domain::field_overlay::FieldQuantity;

        let mut field_mesh = match quantity {
            FieldQuantity::Jsurf | FieldQuantity::ChargeDistribution => {
                surface_extraction::extract_triangles(msh)
            }
            _ => surface_extraction::extract_surface(msh),
        };

        let node_map = surface_extraction::build_node_to_vertex_map(msh, &field_mesh);

        if let Ok(block) = fld.load_block(0) {
            field_mapping::map_field_to_mesh(
                &mut field_mesh,
                &block,
                self.field_component,
                &node_map,
            );
        }

        self.node_to_vertex = Some(node_map);
        self.loaded_mesh = Some(field_mesh);
        self.loaded_msh = Some(msh.clone());
        self.loaded_field = Some(fld);
        self.mode_dirty = true;
    }

    // -----------------------------------------------------------------------
    // Picking
    // -----------------------------------------------------------------------

    /// Perform a pick at screen coordinates (relative to viewport rect).
    pub fn pick_at(&mut self, x: f32, y: f32, viewport: [f32; 2]) {
        let mesh = if self.has_real_data() {
            self.loaded_mesh.as_ref()
        } else {
            self.sphere_mesh.as_ref()
        };
        let mesh = match mesh {
            Some(m) => m,
            None => return,
        };

        let aspect = viewport[0] / viewport[1].max(1.0);
        let (ray_origin, ray_dir) = self.camera.screen_to_ray(x, y, viewport, aspect);
        self.probe_result = picking::pick_field(mesh, ray_origin, ray_dir);
    }
}
