pub mod animation;
pub mod arrow_pipeline;
pub mod camera;
pub mod colormap;
pub mod far_field;
pub mod field_mapping;
pub mod field_pipeline;
pub mod isosurface;
pub mod mesh_data;
pub mod mesh_quality;
pub mod picking;
pub mod scene;
pub mod screenshot;
pub mod slice;
pub mod surface_extraction;

pub use camera::{OrbitCamera, ViewPreset};
pub use colormap::ColormapType;
pub use mesh_data::{FieldMesh, FieldVertex};
pub use scene::{FieldSceneState, VisMode};

use std::collections::HashMap;

use egui::{Sense, Ui, Vec2};
use rcad_kernel::BRep;
use rcad_render::{
    build_faces_highlight_mesh, build_edges_highlight_mesh, merge_meshes,
    Camera, Mesh, SelectionState, Tessellator, WgpuRenderer,
};

// ─── Render Callback (egui_wgpu) ────────────────────────────────────────────

struct GeometryRenderCallback {
    mesh: Mesh,
    camera: Camera,
    aspect: f32,
    brep: BRep,
    selected_faces: Vec<usize>,
    selected_edges: Vec<usize>,
}

impl egui_wgpu::CallbackTrait for GeometryRenderCallback {
    fn prepare(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        _screen_descriptor: &egui_wgpu::ScreenDescriptor,
        _egui_encoder: &mut wgpu::CommandEncoder,
        callback_resources: &mut egui_wgpu::CallbackResources,
    ) -> Vec<wgpu::CommandBuffer> {
        let Some(renderer) = callback_resources.get::<WgpuRenderer>() else {
            return Vec::new();
        };

        let face_mesh = build_faces_highlight_mesh(&self.brep, &self.selected_faces);
        let edge_mesh = build_edges_highlight_mesh(&self.brep, &self.selected_edges);
        renderer.upload_highlights(device, face_mesh.as_ref(), edge_mesh.as_ref());
        renderer.prepare_scene(device, queue, &self.mesh, &self.camera, self.aspect);
        Vec::new()
    }

    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        render_pass: &mut wgpu::RenderPass<'static>,
        callback_resources: &egui_wgpu::CallbackResources,
    ) {
        let Some(renderer) = callback_resources.get::<WgpuRenderer>() else {
            return;
        };
        renderer.draw_in_render_pass(render_pass, false);
    }
}

// ─── SceneViewport ──────────────────────────────────────────────────────────

pub struct SceneViewport {
    pub title: String,
    camera: Camera,
    selection: SelectionState,
    renderer_ready: bool,
    cached_mesh: Option<Mesh>,
    cached_brep: Option<BRep>,
    brep_generation: u64,
}

impl Default for SceneViewport {
    fn default() -> Self {
        Self {
            title: "3D View".to_string(),
            camera: Camera::new(),
            selection: SelectionState::default(),
            renderer_ready: false,
            cached_mesh: None,
            cached_brep: None,
            brep_generation: 0,
        }
    }
}

impl SceneViewport {
    /// Initialize the WgpuRenderer into egui_wgpu callback resources.
    /// Must be called once from the eframe App's creation context or first frame.
    pub fn init_renderer(&mut self, device: &wgpu::Device, format: wgpu::TextureFormat, resources: &mut egui_wgpu::CallbackResources) {
        if self.renderer_ready {
            return;
        }
        let renderer = WgpuRenderer::new(device, format);
        resources.insert(renderer);
        self.renderer_ready = true;
    }

    /// Main UI entry point. `breps` is the live geometry from GeometryEngine.
    /// `generation` should change whenever the geometry is modified, to trigger
    /// re-tessellation.
    pub fn ui(
        &mut self,
        ui: &mut Ui,
        breps: &HashMap<String, BRep>,
        generation: u64,
    ) {
        // Re-tessellate when geometry changes
        if generation != self.brep_generation || self.cached_mesh.is_none() {
            self.brep_generation = generation;
            self.retessellate(breps);
        }

        let desired_size = Vec2::new(ui.available_width(), ui.available_height().max(180.0));
        let (rect, response) = ui.allocate_exact_size(desired_size, Sense::click_and_drag());

        // ── Camera interaction ──
        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll != 0.0 {
            self.camera.distance -= scroll * 0.005 * self.camera.distance;
            self.camera.distance = self.camera.distance.clamp(0.1, 200.0);
            ui.ctx().request_repaint();
        }

        if response.dragged() {
            let delta = ui.input(|i| i.pointer.delta());
            let pan_with_middle =
                ui.input(|i| i.pointer.button_down(egui::PointerButton::Middle));
            let rotate_with_alt =
                ui.input(|i| i.modifiers.alt && i.pointer.button_down(egui::PointerButton::Primary));
            let rotate_with_right =
                ui.input(|i| i.pointer.button_down(egui::PointerButton::Secondary));

            if pan_with_middle {
                self.camera.pan_pixels(delta.x, delta.y);
                ui.ctx().request_repaint();
            } else if rotate_with_alt || rotate_with_right {
                self.camera.rot_y += delta.x * 0.008;
                self.camera.rot_x += delta.y * 0.008;
                ui.ctx().request_repaint();
            }
        }

        // Click for selection
        if response.clicked() {
            if let Some(pointer) = response.interact_pointer_pos() {
                let local = pointer - rect.min;
                let alt_down = ui.input(|i| i.modifiers.alt);
                if !alt_down {
                    if let Some(brep) = &self.cached_brep {
                        let aspect = rect.width() / rect.height().max(1.0);
                        self.selection.click_at(
                            brep,
                            &self.camera,
                            aspect,
                            [rect.width(), rect.height()],
                            [local.x, local.y],
                            rcad_render::DEFAULT_EDGE_PICK_RADIUS_PX,
                        );
                    }
                }
            }
        }

        // Hover for highlighting
        if let Some(pointer) = ui.input(|i| i.pointer.hover_pos()) {
            if rect.contains(pointer) {
                let local = pointer - rect.min;
                if let Some(brep) = &self.cached_brep {
                    let aspect = rect.width() / rect.height().max(1.0);
                    self.selection.hover_at(
                        brep,
                        &self.camera,
                        aspect,
                        [rect.width(), rect.height()],
                        [local.x, local.y],
                        rcad_render::DEFAULT_EDGE_PICK_RADIUS_PX,
                    );
                }
            } else {
                self.selection.clear_hover();
            }
        }

        // ── Render via egui_wgpu callback ──
        if self.renderer_ready {
            if let Some(mesh) = &self.cached_mesh {
                let aspect = rect.width() / rect.height().max(1.0);
                let cb = GeometryRenderCallback {
                    mesh: mesh.clone(),
                    camera: self.camera,
                    aspect,
                    brep: self.cached_brep.clone().unwrap_or_default(),
                    selected_faces: self.selection.highlighted_faces(),
                    selected_edges: self.selection.highlighted_edges(),
                };
                ui.painter()
                    .add(egui_wgpu::Callback::new_paint_callback(rect, cb));
            }
        } else {
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 6.0, egui::Color32::from_rgb(16, 24, 34));
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "Initializing 3D viewport...",
                egui::TextStyle::Body.resolve(ui.style()),
                egui::Color32::WHITE,
            );
        }
    }

    pub fn camera(&self) -> &Camera {
        &self.camera
    }

    pub fn camera_mut(&mut self) -> &mut Camera {
        &mut self.camera
    }

    pub fn selection(&self) -> &SelectionState {
        &self.selection
    }

    pub fn selection_mut(&mut self) -> &mut SelectionState {
        &mut self.selection
    }

    fn retessellate(&mut self, breps: &HashMap<String, BRep>) {
        if breps.is_empty() {
            self.cached_mesh = Some(Mesh {
                vertices: Vec::new(),
                indices: Vec::new(),
                line_indices: Vec::new(),
            });
            self.cached_brep = None;
            return;
        }

        let meshes: Vec<Mesh> = breps.values().map(|b| Tessellator::tessellate(b)).collect();
        let mesh_refs: Vec<&Mesh> = meshes.iter().collect();
        let merged = merge_meshes(&mesh_refs).unwrap_or(Mesh {
            vertices: Vec::new(),
            indices: Vec::new(),
            line_indices: Vec::new(),
        });
        self.cached_mesh = Some(merged);

        // For picking: use first BRep (multi-body picking needs per-object dispatch later)
        self.cached_brep = breps.values().next().cloned();
    }
}
