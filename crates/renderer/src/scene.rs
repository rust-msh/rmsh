use rcad_render::{Camera, DisplayMode, Tessellator, WgpuRenderer};
use rcad_kernel::{BRep, Shell, Solid, Vertex, Wire, Face};
use rmsh_geo::extract::{PointData, SurfaceData, WireframeData};

/// Re-export rcad-render Camera as the orbit camera.
pub use rcad_render::Camera as OrbitCamera;

/// Extension methods for `Camera` (removed from rcad-render upstream).
pub trait CameraExt {
    fn rotate(&mut self, delta_yaw: f32, delta_pitch: f32);
    fn zoom(&mut self, delta: f32);
    fn pan(&mut self, delta_x: f32, delta_y: f32);
    fn fit_to_bbox(&mut self, center: [f32; 3], diagonal: f32);
    fn set_isometric(&mut self);
    fn toggle_projection(&mut self);
    fn orthographic(&self) -> bool;
}

impl CameraExt for Camera {
    fn rotate(&mut self, delta_yaw: f32, delta_pitch: f32) {
        self.rot_y += delta_yaw;
        self.rot_x = (self.rot_x + delta_pitch).clamp(
            -std::f32::consts::FRAC_PI_2 + 0.01,
            std::f32::consts::FRAC_PI_2 - 0.01,
        );
    }

    fn zoom(&mut self, delta: f32) {
        self.distance = (self.distance * (1.0 - delta)).max(0.01);
    }

    fn pan(&mut self, delta_x: f32, delta_y: f32) {
        self.pan_pixels(delta_x, delta_y);
    }

    fn fit_to_bbox(&mut self, center: [f32; 3], diagonal: f32) {
        self.target = glam::Vec3::new(center[0], center[1], center[2]);
        self.distance = diagonal * 1.5;
    }

    fn set_isometric(&mut self) {
        self.rot_x = 0.6154_8246; // atan(1/sqrt(2))
        self.rot_y = std::f32::consts::FRAC_PI_4;
    }

    fn toggle_projection(&mut self) {
        // Perspective-only camera; no-op.
    }

    fn orthographic(&self) -> bool {
        false
    }
}

/// Rendering configuration — controls visibility and Gmsh-style visual parameters.
#[derive(Debug, Clone)]
pub struct RenderConfig {
    // ── Visibility ─────────────────────────────────────────────────────────────
    pub show_nodes: bool,
    pub show_edges: bool,
    pub show_faces: bool,
    pub show_volumes: bool,
    pub show_axes: bool,
    pub show_scale_ruler: bool,

    // ── Background ─────────────────────────────────────────────────────────────
    /// Top background color (gradient top — used for egui background fill).
    pub bg_color_top: [f32; 3],
    /// Bottom background color (gradient bottom).
    pub bg_color_bottom: [f32; 3],

    // ── Surface ────────────────────────────────────────────────────────────────
    /// Base face/surface color RGB.
    pub face_color: [f32; 3],
    /// Surface opacity [0, 1].
    pub surface_opacity: f32,

    // ── Wireframe ──────────────────────────────────────────────────────────────
    /// Edge/wireframe color RGB.
    pub edge_color: [f32; 3],

    // ── Nodes ──────────────────────────────────────────────────────────────────
    /// Node/point color RGB.
    pub node_color: [f32; 3],

    // ── Highlights ─────────────────────────────────────────────────────────────
    pub face_highlight_color: [f32; 4],
    pub edge_highlight_color: [f32; 4],
}

impl Default for RenderConfig {
    fn default() -> Self {
        // Gmsh dark-theme defaults
        Self {
            show_nodes: false,
            show_edges: true,
            show_faces: true,
            show_volumes: true,
            show_axes: true,
            show_scale_ruler: true,

            // Sky gradient: white at bottom, light blue at top
            bg_color_top: [0.53, 0.81, 0.92],      // Light sky blue
            bg_color_bottom: [1.0, 1.0, 1.0],      // White

            // Gmsh: light steel-blue surfaces
            face_color: [0.75, 0.85, 0.95],
            surface_opacity: 0.92,

            // Dark gray edges for contrast on light background
            edge_color: [0.2, 0.2, 0.25],

            // Gmsh: orange node dots
            node_color: [1.0, 0.60, 0.10],

            face_highlight_color: [1.0, 0.50, 0.05, 0.50],
            edge_highlight_color: [1.0, 0.90, 0.10, 1.0],
        }
    }
}

impl RenderConfig {
    fn to_display_mode(&self) -> DisplayMode {
        match (self.show_faces, self.show_edges) {
            (true, true) => DisplayMode::SolidWithEdges,
            (true, false) => DisplayMode::Solid,
            (false, true) => DisplayMode::Wireframe,
            (false, false) => DisplayMode::Wireframe,
        }
    }

    fn model_color_rgba(&self) -> [f32; 4] {
        [
            self.face_color[0],
            self.face_color[1],
            self.face_color[2],
            self.surface_opacity,
        ]
    }
}

/// The 3D scene — wraps `rcad_render::WgpuRenderer` and bridges rmsh geometry types.
pub struct Scene {
    pub camera: Camera,
    pub config: RenderConfig,
    pub renderer: WgpuRenderer,
    pub show_axis_gizmo: bool,
}

impl Scene {
    pub fn new(device: &wgpu::Device, target_format: wgpu::TextureFormat) -> Self {
        Self {
            camera: Camera::new(),
            config: RenderConfig::default(),
            renderer: WgpuRenderer::new(device, target_format),
            show_axis_gizmo: true,
        }
    }

    /// Upload mesh data to GPU from extracted geometry.
    pub fn upload_mesh(
        &mut self,
        device: &wgpu::Device,
        surface: &SurfaceData,
        wireframe: &WireframeData,
        _points: &PointData,
    ) {
        let brep = surface_wireframe_to_brep(surface, wireframe);
        let mesh = Tessellator::tessellate(&brep);
        self.renderer.upload_mesh(device, &mesh);
    }

    /// Upload highlight geometry.
    pub fn upload_highlight(
        &mut self,
        device: &wgpu::Device,
        surface: Option<&SurfaceData>,
        wireframe: Option<&WireframeData>,
    ) {
        let face_mesh = surface.map(|s| {
            let brep = surface_to_brep(s);
            Tessellator::tessellate(&brep)
        });
        let edge_mesh = wireframe.map(|w| {
            let brep = wireframe_to_brep(w);
            Tessellator::tessellate(&brep)
        });
        self.renderer.upload_highlights(device, face_mesh.as_ref(), edge_mesh.as_ref());
    }

    /// Clear highlight geometry.
    pub fn clear_highlight(&mut self, device: &wgpu::Device) {
        self.renderer.upload_highlights(device, None, None);
    }

    /// Update camera and style uniforms (called from egui prepare callback).
    pub fn update_uniforms(&self, queue: &wgpu::Queue, width: u32, height: u32) {
        if width == 0 || height == 0 {
            return;
        }
        let aspect = width as f32 / height as f32;
        self.renderer.update_camera(queue, &self.camera, aspect);
        self.renderer.update_axis_gizmo_camera(queue, &self.camera);
        self.renderer
            .set_model_color(queue, self.config.model_color_rgba());
    }

    /// Sync draw flags from config to renderer.
    pub fn sync_config(&mut self) {
        self.renderer
            .set_display_mode(self.config.to_display_mode());
        self.renderer.set_show_axes(self.config.show_axes);
        self.renderer.set_show_grid(self.config.show_scale_ruler);
    }

    /// Set the scene axes scale factor (default 0.3).
    pub fn set_scene_axes_scale(&self, _scale: f32) {
        // rcad-render no longer exposes per-scene axes scale mutators.
        // Keep this API for compatibility with viewer UI calls.
    }

    /// Get the scene axes scale factor.
    pub fn scene_axes_scale(&self) -> f32 {
        0.3
    }

    /// Enable or disable corner axis gizmo.
    pub fn set_show_axis_gizmo(&mut self, show: bool) {
        self.show_axis_gizmo = show;
    }

    /// Draw corner axis gizmo when enabled.
    pub fn draw_axis_gizmo_in_render_pass(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        viewport_origin_px: [u32; 2],
        viewport_size_px: [u32; 2],
    ) {
        if !self.show_axis_gizmo {
            return;
        }
        self.renderer.draw_axis_gizmo_in_render_pass(
            render_pass,
            viewport_origin_px,
            viewport_size_px,
            false,
        );
    }

    /// Draw into an active render pass (called from egui paint callback).
    pub fn draw_in_render_pass(&self, render_pass: &mut wgpu::RenderPass<'_>) {
        self.renderer.draw_in_render_pass(render_pass, false);
    }
}

// ── Geometry conversion helpers ───────────────────────────────────────────────

/// Convert SurfaceData + WireframeData into a BRep for rcad-render tessellation.
fn surface_wireframe_to_brep(surface: &SurfaceData, wireframe: &WireframeData) -> BRep {
    let vertices: Vec<Vertex> = surface
        .positions
        .iter()
        .map(|p| Vertex {
            point: glam::DVec3::new(p[0] as f64, p[1] as f64, p[2] as f64),
        })
        .collect();

    let triangles: Vec<[usize; 3]> = surface.indices
        .chunks_exact(3)
        .map(|c| [c[0] as usize, c[1] as usize, c[2] as usize])
        .collect();

    let edges: Vec<rcad_kernel::Edge> = wireframe.indices
        .chunks_exact(2)
        .map(|c| rcad_kernel::Edge { start: c[0] as usize, end: c[1] as usize })
        .collect();

    let face = Face {
        outer_wire: Wire { edges: Vec::new() },
        inner_wires: Vec::new(),
        normal: glam::DVec3::Z,
        triangles,
        mesh_dirty: true,
        sample_point: None,
        surface_idx: Some(0),
    };

    BRep {
        vertices,
        edges,
        solids: vec![Solid {
            shells: vec![Shell { faces: vec![face] }],
        }],
        geom: rcad_kernel::GeomStore::default(),
        compound: None,
        compsolid: None,
    }
}

fn surface_to_brep(surface: &SurfaceData) -> BRep {
    let empty_wireframe = WireframeData { positions: Vec::new(), indices: Vec::new() };
    surface_wireframe_to_brep(surface, &empty_wireframe)
}

fn wireframe_to_brep(wireframe: &WireframeData) -> BRep {
    let vertices: Vec<Vertex> = wireframe
        .positions
        .iter()
        .map(|p| Vertex {
            point: glam::DVec3::new(p[0] as f64, p[1] as f64, p[2] as f64),
        })
        .collect();

    let edges: Vec<rcad_kernel::Edge> = wireframe.indices
        .chunks_exact(2)
        .map(|c| rcad_kernel::Edge { start: c[0] as usize, end: c[1] as usize })
        .collect();

    BRep {
        vertices,
        edges,
        solids: Vec::new(),
        geom: rcad_kernel::GeomStore::default(),
        compound: None,
        compsolid: None,
    }
}
