mod app;
mod io;
mod viewport;

#[cfg(not(target_arch = "wasm32"))]
use std::path::PathBuf;

#[cfg(not(target_arch = "wasm32"))]
use rmsh_model::Mesh;

pub use app::RmshApp;

#[derive(Debug, Clone, Default)]
pub struct ViewerConfig {
    pub show_nodes: Option<bool>,
    pub show_edges: Option<bool>,
    pub show_faces: Option<bool>,
    pub show_volumes: Option<bool>,
}

#[cfg(not(target_arch = "wasm32"))]
pub fn run_native_viewer(
    startup_path: Option<PathBuf>,
    initial_mesh: Option<(Mesh, String)>,
    initial_config: Option<ViewerConfig>,
) -> eframe::Result {
    let _ = env_logger::try_init();

    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1280.0, 800.0])
            .with_title("rmsh - Finite Element Mesh Viewer"),
        renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };

    eframe::run_native(
        "rmsh",
        native_options,
        Box::new(move |cc| {
            Ok(Box::new(RmshApp::new_with_inputs(
                cc,
                startup_path.clone(),
                initial_mesh.clone(),
                initial_config.clone(),
            )))
        }),
    )
}