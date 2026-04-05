use eframe::{egui, NativeOptions};
use emstudio_render::{FieldMesh, FieldSceneState};

// ---------------------------------------------------------------------------
// Pane definition
// ---------------------------------------------------------------------------

enum Pane {
    Viewport,
    Controls,
    Colorbar,
}

// ---------------------------------------------------------------------------
// Behavior: bridges egui_tiles with our FieldSceneState
// ---------------------------------------------------------------------------

struct FieldVisBehavior<'a> {
    scene: &'a mut FieldSceneState,
}

impl<'a> egui_tiles::Behavior<Pane> for FieldVisBehavior<'a> {
    fn tab_title_for_pane(&mut self, pane: &Pane) -> egui::WidgetText {
        match pane {
            Pane::Viewport => "3D Viewport".into(),
            Pane::Controls => "Controls".into(),
            Pane::Colorbar => "Colorbar".into(),
        }
    }

    fn pane_ui(
        &mut self,
        ui: &mut egui::Ui,
        _tile_id: egui_tiles::TileId,
        pane: &mut Pane,
    ) -> egui_tiles::UiResponse {
        match pane {
            Pane::Viewport => self.scene.show_viewport(ui),
            Pane::Controls => self.scene.show_controls(ui),
            Pane::Colorbar => self.scene.show_colorbar(ui),
        }
        Default::default()
    }
}

// ---------------------------------------------------------------------------
// App
// ---------------------------------------------------------------------------

struct FieldVisApp {
    scene: FieldSceneState,
    tree: egui_tiles::Tree<Pane>,
}

impl FieldVisApp {
    fn new(cc: &eframe::CreationContext) -> Self {
        let mesh = FieldMesh::uv_sphere(64, 128, 1.0);
        let mut scene = FieldSceneState::new();
        if let Some(rs) = &cc.wgpu_render_state {
            scene.init_gpu(rs, &mesh);
        }

        let tree = build_tree();
        Self { scene, tree }
    }
}

fn build_tree() -> egui_tiles::Tree<Pane> {
    let mut tiles = egui_tiles::Tiles::default();

    let viewport_id = tiles.insert_pane(Pane::Viewport);
    let controls_id = tiles.insert_pane(Pane::Controls);
    let colorbar_id = tiles.insert_pane(Pane::Colorbar);

    // Right column: Controls on top, Colorbar on bottom
    let right_col = {
        let mut linear = egui_tiles::Linear::new(
            egui_tiles::LinearDir::Vertical,
            vec![controls_id, colorbar_id],
        );
        linear.shares.set_share(controls_id, 0.5);
        linear.shares.set_share(colorbar_id, 0.5);
        tiles.insert_container(egui_tiles::Container::Linear(linear))
    };

    // Root: Viewport (left, 75%) | Right column (25%)
    let root = {
        let mut linear = egui_tiles::Linear::new(
            egui_tiles::LinearDir::Horizontal,
            vec![viewport_id, right_col],
        );
        linear.shares.set_share(viewport_id, 3.0);
        linear.shares.set_share(right_col, 1.0);
        tiles.insert_container(egui_tiles::Container::Linear(linear))
    };

    egui_tiles::Tree::new("field_vis_tree", root, tiles)
}

impl eframe::App for FieldVisApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            let mut behavior = FieldVisBehavior {
                scene: &mut self.scene,
            };
            self.tree.ui(&mut behavior, ui);
        });
    }
}

fn main() -> eframe::Result<()> {
    let options = NativeOptions {
        renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };
    eframe::run_native(
        "EMStudio - Field Visualization",
        options,
        Box::new(|cc| Ok(Box::new(FieldVisApp::new(cc)))),
    )
}
