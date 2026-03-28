use eframe::egui;
use egui_dock::{DockArea, DockState, TabViewer};

use emstudio_components::dock;
use emstudio_components::ribbon::{RibbonAction, show_ribbon};
use emstudio_domain::{Project, SimulationStatus};
use emstudio_infra::{Backend, RunMode, default_backend};
use emstudio_render::SceneViewport;

#[derive(Debug, Clone, PartialEq, Eq)]
enum CenterTab {
    Modeling,
    Result,
    Log,
}

impl CenterTab {
    fn title(&self) -> &'static str {
        match self {
            Self::Modeling => "Modeling",
            Self::Result => "Result",
            Self::Log => "Log",
        }
    }
}

struct CenterTabViewer<'a> {
    project: &'a Project,
    viewport: &'a mut SceneViewport,
    log_text: &'a str,
}

impl TabViewer for CenterTabViewer<'_> {
    type Tab = CenterTab;

    fn title(&mut self, tab: &mut Self::Tab) -> egui::WidgetText {
        tab.title().into()
    }

    fn ui(&mut self, ui: &mut egui::Ui, tab: &mut Self::Tab) {
        match tab {
            CenterTab::Modeling => self.viewport.ui(ui, self.project),
            CenterTab::Result => {
                ui.heading("Result Preview");
                ui.separator();
                if let Some(result) = &self.project.last_result {
                    ui.label(format!("Converged: {}", result.converged));
                    ui.label(result.field_preview.as_str());
                } else {
                    ui.label("No result yet.");
                }
            }
            CenterTab::Log => {
                ui.heading("Run Log");
                ui.separator();
                ui.monospace(self.log_text);
            }
        }
    }
}

pub struct EmStudioApp {
    project: Project,
    backend: Box<dyn Backend>,
    viewport: SceneViewport,
    dock_state: DockState<CenterTab>,
    status_text: String,
    log_text: String,
}

impl EmStudioApp {
    fn mode_label(mode: RunMode) -> &'static str {
        match mode {
            RunMode::Standalone => "standalone",
            RunMode::Cloud => "cloud",
        }
    }

    pub fn new(mode: RunMode) -> Self {
        let project = Project::default();
        let mut dock_state = DockState::new(vec![CenterTab::Modeling]);
        let [main, right] = dock_state.main_surface_mut().split_right(
            egui_dock::NodeIndex::root(),
            0.35,
            vec![CenterTab::Result],
        );
        let _ = main;
        let _bottom = dock_state
            .main_surface_mut()
            .split_below(right, 0.5, vec![CenterTab::Log]);

        Self {
            project,
            backend: default_backend(mode),
            viewport: SceneViewport::default(),
            dock_state,
            status_text: format!("Ready ({})", Self::mode_label(mode)),
            log_text: format!("[boot] EmStudio shell started (mode={})", Self::mode_label(mode)),
        }
    }

    pub fn new_default() -> Self {
        Self::new(RunMode::Standalone)
    }

    fn on_ribbon_action(&mut self, action: RibbonAction) {
        match action {
            RibbonAction::NewProject => {
                self.project = Project::default();
                self.status_text = "New project created".to_string();
                self.log_text.push_str("\n[action] new project");
            }
            RibbonAction::OpenProject => {
                self.status_text = "Open is not wired yet".to_string();
                self.log_text.push_str("\n[action] open project (todo)");
            }
            RibbonAction::SaveProject => match self.backend.save_project(self.project.clone()) {
                Ok(()) => {
                    self.status_text = "Project saved".to_string();
                    self.log_text.push_str("\n[action] save project");
                }
                Err(err) => {
                    self.status_text = format!("Save failed: {err}");
                }
            },
            RibbonAction::Solve => {
                self.project.status = SimulationStatus::Solving;
                match self.backend.solve(&self.project) {
                    Ok(result) => {
                        self.project.status = SimulationStatus::Finished;
                        self.project.last_result = Some(result);
                        self.status_text = "Solve completed".to_string();
                        self.log_text.push_str("\n[action] solve project");
                    }
                    Err(err) => {
                        self.project.status = SimulationStatus::Failed;
                        self.status_text = format!("Solve failed: {err}");
                    }
                }
            }
        }
    }
}

impl eframe::App for EmStudioApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("New Project").clicked() {
                        self.on_ribbon_action(RibbonAction::NewProject);
                        ui.close_menu();
                    }
                    if ui.button("Save Project").clicked() {
                        self.on_ribbon_action(RibbonAction::SaveProject);
                        ui.close_menu();
                    }
                });
            });
        });

        egui::TopBottomPanel::top("ribbon").show(ctx, |ui| {
            if let Some(action) = show_ribbon(ui) {
                self.on_ribbon_action(action);
            }
        });

        egui::SidePanel::left("left_dock")
            .resizable(true)
            .default_width(210.0)
            .show(ctx, |ui| {
                dock::left_panel(ui, &self.project);
            });

        egui::SidePanel::right("right_dock")
            .resizable(true)
            .default_width(240.0)
            .show(ctx, |ui| {
                dock::right_panel(ui, &self.project);
            });

        egui::TopBottomPanel::bottom("bottom_status")
            .resizable(true)
            .default_height(30.0)
            .show(ctx, |ui| {
                dock::bottom_panel(ui, &self.status_text);
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            let mut viewer = CenterTabViewer {
                project: &self.project,
                viewport: &mut self.viewport,
                log_text: &self.log_text,
            };
            DockArea::new(&mut self.dock_state).show_inside(ui, &mut viewer);
        });
    }
}
