use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::mpsc;

use eframe::egui;
use egui_dock::{DockArea, DockState, TabViewer};

use emstudio_components::dock::{self, LeftDockContext};
use emstudio_components::draw_dialog::DrawDialog;
use emstudio_components::message_manager::{self, BottomTab, MessageEntry};
use emstudio_components::report_panel::ReportPanel;
use emstudio_components::ribbon::{
    RibbonAction, RibbonState, RibbonTab, build_default_tabs, show_ribbon,
};
use emstudio_components::status_bar::{self, StatusBarState};
use emstudio_components::{LeftPanelTab, menu_bar, qat};
use emstudio_domain::{Edition, Project, SimulationStatus};
use emstudio_domain::geometry_engine::GeometryEngine;
use emstudio_domain::result_store::ResultDataStore;
use emstudio_domain::variable::Variable;
use emstudio_infra::{Backend, RunMode, default_backend};
#[cfg(not(target_arch = "wasm32"))]
use emstudio_infra::{
    export_hfss_project_to_file, export_q3d_project_to_file, import_hfss_project_from_file,
    import_q3d_project_from_file, load_project_from_file, save_project_to_file,
};
use emstudio_render::SceneViewport;

// ---------------------------------------------------------------------------
// Async file dialog result
// ---------------------------------------------------------------------------

enum FileDialogResult {
    OpenFile(PathBuf),
    SaveFile(PathBuf),
    ImportHfssAedt(PathBuf),
    ImportQ3dAedt(PathBuf),
    ExportHfssPy(PathBuf),
    ExportQ3dPy(PathBuf),
}

fn spawn_future<F: std::future::Future<Output = ()> + 'static>(f: F) {
    #[cfg(target_arch = "wasm32")]
    wasm_bindgen_futures::spawn_local(f);

    #[cfg(not(target_arch = "wasm32"))]
    pollster::block_on(f);
}

// ---------------------------------------------------------------------------
// Center tabs (dock area)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
enum CenterTab {
    Modeling,
    Result,
    /// A named report tab. The String is the report name.
    Report(String),
}

impl CenterTab {
    fn title(&self) -> String {
        match self {
            Self::Modeling => "Modeling".into(),
            Self::Result => "Result".into(),
            Self::Report(name) => name.clone(),
        }
    }
}

struct CenterTabViewer<'a> {
    project: &'a Project,
    viewport: &'a mut SceneViewport,
    engine: &'a GeometryEngine,
    geometry_generation: u64,
    report_panels: &'a mut HashMap<String, ReportPanel>,
}

impl TabViewer for CenterTabViewer<'_> {
    type Tab = CenterTab;

    fn title(&mut self, tab: &mut Self::Tab) -> egui::WidgetText {
        tab.title().into()
    }

    fn ui(&mut self, ui: &mut egui::Ui, tab: &mut Self::Tab) {
        match tab {
            CenterTab::Modeling => {
                self.viewport.ui(ui, self.engine.all_breps(), self.geometry_generation);
            }
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
            CenterTab::Report(name) => {
                if let Some(panel) = self.report_panels.get_mut(name) {
                    panel.ui(ui);
                } else {
                    ui.label(format!("Report '{}' not found", name));
                }
            }
        }
    }

    fn closeable(&mut self, tab: &mut Self::Tab) -> bool {
        matches!(tab, CenterTab::Report(_))
    }

    fn on_close(&mut self, tab: &mut Self::Tab) -> egui_dock::widgets::tab_viewer::OnCloseResponse {
        if let CenterTab::Report(name) = tab {
            self.report_panels.remove(name);
        }
        egui_dock::widgets::tab_viewer::OnCloseResponse::Close
    }
}

// ---------------------------------------------------------------------------
// App
// ---------------------------------------------------------------------------

pub struct App {
    project: Project,
    backend: Box<dyn Backend>,
    edition: Edition,
    mode: RunMode,
    viewport: SceneViewport,
    engine: GeometryEngine,
    geometry_generation: u64,
    dock_state: DockState<CenterTab>,
    ribbon_state: RibbonState,
    ribbon_tabs: Vec<RibbonTab>,
    current_file: Option<PathBuf>,
    unsaved_changes: bool,
    status_text: String,
    log_text: String,
    messages: Vec<MessageEntry>,
    file_dialog_rx: mpsc::Receiver<FileDialogResult>,
    file_dialog_tx: mpsc::Sender<FileDialogResult>,

    // Report system
    report_panels: HashMap<String, ReportPanel>,
    result_store: Option<ResultDataStore>,

    // Geometry modeling UI state
    selected_object: Option<String>,
    design_variables: HashMap<String, Variable>,
    draw_dialog: Option<DrawDialog>,
    variable_edit_buffers: HashMap<String, String>,
    param_edit_buffers: HashMap<String, String>,

    // Layout state
    show_project_manager: bool,
    show_message_manager: bool,
    left_panel_active_tab: LeftPanelTab,
    bottom_tab: BottomTab,
}

// Public accessors for testing and external inspection
impl App {
    pub fn project(&self) -> &Project {
        &self.project
    }

    pub fn edition(&self) -> Edition {
        self.edition
    }

    pub fn current_file(&self) -> Option<&PathBuf> {
        self.current_file.as_ref()
    }

    pub fn unsaved_changes(&self) -> bool {
        self.unsaved_changes
    }

    pub fn status_text(&self) -> &str {
        &self.status_text
    }

    pub fn log_text(&self) -> &str {
        &self.log_text
    }

    pub fn ribbon_state(&self) -> &RibbonState {
        &self.ribbon_state
    }

    pub fn ribbon_state_mut(&mut self) -> &mut RibbonState {
        &mut self.ribbon_state
    }

    /// Dispatch a ribbon action programmatically (used by UI and tests).
    pub fn dispatch_action(&mut self, action: RibbonAction) {
        self.on_ribbon_action(action);
    }

    /// Save the current project to a specific path (native) or via backend (WASM).
    pub fn save_to(&mut self, path: &std::path::Path) {
        #[cfg(not(target_arch = "wasm32"))]
        {
            match save_project_to_file(&self.project, path) {
                Ok(()) => {
                    self.current_file = Some(path.to_path_buf());
                    self.unsaved_changes = false;
                    self.status_text = format!("Saved: {}", path.display());
                    self.log_text
                        .push_str(&format!("\n[file] saved to {}", path.display()));
                }
                Err(e) => {
                    self.status_text = format!("Save failed: {e}");
                    self.messages
                        .push(MessageEntry::error(format!("Save failed: {e}")));
                }
            }
        }
        #[cfg(target_arch = "wasm32")]
        {
            let _ = path; // unused on WASM
            self.save_to_backend();
        }
    }

    /// Save the current project through the backend (used on WASM / LocalFirst).
    pub fn save_to_backend(&mut self) {
        match self.backend.save_project(self.project.clone()) {
            Ok(()) => {
                self.unsaved_changes = false;
                self.status_text = "Saved to OPFS".into();
                self.log_text.push_str("\n[file] saved to OPFS");
            }
            Err(e) => {
                self.status_text = format!("Save failed: {e}");
                self.messages
                    .push(MessageEntry::error(format!("Save failed: {e}")));
            }
        }
    }

    /// Load a project from a specific path (native) or via backend (WASM).
    pub fn open_from(&mut self, path: &std::path::Path) {
        #[cfg(not(target_arch = "wasm32"))]
        {
            match load_project_from_file(path) {
                Ok(project) => {
                    self.project = project;
                    self.current_file = Some(path.to_path_buf());
                    self.unsaved_changes = false;
                    self.status_text = format!("Opened: {}", path.display());
                    self.log_text
                        .push_str(&format!("\n[file] opened {}", path.display()));
                    // Initialize result store from project path
                    self.result_store = Some(ResultDataStore::from_project_path(path));
                    // Clear existing report panels on new project load
                    self.report_panels.clear();
                }
                Err(e) => {
                    self.status_text = format!("Open failed: {e}");
                    self.messages
                        .push(MessageEntry::error(format!("Open failed: {e}")));
                    self.log_text
                        .push_str(&format!("\n[file] open failed: {e}"));
                }
            }
        }
        #[cfg(target_arch = "wasm32")]
        {
            let _ = path;
            // On WASM, projects are loaded through the backend
            self.status_text = "Use project list to open projects on web".into();
        }
    }

    /// Load a project by id through the backend (WASM / LocalFirst).
    pub fn load_from_backend(&mut self, id: &str) {
        match self.backend.load_project(id) {
            Ok(project) => {
                self.project = project;
                self.unsaved_changes = false;
                self.status_text = format!("Opened project: {id}");
                self.log_text
                    .push_str(&format!("\n[file] opened project {id}"));
            }
            Err(_) => {
                // On WASM with LocalFirst, this means the load is pending (async).
                // The project will arrive via poll().
                self.status_text = format!("Loading project: {id}...");
            }
        }
    }
}

impl App {
    fn mode_label(mode: RunMode) -> &'static str {
        match mode {
            RunMode::Standalone => "standalone",
            RunMode::Cloud => "cloud",
            RunMode::LocalFirst => "local-first",
        }
    }

    pub fn new(mode: RunMode, edition: Edition, cc: &eframe::CreationContext<'_>) -> Self {
        let mut app = Self::new_headless(mode, edition);
        if let Some(rs) = &cc.wgpu_render_state {
            app.viewport.init_renderer(&rs.device, rs.target_format, &mut rs.renderer.write().callback_resources);
        }
        app
    }

    /// Create an App without GPU renderer (for tests and WASM fallback).
    pub fn new_headless(mode: RunMode, edition: Edition) -> Self {
        let project = Project::default();

        // Center dock: Modeling and Result as sibling tabs (not split)
        let dock_state = DockState::new(vec![CenterTab::Modeling, CenterTab::Result]);

        let (tx, rx) = mpsc::channel();

        let mut app = Self {
            project,
            backend: default_backend(mode),
            edition,
            mode,
            viewport: SceneViewport::default(),
            engine: GeometryEngine::new(),
            geometry_generation: 0,
            dock_state,
            ribbon_state: RibbonState::default(),
            ribbon_tabs: build_default_tabs(edition, mode == RunMode::LocalFirst),
            current_file: None,
            unsaved_changes: false,
            status_text: format!("Ready ({}, {})", Self::mode_label(mode), edition.display_name()),
            log_text: format!(
                "[boot] EmStudio shell started (mode={}, edition={})",
                Self::mode_label(mode),
                edition.display_name(),
            ),
            messages: vec![MessageEntry::info("EmStudio started.")],
            file_dialog_rx: rx,
            file_dialog_tx: tx,

            // Report system
            report_panels: HashMap::new(),
            result_store: None,

            // Geometry modeling UI state
            selected_object: None,
            design_variables: HashMap::new(),
            draw_dialog: None,
            variable_edit_buffers: HashMap::new(),
            param_edit_buffers: HashMap::new(),

            // Layout defaults
            show_project_manager: true,
            show_message_manager: true,
            left_panel_active_tab: LeftPanelTab::ProjectManager,
            bottom_tab: BottomTab::Messages,
        };

        if mode == RunMode::LocalFirst {
            app.load_from_backend("default");
            app.status_text = "Loading default web project...".into();
        }

        app
    }

    pub fn new_default(cc: &eframe::CreationContext<'_>) -> Self {
        Self::new(RunMode::Standalone, Edition::Professional, cc)
    }

    pub fn new_default_headless() -> Self {
        Self::new_headless(RunMode::Standalone, Edition::Professional)
    }

    // -----------------------------------------------------------------------
    // File title helper
    // -----------------------------------------------------------------------

    pub fn file_display_name(&self) -> String {
        let name = match &self.current_file {
            Some(p) => p
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_else(|| "Unknown".into()),
            None => "Untitled".into(),
        };
        if self.unsaved_changes {
            format!("{}*", name)
        } else {
            name
        }
    }

    // -----------------------------------------------------------------------
    // Async file dialog helpers
    // -----------------------------------------------------------------------

    fn spawn_open_dialog(&self, ctx: &egui::Context) {
        let tx = self.file_dialog_tx.clone();
        let ctx = ctx.clone();
        spawn_future(async move {
            let file = rfd::AsyncFileDialog::new()
                .add_filter("EmStudio Project", &["emsp"])
                .pick_file()
                .await;
            if let Some(handle) = file {
                #[cfg(not(target_arch = "wasm32"))]
                let path = handle.path().to_path_buf();
                #[cfg(target_arch = "wasm32")]
                let path = PathBuf::from(handle.file_name());

                let _ = tx.send(FileDialogResult::OpenFile(path));
                ctx.request_repaint();
            }
        });
    }

    fn spawn_save_dialog(&self, ctx: &egui::Context) {
        let tx = self.file_dialog_tx.clone();
        let ctx = ctx.clone();
        spawn_future(async move {
            let file = rfd::AsyncFileDialog::new()
                .add_filter("EmStudio Project", &["emsp"])
                .set_file_name("project.emsp")
                .save_file()
                .await;
            if let Some(handle) = file {
                #[cfg(not(target_arch = "wasm32"))]
                let path = handle.path().to_path_buf();
                #[cfg(target_arch = "wasm32")]
                let path = PathBuf::from(handle.file_name());

                let _ = tx.send(FileDialogResult::SaveFile(path));
                ctx.request_repaint();
            }
        });
    }

    fn spawn_import_hfss_dialog(&self, ctx: &egui::Context) {
        let tx = self.file_dialog_tx.clone();
        let ctx = ctx.clone();
        spawn_future(async move {
            let file = rfd::AsyncFileDialog::new()
                .add_filter("Ansys Project", &["aedt", "py", "txt"])
                .pick_file()
                .await;
            if let Some(handle) = file {
                #[cfg(not(target_arch = "wasm32"))]
                let path = handle.path().to_path_buf();
                #[cfg(target_arch = "wasm32")]
                let path = PathBuf::from(handle.file_name());

                let _ = tx.send(FileDialogResult::ImportHfssAedt(path));
                ctx.request_repaint();
            }
        });
    }

    fn spawn_import_q3d_dialog(&self, ctx: &egui::Context) {
        let tx = self.file_dialog_tx.clone();
        let ctx = ctx.clone();
        spawn_future(async move {
            let file = rfd::AsyncFileDialog::new()
                .add_filter("Ansys Project", &["aedt", "py", "txt"])
                .pick_file()
                .await;
            if let Some(handle) = file {
                #[cfg(not(target_arch = "wasm32"))]
                let path = handle.path().to_path_buf();
                #[cfg(target_arch = "wasm32")]
                let path = PathBuf::from(handle.file_name());

                let _ = tx.send(FileDialogResult::ImportQ3dAedt(path));
                ctx.request_repaint();
            }
        });
    }

    fn spawn_export_hfss_dialog(&self, ctx: &egui::Context) {
        let tx = self.file_dialog_tx.clone();
        let ctx = ctx.clone();
        spawn_future(async move {
            let file = rfd::AsyncFileDialog::new()
                .add_filter("Python Script", &["py"])
                .set_file_name("hfss_export.py")
                .save_file()
                .await;
            if let Some(handle) = file {
                #[cfg(not(target_arch = "wasm32"))]
                let path = handle.path().to_path_buf();
                #[cfg(target_arch = "wasm32")]
                let path = PathBuf::from(handle.file_name());

                let _ = tx.send(FileDialogResult::ExportHfssPy(path));
                ctx.request_repaint();
            }
        });
    }

    fn spawn_export_q3d_dialog(&self, ctx: &egui::Context) {
        let tx = self.file_dialog_tx.clone();
        let ctx = ctx.clone();
        spawn_future(async move {
            let file = rfd::AsyncFileDialog::new()
                .add_filter("Python Script", &["py"])
                .set_file_name("q3d_export.py")
                .save_file()
                .await;
            if let Some(handle) = file {
                #[cfg(not(target_arch = "wasm32"))]
                let path = handle.path().to_path_buf();
                #[cfg(target_arch = "wasm32")]
                let path = PathBuf::from(handle.file_name());

                let _ = tx.send(FileDialogResult::ExportQ3dPy(path));
                ctx.request_repaint();
            }
        });
    }

    /// Poll for completed file dialog results. Called each frame in update().
    fn poll_file_dialogs(&mut self) {
        while let Ok(result) = self.file_dialog_rx.try_recv() {
            match result {
                FileDialogResult::OpenFile(path) => {
                    self.open_from(&path);
                }
                FileDialogResult::SaveFile(path) => {
                    self.save_to(&path);
                }
                FileDialogResult::ImportHfssAedt(path) => {
                    #[cfg(not(target_arch = "wasm32"))]
                    {
                        match import_hfss_project_from_file(&path) {
                            Ok(project) => {
                                self.project = project;
                                self.current_file = None;
                                self.unsaved_changes = true;
                                self.status_text = format!("Imported HFSS: {}", path.display());
                                self.messages.push(MessageEntry::info(format!(
                                    "Imported HFSS project from {}",
                                    path.display()
                                )));
                            }
                            Err(e) => {
                                self.status_text = format!("HFSS import failed: {e}");
                                self.messages
                                    .push(MessageEntry::error(format!("HFSS import failed: {e}")));
                            }
                        }
                    }
                }
                FileDialogResult::ImportQ3dAedt(path) => {
                    #[cfg(not(target_arch = "wasm32"))]
                    {
                        match import_q3d_project_from_file(&path) {
                            Ok(project) => {
                                self.project = project;
                                self.current_file = None;
                                self.unsaved_changes = true;
                                self.status_text = format!("Imported Q3D: {}", path.display());
                                self.messages.push(MessageEntry::info(format!(
                                    "Imported Q3D project from {}",
                                    path.display()
                                )));
                            }
                            Err(e) => {
                                self.status_text = format!("Q3D import failed: {e}");
                                self.messages
                                    .push(MessageEntry::error(format!("Q3D import failed: {e}")));
                            }
                        }
                    }
                }
                FileDialogResult::ExportHfssPy(path) => {
                    #[cfg(not(target_arch = "wasm32"))]
                    {
                        match export_hfss_project_to_file(&self.project, &path) {
                            Ok(()) => {
                                self.status_text = format!("Exported HFSS script: {}", path.display());
                                self.messages.push(MessageEntry::info(format!(
                                    "Exported HFSS script to {}",
                                    path.display()
                                )));
                            }
                            Err(e) => {
                                self.status_text = format!("HFSS export failed: {e}");
                                self.messages
                                    .push(MessageEntry::error(format!("HFSS export failed: {e}")));
                            }
                        }
                    }
                }
                FileDialogResult::ExportQ3dPy(path) => {
                    #[cfg(not(target_arch = "wasm32"))]
                    {
                        match export_q3d_project_to_file(&self.project, &path) {
                            Ok(()) => {
                                self.status_text = format!("Exported Q3D script: {}", path.display());
                                self.messages.push(MessageEntry::info(format!(
                                    "Exported Q3D script to {}",
                                    path.display()
                                )));
                            }
                            Err(e) => {
                                self.status_text = format!("Q3D export failed: {e}");
                                self.messages
                                    .push(MessageEntry::error(format!("Q3D export failed: {e}")));
                            }
                        }
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Core file operations
    // -----------------------------------------------------------------------

    fn do_save(&mut self, ctx: &egui::Context) {
        if let Some(path) = self.current_file.clone() {
            self.save_to(&path);
        } else {
            self.spawn_save_dialog(ctx);
        }
    }

    fn do_save_as(&self, ctx: &egui::Context) {
        self.spawn_save_dialog(ctx);
    }

    fn do_open(&self, ctx: &egui::Context) {
        self.spawn_open_dialog(ctx);
    }

    fn do_import_hfss(&self, ctx: &egui::Context) {
        self.spawn_import_hfss_dialog(ctx);
    }

    fn do_import_q3d(&self, ctx: &egui::Context) {
        self.spawn_import_q3d_dialog(ctx);
    }

    fn do_export_hfss(&self, ctx: &egui::Context) {
        self.spawn_export_hfss_dialog(ctx);
    }

    fn do_export_q3d(&self, ctx: &egui::Context) {
        self.spawn_export_q3d_dialog(ctx);
    }

    fn do_new(&mut self) {
        self.project = Project::default();
        self.current_file = None;
        self.unsaved_changes = false;
        self.status_text = "New project created".into();
        self.log_text.push_str("\n[file] new project");
        self.messages.push(MessageEntry::info("New project created."));
    }

    // -----------------------------------------------------------------------
    // Ribbon action dispatch
    // -----------------------------------------------------------------------

    fn on_ribbon_action(&mut self, action: RibbonAction) {
        // File operations that need ctx are handled in on_ribbon_action_with_ctx
        match action {
            RibbonAction::NewProject => {
                self.do_new();
            }
            RibbonAction::CloseProject => {
                self.do_new();
            }

            // -- Solve --
            RibbonAction::Solve => {
                self.project.status = SimulationStatus::Solving;
                match self.backend.solve(&self.project) {
                    Ok(result) => {
                        self.project.status = SimulationStatus::Finished;
                        self.project.last_result = Some(result);
                        self.status_text = "Solve completed".into();
                        self.log_text.push_str("\n[action] solve");
                        self.messages
                            .push(MessageEntry::info("Solve completed successfully."));
                        self.unsaved_changes = true;
                    }
                    Err(err) => {
                        self.project.status = SimulationStatus::Failed;
                        self.status_text = format!("Solve failed: {err}");
                        self.messages
                            .push(MessageEntry::error(format!("Solve failed: {err}")));
                    }
                }
            }
            RibbonAction::SolveAll => {
                if !self.edition.allows_solve_all() {
                    self.status_text = "Solve All requires Professional edition or higher".into();
                    self.messages.push(MessageEntry::warning(
                        "Solve All is not available in Basic edition.",
                    ));
                    return;
                }
                self.project.status = SimulationStatus::Solving;
                match self.backend.solve(&self.project) {
                    Ok(result) => {
                        self.project.status = SimulationStatus::Finished;
                        self.project.last_result = Some(result);
                        self.status_text = "Solve All completed".into();
                        self.log_text.push_str("\n[action] solve all");
                        self.messages
                            .push(MessageEntry::info("Solve All completed successfully."));
                        self.unsaved_changes = true;
                    }
                    Err(err) => {
                        self.project.status = SimulationStatus::Failed;
                        self.status_text = format!("Solve All failed: {err}");
                        self.messages
                            .push(MessageEntry::error(format!("Solve All failed: {err}")));
                    }
                }
            }

            // -- Toggle actions (state already toggled in ribbon) --
            RibbonAction::ToggleGrid => {
                let on = self
                    .ribbon_state
                    .toggles
                    .get("grid")
                    .copied()
                    .unwrap_or(false);
                self.status_text = format!("Grid: {}", if on { "ON" } else { "OFF" });
            }
            RibbonAction::ToggleRuler => {
                let on = self
                    .ribbon_state
                    .toggles
                    .get("ruler")
                    .copied()
                    .unwrap_or(false);
                self.status_text = format!("Ruler: {}", if on { "ON" } else { "OFF" });
            }
            RibbonAction::ToggleCoordSystem => {
                let on = self
                    .ribbon_state
                    .toggles
                    .get("coord_system")
                    .copied()
                    .unwrap_or(false);
                self.status_text = format!("Coord System: {}", if on { "ON" } else { "OFF" });
            }
            RibbonAction::RenderShaded | RibbonAction::RenderWireframe => {
                self.status_text = format!("{:?} mode selected", action);
            }

            // -- Layout toggles --
            RibbonAction::ToggleProjectManager => {
                self.show_project_manager = !self.show_project_manager;
            }
            RibbonAction::ToggleMessageManager => {
                self.show_message_manager = !self.show_message_manager;
            }

            // File dialog actions are no-ops here — handled below
            RibbonAction::SaveProject | RibbonAction::SaveAs | RibbonAction::OpenProject => {}

            // -- Reports --
            RibbonAction::CreateReport => {
                self.create_default_report();
            }

            // -- Draw commands --
            RibbonAction::DrawBox => {
                self.draw_dialog = Some(DrawDialog::new_box());
            }
            RibbonAction::DrawCylinder => {
                self.draw_dialog = Some(DrawDialog::new_cylinder());
            }
            RibbonAction::DrawSphere => {
                self.draw_dialog = Some(DrawDialog::new_sphere());
            }
            RibbonAction::DrawCone => {
                self.draw_dialog = Some(DrawDialog::new_cone());
            }
            RibbonAction::DrawTorus => {
                self.draw_dialog = Some(DrawDialog::new_torus());
            }

            // -- All other actions: mark dirty + stub --
            other => {
                let name = format!("{:?}", other);
                self.status_text = format!("{} (not implemented yet)", name);
                self.log_text
                    .push_str(&format!("\n[action] {} (stub)", name));
                self.unsaved_changes = true;
            }
        }
    }

    /// Handle actions that require egui::Context (for spawning async dialogs).
    fn on_ribbon_action_with_ctx(&mut self, action: RibbonAction, ctx: &egui::Context) {
        match action {
            RibbonAction::SaveProject => self.do_save(ctx),
            RibbonAction::SaveAs => self.do_save_as(ctx),
            RibbonAction::OpenProject => self.do_open(ctx),
            RibbonAction::ImportHfssAedt => self.do_import_hfss(ctx),
            RibbonAction::ImportQ3dAedt => self.do_import_q3d(ctx),
            RibbonAction::ExportHfssPyAedt => self.do_export_hfss(ctx),
            RibbonAction::ExportQ3dPyAedt => self.do_export_q3d(ctx),
            _ => self.on_ribbon_action(action),
        }
    }

    // -----------------------------------------------------------------------
    // Report creation
    // -----------------------------------------------------------------------

    fn create_default_report(&mut self) {
        use emstudio_domain::report::*;

        // Generate unique report name
        let report_num = self.report_panels.len() + 1;
        let name = format!("S Parameter Plot {}", report_num);

        let report = Report {
            name: name.clone(),
            category: ReportCategory::SParameter,
            chart_type: ChartType::Rectangular,
            solution: "Setup1".into(),
            domain: ReportDomain {
                domain_type: "Frequency".into(),
                primary_sweep: "Freq".into(),
                fixed_values: None,
            },
            traces: vec![ReportTrace {
                name: "dB(S(1,1))".into(),
                expression: "dB(S(1,1))".into(),
                style: Some(TraceStyle {
                    color: [0, 0, 255],
                    line_width: 2,
                    line_style: "Solid".into(),
                }),
                parametric_values: None,
                fixed_values: None,
            }],
            x_axis: Some(AxisConfig {
                label: "Frequency".into(),
                unit: "GHz".into(),
                min: None,
                max: None,
                auto_range: Some(true),
                scale: None,
            }),
            y_axis: Some(AxisConfig {
                label: "Magnitude".into(),
                unit: "dB".into(),
                min: None,
                max: None,
                auto_range: Some(true),
                scale: None,
            }),
            markers: Vec::new(),
            limit_lines: Vec::new(),
            far_field_setup: None,
            matrix_type: None,
            display_options: None,
        };

        let mut panel = ReportPanel::new(report);

        // If we have a result store with S-parameter data, try to load it
        self.try_load_report_traces(&mut panel);

        self.report_panels.insert(name.clone(), panel);

        // Add tab to dock
        self.dock_state
            .main_surface_mut()
            .push_to_first_leaf(CenterTab::Report(name.clone()));

        self.status_text = format!("Created report: {}", name);
        self.log_text
            .push_str(&format!("\n[report] created: {}", name));
        self.messages
            .push(MessageEntry::info(format!("Report '{}' created", name)));
    }

    fn try_load_report_traces(&self, panel: &mut ReportPanel) {
        use emstudio_domain::quantity_expr::{EvalContext, QuantityExpression};

        // Try to load from result store if available
        if let Some(store) = &self.result_store {
            let s_param_path = store.base_path().join("s_parameters.json");
            if s_param_path.exists() {
                if let Ok(data) =
                    emstudio_domain::result_store::load_s_parameters_from_file(&s_param_path)
                {
                    // Collect trace info first to avoid borrow conflict
                    let trace_info: Vec<_> = panel
                        .report
                        .traces
                        .iter()
                        .map(|t| (t.name.clone(), t.expression.clone(), t.style.clone()))
                        .collect();
                    for (name, expression, style) in &trace_info {
                        if let Ok(expr) = QuantityExpression::parse(expression) {
                            let ctx = EvalContext::SParameter { data: &data };
                            if let Ok(points) = expr.evaluate(&ctx) {
                                panel.set_trace_data(
                                    name.clone(),
                                    points,
                                    style.as_ref(),
                                );
                            }
                        }
                    }
                }
            }
        }

        // If no data was loaded, provide demo data so the chart isn't empty
        if panel.trace_data.is_empty() {
            self.load_demo_data(panel);
        }
    }

    fn load_demo_data(&self, panel: &mut ReportPanel) {
        // Collect trace info first to avoid borrow conflict
        let trace_info: Vec<_> = panel
            .report
            .traces
            .iter()
            .map(|t| (t.name.clone(), t.style.clone()))
            .collect();
        for (name, style) in &trace_info {
            // Generate synthetic S11 data (typical patch antenna response)
            let mut points = Vec::with_capacity(201);
            for i in 0..=200 {
                let freq = 1.0 + i as f64 * 0.02; // 1-5 GHz
                let f0 = 2.45; // Resonant frequency
                let bw = 0.15; // Bandwidth
                let delta = (freq - f0) / bw;
                let s11_db = -5.0 - 25.0 / (1.0 + delta * delta);
                points.push([freq, s11_db]);
            }
            panel.set_trace_data(name.clone(), points, style.as_ref());
        }
    }
}

// ---------------------------------------------------------------------------
// eframe::App — AEDT Standard Layout
// ---------------------------------------------------------------------------

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Poll async file dialog results
        self.poll_file_dialogs();

        // Poll backend for async results (WASM worker responses)
        self.backend.poll();

        // Check for async solve results
        if let Some(result) = self.backend.take_solve_result() {
            self.project.status = if result.converged {
                SimulationStatus::Finished
            } else {
                SimulationStatus::Failed
            };
            self.project.last_result = Some(result);
            self.unsaved_changes = true;
            self.status_text = "Solve completed".into();
            self.log_text.push_str("\n[solver] async solve completed");
            self.messages
                .push(MessageEntry::info("Solve completed (async)."));
        }

        // Check for async project load results (WASM worker responses)
        if let Some(project) = self.backend.take_loaded_project() {
            let project_id = project.id.clone();
            self.project = project;
            self.unsaved_changes = false;
            self.status_text = format!("Opened project: {project_id}");
            self.log_text
                .push_str(&format!("\n[file] opened project {project_id} (async)"));
            self.messages
                .push(MessageEntry::info(format!("Project loaded: {project_id}")));
            self.report_panels.clear();
        }

        // =================================================================
        // TOP PANELS (menu bar → QAT → ribbon)
        // =================================================================

        // 1. Menu Bar
        egui::TopBottomPanel::top("menu_bar")
            .exact_height(24.0)
            .show(ctx, |ui| {
                if let Some(action) = menu_bar::show_menu_bar(ui, self.edition, self.mode == RunMode::LocalFirst) {
                    self.on_ribbon_action_with_ctx(action, ctx);
                }
            });

        // 2. Quick Access Toolbar
        egui::TopBottomPanel::top("qat")
            .exact_height(26.0)
            .show(ctx, |ui| {
                if let Some(action) = qat::show_qat(ui) {
                    self.on_ribbon_action_with_ctx(action, ctx);
                }
            });

        // 3. Ribbon
        egui::TopBottomPanel::top("ribbon").show(ctx, |ui| {
            if let Some(action) = show_ribbon(ui, &mut self.ribbon_state, &self.ribbon_tabs) {
                self.on_ribbon_action_with_ctx(action, ctx);
            }
        });

        // =================================================================
        // BOTTOM PANELS (status bar, then message manager)
        // =================================================================

        // 4. Status Bar (always visible, bottom-most)
        egui::TopBottomPanel::bottom("status_bar")
            .exact_height(24.0)
            .show(ctx, |ui| {
                let state = StatusBarState {
                    file_name: self.file_display_name(),
                    unsaved: self.unsaved_changes,
                    status_text: self.status_text.clone(),
                    coordinates: None,
                    units: "mm".into(),
                    message_manager_visible: self.show_message_manager,
                };
                let response = status_bar::show_status_bar(ui, &state);
                if response.toggle_message_manager {
                    self.show_message_manager = !self.show_message_manager;
                }
            });

        // 5. Bottom Dock: Message Manager / Log (conditional)
        if self.show_message_manager {
            egui::TopBottomPanel::bottom("bottom_dock")
                .resizable(true)
                .default_height(150.0)
                .min_height(80.0)
                .show(ctx, |ui| {
                    message_manager::show_message_manager(
                        ui,
                        &mut self.bottom_tab,
                        &self.messages,
                        &self.log_text,
                    );
                });
        }

        // =================================================================
        // SIDE PANELS
        // =================================================================

        // 6. Left Panel: Project Manager + Properties (tabbed)
        if self.show_project_manager {
            egui::SidePanel::left("left_dock")
                .resizable(true)
                .default_width(240.0)
                .min_width(180.0)
                .show(ctx, |ui| {
                    let mut dock_ctx = LeftDockContext {
                        project: &self.project,
                        selected_object: self.selected_object.as_deref(),
                        selected_operation: None, // TODO: lookup from geometry operations
                        design_variables: &self.design_variables,
                        variable_edit_buffers: &mut self.variable_edit_buffers,
                        param_edit_buffers: &mut self.param_edit_buffers,
                    };
                    let dock_response = dock::left_dock_panel(
                        ui,
                        &mut self.left_panel_active_tab,
                        &mut dock_ctx,
                    );

                    // Handle tree selection
                    if let Some(obj_name) = dock_response.tree.selected_object {
                        self.selected_object = Some(obj_name);
                    }

                    // Handle variable add
                    if dock_response.tree.add_variable {
                        let name = format!("var{}", self.design_variables.len() + 1);
                        self.design_variables.insert(
                            name,
                            Variable {
                                value: Some("0".into()),
                                expression: None,
                                description: String::new(),
                                unit_type: None,
                            },
                        );
                        self.unsaved_changes = true;
                    }

                    // Handle variable edit
                    if let Some((name, new_value)) = dock_response.tree.variable_edited {
                        if let Some(var) = self.design_variables.get_mut(&name) {
                            var.value = Some(new_value);
                            self.unsaved_changes = true;
                        }
                    }
                });
        }

        // Draw dialog
        if let Some(dialog) = &mut self.draw_dialog {
            let var_names: Vec<String> = self.design_variables.keys().cloned().collect();
            if let Some(result) = dialog.show(ctx, &var_names) {
                // Assign step number
                let mut op = result.operation;
                let next_step = self
                    .project
                    .model
                    .objects
                    .len() as u32
                    + 1;
                op.step = next_step;

                // Execute the operation
                let vars: HashMap<String, f64> = self
                    .design_variables
                    .iter()
                    .filter_map(|(k, v)| {
                        v.value
                            .as_ref()
                            .and_then(|val| emstudio_domain::expression::parse_value_with_unit(val).ok())
                            .map(|f| (k.clone(), f))
                    })
                    .collect();

                match self.engine.execute(&op, &vars) {
                    Ok(Some(obj)) => {
                        // Add object to project (legacy model)
                        self.project.model.objects.push(emstudio_domain::GeometryObject {
                            id: obj.id,
                            name: obj.name.clone(),
                            mesh_hint: "auto".into(),
                        });
                        self.geometry_generation += 1;
                        self.unsaved_changes = true;
                        self.status_text = format!("Created: {}", obj.name);
                        self.messages.push(MessageEntry::info(format!("Created object '{}'", obj.name)));
                    }
                    Ok(None) => {}
                    Err(e) => {
                        self.status_text = format!("Create failed: {e}");
                        self.messages.push(MessageEntry::error(format!("Create failed: {e}")));
                    }
                }

                self.draw_dialog = None;
            } else if !dialog.open {
                self.draw_dialog = None;
            }
        }

        // =================================================================
        // CENTRAL PANEL
        // =================================================================

        // 7. Central Workspace: 3D Modeler + Result + Report dock tabs
        egui::CentralPanel::default().show(ctx, |ui| {
            let mut viewer = CenterTabViewer {
                project: &self.project,
                viewport: &mut self.viewport,
                engine: &self.engine,
                geometry_generation: self.geometry_generation,
                report_panels: &mut self.report_panels,
            };
            DockArea::new(&mut self.dock_state)
                .show_leaf_collapse_buttons(false)
                .show_close_buttons(true)
                .show_leaf_close_all_buttons(false)
                .show_inside(ui, &mut viewer);
                
        });
    }
}
