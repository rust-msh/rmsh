use std::path::PathBuf;
use std::sync::mpsc;

use eframe::egui;
use egui_dock::{DockArea, DockState, TabViewer};

use emstudio_components::dock;
use emstudio_components::message_manager::{self, BottomTab, MessageEntry};
use emstudio_components::ribbon::{
    RibbonAction, RibbonState, RibbonTab, build_default_tabs, show_ribbon,
};
use emstudio_components::status_bar::{self, StatusBarState};
use emstudio_components::{LeftPanelTab, menu_bar, qat};
use emstudio_domain::{Project, SimulationStatus};
use emstudio_domain::geometry_engine::GeometryEngine;
use emstudio_infra::{Backend, RunMode, default_backend};
#[cfg(not(target_arch = "wasm32"))]
use emstudio_infra::{load_project_from_file, save_project_to_file};
use emstudio_render::SceneViewport;

// ---------------------------------------------------------------------------
// Async file dialog result
// ---------------------------------------------------------------------------

enum FileDialogResult {
    OpenFile(PathBuf),
    SaveFile(PathBuf),
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
}

impl CenterTab {
    fn title(&self) -> &'static str {
        match self {
            Self::Modeling => "Modeling",
            Self::Result => "Result",
        }
    }
}

struct CenterTabViewer<'a> {
    project: &'a Project,
    viewport: &'a mut SceneViewport,
    engine: &'a GeometryEngine,
    geometry_generation: u64,
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
        }
    }
}

// ---------------------------------------------------------------------------
// App
// ---------------------------------------------------------------------------

pub struct App {
    project: Project,
    backend: Box<dyn Backend>,
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

    pub fn new(mode: RunMode, cc: &eframe::CreationContext<'_>) -> Self {
        let mut app = Self::new_headless(mode);
        if let Some(rs) = &cc.wgpu_render_state {
            app.viewport.init_renderer(&rs.device, rs.target_format, &mut rs.renderer.write().callback_resources);
        }
        app
    }

    /// Create an App without GPU renderer (for tests and WASM fallback).
    pub fn new_headless(mode: RunMode) -> Self {
        let project = Project::default();

        // Center dock: Modeling and Result as sibling tabs (not split)
        let dock_state = DockState::new(vec![CenterTab::Modeling, CenterTab::Result]);

        let (tx, rx) = mpsc::channel();

        Self {
            project,
            backend: default_backend(mode),
            viewport: SceneViewport::default(),
            engine: GeometryEngine::new(),
            geometry_generation: 0,
            dock_state,
            ribbon_state: RibbonState::default(),
            ribbon_tabs: build_default_tabs(),
            current_file: None,
            unsaved_changes: false,
            status_text: format!("Ready ({})", Self::mode_label(mode)),
            log_text: format!(
                "[boot] EmStudio shell started (mode={})",
                Self::mode_label(mode)
            ),
            messages: vec![MessageEntry::info("EmStudio started.")],
            file_dialog_rx: rx,
            file_dialog_tx: tx,

            // Layout defaults
            show_project_manager: true,
            show_message_manager: true,
            left_panel_active_tab: LeftPanelTab::ProjectManager,
            bottom_tab: BottomTab::Messages,
        }
    }

    pub fn new_default(cc: &eframe::CreationContext<'_>) -> Self {
        Self::new(RunMode::Standalone, cc)
    }

    pub fn new_default_headless() -> Self {
        Self::new_headless(RunMode::Standalone)
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
            RibbonAction::Solve | RibbonAction::SolveAll => {
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
            _ => self.on_ribbon_action(action),
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

        // =================================================================
        // TOP PANELS (menu bar → QAT → ribbon)
        // =================================================================

        // 1. Menu Bar
        egui::TopBottomPanel::top("menu_bar")
            .exact_height(24.0)
            .show(ctx, |ui| {
                if let Some(action) = menu_bar::show_menu_bar(ui) {
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
                    dock::left_dock_panel(ui, &self.project, &mut self.left_panel_active_tab);
                });
        }

        // =================================================================
        // CENTRAL PANEL
        // =================================================================

        // 7. Central Workspace: 3D Modeler + Result dock tabs
        egui::CentralPanel::default().show(ctx, |ui| {
            let mut viewer = CenterTabViewer {
                project: &self.project,
                viewport: &mut self.viewport,
                engine: &self.engine,
                geometry_generation: self.geometry_generation,
            };
            DockArea::new(&mut self.dock_state).show_inside(ui, &mut viewer);
        });
    }
}
