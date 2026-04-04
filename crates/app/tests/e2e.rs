//! End-to-end integration tests for the App.
//!
//! These tests exercise the full app lifecycle: project creation, modification,
//! save/load roundtrip, dirty flag tracking, and ribbon action dispatch — without
//! launching a GUI window.

use std::path::Path;

use emstudio_app::App;
use emstudio_components::ribbon::RibbonAction;
use emstudio_domain::SimulationStatus;
use emstudio_infra::RunMode;

fn new_app() -> App {
    App::new_headless(RunMode::Standalone)
}

// ---------------------------------------------------------------------------
// Initial state
// ---------------------------------------------------------------------------

#[test]
fn app_starts_with_default_project() {
    let app = new_app();
    assert_eq!(app.project().title, "New Project");
    assert_eq!(app.project().status, SimulationStatus::Idle);
    assert!(app.project().last_result.is_none());
    assert!(app.project().model.objects.is_empty());
}

#[test]
fn app_starts_with_no_file_and_clean_state() {
    let app = new_app();
    assert!(app.current_file().is_none());
    assert!(!app.unsaved_changes());
    assert!(app.status_text().contains("Ready"));
}

#[test]
fn initial_display_name_is_untitled() {
    let app = new_app();
    assert_eq!(app.file_display_name(), "Untitled");
}

// ---------------------------------------------------------------------------
// New project
// ---------------------------------------------------------------------------

#[test]
fn new_project_resets_state() {
    let mut app = new_app();

    // Dirty the project first
    app.dispatch_action(RibbonAction::DrawBox);
    assert!(app.unsaved_changes());

    // Save to a file so current_file is set
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("test.emsp");
    app.save_to(&path);
    assert!(app.current_file().is_some());

    // New project (bypassing dialog — no unsaved prompt in dispatch_action
    // because we just saved)
    app.dispatch_action(RibbonAction::NewProject);

    assert_eq!(app.project().title, "New Project");
    assert!(app.current_file().is_none());
    assert!(!app.unsaved_changes());
    assert!(app.status_text().contains("New project"));
}

// ---------------------------------------------------------------------------
// Save / Load roundtrip
// ---------------------------------------------------------------------------

#[test]
fn save_creates_file_and_clears_dirty() {
    let mut app = new_app();

    // Make a change
    app.dispatch_action(RibbonAction::Solve);
    assert!(app.unsaved_changes());

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("roundtrip.emsp");

    app.save_to(&path);

    assert!(path.exists());
    assert!(!app.unsaved_changes());
    assert_eq!(app.current_file().unwrap(), &path);
    assert!(app.status_text().contains("Saved"));
}

#[test]
fn save_then_load_preserves_project_data() {
    let mut app = new_app();

    // Solve to generate result data
    app.dispatch_action(RibbonAction::Solve);
    let original_title = app.project().title.clone();
    let original_converged = app
        .project()
        .last_result
        .as_ref()
        .map(|r| r.converged)
        .unwrap_or(false);

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("persist.emsp");
    app.save_to(&path);

    // Open in a fresh app
    let mut app2 = new_app();
    app2.open_from(&path);

    assert_eq!(app2.project().title, original_title);
    assert_eq!(
        app2.project()
            .last_result
            .as_ref()
            .map(|r| r.converged)
            .unwrap_or(false),
        original_converged
    );
    assert_eq!(app2.project().status, SimulationStatus::Finished);
    assert!(!app2.unsaved_changes());
    assert!(app2.status_text().contains("Opened"));
}

#[test]
fn open_nonexistent_file_shows_error() {
    let mut app = new_app();
    app.open_from(Path::new("/tmp/nonexistent_emstudio_file.emsp"));

    assert!(app.status_text().contains("Open failed"));
    assert!(app.current_file().is_none()); // should not change
}

#[test]
fn save_to_updates_current_file() {
    let mut app = new_app();
    assert!(app.current_file().is_none());

    let dir = tempfile::tempdir().unwrap();
    let path1 = dir.path().join("file1.emsp");
    let path2 = dir.path().join("file2.emsp");

    app.save_to(&path1);
    assert_eq!(app.current_file().unwrap(), &path1);
    assert_eq!(app.file_display_name(), "file1.emsp");

    app.save_to(&path2);
    assert_eq!(app.current_file().unwrap(), &path2);
    assert_eq!(app.file_display_name(), "file2.emsp");
}

// ---------------------------------------------------------------------------
// Dirty flag tracking
// ---------------------------------------------------------------------------

#[test]
fn solve_marks_project_dirty() {
    let mut app = new_app();
    assert!(!app.unsaved_changes());

    app.dispatch_action(RibbonAction::Solve);

    assert!(app.unsaved_changes());
    assert_eq!(app.project().status, SimulationStatus::Finished);
    assert!(app.project().last_result.is_some());
}

#[test]
fn stub_actions_mark_project_dirty() {
    let stub_actions = [
        RibbonAction::DrawBox,
        RibbonAction::DrawCylinder,
        RibbonAction::BoolUnite,
        RibbonAction::AssignMaterial,
    ];

    for action in stub_actions {
        let mut app = new_app();
        app.dispatch_action(action);
        assert!(
            app.unsaved_changes(),
            "{:?} should mark project dirty",
            action
        );
        assert!(
            app.status_text().contains("not implemented"),
            "{:?} should show stub status",
            action
        );
    }
}

#[test]
fn toggle_actions_do_not_mark_dirty() {
    let mut app = new_app();

    app.ribbon_state_mut()
        .toggles
        .insert("grid".to_string(), false);
    app.dispatch_action(RibbonAction::ToggleGrid);

    assert!(!app.unsaved_changes());
    assert!(app.status_text().contains("Grid:"));
}

#[test]
fn save_clears_dirty_flag() {
    let mut app = new_app();
    app.dispatch_action(RibbonAction::Solve);
    assert!(app.unsaved_changes());

    let dir = tempfile::tempdir().unwrap();
    app.save_to(&dir.path().join("clean.emsp"));
    assert!(!app.unsaved_changes());
}

#[test]
fn dirty_flag_shown_in_display_name() {
    let mut app = new_app();
    assert_eq!(app.file_display_name(), "Untitled");

    app.dispatch_action(RibbonAction::DrawBox);
    assert_eq!(app.file_display_name(), "Untitled*");

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("named.emsp");
    app.save_to(&path);
    assert_eq!(app.file_display_name(), "named.emsp");

    app.dispatch_action(RibbonAction::DrawSphere);
    assert_eq!(app.file_display_name(), "named.emsp*");
}

// ---------------------------------------------------------------------------
// Solve action
// ---------------------------------------------------------------------------

#[test]
fn solve_all_works_same_as_solve() {
    let mut app = new_app();
    app.dispatch_action(RibbonAction::SolveAll);

    assert_eq!(app.project().status, SimulationStatus::Finished);
    assert!(app.project().last_result.is_some());
    assert!(app.status_text().contains("Solve completed"));
    assert!(app.unsaved_changes());
}

// ---------------------------------------------------------------------------
// View / render toggle actions
// ---------------------------------------------------------------------------

#[test]
fn render_mode_action_updates_status() {
    let mut app = new_app();
    app.dispatch_action(RibbonAction::RenderShaded);
    assert!(app.status_text().contains("RenderShaded"));

    app.dispatch_action(RibbonAction::RenderWireframe);
    assert!(app.status_text().contains("RenderWireframe"));
}

// ---------------------------------------------------------------------------
// Full lifecycle: create → modify → save → new → open → verify
// ---------------------------------------------------------------------------

#[test]
fn full_project_lifecycle() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("lifecycle.emsp");

    // 1. Create app, solve (generates result data)
    let mut app = new_app();
    app.dispatch_action(RibbonAction::Solve);
    assert_eq!(app.project().status, SimulationStatus::Finished);
    assert!(app.unsaved_changes());

    // 2. Save
    app.save_to(&path);
    assert!(!app.unsaved_changes());
    assert!(path.exists());

    // 3. New project (clears everything)
    app.dispatch_action(RibbonAction::NewProject);
    assert_eq!(app.project().title, "New Project");
    assert!(app.project().last_result.is_none());
    assert!(app.current_file().is_none());

    // 4. Open the saved file
    app.open_from(&path);
    assert_eq!(app.project().status, SimulationStatus::Finished);
    assert!(app.project().last_result.is_some());
    assert!(!app.unsaved_changes());
    assert_eq!(app.current_file().unwrap(), &path);

    // 5. Make another change, verify dirty
    app.dispatch_action(RibbonAction::DrawBox);
    assert!(app.unsaved_changes());
    assert!(app.file_display_name().ends_with('*'));
}
