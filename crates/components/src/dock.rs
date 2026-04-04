use egui::Ui;

use emstudio_domain::Project;

use crate::project_tree;
use crate::properties_panel;

// ---------------------------------------------------------------------------
// Left panel tab enum
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LeftPanelTab {
    ProjectManager,
    Properties,
}

// ---------------------------------------------------------------------------
// Left dock panel (Project Manager + Properties as tabs)
// ---------------------------------------------------------------------------

/// Renders the left dock panel with a tab group:
/// - Project Manager (tree view)
/// - Properties (attributes/definition)
pub fn left_dock_panel(ui: &mut Ui, project: &Project, active_tab: &mut LeftPanelTab) {
    // Tab strip at top
    ui.horizontal(|ui| {
        if ui
            .selectable_label(*active_tab == LeftPanelTab::ProjectManager, "Project Manager")
            .clicked()
        {
            *active_tab = LeftPanelTab::ProjectManager;
        }
        if ui
            .selectable_label(*active_tab == LeftPanelTab::Properties, "Properties")
            .clicked()
        {
            *active_tab = LeftPanelTab::Properties;
        }
    });
    ui.separator();

    // Tab content
    egui::ScrollArea::vertical()
        .auto_shrink(false)
        .show(ui, |ui| match active_tab {
            LeftPanelTab::ProjectManager => {
                project_tree::show_project_tree(ui, project);
            }
            LeftPanelTab::Properties => {
                properties_panel::show_properties_panel(ui, project);
            }
        });
}
