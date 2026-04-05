use std::collections::HashMap;

use egui::Ui;

use emstudio_domain::geometry::GeometryOperation;
use emstudio_domain::variable::Variable;
use emstudio_domain::Project;

use crate::project_tree::{self, ProjectTreeResponse};
use crate::properties_panel::{self, PropertiesResponse};

// ---------------------------------------------------------------------------
// Left panel tab enum
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LeftPanelTab {
    ProjectManager,
    Properties,
}

// ---------------------------------------------------------------------------
// Left dock panel context
// ---------------------------------------------------------------------------

/// All state needed by the left dock panel.
pub struct LeftDockContext<'a> {
    pub project: &'a Project,
    pub selected_object: Option<&'a str>,
    pub selected_operation: Option<&'a GeometryOperation>,
    pub design_variables: &'a HashMap<String, Variable>,
    pub variable_edit_buffers: &'a mut HashMap<String, String>,
    pub param_edit_buffers: &'a mut HashMap<String, String>,
}

/// Combined response from left dock panel.
pub struct LeftDockResponse {
    pub tree: ProjectTreeResponse,
    pub properties: PropertiesResponse,
}

// ---------------------------------------------------------------------------
// Left dock panel (Project Manager + Properties as tabs)
// ---------------------------------------------------------------------------

/// Renders the left dock panel with a tab group:
/// - Project Manager (tree view)
/// - Properties (attributes/definition)
pub fn left_dock_panel(
    ui: &mut Ui,
    active_tab: &mut LeftPanelTab,
    ctx: &mut LeftDockContext<'_>,
) -> LeftDockResponse {
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

    let mut tree_response = ProjectTreeResponse {
        selected_object: None,
        add_variable: false,
        variable_edited: None,
    };
    let mut props_response = PropertiesResponse {
        parameter_changed: None,
    };

    // Tab content
    egui::ScrollArea::vertical()
        .auto_shrink(false)
        .show(ui, |ui| match active_tab {
            LeftPanelTab::ProjectManager => {
                tree_response = project_tree::show_project_tree(
                    ui,
                    ctx.project,
                    ctx.selected_object,
                    ctx.design_variables,
                    ctx.variable_edit_buffers,
                );
            }
            LeftPanelTab::Properties => {
                props_response = properties_panel::show_properties_panel(
                    ui,
                    ctx.project,
                    ctx.selected_object,
                    ctx.selected_operation,
                    ctx.param_edit_buffers,
                );
            }
        });

    LeftDockResponse {
        tree: tree_response,
        properties: props_response,
    }
}
