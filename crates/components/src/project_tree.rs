// ---------------------------------------------------------------------------
// Project tree — hierarchical AEDT-style project manager
// ---------------------------------------------------------------------------

use std::collections::HashMap;

use egui::Ui;

use emstudio_domain::variable::Variable;
use emstudio_domain::Project;

/// Response from the project tree when user interacts.
pub struct ProjectTreeResponse {
    /// User clicked on a geometry object.
    pub selected_object: Option<String>,
    /// User wants to add a new variable.
    pub add_variable: bool,
    /// User edited a variable value.
    pub variable_edited: Option<(String, String)>,
}

/// Renders a hierarchical project tree using collapsing headers.
pub fn show_project_tree(
    ui: &mut Ui,
    project: &Project,
    selected_object: Option<&str>,
    design_variables: &HashMap<String, Variable>,
    variable_edit_buffers: &mut HashMap<String, String>,
) -> ProjectTreeResponse {
    let mut response = ProjectTreeResponse {
        selected_object: None,
        add_variable: false,
        variable_edited: None,
    };

    egui::CollapsingHeader::new(
        egui::RichText::new(format!("\u{1F4C1} {}", project.title)).strong(),
    )
    .default_open(true)
    .show(ui, |ui| {
        // Variables
        show_variables_node(ui, design_variables, variable_edit_buffers, &mut response);

        // Geometry / 3D Objects
        show_geometry_node(ui, project, selected_object, &mut response);

        // Materials
        egui::CollapsingHeader::new("\u{1F3A8} Materials")
            .default_open(false)
            .show(ui, |ui| {
                if project.model.materials.is_empty() {
                    ui.label(
                        egui::RichText::new("  (no materials)")
                            .italics()
                            .color(egui::Color32::from_rgb(140, 140, 140)),
                    );
                } else {
                    for mat in &project.model.materials {
                        ui.horizontal(|ui| {
                            ui.add_space(8.0);
                            ui.label(format!("\u{25C9} {}", mat.name));
                        });
                    }
                }
            });

        // Analysis
        egui::CollapsingHeader::new("\u{2699} Analysis")
            .default_open(false)
            .show(ui, |ui| {
                ui.label(
                    egui::RichText::new("  (no setups)")
                        .italics()
                        .color(egui::Color32::from_rgb(140, 140, 140)),
                );
            });

        // Results
        egui::CollapsingHeader::new("\u{1F4CA} Results")
            .default_open(false)
            .show(ui, |ui| {
                if let Some(result) = &project.last_result {
                    ui.label(format!(
                        "  Converged: {}",
                        if result.converged { "Yes" } else { "No" }
                    ));
                } else {
                    ui.label(
                        egui::RichText::new("  (no results)")
                            .italics()
                            .color(egui::Color32::from_rgb(140, 140, 140)),
                    );
                }
            });
    });

    response
}

fn show_variables_node(
    ui: &mut Ui,
    variables: &HashMap<String, Variable>,
    edit_buffers: &mut HashMap<String, String>,
    response: &mut ProjectTreeResponse,
) {
    egui::CollapsingHeader::new("\u{1D465} Variables")
        .default_open(true)
        .show(ui, |ui| {
            if variables.is_empty() {
                ui.horizontal(|ui| {
                    ui.add_space(8.0);
                    ui.label(
                        egui::RichText::new("(no variables)")
                            .italics()
                            .color(egui::Color32::from_rgb(140, 140, 140)),
                    );
                });
            } else {
                egui::Grid::new("var_grid")
                    .num_columns(2)
                    .spacing([8.0, 2.0])
                    .show(ui, |ui| {
                        for (name, var) in variables {
                            ui.horizontal(|ui| {
                                ui.add_space(8.0);
                                ui.label(
                                    egui::RichText::new(name)
                                        .color(egui::Color32::from_rgb(100, 200, 100)),
                                );
                            });

                            let display_value = var
                                .value
                                .as_deref()
                                .or(var.expression.as_deref())
                                .unwrap_or("0");
                            let buf = edit_buffers
                                .entry(name.clone())
                                .or_insert_with(|| display_value.to_string());

                            let text_edit = egui::TextEdit::singleline(buf)
                                .desired_width(100.0);
                            let r = ui.add(text_edit);

                            if r.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                                if buf.as_str() != display_value {
                                    response.variable_edited =
                                        Some((name.clone(), buf.clone()));
                                }
                            }
                            ui.end_row();
                        }
                    });
            }

            ui.horizontal(|ui| {
                ui.add_space(8.0);
                if ui.small_button("+ Add Variable").clicked() {
                    response.add_variable = true;
                }
            });
        });
}

fn show_geometry_node(
    ui: &mut Ui,
    project: &Project,
    selected_object: Option<&str>,
    response: &mut ProjectTreeResponse,
) {
    egui::CollapsingHeader::new("\u{1F4D0} Geometry")
        .default_open(true)
        .show(ui, |ui| {
            if project.model.objects.is_empty() {
                ui.label(
                    egui::RichText::new("  (no objects)")
                        .italics()
                        .color(egui::Color32::from_rgb(140, 140, 140)),
                );
            } else {
                for obj in &project.model.objects {
                    let is_selected = selected_object == Some(obj.name.as_str());
                    ui.horizontal(|ui| {
                        ui.add_space(8.0);
                        let label = if is_selected {
                            egui::RichText::new(format!("\u{25A3} {}", obj.name))
                                .strong()
                                .color(egui::Color32::from_rgb(100, 180, 255))
                        } else {
                            egui::RichText::new(format!("\u{25A3} {}", obj.name))
                        };
                        if ui.selectable_label(is_selected, label).clicked() {
                            response.selected_object = Some(obj.name.clone());
                        }
                    });
                }
            }
        });
}
