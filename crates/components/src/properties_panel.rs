// ---------------------------------------------------------------------------
// Properties panel — shows attributes of the selected geometry object
// ---------------------------------------------------------------------------

use egui::Ui;

use emstudio_domain::geometry::GeometryOperation;
use emstudio_domain::{GeometryObject, Project};

/// Response from the properties panel when user edits a parameter.
pub struct PropertiesResponse {
    /// If set, the user modified a geometry parameter and a rebuild is needed.
    pub parameter_changed: Option<ParameterEdit>,
}

pub struct ParameterEdit {
    pub object_name: String,
    pub key: String,
    pub new_value: String,
}

/// Renders the Properties panel showing details of the current model/selection.
pub fn show_properties_panel(
    ui: &mut Ui,
    project: &Project,
    selected_object: Option<&str>,
    selected_operation: Option<&GeometryOperation>,
    param_edit_buffers: &mut std::collections::HashMap<String, String>,
) -> PropertiesResponse {
    let mut response = PropertiesResponse {
        parameter_changed: None,
    };

    if let Some(obj_name) = selected_object {
        // Find the object in the model
        let obj = project.model.objects.iter().find(|o| o.name == obj_name);

        if let Some(obj) = obj {
            show_object_properties(ui, obj);
            ui.add_space(8.0);

            // Show editable parameters from the operation
            if let Some(op) = selected_operation {
                if let Some(edit) = show_operation_parameters(ui, obj, op, param_edit_buffers) {
                    response.parameter_changed = Some(edit);
                }
            }
        } else {
            ui.label(format!("Object '{}' not found", obj_name));
        }
    } else {
        show_model_summary(ui, project);
    }

    response
}

fn show_model_summary(ui: &mut Ui, project: &Project) {
    ui.heading("Attributes");
    ui.separator();

    egui::Grid::new("props_grid")
        .num_columns(2)
        .spacing([8.0, 4.0])
        .show(ui, |ui| {
            ui.label("Model:");
            ui.label(&project.model.name);
            ui.end_row();

            ui.label("Objects:");
            ui.label(format!("{}", project.model.objects.len()));
            ui.end_row();

            ui.label("Materials:");
            ui.label(format!("{}", project.model.materials.len()));
            ui.end_row();

            ui.label("Status:");
            ui.label(format!("{:?}", project.status));
            ui.end_row();
        });

    ui.add_space(12.0);
    ui.heading("Definition");
    ui.separator();
    ui.label(
        egui::RichText::new("Select an object to view its parameters.")
            .italics()
            .color(egui::Color32::from_rgb(140, 140, 140)),
    );
}

fn show_object_properties(ui: &mut Ui, obj: &GeometryObject) {
    ui.heading("Object Properties");
    ui.separator();

    egui::Grid::new("obj_props_grid")
        .num_columns(2)
        .spacing([8.0, 4.0])
        .show(ui, |ui| {
            ui.label("Name:");
            ui.strong(&obj.name);
            ui.end_row();

            ui.label("ID:");
            ui.label(format!("{}", obj.id));
            ui.end_row();

            ui.label("Mesh hint:");
            ui.label(&obj.mesh_hint);
            ui.end_row();
        });
}

fn show_operation_parameters(
    ui: &mut Ui,
    obj: &GeometryObject,
    op: &GeometryOperation,
    param_buffers: &mut std::collections::HashMap<String, String>,
) -> Option<ParameterEdit> {
    ui.heading("Parameters");
    ui.separator();

    let mut edit_result = None;

    ui.label(
        egui::RichText::new(format!("Command: {:?}", op.command))
            .color(egui::Color32::from_rgb(100, 180, 255)),
    );
    ui.add_space(4.0);

    if let Some(map) = op.parameters.as_object() {
        egui::Grid::new("param_edit_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                for (key, value) in map {
                    ui.label(format!("{}:", key));

                    let buf_key = format!("{}:{}", obj.name, key);
                    let display_value = format_json_value(value);

                    let buf = param_buffers
                        .entry(buf_key.clone())
                        .or_insert_with(|| display_value.clone());

                    let text_edit = egui::TextEdit::singleline(buf)
                        .desired_width(120.0);
                    let response = ui.add(text_edit);

                    if response.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                        if *buf != display_value {
                            edit_result = Some(ParameterEdit {
                                object_name: obj.name.clone(),
                                key: key.clone(),
                                new_value: buf.clone(),
                            });
                        }
                    }
                    ui.end_row();
                }
            });
    }

    edit_result
}

fn format_json_value(value: &serde_json::Value) -> String {
    match value {
        serde_json::Value::Number(n) => n.to_string(),
        serde_json::Value::String(s) => s.clone(),
        serde_json::Value::Array(arr) => {
            let parts: Vec<String> = arr.iter().map(|v| format_json_value(v)).collect();
            format!("[{}]", parts.join(", "))
        }
        other => other.to_string(),
    }
}
