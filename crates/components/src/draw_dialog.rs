// ---------------------------------------------------------------------------
// Draw dialog — parameter input for geometry creation commands
// ---------------------------------------------------------------------------

use std::collections::HashMap;

use egui::Ui;

use emstudio_domain::geometry::{GeometryOperation, OperationCommand};

/// State for an active draw dialog.
pub struct DrawDialog {
    pub command: OperationCommand,
    pub title: String,
    pub fields: Vec<DialogField>,
    pub name_field: String,
    pub open: bool,
}

pub struct DialogField {
    pub key: String,
    pub label: String,
    pub value: String,
    pub placeholder: String,
}

/// Result when dialog is confirmed.
pub struct DrawDialogResult {
    pub operation: GeometryOperation,
}

impl DrawDialog {
    pub fn new_box() -> Self {
        Self {
            command: OperationCommand::CreateBox,
            title: "Create Box".into(),
            fields: vec![
                DialogField { key: "origin_x".into(), label: "Origin X".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "origin_y".into(), label: "Origin Y".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "origin_z".into(), label: "Origin Z".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "width".into(), label: "Width (X)".into(), value: "10".into(), placeholder: "10mm".into() },
                DialogField { key: "height".into(), label: "Height (Y)".into(), value: "10".into(), placeholder: "10mm".into() },
                DialogField { key: "depth".into(), label: "Depth (Z)".into(), value: "10".into(), placeholder: "10mm".into() },
            ],
            name_field: "Box1".into(),
            open: true,
        }
    }

    pub fn new_cylinder() -> Self {
        Self {
            command: OperationCommand::CreateCylinder,
            title: "Create Cylinder".into(),
            fields: vec![
                DialogField { key: "center_x".into(), label: "Center X".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "center_y".into(), label: "Center Y".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "center_z".into(), label: "Center Z".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "radius".into(), label: "Radius".into(), value: "5".into(), placeholder: "5mm".into() },
                DialogField { key: "height".into(), label: "Height".into(), value: "10".into(), placeholder: "10mm".into() },
            ],
            name_field: "Cylinder1".into(),
            open: true,
        }
    }

    pub fn new_sphere() -> Self {
        Self {
            command: OperationCommand::CreateSphere,
            title: "Create Sphere".into(),
            fields: vec![
                DialogField { key: "center_x".into(), label: "Center X".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "center_y".into(), label: "Center Y".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "center_z".into(), label: "Center Z".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "radius".into(), label: "Radius".into(), value: "5".into(), placeholder: "5mm".into() },
            ],
            name_field: "Sphere1".into(),
            open: true,
        }
    }

    pub fn new_cone() -> Self {
        Self {
            command: OperationCommand::CreateCone,
            title: "Create Cone".into(),
            fields: vec![
                DialogField { key: "center_x".into(), label: "Center X".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "center_y".into(), label: "Center Y".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "center_z".into(), label: "Center Z".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "radius".into(), label: "Base Radius".into(), value: "5".into(), placeholder: "5mm".into() },
                DialogField { key: "height".into(), label: "Height".into(), value: "10".into(), placeholder: "10mm".into() },
            ],
            name_field: "Cone1".into(),
            open: true,
        }
    }

    pub fn new_torus() -> Self {
        Self {
            command: OperationCommand::CreateTorus,
            title: "Create Torus".into(),
            fields: vec![
                DialogField { key: "center_x".into(), label: "Center X".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "center_y".into(), label: "Center Y".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "center_z".into(), label: "Center Z".into(), value: "0".into(), placeholder: "0".into() },
                DialogField { key: "major_radius".into(), label: "Major Radius".into(), value: "10".into(), placeholder: "10mm".into() },
                DialogField { key: "minor_radius".into(), label: "Minor Radius".into(), value: "3".into(), placeholder: "3mm".into() },
            ],
            name_field: "Torus1".into(),
            open: true,
        }
    }

    /// Show the dialog as an egui::Window. Returns Some(result) if confirmed.
    pub fn show(
        &mut self,
        ctx: &egui::Context,
        variable_names: &[String],
    ) -> Option<DrawDialogResult> {
        if !self.open {
            return None;
        }

        let mut result = None;
        let title = self.title.clone();
        let mut close = false;

        egui::Window::new(title)
            .resizable(false)
            .collapsible(false)
            .show(ctx, |ui| {
                // Name field
                ui.horizontal(|ui| {
                    ui.label("Name:");
                    ui.text_edit_singleline(&mut self.name_field);
                });
                ui.separator();

                // Parameter fields
                egui::Grid::new("draw_dialog_grid")
                    .num_columns(2)
                    .spacing([8.0, 4.0])
                    .show(ui, |ui| {
                        for field in &mut self.fields {
                            ui.label(&field.label);
                            ui.horizontal(|ui| {
                                let text_edit = egui::TextEdit::singleline(&mut field.value)
                                    .desired_width(100.0)
                                    .hint_text(&field.placeholder);
                                ui.add(text_edit);

                                // Variable quick-insert dropdown
                                if !variable_names.is_empty() {
                                    egui::ComboBox::from_id_salt(format!("var_{}", field.key))
                                        .selected_text("")
                                        .width(20.0)
                                        .show_ui(ui, |ui| {
                                            for var_name in variable_names {
                                                if ui.selectable_label(false, var_name).clicked() {
                                                    field.value = var_name.clone();
                                                }
                                            }
                                        });
                                }
                            });
                            ui.end_row();
                        }
                    });

                ui.separator();
                ui.horizontal(|ui| {
                    if ui.button("OK").clicked() {
                        result = Some(self.build_operation());
                        close = true;
                    }
                    if ui.button("Cancel").clicked() {
                        close = true;
                    }
                });
            });

        if close {
            self.open = false;
        }

        if !self.open && result.is_none() {
            // Dialog was closed without confirming
            return None;
        }

        result
    }

    fn build_operation(&self) -> DrawDialogResult {
        let mut params = serde_json::Map::new();

        // Build origin/center from _x/_y/_z fields
        match self.command {
            OperationCommand::CreateBox => {
                let ox = self.get_field("origin_x");
                let oy = self.get_field("origin_y");
                let oz = self.get_field("origin_z");
                params.insert("origin".into(), build_vec3_value(&ox, &oy, &oz));
                self.insert_param(&mut params, "width");
                self.insert_param(&mut params, "height");
                self.insert_param(&mut params, "depth");
            }
            OperationCommand::CreateCylinder | OperationCommand::CreateSphere | OperationCommand::CreateCone => {
                let cx = self.get_field("center_x");
                let cy = self.get_field("center_y");
                let cz = self.get_field("center_z");
                params.insert("center".into(), build_vec3_value(&cx, &cy, &cz));
                self.insert_param(&mut params, "radius");
                if self.command != OperationCommand::CreateSphere {
                    self.insert_param(&mut params, "height");
                }
            }
            OperationCommand::CreateTorus => {
                let cx = self.get_field("center_x");
                let cy = self.get_field("center_y");
                let cz = self.get_field("center_z");
                params.insert("center".into(), build_vec3_value(&cx, &cy, &cz));
                self.insert_param(&mut params, "major_radius");
                self.insert_param(&mut params, "minor_radius");
            }
            _ => {
                // Generic: just pass all fields
                for field in &self.fields {
                    insert_smart_value(&mut params, &field.key, &field.value);
                }
            }
        }

        let step = 0; // Will be assigned by the app

        DrawDialogResult {
            operation: GeometryOperation {
                step,
                command: self.command,
                result_object: Some(self.name_field.clone()),
                parameters: serde_json::Value::Object(params),
                attributes: None,
            },
        }
    }

    fn get_field(&self, key: &str) -> String {
        self.fields
            .iter()
            .find(|f| f.key == key)
            .map(|f| f.value.clone())
            .unwrap_or_default()
    }

    fn insert_param(&self, params: &mut serde_json::Map<String, serde_json::Value>, key: &str) {
        if let Some(field) = self.fields.iter().find(|f| f.key == key) {
            insert_smart_value(params, key, &field.value);
        }
    }
}

/// Insert a value as number if parseable, otherwise as string (for expression support).
fn insert_smart_value(params: &mut serde_json::Map<String, serde_json::Value>, key: &str, value: &str) {
    if let Ok(n) = value.parse::<f64>() {
        params.insert(key.into(), serde_json::Value::from(n));
    } else {
        params.insert(key.into(), serde_json::Value::from(value));
    }
}

/// Build a vec3 JSON value from three component strings.
fn build_vec3_value(x: &str, y: &str, z: &str) -> serde_json::Value {
    let to_json = |s: &str| -> serde_json::Value {
        if let Ok(n) = s.parse::<f64>() {
            serde_json::Value::from(n)
        } else {
            serde_json::Value::from(s)
        }
    };
    serde_json::Value::Array(vec![to_json(x), to_json(y), to_json(z)])
}
