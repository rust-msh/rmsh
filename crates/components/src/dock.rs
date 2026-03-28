use egui::Ui;

use emstudio_domain::Project;

pub fn left_panel(ui: &mut Ui, project: &Project) {
    ui.heading("Model Tree");
    ui.separator();
    ui.label(format!("Project: {}", project.title));
    ui.label(format!("Objects: {}", project.model.objects.len()));
    for object in &project.model.objects {
        ui.label(format!("- {}", object.name));
    }
}

pub fn right_panel(ui: &mut Ui, project: &Project) {
    ui.heading("Properties");
    ui.separator();
    ui.label(format!("Model: {}", project.model.name));
    ui.label(format!("Materials: {}", project.model.materials.len()));
}

pub fn bottom_panel(ui: &mut Ui, status: &str) {
    ui.horizontal(|ui| {
        ui.label("Status:");
        ui.monospace(status);
    });
}
