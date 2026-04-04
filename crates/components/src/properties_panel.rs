use egui::Ui;

use emstudio_domain::Project;

/// Renders the Properties panel showing details of the current model/selection.
/// Mirrors the AEDT Properties Window.
pub fn show_properties_panel(ui: &mut Ui, project: &Project) {
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
