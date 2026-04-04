use egui::Ui;

use emstudio_domain::Project;

/// Renders a hierarchical project tree using collapsing headers.
/// Mirrors the AEDT Project Manager tree structure.
pub fn show_project_tree(ui: &mut Ui, project: &Project) {
    egui::CollapsingHeader::new(
        egui::RichText::new(format!("\u{1F4C1} {}", project.title)).strong(),
    )
    .default_open(true)
    .show(ui, |ui| {
        // Geometry / 3D Objects
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
                        ui.horizontal(|ui| {
                            ui.add_space(8.0);
                            ui.label(format!("\u{25A3} {}", obj.name));
                        });
                    }
                }
            });

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

        // Simulation Setup
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
}
