use egui::Ui;

use emstudio_domain::Edition;

use crate::ribbon::RibbonAction;

/// Traditional menu bar replicating AEDT's menu structure.
/// Returns an action if a menu item was clicked.
///
/// Items are gated by `edition` (feature tier) and `is_web` (web vs native).
pub fn show_menu_bar(ui: &mut Ui, edition: Edition, is_web: bool) -> Option<RibbonAction> {
    let mut action = None;
    let file_ops = !is_web;

    egui::MenuBar::new().ui(ui, |ui| {
        ui.menu_button("File", |ui| {
            if ui.add_enabled(file_ops, egui::Button::new("New Project")).clicked() {
                action = Some(RibbonAction::NewProject);
                ui.close();
            }
            if ui.add_enabled(file_ops, egui::Button::new("Open...")).clicked() {
                action = Some(RibbonAction::OpenProject);
                ui.close();
            }
            if ui.add_enabled(file_ops, egui::Button::new("Save")).clicked() {
                action = Some(RibbonAction::SaveProject);
                ui.close();
            }
            if ui.add_enabled(file_ops, egui::Button::new("Save As...")).clicked() {
                action = Some(RibbonAction::SaveAs);
                ui.close();
            }
            ui.separator();
            if ui.button("Close Project").clicked() {
                action = Some(RibbonAction::CloseProject);
                ui.close();
            }
            ui.separator();
            ui.menu_button("Import", |ui| {
                if ui.add_enabled(file_ops, egui::Button::new("Import STEP")).clicked() {
                    action = Some(RibbonAction::ImportStep);
                    ui.close();
                }
                if ui.add_enabled(file_ops, egui::Button::new("Import SAT")).clicked() {
                    action = Some(RibbonAction::ImportSat);
                    ui.close();
                }
            });
            ui.menu_button("Export", |ui| {
                if ui.add_enabled(file_ops, egui::Button::new("Export STEP")).clicked() {
                    action = Some(RibbonAction::ExportStep);
                    ui.close();
                }
                if ui.add_enabled(file_ops, egui::Button::new("Export SAT")).clicked() {
                    action = Some(RibbonAction::ExportSat);
                    ui.close();
                }
            });
        });

        ui.menu_button("Edit", |ui| {
            ui.add_enabled(false, egui::Button::new("Undo"));
            ui.add_enabled(false, egui::Button::new("Redo"));
        });

        ui.menu_button("View", |ui| {
            if ui.button("Toggle Grid").clicked() {
                action = Some(RibbonAction::ToggleGrid);
                ui.close();
            }
            if ui.button("Toggle Ruler").clicked() {
                action = Some(RibbonAction::ToggleRuler);
                ui.close();
            }
            if ui.button("Toggle Axes").clicked() {
                action = Some(RibbonAction::ToggleCoordSystem);
                ui.close();
            }
            ui.separator();
            if ui.button("Shaded").clicked() {
                action = Some(RibbonAction::RenderShaded);
                ui.close();
            }
            if ui.button("Wireframe").clicked() {
                action = Some(RibbonAction::RenderWireframe);
                ui.close();
            }
            ui.separator();
            if ui.button("Fit All").clicked() {
                action = Some(RibbonAction::FitAll);
                ui.close();
            }
            if ui.button("Zoom In").clicked() {
                action = Some(RibbonAction::ZoomIn);
                ui.close();
            }
            if ui.button("Zoom Out").clicked() {
                action = Some(RibbonAction::ZoomOut);
                ui.close();
            }
            ui.separator();
            if ui.button("Toggle Project Manager").clicked() {
                action = Some(RibbonAction::ToggleProjectManager);
                ui.close();
            }
            if ui.button("Toggle Message Manager").clicked() {
                action = Some(RibbonAction::ToggleMessageManager);
                ui.close();
            }
        });

        ui.menu_button("Project", |ui| {
            if ui.button("Validate").clicked() {
                action = Some(RibbonAction::Validate);
                ui.close();
            }
            if ui.button("Add Setup").clicked() {
                action = Some(RibbonAction::AddSetup);
                ui.close();
            }
            if ui.button("Add Sweep").clicked() {
                action = Some(RibbonAction::AddSweep);
                ui.close();
            }
        });

        ui.menu_button("Draw", |ui| {
            if ui.button("Box").clicked() {
                action = Some(RibbonAction::DrawBox);
                ui.close();
            }
            if ui.button("Cylinder").clicked() {
                action = Some(RibbonAction::DrawCylinder);
                ui.close();
            }
            if ui.button("Sphere").clicked() {
                action = Some(RibbonAction::DrawSphere);
                ui.close();
            }
            if ui.button("Cone").clicked() {
                action = Some(RibbonAction::DrawCone);
                ui.close();
            }
            if ui.button("Torus").clicked() {
                action = Some(RibbonAction::DrawTorus);
                ui.close();
            }
            ui.separator();
            if ui.button("Rectangle").clicked() {
                action = Some(RibbonAction::DrawRectangle);
                ui.close();
            }
            if ui.button("Ellipse").clicked() {
                action = Some(RibbonAction::DrawEllipse);
                ui.close();
            }
            if ui.button("Circle").clicked() {
                action = Some(RibbonAction::DrawCircle);
                ui.close();
            }
            if ui.button("Polygon").clicked() {
                action = Some(RibbonAction::DrawPolygon);
                ui.close();
            }
            if ui.button("Polyline").clicked() {
                action = Some(RibbonAction::DrawPolyline);
                ui.close();
            }
            if ui.button("Arc").clicked() {
                action = Some(RibbonAction::DrawArc);
                ui.close();
            }
            if ui.button("Spline").clicked() {
                action = Some(RibbonAction::DrawSpline);
                ui.close();
            }
        });

        ui.menu_button("Simulation", |ui| {
            if ui.button("Solve").clicked() {
                action = Some(RibbonAction::Solve);
                ui.close();
            }
            if ui.add_enabled(edition.allows_solve_all(), egui::Button::new("Solve All")).clicked() {
                action = Some(RibbonAction::SolveAll);
                ui.close();
            }
            ui.separator();
            ui.add_enabled_ui(false, |ui| {
            let _ = ui.button("Abort");
            });
        });

        ui.menu_button("Results", |ui| {
            if ui.button("Create Report").clicked() {
                action = Some(RibbonAction::CreateReport);
                ui.close();
            }
            if ui.button("Solution Data").clicked() {
                action = Some(RibbonAction::SolutionData);
                ui.close();
            }
            ui.separator();
            ui.menu_button("Plot Fields", |ui| {
                if ui.button("E-Field").clicked() {
                    action = Some(RibbonAction::PlotEField);
                    ui.close();
                }
                if ui.button("H-Field").clicked() {
                    action = Some(RibbonAction::PlotHField);
                    ui.close();
                }
                if ui.button("SAR").clicked() {
                    action = Some(RibbonAction::PlotSAR);
                    ui.close();
                }
            });
            if ui.button("Animate").clicked() {
                action = Some(RibbonAction::Animate);
                ui.close();
            }
        });

        ui.menu_button("Help", |ui| {
            ui.add_enabled(false, egui::Button::new("About EmStudio"));
        });
    });

    action
}
