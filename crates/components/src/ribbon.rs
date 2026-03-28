use egui::{Button, RichText, Ui};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RibbonAction {
    NewProject,
    OpenProject,
    SaveProject,
    Solve,
}

pub fn show_ribbon(ui: &mut Ui) -> Option<RibbonAction> {
    let mut action = None;
    ui.horizontal_wrapped(|ui| {
        ui.spacing_mut().item_spacing.x = 10.0;

        if ui
            .add(Button::new(RichText::new("New").strong()).min_size(egui::vec2(56.0, 36.0)))
            .clicked()
        {
            action = Some(RibbonAction::NewProject);
        }

        if ui
            .add(Button::new("Open").min_size(egui::vec2(56.0, 36.0)))
            .clicked()
        {
            action = Some(RibbonAction::OpenProject);
        }

        if ui
            .add(Button::new("Save").min_size(egui::vec2(56.0, 36.0)))
            .clicked()
        {
            action = Some(RibbonAction::SaveProject);
        }

        ui.separator();

        if ui
            .add(
                Button::new(RichText::new("Solve").strong())
                    .min_size(egui::vec2(72.0, 36.0))
                    .fill(egui::Color32::from_rgb(16, 103, 72)),
            )
            .clicked()
        {
            action = Some(RibbonAction::Solve);
        }
    });

    action
}
