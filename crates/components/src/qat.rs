use egui::{
    Color32, CornerRadius, FontId, Rect, Sense, Stroke, StrokeKind, Ui, vec2,
};

use crate::ribbon::RibbonAction;
use crate::theme;

const QAT_HEIGHT: f32 = 24.0;
const QAT_BTN_SIZE: f32 = 22.0;
const QAT_ICON_FONT: f32 = 14.0;

struct QatButton {
    icon: char,
    tooltip: &'static str,
    action: Option<RibbonAction>,
    enabled: bool,
}

/// Quick Access Toolbar - thin strip of icon-only buttons for frequent actions.
pub fn show_qat(ui: &mut Ui) -> Option<RibbonAction> {
    let mut result = None;

    // Background
    let rect = ui.available_rect_before_wrap();
    ui.painter().rect_filled(
        Rect::from_min_size(rect.min, vec2(rect.width(), QAT_HEIGHT)),
        CornerRadius::ZERO,
        theme::QAT_BG,
    );

    let buttons = [
        QatButton {
            icon: '\u{1F4C4}',
            tooltip: "New Project",
            action: Some(RibbonAction::NewProject),
            enabled: true,
        },
        QatButton {
            icon: '\u{1F4C2}',
            tooltip: "Open Project",
            action: Some(RibbonAction::OpenProject),
            enabled: true,
        },
        QatButton {
            icon: '\u{1F4BE}',
            tooltip: "Save Project",
            action: Some(RibbonAction::SaveProject),
            enabled: true,
        },
    ];

    ui.horizontal(|ui| {
        ui.set_height(QAT_HEIGHT);
        ui.add_space(4.0);

        for btn in &buttons {
            let (btn_rect, response) = ui.allocate_exact_size(
                vec2(QAT_BTN_SIZE, QAT_BTN_SIZE),
                if btn.enabled {
                    Sense::click()
                } else {
                    Sense::hover()
                },
            );

            // Paint hover/pressed bg
            if btn.enabled {
                if response.is_pointer_button_down_on() {
                    ui.painter().rect(
                        btn_rect,
                        CornerRadius::same(2),
                        theme::PRESSED_FILL,
                        Stroke::new(1.0, theme::HOVER_STROKE),
                        StrokeKind::Outside,
                    );
                } else if response.hovered() {
                    ui.painter().rect(
                        btn_rect,
                        CornerRadius::same(2),
                        theme::HOVER_FILL,
                        Stroke::new(1.0, theme::HOVER_STROKE),
                        StrokeKind::Outside,
                    );
                }
            }

            let text_color = if btn.enabled {
                Color32::from_rgb(30, 30, 30)
            } else {
                theme::DISABLED_TEXT
            };

            let galley = ui.painter().layout_no_wrap(
                btn.icon.to_string(),
                FontId::proportional(QAT_ICON_FONT),
                text_color,
            );
            let pos = btn_rect.center() - galley.size() / 2.0;
            ui.painter().galley(pos, galley, text_color);

            if btn.enabled && response.clicked() {
                if let Some(action) = btn.action {
                    result = Some(action);
                }
            }

            response.on_hover_text(btn.tooltip);
        }
    });

    result
}
