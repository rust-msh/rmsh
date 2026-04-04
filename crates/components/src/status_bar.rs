use egui::{Color32, RichText, Ui};

// ---------------------------------------------------------------------------
// State / Response types
// ---------------------------------------------------------------------------

pub struct StatusBarState {
    pub file_name: String,
    pub unsaved: bool,
    pub status_text: String,
    pub coordinates: Option<(f32, f32, f32)>,
    pub units: String,
    pub message_manager_visible: bool,
}

pub struct StatusBarResponse {
    pub toggle_message_manager: bool,
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

/// Renders the AEDT-style status bar at the bottom of the window.
pub fn show_status_bar(ui: &mut Ui, state: &StatusBarState) -> StatusBarResponse {
    let mut response = StatusBarResponse {
        toggle_message_manager: false,
    };

    ui.horizontal(|ui| {
        // File name with unsaved indicator
        let file_text = if state.unsaved {
            format!("{}*", state.file_name)
        } else {
            state.file_name.clone()
        };
        let file_color = if state.unsaved {
            Color32::from_rgb(200, 120, 0)
        } else {
            Color32::from_rgb(100, 100, 100)
        };
        ui.label(RichText::new(format!("[{}]", file_text)).strong().color(file_color));

        ui.separator();

        // Status text
        ui.label(&state.status_text);

        // Right-aligned section
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            // Message Manager toggle button
            let msg_label = if state.message_manager_visible {
                "Messages \u{25BC}"
            } else {
                "Messages \u{25B6}"
            };
            if ui
                .small_button(msg_label)
                .on_hover_text("Toggle Message Manager")
                .clicked()
            {
                response.toggle_message_manager = true;
            }

            ui.separator();

            // Units
            ui.label(
                RichText::new(&state.units)
                    .color(Color32::from_rgb(80, 80, 80))
                    .strong(),
            );

            ui.separator();

            // Coordinates
            if let Some((x, y, z)) = state.coordinates {
                ui.label(
                    RichText::new(format!("X:{x:.2}  Y:{y:.2}  Z:{z:.2}"))
                        .color(Color32::from_rgb(80, 80, 80))
                        .monospace(),
                );
            } else {
                ui.label(
                    RichText::new("X:--  Y:--  Z:--")
                        .color(Color32::from_rgb(160, 160, 160))
                        .monospace(),
                );
            }
        });
    });

    response
}
