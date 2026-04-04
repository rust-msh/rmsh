use egui::{Color32, RichText, ScrollArea, Ui};

use crate::theme;

// ---------------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    Info,
    Warning,
    Error,
}

#[derive(Debug, Clone)]
pub struct MessageEntry {
    pub severity: Severity,
    pub text: String,
}

impl MessageEntry {
    pub fn info(text: impl Into<String>) -> Self {
        Self {
            severity: Severity::Info,
            text: text.into(),
        }
    }

    pub fn warning(text: impl Into<String>) -> Self {
        Self {
            severity: Severity::Warning,
            text: text.into(),
        }
    }

    pub fn error(text: impl Into<String>) -> Self {
        Self {
            severity: Severity::Error,
            text: text.into(),
        }
    }
}

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BottomTab {
    Messages,
    Log,
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

/// Renders the bottom message manager / log area.
pub fn show_message_manager(
    ui: &mut Ui,
    active_tab: &mut BottomTab,
    messages: &[MessageEntry],
    log_text: &str,
) {
    // Tab bar
    ui.horizontal(|ui| {
        if ui
            .selectable_label(*active_tab == BottomTab::Messages, "Messages")
            .clicked()
        {
            *active_tab = BottomTab::Messages;
        }
        if ui
            .selectable_label(*active_tab == BottomTab::Log, "Log")
            .clicked()
        {
            *active_tab = BottomTab::Log;
        }
    });
    ui.separator();

    match active_tab {
        BottomTab::Messages => {
            ScrollArea::vertical().auto_shrink(false).show(ui, |ui| {
                if messages.is_empty() {
                    ui.label(
                        RichText::new("No messages.")
                            .italics()
                            .color(Color32::from_rgb(140, 140, 140)),
                    );
                } else {
                    for msg in messages {
                        let (icon, color) = match msg.severity {
                            Severity::Error => ("\u{274C}", theme::SEVERITY_ERROR),
                            Severity::Warning => ("\u{26A0}", theme::SEVERITY_WARNING),
                            Severity::Info => ("\u{2139}", theme::SEVERITY_INFO),
                        };
                        ui.horizontal(|ui| {
                            ui.label(RichText::new(icon).color(color));
                            ui.label(&msg.text);
                        });
                    }
                }
            });
        }
        BottomTab::Log => {
            ScrollArea::vertical().auto_shrink(false).show(ui, |ui| {
                ui.monospace(log_text);
            });
        }
    }
}
