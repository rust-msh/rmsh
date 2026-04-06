use std::collections::HashMap;

use egui::{
    Align, Color32, CornerRadius, FontId, Id, Layout, Rect, Response, RichText, Sense, Stroke,
    StrokeKind, Ui, vec2,
};

use crate::theme;

// ---------------------------------------------------------------------------
// Layout constants (AEDT-inspired, scaled for egui)
// ---------------------------------------------------------------------------

const TAB_STRIP_HEIGHT: f32 = 26.0;
const COMMAND_AREA_HEIGHT: f32 = 74.0;
const GROUP_LABEL_HEIGHT: f32 = 16.0;
const LARGE_BTN_WIDTH: f32 = 52.0;
const LARGE_BTN_HEIGHT: f32 = 58.0; // icon + 2-line label
const SMALL_ROW_HEIGHT: f32 = 20.0;
const GROUP_H_PAD: f32 = 6.0;
const ICON_FONT_LARGE: f32 = 22.0;
const ICON_FONT_SMALL: f32 = 14.0;
const LABEL_FONT: f32 = 11.0;
const GROUP_LABEL_FONT: f32 = 10.0;

// Re-import shared theme colors for local use
use theme::{
    CHECKED_FILL, CHECKED_STROKE, DISABLED_TEXT, GROUP_LABEL_COLOR, GROUP_SEP_COLOR, HOVER_FILL,
    HOVER_STROKE, MENU_HOVER, PRESSED_FILL, RIBBON_BG, TAB_ACTIVE_TEXT, TAB_NORMAL_TEXT,
};

// ---------------------------------------------------------------------------
// Actions
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RibbonAction {
    // Desktop
    NewProject,
    OpenProject,
    SaveProject,
    SaveAs,
    CloseProject,
    ImportStep,
    ImportSat,
    ImportHfssAedt,
    ImportQ3dAedt,
    ExportStep,
    ExportSat,
    ExportHfssPyAedt,
    ExportQ3dPyAedt,
    // View
    ToggleGrid,
    ToggleRuler,
    ToggleCoordSystem,
    RenderShaded,
    RenderWireframe,
    FitAll,
    ZoomIn,
    ZoomOut,
    // Simulation
    Validate,
    AddSetup,
    AddSweep,
    Solve,
    SolveAll,
    Abort,
    // Draw
    DrawBox,
    DrawCylinder,
    DrawSphere,
    DrawCone,
    DrawTorus,
    DrawRectangle,
    DrawEllipse,
    DrawCircle,
    DrawPolygon,
    DrawPolyline,
    DrawArc,
    DrawSpline,
    SetPlane,
    SetUnits,
    AssignMaterial,
    // Model
    BoolUnite,
    BoolSubtract,
    BoolIntersect,
    BoolSplit,
    GroupObjects,
    UngroupObjects,
    AssignColor,
    SetTransparency,
    // Results
    CreateReport,
    SolutionData,
    PlotFields,
    PlotEField,
    PlotHField,
    PlotSAR,
    Animate,
    // View layout toggles
    ToggleProjectManager,
    ToggleMessageManager,
}

// ---------------------------------------------------------------------------
// Data model
// ---------------------------------------------------------------------------

pub struct RibbonState {
    pub active_tab: usize,
    pub toggles: HashMap<String, bool>,
    open_popup: Option<Id>,
}

impl Default for RibbonState {
    fn default() -> Self {
        Self {
            active_tab: 0,
            toggles: HashMap::new(),
            open_popup: None,
        }
    }
}

impl RibbonState {
    fn is_toggled(&self, key: &str) -> bool {
        self.toggles.get(key).copied().unwrap_or(false)
    }

    fn toggle(&mut self, key: &str) {
        let v = self.toggles.entry(key.to_string()).or_insert(false);
        *v = !*v;
    }
}

#[derive(Clone)]
pub struct MenuItem {
    pub label: String,
    pub action: RibbonAction,
    pub enabled: bool,
}

#[derive(Clone)]
pub enum RibbonItem {
    LargeButton {
        label: String,
        icon: char,
        action: RibbonAction,
        enabled: bool,
        tooltip: String,
    },
    SmallButton {
        label: String,
        icon: char,
        action: RibbonAction,
        enabled: bool,
        tooltip: String,
    },
    LargeSplitButton {
        label: String,
        icon: char,
        action: RibbonAction,
        enabled: bool,
        tooltip: String,
        menu_items: Vec<MenuItem>,
    },
    LargeDropdown {
        label: String,
        icon: char,
        enabled: bool,
        tooltip: String,
        menu_items: Vec<MenuItem>,
    },
    ToggleButton {
        label: String,
        icon: char,
        state_key: String,
        action: RibbonAction,
        tooltip: String,
    },
    ComboBox {
        label: String,
        options: Vec<String>,
        selected: usize,
        action: RibbonAction,
    },
}

#[derive(Clone)]
pub struct RibbonSubGroup {
    pub items: Vec<RibbonItem>,
}

#[derive(Clone)]
pub struct RibbonGroup {
    pub label: String,
    pub sub_groups: Vec<RibbonSubGroup>,
}

#[derive(Clone)]
pub struct RibbonTab {
    pub label: String,
    pub accent_color: Option<Color32>,
    pub groups: Vec<RibbonGroup>,
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

pub fn show_ribbon(ui: &mut Ui, state: &mut RibbonState, tabs: &[RibbonTab]) -> Option<RibbonAction> {
    let mut action = None;

    // Background
    let painter = ui.painter();
    let rect = ui.available_rect_before_wrap();
    painter.rect_filled(
        Rect::from_min_size(rect.min, vec2(rect.width(), TAB_STRIP_HEIGHT + COMMAND_AREA_HEIGHT + GROUP_LABEL_HEIGHT + 4.0)),
        CornerRadius::ZERO,
        RIBBON_BG,
    );

    ui.vertical(|ui| {
        ui.set_min_height(TAB_STRIP_HEIGHT + COMMAND_AREA_HEIGHT + GROUP_LABEL_HEIGHT);

        // Tab strip
        action = render_tab_strip(ui, state, tabs);

        // Thin separator between tab strip and command area
        let sep_rect = Rect::from_min_size(
            ui.cursor().min,
            vec2(ui.available_width(), 1.0),
        );
        ui.painter().rect_filled(sep_rect, CornerRadius::ZERO, GROUP_SEP_COLOR);
        ui.advance_cursor_after_rect(sep_rect);

        // Command area
        if let Some(tab) = tabs.get(state.active_tab) {
            let cmd_action = render_command_area(ui, state, tab);
            if action.is_none() {
                action = cmd_action;
            }
        }
    });

    action
}

// ---------------------------------------------------------------------------
// Tab strip
// ---------------------------------------------------------------------------

fn render_tab_strip(ui: &mut Ui, state: &mut RibbonState, tabs: &[RibbonTab]) -> Option<RibbonAction> {
    ui.horizontal(|ui| {
        ui.set_height(TAB_STRIP_HEIGHT);
        ui.spacing_mut().item_spacing.x = 0.0;

        for (i, tab) in tabs.iter().enumerate() {
            let is_active = state.active_tab == i;

            // Contextual tab accent bar
            if let Some(accent) = tab.accent_color {
                let tab_rect = Rect::from_min_size(
                    ui.cursor().min,
                    vec2(60.0, 3.0),
                );
                ui.painter().rect_filled(tab_rect, CornerRadius::ZERO, accent);
            }

            let text_color = if is_active { TAB_ACTIVE_TEXT } else { TAB_NORMAL_TEXT };
            let font = if is_active {
                FontId::proportional(12.0)
            } else {
                FontId::proportional(11.5)
            };

            let text = RichText::new(&tab.label).font(font.clone()).color(text_color);
            let padding = vec2(12.0, 0.0);

            let (rect, response) = ui.allocate_exact_size(
                vec2(text.text().len() as f32 * 7.5 + padding.x * 2.0, TAB_STRIP_HEIGHT),
                Sense::click(),
            );

            // Paint background
            if is_active {
                ui.painter().rect_filled(
                    rect,
                    CornerRadius::same(2),
                    RIBBON_BG,
                );
                // Bottom edge blends with command area (no border)
            } else if response.hovered() {
                ui.painter().rect_filled(
                    rect,
                    CornerRadius::same(2),
                    Color32::from_rgb(230, 230, 230),
                );
            }

            // Paint text centered
            let galley = ui.painter().layout_no_wrap(
                tab.label.clone(),
                font.clone(),
                text_color,
            );
            let text_pos = rect.center() - galley.size() / 2.0;
            ui.painter().galley(text_pos, galley, text_color);

            if response.clicked() {
                state.active_tab = i;
            }
        }
    })
    .inner;
    None
}

// ---------------------------------------------------------------------------
// Command area
// ---------------------------------------------------------------------------

fn render_command_area(ui: &mut Ui, state: &mut RibbonState, tab: &RibbonTab) -> Option<RibbonAction> {
    let mut action = None;

    ui.horizontal(|ui| {
        ui.set_min_height(COMMAND_AREA_HEIGHT + GROUP_LABEL_HEIGHT);
        ui.spacing_mut().item_spacing.x = 0.0;

        for (gi, group) in tab.groups.iter().enumerate() {
            ui.add_space(GROUP_H_PAD);

            let group_response = ui.vertical(|ui| {
                ui.set_min_height(COMMAND_AREA_HEIGHT + GROUP_LABEL_HEIGHT);

                // Button content area: sub_groups separated by vertical lines
                ui.with_layout(Layout::left_to_right(Align::Min), |ui| {
                    ui.set_min_height(COMMAND_AREA_HEIGHT);
                    ui.spacing_mut().item_spacing.x = 2.0;

                    for (si, sub) in group.sub_groups.iter().enumerate() {
                        let a = render_group_items(ui, state, &sub.items);
                        if action.is_none() {
                            action = a;
                        }

                        // Vertical separator between sub-groups
                        if si < group.sub_groups.len() - 1 {
                            ui.add_space(3.0);
                            let h = COMMAND_AREA_HEIGHT;
                            let (sep_rect, _) = ui.allocate_exact_size(vec2(1.0, h), Sense::hover());
                            ui.painter().rect_filled(
                                Rect::from_center_size(
                                    sep_rect.center(),
                                    vec2(1.0, h - 8.0),
                                ),
                                CornerRadius::ZERO,
                                GROUP_SEP_COLOR,
                            );
                            ui.add_space(3.0);
                        }
                    }
                });

                // Group label at bottom, centered
                ui.with_layout(Layout::top_down(Align::Center), |ui| {
                    ui.set_height(GROUP_LABEL_HEIGHT);
                    ui.label(
                        RichText::new(&group.label)
                            .font(FontId::proportional(GROUP_LABEL_FONT))
                            .color(GROUP_LABEL_COLOR),
                    );
                });
            });

            ui.add_space(GROUP_H_PAD);

            // Draw vertical separator line between groups (not after last)
            if gi < tab.groups.len() - 1 {
                let grect = group_response.response.rect;
                // Allocate a narrow rect for the separator and paint it filled
                let (sep_rect, _) = ui.allocate_exact_size(
                    vec2(1.0, grect.height()),
                    Sense::hover(),
                );
                ui.painter().rect_filled(
                    Rect::from_center_size(
                        sep_rect.center(),
                        vec2(1.0, grect.height() - 8.0),
                    ),
                    CornerRadius::ZERO,
                    GROUP_SEP_COLOR,
                );
            }
        }
    });

    action
}

// ---------------------------------------------------------------------------
// Group items rendering
// ---------------------------------------------------------------------------

fn render_group_items(ui: &mut Ui, state: &mut RibbonState, items: &[RibbonItem]) -> Option<RibbonAction> {
    let mut action = None;

    // Partition items: render large items directly, stack small items in groups of 3
    let mut small_batch: Vec<&RibbonItem> = Vec::new();

    for item in items {
        match item {
            RibbonItem::LargeButton { .. }
            | RibbonItem::LargeSplitButton { .. }
            | RibbonItem::LargeDropdown { .. } => {
                // Flush small batch first
                if !small_batch.is_empty() {
                    let a = render_small_stack(ui, state, &small_batch);
                    if action.is_none() { action = a; }
                    small_batch.clear();
                }
                let a = render_item(ui, state, item);
                if action.is_none() { action = a; }
            }
            RibbonItem::SmallButton { .. }
            | RibbonItem::ToggleButton { .. } => {
                small_batch.push(item);
                if small_batch.len() == 3 {
                    let a = render_small_stack(ui, state, &small_batch);
                    if action.is_none() { action = a; }
                    small_batch.clear();
                }
            }
            RibbonItem::ComboBox { .. } => {
                // Flush and render inline
                if !small_batch.is_empty() {
                    let a = render_small_stack(ui, state, &small_batch);
                    if action.is_none() { action = a; }
                    small_batch.clear();
                }
                let a = render_item(ui, state, item);
                if action.is_none() { action = a; }
            }
        }
    }

    // Flush remaining small items
    if !small_batch.is_empty() {
        let a = render_small_stack(ui, state, &small_batch);
        if action.is_none() { action = a; }
    }

    action
}

fn render_small_stack(ui: &mut Ui, state: &mut RibbonState, items: &[&RibbonItem]) -> Option<RibbonAction> {
    let mut action = None;

    ui.vertical(|ui| {
        ui.set_min_height(COMMAND_AREA_HEIGHT);
        ui.spacing_mut().item_spacing.y = 1.0;

        for item in items {
            let a = render_item(ui, state, item);
            if action.is_none() { action = a; }
        }
    });

    action
}

// ---------------------------------------------------------------------------
// Individual item rendering
// ---------------------------------------------------------------------------

fn render_item(ui: &mut Ui, state: &mut RibbonState, item: &RibbonItem) -> Option<RibbonAction> {
    match item {
        RibbonItem::LargeButton { label, icon, action, enabled, tooltip } => {
            render_large_button(ui, label, *icon, *action, *enabled, tooltip)
        }
        RibbonItem::SmallButton { label, icon, action, enabled, tooltip } => {
            render_small_button(ui, label, *icon, *action, *enabled, tooltip)
        }
        RibbonItem::LargeSplitButton { label, icon, action, enabled, tooltip, menu_items } => {
            render_large_split(ui, state, label, *icon, *action, *enabled, tooltip, menu_items)
        }
        RibbonItem::LargeDropdown { label, icon, enabled, tooltip, menu_items } => {
            render_large_dropdown(ui, state, label, *icon, *enabled, tooltip, menu_items)
        }
        RibbonItem::ToggleButton { label, icon, state_key, action, tooltip } => {
            render_toggle_button(ui, state, label, *icon, state_key, *action, tooltip)
        }
        RibbonItem::ComboBox { label, options, selected, action } => {
            render_combo_box(ui, label, options, *selected, *action)
        }
    }
}

// ---------------------------------------------------------------------------
// Large button
// ---------------------------------------------------------------------------

fn render_large_button(
    ui: &mut Ui,
    label: &str,
    icon: char,
    action: RibbonAction,
    enabled: bool,
    tooltip: &str,
) -> Option<RibbonAction> {
    let mut result = None;

    let (rect, response) = ui.allocate_exact_size(
        vec2(LARGE_BTN_WIDTH, LARGE_BTN_HEIGHT),
        if enabled { Sense::click() } else { Sense::hover() },
    );

    paint_button_bg(ui, &response, rect, false, enabled);

    let text_color = if enabled { Color32::from_rgb(30, 30, 30) } else { DISABLED_TEXT };

    // Icon
    let icon_galley = ui.painter().layout_no_wrap(
        icon.to_string(),
        FontId::proportional(ICON_FONT_LARGE),
        text_color,
    );
    let icon_pos = egui::pos2(
        rect.center().x - icon_galley.size().x / 2.0,
        rect.min.y + 4.0,
    );
    ui.painter().galley(icon_pos, icon_galley, text_color);

    // Label (up to 2 lines, centered)
    let label_galley = ui.painter().layout(
        label.to_string(),
        FontId::proportional(LABEL_FONT),
        text_color,
        LARGE_BTN_WIDTH - 2.0,
    );
    let label_pos = egui::pos2(
        rect.center().x - label_galley.size().x / 2.0,
        rect.min.y + 4.0 + ICON_FONT_LARGE + 4.0,
    );
    ui.painter().galley(label_pos, label_galley, text_color);

    if enabled && response.clicked() {
        result = Some(action);
    }

    if !tooltip.is_empty() {
        response.on_hover_text(tooltip);
    }

    result
}

// ---------------------------------------------------------------------------
// Small button
// ---------------------------------------------------------------------------

fn render_small_button(
    ui: &mut Ui,
    label: &str,
    icon: char,
    action: RibbonAction,
    enabled: bool,
    tooltip: &str,
) -> Option<RibbonAction> {
    let mut result = None;

    let desired_width = ICON_FONT_SMALL + 4.0 + label.len() as f32 * 6.5 + 8.0;
    let (rect, response) = ui.allocate_exact_size(
        vec2(desired_width.max(60.0), SMALL_ROW_HEIGHT),
        if enabled { Sense::click() } else { Sense::hover() },
    );

    paint_button_bg(ui, &response, rect, false, enabled);

    let text_color = if enabled { Color32::from_rgb(30, 30, 30) } else { DISABLED_TEXT };

    // Icon
    let icon_galley = ui.painter().layout_no_wrap(
        icon.to_string(),
        FontId::proportional(ICON_FONT_SMALL),
        text_color,
    );
    let icon_pos = egui::pos2(rect.min.x + 4.0, rect.center().y - icon_galley.size().y / 2.0);
    ui.painter().galley(icon_pos, icon_galley, text_color);

    // Label
    let label_galley = ui.painter().layout_no_wrap(
        label.to_string(),
        FontId::proportional(LABEL_FONT),
        text_color,
    );
    let label_pos = egui::pos2(
        rect.min.x + 4.0 + ICON_FONT_SMALL + 4.0,
        rect.center().y - label_galley.size().y / 2.0,
    );
    ui.painter().galley(label_pos, label_galley, text_color);

    if enabled && response.clicked() {
        result = Some(action);
    }

    if !tooltip.is_empty() {
        response.on_hover_text(tooltip);
    }

    result
}

// ---------------------------------------------------------------------------
// Large split button
// ---------------------------------------------------------------------------

fn render_large_split(
    ui: &mut Ui,
    state: &mut RibbonState,
    label: &str,
    icon: char,
    action: RibbonAction,
    enabled: bool,
    tooltip: &str,
    menu_items: &[MenuItem],
) -> Option<RibbonAction> {
    let mut result = None;
    let btn_id = Id::new(("split", label));

    let total_rect_size = vec2(LARGE_BTN_WIDTH, LARGE_BTN_HEIGHT);
    let (total_rect, _) = ui.allocate_exact_size(total_rect_size, Sense::hover());

    // Top zone: icon area (click = default action)
    let top_rect = Rect::from_min_size(total_rect.min, vec2(LARGE_BTN_WIDTH, 32.0));
    let top_response = ui.interact(top_rect, btn_id.with("top"), if enabled { Sense::click() } else { Sense::hover() });

    // Bottom zone: label + arrow (click = open menu)
    let bot_rect = Rect::from_min_max(
        egui::pos2(total_rect.min.x, total_rect.min.y + 32.0),
        total_rect.max,
    );
    let bot_response = ui.interact(bot_rect, btn_id.with("bot"), if enabled { Sense::click() } else { Sense::hover() });

    // Paint backgrounds
    paint_button_bg(ui, &top_response, top_rect, false, enabled);
    paint_button_bg(ui, &bot_response, bot_rect, false, enabled);

    // Separator line between zones (on hover)
    if top_response.hovered() || bot_response.hovered() {
        ui.painter().hline(
            top_rect.min.x + 2.0..=top_rect.max.x - 2.0,
            bot_rect.min.y,
            Stroke::new(1.0, HOVER_STROKE),
        );
    }

    let text_color = if enabled { Color32::from_rgb(30, 30, 30) } else { DISABLED_TEXT };

    // Icon
    let icon_galley = ui.painter().layout_no_wrap(
        icon.to_string(),
        FontId::proportional(ICON_FONT_LARGE),
        text_color,
    );
    let icon_pos = egui::pos2(
        top_rect.center().x - icon_galley.size().x / 2.0,
        top_rect.min.y + 4.0,
    );
    ui.painter().galley(icon_pos, icon_galley, text_color);

    // Label + arrow
    let arrow_text = format!("{} \u{25BC}", label);
    let label_galley = ui.painter().layout(
        arrow_text,
        FontId::proportional(LABEL_FONT),
        text_color,
        LARGE_BTN_WIDTH - 2.0,
    );
    let label_pos = egui::pos2(
        bot_rect.center().x - label_galley.size().x / 2.0,
        bot_rect.min.y + 2.0,
    );
    ui.painter().galley(label_pos, label_galley, text_color);

    if enabled && top_response.clicked() {
        result = Some(action);
    }

    if enabled && bot_response.clicked() {
        let is_open = state.open_popup == Some(btn_id);
        state.open_popup = if is_open { None } else { Some(btn_id) };
    }

    let top_hovered = top_response.hovered();
    let bot_hovered = bot_response.hovered();

    if !tooltip.is_empty() {
        top_response.on_hover_text(tooltip);
    }

    // Popup menu
    if state.open_popup == Some(btn_id) {
        let menu_action = show_popup_menu(ui, btn_id, total_rect, menu_items);
        if let Some(a) = menu_action {
            result = Some(a);
            state.open_popup = None;
        }
        // Close if clicked elsewhere
        if ui.input(|i| i.pointer.any_click()) && !top_hovered && !bot_hovered {
            // Will close next frame if no menu item clicked
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Large dropdown
// ---------------------------------------------------------------------------

fn render_large_dropdown(
    ui: &mut Ui,
    state: &mut RibbonState,
    label: &str,
    icon: char,
    enabled: bool,
    tooltip: &str,
    menu_items: &[MenuItem],
) -> Option<RibbonAction> {
    let mut result = None;
    let btn_id = Id::new(("dropdown", label));

    let (rect, response) = ui.allocate_exact_size(
        vec2(LARGE_BTN_WIDTH, LARGE_BTN_HEIGHT),
        if enabled { Sense::click() } else { Sense::hover() },
    );

    paint_button_bg(ui, &response, rect, false, enabled);

    let text_color = if enabled { Color32::from_rgb(30, 30, 30) } else { DISABLED_TEXT };

    // Icon
    let icon_galley = ui.painter().layout_no_wrap(
        icon.to_string(),
        FontId::proportional(ICON_FONT_LARGE),
        text_color,
    );
    let icon_pos = egui::pos2(
        rect.center().x - icon_galley.size().x / 2.0,
        rect.min.y + 4.0,
    );
    ui.painter().galley(icon_pos, icon_galley, text_color);

    // Label + arrow
    let arrow_text = format!("{} \u{25BC}", label);
    let label_galley = ui.painter().layout(
        arrow_text,
        FontId::proportional(LABEL_FONT),
        text_color,
        LARGE_BTN_WIDTH - 2.0,
    );
    let label_pos = egui::pos2(
        rect.center().x - label_galley.size().x / 2.0,
        rect.min.y + 4.0 + ICON_FONT_LARGE + 4.0,
    );
    ui.painter().galley(label_pos, label_galley, text_color);

    if enabled && response.clicked() {
        let is_open = state.open_popup == Some(btn_id);
        state.open_popup = if is_open { None } else { Some(btn_id) };
    }

    if !tooltip.is_empty() {
        response.on_hover_text(tooltip);
    }

    if state.open_popup == Some(btn_id) {
        let menu_action = show_popup_menu(ui, btn_id, rect, menu_items);
        if let Some(a) = menu_action {
            result = Some(a);
            state.open_popup = None;
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Toggle button (small)
// ---------------------------------------------------------------------------

fn render_toggle_button(
    ui: &mut Ui,
    state: &mut RibbonState,
    label: &str,
    icon: char,
    state_key: &str,
    action: RibbonAction,
    tooltip: &str,
) -> Option<RibbonAction> {
    let mut result = None;
    let checked = state.is_toggled(state_key);

    let desired_width = ICON_FONT_SMALL + 4.0 + label.len() as f32 * 6.5 + 8.0;
    let (rect, response) = ui.allocate_exact_size(
        vec2(desired_width.max(60.0), SMALL_ROW_HEIGHT),
        Sense::click(),
    );

    paint_button_bg(ui, &response, rect, checked, true);

    let text_color = Color32::from_rgb(30, 30, 30);

    // Icon
    let icon_galley = ui.painter().layout_no_wrap(
        icon.to_string(),
        FontId::proportional(ICON_FONT_SMALL),
        text_color,
    );
    let icon_pos = egui::pos2(rect.min.x + 4.0, rect.center().y - icon_galley.size().y / 2.0);
    ui.painter().galley(icon_pos, icon_galley, text_color);

    // Label
    let label_galley = ui.painter().layout_no_wrap(
        label.to_string(),
        FontId::proportional(LABEL_FONT),
        text_color,
    );
    let label_pos = egui::pos2(
        rect.min.x + 4.0 + ICON_FONT_SMALL + 4.0,
        rect.center().y - label_galley.size().y / 2.0,
    );
    ui.painter().galley(label_pos, label_galley, text_color);

    if response.clicked() {
        state.toggle(state_key);
        result = Some(action);
    }

    if !tooltip.is_empty() {
        response.on_hover_text(tooltip);
    }

    result
}

// ---------------------------------------------------------------------------
// Combo box
// ---------------------------------------------------------------------------

fn render_combo_box(
    ui: &mut Ui,
    label: &str,
    options: &[String],
    selected: usize,
    _action: RibbonAction,
) -> Option<RibbonAction> {
    ui.vertical(|ui| {
        ui.set_min_height(COMMAND_AREA_HEIGHT);
        ui.add_space(4.0);
        ui.label(
            RichText::new(label)
                .font(FontId::proportional(LABEL_FONT))
                .color(GROUP_LABEL_COLOR),
        );

        let current = options.get(selected).map(|s| s.as_str()).unwrap_or("");
        let mut _sel = selected;
        egui::ComboBox::from_id_salt(label)
            .width(90.0)
            .selected_text(current)
            .show_ui(ui, |ui| {
                for (i, opt) in options.iter().enumerate() {
                    ui.selectable_value(&mut _sel, i, opt);
                }
            });
    });

    None // ComboBox action handled separately via state
}

// ---------------------------------------------------------------------------
// Popup menu
// ---------------------------------------------------------------------------

fn show_popup_menu(
    ui: &mut Ui,
    popup_id: Id,
    anchor_rect: Rect,
    items: &[MenuItem],
) -> Option<RibbonAction> {
    let mut result = None;

    let area_resp = egui::Area::new(popup_id.with("menu"))
        .order(egui::Order::Foreground)
        .fixed_pos(egui::pos2(anchor_rect.min.x, anchor_rect.max.y + 1.0))
        .show(ui.ctx(), |ui| {
            egui::Frame::popup(ui.style()).show(ui, |ui| {
                ui.set_min_width(140.0);
                for item in items {
                    let text_color = if item.enabled {
                        Color32::from_rgb(30, 30, 30)
                    } else {
                        DISABLED_TEXT
                    };

                    let response = ui.add_sized(
                        vec2(140.0, 24.0),
                        egui::Label::new(
                            RichText::new(&item.label)
                                .font(FontId::proportional(LABEL_FONT))
                                .color(text_color),
                        )
                        .sense(if item.enabled { Sense::click() } else { Sense::hover() }),
                    );

                    // Hover highlight
                    if item.enabled && response.hovered() {
                        ui.painter().rect_filled(
                            response.rect,
                            CornerRadius::same(2),
                            MENU_HOVER,
                        );
                    }

                    if item.enabled && response.clicked() {
                        result = Some(item.action);
                    }
                }
            });
        });

    // Close popup on click outside
    if ui.input(|i| i.pointer.any_click()) && result.is_none() {
        if let Some(pos) = ui.input(|i| i.pointer.latest_pos()) {
            if !area_resp.response.rect.contains(pos) && !anchor_rect.contains(pos) {
                // Signal close by returning None — caller checks open_popup
            }
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Shared paint helpers
// ---------------------------------------------------------------------------

fn paint_button_bg(ui: &Ui, response: &Response, rect: Rect, checked: bool, enabled: bool) {
    if !enabled {
        return;
    }

    if checked {
        ui.painter().rect(
            rect,
            CornerRadius::same(2),
            CHECKED_FILL,
            Stroke::new(1.0, CHECKED_STROKE),
            StrokeKind::Outside,
        );
    } else if response.is_pointer_button_down_on() {
        ui.painter().rect(
            rect,
            CornerRadius::same(2),
            PRESSED_FILL,
            Stroke::new(1.0, HOVER_STROKE),
            StrokeKind::Outside,
        );
    } else if response.hovered() {
        ui.painter().rect(
            rect,
            CornerRadius::same(2),
            HOVER_FILL,
            Stroke::new(1.0, HOVER_STROKE),
            StrokeKind::Outside,
        );
    }
}

// ---------------------------------------------------------------------------
// Tab definition builders (convenience)
// ---------------------------------------------------------------------------

use emstudio_domain::Edition;

pub fn build_default_tabs(edition: Edition, is_web: bool) -> Vec<RibbonTab> {
    vec![
        build_desktop_tab(is_web),
        build_view_tab(),
        build_simulation_tab(edition),
        build_draw_tab(),
        build_model_tab(),
        build_results_tab(),
    ]
}

fn build_desktop_tab(is_web: bool) -> RibbonTab {
    let file_ops = !is_web;
    RibbonTab {
        label: "Desktop".into(),
        accent_color: None,
        groups: vec![
            RibbonGroup {
                label: "Project".into(),
                sub_groups: vec![
                    RibbonSubGroup {
                        items: vec![
                            RibbonItem::LargeButton {
                                label: "New".into(),
                                icon: '\u{1F4C4}',
                                action: RibbonAction::NewProject,
                                enabled: file_ops,
                                tooltip: if is_web { "Not available in web version".into() } else { "Create a new project".into() },
                            },
                            RibbonItem::LargeButton {
                                label: "Open".into(),
                                icon: '\u{1F4C2}',
                                action: RibbonAction::OpenProject,
                                enabled: file_ops,
                                tooltip: if is_web { "Not available in web version".into() } else { "Open an existing project".into() },
                            },
                            RibbonItem::LargeButton {
                                label: "Save".into(),
                                icon: '\u{1F4BE}',
                                action: RibbonAction::SaveProject,
                                enabled: file_ops,
                                tooltip: if is_web { "Not available in web version".into() } else { "Save current project".into() },
                            },
                        ],
                    },
                    RibbonSubGroup {
                        items: vec![
                            RibbonItem::SmallButton {
                                label: "Save As".into(),
                                icon: '\u{1F4CB}',
                                action: RibbonAction::SaveAs,
                                enabled: file_ops,
                                tooltip: if is_web { "Not available in web version".into() } else { "Save project as a new file".into() },
                            },
                            RibbonItem::SmallButton {
                                label: "Close".into(),
                                icon: '\u{2716}',
                                action: RibbonAction::CloseProject,
                                enabled: true,
                                tooltip: "Close current project".into(),
                            },
                        ],
                    },
                ],
            },
            RibbonGroup {
                label: "Import / Export".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::LargeSplitButton {
                            label: "Import".into(),
                            icon: '\u{1F4E5}',
                            action: RibbonAction::ImportStep,
                            enabled: file_ops,
                            tooltip: if is_web { "Not available in web version".into() } else { "Import geometry".into() },
                            menu_items: vec![
                                MenuItem { label: "Import STEP".into(), action: RibbonAction::ImportStep, enabled: file_ops },
                                MenuItem { label: "Import SAT".into(), action: RibbonAction::ImportSat, enabled: file_ops },
                                MenuItem { label: "Import HFSS (.aedt/.py)".into(), action: RibbonAction::ImportHfssAedt, enabled: file_ops },
                                MenuItem { label: "Import Q3D (.aedt/.py)".into(), action: RibbonAction::ImportQ3dAedt, enabled: file_ops },
                            ],
                        },
                        RibbonItem::LargeSplitButton {
                            label: "Export".into(),
                            icon: '\u{1F4E4}',
                            action: RibbonAction::ExportStep,
                            enabled: file_ops,
                            tooltip: if is_web { "Not available in web version".into() } else { "Export geometry".into() },
                            menu_items: vec![
                                MenuItem { label: "Export STEP".into(), action: RibbonAction::ExportStep, enabled: file_ops },
                                MenuItem { label: "Export SAT".into(), action: RibbonAction::ExportSat, enabled: file_ops },
                                MenuItem { label: "Export HFSS Script (.py)".into(), action: RibbonAction::ExportHfssPyAedt, enabled: file_ops },
                                MenuItem { label: "Export Q3D Script (.py)".into(), action: RibbonAction::ExportQ3dPyAedt, enabled: file_ops },
                            ],
                        },
                    ],
                }],
            },
        ],
    }
}

fn build_view_tab() -> RibbonTab {
    RibbonTab {
        label: "View".into(),
        accent_color: None,
        groups: vec![
            RibbonGroup {
                label: "Visibility".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::ToggleButton {
                            label: "Grid".into(),
                            icon: '\u{2630}',
                            state_key: "grid".into(),
                            action: RibbonAction::ToggleGrid,
                            tooltip: "Toggle grid display".into(),
                        },
                        RibbonItem::ToggleButton {
                            label: "Ruler".into(),
                            icon: '\u{1F4CF}',
                            state_key: "ruler".into(),
                            action: RibbonAction::ToggleRuler,
                            tooltip: "Toggle ruler display".into(),
                        },
                        RibbonItem::ToggleButton {
                            label: "Axes".into(),
                            icon: '\u{2742}',
                            state_key: "coord_system".into(),
                            action: RibbonAction::ToggleCoordSystem,
                            tooltip: "Toggle coordinate system display".into(),
                        },
                    ],
                }],
            },
            RibbonGroup {
                label: "Render".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::ToggleButton {
                            label: "Shaded".into(),
                            icon: '\u{25A3}',
                            state_key: "shaded".into(),
                            action: RibbonAction::RenderShaded,
                            tooltip: "Shaded render mode".into(),
                        },
                        RibbonItem::ToggleButton {
                            label: "Wireframe".into(),
                            icon: '\u{25A1}',
                            state_key: "wireframe".into(),
                            action: RibbonAction::RenderWireframe,
                            tooltip: "Wireframe render mode".into(),
                        },
                    ],
                }],
            },
            RibbonGroup {
                label: "Zoom".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::LargeButton {
                            label: "Fit All".into(),
                            icon: '\u{26F6}',
                            action: RibbonAction::FitAll,
                            enabled: true,
                            tooltip: "Fit all objects in view".into(),
                        },
                        RibbonItem::SmallButton {
                            label: "Zoom In".into(),
                            icon: '\u{1F50D}',
                            action: RibbonAction::ZoomIn,
                            enabled: true,
                            tooltip: "Zoom in".into(),
                        },
                        RibbonItem::SmallButton {
                            label: "Zoom Out".into(),
                            icon: '\u{1F50E}',
                            action: RibbonAction::ZoomOut,
                            enabled: true,
                            tooltip: "Zoom out".into(),
                        },
                    ],
                }],
            },
        ],
    }
}

fn build_simulation_tab(edition: Edition) -> RibbonTab {
    RibbonTab {
        label: "Simulation".into(),
        accent_color: None,
        groups: vec![
            RibbonGroup {
                label: "Validate".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::LargeButton {
                            label: "Validate".into(),
                            icon: '\u{2714}',
                            action: RibbonAction::Validate,
                            enabled: true,
                            tooltip: "Validate current design".into(),
                        },
                    ],
                }],
            },
            RibbonGroup {
                label: "Analysis".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::LargeButton {
                            label: "Add\nSetup".into(),
                            icon: '\u{2699}',
                            action: RibbonAction::AddSetup,
                            enabled: true,
                            tooltip: "Add solution setup".into(),
                        },
                        RibbonItem::SmallButton {
                            label: "Add Sweep".into(),
                            icon: '\u{27A1}',
                            action: RibbonAction::AddSweep,
                            enabled: true,
                            tooltip: "Add frequency sweep".into(),
                        },
                    ],
                }],
            },
            RibbonGroup {
                label: "Solve".into(),
                sub_groups: vec![
                    RibbonSubGroup {
                        items: vec![
                            RibbonItem::LargeButton {
                                label: "Analyze\nAll".into(),
                                icon: '\u{25B6}',
                                action: RibbonAction::SolveAll,
                                enabled: edition.allows_solve_all(),
                                tooltip: if edition.allows_solve_all() {
                                    "Analyze all setups".into()
                                } else {
                                    "Requires Professional edition".into()
                                },
                            },
                            RibbonItem::LargeButton {
                                label: "Solve".into(),
                                icon: '\u{23EF}',
                                action: RibbonAction::Solve,
                                enabled: true,
                                tooltip: "Solve current setup".into(),
                            },
                        ],
                    },
                    RibbonSubGroup {
                        items: vec![
                            RibbonItem::SmallButton {
                                label: "Abort".into(),
                                icon: '\u{23F9}',
                                action: RibbonAction::Abort,
                                enabled: false,
                                tooltip: "Abort current simulation".into(),
                            },
                        ],
                    },
                ],
            },
        ],
    }
}

fn build_draw_tab() -> RibbonTab {
    RibbonTab {
        label: "Draw".into(),
        accent_color: Some(Color32::from_rgb(0, 120, 200)),
        groups: vec![
            RibbonGroup {
                label: "Primitives".into(),
                sub_groups: vec![
                    RibbonSubGroup {
                        items: vec![
                            RibbonItem::LargeSplitButton {
                                label: "Box".into(),
                                icon: '\u{25A3}',
                                action: RibbonAction::DrawBox,
                                enabled: true,
                                tooltip: "Draw a box".into(),
                                menu_items: vec![
                                    MenuItem { label: "Box".into(), action: RibbonAction::DrawBox, enabled: true },
                                    MenuItem { label: "Box by Corner".into(), action: RibbonAction::DrawBox, enabled: true },
                                ],
                            },
                            RibbonItem::LargeSplitButton {
                                label: "Cylinder".into(),
                                icon: '\u{25CB}',
                                action: RibbonAction::DrawCylinder,
                                enabled: true,
                                tooltip: "Draw a cylinder".into(),
                                menu_items: vec![
                                    MenuItem { label: "Cylinder".into(), action: RibbonAction::DrawCylinder, enabled: true },
                                    MenuItem { label: "Cylinder by Axis".into(), action: RibbonAction::DrawCylinder, enabled: true },
                                ],
                            },
                            RibbonItem::LargeButton {
                                label: "Sphere".into(),
                                icon: '\u{25CF}',
                                action: RibbonAction::DrawSphere,
                                enabled: true,
                                tooltip: "Draw a sphere".into(),
                            },
                            RibbonItem::LargeButton {
                                label: "Polyline".into(),
                                icon: '\u{2571}',
                                action: RibbonAction::DrawPolyline,
                                enabled: true,
                                tooltip: "Draw a polyline".into(),
                            },
                        ],
                    },
                    RibbonSubGroup {
                        items: vec![
                            RibbonItem::SmallButton {
                                label: "Cone".into(),
                                icon: '\u{25B3}',
                                action: RibbonAction::DrawCone,
                                enabled: true,
                                tooltip: "Draw a cone".into(),
                            },
                            RibbonItem::SmallButton {
                                label: "Torus".into(),
                                icon: '\u{25CE}',
                                action: RibbonAction::DrawTorus,
                                enabled: true,
                                tooltip: "Draw a torus".into(),
                            },
                            RibbonItem::SmallButton {
                                label: "Rectangle".into(),
                                icon: '\u{25AD}',
                                action: RibbonAction::DrawRectangle,
                                enabled: true,
                                tooltip: "Draw a rectangle".into(),
                            },
                            RibbonItem::SmallButton {
                                label: "Ellipse".into(),
                                icon: '\u{2B2D}',
                                action: RibbonAction::DrawEllipse,
                                enabled: true,
                                tooltip: "Draw an ellipse".into(),
                            },
                            RibbonItem::SmallButton {
                                label: "Circle".into(),
                                icon: '\u{25EF}',
                                action: RibbonAction::DrawCircle,
                                enabled: true,
                                tooltip: "Draw a circle".into(),
                            },
                            RibbonItem::SmallButton {
                                label: "Polygon".into(),
                                icon: '\u{2B23}',
                                action: RibbonAction::DrawPolygon,
                                enabled: true,
                                tooltip: "Draw a regular polygon".into(),
                            },
                            RibbonItem::SmallButton {
                                label: "Arc".into(),
                                icon: '\u{25DC}',
                                action: RibbonAction::DrawArc,
                                enabled: true,
                                tooltip: "Draw an arc".into(),
                            },
                            RibbonItem::SmallButton {
                                label: "Spline".into(),
                                icon: '\u{223F}',
                                action: RibbonAction::DrawSpline,
                                enabled: true,
                                tooltip: "Draw a spline".into(),
                            },
                        ],
                    },
                ],
            },
            RibbonGroup {
                label: "Plane".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::ComboBox {
                            label: "Plane".into(),
                            options: vec!["XY".into(), "YZ".into(), "XZ".into()],
                            selected: 0,
                            action: RibbonAction::SetPlane,
                        },
                    ],
                }],
            },
            RibbonGroup {
                label: "Units".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::LargeDropdown {
                            label: "Units".into(),
                            icon: '\u{1F4D0}',
                            enabled: true,
                            tooltip: "Set model units".into(),
                            menu_items: vec![
                                MenuItem { label: "mm".into(), action: RibbonAction::SetUnits, enabled: true },
                                MenuItem { label: "cm".into(), action: RibbonAction::SetUnits, enabled: true },
                                MenuItem { label: "m".into(), action: RibbonAction::SetUnits, enabled: true },
                                MenuItem { label: "mil".into(), action: RibbonAction::SetUnits, enabled: true },
                                MenuItem { label: "in".into(), action: RibbonAction::SetUnits, enabled: true },
                            ],
                        },
                    ],
                }],
            },
            RibbonGroup {
                label: "Material".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::LargeButton {
                            label: "Assign\nMaterial".into(),
                            icon: '\u{1F9F1}',
                            action: RibbonAction::AssignMaterial,
                            enabled: true,
                            tooltip: "Assign material to selected object".into(),
                        },
                    ],
                }],
            },
        ],
    }
}

fn build_model_tab() -> RibbonTab {
    RibbonTab {
        label: "Model".into(),
        accent_color: Some(Color32::from_rgb(0, 120, 200)),
        groups: vec![
            RibbonGroup {
                label: "Boolean".into(),
                sub_groups: vec![
                    RibbonSubGroup {
                        items: vec![
                            RibbonItem::LargeButton {
                                label: "Unite".into(),
                                icon: '\u{222A}',
                                action: RibbonAction::BoolUnite,
                                enabled: true,
                                tooltip: "Unite selected objects".into(),
                            },
                            RibbonItem::LargeButton {
                                label: "Subtract".into(),
                                icon: '\u{2212}',
                                action: RibbonAction::BoolSubtract,
                                enabled: true,
                                tooltip: "Subtract objects".into(),
                            },
                            RibbonItem::LargeButton {
                                label: "Intersect".into(),
                                icon: '\u{2229}',
                                action: RibbonAction::BoolIntersect,
                                enabled: true,
                                tooltip: "Intersect objects".into(),
                            },
                        ],
                    },
                    RibbonSubGroup {
                        items: vec![
                            RibbonItem::SmallButton {
                                label: "Split".into(),
                                icon: '\u{2702}',
                                action: RibbonAction::BoolSplit,
                                enabled: true,
                                tooltip: "Split object".into(),
                            },
                        ],
                    },
                ],
            },
            RibbonGroup {
                label: "Object".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::SmallButton {
                            label: "Group".into(),
                            icon: '\u{1F4E6}',
                            action: RibbonAction::GroupObjects,
                            enabled: true,
                            tooltip: "Group selected objects".into(),
                        },
                        RibbonItem::SmallButton {
                            label: "Ungroup".into(),
                            icon: '\u{1F4E4}',
                            action: RibbonAction::UngroupObjects,
                            enabled: true,
                            tooltip: "Ungroup selected objects".into(),
                        },
                        RibbonItem::SmallButton {
                            label: "Color".into(),
                            icon: '\u{1F3A8}',
                            action: RibbonAction::AssignColor,
                            enabled: true,
                            tooltip: "Assign color to object".into(),
                        },
                        RibbonItem::SmallButton {
                            label: "Transparency".into(),
                            icon: '\u{1F4A7}',
                            action: RibbonAction::SetTransparency,
                            enabled: true,
                            tooltip: "Set object transparency".into(),
                        },
                    ],
                }],
            },
        ],
    }
}

fn build_results_tab() -> RibbonTab {
    RibbonTab {
        label: "Results".into(),
        accent_color: Some(Color32::from_rgb(50, 160, 80)),
        groups: vec![
            RibbonGroup {
                label: "Report".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::LargeSplitButton {
                            label: "Create\nReport".into(),
                            icon: '\u{1F4CA}',
                            action: RibbonAction::CreateReport,
                            enabled: true,
                            tooltip: "Create a new report".into(),
                            menu_items: vec![
                                MenuItem { label: "Rectangular Plot".into(), action: RibbonAction::CreateReport, enabled: true },
                                MenuItem { label: "Smith Chart".into(), action: RibbonAction::CreateReport, enabled: true },
                                MenuItem { label: "Polar Plot".into(), action: RibbonAction::CreateReport, enabled: true },
                                MenuItem { label: "Data Table".into(), action: RibbonAction::CreateReport, enabled: true },
                            ],
                        },
                        RibbonItem::LargeButton {
                            label: "Solution\nData".into(),
                            icon: '\u{1F4C8}',
                            action: RibbonAction::SolutionData,
                            enabled: true,
                            tooltip: "View solution data".into(),
                        },
                    ],
                }],
            },
            RibbonGroup {
                label: "Field Overlays".into(),
                sub_groups: vec![RibbonSubGroup {
                    items: vec![
                        RibbonItem::LargeSplitButton {
                            label: "Plot\nFields".into(),
                            icon: '\u{1F30A}',
                            action: RibbonAction::PlotFields,
                            enabled: true,
                            tooltip: "Plot field overlays".into(),
                            menu_items: vec![
                                MenuItem { label: "E-Field".into(), action: RibbonAction::PlotEField, enabled: true },
                                MenuItem { label: "H-Field".into(), action: RibbonAction::PlotHField, enabled: true },
                                MenuItem { label: "SAR".into(), action: RibbonAction::PlotSAR, enabled: true },
                            ],
                        },
                        RibbonItem::SmallButton {
                            label: "Animate".into(),
                            icon: '\u{23EF}',
                            action: RibbonAction::Animate,
                            enabled: true,
                            tooltip: "Animate field results".into(),
                        },
                    ],
                }],
            },
        ],
    }
}
