// ---------------------------------------------------------------------------
// Q3D Net/Terminal Panel — UI for editing Q3D electrical networks
// ---------------------------------------------------------------------------

use egui::Ui;

use emstudio_domain::boundary::{Assignment, AssignmentTarget};
use emstudio_domain::excitation::ExcitationType;
use emstudio_domain::net::{Net, Terminal};

/// Response from the Q3D net panel.
pub struct Q3dNetResponse {
    /// Whether any nets were modified.
    pub nets_modified: bool,
    /// Currently selected net name.
    pub selected_net: Option<String>,
}

/// State for the Q3D net panel.
pub struct Q3dNetPanelState {
    pub selected_net_idx: Option<usize>,
    pub new_net_name: String,
    pub new_terminal_name: String,
}

impl Default for Q3dNetPanelState {
    fn default() -> Self {
        Self {
            selected_net_idx: None,
            new_net_name: String::new(),
            new_terminal_name: String::new(),
        }
    }
}

/// Render the Q3D Net/Terminal editing panel.
pub fn show_q3d_net_panel(
    ui: &mut Ui,
    nets: &mut Vec<Net>,
    geo_objects: &[String],
    state: &mut Q3dNetPanelState,
) -> Q3dNetResponse {
    let mut modified = false;

    ui.heading("Q3D Networks");
    ui.separator();

    // Net list with add/remove buttons
    ui.horizontal(|ui| {
        ui.label("Nets:");
        if ui.button("+").clicked() {
            let name = if state.new_net_name.is_empty() {
                format!("Net{}", nets.len() + 1)
            } else {
                state.new_net_name.clone()
            };
            nets.push(Net {
                name,
                objects: Vec::new(),
                is_ground_reference: false,
                terminals: Vec::new(),
            });
            state.selected_net_idx = Some(nets.len() - 1);
            state.new_net_name.clear();
            modified = true;
        }
        if let Some(idx) = state.selected_net_idx {
            if idx < nets.len() && ui.button("-").clicked() {
                nets.remove(idx);
                state.selected_net_idx = if nets.is_empty() {
                    None
                } else {
                    Some(idx.min(nets.len() - 1))
                };
                modified = true;
            }
        }
    });

    // Net selection list
    for (i, net) in nets.iter().enumerate() {
        let label = if net.is_ground_reference {
            format!("{} (GND)", net.name)
        } else {
            net.name.clone()
        };
        let selected = state.selected_net_idx == Some(i);
        if ui.selectable_label(selected, &label).clicked() {
            state.selected_net_idx = Some(i);
        }
    }

    ui.separator();

    // Selected net details
    let selected_net_name = if let Some(idx) = state.selected_net_idx {
        if idx < nets.len() {
            let net = &mut nets[idx];

            ui.label("Net Name:");
            let prev_name = net.name.clone();
            ui.text_edit_singleline(&mut net.name);
            if net.name != prev_name {
                modified = true;
            }

            // Ground reference toggle
            let prev_gnd = net.is_ground_reference;
            ui.checkbox(&mut net.is_ground_reference, "Ground Reference");
            if net.is_ground_reference != prev_gnd {
                modified = true;
            }

            // Object assignment
            ui.separator();
            ui.label("Assigned Objects:");
            let mut objects_changed = false;
            let mut removal_idx = None;
            for (oi, obj_name) in net.objects.iter().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(obj_name);
                    if ui.small_button("x").clicked() {
                        removal_idx = Some(oi);
                    }
                });
            }
            if let Some(ri) = removal_idx {
                net.objects.remove(ri);
                objects_changed = true;
            }

            // Add object from available geometry objects
            let available: Vec<&String> = geo_objects
                .iter()
                .filter(|o| !net.objects.contains(o))
                .collect();
            if !available.is_empty() {
                ui.horizontal(|ui| {
                    ui.label("Add:");
                    for obj in available {
                        if ui.small_button(obj).clicked() {
                            net.objects.push(obj.clone());
                            objects_changed = true;
                        }
                    }
                });
            }
            if objects_changed {
                modified = true;
            }

            // Terminal editing
            ui.separator();
            ui.label("Terminals:");

            let mut terminal_removal = None;
            for (ti, terminal) in net.terminals.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.text_edit_singleline(&mut terminal.name);
                    // Terminal type is read-only display for now
                    ui.label(format!("{:?}", terminal.terminal_type));
                    if ui.small_button("x").clicked() {
                        terminal_removal = Some(ti);
                    }
                });
            }
            if let Some(ri) = terminal_removal {
                net.terminals.remove(ri);
                modified = true;
            }

            // Add terminal
            ui.horizontal(|ui| {
                ui.text_edit_singleline(&mut state.new_terminal_name);
                if ui.button("Add Terminal").clicked() && !state.new_terminal_name.is_empty() {
                    net.terminals.push(Terminal {
                        name: state.new_terminal_name.clone(),
                        terminal_type: ExcitationType::Source,
                        assignment: Assignment {
                            target_type: AssignmentTarget::Face,
                            targets: Vec::new(),
                        },
                    });
                    state.new_terminal_name.clear();
                    modified = true;
                }
            });

            Some(net.name.clone())
        } else {
            None
        }
    } else {
        ui.label("Select a net to edit.");
        None
    };

    Q3dNetResponse {
        nets_modified: modified,
        selected_net: selected_net_name,
    }
}
