// ---------------------------------------------------------------------------
// Optimetrics Setup Panel — UI for configuring parametric sweeps, optimization,
// sensitivity, statistical analysis, and tuning setups
// ---------------------------------------------------------------------------

use egui::Ui;

use emstudio_domain::optimetrics::OptimetricsSetup;

/// UI panel for managing optimetrics setups.
pub struct OptimetricsPanel {
    pub setups: Vec<OptimetricsSetup>,
    pub selected_idx: Option<usize>,
    /// Type template for the "Add" dropdown.
    add_type: OptimetricsType,
    /// Name for new setup.
    new_name: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OptimetricsType {
    ParametricSweep,
    Optimization,
    Sensitivity,
    Statistical,
    Tuning,
}

impl OptimetricsType {
    const ALL: &[OptimetricsType] = &[
        OptimetricsType::ParametricSweep,
        OptimetricsType::Optimization,
        OptimetricsType::Sensitivity,
        OptimetricsType::Statistical,
        OptimetricsType::Tuning,
    ];

    fn name(&self) -> &'static str {
        match self {
            Self::ParametricSweep => "Parametric Sweep",
            Self::Optimization => "Optimization",
            Self::Sensitivity => "Sensitivity",
            Self::Statistical => "Statistical",
            Self::Tuning => "Tuning",
        }
    }
}

impl OptimetricsPanel {
    pub fn new(setups: Vec<OptimetricsSetup>) -> Self {
        Self {
            setups,
            selected_idx: None,
            add_type: OptimetricsType::ParametricSweep,
            new_name: "NewSetup".into(),
        }
    }

    /// Main UI rendering.
    pub fn ui(&mut self, ui: &mut Ui) {
        ui.heading("Optimetrics");
        ui.separator();

        // Add controls
        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt("optim_type")
                .selected_text(self.add_type.name())
                .show_ui(ui, |ui| {
                    for &t in OptimetricsType::ALL {
                        ui.selectable_value(&mut self.add_type, t, t.name());
                    }
                });
            ui.text_edit_singleline(&mut self.new_name);
            if ui.button("Add").clicked() && !self.new_name.is_empty() {
                let setup = create_default_setup(&self.new_name, self.add_type);
                self.setups.push(setup);
                self.selected_idx = Some(self.setups.len() - 1);
                self.new_name = format!("Setup{}", self.setups.len() + 1);
            }
        });

        if ui.button("Delete Selected").clicked() {
            if let Some(idx) = self.selected_idx {
                if idx < self.setups.len() {
                    self.setups.remove(idx);
                    self.selected_idx = if self.setups.is_empty() {
                        None
                    } else {
                        Some(idx.min(self.setups.len() - 1))
                    };
                }
            }
        }

        ui.separator();

        // Setup list
        egui::ScrollArea::vertical()
            .max_height(150.0)
            .show(ui, |ui| {
                for (i, setup) in self.setups.iter().enumerate() {
                    let (name, type_name, enabled) = setup_info(setup);
                    let label = format!("[{}] {} ({})", if enabled { "x" } else { " " }, name, type_name);
                    let selected = self.selected_idx == Some(i);
                    if ui.selectable_label(selected, &label).clicked() {
                        self.selected_idx = Some(i);
                    }
                }
            });

        ui.separator();

        // Detail editor for selected setup
        if let Some(idx) = self.selected_idx {
            if idx < self.setups.len() {
                let setup = &mut self.setups[idx];
                render_setup_editor(ui, setup);
            }
        }
    }
}

fn setup_info(setup: &OptimetricsSetup) -> (&str, &str, bool) {
    match setup {
        OptimetricsSetup::ParametricSweep { name, enabled, .. } => (name, "Sweep", *enabled),
        OptimetricsSetup::Optimization { name, enabled, .. } => (name, "Optimize", *enabled),
        OptimetricsSetup::Sensitivity { name, enabled, .. } => (name, "Sensitivity", *enabled),
        OptimetricsSetup::Statistical { name, enabled, .. } => (name, "Statistical", *enabled),
        OptimetricsSetup::Tuning { name, enabled, .. } => (name, "Tuning", *enabled),
    }
}

fn render_setup_editor(ui: &mut Ui, setup: &mut OptimetricsSetup) {
    match setup {
        OptimetricsSetup::ParametricSweep {
            name,
            enabled,
            setup: analysis_setup,
            sweep_definitions,
            ..
        } => {
            ui.heading("Parametric Sweep");
            ui.horizontal(|ui| {
                ui.label("Name:");
                ui.text_edit_singleline(name);
            });
            ui.checkbox(enabled, "Enabled");
            ui.horizontal(|ui| {
                ui.label("Analysis Setup:");
                ui.text_edit_singleline(analysis_setup);
            });

            ui.separator();
            ui.label(format!("Sweep Definitions: {}", sweep_definitions.len()));

            for (i, def) in sweep_definitions.iter_mut().enumerate() {
                ui.group(|ui| {
                    ui.label(format!("Variable #{}", i + 1));
                    ui.horizontal(|ui| {
                        ui.label("Variable:");
                        ui.text_edit_singleline(&mut def.variable);
                    });
                    ui.horizontal(|ui| {
                        ui.label("Type:");
                        ui.text_edit_singleline(&mut def.sweep_type);
                    });
                    ui.horizontal(|ui| {
                        ui.label("Params:");
                        let mut params_str = def.params.to_string();
                        if ui.text_edit_singleline(&mut params_str).changed() {
                            if let Ok(v) = serde_json::from_str(&params_str) {
                                def.params = v;
                            }
                        }
                    });
                });
            }
        }

        OptimetricsSetup::Optimization {
            name,
            enabled,
            setup: analysis_setup,
            algorithm,
            max_iterations,
            variables,
            goals,
            ..
        } => {
            ui.heading("Optimization");
            ui.horizontal(|ui| {
                ui.label("Name:");
                ui.text_edit_singleline(name);
            });
            ui.checkbox(enabled, "Enabled");
            ui.horizontal(|ui| {
                ui.label("Setup:");
                ui.text_edit_singleline(analysis_setup);
            });
            ui.horizontal(|ui| {
                ui.label("Algorithm:");
                egui::ComboBox::from_id_salt("opt_algo")
                    .selected_text(algorithm.as_str())
                    .show_ui(ui, |ui| {
                        for &alg in &["QuasiNewton", "PatternSearch", "GeneticAlgorithm", "SNLP"]
                        {
                            ui.selectable_value(algorithm, alg.to_string(), alg);
                        }
                    });
            });
            ui.add(egui::Slider::new(max_iterations, 10..=1000).text("Max Iterations"));

            ui.separator();
            ui.label(format!("Variables: {}", variables.len()));
            for var in variables.iter_mut() {
                ui.horizontal(|ui| {
                    ui.text_edit_singleline(&mut var.variable);
                    ui.label("min:");
                    ui.text_edit_singleline(&mut var.min);
                    ui.label("max:");
                    ui.text_edit_singleline(&mut var.max);
                });
            }

            ui.separator();
            ui.label(format!("Goals: {}", goals.len()));
            for (i, goal) in goals.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(format!("Goal #{}", i + 1));
                    ui.horizontal(|ui| {
                        ui.label("Name:");
                        ui.text_edit_singleline(&mut goal.name);
                    });
                    ui.horizontal(|ui| {
                        ui.label("Expr:");
                        ui.text_edit_singleline(&mut goal.expression);
                    });
                    ui.horizontal(|ui| {
                        ui.label("Condition:");
                        ui.text_edit_singleline(&mut goal.condition);
                    });
                });
            }
        }

        OptimetricsSetup::Sensitivity {
            name,
            enabled,
            setup: analysis_setup,
            variables,
            output,
            num_samples,
        } => {
            ui.heading("Sensitivity Analysis");
            ui.horizontal(|ui| {
                ui.label("Name:");
                ui.text_edit_singleline(name);
            });
            ui.checkbox(enabled, "Enabled");
            ui.horizontal(|ui| {
                ui.label("Setup:");
                ui.text_edit_singleline(analysis_setup);
            });
            ui.horizontal(|ui| {
                ui.label("Output:");
                ui.text_edit_singleline(output);
            });
            ui.add(egui::Slider::new(num_samples, 10..=500).text("Samples"));

            ui.separator();
            ui.label(format!("Variables: {}", variables.len()));
            for var in variables.iter_mut() {
                ui.horizontal(|ui| {
                    ui.text_edit_singleline(&mut var.variable);
                    ui.label("variation:");
                    ui.text_edit_singleline(&mut var.variation);
                    ui.label("dist:");
                    ui.text_edit_singleline(&mut var.distribution);
                });
            }
        }

        OptimetricsSetup::Statistical {
            name,
            enabled,
            setup: analysis_setup,
            variables,
            num_trials,
        } => {
            ui.heading("Statistical Analysis (Monte Carlo)");
            ui.horizontal(|ui| {
                ui.label("Name:");
                ui.text_edit_singleline(name);
            });
            ui.checkbox(enabled, "Enabled");
            ui.horizontal(|ui| {
                ui.label("Setup:");
                ui.text_edit_singleline(analysis_setup);
            });
            ui.add(egui::Slider::new(num_trials, 100..=10000).text("Trials"));

            ui.separator();
            for var in variables.iter_mut() {
                ui.horizontal(|ui| {
                    ui.text_edit_singleline(&mut var.variable);
                    ui.label("±");
                    ui.text_edit_singleline(&mut var.variation);
                    ui.text_edit_singleline(&mut var.distribution);
                });
            }
        }

        OptimetricsSetup::Tuning {
            name,
            enabled,
            setup: analysis_setup,
            variables,
        } => {
            ui.heading("Interactive Tuning");
            ui.horizontal(|ui| {
                ui.label("Name:");
                ui.text_edit_singleline(name);
            });
            ui.checkbox(enabled, "Enabled");
            ui.horizontal(|ui| {
                ui.label("Setup:");
                ui.text_edit_singleline(analysis_setup);
            });

            ui.separator();
            ui.label("Tuning Variables:");
            for var in variables.iter_mut() {
                ui.text_edit_singleline(var);
            }
        }
    }
}

fn create_default_setup(name: &str, typ: OptimetricsType) -> OptimetricsSetup {
    match typ {
        OptimetricsType::ParametricSweep => OptimetricsSetup::ParametricSweep {
            name: name.into(),
            enabled: true,
            setup: "Setup1".into(),
            sweep_definitions: Vec::new(),
            constraints: Vec::new(),
            goals: Vec::new(),
        },
        OptimetricsType::Optimization => OptimetricsSetup::Optimization {
            name: name.into(),
            enabled: true,
            setup: "Setup1".into(),
            algorithm: "PatternSearch".into(),
            max_iterations: 100,
            variables: Vec::new(),
            goals: Vec::new(),
            constraints: Vec::new(),
        },
        OptimetricsType::Sensitivity => OptimetricsSetup::Sensitivity {
            name: name.into(),
            enabled: true,
            setup: "Setup1".into(),
            variables: Vec::new(),
            output: String::new(),
            num_samples: 50,
        },
        OptimetricsType::Statistical => OptimetricsSetup::Statistical {
            name: name.into(),
            enabled: true,
            setup: "Setup1".into(),
            variables: Vec::new(),
            num_trials: 1000,
        },
        OptimetricsType::Tuning => OptimetricsSetup::Tuning {
            name: name.into(),
            enabled: true,
            setup: "Setup1".into(),
            variables: Vec::new(),
        },
    }
}
