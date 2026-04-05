// ---------------------------------------------------------------------------
// Tuning Panel — Interactive variable sliders with real-time feedback
// ---------------------------------------------------------------------------

use std::collections::HashMap;

use egui::Ui;

/// Interactive tuning panel with per-variable sliders.
pub struct TuningPanel {
    pub variable_names: Vec<String>,
    pub values: Vec<f64>,
    pub ranges: Vec<(f64, f64)>,
    /// Display results from the last evaluation.
    pub results: HashMap<String, f64>,
    /// Whether any slider changed since last evaluation.
    pub dirty: bool,
}

impl TuningPanel {
    pub fn new(
        variable_names: Vec<String>,
        initial_values: Vec<f64>,
        ranges: Vec<(f64, f64)>,
    ) -> Self {
        Self {
            variable_names,
            values: initial_values,
            ranges,
            results: HashMap::new(),
            dirty: true,
        }
    }

    /// Main UI rendering. Returns true if any value changed.
    pub fn ui(&mut self, ui: &mut Ui) -> bool {
        ui.heading("Interactive Tuning");
        ui.separator();

        let mut changed = false;

        for i in 0..self.variable_names.len() {
            let name = &self.variable_names[i];
            let (lo, hi) = self.ranges[i];
            ui.horizontal(|ui| {
                ui.label(format!("{}:", name));
                let prev = self.values[i];
                let val_text = format!("{:.4}", self.values[i]);
                ui.add(
                    egui::Slider::new(&mut self.values[i], lo..=hi)
                        .text(val_text)
                        .min_decimals(2)
                        .max_decimals(6),
                );
                if (self.values[i] - prev).abs() > 1e-15 {
                    changed = true;
                }
            });
        }

        if changed {
            self.dirty = true;
        }

        // Show results
        if !self.results.is_empty() {
            ui.separator();
            ui.heading("Results");
            egui::Grid::new("tuning_results")
                .striped(true)
                .show(ui, |ui| {
                    for (name, value) in &self.results {
                        ui.label(name);
                        ui.label(format!("{:.6}", value));
                        ui.end_row();
                    }
                });
        }

        // Reset button
        ui.separator();
        if ui.button("Reset to Center").clicked() {
            for i in 0..self.values.len() {
                self.values[i] = (self.ranges[i].0 + self.ranges[i].1) / 2.0;
            }
            self.dirty = true;
            changed = true;
        }

        changed
    }

    /// Get current variable values as a HashMap.
    pub fn current_values(&self) -> HashMap<String, f64> {
        self.variable_names
            .iter()
            .zip(self.values.iter())
            .map(|(n, &v)| (n.clone(), v))
            .collect()
    }

    /// Update displayed results.
    pub fn set_results(&mut self, results: HashMap<String, f64>) {
        self.results = results;
        self.dirty = false;
    }
}
