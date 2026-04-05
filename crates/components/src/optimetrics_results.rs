// ---------------------------------------------------------------------------
// Optimetrics Results Visualization — Charts and tables for sweep/optimization
// ---------------------------------------------------------------------------

use std::collections::HashMap;

use egui::Ui;
use egui_plot::{Bar, BarChart, Line, Plot, PlotPoints};

use emstudio_domain::optimetrics_result::OptimetricsSummary;
use emstudio_domain::sensitivity::SensitivityResult;
use emstudio_domain::statistical::MonteCarloResult;

// ---------------------------------------------------------------------------
// Parametric Sweep Summary Table
// ---------------------------------------------------------------------------

/// Render a parametric sweep summary as a table.
pub fn render_sweep_table(ui: &mut Ui, summary: &OptimetricsSummary) {
    ui.heading(&summary.optimetrics_name);
    ui.label(format!(
        "Completed: {}/{}  Failed: {}",
        summary.completed_variations, summary.total_variations, summary.failed_variations
    ));
    ui.separator();

    if summary.variations.is_empty() {
        ui.label("No variation data");
        return;
    }

    egui::ScrollArea::both().show(ui, |ui| {
        egui::Grid::new("sweep_table")
            .striped(true)
            .min_col_width(80.0)
            .show(ui, |ui| {
                // Header
                ui.strong("#");
                for var in &summary.swept_variables {
                    ui.strong(var);
                }
                for out in &summary.output_variables {
                    ui.strong(out);
                }
                ui.strong("Status");
                ui.end_row();

                // Data rows
                for v in &summary.variations {
                    ui.label(format!("{}", v.index));
                    for var_name in &summary.swept_variables {
                        let val = v.variables.get(var_name).map_or("-", |s| s.as_str());
                        ui.label(val);
                    }
                    for out_name in &summary.output_variables {
                        if let Some(ov) = v.outputs.get(out_name) {
                            ui.label(format!("{:.4} {}", ov.value, ov.unit));
                        } else {
                            ui.label("-");
                        }
                    }
                    ui.label(&v.status);
                    ui.end_row();
                }
            });
    });
}

// ---------------------------------------------------------------------------
// Optimization Convergence Curve
// ---------------------------------------------------------------------------

/// Render the optimization cost vs iteration curve.
pub fn render_optimization_convergence(ui: &mut Ui, summary: &OptimetricsSummary) {
    ui.heading(format!("{} — Convergence", summary.optimetrics_name));

    if let Some(alg) = &summary.algorithm {
        ui.label(format!(
            "Algorithm: {}  Converged: {}",
            alg,
            summary.converged.map_or("N/A", |c| if c { "Yes" } else { "No" })
        ));
    }

    let curve = summary.cost_curve();
    if curve.is_empty() {
        ui.label("No convergence data");
        return;
    }

    let plot = Plot::new("optim_convergence")
        .x_axis_label("Iteration")
        .y_axis_label("Cost")
        .legend(egui_plot::Legend::default())
        .allow_zoom(true)
        .show_grid(true);

    plot.show(ui, |plot_ui| {
        let line = Line::new("Cost", PlotPoints::new(curve))
            .color(egui::Color32::LIGHT_GREEN)
            .stroke(egui::Stroke::new(2.0, egui::Color32::LIGHT_GREEN));
        plot_ui.line(line);
    });

    // Best result
    if let Some(best) = &summary.best_result {
        ui.separator();
        ui.label(format!(
            "Best: iteration {} | cost = {:.6}",
            best.iteration, best.cost
        ));
        for (k, v) in &best.variables {
            ui.label(format!("  {} = {}", k, v));
        }
    }
}

// ---------------------------------------------------------------------------
// Parameter Comparison Overlay
// ---------------------------------------------------------------------------

/// Render overlaid traces comparing output across variations.
pub fn render_parameter_comparison(
    ui: &mut Ui,
    summary: &OptimetricsSummary,
    output_variable: &str,
) {
    ui.heading(format!("{} — {}", summary.optimetrics_name, output_variable));

    if summary.variations.is_empty() {
        ui.label("No data");
        return;
    }

    // Collect points: x = variation index, y = output value
    let points: Vec<[f64; 2]> = summary
        .variations
        .iter()
        .filter_map(|v| {
            v.outputs
                .get(output_variable)
                .map(|ov| [v.index as f64, ov.value])
        })
        .collect();

    let plot = Plot::new("param_comparison")
        .x_axis_label("Variation")
        .y_axis_label(output_variable)
        .show_grid(true);

    plot.show(ui, |plot_ui| {
        let line = Line::new(output_variable, PlotPoints::new(points))
            .color(egui::Color32::LIGHT_BLUE)
            .stroke(egui::Stroke::new(2.0, egui::Color32::LIGHT_BLUE));
        plot_ui.line(line);
    });
}

// ---------------------------------------------------------------------------
// Sensitivity Bar Chart
// ---------------------------------------------------------------------------

/// Render sensitivity gradients as a horizontal bar chart.
pub fn render_sensitivity_bars(ui: &mut Ui, result: &SensitivityResult) {
    ui.heading("Sensitivity Analysis");
    ui.label(format!("Base output: {:.6}", result.base_value));
    ui.separator();

    if result.variable_names.is_empty() {
        ui.label("No sensitivity data");
        return;
    }

    let bars: Vec<Bar> = result
        .variable_names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            Bar::new(i as f64, result.absolute_sensitivities[i])
                .name(name)
                .width(0.6)
        })
        .collect();

    let chart = BarChart::new("Sensitivity", bars).color(egui::Color32::from_rgb(70, 130, 180));

    let plot = Plot::new("sensitivity_chart")
        .y_axis_label("Absolute Sensitivity")
        .show_grid(true)
        .legend(egui_plot::Legend::default());

    plot.show(ui, |plot_ui| {
        plot_ui.bar_chart(chart);
    });

    // Table with details
    ui.separator();
    egui::Grid::new("sens_table")
        .striped(true)
        .show(ui, |ui| {
            ui.strong("Variable");
            ui.strong("Gradient");
            ui.strong("|Sensitivity|");
            ui.end_row();

            for (i, name) in result.variable_names.iter().enumerate() {
                ui.label(name);
                ui.label(format!("{:.6}", result.gradients[i]));
                ui.label(format!("{:.6}", result.absolute_sensitivities[i]));
                ui.end_row();
            }
        });
}

// ---------------------------------------------------------------------------
// Monte Carlo Histogram + Statistics
// ---------------------------------------------------------------------------

/// Render Monte Carlo results as histogram + statistics summary.
pub fn render_monte_carlo(ui: &mut Ui, result: &MonteCarloResult) {
    ui.heading("Statistical Analysis (Monte Carlo)");
    ui.label(format!(
        "Trials: {}  Mean: {:.4}  Std: {:.4}",
        result.outputs.len(),
        result.mean,
        result.std_dev
    ));
    ui.label(format!(
        "Min: {:.4}  Max: {:.4}",
        result.min, result.max
    ));

    // Percentiles
    ui.horizontal(|ui| {
        for (label, val) in &result.percentiles {
            ui.label(format!("{}: {:.4}", label, val));
        }
    });

    ui.separator();

    // Histogram
    let (bin_centers, counts) =
        emstudio_domain::statistical::build_histogram(&result.outputs, 25);

    if !bin_centers.is_empty() {
        let bars: Vec<Bar> = bin_centers
            .iter()
            .zip(counts.iter())
            .map(|(&center, &count)| Bar::new(center, count as f64).width(bin_centers[1] - bin_centers[0]))
            .collect();

        let chart = BarChart::new("Histogram", bars).color(egui::Color32::from_rgb(100, 149, 237));

        let plot = Plot::new("mc_histogram")
            .x_axis_label("Output Value")
            .y_axis_label("Count")
            .show_grid(true);

        plot.show(ui, |plot_ui| {
            plot_ui.bar_chart(chart);
        });
    }
}

// ---------------------------------------------------------------------------
// Combined Optimetrics Results View
// ---------------------------------------------------------------------------

/// A wrapper that renders appropriate visualization based on summary type.
pub fn render_optimetrics_results(
    ui: &mut Ui,
    summary: &OptimetricsSummary,
    sensitivity: Option<&SensitivityResult>,
    monte_carlo: Option<&MonteCarloResult>,
) {
    if summary.is_sweep() {
        render_sweep_table(ui, summary);
        if !summary.output_variables.is_empty() {
            ui.separator();
            render_parameter_comparison(ui, summary, &summary.output_variables[0]);
        }
    } else if summary.is_optimization() {
        render_optimization_convergence(ui, summary);
    }

    if let Some(sens) = sensitivity {
        ui.separator();
        render_sensitivity_bars(ui, sens);
    }

    if let Some(mc) = monte_carlo {
        ui.separator();
        render_monte_carlo(ui, mc);
    }
}
