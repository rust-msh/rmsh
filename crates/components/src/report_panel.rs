// ---------------------------------------------------------------------------
// ReportPanel — 2D chart visualization component using egui_plot
// ---------------------------------------------------------------------------

use std::collections::HashMap;
use std::f64::consts::PI;

use egui::Ui;
use egui_plot::{HLine, Line, Plot, PlotPoints, Points, Text};

use emstudio_domain::report::{ChartType, Report, TraceStyle};

// ---------------------------------------------------------------------------
// Trace data
// ---------------------------------------------------------------------------

/// Evaluated trace data ready for rendering.
#[derive(Debug, Clone)]
pub struct TraceData {
    pub name: String,
    pub points: Vec<[f64; 2]>,
    pub color: egui::Color32,
    pub line_width: f32,
}

// ---------------------------------------------------------------------------
// ReportPanel
// ---------------------------------------------------------------------------

pub struct ReportPanel {
    pub report: Report,
    pub trace_data: HashMap<String, TraceData>,
    /// Interactive markers placed by user.
    pub active_markers: Vec<ActiveMarker>,
    /// Whether trace data needs re-evaluation.
    pub dirty: bool,
}

#[derive(Debug, Clone)]
pub struct ActiveMarker {
    pub name: String,
    pub trace_name: String,
    pub x: f64,
    pub y: f64,
}

impl ReportPanel {
    pub fn new(report: Report) -> Self {
        Self {
            report,
            trace_data: HashMap::new(),
            active_markers: Vec::new(),
            dirty: true,
        }
    }

    /// Set evaluated trace data. Called after expression evaluation.
    pub fn set_trace_data(&mut self, name: String, points: Vec<[f64; 2]>, style: Option<&TraceStyle>) {
        let (color, line_width) = style
            .map(|s| {
                (
                    egui::Color32::from_rgb(s.color[0], s.color[1], s.color[2]),
                    s.line_width as f32,
                )
            })
            .unwrap_or((egui::Color32::BLUE, 2.0));

        self.trace_data.insert(
            name.clone(),
            TraceData {
                name,
                points,
                color,
                line_width,
            },
        );
        self.dirty = false;
    }

    /// Main UI rendering.
    pub fn ui(&mut self, ui: &mut Ui) {
        match self.report.chart_type {
            ChartType::Rectangular => self.render_rectangular(ui),
            ChartType::Smith => self.render_smith_chart(ui),
            ChartType::Polar => self.render_polar_chart(ui),
            ChartType::DataTable => self.render_data_table(ui),
            ChartType::MatrixTable => self.render_data_table(ui),
            ChartType::Polar3D => {
                ui.label("3D Polar chart not yet implemented");
            }
        }
    }

    // -----------------------------------------------------------------------
    // Rectangular plot (egui_plot)
    // -----------------------------------------------------------------------

    fn render_rectangular(&self, ui: &mut Ui) {
        let x_label = self
            .report
            .x_axis
            .as_ref()
            .map(|a| format!("{} ({})", a.label, a.unit))
            .unwrap_or_else(|| "X".into());
        let y_label = self
            .report
            .y_axis
            .as_ref()
            .map(|a| format!("{} ({})", a.label, a.unit))
            .unwrap_or_else(|| "Y".into());

        let plot = Plot::new(format!("report_plot_{}", self.report.name))
            .x_axis_label(x_label)
            .y_axis_label(y_label)
            .legend(egui_plot::Legend::default())
            .allow_zoom(true)
            .allow_drag(true)
            .allow_scroll(true)
            .show_axes(true)
            .show_grid(true);

        plot.show(ui, |plot_ui| {
            // Render traces
            for trace in &self.report.traces {
                if let Some(td) = self.trace_data.get(&trace.name) {
                    let line = Line::new(&td.name, PlotPoints::new(td.points.clone()))
                        .color(td.color)
                        .stroke(egui::Stroke::new(td.line_width, td.color));
                    plot_ui.line(line);
                }
            }

            // Render limit lines
            for ll in &self.report.limit_lines {
                let color = ll
                    .style
                    .as_ref()
                    .map(|s| egui::Color32::from_rgb(s.color[0], s.color[1], s.color[2]))
                    .unwrap_or(egui::Color32::RED);
                let hline = HLine::new(&ll.name, ll.y_value)
                    .color(color)
                    .style(egui_plot::LineStyle::Dashed { length: 8.0 });
                plot_ui.hline(hline);
            }

            // Render markers (from report definition)
            for marker in &self.report.markers {
                if let Some(td) = self.trace_data.get(&marker.trace) {
                    if let Ok(x_val) = marker.x_value.parse::<f64>() {
                        if let Some(y_val) = interpolate_at(&td.points, x_val) {
                            let pt = Points::new(&marker.name, vec![[x_val, y_val]])
                                .radius(5.0)
                                .color(egui::Color32::YELLOW);
                            plot_ui.points(pt);

                            let label = Text::new(
                                format!("{}_label", marker.name),
                                egui_plot::PlotPoint::new(x_val, y_val),
                                format!("  {} ({:.4}, {:.4})", marker.name, x_val, y_val),
                            );
                            plot_ui.text(label);
                        }
                    }
                }
            }

            // Render active (user-placed) markers
            for m in &self.active_markers {
                let pt = Points::new(&m.name, vec![[m.x, m.y]])
                    .radius(5.0)
                    .color(egui::Color32::GOLD);
                plot_ui.points(pt);

                let label = Text::new(
                    format!("{}_label", m.name),
                    egui_plot::PlotPoint::new(m.x, m.y),
                    format!("  {} ({:.4}, {:.4})", m.name, m.x, m.y),
                );
                plot_ui.text(label);
            }
        });
    }

    // -----------------------------------------------------------------------
    // Smith chart (custom egui painting)
    // -----------------------------------------------------------------------

    fn render_smith_chart(&self, ui: &mut Ui) {
        let available_size = ui.available_size();
        let size = available_size.x.min(available_size.y) - 20.0;
        let (response, painter) =
            ui.allocate_painter(egui::Vec2::new(size, size), egui::Sense::hover());
        let center = response.rect.center();
        let radius = size * 0.45;

        let grid_color = egui::Color32::from_gray(80);
        let grid_stroke = egui::Stroke::new(0.5, grid_color);

        // Unit circle (|Γ| = 1)
        painter.circle_stroke(center, radius, egui::Stroke::new(1.0, egui::Color32::GRAY));

        // Constant resistance circles
        for &r in &[0.0_f32, 0.2, 0.5, 1.0, 2.0, 5.0] {
            let cx = center.x + radius * r / (1.0 + r);
            let cr = radius / (1.0 + r);
            painter.circle_stroke(egui::Pos2::new(cx, center.y), cr, grid_stroke);
        }

        // Constant reactance arcs (positive and negative)
        for &x in &[0.2_f32, 0.5, 1.0, 2.0, 5.0] {
            let cx = center.x + radius;
            let arc_radius = radius / x;
            painter.circle_stroke(
                egui::Pos2::new(cx, center.y - arc_radius),
                arc_radius,
                grid_stroke,
            );
            painter.circle_stroke(
                egui::Pos2::new(cx, center.y + arc_radius),
                arc_radius,
                grid_stroke,
            );
        }

        // Horizontal axis
        painter.line_segment(
            [
                egui::Pos2::new(center.x - radius, center.y),
                egui::Pos2::new(center.x + radius, center.y),
            ],
            grid_stroke,
        );

        // Plot S-parameter traces on Smith chart (Γ plane: re on x, im on y)
        for trace in &self.report.traces {
            if let Some(td) = self.trace_data.get(&trace.name) {
                let color = td.color;
                let points: Vec<egui::Pos2> = td
                    .points
                    .iter()
                    .map(|&[re, im]| {
                        egui::Pos2::new(
                            center.x + (re as f32) * radius,
                            center.y - (im as f32) * radius,
                        )
                    })
                    .collect();

                for pair in points.windows(2) {
                    painter.line_segment([pair[0], pair[1]], egui::Stroke::new(td.line_width, color));
                }
                for &pt in &points {
                    painter.circle_filled(pt, 2.0, color);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Polar chart (custom egui painting)
    // -----------------------------------------------------------------------

    fn render_polar_chart(&self, ui: &mut Ui) {
        let available_size = ui.available_size();
        let size = available_size.x.min(available_size.y) - 20.0;
        let (response, painter) =
            ui.allocate_painter(egui::Vec2::new(size, size), egui::Sense::hover());
        let center = response.rect.center();
        let max_radius = size * 0.42;

        let grid_color = egui::Color32::from_gray(80);
        let grid_stroke = egui::Stroke::new(0.5, grid_color);

        // Find max value for auto-scaling
        let max_val = self
            .trace_data
            .values()
            .flat_map(|td| td.points.iter().map(|p| p[1]))
            .fold(f64::NEG_INFINITY, f64::max);
        let scale = if max_val > 0.0 {
            (max_radius as f64) / max_val
        } else {
            1.0
        };

        // Concentric circles (radial grid)
        for i in 1..=5 {
            let r = max_radius * i as f32 / 5.0;
            painter.circle_stroke(center, r, grid_stroke);
            let val = max_val * i as f64 / 5.0;
            painter.text(
                egui::Pos2::new(center.x + 3.0, center.y - r - 2.0),
                egui::Align2::LEFT_BOTTOM,
                format!("{:.1}", val),
                egui::FontId::proportional(9.0),
                grid_color,
            );
        }

        // Radial lines (every 30 degrees)
        for deg in (0..360).step_by(30) {
            let angle = (deg as f64) * PI / 180.0;
            let cos_a = angle.cos() as f32;
            let sin_a = angle.sin() as f32;
            let end = egui::Pos2::new(
                center.x + max_radius * cos_a,
                center.y - max_radius * sin_a,
            );
            painter.line_segment([center, end], grid_stroke);
            painter.text(
                egui::Pos2::new(
                    center.x + (max_radius + 12.0) * cos_a,
                    center.y - (max_radius + 12.0) * sin_a,
                ),
                egui::Align2::CENTER_CENTER,
                format!("{}deg", deg),
                egui::FontId::proportional(10.0),
                egui::Color32::from_gray(160),
            );
        }

        // Plot traces: points are [angle_deg, value]
        for trace in &self.report.traces {
            if let Some(td) = self.trace_data.get(&trace.name) {
                let color = td.color;
                let screen_points: Vec<egui::Pos2> = td
                    .points
                    .iter()
                    .map(|&[angle_deg, value]| {
                        let angle_rad = angle_deg * PI / 180.0;
                        let r = (value * scale) as f32;
                        egui::Pos2::new(
                            center.x + r * (angle_rad.cos() as f32),
                            center.y - r * (angle_rad.sin() as f32),
                        )
                    })
                    .collect();

                for pair in screen_points.windows(2) {
                    painter.line_segment(
                        [pair[0], pair[1]],
                        egui::Stroke::new(td.line_width, color),
                    );
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Data table
    // -----------------------------------------------------------------------

    fn render_data_table(&self, ui: &mut Ui) {
        if self.trace_data.is_empty() {
            ui.label("No data loaded");
            return;
        }

        egui::ScrollArea::both().show(ui, |ui| {
            egui::Grid::new(format!("report_table_{}", self.report.name))
                .striped(true)
                .min_col_width(80.0)
                .show(ui, |ui| {
                    // Header row
                    ui.strong("X");
                    for trace in &self.report.traces {
                        ui.strong(&trace.name);
                    }
                    ui.end_row();

                    let max_len = self
                        .trace_data
                        .values()
                        .map(|td| td.points.len())
                        .max()
                        .unwrap_or(0);

                    for i in 0..max_len {
                        let x_val = self
                            .trace_data
                            .values()
                            .next()
                            .and_then(|td| td.points.get(i))
                            .map(|p| p[0]);
                        if let Some(x) = x_val {
                            ui.label(format!("{:.6}", x));
                        } else {
                            ui.label("-");
                        }

                        for trace in &self.report.traces {
                            if let Some(td) = self.trace_data.get(&trace.name) {
                                if let Some(pt) = td.points.get(i) {
                                    ui.label(format!("{:.6}", pt[1]));
                                } else {
                                    ui.label("-");
                                }
                            } else {
                                ui.label("-");
                            }
                        }
                        ui.end_row();
                    }
                });
        });
    }

    // -----------------------------------------------------------------------
    // CSV export
    // -----------------------------------------------------------------------

    /// Export trace data to CSV format string.
    pub fn export_csv(&self) -> String {
        let mut csv = String::new();

        csv.push_str("X");
        for trace in &self.report.traces {
            csv.push(',');
            csv.push_str(&trace.name);
        }
        csv.push('\n');

        let max_len = self
            .trace_data
            .values()
            .map(|td| td.points.len())
            .max()
            .unwrap_or(0);

        for i in 0..max_len {
            let x_val = self
                .trace_data
                .values()
                .next()
                .and_then(|td| td.points.get(i))
                .map(|p| p[0]);
            if let Some(x) = x_val {
                csv.push_str(&format!("{:.10}", x));
            }

            for trace in &self.report.traces {
                csv.push(',');
                if let Some(td) = self.trace_data.get(&trace.name) {
                    if let Some(pt) = td.points.get(i) {
                        csv.push_str(&format!("{:.10}", pt[1]));
                    }
                }
            }
            csv.push('\n');
        }

        csv
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Linear interpolation to find y at a given x in sorted data.
fn interpolate_at(points: &[[f64; 2]], x: f64) -> Option<f64> {
    if points.is_empty() {
        return None;
    }
    if points.len() == 1 {
        return Some(points[0][1]);
    }

    for i in 0..points.len() - 1 {
        let (x0, y0) = (points[i][0], points[i][1]);
        let (x1, y1) = (points[i + 1][0], points[i + 1][1]);
        if x >= x0 && x <= x1 {
            let t = (x - x0) / (x1 - x0);
            return Some(y0 + t * (y1 - y0));
        }
    }

    if x <= points[0][0] {
        Some(points[0][1])
    } else {
        Some(points.last().unwrap()[1])
    }
}
