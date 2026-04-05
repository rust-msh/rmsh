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

/// Matrix data for heatmap display (RLCG MatrixTable).
#[derive(Debug, Clone)]
pub struct MatrixData {
    pub title: String,
    pub row_labels: Vec<String>,
    pub col_labels: Vec<String>,
    pub values: Vec<Vec<f64>>,
    pub unit: String,
}

/// 3D surface data for parameter sweep plots.
#[derive(Debug, Clone)]
pub struct SurfaceData {
    pub name: String,
    /// X-axis values (e.g. frequency).
    pub x_values: Vec<f64>,
    /// Y-axis values (e.g. parameter sweep variable).
    pub y_values: Vec<f64>,
    /// Z values as a 2D grid: z_grid[y_idx][x_idx].
    pub z_grid: Vec<Vec<f64>>,
}

// ---------------------------------------------------------------------------
// ReportPanel
// ---------------------------------------------------------------------------

pub struct ReportPanel {
    pub report: Report,
    pub trace_data: HashMap<String, TraceData>,
    /// Matrix data for heatmap display.
    pub matrix_data: Option<MatrixData>,
    /// Surface data for 3D rectangular plots.
    pub surface_data: Option<SurfaceData>,
    /// Interactive markers placed by user.
    pub active_markers: Vec<ActiveMarker>,
    /// Whether trace data needs re-evaluation.
    pub dirty: bool,
    /// Selected frequency index for matrix heatmap display.
    pub selected_freq_idx: usize,
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
            matrix_data: None,
            surface_data: None,
            active_markers: Vec::new(),
            dirty: true,
            selected_freq_idx: 0,
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
            ChartType::MatrixTable => {
                if self.matrix_data.is_some()
                    && self
                        .report
                        .display_options
                        .as_ref()
                        .and_then(|d| d.heatmap_enabled)
                        .unwrap_or(false)
                {
                    self.render_rlcg_heatmap(ui);
                } else {
                    self.render_data_table(ui);
                }
            }
            ChartType::Rectangular3D => self.render_rectangular_3d(ui),
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
    // RLCG Matrix Heatmap (egui painter)
    // -----------------------------------------------------------------------

    fn render_rlcg_heatmap(&mut self, ui: &mut Ui) {
        let matrix = match &self.matrix_data {
            Some(m) => m,
            None => {
                ui.label("No matrix data loaded");
                return;
            }
        };

        // Title
        ui.heading(&matrix.title);

        // Frequency selector (if multiple frequencies available)
        ui.horizontal(|ui| {
            ui.label("Freq index:");
            ui.add(
                egui::DragValue::new(&mut self.selected_freq_idx)
                    .range(0..=0usize) // Will be updated by caller with actual range
                    .speed(1.0),
            );
            if !matrix.unit.is_empty() {
                ui.label(format!("Unit: {}", matrix.unit));
            }
        });

        ui.separator();

        if matrix.values.is_empty() || matrix.row_labels.is_empty() {
            ui.label("Empty matrix");
            return;
        }

        let nrows = matrix.values.len();
        let ncols = matrix.values.first().map_or(0, |r| r.len());

        // Find value range
        let (mut vmin, mut vmax) = (f64::MAX, f64::MIN);
        for row in &matrix.values {
            for &v in row {
                vmin = vmin.min(v);
                vmax = vmax.max(v);
            }
        }
        if (vmax - vmin).abs() < 1e-20 {
            vmax = vmin + 1.0;
        }

        // Layout parameters
        let label_width = 80.0f32;
        let cell_size = 60.0f32;
        let header_height = 30.0f32;
        let colorbar_width = 60.0f32;

        let total_w = label_width + ncols as f32 * cell_size + colorbar_width + 20.0;
        let total_h = header_height + nrows as f32 * cell_size + 10.0;

        let (response, painter) =
            ui.allocate_painter(egui::Vec2::new(total_w, total_h), egui::Sense::hover());
        let origin = response.rect.min;

        let font = egui::FontId::proportional(10.0);
        let label_font = egui::FontId::proportional(9.0);

        // Draw column headers
        for (j, label) in matrix.col_labels.iter().enumerate().take(ncols) {
            let cx = origin.x + label_width + j as f32 * cell_size + cell_size * 0.5;
            let cy = origin.y + header_height * 0.5;
            painter.text(
                egui::Pos2::new(cx, cy),
                egui::Align2::CENTER_CENTER,
                label,
                label_font.clone(),
                egui::Color32::WHITE,
            );
        }

        // Draw rows
        for (i, row) in matrix.values.iter().enumerate().take(nrows) {
            let y_top = origin.y + header_height + i as f32 * cell_size;

            // Row label
            let label = matrix.row_labels.get(i).map_or("", |s| s.as_str());
            painter.text(
                egui::Pos2::new(origin.x + label_width - 4.0, y_top + cell_size * 0.5),
                egui::Align2::RIGHT_CENTER,
                label,
                label_font.clone(),
                egui::Color32::WHITE,
            );

            // Cells
            for (j, &val) in row.iter().enumerate().take(ncols) {
                let x_left = origin.x + label_width + j as f32 * cell_size;
                let cell_rect = egui::Rect::from_min_size(
                    egui::Pos2::new(x_left + 1.0, y_top + 1.0),
                    egui::Vec2::new(cell_size - 2.0, cell_size - 2.0),
                );

                let color = value_to_heatmap_color(val, vmin, vmax);
                painter.rect_filled(cell_rect, 2.0, color);

                // Value text (use black or white depending on brightness)
                let brightness = color.r() as u32 + color.g() as u32 + color.b() as u32;
                let text_color = if brightness > 384 {
                    egui::Color32::BLACK
                } else {
                    egui::Color32::WHITE
                };
                let decimal_places = self
                    .report
                    .display_options
                    .as_ref()
                    .and_then(|d| d.decimal_places)
                    .unwrap_or(4) as usize;
                painter.text(
                    cell_rect.center(),
                    egui::Align2::CENTER_CENTER,
                    format!("{:.prec$}", val, prec = decimal_places),
                    font.clone(),
                    text_color,
                );
            }
        }

        // Colorbar
        let bar_x = origin.x + label_width + ncols as f32 * cell_size + 10.0;
        let bar_top = origin.y + header_height;
        let bar_height = nrows as f32 * cell_size;
        let bar_width = 15.0f32;

        let steps = 32;
        for s in 0..steps {
            let frac_top = s as f32 / steps as f32;
            let frac_bot = (s + 1) as f32 / steps as f32;
            let t = 1.0 - (frac_top + frac_bot) * 0.5; // top=max, bottom=min
            let v = vmin + t as f64 * (vmax - vmin);
            let color = value_to_heatmap_color(v, vmin, vmax);

            let y0 = bar_top + frac_top * bar_height;
            let y1 = bar_top + frac_bot * bar_height;
            painter.rect_filled(
                egui::Rect::from_min_max(
                    egui::Pos2::new(bar_x, y0),
                    egui::Pos2::new(bar_x + bar_width, y1),
                ),
                0.0,
                color,
            );
        }

        painter.rect_stroke(
            egui::Rect::from_min_size(
                egui::Pos2::new(bar_x, bar_top),
                egui::Vec2::new(bar_width, bar_height),
            ),
            0.0,
            egui::Stroke::new(1.0, egui::Color32::GRAY),
            egui::StrokeKind::Outside,
        );

        // Colorbar labels
        let text_x = bar_x + bar_width + 4.0;
        painter.text(
            egui::Pos2::new(text_x, bar_top),
            egui::Align2::LEFT_TOP,
            format!("{:.3}", vmax),
            label_font.clone(),
            egui::Color32::WHITE,
        );
        painter.text(
            egui::Pos2::new(text_x, bar_top + bar_height * 0.5),
            egui::Align2::LEFT_CENTER,
            format!("{:.3}", (vmin + vmax) / 2.0),
            label_font.clone(),
            egui::Color32::WHITE,
        );
        painter.text(
            egui::Pos2::new(text_x, bar_top + bar_height),
            egui::Align2::LEFT_BOTTOM,
            format!("{:.3}", vmin),
            label_font,
            egui::Color32::WHITE,
        );
    }

    // -----------------------------------------------------------------------
    // 3D Rectangular Plot — 2D heatmap projection
    // -----------------------------------------------------------------------

    fn render_rectangular_3d(&self, ui: &mut Ui) {
        let surface = match &self.surface_data {
            Some(s) => s,
            None => {
                ui.label("No surface data loaded");
                return;
            }
        };

        if surface.x_values.is_empty() || surface.y_values.is_empty() || surface.z_grid.is_empty()
        {
            ui.label("Empty surface data");
            return;
        }

        ui.heading(&surface.name);
        ui.separator();

        let nx = surface.x_values.len();
        let ny = surface.y_values.len();

        // Find Z range
        let (mut zmin, mut zmax) = (f64::MAX, f64::MIN);
        for row in &surface.z_grid {
            for &z in row {
                zmin = zmin.min(z);
                zmax = zmax.max(z);
            }
        }
        if (zmax - zmin).abs() < 1e-20 {
            zmax = zmin + 1.0;
        }

        // Layout
        let margin = 50.0f32;
        let colorbar_space = 70.0f32;
        let available = ui.available_size();
        let plot_w = (available.x - margin * 2.0 - colorbar_space).max(100.0);
        let plot_h = (available.y - margin * 2.0).max(100.0);

        let (response, painter) = ui.allocate_painter(
            egui::Vec2::new(plot_w + margin * 2.0 + colorbar_space, plot_h + margin * 2.0),
            egui::Sense::hover(),
        );
        let origin = response.rect.min;
        let plot_origin = egui::Pos2::new(origin.x + margin, origin.y + margin);

        let cell_w = plot_w / nx as f32;
        let cell_h = plot_h / ny as f32;

        let label_font = egui::FontId::proportional(9.0);
        let grid_color = egui::Color32::from_gray(60);

        // Draw cells
        for iy in 0..ny {
            for ix in 0..nx {
                let z = surface
                    .z_grid
                    .get(iy)
                    .and_then(|row| row.get(ix))
                    .copied()
                    .unwrap_or(0.0);

                let x_left = plot_origin.x + ix as f32 * cell_w;
                let y_top = plot_origin.y + (ny - 1 - iy) as f32 * cell_h; // Y axis: bottom = low

                let cell_rect = egui::Rect::from_min_size(
                    egui::Pos2::new(x_left, y_top),
                    egui::Vec2::new(cell_w, cell_h),
                );

                let color = value_to_heatmap_color(z, zmin, zmax);
                painter.rect_filled(cell_rect, 0.0, color);
            }
        }

        // Plot border
        let plot_rect = egui::Rect::from_min_size(plot_origin, egui::Vec2::new(plot_w, plot_h));
        painter.rect_stroke(
            plot_rect,
            0.0,
            egui::Stroke::new(1.0, grid_color),
            egui::StrokeKind::Outside,
        );

        // X-axis labels (show a few ticks)
        let x_ticks = 5.min(nx);
        for i in 0..=x_ticks {
            let idx = if x_ticks > 0 {
                i * (nx - 1) / x_ticks
            } else {
                0
            };
            let x = plot_origin.x + idx as f32 * cell_w + cell_w * 0.5;
            let val = surface.x_values.get(idx).copied().unwrap_or(0.0);
            painter.text(
                egui::Pos2::new(x, plot_origin.y + plot_h + 4.0),
                egui::Align2::CENTER_TOP,
                format!("{:.2}", val),
                label_font.clone(),
                egui::Color32::from_gray(180),
            );
        }

        // Y-axis labels
        let y_ticks = 5.min(ny);
        for i in 0..=y_ticks {
            let idx = if y_ticks > 0 {
                i * (ny - 1) / y_ticks
            } else {
                0
            };
            let y = plot_origin.y + (ny - 1 - idx) as f32 * cell_h + cell_h * 0.5;
            let val = surface.y_values.get(idx).copied().unwrap_or(0.0);
            painter.text(
                egui::Pos2::new(plot_origin.x - 4.0, y),
                egui::Align2::RIGHT_CENTER,
                format!("{:.2}", val),
                label_font.clone(),
                egui::Color32::from_gray(180),
            );
        }

        // X-axis label
        if let Some(x_axis) = &self.report.x_axis {
            painter.text(
                egui::Pos2::new(
                    plot_origin.x + plot_w * 0.5,
                    plot_origin.y + plot_h + 20.0,
                ),
                egui::Align2::CENTER_TOP,
                format!("{} ({})", x_axis.label, x_axis.unit),
                label_font.clone(),
                egui::Color32::WHITE,
            );
        }

        // Y-axis label
        if let Some(y_axis) = &self.report.y_axis {
            painter.text(
                egui::Pos2::new(origin.x + 10.0, plot_origin.y + plot_h * 0.5),
                egui::Align2::CENTER_CENTER,
                format!("{} ({})", y_axis.label, y_axis.unit),
                label_font.clone(),
                egui::Color32::WHITE,
            );
        }

        // Colorbar
        let bar_x = plot_origin.x + plot_w + 10.0;
        let bar_w = 15.0f32;
        let bar_top = plot_origin.y;
        let bar_h = plot_h;
        let steps = 32;

        for s in 0..steps {
            let frac_top = s as f32 / steps as f32;
            let frac_bot = (s + 1) as f32 / steps as f32;
            let t = 1.0 - (frac_top + frac_bot) * 0.5;
            let v = zmin + t as f64 * (zmax - zmin);
            let color = value_to_heatmap_color(v, zmin, zmax);

            painter.rect_filled(
                egui::Rect::from_min_max(
                    egui::Pos2::new(bar_x, bar_top + frac_top * bar_h),
                    egui::Pos2::new(bar_x + bar_w, bar_top + frac_bot * bar_h),
                ),
                0.0,
                color,
            );
        }

        painter.rect_stroke(
            egui::Rect::from_min_size(
                egui::Pos2::new(bar_x, bar_top),
                egui::Vec2::new(bar_w, bar_h),
            ),
            0.0,
            egui::Stroke::new(1.0, grid_color),
            egui::StrokeKind::Outside,
        );

        let text_x = bar_x + bar_w + 4.0;
        painter.text(
            egui::Pos2::new(text_x, bar_top),
            egui::Align2::LEFT_TOP,
            format!("{:.3}", zmax),
            label_font.clone(),
            egui::Color32::WHITE,
        );
        painter.text(
            egui::Pos2::new(text_x, bar_top + bar_h * 0.5),
            egui::Align2::LEFT_CENTER,
            format!("{:.3}", (zmin + zmax) / 2.0),
            label_font.clone(),
            egui::Color32::WHITE,
        );
        painter.text(
            egui::Pos2::new(text_x, bar_top + bar_h),
            egui::Align2::LEFT_BOTTOM,
            format!("{:.3}", zmin),
            label_font,
            egui::Color32::WHITE,
        );

        // Hover tooltip
        if let Some(hover_pos) = response.hover_pos() {
            let rx = hover_pos.x - plot_origin.x;
            let ry = hover_pos.y - plot_origin.y;
            if rx >= 0.0 && rx < plot_w && ry >= 0.0 && ry < plot_h {
                let ix = (rx / cell_w) as usize;
                let iy = ny - 1 - (ry / cell_h) as usize;
                if ix < nx && iy < ny {
                    let xv = surface.x_values.get(ix).copied().unwrap_or(0.0);
                    let yv = surface.y_values.get(iy).copied().unwrap_or(0.0);
                    let zv = surface
                        .z_grid
                        .get(iy)
                        .and_then(|r| r.get(ix))
                        .copied()
                        .unwrap_or(0.0);
                    egui::show_tooltip_at_pointer(
                        ui.ctx(),
                        response.layer_id,
                        egui::Id::new("surface_tooltip"),
                        |ui| {
                            ui.label(format!("X: {:.4}  Y: {:.4}  Z: {:.6}", xv, yv, zv));
                        },
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

/// Blue-white-red diverging colormap for heatmap display.
fn value_to_heatmap_color(value: f64, min: f64, max: f64) -> egui::Color32 {
    let t = ((value - min) / (max - min)).clamp(0.0, 1.0) as f32;

    // Blue (0) → White (0.5) → Red (1.0)
    let (r, g, b) = if t < 0.5 {
        let s = t * 2.0; // 0..1 within blue→white
        (
            (s * 255.0) as u8,
            (s * 255.0) as u8,
            255u8,
        )
    } else {
        let s = (t - 0.5) * 2.0; // 0..1 within white→red
        (
            255u8,
            ((1.0 - s) * 255.0) as u8,
            ((1.0 - s) * 255.0) as u8,
        )
    };
    egui::Color32::from_rgb(r, g, b)
}
