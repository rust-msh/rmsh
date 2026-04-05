use std::collections::HashMap;

use serde::{Deserialize, Serialize};

/// A post-processing report (2D chart, data table, etc.).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Report {
    pub name: String,
    pub category: ReportCategory,
    pub chart_type: ChartType,
    /// Analysis setup name.
    pub solution: String,
    pub domain: ReportDomain,
    #[serde(default)]
    pub traces: Vec<ReportTrace>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub x_axis: Option<AxisConfig>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub y_axis: Option<AxisConfig>,
    #[serde(default)]
    pub markers: Vec<ReportMarker>,
    #[serde(default)]
    pub limit_lines: Vec<LimitLine>,
    /// Far-field setup name (HFSS FarField reports only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub far_field_setup: Option<String>,
    /// Matrix type for Q3D reports: `"L"`, `"R"`, `"C"`, `"G"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub matrix_type: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub display_options: Option<DisplayOptions>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReportCategory {
    SParameter,
    FarField,
    NearField,
    Fields,
    Eigenmode,
    Emission,
    RLCGMatrix,
    Q3DFields,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum ChartType {
    #[default]
    Rectangular,
    Polar,
    Smith,
    DataTable,
    Polar3D,
    MatrixTable,
    /// 3D rectangular plot for parameter sweep results.
    Rectangular3D,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportDomain {
    /// `"Frequency"`, `"Time"`, `"Theta"`, `"Phi"`, etc.
    pub domain_type: String,
    pub primary_sweep: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fixed_values: Option<HashMap<String, String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportTrace {
    pub name: String,
    /// Expression, e.g. `"dB(S(1,1))"`, `"GainTotal"`.
    pub expression: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub style: Option<TraceStyle>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parametric_values: Option<HashMap<String, String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fixed_values: Option<HashMap<String, String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraceStyle {
    #[serde(default = "default_trace_color")]
    pub color: [u8; 3],
    #[serde(default = "default_line_width")]
    pub line_width: u32,
    #[serde(default = "default_line_style")]
    pub line_style: String,
}

fn default_trace_color() -> [u8; 3] {
    [0, 0, 255]
}
fn default_line_width() -> u32 {
    2
}
fn default_line_style() -> String {
    "Solid".to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AxisConfig {
    pub label: String,
    #[serde(default)]
    pub unit: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub auto_range: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scale: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportMarker {
    pub name: String,
    pub trace: String,
    pub x_value: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LimitLine {
    pub name: String,
    pub y_value: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub style: Option<TraceStyle>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DisplayOptions {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub heatmap_enabled: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub show_unit: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub decimal_places: Option<u32>,
}
