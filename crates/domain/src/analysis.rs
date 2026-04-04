use serde::{Deserialize, Serialize};

/// Simulation analysis setup (shared structure for HFSS and Q3D).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisSetup {
    pub name: String,
    #[serde(default = "default_true")]
    pub enabled: bool,
    /// Solution frequency (can reference a variable, e.g. `"$freq"`).
    pub solution_frequency: String,
    #[serde(default = "default_max_passes")]
    pub max_passes: u32,
    /// HFSS convergence criterion (S-parameter delta).
    #[serde(default = "default_delta_s")]
    pub max_delta_s: f64,
    #[serde(default = "default_min_converged")]
    pub min_converged_passes: u32,
    #[serde(default = "default_order_basis")]
    pub order_basis: String,
    #[serde(default = "default_solver_type")]
    pub solver_type: String,
    #[serde(default)]
    pub frequency_sweeps: Vec<FrequencySweep>,

    // Q3D-specific fields
    /// Energy convergence threshold (Q3D).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_delta_energy: Option<f64>,
    /// Adaptive frequency for MoM refinement (Q3D).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub adaptive_frequency: Option<String>,
    /// Percent of elements to refine per pass (Q3D).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub percent_refinement: Option<u32>,
    /// DC resistance/inductance extraction settings (Q3D).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dc_settings: Option<DcSettings>,
}

fn default_true() -> bool {
    true
}
fn default_max_passes() -> u32 {
    15
}
fn default_delta_s() -> f64 {
    0.02
}
fn default_min_converged() -> u32 {
    2
}
fn default_order_basis() -> String {
    "Mixed".to_string()
}
fn default_solver_type() -> String {
    "Direct".to_string()
}

/// Frequency sweep definition within an analysis setup.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrequencySweep {
    pub name: String,
    pub sweep_type: SweepType,
    pub start: String,
    pub stop: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub step: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub count: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scale: Option<String>,
    #[serde(default)]
    pub save_fields: bool,
    #[serde(default)]
    pub save_rad_fields: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum SweepType {
    #[default]
    Discrete,
    Interpolating,
    Fast,
}

/// DC extraction settings for Q3D.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DcSettings {
    #[serde(default)]
    pub compute_dc_resistance: bool,
    #[serde(default)]
    pub compute_dc_inductance: bool,
}
