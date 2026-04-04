use std::collections::HashMap;

use serde::{Deserialize, Serialize};

/// Tracks the solve state of all analysis setups in a design.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SolutionIndex {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub last_solve_time: Option<String>,
    #[serde(default)]
    pub setups: HashMap<String, SetupSolutionStatus>,
    #[serde(default)]
    pub optimetrics: HashMap<String, OptimetricsSolutionStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exports: Option<HashMap<String, ExportInfo>>,
    /// If set, describes why results are stale.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stale_reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SetupSolutionStatus {
    pub status: SolveStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub solved_at: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub converged_pass: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub num_tetrahedra: Option<u64>,
    /// Q3D: number of surface triangles.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub num_triangles: Option<u64>,
    /// Q3D: energy convergence metric.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub final_delta_energy: Option<f64>,
    #[serde(default)]
    pub is_stale: bool,
    #[serde(default)]
    pub solved_variations: HashMap<String, VariationResult>,
    #[serde(default)]
    pub sweeps: HashMap<String, SweepSolutionResult>,
    /// Q3D: RLCG matrix summary at the nominal point.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rlcg_summary: Option<RlcgSummary>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum SolveStatus {
    #[default]
    NotSolved,
    InProgress,
    Converged,
    NotConverged,
    Failed,
    Completed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariationResult {
    #[serde(default)]
    pub variables: HashMap<String, String>,
    pub result_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SweepSolutionResult {
    pub status: SolveStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub num_frequency_points: Option<u32>,
    pub result_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimetricsSolutionStatus {
    pub status: SolveStatus,
    #[serde(default)]
    pub total_variations: u32,
    #[serde(default)]
    pub completed_variations: u32,
    pub result_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RlcgSummary {
    pub num_nets: u32,
    pub num_terminals: u32,
    pub r_max_ohm: f64,
    pub l_max_n_h: f64,
    pub c_max_p_f: f64,
    pub g_max_m_s: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExportInfo {
    pub exported_at: String,
    pub model_type: String,
    pub file_path: String,
}
