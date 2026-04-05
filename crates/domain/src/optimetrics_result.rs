// ---------------------------------------------------------------------------
// Optimetrics Result I/O — Load/save optimization and sweep summary files
// ---------------------------------------------------------------------------

use std::collections::HashMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::result_store::ResultError;

// ---------------------------------------------------------------------------
// OptimetricsSummary — matches docs/em-result-file-formats.md §2.11
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimetricsSummary {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
    pub optimetrics_name: String,
    #[serde(rename = "type")]
    pub optimetrics_type: String,
    pub setup: String,

    // Timing
    #[serde(default)]
    pub started_at: String,
    #[serde(default)]
    pub finished_at: String,

    // Parametric sweep fields
    #[serde(default)]
    pub swept_variables: Vec<String>,
    #[serde(default)]
    pub output_variables: Vec<String>,
    #[serde(default)]
    pub total_variations: u32,
    #[serde(default)]
    pub completed_variations: u32,
    #[serde(default)]
    pub failed_variations: u32,
    #[serde(default)]
    pub variations: Vec<Variation>,

    // Optimization fields
    #[serde(default)]
    pub algorithm: Option<String>,
    #[serde(default)]
    pub max_iterations: Option<u32>,
    #[serde(default)]
    pub optimized_variables: Vec<String>,
    #[serde(default)]
    pub goals: Vec<serde_json::Value>,
    #[serde(default)]
    pub total_iterations: Option<u32>,
    #[serde(default)]
    pub converged: Option<bool>,
    #[serde(default)]
    pub convergence_history: Vec<ConvergenceStep>,
    #[serde(default)]
    pub best_result: Option<BestResult>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Variation {
    pub index: u32,
    pub variables: HashMap<String, String>,
    pub status: String,
    #[serde(default)]
    pub num_passes: u32,
    #[serde(default)]
    pub outputs: HashMap<String, OutputValue>,
    #[serde(default)]
    pub result_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputValue {
    pub value: f64,
    #[serde(default)]
    pub unit: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceStep {
    pub iteration: u32,
    pub cost: f64,
    #[serde(default)]
    pub variables: HashMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BestResult {
    pub iteration: u32,
    pub variables: HashMap<String, String>,
    pub cost: f64,
    #[serde(default)]
    pub outputs: HashMap<String, OutputValue>,
    #[serde(default)]
    pub result_path: String,
}

// ---------------------------------------------------------------------------
// File I/O
// ---------------------------------------------------------------------------

pub fn load_summary(path: &Path) -> Result<OptimetricsSummary, ResultError> {
    let contents = std::fs::read_to_string(path)?;
    let data: OptimetricsSummary = serde_json::from_str(&contents)?;
    Ok(data)
}

pub fn save_summary(summary: &OptimetricsSummary, path: &Path) -> Result<(), ResultError> {
    let json = serde_json::to_string_pretty(summary).map_err(|e| {
        ResultError::InvalidData(format!("failed to serialize summary: {}", e))
    })?;
    std::fs::write(path, json)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Helper methods
// ---------------------------------------------------------------------------

impl OptimetricsSummary {
    /// Create a new empty parametric sweep summary.
    pub fn new_sweep(
        design_id: &str,
        name: &str,
        setup: &str,
        swept_variables: Vec<String>,
        output_variables: Vec<String>,
    ) -> Self {
        Self {
            format_version: "1.0".into(),
            file_type: "OptimetricsSummary".into(),
            design_id: design_id.into(),
            optimetrics_name: name.into(),
            optimetrics_type: "ParametricSweep".into(),
            setup: setup.into(),
            started_at: String::new(),
            finished_at: String::new(),
            swept_variables,
            output_variables,
            total_variations: 0,
            completed_variations: 0,
            failed_variations: 0,
            variations: Vec::new(),
            algorithm: None,
            max_iterations: None,
            optimized_variables: Vec::new(),
            goals: Vec::new(),
            total_iterations: None,
            converged: None,
            convergence_history: Vec::new(),
            best_result: None,
        }
    }

    /// Create a new empty optimization summary.
    pub fn new_optimization(
        design_id: &str,
        name: &str,
        setup: &str,
        algorithm: &str,
        variables: Vec<String>,
    ) -> Self {
        Self {
            format_version: "1.0".into(),
            file_type: "OptimetricsSummary".into(),
            design_id: design_id.into(),
            optimetrics_name: name.into(),
            optimetrics_type: "Optimization".into(),
            setup: setup.into(),
            started_at: String::new(),
            finished_at: String::new(),
            swept_variables: Vec::new(),
            output_variables: Vec::new(),
            total_variations: 0,
            completed_variations: 0,
            failed_variations: 0,
            variations: Vec::new(),
            algorithm: Some(algorithm.into()),
            max_iterations: Some(100),
            optimized_variables: variables,
            goals: Vec::new(),
            total_iterations: None,
            converged: None,
            convergence_history: Vec::new(),
            best_result: None,
        }
    }

    /// Whether this is a parametric sweep.
    pub fn is_sweep(&self) -> bool {
        self.optimetrics_type == "ParametricSweep"
    }

    /// Whether this is an optimization.
    pub fn is_optimization(&self) -> bool {
        self.optimetrics_type == "Optimization"
    }

    /// Get the cost convergence curve: iteration vs cost.
    pub fn cost_curve(&self) -> Vec<[f64; 2]> {
        self.convergence_history
            .iter()
            .map(|s| [s.iteration as f64, s.cost])
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_sweep_summary() {
        let json = r#"{
            "format_version": "1.0",
            "file_type": "OptimetricsSummary",
            "design_id": "design-001",
            "optimetrics_name": "LengthSweep",
            "type": "ParametricSweep",
            "setup": "Setup1",
            "started_at": "2026-04-04T15:00:00Z",
            "finished_at": "2026-04-04T16:30:00Z",
            "swept_variables": ["patch_l"],
            "output_variables": ["S11_at_center"],
            "total_variations": 2,
            "completed_variations": 2,
            "failed_variations": 0,
            "variations": [
                {
                    "index": 1,
                    "variables": { "patch_l": "25.0mm" },
                    "status": "Converged",
                    "num_passes": 4,
                    "outputs": {
                        "S11_at_center": { "value": -8.2, "unit": "dB" }
                    },
                    "result_path": "variation_001/"
                },
                {
                    "index": 2,
                    "variables": { "patch_l": "25.5mm" },
                    "status": "Converged",
                    "num_passes": 3,
                    "outputs": {
                        "S11_at_center": { "value": -12.5, "unit": "dB" }
                    },
                    "result_path": "variation_002/"
                }
            ]
        }"#;

        let summary: OptimetricsSummary = serde_json::from_str(json).unwrap();
        assert!(summary.is_sweep());
        assert!(!summary.is_optimization());
        assert_eq!(summary.variations.len(), 2);
        assert_eq!(summary.total_variations, 2);
        assert!((summary.variations[0].outputs["S11_at_center"].value - (-8.2)).abs() < 1e-10);
    }

    #[test]
    fn parse_optimization_summary() {
        let json = r#"{
            "format_version": "1.0",
            "file_type": "OptimetricsSummary",
            "design_id": "design-001",
            "optimetrics_name": "MatchOptimize",
            "type": "Optimization",
            "setup": "Setup1",
            "algorithm": "QuasiNewton",
            "max_iterations": 50,
            "optimized_variables": ["patch_l", "patch_w"],
            "total_iterations": 23,
            "converged": true,
            "convergence_history": [
                { "iteration": 1, "cost": -5.2, "variables": { "patch_l": "28.5mm" } },
                { "iteration": 23, "cost": -22.5, "variables": { "patch_l": "29.8mm" } }
            ],
            "best_result": {
                "iteration": 23,
                "variables": { "patch_l": "29.8mm", "patch_w": "35.5mm" },
                "cost": -22.5,
                "outputs": {
                    "S11_at_center": { "value": -22.5, "unit": "dB" }
                },
                "result_path": "iteration_023/"
            }
        }"#;

        let summary: OptimetricsSummary = serde_json::from_str(json).unwrap();
        assert!(summary.is_optimization());
        assert_eq!(summary.converged, Some(true));
        assert_eq!(summary.convergence_history.len(), 2);

        let curve = summary.cost_curve();
        assert_eq!(curve.len(), 2);
        assert!((curve[0][0] - 1.0).abs() < 1e-10);
        assert!((curve[0][1] - (-5.2)).abs() < 1e-10);

        let best = summary.best_result.as_ref().unwrap();
        assert!((best.cost - (-22.5)).abs() < 1e-10);
    }

    #[test]
    fn roundtrip_save_load() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("summary.json");

        let mut summary = OptimetricsSummary::new_sweep(
            "d1",
            "Sweep1",
            "Setup1",
            vec!["x".into()],
            vec!["output".into()],
        );
        summary.total_variations = 1;
        summary.variations.push(Variation {
            index: 1,
            variables: HashMap::from([("x".into(), "5.0".into())]),
            status: "Converged".into(),
            num_passes: 3,
            outputs: HashMap::from([(
                "output".into(),
                OutputValue {
                    value: 42.0,
                    unit: "dB".into(),
                },
            )]),
            result_path: "variation_001/".into(),
        });

        save_summary(&summary, &path).unwrap();
        let loaded = load_summary(&path).unwrap();

        assert_eq!(loaded.optimetrics_name, "Sweep1");
        assert_eq!(loaded.variations.len(), 1);
        assert!((loaded.variations[0].outputs["output"].value - 42.0).abs() < 1e-10);
    }

    #[test]
    fn new_optimization_summary() {
        let summary = OptimetricsSummary::new_optimization(
            "d1",
            "Opt1",
            "Setup1",
            "GeneticAlgorithm",
            vec!["x".into(), "y".into()],
        );
        assert!(summary.is_optimization());
        assert_eq!(summary.algorithm, Some("GeneticAlgorithm".into()));
        assert_eq!(summary.optimized_variables.len(), 2);
    }
}
