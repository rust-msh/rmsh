use serde::{Deserialize, Serialize};

/// Optimetrics setup (parametric sweep, optimization, sensitivity, etc.).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum OptimetricsSetup {
    ParametricSweep {
        name: String,
        #[serde(default = "default_true")]
        enabled: bool,
        setup: String,
        #[serde(default)]
        sweep_definitions: Vec<SweepDefinition>,
        #[serde(default)]
        constraints: Vec<serde_json::Value>,
        #[serde(default)]
        goals: Vec<serde_json::Value>,
    },
    Optimization {
        name: String,
        #[serde(default = "default_true")]
        enabled: bool,
        setup: String,
        /// `"QuasiNewton"`, `"PatternSearch"`, `"GeneticAlgorithm"`, `"SNLP"`.
        #[serde(default = "default_quasi_newton")]
        algorithm: String,
        #[serde(default = "default_max_iterations")]
        max_iterations: u32,
        #[serde(default)]
        variables: Vec<OptimizationVariable>,
        #[serde(default)]
        goals: Vec<OptimizationGoal>,
        #[serde(default)]
        constraints: Vec<serde_json::Value>,
    },
    Sensitivity {
        name: String,
        #[serde(default = "default_true")]
        enabled: bool,
        setup: String,
        #[serde(default)]
        variables: Vec<SensitivityVariable>,
        output: String,
        #[serde(default = "default_num_samples")]
        num_samples: u32,
    },
    Statistical {
        name: String,
        #[serde(default = "default_true")]
        enabled: bool,
        setup: String,
        #[serde(default)]
        variables: Vec<SensitivityVariable>,
        #[serde(default = "default_num_trials")]
        num_trials: u32,
    },
    Tuning {
        name: String,
        #[serde(default = "default_true")]
        enabled: bool,
        setup: String,
        #[serde(default)]
        variables: Vec<String>,
    },
}

fn default_true() -> bool {
    true
}
fn default_quasi_newton() -> String {
    "QuasiNewton".to_string()
}
fn default_max_iterations() -> u32 {
    100
}
fn default_num_samples() -> u32 {
    50
}
fn default_num_trials() -> u32 {
    1000
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SweepDefinition {
    pub variable: String,
    /// `"LinearStep"`, `"LinearCount"`, `"LogScale"`, `"DiscreteList"`.
    pub sweep_type: String,
    #[serde(default)]
    pub params: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationVariable {
    pub variable: String,
    pub min: String,
    pub max: String,
    pub starting: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationGoal {
    pub name: String,
    pub expression: String,
    /// `"Minimize"`, `"Maximize"`, `"LessThan"`, `"GreaterThan"`, `"EqualTo"`.
    pub condition: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub frequency_range: Option<FrequencyRange>,
    #[serde(default = "default_weight")]
    pub weight: f64,
}

fn default_weight() -> f64 {
    1.0
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrequencyRange {
    pub start: String,
    pub stop: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SensitivityVariable {
    pub variable: String,
    /// e.g., `"5%"`, `"±10%"`.
    pub variation: String,
    /// `"Uniform"`, `"Gaussian"`, `"LogNormal"`.
    #[serde(default = "default_uniform")]
    pub distribution: String,
}

fn default_uniform() -> String {
    "Uniform".to_string()
}
