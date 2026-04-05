// ---------------------------------------------------------------------------
// ResultDataStore — Unified result data management for Milestone 6
// ---------------------------------------------------------------------------

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use thiserror::Error;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug, Error)]
pub enum ResultError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON parse error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("Result file not found: {0}")]
    NotFound(String),
    #[error("Invalid result data: {0}")]
    InvalidData(String),
}

// ---------------------------------------------------------------------------
// S-Parameter data (from s_parameters.json)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SParameterData {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
    pub setup: String,
    #[serde(default)]
    pub sweep: String,
    pub solution_type: String,
    pub reference_impedance_ohm: f64,
    pub num_ports: usize,
    pub port_names: Vec<String>,
    pub num_frequencies: usize,
    pub frequency_unit: String,
    pub data_format: String,
    pub frequencies: Vec<f64>,
    /// Keyed by e.g. "S11", "S21" etc. Each contains real/imag arrays.
    pub parameters: HashMap<String, SParamValues>,
    /// Pre-computed derived quantities (optional).
    #[serde(default)]
    pub derived: HashMap<String, Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SParamValues {
    pub real: Vec<f64>,
    pub imag: Vec<f64>,
}

impl SParameterData {
    /// Get complex S-parameter values for a given key (e.g. "S11").
    pub fn get_complex(&self, key: &str) -> Option<Vec<[f64; 2]>> {
        self.parameters.get(key).map(|v| {
            v.real
                .iter()
                .zip(v.imag.iter())
                .map(|(&re, &im)| [re, im])
                .collect()
        })
    }

    /// Build the parameter key string for S(row, col), 1-based indexing.
    pub fn s_key(row: usize, col: usize) -> String {
        format!("S{}{}", row, col)
    }

    /// Get frequency values converted to Hz.
    pub fn frequencies_hz(&self) -> Vec<f64> {
        let mult = match self.frequency_unit.to_uppercase().as_str() {
            "GHZ" => 1e9,
            "MHZ" => 1e6,
            "KHZ" => 1e3,
            _ => 1.0,
        };
        self.frequencies.iter().map(|&f| f * mult).collect()
    }

    /// Get frequency values in the original unit.
    pub fn frequencies_original(&self) -> &[f64] {
        &self.frequencies
    }
}

// ---------------------------------------------------------------------------
// Convergence data (from convergence.json)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceData {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
    pub setup: String,
    #[serde(default)]
    pub solution_frequency: String,
    #[serde(default)]
    pub target_max_delta_s: Option<f64>,
    #[serde(default)]
    pub target_max_delta_energy: Option<f64>,
    pub passes: Vec<ConvergencePass>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergencePass {
    pub pass_number: u32,
    #[serde(default)]
    pub timestamp: String,
    pub mesh: MeshStats,
    pub solution: SolutionStats,
    #[serde(default)]
    pub performance: Option<PerformanceStats>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshStats {
    pub num_tetrahedra: u64,
    #[serde(default)]
    pub num_nodes: u64,
    #[serde(default)]
    pub mean_edge_length_mm: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolutionStats {
    pub max_delta_s: Option<f64>,
    #[serde(default)]
    pub max_delta_energy: Option<f64>,
    #[serde(default)]
    pub matrix_size: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceStats {
    #[serde(default)]
    pub mesh_time_sec: Option<f64>,
    #[serde(default)]
    pub solve_time_sec: Option<f64>,
    #[serde(default)]
    pub peak_memory_mb: Option<f64>,
}

impl ConvergenceData {
    /// Extract delta S convergence curve: pass_number vs max_delta_s.
    pub fn delta_s_curve(&self) -> Vec<[f64; 2]> {
        self.passes
            .iter()
            .filter_map(|p| {
                p.solution
                    .max_delta_s
                    .map(|ds| [p.pass_number as f64, ds])
            })
            .collect()
    }

    /// Extract delta energy convergence curve.
    pub fn delta_energy_curve(&self) -> Vec<[f64; 2]> {
        self.passes
            .iter()
            .filter_map(|p| {
                p.solution
                    .max_delta_energy
                    .map(|de| [p.pass_number as f64, de])
            })
            .collect()
    }

    /// Extract mesh tetrahedra count curve: pass_number vs num_tetrahedra.
    pub fn tetrahedra_curve(&self) -> Vec<[f64; 2]> {
        self.passes
            .iter()
            .map(|p| [p.pass_number as f64, p.mesh.num_tetrahedra as f64])
            .collect()
    }
}

// ---------------------------------------------------------------------------
// RLCG matrix data (from rlcg_matrix.json)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RlcgMatrixData {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
    pub setup: String,
    pub solution_type: String,
    pub num_nets: usize,
    pub net_names: Vec<String>,
    #[serde(default)]
    pub terminal_names: Vec<String>,
    pub frequencies: Vec<f64>,
    pub matrices: HashMap<String, RlcgMatrix>,
    #[serde(default)]
    pub dc_data: Option<HashMap<String, RlcgDcMatrix>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RlcgMatrix {
    #[serde(default)]
    pub description: String,
    #[serde(default)]
    pub unit: String,
    pub data_per_frequency: Vec<RlcgFrequencyPoint>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RlcgFrequencyPoint {
    pub frequency: f64,
    pub matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RlcgDcMatrix {
    pub matrix: Vec<Vec<f64>>,
}

impl RlcgMatrixData {
    /// Extract a specific matrix element vs frequency curve.
    /// `matrix_type` is "R", "L", "C", or "G"; row/col are 0-based.
    pub fn element_curve(
        &self,
        matrix_type: &str,
        row: usize,
        col: usize,
    ) -> Option<Vec<[f64; 2]>> {
        let mat = self.matrices.get(matrix_type)?;
        Some(
            mat.data_per_frequency
                .iter()
                .filter_map(|fp| {
                    fp.matrix
                        .get(row)
                        .and_then(|r| r.get(col))
                        .map(|&val| [fp.frequency, val])
                })
                .collect(),
        )
    }

    /// Get matrix at a specific frequency index.
    pub fn matrix_at_frequency(
        &self,
        matrix_type: &str,
        freq_idx: usize,
    ) -> Option<&Vec<Vec<f64>>> {
        let mat = self.matrices.get(matrix_type)?;
        mat.data_per_frequency
            .get(freq_idx)
            .map(|fp| &fp.matrix)
    }
}

// ---------------------------------------------------------------------------
// Far-field data (from far_field_*.json)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FarFieldData {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
    pub setup: String,
    #[serde(default)]
    pub far_field_setup: String,
    #[serde(default)]
    pub frequency: String,
    pub theta: AngleRange,
    pub phi: AngleRange,
    #[serde(default)]
    pub fields: HashMap<String, ComplexFieldData>,
    #[serde(default)]
    pub derived_quantities: HashMap<String, DerivedQuantity>,
    #[serde(default)]
    pub antenna_parameters: Option<AntennaParameters>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AngleRange {
    pub start_deg: f64,
    pub stop_deg: f64,
    pub step_deg: f64,
    pub num_points: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComplexFieldData {
    pub data_real: Vec<f64>,
    pub data_imag: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DerivedQuantity {
    #[serde(default)]
    pub unit: String,
    pub data: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AntennaParameters {
    #[serde(default)]
    pub peak_gain_dbi: Option<f64>,
    #[serde(default)]
    pub peak_gain_theta_deg: Option<f64>,
    #[serde(default)]
    pub peak_gain_phi_deg: Option<f64>,
    #[serde(default)]
    pub radiation_efficiency: Option<f64>,
    #[serde(default)]
    pub beamwidth_e_plane_deg: Option<f64>,
    #[serde(default)]
    pub beamwidth_h_plane_deg: Option<f64>,
}

impl FarFieldData {
    /// Get a derived quantity at a fixed phi, varying theta.
    /// Returns Vec<[theta_deg, value]>.
    pub fn theta_cut(&self, quantity: &str, phi_idx: usize) -> Option<Vec<[f64; 2]>> {
        let dq = self.derived_quantities.get(quantity)?;
        let n_theta = self.theta.num_points;
        let n_phi = self.phi.num_points;

        if dq.data.len() != n_theta * n_phi {
            return None;
        }

        let mut result = Vec::with_capacity(n_theta);
        for i_theta in 0..n_theta {
            let theta = self.theta.start_deg + i_theta as f64 * self.theta.step_deg;
            let idx = i_theta * n_phi + phi_idx;
            if let Some(&val) = dq.data.get(idx) {
                result.push([theta, val]);
            }
        }
        Some(result)
    }

    /// Get a derived quantity at a fixed theta, varying phi.
    /// Returns Vec<[phi_deg, value]>.
    pub fn phi_cut(&self, quantity: &str, theta_idx: usize) -> Option<Vec<[f64; 2]>> {
        let dq = self.derived_quantities.get(quantity)?;
        let n_phi = self.phi.num_points;
        let n_theta = self.theta.num_points;

        if dq.data.len() != n_theta * n_phi {
            return None;
        }

        let mut result = Vec::with_capacity(n_phi);
        for i_phi in 0..n_phi {
            let phi = self.phi.start_deg + i_phi as f64 * self.phi.step_deg;
            let idx = theta_idx * n_phi + i_phi;
            if let Some(&val) = dq.data.get(idx) {
                result.push([phi, val]);
            }
        }
        Some(result)
    }
}

// ---------------------------------------------------------------------------
// Cached result enum
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub enum CachedResult {
    SParameter(SParameterData),
    Convergence(ConvergenceData),
    RlcgMatrix(RlcgMatrixData),
    FarField(FarFieldData),
}

// ---------------------------------------------------------------------------
// ResultDataStore
// ---------------------------------------------------------------------------

/// Unified result data store. Loads and caches simulation results from the
/// `.emsp.results/` directory alongside the project file.
pub struct ResultDataStore {
    base_path: PathBuf,
    cache: HashMap<String, CachedResult>,
}

impl ResultDataStore {
    /// Create a new store rooted at the given results directory.
    pub fn new(base_path: PathBuf) -> Self {
        Self {
            base_path,
            cache: HashMap::new(),
        }
    }

    /// Create a store from a project file path (appends `.results/`).
    pub fn from_project_path(project_path: &Path) -> Self {
        let results_dir = project_path.with_extension("emsp.results");
        Self::new(results_dir)
    }

    /// Clear all cached results.
    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }

    /// Base results directory path.
    pub fn base_path(&self) -> &Path {
        &self.base_path
    }

    // -----------------------------------------------------------------------
    // S-Parameter loading
    // -----------------------------------------------------------------------

    /// Load S-parameter data from a JSON file.
    pub fn load_s_parameters(&mut self, relative_path: &str) -> Result<&SParameterData, ResultError> {
        if !self.cache.contains_key(relative_path) {
            let full_path = self.base_path.join(relative_path);
            let data = load_s_parameters_from_file(&full_path)?;
            self.cache
                .insert(relative_path.to_string(), CachedResult::SParameter(data));
        }
        match self.cache.get(relative_path) {
            Some(CachedResult::SParameter(d)) => Ok(d),
            _ => Err(ResultError::NotFound(relative_path.to_string())),
        }
    }

    // -----------------------------------------------------------------------
    // Convergence loading
    // -----------------------------------------------------------------------

    /// Load convergence history from a JSON file.
    pub fn load_convergence(&mut self, relative_path: &str) -> Result<&ConvergenceData, ResultError> {
        if !self.cache.contains_key(relative_path) {
            let full_path = self.base_path.join(relative_path);
            let data = load_convergence_from_file(&full_path)?;
            self.cache
                .insert(relative_path.to_string(), CachedResult::Convergence(data));
        }
        match self.cache.get(relative_path) {
            Some(CachedResult::Convergence(d)) => Ok(d),
            _ => Err(ResultError::NotFound(relative_path.to_string())),
        }
    }

    // -----------------------------------------------------------------------
    // RLCG matrix loading
    // -----------------------------------------------------------------------

    /// Load RLCG matrix data from a JSON file.
    pub fn load_rlcg_matrix(&mut self, relative_path: &str) -> Result<&RlcgMatrixData, ResultError> {
        if !self.cache.contains_key(relative_path) {
            let full_path = self.base_path.join(relative_path);
            let data = load_rlcg_matrix_from_file(&full_path)?;
            self.cache
                .insert(relative_path.to_string(), CachedResult::RlcgMatrix(data));
        }
        match self.cache.get(relative_path) {
            Some(CachedResult::RlcgMatrix(d)) => Ok(d),
            _ => Err(ResultError::NotFound(relative_path.to_string())),
        }
    }

    // -----------------------------------------------------------------------
    // Far-field loading
    // -----------------------------------------------------------------------

    /// Load far-field data from a JSON file.
    pub fn load_far_field(&mut self, relative_path: &str) -> Result<&FarFieldData, ResultError> {
        if !self.cache.contains_key(relative_path) {
            let full_path = self.base_path.join(relative_path);
            let data = load_far_field_from_file(&full_path)?;
            self.cache
                .insert(relative_path.to_string(), CachedResult::FarField(data));
        }
        match self.cache.get(relative_path) {
            Some(CachedResult::FarField(d)) => Ok(d),
            _ => Err(ResultError::NotFound(relative_path.to_string())),
        }
    }

    // -----------------------------------------------------------------------
    // Direct data insertion (for testing or in-memory results)
    // -----------------------------------------------------------------------

    /// Insert a cached result directly.
    pub fn insert(&mut self, key: String, result: CachedResult) {
        self.cache.insert(key, result);
    }

    /// Check if a result is cached.
    pub fn contains(&self, key: &str) -> bool {
        self.cache.contains_key(key)
    }
}

// ---------------------------------------------------------------------------
// Standalone file loaders
// ---------------------------------------------------------------------------

pub fn load_s_parameters_from_file(path: &Path) -> Result<SParameterData, ResultError> {
    let contents = std::fs::read_to_string(path)?;
    let data: SParameterData = serde_json::from_str(&contents)?;
    Ok(data)
}

pub fn load_convergence_from_file(path: &Path) -> Result<ConvergenceData, ResultError> {
    let contents = std::fs::read_to_string(path)?;
    let data: ConvergenceData = serde_json::from_str(&contents)?;
    Ok(data)
}

pub fn load_rlcg_matrix_from_file(path: &Path) -> Result<RlcgMatrixData, ResultError> {
    let contents = std::fs::read_to_string(path)?;
    let data: RlcgMatrixData = serde_json::from_str(&contents)?;
    Ok(data)
}

pub fn load_far_field_from_file(path: &Path) -> Result<FarFieldData, ResultError> {
    let contents = std::fs::read_to_string(path)?;
    let data: FarFieldData = serde_json::from_str(&contents)?;
    Ok(data)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_s_parameter_json() {
        let json = r#"{
            "format_version": "1.0",
            "file_type": "SParameterData",
            "design_id": "design-001",
            "setup": "Setup1",
            "sweep": "Sweep1",
            "solution_type": "DrivenModal",
            "reference_impedance_ohm": 50.0,
            "num_ports": 1,
            "port_names": ["Port1"],
            "num_frequencies": 3,
            "frequency_unit": "GHz",
            "data_format": "RealImaginary",
            "frequencies": [1.0, 2.0, 3.0],
            "parameters": {
                "S11": {
                    "real": [-0.85, -0.72, -0.50],
                    "imag": [0.12, 0.25, 0.40]
                }
            },
            "derived": {}
        }"#;

        let data: SParameterData = serde_json::from_str(json).unwrap();
        assert_eq!(data.num_ports, 1);
        assert_eq!(data.num_frequencies, 3);
        assert_eq!(data.frequencies.len(), 3);

        let complex = data.get_complex("S11").unwrap();
        assert_eq!(complex.len(), 3);
        assert!((complex[0][0] - (-0.85)).abs() < 1e-10);
        assert!((complex[0][1] - 0.12).abs() < 1e-10);

        let freqs_hz = data.frequencies_hz();
        assert!((freqs_hz[0] - 1e9).abs() < 1.0);
        assert!((freqs_hz[2] - 3e9).abs() < 1.0);
    }

    #[test]
    fn parse_convergence_json() {
        let json = r#"{
            "format_version": "1.0",
            "file_type": "ConvergenceHistory",
            "design_id": "design-001",
            "setup": "Setup1",
            "solution_frequency": "2.4GHz",
            "target_max_delta_s": 0.02,
            "passes": [
                {
                    "pass_number": 1,
                    "timestamp": "2026-04-04T14:25:15Z",
                    "mesh": {
                        "num_tetrahedra": 5420,
                        "num_nodes": 1210
                    },
                    "solution": {
                        "max_delta_s": null,
                        "matrix_size": 16260
                    }
                },
                {
                    "pass_number": 2,
                    "timestamp": "2026-04-04T14:25:30Z",
                    "mesh": {
                        "num_tetrahedra": 8100,
                        "num_nodes": 1800
                    },
                    "solution": {
                        "max_delta_s": 0.15,
                        "matrix_size": 24300
                    }
                }
            ]
        }"#;

        let data: ConvergenceData = serde_json::from_str(json).unwrap();
        assert_eq!(data.passes.len(), 2);

        let delta_s = data.delta_s_curve();
        // First pass has null delta_s, so only pass 2 is included
        assert_eq!(delta_s.len(), 1);
        assert!((delta_s[0][0] - 2.0).abs() < 1e-10);
        assert!((delta_s[0][1] - 0.15).abs() < 1e-10);

        let tet_curve = data.tetrahedra_curve();
        assert_eq!(tet_curve.len(), 2);
        assert!((tet_curve[0][1] - 5420.0).abs() < 1e-10);
    }

    #[test]
    fn parse_rlcg_matrix_json() {
        let json = r#"{
            "format_version": "1.0",
            "file_type": "RLCGMatrixData",
            "design_id": "design-002",
            "setup": "Q3D_Setup1",
            "solution_type": "Q3D_ACRL",
            "num_nets": 2,
            "net_names": ["Signal1", "Ground"],
            "frequencies": [0.01, 1.0],
            "matrices": {
                "R": {
                    "description": "Resistance matrix",
                    "unit": "ohm",
                    "data_per_frequency": [
                        { "frequency": 0.01, "matrix": [[0.125, 0.003], [0.003, 0.125]] },
                        { "frequency": 1.0, "matrix": [[0.250, 0.006], [0.006, 0.250]] }
                    ]
                }
            }
        }"#;

        let data: RlcgMatrixData = serde_json::from_str(json).unwrap();
        assert_eq!(data.num_nets, 2);

        let curve = data.element_curve("R", 0, 0).unwrap();
        assert_eq!(curve.len(), 2);
        assert!((curve[0][1] - 0.125).abs() < 1e-10);
        assert!((curve[1][1] - 0.250).abs() < 1e-10);

        let mat = data.matrix_at_frequency("R", 0).unwrap();
        assert_eq!(mat.len(), 2);
        assert!((mat[0][0] - 0.125).abs() < 1e-10);
    }

    #[test]
    fn s_param_key_formatting() {
        assert_eq!(SParameterData::s_key(1, 1), "S11");
        assert_eq!(SParameterData::s_key(2, 1), "S21");
    }

    #[test]
    fn result_data_store_insert_and_retrieve() {
        let mut store = ResultDataStore::new(PathBuf::from("/tmp/test"));
        let data = ConvergenceData {
            format_version: "1.0".into(),
            file_type: "ConvergenceHistory".into(),
            design_id: "d1".into(),
            setup: "Setup1".into(),
            solution_frequency: String::new(),
            target_max_delta_s: Some(0.02),
            target_max_delta_energy: None,
            passes: vec![],
        };
        store.insert("conv".into(), CachedResult::Convergence(data));
        assert!(store.contains("conv"));
    }
}
