use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::analysis::AnalysisSetup;
use crate::boundary::Boundary;
use crate::excitation::Excitation;
use crate::field_overlay::FieldOverlay;
use crate::geometry::Geometry;
use crate::material::MaterialDef;
use crate::mesh::MeshOperation;
use crate::net::Net;
use crate::optimetrics::OptimetricsSetup;
use crate::output_variable::OutputVariable;
use crate::radiation::RadiationSetup;
use crate::report::Report;
use crate::solution_index::SolutionIndex;
use crate::solution_type::SolutionType;
use crate::variable::Variable;

/// A single simulation design within a project.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Design {
    pub id: String,
    pub name: String,
    pub solution_type: SolutionType,
    #[serde(default = "default_units")]
    pub units: String,
    #[serde(default)]
    pub design_settings: DesignSettings,
    /// Design-level variables (no `$` prefix).
    #[serde(default)]
    pub local_variables: HashMap<String, Variable>,
    #[serde(default)]
    pub definitions: Definitions,
    #[serde(default)]
    pub geometry: Geometry,
    #[serde(default)]
    pub boundaries: Vec<Boundary>,
    #[serde(default)]
    pub excitations: Vec<Excitation>,
    /// Q3D network definitions.
    #[serde(default)]
    pub nets: Vec<Net>,
    #[serde(default)]
    pub mesh_operations: Vec<MeshOperation>,
    #[serde(default)]
    pub analysis_setups: Vec<AnalysisSetup>,
    /// HFSS radiation settings (far-field / near-field).
    #[serde(default)]
    pub radiation: RadiationSetup,
    #[serde(default)]
    pub output_variables: Vec<OutputVariable>,
    #[serde(default)]
    pub field_overlays: Vec<FieldOverlay>,
    #[serde(default)]
    pub optimetrics: Vec<OptimetricsSetup>,
    #[serde(default)]
    pub reports: Vec<Report>,
    #[serde(default)]
    pub solution_index: SolutionIndex,
}

fn default_units() -> String {
    "mm".to_string()
}

/// Centralized reusable definitions (definition-reference architecture).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Definitions {
    #[serde(default)]
    pub materials: Vec<MaterialDef>,
    #[serde(default)]
    pub coordinate_systems: Vec<CoordinateSystem>,
    #[serde(default)]
    pub named_selections: Vec<NamedSelection>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoordinateSystem {
    pub name: String,
    /// `"Cartesian"`, `"Cylindrical"`, `"Spherical"`.
    #[serde(default = "default_cartesian")]
    pub coordinate_type: String,
    #[serde(default)]
    pub origin: [f64; 3],
    #[serde(default = "default_x_axis")]
    pub x_axis: [f64; 3],
    #[serde(default = "default_y_axis")]
    pub y_axis: [f64; 3],
}

fn default_cartesian() -> String {
    "Cartesian".to_string()
}
fn default_x_axis() -> [f64; 3] {
    [1.0, 0.0, 0.0]
}
fn default_y_axis() -> [f64; 3] {
    [0.0, 1.0, 0.0]
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NamedSelection {
    pub name: String,
    pub selection_type: SelectionType,
    /// References to objects, faces, edges, or vertices.
    #[serde(default)]
    pub selection: Vec<serde_json::Value>,
    #[serde(default)]
    pub description: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum SelectionType {
    Face,
    Edge,
    Vertex,
    #[default]
    Object,
}

/// Design-level settings (port normalization, solver options, etc.).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DesignSettings {
    #[serde(default)]
    pub port_impedance_normalization: PortNormalization,
    #[serde(default)]
    pub deembedding: DeembeddingSettings,
    #[serde(default = "default_s_matrix_type")]
    pub s_matrix_type: String,
    #[serde(default = "default_temperature")]
    pub environment_temperature: String,
    #[serde(default)]
    pub model_validation: ValidationSettings,
    #[serde(default)]
    pub solver_options: SolverOptions,
}

fn default_s_matrix_type() -> String {
    "Modal".to_string()
}
fn default_temperature() -> String {
    "22cel".to_string()
}

impl Default for DesignSettings {
    fn default() -> Self {
        Self {
            port_impedance_normalization: PortNormalization::default(),
            deembedding: DeembeddingSettings::default(),
            s_matrix_type: default_s_matrix_type(),
            environment_temperature: default_temperature(),
            model_validation: ValidationSettings::default(),
            solver_options: SolverOptions::default(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PortNormalization {
    #[serde(default = "default_true")]
    pub enabled: bool,
    #[serde(default = "default_ref_impedance")]
    pub reference_impedance: String,
}

fn default_true() -> bool {
    true
}
fn default_ref_impedance() -> String {
    "50ohm".to_string()
}

impl Default for PortNormalization {
    fn default() -> Self {
        Self {
            enabled: true,
            reference_impedance: default_ref_impedance(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DeembeddingSettings {
    #[serde(default)]
    pub enabled: bool,
    #[serde(default = "default_zero_mm")]
    pub default_distance: String,
}

fn default_zero_mm() -> String {
    "0mm".to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationSettings {
    #[serde(default = "default_true")]
    pub validate_before_solve: bool,
    #[serde(default = "default_true")]
    pub check_intersections: bool,
    #[serde(default = "default_true")]
    pub check_duplicate_boundaries: bool,
    #[serde(default = "default_true")]
    pub check_port_on_boundary: bool,
}

impl Default for ValidationSettings {
    fn default() -> Self {
        Self {
            validate_before_solve: true,
            check_intersections: true,
            check_duplicate_boundaries: true,
            check_port_on_boundary: true,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolverOptions {
    #[serde(default)]
    pub use_shell_elements: bool,
    #[serde(default = "default_curved_order")]
    pub curved_elements_order: String,
    #[serde(default = "default_true")]
    pub allow_solver_fallback: bool,
}

fn default_curved_order() -> String {
    "1st".to_string()
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            use_shell_elements: false,
            curved_elements_order: default_curved_order(),
            allow_solver_fallback: true,
        }
    }
}
