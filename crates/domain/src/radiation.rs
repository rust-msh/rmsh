use serde::{Deserialize, Serialize};

/// HFSS radiation setup (far-field and near-field definitions).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct RadiationSetup {
    #[serde(default)]
    pub far_field_setups: Vec<FarFieldSetup>,
    #[serde(default)]
    pub near_field_setups: Vec<NearFieldSetup>,
    #[serde(default)]
    pub antenna_parameters: AntennaParameters,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FarFieldSetup {
    pub name: String,
    /// `"InfiniteSphere"` or `"InfinitePlane"`.
    #[serde(default = "default_infinite_sphere")]
    pub setup_type: String,
    #[serde(default = "default_global_cs")]
    pub coordinate_system: String,
    pub theta: AngleRange,
    pub phi: AngleRange,
    #[serde(default)]
    pub use_custom_radiation_surface: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub radiation_surface: Option<String>,
}

fn default_infinite_sphere() -> String {
    "InfiniteSphere".to_string()
}
fn default_global_cs() -> String {
    "Global".to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AngleRange {
    pub start: String,
    pub stop: String,
    pub step: String,
}

/// Near-field sampling setup (one of several geometric types).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum NearFieldSetup {
    Line {
        name: String,
        start_point: [f64; 3],
        end_point: [f64; 3],
        num_points: u32,
    },
    Rectangle {
        name: String,
        center: [f64; 3],
        width: f64,
        height: f64,
        axis: String,
        num_points_u: u32,
        num_points_v: u32,
    },
    Sphere {
        name: String,
        center: [f64; 3],
        radius: f64,
        num_theta: u32,
        num_phi: u32,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AntennaParameters {
    #[serde(default = "default_ref_impedance")]
    pub reference_impedance: String,
    #[serde(default = "default_true")]
    pub calculate_antenna_params: bool,
}

fn default_ref_impedance() -> String {
    "50ohm".to_string()
}
fn default_true() -> bool {
    true
}

impl Default for AntennaParameters {
    fn default() -> Self {
        Self {
            reference_impedance: default_ref_impedance(),
            calculate_antenna_params: true,
        }
    }
}
