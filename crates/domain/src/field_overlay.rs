use serde::{Deserialize, Serialize};

/// A 3D field overlay plot definition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldOverlay {
    pub name: String,
    pub quantity: FieldQuantity,
    /// Component selector: `"Mag"`, `"MagX"`, `"Real"`, `"Vector"`, etc.
    #[serde(default = "default_mag")]
    pub component: String,
    pub plot_type: FieldPlotType,
    /// Geometry assignment (faces, objects, cut planes).
    #[serde(default)]
    pub assignment: serde_json::Value,
    /// Analysis setup name.
    pub solution: String,
    /// Frequency point.
    pub frequency: String,
    #[serde(default = "default_zero_deg")]
    pub phase: String,
    #[serde(default)]
    pub scale: FieldScale,
    #[serde(default)]
    pub display: FieldDisplay,
}

fn default_mag() -> String {
    "Mag".to_string()
}
fn default_zero_deg() -> String {
    "0deg".to_string()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FieldQuantity {
    E,
    H,
    Jvol,
    Jsurf,
    SAR,
    Poynting,
    /// Q3D: charge density.
    ChargeDistribution,
    /// Q3D: ohmic loss density.
    OhmicLoss,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum FieldPlotType {
    #[default]
    Surface,
    CutPlane,
    Volume,
    Line,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldScale {
    #[serde(default = "default_linear")]
    pub scale_type: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max: Option<f64>,
    #[serde(default)]
    pub unit: String,
}

fn default_linear() -> String {
    "Linear".to_string()
}

impl Default for FieldScale {
    fn default() -> Self {
        Self {
            scale_type: default_linear(),
            min: None,
            max: None,
            unit: String::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldDisplay {
    /// `"Shaded"`, `"Arrow"`, `"Contour"`.
    #[serde(default = "default_shaded")]
    pub plot_style: String,
    #[serde(default)]
    pub show_arrows: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub arrow_spacing: Option<u32>,
    #[serde(default = "default_num_colors")]
    pub num_colors: u32,
    #[serde(default = "default_opacity")]
    pub opacity: f32,
}

fn default_shaded() -> String {
    "Shaded".to_string()
}
fn default_num_colors() -> u32 {
    256
}
fn default_opacity() -> f32 {
    1.0
}

impl Default for FieldDisplay {
    fn default() -> Self {
        Self {
            plot_style: default_shaded(),
            show_arrows: false,
            arrow_spacing: None,
            num_colors: default_num_colors(),
            opacity: default_opacity(),
        }
    }
}
