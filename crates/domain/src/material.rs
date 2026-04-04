use serde::{Deserialize, Serialize};

use crate::variable::PropertyValue;

/// Full material definition with category and multi-source properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaterialDef {
    pub name: String,
    #[serde(default)]
    pub category: MaterialCategory,
    #[serde(default)]
    pub properties: MaterialProperties,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub appearance: Option<MaterialAppearance>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum MaterialCategory {
    Conductor,
    #[default]
    Dielectric,
    Magnetic,
    Composite,
    Gas,
}

/// Electromagnetic material properties. Each property supports constant /
/// expression / dataset sourcing via `PropertyValue`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaterialProperties {
    #[serde(default = "default_one")]
    pub permittivity: PropertyValue,
    #[serde(default = "default_one")]
    pub permeability: PropertyValue,
    #[serde(default)]
    pub conductivity: PropertyValue,
    #[serde(default)]
    pub dielectric_loss_tangent: PropertyValue,
    #[serde(default)]
    pub magnetic_loss_tangent: PropertyValue,
    #[serde(default)]
    pub mass_density: PropertyValue,
}

impl Default for MaterialProperties {
    fn default() -> Self {
        Self {
            permittivity: PropertyValue::Constant { value: 1.0 },
            permeability: PropertyValue::Constant { value: 1.0 },
            conductivity: PropertyValue::Constant { value: 0.0 },
            dielectric_loss_tangent: PropertyValue::Constant { value: 0.0 },
            magnetic_loss_tangent: PropertyValue::Constant { value: 0.0 },
            mass_density: PropertyValue::Constant { value: 0.0 },
        }
    }
}

fn default_one() -> PropertyValue {
    PropertyValue::Constant { value: 1.0 }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaterialAppearance {
    #[serde(default = "default_color")]
    pub color: [u8; 3],
    #[serde(default)]
    pub transparency: f32,
}

fn default_color() -> [u8; 3] {
    [128, 128, 128]
}

impl Default for MaterialAppearance {
    fn default() -> Self {
        Self {
            color: default_color(),
            transparency: 0.0,
        }
    }
}
