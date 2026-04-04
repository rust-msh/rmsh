use std::collections::HashMap;

use serde::{Deserialize, Serialize};

/// A project- or design-level variable that can hold a constant, expression, or
/// dataset reference.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Variable {
    /// Constant value string, e.g. `"2.4GHz"`, `"28.5mm"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub value: Option<String>,
    /// Mathematical expression, e.g. `"($freq - 0.5GHz) * 2"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expression: Option<String>,
    /// Human-readable description.
    #[serde(default)]
    pub description: String,
    /// Unit category hint: `"Frequency"`, `"Length"`, `"None"`, etc.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit_type: Option<String>,
}

/// How a material property value is specified.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum PropertyValue {
    /// A single constant.
    Constant { value: f64 },
    /// A mathematical expression referencing variables.
    Expression { expression: String },
    /// A dataset lookup (frequency/temperature dependent).
    Dataset {
        dataset: String,
        independent_variable: String,
    },
}

impl Default for PropertyValue {
    fn default() -> Self {
        Self::Constant { value: 0.0 }
    }
}

/// A lookup table for frequency- or temperature-dependent data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetDefinition {
    #[serde(default)]
    pub description: String,
    /// Name of the independent variable, e.g. `"Freq"`, `"Temp"`.
    pub independent_variable: String,
    pub independent_unit: String,
    pub dependent_unit: String,
    pub data: Vec<DataPoint>,
    #[serde(default)]
    pub interpolation: InterpolationType,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DataPoint {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum InterpolationType {
    #[default]
    PiecewiseLinear,
    CubicSpline,
    Debye,
    DjordjevicSarkar,
}

/// Resolved variable context used during expression evaluation.
#[derive(Debug, Clone, Default)]
pub struct VariableContext {
    /// Merged project + design variables with resolved numeric values.
    pub resolved: HashMap<String, f64>,
}
