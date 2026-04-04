use serde::{Deserialize, Serialize};

use crate::boundary::Assignment;

/// An electromagnetic excitation (port or source).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Excitation {
    pub name: String,
    pub excitation_type: ExcitationType,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assignment: Option<Assignment>,
    /// Type-specific properties (impedance, mode count, etc.).
    #[serde(default)]
    pub properties: serde_json::Value,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExcitationType {
    // HFSS
    WavePort,
    LumpedPort,
    FloquetPort,
    IncidentWave,
    VoltageDrop,

    // Q3D
    Source,
    Sink,
}
