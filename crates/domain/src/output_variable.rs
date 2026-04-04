use serde::{Deserialize, Serialize};

/// A derived mathematical expression computed from simulation results.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputVariable {
    pub name: String,
    /// Expression referencing simulation quantities, e.g. `"dB(S(Port1,Port1))"`.
    pub expression: String,
    #[serde(default)]
    pub description: String,
}
