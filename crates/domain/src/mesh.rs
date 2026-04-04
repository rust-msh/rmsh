use serde::{Deserialize, Serialize};

use crate::boundary::Assignment;

/// A mesh refinement operation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshOperation {
    pub name: String,
    pub mesh_type: MeshOperationType,
    pub assignment: Assignment,
    /// Type-specific parameters (max element size, skin depth layers, etc.).
    #[serde(default)]
    pub properties: serde_json::Value,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MeshOperationType {
    /// Maximum element length constraint.
    LengthBased,
    /// Skin depth layer count control.
    SkinDepth,
    /// Surface normal deviation control.
    CurvatureBased,
    /// Global minimum feature size.
    ModelResolution,
}
