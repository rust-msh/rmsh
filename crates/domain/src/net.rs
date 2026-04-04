use serde::{Deserialize, Serialize};

use crate::boundary::Assignment;
use crate::excitation::ExcitationType;

/// A Q3D electrical network definition.
///
/// Groups conductor objects into named nets with source/sink terminal pairs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Net {
    pub name: String,
    /// Geometry object names belonging to this net.
    #[serde(default)]
    pub objects: Vec<String>,
    /// Whether this net is the ground reference.
    #[serde(default)]
    pub is_ground_reference: bool,
    #[serde(default)]
    pub terminals: Vec<Terminal>,
}

/// A source or sink terminal within a Q3D net.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Terminal {
    pub name: String,
    pub terminal_type: ExcitationType,
    pub assignment: Assignment,
}
