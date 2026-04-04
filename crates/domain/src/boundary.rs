use serde::{Deserialize, Serialize};

/// A boundary condition assigned to geometry faces/edges/objects.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Boundary {
    pub name: String,
    pub boundary_type: BoundaryType,
    pub assignment: Assignment,
    /// Type-specific properties (e.g., impedance value, conductivity).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub properties: Option<serde_json::Value>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoundaryType {
    // HFSS & Q3D shared
    PerfectE,
    PerfectH,
    Impedance,
    FiniteConductivity,
    Symmetry,
    MasterSlave,

    // HFSS-specific
    Radiation,
    PML,

    // Q3D-specific
    ThinConductor,
    InfiniteGroundPlane,
    OpenBoundary,
}

/// Target assignment for boundaries, excitations, and mesh operations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Assignment {
    pub target_type: AssignmentTarget,
    /// References to geometry objects, faces, or edges.
    #[serde(default)]
    pub targets: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum AssignmentTarget {
    #[default]
    Object,
    Face,
    Edge,
}
