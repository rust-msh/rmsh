use serde::{Deserialize, Serialize};

/// History-based parametric geometry container.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Geometry {
    /// Ordered list of modeling operations (the history).
    #[serde(default)]
    pub operations: Vec<GeometryOperation>,
    /// Current geometry state snapshot (regeneratable from operations).
    #[serde(default)]
    pub objects: Vec<GeoObject>,
}

/// A single step in the geometry operation history.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometryOperation {
    /// Step number (1-based, determines replay order).
    pub step: u32,
    /// The command performed.
    pub command: OperationCommand,
    /// Name of the object created or modified (if any).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub result_object: Option<String>,
    /// Command-specific parameters (parametric — may contain variable refs).
    #[serde(default)]
    pub parameters: serde_json::Value,
    /// Optional attributes set on the result object.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub attributes: Option<ObjectAttributes>,
}

/// Geometry modeling command type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OperationCommand {
    // Primitives
    CreateBox,
    CreateCylinder,
    CreateSphere,
    CreateCone,
    CreateTorus,
    CreatePolyline,
    CreateRectangle,
    CreateCircle,

    // Boolean operations
    Unite,
    Subtract,
    Intersect,

    // Transformations
    Move,
    Rotate,
    Mirror,
    Scale,
    DuplicateAlongLine,
    DuplicateAroundAxis,

    // Sweeps
    SweepAlongVector,
    SweepAlongPath,
    SweepAroundAxis,

    // Property modifications
    SetMaterial,
    SetColor,
    Rename,
    SetGroup,
    SetSolveInside,

    // Advanced
    Fillet,
    Chamfer,
    Section,
    Import,
}

/// Attributes that can be set on a geometry object.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ObjectAttributes {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub material: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub solve_inside: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub color: Option<[u8; 3]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transparency: Option<f32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coordinate_system: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group: Option<String>,
}

/// A geometry object in the current snapshot (derived from operation history).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeoObject {
    pub id: u64,
    pub name: String,
    /// Step number of the last operation that produced this object.
    pub derived_from_step: u32,
    #[serde(default = "default_material")]
    pub material: String,
    #[serde(default)]
    pub solve_inside: bool,
    #[serde(default = "default_color")]
    pub color: [u8; 3],
    #[serde(default)]
    pub transparency: f32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub group: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bounding_box: Option<BoundingBox>,
}

fn default_material() -> String {
    "vacuum".to_string()
}

fn default_color() -> [u8; 3] {
    [128, 128, 128]
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BoundingBox {
    pub min: [f64; 3],
    pub max: [f64; 3],
}
