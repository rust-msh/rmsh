use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SimulationStatus {
    Idle,
    Solving,
    Finished,
    Failed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Material {
    pub name: String,
    pub relative_permittivity: f32,
    pub conductivity: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometryObject {
    pub id: u64,
    pub name: String,
    pub mesh_hint: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EmModel {
    pub name: String,
    pub objects: Vec<GeometryObject>,
    pub materials: Vec<Material>,
}

impl Default for EmModel {
    fn default() -> Self {
        Self {
            name: "Untitled Model".to_string(),
            objects: Vec::new(),
            materials: vec![Material {
                name: "Vacuum".to_string(),
                relative_permittivity: 1.0,
                conductivity: 0.0,
            }],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolveResult {
    pub field_preview: String,
    pub converged: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Project {
    pub id: String,
    pub title: String,
    pub model: EmModel,
    pub status: SimulationStatus,
    pub last_result: Option<SolveResult>,
}

impl Default for Project {
    fn default() -> Self {
        Self {
            id: "local-default".to_string(),
            title: "New Project".to_string(),
            model: EmModel::default(),
            status: SimulationStatus::Idle,
            last_result: None,
        }
    }
}
