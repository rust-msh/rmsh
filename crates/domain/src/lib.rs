pub mod worker_protocol;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn project_default_values() {
        let p = Project::default();
        assert_eq!(p.title, "New Project");
        assert_eq!(p.status, SimulationStatus::Idle);
        assert!(p.last_result.is_none());
        assert_eq!(p.model.materials.len(), 1);
        assert_eq!(p.model.materials[0].name, "Vacuum");
    }

    #[test]
    fn project_json_roundtrip() {
        let mut p = Project::default();
        p.title = "Test Project".into();
        p.model.objects.push(GeometryObject {
            id: 42,
            name: "Box1".into(),
            mesh_hint: "auto".into(),
        });
        p.last_result = Some(SolveResult {
            field_preview: "E-field preview".into(),
            converged: true,
        });

        let json = serde_json::to_string(&p).unwrap();
        let loaded: Project = serde_json::from_str(&json).unwrap();

        assert_eq!(loaded.title, "Test Project");
        assert_eq!(loaded.model.objects.len(), 1);
        assert_eq!(loaded.model.objects[0].name, "Box1");
        assert!(loaded.last_result.as_ref().unwrap().converged);
    }

    #[test]
    fn project_msgpack_roundtrip() {
        let mut p = Project::default();
        p.title = "MsgPack Test".into();
        p.model.materials.push(Material {
            name: "Copper".into(),
            relative_permittivity: 1.0,
            conductivity: 5.8e7,
        });

        let data = rmp_serde::to_vec(&p).unwrap();
        let loaded: Project = rmp_serde::from_slice(&data).unwrap();

        assert_eq!(loaded.title, "MsgPack Test");
        assert_eq!(loaded.model.materials.len(), 2);
        assert_eq!(loaded.model.materials[1].name, "Copper");
    }

    #[test]
    fn material_serialization() {
        let m = Material {
            name: "FR4".into(),
            relative_permittivity: 4.4,
            conductivity: 0.0,
        };
        let data = rmp_serde::to_vec(&m).unwrap();
        let loaded: Material = rmp_serde::from_slice(&data).unwrap();
        assert_eq!(loaded.name, "FR4");
        assert!((loaded.relative_permittivity - 4.4).abs() < 1e-6);
    }
}
