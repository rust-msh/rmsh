// ---------------------------------------------------------------------------
// New domain model modules (Milestone 3)
// ---------------------------------------------------------------------------
pub mod solution_type;
pub mod variable;
pub mod material;
pub mod geometry;
pub mod geometry_engine;
pub mod boundary;
pub mod excitation;
pub mod net;
pub mod mesh;
pub mod analysis;
pub mod radiation;
pub mod output_variable;
pub mod field_overlay;
pub mod optimetrics;
pub mod report;
pub mod result_store;
pub mod quantity_expr;
pub mod solution_index;
pub mod design;
pub mod project;
pub mod expression;
pub mod validation;
pub mod dependency;
pub mod file_io;
pub mod worker_protocol;

// Re-export top-level types for convenience
pub use project::EmProject;
pub use design::Design;
pub use solution_type::SolutionType;

// ---------------------------------------------------------------------------
// Legacy types (kept for backward compatibility with app/infra/solver/worker)
// ---------------------------------------------------------------------------
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

    // -- Legacy type tests (kept for backward compat) --

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

    // -- New EmProject tests --

    #[test]
    fn em_project_default_json_roundtrip() {
        let proj = EmProject::default();
        let json = proj.to_json_string().unwrap();
        let loaded = EmProject::from_json_str(&json).unwrap();
        assert_eq!(loaded.metadata.version, "1.0.0");
        assert_eq!(loaded.metadata.application, "EMStudio");
        assert!(loaded.designs.is_empty());
    }

    #[test]
    fn em_project_with_design_roundtrip() {
        use crate::design::Design;
        use crate::material::{MaterialDef, MaterialCategory};
        use crate::solution_type::SolutionType;
        use std::collections::HashMap;

        let mut proj = EmProject::default();
        let mut design = Design {
            id: "design-001".into(),
            name: "Patch Antenna".into(),
            solution_type: SolutionType::DrivenModal,
            units: "mm".into(),
            design_settings: Default::default(),
            local_variables: HashMap::new(),
            definitions: Default::default(),
            geometry: Default::default(),
            boundaries: Vec::new(),
            excitations: Vec::new(),
            nets: Vec::new(),
            mesh_operations: Vec::new(),
            analysis_setups: Vec::new(),
            radiation: Default::default(),
            output_variables: Vec::new(),
            field_overlays: Vec::new(),
            optimetrics: Vec::new(),
            reports: Vec::new(),
            solution_index: Default::default(),
        };

        design.definitions.materials.push(MaterialDef {
            name: "copper".into(),
            category: MaterialCategory::Conductor,
            properties: Default::default(),
            appearance: None,
        });

        proj.designs.push(design);

        let json = proj.to_json_string().unwrap();
        let loaded = EmProject::from_json_str(&json).unwrap();

        assert_eq!(loaded.designs.len(), 1);
        assert_eq!(loaded.designs[0].name, "Patch Antenna");
        assert_eq!(loaded.designs[0].solution_type, SolutionType::DrivenModal);
        assert_eq!(loaded.designs[0].definitions.materials[0].name, "copper");
    }

    #[test]
    fn em_project_validation_catches_dangling_material() {
        use crate::design::Design;
        use crate::geometry::{GeoObject, Geometry};
        use crate::solution_type::SolutionType;
        use std::collections::HashMap;

        let mut proj = EmProject::default();
        let design = Design {
            id: "d1".into(),
            name: "Test".into(),
            solution_type: SolutionType::DrivenModal,
            units: "mm".into(),
            design_settings: Default::default(),
            local_variables: HashMap::new(),
            definitions: Default::default(),
            geometry: Geometry {
                operations: Vec::new(),
                objects: vec![GeoObject {
                    id: 1,
                    name: "Box1".into(),
                    derived_from_step: 1,
                    material: "nonexistent_material".into(),
                    solve_inside: false,
                    color: [128, 128, 128],
                    transparency: 0.0,
                    group: None,
                    bounding_box: None,
                }],
            },
            boundaries: Vec::new(),
            excitations: Vec::new(),
            nets: Vec::new(),
            mesh_operations: Vec::new(),
            analysis_setups: Vec::new(),
            radiation: Default::default(),
            output_variables: Vec::new(),
            field_overlays: Vec::new(),
            optimetrics: Vec::new(),
            reports: Vec::new(),
            solution_index: Default::default(),
        };
        proj.designs.push(design);

        let errors = proj.validate();
        assert!(!errors.is_empty());
        assert!(errors[0].contains("nonexistent_material"));
    }
}
