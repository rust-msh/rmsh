use std::collections::HashMap;
#[cfg(not(target_arch = "wasm32"))]
use std::path::Path;

use emstudio_domain::{Project, SolveResult};
use emstudio_solver::{PlaceholderSolver, Solver};
use thiserror::Error;

#[cfg(target_arch = "wasm32")]
mod wasm_backend;
#[cfg(target_arch = "wasm32")]
pub use wasm_backend::WasmBackend;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunMode {
    Standalone,
    Cloud,
    /// Web Local-First mode: OPFS storage + Web Worker solver.
    LocalFirst,
}

impl Default for RunMode {
    fn default() -> Self {
        Self::Standalone
    }
}

#[derive(Debug, Error)]
pub enum BackendError {
    #[error("project not found: {0}")]
    ProjectNotFound(String),
    #[error("IO error: {0}")]
    IoError(String),
    #[error("serialize error: {0}")]
    SerializeError(String),
    #[error("deserialize error: {0}")]
    DeserializeError(String),
}

// ---------------------------------------------------------------------------
// File I/O (MessagePack .emsp format) — native only
// ---------------------------------------------------------------------------

#[cfg(not(target_arch = "wasm32"))]
pub fn save_project_to_file(project: &Project, path: &Path) -> Result<(), BackendError> {
    let data =
        rmp_serde::to_vec(project).map_err(|e| BackendError::SerializeError(e.to_string()))?;
    std::fs::write(path, data).map_err(|e| BackendError::IoError(e.to_string()))?;
    Ok(())
}

#[cfg(not(target_arch = "wasm32"))]
pub fn load_project_from_file(path: &Path) -> Result<Project, BackendError> {
    let data = std::fs::read(path).map_err(|e| BackendError::IoError(e.to_string()))?;
    let project: Project =
        rmp_serde::from_slice(&data).map_err(|e| BackendError::DeserializeError(e.to_string()))?;
    Ok(project)
}

// ---------------------------------------------------------------------------
// Backend trait
// ---------------------------------------------------------------------------

pub trait Backend {
    fn save_project(&mut self, project: Project) -> Result<(), BackendError>;
    fn load_project(&self, id: &str) -> Result<Project, BackendError>;
    fn solve(&self, project: &Project) -> Result<SolveResult, BackendError>;
    fn mode(&self) -> RunMode;

    /// Poll for async results (WASM: drains worker responses). Default no-op.
    fn poll(&mut self) {}

    /// Take a pending solve result, if one has arrived. Default returns `None`.
    fn take_solve_result(&mut self) -> Option<SolveResult> {
        None
    }
}

pub struct StandaloneBackend {
    projects: HashMap<String, Project>,
    solver: PlaceholderSolver,
}

impl Default for StandaloneBackend {
    fn default() -> Self {
        Self {
            projects: HashMap::new(),
            solver: PlaceholderSolver,
        }
    }
}

impl Backend for StandaloneBackend {
    fn save_project(&mut self, project: Project) -> Result<(), BackendError> {
        self.projects.insert(project.id.clone(), project);
        Ok(())
    }

    fn load_project(&self, id: &str) -> Result<Project, BackendError> {
        self.projects
            .get(id)
            .cloned()
            .ok_or_else(|| BackendError::ProjectNotFound(id.to_string()))
    }

    fn solve(&self, project: &Project) -> Result<SolveResult, BackendError> {
        Ok(self.solver.solve(&project.model))
    }

    fn mode(&self) -> RunMode {
        RunMode::Standalone
    }
}

pub struct CloudBackend {
    endpoint: String,
}

impl CloudBackend {
    pub fn new(endpoint: impl Into<String>) -> Self {
        Self {
            endpoint: endpoint.into(),
        }
    }

    pub fn endpoint(&self) -> &str {
        &self.endpoint
    }
}

impl Backend for CloudBackend {
    fn save_project(&mut self, _project: Project) -> Result<(), BackendError> {
        Ok(())
    }

    fn load_project(&self, id: &str) -> Result<Project, BackendError> {
        Ok(Project {
            id: id.to_string(),
            title: format!("Cloud Project @ {}", self.endpoint),
            ..Project::default()
        })
    }

    fn solve(&self, project: &Project) -> Result<SolveResult, BackendError> {
        Ok(SolveResult {
            field_preview: format!(
                "Cloud placeholder solve for '{}' via {}",
                project.title, self.endpoint
            ),
            converged: true,
        })
    }

    fn mode(&self) -> RunMode {
        RunMode::Cloud
    }
}

pub fn default_backend(mode: RunMode) -> Box<dyn Backend> {
    match mode {
        RunMode::Standalone => Box::<StandaloneBackend>::default(),
        RunMode::Cloud => Box::new(CloudBackend::new("https://api.example.local")),
        #[cfg(target_arch = "wasm32")]
        RunMode::LocalFirst => {
            Box::new(WasmBackend::new("./worker_bootstrap.js").expect("failed to spawn worker"))
        }
        #[cfg(not(target_arch = "wasm32"))]
        RunMode::LocalFirst => {
            // LocalFirst mode is only available on WASM; fall back to Standalone on native.
            Box::<StandaloneBackend>::default()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use emstudio_domain::{GeometryObject, Material, SolveResult};
    use std::fs;

    fn sample_project() -> Project {
        let mut p = Project::default();
        p.title = "File IO Test".into();
        p.model.objects.push(GeometryObject {
            id: 1,
            name: "Cylinder1".into(),
            mesh_hint: "fine".into(),
        });
        p.model.materials.push(Material {
            name: "Gold".into(),
            relative_permittivity: 1.0,
            conductivity: 4.1e7,
        });
        p.last_result = Some(SolveResult {
            field_preview: "test field data".into(),
            converged: true,
        });
        p
    }

    #[test]
    fn save_and_load_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.emsp");

        let original = sample_project();
        save_project_to_file(&original, &path).unwrap();

        assert!(path.exists());
        assert!(fs::metadata(&path).unwrap().len() > 0);

        let loaded = load_project_from_file(&path).unwrap();
        assert_eq!(loaded.title, "File IO Test");
        assert_eq!(loaded.model.objects.len(), 1);
        assert_eq!(loaded.model.objects[0].name, "Cylinder1");
        assert_eq!(loaded.model.materials.len(), 2);
        assert_eq!(loaded.model.materials[1].name, "Gold");
        assert!(loaded.last_result.as_ref().unwrap().converged);
    }

    #[test]
    fn save_overwrites_existing_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("overwrite.emsp");

        let mut p = Project::default();
        p.title = "Version 1".into();
        save_project_to_file(&p, &path).unwrap();

        p.title = "Version 2".into();
        save_project_to_file(&p, &path).unwrap();

        let loaded = load_project_from_file(&path).unwrap();
        assert_eq!(loaded.title, "Version 2");
    }

    #[test]
    fn load_nonexistent_file_returns_error() {
        let result = load_project_from_file(Path::new("/nonexistent/path/file.emsp"));
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(err, BackendError::IoError(_)));
    }

    #[test]
    fn load_corrupt_file_returns_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("corrupt.emsp");
        fs::write(&path, b"this is not valid msgpack data").unwrap();

        let result = load_project_from_file(&path);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), BackendError::DeserializeError(_)));
    }

    #[test]
    fn standalone_backend_in_memory() {
        let mut backend = StandaloneBackend::default();
        let p = sample_project();

        backend.save_project(p.clone()).unwrap();
        let loaded = backend.load_project(&p.id).unwrap();
        assert_eq!(loaded.title, p.title);
    }

    #[test]
    fn standalone_backend_not_found() {
        let backend = StandaloneBackend::default();
        let result = backend.load_project("nonexistent");
        assert!(matches!(result.unwrap_err(), BackendError::ProjectNotFound(_)));
    }

    #[test]
    fn standalone_backend_solve() {
        let backend = StandaloneBackend::default();
        let p = Project::default();
        let result = backend.solve(&p).unwrap();
        // PlaceholderSolver returns something
        assert!(!result.field_preview.is_empty());
    }
}
