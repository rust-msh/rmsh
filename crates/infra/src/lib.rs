use std::collections::HashMap;

use emstudio_domain::{Project, SolveResult};
use emstudio_solver::{PlaceholderSolver, Solver};
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunMode {
    Standalone,
    Cloud,
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
}

pub trait Backend {
    fn save_project(&mut self, project: Project) -> Result<(), BackendError>;
    fn load_project(&self, id: &str) -> Result<Project, BackendError>;
    fn solve(&self, project: &Project) -> Result<SolveResult, BackendError>;
    fn mode(&self) -> RunMode;
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
    }
}
