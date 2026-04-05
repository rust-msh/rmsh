//! Solver error types.

use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SolverError {
    #[error("validation failed: {0}")]
    Validation(String),

    #[error("mesh generation failed: {0}")]
    MeshGeneration(String),

    #[error("solver execution failed: {0}")]
    SolverExecution(String),

    #[error("config generation failed: {0}")]
    ConfigGeneration(String),

    #[error("result extraction failed: {0}")]
    ResultExtraction(String),

    #[error("I/O error at {path}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("simulation cancelled by user")]
    Cancelled,
}

impl SolverError {
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        SolverError::Io {
            path: path.into(),
            source,
        }
    }
}
