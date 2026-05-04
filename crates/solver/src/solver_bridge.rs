//! Dispatch rem electromagnetic solvers based on problem type.

use std::path::Path;

use rem_config::{load_config, ProblemType};
use rem_parallel::NoComm;

use crate::error::SolverError;
use crate::progress::ProgressCallback;

/// Result returned by a solver run.
#[derive(Debug, Clone)]
pub struct SolverResult {
    pub converged: bool,
    pub message: String,
}

/// Load the Palace config from `config_path` and dispatch to the appropriate
/// rem solver module.
pub fn dispatch_solver(
    config_path: &Path,
    progress: &dyn ProgressCallback,
) -> Result<SolverResult, SolverError> {
    if progress.is_cancelled() {
        return Err(SolverError::Cancelled);
    }

    progress.on_phase("Loading solver configuration");

    let config = load_config(config_path)
        .map_err(|e| SolverError::SolverExecution(format!("Failed to load config: {e}")))?;

    progress.on_phase(&format!(
        "Starting {:?} solver",
        config.problem.problem_type
    ));

    let result = match config.problem.problem_type {
        ProblemType::Electrostatic => {
            let comm = NoComm;
            rem_electrostatic::run(&config, &comm)
                .map(|()| SolverResult {
                    converged: true,
                    message: "Electrostatic solve completed".into(),
                })
        }
        ProblemType::Magnetostatic => {
            let comm = NoComm;
            rem_magnetostatic::run(&config, &comm)
                .map(|()| SolverResult {
                    converged: true,
                    message: "Magnetostatic solve completed".into(),
                })
        }
        ProblemType::Eigenmode => {
            let comm = NoComm;
            rem_eigenmode::run(&config, &comm)
                .map(|()| SolverResult {
                    converged: true,
                    message: "Eigenmode solve completed".into(),
                })
        }
        ProblemType::Driven => {
            let comm = NoComm;
            rem_driven::run(&config, &comm)
                .map(|()| SolverResult {
                    converged: true,
                    message: "Driven frequency solve completed".into(),
                })
        }
        ProblemType::MoM => {
            rem_mom::run(&config)
                .map(|()| SolverResult {
                    converged: true,
                    message: "MoM solve completed".into(),
                })
        }
        ProblemType::SBR => {
            rem_sbr::run(&config)
                .map(|()| SolverResult {
                    converged: true,
                    message: "SBR+ solve completed".into(),
                })
        }
        ProblemType::BEM => {
            rem_bem::run(&config)
                .map(|()| SolverResult {
                    converged: true,
                    message: "BEM solve completed".into(),
                })
        }
        ProblemType::Planar => {
            rem_planar::run(&config)
                .map(|()| SolverResult {
                    converged: true,
                    message: "Planar MoM solve completed".into(),
                })
        }
        ProblemType::Transient => {
            Err(rem_core::RemError::NotImplemented(
                "Transient solver not yet available".into(),
            ))
        }
        ProblemType::FEBI => {
            Err(rem_core::RemError::NotImplemented(
                "FEBI solver not yet available via emstudio bridge".into(),
            ))
        }
    };

    result.map_err(|e: rem_core::RemError| {
        let msg = e.to_string();
        progress.on_log(&format!("Solver error: {msg}"));
        SolverError::SolverExecution(msg)
    })
}
