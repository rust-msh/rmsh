//! EMStudio solver crate — orchestrates rmsh mesh generation and rem EM solvers.

pub mod error;
pub mod config;
pub mod mesh_bridge;
pub mod solver_bridge;
pub mod pipeline;
pub mod progress;
pub mod validate;
pub mod adaptive;
pub mod sweep;
pub mod results;
pub mod solver_log;
pub mod wasm_compat;

use emstudio_domain::{EmModel, SolveResult};

pub use error::SolverError;
pub use pipeline::{SolverPipeline, SolverOutput};
pub use progress::{ProgressCallback, SolverProgress, NoOpProgress};

// Route to platform-specific solver implementation
pub use wasm_compat::dispatch_solver;

// ---------------------------------------------------------------------------
// Legacy interface (kept for backward compatibility with worker crate)
// ---------------------------------------------------------------------------

pub trait Solver {
    fn solve(&self, model: &EmModel) -> SolveResult;
}

#[derive(Default)]
pub struct PlaceholderSolver;

impl Solver for PlaceholderSolver {
    fn solve(&self, model: &EmModel) -> SolveResult {
        SolveResult {
            field_preview: format!(
                "Placeholder result for model '{}' with {} objects",
                model.name,
                model.objects.len()
            ),
            converged: true,
        }
    }
}
