//! Simulation progress reporting.

/// Progress updates from the solver pipeline.
#[derive(Debug, Clone)]
pub enum SolverProgress {
    /// Pre-solve model validation.
    Validating,
    /// Mesh generation in progress.
    Meshing { percent: f32 },
    /// Solver running an adaptive pass.
    Solving {
        pass: u32,
        max_passes: u32,
        delta_s: f64,
    },
    /// Frequency sweep in progress.
    Sweeping {
        freq_idx: usize,
        total_freqs: usize,
        freq_hz: f64,
    },
    /// Extracting and writing results.
    ExtractingResults,
    /// Simulation completed.
    Completed { converged: bool },
    /// Simulation failed.
    Failed { message: String },
}

/// Callback trait for receiving solver progress updates.
pub trait ProgressCallback: Send {
    /// Called with progress updates during the simulation.
    fn on_progress(&self, progress: &SolverProgress);

    /// Called when a named phase begins (e.g., "Mesh generation", "Solving pass 3").
    fn on_phase(&self, phase: &str);

    /// Called with log messages.
    fn on_log(&self, message: &str);

    /// Returns true if the user has requested cancellation.
    fn is_cancelled(&self) -> bool;
}

/// No-op progress callback for batch/test usage.
pub struct NoOpProgress;

impl ProgressCallback for NoOpProgress {
    fn on_progress(&self, _progress: &SolverProgress) {}
    fn on_phase(&self, _phase: &str) {}
    fn on_log(&self, _message: &str) {}
    fn is_cancelled(&self) -> bool {
        false
    }
}
