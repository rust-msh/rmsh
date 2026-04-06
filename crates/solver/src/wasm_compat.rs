//! WASM compatibility layer for EMStudio solver.
//!
//! In WASM environments, we cannot use the full rem electromagnetic solver
//! (which includes C++ bindings and parallel libraries). Instead, we provide
//! stub implementations that:
//! 1. Accept the same API as Native solver
//! 2. Return simulated/synthetic results for preview
//! 3. Can optionally send solve requests to a backend worker
//!
//! This allows WASM builds to compile and run, with actual solving deferred
//! to Cloud services or native backends.

#[cfg(target_arch = "wasm32")]
pub mod wasm_stubs {
    use std::path::Path;
    use crate::progress::ProgressCallback;
    use crate::error::SolverError;
    use crate::solver_bridge::SolverResult;

    /// WASM stub: Pretend to solve by returning synthetic results.
    /// 
    /// This allows Local-First WASM mode to complete a solve cycle
    /// for demonstration. Real solves should be routed to:
    /// 1. Cloud backend (if professional/enterprise edition)
    /// 2. Native Desktop app via IPC
    pub fn dispatch_solver_stub(
        _config_path: &Path,
        progress: &dyn ProgressCallback,
    ) -> Result<SolverResult, SolverError> {
        if progress.is_cancelled() {
            return Err(SolverError::Cancelled);
        }

        // Simulate solve phases
        progress.on_phase("Loading configuration (WASM stub)");

        progress.on_phase("Assembling system (WASM stub)");
        progress.on_progress(0.3);

        progress.on_phase("Solving (WASM stub)");
        progress.on_progress(0.7);

        progress.on_phase("Post-processing (WASM stub)");
        progress.on_progress(0.95);

        Ok(SolverResult {
            converged: true,
            message: "Solve completed (WASM stub - use Desktop or Cloud for real results)"
                .to_string(),
        })
    }
}

/// Compiles to native solver bridge when target is NOT wasm32
#[cfg(not(target_arch = "wasm32"))]
pub use crate::solver_bridge::dispatch_solver as dispatch_solver;

/// Routes to appropriate solver based on target architecture
#[cfg(target_arch = "wasm32")]
pub use wasm_stubs::dispatch_solver_stub as dispatch_solver;

#[cfg(test)]
mod tests {
    #[test]
    #[cfg(target_arch = "wasm32")]
    fn test_wasm_stub_solver() {
        struct TestProgress;
        impl crate::progress::ProgressCallback for TestProgress {
            fn on_phase(&self, msg: &str) {
                println!("[Phase] {}", msg);
            }
            fn on_progress(&self, pct: f32) {
                println!("[Progress] {}%", (pct * 100.0) as u32);
            }
            fn is_cancelled(&self) -> bool {
                false
            }
        }

        let progress = TestProgress;
        let config_path = std::path::Path::new("/tmp/dummy.json");
        let result = dispatch_solver(config_path, &progress);

        assert!(result.is_ok());
        let res = result.unwrap();
        assert!(res.converged);
        assert!(res.message.contains("stub"));
    }
}
