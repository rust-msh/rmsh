//! End-to-end solver pipeline: validate → mesh → solve → collect results.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use emstudio_domain::Design;

use crate::config;
use crate::error::SolverError;
use crate::mesh_bridge::{self, MeshStats};
use crate::progress::{ProgressCallback, SolverProgress};
use crate::solver_bridge::{self, SolverResult};
use crate::validate;

/// Output produced by a complete solver pipeline run.
#[derive(Debug)]
pub struct SolverOutput {
    /// Path to the generated .msh mesh file.
    pub mesh_path: PathBuf,
    /// Mesh statistics.
    pub mesh_stats: MeshStats,
    /// Solver result (convergence, message).
    pub solver_result: SolverResult,
    /// Path to the rem output directory (contains CSV, VTK, etc.).
    pub output_dir: PathBuf,
}

/// Orchestrates the full simulation pipeline.
pub struct SolverPipeline {
    /// Working directory for mesh files, config, and solver output.
    work_dir: PathBuf,
}

impl SolverPipeline {
    /// Create a pipeline using the given working directory.
    pub fn new(work_dir: impl Into<PathBuf>) -> Self {
        SolverPipeline {
            work_dir: work_dir.into(),
        }
    }

    /// Create a pipeline using a temporary directory.
    pub fn with_temp_dir() -> Result<Self, SolverError> {
        let dir = tempfile::tempdir()
            .map_err(|e| SolverError::io("tempdir", e))?;
        #[allow(deprecated)]
        let path = dir.into_path();
        Ok(SolverPipeline { work_dir: path })
    }

    /// Run the complete simulation pipeline.
    ///
    /// `vars` provides resolved variable values for parametric geometry.
    pub fn run(
        &self,
        design: &Design,
        vars: &HashMap<String, f64>,
        progress: &dyn ProgressCallback,
    ) -> Result<SolverOutput, SolverError> {
        // Ensure work directory exists.
        std::fs::create_dir_all(&self.work_dir)
            .map_err(|e| SolverError::io(&self.work_dir, e))?;

        // ── Step 1: Validate ─────────────────────────────────────
        progress.on_progress(&SolverProgress::Validating);
        progress.on_phase("Validating design");
        validate::validate_for_solve(design)?;

        if progress.is_cancelled() {
            return Err(SolverError::Cancelled);
        }

        // ── Step 2: Rebuild geometry ─────────────────────────────
        progress.on_phase("Rebuilding geometry");
        let mut engine = emstudio_domain::geometry_engine::GeometryEngine::new();
        engine
            .rebuild(&design.geometry.operations, vars)
            .map_err(|e| SolverError::MeshGeneration(format!("Geometry rebuild failed: {e}")))?;

        if progress.is_cancelled() {
            return Err(SolverError::Cancelled);
        }

        // ── Step 3: Extract surfaces and generate mesh ───────────
        progress.on_progress(&SolverProgress::Meshing { percent: 0.0 });
        progress.on_phase("Generating mesh");

        let surfaces = mesh_bridge::extract_brep_surfaces(
            engine.all_breps(),
            &design.geometry.objects,
        );
        let mesh_config = mesh_bridge::mesh_config_from_design(design);

        progress.on_progress(&SolverProgress::Meshing { percent: 30.0 });

        let (mesh_path, mesh_stats) =
            mesh_bridge::generate_mesh(&surfaces, &mesh_config, &self.work_dir)?;

        progress.on_progress(&SolverProgress::Meshing { percent: 100.0 });
        progress.on_log(&format!(
            "Mesh generated: {} nodes, {} tetrahedra",
            mesh_stats.num_nodes, mesh_stats.num_tetrahedra
        ));

        if progress.is_cancelled() {
            return Err(SolverError::Cancelled);
        }

        // ── Step 4: Generate rem config ──────────────────────────
        progress.on_phase("Generating solver configuration");
        let output_dir = self.work_dir.join("output");
        std::fs::create_dir_all(&output_dir)
            .map_err(|e| SolverError::io(&output_dir, e))?;

        let config_path = config::write_palace_config(design, &mesh_path, &output_dir)?;

        if progress.is_cancelled() {
            return Err(SolverError::Cancelled);
        }

        // ── Step 5: Run solver ───────────────────────────────────
        progress.on_progress(&SolverProgress::Solving {
            pass: 1,
            max_passes: design
                .analysis_setups
                .iter()
                .find(|s| s.enabled)
                .map(|s| s.max_passes)
                .unwrap_or(1),
            delta_s: 1.0,
        });

        let solver_result = solver_bridge::dispatch_solver(&config_path, progress)?;

        // ── Step 6: Report completion ────────────────────────────
        progress.on_progress(&SolverProgress::Completed {
            converged: solver_result.converged,
        });

        Ok(SolverOutput {
            mesh_path,
            mesh_stats,
            solver_result,
            output_dir,
        })
    }

    /// Returns the working directory path.
    pub fn work_dir(&self) -> &Path {
        &self.work_dir
    }
}
