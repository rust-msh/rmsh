//! End-to-end solver pipeline: validate → mesh → adaptive solve → sweep → extract results.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use emstudio_domain::Design;
use emstudio_domain::result_store::{ConvergenceData, SParameterData, RlcgMatrixData, FarFieldData};
use emstudio_domain::solution_index::{SetupSolutionStatus, SolveStatus, SweepSolutionResult};
use emstudio_domain::solution_type::SolutionType;

use crate::adaptive;
use crate::error::SolverError;
use crate::mesh_bridge::{self, MeshStats};
use crate::progress::{ProgressCallback, SolverProgress};
use crate::results;
use crate::solver_bridge::SolverResult;
use crate::solver_log::SolverLog;
use crate::sweep;
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
    /// Convergence history.
    pub convergence: Option<ConvergenceData>,
    /// S-parameter data (from sweep).
    pub s_parameters: Vec<SParameterData>,
    /// RLCG matrix data (Q3D).
    pub rlcg_matrix: Option<RlcgMatrixData>,
    /// Far-field data (MoM/SBR).
    pub far_field: Option<FarFieldData>,
    /// Field data file path (.emsfld).
    pub field_data_path: Option<PathBuf>,
    /// Solution status for updating SolutionIndex.
    pub solution_status: SetupSolutionStatus,
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
        let log = SolverLog::new(&self.work_dir);

        // Ensure work directory exists.
        std::fs::create_dir_all(&self.work_dir)
            .map_err(|e| SolverError::io(&self.work_dir, e))?;

        let setup_name = design
            .analysis_setups
            .iter()
            .find(|s| s.enabled)
            .map(|s| s.name.clone())
            .unwrap_or_default();

        // ── Step 1: Validate ─────────────────────────────────────
        progress.on_progress(&SolverProgress::Validating);
        progress.on_phase("Validating design");
        log.info("Starting design validation");
        validate::validate_for_solve(design)?;
        log.info("Design validation passed");

        if progress.is_cancelled() {
            return Err(SolverError::Cancelled);
        }

        // ── Step 2: Rebuild geometry ─────────────────────────────
        progress.on_phase("Rebuilding geometry");
        log.info("Rebuilding geometry from operation history");
        let mut engine = emstudio_domain::geometry_engine::GeometryEngine::new();
        engine
            .rebuild(&design.geometry.operations, vars)
            .map_err(|e| SolverError::MeshGeneration(format!("Geometry rebuild failed: {e}")))?;
        log.info(format!("Geometry rebuilt: {} objects", engine.all_breps().len()));

        if progress.is_cancelled() {
            return Err(SolverError::Cancelled);
        }

        // ── Step 3: Extract surfaces ─────────────────────────────
        progress.on_progress(&SolverProgress::Meshing { percent: 0.0 });
        progress.on_phase("Extracting surface meshes");

        let surfaces = mesh_bridge::extract_brep_surfaces(
            engine.all_breps(),
            &design.geometry.objects,
        );
        let mesh_config = mesh_bridge::mesh_config_from_design(design);
        log.info(format!(
            "Extracted {} surface objects, max element size = {:.4}",
            surfaces.len(),
            mesh_config.max_element_size
        ));

        // ── Step 4: Adaptive solve loop ──────────────────────────
        log.info("Starting adaptive mesh refinement loop");
        let (convergence, converged) =
            adaptive::run_adaptive_loop(design, &surfaces, &mesh_config, &self.work_dir, progress)?;

        // Find the final mesh path (last pass).
        let final_pass = convergence.passes.len() as u32;
        let final_pass_dir = self.work_dir.join(format!("pass_{}", final_pass));
        let mesh_path = final_pass_dir.join("mesh.msh");
        let output_dir = final_pass_dir.join("output");

        let mesh_stats = MeshStats {
            num_nodes: convergence
                .passes
                .last()
                .map(|p| p.mesh.num_nodes as usize)
                .unwrap_or(0),
            num_tetrahedra: convergence
                .passes
                .last()
                .map(|p| p.mesh.num_tetrahedra as usize)
                .unwrap_or(0),
            num_triangles: 0,
        };

        log.info(format!(
            "Adaptive loop finished: {} passes, converged={}",
            final_pass, converged
        ));

        // ── Step 5: Write convergence data ───────────────────────
        let results_dir = self.work_dir.join("results");
        std::fs::create_dir_all(&results_dir)
            .map_err(|e| SolverError::io(&results_dir, e))?;

        results::write_convergence_json(&convergence, &results_dir)?;
        log.info("Convergence history written");

        // ── Step 6: Frequency sweeps ─────────────────────────────
        let s_parameters = if !design
            .analysis_setups
            .iter()
            .filter(|s| s.enabled)
            .flat_map(|s| &s.frequency_sweeps)
            .next()
            .is_none()
        {
            // No sweeps defined — try to extract S-params from single-frequency solve.
            if output_dir.join("port-S.csv").exists() {
                match results::extract_s_parameters(&output_dir, design, &setup_name) {
                    Ok(sp) => {
                        results::write_s_parameters_json(&sp, &results_dir)?;
                        log.info("S-parameters extracted from single-frequency solve");
                        vec![sp]
                    }
                    Err(_) => Vec::new(),
                }
            } else {
                Vec::new()
            }
        } else {
            log.info("Starting frequency sweeps");
            let sweep_results =
                sweep::run_frequency_sweeps(design, &mesh_path, &self.work_dir, progress)?;
            for (i, sp) in sweep_results.iter().enumerate() {
                let sweep_results_dir = results_dir.join(format!("sweep_{}", i));
                std::fs::create_dir_all(&sweep_results_dir)
                    .map_err(|e| SolverError::io(&sweep_results_dir, e))?;
                results::write_s_parameters_json(sp, &sweep_results_dir)?;
            }
            log.info(format!("{} frequency sweeps completed", sweep_results.len()));
            sweep_results
        };

        // ── Step 7: Extract additional results ───────────────────
        progress.on_progress(&SolverProgress::ExtractingResults);
        progress.on_phase("Extracting results");

        // RLCG matrix (Q3D electrostatic).
        let rlcg_matrix = if design.solution_type.is_q3d() {
            match results::extract_rlcg_matrix(&output_dir, design, &setup_name) {
                Ok(rlcg) => {
                    results::write_rlcg_json(&rlcg, &results_dir)?;
                    log.info("RLCG matrix extracted");
                    Some(rlcg)
                }
                Err(_) => None,
            }
        } else {
            None
        };

        // Far-field (MoM/SBR).
        let far_field = if matches!(
            design.solution_type,
            SolutionType::SBRPlus
        ) {
            match results::extract_far_field(&output_dir, design, &setup_name) {
                Ok(ff) => {
                    results::write_far_field_json(&ff, &results_dir)?;
                    log.info("Far-field data extracted");
                    Some(ff)
                }
                Err(_) => None,
            }
        } else {
            None
        };

        // Field data (.emsfld).
        let field_data_path =
            results::export_field_data(&output_dir, &results_dir, design)?;
        if field_data_path.is_some() {
            log.info("Field data exported to .emsfld");
        }

        // ── Step 8: Build solution status ────────────────────────
        let solution_status = SetupSolutionStatus {
            status: if converged {
                SolveStatus::Converged
            } else {
                SolveStatus::NotConverged
            },
            solved_at: Some(timestamp_iso8601()),
            converged_pass: if converged { Some(final_pass) } else { None },
            num_tetrahedra: Some(mesh_stats.num_tetrahedra as u64),
            num_triangles: None,
            final_delta_energy: convergence
                .passes
                .last()
                .and_then(|p| p.solution.max_delta_energy),
            is_stale: false,
            solved_variations: HashMap::new(),
            sweeps: s_parameters
                .iter()
                .map(|sp| {
                    (
                        sp.sweep.clone(),
                        SweepSolutionResult {
                            status: SolveStatus::Completed,
                            num_frequency_points: Some(sp.num_frequencies as u32),
                            result_path: format!("results/s_parameters.json"),
                        },
                    )
                })
                .collect(),
            rlcg_summary: None,
        };

        let solver_result = SolverResult {
            converged,
            message: if converged {
                format!("Converged after {} passes", final_pass)
            } else {
                format!("Did not converge after {} passes", final_pass)
            },
        };

        // ── Step 9: Flush log and report completion ──────────────
        log.info("Simulation pipeline completed");
        let _ = log.flush();

        progress.on_progress(&SolverProgress::Completed { converged });

        Ok(SolverOutput {
            mesh_path,
            mesh_stats,
            solver_result,
            output_dir,
            convergence: Some(convergence),
            s_parameters,
            rlcg_matrix,
            far_field,
            field_data_path,
            solution_status,
        })
    }

    /// Returns the working directory path.
    pub fn work_dir(&self) -> &Path {
        &self.work_dir
    }
}

fn timestamp_iso8601() -> String {
    use std::time::SystemTime;
    let secs = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    format!("{}Z", secs)
}
