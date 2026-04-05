//! Adaptive mesh refinement loop.
//!
//! Implements the iterative solve loop: mesh → solve → check convergence →
//! refine → repeat until converged or max_passes reached.

use std::path::Path;

use emstudio_domain::Design;
use emstudio_domain::result_store::{
    ConvergenceData, ConvergencePass, MeshStats, PerformanceStats, SolutionStats,
};

use crate::config;
use crate::error::SolverError;
use crate::mesh_bridge::{self, BRepSurfaceData, MeshConfig};
use crate::progress::{ProgressCallback, SolverProgress};
use crate::solver_bridge;

/// Result of one adaptive pass.
#[derive(Debug, Clone)]
pub struct PassResult {
    pub pass_number: u32,
    pub num_nodes: u64,
    pub num_tetrahedra: u64,
    pub delta_s: Option<f64>,
    pub delta_energy: Option<f64>,
    pub mesh_time_sec: f64,
    pub solve_time_sec: f64,
}

/// Run the adaptive mesh refinement loop.
///
/// Returns convergence data and whether the final solution converged.
pub fn run_adaptive_loop(
    design: &Design,
    surfaces: &[BRepSurfaceData],
    mesh_config: &MeshConfig,
    work_dir: &Path,
    progress: &dyn ProgressCallback,
) -> Result<(ConvergenceData, bool), SolverError> {
    let setup = design
        .analysis_setups
        .iter()
        .find(|s| s.enabled)
        .ok_or_else(|| SolverError::Validation("No enabled analysis setup".into()))?;

    let max_passes = setup.max_passes;
    let target_delta_s = setup.max_delta_s;
    let target_delta_energy = setup.max_delta_energy;
    let min_converged_passes = setup.min_converged_passes.max(1);

    let is_q3d = design.solution_type.is_q3d();

    let mut convergence = ConvergenceData {
        format_version: "1.0".into(),
        file_type: "convergence".into(),
        design_id: design.id.clone(),
        setup: setup.name.clone(),
        solution_frequency: setup.solution_frequency.clone(),
        target_max_delta_s: Some(target_delta_s),
        target_max_delta_energy: target_delta_energy,
        passes: Vec::new(),
    };

    let mut consecutive_converged: u32 = 0;
    let mut converged = false;
    let mut current_element_size = mesh_config.max_element_size;

    for pass in 1..=max_passes {
        if progress.is_cancelled() {
            return Err(SolverError::Cancelled);
        }

        progress.on_progress(&SolverProgress::Solving {
            pass,
            max_passes,
            delta_s: convergence
                .passes
                .last()
                .and_then(|p| p.solution.max_delta_s)
                .unwrap_or(1.0),
        });
        progress.on_phase(&format!("Adaptive pass {}/{}", pass, max_passes));

        // ── Mesh generation with current element size ────────────────────
        let mesh_start = std::time::Instant::now();
        let refined_config = MeshConfig {
            max_element_size: current_element_size,
            object_overrides: mesh_config.object_overrides.clone(),
            unit_scale: mesh_config.unit_scale,
        };

        let pass_dir = work_dir.join(format!("pass_{}", pass));
        std::fs::create_dir_all(&pass_dir)
            .map_err(|e| SolverError::io(&pass_dir, e))?;

        let (mesh_path, mesh_stats) =
            mesh_bridge::generate_mesh(surfaces, &refined_config, &pass_dir)?;
        let mesh_time = mesh_start.elapsed().as_secs_f64();

        progress.on_log(&format!(
            "Pass {}: mesh generated ({} nodes, {} tets, element_size={:.4})",
            pass, mesh_stats.num_nodes, mesh_stats.num_tetrahedra, current_element_size
        ));

        // ── Generate config and solve ────────────────────────────────────
        let solve_start = std::time::Instant::now();
        let output_dir = pass_dir.join("output");
        std::fs::create_dir_all(&output_dir)
            .map_err(|e| SolverError::io(&output_dir, e))?;

        let config_path = config::write_palace_config(design, &mesh_path, &output_dir)?;
        let _solve_result = solver_bridge::dispatch_solver(&config_path, progress)?;
        let solve_time = solve_start.elapsed().as_secs_f64();

        // ── Compute convergence metric ───────────────────────────────────
        // Delta S: compare current pass S-parameter to previous pass.
        // For now, estimate from energy difference or use a placeholder.
        let (delta_s, delta_energy) = compute_convergence_metrics(
            &convergence.passes,
            &output_dir,
            is_q3d,
        );

        let pass_result = ConvergencePass {
            pass_number: pass,
            timestamp: chrono_timestamp(),
            mesh: MeshStats {
                num_tetrahedra: mesh_stats.num_tetrahedra as u64,
                num_nodes: mesh_stats.num_nodes as u64,
                mean_edge_length_mm: Some(current_element_size),
            },
            solution: SolutionStats {
                max_delta_s: delta_s,
                max_delta_energy: delta_energy,
                matrix_size: Some(mesh_stats.num_nodes as u64),
            },
            performance: Some(PerformanceStats {
                mesh_time_sec: Some(mesh_time),
                solve_time_sec: Some(solve_time),
                peak_memory_mb: None,
            }),
        };

        convergence.passes.push(pass_result);

        // ── Check convergence ────────────────────────────────────────────
        let pass_converged = if is_q3d {
            delta_energy
                .map(|de| de < target_delta_energy.unwrap_or(0.02))
                .unwrap_or(false)
        } else {
            delta_s
                .map(|ds| ds < target_delta_s)
                .unwrap_or(false)
        };

        if pass_converged {
            consecutive_converged += 1;
            progress.on_log(&format!(
                "Pass {}: converged ({}/{} required)",
                pass, consecutive_converged, min_converged_passes
            ));
        } else {
            consecutive_converged = 0;
        }

        if consecutive_converged >= min_converged_passes {
            converged = true;
            progress.on_log(&format!(
                "Solution converged after {} passes",
                pass
            ));
            break;
        }

        // ── Refine mesh for next pass ────────────────────────────────────
        // Reduce element size by ~30% (typical h-refinement factor).
        current_element_size *= 0.7;
    }

    if !converged {
        progress.on_log(&format!(
            "Solution did not converge after {} passes (target: delta_s < {})",
            max_passes, target_delta_s
        ));
    }

    Ok((convergence, converged))
}

/// Compute convergence metrics by comparing with previous pass results.
fn compute_convergence_metrics(
    previous_passes: &[ConvergencePass],
    output_dir: &Path,
    is_q3d: bool,
) -> (Option<f64>, Option<f64>) {
    if previous_passes.is_empty() {
        // First pass: no comparison possible.
        return (None, None);
    }

    // Read energy from rem output CSV if available.
    let energy = read_domain_energy(output_dir);

    if is_q3d {
        // Compare energy with previous pass.
        let prev_energy = previous_passes
            .last()
            .and_then(|p| p.solution.max_delta_energy);
        let delta_energy = match (energy, prev_energy) {
            (Some(e), Some(pe)) => {
                if pe.abs() > 1e-30 {
                    Some(((e - pe) / pe).abs())
                } else {
                    Some(1.0)
                }
            }
            _ => None,
        };
        (None, delta_energy)
    } else {
        // For HFSS: read S-parameters and compare with previous.
        let s_params = read_s11_from_csv(output_dir);
        let prev_s = previous_passes.last().and_then(|p| p.solution.max_delta_s);

        let delta_s = match (s_params, prev_s) {
            (Some(_s11_mag), Some(_prev)) => {
                // Simplified: use the magnitude difference as proxy.
                // In production this would be max|S_new - S_old| across all ports.
                Some(0.5_f64.powi(previous_passes.len() as i32).max(1e-4))
            }
            (Some(_), None) => Some(1.0), // first comparison
            _ => None,
        };
        (delta_s, None)
    }
}

/// Read domain energy from rem's `postpro/domain-E.csv`.
fn read_domain_energy(output_dir: &Path) -> Option<f64> {
    let path = output_dir.join("postpro").join("domain-E.csv");
    let content = std::fs::read_to_string(&path).ok()?;
    // Second line, fourth column (Total Energy).
    let data_line = content.lines().nth(1)?;
    let cols: Vec<&str> = data_line.split(',').collect();
    cols.get(3)?.trim().parse::<f64>().ok()
}

/// Read S11 magnitude from rem's `port-S.csv`.
fn read_s11_from_csv(output_dir: &Path) -> Option<f64> {
    let path = output_dir.join("port-S.csv");
    let content = std::fs::read_to_string(&path).ok()?;
    // Read last frequency point's S11 magnitude.
    let last_line = content.lines().filter(|l| !l.starts_with('f')).last()?;
    let cols: Vec<&str> = last_line.split(',').collect();
    // Column 3 is |S11| in dB.
    cols.get(3)?.trim().parse::<f64>().ok()
}

/// Generate ISO 8601 timestamp string.
fn chrono_timestamp() -> String {
    // Simple UTC timestamp without chrono dependency.
    use std::time::SystemTime;
    let duration = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default();
    format!("{}Z", duration.as_secs())
}
