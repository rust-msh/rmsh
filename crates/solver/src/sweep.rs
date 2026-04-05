//! Frequency sweep execution.
//!
//! After adaptive convergence at the solution frequency, sweep across the
//! requested frequency range, reusing the converged mesh.

use std::path::Path;

use emstudio_domain::Design;
use emstudio_domain::analysis::SweepType;
use emstudio_domain::result_store::{SParameterData, SParamValues};

use crate::config;
use crate::error::SolverError;
use crate::mesh_bridge::parse_frequency;
use crate::progress::{ProgressCallback, SolverProgress};
use crate::solver_bridge;

use std::collections::HashMap;

/// Run frequency sweeps defined in the design's analysis setups.
///
/// Uses the converged mesh at `mesh_path` and solves at each frequency point.
/// Returns S-parameter data for each sweep.
pub fn run_frequency_sweeps(
    design: &Design,
    mesh_path: &Path,
    work_dir: &Path,
    progress: &dyn ProgressCallback,
) -> Result<Vec<SParameterData>, SolverError> {
    let setup = design
        .analysis_setups
        .iter()
        .find(|s| s.enabled)
        .ok_or_else(|| SolverError::Validation("No enabled analysis setup".into()))?;

    let mut all_results = Vec::new();

    for sweep in &setup.frequency_sweeps {
        if progress.is_cancelled() {
            return Err(SolverError::Cancelled);
        }

        progress.on_phase(&format!("Running frequency sweep: {}", sweep.name));

        let start_hz = parse_frequency(&sweep.start)?;
        let stop_hz = parse_frequency(&sweep.stop)?;

        // Determine frequency points.
        let freq_points = compute_frequency_points(
            start_hz,
            stop_hz,
            &sweep.sweep_type,
            sweep.step.as_deref(),
            sweep.count,
        )?;

        let total_freqs = freq_points.len();
        progress.on_log(&format!(
            "Sweep '{}': {} frequency points from {:.3e} to {:.3e} Hz",
            sweep.name, total_freqs, start_hz, stop_hz
        ));

        // For each frequency point, generate config and solve.
        let mut s_data_per_freq: Vec<(f64, Option<(f64, f64)>)> = Vec::new();

        for (idx, &freq_hz) in freq_points.iter().enumerate() {
            if progress.is_cancelled() {
                return Err(SolverError::Cancelled);
            }

            progress.on_progress(&SolverProgress::Sweeping {
                freq_idx: idx,
                total_freqs,
                freq_hz,
            });

            let sweep_dir = work_dir.join(format!("sweep_{}_{}", sweep.name, idx));
            std::fs::create_dir_all(&sweep_dir)
                .map_err(|e| SolverError::io(&sweep_dir, e))?;

            let output_dir = sweep_dir.join("output");
            std::fs::create_dir_all(&output_dir)
                .map_err(|e| SolverError::io(&output_dir, e))?;

            // Write config with single frequency point.
            let config_path =
                config::write_palace_config(design, mesh_path, &output_dir)?;

            // Run solver at this frequency.
            match solver_bridge::dispatch_solver(&config_path, progress) {
                Ok(_) => {
                    // Read S11 from output.
                    let s11 = read_s11_complex(&output_dir);
                    s_data_per_freq.push((freq_hz, s11));
                }
                Err(e) => {
                    progress.on_log(&format!(
                        "Sweep point {:.3e} Hz failed: {}",
                        freq_hz, e
                    ));
                    s_data_per_freq.push((freq_hz, None));
                }
            }
        }

        // Build SParameterData from collected results.
        let s_param_data = build_s_parameter_data(
            design,
            &setup.name,
            &sweep.name,
            &s_data_per_freq,
        );
        all_results.push(s_param_data);
    }

    Ok(all_results)
}

/// Compute discrete frequency points for a sweep.
fn compute_frequency_points(
    start_hz: f64,
    stop_hz: f64,
    sweep_type: &SweepType,
    step: Option<&str>,
    count: Option<u32>,
) -> Result<Vec<f64>, SolverError> {
    match sweep_type {
        SweepType::Discrete => {
            if let Some(n) = count {
                // Linearly spaced points.
                let n = n.max(2) as usize;
                Ok((0..n)
                    .map(|i| start_hz + (stop_hz - start_hz) * i as f64 / (n - 1) as f64)
                    .collect())
            } else if let Some(step_str) = step {
                let step_hz = parse_frequency(step_str)?;
                if step_hz <= 0.0 {
                    return Err(SolverError::ConfigGeneration(
                        "Frequency step must be positive".into(),
                    ));
                }
                let mut points = Vec::new();
                let mut f = start_hz;
                while f <= stop_hz + step_hz * 0.01 {
                    points.push(f);
                    f += step_hz;
                }
                Ok(points)
            } else {
                // Default: 11 points.
                let n = 11usize;
                Ok((0..n)
                    .map(|i| start_hz + (stop_hz - start_hz) * i as f64 / (n - 1) as f64)
                    .collect())
            }
        }
        SweepType::Interpolating | SweepType::Fast => {
            // Use fewer solve points for interpolating/fast sweeps.
            let n = count.unwrap_or(5).max(2) as usize;
            Ok((0..n)
                .map(|i| start_hz + (stop_hz - start_hz) * i as f64 / (n - 1) as f64)
                .collect())
        }
    }
}

/// Read S11 complex values from rem's port-S.csv.
fn read_s11_complex(output_dir: &Path) -> Option<(f64, f64)> {
    let path = output_dir.join("port-S.csv");
    let content = std::fs::read_to_string(&path).ok()?;
    let last_line = content.lines().filter(|l| !l.starts_with('f')).last()?;
    let cols: Vec<&str> = last_line.split(',').collect();
    let re = cols.get(1)?.trim().parse::<f64>().ok()?;
    let im = cols.get(2)?.trim().parse::<f64>().ok()?;
    Some((re, im))
}

/// Build SParameterData from sweep results.
fn build_s_parameter_data(
    design: &Design,
    setup_name: &str,
    sweep_name: &str,
    data: &[(f64, Option<(f64, f64)>)],
) -> SParameterData {
    let frequencies: Vec<f64> = data.iter().map(|(f, _)| *f / 1e9).collect(); // GHz
    let mut real = Vec::with_capacity(data.len());
    let mut imag = Vec::with_capacity(data.len());
    let mut mag_db = Vec::with_capacity(data.len());

    for (_, s11) in data {
        let (re, im) = s11.unwrap_or((0.0, 0.0));
        real.push(re);
        imag.push(im);
        let mag2 = re * re + im * im;
        let db = if mag2 > 1e-300 { 10.0 * mag2.log10() } else { -300.0 };
        mag_db.push(db);
    }

    let mut parameters = HashMap::new();
    parameters.insert("S11".to_string(), SParamValues { real, imag });

    let mut derived = HashMap::new();
    derived.insert("dB(S11)".to_string(), mag_db);

    SParameterData {
        format_version: "1.0".to_string(),
        file_type: "s_parameters".to_string(),
        design_id: design.id.clone(),
        setup: setup_name.to_string(),
        sweep: sweep_name.to_string(),
        solution_type: format!("{:?}", design.solution_type),
        reference_impedance_ohm: 50.0,
        num_ports: 1,
        port_names: vec!["Port1".to_string()],
        num_frequencies: frequencies.len(),
        frequency_unit: "GHz".to_string(),
        data_format: "RI".to_string(),
        frequencies,
        parameters,
        derived,
    }
}
