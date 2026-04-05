//! Result extraction from rem solver outputs.
//!
//! Reads rem's CSV/VTK output files and converts them to emstudio's domain
//! result types (SParameterData, ConvergenceData, RlcgMatrixData, FarFieldData).

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use emstudio_domain::Design;
use emstudio_domain::result_store::{
    AngleRange, AntennaParameters, ConvergenceData,
    DerivedQuantity, FarFieldData, RlcgFrequencyPoint, RlcgMatrix,
    RlcgMatrixData, SParamValues, SParameterData,
};

use crate::error::SolverError;

// ---------------------------------------------------------------------------
// S-parameter extraction from port-S.csv
// ---------------------------------------------------------------------------

/// Read S-parameter data from rem's `port-S.csv` output.
///
/// CSV format: `f (Hz),Re(S[1][1]),Im(S[1][1]),|S[1][1]| (dB)`
pub fn extract_s_parameters(
    output_dir: &Path,
    design: &Design,
    setup_name: &str,
) -> Result<SParameterData, SolverError> {
    let csv_path = output_dir.join("port-S.csv");
    let content = std::fs::read_to_string(&csv_path)
        .map_err(|e| SolverError::io(&csv_path, e))?;

    let mut frequencies = Vec::new();
    let mut s11_real = Vec::new();
    let mut s11_imag = Vec::new();
    let mut s11_db = Vec::new();

    for line in content.lines().skip(1) {
        // skip header
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() < 4 {
            continue;
        }
        let f: f64 = cols[0].trim().parse().unwrap_or(0.0);
        let re: f64 = cols[1].trim().parse().unwrap_or(0.0);
        let im: f64 = cols[2].trim().parse().unwrap_or(0.0);
        let db: f64 = cols[3].trim().parse().unwrap_or(-300.0);

        frequencies.push(f / 1e9); // Convert Hz → GHz
        s11_real.push(re);
        s11_imag.push(im);
        s11_db.push(db);
    }

    let mut parameters = HashMap::new();
    parameters.insert(
        "S11".to_string(),
        SParamValues {
            real: s11_real,
            imag: s11_imag,
        },
    );

    let mut derived = HashMap::new();
    derived.insert("dB(S11)".to_string(), s11_db);

    Ok(SParameterData {
        format_version: "1.0".to_string(),
        file_type: "s_parameters".to_string(),
        design_id: design.id.clone(),
        setup: setup_name.to_string(),
        sweep: String::new(),
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
    })
}

// ---------------------------------------------------------------------------
// RLCG matrix extraction from terminal-C.csv
// ---------------------------------------------------------------------------

/// Read capacitance matrix from rem's `postpro/terminal-C.csv`.
///
/// CSV format: `"Frequency (GHz)","C[1][1] (F)","C[1][2] (F)",...`
pub fn extract_rlcg_matrix(
    output_dir: &Path,
    design: &Design,
    setup_name: &str,
) -> Result<RlcgMatrixData, SolverError> {
    let csv_path = output_dir.join("postpro").join("terminal-C.csv");
    let content = std::fs::read_to_string(&csv_path)
        .map_err(|e| SolverError::io(&csv_path, e))?;

    let lines: Vec<&str> = content.lines().collect();
    if lines.len() < 2 {
        return Err(SolverError::ResultExtraction(
            "terminal-C.csv has no data rows".into(),
        ));
    }

    // Parse header to determine matrix size.
    let header = lines[0];
    let header_cols: Vec<&str> = header.split(',').collect();
    // Count C[i][j] entries to determine N.
    let n_entries = header_cols.len() - 1; // subtract frequency column
    let n = (n_entries as f64).sqrt() as usize;

    // Parse data rows.
    let mut frequencies = Vec::new();
    let mut c_per_freq = Vec::new();

    for line in &lines[1..] {
        let cols: Vec<&str> = line.split(',').collect();
        if cols.is_empty() {
            continue;
        }
        let freq: f64 = cols[0].trim().parse().unwrap_or(0.0);
        frequencies.push(freq);

        let mut matrix = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                let idx = 1 + i * n + j;
                if let Some(val) = cols.get(idx) {
                    matrix[i][j] = val.trim().parse().unwrap_or(0.0);
                }
            }
        }
        c_per_freq.push(RlcgFrequencyPoint {
            frequency: freq,
            matrix,
        });
    }

    let net_names: Vec<String> = (1..=n).map(|i| format!("Net{}", i)).collect();
    let terminal_names: Vec<String> = (1..=n).map(|i| format!("Terminal{}", i)).collect();

    let mut matrices = HashMap::new();
    matrices.insert(
        "C".to_string(),
        RlcgMatrix {
            description: "Capacitance matrix".to_string(),
            unit: "F".to_string(),
            data_per_frequency: c_per_freq,
        },
    );

    Ok(RlcgMatrixData {
        format_version: "1.0".to_string(),
        file_type: "rlcg_matrix".to_string(),
        design_id: design.id.clone(),
        setup: setup_name.to_string(),
        solution_type: format!("{:?}", design.solution_type),
        num_nets: n,
        net_names,
        terminal_names,
        frequencies,
        matrices,
        dc_data: None,
    })
}

// ---------------------------------------------------------------------------
// Far-field / RCS extraction from rcs.csv
// ---------------------------------------------------------------------------

/// Read RCS far-field data from rem's `postpro/rcs.csv`.
///
/// CSV format: `Freq (GHz),Theta (deg),Phi (deg),RCS (dBsm)`
pub fn extract_far_field(
    output_dir: &Path,
    design: &Design,
    setup_name: &str,
) -> Result<FarFieldData, SolverError> {
    let csv_path = output_dir.join("postpro").join("rcs.csv");
    let content = std::fs::read_to_string(&csv_path)
        .map_err(|e| SolverError::io(&csv_path, e))?;

    let mut theta_vals: Vec<f64> = Vec::new();
    let mut phi_vals: Vec<f64> = Vec::new();
    let mut rcs_data: Vec<f64> = Vec::new();
    let mut freq_ghz = 0.0;

    for line in content.lines().skip(1) {
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() < 4 {
            continue;
        }
        freq_ghz = cols[0].trim().parse().unwrap_or(0.0);
        let theta: f64 = cols[1].trim().parse().unwrap_or(0.0);
        let phi: f64 = cols[2].trim().parse().unwrap_or(0.0);
        let rcs: f64 = cols[3].trim().parse().unwrap_or(-999.0);

        if !theta_vals.contains(&theta) {
            theta_vals.push(theta);
        }
        if !phi_vals.contains(&phi) {
            phi_vals.push(phi);
        }
        rcs_data.push(rcs);
    }

    theta_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    phi_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let theta_step = if theta_vals.len() > 1 {
        theta_vals[1] - theta_vals[0]
    } else {
        1.0
    };
    let phi_step = if phi_vals.len() > 1 {
        phi_vals[1] - phi_vals[0]
    } else {
        1.0
    };

    let mut derived_quantities = HashMap::new();
    derived_quantities.insert(
        "RCS_dBsm".to_string(),
        DerivedQuantity {
            unit: "dBsm".to_string(),
            data: rcs_data,
        },
    );

    // Compute peak gain from RCS (simplified).
    let peak_rcs = derived_quantities
        .get("RCS_dBsm")
        .map(|dq| dq.data.iter().copied().fold(f64::NEG_INFINITY, f64::max))
        .unwrap_or(-999.0);

    Ok(FarFieldData {
        format_version: "1.0".to_string(),
        file_type: "far_field".to_string(),
        design_id: design.id.clone(),
        setup: setup_name.to_string(),
        far_field_setup: "default".to_string(),
        frequency: format!("{:.6}GHz", freq_ghz),
        theta: AngleRange {
            start_deg: *theta_vals.first().unwrap_or(&0.0),
            stop_deg: *theta_vals.last().unwrap_or(&180.0),
            step_deg: theta_step,
            num_points: theta_vals.len(),
        },
        phi: AngleRange {
            start_deg: *phi_vals.first().unwrap_or(&0.0),
            stop_deg: *phi_vals.last().unwrap_or(&360.0),
            step_deg: phi_step,
            num_points: phi_vals.len(),
        },
        fields: HashMap::new(),
        derived_quantities,
        antenna_parameters: Some(AntennaParameters {
            peak_gain_dbi: Some(peak_rcs),
            peak_gain_theta_deg: None,
            peak_gain_phi_deg: None,
            radiation_efficiency: None,
            beamwidth_e_plane_deg: None,
            beamwidth_h_plane_deg: None,
        }),
    })
}

// ---------------------------------------------------------------------------
// Field data export (.emsfld binary format)
// ---------------------------------------------------------------------------

/// Write field data from VTK output to .emsfld binary format.
///
/// Reads rem's VTK output and converts to emstudio's binary field format.
pub fn export_field_data(
    output_dir: &Path,
    work_dir: &Path,
    _design: &Design,
) -> Result<Option<PathBuf>, SolverError> {
    // Look for VTK files in the output directory.
    let vtk_patterns = ["solution.vtk", "driven_0001.vtk"];
    let mut vtk_path = None;

    for pattern in &vtk_patterns {
        let candidate = output_dir.join("paraview").join(pattern);
        if candidate.exists() {
            vtk_path = Some(candidate);
            break;
        }
        let candidate = output_dir.join(pattern);
        if candidate.exists() {
            vtk_path = Some(candidate);
            break;
        }
    }

    let vtk_path = match vtk_path {
        Some(p) => p,
        None => return Ok(None), // No VTK output available.
    };

    // Parse VTK file and extract field data.
    let content = std::fs::read_to_string(&vtk_path)
        .map_err(|e| SolverError::io(&vtk_path, e))?;

    let (nodes, field_values) = parse_vtk_scalar_field(&content)?;

    if nodes.is_empty() || field_values.is_empty() {
        return Ok(None);
    }

    // Write .emsfld binary file.
    let emsfld_path = work_dir.join("field_data.emsfld");
    write_emsfld(&emsfld_path, &nodes, &field_values)?;

    Ok(Some(emsfld_path))
}

/// Parse a VTK legacy ASCII file and extract scalar point data.
fn parse_vtk_scalar_field(content: &str) -> Result<(Vec<[f64; 3]>, Vec<f64>), SolverError> {
    let mut nodes = Vec::new();
    let mut field_values = Vec::new();
    let mut in_points = false;
    let mut in_scalars = false;
    let mut points_remaining = 0usize;
    let mut scalars_remaining = 0usize;
    let mut skip_lookup = false;

    for line in content.lines() {
        let trimmed = line.trim();

        if trimmed.starts_with("POINTS") {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            points_remaining = parts.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
            in_points = true;
            continue;
        }

        if in_points && points_remaining > 0 {
            let vals: Vec<f64> = trimmed
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if vals.len() >= 3 {
                nodes.push([vals[0], vals[1], vals[2]]);
                points_remaining -= 1;
            }
            if points_remaining == 0 {
                in_points = false;
            }
            continue;
        }

        if trimmed.starts_with("SCALARS") {
            scalars_remaining = nodes.len();
            in_scalars = false;
            skip_lookup = true;
            continue;
        }

        if skip_lookup && trimmed.starts_with("LOOKUP_TABLE") {
            skip_lookup = false;
            in_scalars = true;
            continue;
        }

        if in_scalars && scalars_remaining > 0 {
            if let Ok(val) = trimmed.parse::<f64>() {
                field_values.push(val);
                scalars_remaining -= 1;
            }
            if scalars_remaining == 0 {
                in_scalars = false;
            }
            continue;
        }
    }

    Ok((nodes, field_values))
}

/// Write .emsfld binary file.
///
/// Format: 128-byte header + frequency table + block index + field data blocks.
fn write_emsfld(
    path: &Path,
    _nodes: &[[f64; 3]],
    field_values: &[f64],
) -> Result<(), SolverError> {
    let mut file = std::fs::File::create(path).map_err(|e| SolverError::io(path, e))?;

    // Header: 128 bytes
    let mut header = [0u8; 128];
    // Magic: "EMSFLD\0\0"
    header[0..6].copy_from_slice(b"EMSFLD");
    // Version: 1
    header[8] = 1;
    // Field type: 0 = E-field
    header[9] = 0;
    // Data type: 0 = complex f64
    header[10] = 0;
    // Num components: 1 (scalar)
    header[11] = 1;
    // Num frequencies: 1
    let num_freqs: u32 = 1;
    header[12..16].copy_from_slice(&num_freqs.to_le_bytes());
    // Num nodes
    let num_nodes: u64 = field_values.len() as u64;
    header[16..24].copy_from_slice(&num_nodes.to_le_bytes());

    file.write_all(&header).map_err(|e| SolverError::io(path, e))?;

    // Frequency table: one f64
    let freq: f64 = 0.0;
    file.write_all(&freq.to_le_bytes())
        .map_err(|e| SolverError::io(path, e))?;

    // Block index: offset of the single block
    let block_offset: u64 = 128 + 8 + 8; // header + freq_table + index
    file.write_all(&block_offset.to_le_bytes())
        .map_err(|e| SolverError::io(path, e))?;

    // Field data block: real values (imaginary = 0 for real-valued fields).
    for &val in field_values {
        file.write_all(&val.to_le_bytes())
            .map_err(|e| SolverError::io(path, e))?;
        let imag: f64 = 0.0;
        file.write_all(&imag.to_le_bytes())
            .map_err(|e| SolverError::io(path, e))?;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Convergence data serialization
// ---------------------------------------------------------------------------

/// Write convergence data to JSON file.
pub fn write_convergence_json(
    convergence: &ConvergenceData,
    output_dir: &Path,
) -> Result<PathBuf, SolverError> {
    let path = output_dir.join("convergence.json");
    let json = serde_json::to_string_pretty(convergence)
        .map_err(|e| SolverError::ResultExtraction(e.to_string()))?;
    std::fs::write(&path, json).map_err(|e| SolverError::io(&path, e))?;
    Ok(path)
}

/// Write S-parameter data to JSON file.
pub fn write_s_parameters_json(
    data: &SParameterData,
    output_dir: &Path,
) -> Result<PathBuf, SolverError> {
    let path = output_dir.join("s_parameters.json");
    let json = serde_json::to_string_pretty(data)
        .map_err(|e| SolverError::ResultExtraction(e.to_string()))?;
    std::fs::write(&path, json).map_err(|e| SolverError::io(&path, e))?;
    Ok(path)
}

/// Write RLCG matrix data to JSON file.
pub fn write_rlcg_json(
    data: &RlcgMatrixData,
    output_dir: &Path,
) -> Result<PathBuf, SolverError> {
    let path = output_dir.join("rlcg_matrix.json");
    let json = serde_json::to_string_pretty(data)
        .map_err(|e| SolverError::ResultExtraction(e.to_string()))?;
    std::fs::write(&path, json).map_err(|e| SolverError::io(&path, e))?;
    Ok(path)
}

/// Write far-field data to JSON file.
pub fn write_far_field_json(
    data: &FarFieldData,
    output_dir: &Path,
) -> Result<PathBuf, SolverError> {
    let path = output_dir.join("far_field.json");
    let json = serde_json::to_string_pretty(data)
        .map_err(|e| SolverError::ResultExtraction(e.to_string()))?;
    std::fs::write(&path, json).map_err(|e| SolverError::io(&path, e))?;
    Ok(path)
}
