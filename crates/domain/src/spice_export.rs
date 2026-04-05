// ---------------------------------------------------------------------------
// SPICE Export — Generate equivalent circuit netlists from RLCG data
// ---------------------------------------------------------------------------

use std::fmt::Write;
use std::path::Path;

use crate::result_store::{RlcgMatrixData, ResultError};

/// SPICE equivalent circuit model type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpiceModelType {
    /// Pi-model: shunt-series-shunt
    PiModel,
    /// T-model: series-shunt-series
    TModel,
}

/// Options for SPICE netlist export.
#[derive(Debug, Clone)]
pub struct SpiceExportOptions {
    /// Frequency index to use for the lumped model.
    pub frequency_idx: usize,
    /// Circuit model topology.
    pub model_type: SpiceModelType,
    /// Include mutual inductance coupling.
    pub include_mutual: bool,
}

impl Default for SpiceExportOptions {
    fn default() -> Self {
        Self {
            frequency_idx: 0,
            model_type: SpiceModelType::PiModel,
            include_mutual: true,
        }
    }
}

/// Generate a SPICE netlist string from RLCG matrix data.
pub fn export_spice_netlist(
    rlcg: &RlcgMatrixData,
    options: &SpiceExportOptions,
) -> Result<String, ResultError> {
    let n = rlcg.num_nets;
    if n == 0 {
        return Err(ResultError::InvalidData("No nets in RLCG data".to_string()));
    }

    let freq_idx = options.frequency_idx;

    // Get matrices at the specified frequency
    let r_mat = rlcg.matrix_at_frequency("R", freq_idx);
    let l_mat = rlcg.matrix_at_frequency("L", freq_idx);
    let c_mat = rlcg.matrix_at_frequency("C", freq_idx);
    let g_mat = rlcg.matrix_at_frequency("G", freq_idx);

    let freq = rlcg.frequencies.get(freq_idx).copied().unwrap_or(0.0);

    let mut netlist = String::new();
    writeln!(netlist, "* EMStudio Q3D Equivalent Circuit").unwrap();
    writeln!(netlist, "* Design: {}", rlcg.design_id).unwrap();
    writeln!(netlist, "* Setup: {}", rlcg.setup).unwrap();
    writeln!(netlist, "* Frequency: {} GHz", freq).unwrap();
    writeln!(netlist, "* Model: {:?}", options.model_type).unwrap();
    writeln!(netlist, "* Nets: {}", rlcg.net_names.join(", ")).unwrap();
    writeln!(netlist).unwrap();
    writeln!(netlist, ".SUBCKT Q3D_EQUIV {} 0",
        (1..=n).map(|i| format!("net{}", i)).collect::<Vec<_>>().join(" ")).unwrap();
    writeln!(netlist).unwrap();

    // Self-impedance elements (diagonal)
    for i in 0..n {
        let net_name = &rlcg.net_names[i];
        let node_in = format!("net{}", i + 1);
        let node_mid = format!("mid{}", i + 1);
        let node_out = format!("out{}", i + 1);

        writeln!(netlist, "* --- Net: {} ---", net_name).unwrap();

        // Series resistance
        if let Some(r) = r_mat.and_then(|m| m.get(i).and_then(|r| r.get(i))) {
            writeln!(netlist, "R_self_{} {} {} {:.6e}", i + 1, node_in, node_mid, r).unwrap();
        }

        // Series inductance
        if let Some(l) = l_mat.and_then(|m| m.get(i).and_then(|r| r.get(i))) {
            writeln!(netlist, "L_self_{} {} {} {:.6e}", i + 1, node_mid, node_out, l).unwrap();
        }

        match options.model_type {
            SpiceModelType::PiModel => {
                // Shunt capacitance at input and output (C/2 each)
                if let Some(c) = c_mat.and_then(|m| m.get(i).and_then(|r| r.get(i))) {
                    writeln!(netlist, "C_in_{} {} 0 {:.6e}", i + 1, node_in, c / 2.0).unwrap();
                    writeln!(netlist, "C_out_{} {} 0 {:.6e}", i + 1, node_out, c / 2.0).unwrap();
                }
                // Shunt conductance at input (G/2 each)
                if let Some(g) = g_mat.and_then(|m| m.get(i).and_then(|r| r.get(i))) {
                    if g.abs() > 1e-15 {
                        writeln!(netlist, "G_in_{} {} 0 {:.6e}", i + 1, node_in, g / 2.0).unwrap();
                        writeln!(netlist, "G_out_{} {} 0 {:.6e}", i + 1, node_out, g / 2.0).unwrap();
                    }
                }
            }
            SpiceModelType::TModel => {
                // Series elements split R/2 - shunt C - R/2
                // (simplified T-model)
                if let Some(c) = c_mat.and_then(|m| m.get(i).and_then(|r| r.get(i))) {
                    writeln!(netlist, "C_shunt_{} {} 0 {:.6e}", i + 1, node_mid, c).unwrap();
                }
                if let Some(g) = g_mat.and_then(|m| m.get(i).and_then(|r| r.get(i))) {
                    if g.abs() > 1e-15 {
                        writeln!(netlist, "G_shunt_{} {} 0 {:.6e}", i + 1, node_mid, g).unwrap();
                    }
                }
            }
        }

        writeln!(netlist).unwrap();
    }

    // Mutual coupling (off-diagonal)
    if options.include_mutual && n > 1 {
        writeln!(netlist, "* --- Mutual Coupling ---").unwrap();
        if let Some(l_matrix) = l_mat {
            for i in 0..n {
                for j in (i + 1)..n {
                    let l_ii = l_matrix.get(i).and_then(|r| r.get(i)).copied().unwrap_or(0.0);
                    let l_jj = l_matrix.get(j).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                    let l_ij = l_matrix.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);

                    if l_ii > 1e-15 && l_jj > 1e-15 {
                        let k = l_ij / (l_ii * l_jj).sqrt();
                        if k.abs() > 1e-6 {
                            writeln!(
                                netlist,
                                "K_{}{} L_self_{} L_self_{} {:.6e}",
                                i + 1,
                                j + 1,
                                i + 1,
                                j + 1,
                                k
                            )
                            .unwrap();
                        }
                    }
                }
            }
        }

        // Mutual capacitance
        if let Some(c_matrix) = c_mat {
            for i in 0..n {
                for j in (i + 1)..n {
                    let c_ij = c_matrix.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                    if c_ij.abs() > 1e-18 {
                        writeln!(
                            netlist,
                            "C_mutual_{}_{} net{} net{} {:.6e}",
                            i + 1,
                            j + 1,
                            i + 1,
                            j + 1,
                            c_ij.abs()
                        )
                        .unwrap();
                    }
                }
            }
        }

        writeln!(netlist).unwrap();
    }

    writeln!(netlist, ".ENDS Q3D_EQUIV").unwrap();
    writeln!(netlist, ".END").unwrap();

    Ok(netlist)
}

/// Save SPICE netlist to a file.
pub fn save_spice_file(
    path: &Path,
    rlcg: &RlcgMatrixData,
    options: &SpiceExportOptions,
) -> Result<(), ResultError> {
    let netlist = export_spice_netlist(rlcg, options)?;
    std::fs::write(path, netlist).map_err(ResultError::Io)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::result_store::{RlcgMatrix, RlcgFrequencyPoint};
    use std::collections::HashMap;

    fn make_test_rlcg() -> RlcgMatrixData {
        let mut matrices = HashMap::new();
        matrices.insert(
            "R".to_string(),
            RlcgMatrix {
                description: "Resistance".to_string(),
                unit: "ohm".to_string(),
                data_per_frequency: vec![RlcgFrequencyPoint {
                    frequency: 1.0,
                    matrix: vec![vec![0.125, 0.003], vec![0.003, 0.125]],
                }],
            },
        );
        matrices.insert(
            "L".to_string(),
            RlcgMatrix {
                description: "Inductance".to_string(),
                unit: "H".to_string(),
                data_per_frequency: vec![RlcgFrequencyPoint {
                    frequency: 1.0,
                    matrix: vec![vec![1e-9, 0.1e-9], vec![0.1e-9, 1e-9]],
                }],
            },
        );
        matrices.insert(
            "C".to_string(),
            RlcgMatrix {
                description: "Capacitance".to_string(),
                unit: "F".to_string(),
                data_per_frequency: vec![RlcgFrequencyPoint {
                    frequency: 1.0,
                    matrix: vec![vec![1e-12, 0.05e-12], vec![0.05e-12, 1e-12]],
                }],
            },
        );

        RlcgMatrixData {
            format_version: "1.0".to_string(),
            file_type: "RLCGMatrixData".to_string(),
            design_id: "test-design".to_string(),
            setup: "Q3D_Setup1".to_string(),
            solution_type: "Q3D_ACRL".to_string(),
            num_nets: 2,
            net_names: vec!["Signal".to_string(), "Ground".to_string()],
            terminal_names: vec!["T1".to_string(), "T2".to_string()],
            frequencies: vec![1.0],
            matrices,
            dc_data: None,
        }
    }

    #[test]
    fn generate_spice_netlist() {
        let rlcg = make_test_rlcg();
        let options = SpiceExportOptions::default();
        let netlist = export_spice_netlist(&rlcg, &options).unwrap();

        assert!(netlist.contains("Q3D_EQUIV"));
        assert!(netlist.contains("R_self_1"));
        assert!(netlist.contains("L_self_1"));
        assert!(netlist.contains("C_in_1"));
        assert!(netlist.contains("K_12")); // mutual coupling
        assert!(netlist.contains(".END"));
    }

    #[test]
    fn spice_t_model() {
        let rlcg = make_test_rlcg();
        let options = SpiceExportOptions {
            model_type: SpiceModelType::TModel,
            ..Default::default()
        };
        let netlist = export_spice_netlist(&rlcg, &options).unwrap();

        assert!(netlist.contains("C_shunt_1"));
        assert!(!netlist.contains("C_in_1"));
    }
}
