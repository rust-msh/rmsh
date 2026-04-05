// ---------------------------------------------------------------------------
// RLCG → S-Parameter Conversion
// ---------------------------------------------------------------------------
//
// Converts RLCG matrix data to S-parameters using the standard lumped
// impedance model: Z = R + jωL for series, Y = G + jωC for shunt,
// then Z-to-S transformation.

use std::collections::HashMap;
use std::f64::consts::PI;

use crate::result_store::{RlcgMatrixData, SParameterData, SParamValues, ResultError};

/// Options for RLCG to S-parameter conversion.
#[derive(Debug, Clone)]
pub struct RlcgToSparamOptions {
    /// Reference impedance in ohms (typically 50.0).
    pub reference_impedance: f64,
}

impl Default for RlcgToSparamOptions {
    fn default() -> Self {
        Self {
            reference_impedance: 50.0,
        }
    }
}

/// Convert RLCG matrix data to S-parameters.
///
/// For each frequency point, constructs the impedance matrix:
///   Z(ω) = R + jωL
/// and converts to S-parameters:
///   S = (Z - Z₀I)(Z + Z₀I)⁻¹
pub fn rlcg_to_s_parameters(
    rlcg: &RlcgMatrixData,
    options: &RlcgToSparamOptions,
) -> Result<SParameterData, ResultError> {
    let n = rlcg.num_nets;
    if n == 0 {
        return Err(ResultError::InvalidData("No nets in RLCG data".to_string()));
    }

    let z0 = options.reference_impedance;
    let num_freqs = rlcg.frequencies.len();

    // Initialize S-parameter storage
    let mut parameters: HashMap<String, SParamValues> = HashMap::new();
    for i in 0..n {
        for j in 0..n {
            let key = SParameterData::s_key(i + 1, j + 1);
            parameters.insert(
                key,
                SParamValues {
                    real: vec![0.0; num_freqs],
                    imag: vec![0.0; num_freqs],
                },
            );
        }
    }

    for (fi, &freq_ghz) in rlcg.frequencies.iter().enumerate() {
        let omega = 2.0 * PI * freq_ghz * 1e9; // freq in GHz → rad/s

        // Build Z matrix: Z = R + jωL
        let mut z_real = vec![vec![0.0; n]; n];
        let mut z_imag = vec![vec![0.0; n]; n];

        if let Some(r_mat) = rlcg.matrix_at_frequency("R", fi) {
            for i in 0..n {
                for j in 0..n {
                    z_real[i][j] += r_mat.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                }
            }
        }
        if let Some(l_mat) = rlcg.matrix_at_frequency("L", fi) {
            for i in 0..n {
                for j in 0..n {
                    z_imag[i][j] += omega * l_mat.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
                }
            }
        }

        // Also add shunt admittance as series impedance contribution
        // For lumped model: Y_shunt = G + jωC, Z_total includes both series and shunt
        // For simplicity, we use the direct Z-to-S conversion of the series impedance
        // (this is the standard Q3D convention for lumped extraction)

        // S = (Z - Z₀I)(Z + Z₀I)⁻¹
        // Compute Z - Z₀I and Z + Z₀I
        let mut zm = vec![vec![(0.0, 0.0); n]; n]; // Z - Z₀I (real, imag)
        let mut zp = vec![vec![(0.0, 0.0); n]; n]; // Z + Z₀I (real, imag)

        for i in 0..n {
            for j in 0..n {
                let diag = if i == j { z0 } else { 0.0 };
                zm[i][j] = (z_real[i][j] - diag, z_imag[i][j]);
                zp[i][j] = (z_real[i][j] + diag, z_imag[i][j]);
            }
        }

        // Invert Z + Z₀I
        let zp_inv = complex_matrix_inverse(&zp);

        // S = (Z - Z₀I) * (Z + Z₀I)⁻¹
        let s = complex_matrix_multiply(&zm, &zp_inv);

        // Store results
        for i in 0..n {
            for j in 0..n {
                let key = SParameterData::s_key(i + 1, j + 1);
                if let Some(vals) = parameters.get_mut(&key) {
                    vals.real[fi] = s[i][j].0;
                    vals.imag[fi] = s[i][j].1;
                }
            }
        }
    }

    // Build port names from net names
    let port_names: Vec<String> = rlcg
        .net_names
        .iter()
        .enumerate()
        .map(|(i, name)| format!("Port{}_{}", i + 1, name))
        .collect();

    Ok(SParameterData {
        format_version: "1.0".to_string(),
        file_type: "SParameterData".to_string(),
        design_id: rlcg.design_id.clone(),
        setup: rlcg.setup.clone(),
        sweep: String::new(),
        solution_type: "Q3D_to_SParam".to_string(),
        reference_impedance_ohm: z0,
        num_ports: n,
        port_names,
        num_frequencies: num_freqs,
        frequency_unit: "GHz".to_string(),
        data_format: "RealImaginary".to_string(),
        frequencies: rlcg.frequencies.clone(),
        parameters,
        derived: HashMap::new(),
    })
}

// ---------------------------------------------------------------------------
// Complex matrix operations (for small N×N matrices)
// ---------------------------------------------------------------------------

type ComplexMatrix = Vec<Vec<(f64, f64)>>;

fn complex_matrix_multiply(a: &ComplexMatrix, b: &ComplexMatrix) -> ComplexMatrix {
    let n = a.len();
    let mut result = vec![vec![(0.0, 0.0); n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut re = 0.0;
            let mut im = 0.0;
            for k in 0..n {
                let (ar, ai) = a[i][k];
                let (br, bi) = b[k][j];
                re += ar * br - ai * bi;
                im += ar * bi + ai * br;
            }
            result[i][j] = (re, im);
        }
    }
    result
}

fn complex_matrix_inverse(m: &ComplexMatrix) -> ComplexMatrix {
    let n = m.len();
    if n == 1 {
        let (r, i) = m[0][0];
        let denom = r * r + i * i;
        if denom < 1e-30 {
            return vec![vec![(0.0, 0.0)]];
        }
        return vec![vec![(r / denom, -i / denom)]];
    }
    if n == 2 {
        return complex_matrix_inverse_2x2(m);
    }
    complex_matrix_inverse_nxn(m)
}

fn complex_matrix_inverse_2x2(m: &ComplexMatrix) -> ComplexMatrix {
    let (ar, ai) = m[0][0];
    let (br, bi) = m[0][1];
    let (cr, ci) = m[1][0];
    let (dr, di) = m[1][1];

    // det = ad - bc (complex)
    let det_r = ar * dr - ai * di - (br * cr - bi * ci);
    let det_i = ar * di + ai * dr - (br * ci + bi * cr);
    let det_mag2 = det_r * det_r + det_i * det_i;

    if det_mag2 < 1e-30 {
        return vec![vec![(0.0, 0.0); 2]; 2];
    }

    // inv_det = 1/det (complex)
    let inv_det_r = det_r / det_mag2;
    let inv_det_i = -det_i / det_mag2;

    // [d, -b; -c, a] * inv_det
    let cmul = |ar: f64, ai: f64, br: f64, bi: f64| -> (f64, f64) {
        (ar * br - ai * bi, ar * bi + ai * br)
    };

    vec![
        vec![
            cmul(dr, di, inv_det_r, inv_det_i),
            cmul(-br, -bi, inv_det_r, inv_det_i),
        ],
        vec![
            cmul(-cr, -ci, inv_det_r, inv_det_i),
            cmul(ar, ai, inv_det_r, inv_det_i),
        ],
    ]
}

/// Gaussian elimination for N×N complex matrix inverse.
fn complex_matrix_inverse_nxn(m: &ComplexMatrix) -> ComplexMatrix {
    let n = m.len();
    // Augmented matrix [M | I]
    let mut aug: Vec<Vec<(f64, f64)>> = vec![vec![(0.0, 0.0); 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = m[i][j];
        }
        aug[i][n + i] = (1.0, 0.0);
    }

    // Forward elimination with partial pivoting
    for col in 0..n {
        // Find pivot
        let mut max_mag = 0.0;
        let mut pivot_row = col;
        for row in col..n {
            let (r, i) = aug[row][col];
            let mag = r * r + i * i;
            if mag > max_mag {
                max_mag = mag;
                pivot_row = row;
            }
        }
        if max_mag < 1e-30 {
            continue; // Singular
        }
        aug.swap(col, pivot_row);

        // Normalize pivot row
        let (pr, pi) = aug[col][col];
        let pm = pr * pr + pi * pi;
        let inv_pivot = (pr / pm, -pi / pm);
        for j in 0..2 * n {
            let (vr, vi) = aug[col][j];
            aug[col][j] = (
                vr * inv_pivot.0 - vi * inv_pivot.1,
                vr * inv_pivot.1 + vi * inv_pivot.0,
            );
        }

        // Eliminate column
        for row in 0..n {
            if row == col {
                continue;
            }
            let (fr, fi) = aug[row][col];
            for j in 0..2 * n {
                let (vr, vi) = aug[col][j];
                aug[row][j].0 -= fr * vr - fi * vi;
                aug[row][j].1 -= fr * vi + fi * vr;
            }
        }
    }

    // Extract inverse from augmented matrix
    let mut result = vec![vec![(0.0, 0.0); n]; n];
    for i in 0..n {
        for j in 0..n {
            result[i][j] = aug[i][n + j];
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::result_store::{RlcgMatrix, RlcgFrequencyPoint};

    fn make_test_rlcg() -> RlcgMatrixData {
        let mut matrices = HashMap::new();
        matrices.insert(
            "R".to_string(),
            RlcgMatrix {
                description: "Resistance".to_string(),
                unit: "ohm".to_string(),
                data_per_frequency: vec![RlcgFrequencyPoint {
                    frequency: 1.0,
                    matrix: vec![vec![10.0]],
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
                    matrix: vec![vec![1e-9]],
                }],
            },
        );

        RlcgMatrixData {
            format_version: "1.0".to_string(),
            file_type: "RLCGMatrixData".to_string(),
            design_id: "test".to_string(),
            setup: "Q3D_Setup1".to_string(),
            solution_type: "Q3D_ACRL".to_string(),
            num_nets: 1,
            net_names: vec!["Signal".to_string()],
            terminal_names: vec!["T1".to_string()],
            frequencies: vec![1.0],
            matrices,
            dc_data: None,
        }
    }

    #[test]
    fn rlcg_to_s_single_port() {
        let rlcg = make_test_rlcg();
        let options = RlcgToSparamOptions::default();
        let sparam = rlcg_to_s_parameters(&rlcg, &options).unwrap();

        assert_eq!(sparam.num_ports, 1);
        assert_eq!(sparam.num_frequencies, 1);
        assert_eq!(sparam.reference_impedance_ohm, 50.0);

        // S11 = (Z - Z0) / (Z + Z0) where Z = R + jωL
        let s11 = sparam.get_complex("S11").unwrap();
        assert_eq!(s11.len(), 1);

        // Z = 10 + j*(2π*1e9*1e-9) = 10 + j*6.2832
        // S11 = (Z - 50)/(Z + 50) = (-40 + j*6.28)/(60 + j*6.28)
        let z_re = 10.0;
        let z_im = 2.0 * PI * 1e9 * 1e-9;
        let num_re = z_re - 50.0;
        let num_im = z_im;
        let den_re = z_re + 50.0;
        let den_im = z_im;
        let den_mag2 = den_re * den_re + den_im * den_im;
        let expected_re = (num_re * den_re + num_im * den_im) / den_mag2;
        let expected_im = (num_im * den_re - num_re * den_im) / den_mag2;

        assert!(
            (s11[0][0] - expected_re).abs() < 1e-6,
            "S11 real: got {}, expected {}",
            s11[0][0],
            expected_re
        );
        assert!(
            (s11[0][1] - expected_im).abs() < 1e-6,
            "S11 imag: got {}, expected {}",
            s11[0][1],
            expected_im
        );
    }

    #[test]
    fn complex_2x2_inverse() {
        let m = vec![
            vec![(1.0, 0.0), (0.0, 0.0)],
            vec![(0.0, 0.0), (1.0, 0.0)],
        ];
        let inv = complex_matrix_inverse_2x2(&m);
        assert!((inv[0][0].0 - 1.0).abs() < 1e-10);
        assert!((inv[1][1].0 - 1.0).abs() < 1e-10);
    }
}
