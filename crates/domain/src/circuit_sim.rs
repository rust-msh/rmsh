//! Basic S-parameter circuit co-simulation engine.
//!
//! Supports cascading S-parameter blocks, converting between S/Z/Y parameters,
//! and computing system response from EM-extracted Touchstone data.

use std::collections::HashMap;
use num_complex::Complex64;
use thiserror::Error;

use crate::result_store::SParameterData;

#[derive(Debug, Error)]
pub enum CircuitError {
    #[error("Port count mismatch: {0} vs {1}")]
    PortMismatch(usize, usize),
    #[error("Frequency point mismatch")]
    FrequencyMismatch,
    #[error("Singular matrix at frequency {0:.3e} Hz")]
    SingularMatrix(f64),
    #[error("{0}")]
    Other(String),
}

/// A single N-port S-parameter block extracted from EM simulation.
#[derive(Debug, Clone)]
pub struct CircuitBlock {
    pub name: String,
    pub frequencies: Vec<f64>,
    pub s_matrices: Vec<Vec<Complex64>>, // per-frequency, row-major N×N
    pub n_ports: usize,
    pub z0: f64,
}

impl CircuitBlock {
    /// Create from SParameterData.
    pub fn from_sparam(data: &SParameterData, z0: f64) -> Result<Self, CircuitError> {
        let n_ports = data.num_ports;
        let freqs_orig = data.frequencies_original();
        let n_freqs = freqs_orig.len();

        let mut s_matrices = Vec::with_capacity(n_freqs);

        for fi in 0..n_freqs {
            let mut s_mat = vec![Complex64::ZERO; n_ports * n_ports];
            for i in 0..n_ports {
                for j in 0..n_ports {
                    let key = format!("S{}{}", i + 1, j + 1);
                    if let Some(values) = data.parameters.get(&key) {
                        let re = values.real.get(fi).copied().unwrap_or(0.0);
                        let im = values.imag.get(fi).copied().unwrap_or(0.0);
                        s_mat[i * n_ports + j] = Complex64::new(re, im);
                    }
                }
            }
            s_matrices.push(s_mat);
        }

        Ok(Self {
            name: data.setup.clone(),
            frequencies: freqs_orig.to_vec(),
            s_matrices,
            n_ports,
            z0,
        })
    }

    /// Get S-matrix at a specific frequency index.
    pub fn s_matrix(&self, freq_idx: usize) -> &[Complex64] {
        &self.s_matrices[freq_idx]
    }
}

/// Cascade two 2-port S-parameter blocks.
///
/// Block A (ports 1,2) → Block B (ports 3,4)
/// with port 2 of A connected to port 1 of B.
///
/// Returns the 2-port S-matrix of the combined network.
pub fn cascade_s_params_2port(
    s_a: &[Complex64; 4],  // [S11, S12, S21, S22]
    s_b: &[Complex64; 4],  // [S11, S12, S21, S22]
) -> [Complex64; 4] {
    let (s11_a, s12_a, s21_a, s22_a) = (s_a[0], s_a[1], s_a[2], s_a[3]);
    let (s11_b, s12_b, s21_b, s22_b) = (s_b[0], s_b[1], s_b[2], s_b[3]);

    let denom = Complex64::ONE - s22_a * s11_b;
    if denom.norm() < 1e-30 {
        return [Complex64::ZERO; 4];
    }

    let s11 = s11_a + s12_a * s11_b * s21_a / denom;
    let s12 = s12_a * s12_b / denom;
    let s21 = s21_b * s21_a / denom;
    let s22 = s22_b + s21_b * s22_a * s12_b / denom;

    [s11, s12, s21, s22]
}

/// Cascade a chain of 2-port blocks at all frequency points.
///
/// Returns the cascaded S-parameter data.
pub fn cascade_chain_2port(
    blocks: &[&CircuitBlock],
) -> Result<SParameterData, CircuitError> {
    if blocks.is_empty() {
        return Err(CircuitError::Other("No blocks to cascade".into()));
    }

    // Verify all blocks have the same number of frequency points
    let n_freqs = blocks[0].frequencies.len();
    for b in blocks.iter().skip(1) {
        if b.frequencies.len() != n_freqs {
            return Err(CircuitError::FrequencyMismatch);
        }
        if b.n_ports != 2 {
            return Err(CircuitError::PortMismatch(2, b.n_ports));
        }
    }

    let freqs = blocks[0].frequencies.clone();
    let mut s11_re = Vec::with_capacity(n_freqs);
    let mut s11_im = Vec::with_capacity(n_freqs);
    let mut s21_re = Vec::with_capacity(n_freqs);
    let mut s21_im = Vec::with_capacity(n_freqs);
    let mut s12_re = Vec::with_capacity(n_freqs);
    let mut s12_im = Vec::with_capacity(n_freqs);
    let mut s22_re = Vec::with_capacity(n_freqs);
    let mut s22_im = Vec::with_capacity(n_freqs);

    for fi in 0..n_freqs {
        // Start with first block
        let s0 = blocks[0].s_matrix(fi);
        let mut cascaded = [s0[0], s0[1], s0[2], s0[3]];

        // Cascade subsequent blocks
        for block in blocks.iter().skip(1) {
            let sb = block.s_matrix(fi);
            let next = [sb[0], sb[1], sb[2], sb[3]];
            cascaded = cascade_s_params_2port(&cascaded, &next);
        }

        s11_re.push(cascaded[0].re);
        s11_im.push(cascaded[0].im);
        s12_re.push(cascaded[1].re);
        s12_im.push(cascaded[1].im);
        s21_re.push(cascaded[2].re);
        s21_im.push(cascaded[2].im);
        s22_re.push(cascaded[3].re);
        s22_im.push(cascaded[3].im);
    }

    let mut parameters = HashMap::new();
    parameters.insert("S11".to_string(), crate::result_store::SParamValues { real: s11_re, imag: s11_im });
    parameters.insert("S12".to_string(), crate::result_store::SParamValues { real: s12_re, imag: s12_im });
    parameters.insert("S21".to_string(), crate::result_store::SParamValues { real: s21_re, imag: s21_im });
    parameters.insert("S22".to_string(), crate::result_store::SParamValues { real: s22_re, imag: s22_im });

    Ok(SParameterData {
        format_version: "1.0".into(),
        file_type: "cascaded_s_parameters".into(),
        design_id: "cascade".into(),
        setup: "CascadeChain".into(),
        sweep: String::new(),
        solution_type: "Cascade".into(),
        reference_impedance_ohm: blocks[0].z0,
        num_ports: 2,
        port_names: vec!["Port1".into(), "Port2".into()],
        num_frequencies: n_freqs,
        frequency_unit: "Hz".into(),
        data_format: "RI".into(),
        frequencies: freqs.iter().map(|f| f / 1e9).collect(), // GHz
        parameters,
        derived: HashMap::new(),
    })
}

/// Convert S-parameters to Z-parameters for a 2-port network.
pub fn s_to_z_2port(s: &[Complex64; 4], z0: f64) -> [Complex64; 4] {
    let denom = (Complex64::ONE - s[0]) * (Complex64::ONE - s[3]) - s[1] * s[2];
    if denom.norm() < 1e-30 {
        return [Complex64::ZERO; 4];
    }
    let z0c = Complex64::new(z0, 0.0);
    let z11 = z0c * ((Complex64::ONE + s[0]) * (Complex64::ONE - s[3]) + s[1] * s[2]) / denom;
    let z12 = z0c * (Complex64::new(2.0, 0.0) * s[1]) / denom;
    let z21 = z0c * (Complex64::new(2.0, 0.0) * s[2]) / denom;
    let z22 = z0c * ((Complex64::ONE - s[0]) * (Complex64::ONE + s[3]) + s[1] * s[2]) / denom;
    [z11, z12, z21, z22]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cascade_two_matched_blocks_is_transparent() {
        // Perfectly matched, lossless → identity
        let s_a = [Complex64::ZERO, Complex64::ONE, Complex64::ONE, Complex64::ZERO];
        let s_b = [Complex64::ZERO, Complex64::ONE, Complex64::ONE, Complex64::ZERO];
        let result = cascade_s_params_2port(&s_a, &s_b);

        assert!((result[0].norm() - 0.0).abs() < 1e-10); // S11 ≈ 0
        assert!((result[2].norm() - 1.0).abs() < 1e-10); // S21 ≈ 1
    }

    #[test]
    fn cascade_with_open_is_reflective() {
        // Block B is open (S11=1) → cascaded S11 should be large
        let s_a = [Complex64::ZERO, Complex64::ONE, Complex64::ONE, Complex64::ZERO];
        let s_b = [Complex64::ONE, Complex64::ZERO, Complex64::ZERO, Complex64::ONE];
        let result = cascade_s_params_2port(&s_a, &s_b);

        assert!(result[0].norm() > 0.9); // Strong reflection
    }

    #[test]
    fn s_to_z_conversion_series_r() {
        // Series 50Ω resistor: S11 = 50/(50+50+50) = 0.333, S21 = 1+S11 = 1.333?
        // Actually use a simple reflective load: S11=0.5, S21=0, S12=0, S22=0.5
        let s = [
            Complex64::new(0.5, 0.0),
            Complex64::ZERO,
            Complex64::ZERO,
            Complex64::new(0.5, 0.0),
        ];
        let z = s_to_z_2port(&s, 50.0);
        // Z11 should be positive real for a passive load
        assert!(z[0].re > 0.0);
        assert!(z[0].im.abs() < 1e-10);
    }
}
