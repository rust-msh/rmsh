// ---------------------------------------------------------------------------
// Sensitivity Analysis — One-at-a-time perturbation
// ---------------------------------------------------------------------------

use std::collections::HashMap;

use thiserror::Error;

use crate::optimetrics::SensitivityVariable;

#[derive(Debug, Error)]
pub enum SensitivityError {
    #[error("Evaluation failed: {0}")]
    EvalFailed(String),
    #[error("Variable not found: {0}")]
    VariableNotFound(String),
    #[error("Invalid variation: {0}")]
    InvalidVariation(String),
}

#[derive(Debug, Clone)]
pub struct SensitivityResult {
    pub variable_names: Vec<String>,
    /// df/dx_i normalized gradient per variable.
    pub gradients: Vec<f64>,
    /// Absolute sensitivity: |df/dx_i * delta_x_i|.
    pub absolute_sensitivities: Vec<f64>,
    pub base_value: f64,
    /// Per-variable delta used.
    pub deltas: Vec<f64>,
}

/// Perform one-at-a-time sensitivity analysis using central differences.
pub fn one_at_a_time<F>(
    variables: &[SensitivityVariable],
    base_values: &HashMap<String, f64>,
    eval_fn: F,
) -> Result<SensitivityResult, SensitivityError>
where
    F: Fn(&HashMap<String, f64>) -> Result<f64, SensitivityError>,
{
    let base_output = eval_fn(base_values)
        .map_err(|e| SensitivityError::EvalFailed(e.to_string()))?;

    let mut gradients = Vec::with_capacity(variables.len());
    let mut absolute_sensitivities = Vec::with_capacity(variables.len());
    let mut variable_names = Vec::with_capacity(variables.len());
    let mut deltas = Vec::with_capacity(variables.len());

    for var in variables {
        let base_val = base_values
            .get(&var.variable)
            .copied()
            .ok_or_else(|| SensitivityError::VariableNotFound(var.variable.clone()))?;

        let delta = parse_variation(&var.variation, base_val)?;

        // Central difference: (f(x+δ) - f(x-δ)) / 2δ
        let mut vals_plus = base_values.clone();
        vals_plus.insert(var.variable.clone(), base_val + delta);
        let f_plus = eval_fn(&vals_plus)
            .map_err(|e| SensitivityError::EvalFailed(e.to_string()))?;

        let mut vals_minus = base_values.clone();
        vals_minus.insert(var.variable.clone(), base_val - delta);
        let f_minus = eval_fn(&vals_minus)
            .map_err(|e| SensitivityError::EvalFailed(e.to_string()))?;

        let gradient = if delta.abs() > 1e-30 {
            (f_plus - f_minus) / (2.0 * delta)
        } else {
            0.0
        };

        variable_names.push(var.variable.clone());
        gradients.push(gradient);
        absolute_sensitivities.push((gradient * delta).abs());
        deltas.push(delta);
    }

    Ok(SensitivityResult {
        variable_names,
        gradients,
        absolute_sensitivities,
        base_value: base_output,
        deltas,
    })
}

/// Parse a variation string like "5%", "±10%", "0.5mm" into an absolute delta.
fn parse_variation(variation: &str, base_val: f64) -> Result<f64, SensitivityError> {
    let s = variation.trim().trim_start_matches('±').trim_start_matches('+').trim_start_matches('-');
    if s.ends_with('%') {
        let pct: f64 = s.trim_end_matches('%').parse().map_err(|_| {
            SensitivityError::InvalidVariation(format!("cannot parse '{}'", variation))
        })?;
        Ok(base_val.abs() * pct / 100.0)
    } else {
        // Absolute value (possibly with unit suffix — strip common suffixes)
        let numeric = s
            .trim_end_matches(|c: char| c.is_alphabetic())
            .parse::<f64>()
            .map_err(|_| {
                SensitivityError::InvalidVariation(format!("cannot parse '{}'", variation))
            })?;
        Ok(numeric)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadratic_sensitivity() {
        // f(x,y) = x^2 + 2*y → df/dx = 2x, df/dy = 2
        let variables = vec![
            SensitivityVariable {
                variable: "x".into(),
                variation: "5%".into(),
                distribution: "Uniform".into(),
            },
            SensitivityVariable {
                variable: "y".into(),
                variation: "5%".into(),
                distribution: "Uniform".into(),
            },
        ];

        let base = HashMap::from([("x".into(), 3.0), ("y".into(), 4.0)]);

        let result = one_at_a_time(&variables, &base, |vals| {
            let x = vals["x"];
            let y = vals["y"];
            Ok(x * x + 2.0 * y)
        })
        .unwrap();

        // df/dx at x=3 should be ~6.0
        assert!((result.gradients[0] - 6.0).abs() < 0.1, "got {}", result.gradients[0]);
        // df/dy should be ~2.0
        assert!((result.gradients[1] - 2.0).abs() < 0.1, "got {}", result.gradients[1]);
        assert!((result.base_value - 17.0).abs() < 1e-10); // 9 + 8
    }

    #[test]
    fn parse_percent_variation() {
        let delta = parse_variation("5%", 100.0).unwrap();
        assert!((delta - 5.0).abs() < 1e-10);

        let delta2 = parse_variation("±10%", 50.0).unwrap();
        assert!((delta2 - 5.0).abs() < 1e-10);
    }

    #[test]
    fn parse_absolute_variation() {
        let delta = parse_variation("0.5mm", 10.0).unwrap();
        assert!((delta - 0.5).abs() < 1e-10);
    }
}
