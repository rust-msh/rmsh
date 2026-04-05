// ---------------------------------------------------------------------------
// Statistical Analysis — Monte Carlo simulation
// ---------------------------------------------------------------------------

use std::collections::HashMap;

use thiserror::Error;

use crate::optimetrics::SensitivityVariable;
use crate::optimization::SimpleRng;

#[derive(Debug, Error)]
pub enum StatError {
    #[error("Evaluation failed: {0}")]
    EvalFailed(String),
    #[error("Variable not found: {0}")]
    VariableNotFound(String),
    #[error("Invalid distribution: {0}")]
    InvalidDistribution(String),
}

#[derive(Debug, Clone)]
pub struct MonteCarloResult {
    pub variable_names: Vec<String>,
    /// [trial][variable] — sampled variable values per trial.
    pub samples: Vec<Vec<f64>>,
    /// Output value per trial.
    pub outputs: Vec<f64>,
    pub mean: f64,
    pub std_dev: f64,
    pub min: f64,
    pub max: f64,
    /// Percentiles: "p5", "p25", "p50", "p75", "p95".
    pub percentiles: HashMap<String, f64>,
}

/// Run Monte Carlo analysis with random sampling.
pub fn run_monte_carlo<F>(
    variables: &[SensitivityVariable],
    base_values: &HashMap<String, f64>,
    num_trials: u32,
    eval_fn: F,
    seed: u64,
) -> Result<MonteCarloResult, StatError>
where
    F: Fn(&HashMap<String, f64>) -> Result<f64, StatError>,
{
    let mut rng = SimpleRng::new(seed);

    let variable_names: Vec<String> = variables.iter().map(|v| v.variable.clone()).collect();
    let mut samples = Vec::with_capacity(num_trials as usize);
    let mut outputs = Vec::with_capacity(num_trials as usize);

    for _ in 0..num_trials {
        let mut trial_values = base_values.clone();
        let mut sample_row = Vec::with_capacity(variables.len());

        for var in variables {
            let base_val = base_values
                .get(&var.variable)
                .copied()
                .ok_or_else(|| StatError::VariableNotFound(var.variable.clone()))?;

            let delta = parse_variation_range(&var.variation, base_val);
            let sampled = sample_distribution(&var.distribution, base_val, delta, &mut rng)?;

            trial_values.insert(var.variable.clone(), sampled);
            sample_row.push(sampled);
        }

        let output = eval_fn(&trial_values)
            .map_err(|e| StatError::EvalFailed(e.to_string()))?;

        samples.push(sample_row);
        outputs.push(output);
    }

    // Compute statistics
    let n = outputs.len() as f64;
    let mean = outputs.iter().sum::<f64>() / n;
    let variance = outputs.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
    let std_dev = variance.sqrt();
    let min = outputs.iter().fold(f64::MAX, |a, &b| a.min(b));
    let max = outputs.iter().fold(f64::MIN, |a, &b| a.max(b));

    // Percentiles
    let mut sorted = outputs.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let percentiles = HashMap::from([
        ("p5".into(), percentile(&sorted, 0.05)),
        ("p25".into(), percentile(&sorted, 0.25)),
        ("p50".into(), percentile(&sorted, 0.50)),
        ("p75".into(), percentile(&sorted, 0.75)),
        ("p95".into(), percentile(&sorted, 0.95)),
    ]);

    Ok(MonteCarloResult {
        variable_names,
        samples,
        outputs,
        mean,
        std_dev,
        min,
        max,
        percentiles,
    })
}

fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = (p * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

fn parse_variation_range(variation: &str, base_val: f64) -> f64 {
    let s = variation.trim().trim_start_matches('±').trim_start_matches('+').trim_start_matches('-');
    if s.ends_with('%') {
        if let Ok(pct) = s.trim_end_matches('%').parse::<f64>() {
            return base_val.abs() * pct / 100.0;
        }
    }
    let numeric = s.trim_end_matches(|c: char| c.is_alphabetic());
    numeric.parse::<f64>().unwrap_or(base_val * 0.1)
}

fn sample_distribution(
    distribution: &str,
    center: f64,
    half_range: f64,
    rng: &mut SimpleRng,
) -> Result<f64, StatError> {
    match distribution {
        "Uniform" => {
            let v = center - half_range + rng.next_f64() * 2.0 * half_range;
            Ok(v)
        }
        "Gaussian" | "Normal" => {
            let v = rng.next_normal(center, half_range / 3.0); // 3σ = half_range
            Ok(v)
        }
        "LogNormal" => {
            // LogNormal: log(X) ~ N(mu, sigma)
            let sigma = (1.0 + (half_range / center).powi(2)).ln().sqrt();
            let mu = center.ln() - sigma * sigma / 2.0;
            let z = rng.next_normal(mu, sigma);
            Ok(z.exp())
        }
        other => Err(StatError::InvalidDistribution(other.into())),
    }
}

/// Build histogram bins from output data.
pub fn build_histogram(data: &[f64], num_bins: usize) -> (Vec<f64>, Vec<usize>) {
    if data.is_empty() || num_bins == 0 {
        return (Vec::new(), Vec::new());
    }

    let min = data.iter().fold(f64::MAX, |a, &b| a.min(b));
    let max = data.iter().fold(f64::MIN, |a, &b| a.max(b));
    let range = max - min;
    let bin_width = if range > 1e-30 { range / num_bins as f64 } else { 1.0 };

    let bin_centers: Vec<f64> = (0..num_bins)
        .map(|i| min + (i as f64 + 0.5) * bin_width)
        .collect();

    let mut counts = vec![0usize; num_bins];
    for &v in data {
        let idx = ((v - min) / bin_width) as usize;
        let idx = idx.min(num_bins - 1);
        counts[idx] += 1;
    }

    (bin_centers, counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn monte_carlo_identity() {
        let variables = vec![SensitivityVariable {
            variable: "x".into(),
            variation: "10%".into(),
            distribution: "Uniform".into(),
        }];

        let base = HashMap::from([("x".into(), 100.0)]);

        let result = run_monte_carlo(
            &variables,
            &base,
            1000,
            |vals| Ok(vals["x"]),
            42,
        )
        .unwrap();

        // Mean should be close to 100.0
        assert!(
            (result.mean - 100.0).abs() < 5.0,
            "mean={}",
            result.mean
        );
        // Std dev should be > 0 (there's variance)
        assert!(result.std_dev > 0.0);
        // Outputs should span the range
        assert!(result.min < 100.0);
        assert!(result.max > 100.0);
        assert_eq!(result.outputs.len(), 1000);
    }

    #[test]
    fn monte_carlo_gaussian() {
        let variables = vec![SensitivityVariable {
            variable: "x".into(),
            variation: "5%".into(),
            distribution: "Gaussian".into(),
        }];

        let base = HashMap::from([("x".into(), 50.0)]);

        let result = run_monte_carlo(&variables, &base, 500, |vals| Ok(vals["x"]), 123).unwrap();
        assert!((result.mean - 50.0).abs() < 3.0, "mean={}", result.mean);
    }

    #[test]
    fn histogram_basic() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let (centers, counts) = build_histogram(&data, 5);
        assert_eq!(centers.len(), 5);
        assert_eq!(counts.len(), 5);
        let total: usize = counts.iter().sum();
        assert_eq!(total, 10);
    }

    #[test]
    fn percentiles() {
        let sorted: Vec<f64> = (0..100).map(|i| i as f64).collect();
        assert!((percentile(&sorted, 0.5) - 50.0).abs() < 1.0);
        assert!((percentile(&sorted, 0.0) - 0.0).abs() < 1.0);
        assert!((percentile(&sorted, 1.0) - 99.0).abs() < 1.0);
    }
}
