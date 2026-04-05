// ---------------------------------------------------------------------------
// Sweep Engine — Generate parameter combinations for parametric sweeps
// ---------------------------------------------------------------------------

use std::collections::HashMap;

use thiserror::Error;

use crate::optimetrics::SweepDefinition;

#[derive(Debug, Error)]
pub enum SweepError {
    #[error("Unknown sweep type: {0}")]
    UnknownType(String),
    #[error("Invalid parameter: {0}")]
    InvalidParam(String),
    #[error("Empty sweep definition")]
    Empty,
}

/// A fully expanded sweep plan with all variable combinations.
#[derive(Debug, Clone)]
pub struct SweepPlan {
    pub variable_names: Vec<String>,
    pub variations: Vec<HashMap<String, f64>>,
}

impl SweepPlan {
    pub fn num_variations(&self) -> usize {
        self.variations.len()
    }
}

/// Expand a single sweep definition into a list of values.
pub fn expand_sweep(def: &SweepDefinition) -> Result<Vec<f64>, SweepError> {
    match def.sweep_type.as_str() {
        "LinearStep" => {
            let start = get_param_f64(&def.params, "start")?;
            let stop = get_param_f64(&def.params, "stop")?;
            let step = get_param_f64(&def.params, "step")?;
            if step.abs() < 1e-30 {
                return Err(SweepError::InvalidParam("step cannot be zero".into()));
            }
            let mut values = Vec::new();
            let mut v = start;
            if step > 0.0 {
                while v <= stop + step * 1e-10 {
                    values.push(v);
                    v += step;
                }
            } else {
                while v >= stop + step * 1e-10 {
                    values.push(v);
                    v += step;
                }
            }
            Ok(values)
        }
        "LinearCount" => {
            let start = get_param_f64(&def.params, "start")?;
            let stop = get_param_f64(&def.params, "stop")?;
            let count = get_param_usize(&def.params, "count")?;
            if count == 0 {
                return Err(SweepError::InvalidParam("count must be > 0".into()));
            }
            if count == 1 {
                return Ok(vec![start]);
            }
            let step = (stop - start) / (count - 1) as f64;
            let values: Vec<f64> = (0..count).map(|i| start + i as f64 * step).collect();
            Ok(values)
        }
        "LogScale" => {
            let start = get_param_f64(&def.params, "start")?;
            let stop = get_param_f64(&def.params, "stop")?;
            let count = get_param_usize(&def.params, "count")?;
            if start <= 0.0 || stop <= 0.0 {
                return Err(SweepError::InvalidParam(
                    "LogScale requires positive start/stop".into(),
                ));
            }
            if count == 0 {
                return Err(SweepError::InvalidParam("count must be > 0".into()));
            }
            if count == 1 {
                return Ok(vec![start]);
            }
            let log_start = start.log10();
            let log_stop = stop.log10();
            let step = (log_stop - log_start) / (count - 1) as f64;
            let values: Vec<f64> = (0..count)
                .map(|i| 10.0f64.powf(log_start + i as f64 * step))
                .collect();
            Ok(values)
        }
        "DiscreteList" => {
            let list = def
                .params
                .as_array()
                .or_else(|| def.params.get("values").and_then(|v| v.as_array()))
                .ok_or_else(|| {
                    SweepError::InvalidParam("DiscreteList requires an array of values".into())
                })?;
            let values: Vec<f64> = list
                .iter()
                .filter_map(|v| v.as_f64())
                .collect();
            if values.is_empty() {
                return Err(SweepError::InvalidParam("empty discrete list".into()));
            }
            Ok(values)
        }
        other => Err(SweepError::UnknownType(other.into())),
    }
}

/// Generate the full Cartesian product of all sweep definitions.
pub fn generate_sweep_plan(defs: &[SweepDefinition]) -> Result<SweepPlan, SweepError> {
    if defs.is_empty() {
        return Err(SweepError::Empty);
    }

    let mut variable_names = Vec::new();
    let mut value_lists: Vec<Vec<f64>> = Vec::new();

    for def in defs {
        variable_names.push(def.variable.clone());
        value_lists.push(expand_sweep(def)?);
    }

    let variations = cartesian_product(&variable_names, &value_lists);

    Ok(SweepPlan {
        variable_names,
        variations,
    })
}

/// Generate a variation directory name: "variation_001", "variation_002", etc.
pub fn variation_dir_name(index: usize) -> String {
    format!("variation_{:03}", index)
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn get_param_f64(params: &serde_json::Value, key: &str) -> Result<f64, SweepError> {
    params
        .get(key)
        .and_then(|v| v.as_f64())
        .ok_or_else(|| SweepError::InvalidParam(format!("missing or invalid '{}'", key)))
}

fn get_param_usize(params: &serde_json::Value, key: &str) -> Result<usize, SweepError> {
    params
        .get(key)
        .and_then(|v| v.as_u64())
        .map(|v| v as usize)
        .ok_or_else(|| SweepError::InvalidParam(format!("missing or invalid '{}'", key)))
}

fn cartesian_product(
    names: &[String],
    value_lists: &[Vec<f64>],
) -> Vec<HashMap<String, f64>> {
    if names.is_empty() || value_lists.is_empty() {
        return Vec::new();
    }

    let total: usize = value_lists.iter().map(|v| v.len()).product();
    let mut result = Vec::with_capacity(total);

    let mut indices = vec![0usize; names.len()];
    for _ in 0..total {
        let mut map = HashMap::new();
        for (i, name) in names.iter().enumerate() {
            map.insert(name.clone(), value_lists[i][indices[i]]);
        }
        result.push(map);

        // Increment indices (rightmost first, like an odometer)
        let mut carry = true;
        for i in (0..names.len()).rev() {
            if carry {
                indices[i] += 1;
                if indices[i] >= value_lists[i].len() {
                    indices[i] = 0;
                } else {
                    carry = false;
                }
            }
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn sweep_def(variable: &str, sweep_type: &str, params: serde_json::Value) -> SweepDefinition {
        SweepDefinition {
            variable: variable.into(),
            sweep_type: sweep_type.into(),
            params,
        }
    }

    #[test]
    fn linear_step() {
        let def = sweep_def(
            "x",
            "LinearStep",
            serde_json::json!({"start": 0.0, "stop": 10.0, "step": 2.0}),
        );
        let values = expand_sweep(&def).unwrap();
        assert_eq!(values, vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
    }

    #[test]
    fn linear_count() {
        let def = sweep_def(
            "x",
            "LinearCount",
            serde_json::json!({"start": 0.0, "stop": 1.0, "count": 5}),
        );
        let values = expand_sweep(&def).unwrap();
        assert_eq!(values.len(), 5);
        assert!((values[0] - 0.0).abs() < 1e-10);
        assert!((values[4] - 1.0).abs() < 1e-10);
        assert!((values[2] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn log_scale() {
        let def = sweep_def(
            "freq",
            "LogScale",
            serde_json::json!({"start": 1.0, "stop": 1000.0, "count": 4}),
        );
        let values = expand_sweep(&def).unwrap();
        assert_eq!(values.len(), 4);
        assert!((values[0] - 1.0).abs() < 1e-10);
        assert!((values[1] - 10.0).abs() < 1e-6);
        assert!((values[2] - 100.0).abs() < 1e-4);
        assert!((values[3] - 1000.0).abs() < 1e-2);
    }

    #[test]
    fn discrete_list() {
        let def = sweep_def("mat", "DiscreteList", serde_json::json!([1.0, 2.5, 4.2]));
        let values = expand_sweep(&def).unwrap();
        assert_eq!(values, vec![1.0, 2.5, 4.2]);
    }

    #[test]
    fn cartesian_product_2_vars() {
        let defs = vec![
            sweep_def(
                "x",
                "LinearCount",
                serde_json::json!({"start": 1.0, "stop": 3.0, "count": 3}),
            ),
            sweep_def(
                "y",
                "LinearCount",
                serde_json::json!({"start": 10.0, "stop": 20.0, "count": 2}),
            ),
        ];
        let plan = generate_sweep_plan(&defs).unwrap();
        assert_eq!(plan.num_variations(), 6); // 3 × 2
        assert_eq!(plan.variable_names, vec!["x", "y"]);

        // First variation: x=1.0, y=10.0
        assert!((plan.variations[0]["x"] - 1.0).abs() < 1e-10);
        assert!((plan.variations[0]["y"] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn variation_dir_names() {
        assert_eq!(variation_dir_name(1), "variation_001");
        assert_eq!(variation_dir_name(42), "variation_042");
        assert_eq!(variation_dir_name(100), "variation_100");
    }

    #[test]
    fn unknown_sweep_type() {
        let def = sweep_def("x", "Invalid", serde_json::json!({}));
        assert!(expand_sweep(&def).is_err());
    }

    #[test]
    fn zero_step_error() {
        let def = sweep_def(
            "x",
            "LinearStep",
            serde_json::json!({"start": 0.0, "stop": 10.0, "step": 0.0}),
        );
        assert!(expand_sweep(&def).is_err());
    }

    #[test]
    fn single_value_count() {
        let def = sweep_def(
            "x",
            "LinearCount",
            serde_json::json!({"start": 5.0, "stop": 5.0, "count": 1}),
        );
        let values = expand_sweep(&def).unwrap();
        assert_eq!(values, vec![5.0]);
    }
}
