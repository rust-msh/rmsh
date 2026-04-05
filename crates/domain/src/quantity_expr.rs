// ---------------------------------------------------------------------------
// QuantityExpression — Parse and evaluate report trace expressions
// ---------------------------------------------------------------------------
//
// Supported expression syntax (HFSS/Q3D-style):
//   dB(S(1,1))       → S11 magnitude in dB
//   ang_deg(S(2,1))  → S21 phase in degrees
//   ang_rad(S(1,1))  → S11 phase in radians
//   mag(S(1,1))      → S11 linear magnitude
//   re(S(1,1))       → S11 real part
//   im(S(1,1))       → S11 imaginary part
//   VSWR(S(1,1))     → VSWR from S11
//   S(1,1)           → bare (re+j*im), used for Smith chart
//   GainTotal         → far-field gain
//   DeltaS            → convergence delta S
//   Tetrahedra        → convergence mesh count
//   R(1,1), L(1,1), C(1,1), G(1,1) → RLCG matrix elements

use std::f64::consts::PI;

use thiserror::Error;

use crate::result_store::{
    ConvergenceData, FarFieldData, RlcgMatrixData, SParameterData,
};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub enum BaseQuantity {
    /// S-parameter S(row, col), 1-based indexing as in UI expressions.
    SParameter { row: usize, col: usize },
    /// Z-parameter Z(row, col), 1-based.
    ZParameter { row: usize, col: usize },
    /// Far-field derived quantity by name (e.g. "GainTotal", "GainTheta").
    FarFieldQuantity(String),
    /// Convergence: max delta S per pass.
    ConvergenceDeltaS,
    /// Convergence: max delta energy per pass.
    ConvergenceDeltaEnergy,
    /// Convergence: number of tetrahedra per pass.
    ConvergenceTetrahedra,
    /// RLCG matrix element: matrix_type ('R','L','C','G'), row, col (0-based).
    RlcgElement {
        matrix_type: char,
        row: usize,
        col: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Transform {
    /// Linear magnitude: |z|
    Magnitude,
    /// Magnitude in dB: 20*log10(|z|)
    MagnitudeDB,
    /// Phase in degrees: arg(z) * 180/π
    AngleDeg,
    /// Phase in radians: arg(z)
    AngleRad,
    /// Real part
    Real,
    /// Imaginary part
    Imaginary,
    /// VSWR: (1+|Γ|)/(1-|Γ|)
    VSWR,
}

#[derive(Debug, Clone, PartialEq)]
pub struct QuantityExpression {
    pub base: BaseQuantity,
    pub transforms: Vec<Transform>,
}

#[derive(Debug, Error)]
pub enum ExprError {
    #[error("Parse error: {0}")]
    Parse(String),
    #[error("Evaluation error: {0}")]
    Eval(String),
    #[error("Data not available: {0}")]
    DataNotAvailable(String),
}

// ---------------------------------------------------------------------------
// Data source abstraction for evaluation
// ---------------------------------------------------------------------------

/// Evaluation context providing access to different result data sources.
pub enum EvalContext<'a> {
    SParameter {
        data: &'a SParameterData,
    },
    Convergence {
        data: &'a ConvergenceData,
    },
    FarField {
        data: &'a FarFieldData,
        phi_idx: usize,
    },
    Rlcg {
        data: &'a RlcgMatrixData,
    },
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

impl QuantityExpression {
    /// Parse a trace expression string into a QuantityExpression.
    pub fn parse(expr: &str) -> Result<Self, ExprError> {
        let expr = expr.trim();

        // Try convergence quantities first (simple names)
        match expr {
            "DeltaS" => {
                return Ok(Self {
                    base: BaseQuantity::ConvergenceDeltaS,
                    transforms: vec![],
                });
            }
            "DeltaEnergy" => {
                return Ok(Self {
                    base: BaseQuantity::ConvergenceDeltaEnergy,
                    transforms: vec![],
                });
            }
            "Tetrahedra" => {
                return Ok(Self {
                    base: BaseQuantity::ConvergenceTetrahedra,
                    transforms: vec![],
                });
            }
            _ => {}
        }

        // Try far-field named quantities (no parentheses)
        if matches!(
            expr,
            "GainTotal"
                | "GainTheta"
                | "GainPhi"
                | "DirectivityTotal"
                | "AxialRatio"
        ) {
            return Ok(Self {
                base: BaseQuantity::FarFieldQuantity(expr.to_string()),
                transforms: vec![],
            });
        }

        // Parse transform-wrapped expressions: transform(inner)
        if let Some((transform, inner)) = parse_outer_function(expr) {
            let t = match transform {
                "dB" => Transform::MagnitudeDB,
                "mag" => Transform::Magnitude,
                "ang_deg" => Transform::AngleDeg,
                "ang_rad" => Transform::AngleRad,
                "re" => Transform::Real,
                "im" => Transform::Imaginary,
                "VSWR" => Transform::VSWR,
                other => {
                    return Err(ExprError::Parse(format!(
                        "Unknown transform: '{other}'"
                    )));
                }
            };

            let mut inner_expr = Self::parse(inner)?;
            // Prepend transform (outermost is applied last, so we push to end)
            inner_expr.transforms.push(t);
            return Ok(inner_expr);
        }

        // Parse base quantities: S(row,col), Z(row,col), R(row,col), etc.
        if let Some((name, args)) = parse_function_call(expr) {
            match name {
                "S" => {
                    let (row, col) = parse_two_args(args)?;
                    return Ok(Self {
                        base: BaseQuantity::SParameter { row, col },
                        transforms: vec![],
                    });
                }
                "Z" => {
                    let (row, col) = parse_two_args(args)?;
                    return Ok(Self {
                        base: BaseQuantity::ZParameter { row, col },
                        transforms: vec![],
                    });
                }
                "R" | "L" | "C" | "G" => {
                    let (row, col) = parse_two_args(args)?;
                    let matrix_type = name.chars().next().unwrap();
                    return Ok(Self {
                        base: BaseQuantity::RlcgElement {
                            matrix_type,
                            row: row - 1, // UI uses 1-based, store 0-based
                            col: col - 1,
                        },
                        transforms: vec![],
                    });
                }
                _ => {
                    return Err(ExprError::Parse(format!(
                        "Unknown function: '{name}'"
                    )));
                }
            }
        }

        Err(ExprError::Parse(format!(
            "Cannot parse expression: '{expr}'"
        )))
    }

    /// Evaluate the expression against a data context.
    /// Returns Vec<[x, y]> where x is the sweep variable (frequency, pass number, etc.).
    pub fn evaluate(&self, ctx: &EvalContext<'_>) -> Result<Vec<[f64; 2]>, ExprError> {
        match ctx {
            EvalContext::SParameter { data } => self.eval_s_parameter(data),
            EvalContext::Convergence { data } => self.eval_convergence(data),
            EvalContext::FarField { data, phi_idx } => self.eval_far_field(data, *phi_idx),
            EvalContext::Rlcg { data } => self.eval_rlcg(data),
        }
    }

    fn eval_s_parameter(&self, data: &SParameterData) -> Result<Vec<[f64; 2]>, ExprError> {
        let (row, col) = match &self.base {
            BaseQuantity::SParameter { row, col } => (*row, *col),
            other => {
                return Err(ExprError::Eval(format!(
                    "Cannot evaluate {other:?} with S-parameter context"
                )));
            }
        };

        let key = SParameterData::s_key(row, col);
        let complex_data = data.get_complex(&key).ok_or_else(|| {
            ExprError::DataNotAvailable(format!("S-parameter '{key}' not found"))
        })?;

        let freqs = data.frequencies_original();
        if freqs.len() != complex_data.len() {
            return Err(ExprError::Eval(
                "Frequency and parameter data length mismatch".into(),
            ));
        }

        let result: Vec<[f64; 2]> = freqs
            .iter()
            .zip(complex_data.iter())
            .map(|(&freq, &[re, im])| {
                let y = apply_transforms(re, im, &self.transforms);
                [freq, y]
            })
            .collect();

        Ok(result)
    }

    fn eval_convergence(&self, data: &ConvergenceData) -> Result<Vec<[f64; 2]>, ExprError> {
        match &self.base {
            BaseQuantity::ConvergenceDeltaS => Ok(data.delta_s_curve()),
            BaseQuantity::ConvergenceDeltaEnergy => Ok(data.delta_energy_curve()),
            BaseQuantity::ConvergenceTetrahedra => Ok(data.tetrahedra_curve()),
            other => Err(ExprError::Eval(format!(
                "Cannot evaluate {other:?} with convergence context"
            ))),
        }
    }

    fn eval_far_field(
        &self,
        data: &FarFieldData,
        phi_idx: usize,
    ) -> Result<Vec<[f64; 2]>, ExprError> {
        match &self.base {
            BaseQuantity::FarFieldQuantity(name) => {
                data.theta_cut(name, phi_idx).ok_or_else(|| {
                    ExprError::DataNotAvailable(format!(
                        "Far-field quantity '{name}' not found"
                    ))
                })
            }
            other => Err(ExprError::Eval(format!(
                "Cannot evaluate {other:?} with far-field context"
            ))),
        }
    }

    fn eval_rlcg(&self, data: &RlcgMatrixData) -> Result<Vec<[f64; 2]>, ExprError> {
        match &self.base {
            BaseQuantity::RlcgElement {
                matrix_type,
                row,
                col,
            } => {
                let mt = matrix_type.to_string();
                data.element_curve(&mt, *row, *col).ok_or_else(|| {
                    ExprError::DataNotAvailable(format!(
                        "RLCG element {mt}({},{}) not found",
                        row + 1,
                        col + 1,
                    ))
                })
            }
            other => Err(ExprError::Eval(format!(
                "Cannot evaluate {other:?} with RLCG context"
            ))),
        }
    }
}

// ---------------------------------------------------------------------------
// Transform application
// ---------------------------------------------------------------------------

fn apply_transforms(re: f64, im: f64, transforms: &[Transform]) -> f64 {
    if transforms.is_empty() {
        // No transforms: return magnitude by default for numeric display
        return (re * re + im * im).sqrt();
    }

    // Apply transforms in order (innermost first, outermost last).
    // The last transform in the list is the outermost one.
    // For a single complex value, we apply them sequentially.
    let mut val_re = re;
    let mut val_im = im;
    let mut is_real = false;

    for t in transforms {
        if is_real {
            // After a transform produces a real scalar, further transforms
            // operate on (val_re, 0.0).
            val_im = 0.0;
        }
        match t {
            Transform::MagnitudeDB => {
                let mag = (val_re * val_re + val_im * val_im).sqrt();
                val_re = 20.0 * mag.log10();
                val_im = 0.0;
                is_real = true;
            }
            Transform::Magnitude => {
                val_re = (val_re * val_re + val_im * val_im).sqrt();
                val_im = 0.0;
                is_real = true;
            }
            Transform::AngleDeg => {
                val_re = val_im.atan2(val_re) * 180.0 / PI;
                val_im = 0.0;
                is_real = true;
            }
            Transform::AngleRad => {
                val_re = val_im.atan2(val_re);
                val_im = 0.0;
                is_real = true;
            }
            Transform::Real => {
                // val_re stays as is
                val_im = 0.0;
                is_real = true;
            }
            Transform::Imaginary => {
                val_re = val_im;
                val_im = 0.0;
                is_real = true;
            }
            Transform::VSWR => {
                let gamma = (val_re * val_re + val_im * val_im).sqrt();
                val_re = if gamma < 1.0 {
                    (1.0 + gamma) / (1.0 - gamma)
                } else {
                    f64::INFINITY
                };
                val_im = 0.0;
                is_real = true;
            }
        }
    }

    val_re
}

// ---------------------------------------------------------------------------
// Parsing helpers
// ---------------------------------------------------------------------------

/// Parse "funcname(inner)" → Some(("funcname", "inner")).
fn parse_outer_function(expr: &str) -> Option<(&str, &str)> {
    let open = expr.find('(')?;
    if !expr.ends_with(')') {
        return None;
    }
    let name = &expr[..open];

    // Only match known transform names to avoid capturing base functions like S(1,1)
    if !matches!(
        name,
        "dB" | "mag" | "ang_deg" | "ang_rad" | "re" | "im" | "VSWR"
    ) {
        return None;
    }

    let inner = &expr[open + 1..expr.len() - 1];
    Some((name, inner))
}

/// Parse "funcname(args)" → Some(("funcname", "args")).
fn parse_function_call(expr: &str) -> Option<(&str, &str)> {
    let open = expr.find('(')?;
    if !expr.ends_with(')') {
        return None;
    }
    let name = &expr[..open];
    let args = &expr[open + 1..expr.len() - 1];
    Some((name, args))
}

/// Parse "row,col" → (row, col) as usize.
fn parse_two_args(args: &str) -> Result<(usize, usize), ExprError> {
    let parts: Vec<&str> = args.split(',').map(|s| s.trim()).collect();
    if parts.len() != 2 {
        return Err(ExprError::Parse(format!(
            "Expected 2 arguments, got {}",
            parts.len()
        )));
    }
    let row: usize = parts[0]
        .parse()
        .map_err(|_| ExprError::Parse(format!("Invalid row: '{}'", parts[0])))?;
    let col: usize = parts[1]
        .parse()
        .map_err(|_| ExprError::Parse(format!("Invalid col: '{}'", parts[1])))?;
    Ok((row, col))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_db_s11() {
        let expr = QuantityExpression::parse("dB(S(1,1))").unwrap();
        assert_eq!(expr.base, BaseQuantity::SParameter { row: 1, col: 1 });
        assert_eq!(expr.transforms, vec![Transform::MagnitudeDB]);
    }

    #[test]
    fn parse_ang_deg_s21() {
        let expr = QuantityExpression::parse("ang_deg(S(2,1))").unwrap();
        assert_eq!(expr.base, BaseQuantity::SParameter { row: 2, col: 1 });
        assert_eq!(expr.transforms, vec![Transform::AngleDeg]);
    }

    #[test]
    fn parse_bare_s_param() {
        let expr = QuantityExpression::parse("S(1,1)").unwrap();
        assert_eq!(expr.base, BaseQuantity::SParameter { row: 1, col: 1 });
        assert!(expr.transforms.is_empty());
    }

    #[test]
    fn parse_vswr() {
        let expr = QuantityExpression::parse("VSWR(S(1,1))").unwrap();
        assert_eq!(expr.base, BaseQuantity::SParameter { row: 1, col: 1 });
        assert_eq!(expr.transforms, vec![Transform::VSWR]);
    }

    #[test]
    fn parse_real_imaginary() {
        let re = QuantityExpression::parse("re(S(1,1))").unwrap();
        assert_eq!(re.transforms, vec![Transform::Real]);

        let im = QuantityExpression::parse("im(S(1,1))").unwrap();
        assert_eq!(im.transforms, vec![Transform::Imaginary]);
    }

    #[test]
    fn parse_convergence_quantities() {
        let ds = QuantityExpression::parse("DeltaS").unwrap();
        assert_eq!(ds.base, BaseQuantity::ConvergenceDeltaS);

        let tet = QuantityExpression::parse("Tetrahedra").unwrap();
        assert_eq!(tet.base, BaseQuantity::ConvergenceTetrahedra);
    }

    #[test]
    fn parse_far_field() {
        let expr = QuantityExpression::parse("GainTotal").unwrap();
        assert_eq!(
            expr.base,
            BaseQuantity::FarFieldQuantity("GainTotal".into())
        );
    }

    #[test]
    fn parse_rlcg() {
        let expr = QuantityExpression::parse("R(1,2)").unwrap();
        assert_eq!(
            expr.base,
            BaseQuantity::RlcgElement {
                matrix_type: 'R',
                row: 0,
                col: 1,
            }
        );
    }

    #[test]
    fn parse_error_unknown_function() {
        let result = QuantityExpression::parse("foo(S(1,1))");
        assert!(result.is_err());
    }

    #[test]
    fn transform_db() {
        // |0.5 + 0j| = 0.5; dB = 20*log10(0.5) ≈ -6.02
        let val = apply_transforms(0.5, 0.0, &[Transform::MagnitudeDB]);
        assert!((val - (-6.0206)).abs() < 0.001);
    }

    #[test]
    fn transform_phase() {
        // 0 + 1j → phase = 90°
        let val = apply_transforms(0.0, 1.0, &[Transform::AngleDeg]);
        assert!((val - 90.0).abs() < 0.001);
    }

    #[test]
    fn transform_vswr() {
        // |0.5 + 0j| = 0.5; VSWR = (1+0.5)/(1-0.5) = 3.0
        let val = apply_transforms(0.5, 0.0, &[Transform::VSWR]);
        assert!((val - 3.0).abs() < 0.001);
    }

    #[test]
    fn eval_s_parameter_db() {
        let data = SParameterData {
            format_version: "1.0".into(),
            file_type: "SParameterData".into(),
            design_id: "d1".into(),
            setup: "Setup1".into(),
            sweep: String::new(),
            solution_type: "DrivenModal".into(),
            reference_impedance_ohm: 50.0,
            num_ports: 1,
            port_names: vec!["Port1".into()],
            num_frequencies: 2,
            frequency_unit: "GHz".into(),
            data_format: "RealImaginary".into(),
            frequencies: vec![1.0, 2.0],
            parameters: {
                let mut m = std::collections::HashMap::new();
                m.insert(
                    "S11".into(),
                    crate::result_store::SParamValues {
                        real: vec![-0.5, -0.3],
                        imag: vec![0.0, 0.0],
                    },
                );
                m
            },
            derived: std::collections::HashMap::new(),
        };

        let expr = QuantityExpression::parse("dB(S(1,1))").unwrap();
        let ctx = EvalContext::SParameter { data: &data };
        let result = expr.evaluate(&ctx).unwrap();

        assert_eq!(result.len(), 2);
        assert!((result[0][0] - 1.0).abs() < 1e-10);
        // dB(0.5) ≈ -6.02
        assert!((result[0][1] - (-6.0206)).abs() < 0.01);
    }
}
