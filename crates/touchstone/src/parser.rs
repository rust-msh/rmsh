use crate::error::TouchstoneError;
use crate::types::*;

/// Parse a Touchstone file from a string.
///
/// Supports Touchstone v1.0 and v2.0 formats.
/// All parameter data is converted to Real/Imaginary representation internally.
pub fn parse(input: &str) -> Result<Touchstone, TouchstoneError> {
    let mut comments = Vec::new();
    let mut option_line: Option<(usize, &str)> = None;
    let mut version: Option<String> = None;
    let mut num_ports_v2: Option<usize> = None;
    let mut two_port_order = None;
    let mut reference_impedances: Option<Vec<f64>> = None;
    let mut data_lines: Vec<(usize, &str)> = Vec::new();
    let mut in_network_data = false;
    let mut past_network_data = false;
    let mut in_reference = false;
    let mut ref_values = Vec::new();

    // First pass: classify lines
    for (line_idx, raw_line) in input.lines().enumerate() {
        let line_num = line_idx + 1;
        let line = raw_line.trim();

        if line.is_empty() {
            continue;
        }

        // Comments
        if line.starts_with('!') {
            comments.push(line[1..].trim().to_string());
            continue;
        }

        // V2.0 keywords
        if line.starts_with('[') {
            if let Some(keyword) = parse_keyword(line) {
                match keyword.0.to_lowercase().as_str() {
                    "version" => {
                        version = Some(keyword.1.to_string());
                    }
                    "number of ports" => {
                        num_ports_v2 = keyword
                            .1
                            .parse::<usize>()
                            .ok();
                    }
                    "two-port data order" => {
                        two_port_order = match keyword.1.to_lowercase().as_str() {
                            "21_12" => Some(TwoPortOrder::Order21_12),
                            "12_21" => Some(TwoPortOrder::Order12_21),
                            _ => None,
                        };
                    }
                    "reference" => {
                        in_reference = true;
                    }
                    "network data" => {
                        in_network_data = true;
                        in_reference = false;
                    }
                    "end" => {
                        in_network_data = false;
                        past_network_data = true;
                    }
                    _ => {}
                }
                continue;
            }
        }

        // Option line (check before reference block, since '#' terminates [Reference])
        if line.starts_with('#') {
            in_reference = false;
            option_line = Some((line_num, raw_line.trim()));
            continue;
        }

        // Reference impedance values (v2.0)
        if in_reference {
            for token in line.split_whitespace() {
                if let Ok(v) = token.parse::<f64>() {
                    ref_values.push(v);
                }
            }
            continue;
        }

        // Data lines: skip anything after [End] in v2
        if past_network_data {
            continue;
        }

        // For v2 with [Network Data], only collect lines inside that block.
        // For v1, collect all non-keyword non-comment non-option lines.
        if version.is_some() && !in_network_data {
            continue;
        }

        data_lines.push((line_num, raw_line.trim()));
    }

    if !ref_values.is_empty() {
        reference_impedances = Some(ref_values);
    }

    // Parse option line
    let options = match option_line {
        Some((line_num, line)) => parse_option_line(line_num, line)?,
        None => return Err(TouchstoneError::NoOptionLine),
    };

    if data_lines.is_empty() {
        return Err(TouchstoneError::NoData);
    }

    // Collect all numeric values from data lines (handle inline ! comments)
    let mut all_values: Vec<(usize, Vec<f64>)> = Vec::new();
    for (line_num, line) in &data_lines {
        // Strip inline comments
        let data_part = if let Some(idx) = line.find('!') {
            &line[..idx]
        } else {
            line
        };
        let tokens: Vec<&str> = data_part.split_whitespace().collect();
        if tokens.is_empty() {
            continue;
        }
        let mut values = Vec::with_capacity(tokens.len());
        for token in &tokens {
            let v = token.parse::<f64>().map_err(|_| TouchstoneError::InvalidNumber {
                line: *line_num,
                value: token.to_string(),
            })?;
            values.push(v);
        }
        all_values.push((*line_num, values));
    }

    // Determine port count
    let num_ports = determine_num_ports(&all_values, num_ports_v2)?;

    // Parse frequency points
    let data = parse_data_points(num_ports, &options, &all_values, two_port_order)?;

    // Validate against v2 port count
    if let Some(np) = num_ports_v2 {
        if np != num_ports {
            return Err(TouchstoneError::InconsistentPortCount {
                expected: np,
                got: num_ports,
            });
        }
    }

    Ok(Touchstone {
        num_ports,
        options,
        comments,
        data,
        version,
        two_port_order,
        reference_impedances,
    })
}

/// Parse a v2.0 keyword line like `[Number of Ports] 2`
fn parse_keyword(line: &str) -> Option<(&str, &str)> {
    let end = line.find(']')?;
    let key = &line[1..end];
    let value = line[end + 1..].trim();
    Some((key, value))
}

/// Parse the option line: `# GHz S RI R 50`
fn parse_option_line(line_num: usize, line: &str) -> Result<TouchstoneOptions, TouchstoneError> {
    // Strip the leading '#' and split
    let content = line[1..].trim();

    // Defaults per spec
    let mut freq_unit = FrequencyUnit::GHz;
    let mut param_type = ParameterType::S;
    let mut data_format = DataFormat::MagnitudeAngle;
    let mut ref_impedance = 50.0;

    let tokens: Vec<&str> = content.split_whitespace().collect();
    let mut i = 0;
    while i < tokens.len() {
        let token = tokens[i].to_uppercase();
        match token.as_str() {
            "HZ" => freq_unit = FrequencyUnit::Hz,
            "KHZ" => freq_unit = FrequencyUnit::KHz,
            "MHZ" => freq_unit = FrequencyUnit::MHz,
            "GHZ" => freq_unit = FrequencyUnit::GHz,
            "S" => param_type = ParameterType::S,
            "Y" => param_type = ParameterType::Y,
            "Z" => param_type = ParameterType::Z,
            "H" => param_type = ParameterType::H,
            "G" => param_type = ParameterType::G,
            "RI" => data_format = DataFormat::RealImaginary,
            "MA" => data_format = DataFormat::MagnitudeAngle,
            "DB" => data_format = DataFormat::DecibelAngle,
            "R" => {
                // Next token is the reference impedance
                i += 1;
                if i < tokens.len() {
                    ref_impedance = tokens[i].parse::<f64>().map_err(|_| {
                        TouchstoneError::ParseError {
                            line: line_num,
                            message: format!(
                                "invalid reference impedance: {}",
                                tokens[i]
                            ),
                        }
                    })?;
                }
            }
            _ => {
                // Unknown token — try to be lenient, skip it
            }
        }
        i += 1;
    }

    Ok(TouchstoneOptions {
        frequency_unit: freq_unit,
        parameter_type: param_type,
        data_format: data_format,
        reference_impedance: ref_impedance,
    })
}

/// Determine port count from data. For 1-port: 3 values/freq. For 2-port: 9. For N-port: 1 + 2*N*N.
fn determine_num_ports(
    all_values: &[(usize, Vec<f64>)],
    v2_ports: Option<usize>,
) -> Result<usize, TouchstoneError> {
    if let Some(np) = v2_ports {
        return Ok(np);
    }

    // Use the first data line's value count to guess port count
    let first_count = all_values
        .first()
        .map(|(_, v)| v.len())
        .ok_or(TouchstoneError::NoData)?;

    // 1-port: freq + 2 values = 3
    // 2-port: freq + 8 values = 9
    // 3-port: freq + 18 values = 19 (but often split across multiple lines)
    // 4-port: freq + 32 values = 33 (always multi-line)
    match first_count {
        3 => Ok(1),
        9 => Ok(2),
        _ => {
            // For N >= 3, first line of each freq point has: freq + 2*N values (first row)
            // Subsequent N-1 lines have: 2*N values each.
            // So first line has odd count (freq + even), the num_params = first_count - 1
            // num_params = 2*N (first row only), so N = (first_count - 1) / 2
            let n = (first_count - 1) / 2;
            if n >= 3 && (first_count - 1) == 2 * n {
                Ok(n)
            } else {
                // Fallback: assume 1-port if unclear
                Err(TouchstoneError::ParseError {
                    line: all_values.first().map(|(l, _)| *l).unwrap_or(0),
                    message: format!(
                        "cannot determine port count from {} values on first data line",
                        first_count
                    ),
                })
            }
        }
    }
}

/// Convert a pair of raw values to a Complex number based on the source format.
fn to_complex(v1: f64, v2: f64, format: DataFormat) -> Complex {
    match format {
        DataFormat::RealImaginary => Complex::new(v1, v2),
        DataFormat::MagnitudeAngle => Complex::from_mag_angle(v1, v2),
        DataFormat::DecibelAngle => Complex::from_db_angle(v1, v2),
    }
}

/// Parse all data points from collected numeric values.
fn parse_data_points(
    num_ports: usize,
    options: &TouchstoneOptions,
    all_values: &[(usize, Vec<f64>)],
    two_port_order: Option<TwoPortOrder>,
) -> Result<Vec<FrequencyPoint>, TouchstoneError> {
    let n = num_ports;
    let num_params = n * n; // Total S-parameters per frequency point
    let values_per_freq = 1 + 2 * num_params; // freq + pairs of (v1, v2)

    let mut result = Vec::new();

    if n <= 2 {
        // For 1-port and 2-port: all values on one line per frequency
        for (line_num, vals) in all_values {
            if vals.len() != values_per_freq {
                return Err(TouchstoneError::InvalidDataLine {
                    line: *line_num,
                    expected: values_per_freq,
                    got: vals.len(),
                });
            }
            let freq = vals[0];
            let params = build_matrix(n, &vals[1..], options.data_format, two_port_order);
            result.push(FrequencyPoint { frequency: freq, params });
        }
    } else {
        // For N >= 3: multi-line format
        // First line: freq + 2*N values (first row of matrix)
        // Subsequent lines: 2*N values each (remaining rows)
        // May be further split if 2*N > some threshold, but typically each row is one line
        let mut idx = 0;
        while idx < all_values.len() {
            let (line_num, first_vals) = &all_values[idx];

            // Gather all values for this frequency point
            let mut flat_values = first_vals.clone();
            idx += 1;

            // Collect continuation lines until we have enough values
            while flat_values.len() < values_per_freq && idx < all_values.len() {
                flat_values.extend_from_slice(&all_values[idx].1);
                idx += 1;
            }

            if flat_values.len() < values_per_freq {
                return Err(TouchstoneError::InvalidDataLine {
                    line: *line_num,
                    expected: values_per_freq,
                    got: flat_values.len(),
                });
            }

            let freq = flat_values[0];
            let params = build_matrix(n, &flat_values[1..], options.data_format, None);
            result.push(FrequencyPoint { frequency: freq, params });
        }
    }

    Ok(result)
}

/// Build the NxN parameter matrix from flat value pairs.
fn build_matrix(
    n: usize,
    pair_values: &[f64],
    format: DataFormat,
    two_port_order: Option<TwoPortOrder>,
) -> Vec<Vec<Complex>> {
    let mut matrix = vec![vec![Complex::zero(); n]; n];

    // Read pairs in order
    let mut pairs: Vec<Complex> = Vec::with_capacity(n * n);
    for i in 0..n * n {
        let v1 = pair_values[2 * i];
        let v2 = pair_values[2 * i + 1];
        pairs.push(to_complex(v1, v2, format));
    }

    if n == 2 {
        // 2-port Touchstone v1 default order: S11, S21, S12, S22
        // 2-port Touchstone v2 with 21_12: S11, S21, S12, S22 (same)
        // 2-port Touchstone v2 with 12_21: S11, S12, S21, S22
        let order = two_port_order.unwrap_or(TwoPortOrder::Order21_12);
        match order {
            TwoPortOrder::Order21_12 => {
                matrix[0][0] = pairs[0]; // S11
                matrix[1][0] = pairs[1]; // S21
                matrix[0][1] = pairs[2]; // S12
                matrix[1][1] = pairs[3]; // S22
            }
            TwoPortOrder::Order12_21 => {
                matrix[0][0] = pairs[0]; // S11
                matrix[0][1] = pairs[1]; // S12
                matrix[1][0] = pairs[2]; // S21
                matrix[1][1] = pairs[3]; // S22
            }
        }
    } else {
        // N-port: row-major order — S11, S12, ..., S1N, S21, S22, ..., SNN
        for row in 0..n {
            for col in 0..n {
                matrix[row][col] = pairs[row * n + col];
            }
        }
    }

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_1port_ri() {
        let input = "\
! 1-port S-parameter data
# MHz S RI R 50
100  0.9  -0.1
200  0.8  -0.2
500  0.5  -0.5
";
        let ts = parse(input).unwrap();
        assert_eq!(ts.num_ports, 1);
        assert_eq!(ts.options.frequency_unit, FrequencyUnit::MHz);
        assert_eq!(ts.options.parameter_type, ParameterType::S);
        assert_eq!(ts.options.data_format, DataFormat::RealImaginary);
        assert_eq!(ts.options.reference_impedance, 50.0);
        assert_eq!(ts.data.len(), 3);
        assert_eq!(ts.comments.len(), 1);

        let p0 = &ts.data[0];
        assert_eq!(p0.frequency, 100.0);
        assert!((p0.params[0][0].re - 0.9).abs() < 1e-10);
        assert!((p0.params[0][0].im - (-0.1)).abs() < 1e-10);
    }

    #[test]
    fn parse_1port_ma() {
        let input = "\
# GHz S MA R 50
2.4  0.5  -30.0
";
        let ts = parse(input).unwrap();
        let c = ts.data[0].params[0][0];
        let expected = Complex::from_mag_angle(0.5, -30.0);
        assert!((c.re - expected.re).abs() < 1e-10);
        assert!((c.im - expected.im).abs() < 1e-10);
    }

    #[test]
    fn parse_1port_db() {
        let input = "\
# GHz S DB R 50
2.4  -6.0  -30.0
";
        let ts = parse(input).unwrap();
        let c = ts.data[0].params[0][0];
        let expected = Complex::from_db_angle(-6.0, -30.0);
        assert!((c.re - expected.re).abs() < 1e-10);
        assert!((c.im - expected.im).abs() < 1e-10);
    }

    #[test]
    fn parse_2port_ri() {
        let input = "\
! 2-port network
# GHz S RI R 50
1.0  0.9 -0.1  0.05 0.01  0.05 0.01  0.85 -0.12
2.0  0.8 -0.2  0.10 0.02  0.10 0.02  0.78 -0.25
";
        let ts = parse(input).unwrap();
        assert_eq!(ts.num_ports, 2);
        assert_eq!(ts.data.len(), 2);

        // S11 of first point
        let s11 = ts.data[0].params[0][0];
        assert!((s11.re - 0.9).abs() < 1e-10);
        assert!((s11.im - (-0.1)).abs() < 1e-10);

        // S21 of first point (second value pair in v1 order)
        let s21 = ts.data[0].params[1][0];
        assert!((s21.re - 0.05).abs() < 1e-10);
        assert!((s21.im - 0.01).abs() < 1e-10);

        // S12
        let s12 = ts.data[0].params[0][1];
        assert!((s12.re - 0.05).abs() < 1e-10);

        // S22
        let s22 = ts.data[0].params[1][1];
        assert!((s22.re - 0.85).abs() < 1e-10);
    }

    #[test]
    fn parse_v2_with_keywords() {
        let input = "\
[Version] 2.0
! v2 example
[Number of Ports] 1
# GHz S RI R 50
[Network Data]
1.0  0.9  -0.1
2.0  0.8  -0.2
[End]
";
        let ts = parse(input).unwrap();
        assert_eq!(ts.version.as_deref(), Some("2.0"));
        assert_eq!(ts.num_ports, 1);
        assert_eq!(ts.data.len(), 2);
    }

    #[test]
    fn parse_v2_reference_impedances() {
        let input = "\
[Version] 2.0
[Number of Ports] 2
[Two-Port Data Order] 12_21
[Reference]
50.0 75.0
# GHz S RI R 50
[Network Data]
1.0  0.9 -0.1  0.05 0.01  0.05 0.01  0.85 -0.12
[End]
";
        let ts = parse(input).unwrap();
        assert_eq!(ts.two_port_order, Some(TwoPortOrder::Order12_21));
        assert_eq!(ts.reference_impedances, Some(vec![50.0, 75.0]));

        // With 12_21 order: S11, S12, S21, S22
        let s12 = ts.data[0].params[0][1];
        assert!((s12.re - 0.05).abs() < 1e-10);

        let s21 = ts.data[0].params[1][0];
        assert!((s21.re - 0.05).abs() < 1e-10);
    }

    #[test]
    fn parse_inline_comments() {
        let input = "\
# GHz S RI R 50
1.0  0.9  -0.1 ! this is inline
2.0  0.8  -0.2
";
        let ts = parse(input).unwrap();
        assert_eq!(ts.data.len(), 2);
    }

    #[test]
    fn error_no_option_line() {
        let input = "1.0  0.9  -0.1\n";
        let err = parse(input).unwrap_err();
        assert!(matches!(err, TouchstoneError::NoOptionLine));
    }

    #[test]
    fn error_no_data() {
        let input = "# GHz S RI R 50\n";
        let err = parse(input).unwrap_err();
        assert!(matches!(err, TouchstoneError::NoData));
    }

    #[test]
    fn error_invalid_number() {
        let input = "\
# GHz S RI R 50
1.0  abc  -0.1
";
        let err = parse(input).unwrap_err();
        assert!(matches!(err, TouchstoneError::InvalidNumber { .. }));
    }

    #[test]
    fn error_wrong_column_count() {
        let input = "\
# GHz S RI R 50
1.0  0.9
";
        let err = parse(input).unwrap_err();
        assert!(
            matches!(err, TouchstoneError::InvalidDataLine { .. } | TouchstoneError::ParseError { .. })
        );
    }

    #[test]
    fn frequencies_hz_conversion() {
        let input = "\
# GHz S RI R 50
1.0  0.9  -0.1
2.4  0.5  -0.5
";
        let ts = parse(input).unwrap();
        let freqs = ts.frequencies_hz();
        assert!((freqs[0] - 1e9).abs() < 1.0);
        assert!((freqs[1] - 2.4e9).abs() < 1.0);
    }

    #[test]
    fn convenience_methods() {
        let input = "\
# GHz S RI R 50
1.0  0.6  -0.8
";
        let ts = parse(input).unwrap();
        let mag_db = ts.magnitude_db(0, 0).unwrap();
        let phase = ts.phase_deg(0, 0).unwrap();

        let c = Complex::new(0.6, -0.8);
        assert!((mag_db[0] - c.magnitude_db()).abs() < 1e-10);
        assert!((phase[0] - c.phase_deg()).abs() < 1e-10);
    }
}
