use crate::types::*;
use std::fmt::Write;

/// Write a Touchstone struct to a string in the Touchstone file format.
///
/// If `version` is set to "2.0", writes in Touchstone v2.0 format with keywords.
/// Otherwise writes in Touchstone v1.0 format.
pub fn write(ts: &Touchstone) -> String {
    let is_v2 = ts.version.as_deref() == Some("2.0");

    let mut out = String::new();

    if is_v2 {
        write_v2(&mut out, ts);
    } else {
        write_v1(&mut out, ts);
    }

    out
}

fn write_v1(out: &mut String, ts: &Touchstone) {
    // Comments
    for comment in &ts.comments {
        writeln!(out, "! {comment}").unwrap();
    }

    // Option line
    write_option_line(out, &ts.options);

    // Data
    write_data(out, ts);
}

fn write_v2(out: &mut String, ts: &Touchstone) {
    writeln!(out, "[Version] 2.0").unwrap();

    // Comments
    for comment in &ts.comments {
        writeln!(out, "! {comment}").unwrap();
    }

    writeln!(out, "[Number of Ports] {}", ts.num_ports).unwrap();

    if let Some(order) = ts.two_port_order {
        let order_str = match order {
            TwoPortOrder::Order21_12 => "21_12",
            TwoPortOrder::Order12_21 => "12_21",
        };
        writeln!(out, "[Two-Port Data Order] {order_str}").unwrap();
    }

    if let Some(ref impedances) = ts.reference_impedances {
        writeln!(out, "[Reference]").unwrap();
        let vals: Vec<String> = impedances.iter().map(|z| format_value(*z)).collect();
        writeln!(out, "{}", vals.join(" ")).unwrap();
    }

    // Option line
    write_option_line(out, &ts.options);

    writeln!(out, "[Network Data]").unwrap();
    write_data(out, ts);
    writeln!(out, "[End]").unwrap();
}

fn write_option_line(out: &mut String, opts: &TouchstoneOptions) {
    write!(
        out,
        "# {} {} {}",
        opts.frequency_unit.as_str(),
        opts.parameter_type.as_str(),
        opts.data_format.as_str(),
    )
    .unwrap();
    write!(out, " R {}", format_value(opts.reference_impedance)).unwrap();
    writeln!(out).unwrap();
}

fn write_data(out: &mut String, ts: &Touchstone) {
    let n = ts.num_ports;
    let format = ts.options.data_format;

    for fp in &ts.data {
        if n <= 2 {
            write_single_line(out, n, fp, format, ts.two_port_order);
        } else {
            write_multi_line(out, n, fp, format);
        }
    }
}

/// Write a frequency point on a single line (1-port and 2-port).
fn write_single_line(
    out: &mut String,
    n: usize,
    fp: &FrequencyPoint,
    format: DataFormat,
    two_port_order: Option<TwoPortOrder>,
) {
    write!(out, "{}", format_value(fp.frequency)).unwrap();

    if n == 1 {
        let (v1, v2) = from_complex(fp.params[0][0], format);
        write!(out, " {} {}", format_value(v1), format_value(v2)).unwrap();
    } else {
        // 2-port: order depends on two_port_order
        let order = two_port_order.unwrap_or(TwoPortOrder::Order21_12);
        let pairs: [(usize, usize); 4] = match order {
            TwoPortOrder::Order21_12 => [(0, 0), (1, 0), (0, 1), (1, 1)],
            TwoPortOrder::Order12_21 => [(0, 0), (0, 1), (1, 0), (1, 1)],
        };
        for (r, c) in pairs {
            let (v1, v2) = from_complex(fp.params[r][c], format);
            write!(out, " {} {}", format_value(v1), format_value(v2)).unwrap();
        }
    }
    writeln!(out).unwrap();
}

/// Write a frequency point in multi-line format (N >= 3 ports).
/// Row-major order: each row on its own line, first line includes frequency.
fn write_multi_line(
    out: &mut String,
    n: usize,
    fp: &FrequencyPoint,
    format: DataFormat,
) {
    for row in 0..n {
        if row == 0 {
            write!(out, "{}", format_value(fp.frequency)).unwrap();
        }
        for col in 0..n {
            let (v1, v2) = from_complex(fp.params[row][col], format);
            write!(out, " {} {}", format_value(v1), format_value(v2)).unwrap();
        }
        writeln!(out).unwrap();
    }
}

/// Convert a Complex value to a pair of f64 values in the given format.
fn from_complex(c: Complex, format: DataFormat) -> (f64, f64) {
    match format {
        DataFormat::RealImaginary => (c.re, c.im),
        DataFormat::MagnitudeAngle => (c.magnitude(), c.phase_deg()),
        DataFormat::DecibelAngle => (c.magnitude_db(), c.phase_deg()),
    }
}

/// Format a floating-point value, stripping unnecessary trailing zeros.
fn format_value(v: f64) -> String {
    if v == v.trunc() && v.abs() < 1e15 {
        // Integer-valued
        format!("{v:.1}")
    } else {
        format!("{v}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    #[test]
    fn roundtrip_1port_ri() {
        let input = "\
! test comment
# MHz S RI R 50.0
100.0 0.9 -0.1
200.0 0.8 -0.2
";
        let ts = parse(input).unwrap();
        let output = write(&ts);
        let ts2 = parse(&output).unwrap();

        assert_eq!(ts2.num_ports, 1);
        assert_eq!(ts2.data.len(), 2);
        assert!((ts2.data[0].params[0][0].re - 0.9).abs() < 1e-10);
        assert!((ts2.data[0].params[0][0].im - (-0.1)).abs() < 1e-10);
    }

    #[test]
    fn roundtrip_2port_ri() {
        let input = "\
# GHz S RI R 50.0
1.0 0.9 -0.1 0.05 0.01 0.05 0.01 0.85 -0.12
";
        let ts = parse(input).unwrap();
        let output = write(&ts);
        let ts2 = parse(&output).unwrap();

        assert_eq!(ts2.num_ports, 2);
        let s11 = ts2.data[0].params[0][0];
        assert!((s11.re - 0.9).abs() < 1e-10);
        let s21 = ts2.data[0].params[1][0];
        assert!((s21.re - 0.05).abs() < 1e-10);
    }

    #[test]
    fn roundtrip_v2() {
        let input = "\
[Version] 2.0
! v2 test
[Number of Ports] 1
# GHz S RI R 50.0
[Network Data]
1.0 0.9 -0.1
2.0 0.8 -0.2
[End]
";
        let ts = parse(input).unwrap();
        let output = write(&ts);
        let ts2 = parse(&output).unwrap();

        assert_eq!(ts2.version.as_deref(), Some("2.0"));
        assert_eq!(ts2.num_ports, 1);
        assert_eq!(ts2.data.len(), 2);
    }

    #[test]
    fn roundtrip_v2_reference_impedances() {
        let input = "\
[Version] 2.0
[Number of Ports] 2
[Two-Port Data Order] 12_21
[Reference]
50.0 75.0
# GHz S RI R 50.0
[Network Data]
1.0 0.9 -0.1 0.05 0.01 0.05 0.01 0.85 -0.12
[End]
";
        let ts = parse(input).unwrap();
        let output = write(&ts);
        let ts2 = parse(&output).unwrap();

        assert_eq!(ts2.two_port_order, Some(TwoPortOrder::Order12_21));
        assert_eq!(ts2.reference_impedances, Some(vec![50.0, 75.0]));
    }

    #[test]
    fn write_preserves_format() {
        let ts = Touchstone {
            num_ports: 1,
            options: TouchstoneOptions {
                frequency_unit: FrequencyUnit::MHz,
                parameter_type: ParameterType::S,
                data_format: DataFormat::MagnitudeAngle,
                reference_impedance: 50.0,
            },
            comments: vec!["antenna measurement".to_string()],
            data: vec![FrequencyPoint {
                frequency: 900.0,
                params: vec![vec![Complex::from_mag_angle(0.5, -30.0)]],
            }],
            version: None,
            two_port_order: None,
            reference_impedances: None,
        };

        let output = write(&ts);
        assert!(output.contains("# MHZ S MA"));
        assert!(output.contains("900.0"));
        assert!(output.contains("0.5"));
    }
}
