use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

// ========================
// Enums
// ========================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FrequencyUnit {
    Hz,
    KHz,
    MHz,
    GHz,
}

impl FrequencyUnit {
    pub fn multiplier(self) -> f64 {
        match self {
            FrequencyUnit::Hz => 1.0,
            FrequencyUnit::KHz => 1e3,
            FrequencyUnit::MHz => 1e6,
            FrequencyUnit::GHz => 1e9,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            FrequencyUnit::Hz => "HZ",
            FrequencyUnit::KHz => "KHZ",
            FrequencyUnit::MHz => "MHZ",
            FrequencyUnit::GHz => "GHZ",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ParameterType {
    S,
    Y,
    Z,
    H,
    G,
}

impl ParameterType {
    pub fn as_str(self) -> &'static str {
        match self {
            ParameterType::S => "S",
            ParameterType::Y => "Y",
            ParameterType::Z => "Z",
            ParameterType::H => "H",
            ParameterType::G => "G",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DataFormat {
    /// Real and imaginary parts
    RealImaginary,
    /// Magnitude and angle (degrees)
    MagnitudeAngle,
    /// Magnitude (dB) and angle (degrees)
    DecibelAngle,
}

impl DataFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            DataFormat::RealImaginary => "RI",
            DataFormat::MagnitudeAngle => "MA",
            DataFormat::DecibelAngle => "DB",
        }
    }
}

/// Two-port data ordering (Touchstone v2.0)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TwoPortOrder {
    /// S11, S21, S12, S22 (default for v2.0)
    Order21_12,
    /// S11, S12, S21, S22
    Order12_21,
}

// ========================
// Complex number
// ========================

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub fn zero() -> Self {
        Self { re: 0.0, im: 0.0 }
    }

    pub fn magnitude(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }

    pub fn magnitude_db(&self) -> f64 {
        20.0 * self.magnitude().log10()
    }

    pub fn phase_rad(&self) -> f64 {
        self.im.atan2(self.re)
    }

    pub fn phase_deg(&self) -> f64 {
        self.phase_rad() * 180.0 / PI
    }

    /// Create from magnitude and angle (degrees)
    pub fn from_mag_angle(mag: f64, angle_deg: f64) -> Self {
        let angle_rad = angle_deg * PI / 180.0;
        Self {
            re: mag * angle_rad.cos(),
            im: mag * angle_rad.sin(),
        }
    }

    /// Create from dB magnitude and angle (degrees)
    pub fn from_db_angle(db: f64, angle_deg: f64) -> Self {
        let mag = 10.0_f64.powf(db / 20.0);
        Self::from_mag_angle(mag, angle_deg)
    }
}

// ========================
// Options
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TouchstoneOptions {
    pub frequency_unit: FrequencyUnit,
    pub parameter_type: ParameterType,
    pub data_format: DataFormat,
    pub reference_impedance: f64,
}

impl Default for TouchstoneOptions {
    fn default() -> Self {
        Self {
            frequency_unit: FrequencyUnit::GHz,
            parameter_type: ParameterType::S,
            data_format: DataFormat::MagnitudeAngle,
            reference_impedance: 50.0,
        }
    }
}

// ========================
// Frequency point
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrequencyPoint {
    /// Frequency value in the unit specified by options
    pub frequency: f64,
    /// NxN parameter matrix, stored as Real/Imaginary internally.
    /// Indexed as params\[row\]\[col\] (0-based).
    pub params: Vec<Vec<Complex>>,
}

// ========================
// Touchstone dataset
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Touchstone {
    pub num_ports: usize,
    pub options: TouchstoneOptions,
    pub comments: Vec<String>,
    pub data: Vec<FrequencyPoint>,
    /// Touchstone version if v2.0+ keywords detected (e.g. "2.0")
    pub version: Option<String>,
    /// Two-port data ordering (v2.0 only)
    pub two_port_order: Option<TwoPortOrder>,
    /// Per-port reference impedances (v2.0 [Reference] keyword)
    pub reference_impedances: Option<Vec<f64>>,
}

impl Touchstone {
    /// Get the parameter at (row, col) across all frequencies.
    pub fn get_parameter(&self, row: usize, col: usize) -> Option<Vec<Complex>> {
        if row >= self.num_ports || col >= self.num_ports {
            return None;
        }
        Some(self.data.iter().map(|fp| fp.params[row][col]).collect())
    }

    /// Get all frequencies converted to Hz.
    pub fn frequencies_hz(&self) -> Vec<f64> {
        let mult = self.options.frequency_unit.multiplier();
        self.data.iter().map(|fp| fp.frequency * mult).collect()
    }

    /// Get magnitude in dB for parameter (row, col) at all frequencies.
    pub fn magnitude_db(&self, row: usize, col: usize) -> Option<Vec<f64>> {
        self.get_parameter(row, col)
            .map(|v| v.iter().map(|c| c.magnitude_db()).collect())
    }

    /// Get phase in degrees for parameter (row, col) at all frequencies.
    pub fn phase_deg(&self, row: usize, col: usize) -> Option<Vec<f64>> {
        self.get_parameter(row, col)
            .map(|v| v.iter().map(|c| c.phase_deg()).collect())
    }

    /// Number of frequency points.
    pub fn num_frequencies(&self) -> usize {
        self.data.len()
    }
}
