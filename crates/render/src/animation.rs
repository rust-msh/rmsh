/// Phase animator for complex field data.
///
/// Computes `E(t) = Re(E) * cos(phase) - Im(E) * sin(phase)` to produce
/// time-domain field values from frequency-domain complex data.
pub struct PhaseAnimator {
    pub phase_deg: f32,
    pub playing: bool,
    pub speed_deg_per_sec: f32,
}

impl Default for PhaseAnimator {
    fn default() -> Self {
        Self {
            phase_deg: 0.0,
            playing: false,
            speed_deg_per_sec: 90.0,
        }
    }
}

impl PhaseAnimator {
    /// Advance the phase by `dt` seconds.
    pub fn tick(&mut self, dt: f32) {
        if self.playing {
            self.phase_deg = (self.phase_deg + self.speed_deg_per_sec * dt) % 360.0;
        }
    }

    /// Compute real field values at the current phase.
    pub fn apply(&self, field_real: &[f32], field_imag: &[f32]) -> Vec<f32> {
        let phase_rad = self.phase_deg.to_radians();
        let cos_p = phase_rad.cos();
        let sin_p = phase_rad.sin();
        field_real
            .iter()
            .zip(field_imag.iter())
            .map(|(&re, &im)| re * cos_p - im * sin_p)
            .collect()
    }

    /// Compute the field range over all possible phases (conservative envelope).
    pub fn envelope_range(field_real: &[f32], field_imag: &[f32]) -> [f32; 2] {
        let mut min_val = f32::MAX;
        let mut max_val = f32::MIN;
        for (&re, &im) in field_real.iter().zip(field_imag.iter()) {
            let amp = (re * re + im * im).sqrt();
            min_val = min_val.min(-amp);
            max_val = max_val.max(amp);
        }
        [min_val, max_val]
    }
}
