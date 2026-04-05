#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColormapType {
    Rainbow,
    Viridis,
    CoolWarm,
    Grayscale,
}

impl ColormapType {
    pub const ALL: &[ColormapType] = &[
        ColormapType::Rainbow,
        ColormapType::Viridis,
        ColormapType::CoolWarm,
        ColormapType::Grayscale,
    ];

    pub fn name(&self) -> &'static str {
        match self {
            Self::Rainbow => "Rainbow",
            Self::Viridis => "Viridis",
            Self::CoolWarm => "Cool-Warm",
            Self::Grayscale => "Grayscale",
        }
    }

    /// Generate an RGBA8 lookup table with `n` entries.
    pub fn generate_lut(&self, n: usize) -> Vec<[u8; 4]> {
        (0..n)
            .map(|i| {
                let t = i as f32 / (n - 1).max(1) as f32;
                self.sample(t)
            })
            .collect()
    }

    fn sample(&self, t: f32) -> [u8; 4] {
        let t = t.clamp(0.0, 1.0);
        match self {
            Self::Rainbow => sample_rainbow(t),
            Self::Viridis => sample_piecewise(t, &VIRIDIS_STOPS),
            Self::CoolWarm => sample_piecewise(t, &COOLWARM_STOPS),
            Self::Grayscale => {
                let v = (t * 255.0) as u8;
                [v, v, v, 255]
            }
        }
    }
}

fn sample_rainbow(t: f32) -> [u8; 4] {
    // HSV sweep: hue from 240° (blue) to 0° (red)
    let h = (1.0 - t) * 240.0;
    let s = 1.0f32;
    let v = 1.0f32;

    let c = v * s;
    let h_prime = h / 60.0;
    let x = c * (1.0 - (h_prime % 2.0 - 1.0).abs());
    let m = v - c;

    let (r, g, b) = if h_prime < 1.0 {
        (c, x, 0.0)
    } else if h_prime < 2.0 {
        (x, c, 0.0)
    } else if h_prime < 3.0 {
        (0.0, c, x)
    } else if h_prime < 4.0 {
        (0.0, x, c)
    } else if h_prime < 5.0 {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };

    [
        ((r + m) * 255.0) as u8,
        ((g + m) * 255.0) as u8,
        ((b + m) * 255.0) as u8,
        255,
    ]
}

type ColorStop = (f32, [f32; 3]);

const VIRIDIS_STOPS: [ColorStop; 9] = [
    (0.000, [0.267, 0.004, 0.329]),
    (0.125, [0.282, 0.141, 0.458]),
    (0.250, [0.253, 0.265, 0.530]),
    (0.375, [0.206, 0.372, 0.553]),
    (0.500, [0.163, 0.471, 0.558]),
    (0.625, [0.128, 0.567, 0.551]),
    (0.750, [0.135, 0.659, 0.518]),
    (0.875, [0.267, 0.749, 0.441]),
    (1.000, [0.993, 0.906, 0.144]),
];

const COOLWARM_STOPS: [ColorStop; 5] = [
    (0.00, [0.230, 0.299, 0.754]),
    (0.25, [0.552, 0.691, 0.996]),
    (0.50, [0.866, 0.866, 0.866]),
    (0.75, [0.956, 0.604, 0.486]),
    (1.00, [0.706, 0.016, 0.150]),
];

fn sample_piecewise(t: f32, stops: &[ColorStop]) -> [u8; 4] {
    if t <= stops[0].0 {
        return float_to_u8(stops[0].1);
    }
    for i in 1..stops.len() {
        if t <= stops[i].0 {
            let (t0, c0) = stops[i - 1];
            let (t1, c1) = stops[i];
            let f = (t - t0) / (t1 - t0);
            return float_to_u8([
                c0[0] + f * (c1[0] - c0[0]),
                c0[1] + f * (c1[1] - c0[1]),
                c0[2] + f * (c1[2] - c0[2]),
            ]);
        }
    }
    float_to_u8(stops[stops.len() - 1].1)
}

fn float_to_u8(c: [f32; 3]) -> [u8; 4] {
    [
        (c[0] * 255.0) as u8,
        (c[1] * 255.0) as u8,
        (c[2] * 255.0) as u8,
        255,
    ]
}
