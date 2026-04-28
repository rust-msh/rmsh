//! Visual appearance: per-face / per-solid RGB color and surface finish.
//!
//! Analogous to OCCT `XCAFDoc_ColorTool` + `XCAFDoc_VisMaterial`.
//!
//! Colors are stored separately from the BRep topology so that geometry
//! remains pure and color information can be added or removed independently.

use serde::{Deserialize, Serialize};

/// An sRGB color with components in [0.0, 1.0].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Color {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl Color {
    pub const fn new(r: f64, g: f64, b: f64) -> Self {
        Self { r, g, b }
    }

    /// Create from 8-bit RGB values.
    pub fn from_rgb8(r: u8, g: u8, b: u8) -> Self {
        Self {
            r: r as f64 / 255.0,
            g: g as f64 / 255.0,
            b: b as f64 / 255.0,
        }
    }

    // ── Preset colors ─────────────────────────────────────────────────
    pub const RED: Self = Self::new(1.0, 0.0, 0.0);
    pub const GREEN: Self = Self::new(0.0, 0.8, 0.0);
    pub const BLUE: Self = Self::new(0.0, 0.4, 1.0);
    pub const YELLOW: Self = Self::new(1.0, 0.9, 0.0);
    pub const CYAN: Self = Self::new(0.0, 0.9, 0.9);
    pub const MAGENTA: Self = Self::new(0.9, 0.0, 0.9);
    pub const WHITE: Self = Self::new(1.0, 1.0, 1.0);
    pub const GRAY: Self = Self::new(0.6, 0.6, 0.6);
    pub const SILVER: Self = Self::new(0.75, 0.75, 0.78);
    pub const GOLD: Self = Self::new(1.0, 0.84, 0.0);
    pub const ORANGE: Self = Self::new(1.0, 0.5, 0.0);
    pub const BLACK: Self = Self::new(0.05, 0.05, 0.05);
}

/// Per-face color assignment.
///
/// `face_index` is the flat face index across all shells of the BRep
/// (same indexing as `GeomStore::face_surface`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FaceColor {
    pub face_index: usize,
    pub color: Color,
}

/// A collection of color assignments for a single BRep.
///
/// Stores an optional solid-level (default) color and per-face overrides.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct StepColor {
    /// Fallback color applied to any face without an explicit override.
    pub solid_color: Option<Color>,
    /// Per-face color overrides.
    pub face_colors: Vec<FaceColor>,
}

impl StepColor {
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the default solid-level color.
    pub fn with_solid_color(mut self, color: Color) -> Self {
        self.solid_color = Some(color);
        self
    }

    /// Assign a color to a specific face index.
    pub fn with_face_color(mut self, face_index: usize, color: Color) -> Self {
        self.face_colors.push(FaceColor { face_index, color });
        self
    }

    /// Look up the color for a given face index.
    /// Returns the face override if present, otherwise the solid color.
    pub fn color_for_face(&self, face_index: usize) -> Option<Color> {
        self.face_colors
            .iter()
            .find(|fc| fc.face_index == face_index)
            .map(|fc| fc.color)
            .or(self.solid_color)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn color_from_rgb8_round_trip() {
        let c = Color::from_rgb8(255, 128, 0);
        assert_eq!(c.r, 1.0);
        assert!((c.g - 128.0 / 255.0).abs() < 1e-12);
        assert_eq!(c.b, 0.0);
    }

    #[test]
    fn color_preset_values() {
        assert_eq!(Color::RED, Color::new(1.0, 0.0, 0.0));
        assert_eq!(Color::WHITE, Color::new(1.0, 1.0, 1.0));
        assert_eq!(Color::BLACK, Color::new(0.05, 0.05, 0.05));
    }

    #[test]
    fn step_color_solid_fallback() {
        let sc = StepColor::new().with_solid_color(Color::BLUE);
        // Any face without override gets the solid color.
        assert_eq!(sc.color_for_face(0), Some(Color::BLUE));
        assert_eq!(sc.color_for_face(99), Some(Color::BLUE));
    }

    #[test]
    fn step_color_face_override() {
        let sc = StepColor::new()
            .with_solid_color(Color::GRAY)
            .with_face_color(2, Color::RED);
        assert_eq!(sc.color_for_face(2), Some(Color::RED));
        // Other faces fall back to solid.
        assert_eq!(sc.color_for_face(0), Some(Color::GRAY));
    }

    #[test]
    fn step_color_no_color_returns_none() {
        let sc = StepColor::new();
        assert_eq!(sc.color_for_face(0), None);
    }

    #[test]
    fn step_color_face_color_without_solid() {
        let sc = StepColor::new().with_face_color(1, Color::GREEN);
        assert_eq!(sc.color_for_face(1), Some(Color::GREEN));
        // Face 2 has no override and no solid color.
        assert_eq!(sc.color_for_face(2), None);
    }
}
