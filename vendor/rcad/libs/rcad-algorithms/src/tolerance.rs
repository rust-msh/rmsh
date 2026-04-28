use glam::DVec3;
use rcad_kernel::BRep;

/// Absolute tolerance for point coincidence.
///
/// Matches `rcad_kernel::tolerance::CONFUSION` = `Precision::Confusion()` in OCCT.
/// Two points are considered coincident when their distance is below this value.
pub const TOLERANCE_ABS: f64 = 1e-7;

/// Angular tolerance for parallel/perpendicular checks (radians, as cross-product magnitude).
///
/// This is intentionally **looser** than `rcad_kernel::tolerance::ANGULAR` (1e-12):
/// the algorithms layer needs to tolerate slightly imperfect parallelism that
/// arises from floating-point accumulation during intersection computation.
/// Used in [`vectors_parallel`] as `cross(a,b).length_squared() < TOLERANCE_ANG²`.
pub const TOLERANCE_ANG: f64 = 1e-9;

/// Tolerance squared — avoids `sqrt` in distance checks.
pub const TOLERANCE_ABS_SQ: f64 = TOLERANCE_ABS * TOLERANCE_ABS;

/// Tolerance level for different geometric operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceLevel {
    /// Strict tolerance for high-precision operations (e.g., intersection points).
    /// Scale factor: 1.0
    Strict,
    /// Normal tolerance for general operations (e.g., point classification).
    /// Scale factor: 10.0
    Normal,
    /// Relaxed tolerance for approximate operations (e.g., bounding box checks).
    /// Scale factor: 100.0
    Relaxed,
    /// Very relaxed tolerance for coarse operations (e.g., AABB pre-filter).
    /// Scale factor: 1000.0
    Coarse,
}

impl ToleranceLevel {
    /// Get the scale factor for this tolerance level.
    pub fn scale_factor(self) -> f64 {
        match self {
            ToleranceLevel::Strict => 1.0,
            ToleranceLevel::Normal => 10.0,
            ToleranceLevel::Relaxed => 100.0,
            ToleranceLevel::Coarse => 1000.0,
        }
    }
}

/// Adaptive tolerance context based on model scale.
///
/// Instead of using hard-coded absolute tolerances, this context computes
/// tolerances relative to the model's bounding box size. This ensures that
/// models at different scales (e.g., nanometer vs kilometer) are handled
/// appropriately.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveTolerance {
    /// Base tolerance (typically TOLERANCE_ABS for unit-scale models).
    pub base_tolerance: f64,
    /// Model scale factor (e.g., bounding box diagonal).
    pub model_scale: f64,
    /// Minimum tolerance to prevent excessive precision requirements.
    pub min_tolerance: f64,
    /// Maximum tolerance to prevent excessive looseness.
    pub max_tolerance: f64,
}

impl Default for AdaptiveTolerance {
    fn default() -> Self {
        Self {
            base_tolerance: TOLERANCE_ABS,
            model_scale: 1.0,
            min_tolerance: 1e-12,
            max_tolerance: 1e-3,
        }
    }
}

impl AdaptiveTolerance {
    /// Create a new adaptive tolerance with default base tolerance.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a new adaptive tolerance from a BRep's bounding box.
    pub fn from_brep(brep: &BRep) -> Self {
        let scale = compute_model_scale(brep);
        Self::from_scale(scale)
    }

    /// Create a new adaptive tolerance from two BReps' combined bounding box.
    pub fn from_two_breps(a: &BRep, b: &BRep) -> Self {
        let scale_a = compute_model_scale(a);
        let scale_b = compute_model_scale(b);
        Self::from_scale(scale_a.max(scale_b))
    }

    /// Create a new adaptive tolerance from a known scale.
    pub fn from_scale(model_scale: f64) -> Self {
        let mut ctx = Self::default();
        ctx.model_scale = model_scale.max(1e-10);
        ctx
    }

    /// Get the effective tolerance for a specific level.
    pub fn tolerance(self, level: ToleranceLevel) -> f64 {
        let raw = self.base_tolerance * level.scale_factor() * self.model_scale;
        raw.clamp(self.min_tolerance, self.max_tolerance)
    }

    /// Get the squared tolerance for a specific level.
    pub fn tolerance_sq(self, level: ToleranceLevel) -> f64 {
        let t = self.tolerance(level);
        t * t
    }

    /// Get the angular tolerance (not affected by model scale).
    pub fn angular_tolerance(self, level: ToleranceLevel) -> f64 {
        TOLERANCE_ANG * level.scale_factor()
    }

    /// Check if two points coincide at the given tolerance level.
    pub fn points_coincide_at(self, a: DVec3, b: DVec3, level: ToleranceLevel) -> bool {
        (a - b).length_squared() < self.tolerance_sq(level)
    }

    /// Check if a vector is zero at the given tolerance level.
    pub fn is_zero_vec_at(self, v: DVec3, level: ToleranceLevel) -> bool {
        v.length_squared() < self.tolerance_sq(level)
    }

    /// Check if two parameters are equal at the given tolerance level.
    pub fn params_equal_at(self, a: f64, b: f64, level: ToleranceLevel) -> bool {
        (a - b).abs() < self.tolerance(level)
    }

    /// Get tolerance for point coincidence (strict level).
    pub fn coincidence(self) -> f64 {
        self.tolerance(ToleranceLevel::Strict)
    }

    /// Get tolerance for classification operations (normal level).
    pub fn classification(self) -> f64 {
        self.tolerance(ToleranceLevel::Normal)
    }

    /// Get tolerance for boundary checks (relaxed level).
    pub fn boundary(self) -> f64 {
        self.tolerance(ToleranceLevel::Relaxed)
    }

    /// Get tolerance for AABB pre-filtering (coarse level).
    pub fn coarse(self) -> f64 {
        self.tolerance(ToleranceLevel::Coarse)
    }
}

/// Compute the characteristic scale of a model from its bounding box.
/// Returns the diagonal of the bounding box, or 1.0 if the model is empty.
pub fn compute_model_scale(brep: &BRep) -> f64 {
    let mut min_pt = DVec3::splat(f64::INFINITY);
    let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
    let mut has_vertices = false;

    for vertex in &brep.vertices {
        min_pt = min_pt.min(vertex.point);
        max_pt = max_pt.max(vertex.point);
        has_vertices = true;
    }

    if !has_vertices {
        return 1.0;
    }

    let diagonal = (max_pt - min_pt).length();
    diagonal.max(1e-10)
}

/// Compute the characteristic scale from a collection of points.
pub fn compute_scale_from_points(points: &[DVec3]) -> f64 {
    if points.is_empty() {
        return 1.0;
    }

    let mut min_pt = DVec3::splat(f64::INFINITY);
    let mut max_pt = DVec3::splat(f64::NEG_INFINITY);

    for &p in points {
        min_pt = min_pt.min(p);
        max_pt = max_pt.max(p);
    }

    let diagonal = (max_pt - min_pt).length();
    diagonal.max(1e-10)
}

#[inline]
pub fn points_coincide(a: DVec3, b: DVec3) -> bool {
    (a - b).length_squared() < TOLERANCE_ABS_SQ
}

#[inline]
pub fn is_zero_vec(v: DVec3) -> bool {
    v.length_squared() < TOLERANCE_ABS_SQ
}

/// Returns true if two unit vectors are parallel (or anti-parallel).
#[inline]
pub fn vectors_parallel(a: DVec3, b: DVec3) -> bool {
    a.cross(b).length_squared() < TOLERANCE_ANG * TOLERANCE_ANG
}

#[inline]
pub fn params_equal(a: f64, b: f64) -> bool {
    (a - b).abs() < TOLERANCE_ABS
}

/// Check if two vectors are parallel using adaptive tolerance.
pub fn vectors_parallel_adaptive(a: DVec3, b: DVec3, tol: AdaptiveTolerance) -> bool {
    let ang_tol = tol.angular_tolerance(ToleranceLevel::Normal);
    a.cross(b).length_squared() < ang_tol * ang_tol
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coincident_points() {
        let a = DVec3::new(1.0, 2.0, 3.0);
        let b = DVec3::new(1.0, 2.0, 3.0 + 1e-8);
        assert!(points_coincide(a, b));
    }

    #[test]
    fn non_coincident_points() {
        let a = DVec3::ZERO;
        let b = DVec3::new(1e-6, 0.0, 0.0);
        assert!(!points_coincide(a, b));
    }

    #[test]
    fn parallel_vectors() {
        assert!(vectors_parallel(DVec3::X, DVec3::X));
        assert!(vectors_parallel(DVec3::X, -DVec3::X));
        assert!(!vectors_parallel(DVec3::X, DVec3::Y));
    }

    #[test]
    fn tolerance_level_scale_factors() {
        assert_eq!(ToleranceLevel::Strict.scale_factor(), 1.0);
        assert_eq!(ToleranceLevel::Normal.scale_factor(), 10.0);
        assert_eq!(ToleranceLevel::Relaxed.scale_factor(), 100.0);
        assert_eq!(ToleranceLevel::Coarse.scale_factor(), 1000.0);
    }

    #[test]
    fn adaptive_tolerance_default() {
        let tol = AdaptiveTolerance::default();
        assert_eq!(tol.tolerance(ToleranceLevel::Strict), TOLERANCE_ABS);
        assert_eq!(tol.tolerance(ToleranceLevel::Normal), TOLERANCE_ABS * 10.0);
    }

    #[test]
    fn adaptive_tolerance_with_scale() {
        let tol = AdaptiveTolerance::from_scale(100.0);
        // At 100x scale, tolerance should be 100x larger
        assert!((tol.tolerance(ToleranceLevel::Strict) - TOLERANCE_ABS * 100.0).abs() < 1e-20);
    }

    #[test]
    fn adaptive_tolerance_clamping() {
        // Very large scale should be clamped to max_tolerance
        let tol = AdaptiveTolerance::from_scale(1e10);
        assert_eq!(tol.tolerance(ToleranceLevel::Strict), tol.max_tolerance);

        // Very small scale should be clamped to min_tolerance
        let mut small_tol = AdaptiveTolerance::from_scale(1e-20);
        small_tol.min_tolerance = 1e-15;
        assert_eq!(small_tol.tolerance(ToleranceLevel::Strict), small_tol.min_tolerance);
    }

    #[test]
    fn compute_model_scale_from_points() {
        let points = vec![
            DVec3::ZERO,
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let scale = compute_scale_from_points(&points);
        // Bounding box is [0,1] x [0,1] x [0,0], diagonal is sqrt(2)
        assert!((scale - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn compute_model_scale_empty() {
        let scale = compute_scale_from_points(&[]);
        assert_eq!(scale, 1.0);
    }

    #[test]
    fn adaptive_points_coincide() {
        let tol = AdaptiveTolerance::default();
        let a = DVec3::new(1.0, 2.0, 3.0);
        let b = DVec3::new(1.0, 2.0, 3.0 + 1e-8);
        assert!(tol.points_coincide_at(a, b, ToleranceLevel::Strict));
    }

    #[test]
    fn convenience_methods() {
        let tol = AdaptiveTolerance::default();
        // Just verify these don't panic and return reasonable values
        assert!(tol.coincidence() > 0.0);
        assert!(tol.classification() > tol.coincidence());
        assert!(tol.boundary() > tol.classification());
        assert!(tol.coarse() > tol.boundary());
    }
}
