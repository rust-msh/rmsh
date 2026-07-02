//! Scalar size field system — Gmsh `Field` equivalent.
//!
//! Provides per-point characteristic length `lc(x,y,z)` that meshing algorithms
//! use to determine local element size.  This is distinct from the Riemannian
//! metric fields (`MetricField2D`/`MetricField3D`) used by BAMG/MMG3D for
//! anisotropic meshing.
//!
//! # Architecture
//!
//! [`SizeField`] is the trait.  Concrete field types (Distance, Threshold, …)
//! implement it.  [`FieldManager`] holds all registered fields by tag ID and
//! provides the evaluation entry point used by meshers.
//!
//! # Gmsh field type mapping
//!
//! | Gmsh type | rmsh field | Description |
//! |---|---|---|
//! | `Distance` | [`DistanceField`] | lc = distance to nearest entity |
//! | `Threshold` | [`ThresholdField`] | smooth step between min/max lc vs distance |
//! | `MathEval` | [`MathEvalField`] | lc = f(x,y,z) via expression |
//! | `Min` | [`MinField`] | lc = min of child fields |
//! | `Max` | [`MaxField`] | lc = max of child fields |
//! | `Box` | [`BoxField`] | lc = constant inside box, blend at boundary |
//! | `Ball` | [`BallField`] | lc = constant inside ball, blend at boundary |
//! | `Restrict` | [`RestrictField`] | lc = child field inside, ∞ outside |
//! | `Constant` | [`ConstantField`] | lc = constant everywhere |

use std::collections::HashMap;
use std::f64;

// ─── Error type ───────────────────────────────────────────────────────────────

/// Errors produced by size field operations.
#[derive(Debug, Clone)]
pub enum FieldError {
    /// A referenced field tag does not exist.
    UnknownField(i32),
    /// Invalid or missing parameter for a field.
    InvalidParameter(String),
}

impl std::fmt::Display for FieldError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FieldError::UnknownField(tag) => write!(f, "unknown field tag {tag}"),
            FieldError::InvalidParameter(msg) => write!(f, "invalid field parameter: {msg}"),
        }
    }
}

// ─── Trait ────────────────────────────────────────────────────────────────────

/// A scalar size field that returns the characteristic length at a point.
///
/// All coordinates are in model-space length units.  A return value of
/// `f64::INFINITY` means "no constraint" (the mesher falls back to its
/// native size).
pub trait SizeField: Send + Sync {
    /// Human-readable type name, e.g. `"Distance"`, `"Threshold"`.
    fn type_name(&self) -> &'static str;

    /// Evaluate the characteristic length `lc` at `(x, y, z)`.
    fn lc_at(&self, x: f64, y: f64, z: f64, deps: &dyn FieldLookup) -> f64;
}

/// Lookup interface for a field to query other fields by tag.
pub trait FieldLookup {
    fn evaluate_field(&self, tag: i32, x: f64, y: f64, z: f64) -> f64;
}

// ─── Field types ──────────────────────────────────────────────────────────────

/// Constant field: returns the same `lc` everywhere.
#[derive(Debug, Clone)]
pub struct ConstantField {
    pub lc: f64,
}

impl SizeField for ConstantField {
    fn type_name(&self) -> &'static str { "Constant" }
    fn lc_at(&self, _x: f64, _y: f64, _z: f64, _deps: &dyn FieldLookup) -> f64 {
        self.lc
    }
}

/// Distance field: lc = distance to a model entity (currently stub: uses
/// distance to a user-specified point).
///
/// Gmsh `Field.Distance` computes exact distance to a curve/surface entity.
/// This simplified version computes `dist(cx,cy,cz, x,y,z)` when `source_type=Point`.
/// For full gmsh compatibility, wire to rcad2 BRep distance evaluation.
#[derive(Debug, Clone)]
pub struct DistanceField {
    /// Source type: "Point", "Curve", "Surface", "Field"
    pub source_type: String,
    /// Source point (used when source_type = "Point")
    pub cx: f64,
    pub cy: f64,
    pub cz: f64,
    /// Source field tag (used when source_type = "Field")
    pub source_tag: i32,
}

impl SizeField for DistanceField {
    fn type_name(&self) -> &'static str { "Distance" }
    fn lc_at(&self, x: f64, y: f64, z: f64, deps: &dyn FieldLookup) -> f64 {
        if self.source_type == "Field" {
            // Wrapper: distance to a scalar field's isosurface (approximated
            // via the field value at the point itself — not true distance).
            deps.evaluate_field(self.source_tag, x, y, z)
        } else {
            // Point distance (Euclidean)
            let dx = x - self.cx;
            let dy = y - self.cy;
            let dz = z - self.cz;
            (dx * dx + dy * dy + dz * dz).sqrt()
        }
    }
}

/// Threshold field: smooth step between `lc_min` and `lc_max` based on
/// distance from `in_field`.
///
/// - If `dist < dist_min`: returns `lc_min`
/// - If `dist > dist_max`: returns `lc_max`
/// - Otherwise: cubic Hermite smooth step between `lc_min` and `lc_max`
#[derive(Debug, Clone)]
pub struct ThresholdField {
    /// Field tag providing the distance values.
    pub in_field: i32,
    /// Characteristic length at `dist_min` (inside).
    pub lc_min: f64,
    /// Characteristic length at `dist_max` (outside).
    pub lc_max: f64,
    /// Distance at which lc reaches lc_min.
    pub dist_min: f64,
    /// Distance at which lc reaches lc_max.
    pub dist_max: f64,
}

impl SizeField for ThresholdField {
    fn type_name(&self) -> &'static str { "Threshold" }
    fn lc_at(&self, x: f64, y: f64, z: f64, deps: &dyn FieldLookup) -> f64 {
        let d = deps.evaluate_field(self.in_field, x, y, z);
        if d <= self.dist_min {
            self.lc_min
        } else if d >= self.dist_max {
            self.lc_max
        } else {
            // Cubic Hermite smooth step: t = (d - dist_min) / (dist_max - dist_min)
            let t = (d - self.dist_min) / (self.dist_max - self.dist_min);
            let t2 = t * t;
            let t3 = t2 * t;
            let s = 3.0 * t2 - 2.0 * t3; // smoothstep 0→1
            self.lc_min + (self.lc_max - self.lc_min) * s
        }
    }
}

/// MathEval field: lc = a simple math expression (currently stub: constant).
///
/// Gmsh evaluates arbitrary `f(x,y,z)` expressions via its built-in parser.
/// rmsh will wire to a tiny expression evaluator in a future PR.
#[derive(Debug, Clone)]
pub struct MathEvalField {
    /// The expression string (stored for round-trip fidelity).
    pub expression: String,
    /// Hard-coded override for now (expression evaluator TBD).
    pub fallback_lc: f64,
}

impl SizeField for MathEvalField {
    fn type_name(&self) -> &'static str { "MathEval" }
    fn lc_at(&self, _x: f64, _y: f64, _z: f64, _deps: &dyn FieldLookup) -> f64 {
        // TODO: wire to a tiny expression evaluator (meval or fasteval)
        self.fallback_lc
    }
}

/// Min field: lc = min of all child fields.
#[derive(Debug, Clone)]
pub struct MinField {
    pub fields: Vec<i32>,
}

impl SizeField for MinField {
    fn type_name(&self) -> &'static str { "Min" }
    fn lc_at(&self, x: f64, y: f64, z: f64, deps: &dyn FieldLookup) -> f64 {
        self.fields
            .iter()
            .map(|&tag| deps.evaluate_field(tag, x, y, z))
            .fold(f64::INFINITY, f64::min)
    }
}

/// Max field: lc = max of all child fields.
#[derive(Debug, Clone)]
pub struct MaxField {
    pub fields: Vec<i32>,
}

impl SizeField for MaxField {
    fn type_name(&self) -> &'static str { "Max" }
    fn lc_at(&self, x: f64, y: f64, z: f64, deps: &dyn FieldLookup) -> f64 {
        self.fields
            .iter()
            .map(|&tag| deps.evaluate_field(tag, x, y, z))
            .fold(0.0_f64, f64::max)
    }
}

/// Box field: constant `lc` inside an axis-aligned box, interpolated at
/// the boundary region.
#[derive(Debug, Clone)]
pub struct BoxField {
    pub lc_inside: f64,
    pub lc_outside: f64,
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
    pub z_min: f64,
    pub z_max: f64,
}

impl SizeField for BoxField {
    fn type_name(&self) -> &'static str { "Box" }
    fn lc_at(&self, x: f64, y: f64, z: f64, _deps: &dyn FieldLookup) -> f64 {
        // Normalised distance from the box in each axis: 0 inside, 1 at ±thickness.
        let d = |v, vmin, vmax| {
            if v < vmin { vmin - v } else if v > vmax { v - vmax } else { 0.0 }
        };
        let dx = d(x, self.x_min, self.x_max);
        let dy = d(y, self.y_min, self.y_max);
        let dz = d(z, self.z_min, self.z_max);
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        // Sharp box: no blend region for now.
        if dist <= 0.0 { self.lc_inside } else { self.lc_outside }
    }
}

/// Restrict field: delegates to a child field inside a region, returns
/// `INFINITY` (no constraint) outside.
#[derive(Debug, Clone)]
pub struct RestrictField {
    pub in_field: i32,
    /// If non-empty, restrict to inside the box defined by these bounds.
    /// Format: `[xmin, xmax, ymin, ymax, zmin, zmax]`
    pub box_bounds: Option<[f64; 6]>,
}

impl SizeField for RestrictField {
    fn type_name(&self) -> &'static str { "Restrict" }
    fn lc_at(&self, x: f64, y: f64, z: f64, deps: &dyn FieldLookup) -> f64 {
        if let Some([x1, x2, y1, y2, z1, z2]) = self.box_bounds {
            if x < x1 || x > x2 || y < y1 || y > y2 || z < z1 || z > z2 {
                return f64::INFINITY;
            }
        }
        deps.evaluate_field(self.in_field, x, y, z)
    }
}

// ─── Field manager ────────────────────────────────────────────────────────────

/// Registry of all size fields with a background field pointer for meshers.
///
/// ## Usage
///
/// ```rust,ignore
/// let mut mgr = FieldManager::new();
/// let d_tag = mgr.add(Box::new(DistanceField { cx: 0.0, cy: 0.0, cz: 0.0,
///     source_type: "Point".into(), source_tag: 0 }));
/// let t_tag = mgr.add(Box::new(ThresholdField {
///     in_field: d_tag, lc_min: 0.1, lc_max: 1.0, dist_min: 0.0, dist_max: 5.0 }));
/// mgr.set_background(t_tag);
///
/// // Mesher queries:
/// let lc = mgr.evaluate(x, y, z);  // uses background field
/// ```
pub struct FieldManager {
    fields: HashMap<i32, Box<dyn SizeField>>,
    next_tag: i32,
    background_tag: Option<i32>,
}

impl FieldManager {
    pub fn new() -> Self {
        Self {
            fields: HashMap::new(),
            next_tag: 1,
            background_tag: None,
        }
    }

    /// Register a new field and return its tag ID.
    pub fn add(&mut self, field: Box<dyn SizeField>) -> i32 {
        let tag = self.next_tag;
        self.fields.insert(tag, field);
        self.next_tag += 1;
        tag
    }

    /// Replace an existing field at `tag`.
    pub fn set(&mut self, tag: i32, field: Box<dyn SizeField>) -> Result<(), FieldError> {
        if !self.fields.contains_key(&tag) {
            return Err(FieldError::UnknownField(tag));
        }
        self.fields.insert(tag, field);
        Ok(())
    }

    /// Set the background mesh field tag (the one meshers evaluate).
    pub fn set_background(&mut self, tag: i32) -> Result<(), FieldError> {
        if !self.fields.contains_key(&tag) {
            return Err(FieldError::UnknownField(tag));
        }
        self.background_tag = Some(tag);
        Ok(())
    }

    /// Evaluate the background field at `(x, y, z)`.  Returns `None` when no
    /// background field is set.
    pub fn evaluate(&self, x: f64, y: f64, z: f64) -> Option<f64> {
        let tag = self.background_tag?;
        Some(self.evaluate_field(tag, x, y, z))
    }

    /// Whether a background field is set.
    pub fn has_background(&self) -> bool {
        self.background_tag.is_some()
    }

    /// Number of registered fields.
    pub fn len(&self) -> usize {
        self.fields.len()
    }

    /// Whether any field is registered.
    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    /// List all field tags.
    pub fn tags(&self) -> Vec<i32> {
        self.fields.keys().copied().collect()
    }
}

impl FieldLookup for FieldManager {
    fn evaluate_field(&self, tag: i32, x: f64, y: f64, z: f64) -> f64 {
        match self.fields.get(&tag) {
            Some(field) => field.lc_at(x, y, z, self),
            None => f64::INFINITY,
        }
    }
}

// ─── Wrapper: make SizeField usable as a MetricField2D/3D bridge ──────────

/// Bridge: wraps any `SizeField` as a uniform-isotropic `MetricField2D`.
///
/// The scalar lc is mapped to an isotropic metric `M = (1/lc²)·I`.
pub struct SizeAsMetric2D<'a> {
    field_tag: i32,
    manager: &'a FieldManager,
}

impl<'a> SizeAsMetric2D<'a> {
    pub fn new(manager: &'a FieldManager, field_tag: i32) -> Self {
        Self { field_tag, manager }
    }
}

impl crate::bamg_2d::MetricField2D for SizeAsMetric2D<'_> {
    fn metric_at(&self, x: f64, y: f64) -> crate::bamg_2d::Metric2 {
        let lc = self.manager.evaluate_field(self.field_tag, x, y, 0.0);
        if lc.is_finite() && lc > 0.0 {
            crate::bamg_2d::Metric2::isotropic(lc)
        } else {
            crate::bamg_2d::Metric2::isotropic(1.0)
        }
    }
}

/// Bridge: wraps any `SizeField` as a uniform-isotropic `MetricField3D`.
pub struct SizeAsMetric3D<'a> {
    field_tag: i32,
    manager: &'a FieldManager,
}

impl<'a> SizeAsMetric3D<'a> {
    pub fn new(manager: &'a FieldManager, field_tag: i32) -> Self {
        Self { field_tag, manager }
    }
}

impl crate::mmg_remesh::MetricField3D for SizeAsMetric3D<'_> {
    fn metric_at(&self, x: f64, y: f64, z: f64) -> crate::mmg_remesh::Metric3 {
        let lc = self.manager.evaluate_field(self.field_tag, x, y, z);
        if lc.is_finite() && lc > 0.0 {
            crate::mmg_remesh::Metric3::isotropic(lc)
        } else {
            crate::mmg_remesh::Metric3::isotropic(1.0)
        }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn empty_lookup() -> impl FieldLookup {
        struct Empty;
        impl FieldLookup for Empty {
            fn evaluate_field(&self, _tag: i32, _x: f64, _y: f64, _z: f64) -> f64 {
                f64::INFINITY
            }
        }
        Empty
    }

    #[test]
    fn constant_field_returns_fixed_value() {
        let f = ConstantField { lc: 0.5 };
        let deps = empty_lookup();
        assert_eq!(f.lc_at(0.0, 0.0, 0.0, &deps), 0.5);
        assert_eq!(f.lc_at(100.0, -3.0, 7.0, &deps), 0.5);
    }

    #[test]
    fn distance_field_computes_euclidean() {
        let f = DistanceField {
            source_type: "Point".into(),
            cx: 0.0, cy: 0.0, cz: 0.0,
            source_tag: 0,
        };
        let deps = empty_lookup();
        assert!((f.lc_at(3.0, 4.0, 0.0, &deps) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn threshold_field_smooth_step() {
        let d = DistanceField {
            source_type: "Point".into(),
            cx: 0.0, cy: 0.0, cz: 0.0,
            source_tag: 0,
        };
        let mut mgr = FieldManager::new();
        let d_tag = mgr.add(Box::new(d));
        let t_tag = mgr.add(Box::new(ThresholdField {
            in_field: d_tag,
            lc_min: 0.1,
            lc_max: 1.0,
            dist_min: 0.0,
            dist_max: 10.0,
        }));

        // At origin (dist=0) → lc_min
        let lc0 = mgr.evaluate_field(t_tag, 0.0, 0.0, 0.0);
        assert!((lc0 - 0.1).abs() < 1e-12);

        // Far away → lc_max
        let lc_far = mgr.evaluate_field(t_tag, 100.0, 0.0, 0.0);
        assert!((lc_far - 1.0).abs() < 1e-12);
    }

    #[test]
    fn min_field_returns_minimum() {
        let mut mgr = FieldManager::new();
        let a = mgr.add(Box::new(ConstantField { lc: 0.5 }));
        let b = mgr.add(Box::new(ConstantField { lc: 0.2 }));
        let m = mgr.add(Box::new(MinField { fields: vec![a, b] }));
        let lc = mgr.evaluate_field(m, 0.0, 0.0, 0.0);
        assert!((lc - 0.2).abs() < 1e-12);
    }

    #[test]
    fn max_field_returns_maximum() {
        let mut mgr = FieldManager::new();
        let a = mgr.add(Box::new(ConstantField { lc: 0.5 }));
        let b = mgr.add(Box::new(ConstantField { lc: 0.2 }));
        let m = mgr.add(Box::new(MaxField { fields: vec![a, b] }));
        let lc = mgr.evaluate_field(m, 0.0, 0.0, 0.0);
        assert!((lc - 0.5).abs() < 1e-12);
    }

    #[test]
    fn box_field_inside_constant() {
        let f = BoxField {
            lc_inside: 0.1,
            lc_outside: 1.0,
            x_min: -1.0, x_max: 1.0,
            y_min: -1.0, y_max: 1.0,
            z_min: -1.0, z_max: 1.0,
        };
        let deps = empty_lookup();
        assert!((f.lc_at(0.0, 0.0, 0.0, &deps) - 0.1).abs() < 1e-12);
        assert!((f.lc_at(10.0, 0.0, 0.0, &deps) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn field_manager_background() {
        let mut mgr = FieldManager::new();
        let tag = mgr.add(Box::new(ConstantField { lc: 0.3 }));
        mgr.set_background(tag).unwrap();
        let lc = mgr.evaluate(1.0, 2.0, 3.0).expect("should have background");
        assert!((lc - 0.3).abs() < 1e-12);
    }

    #[test]
    fn field_manager_no_background_returns_none() {
        let mgr = FieldManager::new();
        assert!(mgr.evaluate(0.0, 0.0, 0.0).is_none());
    }

    #[test]
    fn field_manager_unknown_background_errors() {
        let mut mgr = FieldManager::new();
        let err = mgr.set_background(99);
        assert!(err.is_err());
    }

    #[test]
    fn restrict_field_delegates_inside_box() {
        let mut mgr = FieldManager::new();
        let inner = mgr.add(Box::new(ConstantField { lc: 0.05 }));
        let restrict = mgr.add(Box::new(RestrictField {
            in_field: inner,
            box_bounds: Some([-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]),
        }));
        assert!((mgr.evaluate_field(restrict, 0.0, 0.0, 0.0) - 0.05).abs() < 1e-12);
        assert!(mgr.evaluate_field(restrict, 10.0, 0.0, 0.0).is_infinite());
    }

    #[test]
    fn size_as_metric_bridge_isotropic() {
        use crate::bamg_2d::MetricField2D;
        let mut mgr = FieldManager::new();
        let tag = mgr.add(Box::new(ConstantField { lc: 0.5 }));
        let bridge = SizeAsMetric2D::new(&mgr, tag);
        let m = bridge.metric_at(0.0, 0.0);
        // Isotropic metric: M = (1/lc²)·I = 4·I
        assert!((m.m11 - 4.0).abs() < 1e-12);
        assert!((m.m22 - 4.0).abs() < 1e-12);
        assert!((m.m12 - 0.0).abs() < 1e-12);
    }
}
