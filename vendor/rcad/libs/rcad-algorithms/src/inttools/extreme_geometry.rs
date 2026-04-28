//! Extreme geometry detection and handling for robust boolean operations.
//!
//! This module provides detection and specialized handling for geometric configurations
//! that are challenging for boolean operations:
//!
//! - **Near-tangent geometry**: Surfaces that are almost tangent (nearly parallel at contact)
//! - **High aspect ratio geometry**: Very long thin edges or faces
//! - **Near-degenerate geometry**: Geometry approaching degeneracy (collinear, coplanar, etc.)
//! - **Extreme size differences**: Inputs with vastly different scales

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::{BRep, Edge, Face};

use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG, AdaptiveTolerance, ToleranceLevel};

// ──────────────────────────────────────────────────────────────────────────────
// Near-Tangent Geometry
// ──────────────────────────────────────────────────────────────────────────────

/// Classification of near-tangent configuration severity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NearTangentSeverity {
    /// Configurations clearly beyond tangent threshold.
    NotTangent,
    /// Close to tangent but within normal tolerance handling.
    Marginal,
    /// Very close to tangent, requires special handling.
    NearTangent,
    /// Essentially tangent, may need fuzzy tolerance adjustment.
    Critical,
}

/// Result of near-tangent configuration detection.
#[derive(Debug, Clone)]
pub struct NearTangentConfig {
    /// Point where near-tangency was detected.
    pub point: DVec3,
    /// Normal of surface A at the contact point.
    pub normal_a: DVec3,
    /// Normal of surface B at the contact point.
    pub normal_b: DVec3,
    /// Angle between normals (radians).
    pub angle: f64,
    /// Severity classification.
    pub severity: NearTangentSeverity,
    /// Suggested fuzzy tolerance adjustment.
    pub suggested_fuzzy_adjustment: f64,
}

/// Handler for near-tangent geometry configurations.
#[derive(Debug, Clone)]
pub struct NearTangentHandler {
    /// Base tolerance for tangency detection.
    pub base_tolerance: f64,
    /// Angular threshold for near-tangent detection (radians).
    pub angular_threshold: f64,
    /// Fuzzy tolerance multiplier for near-tangent cases.
    pub fuzzy_multiplier: f64,
    /// Maximum fuzzy tolerance cap.
    pub max_fuzzy: f64,
}

impl Default for NearTangentHandler {
    fn default() -> Self {
        Self {
            base_tolerance: TOLERANCE_ABS,
            angular_threshold: TOLERANCE_ANG.sqrt(), // ~1e-4.5 radians
            fuzzy_multiplier: 10.0,
            max_fuzzy: TOLERANCE_ABS * 1000.0,
        }
    }
}

impl NearTangentHandler {
    /// Create a new handler with custom tolerances.
    pub fn new(base_tolerance: f64, angular_threshold: f64) -> Self {
        Self {
            base_tolerance,
            angular_threshold,
            fuzzy_multiplier: 10.0,
            max_fuzzy: base_tolerance * 1000.0,
        }
    }

    /// Create handler from adaptive tolerance context.
    pub fn from_adaptive(tol: AdaptiveTolerance) -> Self {
        Self {
            base_tolerance: tol.coincidence(),
            angular_threshold: tol.angular_tolerance(ToleranceLevel::Normal),
            fuzzy_multiplier: 10.0,
            max_fuzzy: tol.max_tolerance,
        }
    }

    /// Classify near-tangent severity based on angle.
    pub fn classify_severity(&self, angle: f64) -> NearTangentSeverity {
        let ang_deg = angle.to_degrees();
        if ang_deg < 0.001 {
            NearTangentSeverity::Critical
        } else if ang_deg < 0.01 {
            NearTangentSeverity::NearTangent
        } else if ang_deg < 0.1 {
            NearTangentSeverity::Marginal
        } else {
            NearTangentSeverity::NotTangent
        }
    }

    /// Detect near-tangent configurations between two surfaces.
    pub fn detect_near_tangent_configurations(
        &self,
        s1: &Surface3,
        s2: &Surface3,
        contact_points: &[DVec3],
    ) -> Vec<NearTangentConfig> {
        let mut configs = Vec::new();

        for &point in contact_points {
            // Get surface normals at contact point
            let (n1, n2) = match self.get_surface_normals_at_point(s1, s2, point) {
                Some(normals) => normals,
                None => continue,
            };

            // Compute angle between normals
            let dot = n1.dot(n2).clamp(-1.0, 1.0);
            let angle = dot.acos();

            let severity = self.classify_severity(angle);
            if severity == NearTangentSeverity::NotTangent {
                continue;
            }

            let suggested_adjustment = self.compute_fuzzy_adjustment(angle);

            configs.push(NearTangentConfig {
                point,
                normal_a: n1,
                normal_b: n2,
                angle,
                severity,
                suggested_fuzzy_adjustment: suggested_adjustment,
            });
        }

        configs
    }

    /// Get normals from both surfaces at a point.
    fn get_surface_normals_at_point(&self, s1: &Surface3, s2: &Surface3, point: DVec3) -> Option<(DVec3, DVec3)> {
        // Project point to surface to get UV parameters
        let (u1, v1) = self.project_point_to_surface(s1, point)?;
        let (u2, v2) = self.project_point_to_surface(s2, point)?;

        let n1 = s1.normal_at(u1, v1);
        let n2 = s2.normal_at(u2, v2);

        Some((n1.normalize_or_zero(), n2.normalize_or_zero()))
    }

    /// Project a 3D point onto a surface to get UV parameters.
    fn project_point_to_surface(&self, surface: &Surface3, point: DVec3) -> Option<(f64, f64)> {
        // Simple Newton iteration for projection
        let mut u = 0.5;
        let mut v = 0.5;

        for _ in 0..20 {
            let surf_point = surface.point_at(u, v);
            let diff = point - surf_point;
            if diff.length() < self.base_tolerance {
                return Some((u, v));
            }

            // Simple gradient descent step using normal as approximation
            let normal = surface.normal_at(u, v);
            let step = diff.dot(normal);
            u += step * 0.1;
            v += step * 0.1;
        }

        Some((u, v))
    }

    /// Compute suggested fuzzy tolerance adjustment based on angle.
    pub fn compute_fuzzy_adjustment(&self, angle: f64) -> f64 {
        let ang_deg = angle.to_degrees();
        let factor = if ang_deg < 0.001 {
            1000.0
        } else if ang_deg < 0.01 {
            100.0
        } else if ang_deg < 0.1 {
            10.0
        } else {
            1.0
        };

        (self.base_tolerance * factor * self.fuzzy_multiplier).min(self.max_fuzzy)
    }

    /// Adjust tolerance based on near-tangent configurations.
    pub fn adjust_tolerance_for_tangency(&self, base_fuzzy: f64, configs: &[NearTangentConfig]) -> f64 {
        if configs.is_empty() {
            return base_fuzzy;
        }

        // Find the most severe configuration
        let max_adjustment = configs
            .iter()
            .map(|c| c.suggested_fuzzy_adjustment)
            .fold(0.0, f64::max);

        base_fuzzy.max(max_adjustment)
    }
}

/// Detect near-tangent configurations at contact points between two surfaces.
pub fn detect_near_tangent_configurations(
    s1: &Surface3,
    s2: &Surface3,
    contact_points: &[DVec3],
    tolerance: f64,
) -> Vec<NearTangentConfig> {
    let handler = NearTangentHandler::new(tolerance, TOLERANCE_ANG.sqrt());
    handler.detect_near_tangent_configurations(s1, s2, contact_points)
}

// ──────────────────────────────────────────────────────────────────────────────
// High Aspect Ratio Geometry
// ──────────────────────────────────────────────────────────────────────────────

/// Threshold for high aspect ratio detection.
pub const ASPECT_RATIO_THRESHOLD: f64 = 100.0;

/// Very high aspect ratio threshold requiring special handling.
pub const ASPECT_RATIO_VERY_HIGH: f64 = 1000.0;

/// Result of high aspect ratio edge detection.
#[derive(Debug, Clone)]
pub struct HighAspectRatioEdge {
    /// Edge index in the BRep.
    pub edge_index: usize,
    /// Chord length (end-to-end distance).
    pub chord_length: f64,
    /// Approximate arc length.
    pub arc_length: f64,
    /// Aspect ratio (arc/chord or max/min dimension).
    pub aspect_ratio: f64,
    /// Whether the edge is considered problematic.
    pub is_problematic: bool,
    /// Suggested tolerance multiplier.
    pub suggested_tolerance_multiplier: f64,
}

/// Result of high aspect ratio face detection.
#[derive(Debug, Clone)]
pub struct HighAspectRatioFace {
    /// Face index in the BRep.
    pub face_index: usize,
    /// Characteristic length (e.g., bounding box diagonal).
    pub characteristic_length: f64,
    /// Minimum dimension (e.g., shortest edge or thickness).
    pub min_dimension: f64,
    /// Aspect ratio.
    pub aspect_ratio: f64,
    /// Whether the face is considered problematic.
    pub is_problematic: bool,
    /// Suggested tolerance multiplier.
    pub suggested_tolerance_multiplier: f64,
}

/// Adaptive tolerance for high aspect ratio geometry.
#[derive(Debug, Clone)]
pub struct AspectRatioAdaptiveTolerance {
    /// Base tolerance.
    pub base_tolerance: f64,
    /// Threshold for high aspect ratio detection.
    pub aspect_ratio_threshold: f64,
    /// Maximum tolerance multiplier.
    pub max_multiplier: f64,
}

impl Default for AspectRatioAdaptiveTolerance {
    fn default() -> Self {
        Self {
            base_tolerance: TOLERANCE_ABS,
            aspect_ratio_threshold: ASPECT_RATIO_THRESHOLD,
            max_multiplier: 100.0,
        }
    }
}

impl AspectRatioAdaptiveTolerance {
    /// Create from adaptive tolerance context.
    pub fn from_adaptive(tol: AdaptiveTolerance) -> Self {
        Self {
            base_tolerance: tol.coincidence(),
            aspect_ratio_threshold: ASPECT_RATIO_THRESHOLD,
            max_multiplier: 100.0,
        }
    }

    /// Detect high aspect ratio edges in a BRep.
    pub fn detect_high_aspect_ratio_edges(&self, brep: &BRep) -> Vec<HighAspectRatioEdge> {
        let mut results = Vec::new();

        for (idx, edge) in brep.edges.iter().enumerate() {
            let info = self.analyze_edge_aspect_ratio(brep, idx, edge);
            if info.aspect_ratio >= self.aspect_ratio_threshold {
                results.push(info);
            }
        }

        results
    }

    /// Analyze aspect ratio of a single edge.
    fn analyze_edge_aspect_ratio(&self, brep: &BRep, edge_idx: usize, edge: &Edge) -> HighAspectRatioEdge {
        // Get endpoints from BRep vertices
        let p_start = brep.vertices.get(edge.start).map(|v| v.point).unwrap_or(DVec3::ZERO);
        let p_end = brep.vertices.get(edge.end).map(|v| v.point).unwrap_or(DVec3::ZERO);

        let chord_length = (p_end - p_start).length();

        // Get curve if available
        let (arc_length, aspect_ratio) = if let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) {
            if let Some(curve) = brep.geom.curves.get(curve_idx) {
                // Get curve range
                let range = brep.geom.edge_curve_range.get(edge_idx)
                    .and_then(|v| *v)
                    .unwrap_or_else(|| curve.default_domain());
                let t0 = range[0];
                let t1 = range[1];

                // Approximate arc length via sampling
                let arc_length = self.approximate_arc_length(curve, t0, t1);

                let aspect_ratio = if chord_length > self.base_tolerance {
                    arc_length / chord_length
                } else {
                    1.0
                };

                (arc_length, aspect_ratio)
            } else {
                (chord_length, 1.0)
            }
        } else {
            (chord_length, 1.0)
        };

        let is_problematic = aspect_ratio >= ASPECT_RATIO_VERY_HIGH;
        let suggested_multiplier = self.compute_tolerance_multiplier(aspect_ratio);

        HighAspectRatioEdge {
            edge_index: edge_idx,
            chord_length,
            arc_length,
            aspect_ratio,
            is_problematic,
            suggested_tolerance_multiplier: suggested_multiplier,
        }
    }

    /// Approximate arc length of a curve via numerical integration.
    fn approximate_arc_length(&self, curve: &Curve3, t0: f64, t1: f64) -> f64 {
        let n = 100;
        let dt = (t1 - t0) / n as f64;
        let mut length = 0.0;

        let mut prev_point = curve.point_at(t0);
        for i in 1..=n {
            let t = t0 + dt * i as f64;
            let point = curve.point_at(t);
            length += (point - prev_point).length();
            prev_point = point;
        }

        length
    }

    /// Compute tolerance multiplier based on aspect ratio.
    pub fn compute_tolerance_multiplier(&self, aspect_ratio: f64) -> f64 {
        if aspect_ratio <= self.aspect_ratio_threshold {
            return 1.0;
        }

        // Logarithmic scaling
        let log_ratio = (aspect_ratio / self.aspect_ratio_threshold).ln();
        let multiplier = 1.0 + log_ratio * 10.0;

        multiplier.min(self.max_multiplier)
    }

    /// Get effective tolerance for a given aspect ratio.
    pub fn effective_tolerance(&self, aspect_ratio: f64) -> f64 {
        let multiplier = self.compute_tolerance_multiplier(aspect_ratio);
        self.base_tolerance * multiplier
    }
}

/// Detect high aspect ratio edges in a BRep.
pub fn detect_high_aspect_ratio_edges(brep: &BRep, tolerance: f64) -> Vec<HighAspectRatioEdge> {
    let aat = AspectRatioAdaptiveTolerance {
        base_tolerance: tolerance,
        ..Default::default()
    };
    aat.detect_high_aspect_ratio_edges(brep)
}

// ──────────────────────────────────────────────────────────────────────────────
// Near-Degenerate Geometry
// ──────────────────────────────────────────────────────────────────────────────

/// Classification of near-degenerate geometry types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegenerateType {
    /// Edge with near-zero length.
    NearZeroLengthEdge,
    /// Edge with near-zero curvature (almost straight).
    NearZeroCurvature,
    /// Face with near-zero area.
    NearZeroAreaFace,
    /// Face with near-collinear boundary.
    NearCollinearBoundary,
    /// Surface with near-singular parameterization.
    NearSingularSurface,
    /// Vertex very close to edge (not at endpoint).
    VertexNearEdge,
}

/// Result of near-degenerate geometry detection.
#[derive(Debug, Clone)]
pub struct NearDegenerateGeometry {
    /// Type of degeneracy.
    pub degenerate_type: DegenerateType,
    /// Index of the affected element.
    pub element_index: usize,
    /// Severity measure (0 = degenerate, 1 = safe).
    pub severity: f64,
    /// Description of the issue.
    pub description: String,
    /// Suggested repair action.
    pub suggested_repair: String,
}

/// Handler for near-degenerate geometry.
#[derive(Debug, Clone)]
pub struct DegenerateGeometryHandler {
    /// Tolerance for near-zero detection.
    pub zero_tolerance: f64,
    /// Tolerance for collinearity detection.
    pub collinear_tolerance: f64,
    /// Minimum area threshold for faces.
    pub min_area: f64,
    /// Minimum edge length threshold.
    pub min_edge_length: f64,
}

impl Default for DegenerateGeometryHandler {
    fn default() -> Self {
        Self {
            zero_tolerance: TOLERANCE_ABS,
            collinear_tolerance: TOLERANCE_ANG.sqrt(),
            min_area: TOLERANCE_ABS * TOLERANCE_ABS,
            min_edge_length: TOLERANCE_ABS,
        }
    }
}

impl DegenerateGeometryHandler {
    /// Create from adaptive tolerance context.
    pub fn from_adaptive(tol: AdaptiveTolerance) -> Self {
        Self {
            zero_tolerance: tol.coincidence(),
            collinear_tolerance: tol.angular_tolerance(ToleranceLevel::Normal),
            min_area: tol.coincidence() * tol.coincidence(),
            min_edge_length: tol.coincidence(),
        }
    }

    /// Detect all near-degenerate geometry in a BRep.
    pub fn detect_near_degenerate_geometry(&self, brep: &BRep) -> Vec<NearDegenerateGeometry> {
        let mut results = Vec::new();

        // Check edges
        for (idx, edge) in brep.edges.iter().enumerate() {
            if let Some(issue) = self.check_edge_degeneracy(brep, idx, edge) {
                results.push(issue);
            }
        }

        // Check faces - iterate through solids/shells/faces
        let mut face_idx = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    if let Some(issue) = self.check_face_degeneracy(brep, face_idx, face) {
                        results.push(issue);
                    }
                    face_idx += 1;
                }
            }
        }

        results
    }

    /// Check an edge for degeneracy.
    fn check_edge_degeneracy(&self, brep: &BRep, idx: usize, edge: &Edge) -> Option<NearDegenerateGeometry> {
        // Get endpoints from BRep vertices
        let p_start = brep.vertices.get(edge.start).map(|v| v.point)?;
        let p_end = brep.vertices.get(edge.end).map(|v| v.point)?;

        let length = (p_end - p_start).length();

        // Check for near-zero length
        if length < self.min_edge_length {
            let severity = length / self.min_edge_length;
            return Some(NearDegenerateGeometry {
                degenerate_type: DegenerateType::NearZeroLengthEdge,
                element_index: idx,
                severity,
                description: format!("Edge {} has near-zero length ({:.2e})", idx, length),
                suggested_repair: "Remove edge or merge with adjacent geometry".to_string(),
            });
        }

        // Check for near-zero curvature (almost straight curve that should be a line)
        if self.is_near_zero_curvature(brep, idx) {
            return Some(NearDegenerateGeometry {
                degenerate_type: DegenerateType::NearZeroCurvature,
                element_index: idx,
                severity: 0.5,
                description: format!("Edge {} has near-zero curvature", idx),
                suggested_repair: "Consider simplifying to a line segment".to_string(),
            });
        }

        None
    }

    /// Check if edge has near-zero curvature.
    fn is_near_zero_curvature(&self, brep: &BRep, edge_idx: usize) -> bool {
        // Get curve if available
        let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) else {
            // No curve - assume it's a line (zero curvature)
            return true;
        };
        let Some(curve) = brep.geom.curves.get(curve_idx) else {
            return true;
        };

        // Get curve range
        let range = brep.geom.edge_curve_range.get(edge_idx)
            .and_then(|v| *v)
            .unwrap_or_else(|| curve.default_domain());
        let t0 = range[0];
        let t1 = range[1];

        // Sample curvature at several points
        // For now, check if curve is a line
        match curve {
            Curve3::Line(_) => true,
            _ => {
                // For non-line curves, check if they're nearly straight
                let n = 10;
                for i in 0..n {
                    let t = t0 + (t1 - t0) * i as f64 / (n - 1).max(1) as f64;
                    // Simple check: compare point distance to chord
                    let mid_t = (t0 + t1) * 0.5;
                    let mid_point = curve.point_at(mid_t);
                    let start_point = curve.point_at(t0);
                    let end_point = curve.point_at(t1);
                    let chord_mid = (start_point + end_point) * 0.5;
                    let deviation = (mid_point - chord_mid).length();
                    if deviation > self.collinear_tolerance {
                        return false;
                    }
                }
                true
            }
        }
    }

    /// Check a face for degeneracy.
    fn check_face_degeneracy(&self, _brep: &BRep, idx: usize, _face: &Face) -> Option<NearDegenerateGeometry> {
        // This would require more complex face analysis
        // For now, we return None as placeholder
        // Full implementation would compute face area, check for sliver faces, etc.
        let _ = (idx, _face);
        None
    }

    /// Generate repair recommendations for all detected issues.
    pub fn generate_repair_recommendations(&self, issues: &[NearDegenerateGeometry]) -> Vec<String> {
        issues
            .iter()
            .map(|issue| {
                format!(
                    "Element {}: {} - {}",
                    issue.element_index,
                    issue.description,
                    issue.suggested_repair
                )
            })
            .collect()
    }
}

/// Detect near-degenerate geometry in a BRep.
pub fn detect_near_degenerate_geometry(brep: &BRep, tolerance: f64) -> Vec<NearDegenerateGeometry> {
    let handler = DegenerateGeometryHandler {
        zero_tolerance: tolerance,
        ..Default::default()
    };
    handler.detect_near_degenerate_geometry(brep)
}

// ──────────────────────────────────────────────────────────────────────────────
// Extreme Size Difference
// ──────────────────────────────────────────────────────────────────────────────

/// Threshold for size ratio detection.
pub const SIZE_RATIO_THRESHOLD: f64 = 1000.0;

/// Result of size difference analysis.
#[derive(Debug, Clone)]
pub struct SizeDifferenceAnalysis {
    /// Characteristic size of shape A.
    pub size_a: f64,
    /// Characteristic size of shape B.
    pub size_b: f64,
    /// Size ratio (larger/smaller).
    pub size_ratio: f64,
    /// Whether the size difference is considered extreme.
    pub is_extreme: bool,
    /// Suggested tolerance multiplier.
    pub suggested_tolerance_multiplier: f64,
    /// Whether to use relative tolerances.
    pub use_relative_tolerances: bool,
}

/// Handler for extreme size differences between shapes.
#[derive(Debug, Clone)]
pub struct SizeDifferenceHandler {
    /// Base tolerance.
    pub base_tolerance: f64,
    /// Threshold for extreme size ratio.
    pub size_ratio_threshold: f64,
    /// Maximum tolerance multiplier.
    pub max_multiplier: f64,
}

impl Default for SizeDifferenceHandler {
    fn default() -> Self {
        Self {
            base_tolerance: TOLERANCE_ABS,
            size_ratio_threshold: SIZE_RATIO_THRESHOLD,
            max_multiplier: 1000.0,
        }
    }
}

impl SizeDifferenceHandler {
    /// Create from adaptive tolerance context.
    pub fn from_adaptive(tol: AdaptiveTolerance) -> Self {
        Self {
            base_tolerance: tol.coincidence(),
            size_ratio_threshold: SIZE_RATIO_THRESHOLD,
            max_multiplier: 100.0,
        }
    }

    /// Compute characteristic size of a BRep.
    pub fn compute_characteristic_size(&self, brep: &BRep) -> f64 {
        if brep.vertices.is_empty() {
            return 1.0;
        }

        let mut min_pt = DVec3::splat(f64::INFINITY);
        let mut max_pt = DVec3::splat(f64::NEG_INFINITY);

        for vertex in &brep.vertices {
            min_pt = min_pt.min(vertex.point);
            max_pt = max_pt.max(vertex.point);
        }

        (max_pt - min_pt).length().max(1e-10)
    }

    /// Analyze size difference between two BReps.
    pub fn analyze_size_difference(&self, a: &BRep, b: &BRep) -> SizeDifferenceAnalysis {
        let size_a = self.compute_characteristic_size(a);
        let size_b = self.compute_characteristic_size(b);

        let (larger, smaller) = if size_a >= size_b {
            (size_a, size_b)
        } else {
            (size_b, size_a)
        };

        let size_ratio = larger / smaller.max(self.base_tolerance);
        let is_extreme = size_ratio >= self.size_ratio_threshold;

        // Compute suggested tolerance multiplier
        let suggested_multiplier = if is_extreme {
            let log_ratio = (size_ratio / self.size_ratio_threshold).ln();
            (10.0 * (1.0 + log_ratio)).min(self.max_multiplier)
        } else {
            1.0
        };

        // For extreme ratios, use relative tolerances based on the smaller shape
        let use_relative_tolerances = is_extreme;

        SizeDifferenceAnalysis {
            size_a,
            size_b,
            size_ratio,
            is_extreme,
            suggested_tolerance_multiplier: suggested_multiplier,
            use_relative_tolerances,
        }
    }

    /// Get effective tolerance for the smaller shape.
    pub fn effective_tolerance_for_small_shape(&self, analysis: &SizeDifferenceAnalysis) -> f64 {
        let smaller_size = analysis.size_a.min(analysis.size_b);
        // Use a relative tolerance of 1e-7 of the smaller dimension
        (smaller_size * 1e-7).max(self.base_tolerance)
    }
}

/// Analyze size difference between two BReps.
pub fn analyze_size_difference(a: &BRep, b: &BRep, tolerance: f64) -> SizeDifferenceAnalysis {
    let handler = SizeDifferenceHandler {
        base_tolerance: tolerance,
        ..Default::default()
    };
    handler.analyze_size_difference(a, b)
}

// ──────────────────────────────────────────────────────────────────────────────
// Comprehensive Extreme Geometry Analysis
// ──────────────────────────────────────────────────────────────────────────────

/// Comprehensive analysis of extreme geometry conditions.
#[derive(Debug, Clone, Default)]
pub struct ExtremeGeometryAnalysis {
    /// Near-tangent configurations detected.
    pub near_tangent_configs: Vec<NearTangentConfig>,
    /// High aspect ratio edges detected.
    pub high_aspect_ratio_edges: Vec<HighAspectRatioEdge>,
    /// Near-degenerate geometry detected.
    pub degenerate_geometry: Vec<NearDegenerateGeometry>,
    /// Size difference analysis (for two-shape operations).
    pub size_difference: Option<SizeDifferenceAnalysis>,
    /// Overall recommended fuzzy tolerance.
    pub recommended_fuzzy_tolerance: f64,
    /// Whether any extreme geometry conditions were detected.
    pub has_extreme_geometry: bool,
    /// Summary of detected issues.
    pub issues_summary: Vec<String>,
}

/// Options for extreme geometry analysis.
#[derive(Debug, Clone)]
pub struct ExtremeGeometryAnalysisOptions {
    /// Base tolerance for detection.
    pub tolerance: f64,
    /// Whether to check for near-tangent configurations.
    pub check_near_tangent: bool,
    /// Whether to check for high aspect ratio.
    pub check_aspect_ratio: bool,
    /// Whether to check for degenerate geometry.
    pub check_degenerate: bool,
    /// Whether to check for size differences.
    pub check_size_difference: bool,
    /// Contact points for near-tangent detection (if known).
    pub contact_points: Vec<DVec3>,
    /// Surfaces for near-tangent detection (if known).
    pub surfaces: Option<(Surface3, Surface3)>,
}

impl Default for ExtremeGeometryAnalysisOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            check_near_tangent: true,
            check_aspect_ratio: true,
            check_degenerate: true,
            check_size_difference: true,
            contact_points: Vec::new(),
            surfaces: None,
        }
    }
}

/// Perform comprehensive extreme geometry analysis.
pub fn analyze_extreme_geometry(
    a: &BRep,
    b: Option<&BRep>,
    options: &ExtremeGeometryAnalysisOptions,
) -> ExtremeGeometryAnalysis {
    let mut analysis = ExtremeGeometryAnalysis::default();
    let mut max_fuzzy = options.tolerance;

    // Check for near-tangent configurations
    if options.check_near_tangent {
        if let Some((ref s1, ref s2)) = options.surfaces {
            let handler = NearTangentHandler::new(options.tolerance, TOLERANCE_ANG.sqrt());
            analysis.near_tangent_configs = handler
                .detect_near_tangent_configurations(s1, s2, &options.contact_points);

            for config in &analysis.near_tangent_configs {
                max_fuzzy = max_fuzzy.max(config.suggested_fuzzy_adjustment);
            }
        }
    }

    // Check for high aspect ratio edges
    if options.check_aspect_ratio {
        let aat = AspectRatioAdaptiveTolerance::from_adaptive(AdaptiveTolerance::from_scale(
            compute_brep_scale(a),
        ));
        analysis.high_aspect_ratio_edges = aat.detect_high_aspect_ratio_edges(a);

        for edge in &analysis.high_aspect_ratio_edges {
            max_fuzzy = max_fuzzy.max(options.tolerance * edge.suggested_tolerance_multiplier);
        }
    }

    // Check for degenerate geometry
    if options.check_degenerate {
        let handler = DegenerateGeometryHandler::from_adaptive(AdaptiveTolerance::from_scale(
            compute_brep_scale(a),
        ));
        analysis.degenerate_geometry = handler.detect_near_degenerate_geometry(a);
    }

    // Check for size differences
    if options.check_size_difference {
        if let Some(b_rep) = b {
            let handler = SizeDifferenceHandler::from_adaptive(AdaptiveTolerance::from_scale(
                compute_brep_scale(a).max(compute_brep_scale(b_rep)),
            ));
            analysis.size_difference = Some(handler.analyze_size_difference(a, b_rep));

            if let Some(ref sd) = analysis.size_difference {
                max_fuzzy = max_fuzzy.max(options.tolerance * sd.suggested_tolerance_multiplier);
            }
        }
    }

    // Determine if any extreme geometry was detected
    analysis.has_extreme_geometry = !analysis.near_tangent_configs.is_empty()
        || analysis.high_aspect_ratio_edges.iter().any(|e| e.is_problematic)
        || !analysis.degenerate_geometry.is_empty()
        || analysis.size_difference.as_ref().map_or(false, |sd| sd.is_extreme);

    // Build issues summary
    if !analysis.near_tangent_configs.is_empty() {
        analysis.issues_summary.push(format!(
            "Near-tangent configurations: {}",
            analysis.near_tangent_configs.len()
        ));
    }
    if !analysis.high_aspect_ratio_edges.is_empty() {
        analysis.issues_summary.push(format!(
            "High aspect ratio edges: {}",
            analysis.high_aspect_ratio_edges.len()
        ));
    }
    if !analysis.degenerate_geometry.is_empty() {
        analysis.issues_summary.push(format!(
            "Near-degenerate geometry: {}",
            analysis.degenerate_geometry.len()
        ));
    }
    if let Some(ref sd) = analysis.size_difference {
        if sd.is_extreme {
            analysis.issues_summary.push(format!(
                "Extreme size ratio: {:.1}",
                sd.size_ratio
            ));
        }
    }

    analysis.recommended_fuzzy_tolerance = max_fuzzy;
    analysis
}

/// Compute characteristic scale of a BRep.
fn compute_brep_scale(brep: &BRep) -> f64 {
    if brep.vertices.is_empty() {
        return 1.0;
    }

    let mut min_pt = DVec3::splat(f64::INFINITY);
    let mut max_pt = DVec3::splat(f64::NEG_INFINITY);

    for vertex in &brep.vertices {
        min_pt = min_pt.min(vertex.point);
        max_pt = max_pt.max(vertex.point);
    }

    (max_pt - min_pt).length().max(1e-10)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_near_tangent_severity_classification() {
        let handler = NearTangentHandler::default();

        // 0.0 degrees -> Critical
        assert_eq!(
            handler.classify_severity(0.0),
            NearTangentSeverity::Critical
        );
        // 0.001 degrees -> NearTangent (just above Critical threshold)
        assert_eq!(
            handler.classify_severity(0.001_f64.to_radians()),
            NearTangentSeverity::NearTangent
        );
        // 0.05 degrees -> Marginal
        assert_eq!(
            handler.classify_severity(0.05_f64.to_radians()),
            NearTangentSeverity::Marginal
        );
        // 1.0 degrees -> NotTangent
        assert_eq!(
            handler.classify_severity(1.0_f64.to_radians()),
            NearTangentSeverity::NotTangent
        );
    }

    #[test]
    fn test_aspect_ratio_tolerance_multiplier() {
        let aat = AspectRatioAdaptiveTolerance::default();

        // Normal aspect ratio
        assert_eq!(aat.compute_tolerance_multiplier(10.0), 1.0);

        // High aspect ratio
        let mult = aat.compute_tolerance_multiplier(ASPECT_RATIO_THRESHOLD * 10.0);
        assert!(mult > 1.0);
        assert!(mult < aat.max_multiplier);
    }

    #[test]
    fn test_size_difference_analysis() {
        let handler = SizeDifferenceHandler::default();

        // Create two simple BReps with different sizes
        // This is a placeholder test - full test would need actual BReps
        assert!(handler.size_ratio_threshold > 0.0);
    }

    #[test]
    fn test_degenerate_geometry_handler_default() {
        let handler = DegenerateGeometryHandler::default();
        assert!(handler.zero_tolerance > 0.0);
        assert!(handler.collinear_tolerance > 0.0);
    }
}
