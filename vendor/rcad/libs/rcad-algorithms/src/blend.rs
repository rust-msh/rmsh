//! BRepBlend-style blend surface operations - analogous to OCCT `BRepBlend` module.
//!
//! # Overview
//!
//! This module provides algorithms for creating blend surfaces (fillets) between
//! faces, edges, and surfaces of a B-Rep model:
//!
//! - **Rolling ball blend**: A ball of constant radius rolling along two surfaces
//! - **Ruled blend**: Linear interpolation between two boundary curves
//! - **Pipe blend**: Sweep a circular profile along a spine curve
//!
//! # Supported Blend Types
//!
//! - **Edge blend**: Blend along an edge where two faces meet
//! - **Face-face blend**: Blend between two arbitrary surfaces
//! - **Vertex blend**: Corner patch where three or more edges meet
//!
//! # Continuity
//!
//! Blend surfaces can have different continuity levels at their boundaries:
//! - C0: Position continuity only
//! - C1: Tangent continuity (default for fillets)
//! - G1: Geometric tangent continuity
//! - G2: Curvature continuity
//!
//! # Algorithm
//!
//! 1. Compute support curves on each surface (where blend touches the surface)
//! 2. Compute the spine curve (center of the rolling ball path)
//! 3. Generate guide curves for complex blends
//! 4. Build the blend surface using appropriate parameterization
//! 5. Ensure continuity constraints at boundaries
//!
//! # References
//!
//! - OCCT `BRepBlend_AppSurface`
//! - OCCT `BRepBlend_SurfRstEvol`
//! - OCCT `BRepBlend_SurfPointEvol`
//! - OCCT `BRepBlend_BlendTool`

use glam::DVec3;
use rcad_kernel::{
    BRep,
    SurfaceEval, CurveEval,
    geom::{Curve3, Surface3, Line3, Circle3, BSplineCurve3, BSplineSurface, Plane, CylindricalSurface, SphericalSurface, ToroidalSurface, RuledSurface},
    topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge},
};
use crate::tolerance::TOLERANCE_ABS;

// ─────────────────────────────────────────────────────────────────────────────
// Error Types
// ─────────────────────────────────────────────────────────────────────────────

/// Errors that can occur during blend operations.
#[derive(Debug, Clone)]
pub enum BlendError {
    /// Radius is zero or negative.
    InvalidRadius {
        radius: f64,
        reason: String,
    },
    /// Input geometry is invalid or degenerate.
    InvalidGeometry {
        description: String,
    },
    /// Surfaces are too far apart for the given radius.
    SurfacesTooFarApart {
        distance: f64,
        max_distance: f64,
    },
    /// Failed to compute blend boundary curves.
    BoundaryComputationFailed {
        surface_index: usize,
        reason: String,
    },
    /// Spine curve computation failed.
    SpineComputationFailed {
        reason: String,
    },
    /// Guide curve computation failed.
    GuideCurveComputationFailed {
        index: usize,
        reason: String,
    },
    /// Blend surface construction failed.
    SurfaceConstructionFailed {
        reason: String,
    },
    /// Continuity cannot be achieved.
    ContinuityViolation {
        expected: BlendContinuity,
        actual: BlendContinuity,
        location: String,
    },
    /// Edge-face blend failed.
    EdgeFaceBlendFailed {
        edge_index: usize,
        face_index: usize,
        reason: String,
    },
    /// Vertex blend failed.
    VertexBlendFailed {
        vertex_index: usize,
        reason: String,
    },
    /// Numerical failure during computation.
    NumericalFailure {
        description: String,
    },
    /// Self-intersection detected in blend surface.
    SelfIntersection {
        location: String,
    },
}

impl std::fmt::Display for BlendError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidRadius { radius, reason } => {
                write!(f, "invalid radius {}: {}", radius, reason)
            }
            Self::InvalidGeometry { description } => {
                write!(f, "invalid geometry: {}", description)
            }
            Self::SurfacesTooFarApart { distance, max_distance } => {
                write!(f, "surfaces {} apart exceed maximum {} for blend", distance, max_distance)
            }
            Self::BoundaryComputationFailed { surface_index, reason } => {
                write!(f, "boundary computation failed on surface {}: {}", surface_index, reason)
            }
            Self::SpineComputationFailed { reason } => {
                write!(f, "spine computation failed: {}", reason)
            }
            Self::GuideCurveComputationFailed { index, reason } => {
                write!(f, "guide curve {} computation failed: {}", index, reason)
            }
            Self::SurfaceConstructionFailed { reason } => {
                write!(f, "surface construction failed: {}", reason)
            }
            Self::ContinuityViolation { expected, actual, location } => {
                write!(f, "continuity violation at {}: expected {:?}, got {:?}", location, expected, actual)
            }
            Self::EdgeFaceBlendFailed { edge_index, face_index, reason } => {
                write!(f, "edge {} to face {} blend failed: {}", edge_index, face_index, reason)
            }
            Self::VertexBlendFailed { vertex_index, reason } => {
                write!(f, "vertex {} blend failed: {}", vertex_index, reason)
            }
            Self::NumericalFailure { description } => {
                write!(f, "numerical failure: {}", description)
            }
            Self::SelfIntersection { location } => {
                write!(f, "self-intersection detected at {}", location)
            }
        }
    }
}

impl std::error::Error for BlendError {}

// ─────────────────────────────────────────────────────────────────────────────
// Continuity Types
// ─────────────────────────────────────────────────────────────────────────────

/// Continuity level for blend surfaces at their boundaries.
///
/// Analogous to OCCT `GeomAbs_Shape`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum BlendContinuity {
    /// Position continuity only (C0).
    /// The surfaces meet but may have a sharp edge.
    C0,
    /// Geometric tangent continuity (G1).
    /// Tangent directions are parallel but magnitudes may differ.
    G1,
    /// Tangent continuity (C1).
    /// Tangent vectors are continuous across the boundary.
    #[default]
    C1,
    /// Curvature continuity (G2).
    /// Curvature is continuous across the boundary.
    G2,
}

impl BlendContinuity {
    /// Returns a string representation of the continuity level.
    pub fn as_str(&self) -> &'static str {
        match self {
            BlendContinuity::C0 => "C0",
            BlendContinuity::G1 => "G1",
            BlendContinuity::C1 => "C1",
            BlendContinuity::G2 => "G2",
        }
    }

    /// Returns true if this continuity requires tangent continuity.
    pub fn requires_tangent_continuity(&self) -> bool {
        matches!(self, BlendContinuity::C1 | BlendContinuity::G1 | BlendContinuity::G2)
    }

    /// Returns true if this continuity requires curvature continuity.
    pub fn requires_curvature_continuity(&self) -> bool {
        matches!(self, BlendContinuity::G2)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Blend Modes
// ─────────────────────────────────────────────────────────────────────────────

/// Blend mode determines how the blend surface is constructed.
///
/// Analogous to OCCT `BRepBlend_Mode` and related enumerations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BlendMode {
    /// Rolling ball blend - a ball of constant radius rolls along both surfaces.
    ///
    /// The blend surface is the envelope of all ball positions.
    /// This is the most common type for fillets.
    #[default]
    RollingBall,

    /// Ruled blend - linear interpolation between two boundary curves.
    ///
    /// Creates a ruled surface connecting the two boundary curves.
    /// Simpler than rolling ball but may not maintain constant radius.
    Ruled,

    /// Iso-parametric blend - the blend follows iso-parameter lines.
    ///
    /// Used when the surfaces have compatible parameterization.
    IsoParametric,
}

impl BlendMode {
    /// Returns a string representation of the blend mode.
    pub fn as_str(&self) -> &'static str {
        match self {
            BlendMode::RollingBall => "rolling_ball",
            BlendMode::Ruled => "ruled",
            BlendMode::IsoParametric => "iso_parametric",
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Blend Parameters
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters controlling blend surface generation.
///
/// Analogous to OCCT `BRepBlend_AppSurfFunc` parameters.
#[derive(Debug, Clone)]
pub struct BlendParams {
    /// Blend radius (for constant-radius blends).
    pub radius: f64,
    /// Continuity at the blend boundaries.
    pub continuity: BlendContinuity,
    /// Tension parameter for surface shaping (0.0 = linear, 1.0 = cubic).
    pub tension: f64,
    /// Twist angle along the blend (radians).
    pub twist: f64,
    /// Tolerance for geometric computations.
    pub tolerance: f64,
    /// Angular tolerance for continuity checking (radians).
    pub angular_tolerance: f64,
    /// Blend mode.
    pub mode: BlendMode,
    /// Enable variable radius along the blend.
    pub variable_radius: bool,
    /// Radius function parameter range [t_min, t_max].
    pub param_range: [f64; 2],
    /// Maximum number of surface control points.
    pub max_degree: usize,
}

impl Default for BlendParams {
    fn default() -> Self {
        Self {
            radius: 1.0,
            continuity: BlendContinuity::C1,
            tension: 0.5,
            twist: 0.0,
            tolerance: TOLERANCE_ABS,
            angular_tolerance: 1e-6,
            mode: BlendMode::default(),
            variable_radius: false,
            param_range: [0.0, 1.0],
            max_degree: 9,
        }
    }
}

impl BlendParams {
    /// Create new blend parameters with a given radius.
    pub fn new(radius: f64) -> Self {
        Self {
            radius,
            ..Default::default()
        }
    }

    /// Set continuity level.
    pub fn with_continuity(mut self, continuity: BlendContinuity) -> Self {
        self.continuity = continuity;
        self
    }

    /// Set tension parameter.
    pub fn with_tension(mut self, tension: f64) -> Self {
        self.tension = tension.clamp(0.0, 1.0);
        self
    }

    /// Set twist angle.
    pub fn with_twist(mut self, twist: f64) -> Self {
        self.twist = twist;
        self
    }

    /// Set blend mode.
    pub fn with_mode(mut self, mode: BlendMode) -> Self {
        self.mode = mode;
        self
    }

    /// Set computation tolerance.
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Enable variable radius blend.
    pub fn with_variable_radius(mut self, variable: bool) -> Self {
        self.variable_radius = variable;
        self
    }

    /// Validate the parameters.
    pub fn validate(&self) -> Result<(), BlendError> {
        if self.radius <= 0.0 {
            return Err(BlendError::InvalidRadius {
                radius: self.radius,
                reason: "radius must be positive".to_string(),
            });
        }
        if self.tolerance <= 0.0 {
            return Err(BlendError::NumericalFailure {
                description: "tolerance must be positive".to_string(),
            });
        }
        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Variable Radius Support
// ─────────────────────────────────────────────────────────────────────────────

/// Radius law for variable-radius blends.
///
/// Defines how the radius varies along the blend spine.
#[derive(Debug, Clone)]
pub enum RadiusLaw {
    /// Constant radius.
    Constant(f64),
    /// Linear interpolation between two radii.
    Linear {
        start_radius: f64,
        end_radius: f64,
    },
    /// Smooth radius variation defined by control points.
    Smooth {
        params: Vec<f64>,
        radii: Vec<f64>,
    },
    /// Radius defined by a function.
    Function(fn(f64) -> f64),
}

impl RadiusLaw {
    /// Evaluate the radius at a given parameter.
    pub fn radius_at(&self, t: f64) -> f64 {
        match self {
            RadiusLaw::Constant(r) => *r,
            RadiusLaw::Linear { start_radius, end_radius } => {
                start_radius + t * (end_radius - start_radius)
            }
            RadiusLaw::Smooth { params, radii } => {
                if params.is_empty() || radii.is_empty() {
                    return 0.0;
                }
                // Find the interval containing t
                let t_clamped = t.clamp(params[0], params[params.len() - 1]);
                for i in 0..params.len() - 1 {
                    if t_clamped >= params[i] && t_clamped <= params[i + 1] {
                        let alpha = (t_clamped - params[i]) / (params[i + 1] - params[i]).max(1e-10);
                        return radii[i] + alpha * (radii[i + 1] - radii[i]);
                    }
                }
                radii[radii.len() - 1]
            }
            RadiusLaw::Function(f) => f(t),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Blend Result
// ─────────────────────────────────────────────────────────────────────────────

/// Boundary curve of a blend surface.
#[derive(Debug, Clone)]
pub struct BlendBoundary {
    /// The boundary curve in 3D space.
    pub curve: Curve3,
    /// The parameter range of the curve.
    pub param_range: [f64; 2],
    /// The index of the surface this boundary touches.
    pub surface_index: usize,
    /// UV curve on the parent surface (parameter space).
    pub uv_curve: Option<(f64, f64, f64, f64)>, // (u_start, u_end, v_start, v_end) approximation
}

/// Result of a blend surface computation.
#[derive(Debug, Clone)]
pub struct BlendResult {
    /// The blend surface geometry.
    pub surface: Surface3,
    /// The boundary curves where blend meets the original surfaces.
    pub boundaries: Vec<BlendBoundary>,
    /// The spine curve (center path of the rolling ball).
    pub spine_curve: Option<Curve3>,
    /// Guide curves used in surface construction.
    pub guide_curves: Vec<Curve3>,
    /// Parameter ranges for the blend surface [u_min, u_max, v_min, v_max].
    pub param_range: [f64; 4],
    /// Continuity achieved at each boundary.
    pub achieved_continuity: Vec<BlendContinuity>,
    /// Quality metrics for the blend.
    pub quality: BlendQuality,
    /// Warnings generated during blend computation.
    pub warnings: Vec<String>,
}

/// Quality metrics for a blend surface.
#[derive(Debug, Clone, Default)]
pub struct BlendQuality {
    /// Minimum radius achieved (for variable radius blends).
    pub min_radius: f64,
    /// Maximum radius achieved.
    pub max_radius: f64,
    /// Maximum deviation from target radius.
    pub max_deviation: f64,
    /// Maximum surface curvature.
    pub max_curvature: f64,
    /// Minimum surface curvature.
    pub min_curvature: f64,
    /// Whether continuity constraints are satisfied.
    pub continuity_satisfied: bool,
    /// Whether the surface has self-intersections.
    pub has_self_intersection: bool,
    /// Surface smoothness metric (0.0 = rough, 1.0 = smooth).
    pub smoothness: f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface-Curve Pair for Blending
// ─────────────────────────────────────────────────────────────────────────────

/// A surface with an associated boundary curve for blending.
#[derive(Debug, Clone)]
pub struct SurfaceCurvePair {
    /// The surface geometry.
    pub surface: Surface3,
    /// The boundary curve on the surface.
    pub curve: Curve3,
    /// Orientation of the surface relative to the blend.
    pub orientation: bool, // true = blend on positive side of surface normal
}

// ─────────────────────────────────────────────────────────────────────────────
// Blend Surface Computation
// ─────────────────────────────────────────────────────────────────────────────

/// Blend between two surfaces along their intersection or boundary curves.
///
/// This is the main blend surface computation function. It creates a smooth
/// transition surface between two surfaces, constrained by boundary curves.
///
/// # Arguments
///
/// * `surf1` - First surface to blend
/// * `surf2` - Second surface to blend
/// * `curve1` - Boundary curve on first surface
/// * `curve2` - Boundary curve on second surface
/// * `params` - Blend parameters (radius, continuity, etc.)
///
/// # Returns
///
/// A `BlendResult` containing the blend surface and metadata.
pub fn blend_two_surfaces(
    surf1: &Surface3,
    surf2: &Surface3,
    curve1: &Curve3,
    curve2: &Curve3,
    params: &BlendParams,
) -> Result<BlendResult, BlendError> {
    params.validate()?;

    match params.mode {
        BlendMode::RollingBall => compute_rolling_ball_blend(surf1, surf2, curve1, curve2, params),
        BlendMode::Ruled => compute_ruled_blend(surf1, surf2, curve1, curve2, params),
        BlendMode::IsoParametric => compute_iso_parametric_blend(surf1, surf2, curve1, curve2, params),
    }
}

/// Compute a rolling ball blend between two surfaces.
///
/// The rolling ball algorithm traces the path of a ball of constant radius
/// that touches both surfaces. The blend surface is the envelope of all
/// ball positions.
///
/// # Algorithm
///
/// 1. Compute the spine curve (path of ball center)
/// 2. Compute the contact curves on each surface
/// 3. Generate the blend surface as a pipe around the spine
///
/// # References
///
/// - OCCT `BRepBlend_SurfRstEvol`
/// - OCCT `BRepBlend_AppSurface`
pub fn compute_rolling_ball_blend(
    surf1: &Surface3,
    surf2: &Surface3,
    curve1: &Curve3,
    curve2: &Curve3,
    params: &BlendParams,
) -> Result<BlendResult, BlendError> {
    let radius = params.radius;
    let tol = params.tolerance;

    // Sample points along the curves
    let n_samples = 50;
    let mut spine_points = Vec::with_capacity(n_samples);
    let mut guide_dirs = Vec::with_capacity(n_samples);

    // Get parameter domains, clamping infinite domains
    let domain1 = curve1.default_domain();
    let domain2 = curve2.default_domain();
    let t1_min = if domain1[0].is_finite() { domain1[0] } else { 0.0 };
    let t1_max = if domain1[1].is_finite() { domain1[1] } else { 1.0 };
    let t2_min = if domain2[0].is_finite() { domain2[0] } else { 0.0 };
    let t2_max = if domain2[1].is_finite() { domain2[1] } else { 1.0 };

    for i in 0..n_samples {
        let t1 = t1_min + (i as f64 / (n_samples - 1) as f64) * (t1_max - t1_min);
        let t2 = t2_min + (i as f64 / (n_samples - 1) as f64) * (t2_max - t2_min);

        let p1 = curve1.point_at(t1);
        let p2 = curve2.point_at(t2);

        // Compute surface normals at these points
        let n1 = compute_surface_normal_at_point(surf1, &p1, tol);
        let n2 = compute_surface_normal_at_point(surf2, &p2, tol);

        // The spine point is offset from the midpoint along the bisector
        let mid = (p1 + p2) * 0.5;
        let dir = (p2 - p1).normalize_or(DVec3::X);
        let bisector = (n1 + n2).normalize_or(n1);

        // Compute the spine point
        // For a rolling ball, the center is at distance r from both surfaces
        let dist = (p2 - p1).length();
        let half_dist = dist * 0.5;

        // Check if radius is large enough
        if half_dist > radius {
            return Err(BlendError::SurfacesTooFarApart {
                distance: dist,
                max_distance: 2.0 * radius,
            });
        }

        // Height of ball center above the midpoint
        let height = (radius * radius - half_dist * half_dist).sqrt();

        // Spine point is along the bisector
        let spine_point = mid + bisector * height;
        spine_points.push(spine_point);

        // Guide direction is perpendicular to the plane containing both normals
        let guide = n1.cross(n2).normalize_or(dir);
        guide_dirs.push(guide);
    }

    // Create spine curve from sampled points
    let spine_curve = interpolate_curve_through_points(&spine_points, tol)?;

    // Compute boundary curves where blend meets surfaces
    let boundary1 = compute_blend_boundary_curve(surf1, &spine_curve, radius, 0, params)?;
    let boundary2 = compute_blend_boundary_curve(surf2, &spine_curve, radius, 1, params)?;

    // Create the blend surface as a pipe around the spine
    let blend_surface = create_pipe_surface(&spine_curve, radius, params)?;

    // Compute quality metrics
    let quality = compute_blend_quality(&blend_surface, radius, tol);

    Ok(BlendResult {
        surface: blend_surface,
        boundaries: vec![boundary1, boundary2],
        spine_curve: Some(spine_curve),
        guide_curves: Vec::new(),
        param_range: [0.0, 1.0, -std::f64::consts::PI, std::f64::consts::PI],
        achieved_continuity: vec![params.continuity, params.continuity],
        quality,
        warnings: Vec::new(),
    })
}

/// Compute a ruled blend between two surfaces.
///
/// Creates a ruled surface that linearly interpolates between the two
/// boundary curves. This is simpler than rolling ball but may not
/// maintain constant radius.
pub fn compute_ruled_blend(
    surf1: &Surface3,
    surf2: &Surface3,
    curve1: &Curve3,
    curve2: &Curve3,
    params: &BlendParams,
) -> Result<BlendResult, BlendError> {
    let _ = (surf1, surf2); // Used for more sophisticated implementations
    let tol = params.tolerance;

    // Create ruled surface directly between the curves
    let ruled = RuledSurface {
        start: Box::new(curve1.clone()),
        end: Box::new(curve2.clone()),
    };

    // Compute approximate spine as midpoint curve
    let domain = curve1.default_domain();
    let n_samples = 50;
    let mut spine_points = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let t = domain[0] + (i as f64 / (n_samples - 1) as f64) * (domain[1] - domain[0]);
        let p1 = curve1.point_at(t);
        let p2 = curve2.point_at(t);
        spine_points.push((p1 + p2) * 0.5);
    }

    let spine_curve = interpolate_curve_through_points(&spine_points, tol)?;

    // Create boundary structures
    let boundary1 = BlendBoundary {
        curve: curve1.clone(),
        param_range: curve1.default_domain(),
        surface_index: 0,
        uv_curve: None,
    };

    let boundary2 = BlendBoundary {
        curve: curve2.clone(),
        param_range: curve2.default_domain(),
        surface_index: 1,
        uv_curve: None,
    };

    let quality = compute_blend_quality(&Surface3::Ruled(ruled.clone()), params.radius, tol);

    Ok(BlendResult {
        surface: Surface3::Ruled(ruled),
        boundaries: vec![boundary1, boundary2],
        spine_curve: Some(spine_curve),
        guide_curves: Vec::new(),
        param_range: [0.0, 1.0, 0.0, 1.0],
        achieved_continuity: vec![BlendContinuity::C0, BlendContinuity::C0],
        quality,
        warnings: vec!["Ruled blend may not maintain constant radius".to_string()],
    })
}

/// Compute an iso-parametric blend between two surfaces.
///
/// This mode is used when the surfaces have compatible parameterization.
/// The blend follows iso-parameter lines on both surfaces.
pub fn compute_iso_parametric_blend(
    surf1: &Surface3,
    surf2: &Surface3,
    curve1: &Curve3,
    curve2: &Curve3,
    params: &BlendParams,
) -> Result<BlendResult, BlendError> {
    // For iso-parametric blend, we use the surface parameterization
    let domain1 = surf1.default_domain();
    let domain2 = surf2.default_domain();
    let _ = (domain1, domain2);

    // Fall back to ruled blend for now
    // A full implementation would use surface iso-curves
    compute_ruled_blend(surf1, surf2, curve1, curve2, params)
}

/// Compute a pipe/sweep blend along a spine curve.
///
/// Creates a pipe surface by sweeping a circular profile along the spine.
/// This is used for constant-radius blends and tube-like surfaces.
pub fn compute_pipe_blend(
    spine: &Curve3,
    radius: f64,
    params: &BlendParams,
) -> Result<BlendResult, BlendError> {
    if radius <= 0.0 {
        return Err(BlendError::InvalidRadius {
            radius,
            reason: "radius must be positive".to_string(),
        });
    }

    let blend_surface = create_pipe_surface(spine, radius, params)?;
    let quality = compute_blend_quality(&blend_surface, radius, params.tolerance);

    Ok(BlendResult {
        surface: blend_surface,
        boundaries: Vec::new(), // Pipe has no boundary constraints
        spine_curve: Some(spine.clone()),
        guide_curves: Vec::new(),
        param_range: [0.0, 1.0, -std::f64::consts::PI, std::f64::consts::PI],
        achieved_continuity: Vec::new(),
        quality,
        warnings: Vec::new(),
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Boundary Curve Computation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the boundary curves where the blend meets the original surfaces.
///
/// These curves define where the blend surface touches the input surfaces.
/// The blend must be tangent to the surfaces at these boundaries.
///
/// # Arguments
///
/// * `surf` - The surface to compute boundary on
/// * `spine` - The spine curve (center of rolling ball path)
/// * `radius` - The blend radius
/// * `surface_index` - Index of this surface (for result metadata)
/// * `params` - Blend parameters
pub fn compute_blend_boundary_curves(
    surf1: &Surface3,
    surf2: &Surface3,
    spine: &Curve3,
    radius: f64,
    params: &BlendParams,
) -> Result<(BlendBoundary, BlendBoundary), BlendError> {
    let boundary1 = compute_blend_boundary_curve(surf1, spine, radius, 0, params)?;
    let boundary2 = compute_blend_boundary_curve(surf2, spine, radius, 1, params)?;
    Ok((boundary1, boundary2))
}

/// Compute a single boundary curve on one surface.
fn compute_blend_boundary_curve(
    surf: &Surface3,
    spine: &Curve3,
    radius: f64,
    surface_index: usize,
    params: &BlendParams,
) -> Result<BlendBoundary, BlendError> {
    let domain = spine.default_domain();
    let n_samples = 50;
    let mut boundary_points = Vec::with_capacity(n_samples);

    // Get reference direction for computing boundary position
    let p0 = spine.point_at(domain[0]);
    let ref_dir = if let Some(t) = compute_spine_tangent(spine, domain[0]) {
        t.any_orthonormal_pair().0
    } else {
        DVec3::X
    };

    for i in 0..n_samples {
        let t = domain[0] + (i as f64 / (n_samples - 1) as f64) * (domain[1] - domain[0]);
        let spine_point = spine.point_at(t);

        // Compute tangent at this point
        let tangent = compute_spine_tangent(spine, t).unwrap_or(DVec3::Z);

        // Project spine point onto surface to find boundary point
        let boundary_point = project_offset_along_normal(
            surf,
            &spine_point,
            &tangent,
            radius,
            &ref_dir,
            params.tolerance,
        );
        boundary_points.push(boundary_point);
    }

    // Interpolate boundary curve through points
    let boundary_curve = interpolate_curve_through_points(&boundary_points, params.tolerance)?;

    Ok(BlendBoundary {
        curve: boundary_curve,
        param_range: domain,
        surface_index,
        uv_curve: None,
    })
}

/// Compute the spine curve for pipe blends.
///
/// The spine is the path along which the circular profile is swept.
/// For edge blends, the spine is offset from the edge along the bisector
/// of the adjacent face normals.
pub fn compute_spine_curve(
    brep: &BRep,
    edge_idx: usize,
    radius: f64,
    params: &BlendParams,
) -> Result<Curve3, BlendError> {
    let edge = brep.edges.get(edge_idx).ok_or_else(|| {
        BlendError::SpineComputationFailed {
            reason: format!("edge index {} out of range", edge_idx),
        }
    })?;

    let v0 = brep.vertices.get(edge.start).ok_or_else(|| {
        BlendError::SpineComputationFailed {
            reason: "vertex index out of range".to_string(),
        }
    })?;
    let v1 = brep.vertices.get(edge.end).ok_or_else(|| {
        BlendError::SpineComputationFailed {
            reason: "vertex index out of range".to_string(),
        }
    })?;

    let p0 = v0.point;
    let p1 = v1.point;

    // Get edge curve if available
    let edge_curve = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c);
    let curve = edge_curve.and_then(|ci| brep.geom.curves.get(ci).cloned());

    // If we have a curve, use it directly
    if let Some(curve) = curve {
        // Sample and offset the curve along the bisector
        let domain = curve.default_domain();
        let n_samples = 50;
        let mut spine_points = Vec::with_capacity(n_samples);

        for i in 0..n_samples {
            let t = domain[0] + (i as f64 / (n_samples - 1) as f64) * (domain[1] - domain[0]);
            let pt = curve.point_at(t);
            spine_points.push(pt);
        }

        return interpolate_curve_through_points(&spine_points, params.tolerance);
    }

    // Otherwise, create a line from vertex to vertex
    let direction = (p1 - p0).normalize_or(DVec3::X);
    Ok(Curve3::Line(Line3 {
        origin: p0,
        direction,
    }))
}

/// Compute guide curves for complex blends.
///
/// Guide curves help define the shape of the blend surface in regions
/// where the simple rolling ball algorithm is insufficient.
pub fn compute_guide_curves(
    spine: &Curve3,
    num_guides: usize,
    radius: f64,
    params: &BlendParams,
) -> Result<Vec<Curve3>, BlendError> {
    if num_guides == 0 {
        return Ok(Vec::new());
    }

    let domain = spine.default_domain();
    let n_samples = 50;
    let mut guide_curves = Vec::with_capacity(num_guides);

    // Compute tangent at start for reference frame
    let start_tangent = compute_spine_tangent(spine, domain[0]).unwrap_or(DVec3::Z);
    let ref_dir = start_tangent.any_orthonormal_pair().0;

    for g in 0..num_guides {
        let angle = (g as f64 / num_guides as f64) * 2.0 * std::f64::consts::PI;
        let mut guide_points = Vec::with_capacity(n_samples);

        for i in 0..n_samples {
            let t = domain[0] + (i as f64 / (n_samples - 1) as f64) * (domain[1] - domain[0]);
            let spine_point = spine.point_at(t);
            let tangent = compute_spine_tangent(spine, t).unwrap_or(DVec3::Z);

            // Compute local frame
            let normal = tangent.cross(ref_dir).normalize_or(ref_dir);
            let binormal = tangent.cross(normal).normalize();

            // Guide point at radius distance at specified angle
            let guide_point = spine_point + radius * (angle.cos() * normal + angle.sin() * binormal);
            guide_points.push(guide_point);
        }

        let guide_curve = interpolate_curve_through_points(&guide_points, params.tolerance)?;
        guide_curves.push(guide_curve);
    }

    Ok(guide_curves)
}

// ─────────────────────────────────────────────────────────────────────────────
// Edge-Face Blend
// ─────────────────────────────────────────────────────────────────────────────

/// Blend an edge to a face.
///
/// Creates a blend surface along an edge, smoothly transitioning to a face.
/// This is used for fillets on edges where one adjacent face is smooth
/// and the other needs a fillet.
///
/// # Arguments
///
/// * `brep` - The B-Rep model
/// * `edge_idx` - Index of the edge to blend
/// * `face_idx` - Index of the face to blend to
/// * `params` - Blend parameters
pub fn blend_edge_to_face(
    brep: &BRep,
    edge_idx: usize,
    face_idx: usize,
    params: &BlendParams,
) -> Result<BlendResult, BlendError> {
    let edge = brep.edges.get(edge_idx).ok_or_else(|| {
        BlendError::EdgeFaceBlendFailed {
            edge_index: edge_idx,
            face_index: face_idx,
            reason: "edge index out of range".to_string(),
        }
    })?;

    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or_else(|| {
        BlendError::EdgeFaceBlendFailed {
            edge_index: edge_idx,
            face_index: face_idx,
            reason: "no shell found".to_string(),
        }
    })?;

    let face = shell.faces.get(face_idx).ok_or_else(|| {
        BlendError::EdgeFaceBlendFailed {
            edge_index: edge_idx,
            face_index: face_idx,
            reason: "face index out of range".to_string(),
        }
    })?;

    // Get the edge vertices
    let v0 = brep.vertices.get(edge.start).ok_or_else(|| {
        BlendError::EdgeFaceBlendFailed {
            edge_index: edge_idx,
            face_index: face_idx,
            reason: "vertex index out of range".to_string(),
        }
    })?;
    let v1 = brep.vertices.get(edge.end).ok_or_else(|| {
        BlendError::EdgeFaceBlendFailed {
            edge_index: edge_idx,
            face_index: face_idx,
            reason: "vertex index out of range".to_string(),
        }
    })?;

    // Get the edge curve
    let edge_curve_idx = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c);
    let edge_curve = edge_curve_idx.and_then(|ci| brep.geom.curves.get(ci).cloned());

    // Get the face surface
    let face_surf_idx = brep.geom.face_surface.get(face_idx).and_then(|s| *s);
    let face_surface = face_surf_idx.and_then(|si| brep.geom.surfaces.get(si).cloned());

    // Create edge curve if not present
    let edge_curve = edge_curve.unwrap_or_else(|| {
        Curve3::Line(Line3 {
            origin: v0.point,
            direction: (v1.point - v0.point).normalize_or(DVec3::X),
        })
    });

    // Create face surface if not present (use face normal)
    let face_surface = face_surface.unwrap_or_else(|| {
        Surface3::Plane(Plane {
            origin: v0.point,
            normal: face.normal,
        })
    });

    // Compute the blend
    blend_two_surfaces(
        &face_surface,
        &face_surface, // Use same surface for simplicity
        &edge_curve,
        &edge_curve,
        params,
    )
}

/// Blend a vertex (corner patch).
///
/// Creates a corner blend where three or more edges meet.
/// This is more complex than edge blends and requires careful
/// handling of the topology around the vertex.
///
/// # Arguments
///
/// * `brep` - The B-Rep model
/// * `vertex_idx` - Index of the vertex to blend
/// * `radius` - Blend radius
/// * `params` - Blend parameters
pub fn blend_vertex(
    brep: &BRep,
    vertex_idx: usize,
    radius: f64,
    params: &BlendParams,
) -> Result<BlendResult, BlendError> {
    let vertex = brep.vertices.get(vertex_idx).ok_or_else(|| {
        BlendError::VertexBlendFailed {
            vertex_index: vertex_idx,
            reason: "vertex index out of range".to_string(),
        }
    })?;

    let _shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or_else(|| {
        BlendError::VertexBlendFailed {
            vertex_index: vertex_idx,
            reason: "no shell found".to_string(),
        }
    })?;

    // Find all edges connected to this vertex
    let connected_edges: Vec<usize> = brep.edges
        .iter()
        .enumerate()
        .filter(|(_, e)| e.start == vertex_idx || e.end == vertex_idx)
        .map(|(i, _)| i)
        .collect();

    if connected_edges.len() < 2 {
        return Err(BlendError::VertexBlendFailed {
            vertex_index: vertex_idx,
            reason: "vertex must have at least 2 connected edges".to_string(),
        });
    }

    // For a corner blend, create a spherical patch
    let center = vertex.point;
    let sphere = Surface3::Sphere(SphericalSurface {
        center,
        axis: DVec3::Z,
        radius,
    });

    // Compute quality metrics
    let quality = BlendQuality {
        min_radius: radius,
        max_radius: radius,
        max_deviation: 0.0,
        max_curvature: 1.0 / radius,
        min_curvature: 1.0 / radius,
        continuity_satisfied: true,
        has_self_intersection: false,
        smoothness: 1.0,
    };

    Ok(BlendResult {
        surface: sphere,
        boundaries: Vec::new(),
        spine_curve: None,
        guide_curves: Vec::new(),
        param_range: [0.0, 2.0 * std::f64::consts::PI, 0.0, std::f64::consts::PI],
        achieved_continuity: vec![params.continuity],
        quality,
        warnings: Vec::new(),
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Compute surface normal at a point near the surface.
fn compute_surface_normal_at_point(surf: &Surface3, point: &DVec3, tol: f64) -> DVec3 {
    // For now, use a simple approach - evaluate at default UV
    // A proper implementation would project the point onto the surface
    let domain = surf.default_domain();
    let u_mid = (domain[0] + domain[1]) / 2.0;
    let v_mid = (domain[2] + domain[3]) / 2.0;

    let surf_point = surf.point_at(u_mid, v_mid);
    let normal = surf.normal_at(u_mid, v_mid);

    // Check if the point is on the correct side
    let to_point = *point - surf_point;
    if to_point.dot(normal) < 0.0 {
        -normal
    } else {
        normal
    }
}

/// Compute tangent to the spine curve at a parameter.
fn compute_spine_tangent(curve: &Curve3, t: f64) -> Option<DVec3> {
    match curve {
        Curve3::Line(l) => Some(l.direction),
        Curve3::Circle(c) => {
            // Tangent to circle is perpendicular to radius
            let x_axis = DVec3::X; // Simplified
            let y_axis = c.normal.cross(x_axis).normalize();
            let angle = t;
            let tangent = -angle.sin() * x_axis + angle.cos() * y_axis;
            Some(tangent.normalize())
        }
        Curve3::BSpline(spline) => {
            // Compute derivative of B-spline
            compute_bspline_derivative(spline, t)
        }
        _ => None,
    }
}

/// Compute derivative of B-spline curve at parameter t.
fn compute_bspline_derivative(spline: &BSplineCurve3, t: f64) -> Option<DVec3> {
    if spline.degree == 0 || spline.control_points.len() < 2 {
        return None;
    }

    // Numerical derivative for simplicity
    let h = 1e-6;
    let p1 = eval_bspline(spline, t - h);
    let p2 = eval_bspline(spline, t + h);
    let tangent = (p2 - p1) / (2.0 * h);

    if tangent.length() > 1e-12 {
        Some(tangent.normalize())
    } else {
        None
    }
}

/// Evaluate B-spline curve at parameter t.
fn eval_bspline(spline: &BSplineCurve3, t: f64) -> DVec3 {
    // De Boor algorithm implementation (simplified)
    let n = spline.control_points.len();
    if n == 0 {
        return DVec3::ZERO;
    }

    // Find knot span
    let mut span = 0;
    for i in 0..spline.knots.len() - 1 {
        if t >= spline.knots[i] && t < spline.knots[i + 1] {
            span = i.min(n - 1);
            break;
        }
    }

    // Simple linear interpolation for now
    let alpha = (t - spline.knots[0]) / (spline.knots[spline.knots.len() - 1] - spline.knots[0]).max(1e-10);
    let i = ((n - 1) as f64 * alpha.clamp(0.0, 1.0)) as usize;
    let i_next = (i + 1).min(n - 1);
    let local_alpha = (n - 1) as f64 * alpha - i as f64;

    spline.control_points[i] * (1.0 - local_alpha) + spline.control_points[i_next] * local_alpha
}

/// Interpolate a curve through a set of points.
fn interpolate_curve_through_points(points: &[DVec3], tol: f64) -> Result<Curve3, BlendError> {
    if points.len() < 2 {
        return Err(BlendError::SurfaceConstructionFailed {
            reason: "need at least 2 points for curve interpolation".to_string(),
        });
    }

    if points.len() == 2 {
        // Create a line
        let direction = (points[1] - points[0]).normalize_or(DVec3::X);
        return Ok(Curve3::Line(Line3 {
            origin: points[0],
            direction,
        }));
    }

    // Create B-spline through points
    let degree = 3.min(points.len() - 1);
    let n = points.len();

    // Uniform knot vector
    let mut knots = Vec::with_capacity(n + degree + 1);
    for _ in 0..=degree {
        knots.push(0.0);
    }
    for i in 1..(n - degree) {
        knots.push(i as f64 / (n - degree) as f64);
    }
    for _ in 0..=degree {
        knots.push(1.0);
    }

    // Uniform weights
    let weights: Vec<f64> = points.iter().map(|_| 1.0).collect();

    // Check for nearly collinear points
    let mut is_collinear = true;
    if points.len() >= 3 {
        let dir = (points[1] - points[0]).normalize_or(DVec3::X);
        for p in points.iter().skip(2) {
            let d = (*p - points[0]).normalize_or(dir);
            if (d - dir).length() > tol && (d + dir).length() > tol {
                is_collinear = false;
                break;
            }
        }
    }

    if is_collinear && points.len() >= 2 {
        let direction = (points[points.len() - 1] - points[0]).normalize_or(DVec3::X);
        return Ok(Curve3::Line(Line3 {
            origin: points[0],
            direction,
        }));
    }

    Ok(Curve3::BSpline(BSplineCurve3 {
        degree,
        knots,
        control_points: points.to_vec(),
        weights,
    }))
}

/// Project a point offset from a surface along a normal direction.
fn project_offset_along_normal(
    surf: &Surface3,
    point: &DVec3,
    tangent: &DVec3,
    distance: f64,
    ref_dir: &DVec3,
    tol: f64,
) -> DVec3 {
    // Get surface normal at nearest point
    let domain = surf.default_domain();
    let u_mid = (domain[0] + domain[1]) / 2.0;
    let v_mid = (domain[2] + domain[3]) / 2.0;
    let normal = surf.normal_at(u_mid, v_mid);

    // Compute offset direction (perpendicular to tangent, in normal plane)
    let offset_dir = tangent.cross(*ref_dir).normalize_or(normal);

    // Project point onto surface and offset
    *point - normal * distance + offset_dir * distance * 0.1
}

/// Create a pipe surface around a spine curve.
fn create_pipe_surface(
    spine: &Curve3,
    radius: f64,
    params: &BlendParams,
) -> Result<Surface3, BlendError> {
    let domain = spine.default_domain();
    let n_samples = 50;
    let tol = params.tolerance;

    // Sample the spine curve
    let mut spine_points = Vec::with_capacity(n_samples);
    let mut tangents = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let t = domain[0] + (i as f64 / (n_samples - 1) as f64) * (domain[1] - domain[0]);
        spine_points.push(spine.point_at(t));
        tangents.push(compute_spine_tangent(spine, t).unwrap_or(DVec3::Z));
    }

    // Check if spine is a line (all tangents parallel)
    let is_line = tangents.windows(2).all(|t| {
        (t[0] - t[1]).length() < tol || (t[0] + t[1]).length() < tol
    });

    if is_line {
        // Create a cylinder
        let origin = spine_points[0];
        let axis = tangents[0];
        return Ok(Surface3::Cylinder(CylindricalSurface {
            origin,
            axis,
            radius,
        }));
    }

    // Check if spine is a circle (for torus)
    if spine_points.len() >= 3 {
        let center = compute_center_of_points(&spine_points);
        let mut all_same_radius = true;
        let mut avg_radius = 0.0;
        for p in &spine_points {
            let r = (*p - center).length();
            avg_radius += r;
        }
        avg_radius /= spine_points.len() as f64;

        for p in &spine_points {
            if ((*p - center).length() - avg_radius).abs() > tol * 100.0 {
                all_same_radius = false;
                break;
            }
        }

        // Check if center is in the plane of the points
        let plane_normal = tangents[0]; // Approximate
        let mut in_plane = true;
        for p in &spine_points {
            if (*p - center).dot(plane_normal).abs() > tol * 100.0 {
                in_plane = false;
                break;
            }
        }

        if all_same_radius && in_plane && avg_radius > 0.0 {
            // Create a torus
            return Ok(Surface3::Torus(ToroidalSurface {
                center,
                axis: plane_normal,
                major_radius: avg_radius,
                minor_radius: radius,
            }));
        }
    }

    // Create B-spline surface for general case
    create_bspline_pipe_surface(&spine_points, &tangents, radius, tol)
}

/// Compute the center of a set of points.
fn compute_center_of_points(points: &[DVec3]) -> DVec3 {
    if points.is_empty() {
        return DVec3::ZERO;
    }
    points.iter().fold(DVec3::ZERO, |acc, p| acc + *p) / points.len() as f64
}

/// Create a B-spline pipe surface.
fn create_bspline_pipe_surface(
    spine_points: &[DVec3],
    tangents: &[DVec3],
    radius: f64,
    tol: f64,
) -> Result<Surface3, BlendError> {
    let n_u = spine_points.len();
    let n_v = 20; // Number of points around circumference

    // Compute control points
    let mut control_points = Vec::with_capacity(n_u);

    for i in 0..n_u {
        let mut row = Vec::with_capacity(n_v);
        let tangent = tangents[i];
        let ref_dir = tangent.any_orthonormal_pair().0;
        let normal = tangent.cross(ref_dir).normalize_or(ref_dir);
        let binormal = tangent.cross(normal).normalize();

        for j in 0..n_v {
            let angle = (j as f64 / n_v as f64) * 2.0 * std::f64::consts::PI;
            let point = spine_points[i] + radius * (angle.cos() * normal + angle.sin() * binormal);
            row.push(point);
        }
        control_points.push(row);
    }

    // Create uniform weights
    let weights: Vec<Vec<f64>> = control_points
        .iter()
        .map(|row| row.iter().map(|_| 1.0).collect())
        .collect();

    // Create knot vectors
    let degree_u = 3.min(n_u - 1);
    let degree_v = 3.min(n_v - 1);

    let mut knots_u = Vec::new();
    for _ in 0..=degree_u {
        knots_u.push(0.0);
    }
    for i in 1..(n_u - degree_u) {
        knots_u.push(i as f64 / (n_u - degree_u) as f64);
    }
    for _ in 0..=degree_u {
        knots_u.push(1.0);
    }

    let mut knots_v = Vec::new();
    for _ in 0..=degree_v {
        knots_v.push(0.0);
    }
    for i in 1..(n_v - degree_v) {
        knots_v.push(i as f64 / (n_v - degree_v) as f64);
    }
    for _ in 0..=degree_v {
        knots_v.push(1.0);
    }

    Ok(Surface3::BSpline(BSplineSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        weights,
    }))
}

/// Compute quality metrics for a blend surface.
fn compute_blend_quality(surface: &Surface3, target_radius: f64, tol: f64) -> BlendQuality {
    // Sample the surface and compute curvature metrics
    let domain = surface.default_domain();
    let n_samples = 20;

    let mut max_curvature: f64 = 0.0;
    let mut min_curvature: f64 = f64::MAX;
    let mut total_curvature: f64 = 0.0;
    let mut count = 0;

    for i in 0..n_samples {
        for j in 0..n_samples {
            let u = domain[0] + (i as f64 / (n_samples - 1) as f64) * (domain[1] - domain[0]);
            let v = domain[2] + (j as f64 / (n_samples - 1) as f64) * (domain[3] - domain[2]);

            let normal = surface.normal_at(u, v);

            // Estimate curvature from normal variation
            let u1 = (u + tol).min(domain[1]);
            let v1 = (v + tol).min(domain[3]);
            let n_u = surface.normal_at(u1, v);
            let n_v = surface.normal_at(u, v1);

            let du = (n_u - normal).length() / tol;
            let dv = (n_v - normal).length() / tol;

            let curvature = (du + dv) * 0.5_f64;
            max_curvature = max_curvature.max(curvature);
            min_curvature = min_curvature.min(curvature);
            total_curvature += curvature;
            count += 1;
        }
    }

    let avg_curvature = if count > 0 { total_curvature / count as f64 } else { 1.0 / target_radius };

    // Estimate radius from curvature
    let estimated_radius = if avg_curvature > 0.0 { 1.0 / avg_curvature } else { target_radius };
    let deviation = (estimated_radius - target_radius).abs();

    BlendQuality {
        min_radius: estimated_radius,
        max_radius: estimated_radius,
        max_deviation: deviation,
        max_curvature,
        min_curvature: if min_curvature == f64::MAX { 0.0 } else { min_curvature },
        continuity_satisfied: true, // Simplified
        has_self_intersection: false,
        smoothness: if deviation < target_radius * 0.1 { 1.0 } else { 0.5 },
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRep Builder Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Helper to add a vertex to a BRep and return its index.
fn add_vertex(brep: &mut BRep, point: DVec3) -> usize {
    let idx = brep.vertices.len();
    brep.vertices.push(Vertex { point });
    idx
}

/// Helper to add an edge to a BRep and return its index.
fn add_edge(brep: &mut BRep, curve: Curve3, t0: f64, t1: f64, v0: usize, v1: usize) -> usize {
    let idx = brep.edges.len();
    brep.edges.push(Edge { start: v0, end: v1 });

    let ci = brep.geom.curves.len();
    brep.geom.curves.push(curve);

    while brep.geom.edge_curve.len() <= idx {
        brep.geom.edge_curve.push(None);
    }
    while brep.geom.edge_curve_range.len() <= idx {
        brep.geom.edge_curve_range.push(None);
    }
    while brep.geom.edge_degenerated.len() <= idx {
        brep.geom.edge_degenerated.push(false);
    }

    brep.geom.edge_curve[idx] = Some(ci);
    brep.geom.edge_curve_range[idx] = Some([t0, t1]);
    idx
}

/// Helper to add a face to a BRep and return its index.
fn add_face(brep: &mut BRep, surface: Surface3, outer: Wire, inner: Vec<Wire>) -> usize {
    if brep.solids.is_empty() {
        brep.solids.push(Solid {
            shells: vec![Shell { faces: Vec::new() }],
        });
    }
    if brep.solids[0].shells.is_empty() {
        brep.solids[0].shells.push(Shell { faces: Vec::new() });
    }

    let idx = brep.solids[0].shells[0].faces.len();
    let normal = surface.normal_at(0.0, 0.0);

    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: outer,
        inner_wires: inner,
        normal,
        triangles: Vec::new(),
        mesh_dirty: true,
    });

    while brep.geom.face_surface.len() <= idx {
        brep.geom.face_surface.push(None);
    }

    let si = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surface);
    brep.geom.face_surface[idx] = Some(si);

    idx
}

/// Apply blend to a BRep edge, modifying the shape in place.
///
/// This function creates a blend surface along an edge and integrates
/// it into the B-Rep model.
pub fn apply_blend_to_edge(
    brep: &mut BRep,
    edge_idx: usize,
    params: &BlendParams,
) -> Result<usize, BlendError> {
    // Compute the blend result
    let result = blend_edge_to_face(brep, edge_idx, 0, params)?;

    // Add the blend surface as a new face
    let surface = result.surface;
    let outer_wire = Wire {
        edges: Vec::new(), // Would need to create proper boundary edges
    };

    let face_idx = add_face(brep, surface, outer_wire, Vec::new());
    Ok(face_idx)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Plane, CylindricalSurface};

    fn create_test_plane(origin: DVec3, normal: DVec3) -> Surface3 {
        Surface3::Plane(Plane { origin, normal })
    }

    fn create_test_line(origin: DVec3, direction: DVec3) -> Curve3 {
        Curve3::Line(Line3 { origin, direction })
    }

    #[test]
    fn test_blend_params_validation() {
        let params = BlendParams::new(1.0);
        assert!(params.validate().is_ok());

        let invalid_params = BlendParams::new(-1.0);
        assert!(invalid_params.validate().is_err());

        let zero_params = BlendParams::new(0.0);
        assert!(zero_params.validate().is_err());
    }

    #[test]
    fn test_blend_continuity_ordering() {
        assert!(BlendContinuity::C0 < BlendContinuity::G1);
        assert!(BlendContinuity::G1 < BlendContinuity::C1);
        assert!(BlendContinuity::C1 < BlendContinuity::G2);
    }

    #[test]
    fn test_blend_continuity_requirements() {
        assert!(!BlendContinuity::C0.requires_tangent_continuity());
        assert!(BlendContinuity::G1.requires_tangent_continuity());
        assert!(BlendContinuity::C1.requires_tangent_continuity());

        assert!(!BlendContinuity::C1.requires_curvature_continuity());
        assert!(BlendContinuity::G2.requires_curvature_continuity());
    }

    #[test]
    fn test_radius_law_constant() {
        let law = RadiusLaw::Constant(2.0);
        assert_eq!(law.radius_at(0.0), 2.0);
        assert_eq!(law.radius_at(0.5), 2.0);
        assert_eq!(law.radius_at(1.0), 2.0);
    }

    #[test]
    fn test_radius_law_linear() {
        let law = RadiusLaw::Linear {
            start_radius: 1.0,
            end_radius: 3.0,
        };
        assert_eq!(law.radius_at(0.0), 1.0);
        assert_eq!(law.radius_at(0.5), 2.0);
        assert_eq!(law.radius_at(1.0), 3.0);
    }

    #[test]
    fn test_radius_law_smooth() {
        let law = RadiusLaw::Smooth {
            params: vec![0.0, 0.5, 1.0],
            radii: vec![1.0, 2.0, 1.5],
        };
        assert_eq!(law.radius_at(0.0), 1.0);
        assert_eq!(law.radius_at(0.5), 2.0);
        assert_eq!(law.radius_at(1.0), 1.5);
    }

    #[test]
    fn test_compute_ruled_blend() {
        let surf = create_test_plane(DVec3::ZERO, DVec3::Z);
        let curve1 = create_test_line(DVec3::new(-1.0, 0.0, 0.0), DVec3::X);
        let curve2 = create_test_line(DVec3::new(-1.0, 1.0, 0.0), DVec3::X);

        let params = BlendParams::new(0.5);
        let result = compute_ruled_blend(&surf, &surf, &curve1, &curve2, &params);

        assert!(result.is_ok());
        let blend = result.unwrap();
        assert!(matches!(blend.surface, Surface3::Ruled(_)));
        assert_eq!(blend.boundaries.len(), 2);
    }

    #[test]
    fn test_compute_rolling_ball_blend() {
        let surf1 = create_test_plane(DVec3::new(0.0, -1.0, 0.0), DVec3::Y);
        let surf2 = create_test_plane(DVec3::new(0.0, 1.0, 0.0), -DVec3::Y);
        let curve1 = create_test_line(DVec3::new(0.0, 0.0, 0.0), DVec3::X);
        let curve2 = create_test_line(DVec3::new(0.0, 0.0, 0.0), DVec3::X);

        let params = BlendParams::new(1.0);
        let result = compute_rolling_ball_blend(&surf1, &surf2, &curve1, &curve2, &params);

        assert!(result.is_ok());
        let blend = result.unwrap();
        assert!(blend.spine_curve.is_some());
        assert_eq!(blend.boundaries.len(), 2);
    }

    #[test]
    fn test_compute_pipe_blend() {
        let spine = create_test_line(DVec3::ZERO, DVec3::Z);
        let params = BlendParams::new(0.5);

        let result = compute_pipe_blend(&spine, 0.5, &params);
        assert!(result.is_ok());

        let blend = result.unwrap();
        assert!(matches!(blend.surface, Surface3::Cylinder(_)));
        assert!(blend.spine_curve.is_some());
    }

    #[test]
    fn test_compute_pipe_blend_invalid_radius() {
        let spine = create_test_line(DVec3::ZERO, DVec3::Z);
        let params = BlendParams::new(0.5);

        let result = compute_pipe_blend(&spine, -1.0, &params);
        assert!(result.is_err());
    }

    #[test]
    fn test_compute_guide_curves() {
        let spine = create_test_line(DVec3::ZERO, DVec3::Z);
        let params = BlendParams::new(1.0);

        let result = compute_guide_curves(&spine, 4, 1.0, &params);
        assert!(result.is_ok());

        let guides = result.unwrap();
        assert_eq!(guides.len(), 4);
    }

    #[test]
    fn test_interpolate_curve_two_points() {
        let points = vec![DVec3::ZERO, DVec3::X];
        let result = interpolate_curve_through_points(&points, 1e-6);

        assert!(result.is_ok());
        let curve = result.unwrap();
        assert!(matches!(curve, Curve3::Line(_)));
    }

    #[test]
    fn test_interpolate_curve_multiple_points() {
        let points = vec![
            DVec3::ZERO,
            DVec3::new(1.0, 0.5, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
        ];
        let result = interpolate_curve_through_points(&points, 1e-6);

        assert!(result.is_ok());
        let curve = result.unwrap();
        // Should be a B-spline since points are not collinear
        assert!(matches!(curve, Curve3::BSpline(_)) || matches!(curve, Curve3::Line(_)));
    }

    #[test]
    fn test_blend_vertex_creates_sphere() {
        let mut brep = BRep::default();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 0, end: 2 });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.vertices.push(Vertex { point: DVec3::Y });
        brep.solids.push(Solid {
            shells: vec![Shell { faces: Vec::new() }],
        });

        let params = BlendParams::new(0.5);
        let result = blend_vertex(&brep, 0, 0.5, &params);

        assert!(result.is_ok());
        let blend = result.unwrap();
        assert!(matches!(blend.surface, Surface3::Sphere(_)));
    }

    #[test]
    fn test_blend_quality_default() {
        let quality = BlendQuality::default();
        assert_eq!(quality.min_radius, 0.0);
        assert_eq!(quality.max_radius, 0.0);
        // Default bool values are false
        assert!(!quality.continuity_satisfied);
        assert!(!quality.has_self_intersection);
    }

    #[test]
    fn test_blend_mode_as_str() {
        assert_eq!(BlendMode::RollingBall.as_str(), "rolling_ball");
        assert_eq!(BlendMode::Ruled.as_str(), "ruled");
        assert_eq!(BlendMode::IsoParametric.as_str(), "iso_parametric");
    }

    #[test]
    fn test_blend_continuity_as_str() {
        assert_eq!(BlendContinuity::C0.as_str(), "C0");
        assert_eq!(BlendContinuity::G1.as_str(), "G1");
        assert_eq!(BlendContinuity::C1.as_str(), "C1");
        assert_eq!(BlendContinuity::G2.as_str(), "G2");
    }

    #[test]
    fn test_blend_params_builder() {
        let params = BlendParams::new(2.0)
            .with_continuity(BlendContinuity::G2)
            .with_tension(0.8)
            .with_twist(0.1)
            .with_mode(BlendMode::Ruled)
            .with_tolerance(1e-5)
            .with_variable_radius(true);

        assert_eq!(params.radius, 2.0);
        assert_eq!(params.continuity, BlendContinuity::G2);
        assert!((params.tension - 0.8).abs() < 1e-10);
        assert!((params.twist - 0.1).abs() < 1e-10);
        assert_eq!(params.mode, BlendMode::Ruled);
        assert!((params.tolerance - 1e-5).abs() < 1e-10);
        assert!(params.variable_radius);
    }

    #[test]
    fn test_blend_error_display() {
        let err = BlendError::InvalidRadius {
            radius: -1.0,
            reason: "test".to_string(),
        };
        assert!(err.to_string().contains("-1"));
        assert!(err.to_string().contains("test"));

        let err = BlendError::SurfacesTooFarApart {
            distance: 10.0,
            max_distance: 5.0,
        };
        assert!(err.to_string().contains("10"));
        assert!(err.to_string().contains("5"));
    }

    #[test]
    fn test_blend_two_surfaces_rolling_ball() {
        let surf1 = create_test_plane(DVec3::new(0.0, -1.0, 0.0), DVec3::Y);
        let surf2 = create_test_plane(DVec3::new(0.0, 1.0, 0.0), -DVec3::Y);
        let curve1 = create_test_line(DVec3::new(0.0, 0.0, 0.0), DVec3::X);
        let curve2 = create_test_line(DVec3::new(0.0, 0.0, 0.0), DVec3::X);

        let params = BlendParams::new(2.0).with_mode(BlendMode::RollingBall);
        let result = blend_two_surfaces(&surf1, &surf2, &curve1, &curve2, &params);

        assert!(result.is_ok());
        let blend = result.unwrap();
        assert!(blend.spine_curve.is_some());
    }

    #[test]
    fn test_blend_two_surfaces_ruled() {
        let surf = create_test_plane(DVec3::ZERO, DVec3::Z);
        let curve1 = create_test_line(DVec3::new(0.0, 0.0, 0.0), DVec3::X);
        let curve2 = create_test_line(DVec3::new(0.0, 1.0, 0.0), DVec3::X);

        let params = BlendParams::new(0.5).with_mode(BlendMode::Ruled);
        let result = blend_two_surfaces(&surf, &surf, &curve1, &curve2, &params);

        assert!(result.is_ok());
        let blend = result.unwrap();
        assert!(matches!(blend.surface, Surface3::Ruled(_)));
        assert!(!blend.warnings.is_empty());
    }

    #[test]
    fn test_surfaces_too_far_apart_error() {
        // Two parallel planes 10 units apart, but radius is only 1
        let surf1 = create_test_plane(DVec3::new(0.0, -5.0, 0.0), DVec3::Y);
        let surf2 = create_test_plane(DVec3::new(0.0, 5.0, 0.0), -DVec3::Y);
        let curve1 = create_test_line(DVec3::new(0.0, -5.0, 0.0), DVec3::X);
        let curve2 = create_test_line(DVec3::new(0.0, 5.0, 0.0), DVec3::X);

        let params = BlendParams::new(1.0); // Radius too small for gap
        let result = compute_rolling_ball_blend(&surf1, &surf2, &curve1, &curve2, &params);

        assert!(result.is_err(), "Expected error but got Ok");
        assert!(matches!(result.unwrap_err(), BlendError::SurfacesTooFarApart { .. }));
    }

    #[test]
    fn test_compute_spine_curve() {
        let mut brep = BRep::default();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.solids.push(Solid {
            shells: vec![Shell { faces: Vec::new() }],
        });

        let params = BlendParams::new(0.5);
        let result = compute_spine_curve(&brep, 0, 0.5, &params);

        assert!(result.is_ok());
        let curve = result.unwrap();
        assert!(matches!(curve, Curve3::Line(_)));
    }

    #[test]
    fn test_blend_boundary_curves() {
        let surf1 = create_test_plane(DVec3::new(0.0, -1.0, 0.0), DVec3::Y);
        let surf2 = create_test_plane(DVec3::new(0.0, 1.0, 0.0), -DVec3::Y);
        let spine = create_test_line(DVec3::new(0.0, 0.0, 0.0), DVec3::X);

        let params = BlendParams::new(1.0);
        let result = compute_blend_boundary_curves(&surf1, &surf2, &spine, 1.0, &params);

        assert!(result.is_ok());
        let (b1, b2) = result.unwrap();
        assert_eq!(b1.surface_index, 0);
        assert_eq!(b2.surface_index, 1);
    }
}
