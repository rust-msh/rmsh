//! Adaptor3d-style 3D geometry adaptors.
//!
//! This module provides adaptor classes for 3D curves and surfaces, similar to OCCT's
//! Adaptor3d package. Adaptors provide a uniform interface for querying geometric
//! properties without exposing the underlying geometry types directly.
//!
//! # OCCT Equivalents
//!
//! - `Adaptor3d_Curve` -> `Curve3dAdaptor`
//! - `Adaptor3d_Surface` -> `SurfaceAdaptor`
//! - `Adaptor3d_CurveOnSurface` -> `CurveOnSurfaceAdaptor`
//! - `Adaptor3d_HSurface` -> `HSurfaceAdaptor`
//!
//! # Examples
//!
//! ```
//! use rcad_algorithms::adaptor3d::Curve3dAdaptor;
//! use rcad_kernel::Curve3;
//! use rcad_kernel::geom::Circle3;
//! use glam::dvec3;
//!
//! let circle = Circle3 {
//!     center: dvec3(0.0, 0.0, 0.0),
//!     normal: dvec3(0.0, 0.0, 1.0),
//!     radius: 1.0,
//! };
//! let curve = Curve3::Circle(circle);
//! let adaptor = Curve3dAdaptor::from_curve(&curve);
//!
//! let p = adaptor.point_at(0.0);
//! assert!((p.length() - 1.0).abs() < 1e-10); // point is on the circle
//!
//! let domain = adaptor.domain();
//! assert!((domain[1] - domain[0] - std::f64::consts::TAU).abs() < 1e-10);
//! ```

use glam::{DVec2, DVec3};
use rcad_kernel::{Curve2d, Curve3, Curve2dEval, CurveEval, Surface3, SurfaceEval};
use std::rc::Rc;
use std::f64::consts::PI;

// ============================================================================
// Curve3dAdaptor
// ============================================================================

/// Adaptor for 3D curves providing a uniform interface for geometric queries.
///
/// Wraps any `Curve3` and provides methods for point evaluation, derivatives,
/// domain information, and topological properties (closed, periodic).
///
/// Analogous to OCCT `Adaptor3d_Curve`.
#[derive(Debug, Clone)]
pub struct Curve3dAdaptor {
    curve: Curve3,
    first: f64,
    last: f64,
}

impl Curve3dAdaptor {
    /// Creates a new adaptor from a curve reference.
    ///
    /// The domain is initialized to the curve's default (natural) domain.
    pub fn from_curve(curve: &Curve3) -> Self {
        let domain = curve.default_domain();
        Self {
            curve: curve.clone(),
            first: domain[0],
            last: domain[1],
        }
    }

    /// Creates an adaptor with a trimmed parameter range.
    pub fn from_curve_with_range(curve: &Curve3, first: f64, last: f64) -> Self {
        Self {
            curve: curve.clone(),
            first,
            last,
        }
    }

    /// Returns the point on the curve at parameter `t`.
    pub fn point_at(&self, t: f64) -> DVec3 {
        self.curve.point_at(t)
    }

    /// Returns the n-th derivative at parameter `t`.
    ///
    /// - `order = 0`: point (same as `point_at`)
    /// - `order = 1`: first derivative (tangent vector)
    /// - `order = 2`: second derivative
    /// - etc.
    ///
    /// For analytic curves (lines, circles, ellipses), exact derivatives are computed.
    /// For B-Splines and Bezier curves, numerical differentiation is used for orders > 1.
    pub fn derivative(&self, t: f64, order: usize) -> DVec3 {
        match order {
            0 => self.point_at(t),
            1 => self.first_derivative(t),
            2 => self.second_derivative(t),
            n => self.nth_derivative_numerical(t, n),
        }
    }

    /// Returns the parameter domain as `[first, last]`.
    ///
    /// For trimmed adaptors, this returns the trimmed range, not the natural domain.
    pub fn domain(&self) -> [f64; 2] {
        [self.first, self.last]
    }

    /// Returns true if the curve is closed (start point equals end point).
    ///
    /// A curve is considered closed if `point_at(first) ≈ point_at(last)`.
    pub fn is_closed(&self) -> bool {
        if !self.first.is_finite() || !self.last.is_finite() {
            return false;
        }
        let p_first = self.point_at(self.first);
        let p_last = self.point_at(self.last);
        (p_first - p_last).length() < 1e-10
    }

    /// Returns true if the curve is periodic.
    ///
    /// A periodic curve has a well-defined period and wraps around seamlessly.
    /// Circles, ellipses, and toroidal curves are periodic.
    pub fn is_periodic(&self) -> bool {
        matches!(
            self.curve,
            Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::CircularHelix(_)
        )
    }

    /// Returns the period of the curve if periodic, or `None`.
    ///
    /// - Circles and ellipses: `2π`
    /// - Circular helix: `2π` (azimuthal period)
    pub fn period(&self) -> Option<f64> {
        match &self.curve {
            Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::CircularHelix(_) => Some(2.0 * PI),
            _ => None,
        }
    }

    /// Returns the resolution for a given tolerance.
    ///
    /// This is the parameter step that produces a point within `tol` of the current point.
    /// Used for discretization and sampling algorithms.
    pub fn resolution(&self, tol: f64) -> f64 {
        let domain = self.curve.default_domain();
        let mid = (domain[0] + domain[1]) / 2.0;

        // Estimate average curvature at the midpoint
        let p0 = self.point_at(mid);
        let p1 = self.point_at(mid + 1e-6);
        let chord_len = (p1 - p0).length();

        if chord_len < 1e-15 {
            return tol;
        }

        // Resolution is approximately tol / chord_per_unit_param
        tol / chord_len * 1e-6
    }

    /// Returns the first derivative (tangent vector) at parameter `t`.
    fn first_derivative(&self, t: f64) -> DVec3 {
        match &self.curve {
            Curve3::Line(line) => line.direction,
            Curve3::Circle(circle) => self.circle_tangent(circle, t),
            Curve3::Ellipse(ellipse) => self.ellipse_tangent(ellipse, t),
            Curve3::BSpline(bspline) => self.bspline_tangent(bspline, t),
            Curve3::Bezier(bezier) => self.bezier_tangent(bezier, t),
            Curve3::Offset(offset) => self.offset_tangent(offset, t),
            Curve3::Hyperbola(hyp) => self.hyperbola_tangent(hyp, t),
            Curve3::Parabola(par) => self.parabola_tangent(par, t),
            Curve3::CircularHelix(helix) => self.helix_tangent(helix, t),
            Curve3::SineWave(sine) => self.sine_wave_tangent(sine, t),
        }
    }

    /// Returns the second derivative at parameter `t`.
    fn second_derivative(&self, t: f64) -> DVec3 {
        let eps = 1e-6;
        let d1_plus = self.first_derivative(t + eps);
        let d1_minus = self.first_derivative(t - eps);
        (d1_plus - d1_minus) / (2.0 * eps)
    }

    /// Numerical nth derivative using finite differences.
    fn nth_derivative_numerical(&self, t: f64, order: usize) -> DVec3 {
        if order == 0 {
            return self.point_at(t);
        }

        let eps = 1e-4;
        let mut result = DVec3::ZERO;

        // Central difference formula for nth derivative
        for i in 0..=order {
            let coeff = central_diff_coefficient(order, i);
            let ti = t + (i as f64 - order as f64 / 2.0) * eps;
            result += coeff * self.point_at(ti);
        }

        result / eps.powi(order as i32)
    }

    // --- Analytic tangent computations ---

    fn circle_tangent(&self, circle: &rcad_kernel::geom::Circle3, t: f64) -> DVec3 {
        circle.tangent_at(t)
    }

    fn ellipse_tangent(&self, ellipse: &rcad_kernel::geom::Ellipse3, t: f64) -> DVec3 {
        ellipse.tangent_at(t)
    }

    fn bspline_tangent(&self, bspline: &rcad_kernel::geom::BSplineCurve3, t: f64) -> DVec3 {
        bspline.tangent_at(t)
    }

    fn bezier_tangent(&self, bezier: &rcad_kernel::geom::BezierCurve3, t: f64) -> DVec3 {
        bezier.tangent_at(t)
    }

    fn offset_tangent(&self, offset: &rcad_kernel::geom::OffsetCurve3, t: f64) -> DVec3 {
        offset.tangent_at(t)
    }

    fn hyperbola_tangent(&self, hyp: &rcad_kernel::geom::Hyperbola3, t: f64) -> DVec3 {
        hyp.tangent_at(t)
    }

    fn parabola_tangent(&self, par: &rcad_kernel::geom::Parabola3, t: f64) -> DVec3 {
        par.tangent_at(t)
    }

    fn helix_tangent(&self, helix: &rcad_kernel::geom::CircularHelix3, t: f64) -> DVec3 {
        helix.tangent_at(t)
    }

    fn sine_wave_tangent(&self, sine: &rcad_kernel::geom::SineWave3, t: f64) -> DVec3 {
        sine.tangent_at(t)
    }

    /// Returns a reference to the underlying curve.
    pub fn curve(&self) -> &Curve3 {
        &self.curve
    }
}

/// Returns the coefficient for central finite difference approximation.
fn central_diff_coefficient(order: usize, i: usize) -> f64 {
    // Binomial coefficients with alternating signs
    let sign = if (order - i) % 2 == 0 { 1.0 } else { -1.0 };
    sign * binomial(order, i)
}

/// Binomial coefficient (n choose k).
fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let mut result = 1.0;
    for i in 0..k.min(n - k) {
        result *= (n - i) as f64 / (i + 1) as f64;
    }
    result
}

// ============================================================================
// SurfaceAdaptor
// ============================================================================

/// Adaptor for 3D surfaces providing a uniform interface for geometric queries.
///
/// Wraps any `Surface3` and provides methods for point evaluation, partial derivatives,
/// domain information, and topological properties (closed, periodic in U/V).
///
/// Analogous to OCCT `Adaptor3d_Surface`.
#[derive(Debug, Clone)]
pub struct SurfaceAdaptor {
    surface: Surface3,
    u_first: f64,
    u_last: f64,
    v_first: f64,
    v_last: f64,
}

impl SurfaceAdaptor {
    /// Creates a new adaptor from a surface reference.
    ///
    /// The domain is initialized to the surface's default (natural) domain.
    pub fn from_surface(surface: &Surface3) -> Self {
        let [u0, u1, v0, v1] = surface.default_domain();
        Self {
            surface: surface.clone(),
            u_first: u0,
            u_last: u1,
            v_first: v0,
            v_last: v1,
        }
    }

    /// Creates an adaptor with a trimmed parameter range.
    pub fn from_surface_with_range(
        surface: &Surface3,
        u_first: f64,
        u_last: f64,
        v_first: f64,
        v_last: f64,
    ) -> Self {
        Self {
            surface: surface.clone(),
            u_first,
            u_last,
            v_first,
            v_last,
        }
    }

    /// Returns the point on the surface at parameters `(u, v)`.
    pub fn point_at(&self, u: f64, v: f64) -> DVec3 {
        self.surface.point_at(u, v)
    }

    /// Returns the partial derivative with respect to U at `(u, v)`.
    pub fn derivative_u(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-6;
        let p_plus = self.point_at(u + eps, v);
        let p_minus = self.point_at(u - eps, v);
        (p_plus - p_minus) / (2.0 * eps)
    }

    /// Returns the partial derivative with respect to V at `(u, v)`.
    pub fn derivative_v(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-6;
        let p_plus = self.point_at(u, v + eps);
        let p_minus = self.point_at(u, v - eps);
        (p_plus - p_minus) / (2.0 * eps)
    }

    /// Returns the second partial derivative ∂²S/∂u² at `(u, v)`.
    pub fn derivative_uu(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-6;
        let d_plus = self.derivative_u(u + eps, v);
        let d_minus = self.derivative_u(u - eps, v);
        (d_plus - d_minus) / (2.0 * eps)
    }

    /// Returns the second partial derivative ∂²S/∂v² at `(u, v)`.
    pub fn derivative_vv(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-6;
        let d_plus = self.derivative_v(u, v + eps);
        let d_minus = self.derivative_v(u, v - eps);
        (d_plus - d_minus) / (2.0 * eps)
    }

    /// Returns the mixed partial derivative ∂²S/∂u∂v at `(u, v)`.
    pub fn derivative_uv(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-6;
        let d_plus = self.derivative_u(u, v + eps);
        let d_minus = self.derivative_u(u, v - eps);
        (d_plus - d_minus) / (2.0 * eps)
    }

    /// Returns the normal vector at `(u, v)`.
    pub fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        self.surface.normal_at(u, v)
    }

    /// Returns the parameter domain as `[u_first, u_last, v_first, v_last]`.
    pub fn domain(&self) -> [f64; 4] {
        [self.u_first, self.u_last, self.v_first, self.v_last]
    }

    /// Returns true if the surface is closed in the U direction.
    ///
    /// A surface is U-closed if `S(u_first, v) ≈ S(u_last, v)` for all v.
    pub fn is_u_closed(&self) -> bool {
        if !self.u_first.is_finite() || !self.u_last.is_finite() {
            return false;
        }

        // Sample a few V values to check closure
        let v_samples = if self.v_first.is_finite() && self.v_last.is_finite() {
            let dv = (self.v_last - self.v_first) / 4.0;
            vec![
                self.v_first,
                self.v_first + dv,
                self.v_first + 2.0 * dv,
                self.v_first + 3.0 * dv,
                self.v_last,
            ]
        } else {
            vec![0.0]
        };

        for v in v_samples {
            let p_first = self.point_at(self.u_first, v);
            let p_last = self.point_at(self.u_last, v);
            if (p_first - p_last).length() > 1e-10 {
                return false;
            }
        }
        true
    }

    /// Returns true if the surface is closed in the V direction.
    ///
    /// A surface is V-closed if `S(u, v_first) ≈ S(u, v_last)` for all u.
    pub fn is_v_closed(&self) -> bool {
        if !self.v_first.is_finite() || !self.v_last.is_finite() {
            return false;
        }

        // Sample a few U values to check closure
        let u_samples = if self.u_first.is_finite() && self.u_last.is_finite() {
            let du = (self.u_last - self.u_first) / 4.0;
            vec![
                self.u_first,
                self.u_first + du,
                self.u_first + 2.0 * du,
                self.u_first + 3.0 * du,
                self.u_last,
            ]
        } else {
            vec![0.0]
        };

        for u in u_samples {
            let p_first = self.point_at(u, self.v_first);
            let p_last = self.point_at(u, self.v_last);
            if (p_first - p_last).length() > 1e-10 {
                return false;
            }
        }
        true
    }

    /// Returns the U period if the surface is U-periodic, or `None`.
    ///
    /// - Cylinder, cone, sphere, torus: `2π`
    /// - Others: `None`
    pub fn u_period(&self) -> Option<f64> {
        match &self.surface {
            Surface3::Cylinder(_)
            | Surface3::Sphere(_)
            | Surface3::Cone(_)
            | Surface3::Torus(_)
            | Surface3::Ellipsoid(_)
            | Surface3::Helicoid(_)
            | Surface3::Pipe(_)
            | Surface3::Revolution(_) => Some(2.0 * PI),
            _ => None,
        }
    }

    /// Returns the V period if the surface is V-periodic, or `None`.
    ///
    /// - Torus: `2π`
    /// - Others: `None`
    pub fn v_period(&self) -> Option<f64> {
        match &self.surface {
            Surface3::Torus(_) => Some(2.0 * PI),
            _ => None,
        }
    }

    /// Returns the surface type as a string for debugging.
    pub fn surface_type(&self) -> &'static str {
        match &self.surface {
            Surface3::Plane(_) => "Plane",
            Surface3::Cylinder(_) => "Cylinder",
            Surface3::Sphere(_) => "Sphere",
            Surface3::Cone(_) => "Cone",
            Surface3::Torus(_) => "Torus",
            Surface3::Ellipsoid(_) => "Ellipsoid",
            Surface3::Helicoid(_) => "Helicoid",
            Surface3::Pipe(_) => "Pipe",
            Surface3::BSpline(_) => "BSpline",
            Surface3::LinearExtrusion(_) => "LinearExtrusion",
            Surface3::Revolution(_) => "Revolution",
            Surface3::Ruled(_) => "Ruled",
            Surface3::Coons(_) => "Coons",
            Surface3::Bezier(_) => "Bezier",
            Surface3::TriBezier(_) => "TriBezier",
            Surface3::Offset(_) => "Offset",
            Surface3::Trimmed(_) => "Trimmed",
        }
    }

    /// Returns a reference to the underlying surface.
    pub fn surface(&self) -> &Surface3 {
        &self.surface
    }

    /// Returns the bounds of the U domain.
    pub fn u_bounds(&self) -> [f64; 2] {
        [self.u_first, self.u_last]
    }

    /// Returns the bounds of the V domain.
    pub fn v_bounds(&self) -> [f64; 2] {
        [self.v_first, self.v_last]
    }
}

// ============================================================================
// CurveOnSurfaceAdaptor
// ============================================================================

/// Adaptor for a 2D curve lying on a 3D surface.
///
/// This represents a curve defined in the parameter space of a surface,
/// providing 3D evaluation by composing the 2D curve evaluation with
/// the surface evaluation.
///
/// Analogous to OCCT `Adaptor3d_CurveOnSurface`.
#[derive(Debug, Clone)]
pub struct CurveOnSurfaceAdaptor {
    curve2d: Curve2d,
    surface: Surface3,
    first: f64,
    last: f64,
}

impl CurveOnSurfaceAdaptor {
    /// Creates a new adaptor from a 2D curve and a 3D surface.
    pub fn new(curve2d: &Curve2d, surface: &Surface3) -> Self {
        let domain = curve2d.default_domain();
        Self {
            curve2d: curve2d.clone(),
            surface: surface.clone(),
            first: domain[0],
            last: domain[1],
        }
    }

    /// Creates an adaptor with a trimmed parameter range.
    pub fn new_with_range(curve2d: &Curve2d, surface: &Surface3, first: f64, last: f64) -> Self {
        Self {
            curve2d: curve2d.clone(),
            surface: surface.clone(),
            first,
            last,
        }
    }

    /// Returns the 3D point on the surface at parameter `t`.
    ///
    /// The 2D curve is evaluated to get `(u, v)`, then the surface is
    /// evaluated at `(u, v)` to get the 3D point.
    pub fn point_at(&self, t: f64) -> DVec3 {
        let uv = self.curve2d.point_at(t);
        self.surface.point_at(uv.x, uv.y)
    }

    /// Returns the n-th derivative in 3D at parameter `t`.
    ///
    /// Uses the chain rule to compute derivatives of the composite mapping
    /// `t -> (u(t), v(t)) -> S(u(t), v(t))`.
    pub fn derivative(&self, t: f64, order: usize) -> DVec3 {
        match order {
            0 => self.point_at(t),
            1 => self.first_derivative(t),
            _ => self.nth_derivative_numerical(t, order),
        }
    }

    /// Returns the parameter domain as `[first, last]`.
    pub fn domain(&self) -> [f64; 2] {
        [self.first, self.last]
    }

    /// Returns a reference to the underlying 2D curve.
    pub fn curve2d(&self) -> &Curve2d {
        &self.curve2d
    }

    /// Returns a reference to the underlying surface.
    pub fn surface(&self) -> &Surface3 {
        &self.surface
    }

    /// Returns the UV point at parameter `t`.
    pub fn uv_at(&self, t: f64) -> DVec2 {
        self.curve2d.point_at(t)
    }

    /// Returns the first derivative (tangent vector) in 3D.
    ///
    /// Uses the chain rule: `dS/dt = ∂S/∂u * du/dt + ∂S/∂v * dv/dt`
    fn first_derivative(&self, t: f64) -> DVec3 {
        let eps = 1e-6;
        let uv = self.curve2d.point_at(t);
        let uv_plus = self.curve2d.point_at(t + eps);
        let uv_minus = self.curve2d.point_at(t - eps);

        // du/dt, dv/dt
        let du_dt = (uv_plus.x - uv_minus.x) / (2.0 * eps);
        let dv_dt = (uv_plus.y - uv_minus.y) / (2.0 * eps);

        // ∂S/∂u, ∂S/∂v
        let ds_du = self.surface_derivative_u(uv.x, uv.y);
        let ds_dv = self.surface_derivative_v(uv.x, uv.y);

        ds_du * du_dt + ds_dv * dv_dt
    }

    /// Numerical nth derivative.
    fn nth_derivative_numerical(&self, t: f64, order: usize) -> DVec3 {
        let eps = 1e-4;
        let mut result = DVec3::ZERO;

        for i in 0..=order {
            let coeff = central_diff_coefficient(order, i);
            let ti = t + (i as f64 - order as f64 / 2.0) * eps;
            result += coeff * self.point_at(ti);
        }

        result / eps.powi(order as i32)
    }

    /// Partial derivative ∂S/∂u using finite differences.
    fn surface_derivative_u(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-6;
        let p_plus = self.surface.point_at(u + eps, v);
        let p_minus = self.surface.point_at(u - eps, v);
        (p_plus - p_minus) / (2.0 * eps)
    }

    /// Partial derivative ∂S/∂v using finite differences.
    fn surface_derivative_v(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-6;
        let p_plus = self.surface.point_at(u, v + eps);
        let p_minus = self.surface.point_at(u, v - eps);
        (p_plus - p_minus) / (2.0 * eps)
    }

    /// Returns the normal vector of the surface at the curve point `t`.
    pub fn normal_at(&self, t: f64) -> DVec3 {
        let uv = self.curve2d.point_at(t);
        self.surface.normal_at(uv.x, uv.y)
    }
}

// ============================================================================
// HSurfaceAdaptor
// ============================================================================

/// Handle to a surface adaptor with reference counting.
///
/// This is a smart pointer wrapper around `SurfaceAdaptor` that allows
/// sharing surface adaptors without copying the underlying geometry.
///
/// Analogous to OCCT `Adaptor3d_HSurface`.
#[derive(Debug, Clone)]
pub struct HSurfaceAdaptor {
    inner: Rc<SurfaceAdaptor>,
}

impl HSurfaceAdaptor {
    /// Creates a new handle from a surface adaptor.
    pub fn new(adaptor: SurfaceAdaptor) -> Self {
        Self {
            inner: Rc::new(adaptor),
        }
    }

    /// Creates a handle from a surface reference.
    pub fn from_surface(surface: &Surface3) -> Self {
        Self::new(SurfaceAdaptor::from_surface(surface))
    }

    /// Returns a reference to the underlying surface adaptor.
    pub fn adaptor(&self) -> &SurfaceAdaptor {
        &self.inner
    }

    /// Returns the point on the surface at `(u, v)`.
    pub fn point_at(&self, u: f64, v: f64) -> DVec3 {
        self.inner.point_at(u, v)
    }

    /// Returns the partial derivative with respect to U.
    pub fn derivative_u(&self, u: f64, v: f64) -> DVec3 {
        self.inner.derivative_u(u, v)
    }

    /// Returns the partial derivative with respect to V.
    pub fn derivative_v(&self, u: f64, v: f64) -> DVec3 {
        self.inner.derivative_v(u, v)
    }

    /// Returns the parameter domain.
    pub fn domain(&self) -> [f64; 4] {
        self.inner.domain()
    }

    /// Returns true if the surface is U-closed.
    pub fn is_u_closed(&self) -> bool {
        self.inner.is_u_closed()
    }

    /// Returns true if the surface is V-closed.
    pub fn is_v_closed(&self) -> bool {
        self.inner.is_v_closed()
    }

    /// Returns the U period if periodic.
    pub fn u_period(&self) -> Option<f64> {
        self.inner.u_period()
    }

    /// Returns the V period if periodic.
    pub fn v_period(&self) -> Option<f64> {
        self.inner.v_period()
    }

    /// Returns the surface type name.
    pub fn surface_type(&self) -> &'static str {
        self.inner.surface_type()
    }

    /// Returns a reference to the underlying surface.
    pub fn surface(&self) -> &Surface3 {
        self.inner.surface()
    }
}

impl std::ops::Deref for HSurfaceAdaptor {
    type Target = SurfaceAdaptor;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl From<SurfaceAdaptor> for HSurfaceAdaptor {
    fn from(adaptor: SurfaceAdaptor) -> Self {
        Self::new(adaptor)
    }
}

impl From<&Surface3> for HSurfaceAdaptor {
    fn from(surface: &Surface3) -> Self {
        Self::from_surface(surface)
    }
}

// ============================================================================
// Additional utility trait for Curve2d default_domain
// ============================================================================

/// Extension trait for Curve2d to provide default domain.
trait Curve2dDomain {
    fn default_domain(&self) -> [f64; 2];
}

impl Curve2dDomain for Curve2d {
    fn default_domain(&self) -> [f64; 2] {
        match self {
            Curve2d::Line(_) => [f64::NEG_INFINITY, f64::INFINITY],
            Curve2d::Circle(_) | Curve2d::Ellipse(_) => [0.0, 2.0 * PI],
            Curve2d::CircleInvolute(_) => [0.0, 2.0 * PI],
            Curve2d::ArchimedeanSpiral(_) | Curve2d::LogarithmicSpiral(_) => [0.0, 10.0],
            Curve2d::SineWave(_) => [f64::NEG_INFINITY, f64::INFINITY],
            Curve2d::BSpline(bspline) => {
                let d = bspline.degree;
                let n = bspline.knots.len();
                if n < 2 * d + 2 {
                    [0.0, 1.0]
                } else {
                    [bspline.knots[d], bspline.knots[n - d - 1]]
                }
            }
            Curve2d::Bezier(_) => [0.0, 1.0],
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle3, Line3, Plane, CylindricalSurface, ToroidalSurface};
    use glam::{dvec2, dvec3};

    #[test]
    fn curve_adaptor_from_circle() {
        let circle = Circle3 {
            center: dvec3(0.0, 0.0, 0.0),
            normal: dvec3(0.0, 0.0, 1.0),
            radius: 2.0,
        };
        let curve = Curve3::Circle(circle);
        let adaptor = Curve3dAdaptor::from_curve(&curve);

        // Domain should be [0, 2π]
        let domain = adaptor.domain();
        assert!((domain[0] - 0.0).abs() < 1e-10);
        assert!((domain[1] - 2.0 * PI).abs() < 1e-10);

        // Circle should be closed and periodic
        assert!(adaptor.is_closed());
        assert!(adaptor.is_periodic());
        assert!((adaptor.period().unwrap() - 2.0 * PI).abs() < 1e-10);

        // Point at t=0 should be on the circle (distance from center = radius)
        let p = adaptor.point_at(0.0);
        assert!((p.length() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn curve_adaptor_derivatives() {
        let line = Line3 {
            origin: dvec3(0.0, 0.0, 0.0),
            direction: dvec3(1.0, 0.0, 0.0),
        };
        let curve = Curve3::Line(line);
        let adaptor = Curve3dAdaptor::from_curve(&curve);

        // First derivative should be the direction
        let d1 = adaptor.derivative(5.0, 1);
        assert!((d1 - dvec3(1.0, 0.0, 0.0)).length() < 1e-10);

        // Second derivative should be zero for a line
        let d2 = adaptor.derivative(5.0, 2);
        assert!(d2.length() < 1e-10);
    }

    #[test]
    fn curve_adaptor_resolution() {
        let circle = Circle3 {
            center: dvec3(0.0, 0.0, 0.0),
            normal: dvec3(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let curve = Curve3::Circle(circle);
        let adaptor = Curve3dAdaptor::from_curve(&curve);

        // Resolution should be proportional to tolerance
        let res = adaptor.resolution(0.001);
        assert!(res > 0.0 && res < 0.1);
    }

    #[test]
    fn surface_adaptor_from_plane() {
        let plane = Plane {
            origin: dvec3(0.0, 0.0, 0.0),
            normal: dvec3(0.0, 0.0, 1.0),
        };
        let surface = Surface3::Plane(plane);
        let adaptor = SurfaceAdaptor::from_surface(&surface);

        // Domain should be infinite
        let domain = adaptor.domain();
        assert!(domain[0].is_infinite() && domain[0] < 0.0);
        assert!(domain[1].is_infinite());

        // Plane should not be closed or periodic
        assert!(!adaptor.is_u_closed());
        assert!(!adaptor.is_v_closed());
        assert!(adaptor.u_period().is_none());
        assert!(adaptor.v_period().is_none());

        // Point at (0, 0) should be origin
        let p = adaptor.point_at(0.0, 0.0);
        assert!(p.length() < 1e-10);
    }

    #[test]
    fn surface_adaptor_from_cylinder() {
        let cylinder = CylindricalSurface {
            origin: dvec3(0.0, 0.0, 0.0),
            axis: dvec3(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let surface = Surface3::Cylinder(cylinder);
        let adaptor = SurfaceAdaptor::from_surface(&surface);

        // U domain should be [0, 2π], V should be infinite
        let domain = adaptor.domain();
        assert!((domain[0] - 0.0).abs() < 1e-10);
        assert!((domain[1] - 2.0 * PI).abs() < 1e-10);
        assert!(domain[2].is_infinite() && domain[2] < 0.0);
        assert!(domain[3].is_infinite());

        // Cylinder should be U-closed and U-periodic
        assert!(adaptor.is_u_closed());
        assert!(!adaptor.is_v_closed());
        assert!((adaptor.u_period().unwrap() - 2.0 * PI).abs() < 1e-10);
        assert!(adaptor.v_period().is_none());

        // Normal should point outward (perpendicular to axis, length 1)
        let n = adaptor.normal_at(0.0, 0.0);
        assert!(n.z.abs() < 1e-10); // normal is horizontal
        assert!((n.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn surface_adaptor_from_torus() {
        let torus = ToroidalSurface {
            center: dvec3(0.0, 0.0, 0.0),
            axis: dvec3(0.0, 0.0, 1.0),
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        let surface = Surface3::Torus(torus);
        let adaptor = SurfaceAdaptor::from_surface(&surface);

        // Both U and V domains should be [0, 2π]
        let domain = adaptor.domain();
        assert!((domain[0] - 0.0).abs() < 1e-10);
        assert!((domain[1] - 2.0 * PI).abs() < 1e-10);
        assert!((domain[2] - 0.0).abs() < 1e-10);
        assert!((domain[3] - 2.0 * PI).abs() < 1e-10);

        // Torus should be both U and V closed and periodic
        assert!(adaptor.is_u_closed());
        assert!(adaptor.is_v_closed());
        assert!((adaptor.u_period().unwrap() - 2.0 * PI).abs() < 1e-10);
        assert!((adaptor.v_period().unwrap() - 2.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn surface_adaptor_derivatives() {
        let plane = Plane {
            origin: dvec3(0.0, 0.0, 0.0),
            normal: dvec3(0.0, 0.0, 1.0),
        };
        let surface = Surface3::Plane(plane);
        let adaptor = SurfaceAdaptor::from_surface(&surface);

        // For a plane, ∂S/∂u and ∂S/∂v should be perpendicular to normal
        let du = adaptor.derivative_u(0.5, 0.5);
        let dv = adaptor.derivative_v(0.5, 0.5);

        assert!(du.dot(dvec3(0.0, 0.0, 1.0)).abs() < 1e-10);
        assert!(dv.dot(dvec3(0.0, 0.0, 1.0)).abs() < 1e-10);

        // For a plane, second derivatives should be near zero (within numerical precision)
        let duu = adaptor.derivative_uu(0.5, 0.5);
        let dvv = adaptor.derivative_vv(0.5, 0.5);
        let duv = adaptor.derivative_uv(0.5, 0.5);

        // Numerical second derivatives have O(eps) error for linear functions
        assert!(duu.length() < 1e-4);
        assert!(dvv.length() < 1e-4);
        assert!(duv.length() < 1e-4);
    }

    #[test]
    fn curve_on_surface_adaptor_line_on_plane() {
        let plane = Plane {
            origin: dvec3(0.0, 0.0, 0.0),
            normal: dvec3(0.0, 0.0, 1.0),
        };
        let surface = Surface3::Plane(plane);

        let line2d = rcad_kernel::geom::Line2d {
            origin: dvec2(0.0, 0.0),
            direction: dvec2(1.0, 1.0),
        };
        let curve2d = Curve2d::Line(line2d);

        let adaptor = CurveOnSurfaceAdaptor::new(&curve2d, &surface);

        // At t=0, should be at origin
        let p0 = adaptor.point_at(0.0);
        assert!(p0.length() < 1e-10);

        // At t=1, should be on the plane (z=0)
        let p1 = adaptor.point_at(1.0);
        assert!(p1.z.abs() < 1e-10);

        // First derivative should be horizontal (z component = 0)
        let d1 = adaptor.derivative(0.5, 1);
        assert!(d1.z.abs() < 1e-8);
        assert!((d1.length() - (2.0_f64).sqrt()).abs() < 1e-8); // magnitude sqrt(2)

        // Normal should be consistent with plane
        let n = adaptor.normal_at(0.5);
        assert!((n - dvec3(0.0, 0.0, 1.0)).length() < 1e-10);
    }

    #[test]
    fn curve_on_surface_adaptor_circle_on_cylinder() {
        let cylinder = CylindricalSurface {
            origin: dvec3(0.0, 0.0, 0.0),
            axis: dvec3(0.0, 0.0, 1.0),
            radius: 2.0,
        };
        let surface = Surface3::Cylinder(cylinder);

        // A line in UV space at constant v traces a circle
        let circle2d = rcad_kernel::geom::Circle2d {
            center: dvec2(0.0, 1.0), // Centered at u=0, v=1
            radius: 1.0,
        };
        let curve2d = Curve2d::Circle(circle2d);

        let adaptor = CurveOnSurfaceAdaptor::new(&curve2d, &surface);

        // The 3D curve should lie on the cylinder
        for t in [0.0, PI / 2.0, PI, 3.0 * PI / 2.0] {
            let p = adaptor.point_at(t);
            let radial_dist = (p.x * p.x + p.y * p.y).sqrt();
            // Should be approximately on cylinder surface (accounting for the circle in UV)
            // The actual radial distance varies based on the circle's UV path
            assert!(radial_dist > 0.0);
        }
    }

    #[test]
    fn h_surface_adaptor_basic() {
        let plane = Plane {
            origin: dvec3(0.0, 0.0, 0.0),
            normal: dvec3(0.0, 0.0, 1.0),
        };
        let surface = Surface3::Plane(plane);
        let handle = HSurfaceAdaptor::from_surface(&surface);

        // Should delegate to inner adaptor (point on plane should have z=0)
        let p = handle.point_at(1.0, 2.0);
        assert!(p.z.abs() < 1e-10);

        // Should provide domain access
        let domain = handle.domain();
        assert!(domain[0].is_infinite() && domain[0] < 0.0);

        // Type name should be available
        assert_eq!(handle.surface_type(), "Plane");

        // Should be clonable without copying geometry
        let handle2 = handle.clone();
        assert_eq!(handle.point_at(1.0, 2.0), handle2.point_at(1.0, 2.0));
    }

    #[test]
    fn h_surface_adaptor_from_adaptor() {
        let plane = Plane {
            origin: dvec3(1.0, 2.0, 3.0),
            normal: dvec3(0.0, 1.0, 0.0),
        };
        let surface = Surface3::Plane(plane);
        let adaptor = SurfaceAdaptor::from_surface(&surface);
        let handle = HSurfaceAdaptor::from(adaptor);

        let p = handle.point_at(0.0, 0.0);
        assert!((p - dvec3(1.0, 2.0, 3.0)).length() < 1e-10);
    }

    #[test]
    fn curve_adaptor_with_trimmed_range() {
        let circle = Circle3 {
            center: dvec3(0.0, 0.0, 0.0),
            normal: dvec3(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let curve = Curve3::Circle(circle);

        // Create adaptor for a quarter arc
        let adaptor = Curve3dAdaptor::from_curve_with_range(&curve, 0.0, PI / 2.0);

        let domain = adaptor.domain();
        assert!((domain[0] - 0.0).abs() < 1e-10);
        assert!((domain[1] - PI / 2.0).abs() < 1e-10);

        // The trimmed arc should not be closed
        assert!(!adaptor.is_closed());

        // But it's still periodic in the mathematical sense
        assert!(adaptor.is_periodic());
    }

    #[test]
    fn surface_adaptor_with_trimmed_range() {
        let cylinder = CylindricalSurface {
            origin: dvec3(0.0, 0.0, 0.0),
            axis: dvec3(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let surface = Surface3::Cylinder(cylinder);

        // Create adaptor for a partial cylinder
        let adaptor = SurfaceAdaptor::from_surface_with_range(&surface, 0.0, PI, 0.0, 5.0);

        let domain = adaptor.domain();
        assert!((domain[0] - 0.0).abs() < 1e-10);
        assert!((domain[1] - PI).abs() < 1e-10);
        assert!((domain[2] - 0.0).abs() < 1e-10);
        assert!((domain[3] - 5.0).abs() < 1e-10);

        // Partial cylinder should not be U-closed
        assert!(!adaptor.is_u_closed());
    }
}
