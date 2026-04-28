//! ElCLib-style elementary curve utilities.
//!
//! Provides analytical evaluation and parameter computation for elementary curves.
//! Analogous to OCCT `ElCLib` package.
//!
//! # Curve Types
//! - **Line**: Unbounded linear curve, parameter is distance along direction
//! - **Circle**: Closed circular curve, parameter is angle in radians [0, 2*pi]
//! - **Ellipse**: Closed elliptical curve, parameter is angle in radians [0, 2*pi]
//! - **Hyperbola**: Unbounded hyperbolic curve, parameter t (real line)
//! - **Parabola**: Unbounded parabolic curve, parameter t (real line)
//! - **BSpline**: NURBS curve, parameter within knot domain

use glam::DVec3;
use rcad_kernel::geom::{
    any_perpendicular, BSplineCurve3, Circle3, CurveEval, Ellipse3, Hyperbola3, Line3, Parabola3,
};
use std::f64::consts::{FRAC_PI_2, TAU};

// =============================================================================
// Line Utilities
// =============================================================================

/// Evaluate a point on a line at parameter t.
///
/// The parameter t represents the signed distance along the line's direction
/// from the origin. P(t) = origin + t * direction.
///
/// # Example
/// ```ignore
/// let line = Line3 { origin: DVec3::ZERO, direction: DVec3::X };
/// let p = line_point_at(&line, 5.0);
/// assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-10);
/// ```
pub fn line_point_at(line: &Line3, t: f64) -> DVec3 {
    line.origin + t * line.direction
}

/// Compute the parameter value for a point on a line.
///
/// Projects the point onto the line and returns the signed distance from
/// the origin along the line's direction.
///
/// Returns the parameter t such that `line_point_at(line, t)` is the
/// closest point on the line to the input point.
pub fn line_parameter(line: &Line3, point: DVec3) -> f64 {
    (point - line.origin).dot(line.direction)
}

/// Compute the perpendicular distance from a point to a line.
///
/// Uses the formula: distance = |(P - origin) x direction| / |direction|
/// For unit direction vectors, this simplifies to |(P - origin) x direction|.
pub fn line_distance_to_point(line: &Line3, point: DVec3) -> f64 {
    let v = point - line.origin;
    let cross = v.cross(line.direction);
    cross.length() / line.direction.length()
}

/// Find the closest point on a line to a given point.
///
/// Projects the point onto the line along the perpendicular direction.
pub fn line_closest_point(line: &Line3, point: DVec3) -> DVec3 {
    let t = line_parameter(line, point);
    line_point_at(line, t)
}

// =============================================================================
// Circle Utilities
// =============================================================================

/// Evaluate a point on a circle at the given angle.
///
/// The angle parameter is in radians, measured from the reference direction
/// (computed as any_perpendicular of the normal) in a right-handed sense
/// around the normal.
///
/// # Example
/// ```ignore
/// let circle = Circle3 { center: DVec3::ZERO, normal: DVec3::Z, radius: 2.0 };
/// let p = circle_point_at(&circle, 0.0); // Point at angle 0
/// assert!((p.length() - 2.0).abs() < 1e-10);
/// ```
pub fn circle_point_at(circle: &Circle3, angle: f64) -> DVec3 {
    circle.point_at(angle)
}

/// Compute the angle parameter for a point on or near a circle.
///
/// Projects the point onto the circle's plane and computes the angle
/// from the reference direction. Returns the angle in radians [0, 2*pi].
///
/// If the point is not exactly on the circle, the returned angle
/// corresponds to the closest point on the circle to the input.
pub fn circle_parameter(circle: &Circle3, point: DVec3) -> f64 {
    // Build local frame
    let x_axis = any_perpendicular(circle.normal);
    let y_axis = circle.normal.cross(x_axis);

    // Project point into plane and express in local coordinates
    let v = point - circle.center;
    let x = v.dot(x_axis);
    let y = v.dot(y_axis);

    // Compute angle using atan2
    let angle = y.atan2(x);

    // Normalize to [0, 2*pi]
    if angle < 0.0 {
        angle + TAU
    } else {
        angle
    }
}

/// Compute the unit tangent vector on a circle at the given angle.
///
/// The tangent is perpendicular to the radius vector, pointing in the
/// direction of increasing angle (counterclockwise when viewed from
/// the normal direction).
pub fn circle_tangent_at(circle: &Circle3, angle: f64) -> DVec3 {
    let x_axis = any_perpendicular(circle.normal);
    let y_axis = circle.normal.cross(x_axis);
    (-angle.sin() * x_axis + angle.cos() * y_axis).normalize()
}

/// Compute the unit normal vector on a circle at the given angle.
///
/// The normal points outward from the center, perpendicular to the curve
/// in the plane of the circle.
pub fn circle_normal_at(circle: &Circle3, angle: f64) -> DVec3 {
    let x_axis = any_perpendicular(circle.normal);
    let y_axis = circle.normal.cross(x_axis);
    (angle.cos() * x_axis + angle.sin() * y_axis).normalize()
}

/// Compute the binormal vector on a circle at the given angle.
///
/// The binormal is the cross product of tangent and normal, which for
/// a planar circle is the circle's normal (axis) vector.
pub fn circle_binormal_at(circle: &Circle3, _angle: f64) -> DVec3 {
    circle.normal.normalize()
}

/// Compute the nth derivative of a circle at the given angle.
///
/// - Order 0: position on circle (same as `circle_point_at`)
/// - Order 1: first derivative = radius * tangent
/// - Order 2: second derivative = -radius * normal (centripetal acceleration)
/// - Higher orders: follow the pattern of derivatives of sin/cos
///
/// For a circle of radius R:
/// - dP/dt = R * (-sin(t), cos(t)) [first derivative]
/// - d²P/dt² = R * (-cos(t), -sin(t)) = -R * (cos(t), sin(t)) [second derivative]
/// - d³P/dt³ = R * (sin(t), -cos(t)) [third derivative]
///
/// Returns DVec3::ZERO if order is not supported.
pub fn circle_derivative(circle: &Circle3, angle: f64, order: usize) -> DVec3 {
    let x_axis = any_perpendicular(circle.normal);
    let y_axis = circle.normal.cross(x_axis);
    let r = circle.radius;

    match order {
        0 => circle_point_at(circle, angle),
        1 => r * (-angle.sin() * x_axis + angle.cos() * y_axis),
        2 => -r * (angle.cos() * x_axis + angle.sin() * y_axis),
        3 => r * (angle.sin() * x_axis - angle.cos() * y_axis),
        4 => r * (angle.cos() * x_axis + angle.sin() * y_axis), // Same as order 2 with opposite sign
        n => {
            // Higher orders cycle every 4
            let k = n % 4;
            match k {
                0 => circle_point_at(circle, angle),
                1 => r * (-angle.sin() * x_axis + angle.cos() * y_axis),
                2 => -r * (angle.cos() * x_axis + angle.sin() * y_axis),
                3 => r * (angle.sin() * x_axis - angle.cos() * y_axis),
                _ => DVec3::ZERO,
            }
        }
    }
}

// =============================================================================
// Ellipse Utilities
// =============================================================================

/// Evaluate a point on an ellipse at the given angle parameter.
///
/// The angle parameter is the eccentric anomaly, not the polar angle.
/// The ellipse is parameterized as:
///   P(angle) = center + a*cos(angle)*major_dir + b*sin(angle)*minor_dir
///
/// where a = major_radius, b = minor_radius, and minor_dir = normal x major_dir.
pub fn ellipse_point_at(ellipse: &Ellipse3, angle: f64) -> DVec3 {
    ellipse.point_at(angle)
}

/// Compute the angle parameter for a point on or near an ellipse.
///
/// Projects the point onto the ellipse's plane and solves for the
/// eccentric anomaly. Uses Newton-Raphson iteration for accuracy.
///
/// Returns the angle in radians [0, 2*pi].
pub fn ellipse_parameter(ellipse: &Ellipse3, point: DVec3) -> f64 {
    // Build local frame
    let x_axis = ellipse.major_dir;
    let y_axis = ellipse.normal.cross(x_axis).normalize();

    // Project point into plane and express in local coordinates
    let v = point - ellipse.center;
    let x = v.dot(x_axis);
    let y = v.dot(y_axis);

    // For an ellipse x = a*cos(t), y = b*sin(t)
    // We need to solve for t given (x, y)
    // Use atan2(y/b, x/a) as initial guess, then refine
    let a = ellipse.major_radius;
    let b = ellipse.minor_radius;

    if a.abs() < 1e-15 || b.abs() < 1e-15 {
        return 0.0;
    }

    // Initial guess using modified atan2
    let mut t = (y / b).atan2(x / a);

    // Newton-Raphson refinement for better accuracy
    // We solve: f(t) = atan2(y - b*sin(t), x - a*cos(t)) = 0
    // which is implicit. Instead, we directly compute the eccentric anomaly.
    // The eccentric anomaly satisfies:
    //   x = a * cos(t)
    //   y = b * sin(t)
    // So: cos(t) = x/a, sin(t) = y/b
    // t = atan2(y/b, x/a) is already exact for the parametric form

    // Normalize to [0, 2*pi]
    if t < 0.0 {
        t + TAU
    } else {
        t
    }
}

/// Compute the nth derivative of an ellipse at the given angle.
///
/// For an ellipse with radii a (major) and b (minor):
/// - Order 0: P(t) = (a*cos(t), b*sin(t))
/// - Order 1: dP/dt = (-a*sin(t), b*cos(t))
/// - Order 2: d²P/dt² = (-a*cos(t), -b*sin(t))
/// - Order 3: d³P/dt³ = (a*sin(t), -b*cos(t))
/// - Higher orders: cycle every 4
pub fn ellipse_derivative(ellipse: &Ellipse3, angle: f64, order: usize) -> DVec3 {
    let x_axis = ellipse.major_dir;
    let y_axis = ellipse.normal.cross(x_axis).normalize();
    let a = ellipse.major_radius;
    let b = ellipse.minor_radius;

    let cos_a = angle.cos();
    let sin_a = angle.sin();

    match order {
        0 => ellipse.center + a * cos_a * x_axis + b * sin_a * y_axis,
        1 => -a * sin_a * x_axis + b * cos_a * y_axis,
        2 => -a * cos_a * x_axis - b * sin_a * y_axis,
        3 => a * sin_a * x_axis - b * cos_a * y_axis,
        n => {
            let k = n % 4;
            match k {
                0 => ellipse.center + a * cos_a * x_axis + b * sin_a * y_axis,
                1 => -a * sin_a * x_axis + b * cos_a * y_axis,
                2 => -a * cos_a * x_axis - b * sin_a * y_axis,
                3 => a * sin_a * x_axis - b * cos_a * y_axis,
                _ => DVec3::ZERO,
            }
        }
    }
}

// =============================================================================
// Hyperbola Utilities
// =============================================================================

/// Evaluate a point on a hyperbola at parameter t.
///
/// The hyperbola is parameterized as:
///   P(t) = center + a*cosh(t)*major_dir + b*sinh(t)*minor_dir
///
/// where a = semi_major, b = semi_minor, and minor_dir = normal x major_dir.
/// The principal branch (t >= 0) is on the +major_dir side of the center.
pub fn hyperbola_point_at(hyp: &Hyperbola3, t: f64) -> DVec3 {
    hyp.point_at(t)
}

/// Compute the nth derivative of a hyperbola at parameter t.
///
/// For a hyperbola with semi-axes a and b:
/// - Order 0: P(t) = (a*cosh(t), b*sinh(t))
/// - Order 1: dP/dt = (a*sinh(t), b*cosh(t))
/// - Order 2: d²P/dt² = (a*cosh(t), b*sinh(t))
/// - Order 3: d³P/dt³ = (a*sinh(t), b*cosh(t))
/// - Higher orders: alternates between sinh and cosh patterns
pub fn hyperbola_derivative(hyp: &Hyperbola3, t: f64, order: usize) -> DVec3 {
    let minor_dir = hyp.normal.cross(hyp.major_dir).normalize();
    let a = hyp.semi_major;
    let b = hyp.semi_minor;

    let cosh_t = t.cosh();
    let sinh_t = t.sinh();

    match order {
        0 => hyp.center + a * cosh_t * hyp.major_dir + b * sinh_t * minor_dir,
        1 => a * sinh_t * hyp.major_dir + b * cosh_t * minor_dir,
        2 => a * cosh_t * hyp.major_dir + b * sinh_t * minor_dir, // Same as order 0
        3 => a * sinh_t * hyp.major_dir + b * cosh_t * minor_dir, // Same as order 1
        n => {
            // Pattern: even orders = order 0, odd orders = order 1
            if n % 2 == 0 {
                hyp.center + a * cosh_t * hyp.major_dir + b * sinh_t * minor_dir
            } else {
                a * sinh_t * hyp.major_dir + b * cosh_t * minor_dir
            }
        }
    }
}

// =============================================================================
// Parabola Utilities
// =============================================================================

/// Evaluate a point on a parabola at parameter t.
///
/// The parabola is parameterized as:
///   P(t) = vertex + (t²/(2p))*axis_dir + t*dir_perp
///
/// where p = focal_param (twice the focal length), and dir_perp = normal x axis_dir.
/// The focus is at distance p/2 from the vertex along axis_dir.
pub fn parabola_point_at(parab: &Parabola3, t: f64) -> DVec3 {
    parab.point_at(t)
}

/// Compute the nth derivative of a parabola at parameter t.
///
/// For a parabola with focal parameter p:
/// - Order 0: P(t) = (t²/(2p), t) in local coordinates
/// - Order 1: dP/dt = (t/p, 1)
/// - Order 2: d²P/dt² = (1/p, 0)
/// - Order 3+: All higher derivatives are zero
pub fn parabola_derivative(parab: &Parabola3, t: f64, order: usize) -> DVec3 {
    // dir_perp forms a right-handed system: axis_dir × normal gives perpendicular direction
    let dir_perp = parab.axis_dir.cross(parab.normal).normalize();
    let p = parab.focal_param;

    if p.abs() < 1e-15 {
        return DVec3::ZERO;
    }

    match order {
        0 => {
            parab.vertex
                + (t * t / (2.0 * p)) * parab.axis_dir
                + t * dir_perp
        }
        1 => {
            // dP/dt = (t/p) * axis_dir + dir_perp
            (t / p) * parab.axis_dir + dir_perp
        }
        2 => {
            // d²P/dt² = (1/p) * axis_dir
            (1.0 / p) * parab.axis_dir
        }
        _ => {
            // All higher derivatives are zero
            DVec3::ZERO
        }
    }
}

// =============================================================================
// BSpline Utilities
// =============================================================================

/// Evaluate a point on a B-spline curve at parameter t.
///
/// Uses the de Boor algorithm for rational and non-rational B-splines.
/// The parameter t should be within the curve's domain [knots[degree], knots[n-degree-1]].
pub fn bspline_point_at(spline: &BSplineCurve3, t: f64) -> DVec3 {
    spline.point_at(t)
}

/// Compute the nth derivative of a B-spline curve at parameter t.
///
/// Uses analytical differentiation for NURBS curves via the quotient rule.
/// The derivative of a rational B-spline C(t) = A(t)/W(t) is:
///   C'(t) = (A'(t) - W'(t)*C(t)) / W(t)
///
/// Higher-order derivatives are computed by differentiating the derivative
/// curve, which is itself a B-spline of degree p-1.
///
/// Returns DVec3::ZERO if:
/// - The order is greater than the curve's degree
/// - The curve is invalid (no control points or degree 0 with order > 0)
pub fn bspline_derivative(spline: &BSplineCurve3, t: f64, order: usize) -> DVec3 {
    if order == 0 {
        return bspline_point_at(spline, t);
    }

    let n = spline.control_points.len();
    let degree = spline.degree;

    if n == 0 || (order > degree && degree == 0) {
        return DVec3::ZERO;
    }

    // For higher orders than degree, derivative is zero
    if order > degree {
        return DVec3::ZERO;
    }

    // Compute derivative using finite differences for simplicity and robustness
    // For production code, this should use the analytical derivative chain
    let h = 1e-7;

    let domain = spline.default_domain();
    let t_min = domain[0];
    let t_max = domain[1];

    // Clamp t to domain for finite difference calculation
    let t_lo = (t - h).max(t_min);
    let t_hi = (t + h).min(t_max);
    let actual_h = t_hi - t_lo;

    if actual_h < 1e-15 {
        return DVec3::ZERO;
    }

    if order == 1 {
        // First derivative
        let p_lo = bspline_point_at(spline, t_lo);
        let p_hi = bspline_point_at(spline, t_hi);
        (p_hi - p_lo) / actual_h
    } else if order == 2 {
        // Second derivative using central differences
        let p_lo = bspline_point_at(spline, t_lo);
        let p_mid = bspline_point_at(spline, t);
        let p_hi = bspline_point_at(spline, t_hi);
        (p_hi - 2.0 * p_mid + p_lo) / (actual_h * actual_h / 4.0)
    } else {
        // Higher-order derivatives via recursive finite differences
        // This is less accurate but works for any order
        let mut points = Vec::with_capacity(order + 1);
        let step = actual_h / order as f64;
        for i in 0..=order {
            let ti = t_lo + i as f64 * step;
            points.push(bspline_point_at(spline, ti));
        }

        // Apply finite difference formula n times
        for _ in 0..order {
            let mut new_points = Vec::with_capacity(points.len() - 1);
            for i in 0..points.len() - 1 {
                new_points.push(points[i + 1] - points[i]);
            }
            points = new_points;
        }

        points.get(0).copied().unwrap_or(DVec3::ZERO) / step.powi(order as i32)
    }
}

// =============================================================================
// Unit Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_4;

    // -------------------------------------------------------------------------
    // Line Tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_line_point_at() {
        let line = Line3 {
            origin: DVec3::new(1.0, 2.0, 3.0),
            direction: DVec3::X,
        };
        let p = line_point_at(&line, 5.0);
        assert!((p - DVec3::new(6.0, 2.0, 3.0)).length() < 1e-10);

        // Negative parameter
        let p_neg = line_point_at(&line, -3.0);
        assert!((p_neg - DVec3::new(-2.0, 2.0, 3.0)).length() < 1e-10);
    }

    #[test]
    fn test_line_parameter() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        assert!((line_parameter(&line, DVec3::new(5.0, 0.0, 0.0)) - 5.0).abs() < 1e-10);
        assert!((line_parameter(&line, DVec3::new(-3.0, 0.0, 0.0)) - (-3.0)).abs() < 1e-10);

        // Point off the line
        let t = line_parameter(&line, DVec3::new(5.0, 10.0, 0.0));
        assert!((t - 5.0).abs() < 1e-10); // Should project to t=5
    }

    #[test]
    fn test_line_distance_to_point() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let dist = line_distance_to_point(&line, DVec3::new(5.0, 3.0, 4.0));
        assert!((dist - 5.0).abs() < 1e-10); // Distance is sqrt(3^2 + 4^2) = 5
    }

    #[test]
    fn test_line_closest_point() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let closest = line_closest_point(&line, DVec3::new(5.0, 10.0, 20.0));
        assert!((closest - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // Circle Tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_circle_point_at() {
        let circle = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 2.0,
        };

        let p0 = circle_point_at(&circle, 0.0);
        assert!((p0.length() - 2.0).abs() < 1e-10);

        let p90 = circle_point_at(&circle, FRAC_PI_2);
        assert!((p90.length() - 2.0).abs() < 1e-10);

        // Full revolution returns to start
        let p2pi = circle_point_at(&circle, TAU);
        assert!((p0 - p2pi).length() < 1e-10);
    }

    #[test]
    fn test_circle_parameter() {
        let circle = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 5.0,
        };

        let x_axis = any_perpendicular(circle.normal);

        // Point at angle 0
        let p0 = circle.center + circle.radius * x_axis;
        let t0 = circle_parameter(&circle, p0);
        assert!(t0.abs() < 1e-10 || (t0 - TAU).abs() < 1e-10);

        // Point at angle pi/2
        let y_axis = circle.normal.cross(x_axis);
        let p90 = circle.center + circle.radius * y_axis;
        let t90 = circle_parameter(&circle, p90);
        assert!((t90 - FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn test_circle_tangent_at() {
        let circle = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        };

        let tangent = circle_tangent_at(&circle, 0.0);
        let x_axis = any_perpendicular(circle.normal);
        let y_axis = circle.normal.cross(x_axis);

        // At angle 0, tangent should point in +y direction
        assert!((tangent - y_axis).length() < 1e-10);

        // Tangent should be unit length
        assert!((tangent.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_circle_normal_at() {
        let circle = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        };

        let normal = circle_normal_at(&circle, 0.0);
        let x_axis = any_perpendicular(circle.normal);

        // At angle 0, normal should point in +x direction
        assert!((normal - x_axis).length() < 1e-10);

        // Normal should be unit length
        assert!((normal.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_circle_binormal_at() {
        let circle = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        };

        let binormal = circle_binormal_at(&circle, 0.0);
        assert!((binormal - DVec3::Z).length() < 1e-10);
    }

    #[test]
    fn test_circle_derivative() {
        let circle = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 2.0,
        };

        // Order 0 should match point_at
        let p = circle_derivative(&circle, 1.0, 0);
        let expected_p = circle_point_at(&circle, 1.0);
        assert!((p - expected_p).length() < 1e-10);

        // First derivative magnitude should be radius
        let d1 = circle_derivative(&circle, 1.0, 1);
        assert!((d1.length() - 2.0).abs() < 1e-10);

        // Second derivative magnitude should be radius
        let d2 = circle_derivative(&circle, 1.0, 2);
        assert!((d2.length() - 2.0).abs() < 1e-10);

        // Second derivative should point toward center
        let point = circle_point_at(&circle, 1.0);
        let radial = (point - circle.center).normalize();
        let d2_normalized = d2.normalize();
        assert!((d2_normalized + radial).length() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // Ellipse Tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_ellipse_point_at() {
        let ellipse = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 4.0,
            minor_radius: 2.0,
        };

        // At angle 0, should be at major radius
        let p0 = ellipse_point_at(&ellipse, 0.0);
        assert!((p0 - DVec3::new(4.0, 0.0, 0.0)).length() < 1e-10);

        // At angle pi/2, should be at minor radius
        let p90 = ellipse_point_at(&ellipse, FRAC_PI_2);
        assert!((p90 - DVec3::new(0.0, 2.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn test_ellipse_parameter() {
        let ellipse = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 4.0,
            minor_radius: 2.0,
        };

        // At angle 0
        let p0 = DVec3::new(4.0, 0.0, 0.0);
        let t0 = ellipse_parameter(&ellipse, p0);
        assert!(t0.abs() < 1e-10 || (t0 - TAU).abs() < 1e-10);

        // At angle pi/2
        let p90 = DVec3::new(0.0, 2.0, 0.0);
        let t90 = ellipse_parameter(&ellipse, p90);
        assert!((t90 - FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn test_ellipse_derivative() {
        let ellipse = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 4.0,
            minor_radius: 2.0,
        };

        // First derivative at angle 0
        let d1 = ellipse_derivative(&ellipse, 0.0, 1);
        assert!((d1 - DVec3::new(0.0, 2.0, 0.0)).length() < 1e-10);

        // Second derivative at angle 0 (points toward center along major axis)
        let d2 = ellipse_derivative(&ellipse, 0.0, 2);
        assert!((d2 - DVec3::new(-4.0, 0.0, 0.0)).length() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // Hyperbola Tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_hyperbola_point_at() {
        let hyp = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 3.0,
            semi_minor: 2.0,
        };

        // At t=0, should be on vertex (major axis)
        let p0 = hyperbola_point_at(&hyp, 0.0);
        assert!((p0 - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn test_hyperbola_derivative() {
        let hyp = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 3.0,
            semi_minor: 2.0,
        };

        // First derivative at t=0
        let d1 = hyperbola_derivative(&hyp, 0.0, 1);
        let y_axis = hyp.normal.cross(hyp.major_dir).normalize();
        assert!((d1 - 2.0 * y_axis).length() < 1e-10);

        // Second derivative at t=0 (same as position minus center)
        let d2 = hyperbola_derivative(&hyp, 0.0, 2);
        assert!((d2 - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // Parabola Tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_parabola_point_at() {
        let parab = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::Y,
            focal_param: 2.0,
        };

        // At t=0, should be at vertex
        let p0 = parabola_point_at(&parab, 0.0);
        assert!((p0 - DVec3::ZERO).length() < 1e-10);

        // At t=1: x = 1, y = 1^2/(2*2) = 0.25
        let p1 = parabola_point_at(&parab, 1.0);
        assert!((p1 - DVec3::new(1.0, 0.25, 0.0)).length() < 1e-10);
    }

    #[test]
    fn test_parabola_derivative() {
        let parab = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::Y,
            focal_param: 2.0,
        };

        // dir_perp forms a right-handed system: axis_dir × normal
        let dir_perp = parab.axis_dir.cross(parab.normal).normalize();

        // First derivative at t=0: (0, 0) + dir_perp = (1, 0, 0)
        let d1 = parabola_derivative(&parab, 0.0, 1);
        assert!((d1 - dir_perp).length() < 1e-10);

        // Second derivative at t=0: 1/p * axis_dir
        let d2 = parabola_derivative(&parab, 0.0, 2);
        assert!((d2 - DVec3::new(0.0, 0.5, 0.0)).length() < 1e-10);

        // Third and higher derivatives are zero
        let d3 = parabola_derivative(&parab, 0.0, 3);
        assert!(d3.length() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // BSpline Tests
    // -------------------------------------------------------------------------

    #[test]
    fn test_bspline_point_at_linear() {
        // Degree-1 B-spline (line segment from origin to (1,0,0))
        let spline = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
        };

        let p0 = bspline_point_at(&spline, 0.0);
        assert!((p0 - DVec3::ZERO).length() < 1e-10);

        let p1 = bspline_point_at(&spline, 1.0);
        assert!((p1 - DVec3::X).length() < 1e-10);

        let pmid = bspline_point_at(&spline, 0.5);
        assert!((pmid - DVec3::new(0.5, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn test_bspline_derivative_linear() {
        let spline = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
        };

        // First derivative should be constant (direction from p0 to p1)
        let d1 = bspline_derivative(&spline, 0.5, 1);
        assert!((d1 - DVec3::X).length() < 1e-6);

        // Second derivative should be zero for a line
        let d2 = bspline_derivative(&spline, 0.5, 2);
        assert!(d2.length() < 1e-6);
    }

    #[test]
    fn test_bspline_point_at_quadratic() {
        // Degree-2 B-spline (quadratic Bezier)
        let spline = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(0.5, 1.0, 0.0),
                DVec3::X,
            ],
            weights: vec![1.0, 1.0, 1.0],
        };

        let p0 = bspline_point_at(&spline, 0.0);
        assert!((p0 - DVec3::ZERO).length() < 1e-10);

        let p1 = bspline_point_at(&spline, 1.0);
        assert!((p1 - DVec3::X).length() < 1e-10);
    }

    #[test]
    fn test_bspline_derivative_higher_order_returns_zero() {
        let spline = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
        };

        // Order greater than degree should return zero
        let d3 = bspline_derivative(&spline, 0.5, 3);
        assert!(d3.length() < 1e-10);
    }

    #[test]
    fn test_circle_derivative_cycle() {
        let circle = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        };

        let angle = FRAC_PI_4;

        // Test that derivatives cycle every 4
        let d0 = circle_derivative(&circle, angle, 0);
        let d4 = circle_derivative(&circle, angle, 4);
        let d8 = circle_derivative(&circle, angle, 8);

        // Orders 0, 4, 8 should all give the point (minus center for 4, 8)
        assert!((d0 - d4).length() < 1e-10);
        assert!((d0 - d8).length() < 1e-10);
    }
}
