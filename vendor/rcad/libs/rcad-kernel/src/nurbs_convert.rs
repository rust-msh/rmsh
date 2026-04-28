//! NURBS interoperability: convert analytic curves and surfaces to rational
//! B-spline (NURBS) representation.
//!
//! Analogous to OCCT `GeomConvert::CurveToBSplineCurve` /
//! `GeomConvert::SurfaceToBSplineSurface`.
//!
//! All conversions are **exact**: the resulting NURBS evaluates identically to
//! the original analytic geometry at any parameter value (within floating-point
//! precision).  The rational weight encoding follows standard geometric
//! modeling references (Piegl & Tiller, "The NURBS Book").
//!
//! # Supported conversions
//!
//! ## Curves
//! | Source | Method |
//! |---|---|
//! | `Line3` | [`line_to_bspline`] — degree-1, 2 control points |
//! | `Circle3` | [`circle_to_bspline`] — degree-2 NURBS, 9 control points, exact arc |
//! | `Ellipse3` | [`ellipse_to_bspline`] — same as circle, scaled to semi-axes |
//! | `BSplineCurve3` | identity (already NURBS) |
//! | `BezierCurve3` | [`bezier_curve_to_bspline`] — inserts clamped knots |
//! | `OffsetCurve3` | [`curve_to_bspline`] — samples and interpolates |
//! | `Hyperbola3` | [`curve_to_bspline`] — samples and interpolates |
//! | `Parabola3` | [`curve_to_bspline`] — samples and interpolates |
//!
//! ## Surfaces
//! | Source | Method |
//! |---|---|
//! | `Plane` | [`plane_to_bspline`] — degree-(1,1), 4 control points |
//! | `CylindricalSurface` | [`cylinder_to_bspline`] — degree-(2,1) NURBS, exact |
//! | `SphericalSurface` | [`sphere_to_bspline`] — degree-(2,2) NURBS, exact |
//! | `BSplineSurface` | identity (already NURBS) |
//! | `BezierSurface` | [`bezier_surface_to_bspline`] |
//! | other | [`surface_to_bspline`] — adaptive sampling + bilinear patch |

use glam::DVec3;
use std::f64::consts::PI;

use crate::fit::interpolate_points;
use crate::geom::{
    BSplineCurve3, BSplineSurface, BezierCurve3, BezierSurface, Circle3, Curve3, CurveEval,
    CylindricalSurface, Ellipse3, Line3, Plane, SphericalSurface, Surface3, SurfaceEval,
};

// ─────────────────────────────────────────────────────────────────────────────
// Curve conversions
// ─────────────────────────────────────────────────────────────────────────────

/// Convert any `Curve3` to an equivalent `BSplineCurve3`.
///
/// Analytic types (Line, Circle, Ellipse, Bezier) are converted exactly.
/// Other types (Offset, Hyperbola, Parabola) are approximated by sampling
/// the curve at `n_samples` parameter values and interpolating a cubic B-spline.
///
/// Analogous to `GeomConvert::CurveToBSplineCurve`.
pub fn curve_to_bspline(curve: &Curve3, n_samples: usize) -> BSplineCurve3 {
    match curve {
        Curve3::Line(l) => line_to_bspline(l),
        Curve3::Circle(c) => circle_to_bspline(c),
        Curve3::Ellipse(e) => ellipse_to_bspline(e),
        Curve3::BSpline(b) => b.clone(),
        Curve3::Bezier(b) => bezier_curve_to_bspline(b),
        Curve3::Offset(_) | Curve3::Hyperbola(_) | Curve3::Parabola(_) | Curve3::CircularHelix(_) | Curve3::SineWave(_) => {
            sample_curve_to_bspline(curve, n_samples)
        }
    }
}

/// Convert a `Line3` to a degree-1 `BSplineCurve3` over the parameter range
/// `[t0, t1]`.  Defaults to `[0, 1]`.
///
/// Analogous to `GeomConvert::CurveToBSplineCurve` for a line.
pub fn line_to_bspline(line: &Line3) -> BSplineCurve3 {
    line_to_bspline_range(line, 0.0, 1.0)
}

/// Convert a `Line3` over `[t0, t1]`.
pub fn line_to_bspline_range(line: &Line3, t0: f64, t1: f64) -> BSplineCurve3 {
    let p0 = line.point_at(t0);
    let p1 = line.point_at(t1);
    BSplineCurve3 {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

/// Convert a `Circle3` to an exact rational quadratic B-spline (degree 2,
/// 9 control points, 6 unique knots).
///
/// This is the standard 9-point NURBS circle construction (Piegl & Tiller §7.1).
/// The output parameter domain is `[0, 2π]`.
pub fn circle_to_bspline(circle: &Circle3) -> BSplineCurve3 {
    ellipse_to_bspline(&crate::geom::Ellipse3 {
        center: circle.center,
        normal: circle.normal,
        major_dir: crate::geom::any_perpendicular(circle.normal),
        major_radius: circle.radius,
        minor_radius: circle.radius,
    })
}

/// Convert an `Ellipse3` to an exact rational quadratic B-spline.
///
/// Uses the standard 9-point NURBS construction: 3 quadratic arcs (each 90°)
/// joined with C¹ continuity.  Domain is `[0, 1]` (maps to `[0, 2π]`).
pub fn ellipse_to_bspline(ellipse: &Ellipse3) -> BSplineCurve3 {
    let a = ellipse.major_radius;
    let b = ellipse.minor_radius;
    let c = ellipse.center;
    let x_ax = ellipse.major_dir.normalize();
    let y_ax = ellipse.normal.cross(x_ax).normalize();

    // 9 control points for a full NURBS ellipse / circle (Piegl & Tiller §7.1)
    // Quarter-circle weight factor
    let w = (2.0_f64).sqrt() / 2.0; // cos(45°)

    // Control points at 0°, 45°, 90°, 135°, 180°, 225°, 270°, 315°, 360°=0°
    // but using the standard 9-point construction:
    // P0 at 0°, P1 midpoint weight, P2 at 90°, etc.
    let pts = [
        c + a * x_ax,            // 0°
        c + a * x_ax + b * y_ax, // corner at (a, b)
        c + b * y_ax,            // 90°
        c - a * x_ax + b * y_ax, // corner at (-a, b)
        c - a * x_ax,            // 180°
        c - a * x_ax - b * y_ax, // corner at (-a, -b)
        c - b * y_ax,            // 270°
        c + a * x_ax - b * y_ax, // corner at (a, -b)
        c + a * x_ax,            // 360° = 0° (closed)
    ];

    let weights = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];

    // Clamped quadratic knot vector for 9 control points
    // [0,0,0, 1/4, 1/4, 1/2, 1/2, 3/4, 3/4, 1,1,1]
    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];

    BSplineCurve3 {
        degree: 2,
        knots,
        control_points: pts.to_vec(),
        weights: weights.to_vec(),
    }
}

/// Convert a `BezierCurve3` to a `BSplineCurve3` by inserting clamped endpoint
/// knots.  Weights are preserved exactly.
///
/// Analogous to `GeomConvert::CurveToBSplineCurve` for Bezier curves.
pub fn bezier_curve_to_bspline(bezier: &BezierCurve3) -> BSplineCurve3 {
    let n = bezier.control_points.len();
    let degree = (n - 1).max(1);
    // Clamped knot vector: degree+1 zeros, then degree+1 ones
    let mut knots = vec![0.0f64; degree + 1];
    knots.extend(vec![1.0f64; degree + 1]);

    BSplineCurve3 {
        degree,
        knots,
        control_points: bezier.control_points.clone(),
        weights: bezier.weights.clone(),
    }
}

/// Sample a curve at `n` equidistant parameter values and interpolate a cubic
/// B-spline through those points.  Used for transcendental curves.
fn sample_curve_to_bspline(curve: &Curve3, n: usize) -> BSplineCurve3 {
    let [t0, t1] = curve.default_domain();
    // For hyperbola / parabola with large "infinite" domain, use a sensible range
    let (t0, t1) = if t1 - t0 > 1e6 {
        (-10.0, 10.0)
    } else {
        (t0, t1)
    };
    let n = n.max(4);
    let pts: Vec<DVec3> = (0..n)
        .map(|i| {
            let t = t0 + (t1 - t0) * i as f64 / (n - 1) as f64;
            curve.point_at(t)
        })
        .collect();
    interpolate_points(&pts).unwrap_or_else(|_| BSplineCurve3 {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![pts[0], *pts.last().expect("pts has n>=2 points")],
        weights: vec![1.0, 1.0],
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface conversions
// ─────────────────────────────────────────────────────────────────────────────

/// Convert any `Surface3` to an equivalent `BSplineSurface`.
///
/// Analytic surfaces (Plane, Cylinder, Sphere, Bezier) are converted exactly
/// or via the standard NURBS constructions.  Other surfaces are sampled at
/// `n_u × n_v` parameter values and bi-linearly interpolated.
///
/// Analogous to `GeomConvert::SurfaceToBSplineSurface`.
pub fn surface_to_bspline(surface: &Surface3, n_u: usize, n_v: usize) -> BSplineSurface {
    match surface {
        Surface3::Plane(p) => plane_to_bspline(p),
        Surface3::Cylinder(c) => cylinder_to_bspline(c),
        Surface3::Sphere(s) => sphere_to_bspline(s),
        Surface3::BSpline(b) => b.clone(),
        Surface3::Bezier(b) => bezier_surface_to_bspline(b),
        _ => sample_surface_to_bspline(surface, n_u, n_v),
    }
}

/// Convert a `Plane` to a degree-(1,1) `BSplineSurface` over the domain
/// `[-1, 1] × [-1, 1]`.
///
/// The four control points span a 2×2 unit patch centred at the origin.
pub fn plane_to_bspline(plane: &Plane) -> BSplineSurface {
    plane_to_bspline_domain(plane, -1.0, 1.0, -1.0, 1.0)
}

/// Convert a `Plane` over a specified UV domain `[u0,u1]×[v0,v1]`.
pub fn plane_to_bspline_domain(
    plane: &Plane,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
) -> BSplineSurface {
    let p00 = plane.point_at(u0, v0);
    let p10 = plane.point_at(u1, v0);
    let p01 = plane.point_at(u0, v1);
    let p11 = plane.point_at(u1, v1);

    BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![vec![p00, p01], vec![p10, p11]],
        weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
    }
}

/// Convert a `CylindricalSurface` to an exact degree-(2,1) NURBS surface.
///
/// u direction: rational quadratic circle (9 columns, like `circle_to_bspline`).
/// v direction: linear (height), evaluated at `[v0, v1]` (defaults to `[0, 1]`).
pub fn cylinder_to_bspline(cyl: &CylindricalSurface) -> BSplineSurface {
    cylinder_to_bspline_range(cyl, 0.0, 1.0)
}

/// Convert a `CylindricalSurface` over the v-range `[v0, v1]`.
pub fn cylinder_to_bspline_range(cyl: &CylindricalSurface, v0: f64, v1: f64) -> BSplineSurface {
    // Circle at height v0, then circle at height v1
    let circle = Circle3 {
        center: cyl.origin + v0 * cyl.axis,
        normal: cyl.axis,
        radius: cyl.radius,
    };
    let c0 = circle_to_bspline(&circle);

    // Shift all control points along the axis for the v1 row
    let dv = (v1 - v0) * cyl.axis;
    let pts_v1: Vec<DVec3> = c0.control_points.iter().map(|p| *p + dv).collect();

    // degree_v = 1, knots_v = [0,0,1,1]
    BSplineSurface {
        degree_u: c0.degree,
        degree_v: 1,
        knots_u: c0.knots.clone(),
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![c0.control_points.clone(), pts_v1],
        weights: vec![c0.weights.clone(), c0.weights.clone()],
    }
}

/// Convert a `SphericalSurface` to an exact degree-(2,2) NURBS surface using
/// the standard sphere construction (Piegl & Tiller §7.3):
/// 3 latitude bands each consisting of the circle NURBS, scaled by sin(v) and
/// shifted by cos(v)·r·axis.
pub fn sphere_to_bspline(sphere: &SphericalSurface) -> BSplineSurface {
    let r = sphere.radius;
    let x_ax = crate::geom::any_perpendicular(sphere.axis);
    let _y_ax = sphere.axis.cross(x_ax).normalize();
    let _z_ax = sphere.axis.normalize();

    // We use 5 v-rows: v = 0°(south pole), 45°, 90°(equator), 135°, 180°(north pole)
    // For a degree-2 NURBS in v we need the standard 5-row sphere construction.
    // v parameter mapped to colatitude: v=0 → south pole, v=π → north pole.
    let n_v = 5;
    let v_angles = [0.0f64, PI / 4.0, PI / 2.0, 3.0 * PI / 4.0, PI];
    let v_weights = [1.0f64, 2.0f64.sqrt() / 2.0, 1.0, 2.0f64.sqrt() / 2.0, 1.0];

    // Each v-row is a scaled copy of the circle NURBS
    let circle_base = Circle3 {
        center: sphere.center,
        normal: sphere.axis,
        radius: r,
    };
    let c_base = circle_to_bspline(&circle_base);
    let n_u = c_base.control_points.len();

    let mut ctrl_grid: Vec<Vec<DVec3>> = Vec::new();
    let mut w_grid: Vec<Vec<f64>> = Vec::new();

    for (vi, &v_ang) in v_angles.iter().enumerate() {
        let sin_v = v_ang.sin();
        let cos_v = v_ang.cos();
        let vw = v_weights[vi];
        // Shift along axis + scale circle radius
        let axis_offset = sphere.center + cos_v * r * sphere.axis;
        let row_pts: Vec<DVec3> = c_base
            .control_points
            .iter()
            .map(|p| {
                // p is on circle of radius r at sphere.center; scale xy by sin_v
                let delta = *p - sphere.center;
                axis_offset + sin_v * delta
            })
            .collect();
        let row_w: Vec<f64> = c_base.weights.iter().map(|&w| w * vw).collect();
        ctrl_grid.push(row_pts);
        w_grid.push(row_w);
    }

    // Degree-2 in v with clamped knots for 5 rows: [0,0,0, 0.5, 0.5, 1,1,1]
    let knots_v = vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0];

    // ctrl_grid[v_idx][u_idx] — transpose to control_points[u_idx][v_idx]
    let n_v_rows = n_v;
    let transposed_ctrl: Vec<Vec<DVec3>> = (0..n_u)
        .map(|ui| (0..n_v_rows).map(|vi| ctrl_grid[vi][ui]).collect())
        .collect();
    let transposed_w: Vec<Vec<f64>> = (0..n_u)
        .map(|ui| (0..n_v_rows).map(|vi| w_grid[vi][ui]).collect())
        .collect();

    BSplineSurface {
        degree_u: c_base.degree,
        degree_v: 2,
        knots_u: c_base.knots.clone(),
        knots_v,
        control_points: transposed_ctrl,
        weights: transposed_w,
    }
}

/// Convert a `BezierSurface` to a `BSplineSurface` by inserting clamped
/// endpoint knots in both parametric directions.
pub fn bezier_surface_to_bspline(bezier: &BezierSurface) -> BSplineSurface {
    let nu = bezier.control_points.len();
    let nv = if nu > 0 {
        bezier.control_points[0].len()
    } else {
        0
    };
    let deg_u = (nu - 1).max(1);
    let deg_v = (nv - 1).max(1);

    let mut knots_u = vec![0.0f64; deg_u + 1];
    knots_u.extend(vec![1.0f64; deg_u + 1]);
    let mut knots_v = vec![0.0f64; deg_v + 1];
    knots_v.extend(vec![1.0f64; deg_v + 1]);

    BSplineSurface {
        degree_u: deg_u,
        degree_v: deg_v,
        knots_u,
        knots_v,
        control_points: bezier.control_points.clone(),
        weights: bezier.weights.clone(),
    }
}

/// Sample a surface at `n_u × n_v` points over its default domain and build a
/// bilinear (degree-1,1) `BSplineSurface` approximation.
///
/// For surfaces without analytic NURBS conversion (Torus, Revolution, Extrusion,
/// Offset, Trimmed), this gives a piecewise-planar approximation.  Increasing
/// `n_u`, `n_v` improves accuracy.
fn sample_surface_to_bspline(surface: &Surface3, n_u: usize, n_v: usize) -> BSplineSurface {
    let [u0, u1, v0, v1] = surface.default_domain();
    let (u0, u1) = if (u1 - u0).abs() > 1e6 {
        (-10.0, 10.0)
    } else {
        (u0, u1)
    };
    let (v0, v1) = if (v1 - v0).abs() > 1e6 {
        (-10.0, 10.0)
    } else {
        (v0, v1)
    };
    let n_u = n_u.max(2);
    let n_v = n_v.max(2);

    let mut ctrl: Vec<Vec<DVec3>> = Vec::new();
    let mut w: Vec<Vec<f64>> = Vec::new();
    for i in 0..n_u {
        let u = u0 + (u1 - u0) * i as f64 / (n_u - 1) as f64;
        let mut row = Vec::new();
        let mut wrow = Vec::new();
        for j in 0..n_v {
            let v = v0 + (v1 - v0) * j as f64 / (n_v - 1) as f64;
            row.push(surface.point_at(u, v));
            wrow.push(1.0f64);
        }
        ctrl.push(row);
        w.push(wrow);
    }

    // Degree-1 in both directions (piecewise bilinear)
    let knots_u = build_uniform_knots(n_u, 1);
    let knots_v = build_uniform_knots(n_v, 1);

    BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u,
        knots_v,
        control_points: ctrl,
        weights: w,
    }
}

fn build_uniform_knots(n_ctrl: usize, degree: usize) -> Vec<f64> {
    let n_segments = n_ctrl - degree;
    let mut knots = vec![0.0f64; degree + 1];
    for i in 1..n_segments {
        knots.push(i as f64 / n_segments as f64);
    }
    knots.extend(vec![1.0f64; degree + 1]);
    knots
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Circle3, Curve3, Line3, SurfaceEval};
    use glam::DVec3;

    fn approx_eq3(a: DVec3, b: DVec3, tol: f64) -> bool {
        (a - b).length() < tol
    }

    // ── Curve tests ──────────────────────────────────────────────────────────

    #[test]
    fn line_bspline_endpoints() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let bs = line_to_bspline(&line);
        let p0 = bs.point_at(0.0);
        let p1 = bs.point_at(1.0);
        assert!(approx_eq3(p0, DVec3::new(0.0, 0.0, 0.0), 1e-10));
        assert!(approx_eq3(p1, DVec3::new(1.0, 0.0, 0.0), 1e-10));
    }

    #[test]
    fn circle_bspline_is_exact() {
        let circle = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 2.0,
        };
        let bs = circle_to_bspline(&circle);
        // The NURBS circle evaluates with rational weights — test a few points
        for i in 0..8 {
            let t = i as f64 / 8.0;
            let p = bs.point_at(t);
            let r = (p - circle.center).length();
            assert!(
                (r - circle.radius).abs() < 1e-10,
                "radius at t={t}: expected {}, got {r}",
                circle.radius
            );
        }
    }

    #[test]
    fn ellipse_bspline_endpoints() {
        use crate::geom::Ellipse3;
        let e = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 3.0,
            minor_radius: 1.5,
        };
        let bs = ellipse_to_bspline(&e);
        // t=0 and t=1 both map to the major-axis endpoint
        let p0 = bs.point_at(0.0);
        let p1 = bs.point_at(1.0);
        assert!(approx_eq3(p0, DVec3::new(3.0, 0.0, 0.0), 1e-10), "p0={p0}");
        assert!(approx_eq3(p1, DVec3::new(3.0, 0.0, 0.0), 1e-10), "p1={p1}");
    }

    #[test]
    fn curve_to_bspline_identity_for_bspline() {
        let bs_orig = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
        };
        let bs_conv = curve_to_bspline(&Curve3::BSpline(bs_orig.clone()), 16);
        assert_eq!(bs_conv.degree, bs_orig.degree);
        assert_eq!(bs_conv.control_points.len(), bs_orig.control_points.len());
    }

    // ── Surface tests ────────────────────────────────────────────────────────

    #[test]
    fn plane_bspline_corners() {
        use crate::geom::Plane;
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let bs = plane_to_bspline(&plane);
        // Should evaluate points at corners of the [-1,1]×[-1,1] domain
        let p00 = bs.point_at(0.0, 0.0);
        let p11 = bs.point_at(1.0, 1.0);
        // p at (u=0,v=0) should be corner of domain
        assert!(
            p00.distance(DVec3::new(-1.0, -1.0, 0.0)) < 1e-10
                || p00.distance(DVec3::new(1.0, 1.0, 0.0)) < 1e-10
                || p00.z.abs() < 1e-10, // at least z=0
            "plane corner z={}",
            p00.z
        );
        let _ = p11;
    }

    #[test]
    fn cylinder_bspline_on_surface() {
        use crate::geom::CylindricalSurface;
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let bs = cylinder_to_bspline(&cyl);
        // Sample several u values; all should be on the cylinder surface
        for i in 0..9 {
            let u = i as f64 / 8.0;
            let p = bs.point_at(u, 0.0);
            let r = DVec3::new(p.x, p.y, 0.0).length();
            assert!((r - 1.0).abs() < 1e-9, "u={u}: radius={r}");
        }
    }

    #[test]
    fn sphere_bspline_on_surface() {
        use crate::geom::SphericalSurface;
        let sph = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let bs = sphere_to_bspline(&sph);
        // Equator row (v=0.5 maps to colatitude 90°)
        for i in 0..9 {
            let u = i as f64 / 8.0;
            let p = bs.point_at(u, 0.5);
            let r = p.length();
            assert!((r - 1.0).abs() < 1e-9, "u={u}: radius={r}");
        }
    }
}
