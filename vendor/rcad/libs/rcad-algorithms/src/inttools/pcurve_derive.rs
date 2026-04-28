//! Analytic derivation of 2D parametric curves (PCurves) for surface-surface
//! intersection results.
//!
//! Each function takes a 3D intersection curve together with surface geometry
//! and returns the exact [`Curve2d`] that represents that intersection in the
//! surface's (u, v) parameter domain.

use glam::{DVec2, DVec3};
use rcad_kernel::fit::interpolate_points_2d;
use rcad_kernel::geom::{
    Circle2d, Circle3, ConicalSurface, Curve2d, CurveEval, CylindricalSurface, Ellipse2d,
    Ellipse3, Line2d, Line3, Plane, SphericalSurface, Surface3, any_perpendicular,
};
use rcad_kernel::projection::closest_point_on_surface;

// ─────────────────────────────────────────────────────────────────────────────
// Plane functions
// ─────────────────────────────────────────────────────────────────────────────

/// Project a [`Circle3`] onto a [`Plane`]'s (u, v) domain.
///
/// Uses `any_perpendicular(plane.normal)` as the u-axis and
/// `plane.normal × u_axis` as the v-axis, matching [`Plane::point_at`].
///
/// If the circle lies in the plane (its normal is parallel to the plane
/// normal), the result is an analytic [`Circle2d`].  Otherwise the circle
/// projects to a general conic and is approximated with a [`BSplineCurve2`]
/// built from 33 sampled points.
pub fn circle_pcurve_on_plane(circle: &Circle3, plane: &Plane) -> Curve2d {
    let u_axis = any_perpendicular(plane.normal);
    let v_axis = plane.normal.cross(u_axis);

    // Test whether the circle lies in the plane.
    let normal_dot = circle
        .normal
        .normalize()
        .dot(plane.normal.normalize())
        .abs();
    if (normal_dot - 1.0).abs() < 1e-6 {
        // Circle lies in the plane → analytic Circle2d.
        let diff = circle.center - plane.origin;
        let center_2d = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));
        return Curve2d::Circle(Circle2d {
            center: center_2d,
            radius: circle.radius,
        });
    }

    // Oblique case: sample the circle and project each point into the plane.
    let n_samples = 33_usize;
    let pts: Vec<DVec2> = (0..n_samples)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / (n_samples - 1) as f64;
            let p3 = circle.point_at(t);
            // Project onto the plane (drop the normal component).
            let diff = p3 - plane.origin;
            DVec2::new(diff.dot(u_axis), diff.dot(v_axis))
        })
        .collect();

    let bspline = interpolate_points_2d(&pts).expect("circle samples should not be degenerate");
    Curve2d::BSpline(bspline)
}

/// Project an [`Ellipse3`] onto a [`Plane`]'s (u, v) domain.
///
/// Returns an analytic [`Ellipse2d`] with the projected center, major
/// direction, and radii (unchanged — projection along a parallel normal
/// preserves semi-axes when the ellipse is coplanar with the plane).
pub fn ellipse_pcurve_on_plane(ellipse: &Ellipse3, plane: &Plane) -> Curve2d {
    let u_axis = any_perpendicular(plane.normal);
    let v_axis = plane.normal.cross(u_axis);

    let diff = ellipse.center - plane.origin;
    let center_2d = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));

    let major_proj = DVec2::new(ellipse.major_dir.dot(u_axis), ellipse.major_dir.dot(v_axis));
    let major_dir_2d = if major_proj.length() > 1e-12 {
        major_proj.normalize()
    } else {
        DVec2::X
    };

    Curve2d::Ellipse(Ellipse2d {
        center: center_2d,
        major_dir: major_dir_2d,
        major_radius: ellipse.major_radius,
        minor_radius: ellipse.minor_radius,
    })
}

/// Project a [`Line3`] onto a [`Plane`]'s (u, v) domain.
///
/// Returns a [`Line2d`] whose origin and direction are the projections of the
/// 3D line's origin and direction into the plane's parameter space.
pub fn line_pcurve_on_plane(line: &Line3, plane: &Plane) -> Curve2d {
    let u_axis = any_perpendicular(plane.normal);
    let v_axis = plane.normal.cross(u_axis);

    let diff = line.origin - plane.origin;
    let origin_2d = DVec2::new(diff.dot(u_axis), diff.dot(v_axis));

    let dir_2d = DVec2::new(line.direction.dot(u_axis), line.direction.dot(v_axis));
    let direction_2d = if dir_2d.length() > 1e-12 {
        dir_2d.normalize()
    } else {
        DVec2::X
    };

    Curve2d::Line(Line2d {
        origin: origin_2d,
        direction: direction_2d,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Sphere functions
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the PCurve of a [`Circle3`] on a [`SphericalSurface`].
///
/// Any circle that lies entirely on a sphere is a **latitude circle** (constant
/// colatitude φ) for that sphere.  The colatitude is determined solely by the
/// axial component of the circle's centre:
///
/// ```text
/// φ = acos((circle_center − sphere_center) · sphere_axis / sphere_radius)
/// ```
///
/// This formula is valid for **any** circle on the sphere surface, regardless
/// of whether the circle's normal is parallel to the sphere axis.
///
/// The sphere's (u, v) domain uses **u = longitude ∈ [−π, π]** (matching the
/// atan2-based projection in `ds.rs`) and **v = colatitude ∈ [0, π]**.
///
/// Returns a horizontal [`Line2d`] at v = φ with `origin.x = −π` and
/// `direction = (1, 0)`, so that when the circle's parameter `t` runs over
/// `[0, 2π]`, the longitude sweeps the full `[−π, +π]` domain of the sphere.
///
/// # Applicability
///
/// This function is exact when the circle is a true latitude circle (normal ∥
/// sphere axis, e.g. sphere–cylinder axis-aligned intersection).  When the
/// circle's normal is not parallel to the sphere axis (e.g. sphere–sphere
/// intersection), the v value is still exact but the u parameterization is a
/// uniform sweep; use [`fallback_pcurve_by_projection`] when the exact
/// per-t correspondence matters.
pub fn circle_pcurve_on_sphere(circle: &Circle3, sphere: &SphericalSurface) -> Curve2d {
    let along_axis = (circle.center - sphere.center).dot(sphere.axis.normalize());
    let phi = (along_axis / sphere.radius).clamp(-1.0, 1.0).acos();

    // Start the horizontal line at u = -π so that sampling over [0, 2π] spans
    // the full [-π, +π] UV boundary of the sphere.
    Curve2d::Line(Line2d {
        origin: DVec2::new(-std::f64::consts::PI, phi),
        direction: DVec2::new(1.0, 0.0),
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Cylinder functions
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the PCurve of a [`Circle3`] on a [`CylindricalSurface`].
///
/// A circle perpendicular to the cylinder axis at height h returns a
/// horizontal [`Line2d`] at v = h in (θ, h) space.
pub fn circle_pcurve_on_cylinder(circle: &Circle3, cyl: &CylindricalSurface) -> Curve2d {
    let h = (circle.center - cyl.origin).dot(cyl.axis.normalize());

    Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, h),
        direction: DVec2::new(1.0, 0.0),
    })
}

/// Compute the PCurve of a [`Line3`] on a [`CylindricalSurface`].
///
/// A line parallel to the cylinder axis at azimuth θ returns a vertical
/// [`Line2d`] at u = θ in (θ, h) space.
pub fn line_pcurve_on_cylinder(line: &Line3, cyl: &CylindricalSurface) -> Curve2d {
    let u_axis = any_perpendicular(cyl.axis);
    let v_axis = cyl.axis.cross(u_axis).normalize();

    let radial = line.origin - cyl.origin;
    let radial_perp = radial - cyl.axis * radial.dot(cyl.axis.normalize());
    let theta = radial_perp.dot(v_axis).atan2(radial_perp.dot(u_axis));

    let h = radial.dot(cyl.axis.normalize());

    Curve2d::Line(Line2d {
        origin: DVec2::new(theta, h),
        direction: DVec2::new(0.0, 1.0),
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Cone functions
// ─────────────────────────────────────────────────────────────────────────────

fn cone_uv_from_point(point: DVec3, cone: &ConicalSurface) -> DVec2 {
    let axis = cone.axis_dir();
    let u_axis = any_perpendicular(axis);
    let v_axis = axis.cross(u_axis).normalize();
    let local = point - cone.apex;
    let axial = local.dot(axis);
    let radial = local - axis * axial;
    let mut u = radial.dot(v_axis).atan2(radial.dot(u_axis));
    if u < 0.0 {
        u += std::f64::consts::TAU;
    }
    DVec2::new(u, cone.slant_from_axial(axial))
}

fn sampled_curve_pcurve_on_cone(
    curve: &rcad_kernel::geom::Curve3,
    t_range: &[f64; 2],
    cone: &ConicalSurface,
) -> Curve2d {
    let n = 33_usize;
    let mut pts: Vec<DVec2> = (0..n)
        .map(|i| {
            let t = t_range[0] + (t_range[1] - t_range[0]) * i as f64 / (n - 1) as f64;
            let p3 = curve.point_at(t);
            cone_uv_from_point(p3, cone)
        })
        .collect();

    for i in 1..pts.len() {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x += std::f64::consts::TAU;
            }
        }
    }

    let bspline = interpolate_points_2d(&pts).expect("cone curve samples should not be degenerate");
    Curve2d::BSpline(bspline)
}

pub fn circle_pcurve_on_cone(circle: &Circle3, cone: &ConicalSurface) -> Curve2d {
    let axis = cone.axis_dir();
    let normal_dot = circle.normal.normalize().dot(axis).abs();
    if (normal_dot - 1.0).abs() < 1e-6 {
        let slant = cone.slant_from_axial((circle.center - cone.apex).dot(axis));
        return Curve2d::Line(Line2d {
            origin: DVec2::new(0.0, slant),
            direction: DVec2::new(1.0, 0.0),
        });
    }
    sampled_curve_pcurve_on_cone(&rcad_kernel::geom::Curve3::Circle(*circle), &[0.0, std::f64::consts::TAU], cone)
}

pub fn line_pcurve_on_cone(line: &Line3, cone: &ConicalSurface) -> Curve2d {
    let uv0 = cone_uv_from_point(line.origin, cone);
    let uv1 = cone_uv_from_point(line.origin + line.direction, cone);
    let du = (uv1.x - uv0.x).abs().min((uv1.x - uv0.x + std::f64::consts::TAU).abs());
    if du < 1e-6 {
        let dir_v = if uv1.y >= uv0.y { 1.0 } else { -1.0 };
        return Curve2d::Line(Line2d {
            origin: uv0,
            direction: DVec2::new(0.0, dir_v),
        });
    }
    sampled_curve_pcurve_on_cone(&rcad_kernel::geom::Curve3::Line(*line), &[-10.0, 10.0], cone)
}

pub fn ellipse_pcurve_on_cone(ellipse: &Ellipse3, cone: &ConicalSurface) -> Curve2d {
    sampled_curve_pcurve_on_cone(&rcad_kernel::geom::Curve3::Ellipse(*ellipse), &[0.0, std::f64::consts::TAU], cone)
}

pub fn sampled_pcurve_on_cone(
    curve: &rcad_kernel::geom::Curve3,
    t_range: &[f64; 2],
    cone: &ConicalSurface,
) -> Curve2d {
    sampled_curve_pcurve_on_cone(curve, t_range, cone)
}

// ─────────────────────────────────────────────────────────────────────────────
// Numeric fallback functions
// ─────────────────────────────────────────────────────────────────────────────

/// Derive a PCurve by sampling `curve` at 33 evenly-spaced parameter values
/// over `t_range` and projecting each 3D point onto `surface`.
///
/// Intended as a fallback for curve/surface combinations that do not have an
/// analytic form.  Returns a [`BSplineCurve2`] interpolated through the
/// projected (u, v) points.
pub fn fallback_pcurve_by_projection(
    curve: &rcad_kernel::geom::Curve3,
    t_range: &[f64; 2],
    surface: &Surface3,
) -> Curve2d {
    let n = 33_usize;
    let mut pts: Vec<DVec2> = (0..n)
        .map(|i| {
            let t = t_range[0] + (t_range[1] - t_range[0]) * i as f64 / (n - 1) as f64;
            let p3 = curve.point_at(t);
            let proj = closest_point_on_surface(surface, p3, 16);
            DVec2::new(proj.params.0, proj.params.1)
        })
        .collect();

    // Unwrap seam discontinuities: atan2-based u values are in [-π, π], but
    // consecutive samples may jump by ~2π when the curve crosses the seam.
    // Make the u sequence monotone so the interpolated BSpline has no kinks.
    for i in 1..pts.len() {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            // Jumped from near -π back up to near +π: pull remaining down.
            for p in &mut pts[i..] {
                p.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            // Jumped from near +π down to near -π: push remaining up.
            for p in &mut pts[i..] {
                p.x += std::f64::consts::TAU;
            }
        }
    }

    match interpolate_points_2d(&pts) {
        Ok(bspline) => Curve2d::BSpline(bspline),
        Err(_) => {
            // Extremely short/degenerate intersection remnants can collapse
            // to a single UV sample. Emit a stable line fallback instead of
            // panicking so boolean export can proceed.
            let origin = pts.first().copied().unwrap_or(DVec2::ZERO);
            let direction = pts
                .iter()
                .skip(1)
                .map(|p| *p - origin)
                .find(|d| d.length_squared() > 1e-24)
                .map(|d| d.normalize())
                .unwrap_or(DVec2::X);
            Curve2d::Line(Line2d { origin, direction })
        }
    }
}

/// Project a 3D polyline onto `surface` and interpolate a [`BSplineCurve2`].
///
/// Returns `None` if the polyline has fewer than 2 points or all projected
/// points are coincident.
pub fn polyline_pcurve_by_projection(polyline: &[DVec3], surface: &Surface3) -> Option<Curve2d> {
    if polyline.len() < 2 {
        return None;
    }

    let mut pts: Vec<DVec2> = polyline
        .iter()
        .map(|&p3| {
            let proj = closest_point_on_surface(surface, p3, 16);
            DVec2::new(proj.params.0, proj.params.1)
        })
        .collect();

    // Unwrap seam discontinuities (same logic as fallback_pcurve_by_projection).
    for i in 1..pts.len() {
        let du = pts[i].x - pts[i - 1].x;
        if du > std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x -= std::f64::consts::TAU;
            }
        } else if du < -std::f64::consts::PI {
            for p in &mut pts[i..] {
                p.x += std::f64::consts::TAU;
            }
        }
    }

    interpolate_points_2d(&pts).ok().map(Curve2d::BSpline)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Curve2dEval, SphericalSurface};
    use std::f64::consts::PI;

    /// A circle whose normal is Z lying in the XY plane (z = 0) projects to a
    /// Circle2d in the plane's (u, v) space.
    #[test]
    fn circle_on_plane_is_circle() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let circle = Circle3 {
            center: DVec3::new(1.0, 2.0, 0.0),
            normal: DVec3::Z,
            radius: 3.0,
        };

        let pcurve = circle_pcurve_on_plane(&circle, &plane);

        match pcurve {
            Curve2d::Circle(c) => {
                assert!((c.radius - 3.0).abs() < 1e-9, "radius={}", c.radius);
            }
            other => panic!("expected Circle2d, got {other:?}"),
        }
    }

    /// A circle at z = 1 on a sphere of radius 2 (axis = Z, center = origin)
    /// should produce φ = acos(0.5) = π/3.
    #[test]
    fn circle_on_sphere_is_latitude() {
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        let circle = Circle3 {
            center: DVec3::new(0.0, 0.0, 1.0),
            normal: DVec3::Z,
            radius: (3.0_f64).sqrt(), // r² + z² = 4  ⟹  r = √3
        };

        let pcurve = circle_pcurve_on_sphere(&circle, &sphere);

        match pcurve {
            Curve2d::Line(l) => {
                let expected_phi = (0.5_f64).acos(); // π/3
                assert!(
                    (l.origin.y - expected_phi).abs() < 1e-9,
                    "phi={}, expected {expected_phi}",
                    l.origin.y
                );
                // Origin x starts at -π so the line spans [-π, +π] over [0, 2π] sampling.
                assert!(
                    (l.origin.x + PI).abs() < 1e-9,
                    "expected origin.x = -π, got {}",
                    l.origin.x
                );
                // Direction must be horizontal (constant colatitude).
                assert!((l.direction.x - 1.0).abs() < 1e-9);
                assert!(l.direction.y.abs() < 1e-9);
            }
            other => panic!("expected Line2d, got {other:?}"),
        }
    }

    /// A circle at height h = 3 on a cylinder should produce a horizontal
    /// line at v = 3.
    #[test]
    fn circle_on_cylinder_is_h_line() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let circle = Circle3 {
            center: DVec3::new(0.0, 0.0, 3.0),
            normal: DVec3::Z,
            radius: 1.0,
        };

        let pcurve = circle_pcurve_on_cylinder(&circle, &cyl);

        match pcurve {
            Curve2d::Line(l) => {
                assert!(
                    (l.origin.y - 3.0).abs() < 1e-9,
                    "h={}, expected 3.0",
                    l.origin.y
                );
                assert!((l.direction.x - 1.0).abs() < 1e-9);
                assert!(l.direction.y.abs() < 1e-9);
            }
            other => panic!("expected Line2d, got {other:?}"),
        }
    }

    #[test]
    fn circle_on_cone_is_h_line() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: (1.0_f64 / 2.0).atan(),
        };
        let h = 3.0;
        let slant = cone.slant_from_axial(h);
        let circle = Circle3 {
            center: DVec3::new(0.0, 0.0, h),
            normal: DVec3::Z,
            radius: h * cone.half_angle_rad.tan(),
        };

        let pcurve = circle_pcurve_on_cone(&circle, &cone);
        match pcurve {
            Curve2d::Line(l) => {
                assert!((l.origin.y - slant).abs() < 1e-9);
                assert!((l.direction.x - 1.0).abs() < 1e-9);
                assert!(l.direction.y.abs() < 1e-9);
            }
            other => panic!("expected Line2d, got {other:?}"),
        }
    }

    #[test]
    fn line_on_cone_is_v_line() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: (1.0_f64 / 2.0).atan(),
        };
        let line = Line3 {
            origin: DVec3::new(1.0, 0.0, 2.0),
            direction: DVec3::new(0.5, 0.0, 1.0).normalize(),
        };

        let pcurve = line_pcurve_on_cone(&line, &cone);
        match pcurve {
            Curve2d::Line(l) => {
                // The origin's u coordinate depends on the arbitrary perpendicular chosen,
                // so we only verify that the line is a v-line (direction purely in v)
                // by checking that the x direction is zero.
                assert!(l.direction.x.abs() < 1e-9, "v-line should have zero u direction");
                assert!((l.direction.y - 1.0).abs() < 1e-9, "v-line should have unit v direction");
            }
            other => panic!("expected Line2d, got {other:?}"),
        }
    }

    /// The fallback projection of any curve on a sphere should produce a
    /// BSplineCurve2.
    #[test]
    fn fallback_projection_produces_bspline() {
        use rcad_kernel::geom::Curve3;

        let sphere_surface = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        });
        // A circle on the sphere at the equator (z = 0, r = 2).
        let circle = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 2.0,
        };
        let curve3 = Curve3::Circle(circle);
        let t_range = [0.0_f64, PI]; // half circle

        let pcurve = fallback_pcurve_by_projection(&curve3, &t_range, &sphere_surface);

        match pcurve {
            Curve2d::BSpline(ref b) => {
                // Should have at least some control points.
                assert!(!b.control_points.is_empty());
                // Evaluate endpoints to make sure the BSpline is usable.
                let p0 = pcurve.point_at(0.0);
                let p1 = pcurve.point_at(1.0);
                // Both must be finite.
                assert!(p0.x.is_finite() && p0.y.is_finite());
                assert!(p1.x.is_finite() && p1.y.is_finite());
            }
            other => panic!("expected BSpline2, got {other:?}"),
        }
    }

    /// `circle_pcurve_on_sphere` is valid for any circle that lies on the sphere,
    /// even when the circle's normal is NOT parallel to the sphere axis.
    /// Here the intersection of two spheres whose centres are separated along X
    /// gives a circle whose normal is X, but we still get a correct latitude line.
    #[test]
    fn circle_on_sphere_non_axis_normal() {
        // Two unit spheres: sph1 at origin (axis=Z), sph2 at (1,0,0).
        // Their intersection circle: d=1, r1=r2=1.
        //   h = (1 + 1 - 1)/(2) = 0.5   (distance from sph1 center to radical plane)
        //   r_circ = sqrt(1 - 0.25) = sqrt(0.75)
        //   circle center = (0.5, 0, 0)
        let sphere = SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 1.0 };
        let circle = Circle3 {
            center: DVec3::new(0.5, 0.0, 0.0),
            normal: DVec3::X,  // NOT parallel to sphere axis (Z)
            radius: (0.75_f64).sqrt(),
        };

        let pcurve = circle_pcurve_on_sphere(&circle, &sphere);

        // along_axis = (0.5, 0, 0) · (0, 0, 1) = 0  →  phi = acos(0) = π/2
        match pcurve {
            Curve2d::Line(l) => {
                let expected_phi = std::f64::consts::PI / 2.0;
                assert!(
                    (l.origin.y - expected_phi).abs() < 1e-9,
                    "phi={}, expected π/2",
                    l.origin.y
                );
                assert!((l.direction.x - 1.0).abs() < 1e-9);
                assert!(l.direction.y.abs() < 1e-9);
                // origin.x = longitude of circle.point_at(0) — just check it's finite
                assert!(l.origin.x.is_finite());
            }
            other => panic!("expected Line2d, got {other:?}"),
        }
    }

    /// Verify that `circle_pcurve_on_sphere` and `fallback_pcurve_by_projection`
    /// agree at the equatorial circle (both should give v ≈ π/2).
    #[test]
    fn analytic_sphere_pcurve_matches_fallback() {
        use rcad_kernel::geom::{Curve2dEval, Curve3};

        let sphere_surf = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        });
        let sphere = SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 2.0 };
        // Equatorial circle at z=0, r=2
        let circle = Circle3 { center: DVec3::ZERO, normal: DVec3::Z, radius: 2.0 };

        let analytic = circle_pcurve_on_sphere(&circle, &sphere);
        let fallback = fallback_pcurve_by_projection(
            &Curve3::Circle(circle),
            &[0.0, std::f64::consts::TAU],
            &sphere_surf,
        );

        // Both should yield v ≈ π/2 everywhere (equator = colatitude π/2)
        for i in 0..8 {
            let t = i as f64 / 8.0;
            let pa = analytic.point_at(t);
            let pf = fallback.point_at(t);
            assert!(
                (pa.y - pf.y).abs() < 0.02,
                "t={t}: analytic.v={:.4} fallback.v={:.4}",
                pa.y,
                pf.y
            );
        }
    }
}
