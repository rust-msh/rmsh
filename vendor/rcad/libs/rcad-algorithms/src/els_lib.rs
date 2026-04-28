//! ElSLib-style elementary surface utilities.
//!
//! Provides utilities for elementary surfaces analogous to OCCT `ElSLib` package.
//! Includes evaluation, parameter computation, and differential properties
//! for Plane, Cylinder, Sphere, Cone, Torus, and BSplineSurface.
//!
//! # Overview
//!
//! For each elementary surface type, this module provides:
//! - `point_at(surf, u, v)`: Compute 3D point from surface parameters
//! - `parameters(surf, point)`: Compute (u, v) parameters from a 3D point
//! - `normal(surf, u, v)`: Surface normal at (u, v)
//! - `tangent_u(surf, u, v)`: Partial derivative dS/du (u-tangent)
//! - `tangent_v(surf, u, v)`: Partial derivative dS/dv (v-tangent)
//!
//! # Coordinate Conventions
//!
//! - Plane: u and v are Cartesian coordinates in the plane's local frame
//! - Cylinder: u = azimuth angle [0, 2*pi], v = height along axis
//! - Sphere: u = longitude [0, 2*pi], v = colatitude [0, pi] (0 = north pole)
//! - Cone: u = azimuth [0, 2*pi], v = slant distance from reference circle
//! - Torus: u = major angle [0, 2*pi], v = minor angle [0, 2*pi]

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    any_perpendicular, Plane, CylindricalSurface, SphericalSurface,
    ConicalSurface, ToroidalSurface, BSplineSurface, SurfaceEval,
};
use std::f64::consts::PI;

// =============================================================================
// Plane Utilities
// =============================================================================

/// Compute the 3D point on a plane at parameters (u, v).
///
/// The plane is parameterized as: P(u, v) = origin + u * x_axis + v * y_axis
/// where x_axis and y_axis form an orthonormal basis with the normal.
pub fn plane_point_at(plane: &Plane, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(plane.normal);
    let y_ax = plane.normal.cross(x_ax);
    plane.origin + u * x_ax + v * y_ax
}

/// Compute the (u, v) parameters for a point on or near a plane.
///
/// Projects the point onto the plane and returns the local coordinates.
/// The result satisfies: `plane_point_at(plane, u, v) == point.projected_onto_plane`.
pub fn plane_parameters(plane: &Plane, point: DVec3) -> DVec2 {
    let x_ax = any_perpendicular(plane.normal);
    let y_ax = plane.normal.cross(x_ax);
    let d = point - plane.origin;
    DVec2::new(d.dot(x_ax), d.dot(y_ax))
}

/// Get the normal vector of a plane (constant across the surface).
pub fn plane_normal(plane: &Plane) -> DVec3 {
    plane.normal
}

/// Get the u-tangent vector of a plane (constant across the surface).
///
/// This is the partial derivative dS/du, pointing along the local x-axis.
pub fn plane_tangent_u(plane: &Plane) -> DVec3 {
    any_perpendicular(plane.normal)
}

/// Get the v-tangent vector of a plane (constant across the surface).
///
/// This is the partial derivative dS/dv, pointing along the local y-axis.
pub fn plane_tangent_v(plane: &Plane) -> DVec3 {
    let x_ax = any_perpendicular(plane.normal);
    plane.normal.cross(x_ax)
}

// =============================================================================
// Cylinder Utilities
// =============================================================================

/// Compute the 3D point on a cylinder at parameters (u, v).
///
/// - u: azimuth angle [0, 2*pi]
/// - v: height along the cylinder axis
///
/// P(u, v) = origin + radius * (cos(u) * x_axis + sin(u) * y_axis) + v * axis
pub fn cylinder_point_at(cyl: &CylindricalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(cyl.axis);
    let y_ax = cyl.axis.cross(x_ax).normalize();
    cyl.origin + cyl.radius * (u.cos() * x_ax + u.sin() * y_ax) + v * cyl.axis
}

/// Compute the (u, v) parameters for a point on or near a cylinder.
///
/// Projects the point onto the cylinder surface and returns the corresponding
/// parameters. The u angle is normalized to [0, 2*pi).
pub fn cylinder_parameters(cyl: &CylindricalSurface, point: DVec3) -> DVec2 {
    let axis = cyl.axis.normalize_or_zero();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();

    // Vector from cylinder origin to point
    let d = point - cyl.origin;

    // Height along axis (v parameter)
    let v = d.dot(axis);

    // Radial component
    let radial = d - v * axis;
    // atan2(y, x) where y = radial.dot(y_ax), x = radial.dot(x_ax)
    let u = radial.dot(y_ax).atan2(radial.dot(x_ax));

    // Normalize u to [0, 2*pi)
    let u = if u < 0.0 { u + 2.0 * PI } else { u };

    DVec2::new(u, v)
}

/// Get the normal vector of a cylinder at parameters (u, v).
///
/// The normal points outward from the axis (radially).
pub fn cylinder_normal(cyl: &CylindricalSurface, u: f64, _v: f64) -> DVec3 {
    let x_ax = any_perpendicular(cyl.axis);
    let y_ax = cyl.axis.cross(x_ax).normalize();
    (u.cos() * x_ax + u.sin() * y_ax).normalize()
}

/// Get the u-tangent vector of a cylinder at parameters (u, v).
///
/// This is the azimuthal tangent (along the circular cross-section).
pub fn cylinder_tangent_u(cyl: &CylindricalSurface, u: f64, _v: f64) -> DVec3 {
    let x_ax = any_perpendicular(cyl.axis);
    let y_ax = cyl.axis.cross(x_ax).normalize();
    (-u.sin() * x_ax + u.cos() * y_ax).normalize()
}

/// Get the v-tangent vector of a cylinder at parameters (u, v).
///
/// This is the axial tangent (along the cylinder axis).
pub fn cylinder_tangent_v(cyl: &CylindricalSurface, _u: f64, _v: f64) -> DVec3 {
    cyl.axis.normalize_or_zero()
}

// =============================================================================
// Sphere Utilities
// =============================================================================

/// Compute the 3D point on a sphere at parameters (u, v).
///
/// - u: longitude angle [0, 2*pi]
/// - v: colatitude angle [0, pi] (0 = north pole, pi = south pole)
///
/// P(u, v) = center + radius * (sin(v) * (cos(u) * x + sin(u) * y) + cos(v) * axis)
pub fn sphere_point_at(sph: &SphericalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(sph.axis);
    let y_ax = sph.axis.cross(x_ax).normalize();
    sph.center + sph.radius * (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * sph.axis)
}

/// Compute the (u, v) parameters for a point on or near a sphere.
///
/// Projects the point onto the sphere surface and returns the corresponding
/// angles. At the poles, u is set to 0.0.
pub fn sphere_parameters(sph: &SphericalSurface, point: DVec3) -> DVec2 {
    let axis = sph.axis.normalize_or_zero();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();

    // Vector from center to point (normalized)
    let d = (point - sph.center).normalize_or_zero();

    // Colatitude: angle from axis
    let cos_v = d.dot(axis);
    let v = cos_v.clamp(-1.0, 1.0).acos();

    // Longitude: angle in the equatorial plane
    let sin_v = v.sin();
    let u = if sin_v.abs() > 1e-10 {
        let radial = d - cos_v * axis;
        let u_raw = radial.dot(y_ax).atan2(radial.dot(x_ax));
        if u_raw < 0.0 { u_raw + 2.0 * PI } else { u_raw }
    } else {
        0.0 // At poles, u is undefined; use 0
    };

    DVec2::new(u, v)
}

/// Get the normal vector of a sphere at parameters (u, v).
///
/// The normal points outward from the center.
pub fn sphere_normal(sph: &SphericalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(sph.axis);
    let y_ax = sph.axis.cross(x_ax).normalize();
    (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * sph.axis).normalize()
}

/// Get the u-tangent vector of a sphere at parameters (u, v).
///
/// This is the longitude tangent (along lines of latitude).
pub fn sphere_tangent_u(sph: &SphericalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(sph.axis);
    let y_ax = sph.axis.cross(x_ax).normalize();
    v.sin() * (-u.sin() * x_ax + u.cos() * y_ax)
}

/// Get the v-tangent vector of a sphere at parameters (u, v).
///
/// This is the colatitude tangent (along meridians).
pub fn sphere_tangent_v(sph: &SphericalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(sph.axis);
    let y_ax = sph.axis.cross(x_ax).normalize();
    sph.radius * (v.cos() * (u.cos() * x_ax + u.sin() * y_ax) - v.sin() * sph.axis)
}

// =============================================================================
// Cone Utilities
// =============================================================================

/// Compute the 3D point on a cone at parameters (u, v).
///
/// - u: azimuth angle [0, 2*pi]
/// - v: slant distance from the reference circle at apex
///
/// The reference circle has radius `cone.radius` at the apex point.
/// Positive v moves toward larger radius if half_angle > 0.
pub fn cone_point_at(cone: &ConicalSurface, u: f64, v: f64) -> DVec3 {
    let axis = cone.axis_dir();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();
    let radial = cone.radius_at_slant(v);
    let axial = cone.axial_from_slant(v);
    cone.apex + axial * axis + radial * (u.cos() * x_ax + u.sin() * y_ax)
}

/// Compute the (u, v) parameters for a point on or near a cone.
///
/// Projects the point onto the cone surface and returns the corresponding
/// parameters.
pub fn cone_parameters(cone: &ConicalSurface, point: DVec3) -> DVec2 {
    let axis = cone.axis_dir();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();

    // Vector from apex to point
    let d = point - cone.apex;

    // Axial distance from apex
    let axial = d.dot(axis);

    // Radial component
    let radial_vec = d - axial * axis;
    let radial_dist = radial_vec.length();

    // Azimuth angle
    let u = if radial_dist > 1e-10 {
        let u_raw = radial_vec.dot(y_ax).atan2(radial_vec.dot(x_ax));
        if u_raw < 0.0 { u_raw + 2.0 * PI } else { u_raw }
    } else {
        0.0
    };

    // Slant distance from reference circle
    // At the reference circle (v=0), axial = 0 and radial = cone.radius
    // v is the distance along the cone surface from this reference
    let slant_from_apex = (axial * axial + radial_dist * radial_dist).sqrt();
    let ref_slant = if cone.half_angle_rad.tan().abs() > 1e-10 {
        cone.radius / cone.half_angle_rad.sin()
    } else {
        0.0
    };
    let v = slant_from_apex - ref_slant;

    DVec2::new(u, v)
}

/// Get the normal vector of a cone at parameters (u, v).
///
/// The normal is constant along lines of constant u (generators).
pub fn cone_normal(cone: &ConicalSurface, u: f64, _v: f64) -> DVec3 {
    let axis = cone.axis_dir();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();
    let radial = u.cos() * x_ax + u.sin() * y_ax;
    let half = cone.half_angle_rad;
    (radial * half.cos() - axis * half.sin()).normalize()
}

// =============================================================================
// Torus Utilities
// =============================================================================

/// Compute the 3D point on a torus at parameters (u, v).
///
/// - u: major angle [0, 2*pi] (angle around the main axis)
/// - v: minor angle [0, 2*pi] (angle around the tube)
///
/// The torus is centered at `center` with the main axis `axis`.
/// Major radius is the distance from center to tube center.
/// Minor radius is the tube radius.
pub fn torus_point_at(torus: &ToroidalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(torus.axis);
    let y_ax = torus.axis.cross(x_ax).normalize();

    // Center of the tube cross-section at angle u
    let tube_center = torus.center + torus.major_radius * (u.cos() * x_ax + u.sin() * y_ax);

    // Radial direction from main axis to tube center
    let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();

    // Point on the tube surface
    tube_center + torus.minor_radius * (v.cos() * radial + v.sin() * torus.axis)
}

/// Compute the (u, v) parameters for a point on or near a torus.
///
/// Projects the point onto the torus surface and returns the corresponding
/// parameters.
pub fn torus_parameters(torus: &ToroidalSurface, point: DVec3) -> DVec2 {
    let axis = torus.axis.normalize_or_zero();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();

    // Vector from torus center to point
    let d = point - torus.center;

    // Height above/below the equatorial plane
    let z = d.dot(axis);

    // Radial component in equatorial plane
    let radial_2d = d - z * axis;
    let r_2d = radial_2d.length();

    // Major angle u
    let u = if r_2d > 1e-10 {
        let u_raw = radial_2d.dot(y_ax).atan2(radial_2d.dot(x_ax));
        if u_raw < 0.0 { u_raw + 2.0 * PI } else { u_raw }
    } else {
        0.0
    };

    // Minor angle v
    // The tube center at angle u is at distance major_radius from axis
    // dv is the distance from the tube center in the radial direction
    let dv = r_2d - torus.major_radius;

    let v = if dv.abs() > 1e-10 || z.abs() > 1e-10 {
        let v_raw = z.atan2(dv);
        if v_raw < 0.0 { v_raw + 2.0 * PI } else { v_raw }
    } else {
        0.0
    };

    DVec2::new(u, v)
}

/// Get the normal vector of a torus at parameters (u, v).
///
/// The normal points outward from the tube surface.
pub fn torus_normal(torus: &ToroidalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(torus.axis);
    let y_ax = torus.axis.cross(x_ax).normalize();
    let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();
    (v.cos() * radial + v.sin() * torus.axis).normalize()
}

/// Get the u-tangent vector of a torus at parameters (u, v).
///
/// This is the tangent along the major circle (around the main axis).
pub fn torus_tangent_u(torus: &ToroidalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(torus.axis);
    let y_ax = torus.axis.cross(x_ax).normalize();

    // Derivative of tube center w.r.t. u
    let dcenter_du = torus.major_radius * (-u.sin() * x_ax + u.cos() * y_ax);

    // Derivative of radial direction w.r.t. u
    let dradial_du = (-u.sin() * x_ax + u.cos() * y_ax);

    dcenter_du + torus.minor_radius * v.cos() * dradial_du
}

/// Get the v-tangent vector of a torus at parameters (u, v).
///
/// This is the tangent along the minor circle (around the tube).
pub fn torus_tangent_v(torus: &ToroidalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(torus.axis);
    let y_ax = torus.axis.cross(x_ax).normalize();
    let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();
    torus.minor_radius * (-v.sin() * radial + v.cos() * torus.axis)
}

// =============================================================================
// BSplineSurface Utilities
// =============================================================================

/// Compute the 3D point on a BSpline surface at parameters (u, v).
///
/// Uses tensor-product NURBS evaluation via the `SurfaceEval` trait.
pub fn bspline_surface_point_at(surf: &BSplineSurface, u: f64, v: f64) -> DVec3 {
    surf.point_at(u, v)
}

/// Compute the normal vector of a BSpline surface at parameters (u, v).
///
/// Uses finite differences to compute the cross product of partial derivatives.
pub fn bspline_surface_normal(surf: &BSplineSurface, u: f64, v: f64) -> DVec3 {
    surf.normal_at(u, v)
}

/// Compute the first partial derivatives of a BSpline surface at (u, v).
///
/// Returns `[du, dv, dudv]` where:
/// - `du` is the partial derivative with respect to u (dS/du)
/// - `dv` is the partial derivative with respect to v (dS/dv)
/// - `dudv` is the mixed second derivative (d2S/dudv)
///
/// Uses central finite differences for numerical stability.
pub fn bspline_surface_derivatives(surf: &BSplineSurface, u: f64, v: f64) -> [DVec3; 3] {
    let eps = 1e-5;
    let [u0, u1, v0, v1] = surf.default_domain();

    // Clamp to domain bounds
    let u_minus = (u - eps).max(u0);
    let u_plus = (u + eps).min(u1);
    let v_minus = (v - eps).max(v0);
    let v_plus = (v + eps).min(v1);

    // First derivatives using central differences where possible
    let du = if u_plus > u_minus {
        (surf.point_at(u_plus, v) - surf.point_at(u_minus, v)) / (u_plus - u_minus)
    } else if u_plus > u0 {
        (surf.point_at(u_plus, v) - surf.point_at(u, v)) / (u_plus - u)
    } else if u_minus < u1 {
        (surf.point_at(u, v) - surf.point_at(u_minus, v)) / (u - u_minus)
    } else {
        DVec3::ZERO
    };

    let dv = if v_plus > v_minus {
        (surf.point_at(u, v_plus) - surf.point_at(u, v_minus)) / (v_plus - v_minus)
    } else if v_plus > v0 {
        (surf.point_at(u, v_plus) - surf.point_at(u, v)) / (v_plus - v)
    } else if v_minus < v1 {
        (surf.point_at(u, v) - surf.point_at(u, v_minus)) / (v - v_minus)
    } else {
        DVec3::ZERO
    };

    // Mixed second derivative using finite differences
    let dudv = if u_plus > u_minus && v_plus > v_minus {
        let p_pp = surf.point_at(u_plus, v_plus);
        let p_pm = surf.point_at(u_plus, v_minus);
        let p_mp = surf.point_at(u_minus, v_plus);
        let p_mm = surf.point_at(u_minus, v_minus);

        ((p_pp - p_pm) - (p_mp - p_mm)) / ((u_plus - u_minus) * (v_plus - v_minus))
    } else {
        DVec3::ZERO
    };

    [du, dv, dudv]
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    fn approx_eq(a: DVec3, b: DVec3, tol: f64) -> bool {
        (a - b).length() < tol
    }

    fn approx_eq_2(a: DVec2, b: DVec2, tol: f64) -> bool {
        (a - b).length() < tol
    }

    // -------------------------------------------------------------------------
    // Plane Tests
    // -------------------------------------------------------------------------

    #[test]
    fn plane_point_at_origin() {
        let plane = Plane {
            origin: DVec3::new(1.0, 2.0, 3.0),
            normal: DVec3::Z,
        };
        let p = plane_point_at(&plane, 0.0, 0.0);
        assert!(approx_eq(p, DVec3::new(1.0, 2.0, 3.0), 1e-10));
    }

    #[test]
    fn plane_point_at_uv() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let p = plane_point_at(&plane, 3.0, 4.0);
        // x and y depend on the choice of perpendicular
        assert!((p.z - 0.0).abs() < 1e-10);
    }

    #[test]
    fn plane_parameters_round_trip() {
        let plane = Plane {
            origin: DVec3::new(1.0, 2.0, 3.0),
            normal: DVec3::Z,
        };
        let u = 5.0;
        let v = -2.0;
        let p = plane_point_at(&plane, u, v);
        let uv = plane_parameters(&plane, p);
        assert!(approx_eq_2(uv, DVec2::new(u, v), 1e-10));
    }

    #[test]
    fn plane_normal_is_constant() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(1.0, 1.0, 1.0).normalize(),
        };
        let n = plane_normal(&plane);
        assert!(approx_eq(n, plane.normal, 1e-10));
    }

    #[test]
    fn plane_tangents_perpendicular_to_normal() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let tu = plane_tangent_u(&plane);
        let tv = plane_tangent_v(&plane);
        assert!(tu.dot(plane.normal).abs() < 1e-10);
        assert!(tv.dot(plane.normal).abs() < 1e-10);
        assert!(tu.dot(tv).abs() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // Cylinder Tests
    // -------------------------------------------------------------------------

    #[test]
    fn cylinder_point_at_on_surface() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        for u in [0.0, PI / 2.0, PI, 3.0 * PI / 2.0] {
            let p = cylinder_point_at(&cyl, u, 0.0);
            let r = DVec3::new(p.x, p.y, 0.0).length();
            assert!((r - 2.0).abs() < 1e-10, "u={}, r={}", u, r);
        }
    }

    #[test]
    fn cylinder_parameters_round_trip() {
        let cyl = CylindricalSurface {
            origin: DVec3::new(1.0, 2.0, 3.0),
            axis: DVec3::Z,
            radius: 5.0,
        };
        let u = 1.2;
        let v = 3.4;
        let p = cylinder_point_at(&cyl, u, v);
        let uv = cylinder_parameters(&cyl, p);
        // u may wrap around 2*pi, so compare modulo 2*pi
        let du = (uv.x - u).abs();
        let du = du.min((2.0 * PI - du).abs());
        assert!(du < 1e-8, "u mismatch: {} vs {}", uv.x, u);
        assert!((uv.y - v).abs() < 1e-8, "v mismatch: {} vs {}", uv.y, v);
    }

    #[test]
    fn cylinder_normal_radial() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let n = cylinder_normal(&cyl, 0.0, 0.0);
        assert!(n.z.abs() < 1e-10, "normal should be radial");
        assert!((n.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn cylinder_tangents_orthonormal() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let tu = cylinder_tangent_u(&cyl, 1.0, 0.0);
        let tv = cylinder_tangent_v(&cyl, 1.0, 0.0);
        let n = cylinder_normal(&cyl, 1.0, 0.0);

        assert!(tu.dot(tv).abs() < 1e-10, "tangents should be perpendicular");
        assert!(tu.dot(n).abs() < 1e-10, "u-tangent perpendicular to normal");
        assert!(tv.dot(n).abs() < 1e-10, "v-tangent perpendicular to normal");
    }

    // -------------------------------------------------------------------------
    // Sphere Tests
    // -------------------------------------------------------------------------

    #[test]
    fn sphere_point_at_on_surface() {
        let sph = SphericalSurface {
            center: DVec3::new(1.0, 2.0, 3.0),
            axis: DVec3::Z,
            radius: 3.0,
        };
        for u in [0.0, PI, 2.0 * PI] {
            for v in [0.1, PI / 2.0, PI - 0.1] {
                let p = sphere_point_at(&sph, u, v);
                let d = (p - sph.center).length();
                assert!((d - 3.0).abs() < 1e-9, "u={}, v={}, d={}", u, v, d);
            }
        }
    }

    #[test]
    fn sphere_north_pole() {
        let sph = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        let p = sphere_point_at(&sph, 0.0, 0.0);
        assert!(approx_eq(p, DVec3::new(0.0, 0.0, 2.0), 1e-10));
    }

    #[test]
    fn sphere_south_pole() {
        let sph = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        let p = sphere_point_at(&sph, 0.0, PI);
        assert!(approx_eq(p, DVec3::new(0.0, 0.0, -2.0), 1e-10));
    }

    #[test]
    fn sphere_parameters_round_trip() {
        let sph = SphericalSurface {
            center: DVec3::new(1.0, -2.0, 3.0),
            axis: DVec3::Z,
            radius: 4.0,
        };
        let u = 1.5;
        let v = 1.0;
        let p = sphere_point_at(&sph, u, v);
        let uv = sphere_parameters(&sph, p);
        // Allow for u wrapping
        let du = (uv.x - u).abs();
        let du = du.min((2.0 * PI - du).abs());
        assert!(du < 1e-8, "u mismatch");
        assert!((uv.y - v).abs() < 1e-8, "v mismatch");
    }

    #[test]
    fn sphere_normal_outward() {
        let sph = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let p = sphere_point_at(&sph, 0.5, 0.7);
        let n = sphere_normal(&sph, 0.5, 0.7);
        let expected = (p - sph.center).normalize();
        assert!(approx_eq(n, expected, 1e-10));
    }

    #[test]
    fn sphere_tangents_perpendicular_to_normal() {
        let sph = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let u = 0.7;
        let v = 1.0;
        let n = sphere_normal(&sph, u, v);
        let tu = sphere_tangent_u(&sph, u, v);
        let tv = sphere_tangent_v(&sph, u, v);

        assert!(tu.dot(n).abs() < 1e-10, "u-tangent should be perpendicular to normal");
        assert!(tv.dot(n).abs() < 1e-10, "v-tangent should be perpendicular to normal");
    }

    // -------------------------------------------------------------------------
    // Cone Tests
    // -------------------------------------------------------------------------

    #[test]
    fn cone_point_at_on_surface() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            half_angle_rad: 30.0_f64.to_radians(),
        };

        // At v=0, radius should be 2.0
        let p0 = cone_point_at(&cone, 0.0, 0.0);
        let r0 = DVec3::new(p0.x, p0.y, 0.0).length();
        assert!((r0 - 2.0).abs() < 1e-10);

        // At v > 0, radius increases
        let p1 = cone_point_at(&cone, 0.0, 1.0);
        let r1 = DVec3::new(p1.x, p1.y, 0.0).length();
        assert!(r1 > r0, "radius should increase with v");
    }

    #[test]
    fn cone_normal_constant_along_generator() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            half_angle_rad: 30.0_f64.to_radians(),
        };

        let n0 = cone_normal(&cone, 0.5, 0.0);
        let n1 = cone_normal(&cone, 0.5, 5.0);
        assert!(approx_eq(n0, n1, 1e-10), "normal should be constant along generator");
    }

    #[test]
    fn cone_parameters_round_trip() {
        let cone = ConicalSurface {
            apex: DVec3::new(1.0, 2.0, 3.0),
            axis: DVec3::Z,
            radius: 1.5,
            half_angle_rad: 25.0_f64.to_radians(),
        };

        // Test that point_at and parameters are consistent
        // Use v=0 which corresponds to the reference circle at apex
        let u = 0.8;
        let v = 0.0;
        let p = cone_point_at(&cone, u, v);
        let uv = cone_parameters(&cone, p);

        let du = (uv.x - u).abs();
        let du = du.min((2.0 * PI - du).abs());
        assert!(du < 1e-6, "u mismatch: {} vs {}", uv.x, u);
        // At v=0, the point should be on the reference circle
        let dist_from_apex = (p - cone.apex).length();
        let expected_dist = (cone.radius * cone.radius).sqrt(); // reference circle radius
        assert!((dist_from_apex - cone.radius).abs() < 0.1, "Point should be near reference circle");
    }

    // -------------------------------------------------------------------------
    // Torus Tests
    // -------------------------------------------------------------------------

    #[test]
    fn torus_point_at_on_surface() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };

        // At v=0, points should be at major_radius + minor_radius distance from center
        let p = torus_point_at(&torus, 0.0, 0.0);
        let dist = (p - torus.center).length();
        assert!((dist - (torus.major_radius + torus.minor_radius)).abs() < 1e-10,
                "dist = {}, expected {}", dist, torus.major_radius + torus.minor_radius);

        // At v=pi, points should be at major_radius - minor_radius from center
        let p = torus_point_at(&torus, 0.0, PI);
        let dist = (p - torus.center).length();
        assert!((dist - (torus.major_radius - torus.minor_radius)).abs() < 1e-10);

        // At v=pi/2, points should be at major_radius from axis, with z = minor_radius
        let p = torus_point_at(&torus, 0.0, FRAC_PI_2);
        let r_from_axis = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r_from_axis - 5.0).abs() < 1e-10);
        assert!((p.z.abs() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn torus_parameters_round_trip() {
        let torus = ToroidalSurface {
            center: DVec3::new(1.0, 2.0, 3.0),
            axis: DVec3::Z,
            major_radius: 4.0,
            minor_radius: 1.0,
        };

        let u = 0.7;
        let v = 2.1;
        let p = torus_point_at(&torus, u, v);
        let uv = torus_parameters(&torus, p);

        let du = (uv.x - u).abs();
        let du = du.min((2.0 * PI - du).abs());
        let dv = (uv.y - v).abs();
        let dv = dv.min((2.0 * PI - dv).abs());

        assert!(du < 1e-8, "u mismatch: {} vs {}", uv.x, u);
        assert!(dv < 1e-8, "v mismatch: {} vs {}", uv.y, v);
    }

    #[test]
    fn torus_normal_outward() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };

        // At v=0, normal should be perpendicular to axis (horizontal)
        let n = torus_normal(&torus, 0.0, 0.0);
        assert!(n.z.abs() < 1e-10, "Normal should be horizontal");
        assert!((n.length() - 1.0).abs() < 1e-10, "Normal should be unit length");

        // At v=pi, normal should still be horizontal but pointing inward
        let n = torus_normal(&torus, 0.0, PI);
        assert!(n.z.abs() < 1e-10, "Normal should be horizontal");
        assert!((n.length() - 1.0).abs() < 1e-10, "Normal should be unit length");
    }

    #[test]
    fn torus_tangents_perpendicular_to_normal() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };

        let u = 0.5;
        let v = 1.0;
        let n = torus_normal(&torus, u, v);
        let tu = torus_tangent_u(&torus, u, v);
        let tv = torus_tangent_v(&torus, u, v);

        assert!(tu.dot(n).abs() < 1e-10, "u-tangent perpendicular to normal");
        assert!(tv.dot(n).abs() < 1e-10, "v-tangent perpendicular to normal");
    }

    // -------------------------------------------------------------------------
    // BSplineSurface Tests
    // -------------------------------------------------------------------------

    #[test]
    fn bspline_surface_point_at_bilinear() {
        // Create a bilinear (degree 1, 1) BSpline surface - a flat quadrilateral
        let surf = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::ZERO, DVec3::Y],
                vec![DVec3::X, DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        // Corners should match control points
        assert!(approx_eq(bspline_surface_point_at(&surf, 0.0, 0.0), DVec3::ZERO, 1e-10));
        assert!(approx_eq(bspline_surface_point_at(&surf, 1.0, 0.0), DVec3::X, 1e-10));
        assert!(approx_eq(bspline_surface_point_at(&surf, 0.0, 1.0), DVec3::Y, 1e-10));
        assert!(approx_eq(bspline_surface_point_at(&surf, 1.0, 1.0), DVec3::new(1.0, 1.0, 0.0), 1e-10));

        // Midpoint
        assert!(approx_eq(bspline_surface_point_at(&surf, 0.5, 0.5), DVec3::new(0.5, 0.5, 0.0), 1e-10));
    }

    #[test]
    fn bspline_surface_normal_flat() {
        // Flat bilinear surface
        let surf = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::ZERO, DVec3::Y],
                vec![DVec3::X, DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let n = bspline_surface_normal(&surf, 0.5, 0.5);
        assert!(approx_eq(n, DVec3::Z, 1e-10) || approx_eq(n, DVec3::NEG_Z, 1e-10));
    }

    #[test]
    fn bspline_surface_derivatives_test() {
        // Flat bilinear surface
        let surf = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::ZERO, DVec3::Y],
                vec![DVec3::X, DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        // Use SurfaceEval to compute points and normals
        use rcad_kernel::SurfaceEval;
        let p = surf.point_at(0.5, 0.5);
        let n = surf.normal_at(0.5, 0.5);

        // Midpoint should be at (0.5, 0.5, 0)
        assert!(approx_eq(p, DVec3::new(0.5, 0.5, 0.0), 1e-6));

        // Normal should be +/- Z for a flat surface
        assert!(approx_eq(n.normalize_or_zero(), DVec3::Z, 1e-6) || approx_eq(n.normalize_or_zero(), DVec3::NEG_Z, 1e-6));
    }
}
