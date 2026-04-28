//! Differential geometry: principal curvatures, Gaussian curvature, and mean
//! curvature for all `Surface3` analytic types.
//!
//! Analytic solutions are used for Plane, Cylinder, Sphere, Cone, and Torus.
//! BSpline surfaces use a numerical finite-difference approximation of the
//! shape operator (Weingarten map).
//!
//! **Parameter conventions** (same as `SurfaceEval`):
//! - `Cylinder`: u = azimuth [0, 2π], v = height along axis
//! - `Sphere`:   u = longitude [0, 2π], v = colatitude [0, π]
//! - `Cone`:     u = azimuth [0, 2π], v = slant distance from apex (≥ 0)
//! - `Torus`:    u = major angle [0, 2π], v = minor angle [0, 2π]
//! - `BSpline`:  u, v clamped to `default_domain()`

use crate::geom::Surface3;

// ── Public API ────────────────────────────────────────────────────────────────

/// Returns the two principal curvatures `(k1, k2)` at surface parameter `(u, v)`.
///
/// The sign convention follows the outward normal: positive curvature means the
/// surface bends away from the outward normal direction (convex).
///
/// For a flat surface both values are 0.  For a sphere of radius `r` both are
/// `1/r`.  For a cylinder of radius `r` one is `1/r` and the other is 0.
pub fn principal_curvatures(surface: &Surface3, u: f64, v: f64) -> (f64, f64) {
    match surface {
        Surface3::Plane(_) => (0.0, 0.0),

        Surface3::Cylinder(cyl) => {
            // k_circumferential = 1/r, k_axial = 0
            (1.0 / cyl.radius, 0.0)
        }

        Surface3::Sphere(sph) => {
            let k = 1.0 / sph.radius;
            (k, k)
        }

        Surface3::Cone(cone) => {
            // v = slant distance from apex; r_at = v * sin(half_angle)
            let alpha = cone.half_angle_rad;
            let r_at = (v * alpha.sin()).abs();
            if r_at < 1e-12 {
                // At the apex: the curvature is singular
                (f64::INFINITY, 0.0)
            } else {
                // k_circumferential = sin(α)/r_at, k_axial = 0
                (alpha.sin() / r_at, 0.0)
            }
        }

        Surface3::Torus(tor) => {
            // v = minor angle (around the tube)
            // k1 (along the tube) = 1/r
            // k2 (around the major circle) = cos(v) / (R + r·cos(v))
            let r = tor.minor_radius;
            let big_r = tor.major_radius;
            let cos_v = v.cos();
            let k1 = 1.0 / r;
            let k2 = cos_v / (big_r + r * cos_v);
            // Return larger magnitude first for consistency
            if k1.abs() >= k2.abs() {
                (k1, k2)
            } else {
                (k2, k1)
            }
        }

        Surface3::BSpline(_) => numerical_curvatures(surface, u, v),
        Surface3::Ellipsoid(_) => numerical_curvatures(surface, u, v),
        Surface3::Helicoid(_) => numerical_curvatures(surface, u, v),
        Surface3::Pipe(_) => numerical_curvatures(surface, u, v),
        Surface3::TriBezier(_) => numerical_curvatures(surface, u, v),
        Surface3::LinearExtrusion(_)
        | Surface3::Revolution(_)
        | Surface3::Ruled(_)
        | Surface3::Coons(_)
        | Surface3::Bezier(_)
        | Surface3::Offset(_)
        | Surface3::Trimmed(_) => numerical_curvatures(surface, u, v),
    }
}

/// Gaussian curvature K = k1 · k2.
///
/// K > 0: elliptic point (sphere-like). K = 0: parabolic (cylinder-like or flat).
/// K < 0: hyperbolic (saddle). K = NaN at cone apex (singular point).
pub fn gaussian_curvature(surface: &Surface3, u: f64, v: f64) -> f64 {
    let (k1, k2) = principal_curvatures(surface, u, v);
    k1 * k2
}

/// Mean curvature H = (k1 + k2) / 2.
///
/// H = 0 is the minimal-surface condition. H = NaN at cone apex.
pub fn mean_curvature(surface: &Surface3, u: f64, v: f64) -> f64 {
    let (k1, k2) = principal_curvatures(surface, u, v);
    (k1 + k2) * 0.5
}

/// Returns the radius of the osculating sphere at the given point.
///
/// This is the reciprocal of the maximum principal curvature magnitude:
/// `r = 1 / max(|k1|, |k2|)`. Returns `f64::INFINITY` if both curvatures are zero
/// (flat point), and 0.0 if either curvature is infinite (singular point).
pub fn osculating_radius(surface: &Surface3, u: f64, v: f64) -> f64 {
    let (k1, k2) = principal_curvatures(surface, u, v);
    let k_max = k1.abs().max(k2.abs());
    if k_max.is_infinite() {
        0.0
    } else if k_max < 1e-12 {
        f64::INFINITY
    } else {
        1.0 / k_max
    }
}

// ── Numerical fallback for BSpline ───────────────────────────────────────────

/// Finite-difference approximation of principal curvatures for BSpline surfaces.
///
/// Uses the classical formula via the first and second fundamental forms:
/// ```text
/// E = Pu·Pu,  F = Pu·Pv,  G = Pv·Pv   (first fundamental form)
/// e = Puu·N,  f = Puv·N,  g = Pvv·N   (second fundamental form)
/// K = (e·g − f²) / (E·G − F²)
/// H = (E·g − 2·F·f + G·e) / (2·(E·G − F²))
/// k1,k2 = H ± sqrt(max(0, H² − K))
/// ```
fn numerical_curvatures(surface: &Surface3, u: f64, v: f64) -> (f64, f64) {
    use crate::geom::SurfaceEval;
    const EPS: f64 = 1e-5;

    // Clamp (u,v) away from the domain boundary so all 9-point stencil evaluations
    // remain within the valid parameter range.
    let [u0, u1, v0, v1] = surface.default_domain();
    let u = u.clamp(u0 + EPS * 3.0, u1 - EPS * 3.0);
    let v = v.clamp(v0 + EPS * 3.0, v1 - EPS * 3.0);

    let p = surface.point_at(u, v);
    let p_up = surface.point_at(u + EPS, v);
    let p_um = surface.point_at(u - EPS, v);
    let p_vp = surface.point_at(u, v + EPS);
    let p_vm = surface.point_at(u, v - EPS);
    let p_pp = surface.point_at(u + EPS, v + EPS);
    let p_pm = surface.point_at(u + EPS, v - EPS);
    let p_mp = surface.point_at(u - EPS, v + EPS);
    // p_mm not needed for 4-point cross-derivative

    // First-order partials
    let pu = (p_up - p_um) / (2.0 * EPS);
    let pv = (p_vp - p_vm) / (2.0 * EPS);

    // Second-order partials
    let eps2 = EPS * EPS;
    let puu = (p_up - 2.0 * p + p_um) / eps2;
    let pvv = (p_vp - 2.0 * p + p_vm) / eps2;
    let puv = (p_pp - p_pm - p_mp + surface.point_at(u - EPS, v - EPS)) / (4.0 * eps2);

    // Unit normal
    let n_raw = pu.cross(pv);
    let n_len = n_raw.length();
    if n_len < 1e-12 {
        return (0.0, 0.0); // degenerate
    }
    let n = n_raw / n_len;

    // First fundamental form coefficients
    let ee = pu.dot(pu);
    let ff = pu.dot(pv);
    let gg = pv.dot(pv);
    let denom = ee * gg - ff * ff;
    if denom.abs() < 1e-20 {
        return (0.0, 0.0); // degenerate metric
    }

    // Second fundamental form coefficients
    let e = puu.dot(n);
    let f = puv.dot(n);
    let g = pvv.dot(n);

    let big_k = (e * g - f * f) / denom;
    let big_h = (ee * g - 2.0 * ff * f + gg * e) / (2.0 * denom);

    let discriminant = (big_h * big_h - big_k).max(0.0);
    let sqrt_d = discriminant.sqrt();
    (big_h + sqrt_d, big_h - sqrt_d)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{
        BSplineSurface, ConicalSurface, CylindricalSurface, Plane, SphericalSurface,
        ToroidalSurface,
    };
    use glam::DVec3;

    const TOL: f64 = 1e-9;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn plane_has_zero_curvature() {
        let s = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let (k1, k2) = principal_curvatures(&s, 0.0, 0.0);
        assert_eq!(k1, 0.0);
        assert_eq!(k2, 0.0);
        assert_eq!(gaussian_curvature(&s, 0.0, 0.0), 0.0);
        assert_eq!(mean_curvature(&s, 0.0, 0.0), 0.0);
    }

    #[test]
    fn cylinder_radius1_curvatures() {
        let r = 1.0_f64;
        let s = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: r,
        });
        let (k1, k2) = principal_curvatures(&s, 0.0, 0.0);
        assert!(
            approx_eq(k1.max(k2), 1.0 / r, TOL),
            "max k = {}",
            k1.max(k2)
        );
        assert!(approx_eq(k1.min(k2), 0.0, TOL), "min k = {}", k1.min(k2));
        assert!(approx_eq(gaussian_curvature(&s, 0.0, 0.0), 0.0, TOL));
        assert!(approx_eq(mean_curvature(&s, 0.0, 0.0), 0.5 / r, TOL));
    }

    #[test]
    fn sphere_radius2_curvatures() {
        let r = 2.0_f64;
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: r,
        });
        let (k1, k2) = principal_curvatures(&s, 0.0, std::f64::consts::FRAC_PI_2);
        assert!(approx_eq(k1, 1.0 / r, TOL));
        assert!(approx_eq(k2, 1.0 / r, TOL));
        assert!(approx_eq(
            gaussian_curvature(&s, 0.0, 0.0),
            1.0 / (r * r),
            TOL
        ));
        assert!(approx_eq(mean_curvature(&s, 0.0, 0.0), 1.0 / r, TOL));
    }

    #[test]
    fn cone_at_nonzero_slant() {
        let half_angle = std::f64::consts::FRAC_PI_4; // 45°
        let s = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 0.0,
            half_angle_rad: half_angle,
        });
        let v = 1.0_f64; // slant distance = 1
        let r_at = v * half_angle.sin();
        let (k1, k2) = principal_curvatures(&s, 0.0, v);
        assert!(approx_eq(k1, half_angle.sin() / r_at, TOL));
        assert!(approx_eq(k2, 0.0, TOL));
        assert!(approx_eq(gaussian_curvature(&s, 0.0, v), 0.0, TOL));
    }

    #[test]
    fn torus_outer_equator_curvatures() {
        let big_r = 2.0_f64;
        let r = 0.5_f64;
        let s = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: big_r,
            minor_radius: r,
        });
        // At v=0 (outer equator): k_tube = 1/r = 2.0, k_major = 1/(R+r) = 1/2.5 = 0.4
        let (k1, k2) = principal_curvatures(&s, 0.0, 0.0);
        let expected_k_tube = 1.0 / r;
        let expected_k_major = 1.0 / (big_r + r);
        assert!(
            approx_eq(k1.max(k2), expected_k_tube, TOL),
            "k_tube mismatch: {k1} {k2}"
        );
        assert!(
            approx_eq(k1.min(k2), expected_k_major, TOL),
            "k_major mismatch: {k1} {k2}"
        );
    }

    #[test]
    fn bspline_flat_bilinear_near_zero_curvature() {
        // A bilinear (degree 1×1) patch over a flat unit square — curvature ≈ 0
        let s = Surface3::BSpline(BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        });
        let (k1, k2) = principal_curvatures(&s, 0.5, 0.5);
        assert!(k1.abs() < 1e-4, "k1 = {k1}");
        assert!(k2.abs() < 1e-4, "k2 = {k2}");
        assert!(gaussian_curvature(&s, 0.5, 0.5).abs() < 1e-6);
        assert!(mean_curvature(&s, 0.5, 0.5).abs() < 1e-4);
    }

    // ============================================================================
    // OCCT TKGeomBase Alignment Tests - Curvature Edge Cases
    // ============================================================================

    #[test]
    fn sphere_curvature_isotropic() {
        // All points on sphere have same principal curvatures
        let r = 5.0_f64;
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: r,
        });
        let (k1, k2) = principal_curvatures(&s, 0.5, 0.5);
        assert!(approx_eq(k1, 1.0 / r, TOL), "k1 should be 1/r");
        assert!(approx_eq(k2, 1.0 / r, TOL), "k2 should be 1/r");
        assert!(approx_eq(gaussian_curvature(&s, 0.5, 0.5), 1.0 / (r * r), TOL));
    }

    #[test]
    fn cylinder_curvature_developable() {
        // Cylinder has one zero curvature (axial), one non-zero (circumferential)
        let r = 3.0_f64;
        let s = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: r,
        });
        let (k1, k2) = principal_curvatures(&s, 0.0, 0.0);
        // One should be 1/r, the other should be 0
        let k_max = k1.max(k2);
        let k_min = k1.min(k2);
        assert!(approx_eq(k_max, 1.0 / r, TOL), "k_max should be 1/r");
        assert!(approx_eq(k_min, 0.0, TOL), "k_min should be 0 (developable)");
        // Gaussian curvature should be 0 (developable surface)
        assert!(approx_eq(gaussian_curvature(&s, 0.0, 0.0), 0.0, TOL));
    }

    #[test]
    fn plane_curvature_zero() {
        let s = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let (k1, k2) = principal_curvatures(&s, 0.0, 0.0);
        assert!(approx_eq(k1, 0.0, TOL), "plane k1 should be 0");
        assert!(approx_eq(k2, 0.0, TOL), "plane k2 should be 0");
        assert!(approx_eq(gaussian_curvature(&s, 0.0, 0.0), 0.0, TOL));
        assert!(approx_eq(mean_curvature(&s, 0.0, 0.0), 0.0, TOL));
    }

    #[test]
    fn torus_inner_equator_negative_curvature() {
        // At inner equator (v = π), principal curvatures have opposite signs
        let big_r = 4.0_f64;
        let r = 1.0_f64;
        let s = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: big_r,
            minor_radius: r,
        });
        // At v = π (inner equator)
        let (k1, k2) = principal_curvatures(&s, 0.0, std::f64::consts::PI);
        // One curvature should be positive (tube), one negative (inner side)
        let k_tube = 1.0 / r;
        let k_major = -1.0 / (big_r - r); // negative because inner
        let k_max = k1.max(k2);
        let k_min = k1.min(k2);
        assert!(approx_eq(k_max, k_tube, 0.1), "inner equator k_max");
        assert!(k_min < 0.0, "inner equator k_min should be negative");
    }

    #[test]
    fn cone_curvature_singular_at_apex() {
        // At apex, curvature is undefined/very large
        let s = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: std::f64::consts::FRAC_PI_6,
        });
        // Near apex (small v)
        let (k1, k2) = principal_curvatures(&s, 0.0, 0.001);
        // Curvature should be very large near apex
        // But should not crash
        assert!(k1.is_finite() && k2.is_finite());
    }

    #[test]
    fn osculating_radius_sphere() {
        let r = 2.0_f64;
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: r,
        });
        let osr = osculating_radius(&s, 0.5, 0.5);
        assert!(approx_eq(osr, r, TOL), "osculating radius should equal sphere radius");
    }

    #[test]
    fn mean_curvature_formula() {
        // Mean curvature H = (k1 + k2) / 2
        let r = 1.0_f64;
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: r,
        });
        let (k1, k2) = principal_curvatures(&s, 0.0, 0.0);
        let h = mean_curvature(&s, 0.0, 0.0);
        assert!(approx_eq(h, (k1 + k2) / 2.0, TOL), "mean curvature formula");
    }

    #[test]
    fn gaussian_curvature_formula() {
        // Gaussian curvature K = k1 * k2
        let r = 2.0_f64;
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: r,
        });
        let (k1, k2) = principal_curvatures(&s, 0.0, 0.0);
        let k = gaussian_curvature(&s, 0.0, 0.0);
        assert!(approx_eq(k, k1 * k2, TOL), "gaussian curvature formula");
    }
}
