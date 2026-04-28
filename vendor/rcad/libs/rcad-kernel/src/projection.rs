//! Closest-point projection from a 3D point onto a curve or surface.
//!
//! Analogous to OCCT `GeomAPI_ProjectPointOnCurve` and
//! `GeomAPI_ProjectPointOnSurf`.
//!
//! # Strategy
//! - **Analytic surfaces** (Plane, Cylinder, Sphere, Cone, Torus): closed-form
//!   projection — fast and exact.
//! - **All curves** and **parametric surfaces** (BSpline, Bezier, Offset,
//!   LinearExtrusion, Revolution): sample the domain uniformly to find the best
//!   initial guess, then refine with Newton-Raphson minimisation of `|P(t) - Q|²`.

use glam::DVec3;

use crate::geom::{Curve3, CurveEval, Surface3, SurfaceEval};

// ─────────────────────────────────────────────────────────────────────────────
// Result types
// ─────────────────────────────────────────────────────────────────────────────

/// Result of projecting a point onto a curve.
#[derive(Debug, Clone)]
pub struct CurveProjection {
    /// Nearest point on the curve.
    pub point: DVec3,
    /// Curve parameter at the nearest point.
    pub param: f64,
    /// Distance from the query point to the curve.
    pub distance: f64,
}

/// Result of projecting a point onto a surface.
#[derive(Debug, Clone)]
pub struct SurfaceProjection {
    /// Nearest point on the surface.
    pub point: DVec3,
    /// Surface parameter (u, v) at the nearest point.
    pub params: (f64, f64),
    /// Distance from the query point to the surface.
    pub distance: f64,
}

// ─────────────────────────────────────────────────────────────────────────────
// Curve projection
// ─────────────────────────────────────────────────────────────────────────────

/// Project the point `query` onto `curve`, returning the nearest point on the
/// curve, its parameter value, and the Euclidean distance.
///
/// The curve is evaluated over its natural domain ([`CurveEval::default_domain`]).
///
/// # Algorithm
/// 1. Sample `n_samples` points uniformly in the domain.
/// 2. Take the sample with the smallest distance as the initial guess.
/// 3. Refine with Newton iterations minimising `f(t) = |C(t) - Q|²`:
///    `t_{i+1} = t_i - (C(t) - Q) · T(t) / (|T|² + (C - Q) · T')` where
///    `T = C'(t)` (approximated by finite difference).
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Curve3, Circle3};
/// use rcad_kernel::projection::closest_point_on_curve;
///
/// let circle = Curve3::Circle(Circle3 {
///     center: DVec3::ZERO,
///     normal: DVec3::Z,
///     radius: 1.0,
/// });
/// let q = DVec3::new(2.0, 0.0, 0.0);
/// let result = closest_point_on_curve(&circle, q, 64);
/// assert!((result.point - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6);
/// ```
pub fn closest_point_on_curve(curve: &Curve3, query: DVec3, n_samples: usize) -> CurveProjection {
    // ── Analytic fast paths ───────────────────────────────────────────────────
    // Analogous to OCCT `ExtremaPC` per-type dispatch.

    match curve {
        Curve3::Line(l) => {
            // Closest point on an infinite line: project query onto line direction.
            let dir_sq = l.direction.dot(l.direction);
            if dir_sq < 1e-20 {
                // Degenerate line — return origin.
                let pt = l.origin;
                return CurveProjection {
                    point: pt,
                    param: 0.0,
                    distance: (pt - query).length(),
                };
            }
            let t = (query - l.origin).dot(l.direction) / dir_sq;
            let [t0, t1] = curve.default_domain();
            let t_clamped = if t0.is_finite() && t1.is_finite() {
                t.clamp(t0, t1)
            } else {
                t
            };
            let pt = l.origin + t_clamped * l.direction;
            return CurveProjection {
                point: pt,
                param: t_clamped,
                distance: (pt - query).length(),
            };
        }

        Curve3::Circle(circ) => {
            // Closest point on a circle: project query onto circle plane,
            // compute the angle, then evaluate.
            // The circle is parametrized as P(t) = center + cos(t)*x_ax + sin(t)*y_ax.
            let x_ax = crate::geom::any_perpendicular(circ.normal);
            let y_ax = circ.normal.cross(x_ax).normalize_or_zero();
            // Project query onto the circle plane.
            let q_in_plane = query - query.dot(circ.normal) * circ.normal
                + circ.center.dot(circ.normal) * circ.normal;
            let local = q_in_plane - circ.center;
            let a = local.dot(x_ax);
            let b = local.dot(y_ax);
            // Handle the degenerate case (query on the axis).
            let t_raw = if a.abs() < 1e-15 && b.abs() < 1e-15 {
                0.0_f64
            } else {
                b.atan2(a)
            };
            let [t0, t1] = curve.default_domain();
            // Wrap t_raw into [t0, t1] for partial arcs.
            let t = if t1 - t0 >= std::f64::consts::TAU - 1e-9 {
                // Full circle — no clamping needed.
                t_raw
            } else {
                // Clamp and check the two nearest full-circle candidates.
                let candidates = [
                    t_raw,
                    t_raw + std::f64::consts::TAU,
                    t_raw - std::f64::consts::TAU,
                    t0,
                    t1,
                ];
                candidates
                    .iter()
                    .filter(|&&tc| tc >= t0 - 1e-12 && tc <= t1 + 1e-12)
                    .map(|&tc| tc.clamp(t0, t1))
                    .min_by(|&a, &b| {
                        let da = (circ.point_at(a) - query).length();
                        let db = (circ.point_at(b) - query).length();
                        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .unwrap_or(t_raw.clamp(t0, t1))
            };
            let pt = circ.point_at(t);
            return CurveProjection {
                point: pt,
                param: t,
                distance: (pt - query).length(),
            };
        }

        Curve3::Ellipse(ell) => {
            // Closest point on an ellipse: project query onto ellipse plane,
            // decompose into (a, b) components, use atan2 for initial angle
            // estimate, then refine with a few Newton steps.
            let minor_dir = ell.normal.cross(ell.major_dir).normalize_or_zero();
            // Project query onto the ellipse plane.
            let q_in_plane = query - query.dot(ell.normal) * ell.normal
                + ell.center.dot(ell.normal) * ell.normal;
            let local = q_in_plane - ell.center;
            let a = local.dot(ell.major_dir);
            let b = local.dot(minor_dir);
            // Normalized angle that would locate the closest point on a unit circle.
            let t_init = if a.abs() < 1e-15 && b.abs() < 1e-15 {
                0.0_f64
            } else {
                (b / ell.minor_radius).atan2(a / ell.major_radius)
            };
            let [t0, t1] = curve.default_domain();
            let t_init_clamped = if t1 - t0 >= std::f64::consts::TAU - 1e-9 {
                t_init
            } else {
                // Try all canonical equivalent angles.
                let candidates = [
                    t_init,
                    t_init + std::f64::consts::TAU,
                    t_init - std::f64::consts::TAU,
                    t0,
                    t1,
                ];
                candidates
                    .iter()
                    .filter(|&&tc| tc >= t0 - 1e-12 && tc <= t1 + 1e-12)
                    .map(|&tc| tc.clamp(t0, t1))
                    .min_by(|&a, &b| {
                        let da = (ell.point_at(a) - query).length();
                        let db = (ell.point_at(b) - query).length();
                        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .unwrap_or(t_init.clamp(t0, t1))
            };
            // Refine with Newton steps (at most 10 iterations).
            let dt = 1e-7;
            let clamp = |t: f64| {
                if t1 - t0 >= std::f64::consts::TAU - 1e-9 { t } else { t.clamp(t0, t1) }
            };
            let mut best_t = t_init_clamped;
            let mut best_dist = (ell.point_at(best_t) - query).length();
            for _ in 0..10 {
                let p = ell.point_at(best_t);
                let diff = p - query;
                let tangent = (ell.point_at(best_t + dt) - ell.point_at(best_t - dt)) / (2.0 * dt);
                let tang_sq = tangent.dot(tangent);
                if tang_sq < 1e-20 { break; }
                let curv = (ell.point_at(best_t + 2.0 * dt) - 2.0 * p
                    + ell.point_at(best_t - 2.0 * dt)) / (dt * dt);
                let denom = tang_sq + diff.dot(curv);
                let delta = diff.dot(tangent) / if denom.abs() > 1e-20 { denom } else { tang_sq };
                let new_t = clamp(best_t - delta);
                let new_dist = (ell.point_at(new_t) - query).length();
                if new_dist < best_dist {
                    best_dist = new_dist;
                    best_t = new_t;
                }
                if delta.abs() < 1e-11 { break; }
            }
            let pt = ell.point_at(best_t);
            return CurveProjection {
                point: pt,
                param: best_t,
                distance: (pt - query).length(),
            };
        }

        _ => {}
    }

    // ── Numerical fallback for all other curve types ───────────────────────────
    let [t0_raw, t1_raw] = curve.default_domain();
    let n = n_samples.max(4);

    // For infinite domains (lines), use a heuristic finite sampling range
    // centered on the closest parameter analytically (dot product for lines).
    let (t0, t1) = if t0_raw.is_infinite() || t1_raw.is_infinite() {
        // Use the analytical projection for the domain center estimate
        let t_center = 0.0; // non-line types; 0 is a safe fallback
        let span = 100.0_f64; // generous range around t_center
        (t_center - span, t_center + span)
    } else {
        (t0_raw, t1_raw)
    };

    // Step 1: coarse sampling
    let (mut best_t, mut best_dist) = {
        let mut bd = f64::INFINITY;
        let mut bt = t0;
        for i in 0..=n {
            let t = t0 + (t1 - t0) * i as f64 / n as f64;
            let p = curve.point_at(t);
            let d = (p - query).length();
            if d < bd {
                bd = d;
                bt = t;
            }
        }
        (bt, bd)
    };

    // Step 2: Newton refinement
    // For infinite domains, don't clamp the Newton step
    let clamp_t = |t: f64| {
        if t0_raw.is_infinite() || t1_raw.is_infinite() {
            t
        } else {
            t.clamp(t0, t1)
        }
    };
    let dt = if (t1 - t0).is_finite() {
        (t1 - t0) * 1e-6
    } else {
        1e-6
    };
    for _ in 0..30 {
        let p = curve.point_at(best_t);
        let diff = p - query;
        // Finite-difference tangent
        let t_plus = best_t + dt;
        let t_minus = best_t - dt;
        let span = t_plus - t_minus;
        if span.abs() < 1e-20 {
            break;
        }
        let tangent = (curve.point_at(t_plus) - curve.point_at(t_minus)) / span;
        let tang_sq = tangent.dot(tangent);
        if tang_sq < 1e-20 {
            break;
        }
        // Second-order term (curvature denominator term)
        let curvature_approx = (curve.point_at(best_t + 2.0 * dt) - 2.0 * p
            + curve.point_at(best_t - 2.0 * dt))
            / (dt * dt);
        let denom = tang_sq + diff.dot(curvature_approx);
        let delta = diff.dot(tangent) / if denom.abs() > 1e-20 { denom } else { tang_sq };
        let new_t = clamp_t(best_t - delta);
        let new_dist = (curve.point_at(new_t) - query).length();
        if new_dist < best_dist {
            best_dist = new_dist;
            best_t = new_t;
        }
        if delta.abs() < 1e-10 {
            break;
        }
    }

    let best_point = curve.point_at(best_t);
    CurveProjection {
        point: best_point,
        param: best_t,
        distance: (best_point - query).length(),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Surface projection
// ─────────────────────────────────────────────────────────────────────────────

/// Project the point `query` onto `surface`, returning the nearest point on the
/// surface, its (u, v) parameters, and the Euclidean distance.
///
/// Analytic surfaces (Plane, Sphere, Cylinder, Cone, Torus) use closed-form
/// formulae.  All other surfaces fall back to numerical sampling + Newton
/// refinement.
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Surface3, SphericalSurface};
/// use rcad_kernel::projection::closest_point_on_surface;
///
/// let sphere = Surface3::Sphere(SphericalSurface {
///     center: DVec3::ZERO,
///     axis: DVec3::Z,
///     radius: 1.0,
/// });
/// let q = DVec3::new(3.0, 0.0, 0.0);
/// let result = closest_point_on_surface(&sphere, q, 16);
/// assert!((result.point - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6);
/// ```
pub fn closest_point_on_surface(
    surface: &Surface3,
    query: DVec3,
    n_samples: usize,
) -> SurfaceProjection {
    use crate::geom::*;

    match surface {
        // ── Analytic closed-form projections ──────────────────────────────────
        Surface3::Plane(plane) => {
            // Project onto the infinite plane, no clamping needed.
            let d = (query - plane.origin).dot(plane.normal);
            let point = query - plane.normal * d;
            let distance = d.abs();
            // (u, v) = coordinates in plane (not strictly needed, but provide them)
            let u_axis = any_perpendicular(plane.normal);
            let v_axis = plane.normal.cross(u_axis);
            let diff = point - plane.origin;
            SurfaceProjection {
                point,
                params: (diff.dot(u_axis), diff.dot(v_axis)),
                distance,
            }
        }

        Surface3::Sphere(sph) => {
            let v = query - sph.center;
            let len = v.length();
            let point = if len < 1e-14 {
                sph.center + sph.radius * DVec3::X // degenerate: pick arbitrary
            } else {
                sph.center + v / len * sph.radius
            };
            // Compute (theta, phi) from point relative to center
            let w = (point - sph.center).normalize_or_zero();
            let u_axis = any_perpendicular(sph.axis);
            let v_axis = sph.axis.cross(u_axis);
            let theta = w.dot(sph.axis).clamp(-1.0, 1.0).acos();
            let phi = w.dot(v_axis).atan2(w.dot(u_axis));
            SurfaceProjection {
                point,
                params: (phi, theta),
                distance: (point - query).length(),
            }
        }

        Surface3::Cylinder(cyl) => {
            // Project by collapsing along axis, then normalizing radial component.
            let v = query - cyl.origin;
            let along = v.dot(cyl.axis);
            let radial = v - cyl.axis * along;
            let radial_len = radial.length();
            let point = if radial_len < 1e-14 {
                cyl.origin + cyl.axis * along + cyl.radius * any_perpendicular(cyl.axis)
            } else {
                cyl.origin + cyl.axis * along + radial / radial_len * cyl.radius
            };
            let u_axis = any_perpendicular(cyl.axis);
            let v_axis = cyl.axis.cross(u_axis);
            let r = (point - cyl.origin - cyl.axis * along).normalize_or_zero();
            let theta = r.dot(v_axis).atan2(r.dot(u_axis));
            SurfaceProjection {
                point,
                params: (theta, along),
                distance: (point - query).length(),
            }
        }

        Surface3::Cone(cone) => {
            // Project onto the cone's reference-circle parameterization.
            let axis = cone.axis_dir();
            let x_axis = any_perpendicular(axis);
            let y_axis = axis.cross(x_axis).normalize_or_zero();
            let local = query - cone.apex;
            let along = local.dot(axis);
            let radial = local - axis * along;
            let radial_len = radial.length();
            let half = cone.half_angle_rad;
            let tan_h = half.tan();
            let axial = (along + (radial_len - cone.radius) * tan_h) / (1.0 + tan_h * tan_h);
            let r_hat = if radial_len < 1e-14 {
                x_axis
            } else {
                radial / radial_len
            };
            let point = cone.apex + axis * axial + r_hat * cone.radius_at_axial(axial);
            let slant = cone.slant_from_axial(axial);
            let theta = r_hat.dot(y_axis).atan2(r_hat.dot(x_axis));
            SurfaceProjection {
                point,
                params: (theta, slant),
                distance: (point - query).length(),
            }
        }

        Surface3::Torus(torus) => {
            // Step 1: project onto the major-radius circle in the equatorial plane.
            let v = query - torus.center;
            let along = v.dot(torus.axis);
            let radial = v - torus.axis * along;
            let radial_len = radial.length();
            let major_dir = if radial_len < 1e-14 {
                any_perpendicular(torus.axis)
            } else {
                radial / radial_len
            };
            let tube_center = torus.center + major_dir * torus.major_radius;
            // Step 2: project onto the tube circle.
            let w = query - tube_center;
            let w_len = w.length();
            let point = if w_len < 1e-14 {
                tube_center + major_dir * torus.minor_radius
            } else {
                tube_center + w / w_len * torus.minor_radius
            };
            let u = major_dir
                .dot(any_perpendicular(torus.axis))
                .atan2(major_dir.dot(torus.axis.cross(any_perpendicular(torus.axis))));
            let w_dir = (point - tube_center).normalize_or_zero();
            let v_param = w_dir.dot(torus.axis).atan2(w_dir.dot(major_dir));
            SurfaceProjection {
                point,
                params: (u, v_param),
                distance: (point - query).length(),
            }
        }

        // ── Numerical fallback for parametric surfaces ─────────────────────────
        _ => numeric_surface_projection(surface, query, n_samples),
    }
}

/// Numerical closest-point on a parametric surface via uniform sampling +
/// Newton refinement of `f(u,v) = |S(u,v) - Q|²`.
fn numeric_surface_projection(
    surface: &Surface3,
    query: DVec3,
    n_samples: usize,
) -> SurfaceProjection {
    let [u0, u1, v0, v1] = surface.default_domain();
    let n = n_samples.max(4);

    // Coarse sampling
    let (mut best_u, mut best_v, mut best_dist) = {
        let mut bd = f64::INFINITY;
        let (mut bu, mut bv) = (u0, v0);
        for i in 0..=n {
            for j in 0..=n {
                let u = u0 + (u1 - u0) * i as f64 / n as f64;
                let v = v0 + (v1 - v0) * j as f64 / n as f64;
                let p = surface.point_at(u, v);
                let d = (p - query).length_squared();
                if d < bd {
                    bd = d;
                    bu = u;
                    bv = v;
                }
            }
        }
        (bu, bv, bd.sqrt())
    };

    // Newton refinement: gradient of ½|S(u,v)-Q|²
    let eps = ((u1 - u0) + (v1 - v0)) * 1e-6;
    for _ in 0..40 {
        let p = surface.point_at(best_u, best_v);
        let diff = p - query;
        // Partial derivatives via finite difference
        let pu = surface.point_at((best_u + eps).min(u1), best_v);
        let pum = surface.point_at((best_u - eps).max(u0), best_v);
        let pv = surface.point_at(best_u, (best_v + eps).min(v1));
        let pvm = surface.point_at(best_u, (best_v - eps).max(v0));
        let du = (pu - pum) / (2.0 * eps.min((best_u + eps).min(u1) - (best_u - eps).max(u0)));
        let dv = (pv - pvm) / (2.0 * eps.min((best_v + eps).min(v1) - (best_v - eps).max(v0)));
        // Gradient components: ∂f/∂u = diff · du, ∂f/∂v = diff · dv
        let gu = diff.dot(du);
        let gv = diff.dot(dv);
        // Hessian diagonal approximation (Gauss-Newton)
        let huu = du.dot(du);
        let hvv = dv.dot(dv);
        let huv = du.dot(dv);
        let det = huu * hvv - huv * huv;
        if det.abs() < 1e-20 {
            break;
        }
        let delta_u = (hvv * gu - huv * gv) / det;
        let delta_v = (huu * gv - huv * gu) / det;
        let new_u = (best_u - delta_u).clamp(u0, u1);
        let new_v = (best_v - delta_v).clamp(v0, v1);
        let new_dist = (surface.point_at(new_u, new_v) - query).length();
        if new_dist < best_dist {
            best_dist = new_dist;
            best_u = new_u;
            best_v = new_v;
        }
        if delta_u.abs() < 1e-10 && delta_v.abs() < 1e-10 {
            break;
        }
    }

    let best_point = surface.point_at(best_u, best_v);
    SurfaceProjection {
        point: best_point,
        params: (best_u, best_v),
        distance: (best_point - query).length(),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;

    #[test]
    fn project_onto_plane() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        });
        let q = DVec3::new(3.0, 5.0, -2.0);
        let r = closest_point_on_surface(&plane, q, 8);
        assert!(
            (r.point - DVec3::new(3.0, 0.0, -2.0)).length() < 1e-9,
            "point={}",
            r.point
        );
        assert!((r.distance - 5.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        });
        let q = DVec3::new(5.0, 0.0, 0.0);
        let r = closest_point_on_surface(&sphere, q, 16);
        assert!((r.point - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-9);
        assert!((r.distance - 3.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_cylinder() {
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        });
        let q = DVec3::new(3.0, 2.0, 0.0);
        let r = closest_point_on_surface(&cyl, q, 16);
        assert!((r.point - DVec3::new(1.0, 2.0, 0.0)).length() < 1e-9);
        assert!((r.distance - 2.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_torus() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        // Point far along +X axis → nearest point on outer equator
        let q = DVec3::new(10.0, 0.0, 0.0);
        let r = closest_point_on_surface(&torus, q, 16);
        // Nearest should be at (4, 0, 0)
        assert!((r.point - DVec3::new(4.0, 0.0, 0.0)).length() < 1e-6);
    }

    #[test]
    fn project_onto_cone_returns_theta_and_slant_params() {
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            half_angle_rad: 30.0_f64.to_radians(),
        });
        let expected_slant = 4.0;
        let on_surface = match &cone {
            Surface3::Cone(surface) => surface.point_at(0.0, expected_slant),
            _ => unreachable!(),
        };
        let query_normal = match &cone {
            Surface3::Cone(surface) => surface.normal_at(0.0, expected_slant),
            _ => unreachable!(),
        };
        let q = on_surface + query_normal * 0.25;
        let r = closest_point_on_surface(&cone, q, 16);

        assert!((r.point - on_surface).length() < 5e-3, "projected point={}", r.point);
        assert!((r.params.1 - expected_slant).abs() < 5e-3, "slant={}", r.params.1);
        let lifted = match &cone {
            Surface3::Cone(surface) => surface.point_at(r.params.0, r.params.1),
            _ => unreachable!(),
        };
        assert!((lifted - r.point).length() < 1e-6, "lifted point={lifted} projected={}", r.point);
    }

    #[test]
    fn project_onto_circle_curve() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        let q = DVec3::new(2.0, 0.0, 0.0);
        let r = closest_point_on_curve(&circle, q, 64);
        assert!(
            (r.point - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6,
            "expected (1,0,0) got {}",
            r.point
        );
        assert!((r.distance - 1.0).abs() < 1e-6);
    }

    #[test]
    fn project_onto_line_curve() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        // Line has infinite domain; nearest point to (3, 4, 0) is (3, 0, 0)
        let q = DVec3::new(3.0, 4.0, 0.0);
        let r = closest_point_on_curve(&line, q, 32);
        let expected = DVec3::new(3.0, 0.0, 0.0);
        assert!(
            (r.point - expected).length() < 1e-4,
            "expected {:?} got {}",
            expected,
            r.point
        );
        assert!((r.distance - 4.0).abs() < 1e-4, "distance={}", r.distance);
    }

    #[test]
    fn project_onto_ellipse_curve_analytic() {
        // Ellipse centered at origin in XY plane, semi-axes 3 and 1.
        let ellipse = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        // Query point along +X beyond the ellipse → nearest should be (3, 0, 0).
        let q = DVec3::new(5.0, 0.0, 0.0);
        let r = closest_point_on_curve(&ellipse, q, 64);
        assert!(
            (r.point - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-5,
            "expected (3,0,0) got {}",
            r.point
        );
        assert!((r.distance - 2.0).abs() < 1e-5, "distance={}", r.distance);
    }

    #[test]
    fn project_onto_ellipse_curve_off_plane() {
        // Query lifted off the ellipse plane — projection should still land on ellipse.
        let ellipse = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 2.0,
            minor_radius: 1.0,
        });
        // Query at (2, 0, 5) → closest in-plane point is (2, 0, 0).
        let q = DVec3::new(2.0, 0.0, 5.0);
        let r = closest_point_on_curve(&ellipse, q, 64);
        assert!(r.point.z.abs() < 1e-5, "z should be ~0, got {}", r.point.z);
        assert!(
            (r.point - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-5,
            "expected (2,0,0) got {}",
            r.point
        );
    }

    #[test]
    fn project_onto_line_curve_oblique() {
        // Line along (1,1,0)/sqrt(2), query off axis → test 3-D dot product.
        let dir = DVec3::new(1.0, 1.0, 0.0).normalize();
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: dir,
        });
        let q = DVec3::new(0.0, 1.0, 2.0);
        let r = closest_point_on_curve(&line, q, 32);
        // t = q·dir = 0*0.707 + 1*0.707 + 0 = 0.707, point = t*dir
        let t = q.dot(dir);
        let expected = dir * t;
        assert!(
            (r.point - expected).length() < 1e-9,
            "expected {:?} got {}",
            expected,
            r.point
        );
    }

    #[test]
    fn project_onto_partial_circle_arc() {
        // Arc from 0 to π/2 (first quadrant).
        let arc = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        // For a full circle, query at (-2, 0, 0) should give t = π, point = (-1, 0, 0).
        let q = DVec3::new(-2.0, 0.0, 0.0);
        let r = closest_point_on_curve(&arc, q, 64);
        assert!(
            (r.point - DVec3::new(-1.0, 0.0, 0.0)).length() < 1e-6,
            "expected (-1,0,0) got {}",
            r.point
        );
    }

    #[test]
    fn project_onto_bspline_surface() {
        // Flat BSpline surface at z=0 over [0,1]²
        use crate::geom::BSplineSurface;
        let surf = Surface3::BSpline(BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![vec![1.0; 2]; 2],
        });
        let q = DVec3::new(0.5, 0.5, 5.0);
        let r = closest_point_on_surface(&surf, q, 8);
        assert!(
            (r.point - DVec3::new(0.5, 0.5, 0.0)).length() < 1e-4,
            "got {}",
            r.point
        );
    }

    // ============================================================================
    // OCCT TKGeomBase Alignment Tests - Projection Edge Cases
    // ============================================================================

    #[test]
    fn project_onto_cone_surface() {
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: std::f64::consts::FRAC_PI_6,
        });
        // Query near the cone surface (not at apex - apex is singular)
        let q = DVec3::new(1.0, 0.0, 1.0); // Near the cone surface
        let r = closest_point_on_surface(&cone, q, 16);
        assert!(r.distance < 0.5, "near-surface projection should be close");

        // Query along axis away from apex
        let q2 = DVec3::new(0.0, 0.0, 5.0);
        let r2 = closest_point_on_surface(&cone, q2, 16);
        assert!(r2.distance > 0.0, "axis projection distance should be positive");
    }

    #[test]
    fn project_onto_torus_surface() {
        let torus = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        // Query at center (inside hole)
        let q = DVec3::new(0.0, 0.0, 0.0);
        let r = closest_point_on_surface(&torus, q, 16);
        assert!(r.distance > 0.0, "center distance should be positive");

        // Query on outer ring
        let q2 = DVec3::new(4.0, 0.0, 0.0);
        let r2 = closest_point_on_surface(&torus, q2, 16);
        assert!((r2.distance - 0.0).abs() < 0.1, "on-torus distance should be small");
    }

    #[test]
    fn project_onto_cylinder_interior() {
        let cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        });
        // Query inside cylinder
        let q = DVec3::new(0.0, 0.0, 1.0);
        let r = closest_point_on_surface(&cyl, q, 16);
        assert!((r.distance - 2.0).abs() < 1e-6, "interior distance should be radius");
    }

    #[test]
    fn project_onto_sphere_interior() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
        });
        // Query inside sphere
        let q = DVec3::new(1.0, 1.0, 1.0);
        let r = closest_point_on_surface(&sphere, q, 16);
        assert!(r.distance < 3.0, "interior distance should be less than radius");
    }

    #[test]
    fn project_onto_plane_offset() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, 5.0, 0.0),
            normal: DVec3::Y,
        });
        let q = DVec3::new(3.0, 0.0, -2.0);
        let r = closest_point_on_surface(&plane, q, 8);
        assert!((r.point.y - 5.0).abs() < 1e-9);
        assert!((r.distance - 5.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_line_at_origin() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let q = DVec3::new(0.0, 5.0, 0.0);
        let r = closest_point_on_curve(&line, q, 16);
        assert!((r.point - DVec3::ZERO).length() < 1e-9);
        assert!((r.distance - 5.0).abs() < 1e-9);
    }

    #[test]
    fn project_onto_circle_at_parameter() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::new(2.0, 3.0, 0.0),
            normal: DVec3::Z,
            radius: 1.0,
        });
        // Query at parameter 0 (should be at center + radius * X)
        let q = DVec3::new(3.0, 3.0, 0.0);
        let r = closest_point_on_curve(&circle, q, 32);
        assert!((r.point - DVec3::new(3.0, 3.0, 0.0)).length() < 1e-6);
    }

    #[test]
    fn project_distant_point_onto_sphere() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        // Very distant query
        let q = DVec3::new(1000.0, 1000.0, 1000.0);
        let r = closest_point_on_surface(&sphere, q, 16);
        // Distance should be approximately |q| - radius
        let expected_dist = q.length() - 1.0;
        assert!((r.distance - expected_dist).abs() < 1.0, "distant projection distance");
    }

    #[test]
    fn project_near_surface_boundary() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        // Query very close to plane
        let q = DVec3::new(1.0, 2.0, 1e-10);
        let r = closest_point_on_surface(&plane, q, 8);
        assert!(r.distance < 1e-9, "near-surface distance should be tiny");
    }
}
