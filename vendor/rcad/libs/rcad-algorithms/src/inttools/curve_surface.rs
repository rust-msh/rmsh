use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

/// Hit from a curve-surface intersection.
pub struct CurveSurfaceHit {
    pub point: DVec3,
    /// Parametric value on the curve.
    pub curve_param: f64,
}

/// Intersect a line with a cylindrical surface (infinite cylinder).
/// Returns 0, 1, or 2 hit points within the parameter range.
pub fn intersect_line_cylinder(
    line: &Line3,
    t_range: [f64; 2],
    cyl: &CylindricalSurface,
) -> Vec<CurveSurfaceHit> {
    // Project to 2D perpendicular to cylinder axis
    let oc = line.origin - cyl.origin;

    // Component of direction perpendicular to axis
    let d_perp = line.direction - cyl.axis * line.direction.dot(cyl.axis);
    let oc_perp = oc - cyl.axis * oc.dot(cyl.axis);

    let a = d_perp.dot(d_perp);
    let b = 2.0 * oc_perp.dot(d_perp);
    let c = oc_perp.dot(oc_perp) - cyl.radius * cyl.radius;

    solve_quadratic_hits(a, b, c, line, t_range)
}

/// Intersect a line with a sphere.
pub fn intersect_line_sphere(
    line: &Line3,
    t_range: [f64; 2],
    sphere: &SphericalSurface,
) -> Vec<CurveSurfaceHit> {
    let oc = line.origin - sphere.center;
    let a = line.direction.dot(line.direction);
    let b = 2.0 * oc.dot(line.direction);
    let c = oc.length_squared() - sphere.radius * sphere.radius;

    solve_quadratic_hits(a, b, c, line, t_range)
}

/// Intersect a line with a conical surface (infinite cone).
pub fn intersect_line_cone(
    line: &Line3,
    t_range: [f64; 2],
    cone: &ConicalSurface,
) -> Vec<CurveSurfaceHit> {
    let apex = cone.apex_point();
    let axis = cone.axis_dir();
    let co = line.origin - apex;
    let cos2 = cone.half_angle_rad.cos().powi(2);

    let d_dot_a = line.direction.dot(axis);
    let co_dot_a = co.dot(axis);

    // Point P on cone satisfies:
    //   ((P-apex)·axis)² = cos²(half_angle) * |P-apex|²
    // Substituting P = O + t*D:
    let d_d = line.direction.dot(line.direction);
    let co_d = co.dot(line.direction);
    let co_co = co.dot(co);

    let a = d_dot_a * d_dot_a - cos2 * d_d;
    let b = 2.0 * (d_dot_a * co_dot_a - cos2 * co_d);
    let c = co_dot_a * co_dot_a - cos2 * co_co;

    // Filter hits: only keep those where the point is on the correct nappe
    // (same side as cone.axis from apex), and the cone has positive height direction.
    let mut hits = solve_quadratic_hits(a, b, c, line, t_range);
    hits.retain(|hit| {
        let v = hit.point - apex;
        v.dot(axis) > -TOLERANCE_ABS // on the positive nappe
    });
    hits
}

/// Intersect a circle arc with a plane. Returns 0, 1, or 2 hits.
pub fn intersect_circle_plane(
    circle: &Circle3,
    t_range: [f64; 2], // angular range in radians
    plane: &Plane,
) -> Vec<CurveSurfaceHit> {
    // Circle parametric: P(θ) = center + radius*(u*cos(θ) + v*sin(θ))
    // where u,v are orthonormal vectors in the circle plane
    let u = if circle.normal.x.abs() < 0.9 {
        circle.normal.cross(DVec3::X).normalize()
    } else {
        circle.normal.cross(DVec3::Y).normalize()
    };
    let v = circle.normal.cross(u);

    // Plane equation: (P - plane.origin) · plane.normal = 0
    // Substituting: (center - plane.origin)·n + r*(u·n*cos(θ) + v·n*sin(θ)) = 0
    let d = (circle.center - plane.origin).dot(plane.normal);
    let a_coeff = circle.radius * u.dot(plane.normal);
    let b_coeff = circle.radius * v.dot(plane.normal);

    // d + a_coeff*cos(θ) + b_coeff*sin(θ) = 0
    // A*cos(θ) + B*sin(θ) = -d
    // R*cos(θ - φ) = -d  where R = sqrt(A²+B²)
    let r_amp = (a_coeff * a_coeff + b_coeff * b_coeff).sqrt();
    if r_amp < TOLERANCE_ABS {
        return vec![]; // circle parallel to plane
    }

    let ratio = -d / r_amp;
    if ratio.abs() > 1.0 + TOLERANCE_ABS {
        return vec![];
    }
    let ratio = ratio.clamp(-1.0, 1.0);

    let phi = b_coeff.atan2(a_coeff);
    let alpha = ratio.acos();

    let mut hits = Vec::new();
    for theta in [phi + alpha, phi - alpha] {
        // Normalize theta to [0, 2π)
        let theta = ((theta % (2.0 * std::f64::consts::PI)) + 2.0 * std::f64::consts::PI)
            % (2.0 * std::f64::consts::PI);

        if theta >= t_range[0] - TOLERANCE_ABS && theta <= t_range[1] + TOLERANCE_ABS {
            let point = circle.center
                + u * (circle.radius * theta.cos())
                + v * (circle.radius * theta.sin());
            hits.push(CurveSurfaceHit {
                point,
                curve_param: theta,
            });
        }
    }
    hits
}

/// Intersect a circle arc with a cylindrical surface. Returns 0–4 hits.
///
/// Strategy: substitute circle parametric `P(θ) = C + r*(u*cosθ + v*sinθ)` into
/// the cylinder implicit `|P⊥ - cyl.origin⊥|² = R²` (projecting out the axis).
/// This yields `A + B*cosθ + C*sinθ + D*cos2θ + E*sin2θ = 0`, which we solve
/// numerically via Newton refinement seeded from a coarse angle grid.
pub fn intersect_circle_cylinder(
    circle: &Circle3,
    t_range: [f64; 2],
    cyl: &CylindricalSurface,
) -> Vec<CurveSurfaceHit> {
    circle_vs_implicit_surface(
        circle,
        t_range,
        |p: DVec3| -> f64 {
            let v = p - cyl.origin;
            let along = v.dot(cyl.axis);
            let perp = v - cyl.axis * along;
            perp.length_squared() - cyl.radius * cyl.radius
        },
    )
}

/// Intersect a circle arc with a spherical surface. Returns 0–2 hits.
///
/// Substitutes circle parametric into sphere implicit `|P - center|² = R²`.
pub fn intersect_circle_sphere(
    circle: &Circle3,
    t_range: [f64; 2],
    sph: &SphericalSurface,
) -> Vec<CurveSurfaceHit> {
    circle_vs_implicit_surface(
        circle,
        t_range,
        |p: DVec3| -> f64 { (p - sph.center).length_squared() - sph.radius * sph.radius },
    )
}

/// Intersect a circle arc with a conical surface. Returns 0–4 hits.
pub fn intersect_circle_cone(
    circle: &Circle3,
    t_range: [f64; 2],
    cone: &ConicalSurface,
) -> Vec<CurveSurfaceHit> {
    let cos2 = cone.half_angle_rad.cos().powi(2);
    let apex = cone.apex_point();
    let axis = cone.axis_dir();
    circle_vs_implicit_surface(
        circle,
        t_range,
        |p: DVec3| -> f64 {
            let v = p - apex;
            let along = v.dot(axis);
            let along2 = along * along;
            let len2 = v.length_squared();
            // Cone implicit: (v·axis)² = cos²(half) * |v|²
            along2 - cos2 * len2
        },
    )
}

/// Generic circle-vs-implicit-surface intersection via Newton refinement.
///
/// Evaluates `f(P(θ)) = 0` on a circle arc.  Seeds are found by coarse sampling;
/// each seed is refined with Newton's method.  Duplicate roots within TOLERANCE_ABS
/// are deduplicated.
fn circle_vs_implicit_surface(
    circle: &Circle3,
    t_range: [f64; 2],
    f: impl Fn(DVec3) -> f64,
) -> Vec<CurveSurfaceHit> {
    use std::f64::consts::TAU;

    // Build a local orthonormal frame for the circle
    let cn = circle.normal.normalize();
    let cu = if cn.x.abs() < 0.9 {
        cn.cross(DVec3::X).normalize()
    } else {
        cn.cross(DVec3::Y).normalize()
    };
    let cv = cn.cross(cu);

    let pt = |theta: f64| -> DVec3 {
        circle.center + circle.radius * (theta.cos() * cu + theta.sin() * cv)
    };

    const N_SEEDS: usize = 64;
    let [t0, t1] = t_range;
    let span = (t1 - t0).min(TAU);

    // Sign-change detection over coarse grid
    let mut seeds: Vec<f64> = Vec::new();
    let mut prev_val = f(pt(t0));
    for i in 1..=N_SEEDS {
        let theta = t0 + span * i as f64 / N_SEEDS as f64;
        let val = f(pt(theta));
        if prev_val * val <= 0.0 {
            // Sign change — midpoint as seed
            seeds.push(theta - span * 0.5 / N_SEEDS as f64);
        }
        prev_val = val;
    }

    // Newton refinement
    let mut hits: Vec<CurveSurfaceHit> = Vec::new();
    const MAX_ITER: usize = 20;
    const H: f64 = 1e-7;
    for seed in seeds {
        let mut theta = seed;
        for _ in 0..MAX_ITER {
            let fv = f(pt(theta));
            let dfdtheta = (f(pt(theta + H)) - f(pt(theta - H))) / (2.0 * H);
            if dfdtheta.abs() < 1e-30 {
                break;
            }
            let delta = -fv / dfdtheta;
            theta += delta;
            if delta.abs() < TOLERANCE_ABS * 0.01 {
                break;
            }
        }

        // Validate within t_range and on the surface
        if theta < t0 - TOLERANCE_ABS || theta > t1 + TOLERANCE_ABS {
            continue;
        }
        let point = pt(theta);
        if f(point).abs() > TOLERANCE_ABS * 10.0 {
            continue;
        }

        // Deduplicate
        let duplicate = hits.iter().any(|h: &CurveSurfaceHit| {
            (h.curve_param - theta).abs() < TOLERANCE_ABS * 5.0
        });
        if !duplicate {
            hits.push(CurveSurfaceHit {
                point,
                curve_param: theta,
            });
        }
    }
    hits
}

fn solve_quadratic_hits(
    a: f64,
    b: f64,
    c: f64,
    line: &Line3,
    t_range: [f64; 2],
) -> Vec<CurveSurfaceHit> {
    let mut hits = Vec::new();

    if a.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
        // Linear
        if b.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
            return hits;
        }
        let t = -c / b;
        if t >= t_range[0] - TOLERANCE_ABS && t <= t_range[1] + TOLERANCE_ABS {
            hits.push(CurveSurfaceHit {
                point: line.origin + line.direction * t,
                curve_param: t,
            });
        }
        return hits;
    }

    let discriminant = b * b - 4.0 * a * c;
    if discriminant < -TOLERANCE_ABS {
        return hits;
    }

    if discriminant.abs() < TOLERANCE_ABS {
        let t = -b / (2.0 * a);
        if t >= t_range[0] - TOLERANCE_ABS && t <= t_range[1] + TOLERANCE_ABS {
            hits.push(CurveSurfaceHit {
                point: line.origin + line.direction * t,
                curve_param: t,
            });
        }
    } else {
        let sqrt_d = discriminant.sqrt();
        for t in [(-b - sqrt_d) / (2.0 * a), (-b + sqrt_d) / (2.0 * a)] {
            if t >= t_range[0] - TOLERANCE_ABS && t <= t_range[1] + TOLERANCE_ABS {
                hits.push(CurveSurfaceHit {
                    point: line.origin + line.direction * t,
                    curve_param: t,
                });
            }
        }
    }
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_through_sphere() {
        let line = Line3 {
            origin: DVec3::new(-5.0, 0.0, 0.0),
            direction: DVec3::X,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        };
        let hits = intersect_line_sphere(&line, [-10.0, 10.0], &sphere);
        assert_eq!(hits.len(), 2);
        // Should hit at x = -1 and x = 1
        let xs: Vec<f64> = hits.iter().map(|h| h.point.x).collect();
        assert!((xs[0] - (-1.0)).abs() < TOLERANCE_ABS);
        assert!((xs[1] - 1.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn line_tangent_to_sphere() {
        let line = Line3 {
            origin: DVec3::new(-5.0, 1.0, 0.0),
            direction: DVec3::X,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        };
        let hits = intersect_line_sphere(&line, [-10.0, 10.0], &sphere);
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn line_misses_sphere() {
        let line = Line3 {
            origin: DVec3::new(-5.0, 2.0, 0.0),
            direction: DVec3::X,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        };
        let hits = intersect_line_sphere(&line, [-10.0, 10.0], &sphere);
        assert!(hits.is_empty());
    }

    #[test]
    fn line_through_cylinder() {
        let line = Line3 {
            origin: DVec3::new(-5.0, 0.5, 0.0),
            direction: DVec3::X,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        };
        let hits = intersect_line_cylinder(&line, [-10.0, 10.0], &cyl);
        assert_eq!(hits.len(), 2);
        // Hits at x = ±1, z = 0
        for h in &hits {
            let dist = ((h.point.x * h.point.x + h.point.z * h.point.z).sqrt() - 1.0).abs();
            assert!(dist < TOLERANCE_ABS, "point not on cylinder: dist={dist}");
        }
    }

    #[test]
    fn line_misses_cylinder() {
        let line = Line3 {
            origin: DVec3::new(-5.0, 0.0, 3.0),
            direction: DVec3::X,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        };
        let hits = intersect_line_cylinder(&line, [-10.0, 10.0], &cyl);
        assert!(hits.is_empty());
    }

    #[test]
    fn line_through_cone() {
        let line = Line3 {
            origin: DVec3::new(-5.0, 1.0, 0.0),
            direction: DVec3::X,
        };
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
            half_angle_rad: std::f64::consts::FRAC_PI_4,
        };
        let hits = intersect_line_cone(&line, [-10.0, 10.0], &cone);
        // At y=1, the cone radius = tan(45°) * 1 = 1, so hits at x = ±1
        assert_eq!(hits.len(), 2);
        for h in &hits {
            assert!((h.point.y - 1.0).abs() < TOLERANCE_ABS);
        }
    }
}
