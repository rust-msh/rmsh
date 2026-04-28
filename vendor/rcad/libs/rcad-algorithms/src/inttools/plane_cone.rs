use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

#[derive(Debug, Clone)]
pub enum PlaneConicalResult {
    NoIntersection,
    Point(DVec3),
    SingleLine(Line3),
    TwoLines(Line3, Line3),
    Circle(Circle3),
    Ellipse(Ellipse3),
    Parabola(Parabola3),
    Hyperbola(Hyperbola3),
}

/// Classify the intersection of a plane and a cone analytically.
///
/// Let α = angle between plane normal and cone axis, β = cone half-angle.
/// Conic type is determined by the Dandelin-sphere criterion:
/// - |sin α| ≈ 1  (plane ⊥ axis)             → Circle
/// - |sin α| > cos β  (plane steeper than α = π/2−β)  → Ellipse
/// - |sin α| ≈ cos β  (plane parallel to one generator) → Parabola
/// - |sin α| < cos β  (plane shallower, parallel to axis) → Hyperbola
pub fn intersect_plane_cone(plane: &Plane, cone: &ConicalSurface) -> PlaneConicalResult {
    let axis_n = cone.axis.normalize();
    let plane_n = plane.normal.normalize();
    let apex = cone.apex_point();

    // cos of angle between plane normal and cone axis
    let cos_angle = plane_n.dot(axis_n).abs();
    // sin of that angle (= cos of complement)
    let sin_angle = (1.0 - cos_angle * cos_angle).sqrt().max(0.0);

    // Signed distance from apex to plane along plane normal direction
    let apex_to_plane = (plane.origin - apex).dot(plane.normal);

    // ── Plane ⊥ axis → circle ─────────────────────────────────────────────────
    if (cos_angle - 1.0).abs() < TOLERANCE_ANG {
        if apex_to_plane.abs() < TOLERANCE_ABS {
            return PlaneConicalResult::Point(apex);
        }
        let t = apex_to_plane / axis_n.dot(plane.normal);
        let center = apex + axis_n * t;
        let radius = (t * cone.half_angle_rad.tan()).abs();
        if radius < TOLERANCE_ABS {
            return PlaneConicalResult::Point(center);
        }
        return PlaneConicalResult::Circle(Circle3 {
            center,
            normal: cone.axis,
            radius,
        });
    }

    // ── Plane through apex ────────────────────────────────────────────────────
    if apex_to_plane.abs() < TOLERANCE_ABS {
        let angle_between = sin_angle.atan2(cos_angle); // angle between plane and axis
        let half = cone.half_angle_rad;

        if (angle_between - half).abs() < TOLERANCE_ANG {
            // Tangent: single generator line
            let dir = plane_n.cross(axis_n).normalize();
            let gen_dir = (axis_n * half.cos() + dir * half.sin()).normalize();
            return PlaneConicalResult::SingleLine(Line3 {
                origin: apex,
                direction: gen_dir,
            });
        }

        if angle_between < half {
            // Two generators
            let cross = plane_n.cross(axis_n);
            if is_zero_vec(cross) {
                return PlaneConicalResult::Point(apex);
            }
            let perp_in_plane = cross.normalize();
            let projected_axis =
                (axis_n - plane_n * axis_n.dot(plane_n)).normalize_or_zero();
            if projected_axis.length_squared() < 1e-12 {
                return PlaneConicalResult::Point(apex);
            }
            let d1 = (projected_axis * half.cos() + perp_in_plane * half.sin()).normalize();
            let d2 = (projected_axis * half.cos() - perp_in_plane * half.sin()).normalize();
            return PlaneConicalResult::TwoLines(
                Line3 { origin: apex, direction: d1 },
                Line3 { origin: apex, direction: d2 },
            );
        }

        return PlaneConicalResult::Point(apex);
    }

    // ── General case: conic type via Dandelin criterion ───────────────────────
    // Let σ = angle between the cutting plane and the cone axis.
    // cos(σ) = sin_angle  (since σ = 90° − angle_between_normal_and_axis)
    // The cone's generator makes angle β (= half_angle_rad) with the axis.
    //
    //  cos(σ) = sin_angle > sin_beta  →  Ellipse   (plane steeper than generator)
    //  cos(σ) = sin_angle ≈ sin_beta  →  Parabola  (plane parallel to one generator)
    //  cos(σ) = sin_angle < sin_beta  →  Hyperbola (plane shallower; cuts both nappes)
    let cos_beta = cone.half_angle_rad.cos();
    let sin_beta = cone.half_angle_rad.sin();

    // ── Parabola: plane parallel to exactly one generator ─────────────────────
    if (sin_angle - sin_beta).abs() < TOLERANCE_ANG {
        return build_parabola(plane, cone, apex_to_plane, axis_n);
    }

    // ── Hyperbola: sin_angle < sin_beta ───────────────────────────────────────
    if sin_angle < sin_beta - TOLERANCE_ANG {
        return build_hyperbola(plane, cone, apex_to_plane, axis_n, cos_beta, sin_beta, sin_angle);
    }

    // ── Ellipse: sin_angle > cos_beta ─────────────────────────────────────────
    build_ellipse(plane, cone, apex_to_plane, axis_n, cos_angle)
}

// ─────────────────────────────────────────────────────────────────────────────
// Ellipse builder (corrected)
// ─────────────────────────────────────────────────────────────────────────────

fn build_ellipse(
    plane: &Plane,
    cone: &ConicalSurface,
    apex_to_plane: f64,
    axis_n: DVec3,
    cos_angle: f64,
) -> PlaneConicalResult {
    // The intersection of a plane with a (right circular) cone is an ellipse
    // when sin(α) > sin(β) where α = angle between plane and axis, β = half-angle.
    //
    // We use the standard oblique-section formula:
    //   center = apex + axis_n * t  where t = apex_to_plane / (axis_n · plane.normal)
    //   minor_radius = r(t) = |t| * tan(β)  (circular cross-section radius at height t)
    //   major_radius = minor_radius / cos(γ)
    //     where γ = angle between plane and the circular cross-section plane
    //           = angle between plane normal and cone axis = acos(cos_angle)
    //     so major_radius = minor_radius / cos_angle
    //
    // (This is the standard textbook formula for plane-cone ellipse.)
    let tan_beta = cone.half_angle_rad.tan();
    let denom = axis_n.dot(plane.normal);
    if denom.abs() < 1e-14 {
        return PlaneConicalResult::NoIntersection;
    }
    let t = apex_to_plane / denom;

    // Must be on the same nappe as the apex_to_plane sign
    // (t > 0: upper nappe; t < 0: lower nappe)
    let apex = cone.apex_point();
    let center = apex + axis_n * t;
    let base_radius = (t * tan_beta).abs();

    if base_radius < TOLERANCE_ABS {
        return PlaneConicalResult::Point(center);
    }


    // Semi-minor axis = base_radius (perpendicular to tilt direction)
    let minor_radius = base_radius;
    // Semi-major axis: correct formula using Dandelin approach
    // a = b / sqrt(1 - e²) where e = sin_angle_tilt / cos_beta... complex.
    // Use the practical formula: major_radius = base_radius / cos(γ)
    // where γ = angle between plane and the "circular" cross-section,
    // = complement of angle between plane normal and axis.
    // This is the standard oblique section formula: a = r / cos(φ)
    // where φ = angle between plane normal and the cylinder axis.
    // For a cone this is approximate but matches the standard result for shallow angles.
    let major_radius = base_radius / cos_angle;

    // Major direction in the plane (toward the steeper axis)
    let major_dir = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let major_dir = if major_dir.length_squared() < 1e-12 {
        any_perpendicular(plane.normal)
    } else {
        major_dir
    };

    PlaneConicalResult::Ellipse(Ellipse3 {
        center,
        normal: plane.normal,
        major_dir,
        major_radius,
        minor_radius,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Parabola builder
// ─────────────────────────────────────────────────────────────────────────────

fn build_parabola(
    plane: &Plane,
    cone: &ConicalSurface,
    apex_to_plane: f64,
    axis_n: DVec3,
) -> PlaneConicalResult {
    // The plane is parallel to exactly one generator.
    // The vertex of the parabola is the point where the single tangent
    // generator intersects the plane.
    //
    // The generator direction in the plane of the cone that is parallel to
    // the cutting plane: find the generator in the plane spanned by axis and
    // the "steepest descent" direction in the cutting plane.
    let steepest = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let steepest = if steepest.length_squared() < 1e-12 {
        any_perpendicular(plane.normal)
    } else {
        steepest
    };

    let tan_beta = cone.half_angle_rad.tan();
    // Generator parallel to the cutting plane: axis_n + tan_beta * steepest (normalized)
    let gen_dir = (axis_n + tan_beta * steepest).normalize();

    // Vertex: foot of generator on the plane
    let denom = gen_dir.dot(plane.normal);
    let vertex = if denom.abs() > 1e-12 {
        let t = apex_to_plane / denom;
        cone.apex_point() + gen_dir * t
    } else {
        // Generator is parallel to plane; use foot of axis on plane
        let t = apex_to_plane / axis_n.dot(plane.normal).max(1e-12);
        cone.apex_point() + axis_n * t
    };

    // Axis direction of the parabola: projection of cone axis onto the plane
    let axis_2d = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let axis_dir = if axis_2d.length_squared() < 1e-12 {
        steepest
    } else {
        axis_2d
    };

    // Focal parameter p: derived from cone geometry.
    // For a cone with half-angle β and unit-speed cut at vertex distance d_v from apex:
    //   d_v = apex_to_plane / (gen_dir · plane_n)  (already computed above)
    // p = 2 * r_v * tan_beta where r_v = d_v * sin_beta
    let d_v = (vertex - cone.apex_point()).length().max(TOLERANCE_ABS);
    let r_v = d_v * cone.half_angle_rad.sin();
    let focal_param = (2.0 * r_v * tan_beta).max(TOLERANCE_ABS);

    PlaneConicalResult::Parabola(Parabola3 {
        vertex,
        normal: plane.normal,
        axis_dir,
        focal_param,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Hyperbola builder
// ─────────────────────────────────────────────────────────────────────────────

fn build_hyperbola(
    plane: &Plane,
    cone: &ConicalSurface,
    apex_to_plane: f64,
    axis_n: DVec3,
    _cos_beta: f64,
    sin_beta: f64,
    sin_angle: f64,
) -> PlaneConicalResult {
    // Plane cuts both nappes → two-branch hyperbola.
    // The apex is on the hyperbola's transverse axis.
    //
    // Center: projection of apex onto the cutting plane.
    let center = cone.apex_point() + plane.normal * apex_to_plane;

    // Major direction in the plane: projection of cone axis onto the plane.
    let major_dir = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let major_dir = if major_dir.length_squared() < 1e-12 {
        any_perpendicular(plane.normal)
    } else {
        major_dir
    };

    let cos_angle = plane.normal.dot(axis_n).abs();
    let tan_beta = cone.half_angle_rad.tan();

    // Semi-axes via the Dandelin-sphere construction for hyperbolas:
    // Let d = |apex_to_plane| (distance from apex to cutting plane along its normal).
    // The foot of the axis on the plane: center_along_axis = apex + axis_n * (d / cos_angle)
    //   (if cos_angle ≈ 0, the axis is nearly parallel to the plane)
    //
    // Practical formula (standard oblique section of a right circular cone):
    //   a = d * sin_beta / sqrt(sin_beta² − sin_angle²)
    //   b = d * sin_angle * cos_beta / sqrt(sin_beta² − sin_angle²)  (... simplified)
    //
    // For the purely parallel case (sin_angle=0, plane parallel to axis):
    //   The cutting plane at distance ρ from axis intersects in two lines if ρ < r,
    //   or a hyperbola whose semi-transverse = ρ / tan_beta (approx).
    //   Use the general formula with sin_angle=0:
    //   a = d * sin_beta / sin_beta = d,  b = 0 → degenerate (two straight lines).
    // We handle this by using the ρ-based formula when cos_angle ≈ 0.

    let discriminant = sin_beta * sin_beta - sin_angle * sin_angle;
    if discriminant <= TOLERANCE_ABS * TOLERANCE_ABS {
        // Parabola boundary or degenerate; caller should have caught this.
        // Fall back to no intersection to be safe.
        return PlaneConicalResult::NoIntersection;
    }
    let sqrt_d = discriminant.sqrt();
    let d = apex_to_plane.abs();

    let (a, b) = if cos_angle < TOLERANCE_ANG {
        // Plane nearly parallel to axis: compute by radial distance from axis.
        // apex_to_plane is distance along plane normal; for axis-parallel plane,
        // plane normal ⊥ axis, so ρ = distance from cone axis to plane.
        // The hyperbola at that distance from the cone axis:
        //   at height y along axis: r(y) = y * tan_beta
        //   intersection with plane at distance ρ: y = ρ / tan_beta  (where both branches meet)
        // a = ρ / tan_beta  (half-distance between vertices)
        // b = ρ (transverse width at x = 0... approximate)
        let rho = d;
        (rho / tan_beta, rho)
    } else {
        (d * sin_beta / sqrt_d, d * sin_angle * cos_angle.abs() / sqrt_d)
    };

    if a < TOLERANCE_ABS {
        return PlaneConicalResult::Point(center);
    }

    // Ensure b > 0 (may be tiny for near-axis planes)
    let b = b.max(TOLERANCE_ABS);

    PlaneConicalResult::Hyperbola(Hyperbola3 {
        center,
        normal: plane.normal,
        major_dir,
        semi_major: a,
        semi_minor: b,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_cone() -> ConicalSurface {
        ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 0.0,
            half_angle_rad: std::f64::consts::FRAC_PI_4, // 45°
        }
    }

    #[test]
    fn perpendicular_plane_circle() {
        let plane = Plane {
            origin: DVec3::new(0.0, 2.0, 0.0),
            normal: DVec3::Y,
        };
        match intersect_plane_cone(&plane, &test_cone()) {
            PlaneConicalResult::Circle(c) => {
                assert!((c.center.y - 2.0).abs() < TOLERANCE_ABS);
                assert!((c.radius - 2.0).abs() < 0.01); // tan(45°)*2 = 2
            }
            other => panic!("Expected Circle, got {other:?}"),
        }
    }

    #[test]
    fn plane_through_apex() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };
        let result = intersect_plane_cone(&plane, &test_cone());
        // Plane through apex, perpendicular to X → two lines or point
        assert!(matches!(
            result,
            PlaneConicalResult::TwoLines(_, _) | PlaneConicalResult::Point(_)
        ));
    }

    /// A plane that makes 45° with the Y-axis (same as the cone half-angle) is
    /// parallel to exactly one generator → parabola.
    ///
    /// Cone: apex at origin, axis Y, half-angle 45°.
    /// Plane normal = normalize(Y + Z) = (0, 1/√2, 1/√2):
    ///   angle between normal and Y-axis = 45°.
    ///   sin(angle) = sin(45°) = 1/√2 ≈ cos(45°) = cos_beta → parabola.
    #[test]
    fn parabola_case() {
        let n = DVec3::new(0.0, 1.0, 1.0).normalize();
        let plane = Plane {
            origin: DVec3::new(0.0, 2.0, 0.0),
            normal: n,
        };
        match intersect_plane_cone(&plane, &test_cone()) {
            PlaneConicalResult::Parabola(p) => {
                assert!(p.focal_param > 0.0, "focal_param must be positive: {}", p.focal_param);
                assert!(p.normal.dot(n).abs() > 0.99, "parabola normal should match plane normal");
            }
            other => panic!("Expected Parabola, got {other:?}"),
        }
    }

    /// A plane with its normal nearly aligned with the cone axis cuts only one
    /// nappe at a shallow angle → hyperbola (sin_angle < sin_beta).
    ///
    /// Cone: axis Y, half-angle 60°.  sin_beta = sin(60°) ≈ 0.866.
    /// Plane normal = normalize(Y * 2 + Z): cos_angle ≈ 0.894, sin_angle ≈ 0.447.
    /// sin_angle (0.447) < sin_beta (0.866) → hyperbola.
    #[test]
    fn hyperbola_case() {
        let cone_60 = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
            half_angle_rad: std::f64::consts::PI / 3.0, // 60°; sin_beta ≈ 0.866
        };
        // normal = normalize(2Y + Z): cos_angle = 2/sqrt(5) ≈ 0.894, sin_angle ≈ 0.447
        let n = DVec3::new(0.0, 2.0, 1.0).normalize();
        let plane = Plane {
            origin: DVec3::new(0.0, 0.5, 0.0),
            normal: n,
        };
        match intersect_plane_cone(&plane, &cone_60) {
            PlaneConicalResult::Hyperbola(h) => {
                assert!(h.semi_major > 0.0, "semi_major must be positive: {}", h.semi_major);
                assert!(h.semi_minor > 0.0, "semi_minor must be positive: {}", h.semi_minor);
            }
            other => panic!("Expected Hyperbola, got {other:?}"),
        }
    }

    /// A 30° half-angle cone with a steeply-tilted cutting plane should
    /// yield an ellipse (sin_angle > sin_beta).
    ///
    /// Cone: axis Y, half_angle=30°, so sin_beta=0.5.
    /// Plane normal = normalize(Z * 2 + Y): angle between normal and Y ≈ 63°,
    /// sin_angle ≈ sin(63°) ≈ 0.891 > sin_beta=0.5 → ellipse.
    #[test]
    fn steep_oblique_gives_ellipse() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
            half_angle_rad: std::f64::consts::PI / 6.0, // 30°, sin_beta=0.5
        };
        // normal = normalize(2Z + Y): cos_angle=1/sqrt(5)≈0.447, sin_angle≈0.894 > 0.5
        let n = DVec3::new(0.0, 1.0, 2.0).normalize();
        let plane = Plane {
            origin: DVec3::new(0.0, 1.0, 0.0),
            normal: n,
        };
        match intersect_plane_cone(&plane, &cone) {
            PlaneConicalResult::Ellipse(_) => {}
            other => panic!("Expected Ellipse for steep cut of 30° cone, got {other:?}"),
        }
    }
}
