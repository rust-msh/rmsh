//! IntAna-style analytical intersection algorithms.
//!
//! This module provides analytical solutions for geometric intersections between
//! curves (lines) and surfaces (planes, cylinders, spheres, cones, tori), as well
//! as surface-surface intersections.
//!
//! Named after OCCT's IntAna package which provides similar functionality.

use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

// =============================================================================
// Line-Surface Intersections
// =============================================================================

/// Result of a line-plane intersection.
///
/// Returns the intersection point and the parameter value on the line where
/// the intersection occurs: `P = line.origin + t * line.direction`.
#[derive(Debug, Clone, Copy)]
pub struct LinPlnIntersection {
    /// The intersection point in 3D space.
    pub point: DVec3,
    /// The parameter on the line where intersection occurs.
    pub param: f64,
}

/// Intersect a line with a plane.
///
/// # Returns
/// - `Some((point, param))` if the line intersects the plane
/// - `None` if the line is parallel to the plane
///
/// # Example
/// ```
/// use glam::DVec3;
/// use rcad_kernel::geom::{Line3, Plane};
/// use rcad_algorithms::int_ana::intersect_line_plane;
///
/// let line = Line3 { origin: DVec3::ZERO, direction: DVec3::Z };
/// let plane = Plane { origin: DVec3::new(0.0, 0.0, 5.0), normal: DVec3::Z };
/// let result = intersect_line_plane(&line, &plane);
/// assert!(result.is_some());
/// let intersection = result.unwrap();
/// assert!((intersection.param - 5.0).abs() < 1e-10);
/// ```
pub fn intersect_line_plane(line: &Line3, plane: &Plane) -> Option<LinPlnIntersection> {
    let denom = line.direction.dot(plane.normal);

    if denom.abs() < TOLERANCE_ANG {
        // Line is parallel to plane
        return None;
    }

    let t = (plane.origin - line.origin).dot(plane.normal) / denom;
    let point = line.origin + t * line.direction;

    Some(LinPlnIntersection { point, param: t })
}

// -----------------------------------------------------------------------------

/// Intersect a line with a cylinder.
///
/// A line can intersect a cylinder at 0, 1 (tangent), or 2 points.
///
/// # Returns
/// A vector of intersection points, each with the line parameter.
/// The vector is sorted by parameter value.
pub fn intersect_line_cylinder(line: &Line3, cyl: &CylindricalSurface) -> Vec<(DVec3, f64)> {
    // Transform line into cylinder's local coordinate system where:
    // - Cylinder axis is along Z
    // - Cylinder origin is at (0, 0, 0)
    // - Cylinder equation: x^2 + y^2 = r^2

    let axis = cyl.axis.normalize();
    let origin = cyl.origin;

    // Build orthonormal basis for cylinder
    let x_axis = any_perpendicular(axis);
    let y_axis = axis.cross(x_axis).normalize();

    // Transform line origin and direction to cylinder frame
    let rel_origin = line.origin - origin;
    let o = DVec3::new(
        rel_origin.dot(x_axis),
        rel_origin.dot(y_axis),
        rel_origin.dot(axis),
    );
    let d = DVec3::new(
        line.direction.dot(x_axis),
        line.direction.dot(y_axis),
        line.direction.dot(axis),
    );

    // Solve: (ox + t*dx)^2 + (oy + t*dy)^2 = r^2
    // => (dx^2 + dy^2) t^2 + 2(ox*dx + oy*dy) t + (ox^2 + oy^2 - r^2) = 0
    let a = d.x * d.x + d.y * d.y;
    let b = 2.0 * (o.x * d.x + o.y * d.y);
    let c = o.x * o.x + o.y * o.y - cyl.radius * cyl.radius;

    // If a is nearly zero, line is parallel to cylinder axis
    if a.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
        // Check if line is on the cylinder surface
        if c.abs() < TOLERANCE_ABS {
            // Line is on cylinder surface (coincident) - return no discrete points
            return vec![];
        }
        // Line is parallel but not on surface
        return vec![];
    }

    let discriminant = b * b - 4.0 * a * c;

    if discriminant < -TOLERANCE_ABS {
        // No intersection
        return vec![];
    }

    if discriminant.abs() < TOLERANCE_ABS {
        // Tangent (one intersection point)
        let t = -b / (2.0 * a);
        let point = line.origin + t * line.direction;
        return vec![(point, t)];
    }

    // Two intersection points
    let sqrt_disc = discriminant.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);

    let p1 = line.origin + t1 * line.direction;
    let p2 = line.origin + t2 * line.direction;

    if t1 <= t2 {
        vec![(p1, t1), (p2, t2)]
    } else {
        vec![(p2, t2), (p1, t1)]
    }
}

// -----------------------------------------------------------------------------

/// Intersect a line with a sphere.
///
/// A line can intersect a sphere at 0, 1 (tangent), or 2 points.
///
/// # Returns
/// A vector of intersection points, each with the line parameter.
/// The vector is sorted by parameter value.
pub fn intersect_line_sphere(line: &Line3, sphere: &SphericalSurface) -> Vec<(DVec3, f64)> {
    // Solve: |P - C|^2 = r^2
    // |O + t*D - C|^2 = r^2
    // Let L = O - C
    // |L + t*D|^2 = r^2
    // L^2 + 2t(L.D) + t^2*D^2 = r^2
    // t^2*D^2 + 2t(L.D) + (L^2 - r^2) = 0

    let l = line.origin - sphere.center;
    let d = line.direction;

    let a = d.dot(d);
    let b = 2.0 * l.dot(d);
    let c = l.dot(l) - sphere.radius * sphere.radius;

    let discriminant = b * b - 4.0 * a * c;

    if discriminant < -TOLERANCE_ABS {
        return vec![];
    }

    if discriminant.abs() < TOLERANCE_ABS {
        // Tangent point
        let t = -b / (2.0 * a);
        let point = line.origin + t * line.direction;
        return vec![(point, t)];
    }

    // Two intersection points
    let sqrt_disc = discriminant.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);

    let p1 = line.origin + t1 * line.direction;
    let p2 = line.origin + t2 * line.direction;

    if t1 <= t2 {
        vec![(p1, t1), (p2, t2)]
    } else {
        vec![(p2, t2), (p1, t1)]
    }
}

// -----------------------------------------------------------------------------

/// Intersect a line with a cone.
///
/// A line can intersect a cone at 0, 1, or 2 points.
///
/// The cone is defined by its apex, axis, radius at reference point, and half-angle.
///
/// # Returns
/// A vector of intersection points, each with the line parameter.
pub fn intersect_line_cone(line: &Line3, cone: &ConicalSurface) -> Vec<(DVec3, f64)> {
    // Cone equation in local coordinates (axis along Z, apex at origin):
    // (x^2 + y^2) = (z * tan(half_angle))^2
    //
    // We transform the line to the cone's local frame and solve the resulting
    // quadratic equation.

    let axis = cone.axis_dir();
    let apex = cone.apex_point();
    let tan_half = cone.half_angle_rad.tan();

    // Build orthonormal basis
    let x_axis = any_perpendicular(axis);
    let y_axis = axis.cross(x_axis).normalize();

    // Transform line to cone's local frame (apex at origin, axis along Z)
    let rel_origin = line.origin - apex;
    let o = DVec3::new(
        rel_origin.dot(x_axis),
        rel_origin.dot(y_axis),
        rel_origin.dot(axis),
    );
    let d = DVec3::new(
        line.direction.dot(x_axis),
        line.direction.dot(y_axis),
        line.direction.dot(axis),
    );

    // Cone equation: x^2 + y^2 = (z * tan_half)^2
    // Line: (ox + t*dx, oy + t*dy, oz + t*dz)
    //
    // Substitute and expand:
    // (ox + t*dx)^2 + (oy + t*dy)^2 = tan_half^2 * (oz + t*dz)^2
    //
    // Let T = tan_half
    // (dx^2 + dy^2 - T^2*dz^2) t^2 + 2(ox*dx + oy*dy - T^2*oz*dz) t + (ox^2 + oy^2 - T^2*oz^2) = 0

    let t2 = tan_half * tan_half;

    let a = d.x * d.x + d.y * d.y - t2 * d.z * d.z;
    let b = 2.0 * (o.x * d.x + o.y * d.y - t2 * o.z * d.z);
    let c = o.x * o.x + o.y * o.y - t2 * o.z * o.z;

    // Check if the line is parallel to the cone axis
    if a.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
        // Line is parallel to cone axis or nearly so
        // Check if it intersects
        if c.abs() < TOLERANCE_ABS {
            // Line passes through apex
            return vec![(apex, -o.z / d.z.max(1e-12))];
        }
        // Check if the line hits the cone surface
        // At height z from apex, cone radius is |z| * tan_half
        // Line at height z has radial distance sqrt((ox + (z-oz)/dz*dx)^2 + (oy + (z-oz)/dz*dy)^2)
        // This is complex; fall back to sampling for now
        if d.z.abs() > TOLERANCE_ABS {
            // Line is along axis direction, check if it's within the cone
            let radial_dist_sq = o.x * o.x + o.y * o.y;
            let z_at_origin = o.z;
            let cone_radius_at_origin = (z_at_origin * tan_half).abs();
            if radial_dist_sq < (cone_radius_at_origin * cone_radius_at_origin) + TOLERANCE_ABS {
                // Line is inside cone - no intersection (passes through apex if radial = 0)
                return vec![];
            }
        }
        return vec![];
    }

    let discriminant = b * b - 4.0 * a * c;

    if discriminant < -TOLERANCE_ABS {
        return vec![];
    }

    if discriminant.abs() < TOLERANCE_ABS {
        let t = -b / (2.0 * a);
        let point = line.origin + t * line.direction;
        return vec![(point, t)];
    }

    let sqrt_disc = discriminant.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2_param = (-b + sqrt_disc) / (2.0 * a);

    let p1 = line.origin + t1 * line.direction;
    let p2 = line.origin + t2_param * line.direction;

    // Verify points are on the correct nappe (z >= 0 or z <= 0 depending on sign)
    // by checking that the point satisfies the cone equation with consistent sign
    let mut results = Vec::new();

    for (p, t) in [(p1, t1), (p2, t2_param)] {
        let rel_p = p - apex;
        let z = rel_p.dot(axis);
        let radial = (rel_p - z * axis).length();
        let expected_radius = (z * tan_half).abs();

        // Accept point if it's close to the cone surface
        if (radial - expected_radius).abs() < TOLERANCE_ABS * 100.0 {
            results.push((p, t));
        }
    }

    results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    results
}

// -----------------------------------------------------------------------------

/// Intersect a line with a torus.
///
/// A line can intersect a torus at 0, 1, 2, 3, or 4 points.
///
/// The torus is defined by its center, axis, major radius (distance from center
/// to tube center), and minor radius (tube radius).
///
/// # Returns
/// A vector of intersection points, each with the line parameter.
pub fn intersect_line_torus(line: &Line3, torus: &ToroidalSurface) -> Vec<(DVec3, f64)> {
    // Torus equation in local coordinates (axis along Z):
    // (sqrt(x^2 + y^2) - R)^2 + z^2 = r^2
    // where R = major_radius, r = minor_radius
    //
    // Expanding: (x^2 + y^2 + z^2 + R^2 - r^2)^2 = 4R^2(x^2 + y^2)
    //
    // This is a quartic equation in t when we substitute the line equation.

    let axis = torus.axis.normalize();
    let center = torus.center;
    let r_major = torus.major_radius;
    let r_minor = torus.minor_radius;

    // Build orthonormal basis
    let x_axis = any_perpendicular(axis);
    let y_axis = axis.cross(x_axis).normalize();

    // Transform line to torus's local frame
    let rel_origin = line.origin - center;
    let o = DVec3::new(
        rel_origin.dot(x_axis),
        rel_origin.dot(y_axis),
        rel_origin.dot(axis),
    );
    let d = DVec3::new(
        line.direction.dot(x_axis),
        line.direction.dot(y_axis),
        line.direction.dot(axis),
    );

    // Quartic equation: a4*t^4 + a3*t^3 + a2*t^2 + a1*t + a0 = 0
    // From the torus equation: (x^2 + y^2 + z^2 + R^2 - r^2)^2 - 4R^2(x^2 + y^2) = 0
    //
    // Let f = o + t*d, and let:
    //   s = x^2 + y^2 + z^2 = |f|^2
    //   p = x^2 + y^2 (distance from axis squared)
    //
    // Then: (s + R^2 - r^2)^2 - 4R^2*p = 0

    let r2 = r_minor * r_minor;
    let r_major2 = r_major * r_major;
    let alpha = r_major2 - r2; // R^2 - r^2

    // s(t) = |o + t*d|^2 = o^2 + 2t(o.d) + t^2*d^2
    // p(t) = (ox + t*dx)^2 + (oy + t*dy)^2

    let o_sq = o.dot(o);
    let o_dot_d = o.dot(d);
    let d_sq = d.dot(d);

    let ox_sq = o.x * o.x;
    let oy_sq = o.y * o.y;

    let p0 = ox_sq + oy_sq; // p(0)
    let p1 = 2.0 * (o.x * d.x + o.y * d.y); // coefficient of t in p(t)
    let p2 = d.x * d.x + d.y * d.y; // coefficient of t^2 in p(t)

    // s(t) = o_sq + 2*o_dot_d*t + d_sq*t^2
    // s(t) + alpha = (o_sq + alpha) + 2*o_dot_d*t + d_sq*t^2
    //
    // (s + alpha)^2 = [(o_sq + alpha) + 2*o_dot_d*t + d_sq*t^2]^2
    //
    // 4*R^2*p(t) = 4*R^2*(p0 + p1*t + p2*t^2)

    let s0 = o_sq + alpha;
    let s1 = 2.0 * o_dot_d;
    let s2 = d_sq;

    // Expand (s0 + s1*t + s2*t^2)^2:
    // = s0^2 + 2*s0*s1*t + (2*s0*s2 + s1^2)*t^2 + 2*s1*s2*t^3 + s2^2*t^4

    let a4 = s2 * s2;
    let a3 = 2.0 * s1 * s2;
    let a2 = 2.0 * s0 * s2 + s1 * s1;
    let a1 = 2.0 * s0 * s1;
    let a0 = s0 * s0;

    // Subtract 4*R^2*p(t):
    let four_r_major2 = 4.0 * r_major2;
    let a2_final = a2 - four_r_major2 * p2;
    let a1_final = a1 - four_r_major2 * p1;
    let a0_final = a0 - four_r_major2 * p0;

    // Solve the quartic
    let roots = solve_quartic(a4, a3, a2_final, a1_final, a0_final);

    // Convert roots to points
    let mut results: Vec<(DVec3, f64)> = roots
        .into_iter()
        .map(|t| {
            let point = line.origin + t * line.direction;
            (point, t)
        })
        .collect();

    results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    results
}

// =============================================================================
// Plane-Surface Intersections
// =============================================================================

/// Result of plane-plane intersection.
#[derive(Debug, Clone)]
pub enum PlnPlnResult {
    /// Planes are parallel and do not intersect.
    Parallel,
    /// Planes are coincident (same plane).
    Coincident,
    /// Planes intersect in a line.
    Line(Line3),
}

/// Intersect two planes.
///
/// Two planes either:
/// - Are parallel (no intersection)
/// - Are coincident (same plane)
/// - Intersect in a line
pub fn intersect_plane_plane_intana(p1: &Plane, p2: &Plane) -> PlnPlnResult {
    let n1 = p1.normal;
    let n2 = p2.normal;
    let cross = n1.cross(n2);

    if is_zero_vec(cross) {
        let d = (p2.origin - p1.origin).dot(n1);
        if d.abs() < TOLERANCE_ABS {
            return PlnPlnResult::Coincident;
        }
        return PlnPlnResult::Parallel;
    }

    let direction = cross.normalize();

    // Find a point on the intersection line
    // The point satisfies both plane equations: n1.P = d1, n2.P = d2
    // where d = n.origin (distance from origin along normal)
    let d1 = n1.dot(p1.origin);
    let d2 = n2.dot(p2.origin);

    let origin = solve_two_plane_point(n1, d1, n2, d2, direction);

    PlnPlnResult::Line(Line3 { origin, direction })
}

/// Find a point on the intersection line of two planes by zeroing the largest
/// component of the line direction and solving the resulting 2x2 system.
fn solve_two_plane_point(n1: DVec3, d1: f64, n2: DVec3, d2: f64, dir: DVec3) -> DVec3 {
    let abs_dir = DVec3::new(dir.x.abs(), dir.y.abs(), dir.z.abs());

    if abs_dir.x >= abs_dir.y && abs_dir.x >= abs_dir.z {
        // Set x = 0
        let det = n1.y * n2.z - n1.z * n2.y;
        let y = (d1 * n2.z - d2 * n1.z) / det;
        let z = (n1.y * d2 - n2.y * d1) / det;
        DVec3::new(0.0, y, z)
    } else if abs_dir.y >= abs_dir.z {
        // Set y = 0
        let det = n1.x * n2.z - n1.z * n2.x;
        let x = (d1 * n2.z - d2 * n1.z) / det;
        let z = (n1.x * d2 - n2.x * d1) / det;
        DVec3::new(x, 0.0, z)
    } else {
        // Set z = 0
        let det = n1.x * n2.y - n1.y * n2.x;
        let x = (d1 * n2.y - d2 * n1.y) / det;
        let y = (n1.x * d2 - n2.x * d1) / det;
        DVec3::new(x, y, 0.0)
    }
}

// -----------------------------------------------------------------------------

/// Result of plane-cylinder intersection.
#[derive(Debug, Clone)]
pub enum PlnCylResult {
    /// No intersection (plane outside cylinder).
    NoIntersection,
    /// Single tangent line.
    TangentLine(Line3),
    /// Two parallel lines (plane parallel to axis, cutting through).
    TwoLines(Line3, Line3),
    /// Circle (plane perpendicular to axis).
    Circle(Circle3),
    /// Ellipse (plane at angle to axis).
    Ellipse(Ellipse3),
}

/// Intersect a plane with a cylinder.
pub fn intersect_plane_cylinder_intana(plane: &Plane, cyl: &CylindricalSurface) -> PlnCylResult {
    let cos_angle = plane.normal.dot(cyl.axis).abs();

    if cos_angle < TOLERANCE_ANG {
        // Plane parallel to cylinder axis
        let axis_to_plane = (plane.origin - cyl.origin).dot(plane.normal);
        let dist = axis_to_plane.abs();

        if dist > cyl.radius + TOLERANCE_ABS {
            return PlnCylResult::NoIntersection;
        }
        if (dist - cyl.radius).abs() < TOLERANCE_ABS {
            let tang_point = cyl.origin + plane.normal * (-axis_to_plane);
            return PlnCylResult::TangentLine(Line3 {
                origin: tang_point,
                direction: cyl.axis,
            });
        }
        let offset_dir = plane.normal.cross(cyl.axis).normalize();
        let half_chord = (cyl.radius * cyl.radius - dist * dist).sqrt();
        let center_on_plane = cyl.origin - plane.normal * axis_to_plane;

        let l1_origin = center_on_plane + offset_dir * half_chord;
        let l2_origin = center_on_plane - offset_dir * half_chord;

        return PlnCylResult::TwoLines(
            Line3 {
                origin: l1_origin,
                direction: cyl.axis,
            },
            Line3 {
                origin: l2_origin,
                direction: cyl.axis,
            },
        );
    }

    if (cos_angle - 1.0).abs() < TOLERANCE_ANG {
        // Plane perpendicular to cylinder axis -> circle
        let t = (plane.origin - cyl.origin).dot(cyl.axis);
        let center = cyl.origin + cyl.axis * t;
        return PlnCylResult::Circle(Circle3 {
            center,
            normal: cyl.axis,
            radius: cyl.radius,
        });
    }

    // General oblique case -> ellipse
    let major_radius = cyl.radius / cos_angle;
    let minor_radius = cyl.radius;

    let t = (plane.origin - cyl.origin).dot(plane.normal) / cyl.axis.dot(plane.normal);
    let center = cyl.origin + cyl.axis * t;

    let major_dir = (cyl.axis - plane.normal * cyl.axis.dot(plane.normal)).normalize();

    PlnCylResult::Ellipse(Ellipse3 {
        center,
        normal: plane.normal,
        major_dir,
        major_radius,
        minor_radius,
    })
}

// -----------------------------------------------------------------------------

/// Result of plane-sphere intersection.
#[derive(Debug, Clone)]
pub enum PlnSphResult {
    /// No intersection (plane outside sphere).
    NoIntersection,
    /// Single tangent point.
    TangentPoint(DVec3),
    /// Circle intersection.
    Circle(Circle3),
}

/// Intersect a plane with a sphere.
pub fn intersect_plane_sphere_intana(plane: &Plane, sphere: &SphericalSurface) -> PlnSphResult {
    let signed_dist = (sphere.center - plane.origin).dot(plane.normal);
    let abs_dist = signed_dist.abs();

    if abs_dist > sphere.radius + TOLERANCE_ABS {
        return PlnSphResult::NoIntersection;
    }
    if (abs_dist - sphere.radius).abs() < TOLERANCE_ABS {
        let point = sphere.center - plane.normal * signed_dist;
        return PlnSphResult::TangentPoint(point);
    }

    let circle_radius = (sphere.radius * sphere.radius - signed_dist * signed_dist).sqrt();
    let center = sphere.center - plane.normal * signed_dist;

    PlnSphResult::Circle(Circle3 {
        center,
        normal: plane.normal,
        radius: circle_radius,
    })
}

// -----------------------------------------------------------------------------

/// Result of plane-cone intersection.
#[derive(Debug, Clone)]
pub enum PlnConResult {
    /// No intersection.
    NoIntersection,
    /// Single point (plane through apex at tangent angle).
    Point(DVec3),
    /// Single generator line (plane tangent to cone).
    SingleLine(Line3),
    /// Two generator lines (plane through apex).
    TwoLines(Line3, Line3),
    /// Circle (plane perpendicular to axis).
    Circle(Circle3),
    /// Ellipse (plane cutting at angle, sin_angle > sin_beta).
    Ellipse(Ellipse3),
    /// Parabola (plane parallel to one generator).
    Parabola(Parabola3),
    /// Hyperbola (plane cutting both nappes).
    Hyperbola(Hyperbola3),
}

/// Intersect a plane with a cone.
///
/// The intersection type depends on the angle between the plane normal and
/// the cone axis relative to the cone's half-angle:
/// - Perpendicular to axis: circle
/// - Steep angle (sin_angle > sin_beta): ellipse
/// - Parallel to generator (sin_angle = sin_beta): parabola
/// - Shallow angle (sin_angle < sin_beta): hyperbola
pub fn intersect_plane_cone_intana(plane: &Plane, cone: &ConicalSurface) -> PlnConResult {
    let axis_n = cone.axis_dir();
    let plane_n = plane.normal.normalize();
    let apex = cone.apex_point();

    // cos of angle between plane normal and cone axis
    let cos_angle = plane_n.dot(axis_n).abs();
    // sin of that angle
    let sin_angle = (1.0 - cos_angle * cos_angle).sqrt().max(0.0);

    // Signed distance from apex to plane along plane normal direction
    let apex_to_plane = (plane.origin - apex).dot(plane.normal);

    // Plane perpendicular to axis -> circle
    if (cos_angle - 1.0).abs() < TOLERANCE_ANG {
        if apex_to_plane.abs() < TOLERANCE_ABS {
            return PlnConResult::Point(apex);
        }
        let t = apex_to_plane / axis_n.dot(plane.normal);
        let center = apex + axis_n * t;
        let radius = (t * cone.half_angle_rad.tan()).abs();
        if radius < TOLERANCE_ABS {
            return PlnConResult::Point(center);
        }
        return PlnConResult::Circle(Circle3 {
            center,
            normal: cone.axis,
            radius,
        });
    }

    // Plane through apex
    if apex_to_plane.abs() < TOLERANCE_ABS {
        let angle_between = sin_angle.atan2(cos_angle);
        let half = cone.half_angle_rad;

        if (angle_between - half).abs() < TOLERANCE_ANG {
            // Tangent: single generator line
            let dir = plane_n.cross(axis_n).normalize();
            let gen_dir = (axis_n * half.cos() + dir * half.sin()).normalize();
            return PlnConResult::SingleLine(Line3 {
                origin: apex,
                direction: gen_dir,
            });
        }

        if angle_between < half {
            // Two generators
            let cross = plane_n.cross(axis_n);
            if is_zero_vec(cross) {
                return PlnConResult::Point(apex);
            }
            let perp_in_plane = cross.normalize();
            let projected_axis = (axis_n - plane_n * axis_n.dot(plane_n)).normalize_or_zero();
            if projected_axis.length_squared() < 1e-12 {
                return PlnConResult::Point(apex);
            }
            let d1 = (projected_axis * half.cos() + perp_in_plane * half.sin()).normalize();
            let d2 = (projected_axis * half.cos() - perp_in_plane * half.sin()).normalize();
            return PlnConResult::TwoLines(
                Line3 { origin: apex, direction: d1 },
                Line3 { origin: apex, direction: d2 },
            );
        }

        return PlnConResult::Point(apex);
    }

    // General case: conic type via Dandelin criterion
    let sin_beta = cone.half_angle_rad.sin();

    // Parabola: plane parallel to exactly one generator
    if (sin_angle - sin_beta).abs() < TOLERANCE_ANG {
        return build_parabola_result(plane, cone, apex_to_plane, axis_n);
    }

    // Hyperbola: sin_angle < sin_beta
    if sin_angle < sin_beta - TOLERANCE_ANG {
        return build_hyperbola_result(plane, cone, apex_to_plane, axis_n, sin_angle);
    }

    // Ellipse: sin_angle > sin_beta
    build_ellipse_result(plane, cone, apex_to_plane, axis_n, cos_angle)
}

fn build_ellipse_result(
    plane: &Plane,
    cone: &ConicalSurface,
    apex_to_plane: f64,
    axis_n: DVec3,
    cos_angle: f64,
) -> PlnConResult {
    let tan_beta = cone.half_angle_rad.tan();
    let denom = axis_n.dot(plane.normal);
    if denom.abs() < 1e-14 {
        return PlnConResult::NoIntersection;
    }
    let t = apex_to_plane / denom;

    let apex = cone.apex_point();
    let center = apex + axis_n * t;
    let base_radius = (t * tan_beta).abs();

    if base_radius < TOLERANCE_ABS {
        return PlnConResult::Point(center);
    }

    let minor_radius = base_radius;
    let major_radius = base_radius / cos_angle;

    let major_dir = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let major_dir = if major_dir.length_squared() < 1e-12 {
        any_perpendicular(plane.normal)
    } else {
        major_dir
    };

    PlnConResult::Ellipse(Ellipse3 {
        center,
        normal: plane.normal,
        major_dir,
        major_radius,
        minor_radius,
    })
}

fn build_parabola_result(
    plane: &Plane,
    cone: &ConicalSurface,
    apex_to_plane: f64,
    axis_n: DVec3,
) -> PlnConResult {
    let steepest = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let steepest = if steepest.length_squared() < 1e-12 {
        any_perpendicular(plane.normal)
    } else {
        steepest
    };

    let tan_beta = cone.half_angle_rad.tan();
    let gen_dir = (axis_n + tan_beta * steepest).normalize();

    let denom = gen_dir.dot(plane.normal);
    let vertex = if denom.abs() > 1e-12 {
        let t = apex_to_plane / denom;
        cone.apex_point() + gen_dir * t
    } else {
        let t = apex_to_plane / axis_n.dot(plane.normal).max(1e-12);
        cone.apex_point() + axis_n * t
    };

    let axis_2d = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let axis_dir = if axis_2d.length_squared() < 1e-12 {
        steepest
    } else {
        axis_2d
    };

    let d_v = (vertex - cone.apex_point()).length().max(TOLERANCE_ABS);
    let r_v = d_v * cone.half_angle_rad.sin();
    let focal_param = (2.0 * r_v * tan_beta).max(TOLERANCE_ABS);

    PlnConResult::Parabola(Parabola3 {
        vertex,
        normal: plane.normal,
        axis_dir,
        focal_param,
    })
}

fn build_hyperbola_result(
    plane: &Plane,
    cone: &ConicalSurface,
    apex_to_plane: f64,
    axis_n: DVec3,
    sin_angle: f64,
) -> PlnConResult {
    let center = cone.apex_point() + plane.normal * apex_to_plane;

    let major_dir = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let major_dir = if major_dir.length_squared() < 1e-12 {
        any_perpendicular(plane.normal)
    } else {
        major_dir
    };

    let cos_angle = plane.normal.dot(axis_n).abs();
    let sin_beta = cone.half_angle_rad.sin();
    let tan_beta = cone.half_angle_rad.tan();

    let discriminant = sin_beta * sin_beta - sin_angle * sin_angle;
    if discriminant <= TOLERANCE_ABS * TOLERANCE_ABS {
        return PlnConResult::NoIntersection;
    }
    let sqrt_d = discriminant.sqrt();
    let d = apex_to_plane.abs();

    let (a, b) = if cos_angle < TOLERANCE_ANG {
        let rho = d;
        (rho / tan_beta, rho)
    } else {
        (d * sin_beta / sqrt_d, d * sin_angle * cos_angle.abs() / sqrt_d)
    };

    if a < TOLERANCE_ABS {
        return PlnConResult::Point(center);
    }

    let b = b.max(TOLERANCE_ABS);

    PlnConResult::Hyperbola(Hyperbola3 {
        center,
        normal: plane.normal,
        major_dir,
        semi_major: a,
        semi_minor: b,
    })
}

// =============================================================================
// Cylinder-Cylinder Intersection
// =============================================================================

/// Result of cylinder-cylinder intersection.
#[derive(Debug, Clone)]
pub enum CylCylResult {
    /// No intersection (cylinders too far apart).
    NoIntersection,
    /// Single curve (tangent case).
    SingleCurve(Curve3),
    /// Two intersection curves.
    TwoCurves(Curve3, Curve3),
    /// Complex intersection requiring numerical methods.
    Complex(Vec<Curve3>),
}

/// Intersect two cylinders.
///
/// Cylinder-cylinder intersection is complex and can produce:
/// - Circles (coaxial, same radius)
/// - Ellipses (rare, special cases)
/// - Bicurves (general case, typically 3D curves)
pub fn intersect_cylinder_cylinder(cyl1: &CylindricalSurface, cyl2: &CylindricalSurface) -> CylCylResult {
    // Check for parallel axes
    let axis1 = cyl1.axis.normalize();
    let axis2 = cyl2.axis.normalize();

    let axes_parallel = vectors_parallel(axis1, axis2);

    if axes_parallel {
        return intersect_parallel_cylinders(cyl1, cyl2, axis1);
    }

    // Check for perpendicular axes (special case)
    let dot_axes = axis1.dot(axis2);
    if dot_axes.abs() < TOLERANCE_ANG {
        return intersect_perpendicular_cylinders(cyl1, cyl2, axis1, axis2);
    }

    // General case: axes are skew or intersecting at an angle
    // The intersection is typically a spatial curve (not a simple conic)
    intersect_general_cylinders(cyl1, cyl2, axis1, axis2)
}

fn intersect_parallel_cylinders(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    axis: DVec3,
) -> CylCylResult {
    // Distance between axes
    let diff = cyl2.origin - cyl1.origin;
    let dist_along_axis = diff.dot(axis);
    let radial_vec = diff - dist_along_axis * axis;
    let dist_between_axes = radial_vec.length();

    let r1 = cyl1.radius;
    let r2 = cyl2.radius;

    // Check for intersection
    if dist_between_axes > r1 + r2 + TOLERANCE_ABS {
        return CylCylResult::NoIntersection;
    }

    if (dist_between_axes - (r1 + r2)).abs() < TOLERANCE_ABS {
        // External tangent - single line of contact (rare)
        // This would be a line, not a curve
        let contact_point = cyl1.origin + radial_vec.normalize() * r1;
        return CylCylResult::SingleCurve(Curve3::Line(Line3 {
            origin: contact_point,
            direction: axis,
        }));
    }

    if dist_between_axes < (r1 - r2).abs() - TOLERANCE_ABS {
        // One cylinder is inside the other
        if dist_between_axes + r2 < r1 {
            return CylCylResult::NoIntersection; // cyl2 entirely inside cyl1
        }
    }

    // Concentric same-radius cylinders
    if dist_between_axes < TOLERANCE_ABS && (r1 - r2).abs() < TOLERANCE_ABS {
        // Same cylinder - return a representative circle
        return CylCylResult::SingleCurve(Curve3::Circle(Circle3 {
            center: cyl1.origin,
            normal: axis,
            radius: r1,
        }));
    }

    // Two parallel cylinders with overlapping cross-sections
    // Intersection consists of two lines parallel to the axis
    if dist_between_axes < TOLERANCE_ABS {
        // Coaxial different radii - no intersection if radii differ
        return CylCylResult::NoIntersection;
    }

    // Compute the two intersection lines
    let dir = radial_vec.normalize();
    let d = dist_between_axes;
    let x = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
    let h_sq = r1 * r1 - x * x;

    if h_sq < -TOLERANCE_ABS {
        return CylCylResult::NoIntersection;
    }

    if h_sq < TOLERANCE_ABS {
        // Single tangent line
        let contact = cyl1.origin + dir * x;
        return CylCylResult::SingleCurve(Curve3::Line(Line3 {
            origin: contact,
            direction: axis,
        }));
    }

    let h = h_sq.sqrt();
    let perp = axis.cross(dir).normalize();

    let p1 = cyl1.origin + dir * x + perp * h;
    let p2 = cyl1.origin + dir * x - perp * h;

    CylCylResult::TwoCurves(
        Curve3::Line(Line3 { origin: p1, direction: axis }),
        Curve3::Line(Line3 { origin: p2, direction: axis }),
    )
}

fn intersect_perpendicular_cylinders(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    axis1: DVec3,
    axis2: DVec3,
) -> CylCylResult {
    // Build a local frame: axis1 along Z, axis2 along X (after projection)
    let x_axis = any_perpendicular(axis1);
    let y_axis = axis1.cross(x_axis).normalize();

    // Transform cyl2's origin to cyl1's frame
    let rel_origin = cyl2.origin - cyl1.origin;
    let o2 = DVec3::new(
        rel_origin.dot(x_axis),
        rel_origin.dot(y_axis),
        rel_origin.dot(axis1),
    );

    // For perpendicular cylinders with offset d in the xy-plane,
    // the intersection curve is symmetric about both axes.
    // The intersection can be computed numerically or approximated.

    // Check if the cylinders actually intersect
    // The distance from cyl2's axis to cyl1's axis in the xy-plane
    let d_xy = DVec3::new(o2.x, o2.y, 0.0).length();

    if d_xy > cyl1.radius + cyl2.radius + TOLERANCE_ABS {
        return CylCylResult::NoIntersection;
    }

    // For the general perpendicular case, we'd need to compute the 3D intersection curve
    // This is complex and typically requires numerical methods or B-spline approximation
    intersect_general_cylinders(cyl1, cyl2, axis1, axis2)
}

fn intersect_general_cylinders(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    axis1: DVec3,
    axis2: DVec3,
) -> CylCylResult {
    // Find the closest points between the two axes
    let diff = cyl2.origin - cyl1.origin;

    // Direction of the line perpendicular to both axes
    let cross = axis1.cross(axis2);
    let cross_len = cross.length();

    if cross_len < TOLERANCE_ANG {
        // Parallel axes (should have been handled above)
        return intersect_parallel_cylinders(cyl1, cyl2, axis1);
    }

    let cross_dir = cross / cross_len;

    // Distance between axes at their closest approach
    let d = diff.dot(cross_dir).abs();

    // Check if cylinders intersect at all
    // This is a simplified check - actual intersection depends on the geometry
    if d > cyl1.radius + cyl2.radius + TOLERANCE_ABS {
        return CylCylResult::NoIntersection;
    }

    // The intersection of two general cylinders is a complex 3D curve
    // that cannot be expressed as a simple conic. It typically requires
    // numerical approximation or B-spline representation.
    //
    // For now, we return Complex with an empty vector, indicating that
    // numerical methods are required. A full implementation would sample
    // the intersection curve and fit a B-spline.
    CylCylResult::Complex(vec![])
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Solve a quartic equation a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 = 0
/// Returns real roots.
fn solve_quartic(a4: f64, a3: f64, a2: f64, a1: f64, a0: f64) -> Vec<f64> {
    // Normalize to monic polynomial
    if a4.abs() < 1e-15 {
        return solve_cubic(a3, a2, a1, a0);
    }

    let a = a3 / a4;
    let b = a2 / a4;
    let c = a1 / a4;
    let d = a0 / a4;

    // Use Ferrari's method
    // Substitute x = y - a/4 to eliminate cubic term
    // y^4 + p*y^2 + q*y + r = 0
    let p = b - 3.0 * a * a / 8.0;
    let q = c + a * a * a / 8.0 - a * b / 2.0;
    let r = d - 3.0 * a * a * a * a / 256.0 + a * a * b / 16.0 - a * c / 4.0;

    // If q = 0, we have a quadratic in y^2
    if q.abs() < 1e-12 {
        let disc = p * p - 4.0 * r;
        if disc < -TOLERANCE_ABS {
            return vec![];
        }
        if disc.abs() < TOLERANCE_ABS {
            let y = (-p / 2.0).sqrt();
            let x = y - a / 4.0;
            let mut roots = vec![x, -x];
            roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            return roots;
        }
        let sqrt_disc = disc.sqrt();
        let y1_sq = (-p + sqrt_disc) / 2.0;
        let y2_sq = (-p - sqrt_disc) / 2.0;

        let mut roots = Vec::new();
        if y1_sq >= -TOLERANCE_ABS {
            let y1 = y1_sq.max(0.0).sqrt();
            roots.push(y1 - a / 4.0);
            roots.push(-y1 - a / 4.0);
        }
        if y2_sq >= -TOLERANCE_ABS {
            let y2 = y2_sq.max(0.0).sqrt();
            roots.push(y2 - a / 4.0);
            roots.push(-y2 - a / 4.0);
        }
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        return roots;
    }

    // Solve the resolvent cubic: z^3 + (p/2)*z^2 + ((p^2-4r)/16)*z - q^2/64 = 0
    // Let u = z + p/2: u^3 - (p/2)*u^2 - r*u + (p*r - q^2/4)/2 = 0
    // Simplified resolvent: t^3 + 2p*t^2 + (p^2 - 4r)*t - q^2 = 0
    let cubic_a = 1.0;
    let cubic_b = 2.0 * p;
    let cubic_c = p * p - 4.0 * r;
    let cubic_d = -q * q;

    let cubic_roots = solve_cubic(cubic_a, cubic_b, cubic_c, cubic_d);

    // Find a positive real root for t
    let t = cubic_roots
        .iter()
        .find(|&&r| r > TOLERANCE_ABS)
        .copied()
        .unwrap_or_else(|| {
            cubic_roots
                .iter()
                .find(|&&r| r > -TOLERANCE_ABS)
                .copied()
                .unwrap_or(0.0)
        });

    let sqrt_t = t.max(0.0).sqrt();

    // The four roots are given by:
    // y = (sqrt_t + sqrt(-(p + t + q/sqrt_t))) / 2
    // y = (sqrt_t - sqrt(-(p + t + q/sqrt_t))) / 2
    // y = (-sqrt_t + sqrt(-(p + t - q/sqrt_t))) / 2
    // y = (-sqrt_t - sqrt(-(p + t - q/sqrt_t))) / 2

    let mut roots = Vec::new();

    if sqrt_t > TOLERANCE_ABS {
        let inner1 = -(p + t + q / sqrt_t);
        let inner2 = -(p + t - q / sqrt_t);

        if inner1 >= -TOLERANCE_ABS {
            let s1 = inner1.max(0.0).sqrt();
            roots.push((sqrt_t + s1) / 2.0 - a / 4.0);
            roots.push((sqrt_t - s1) / 2.0 - a / 4.0);
        }
        if inner2 >= -TOLERANCE_ABS {
            let s2 = inner2.max(0.0).sqrt();
            roots.push((-sqrt_t + s2) / 2.0 - a / 4.0);
            roots.push((-sqrt_t - s2) / 2.0 - a / 4.0);
        }
    } else {
        // t is nearly zero, use alternative formula
        let inner = -(p + t);
        if inner >= -TOLERANCE_ABS {
            let s = inner.max(0.0).sqrt();
            roots.push(s / 2.0 - a / 4.0);
            roots.push(-s / 2.0 - a / 4.0);
        }
    }

    // Filter out duplicates and near-duplicates
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    roots.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE_ABS);

    roots
}

/// Solve a cubic equation a*x^3 + b*x^2 + c*x + d = 0
/// Returns real roots.
fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    if a.abs() < 1e-15 {
        return solve_quadratic(b, c, d);
    }

    // Normalize to monic: x^3 + px^2 + qx + r = 0
    let p = b / a;
    let q = c / a;
    let r = d / a;

    // Substitute x = t - p/3 to get depressed cubic: t^3 + at + b = 0
    let a_depressed = q - p * p / 3.0;
    let b_depressed = 2.0 * p * p * p / 27.0 - p * q / 3.0 + r;

    // Discriminant
    let discriminant = b_depressed * b_depressed / 4.0 + a_depressed * a_depressed * a_depressed / 27.0;

    let offset = p / 3.0;

    if discriminant > TOLERANCE_ABS {
        // One real root
        let sqrt_disc = discriminant.sqrt();
        let u = (-b_depressed / 2.0 + sqrt_disc).cbrt();
        let v = (-b_depressed / 2.0 - sqrt_disc).cbrt();
        vec![u + v - offset]
    } else if discriminant < -TOLERANCE_ABS {
        // Three distinct real roots
        let radius = (-a_depressed * a_depressed * a_depressed / 27.0).sqrt();
        let theta = (-b_depressed / (2.0 * radius)).acos() / 3.0;
        let cbrt_r = radius.cbrt();

        vec![
            2.0 * cbrt_r * theta.cos() - offset,
            2.0 * cbrt_r * (theta + std::f64::consts::FRAC_PI_2 * (4.0 / 3.0)).cos() - offset,
            2.0 * cbrt_r * (theta + std::f64::consts::FRAC_PI_2 * (8.0 / 3.0)).cos() - offset,
        ]
    } else {
        // One or two real roots (repeated)
        if a_depressed.abs() < TOLERANCE_ABS {
            vec![-offset]
        } else {
            let u = (-b_depressed / 2.0).cbrt();
            vec![2.0 * u - offset, -u - offset]
        }
    }
}

/// Solve a quadratic equation a*x^2 + b*x + c = 0
/// Returns real roots.
fn solve_quadratic(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() < 1e-15 {
        if b.abs() < 1e-15 {
            return vec![];
        }
        return vec![-c / b];
    }

    let discriminant = b * b - 4.0 * a * c;

    if discriminant < -TOLERANCE_ABS {
        return vec![];
    }

    if discriminant.abs() < TOLERANCE_ABS {
        return vec![-b / (2.0 * a)];
    }

    let sqrt_disc = discriminant.sqrt();
    vec![
        (-b - sqrt_disc) / (2.0 * a),
        (-b + sqrt_disc) / (2.0 * a),
    ]
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------------------------------------------------------
    // Line-Plane Tests
    // -------------------------------------------------------------------------

    #[test]
    fn line_plane_simple_intersection() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::Z,
        };
        let plane = Plane {
            origin: DVec3::new(0.0, 0.0, 5.0),
            normal: DVec3::Z,
        };
        let result = intersect_line_plane(&line, &plane).unwrap();
        assert!((result.param - 5.0).abs() < TOLERANCE_ABS);
        assert!((result.point - DVec3::new(0.0, 0.0, 5.0)).length() < TOLERANCE_ABS);
    }

    #[test]
    fn line_plane_parallel() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let plane = Plane {
            origin: DVec3::new(0.0, 0.0, 5.0),
            normal: DVec3::Z,
        };
        assert!(intersect_line_plane(&line, &plane).is_none());
    }

    #[test]
    fn line_plane_oblique() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::new(1.0, 1.0, 1.0).normalize(),
        };
        let plane = Plane {
            origin: DVec3::new(0.0, 0.0, 3.0),
            normal: DVec3::Z,
        };
        let result = intersect_line_plane(&line, &plane).unwrap();
        assert!((result.point.z - 3.0).abs() < TOLERANCE_ABS);
    }

    // -------------------------------------------------------------------------
    // Line-Cylinder Tests
    // -------------------------------------------------------------------------

    #[test]
    fn line_cylinder_through_axis() {
        let line = Line3 {
            origin: DVec3::new(-10.0, 0.0, 0.0),
            direction: DVec3::X,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        let results = intersect_line_cylinder(&line, &cyl);
        assert_eq!(results.len(), 2);
        assert!((results[0].0.x + 2.0).abs() < TOLERANCE_ABS);
        assert!((results[1].0.x - 2.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn line_cylinder_tangent() {
        let line = Line3 {
            origin: DVec3::new(-10.0, 2.0, 0.0),
            direction: DVec3::X,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        let results = intersect_line_cylinder(&line, &cyl);
        assert_eq!(results.len(), 1);
        assert!((results[0].0.y - 2.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn line_cylinder_no_intersection() {
        let line = Line3 {
            origin: DVec3::new(-10.0, 5.0, 0.0),
            direction: DVec3::X,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        let results = intersect_line_cylinder(&line, &cyl);
        assert!(results.is_empty());
    }

    // -------------------------------------------------------------------------
    // Line-Sphere Tests
    // -------------------------------------------------------------------------

    #[test]
    fn line_sphere_through_center() {
        let line = Line3 {
            origin: DVec3::new(-10.0, 0.0, 0.0),
            direction: DVec3::X,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        };
        let results = intersect_line_sphere(&line, &sphere);
        assert_eq!(results.len(), 2);
        assert!((results[0].0.x + 3.0).abs() < TOLERANCE_ABS);
        assert!((results[1].0.x - 3.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn line_sphere_tangent() {
        let line = Line3 {
            origin: DVec3::new(-10.0, 3.0, 0.0),
            direction: DVec3::X,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        };
        let results = intersect_line_sphere(&line, &sphere);
        assert_eq!(results.len(), 1);
        assert!((results[0].0.y - 3.0).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn line_sphere_no_intersection() {
        let line = Line3 {
            origin: DVec3::new(-10.0, 10.0, 0.0),
            direction: DVec3::X,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        };
        let results = intersect_line_sphere(&line, &sphere);
        assert!(results.is_empty());
    }

    // -------------------------------------------------------------------------
    // Line-Cone Tests
    // -------------------------------------------------------------------------

    #[test]
    fn line_cone_through_surface() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 0.0,
            half_angle_rad: std::f64::consts::FRAC_PI_4,
        };
        let line = Line3 {
            origin: DVec3::new(-10.0, 5.0, 0.0),
            direction: DVec3::X,
        };
        let results = intersect_line_cone(&line, &cone);
        // At y=5, cone radius = 5*tan(45) = 5, so line at x=+/-5
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn line_cone_through_apex() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 0.0,
            half_angle_rad: std::f64::consts::FRAC_PI_4,
        };
        let line = Line3 {
            origin: DVec3::new(-10.0, -10.0, 0.0),
            direction: DVec3::new(1.0, 1.0, 0.0).normalize(),
        };
        let results = intersect_line_cone(&line, &cone);
        // Line passes through apex at origin
        assert!(!results.is_empty());
    }

    // -------------------------------------------------------------------------
    // Line-Torus Tests
    // -------------------------------------------------------------------------

    #[test]
    fn line_torus_through_hole() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        let line = Line3 {
            origin: DVec3::new(-10.0, 0.0, 0.0),
            direction: DVec3::X,
        };
        let results = intersect_line_torus(&line, &torus);
        // Line passes through torus center
        // At z=0, the torus has a hole of radius 5-1=4, so line at |x| < 4 doesn't intersect
        // At z=0, the torus surface is at |x| = 4 (inner) or |x| = 6 (outer)
        // So we expect 4 intersections: at x = -6, -4, 4, 6
        assert!(results.len() >= 2, "Expected at least 2 intersections, got {}", results.len());
    }

    #[test]
    fn line_torus_axial() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        let line = Line3 {
            origin: DVec3::new(5.0, 0.0, -10.0),
            direction: DVec3::Z,
        };
        let results = intersect_line_torus(&line, &torus);
        // Line at radius 5 (tube center), passes through the tube
        // Expect 2 intersections at z = -1 and z = 1
        assert_eq!(results.len(), 2);
        assert!((results[0].0.z + 1.0).abs() < TOLERANCE_ABS);
        assert!((results[1].0.z - 1.0).abs() < TOLERANCE_ABS);
    }

    // -------------------------------------------------------------------------
    // Plane-Plane Tests
    // -------------------------------------------------------------------------

    #[test]
    fn plane_plane_intersection() {
        let p1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let p2 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        };
        match intersect_plane_plane_intana(&p1, &p2) {
            PlnPlnResult::Line(line) => {
                assert!(vectors_parallel(line.direction, DVec3::X));
                assert!(line.origin.y.abs() < TOLERANCE_ABS);
                assert!(line.origin.z.abs() < TOLERANCE_ABS);
            }
            _ => panic!("Expected Line"),
        }
    }

    #[test]
    fn plane_plane_parallel() {
        let p1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let p2 = Plane {
            origin: DVec3::new(0.0, 0.0, 5.0),
            normal: DVec3::Z,
        };
        assert!(matches!(
            intersect_plane_plane_intana(&p1, &p2),
            PlnPlnResult::Parallel
        ));
    }

    #[test]
    fn plane_plane_coincident() {
        let p1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let p2 = Plane {
            origin: DVec3::new(5.0, 3.0, 0.0),
            normal: DVec3::Z,
        };
        assert!(matches!(
            intersect_plane_plane_intana(&p1, &p2),
            PlnPlnResult::Coincident
        ));
    }

    // -------------------------------------------------------------------------
    // Plane-Cylinder Tests
    // -------------------------------------------------------------------------

    #[test]
    fn plane_cylinder_circle() {
        let plane = Plane {
            origin: DVec3::new(0.0, 5.0, 0.0),
            normal: DVec3::Y,
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        };
        match intersect_plane_cylinder_intana(&plane, &cyl) {
            PlnCylResult::Circle(c) => {
                assert!((c.radius - 3.0).abs() < TOLERANCE_ABS);
                assert!((c.center.y - 5.0).abs() < TOLERANCE_ABS);
            }
            _ => panic!("Expected Circle"),
        }
    }

    #[test]
    fn plane_cylinder_ellipse() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(0.0, 1.0, 1.0).normalize(),
        };
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        };
        match intersect_plane_cylinder_intana(&plane, &cyl) {
            PlnCylResult::Ellipse(e) => {
                assert!((e.minor_radius - 2.0).abs() < TOLERANCE_ABS);
                assert!(e.major_radius > 2.0);
            }
            _ => panic!("Expected Ellipse"),
        }
    }

    // -------------------------------------------------------------------------
    // Plane-Sphere Tests
    // -------------------------------------------------------------------------

    #[test]
    fn plane_sphere_circle() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 5.0,
        };
        match intersect_plane_sphere_intana(&plane, &sphere) {
            PlnSphResult::Circle(c) => {
                assert!((c.radius - 5.0).abs() < TOLERANCE_ABS);
            }
            _ => panic!("Expected Circle"),
        }
    }

    #[test]
    fn plane_sphere_tangent() {
        let plane = Plane {
            origin: DVec3::new(0.0, 5.0, 0.0),
            normal: DVec3::Y,
        };
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 5.0,
        };
        assert!(matches!(
            intersect_plane_sphere_intana(&plane, &sphere),
            PlnSphResult::TangentPoint(_)
        ));
    }

    // -------------------------------------------------------------------------
    // Plane-Cone Tests
    // -------------------------------------------------------------------------

    #[test]
    fn plane_cone_circle() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 0.0,
            half_angle_rad: std::f64::consts::FRAC_PI_4,
        };
        let plane = Plane {
            origin: DVec3::new(0.0, 4.0, 0.0),
            normal: DVec3::Y,
        };
        match intersect_plane_cone_intana(&plane, &cone) {
            PlnConResult::Circle(c) => {
                assert!((c.center.y - 4.0).abs() < TOLERANCE_ABS);
                assert!((c.radius - 4.0).abs() < 0.01); // tan(45)*4 = 4
            }
            _ => panic!("Expected Circle"),
        }
    }

    #[test]
    fn plane_cone_ellipse() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 0.0,
            half_angle_rad: std::f64::consts::PI / 6.0, // 30 degrees
        };
        let plane = Plane {
            origin: DVec3::new(0.0, 3.0, 0.0),
            normal: DVec3::new(0.0, 1.0, 2.0).normalize(),
        };
        match intersect_plane_cone_intana(&plane, &cone) {
            PlnConResult::Ellipse(_) => {}
            other => panic!("Expected Ellipse, got {:?}", other),
        }
    }

    // -------------------------------------------------------------------------
    // Cylinder-Cylinder Tests
    // -------------------------------------------------------------------------

    #[test]
    fn cylinder_cylinder_parallel_no_intersection() {
        let cyl1 = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        let cyl2 = CylindricalSurface {
            origin: DVec3::new(10.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 2.0,
        };
        assert!(matches!(
            intersect_cylinder_cylinder(&cyl1, &cyl2),
            CylCylResult::NoIntersection
        ));
    }

    #[test]
    fn cylinder_cylinder_parallel_two_lines() {
        let cyl1 = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
        };
        let cyl2 = CylindricalSurface {
            origin: DVec3::new(2.0, 0.0, 0.0),
            axis: DVec3::Z,
            radius: 3.0,
        };
        match intersect_cylinder_cylinder(&cyl1, &cyl2) {
            CylCylResult::TwoCurves(_, _) => {}
            other => panic!("Expected TwoCurves, got {:?}", other),
        }
    }

    #[test]
    fn cylinder_cylinder_perpendicular() {
        let cyl1 = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        let cyl2 = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::X,
            radius: 2.0,
        };
        // Perpendicular intersecting cylinders
        let result = intersect_cylinder_cylinder(&cyl1, &cyl2);
        // Should have some intersection
        match result {
            CylCylResult::Complex(_) | CylCylResult::TwoCurves(_, _) => {}
            other => panic!("Expected intersection, got {:?}", other),
        }
    }

    // -------------------------------------------------------------------------
    // Quartic Solver Tests
    // -------------------------------------------------------------------------

    #[test]
    fn quartic_four_real_roots() {
        // (x-1)(x-2)(x-3)(x-4) = x^4 - 10x^3 + 35x^2 - 50x + 24
        let roots = solve_quartic(1.0, -10.0, 35.0, -50.0, 24.0);
        assert_eq!(roots.len(), 4);
        assert!((roots[0] - 1.0).abs() < 1e-6);
        assert!((roots[1] - 2.0).abs() < 1e-6);
        assert!((roots[2] - 3.0).abs() < 1e-6);
        assert!((roots[3] - 4.0).abs() < 1e-6);
    }

    #[test]
    fn quartic_two_real_roots() {
        // (x^2 - 1)(x^2 + 1) = x^4 - 1
        let roots = solve_quartic(1.0, 0.0, 0.0, 0.0, -1.0);
        assert_eq!(roots.len(), 2);
        assert!((roots[0] + 1.0).abs() < 1e-6);
        assert!((roots[1] - 1.0).abs() < 1e-6);
    }
}
