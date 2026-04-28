//! ShapeConstruct-style low-level shape construction utilities.
//!
//! This module provides low-level construction functions for geometric primitives
//! and topological shapes. Analogous to OCCT's `ShapeConstruct` package.
//!
//! # Modules
//!
//! | Function | Description | OCCT equivalent |
//! |---|---|---|
//! | `construct_line` | Construct a line from two points | `GC_MakeLine` |
//! | `construct_circle_from_3_points` | Circle through 3 points | `GC_MakeCircle` |
//! | `construct_circle_center_normal` | Circle from center and normal | `GC_MakeCircle` |
//! | `construct_ellipse_from_points` | Ellipse from 3 points | `GC_MakeEllipse` |
//! | `construct_plane_from_3_points` | Plane through 3 points | `GC_MakePlane` |
//! | `construct_plane_from_point_normal` | Plane from point and normal | `GC_MakePlane` |
//! | `construct_cylinder_from_axis` | Cylindrical surface | `GC_MakeCylindricalSurface` |
//! | `construct_cone_from_axis` | Conical surface | `GC_MakeConicalSurface` |
//! | `construct_sphere_from_center_radius` | Spherical surface | `GC_MakeSphere` |
//! | `construct_torus_from_center_radii` | Toroidal surface | `GC_MakeTorus` |
//! | `construct_bspline_curve` | BSpline curve from control points | `GeomAPI_PointsToBSpline` |
//! | `construct_bspline_surface` | BSpline surface from control grid | `GeomAPI_PointsToBSplineSurface` |
//! | `construct_polygon_wire` | Wire from polygon points | `BRepBuilderAPI_MakePolygon` |
//! | `construct_circle_wire` | Wire from circle discretization | `BRepBuilderAPI_MakePolygon` |
//! | `construct_planar_face_from_wire` | Planar face from wire | `BRepBuilderAPI_MakeFace` |
//! | `construct_face_from_boundary` | Face with inner wires | `BRepBuilderAPI_MakeFace` |

use glam::DVec3;
use rcad_kernel::geom::{
    BSplineCurve3, BSplineSurface, Circle3, ConicalSurface, Curve3, CylindricalSurface,
    Ellipse3, Line3, Plane, SphericalSurface, Surface3, ToroidalSurface,
};
use rcad_kernel::topology::{Edge, Face, Wire, WireEdge};

use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG};

// =============================================================================
// Curve Construction
// =============================================================================

/// Construct a line from two points.
///
/// Returns `None` if the points are coincident.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_algorithms::shape_construct::construct_line;
///
/// let p1 = DVec3::ZERO;
/// let p2 = DVec3::X;
/// let line = construct_line(p1, p2).unwrap();
/// assert!((line.origin - p1).length() < 1e-10);
/// ```
pub fn construct_line(p1: DVec3, p2: DVec3) -> Option<Line3> {
    let delta = p2 - p1;
    let len = delta.length();
    if len < TOLERANCE_ABS {
        return None;
    }
    Some(Line3 {
        origin: p1,
        direction: delta / len,
    })
}

/// Construct a circle from three points.
///
/// The three points must not be collinear and must be distinct.
/// Returns `None` if the construction fails.
///
/// # Algorithm
/// The circle center is the intersection of the perpendicular bisectors
/// of the chords (p1, p2) and (p2, p3).
pub fn construct_circle_from_3_points(p1: DVec3, p2: DVec3, p3: DVec3) -> Option<Circle3> {
    // Check for coincident points
    if (p1 - p2).length() < TOLERANCE_ABS || (p2 - p3).length() < TOLERANCE_ABS || (p1 - p3).length() < TOLERANCE_ABS {
        return None;
    }

    // Compute vectors
    let v12 = p2 - p1;
    let v23 = p3 - p2;
    let v13 = p3 - p1;

    // Compute the plane normal
    let normal = v12.cross(v13);
    let normal_len_sq = normal.length_squared();
    if normal_len_sq < TOLERANCE_ANG * TOLERANCE_ANG {
        // Points are collinear
        return None;
    }
    let normal = normal / normal_len_sq.sqrt();

    // Find the circle center using perpendicular bisector intersection
    // Midpoints
    let m12 = (p1 + p2) * 0.5;
    let m23 = (p2 + p3) * 0.5;

    // Perpendicular bisector directions (in the plane)
    let d12 = v12.cross(normal).normalize_or(DVec3::ZERO);
    let d23 = v23.cross(normal).normalize_or(DVec3::ZERO);

    // Solve for intersection: m12 + t * d12 = m23 + s * d23
    // Rearranged: t * d12 - s * d23 = m23 - m12
    // In matrix form: [d12, -d23] * [t, s]^T = m23 - m12
    let rhs = m23 - m12;
    let det = d12.x * (-d23.y) - d12.y * (-d23.x);

    // If determinant is too small, use different approach
    let center = if det.abs() > TOLERANCE_ANG {
        let t = (rhs.x * (-d23.y) - rhs.y * (-d23.x)) / det;
        m12 + d12 * t
    } else {
        // Try solving in different coordinate plane
        let det = d12.y * (-d23.z) - d12.z * (-d23.y);
        if det.abs() > TOLERANCE_ANG {
            let t = (rhs.y * (-d23.z) - rhs.z * (-d23.y)) / det;
            m12 + d12 * t
        } else {
            let det = d12.z * (-d23.x) - d12.x * (-d23.z);
            if det.abs() > TOLERANCE_ANG {
                let t = (rhs.z * (-d23.x) - rhs.x * (-d23.z)) / det;
                m12 + d12 * t
            } else {
                // Fall back to direct formula
                // Using the circumcenter formula
                let a_sq = v12.length_squared();
                let b_sq = v23.length_squared();
                let c_sq = v13.length_squared();
                let denom = 2.0 * (a_sq * b_sq + b_sq * c_sq + c_sq * a_sq
                    - a_sq * a_sq - b_sq * b_sq - c_sq * c_sq);
                if denom.abs() < TOLERANCE_ABS {
                    return None;
                }
                let alpha = a_sq * (b_sq + c_sq - a_sq) / denom;
                let beta = b_sq * (c_sq + a_sq - b_sq) / denom;
                let gamma = 1.0 - alpha - beta;
                alpha * p1 + beta * p2 + gamma * p3
            }
        }
    };

    let radius = (center - p1).length();
    if radius < TOLERANCE_ABS {
        return None;
    }

    Some(Circle3 {
        center,
        normal,
        radius,
    })
}

/// Construct a circle from center, normal, and radius.
///
/// The normal vector will be normalized. The radius must be positive.
pub fn construct_circle_center_normal(center: DVec3, normal: DVec3, radius: f64) -> Circle3 {
    Circle3 {
        center,
        normal: normal.normalize_or(DVec3::Z),
        radius: radius.abs().max(TOLERANCE_ABS),
    }
}

/// Construct an ellipse from three points.
///
/// The three points define the ellipse with:
/// - p1 and p2: endpoints of the major axis
/// - p3: a point on the ellipse (defines minor axis length)
///
/// Returns `None` if the construction fails (e.g., collinear points).
pub fn construct_ellipse_from_points(p1: DVec3, p2: DVec3, p3: DVec3) -> Option<Ellipse3> {
    // Check for coincident points
    if (p1 - p2).length() < TOLERANCE_ABS {
        return None;
    }

    // Center is midpoint of major axis
    let center = (p1 + p2) * 0.5;

    // Major axis direction and length
    let major_axis = p2 - p1;
    let major_radius = major_axis.length() * 0.5;
    let major_dir = major_axis.normalize_or(DVec3::X);

    // Compute plane normal from major axis and third point
    let v3 = p3 - center;
    let normal = major_dir.cross(v3);
    let normal_len_sq = normal.length_squared();
    if normal_len_sq < TOLERANCE_ANG * TOLERANCE_ANG {
        // Third point is collinear with major axis
        return None;
    }
    let normal = normal / normal_len_sq.sqrt();

    // Project third point onto the plane perpendicular to major axis
    let v3_proj = v3 - v3.dot(major_dir) * major_dir;
    let minor_radius = v3_proj.length();

    if minor_radius < TOLERANCE_ABS {
        // Degenerate ellipse (would be a line)
        return None;
    }

    Some(Ellipse3 {
        center,
        normal,
        major_dir,
        major_radius: major_radius.max(TOLERANCE_ABS),
        minor_radius: minor_radius.max(TOLERANCE_ABS).min(major_radius),
    })
}

// =============================================================================
// Surface Construction
// =============================================================================

/// Construct a plane from three points.
///
/// Returns `None` if the points are collinear.
///
/// # Algorithm
/// The plane origin is p1, and the normal is computed from the cross product
/// of vectors (p2-p1) and (p3-p1).
pub fn construct_plane_from_3_points(p1: DVec3, p2: DVec3, p3: DVec3) -> Option<Plane> {
    let v1 = p2 - p1;
    let v2 = p3 - p1;
    let normal = v1.cross(v2);
    let len = normal.length();
    if len < TOLERANCE_ABS {
        return None;
    }
    Some(Plane {
        origin: p1,
        normal: normal / len,
    })
}

/// Construct a plane from a point and normal vector.
///
/// The normal will be normalized.
pub fn construct_plane_from_point_normal(point: DVec3, normal: DVec3) -> Plane {
    Plane {
        origin: point,
        normal: normal.normalize_or(DVec3::Z),
    }
}

/// Construct a cylindrical surface from axis point, direction, and radius.
///
/// # Arguments
/// * `axis_point` - A point on the cylinder axis
/// * `axis_dir` - Direction of the cylinder axis (will be normalized)
/// * `radius` - Radius of the cylinder (must be positive)
pub fn construct_cylinder_from_axis(axis_point: DVec3, axis_dir: DVec3, radius: f64) -> CylindricalSurface {
    CylindricalSurface {
        origin: axis_point,
        axis: axis_dir.normalize_or(DVec3::Z),
        radius: radius.abs().max(TOLERANCE_ABS),
    }
}

/// Construct a conical surface from apex, axis direction, and half-angle.
///
/// # Arguments
/// * `apex` - The apex point of the cone
/// * `axis_dir` - Direction of the cone axis (will be normalized)
/// * `angle` - Half-angle of the cone in radians (angle between axis and surface)
pub fn construct_cone_from_axis(apex: DVec3, axis_dir: DVec3, angle: f64) -> ConicalSurface {
    ConicalSurface {
        apex,
        axis: axis_dir.normalize_or(DVec3::Z),
        radius: 0.0,
        half_angle_rad: angle.abs().clamp(TOLERANCE_ANG, std::f64::consts::FRAC_PI_2 - TOLERANCE_ANG),
    }
}

/// Construct a spherical surface from center and radius.
///
/// # Arguments
/// * `center` - Center of the sphere
/// * `radius` - Radius of the sphere (must be positive)
pub fn construct_sphere_from_center_radius(center: DVec3, radius: f64) -> SphericalSurface {
    SphericalSurface {
        center,
        axis: DVec3::Z, // Default axis
        radius: radius.abs().max(TOLERANCE_ABS),
    }
}

/// Construct a toroidal surface from center, axis, and radii.
///
/// # Arguments
/// * `center` - Center of the torus
/// * `axis` - Axis direction (will be normalized)
/// * `major` - Major radius (distance from center to tube center)
/// * `minor` - Minor radius (tube radius)
pub fn construct_torus_from_center_radii(
    center: DVec3,
    axis: DVec3,
    major: f64,
    minor: f64,
) -> ToroidalSurface {
    ToroidalSurface {
        center,
        axis: axis.normalize_or(DVec3::Z),
        major_radius: major.abs().max(TOLERANCE_ABS),
        minor_radius: minor.abs().max(TOLERANCE_ABS),
    }
}

// =============================================================================
// BSpline Construction
// =============================================================================

/// Construct a BSpline curve from control points, knots, and degree.
///
/// # Arguments
/// * `control_points` - Control point positions
/// * `knots` - Knot vector (must have `control_points.len() + degree + 1` knots)
/// * `degree` - Degree of the BSpline curve
///
/// # Panics
/// Does not panic; returns a default curve if inputs are invalid.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_algorithms::shape_construct::construct_bspline_curve;
///
/// let pts = vec![
///     DVec3::ZERO,
///     DVec3::new(1.0, 1.0, 0.0),
///     DVec3::new(2.0, 0.0, 0.0),
/// ];
/// let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
/// let curve = construct_bspline_curve(&pts, &knots, 2);
/// ```
pub fn construct_bspline_curve(control_points: &[DVec3], knots: &[f64], degree: usize) -> BSplineCurve3 {
    let n_cp = control_points.len();
    let expected_knots = n_cp + degree + 1;

    // Validate inputs
    if n_cp < 2 || degree == 0 || knots.len() != expected_knots {
        // Return a degenerate curve
        return BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
        };
    }

    // Normalize knots to [0, 1]
    let t_min = knots[0];
    let t_max = knots[knots.len() - 1];
    let range = t_max - t_min;
    let normalized_knots: Vec<f64> = if range > TOLERANCE_ABS {
        knots.iter().map(|&t| (t - t_min) / range).collect()
    } else {
        knots.to_vec()
    };

    // Default weights (non-rational curve)
    let weights = vec![1.0; n_cp];

    BSplineCurve3 {
        degree,
        knots: normalized_knots,
        control_points: control_points.to_vec(),
        weights,
    }
}

/// Construct a BSpline surface from control point grid, knots, and degrees.
///
/// # Arguments
/// * `control_points` - Control point grid [u_index][v_index]
/// * `u_knots` - Knot vector for U direction
/// * `v_knots` - Knot vector for V direction
/// * `u_deg` - Degree in U direction
/// * `v_deg` - Degree in V direction
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_algorithms::shape_construct::construct_bspline_surface;
///
/// let pts = vec![
///     vec![DVec3::ZERO, DVec3::Y],
///     vec![DVec3::X, DVec3::new(1.0, 1.0, 0.0)],
/// ];
/// let u_knots = vec![0.0, 0.0, 1.0, 1.0];
/// let v_knots = vec![0.0, 0.0, 1.0, 1.0];
/// let surf = construct_bspline_surface(&pts, &u_knots, &v_knots, 1, 1);
/// ```
pub fn construct_bspline_surface(
    control_points: &[Vec<DVec3>],
    u_knots: &[f64],
    v_knots: &[f64],
    u_deg: usize,
    v_deg: usize,
) -> BSplineSurface {
    let n_u = control_points.len();
    let n_v = control_points.first().map(|row| row.len()).unwrap_or(0);

    // Validate inputs
    let expected_u_knots = n_u + u_deg + 1;
    let expected_v_knots = n_v + v_deg + 1;

    if n_u < 2 || n_v < 2 || u_deg == 0 || v_deg == 0
        || u_knots.len() != expected_u_knots || v_knots.len() != expected_v_knots
        || control_points.iter().any(|row| row.len() != n_v)
    {
        // Return a degenerate surface
        return BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![vec![DVec3::ZERO, DVec3::Y], vec![DVec3::X, DVec3::X + DVec3::Y]],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        };
    }

    // Normalize knots
    let u_min = u_knots[0];
    let u_max = u_knots[u_knots.len() - 1];
    let v_min = v_knots[0];
    let v_max = v_knots[v_knots.len() - 1];

    let u_range = u_max - u_min;
    let v_range = v_max - v_min;

    let normalized_u_knots: Vec<f64> = if u_range > TOLERANCE_ABS {
        u_knots.iter().map(|&t| (t - u_min) / u_range).collect()
    } else {
        u_knots.to_vec()
    };

    let normalized_v_knots: Vec<f64> = if v_range > TOLERANCE_ABS {
        v_knots.iter().map(|&t| (t - v_min) / v_range).collect()
    } else {
        v_knots.to_vec()
    };

    // Default weights (non-rational surface)
    let weights: Vec<Vec<f64>> = control_points.iter().map(|row| vec![1.0; row.len()]).collect();

    BSplineSurface {
        degree_u: u_deg,
        degree_v: v_deg,
        knots_u: normalized_u_knots,
        knots_v: normalized_v_knots,
        control_points: control_points.to_vec(),
        weights,
    }
}

// =============================================================================
// Wire Construction
// =============================================================================

/// Construct a wire from a polygon of points.
///
/// # Arguments
/// * `points` - Ordered points defining the polygon vertices
/// * `closed` - If true, adds an edge from the last point back to the first
///
/// # Returns
/// A tuple `(vertices, edges, wire)` where:
/// - `vertices` are the vertex positions
/// - `edges` connect consecutive vertices (and last to first if closed)
/// - `wire` references the edges
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_algorithms::shape_construct::construct_polygon_wire;
///
/// let pts = vec![DVec3::ZERO, DVec3::X, DVec3::X + DVec3::Y, DVec3::Y];
/// let (vertices, edges, wire) = construct_polygon_wire(&pts, true);
/// assert_eq!(vertices.len(), 4);
/// assert_eq!(edges.len(), 4); // Closed: 4 edges
/// ```
pub fn construct_polygon_wire(points: &[DVec3], closed: bool) -> (Vec<DVec3>, Vec<Edge>, Wire) {
    if points.is_empty() {
        return (vec![], vec![], Wire { edges: vec![] });
    }

    let n = points.len();
    let edge_count = if closed && n >= 3 { n } else { n.saturating_sub(1) };

    // Build edges
    let mut edges = Vec::with_capacity(edge_count);
    for i in 0..n.saturating_sub(1) {
        edges.push(Edge { start: i, end: i + 1 });
    }
    if closed && n >= 3 {
        edges.push(Edge { start: n - 1, end: 0 });
    }

    // Build wire with forward-oriented edges
    let wire = Wire {
        edges: edges.iter().enumerate().map(|(i, _)| WireEdge::fwd(i)).collect(),
    };

    (points.to_vec(), edges, wire)
}

/// Construct a wire approximating a circle by polygonal segments.
///
/// # Arguments
/// * `center` - Center of the circle
/// * `normal` - Normal to the plane of the circle (will be normalized)
/// * `radius` - Radius of the circle
/// * `n_segments` - Number of line segments to approximate the circle
///
/// # Returns
/// A tuple `(vertices, edges, wire)` representing the discretized circle.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_algorithms::shape_construct::construct_circle_wire;
///
/// let (vertices, edges, wire) = construct_circle_wire(DVec3::ZERO, DVec3::Z, 1.0, 8);
/// assert_eq!(vertices.len(), 8);
/// assert_eq!(edges.len(), 8);
/// ```
pub fn construct_circle_wire(
    center: DVec3,
    normal: DVec3,
    radius: f64,
    n_segments: usize,
) -> (Vec<DVec3>, Vec<Edge>, Wire) {
    if n_segments < 3 {
        return (vec![], vec![], Wire { edges: vec![] });
    }

    let normal = normal.normalize_or(DVec3::Z);
    let radius = radius.abs().max(TOLERANCE_ABS);

    // Build orthonormal basis in the plane
    let x_dir = any_perpendicular(normal);
    let y_dir = normal.cross(x_dir).normalize();

    // Generate vertices
    let mut vertices = Vec::with_capacity(n_segments);
    for i in 0..n_segments {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n_segments as f64);
        let x = radius * angle.cos();
        let y = radius * angle.sin();
        let point = center + x * x_dir + y * y_dir;
        vertices.push(point);
    }

    // Build edges (closed loop)
    let mut edges = Vec::with_capacity(n_segments);
    for i in 0..n_segments {
        edges.push(Edge {
            start: i,
            end: (i + 1) % n_segments,
        });
    }

    // Build wire
    let wire = Wire {
        edges: edges.iter().enumerate().map(|(i, _)| WireEdge::fwd(i)).collect(),
    };

    (vertices, edges, wire)
}

// =============================================================================
// Face Construction
// =============================================================================

/// Construct a planar face from a wire and surface.
///
/// # Arguments
/// * `wire` - The boundary wire
/// * `surface` - The underlying surface (should be a plane)
///
/// # Returns
/// A face with the given wire as outer boundary.
pub fn construct_planar_face_from_wire(wire: &Wire, surface: &Surface3) -> Face {
    let normal = match surface {
        Surface3::Plane(p) => p.normal,
        _ => DVec3::Z,
    };

    Face {
        outer_wire: wire.clone(),
        inner_wires: vec![],
        normal,
        triangles: vec![],
        mesh_dirty: true,
    }
}

/// Construct a face from boundary wires and surface.
///
/// # Arguments
/// * `outer_wire` - The outer boundary wire
/// * `inner_wires` - Inner boundary wires (holes)
/// * `surface` - The underlying surface
///
/// # Returns
/// A face with outer and inner boundaries.
///
/// # Example
/// ```rust
/// use glam::DVec3;
/// use rcad_algorithms::shape_construct::{construct_polygon_wire, construct_face_from_boundary};
/// use rcad_kernel::geom::{Surface3, Plane};
///
/// let (verts, edges, outer) = construct_polygon_wire(
///     &[DVec3::ZERO, DVec3::new(10.0, 0.0, 0.0), DVec3::new(10.0, 10.0, 0.0), DVec3::new(0.0, 10.0, 0.0)],
///     true
/// );
/// let surface = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
/// let face = construct_face_from_boundary(outer, vec![], surface);
/// assert!(face.inner_wires.is_empty());
/// ```
pub fn construct_face_from_boundary(outer_wire: Wire, inner_wires: Vec<Wire>, surface: Surface3) -> Face {
    let normal = match &surface {
        Surface3::Plane(p) => p.normal,
        Surface3::Cylinder(c) => {
            // Use cylinder axis as reference
            c.axis
        }
        Surface3::Sphere(s) => {
            // Use Z as default for spherical faces
            s.axis
        }
        Surface3::Cone(c) => c.axis,
        Surface3::Torus(t) => t.axis,
        _ => DVec3::Z,
    };

    Face {
        outer_wire,
        inner_wires,
        normal,
        triangles: vec![],
        mesh_dirty: true,
    }
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Return any unit vector perpendicular to the given vector.
fn any_perpendicular(v: DVec3) -> DVec3 {
    let v = v.normalize_or(DVec3::Z);
    let perp = if v.x.abs() > 0.5 {
        DVec3::new(-v.y, v.x, 0.0)
    } else {
        DVec3::new(0.0, -v.z, v.y)
    };
    perp.normalize_or(DVec3::X)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{CurveEval, SurfaceEval};

    // -------------------------------------------------------------------------
    // Curve Construction Tests
    // -------------------------------------------------------------------------

    #[test]
    fn construct_line_basic() {
        let p1 = DVec3::ZERO;
        let p2 = DVec3::X;
        let line = construct_line(p1, p2).unwrap();

        assert!((line.origin - p1).length() < TOLERANCE_ABS);
        assert!((line.direction - DVec3::X).length() < TOLERANCE_ANG);
    }

    #[test]
    fn construct_line_coincident_points() {
        let p = DVec3::new(1.0, 2.0, 3.0);
        assert!(construct_line(p, p).is_none());
    }

    #[test]
    fn construct_circle_from_3_points_basic() {
        // Unit circle in XY plane
        let p1 = DVec3::new(1.0, 0.0, 0.0);
        let p2 = DVec3::new(0.0, 1.0, 0.0);
        let p3 = DVec3::new(-1.0, 0.0, 0.0);

        let circle = construct_circle_from_3_points(p1, p2, p3).unwrap();

        assert!((circle.center - DVec3::ZERO).length() < TOLERANCE_ABS);
        assert!((circle.radius - 1.0).abs() < 1e-6);
        assert!((circle.normal - DVec3::Z).length() < TOLERANCE_ANG
            || (circle.normal + DVec3::Z).length() < TOLERANCE_ANG);
    }

    #[test]
    fn construct_circle_from_3_points_collinear() {
        let p1 = DVec3::ZERO;
        let p2 = DVec3::X;
        let p3 = DVec3::new(2.0, 0.0, 0.0);

        assert!(construct_circle_from_3_points(p1, p2, p3).is_none());
    }

    #[test]
    fn construct_circle_from_3_points_coincident() {
        let p1 = DVec3::ZERO;
        let p2 = DVec3::X;
        let p3 = DVec3::X;

        assert!(construct_circle_from_3_points(p1, p2, p3).is_none());
    }

    #[test]
    fn construct_circle_center_normal_basic() {
        let center = DVec3::new(1.0, 2.0, 3.0);
        let normal = DVec3::Z;
        let radius = 5.0;

        let circle = construct_circle_center_normal(center, normal, radius);

        assert!((circle.center - center).length() < TOLERANCE_ABS);
        assert!((circle.radius - radius).abs() < TOLERANCE_ABS);
        assert!((circle.normal - DVec3::Z).length() < TOLERANCE_ANG);
    }

    #[test]
    fn construct_ellipse_from_points_basic() {
        // Ellipse with major axis along X, minor along Y
        let p1 = DVec3::new(-2.0, 0.0, 0.0);
        let p2 = DVec3::new(2.0, 0.0, 0.0);
        let p3 = DVec3::new(0.0, 1.0, 0.0);

        let ellipse = construct_ellipse_from_points(p1, p2, p3).unwrap();

        assert!((ellipse.center - DVec3::ZERO).length() < TOLERANCE_ABS);
        assert!((ellipse.major_radius - 2.0).abs() < 1e-6);
        assert!((ellipse.minor_radius - 1.0).abs() < 1e-6);
    }

    #[test]
    fn construct_ellipse_from_points_coincident() {
        let p1 = DVec3::ZERO;
        let p2 = DVec3::ZERO;
        let p3 = DVec3::Y;

        assert!(construct_ellipse_from_points(p1, p2, p3).is_none());
    }

    // -------------------------------------------------------------------------
    // Surface Construction Tests
    // -------------------------------------------------------------------------

    #[test]
    fn construct_plane_from_3_points_basic() {
        let p1 = DVec3::ZERO;
        let p2 = DVec3::X;
        let p3 = DVec3::Y;

        let plane = construct_plane_from_3_points(p1, p2, p3).unwrap();

        assert!((plane.origin - p1).length() < TOLERANCE_ABS);
        assert!((plane.normal - DVec3::Z).length() < TOLERANCE_ANG
            || (plane.normal + DVec3::Z).length() < TOLERANCE_ANG);
    }

    #[test]
    fn construct_plane_from_3_points_collinear() {
        let p1 = DVec3::ZERO;
        let p2 = DVec3::X;
        let p3 = DVec3::new(2.0, 0.0, 0.0);

        assert!(construct_plane_from_3_points(p1, p2, p3).is_none());
    }

    #[test]
    fn construct_plane_from_point_normal_basic() {
        let point = DVec3::new(1.0, 2.0, 3.0);
        let normal = DVec3::new(1.0, 1.0, 0.0);

        let plane = construct_plane_from_point_normal(point, normal);

        assert!((plane.origin - point).length() < TOLERANCE_ABS);
        let expected_normal = normal.normalize();
        assert!((plane.normal - expected_normal).length() < TOLERANCE_ANG);
    }

    #[test]
    fn construct_cylinder_from_axis_basic() {
        let axis_point = DVec3::new(1.0, 2.0, 3.0);
        let axis_dir = DVec3::Z;
        let radius = 2.5;

        let cyl = construct_cylinder_from_axis(axis_point, axis_dir, radius);

        assert!((cyl.origin - axis_point).length() < TOLERANCE_ABS);
        assert!((cyl.radius - radius).abs() < TOLERANCE_ABS);
        assert!((cyl.axis - DVec3::Z).length() < TOLERANCE_ANG);
    }

    #[test]
    fn construct_cone_from_axis_basic() {
        let apex = DVec3::ZERO;
        let axis_dir = DVec3::Y;
        let angle = std::f64::consts::FRAC_PI_4;

        let cone = construct_cone_from_axis(apex, axis_dir, angle);

        assert!((cone.apex - apex).length() < TOLERANCE_ABS);
        assert!((cone.half_angle_rad - angle).abs() < TOLERANCE_ANG);
        assert!((cone.axis - DVec3::Y).length() < TOLERANCE_ANG);
    }

    #[test]
    fn construct_sphere_from_center_radius_basic() {
        let center = DVec3::new(1.0, 2.0, 3.0);
        let radius = 5.0;

        let sphere = construct_sphere_from_center_radius(center, radius);

        assert!((sphere.center - center).length() < TOLERANCE_ABS);
        assert!((sphere.radius - radius).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn construct_torus_from_center_radii_basic() {
        let center = DVec3::ZERO;
        let axis = DVec3::Z;
        let major = 5.0;
        let minor = 1.0;

        let torus = construct_torus_from_center_radii(center, axis, major, minor);

        assert!((torus.center - center).length() < TOLERANCE_ABS);
        assert!((torus.major_radius - major).abs() < TOLERANCE_ABS);
        assert!((torus.minor_radius - minor).abs() < TOLERANCE_ABS);
        assert!((torus.axis - DVec3::Z).length() < TOLERANCE_ANG);
    }

    // -------------------------------------------------------------------------
    // BSpline Construction Tests
    // -------------------------------------------------------------------------

    #[test]
    fn construct_bspline_curve_quadratic() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(0.5, 1.0, 0.0),
            DVec3::X,
        ];
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let curve = construct_bspline_curve(&pts, &knots, 2);

        assert_eq!(curve.degree, 2);
        assert_eq!(curve.control_points.len(), 3);
        assert_eq!(curve.knots.len(), 6);

        // Curve should start at first control point
        let p_start = Curve3::BSpline(curve.clone()).point_at(0.0);
        assert!((p_start - pts[0]).length() < TOLERANCE_ABS);

        // Curve should end at last control point
        let p_end = Curve3::BSpline(curve).point_at(1.0);
        assert!((p_end - pts[2]).length() < TOLERANCE_ABS);
    }

    #[test]
    fn construct_bspline_curve_linear() {
        let pts = vec![DVec3::ZERO, DVec3::X];
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        let curve = construct_bspline_curve(&pts, &knots, 1);

        assert_eq!(curve.degree, 1);
        assert_eq!(curve.control_points.len(), 2);

        // Linear B-spline should be a straight line
        let curve_geom = Curve3::BSpline(curve);
        let p_mid = curve_geom.point_at(0.5);
        assert!((p_mid - DVec3::new(0.5, 0.0, 0.0)).length() < 1e-6);
    }

    #[test]
    fn construct_bspline_curve_invalid_inputs() {
        // Too few control points
        let pts = vec![DVec3::ZERO];
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        let curve = construct_bspline_curve(&pts, &knots, 1);
        // Should return default curve
        assert_eq!(curve.degree, 1);
        assert_eq!(curve.control_points.len(), 2);
    }

    #[test]
    fn construct_bspline_surface_bilinear() {
        let pts = vec![
            vec![DVec3::ZERO, DVec3::Y],
            vec![DVec3::X, DVec3::new(1.0, 1.0, 0.0)],
        ];
        let u_knots = vec![0.0, 0.0, 1.0, 1.0];
        let v_knots = vec![0.0, 0.0, 1.0, 1.0];
        let surf = construct_bspline_surface(&pts, &u_knots, &v_knots, 1, 1);

        assert_eq!(surf.degree_u, 1);
        assert_eq!(surf.degree_v, 1);
        assert_eq!(surf.control_points.len(), 2);

        // Check corners
        let surf_geom = Surface3::BSpline(surf);
        let p00 = surf_geom.point_at(0.0, 0.0);
        let p11 = surf_geom.point_at(1.0, 1.0);

        assert!((p00 - DVec3::ZERO).length() < TOLERANCE_ABS);
        assert!((p11 - DVec3::new(1.0, 1.0, 0.0)).length() < TOLERANCE_ABS);
    }

    // -------------------------------------------------------------------------
    // Wire Construction Tests
    // -------------------------------------------------------------------------

    #[test]
    fn construct_polygon_wire_open() {
        let pts = vec![DVec3::ZERO, DVec3::X, DVec3::X + DVec3::Y];
        let (vertices, edges, wire) = construct_polygon_wire(&pts, false);

        assert_eq!(vertices.len(), 3);
        assert_eq!(edges.len(), 2); // Open: 2 edges
        assert_eq!(wire.edges.len(), 2);
    }

    #[test]
    fn construct_polygon_wire_closed() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::X,
            DVec3::X + DVec3::Y,
            DVec3::Y,
        ];
        let (vertices, edges, wire) = construct_polygon_wire(&pts, true);

        assert_eq!(vertices.len(), 4);
        assert_eq!(edges.len(), 4); // Closed: 4 edges
        assert_eq!(wire.edges.len(), 4);

        // Check that last edge closes the loop
        let last_edge = &edges[3];
        assert_eq!(last_edge.start, 3);
        assert_eq!(last_edge.end, 0);
    }

    #[test]
    fn construct_polygon_wire_empty() {
        let (vertices, edges, wire) = construct_polygon_wire(&[], false);

        assert!(vertices.is_empty());
        assert!(edges.is_empty());
        assert!(wire.edges.is_empty());
    }

    #[test]
    fn construct_circle_wire_basic() {
        let (vertices, edges, wire) = construct_circle_wire(DVec3::ZERO, DVec3::Z, 1.0, 8);

        assert_eq!(vertices.len(), 8);
        assert_eq!(edges.len(), 8);
        assert_eq!(wire.edges.len(), 8);

        // Check that vertices lie on a circle
        for v in &vertices {
            let r = (v.x * v.x + v.y * v.y).sqrt();
            assert!((r - 1.0).abs() < 1e-10);
        }

        // Check that wire is closed
        let first_edge = &edges[0];
        let last_edge = &edges[7];
        assert_eq!(first_edge.start, 0);
        assert_eq!(last_edge.end, 0);
    }

    #[test]
    fn construct_circle_wire_too_few_segments() {
        let (vertices, edges, wire) = construct_circle_wire(DVec3::ZERO, DVec3::Z, 1.0, 2);

        assert!(vertices.is_empty());
        assert!(edges.is_empty());
    }

    #[test]
    fn construct_circle_wire_offset_center() {
        let center = DVec3::new(5.0, 5.0, 5.0);
        let (vertices, _, _) = construct_circle_wire(center, DVec3::Z, 2.0, 16);

        // Check that vertices are centered correctly
        let avg: DVec3 = vertices.iter().sum::<DVec3>() / vertices.len() as f64;
        assert!((avg - center).length() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // Face Construction Tests
    // -------------------------------------------------------------------------

    #[test]
    fn construct_planar_face_from_wire_basic() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::X,
            DVec3::X + DVec3::Y,
            DVec3::Y,
        ];
        let (_, edges, wire) = construct_polygon_wire(&pts, true);
        let surface = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let face = construct_planar_face_from_wire(&wire, &surface);

        assert_eq!(face.outer_wire.edges.len(), 4);
        assert!((face.normal - DVec3::Z).length() < TOLERANCE_ANG);
        assert!(face.inner_wires.is_empty());
    }

    #[test]
    fn construct_face_from_boundary_basic() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(10.0, 0.0, 0.0),
            DVec3::new(10.0, 10.0, 0.0),
            DVec3::new(0.0, 10.0, 0.0),
        ];
        let (_, _, outer) = construct_polygon_wire(&pts, true);

        let surface = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let face = construct_face_from_boundary(outer.clone(), vec![], surface);

        assert_eq!(face.outer_wire.edges.len(), 4);
        assert!(face.inner_wires.is_empty());
    }

    #[test]
    fn construct_face_from_boundary_with_hole() {
        let outer_pts = vec![
            DVec3::ZERO,
            DVec3::new(10.0, 0.0, 0.0),
            DVec3::new(10.0, 10.0, 0.0),
            DVec3::new(0.0, 10.0, 0.0),
        ];
        let (_, _, outer) = construct_polygon_wire(&outer_pts, true);

        let inner_pts = vec![
            DVec3::new(3.0, 3.0, 0.0),
            DVec3::new(7.0, 3.0, 0.0),
            DVec3::new(7.0, 7.0, 0.0),
            DVec3::new(3.0, 7.0, 0.0),
        ];
        let (_, _, inner) = construct_polygon_wire(&inner_pts, true);

        let surface = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let face = construct_face_from_boundary(outer, vec![inner], surface);

        assert_eq!(face.outer_wire.edges.len(), 4);
        assert_eq!(face.inner_wires.len(), 1);
        assert_eq!(face.inner_wires[0].edges.len(), 4);
    }

    // -------------------------------------------------------------------------
    // Helper Function Tests
    // -------------------------------------------------------------------------

    #[test]
    fn any_perpendicular_basic() {
        let v = DVec3::Z;
        let perp = any_perpendicular(v);

        // Should be perpendicular
        assert!((v.cross(perp)).length() > 0.9);
        assert!(v.dot(perp).abs() < TOLERANCE_ANG);
    }

    #[test]
    fn any_perpendicular_x_axis() {
        let v = DVec3::X;
        let perp = any_perpendicular(v);

        assert!(v.dot(perp).abs() < TOLERANCE_ANG);
        assert!((perp.length() - 1.0).abs() < TOLERANCE_ANG);
    }

    #[test]
    fn any_perpendicular_arbitrary() {
        let v = DVec3::new(1.0, 2.0, 3.0).normalize();
        let perp = any_perpendicular(v);

        assert!(v.dot(perp).abs() < TOLERANCE_ANG);
        assert!((perp.length() - 1.0).abs() < TOLERANCE_ANG);
    }
}
