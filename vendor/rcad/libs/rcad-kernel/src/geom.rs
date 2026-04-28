use glam::{DVec2, DVec3};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

pub type Point3 = DVec3;
pub type Vec3 = DVec3;
pub type Point2 = DVec2;
pub type Vec2 = DVec2;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Line3 {
    pub origin: Point3,
    pub direction: Vec3,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Circle3 {
    pub center: Point3,
    pub normal: Vec3,
    pub radius: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Ellipse3 {
    pub center: Point3,
    pub normal: Vec3,
    pub major_dir: Vec3,
    pub major_radius: f64,
    pub minor_radius: f64,
}

/// A non-uniform rational B-spline curve in 3D.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BSplineCurve3 {
    pub degree: usize,
    /// Full knot vector (with multiplicities expanded).
    pub knots: Vec<f64>,
    pub control_points: Vec<DVec3>,
    /// Homogeneous weights; 1.0 for non-rational.
    pub weights: Vec<f64>,
}

/// A rational or non-rational Bezier curve in 3D.
///
/// Evaluated via de Casteljau's algorithm. Domain is always `[0.0, 1.0]`.
/// Analogous to OCCT `Geom_BezierCurve`.
///
/// Note: a Bezier curve of degree n is equivalent to a B-spline of degree n
/// with knot vector `[0, 鈥? 0, 1, 鈥? 1]` (n+1 times each).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BezierCurve3 {
    pub control_points: Vec<DVec3>,
    /// Homogeneous weights; 1.0 for non-rational (polynomial Bezier).
    pub weights: Vec<f64>,
}

/// A 3D hyperbola defined by center, normal, semi-transverse axis `a`, and
/// semi-conjugate axis `b`.  Parametric form:
///
///   P(t) = center + a路cosh(t)路major_dir + b路sinh(t)路minor_dir
///
/// where `minor_dir = normal 脳 major_dir`.  Domain is `(鈭掆垶, +鈭?`;
/// the principal branch (t 鈮?0) is on the `+major_dir` side.
/// Analogous to OCCT `Geom_Hyperbola`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Hyperbola3 {
    pub center: Point3,
    pub normal: Vec3,
    pub major_dir: Vec3,
    pub semi_major: f64, // a  (transverse semi-axis)
    pub semi_minor: f64, // b  (conjugate semi-axis)
}

/// A 3D parabola defined by its vertex, axis, and focal parameter `p`
/// (where the focus is at distance `p/2` from the vertex along the axis).
///
///   P(t) = vertex + (t虏/(2p))路axis_dir + t路dir_perp
///
/// where `dir_perp = normal 脳 axis_dir` is the cross-axis direction.
/// Domain is `(鈭掆垶, +鈭?`.  Analogous to OCCT `Geom_Parabola`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Parabola3 {
    pub vertex: Point3,
    pub normal: Vec3,
    pub axis_dir: Vec3,   // direction from vertex toward focus
    pub focal_param: f64, // p  (= 2 脳 focal_length)
}

/// A circular helix curve around an axis.
///
/// Parameterization:
/// `P(t) = origin + radius*(cos t * x_axis + sin t * y_axis) + (pitch/(2*pi))*t * axis`
///
/// Analogous to OCCT TKHelix circular helix primitives.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CircularHelix3 {
    pub origin: Point3,
    pub axis: Vec3,
    /// A reference direction orthogonalized against `axis` at evaluation time.
    pub ref_dir: Vec3,
    pub radius: f64,
    /// Axial advance per full revolution (2*pi in parameter).
    pub pitch: f64,
}

/// A 3D sine-wave curve traveling along a baseline direction with amplitude
/// in a perpendicular `amplitude_dir`.
///
/// Parameterization:
/// `P(t) = origin + t * baseline_dir + amplitude * sin(frequency * t + phase) * amplitude_dir`
///
/// `baseline_dir` and `amplitude_dir` should be orthogonal unit vectors.
/// Analogous to OCCT `GeomEval_SineWaveCurve`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SineWave3 {
    pub origin: Point3,
    /// Unit direction along which the parameter `t` advances.
    pub baseline_dir: Vec3,
    /// Unit direction of the sine-wave displacement (orthogonal to `baseline_dir`).
    pub amplitude_dir: Vec3,
    pub amplitude: f64,
    pub frequency: f64,
    pub phase: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Curve3 {
    Line(Line3),
    Circle(Circle3),
    Ellipse(Ellipse3),
    BSpline(BSplineCurve3),
    Bezier(BezierCurve3),  // Phase M
    Offset(OffsetCurve3),  // Phase M
    Hyperbola(Hyperbola3), // Phase S
    Parabola(Parabola3),   // Phase S
    CircularHelix(CircularHelix3), // Phase P7
    SineWave(SineWave3),   // Phase P7
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Plane {
    pub origin: Point3,
    pub normal: Vec3,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CylindricalSurface {
    pub origin: Point3,
    pub axis: Vec3,
    pub radius: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SphericalSurface {
    pub center: Point3,
    pub axis: Vec3,
    pub radius: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ConicalSurface {
    /// Point on the cone axis where the surface radius equals `radius`.
    ///
    /// Historically this field was used as an apex for zero-radius primitive
    /// cones. For general conical surfaces, the true apex is derived from this
    /// reference point, `radius`, and `half_angle_rad`.
    pub apex: Point3,
    pub axis: Vec3,
    /// Radius of the reference circle at `apex`.
    pub radius: f64,
    pub half_angle_rad: f64,
}

impl ConicalSurface {
    pub fn axis_dir(&self) -> DVec3 {
        self.axis.normalize_or_zero()
    }

    pub fn apex_point(&self) -> DVec3 {
        let tan_half = self.half_angle_rad.tan();
        if tan_half.abs() < 1e-12 {
            self.apex
        } else {
            self.apex - self.axis_dir() * (self.radius / tan_half)
        }
    }

    pub fn axial_from_slant(&self, slant: f64) -> f64 {
        slant * self.half_angle_rad.cos()
    }

    pub fn slant_from_axial(&self, axial: f64) -> f64 {
        let cos_half = self.half_angle_rad.cos();
        if cos_half.abs() < 1e-12 {
            0.0
        } else {
            axial / cos_half
        }
    }

    pub fn radius_at_slant(&self, slant: f64) -> f64 {
        self.radius + slant * self.half_angle_rad.sin()
    }

    pub fn radius_at_axial(&self, axial: f64) -> f64 {
        self.radius + axial * self.half_angle_rad.tan()
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ToroidalSurface {
    pub center: Point3,
    pub axis: Vec3,
    pub major_radius: f64,
    pub minor_radius: f64,
}

/// An ellipsoidal surface aligned to a local orthonormal frame.
///
/// Parameterization matches sphere-like angles:
/// - `u` = longitude `[0, 2π]`
/// - `v` = colatitude `[0, π]` (0 at +axis pole)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct EllipsoidalSurface {
    pub center: Point3,
    pub axis: Vec3,
    /// Reference direction used to derive the local X axis.
    pub ref_dir: Vec3,
    pub radius_x: f64,
    pub radius_y: f64,
    pub radius_z: f64,
}

/// A classical helicoid surface around an axis.
///
/// Parameterization:
/// `S(u, v) = origin + v * (cos(u) * x_axis + sin(u) * y_axis) + (pitch/(2*pi))*u * axis`
///
/// `u` is the azimuth / screw parameter and `v` is the signed radial distance.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct HelicoidSurface {
    pub origin: Point3,
    pub axis: Vec3,
    /// Reference direction used to derive the local X axis.
    pub ref_dir: Vec3,
    /// Axial advance per full revolution.
    pub pitch: f64,
}

/// A circular pipe/tube surface around a spine curve.
///
/// `u` is the azimuth angle around the local section frame and `v` follows the
/// natural parameter of the spine curve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipeSurface {
    pub spine: Box<Curve3>,
    /// Initial/reference direction projected onto the normal plane of the
    /// spine tangent at evaluation time.
    pub ref_dir: Vec3,
    pub radius: f64,
}

/// A non-uniform rational B-spline surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BSplineSurface {
    pub degree_u: usize,
    pub degree_v: usize,
    /// Full knot vector for u (with multiplicities expanded).
    pub knots_u: Vec<f64>,
    /// Full knot vector for v (with multiplicities expanded).
    pub knots_v: Vec<f64>,
    /// Control point grid [u_index][v_index].
    pub control_points: Vec<Vec<DVec3>>,
    /// Weight grid [u_index][v_index]; 1.0 for non-rational.
    pub weights: Vec<Vec<f64>>,
}

/// A rational or non-rational Bezier surface (tensor-product bicubic patch).
///
/// Evaluated by applying de Casteljau in u, then in v. Domain is `[0, 1] 脳 [0, 1]`.
/// Analogous to OCCT `Geom_BezierSurface`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BezierSurface {
    /// Control point grid [u_count][v_count].
    pub control_points: Vec<Vec<DVec3>>,
    /// Weight grid [u_count][v_count]; 1.0 for non-rational.
    pub weights: Vec<Vec<f64>>,
}

/// A triangular rational Bezier surface using barycentric coordinates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriBezierSurface {
    /// Triangular control net rows. Row `i` has `degree + 1 - i` points.
    pub control_points: Vec<Vec<DVec3>>,
    /// Weight rows with the same triangular layout as `control_points`.
    pub weights: Vec<Vec<f64>>,
}

/// A curve offset from a base curve by a fixed distance in a reference plane.
///
/// `S(t) = basis.point_at(t) + offset_distance * (tangent(t) 脳 offset_dir).normalize()`
///
/// Analogous to OCCT `Geom_OffsetCurve`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OffsetCurve3 {
    pub basis: Box<Curve3>,
    /// Offset distance (positive = outward from the curve's "left" side).
    pub offset_distance: f64,
    /// Fixed reference direction (normal to the offset plane).
    /// The offset direction at each point is `(tangent 脳 offset_dir).normalize()`.
    pub offset_dir: Vec3,
}

/// A surface offset from a base surface by a fixed distance along the normal.
///
/// `S(u,v) = basis.point_at(u,v) + offset_distance * basis.normal_at(u,v)`
///
/// The offset normal is the same as the basis normal. Domain equals the basis domain.
/// Analogous to OCCT `Geom_OffsetSurface`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OffsetSurface {
    pub basis: Box<Surface3>,
    /// Offset distance along the outward normal (positive = outward).
    pub offset_distance: f64,
}

/// A rectangular trimmed surface 鈥?a base surface restricted to the UV box
/// `[u1, u2] 脳 [v1, v2]`.
///
/// Evaluation delegates fully to the basis surface; only the reported domain
/// changes. Analogous to OCCT `Geom_RectangularTrimmedSurface`.
///
/// Appears in STEP as `RECTANGULAR_TRIMMED_SURFACE`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrimmedSurface {
    /// The underlying surface being trimmed.
    pub basis: Box<Surface3>,
    /// Trim bounds `[u1, u2, v1, v2]`.
    pub trim: [f64; 4],
}

impl TrimmedSurface {
    pub fn new(basis: Surface3, u1: f64, u2: f64, v1: f64, v2: f64) -> Self {
        Self {
            basis: Box::new(basis),
            trim: [u1, u2, v1, v2],
        }
    }
}

/// Surface formed by translating a 3D profile curve along a direction.
/// S(u,v) = profile.point_at(u) + v * direction
/// Analogous to OCCT Geom_SurfaceOfLinearExtrusion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinearExtrusionSurface {
    pub profile: Box<Curve3>,
    /// Normalized extrusion direction.
    pub direction: Vec3,
}

/// Surface formed by rotating a 3D profile curve around an axis.
/// S(u,v) = rotate(profile.point_at(v), axis_origin, axis_dir, angle=u)
/// Analogous to OCCT Geom_SurfaceOfRevolution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RevolutionSurface {
    pub profile: Box<Curve3>,
    pub axis_origin: Point3,
    /// Normalized rotation axis direction.
    pub axis_dir: Vec3,
}

/// Surface linearly interpolating between two 3D curves with a shared parameter domain.
/// S(u,v) = lerp(start.point_at(u), end.point_at(u), v)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RuledSurface {
    pub start: Box<Curve3>,
    pub end: Box<Curve3>,
}

/// A Coons patch blending four boundary curves over `[0,1] x [0,1]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoonsSurface {
    /// Boundary curve at `v = 0`, parameterized by `u`.
    pub south: Box<Curve3>,
    /// Boundary curve at `v = 1`, parameterized by `u`.
    pub north: Box<Curve3>,
    /// Boundary curve at `u = 0`, parameterized by `v`.
    pub west: Box<Curve3>,
    /// Boundary curve at `u = 1`, parameterized by `v`.
    pub east: Box<Curve3>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Surface3 {
    Plane(Plane),
    Cylinder(CylindricalSurface),
    Sphere(SphericalSurface),
    Cone(ConicalSurface),
    Torus(ToroidalSurface),
    Ellipsoid(EllipsoidalSurface),
    Helicoid(HelicoidSurface),
    Pipe(PipeSurface),
    BSpline(BSplineSurface),
    LinearExtrusion(LinearExtrusionSurface), // Phase K
    Revolution(RevolutionSurface),           // Phase K
    Ruled(RuledSurface),                     // Phase U
    Coons(CoonsSurface),                     // Phase V
    Bezier(BezierSurface),                   // Phase M
    TriBezier(TriBezierSurface),             // Phase T
    Offset(OffsetSurface),                   // Phase M
    Trimmed(TrimmedSurface),                 // Phase Q
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum PrimitiveSolid {
    Box {
        width: f64,
        height: f64,
        depth: f64,
    },
    Sphere {
        radius: f64,
    },
    Cylinder {
        radius: f64,
        height: f64,
    },
    Cone {
        base_radius: f64,
        height: f64,
    },
    Torus {
        major_radius: f64,
        minor_radius: f64,
    },
}

// 鈹€鈹€ 2D Geometry (parameter-space / PCurve types) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// A line in 2D parameter space: point + direction.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Line2d {
    pub origin: Point2,
    pub direction: Vec2,
}

/// A circle in 2D parameter space.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Circle2d {
    pub center: Point2,
    pub radius: f64,
}

/// An ellipse in 2D parameter space.
///
/// Analogous to OCCT `Geom2d_Ellipse`. Used as a PCurve when an edge traces
/// an elliptical path on the parameter domain of an adjacent surface.
///
/// Parametric form: `center + major_dir * a*cos(t) + minor_dir * b*sin(t)`
/// where `minor_dir = rotate_ccw_90(major_dir)`.  Default domain: `[0, 2蟺]`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Ellipse2d {
    pub center: Point2,
    /// Normalized major-axis direction in (u, v) space.
    pub major_dir: Vec2,
    pub major_radius: f64,
    pub minor_radius: f64,
}

/// A 2D involute of a base circle in parameter space.
///
/// Parametric form around the local x-axis:
/// `x(t) = r * (cos t + t sin t)`
/// `y(t) = r * (sin t - t cos t)`
///
/// The local frame is then rotated by `start_angle` and translated by `center`.
/// This curve is commonly used for gear-tooth flank profiles.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CircleInvolute2d {
    pub center: Point2,
    pub base_radius: f64,
    /// Rotation of the local involute frame in radians.
    pub start_angle: f64,
}

/// A 2D Archimedean spiral in parameter space.
///
/// `r(t) = a + b*t`, `theta(t) = start_angle + t`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ArchimedeanSpiral2d {
    pub center: Point2,
    pub a: f64,
    pub b: f64,
    pub start_angle: f64,
}

/// A 2D logarithmic spiral in parameter space.
///
/// `r(t) = a * exp(b*t)`, `theta(t) = start_angle + t`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct LogarithmicSpiral2d {
    pub center: Point2,
    pub a: f64,
    pub b: f64,
    pub start_angle: f64,
}

/// A 2D sine-wave curve in parameter space.
///
/// Parametric form:
/// `x(t) = t`
/// `y(t) = amplitude * sin(frequency * t + phase)`
///
/// Useful for procedural sketching and for matching OCCT's sine-wave evaluator
/// family in a lightweight form.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SineWave2d {
    pub amplitude: f64,
    pub frequency: f64,
    pub phase: f64,
}

/// A non-uniform rational B-spline curve in 2D parameter space.
///
/// Analogous to OCCT `Geom2d_BSplineCurve`. Used for PCurves: the image of
/// a 3D edge in the (u, v) domain of an adjacent surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BSplineCurve2 {
    pub degree: usize,
    /// Full knot vector (with multiplicities expanded).
    pub knots: Vec<f64>,
    pub control_points: Vec<DVec2>,
    /// Homogeneous weights; 1.0 for non-rational.
    pub weights: Vec<f64>,
}

/// A rational or non-rational Bezier curve in 2D parameter space.
///
/// Analogous to OCCT `Geom2d_BezierCurve`. Domain is `[0, 1]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BezierCurve2 {
    pub control_points: Vec<DVec2>,
    /// Homogeneous weights; 1.0 for non-rational.
    pub weights: Vec<f64>,
}

/// A curve defined in the 2D parameter space (u, v) of a surface.
///
/// Used for PCurves: the image of a 3D edge on the parameter domain of an
/// adjacent face surface. Analogous to OCCT `Geom2d_Curve`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Curve2d {
    Line(Line2d),
    Circle(Circle2d),
    Ellipse(Ellipse2d), // Phase J
    CircleInvolute(CircleInvolute2d), // Phase P7
    ArchimedeanSpiral(ArchimedeanSpiral2d), // Phase P7
    LogarithmicSpiral(LogarithmicSpiral2d), // Phase P7
    SineWave(SineWave2d), // Phase P7
    BSpline(BSplineCurve2),
    Bezier(BezierCurve2), // Phase M
}

// 鈹€鈹€ Geometric evaluation traits 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Returns a vector perpendicular to `v`. Stable for any non-zero input.
pub fn any_perpendicular(v: DVec3) -> DVec3 {
    // Pick the axis least aligned with v, then cross.
    let abs = v.abs();
    let candidate = if abs.x <= abs.y && abs.x <= abs.z {
        DVec3::X
    } else if abs.y <= abs.z {
        DVec3::Y
    } else {
        DVec3::Z
    };
    v.cross(candidate).normalize()
}

fn orthonormal_frame(axis: DVec3, ref_dir: DVec3) -> (DVec3, DVec3, DVec3) {
    let axis = axis.normalize_or_zero();
    let mut x_axis = ref_dir - axis * ref_dir.dot(axis);
    if x_axis.length_squared() <= 1e-24 {
        x_axis = any_perpendicular(axis);
    } else {
        x_axis = x_axis.normalize();
    }
    let y_axis = axis.cross(x_axis).normalize_or_zero();
    (axis, x_axis, y_axis)
}

/// Parametric evaluation of a 3D curve: `t 鈫?Point3`.
///
/// Mirrors OCCT `Geom_Curve::Value(t)` / `D1(t)`.
pub trait CurveEval {
    /// Point on the curve at parameter `t`.
    fn point_at(&self, t: f64) -> DVec3;
    /// Unit tangent vector at parameter `t`.
    fn tangent_at(&self, t: f64) -> DVec3;
    /// Natural parameter domain `[t_min, t_max]`.
    /// Lines use `[NEG_INFINITY, INFINITY]`; circles/ellipses use `[0, 2蟺]`.
    fn default_domain(&self) -> [f64; 2];
}

/// Parametric evaluation of a 3D surface: `(u, v) 鈫?Point3`.
///
/// Mirrors OCCT `Geom_Surface::Value(u, v)`.
pub trait SurfaceEval {
    /// Point on the surface at parameter `(u, v)`.
    fn point_at(&self, u: f64, v: f64) -> DVec3;
    /// Outward unit normal at parameter `(u, v)`.
    fn normal_at(&self, u: f64, v: f64) -> DVec3;
    /// Natural parameter domain `[u_min, u_max, v_min, v_max]`.
    fn default_domain(&self) -> [f64; 4];
}

/// Parametric evaluation of a 2D curve (PCurve): `t 鈫?Point2`.
pub trait Curve2dEval {
    fn point_at(&self, t: f64) -> DVec2;
}

// 鈹€鈹€ CurveEval implementations 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

impl CurveEval for Line3 {
    fn point_at(&self, t: f64) -> DVec3 {
        self.origin + t * self.direction
    }
    fn tangent_at(&self, _t: f64) -> DVec3 {
        self.direction
    }
    fn default_domain(&self) -> [f64; 2] {
        [f64::NEG_INFINITY, f64::INFINITY]
    }
}

impl CurveEval for Circle3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.normal);
        let y_ax = self.normal.cross(x_ax);
        self.center + self.radius * (t.cos() * x_ax + t.sin() * y_ax)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.normal);
        let y_ax = self.normal.cross(x_ax);
        (-t.sin() * x_ax + t.cos() * y_ax).normalize()
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 2.0 * PI]
    }
}

impl CurveEval for Ellipse3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let x_ax = self.major_dir;
        let y_ax = self.normal.cross(x_ax).normalize();
        self.center + self.major_radius * t.cos() * x_ax + self.minor_radius * t.sin() * y_ax
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let x_ax = self.major_dir;
        let y_ax = self.normal.cross(x_ax).normalize();
        (-self.major_radius * t.sin() * x_ax + self.minor_radius * t.cos() * y_ax).normalize()
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 2.0 * PI]
    }
}

impl CurveEval for Hyperbola3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let minor_dir = self.normal.cross(self.major_dir).normalize();
        self.center
            + self.semi_major * t.cosh() * self.major_dir
            + self.semi_minor * t.sinh() * minor_dir
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let minor_dir = self.normal.cross(self.major_dir).normalize();
        let v =
            self.semi_major * t.sinh() * self.major_dir + self.semi_minor * t.cosh() * minor_dir;
        v.normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        [-1e4, 1e4] // unbounded; caller trims as needed
    }
}

impl CurveEval for Parabola3 {
    fn point_at(&self, t: f64) -> DVec3 {
        // dir_perp forms a right-handed system: axis_dir × normal gives perpendicular direction
        let dir_perp = self.axis_dir.cross(self.normal).normalize();
        self.vertex + (t * t / (2.0 * self.focal_param)) * self.axis_dir + t * dir_perp
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let dir_perp = self.axis_dir.cross(self.normal).normalize();
        let v = (t / self.focal_param) * self.axis_dir + dir_perp;
        v.normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        [-1e4, 1e4] // unbounded
    }
}

impl CurveEval for CircularHelix3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let axis = self.axis.normalize_or_zero();
        let mut x_axis = self.ref_dir - axis * self.ref_dir.dot(axis);
        if x_axis.length_squared() <= 1e-24 {
            x_axis = any_perpendicular(axis);
        } else {
            x_axis = x_axis.normalize();
        }
        let y_axis = axis.cross(x_axis).normalize_or_zero();
        let lead = self.pitch / (2.0 * PI);
        self.origin + self.radius * (t.cos() * x_axis + t.sin() * y_axis) + (lead * t) * axis
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let axis = self.axis.normalize_or_zero();
        let mut x_axis = self.ref_dir - axis * self.ref_dir.dot(axis);
        if x_axis.length_squared() <= 1e-24 {
            x_axis = any_perpendicular(axis);
        } else {
            x_axis = x_axis.normalize();
        }
        let y_axis = axis.cross(x_axis).normalize_or_zero();
        let lead = self.pitch / (2.0 * PI);
        (-self.radius * t.sin() * x_axis + self.radius * t.cos() * y_axis + lead * axis)
            .normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        [-1e4, 1e4]
    }
}

impl CurveEval for SineWave3 {
    fn point_at(&self, t: f64) -> DVec3 {
        self.origin
            + t * self.baseline_dir
            + self.amplitude * (self.frequency * t + self.phase).sin() * self.amplitude_dir
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let v = self.baseline_dir
            + self.amplitude * self.frequency * (self.frequency * t + self.phase).cos()
                * self.amplitude_dir;
        v.normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        [-1e4, 1e4]
    }
}

impl CurveEval for Curve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line(c) => c.point_at(t),
            Curve3::Circle(c) => c.point_at(t),
            Curve3::Ellipse(c) => c.point_at(t),
            Curve3::BSpline(c) => c.point_at(t),
            Curve3::Bezier(c) => c.point_at(t),
            Curve3::Offset(c) => c.point_at(t),
            Curve3::Hyperbola(c) => c.point_at(t),
            Curve3::Parabola(c) => c.point_at(t),
            Curve3::CircularHelix(c) => c.point_at(t),
            Curve3::SineWave(c) => c.point_at(t),
        }
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line(c) => c.tangent_at(t),
            Curve3::Circle(c) => c.tangent_at(t),
            Curve3::Ellipse(c) => c.tangent_at(t),
            Curve3::BSpline(c) => c.tangent_at(t),
            Curve3::Bezier(c) => c.tangent_at(t),
            Curve3::Offset(c) => c.tangent_at(t),
            Curve3::Hyperbola(c) => c.tangent_at(t),
            Curve3::Parabola(c) => c.tangent_at(t),
            Curve3::CircularHelix(c) => c.tangent_at(t),
            Curve3::SineWave(c) => c.tangent_at(t),
        }
    }
    fn default_domain(&self) -> [f64; 2] {
        match self {
            Curve3::Line(c) => c.default_domain(),
            Curve3::Circle(c) => c.default_domain(),
            Curve3::Ellipse(c) => c.default_domain(),
            Curve3::BSpline(c) => c.default_domain(),
            Curve3::Bezier(c) => c.default_domain(),
            Curve3::Offset(c) => c.default_domain(),
            Curve3::Hyperbola(c) => c.default_domain(),
            Curve3::Parabola(c) => c.default_domain(),
            Curve3::CircularHelix(c) => c.default_domain(),
            Curve3::SineWave(c) => c.default_domain(),
        }
    }
}

// 鈹€鈹€ SurfaceEval implementations 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

impl SurfaceEval for Plane {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.normal);
        let y_ax = self.normal.cross(x_ax);
        self.origin + u * x_ax + v * y_ax
    }
    fn normal_at(&self, _u: f64, _v: f64) -> DVec3 {
        self.normal
    }
    fn default_domain(&self) -> [f64; 4] {
        [
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
        ]
    }
}

impl SurfaceEval for CylindricalSurface {
    /// u = azimuth angle [0, 2蟺], v = height along axis.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        self.origin + self.radius * (u.cos() * x_ax + u.sin() * y_ax) + v * self.axis
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        (u.cos() * x_ax + u.sin() * y_ax).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, f64::NEG_INFINITY, f64::INFINITY]
    }
}

impl SurfaceEval for SphericalSurface {
    /// u = longitude [0, 2蟺], v = colatitude [0, 蟺] (0 = north pole).
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        self.center
            + self.radius * (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * self.axis)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * self.axis).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, PI]
    }
}

impl SurfaceEval for ConicalSurface {
    /// u = azimuth [0, 2π], v = distance along the cone generatrix from the
    /// reference circle at `self.apex`.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let axis = self.axis_dir();
        let x_ax = any_perpendicular(axis);
        let y_ax = axis.cross(x_ax).normalize();
        let radial = self.radius_at_slant(v);
        let axial = self.axial_from_slant(v);
        self.apex + axial * axis + radial * (u.cos() * x_ax + u.sin() * y_ax)
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let axis = self.axis_dir();
        let x_ax = any_perpendicular(axis);
        let y_ax = axis.cross(x_ax).normalize();
        let radial = u.cos() * x_ax + u.sin() * y_ax;
        let half = self.half_angle_rad;
        (radial * half.cos() - axis * half.sin()).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, f64::INFINITY]
    }
}

impl SurfaceEval for ToroidalSurface {
    /// u = major angle [0, 2蟺], v = minor angle [0, 2蟺].
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        let tube_center = self.center + self.major_radius * (u.cos() * x_ax + u.sin() * y_ax);
        let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();
        tube_center + self.minor_radius * (v.cos() * radial + v.sin() * self.axis)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = any_perpendicular(self.axis);
        let y_ax = self.axis.cross(x_ax).normalize();
        let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();
        (v.cos() * radial + v.sin() * self.axis).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, 2.0 * PI]
    }
}

impl SurfaceEval for EllipsoidalSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        self.center
            + self.radius_x * v.sin() * u.cos() * x_axis
            + self.radius_y * v.sin() * u.sin() * y_axis
            + self.radius_z * v.cos() * axis
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        let p = self.point_at(u, v) - self.center;
        let x = p.dot(x_axis);
        let y = p.dot(y_axis);
        let z = p.dot(axis);
        let grad = (x / (self.radius_x * self.radius_x)) * x_axis
            + (y / (self.radius_y * self.radius_y)) * y_axis
            + (z / (self.radius_z * self.radius_z)) * axis;
        grad.normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, PI]
    }
}

impl SurfaceEval for HelicoidSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        let lead = self.pitch / (2.0 * PI);
        self.origin + v * (u.cos() * x_axis + u.sin() * y_axis) + (lead * u) * axis
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        let lead = self.pitch / (2.0 * PI);
        let du = v * (-u.sin() * x_axis + u.cos() * y_axis) + lead * axis;
        let dv = u.cos() * x_axis + u.sin() * y_axis;
        du.cross(dv).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 4] {
        [-2.0 * PI, 2.0 * PI, -10.0, 10.0]
    }
}

impl SurfaceEval for LinearExtrusionSurface {
    /// u = profile parameter, v = extrusion distance along direction.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        self.profile.point_at(u) + v * self.direction
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let tangent = self.profile.tangent_at(u);
        let n = tangent.cross(self.direction);
        if n.length_squared() < 1e-20 {
            return DVec3::Z;
        }
        n.normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        let [t1, t2] = self.profile.default_domain();
        [t1, t2, f64::NEG_INFINITY, f64::INFINITY]
    }
}

impl SurfaceEval for RevolutionSurface {
    /// u = azimuth angle [0, 2蟺], v = profile parameter.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let p = self.profile.point_at(v);
        let d = p - self.axis_origin;
        let d_par = self.axis_dir * d.dot(self.axis_dir);
        let d_perp = d - d_par;
        self.axis_origin + d_par + d_perp * u.cos() + self.axis_dir.cross(d_perp) * u.sin()
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-6;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (self.point_at(u, v + eps) - self.point_at(u, v - eps)) / (2.0 * eps);
        let n = du.cross(dv);
        if n.length_squared() < 1e-20 {
            return DVec3::Z;
        }
        n.normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        let [t1, t2] = self.profile.default_domain();
        [0.0, 2.0 * PI, t1, t2]
    }
}

impl SurfaceEval for TrimmedSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        self.basis.point_at(u, v)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        self.basis.normal_at(u, v)
    }
    fn default_domain(&self) -> [f64; 4] {
        self.trim
    }
}

impl SurfaceEval for Surface3 {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        match self {
            Surface3::Plane(s) => s.point_at(u, v),
            Surface3::Cylinder(s) => s.point_at(u, v),
            Surface3::Sphere(s) => s.point_at(u, v),
            Surface3::Cone(s) => s.point_at(u, v),
            Surface3::Torus(s) => s.point_at(u, v),
            Surface3::Ellipsoid(s) => s.point_at(u, v),
            Surface3::Helicoid(s) => s.point_at(u, v),
            Surface3::Pipe(s) => s.point_at(u, v),
            Surface3::BSpline(s) => s.point_at(u, v),
            Surface3::LinearExtrusion(s) => s.point_at(u, v),
            Surface3::Revolution(s) => s.point_at(u, v),
            Surface3::Ruled(s) => s.point_at(u, v),
            Surface3::Coons(s) => s.point_at(u, v),
            Surface3::Bezier(s) => s.point_at(u, v),
            Surface3::TriBezier(s) => s.point_at(u, v),
            Surface3::Offset(s) => s.point_at(u, v),
            Surface3::Trimmed(s) => s.point_at(u, v),
        }
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        match self {
            Surface3::Plane(s) => s.normal_at(u, v),
            Surface3::Cylinder(s) => s.normal_at(u, v),
            Surface3::Sphere(s) => s.normal_at(u, v),
            Surface3::Cone(s) => s.normal_at(u, v),
            Surface3::Torus(s) => s.normal_at(u, v),
            Surface3::Ellipsoid(s) => s.normal_at(u, v),
            Surface3::Helicoid(s) => s.normal_at(u, v),
            Surface3::Pipe(s) => s.normal_at(u, v),
            Surface3::BSpline(s) => s.normal_at(u, v),
            Surface3::LinearExtrusion(s) => s.normal_at(u, v),
            Surface3::Revolution(s) => s.normal_at(u, v),
            Surface3::Ruled(s) => s.normal_at(u, v),
            Surface3::Coons(s) => s.normal_at(u, v),
            Surface3::Bezier(s) => s.normal_at(u, v),
            Surface3::TriBezier(s) => s.normal_at(u, v),
            Surface3::Offset(s) => s.normal_at(u, v),
            Surface3::Trimmed(s) => s.normal_at(u, v),
        }
    }
    fn default_domain(&self) -> [f64; 4] {
        match self {
            Surface3::Plane(s) => s.default_domain(),
            Surface3::Cylinder(s) => s.default_domain(),
            Surface3::Sphere(s) => s.default_domain(),
            Surface3::Cone(s) => s.default_domain(),
            Surface3::Torus(s) => s.default_domain(),
            Surface3::Ellipsoid(s) => s.default_domain(),
            Surface3::Helicoid(s) => s.default_domain(),
            Surface3::Pipe(s) => s.default_domain(),
            Surface3::BSpline(s) => s.default_domain(),
            Surface3::LinearExtrusion(s) => s.default_domain(),
            Surface3::Revolution(s) => s.default_domain(),
            Surface3::Ruled(s) => s.default_domain(),
            Surface3::Coons(s) => s.default_domain(),
            Surface3::Bezier(s) => s.default_domain(),
            Surface3::TriBezier(s) => s.default_domain(),
            Surface3::Offset(s) => s.default_domain(),
            Surface3::Trimmed(s) => s.default_domain(),
        }
    }
}

/// Evaluate all Lagrange basis functions for the given nodes at t.
///
/// Uses safer numerical handling with explicit tolerance for near-singular cases.
/// Returns basis functions that satisfy partition of unity (sum = 1).
fn lagrange_basis(nodes: &[f64], t: f64) -> Vec<f64> {
    let n = nodes.len();
    if n == 0 {
        return vec![];
    }

    let mut basis = vec![1.0; n];
    let tol = 1e-14;

    for i in 0..n {
        for j in 0..n {
            if i != j {
                let denom = nodes[i] - nodes[j];
                if denom.abs() > tol {
                    basis[i] *= (t - nodes[j]) / denom;
                } else {
                    // Nodes too close - this indicates invalid input,
                    // but we handle gracefully by setting to 0
                    basis[i] = 0.0;
                }
            }
        }

        // Guard against NaN/Inf
        if !basis[i].is_finite() {
            basis[i] = 0.0;
        }
    }

    // Ensure partition of unity for stability
    let sum: f64 = basis.iter().sum();
    if sum.abs() > tol {
        for b in &mut basis {
            *b /= sum;
        }
    }

    basis
}

fn remap_unit_to_curve_domain(curve: &Curve3, t: f64) -> f64 {
    let [t0, t1] = curve.default_domain();
    if !t0.is_finite() || !t1.is_finite() {
        return t;
    }
    t0 + (t1 - t0) * t
}

fn projected_frame_from_tangent(tangent: DVec3, ref_dir: DVec3) -> (DVec3, DVec3) {
    let tangent = tangent.normalize_or_zero();
    let mut x_axis = ref_dir - tangent * ref_dir.dot(tangent);
    if x_axis.length_squared() <= 1e-24 {
        x_axis = any_perpendicular(tangent);
    } else {
        x_axis = x_axis.normalize();
    }
    let y_axis = tangent.cross(x_axis).normalize_or_zero();
    (x_axis, y_axis)
}

impl SurfaceEval for PipeSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let center = self.spine.point_at(v);
        let tangent = self.spine.tangent_at(v);
        let (x_axis, y_axis) = projected_frame_from_tangent(tangent, self.ref_dir);
        center + self.radius * (u.cos() * x_axis + u.sin() * y_axis)
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = self.point_at(u + eps, v) - self.point_at(u - eps, v);
        let dv = self.point_at(u, v + eps) - self.point_at(u, v - eps);
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        let [v0, v1] = self.spine.default_domain();
        [0.0, 2.0 * PI, v0, v1]
    }
}

impl SurfaceEval for CoonsSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let south = self.south.point_at(remap_unit_to_curve_domain(&self.south, u));
        let north = self.north.point_at(remap_unit_to_curve_domain(&self.north, u));
        let west = self.west.point_at(remap_unit_to_curve_domain(&self.west, v));
        let east = self.east.point_at(remap_unit_to_curve_domain(&self.east, v));

        let p00 = self.south.point_at(remap_unit_to_curve_domain(&self.south, 0.0));
        let p10 = self.south.point_at(remap_unit_to_curve_domain(&self.south, 1.0));
        let p01 = self.north.point_at(remap_unit_to_curve_domain(&self.north, 0.0));
        let p11 = self.north.point_at(remap_unit_to_curve_domain(&self.north, 1.0));

        let linear_u = south * (1.0 - v) + north * v;
        let linear_v = west * (1.0 - u) + east * u;
        let bilinear = p00 * ((1.0 - u) * (1.0 - v))
            + p10 * (u * (1.0 - v))
            + p01 * ((1.0 - u) * v)
            + p11 * (u * v);
        linear_u + linear_v - bilinear
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = self.point_at((u + eps).clamp(0.0, 1.0), v)
            - self.point_at((u - eps).clamp(0.0, 1.0), v);
        let dv = self.point_at(u, (v + eps).clamp(0.0, 1.0))
            - self.point_at(u, (v - eps).clamp(0.0, 1.0));
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        [0.0, 1.0, 0.0, 1.0]
    }
}

// 鈹€鈹€ BSpline evaluation 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// De Boor's algorithm in homogeneous 4D space.
/// Returns `[wx, wy, wz, w]` (not divided by w yet).
fn de_boor_homo(
    degree: usize,
    knots: &[f64],
    points: &[DVec3],
    weights: &[f64],
    t: f64,
) -> [f64; 4] {
    let n = points.len();
    if n == 0 {
        return [0.0; 4];
    }
    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };
    let mut d: Vec<[f64; 4]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
            let w = weights[idx];
            [points[idx].x * w, points[idx].y * w, points[idx].z * w, w]
        })
        .collect();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }
    d[degree]
}

/// De Boor's algorithm for rational B-spline evaluation.
/// Returns the 3D point at parameter `t`.
fn de_boor(degree: usize, knots: &[f64], points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n == 0 {
        return DVec3::ZERO;
    }

    // Find knot span index k such that knots[k] <= t < knots[k+1]
    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };

    // Initialize homogeneous control points for the span
    let mut d: Vec<[f64; 4]> = (0..=degree)
        .map(|j| {
            let idx = k - degree + j;
            let idx = idx.min(n - 1);
            let w = weights[idx];
            [points[idx].x * w, points[idx].y * w, points[idx].z * w, w]
        })
        .collect();

    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }

    let w = d[degree][3];
    if w.abs() < 1e-15 {
        DVec3::ZERO
    } else {
        DVec3::new(d[degree][0] / w, d[degree][1] / w, d[degree][2] / w)
    }
}

/// De Boor's algorithm for rational B-spline evaluation in 2D parameter space.
/// Returns the 2D point at parameter `t`. Identical logic to `de_boor` with DVec2.
fn de_boor_2d(degree: usize, knots: &[f64], points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n == 0 {
        return DVec2::ZERO;
    }

    let k = {
        let t_min = knots[degree];
        let t_max = knots[knots.len() - degree - 1];
        let t_clamped = t.clamp(t_min, t_max);
        let mut span = degree;
        for (i, &knot) in knots
            .iter()
            .enumerate()
            .take(knots.len() - degree - 1)
            .skip(degree)
        {
            if knot <= t_clamped {
                span = i;
            } else {
                break;
            }
        }
        span
    };

    // Homogeneous control points [x*w, y*w, w]
    let mut d: Vec<[f64; 3]> = (0..=degree)
        .map(|j| {
            let idx = (k - degree + j).min(n - 1);
            let w = weights[idx];
            [points[idx].x * w, points[idx].y * w, w]
        })
        .collect();

    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = k - degree + j;
            let denom = knots[i + degree - r + 1] - knots[i];
            let alpha = if denom.abs() < 1e-15 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            let prev = d[j - 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(prev.iter()) {
                *elem = (1.0 - alpha) * p + alpha * *elem;
            }
        }
    }

    let w = d[degree][2];
    if w.abs() < 1e-15 {
        DVec2::ZERO
    } else {
        DVec2::new(d[degree][0] / w, d[degree][1] / w)
    }
}

// 鈹€鈹€ Analytic curve derivative helpers 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Analytic tangent for a rational B-Spline curve (NURBS) using the quotient rule.
///
/// The derivative of C(t) = A(t)/W(t) is:
///   C'(t) = (A'(t) 鈭?W'(t)路C(t)) / W(t)
///
/// A'(t) and W'(t) are degree-(p鈭?) B-Splines with control points:
///   A'_i = p 路 (w_{i+1}路P_{i+1} 鈭?w_i路P_i) / (t_{i+p+1} 鈭?t_{i+1})
///   W'_i = p 路 (w_{i+1} 鈭?w_i)              / (t_{i+p+1} 鈭?t_{i+1})
///
/// Returns the unnormalised derivative vector (caller normalises if needed).
fn bspline_tangent_analytic(
    degree: usize,
    knots: &[f64],
    points: &[DVec3],
    weights: &[f64],
    t: f64,
) -> DVec3 {
    let n = points.len();
    if n < 2 || degree == 0 {
        return DVec3::ZERO;
    }

    let p = degree as f64;
    let m = n - 1; // number of derivative control points

    let mut a_prime: Vec<DVec3> = Vec::with_capacity(m);
    let mut w_prime: Vec<DVec3> = Vec::with_capacity(m); // scalar stored in .x
    for i in 0..m {
        let denom = knots[i + degree + 1] - knots[i + 1];
        if denom.abs() < 1e-15 {
            a_prime.push(DVec3::ZERO);
            w_prime.push(DVec3::ZERO);
        } else {
            let s = p / denom;
            a_prime.push(s * (weights[i + 1] * points[i + 1] - weights[i] * points[i]));
            w_prime.push(DVec3::new(s * (weights[i + 1] - weights[i]), 0.0, 0.0));
        }
    }

    let deriv_knots = &knots[1..knots.len() - 1];
    let unit = vec![1.0f64; m];

    // A'(t): non-rational B-Spline of degree p-1
    let a_prime_t = de_boor(degree - 1, deriv_knots, &a_prime, &unit, t);
    // W'(t): scalar B-Spline of degree p-1 (embedded in .x)
    let w_prime_t = de_boor(degree - 1, deriv_knots, &w_prime, &unit, t).x;

    // W(t) and C(t) from the homogeneous evaluation
    let h = de_boor_homo(degree, knots, points, weights, t);
    let w_t = h[3];
    if w_t.abs() < 1e-15 {
        return DVec3::ZERO;
    }
    let c_t = DVec3::new(h[0] / w_t, h[1] / w_t, h[2] / w_t);

    (a_prime_t - w_prime_t * c_t) / w_t
}

/// Analytic tangent for a rational Bezier curve using the quotient rule.
///
/// The derivative of a degree-n Bezier is a degree-(n鈭?) Bezier with:
///   A'_i = n路(w_{i+1}路P_{i+1} 鈭?w_i路P_i)
///   W'_i = n路(w_{i+1} 鈭?w_i)
fn bezier_tangent_analytic(points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n < 2 {
        return DVec3::ZERO;
    }
    let deg = (n - 1) as f64;

    let mut a_prime: Vec<DVec3> = Vec::with_capacity(n - 1);
    let mut w_prime: Vec<DVec3> = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        a_prime.push(deg * (weights[i + 1] * points[i + 1] - weights[i] * points[i]));
        w_prime.push(DVec3::new(deg * (weights[i + 1] - weights[i]), 0.0, 0.0));
    }

    let unit = vec![1.0f64; n - 1];
    let a_prime_t = de_casteljau_3d(&a_prime, &unit, t);
    let w_prime_t = de_casteljau_3d(&w_prime, &unit, t).x;

    // W(t): evaluate weights as scalar Bezier (embed in .x with unit weights)
    let w_pts: Vec<DVec3> = weights.iter().map(|&w| DVec3::new(w, 0.0, 0.0)).collect();
    let w_unit = vec![1.0f64; n]; // n elements to match w_pts
    let w_t = de_casteljau_3d(&w_pts, &w_unit, t).x;
    if w_t.abs() < 1e-15 {
        return DVec3::ZERO;
    }

    // C(t) from the standard rational evaluation
    let c_t = de_casteljau_3d(points, weights, t);

    (a_prime_t - w_prime_t * c_t) / w_t
}

impl CurveEval for BSplineCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        de_boor(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        bspline_tangent_analytic(self.degree, &self.knots, &self.control_points, &self.weights, t)
            .normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        let d = self.degree;
        let n = self.knots.len();
        if n < 2 * d + 2 {
            return [0.0, 1.0];
        }
        [self.knots[d], self.knots[n - d - 1]]
    }
}

impl SurfaceEval for BSplineSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        // Tensor product rational evaluation (NURBS):
        // 1. For each v-column, evaluate the u-direction NURBS in homogeneous coords
        //    鈫?get (wx, wy, wz, w) for each column index.
        // 2. Collect column weights and weighted positions.
        // 3. Run de Boor in v on the homogeneous results, then divide by weight.
        let n_u = self.control_points.len();
        if n_u == 0 {
            return DVec3::ZERO;
        }
        let n_v = self.control_points[0].len();
        if n_v == 0 {
            return DVec3::ZERO;
        }
        // Step 1: evaluate each v-column in the u direction 鈫?homogeneous 4-vector
        let col_homo: Vec<[f64; 4]> = (0..n_v)
            .map(|j| {
                let pts: Vec<DVec3> = (0..n_u).map(|i| self.control_points[i][j]).collect();
                let wts: Vec<f64> = (0..n_u).map(|i| self.weights[i][j]).collect();
                de_boor_homo(self.degree_u, &self.knots_u, &pts, &wts, u)
            })
            .collect();
        // Step 2: build the v-direction "control points" and "weights" from col_homo
        let v_pts: Vec<DVec3> = col_homo
            .iter()
            .map(|h| {
                let w = h[3];
                if w.abs() < 1e-15 {
                    DVec3::ZERO
                } else {
                    DVec3::new(h[0] / w, h[1] / w, h[2] / w)
                }
            })
            .collect();
        let v_wts: Vec<f64> = col_homo.iter().map(|h| h[3]).collect();
        // Step 3: rational de Boor in v
        de_boor(self.degree_v, &self.knots_v, &v_pts, &v_wts, v)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let [_u0, u1, _v0, v1] = self.default_domain();
        let du = if u + eps <= u1 {
            self.point_at(u + eps, v) - self.point_at(u, v)
        } else {
            self.point_at(u, v) - self.point_at(u - eps, v)
        };
        let dv = if v + eps <= v1 {
            self.point_at(u, v + eps) - self.point_at(u, v)
        } else {
            self.point_at(u, v) - self.point_at(u, v - eps)
        };
        let n = du.cross(dv);
        let len = n.length();
        if len < 1e-15 { DVec3::Z } else { n / len }
    }
    fn default_domain(&self) -> [f64; 4] {
        let du = self.degree_u;
        let dv = self.degree_v;
        let nu = self.knots_u.len();
        let nv = self.knots_v.len();
        let u0 = if nu > du { self.knots_u[du] } else { 0.0 };
        let u1 = if nu > du + 1 {
            self.knots_u[nu - du - 1]
        } else {
            1.0
        };
        let v0 = if nv > dv { self.knots_v[dv] } else { 0.0 };
        let v1 = if nv > dv + 1 {
            self.knots_v[nv - dv - 1]
        } else {
            1.0
        };
        [u0, u1, v0, v1]
    }
}

// 鈹€鈹€ Curve2dEval implementations 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

impl Curve2dEval for Line2d {
    fn point_at(&self, t: f64) -> DVec2 {
        self.origin + t * self.direction
    }
}

impl Curve2dEval for Circle2d {
    fn point_at(&self, t: f64) -> DVec2 {
        self.center + self.radius * DVec2::new(t.cos(), t.sin())
    }
}

impl Curve2dEval for Ellipse2d {
    fn point_at(&self, t: f64) -> DVec2 {
        // minor_dir = rotate major_dir by 90掳 counter-clockwise
        let minor_dir = DVec2::new(-self.major_dir.y, self.major_dir.x);
        self.center
            + self.major_dir * (self.major_radius * t.cos())
            + minor_dir * (self.minor_radius * t.sin())
    }
}

impl Curve2dEval for CircleInvolute2d {
    fn point_at(&self, t: f64) -> DVec2 {
        let r = self.base_radius.max(0.0);
        let x = r * (t.cos() + t * t.sin());
        let y = r * (t.sin() - t * t.cos());

        let ca = self.start_angle.cos();
        let sa = self.start_angle.sin();
        let xr = x * ca - y * sa;
        let yr = x * sa + y * ca;
        self.center + DVec2::new(xr, yr)
    }
}

impl Curve2dEval for ArchimedeanSpiral2d {
    fn point_at(&self, t: f64) -> DVec2 {
        let r = self.a + self.b * t;
        let th = self.start_angle + t;
        self.center + DVec2::new(r * th.cos(), r * th.sin())
    }
}

impl Curve2dEval for LogarithmicSpiral2d {
    fn point_at(&self, t: f64) -> DVec2 {
        let r = self.a * (self.b * t).exp();
        let th = self.start_angle + t;
        self.center + DVec2::new(r * th.cos(), r * th.sin())
    }
}

impl Curve2dEval for SineWave2d {
    fn point_at(&self, t: f64) -> DVec2 {
        DVec2::new(t, self.amplitude * (self.frequency * t + self.phase).sin())
    }
}

impl Curve2dEval for BSplineCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        de_boor_2d(
            self.degree,
            &self.knots,
            &self.control_points,
            &self.weights,
            t,
        )
    }
}

impl Curve2dEval for Curve2d {
    fn point_at(&self, t: f64) -> DVec2 {
        match self {
            Curve2d::Line(c) => c.point_at(t),
            Curve2d::Circle(c) => c.point_at(t),
            Curve2d::Ellipse(c) => c.point_at(t),
            Curve2d::CircleInvolute(c) => c.point_at(t),
            Curve2d::ArchimedeanSpiral(c) => c.point_at(t),
            Curve2d::LogarithmicSpiral(c) => c.point_at(t),
            Curve2d::SineWave(c) => c.point_at(t),
            Curve2d::BSpline(c) => c.point_at(t),
            Curve2d::Bezier(c) => c.point_at(t),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circle_involute_starts_on_base_circle() {
        let inv = CircleInvolute2d {
            center: DVec2::new(2.0, -1.0),
            base_radius: 3.0,
            start_angle: 0.0,
        };
        let p0 = inv.point_at(0.0);
        assert!((p0.x - 5.0).abs() < 1e-12);
        assert!((p0.y + 1.0).abs() < 1e-12);
    }

    #[test]
    fn archimedean_spiral_point_progresses_radially() {
        let s = ArchimedeanSpiral2d {
            center: DVec2::ZERO,
            a: 1.0,
            b: 0.5,
            start_angle: 0.0,
        };
        let p0 = s.point_at(0.0);
        let p1 = s.point_at(2.0);
        assert!(p1.length() > p0.length(), "spiral radius should increase with t");
    }

    #[test]
    fn logarithmic_spiral_grows_exponentially() {
        let s = LogarithmicSpiral2d {
            center: DVec2::ZERO,
            a: 1.0,
            b: 0.4,
            start_angle: 0.0,
        };
        let p0 = s.point_at(0.0);
        let p1 = s.point_at(2.0);
        assert!(p1.length() > p0.length() * 1.5, "log spiral should grow faster than linear at this sample");
    }

    #[test]
    fn sine_wave_samples_match_expected_values() {
        let s = SineWave2d {
            amplitude: 2.0,
            frequency: 1.0,
            phase: 0.0,
        };
        let p0 = s.point_at(0.0);
        let p90 = s.point_at(std::f64::consts::FRAC_PI_2);
        assert!((p0.x - 0.0).abs() < 1e-12);
        assert!((p0.y - 0.0).abs() < 1e-12);
        assert!((p90.y - 2.0).abs() < 1e-12);
    }

    #[test]
    fn curve2d_sine_wave_variant_dispatches_evaluator() {
        let c = Curve2d::SineWave(SineWave2d {
            amplitude: 1.5,
            frequency: 2.0,
            phase: 0.25,
        });
        let t = 0.3;
        let p = c.point_at(t);
        let expected_y = 1.5 * (2.0 * t + 0.25).sin();
        assert!((p.x - t).abs() < 1e-12);
        assert!((p.y - expected_y).abs() < 1e-12);
    }

    #[test]
    fn sine_wave3_origin_phase_zero_evaluates_at_zero_offset() {
        let c = SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 3.0,
            frequency: 1.0,
            phase: 0.0,
        };
        // At t=0, sin(0)=0 → point should be at origin.
        let p = c.point_at(0.0);
        assert!(p.length() < 1e-12, "phase-zero at t=0 should be at origin: {p:?}");
        // At t=pi/2, sin(pi/2)=1 → y should equal amplitude.
        let p2 = c.point_at(std::f64::consts::FRAC_PI_2);
        assert!((p2.y - 3.0).abs() < 1e-9, "y at t=pi/2 should be amplitude=3: {p2:?}");
    }

    #[test]
    fn curve3_sine_wave_variant_dispatches_evaluator() {
        let c = Curve3::SineWave(SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 1.0,
            frequency: 2.0,
            phase: 0.0,
        });
        let t = 0.5;
        let p = c.point_at(t);
        let expected = DVec3::new(0.5, (2.0_f64 * t).sin(), 0.0);
        assert!((p - expected).length() < 1e-12);
        // Tangent should be non-zero
        let tan = c.tangent_at(t);
        assert!(tan.length() > 0.9, "tangent should be roughly unit-length: {tan:?}");
    }
}

// 鈹€鈹€ Bezier (de Casteljau) implementations 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// De Casteljau algorithm for rational Bezier curve evaluation in 3D.
/// `t` 鈭?[0, 1].
fn de_casteljau_3d(points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n == 0 {
        return DVec3::ZERO;
    }
    // Work in homogeneous coordinates [x*w, y*w, z*w, w]
    let mut d: Vec<[f64; 4]> = points
        .iter()
        .zip(weights)
        .map(|(p, &w)| [p.x * w, p.y * w, p.z * w, w])
        .collect();
    for r in 1..n {
        for j in 0..n - r {
            let next = d[j + 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(next.iter()) {
                *elem = (1.0 - t) * *elem + t * p;
            }
        }
    }
    let w = d[0][3];
    if w.abs() < 1e-15 {
        DVec3::ZERO
    } else {
        DVec3::new(d[0][0] / w, d[0][1] / w, d[0][2] / w)
    }
}

/// De Casteljau algorithm for rational Bezier curve evaluation in 2D.
fn de_casteljau_2d(points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n == 0 {
        return DVec2::ZERO;
    }
    let mut d: Vec<[f64; 3]> = points
        .iter()
        .zip(weights)
        .map(|(p, &w)| [p.x * w, p.y * w, w])
        .collect();
    for r in 1..n {
        for j in 0..n - r {
            let next = d[j + 1];
            let cur = &mut d[j];
            for (elem, p) in cur.iter_mut().zip(next.iter()) {
                *elem = (1.0 - t) * *elem + t * p;
            }
        }
    }
    let w = d[0][2];
    if w.abs() < 1e-15 {
        DVec2::ZERO
    } else {
        DVec2::new(d[0][0] / w, d[0][1] / w)
    }
}

impl CurveEval for BezierCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        de_casteljau_3d(&self.control_points, &self.weights, t)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        bezier_tangent_analytic(&self.control_points, &self.weights, t).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }
}

impl SurfaceEval for BezierSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let n_u = self.control_points.len();
        if n_u == 0 {
            return DVec3::ZERO;
        }
        let n_v = self.control_points[0].len();
        if n_v == 0 {
            return DVec3::ZERO;
        }
        // Apply de Casteljau in u for each v-column, producing n_v intermediate points
        let row_points: Vec<DVec3> = (0..n_v)
            .map(|j| {
                let col_pts: Vec<DVec3> = (0..n_u).map(|i| self.control_points[i][j]).collect();
                let col_wts: Vec<f64> = (0..n_u).map(|i| self.weights[i][j]).collect();
                de_casteljau_3d(&col_pts, &col_wts, u)
            })
            .collect();
        let unit_wts = vec![1.0; n_v];
        de_casteljau_3d(&row_points, &unit_wts, v)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (self.point_at(u, v + eps) - self.point_at(u, v - eps)) / (2.0 * eps);
        let n = du.cross(dv);
        let len = n.length();
        if len < 1e-15 { DVec3::Z } else { n / len }
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 1.0, 0.0, 1.0]
    }
}

fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |acc, v| acc * v as f64)
}

fn trinomial_coeff(n: usize, i: usize, j: usize, k: usize) -> f64 {
    factorial(n) / (factorial(i) * factorial(j) * factorial(k))
}

impl SurfaceEval for TriBezierSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let degree = self.control_points.len().saturating_sub(1);
        if self.control_points.is_empty() || self.weights.len() != self.control_points.len() {
            return DVec3::ZERO;
        }

        let w = 1.0 - u - v;
        let mut homo = [0.0; 4];
        for (i, row) in self.control_points.iter().enumerate() {
            if row.len() != degree + 1 - i || self.weights.get(i).map(|r| r.len()) != Some(row.len()) {
                return DVec3::ZERO;
            }
            for (j, point) in row.iter().enumerate() {
                let k = degree - i - j;
                let basis = trinomial_coeff(degree, i, j, k)
                    * u.powi(i as i32)
                    * v.powi(j as i32)
                    * w.powi(k as i32);
                let weight = self.weights[i][j];
                homo[0] += basis * weight * point.x;
                homo[1] += basis * weight * point.y;
                homo[2] += basis * weight * point.z;
                homo[3] += basis * weight;
            }
        }

        if homo[3].abs() < 1e-15 {
            DVec3::ZERO
        } else {
            DVec3::new(homo[0] / homo[3], homo[1] / homo[3], homo[2] / homo[3])
        }
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (self.point_at(u, v + eps) - self.point_at(u, v - eps)) / (2.0 * eps);
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        [0.0, 1.0, 0.0, 1.0]
    }
}

impl SurfaceEval for RuledSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let start = self.start.point_at(u);
        let end = self.end.point_at(u);
        start.lerp(end, v)
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = self.end.point_at(u) - self.start.point_at(u);
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        let [u0, u1] = self.start.default_domain();
        [u0, u1, 0.0, 1.0]
    }
}

impl Curve2dEval for BezierCurve2 {
    fn point_at(&self, t: f64) -> DVec2 {
        de_casteljau_2d(&self.control_points, &self.weights, t)
    }
}

// 鈹€鈹€ Offset Curve / Surface implementations 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

impl CurveEval for OffsetCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let base_pt = self.basis.point_at(t);
        let tangent = self.basis.tangent_at(t);
        let perp = tangent.cross(self.offset_dir);
        let perp_len = perp.length();
        if perp_len < 1e-15 {
            return base_pt;
        }
        base_pt + self.offset_distance * (perp / perp_len)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let eps = 1e-6;
        let [t0, t1] = self.basis.default_domain();
        let t_lo = (t - eps).max(t0);
        let t_hi = (t + eps).min(t1);
        let dp = self.point_at(t_hi) - self.point_at(t_lo);
        let len = dp.length();
        if len < 1e-15 { DVec3::X } else { dp / len }
    }
    fn default_domain(&self) -> [f64; 2] {
        self.basis.default_domain()
    }
}

impl SurfaceEval for OffsetSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let base_pt = self.basis.point_at(u, v);
        let n = self.basis.normal_at(u, v);
        base_pt + self.offset_distance * n
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        // Offset preserves the normal direction (first-order approximation)
        self.basis.normal_at(u, v)
    }
    fn default_domain(&self) -> [f64; 4] {
        self.basis.default_domain()
    }
}

#[cfg(test)]
mod eval_tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn line3_point_at() {
        let l = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        assert!((l.point_at(3.0) - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn circle3_point_at_zero_is_on_circle() {
        // Circle in XY plane, normal = Z
        let c = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 2.0,
        };
        let p0 = c.point_at(0.0);
        assert!((p0.length() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn circle3_full_revolution_closes() {
        let c = Circle3 {
            center: DVec3::new(1.0, 2.0, 3.0),
            normal: DVec3::Y,
            radius: 5.0,
        };
        let p0 = c.point_at(0.0);
        let p2pi = c.point_at(2.0 * PI);
        assert!((p0 - p2pi).length() < 1e-10);
    }

    #[test]
    fn circle3_quarter_turn() {
        let c = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        };
        let p0 = c.point_at(0.0);
        let p90 = c.point_at(FRAC_PI_2);
        // 90掳 rotation: p0 and p90 should be perpendicular from center
        assert!((p0.dot(p90)).abs() < 1e-10);
        assert!((p90.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sphere_surface_north_pole() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        };
        // v=0 is north pole regardless of u
        let p = s.point_at(0.0, 0.0);
        // Should be at (0, 3, 0)
        assert!((p - DVec3::new(0.0, 3.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn sphere_surface_point_on_sphere() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        };
        for u in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0] {
            for v in [0.1, 0.5, 1.0, PI / 2.0, PI - 0.1] {
                let p = s.point_at(u, v);
                assert!(
                    (p.length() - 2.0).abs() < 1e-9,
                    "u={u} v={v} |p|={}",
                    p.length()
                );
            }
        }
    }

    #[test]
    fn cylinder_surface_point_on_cylinder() {
        let c = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
        };
        for u in [0.0, 1.0, PI, 2.0 * PI - 0.1] {
            let p = c.point_at(u, 0.0);
            let radial = DVec3::new(p.x, 0.0, p.z).length();
            assert!((radial - 3.0).abs() < 1e-9, "u={u} radial={radial}");
        }
    }

    #[test]
    fn bspline_degree1_linear_interpolation() {
        // Degree-1 BSpline with 2 control points = straight line
        let c = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
        };
        let p0 = c.point_at(0.0);
        let p1 = c.point_at(1.0);
        let pmid = c.point_at(0.5);
        assert!((p0 - DVec3::ZERO).length() < 1e-10);
        assert!((p1 - DVec3::X).length() < 1e-10);
        assert!((pmid - DVec3::new(0.5, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn bspline_degree2_quadratic() {
        // Degree-2 quadratic arc through 3 control points
        let c = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::new(0.5, 1.0, 0.0), DVec3::X],
            weights: vec![1.0, 1.0, 1.0],
        };
        let p0 = c.point_at(0.0);
        let p1 = c.point_at(1.0);
        assert!((p0 - DVec3::ZERO).length() < 1e-10);
        assert!((p1 - DVec3::X).length() < 1e-10);
    }

    #[test]
    fn torus_surface_point_on_torus() {
        let t = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        for u in [0.0, PI / 2.0, PI] {
            for v in [0.0, PI / 2.0, PI] {
                let p = t.point_at(u, v);
                // Distance from the tube center circle should be minor_radius
                let x_ax = any_perpendicular(DVec3::Y);
                let y_ax = DVec3::Y.cross(x_ax).normalize();
                let tube_center = t.center + t.major_radius * (u.cos() * x_ax + u.sin() * y_ax);
                assert!((p - tube_center).length() - 1.0 < 1e-9, "u={u} v={v}");
            }
        }
    }

    #[test]
    fn ellipsoid_surface_satisfies_implicit_equation() {
        let s = EllipsoidalSurface {
            center: DVec3::new(1.0, -2.0, 0.5),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 4.0,
            radius_y: 2.0,
            radius_z: 1.5,
        };
        let p = s.point_at(0.7, 1.2) - s.center;
        let value = (p.x / s.radius_x).powi(2)
            + (p.y / s.radius_y).powi(2)
            + (p.z / s.radius_z).powi(2);
        assert!((value - 1.0).abs() < 1e-9, "implicit value should be 1, got {value}");
    }

    #[test]
    fn ellipsoid_surface_normal_matches_gradient_direction() {
        let s = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };
        let u = 0.9;
        let v = 1.1;
        let p = s.point_at(u, v);
        let expected = DVec3::new(
            p.x / (s.radius_x * s.radius_x),
            p.y / (s.radius_y * s.radius_y),
            p.z / (s.radius_z * s.radius_z),
        )
        .normalize();
        let n = s.normal_at(u, v);
        assert!((n - expected).length() < 1e-9, "n={n:?} expected={expected:?}");
    }

    #[test]
    fn helicoid_surface_advances_by_pitch_per_turn() {
        let s = HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 6.0,
        };
        let p0 = s.point_at(0.0, 2.0);
        let p1 = s.point_at(2.0 * PI, 2.0);
        let delta = p1 - p0;
        assert!((delta - DVec3::new(0.0, 0.0, 6.0)).length() < 1e-9, "delta={delta:?}");
    }

    #[test]
    fn helicoid_surface_normal_is_perpendicular_to_parametric_tangents() {
        let s = HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 4.0,
        };
        let u = 0.6;
        let v = 1.75;
        let n = s.normal_at(u, v);
        let eps = 1e-6;
        let du = (s.point_at(u + eps, v) - s.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (s.point_at(u, v + eps) - s.point_at(u, v - eps)) / (2.0 * eps);
        assert!(n.dot(du).abs() < 1e-6, "n·du={} should be near 0", n.dot(du));
        assert!(n.dot(dv).abs() < 1e-6, "n·dv={} should be near 0", n.dot(dv));
        assert!(n.length() > 0.99, "normal should be unit-length: {n:?}");
    }

    #[test]
    fn pipe_surface_with_line_spine_matches_cylindrical_section() {
        let surface = PipeSurface {
            spine: Box::new(Curve3::Line(Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::Z,
            })),
            ref_dir: DVec3::X,
            radius: 2.0,
        };

        assert!((surface.point_at(0.0, 0.0) - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-9);
        assert!((surface.point_at(PI * 0.5, 0.5) - DVec3::new(0.0, 2.0, 0.5)).length() < 1e-9);
        assert!((surface.default_domain()[0] - 0.0).abs() < 1e-12);
        assert!((surface.default_domain()[1] - 2.0 * PI).abs() < 1e-12);
    }

    #[test]
    fn tri_bezier_surface_hits_triangle_corners() {
        let surface = TriBezierSurface {
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0]],
        };
        assert!((surface.point_at(0.0, 0.0) - DVec3::new(0.0, 0.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(1.0, 0.0) - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(0.0, 1.0) - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn tri_bezier_surface_dispatches_through_surface3() {
        let surface = Surface3::TriBezier(TriBezierSurface {
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0]],
        });
        let p = surface.point_at(0.25, 0.5);
        assert!(p.x >= -1e-12 && p.y >= -1e-12);
        assert!(surface.normal_at(0.2, 0.2).length() > 0.99);
    }

    #[test]
    fn ruled_surface_interpolates_between_curves() {
        let surface = RuledSurface {
            start: Box::new(Curve3::Line(Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::X,
            })),
            end: Box::new(Curve3::Line(Line3 {
                origin: DVec3::Y,
                direction: DVec3::X,
            })),
        };
        assert!((surface.point_at(0.25, 0.0) - DVec3::new(0.25, 0.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(0.25, 1.0) - DVec3::new(0.25, 1.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(0.25, 0.5) - DVec3::new(0.25, 0.5, 0.0)).length() < 1e-12);
        assert!(surface.normal_at(0.25, 0.5).length() > 0.99);
    }

    #[test]
    fn coons_surface_interpolates_all_four_boundaries() {
        let surface = CoonsSurface {
            south: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::X,
            })),
            north: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 1.0, 1.0),
                direction: DVec3::X,
            })),
            west: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::new(0.0, 1.0, 1.0),
            })),
            east: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(1.0, 0.0, 0.0),
                direction: DVec3::new(0.0, 1.0, 1.0),
            })),
        };

        assert!((surface.point_at(0.3, 0.0) - DVec3::new(0.3, 0.0, 0.0)).length() < 1e-9);
        assert!((surface.point_at(0.3, 1.0) - DVec3::new(0.3, 1.0, 1.0)).length() < 1e-9);
        assert!((surface.point_at(0.0, 0.4) - DVec3::new(0.0, 0.4, 0.4)).length() < 1e-9);
        assert!((surface.point_at(1.0, 0.4) - DVec3::new(1.0, 0.4, 0.4)).length() < 1e-9);
        assert!((surface.point_at(0.5, 0.5) - DVec3::new(0.5, 0.5, 0.5)).length() < 1e-9);
    }

    #[test]
    fn conical_surface_uses_slant_distance_from_reference_circle() {
        let surface = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            half_angle_rad: 30.0_f64.to_radians(),
        };

        let p0 = surface.point_at(0.0, 0.0);
        assert!(p0.dot(surface.axis_dir()).abs() < 1e-9);
        assert!((p0.length() - 2.0).abs() < 1e-9);

        let slant = 4.0;
        let p1 = surface.point_at(0.0, slant);
        assert!((p1.z - slant * surface.half_angle_rad.cos()).abs() < 1e-9);
        let radial = p1 - surface.axis_dir() * p1.dot(surface.axis_dir());
        assert!((radial.length() - (2.0 + slant * surface.half_angle_rad.sin())).abs() < 1e-9);
    }

    #[test]
    fn conical_surface_derives_true_apex_from_reference_circle() {
        let surface = ConicalSurface {
            apex: DVec3::new(0.0, 0.0, 5.0),
            axis: DVec3::Z,
            radius: 2.0,
            half_angle_rad: 45.0_f64.to_radians(),
        };

        assert!((surface.apex_point() - DVec3::new(0.0, 0.0, 3.0)).length() < 1e-9);
    }

    // 鈹€鈹€ Analytic derivative tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    /// Quadratic Bezier: P0=(0,0,0), P1=(0.5,1,0), P2=(1,0,0), unit weights.
    /// Analytic tangent at t=0 should be (0.5,1,0).normalize() = (1,2,0)/鈭?.
    #[test]
    fn bezier_tangent_at_endpoint_analytic() {
        let pts = vec![DVec3::ZERO, DVec3::new(0.5, 1.0, 0.0), DVec3::new(1.0, 0.0, 0.0)];
        let wts = vec![1.0, 1.0, 1.0];
        let c = BezierCurve3 { control_points: pts, weights: wts };
        let tan = c.tangent_at(0.0);
        let expected = DVec3::new(1.0, 2.0, 0.0).normalize();
        assert!((tan - expected).length() < 1e-10, "tan={tan:?} expected={expected:?}");
    }

    /// Quadratic Bezier tangent at t=1 should be (1,-2,0)/鈭?.
    #[test]
    fn bezier_tangent_at_end_analytic() {
        let pts = vec![DVec3::ZERO, DVec3::new(0.5, 1.0, 0.0), DVec3::new(1.0, 0.0, 0.0)];
        let wts = vec![1.0, 1.0, 1.0];
        let c = BezierCurve3 { control_points: pts, weights: wts };
        let tan = c.tangent_at(1.0);
        let expected = DVec3::new(1.0, -2.0, 0.0).normalize();
        assert!((tan - expected).length() < 1e-10, "tan={tan:?} expected={expected:?}");
    }

    /// Degree-1 B-Spline (polyline): tangent should be constant along each segment.
    #[test]
    fn bspline_degree1_tangent_is_segment_direction() {
        // Two-segment polyline: (0,0,0)鈫?1,0,0)鈫?1,1,0)
        let pts = vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)];
        let wts = vec![1.0, 1.0, 1.0];
        let knots = vec![0.0, 0.0, 0.5, 1.0, 1.0];
        let c = BSplineCurve3 { degree: 1, knots, control_points: pts, weights: wts };
        let tan0 = c.tangent_at(0.1);
        assert!((tan0 - DVec3::X).length() < 1e-10, "first segment should be +X, got {tan0:?}");
        let tan1 = c.tangent_at(0.9);
        assert!((tan1 - DVec3::Y).length() < 1e-10, "second segment should be +Y, got {tan1:?}");
    }

    /// Degree-2 B-Spline circle arc: tangent should be perpendicular to radius.
    #[test]
    fn bspline_circle_tangent_perpendicular_to_radius() {
        // Use circle_to_bspline to get an exact NURBS circle, then check tangents.
        let circle = Circle3 { center: DVec3::ZERO, normal: DVec3::Z, radius: 1.0 };
        let c = crate::nurbs_convert::circle_to_bspline(&circle);
        for &t in &[0.0, 0.5, 1.0, 1.5, 2.0] {
            let pt = c.point_at(t);
            let tan = c.tangent_at(t);
            // Tangent must be perpendicular to the radius vector
            let dot = pt.normalize_or_zero().dot(tan);
            assert!(dot.abs() < 1e-8, "t={t}: radius路tangent={dot} (should be 0)");
            // Tangent must be a unit vector
            assert!((tan.length() - 1.0).abs() < 1e-10, "t={t}: |tan|={}", tan.length());
        }
    }
}
