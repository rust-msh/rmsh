use super::{BuildError, normalize_rejection, normalize_vector, validate_point};
use glam::DVec3;
use rcad_kernel::Curve3;
use rcad_kernel::geom::{Circle3, Ellipse3, Line3};

pub fn line(origin: DVec3, direction: DVec3) -> Result<Curve3, BuildError> {
    let origin = validate_point("origin", origin)?;
    let direction = normalize_vector("direction", direction)?;
    Ok(Curve3::Line(Line3 { origin, direction }))
}

pub fn make_line(origin: DVec3, direction: DVec3) -> Result<Curve3, BuildError> {
    line(origin, direction)
}

pub fn circle(center: DVec3, normal: DVec3, radius: f64) -> Result<Curve3, BuildError> {
    let center = validate_point("center", center)?;
    let normal = normalize_vector("normal", normal)?;
    let radius = super::validate_positive("radius", radius)?;
    Ok(Curve3::Circle(Circle3 {
        center,
        normal,
        radius,
    }))
}

pub fn make_circle(center: DVec3, normal: DVec3, radius: f64) -> Result<Curve3, BuildError> {
    circle(center, normal, radius)
}

pub fn ellipse(
    center: DVec3,
    normal: DVec3,
    major_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
) -> Result<Curve3, BuildError> {
    let center = validate_point("center", center)?;
    let normal = normalize_vector("normal", normal)?;
    let major_dir = normalize_rejection("major_dir", major_dir, "normal", normal)?;
    let major_radius = super::validate_positive("major_radius", major_radius)?;
    let minor_radius = super::validate_positive("minor_radius", minor_radius)?;
    Ok(Curve3::Ellipse(Ellipse3 {
        center,
        normal,
        major_dir,
        major_radius,
        minor_radius,
    }))
}

pub fn make_ellipse(
    center: DVec3,
    normal: DVec3,
    major_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
) -> Result<Curve3, BuildError> {
    ellipse(center, normal, major_dir, major_radius, minor_radius)
}
