use super::{BuildError, normalize_vector, validate_point, validate_positive};
use glam::DVec3;
use rcad_kernel::Surface3;
use rcad_kernel::geom::{
    ConicalSurface, CylindricalSurface, Plane, SphericalSurface, ToroidalSurface,
};

pub fn plane(origin: DVec3, normal: DVec3) -> Result<Surface3, BuildError> {
    let origin = validate_point("origin", origin)?;
    let normal = normalize_vector("normal", normal)?;
    Ok(Surface3::Plane(Plane { origin, normal }))
}

pub fn make_plane(origin: DVec3, normal: DVec3) -> Result<Surface3, BuildError> {
    plane(origin, normal)
}

pub fn cylindrical_surface(
    origin: DVec3,
    axis: DVec3,
    radius: f64,
) -> Result<Surface3, BuildError> {
    let origin = validate_point("origin", origin)?;
    let axis = normalize_vector("axis", axis)?;
    let radius = validate_positive("radius", radius)?;
    Ok(Surface3::Cylinder(CylindricalSurface {
        origin,
        axis,
        radius,
    }))
}

pub fn make_cylindrical_surface(
    origin: DVec3,
    axis: DVec3,
    radius: f64,
) -> Result<Surface3, BuildError> {
    cylindrical_surface(origin, axis, radius)
}

pub fn spherical_surface(center: DVec3, radius: f64) -> Result<Surface3, BuildError> {
    let center = validate_point("center", center)?;
    let radius = validate_positive("radius", radius)?;
    Ok(Surface3::Sphere(SphericalSurface {
        center,
        axis: DVec3::Z,
        radius,
    }))
}

pub fn make_spherical_surface(center: DVec3, radius: f64) -> Result<Surface3, BuildError> {
    spherical_surface(center, radius)
}

pub fn conical_surface(
    apex: DVec3,
    axis: DVec3,
    half_angle_rad: f64,
) -> Result<Surface3, BuildError> {
    let apex = validate_point("apex", apex)?;
    let axis = normalize_vector("axis", axis)?;
    let half_angle_rad = validate_positive("half_angle_rad", half_angle_rad)?;
    Ok(Surface3::Cone(ConicalSurface {
        apex,
        axis,
        radius: 0.0,
        half_angle_rad,
    }))
}

pub fn make_conical_surface(
    apex: DVec3,
    axis: DVec3,
    half_angle_rad: f64,
) -> Result<Surface3, BuildError> {
    conical_surface(apex, axis, half_angle_rad)
}

pub fn toroidal_surface(
    center: DVec3,
    axis: DVec3,
    major_radius: f64,
    minor_radius: f64,
) -> Result<Surface3, BuildError> {
    let center = validate_point("center", center)?;
    let axis = normalize_vector("axis", axis)?;
    let major_radius = validate_positive("major_radius", major_radius)?;
    let minor_radius = validate_positive("minor_radius", minor_radius)?;
    Ok(Surface3::Torus(ToroidalSurface {
        center,
        axis,
        major_radius,
        minor_radius,
    }))
}

pub fn make_toroidal_surface(
    center: DVec3,
    axis: DVec3,
    major_radius: f64,
    minor_radius: f64,
) -> Result<Surface3, BuildError> {
    toroidal_surface(center, axis, major_radius, minor_radius)
}
