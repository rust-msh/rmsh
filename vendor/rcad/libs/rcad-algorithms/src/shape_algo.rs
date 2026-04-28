//! ShapeAlgo-style additional shape algorithms.
//!
//! Provides utilities for shape analysis and geometry extraction, analogous to
//! OCCT `ShapeAlgo` package. This module includes:
//!
//! - `AlgoContainer`: Container for pluggable shape algorithms
//! - `GetBoxGeometry`: Extract box dimensions from a BRep
//! - `GetCylinderGeometry`: Extract cylinder parameters from a BRep
//! - `GetSphereGeometry`: Extract sphere parameters from a BRep
//! - `GetConeGeometry`: Extract cone parameters from a BRep
//! - `GetTorusGeometry`: Extract torus parameters from a BRep
//! - `IsPrimitive`: Check if a shape matches a primitive type

use glam::DVec3;
use rcad_kernel::geom::{
    ConicalSurface, CylindricalSurface, Plane, SphericalSurface, Surface3, ToroidalSurface,
};
use rcad_kernel::BRep;
use std::collections::HashMap;

// =============================================================================
// Geometry Extraction Structures
// =============================================================================

/// Extracted box geometry parameters.
///
/// Represents an axis-aligned box with origin at the minimum corner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoxGeometry {
    /// Origin (minimum corner) of the box.
    pub origin: DVec3,
    /// Dimension along the X axis.
    pub dx: f64,
    /// Dimension along the Y axis.
    pub dy: f64,
    /// Dimension along the Z axis.
    pub dz: f64,
}

/// Extracted cylinder geometry parameters.
///
/// Represents a cylinder defined by origin (center of bottom), axis, radius, and height.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderGeometry {
    /// Origin at the center of the bottom face.
    pub origin: DVec3,
    /// Cylinder axis direction (normalized).
    pub axis: DVec3,
    /// Cylinder radius.
    pub radius: f64,
    /// Cylinder height along the axis.
    pub height: f64,
}

/// Extracted sphere geometry parameters.
///
/// Represents a sphere defined by center and radius.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereGeometry {
    /// Center of the sphere.
    pub center: DVec3,
    /// Sphere radius.
    pub radius: f64,
}

/// Extracted cone geometry parameters.
///
/// Represents a cone defined by apex, axis, and half-angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConeGeometry {
    /// Apex point of the cone.
    pub apex: DVec3,
    /// Cone axis direction (normalized).
    pub axis: DVec3,
    /// Half-angle of the cone (radians).
    pub angle: f64,
}

/// Extracted torus geometry parameters.
///
/// Represents a torus defined by center, axis, and two radii.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorusGeometry {
    /// Center of the torus.
    pub center: DVec3,
    /// Torus axis direction (normalized).
    pub axis: DVec3,
    /// Distance from center to the center of the tube.
    pub major_radius: f64,
    /// Radius of the tube.
    pub minor_radius: f64,
}

// =============================================================================
// ShapeAlgorithm Trait
// =============================================================================

/// Trait for pluggable shape algorithms.
///
/// Algorithms implementing this trait can be registered with an `AlgoContainer`
/// and executed on BRep shapes.
pub trait ShapeAlgorithm: Send + Sync {
    /// Get the name of this algorithm.
    fn name(&self) -> &str;

    /// Execute the algorithm on the given BRep.
    ///
    /// Returns `true` if the algorithm succeeded, `false` otherwise.
    fn execute(&self, brep: &BRep) -> bool;
}

// =============================================================================
// AlgoContainer
// =============================================================================

/// Container for pluggable shape algorithms.
///
/// Provides a registry for algorithms that can be looked up by name and
/// executed on BRep shapes. Analogous to OCCT `ShapeAlgo_AlgoContainer`.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::shape_algo::{AlgoContainer, ShapeAlgorithm};
///
/// let mut container = AlgoContainer::new();
/// // Algorithms can be added via add_algorithm
/// ```
pub struct AlgoContainer {
    /// Registered algorithms indexed by name.
    algorithms: HashMap<String, Box<dyn ShapeAlgorithm>>,
}

impl AlgoContainer {
    /// Create a new empty algorithm container.
    pub fn new() -> Self {
        Self {
            algorithms: HashMap::new(),
        }
    }

    /// Add an algorithm to the container.
    ///
    /// If an algorithm with the same name already exists, it will be replaced.
    pub fn add_algorithm(&mut self, name: &str, algorithm: Box<dyn ShapeAlgorithm>) {
        self.algorithms.insert(name.to_string(), algorithm);
    }

    /// Get an algorithm by name.
    pub fn get_algorithm(&self, name: &str) -> Option<&dyn ShapeAlgorithm> {
        self.algorithms.get(name).map(|b| b.as_ref())
    }

    /// Check if an algorithm exists.
    pub fn has_algorithm(&self, name: &str) -> bool {
        self.algorithms.contains_key(name)
    }

    /// Remove an algorithm by name.
    pub fn remove_algorithm(&mut self, name: &str) -> bool {
        self.algorithms.remove(name).is_some()
    }

    /// Get the number of registered algorithms.
    pub fn len(&self) -> usize {
        self.algorithms.len()
    }

    /// Check if the container is empty.
    pub fn is_empty(&self) -> bool {
        self.algorithms.is_empty()
    }

    /// Get all algorithm names.
    pub fn algorithm_names(&self) -> Vec<&str> {
        self.algorithms.keys().map(|s| s.as_str()).collect()
    }
}

impl Default for AlgoContainer {
    fn default() -> Self {
        Self::new()
    }
}

// =============================================================================
// Geometry Extraction Functions
// =============================================================================

/// Extract box geometry from a BRep.
///
/// A box is recognized as a solid with 6 planar faces arranged in 3 pairs
/// of parallel faces with appropriate normals.
///
/// Returns `None` if the shape is not a valid box.
pub fn get_box_geometry(brep: &BRep) -> Option<BoxGeometry> {
    // Must have exactly one solid with one shell
    if brep.solids.len() != 1 {
        return None;
    }
    let solid = &brep.solids[0];
    if solid.shells.len() != 1 {
        return None;
    }
    let shell = &solid.shells[0];

    // A box has exactly 6 faces
    if shell.faces.len() != 6 {
        return None;
    }

    // Get all face surfaces
    let surfaces = get_face_surfaces(brep, shell);
    if surfaces.len() != 6 {
        return None;
    }

    // All surfaces must be planes
    let planes: Vec<&Plane> = surfaces
        .iter()
        .filter_map(|s| match s {
            Surface3::Plane(p) => Some(p),
            _ => None,
        })
        .collect();

    if planes.len() != 6 {
        return None;
    }

    // Group planes by normal direction (allowing for opposite directions)
    let mut normal_groups: Vec<(DVec3, Vec<&Plane>)> = Vec::new();

    for plane in &planes {
        let normal = plane.normal.normalize_or_zero();
        let found = normal_groups.iter_mut().find(|(n, _)| {
            let dot = n.dot(normal).abs();
            dot > 0.999 // Nearly parallel or anti-parallel
        });

        if let Some((_, group)) = found {
            group.push(*plane);
        } else {
            normal_groups.push((normal, vec![*plane]));
        }
    }

    // A box has 3 pairs of parallel faces
    if normal_groups.len() != 3 {
        return None;
    }

    for (_, group) in &normal_groups {
        if group.len() != 2 {
            return None;
        }
    }

    // Compute box dimensions and origin
    // The origin is the minimum corner, dimensions are the extent along each axis
    let bbox = brep.bounding_box()?;
    let min_pt = bbox[0];
    let max_pt = bbox[1];

    let dx = (max_pt.x - min_pt.x).abs();
    let dy = (max_pt.y - min_pt.y).abs();
    let dz = (max_pt.z - min_pt.z).abs();

    // Verify the planes correspond to the bounding box
    let tolerance = 1e-6;
    for plane in &planes {
        // Each plane should pass through one of the bbox corners
        let d = plane.normal.dot(plane.origin);
        let d_min = plane.normal.dot(min_pt);
        let d_max = plane.normal.dot(max_pt);

        let near_min = (d - d_min).abs() < tolerance;
        let near_max = (d - d_max).abs() < tolerance;

        if !near_min && !near_max {
            return None;
        }
    }

    Some(BoxGeometry {
        origin: min_pt,
        dx,
        dy,
        dz,
    })
}

/// Extract cylinder geometry from a BRep.
///
/// A cylinder is recognized as a solid with:
/// - One cylindrical lateral face
/// - Two planar end caps (optional for partial cylinders)
///
/// Returns `None` if the shape is not a valid cylinder.
pub fn get_cylinder_geometry(brep: &BRep) -> Option<CylinderGeometry> {
    // Must have exactly one solid with one shell
    if brep.solids.len() != 1 {
        return None;
    }
    let solid = &brep.solids[0];
    if solid.shells.len() != 1 {
        return None;
    }
    let shell = &solid.shells[0];

    // Get all face surfaces
    let surfaces = get_face_surfaces(brep, shell);

    // Find the cylindrical surface
    let mut cyl_surf: Option<&CylindricalSurface> = None;
    for surf in &surfaces {
        if let Surface3::Cylinder(c) = surf {
            if cyl_surf.is_some() {
                // Multiple cylinders - not a simple cylinder
                return None;
            }
            cyl_surf = Some(c);
        }
    }

    let cyl = cyl_surf?;

    // Check that other surfaces are planes (caps)
    for surf in &surfaces {
        match surf {
            Surface3::Cylinder(_) => {}
            Surface3::Plane(_) => {}
            _ => return None, // Non-cylindrical, non-planar surface
        }
    }

    // Compute height from the bounding box along the cylinder axis
    let axis = cyl.axis.normalize_or_zero();
    let bbox = brep.bounding_box()?;
    let min_pt = bbox[0];
    let max_pt = bbox[1];

    // Project bbox corners onto the axis to find extent
    let mut min_proj = f64::INFINITY;
    let mut max_proj = f64::NEG_INFINITY;

    for corner in &[min_pt, max_pt, DVec3::new(min_pt.x, min_pt.y, max_pt.z)] {
        let proj = corner.dot(axis);
        min_proj = min_proj.min(proj);
        max_proj = max_proj.max(proj);
    }

    // Also check all combinations of bbox corners
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                let corner = DVec3::new(
                    if i == 0 { min_pt.x } else { max_pt.x },
                    if j == 0 { min_pt.y } else { max_pt.y },
                    if k == 0 { min_pt.z } else { max_pt.z },
                );
                let proj = corner.dot(axis);
                min_proj = min_proj.min(proj);
                max_proj = max_proj.max(proj);
            }
        }
    }

    let height = (max_proj - min_proj).abs();

    // Determine the actual bottom of the cylinder
    let origin_proj = cyl.origin.dot(axis);
    let origin = if (origin_proj - min_proj).abs() < 1e-6 {
        // Cylinder origin is at the bottom
        cyl.origin
    } else {
        // Compute bottom center
        cyl.origin + axis * (min_proj - origin_proj)
    };

    Some(CylinderGeometry {
        origin,
        axis,
        radius: cyl.radius,
        height,
    })
}

/// Extract sphere geometry from a BRep.
///
/// A sphere is recognized as a solid with a single spherical face.
///
/// Returns `None` if the shape is not a valid sphere.
pub fn get_sphere_geometry(brep: &BRep) -> Option<SphereGeometry> {
    // Must have exactly one solid with one shell
    if brep.solids.len() != 1 {
        return None;
    }
    let solid = &brep.solids[0];
    if solid.shells.len() != 1 {
        return None;
    }
    let shell = &solid.shells[0];

    // Get all face surfaces
    let surfaces = get_face_surfaces(brep, shell);

    // For a sphere, we expect exactly one spherical surface
    let mut sphere_surf: Option<&SphericalSurface> = None;
    for surf in &surfaces {
        match surf {
            Surface3::Sphere(s) => {
                if sphere_surf.is_some() {
                    // Multiple spheres - not a simple sphere
                    return None;
                }
                sphere_surf = Some(s);
            }
            _ => return None, // Non-spherical surface
        }
    }

    let sphere = sphere_surf?;

    Some(SphereGeometry {
        center: sphere.center,
        radius: sphere.radius,
    })
}

/// Extract cone geometry from a BRep.
///
/// A cone is recognized as a solid with:
/// - One conical lateral face
/// - Optionally one planar base cap
///
/// Returns `None` if the shape is not a valid cone.
pub fn get_cone_geometry(brep: &BRep) -> Option<ConeGeometry> {
    // Must have exactly one solid with one shell
    if brep.solids.len() != 1 {
        return None;
    }
    let solid = &brep.solids[0];
    if solid.shells.len() != 1 {
        return None;
    }
    let shell = &solid.shells[0];

    // Get all face surfaces
    let surfaces = get_face_surfaces(brep, shell);

    // Find the conical surface
    let mut cone_surf: Option<&ConicalSurface> = None;
    for surf in &surfaces {
        if let Surface3::Cone(c) = surf {
            if cone_surf.is_some() {
                // Multiple cones - not a simple cone
                return None;
            }
            cone_surf = Some(c);
        }
    }

    let cone = cone_surf?;

    // Check that other surfaces are planes (caps)
    for surf in &surfaces {
        match surf {
            Surface3::Cone(_) => {}
            Surface3::Plane(_) => {}
            _ => return None, // Non-conical, non-planar surface
        }
    }

    // Get the apex and angle from the cone
    let apex = cone.apex_point();
    let axis = cone.axis.normalize_or_zero();
    let angle = cone.half_angle_rad;

    Some(ConeGeometry { apex, axis, angle })
}

/// Extract torus geometry from a BRep.
///
/// A torus is recognized as a solid with a single toroidal face.
///
/// Returns `None` if the shape is not a valid torus.
pub fn get_torus_geometry(brep: &BRep) -> Option<TorusGeometry> {
    // Must have exactly one solid with one shell
    if brep.solids.len() != 1 {
        return None;
    }
    let solid = &brep.solids[0];
    if solid.shells.len() != 1 {
        return None;
    }
    let shell = &solid.shells[0];

    // Get all face surfaces
    let surfaces = get_face_surfaces(brep, shell);

    // For a torus, we expect exactly one toroidal surface
    let mut torus_surf: Option<&ToroidalSurface> = None;
    for surf in &surfaces {
        match surf {
            Surface3::Torus(t) => {
                if torus_surf.is_some() {
                    // Multiple tori - not a simple torus
                    return None;
                }
                torus_surf = Some(t);
            }
            _ => return None, // Non-toroidal surface
        }
    }

    let torus = torus_surf?;

    Some(TorusGeometry {
        center: torus.center,
        axis: torus.axis.normalize_or_zero(),
        major_radius: torus.major_radius,
        minor_radius: torus.minor_radius,
    })
}

// =============================================================================
// Primitive Detection Functions
// =============================================================================

/// Check if a BRep represents a box.
pub fn is_box(brep: &BRep) -> bool {
    get_box_geometry(brep).is_some()
}

/// Check if a BRep represents a cylinder.
pub fn is_cylinder(brep: &BRep) -> bool {
    get_cylinder_geometry(brep).is_some()
}

/// Check if a BRep represents a sphere.
pub fn is_sphere(brep: &BRep) -> bool {
    get_sphere_geometry(brep).is_some()
}

/// Check if a BRep represents a cone.
pub fn is_cone(brep: &BRep) -> bool {
    get_cone_geometry(brep).is_some()
}

/// Check if a BRep represents a torus.
pub fn is_torus(brep: &BRep) -> bool {
    get_torus_geometry(brep).is_some()
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Get all surface references for faces in a shell.
fn get_face_surfaces<'a>(brep: &'a BRep, shell: &rcad_kernel::topology::Shell) -> Vec<&'a Surface3> {
    let mut surfaces = Vec::new();

    // Count faces before this shell
    let mut face_offset = 0;
    for solid in &brep.solids {
        for s in &solid.shells {
            if std::ptr::eq(s, shell) {
                break;
            }
            face_offset += s.faces.len();
        }
    }

    for (i, _face) in shell.faces.iter().enumerate() {
        let face_idx = face_offset + i;
        if let Some(&Some(surf_idx)) = brep.geom.face_surface.get(face_idx) {
            if let Some(surf) = brep.geom.surfaces.get(surf_idx) {
                surfaces.push(surf);
            }
        }
    }

    surfaces
}

/// Check if a BRep has valid geometry data.
fn has_geometry(brep: &BRep) -> bool {
    !brep.geom.surfaces.is_empty() || !brep.geom.curves.is_empty()
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{
        ConicalSurface, CylindricalSurface, Plane, SphericalSurface, ToroidalSurface,
    };
    use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
    use rcad_kernel::{BRep, GeomStore, PrimitiveSolid};
    use std::f64::consts::PI;

    fn make_test_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        })
    }

    fn make_test_sphere() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Sphere { radius: 5.0 })
    }

    fn make_test_cylinder() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 3.0,
            height: 10.0,
        })
    }

    fn make_test_cone() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Cone {
            base_radius: 4.0,
            height: 8.0,
        })
    }

    fn make_test_torus() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 5.0,
            minor_radius: 1.5,
        })
    }

    // ─────────────────────────────────────────────────────────────────────────
    // AlgoContainer Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_algo_container_new() {
        let container = AlgoContainer::new();
        assert!(container.is_empty());
        assert_eq!(container.len(), 0);
    }

    #[test]
    fn test_algo_container_add_algorithm() {
        struct TestAlgo;
        impl ShapeAlgorithm for TestAlgo {
            fn name(&self) -> &str {
                "test"
            }
            fn execute(&self, _brep: &BRep) -> bool {
                true
            }
        }

        let mut container = AlgoContainer::new();
        container.add_algorithm("test", Box::new(TestAlgo));

        assert!(!container.is_empty());
        assert_eq!(container.len(), 1);
        assert!(container.has_algorithm("test"));
    }

    #[test]
    fn test_algo_container_remove_algorithm() {
        struct TestAlgo;
        impl ShapeAlgorithm for TestAlgo {
            fn name(&self) -> &str {
                "test"
            }
            fn execute(&self, _brep: &BRep) -> bool {
                true
            }
        }

        let mut container = AlgoContainer::new();
        container.add_algorithm("test", Box::new(TestAlgo));
        assert!(container.remove_algorithm("test"));
        assert!(container.is_empty());
        assert!(!container.remove_algorithm("nonexistent"));
    }

    #[test]
    fn test_algo_container_algorithm_names() {
        struct AlgoA;
        struct AlgoB;

        impl ShapeAlgorithm for AlgoA {
            fn name(&self) -> &str {
                "a"
            }
            fn execute(&self, _brep: &BRep) -> bool {
                true
            }
        }

        impl ShapeAlgorithm for AlgoB {
            fn name(&self) -> &str {
                "b"
            }
            fn execute(&self, _brep: &BRep) -> bool {
                true
            }
        }

        let mut container = AlgoContainer::new();
        container.add_algorithm("a", Box::new(AlgoA));
        container.add_algorithm("b", Box::new(AlgoB));

        let mut names = container.algorithm_names();
        names.sort();
        assert_eq!(names, vec!["a", "b"]);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Box Geometry Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_get_box_geometry_valid() {
        let brep = make_test_box();
        let geom = get_box_geometry(&brep);

        // Box should be centered at origin
        assert!(geom.is_some());
        let g = geom.unwrap();
        // The box is centered at origin, so dimensions should match
        assert!((g.dx - 2.0).abs() < 1e-6);
        assert!((g.dy - 3.0).abs() < 1e-6);
        assert!((g.dz - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_is_box() {
        let box_brep = make_test_box();
        let sphere_brep = make_test_sphere();

        assert!(is_box(&box_brep));
        assert!(!is_box(&sphere_brep));
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Sphere Geometry Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_get_sphere_geometry_valid() {
        let brep = make_test_sphere();
        let geom = get_sphere_geometry(&brep);

        assert!(geom.is_some());
        let g = geom.unwrap();
        assert_eq!(g.center, DVec3::ZERO);
        assert!((g.radius - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_is_sphere() {
        let sphere_brep = make_test_sphere();
        let box_brep = make_test_box();

        assert!(is_sphere(&sphere_brep));
        assert!(!is_sphere(&box_brep));
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Cylinder Geometry Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_get_cylinder_geometry_valid() {
        let brep = make_test_cylinder();
        let geom = get_cylinder_geometry(&brep);

        assert!(geom.is_some());
        let g = geom.unwrap();
        assert!((g.radius - 3.0).abs() < 1e-6);
        assert!((g.height - 10.0).abs() < 1e-6);
        // Axis should be normalized
        assert!((g.axis.length() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_is_cylinder() {
        let cyl_brep = make_test_cylinder();
        let sphere_brep = make_test_sphere();

        assert!(is_cylinder(&cyl_brep));
        assert!(!is_cylinder(&sphere_brep));
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Cone Geometry Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_get_cone_geometry_valid() {
        let brep = make_test_cone();
        let geom = get_cone_geometry(&brep);

        assert!(geom.is_some());
        let g = geom.unwrap();
        // Axis should be normalized
        assert!((g.axis.length() - 1.0).abs() < 1e-6);
        // Angle should be positive
        assert!(g.angle > 0.0);
    }

    #[test]
    fn test_is_cone() {
        let cone_brep = make_test_cone();
        let sphere_brep = make_test_sphere();

        assert!(is_cone(&cone_brep));
        assert!(!is_cone(&sphere_brep));
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Torus Geometry Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_get_torus_geometry_valid() {
        let brep = make_test_torus();
        let geom = get_torus_geometry(&brep);

        assert!(geom.is_some());
        let g = geom.unwrap();
        assert_eq!(g.center, DVec3::ZERO);
        assert!((g.major_radius - 5.0).abs() < 1e-6);
        assert!((g.minor_radius - 1.5).abs() < 1e-6);
        // Axis should be normalized
        assert!((g.axis.length() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_is_torus() {
        let torus_brep = make_test_torus();
        let sphere_brep = make_test_sphere();

        assert!(is_torus(&torus_brep));
        assert!(!is_torus(&sphere_brep));
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Edge Cases
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_empty_brep() {
        let brep = BRep::new();

        assert!(!is_box(&brep));
        assert!(!is_sphere(&brep));
        assert!(!is_cylinder(&brep));
        assert!(!is_cone(&brep));
        assert!(!is_torus(&brep));
    }

    #[test]
    fn test_multi_solid_brep() {
        let mut brep = BRep::new();
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![] }],
        });
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![] }],
        });

        assert!(!is_box(&brep));
        assert!(!is_sphere(&brep));
    }

    #[test]
    fn test_primitive_detection_comprehensive() {
        let box_brep = make_test_box();
        let sphere_brep = make_test_sphere();
        let cyl_brep = make_test_cylinder();
        let cone_brep = make_test_cone();
        let torus_brep = make_test_torus();

        // Box
        assert!(is_box(&box_brep));
        assert!(!is_sphere(&box_brep));
        assert!(!is_cylinder(&box_brep));
        assert!(!is_cone(&box_brep));
        assert!(!is_torus(&box_brep));

        // Sphere
        assert!(!is_box(&sphere_brep));
        assert!(is_sphere(&sphere_brep));
        assert!(!is_cylinder(&sphere_brep));
        assert!(!is_cone(&sphere_brep));
        assert!(!is_torus(&sphere_brep));

        // Cylinder
        assert!(!is_box(&cyl_brep));
        assert!(!is_sphere(&cyl_brep));
        assert!(is_cylinder(&cyl_brep));
        assert!(!is_cone(&cyl_brep));
        assert!(!is_torus(&cyl_brep));

        // Cone
        assert!(!is_box(&cone_brep));
        assert!(!is_sphere(&cone_brep));
        assert!(!is_cylinder(&cone_brep));
        assert!(is_cone(&cone_brep));
        assert!(!is_torus(&cone_brep));

        // Torus
        assert!(!is_box(&torus_brep));
        assert!(!is_sphere(&torus_brep));
        assert!(!is_cylinder(&torus_brep));
        assert!(!is_cone(&torus_brep));
        assert!(is_torus(&torus_brep));
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Geometry Structure Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_box_geometry_fields() {
        let geom = BoxGeometry {
            origin: DVec3::new(1.0, 2.0, 3.0),
            dx: 4.0,
            dy: 5.0,
            dz: 6.0,
        };

        assert_eq!(geom.origin, DVec3::new(1.0, 2.0, 3.0));
        assert!((geom.dx - 4.0).abs() < 1e-10);
        assert!((geom.dy - 5.0).abs() < 1e-10);
        assert!((geom.dz - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_cylinder_geometry_fields() {
        let geom = CylinderGeometry {
            origin: DVec3::new(0.0, 0.0, 0.0),
            axis: DVec3::Y,
            radius: 5.0,
            height: 10.0,
        };

        assert_eq!(geom.origin, DVec3::ZERO);
        assert_eq!(geom.axis, DVec3::Y);
        assert!((geom.radius - 5.0).abs() < 1e-10);
        assert!((geom.height - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_sphere_geometry_fields() {
        let geom = SphereGeometry {
            center: DVec3::new(1.0, 2.0, 3.0),
            radius: 7.0,
        };

        assert_eq!(geom.center, DVec3::new(1.0, 2.0, 3.0));
        assert!((geom.radius - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_cone_geometry_fields() {
        let geom = ConeGeometry {
            apex: DVec3::new(0.0, 5.0, 0.0),
            axis: DVec3::Y,
            angle: PI / 6.0, // 30 degrees
        };

        assert_eq!(geom.apex, DVec3::new(0.0, 5.0, 0.0));
        assert_eq!(geom.axis, DVec3::Y);
        assert!((geom.angle - PI / 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_torus_geometry_fields() {
        let geom = TorusGeometry {
            center: DVec3::new(1.0, 2.0, 3.0),
            axis: DVec3::Z,
            major_radius: 10.0,
            minor_radius: 2.0,
        };

        assert_eq!(geom.center, DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(geom.axis, DVec3::Z);
        assert!((geom.major_radius - 10.0).abs() < 1e-10);
        assert!((geom.minor_radius - 2.0).abs() < 1e-10);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Clone and PartialEq Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn test_geometry_clone() {
        let geom = SphereGeometry {
            center: DVec3::new(1.0, 2.0, 3.0),
            radius: 5.0,
        };
        let cloned = geom.clone();
        assert_eq!(geom, cloned);
    }

    #[test]
    fn test_geometry_partial_eq() {
        let g1 = BoxGeometry {
            origin: DVec3::ZERO,
            dx: 1.0,
            dy: 2.0,
            dz: 3.0,
        };
        let g2 = BoxGeometry {
            origin: DVec3::ZERO,
            dx: 1.0,
            dy: 2.0,
            dz: 3.0,
        };
        let g3 = BoxGeometry {
            origin: DVec3::ZERO,
            dx: 1.0,
            dy: 2.0,
            dz: 4.0,
        };

        assert_eq!(g1, g2);
        assert_ne!(g1, g3);
    }
}
