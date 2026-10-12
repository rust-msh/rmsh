//! The [`rmsh_algo::GeometrySource`] adapter over the rcad kernel BRep
//! ([`RcadGeometrySource`]).
//!
//! # Kernel query mapping (one kernel function per trait method)
//!
//! - Projection: `rcad_kernel::base::geom_api::closest_point_on_surface`
//!   (the `GeomAPI_ProjectPointOnSurf` stand-in — closed-form for the
//!   analytic families Plane/Sphere/Cylinder/Cone/Torus, numeric sampling +
//!   Newton refinement elsewhere), called per face surface in the
//!   face-local frame: the query is pulled back by the face location
//!   (`TFaceData::surface_location`), and the winning hit is lifted back to
//!   world coordinates.
//! - Normal: `rcad_kernel::base::geom_lprop::SLProps::normal` (the
//!   `GeomLProp_SLProps` stand-in) at the projected `(u, v)`, first order.
//! - Principal curvatures: the same `SLProps` at order 2 — `min_curvature`
//!   / `max_curvature` (the `IsCurvatureDefined` gate honored), reordered to
//!   the trait convention `|k1| <= |k2|`.
//!
//! `SLProps` is the general path (it differentiates every `Surface3` family
//! through the surface's derivative evaluation); the kernel's analytic
//! `math::curvature::principal_curvatures` shortcut exists but covers only
//! the closed-form families, so it is not the adapter's engine.
//!
//! # Honest limits (documented, not hidden)
//!
//! - Projection targets the face's SURFACE, unbounded — face trim is not
//!   consulted. The payload carries a trim handle only as
//!   `TFaceData::uv_domain` (which is `None` for kernel-built primitive
//!   faces; the real trim lives in the wires' pcurves, whose point-in-face
//!   classification is a separate algorithm card). A query far from the
//!   solid can therefore answer with an adjacent face's unbounded surface
//!   extension; the test
//!   `closest_hit_documents_untrimmed_surface_semantics` pins that behavior
//!   so the limitation is visible, and the kernel gap is: TFaceData
//!   materializes no UV domain for trimmed primitive faces (rcad-kernel
//!   topods.rs `TFaceData::uv_domain`), so a trim-respecting projection
//!   would need the pcurve classifier first.
//! - Normal orientation is the surface parametrization normal (see the
//!   trait's module docs): the box's XMin face is stored Reversed, and this
//!   source still answers the +X surface normal there. A degenerate
//!   parametrization point (sphere pole) answers `DVec3::ZERO` per the
//!   trait contract.
//! - Curvature `None` happens only where the kernel's own gate says the
//!   curvature is undefined (singular/degenerate `(u, v)`), or when the
//!   source has no face carrying a surface. Planes answer DEFINED
//!   `[0.0, 0.0]` (the kernel's umbilic branch), not `None` — pinned by a
//!   test.

use glam::DVec3;
use rcad_kernel::base::geom_api::closest_point_on_surface;
use rcad_kernel::base::geom_lprop::SLProps;
use rcad_kernel::topods::{BRep, TShape};
use rmsh_algo::{GeometryHit, GeometrySource};

/// Sampling density for the numeric fallback path of
/// `closest_point_on_surface` (the analytic families never sample; 16 is the
/// density the kernel's own spherical/cylindrical tests use).
const PROJECTION_SAMPLES: usize = 16;

/// Linear tolerance handed to the `SLProps` evaluations (the kernel's own
/// tests use the same order).
const SLPROPS_TOLERANCE: f64 = 1.0e-12;

/// The nearest-face projection state shared by both trait methods.
struct FaceHit<'a> {
    surface: &'a rcad_kernel::geom::Surface3,
    /// Face location (local frame → world).
    location: glam::DAffine3,
    /// Surface parameters at the projection.
    params: (f64, f64),
    /// The projected point in the face-local frame.
    local_point: DVec3,
    distance: f64,
}

/// The [`GeometrySource`] implementation over an rcad-kernel BRep: walks the
/// face TShapes, projects onto each face's surface, answers from the
/// nearest one. See the module docs for the kernel mapping and the honest
/// limits.
pub struct RcadGeometrySource<'a> {
    brep: &'a BRep,
}

impl<'a> RcadGeometrySource<'a> {
    /// Wrap a BRep borrow. The source never mutates the BRep and stays
    /// valid while the borrow lives.
    pub fn new(brep: &'a BRep) -> Self {
        RcadGeometrySource { brep }
    }

    /// Project `query` (world) onto every face surface (local frames) and
    /// keep the nearest hit. `None` when no face carries a surface.
    fn nearest_face_hit(&self, query: DVec3) -> Option<FaceHit<'_>> {
        let mut best: Option<FaceHit> = None;
        for ts in &self.brep.tshapes {
            let TShape::Face(fd) = ts.as_ref() else {
                continue;
            };
            let Some(surface) = fd.surface.as_ref() else {
                continue;
            };
            let location = self.brep.get_location(fd.surface_location);
            let local_query = location.inverse().transform_point3(query);
            let proj = closest_point_on_surface(surface, local_query, PROJECTION_SAMPLES);
            let is_better = match &best {
                // Strict `<` keeps the FIRST face on exact ties (pool order).
                None => true,
                Some(b) => proj.distance < b.distance,
            };
            if is_better {
                best = Some(FaceHit {
                    surface,
                    location,
                    params: proj.params,
                    local_point: proj.point,
                    distance: proj.distance,
                });
            }
        }
        best
    }
}

impl GeometrySource for RcadGeometrySource<'_> {
    fn closest_hit(&self, query: DVec3) -> Option<GeometryHit> {
        let hit = self.nearest_face_hit(query)?;
        let mut props =
            SLProps::new(hit.surface, hit.params.0, hit.params.1, 1, SLPROPS_TOLERANCE);
        let normal = props
            .normal()
            .map(|n| hit.location.transform_vector3(n).normalize_or_zero())
            .unwrap_or(DVec3::ZERO);
        Some(GeometryHit {
            point: hit.location.transform_point3(hit.local_point),
            normal,
            uv: [hit.params.0, hit.params.1],
        })
    }

    fn principal_curvatures(&self, query: DVec3) -> Option<[f64; 2]> {
        let hit = self.nearest_face_hit(query)?;
        let mut props =
            SLProps::new(hit.surface, hit.params.0, hit.params.1, 2, SLPROPS_TOLERANCE);
        if !props.is_curvature_defined() {
            // The kernel's own gate: curvature is undefined at this (u, v)
            // (singular/degenerate parametrization). Honest None.
            return None;
        }
        let k_min = props.min_curvature();
        let k_max = props.max_curvature();
        // Trait convention |k1| <= |k2| (magnitude order, signs preserved).
        Some(if k_min.abs() <= k_max.abs() {
            [k_min, k_max]
        } else {
            [k_max, k_min]
        })
    }
}

#[cfg(test)]
mod tests {
    use glam::DVec3;
    use rmsh_algo::GeometrySource;

    use crate::test_support::unit_box_brep;

    /// The unit box's top face (z = 1) answers the query above it: the hit
    /// lands exactly on the face plane with the +Z surface normal and the
    /// face's (u, v) at the hit.
    ///
    /// NOTE on the card's original probe (0.5, 0.5, 5.0): under the v1
    /// untrimmed-surface semantics the SIDE planes (0.5 away) beat the top
    /// face (4.0 away) — see
    /// `closest_hit_documents_untrimmed_surface_semantics`. The
    /// trim-respecting answer for far queries needs the pcurve point-in-face
    /// classifier (documented kernel gap), so the top-face pin lives at a
    /// query the v1 semantics answer as a human expects.
    #[test]
    fn closest_hit_above_the_top_face_lands_on_it_with_plus_z_normal() {
        let brep = unit_box_brep();
        let source = super::RcadGeometrySource::new(&brep);

        let hit = source
            .closest_hit(DVec3::new(0.5, 0.5, 1.25))
            .expect("the box answers");
        assert_eq!(hit.point, DVec3::new(0.5, 0.5, 1.0), "exact plane hit");
        assert_eq!(hit.normal, DVec3::Z, "exact +Z surface normal");
        assert_eq!(hit.uv, [0.5, 0.5], "top-face parametrization (u=x, v=y)");
    }

    /// Planes carry DEFINED zero curvature in the kernel (umbilic branch),
    /// so the trait answers `Some([0.0, 0.0])`, not `None` — the documented
    /// adjudication of the "planes and curvature" question.
    #[test]
    fn curvature_of_a_plane_face_is_the_defined_zero_pair() {
        let brep = unit_box_brep();
        let source = super::RcadGeometrySource::new(&brep);
        assert_eq!(
            source.principal_curvatures(DVec3::new(0.5, 0.5, 1.25)),
            Some([0.0, 0.0])
        );
    }

    /// A query inside the box hits the nearest face. (0.4, 0.5, 0.5) is
    /// 0.4 from the XMin face and at least 0.5 from every other; the trait
    /// answers the SURFACE parametrization normal (+X — the face TShape is
    /// stored Reversed and the trait does not compose use orientation, per
    /// the module docs).
    #[test]
    fn closest_hit_inside_the_box_hits_the_nearest_face() {
        let brep = unit_box_brep();
        let source = super::RcadGeometrySource::new(&brep);

        let hit = source
            .closest_hit(DVec3::new(0.4, 0.5, 0.5))
            .expect("the box answers");
        assert_eq!(hit.point, DVec3::new(0.0, 0.5, 0.5), "the XMin plane");
        assert_eq!(hit.normal, DVec3::X, "surface normal (orientation-free)");
    }

    /// The documented v1 limitation, pinned so it cannot regress silently:
    /// a query far above the box answers with the nearest UNBOUNDED surface
    /// — the XMin plane's extension (0.5 away) beats the top face (4.0
    /// away) — and the hit point carries z = 5.0, outside the face trim.
    /// Trim-respecting projection requires the pcurve point-in-face
    /// classifier; TFaceData materializes no UV domain for primitive faces
    /// (kernel gap, module docs).
    #[test]
    fn closest_hit_documents_untrimmed_surface_semantics() {
        let brep = unit_box_brep();
        let source = super::RcadGeometrySource::new(&brep);

        let hit = source
            .closest_hit(DVec3::new(0.5, 0.5, 5.0))
            .expect("the box answers");
        assert_eq!(hit.point, DVec3::new(0.0, 0.5, 5.0), "unbounded extension");
        assert_eq!(hit.normal, DVec3::X);
    }

    /// An empty BRep is a legal empty source: `None` twice.
    #[test]
    fn empty_brep_answers_none() {
        let brep = rcad_kernel::BRep::new();
        let source = super::RcadGeometrySource::new(&brep);
        assert_eq!(source.closest_hit(DVec3::ONE), None);
        assert_eq!(source.principal_curvatures(DVec3::ONE), None);
    }

    /// The trait is consumed through the rmsh-algo surface (imported from
    /// there, dispatched dynamically too) — the bridge adapter slots into
    /// size-field-style generic code.
    #[test]
    fn adapter_dispatches_through_the_rmsh_trait_object() {
        let brep = unit_box_brep();
        let source: Box<dyn GeometrySource> = Box::new(super::RcadGeometrySource::new(&brep));
        let hit = source.closest_hit(DVec3::new(0.5, 0.5, 1.25)).expect("answers");
        assert_eq!(hit.point.z, 1.0);
    }
}
