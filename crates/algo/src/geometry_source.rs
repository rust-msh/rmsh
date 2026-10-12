//! Geometry-source abstraction: the dependency-inversion seam through which
//! rmsh meshing algorithms query geometry truth (closest point, normal,
//! curvature) without depending on any CAD kernel.
//!
//! # Layering (charter, ROADMAP §0; adjudication in the landscape doc §C)
//!
//! The trait is DEFINED here in rmsh-algo — a zero-rcad crate — and
//! IMPLEMENTED by rcad-consuming adapters (the `rmsh-bridge` crate wraps the
//! open-source rcad kernel). The dependency direction is always
//! rcad → rmsh: rmsh never implements this trait against a concrete kernel,
//! and the core rmsh crates stay rcad-free. Purely mathematical sources
//! (distance to a point cloud, analytic fields) live in rmsh itself and do
//! not need this trait's kernel half.
//!
//! # v1 surface (minimal, P2 card P2a)
//!
//! Two queries, both optional-answer:
//!
//! - [`GeometrySource::closest_hit`] — the projection of a query point onto
//!   the source's geometry, with the surface normal and surface parameters
//!   at the hit. `None` means "this source has no geometry to answer"; an
//!   empty source is legal, not an error.
//! - [`GeometrySource::principal_curvatures`] — the principal curvature
//!   magnitudes at the hit nearest `query`, ordered `|k1| <= |k2|`. `None`
//!   means "not answerable for that region or surface family" (a kernel that
//!   cannot differentiate a family answers `None` for it — honest absence,
//!   not a zero).
//!
//! The size-field framework (next P2 card: curvature field, distance field,
//! composition, gradient limiting) consumes exactly this surface and nothing
//! else.
//!
//! # Adjudications recorded for the main line
//!
//! - Points are [`glam::DVec3`]: the meshing-facing query layer aligns with
//!   the rcad kernel's math type so adapters convert nothing. The rmsh core
//!   model/geo crates keep nalgebra and are untouched by this (rmsh-io
//!   already uses glam; glam is plain open-source math, not an rcad crate).
//! - `uv` IS carried in v1: (point, normal, uv) is the triple a surface
//!   local-property evaluation answers in one shot (GeomLProp_SLProps
//!   semantics), and the size-field framework needs the hit parameters to
//!   re-query in the neighbourhood of the hit.
//! - Normals are the answering SURFACE's natural normal (the parametrization
//!   normal dU × dV), not the use-oriented normal of a topology occurrence:
//!   a shared face TShape can appear several times with different
//!   orientations, and an occurrence-independent source cannot compose them.
//!   Callers needing face-oriented normals compose the orientation where
//!   they know it.
//! - There is deliberately NO feature-size / edge query in v1. The ROADMAP
//!   names feature-size queries under this trait's umbrella; they land with
//!   the feature-capture card that first needs them, not speculatively.
//!
//! # Contract
//!
//! - `closest_hit` / `principal_curvatures` must be pure (no observable
//!   mutation): size fields call them per candidate point and rely on
//!   repeatable answers.
//! - A source answers in world coordinates (the same frame its geometry is
//!   published in); adapters owning located geometry do the lifting.
//! - Curvature ordering is by MAGNITUDE (`|k1| <= |k2|`), signs preserved:
//!   adaptive sizing keys on the largest bend, which is `k2`.

use glam::DVec3;

/// A hit from a geometry query: the projected point, the surface normal
/// there, and the surface parameters at the hit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeometryHit {
    /// The projected point, in world coordinates.
    pub point: DVec3,
    /// The unit surface normal at the hit (see the module docs: the
    /// parametrization normal, not a use-oriented occurrence normal). A
    /// degenerate parametrization point may answer `DVec3::ZERO`; callers
    /// must tolerate that rather than assume unit length blindly.
    pub normal: DVec3,
    /// The surface parameters `(u, v)` at the hit, in the answering face's
    /// parametrization.
    pub uv: [f64; 2],
}

/// The meshing-oriented geometry abstraction (module docs carry the full
/// layering + contract).
///
/// Implemented by rcad adapters (see the `rmsh-bridge` crate), never by
/// rmsh itself against a concrete kernel.
pub trait GeometrySource {
    /// Closest point on the geometry to `query`, with the surface normal
    /// and parameters at the hit. `None` = the source has no geometry to
    /// answer (an empty source is legal, not an error).
    fn closest_hit(&self, query: DVec3) -> Option<GeometryHit>;

    /// The principal curvature magnitudes at the hit nearest `query`
    /// (`k1`, `k2` by convention `|k1| <= |k2|`); `None` = not answerable.
    fn principal_curvatures(&self, query: DVec3) -> Option<[f64; 2]>;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A known analytic answer: the plane z = 0.
    ///
    /// `closest_hit` returns the query's own (x, y) with z clamped to 0 and
    /// the +Z normal; curvature is exactly `[0, 0]` (defined, not `None` —
    /// this is what the rcad kernel answers for planes and the contract the
    /// curvature field will rely on).
    struct FakePlaneSource;

    impl GeometrySource for FakePlaneSource {
        fn closest_hit(&self, query: DVec3) -> Option<GeometryHit> {
            Some(GeometryHit {
                point: DVec3::new(query.x, query.y, 0.0),
                normal: DVec3::Z,
                uv: [query.x, query.y],
            })
        }

        fn principal_curvatures(&self, _query: DVec3) -> Option<[f64; 2]> {
            Some([0.0, 0.0])
        }
    }

    /// A sphere-approximating fake for the curvature half of the contract:
    /// the sphere of radius `r` at the origin answers the radial projection
    /// with the outward normal and both principal curvatures `1/r`.
    struct FakeSphereSource {
        radius: f64,
    }

    impl GeometrySource for FakeSphereSource {
        fn closest_hit(&self, query: DVec3) -> Option<GeometryHit> {
            let dir = query.normalize_or_zero();
            Some(GeometryHit {
                point: dir * self.radius,
                normal: dir,
                uv: [0.0, 0.0],
            })
        }

        fn principal_curvatures(&self, _query: DVec3) -> Option<[f64; 2]> {
            Some([1.0 / self.radius, 1.0 / self.radius])
        }
    }

    /// An empty source is legal: `None` twice, not an error.
    struct EmptySource;

    impl GeometrySource for EmptySource {
        fn closest_hit(&self, _query: DVec3) -> Option<GeometryHit> {
            None
        }

        fn principal_curvatures(&self, _query: DVec3) -> Option<[f64; 2]> {
            None
        }
    }

    #[test]
    fn fake_plane_closest_hit_is_the_vertical_projection() {
        let hit = FakePlaneSource
            .closest_hit(DVec3::new(0.3, -1.25, 7.0))
            .expect("the plane answers every query");
        assert_eq!(hit.point, DVec3::new(0.3, -1.25, 0.0));
        assert_eq!(hit.normal, DVec3::Z);
        assert_eq!(hit.uv, [0.3, -1.25]);
    }

    #[test]
    fn fake_plane_curvature_is_defined_zero() {
        assert_eq!(
            FakePlaneSource.principal_curvatures(DVec3::new(9.0, 9.0, 9.0)),
            Some([0.0, 0.0])
        );
    }

    #[test]
    fn fake_sphere_hits_radially_with_uniform_curvature() {
        let sphere = FakeSphereSource { radius: 2.0 };
        let hit = sphere
            .closest_hit(DVec3::new(0.0, 0.0, 5.0))
            .expect("the sphere answers every query");
        assert_eq!(hit.point, DVec3::new(0.0, 0.0, 2.0));
        assert_eq!(hit.normal, DVec3::Z);
        assert_eq!(
            sphere.principal_curvatures(hit.point),
            Some([0.5, 0.5]),
            "|k1| <= |k2| with both equal 1/r"
        );
    }

    #[test]
    fn empty_source_answers_none_and_that_is_legal() {
        assert_eq!(EmptySource.closest_hit(DVec3::ONE), None);
        assert_eq!(EmptySource.principal_curvatures(DVec3::ONE), None);
    }

    /// The size-field framework consumes the trait generically (monomorphic
    /// call) and through trait objects (mixed source lists) — pin both call
    /// forms compile and dispatch.
    fn field_probe<G: GeometrySource + ?Sized>(source: &G, query: DVec3) -> Option<f64> {
        let hit = source.closest_hit(query)?;
        let [k1, k2] = source.principal_curvatures(query)?;
        Some(hit.point.distance(query) / (1.0 + k2.abs() + k1.abs()))
    }

    #[test]
    fn trait_consumed_generically_and_as_trait_object() {
        let plane = FakePlaneSource;
        assert!(field_probe(&plane, DVec3::new(1.0, 2.0, 3.0)).is_some());

        let objects: Vec<Box<dyn GeometrySource>> =
            vec![Box::new(FakePlaneSource), Box::new(FakeSphereSource { radius: 2.0 })];
        let via_object = field_probe(objects[1].as_ref(), DVec3::new(0.0, 0.0, 5.0));
        assert!(via_object.is_some());
    }
}
