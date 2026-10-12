//! rmsh bridge to the open-source rcad kernel: the BRep tessellation walk
//! ([`brep_surface_mesh`]) and the [`rmsh_algo::GeometrySource`] adapter
//! ([`RcadGeometrySource`]).
//!
//! # Bridge adoption (P2 card P2a, ROADMAP §2)
//!
//! The behavioral reference for [`brep_surface_mesh`] was the pro-side
//! implementation (`rcad-pro/libs/rcad-io/src/brep_to_mesh.rs`, the
//! `brep_to_surface_mesh` walk over the `BRep::triangulations` pool + face
//! TShapes). This crate RE-AUTHORS that walk as open-source rmsh code, and
//! the pro side is scheduled to consume THIS bridge afterward (the pro
//! rcad-io facade retires in favor of rmsh's format/bridge layer, the same
//! adjudication that retired the earlier pro msh facade). Two deliberate
//! departures from the pro reference are documented where they happen: the
//! triangulation-pool indexing (module `brep_surface_mesh`) and the node-id
//! policy (same module).
//!
//! # Layering (charter, ROADMAP §0)
//!
//! Dependency direction is rcad → rmsh, never the reverse: this crate
//! depends on rmsh-algo (the [`rmsh_algo::GeometrySource`] trait),
//! rmsh-model (the mesh data model) and rcad-kernel (the geometry); rmsh-algo
//! itself stays zero-rcad. This crate is the ONLY sanctioned place where an
//! rmsh-adjacent implementation touches rcad types.
//!
//! # Tessellation entry point (tests)
//!
//! The walk reads whatever triangulations the BRep already carries and never
//! meshes. The tests mount triangulations through the open-source path
//! `rcad_mesh::b_rep_mesh::facade::IncrementalMesh` (the OCCT
//! `BRepMesh_IncrementalMesh` facade — the `incmesh shape defl` DRAW
//! mapping) over primitives built with `rcad_modeling::make_box_brep`; both
//! are dev-dependencies only.

pub mod brep_surface_mesh;
pub mod geometry_source;

#[cfg(test)]
pub(crate) mod test_support;

pub use brep_surface_mesh::brep_surface_mesh;
pub use geometry_source::RcadGeometrySource;
