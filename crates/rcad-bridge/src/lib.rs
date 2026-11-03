//! rmsh bridge to the open-source rcad kernel: the [`rmsh_algo::GeometrySource`]
//! adapter ([`RcadGeometrySource`]).
//!
//! # Layering (charter, ROADMAP §0)
//!
//! The dependency direction is rmsh → rcad: this crate depends on rmsh-algo
//! (the [`rmsh_algo::GeometrySource`] trait), rmsh-model (the mesh data model,
//! re-exported from `rcad-mesh`) and rcad-kernel (the geometry); rmsh-algo
//! itself stays zero-rcad. This crate is the ONLY sanctioned place where an
//! rmsh-adjacent implementation touches rcad types.
//!
//! The BRep → Mesh tessellation walk now lives in `rcad-mesh` (the shared mesh
//! data model crate), which both the kernel side and the mesh product consume.

pub mod geometry_source;

#[cfg(test)]
mod brep_walk_consumability;
#[cfg(test)]
pub(crate) mod test_support;

pub use geometry_source::RcadGeometrySource;
