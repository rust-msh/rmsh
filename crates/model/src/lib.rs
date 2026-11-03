//! rmsh mesh-product model surface.
//!
//! The shared finite-element mesh data model ([`Mesh`], [`Node`],
//! [`Element`], [`ElementType`], the MSH partition semantics) is defined
//! ONCE in `rcad-mesh` and re-exported here, so the mesh product and the
//! geometry kernel import the same types. What stays in this crate is the
//! Gmsh-style geometry topology model ([`GModel`] and its entities), which
//! is rmsh-specific and has no rcad counterpart.

pub mod topology;

pub use rcad_mesh::{
    brep_surface_mesh, Element, ElementType, GhostElement, GhostEntity, Mesh, Node, PartitionInfo,
    PartitionedEntity, PeriodicLink, Point3, Vector3, gmsh_node_count,
};

// Backward compatibility aliases
pub use topology::{TopoEdge, TopoFace, TopoSelection, TopoVertex, TopoVolume, Topology};
pub use topology::{GEdge, GFace, GModel, GRegion, GSelection, GVertex};
