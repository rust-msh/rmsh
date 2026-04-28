//! 3D geometric constraint solver.
//!
//! Analogous to the 2D [`crate::Sketch`] but operates in 3D space with
//! [`SpacePoint`], [`SpaceLine`], [`Plane`], and [`Sphere`] entities.
//!
//! # Quick start
//!
//! ```
//! use rcad_constraints::space3d::SpaceSketch;
//! use rcad_constraints::space3d::constraint::SpaceConstraint;
//!
//! let mut sk = SpaceSketch::new();
//! let p1 = sk.add_point(0.0, 0.0, 0.0);
//! let p2 = sk.add_point(3.0, 4.0, 0.0);
//!
//! // Fix p1 at the origin
//! sk.add_constraint(SpaceConstraint::fix_point(p1, 0.0, 0.0, 0.0));
//! // Constrain the distance p1→p2 to exactly 5
//! sk.add_constraint(SpaceConstraint::point_distance(p1, p2, 5.0));
//!
//! let result = sk.solve();
//! assert!(result.converged);
//! ```

pub mod constraint;
pub mod entity;
pub mod sketch;
pub mod solver;

pub use sketch::SpaceSketch;
pub use entity::{SpaceEntityId, SpaceEntityKind, SpacePointRef};
pub use constraint::SpaceConstraint;
pub use solver::SpaceSolveResult;
