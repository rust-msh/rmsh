//! Parametric 2D geometric constraint solver (GCS).
//!
//! # Overview
//!
//! A [`Sketch`] holds a set of 2D geometric entities (points, lines, circles,
//! arcs) and a set of geometric constraints between them.  Calling
//! [`Sketch::solve`] runs a Newton-Raphson iteration that adjusts entity
//! parameters until all constraints are satisfied (or reports failure).
//!
//! # Quick start
//!
//! ```
//! use rcad_constraints::Sketch;
//! use rcad_constraints::constraint::Constraint;
//!
//! let mut sk = Sketch::new();
//! let p1 = sk.add_point(0.0, 0.0);
//! let p2 = sk.add_point(3.0, 4.0);
//!
//! // Fix p1 at the origin
//! sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
//! // Constrain the distance p1→p2 to exactly 5
//! sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
//!
//! let result = sk.solve();
//! assert!(result.converged);
//! ```

pub mod constraint;
pub mod entity;
pub mod sketch;
pub mod solver;
pub mod to_brep;
pub mod space3d;

pub use sketch::Sketch;
pub use entity::{EntityId, EntityKind};
pub use constraint::Constraint;
pub use solver::SolveResult;
