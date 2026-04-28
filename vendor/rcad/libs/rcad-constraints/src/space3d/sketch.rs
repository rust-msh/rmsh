//! 3D parametric sketch — container for space entities and constraints.
//!
//! Analogous to the 2D [`crate::Sketch`] but operates in 3D space with
//! [`SpacePoint`], [`SpaceLine`], [`Plane`], and [`Sphere`] entities.

use crate::space3d::constraint::SpaceConstraint;
use crate::space3d::entity::{SpaceEntity, SpaceEntityId, SpaceEntityKind};
use crate::space3d::solver::{SpaceSolveResult, solve_space};

/// Detailed DOF analysis result for 3D sketches.
#[derive(Debug, Clone)]
pub struct SpaceDofReport {
    /// Number of free (non-fixed) parameters.
    pub free_params: usize,
    /// Total number of constraint equations.
    pub equations: usize,
    /// Net degrees of freedom = free_params - equations.
    pub dof: i64,
    /// Indices of free parameters in the sketch's parameter vector.
    pub free_param_indices: Vec<usize>,
}

/// A 3D geometric constraint sketch.
///
/// Holds a flat parameter vector and a list of 3D entities and constraints.
/// Call [`SpaceSketch::solve`] to run the Newton-Raphson solver.
#[derive(Debug, Clone)]
pub struct SpaceSketch {
    pub params: Vec<f64>,
    pub entities: Vec<SpaceEntity>,
    pub constraints: Vec<SpaceConstraint>,
    fixed: Vec<bool>,
}

impl Default for SpaceSketch {
    fn default() -> Self {
        Self::new()
    }
}

impl SpaceSketch {
    /// Create a new empty 3D sketch.
    pub fn new() -> Self {
        Self {
            params: Vec::new(),
            entities: Vec::new(),
            constraints: Vec::new(),
            fixed: Vec::new(),
        }
    }

    /// Add a 3D point entity at the given coordinates.
    pub fn add_point(&mut self, x: f64, y: f64, z: f64) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        self.params.push(x);
        self.params.push(y);
        self.params.push(z);
        self.fixed.push(false);
        self.fixed.push(false);
        self.fixed.push(false);
        self.entities.push(SpaceEntity::new(SpaceEntityKind::SpacePoint, param_start));
        id
    }

    /// Add a 3D line entity from (x1,y1,z1) to (x2,y2,z2).
    pub fn add_line(&mut self, x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        for &v in &[x1, y1, z1, x2, y2, z2] {
            self.params.push(v);
            self.fixed.push(false);
        }
        self.entities.push(SpaceEntity::new(SpaceEntityKind::SpaceLine, param_start));
        id
    }

    /// Add a plane entity with normal (nx,ny,nz) and distance d from origin.
    pub fn add_plane(&mut self, nx: f64, ny: f64, nz: f64, d: f64) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        for &v in &[nx, ny, nz, d] {
            self.params.push(v);
            self.fixed.push(false);
        }
        self.entities.push(SpaceEntity::new(SpaceEntityKind::Plane, param_start));
        id
    }

    /// Add a sphere entity at (cx,cy,cz) with radius r.
    pub fn add_sphere(&mut self, cx: f64, cy: f64, cz: f64, r: f64) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        for &v in &[cx, cy, cz, r] {
            self.params.push(v);
            self.fixed.push(false);
        }
        self.entities.push(SpaceEntity::new(SpaceEntityKind::Sphere, param_start));
        id
    }

    /// Add a cylinder entity with axis point (cx,cy,cz), axis direction (ax,ay,az), and radius r.
    pub fn add_cylinder(&mut self, cx: f64, cy: f64, cz: f64, ax: f64, ay: f64, az: f64, r: f64) -> SpaceEntityId {
        let id = self.entities.len();
        let param_start = self.params.len();
        for &v in &[cx, cy, cz, ax, ay, az, r] {
            self.params.push(v);
            self.fixed.push(false);
        }
        self.entities.push(SpaceEntity::new(SpaceEntityKind::Cylinder, param_start));
        id
    }

    /// Add a constraint.
    pub fn add_constraint(&mut self, constraint: SpaceConstraint) {
        self.constraints.push(constraint);
    }

    /// Fix a specific parameter so the solver won't change it.
    pub fn fix_param(&mut self, param_idx: usize) {
        if param_idx < self.fixed.len() {
            self.fixed[param_idx] = true;
        }
    }

    /// Fix all parameters of an entity.
    pub fn fix_entity(&mut self, entity_id: SpaceEntityId) {
        if let Some(entity) = self.entities.get(entity_id) {
            let count = entity.kind.param_count();
            for i in 0..count {
                let idx = entity.param_start + i;
                if idx < self.fixed.len() {
                    self.fixed[idx] = true;
                }
            }
        }
    }

    /// Get the coordinates of a SpacePoint entity.
    pub fn point_coords(&self, id: SpaceEntityId) -> Option<(f64, f64, f64)> {
        let e = self.entities.get(id)?;
        if e.kind != SpaceEntityKind::SpacePoint {
            return None;
        }
        Some((
            self.params[e.param(0)],
            self.params[e.param(1)],
            self.params[e.param(2)],
        ))
    }

    /// Get the endpoints of a SpaceLine entity.
    pub fn line_endpoints(&self, id: SpaceEntityId) -> Option<(f64, f64, f64, f64, f64, f64)> {
        let e = self.entities.get(id)?;
        if e.kind != SpaceEntityKind::SpaceLine {
            return None;
        }
        Some((
            self.params[e.param(0)],
            self.params[e.param(1)],
            self.params[e.param(2)],
            self.params[e.param(3)],
            self.params[e.param(4)],
            self.params[e.param(5)],
        ))
    }

    /// Run the Newton-Raphson solver on this sketch.
    pub fn solve(&mut self) -> SpaceSolveResult {
        solve_space(&mut self.params, &self.fixed, &self.entities, &self.constraints)
    }

    /// Degrees of freedom = (free parameters) − (constraint equations).
    pub fn dof(&self) -> i64 {
        let free = self.fixed.iter().filter(|&&f| !f).count() as i64;
        let eqs: i64 = self.constraints.iter().map(|c| c.equation_count() as i64).sum();
        free - eqs
    }

    /// Detailed DOF analysis.
    pub fn dof_analysis(&self) -> SpaceDofReport {
        let free_param_indices: Vec<usize> = (0..self.params.len())
            .filter(|&i| !self.fixed[i])
            .collect();
        let free = free_param_indices.len() as i64;
        let eqs: i64 = self.constraints.iter().map(|c| c.equation_count() as i64).sum();
        SpaceDofReport {
            free_params: free_param_indices.len(),
            equations: eqs as usize,
            dof: free - eqs,
            free_param_indices,
        }
    }
}
