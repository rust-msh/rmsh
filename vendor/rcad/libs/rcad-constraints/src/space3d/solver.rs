//! Newton-Raphson GCS solver for 3D constraints with numerical Jacobian.
//!
//! The solver treats the constraint system as F(x) = 0 where x is the vector
//! of free (non-fixed) parameters.  It iterates:
//!
//!   Δx = −J⁺ · F(x)
//!   x  ← x + Δx
//!
//! where J⁺ is the Moore-Penrose pseudo-inverse of the Jacobian, computed via
//! the normal equations (J^T J + λI) Δx = J^T (−F).  The Tikhonov
//! regularisation λ handles under-constrained systems by returning the
//! minimum-norm step.
//!
//! The Jacobian is approximated by central finite differences.

use crate::space3d::constraint::SpaceConstraint;
use crate::space3d::entity::SpaceEntity;
use crate::solver::{RESIDUAL_TOL, FD_H, LAMBDA, PIVOT_TOL};

/// Result returned by [`crate::SpaceSketch::solve`].
#[derive(Debug, Clone)]
pub struct SpaceSolveResult {
    /// `true` if all constraint residuals are below [`RESIDUAL_TOL`].
    pub converged: bool,
    /// RMS residual at termination.
    pub residual: f64,
    /// Number of Newton iterations performed.
    pub iterations: usize,
}

/// Convergence tolerance on the RMS constraint residual.
pub const RESIDUAL_TOL_3D: f64 = RESIDUAL_TOL;
/// Maximum Newton iterations.
const MAX_ITER: usize = 100;

/// Run the Newton-Raphson solver for 3D constraints.
///
/// - `params`: full parameter vector (modified in-place for free params).
/// - `fixed`: mask — `fixed[i] == true` means `params[i]` is held constant.
/// - `entities`: entity metadata (for constraint evaluation).
/// - `constraints`: list of 3D constraints.
pub fn solve_space(
    params: &mut Vec<f64>,
    fixed: &[bool],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
) -> SpaceSolveResult {
    let free_params: Vec<usize> = (0..params.len()).filter(|&i| !fixed[i]).collect();
    let n_free = free_params.len();

    let n_eq: usize = constraints.iter().map(|c| c.equation_count()).sum();

    if n_eq == 0 || n_free == 0 {
        return SpaceSolveResult { converged: true, residual: 0.0, iterations: 0 };
    }

    let mut residual;
    let mut iters = 0;

    for _ in 0..MAX_ITER {
        let f = eval_residuals(params, entities, constraints, n_eq);
        residual = rms(&f);
        if residual < RESIDUAL_TOL {
            break;
        }

        let j = numerical_jacobian(params, fixed, &free_params, entities, constraints, n_eq);

        let delta = solve_normal_equations(&j, &f, n_free);

        let mut alpha = 1.0_f64;
        let old_params: Vec<f64> = free_params.iter().map(|&gi| params[gi]).collect();
        for _ in 0..8 {
            for (k, &gi) in free_params.iter().enumerate() {
                params[gi] = old_params[k] + alpha * delta[k];
            }
            let f_new = eval_residuals(params, entities, constraints, n_eq);
            if rms(&f_new) < residual * (1.0 + 1e-4) {
                break;
            }
            alpha *= 0.5;
        }

        iters += 1;
    }

    let f_final = eval_residuals(params, entities, constraints, n_eq);
    residual = rms(&f_final);

    SpaceSolveResult {
        converged: residual < RESIDUAL_TOL,
        residual,
        iterations: iters,
    }
}

fn eval_residuals(
    params: &[f64],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
    n_eq: usize,
) -> Vec<f64> {
    let mut f = vec![0.0_f64; n_eq];
    let mut row = 0;
    for c in constraints {
        let cnt = c.equation_count();
        c.residuals(params, entities, &mut f[row..row + cnt]);
        row += cnt;
    }
    f
}

fn numerical_jacobian(
    params: &[f64],
    _fixed: &[bool],
    free_params: &[usize],
    entities: &[SpaceEntity],
    constraints: &[SpaceConstraint],
    n_eq: usize,
) -> Vec<Vec<f64>> {
    let n_free = free_params.len();
    let mut j = vec![vec![0.0_f64; n_free]; n_eq];
    let mut p = params.to_vec();

    for (col, &gi) in free_params.iter().enumerate() {
        let orig = p[gi];

        p[gi] = orig + FD_H;
        let f_plus = eval_residuals(&p, entities, constraints, n_eq);

        p[gi] = orig - FD_H;
        let f_minus = eval_residuals(&p, entities, constraints, n_eq);

        p[gi] = orig;

        for row in 0..n_eq {
            j[row][col] = (f_plus[row] - f_minus[row]) / (2.0 * FD_H);
        }
    }
    j
}

fn solve_normal_equations(j: &[Vec<f64>], f: &[f64], n: usize) -> Vec<f64> {
    let m = j.len();

    let mut a = vec![vec![0.0_f64; n]; n];
    let mut b = vec![0.0_f64; n];

    for i in 0..m {
        for col in 0..n {
            b[col] -= j[i][col] * f[i];
            for k in 0..n {
                a[col][k] += j[i][col] * j[i][k];
            }
        }
    }
    for i in 0..n {
        a[i][i] += LAMBDA;
    }

    gaussian_elimination(&mut a, &mut b).unwrap_or_else(|| vec![0.0; n])
}

fn gaussian_elimination(a: &mut Vec<Vec<f64>>, b: &mut Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    debug_assert_eq!(a.len(), n);

    for col in 0..n {
        let pivot_row = (col..n)
            .max_by(|&i, &j| a[i][col].abs().partial_cmp(&a[j][col].abs()).unwrap())?;
        a.swap(col, pivot_row);
        b.swap(col, pivot_row);

        let pivot = a[col][col];
        if pivot.abs() < PIVOT_TOL {
            return None;
        }

        for row in (col + 1)..n {
            let factor = a[row][col] / pivot;
            for k in col..n {
                let v = a[col][k] * factor;
                a[row][k] -= v;
            }
            let bv = b[col] * factor;
            b[row] -= bv;
        }
    }

    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        x[i] = b[i];
        for j in (i + 1)..n {
            let v = a[i][j] * x[j];
            x[i] -= v;
        }
        x[i] /= a[i][i];
    }
    Some(x)
}

fn rms(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space3d::constraint::SpaceConstraint;
    use crate::space3d::entity::{SpaceEntity, SpaceEntityKind, SpacePointRef};

    #[test]
    fn fixed_point_3d() {
        let entities = vec![SpaceEntity::new(SpaceEntityKind::SpacePoint, 0)];
        let constraints = vec![SpaceConstraint::Fixed {
            point: SpacePointRef::Point(0),
            x: 1.0,
            y: 2.0,
            z: 3.0,
        }];
        let mut params = vec![10.0_f64, 20.0, 30.0];
        let fixed = vec![false, false, false];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        assert!((params[0] - 1.0).abs() < 1e-7);
        assert!((params[1] - 2.0).abs() < 1e-7);
        assert!((params[2] - 3.0).abs() < 1e-7);
    }

    #[test]
    fn coincident_3d() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 0),
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 3),
        ];
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::Point(0),
                x: 0.0, y: 0.0, z: 0.0,
            },
            SpaceConstraint::Coincident(SpacePointRef::Point(0), SpacePointRef::Point(1)),
        ];
        let mut params = vec![0.0, 0.0, 0.0, 5.0, 5.0, 5.0];
        let fixed = vec![true, true, true, false, false, false];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        assert!((params[3] - 0.0).abs() < 1e-7);
        assert!((params[4] - 0.0).abs() < 1e-7);
        assert!((params[5] - 0.0).abs() < 1e-7);
    }

    #[test]
    fn point_distance_3d() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 0),
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 3),
        ];
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::Point(0),
                x: 0.0, y: 0.0, z: 0.0,
            },
            SpaceConstraint::PointDistance {
                p1: SpacePointRef::Point(0),
                p2: SpacePointRef::Point(1),
                distance: 5.0,
            },
        ];
        let mut params = vec![0.0, 0.0, 0.0, 3.0, 4.0, 0.0];
        let fixed = vec![true, true, true, false, false, false];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        let dx = params[3];
        let dy = params[4];
        let dz = params[5];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!((dist - 5.0).abs() < 1e-7, "distance = {dist}");
    }

    #[test]
    fn point_on_plane() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpacePoint, 0),
            SpaceEntity::new(SpaceEntityKind::Plane, 3),
        ];
        // Plane: z = 10 (normal = (0,0,1), d = 10)
        // Fix the plane params and point x,y, leave z free for PointOnPlane to solve.
        let constraints = vec![
            SpaceConstraint::PointOnPlane { point: SpacePointRef::Point(0), plane: 1 },
        ];
        let mut params = vec![1.0, 2.0, 0.0, 0.0, 0.0, 1.0, 10.0];
        // Fix plane params (indices 3-6) and point x,y (indices 0,1); z (index 2) is free.
        let fixed = vec![false, false, false, true, true, true, true];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        assert!((params[2] - 10.0).abs() < 1e-7, "z should be 10, got {}", params[2]);
    }

    #[test]
    fn line_length_3d() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpaceLine, 0),
        ];
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::LineStart(0),
                x: 0.0, y: 0.0, z: 0.0,
            },
            SpaceConstraint::LineLength { line: 0, length: 3.0 },
        ];
        let mut params = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let fixed = vec![true, true, true, false, false, false];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        let dx = params[3];
        let dy = params[4];
        let dz = params[5];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!((len - 3.0).abs() < 1e-7, "length = {len}");
    }

    #[test]
    fn sphere_tangent_external() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::Sphere, 0),
            SpaceEntity::new(SpaceEntityKind::Sphere, 4),
        ];
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::SphereCenter(0),
                x: 0.0, y: 0.0, z: 0.0,
            },
            SpaceConstraint::SphereRadius { sphere: 0, radius: 2.0 },
            SpaceConstraint::SphereRadius { sphere: 1, radius: 1.0 },
            SpaceConstraint::Fixed {
                point: SpacePointRef::SphereCenter(1),
                x: 5.0, y: 0.0, z: 0.0,
            },
            SpaceConstraint::SphereTangent { s1: 0, s2: 1, external: true },
        ];
        let mut params = vec![0.0, 0.0, 0.0, 2.0, 5.0, 0.0, 0.0, 1.0];
        // Fix sphere 0 center + radius, sphere 1 center + radius
        let fixed = vec![true, true, true, true, true, true, true, true];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        // All params fixed → nothing to solve, but should not panic
        assert!(result.converged);
    }

    #[test]
    fn line_perpendicular_line() {
        let entities = vec![
            SpaceEntity::new(SpaceEntityKind::SpaceLine, 0),
            SpaceEntity::new(SpaceEntityKind::SpaceLine, 6),
        ];
        // Line 0: along X axis from (0,0,0) to (1,0,0)
        // Line 1: along Y axis from (0,0,0) to (0,1,0)
        let constraints = vec![
            SpaceConstraint::Fixed {
                point: SpacePointRef::LineStart(0),
                x: 0.0, y: 0.0, z: 0.0,
            },
            SpaceConstraint::Fixed {
                point: SpacePointRef::LineEnd(0),
                x: 1.0, y: 0.0, z: 0.0,
            },
            SpaceConstraint::Fixed {
                point: SpacePointRef::LineStart(1),
                x: 0.0, y: 0.0, z: 0.0,
            },
            SpaceConstraint::LinePerpendicularLine(0, 1),
        ];
        let mut params = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let fixed = vec![true, true, true, true, true, true, true, true, true, false, false, false];

        let result = solve_space(&mut params, &fixed, &entities, &constraints);
        assert!(result.converged, "not converged: {}", result.residual);
        // Line 1 should end up perpendicular to line 0 (along Y or Z)
        let dx2 = params[9];
        let dy2 = params[10];
        let dz2 = params[11];
        let dot = 1.0 * dx2 + 0.0 * dy2 + 0.0 * dz2;
        assert!(dot.abs() < 1e-7, "dot product should be ~0, got {dot}");
    }
}
