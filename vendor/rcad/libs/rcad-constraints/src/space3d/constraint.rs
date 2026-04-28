//! 3D geometric constraints and their residual equations.
//!
//! Each [`SpaceConstraint`] variant encodes one or more scalar equations f(x)=0.
//! The [`SpaceConstraint::residuals`] method evaluates those equations given the
//! current parameter vector, and [`SpaceConstraint::equation_count`] returns how
//! many equations the constraint contributes.

use crate::space3d::entity::{SpaceEntity, SpaceEntityId, SpacePointRef};
use crate::solver::SOLVER_NORM_TOL;

/// A geometric constraint between 3D space entities.
#[derive(Debug, Clone)]
pub enum SpaceConstraint {
    // ── Point constraints ────────────────────────────────────────────────────
    /// Two points occupy the same location in 3D. (3 equations)
    Coincident(SpacePointRef, SpacePointRef),

    /// A point is fixed at a specific (x, y, z) location. (3 equations)
    Fixed { point: SpacePointRef, x: f64, y: f64, z: f64 },

    // ── Distance ─────────────────────────────────────────────────────────────
    /// Euclidean distance between two points equals `distance`. (1 equation)
    PointDistance { p1: SpacePointRef, p2: SpacePointRef, distance: f64 },

    /// Length of a SpaceLine entity equals `length`. (1 equation)
    LineLength { line: SpaceEntityId, length: f64 },

    // ── Point-on-entity constraints ──────────────────────────────────────────
    /// A point lies on a line (infinite extension). (2 equations)
    PointOnLine { point: SpacePointRef, line: SpaceEntityId },

    /// A point lies on a plane. (1 equation)
    PointOnPlane { point: SpacePointRef, plane: SpaceEntityId },

    /// A point lies on a sphere surface. (1 equation)
    PointOnSphere { point: SpacePointRef, sphere: SpaceEntityId },

    // ── Line orientation ─────────────────────────────────────────────────────
    /// A SpaceLine is parallel to a given direction vector (dx, dy, dz). (2 equations)
    /// The direction is normalized internally.
    LineParallelToDir { line: SpaceEntityId, dx: f64, dy: f64, dz: f64 },

    /// A SpaceLine is parallel to another SpaceLine. (2 equations)
    LineParallelLine(SpaceEntityId, SpaceEntityId),

    /// A SpaceLine is perpendicular to another SpaceLine. (1 equation)
    LinePerpendicularLine(SpaceEntityId, SpaceEntityId),

    /// Two SpaceLines have equal length. (1 equation)
    EqualLineLength(SpaceEntityId, SpaceEntityId),

    // ── Plane constraints ────────────────────────────────────────────────────
    /// A plane has a fixed normal direction (nx, ny, nz). (2 equations)
    /// The normal is normalized internally.
    PlaneNormal { plane: SpaceEntityId, nx: f64, ny: f64, nz: f64 },

    /// Two planes are parallel (normals are parallel). (2 equations)
    PlaneParallel(SpaceEntityId, SpaceEntityId),

    /// Two planes are perpendicular (normals are orthogonal). (1 equation)
    PlanePerpendicular(SpaceEntityId, SpaceEntityId),

    /// A line lies in a plane (both endpoints on the plane). (2 equations)
    LineInPlane { line: SpaceEntityId, plane: SpaceEntityId },

    /// A line is perpendicular to a plane (line direction parallel to plane normal). (2 equations)
    LinePerpendicularPlane { line: SpaceEntityId, plane: SpaceEntityId },

    /// A line is parallel to a plane (line direction perpendicular to plane normal). (1 equation)
    LineParallelPlane { line: SpaceEntityId, plane: SpaceEntityId },

    // ── Sphere constraints ───────────────────────────────────────────────────
    /// Sphere radius equals a fixed value. (1 equation)
    SphereRadius { sphere: SpaceEntityId, radius: f64 },

    /// Two spheres are tangent (distance between centers = sum/diff of radii). (1 equation)
    SphereTangent { s1: SpaceEntityId, s2: SpaceEntityId, external: bool },

    /// Two spheres have equal radius. (1 equation)
    EqualSphereRadius(SpaceEntityId, SpaceEntityId),

    // ── Angle ────────────────────────────────────────────────────────────────
    /// Angle between two lines equals `angle_rad`. (1 equation)
    AngleLineLine { l1: SpaceEntityId, l2: SpaceEntityId, angle_rad: f64 },

    /// Angle between two planes equals `angle_rad`. (1 equation)
    AnglePlanePlane { p1: SpaceEntityId, p2: SpaceEntityId, angle_rad: f64 },

    /// Angle between a line and a plane equals `angle_rad`. (1 equation)
    AngleLinePlane { line: SpaceEntityId, plane: SpaceEntityId, angle_rad: f64 },

    // ── Distance (entity to entity) ──────────────────────────────────────────
    /// Distance from a point to a plane equals `distance`. (1 equation)
    PointPlaneDistance { point: SpacePointRef, plane: SpaceEntityId, distance: f64 },

    // ── Cylinder constraints ─────────────────────────────────────────────────
    /// Two Cylinder/Sphere entities share the same center. (3 equations)
    Concentric3D(SpaceEntityId, SpaceEntityId),

    /// Cylinder radius equals a fixed value. (1 equation)
    CylinderRadius { cylinder: SpaceEntityId, radius: f64 },

    /// A point lies on a cylinder surface. (1 equation)
    PointOnCylinder { point: SpacePointRef, cylinder: SpaceEntityId },

    /// Two cylinders with parallel axes are tangent. (1 equation)
    CylinderTangent { c1: SpaceEntityId, c2: SpaceEntityId, external: bool },

    /// A plane is tangent to a sphere. (1 equation)
    PlaneTangentToSphere { plane: SpaceEntityId, sphere: SpaceEntityId },
}

impl SpaceConstraint {
    // ── Convenience constructors ─────────────────────────────────────────────

    pub fn fix_point(point: impl Into<SpacePointRef>, x: f64, y: f64, z: f64) -> Self {
        SpaceConstraint::Fixed { point: point.into(), x, y, z }
    }

    pub fn point_distance(p1: impl Into<SpacePointRef>, p2: impl Into<SpacePointRef>, distance: f64) -> Self {
        SpaceConstraint::PointDistance { p1: p1.into(), p2: p2.into(), distance }
    }

    pub fn coincident(p1: impl Into<SpacePointRef>, p2: impl Into<SpacePointRef>) -> Self {
        SpaceConstraint::Coincident(p1.into(), p2.into())
    }

    // ── Equation count ───────────────────────────────────────────────────────

    /// Number of scalar equations this constraint contributes.
    pub fn equation_count(&self) -> usize {
        match self {
            SpaceConstraint::Coincident(..) | SpaceConstraint::Fixed { .. }
            | SpaceConstraint::Concentric3D(..) => 3,
            SpaceConstraint::PointOnLine { .. } | SpaceConstraint::PlaneNormal { .. }
            | SpaceConstraint::PlaneParallel(..) | SpaceConstraint::LineInPlane { .. }
            | SpaceConstraint::LinePerpendicularPlane { .. } => 2,
            _ => 1,
        }
    }

    // ── Residual evaluation ──────────────────────────────────────────────────

    /// Evaluate the constraint residuals f(x) into `out`.
    ///
    /// `out` must have length == `equation_count()`.
    /// Returns `false` if the constraint references an entity of the wrong kind.
    pub fn residuals(&self, params: &[f64], entities: &[SpaceEntity], out: &mut [f64]) -> bool {
        match self {
            // ── Coincident ───────────────────────────────────────────────────
            SpaceConstraint::Coincident(p1, p2) => {
                let (x1, y1, z1) = p1.param_indices(entities);
                let (x2, y2, z2) = p2.param_indices(entities);
                out[0] = params[x1] - params[x2];
                out[1] = params[y1] - params[y2];
                out[2] = params[z1] - params[z2];
            }

            // ── Fixed ────────────────────────────────────────────────────────
            SpaceConstraint::Fixed { point, x, y, z } => {
                let (xi, yi, zi) = point.param_indices(entities);
                out[0] = params[xi] - x;
                out[1] = params[yi] - y;
                out[2] = params[zi] - z;
            }

            // ── PointDistance ────────────────────────────────────────────────
            SpaceConstraint::PointDistance { p1, p2, distance } => {
                let (x1, y1, z1) = p1.param_indices(entities);
                let (x2, y2, z2) = p2.param_indices(entities);
                let dx = params[x1] - params[x2];
                let dy = params[y1] - params[y2];
                let dz = params[z1] - params[z2];
                out[0] = (dx * dx + dy * dy + dz * dz).sqrt() - distance;
            }

            // ── LineLength ───────────────────────────────────────────────────
            SpaceConstraint::LineLength { line, length } => {
                let e = &entities[*line];
                let dx = params[e.param(3)] - params[e.param(0)];
                let dy = params[e.param(4)] - params[e.param(1)];
                let dz = params[e.param(5)] - params[e.param(2)];
                out[0] = (dx * dx + dy * dy + dz * dz).sqrt() - length;
            }

            // ── PointOnLine ──────────────────────────────────────────────────
            SpaceConstraint::PointOnLine { point, line } => {
                let (px, py, pz) = point.param_indices(entities);
                let el = &entities[*line];
                let x1 = params[el.param(0)];
                let y1 = params[el.param(1)];
                let z1 = params[el.param(2)];
                let x2 = params[el.param(3)];
                let y2 = params[el.param(4)];
                let z2 = params[el.param(5)];
                let ldx = x2 - x1;
                let ldy = y2 - y1;
                let ldz = z2 - z1;
                let len = (ldx * ldx + ldy * ldy + ldz * ldz).sqrt().max(SOLVER_NORM_TOL);
                // Cross product of (point - p1) with line direction = 0
                let dx = params[px] - x1;
                let dy = params[py] - y1;
                let dz = params[pz] - z1;
                out[0] = (dy * ldz - dz * ldy) / len;
                out[1] = (dz * ldx - dx * ldz) / len;
            }

            // ── PointOnPlane ─────────────────────────────────────────────────
            SpaceConstraint::PointOnPlane { point, plane } => {
                let (px, py, pz) = point.param_indices(entities);
                let ep = &entities[*plane];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let d = params[ep.param(3)];
                out[0] = nx * params[px] + ny * params[py] + nz * params[pz] - d;
            }

            // ── PointOnSphere ────────────────────────────────────────────────
            SpaceConstraint::PointOnSphere { point, sphere } => {
                let (px, py, pz) = point.param_indices(entities);
                let es = &entities[*sphere];
                let cx = params[es.param(0)];
                let cy = params[es.param(1)];
                let cz = params[es.param(2)];
                let r = params[es.param(3)];
                let dx = params[px] - cx;
                let dy = params[py] - cy;
                let dz = params[pz] - cz;
                out[0] = (dx * dx + dy * dy + dz * dz).sqrt() - r;
            }

            // ── LineParallelToDir ────────────────────────────────────────────
            SpaceConstraint::LineParallelToDir { line, dx, dy, dz } => {
                let e = &entities[*line];
                let ldx = params[e.param(3)] - params[e.param(0)];
                let ldy = params[e.param(4)] - params[e.param(1)];
                let ldz = params[e.param(5)] - params[e.param(2)];
                let len_l = (ldx * ldx + ldy * ldy + ldz * ldz).sqrt().max(SOLVER_NORM_TOL);
                let len_d = (dx * dx + dy * dy + dz * dz).sqrt().max(SOLVER_NORM_TOL);
                // Cross product = 0 (normalized)
                let ndx = dx / len_d;
                let ndy = dy / len_d;
                let ndz = dz / len_d;
                out[0] = (ldy * ndz - ldz * ndy) / len_l;
                out[1] = (ldz * ndx - ldx * ndz) / len_l;
            }

            // ── LineParallelLine ─────────────────────────────────────────────
            SpaceConstraint::LineParallelLine(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(3)] - params[e1.param(0)];
                let dy1 = params[e1.param(4)] - params[e1.param(1)];
                let dz1 = params[e1.param(5)] - params[e1.param(2)];
                let dx2 = params[e2.param(3)] - params[e2.param(0)];
                let dy2 = params[e2.param(4)] - params[e2.param(1)];
                let dz2 = params[e2.param(5)] - params[e2.param(2)];
                let len1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1).sqrt().max(SOLVER_NORM_TOL);
                let len2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2).sqrt().max(SOLVER_NORM_TOL);
                // Cross product = 0 (normalized)
                out[0] = (dy1 * dz2 - dz1 * dy2) / (len1 * len2);
                out[1] = (dz1 * dx2 - dx1 * dz2) / (len1 * len2);
            }

            // ── LinePerpendicularLine ────────────────────────────────────────
            SpaceConstraint::LinePerpendicularLine(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(3)] - params[e1.param(0)];
                let dy1 = params[e1.param(4)] - params[e1.param(1)];
                let dz1 = params[e1.param(5)] - params[e1.param(2)];
                let dx2 = params[e2.param(3)] - params[e2.param(0)];
                let dy2 = params[e2.param(4)] - params[e2.param(1)];
                let dz2 = params[e2.param(5)] - params[e2.param(2)];
                let len1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1).sqrt().max(SOLVER_NORM_TOL);
                let len2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2).sqrt().max(SOLVER_NORM_TOL);
                // Dot product = 0 (normalized)
                out[0] = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (len1 * len2);
            }

            // ── EqualLineLength ──────────────────────────────────────────────
            SpaceConstraint::EqualLineLength(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(3)] - params[e1.param(0)];
                let dy1 = params[e1.param(4)] - params[e1.param(1)];
                let dz1 = params[e1.param(5)] - params[e1.param(2)];
                let dx2 = params[e2.param(3)] - params[e2.param(0)];
                let dy2 = params[e2.param(4)] - params[e2.param(1)];
                let dz2 = params[e2.param(5)] - params[e2.param(2)];
                out[0] = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1) - (dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
            }

            // ── PlaneNormal ──────────────────────────────────────────────────
            SpaceConstraint::PlaneNormal { plane, nx, ny, nz } => {
                let e = &entities[*plane];
                let enx = params[e.param(0)];
                let eny = params[e.param(1)];
                let enz = params[e.param(2)];
                let len = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                let nnx = nx / len;
                let nny = ny / len;
                let nnz = nz / len;
                // Cross product of current normal and target normal = 0
                out[0] = eny * nnz - enz * nny;
                out[1] = enz * nnx - enx * nnz;
            }

            // ── PlaneParallel ────────────────────────────────────────────────
            SpaceConstraint::PlaneParallel(p1, p2) => {
                let e1 = &entities[*p1];
                let e2 = &entities[*p2];
                let nx1 = params[e1.param(0)];
                let ny1 = params[e1.param(1)];
                let nz1 = params[e1.param(2)];
                let nx2 = params[e2.param(0)];
                let ny2 = params[e2.param(1)];
                let nz2 = params[e2.param(2)];
                let len1 = (nx1 * nx1 + ny1 * ny1 + nz1 * nz1).sqrt().max(SOLVER_NORM_TOL);
                let len2 = (nx2 * nx2 + ny2 * ny2 + nz2 * nz2).sqrt().max(SOLVER_NORM_TOL);
                // Cross product = 0 (normalized)
                out[0] = (ny1 * nz2 - nz1 * ny2) / (len1 * len2);
                out[1] = (nz1 * nx2 - nx1 * nz2) / (len1 * len2);
            }

            // ── PlanePerpendicular ───────────────────────────────────────────
            SpaceConstraint::PlanePerpendicular(p1, p2) => {
                let e1 = &entities[*p1];
                let e2 = &entities[*p2];
                let nx1 = params[e1.param(0)];
                let ny1 = params[e1.param(1)];
                let nz1 = params[e1.param(2)];
                let nx2 = params[e2.param(0)];
                let ny2 = params[e2.param(1)];
                let nz2 = params[e2.param(2)];
                let len1 = (nx1 * nx1 + ny1 * ny1 + nz1 * nz1).sqrt().max(SOLVER_NORM_TOL);
                let len2 = (nx2 * nx2 + ny2 * ny2 + nz2 * nz2).sqrt().max(SOLVER_NORM_TOL);
                out[0] = (nx1 * nx2 + ny1 * ny2 + nz1 * nz2) / (len1 * len2);
            }

            // ── LineInPlane ──────────────────────────────────────────────────
            SpaceConstraint::LineInPlane { line, plane } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let d = params[ep.param(3)];
                // Both endpoints must satisfy plane equation
                let x1 = params[el.param(0)];
                let y1 = params[el.param(1)];
                let z1 = params[el.param(2)];
                let x2 = params[el.param(3)];
                let y2 = params[el.param(4)];
                let z2 = params[el.param(5)];
                out[0] = nx * x1 + ny * y1 + nz * z1 - d;
                out[1] = nx * x2 + ny * y2 + nz * z2 - d;
            }

            // ── LinePerpendicularPlane ───────────────────────────────────────
            SpaceConstraint::LinePerpendicularPlane { line, plane } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let ldx = params[el.param(3)] - params[el.param(0)];
                let ldy = params[el.param(4)] - params[el.param(1)];
                let ldz = params[el.param(5)] - params[el.param(2)];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let len_l = (ldx * ldx + ldy * ldy + ldz * ldz).sqrt().max(SOLVER_NORM_TOL);
                let len_n = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                // Cross product of line direction and plane normal = 0
                out[0] = (ldy * nz - ldz * ny) / (len_l * len_n);
                out[1] = (ldz * nx - ldx * nz) / (len_l * len_n);
            }

            // ── LineParallelPlane ────────────────────────────────────────────
            SpaceConstraint::LineParallelPlane { line, plane } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let ldx = params[el.param(3)] - params[el.param(0)];
                let ldy = params[el.param(4)] - params[el.param(1)];
                let ldz = params[el.param(5)] - params[el.param(2)];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let len_l = (ldx * ldx + ldy * ldy + ldz * ldz).sqrt().max(SOLVER_NORM_TOL);
                let len_n = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                // Dot product = 0 (line direction perpendicular to plane normal)
                out[0] = (ldx * nx + ldy * ny + ldz * nz) / (len_l * len_n);
            }

            // ── SphereRadius ─────────────────────────────────────────────────
            SpaceConstraint::SphereRadius { sphere, radius } => {
                let e = &entities[*sphere];
                out[0] = params[e.param(3)] - radius;
            }

            // ── SphereTangent ────────────────────────────────────────────────
            SpaceConstraint::SphereTangent { s1, s2, external } => {
                let e1 = &entities[*s1];
                let e2 = &entities[*s2];
                let cx1 = params[e1.param(0)];
                let cy1 = params[e1.param(1)];
                let cz1 = params[e1.param(2)];
                let r1 = params[e1.param(3)];
                let cx2 = params[e2.param(0)];
                let cy2 = params[e2.param(1)];
                let cz2 = params[e2.param(2)];
                let r2 = params[e2.param(3)];
                let dx = cx2 - cx1;
                let dy = cy2 - cy1;
                let dz = cz2 - cz1;
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if *external {
                    out[0] = dist - (r1 + r2);
                } else {
                    out[0] = dist - (r1 - r2).abs();
                }
            }

            // ── EqualSphereRadius ────────────────────────────────────────────
            SpaceConstraint::EqualSphereRadius(s1, s2) => {
                let e1 = &entities[*s1];
                let e2 = &entities[*s2];
                out[0] = params[e1.param(3)] - params[e2.param(3)];
            }

            // ── AngleLineLine ────────────────────────────────────────────────
            SpaceConstraint::AngleLineLine { l1, l2, angle_rad } => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(3)] - params[e1.param(0)];
                let dy1 = params[e1.param(4)] - params[e1.param(1)];
                let dz1 = params[e1.param(5)] - params[e1.param(2)];
                let dx2 = params[e2.param(3)] - params[e2.param(0)];
                let dy2 = params[e2.param(4)] - params[e2.param(1)];
                let dz2 = params[e2.param(5)] - params[e2.param(2)];
                let len1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1).sqrt().max(SOLVER_NORM_TOL);
                let len2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2).sqrt().max(SOLVER_NORM_TOL);
                let dot = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (len1 * len2);
                // cos(angle) - dot = 0
                out[0] = angle_rad.cos() - dot.clamp(-1.0, 1.0);
            }

            // ── AnglePlanePlane ──────────────────────────────────────────────
            SpaceConstraint::AnglePlanePlane { p1, p2, angle_rad } => {
                let e1 = &entities[*p1];
                let e2 = &entities[*p2];
                let nx1 = params[e1.param(0)];
                let ny1 = params[e1.param(1)];
                let nz1 = params[e1.param(2)];
                let nx2 = params[e2.param(0)];
                let ny2 = params[e2.param(1)];
                let nz2 = params[e2.param(2)];
                let len1 = (nx1 * nx1 + ny1 * ny1 + nz1 * nz1).sqrt().max(SOLVER_NORM_TOL);
                let len2 = (nx2 * nx2 + ny2 * ny2 + nz2 * nz2).sqrt().max(SOLVER_NORM_TOL);
                let dot = (nx1 * nx2 + ny1 * ny2 + nz1 * nz2) / (len1 * len2);
                out[0] = angle_rad.cos() - dot.clamp(-1.0, 1.0);
            }

            // ── AngleLinePlane ───────────────────────────────────────────────
            SpaceConstraint::AngleLinePlane { line, plane, angle_rad } => {
                let el = &entities[*line];
                let ep = &entities[*plane];
                let ldx = params[el.param(3)] - params[el.param(0)];
                let ldy = params[el.param(4)] - params[el.param(1)];
                let ldz = params[el.param(5)] - params[el.param(2)];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let len_l = (ldx * ldx + ldy * ldy + ldz * ldz).sqrt().max(SOLVER_NORM_TOL);
                let len_n = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                // Angle between line and plane: sin(angle) = |dir·n| / (|dir|*|n|)
                let sin_val = (ldx * nx + ldy * ny + ldz * nz).abs() / (len_l * len_n);
                out[0] = angle_rad.sin() - sin_val.clamp(0.0, 1.0);
            }

            // ── PointPlaneDistance ───────────────────────────────────────────
            SpaceConstraint::PointPlaneDistance { point, plane, distance } => {
                let (px, py, pz) = point.param_indices(entities);
                let ep = &entities[*plane];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let d = params[ep.param(3)];
                let len_n = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                let num = (nx * params[px] + ny * params[py] + nz * params[pz] - d).abs();
                out[0] = num / len_n - distance;
            }

            // ── Concentric3D ─────────────────────────────────────────────────
            SpaceConstraint::Concentric3D(e1, e2) => {
                let a = &entities[*e1];
                let b = &entities[*e2];
                out[0] = params[a.param(0)] - params[b.param(0)];
                out[1] = params[a.param(1)] - params[b.param(1)];
                out[2] = params[a.param(2)] - params[b.param(2)];
            }

            // ── CylinderRadius ───────────────────────────────────────────────
            SpaceConstraint::CylinderRadius { cylinder, radius } => {
                let e = &entities[*cylinder];
                out[0] = params[e.param(6)] - radius;
            }

            // ── PointOnCylinder ──────────────────────────────────────────────
            SpaceConstraint::PointOnCylinder { point, cylinder } => {
                let (px, py, pz) = point.param_indices(entities);
                let ec = &entities[*cylinder];
                let cx = params[ec.param(0)];
                let cy = params[ec.param(1)];
                let cz = params[ec.param(2)];
                let ax = params[ec.param(3)];
                let ay = params[ec.param(4)];
                let az = params[ec.param(5)];
                let r = params[ec.param(6)];
                let len_a = (ax * ax + ay * ay + az * az).sqrt().max(SOLVER_NORM_TOL);
                let nax = ax / len_a;
                let nay = ay / len_a;
                let naz = az / len_a;
                // Vector from axis point to point
                let vx = params[px] - cx;
                let vy = params[py] - cy;
                let vz = params[pz] - cz;
                // Project onto axis
                let proj = vx * nax + vy * nay + vz * naz;
                // Perpendicular component
                let perp_x = vx - proj * nax;
                let perp_y = vy - proj * nay;
                let perp_z = vz - proj * naz;
                let dist = (perp_x * perp_x + perp_y * perp_y + perp_z * perp_z).sqrt();
                out[0] = dist - r;
            }

            // ── CylinderTangent ──────────────────────────────────────────────
            SpaceConstraint::CylinderTangent { c1, c2, external } => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                let cx1 = params[e1.param(0)];
                let cy1 = params[e1.param(1)];
                let cz1 = params[e1.param(2)];
                let ax1 = params[e1.param(3)];
                let ay1 = params[e1.param(4)];
                let az1 = params[e1.param(5)];
                let r1 = params[e1.param(6)];
                let cx2 = params[e2.param(0)];
                let cy2 = params[e2.param(1)];
                let cz2 = params[e2.param(2)];
                let _ax2 = params[e2.param(3)];
                let _ay2 = params[e2.param(4)];
                let _az2 = params[e2.param(5)];
                let r2 = params[e2.param(6)];
                // Normalize axis 1
                let len_a1 = (ax1 * ax1 + ay1 * ay1 + az1 * az1).sqrt().max(SOLVER_NORM_TOL);
                let nax1 = ax1 / len_a1;
                let nay1 = ay1 / len_a1;
                let naz1 = az1 / len_a1;
                // Vector between axis points
                let dx = cx2 - cx1;
                let dy = cy2 - cy1;
                let dz = cz2 - cz1;
                // Perpendicular distance from c2 axis point to c1 axis line
                let proj = dx * nax1 + dy * nay1 + dz * naz1;
                let perp_x = dx - proj * nax1;
                let perp_y = dy - proj * nay1;
                let perp_z = dz - proj * naz1;
                let perp_dist = (perp_x * perp_x + perp_y * perp_y + perp_z * perp_z).sqrt();
                if *external {
                    out[0] = perp_dist - (r1 + r2);
                } else {
                    out[0] = perp_dist - (r1 - r2).abs();
                }
            }

            // ── PlaneTangentToSphere ─────────────────────────────────────────
            SpaceConstraint::PlaneTangentToSphere { plane, sphere } => {
                let ep = &entities[*plane];
                let nx = params[ep.param(0)];
                let ny = params[ep.param(1)];
                let nz = params[ep.param(2)];
                let d = params[ep.param(3)];
                let es = &entities[*sphere];
                let cx = params[es.param(0)];
                let cy = params[es.param(1)];
                let cz = params[es.param(2)];
                let r = params[es.param(3)];
                let len_n = (nx * nx + ny * ny + nz * nz).sqrt().max(SOLVER_NORM_TOL);
                let dist = (nx * cx + ny * cy + nz * cz - d).abs() / len_n;
                out[0] = dist - r;
            }
        }
        true
    }
}
