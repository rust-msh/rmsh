//! Geometric constraints and their residual equations.
//!
//! Each [`Constraint`] variant encodes one or more scalar equations f(x)=0.
//! The [`Constraint::residuals`] method evaluates those equations given the
//! current parameter vector, and [`Constraint::equation_count`] returns how
//! many equations the constraint contributes.

use crate::entity::{Entity, EntityId, PointRef};
use crate::solver::SOLVER_NORM_TOL;

/// A geometric constraint between sketch entities.
#[derive(Debug, Clone)]
pub enum Constraint {
    // ── Point constraints ────────────────────────────────────────────────────
    /// Two points occupy the same location. (2 equations)
    Coincident(PointRef, PointRef),

    /// A point is fixed at a specific (x, y) location. (2 equations)
    Fixed { point: PointRef, x: f64, y: f64 },

    // ── Distance / length ────────────────────────────────────────────────────
    /// Euclidean distance between two points equals `distance`. (1 equation)
    PointDistance { p1: PointRef, p2: PointRef, distance: f64 },

    /// Length of a Line entity equals `length`. (1 equation)
    LineLength { line: EntityId, length: f64 },

    // ── Line orientation ─────────────────────────────────────────────────────
    /// A Line entity is horizontal (y1 == y2). (1 equation)
    Horizontal(EntityId),

    /// A Line entity is vertical (x1 == x2). (1 equation)
    Vertical(EntityId),

    // ── Line–line relations ──────────────────────────────────────────────────
    /// Two Line entities are parallel. (1 equation)
    Parallel(EntityId, EntityId),

    /// Two Line entities are perpendicular. (1 equation)
    Perpendicular(EntityId, EntityId),

    /// Two Line entities have equal length. (1 equation)
    EqualLength(EntityId, EntityId),

    /// Angle from line `l1` to line `l2` equals `angle_rad`. (1 equation)
    Angle { l1: EntityId, l2: EntityId, angle_rad: f64 },

    // ── Circle / arc constraints ─────────────────────────────────────────────
    /// Two Circle or Arc entities have equal radius. (1 equation)
    EqualRadius(EntityId, EntityId),

    /// A point lies on a circle (or arc circle). (1 equation)
    PointOnCircle { point: PointRef, circle: EntityId },

    /// A point lies on a line (infinite extension). (1 equation)
    PointOnLine { point: PointRef, line: EntityId },

    /// A circle is tangent to a line. (1 equation)
    Tangent { circle: EntityId, line: EntityId },

    /// Two circles are externally or internally tangent. (1 equation)
    ///
    /// External tangency: dist(c1, c2) = r1 + r2.
    /// Internal tangency: dist(c1, c2) = |r1 - r2|.
    /// The solver finds whichever solution is closest to the initial guess.
    /// Use `external: true` for external tangency, `false` for internal.
    CircleCircleTangent { c1: EntityId, c2: EntityId, external: bool },

    /// Circle radius equals a fixed value. (1 equation)
    Radius { circle: EntityId, radius: f64 },

    /// Two arcs (or arc and circle) are tangent. (1 equation)
    ///
    /// Uses the same residual as `CircleCircleTangent` — arcs share the same
    /// center/radius parameter layout as circles.
    ArcArcTangent { a1: EntityId, a2: EntityId, external: bool },

    /// Two points are symmetric about a line. (2 equations)
    ///
    /// The midpoint of p1–p2 lies on the line, and the segment p1–p2 is
    /// perpendicular to the line direction.
    Symmetric { p1: PointRef, p2: PointRef, line: EntityId },

    // ── Circle / arc relations ────────────────────────────────────────────────
    /// Two Circle/Arc entities share the same center. (2 equations)
    Concentric(EntityId, EntityId),

    /// A point is the midpoint of a Line segment. (2 equations)
    Midpoint { point: PointRef, line: EntityId },

    /// Circle/Arc diameter equals a fixed value. (1 equation)
    Diameter { circle: EntityId, diameter: f64 },
}

impl Constraint {
    // ── Convenience constructors ─────────────────────────────────────────────

    pub fn fix_point(point: impl Into<PointRef>, x: f64, y: f64) -> Self {
        Constraint::Fixed { point: point.into(), x, y }
    }

    pub fn point_distance(p1: impl Into<PointRef>, p2: impl Into<PointRef>, distance: f64) -> Self {
        Constraint::PointDistance { p1: p1.into(), p2: p2.into(), distance }
    }

    pub fn coincident(p1: impl Into<PointRef>, p2: impl Into<PointRef>) -> Self {
        Constraint::Coincident(p1.into(), p2.into())
    }

    // ── Equation count ───────────────────────────────────────────────────────

    /// Number of scalar equations this constraint contributes.
    pub fn equation_count(&self) -> usize {
        match self {
            Constraint::Coincident(..) | Constraint::Fixed { .. } | Constraint::Symmetric { .. }
            | Constraint::Concentric(..) | Constraint::Midpoint { .. } => 2,
            _ => 1,
        }
    }

    // ── Residual evaluation ──────────────────────────────────────────────────

    /// Evaluate the constraint residuals f(x) into `out`.
    ///
    /// `out` must have length == `equation_count()`.
    /// Returns `false` if the constraint references an entity of the wrong kind.
    pub fn residuals(&self, params: &[f64], entities: &[Entity], out: &mut [f64]) -> bool {
        match self {
            // ── Coincident ───────────────────────────────────────────────────
            Constraint::Coincident(p1, p2) => {
                let (x1, y1) = p1.param_indices(entities);
                let (x2, y2) = p2.param_indices(entities);
                out[0] = params[x1] - params[x2];
                out[1] = params[y1] - params[y2];
            }

            // ── Fixed ────────────────────────────────────────────────────────
            Constraint::Fixed { point, x, y } => {
                let (xi, yi) = point.param_indices(entities);
                out[0] = params[xi] - x;
                out[1] = params[yi] - y;
            }

            // ── PointDistance ────────────────────────────────────────────────
            Constraint::PointDistance { p1, p2, distance } => {
                let (x1, y1) = p1.param_indices(entities);
                let (x2, y2) = p2.param_indices(entities);
                let dx = params[x1] - params[x2];
                let dy = params[y1] - params[y2];
                // Use sqrt form for better conditioning near solution.
                out[0] = (dx * dx + dy * dy).sqrt() - distance;
            }

            // ── LineLength ───────────────────────────────────────────────────
            Constraint::LineLength { line, length } => {
                let e = &entities[*line];
                let dx = params[e.param(2)] - params[e.param(0)];
                let dy = params[e.param(3)] - params[e.param(1)];
                out[0] = (dx * dx + dy * dy).sqrt() - length;
            }

            // ── Horizontal ───────────────────────────────────────────────────
            Constraint::Horizontal(id) => {
                let e = &entities[*id];
                out[0] = params[e.param(1)] - params[e.param(3)]; // y1 - y2
            }

            // ── Vertical ─────────────────────────────────────────────────────
            Constraint::Vertical(id) => {
                let e = &entities[*id];
                out[0] = params[e.param(0)] - params[e.param(2)]; // x1 - x2
            }

            // ── Parallel ─────────────────────────────────────────────────────
            Constraint::Parallel(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(2)] - params[e1.param(0)];
                let dy1 = params[e1.param(3)] - params[e1.param(1)];
                let dx2 = params[e2.param(2)] - params[e2.param(0)];
                let dy2 = params[e2.param(3)] - params[e2.param(1)];
                // cross product of direction vectors = 0
                out[0] = dx1 * dy2 - dy1 * dx2;
            }

            // ── Perpendicular ─────────────────────────────────────────────────
            Constraint::Perpendicular(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(2)] - params[e1.param(0)];
                let dy1 = params[e1.param(3)] - params[e1.param(1)];
                let dx2 = params[e2.param(2)] - params[e2.param(0)];
                let dy2 = params[e2.param(3)] - params[e2.param(1)];
                // dot product = 0
                out[0] = dx1 * dx2 + dy1 * dy2;
            }

            // ── EqualLength ───────────────────────────────────────────────────
            Constraint::EqualLength(l1, l2) => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(2)] - params[e1.param(0)];
                let dy1 = params[e1.param(3)] - params[e1.param(1)];
                let dx2 = params[e2.param(2)] - params[e2.param(0)];
                let dy2 = params[e2.param(3)] - params[e2.param(1)];
                out[0] = (dx1 * dx1 + dy1 * dy1) - (dx2 * dx2 + dy2 * dy2);
            }

            // ── Angle ─────────────────────────────────────────────────────────
            Constraint::Angle { l1, l2, angle_rad } => {
                let e1 = &entities[*l1];
                let e2 = &entities[*l2];
                let dx1 = params[e1.param(2)] - params[e1.param(0)];
                let dy1 = params[e1.param(3)] - params[e1.param(1)];
                let dx2 = params[e2.param(2)] - params[e2.param(0)];
                let dy2 = params[e2.param(3)] - params[e2.param(1)];
                // sin(angle)*(d1·d2) - cos(angle)*(d1×d2) = 0
                let dot = dx1 * dx2 + dy1 * dy2;
                let cross = dx1 * dy2 - dy1 * dx2;
                out[0] = angle_rad.sin() * dot - angle_rad.cos() * cross;
            }

            // ── EqualRadius ───────────────────────────────────────────────────
            Constraint::EqualRadius(c1, c2) => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                // radius is param[2] for both Circle and Arc
                out[0] = params[e1.param(2)] - params[e2.param(2)];
            }

            // ── Radius ────────────────────────────────────────────────────────
            Constraint::Radius { circle, radius } => {
                let e = &entities[*circle];
                out[0] = params[e.param(2)] - radius;
            }

            // ── PointOnCircle ─────────────────────────────────────────────────
            Constraint::PointOnCircle { point, circle } => {
                let (px, py) = point.param_indices(entities);
                let ec = &entities[*circle];
                let cx = params[ec.param(0)];
                let cy = params[ec.param(1)];
                let r = params[ec.param(2)];
                let dx = params[px] - cx;
                let dy = params[py] - cy;
                out[0] = (dx * dx + dy * dy).sqrt() - r;
            }

            // ── PointOnLine ───────────────────────────────────────────────────
            Constraint::PointOnLine { point, line } => {
                let (px, py) = point.param_indices(entities);
                let el = &entities[*line];
                let x1 = params[el.param(0)];
                let y1 = params[el.param(1)];
                let x2 = params[el.param(2)];
                let y2 = params[el.param(3)];
                // Normalize by line length to improve conditioning.
                let len = ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt().max(SOLVER_NORM_TOL);
                out[0] = ((x2 - x1) * (params[py] - y1) - (y2 - y1) * (params[px] - x1)) / len;
            }

            // ── Tangent (circle–line) ─────────────────────────────────────────
            Constraint::Tangent { circle, line } => {
                let ec = &entities[*circle];
                let el = &entities[*line];
                let cx = params[ec.param(0)];
                let cy = params[ec.param(1)];
                let r = params[ec.param(2)];
                let x1 = params[el.param(0)];
                let y1 = params[el.param(1)];
                let x2 = params[el.param(2)];
                let y2 = params[el.param(3)];
                let dx = x2 - x1;
                let dy = y2 - y1;
                let len = (dx * dx + dy * dy).sqrt().max(SOLVER_NORM_TOL);
                // dist(center, line) - r = 0
                let num = (cx - x1) * dy - (cy - y1) * dx;
                out[0] = num / len - r;
            }

            // ── CircleCircleTangent ───────────────────────────────────────────
            Constraint::CircleCircleTangent { c1, c2, external } => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                let cx1 = params[e1.param(0)];
                let cy1 = params[e1.param(1)];
                let r1 = params[e1.param(2)];
                let cx2 = params[e2.param(0)];
                let cy2 = params[e2.param(1)];
                let r2 = params[e2.param(2)];
                let dx = cx2 - cx1;
                let dy = cy2 - cy1;
                let dist = (dx * dx + dy * dy).sqrt();
                if *external {
                    // External tangency: dist = r1 + r2
                    out[0] = dist - (r1 + r2);
                } else {
                    // Internal tangency: dist = |r1 - r2|
                    out[0] = dist - (r1 - r2).abs();
                }
            }

            // ── ArcArcTangent ─────────────────────────────────────────────────
            Constraint::ArcArcTangent { a1, a2, external } => {
                let e1 = &entities[*a1];
                let e2 = &entities[*a2];
                let cx1 = params[e1.param(0)];
                let cy1 = params[e1.param(1)];
                let r1 = params[e1.param(2)];
                let cx2 = params[e2.param(0)];
                let cy2 = params[e2.param(1)];
                let r2 = params[e2.param(2)];
                let dx = cx2 - cx1;
                let dy = cy2 - cy1;
                let dist = (dx * dx + dy * dy).sqrt();
                if *external {
                    out[0] = dist - (r1 + r2);
                } else {
                    out[0] = dist - (r1 - r2).abs();
                }
            }

            // ── Symmetric ─────────────────────────────────────────────────────
            Constraint::Symmetric { p1, p2, line } => {
                let (x1, y1) = p1.param_indices(entities);
                let (x2, y2) = p2.param_indices(entities);
                let el = &entities[*line];
                let lx1 = params[el.param(0)];
                let ly1 = params[el.param(1)];
                let lx2 = params[el.param(2)];
                let ly2 = params[el.param(3)];
                let ldx = lx2 - lx1;
                let ldy = ly2 - ly1;
                let len = (ldx * ldx + ldy * ldy).sqrt().max(SOLVER_NORM_TOL);
                // Midpoint of p1–p2
                let mx = (params[x1] + params[x2]) * 0.5;
                let my = (params[y1] + params[y2]) * 0.5;
                // Eq 1: midpoint lies on the line (cross product = 0)
                out[0] = ((lx2 - lx1) * (my - ly1) - (ly2 - ly1) * (mx - lx1)) / len;
                // Eq 2: segment p1–p2 is perpendicular to line direction (dot = 0)
                let sdx = params[x2] - params[x1];
                let sdy = params[y2] - params[y1];
                out[1] = (sdx * ldx + sdy * ldy) / len;
            }

            // ── Concentric ────────────────────────────────────────────────────
            Constraint::Concentric(c1, c2) => {
                let e1 = &entities[*c1];
                let e2 = &entities[*c2];
                out[0] = params[e1.param(0)] - params[e2.param(0)];
                out[1] = params[e1.param(1)] - params[e2.param(1)];
            }

            // ── Midpoint ──────────────────────────────────────────────────────
            Constraint::Midpoint { point, line } => {
                let (px, py) = point.param_indices(entities);
                let el = &entities[*line];
                let mx = (params[el.param(0)] + params[el.param(2)]) * 0.5;
                let my = (params[el.param(1)] + params[el.param(3)]) * 0.5;
                out[0] = params[px] - mx;
                out[1] = params[py] - my;
            }

            // ── Diameter ──────────────────────────────────────────────────────
            Constraint::Diameter { circle, diameter } => {
                let e = &entities[*circle];
                out[0] = params[e.param(2)] - diameter / 2.0;
            }
        }
        true
    }
}

// ── From<EntityId> for PointRef (Point entities only) ────────────────────────

impl From<EntityId> for PointRef {
    fn from(id: EntityId) -> Self {
        PointRef::Point(id)
    }
}
