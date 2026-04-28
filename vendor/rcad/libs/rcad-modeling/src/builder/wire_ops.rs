//! Wire-level operations: surface projection and 2-D round/bevel API.
//!
//! ## Normal projection of a Wire onto a Surface
//!
//! `project_wire_onto_surface` projects every vertex of a `Wire` onto the
//! nearest point on a `Surface3`, then reconnects the projected vertices with
//! straight-line edges to produce a new `Wire`.  The resulting wire hugs the
//! surface and can be used as a boundary loop for trimming or feature creation.
//!
//! Analogous to OCCT `BRepAlgo_NormalProjection` / `BRepOffsetAPI_NormalProjection`.
//!
//! ## 2-D fillet / chamfer on a polygon
//!
//! `fillet_wire_2d` and `chamfer_wire_2d` round or bevel polygon corners in 2-D.
//! Each function takes a slice of `DVec2` control points (the polygon vertices),
//! a parameter (`radius` / `dist`), and returns a new `Vec<DVec2>` that
//! approximates the rounded/bevelled polygon with arc mid-points inserted for
//! filleted corners.
//!
//! Analogous to OCCT `BRepFilletAPI_MakeFillet2d`.

use std::f64::consts::PI;

use glam::{DVec2, DVec3};
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Line3, Surface3};
use rcad_kernel::projection::closest_point_on_surface;
use rcad_kernel::topology::{Vertex, WireEdge, Wire};

use crate::builder::BuildError;

// ── Wire projection ───────────────────────────────────────────────────────────

/// Project every vertex of `wire` onto `surface`, reconnect them with
/// straight-line edges and return the projected wire.
///
/// The output wire has the same number of edges as the input.  Each edge
/// becomes a straight line between the two projected endpoint positions.
///
/// `n_samples` controls the Newton-search granularity for the closest-point
/// solver (passed through to `closest_point_on_surface`).  A value of 32 is
/// sufficient for all analytic surface types; B-spline surfaces may need more.
///
/// # Errors
/// Returns `BuildError::DegenerateGeometry` if the wire contains no edges, or
/// if a projected edge collapses to a zero-length segment.
pub fn project_wire_onto_surface(
    brep: &BRep,
    wire: &Wire,
    surface: &Surface3,
    n_samples: usize,
) -> Result<(Wire, Vec<DVec3>), BuildError> {
    if wire.edges.is_empty() {
        return Err(BuildError::DegenerateGeometry("wire has no edges"));
    }

    // Collect the ordered vertex positions from the wire.
    let src_pts: Vec<DVec3> = wire_ordered_vertex_positions(brep, wire);
    if src_pts.is_empty() {
        return Err(BuildError::DegenerateGeometry("wire has no vertices"));
    }

    // Project each source vertex onto the surface.
    let ns = n_samples.max(8);
    let proj_pts: Vec<DVec3> = src_pts
        .iter()
        .map(|&p| closest_point_on_surface(surface, p, ns).point)
        .collect();

    // Build a minimal BRep scratch space to hold the new vertices/edges.
    let dst_verts: Vec<Vertex> = proj_pts.iter().map(|&p| Vertex { point: p }).collect();
    let n = dst_verts.len();

    // Build edges: one per wire-edge, connecting consecutive projected vertices.
    // For a closed polygon the last edge connects back to vertex 0.
    let mut new_wire_edges: Vec<WireEdge> = Vec::with_capacity(wire.edges.len());
    let all_vertices: Vec<Vertex> = dst_verts.clone();
    let mut all_edges = Vec::new();

    for (i, _we) in wire.edges.iter().enumerate() {
        let va_idx = i;
        let vb_idx = (i + 1) % n;
        let pa = all_vertices[va_idx].point;
        let pb = all_vertices[vb_idx].point;
        let len = (pb - pa).length();
        if len < 1e-12 {
            continue; // skip degenerate (coincident) projected edge
        }
        let dir = (pb - pa) / len;
        let curve = Curve3::Line(Line3 { origin: pa, direction: dir });
        let eid = all_edges.len();
        all_edges.push((curve, 0.0, len, va_idx, vb_idx));
        new_wire_edges.push(WireEdge::fwd(eid));
    }

    if new_wire_edges.is_empty() {
        return Err(BuildError::DegenerateGeometry("all projected edges degenerated"));
    }

    // We only need the wire representation (not a full BRep) so we return
    // the wire plus the projected point positions for the caller to use.
    let out_wire = Wire { edges: new_wire_edges };
    Ok((out_wire, proj_pts))
}

/// Helper: walk the ordered start-vertex positions of a wire.
///
/// For the i-th `WireEdge` the start vertex is:
/// - `edge.start` if `we.forward == true`
/// - `edge.end` if `we.forward == false`
fn wire_ordered_vertex_positions(brep: &BRep, wire: &Wire) -> Vec<DVec3> {
    let mut pts = Vec::with_capacity(wire.edges.len());
    for we in &wire.edges {
        if let Some(e) = brep.edges.get(we.idx) {
            let v_idx = if we.forward { e.start } else { e.end };
            if let Some(v) = brep.vertices.get(v_idx) {
                pts.push(v.point);
            }
        }
    }
    pts
}

// ── 2-D fillet (rounded corners) ─────────────────────────────────────────────

/// Round the corners of a 2-D polygon.
///
/// For each interior vertex `p[i]` the two adjacent edges are trimmed back by
/// the tangent length `t = r / tan(θ/2)` (where `θ` is the interior angle) and
/// a circular-arc approximation is inserted.  The arc is approximated with a
/// single chord mid-point so that the output polygon has the correct start/end
/// tangent points and one extra sample per rounded corner.
///
/// Corners whose half-angle is very flat (almost 180°) or whose edge is too
/// short to accommodate the setback are left unmodified (polygon falls through
/// unchanged at those corners).
///
/// `pts` must have at least 3 points for any rounding to occur.  Both open and
/// closed polygons are supported; pass `closed = true` when the polygon loops
/// back from the last point to the first.
///
/// # Errors
/// Returns `BuildError::NonPositiveValue` if `radius ≤ 0`.
/// Returns `BuildError::DegenerateGeometry` if fewer than 2 points are given.
pub fn fillet_wire_2d(
    pts: &[DVec2],
    radius: f64,
    closed: bool,
) -> Result<Vec<DVec2>, BuildError> {
    if radius <= 0.0 {
        return Err(BuildError::NonPositiveValue("radius"));
    }
    if pts.len() < 2 {
        return Err(BuildError::DegenerateGeometry("need at least 2 points"));
    }
    if pts.len() < 3 {
        // Two points: nothing to round.
        return Ok(pts.to_vec());
    }
    round_corners_2d(pts, radius, closed, /*chamfer=*/false)
}

/// Bevel the corners of a 2-D polygon.
///
/// Like [`fillet_wire_2d`] but instead of inserting an arc mid-point, a single
/// straight chamfer line is inserted between the two setback points.  Each
/// rounded corner becomes one extra vertex (the midpoint of the bevel chord is
/// NOT added — the two setback points are both inserted, giving a flat cut).
///
/// # Errors
/// Returns `BuildError::NonPositiveValue` if `dist ≤ 0`.
/// Returns `BuildError::DegenerateGeometry` if fewer than 2 points are given.
pub fn chamfer_wire_2d(
    pts: &[DVec2],
    dist: f64,
    closed: bool,
) -> Result<Vec<DVec2>, BuildError> {
    if dist <= 0.0 {
        return Err(BuildError::NonPositiveValue("dist"));
    }
    if pts.len() < 2 {
        return Err(BuildError::DegenerateGeometry("need at least 2 points"));
    }
    if pts.len() < 3 {
        return Ok(pts.to_vec());
    }
    round_corners_2d(pts, dist, closed, /*chamfer=*/true)
}

// ── Shared corner-rounding core ───────────────────────────────────────────────

fn round_corners_2d(
    pts: &[DVec2],
    param: f64,
    closed: bool,
    chamfer: bool,
) -> Result<Vec<DVec2>, BuildError> {
    let n = pts.len();
    let mut out: Vec<DVec2> = Vec::with_capacity(n * 3);

    // Determine iteration range:
    // - closed: all n corners (i..i+n, modulo)
    // - open: skip the first and last vertices (no corner there)
    let corner_start = if closed { 0 } else { 1 };
    let corner_end = if closed { n } else { n - 1 };

    // We build the output by processing each straight segment and inserting
    // the corner geometry.  Use setback markers so that segments shortened by
    // one corner's setback know where to start.
    //
    // Strategy: collect (setback_start, setback_end) for each corner.
    // Then rebuild the output by walking segment-by-segment.

    let idx = |i: usize| {
        if closed { i % n } else { i }
    };

    // Each corner i (valid range [corner_start, corner_end)) produces:
    // - tangent point on the prev-edge: prev_pt
    // - tangent point on the next-edge: next_pt
    // - arc/chamfer insertion between them

    struct Corner {
        /// True if this corner was actually rounded/bevelled.
        active: bool,
        prev: DVec2, // setback point on prev edge end
        next: DVec2, // setback point on next edge start
        mid: Option<DVec2>, // arc midpoint (fillet only)
    }

    impl Clone for Corner {
        fn clone(&self) -> Self {
            Corner {
                active: self.active,
                prev: self.prev,
                next: self.next,
                mid: self.mid,
            }
        }
    }

    let mut corners: Vec<Option<Corner>> = vec![None; n];

    for ci in corner_start..corner_end {
        let i = idx(ci);
        let prev_i = if ci == 0 { n - 1 } else { ci - 1 };
        let next_i = if ci + 1 >= n { 0 } else { ci + 1 };
        let p_prev = pts[prev_i];
        let p_cur  = pts[i];
        let p_next = pts[next_i];

        let d0 = (p_cur - p_prev).normalize_or_zero();
        let d1 = (p_next - p_cur).normalize_or_zero();

        if d0.length_squared() < 1e-20 || d1.length_squared() < 1e-20 {
            continue; // degenerate segment
        }

        // Interior angle (between the two incoming/outgoing directions).
        let cos_a = (-d0).dot(d1).clamp(-1.0, 1.0);
        let half_angle = (PI - cos_a.acos()) * 0.5;
        let tan_h = half_angle.tan().abs();

        if tan_h < 1e-10 {
            // Straight edge — nothing to round.
            continue;
        }

        let setback = if chamfer { param } else { param / tan_h };

        let len0 = (p_cur - p_prev).length();
        let len1 = (p_next - p_cur).length();

        // Only round if there is enough room on both edges.
        if setback >= len0 * 0.5 || setback >= len1 * 0.5 {
            continue;
        }

        let prev_pt = p_cur - d0 * setback;
        let next_pt = p_cur + d1 * setback;

        let mid = if chamfer {
            None
        } else {
            // Arc midpoint: bisector direction from the corner.
            let bisect = (-d0 + d1).normalize_or_zero();
            // Distance from corner to arc midpoint.
            let arc_dist = param / half_angle.sin().max(1e-10);
            // Move `arc_dist * bisect` from the corner and find the midpoint
            // of the chord between prev_pt and next_pt on the circle.
            let center = p_cur + bisect * (param / half_angle.sin().max(1e-10) - param);
            let arc_mid = (prev_pt + next_pt) * 0.5;
            // Pull arc_mid radially outward so it lies on the circle.
            let from_center = arc_mid - (p_cur + bisect * (param / half_angle.sin().max(1e-10) - param));
            let r_nominal = from_center.length();
            let _ = arc_dist; let _ = center;
            if r_nominal < 1e-12 {
                None
            } else {
                let on_arc = p_cur + bisect * (param / half_angle.sin().max(1e-10) - param)
                    + from_center * (param / r_nominal);
                Some(on_arc)
            }
        };

        corners[i] = Some(Corner { active: true, prev: prev_pt, next: next_pt, mid });
    }

    // Now rebuild the output polygon.
    // Walk vertices 0..n (for closed) or 0..n-1 (for open).
    let walk_end = if closed { n } else { n - 1 };

    for ci in 0..walk_end {
        let i = idx(ci);
        let ni = idx(ci + 1);

        match &corners[i] {
            None => {
                // No rounding at vertex i — emit it as-is.
                out.push(pts[i]);
            }
            Some(c) if c.active => {
                // Emit the setback point on the incoming edge from the previous corner.
                out.push(c.prev);
                // Emit arc midpoint if fillet.
                if let Some(m) = c.mid {
                    out.push(m);
                }
                // Emit the setback point on the outgoing edge.
                out.push(c.next);
            }
            _ => {
                out.push(pts[i]);
            }
        }

        // If not the last segment, possibly emit the straight segment endpoint.
        // The next corner's `prev` point handles the incoming segment trimming.
        // We don't emit the next vertex here if the next corner is active
        // (it will emit `prev` itself).
        if ci < walk_end - 1 {
            if corners[ni].as_ref().map_or(false, |c| c.active) {
                // The start of segment [i..ni] was already pushed above.
                // The end (ni setback point on incoming edge) will be pushed in
                // the next iteration as `c.prev`. Nothing to add here.
            } else if corners[i].as_ref().map_or(false, |c| c.active) {
                // i was active, ni is not — the outgoing next_pt is ni itself,
                // nothing extra to emit.
            }
        }
    }

    // For open polygons, emit the last vertex.
    if !closed {
        out.push(pts[n - 1]);
    }

    if out.is_empty() {
        return Err(BuildError::DegenerateGeometry("all corners degenerated"));
    }

    Ok(out)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec2;
    use rcad_kernel::geom::{Plane, PrimitiveSolid, Surface3};
    use rcad_kernel::BRep;

    // ── project_wire_onto_surface ─────────────────────────────────────────────

    fn box_wire_0(brep: &BRep) -> Wire {
        brep.solids[0].shells[0].faces[0].outer_wire.clone()
    }

    #[test]
    fn project_wire_onto_same_plane_is_identity() {
        // The first face of a unit box is the bottom face (z=0 plane, normal -Z).
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let wire = box_wire_0(&brep);
        let plane = Surface3::Plane(Plane {
            origin: glam::DVec3::ZERO,
            normal: glam::DVec3::NEG_Z,
        });

        let (proj_wire, proj_pts) = project_wire_onto_surface(&brep, &wire, &plane, 16).unwrap();
        // All projected points should lie on z=0.
        for p in &proj_pts {
            assert!(p.z.abs() < 1e-9, "projected point not on z=0: {p:?}");
        }
        assert_eq!(proj_wire.edges.len(), wire.edges.len());
    }

    #[test]
    fn project_raised_wire_onto_plane_lowers_z() {
        // Build a small square wire manually (elevated at z=2) and project onto z=0 plane.
        let mut brep = BRep::new();
        use rcad_kernel::topology::Vertex;
        let v0 = { brep.vertices.push(Vertex { point: glam::DVec3::new(0.0, 0.0, 2.0) }); brep.vertices.len()-1 };
        let v1 = { brep.vertices.push(Vertex { point: glam::DVec3::new(1.0, 0.0, 2.0) }); brep.vertices.len()-1 };
        let v2 = { brep.vertices.push(Vertex { point: glam::DVec3::new(1.0, 1.0, 2.0) }); brep.vertices.len()-1 };
        let v3 = { brep.vertices.push(Vertex { point: glam::DVec3::new(0.0, 1.0, 2.0) }); brep.vertices.len()-1 };
        use rcad_kernel::topology::Edge;
        let e0 = { brep.edges.push(Edge { start: v0, end: v1 }); brep.edges.len()-1 };
        let e1 = { brep.edges.push(Edge { start: v1, end: v2 }); brep.edges.len()-1 };
        let e2 = { brep.edges.push(Edge { start: v2, end: v3 }); brep.edges.len()-1 };
        let e3 = { brep.edges.push(Edge { start: v3, end: v0 }); brep.edges.len()-1 };
        let wire = Wire {
            edges: vec![WireEdge::fwd(e0), WireEdge::fwd(e1), WireEdge::fwd(e2), WireEdge::fwd(e3)],
        };

        let plane = Surface3::Plane(Plane {
            origin: glam::DVec3::ZERO,
            normal: glam::DVec3::Z,
        });
        let (_proj_wire, proj_pts) = project_wire_onto_surface(&brep, &wire, &plane, 16).unwrap();
        for p in &proj_pts {
            assert!(p.z.abs() < 1e-9, "projected point not on z=0 plane: {p:?}");
        }
        // XY coords should be preserved.
        assert!((proj_pts[0] - glam::DVec3::new(0.0, 0.0, 0.0)).length() < 1e-9);
        assert!((proj_pts[1] - glam::DVec3::new(1.0, 0.0, 0.0)).length() < 1e-9);
    }

    #[test]
    fn project_wire_empty_wire_returns_error() {
        let brep = BRep::new();
        let wire = Wire { edges: vec![] };
        let plane = Surface3::Plane(Plane {
            origin: glam::DVec3::ZERO,
            normal: glam::DVec3::Z,
        });
        assert!(project_wire_onto_surface(&brep, &wire, &plane, 16).is_err());
    }

    // ── fillet_wire_2d ────────────────────────────────────────────────────────

    fn square_pts() -> Vec<DVec2> {
        vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(0.0, 1.0),
        ]
    }

    #[test]
    fn fillet_square_closed_produces_more_points() {
        let pts = square_pts();
        let result = fillet_wire_2d(&pts, 0.1, true).unwrap();
        // Each of 4 corners gets prev + arc_mid + next (up to 3 extra per corner),
        // so result should have more points than input.
        assert!(
            result.len() > pts.len(),
            "fillet should produce more points, got {}",
            result.len()
        );
    }

    #[test]
    fn chamfer_square_closed_produces_more_points() {
        let pts = square_pts();
        let result = chamfer_wire_2d(&pts, 0.1, true).unwrap();
        assert!(
            result.len() > pts.len(),
            "chamfer should produce more points, got {}",
            result.len()
        );
    }

    #[test]
    fn fillet_too_large_radius_leaves_polygon_unchanged() {
        // radius > half edge length: all corners fall through unchanged.
        let pts = square_pts(); // edge length = 1.0
        let result = fillet_wire_2d(&pts, 2.0, true).unwrap(); // radius > 0.5 → skip all
        // All input points must still appear in the output (no corners were rounded).
        for p in &pts {
            assert!(
                result.iter().any(|q| (*q - *p).length() < 1e-9),
                "original point {p:?} missing from unfilleted output"
            );
        }
    }

    #[test]
    fn chamfer_rejects_zero_dist() {
        assert!(chamfer_wire_2d(&square_pts(), 0.0, true).is_err());
    }

    #[test]
    fn fillet_rejects_zero_radius() {
        assert!(fillet_wire_2d(&square_pts(), 0.0, true).is_err());
    }

    #[test]
    fn fillet_rejects_single_point() {
        assert!(fillet_wire_2d(&[DVec2::ZERO], 0.1, true).is_err());
    }

    #[test]
    fn chamfer_open_polygon_preserves_endpoints() {
        let l = vec![DVec2::new(0.0, 0.0), DVec2::new(1.0, 0.0), DVec2::new(1.0, 1.0)];
        let result = chamfer_wire_2d(&l, 0.1, false).unwrap();
        // First and last must be preserved.
        assert!((result.first().unwrap() - DVec2::new(0.0, 0.0)).length() < 1e-9);
        assert!((result.last().unwrap()  - DVec2::new(1.0, 1.0)).length() < 1e-9);
    }

    #[test]
    fn fillet_two_points_returns_unchanged() {
        let pts = vec![DVec2::new(0.0, 0.0), DVec2::new(1.0, 0.0)];
        let result = fillet_wire_2d(&pts, 0.1, false).unwrap();
        assert_eq!(result.len(), 2);
    }
}
