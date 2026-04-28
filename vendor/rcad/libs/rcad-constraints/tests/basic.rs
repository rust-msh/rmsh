//! Integration tests for the GCS solver.

use rcad_constraints::{Constraint, Sketch};
use rcad_constraints::constraint::Constraint::*;
use rcad_constraints::entity::PointRef;

const TOL: f64 = 1e-7;

// ─── helpers ──────────────────────────────────────────────────────────────────

fn assert_near(a: f64, b: f64, label: &str) {
    assert!((a - b).abs() < TOL, "{label}: expected {b}, got {a}");
}

// ─── point constraints ────────────────────────────────────────────────────────

/// Fix a point at the origin.
#[test]
fn fixed_point() {
    let mut sk = Sketch::new();
    let p = sk.add_point(3.0, 7.0);
    sk.add_constraint(Constraint::fix_point(p, 0.0, 0.0));
    let res = sk.solve();
    assert!(res.converged, "not converged: residual={}", res.residual);
    let c = sk.point_coords(PointRef::Point(p));
    assert_near(c.x, 0.0, "x");
    assert_near(c.y, 0.0, "y");
}

/// Two points must be coincident.
#[test]
fn coincident_points() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(5.0, 3.0);
    sk.add_constraint(Constraint::fix_point(p1, 1.0, 2.0));
    sk.add_constraint(Constraint::coincident(p1, p2));
    let res = sk.solve();
    assert!(res.converged);
    let c1 = sk.point_coords(PointRef::Point(p1));
    let c2 = sk.point_coords(PointRef::Point(p2));
    assert_near(c1.x, c2.x, "coincident x");
    assert_near(c1.y, c2.y, "coincident y");
}

/// Distance constraint: p1 fixed at origin, p2 at distance 5.
#[test]
fn point_distance() {
    let mut sk = Sketch::new();
    let p1 = sk.add_point(0.0, 0.0);
    let p2 = sk.add_point(3.0, 4.0); // initial distance 5 — already satisfied
    sk.add_constraint(Constraint::fix_point(p1, 0.0, 0.0));
    // fix y of p2 so the system is fully constrained
    let p2_y_param = sk.entities[p2].param_start + 1;
    sk.fix_param(p2_y_param);
    sk.add_constraint(Constraint::point_distance(p1, p2, 5.0));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c2 = sk.point_coords(PointRef::Point(p2));
    let dist = (c2.x * c2.x + c2.y * c2.y).sqrt();
    assert_near(dist, 5.0, "distance");
}

// ─── line constraints ─────────────────────────────────────────────────────────

/// Horizontal constraint: line becomes horizontal.
#[test]
fn horizontal_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 2.0, 1.5); // not horizontal
    sk.add_constraint(Horizontal(l));
    // Fix start point to remove translation DOF
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(l);
    assert_near(p[1], p[3], "y1 == y2 (horizontal)");
}

/// Vertical constraint.
#[test]
fn vertical_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 1.5, 2.0);
    sk.add_constraint(Vertical(l));
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(l);
    assert_near(p[0], p[2], "x1 == x2 (vertical)");
}

/// Line length constraint.
#[test]
fn line_length() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 3.0, 0.0); // length 3, want 5
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    // fix direction: keep y2=0 (horizontal)
    sk.add_constraint(Horizontal(l));
    sk.add_constraint(LineLength { line: l, length: 5.0 });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(l);
    let len = ((p[2] - p[0]).powi(2) + (p[3] - p[1]).powi(2)).sqrt();
    assert_near(len, 5.0, "length");
}

/// Parallel constraint: two lines become parallel.
#[test]
fn parallel_lines() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 2.0, 0.0); // horizontal
    let l2 = sk.add_line(0.0, 1.0, 2.0, 1.5); // slightly tilted
    sk.fix_entity(l1); // fix l1 as reference
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l2), 0.0, 1.0));
    sk.add_constraint(Parallel(l1, l2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(l1);
    let p2 = sk.entity_params(l2);
    let dy1 = p1[3] - p1[1];
    let dx1 = p1[2] - p1[0];
    let dy2 = p2[3] - p2[1];
    let dx2 = p2[2] - p2[0];
    // cross product ≈ 0
    assert_near(dx1 * dy2 - dy1 * dx2, 0.0, "parallel cross product");
}

/// Perpendicular constraint.
#[test]
fn perpendicular_lines() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 2.0, 0.0); // horizontal
    let l2 = sk.add_line(1.0, 0.0, 1.5, 1.5); // almost vertical
    sk.fix_entity(l1);
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l2), 1.0, 0.0));
    sk.add_constraint(Perpendicular(l1, l2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(l1);
    let p2 = sk.entity_params(l2);
    let dot = (p1[2] - p1[0]) * (p2[2] - p2[0])
            + (p1[3] - p1[1]) * (p2[3] - p2[1]);
    assert_near(dot, 0.0, "perpendicular dot product");
}

/// Equal length constraint.
#[test]
fn equal_length() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 4.0, 0.0); // length 4
    let l2 = sk.add_line(0.0, 1.0, 2.0, 1.0); // length 2 → should become 4
    sk.fix_entity(l1);
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l2), 0.0, 1.0));
    sk.add_constraint(Horizontal(l2));
    sk.add_constraint(EqualLength(l1, l2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p2 = sk.entity_params(l2);
    let len2 = ((p2[2] - p2[0]).powi(2) + (p2[3] - p2[1]).powi(2)).sqrt();
    assert_near(len2, 4.0, "equal length");
}

// ─── circle constraints ───────────────────────────────────────────────────────

/// Radius constraint.
#[test]
fn circle_radius() {
    let mut sk = Sketch::new();
    let c = sk.add_circle(0.0, 0.0, 1.0);
    sk.add_constraint(Constraint::fix_point(PointRef::Center(c), 0.0, 0.0));
    sk.add_constraint(Radius { circle: c, radius: 3.0 });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(c);
    assert_near(p[2], 3.0, "radius");
}

/// Point on circle.
#[test]
fn point_on_circle() {
    let mut sk = Sketch::new();
    let c = sk.add_circle(0.0, 0.0, 2.0);
    sk.fix_entity(c);
    let p = sk.add_point(3.0, 0.0); // outside circle
    sk.add_constraint(Constraint::fix_point(p, 0.0, 3.0)); // fix to (0,3)
    // The point at (0,3) is at distance 3 from origin, but circle radius is 2 — can't both hold.
    // Instead: just fix y and let PointOnCircle pull x to satisfy radius.
    // Re-do: fix only y of p at 0, constrain PointOnCircle.
    let mut sk2 = Sketch::new();
    let c2 = sk2.add_circle(0.0, 0.0, 2.0);
    sk2.fix_entity(c2);
    let p2 = sk2.add_point(1.5, 1.0);
    // fix y at 0 so x can adjust to land on circle
    let p2_y = sk2.entities[p2].param_start + 1;
    sk2.fix_param(p2_y);
    sk2.add_constraint(PointOnCircle { point: PointRef::Point(p2), circle: c2 });
    let res2 = sk2.solve();
    assert!(res2.converged, "residual={}", res2.residual);
    let coords = sk2.point_coords(PointRef::Point(p2));
    let r = (coords.x * coords.x + coords.y * coords.y).sqrt();
    assert_near(r, 2.0, "point on circle radius");
}

// ─── DOF analysis ─────────────────────────────────────────────────────────────

#[test]
fn dof_fully_constrained_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 1.0, 0.0);
    // 4 params; add: fix_start(2) + horizontal(1) + line_length(1) = 4 equations
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l), 0.0, 0.0));
    sk.add_constraint(Horizontal(l));
    sk.add_constraint(LineLength { line: l, length: 5.0 });
    assert_eq!(sk.dof(), 0, "fully constrained");
}

#[test]
fn dof_under_constrained() {
    let mut sk = Sketch::new();
    let _p = sk.add_point(1.0, 2.0); // 2 params, 0 constraints → DOF=2
    assert_eq!(sk.dof(), 2);
}

// ─── wire BRep ───────────────────────────────────────────────────────────────

#[test]
fn to_wire_brep_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 3.0, 4.0);
    sk.fix_entity(l);
    sk.solve();
    let brep = sk.to_wire_brep();
    // 1 line → 2 vertices, 1 edge
    assert_eq!(brep.vertices.len(), 2);
    assert_eq!(brep.edges.len(), 1);
    assert_eq!(brep.geom.curves.len(), 1);
}

#[test]
fn to_wire_brep_circle() {
    let mut sk = Sketch::new();
    let _c = sk.add_circle(1.0, 2.0, 3.0);
    sk.solve();
    let brep = sk.to_wire_brep();
    assert_eq!(brep.vertices.len(), 1);  // one vertex (start==end)
    assert_eq!(brep.edges.len(), 1);
    assert_eq!(brep.geom.curves.len(), 1);
}

// ─── additional constraint tests ─────────────────────────────────────────────

/// Angle constraint: two lines meet at a 90° angle.
#[test]
fn angle_between_lines() {
    use std::f64::consts::FRAC_PI_2;
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 2.0, 0.0);
    let l2 = sk.add_line(0.0, 0.0, 1.5, 1.5);
    sk.fix_entity(l1);
    sk.add_constraint(Constraint::fix_point(PointRef::LineStart(l2), 0.0, 0.0));
    sk.add_constraint(Angle { l1, l2, angle_rad: FRAC_PI_2 });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(l1);
    let p2 = sk.entity_params(l2);
    let dot = (p1[2] - p1[0]) * (p2[2] - p2[0]) + (p1[3] - p1[1]) * (p2[3] - p2[1]);
    assert_near(dot, 0.0, "90° angle => dot product = 0");
}

/// EqualRadius: two circles end up with the same radius.
#[test]
fn equal_radius_circles() {
    let mut sk = Sketch::new();
    let c1 = sk.add_circle(0.0, 0.0, 3.0);
    let c2 = sk.add_circle(5.0, 0.0, 1.0);
    sk.fix_entity(c1);
    sk.add_constraint(Constraint::fix_point(PointRef::Center(c2), 5.0, 0.0));
    sk.add_constraint(EqualRadius(c1, c2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(c1);
    let p2 = sk.entity_params(c2);
    assert_near(p1[2], p2[2], "equal radius");
}

/// PointOnLine: a free point is pulled onto a fixed line.
#[test]
fn point_on_line_constraint() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 4.0, 0.0);
    sk.fix_entity(l);
    let p = sk.add_point(2.0, 3.0);
    let px_idx = sk.entities[p].param_start;
    sk.fix_param(px_idx);
    sk.add_constraint(PointOnLine { point: PointRef::Point(p), line: l });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let coords = sk.point_coords(PointRef::Point(p));
    assert_near(coords.y, 0.0, "point y should be 0 (on horizontal line)");
}

/// Tangent (circle-line): circle with fixed radius moves so it is tangent to
/// a horizontal line at y=0.  The center y-coordinate should converge to r.
#[test]
fn circle_tangent_to_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(-5.0, 0.0, 5.0, 0.0);
    sk.fix_entity(l);
    // Circle at (0, 2), radius 1.  Fix x and radius; let y move to satisfy tangency.
    let c = sk.add_circle(0.0, 2.0, 1.0);
    let c_params = sk.entities[c].param_start;
    sk.fix_param(c_params);     // fix cx
    sk.fix_param(c_params + 2); // fix r
    sk.add_constraint(Tangent { circle: c, line: l });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(c);
    // dist(center, line y=0) = |cy| = r = 1
    assert_near(p[1].abs(), p[2], "|center y| should equal radius for tangency");
}

/// CircleCircleTangent (external): c2 slides along X until externally tangent
/// to c1 (r1=2, r2=1 → expected dist = 3).
#[test]
fn circle_circle_external_tangent() {
    use rcad_constraints::constraint::Constraint::CircleCircleTangent;
    let mut sk = Sketch::new();
    let c1 = sk.add_circle(0.0, 0.0, 2.0);
    sk.fix_entity(c1);
    // c2 starts at (6, 0) with r=1.  Fix y and r; let x move.
    let c2 = sk.add_circle(6.0, 0.0, 1.0);
    let c2_params = sk.entities[c2].param_start;
    sk.fix_param(c2_params + 1); // fix cy
    sk.fix_param(c2_params + 2); // fix r
    sk.add_constraint(CircleCircleTangent { c1, c2, external: true });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(c1);
    let p2 = sk.entity_params(c2);
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let dist = (dx * dx + dy * dy).sqrt();
    assert_near(dist, p1[2] + p2[2], "external tangency: dist = r1 + r2");
}

/// ArcArcTangent (external): two arcs become externally tangent.
#[test]
fn arc_arc_external_tangent() {
    use rcad_constraints::constraint::Constraint::ArcArcTangent;
    let mut sk = Sketch::new();
    let a1 = sk.add_arc(0.0, 0.0, 2.0, 0.0, std::f64::consts::PI);
    sk.fix_entity(a1);
    // a2 starts at (8, 0) with r=1. Fix cy and r; let cx move.
    let a2 = sk.add_arc(8.0, 0.0, 1.0, 0.0, std::f64::consts::PI);
    let a2_params = sk.entities[a2].param_start;
    sk.fix_param(a2_params + 1); // fix cy
    sk.fix_param(a2_params + 2); // fix r
    sk.add_constraint(ArcArcTangent { a1, a2, external: true });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(a1);
    let p2 = sk.entity_params(a2);
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let dist = (dx * dx + dy * dy).sqrt();
    assert_near(dist, p1[2] + p2[2], "arc external tangency: dist = r1 + r2");
}

/// Symmetric: two points are mirror images across a horizontal line.
#[test]
fn symmetric_across_line() {
    use rcad_constraints::constraint::Constraint::Symmetric;
    let mut sk = Sketch::new();
    // Mirror line: y=2 (horizontal)
    let mirror = sk.add_line(-5.0, 2.0, 5.0, 2.0);
    sk.fix_entity(mirror);
    // p1 fixed at (3, 0); p2 free — should end up at (3, 4)
    let p1 = sk.add_point(3.0, 0.0);
    sk.add_constraint(Constraint::fix_point(p1, 3.0, 0.0));
    let p2 = sk.add_point(3.0, 1.0); // wrong y, should move to 4
    // Fix x of p2 so only y is free
    let p2_x = sk.entities[p2].param_start;
    sk.fix_param(p2_x);
    sk.add_constraint(Symmetric {
        p1: PointRef::Point(p1),
        p2: PointRef::Point(p2),
        line: mirror,
    });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let c2 = sk.point_coords(PointRef::Point(p2));
    assert_near(c2.x, 3.0, "symmetric x");
    assert_near(c2.y, 4.0, "symmetric y (mirror of y=0 across y=2)");
}

/// to_solid_brep: a square sketch extruded to height 2 should produce 6 faces.
#[test]
fn sketch_to_solid_brep_square() {
    let mut sk = Sketch::new();
    // Unit square: 4 lines forming a closed polygon
    sk.add_line(0.0, 0.0, 1.0, 0.0);
    sk.add_line(1.0, 0.0, 1.0, 1.0);
    sk.add_line(1.0, 1.0, 0.0, 1.0);
    sk.add_line(0.0, 1.0, 0.0, 0.0);
    sk.solve();

    let solid = sk.to_solid_brep(2.0).expect("to_solid_brep should succeed for a closed square");
    let n_faces = solid.solids[0].shells[0].faces.len();
    assert_eq!(n_faces, 6, "extruded square should have 6 faces, got {n_faces}");
}

// ─── new constraint tests ─────────────────────────────────────────────────────

/// Concentric: two circles share the same center.
#[test]
fn concentric_circles() {
    let mut sk = Sketch::new();
    let c1 = sk.add_circle(0.0, 0.0, 3.0);
    sk.fix_entity(c1);
    let c2 = sk.add_circle(5.0, 4.0, 1.0); // wrong center
    sk.add_constraint(Constraint::Concentric(c1, c2));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p1 = sk.entity_params(c1);
    let p2 = sk.entity_params(c2);
    assert_near(p1[0], p2[0], "concentric cx");
    assert_near(p1[1], p2[1], "concentric cy");
}

/// Midpoint: a point is pulled to the midpoint of a line.
#[test]
fn midpoint_on_line() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 4.0, 2.0);
    sk.fix_entity(l);
    let p = sk.add_point(0.0, 0.0); // wrong position
    sk.add_constraint(Constraint::Midpoint { point: PointRef::Point(p), line: l });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let coords = sk.point_coords(PointRef::Point(p));
    assert_near(coords.x, 2.0, "midpoint x");
    assert_near(coords.y, 1.0, "midpoint y");
}

/// Diameter: circle diameter constrained to a value.
#[test]
fn diameter_constraint() {
    let mut sk = Sketch::new();
    let c = sk.add_circle(0.0, 0.0, 1.0); // radius 1, want diameter 6 → radius 3
    sk.add_constraint(Constraint::fix_point(PointRef::Center(c), 0.0, 0.0));
    sk.add_constraint(Constraint::Diameter { circle: c, diameter: 6.0 });
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
    let p = sk.entity_params(c);
    assert_near(p[2], 3.0, "radius should be diameter/2");
}

/// Rectangle with diagonal: 4 lines forming a rectangle + diagonal line.
#[test]
fn rectangle_with_diagonal() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 3.0, 0.0); // bottom
    let l2 = sk.add_line(3.0, 0.0, 3.0, 2.0); // right
    let l3 = sk.add_line(3.0, 2.0, 0.0, 2.0); // top
    let l4 = sk.add_line(0.0, 2.0, 0.0, 0.0); // left
    let d = sk.add_line(0.0, 0.0, 3.0, 2.0);  // diagonal
    // Fix the rectangle
    sk.fix_entity(l1);
    // Coincident constraints for corners
    sk.add_constraint(Constraint::coincident(PointRef::LineEnd(l1), PointRef::LineStart(l2)));
    sk.add_constraint(Constraint::coincident(PointRef::LineEnd(l2), PointRef::LineStart(l3)));
    sk.add_constraint(Constraint::coincident(PointRef::LineEnd(l3), PointRef::LineStart(l4)));
    sk.add_constraint(Constraint::coincident(PointRef::LineEnd(l4), PointRef::LineStart(l1)));
    // Diagonal connects l1 start to l3 start
    sk.add_constraint(Constraint::coincident(PointRef::LineStart(d), PointRef::LineStart(l1)));
    sk.add_constraint(Constraint::coincident(PointRef::LineEnd(d), PointRef::LineStart(l3)));
    // Opposite sides equal
    sk.add_constraint(EqualLength(l1, l3));
    sk.add_constraint(EqualLength(l2, l4));
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
}

/// Four-bar linkage: 4 connected lines, fix one bar, verify DOF.
#[test]
fn four_bar_linkage() {
    let mut sk = Sketch::new();
    let l1 = sk.add_line(0.0, 0.0, 4.0, 0.0); // ground (fixed)
    let l2 = sk.add_line(0.0, 0.0, 1.5, 2.5); // left arm
    let l3 = sk.add_line(1.5, 2.5, 5.5, 2.5); // coupler
    let l4 = sk.add_line(5.5, 2.5, 4.0, 0.0); // right arm
    sk.fix_entity(l1);
    // Pin joints
    sk.add_constraint(Constraint::coincident(PointRef::LineStart(l2), PointRef::LineStart(l1)));
    sk.add_constraint(Constraint::coincident(PointRef::LineEnd(l2), PointRef::LineStart(l3)));
    sk.add_constraint(Constraint::coincident(PointRef::LineEnd(l3), PointRef::LineStart(l4)));
    sk.add_constraint(Constraint::coincident(PointRef::LineEnd(l4), PointRef::LineEnd(l1)));
    // Fix lengths
    sk.add_constraint(LineLength { line: l2, length: 3.0 });
    sk.add_constraint(LineLength { line: l3, length: 4.0 });
    sk.add_constraint(LineLength { line: l4, length: 3.0 });
    // Four-bar has 1 DOF (input angle free)
    assert_eq!(sk.dof(), 1, "four-bar should have 1 DOF");
    let res = sk.solve();
    assert!(res.converged, "residual={}", res.residual);
}
