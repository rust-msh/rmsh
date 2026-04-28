//! First-stage feature operations (TKFeat-like APIs).
//!
//! This module builds practical feature workflows on top of the existing
//! boolean kernel. The first shipped feature is a cylindrical hole.

use glam::{DAffine3, DMat3, DVec3};
use rcad_kernel::{BRep, GeomStore, PrimitiveSolid};
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

use crate::{BooleanError, BooleanOpType, boolean_op};

/// Errors returned by feature operations.
#[derive(Debug)]
pub enum FeatureError {
    NonFiniteInput(&'static str),
    NonPositiveInput(&'static str),
    InvalidInput(&'static str),
    ZeroVector(&'static str),
    ParallelVectors(&'static str, &'static str),
    Modeling(String),
    Boolean(BooleanError),
}

impl std::fmt::Display for FeatureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFiniteInput(name) => write!(f, "{name} must be finite"),
            Self::NonPositiveInput(name) => write!(f, "{name} must be > 0"),
            Self::InvalidInput(name) => write!(f, "{name} is invalid"),
            Self::ZeroVector(name) => write!(f, "{name} must be non-zero"),
            Self::ParallelVectors(a, b) => write!(f, "{a} must not be parallel to {b}"),
            Self::Modeling(msg) => write!(f, "modeling operation failed: {msg}"),
            Self::Boolean(err) => write!(f, "boolean operation failed: {err}"),
        }
    }
}

impl std::error::Error for FeatureError {}

impl From<BooleanError> for FeatureError {
    fn from(value: BooleanError) -> Self {
        Self::Boolean(value)
    }
}

impl From<rcad_modeling::BuildError> for FeatureError {
    fn from(value: rcad_modeling::BuildError) -> Self {
        Self::Modeling(value.to_string())
    }
}

const EPS: f64 = 1e-12;

fn validate_finite(name: &'static str, v: f64) -> Result<f64, FeatureError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(FeatureError::NonFiniteInput(name))
    }
}

fn validate_positive(name: &'static str, v: f64) -> Result<f64, FeatureError> {
    let v = validate_finite(name, v)?;
    if v > 0.0 {
        Ok(v)
    } else {
        Err(FeatureError::NonPositiveInput(name))
    }
}

fn normalize(name: &'static str, v: DVec3) -> Result<DVec3, FeatureError> {
    if !v.is_finite() {
        return Err(FeatureError::NonFiniteInput(name));
    }
    if v.length_squared() <= EPS {
        return Err(FeatureError::ZeroVector(name));
    }
    Ok(v.normalize())
}

fn axis_ref_basis(axis: DVec3, ref_dir: DVec3) -> Result<(DVec3, DVec3, DVec3), FeatureError> {
    let y_axis = normalize("axis", axis)?;
    let ref_dir = normalize("ref_dir", ref_dir)?;
    let x_reject = ref_dir - y_axis * ref_dir.dot(y_axis);
    if x_reject.length_squared() <= EPS {
        return Err(FeatureError::ParallelVectors("ref_dir", "axis"));
    }
    let x_axis = x_reject.normalize();
    let z_axis = x_axis.cross(y_axis).normalize();
    Ok((x_axis, y_axis, z_axis))
}

/// Create a cylindrical through/blind hole by subtracting an oriented cylinder
/// from `target`.
///
/// - `center`: center of the tool cylinder.
/// - `axis`: cylinder axis direction.
/// - `ref_dir`: reference direction used to build local orientation.
/// - `radius`: hole radius.
/// - `depth`: tool cylinder height.
///
/// For through holes, pass a `depth` larger than the part thickness along
/// `axis`.
pub fn make_cylindrical_hole(
    target: &BRep,
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    depth: f64,
) -> Result<BRep, FeatureError> {
    if !center.is_finite() {
        return Err(FeatureError::NonFiniteInput("center"));
    }
    let radius = validate_positive("radius", radius)?;
    let depth = validate_positive("depth", depth)?;

    let (x_axis, y_axis, z_axis) = axis_ref_basis(axis, ref_dir)?;

    let mut tool = BRep::from_primitive(PrimitiveSolid::Cylinder { radius, height: depth });
    let rot = DMat3::from_cols(x_axis, y_axis, z_axis);
    tool.apply_transform(DAffine3::from_mat3_translation(rot, center));

    Ok(boolean_op(BooleanOpType::Difference, target, &tool)?)
}

/// Create a prismatic boss or pocket by extruding a polygon profile and
/// performing a boolean union (boss) or difference (pocket) with `target`.
///
/// - `profile_verts`: 3D coplanar polygon vertices in CCW order when viewed
///   along the extrusion direction.  Minimum 3 vertices required.
/// - `direction`: extrusion direction (unit vector is computed internally).
/// - `depth`: extrusion length (must be > 0).
/// - `op`: [`BooleanOpType::Union`] = boss; [`BooleanOpType::Difference`] = pocket.
///
/// Analogous to OCCT `BRepFeat_MakePrism` for linear boss/pocket features.
pub fn make_prism(
    target: &BRep,
    profile_verts: &[DVec3],
    direction: DVec3,
    depth: f64,
    op: BooleanOpType,
) -> Result<BRep, FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::NonPositiveInput("profile_verts needs >= 3 vertices"));
    }
    let dir = normalize("direction", direction)?;
    let depth = validate_positive("depth", depth)?;

    let tool = build_polygon_prism(profile_verts, dir, depth)?;
    Ok(boolean_op(op, target, &tool)?)
}

/// Create a drafted prismatic boss or pocket by extruding a polygon profile
/// with radial taper and applying a boolean operation.
///
/// Positive `draft_angle_rad` expands the top profile outward; negative values
/// shrink it inward.
///
/// Analogous to OCCT `BRepFeat_MakeDPrism` (linear draft prism).
pub fn make_draft_prism(
    target: &BRep,
    profile_verts: &[DVec3],
    direction: DVec3,
    depth: f64,
    draft_angle_rad: f64,
    op: BooleanOpType,
) -> Result<BRep, FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::InvalidInput("profile_verts needs >= 3 vertices"));
    }
    let dir = normalize("direction", direction)?;
    let depth = validate_positive("depth", depth)?;
    let angle = validate_finite("draft_angle_rad", draft_angle_rad)?;
    if angle.abs() >= std::f64::consts::FRAC_PI_2 - 1e-6 {
        return Err(FeatureError::InvalidInput("draft_angle_rad must be in (-pi/2, pi/2)"));
    }

    let bot: Vec<DVec3> = profile_verts.to_vec();
    let centroid = bot.iter().copied().fold(DVec3::ZERO, |acc, p| acc + p) / bot.len() as f64;
    let axial = dir * depth;
    let taper = depth * angle.tan();

    let top: Vec<DVec3> = bot
        .iter()
        .map(|&p| {
            let v = p - centroid;
            let radial = v - dir * v.dot(dir);
            let radial_dir = if radial.length_squared() > EPS {
                radial.normalize()
            } else {
                DVec3::ZERO
            };
            p + axial + radial_dir * taper
        })
        .collect();

    let tool = build_prism_from_sections(&bot, &top, dir)?;
    Ok(boolean_op(op, target, &tool)?)
}

/// Create a revolution boss/pocket feature from a planar profile.
///
/// The profile polygon is revolved around `axis_origin + t * axis_dir` by
/// `angle_rad`, then combined with `target` by boolean `op`.
///
/// Analogous to OCCT `BRepFeat_MakeRevol` for linear profile faces.
pub fn make_revolution(
    target: &BRep,
    profile_verts: &[DVec3],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle_rad: f64,
    op: BooleanOpType,
) -> Result<BRep, FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::InvalidInput("profile_verts needs >= 3 vertices"));
    }
    if !axis_origin.is_finite() {
        return Err(FeatureError::NonFiniteInput("axis_origin"));
    }
    let axis_dir = normalize("axis_dir", axis_dir)?;
    let angle_rad = validate_positive("angle_rad", angle_rad)?;

    let profile = build_polygon_face_brep(profile_verts)?;
    let tool = rcad_modeling::revolve(&profile, 0, axis_origin, axis_dir, angle_rad)?;

    Ok(boolean_op(op, target, &tool)?)
}

fn build_polygon_face_brep(profile_verts: &[DVec3]) -> Result<BRep, FeatureError> {
    if profile_verts.len() < 3 {
        return Err(FeatureError::InvalidInput("profile_verts needs >= 3 vertices"));
    }

    let n = profile_verts.len();
    let mut brep = BRep {
        vertices: Vec::with_capacity(n),
        edges: Vec::with_capacity(n),
        solids: Vec::new(),
        geom: GeomStore::default(),
        compound: None,
        compsolid: None,
    };

    for &p in profile_verts {
        brep.vertices.push(Vertex { point: p });
    }

    for i in 0..n {
        let j = (i + 1) % n;
        brep.edges.push(Edge { start: i, end: j });
    }

    let normal = {
        let a = profile_verts[0];
        let b = profile_verts[1];
        let c = profile_verts[2];
        let n = (b - a).cross(c - a);
        if n.length_squared() <= EPS {
            return Err(FeatureError::InvalidInput("profile_verts are degenerate"));
        }
        n.normalize()
    };

    let face = Face {
        outer_wire: Wire {
            edges: (0..n).map(WireEdge::fwd).collect(),
        },
        inner_wires: vec![],
        normal,
        triangles: vec![],
        mesh_dirty: true,
    };

    brep.solids.push(Solid {
        shells: vec![Shell { faces: vec![face] }],
    });

    Ok(brep)
}

/// Build a solid BRep prism from a polygon profile (n vertices, coplanar) extruded
fn build_polygon_prism(profile_verts: &[DVec3], dir: DVec3, depth: f64) -> Result<BRep, FeatureError> {
    let bot: Vec<DVec3> = profile_verts.to_vec();
    let top: Vec<DVec3> = bot.iter().map(|&p| p + dir * depth).collect();
    build_prism_from_sections(&bot, &top, dir)
}

fn build_prism_from_sections(bot: &[DVec3], top: &[DVec3], dir: DVec3) -> Result<BRep, FeatureError> {
    let n = bot.len();
    if n < 3 || top.len() != n {
        return Err(FeatureError::InvalidInput("section vertex count mismatch"));
    }

    let mut brep = BRep {
        vertices: Vec::with_capacity(2 * n),
        edges: Vec::new(),
        solids: Vec::new(),
        geom: GeomStore::default(),
        compound: None,
        compsolid: None,
    };

    // Add vertices: bot[0..n] then top[0..n]
    // bot vertex index: i; top vertex index: n + i
    for &p in bot { brep.vertices.push(Vertex { point: p }); }
    for &p in top { brep.vertices.push(Vertex { point: p }); }

    /// Add a line edge from start to end and return its index.
    fn add_line_edge(brep: &mut BRep, start: usize, end: usize) -> usize {
        let p0 = brep.vertices[start].point;
        let p1 = brep.vertices[end].point;
        let d = p1 - p0;
        let len = d.length();
        let dir = if len > 0.0 { d / len } else { DVec3::X };
        let ei = brep.edges.len();
        brep.edges.push(Edge { start, end });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
        brep.geom.edge_curve.push(Some(ci));
        brep.geom.edge_curve_range.push(Some([0.0, len]));
        brep.geom.edge_degenerated.push(false);
        ei
    }

    // Bottom-cap edges: bot[i] -> bot[(i+1)%n]
    let bot_edges: Vec<usize> = (0..n).map(|i| add_line_edge(&mut brep, i, (i + 1) % n)).collect();
    // Top-cap edges: top[i] -> top[(i+1)%n]
    let top_edges: Vec<usize> = (0..n).map(|i| add_line_edge(&mut brep, n + i, n + (i + 1) % n)).collect();
    // Vertical edges: bot[i] -> top[i]
    let vert_edges: Vec<usize> = (0..n).map(|i| add_line_edge(&mut brep, i, n + i)).collect();

    let mut faces = Vec::with_capacity(n + 2);

    // Bottom cap (outward normal = -dir): reverse traversal of bot edges
    {
        let wire_edges: Vec<WireEdge> = (0..n)
            .map(|i| WireEdge { idx: bot_edges[n - 1 - i], forward: false })
            .collect();
        faces.push(Face { outer_wire: Wire { edges: wire_edges }, inner_wires: vec![],
            normal: -dir, triangles: vec![], mesh_dirty: true });
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane { origin: bot[0], normal: -dir }));
        brep.geom.face_surface.push(Some(si));
    }

    // Top cap (outward normal = +dir): forward traversal of top edges
    {
        let wire_edges: Vec<WireEdge> = (0..n)
            .map(|i| WireEdge { idx: top_edges[i], forward: true })
            .collect();
        faces.push(Face { outer_wire: Wire { edges: wire_edges }, inner_wires: vec![],
            normal: dir, triangles: vec![], mesh_dirty: true });
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane { origin: top[0], normal: dir }));
        brep.geom.face_surface.push(Some(si));
    }

    // Lateral quad faces: quad bot[i] -> bot[j] -> top[j] -> top[i] for each edge i
    for i in 0..n {
        let j = (i + 1) % n;
        let a = bot[i];
        let b = bot[j];
        let c = top[j];
        let face_normal = {
            let ab = b - a;
            let ac = c - a;
            let nv = ab.cross(ac);
            if nv.length_squared() > 1e-24 { nv.normalize() } else { -dir.cross(ab).normalize() }
        };
        // wire: bot[i]->bot[j], vert bot[j]->top[j], top[j]->top[i] (reversed), vert top[i]->bot[i] (reversed)
        let wire_edges = vec![
            WireEdge { idx: bot_edges[i],  forward: true },
            WireEdge { idx: vert_edges[j], forward: true },
            WireEdge { idx: top_edges[i],  forward: false },
            WireEdge { idx: vert_edges[i], forward: false },
        ];
        faces.push(Face { outer_wire: Wire { edges: wire_edges }, inner_wires: vec![],
            normal: face_normal, triangles: vec![], mesh_dirty: true });
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane { origin: a, normal: face_normal }));
        brep.geom.face_surface.push(Some(si));
    }

    brep.solids.push(Solid { shells: vec![Shell { faces }] });
    Ok(brep)
}

#[cfg(test)]
mod tests {
    use glam::DVec3;
    use rcad_kernel::{BRep, PrimitiveSolid};

    use super::*;
    #[test]
    fn cylindrical_hole_subtracts_from_box() {
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        });

        let result = make_cylindrical_hole(
            &target,
            DVec3::ZERO,
            DVec3::Y,
            DVec3::X,
            0.6,
            6.0,
        )
        .expect("cylindrical hole should succeed");

        assert!(
            result.solids[0].shells[0].faces.len() >= target.solids[0].shells[0].faces.len(),
            "hole operation should keep or increase face count"
        );
        assert!(!result.edges.is_empty(), "hole result should keep edge topology");
    }

    #[test]
    fn cylindrical_hole_rejects_parallel_axis_ref_dir() {
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        let err = make_cylindrical_hole(
            &target,
            DVec3::ZERO,
            DVec3::Y,
            DVec3::Y,
            0.3,
            3.0,
        )
        .expect_err("parallel axis/ref_dir must be rejected");

        assert!(matches!(err, FeatureError::ParallelVectors(_, _)));
    }

    #[test]
    fn make_prism_boss_adds_material() {
        // Start with a box, add a smaller rectangular prism on top.
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 0.5,
            depth: 4.0,
        });

        // Profile: a 1脳1 square centred at origin on the Z=0 plane
        let profile = vec![
            DVec3::new(-0.5, 0.0, -0.5),
            DVec3::new( 0.5, 0.0, -0.5),
            DVec3::new( 0.5, 0.0,  0.5),
            DVec3::new(-0.5, 0.0,  0.5),
        ];

        let result = make_prism(&target, &profile, DVec3::Y, 1.0, BooleanOpType::Union)
            .expect("make_prism boss should succeed");

        assert!(!result.edges.is_empty(), "boss result must have edges");
        assert!(
            result.solids[0].shells[0].faces.len() >= target.solids[0].shells[0].faces.len(),
            "boss should keep or increase face count"
        );
    }

    #[test]
    fn make_prism_pocket_removes_material() {
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 2.0,
            depth: 4.0,
        });

        // Profile: a 0.5脳0.5 square, extruded 3.0 through the 2.0-tall box
        let profile = vec![
            DVec3::new(-0.25, 0.0, -0.25),
            DVec3::new( 0.25, 0.0, -0.25),
            DVec3::new( 0.25, 0.0,  0.25),
            DVec3::new(-0.25, 0.0,  0.25),
        ];

        let result = make_prism(&target, &profile, DVec3::Y, 3.0, BooleanOpType::Difference)
            .expect("make_prism pocket should succeed");

        assert!(!result.edges.is_empty(), "pocket result must have edges");
    }

    #[test]
    fn make_prism_rejects_degenerate_profile() {
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0, height: 2.0, depth: 2.0,
        });
        let err = make_prism(&target, &[DVec3::ZERO, DVec3::X], DVec3::Y, 1.0, BooleanOpType::Union)
            .expect_err("profile with 2 verts must be rejected");
        assert!(matches!(err, FeatureError::NonPositiveInput(_)));
    }

    #[test]
    fn make_draft_prism_boss_adds_material() {
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 2.0,
            depth: 4.0,
        });
        let profile = vec![
            DVec3::new(-0.4, 1.0, -0.4),
            DVec3::new(0.4, 1.0, -0.4),
            DVec3::new(0.4, 1.0, 0.4),
            DVec3::new(-0.4, 1.0, 0.4),
        ];

        let out = make_draft_prism(
            &target,
            &profile,
            DVec3::Y,
            0.8,
            8.0_f64.to_radians(),
            BooleanOpType::Union,
        )
        .expect("draft prism boss should succeed");

        assert!(!out.solids.is_empty());
        assert!(out.edges.len() >= target.edges.len());
    }

    #[test]
    fn make_revolution_boss_adds_material() {
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 5.0,
            height: 5.0,
            depth: 5.0,
        });
        let profile = vec![
            DVec3::new(1.5, -0.3, 0.0),
            DVec3::new(2.0, -0.3, 0.0),
            DVec3::new(2.0, 0.3, 0.0),
            DVec3::new(1.5, 0.3, 0.0),
        ];

        let out = make_revolution(
            &target,
            &profile,
            DVec3::ZERO,
            DVec3::Z,
            std::f64::consts::TAU,
            BooleanOpType::Union,
        )
        .expect("make_revolution boss should succeed");

        assert!(!out.solids.is_empty());
        assert!(!out.edges.is_empty());
    }

    #[test]
    fn make_revolution_rejects_invalid_profile() {
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        let err = make_revolution(
            &target,
            &[DVec3::ZERO, DVec3::X],
            DVec3::ZERO,
            DVec3::Z,
            1.0,
            BooleanOpType::Union,
        )
        .expect_err("profile with <3 verts should fail");
        assert!(matches!(err, FeatureError::InvalidInput(_)));
    }
}

// ─── SplitShape: split a face by a cutting wire ──────────────────────────────

/// Error returned by [`split_face_by_wire`].
#[derive(Debug)]
pub enum SplitShapeError {
    FaceNotFound,
    CutPathTooShort,
    CutVertexNotOnWire { vertex_idx: usize },
    CutPathClosedLoop,
    DegenerateResult,
}

impl std::fmt::Display for SplitShapeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::FaceNotFound => write!(f, "face index out of range"),
            Self::CutPathTooShort => write!(f, "cut path needs at least 2 vertices"),
            Self::CutVertexNotOnWire { vertex_idx } => {
                write!(f, "cut vertex {vertex_idx} is not on the face outer wire")
            }
            Self::CutPathClosedLoop => write!(f, "cut path start and end are the same wire vertex"),
            Self::DegenerateResult => write!(f, "split produced a degenerate wire"),
        }
    }
}

impl std::error::Error for SplitShapeError {}

/// Split a face by a cutting wire (path of vertex indices already in `brep.vertices`).
///
/// - `cut_path`: at least 2 vertex indices; first and last must appear as the
///   *start* vertex of some edge in the face's outer wire.
/// - New line edges are inserted for each segment of `cut_path`.
/// - The face is replaced by two sub-faces.
///
/// Analogous to OCCT `BRepFeat_SplitShape`.
pub fn split_face_by_wire(
    brep: &mut BRep,
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    cut_path: &[usize],
) -> Result<usize, SplitShapeError> {
    if cut_path.len() < 2 {
        return Err(SplitShapeError::CutPathTooShort);
    }
    if solid_idx >= brep.solids.len()
        || shell_idx >= brep.solids[solid_idx].shells.len()
        || face_idx >= brep.solids[solid_idx].shells[shell_idx].faces.len()
    {
        return Err(SplitShapeError::FaceNotFound);
    }

    let start_v = cut_path[0];
    let end_v = *cut_path.last().unwrap();

    let outer_edges = brep.solids[solid_idx].shells[shell_idx].faces[face_idx]
        .outer_wire.edges.clone();

    // Build ordered vertex sequence of the outer wire.
    let wire_verts: Vec<usize> = outer_edges.iter().map(|we| {
        let e = &brep.edges[we.idx];
        if we.forward { e.start } else { e.end }
    }).collect();

    let pos_start = wire_verts.iter().position(|&v| v == start_v)
        .ok_or(SplitShapeError::CutVertexNotOnWire { vertex_idx: start_v })?;
    let pos_end = wire_verts.iter().position(|&v| v == end_v)
        .ok_or(SplitShapeError::CutVertexNotOnWire { vertex_idx: end_v })?;
    if pos_start == pos_end {
        return Err(SplitShapeError::CutPathClosedLoop);
    }

    let n = outer_edges.len();

    // Add line edges for the cut path segments.
    let cut_edge_indices: Vec<usize> = cut_path.windows(2).map(|w| {
        let (sv, ev) = (w[0], w[1]);
        let ei = brep.edges.len();
        let p0 = brep.vertices[sv].point;
        let p1 = brep.vertices[ev].point;
        let d = p1 - p0;
        let len = d.length();
        let dir = if len > 1e-30 { d / len } else { DVec3::X };
        brep.edges.push(Edge { start: sv, end: ev });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
        brep.geom.edge_curve.push(Some(ci));
        brep.geom.edge_curve_range.push(Some([0.0, len]));
        brep.geom.edge_degenerated.push(false);
        ei
    }).collect();

    // Half A: outer[pos_start..pos_end] + cut forward.
    let half_a: Vec<WireEdge> = (0..(pos_end - pos_start))
        .map(|i| outer_edges[(pos_start + i) % n]).collect();
    let cut_fwd: Vec<WireEdge> = cut_edge_indices.iter().map(|&ei| WireEdge::fwd(ei)).collect();
    let mut wire_a = half_a;
    wire_a.extend_from_slice(&cut_fwd);

    // Half B: outer[pos_end..] + outer[..pos_start] + cut reversed.
    let half_b_len = n - (pos_end - pos_start);
    let half_b: Vec<WireEdge> = (0..half_b_len)
        .map(|i| outer_edges[(pos_end + i) % n]).collect();
    let cut_rev: Vec<WireEdge> = cut_edge_indices.iter().rev().map(|&ei| WireEdge::rev(ei)).collect();
    let mut wire_b = half_b;
    wire_b.extend_from_slice(&cut_rev);

    if wire_a.len() < 3 || wire_b.len() < 3 {
        return Err(SplitShapeError::DegenerateResult);
    }

    let orig_normal = brep.solids[solid_idx].shells[shell_idx].faces[face_idx].normal;
    let orig_inner = brep.solids[solid_idx].shells[shell_idx].faces[face_idx].inner_wires.clone();

    let face_a = Face {
        outer_wire: Wire { edges: wire_a },
        inner_wires: orig_inner.clone(),
        normal: orig_normal,
        triangles: vec![],
        mesh_dirty: true,
    };
    let face_b = Face {
        outer_wire: Wire { edges: wire_b },
        inner_wires: orig_inner,
        normal: orig_normal,
        triangles: vec![],
        mesh_dirty: true,
    };

    // Update GeomStore face_surface flat index.
    let flat_idx: usize = brep.solids[..solid_idx].iter()
        .flat_map(|s| s.shells.iter()).map(|sh| sh.faces.len()).sum::<usize>()
        + brep.solids[solid_idx].shells[..shell_idx].iter()
            .map(|sh| sh.faces.len()).sum::<usize>()
        + face_idx;
    let orig_surf = brep.geom.face_surface.get(flat_idx).copied().flatten();
    if flat_idx + 1 <= brep.geom.face_surface.len() {
        brep.geom.face_surface.insert(flat_idx + 1, orig_surf);
    }
    if flat_idx + 1 <= brep.geom.face_tolerance.len() {
        let ft = brep.geom.face_tolerance[flat_idx];
        brep.geom.face_tolerance.insert(flat_idx + 1, ft);
    }

    brep.solids[solid_idx].shells[shell_idx].faces[face_idx] = face_a;
    brep.solids[solid_idx].shells[shell_idx].faces.insert(face_idx + 1, face_b);
    Ok(1)
}

// ─── Linear rib / slot ───────────────────────────────────────────────────────

/// Create a linear rib (or slot) feature via prism boolean.
///
/// Analogous to OCCT `BRepFeat_MakeLinearForm`.
pub fn make_linear_rib(
    target: &BRep,
    profile_verts: &[DVec3],
    direction: DVec3,
    depth: f64,
    op: BooleanOpType,
) -> Result<BRep, FeatureError> {
    make_prism(target, profile_verts, direction, depth, op)
}

/// Create a revolution rib/slot feature via revolve boolean.
///
/// Analogous to OCCT `BRepFeat_MakeRevolutionForm`.
pub fn make_revolution_rib(
    target: &BRep,
    profile_verts: &[DVec3],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle_rad: f64,
    op: BooleanOpType,
) -> Result<BRep, FeatureError> {
    make_revolution(target, profile_verts, axis_origin, axis_dir, angle_rad, op)
}

#[cfg(test)]
mod split_tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

    fn make_square_brep() -> BRep {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 3
        brep.edges.push(Edge { start: 0, end: 1 }); // e0
        brep.edges.push(Edge { start: 1, end: 2 }); // e1
        brep.edges.push(Edge { start: 2, end: 3 }); // e2
        brep.edges.push(Edge { start: 3, end: 0 }); // e3
        brep.solids.push(Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![
                            WireEdge::fwd(0), WireEdge::fwd(1),
                            WireEdge::fwd(2), WireEdge::fwd(3),
                        ],
                    },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                }],
            }],
        });
        brep
    }

    #[test]
    fn split_face_by_wire_splits_square_into_two_triangles() {
        let mut brep = make_square_brep();
        // Cut from v0 (wire pos 0) to v2 (wire pos 2).
        let result = split_face_by_wire(&mut brep, 0, 0, 0, &[0, 2]);
        assert!(result.is_ok(), "split should succeed: {:?}", result);
        let shell = &brep.solids[0].shells[0];
        assert_eq!(shell.faces.len(), 2, "face should be split into 2");
        for face in &shell.faces {
            assert_eq!(face.outer_wire.edges.len(), 3, "each sub-face should be a triangle");
        }
    }

    #[test]
    fn split_face_by_wire_rejects_missing_vertex() {
        let mut brep = make_square_brep();
        let result = split_face_by_wire(&mut brep, 0, 0, 0, &[0, 99]);
        assert!(matches!(result, Err(SplitShapeError::CutVertexNotOnWire { .. })));
    }

    #[test]
    fn make_linear_rib_creates_boss() {
        let target = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0, height: 2.0, depth: 0.5,
        });
        let profile = vec![
            DVec3::new(0.5, 0.0, 0.0),
            DVec3::new(1.5, 0.0, 0.0),
            DVec3::new(1.5, 0.0, 1.0),
            DVec3::new(0.5, 0.0, 1.0),
        ];
        let result = make_linear_rib(&target, &profile, DVec3::Y, 0.2, BooleanOpType::Union);
        assert!(result.is_ok(), "linear rib should succeed: {:?}", result);
    }
}
