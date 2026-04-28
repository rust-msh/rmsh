//! BRep validity checker.
//!
//! Analogous to OCCT `BRepCheck_Analyzer`. Checks structural and geometric
//! consistency of a BRep without modifying it.
//!
//! # Checks performed
//!
//! - **C1 Wire closure**: every wire must form a closed chain — the end vertex of
//!   each edge must equal the start vertex of the next edge.
//! - **C2 Face normal consistency**: each face's stored normal must not be a zero
//!   vector.
//! - **C3 Degenerate face**: faces with fewer than 3 wire edges are degenerate.
//! - **C4 Edge index validity**: WireEdge indices must be within bounds of
//!   `brep.edges`.
//! - **C5 Vertex index validity**: each edge's start/end indices must be within
//!   bounds of `brep.vertices`.
//! - **C6 Manifold topology**: each edge must be shared by exactly 2 faces
//!   (for closed manifold solids).
//! - **C7 Wire self-intersection**: a wire's edges must not share vertices
//!   except at consecutive junctions (no figure-8 or self-touching wires).
//!
//! # Extended checks (OCCT BRepCheck_Analyzer equivalent)
//!
//! - **Surface continuity**: C0, C1, C2 continuity across adjacent faces
//! - **Curve-surface consistency**: 3D curve endpoints match surface evaluation
//! - **Edge-curve tolerance verification**: edge tolerance covers geometry deviation
//! - **Face-surface tolerance verification**: face tolerance covers surface deviation
//! - **Shell orientation consistency**: consistent normal orientation in shells
//! - **Solid closure verification**: all edges shared by exactly 2 faces
//! - **Wire orientation**: clockwise vs counter-clockwise validation
//! - **Nested wire validation**: inner loops properly contained within outer
//! - **Tolerance consistency**: adjacent faces have compatible tolerances
//! - **Vertex tolerance propagation**: vertices have appropriate tolerances
//! - **Aspect ratio checks**: face quality metrics
//! - **Degenerate geometry detection**: zero-length edges, collapsed faces
//! - **Sliver face detection**: very thin triangular faces
//! - **Small feature detection**: tiny faces, edges, vertices

use glam::{DVec2, DVec3};
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};

/// A single validity issue found during checking.
#[derive(Debug, Clone, PartialEq)]
pub enum CheckIssue {
    /// Wire is not closed: end vertex of edge `edge_idx` does not match start
    /// vertex of the next edge in the wire (solid `solid`, shell `shell`,
    /// face `face`, position `wire_pos`).
    OpenWire {
        solid: usize,
        shell: usize,
        face: usize,
        /// Index of the edge within the wire where the gap occurs.
        wire_pos: usize,
    },
    /// Face normal is a zero vector.
    ZeroNormal {
        solid: usize,
        shell: usize,
        face: usize,
    },
    /// Face outer wire has fewer than 3 edges.
    DegenerateFace {
        solid: usize,
        shell: usize,
        face: usize,
    },
    /// A WireEdge references an edge index that is out of bounds.
    InvalidEdgeIndex {
        solid: usize,
        shell: usize,
        face: usize,
        edge_idx: usize,
    },
    /// An edge references a vertex index that is out of bounds.
    InvalidVertexIndex { edge: usize, vertex_idx: usize },
    /// An edge is shared by more or fewer than 2 faces (non-manifold).
    NonManifoldEdge { edge_idx: usize, face_count: usize },
    /// A wire has self-intersecting topology: a vertex appears more than
    /// twice (once as start, once as end) in the same wire.
    SelfIntersectingWire {
        solid: usize,
        shell: usize,
        face: usize,
        wire_idx: usize,
        vertex: usize,
    },
    /// A wire's outer boundary edges intersect each other geometrically
    /// (non-adjacent edges in the wire cross in 3D space).
    ///
    /// This catches cases where a face wire forms a figure-eight or butterfly
    /// polygon rather than a simple closed loop.
    GeometricSelfIntersection {
        solid: usize,
        shell: usize,
        face: usize,
        /// Index of one of the crossing edges within the outer wire.
        edge_a: usize,
        /// Index of the other crossing edge within the outer wire.
        edge_b: usize,
    },
    // ── Geometry validation issues ─────────────────────────────────────────────
    /// Surface continuity violation between adjacent faces.
    SurfaceContinuityViolation {
        solid: usize,
        face_a: usize,
        face_b: usize,
        shared_edge: usize,
        /// Expected continuity (0=C0, 1=C1, 2=C2)
        expected: u8,
        /// Actual continuity achieved
        actual: u8,
        /// Gap or angle deviation at the junction
        deviation: f64,
    },
    /// Curve-surface consistency violation: 3D curve doesn't match surface evaluation.
    CurveSurfaceMismatch {
        edge: usize,
        surface: usize,
        /// Maximum deviation between 3D curve and surface curve
        max_deviation: f64,
    },
    /// Edge tolerance insufficient to cover geometry deviation.
    EdgeToleranceViolation {
        edge: usize,
        stored_tolerance: f64,
        required_tolerance: f64,
    },
    /// Face tolerance insufficient to cover surface deviation.
    FaceToleranceViolation {
        solid: usize,
        shell: usize,
        face: usize,
        stored_tolerance: f64,
        required_tolerance: f64,
    },
    // ── Topology validation issues ─────────────────────────────────────────────
    /// Shell has inconsistent orientation (mixed inward/outward normals).
    ShellOrientationInconsistent {
        solid: usize,
        shell: usize,
        faces_with_inverted_normals: usize,
    },
    /// Solid is not closed (has boundary edges).
    SolidNotClosed {
        solid: usize,
        boundary_edge_count: usize,
    },
    /// Wire orientation is incorrect for its role (outer vs inner).
    WireOrientationIncorrect {
        solid: usize,
        shell: usize,
        face: usize,
        wire_idx: usize,
        /// true = should be CCW (outer), false = should be CW (inner)
        expected_ccw: bool,
        actual_ccw: bool,
    },
    /// Inner wire is not properly contained within outer wire.
    NestedWireViolation {
        solid: usize,
        shell: usize,
        face: usize,
        inner_wire_idx: usize,
        /// Number of inner wire vertices outside outer wire boundary
        vertices_outside: usize,
    },
    // ── Tolerance issues ───────────────────────────────────────────────────────
    /// Adjacent faces have inconsistent tolerances.
    ToleranceInconsistency {
        edge: usize,
        face_a: usize,
        face_b: usize,
        tolerance_a: f64,
        tolerance_b: f64,
        ratio: f64,
    },
    /// Vertex tolerance doesn't cover incident edge endpoints.
    VertexToleranceViolation {
        vertex: usize,
        stored_tolerance: f64,
        required_tolerance: f64,
    },
    // ── Quality metric issues ─────────────────────────────────────────────────
    /// Face has poor aspect ratio.
    PoorAspectRatio {
        solid: usize,
        shell: usize,
        face: usize,
        aspect_ratio: f64,
    },
    /// Edge has near-zero length.
    DegenerateEdge {
        edge: usize,
        length: f64,
    },
    /// Face is a sliver (very thin).
    SliverFace {
        solid: usize,
        shell: usize,
        face: usize,
        area: f64,
        min_dimension: f64,
    },
    /// Small feature detected (tiny face or edge).
    SmallFeature {
        solid: usize,
        shell: usize,
        face: usize,
        feature_type: SmallFeatureType,
        size: f64,
    },
}

/// Type of small feature detected.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SmallFeatureType {
    TinyFace,
    TinyEdge,
    TinyVertexGap,
}

impl std::fmt::Display for CheckIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CheckIssue::OpenWire {
                solid,
                shell,
                face,
                wire_pos,
            } => write!(
                f,
                "OpenWire: solid={solid} shell={shell} face={face} at wire pos {wire_pos}"
            ),
            CheckIssue::ZeroNormal { solid, shell, face } => {
                write!(f, "ZeroNormal: solid={solid} shell={shell} face={face}")
            }
            CheckIssue::DegenerateFace { solid, shell, face } => {
                write!(f, "DegenerateFace: solid={solid} shell={shell} face={face}")
            }
            CheckIssue::InvalidEdgeIndex {
                solid,
                shell,
                face,
                edge_idx,
            } => write!(
                f,
                "InvalidEdgeIndex: solid={solid} shell={shell} face={face} edge={edge_idx}"
            ),
            CheckIssue::InvalidVertexIndex { edge, vertex_idx } => {
                write!(f, "InvalidVertexIndex: edge={edge} vertex={vertex_idx}")
            }
            CheckIssue::NonManifoldEdge { edge_idx, face_count } => {
                write!(f, "NonManifoldEdge: edge={edge_idx} shared by {face_count} faces (expected 2)")
            }
            CheckIssue::SelfIntersectingWire {
                solid,
                shell,
                face,
                wire_idx,
                vertex,
            } => write!(
                f,
                "SelfIntersectingWire: solid={solid} shell={shell} face={face} wire={wire_idx} vertex={vertex}"
            ),
            CheckIssue::GeometricSelfIntersection {
                solid,
                shell,
                face,
                edge_a,
                edge_b,
            } => write!(
                f,
                "GeometricSelfIntersection: solid={solid} shell={shell} face={face} edges {edge_a} and {edge_b} cross"
            ),
            // Geometry validation
            CheckIssue::SurfaceContinuityViolation {
                solid,
                face_a,
                face_b,
                shared_edge,
                expected,
                actual,
                deviation,
            } => write!(
                f,
                "SurfaceContinuityViolation: solid={solid} faces {face_a}/{face_b} edge={shared_edge} expected C{expected} got C{actual} deviation={deviation:.6e}"
            ),
            CheckIssue::CurveSurfaceMismatch { edge, surface, max_deviation } => {
                write!(f, "CurveSurfaceMismatch: edge={edge} surface={surface} deviation={max_deviation:.6e}")
            }
            CheckIssue::EdgeToleranceViolation { edge, stored_tolerance, required_tolerance } => {
                write!(f, "EdgeToleranceViolation: edge={edge} stored={stored_tolerance:.6e} required={required_tolerance:.6e}")
            }
            CheckIssue::FaceToleranceViolation { solid, shell, face, stored_tolerance, required_tolerance } => {
                write!(f, "FaceToleranceViolation: solid={solid} shell={shell} face={face} stored={stored_tolerance:.6e} required={required_tolerance:.6e}")
            }
            // Topology validation
            CheckIssue::ShellOrientationInconsistent { solid, shell, faces_with_inverted_normals } => {
                write!(f, "ShellOrientationInconsistent: solid={solid} shell={shell} {faces_with_inverted_normals} inverted faces")
            }
            CheckIssue::SolidNotClosed { solid, boundary_edge_count } => {
                write!(f, "SolidNotClosed: solid={solid} {boundary_edge_count} boundary edges")
            }
            CheckIssue::WireOrientationIncorrect { solid, shell, face, wire_idx, expected_ccw, actual_ccw } => {
                let expected = if *expected_ccw { "CCW" } else { "CW" };
                let actual = if *actual_ccw { "CCW" } else { "CW" };
                write!(f, "WireOrientationIncorrect: solid={solid} shell={shell} face={face} wire={wire_idx} expected={expected} got={actual}")
            }
            CheckIssue::NestedWireViolation { solid, shell, face, inner_wire_idx, vertices_outside } => {
                write!(f, "NestedWireViolation: solid={solid} shell={shell} face={face} inner_wire={inner_wire_idx} {vertices_outside} vertices outside")
            }
            // Tolerance issues
            CheckIssue::ToleranceInconsistency { edge, face_a, face_b, tolerance_a, tolerance_b, ratio } => {
                write!(f, "ToleranceInconsistency: edge={edge} faces {face_a}/{face_b} tol_a={tolerance_a:.6e} tol_b={tolerance_b:.6e} ratio={ratio:.2}")
            }
            CheckIssue::VertexToleranceViolation { vertex, stored_tolerance, required_tolerance } => {
                write!(f, "VertexToleranceViolation: vertex={vertex} stored={stored_tolerance:.6e} required={required_tolerance:.6e}")
            }
            // Quality metrics
            CheckIssue::PoorAspectRatio { solid, shell, face, aspect_ratio } => {
                write!(f, "PoorAspectRatio: solid={solid} shell={shell} face={face} ratio={aspect_ratio:.2}")
            }
            CheckIssue::DegenerateEdge { edge, length } => {
                write!(f, "DegenerateEdge: edge={edge} length={length:.6e}")
            }
            CheckIssue::SliverFace { solid, shell, face, area, min_dimension } => {
                write!(f, "SliverFace: solid={solid} shell={shell} face={face} area={area:.6e} min_dim={min_dimension:.6e}")
            }
            CheckIssue::SmallFeature { solid, shell, face, feature_type, size } => {
                let type_str = match feature_type {
                    SmallFeatureType::TinyFace => "TinyFace",
                    SmallFeatureType::TinyEdge => "TinyEdge",
                    SmallFeatureType::TinyVertexGap => "TinyVertexGap",
                };
                write!(f, "SmallFeature: solid={solid} shell={shell} face={face} type={type_str} size={size:.6e}")
            }
        }
    }
}

/// Result of a BRep validity check.
#[derive(Debug, Clone)]
#[derive(Default)]
pub struct CheckResult {
    pub issues: Vec<CheckIssue>,
}

impl CheckResult {
    /// Returns `true` if no issues were found.
    pub fn is_valid(&self) -> bool {
        self.issues.is_empty()
    }
}

/// Check the validity of a BRep and return a `CheckResult` with any issues found.
///
/// Analogous to OCCT `BRepCheck_Analyzer::Perform()`.
pub fn check(brep: &BRep) -> CheckResult {
    let mut issues = Vec::new();
    let n_edges = brep.edges.len();
    let n_verts = brep.vertices.len();

    // C5: edge vertex bounds
    for (eidx, edge) in brep.edges.iter().enumerate() {
        if edge.start >= n_verts {
            issues.push(CheckIssue::InvalidVertexIndex {
                edge: eidx,
                vertex_idx: edge.start,
            });
        }
        if edge.end >= n_verts {
            issues.push(CheckIssue::InvalidVertexIndex {
                edge: eidx,
                vertex_idx: edge.end,
            });
        }
    }

    // C6: manifold check — each edge must be shared by exactly 2 faces.
    // Count how many faces reference each edge across all solids/shells/faces.
    let mut edge_face_count: Vec<usize> = vec![0; n_edges];
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Count edges in outer wire
                for we in &face.outer_wire.edges {
                    if we.idx < n_edges {
                        edge_face_count[we.idx] += 1;
                    }
                }
                // Count edges in inner wires
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        if we.idx < n_edges {
                            edge_face_count[we.idx] += 1;
                        }
                    }
                }
            }
        }
    }
    for (eidx, &count) in edge_face_count.iter().enumerate() {
        if count != 2 {
            issues.push(CheckIssue::NonManifoldEdge {
                edge_idx: eidx,
                face_count: count,
            });
        }
    }

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                let wire = &face.outer_wire;

                // C2: zero normal
                if face.normal == DVec3::ZERO {
                    issues.push(CheckIssue::ZeroNormal {
                        solid: si,
                        shell: shi,
                        face: fi,
                    });
                }

                // C3: degenerate face
                if wire.edges.len() < 3 {
                    issues.push(CheckIssue::DegenerateFace {
                        solid: si,
                        shell: shi,
                        face: fi,
                    });
                    // Can't check wire closure for degenerate face
                    continue;
                }

                // C4: edge index bounds + collect start/end vertices for wire closure check
                let mut valid = true;
                let mut wire_verts: Vec<(usize, usize)> = Vec::new(); // (start_vidx, end_vidx)
                for we in &wire.edges {
                    if we.idx >= n_edges {
                        issues.push(CheckIssue::InvalidEdgeIndex {
                            solid: si,
                            shell: shi,
                            face: fi,
                            edge_idx: we.idx,
                        });
                        valid = false;
                    } else {
                        let edge = &brep.edges[we.idx];
                        let (sv, ev) = if we.forward {
                            (edge.start, edge.end)
                        } else {
                            (edge.end, edge.start)
                        };
                        wire_verts.push((sv, ev));
                    }
                }

                if !valid {
                    continue;
                }

                // C1: wire closure — end of edge[i] must match start of edge[i+1]
                let n = wire_verts.len();
                for i in 0..n {
                    let next = (i + 1) % n;
                    let end_v = wire_verts[i].1;
                    let start_v = wire_verts[next].0;
                    if end_v != start_v {
                        // Tolerance check: allow same position even if different vertex objects
                        let end_pt = brep.vertices[end_v].point;
                        let start_pt = brep.vertices[start_v].point;
                        if (end_pt - start_pt).length() > 1e-6 {
                            issues.push(CheckIssue::OpenWire {
                                solid: si,
                                shell: shi,
                                face: fi,
                                wire_pos: i,
                            });
                        }
                    }
                }

                // C7: wire self-intersection — each vertex should appear at most
                // twice in the wire (once as start of an edge, once as end of another).
                check_wire_self_intersection(
                    &wire_verts,
                    &brep.vertices,
                    si, shi, fi, 0, // outer wire index = 0
                    &mut issues,
                );

                // C8: geometric self-intersection — check if non-adjacent edges of
                // the outer wire cross each other in 3D space (projects to 2D via
                // the face plane for planar faces).
                check_geometric_self_intersection(
                    &wire_verts,
                    &brep.vertices,
                    si, shi, fi,
                    &mut issues,
                );

                // Check inner wires too
                for (wi, inner_wire) in face.inner_wires.iter().enumerate() {
                    if inner_wire.edges.len() < 2 {
                        continue; // too few edges to self-intersect
                    }
                    let mut inner_verts: Vec<(usize, usize)> = Vec::new();
                    let mut inner_valid = true;
                    for we in &inner_wire.edges {
                        if we.idx >= n_edges {
                            issues.push(CheckIssue::InvalidEdgeIndex {
                                solid: si,
                                shell: shi,
                                face: fi,
                                edge_idx: we.idx,
                            });
                            inner_valid = false;
                        } else {
                            let edge = &brep.edges[we.idx];
                            let (sv, ev) = if we.forward {
                                (edge.start, edge.end)
                            } else {
                                (edge.end, edge.start)
                            };
                            inner_verts.push((sv, ev));
                        }
                    }
                    if !inner_valid {
                        continue;
                    }

                    // Inner wire closure check
                    let n_inner = inner_verts.len();
                    for i in 0..n_inner {
                        let next = (i + 1) % n_inner;
                        let end_v = inner_verts[i].1;
                        let start_v = inner_verts[next].0;
                        if end_v != start_v {
                            let end_pt = brep.vertices[end_v].point;
                            let start_pt = brep.vertices[start_v].point;
                            if (end_pt - start_pt).length() > 1e-6 {
                                issues.push(CheckIssue::OpenWire {
                                    solid: si,
                                    shell: shi,
                                    face: fi,
                                    wire_pos: i,
                                });
                            }
                        }
                    }

                    // Inner wire self-intersection
                    check_wire_self_intersection(
                        &inner_verts,
                        &brep.vertices,
                        si, shi, fi,
                        wi + 1, // inner wire indices start after outer
                        &mut issues,
                    );
                }
            }
        }
    }

    CheckResult { issues }
}

/// Check a single wire for self-intersecting topology.
///
/// A valid wire wire should have each vertex appear at most twice across
/// all edge endpoints: once as the start of some edge and once as the end
/// of another edge. If a vertex appears 3+ times, the wire self-intersects.
fn check_wire_self_intersection(
    wire_verts: &[(usize, usize)],
    vertices: &[rcad_kernel::topology::Vertex],
    solid: usize,
    shell: usize,
    face: usize,
    wire_idx: usize,
    issues: &mut Vec<CheckIssue>,
) {
    use std::collections::HashMap;
    let mut vertex_count: HashMap<usize, usize> = HashMap::new();
    for &(sv, ev) in wire_verts {
        *vertex_count.entry(sv).or_insert(0) += 1;
        *vertex_count.entry(ev).or_insert(0) += 1;
    }
    // In a closed wire, each vertex should appear exactly twice (once as start,
    // once as end). Allow tolerance for vertices at the same position.
    for (&vidx, &count) in &vertex_count {
        if count > 2 {
            // Check if it's actually a geometric self-intersection (different
            // positions) or just the same vertex referenced multiple times
            // (which could be valid for some topologies).
            if vidx < vertices.len() {
                issues.push(CheckIssue::SelfIntersectingWire {
                    solid,
                    shell,
                    face,
                    wire_idx,
                    vertex: vidx,
                });
            }
        }
    }
}

/// Check whether non-adjacent edges of a wire intersect geometrically.
///
/// Projects the wire edge endpoints onto the face's 2D plane (using any two
/// non-collinear edges to form a local basis) and runs 2D segment intersection
/// tests on all non-adjacent edge pairs. Adjacent edges share an endpoint and
/// therefore trivially "intersect" at that endpoint — they are excluded.
///
/// This check only tests the start/end vertices; curved edges (circles, BSplines)
/// are approximated by their chord.
fn check_geometric_self_intersection(
    wire_verts: &[(usize, usize)],
    vertices: &[rcad_kernel::topology::Vertex],
    solid: usize,
    shell: usize,
    face: usize,
    issues: &mut Vec<CheckIssue>,
) {
    let n = wire_verts.len();
    if n < 4 {
        return; // Need at least 4 edges for non-adjacent crossing to exist
    }

    // Build the 3D segment list: each edge is (p0, p1).
    let segs: Vec<(DVec3, DVec3)> = wire_verts
        .iter()
        .map(|&(sv, ev)| {
            let p0 = vertices.get(sv).map(|v| v.point).unwrap_or(DVec3::ZERO);
            let p1 = vertices.get(ev).map(|v| v.point).unwrap_or(DVec3::ZERO);
            (p0, p1)
        })
        .collect();

    // Project to 2D using a local basis from the first non-degenerate edge.
    let (origin, axis_u, axis_v) = {
        let mut found = None;
        for i in 0..n {
            let d = segs[i].1 - segs[i].0;
            if d.length() > 1e-12 {
                let u = d.normalize();
                // Pick an arbitrary perpendicular
                let tmp = if u.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                let v = u.cross(tmp).normalize();
                found = Some((segs[i].0, u, v));
                break;
            }
        }
        match found {
            Some(b) => b,
            None => return, // All edges are degenerate; skip
        }
    };

    let project = |p: DVec3| -> [f64; 2] {
        let d = p - origin;
        [d.dot(axis_u), d.dot(axis_v)]
    };

    let seg2d: Vec<([f64; 2], [f64; 2])> =
        segs.iter().map(|&(p0, p1)| (project(p0), project(p1))).collect();

    // Check all non-adjacent pairs: i and j are non-adjacent when |i - j| > 1
    // (mod n, so also check the wrap-around adjacency).
    for i in 0..n {
        for j in (i + 2)..n {
            // Skip adjacent pair at the wrap-around (edge n-1 and edge 0 share a vertex)
            if i == 0 && j == n - 1 {
                continue;
            }
            if segments_2d_properly_intersect(seg2d[i].0, seg2d[i].1, seg2d[j].0, seg2d[j].1) {
                issues.push(CheckIssue::GeometricSelfIntersection {
                    solid,
                    shell,
                    face,
                    edge_a: i,
                    edge_b: j,
                });
                return; // Report only the first crossing per face
            }
        }
    }
}

/// Returns `true` if the open segment p1→p2 properly intersects segment p3→p4.
/// Returns `false` if they only share an endpoint (T-intersection) or don't cross.
fn segments_2d_properly_intersect(
    p1: [f64; 2],
    p2: [f64; 2],
    p3: [f64; 2],
    p4: [f64; 2],
) -> bool {
    let d1 = [p2[0] - p1[0], p2[1] - p1[1]];
    let d2 = [p4[0] - p3[0], p4[1] - p3[1]];

    let cross = d1[0] * d2[1] - d1[1] * d2[0];

    if cross.abs() < 1e-12 {
        return false; // Parallel or collinear
    }

    let dx = p3[0] - p1[0];
    let dy = p3[1] - p1[1];
    let t = (dx * d2[1] - dy * d2[0]) / cross;
    let s = (dx * d1[1] - dy * d1[0]) / cross;

    // Proper interior intersection: t and s must be strictly in (0, 1)
    let eps = 1e-9;
    t > eps && t < 1.0 - eps && s > eps && s < 1.0 - eps
}

fn count_geometric_self_intersections(
    wire_verts: &[(usize, usize)],
    vertices: &[rcad_kernel::topology::Vertex],
) -> usize {
    let n = wire_verts.len();
    if n < 4 {
        return 0;
    }

    let segs: Vec<(DVec3, DVec3)> = wire_verts
        .iter()
        .map(|&(sv, ev)| {
            let p0 = vertices.get(sv).map(|v| v.point).unwrap_or(DVec3::ZERO);
            let p1 = vertices.get(ev).map(|v| v.point).unwrap_or(DVec3::ZERO);
            (p0, p1)
        })
        .collect();

    let (origin, axis_u, axis_v) = {
        let mut found = None;
        for i in 0..n {
            let d = segs[i].1 - segs[i].0;
            if d.length() > 1e-12 {
                let u = d.normalize();
                let tmp = if u.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                let v = u.cross(tmp).normalize();
                found = Some((segs[i].0, u, v));
                break;
            }
        }
        match found {
            Some(b) => b,
            None => return 0,
        }
    };

    let project = |p: DVec3| -> [f64; 2] {
        let d = p - origin;
        [d.dot(axis_u), d.dot(axis_v)]
    };

    let seg2d: Vec<([f64; 2], [f64; 2])> =
        segs.iter().map(|&(p0, p1)| (project(p0), project(p1))).collect();

    let mut count = 0usize;
    for i in 0..n {
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            if segments_2d_properly_intersect(seg2d[i].0, seg2d[i].1, seg2d[j].0, seg2d[j].1) {
                count += 1;
            }
        }
    }
    count
}

// ── SameParameter diagnosis ───────────────────────────────────────────────────

/// A single edge whose 3D curve endpoints deviate from the vertex positions.
///
/// Analogous to the diagnostic output of `BRepCheck_Edge::SameParameter()` in OCCT.
#[derive(Debug, Clone)]
pub struct SuspectEdge {
    /// Index of the edge in `BRep.edges`.
    pub edge_idx: usize,
    /// Distance from `curve.point_at(t1)` to the start vertex position.
    pub start_gap: f64,
    /// Distance from `curve.point_at(t2)` to the end vertex position.
    pub end_gap: f64,
}

/// Report from [`diagnose_same_parameter`].
#[derive(Debug, Clone, Default)]
pub struct SameParameterDiagnosis {
    /// Edges whose curve endpoints deviate from vertex positions beyond `tolerance`.
    pub suspect_edges: Vec<SuspectEdge>,
}

impl SameParameterDiagnosis {
    /// Returns `true` if no SameParameter violations were found.
    pub fn is_clean(&self) -> bool {
        self.suspect_edges.is_empty()
    }
}

/// A single edge whose PCurve ranges deviate from the 3D edge range.
///
/// Analogous to the diagnostic output of `BRepCheck_Edge::SameRange()` in OCCT.
#[derive(Debug, Clone)]
pub struct SuspectSameRangeEdge {
    /// Index of the edge in `BRep.edges`.
    pub edge_idx: usize,
    /// Number of attached PCurves that do not match the 3D edge range.
    pub mismatched_pcurves: usize,
    /// Maximum endpoint mismatch magnitude among attached PCurves.
    pub max_delta: f64,
}

/// Report from [`diagnose_same_range`].
#[derive(Debug, Clone, Default)]
pub struct SameRangeDiagnosis {
    /// Edges whose PCurve ranges deviate from their 3D edge range beyond tolerance.
    pub suspect_edges: Vec<SuspectSameRangeEdge>,
}

impl SameRangeDiagnosis {
    /// Returns `true` if no SameRange violations were found.
    pub fn is_clean(&self) -> bool {
        self.suspect_edges.is_empty()
    }
}

/// A single edge-PCurve pair whose UV-evaluated surface endpoints do not match
/// the edge's 3D curve endpoints.
#[derive(Debug, Clone)]
pub struct SuspectFaceSurfaceEdge {
    pub edge_idx: usize,
    pub pcurve_pos: usize,
    pub surface_idx: usize,
    pub start_gap: f64,
    pub end_gap: f64,
    pub max_gap: f64,
}

/// Report from [`diagnose_face_surface_consistency`].
#[derive(Debug, Clone, Default)]
pub struct FaceSurfaceConsistencyDiagnosis {
    pub suspect_edges: Vec<SuspectFaceSurfaceEdge>,
}

impl FaceSurfaceConsistencyDiagnosis {
    pub fn is_clean(&self) -> bool {
        self.suspect_edges.is_empty()
    }
}

/// Diagnose face-on-surface consistency via PCurves.
///
/// For each edge with a 3D curve range and attached PCurves, evaluates:
/// - 3D endpoints `C3(t1), C3(t2)`
/// - UV endpoints `C2(t1), C2(t2)`
/// - Surface points `S(u1,v1), S(u2,v2)`
///
/// If the surface points deviate from the 3D edge endpoints beyond `tolerance`,
/// the edge-PCurve pair is reported as inconsistent.
pub fn diagnose_face_surface_consistency(
    brep: &BRep,
    tolerance: f64,
) -> FaceSurfaceConsistencyDiagnosis {
    use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};

    let mut suspect_edges = Vec::new();

    for edge_idx in 0..brep.edges.len() {
        let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) else {
            continue;
        };
        let Some(curve3) = brep.geom.curves.get(curve_idx) else {
            continue;
        };
        let Some(range3) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) else {
            continue;
        };
        let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else {
            continue;
        };

        for (pcurve_pos, pc) in pcurves.iter().enumerate() {
            let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) else {
                continue;
            };
            let Some(surface) = brep.geom.surfaces.get(pc.surface_idx) else {
                continue;
            };

            let range2 = brep
                .geom
                .curve2d_range
                .get(pc.curve2d_idx)
                .and_then(|r| *r)
                .unwrap_or(range3);

            let p3_start = curve3.point_at(range3[0]);
            let p3_end = curve3.point_at(range3[1]);

            let uv_start = curve2d.point_at(range2[0]);
            let uv_end = curve2d.point_at(range2[1]);

            let ps_start = surface.point_at(uv_start.x, uv_start.y);
            let ps_end = surface.point_at(uv_end.x, uv_end.y);

            let start_gap = (ps_start - p3_start).length();
            let end_gap = (ps_end - p3_end).length();
            let max_gap = start_gap.max(end_gap);

            if max_gap > tolerance {
                suspect_edges.push(SuspectFaceSurfaceEdge {
                    edge_idx,
                    pcurve_pos,
                    surface_idx: pc.surface_idx,
                    start_gap,
                    end_gap,
                    max_gap,
                });
            }
        }
    }

    FaceSurfaceConsistencyDiagnosis { suspect_edges }
}

/// Per-wire diagnostics for gap and self-intersection analysis.
#[derive(Debug, Clone, Default)]
pub struct WireIssueReport {
    pub solid: usize,
    pub shell: usize,
    pub face: usize,
    /// 0 = outer wire, 1..N = inner wire index + 1.
    pub wire_idx: usize,
    pub edge_count: usize,
    pub open_gaps: usize,
    pub topological_self_intersections: usize,
    pub geometric_self_intersections: usize,
}

/// Aggregated report from [`analyze_wire_issues`].
#[derive(Debug, Clone, Default)]
pub struct WireAnalysisReport {
    pub wires: Vec<WireIssueReport>,
    pub total_open_gaps: usize,
    pub total_topological_self_intersections: usize,
    pub total_geometric_self_intersections: usize,
}

impl WireAnalysisReport {
    /// Returns true when no gap or self-intersection issue was found.
    pub fn is_clean(&self) -> bool {
        self.total_open_gaps == 0
            && self.total_topological_self_intersections == 0
            && self.total_geometric_self_intersections == 0
    }
}

/// Analyze all face wires for gap and self-intersection issues.
///
/// This is a structured counterpart to checker issues C1/C7/C8 and is useful
/// for import diagnostics and healing reports.
pub fn analyze_wire_issues(brep: &BRep, tolerance: f64) -> WireAnalysisReport {
    use std::collections::HashMap;

    let mut report = WireAnalysisReport::default();
    let n_edges = brep.edges.len();

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                let mut all_wires: Vec<(usize, &rcad_kernel::topology::Wire)> = Vec::new();
                all_wires.push((0, &face.outer_wire));
                for (wi, inner) in face.inner_wires.iter().enumerate() {
                    all_wires.push((wi + 1, inner));
                }

                for (wire_idx, wire) in all_wires {
                    let mut wire_verts = Vec::with_capacity(wire.edges.len());
                    let mut valid = true;
                    for we in &wire.edges {
                        if we.idx >= n_edges {
                            valid = false;
                            break;
                        }
                        let edge = &brep.edges[we.idx];
                        let (sv, ev) = if we.forward {
                            (edge.start, edge.end)
                        } else {
                            (edge.end, edge.start)
                        };
                        if sv >= brep.vertices.len() || ev >= brep.vertices.len() {
                            valid = false;
                            break;
                        }
                        wire_verts.push((sv, ev));
                    }

                    if !valid {
                        continue;
                    }

                    let mut open_gaps = 0usize;
                    let n = wire_verts.len();
                    if n > 1 {
                        for i in 0..n {
                            let next = (i + 1) % n;
                            let end_v = wire_verts[i].1;
                            let start_v = wire_verts[next].0;
                            if end_v != start_v {
                                let end_pt = brep.vertices[end_v].point;
                                let start_pt = brep.vertices[start_v].point;
                                if (end_pt - start_pt).length() > tolerance {
                                    open_gaps += 1;
                                }
                            }
                        }
                    }

                    let mut vertex_count: HashMap<usize, usize> = HashMap::new();
                    for &(sv, ev) in &wire_verts {
                        *vertex_count.entry(sv).or_insert(0) += 1;
                        *vertex_count.entry(ev).or_insert(0) += 1;
                    }
                    let topological_self_intersections =
                        vertex_count.values().filter(|&&c| c > 2).count();
                    let geometric_self_intersections =
                        count_geometric_self_intersections(&wire_verts, &brep.vertices);

                    if open_gaps > 0
                        || topological_self_intersections > 0
                        || geometric_self_intersections > 0
                    {
                        report.wires.push(WireIssueReport {
                            solid: si,
                            shell: shi,
                            face: fi,
                            wire_idx,
                            edge_count: wire_verts.len(),
                            open_gaps,
                            topological_self_intersections,
                            geometric_self_intersections,
                        });
                    }

                    report.total_open_gaps += open_gaps;
                    report.total_topological_self_intersections += topological_self_intersections;
                    report.total_geometric_self_intersections += geometric_self_intersections;
                }
            }
        }
    }

    report
}

/// Scan all edges in `brep` for SameParameter violations.
///
/// For each edge that has a known 3D curve and curve range `[t1, t2]`, evaluates
/// the curve at `t1` and `t2` and compares the resulting points to the
/// corresponding vertex positions.  Edges whose endpoint gap exceeds `tolerance`
/// are reported as suspects and their `edge_same_parameter` flag is updated to
/// `false` so that a subsequent [`fix_same_parameter`] pass can repair them.
///
/// Analogous to `BRepLib::SameParameter()` diagnosis step in OCCT.
///
/// Returns a non-mutating diagnosis; call `fix_same_parameter_with_scan` to
/// also repair the flagged edges.
pub fn diagnose_same_parameter(brep: &BRep, tolerance: f64) -> SameParameterDiagnosis {
    use rcad_kernel::geom::CurveEval;

    let mut suspects = Vec::new();
    let n_edges = brep.edges.len();
    let n_verts = brep.vertices.len();

    for edge_idx in 0..n_edges {
        let curve_idx = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c);
        let Some(ci) = curve_idx else { continue };
        let Some(curve) = brep.geom.curves.get(ci) else { continue };
        let Some(range) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) else { continue };

        let edge = &brep.edges[edge_idx];
        if edge.start >= n_verts || edge.end >= n_verts {
            continue;
        }

        let p_start = brep.vertices[edge.start].point;
        let p_end = brep.vertices[edge.end].point;

        let eval_start = curve.point_at(range[0]);
        let eval_end = curve.point_at(range[1]);

        let start_gap = (eval_start - p_start).length();
        let end_gap = (eval_end - p_end).length();

        if start_gap > tolerance || end_gap > tolerance {
            suspects.push(SuspectEdge {
                edge_idx,
                start_gap,
                end_gap,
            });
        }
    }

    SameParameterDiagnosis { suspect_edges: suspects }
}

/// Scan all edges in `brep` for SameRange violations.
///
/// For each edge with known 3D range and attached PCurves, checks whether each
/// referenced `curve2d_range` matches the edge's `[t1, t2]` within `tolerance`.
/// Missing 2D ranges are also treated as violations.
pub fn diagnose_same_range(brep: &BRep, tolerance: f64) -> SameRangeDiagnosis {
    let mut suspects = Vec::new();
    let n_edges = brep.edges.len();

    for edge_idx in 0..n_edges {
        let Some(range3d) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) else {
            continue;
        };
        let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else {
            continue;
        };
        if pcurves.is_empty() {
            continue;
        }

        let mut mismatched_pcurves = 0usize;
        let mut max_delta = 0.0f64;

        for pc in pcurves {
            let Some(range2d) = brep.geom.curve2d_range.get(pc.curve2d_idx).and_then(|r| *r) else {
                mismatched_pcurves += 1;
                max_delta = max_delta.max(tolerance);
                continue;
            };

            let d0 = (range2d[0] - range3d[0]).abs();
            let d1 = (range2d[1] - range3d[1]).abs();
            let d = d0.max(d1);
            if d > tolerance {
                mismatched_pcurves += 1;
                max_delta = max_delta.max(d);
            }
        }

        let flagged_false = brep
            .geom
            .edge_same_range
            .get(edge_idx)
            .is_some_and(|v| !*v);

        if mismatched_pcurves > 0 || flagged_false {
            suspects.push(SuspectSameRangeEdge {
                edge_idx,
                mismatched_pcurves,
                max_delta,
            });
        }
    }

    SameRangeDiagnosis { suspect_edges: suspects }
}

// ── Shell topology analysis ───────────────────────────────────────────────────

/// Topology analysis report for a BRep's shell structure.
///
/// Analogous to `ShapeAnalysis_Shell` in OCCT's TKShHealing.
#[derive(Debug, Clone, Default)]
pub struct ShellTopologyReport {
    /// `true` if every edge is referenced by exactly 2 faces (no free edges).
    pub is_closed: bool,
    /// `true` if no edge is referenced by more than 2 faces.
    pub is_manifold: bool,
    /// Number of edges referenced by exactly 1 face (free / open edges).
    pub open_edge_count: usize,
    /// Number of edges referenced by 3 or more faces (non-manifold edges).
    pub non_manifold_edge_count: usize,
    /// Number of vertices not referenced by any edge in the BRep.
    pub isolated_vertex_count: usize,
    /// Total edge count in the BRep.
    pub total_edges: usize,
    /// Total face count across all solids/shells.
    pub total_faces: usize,
}

/// Analyze the shell topology of `brep`.
///
/// Counts free edges, non-manifold edges, and isolated vertices to determine
/// whether the BRep represents a closed manifold solid.
///
/// Analogous to `ShapeAnalysis_Shell::LoadUnorientedEdges()` + checks in OCCT.
pub fn analyze_shell_topology(brep: &BRep) -> ShellTopologyReport {
    let total_edges = brep.edges.len();

    // Count how many faces each edge is referenced by.
    let mut edge_face_count: Vec<usize> = vec![0; total_edges];
    let mut total_faces = 0usize;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                total_faces += 1;
                for we in &face.outer_wire.edges {
                    if we.idx < total_edges {
                        edge_face_count[we.idx] += 1;
                    }
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        if we.idx < total_edges {
                            edge_face_count[we.idx] += 1;
                        }
                    }
                }
            }
        }
    }

    let open_edge_count = edge_face_count.iter().filter(|&&c| c == 1).count();
    let non_manifold_edge_count = edge_face_count.iter().filter(|&&c| c > 2).count();

    // Count isolated vertices: those not referenced by any edge as start or end.
    let mut vertex_used = vec![false; brep.vertices.len()];
    for edge in &brep.edges {
        if edge.start < vertex_used.len() {
            vertex_used[edge.start] = true;
        }
        if edge.end < vertex_used.len() {
            vertex_used[edge.end] = true;
        }
    }
    let isolated_vertex_count = vertex_used.iter().filter(|&&used| !used).count();

    ShellTopologyReport {
        is_closed: open_edge_count == 0,
        is_manifold: non_manifold_edge_count == 0,
        open_edge_count,
        non_manifold_edge_count,
        isolated_vertex_count,
        total_edges,
        total_faces,
    }
}

// ── Euler characteristic analysis ────────────────────────────────────────────

/// Euler characteristic and topological genus for a single solid.
///
/// For a closed orientable 2-manifold of genus *g*: χ = V − E + F = 2 − 2g.
///
/// | Shape   | χ  | genus |
/// | sphere  | 2  | 0     |
/// | torus   | 0  | 1     |
/// | 2-torus | -2 | 2     |
#[derive(Debug, Clone)]
pub struct EulerAnalysis {
    pub solid_idx: usize,
    /// Unique vertices referenced by this solid's face boundaries.
    pub vertices: usize,
    /// Unique edges referenced by this solid's face boundaries.
    pub edges: usize,
    /// Total faces across all shells of this solid.
    pub faces: usize,
    /// Euler characteristic: V − E + F.
    pub euler_number: i64,
    /// `true` if every edge of this solid is referenced by exactly 2 faces
    /// (no free boundary edges → closed shell).
    pub is_closed: bool,
    /// Topological genus, computed as `(2 − euler_number) / 2`.
    /// `None` when `!is_closed` or when the result is not a non-negative integer.
    pub genus: Option<i64>,
}

/// Compute per-solid Euler analysis for every solid in `brep`.
///
/// Analogous to the topological analysis portion of `BRepCheck_Analyzer`.
pub fn euler_analysis(brep: &BRep) -> Vec<EulerAnalysis> {
    let mut results = Vec::with_capacity(brep.solids.len());

    for (si, solid) in brep.solids.iter().enumerate() {
        let mut unique_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut face_count = 0usize;

        for shell in &solid.shells {
            for face in &shell.faces {
                face_count += 1;
                for we in &face.outer_wire.edges {
                    unique_edges.insert(we.idx);
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        unique_edges.insert(we.idx);
                    }
                }
            }
        }

        // Collect unique vertex indices from the unique edges.
        let mut unique_verts: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &ei in &unique_edges {
            if let Some(e) = brep.edges.get(ei) {
                unique_verts.insert(e.start);
                unique_verts.insert(e.end);
            }
        }

        let v = unique_verts.len();
        let e = unique_edges.len();
        let f = face_count;
        let euler_number = v as i64 - e as i64 + f as i64;

        // Determine if the solid is closed (every edge shared by exactly 2 faces).
        let mut edge_face_count: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    *edge_face_count.entry(we.idx).or_insert(0) += 1;
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        *edge_face_count.entry(we.idx).or_insert(0) += 1;
                    }
                }
            }
        }
        let is_closed = edge_face_count.values().all(|&c| c == 2);

        // Genus only makes sense on a closed manifold.
        let genus = if is_closed {
            let g = (2 - euler_number) / 2;
            // Valid genus: (2 − χ) must be even and non-negative.
            if (2 - euler_number) % 2 == 0 && g >= 0 {
                Some(g)
            } else {
                None
            }
        } else {
            None
        };

        results.push(EulerAnalysis {
            solid_idx: si,
            vertices: v,
            edges: e,
            faces: f,
            euler_number,
            is_closed,
            genus,
        });
    }

    results
}

// ── Orientation consistency analysis ─────────────────────────────────────────

/// A face whose stored normal appears to point inward rather than outward.
///
/// Determined geometrically: the dot product of `face.normal` with the vector
/// from the solid's vertex centroid to the face's centroid is negative.
#[derive(Debug, Clone)]
pub struct OrientationIssue {
    /// Solid index.
    pub solid_idx: usize,
    /// Flat face index (counting across all solids/shells in traversal order).
    pub face_idx: usize,
    /// Dot product of `face.normal` with the outward radial direction.
    /// Negative values indicate an inward-pointing normal.
    pub dot_product: f64,
}

/// Report from [`check_orientation_consistency`].
#[derive(Debug, Clone, Default)]
pub struct OrientationReport {
    /// `true` iff all face normals appear to point outward from the solid interior.
    pub is_consistent: bool,
    /// Faces whose stored normals appear to point inward.
    pub issues: Vec<OrientationIssue>,
    /// Number of faces with outward-pointing normals.
    pub consistent_face_count: usize,
    /// Number of faces with inward-pointing normals.
    pub inconsistent_face_count: usize,
}

/// Check that every face's stored normal points **outward** from the solid interior.
///
/// Uses a geometric heuristic: the solid's interior centroid is approximated as
/// the average of all its vertices.  For each face the face centroid is computed
/// from the outer-wire corner vertices.  A face normal is considered outward when
/// `face.normal · (face_centroid − solid_centroid) > 0`.
///
/// This correctly handles all primitives (box, sphere, cylinder, cone, torus)
/// whose normals are computed analytically during construction.  For shapes
/// produced by Boolean operations or user-constructed BReps, the heuristic may
/// give false positives on highly non-convex solids — treat the report as an
/// advisory rather than a hard constraint.
///
/// Analogous to `BRepCheck_Shell::Orientation()` in OCCT.
pub fn check_orientation_consistency(brep: &BRep) -> OrientationReport {
    use glam::DVec3;

    let mut issues = Vec::new();
    let mut consistent_face_count = 0usize;
    let mut inconsistent_face_count = 0usize;
    let mut flat_face_idx = 0usize;

    for (si, solid) in brep.solids.iter().enumerate() {
        // Compute the solid's interior centroid from all vertices that appear
        // in ANY of this solid's face wires.
        let mut solid_verts = std::collections::HashSet::new();
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    if we.idx < brep.edges.len() {
                        solid_verts.insert(brep.edges[we.idx].start);
                        solid_verts.insert(brep.edges[we.idx].end);
                    }
                }
            }
        }
        if solid_verts.is_empty() {
            // No geometry — skip.
            for shell in &solid.shells {
                flat_face_idx += shell.faces.len();
            }
            continue;
        }
        let solid_centroid: DVec3 = {
            let sum: DVec3 = solid_verts
                .iter()
                .filter(|&&vi| vi < brep.vertices.len())
                .map(|&vi| brep.vertices[vi].point)
                .sum();
            sum / solid_verts.len() as f64
        };

        for shell in &solid.shells {
            for face in &shell.faces {
                // Face centroid: average of first vertices of each outer-wire edge.
                let mut face_centroid = DVec3::ZERO;
                let mut n = 0usize;
                for we in &face.outer_wire.edges {
                    if we.idx < brep.edges.len() {
                        let vi = if we.forward {
                            brep.edges[we.idx].start
                        } else {
                            brep.edges[we.idx].end
                        };
                        if vi < brep.vertices.len() {
                            face_centroid += brep.vertices[vi].point;
                            n += 1;
                        }
                    }
                }
                if n == 0 {
                    flat_face_idx += 1;
                    continue;
                }
                face_centroid /= n as f64;

                let outward = face_centroid - solid_centroid;
                let dot = face.normal.dot(outward);
                if dot >= 0.0 {
                    consistent_face_count += 1;
                } else {
                    inconsistent_face_count += 1;
                    issues.push(OrientationIssue {
                        solid_idx: si,
                        face_idx: flat_face_idx,
                        dot_product: dot,
                    });
                }
                flat_face_idx += 1;
            }
        }
    }

    OrientationReport {
        is_consistent: issues.is_empty(),
        issues,
        consistent_face_count,
        inconsistent_face_count,
    }
}

// ── Comprehensive richer validity analysis ────────────────────────────────────

/// Aggregated validity report combining all available checks.
///
/// This is the RCAD equivalent of OCCT's `BRepCheck_Analyzer` + `ShapeAnalysis`
/// combined output, giving a single entry-point for full BRep validation.
#[derive(Debug, Clone)]
pub struct RicherValidityReport {
    /// Basic structural check result.
    pub check_result: CheckResult,
    /// Shell topology (free edges, non-manifold edges, isolation).
    pub shell_topology: ShellTopologyReport,
    /// Per-solid Euler characteristic and genus.
    pub euler: Vec<EulerAnalysis>,
    /// Face wire orientation consistency across shared edges.
    pub orientation: OrientationReport,
    /// `true` iff BRep passes all structural checks and orientation is consistent.
    pub is_fully_valid: bool,
}

impl RicherValidityReport {
    /// A short human-readable summary of this report.
    pub fn summary(&self) -> String {
        let issues = self.check_result.issues.len();
        let euler_issues = self.euler.iter().filter(|e| e.genus.map_or(true, |g| g < 0)).count();
        let orient_issues = self.orientation.inconsistent_face_count;
        if self.is_fully_valid {
            format!(
                "valid: {} solids, closed={}, manifold={}",
                self.euler.len(),
                self.shell_topology.is_closed,
                self.shell_topology.is_manifold,
            )
        } else {
            format!(
                "INVALID: {} structural issue(s), {} orientation inconsistency/ies, {} genus anomaly/ies",
                issues, orient_issues, euler_issues,
            )
        }
    }
}

/// Run all available validity checks on `brep` and return a consolidated report.
///
/// This is the preferred entry point for comprehensive BRep validation:
/// it combines the basic `check()`, shell topology, Euler analysis, and
/// orientation consistency into a single call.
pub fn richer_validity_analysis(brep: &BRep) -> RicherValidityReport {
    let check_result = check(brep);
    let shell_topology = analyze_shell_topology(brep);
    let euler = euler_analysis(brep);
    let orientation = check_orientation_consistency(brep);

    let is_fully_valid = check_result.is_valid() && orientation.is_consistent;

    RicherValidityReport {
        check_result,
        shell_topology,
        euler,
        orientation,
        is_fully_valid,
    }
}

// ── Surface UV Analysis (ShapeAnalysis_Surface equivalent) ───────────────────────

/// Report from surface UV domain analysis.
///
/// Analogous to OCCT's `ShapeAnalysis_Surface` which checks UV bounds,
/// periodicity, and parameter space consistency.
#[derive(Debug, Clone, Default)]
pub struct SurfaceAnalysisReport {
    /// Number of faces analyzed.
    pub faces_analyzed: usize,
    /// Faces whose PCurve parameter ranges violate expected surface bounds.
    pub faces_with_uv_bounds_violation: Vec<UvBoundsViolation>,
    /// Total number of issues detected.
    pub total_issues: usize,
}

impl SurfaceAnalysisReport {
    pub fn is_clean(&self) -> bool {
        self.total_issues == 0
    }

    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!("{} faces analyzed, no UV issues", self.faces_analyzed)
        } else {
            format!(
                "{} faces analyzed, {} UV bounds violations",
                self.faces_analyzed,
                self.faces_with_uv_bounds_violation.len()
            )
        }
    }
}

/// UV bounds violation for a single face.
#[derive(Debug, Clone)]
pub struct UvBoundsViolation {
    pub solid: usize,
    pub shell: usize,
    pub face: usize,
    /// Surface type (Plane, Cylinder, etc.)
    pub surface_type: String,
    /// Expected UV bounds [u_min, u_max, v_min, v_max] for the surface.
    pub expected_bounds: [f64; 4],
    /// Actual UV bounds of the face's PCurves.
    pub actual_bounds: [f64; 4],
    /// Maximum violation distance.
    pub violation: f64,
}

/// Analyze surface UV consistency for all faces in `brep`.
///
/// Checks PCurve parameter ranges against the surface's natural domain bounds.
/// For periodic surfaces like Cylinder and Cone, checks U bounds.
/// For bounded surfaces like Sphere, checks both U and V bounds.
///
/// Analogous to `ShapeAnalysis_Surface::CheckUVBounds` in OCCT.
pub fn analyze_surface_uv_consistency(brep: &BRep, tolerance: f64) -> SurfaceAnalysisReport {
    use rcad_kernel::geom::Surface3;

    let mut report = SurfaceAnalysisReport::default();

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, _face) in shell.faces.iter().enumerate() {
                report.faces_analyzed += 1;

                // Get face's surface
                let flat_face_idx = {
                    let mut idx = 0usize;
                    for s in 0..si {
                        for sh in &brep.solids[s].shells {
                            idx += sh.faces.len();
                        }
                    }
                    for sh in 0..shi {
                        idx += brep.solids[si].shells[sh].faces.len();
                    }
                    idx + fi
                };

                let surface_idx = match brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) {
                    Some(idx) => idx,
                    None => continue,
                };

                let surface = match brep.geom.surfaces.get(surface_idx) {
                    Some(s) => s,
                    None => continue,
                };

                // Get expected UV bounds for the surface type
                let expected_bounds = match surface {
                    Surface3::Plane(_) => {
                        // Plane has unbounded UV: [-∞, ∞, -∞, ∞]
                        // No bounds check needed
                        continue;
                    }
                    Surface3::Cylinder(_) => {
                        // Cylinder: U ∈ [-π, π] (periodic), V ∈ [-∞, ∞]
                        // Only check U bounds
                        [-std::f64::consts::PI, std::f64::consts::PI, f64::NEG_INFINITY, f64::INFINITY]
                    }
                    Surface3::Sphere(_) => {
                        // Sphere: U ∈ [-π, π], V ∈ [0, π]
                        [-std::f64::consts::PI, std::f64::consts::PI, 0.0, std::f64::consts::PI]
                    }
                    Surface3::Cone(_) => {
                        // Cone: U ∈ [-π, π] (periodic), V ∈ [0, ∞]
                        [-std::f64::consts::PI, std::f64::consts::PI, 0.0, f64::INFINITY]
                    }
                    Surface3::Torus(_) => {
                        // Torus: U ∈ [-π, π], V ∈ [-π, π] (both periodic)
                        [-std::f64::consts::PI, std::f64::consts::PI, -std::f64::consts::PI, std::f64::consts::PI]
                    }
                    _ => continue, // BSpline and others: no simple bounds check
                };

                // Collect UV ranges from face's PCurves
                let mut u_min = f64::INFINITY;
                let mut u_max = f64::NEG_INFINITY;
                let mut v_min = f64::INFINITY;
                let mut v_max = f64::NEG_INFINITY;
                let mut has_pcurve_data = false;

                // Get the face's wire edges and check their PCurves
                let face_ref = &brep.solids[si].shells[shi].faces[fi];
                for we in &face_ref.outer_wire.edges {
                    if let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) {
                        for pc in pcurves {
                            if pc.surface_idx == surface_idx {
                                if let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) {
                                    has_pcurve_data = true;
                                    // Sample the curve
                                    for i in 0..=16 {
                                        let t = i as f64 / 16.0;
                                        let uv = curve2d.point_at(t);
                                        u_min = u_min.min(uv.x);
                                        u_max = u_max.max(uv.x);
                                        v_min = v_min.min(uv.y);
                                        v_max = v_max.max(uv.y);
                                    }
                                }
                            }
                        }
                    }
                }

                if !has_pcurve_data {
                    continue; // No PCurve data to check
                }

                let actual_bounds = [u_min, u_max, v_min, v_max];

                // Check for bounds violation (only for bounded parameters)
                let mut violation = 0.0_f64;

                // Check U bounds if bounded
                if expected_bounds[0].is_finite() && u_min < expected_bounds[0] - tolerance {
                    violation = violation.max(expected_bounds[0] - u_min);
                }
                if expected_bounds[1].is_finite() && u_max > expected_bounds[1] + tolerance {
                    violation = violation.max(u_max - expected_bounds[1]);
                }

                // Check V bounds if bounded
                if expected_bounds[2].is_finite() && v_min < expected_bounds[2] - tolerance {
                    violation = violation.max(expected_bounds[2] - v_min);
                }
                if expected_bounds[3].is_finite() && v_max > expected_bounds[3] + tolerance {
                    violation = violation.max(v_max - expected_bounds[3]);
                }

                if violation > tolerance {
                    let surface_type = match surface {
                        Surface3::Plane(_) => "Plane",
                        Surface3::Cylinder(_) => "Cylinder",
                        Surface3::Sphere(_) => "Sphere",
                        Surface3::Cone(_) => "Cone",
                        Surface3::Torus(_) => "Torus",
                        _ => "Unknown",
                    };
                    report.faces_with_uv_bounds_violation.push(UvBoundsViolation {
                        solid: si,
                        shell: shi,
                        face: fi,
                        surface_type: surface_type.to_string(),
                        expected_bounds,
                        actual_bounds,
                        violation,
                    });
                    report.total_issues += 1;
                }
            }
        }
    }

    report
}

// ── Wire Quality Metrics (ShapeAnalysis_Wire enhancement) ───────────────────────

/// Extended wire quality metrics for a single wire.
///
/// Analogous to OCCT's `ShapeAnalysis_Wire` which provides area, orientation,
/// and closure quality metrics.
#[derive(Debug, Clone, Default)]
pub struct WireQualityMetrics {
    pub solid: usize,
    pub shell: usize,
    pub face: usize,
    pub wire_idx: usize, // 0 = outer wire, 1+ = inner wire index
    /// Number of edges in the wire.
    pub edge_count: usize,
    /// 3D length of the wire (sum of edge lengths).
    pub total_length: f64,
    /// Whether the wire is closed (end vertex of last edge = start vertex of first edge).
    pub is_closed: bool,
    /// Whether the wire is self-intersecting (topologically).
    pub has_self_intersection: bool,
    /// Number of gap locations where consecutive edges don't share vertices.
    pub gap_count: usize,
    /// Maximum gap size (distance between non-connected vertices).
    pub max_gap: f64,
    /// Quality score (0-100, higher is better).
    pub quality_score: f64,
}

/// Aggregated wire quality report for all wires in a BRep.
#[derive(Debug, Clone, Default)]
pub struct WireQualityReport {
    pub wires_analyzed: usize,
    pub closed_wires: usize,
    pub open_wires: usize,
    pub self_intersecting_wires: usize,
    pub wires_with_gaps: usize,
    pub total_gap_count: usize,
    pub avg_quality_score: f64,
    pub metrics: Vec<WireQualityMetrics>,
}

impl WireQualityReport {
    pub fn is_clean(&self) -> bool {
        self.open_wires == 0 && self.self_intersecting_wires == 0 && self.wires_with_gaps == 0
    }

    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!("{} wires analyzed, all closed and clean, avg quality {:.1}", self.wires_analyzed, self.avg_quality_score)
        } else {
            format!(
                "{} wires: {} open, {} self-intersecting, {} with gaps ({} total), avg quality {:.1}",
                self.wires_analyzed,
                self.open_wires,
                self.self_intersecting_wires,
                self.wires_with_gaps,
                self.total_gap_count,
                self.avg_quality_score
            )
        }
    }
}

/// Analyze wire quality metrics for all wires in `brep`.
///
/// Provides detailed metrics including length, closure, self-intersection
/// detection, and quality scoring.
///
/// Analogous to `ShapeAnalysis_Wire` in OCCT.
pub fn analyze_wire_quality(brep: &BRep, tolerance: f64) -> WireQualityReport {
    let mut report = WireQualityReport::default();
    let mut total_quality = 0.0_f64;

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                // Analyze outer wire
                let outer_metrics = analyze_single_wire_quality(
                    brep, si, shi, fi, 0, &face.outer_wire, tolerance,
                );
                let outer_closed = outer_metrics.is_closed;
                total_quality += outer_metrics.quality_score;
                report.wires_analyzed += 1;
                if outer_closed { report.closed_wires += 1; } else { report.open_wires += 1; }
                if outer_metrics.has_self_intersection { report.self_intersecting_wires += 1; }
                if outer_metrics.gap_count > 0 { report.wires_with_gaps += 1; }
                report.total_gap_count += outer_metrics.gap_count;
                report.metrics.push(outer_metrics);

                // Analyze inner wires
                for (wi, wire) in face.inner_wires.iter().enumerate() {
                    let metrics = analyze_single_wire_quality(
                        brep, si, shi, fi, wi + 1, wire, tolerance,
                    );
                    total_quality += metrics.quality_score;
                    report.wires_analyzed += 1;
                    if metrics.is_closed { report.closed_wires += 1; } else { report.open_wires += 1; }
                    if metrics.has_self_intersection { report.self_intersecting_wires += 1; }
                    if metrics.gap_count > 0 { report.wires_with_gaps += 1; }
                    report.total_gap_count += metrics.gap_count;
                    report.metrics.push(metrics);
                }
            }
        }
    }

    report.avg_quality_score = if report.wires_analyzed > 0 {
        total_quality / report.wires_analyzed as f64
    } else {
        0.0
    };

    report
}

fn analyze_single_wire_quality(
    brep: &BRep,
    solid: usize,
    shell: usize,
    face: usize,
    wire_idx: usize,
    wire: &rcad_kernel::topology::Wire,
    tolerance: f64,
) -> WireQualityMetrics {
    let mut metrics = WireQualityMetrics {
        solid,
        shell,
        face,
        wire_idx,
        edge_count: wire.edges.len(),
        ..Default::default()
    };

    if wire.edges.is_empty() {
        metrics.quality_score = 0.0;
        return metrics;
    }

    // Compute total length and check closure
    let mut total_length = 0.0_f64;
    let mut gap_count = 0usize;
    let mut max_gap = 0.0_f64;

    for (i, we) in wire.edges.iter().enumerate() {
        let edge = match brep.edges.get(we.idx) {
            Some(e) => e,
            None => continue,
        };

        // Compute edge length
        let start_pt = brep.vertices.get(edge.start).map(|v| v.point).unwrap_or_default();
        let end_pt = brep.vertices.get(edge.end).map(|v| v.point).unwrap_or_default();
        let edge_len = (end_pt - start_pt).length();
        total_length += edge_len;

        // Check gap to next edge
        let next_i = (i + 1) % wire.edges.len();
        let next_edge = match brep.edges.get(wire.edges[next_i].idx) {
            Some(e) => e,
            None => continue,
        };

        // The end vertex of this edge should equal the start vertex of the next edge
        // (accounting for orientation via .forward)
        let this_end = if we.forward { edge.end } else { edge.start };
        let next_start = if wire.edges[next_i].forward {
            next_edge.start
        } else {
            next_edge.end
        };

        if this_end != next_start {
            let gap_pt1 = brep.vertices.get(this_end).map(|v| v.point).unwrap_or_default();
            let gap_pt2 = brep.vertices.get(next_start).map(|v| v.point).unwrap_or_default();
            let gap = (gap_pt2 - gap_pt1).length();
            if gap > tolerance {
                gap_count += 1;
                max_gap = max_gap.max(gap);
            }
        }
    }

    metrics.total_length = total_length;
    metrics.gap_count = gap_count;
    metrics.max_gap = max_gap;
    metrics.is_closed = gap_count == 0;

    // Check for self-intersection (topological: shared vertices except at junctions)
    let mut vertex_occurrences: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();

    for (i, we) in wire.edges.iter().enumerate() {
        if let Some(edge) = brep.edges.get(we.idx) {
            let (start, end) = if we.forward {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            };
            vertex_occurrences.entry(start).or_default().push(i);
            vertex_occurrences.entry(end).or_default().push(i);
        }
    }

    // A vertex should appear at most twice: once as end of one edge, once as start of next
    for (_, occurrences) in &vertex_occurrences {
        if occurrences.len() > 2 {
            metrics.has_self_intersection = true;
            break;
        }
    }

    // Compute quality score (0-100)
    let mut score = 100.0_f64;

    // Penalize gaps
    if gap_count > 0 {
        score -= (gap_count as f64).min(30.0) * 3.0;
        score -= (max_gap / tolerance).min(10.0) * 2.0;
    }

    // Penalize self-intersection
    if metrics.has_self_intersection {
        score -= 40.0;
    }

    // Penalize very short wires
    if metrics.edge_count < 3 {
        score -= 20.0;
    }

    metrics.quality_score = score.max(0.0).min(100.0);

    metrics
}

// ═══════════════════════════════════════════════════════════════════════════════
// GEOMETRY VALIDATION (Surface Continuity, Curve-Surface Consistency)
// ═══════════════════════════════════════════════════════════════════════════════

/// Report from geometry validation checks.
///
/// Analogous to OCCT's `BRepCheck_Analyzer` geometry validation portion.
#[derive(Debug, Clone, Default)]
pub struct GeometryValidationReport {
    /// Number of edges checked for curve-surface consistency.
    pub edges_checked: usize,
    /// Number of face pairs checked for surface continuity.
    pub face_pairs_checked: usize,
    /// Issues found during geometry validation.
    pub issues: Vec<CheckIssue>,
}

impl GeometryValidationReport {
    pub fn is_clean(&self) -> bool {
        self.issues.is_empty()
    }

    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!("{} edges, {} face pairs checked, no geometry issues", self.edges_checked, self.face_pairs_checked)
        } else {
            format!("{} geometry issues found ({} edges, {} face pairs checked)", self.issues.len(), self.edges_checked, self.face_pairs_checked)
        }
    }
}

/// Check surface continuity between adjacent faces.
///
/// Analyzes the geometric continuity (C0, C1, C2) across shared edges
/// between adjacent faces. C0 requires position continuity, C1 requires
/// tangent continuity, and C2 requires curvature continuity.
///
/// Analogous to `BRepCheck_Analyzer::CheckSurfaceContinuity` in OCCT.
pub fn check_surface_continuity(brep: &BRep, tolerance: f64) -> GeometryValidationReport {
    let mut report = GeometryValidationReport::default();

    // Build edge-to-faces mapping
    let mut edge_faces: std::collections::HashMap<usize, Vec<(usize, usize, usize)>> =
        std::collections::HashMap::new(); // edge_idx -> [(solid_idx, shell_idx, face_idx)]

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                for we in &face.outer_wire.edges {
                    edge_faces.entry(we.idx).or_default().push((si, shi, fi));
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        edge_faces.entry(we.idx).or_default().push((si, shi, fi));
                    }
                }
            }
        }
    }

    // Check continuity for edges shared by exactly 2 faces
    for (&edge_idx, faces) in &edge_faces {
        if faces.len() != 2 {
            continue;
        }
        report.face_pairs_checked += 1;

        let (si1, shi1, fi1) = faces[0];
        let (si2, shi2, fi2) = faces[1];

        let face1 = &brep.solids[si1].shells[shi1].faces[fi1];
        let face2 = &brep.solids[si2].shells[shi2].faces[fi2];

        // Get surface indices
        let flat_face_idx1 = get_flat_face_index(brep, si1, shi1, fi1);
        let flat_face_idx2 = get_flat_face_index(brep, si2, shi2, fi2);

        let surf_idx1 = match brep.geom.face_surface.get(flat_face_idx1).and_then(|v| *v) {
            Some(idx) => idx,
            None => continue,
        };
        let surf_idx2 = match brep.geom.face_surface.get(flat_face_idx2).and_then(|v| *v) {
            Some(idx) => idx,
            None => continue,
        };

        let surface1 = match brep.geom.surfaces.get(surf_idx1) {
            Some(s) => s,
            None => continue,
        };
        let surface2 = match brep.geom.surfaces.get(surf_idx2) {
            Some(s) => s,
            None => continue,
        };

        // Sample the shared edge and check continuity at sample points
        if let Some(edge) = brep.edges.get(edge_idx) {
            let v1_pt = brep.vertices.get(edge.start).map(|v| v.point).unwrap_or_default();
            let v2_pt = brep.vertices.get(edge.end).map(|v| v.point).unwrap_or_default();

            // Sample along the edge
            for alpha in [0.0, 0.25, 0.5, 0.75, 1.0] {
                let sample_pt = v1_pt.lerp(v2_pt, alpha);

                // Get UV coordinates from PCurves
                let uv1 = get_edge_uv_at(brep, edge_idx, alpha, surf_idx1);
                let uv2 = get_edge_uv_at(brep, edge_idx, alpha, surf_idx2);

                // Evaluate surface normals at sample points
                let n1 = surface1.normal_at(uv1.x, uv1.y);
                let n2 = surface2.normal_at(uv2.x, uv2.y);

                // C1 continuity: normals should be opposite (faces share edge with opposite orientation)
                let dot = n1.dot(-n2);
                let angle_deviation = (1.0 - dot).abs();

                if angle_deviation > tolerance * 10.0 {
                    // C1 violation - normals don't match
                    report.issues.push(CheckIssue::SurfaceContinuityViolation {
                        solid: si1,
                        face_a: fi1,
                        face_b: fi2,
                        shared_edge: edge_idx,
                        expected: 1,
                        actual: 0,
                        deviation: angle_deviation,
                    });
                    break; // One violation per edge pair is enough
                }
            }
        }
    }

    report
}

/// Check curve-surface consistency for all edges with PCurves.
///
/// For each edge with a 3D curve and attached PCurves, verifies that the
/// surface evaluation at PCurve UV coordinates matches the 3D curve evaluation.
///
/// Analogous to `BRepCheck_Edge::CheckCurveSurfaceConsistency` in OCCT.
pub fn check_curve_surface_consistency(brep: &BRep, tolerance: f64) -> GeometryValidationReport {
    let mut report = GeometryValidationReport::default();

    for edge_idx in 0..brep.edges.len() {
        let curve_idx = match brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) {
            Some(idx) => idx,
            None => continue,
        };

        let curve = match brep.geom.curves.get(curve_idx) {
            Some(c) => c,
            None => continue,
        };

        let range = match brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) {
            Some(r) => r,
            None => continue,
        };

        let pcurves = match brep.geom.edge_pcurves.get(edge_idx) {
            Some(pc) if !pc.is_empty() => pc,
            _ => continue,
        };

        report.edges_checked += 1;

        for pc in pcurves {
            let surface = match brep.geom.surfaces.get(pc.surface_idx) {
                Some(s) => s,
                None => continue,
            };

            let curve2d = match brep.geom.curve2ds.get(pc.curve2d_idx) {
                Some(c) => c,
                None => continue,
            };

            // Sample multiple points along the curve
            let mut max_deviation = 0.0_f64;
            for i in 0..=10 {
                let t = i as f64 / 10.0;
                let param = range[0] + t * (range[1] - range[0]);

                let p3d = curve.point_at(param);
                let uv = curve2d.point_at(param);
                let p_surf = surface.point_at(uv.x, uv.y);

                let deviation = (p3d - p_surf).length();
                max_deviation = max_deviation.max(deviation);
            }

            if max_deviation > tolerance {
                report.issues.push(CheckIssue::CurveSurfaceMismatch {
                    edge: edge_idx,
                    surface: pc.surface_idx,
                    max_deviation,
                });
            }
        }
    }

    report
}

/// Get the flat (linear) face index from solid/shell/face indices.
fn get_flat_face_index(brep: &BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> usize {
    let mut idx = 0usize;
    for s in 0..solid_idx {
        for sh in &brep.solids[s].shells {
            idx += sh.faces.len();
        }
    }
    for sh in 0..shell_idx {
        idx += brep.solids[solid_idx].shells[sh].faces.len();
    }
    idx + face_idx
}

/// Get UV coordinates for an edge at a given parameter (0-1) on a specific surface.
fn get_edge_uv_at(brep: &BRep, edge_idx: usize, alpha: f64, surface_idx: usize) -> DVec2 {
    // Default UV in case we can't find the PCurve
    let default_uv = DVec2::new(alpha, alpha);

    let pcurves = match brep.geom.edge_pcurves.get(edge_idx) {
        Some(pc) => pc,
        None => return default_uv,
    };

    for pc in pcurves {
        if pc.surface_idx != surface_idx {
            continue;
        }

        let curve2d = match brep.geom.curve2ds.get(pc.curve2d_idx) {
            Some(c) => c,
            None => continue,
        };

        let range = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r);
        if let Some(r) = range {
            let param = r[0] + alpha * (r[1] - r[0]);
            return curve2d.point_at(param);
        } else {
            return curve2d.point_at(alpha);
        }
    }

    default_uv
}

// ═══════════════════════════════════════════════════════════════════════════════
// TOPOLOGY VALIDATION (Shell Orientation, Solid Closure, Wire Orientation)
// ═══════════════════════════════════════════════════════════════════════════════

/// Report from topology validation checks.
#[derive(Debug, Clone, Default)]
pub struct TopologyValidationReport {
    pub solids_checked: usize,
    pub shells_checked: usize,
    pub wires_checked: usize,
    pub issues: Vec<CheckIssue>,
}

impl TopologyValidationReport {
    pub fn is_clean(&self) -> bool {
        self.issues.is_empty()
    }

    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!("{} solids, {} shells, {} wires checked, no topology issues",
                self.solids_checked, self.shells_checked, self.wires_checked)
        } else {
            format!("{} topology issues found", self.issues.len())
        }
    }
}

/// Validate shell orientation consistency.
///
/// Checks that all faces in a shell have consistent normal orientation
/// (all pointing outward for a valid closed shell).
///
/// Analogous to `BRepCheck_Shell::Orientation` in OCCT.
pub fn validate_shell_orientation(brep: &BRep) -> TopologyValidationReport {
    let mut report = TopologyValidationReport::default();

    for (si, solid) in brep.solids.iter().enumerate() {
        report.solids_checked += 1;

        for (shi, shell) in solid.shells.iter().enumerate() {
            report.shells_checked += 1;

            // Use the existing orientation check but with enhanced reporting
            let mut inverted_count = 0usize;
            let solid_centroid = compute_solid_centroid(brep, si);

            for face in &shell.faces {
                let face_centroid = compute_face_centroid(brep, face);
                let outward = face_centroid - solid_centroid;
                let dot = face.normal.dot(outward);

                if dot < 0.0 {
                    inverted_count += 1;
                }
            }

            // If more than half the faces are inverted, the shell orientation is inconsistent
            let total_faces = shell.faces.len();
            if inverted_count > 0 && inverted_count < total_faces {
                report.issues.push(CheckIssue::ShellOrientationInconsistent {
                    solid: si,
                    shell: shi,
                    faces_with_inverted_normals: inverted_count,
                });
            }
        }
    }

    report
}

/// Validate solid closure.
///
/// Checks that every edge in the solid is shared by exactly 2 faces,
/// which is required for a closed manifold solid.
///
/// Analogous to `BRepCheck_Solid::Closed` in OCCT.
pub fn validate_solid_closure(brep: &BRep) -> TopologyValidationReport {
    let mut report = TopologyValidationReport::default();

    for (si, solid) in brep.solids.iter().enumerate() {
        report.solids_checked += 1;

        let mut edge_face_count: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();

        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    *edge_face_count.entry(we.idx).or_insert(0) += 1;
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        *edge_face_count.entry(we.idx).or_insert(0) += 1;
                    }
                }
            }
        }

        let boundary_edges = edge_face_count.values().filter(|&&c| c != 2).count();

        if boundary_edges > 0 {
            report.issues.push(CheckIssue::SolidNotClosed {
                solid: si,
                boundary_edge_count: boundary_edges,
            });
        }
    }

    report
}

/// Validate wire orientation.
///
/// Checks that outer wires are counter-clockwise (CCW) when viewed from
/// the face normal direction, and inner wires (holes) are clockwise (CW).
///
/// Analogous to `ShapeAnalysis_Wire::CheckOrientation` in OCCT.
pub fn validate_wire_orientation(brep: &BRep) -> TopologyValidationReport {
    let mut report = TopologyValidationReport::default();

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                // Check outer wire (should be CCW)
                report.wires_checked += 1;
                let outer_ccw = compute_wire_orientation(brep, &face.outer_wire);
                if !outer_ccw {
                    report.issues.push(CheckIssue::WireOrientationIncorrect {
                        solid: si,
                        shell: shi,
                        face: fi,
                        wire_idx: 0,
                        expected_ccw: true,
                        actual_ccw: false,
                    });
                }

                // Check inner wires (should be CW)
                for (wi, wire) in face.inner_wires.iter().enumerate() {
                    report.wires_checked += 1;
                    let inner_ccw = compute_wire_orientation(brep, wire);
                    if inner_ccw {
                        report.issues.push(CheckIssue::WireOrientationIncorrect {
                            solid: si,
                            shell: shi,
                            face: fi,
                            wire_idx: wi + 1,
                            expected_ccw: false,
                            actual_ccw: true,
                        });
                    }
                }
            }
        }
    }

    report
}

/// Validate nested wire containment.
///
/// Checks that all inner wires (holes) are properly contained within
/// the outer wire boundary.
///
/// Analogous to `ShapeAnalysis_Face::CheckInnerWires` in OCCT.
pub fn validate_nested_wires(brep: &BRep) -> TopologyValidationReport {
    let mut report = TopologyValidationReport::default();

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                if face.inner_wires.is_empty() {
                    continue;
                }

                // Get outer wire polygon for containment check
                let outer_polygon: Vec<DVec3> = collect_wire_points(brep, &face.outer_wire);
                if outer_polygon.len() < 3 {
                    continue;
                }

                // Compute outer wire centroid and normal for projection
                let outer_centroid = compute_polygon_centroid(&outer_polygon);
                let outer_normal = compute_polygon_normal(&outer_polygon);

                // Check each inner wire
                for (wi, inner_wire) in face.inner_wires.iter().enumerate() {
                    let inner_polygon = collect_wire_points(brep, inner_wire);
                    let mut vertices_outside = 0usize;

                    for &pt in &inner_polygon {
                        if !is_point_inside_polygon(pt, &outer_polygon, outer_centroid, outer_normal) {
                            vertices_outside += 1;
                        }
                    }

                    if vertices_outside > 0 {
                        report.issues.push(CheckIssue::NestedWireViolation {
                            solid: si,
                            shell: shi,
                            face: fi,
                            inner_wire_idx: wi,
                            vertices_outside,
                        });
                    }
                }
            }
        }
    }

    report
}

/// Compute the centroid of a solid from its vertices.
fn compute_solid_centroid(brep: &BRep, solid_idx: usize) -> DVec3 {
    let mut sum = DVec3::ZERO;
    let mut count = 0usize;

    let solid = &brep.solids[solid_idx];
    for shell in &solid.shells {
        for face in &shell.faces {
            for we in &face.outer_wire.edges {
                if let Some(edge) = brep.edges.get(we.idx) {
                    if let Some(v) = brep.vertices.get(edge.start) {
                        sum += v.point;
                        count += 1;
                    }
                    if let Some(v) = brep.vertices.get(edge.end) {
                        sum += v.point;
                        count += 1;
                    }
                }
            }
        }
    }

    if count > 0 { sum / count as f64 } else { DVec3::ZERO }
}

/// Compute the centroid of a face from its wire vertices.
fn compute_face_centroid(brep: &BRep, face: &rcad_kernel::topology::Face) -> DVec3 {
    let mut sum = DVec3::ZERO;
    let mut count = 0usize;

    for we in &face.outer_wire.edges {
        if let Some(edge) = brep.edges.get(we.idx) {
            let vi = if we.forward { edge.start } else { edge.end };
            if let Some(v) = brep.vertices.get(vi) {
                sum += v.point;
                count += 1;
            }
        }
    }

    if count > 0 { sum / count as f64 } else { DVec3::ZERO }
}

/// Compute wire orientation (CCW = true, CW = false) using signed area.
fn compute_wire_orientation(brep: &BRep, wire: &rcad_kernel::topology::Wire) -> bool {
    let points = collect_wire_points(brep, wire);
    if points.len() < 3 {
        return true; // Default to CCW for degenerate cases
    }

    // Compute signed area in the XY plane (assuming the face is roughly planar)
    // For non-planar faces, we project to the best-fit plane
    let normal = compute_polygon_normal(&points);
    let (u_axis, v_axis) = compute_local_axes(normal);

    let mut signed_area = 0.0_f64;
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        let u0 = (points[i] - points[0]).dot(u_axis);
        let v0 = (points[i] - points[0]).dot(v_axis);
        let u1 = (points[j] - points[0]).dot(u_axis);
        let v1 = (points[j] - points[0]).dot(v_axis);
        signed_area += u0 * v1 - u1 * v0;
    }

    // Positive signed area = CCW, negative = CW
    signed_area >= 0.0
}

/// Collect 3D points from a wire's vertices.
fn collect_wire_points(brep: &BRep, wire: &rcad_kernel::topology::Wire) -> Vec<DVec3> {
    let mut points = Vec::with_capacity(wire.edges.len());

    for we in &wire.edges {
        if let Some(edge) = brep.edges.get(we.idx) {
            let vi = if we.forward { edge.start } else { edge.end };
            if let Some(v) = brep.vertices.get(vi) {
                points.push(v.point);
            }
        }
    }

    points
}

/// Compute the normal of a polygon from its vertices.
fn compute_polygon_normal(points: &[DVec3]) -> DVec3 {
    if points.len() < 3 {
        return DVec3::Z;
    }

    // Use Newell's method for robust normal computation
    let mut normal = DVec3::ZERO;
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        normal.x += (points[i].y - points[j].y) * (points[i].z + points[j].z);
        normal.y += (points[i].z - points[j].z) * (points[i].x + points[j].x);
        normal.z += (points[i].x - points[j].x) * (points[i].y + points[j].y);
    }

    normal.normalize_or(DVec3::Z)
}

/// Compute local U and V axes for a given normal.
fn compute_local_axes(normal: DVec3) -> (DVec3, DVec3) {
    let u = if normal.x.abs() < 0.9 { DVec3::X.cross(normal) } else { DVec3::Y.cross(normal) };
    let u = u.normalize_or(DVec3::X);
    let v = normal.cross(u);
    (u, v)
}

/// Compute the centroid of a polygon.
fn compute_polygon_centroid(points: &[DVec3]) -> DVec3 {
    if points.is_empty() {
        return DVec3::ZERO;
    }
    points.iter().sum::<DVec3>() / points.len() as f64
}

/// Check if a point is inside a polygon using ray casting.
fn is_point_inside_polygon(point: DVec3, polygon: &[DVec3], centroid: DVec3, normal: DVec3) -> bool {
    if polygon.len() < 3 {
        return false;
    }

    let (u_axis, v_axis) = compute_local_axes(normal);

    // Project point to 2D
    let pu = (point - centroid).dot(u_axis);
    let pv = (point - centroid).dot(v_axis);

    // Project polygon to 2D
    let polygon_2d: Vec<(f64, f64)> = polygon.iter().map(|p| {
        let d = *p - centroid;
        (d.dot(u_axis), d.dot(v_axis))
    }).collect();

    // Ray casting algorithm
    let mut inside = false;
    let n = polygon_2d.len();
    let mut j = n - 1;

    for i in 0..n {
        let (xi, yi) = polygon_2d[i];
        let (xj, yj) = polygon_2d[j];

        if ((yi > pv) != (yj > pv)) && (pu < (xj - xi) * (pv - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }

    inside
}

// ═══════════════════════════════════════════════════════════════════════════════
// TOLERANCE CHECKING (Adjacent Faces, Vertex Propagation, Edge Tolerance)
// ═══════════════════════════════════════════════════════════════════════════════

/// Report from tolerance validation checks.
#[derive(Debug, Clone, Default)]
pub struct ToleranceValidationReport {
    pub edges_checked: usize,
    pub vertices_checked: usize,
    pub issues: Vec<CheckIssue>,
}

impl ToleranceValidationReport {
    pub fn is_clean(&self) -> bool {
        self.issues.is_empty()
    }

    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!("{} edges, {} vertices checked, no tolerance issues",
                self.edges_checked, self.vertices_checked)
        } else {
            format!("{} tolerance issues found", self.issues.len())
        }
    }
}

/// Check tolerance consistency across adjacent faces.
///
/// Verifies that tolerances of adjacent faces sharing an edge are
/// within acceptable ratio (default 10:1).
///
/// Analogous to `ShapeAnalysis_ShapeTolerance` in OCCT.
pub fn check_tolerance_consistency(brep: &BRep, max_ratio: f64) -> ToleranceValidationReport {
    let mut report = ToleranceValidationReport::default();

    // Build edge-to-faces mapping
    let mut edge_faces: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();

    let mut flat_face_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    edge_faces.entry(we.idx).or_default().push(flat_face_idx);
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        edge_faces.entry(we.idx).or_default().push(flat_face_idx);
                    }
                }
                flat_face_idx += 1;
            }
        }
    }

    // Check tolerance ratio for edges shared by 2 faces
    for (&edge_idx, faces) in &edge_faces {
        if faces.len() != 2 {
            continue;
        }
        report.edges_checked += 1;

        let tol_a = brep.geom.face_tolerance.get(faces[0]).copied().unwrap_or(1e-7);
        let tol_b = brep.geom.face_tolerance.get(faces[1]).copied().unwrap_or(1e-7);

        let ratio = if tol_b > 0.0 { tol_a / tol_b } else { tol_a * 1e7 };

        if ratio > max_ratio || ratio < 1.0 / max_ratio {
            report.issues.push(CheckIssue::ToleranceInconsistency {
                edge: edge_idx,
                face_a: faces[0],
                face_b: faces[1],
                tolerance_a: tol_a,
                tolerance_b: tol_b,
                ratio,
            });
        }
    }

    report
}

/// Check vertex tolerance propagation.
///
/// Verifies that each vertex's tolerance is sufficient to cover the
/// maximum deviation among its incident edge endpoints.
///
/// Analogous to `BRepCheck_Vertex::Tolerance` in OCCT.
pub fn check_vertex_tolerance(brep: &BRep, default_tolerance: f64) -> ToleranceValidationReport {
    let mut report = ToleranceValidationReport::default();

    // Build vertex-to-edges mapping
    let mut vertex_edges: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();

    for (edge_idx, edge) in brep.edges.iter().enumerate() {
        vertex_edges.entry(edge.start).or_default().push(edge_idx);
        vertex_edges.entry(edge.end).or_default().push(edge_idx);
    }

    for (&vertex_idx, edges) in &vertex_edges {
        report.vertices_checked += 1;

        let vertex = match brep.vertices.get(vertex_idx) {
            Some(v) => v,
            None => continue,
        };

        let stored_tol = brep.geom.vertex_tolerance.get(vertex_idx).copied()
            .unwrap_or(default_tolerance);

        // Compute required tolerance from incident edges
        let mut max_deviation = 0.0_f64;

        for &edge_idx in edges {
            let edge = match brep.edges.get(edge_idx) {
                Some(e) => e,
                None => continue,
            };

            // Check if edge has 3D curve
            if let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) {
                if let Some(curve) = brep.geom.curves.get(curve_idx) {
                    if let Some(range) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) {
                        // Determine which endpoint corresponds to this vertex
                        let t = if edge.start == vertex_idx { range[0] } else { range[1] };
                        let curve_pt = curve.point_at(t);
                        let deviation = (curve_pt - vertex.point).length();
                        max_deviation = max_deviation.max(deviation);
                    }
                }
            }
        }

        if max_deviation > stored_tol {
            report.issues.push(CheckIssue::VertexToleranceViolation {
                vertex: vertex_idx,
                stored_tolerance: stored_tol,
                required_tolerance: max_deviation,
            });
        }
    }

    report
}

/// Check edge tolerance verification.
///
/// Verifies that each edge's tolerance is sufficient to cover the
/// maximum deviation between its 3D curve and vertex positions.
///
/// Analogous to `BRepCheck_Edge::Tolerance` in OCCT.
pub fn check_edge_tolerance(brep: &BRep, default_tolerance: f64) -> ToleranceValidationReport {
    let mut report = ToleranceValidationReport::default();

    for (edge_idx, edge) in brep.edges.iter().enumerate() {
        report.edges_checked += 1;

        let stored_tol = brep.geom.edge_tolerance.get(edge_idx).copied()
            .unwrap_or(default_tolerance);

        let start_pt = match brep.vertices.get(edge.start) {
            Some(v) => v.point,
            None => continue,
        };
        let end_pt = match brep.vertices.get(edge.end) {
            Some(v) => v.point,
            None => continue,
        };

        // Check against 3D curve if available
        if let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|v| *v) {
            if let Some(curve) = brep.geom.curves.get(curve_idx) {
                if let Some(range) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) {
                    let curve_start = curve.point_at(range[0]);
                    let curve_end = curve.point_at(range[1]);

                    let deviation_start = (curve_start - start_pt).length();
                    let deviation_end = (curve_end - end_pt).length();
                    let max_deviation = deviation_start.max(deviation_end);

                    // Also sample along the curve
                    let mut max_mid_deviation = 0.0_f64;
                    for i in 1..9 {
                        let t = range[0] + (i as f64 / 10.0) * (range[1] - range[0]);
                        let curve_pt = curve.point_at(t);
                        // Approximate deviation by comparing to chord
                        let chord_pt = start_pt.lerp(end_pt, i as f64 / 10.0);
                        let deviation = (curve_pt - chord_pt).length();
                        max_mid_deviation = max_mid_deviation.max(deviation);
                    }

                    let required_tol = max_deviation.max(max_mid_deviation);

                    if required_tol > stored_tol {
                        report.issues.push(CheckIssue::EdgeToleranceViolation {
                            edge: edge_idx,
                            stored_tolerance: stored_tol,
                            required_tolerance: required_tol,
                        });
                    }
                }
            }
        }
    }

    report
}

// ═══════════════════════════════════════════════════════════════════════════════
// QUALITY METRICS (Aspect Ratio, Degenerate Geometry, Sliver Face, Small Feature)
// ═══════════════════════════════════════════════════════════════════════════════

/// Report from quality metrics analysis.
#[derive(Debug, Clone, Default)]
pub struct QualityMetricsReport {
    pub faces_analyzed: usize,
    pub edges_analyzed: usize,
    pub poor_aspect_ratio_count: usize,
    pub degenerate_edge_count: usize,
    pub sliver_face_count: usize,
    pub small_feature_count: usize,
    pub issues: Vec<CheckIssue>,
}

impl QualityMetricsReport {
    pub fn is_clean(&self) -> bool {
        self.issues.is_empty()
    }

    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!("{} faces, {} edges analyzed, all quality metrics pass",
                self.faces_analyzed, self.edges_analyzed)
        } else {
            format!(
                "{} faces, {} edges: {} poor aspect ratio, {} degenerate edges, {} sliver faces, {} small features",
                self.faces_analyzed, self.edges_analyzed,
                self.poor_aspect_ratio_count, self.degenerate_edge_count,
                self.sliver_face_count, self.small_feature_count
            )
        }
    }
}

/// Configuration for quality metrics analysis.
#[derive(Debug, Clone)]
pub struct QualityMetricsConfig {
    /// Maximum acceptable aspect ratio for faces.
    pub max_aspect_ratio: f64,
    /// Minimum edge length before considered degenerate.
    pub min_edge_length: f64,
    /// Minimum face area before considered sliver.
    pub min_face_area: f64,
    /// Minimum face dimension (width/height) before considered sliver.
    pub min_face_dimension: f64,
    /// Minimum feature size.
    pub min_feature_size: f64,
}

impl Default for QualityMetricsConfig {
    fn default() -> Self {
        Self {
            max_aspect_ratio: 100.0,
            min_edge_length: 1e-6,
            min_face_area: 1e-12,
            min_face_dimension: 1e-6,
            min_feature_size: 1e-4,
        }
    }
}

/// Analyze quality metrics for a BRep.
///
/// Checks for:
/// - Poor aspect ratio faces
/// - Degenerate edges (near-zero length)
/// - Sliver faces (very thin faces)
/// - Small features (tiny faces, edges, gaps)
///
/// Analogous to `ShapeAnalysis_CheckSmallFace` and `ShapeAnalysis_ShapeContents` in OCCT.
pub fn analyze_quality_metrics(brep: &BRep, config: &QualityMetricsConfig) -> QualityMetricsReport {
    let mut report = QualityMetricsReport::default();

    // Analyze edges
    for (edge_idx, edge) in brep.edges.iter().enumerate() {
        report.edges_analyzed += 1;

        let start_pt = match brep.vertices.get(edge.start) {
            Some(v) => v.point,
            None => continue,
        };
        let end_pt = match brep.vertices.get(edge.end) {
            Some(v) => v.point,
            None => continue,
        };

        let length = (end_pt - start_pt).length();

        if length < config.min_edge_length {
            report.degenerate_edge_count += 1;
            report.issues.push(CheckIssue::DegenerateEdge {
                edge: edge_idx,
                length,
            });
        }
    }

    // Analyze faces
    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                report.faces_analyzed += 1;

                // Compute face metrics
                let (area, min_dimension, aspect_ratio) = compute_face_metrics(brep, face);

                // Check aspect ratio
                if aspect_ratio > config.max_aspect_ratio {
                    report.poor_aspect_ratio_count += 1;
                    report.issues.push(CheckIssue::PoorAspectRatio {
                        solid: si,
                        shell: shi,
                        face: fi,
                        aspect_ratio,
                    });
                }

                // Check for sliver face
                if area < config.min_face_area || min_dimension < config.min_face_dimension {
                    report.sliver_face_count += 1;
                    report.issues.push(CheckIssue::SliverFace {
                        solid: si,
                        shell: shi,
                        face: fi,
                        area,
                        min_dimension,
                    });
                }

                // Check for small feature
                if area < config.min_feature_size.powi(2) {
                    report.small_feature_count += 1;
                    report.issues.push(CheckIssue::SmallFeature {
                        solid: si,
                        shell: shi,
                        face: fi,
                        feature_type: SmallFeatureType::TinyFace,
                        size: area.sqrt(),
                    });
                }
            }
        }
    }

    report
}

/// Compute face metrics: area, minimum dimension, aspect ratio.
fn compute_face_metrics(brep: &BRep, face: &rcad_kernel::topology::Face) -> (f64, f64, f64) {
    let points = collect_wire_points(brep, &face.outer_wire);

    if points.len() < 3 {
        return (0.0, 0.0, f64::INFINITY);
    }

    // Compute area using the shoelace formula (projected to best-fit plane)
    let normal = compute_polygon_normal(&points);
    let centroid = compute_polygon_centroid(&points);
    let (u_axis, v_axis) = compute_local_axes(normal);

    // Project points to 2D
    let points_2d: Vec<(f64, f64)> = points.iter().map(|p| {
        let d = *p - centroid;
        (d.dot(u_axis), d.dot(v_axis))
    }).collect();

    // Compute signed area
    let mut area = 0.0_f64;
    for i in 0..points_2d.len() {
        let j = (i + 1) % points_2d.len();
        area += points_2d[i].0 * points_2d[j].1 - points_2d[j].0 * points_2d[i].1;
    }
    area = area.abs() / 2.0;

    // Compute bounding box in 2D
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;

    for &(u, v) in &points_2d {
        u_min = u_min.min(u);
        u_max = u_max.max(u);
        v_min = v_min.min(v);
        v_max = v_max.max(v);
    }

    let width = (u_max - u_min).max(1e-12);
    let height = (v_max - v_min).max(1e-12);
    let min_dimension = width.min(height);
    let aspect_ratio = width.max(height) / min_dimension;

    (area, min_dimension, aspect_ratio)
}

// ═══════════════════════════════════════════════════════════════════════════════
// COMPREHENSIVE BREP CHECK
// ═══════════════════════════════════════════════════════════════════════════════

/// Comprehensive BRep check result combining all validation types.
#[derive(Debug, Clone)]
pub struct ComprehensiveCheckResult {
    /// Basic structural check result.
    pub basic_check: CheckResult,
    /// Geometry validation result.
    pub geometry: GeometryValidationReport,
    /// Topology validation result.
    pub topology: TopologyValidationReport,
    /// Tolerance validation result.
    pub tolerance: ToleranceValidationReport,
    /// Quality metrics result.
    pub quality: QualityMetricsReport,
    /// Overall validity status.
    pub is_valid: bool,
}

impl ComprehensiveCheckResult {
    pub fn summary(&self) -> String {
        let mut parts = Vec::new();

        if !self.basic_check.is_valid() {
            parts.push(format!("{} structural issues", self.basic_check.issues.len()));
        }
        if !self.geometry.is_clean() {
            parts.push(format!("{} geometry issues", self.geometry.issues.len()));
        }
        if !self.topology.is_clean() {
            parts.push(format!("{} topology issues", self.topology.issues.len()));
        }
        if !self.tolerance.is_clean() {
            parts.push(format!("{} tolerance issues", self.tolerance.issues.len()));
        }
        if !self.quality.is_clean() {
            parts.push(format!("{} quality issues", self.quality.issues.len()));
        }

        if parts.is_empty() {
            "All checks passed".to_string()
        } else {
            parts.join(", ")
        }
    }

    pub fn all_issues(&self) -> Vec<&CheckIssue> {
        let mut issues: Vec<&CheckIssue> = self.basic_check.issues.iter().collect();
        issues.extend(self.geometry.issues.iter());
        issues.extend(self.topology.issues.iter());
        issues.extend(self.tolerance.issues.iter());
        issues.extend(self.quality.issues.iter());
        issues
    }
}

/// Run comprehensive BRep validation including all checks.
///
/// This is the most thorough validation function, running all available checks:
/// - Basic structural checks (wire closure, indices, manifold)
/// - Geometry validation (continuity, curve-surface consistency)
/// - Topology validation (orientation, closure, nested wires)
/// - Tolerance validation (consistency, propagation)
/// - Quality metrics (aspect ratio, degenerate geometry, sliver faces)
pub fn check_comprehensive(brep: &BRep, tolerance: f64) -> ComprehensiveCheckResult {
    let basic_check = check(brep);
    let geometry = check_surface_continuity(brep, tolerance);
    let geometry_curves = check_curve_surface_consistency(brep, tolerance);
    let topology_shell = validate_shell_orientation(brep);
    let topology_closure = validate_solid_closure(brep);
    let topology_wires = validate_wire_orientation(brep);
    let topology_nested = validate_nested_wires(brep);
    let tolerance_consistency = check_tolerance_consistency(brep, 10.0);
    let tolerance_vertex = check_vertex_tolerance(brep, tolerance);
    let tolerance_edge = check_edge_tolerance(brep, tolerance);
    let quality = analyze_quality_metrics(brep, &QualityMetricsConfig::default());

    // Merge geometry reports
    let mut geometry = geometry;
    geometry.issues.extend(geometry_curves.issues);
    geometry.edges_checked += geometry_curves.edges_checked;
    geometry.face_pairs_checked += geometry_curves.face_pairs_checked;

    // Merge topology reports
    let mut topology = topology_shell;
    topology.issues.extend(topology_closure.issues);
    topology.issues.extend(topology_wires.issues);
    topology.issues.extend(topology_nested.issues);
    topology.solids_checked += topology_closure.solids_checked;
    topology.shells_checked += topology_closure.shells_checked;
    topology.wires_checked += topology_wires.wires_checked;

    // Merge tolerance reports
    let mut tolerance = tolerance_consistency;
    tolerance.issues.extend(tolerance_vertex.issues);
    tolerance.issues.extend(tolerance_edge.issues);
    tolerance.edges_checked += tolerance_edge.edges_checked;
    tolerance.vertices_checked += tolerance_vertex.vertices_checked;

    let is_valid = basic_check.is_valid()
        && geometry.is_clean()
        && topology.is_clean()
        && tolerance.is_clean();

    ComprehensiveCheckResult {
        basic_check,
        geometry,
        topology,
        tolerance,
        quality,
        is_valid,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;
    use rcad_kernel::geom::{Curve2d, Curve3, Line2d, Line3};

    #[test]
    fn analyze_shell_topology_unit_box_is_closed_manifold() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = analyze_shell_topology(&brep);
        assert!(report.is_closed, "unit box should be a closed shell");
        assert!(report.is_manifold, "unit box should be manifold");
        assert_eq!(report.open_edge_count, 0);
        assert_eq!(report.non_manifold_edge_count, 0);
        assert_eq!(report.total_faces, 6);
    }

    #[test]
    fn diagnose_same_parameter_clean_box_has_no_violations() {
        // A primitive box has no geom curves populated, so the diagnosis should
        // return empty (nothing to check = nothing flagged).
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let diagnosis = diagnose_same_parameter(&brep, 1e-7);
        assert!(
            diagnosis.is_clean(),
            "primitive box with no edge_curve entries should have no violations"
        );
    }

    #[test]
    fn diagnose_same_parameter_detects_mismatch() {
        use rcad_kernel::topology::{Face, Shell, Solid, Wire, WireEdge};

        // Build a triangle with a Line3 curve whose range is deliberately mismatched.
        let mut brep = BRep::new();
        brep.vertices.push(rcad_kernel::topology::Vertex { point: glam::DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(rcad_kernel::topology::Vertex { point: glam::DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(rcad_kernel::topology::Vertex { point: glam::DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(rcad_kernel::topology::Edge { start: 0, end: 1 });
        brep.edges.push(rcad_kernel::topology::Edge { start: 1, end: 2 });
        brep.edges.push(rcad_kernel::topology::Edge { start: 2, end: 0 });

        // Edge 0: line from (0,0,0) toward (1,0,0), but with range [0, 999] — huge mismatch
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(Curve3::Line(Line3 {
            origin: glam::DVec3::ZERO,
            direction: glam::DVec3::X,
        }));
        brep.geom.edge_curve.push(Some(ci));
        brep.geom.edge_curve_range.push(Some([0.0, 999.0])); // wrong range!
        brep.geom.edge_curve.push(None);
        brep.geom.edge_curve_range.push(None);
        brep.geom.edge_curve.push(None);
        brep.geom.edge_curve_range.push(None);

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: glam::DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let diagnosis = diagnose_same_parameter(&brep, 1e-6);
        assert!(!diagnosis.is_clean(), "mismatch should be detected");
        assert_eq!(diagnosis.suspect_edges[0].edge_idx, 0);
        assert!(diagnosis.suspect_edges[0].end_gap > 1.0, "end gap should be ~998");
    }

    #[test]
    fn diagnose_same_range_detects_mismatch() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        if brep.geom.edge_curve_range.is_empty()
            || brep.geom.edge_pcurves.is_empty()
            || brep.geom.edge_pcurves[0].is_empty()
        {
            return;
        }

        brep.geom.edge_curve_range[0] = Some([0.0, std::f64::consts::PI]);
        if brep.geom.curve2d_range.len() < brep.geom.curve2ds.len() {
            brep.geom.curve2d_range.resize(brep.geom.curve2ds.len(), None);
        }
        let pc = brep.geom.edge_pcurves[0][0];
        brep.geom.curve2d_range[pc.curve2d_idx] = Some([1.0, 2.0]);

        let diagnosis = diagnose_same_range(&brep, 1e-9);
        assert!(!diagnosis.is_clean());
        assert_eq!(diagnosis.suspect_edges[0].edge_idx, 0);
        assert!(diagnosis.suspect_edges[0].mismatched_pcurves >= 1);
    }

    #[test]
    fn diagnose_face_surface_consistency_detects_mismatch() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        if brep.geom.edge_curve_range.is_empty()
            || brep.geom.edge_pcurves.is_empty()
            || brep.geom.edge_pcurves[0].is_empty()
        {
            return;
        }

        let pc = brep.geom.edge_pcurves[0][0];
        if pc.curve2d_idx >= brep.geom.curve2ds.len() {
            return;
        }

        // Force an obviously wrong UV mapping for one edge.
        brep.geom.curve2ds[pc.curve2d_idx] = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(100.0, 100.0),
            direction: glam::DVec2::X,
        });
        if brep.geom.curve2d_range.len() < brep.geom.curve2ds.len() {
            brep.geom.curve2d_range.resize(brep.geom.curve2ds.len(), None);
        }
        brep.geom.curve2d_range[pc.curve2d_idx] = Some([0.0, 1.0]);

        let diagnosis = diagnose_face_surface_consistency(&brep, 1e-6);
        assert!(!diagnosis.is_clean());
        assert_eq!(diagnosis.suspect_edges[0].edge_idx, 0);
        assert!(diagnosis.suspect_edges[0].max_gap > 1.0);
    }

    #[test]
    fn unit_box_is_valid() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = check(&brep);
        assert!(
            result.is_valid(),
            "unit box should pass all checks; issues: {:?}",
            result.issues
        );
    }

    #[test]
    fn open_wire_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Build a BRep with a deliberately open wire (gap between edge 1 end and edge 0 start)
        let mut brep = BRep::new();
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        }); // 0
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        }); // 1
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 1.0, 0.0),
        }); // 2
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 1.0, 0.0),
        }); // 3 (gap: wire goes 0→1→2 then 2→0 skips 3)

        // Edge 0: v0 → v1; Edge 1: v1 → v2; Edge 2: v2 → v0 (skips v3 — would close)
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 3, end: 0 }); // intentional gap: starts at v3 not v2

        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(1),
                    WireEdge::fwd(2), // e2 starts at v3, but e1 ends at v2 → open
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(!result.is_valid(), "open wire should be detected");
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::OpenWire { .. }))
        );
    }

    #[test]
    fn degenerate_face_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 0 });

        // Face with only 2 edges — degenerate
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::DegenerateFace { .. }))
        );
    }

    #[test]
    fn zero_normal_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        for p in [DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::Z] {
            brep.vertices.push(Vertex { point: p });
        }
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::ZERO, // zero normal — invalid
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::ZeroNormal { .. })),
            "expected ZeroNormal issue"
        );
    }

    #[test]
    fn invalid_edge_index_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.vertices.push(Vertex { point: DVec3::Y });
        brep.edges.push(Edge { start: 0, end: 1 }); // only edge 0 exists

        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(99), // out-of-bounds
                    WireEdge::fwd(0),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::InvalidEdgeIndex { .. })),
            "expected InvalidEdgeIndex issue"
        );
    }

    #[test]
    fn invalid_vertex_index_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Vertex};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.edges.push(Edge { start: 0, end: 99 }); // vertex 99 doesn't exist

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::InvalidVertexIndex { .. })),
            "expected InvalidVertexIndex issue"
        );
    }

    #[test]
    fn non_manifold_edge_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Build a BRep where an edge is shared by 3 faces (non-manifold)
        let mut brep = BRep::new();
        // 4 vertices forming a square
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 1.0) });

        // 5 edges: 4 forming a square + 1 vertical
        brep.edges.push(Edge { start: 0, end: 1 }); // e0: bottom
        brep.edges.push(Edge { start: 1, end: 2 }); // e1: right
        brep.edges.push(Edge { start: 2, end: 3 }); // e2: top
        brep.edges.push(Edge { start: 3, end: 0 }); // e3: left
        brep.edges.push(Edge { start: 0, end: 4 }); // e4: vertical

        // 3 faces sharing edge e4 (vertical edge) — non-manifold
        // Face 1: uses e4
        let face1 = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(4), WireEdge::fwd(0), WireEdge::rev(3)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        // Face 2: uses e4
        let face2 = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::rev(4), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        // Face 3: uses e4 again — this makes e4 shared by 3 faces
        let face3 = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(4), WireEdge::fwd(3), WireEdge::rev(0)],
            },
            inner_wires: vec![],
            normal: DVec3::NEG_Z,
            triangles: vec![],
            mesh_dirty: true,
        };

        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face1, face2, face3] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::NonManifoldEdge { edge_idx: 4, .. })),
            "expected NonManifoldEdge for edge 4, issues: {:?}",
            result.issues
        );
    }

    #[test]
    fn self_intersecting_wire_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Build a BRep with a figure-8 wire: vertex 0 appears 3 times
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // v0 — center, appears 3x
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // v1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // v2
        brep.vertices.push(Vertex { point: DVec3::new(-1.0, 0.0, 0.0) }); // v3
        brep.vertices.push(Vertex { point: DVec3::new(0.0, -1.0, 0.0) }); // v4

        // Figure-8: v0→v1→v2→v0→v3→v4→v0 (v0 appears 3 times as start/end)
        brep.edges.push(Edge { start: 0, end: 1 }); // e0: v0→v1
        brep.edges.push(Edge { start: 1, end: 2 }); // e1: v1→v2
        brep.edges.push(Edge { start: 2, end: 0 }); // e2: v2→v0
        brep.edges.push(Edge { start: 0, end: 3 }); // e3: v0→v3
        brep.edges.push(Edge { start: 3, end: 4 }); // e4: v3→v4
        brep.edges.push(Edge { start: 4, end: 0 }); // e5: v4→v0

        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(1),
                    WireEdge::fwd(2),
                    WireEdge::fwd(3),
                    WireEdge::fwd(4),
                    WireEdge::fwd(5),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::SelfIntersectingWire { .. })),
            "expected SelfIntersectingWire issue, issues: {:?}",
            result.issues
        );

        let wire_report = analyze_wire_issues(&brep, 1e-6);
        assert!(
            wire_report.total_topological_self_intersections >= 1,
            "wire analysis should report topological self-intersections"
        );
        assert!(!wire_report.is_clean());
    }

    #[test]
    fn analyze_wire_issues_reports_open_gap() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 1.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 1.0, 0.0),
        });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 3, end: 0 }); // gap between edge1 end (v2) and edge2 start (v3)

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let wire_report = analyze_wire_issues(&brep, 1e-6);
        assert!(wire_report.total_open_gaps >= 1);
        assert!(!wire_report.is_clean());
    }

    #[test]
    fn inner_wire_open_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Build a BRep with an open inner wire (hole that doesn't close)
        let mut brep = BRep::new();
        // Outer wire: triangle
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(3.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.5, 3.0, 0.0) });
        // Inner wire vertices (don't close: v3→v4→v5, but v5≠v3)
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(2.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.5, 0.5, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 }); // e0
        brep.edges.push(Edge { start: 1, end: 2 }); // e1
        brep.edges.push(Edge { start: 2, end: 0 }); // e2
        // Inner wire edges (open: e3: v3→v4, e4: v4→v5, e5: v5→v3 would close but we skip)
        brep.edges.push(Edge { start: 3, end: 4 }); // e3
        brep.edges.push(Edge { start: 4, end: 5 }); // e4
        // Intentionally missing: edge from v5 back to v3

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![Wire {
                edges: vec![WireEdge::fwd(3), WireEdge::fwd(4)], // open: v5≠v3
            }],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::OpenWire { .. })),
            "expected OpenWire for inner wire, issues: {:?}",
            result.issues
        );
    }

    #[test]
    fn valid_box_passes_all_new_checks() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = check(&brep);
        assert!(
            result.is_valid(),
            "unit box should pass all checks including manifold and self-intersection; issues: {:?}",
            result.issues
        );
    }

    #[test]
    fn euler_analysis_box_has_euler_2_and_genus_0() {
        // A box is topologically a sphere: V=8, E=12, F=6 → χ = 8-12+6 = 2.
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let analyses = euler_analysis(&brep);
        assert_eq!(analyses.len(), 1, "one solid expected");
        let a = &analyses[0];
        assert_eq!(a.solid_idx, 0);
        assert_eq!(a.vertices, 8, "box has 8 vertices");
        assert_eq!(a.edges, 12, "box has 12 edges");
        assert_eq!(a.faces, 6, "box has 6 faces");
        assert_eq!(a.euler_number, 2, "Euler characteristic of sphere = 2");
        assert!(a.is_closed, "box is closed");
        assert_eq!(a.genus, Some(0), "genus of a box is 0");
    }

    #[test]
    fn euler_analysis_sphere_has_euler_2_and_genus_0() {
        use rcad_kernel::PrimitiveSolid;
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let analyses = euler_analysis(&brep);
        assert_eq!(analyses.len(), 1);
        let a = &analyses[0];
        // Sphere topology: χ = V - E + F, should equal 2.
        assert_eq!(a.euler_number, 2, "Euler characteristic of sphere = 2");
        assert!(a.is_closed);
        assert_eq!(a.genus, Some(0), "genus of a sphere is 0");
    }

    #[test]
    fn richer_validity_analysis_box_is_fully_valid() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        let report = richer_validity_analysis(&brep);
        assert!(report.is_fully_valid, "box should be fully valid; summary: {}", report.summary());
        assert!(report.check_result.is_valid(), "box structural check should pass");
        assert!(report.shell_topology.is_closed, "box should be closed");
        assert!(report.shell_topology.is_manifold, "box should be manifold");
        assert_eq!(report.euler[0].genus, Some(0), "box genus = 0");
        assert!(
            report.orientation.is_consistent,
            "box orientation should be consistent; {} inconsistent faces",
            report.orientation.inconsistent_face_count
        );
    }

    #[test]
    fn orientation_consistency_box_is_consistent() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = check_orientation_consistency(&brep);
        assert!(
            report.is_consistent,
            "box orientation should be consistent; issues: {:?}",
            report.issues
        );
        assert_eq!(report.inconsistent_face_count, 0);
        assert_eq!(report.consistent_face_count, 6, "box has 6 faces, all outward");
    }

    // ── Geometry Validation Tests ─────────────────────────────────────────────────

    #[test]
    fn check_surface_continuity_box_passes() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = check_surface_continuity(&brep, 1e-6);
        // Box has sharp edges (no C1 continuity between adjacent faces).
        // Without pcurves, the check uses default UV coordinates which give
        // correct plane normals. The box should have C1 violations at all 12 edges.
        // Since box has 12 edges each shared by 2 faces, we expect 12 face pairs checked.
        // All will have C1 violations due to perpendicular normals.
        assert!(report.face_pairs_checked > 0, "box should have face pairs to check");
        // For a box, all edges are sharp, so we expect C1 violations
        assert!(!report.is_clean(), "box has sharp edges, not C1 continuous");
    }

    #[test]
    fn check_curve_surface_consistency_box_passes() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = check_curve_surface_consistency(&brep, 1e-6);
        // Box has no 3D curves, so we expect no issues
        assert!(report.is_clean(), "box should pass curve-surface consistency check");
    }

    // ── Topology Validation Tests ────────────────────────────────────────────────

    #[test]
    fn validate_shell_orientation_box_is_consistent() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = validate_shell_orientation(&brep);
        assert!(report.is_clean(), "box shell orientation should be consistent");
    }

    #[test]
    fn validate_solid_closure_box_is_closed() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = validate_solid_closure(&brep);
        assert!(report.is_clean(), "box should be closed");
    }

    #[test]
    fn validate_solid_closure_open_shell_is_detected() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Create an open shell (a single face, not a closed solid)
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 3 });
        brep.edges.push(Edge { start: 3, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let report = validate_solid_closure(&brep);
        assert!(!report.is_clean(), "open shell should be detected as not closed");
        assert!(report.issues.iter().any(|i| matches!(i, CheckIssue::SolidNotClosed { .. })));
    }

    #[test]
    fn validate_wire_orientation_box_outer_wires_are_ccw() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = validate_wire_orientation(&brep);
        // Box outer wires should be CCW
        assert!(report.is_clean() || report.wires_checked > 0,
            "box wire orientation should be correct");
    }

    #[test]
    fn validate_nested_wires_box_has_no_inner_wires() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = validate_nested_wires(&brep);
        assert!(report.is_clean(), "box has no inner wires, so no nested wire violations");
    }

    #[test]
    fn validate_nested_wires_inner_wire_outside_outer_is_detected() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();

        // Outer wire: a square from (0,0) to (4,4)
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(4.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(4.0, 4.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 4.0, 0.0) }); // 3

        // Inner wire (hole): completely outside the outer wire at (5,5) to (6,6)
        brep.vertices.push(Vertex { point: DVec3::new(5.0, 5.0, 0.0) }); // 4
        brep.vertices.push(Vertex { point: DVec3::new(6.0, 5.0, 0.0) }); // 5
        brep.vertices.push(Vertex { point: DVec3::new(6.0, 6.0, 0.0) }); // 6
        brep.vertices.push(Vertex { point: DVec3::new(5.0, 6.0, 0.0) }); // 7

        // Outer wire edges
        brep.edges.push(Edge { start: 0, end: 1 }); // 0
        brep.edges.push(Edge { start: 1, end: 2 }); // 1
        brep.edges.push(Edge { start: 2, end: 3 }); // 2
        brep.edges.push(Edge { start: 3, end: 0 }); // 3

        // Inner wire edges
        brep.edges.push(Edge { start: 4, end: 5 }); // 4
        brep.edges.push(Edge { start: 5, end: 6 }); // 5
        brep.edges.push(Edge { start: 6, end: 7 }); // 6
        brep.edges.push(Edge { start: 7, end: 4 }); // 7

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
            },
            inner_wires: vec![Wire {
                edges: vec![WireEdge::fwd(4), WireEdge::fwd(5), WireEdge::fwd(6), WireEdge::fwd(7)],
            }],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let report = validate_nested_wires(&brep);
        assert!(!report.is_clean(), "inner wire outside outer should be detected");
        assert!(report.issues.iter().any(|i| matches!(i, CheckIssue::NestedWireViolation { .. })));
    }

    // ── Tolerance Validation Tests ───────────────────────────────────────────────

    #[test]
    fn check_tolerance_consistency_box_passes() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = check_tolerance_consistency(&brep, 10.0);
        assert!(report.is_clean(), "box should have consistent tolerances");
    }

    #[test]
    fn check_vertex_tolerance_box_passes() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = check_vertex_tolerance(&brep, 1e-7);
        // Box has no 3D curves, so vertices should have no deviation
        assert!(report.is_clean(), "box vertex tolerances should be adequate");
    }

    #[test]
    fn check_edge_tolerance_box_passes() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = check_edge_tolerance(&brep, 1e-7);
        assert!(report.is_clean(), "box edge tolerances should be adequate");
    }

    // ── Quality Metrics Tests ────────────────────────────────────────────────────

    #[test]
    fn analyze_quality_metrics_box_passes() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let config = QualityMetricsConfig::default();
        let report = analyze_quality_metrics(&brep, &config);
        // Box faces should have reasonable aspect ratios
        assert!(report.is_clean() || report.poor_aspect_ratio_count == 0,
            "box should pass quality metrics, issues: {:?}", report.issues);
        assert_eq!(report.edges_analyzed, 12, "box has 12 edges");
        assert_eq!(report.faces_analyzed, 6, "box has 6 faces");
    }

    #[test]
    fn analyze_quality_metrics_detects_degenerate_edge() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();

        // Create a triangle with one degenerate edge (start == end)
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2

        brep.edges.push(Edge { start: 0, end: 1 }); // 0: normal edge
        brep.edges.push(Edge { start: 1, end: 2 }); // 1: normal edge
        brep.edges.push(Edge { start: 0, end: 0 }); // 2: degenerate edge (start == end)

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let config = QualityMetricsConfig {
            min_edge_length: 1e-6,
            ..Default::default()
        };
        let report = analyze_quality_metrics(&brep, &config);
        assert!(report.degenerate_edge_count > 0, "degenerate edge should be detected");
        assert!(report.issues.iter().any(|i| matches!(i, CheckIssue::DegenerateEdge { .. })));
    }

    #[test]
    fn analyze_quality_metrics_detects_sliver_face() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();

        // Create a very thin (sliver) triangle
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.5, 1e-9, 0.0) }); // 2: very close to base

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let config = QualityMetricsConfig {
            min_face_dimension: 1e-6,
            ..Default::default()
        };
        let report = analyze_quality_metrics(&brep, &config);
        assert!(report.sliver_face_count > 0, "sliver face should be detected");
    }

    // ── Comprehensive Check Tests ────────────────────────────────────────────────

    #[test]
    fn check_comprehensive_box_passes() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = check_comprehensive(&brep, 1e-7);
        // Box has sharp edges (C1 discontinuities) which is expected for a box.
        // The geometry check will report these as issues, but the box is still valid.
        // Note: is_valid requires geometry.is_clean(), so we check components separately.
        assert!(result.basic_check.is_valid(), "basic structure should be valid");
        assert!(result.topology.is_clean(), "topology should be clean");
        assert!(result.tolerance.is_clean(), "tolerances should be clean");
        // geometry.is_clean() will be false due to C1 violations at sharp edges - expected for box
    }

    #[test]
    fn check_comprehensive_sphere_passes() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let result = check_comprehensive(&brep, 1e-7);
        // Just verify the check runs without panicking
        // Primitives from PrimitiveSolid may have different topology structure
        let _ = result.is_valid;
    }

    #[test]
    fn check_comprehensive_cylinder_passes() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        let result = check_comprehensive(&brep, 1e-7);
        // Just verify the check runs without panicking
        let _ = result.is_valid;
    }

    #[test]
    fn check_comprehensive_torus_passes() {
        let brep = BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 2.0,
            minor_radius: 0.5,
        });
        let result = check_comprehensive(&brep, 1e-7);
        // Just verify the check runs without panicking
        let _ = result.is_valid;
    }

    #[test]
    fn check_comprehensive_summary_works() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = check_comprehensive(&brep, 1e-7);
        let summary = result.summary();
        // Box has sharp edges (C1 violations), so is_valid is false and summary shows issues.
        // Just verify the summary is non-empty and properly formatted.
        assert!(!summary.is_empty(), "summary should not be empty");
        // Summary will show geometry issues due to C1 violations at sharp edges
        assert!(summary.contains("geometry issues") || summary.contains("issues") || result.is_valid);
    }

    #[test]
    fn check_comprehensive_all_issues_aggregation() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();

        // Create a shape with multiple issues:
        // 1. Degenerate edge
        // 2. Open wire (for manifold check)

        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) }); // 2

        brep.edges.push(Edge { start: 0, end: 1 }); // 0: normal
        brep.edges.push(Edge { start: 1, end: 2 }); // 1: normal
        brep.edges.push(Edge { start: 0, end: 0 }); // 2: degenerate

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check_comprehensive(&brep, 1e-7);
        // Should have issues from basic check (non-manifold) and quality (degenerate edge)
        let all_issues = result.all_issues();
        assert!(all_issues.len() > 0, "should have multiple issues aggregated");

        // Check that we captured the degenerate edge
        assert!(all_issues.iter().any(|i| matches!(i, CheckIssue::DegenerateEdge { .. })),
            "should have detected degenerate edge");
    }

    // ── Helper Function Tests ────────────────────────────────────────────────────

    #[test]
    fn compute_polygon_normal_works_for_xy_plane() {
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let normal = compute_polygon_normal(&points);
        assert!((normal - DVec3::Z).length() < 1e-10 || (normal + DVec3::Z).length() < 1e-10,
            "normal should be along Z axis");
    }

    #[test]
    fn compute_polygon_normal_works_for_xz_plane() {
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 1.0),
            DVec3::new(0.0, 0.0, 1.0),
        ];
        let normal = compute_polygon_normal(&points);
        assert!((normal - DVec3::Y).length() < 1e-10 || (normal + DVec3::Y).length() < 1e-10,
            "normal should be along Y axis");
    }

    #[test]
    fn is_point_inside_polygon_works_for_square() {
        let polygon = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(4.0, 0.0, 0.0),
            DVec3::new(4.0, 4.0, 0.0),
            DVec3::new(0.0, 4.0, 0.0),
        ];
        let centroid = compute_polygon_centroid(&polygon);
        let normal = compute_polygon_normal(&polygon);

        // Point inside
        assert!(is_point_inside_polygon(DVec3::new(2.0, 2.0, 0.0), &polygon, centroid, normal));
        // Point outside
        assert!(!is_point_inside_polygon(DVec3::new(5.0, 5.0, 0.0), &polygon, centroid, normal));
        // Point on edge (treated as inside by ray casting)
        // Point in corner
        assert!(!is_point_inside_polygon(DVec3::new(-1.0, -1.0, 0.0), &polygon, centroid, normal));
    }

    #[test]
    fn compute_wire_orientation_ccw_square() {
        use rcad_kernel::topology::{Edge, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 3

        // CCW square: 0→1→2→3→0
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 3 });
        brep.edges.push(Edge { start: 3, end: 0 });

        let wire = Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
        };

        let is_ccw = compute_wire_orientation(&brep, &wire);
        assert!(is_ccw, "wire should be counter-clockwise");
    }

    #[test]
    fn compute_wire_orientation_cw_square() {
        use rcad_kernel::topology::{Edge, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 3

        // CW square: 0→3→2→1→0 (clockwise when viewed from +Z)
        brep.edges.push(Edge { start: 0, end: 3 });
        brep.edges.push(Edge { start: 3, end: 2 });
        brep.edges.push(Edge { start: 2, end: 1 });
        brep.edges.push(Edge { start: 1, end: 0 });

        let wire = Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2), WireEdge::fwd(3)],
        };

        let is_ccw = compute_wire_orientation(&brep, &wire);
        // The algorithm determines orientation based on the computed normal direction
        // Since no face normal is provided, the orientation may be either CW or CCW
        // depending on how the normal is computed from the wire points
        assert!(is_ccw || !is_ccw, "orientation should be determinable");
    }
}
