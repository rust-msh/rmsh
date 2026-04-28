use std::collections::HashMap;
use std::collections::HashSet;

use glam::{DVec2, DVec3};
use rayon::prelude::*;
use rcad_kernel::closest_point_on_curve;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve2dEval, SurfaceEval, *};
use rcad_kernel::topology::*;

use crate::bopds::ds::*;
use crate::classify::{Classification, classify_point};
use crate::history::{
    BooleanHistory, EdgeOrigin, FaceOrigin, HistoryTracker, ShellOrigin, SolidOrigin, VertexOrigin,
};
use crate::inttools::edge_face::plane_local_basis;
use crate::tolerance::*;
use crate::triangulate::triangulate_polygon;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOpType {
    Union,
    Intersection,
    Difference,
}

#[derive(Debug)]
pub enum BooleanError {
    EmptyInput,
    MissingGeometry(&'static str),
    DegenerateResult,
    /// A numeric operation produced a non-finite or NaN value.
    NumericalFailure(&'static str),
    /// An expected non-empty collection was empty (e.g. polyline with no points).
    EmptyCollection(&'static str),
    /// Result fails validity checks (non-manifold, open shells, invalid orientation).
    InvalidResult(&'static str),
    /// Missing intersection curves between surfaces that should intersect.
    IncompleteIntersection(&'static str),
    /// Result contains self-intersecting geometry.
    SelfIntersection(&'static str),
}

impl std::fmt::Display for BooleanError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "empty input"),
            Self::MissingGeometry(msg) => write!(f, "missing geometry: {msg}"),
            Self::DegenerateResult => write!(f, "degenerate result"),
            Self::NumericalFailure(msg) => write!(f, "numerical failure: {msg}"),
            Self::EmptyCollection(msg) => write!(f, "unexpected empty collection: {msg}"),
            Self::InvalidResult(msg) => write!(f, "invalid result: {msg}"),
            Self::IncompleteIntersection(msg) => write!(f, "incomplete intersection: {msg}"),
            Self::SelfIntersection(msg) => write!(f, "self-intersection: {msg}"),
        }
    }
}

impl std::error::Error for BooleanError {}

/// A sub-region of an original face after splitting by intersection curves.
#[derive(Debug, Clone)]
pub struct SubFace {
    /// Boundary vertex positions in 3D (ordered polygon).
    pub boundary: Vec<DVec3>,
    /// The surface this lies on.
    pub surface: Surface3,
    /// Normal direction.
    pub normal: DVec3,
    /// UV centroid of this sub-face's parameter-space polygon (for curved surfaces).
    /// Used by `sample_point` to produce a geometrically representative interior point.
    pub uv_centroid: Option<DVec2>,
    /// Explicit override for the sample point. When set, `sample_point()` uses this
    /// instead of computing it from the boundary centroid. Used when the centroid would
    /// fall in a different classification region (e.g. the outer annular region around
    /// an embedded circle, whose centroid falls inside the circle).
    pub sample_override: Option<DVec3>,
    /// UV domain [u0, u1, v0, v1] of this sub-face's parameter-space region.
    /// Propagated to `GeomStore.face_surface_range` in the result BRep so that
    /// `tessellate_curved_face` uses the correct sub-domain instead of the full
    /// surface domain.
    pub uv_domain: Option<[f64; 4]>,
    /// Inner wire boundaries (holes) in 3D. Each inner wire is an ordered polygon
    /// representing a closed trim curve that forms a hole in the face.
    pub inner_wires: Vec<Vec<DVec3>>,
    /// Optional analytic circular outer boundary. When present, the result
    /// builder emits a single closed edge with circle geometry instead of
    /// polyline segments.
    pub analytic_outer_circle: Option<Circle3>,
}

impl SubFace {
    fn sample_point(&self) -> DVec3 {
        // Returns a point slightly INSIDE the surface (toward the interior of the solid),
        // so classify_point can tell whether this sub-face is inside or outside
        // the other solid.
        //
        // For sphere sub-faces the outward normal points AWAY from the sphere center,
        // so we must offset toward the center to stay inside the sphere's volume.
        // We use the UV centroid to get a point in the middle of the spherical cap.
        if let Some(pt) = self.sample_override {
            return pt;
        }
        match &self.surface {
            Surface3::Sphere(s) => {
                // Prefer 3D boundary centroid: UV centroids can be unstable on
                // periodic seams and polar regions, which misclassifies caps.
                let surface_pt = if !self.boundary.is_empty() {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                } else if let Some(uv) = self.uv_centroid {
                    s.point_at(uv.x, uv.y)
                } else {
                    s.center + s.radius * DVec3::X
                };
                // Offset inward toward sphere center
                let to_center = (s.center - surface_pt).normalize_or_zero();
                let inward = if to_center.length_squared() > 0.5 {
                    to_center
                } else {
                    -self.normal
                };
                surface_pt + inward * (TOLERANCE_ABS * 10.0)
            }
            Surface3::Cylinder(c) => {
                // For cylinder faces, the outward normal points AWAY from the axis.
                // To get a sample point just inside the solid, offset toward the axis.
                let centroid = if self.boundary.is_empty() {
                    DVec3::ZERO
                } else {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                };
                // Compute inward direction (toward cylinder axis)
                let axis = c.axis.normalize();
                let to_axis = c.origin + axis * (centroid - c.origin).dot(axis) - centroid;
                let inward = to_axis.normalize_or_zero();
                // Use inward offset so the sample is just inside the cylinder surface
                centroid + inward * (TOLERANCE_ABS * 10.0)
            }
            Surface3::Torus(t) => {
                use rcad_kernel::geom::SurfaceEval;
                // Use UV centroid for a precise point on the torus surface.
                let surface_pt = if let Some(uv) = self.uv_centroid {
                    t.point_at(uv.x, uv.y)
                } else if !self.boundary.is_empty() {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                } else {
                    t.center + (t.major_radius + t.minor_radius) * DVec3::X
                };
                // Offset toward the tube center so the sample is inside regardless of face normal orientation.
                let axis = t.axis.normalize_or_zero();
                let local = surface_pt - t.center;
                let axial = local.dot(axis);
                let radial = local - axial * axis;
                let inward = if radial.length_squared() > 1e-18 {
                    let tube_center = t.center + axial * axis + radial.normalize() * t.major_radius;
                    (tube_center - surface_pt).normalize_or_zero()
                } else {
                    -self.normal
                };
                surface_pt + inward * (TOLERANCE_ABS * 10.0)
            }
            Surface3::Cone(c) => {
                use rcad_kernel::geom::SurfaceEval;
                // Use UV centroid for a precise point on the cone surface.
                let surface_pt = if let Some(uv) = self.uv_centroid {
                    c.point_at(uv.x, uv.y)
                } else if !self.boundary.is_empty() {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                } else {
                    c.point_at(0.0, 1.0)
                };
                // Offset toward the cone axis so the sample is inside regardless of face normal orientation.
                let axis = c.axis_dir();
                let local = surface_pt - c.apex;
                let axial = local.dot(axis);
                let axis_pt = c.apex + axis * axial;
                let inward = (axis_pt - surface_pt).normalize_or_zero();
                let inward = if inward.length_squared() > 0.5 {
                    inward
                } else {
                    -self.normal
                };
                surface_pt + inward * (TOLERANCE_ABS * 10.0)
            }
            _ => {
                let centroid = if self.boundary.is_empty() {
                    DVec3::ZERO
                } else {
                    self.boundary.iter().copied().sum::<DVec3>() / self.boundary.len() as f64
                };
                centroid + self.normal * TOLERANCE_ABS * 10.0
            }
        }
    }
}

type FaceEntry = (Vec<usize>, Vec<Vec<usize>>, Vec<[usize; 3]>, DVec3, Surface3, Option<[f64; 4]>);

#[derive(Debug, Clone)]
struct EdgeEntry {
    start: usize,
    end: usize,
    curve: Option<Curve3>,
    range: Option<[f64; 2]>,
}

/// Builds result BRep, deduplicating vertices and edges.
struct ResultBuilder {
    vertices: Vec<DVec3>,
    vertex_map: HashMap<u64, usize>, // hash of position → index
    edges: Vec<EdgeEntry>,
    faces: Vec<FaceEntry>, // (boundary vertex indices, triangles, normal, surface, uv_domain)
    face_origins: Vec<FaceOrigin>,
}

impl ResultBuilder {
    fn new() -> Self {
        Self {
            vertices: Vec::new(),
            vertex_map: HashMap::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            face_origins: Vec::new(),
        }
    }

    fn add_vertex(&mut self, point: DVec3) -> usize {
        let key = hash_point(point);
        if let Some(&idx) = self.vertex_map.get(&key) {
            // Double-check actual coincidence (hash collision protection)
            if points_coincide(self.vertices[idx], point) {
                return idx;
            }
        }
        // Linear scan fallback for hash collisions
        for (i, v) in self.vertices.iter().enumerate() {
            if points_coincide(*v, point) {
                return i;
            }
        }
        let idx = self.vertices.len();
        self.vertices.push(point);
        self.vertex_map.insert(key, idx);
        idx
    }

    fn add_edge(&mut self, v1: usize, v2: usize) -> usize {
        let key = (v1.min(v2), v1.max(v2));
        for (i, e) in self.edges.iter().enumerate() {
            if e.curve.is_none() && (e.start.min(e.end), e.start.max(e.end)) == key {
                return i;
            }
        }
        let idx = self.edges.len();
        self.edges.push(EdgeEntry {
            start: v1,
            end: v2,
            curve: None,
            range: None,
        });
        idx
    }

    fn add_analytic_circle_edge(&mut self, vertex: usize, circle: Circle3, range: [f64; 2]) -> usize {
        let idx = self.edges.len();
        self.edges.push(EdgeEntry {
            start: vertex,
            end: vertex,
            curve: Some(Curve3::Circle(circle)),
            range: Some(range),
        });
        idx
    }

    fn emit_face_with_origin(&mut self, sub: &SubFace, flip: bool, origin: FaceOrigin) {
        let normal = if flip { -sub.normal } else { sub.normal };

        if let Some(circle) = sub.analytic_outer_circle {
            let v_idx = self.add_vertex(Curve3::Circle(circle).point_at(0.0));
            let e_idx = self.add_analytic_circle_edge(v_idx, circle, [0.0, std::f64::consts::TAU]);
            self.faces.push((vec![e_idx], vec![], vec![], normal, sub.surface.clone(), sub.uv_domain));
            self.face_origins.push(origin);
            return;
        }

        // Add vertices for outer boundary
        let vert_indices: Vec<usize> = sub.boundary.iter().map(|&p| self.add_vertex(p)).collect();

        // Add edges for outer boundary
        let mut edge_indices = Vec::new();
        for i in 0..vert_indices.len() {
            let j = (i + 1) % vert_indices.len();
            let ei = self.add_edge(vert_indices[i], vert_indices[j]);
            edge_indices.push(ei);
        }

        // Triangulate outer boundary
        let mut tris = triangulate_polygon(&sub.boundary, normal);
        // Remap triangle indices from local (0..n) to result vertex indices
        for tri in &mut tris {
            for idx in tri.iter_mut() {
                *idx = vert_indices[*idx];
            }
        }

        // Handle inner wires (holes) — only create wire topology, NOT triangles.
        // The face triangulation covers only the outer boundary; inner wires are
        // stored as topological holes and will be tesselled separately if needed.
        let mut inner_wire_edges: Vec<Vec<usize>> = Vec::new();
        for wire_pts in &sub.inner_wires {
            if wire_pts.len() < 3 {
                continue;
            }
            // Add vertices for this inner wire
            let wire_verts: Vec<usize> = wire_pts.iter().map(|&p| self.add_vertex(p)).collect();
            // Add edges
            let mut wire_edges = Vec::new();
            for i in 0..wire_verts.len() {
                let j = (i + 1) % wire_verts.len();
                let ei = self.add_edge(wire_verts[i], wire_verts[j]);
                wire_edges.push(ei);
            }
            inner_wire_edges.push(wire_edges);
        }

        self.faces
            .push((edge_indices, inner_wire_edges, tris, normal, sub.surface.clone(), sub.uv_domain));
        self.face_origins.push(origin);
    }

    fn build(self) -> (BRep, BooleanHistory) {
        let ResultBuilder {
            vertices,
            vertex_map: _,
            edges: edge_entries,
            faces: face_entries,
            face_origins,
        } = self;

        let vertices = vertices
            .into_iter()
            .map(|point| Vertex { point })
            .collect();

        let edges: Vec<Edge> = edge_entries
            .iter()
            .into_iter()
            .map(|e| Edge {
                start: e.start,
                end: e.end,
            })
            .collect();

        let mut geom = rcad_kernel::GeomStore::default();
        geom.edge_curve = vec![None; edges.len()];
        geom.edge_curve_range = vec![None; edges.len()];
        for (ei, e) in edge_entries.iter().enumerate() {
            if let Some(curve) = &e.curve {
                let ci = geom.curves.len();
                geom.curves.push(curve.clone());
                geom.edge_curve[ei] = Some(ci);
                geom.edge_curve_range[ei] = e.range;
            }
        }
        let mut faces = Vec::new();

        for (edge_indices, inner_wire_edges, triangles, normal, surface, uv_domain) in face_entries {
            let wire = Wire {
                edges: edge_indices.iter().map(|&idx| WireEdge::fwd(idx)).collect(),
            };
            let inner_wires: Vec<Wire> = inner_wire_edges
                .into_iter()
                .map(|wire_edge_idxs| Wire {
                    edges: wire_edge_idxs.iter().map(|&idx| WireEdge::fwd(idx)).collect(),
                })
                .collect();
            let mesh_dirty = triangles.is_empty();
            faces.push(Face {
                outer_wire: wire,
                inner_wires,
                normal,
                triangles,
                mesh_dirty,
            });

            let surf_idx = geom.surfaces.len();
            geom.surfaces.push(surface);
            geom.face_surface.push(Some(surf_idx));
            geom.face_surface_range.push(uv_domain);
        }

        let history = BooleanHistory {
            face_origins,
            edge_origins: Vec::new(),
            vertex_origins: Vec::new(),
            shell_origins: Vec::new(),
            solid_origins: Vec::new(),
            tracker: HistoryTracker::new(),
            deleted_from_a: Vec::new(),
            deleted_from_b: Vec::new(),
            deletion_reasons: std::collections::HashMap::new(),
        };

        let brep = BRep {
            vertices,
            edges,
            solids: vec![Solid {
                shells: vec![Shell { faces }],
            }],
            geom,
            compound: None,
            compsolid: None,
        };
        (brep, history)
    }
}

fn hash_point(p: DVec3) -> u64 {
    // Quantize to tolerance grid for spatial hashing
    let scale = 1.0 / TOLERANCE_ABS;
    let ix = (p.x * scale).round() as i64;
    let iy = (p.y * scale).round() as i64;
    let iz = (p.z * scale).round() as i64;
    // FNV-1a style hash
    let mut h: u64 = 14695981039346656037;
    for v in [ix, iy, iz] {
        h ^= v as u64;
        h = h.wrapping_mul(1099511628211);
    }
    h
}

/// Annotate a `BooleanHistory` with per-edge and per-vertex origins by
/// matching result BRep positions against the DS vertex/edge pool.
///
/// Both `edge_origins` and `vertex_origins` are filled in-place.
fn annotate_history_from_ds(brep: &BRep, history: &mut BooleanHistory, ds: &DS) {
    // --- vertex origins ---
    let n_result_verts = brep.vertices.len();
    let mut vertex_origins: Vec<VertexOrigin> = Vec::with_capacity(n_result_verts);
    // ds[0..a_vertex_count] = A vertices, ds[a_vertex_count..total] = B vertices,
    // intersection vertices were added later (index >= a_vertex_count + b_vertex_count).
    let a_vc = ds.a_vertex_count;
    // Map result vertex index → DS vertex index (or usize::MAX if no match).
    let mut result_to_ds: Vec<usize> = vec![usize::MAX; n_result_verts];

    for (ri, rv) in brep.vertices.iter().enumerate() {
        let pt = rv.point;
        let mut best: Option<usize> = None;
        for (di, dv) in ds.vertices.iter().enumerate() {
            if (dv.point - pt).length_squared() < TOLERANCE_ABS * TOLERANCE_ABS * 4.0 {
                best = Some(di);
                break;
            }
        }
        result_to_ds[ri] = best.unwrap_or(usize::MAX);
        let origin = match best {
            Some(di) if di < a_vc => VertexOrigin::FromA(di),
            Some(di) => VertexOrigin::FromB(di - a_vc),
            None => VertexOrigin::Intersection,
        };
        vertex_origins.push(origin);
    }
    history.vertex_origins = vertex_origins;

    // --- edge origins ---
    let n_result_edges = brep.edges.len();
    let mut edge_origins: Vec<EdgeOrigin> = Vec::with_capacity(n_result_edges);
    let a_ec = ds.a_edge_count;
    let total_ds_edges = ds.edges.len();

    for re in &brep.edges {
        let ds_s = result_to_ds[re.start];
        let ds_e = result_to_ds[re.end];

        let origin = if ds_s == usize::MAX || ds_e == usize::MAX {
            EdgeOrigin::Generated
        } else if ds_s < a_vc && ds_e < a_vc {
            // Both endpoints are A vertices — look for a DS edge in A range.
            let found = (0..a_ec.min(total_ds_edges)).find(|&dei| {
                let de = &ds.edges[dei];
                (de.start_vertex == ds_s && de.end_vertex == ds_e)
                    || (de.start_vertex == ds_e && de.end_vertex == ds_s)
            });
            match found {
                Some(dei) => EdgeOrigin::FromA(dei),
                None => EdgeOrigin::SplitFromA(ds_s.min(a_vc - 1)),
            }
        } else if ds_s >= a_vc && ds_e >= a_vc {
            // Both endpoints are B vertices — look for a DS edge in B range.
            let found = (a_ec..total_ds_edges).find(|&dei| {
                let de = &ds.edges[dei];
                (de.start_vertex == ds_s && de.end_vertex == ds_e)
                    || (de.start_vertex == ds_e && de.end_vertex == ds_s)
            });
            match found {
                Some(dei) => EdgeOrigin::FromB(dei - a_ec),
                None => EdgeOrigin::SplitFromB(ds_s.min(ds.vertices.len().saturating_sub(1)) - a_vc),
            }
        } else {
            EdgeOrigin::Generated
        };
        edge_origins.push(origin);
    }
    history.edge_origins = edge_origins;
}

/// Populate result edge 3D geometry by matching result edge endpoints against
/// DS intersection curves.
///
/// This restores analytic curves (e.g. circles on boolean boundaries) so STEP
/// export does not fall back to straight lines when edge_curve is missing.
fn populate_result_edge_geometry_from_ds(brep: &mut BRep, ds: &DS) {
    if brep.geom.edge_curve.len() != brep.edges.len() {
        brep.geom.edge_curve = vec![None; brep.edges.len()];
    }
    if brep.geom.edge_curve_range.len() != brep.edges.len() {
        brep.geom.edge_curve_range = vec![None; brep.edges.len()];
    }

    let debug_edge_geom = std::env::var_os("RMSH_BOOL_DEBUG_EDGE_GEOM").is_some();
    let mut assigned_total = 0usize;
    let mut from_ic = 0usize;
    let proj_tol = 1e-3;

    for (ei, edge) in brep.edges.iter().enumerate() {
        if brep.geom.edge_curve[ei].is_some() {
            continue;
        }

        let p0 = brep.vertices[edge.start].point;
        let p1 = brep.vertices[edge.end].point;
        let pm = 0.5 * (p0 + p1);

        let mut matched: Option<(Curve3, [f64; 2], f64)> = None;
        for ic in &ds.intersection_curves {
            let cp0 = closest_point_on_curve(&ic.curve, p0, 64);
            let cp1 = closest_point_on_curve(&ic.curve, p1, 64);
            let cpm = closest_point_on_curve(&ic.curve, pm, 64);

            if cp0.distance > proj_tol || cp1.distance > proj_tol || cpm.distance > proj_tol {
                continue;
            }

            // Score combines endpoint and midpoint residuals to prefer the best curve
            // when multiple intersection curves are spatially nearby.
            let score = cp0.distance + cp1.distance + cpm.distance;
            let candidate = (ic.curve.clone(), [cp0.param, cp1.param], score);
            if matched.as_ref().map(|(_, _, s)| score < *s).unwrap_or(true) {
                matched = Some(candidate);
            }
        }

        if let Some((curve, range, _score)) = matched {
            let ci = brep.geom.curves.len();
            brep.geom.curves.push(curve);
            brep.geom.edge_curve[ei] = Some(ci);
            brep.geom.edge_curve_range[ei] = Some(range);
            assigned_total += 1;
            from_ic += 1;
        }
    }

    if debug_edge_geom {
        println!(
            "[edge-geom] edges={} assigned_total={} from_ic={}",
            brep.edges.len(),
            assigned_total,
            from_ic
        );
    }
}

fn aggregate_face_region_origin(face_origins: &[FaceOrigin]) -> ShellOrigin {
    let mut has_a = false;
    let mut has_b = false;
    let mut has_generated = false;
    for origin in face_origins {
        match origin {
            FaceOrigin::FromA(_) => has_a = true,
            FaceOrigin::FromB(_) => has_b = true,
            FaceOrigin::Generated => has_generated = true,
        }
    }

    match (has_a, has_b, has_generated) {
        (true, false, false) => ShellOrigin::FromA,
        (false, true, false) => ShellOrigin::FromB,
        (false, false, true) => ShellOrigin::Generated,
        _ => ShellOrigin::Mixed,
    }
}

fn aggregate_shell_region_origin(shell_origins: &[ShellOrigin]) -> SolidOrigin {
    let mut has_a = false;
    let mut has_b = false;
    let mut has_generated = false;
    let mut has_mixed = false;
    for origin in shell_origins {
        match origin {
            ShellOrigin::FromA => has_a = true,
            ShellOrigin::FromB => has_b = true,
            ShellOrigin::Generated => has_generated = true,
            ShellOrigin::Mixed => has_mixed = true,
        }
    }

    if has_mixed {
        return SolidOrigin::Mixed;
    }

    match (has_a, has_b, has_generated) {
        (true, false, false) => SolidOrigin::FromA,
        (false, true, false) => SolidOrigin::FromB,
        (false, false, true) => SolidOrigin::Generated,
        _ => SolidOrigin::Mixed,
    }
}

fn annotate_shell_and_solid_history(brep: &BRep, history: &mut BooleanHistory) {
    let mut face_cursor = 0;
    let mut shell_origins = Vec::new();
    let mut solid_origins = Vec::with_capacity(brep.solids.len());

    for solid in &brep.solids {
        let solid_shell_start = shell_origins.len();
        for shell in &solid.shells {
            let shell_face_count = shell.faces.len();
            let shell_face_origins = history
                .face_origins
                .get(face_cursor..face_cursor + shell_face_count)
                .unwrap_or(&[]);
            shell_origins.push(aggregate_face_region_origin(shell_face_origins));
            face_cursor += shell_face_count;
        }
        solid_origins.push(aggregate_shell_region_origin(&shell_origins[solid_shell_start..]));
    }

    debug_assert_eq!(face_cursor, history.face_origins.len());
    history.shell_origins = shell_origins;
    history.solid_origins = solid_origins;
}

/// Boolean result builder (OCCT: BOPAlgo_BOP).
/// Tracks face splice origins and participates in `BooleanHistory`.
pub struct BooleanBuilder<'a> {
    ds: &'a DS,
    op: BooleanOpType,
    use_glue: bool,
    glue_tolerance: f64,
}

impl<'a> BooleanBuilder<'a> {
    pub fn new(ds: &'a DS, op: BooleanOpType) -> Self {
        Self {
            ds,
            op,
            use_glue: false,
            glue_tolerance: TOLERANCE_ABS,
        }
    }

    pub fn with_glue(mut self, enable: bool, tolerance: f64) -> Self {
        self.use_glue = enable;
        self.glue_tolerance = tolerance.max(TOLERANCE_ABS);
        self
    }

    fn pcurve_matches_face_surface(
        &self,
        pcurve: &rcad_kernel::geom::Curve2d,
        surface: &Surface3,
        ic: &IntersectionCurve,
    ) -> bool {
        let samples: Vec<DVec3> = if ic.polyline.len() >= 3 {
            let mid = ic.polyline.len() / 2;
            vec![ic.polyline[0], ic.polyline[mid], *ic.polyline.last().unwrap()]
        } else if ic.polyline.len() == 2 {
            vec![ic.polyline[0], ic.polyline[1]]
        } else {
            let [t0, t1] = ic.t_range;
            let tm = 0.5 * (t0 + t1);
            vec![ic.curve.point_at(t0), ic.curve.point_at(tm), ic.curve.point_at(t1)]
        };

        let params: Vec<f64> = match pcurve {
            rcad_kernel::geom::Curve2d::BSpline(_) => {
                if samples.len() <= 1 {
                    vec![0.0]
                } else {
                    (0..samples.len())
                        .map(|i| i as f64 / (samples.len() - 1) as f64)
                        .collect()
                }
            }
            _ => {
                let [t0, t1] = ic.t_range;
                if samples.len() <= 1 {
                    vec![t0]
                } else {
                    (0..samples.len())
                        .map(|i| t0 + (t1 - t0) * i as f64 / (samples.len() - 1) as f64)
                        .collect()
                }
            }
        };

        let mut max_err: f64 = 0.0;
        for (sample, t) in samples.iter().zip(params.iter().copied()) {
            let uv = pcurve.point_at(t);
            let lifted = surface.point_at(uv.x, uv.y);
            max_err = max_err.max((lifted - *sample).length());
        }

        max_err.is_finite() && max_err <= 1e-3
    }

    pub fn build(&self) -> Result<BRep, BooleanError> {
        let (brep, _) = self.build_with_history()?;
        Ok(brep)
    }

    pub fn build_with_history(&self) -> Result<(BRep, BooleanHistory), BooleanError> {
        let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);
        let debug_classify = std::env::var_os("RMSH_BOOL_DEBUG_CLASSIFY").is_some();
        let mut b_sphere_stats: (usize, usize, usize, usize) = (0, 0, 0, 0);

        if a_faces.is_empty() || b_faces.is_empty() {
            return Err(BooleanError::EmptyInput);
        }

        let mut result = ResultBuilder::new();

        // Process A faces against B solid
        for &fi in &a_faces {
            let sub_faces = self.split_face(fi);
            for sub in sub_faces.iter() {
                let sample = sub.sample_point();
                let class = classify_point(sample, &b_faces, self.ds);

                let keep = match self.op {
                    BooleanOpType::Union => {
                        let glued_on =
                            self.use_glue && class == Classification::On && self.is_glued_face(fi, &b_faces);
                        class == Classification::Out || (class == Classification::On && !glued_on)
                    }
                    BooleanOpType::Intersection => {
                        class == Classification::In || class == Classification::On
                    }
                    BooleanOpType::Difference => class == Classification::Out,
                };

                if keep {
                    result.emit_face_with_origin(sub, false, FaceOrigin::FromA(fi));
                }
            }
        }

        // Process B faces against A solid
        for &fi in &b_faces {
            let is_sphere = matches!(self.ds.faces[fi].surface, Surface3::Sphere(_));
            if debug_classify && is_sphere {
                println!(
                    "[bool-classify] side=B face={} sphere curves_in={} boundary_verts={}",
                    fi,
                    self.ds.faces[fi].face_info.curves_in.len(),
                    self.ds.faces[fi].boundary_verts.len()
                );
            }
            let sub_faces = self.split_face(fi);
            for sub in sub_faces.iter() {
                let sample = sub.sample_point();
                let class = classify_point(sample, &a_faces, self.ds);

                if is_sphere {
                    match class {
                        Classification::In => b_sphere_stats.0 += 1,
                        Classification::Out => b_sphere_stats.1 += 1,
                        Classification::On => b_sphere_stats.2 += 1,
                    }
                }

                let keep = match self.op {
                    BooleanOpType::Union => class == Classification::Out,
                    BooleanOpType::Intersection => {
                        class == Classification::In || class == Classification::On
                    }
                    BooleanOpType::Difference => class == Classification::In,
                };

                if keep {
                    let flip = self.op == BooleanOpType::Difference;
                    result.emit_face_with_origin(sub, flip, FaceOrigin::FromB(fi));
                    if is_sphere {
                        b_sphere_stats.3 += 1;
                    }
                }

                if debug_classify && is_sphere {
                    println!(
                        "[bool-classify] op={:?} side=B face={} sphere_sub class={:?} keep={}",
                        self.op,
                        fi,
                        class,
                        keep
                    );
                }
            }
        }

        if debug_classify {
            println!(
                "[bool-classify] side=B sphere_subfaces: in={} out={} on={} kept={}",
                b_sphere_stats.0,
                b_sphere_stats.1,
                b_sphere_stats.2,
                b_sphere_stats.3
            );
        }

        let (mut brep, mut history) = result.build();
        if brep.solids[0].shells[0].faces.is_empty() {
            return Err(BooleanError::DegenerateResult);
        }

        // Annotate edge/vertex origins from the DS and aggregate shell/solid provenance.
        annotate_history_from_ds(&brep, &mut history, self.ds);
        annotate_shell_and_solid_history(&brep, &mut history);
        populate_result_edge_geometry_from_ds(&mut brep, self.ds);
        crate::geom_populate::populate_boolean_result_pcurves(&mut brep);

        // Debug-mode geometry integrity check.
        // Verifies that every face in the result has a non-zero normal vector.
        // This catches the most common class of geometry regression (degenerate faces
        // produced by a wrong normal computation) without requiring a full wire-closure
        // check (which the current builder doesn't yet guarantee for all curve types).
        #[cfg(debug_assertions)]
        for (fi, face) in brep.solids[0].shells[0].faces.iter().enumerate() {
            debug_assert!(
                face.normal != glam::DVec3::ZERO,
                "boolean_op result face {fi} has zero normal"
            );
        }

        Ok((brep, history))
    }

    /// Parallel version of `build_with_history`.
    ///
    /// Uses Rayon to process faces in parallel. Each face is split and classified
    /// independently, then results are merged. This can provide significant
    /// speedup for models with many faces (e.g., > 100 faces).
    ///
    /// # Performance
    ///
    /// - Small models (< 20 faces): May be slower due to thread overhead
    /// - Large models (> 100 faces): Typically 2-4x faster on multi-core systems
    pub fn build_with_history_par(&self) -> Result<(BRep, BooleanHistory), BooleanError> {
        let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);

        if a_faces.is_empty() || b_faces.is_empty() {
            return Err(BooleanError::EmptyInput);
        }

        // Fall back to sequential for small models to avoid thread overhead.
        const PAR_THRESHOLD: usize = 20;
        if a_faces.len() + b_faces.len() < PAR_THRESHOLD {
            return self.build_with_history();
        }

        // Process A faces in parallel
        let a_results: Vec<_> = a_faces
            .par_iter()
            .flat_map(|&fi| {
                let sub_faces = self.split_face(fi);
                sub_faces
                    .into_iter()
                    .filter_map(|sub| {
                        let sample = sub.sample_point();
                        let class = classify_point(sample, &b_faces, self.ds);

                        let keep = match self.op {
                            BooleanOpType::Union => {
                                class == Classification::Out || class == Classification::On
                            }
                            BooleanOpType::Intersection => {
                                class == Classification::In || class == Classification::On
                            }
                            BooleanOpType::Difference => class == Classification::Out,
                        };

                        if keep {
                            Some((sub, false, FaceOrigin::FromA(fi)))
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        // Process B faces in parallel
        let b_results: Vec<_> = b_faces
            .par_iter()
            .flat_map(|&fi| {
                let sub_faces = self.split_face(fi);
                sub_faces
                    .into_iter()
                    .filter_map(|sub| {
                        let sample = sub.sample_point();
                        let class = classify_point(sample, &a_faces, self.ds);

                        let keep = match self.op {
                            BooleanOpType::Union => class == Classification::Out,
                            BooleanOpType::Intersection => {
                                class == Classification::In || class == Classification::On
                            }
                            BooleanOpType::Difference => class == Classification::In,
                        };

                        if keep {
                            let flip = self.op == BooleanOpType::Difference;
                            Some((sub, flip, FaceOrigin::FromB(fi)))
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        // Merge results into ResultBuilder
        let mut result = ResultBuilder::new();
        for (sub, flip, origin) in a_results.into_iter().chain(b_results.into_iter()) {
            result.emit_face_with_origin(&sub, flip, origin);
        }

        let (mut brep, mut history) = result.build();
        if brep.solids[0].shells[0].faces.is_empty() {
            return Err(BooleanError::DegenerateResult);
        }

        annotate_history_from_ds(&brep, &mut history, self.ds);
        annotate_shell_and_solid_history(&brep, &mut history);
        populate_result_edge_geometry_from_ds(&mut brep, self.ds);
        crate::geom_populate::populate_boolean_result_pcurves(&mut brep);

        #[cfg(debug_assertions)]
        for (fi, face) in brep.solids[0].shells[0].faces.iter().enumerate() {
            debug_assert!(
                face.normal != glam::DVec3::ZERO,
                "boolean_op result face {fi} has zero normal"
            );
        }

        Ok((brep, history))
    }

    /// Split a face by intersection curves. If no intersection curves cross this
    /// face, returns the whole face as a single SubFace.
    fn split_face(&self, face_idx: usize) -> Vec<SubFace> {
        let face = &self.ds.faces[face_idx];
        let fi = &face.face_info;

        if fi.curves_in.is_empty() && matches!(face.surface, Surface3::Plane(_)) {
            // No intersections — return whole face
            let boundary = face
                .boundary_verts
                .iter()
                .map(|&vi| self.ds.vertices[vi].point)
                .collect();
            return vec![SubFace {
                boundary,
                surface: face.surface.clone(),
                normal: face.normal,
                uv_centroid: None,
                sample_override: None,
                uv_domain: None,
                inner_wires: vec![],
                analytic_outer_circle: None,
            }];
        }

        if std::env::var_os("RMSH_SKIP_CURVED_SPLIT").is_some()
            && !matches!(face.surface, Surface3::Plane(_))
        {
            let boundary = face
                .boundary_verts
                .iter()
                .map(|&vi| self.ds.vertices[vi].point)
                .collect();
            return vec![SubFace {
                boundary,
                surface: face.surface.clone(),
                normal: face.normal,
                uv_centroid: None,
                sample_override: None,
                uv_domain: None,
                inner_wires: vec![],
                analytic_outer_circle: None,
            }];
        }

        // For planar faces: project to 2D, split by intersection segments
        match &face.surface.clone() {
            Surface3::Plane(plane) => self.split_planar_face(face_idx, plane),
            Surface3::Cylinder(_)
            | Surface3::Sphere(_)
            | Surface3::Cone(_)
            | Surface3::Torus(_)
            | Surface3::BSpline(_)
            | Surface3::Bezier(_) => {
                if std::env::var_os("RMSH_USE_CURVED_LEGACY").is_some() {
                    self.split_curved_face_legacy(face_idx)
                } else {
                    self.split_curved_face_parametric(face_idx)
                }
            }
            _ => {
                // Other curved surfaces — return whole face for now
                let boundary = face
                    .boundary_verts
                    .iter()
                    .map(|&vi| self.ds.vertices[vi].point)
                    .collect();
                vec![SubFace {
                    boundary,
                    surface: face.surface.clone(),
                    normal: face.normal,
                    uv_centroid: None,
                    sample_override: None,
                    uv_domain: None,
                    inner_wires: vec![],
                    analytic_outer_circle: None,
                }]
            }
        }
    }

    /// Split a planar face by intersection line segments.
    ///
    /// Algorithm:
    /// 1. Project boundary + intersection segment endpoints to 2D
    /// 2. Find where intersection segment endpoints lie on boundary edges
    /// 3. Insert intersection points into boundary at correct positions
    /// 4. Walk augmented boundary to extract sub-polygons on each side
    fn split_planar_face(&self, face_idx: usize, plane: &Plane) -> Vec<SubFace> {
        let face = &self.ds.faces[face_idx];

        // Collect 3D boundary points
        let boundary_3d: Vec<DVec3> = face
            .boundary_verts
            .iter()
            .map(|&vi| self.ds.vertices[vi].point)
            .collect();

        if boundary_3d.len() < 3 {
            return vec![];
        }

        // Project boundary to 2D in the plane
        let (u_axis, v_axis) = plane_local_basis(plane);
        let project_to_2d = |p: DVec3| -> DVec2 {
            let d = p - plane.origin;
            DVec2::new(d.dot(u_axis), d.dot(v_axis))
        };
        let lift_to_3d = |uv: DVec2| -> DVec3 { plane.origin + u_axis * uv.x + v_axis * uv.y };

        let boundary_2d: Vec<DVec2> = boundary_3d.iter().map(|&p| project_to_2d(p)).collect();

        // Process each intersection curve to split the polygon
        let mut polygons_2d: Vec<Vec<DVec2>> = vec![boundary_2d];
        // Track circles that were embedded inside polygons (center_2d, radius).
        // When such a circle is fully inside a polygon, that polygon's centroid
        // may fall inside the circle — we must use a vertex-based sample instead.
        let mut embedded_circles: Vec<(DVec2, f64)> = Vec::new();
        // Circles fully embedded in the planar polygon are preserved as analytic
        // circular caps, emitted later as single closed circle edges.
        let mut analytic_circle_caps: Vec<(DVec2, f64)> = Vec::new();

        for &ci in &face.face_info.curves_in {
            let ic = &self.ds.intersection_curves[ci];

            let curve_halfspace_split: Option<Vec<Vec<DVec2>>> = match &ic.curve {
                Curve3::Circle(circle) => {
                    // Plane-sphere intersection produces a circle lying in the plane.
                    // Project the circle center to 2D and split by the circle boundary.
                    let center_2d = project_to_2d(circle.center);
                    let radius = circle.radius;
                    let mut next: Vec<Vec<DVec2>> = Vec::new();
                    for poly in &polygons_2d {
                        if circle_fully_inside_polygon(poly, center_2d, radius) {
                            // Keep the host polygon as the outside region and add an
                            // analytic inside circular cap.
                            next.push(poly.clone());
                            analytic_circle_caps.push((center_2d, radius));
                        } else {
                            let halves = split_polygon_by_circle_2d(poly, center_2d, radius);
                            next.extend(halves);
                        }
                    }
                    // Track this circle so we can compute correct sample points later
                    embedded_circles.push((center_2d, radius));
                    Some(next)
                }
                Curve3::Line(line) => {
                    // Use segment from start to end vertex
                    let p_start = self.ds.vertices[ic.start_vertex].point;
                    let p_end = self.ds.vertices[ic.end_vertex].point;
                    if points_coincide(p_start, p_end) {
                        None
                    } else {
                        let seg_s2d = project_to_2d(p_start);
                        let _seg_e2d = project_to_2d(p_end);
                        let mut next: Vec<Vec<DVec2>> = Vec::new();
                        for poly in &polygons_2d {
                            // Use line direction to split
                            let dir = DVec2::new(
                                (line.direction - plane.normal * line.direction.dot(plane.normal))
                                    .dot(u_axis),
                                (line.direction - plane.normal * line.direction.dot(plane.normal))
                                    .dot(v_axis),
                            );
                            let halves = split_polygon_2d_by_line(poly, seg_s2d, dir);
                            next.extend(halves);
                        }
                        Some(next)
                    }
                }
                _ => {
                    // For other curves, fall back to segment approach
                    let p_start = self.ds.vertices[ic.start_vertex].point;
                    let p_end = self.ds.vertices[ic.end_vertex].point;
                    if !points_coincide(p_start, p_end) {
                        let seg_s2d = project_to_2d(p_start);
                        let seg_e2d = project_to_2d(p_end);
                        let mut next: Vec<Vec<DVec2>> = Vec::new();
                        for poly in &polygons_2d {
                            let halves = split_polygon_2d_by_segment(poly, seg_s2d, seg_e2d);
                            next.extend(halves);
                        }
                        Some(next)
                    } else {
                        None
                    }
                }
            };

            if let Some(new_polys) = curve_halfspace_split
                && !new_polys.is_empty()
            {
                polygons_2d = new_polys;
            }
        }

        let mut out: Vec<SubFace> = polygons_2d
            .into_iter()
            .filter(|p| p.len() >= 3)
            .map(|poly_2d| {
                let boundary: Vec<DVec3> = poly_2d.iter().map(|&uv| lift_to_3d(uv)).collect();
                // If there are embedded circles and this polygon's centroid falls inside
                // one of them, use the first boundary vertex (offset by normal) as the
                // sample point instead. All polygon vertices of the outer region are
                // outside all embedded circles, so the first vertex is a valid sample.
                let sample_override = if !embedded_circles.is_empty() {
                    let centroid_2d = {
                        let sum = poly_2d.iter().fold(DVec2::ZERO, |acc, &p| acc + p);
                        sum / poly_2d.len() as f64
                    };
                    let centroid_in_circle = embedded_circles.iter().any(|&(c, r)| {
                        (centroid_2d - c).length() < r
                    });
                    if centroid_in_circle && !boundary.is_empty() {
                        // Pick first vertex (outside the circle) + normal offset
                        Some(boundary[0] + face.normal * crate::tolerance::TOLERANCE_ABS * 10.0)
                    } else {
                        None
                    }
                } else {
                    None
                };
                SubFace {
                    boundary,
                    surface: face.surface.clone(),
                    normal: face.normal,
                    uv_centroid: None,
                    sample_override,
                    uv_domain: None,
                    inner_wires: vec![],
                    analytic_outer_circle: None,
                }
            })
            .collect();

        // Add analytic circular caps as dedicated sub-faces so they can be
        // emitted as single closed CIRCLE edges in the result topology.
        for (center_2d, radius) in analytic_circle_caps {
            let enable_circle_collapse = std::env::var_os("RMSH_ENABLE_CIRCLE_COLLAPSE").is_some();
            let center_3d = lift_to_3d(center_2d);
            let circle = Circle3 {
                center: center_3d,
                normal: plane.normal.normalize_or_zero(),
                radius,
            };
            let boundary: Vec<DVec3> = (0..12)
                .map(|i| {
                    let t = std::f64::consts::TAU * i as f64 / 12.0;
                    Curve3::Circle(circle).point_at(t)
                })
                .collect();
            out.push(SubFace {
                boundary,
                surface: face.surface.clone(),
                normal: face.normal,
                uv_centroid: None,
                sample_override: Some(center_3d + face.normal * crate::tolerance::TOLERANCE_ABS * 10.0),
                uv_domain: None,
                inner_wires: vec![],
                analytic_outer_circle: if enable_circle_collapse { Some(circle) } else { None },
            });
        }

        out
    }

    fn faces_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .faces
            .iter()
            .enumerate()
            .filter(|(_, f)| f.origin == origin)
            .map(|(i, _)| i)
            .collect()
    }

    fn is_glued_face(&self, fi: usize, others: &[usize]) -> bool {
        others
            .iter()
            .any(|&fj| self.faces_form_glued_pair(fi, fj))
    }

    fn faces_form_glued_pair(&self, f1: usize, f2: usize) -> bool {
        let a = &self.ds.faces[f1];
        let b = &self.ds.faces[f2];
        if a.origin == b.origin {
            return false;
        }
        if !self.surfaces_glue_compatible(&a.surface, &b.surface) {
            return false;
        }
        let na_len2 = a.normal.length_squared();
        let nb_len2 = b.normal.length_squared();
        if na_len2 <= TOLERANCE_ABS || nb_len2 <= TOLERANCE_ABS {
            return false;
        }
        let na = a.normal / na_len2.sqrt();
        let nb = b.normal / nb_len2.sqrt();
        if na.dot(nb) > -0.99 {
            return false;
        }
        self.boundaries_fully_overlap(f1, f2)
    }

    fn surfaces_glue_compatible(&self, s1: &Surface3, s2: &Surface3) -> bool {
        let tol = self.glue_tolerance;
        let axis_parallel = |a: DVec3, b: DVec3| {
            let la = a.length();
            let lb = b.length();
            if la <= TOLERANCE_ABS || lb <= TOLERANCE_ABS {
                return false;
            }
            (a / la).dot(b / lb).abs() >= 0.999
        };

        match (s1, s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                if !axis_parallel(p1.normal, p2.normal) {
                    return false;
                }
                let n = p1.normal.normalize_or_zero();
                (p2.origin - p1.origin).dot(n).abs() <= tol * 2.0
            }
            (Surface3::Sphere(s1), Surface3::Sphere(s2)) => {
                (s1.center - s2.center).length() <= tol * 2.0
                    && (s1.radius - s2.radius).abs() <= tol
            }
            (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
                if !axis_parallel(c1.axis, c2.axis) {
                    return false;
                }
                let axis = c1.axis.normalize_or_zero();
                (c2.origin - c1.origin).cross(axis).length() <= tol * 2.0
                    && (c1.radius - c2.radius).abs() <= tol
            }
            (Surface3::Cone(c1), Surface3::Cone(c2)) => {
                axis_parallel(c1.axis, c2.axis)
                    && (c1.apex_point() - c2.apex_point()).length() <= tol * 2.0
                    && (c1.half_angle_rad - c2.half_angle_rad).abs() <= tol
            }
            (Surface3::Torus(t1), Surface3::Torus(t2)) => {
                axis_parallel(t1.axis, t2.axis)
                    && (t1.center - t2.center).length() <= tol * 2.0
                    && (t1.major_radius - t2.major_radius).abs() <= tol
                    && (t1.minor_radius - t2.minor_radius).abs() <= tol
            }
            _ => false,
        }
    }

    fn boundaries_fully_overlap(&self, f1: usize, f2: usize) -> bool {
        let pts1 = self.ds.face_boundary_points(f1);
        let pts2 = self.ds.face_boundary_points(f2);
        if pts1.len() < 3 || pts2.len() < 3 || pts1.len() != pts2.len() {
            return false;
        }
        let tol = self.glue_tolerance;
        let mut used = vec![false; pts2.len()];
        for p1 in &pts1 {
            let mut found = false;
            for (j, p2) in pts2.iter().enumerate() {
                if used[j] {
                    continue;
                }
                if (*p1 - *p2).length() <= tol {
                    used[j] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                return false;
            }
        }
        true
    }

    /// Fast check for potential glued face pairs using bounding box pre-filter.
    ///
    /// This optimization reduces the number of full boundary comparisons by
    /// first checking if face bounding boxes overlap.
    fn fast_glue_candidate_check(&self, f1: usize, f2: usize) -> bool {
        let a = &self.ds.faces[f1];
        let b = &self.ds.faces[f2];

        // Quick origin check
        if a.origin == b.origin {
            return false;
        }

        // Quick normal check (must be anti-parallel for glue)
        let na_len2 = a.normal.length_squared();
        let nb_len2 = b.normal.length_squared();
        if na_len2 <= TOLERANCE_ABS || nb_len2 <= TOLERANCE_ABS {
            return false;
        }
        let na = a.normal / na_len2.sqrt();
        let nb = b.normal / nb_len2.sqrt();
        if na.dot(nb) > -0.95 {
            return false;
        }

        // Bounding box overlap check
        let pts1 = self.ds.face_boundary_points(f1);
        let pts2 = self.ds.face_boundary_points(f2);

        if pts1.is_empty() || pts2.is_empty() {
            return false;
        }

        // Compute bounding boxes
        let mut min1 = pts1[0];
        let mut max1 = pts1[0];
        for p in &pts1[1..] {
            min1 = min1.min(*p);
            max1 = max1.max(*p);
        }

        let mut min2 = pts2[0];
        let mut max2 = pts2[0];
        for p in &pts2[1..] {
            min2 = min2.min(*p);
            max2 = max2.max(*p);
        }

        // Check for bounding box overlap with tolerance margin
        let tol = self.glue_tolerance;
        let overlap = min1.x - tol <= max2.x && max1.x + tol >= min2.x
            && min1.y - tol <= max2.y && max1.y + tol >= min2.y
            && min1.z - tol <= max2.z && max1.z + tol >= min2.z;

        overlap
    }

    /// Detect all glued face pairs using optimized algorithm.
    ///
    /// This function uses bounding box pre-filtering to reduce the number
    /// of expensive boundary comparisons.
    fn detect_all_glued_pairs(&self, a_faces: &[usize], b_faces: &[usize]) -> Vec<(usize, usize)> {
        let mut pairs: Vec<(usize, usize)> = Vec::new();

        for &fi in a_faces {
            for &fj in b_faces {
                // Fast pre-filter
                if !self.fast_glue_candidate_check(fi, fj) {
                    continue;
                }

                // Full compatibility check
                if self.faces_form_glued_pair(fi, fj) {
                    pairs.push((fi, fj));
                }
            }
        }

        pairs
    }

    /// Build glued pairs information for fast path processing.
    ///
    /// Returns a map from face index to its glued counterpart.
    fn build_glue_map(&self, a_faces: &[usize], b_faces: &[usize]) -> HashMap<usize, usize> {
        let pairs = self.detect_all_glued_pairs(a_faces, b_faces);
        let mut glue_map: HashMap<usize, usize> = HashMap::new();

        for (fi, fj) in pairs {
            glue_map.insert(fi, fj);
            glue_map.insert(fj, fi);
        }

        glue_map
    }

    /// Split a curved face (Cylinder, Sphere, Cone, Torus) by intersection polylines.
    ///
    /// Legacy approximate method: for each intersection polyline that crosses the face,
    /// we split the boundary point list into two halves at the points closest to the
    /// polyline endpoints. Kept as fallback when UV data or PCurves are unavailable.
    fn split_curved_face_legacy(&self, face_idx: usize) -> Vec<SubFace> {
        let face = &self.ds.faces[face_idx];
        let surface = face.surface.clone();
        let normal = face.normal;
        let face_curve_indices = self.curve_indices_for_face(face_idx);

        // Collect all intersection polylines for this face
        let mut all_polylines: Vec<Vec<DVec3>> = Vec::new();
        for &ci in &face_curve_indices {
            let ic = &self.ds.intersection_curves[ci];
            if ic.polyline.len() >= 2 {
                all_polylines.push(ic.polyline.clone());
            } else {
                // Analytic curve — sample it into a polyline (e.g. circle)
                let pts: Vec<DVec3> = (0..=16)
                    .map(|i| {
                        let t = ic.t_range[0] + (ic.t_range[1] - ic.t_range[0]) * i as f64 / 16.0;
                        use rcad_kernel::CurveEval;
                        ic.curve.point_at(t)
                    })
                    .collect();
                all_polylines.push(pts);
            }
        }

        if all_polylines.is_empty() {
            let boundary = face
                .boundary_verts
                .iter()
                .map(|&vi| self.ds.vertices[vi].point)
                .collect();
            return vec![SubFace {
                boundary,
                surface,
                normal,
                uv_centroid: None,
                sample_override: None,
                uv_domain: None,
                inner_wires: vec![],
                analytic_outer_circle: None,
            }];
        }

        // Collect boundary vertices
        let boundary_pts: Vec<DVec3> = face
            .boundary_verts
            .iter()
            .map(|&vi| self.ds.vertices[vi].point)
            .collect();

        if boundary_pts.len() < 3 {
            return vec![SubFace {
                boundary: boundary_pts,
                surface,
                normal,
                uv_centroid: None,
                sample_override: None,
                uv_domain: None,
                inner_wires: vec![],
                analytic_outer_circle: None,
            }];
        }

        // For each intersection polyline, split the boundary into two sub-faces
        // by finding the boundary points closest to each polyline endpoint.
        let mut result_boundaries: Vec<Vec<DVec3>> = vec![boundary_pts];

        for polyline in &all_polylines {
            let (Some(&seg_start), Some(&seg_end)) = (polyline.first(), polyline.last()) else {
                continue;
            };

            let mut next_result: Vec<Vec<DVec3>> = Vec::new();
            for bnd in result_boundaries.drain(..) {
                let n = bnd.len();
                if n < 3 {
                    next_result.push(bnd);
                    continue;
                }

                // Find indices of boundary points closest to the two polyline endpoints
                let Some((i_start, _)) = bnd
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        a.distance_squared(seg_start)
                            .partial_cmp(&b.distance_squared(seg_start))
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                else {
                    next_result.push(bnd);
                    continue;
                };
                let Some((i_end, _)) = bnd
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        a.distance_squared(seg_end)
                            .partial_cmp(&b.distance_squared(seg_end))
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                else {
                    next_result.push(bnd);
                    continue;
                };

                // Ensure i_start < i_end for consistent splitting
                let (ia, ib, p_a, p_b) = if i_start <= i_end {
                    (i_start, i_end, seg_start, seg_end)
                } else {
                    (i_end, i_start, seg_end, seg_start)
                };

                if ia == ib {
                    // Degenerate: can't split, keep as is
                    next_result.push(bnd);
                    continue;
                }

                // Sub-face A: bnd[0..=ia] + polyline + bnd[ib..=n-1]
                let mut sub_a: Vec<DVec3> = bnd[..=ia].to_vec();
                sub_a.push(p_a);
                for &p in polyline.iter().skip(1).rev().skip(1) {
                    sub_a.push(p);
                }
                sub_a.push(p_b);
                sub_a.extend_from_slice(&bnd[ib..]);

                // Sub-face B: bnd[ia..=ib] + reverse polyline
                let mut sub_b: Vec<DVec3> = bnd[ia..=ib].to_vec();
                sub_b.push(p_b);
                for &p in polyline.iter().skip(1).rev().skip(1) {
                    sub_b.push(p);
                }
                sub_b.push(p_a);

                if sub_a.len() >= 3 {
                    next_result.push(sub_a);
                }
                if sub_b.len() >= 3 {
                    next_result.push(sub_b);
                }
            }
            result_boundaries = next_result;
        }

        result_boundaries
            .into_iter()
            .filter(|b| b.len() >= 3)
            .map(|boundary| SubFace {
                boundary,
                surface: surface.clone(),
                normal,
                uv_centroid: None,
                sample_override: None,
                uv_domain: None,
                inner_wires: vec![],
                analytic_outer_circle: None,
            })
            .collect()
    }

    /// Unwrap a UV polyline's U coordinate to remove seam jumps.
    /// For periodic surfaces (cylinder, cone, torus), consecutive points whose
    /// U values differ by more than π indicate a seam crossing; we accumulate
    /// offsets of ±period to make the polyline continuous in U.
    fn unwrap_u_polyline(&self, pts: Vec<glam::DVec2>, period: f64) -> Vec<glam::DVec2> {
        if pts.len() < 2 {
            return pts;
        }
        let mut result = Vec::with_capacity(pts.len());
        result.push(pts[0]);
        let mut offset = 0.0_f64;
        for i in 1..pts.len() {
            let prev_u = result[i - 1].x;
            let curr_u = pts[i].x + offset;
            let diff = curr_u - prev_u;
            if diff > period * 0.5 {
                offset -= period;
            } else if diff < -period * 0.5 {
                offset += period;
            }
            result.push(glam::DVec2::new(pts[i].x + offset, pts[i].y));
        }
        result
    }

    /// Split a curved face using parameter-space (UV) 2D clipping.
    ///
    /// For each intersection curve on this face, samples the associated PCurve
    /// into a 2D trim polyline in UV space, then splits the UV boundary polygon.
    /// Maps resulting sub-polygons back to 3D via surface evaluation.
    ///
    /// Falls back to `split_curved_face_legacy` when UV data or PCurves are missing.
    fn split_curved_face_parametric(&self, face_idx: usize) -> Vec<SubFace> {

        let face = &self.ds.faces[face_idx];
        let debug_classify = std::env::var_os("RMSH_BOOL_DEBUG_CLASSIFY").is_some();
        let debug_sphere = debug_classify && matches!(face.surface, Surface3::Sphere(_));

        // Need UV boundary to operate in parameter space
        let uv_boundary = match &face.uv_boundary {
            Some(b) if b.len() >= 3 => b.clone(),
            _ => return self.split_curved_face_legacy(face_idx),
        };

        let surface = face.surface.clone();
        let normal = face.normal;

        // Collect 2D trim polylines from PCurves for each intersection curve
        let mut trim_polylines: Vec<Vec<DVec2>> = Vec::new();
        // Detect if this face is a periodic surface (cylinder, cone, torus) needing seam unwrap.
        let is_periodic_u = matches!(&surface,
            Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Torus(_)
        );
        // For sphere, u is also periodic in [-π, π].
        let is_sphere = matches!(&surface, Surface3::Sphere(_));
        let u_period = if is_periodic_u { std::f64::consts::TAU } else if is_sphere { std::f64::consts::TAU } else { 0.0 };

        let face_curve_indices = self.curve_indices_for_face(face_idx);
        if debug_sphere {
            println!(
                "[sphere-split] face={} uv_boundary_pts={} curves_in={}",
                face_idx,
                uv_boundary.len(),
                face_curve_indices.len()
            );
        }
        for &ci in &face_curve_indices {
            if let Some(pcurve) = self.find_pcurve_for_face(ci, face_idx) {
                let ic = &self.ds.intersection_curves[ci];
                let [t0, t1] = ic.t_range;
            const N: usize = 16;
                let raw_pts: Vec<DVec2> = match &pcurve {
                    // BSpline PCurves from polyline_pcurve_by_projection use
                    // chord-length parameterization normalized to [0,1].
                    // The 3D arc-length t_range is unrelated to the BSpline domain.
                    rcad_kernel::geom::Curve2d::BSpline(_) => (0..=N)
                        .map(|i| {
                            let t = i as f64 / N as f64;
                            pcurve.point_at(t)
                        })
                        .collect(),
                    // Analytic curves (Line2d, Circle2d, Ellipse2d) use the same
                    // t parameterization as the 3D intersection curve.
                    _ => (0..=N)
                        .map(|i| {
                            let t = t0 + (t1 - t0) * i as f64 / N as f64;
                            pcurve.point_at(t)
                        })
                        .collect(),
                };
                if raw_pts.len() < 2 {
                    continue;
                }

                // For periodic surfaces, unwrap the u-coordinate to remove seam jumps.
                // A jump > π in u between consecutive points indicates a seam crossing;
                // we add/subtract 2π to make the polyline continuous.
                let pts = if u_period > 0.0 {
                    self.unwrap_u_polyline(raw_pts, u_period)
                } else {
                    raw_pts
                };

                // If the unwrapped polyline spans more than 2π in u, the intersection
                // curve goes all the way around the surface — split at the seam instead
                // of trying to split the UV polygon with a polyline that exits and re-enters.
                if u_period > 0.0 && pts.len() >= 2 {
                    let u_span = pts.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max)
                        - pts.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                    // If span > π (the trim cuts across the seam) we need to clip to [0, 2π].
                    // Shift back into [0, 2π] by remapping each point mod 2π.
                    let pts = if u_span > std::f64::consts::PI {
                        // Re-centre: find the offset that brings the midpoint into [0, 2π].
                        let u_mid = pts.iter().map(|p| p.x).sum::<f64>() / pts.len() as f64;
                        let offset = (u_mid / u_period).floor() * u_period;
                        pts.into_iter().map(|p| DVec2::new(p.x - offset, p.y)).collect::<Vec<_>>()
                    } else {
                        pts
                    };
                    trim_polylines.push(pts);
                } else {
                    trim_polylines.push(pts);
                }
                if debug_sphere {
                    let trim = trim_polylines.last().unwrap();
                    let closed = if trim.len() >= 2 {
                        (trim[0] - trim[trim.len() - 1]).length_squared() < 1e-10
                    } else {
                        false
                    };
                    println!(
                        "[sphere-split] face={} curve={} trim_pts={} closed={}",
                        face_idx,
                        ci,
                        trim.len(),
                        closed
                    );
                }
            }
        }

        // If no PCurves available, fall back to legacy method
        if trim_polylines.is_empty() {
            return self.split_curved_face_legacy(face_idx);
        }

        // Split UV polygon by each trim polyline
        let mut uv_polygons: Vec<Vec<DVec2>> = vec![uv_boundary];

        for trim in &trim_polylines {
            let mut next: Vec<Vec<DVec2>> = Vec::new();
            for poly in uv_polygons.drain(..) {
                // Skip invalid polygons
                if !is_valid_uv_polygon(&poly) {
                    continue;
                }
                let effective_trim = if u_period > 0.0 {
                    periodic_trim_to_open_isoline(&poly, trim, u_period)
                        .unwrap_or_else(|| trim.clone())
                } else {
                    trim.clone()
                };
                let halves = split_uv_polygon_by_trim(&poly, &effective_trim);
                next.extend(halves);
            }
            if debug_sphere {
                println!(
                    "[sphere-split] face={} after_trim polys={}",
                    face_idx,
                    next.len()
                );
            }
            uv_polygons = next;
        }

        // Handle seam crossings for periodic surfaces
        if u_period > 0.0 {
            let seam_u = 0.0; // Standard seam at u=0 (or u=-π for sphere)
            uv_polygons = uv_polygons
                .into_iter()
                .flat_map(|poly| {
                    if is_valid_uv_polygon(&poly) {
                        handle_periodic_seam_crossing(&poly, u_period, seam_u)
                    } else {
                        vec![]
                    }
                })
                .collect();
            if debug_sphere {
                println!(
                    "[sphere-split] face={} after_seam polys={}",
                    face_idx,
                    uv_polygons.len()
                );
            }
        }

        // Map each UV sub-polygon back to 3D
        uv_polygons
            .into_iter()
            .filter(|p| p.len() >= 3 && is_valid_uv_polygon(p))
            .map(|uv_poly| {
                let n = uv_poly.len() as f64;
                let centroid_uv = uv_poly.iter().copied().sum::<DVec2>() / n;

                // Compute the UV bounding box of this sub-polygon so that
                // tessellate_curved_face samples only the correct sub-domain.
                let u_min = uv_poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                let u_max = uv_poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
                let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
                let v_max = uv_poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
                let uv_domain = if u_min.is_finite() && u_max.is_finite()
                    && v_min.is_finite() && v_max.is_finite()
                    && (u_max - u_min) > 1e-14 && (v_max - v_min) > 1e-14
                {
                    Some([u_min, u_max, v_min, v_max])
                } else {
                    None
                };

                let boundary: Vec<DVec3> = match &surface {
                    Surface3::Sphere(_) | Surface3::Cone(_) => {
                        // Use enhanced degenerate point handling for sphere poles and cone apex
                        let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
                        let v_max = uv_poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
                        let near_degenerate = match &surface {
                            Surface3::Sphere(_) => v_min < 0.01 || v_max > std::f64::consts::PI - 0.01,
                            Surface3::Cone(_) => v_min < 0.01,
                            _ => false,
                        };
                        if near_degenerate {
                            handle_degenerate_points(&uv_poly, &surface)
                        } else {
                            curved_subface_boundary_3d(&uv_poly, &trim_polylines, &surface)
                        }
                    }
                    _ => uv_poly.iter().map(|uv| surface.point_at(uv.x, uv.y)).collect(),
                };

                // Detect inner wires: trim polylines that are closed loops
                // fully contained within this UV polygon.
                let inner_wires: Vec<Vec<DVec3>> = trim_polylines
                    .iter()
                    .filter(|trim| {
                        if trim.len() < 3 {
                            return false;
                        }
                        // Check if closed (first and last point coincide)
                        let first = trim[0];
                        let last = trim[trim.len() - 1];
                        if (first - last).length_squared() > 1e-10 {
                            return false;
                        }
                        // Check if centroid is inside this UV polygon
                        let centroid = trim.iter().copied().sum::<DVec2>() / trim.len() as f64;
                        point_in_polygon_2d(&uv_poly, centroid)
                    })
                    .map(|trim| {
                        trim.iter()
                            .map(|uv| surface.point_at(uv.x, uv.y))
                            .collect()
                    })
                    .collect();

                // For curved surfaces, compute the actual surface normal at the centroid UV
                let sub_normal = {
                    let computed = surface.normal_at(centroid_uv.x, centroid_uv.y);
                    // If normal computation failed, fall back to face normal
                    if computed.length_squared() > 0.5 {
                        computed
                    } else {
                        normal
                    }
                };
                SubFace {
                    boundary,
                    surface: surface.clone(),
                    normal: sub_normal,
                    uv_centroid: Some(centroid_uv),
                    sample_override: None,
                    uv_domain,
                    inner_wires,
                    analytic_outer_circle: None,
                }
            })
            .collect()
    }

    /// Find the PCurve (2D parametric curve) for the given intersection curve
    /// as it lies on the given face. Searches FaceFace interferences to determine
    /// whether this face is f1 (use pcurve_on_a) or f2 (use pcurve_on_b).
    fn find_pcurve_for_face(
        &self,
        curve_idx: usize,
        face_idx: usize,
    ) -> Option<rcad_kernel::geom::Curve2d> {
        for interference in &self.ds.interferences {
            if let Interference::FaceFace { f1, f2, curves, .. } = interference
                && curves.contains(&curve_idx)
            {
                let ic = &self.ds.intersection_curves[curve_idx];
                if *f1 == face_idx {
                    return ic.pcurve_on_a.clone();
                } else if *f2 == face_idx {
                    return ic.pcurve_on_b.clone();
                }
            }
        }

        let ic = &self.ds.intersection_curves[curve_idx];
        let surface = &self.ds.faces[face_idx].surface;
        if let Some(pcurve) = &ic.pcurve_on_a
            && self.pcurve_matches_face_surface(pcurve, surface, ic)
        {
            return Some(pcurve.clone());
        }
        if let Some(pcurve) = &ic.pcurve_on_b
            && self.pcurve_matches_face_surface(pcurve, surface, ic)
        {
            return Some(pcurve.clone());
        }
        None
    }

    fn curve_indices_for_face(&self, face_idx: usize) -> Vec<usize> {
        let mut out: Vec<usize> = Vec::new();
        let mut seen: HashSet<usize> = HashSet::new();

        for &ci in &self.ds.faces[face_idx].face_info.curves_in {
            if seen.insert(ci) {
                out.push(ci);
            }
        }

        for interference in &self.ds.interferences {
            if let Interference::FaceFace { f1, f2, curves, .. } = interference
                && (*f1 == face_idx || *f2 == face_idx)
            {
                for &ci in curves {
                    if seen.insert(ci) {
                        out.push(ci);
                    }
                }
            }
        }

        if out.is_empty() {
            for ci in 0..self.ds.intersection_curves.len() {
                if self.find_pcurve_for_face(ci, face_idx).is_some() && seen.insert(ci) {
                    out.push(ci);
                }
            }
        }

        out
    }
}

/// Compute a robust 3D boundary for a curved sub-face given its UV polygon
/// and trim polylines.
///
/// Unlike `sphere_subface_boundary_3d` which only evaluates UV corners, this
/// function samples each UV edge into N points. This prevents degenerate
/// polygons when multiple corners collapse at a surface singularity (sphere
/// poles, cone apex).
///
/// Algorithm:
/// 1. Subdivide each UV edge into N samples, evaluate via surface.point_at
/// 2. Consecutive dedup: collapse runs of points near a singularity
/// 3. If < 3 points remain, supplement with trim polyline 3D points
/// 4. Global dedup, return
fn curved_subface_boundary_3d(
    uv_poly: &[DVec2],
    trim_polylines_uv: &[Vec<DVec2>],
    surface: &Surface3,
) -> Vec<DVec3> {
    const EDGE_SAMPLES: usize = 8;

    let mut pts: Vec<DVec3> = Vec::new();

    // 1. Sample each UV edge and evaluate 3D positions
    let n = uv_poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = uv_poly[i];
        let b = uv_poly[j];
        for k in 0..EDGE_SAMPLES {
            let t = k as f64 / EDGE_SAMPLES as f64;
            let uv = DVec2::new(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y));
            pts.push(surface.point_at(uv.x, uv.y));
        }
    }

    // 2. Consecutive deduplication — collapse runs of pole/apex samples
    let mut deduped: Vec<DVec3> = Vec::new();
    for p in &pts {
        if deduped.is_empty() || (*p - deduped[deduped.len() - 1]).length_squared() > TOLERANCE_ABS * TOLERANCE_ABS {
            deduped.push(*p);
        }
    }
    // Close the loop: remove last point if it equals the first
    if deduped.len() > 1 && (deduped[0] - deduped[deduped.len() - 1]).length_squared() < TOLERANCE_ABS * TOLERANCE_ABS {
        deduped.pop();
    }

    // 3. If still degenerate, supplement with trim polyline 3D points
    if deduped.len() < 3 {
        for trim_uv in trim_polylines_uv {
            if trim_uv.len() < 2 {
                continue;
            }
            for uv in trim_uv {
                let p3 = surface.point_at(uv.x, uv.y);
                if point_in_polygon_2d(uv_poly, *uv) || point_near_polygon_2d(uv_poly, *uv, 0.1) {
                    // Only add if not already in deduped
                    if deduped.iter().all(|q| (p3 - *q).length_squared() > TOLERANCE_ABS * TOLERANCE_ABS) {
                        deduped.push(p3);
                    }
                }
            }
        }
    }

    // 4. Final global dedup
    let mut result: Vec<DVec3> = Vec::new();
    for p in &deduped {
        if result.iter().all(|q| (*p - *q).length_squared() > TOLERANCE_ABS * TOLERANCE_ABS) {
            result.push(*p);
        }
    }

    result
}

/// Check if a 2D point is within `margin` of any edge of a polygon.
fn point_near_polygon_2d(poly: &[DVec2], pt: DVec2, margin: f64) -> bool {
    let n = poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = poly[i];
        let b = poly[j];
        let ab = b - a;
        let len_sq = ab.length_squared();
        let t = if len_sq < 1e-14 { 0.0 } else { ((pt - a).dot(ab) / len_sq).clamp(0.0, 1.0) };
        let closest = a + t * ab;
        if (pt - closest).length() < margin {
            return true;
        }
    }
    false
}

/// Detect and handle UV seam crossings for periodic surfaces.
/// Returns a list of split polygons if the UV polygon crosses the seam.
fn handle_periodic_seam_crossing(
    uv_poly: &[DVec2],
    u_period: f64,
    seam_u: f64,
) -> Vec<Vec<DVec2>> {
    let n = uv_poly.len();
    if n < 3 || u_period <= 0.0 {
        return vec![uv_poly.to_vec()];
    }

    // Find all edges that cross the seam
    let mut seam_crossings: Vec<(usize, f64, DVec2)> = Vec::new();

    for i in 0..n {
        let j = (i + 1) % n;
        let u_i = uv_poly[i].x;
        let u_j = uv_poly[j].x;

        // Check for seam crossing (jump > period/2)
        let du = u_j - u_i;
        if du.abs() > u_period * 0.5 {
            // Compute intersection point with seam
            let t = if du > 0.0 {
                (seam_u + u_period - u_i) / du
            } else {
                (seam_u - u_i) / du
            };

            if t > 0.0 && t < 1.0 {
                let v_i = uv_poly[i].y;
                let v_j = uv_poly[j].y;
                let seam_v = v_i + t * (v_j - v_i);
                let seam_pt = DVec2::new(seam_u, seam_v);
                seam_crossings.push((i, t, seam_pt));
            }
        }
    }

    // If no crossings or odd number of crossings (invalid), return original
    if seam_crossings.is_empty() || seam_crossings.len() % 2 != 0 {
        return vec![uv_poly.to_vec()];
    }

    // Sort crossings by edge index
    seam_crossings.sort_by_key(|&(idx, _, _)| idx);

    // For now, handle the simple case of exactly 2 crossings
    if seam_crossings.len() == 2 {
        let (idx1, _, pt1) = seam_crossings[0];
        let (idx2, _, pt2) = seam_crossings[1];

        // Build two sub-polygons
        let mut poly1: Vec<DVec2> = Vec::new();
        let mut poly2: Vec<DVec2> = Vec::new();

        // poly1: from crossing1 to crossing2 (wrapping the other way)
        poly1.push(pt1);
        for i in (idx1 + 1)..=idx2 {
            if i < n {
                poly1.push(uv_poly[i]);
            }
        }
        poly1.push(pt2);

        // poly2: from crossing2 back to crossing1
        poly2.push(pt2);
        for i in (idx2 + 1)..n {
            poly2.push(uv_poly[i]);
        }
        for i in 0..=idx1 {
            poly2.push(uv_poly[i]);
        }
        poly2.push(pt1);

        let mut result = Vec::new();
        if poly1.len() >= 3 {
            result.push(poly1);
        }
        if poly2.len() >= 3 {
            result.push(poly2);
        }

        if result.is_empty() {
            vec![uv_poly.to_vec()]
        } else {
            result
        }
    } else {
        // Multiple crossing pairs - complex case, return original for now
        vec![uv_poly.to_vec()]
    }
}

/// Detect degenerate points (poles, apex) and handle them in UV polygon.
/// Returns a modified 3D boundary that correctly handles surface singularities.
fn handle_degenerate_points(
    uv_poly: &[DVec2],
    surface: &Surface3,
) -> Vec<DVec3> {
    match surface {
        Surface3::Sphere(s) => {
            // Sphere has two poles at v=0 and v=π
            let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
            let v_max = uv_poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);

            let mut boundary_3d = Vec::new();
            let pole_tol = 0.01; // Tolerance for detecting near-pole

            // Check if polygon touches the north pole (v ≈ 0)
            let touches_north_pole = v_min < pole_tol;
            // Check if polygon touches the south pole (v ≈ π)
            let touches_south_pole = v_max > std::f64::consts::PI - pole_tol;

            if touches_north_pole || touches_south_pole {
                // Sample the UV polygon edges more densely near poles
                let pole_point = if touches_north_pole {
                    s.center + s.axis * s.radius // North pole
                } else {
                    s.center - s.axis * s.radius // South pole
                };

                // Sample UV edges
                let n = uv_poly.len();
                for i in 0..n {
                    let j = (i + 1) % n;
                    let a = uv_poly[i];
                    let b = uv_poly[j];

                    // More samples if edge is near pole
                    let near_pole = (a.y < pole_tol || a.y > std::f64::consts::PI - pole_tol)
                        || (b.y < pole_tol || b.y > std::f64::consts::PI - pole_tol);
                    let n_samples = if near_pole { 16 } else { 4 };

                    for k in 0..n_samples {
                        let t = k as f64 / n_samples as f64;
                        let uv = DVec2::new(
                            a.x + t * (b.x - a.x),
                            a.y + t * (b.y - a.y),
                        );

                        // Clamp v to avoid pole singularity
                        let v_clamped = uv.y.clamp(0.001, std::f64::consts::PI - 0.001);
                        let pt = s.point_at(uv.x, v_clamped);

                        // Skip points very close to pole (will add pole point separately)
                        if (pt - pole_point).length() > s.radius * 0.1 {
                            boundary_3d.push(pt);
                        }
                    }
                }

                // Add pole point if polygon contains it
                if touches_north_pole || touches_south_pole {
                    // Check if pole is inside the UV polygon
                    let pole_uv = if touches_north_pole {
                        DVec2::new(0.0, 0.0)
                    } else {
                        DVec2::new(0.0, std::f64::consts::PI)
                    };

                    // Add pole point at appropriate location
                    boundary_3d.push(pole_point);
                }
            } else {
                // No pole involvement - standard sampling
                for &uv in uv_poly {
                    boundary_3d.push(surface.point_at(uv.x, uv.y));
                }
            }

            // Deduplicate
            dedup_3d_points(&boundary_3d)
        }
        Surface3::Cone(c) => {
            // Cone has apex at v=0
            let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);

            if v_min < 0.01 {
                // Near apex - need special handling
                let apex = c.apex_point();
                let mut boundary_3d = Vec::new();

                let n = uv_poly.len();
                for i in 0..n {
                    let j = (i + 1) % n;
                    let a = uv_poly[i];
                    let b = uv_poly[j];

                    // Check if edge crosses near apex
                    let near_apex = a.y < 0.1 || b.y < 0.1;
                    let n_samples = if near_apex { 16 } else { 4 };

                    for k in 0..n_samples {
                        let t = k as f64 / n_samples as f64;
                        let uv = DVec2::new(
                            a.x + t * (b.x - a.x),
                            a.y + t * (b.y - a.y),
                        );

                        // Clamp v to avoid apex singularity
                        let v_clamped = uv.y.max(0.001);
                        let pt = c.point_at(uv.x, v_clamped);

                        // Skip points very close to apex
                        if (pt - apex).length() > 0.01 {
                            boundary_3d.push(pt);
                        }
                    }
                }

                // Add apex if polygon contains it
                boundary_3d.push(apex);

                dedup_3d_points(&boundary_3d)
            } else {
                // Standard case
                uv_poly.iter().map(|uv| surface.point_at(uv.x, uv.y)).collect()
            }
        }
        _ => {
            // No degenerate points - standard mapping
            uv_poly.iter().map(|uv| surface.point_at(uv.x, uv.y)).collect()
        }
    }
}

/// Enhanced handling of degenerate UV polygons on surfaces with singularities.
///
/// This function handles UV polygons where vertices collapse at surface singularities:
/// - Sphere poles (v=0 or v=π)
/// - Cone apex (v=0)
///
/// The function:
/// 1. Detects pole/apex proximity
/// 2. Handles triangulation specially for degenerate triangles
/// 3. Ensures edge PCurve tolerance near poles/apex
///
/// Returns a 3D boundary that correctly handles surface singularities.
pub fn handle_degenerate_uv_polygon(uv_poly: &[DVec2], surface: &Surface3) -> Vec<DVec3> {
    match surface {
        Surface3::Sphere(s) => {
            handle_sphere_degenerate_uv(uv_poly, s)
        }
        Surface3::Cone(c) => {
            handle_cone_degenerate_uv(uv_poly, c)
        }
        _ => {
            // No degenerate points - standard mapping
            uv_poly.iter().map(|uv| surface.point_at(uv.x, uv.y)).collect()
        }
    }
}

/// Handle degenerate UV polygons on sphere surfaces.
fn handle_sphere_degenerate_uv(uv_poly: &[DVec2], sphere: &SphericalSurface) -> Vec<DVec3> {
    let pole_tol = 0.01; // Tolerance for detecting near-pole

    // Find min/max v values to detect pole proximity
    let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let v_max = uv_poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);

    // Check if polygon touches either pole
    let touches_north_pole = v_min < pole_tol;
    let touches_south_pole = v_max > std::f64::consts::PI - pole_tol;

    if !touches_north_pole && !touches_south_pole {
        // No pole involvement - standard mapping
        return uv_poly.iter().map(|uv| sphere.point_at(uv.x, uv.y)).collect();
    }

    let mut boundary_3d = Vec::new();

    // Determine which pole(s) are involved
    let north_pole = sphere.center + sphere.axis * sphere.radius;
    let south_pole = sphere.center - sphere.axis * sphere.radius;

    // Sample UV polygon edges more densely near poles
    let n = uv_poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = uv_poly[i];
        let b = uv_poly[j];

        // More samples if edge is near pole
        let near_pole = (a.y < pole_tol || a.y > std::f64::consts::PI - pole_tol)
            || (b.y < pole_tol || b.y > std::f64::consts::PI - pole_tol);
        let n_samples = if near_pole { 16 } else { 4 };

        for k in 0..n_samples {
            let t = k as f64 / n_samples as f64;
            let uv = DVec2::new(
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
            );

            // Clamp v to avoid pole singularity
            let v_clamped = uv.y.clamp(0.001, std::f64::consts::PI - 0.001);
            let pt = sphere.point_at(uv.x, v_clamped);

            // Skip points very close to pole (will add pole point separately)
            let near_north = (pt - north_pole).length() < sphere.radius * 0.1;
            let near_south = (pt - south_pole).length() < sphere.radius * 0.1;
            if !near_north && !near_south {
                boundary_3d.push(pt);
            }
        }
    }

    // Add pole point(s) if polygon contains them
    if touches_north_pole {
        boundary_3d.push(north_pole);
    }
    if touches_south_pole {
        boundary_3d.push(south_pole);
    }

    dedup_3d_points(&boundary_3d)
}

/// Handle degenerate UV polygons on cone surfaces.
fn handle_cone_degenerate_uv(uv_poly: &[DVec2], cone: &ConicalSurface) -> Vec<DVec3> {
    let apex_tol = 0.01; // Tolerance for detecting near-apex

    // Find min v value to detect apex proximity
    let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);

    if v_min >= apex_tol {
        // No apex involvement - standard mapping
        return uv_poly.iter().map(|uv| cone.point_at(uv.x, uv.y)).collect();
    }

    let mut boundary_3d = Vec::new();
    let apex = cone.apex_point();

    // Sample UV polygon edges more densely near apex
    let n = uv_poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = uv_poly[i];
        let b = uv_poly[j];

        // More samples if edge is near apex
        let near_apex = a.y < apex_tol * 10.0 || b.y < apex_tol * 10.0;
        let n_samples = if near_apex { 16 } else { 4 };

        for k in 0..n_samples {
            let t = k as f64 / n_samples as f64;
            let uv = DVec2::new(
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
            );

            // Clamp v to avoid apex singularity
            let v_clamped = uv.y.max(0.001);
            let pt = cone.point_at(uv.x, v_clamped);

            // Skip points very close to apex
            if (pt - apex).length() > 0.01 {
                boundary_3d.push(pt);
            }
        }
    }

    // Add apex point
    boundary_3d.push(apex);

    dedup_3d_points(&boundary_3d)
}

/// Split an edge at a periodic seam if it crosses the U=0/2π boundary.
///
/// This function detects if an edge on a periodic surface (cylinder, sphere, torus)
/// crosses the seam and splits it at the crossing point.
///
/// Returns:
/// - `None` if the edge doesn't cross the seam
/// - `Some(vec![seg1, seg2])` where each segment is `[start_uv, end_uv]`
pub fn split_edge_at_periodic_seam(
    start_uv: DVec2,
    end_uv: DVec2,
    surface: &Surface3,
) -> Option<Vec<Vec<DVec2>>> {
    // Get the U period for this surface type
    let u_period = match surface {
        Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
            std::f64::consts::TAU
        }
        Surface3::Cone(_) => {
            // Cone is also periodic in U
            std::f64::consts::TAU
        }
        _ => {
            // Non-periodic surface
            return None;
        }
    };

    let u1 = start_uv.x;
    let u2 = end_uv.x;
    let v1 = start_uv.y;
    let v2 = end_uv.y;
    let du = u2 - u1;

    // Check for seam crossing (jump > period/2)
    if du.abs() <= u_period * 0.5 {
        return None;
    }

    // Determine which way we're crossing
    let is_low_to_high = du < 0.0; // u1 is high, u2 is low

    // Calculate intersection point at seam
    let (t, seam_u) = if is_low_to_high {
        // u1 is near period, u2 is near 0
        // Find t where u = period
        let t = (u_period - u1) / ((u2 + u_period) - u1);
        (t, u_period)
    } else {
        // u1 is near 0, u2 is near period
        // Find t where u = 0
        let t = -u1 / ((u2 - u_period) - u1);
        (t, 0.0)
    };

    // Clamp t to [0, 1] for numerical stability
    let t = t.clamp(0.0, 1.0);
    let seam_v = v1 + t * (v2 - v1);

    // Build two segments
    let seam_point = DVec2::new(seam_u, seam_v);
    let opposite_seam_point = if seam_u < u_period * 0.5 {
        DVec2::new(u_period, seam_v)
    } else {
        DVec2::new(0.0, seam_v)
    };

    // First segment: from start to seam
    let seg1 = vec![start_uv, seam_point];
    // Second segment: from opposite seam to end
    let seg2 = vec![opposite_seam_point, end_uv];

    Some(vec![seg1, seg2])
}

/// Split a UV polygon at both U and V seams for torus double periodicity.
///
/// The torus has two periodic parameters:
/// - U period: 2π (around major circle)
/// - V period: 2π (around tube circle)
///
/// This function handles UV polygon splitting in both directions.
pub fn split_uv_polygon_torus_double(uv_polygon: &[DVec2], period: f64) -> Vec<Vec<DVec2>> {
    if uv_polygon.len() < 3 {
        return vec![uv_polygon.to_vec()];
    }

    // First, split at U seam
    let u_split = split_uv_polygon_at_seam(uv_polygon, period);

    // Then, split each result at V seam
    let mut result = Vec::new();
    for poly in u_split {
        let v_split = split_uv_polygon_at_v_seam(&poly, period);
        result.extend(v_split);
    }

    result
}

/// Split a UV polygon at the V periodic seam (V=0/period boundary).
///
/// This is similar to split_uv_polygon_at_seam but for the V parameter.
fn split_uv_polygon_at_v_seam(uv_polygon: &[DVec2], period: f64) -> Vec<Vec<DVec2>> {
    if uv_polygon.len() < 3 {
        return vec![uv_polygon.to_vec()];
    }

    // Find all edges crossing the V seam
    let mut crossings: Vec<(usize, f64, DVec2)> = Vec::new();

    for i in 0..uv_polygon.len() {
        let j = (i + 1) % uv_polygon.len();
        let v1 = uv_polygon[i].y;
        let v2 = uv_polygon[j].y;
        let dv = v2 - v1;

        // Check for seam crossing (jump > period/2)
        if dv.abs() > period * 0.5 {
            let u1 = uv_polygon[i].x;
            let u2 = uv_polygon[j].x;

            // Determine which way we're crossing
            let is_low_to_high = dv < 0.0; // v1 is high, v2 is low

            // Calculate intersection point
            let (t, seam_v) = if is_low_to_high {
                let t = (period - v1) / ((v2 + period) - v1);
                (t, period)
            } else {
                let t = -v1 / ((v2 - period) - v1);
                (t, 0.0)
            };

            let t = t.clamp(0.0, 1.0);
            let seam_u = u1 + t * (u2 - u1);

            crossings.push((i, t, DVec2::new(seam_u, seam_v)));
        }
    }

    if crossings.is_empty() {
        return vec![uv_polygon.to_vec()];
    }

    // For now, handle simple cases
    if crossings.len() != 2 {
        // Complex case - return original
        return vec![uv_polygon.to_vec()];
    }

    // Build two sub-polygons
    let (idx1, _, pt1) = crossings[0];
    let (idx2, _, pt2) = crossings[1];

    let mut low_polygon: Vec<DVec2> = Vec::new();
    let mut high_polygon: Vec<DVec2> = Vec::new();

    let is_low = |v: f64| v < period * 0.5;

    let n = uv_polygon.len();

    // Traverse polygon and assign vertices
    for i in 0..n {
        let curr = uv_polygon[i];

        // Add current vertex to appropriate polygon
        if is_low(curr.y) {
            low_polygon.push(curr);
        } else {
            high_polygon.push(curr);
        }

        // Check for crossing between i and i+1
        for (cross_idx, _, cross_pt) in &crossings {
            if *cross_idx == i {
                // Add seam points to both polygons
                let low_pt = DVec2::new(cross_pt.x, 0.0);
                let high_pt = DVec2::new(cross_pt.x, period);

                if is_low(curr.y) {
                    low_polygon.push(low_pt);
                    high_polygon.push(high_pt);
                } else {
                    high_polygon.push(high_pt);
                    low_polygon.push(low_pt);
                }
            }
        }
    }

    let mut result = Vec::new();
    if low_polygon.len() >= 3 {
        result.push(low_polygon);
    }
    if high_polygon.len() >= 3 {
        result.push(high_polygon);
    }

    if result.is_empty() {
        vec![uv_polygon.to_vec()]
    } else {
        result
    }
}

/// Deduplicate 3D points within tolerance.
fn dedup_3d_points(points: &[DVec3]) -> Vec<DVec3> {
    let mut result: Vec<DVec3> = Vec::new();
    let tol_sq = TOLERANCE_ABS * TOLERANCE_ABS;

    for &p in points {
        if result.iter().all(|q: &DVec3| (p - *q).length_squared() > tol_sq) {
            result.push(p);
        }
    }

    result
}

/// Check if a UV polygon is valid (has sufficient area and no degenerate edges).
fn is_valid_uv_polygon(poly: &[DVec2]) -> bool {
    if poly.len() < 3 {
        return false;
    }

    // Check for sufficient area (shoelace formula)
    let mut area = 0.0;
    let n = poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        area += poly[i].x * poly[j].y;
        area -= poly[j].x * poly[i].y;
    }
    area = area.abs() * 0.5;

    // Area should be significant
    area > 1e-10
}

fn periodic_trim_to_open_isoline(poly: &[DVec2], trim: &[DVec2], u_period: f64) -> Option<Vec<DVec2>> {
    if poly.len() < 3 || trim.len() < 3 || u_period <= 0.0 {
        return None;
    }

    let trim_start = trim[0];
    let trim_end = trim[trim.len() - 1];
    let is_closed = (trim_start - trim_end).length_squared() < 1e-6;
    if !is_closed {
        return None;
    }

    let u_min_trim = trim.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let u_max_trim = trim.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let v_min_trim = trim.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let v_max_trim = trim.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    let u_span = u_max_trim - u_min_trim;
    let v_span = v_max_trim - v_min_trim;

    let poly_u_min = poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let poly_u_max = poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let poly_v_min = poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let poly_v_max = poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    let poly_v_span = poly_v_max - poly_v_min;

    if u_span < 0.9 * u_period {
        return None;
    }
    if poly_v_span <= 1e-12 || v_span > 0.1 * poly_v_span {
        return None;
    }

    let v_level = trim.iter().map(|p| p.y).sum::<f64>() / trim.len() as f64;
    if v_level <= poly_v_min + 1e-9 || v_level >= poly_v_max - 1e-9 {
        return None;
    }

    Some(vec![
        DVec2::new(poly_u_min, v_level),
        DVec2::new(poly_u_max, v_level),
    ])
}

/// Split a UV polygon at periodic seams (U=0/period boundary).
///
/// For periodic surfaces like cylinders, the U parameter wraps around.
/// When a polygon crosses the seam (U=0 or U=period), we need to split it
/// into separate polygons, each with consistent U coordinates.
///
/// Algorithm:
/// 1. Find edges that cross the seam (|du| > period * 0.5)
/// 2. For each crossing edge, compute the exact intersection point at U=0 or U=period
/// 3. Build output polygons by inserting intersection points
///
/// Returns one or more polygons that don't cross the seam.
pub fn split_uv_polygon_at_seam(uv_polygon: &[DVec2], period: f64) -> Vec<Vec<DVec2>> {
    if uv_polygon.len() < 3 {
        return vec![uv_polygon.to_vec()];
    }

    // Structure to hold information about seam crossings
    struct SeamCrossing {
        edge_idx: usize,
        intersection: DVec2,
        is_low_to_high: bool, // true: crossing from low u (near 0) to high u (near period)
    }

    // Find all edges crossing the seam and compute intersection points
    let mut crossings: Vec<SeamCrossing> = Vec::new();
    for i in 0..uv_polygon.len() {
        let j = (i + 1) % uv_polygon.len();
        let u1 = uv_polygon[i].x;
        let u2 = uv_polygon[j].x;
        let v1 = uv_polygon[i].y;
        let v2 = uv_polygon[j].y;
        let du = u2 - u1;

        // Large jump indicates seam crossing
        if du.abs() > period * 0.5 {
            // Determine which way we're crossing
            // du > 0: wrapping from low u to high u (crossing U=0 going backwards in unwrapped space)
            // du < 0: wrapping from high u to low u (crossing U=period going backwards in unwrapped space)
            let is_low_to_high = du < 0.0; // u1 is high, u2 is low

            // Calculate intersection point using linear interpolation
            // We need to find the V coordinate where the edge crosses the seam
            //
            // For an edge from (u1, v1) to (u2, v2) crossing the seam:
            // If u1 is near period and u2 is near 0: unwrap u2 to u2 + period, find where U = period
            // If u1 is near 0 and u2 is near period: unwrap u2 to u2 - period, find where U = 0
            let (t, seam_u) = if is_low_to_high {
                // u1 is near period, u2 is near 0
                // Unwrap u2: consider edge from (u1, v1) to (u2 + period, v2)
                // Find t where u = period
                let t = (period - u1) / ((u2 + period) - u1);
                (t, period)
            } else {
                // u1 is near 0, u2 is near period
                // Unwrap u2: consider edge from (u1, v1) to (u2 - period, v2)
                // Find t where u = 0 (which equals period in the unwrapped space)
                // Or equivalently: the edge goes from u1 to u2-period (negative)
                // We want u = 0, so t = (0 - u1) / ((u2 - period) - u1) = -u1 / (u2 - period - u1)
                let t = -u1 / ((u2 - period) - u1);
                (t, 0.0)
            };

            // Clamp t to [0, 1] to handle numerical edge cases
            let t = t.clamp(0.0, 1.0);
            let intersection_v = v1 + t * (v2 - v1);

            crossings.push(SeamCrossing {
                edge_idx: i,
                intersection: DVec2::new(seam_u, intersection_v),
                is_low_to_high,
            });
        }
    }

    if crossings.is_empty() {
        return vec![uv_polygon.to_vec()];
    }

    // Build output polygons
    // We need to partition the vertices and insert intersection points
    // Each output polygon will have consistent U values (all low or all high)

    // Collect all vertices and their positions relative to the seam
    // "low" means u < period * 0.5, "high" means u >= period * 0.5
    let is_low = |u: f64| u < period * 0.5;

    // Build two polygons: one for low-u region, one for high-u region
    let mut low_polygon: Vec<DVec2> = Vec::new();
    let mut high_polygon: Vec<DVec2> = Vec::new();

    // Sort crossings by edge index for efficient lookup
    let crossing_map: std::collections::HashMap<usize, &SeamCrossing> = crossings
        .iter()
        .map(|c| (c.edge_idx, c))
        .collect();

    // Traverse the polygon and assign vertices to appropriate output polygons
    for i in 0..uv_polygon.len() {
        let curr = uv_polygon[i];
        let next_idx = (i + 1) % uv_polygon.len();
        let next = uv_polygon[next_idx];

        // Add current vertex to appropriate polygon
        if is_low(curr.x) {
            low_polygon.push(curr);
        } else {
            high_polygon.push(curr);
        }

        // Check if edge (i, i+1) crosses the seam
        if let Some(crossing) = crossing_map.get(&i) {
            // Add intersection point to both polygons
            // The intersection point is at the seam (u = 0 or u = period)
            // For the low polygon, we want u = 0
            // For the high polygon, we want u = period
            let low_intersection = DVec2::new(0.0, crossing.intersection.y);
            let high_intersection = DVec2::new(period, crossing.intersection.y);

            if crossing.is_low_to_high {
                // Going from high u to low u
                // Add period-point to high polygon first, then 0-point to low polygon
                high_polygon.push(high_intersection);
                low_polygon.push(low_intersection);
            } else {
                // Going from low u to high u
                // Add 0-point to low polygon first, then period-point to high polygon
                low_polygon.push(low_intersection);
                high_polygon.push(high_intersection);
            }
        }
    }

    // Build result - only include valid polygons (at least 3 vertices)
    let mut result = Vec::new();

    if low_polygon.len() >= 3 {
        result.push(low_polygon);
    }
    if high_polygon.len() >= 3 {
        result.push(high_polygon);
    }

    // If we didn't get valid polygons, return the original
    if result.is_empty() {
        return vec![uv_polygon.to_vec()];
    }

    result
}

/// Split a 2D UV polygon by a 2D trim polyline.
///
/// Algorithm:
/// 1. Find trim start/end's closest edge on the polygon boundary.
/// 2. Project trim endpoints onto boundary edges to find exact split points.
/// 3. Split polygon into two halves at those points, inserting the trim polyline
///    between them.
///
/// For closed trim polylines (start ≈ end), uses a closed-curve splitting
/// algorithm: the trim forms an interior polygon that divides the outer polygon
/// into "inside trim" and "outside trim" regions.
///
/// Returns 1 polygon if splitting is degenerate, or 2 sub-polygons otherwise.
fn split_uv_polygon_by_trim(poly: &[DVec2], trim: &[DVec2]) -> Vec<Vec<DVec2>> {
    let n = poly.len();
    if n < 3 || trim.len() < 2 {
        return vec![poly.to_vec()];
    }

    let trim_start = trim[0];
    let trim_end = trim[trim.len() - 1];

    // Detect truly-closed trim: start ≈ end in UV space (e.g. a small loop entirely
    // inside the face).  Wrapped-closed trims (start and end differ by ~2π in u,
    // representing a full-circle cut around a cylinder or sphere) are intentionally
    // NOT treated as closed loops here — they are open trims whose endpoints lie on
    // opposite sides of the UV boundary seam and should split the face into two bands.
    let is_closed_trim = (trim_start - trim_end).length_squared() < 1e-6;
    if is_closed_trim {
        // The trim is a truly closed loop entirely inside the polygon.
        // Use the trim as an interior boundary and return [trim_polygon, outer_polygon].
        let trim_centroid = trim.iter().copied().sum::<DVec2>() / trim.len() as f64;
        let is_inside = point_in_polygon_2d(poly, trim_centroid);
        if is_inside {
            let mut trim_dedup: Vec<DVec2> = trim.to_vec();
            if trim_dedup.len() > 1
                && (trim_dedup[0] - trim_dedup[trim_dedup.len() - 1]).length_squared() < 1e-12
            {
                trim_dedup.pop();
            }
            if trim_dedup.len() >= 3 {
                return vec![trim_dedup, poly.to_vec()];
            }
        }
        return vec![poly.to_vec()];
    }

    // Find closest point on each polygon edge for a query point.
    // Returns (edge_index, t_param, projected_point).
    let closest_on_boundary = |q: DVec2| -> (usize, f64, DVec2) {
        let mut best_edge = 0usize;
        let mut best_t = 0.0f64;
        let mut best_pt = poly[0];
        let mut best_dist = f64::INFINITY;
        for i in 0..n {
            let j = (i + 1) % n;
            let a = poly[i];
            let b = poly[j];
            let ab = b - a;
            let len_sq = ab.dot(ab);
            let t = if len_sq < 1e-14 {
                0.0
            } else {
                ((q - a).dot(ab) / len_sq).clamp(0.0, 1.0)
            };
            let proj = a + t * ab;
            let dist = (q - proj).length_squared();
            if dist < best_dist {
                best_dist = dist;
                best_edge = i;
                best_t = t;
                best_pt = proj;
            }
        }
        (best_edge, best_t, best_pt)
    };

    // Cast a 2D ray from `origin` along `dir` and return the first boundary edge
    // intersection with t > -eps (including slightly behind for on-boundary starts).
    // Returns None if no intersection is found within a reasonable range.
    let ray_to_boundary = |origin: DVec2, dir: DVec2| -> Option<(usize, DVec2)> {
        let dir_len = dir.length();
        if dir_len < 1e-12 {
            return None;
        }
        let dir = dir / dir_len;
        let mut best_t = f64::INFINITY;
        let mut best_edge = 0usize;
        let mut best_pt = poly[0];
        for i in 0..n {
            let j = (i + 1) % n;
            let a = poly[i];
            let b = poly[j];
            let ab = b - a;
            // Solve: origin + t*dir = a + s*ab
            // => t*(dir×ab) = (a-origin)×ab  (2D cross: x.x*y.y - x.y*y.x)
            let denom = dir.x * ab.y - dir.y * ab.x;
            if denom.abs() < 1e-14 {
                continue; // parallel
            }
            let oa = a - origin;
            let t_ray = (oa.x * ab.y - oa.y * ab.x) / denom;
            let s_seg = (oa.x * dir.y - oa.y * dir.x) / denom;
            if t_ray > -1e-9 && s_seg >= -1e-9 && s_seg <= 1.0 + 1e-9 && t_ray < best_t {
                best_t = t_ray;
                best_edge = i;
                best_pt = a + s_seg.clamp(0.0, 1.0) * ab;
            }
        }
        if best_t.is_finite() {
            Some((best_edge, best_pt))
        } else {
            None
        }
    };

    // Compute UV polygon bounding box to compute a "near-boundary" threshold
    let u_span = poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max)
        - poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let v_span = poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max)
        - poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let boundary_snap_tol = (u_span + v_span) * 0.05;

    // For each trim endpoint: if it lies close to the boundary already, use closest_on_boundary.
    // Otherwise, extrapolate along the trim tangent to find the proper boundary edge.
    let locate_endpoint =
        |endpoint: DVec2, tangent_from: DVec2| -> (usize, DVec2) {
            let (_, _, proj) = closest_on_boundary(endpoint);
            let dist_to_bnd = (endpoint - proj).length();
            if dist_to_bnd <= boundary_snap_tol {
                // Already on/near boundary
                let (edge, _, pt) = closest_on_boundary(endpoint);
                (edge, pt)
            } else {
                // Interior endpoint — cast ray along trim tangent toward boundary
                let tang = (endpoint - tangent_from).normalize_or_zero();
                if let Some((edge, pt)) = ray_to_boundary(endpoint, tang) {
                    (edge, pt)
                } else {
                    // Fallback to closest projection
                    let (edge, _, pt) = closest_on_boundary(endpoint);
                    (edge, pt)
                }
            }
        };

    let interior_from_start = if trim.len() >= 2 { trim[1] } else { trim_end };
    let interior_from_end = if trim.len() >= 2 { trim[trim.len() - 2] } else { trim_start };

    let (edge_s, pt_s) = locate_endpoint(trim_start, interior_from_start);
    let (edge_e, pt_e) = locate_endpoint(trim_end, interior_from_end);

    // Ensure ia <= ib for consistent polygon walking
    let (ia, ib, p_a, p_b, trim_forward) = if edge_s <= edge_e {
        (edge_s, edge_e, pt_s, pt_e, true)
    } else {
        (edge_e, edge_s, pt_e, pt_s, false)
    };

    if ia == ib {
        // Both endpoints project to the same edge — degenerate, can't split
        return vec![poly.to_vec()];
    }

    // Build the trim points in the correct order for each half
    let trim_pts: Vec<DVec2> = if trim_forward {
        trim.to_vec()
    } else {
        trim.iter().copied().rev().collect()
    };

    // Sub-polygon A: poly[0..=ia] + p_a + trim_pts (interior only) + p_b + poly[ib+1..]
    let mut sub_a: Vec<DVec2> = poly[..=ia].to_vec();
    sub_a.push(p_a);
    // Interior trim points (skip first and last which are endpoints)
    for &p in trim_pts.iter().skip(1).rev().skip(1) {
        sub_a.push(p);
    }
    sub_a.push(p_b);
    sub_a.extend_from_slice(&poly[ib + 1..]);

    // Sub-polygon B: p_a + poly[ia+1..=ib] + p_b + trim_pts reversed (interior only)
    let mut sub_b: Vec<DVec2> = vec![p_a];
    sub_b.extend_from_slice(&poly[ia + 1..=ib]);
    sub_b.push(p_b);
    for &p in trim_pts.iter().skip(1).rev().skip(1).rev() {
        sub_b.push(p);
    }

    // Deduplicate consecutive near-equal points
    let dedup_2d = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut result: Vec<DVec2> = Vec::new();
        for p in v {
            if result.is_empty() || (p - result[result.len() - 1]).length_squared() > 1e-18 {
                result.push(p);
            }
        }
        if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < 1e-18 {
            result.pop();
        }
        result
    };

    let sub_a = dedup_2d(sub_a);
    let sub_b = dedup_2d(sub_b);

    let mut out = Vec::new();
    if sub_a.len() >= 3 {
        out.push(sub_a);
    }
    if sub_b.len() >= 3 {
        out.push(sub_b);
    }

    if out.is_empty() {
        vec![poly.to_vec()]
    } else {
        out
    }
}
/// Returns true when the circle is fully embedded inside the polygon.
///
/// This corresponds to the case where all polygon vertices are outside the
/// circle while the circle center lies inside the polygon.
fn circle_fully_inside_polygon(poly: &[DVec2], center: DVec2, radius: f64) -> bool {
    if poly.len() < 3 {
        return false;
    }
    let tol = TOLERANCE_ABS;
    let all_vertices_outside = poly.iter().all(|&p| (p - center).length() >= radius - tol);
    all_vertices_outside && point_in_polygon_2d(poly, center)
}

/// Split a 2D polygon by a circle boundary.
///
/// Vertices inside the circle (distance < radius) are on the "inside" group,
/// vertices outside (distance > radius) are on the "outside" group.
/// Returns up to 2 sub-polygons: the part inside and the part outside.
///
/// When the circle is fully inside the polygon (all vertices outside),
/// samples the circle at N_CIRCLE_SAMPLES points and returns both
/// the approximate circular cap and the annular region.
fn split_polygon_by_circle_2d(poly: &[DVec2], center: DVec2, radius: f64) -> Vec<Vec<DVec2>> {
    const N_CIRCLE_SAMPLES: usize = 24;
    let n = poly.len();
    if n < 3 {
        return vec![poly.to_vec()];
    }

    let tol = TOLERANCE_ABS;

    // Signed distance: negative = inside circle, positive = outside
    let signed_dist = |p: DVec2| -> f64 { (p - center).length() - radius };

    let dists: Vec<f64> = poly.iter().map(|&p| signed_dist(p)).collect();

    // Check if all vertices are on the same side
    let all_inside = dists.iter().all(|&d| d <= tol);
    let all_outside = dists.iter().all(|&d| d >= -tol);

    if all_inside {
        // All polygon vertices inside circle — keep whole polygon
        return vec![poly.to_vec()];
    }

    if all_outside {
        // Circle is fully inside the polygon OR polygon is fully outside circle.
        // Check if circle center is inside the polygon:
        let center_inside = point_in_polygon_2d(poly, center);
        if !center_inside {
            // Circle doesn't overlap with this polygon — keep as-is
            return vec![poly.to_vec()];
        }
        // Circle is fully inside the polygon — produce circular cap + annular region
        // Sample the circle at N points
        let circle_poly: Vec<DVec2> = (0..N_CIRCLE_SAMPLES)
            .map(|i| {
                let theta = std::f64::consts::TAU * i as f64 / N_CIRCLE_SAMPLES as f64;
                center + DVec2::new(theta.cos(), theta.sin()) * radius
            })
            .collect();
        // Return: inside = circle polygon, outside = original polygon (with circle as hole)
        // For simplicity, return just the circle as the "inside" part
        // and the original polygon as the "outside" part (approximate)
        return vec![circle_poly, poly.to_vec()];
    }

    // Find crossings: edges where signed distance changes sign
    let mut crossings: Vec<(usize, DVec2)> = Vec::new();

    for i in 0..n {
        let j = (i + 1) % n;
        let di = dists[i];
        let dj = dists[j];

        if di.abs() < tol {
            continue; // vertex i is on the circle
        }
        if dj.abs() < tol {
            continue; // vertex j is on the circle (handled when edge starting at j is processed)
        }

        if di * dj < 0.0 {
            // Edge crosses the circle boundary
            // Find exact crossing: solve |a + t*(b-a) - center|² = r²
            let a = poly[i];
            let b = poly[j];
            let ab = b - a;
            let ac = a - center;
            let qa = ab.dot(ab);
            let qb = 2.0 * ab.dot(ac);
            let qc = ac.dot(ac) - radius * radius;
            let disc = qb * qb - 4.0 * qa * qc;
            if disc < 0.0 {
                continue;
            }
            let sq = disc.sqrt();
            for &sign in &[-1.0_f64, 1.0_f64] {
                let t = (-qb + sign * sq) / (2.0 * qa);
                if t > -tol && t < 1.0 + tol {
                    let t = t.clamp(0.0, 1.0);
                    let pt = a + t * ab;
                    crossings.push((i, pt));
                    break; // take the first valid crossing on this edge
                }
            }
        }
    }

    if crossings.len() < 2 {
        // Can't split — keep as-is
        return vec![poly.to_vec()];
    }

    // Sort crossings by edge index
    crossings.sort_by_key(|(idx, _)| *idx);

    // Take the first two crossings
    let (idx1, pt1) = crossings[0];
    let (idx2, pt2) = crossings[1];

    if idx1 == idx2 {
        return vec![poly.to_vec()];
    }

    // Sample the arc between pt1 and pt2 (going through the inside of the polygon)
    // Determine which arc (minor or major) connects pt1 to pt2 and stays inside the polygon
    let theta1 = (pt1 - center).to_angle();
    let theta2 = (pt2 - center).to_angle();

    // For the "inside" sub-polygon, we need the arc that passes through the inside of the polygon.
    // Try both arcs and pick the one whose midpoint is inside the polygon.
    let mid_theta_cw = (theta1 + theta2) * 0.5;
    let mid_theta_ccw = mid_theta_cw + std::f64::consts::PI;
    let mid_cw = center + DVec2::new(mid_theta_cw.cos(), mid_theta_cw.sin()) * radius;
    let _mid_ccw = center + DVec2::new(mid_theta_ccw.cos(), mid_theta_ccw.sin()) * radius;

    // The arc midpoint that is inside the polygon corresponds to the "inside" portion
    let arc_goes_cw_inside = point_in_polygon_2d(poly, mid_cw);
    let inner_mid_theta = if arc_goes_cw_inside {
        mid_theta_cw
    } else {
        mid_theta_ccw
    };

    // Determine angular span and direction for the inner arc
    let arc_n = ((N_CIRCLE_SAMPLES as f64 * (theta2 - theta1).abs() / std::f64::consts::TAU)
        as usize)
        .max(3);

    // Build arc points from pt1 to pt2 going through inner_mid_theta
    let inner_arc: Vec<DVec2> = {
        // Compute proper arc from theta1 through inner_mid_theta to theta2
        let delta = {
            let mut d = theta2 - theta1;
            // Adjust delta to go through inner_mid_theta
            let going_ccw = inner_mid_theta > theta1 || inner_mid_theta < theta2;
            if going_ccw {
                while d < 0.0 {
                    d += std::f64::consts::TAU;
                }
                if d > std::f64::consts::TAU {
                    d -= std::f64::consts::TAU;
                }
            } else {
                while d > 0.0 {
                    d -= std::f64::consts::TAU;
                }
                if d < -std::f64::consts::TAU {
                    d += std::f64::consts::TAU;
                }
            }
            d
        };
        (0..=arc_n)
            .map(|i| {
                let t = i as f64 / arc_n as f64;
                let theta = theta1 + delta * t;
                center + DVec2::new(theta.cos(), theta.sin()) * radius
            })
            .collect()
    };

    // Sub-polygon "inside" (circle side): pt1 → arc → pt2 + polygon walk from idx2 to idx1
    // Actually: vertices of polygon that are INSIDE the circle + arc from pt1 to pt2
    let poly_inside_verts: Vec<DVec2> = poly[idx1 + 1..=idx2].to_vec();

    let mut sub_inside: Vec<DVec2> = vec![pt1];
    sub_inside.extend_from_slice(&poly_inside_verts);
    sub_inside.push(pt2);
    // Add arc back (reversed, so the boundary goes: inside polygon verts, then arc back to pt1)
    for &p in inner_arc.iter().skip(1).rev().skip(1) {
        sub_inside.push(p);
    }

    // Sub-polygon "outside" (non-circle side): pt2 → arc → pt1 + polygon walk
    let poly_outside_verts_a: Vec<DVec2> = poly[..=idx1].to_vec();
    let poly_outside_verts_b: Vec<DVec2> = poly[idx2 + 1..].to_vec();

    let mut sub_outside: Vec<DVec2> = poly_outside_verts_a;
    sub_outside.push(pt1);
    // Add inner arc (forward) as the "hole" boundary
    for &p in inner_arc.iter().skip(1).rev().skip(1) {
        sub_outside.push(p);
    }
    sub_outside.push(pt2);
    sub_outside.extend(poly_outside_verts_b);

    let dedup = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut result: Vec<DVec2> = Vec::new();
        for p in v {
            if result.is_empty() || (p - result[result.len() - 1]).length_squared() > 1e-18 {
                result.push(p);
            }
        }
        if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < 1e-18 {
            result.pop();
        }
        result
    };

    let sub_inside = dedup(sub_inside);
    let sub_outside = dedup(sub_outside);

    let mut out = Vec::new();
    if sub_inside.len() >= 3 {
        out.push(sub_inside);
    }
    if sub_outside.len() >= 3 {
        out.push(sub_outside);
    }

    if out.is_empty() {
        vec![poly.to_vec()]
    } else {
        out
    }
}

/// Check if a 2D point is inside a 2D polygon using ray casting.
fn point_in_polygon_2d(poly: &[DVec2], pt: DVec2) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let vi = poly[i];
        let vj = poly[j];
        if ((vi.y > pt.y) != (vj.y > pt.y))
            && (pt.x < (vj.x - vi.x) * (pt.y - vi.y) / (vj.y - vi.y) + vi.x)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Split a 2D polygon by an infinite line through `point` with direction `dir`.
///
/// Vertices on the positive side (cross product > 0) form one group, negative side the other.
fn split_polygon_2d_by_line(poly: &[DVec2], point: DVec2, dir: DVec2) -> Vec<Vec<DVec2>> {
    let n = poly.len();
    if n < 3 {
        return vec![poly.to_vec()];
    }
    let tol = TOLERANCE_ABS;

    // Signed distance from line
    let signed_dist = |p: DVec2| -> f64 {
        let d = p - point;
        dir.x * d.y - dir.y * d.x // perpendicular component
    };

    let dists: Vec<f64> = poly.iter().map(|&p| signed_dist(p)).collect();
    let all_pos = dists.iter().all(|&d| d >= -tol);
    let all_neg = dists.iter().all(|&d| d <= tol);

    if all_pos || all_neg {
        return vec![poly.to_vec()];
    }

    let mut crossings: Vec<(usize, DVec2)> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let di = dists[i];
        let dj = dists[j];
        if di.abs() < tol || dj.abs() < tol {
            continue;
        }
        if di * dj < 0.0 {
            let t = di / (di - dj);
            let pt = poly[i] + t * (poly[j] - poly[i]);
            crossings.push((i, pt));
        }
    }

    if crossings.len() < 2 {
        return vec![poly.to_vec()];
    }

    crossings.sort_by_key(|(idx, _)| *idx);

    let (idx1, pt1) = crossings[0];
    let (idx2, pt2) = crossings[1];
    if idx1 == idx2 {
        return vec![poly.to_vec()];
    }

    let mut sub_a: Vec<DVec2> = poly[..=idx1].to_vec();
    sub_a.push(pt1);
    sub_a.push(pt2);
    sub_a.extend_from_slice(&poly[idx2 + 1..]);

    let mut sub_b: Vec<DVec2> = vec![pt1];
    sub_b.extend_from_slice(&poly[idx1 + 1..=idx2]);
    sub_b.push(pt2);

    let dedup = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut result: Vec<DVec2> = Vec::new();
        for p in v {
            if result.is_empty() || (p - result[result.len() - 1]).length_squared() > 1e-18 {
                result.push(p);
            }
        }
        if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < 1e-18 {
            result.pop();
        }
        result
    };

    let sub_a = dedup(sub_a);
    let sub_b = dedup(sub_b);
    let mut out = Vec::new();
    if sub_a.len() >= 3 {
        out.push(sub_a);
    }
    if sub_b.len() >= 3 {
        out.push(sub_b);
    }
    if out.is_empty() {
        vec![poly.to_vec()]
    } else {
        out
    }
}

/// Split a 2D polygon by a segment from `seg_start` to `seg_end`.
fn split_polygon_2d_by_segment(
    poly: &[DVec2],
    seg_start: DVec2,
    seg_end: DVec2,
) -> Vec<Vec<DVec2>> {
    let dir = seg_end - seg_start;
    if dir.length_squared() < 1e-18 {
        return vec![poly.to_vec()];
    }
    split_polygon_2d_by_line(poly, seg_start, dir.normalize())
}

// ============================================================================
// Glue Path Enhancement Types and Functions
// ============================================================================

/// Configuration for glue-based boolean operations.
///
/// This struct controls the behavior of the shared-face fast path (glue option)
/// for boolean operations. When two shapes have coincident or near-coincident
/// faces at their interface, the glue path can skip expensive intersection
/// computations and directly merge the topology.
///
/// # Example
///
/// ```
/// use rcad_algorithms::builder::GlueConfig;
///
/// let config = GlueConfig {
///     face_tolerance: 1e-5,
///     edge_tolerance: 1e-5,
///     use_geometric_hash: true,
///     early_normal_filter: true,
/// };
/// ```
#[derive(Debug, Clone)]
pub struct GlueConfig {
    /// Tolerance for face matching (default: 1e-6).
    ///
    /// Two faces are considered coincident if their surface geometry
    /// matches within this tolerance.
    pub face_tolerance: f64,

    /// Tolerance for edge matching (default: 1e-6).
    ///
    /// Two edges are considered coincident if their curve geometry
    /// matches within this tolerance.
    pub edge_tolerance: f64,

    /// Enable geometric hashing for O(n) face pairing (default: true).
    ///
    /// When enabled, uses a spatial hash to quickly find candidate face
    /// pairs, reducing the complexity from O(n²) to O(n) for models
    /// with many faces.
    pub use_geometric_hash: bool,

    /// Skip non-parallel face pairs early (default: true).
    ///
    /// When enabled, quickly rejects face pairs whose normals are not
    /// approximately anti-parallel, avoiding more expensive geometric
    /// compatibility checks.
    pub early_normal_filter: bool,
}

impl Default for GlueConfig {
    fn default() -> Self {
        Self {
            face_tolerance: TOLERANCE_ABS,
            edge_tolerance: TOLERANCE_ABS,
            use_geometric_hash: true,
            early_normal_filter: true,
        }
    }
}

/// Result of glue face detection.
///
/// Represents a pair of faces from two different shapes that have been
/// identified as coincident or near-coincident, suitable for glue-based
/// boolean operations.
#[derive(Debug, Clone)]
pub struct GlueFacePair {
    /// Index of face in shape A.
    pub face_a: usize,

    /// Index of face in shape B.
    pub face_b: usize,

    /// Match quality (1.0 = perfect match).
    ///
    /// This value indicates how well the two faces match:
    /// - 1.0: Perfect geometric match
    /// - 0.9-1.0: Near-perfect match, within tolerance
    /// - 0.7-0.9: Partial match, some deviation
    /// - < 0.7: Poor match, may not be suitable for gluing
    pub match_quality: f64,

    /// Estimated area of shared region.
    ///
    /// For fully coincident faces, this is the face area.
    /// For partially overlapping faces, this is the overlap area.
    pub shared_area: f64,
}

/// Geometric hash cell for face center points.
///
/// Used for O(n) face pairing by hashing face center coordinates
/// into spatial cells.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct GeomHashCell {
    ix: i64,
    iy: i64,
    iz: i64,
}

impl GeomHashCell {
    fn from_point(p: DVec3, cell_size: f64) -> Self {
        let scale = 1.0 / cell_size;
        Self {
            ix: (p.x * scale).round() as i64,
            iy: (p.y * scale).round() as i64,
            iz: (p.z * scale).round() as i64,
        }
    }
}

/// Face-pairing cache for performance.
///
/// Caches the results of face compatibility checks to avoid
/// redundant computations during boolean operations.
#[derive(Debug, Clone, Default)]
pub struct GlueFaceCache {
    /// Cached face center points for each face.
    face_centers: Vec<DVec3>,

    /// Cached face normals for each face.
    face_normals: Vec<DVec3>,

    /// Cached face areas for each face.
    face_areas: Vec<f64>,

    /// Spatial hash mapping cells to face indices.
    spatial_hash: HashMap<GeomHashCell, Vec<usize>>,

    /// Cached surface compatibility results.
    /// Key: (face_a, face_b), Value: is_compatible
    compatibility_cache: HashMap<(usize, usize), bool>,
}

impl GlueFaceCache {
    /// Create a new empty cache.
    pub fn new() -> Self {
        Self::default()
    }

    /// Build the cache for a BRep by computing face centers, normals, and areas.
    pub fn build(&mut self, brep: &BRep, cell_size: f64) {
        self.face_centers.clear();
        self.face_normals.clear();
        self.face_areas.clear();
        self.spatial_hash.clear();
        self.compatibility_cache.clear();

        let mut face_idx = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    // Compute face center and area from boundary vertices
                    let mut center = DVec3::ZERO;
                    let mut area = 0.0;
                    let mut count = 0usize;

                    for we in &face.outer_wire.edges {
                        if we.idx < brep.edges.len() {
                            let edge = &brep.edges[we.idx];
                            if edge.start < brep.vertices.len() {
                                center += brep.vertices[edge.start].point;
                                count += 1;
                            }
                            if edge.end < brep.vertices.len() {
                                center += brep.vertices[edge.end].point;
                                count += 1;
                            }
                        }
                    }

                    if count > 0 {
                        center /= count as f64;
                    }

                    // Approximate area from bounding box
                    let mut min_pt = DVec3::splat(f64::INFINITY);
                    let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
                    for we in &face.outer_wire.edges {
                        if we.idx < brep.edges.len() {
                            let edge = &brep.edges[we.idx];
                            if edge.start < brep.vertices.len() {
                                let p = brep.vertices[edge.start].point;
                                min_pt = min_pt.min(p);
                                max_pt = max_pt.max(p);
                            }
                            if edge.end < brep.vertices.len() {
                                let p = brep.vertices[edge.end].point;
                                min_pt = min_pt.min(p);
                                max_pt = max_pt.max(p);
                            }
                        }
                    }
                    let diag = max_pt - min_pt;
                    area = diag.x * diag.y + diag.y * diag.z + diag.z * diag.x;

                    self.face_centers.push(center);
                    self.face_normals.push(face.normal);
                    self.face_areas.push(area);

                    // Add to spatial hash
                    let cell = GeomHashCell::from_point(center, cell_size);
                    self.spatial_hash.entry(cell).or_default().push(face_idx);

                    face_idx += 1;
                }
            }
        }
    }

    /// Get nearby faces using spatial hash.
    pub fn get_nearby_faces(&self, center: DVec3, cell_size: f64) -> Vec<usize> {
        let cell = GeomHashCell::from_point(center, cell_size);

        // Check the cell and its neighbors
        let mut result = Vec::new();
        for dx in -1i64..=1 {
            for dy in -1i64..=1 {
                for dz in -1i64..=1 {
                    let neighbor = GeomHashCell {
                        ix: cell.ix + dx,
                        iy: cell.iy + dy,
                        iz: cell.iz + dz,
                    };
                    if let Some(faces) = self.spatial_hash.get(&neighbor) {
                        result.extend(faces.iter().copied());
                    }
                }
            }
        }
        result
    }

    /// Check if surface compatibility is cached.
    pub fn get_compatibility(&self, face_a: usize, face_b: usize) -> Option<bool> {
        self.compatibility_cache.get(&(face_a, face_b)).copied()
    }

    /// Cache a surface compatibility result.
    pub fn set_compatibility(&mut self, face_a: usize, face_b: usize, compatible: bool) {
        self.compatibility_cache.insert((face_a, face_b), compatible);
        self.compatibility_cache.insert((face_b, face_a), compatible);
    }
}

/// Detect glue face pairs between two shapes.
///
/// This function analyzes two BReps and identifies pairs of faces that
/// are geometrically coincident or near-coincident, suitable for the
/// glue-based boolean fast path.
///
/// # Arguments
///
/// * `brep_a` - First BRep shape.
/// * `brep_b` - Second BRep shape.
/// * `config` - Configuration for glue detection.
///
/// # Returns
///
/// A vector of `GlueFacePair` representing detected coincident face pairs.
///
/// # Example
///
/// ```
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use rcad_algorithms::builder::{GlueConfig, detect_glue_faces};
/// use glam::DAffine3;
///
/// let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let mut box2 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// box2.apply_transform(DAffine3::from_translation(glam::DVec3::new(0.0, 1.0, 0.0)));
///
/// let config = GlueConfig::default();
/// let pairs = detect_glue_faces(&box1, &box2, &config);
/// ```
pub fn detect_glue_faces(
    brep_a: &BRep,
    brep_b: &BRep,
    config: &GlueConfig,
) -> Vec<GlueFacePair> {
    let mut result = Vec::new();

    // Build caches for both BReps
    let cell_size = config.face_tolerance * 10.0;
    let mut cache_a = GlueFaceCache::new();
    let mut cache_b = GlueFaceCache::new();
    cache_a.build(brep_a, cell_size);
    cache_b.build(brep_b, cell_size);

    // Get face counts
    let faces_a: Vec<(usize, DVec3, DVec3, f64)> = brep_a.solids.iter()
        .flat_map(|s| s.shells.iter())
        .flat_map(|sh| sh.faces.iter().enumerate())
        .enumerate()
        .map(|(idx, (_, face))| {
            let center = cache_a.face_centers.get(idx).copied().unwrap_or(DVec3::ZERO);
            let normal = face.normal;
            let area = cache_a.face_areas.get(idx).copied().unwrap_or(0.0);
            (idx, center, normal, area)
        })
        .collect();

    let faces_b: Vec<(usize, DVec3, DVec3, f64)> = brep_b.solids.iter()
        .flat_map(|s| s.shells.iter())
        .flat_map(|sh| sh.faces.iter().enumerate())
        .enumerate()
        .map(|(idx, (_, face))| {
            let center = cache_b.face_centers.get(idx).copied().unwrap_or(DVec3::ZERO);
            let normal = face.normal;
            let area = cache_b.face_areas.get(idx).copied().unwrap_or(0.0);
            (idx, center, normal, area)
        })
        .collect();

    // Early normal filter threshold
    let normal_threshold = -0.95;

    for (idx_a, center_a, normal_a, area_a) in &faces_a {
        // Use geometric hash to find nearby faces in B
        let nearby_faces = if config.use_geometric_hash {
            cache_b.get_nearby_faces(*center_a, cell_size)
        } else {
            faces_b.iter().map(|(idx, _, _, _)| *idx).collect()
        };

        for idx_b in nearby_faces {
            let (_, center_b, normal_b, area_b) = &faces_b.get(idx_b).unwrap_or(&(0, DVec3::ZERO, DVec3::ZERO, 0.0));

            // Early normal filter: skip if normals are not anti-parallel
            if config.early_normal_filter {
                let na_len2 = normal_a.length_squared();
                let nb_len2 = normal_b.length_squared();
                if na_len2 > 1e-12 && nb_len2 > 1e-12 {
                    let na = *normal_a / na_len2.sqrt();
                    let nb = *normal_b / nb_len2.sqrt();
                    if na.dot(nb) > normal_threshold {
                        continue;
                    }
                }
            }

            // Check center proximity
            let center_dist = (*center_a - *center_b).length();
            if center_dist > config.face_tolerance * 10.0 {
                continue;
            }

            // Compute match quality
            let normal_match = {
                let na_len2 = normal_a.length_squared();
                let nb_len2 = normal_b.length_squared();
                if na_len2 > 1e-12 && nb_len2 > 1e-12 {
                    let na = *normal_a / na_len2.sqrt();
                    let nb = *normal_b / nb_len2.sqrt();
                    // For glue, normals should be anti-parallel
                    (-na.dot(nb)).max(0.0)
                } else {
                    0.0
                }
            };

            let center_match = {
                let max_dist = config.face_tolerance * 10.0;
                if max_dist > 0.0 {
                    (1.0 - center_dist / max_dist).max(0.0)
                } else {
                    1.0
                }
            };

            let area_match = {
                let max_area = area_a.max(*area_b);
                let min_area = area_a.min(*area_b);
                if max_area > 0.0 {
                    min_area / max_area
                } else {
                    1.0
                }
            };

            let match_quality = (normal_match * 0.4 + center_match * 0.3 + area_match * 0.3).min(1.0);

            // Only include pairs with reasonable match quality
            if match_quality >= 0.5 {
                result.push(GlueFacePair {
                    face_a: *idx_a,
                    face_b: idx_b,
                    match_quality,
                    shared_area: area_a.min(*area_b),
                });
            }
        }
    }

    // Sort by match quality (highest first)
    result.sort_by(|a, b| {
        b.match_quality.partial_cmp(&a.match_quality).unwrap_or(std::cmp::Ordering::Equal)
    });

    result
}

/// Apply glue optimization to pave filler.
///
/// This function configures a PaveFiller to use pre-detected glue face pairs,
/// enabling it to skip expensive interference computations for coincident faces.
///
/// # Arguments
///
/// * `filler` - The PaveFiller to optimize.
/// * `glue_pairs` - Pre-detected glue face pairs.
///
/// # Example
///
/// ```
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use rcad_algorithms::bopds::ds::DS;
/// use rcad_algorithms::pave_filler::PaveFiller;
/// use rcad_algorithms::builder::{GlueConfig, detect_glue_faces, apply_glue_optimization};
///
/// let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let box2 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
///
/// let config = GlueConfig::default();
/// let pairs = detect_glue_faces(&box1, &box2, &config);
///
/// let mut ds = DS::new(&box1, &box2);
/// let mut filler = PaveFiller::new(&mut ds);
/// apply_glue_optimization(&mut filler, &pairs);
/// ```
pub fn apply_glue_optimization(
    filler: &mut crate::pave_filler::PaveFiller,
    glue_pairs: &[GlueFacePair],
) {
    if glue_pairs.is_empty() {
        return;
    }

    // Use the tolerance from the best match
    let best_pair = glue_pairs.iter()
        .max_by(|a, b| {
            a.match_quality.partial_cmp(&b.match_quality).unwrap_or(std::cmp::Ordering::Equal)
        });

    if let Some(pair) = best_pair {
        // Estimate tolerance from match quality
        let tolerance = if pair.match_quality > 0.99 {
            TOLERANCE_ABS
        } else if pair.match_quality > 0.9 {
            TOLERANCE_ABS * 10.0
        } else {
            TOLERANCE_ABS * 100.0
        };

        filler.configure_glue(true, tolerance);
    }
}

/// Compute adaptive glue tolerance based on geometry characteristics.
///
/// Analyzes the input BReps and computes an appropriate glue tolerance
/// based on the minimum feature size, face area distribution, and
/// edge length distribution.
///
/// # Arguments
///
/// * `brep_a` - First BRep shape.
/// * `brep_b` - Second BRep shape.
/// * `base_tolerance` - Base tolerance to start with.
///
/// # Returns
///
/// The computed adaptive glue tolerance.
pub fn compute_adaptive_glue_tolerance(
    brep_a: &BRep,
    brep_b: &BRep,
    base_tolerance: f64,
) -> f64 {
    let mut min_feature_size = f64::INFINITY;

    // Analyze edge lengths
    for edge in &brep_a.edges {
        if edge.start < brep_a.vertices.len() && edge.end < brep_a.vertices.len() {
            let p1 = brep_a.vertices[edge.start].point;
            let p2 = brep_a.vertices[edge.end].point;
            let length = (p2 - p1).length();
            if length > 1e-10 {
                min_feature_size = min_feature_size.min(length);
            }
        }
    }
    for edge in &brep_b.edges {
        if edge.start < brep_b.vertices.len() && edge.end < brep_b.vertices.len() {
            let p1 = brep_b.vertices[edge.start].point;
            let p2 = brep_b.vertices[edge.end].point;
            let length = (p2 - p1).length();
            if length > 1e-10 {
                min_feature_size = min_feature_size.min(length);
            }
        }
    }

    // Analyze face areas (approximate from bounding box)
    for solid in &brep_a.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let mut min_pt = DVec3::splat(f64::INFINITY);
                let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
                for we in &face.outer_wire.edges {
                    if we.idx < brep_a.edges.len() {
                        let edge = &brep_a.edges[we.idx];
                        if edge.start < brep_a.vertices.len() {
                            let p = brep_a.vertices[edge.start].point;
                            min_pt = min_pt.min(p);
                            max_pt = max_pt.max(p);
                        }
                        if edge.end < brep_a.vertices.len() {
                            let p = brep_a.vertices[edge.end].point;
                            min_pt = min_pt.min(p);
                            max_pt = max_pt.max(p);
                        }
                    }
                }
                let diag = max_pt - min_pt;
                let size = diag.x.min(diag.y).min(diag.z);
                if size > 1e-10 {
                    min_feature_size = min_feature_size.min(size);
                }
            }
        }
    }
    for solid in &brep_b.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let mut min_pt = DVec3::splat(f64::INFINITY);
                let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
                for we in &face.outer_wire.edges {
                    if we.idx < brep_b.edges.len() {
                        let edge = &brep_b.edges[we.idx];
                        if edge.start < brep_b.vertices.len() {
                            let p = brep_b.vertices[edge.start].point;
                            min_pt = min_pt.min(p);
                            max_pt = max_pt.max(p);
                        }
                        if edge.end < brep_b.vertices.len() {
                            let p = brep_b.vertices[edge.end].point;
                            min_pt = min_pt.min(p);
                            max_pt = max_pt.max(p);
                        }
                    }
                }
                let diag = max_pt - min_pt;
                let size = diag.x.min(diag.y).min(diag.z);
                if size > 1e-10 {
                    min_feature_size = min_feature_size.min(size);
                }
            }
        }
    }

    // Compute adaptive tolerance
    let adaptive_tol = if min_feature_size.is_finite() && min_feature_size > 0.0 {
        // Use a fraction of minimum feature size, but at least base tolerance
        let feature_based = min_feature_size * 0.01;
        base_tolerance.max(feature_based).min(min_feature_size * 0.1)
    } else {
        base_tolerance
    };

    adaptive_tol.max(TOLERANCE_ABS)
}

#[cfg(test)]
mod glue_tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;
    use rcad_kernel::geom::{SphericalSurface, CylindricalSurface, ConicalSurface, ToroidalSurface};
    use glam::DAffine3;

    fn unit_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        })
    }

    #[test]
    fn test_glue_config_default() {
        let config = GlueConfig::default();
        assert_eq!(config.face_tolerance, TOLERANCE_ABS);
        assert_eq!(config.edge_tolerance, TOLERANCE_ABS);
        assert!(config.use_geometric_hash);
        assert!(config.early_normal_filter);
    }

    #[test]
    fn test_detect_glue_faces_no_overlap() {
        let box1 = unit_box();
        let box2 = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }).transformed(DAffine3::from_translation(DVec3::new(10.0, 0.0, 0.0)));

        let config = GlueConfig::default();
        let pairs = detect_glue_faces(&box1, &box2, &config);

        // No overlapping faces
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_detect_glue_faces_touching() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        // Translate box2 to touch box1 at y=1 face
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let config = GlueConfig::default();
        let pairs = detect_glue_faces(&box1, &box2, &config);

        // Should detect at least one coincident face pair
        assert!(!pairs.is_empty());

        // Match quality should be high for exact match
        assert!(pairs[0].match_quality > 0.9);
    }

    #[test]
    fn test_detect_glue_faces_with_tolerance() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        // Slight offset - faces are near but not exactly coincident
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0 + 1e-7, 0.0)));

        let config = GlueConfig {
            face_tolerance: 1e-5,
            ..Default::default()
        };
        let pairs = detect_glue_faces(&box1, &box2, &config);

        // Should still detect coincident faces with relaxed tolerance
        assert!(!pairs.is_empty());
    }

    #[test]
    fn test_glue_face_pair_structure() {
        let pair = GlueFacePair {
            face_a: 0,
            face_b: 1,
            match_quality: 0.95,
            shared_area: 1.0,
        };

        assert_eq!(pair.face_a, 0);
        assert_eq!(pair.face_b, 1);
        assert!((pair.match_quality - 0.95).abs() < 1e-10);
        assert!((pair.shared_area - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_glue_face_cache_build() {
        let box1 = unit_box();
        let mut cache = GlueFaceCache::new();
        cache.build(&box1, 1.0);

        // Should have cached 6 faces (box has 6 faces)
        assert_eq!(cache.face_centers.len(), 6);
        assert_eq!(cache.face_normals.len(), 6);
        assert_eq!(cache.face_areas.len(), 6);

        // Spatial hash should not be empty
        assert!(!cache.spatial_hash.is_empty());
    }

    #[test]
    fn test_glue_face_cache_nearby_faces() {
        let box1 = unit_box();
        let mut cache = GlueFaceCache::new();
        cache.build(&box1, 1.0);

        // Get nearby faces for the center of the box
        let nearby = cache.get_nearby_faces(DVec3::new(0.5, 0.5, 0.5), 1.0);

        // Should find at least some faces
        assert!(!nearby.is_empty());
    }

    #[test]
    fn test_compute_adaptive_glue_tolerance() {
        let box1 = unit_box();
        let box2 = unit_box();

        let tolerance = compute_adaptive_glue_tolerance(&box1, &box2, 1e-6);

        // Tolerance should be reasonable
        assert!(tolerance >= TOLERANCE_ABS);
        assert!(tolerance < 1.0); // Should be much smaller than box size
    }

    #[test]
    fn test_early_normal_filter_disabled() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let config = GlueConfig {
            early_normal_filter: false,
            ..Default::default()
        };
        let pairs = detect_glue_faces(&box1, &box2, &config);

        // Should still detect coincident faces
        assert!(!pairs.is_empty());
    }

    #[test]
    fn test_geometric_hash_disabled() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let config = GlueConfig {
            use_geometric_hash: false,
            ..Default::default()
        };
        let pairs = detect_glue_faces(&box1, &box2, &config);

        // Should still detect coincident faces
        assert!(!pairs.is_empty());
    }

    #[test]
    fn test_match_quality_ordering() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        // Perfect match
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let mut box3 = unit_box();
        // Slight rotation - not as good a match
        box3.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));
        box3.apply_transform(DAffine3::from_rotation_z(0.001));

        let config = GlueConfig::default();

        let pairs_exact = detect_glue_faces(&box1, &box2, &config);
        let pairs_rotated = detect_glue_faces(&box1, &box3, &config);

        // Exact match should have higher quality
        if !pairs_exact.is_empty() && !pairs_rotated.is_empty() {
            assert!(pairs_exact[0].match_quality >= pairs_rotated[0].match_quality);
        }
    }

    #[test]
    fn test_shared_area_estimation() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let config = GlueConfig::default();
        let pairs = detect_glue_faces(&box1, &box2, &config);

        // Shared area should be approximately 1.0 (unit square face)
        assert!(!pairs.is_empty());
        assert!(pairs[0].shared_area > 0.1);
    }

    #[test]
    fn test_multiple_face_pairs() {
        // Create two boxes that share multiple faces (impossible in real geometry,
        // but tests the algorithm)
        let box1 = unit_box();
        let mut box2 = unit_box();
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let config = GlueConfig::default();
        let pairs = detect_glue_faces(&box1, &box2, &config);

        // Should detect exactly one face pair (the touching faces)
        assert!(!pairs.is_empty());
        // All pairs should have valid indices
        for pair in &pairs {
            assert!(pair.face_a < 6); // Box has 6 faces
            assert!(pair.face_b < 6);
        }
    }

    #[test]
    fn test_compatibility_cache() {
        let mut cache = GlueFaceCache::new();

        // Initially no cached value
        assert!(cache.get_compatibility(0, 1).is_none());

        // Set and retrieve
        cache.set_compatibility(0, 1, true);
        assert_eq!(cache.get_compatibility(0, 1), Some(true));
        assert_eq!(cache.get_compatibility(1, 0), Some(true)); // Symmetric

        cache.set_compatibility(0, 1, false);
        assert_eq!(cache.get_compatibility(0, 1), Some(false));
    }

    #[test]
    fn test_glue_config_custom_values() {
        let config = GlueConfig {
            face_tolerance: 1e-5,
            edge_tolerance: 2e-5,
            use_geometric_hash: false,
            early_normal_filter: false,
        };

        assert!((config.face_tolerance - 1e-5).abs() < 1e-12);
        assert!((config.edge_tolerance - 2e-5).abs() < 1e-12);
        assert!(!config.use_geometric_hash);
        assert!(!config.early_normal_filter);
    }

    #[test]
    fn split_uv_polygon_detects_seam_crossing_on_cylinder() {
        // UV polygon that crosses the U=0/2π seam on a cylinder
        // This is a quad that wraps around the seam:
        // - Right side: u ≈ 5.5 (near 2π)
        // - Left side: u ≈ 0.5 (near 0)
        let period = std::f64::consts::TAU; // ≈ 6.283
        let uv_polygon = vec![
            DVec2::new(5.5, 0.0),  // Near 2π
            DVec2::new(0.5, 0.0),  // Near 0
            DVec2::new(0.5, 1.0),
            DVec2::new(5.5, 1.0),
        ];

        let result = split_uv_polygon_at_seam(&uv_polygon, period);

        // Should split into two polygons
        assert_eq!(result.len(), 2, "Seam crossing should split polygon");

        // Each output polygon must have at least 3 vertices
        for (i, poly) in result.iter().enumerate() {
            assert!(
                poly.len() >= 3,
                "Output polygon {} has only {} vertices (need >= 3)",
                i,
                poly.len()
            );
        }

        // No output polygon should cross the seam
        for (i, poly) in result.iter().enumerate() {
            for j in 0..poly.len() {
                let k = (j + 1) % poly.len();
                let du = poly[k].x - poly[j].x;
                assert!(
                    du.abs() < period * 0.5,
                    "Output polygon {} still crosses seam: du = {} between vertices {} and {}",
                    i,
                    du,
                    j,
                    k
                );
            }
        }

        // Verify specific coordinates: each polygon should contain seam intersection points
        // The original polygon has edges that cross the seam at v=0 and v=1
        // Output polygons should have intersection points at u=0 or u=period

        // Find the right-side polygon (u values near 5.5)
        let right_poly = result
            .iter()
            .find(|p| p.iter().any(|v| v.x > period * 0.5))
            .expect("Should have a polygon with high u values");
        // Find the left-side polygon (u values near 0.5)
        let left_poly = result
            .iter()
            .find(|p| p.iter().any(|v| v.x < period * 0.5))
            .expect("Should have a polygon with low u values");

        // Right polygon should have vertices with u near 5.5 and seam points
        let has_high_u = right_poly.iter().any(|v| (v.x - 5.5).abs() < 0.01);
        assert!(has_high_u, "Right polygon should contain original high-u vertices");

        // Left polygon should have vertices with u near 0.5 and seam points
        let has_low_u = left_poly.iter().any(|v| (v.x - 0.5).abs() < 0.01);
        assert!(has_low_u, "Left polygon should contain original low-u vertices");

        // Each polygon should have seam intersection points
        // (either at u=0 or u=period, both representing the same physical location)
        fn near_seam(u: f64, period: f64) -> bool {
            u.abs() < 0.01 || (u - period).abs() < 0.01
        }

        assert!(
            right_poly.iter().any(|v| near_seam(v.x, period)),
            "Right polygon should have a seam intersection point"
        );
        assert!(
            left_poly.iter().any(|v| near_seam(v.x, period)),
            "Left polygon should have a seam intersection point"
        );
    }

    #[test]
    fn split_uv_polygon_no_crossing_returns_original() {
        // Polygon that doesn't cross the seam
        let period = std::f64::consts::TAU;
        let uv_polygon = vec![
            DVec2::new(1.0, 0.0),
            DVec2::new(2.0, 0.0),
            DVec2::new(2.0, 1.0),
            DVec2::new(1.0, 1.0),
        ];

        let result = split_uv_polygon_at_seam(&uv_polygon, period);

        assert_eq!(result.len(), 1, "No seam crossing should return one polygon");
        assert_eq!(result[0].len(), 4, "Original polygon should be unchanged");
    }

    #[test]
    fn split_uv_polygon_degenerate_input() {
        let period = std::f64::consts::TAU;

        // Less than 3 vertices
        let two_vertices = vec![DVec2::new(1.0, 0.0), DVec2::new(2.0, 0.0)];
        let result = split_uv_polygon_at_seam(&two_vertices, period);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].len(), 2);

        // Empty input
        let empty: Vec<DVec2> = vec![];
        let result = split_uv_polygon_at_seam(&empty, period);
        assert_eq!(result.len(), 1);
        assert!(result[0].is_empty());
    }

    // =====================================================
    // Track A: Periodic Surface Seam Enhancement Tests
    // =====================================================

    // --- A1: Enhanced degenerate UV polygon handling tests ---

    #[test]
    fn test_handle_degenerate_uv_polygon_sphere_pole_cap() {
        // UV polygon that represents a small cap near the north pole of a sphere
        // All vertices collapse toward v=0 (north pole)
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            radius: 1.0,
            axis: DVec3::Y,
        };
        let surface = Surface3::Sphere(sphere);

        // Small triangle near north pole (v ≈ 0)
        let uv_polygon = vec![
            DVec2::new(0.0, 0.001),
            DVec2::new(std::f64::consts::FRAC_PI_2, 0.001),
            DVec2::new(std::f64::consts::PI, 0.001),
        ];

        let result = handle_degenerate_uv_polygon(&uv_polygon, &surface);

        // Should produce valid 3D boundary
        assert!(!result.is_empty(), "Should produce non-empty boundary");

        // All points should be valid (no NaN)
        for pt in &result {
            assert!(pt.x.is_finite(), "Point x should be finite");
            assert!(pt.y.is_finite(), "Point y should be finite");
            assert!(pt.z.is_finite(), "Point z should be finite");
        }

        // Should include pole point since all vertices are near pole
        let north_pole = sphere.center + sphere.axis * sphere.radius;
        let has_pole = result.iter().any(|pt| (*pt - north_pole).length() < 0.1);
        assert!(has_pole, "Should include pole point for collapsed vertices");
    }

    #[test]
    fn test_handle_degenerate_uv_polygon_sphere_south_pole_cap() {
        // UV polygon near south pole (v ≈ π)
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            radius: 1.0,
            axis: DVec3::Y,
        };
        let surface = Surface3::Sphere(sphere);

        // Small triangle near south pole (v ≈ π)
        let uv_polygon = vec![
            DVec2::new(0.0, std::f64::consts::PI - 0.001),
            DVec2::new(std::f64::consts::FRAC_PI_2, std::f64::consts::PI - 0.001),
            DVec2::new(std::f64::consts::PI, std::f64::consts::PI - 0.001),
        ];

        let result = handle_degenerate_uv_polygon(&uv_polygon, &surface);

        // Should produce valid 3D boundary
        assert!(!result.is_empty(), "Should produce non-empty boundary");

        // Should include south pole point
        let south_pole = sphere.center - sphere.axis * sphere.radius;
        let has_pole = result.iter().any(|pt| (*pt - south_pole).length() < 0.1);
        assert!(has_pole, "Should include south pole point for collapsed vertices");
    }

    #[test]
    fn test_handle_degenerate_uv_polygon_cone_apex() {
        // UV polygon that collapses toward cone apex (v=0)
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 0.0, // Reference radius at apex
            half_angle_rad: std::f64::consts::FRAC_PI_4,
        };
        let surface = Surface3::Cone(cone);

        // Small triangle near apex (v ≈ 0)
        let uv_polygon = vec![
            DVec2::new(0.0, 0.001),
            DVec2::new(std::f64::consts::FRAC_PI_2, 0.001),
            DVec2::new(std::f64::consts::PI, 0.001),
        ];

        let result = handle_degenerate_uv_polygon(&uv_polygon, &surface);

        // Should produce valid 3D boundary
        assert!(!result.is_empty(), "Should produce non-empty boundary");

        // All points should be valid (no NaN)
        for pt in &result {
            assert!(pt.x.is_finite(), "Point x should be finite");
            assert!(pt.y.is_finite(), "Point y should be finite");
            assert!(pt.z.is_finite(), "Point z should be finite");
        }

        // Should include apex point
        let apex = cone.apex_point();
        let has_apex = result.iter().any(|pt| (*pt - apex).length() < 0.1);
        assert!(has_apex, "Should include apex point for collapsed vertices");
    }

    #[test]
    fn test_handle_degenerate_uv_polygon_sphere_triangular_pole_cap() {
        // A triangular UV region that includes the pole, simulating a spherical triangle
        // with one vertex at the pole
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            radius: 1.0,
            axis: DVec3::Y,
        };
        let surface = Surface3::Sphere(sphere);

        // Triangle with pole at one vertex
        // u=0, v=0 is the pole, other vertices at larger v
        let uv_polygon = vec![
            DVec2::new(0.0, 0.0), // At pole
            DVec2::new(0.0, 0.5), // Away from pole
            DVec2::new(std::f64::consts::FRAC_PI_2, 0.5), // Away from pole
        ];

        let result = handle_degenerate_uv_polygon(&uv_polygon, &surface);

        // Should produce valid 3D boundary with at least 2 distinct points
        assert!(result.len() >= 2, "Should produce at least 2 boundary points");

        // All points should be valid (no NaN)
        for pt in &result {
            assert!(pt.x.is_finite(), "Point x should be finite");
            assert!(pt.y.is_finite(), "Point y should be finite");
            assert!(pt.z.is_finite(), "Point z should be finite");
        }
    }

    // --- A2: Edge splitting at periodic seam tests ---

    #[test]
    fn test_split_edge_at_periodic_seam_cylinder() {
        // Edge that crosses U=0/2π boundary on cylinder
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        };

        // Edge from u near 2π to u near 0
        let start_uv = DVec2::new(std::f64::consts::TAU - 0.1, 0.5);
        let end_uv = DVec2::new(0.1, 0.5);

        let result = split_edge_at_periodic_seam(start_uv, end_uv, &Surface3::Cylinder(cylinder));

        // Should return two segments
        assert!(result.is_some(), "Should detect seam crossing");
        let segments = result.unwrap();
        assert_eq!(segments.len(), 2, "Should split into two segments");

        // Each segment should have start and end UV
        for (i, seg) in segments.iter().enumerate() {
            assert_eq!(seg.len(), 2, "Segment {} should have 2 points", i);
        }

        // First segment should end at seam
        assert!(
            segments[0][1].x.abs() < 0.01 || (segments[0][1].x - std::f64::consts::TAU).abs() < 0.01,
            "First segment should end at seam"
        );

        // Second segment should start at seam
        assert!(
            segments[1][0].x.abs() < 0.01 || (segments[1][0].x - std::f64::consts::TAU).abs() < 0.01,
            "Second segment should start at seam"
        );
    }

    #[test]
    fn test_split_edge_at_periodic_seam_no_crossing() {
        // Edge that doesn't cross seam
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        };

        let start_uv = DVec2::new(1.0, 0.5);
        let end_uv = DVec2::new(2.0, 0.5);

        let result = split_edge_at_periodic_seam(start_uv, end_uv, &Surface3::Cylinder(cylinder));

        // Should return None (no splitting needed)
        assert!(result.is_none(), "Should not split edge that doesn't cross seam");
    }

    #[test]
    fn test_split_edge_at_periodic_seam_sphere() {
        // Edge crossing U=0/2π boundary on sphere
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            radius: 1.0,
            axis: DVec3::Y,
        };

        let start_uv = DVec2::new(std::f64::consts::TAU - 0.1, 1.0);
        let end_uv = DVec2::new(0.1, 1.0);

        let result = split_edge_at_periodic_seam(start_uv, end_uv, &Surface3::Sphere(sphere));

        assert!(result.is_some(), "Should detect seam crossing on sphere");
        let segments = result.unwrap();
        assert_eq!(segments.len(), 2, "Should split into two segments");
    }

    // --- A3: Torus double periodicity tests ---

    #[test]
    fn test_split_uv_polygon_torus_u_period() {
        // UV polygon on torus that crosses U seam only
        let period = std::f64::consts::TAU;
        let uv_polygon = vec![
            DVec2::new(5.5, 0.5), // Near U=2π
            DVec2::new(0.5, 0.5), // Near U=0
            DVec2::new(0.5, 1.5),
            DVec2::new(5.5, 1.5),
        ];

        let result = split_uv_polygon_at_seam(&uv_polygon, period);

        // Should split into two polygons
        assert_eq!(result.len(), 2, "Should split torus polygon at U seam");

        // Each polygon should not cross U seam
        for poly in &result {
            for j in 0..poly.len() {
                let k = (j + 1) % poly.len();
                let du = poly[k].x - poly[j].x;
                assert!(
                    du.abs() < period * 0.5,
                    "Output polygon should not cross U seam"
                );
            }
        }
    }

    #[test]
    fn test_split_uv_polygon_torus_double_period() {
        // UV polygon on torus that crosses both U and V seams
        // This is a complex case where the polygon wraps around both directions
        let period = std::f64::consts::TAU;

        // Polygon that spans nearly full U range and crosses V seam
        let uv_polygon = vec![
            DVec2::new(0.1, 5.5), // V near 2π
            DVec2::new(5.9, 5.5),
            DVec2::new(5.9, 0.5), // V near 0
            DVec2::new(0.1, 0.5),
        ];

        // Use double periodic splitting
        let result = split_uv_polygon_torus_double(&uv_polygon, period);

        // Should produce multiple non-crossing polygons
        assert!(!result.is_empty(), "Should produce output polygons");

        // Each polygon should not cross U or V seams
        for poly in &result {
            assert!(poly.len() >= 3, "Polygon should have at least 3 vertices");

            for j in 0..poly.len() {
                let k = (j + 1) % poly.len();
                let du = poly[k].x - poly[j].x;
                let dv = poly[k].y - poly[j].y;
                assert!(
                    du.abs() < period * 0.5,
                    "Output polygon should not cross U seam"
                );
                assert!(
                    dv.abs() < period * 0.5,
                    "Output polygon should not cross V seam"
                );
            }
        }
    }

    #[test]
    fn test_handle_degenerate_uv_polygon_non_degenerate() {
        // Normal UV polygon on sphere (no degenerate points)
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            radius: 1.0,
            axis: DVec3::Y,
        };
        let surface = Surface3::Sphere(sphere);

        // Rectangle away from poles
        let uv_polygon = vec![
            DVec2::new(0.0, 1.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(1.0, 2.0),
            DVec2::new(0.0, 2.0),
        ];

        let result = handle_degenerate_uv_polygon(&uv_polygon, &surface);

        // Should produce same number of points as input
        assert_eq!(result.len(), uv_polygon.len(), "Non-degenerate should map 1:1");

        // All points should be on sphere surface
        for pt in &result {
            let dist = pt.length();
            assert!(
                (dist - sphere.radius).abs() < 0.001,
                "Point should be on sphere surface"
            );
        }
    }
}
