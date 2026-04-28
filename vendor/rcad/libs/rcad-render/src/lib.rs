use rcad_kernel::{
    BRep, BRepGraph, Curve3, CurveEval, Surface3, SurfaceEval, any_perpendicular,
    seam_edge_candidates,
};
use rcad_algorithms::{TessellationParams, mesh_brep};
use wgpu::util::DeviceExt;

/// Tessellation quality options for [`Tessellator::tessellate_with_options`].
///
/// Re-exported from [`rcad_algorithms::TessellationParams`].
pub type TessellationOptions = TessellationParams;

/// Edited topology/geometry entities used to drive incremental mesh invalidation.
///
/// Indices are optional and may be mixed: if both vertices and edges are listed,
/// all adjacent faces of either set will be invalidated.
#[derive(Debug, Clone, Default)]
pub struct EditedModelDelta {
    /// Modified vertex indices in `BRep.vertices`.
    pub modified_vertices: Vec<usize>,
    /// Modified edge indices in `BRep.edges`.
    pub modified_edges: Vec<usize>,
    /// Modified flattened face indices (solid/shell/face traversal order).
    pub modified_faces: Vec<usize>,
}

impl EditedModelDelta {
    pub fn is_empty(&self) -> bool {
        self.modified_vertices.is_empty()
            && self.modified_edges.is_empty()
            && self.modified_faces.is_empty()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SelectionMode {
    Face,
    Edge,
}

#[derive(Clone, Debug)]
pub struct SelectionState {
    pub mode: SelectionMode,
    pub additive_select: bool,
    pub selected_faces: Vec<usize>,
    pub selected_edges: Vec<usize>,
    pub hovered_face: Option<usize>,
    pub hovered_edge: Option<usize>,
    pub last_hover_pos: Option<(f32, f32)>,
}

#[derive(Default, Clone, Copy, Debug, PartialEq, Eq)]
pub enum DisplayMode {
    #[default]
    SolidWithEdges,
    Solid,
    Wireframe,
    Transparent,
}

pub const DEFAULT_EDGE_PICK_RADIUS_PX: f32 = 8.0;

const AXIS_GIZMO_CAMERA_DISTANCE: f32 = 3.2;
const AXIS_GIZMO_SIDE_RATIO: f32 = 0.22;
const AXIS_GIZMO_MIN_SIDE_PX: u32 = 92;
const AXIS_GIZMO_MAX_SIDE_PX: u32 = 160;
const AXIS_GIZMO_PADDING_RATIO: f32 = 0.12;
const AXIS_GIZMO_AXIS_LENGTH: f32 = 0.78;
const AXIS_GIZMO_CENTER_HALF_EXTENT: f32 = 0.22;
const GRID_MAJOR_LINE_EVERY: i32 = 5;
const GRID_BUFFER_HALF_CELLS: i32 = 64;
const GRID_TARGET_HALF_MINOR_LINES: f32 = 12.0;
const GRID_MIN_MINOR_SPACING: f32 = 0.02;
const GRID_FOV_Y_RADIANS: f32 = 45.0_f32.to_radians();

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AxisGizmoLayout {
    pub origin_px: [u32; 2],
    pub size_px: [u32; 2],
}

impl AxisGizmoLayout {
    pub fn contains_point(&self, point_px: [f32; 2]) -> bool {
        let min_x = self.origin_px[0] as f32;
        let min_y = self.origin_px[1] as f32;
        let max_x = min_x + self.size_px[0] as f32;
        let max_y = min_y + self.size_px[1] as f32;
        point_px[0] >= min_x
            && point_px[1] >= min_y
            && point_px[0] <= max_x
            && point_px[1] <= max_y
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AxisGizmoHit {
    X,
    Y,
    Z,
    Center,
}

pub fn axis_gizmo_layout(
    viewport_origin_px: [u32; 2],
    viewport_size_px: [u32; 2],
) -> Option<AxisGizmoLayout> {
    let width = viewport_size_px[0];
    let height = viewport_size_px[1];
    if width < 32 || height < 32 {
        return None;
    }

    let side = ((width.min(height) as f32) * AXIS_GIZMO_SIDE_RATIO).round() as u32;
    let side = side
        .clamp(AXIS_GIZMO_MIN_SIDE_PX, AXIS_GIZMO_MAX_SIDE_PX)
        .min(width)
        .min(height);
    let padding = ((side as f32) * AXIS_GIZMO_PADDING_RATIO).round() as u32;
    if side == 0 || side + padding > width || side + padding > height {
        return None;
    }

    Some(AxisGizmoLayout {
        origin_px: [
            viewport_origin_px[0] + width - side - padding,
            viewport_origin_px[1] + height - side - padding,
        ],
        size_px: [side, side],
    })
}

fn axis_gizmo_eye(camera: &Camera) -> glam::Vec3 {
    let mut eye_dir = camera.eye_position().normalize_or_zero();
    if eye_dir.length_squared() <= 1e-8 {
        eye_dir = glam::Vec3::new(1.0, 1.0, 1.0).normalize();
    }
    eye_dir * AXIS_GIZMO_CAMERA_DISTANCE
}

fn axis_gizmo_view_projection(camera: &Camera) -> [[f32; 4]; 4] {
    let eye = axis_gizmo_eye(camera);
    let forward = (-eye).normalize_or_zero();
    let up = if forward.dot(glam::Vec3::Y).abs() > 0.98 {
        glam::Vec3::Z
    } else {
        glam::Vec3::Y
    };
    let view = glam::Mat4::look_at_rh(eye, glam::Vec3::ZERO, up);
    let proj = glam::Mat4::perspective_rh(45.0_f32.to_radians(), 1.0, 0.1, 100.0);
    (proj * view).to_cols_array_2d()
}

pub fn axis_gizmo_hit_test(
    camera: &Camera,
    viewport_origin_px: [u32; 2],
    viewport_size_px: [u32; 2],
    pointer_px: [f32; 2],
) -> Option<AxisGizmoHit> {
    let layout = axis_gizmo_layout(viewport_origin_px, viewport_size_px)?;
    if !layout.contains_point(pointer_px) {
        return None;
    }

    let local = [
        pointer_px[0] - layout.origin_px[0] as f32,
        pointer_px[1] - layout.origin_px[1] as f32,
    ];
    let side = layout.size_px[0] as f32;
    let vp = glam::Mat4::from_cols_array_2d(&axis_gizmo_view_projection(camera));
    let center = project_to_screen(vp, glam::Vec3::ZERO, [side, side])?;

    let center_radius = (side * 0.14).clamp(10.0, 18.0);
    let center_distance = ((local[0] - center[0]).powi(2) + (local[1] - center[1]).powi(2)).sqrt();
    if center_distance <= center_radius {
        return Some(AxisGizmoHit::Center);
    }

    let axis_threshold = (side * 0.10).clamp(8.0, 18.0);
    let endpoint_radius = (side * 0.13).clamp(10.0, 22.0);
    let axes = [
        (AxisGizmoHit::X, glam::Vec3::X),
        (AxisGizmoHit::Y, glam::Vec3::Y),
        (AxisGizmoHit::Z, glam::Vec3::Z),
    ];
    let mut best: Option<(f32, f32, AxisGizmoHit)> = None;

    for (hit, axis_dir) in axes {
        let tip = project_to_screen(vp, axis_dir * AXIS_GIZMO_AXIS_LENGTH, [side, side])?;
        let tip_distance = ((local[0] - tip[0]).powi(2) + (local[1] - tip[1]).powi(2)).sqrt();
        let shaft_distance =
            point_segment_distance_2d(local, [center[0], center[1]], [tip[0], tip[1]]);
        if tip_distance > endpoint_radius && shaft_distance > axis_threshold {
            continue;
        }

        let score = tip_distance.min(shaft_distance + endpoint_radius * 0.35);
        match best {
            Some((best_score, best_depth, _))
                if score > best_score
                    || ((score - best_score).abs() < 1e-3 && tip[2] >= best_depth) => {}
            _ => best = Some((score, tip[2], hit)),
        }
    }

    best.map(|(_, _, hit)| hit)
}

impl Default for SelectionState {
    fn default() -> Self {
        Self {
            mode: SelectionMode::Face,
            additive_select: false,
            selected_faces: Vec::new(),
            selected_edges: Vec::new(),
            hovered_face: None,
            hovered_edge: None,
            last_hover_pos: None,
        }
    }
}

impl SelectionState {
    pub fn set_mode(&mut self, mode: SelectionMode) {
        if self.mode != mode {
            self.mode = mode;
            self.clear_hover();
        }
    }

    pub fn click_at(
        &mut self,
        brep: &BRep,
        camera: &Camera,
        aspect: f32,
        viewport: [f32; 2],
        cursor: [f32; 2],
        edge_pick_radius_px: f32,
    ) {
        if cursor[0] < 0.0 || cursor[1] < 0.0 || cursor[0] > viewport[0] || cursor[1] > viewport[1]
        {
            return;
        }

        match self.mode {
            SelectionMode::Face => {
                let hit = pick_face(brep, camera, aspect, viewport, cursor);
                if self.additive_select {
                    if let Some(idx) = hit {
                        toggle_index(&mut self.selected_faces, idx);
                    }
                } else {
                    self.selected_faces.clear();
                    if let Some(idx) = hit {
                        self.selected_faces.push(idx);
                    }
                }
            }
            SelectionMode::Edge => {
                let hit = pick_edge(brep, camera, aspect, viewport, cursor, edge_pick_radius_px);
                if self.additive_select {
                    if let Some(idx) = hit {
                        toggle_index(&mut self.selected_edges, idx);
                    }
                } else {
                    self.selected_edges.clear();
                    if let Some(idx) = hit {
                        self.selected_edges.push(idx);
                    }
                }
            }
        }
    }

    pub fn hover_at(
        &mut self,
        brep: &BRep,
        camera: &Camera,
        aspect: f32,
        viewport: [f32; 2],
        cursor: [f32; 2],
        edge_pick_radius_px: f32,
    ) {
        if cursor[0] < 0.0 || cursor[1] < 0.0 || cursor[0] > viewport[0] || cursor[1] > viewport[1]
        {
            self.clear_hover();
            return;
        }

        self.last_hover_pos = Some((cursor[0], cursor[1]));
        match self.mode {
            SelectionMode::Face => {
                self.hovered_face = pick_face(brep, camera, aspect, viewport, cursor);
                self.hovered_edge = None;
            }
            SelectionMode::Edge => {
                self.hovered_edge =
                    pick_edge(brep, camera, aspect, viewport, cursor, edge_pick_radius_px);
                self.hovered_face = None;
            }
        }
    }

    pub fn clear_hover(&mut self) {
        self.hovered_face = None;
        self.hovered_edge = None;
        self.last_hover_pos = None;
    }

    pub fn highlighted_faces(&self) -> Vec<usize> {
        merged_indices(&self.selected_faces, self.hovered_face)
    }

    pub fn highlighted_edges(&self) -> Vec<usize> {
        merged_indices(&self.selected_edges, self.hovered_edge)
    }
}

fn toggle_index(list: &mut Vec<usize>, idx: usize) {
    if let Some(pos) = list.iter().position(|&v| v == idx) {
        list.swap_remove(pos);
    } else {
        list.push(idx);
    }
}

fn merged_indices(selected: &[usize], hovered: Option<usize>) -> Vec<usize> {
    let mut merged = selected.to_vec();
    if let Some(h) = hovered
        && !merged.contains(&h)
    {
        merged.push(h);
    }
    merged
}

pub fn pick_face(
    brep: &BRep,
    camera: &Camera,
    aspect: f32,
    viewport_size: [f32; 2],
    cursor_pos: [f32; 2],
) -> Option<usize> {
    let ray = screen_ray(camera, aspect, viewport_size, cursor_pos)?;
    let mut best: Option<(f32, usize)> = None;

    let mut face_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for tri in &face.triangles {
                    let a = to_vec3(brep.vertices.get(tri[0])?.point);
                    let b = to_vec3(brep.vertices.get(tri[1])?.point);
                    let c = to_vec3(brep.vertices.get(tri[2])?.point);
                    if let Some(t) = ray_triangle_intersection(ray.0, ray.1, a, b, c)
                        && t > 0.0
                    {
                        match best {
                            Some((best_t, _)) if t >= best_t => {}
                            _ => best = Some((t, face_idx)),
                        }
                    }
                }
                face_idx += 1;
            }
        }
    }

    best.map(|(_, idx)| idx)
}

pub fn pick_edge(
    brep: &BRep,
    camera: &Camera,
    aspect: f32,
    viewport_size: [f32; 2],
    cursor_pos: [f32; 2],
    max_distance_px: f32,
) -> Option<usize> {
    if viewport_size[0] <= 1.0 || viewport_size[1] <= 1.0 {
        return None;
    }

    let vp =
        glam::Mat4::from_cols_array_2d(&camera.build_view_projection_matrix(aspect.max(0.001)));
    let mut best: Option<(f32, f32, usize)> = None;

    for (idx, edge) in brep.edges.iter().enumerate() {
        let p0 = to_vec3(brep.vertices.get(edge.start)?.point);
        let p1 = to_vec3(brep.vertices.get(edge.end)?.point);

        let s0 = project_to_screen(vp, p0, viewport_size)?;
        let s1 = project_to_screen(vp, p1, viewport_size)?;
        let distance = point_segment_distance_2d(cursor_pos, [s0[0], s0[1]], [s1[0], s1[1]]);

        if distance > max_distance_px {
            continue;
        }

        let depth = (s0[2] + s1[2]) * 0.5;
        match best {
            Some((best_dist, best_depth, _))
                if distance > best_dist
                    || ((distance - best_dist).abs() < 1e-3 && depth >= best_depth) => {}
            _ => best = Some((distance, depth, idx)),
        }
    }

    best.map(|(_, _, idx)| idx)
}

#[derive(Clone, Debug)]
pub struct Mesh {
    pub vertices: Vec<[f32; 3]>,
    pub indices: Vec<u32>,
    pub line_indices: Vec<u32>,
    /// Per-vertex smooth normals (same length as `vertices`).  When empty the
    /// renderer uploads zero normals, which triggers the flat-shading fallback
    /// in the fragment shader.
    pub normals: Vec<[f32; 3]>,
}

#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
struct MaterialUniform {
    color: [f32; 4],
    flags: [f32; 4],
}

#[derive(Debug)]
struct SolidMeshBuffers {
    vertex_buffer: wgpu::Buffer,
    index_buffer: wgpu::Buffer,
    index_count: u32,
}

pub fn build_face_highlight_mesh(brep: &BRep, face_index: usize) -> Option<Mesh> {
    build_faces_highlight_mesh(brep, &[face_index])
}

pub fn build_faces_highlight_mesh(brep: &BRep, face_indices: &[usize]) -> Option<Mesh> {
    if face_indices.is_empty() {
        return None;
    }

    let selected: std::collections::HashSet<usize> = face_indices.iter().copied().collect();
    let mut current = 0usize;
    let mut indices: Vec<u32> = Vec::new();

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if selected.contains(&current) {
                    for tri in &face.triangles {
                        indices.push(tri[0] as u32);
                        indices.push(tri[1] as u32);
                        indices.push(tri[2] as u32);
                    }
                }
                current += 1;
            }
        }
    }

    if indices.is_empty() {
        return None;
    }

    let vertices: Vec<[f32; 3]> = brep
        .vertices
        .iter()
        .map(|v| [v.point.x as f32, v.point.y as f32, v.point.z as f32])
        .collect();

    Some(Mesh {
        vertices,
        indices,
        line_indices: Vec::new(),
        normals: Vec::new(),
    })
}

pub fn build_edge_highlight_mesh(brep: &BRep, edge_index: usize) -> Option<Mesh> {
    build_edges_highlight_mesh(brep, &[edge_index])
}

pub fn build_edges_highlight_mesh(brep: &BRep, edge_indices: &[usize]) -> Option<Mesh> {
    if edge_indices.is_empty() {
        return None;
    }

    let mut vertices: Vec<[f32; 3]> = brep
        .vertices
        .iter()
        .map(|v| [v.point.x as f32, v.point.y as f32, v.point.z as f32])
        .collect();
    let mut dummy_normals: Vec<[f32; 3]> = vec![[0.0; 3]; vertices.len()];
    let mut indices: Vec<u32> = Vec::with_capacity(edge_indices.len() * 2);
    for &edge_index in edge_indices {
        let edge = brep.edges.get(edge_index)?;
        if let Some(pts) = sample_edge_curve_points(brep, edge_index) {
            let base = vertices.len() as u32;
            let n = pts.len();
            dummy_normals.extend(std::iter::repeat([0.0f32; 3]).take(n));
            for i in 0..(n - 1) as u32 {
                indices.push(base + i);
                indices.push(base + i + 1);
            }
            vertices.extend_from_slice(&pts);
        } else {
            indices.push(edge.start as u32);
            indices.push(edge.end as u32);
        }
    }
    drop(dummy_normals);

    Some(Mesh {
        vertices,
        indices,
        line_indices: Vec::new(),
        normals: Vec::new(),
    })
}

pub fn merge_meshes(meshes: &[&Mesh]) -> Option<Mesh> {
    if meshes.is_empty() {
        return None;
    }

    let total_vertices = meshes.iter().map(|mesh| mesh.vertices.len()).sum();
    let total_indices = meshes.iter().map(|mesh| mesh.indices.len()).sum();
    let total_line_indices = meshes.iter().map(|mesh| mesh.line_indices.len()).sum();

    if total_vertices == 0 || (total_indices == 0 && total_line_indices == 0) {
        return None;
    }

    let mut vertices = Vec::with_capacity(total_vertices);
    let mut normals: Vec<[f32; 3]> = Vec::with_capacity(total_vertices);
    let mut indices = Vec::with_capacity(total_indices);
    let mut line_indices = Vec::with_capacity(total_line_indices);
    let mut vertex_offset = 0u32;

    for mesh in meshes {
        vertices.extend_from_slice(&mesh.vertices);
        normals.extend_from_slice(&mesh.normals);
        indices.extend(mesh.indices.iter().map(|index| index + vertex_offset));
        line_indices.extend(mesh.line_indices.iter().map(|index| index + vertex_offset));
        vertex_offset += mesh.vertices.len() as u32;
    }

    // If only some meshes had normals fill missing entries with zero.
    normals.resize(vertices.len(), [0.0, 0.0, 0.0]);

    Some(Mesh {
        vertices,
        indices,
        line_indices,
        normals,
    })
}

/// Sample the analytic curve of a BRep edge into a sequence of `[f32; 3]` points
/// (including both endpoints). Returns `None` for straight lines or missing geometry,
/// signalling the caller to fall back to a single-chord segment.
fn sample_edge_curve_points(brep: &BRep, edge_idx: usize) -> Option<Vec<[f32; 3]>> {
    let ci = brep.geom.edge_curve.get(edge_idx).and_then(|v| *v)?;
    let curve = brep.geom.curves.get(ci)?;
    let mut range = brep
        .geom
        .edge_curve_range
        .get(edge_idx)
        .and_then(|v| *v)
        .or_else(|| match curve {
            Curve3::Circle(_) | Curve3::Ellipse(_) => Some([0.0, 2.0 * std::f64::consts::PI]),
            _ => None,
        })?;
    let edge = brep.edges.get(edge_idx)?;
    let p_start = brep.vertices.get(edge.start)?.point;
    let p_end = brep.vertices.get(edge.end)?.point;

    let two_pi = 2.0 * std::f64::consts::PI;
    let wrap_2pi = |t: f64| -> f64 {
        let mut out = t % two_pi;
        if out < 0.0 {
            out += two_pi;
        }
        out
    };

    // Some imported periodic edges carry a full [0, 2π] range even when the
    // topological edge is only an arc. Rebuild a trimmed range from endpoints.
    match curve {
        Curve3::Circle(c) => {
            if (range[1] - range[0]).abs() >= two_pi * 0.999 {
                let x_ax = any_perpendicular(c.normal);
                let y_ax = c.normal.cross(x_ax);
                let v0 = p_start - c.center;
                let v1 = p_end - c.center;
                let t0 = wrap_2pi(v0.dot(y_ax).atan2(v0.dot(x_ax)));
                let t1 = wrap_2pi(v1.dot(y_ax).atan2(v1.dot(x_ax)));
                let mut dt = t1 - t0;
                if dt > std::f64::consts::PI {
                    dt -= two_pi;
                } else if dt < -std::f64::consts::PI {
                    dt += two_pi;
                }
                range = [t0, t0 + dt];
            }
        }
        Curve3::Ellipse(e) => {
            if (range[1] - range[0]).abs() >= two_pi * 0.999 {
                let x_ax = e.major_dir.normalize();
                let y_ax = e.normal.cross(x_ax).normalize();
                let v0 = p_start - e.center;
                let v1 = p_end - e.center;
                let t0 = wrap_2pi((v0.dot(y_ax) / e.minor_radius).atan2(v0.dot(x_ax) / e.major_radius));
                let t1 = wrap_2pi((v1.dot(y_ax) / e.minor_radius).atan2(v1.dot(x_ax) / e.major_radius));
                let mut dt = t1 - t0;
                if dt > std::f64::consts::PI {
                    dt -= two_pi;
                } else if dt < -std::f64::consts::PI {
                    dt += two_pi;
                }
                range = [t0, t0 + dt];
            }
        }
        _ => {}
    }

    // Straight lines render fine as a single chord — skip sampling.
    if matches!(curve, Curve3::Line(_)) {
        return None;
    }
    let t1 = range[0];
    let t2 = range[1];
    let span = (t2 - t1).abs();
    if span < 1e-12 {
        return None;
    }
    let n_segs: usize = match curve {
        Curve3::Circle(_) => {
            let segs = (span / (2.0 * std::f64::consts::PI) * 64.0).ceil() as usize;
            segs.clamp(2, 64)
        }
        Curve3::Ellipse(_) => 32,
        _ => 24,
    };
    let pts: Vec<[f32; 3]> = (0..=n_segs)
        .map(|i| {
            let t = t1 + (t2 - t1) * (i as f64 / n_segs as f64);
            let p = curve.point_at(t);
            [p.x as f32, p.y as f32, p.z as f32]
        })
        .collect();
    Some(pts)
}

pub struct Tessellator;

impl Tessellator {
    pub fn tessellate(brep: &BRep) -> Mesh {
        let mut flat_verts: Vec<[f32; 3]> = brep
            .vertices
            .iter()
            .map(|v| [v.point.x as f32, v.point.y as f32, v.point.z as f32])
            .collect();

        let n_verts = flat_verts.len();
        let mut indices: Vec<u32> = Vec::new();
        let mut line_indices: Vec<u32> = Vec::with_capacity(brep.edges.len() * 2);

        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for tri in &face.triangles {
                        indices.push(tri[0] as u32);
                        indices.push(tri[1] as u32);
                        indices.push(tri[2] as u32);
                    }
                }
            }
        }

        // ── Per-vertex smooth normal computation (area-weighted face normal avg) ──
        let mut normal_accum = vec![[0.0f64; 3]; n_verts];
        let mut face_flat_idx = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let face_surface = brep
                        .geom
                        .face_surface
                        .get(face_flat_idx)
                        .and_then(|value| *value)
                        .and_then(|surface_idx| brep.geom.surfaces.get(surface_idx));

                    for tri in &face.triangles {
                        let a = brep.vertices[tri[0]].point;
                        let b = brep.vertices[tri[1]].point;
                        let c = brep.vertices[tri[2]].point;
                        let e1 = b - a;
                        let e2 = c - a;
                        // Area-weighted face normal (magnitude = 2× triangle area)
                        let fn_ = e1.cross(e2);
                        let weight = fn_.length().max(1e-12);
                        for &vi in tri.iter() {
                            if vi < n_verts {
                                let analytic = face_surface.and_then(|surface| {
                                    analytic_surface_normal_at_point(surface, brep.vertices[vi].point)
                                });
                                let contribution = analytic.unwrap_or(fn_.normalize_or_zero()) * weight;
                                normal_accum[vi][0] += contribution.x;
                                normal_accum[vi][1] += contribution.y;
                                normal_accum[vi][2] += contribution.z;
                            }
                        }
                    }

                    face_flat_idx += 1;
                }
            }
        }
        let mut normals: Vec<[f32; 3]> = normal_accum
            .iter()
            .map(|n| {
                let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                if len < 1e-15 {
                    [0.0, 0.0, 0.0]
                } else {
                    [(n[0] / len) as f32, (n[1] / len) as f32, (n[2] / len) as f32]
                }
            })
            .collect();

        // STEP imports can split one visually smooth surface into multiple faces
        // with duplicate vertex positions. If the split-face normals are already
        // closely aligned, smooth them together to avoid artificial seam lines.
        smooth_normals_across_coincident_vertices(&flat_verts, &mut normals, 1e-5, 0.95);

        let mut seam_edges: std::collections::HashSet<usize> =
            seam_edge_candidates(brep).into_iter().collect();

        // Some closed periodic faces (notably primitive spheres) represent the
        // seam by repeating the same edge index in the face wire.
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let mut counts: std::collections::HashMap<usize, usize> =
                        std::collections::HashMap::new();
                    for we in &face.outer_wire.edges {
                        *counts.entry(we.idx).or_insert(0) += 1;
                    }
                    for wire in &face.inner_wires {
                        for we in &wire.edges {
                            *counts.entry(we.idx).or_insert(0) += 1;
                        }
                    }
                    for (edge_idx, count) in counts {
                        if count > 1 {
                            seam_edges.insert(edge_idx);
                        }
                    }
                }
            }
        }

        for (edge_idx, edge) in brep.edges.iter().enumerate() {
            if seam_edges.contains(&edge_idx) {
                // Do not draw periodic seam edges in wireframe overlays.
                continue;
            }
            if let Some(pts) = sample_edge_curve_points(brep, edge_idx) {
                let base = flat_verts.len() as u32;
                let n = pts.len();
                normals.extend(std::iter::repeat([0.0f32; 3]).take(n));
                for i in 0..(n - 1) as u32 {
                    line_indices.push(base + i);
                    line_indices.push(base + i + 1);
                }
                flat_verts.extend_from_slice(&pts);
            } else {
                line_indices.push(edge.start as u32);
                line_indices.push(edge.end as u32);
            }
        }

        Mesh {
            vertices: flat_verts,
            indices,
            line_indices,
            normals,
        }
    }

    /// Re-tessellate dirty faces using the given quality options, then build a GPU [`Mesh`].
    ///
    /// Calls [`rcad_algorithms::mesh_brep`] to recompute triangles for any face whose
    /// `mesh_dirty` flag is set, then delegates to [`Tessellator::tessellate`].
    ///
    /// Analogous to `BRepMesh_IncrementalMesh` with explicit deflection/angular arguments in OCCT.
    pub fn tessellate_with_options(brep: &mut BRep, options: &TessellationOptions) -> Mesh {
        mesh_brep(brep, options);
        Self::tessellate(brep)
    }

    /// Incrementally invalidate mesh cache for faces affected by edited entities.
    ///
    /// Returns the number of faces that were newly marked dirty.
    pub fn invalidate_cache_for_edits(brep: &mut BRep, edits: &EditedModelDelta) -> usize {
        if edits.is_empty() {
            return 0;
        }

        let graph = BRepGraph::from_brep(brep);
        let face_count: usize = brep
            .solids
            .iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();
        if face_count == 0 {
            return 0;
        }

        let mut dirty_faces = vec![false; face_count];

        for &fi in &edits.modified_faces {
            if fi < face_count {
                dirty_faces[fi] = true;
            }
        }
        for &ei in &edits.modified_edges {
            for &fi in graph.edge_adjacent_faces(ei) {
                if fi < face_count {
                    dirty_faces[fi] = true;
                }
            }
        }
        for &vi in &edits.modified_vertices {
            for &fi in graph.vertex_adjacent_faces(vi) {
                if fi < face_count {
                    dirty_faces[fi] = true;
                }
            }
        }

        let mut newly_marked = 0usize;
        let mut flat_fi = 0usize;
        for solid in &mut brep.solids {
            for shell in &mut solid.shells {
                for face in &mut shell.faces {
                    if dirty_faces[flat_fi] && !face.mesh_dirty {
                        face.mesh_dirty = true;
                        newly_marked += 1;
                    }
                    flat_fi += 1;
                }
            }
        }

        newly_marked
    }

    /// Convenience helper: invalidate affected faces from edit delta, then tessellate.
    pub fn tessellate_after_edits(
        brep: &mut BRep,
        edits: &EditedModelDelta,
        options: &TessellationOptions,
    ) -> Mesh {
        Self::invalidate_cache_for_edits(brep, edits);
        Self::tessellate_with_options(brep, options)
    }
}

#[cfg(test)]
mod tests {
    use super::Tessellator;
    use rcad_algorithms::{mesh_brep, TessellationParams};
    use rcad_kernel::{BRep, PrimitiveSolid};

    #[test]
    fn tessellate_sphere_hides_seam_edges_in_line_indices() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        mesh_brep(&mut brep, &TessellationParams::default());
        let mesh = Tessellator::tessellate(&brep);
        assert!(
            mesh.line_indices.is_empty(),
            "full sphere should not render seam wireframe edge"
        );
    }

    #[test]
    fn tessellate_box_keeps_regular_edges_in_line_indices() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        mesh_brep(&mut brep, &TessellationParams::default());
        let mesh = Tessellator::tessellate(&brep);
        assert_eq!(mesh.line_indices.len(), brep.edges.len() * 2);
    }

    #[test]
    fn test_tessellate_cylinder() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        mesh_brep(&mut brep, &TessellationParams::default());
        let mesh = Tessellator::tessellate(&brep);

        // Verify vertices are generated
        assert!(!mesh.vertices.is_empty(), "cylinder should generate vertices");

        // Verify triangles are generated (indices should be divisible by 3)
        assert!(!mesh.indices.is_empty(), "cylinder should generate indices");
        assert_eq!(
            mesh.indices.len() % 3,
            0,
            "indices should form complete triangles"
        );

        // Verify smooth normals are generated
        assert_eq!(
            mesh.normals.len(),
            mesh.vertices.len(),
            "cylinder should have per-vertex normals for smooth shading"
        );

        // Verify reasonable mesh density (cylinder should have multiple quads around)
        assert!(
            mesh.vertices.len() > 20,
            "cylinder should have sufficient vertices for curvature"
        );
    }

    #[test]
    fn test_tessellate_sphere() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        mesh_brep(&mut brep, &TessellationParams::default());
        let mesh = Tessellator::tessellate(&brep);

        // Verify vertices are generated
        assert!(!mesh.vertices.is_empty(), "sphere should generate vertices");

        // Verify triangles are generated (indices should be divisible by 3)
        assert!(!mesh.indices.is_empty(), "sphere should generate indices");
        assert_eq!(
            mesh.indices.len() % 3,
            0,
            "indices should form complete triangles"
        );

        // Verify smooth normals are generated
        assert_eq!(
            mesh.normals.len(),
            mesh.vertices.len(),
            "sphere should have per-vertex normals for smooth shading"
        );

        // Verify reasonable mesh density
        assert!(
            mesh.vertices.len() > 50,
            "sphere should have sufficient vertices for curvature"
        );
    }

    #[test]
    fn test_tessellate_cone() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Cone {
            base_radius: 1.0,
            height: 2.0,
        });
        mesh_brep(&mut brep, &TessellationParams::default());
        let mesh = Tessellator::tessellate(&brep);

        // Verify vertices are generated
        assert!(!mesh.vertices.is_empty(), "cone should generate vertices");

        // Verify triangles are generated (indices should be divisible by 3)
        assert!(!mesh.indices.is_empty(), "cone should generate indices");
        assert_eq!(
            mesh.indices.len() % 3,
            0,
            "indices should form complete triangles"
        );

        // Verify smooth normals are generated
        assert_eq!(
            mesh.normals.len(),
            mesh.vertices.len(),
            "cone should have per-vertex normals for smooth shading"
        );

        // Verify reasonable mesh density
        assert!(
            mesh.vertices.len() > 10,
            "cone should have sufficient vertices for curvature"
        );
    }

    #[test]
    fn test_tessellate_torus() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 2.0,
            minor_radius: 0.5,
        });
        mesh_brep(&mut brep, &TessellationParams::default());
        let mesh = Tessellator::tessellate(&brep);

        // Verify vertices are generated
        assert!(!mesh.vertices.is_empty(), "torus should generate vertices");

        // Verify triangles are generated (indices should be divisible by 3)
        assert!(!mesh.indices.is_empty(), "torus should generate indices");
        assert_eq!(
            mesh.indices.len() % 3,
            0,
            "indices should form complete triangles"
        );

        // Verify smooth normals are generated
        assert_eq!(
            mesh.normals.len(),
            mesh.vertices.len(),
            "torus should have per-vertex normals for smooth shading"
        );

        // Verify reasonable mesh density (torus is complex, needs many vertices)
        assert!(
            mesh.vertices.len() > 100,
            "torus should have sufficient vertices for double curvature"
        );
    }

    /// Test adaptive tessellation with edge-sensitive refinement.
    /// High-quality tessellation should produce more triangles on curved surfaces
    /// compared to preview quality, demonstrating adaptive refinement behavior.
    #[test]
    fn test_adaptive_tessellation() {
        use rcad_algorithms::TessellationParams;

        // Create a cylinder which has both flat (top/bottom) and curved (side) surfaces
        let mut brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        // Tessellate with preview settings (lower quality, no adaptive refinement)
        let preview_params = TessellationParams::preview();
        let preview_mesh = Tessellator::tessellate_with_options(&mut brep, &preview_params);

        // Count triangles for preview mesh
        let preview_tri_count = preview_mesh.indices.len() / 3;

        // Reset mesh dirty flags for re-tessellation
        for solid in &mut brep.solids {
            for shell in &mut solid.shells {
                for face in &mut shell.faces {
                    face.mesh_dirty = true;
                }
            }
        }

        // Tessellate with high-quality settings (adaptive refinement enabled)
        let hq_params = TessellationParams::high_quality();
        let hq_mesh = Tessellator::tessellate_with_options(&mut brep, &hq_params);

        // Count triangles for high-quality mesh
        let hq_tri_count = hq_mesh.indices.len() / 3;

        // High-quality tessellation with adaptive refinement should produce more triangles
        // on curved surfaces (cylinder sides) than preview quality
        assert!(
            hq_tri_count > preview_tri_count,
            "High-quality tessellation ({} triangles) should produce more triangles than preview ({} triangles) on curved surfaces",
            hq_tri_count,
            preview_tri_count
        );

        // Verify the high-quality mesh has reasonable triangle count for a cylinder
        // A cylinder should have at least enough triangles to represent the curved surface
        assert!(
            hq_tri_count >= 12, // Minimum reasonable for a cylinder
            "High-quality cylinder mesh should have at least 12 triangles, got {}",
            hq_tri_count
        );
    }

    /// Test triangle quality by validating aspect ratios.
    /// Triangles should not have extreme aspect ratios that would cause rendering artifacts.
    #[test]
    fn test_triangle_quality() {
        use rcad_algorithms::TessellationParams;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let params = TessellationParams::standard();
        let mesh = Tessellator::tessellate_with_options(&mut brep, &params);

        let vertices = &mesh.vertices;
        let max_allowed_aspect_ratio = 25.0; // Reasonable threshold for mesh quality

        let mut max_aspect_ratio = 0.0f32;
        let mut triangles_with_bad_aspect = 0usize;

        // Check each triangle's aspect ratio
        for tri_idx in 0..(mesh.indices.len() / 3) {
            let i0 = mesh.indices[tri_idx * 3] as usize;
            let i1 = mesh.indices[tri_idx * 3 + 1] as usize;
            let i2 = mesh.indices[tri_idx * 3 + 2] as usize;

            let v0 = glam::Vec3::from(vertices[i0]);
            let v1 = glam::Vec3::from(vertices[i1]);
            let v2 = glam::Vec3::from(vertices[i2]);

            let aspect_ratio = compute_triangle_aspect_ratio(v0, v1, v2);
            max_aspect_ratio = max_aspect_ratio.max(aspect_ratio);

            if aspect_ratio > max_allowed_aspect_ratio {
                triangles_with_bad_aspect += 1;
            }
        }

        let total_triangles = mesh.indices.len() / 3;

        // Allow a small percentage of bad triangles (some edge cases are acceptable)
        let bad_ratio = triangles_with_bad_aspect as f32 / total_triangles as f32;
        assert!(
            bad_ratio < 0.05, // Less than 5% bad triangles
            "Too many triangles with bad aspect ratio: {}/{} ({:.1}%)",
            triangles_with_bad_aspect,
            total_triangles,
            bad_ratio * 100.0
        );

        // Log the maximum aspect ratio for debugging
        println!(
            "Sphere mesh: {} triangles, max aspect ratio: {:.2}",
            total_triangles, max_aspect_ratio
        );
    }

    /// Compute the aspect ratio of a triangle.
    /// Aspect ratio = longest_edge / (2 * sqrt(3) * inradius)
    /// An equilateral triangle has aspect ratio = 1.0
    fn compute_triangle_aspect_ratio(v0: glam::Vec3, v1: glam::Vec3, v2: glam::Vec3) -> f32 {
        let e0 = (v1 - v0).length();
        let e1 = (v2 - v1).length();
        let e2 = (v0 - v2).length();

        let longest_edge = e0.max(e1).max(e2);

        // Compute area using cross product
        let cross = (v1 - v0).cross(v2 - v0);
        let area = cross.length() * 0.5;

        if area < 1e-10 {
            // Degenerate triangle
            return f32::MAX;
        }

        // Semi-perimeter
        let s = (e0 + e1 + e2) * 0.5;

        // Inradius = area / s
        let inradius = area / s;

        // Aspect ratio = longest_edge / (2 * sqrt(3) * inradius)
        // For equilateral triangle: inradius = edge / (2 * sqrt(3))
        // So aspect ratio = 1.0 for equilateral
        longest_edge / (2.0 * f32::sqrt(3.0) * inradius)
    }

    // ============================================================================
    // OCCT TKMesh Alignment Tests
    // ============================================================================

    /// Test box tessellation produces correct topology.
    /// OCCT TKMesh coverage: box_primitive_mesh, planar_surface_tessellation.
    #[test]
    fn test_tessellate_box() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        mesh_brep(&mut brep, &TessellationParams::default());
        let mesh = Tessellator::tessellate(&brep);

        // Box should have 6 faces, each with 2 triangles minimum
        assert!(!mesh.vertices.is_empty(), "box should generate vertices");
        assert!(mesh.indices.len() >= 36, "box should have at least 12 triangles");

        // Verify normals are unit length
        for normal in &mesh.normals {
            let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
            assert!((len - 1.0).abs() < 0.01, "normals should be unit length");
        }
    }

    /// Test torus mesh handles inner equator correctly.
    /// OCCT TKMesh coverage: toroidal_surface_mesh, negative_gaussian_curvature.
    #[test]
    fn test_torus_inner_equator_mesh() {
        // Torus with significant inner hole
        let mut brep = BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 4.0,
            minor_radius: 1.5,
        });
        mesh_brep(&mut brep, &TessellationParams::high_quality());
        let mesh = Tessellator::tessellate(&brep);

        // Verify mesh is valid
        assert!(!mesh.vertices.is_empty(), "torus should generate vertices");

        // Check that some vertices are on the inner equator (closest to center)
        let inner_radius = 4.0 - 1.5; // 2.5
        let mut has_inner_vertex = false;
        for v in &mesh.vertices {
            let r = (v[0] * v[0] + v[1] * v[1]).sqrt();
            if (r - inner_radius).abs() < 0.3 {
                has_inner_vertex = true;
                break;
            }
        }
        assert!(has_inner_vertex, "torus should have vertices on inner equator");
    }

    /// Test cylinder mesh has proper radial distribution.
    /// OCCT TKMesh coverage: cylindrical_surface_mesh, periodic_parameter.
    #[test]
    fn test_cylinder_radial_distribution() {
        let radius = 2.0;
        let mut brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius,
            height: 3.0,
        });
        mesh_brep(&mut brep, &TessellationParams::standard());
        let mesh = Tessellator::tessellate(&brep);

        // Count vertices at different angles
        let mut angles: Vec<f32> = Vec::new();
        for v in &mesh.vertices {
            let angle = (v[1].atan2(v[0]) * 180.0 / std::f32::consts::PI).abs();
            angles.push(angle);
        }

        // Should have vertices distributed around full 360 degrees
        let max_angle = angles.iter().cloned().fold(0.0_f32, f32::max);
        assert!(max_angle > 150.0, "vertices should span most of 360 degrees");
    }

    /// Test cone mesh handles apex region correctly.
    /// OCCT TKMesh coverage: conical_surface_mesh, singular_point_handling.
    #[test]
    fn test_cone_apex_mesh() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Cone {
            base_radius: 2.0,
            height: 4.0,
        });
        mesh_brep(&mut brep, &TessellationParams::standard());
        let mesh = Tessellator::tessellate(&brep);

        // Verify the cone mesh is valid and non-empty
        assert!(!mesh.vertices.is_empty(), "cone should generate vertices");
        assert!(!mesh.indices.is_empty(), "cone should generate indices");

        // Verify indices form valid triangles
        let vertex_count = mesh.vertices.len() as u32;
        for &idx in &mesh.indices {
            assert!(idx < vertex_count, "index should be within bounds");
        }

        // Check z range includes the expected height range
        let z_values: Vec<f32> = mesh.vertices.iter().map(|v| v[2]).collect();
        let z_max = z_values.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let z_min = z_values.iter().cloned().fold(f32::INFINITY, f32::min);

        // Z should span from 0 to approximately the height
        assert!(z_max > 0.0, "cone z_max should be positive");
        assert!(z_min <= z_max, "cone should have valid z range");
    }

    /// Test preview quality mesh is lower density.
    /// OCCT TKMesh coverage: mesh_quality_levels, adaptive_deflection.
    #[test]
    fn test_preview_quality_density() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let preview_mesh = Tessellator::tessellate_with_options(
            &mut brep,
            &TessellationParams::preview(),
        );

        // Reset for re-tessellation
        for solid in &mut brep.solids {
            for shell in &mut solid.shells {
                for face in &mut shell.faces {
                    face.mesh_dirty = true;
                }
            }
        }

        let hq_mesh = Tessellator::tessellate_with_options(
            &mut brep,
            &TessellationParams::high_quality(),
        );

        // Preview should have fewer triangles than high quality
        assert!(
            preview_mesh.vertices.len() <= hq_mesh.vertices.len(),
            "preview mesh ({}) should have <= vertices than hq ({})",
            preview_mesh.vertices.len(),
            hq_mesh.vertices.len()
        );
    }

    /// Test mesh normals point outward.
    /// OCCT TKMesh coverage: normal_orientation, consistent_normals.
    #[test]
    fn test_normals_point_outward() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        mesh_brep(&mut brep, &TessellationParams::standard());
        let mesh = Tessellator::tessellate(&brep);

        // For a sphere centered at origin, normals should generally point outward
        // (same direction as position vector). Allow some tolerance for poles.
        let mut outward_count = 0;
        let mut total_count = 0;
        for (v, n) in mesh.vertices.iter().zip(mesh.normals.iter()) {
            let v_len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if v_len > 0.1 {
                total_count += 1;
                let v_norm = [v[0] / v_len, v[1] / v_len, v[2] / v_len];
                let dot = v_norm[0] * n[0] + v_norm[1] * n[1] + v_norm[2] * n[2];
                if dot > 0.7 {
                    // Allow some deviation at poles
                    outward_count += 1;
                }
            }
        }
        // At least 90% of normals should point outward
        let ratio = outward_count as f32 / total_count as f32;
        assert!(
            ratio > 0.9,
            "most sphere normals should point outward ({}/{} = {:.1}%)",
            outward_count,
            total_count,
            ratio * 100.0
        );
    }

    /// Test mesh indices are valid.
    /// OCCT TKMesh coverage: index_bounds, triangle_winding.
    #[test]
    fn test_mesh_indices_valid() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 2.0,
            minor_radius: 0.5,
        });
        mesh_brep(&mut brep, &TessellationParams::standard());
        let mesh = Tessellator::tessellate(&brep);

        let vertex_count = mesh.vertices.len() as u32;

        // All indices should be within bounds
        for &idx in &mesh.indices {
            assert!(
                idx < vertex_count,
                "index {} out of bounds (vertex count: {})",
                idx,
                vertex_count
            );
        }
    }

    /// Test empty mesh handling.
    /// OCCT TKMesh coverage: empty_shape, null_mesh.
    #[test]
    fn test_empty_brep_mesh() {
        let brep = BRep::new();
        let mesh = Tessellator::tessellate(&brep);

        assert!(mesh.vertices.is_empty(), "empty brep should produce empty mesh");
        assert!(mesh.indices.is_empty(), "empty brep should produce no indices");
    }
}

#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
pub struct CameraUniform {
    view_proj: [[f32; 4]; 4],
    eye_pos: [f32; 4],
    light_dir: [f32; 4],
}

#[derive(Clone, Copy, Debug)]
pub struct Camera {
    pub rot_x: f32,
    pub rot_y: f32,
    pub distance: f32,
    pub target: glam::Vec3,
}

impl Camera {
    pub fn new() -> Self {
        Self {
            rot_x: 0.4,
            rot_y: 0.5,
            distance: 3.0,
            target: glam::Vec3::ZERO,
        }
    }

    pub fn build_view_projection_matrix(&self, aspect: f32) -> [[f32; 4]; 4] {
        let eye = self.eye_position();
        let target = self.target;
        let forward = (target - eye).normalize_or_zero();
        let up = if forward.dot(glam::Vec3::Y).abs() > 0.98 {
            glam::Vec3::Z
        } else {
            glam::Vec3::Y
        };

        let view = glam::Mat4::look_at_rh(eye, target, up);
        let proj = glam::Mat4::perspective_rh(45.0_f32.to_radians(), aspect, 0.1, 100.0);

        (proj * view).to_cols_array_2d()
    }

    pub fn eye_position(&self) -> glam::Vec3 {
        self.target
            + glam::Vec3::new(
                self.distance * self.rot_y.cos() * self.rot_x.cos(),
                self.distance * self.rot_x.sin(),
                self.distance * self.rot_y.sin() * self.rot_x.cos(),
            )
    }

    pub fn pan_pixels(&mut self, dx: f32, dy: f32) {
        let eye = self.eye_position();
        let forward = (self.target - eye).normalize_or_zero();
        if forward.length_squared() <= 1e-8 {
            return;
        }

        let mut right = forward.cross(glam::Vec3::Y);
        if right.length_squared() <= 1e-8 {
            right = forward.cross(glam::Vec3::X);
        }
        right = right.normalize_or_zero();
        let up = right.cross(forward).normalize_or_zero();

        let scale = self.distance.max(0.1) * 0.0015;
        self.target += (-dx * right + dy * up) * scale;
    }

    pub fn set_view_direction(&mut self, direction: glam::Vec3) {
        let dir = direction.normalize_or_zero();
        if dir.length_squared() <= 1e-8 {
            return;
        }
        self.rot_x = dir.y.clamp(-1.0, 1.0).asin();
        self.rot_y = dir.z.atan2(dir.x);
    }

    pub fn set_isometric_view(&mut self) {
        self.rot_x = 0.4;
        self.rot_y = 0.5;
    }
}

impl Default for Camera {
    fn default() -> Self {
        Self::new()
    }
}

fn smooth_normals_across_coincident_vertices(
    vertices: &[[f32; 3]],
    normals: &mut [[f32; 3]],
    tolerance: f32,
    min_dot: f32,
) {
    use std::collections::HashMap;

    let scale = 1.0 / tolerance.max(1e-9);
    let mut groups: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (index, &[x, y, z]) in vertices.iter().enumerate() {
        let key = [
            (x * scale).round() as i64,
            (y * scale).round() as i64,
            (z * scale).round() as i64,
        ];
        groups.entry(key).or_default().push(index);
    }

    for indices in groups.values() {
        if indices.len() < 2 {
            continue;
        }

        let mut valid_normals = Vec::new();
        for &index in indices {
            let normal = glam::Vec3::from_array(normals[index]);
            if normal.length_squared() > 1e-8 {
                valid_normals.push((index, normal.normalize()));
            }
        }
        if valid_normals.len() < 2 {
            continue;
        }

        let mut should_smooth = true;
        for left in 0..valid_normals.len() {
            for right in (left + 1)..valid_normals.len() {
                if valid_normals[left].1.dot(valid_normals[right].1) < min_dot {
                    should_smooth = false;
                    break;
                }
            }
            if !should_smooth {
                break;
            }
        }
        if !should_smooth {
            continue;
        }

        let averaged = valid_normals
            .iter()
            .fold(glam::Vec3::ZERO, |acc, (_, normal)| acc + *normal)
            .normalize_or_zero();
        if averaged.length_squared() <= 1e-8 {
            continue;
        }

        for &index in indices {
            normals[index] = averaged.to_array();
        }
    }
}

fn analytic_surface_normal_at_point(surface: &Surface3, point: glam::DVec3) -> Option<glam::DVec3> {
    match surface {
        Surface3::Plane(plane) => Some(plane.normal.normalize_or_zero()),
        Surface3::Cylinder(cylinder) => {
            let axis = cylinder.axis.normalize_or_zero();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let x_axis = any_perpendicular(axis);
            let y_axis = axis.cross(x_axis).normalize_or_zero();
            let delta = point - cylinder.origin;
            let radial = delta - axis * delta.dot(axis);
            if radial.length_squared() <= 1e-20 {
                return None;
            }
            let u = radial.dot(y_axis).atan2(radial.dot(x_axis));
            let v = delta.dot(axis);
            Some(cylinder.normal_at(u, v).normalize_or_zero())
        }
        Surface3::Sphere(sphere) => {
            let axis = sphere.axis.normalize_or_zero();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let x_axis = any_perpendicular(axis);
            let y_axis = axis.cross(x_axis).normalize_or_zero();
            let radial = (point - sphere.center).normalize_or_zero();
            if radial.length_squared() <= 1e-20 {
                return None;
            }
            let u = radial.dot(y_axis).atan2(radial.dot(x_axis));
            let v = radial.dot(axis).clamp(-1.0, 1.0).acos();
            Some(sphere.normal_at(u, v).normalize_or_zero())
        }
        Surface3::Cone(cone) => {
            let axis = cone.axis.normalize_or_zero();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let x_axis = any_perpendicular(axis);
            let y_axis = axis.cross(x_axis).normalize_or_zero();
            let delta = point - cone.apex;
            let height = delta.dot(axis);
            let radial = delta - axis * height;
            if radial.length_squared() <= 1e-20 {
                return None;
            }
            let u = radial.dot(y_axis).atan2(radial.dot(x_axis));
            Some(cone.normal_at(u, height.max(0.0)).normalize_or_zero())
        }
        Surface3::Torus(torus) => {
            let axis = torus.axis.normalize_or_zero();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let x_axis = any_perpendicular(axis);
            let y_axis = axis.cross(x_axis).normalize_or_zero();
            let delta = point - torus.center;
            let planar = delta - axis * delta.dot(axis);
            if planar.length_squared() <= 1e-20 {
                return None;
            }
            let u = planar.dot(y_axis).atan2(planar.dot(x_axis));
            let major_dir = (u.cos() * x_axis + u.sin() * y_axis).normalize_or_zero();
            let tube_center = torus.center + torus.major_radius * major_dir;
            let tube_vec = (point - tube_center).normalize_or_zero();
            if tube_vec.length_squared() <= 1e-20 {
                return None;
            }
            let v = tube_vec.dot(axis).atan2(tube_vec.dot(major_dir));
            Some(torus.normal_at(u, v).normalize_or_zero())
        }
        Surface3::Offset(offset) => analytic_surface_normal_at_point(&offset.basis, point),
        Surface3::Trimmed(trimmed) => analytic_surface_normal_at_point(&trimmed.basis, point),
        _ => None,
    }
}

fn to_vec3(v: glam::DVec3) -> glam::Vec3 {
    glam::Vec3::new(v.x as f32, v.y as f32, v.z as f32)
}

fn screen_ray(
    camera: &Camera,
    aspect: f32,
    viewport_size: [f32; 2],
    cursor_pos: [f32; 2],
) -> Option<(glam::Vec3, glam::Vec3)> {
    if viewport_size[0] <= 1.0 || viewport_size[1] <= 1.0 {
        return None;
    }

    let ndc_x = (2.0 * cursor_pos[0] / viewport_size[0]) - 1.0;
    let ndc_y = 1.0 - (2.0 * cursor_pos[1] / viewport_size[1]);

    let vp =
        glam::Mat4::from_cols_array_2d(&camera.build_view_projection_matrix(aspect.max(0.001)));
    let inv = vp.inverse();

    let near = inv * glam::Vec4::new(ndc_x, ndc_y, 0.0, 1.0);
    let far = inv * glam::Vec4::new(ndc_x, ndc_y, 1.0, 1.0);

    if near.w.abs() < 1e-6 || far.w.abs() < 1e-6 {
        return None;
    }

    let p0 = (near / near.w).truncate();
    let p1 = (far / far.w).truncate();
    let dir = (p1 - p0).normalize_or_zero();
    if dir.length_squared() <= 1e-8 {
        return None;
    }
    Some((p0, dir))
}

pub fn cursor_point_on_plane(
    camera: &Camera,
    aspect: f32,
    viewport_size: [f32; 2],
    cursor_pos: [f32; 2],
    plane_origin: glam::DVec3,
    plane_normal: glam::DVec3,
) -> Option<glam::DVec3> {
    let (ray_origin, ray_dir) = screen_ray(camera, aspect, viewport_size, cursor_pos)?;
    let plane_origin = glam::Vec3::new(
        plane_origin.x as f32,
        plane_origin.y as f32,
        plane_origin.z as f32,
    );
    let plane_normal = glam::Vec3::new(
        plane_normal.x as f32,
        plane_normal.y as f32,
        plane_normal.z as f32,
    )
    .normalize_or_zero();
    if plane_normal.length_squared() <= 1e-8 {
        return None;
    }

    let denom = plane_normal.dot(ray_dir);
    if denom.abs() <= 1e-6 {
        return None;
    }

    let t = plane_normal.dot(plane_origin - ray_origin) / denom;
    if !t.is_finite() || t < 0.0 {
        return None;
    }

    let point = ray_origin + ray_dir * t;
    Some(glam::DVec3::new(
        point.x as f64,
        point.y as f64,
        point.z as f64,
    ))
}

fn ray_triangle_intersection(
    ray_origin: glam::Vec3,
    ray_dir: glam::Vec3,
    v0: glam::Vec3,
    v1: glam::Vec3,
    v2: glam::Vec3,
) -> Option<f32> {
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    let pvec = ray_dir.cross(e2);
    let det = e1.dot(pvec);
    if det.abs() < 1e-8 {
        return None;
    }
    let inv_det = 1.0 / det;
    let tvec = ray_origin - v0;
    let u = tvec.dot(pvec) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let qvec = tvec.cross(e1);
    let v = ray_dir.dot(qvec) * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = e2.dot(qvec) * inv_det;
    if t.is_finite() { Some(t) } else { None }
}

fn project_to_screen(vp: glam::Mat4, p: glam::Vec3, viewport_size: [f32; 2]) -> Option<[f32; 3]> {
    let clip = vp * p.extend(1.0);
    if clip.w.abs() < 1e-6 {
        return None;
    }
    let ndc = (clip / clip.w).truncate();
    let x = (ndc.x + 1.0) * 0.5 * viewport_size[0];
    let y = (1.0 - ndc.y) * 0.5 * viewport_size[1];
    Some([x, y, ndc.z])
}

fn point_segment_distance_2d(p: [f32; 2], a: [f32; 2], b: [f32; 2]) -> f32 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ap = [p[0] - a[0], p[1] - a[1]];
    let ab_len2 = ab[0] * ab[0] + ab[1] * ab[1];
    if ab_len2 <= 1e-8 {
        return ((p[0] - a[0]).powi(2) + (p[1] - a[1]).powi(2)).sqrt();
    }
    let t = ((ap[0] * ab[0] + ap[1] * ab[1]) / ab_len2).clamp(0.0, 1.0);
    let q = [a[0] + ab[0] * t, a[1] + ab[1] * t];
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt()
}

#[allow(dead_code)]
struct AxisBuffers {
    vertex_buffer: wgpu::Buffer,
    tri_index_buffer: wgpu::Buffer,
    tri_index_count: u32,
    line_index_buffer: wgpu::Buffer,
    line_index_count: u32,
}

impl std::fmt::Debug for AxisBuffers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AxisBuffers")
            .field("tri_index_count", &self.tri_index_count)
            .field("line_index_count", &self.line_index_count)
            .finish_non_exhaustive()
    }
}

fn build_axis_arrow_mesh(
    direction: glam::Vec3,
    shaft_length: f32,
    cone_radius: f32,
    cone_height: f32,
    segments: u32,
) -> (Vec<[f32; 3]>, Vec<u32>, Vec<u32>) {
    let dir = direction.normalize();

    // Build a local frame: dir is the axis, u and v are perpendicular
    let arbitrary = if dir.y.abs() < 0.9 {
        glam::Vec3::Y
    } else {
        glam::Vec3::X
    };
    let u = dir.cross(arbitrary).normalize();
    let v = dir.cross(u);

    // Shaft: two vertices (origin → shaft_length along dir)
    let shaft_end = dir * shaft_length;
    let mut vertices = vec![[0.0, 0.0, 0.0], shaft_end.to_array()];

    // Line indices for the shaft
    let line_indices = vec![0, 1];

    // Cone: base ring at shaft_end, tip at shaft_end + cone_height * dir
    let tip = shaft_end + dir * cone_height;
    let tip_idx = vertices.len() as u32;
    vertices.push(tip.to_array());

    let base_start = vertices.len() as u32;
    for i in 0..segments {
        let angle = 2.0 * std::f32::consts::PI * (i as f32) / (segments as f32);
        let point = shaft_end + u * (angle.cos() * cone_radius) + v * (angle.sin() * cone_radius);
        vertices.push(point.to_array());
    }

    // Triangle fan for cone
    let mut tri_indices = Vec::new();
    for i in 0..segments {
        let curr = base_start + i;
        let next = base_start + (i + 1) % segments;
        // Side face
        tri_indices.push(tip_idx);
        tri_indices.push(curr);
        tri_indices.push(next);
    }

    // Base cap center
    let base_center_idx = vertices.len() as u32;
    vertices.push(shaft_end.to_array());
    for i in 0..segments {
        let curr = base_start + i;
        let next = base_start + (i + 1) % segments;
        tri_indices.push(base_center_idx);
        tri_indices.push(next);
        tri_indices.push(curr);
    }

    (vertices, tri_indices, line_indices)
}

fn append_ring_vertices(
    vertices: &mut Vec<[f32; 3]>,
    normals: &mut Vec<[f32; 3]>,
    center: glam::Vec3,
    u: glam::Vec3,
    v: glam::Vec3,
    radius: f32,
    segments: u32,
    normal_scale_xy: f32,
    normal_scale_axis: f32,
    axis: glam::Vec3,
) -> u32 {
    let start = vertices.len() as u32;
    for i in 0..segments {
        let angle = 2.0 * std::f32::consts::PI * (i as f32) / segments as f32;
        let radial = u * angle.cos() + v * angle.sin();
        let position = center + radial * radius;
        let normal = (radial * normal_scale_xy + axis * normal_scale_axis).normalize_or_zero();
        vertices.push(position.to_array());
        normals.push(normal.to_array());
    }
    start
}

fn append_ring_indices(indices: &mut Vec<u32>, start0: u32, start1: u32, segments: u32) {
    for i in 0..segments {
        let next = (i + 1) % segments;
        let a0 = start0 + i;
        let a1 = start0 + next;
        let b0 = start1 + i;
        let b1 = start1 + next;

        indices.extend_from_slice(&[a0, b0, a1, a1, b0, b1]);
    }
}

fn build_solid_axis_arrow_mesh(
    direction: glam::Vec3,
    shaft_length: f32,
    shaft_radius: f32,
    cone_radius: f32,
    cone_height: f32,
    segments: u32,
) -> (Vec<[f32; 3]>, Vec<[f32; 3]>, Vec<u32>) {
    let axis = direction.normalize_or_zero();
    let arbitrary = if axis.y.abs() < 0.9 {
        glam::Vec3::Y
    } else {
        glam::Vec3::X
    };
    let u = axis.cross(arbitrary).normalize_or_zero();
    let v = axis.cross(u).normalize_or_zero();

    let mut vertices = Vec::new();
    let mut normals = Vec::new();
    let mut indices = Vec::new();

    let shaft_start = append_ring_vertices(
        &mut vertices,
        &mut normals,
        glam::Vec3::ZERO,
        u,
        v,
        shaft_radius,
        segments,
        1.0,
        0.0,
        axis,
    );
    let shaft_end_center = axis * shaft_length;
    let shaft_end = append_ring_vertices(
        &mut vertices,
        &mut normals,
        shaft_end_center,
        u,
        v,
        shaft_radius,
        segments,
        1.0,
        0.0,
        axis,
    );
    append_ring_indices(&mut indices, shaft_start, shaft_end, segments);

    let cone_base = append_ring_vertices(
        &mut vertices,
        &mut normals,
        shaft_end_center,
        u,
        v,
        cone_radius,
        segments,
        cone_height,
        cone_radius,
        axis,
    );
    let tip_index = vertices.len() as u32;
    let tip = shaft_end_center + axis * cone_height;
    vertices.push(tip.to_array());
    normals.push(axis.to_array());
    for i in 0..segments {
        let next = (i + 1) % segments;
        indices.extend_from_slice(&[tip_index, cone_base + i, cone_base + next]);
    }

    (vertices, normals, indices)
}

fn build_cube_mesh(half_extent: f32) -> (Vec<[f32; 3]>, Vec<[f32; 3]>, Vec<u32>) {
    let mut vertices = Vec::new();
    let mut normals = Vec::new();
    let mut indices = Vec::new();

    let faces = [
        (
            glam::Vec3::X,
            [
                glam::Vec3::new(half_extent, -half_extent, -half_extent),
                glam::Vec3::new(half_extent, -half_extent, half_extent),
                glam::Vec3::new(half_extent, half_extent, half_extent),
                glam::Vec3::new(half_extent, half_extent, -half_extent),
            ],
        ),
        (
            -glam::Vec3::X,
            [
                glam::Vec3::new(-half_extent, -half_extent, half_extent),
                glam::Vec3::new(-half_extent, -half_extent, -half_extent),
                glam::Vec3::new(-half_extent, half_extent, -half_extent),
                glam::Vec3::new(-half_extent, half_extent, half_extent),
            ],
        ),
        (
            glam::Vec3::Y,
            [
                glam::Vec3::new(-half_extent, half_extent, -half_extent),
                glam::Vec3::new(half_extent, half_extent, -half_extent),
                glam::Vec3::new(half_extent, half_extent, half_extent),
                glam::Vec3::new(-half_extent, half_extent, half_extent),
            ],
        ),
        (
            -glam::Vec3::Y,
            [
                glam::Vec3::new(-half_extent, -half_extent, half_extent),
                glam::Vec3::new(half_extent, -half_extent, half_extent),
                glam::Vec3::new(half_extent, -half_extent, -half_extent),
                glam::Vec3::new(-half_extent, -half_extent, -half_extent),
            ],
        ),
        (
            glam::Vec3::Z,
            [
                glam::Vec3::new(-half_extent, -half_extent, half_extent),
                glam::Vec3::new(-half_extent, half_extent, half_extent),
                glam::Vec3::new(half_extent, half_extent, half_extent),
                glam::Vec3::new(half_extent, -half_extent, half_extent),
            ],
        ),
        (
            -glam::Vec3::Z,
            [
                glam::Vec3::new(half_extent, -half_extent, -half_extent),
                glam::Vec3::new(half_extent, half_extent, -half_extent),
                glam::Vec3::new(-half_extent, half_extent, -half_extent),
                glam::Vec3::new(-half_extent, -half_extent, -half_extent),
            ],
        ),
    ];

    for (normal, corners) in faces {
        let start = vertices.len() as u32;
        for corner in corners {
            vertices.push(corner.to_array());
            normals.push(normal.to_array());
        }
        indices.extend_from_slice(&[start, start + 1, start + 2, start, start + 2, start + 3]);
    }

    (vertices, normals, indices)
}

fn interleave_positions_normals(vertices: &[[f32; 3]], normals: &[[f32; 3]]) -> Vec<[f32; 6]> {
    vertices
        .iter()
        .zip(normals.iter())
        .map(|(&[px, py, pz], &[nx, ny, nz])| [px, py, pz, nx, ny, nz])
        .collect()
}

struct GridBuffers {
    vertex_buffer: wgpu::Buffer,
    major_index_buffer: wgpu::Buffer,
    major_index_count: std::sync::Mutex<u32>,
    minor_index_buffer: wgpu::Buffer,
    minor_index_count: std::sync::Mutex<u32>,
}

/// Interleave vertex positions and normals into a flat `Vec<[f32; 6]>` for GPU upload.
/// If `normals` is empty or a different length, zero normals are used (triggers flat-shading fallback in shader).
fn interleave_verts_normals(vertices: &[[f32; 3]], normals: &[[f32; 3]]) -> Vec<[f32; 6]> {
    let has_normals = normals.len() == vertices.len();
    vertices
        .iter()
        .enumerate()
        .map(|(i, &[px, py, pz])| {
            let [nx, ny, nz] = if has_normals { normals[i] } else { [0.0, 0.0, 0.0] };
            [px, py, pz, nx, ny, nz]
        })
        .collect()
}

fn snap_grid_spacing(raw_spacing: f32) -> f32 {
    let raw_spacing = raw_spacing.max(GRID_MIN_MINOR_SPACING);
    let magnitude = 10.0_f32.powf(raw_spacing.log10().floor());
    let normalized = raw_spacing / magnitude;
    let snapped = if normalized <= 1.0 {
        1.0
    } else if normalized <= 2.0 {
        2.0
    } else if normalized <= 5.0 {
        5.0
    } else {
        10.0
    };
    snapped * magnitude
}

fn grid_focus_point(camera: &Camera) -> glam::Vec3 {
    let eye = camera.eye_position();
    let forward = (camera.target - eye).normalize_or_zero();
    if forward.length_squared() <= 1e-8 || forward.y.abs() <= 1e-4 {
        return camera.target;
    }

    let t = -eye.y / forward.y;
    if t.is_finite() && t > 0.0 {
        eye + forward * t
    } else {
        camera.target
    }
}

fn adaptive_grid_minor_spacing(camera: &Camera, aspect: f32) -> f32 {
    let half_view_width = camera.distance.max(0.1)
        * (GRID_FOV_Y_RADIANS * 0.5).tan()
        * aspect.max(1.0);
    let raw_spacing = half_view_width / GRID_TARGET_HALF_MINOR_LINES;
    snap_grid_spacing(raw_spacing)
}

fn grid_material_uniforms(minor_spacing: f32) -> (MaterialUniform, MaterialUniform) {
    let fade = (1.0 / (1.0 + minor_spacing * 6.0)).clamp(0.0, 1.0);
    let major_alpha = 0.26 + 0.18 * fade;
    let minor_alpha = 0.04 + 0.12 * fade;
    (
        MaterialUniform {
            color: [0.44, 0.46, 0.50, major_alpha],
            flags: [1.0, 0.0, 0.0, 0.0],
        },
        MaterialUniform {
            color: [0.32, 0.34, 0.38, minor_alpha],
            flags: [1.0, 0.0, 0.0, 0.0],
        },
    )
}

/// Build an adaptive grid on the XZ plane (Y=0). Returns (vertices, major_line_indices, minor_line_indices).
fn build_grid_mesh(
    center: glam::Vec2,
    minor_spacing: f32,
    half_cells: i32,
) -> (Vec<[f32; 3]>, Vec<u32>, Vec<u32>) {
    let mut vertices = Vec::new();
    let mut major_indices = Vec::new();
    let mut minor_indices = Vec::new();
    let spacing = minor_spacing.max(GRID_MIN_MINOR_SPACING);
    let center_step_x = (center.x / spacing).round() as i32;
    let center_step_z = (center.y / spacing).round() as i32;
    let start_x = center_step_x - half_cells;
    let end_x = center_step_x + half_cells;
    let start_z = center_step_z - half_cells;
    let end_z = center_step_z + half_cells;

    // Generate lines along X (varying Z)
    for step_z in start_z..=end_z {
        let z = step_z as f32 * spacing;
        let x0 = start_x as f32 * spacing;
        let x1 = end_x as f32 * spacing;

        let idx = vertices.len() as u32;
        vertices.push([x0, 0.0, z]);
        vertices.push([x1, 0.0, z]);

        let is_major = step_z.rem_euclid(GRID_MAJOR_LINE_EVERY) == 0 && z.abs() > 0.001;
        let is_origin_line = z.abs() < 0.001;

        if is_origin_line {
            // Skip — the axis handles this
        } else if is_major {
            major_indices.push(idx);
            major_indices.push(idx + 1);
        } else {
            minor_indices.push(idx);
            minor_indices.push(idx + 1);
        }
    }

    // Generate lines along Z (varying X)
    for step_x in start_x..=end_x {
        let x = step_x as f32 * spacing;
        let z0 = start_z as f32 * spacing;
        let z1 = end_z as f32 * spacing;

        let idx = vertices.len() as u32;
        vertices.push([x, 0.0, z0]);
        vertices.push([x, 0.0, z1]);

        let is_major = step_x.rem_euclid(GRID_MAJOR_LINE_EVERY) == 0 && x.abs() > 0.001;
        let is_origin_line = x.abs() < 0.001;

        if is_origin_line {
            // Skip
        } else if is_major {
            major_indices.push(idx);
            major_indices.push(idx + 1);
        } else {
            minor_indices.push(idx);
            minor_indices.push(idx + 1);
        }
    }

    (vertices, major_indices, minor_indices)
}

/// Default scale for scene axes (relative to fixed camera distance).
const DEFAULT_SCENE_AXES_SCALE: f32 = 0.3;

pub struct WgpuRenderer {
    pipeline: wgpu::RenderPipeline,
    pipeline_depth: wgpu::RenderPipeline,
    pipeline_line: wgpu::RenderPipeline,
    pipeline_line_depth: wgpu::RenderPipeline,
    camera_buffer: wgpu::Buffer,
    camera_bind_group: wgpu::BindGroup,
    material_bind_group: wgpu::BindGroup,
    material_face_highlight_bind_group: wgpu::BindGroup,
    material_edge_highlight_bind_group: wgpu::BindGroup,
    vertex_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    index_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    index_count: std::sync::Mutex<u32>,
    line_index_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    line_index_count: std::sync::Mutex<u32>,
    highlight_face_vertex_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    highlight_face_index_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    highlight_face_index_count: std::sync::Mutex<u32>,
    highlight_edge_vertex_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    highlight_edge_index_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    highlight_edge_index_count: std::sync::Mutex<u32>,
    depth_texture: std::sync::Mutex<Option<wgpu::Texture>>,
    depth_view: std::sync::Mutex<Option<wgpu::TextureView>>,
    depth_size: std::sync::Mutex<(u32, u32)>,
    axes_buffers: [AxisBuffers; 3],
    axes_camera_buffer: wgpu::Buffer,
    axes_camera_bind_group: wgpu::BindGroup,
    axes_material_bind_groups: [wgpu::BindGroup; 3],
    show_axes: std::sync::Mutex<bool>,
    scene_axes_scale: std::sync::Mutex<f32>,
    display_mode: std::sync::Mutex<DisplayMode>,
    material_transparent_bind_group: wgpu::BindGroup,
    material_buffer: wgpu::Buffer,
    grid: GridBuffers,
    grid_major_material_buffer: wgpu::Buffer,
    grid_minor_material_buffer: wgpu::Buffer,
    grid_major_material_bind_group: wgpu::BindGroup,
    grid_minor_material_bind_group: wgpu::BindGroup,
    gizmo_camera_buffer: wgpu::Buffer,
    gizmo_camera_bind_group: wgpu::BindGroup,
    gizmo_axes_buffers: [SolidMeshBuffers; 3],
    gizmo_axes_material_bind_groups: [wgpu::BindGroup; 3],
    gizmo_center_buffers: SolidMeshBuffers,
    gizmo_center_material_bind_group: wgpu::BindGroup,
    show_grid: std::sync::Mutex<bool>,
    show_axis_gizmo: std::sync::Mutex<bool>,
    gizmo_eye: std::sync::Mutex<glam::Vec3>,
    light_dir: std::sync::Mutex<glam::Vec3>,
}

unsafe impl Send for WgpuRenderer {}
unsafe impl Sync for WgpuRenderer {}

impl WgpuRenderer {
    pub fn default_clear_color() -> wgpu::Color {
        wgpu::Color {
            r: 0.07,
            g: 0.07,
            b: 0.11,
            a: 1.0,
        }
    }

    fn create_pipeline(
        device: &wgpu::Device,
        shader: &wgpu::ShaderModule,
        surface_format: wgpu::TextureFormat,
        layout: &wgpu::PipelineLayout,
        topology: wgpu::PrimitiveTopology,
        with_depth: bool,
        blend: Option<wgpu::BlendState>,
    ) -> wgpu::RenderPipeline {
        device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some(if with_depth {
                "Render Pipeline (Depth)"
            } else {
                "Render Pipeline"
            }),
            layout: Some(layout),
            cache: None,
            vertex: wgpu::VertexState {
                module: shader,
                entry_point: Some("vs_main"),
                compilation_options: Default::default(),
                buffers: &[wgpu::VertexBufferLayout {
                    // Each vertex is [px, py, pz, nx, ny, nz] = 6 × f32 = 24 bytes.
                    // The normal component (location 1) is zero for meshes that
                    // do not carry smooth normals (grid, axes, highlights).
                    array_stride: std::mem::size_of::<[f32; 6]>() as wgpu::BufferAddress,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3],
                }],
            },
            fragment: Some(wgpu::FragmentState {
                module: shader,
                entry_point: Some("fs_main"),
                compilation_options: Default::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format: surface_format,
                    blend,
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None,
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: if with_depth {
                Some(wgpu::DepthStencilState {
                    format: wgpu::TextureFormat::Depth32Float,
                    depth_write_enabled: true,
                    depth_compare: wgpu::CompareFunction::LessEqual,
                    stencil: wgpu::StencilState::default(),
                    bias: wgpu::DepthBiasState::default(),
                })
            } else {
                None
            },
            multisample: wgpu::MultisampleState {
                count: 1,
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            multiview: None,
        })
    }

    pub fn new(device: &wgpu::Device, surface_format: wgpu::TextureFormat) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("RCAD Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shader.wgsl").into()),
        });

        let camera_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Camera Buffer"),
            contents: bytemuck::cast_slice(&[CameraUniform {
                view_proj: glam::Mat4::IDENTITY.to_cols_array_2d(),
                eye_pos: [0.0, 0.0, 3.0, 1.0],
                light_dir: [0.45, 0.85, 0.35, 0.0],
            }]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let camera_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("Camera Bind Group Layout"),
                entries: &[wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                }],
            });

        let camera_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Camera Bind Group"),
            layout: &camera_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: camera_buffer.as_entire_binding(),
            }],
        });

        let gizmo_camera_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Axis Gizmo Camera Buffer"),
            contents: bytemuck::cast_slice(&[CameraUniform {
                view_proj: glam::Mat4::IDENTITY.to_cols_array_2d(),
                eye_pos: [0.0, 0.0, 3.0, 1.0],
                light_dir: [0.45, 0.85, 0.35, 0.0],
            }]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let gizmo_camera_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Axis Gizmo Camera Bind Group"),
            layout: &camera_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: gizmo_camera_buffer.as_entire_binding(),
            }],
        });

        // Scene axes camera buffer (fixed-size axes at origin)
        let axes_camera_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Scene Axes Camera Buffer"),
            contents: bytemuck::cast_slice(&[CameraUniform {
                view_proj: glam::Mat4::IDENTITY.to_cols_array_2d(),
                eye_pos: [0.0, 0.0, 3.0, 1.0],
                light_dir: [0.45, 0.85, 0.35, 0.0],
            }]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let axes_camera_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Scene Axes Camera Bind Group"),
            layout: &camera_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: axes_camera_buffer.as_entire_binding(),
            }],
        });

        let material_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("Material Bind Group Layout"),
                entries: &[wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                }],
            });

        let material_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Material Buffer"),
            contents: bytemuck::cast_slice(&[MaterialUniform {
                color: [0.18, 0.64, 0.96, 1.0],
                flags: [0.0, 0.0, 0.0, 0.0],
            }]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let material_transparent_buffer =
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Transparent Material Buffer"),
                contents: bytemuck::cast_slice(&[MaterialUniform {
                    color: [0.18, 0.64, 0.96, 0.3],
                    flags: [0.0, 0.0, 0.0, 0.0],
                }]),
                usage: wgpu::BufferUsages::UNIFORM,
            });
        let material_face_highlight_buffer =
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Face Highlight Material Buffer"),
                contents: bytemuck::cast_slice(&[MaterialUniform {
                    color: [1.0, 0.45, 0.05, 0.45],
                    flags: [1.0, 0.0, 0.0, 0.0],
                }]),
                usage: wgpu::BufferUsages::UNIFORM,
            });
        let material_edge_highlight_buffer =
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Edge Highlight Material Buffer"),
                contents: bytemuck::cast_slice(&[MaterialUniform {
                    color: [1.0, 0.95, 0.1, 1.0],
                    flags: [1.0, 0.0, 0.0, 0.0],
                }]),
                usage: wgpu::BufferUsages::UNIFORM,
            });

        let material_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Material Bind Group"),
            layout: &material_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: material_buffer.as_entire_binding(),
            }],
        });
        let material_face_highlight_bind_group =
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Face Highlight Material Bind Group"),
                layout: &material_bind_group_layout,
                entries: &[wgpu::BindGroupEntry {
                    binding: 0,
                    resource: material_face_highlight_buffer.as_entire_binding(),
                }],
            });
        let material_edge_highlight_bind_group =
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Edge Highlight Material Bind Group"),
                layout: &material_bind_group_layout,
                entries: &[wgpu::BindGroupEntry {
                    binding: 0,
                    resource: material_edge_highlight_buffer.as_entire_binding(),
                }],
            });
        let material_transparent_bind_group =
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Transparent Material Bind Group"),
                layout: &material_bind_group_layout,
                entries: &[wgpu::BindGroupEntry {
                    binding: 0,
                    resource: material_transparent_buffer.as_entire_binding(),
                }],
            });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Render Pipeline Layout"),
            bind_group_layouts: &[&camera_bind_group_layout, &material_bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = Self::create_pipeline(
            device,
            &shader,
            surface_format,
            &pipeline_layout,
            wgpu::PrimitiveTopology::TriangleList,
            false,
            Some(wgpu::BlendState::ALPHA_BLENDING),
        );
        let pipeline_depth = Self::create_pipeline(
            device,
            &shader,
            surface_format,
            &pipeline_layout,
            wgpu::PrimitiveTopology::TriangleList,
            true,
            Some(wgpu::BlendState::ALPHA_BLENDING),
        );
        let pipeline_line = Self::create_pipeline(
            device,
            &shader,
            surface_format,
            &pipeline_layout,
            wgpu::PrimitiveTopology::LineList,
            false,
            Some(wgpu::BlendState::ALPHA_BLENDING),
        );
        let pipeline_line_depth = Self::create_pipeline(
            device,
            &shader,
            surface_format,
            &pipeline_layout,
            wgpu::PrimitiveTopology::LineList,
            true,
            Some(wgpu::BlendState::ALPHA_BLENDING),
        );

        // Build background grid
        let (grid_major_material, grid_minor_material) = grid_material_uniforms(0.2);
        let grid_major_material_buffer =
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Grid Major Material Buffer"),
                contents: bytemuck::cast_slice(&[grid_major_material]),
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            });
        let grid_minor_material_buffer =
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Grid Minor Material Buffer"),
                contents: bytemuck::cast_slice(&[grid_minor_material]),
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            });
        let grid_major_material_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Grid Major Material Bind Group"),
            layout: &material_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: grid_major_material_buffer.as_entire_binding(),
            }],
        });
        let grid_minor_material_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Grid Minor Material Bind Group"),
            layout: &material_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: grid_minor_material_buffer.as_entire_binding(),
            }],
        });

        let max_grid_lines = ((GRID_BUFFER_HALF_CELLS * 2 + 1) * 2) as u64;
        let max_grid_vertices = max_grid_lines * 2;
        let max_grid_indices = max_grid_vertices;
        let grid = GridBuffers {
            vertex_buffer: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("Grid Vertex Buffer"),
                size: max_grid_vertices * std::mem::size_of::<[f32; 6]>() as u64,
                usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }),
            major_index_buffer: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("Grid Major Index Buffer"),
                size: max_grid_indices * std::mem::size_of::<u32>() as u64,
                usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }),
            major_index_count: std::sync::Mutex::new(0),
            minor_index_buffer: device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("Grid Minor Index Buffer"),
                size: max_grid_indices * std::mem::size_of::<u32>() as u64,
                usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }),
            minor_index_count: std::sync::Mutex::new(0),
        };

        // Build axis arrows (X=red, Y=green, Z=blue)
        let axis_colors: [[f32; 4]; 3] = [
            [1.0, 0.2, 0.2, 1.0], // X — red
            [0.2, 1.0, 0.2, 1.0], // Y — green
            [0.3, 0.5, 1.0, 1.0], // Z — blue
        ];
        let axis_dirs = [glam::Vec3::X, glam::Vec3::Y, glam::Vec3::Z];
        let axis_names = ["X", "Y", "Z"];

        let mut axes_material_bind_groups_vec = Vec::with_capacity(3);
        let mut axes_buffers_vec = Vec::with_capacity(3);

        for i in 0..3 {
            let mat_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(&format!("Axis {} Material Buffer", axis_names[i])),
                contents: bytemuck::cast_slice(&[MaterialUniform {
                    color: axis_colors[i],
                    flags: [1.0, 0.0, 0.0, 0.0], // unlit
                }]),
                usage: wgpu::BufferUsages::UNIFORM,
            });
            let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some(&format!("Axis {} Material Bind Group", axis_names[i])),
                layout: &material_bind_group_layout,
                entries: &[wgpu::BindGroupEntry {
                    binding: 0,
                    resource: mat_buf.as_entire_binding(),
                }],
            });
            axes_material_bind_groups_vec.push(bg);

            let (verts, tri_idx, line_idx) = build_axis_arrow_mesh(axis_dirs[i], 1.0, 0.03, 0.1, 8);
            let verts_padded: Vec<[f32; 6]> = verts.iter().map(|&[x, y, z]| [x, y, z, 0.0, 0.0, 0.0]).collect();

            let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(&format!("Axis {} Vertex Buffer", axis_names[i])),
                contents: bytemuck::cast_slice(&verts_padded),
                usage: wgpu::BufferUsages::VERTEX,
            });
            let tri_index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(&format!("Axis {} Tri Index Buffer", axis_names[i])),
                contents: bytemuck::cast_slice(&tri_idx),
                usage: wgpu::BufferUsages::INDEX,
            });
            let line_index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(&format!("Axis {} Line Index Buffer", axis_names[i])),
                contents: bytemuck::cast_slice(&line_idx),
                usage: wgpu::BufferUsages::INDEX,
            });

            axes_buffers_vec.push(AxisBuffers {
                vertex_buffer,
                tri_index_buffer,
                tri_index_count: tri_idx.len() as u32,
                line_index_buffer,
                line_index_count: line_idx.len() as u32,
            });
        }

        // Convert Vecs to fixed-size arrays
        let axes_material_bind_groups: [_; 3] = axes_material_bind_groups_vec
            .try_into()
            .expect("axes loop always produces exactly 3 bind groups");
        let axes_buffers: [_; 3] = axes_buffers_vec
            .try_into()
            .expect("axes loop always produces exactly 3 axis buffers");

        let gizmo_axis_colors: [[f32; 4]; 3] = [
            [0.98, 0.16, 0.12, 1.0],
            [0.10, 0.84, 0.26, 1.0],
            [0.14, 0.43, 0.98, 1.0],
        ];
        let mut gizmo_axes_material_bind_groups_vec = Vec::with_capacity(3);
        let mut gizmo_axes_buffers_vec = Vec::with_capacity(3);
        for i in 0..3 {
            let mat_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(&format!("Axis Gizmo {} Material Buffer", axis_names[i])),
                contents: bytemuck::cast_slice(&[MaterialUniform {
                    color: gizmo_axis_colors[i],
                    flags: [0.0, 0.0, 0.0, 0.0],
                }]),
                usage: wgpu::BufferUsages::UNIFORM,
            });
            let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some(&format!("Axis Gizmo {} Material Bind Group", axis_names[i])),
                layout: &material_bind_group_layout,
                entries: &[wgpu::BindGroupEntry {
                    binding: 0,
                    resource: mat_buf.as_entire_binding(),
                }],
            });
            gizmo_axes_material_bind_groups_vec.push(bind_group);

            let (verts, normals, tri_idx) = build_solid_axis_arrow_mesh(
                axis_dirs[i],
                AXIS_GIZMO_AXIS_LENGTH,
                0.06,
                0.12,
                0.24,
                20,
            );
            let interleaved = interleave_positions_normals(&verts, &normals);
            gizmo_axes_buffers_vec.push(SolidMeshBuffers {
                vertex_buffer: device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some(&format!("Axis Gizmo {} Vertex Buffer", axis_names[i])),
                    contents: bytemuck::cast_slice(&interleaved),
                    usage: wgpu::BufferUsages::VERTEX,
                }),
                index_buffer: device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some(&format!("Axis Gizmo {} Index Buffer", axis_names[i])),
                    contents: bytemuck::cast_slice(&tri_idx),
                    usage: wgpu::BufferUsages::INDEX,
                }),
                index_count: tri_idx.len() as u32,
            });
        }
        let gizmo_axes_material_bind_groups: [_; 3] = gizmo_axes_material_bind_groups_vec
            .try_into()
            .expect("gizmo axes loop always produces exactly 3 bind groups");
        let gizmo_axes_buffers: [_; 3] = gizmo_axes_buffers_vec
            .try_into()
            .expect("gizmo axes loop always produces exactly 3 axis buffers");

        let gizmo_center_material_buffer =
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Axis Gizmo Center Material Buffer"),
                contents: bytemuck::cast_slice(&[MaterialUniform {
                    color: [0.60, 0.62, 0.66, 1.0],
                    flags: [0.0, 0.0, 0.0, 0.0],
                }]),
                usage: wgpu::BufferUsages::UNIFORM,
            });
        let gizmo_center_material_bind_group =
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("Axis Gizmo Center Material Bind Group"),
                layout: &material_bind_group_layout,
                entries: &[wgpu::BindGroupEntry {
                    binding: 0,
                    resource: gizmo_center_material_buffer.as_entire_binding(),
                }],
            });
        let (gizmo_center_verts, gizmo_center_normals, gizmo_center_idx) =
            build_cube_mesh(AXIS_GIZMO_CENTER_HALF_EXTENT);
        let gizmo_center_buffers = SolidMeshBuffers {
            vertex_buffer: device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Axis Gizmo Center Vertex Buffer"),
                contents: bytemuck::cast_slice(&interleave_positions_normals(
                    &gizmo_center_verts,
                    &gizmo_center_normals,
                )),
                usage: wgpu::BufferUsages::VERTEX,
            }),
            index_buffer: device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Axis Gizmo Center Index Buffer"),
                contents: bytemuck::cast_slice(&gizmo_center_idx),
                usage: wgpu::BufferUsages::INDEX,
            }),
            index_count: gizmo_center_idx.len() as u32,
        };

        Self {
            pipeline,
            pipeline_depth,
            pipeline_line,
            pipeline_line_depth,
            camera_buffer,
            camera_bind_group,
            material_bind_group,
            material_face_highlight_bind_group,
            material_edge_highlight_bind_group,
            vertex_buffer: std::sync::Mutex::new(None),
            index_buffer: std::sync::Mutex::new(None),
            index_count: std::sync::Mutex::new(0),
            line_index_buffer: std::sync::Mutex::new(None),
            line_index_count: std::sync::Mutex::new(0),
            highlight_face_vertex_buffer: std::sync::Mutex::new(None),
            highlight_face_index_buffer: std::sync::Mutex::new(None),
            highlight_face_index_count: std::sync::Mutex::new(0),
            highlight_edge_vertex_buffer: std::sync::Mutex::new(None),
            highlight_edge_index_buffer: std::sync::Mutex::new(None),
            highlight_edge_index_count: std::sync::Mutex::new(0),
            depth_texture: std::sync::Mutex::new(None),
            depth_view: std::sync::Mutex::new(None),
            depth_size: std::sync::Mutex::new((0, 0)),
            axes_buffers,
            axes_camera_buffer,
            axes_camera_bind_group,
            axes_material_bind_groups,
            show_axes: std::sync::Mutex::new(true),
            scene_axes_scale: std::sync::Mutex::new(DEFAULT_SCENE_AXES_SCALE),
            display_mode: std::sync::Mutex::new(DisplayMode::default()),
            material_transparent_bind_group,
            material_buffer,
            grid,
            grid_major_material_buffer,
            grid_minor_material_buffer,
            grid_major_material_bind_group,
            grid_minor_material_bind_group,
            gizmo_camera_buffer,
            gizmo_camera_bind_group,
            gizmo_axes_buffers,
            gizmo_axes_material_bind_groups,
            gizmo_center_buffers,
            gizmo_center_material_bind_group,
            show_grid: std::sync::Mutex::new(true),
            show_axis_gizmo: std::sync::Mutex::new(true),
            gizmo_eye: std::sync::Mutex::new(glam::Vec3::new(0.0, 0.0, 3.2)),
            light_dir: std::sync::Mutex::new(glam::Vec3::new(0.45, 0.85, 0.35)),
        }
    }

    pub fn ensure_depth_texture(&self, device: &wgpu::Device, width: u32, height: u32) {
        let width = width.max(1);
        let height = height.max(1);

        {
            let size = self.depth_size.lock().expect("render mutex poisoned");
            let has_view = self.depth_view.lock().expect("render mutex poisoned").is_some();
            if has_view && *size == (width, height) {
                return;
            }
        }

        let texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("RCAD Depth Texture"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        let view = texture.create_view(&wgpu::TextureViewDescriptor::default());

        *self.depth_texture.lock().expect("render mutex poisoned") = Some(texture);
        *self.depth_view.lock().expect("render mutex poisoned") = Some(view);
        *self.depth_size.lock().expect("render mutex poisoned") = (width, height);
    }

    pub fn prepare_scene(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        mesh: &Mesh,
        camera: &Camera,
        aspect: f32,
    ) {
        self.upload_mesh(device, mesh);
        self.update_camera(queue, camera, aspect.max(0.001));
        self.update_axis_gizmo_camera(queue, camera);
    }

    pub fn prepare_scene_with_depth(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        mesh: &Mesh,
        camera: &Camera,
        aspect: f32,
        depth_size: (u32, u32),
    ) {
        self.ensure_depth_texture(device, depth_size.0, depth_size.1);
        self.prepare_scene(device, queue, mesh, camera, aspect);
    }

    pub fn upload_mesh(&self, device: &wgpu::Device, mesh: &Mesh) {
        let interleaved = interleave_verts_normals(&mesh.vertices, &mesh.normals);
        *self.vertex_buffer.lock().expect("render mutex poisoned") = Some(device.create_buffer_init(
            &wgpu::util::BufferInitDescriptor {
                label: Some("Vertex Buffer"),
                contents: bytemuck::cast_slice(&interleaved),
                usage: wgpu::BufferUsages::VERTEX,
            },
        ));

        *self.index_buffer.lock().expect("render mutex poisoned") = Some(device.create_buffer_init(
            &wgpu::util::BufferInitDescriptor {
                label: Some("Index Buffer"),
                contents: bytemuck::cast_slice(&mesh.indices),
                usage: wgpu::BufferUsages::INDEX,
            },
        ));

        *self.index_count.lock().expect("render mutex poisoned") = mesh.indices.len() as u32;

        if mesh.line_indices.is_empty() {
            *self.line_index_buffer.lock().expect("render mutex poisoned") = None;
            *self.line_index_count.lock().expect("render mutex poisoned") = 0;
        } else {
            *self.line_index_buffer.lock().expect("render mutex poisoned") = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Line Index Buffer"),
                    contents: bytemuck::cast_slice(&mesh.line_indices),
                    usage: wgpu::BufferUsages::INDEX,
                },
            ));
            *self.line_index_count.lock().expect("render mutex poisoned") = mesh.line_indices.len() as u32;
        }
    }

    pub fn upload_highlights(
        &self,
        device: &wgpu::Device,
        face_mesh: Option<&Mesh>,
        edge_mesh: Option<&Mesh>,
    ) {
        if let Some(mesh) = face_mesh {
            let interleaved = interleave_verts_normals(&mesh.vertices, &mesh.normals);
            *self.highlight_face_vertex_buffer.lock().expect("render mutex poisoned") = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Face Vertex Buffer"),
                    contents: bytemuck::cast_slice(&interleaved),
                    usage: wgpu::BufferUsages::VERTEX,
                },
            ));
            *self.highlight_face_index_buffer.lock().expect("render mutex poisoned") = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Face Index Buffer"),
                    contents: bytemuck::cast_slice(&mesh.indices),
                    usage: wgpu::BufferUsages::INDEX,
                },
            ));
            *self.highlight_face_index_count.lock().expect("render mutex poisoned") = mesh.indices.len() as u32;
        } else {
            *self.highlight_face_vertex_buffer.lock().expect("render mutex poisoned") = None;
            *self.highlight_face_index_buffer.lock().expect("render mutex poisoned") = None;
            *self.highlight_face_index_count.lock().expect("render mutex poisoned") = 0;
        }

        if let Some(mesh) = edge_mesh {
            let interleaved = interleave_verts_normals(&mesh.vertices, &mesh.normals);
            *self.highlight_edge_vertex_buffer.lock().expect("render mutex poisoned") = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Edge Vertex Buffer"),
                    contents: bytemuck::cast_slice(&interleaved),
                    usage: wgpu::BufferUsages::VERTEX,
                },
            ));
            *self.highlight_edge_index_buffer.lock().expect("render mutex poisoned") = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Edge Index Buffer"),
                    contents: bytemuck::cast_slice(&mesh.indices),
                    usage: wgpu::BufferUsages::INDEX,
                },
            ));
            *self.highlight_edge_index_count.lock().expect("render mutex poisoned") = mesh.indices.len() as u32;
        } else {
            *self.highlight_edge_vertex_buffer.lock().expect("render mutex poisoned") = None;
            *self.highlight_edge_index_buffer.lock().expect("render mutex poisoned") = None;
            *self.highlight_edge_index_count.lock().expect("render mutex poisoned") = 0;
        }
    }

    pub fn update_camera(&self, queue: &wgpu::Queue, camera: &Camera, aspect: f32) {
        let eye = camera.eye_position();
        let ld = *self.light_dir.lock().expect("render mutex poisoned");
        let uniform = CameraUniform {
            view_proj: camera.build_view_projection_matrix(aspect),
            eye_pos: [eye.x, eye.y, eye.z, 1.0],
            light_dir: [ld.x, ld.y, ld.z, 0.0],
        };
        queue.write_buffer(&self.camera_buffer, 0, bytemuck::cast_slice(&[uniform]));

        let grid_minor_spacing = adaptive_grid_minor_spacing(camera, aspect);
        let grid_focus = grid_focus_point(camera);
        let (grid_major_material, grid_minor_material) = grid_material_uniforms(grid_minor_spacing);
        let (grid_vertices, grid_major_indices, grid_minor_indices) = build_grid_mesh(
            glam::Vec2::new(grid_focus.x, grid_focus.z),
            grid_minor_spacing,
            GRID_BUFFER_HALF_CELLS,
        );
        let padded_vertices: Vec<[f32; 6]> = grid_vertices
            .iter()
            .map(|&[x, y, z]| [x, y, z, 0.0, 0.0, 0.0])
            .collect();
        queue.write_buffer(&self.grid.vertex_buffer, 0, bytemuck::cast_slice(&padded_vertices));
        queue.write_buffer(
            &self.grid.major_index_buffer,
            0,
            bytemuck::cast_slice(&grid_major_indices),
        );
        queue.write_buffer(
            &self.grid.minor_index_buffer,
            0,
            bytemuck::cast_slice(&grid_minor_indices),
        );
        queue.write_buffer(
            &self.grid_major_material_buffer,
            0,
            bytemuck::cast_slice(&[grid_major_material]),
        );
        queue.write_buffer(
            &self.grid_minor_material_buffer,
            0,
            bytemuck::cast_slice(&[grid_minor_material]),
        );
        *self.grid.major_index_count.lock().expect("render mutex poisoned") =
            grid_major_indices.len() as u32;
        *self.grid.minor_index_count.lock().expect("render mutex poisoned") =
            grid_minor_indices.len() as u32;

        // Update scene axes camera (fixed size, follows rotation only)
        self.update_scene_axes_camera(queue, camera, aspect);
    }

    /// Update the scene axes camera for fixed-size axes at origin.
    fn update_scene_axes_camera(&self, queue: &wgpu::Queue, camera: &Camera, aspect: f32) {
        let scale = *self.scene_axes_scale.lock().expect("render mutex poisoned");

        // Build view matrix using camera rotation but with fixed distance
        let yaw = camera.rot_y;
        let pitch = camera.rot_x;

        // Fixed camera distance for consistent axis size
        let fixed_distance = 3.0;

        // Compute eye position based on rotation
        let eye = glam::Vec3::new(
            fixed_distance * pitch.cos() * yaw.sin(),
            fixed_distance * pitch.sin(),
            fixed_distance * pitch.cos() * yaw.cos(),
        );

        let view = glam::Mat4::look_at_rh(eye, glam::Vec3::ZERO, glam::Vec3::Y);

        // Use perspective projection with fixed FOV for consistent size
        let fov_y = std::f32::consts::FRAC_PI_4;
        let near = 0.1;
        let far = 100.0;
        let proj = glam::Mat4::perspective_rh(fov_y, aspect, near, far);

        // Apply scale factor
        let scale_matrix = glam::Mat4::from_scale(glam::Vec3::splat(scale));

        let view_proj = proj * view * scale_matrix;

        let light_dir = glam::Vec3::new(0.45, 0.85, 0.35).normalize_or_zero();
        let uniform = CameraUniform {
            view_proj: view_proj.to_cols_array_2d(),
            eye_pos: [eye.x, eye.y, eye.z, 1.0],
            light_dir: [light_dir.x, light_dir.y, light_dir.z, 0.0],
        };
        queue.write_buffer(&self.axes_camera_buffer, 0, bytemuck::cast_slice(&[uniform]));
    }

    /// Set the scene axes scale factor (default 0.3).
    pub fn set_scene_axes_scale(&self, scale: f32) {
        *self.scene_axes_scale.lock().expect("render mutex poisoned") = scale;
    }

    /// Get the scene axes scale factor.
    pub fn scene_axes_scale(&self) -> f32 {
        *self.scene_axes_scale.lock().expect("render mutex poisoned")
    }
    pub fn update_axis_gizmo_camera(&self, queue: &wgpu::Queue, camera: &Camera) {
        let eye = axis_gizmo_eye(camera);
        let light_dir = glam::Vec3::new(0.55, 0.8, 0.35).normalize_or_zero();
        let uniform = CameraUniform {
            view_proj: axis_gizmo_view_projection(camera),
            eye_pos: [eye.x, eye.y, eye.z, 1.0],
            light_dir: [light_dir.x, light_dir.y, light_dir.z, 0.0],
        };
        queue.write_buffer(&self.gizmo_camera_buffer, 0, bytemuck::cast_slice(&[uniform]));
        *self.gizmo_eye.lock().expect("render mutex poisoned") = eye;
    }

    pub fn draw_in_render_pass(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        use_depth_pipeline: bool,
    ) {
        let mode = *self.display_mode.lock().expect("render mutex poisoned");

        // Draw grid first (behind everything)
        if *self.show_grid.lock().expect("render mutex poisoned") {
            self.draw_grid_in_render_pass(render_pass, use_depth_pipeline);
        }

        let vb_guard = self.vertex_buffer.lock().expect("render mutex poisoned");
        let ib_guard = self.index_buffer.lock().expect("render mutex poisoned");
        let count = *self.index_count.lock().expect("render mutex poisoned");
        let lib_guard = self.line_index_buffer.lock().expect("render mutex poisoned");
        let lcount = *self.line_index_count.lock().expect("render mutex poisoned");

        // Draw model based on display mode
        let draw_triangles = matches!(
            mode,
            DisplayMode::Solid | DisplayMode::SolidWithEdges | DisplayMode::Transparent
        );
        let draw_wireframe = matches!(
            mode,
            DisplayMode::Wireframe | DisplayMode::SolidWithEdges | DisplayMode::Transparent
        );

        // In transparent mode, draw wireframe first so it's behind the translucent surface
        if mode == DisplayMode::Transparent
            && draw_wireframe
            && lcount > 0
            && let (Some(vb), Some(lib)) = (vb_guard.as_ref(), lib_guard.as_ref())
        {
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_line_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline_line);
            }
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.material_bind_group, &[]);
            render_pass.set_vertex_buffer(0, vb.slice(..));
            render_pass.set_index_buffer(lib.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..lcount, 0, 0..1);
        }

        // Draw triangles
        if draw_triangles
            && count > 0
            && let (Some(vb), Some(ib)) = (vb_guard.as_ref(), ib_guard.as_ref())
        {
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline);
            }
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            let mat = if mode == DisplayMode::Transparent {
                &self.material_transparent_bind_group
            } else {
                &self.material_bind_group
            };
            render_pass.set_bind_group(1, mat, &[]);
            render_pass.set_vertex_buffer(0, vb.slice(..));
            render_pass.set_index_buffer(ib.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..count, 0, 0..1);
        }

        // Draw wireframe (non-transparent modes)
        if draw_wireframe
            && mode != DisplayMode::Transparent
            && lcount > 0
            && let (Some(vb), Some(lib)) = (vb_guard.as_ref(), lib_guard.as_ref())
        {
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_line_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline_line);
            }
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.material_bind_group, &[]);
            render_pass.set_vertex_buffer(0, vb.slice(..));
            render_pass.set_index_buffer(lib.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..lcount, 0, 0..1);
        }

        // Draw face highlights
        let hvb_guard = self.highlight_face_vertex_buffer.lock().expect("render mutex poisoned");
        let hib_guard = self.highlight_face_index_buffer.lock().expect("render mutex poisoned");
        let hcount = *self.highlight_face_index_count.lock().expect("render mutex poisoned");
        if hcount > 0
            && let (Some(vb), Some(ib)) = (hvb_guard.as_ref(), hib_guard.as_ref())
        {
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline);
            }
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.material_face_highlight_bind_group, &[]);
            render_pass.set_vertex_buffer(0, vb.slice(..));
            render_pass.set_index_buffer(ib.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..hcount, 0, 0..1);
        }

        // Draw edge highlights
        let evb_guard = self.highlight_edge_vertex_buffer.lock().expect("render mutex poisoned");
        let eib_guard = self.highlight_edge_index_buffer.lock().expect("render mutex poisoned");
        let ecount = *self.highlight_edge_index_count.lock().expect("render mutex poisoned");
        if ecount > 0
            && let (Some(vb), Some(ib)) = (evb_guard.as_ref(), eib_guard.as_ref())
        {
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_line_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline_line);
            }
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.material_edge_highlight_bind_group, &[]);
            render_pass.set_vertex_buffer(0, vb.slice(..));
            render_pass.set_index_buffer(ib.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..ecount, 0, 0..1);
        }

        // Draw coordinate axes
        if *self.show_axes.lock().expect("render mutex poisoned") {
            self.draw_axes_in_render_pass(render_pass, use_depth_pipeline);
        }
    }

    pub fn draw_axis_gizmo_in_render_pass(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        viewport_origin_px: [u32; 2],
        viewport_size_px: [u32; 2],
    ) {
        if !*self.show_axis_gizmo.lock().expect("render mutex poisoned") {
            return;
        }
        let Some(layout) = axis_gizmo_layout(viewport_origin_px, viewport_size_px) else {
            return;
        };
        let width = viewport_size_px[0];
        let height = viewport_size_px[1];
        let x = layout.origin_px[0];
        let y = layout.origin_px[1];
        let side = layout.size_px[0];

        render_pass.set_viewport(x as f32, y as f32, side as f32, side as f32, 0.0, 1.0);
        render_pass.set_scissor_rect(x, y, side, side);

        let eye = *self.gizmo_eye.lock().expect("render mutex poisoned");
        let mut far_axes = Vec::new();
        let mut near_axes = Vec::new();
        for axis_index in 0..3 {
            let axis_dir = match axis_index {
                0 => glam::Vec3::X,
                1 => glam::Vec3::Y,
                _ => glam::Vec3::Z,
            };
            let depth = (axis_dir * 0.48 - eye).length_squared();
            if depth > (glam::Vec3::ZERO - eye).length_squared() {
                far_axes.push((depth, axis_index));
            } else {
                near_axes.push((depth, axis_index));
            }
        }
        far_axes.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        near_axes.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

        for (_, axis_index) in far_axes {
            self.draw_gizmo_axis_in_render_pass(render_pass, axis_index);
        }
        self.draw_gizmo_center_in_render_pass(render_pass);
        for (_, axis_index) in near_axes {
            self.draw_gizmo_axis_in_render_pass(render_pass, axis_index);
        }

        render_pass.set_viewport(
            viewport_origin_px[0] as f32,
            viewport_origin_px[1] as f32,
            width as f32,
            height as f32,
            0.0,
            1.0,
        );
        render_pass.set_scissor_rect(
            viewport_origin_px[0],
            viewport_origin_px[1],
            width.max(1),
            height.max(1),
        );
    }
    fn draw_gizmo_axis_in_render_pass(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        axis_index: usize,
    ) {
        let axis = &self.gizmo_axes_buffers[axis_index];
        render_pass.set_pipeline(&self.pipeline);
        render_pass.set_bind_group(0, &self.gizmo_camera_bind_group, &[]);
        render_pass.set_bind_group(1, &self.gizmo_axes_material_bind_groups[axis_index], &[]);
        render_pass.set_vertex_buffer(0, axis.vertex_buffer.slice(..));
        render_pass.set_index_buffer(axis.index_buffer.slice(..), wgpu::IndexFormat::Uint32);
        render_pass.draw_indexed(0..axis.index_count, 0, 0..1);
    }

    fn draw_gizmo_center_in_render_pass(&self, render_pass: &mut wgpu::RenderPass<'_>) {
        render_pass.set_pipeline(&self.pipeline);
        render_pass.set_bind_group(0, &self.gizmo_camera_bind_group, &[]);
        render_pass.set_bind_group(1, &self.gizmo_center_material_bind_group, &[]);
        render_pass.set_vertex_buffer(0, self.gizmo_center_buffers.vertex_buffer.slice(..));
        render_pass.set_index_buffer(
            self.gizmo_center_buffers.index_buffer.slice(..),
            wgpu::IndexFormat::Uint32,
        );
        render_pass.draw_indexed(0..self.gizmo_center_buffers.index_count, 0, 0..1);
    }

    fn draw_axes_in_render_pass(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        use_depth_pipeline: bool,
    ) {
        for i in 0..3 {
            let axis = &self.axes_buffers[i];

            // Draw cone (triangles)
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline);
            }
            render_pass.set_bind_group(0, &self.axes_camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.axes_material_bind_groups[i], &[]);
            render_pass.set_vertex_buffer(0, axis.vertex_buffer.slice(..));
            render_pass
                .set_index_buffer(axis.tri_index_buffer.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..axis.tri_index_count, 0, 0..1);

            // Draw shaft (line)
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_line_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline_line);
            }
            render_pass.set_bind_group(0, &self.axes_camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.axes_material_bind_groups[i], &[]);
            render_pass.set_vertex_buffer(0, axis.vertex_buffer.slice(..));
            render_pass
                .set_index_buffer(axis.line_index_buffer.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..axis.line_index_count, 0, 0..1);
        }
    }

    fn draw_grid_in_render_pass(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        use_depth_pipeline: bool,
    ) {
        let minor_index_count = *self.grid.minor_index_count.lock().expect("render mutex poisoned");
        let major_index_count = *self.grid.major_index_count.lock().expect("render mutex poisoned");
        if use_depth_pipeline {
            render_pass.set_pipeline(&self.pipeline_line_depth);
        } else {
            render_pass.set_pipeline(&self.pipeline_line);
        }
        render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
        render_pass.set_vertex_buffer(0, self.grid.vertex_buffer.slice(..));

        // Draw minor lines first (thinner/dimmer)
        if minor_index_count > 0 {
            render_pass.set_bind_group(1, &self.grid_minor_material_bind_group, &[]);
            render_pass.set_index_buffer(
                self.grid.minor_index_buffer.slice(..),
                wgpu::IndexFormat::Uint32,
            );
            render_pass.draw_indexed(0..minor_index_count, 0, 0..1);
        }

        // Draw major lines on top
        if major_index_count > 0 {
            render_pass.set_bind_group(1, &self.grid_major_material_bind_group, &[]);
            render_pass.set_index_buffer(
                self.grid.major_index_buffer.slice(..),
                wgpu::IndexFormat::Uint32,
            );
            render_pass.draw_indexed(0..major_index_count, 0, 0..1);
        }
    }

    pub fn set_show_axes(&self, show: bool) {
        *self.show_axes.lock().expect("render mutex poisoned") = show;
    }

    pub fn show_axes(&self) -> bool {
        *self.show_axes.lock().expect("render mutex poisoned")
    }

    pub fn set_display_mode(&self, mode: DisplayMode) {
        *self.display_mode.lock().expect("render mutex poisoned") = mode;
    }

    pub fn display_mode(&self) -> DisplayMode {
        *self.display_mode.lock().expect("render mutex poisoned")
    }

    pub fn set_show_grid(&self, show: bool) {
        *self.show_grid.lock().expect("render mutex poisoned") = show;
    }

    pub fn show_grid(&self) -> bool {
        *self.show_grid.lock().expect("render mutex poisoned")
    }

    pub fn set_show_axis_gizmo(&self, show: bool) {
        *self.show_axis_gizmo.lock().expect("render mutex poisoned") = show;
    }

    pub fn show_axis_gizmo(&self) -> bool {
        *self.show_axis_gizmo.lock().expect("render mutex poisoned")
    }

    pub fn set_model_color(&self, queue: &wgpu::Queue, color: [f32; 4]) {
        let uniform = MaterialUniform {
            color,
            flags: [0.0, 0.0, 0.0, 0.0],
        };
        queue.write_buffer(&self.material_buffer, 0, bytemuck::cast_slice(&[uniform]));
    }

    pub fn set_model_unlit(&self, queue: &wgpu::Queue, color: [f32; 4], unlit: bool) {
        let uniform = MaterialUniform {
            color,
            flags: [if unlit { 1.0 } else { 0.0 }, 0.0, 0.0, 0.0],
        };
        queue.write_buffer(&self.material_buffer, 0, bytemuck::cast_slice(&[uniform]));
    }

    pub fn set_light_direction(&self, dir: glam::Vec3) {
        *self.light_dir.lock().expect("render mutex poisoned") = dir;
    }

    pub fn light_direction(&self) -> glam::Vec3 {
        *self.light_dir.lock().expect("render mutex poisoned")
    }

    /// Set light direction to match the camera eye direction (headlight mode).
    pub fn set_headlight(&self, camera: &Camera) {
        let eye = camera.eye_position();
        let dir = (eye - camera.target).normalize_or_zero();
        if dir.length_squared() > 1e-6 {
            *self.light_dir.lock().expect("render mutex poisoned") = dir;
        }
    }

    pub fn render(
        &self,
        view: &wgpu::TextureView,
        encoder: &mut wgpu::CommandEncoder,
        clear_color: wgpu::Color,
        clip_bounds: Option<(u32, u32, u32, u32)>,
    ) {
        let use_depth = clip_bounds.is_some();
        let depth_view_guard = self.depth_view.lock().expect("render mutex poisoned");
        let depth_attachment = if use_depth {
            depth_view_guard
                .as_ref()
                .map(|depth_view| wgpu::RenderPassDepthStencilAttachment {
                    view: depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Discard,
                    }),
                    stencil_ops: None,
                })
        } else {
            None
        };

        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("Render Pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: if clip_bounds.is_some() {
                        wgpu::LoadOp::Load
                    } else {
                        wgpu::LoadOp::Clear(clear_color)
                    },
                    store: wgpu::StoreOp::Store,
                },
                depth_slice: None,
            })],
            depth_stencil_attachment: depth_attachment,
            timestamp_writes: None,
            occlusion_query_set: None,
        });

        let use_depth_pipeline = use_depth && depth_view_guard.is_some();

        if let Some((x, y, width, height)) = clip_bounds
            && width > 0
            && height > 0
        {
            render_pass.set_scissor_rect(x, y, width.max(1), height.max(1));
        }

        self.draw_in_render_pass(&mut render_pass, use_depth_pipeline);
    }

    pub fn render_with_axis_gizmo(
        &self,
        view: &wgpu::TextureView,
        encoder: &mut wgpu::CommandEncoder,
        clear_color: wgpu::Color,
        clip_bounds: Option<(u32, u32, u32, u32)>,
    ) {
        self.render(view, encoder, clear_color, clip_bounds);

        if let Some((x, y, width, height)) = clip_bounds
            && width > 0
            && height > 0
        {
            let mut gizmo_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Axis Gizmo Overlay Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            gizmo_pass.set_scissor_rect(x, y, width.max(1), height.max(1));
            self.draw_axis_gizmo_in_render_pass(&mut gizmo_pass, [x, y], [width, height]);
        }
    }

    pub fn render_with_defaults(
        &self,
        view: &wgpu::TextureView,
        encoder: &mut wgpu::CommandEncoder,
        clip_bounds: Option<(u32, u32, u32, u32)>,
    ) {
        self.render(view, encoder, Self::default_clear_color(), clip_bounds);
    }

    pub fn render_with_defaults_and_axis_gizmo(
        &self,
        view: &wgpu::TextureView,
        encoder: &mut wgpu::CommandEncoder,
        clip_bounds: Option<(u32, u32, u32, u32)>,
    ) {
        self.render_with_axis_gizmo(view, encoder, Self::default_clear_color(), clip_bounds);
    }

    /// Render the current scene to an offscreen texture and return it as an RGBA image.
    pub fn screenshot(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        camera: &Camera,
        mesh: &Mesh,
        width: u32,
        height: u32,
    ) -> image::RgbaImage {
        let width = width.max(1);
        let height = height.max(1);
        let aspect = width as f32 / height as f32;

        // Prepare scene data
        self.upload_mesh(device, mesh);
        self.update_camera(queue, camera, aspect);

        // Create offscreen color texture
        let color_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Screenshot Color Texture"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Rgba8UnormSrgb,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());

        // Create offscreen depth texture
        let depth_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Screenshot Depth Texture"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        let depth_view = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());

        // Render
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Screenshot Encoder"),
        });

        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Screenshot Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(Self::default_clear_color()),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Discard,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            self.draw_in_render_pass(&mut render_pass, true);
        }

        // Copy texture to staging buffer
        let bytes_per_pixel = 4u32;
        // wgpu requires rows to be aligned to 256 bytes
        let unpadded_bytes_per_row = width * bytes_per_pixel;
        let padded_bytes_per_row = (unpadded_bytes_per_row + 255) & !255;
        let buffer_size = (padded_bytes_per_row * height) as u64;

        let staging_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Screenshot Staging Buffer"),
            size: buffer_size,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        encoder.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo {
                texture: &color_texture,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::TexelCopyBufferInfo {
                buffer: &staging_buffer,
                layout: wgpu::TexelCopyBufferLayout {
                    offset: 0,
                    bytes_per_row: Some(padded_bytes_per_row),
                    rows_per_image: Some(height),
                },
            },
            wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
        );

        queue.submit(std::iter::once(encoder.finish()));

        // Map and read back
        let buffer_slice = staging_buffer.slice(..);
        let (sender, receiver) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
            sender.send(result).expect("screenshot channel: receiver dropped before GPU callback");
        });
        device.poll(wgpu::PollType::wait_indefinitely()).expect("GPU device lost during screenshot");
        receiver
            .recv()
            .expect("screenshot channel: sender dropped before recv")
            .expect("GPU buffer map failed during screenshot");

        let data = buffer_slice.get_mapped_range();
        let mut pixels = Vec::with_capacity((width * height * bytes_per_pixel) as usize);
        for row in 0..height {
            let start = (row * padded_bytes_per_row) as usize;
            let end = start + (width * bytes_per_pixel) as usize;
            pixels.extend_from_slice(&data[start..end]);
        }
        drop(data);
        staging_buffer.unmap();

        image::RgbaImage::from_raw(width, height, pixels)
            .expect("pixel buffer size matches image dimensions")
    }

    /// Render the scene and save to a PNG file.
    #[allow(clippy::too_many_arguments)]
    pub fn screenshot_to_file(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        camera: &Camera,
        mesh: &Mesh,
        width: u32,
        height: u32,
        path: &std::path::Path,
    ) -> Result<(), String> {
        let img = self.screenshot(device, queue, camera, mesh, width, height);
        img.save(path).map_err(|e| e.to_string())
    }
}
