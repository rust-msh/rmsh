pub mod assembly;
pub use assembly::{Assembly, AssemblyMetadata, AssemblyNode, NodeContent, assembly_from_parts};

use glam::DVec3;
use rcad_kernel::{BRep, Edge, Face, Shell, Solid, Wire, WireEdge};
use rcad_modeling::{make_box_brep, make_cone_brep, make_cylinder_brep, make_sphere_brep, make_torus_brep};
use rcad_render::{
    Camera, DEFAULT_EDGE_PICK_RADIUS_PX, SelectionMode, SelectionState, cursor_point_on_plane,
};
use std::collections::HashSet;

const WORK_PLANE_ORIGIN: DVec3 = DVec3::ZERO;

/// The active construction plane for click-based geometry creation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkPlane {
    /// XY plane (normal = +Z).
    XY,
    /// XZ plane (normal = +Y).
    XZ,
    /// YZ plane (normal = +X).
    YZ,
}

impl WorkPlane {
    pub fn normal(self) -> DVec3 {
        match self {
            WorkPlane::XY => DVec3::Z,
            WorkPlane::XZ => DVec3::Y,
            WorkPlane::YZ => DVec3::X,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            WorkPlane::XY => "XY",
            WorkPlane::XZ => "XZ",
            WorkPlane::YZ => "YZ",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tool {
    SelectFace,
    SelectEdge,
    Box,
    Sphere,
    Cylinder,
    Cone,
    Torus,
}

#[derive(Debug, Clone, Copy)]
pub enum CommandState {
    Idle,
    BoxBase {
        first: DVec3,
        current: DVec3,
    },
    BoxHeight {
        first: DVec3,
        second: DVec3,
        anchor_screen_y: f32,
        current_screen_y: f32,
    },
    SphereRadius {
        center: DVec3,
        current: DVec3,
    },
    /// Cylinder: step 1 = pick base center, step 2 = pick radius point, step 3 = pick height
    CylinderRadius {
        center: DVec3,
        current: DVec3,
    },
    CylinderHeight {
        center: DVec3,
        radius_pt: DVec3,
        anchor_screen_y: f32,
        current_screen_y: f32,
    },
    /// Cone: step 1 = pick base center, step 2 = pick radius point, step 3 = pick height
    ConeRadius {
        center: DVec3,
        current: DVec3,
    },
    ConeHeight {
        center: DVec3,
        radius_pt: DVec3,
        anchor_screen_y: f32,
        current_screen_y: f32,
    },
    /// Torus: step 1 = pick center, step 2 = pick major radius, step 3 = pick minor radius
    TorusMajor {
        center: DVec3,
        current: DVec3,
    },
    TorusMinor {
        center: DVec3,
        major_pt: DVec3,
        current: DVec3,
    },
}

#[derive(Debug, Clone)]
pub struct CreationController {
    active_tool: Tool,
    command_state: CommandState,
    work_plane: WorkPlane,
}

impl Default for CreationController {
    fn default() -> Self {
        Self {
            active_tool: Tool::SelectFace,
            command_state: CommandState::Idle,
            work_plane: WorkPlane::XY,
        }
    }
}

impl CreationController {
    pub fn active_tool(&self) -> Tool {
        self.active_tool
    }

    pub fn command_state(&self) -> CommandState {
        self.command_state
    }

    pub fn work_plane(&self) -> WorkPlane {
        self.work_plane
    }

    pub fn set_work_plane(&mut self, plane: WorkPlane) {
        self.work_plane = plane;
        // Reset the active command when the work plane changes.
        self.command_state = CommandState::Idle;
    }

    pub fn is_command_active(&self) -> bool {
        !matches!(self.command_state, CommandState::Idle)
    }

    pub fn set_tool(&mut self, tool: Tool, selection: &mut SelectionState) {
        self.active_tool = tool;
        self.command_state = CommandState::Idle;
        match tool {
            Tool::SelectFace => selection.set_mode(SelectionMode::Face),
            Tool::SelectEdge => selection.set_mode(SelectionMode::Edge),
            Tool::Box
            | Tool::Sphere
            | Tool::Cylinder
            | Tool::Cone
            | Tool::Torus => selection.clear_hover(),
        }
    }

    pub fn set_selection_mode(&mut self, mode: SelectionMode, selection: &mut SelectionState) {
        match mode {
            SelectionMode::Face => self.set_tool(Tool::SelectFace, selection),
            SelectionMode::Edge => self.set_tool(Tool::SelectEdge, selection),
        }
    }

    pub fn set_additive_select(&self, selection: &mut SelectionState, enabled: bool) {
        selection.additive_select = enabled;
    }

    pub fn clear_hover_if_selection_tool(&self, selection: &mut SelectionState) {
        if matches!(self.active_tool, Tool::SelectFace | Tool::SelectEdge) {
            selection.clear_hover();
        }
    }

    pub fn clear_selection(&self, selection: &mut SelectionState) {
        selection.selected_faces.clear();
        selection.selected_edges.clear();
    }

    pub fn select_face_boundary_edges(&self, brep: &BRep, selection: &mut SelectionState) {
        let mut edge_set: HashSet<usize> = HashSet::new();
        for face_idx in &selection.selected_faces {
            if let Some(face) = face_by_index(brep, *face_idx) {
                for edge in &face.outer_wire.edges {
                    edge_set.insert(edge.idx);
                }
                for wire in &face.inner_wires {
                    for edge in &wire.edges {
                        edge_set.insert(edge.idx);
                    }
                }
            }
        }
        selection.selected_edges = edge_set.into_iter().collect();
        selection.selected_edges.sort_unstable();
    }

    pub fn select_edge_incident_faces(&self, brep: &BRep, selection: &mut SelectionState) {
        if selection.selected_edges.is_empty() {
            return;
        }
        let selected_edges: HashSet<usize> = selection.selected_edges.iter().copied().collect();
        let mut faces: Vec<usize> = Vec::new();
        let mut face_idx = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let hit_outer = face
                        .outer_wire
                        .edges
                        .iter()
                        .any(|e| selected_edges.contains(&e.idx));
                    let hit_inner = face
                        .inner_wires
                        .iter()
                        .any(|w| w.edges.iter().any(|e| selected_edges.contains(&e.idx)));
                    if hit_outer || hit_inner {
                        faces.push(face_idx);
                    }
                    face_idx += 1;
                }
            }
        }
        selection.selected_faces = faces;
    }

    pub fn grow_selected_faces(&self, brep: &BRep, selection: &mut SelectionState) {
        if selection.selected_faces.is_empty() {
            return;
        }
        let mut edge_to_faces: Vec<Vec<usize>> = vec![Vec::new(); brep.edges.len()];
        let mut face_idx = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for edge in &face.outer_wire.edges {
                        if edge.idx < edge_to_faces.len() {
                            edge_to_faces[edge.idx].push(face_idx);
                        }
                    }
                    for wire in &face.inner_wires {
                        for edge in &wire.edges {
                            if edge.idx < edge_to_faces.len() {
                                edge_to_faces[edge.idx].push(face_idx);
                            }
                        }
                    }
                    face_idx += 1;
                }
            }
        }

        let mut grown: HashSet<usize> = selection.selected_faces.iter().copied().collect();
        for &f in &selection.selected_faces {
            if let Some(face) = face_by_index(brep, f) {
                for edge in &face.outer_wire.edges {
                    if let Some(adj) = edge_to_faces.get(edge.idx) {
                        for &af in adj {
                            grown.insert(af);
                        }
                    }
                }
                for wire in &face.inner_wires {
                    for edge in &wire.edges {
                        if let Some(adj) = edge_to_faces.get(edge.idx) {
                            for &af in adj {
                                grown.insert(af);
                            }
                        }
                    }
                }
            }
        }

        selection.selected_faces = grown.into_iter().collect();
        selection.selected_faces.sort_unstable();
    }

    pub fn grow_selected_edges(&self, brep: &BRep, selection: &mut SelectionState) {
        if selection.selected_edges.is_empty() {
            return;
        }

        let mut seed_vertices: HashSet<usize> = HashSet::new();
        for &edge_idx in &selection.selected_edges {
            if let Some(edge) = brep.edges.get(edge_idx) {
                seed_vertices.insert(edge.start);
                seed_vertices.insert(edge.end);
            }
        }

        let mut grown: HashSet<usize> = selection.selected_edges.iter().copied().collect();
        for (idx, edge) in brep.edges.iter().enumerate() {
            if seed_vertices.contains(&edge.start) || seed_vertices.contains(&edge.end) {
                grown.insert(idx);
            }
        }

        selection.selected_edges = grown.into_iter().collect();
        selection.selected_edges.sort_unstable();
    }

    pub fn tool_name(&self) -> &'static str {
        match self.active_tool {
            Tool::SelectFace => "Select Face",
            Tool::SelectEdge => "Select Edge",
            Tool::Box => "Box",
            Tool::Sphere => "Sphere",
            Tool::Cylinder => "Cylinder",
            Tool::Cone => "Cone",
            Tool::Torus => "Torus",
        }
    }

    pub fn command_hint(&self) -> &'static str {
        match self.command_state {
            CommandState::Idle => match self.active_tool {
                Tool::SelectFace => "Face selection mode",
                Tool::SelectEdge => "Edge selection mode",
                Tool::Box => "Box: click first corner",
                Tool::Sphere => "Sphere: click center",
                Tool::Cylinder => "Cylinder: click base center",
                Tool::Cone => "Cone: click base center",
                Tool::Torus => "Torus: click center",
            },
            CommandState::BoxBase { .. } => "Box: click opposite base corner",
            CommandState::BoxHeight { .. } => {
                "Box: move mouse for height, click or Enter to finish"
            }
            CommandState::SphereRadius { .. } => "Sphere: click radius point or Enter to finish",
            CommandState::CylinderRadius { .. } => "Cylinder: click radius point",
            CommandState::CylinderHeight { .. } => {
                "Cylinder: move mouse for height, click or Enter to finish"
            }
            CommandState::ConeRadius { .. } => "Cone: click base radius point",
            CommandState::ConeHeight { .. } => {
                "Cone: move mouse for height, click or Enter to finish"
            }
            CommandState::TorusMajor { .. } => "Torus: click major radius point",
            CommandState::TorusMinor { .. } => "Torus: click minor radius point or Enter to finish",
        }
    }

    pub fn handle_primary_click(
        &mut self,
        brep: &BRep,
        camera: &Camera,
        selection: &mut SelectionState,
        cursor: [f32; 2],
        viewport: [f32; 2],
    ) -> Option<BRep> {
        let aspect = viewport[0] / viewport[1].max(1.0);
        match self.active_tool {
            Tool::SelectFace | Tool::SelectEdge => {
                selection.click_at(
                    brep,
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    DEFAULT_EDGE_PICK_RADIUS_PX,
                );
                None
            }
            Tool::Box => {
                let point = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                )?;
                self.command_state = match self.command_state {
                    CommandState::Idle => CommandState::BoxBase {
                        first: point,
                        current: point,
                    },
                    CommandState::BoxBase { first, .. } => {
                        if (point - first).length_squared() < 1e-8 {
                            CommandState::BoxBase {
                                first,
                                current: point,
                            }
                        } else {
                            CommandState::BoxHeight {
                                first,
                                second: point,
                                anchor_screen_y: cursor[1],
                                current_screen_y: cursor[1],
                            }
                        }
                    }
                    CommandState::BoxHeight { .. } => {
                        let finished = self.preview_brep(camera.distance);
                        self.command_state = CommandState::Idle;
                        return finished;
                    }
                    other => other,
                };
                None
            }
            Tool::Sphere => {
                let point = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                )?;
                self.command_state = match self.command_state {
                    CommandState::Idle => CommandState::SphereRadius {
                        center: point,
                        current: point,
                    },
                    CommandState::SphereRadius { .. } => {
                        let finished = self.preview_brep(camera.distance);
                        self.command_state = CommandState::Idle;
                        return finished;
                    }
                    other => other,
                };
                None
            }
            Tool::Cylinder => {
                let point = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                )?;
                self.command_state = match self.command_state {
                    CommandState::Idle => CommandState::CylinderRadius {
                        center: point,
                        current: point,
                    },
                    CommandState::CylinderRadius { center, .. } => {
                        if center.distance_squared(point) < 1e-8 {
                            CommandState::CylinderRadius { center, current: point }
                        } else {
                            CommandState::CylinderHeight {
                                center,
                                radius_pt: point,
                                anchor_screen_y: cursor[1],
                                current_screen_y: cursor[1],
                            }
                        }
                    }
                    CommandState::CylinderHeight { .. } => {
                        let finished = self.preview_brep(camera.distance);
                        self.command_state = CommandState::Idle;
                        return finished;
                    }
                    other => other,
                };
                None
            }
            Tool::Cone => {
                let point = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                )?;
                self.command_state = match self.command_state {
                    CommandState::Idle => CommandState::ConeRadius {
                        center: point,
                        current: point,
                    },
                    CommandState::ConeRadius { center, .. } => {
                        if center.distance_squared(point) < 1e-8 {
                            CommandState::ConeRadius { center, current: point }
                        } else {
                            CommandState::ConeHeight {
                                center,
                                radius_pt: point,
                                anchor_screen_y: cursor[1],
                                current_screen_y: cursor[1],
                            }
                        }
                    }
                    CommandState::ConeHeight { .. } => {
                        let finished = self.preview_brep(camera.distance);
                        self.command_state = CommandState::Idle;
                        return finished;
                    }
                    other => other,
                };
                None
            }
            Tool::Torus => {
                let point = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                )?;
                self.command_state = match self.command_state {
                    CommandState::Idle => CommandState::TorusMajor {
                        center: point,
                        current: point,
                    },
                    CommandState::TorusMajor { center, current } => {
                        let major_r = center.distance(current);
                        if major_r < 1e-6 {
                            CommandState::TorusMajor { center, current: point }
                        } else {
                            CommandState::TorusMinor {
                                center,
                                major_pt: current,
                                current: point,
                            }
                        }
                    }
                    CommandState::TorusMinor { .. } => {
                        let finished = self.preview_brep(camera.distance);
                        self.command_state = CommandState::Idle;
                        return finished;
                    }
                    other => other,
                };
                None
            }
        }
    }

    pub fn handle_pointer_move(
        &mut self,
        brep: &BRep,
        camera: &Camera,
        selection: &mut SelectionState,
        cursor: [f32; 2],
        viewport: [f32; 2],
    ) {
        let aspect = viewport[0] / viewport[1].max(1.0);
        match self.command_state {
            CommandState::BoxBase { first, .. } => {
                if let Some(point) = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                ) {
                    self.command_state = CommandState::BoxBase {
                        first,
                        current: point,
                    };
                }
            }
            CommandState::BoxHeight {
                first,
                second,
                anchor_screen_y,
                ..
            } => {
                self.command_state = CommandState::BoxHeight {
                    first,
                    second,
                    anchor_screen_y,
                    current_screen_y: cursor[1],
                };
            }
            CommandState::SphereRadius { center, .. } => {
                if let Some(point) = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                ) {
                    self.command_state = CommandState::SphereRadius {
                        center,
                        current: point,
                    };
                }
            }
            CommandState::CylinderRadius { center, .. } => {
                if let Some(point) = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                ) {
                    self.command_state = CommandState::CylinderRadius {
                        center,
                        current: point,
                    };
                }
            }
            CommandState::CylinderHeight {
                center,
                radius_pt,
                anchor_screen_y,
                ..
            } => {
                self.command_state = CommandState::CylinderHeight {
                    center,
                    radius_pt,
                    anchor_screen_y,
                    current_screen_y: cursor[1],
                };
            }
            CommandState::ConeRadius { center, .. } => {
                if let Some(point) = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                ) {
                    self.command_state = CommandState::ConeRadius {
                        center,
                        current: point,
                    };
                }
            }
            CommandState::ConeHeight {
                center,
                radius_pt,
                anchor_screen_y,
                ..
            } => {
                self.command_state = CommandState::ConeHeight {
                    center,
                    radius_pt,
                    anchor_screen_y,
                    current_screen_y: cursor[1],
                };
            }
            CommandState::TorusMajor { center, .. } => {
                if let Some(point) = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                ) {
                    self.command_state = CommandState::TorusMajor {
                        center,
                        current: point,
                    };
                }
            }
            CommandState::TorusMinor { center, major_pt, .. } => {
                if let Some(point) = cursor_point_on_plane(
                    camera,
                    aspect,
                    viewport,
                    cursor,
                    WORK_PLANE_ORIGIN,
                    self.work_plane.normal(),
                ) {
                    self.command_state = CommandState::TorusMinor {
                        center,
                        major_pt,
                        current: point,
                    };
                }
            }
            CommandState::Idle => {
                if matches!(self.active_tool, Tool::SelectFace | Tool::SelectEdge) {
                    selection.hover_at(
                        brep,
                        camera,
                        aspect,
                        viewport,
                        cursor,
                        DEFAULT_EDGE_PICK_RADIUS_PX,
                    );
                }
            }
        }
    }

    pub fn cancel_active_command(&mut self) {
        self.command_state = CommandState::Idle;
    }

    pub fn undo_last_step(&mut self) {
        self.command_state = match self.command_state {
            CommandState::Idle => CommandState::Idle,
            CommandState::BoxBase { .. } => CommandState::Idle,
            CommandState::BoxHeight { first, second, .. } => CommandState::BoxBase {
                first,
                current: second,
            },
            CommandState::SphereRadius { .. } => CommandState::Idle,
            CommandState::CylinderRadius { .. } => CommandState::Idle,
            CommandState::CylinderHeight { center, radius_pt, .. } => CommandState::CylinderRadius {
                center,
                current: radius_pt,
            },
            CommandState::ConeRadius { .. } => CommandState::Idle,
            CommandState::ConeHeight { center, radius_pt, .. } => CommandState::ConeRadius {
                center,
                current: radius_pt,
            },
            CommandState::TorusMajor { .. } => CommandState::Idle,
            CommandState::TorusMinor { center, major_pt, .. } => CommandState::TorusMajor {
                center,
                current: major_pt,
            },
        };
    }

    pub fn confirm_active_command(&mut self, camera: &Camera) -> Option<BRep> {
        let finished = self.preview_brep(camera.distance);
        if finished.is_some() {
            self.command_state = CommandState::Idle;
        }
        finished
    }

    pub fn preview_brep(&self, camera_distance: f32) -> Option<BRep> {
        match self.command_state {
            CommandState::Idle => None,
            CommandState::BoxBase { first, current } => build_box_from_points(first, current, 0.02),
            CommandState::BoxHeight {
                first,
                second,
                anchor_screen_y,
                current_screen_y,
            } => build_box_from_points(
                first,
                second,
                box_height_from_screen(anchor_screen_y, current_screen_y, camera_distance),
            ),
            CommandState::SphereRadius { center, current } => {
                build_sphere_from_points(center, current)
            }
            CommandState::CylinderRadius { center, current } => {
                let radius = center.distance(current).max(0.01);
                build_cylinder_from_points(center, radius, 0.1)
            }
            CommandState::CylinderHeight {
                center,
                radius_pt,
                anchor_screen_y,
                current_screen_y,
            } => {
                let radius = center.distance(radius_pt).max(0.01);
                let height = box_height_from_screen(anchor_screen_y, current_screen_y, camera_distance);
                build_cylinder_from_points(center, radius, height)
            }
            CommandState::ConeRadius { center, current } => {
                let radius = center.distance(current).max(0.01);
                build_cone_from_points(center, radius, 0.1)
            }
            CommandState::ConeHeight {
                center,
                radius_pt,
                anchor_screen_y,
                current_screen_y,
            } => {
                let radius = center.distance(radius_pt).max(0.01);
                let height = box_height_from_screen(anchor_screen_y, current_screen_y, camera_distance);
                build_cone_from_points(center, radius, height)
            }
            CommandState::TorusMajor { center, current } => {
                let major_r = center.distance(current).max(0.05);
                build_torus_from_points(center, major_r, major_r * 0.2)
            }
            CommandState::TorusMinor { center, major_pt, current } => {
                let major_r = center.distance(major_pt).max(0.05);
                // Minor radius is the distance from current to the major circle in the plane.
                // Approximate: distance from current to center, clamped relative to major.
                let minor_r = center.distance(current).max(0.01).min(major_r * 0.9);
                build_torus_from_points(center, major_r, minor_r)
            }
        }
    }
}

fn face_by_index(brep: &BRep, face_index: usize) -> Option<&Face> {
    let mut current = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if current == face_index {
                    return Some(face);
                }
                current += 1;
            }
        }
    }
    None
}

fn box_height_from_screen(
    anchor_screen_y: f32,
    current_screen_y: f32,
    camera_distance: f32,
) -> f64 {
    let raw = (anchor_screen_y - current_screen_y) * camera_distance * 0.01;
    if raw.abs() < 0.1 { 1.0 } else { raw as f64 }
}

fn build_box_from_points(first: DVec3, second: DVec3, depth_signed: f64) -> Option<BRep> {
    let min_x = first.x.min(second.x);
    let min_y = first.y.min(second.y);
    let width = (first.x - second.x).abs();
    let height = (first.y - second.y).abs();
    let depth = depth_signed.abs().max(0.02);
    if width < 1e-6 || height < 1e-6 {
        return None;
    }

    let origin_z = if depth_signed < 0.0 {
        depth_signed
    } else {
        0.0
    };
    make_box_brep(
        DVec3::new(min_x, min_y, origin_z),
        DVec3::X,
        DVec3::Y,
        width,
        height,
        depth,
    )
    .ok()
}

fn build_sphere_from_points(center: DVec3, current: DVec3) -> Option<BRep> {
    let radius = center.distance(current);
    if radius < 1e-6 {
        return None;
    }
    make_sphere_brep(center, radius).ok()
}

fn build_cylinder_from_points(center: DVec3, radius: f64, height: f64) -> Option<BRep> {
    let height = height.abs().max(0.01);
    make_cylinder_brep(center, DVec3::Z, DVec3::X, radius, height).ok()
}

fn build_cone_from_points(center: DVec3, base_radius: f64, height: f64) -> Option<BRep> {
    let height = height.abs().max(0.01);
    make_cone_brep(center, DVec3::Z, DVec3::X, base_radius, height).ok()
}

fn build_torus_from_points(center: DVec3, major_radius: f64, minor_radius: f64) -> Option<BRep> {
    make_torus_brep(center, DVec3::Z, DVec3::X, major_radius, minor_radius).ok()
}
pub fn append_brep(dst: &mut BRep, src: BRep) {
    let vertex_offset = dst.vertices.len();
    let edge_offset = dst.edges.len();
    let curve_offset = dst.geom.curves.len();
    let surface_offset = dst.geom.surfaces.len();
    let src_face_surface = src.geom.face_surface.clone();

    dst.vertices.extend(src.vertices.iter().cloned());
    dst.edges.extend(src.edges.iter().map(|edge| Edge {
        start: edge.start + vertex_offset,
        end: edge.end + vertex_offset,
    }));

    dst.geom.curves.extend(src.geom.curves.iter().cloned());
    dst.geom.surfaces.extend(src.geom.surfaces.iter().cloned());
    dst.geom.edge_curve.extend(
        src.geom
            .edge_curve
            .iter()
            .map(|curve| curve.map(|idx| idx + curve_offset)),
    );
    // edge_curve_range has no index offset — parameters are curve-local values
    dst.geom
        .edge_curve_range
        .extend(src.geom.edge_curve_range.iter().cloned());
    dst.geom
        .edge_degenerated
        .extend(src.geom.edge_degenerated.iter().cloned());

    let mut face_counter = 0usize;
    for solid in src.solids {
        let mut new_shells = Vec::with_capacity(solid.shells.len());
        for shell in solid.shells {
            let mut new_faces = Vec::with_capacity(shell.faces.len());
            for face in shell.faces {
                let surface = src_face_surface
                    .get(face_counter)
                    .copied()
                    .flatten()
                    .map(|idx| idx + surface_offset);
                dst.geom.face_surface.push(surface);
                face_counter += 1;

                new_faces.push(Face {
                    outer_wire: Wire {
                        edges: face
                            .outer_wire
                            .edges
                            .into_iter()
                            .map(|we| WireEdge {
                                idx: we.idx + edge_offset,
                                forward: we.forward,
                            })
                            .collect(),
                    },
                    inner_wires: face
                        .inner_wires
                        .into_iter()
                        .map(|wire| Wire {
                            edges: wire
                                .edges
                                .into_iter()
                                .map(|we| WireEdge {
                                    idx: we.idx + edge_offset,
                                    forward: we.forward,
                                })
                                .collect(),
                        })
                        .collect(),
                    normal: face.normal,
                    triangles: face
                        .triangles
                        .into_iter()
                        .map(|tri| {
                            [
                                tri[0] + vertex_offset,
                                tri[1] + vertex_offset,
                                tri[2] + vertex_offset,
                            ]
                        })
                        .collect(),
                    mesh_dirty: face.mesh_dirty,
                });
            }
            new_shells.push(Shell { faces: new_faces });
        }
        dst.solids.push(Solid { shells: new_shells });
    }
}
