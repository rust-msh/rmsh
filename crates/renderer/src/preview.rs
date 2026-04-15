use glam::{Mat4, Vec4};

use crate::{OrbitCamera, RenderConfig};

#[derive(Debug, Clone)]
pub struct PreviewTriangle {
    pub depth: f32,
    pub points: [[f32; 2]; 3],
    pub color_rgb: [u8; 3],
}

#[derive(Debug, Clone)]
pub struct PreviewLine {
    pub start: [f32; 2],
    pub end: [f32; 2],
}

#[derive(Debug, Clone)]
pub struct PreviewOverlayLine {
    pub start: [f32; 2],
    pub end: [f32; 2],
    pub rgba: [u8; 4],
    pub width: f32,
}

#[derive(Debug, Clone)]
pub struct PreviewOverlayText {
    pub position: [f32; 2],
    pub text: String,
    pub rgba: [u8; 4],
    pub font_px: u32,
}

#[derive(Debug, Clone)]
pub struct PreviewFrame {
    pub triangles: Vec<PreviewTriangle>,
    pub lines: Vec<PreviewLine>,
    pub edge_rgb: [u8; 3],
    pub overlay_lines: Vec<PreviewOverlayLine>,
    pub overlay_texts: Vec<PreviewOverlayText>,
}

fn push_background_grid(out: &mut Vec<PreviewOverlayLine>, width: f32, height: f32) {
    let minor_step = 32.0f32;
    let major_every = 4i32;

    let mut i = 0i32;
    let mut x = 0.0f32;
    while x <= width {
        let major = i % major_every == 0;
        out.push(PreviewOverlayLine {
            start: [x, 0.0],
            end: [x, height],
            rgba: if major {
                [180, 195, 220, 28]
            } else {
                [150, 170, 200, 16]
            },
            width: if major { 1.0 } else { 0.8 },
        });
        x += minor_step;
        i += 1;
    }

    i = 0;
    let mut y = 0.0f32;
    while y <= height {
        let major = i % major_every == 0;
        out.push(PreviewOverlayLine {
            start: [0.0, y],
            end: [width, y],
            rgba: if major {
                [180, 195, 220, 28]
            } else {
                [150, 170, 200, 16]
            },
            width: if major { 1.0 } else { 0.8 },
        });
        y += minor_step;
        i += 1;
    }
}

fn push_axes_gizmo(
    lines: &mut Vec<PreviewOverlayLine>,
    texts: &mut Vec<PreviewOverlayText>,
    camera: &OrbitCamera,
    vp: Mat4,
    width: f32,
    height: f32,
) {
    let ox = 18.0f32;
    let oy = height - 22.0;

    let origin = [camera.target.x, camera.target.y, camera.target.z];
    let axis_len_world = camera.distance.max(1e-3) * 0.15;

    let fallback_dirs = [[1.0f32, 0.0f32], [0.0f32, -1.0f32], [-0.75f32, 0.6f32]];
    let axes = [
        ("X", [1.0f32, 0.0f32, 0.0f32], [255, 90, 90, 215], [255, 110, 110, 230], 1.6f32),
        ("Y", [0.0f32, 1.0f32, 0.0f32], [90, 240, 120, 215], [130, 250, 150, 230], 1.6f32),
        ("Z", [0.0f32, 0.0f32, 1.0f32], [120, 170, 255, 215], [145, 190, 255, 230], 1.4f32),
    ];

    let origin_proj = project(vp, origin, width, height);
    let gizmo_len = 34.0f32;

    for (idx, (label, axis, line_rgba, text_rgba, line_width)) in axes.iter().enumerate() {
        let end_world = [
            origin[0] + axis[0] * axis_len_world,
            origin[1] + axis[1] * axis_len_world,
            origin[2] + axis[2] * axis_len_world,
        ];

        let mut dir = fallback_dirs[idx];
        if let (Some(o), Some(e)) = (origin_proj, project(vp, end_world, width, height)) {
            let dx = e.0 - o.0;
            let dy = e.1 - o.1;
            let l = (dx * dx + dy * dy).sqrt();
            if l > 1e-4 {
                dir = [dx / l, dy / l];
            }
        }

        let end = [ox + dir[0] * gizmo_len, oy + dir[1] * gizmo_len];
        lines.push(PreviewOverlayLine {
            start: [ox, oy],
            end,
            rgba: *line_rgba,
            width: *line_width,
        });
        texts.push(PreviewOverlayText {
            position: [end[0] + 4.0, end[1] + 4.0],
            text: (*label).to_string(),
            rgba: *text_rgba,
            font_px: 11,
        });
    }
}

fn nice_step(value: f32) -> f32 {
    if !value.is_finite() || value <= 0.0 {
        return 1.0;
    }
    let exp = value.log10().floor();
    let base = 10f32.powf(exp);
    let n = value / base;
    let m = if n < 1.5 {
        1.0
    } else if n < 3.5 {
        2.0
    } else if n < 7.5 {
        5.0
    } else {
        10.0
    };
    m * base
}

fn format_step(step: f32) -> String {
    if step <= 0.0 || !step.is_finite() {
        return "scale".to_string();
    }
    if step >= 1.0 {
        format!("{:.0}", step)
    } else if step >= 0.1 {
        format!("{:.2}", step)
    } else {
        format!("{:.3}", step)
    }
}

fn push_scale_ruler(
    lines: &mut Vec<PreviewOverlayLine>,
    texts: &mut Vec<PreviewOverlayText>,
    camera: &OrbitCamera,
    vp: Mat4,
    width: f32,
    height: f32,
) {
    let margin = 16.0f32;

    let origin = [camera.target.x, camera.target.y, camera.target.z];
    let unit_x = [origin[0] + 1.0, origin[1], origin[2]];
    let px_per_world = if let (Some(o), Some(x1)) = (
        project(vp, origin, width, height),
        project(vp, unit_x, width, height),
    ) {
        let dx = x1.0 - o.0;
        let dy = x1.1 - o.1;
        (dx * dx + dy * dy).sqrt().max(1e-3)
    } else {
        (36.0 / camera.distance.max(0.2)).max(1e-3)
    };

    let target_px = (width * 0.18).clamp(72.0, 132.0);
    let world_step = nice_step(target_px / px_per_world);
    let ruler_len = (world_step * px_per_world).clamp(56.0, width * 0.35);
    let y = height - margin;
    let x0 = width - margin - ruler_len;
    let x1 = width - margin;

    lines.push(PreviewOverlayLine {
        start: [x0, y],
        end: [x1, y],
        rgba: [245, 245, 245, 210],
        width: 2.0,
    });

    for i in 0..=4 {
        let t = i as f32 / 4.0;
        let x = x0 + (x1 - x0) * t;
        let tick_h = if i % 2 == 0 { 8.0 } else { 5.0 };
        lines.push(PreviewOverlayLine {
            start: [x, y - tick_h],
            end: [x, y + 1.0],
            rgba: [245, 245, 245, 210],
            width: 1.2,
        });
    }

    texts.push(PreviewOverlayText {
        position: [x0, y - 10.0],
        text: format!("{} u", format_step(world_step)),
        rgba: [230, 230, 230, 220],
        font_px: 10,
    });
}

fn to_rgb_u8(color: [f32; 3]) -> [u8; 3] {
    [
        (color[0].clamp(0.0, 1.0) * 255.0) as u8,
        (color[1].clamp(0.0, 1.0) * 255.0) as u8,
        (color[2].clamp(0.0, 1.0) * 255.0) as u8,
    ]
}

fn project(vp: Mat4, p: [f32; 3], width: f32, height: f32) -> Option<(f32, f32, f32)> {
    let clip = vp * Vec4::new(p[0], p[1], p[2], 1.0);
    if clip.w.abs() < 1e-6 {
        return None;
    }

    let ndc = clip / clip.w;
    if ndc.z < -1.5 || ndc.z > 1.5 {
        return None;
    }

    let x = (ndc.x * 0.5 + 0.5) * width;
    let y = (1.0 - (ndc.y * 0.5 + 0.5)) * height;
    Some((x, y, ndc.z))
}

pub fn background_rgb(cfg: &RenderConfig) -> ([u8; 3], [u8; 3]) {
    (to_rgb_u8(cfg.bg_color_top), to_rgb_u8(cfg.bg_color_bottom))
}

pub fn build_preview_frame(
    camera: &OrbitCamera,
    cfg: &RenderConfig,
    surface_positions: &[[f32; 3]],
    surface_colors: &[[f32; 3]],
    surface_indices: &[u32],
    wire_positions: &[[f32; 3]],
    wire_indices: &[u32],
    node_count: usize,
    elem_count: usize,
    vol_count: usize,
    surf_count: usize,
    edge_count: usize,
    status_text: Option<&str>,
    width: u32,
    height: u32,
) -> PreviewFrame {
    let width = width.max(1) as f32;
    let height = height.max(1) as f32;
    let aspect = (width / height).max(0.01);
    let vp = Mat4::from_cols_array_2d(&camera.build_view_projection_matrix(aspect));

    let mut triangles = Vec::new();
    if cfg.show_faces {
        for tri in surface_indices.chunks_exact(3) {
            let p0 = surface_positions.get(tri[0] as usize).copied();
            let p1 = surface_positions.get(tri[1] as usize).copied();
            let p2 = surface_positions.get(tri[2] as usize).copied();
            if let (Some(p0), Some(p1), Some(p2)) = (p0, p1, p2) {
                if let (Some(a), Some(b), Some(c)) = (
                    project(vp, p0, width, height),
                    project(vp, p1, width, height),
                    project(vp, p2, width, height),
                ) {
                    let depth = (a.2 + b.2 + c.2) / 3.0;
                    let c0 = surface_colors.get(tri[0] as usize).copied().unwrap_or(cfg.face_color);
                    let c1 = surface_colors.get(tri[1] as usize).copied().unwrap_or(cfg.face_color);
                    let c2 = surface_colors.get(tri[2] as usize).copied().unwrap_or(cfg.face_color);
                    let avg = [
                        (c0[0] + c1[0] + c2[0]) / 3.0,
                        (c0[1] + c1[1] + c2[1]) / 3.0,
                        (c0[2] + c1[2] + c2[2]) / 3.0,
                    ];

                    triangles.push(PreviewTriangle {
                        depth,
                        points: [[a.0, a.1], [b.0, b.1], [c.0, c.1]],
                        color_rgb: to_rgb_u8(avg),
                    });
                }
            }
        }
        triangles.sort_by(|a, b| {
            a.depth
                .partial_cmp(&b.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    let mut lines = Vec::new();
    if cfg.show_edges {
        for seg in wire_indices.chunks_exact(2) {
            let p0 = wire_positions.get(seg[0] as usize).copied();
            let p1 = wire_positions.get(seg[1] as usize).copied();
            if let (Some(p0), Some(p1)) = (p0, p1) {
                if let (Some(a), Some(b)) = (
                    project(vp, p0, width, height),
                    project(vp, p1, width, height),
                ) {
                    lines.push(PreviewLine {
                        start: [a.0, a.1],
                        end: [b.0, b.1],
                    });
                }
            }
        }
    }

    let mut overlay_lines = Vec::new();
    let mut overlay_texts = Vec::new();

    if cfg.show_scale_ruler {
        push_background_grid(&mut overlay_lines, width, height);
        push_scale_ruler(&mut overlay_lines, &mut overlay_texts, camera, vp, width, height);
    }
    if cfg.show_axes {
        push_axes_gizmo(&mut overlay_lines, &mut overlay_texts, camera, vp, width, height);
    }

    if node_count > 0 || elem_count > 0 || !surface_indices.is_empty() || !wire_indices.is_empty() {
        overlay_texts.push(PreviewOverlayText {
            position: [12.0, 20.0],
            text: format!(
                "nodes={}  elements={} (V:{} S:{} E:{})  tris={}  lines={}  (drag rotate, shift+drag pan, wheel zoom)",
                node_count,
                elem_count,
                vol_count,
                surf_count,
                edge_count,
                surface_indices.len() / 3,
                wire_indices.len() / 2
            ),
            rgba: [255, 255, 255, 217],
            font_px: 12,
        });
    }

    if let Some(status) = status_text {
        overlay_texts.push(PreviewOverlayText {
            position: [12.0, height - 12.0],
            text: status.to_string(),
            rgba: [255, 255, 255, 200],
            font_px: 12,
        });
    }

    PreviewFrame {
        triangles,
        lines,
        edge_rgb: to_rgb_u8(cfg.edge_color),
        overlay_lines,
        overlay_texts,
    }
}
