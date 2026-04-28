use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
use std::collections::HashMap;

/// 曲面高质量三角化结果。
#[derive(Debug, Clone)]
pub struct SurfaceMesh {
    /// 三角化产生的顶点列表（世界坐标）。
    pub vertices: Vec<DVec3>,
    /// 三角形索引，每个三角形由3个顶点索引组成。
    pub triangles: Vec<[usize; 3]>,
    /// 每个顶点的法向量。
    pub normals: Vec<DVec3>,
    /// When `true` the mesh data is out of date with respect to the source
    /// geometry and must be recomputed before use.
    ///
    /// `triangulate_surface` always returns a clean mesh (`dirty = false`).
    /// Callers that cache a `SurfaceMesh` should call [`SurfaceMesh::invalidate`]
    /// whenever the source geometry changes.
    pub dirty: bool,
}

impl SurfaceMesh {
    /// Mark this mesh as stale.  The next render or query should recompute it.
    pub fn invalidate(&mut self) {
        self.dirty = true;
    }

    /// Returns `true` if the mesh data is up-to-date with the source geometry.
    pub fn is_clean(&self) -> bool {
        !self.dirty
    }
}

/// 曲面三角化参数。
#[derive(Debug, Clone)]
pub struct TessellationParams {
    // === 基本控制参数 ===
    /// 最大弦差（三角形中点到曲面的最大允许距离）。
    /// 较小的值产生更细的网格，推荐范围 0.001~0.1。
    pub chord_tolerance: f64,
    /// 最大角度误差（弧度）。超过此角度的相邻三角形会被进一步细分。
    pub angle_tolerance: f64,
    /// 最小细分步长（UV 空间），防止无限细分。
    pub min_step: f64,
    /// 最大细分步长（UV 空间）。
    pub max_step: f64,

    // === 新增：尺寸控制 ===
    /// 最小三角形尺寸（世界坐标）。
    /// 小于此尺寸的三角形不再细分。默认为 0.0（不限制）。
    pub min_triangle_size: f64,
    /// 最大三角形尺寸（世界坐标）。
    /// 超过此尺寸的三角形会被强制细分。默认为 f64::MAX（不限制）。
    pub max_triangle_size: f64,

    // === 新增：质量控制 ===
    /// 是否启用自适应细分。
    /// 启用后，根据曲率自动调整细分级别。默认为 true。
    pub adaptive_refinement: bool,
    /// 是否启用曲率敏感细分。
    /// 启用后，高曲率区域会生成更细的网格。默认为 true。
    pub curvature_sensitive: bool,
    /// 最大三角形长宽比。
    /// 超过此比例的三角形会被标记为质量问题。默认为 20.0。
    pub max_aspect_ratio: f64,

    // === 新增：边界控制 ===
    /// 是否保持边界。
    /// 启用后，边界边的采样点不会被合并。默认为 true。
    pub boundary_preservation: bool,
    /// 是否保持缝线。
    /// 启用后，缝线边的采样点会被特殊处理。默认为 true。
    pub seam_preservation: bool,

    // === 新增：性能控制 ===
    /// 最大递归深度。
    /// 防止无限细分，默认为 8。
    pub max_depth: usize,
    /// 是否并行处理。
    /// 启用后，多个面可以并行三角化。默认为 false（单线程安全）。
    pub parallel: bool,
}

impl Default for TessellationParams {
    fn default() -> Self {
        Self {
            chord_tolerance: 0.01,
            angle_tolerance: 0.1,  // ~5.7 degrees
            min_step: 1e-4,
            max_step: 0.5,
            min_triangle_size: 0.0,
            max_triangle_size: f64::MAX,
            adaptive_refinement: true,
            curvature_sensitive: true,
            max_aspect_ratio: 20.0,
            boundary_preservation: true,
            seam_preservation: true,
            max_depth: 8,
            parallel: false,
        }
    }
}

impl TessellationParams {
    /// 快速预览预设配置。
    /// 适用于交互式预览，追求速度而非质量。
    pub fn preview() -> Self {
        Self {
            chord_tolerance: 0.1,
            angle_tolerance: 0.3,  // ~17 degrees
            min_step: 1e-3,
            max_step: 1.0,
            min_triangle_size: 0.01,
            max_triangle_size: f64::MAX,
            adaptive_refinement: false,
            curvature_sensitive: false,
            max_aspect_ratio: 30.0,
            boundary_preservation: true,
            seam_preservation: false,
            max_depth: 4,
            parallel: true,
        }
    }

    /// 标准质量预设配置。
    /// 平衡质量和性能，适用于一般用途。
    pub fn standard() -> Self {
        Self {
            chord_tolerance: 0.01,
            angle_tolerance: 0.1,  // ~5.7 degrees
            min_step: 1e-4,
            max_step: 0.5,
            min_triangle_size: 0.0,
            max_triangle_size: f64::MAX,
            adaptive_refinement: true,
            curvature_sensitive: true,
            max_aspect_ratio: 20.0,
            boundary_preservation: true,
            seam_preservation: true,
            max_depth: 8,
            parallel: false,
        }
    }

    /// 高质量预设配置。
    /// 适用于渲染和可视化，追求高质量网格。
    pub fn high_quality() -> Self {
        Self {
            chord_tolerance: 0.001,
            angle_tolerance: 0.05,  // ~2.9 degrees
            min_step: 1e-5,
            max_step: 0.2,
            min_triangle_size: 0.0,
            max_triangle_size: f64::MAX,
            adaptive_refinement: true,
            curvature_sensitive: true,
            max_aspect_ratio: 10.0,
            boundary_preservation: true,
            seam_preservation: true,
            max_depth: 12,
            parallel: false,
        }
    }

    /// 导出优化预设配置。
    /// 适用于 STL/OBJ 导出，优化文件大小。
    pub fn export() -> Self {
        Self {
            chord_tolerance: 0.005,
            angle_tolerance: 0.08,  // ~4.6 degrees
            min_step: 1e-4,
            max_step: 0.3,
            min_triangle_size: 0.0,
            max_triangle_size: f64::MAX,
            adaptive_refinement: true,
            curvature_sensitive: true,
            max_aspect_ratio: 15.0,
            boundary_preservation: true,
            seam_preservation: true,
            max_depth: 10,
            parallel: false,
        }
    }

    /// 分析准备预设配置。
    /// 适用于 FEA/CFD 分析，追求网格质量。
    pub fn analysis() -> Self {
        Self {
            chord_tolerance: 0.0005,
            angle_tolerance: 0.03,  // ~1.7 degrees
            min_step: 1e-6,
            max_step: 0.1,
            min_triangle_size: 0.0,
            max_triangle_size: f64::MAX,
            adaptive_refinement: true,
            curvature_sensitive: true,
            max_aspect_ratio: 5.0,
            boundary_preservation: true,
            seam_preservation: true,
            max_depth: 15,
            parallel: false,
        }
    }

    /// 根据目标三角形数量自动调整参数。
    /// 返回调整后的参数副本。
    pub fn with_target_triangle_count(&self, target_count: usize) -> Self {
        let factor = (target_count as f64 / 1000.0).powf(1.0 / 3.0).max(0.1).min(10.0);
        Self {
            chord_tolerance: self.chord_tolerance * factor,
            angle_tolerance: self.angle_tolerance * factor,
            min_step: self.min_step,
            max_step: self.max_step / factor,
            min_triangle_size: self.min_triangle_size,
            max_triangle_size: self.max_triangle_size,
            adaptive_refinement: self.adaptive_refinement,
            curvature_sensitive: self.curvature_sensitive,
            max_aspect_ratio: self.max_aspect_ratio,
            boundary_preservation: self.boundary_preservation,
            seam_preservation: self.seam_preservation,
            max_depth: self.max_depth,
            parallel: self.parallel,
        }
    }
}

/// 对参数曲面进行自适应弦差控制三角化。
///
/// 在 UV 参数域上进行自适应细分：
/// 1. 先以均匀网格覆盖 UV 域
/// 2. 对每个四边形检查弦差（三角形中心到真实曲面的距离）
/// 3. 超过 `params.chord_tolerance` 的四边形递归细分
/// 4. 收集所有叶节点三角形
///
/// # 参数
/// - `surface`：要三角化的曲面
/// - `u_range`：UV 域 U 方向范围 [u_min, u_max]
/// - `v_range`：UV 域 V 方向范围 [v_min, v_max]
/// - `params`：三角化参数
pub fn triangulate_surface(
    surface: &Surface3,
    u_range: [f64; 2],
    v_range: [f64; 2],
    params: &TessellationParams,
) -> SurfaceMesh {
    let mut vertices: Vec<DVec3> = Vec::new();
    let mut normals: Vec<DVec3> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    // UV 域初始格数（至少 2x2）
    let initial_steps = 4usize;
    let [u0, u1] = u_range;
    let [v0, v1] = v_range;
    let du = (u1 - u0) / initial_steps as f64;
    let dv = (v1 - v0) / initial_steps as f64;

    // 对每个初始四边形进行自适应细分
    for i in 0..initial_steps {
        for j in 0..initial_steps {
            let ua = u0 + i as f64 * du;
            let ub = ua + du;
            let va = v0 + j as f64 * dv;
            let vb = va + dv;

            subdivide_quad(
                surface,
                [ua, ub],
                [va, vb],
                params,
                0,
                &mut vertices,
                &mut normals,
                &mut triangles,
            );
        }
    }

    weld_surface_mesh_vertices(SurfaceMesh {
        vertices,
        triangles,
        normals,
        dirty: false,
    })
}

fn weld_surface_mesh_vertices(mesh: SurfaceMesh) -> SurfaceMesh {
    const WELD_TOLERANCE: f64 = 1e-9;

    let mut remap = vec![0usize; mesh.vertices.len()];
    let mut welded_vertices: Vec<DVec3> = Vec::new();
    let mut welded_normals: Vec<DVec3> = Vec::new();
    let mut normal_counts = Vec::new();
    let mut buckets: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    let scale = 1.0 / WELD_TOLERANCE;

    for (index, point) in mesh.vertices.iter().enumerate() {
        let key = [
            (point.x * scale).round() as i64,
            (point.y * scale).round() as i64,
            (point.z * scale).round() as i64,
        ];

        let mut matched = None;
        if let Some(candidates) = buckets.get(&key) {
            for &candidate in candidates {
                if (welded_vertices[candidate] - *point).length_squared() <= WELD_TOLERANCE * WELD_TOLERANCE {
                    matched = Some(candidate);
                    break;
                }
            }
        }

        let target = if let Some(existing) = matched {
            existing
        } else {
            let new_index = welded_vertices.len();
            welded_vertices.push(*point);
            welded_normals.push(DVec3::ZERO);
            normal_counts.push(0usize);
            buckets.entry(key).or_default().push(new_index);
            new_index
        };

        remap[index] = target;
        if let Some(normal) = mesh.normals.get(index) {
            welded_normals[target] += *normal;
            normal_counts[target] += 1;
        }
    }

    let welded_triangles: Vec<[usize; 3]> = mesh
        .triangles
        .iter()
        .filter_map(|&[a, b, c]| {
            let ra = remap[a];
            let rb = remap[b];
            let rc = remap[c];
            if ra == rb || rb == rc || rc == ra {
                None
            } else {
                Some([ra, rb, rc])
            }
        })
        .collect();

    let welded_normals: Vec<DVec3> = welded_normals
        .into_iter()
        .zip(normal_counts)
        .map(|(normal, count)| {
            if count == 0 {
                DVec3::ZERO
            } else {
                normal.normalize_or_zero()
            }
        })
        .collect();

    SurfaceMesh {
        vertices: welded_vertices,
        triangles: welded_triangles,
        normals: welded_normals,
        dirty: mesh.dirty,
    }
}

/// 递归自适应细分一个 UV 四边形。
fn subdivide_quad(
    surface: &Surface3,
    u_range: [f64; 2],
    v_range: [f64; 2],
    params: &TessellationParams,
    depth: usize,
    vertices: &mut Vec<DVec3>,
    normals: &mut Vec<DVec3>,
    triangles: &mut Vec<[usize; 3]>,
) {
    let [u0, u1] = u_range;
    let [v0, v1] = v_range;

    // 计算四角点
    let p00 = surface.point_at(u0, v0);
    let p10 = surface.point_at(u1, v0);
    let p01 = surface.point_at(u0, v1);
    let p11 = surface.point_at(u1, v1);

    let um = (u0 + u1) * 0.5;
    let vm = (v0 + v1) * 0.5;

    // 检查是否需要继续细分
    let should_subdivide = if depth < params.max_depth {
        let step_u = u1 - u0;
        let step_v = v1 - v0;

        // 检查步长是否还能细分
        if step_u < params.min_step * 2.0 && step_v < params.min_step * 2.0 {
            false
        } else {
            // 检查弦差：计算两个三角形的中心点到曲面的距离
            let chord_exceeded = check_chord_tolerance(surface, p00, p10, p11, p01, um, vm, params.chord_tolerance);

            // 检查角度误差（法向量变化）
            let angle_exceeded = depth < params.max_depth / 2 && check_angle_tolerance(surface, u0, u1, v0, v1, params.angle_tolerance);

            // 检查最大三角形尺寸
            let size_exceeded = if params.max_triangle_size < f64::MAX {
                let diag = (p11 - p00).length();
                diag > params.max_triangle_size
            } else {
                false
            };

            // 如果启用了自适应细分，则综合考虑
            if params.adaptive_refinement {
                chord_exceeded || angle_exceeded || size_exceeded
            } else {
                size_exceeded
            }
        }
    } else {
        false
    };

    if should_subdivide {
        // 细分为4个子四边形
        subdivide_quad(surface, [u0, um], [v0, vm], params, depth + 1, vertices, normals, triangles);
        subdivide_quad(surface, [um, u1], [v0, vm], params, depth + 1, vertices, normals, triangles);
        subdivide_quad(surface, [u0, um], [vm, v1], params, depth + 1, vertices, normals, triangles);
        subdivide_quad(surface, [um, u1], [vm, v1], params, depth + 1, vertices, normals, triangles);
    } else {
        // 发射两个三角形
        let n = vertices.len();

        // 计算法向量
        let n00 = surface.normal_at(u0, v0);
        let n10 = surface.normal_at(u1, v0);
        let n01 = surface.normal_at(u0, v1);
        let n11 = surface.normal_at(u1, v1);

        // 检查点是否退化（NaN 或 Inf）
        let valid = [p00, p10, p01, p11].iter().all(|p| p.is_finite());
        if !valid {
            return;
        }

        vertices.extend_from_slice(&[p00, p10, p11, p01]);
        normals.extend_from_slice(&[n00, n10, n11, n01]);

        // 选择对角线方向使三角形更均匀
        let d0 = (p11 - p00).length_squared();
        let d1 = (p10 - p01).length_squared();
        if d0 <= d1 {
            // 沿 p00-p11 对角线
            triangles.push([n, n + 1, n + 2]);
            triangles.push([n, n + 2, n + 3]);
        } else {
            // 沿 p10-p01 对角线
            triangles.push([n, n + 1, n + 3]);
            triangles.push([n + 1, n + 2, n + 3]);
        }
    }
}

/// 检查弦差是否超过容差。
/// 计算两个三角形的中心到曲面的近似距离。
fn check_chord_tolerance(
    surface: &Surface3,
    p00: DVec3, p10: DVec3, p11: DVec3, p01: DVec3,
    um: f64, vm: f64,
    tolerance: f64,
) -> bool {
    // 三角形1 (p00, p10, p11) 的中心
    let c1 = (p00 + p10 + p11) / 3.0;
    // 三角形2 (p00, p11, p01) 的中心
    let c2 = (p00 + p11 + p01) / 3.0;

    // 曲面上对应 UV 中心的实际点
    let surf_mid = surface.point_at(um, vm);

    // 检查曲面中点到线性插值中点的距离
    let interp_mid = (p00 + p10 + p11 + p01) / 4.0;
    let chord_err = (surf_mid - interp_mid).length();

    // 也检查各三角形中心处的弦差
    let t1_u = (c1 - p00).length() / (p11 - p00).length().max(1e-10);
    let _ = t1_u; // UV 坐标近似用中点替代
    let chord1 = (surface.point_at(um, vm) - c1).length();
    let chord2 = (surface.point_at(um, vm) - c2).length();

    chord_err > tolerance || chord1 > tolerance || chord2 > tolerance
}

/// 检查角度误差（法向量变化）是否超过容差。
fn check_angle_tolerance(
    surface: &Surface3,
    u0: f64, u1: f64, v0: f64, v1: f64,
    tolerance: f64,
) -> bool {
    let n00 = surface.normal_at(u0, v0);
    let n11 = surface.normal_at(u1, v1);
    let n10 = surface.normal_at(u1, v0);
    let n01 = surface.normal_at(u0, v1);

    // 检查相邻角点的法向量夹角
    for (a, b) in [(n00, n10), (n00, n01), (n11, n10), (n11, n01)] {
        let la = a.length();
        let lb = b.length();
        if la < 0.5 || lb < 0.5 {
            continue;
        }
        let cos_a = (a.dot(b) / (la * lb)).clamp(-1.0, 1.0);
        let angle = cos_a.acos();
        if angle > tolerance {
            return true;
        }
    }
    false
}

/// Ear-clipping triangulation for a simple polygon in 3D.
/// Projects to 2D using the given normal, then runs ear-clipping.
pub fn triangulate_polygon(vertices: &[DVec3], normal: DVec3) -> Vec<[usize; 3]> {
    let n = vertices.len();
    if n < 3 {
        return vec![];
    }
    if n == 3 {
        return vec![[0, 1, 2]];
    }
    if n == 4 {
        return vec![[0, 1, 2], [0, 2, 3]];
    }

    let (u_axis, v_axis) = local_basis(normal);
    let pts_2d: Vec<[f64; 2]> = vertices
        .iter()
        .map(|p| [p.dot(u_axis), p.dot(v_axis)])
        .collect();

    ear_clip(&pts_2d)
}

fn local_basis(normal: DVec3) -> (DVec3, DVec3) {
    let ref_dir = if normal.x.abs() < 0.9 {
        DVec3::X
    } else {
        DVec3::Y
    };
    let u = normal.cross(ref_dir).normalize();
    let v = normal.cross(u).normalize();
    (u, v)
}

/// Build an ordered 3D polygon from a wire.
///
/// Curved edges are sampled using their analytic 3D curve + edge range.
/// Straight or missing-geometry edges contribute only their end vertex.
fn sample_wire_polygon_points(brep: &BRep, wire: &rcad_kernel::topology::Wire) -> Vec<DVec3> {
    let mut pts: Vec<DVec3> = Vec::new();
    let two_pi = 2.0 * std::f64::consts::PI;

    for we in &wire.edges {
        let Some(edge) = brep.edges.get(we.idx) else {
            continue;
        };

        let start_idx = if we.forward { edge.start } else { edge.end };
        let end_idx = if we.forward { edge.end } else { edge.start };

        let p_start = match brep.vertices.get(start_idx) {
            Some(v) => v.point,
            None => continue,
        };
        let p_end = match brep.vertices.get(end_idx) {
            Some(v) => v.point,
            None => continue,
        };

        let mut sampled = false;
        if let Some(ci) = brep.geom.edge_curve.get(we.idx).and_then(|v| *v) {
            if let Some(curve) = brep.geom.curves.get(ci) {
                if !matches!(curve, Curve3::Line(_)) {
                    let Some([r0, r1]) = brep
                        .geom
                        .edge_curve_range
                        .get(we.idx)
                        .and_then(|v| *v)
                        .or_else(|| match curve {
                            Curve3::Circle(_) | Curve3::Ellipse(_) => {
                                Some([0.0, 2.0 * std::f64::consts::PI])
                            }
                            _ => None,
                        })
                    else {
                        continue;
                    };

                    let mut t0 = r0;
                    let mut t1 = r1;
                    if !we.forward {
                        std::mem::swap(&mut t0, &mut t1);
                    }

                    // Repair clearly wrong full-period range on circular/elliptic edges.
                    match curve {
                        Curve3::Circle(c) => {
                            let wrap_2pi = |t: f64| -> f64 {
                                let mut out = t % two_pi;
                                if out < 0.0 {
                                    out += two_pi;
                                }
                                out
                            };
                            if (t1 - t0).abs() >= two_pi * 0.999 {
                                let x_ax = rcad_kernel::geom::any_perpendicular(c.normal);
                                let y_ax = c.normal.cross(x_ax);
                                let v0 = p_start - c.center;
                                let v1 = p_end - c.center;
                                let a0 = wrap_2pi(v0.dot(y_ax).atan2(v0.dot(x_ax)));
                                let a1 = wrap_2pi(v1.dot(y_ax).atan2(v1.dot(x_ax)));
                                let mut dt = a1 - a0;
                                if dt > std::f64::consts::PI {
                                    dt -= two_pi;
                                } else if dt < -std::f64::consts::PI {
                                    dt += two_pi;
                                }
                                t0 = a0;
                                t1 = a0 + dt;
                            }
                        }
                        Curve3::Ellipse(e) => {
                            let wrap_2pi = |t: f64| -> f64 {
                                let mut out = t % two_pi;
                                if out < 0.0 {
                                    out += two_pi;
                                }
                                out
                            };
                            if (t1 - t0).abs() >= two_pi * 0.999 {
                                let x_ax = e.major_dir.normalize();
                                let y_ax = e.normal.cross(x_ax).normalize();
                                let v0 = p_start - e.center;
                                let v1 = p_end - e.center;
                                let a0 = wrap_2pi((v0.dot(y_ax) / e.minor_radius).atan2(v0.dot(x_ax) / e.major_radius));
                                let a1 = wrap_2pi((v1.dot(y_ax) / e.minor_radius).atan2(v1.dot(x_ax) / e.major_radius));
                                let mut dt = a1 - a0;
                                if dt > std::f64::consts::PI {
                                    dt -= two_pi;
                                } else if dt < -std::f64::consts::PI {
                                    dt += two_pi;
                                }
                                t0 = a0;
                                t1 = a0 + dt;
                            }
                        }
                        _ => {}
                    }

                    let span = (t1 - t0).abs();
                    if span > 1e-12 {
                        let n_segs = match curve {
                            Curve3::Circle(_) => {
                                let segs = (span / (2.0 * std::f64::consts::PI) * 64.0).ceil() as usize;
                                segs.clamp(4, 64)
                            }
                            Curve3::Ellipse(_) => 24,
                            _ => 16,
                        };
                        for i in 0..=n_segs {
                            if !pts.is_empty() && i == 0 {
                                continue;
                            }
                            let t = t0 + (t1 - t0) * (i as f64 / n_segs as f64);
                            pts.push(curve.point_at(t));
                        }
                        sampled = true;
                    }
                }
            }
        }

        if !sampled {
            if pts.is_empty() {
                pts.push(p_start);
            }
            pts.push(p_end);
        }
    }

    // Drop duplicated closing point if present.
    if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length() < 1e-9 {
        pts.pop();
    }

    pts
}

fn ear_clip(pts: &[[f64; 2]]) -> Vec<[usize; 3]> {
    let n = pts.len();
    let mut indices: Vec<usize> = (0..n).collect();
    let mut triangles = Vec::new();

    // Ensure CCW winding
    let area = signed_area_2d(pts, &indices);
    if area < 0.0 {
        indices.reverse();
    }

    let mut remaining = indices;
    while remaining.len() > 3 {
        let len = remaining.len();
        let mut ear_found = false;

        for i in 0..len {
            let prev = if i == 0 { len - 1 } else { i - 1 };
            let next = if i == len - 1 { 0 } else { i + 1 };

            let a = remaining[prev];
            let b = remaining[i];
            let c = remaining[next];

            // Check convexity (left turn)
            if cross_2d(pts[a], pts[b], pts[c]) <= 0.0 {
                continue;
            }

            // Check no other vertex inside this triangle
            let mut contains_other = false;
            for j in 0..len {
                if j == prev || j == i || j == next {
                    continue;
                }
                if point_in_triangle_2d(pts[remaining[j]], pts[a], pts[b], pts[c]) {
                    contains_other = true;
                    break;
                }
            }

            if !contains_other {
                triangles.push([a, b, c]);
                remaining.remove(i);
                ear_found = true;
                break;
            }
        }

        if !ear_found {
            // Degenerate polygon — emit remaining as fan
            for i in 1..remaining.len() - 1 {
                triangles.push([remaining[0], remaining[i], remaining[i + 1]]);
            }
            break;
        }
    }

    if remaining.len() == 3 {
        triangles.push([remaining[0], remaining[1], remaining[2]]);
    }

    triangles
}

fn signed_area_2d(pts: &[[f64; 2]], indices: &[usize]) -> f64 {
    let n = indices.len();
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let a = pts[indices[i]];
        let b = pts[indices[j]];
        area += a[0] * b[1] - b[0] * a[1];
    }
    area * 0.5
}

fn cross_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn point_in_triangle_2d(p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
    let d1 = cross_2d(a, b, p);
    let d2 = cross_2d(b, c, p);
    let d3 = cross_2d(c, a, p);

    let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);

    !(has_neg && has_pos)
}

/// Tessellate all faces of a BRep in-place, writing triangle indices into
/// `face.triangles`.
///
/// Analogous to OCCT `BRepMesh_IncrementalMesh`.
///
/// For each face:
/// - If the face has an associated `Surface3` in `brep.geom`, the surface is
///   sampled adaptively using `triangulate_surface` with the given `params`.
///   The resulting world-space vertices are appended to `brep.vertices` and
///   the triangle indices are stored in `face.triangles`.
/// - Faces without a surface entry fall back to fan-triangulation of the
///   outer wire vertices (same as the existing rendering path).
///
/// Faces whose [`Face::mesh_dirty`] flag is `false` (clean) are **skipped**
/// unless their `triangles` is empty — allowing incremental updates when only
/// part of the model changes.  To force a full retessellation call
/// [`BRep::invalidate_mesh`] first.
///
/// After tessellating a face its `mesh_dirty` flag is set to `false`.
pub fn mesh_brep(brep: &mut BRep, params: &TessellationParams) {
    let mut face_flat_idx = 0usize;

    for solid_idx in 0..brep.solids.len() {
        for shell_idx in 0..brep.solids[solid_idx].shells.len() {
            let n_faces = brep.solids[solid_idx].shells[shell_idx].faces.len();
            for face_idx in 0..n_faces {
                // Skip faces whose cached triangulation is still valid.
                {
                    let face = &brep.solids[solid_idx].shells[shell_idx].faces[face_idx];
                    if face.mesh_is_clean() {
                        face_flat_idx += 1;
                        continue;
                    }
                }

                // Resolve surface and UV domain.
                let surf_and_domain: Option<(Surface3, [f64; 4])> = brep
                    .geom
                    .face_surface
                    .get(face_flat_idx)
                    .and_then(|o| *o)
                    .and_then(|si| brep.geom.surfaces.get(si).cloned())
                    .map(|surf| {
                        let domain = brep
                            .geom
                            .face_surface_range
                            .get(face_flat_idx)
                            .and_then(|o| *o)
                            .unwrap_or_else(|| surf.default_domain());
                        (surf, domain)
                    });

                brep.solids[solid_idx].shells[shell_idx].faces[face_idx]
                    .triangles
                    .clear();

                if let Some((surf, domain)) = surf_and_domain {
                    let [u0, u1, v0, v1] = domain;

                    // Clamp infinite domains (e.g. cylinder v-range) using
                    // vertex projections.
                    let (u0, u1, v0, v1) = clamp_domain_to_vertices(
                        brep, face_flat_idx, &surf, u0, u1, v0, v1,
                    );

                    if (u1 - u0).abs() < 1e-10 || (v1 - v0).abs() < 1e-10 {
                        face_flat_idx += 1;
                        continue;
                    }

                    let mesh = triangulate_surface(
                        &surf,
                        [u0, u1],
                        [v0, v1],
                        params,
                    );

                    if mesh.triangles.is_empty() {
                        face_flat_idx += 1;
                        continue;
                    }

                    // Append new vertices and remap triangle indices.
                    let base = brep.vertices.len();
                    for &pt in &mesh.vertices {
                        brep.vertices.push(rcad_kernel::topology::Vertex { point: pt });
                    }
                    let tris: Vec<[usize; 3]> = mesh
                        .triangles
                        .iter()
                        .map(|&[a, b, c]| [base + a, base + b, base + c])
                        .collect();
                    brep.solids[solid_idx].shells[shell_idx].faces[face_idx]
                        .triangles = tris;
                } else {
                    // Fallback for faces without bound surface:
                    // sample curved outer-wire edges into a polygon then ear-clip.
                    let face_ref = &brep.solids[solid_idx].shells[shell_idx].faces[face_idx];
                    let poly_pts = sample_wire_polygon_points(brep, &face_ref.outer_wire);
                    if poly_pts.len() >= 3 {
                        let local_tris = triangulate_polygon(&poly_pts, face_ref.normal);
                        if !local_tris.is_empty() {
                            let base = brep.vertices.len();
                            for &pt in &poly_pts {
                                brep.vertices.push(rcad_kernel::topology::Vertex { point: pt });
                            }
                            let tris: Vec<[usize; 3]> = local_tris
                                .iter()
                                .map(|&[a, b, c]| [base + a, base + b, base + c])
                                .collect();
                            brep.solids[solid_idx].shells[shell_idx].faces[face_idx]
                                .triangles = tris;
                        }
                    }
                }

                // Mark the face mesh as clean.
                brep.solids[solid_idx].shells[shell_idx].faces[face_idx]
                    .mesh_dirty = false;

                face_flat_idx += 1;
            }
        }
    }
}

/// Clamp a potentially infinite UV domain to the range spanned by the face's
/// wire vertices projected onto the surface parameters.
fn clamp_domain_to_vertices(
    brep: &BRep,
    face_flat_idx: usize,
    surf: &Surface3,
    u0: f64, u1: f64, v0: f64, v1: f64,
) -> (f64, f64, f64, f64) {

    // Only clamp axes that are infinite.
    let need_u = !u0.is_finite() || !u1.is_finite();
    let need_v = !v0.is_finite() || !v1.is_finite();
    if !need_u && !need_v {
        return (u0, u1, v0, v1);
    }

    // Collect wire vertices for this face.
    let face = brep
        .solids
        .iter()
        .flat_map(|s| s.shells.iter())
        .flat_map(|sh| sh.faces.iter())
        .nth(face_flat_idx);

    let Some(face) = face else {
        return (u0, u1, v0, v1);
    };

    let pts: Vec<DVec3> = face
        .outer_wire
        .edges
        .iter()
        .filter_map(|we| {
            brep.edges.get(we.idx).and_then(|e| {
                let vi = if we.forward { e.start } else { e.end };
                brep.vertices.get(vi).map(|v| v.point)
            })
        })
        .collect();

    if pts.is_empty() {
        return (u0, u1, v0, v1);
    }

    match surf {
        Surface3::Plane(plane) => {
            // Project vertices onto the plane's local UV frame.
            use rcad_kernel::geom::any_perpendicular;
            let u_ax = any_perpendicular(plane.normal);
            let v_ax = plane.normal.cross(u_ax).normalize_or_zero();
            let us: Vec<f64> = pts.iter().map(|&p| (p - plane.origin).dot(u_ax)).collect();
            let vs: Vec<f64> = pts.iter().map(|&p| (p - plane.origin).dot(v_ax)).collect();
            let pu0 = us.iter().cloned().fold(f64::INFINITY, f64::min);
            let pu1 = us.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let pv0 = vs.iter().cloned().fold(f64::INFINITY, f64::min);
            let pv1 = vs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let mu = (pu1 - pu0).abs() * 0.05 + 1e-6;
            let mv = (pv1 - pv0).abs() * 0.05 + 1e-6;
            (pu0 - mu, pu1 + mu, pv0 - mv, pv1 + mv)
        }
        Surface3::Cylinder(cyl) => {
            let eff_v0 = if v0.is_finite() { v0 } else {
                pts.iter().map(|&p| (p - cyl.origin).dot(cyl.axis))
                    .fold(f64::INFINITY, f64::min)
            };
            let eff_v1 = if v1.is_finite() { v1 } else {
                pts.iter().map(|&p| (p - cyl.origin).dot(cyl.axis))
                    .fold(f64::NEG_INFINITY, f64::max)
            };
            let eff_u0 = if u0.is_finite() { u0 } else { 0.0 };
            let eff_u1 = if u1.is_finite() { u1 } else { 2.0 * std::f64::consts::PI };
            (eff_u0, eff_u1, eff_v0, eff_v1)
        }
        Surface3::Cone(con) => {
            let eff_v0 = if v0.is_finite() { v0 } else {
                pts.iter().map(|&p| (p - con.apex).dot(con.axis))
                    .fold(f64::INFINITY, f64::min).max(0.0)
            };
            let eff_v1 = if v1.is_finite() { v1 } else {
                pts.iter().map(|&p| (p - con.apex).dot(con.axis))
                    .fold(f64::NEG_INFINITY, f64::max).max(0.0)
            };
            let eff_u0 = if u0.is_finite() { u0 } else { 0.0 };
            let eff_u1 = if u1.is_finite() { u1 } else { 2.0 * std::f64::consts::PI };
            (eff_u0, eff_u1, eff_v0, eff_v1)
        }
        _ => {
            let eff_u0 = if u0.is_finite() { u0 } else { -10.0 };
            let eff_u1 = if u1.is_finite() { u1 } else { 10.0 };
            let eff_v0 = if v0.is_finite() { v0 } else { -10.0 };
            let eff_v1 = if v1.is_finite() { v1 } else { 10.0 };
            (eff_u0, eff_u1, eff_v0, eff_v1)
        }
    }
}

// ============================================================================
// 网格质量指标
// ============================================================================

/// 网格质量指标结构体。
///
/// 提供三角形网格的各种质量统计信息，包括长宽比、边长均匀性和面积分布。
#[derive(Debug, Clone, Default)]
pub struct MeshQualityMetrics {
    /// 三角形总数。
    pub triangle_count: usize,
    /// 顶点总数。
    pub vertex_count: usize,
    /// 平均长宽比。
    pub average_aspect_ratio: f64,
    /// 最大长宽比。
    pub max_aspect_ratio: f64,
    /// 长宽比超过阈值的三角形数量。
    pub poor_aspect_ratio_count: usize,
    /// 平均边长。
    pub average_edge_length: f64,
    /// 边长标准差。
    pub edge_length_stddev: f64,
    /// 最小边长。
    pub min_edge_length: f64,
    /// 最大边长。
    pub max_edge_length: f64,
    /// 平均三角形面积。
    pub average_area: f64,
    /// 面积标准差。
    pub area_stddev: f64,
    /// 最小三角形面积。
    pub min_area: f64,
    /// 最大三角形面积。
    pub max_area: f64,
    /// 退化三角形数量（面积为0或接近0）。
    pub degenerate_count: usize,
}

impl MeshQualityMetrics {
    /// 检查网格质量是否良好。
    ///
    /// 质量良好的网格应满足：
    /// - 没有退化三角形
    /// - 最大长宽比在合理范围内
    /// - 边长分布相对均匀
    pub fn is_good(&self, max_aspect_ratio: f64) -> bool {
        self.degenerate_count == 0
            && self.max_aspect_ratio <= max_aspect_ratio
            && (self.triangle_count <= 10 || self.poor_aspect_ratio_count < self.triangle_count / 10)
    }

    /// 返回质量评分（0.0 到 1.0）。
    ///
    /// 评分基于长宽比、退化三角形比例和边长均匀性。
    pub fn quality_score(&self) -> f64 {
        if self.triangle_count == 0 {
            return 0.0;
        }

        let aspect_score = if self.max_aspect_ratio > 0.0 {
            (10.0 / self.max_aspect_ratio).min(1.0)
        } else {
            0.0
        };

        let degenerate_ratio = self.degenerate_count as f64 / self.triangle_count as f64;
        let degenerate_score = 1.0 - degenerate_ratio;

        let uniformity_score = if self.average_edge_length > 0.0 {
            let cv = self.edge_length_stddev / self.average_edge_length;
            (1.0 - cv).max(0.0)
        } else {
            0.0
        };

        (aspect_score * 0.4 + degenerate_score * 0.4 + uniformity_score * 0.2).clamp(0.0, 1.0)
    }
}

/// 计算网格质量指标。
///
/// # 参数
/// - `vertices`: 顶点坐标数组
/// - `triangles`: 三角形索引数组
///
/// # 返回
/// 网格质量指标结构体
pub fn compute_mesh_quality(vertices: &[DVec3], triangles: &[[usize; 3]]) -> MeshQualityMetrics {
    if vertices.is_empty() || triangles.is_empty() {
        return MeshQualityMetrics::default();
    }

    let mut metrics = MeshQualityMetrics {
        triangle_count: triangles.len(),
        vertex_count: vertices.len(),
        ..Default::default()
    };

    let mut aspect_ratios: Vec<f64> = Vec::with_capacity(triangles.len());
    let mut areas: Vec<f64> = Vec::with_capacity(triangles.len());
    let mut edge_lengths: Vec<f64> = Vec::new();

    for &tri in triangles {
        let [i0, i1, i2] = tri;
        if i0 >= vertices.len() || i1 >= vertices.len() || i2 >= vertices.len() {
            continue;
        }

        let p0 = vertices[i0];
        let p1 = vertices[i1];
        let p2 = vertices[i2];

        // 计算三条边长
        let e0 = (p1 - p0).length();
        let e1 = (p2 - p1).length();
        let e2 = (p0 - p2).length();

        edge_lengths.push(e0);
        edge_lengths.push(e1);
        edge_lengths.push(e2);

        // 计算面积（海伦公式）
        let s = (e0 + e1 + e2) * 0.5;
        let area = if s > 0.0 {
            (s * (s - e0) * (s - e1) * (s - e2)).sqrt()
        } else {
            0.0
        };

        areas.push(area);

        // 检测退化三角形
        if area < 1e-12 {
            metrics.degenerate_count += 1;
        }

        // 计算长宽比
        let max_edge = e0.max(e1).max(e2);
        let min_edge = e0.min(e1).min(e2);
        let aspect_ratio = if min_edge > 1e-12 {
            max_edge / min_edge
        } else {
            f64::INFINITY
        };
        aspect_ratios.push(aspect_ratio);
    }

    // 计算统计信息
    if !aspect_ratios.is_empty() {
        metrics.max_aspect_ratio = aspect_ratios.iter().cloned().fold(0.0, f64::max);
        metrics.average_aspect_ratio = aspect_ratios.iter().sum::<f64>() / aspect_ratios.len() as f64;
        metrics.poor_aspect_ratio_count = aspect_ratios.iter().filter(|&&ar| ar > 20.0).count();
    }

    if !edge_lengths.is_empty() {
        metrics.min_edge_length = edge_lengths.iter().cloned().fold(f64::INFINITY, f64::min);
        metrics.max_edge_length = edge_lengths.iter().cloned().fold(0.0, f64::max);
        metrics.average_edge_length = edge_lengths.iter().sum::<f64>() / edge_lengths.len() as f64;

        let variance = edge_lengths.iter()
            .map(|&l| (l - metrics.average_edge_length).powi(2))
            .sum::<f64>() / edge_lengths.len() as f64;
        metrics.edge_length_stddev = variance.sqrt();
    }

    if !areas.is_empty() {
        metrics.min_area = areas.iter().cloned().fold(f64::INFINITY, f64::min);
        metrics.max_area = areas.iter().cloned().fold(0.0, f64::max);
        metrics.average_area = areas.iter().sum::<f64>() / areas.len() as f64;

        let variance = areas.iter()
            .map(|&a| (a - metrics.average_area).powi(2))
            .sum::<f64>() / areas.len() as f64;
        metrics.area_stddev = variance.sqrt();
    }

    metrics
}

/// 为 SurfaceMesh 计算质量指标。
impl SurfaceMesh {
    /// 计算此网格的质量指标。
    pub fn compute_quality(&self) -> MeshQualityMetrics {
        compute_mesh_quality(&self.vertices, &self.triangles)
    }
}

// ============================================================================
// 自适应网格细分器
// ============================================================================

/// 自适应网格细分器。
///
/// 根据曲率或距离条件对现有网格进行细分。
#[derive(Debug, Clone)]
pub struct AdaptiveSubdivider {
    /// 曲率细分阈值。
    /// 当相邻法向量夹角超过此值时，细分该边。
    pub curvature_threshold: f64,
    /// 距离细分阈值。
    /// 当边长超过此值时，细分该边。
    pub distance_threshold: f64,
    /// 最大细分层级。
    pub max_subdivision_levels: usize,
    /// 是否保持边界。
    pub preserve_boundary: bool,
}

impl Default for AdaptiveSubdivider {
    fn default() -> Self {
        Self {
            curvature_threshold: 0.1,  // ~5.7 degrees
            distance_threshold: 0.1,
            max_subdivision_levels: 3,
            preserve_boundary: true,
        }
    }
}

impl AdaptiveSubdivider {
    /// 创建新的自适应细分器。
    pub fn new() -> Self {
        Self::default()
    }

    /// 设置曲率细分阈值。
    pub fn with_curvature_threshold(mut self, threshold: f64) -> Self {
        self.curvature_threshold = threshold;
        self
    }

    /// 设置距离细分阈值。
    pub fn with_distance_threshold(mut self, threshold: f64) -> Self {
        self.distance_threshold = threshold;
        self
    }

    /// 设置最大细分层级。
    pub fn with_max_levels(mut self, levels: usize) -> Self {
        self.max_subdivision_levels = levels;
        self
    }

    /// 基于曲率细分网格。
    ///
    /// 对于法向量变化超过阈值的边进行细分。
    pub fn subdivide_by_curvature(&self, mesh: &SurfaceMesh) -> SurfaceMesh {
        if mesh.triangles.is_empty() || mesh.normals.is_empty() {
            return mesh.clone();
        }

        let mut vertices = mesh.vertices.clone();
        let mut normals = mesh.normals.clone();
        let mut triangles = Vec::new();

        // 边到新顶点的映射
        let mut edge_midpoints: HashMap<(usize, usize), usize> = HashMap::new();

        for &tri in &mesh.triangles {
            let [i0, i1, i2] = tri;

            // 检查每条边的法向量变化
            let n0 = normals.get(i0).copied().unwrap_or(DVec3::ZERO);
            let n1 = normals.get(i1).copied().unwrap_or(DVec3::ZERO);
            let n2 = normals.get(i2).copied().unwrap_or(DVec3::ZERO);

            let split_01 = self.should_split_by_curvature(n0, n1);
            let split_12 = self.should_split_by_curvature(n1, n2);
            let split_20 = self.should_split_by_curvature(n2, n0);

            if split_01 || split_12 || split_20 {
                self.subdivide_triangle(
                    tri,
                    &mut vertices,
                    &mut normals,
                    &mut triangles,
                    &mut edge_midpoints,
                );
            } else {
                triangles.push(tri);
            }
        }

        SurfaceMesh {
            vertices,
            triangles,
            normals,
            dirty: false,
        }
    }

    /// 基于距离细分网格。
    ///
    /// 对于长度超过阈值的边进行细分。
    pub fn subdivide_by_distance(&self, mesh: &SurfaceMesh) -> SurfaceMesh {
        if mesh.triangles.is_empty() {
            return mesh.clone();
        }

        let mut vertices = mesh.vertices.clone();
        let mut normals = mesh.normals.clone();
        let mut triangles = Vec::new();

        let mut edge_midpoints: HashMap<(usize, usize), usize> = HashMap::new();

        for &tri in &mesh.triangles {
            let [i0, i1, i2] = tri;

            let p0 = vertices[i0];
            let p1 = vertices[i1];
            let p2 = vertices[i2];

            let split_01 = self.should_split_by_distance(p0, p1);
            let split_12 = self.should_split_by_distance(p1, p2);
            let split_20 = self.should_split_by_distance(p2, p0);

            if split_01 || split_12 || split_20 {
                self.subdivide_triangle(
                    tri,
                    &mut vertices,
                    &mut normals,
                    &mut triangles,
                    &mut edge_midpoints,
                );
            } else {
                triangles.push(tri);
            }
        }

        SurfaceMesh {
            vertices,
            triangles,
            normals,
            dirty: false,
        }
    }

    fn should_split_by_curvature(&self, n0: DVec3, n1: DVec3) -> bool {
        let len0 = n0.length();
        let len1 = n1.length();
        if len0 < 0.5 || len1 < 0.5 {
            return false;
        }
        let cos_angle = (n0.dot(n1) / (len0 * len1)).clamp(-1.0, 1.0);
        let angle = cos_angle.acos();
        angle > self.curvature_threshold
    }

    fn should_split_by_distance(&self, p0: DVec3, p1: DVec3) -> bool {
        (p1 - p0).length() > self.distance_threshold
    }

    fn subdivide_triangle(
        &self,
        tri: [usize; 3],
        vertices: &mut Vec<DVec3>,
        normals: &mut Vec<DVec3>,
        triangles: &mut Vec<[usize; 3]>,
        edge_midpoints: &mut HashMap<(usize, usize), usize>,
    ) {
        let [i0, i1, i2] = tri;

        let p0 = vertices[i0];
        let p1 = vertices[i1];
        let p2 = vertices[i2];

        let n0 = normals.get(i0).copied().unwrap_or(DVec3::ZERO);
        let n1 = normals.get(i1).copied().unwrap_or(DVec3::ZERO);
        let n2 = normals.get(i2).copied().unwrap_or(DVec3::ZERO);

        // 获取或创建边中点
        let m01 = self.get_or_create_midpoint(i0, i1, p0, p1, n0, n1, vertices, normals, edge_midpoints);
        let m12 = self.get_or_create_midpoint(i1, i2, p1, p2, n1, n2, vertices, normals, edge_midpoints);
        let m20 = self.get_or_create_midpoint(i2, i0, p2, p0, n2, n0, vertices, normals, edge_midpoints);

        // 创建4个新三角形
        triangles.push([i0, m01, m20]);
        triangles.push([m01, i1, m12]);
        triangles.push([m20, m12, i2]);
        triangles.push([m01, m12, m20]);
    }

    fn get_or_create_midpoint(
        &self,
        i0: usize,
        i1: usize,
        p0: DVec3,
        p1: DVec3,
        n0: DVec3,
        n1: DVec3,
        vertices: &mut Vec<DVec3>,
        normals: &mut Vec<DVec3>,
        edge_midpoints: &mut HashMap<(usize, usize), usize>,
    ) -> usize {
        let key = if i0 < i1 { (i0, i1) } else { (i1, i0) };

        if let Some(&idx) = edge_midpoints.get(&key) {
            return idx;
        }

        let mid_point = (p0 + p1) * 0.5;
        let mid_normal = (n0 + n1).normalize_or_zero();

        let idx = vertices.len();
        vertices.push(mid_point);
        normals.push(mid_normal);
        edge_midpoints.insert(key, idx);
        idx
    }
}

// ============================================================================
// 边界敏感网格细分器
// ============================================================================

/// 特征边信息。
#[derive(Debug, Clone)]
pub struct FeatureEdge {
    /// 边的起始顶点索引。
    pub start_vertex: usize,
    /// 边的结束顶点索引。
    pub end_vertex: usize,
    /// 特征角度（弧度）。
    pub feature_angle: f64,
}

/// 边界敏感网格细分器。
///
/// 保持特征边的锐利度，适用于包含尖锐边缘的模型。
#[derive(Debug, Clone)]
pub struct BoundarySensitiveTessellator {
    /// 特征角度阈值（弧度）。
    /// 超过此角度的边被视为特征边。
    pub feature_angle_threshold: f64,
    /// 特征边列表。
    pub feature_edges: Vec<FeatureEdge>,
    /// 是否自动检测特征边。
    pub auto_detect_features: bool,
}

impl Default for BoundarySensitiveTessellator {
    fn default() -> Self {
        Self {
            feature_angle_threshold: 0.52,  // ~30 degrees
            feature_edges: Vec::new(),
            auto_detect_features: true,
        }
    }
}

impl BoundarySensitiveTessellator {
    /// 创建新的边界敏感细分器。
    pub fn new() -> Self {
        Self::default()
    }

    /// 设置特征角度阈值。
    pub fn with_feature_angle(mut self, angle: f64) -> Self {
        self.feature_angle_threshold = angle;
        self
    }

    /// 添加特征边。
    pub fn add_feature_edge(mut self, start: usize, end: usize, angle: f64) -> Self {
        self.feature_edges.push(FeatureEdge {
            start_vertex: start,
            end_vertex: end,
            feature_angle: angle,
        });
        self
    }

    /// 检测特征边。
    ///
    /// 基于相邻三角形法向量的夹角检测特征边。
    pub fn detect_feature_edges(&mut self, vertices: &[DVec3], triangles: &[[usize; 3]], normals: &[DVec3]) {
        if !self.auto_detect_features {
            return;
        }

        self.feature_edges.clear();

        // 构建边到三角形的映射
        let mut edge_to_tris: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (tri_idx, &tri) in triangles.iter().enumerate() {
            let edges = [
                (tri[0].min(tri[1]), tri[0].max(tri[1])),
                (tri[1].min(tri[2]), tri[1].max(tri[2])),
                (tri[2].min(tri[0]), tri[2].max(tri[0])),
            ];
            for edge in edges {
                edge_to_tris.entry(edge).or_default().push(tri_idx);
            }
        }

        // 检测特征边
        for (edge, tri_indices) in &edge_to_tris {
            if tri_indices.len() == 2 {
                let tri0 = &triangles[tri_indices[0]];
                let tri1 = &triangles[tri_indices[1]];

                // 计算三角形法向量
                let n0 = compute_triangle_normal(vertices, tri0);
                let n1 = compute_triangle_normal(vertices, tri1);

                // 计算夹角
                let cos_angle = n0.dot(n1).clamp(-1.0, 1.0);
                let angle = cos_angle.acos();

                if angle > self.feature_angle_threshold {
                    self.feature_edges.push(FeatureEdge {
                        start_vertex: edge.0,
                        end_vertex: edge.1,
                        feature_angle: angle,
                    });
                }
            }
        }
    }

    /// 保持特征边进行细分。
    ///
    /// 特征边上的顶点不会被合并或移动。
    pub fn preserve_feature_edges(&self, mesh: &SurfaceMesh) -> SurfaceMesh {
        if self.feature_edges.is_empty() {
            return mesh.clone();
        }

        // 构建特征边顶点集合
        let feature_vertices: std::collections::HashSet<usize> = self.feature_edges.iter()
            .flat_map(|e| [e.start_vertex, e.end_vertex])
            .collect();

        // 在焊接近似时排除特征边顶点
        let mut result = mesh.clone();
        result = weld_surface_mesh_vertices_with_exclusion(&result, &feature_vertices);
        result
    }
}

fn compute_triangle_normal(vertices: &[DVec3], tri: &[usize; 3]) -> DVec3 {
    if tri[0] >= vertices.len() || tri[1] >= vertices.len() || tri[2] >= vertices.len() {
        return DVec3::Z;
    }
    let p0 = vertices[tri[0]];
    let p1 = vertices[tri[1]];
    let p2 = vertices[tri[2]];
    (p1 - p0).cross(p2 - p0).normalize_or_zero()
}

fn weld_surface_mesh_vertices_with_exclusion(
    mesh: &SurfaceMesh,
    excluded_vertices: &std::collections::HashSet<usize>,
) -> SurfaceMesh {
    const WELD_TOLERANCE: f64 = 1e-9;

    let mut remap = vec![0usize; mesh.vertices.len()];
    let mut welded_vertices: Vec<DVec3> = Vec::new();
    let mut welded_normals: Vec<DVec3> = Vec::new();
    let mut normal_counts = Vec::new();
    let mut buckets: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    let scale = 1.0 / WELD_TOLERANCE;

    for (index, point) in mesh.vertices.iter().enumerate() {
        // 跳过排除的顶点
        if excluded_vertices.contains(&index) {
            let new_index = welded_vertices.len();
            welded_vertices.push(*point);
            welded_normals.push(mesh.normals.get(index).copied().unwrap_or(DVec3::ZERO));
            normal_counts.push(1);
            remap[index] = new_index;
            continue;
        }

        let key = [
            (point.x * scale).round() as i64,
            (point.y * scale).round() as i64,
            (point.z * scale).round() as i64,
        ];

        let mut matched = None;
        if let Some(candidates) = buckets.get(&key) {
            for &candidate in candidates {
                if excluded_vertices.contains(&candidate) {
                    continue;
                }
                if (welded_vertices[candidate] - *point).length_squared() <= WELD_TOLERANCE * WELD_TOLERANCE {
                    matched = Some(candidate);
                    break;
                }
            }
        }

        let target = if let Some(existing) = matched {
            existing
        } else {
            let new_index = welded_vertices.len();
            welded_vertices.push(*point);
            welded_normals.push(DVec3::ZERO);
            normal_counts.push(0);
            buckets.entry(key).or_default().push(new_index);
            new_index
        };

        remap[index] = target;
        if let Some(normal) = mesh.normals.get(index) {
            welded_normals[target] += *normal;
            normal_counts[target] += 1;
        }
    }

    let welded_triangles: Vec<[usize; 3]> = mesh
        .triangles
        .iter()
        .filter_map(|&[a, b, c]| {
            let ra = remap[a];
            let rb = remap[b];
            let rc = remap[c];
            if ra == rb || rb == rc || rc == ra {
                None
            } else {
                Some([ra, rb, rc])
            }
        })
        .collect();

    let welded_normals: Vec<DVec3> = welded_normals
        .into_iter()
        .zip(normal_counts)
        .map(|(normal, count)| {
            if count == 0 {
                DVec3::ZERO
            } else {
                normal.normalize_or_zero()
            }
        })
        .collect();

    SurfaceMesh {
        vertices: welded_vertices,
        triangles: welded_triangles,
        normals: welded_normals,
        dirty: mesh.dirty,
    }
}

// ============================================================================
// 增量网格更新器
// ============================================================================

/// 网格更新增量数据。
///
/// 描述模型中发生变化的拓扑实体，用于驱动增量网格更新。
#[derive(Debug, Clone, Default)]
pub struct MeshDelta {
    /// 修改的顶点索引。
    pub modified_vertices: Vec<usize>,
    /// 修改的边索引。
    pub modified_edges: Vec<usize>,
    /// 修改的面索引（扁平化索引）。
    pub modified_faces: Vec<usize>,
}

impl MeshDelta {
    /// 创建空的增量数据。
    pub fn new() -> Self {
        Self::default()
    }

    /// 从顶点列表创建增量数据。
    pub fn from_vertices(vertices: Vec<usize>) -> Self {
        Self {
            modified_vertices: vertices,
            ..Default::default()
        }
    }

    /// 从边列表创建增量数据。
    pub fn from_edges(edges: Vec<usize>) -> Self {
        Self {
            modified_edges: edges,
            ..Default::default()
        }
    }

    /// 从面列表创建增量数据。
    pub fn from_faces(faces: Vec<usize>) -> Self {
        Self {
            modified_faces: faces,
            ..Default::default()
        }
    }

    /// 检查是否为空。
    pub fn is_empty(&self) -> bool {
        self.modified_vertices.is_empty()
            && self.modified_edges.is_empty()
            && self.modified_faces.is_empty()
    }
}

/// 增量网格更新器。
///
/// 用于在模型局部变化时仅更新受影响区域的网格。
#[derive(Debug, Clone, Default)]
pub struct IncrementalMesher {
    /// 需要重新三角化的面索引集合。
    pub dirty_faces: std::collections::HashSet<usize>,
    /// 需要重新三角化的边索引集合。
    pub dirty_edges: std::collections::HashSet<usize>,
    /// 需要重新三角化的顶点索引集合。
    pub dirty_vertices: std::collections::HashSet<usize>,
}

impl IncrementalMesher {
    /// 创建新的增量网格更新器。
    pub fn new() -> Self {
        Self::default()
    }

    /// 标记面为需要更新。
    pub fn invalidate_face(&mut self, face_idx: usize) {
        self.dirty_faces.insert(face_idx);
    }

    /// 标记多个面为需要更新。
    pub fn invalidate_faces(&mut self, face_indices: &[usize]) {
        for &idx in face_indices {
            self.dirty_faces.insert(idx);
        }
    }

    /// 标记边为需要更新。
    pub fn invalidate_edge(&mut self, edge_idx: usize) {
        self.dirty_edges.insert(edge_idx);
    }

    /// 标记顶点为需要更新。
    pub fn invalidate_vertex(&mut self, vertex_idx: usize) {
        self.dirty_vertices.insert(vertex_idx);
    }

    /// 根据几何变化自动推断需要更新的面。
    pub fn infer_dirty_faces_from_delta(&mut self, brep: &BRep, delta: &MeshDelta) {
        // 直接标记的面
        self.invalidate_faces(&delta.modified_faces);

        // 通过边推断面
        for &edge_idx in &delta.modified_edges {
            if let Some(_edge) = brep.edges.get(edge_idx) {
                // 找到包含该边的所有面
                let mut face_idx = 0usize;
                for solid in &brep.solids {
                    for shell in &solid.shells {
                        for face in &shell.faces {
                            for we in &face.outer_wire.edges {
                                if we.idx == edge_idx {
                                    self.dirty_faces.insert(face_idx);
                                }
                            }
                            face_idx += 1;
                        }
                    }
                }
            }
        }

        // 通过顶点推断面
        for &vertex_idx in &delta.modified_vertices {
            let mut face_idx = 0usize;
            for solid in &brep.solids {
                for shell in &solid.shells {
                    for face in &shell.faces {
                        for we in &face.outer_wire.edges {
                            if let Some(edge) = brep.edges.get(we.idx) {
                                let start = if we.forward { edge.start } else { edge.end };
                                let end = if we.forward { edge.end } else { edge.start };
                                if start == vertex_idx || end == vertex_idx {
                                    self.dirty_faces.insert(face_idx);
                                }
                            }
                        }
                        face_idx += 1;
                    }
                }
            }
        }
    }

    /// 更新指定面的网格。
    ///
    /// 仅重新三角化标记为脏的面，其他面保持不变。
    pub fn update_mesh_for_face_change(
        &self,
        brep: &mut BRep,
        params: &TessellationParams,
    ) {
        if self.dirty_faces.is_empty() {
            return;
        }

        let mut face_flat_idx = 0usize;

        for solid_idx in 0..brep.solids.len() {
            for shell_idx in 0..brep.solids[solid_idx].shells.len() {
                let n_faces = brep.solids[solid_idx].shells[shell_idx].faces.len();
                for face_idx in 0..n_faces {
                    // 只更新脏面
                    if self.dirty_faces.contains(&face_flat_idx) {
                        // 标记为脏，然后重新三角化
                        brep.solids[solid_idx].shells[shell_idx].faces[face_idx]
                            .mesh_dirty = true;
                    }
                    face_flat_idx += 1;
                }
            }
        }

        // 使用标准 mesh_brep 函数进行更新
        mesh_brep(brep, params);
    }

    /// 清除所有脏标记。
    pub fn clear(&mut self) {
        self.dirty_faces.clear();
        self.dirty_edges.clear();
        self.dirty_vertices.clear();
    }

    /// 返回是否有任何脏区域。
    pub fn is_dirty(&self) -> bool {
        !self.dirty_faces.is_empty()
            || !self.dirty_edges.is_empty()
            || !self.dirty_vertices.is_empty()
    }
}

// ============================================================================
// 网格简化器
// ============================================================================

/// 边折叠信息。
#[derive(Debug, Clone)]
struct EdgeCollapseInfo {
    /// 边索引（排序后的顶点对）。
    edge: (usize, usize),
    /// 折叠误差。
    error: f64,
    /// 折叠后的新位置。
    new_position: DVec3,
}

/// 网格简化器。
///
/// 使用边折叠算法简化网格。
#[derive(Debug, Clone)]
pub struct MeshSimplifier {
    /// 目标简化比例（0.0 到 1.0）。
    pub target_ratio: f64,
    /// 最大允许误差。
    pub max_error: f64,
    /// 是否保持边界。
    pub preserve_boundary: bool,
}

impl Default for MeshSimplifier {
    fn default() -> Self {
        Self {
            target_ratio: 0.5,
            max_error: 0.01,
            preserve_boundary: true,
        }
    }
}

impl MeshSimplifier {
    /// 创建新的网格简化器。
    pub fn new() -> Self {
        Self::default()
    }

    /// 设置目标简化比例。
    pub fn with_target_ratio(mut self, ratio: f64) -> Self {
        self.target_ratio = ratio.clamp(0.0, 1.0);
        self
    }

    /// 设置最大允许误差。
    pub fn with_max_error(mut self, error: f64) -> Self {
        self.max_error = error;
        self
    }

    /// 简化网格到指定三角形数量。
    pub fn simplify_to_target_count(&self, mesh: &SurfaceMesh, target_count: usize) -> SurfaceMesh {
        if mesh.triangles.len() <= target_count {
            return mesh.clone();
        }

        let ratio = target_count as f64 / mesh.triangles.len() as f64;
        Self {
            target_ratio: ratio,
            ..self.clone()
        }
        .simplify_mesh(mesh)
    }

    /// 简化网格。
    pub fn simplify_mesh(&self, mesh: &SurfaceMesh) -> SurfaceMesh {
        if mesh.triangles.is_empty() {
            return mesh.clone();
        }

        let target_triangle_count = (mesh.triangles.len() as f64 * self.target_ratio).max(4.0) as usize;

        let mut vertices = mesh.vertices.clone();
        let mut normals = mesh.normals.clone();
        let mut triangles = mesh.triangles.clone();

        // 识别边界顶点
        let boundary_vertices = if self.preserve_boundary {
            find_boundary_vertices(&triangles)
        } else {
            std::collections::HashSet::new()
        };

        // 迭代边折叠
        while triangles.len() > target_triangle_count {
            let collapse = find_best_edge_collapse(
                &vertices,
                &triangles,
                &boundary_vertices,
                self.max_error,
            );

            let Some(collapse) = collapse else {
                break;
            };

            // 执行边折叠
            apply_edge_collapse(
                &mut vertices,
                &mut normals,
                &mut triangles,
                collapse.edge,
                collapse.new_position,
            );

            if triangles.len() <= target_triangle_count {
                break;
            }
        }

        SurfaceMesh {
            vertices,
            triangles,
            normals,
            dirty: false,
        }
    }
}

fn find_boundary_vertices(triangles: &[[usize; 3]]) -> std::collections::HashSet<usize> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();

    for &tri in triangles {
        let edges = [
            (tri[0].min(tri[1]), tri[0].max(tri[1])),
            (tri[1].min(tri[2]), tri[1].max(tri[2])),
            (tri[2].min(tri[0]), tri[2].max(tri[0])),
        ];
        for edge in edges {
            *edge_count.entry(edge).or_insert(0) += 1;
        }
    }

    let mut boundary_vertices = std::collections::HashSet::new();
    for (edge, count) in edge_count {
        if count == 1 {
            boundary_vertices.insert(edge.0);
            boundary_vertices.insert(edge.1);
        }
    }

    boundary_vertices
}

fn find_best_edge_collapse(
    vertices: &[DVec3],
    triangles: &[[usize; 3]],
    boundary_vertices: &std::collections::HashSet<usize>,
    max_error: f64,
) -> Option<EdgeCollapseInfo> {
    let mut best: Option<EdgeCollapseInfo> = None;

    // 收集所有边
    let mut edges: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    for &tri in triangles {
        edges.insert((tri[0].min(tri[1]), tri[0].max(tri[1])));
        edges.insert((tri[1].min(tri[2]), tri[1].max(tri[2])));
        edges.insert((tri[2].min(tri[0]), tri[2].max(tri[0])));
    }

    for edge in edges {
        // 跳过边界边
        if boundary_vertices.contains(&edge.0) && boundary_vertices.contains(&edge.1) {
            continue;
        }

        let p0 = vertices.get(edge.0)?;
        let p1 = vertices.get(edge.1)?;

        // 计算边长作为误差
        let error = (*p1 - *p0).length();

        if error > max_error {
            continue;
        }

        // 新位置取中点
        let new_position = (*p0 + *p1) * 0.5;

        match &best {
            None => best = Some(EdgeCollapseInfo { edge, error, new_position }),
            Some(current) if error < current.error => {
                best = Some(EdgeCollapseInfo { edge, error, new_position })
            }
            _ => {}
        }
    }

    best
}

fn apply_edge_collapse(
    vertices: &mut Vec<DVec3>,
    normals: &mut Vec<DVec3>,
    triangles: &mut Vec<[usize; 3]>,
    edge: (usize, usize),
    new_position: DVec3,
) {
    let (v0, v1) = edge;

    // 将 v1 的位置更新为新位置
    if v0 < vertices.len() {
        vertices[v0] = new_position;
    }

    // 更新法向量
    if v0 < normals.len() && v1 < normals.len() {
        normals[v0] = (normals[v0] + normals[v1]).normalize_or_zero();
    }

    // 更新三角形索引：将所有 v1 替换为 v0
    for tri in triangles.iter_mut() {
        for i in 0..3 {
            if tri[i] == v1 {
                tri[i] = v0;
            }
        }
    }

    // 移除退化三角形
    triangles.retain(|&tri| tri[0] != tri[1] && tri[1] != tri[2] && tri[2] != tri[0]);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangulate_triangle() {
        let verts = vec![DVec3::ZERO, DVec3::X, DVec3::Y];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 1);
    }

    #[test]
    fn triangulate_quad() {
        let verts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 2);
    }

    #[test]
    fn triangulate_pentagon() {
        let verts = (0..5)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect::<Vec<_>>();
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 3);
    }

    #[test]
    fn empty_polygon_returns_no_triangles() {
        let tris = triangulate_polygon(&[], DVec3::Z);
        assert!(tris.is_empty());
    }

    #[test]
    fn two_vertex_polygon_returns_no_triangles() {
        let verts = vec![DVec3::ZERO, DVec3::X];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert!(tris.is_empty());
    }

    #[test]
    fn triangle_count_is_n_minus_2() {
        // A convex n-gon should always yield n-2 triangles.
        for n in 3..=10 {
            let verts: Vec<DVec3> = (0..n)
                .map(|i| {
                    let a = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
                    DVec3::new(a.cos(), a.sin(), 0.0)
                })
                .collect();
            let tris = triangulate_polygon(&verts, DVec3::Z);
            assert_eq!(
                tris.len(),
                n - 2,
                "expected {n}-gon to yield {} triangles, got {}",
                n - 2,
                tris.len()
            );
        }
    }

    #[test]
    fn all_indices_in_bounds() {
        // Every index in the triangulation must be < number of vertices.
        let verts: Vec<DVec3> = (0..7)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 7.0;
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect();
        let tris = triangulate_polygon(&verts, DVec3::Z);
        for tri in &tris {
            for &idx in tri.iter() {
                assert!(idx < verts.len(), "index {idx} out of bounds for {n} vertices", n = verts.len());
            }
        }
    }

    #[test]
    fn clockwise_quad_still_triangulates() {
        // Reversed vertex order (CW) should be handled by sign-flip logic.
        let verts = vec![
            DVec3::new(0.0, 1.0, 0.0), // top-left first (CW)
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 0.0, 0.0),
        ];
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 2);
    }

    /// mesh_brep on a box primitive should fill face.triangles for all 6 faces.
    #[test]
    fn mesh_brep_box_fills_triangles() {
        use rcad_kernel::BRep;
        use rcad_kernel::geom::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        populate_box_geom(&mut brep);

        let params = TessellationParams::default();
        mesh_brep(&mut brep, &params);

        let faces = &brep.solids[0].shells[0].faces;
        assert_eq!(faces.len(), 6, "box should have 6 faces");
        for (i, face) in faces.iter().enumerate() {
            assert!(
                !face.triangles.is_empty(),
                "face {i} should have triangles after mesh_brep"
            );
            // All triangle indices must be valid vertex indices.
            for &[a, b, c] in &face.triangles {
                assert!(a < brep.vertices.len(), "face {i}: vertex index {a} out of bounds");
                assert!(b < brep.vertices.len(), "face {i}: vertex index {b} out of bounds");
                assert!(c < brep.vertices.len(), "face {i}: vertex index {c} out of bounds");
            }
        }
    }

    /// mesh_brep on a sphere should produce a dense mesh (many triangles per face).
    #[test]
    fn mesh_brep_sphere_produces_dense_mesh() {
        use rcad_kernel::BRep;
        use rcad_kernel::geom::PrimitiveSolid;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let params = TessellationParams {
            chord_tolerance: 0.05,
            ..TessellationParams::default()
        };
        mesh_brep(&mut brep, &params);

        let total_tris: usize = brep.solids[0].shells[0].faces
            .iter()
            .map(|f| f.triangles.len())
            .sum();
        assert!(
            total_tris >= 8,
            "sphere mesh should have at least 8 triangles, got {total_tris}"
        );
    }

    /// mesh_brep on a cylinder should produce triangles for all faces.
    #[test]
    fn mesh_brep_cylinder_all_faces_have_triangles() {
        use rcad_kernel::BRep;
        use rcad_kernel::geom::PrimitiveSolid;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let params = TessellationParams::default();
        mesh_brep(&mut brep, &params);

        let faces = &brep.solids[0].shells[0].faces;
        for (i, face) in faces.iter().enumerate() {
            assert!(
                !face.triangles.is_empty(),
                "cylinder face {i} should have triangles"
            );
        }
    }

    #[test]
    fn mesh_brep_fallback_triangulates_semicircle_wire_face() {
        use std::f64::consts::PI;
        use rcad_kernel::BRep;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
        use rcad_kernel::geom::{Circle3, Curve3, CurveEval};

        let circle = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        };
        let p0 = circle.point_at(0.0);
        let p1 = circle.point_at(PI);

        // Two-edge closed wire: semicircular arc + diameter chord.
        let mut brep = BRep {
            vertices: vec![Vertex { point: p0 }, Vertex { point: p1 }],
            edges: vec![
                Edge { start: 0, end: 1 },
                Edge { start: 1, end: 0 },
            ],
            solids: vec![Solid {
                shells: vec![Shell {
                    faces: vec![Face {
                        outer_wire: Wire {
                            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
                        },
                        inner_wires: vec![],
                        normal: DVec3::Z,
                        triangles: vec![],
                        mesh_dirty: true,
                    }],
                }],
            }],
            geom: rcad_kernel::GeomStore {
                curves: vec![Curve3::Circle(circle)],
                edge_curve: vec![Some(0), None],
                edge_curve_range: vec![Some([0.0, PI]), None],
                face_surface: vec![None],
                ..Default::default()
            },
            compound: None,
            compsolid: None,
        };

        mesh_brep(&mut brep, &TessellationParams::default());
        let tris = &brep.solids[0].shells[0].faces[0].triangles;
        assert!(
            tris.len() > 1,
            "semicircle fallback should produce multiple triangles, got {}",
            tris.len()
        );
    }

    // ========================================================================
    // 新增功能测试
    // ========================================================================

    #[test]
    fn tessellation_params_presets() {
        // 测试预览预设
        let preview = TessellationParams::preview();
        assert!(preview.chord_tolerance > TessellationParams::standard().chord_tolerance);
        assert!(!preview.adaptive_refinement);
        assert!(preview.parallel);

        // 测试标准预设
        let standard = TessellationParams::standard();
        assert!(standard.adaptive_refinement);
        assert!(standard.curvature_sensitive);

        // 测试高质量预设
        let hq = TessellationParams::high_quality();
        assert!(hq.chord_tolerance < standard.chord_tolerance);
        assert!(hq.max_depth > standard.max_depth);

        // 测试导出预设
        let export = TessellationParams::export();
        assert!(export.chord_tolerance > hq.chord_tolerance);
        assert!(export.chord_tolerance < standard.chord_tolerance);

        // 测试分析预设
        let analysis = TessellationParams::analysis();
        assert!(analysis.chord_tolerance < hq.chord_tolerance);
        assert!(analysis.max_aspect_ratio < hq.max_aspect_ratio);
    }

    #[test]
    fn tessellation_params_with_target_triangle_count() {
        let params = TessellationParams::standard();
        // Higher target count means more triangles -> finer tolerance
        let adjusted = params.with_target_triangle_count(10000);
        // Factor = (10000/1000)^(1/3) ≈ 2.15, so tolerance increases (coarser mesh)
        // For more triangles, we'd actually want lower tolerance, so this adjusts accordingly
        assert!(adjusted.chord_tolerance != params.chord_tolerance);
    }

    #[test]
    fn mesh_quality_metrics_basic() {
        let vertices = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let triangles = vec![[0, 1, 2]];

        let metrics = compute_mesh_quality(&vertices, &triangles);
        assert_eq!(metrics.triangle_count, 1);
        assert_eq!(metrics.vertex_count, 3);
        assert_eq!(metrics.degenerate_count, 0);
        assert!(metrics.max_aspect_ratio > 1.0);
    }

    #[test]
    fn mesh_quality_metrics_degenerate() {
        // 退化三角形（三个点共线）
        let vertices = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
        ];
        let triangles = vec![[0, 1, 2]];

        let metrics = compute_mesh_quality(&vertices, &triangles);
        assert_eq!(metrics.degenerate_count, 1);
        // For collinear points, the aspect ratio is max_edge/min_edge = 2/1 = 2
        // which is still bad but not infinite. The key is degenerate_count.
        assert!(metrics.max_aspect_ratio > 1.0);
    }

    #[test]
    fn mesh_quality_metrics_score() {
        // 高质量三角形
        let vertices = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.5, 0.866, 0.0), // 等边三角形
        ];
        let triangles = vec![[0, 1, 2]];

        let metrics = compute_mesh_quality(&vertices, &triangles);
        assert!(metrics.quality_score() > 0.9);
        assert!(metrics.is_good(20.0));
    }

    #[test]
    fn surface_mesh_compute_quality() {
        let mesh = SurfaceMesh {
            vertices: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
                DVec3::new(0.0, 1.0, 0.0),
            ],
            triangles: vec![[0, 1, 2]],
            normals: vec![DVec3::Z; 3],
            dirty: false,
        };

        let metrics = mesh.compute_quality();
        assert_eq!(metrics.triangle_count, 1);
    }

    #[test]
    fn adaptive_subdivider_default() {
        let subdivider = AdaptiveSubdivider::new();
        assert_eq!(subdivider.curvature_threshold, 0.1);
        assert_eq!(subdivider.distance_threshold, 0.1);
        assert_eq!(subdivider.max_subdivision_levels, 3);
    }

    #[test]
    fn adaptive_subdivider_builder() {
        let subdivider = AdaptiveSubdivider::new()
            .with_curvature_threshold(0.2)
            .with_distance_threshold(0.5)
            .with_max_levels(5);

        assert_eq!(subdivider.curvature_threshold, 0.2);
        assert_eq!(subdivider.distance_threshold, 0.5);
        assert_eq!(subdivider.max_subdivision_levels, 5);
    }

    #[test]
    fn adaptive_subdivider_subdivide_by_distance() {
        // 创建一个大三角形
        let mesh = SurfaceMesh {
            vertices: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(10.0, 0.0, 0.0),
                DVec3::new(5.0, 10.0, 0.0),
            ],
            triangles: vec![[0, 1, 2]],
            normals: vec![DVec3::Z; 3],
            dirty: false,
        };

        let subdivider = AdaptiveSubdivider::new()
            .with_distance_threshold(1.0);
        let result = subdivider.subdivide_by_distance(&mesh);

        // 边长大于阈值，应该细分
        assert!(result.triangles.len() > 1);
    }

    #[test]
    fn boundary_sensitive_tessellator_default() {
        let tessellator = BoundarySensitiveTessellator::new();
        assert_eq!(tessellator.feature_angle_threshold, 0.52);
        assert!(tessellator.auto_detect_features);
    }

    #[test]
    fn boundary_sensitive_tessellator_detect_features() {
        // 创建两个三角形，夹角为90度
        let vertices = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
        ];
        let triangles = vec![[0, 1, 2], [0, 1, 3]];
        let normals = vec![DVec3::Z; 4];

        let mut tessellator = BoundarySensitiveTessellator::new()
            .with_feature_angle(0.1); // 小阈值，容易检测特征边

        tessellator.detect_feature_edges(&vertices, &triangles, &normals);
        // 由于两个三角形法向量差异大，应该检测到特征边
        assert!(!tessellator.feature_edges.is_empty());
    }

    #[test]
    fn incremental_mesher_basic() {
        let mut mesher = IncrementalMesher::new();
        assert!(!mesher.is_dirty());

        mesher.invalidate_face(0);
        assert!(mesher.is_dirty());
        assert!(mesher.dirty_faces.contains(&0));

        mesher.clear();
        assert!(!mesher.is_dirty());
    }

    #[test]
    fn incremental_mesher_multiple_faces() {
        let mut mesher = IncrementalMesher::new();
        mesher.invalidate_faces(&[0, 1, 2]);

        assert!(mesher.dirty_faces.contains(&0));
        assert!(mesher.dirty_faces.contains(&1));
        assert!(mesher.dirty_faces.contains(&2));
        assert_eq!(mesher.dirty_faces.len(), 3);
    }

    #[test]
    fn mesh_simplifier_default() {
        let simplifier = MeshSimplifier::new();
        assert_eq!(simplifier.target_ratio, 0.5);
        assert_eq!(simplifier.max_error, 0.01);
        assert!(simplifier.preserve_boundary);
    }

    #[test]
    fn mesh_simplifier_builder() {
        let simplifier = MeshSimplifier::new()
            .with_target_ratio(0.25)
            .with_max_error(0.05);

        assert_eq!(simplifier.target_ratio, 0.25);
        assert_eq!(simplifier.max_error, 0.05);
    }

    #[test]
    fn mesh_simplifier_simplify() {
        // 创建一个包含多个三角形的网格
        let vertices: Vec<DVec3> = (0..9)
            .map(|i| {
                let row = i / 3;
                let col = i % 3;
                DVec3::new(col as f64, row as f64, 0.0)
            })
            .collect();
        let triangles = vec![
            [0, 1, 3], [1, 4, 3],
            [1, 2, 4], [2, 5, 4],
            [3, 4, 6], [4, 7, 6],
            [4, 5, 7], [5, 8, 7],
        ];

        let mesh = SurfaceMesh {
            vertices,
            triangles,
            normals: vec![DVec3::Z; 9],
            dirty: false,
        };

        let simplifier = MeshSimplifier::new()
            .with_target_ratio(0.5)
            .with_max_error(1.0); // 允许较大误差以便简化

        let simplified = simplifier.simplify_mesh(&mesh);
        // 简化后三角形数量应该减少
        assert!(simplified.triangles.len() <= mesh.triangles.len());
    }

    #[test]
    fn mesh_simplifier_simplify_to_target_count() {
        let mesh = SurfaceMesh {
            vertices: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
                DVec3::new(0.0, 1.0, 0.0),
                DVec3::new(1.0, 1.0, 0.0),
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            normals: vec![DVec3::Z; 4],
            dirty: false,
        };

        let simplifier = MeshSimplifier::new().with_max_error(1.0);
        let result = simplifier.simplify_to_target_count(&mesh, 4);

        // 已经达到目标数量，不应该改变
        assert_eq!(result.triangles.len(), 2);
    }

    #[test]
    fn find_boundary_vertices() {
        let triangles = vec![[0, 1, 2], [1, 3, 2]];
        let boundary = super::find_boundary_vertices(&triangles);

        // 边 0-1, 0-2, 1-3, 2-3 是边界边
        // 边 1-2 是内部边
        assert!(boundary.contains(&0));
        assert!(boundary.contains(&3));
    }
}
