//! BVH（包围体层次，Bounding Volume Hierarchy）加速结构。
//!
//! 使用 SAH（Surface Area Heuristic）构建，加速以下查询：
//! - 射线拾取（ray picking）
//! - 形状最小距离（min_distance）
//! - 间隙/重叠检测（detect_gaps_overlaps）
//!
//! 类比 OCCT `BVH_Tree` / `BVH_Builder`。

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::SurfaceEval;

/// 轴对齐包围盒（AABB）。
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    pub min: DVec3,
    pub max: DVec3,
}

impl Aabb {
    /// 空 AABB（min > max，不包含任何点）。
    pub fn empty() -> Self {
        Self {
            min: DVec3::splat(f64::INFINITY),
            max: DVec3::splat(f64::NEG_INFINITY),
        }
    }

    /// 从两个点构建 AABB。
    pub fn from_points(pts: &[DVec3]) -> Self {
        let mut aabb = Self::empty();
        for &p in pts {
            aabb.expand_point(p);
        }
        aabb
    }

    /// 扩展以包含一个点。
    pub fn expand_point(&mut self, p: DVec3) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    /// 扩展以包含另一个 AABB。
    pub fn expand_aabb(&mut self, other: &Aabb) {
        self.min = self.min.min(other.min);
        self.max = self.max.max(other.max);
    }

    /// 返回 AABB 中心点。
    pub fn center(&self) -> DVec3 {
        (self.min + self.max) * 0.5
    }

    /// 返回 AABB 的表面积（用于 SAH）。
    pub fn surface_area(&self) -> f64 {
        let d = self.max - self.min;
        if d.x < 0.0 || d.y < 0.0 || d.z < 0.0 {
            return 0.0; // 空 AABB
        }
        2.0 * (d.x * d.y + d.y * d.z + d.z * d.x)
    }

    /// 检测与另一个 AABB 是否相交。
    pub fn intersects(&self, other: &Aabb) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }

    /// 检测射线与 AABB 是否相交，返回最近交点参数 t（仅正向）。
    pub fn ray_intersect(&self, origin: DVec3, inv_dir: DVec3) -> Option<f64> {
        let t1 = (self.min - origin) * inv_dir;
        let t2 = (self.max - origin) * inv_dir;

        let t_min = t1.min(t2);
        let t_max = t1.max(t2);

        let t_enter = t_min.x.max(t_min.y).max(t_min.z);
        let t_exit = t_max.x.min(t_max.y).min(t_max.z);

        if t_exit >= t_enter.max(0.0) {
            Some(t_enter.max(0.0))
        } else {
            None
        }
    }

    /// 点到 AABB 的最小距离平方。
    pub fn point_dist_sq(&self, p: DVec3) -> f64 {
        let clamped = p.clamp(self.min, self.max);
        (p - clamped).length_squared()
    }

    /// 检测点是否在 AABB 内（含边界）。
    pub fn contains_point(&self, p: DVec3) -> bool {
        p.cmpge(self.min).all() && p.cmple(self.max).all()
    }
}

/// BVH 节点（内部节点或叶节点）。
#[derive(Debug, Clone)]
enum BvhNode {
    /// 叶节点：包含若干面的索引。
    Leaf {
        aabb: Aabb,
        /// 面索引范围（在 `Bvh.face_indices` 中的 [start, end)）。
        start: usize,
        end: usize,
    },
    /// 内部节点：包含左右子节点索引（在 `Bvh.nodes` 中）。
    Internal {
        aabb: Aabb,
        left: usize,
        right: usize,
    },
}

impl BvhNode {
    fn aabb(&self) -> &Aabb {
        match self {
            BvhNode::Leaf { aabb, .. } => aabb,
            BvhNode::Internal { aabb, .. } => aabb,
        }
    }
}

/// BVH 树，绑定到一个 BRep 的面集合。
///
/// 构建后可用于射线拾取、最近面查询等加速操作。
pub struct Bvh {
    /// 节点数组（根节点为 index 0）。
    nodes: Vec<BvhNode>,
    /// 面索引数组（叶节点通过 [start, end) 引用此数组）。
    face_indices: Vec<usize>,
    /// 每个面的 AABB（按原始面索引存储）。
    face_aabbs: Vec<Aabb>,
    /// 每个面的中心点（用于 SAH 排序）。
    face_centers: Vec<DVec3>,
}

/// 每个叶节点最多包含的面数。
const MAX_LEAF_SIZE: usize = 4;

/// SAH 划分候选数（在每个轴上采样的划分位置数）。
const SAH_BUCKETS: usize = 8;

impl Bvh {
    /// 为 BRep 的所有面构建 BVH 树。
    ///
    /// 采样策略：每个面取 boundary 顶点 + 面法向偏移少量采样点，
    /// 保证 AABB 能覆盖整个面（含曲面面片的近似范围）。
    pub fn build(brep: &BRep) -> Self {
        let faces = &brep.solids[0].shells[0].faces;
        let n_faces = faces.len();

        let mut face_aabbs = Vec::with_capacity(n_faces);
        let mut face_centers = Vec::with_capacity(n_faces);

        for (fi, face) in faces.iter().enumerate() {
            let mut aabb = Aabb::empty();

            // 从边界顶点构建 AABB
            for &wire_edge in &face.outer_wire.edges {
                let edge = &brep.edges[wire_edge.idx];
                let v0 = brep.vertices[edge.start].point;
                let v1 = brep.vertices[edge.end].point;
                aabb.expand_point(v0);
                aabb.expand_point(v1);
            }

            // 对于曲面面片，额外采样面内部点扩展 AABB
            if let Some(surf_idx) = brep.geom.face_surface.get(fi).and_then(|s| *s) {
                let surface = &brep.geom.surfaces[surf_idx];
                let domain = surface.default_domain();
                let [u0, u1, v0, v1] = domain;
                // 采样 3x3 网格以覆盖曲面范围
                for i in 0..=2 {
                    for j in 0..=2 {
                        let u = u0 + (u1 - u0) * i as f64 / 2.0;
                        let v = v0 + (v1 - v0) * j as f64 / 2.0;
                        let p = surface.point_at(u, v);
                        if p.is_finite() {
                            aabb.expand_point(p);
                        }
                    }
                }
            }

            // 退化面保护：至少给一个微小体积
            let size = aabb.max - aabb.min;
            if size.x < 1e-10 { aabb.min.x -= 1e-10; aabb.max.x += 1e-10; }
            if size.y < 1e-10 { aabb.min.y -= 1e-10; aabb.max.y += 1e-10; }
            if size.z < 1e-10 { aabb.min.z -= 1e-10; aabb.max.z += 1e-10; }

            let center = aabb.center();
            face_aabbs.push(aabb);
            face_centers.push(center);
        }

        let face_indices: Vec<usize> = (0..n_faces).collect();
        let mut bvh = Bvh {
            nodes: Vec::new(),
            face_indices,
            face_aabbs,
            face_centers,
        };

        if n_faces > 0 {
            bvh.build_recursive(0, n_faces);
        }

        bvh
    }

    /// 递归构建 BVH 节点，返回新节点在 `nodes` 中的索引。
    fn build_recursive(&mut self, start: usize, end: usize) -> usize {
        let count = end - start;

        // 计算当前范围的 AABB
        let mut aabb = Aabb::empty();
        for i in start..end {
            aabb.expand_aabb(&self.face_aabbs[self.face_indices[i]]);
        }

        // 叶节点条件：面数足够少
        if count <= MAX_LEAF_SIZE {
            let node_idx = self.nodes.len();
            self.nodes.push(BvhNode::Leaf { aabb, start, end });
            return node_idx;
        }

        // SAH 选择最佳划分轴和划分位置
        let (split_axis, split_pos) = self.sah_split(start, end, &aabb);

        // 按划分位置重排 face_indices
        let mid = self.partition(start, end, split_axis, split_pos);

        // 防止退化划分（所有面都在同一侧）
        let mid = if mid == start || mid == end {
            (start + end) / 2
        } else {
            mid
        };

        // 占位（先 push 内部节点，再递归）
        let node_idx = self.nodes.len();
        self.nodes.push(BvhNode::Internal { aabb: Aabb::empty(), left: 0, right: 0 });

        let left = self.build_recursive(start, mid);
        let right = self.build_recursive(mid, end);

        // 更新内部节点
        self.nodes[node_idx] = BvhNode::Internal { aabb, left, right };

        node_idx
    }

    /// SAH 选择最佳划分：返回 (轴索引 0/1/2, 划分位置)。
    fn sah_split(&self, start: usize, end: usize, parent_aabb: &Aabb) -> (usize, f64) {
        let parent_sa = parent_aabb.surface_area().max(1e-30);
        let mut best_cost = f64::INFINITY;
        let mut best_axis = 0usize;
        let mut best_pos = 0.0f64;

        for axis in 0..3usize {
            let axis_min = match axis {
                0 => parent_aabb.min.x,
                1 => parent_aabb.min.y,
                _ => parent_aabb.min.z,
            };
            let axis_max = match axis {
                0 => parent_aabb.max.x,
                1 => parent_aabb.max.y,
                _ => parent_aabb.max.z,
            };
            let span = axis_max - axis_min;
            if span < 1e-14 {
                continue;
            }

            for b in 1..SAH_BUCKETS {
                let split = axis_min + span * b as f64 / SAH_BUCKETS as f64;

                let mut left_aabb = Aabb::empty();
                let mut right_aabb = Aabb::empty();
                let mut left_count = 0usize;
                let mut right_count = 0usize;

                for i in start..end {
                    let fi = self.face_indices[i];
                    let center_val = match axis {
                        0 => self.face_centers[fi].x,
                        1 => self.face_centers[fi].y,
                        _ => self.face_centers[fi].z,
                    };
                    if center_val < split {
                        left_aabb.expand_aabb(&self.face_aabbs[fi]);
                        left_count += 1;
                    } else {
                        right_aabb.expand_aabb(&self.face_aabbs[fi]);
                        right_count += 1;
                    }
                }

                if left_count == 0 || right_count == 0 {
                    continue;
                }

                let cost = (left_count as f64 * left_aabb.surface_area()
                    + right_count as f64 * right_aabb.surface_area())
                    / parent_sa;

                if cost < best_cost {
                    best_cost = cost;
                    best_axis = axis;
                    best_pos = split;
                }
            }
        }

        // 如果没有找到有效划分，按最长轴中点划分
        if best_cost.is_infinite() {
            let d = parent_aabb.max - parent_aabb.min;
            best_axis = if d.x >= d.y && d.x >= d.z { 0 } else if d.y >= d.z { 1 } else { 2 };
            best_pos = parent_aabb.center()[best_axis];
        }

        (best_axis, best_pos)
    }

    /// 按轴和位置将 face_indices[start..end] 原地分区，返回分界索引。
    fn partition(&mut self, start: usize, end: usize, axis: usize, split_pos: f64) -> usize {
        let mut mid = start;
        for i in start..end {
            let fi = self.face_indices[i];
            let center_val = match axis {
                0 => self.face_centers[fi].x,
                1 => self.face_centers[fi].y,
                _ => self.face_centers[fi].z,
            };
            if center_val < split_pos {
                self.face_indices.swap(i, mid);
                mid += 1;
            }
        }
        mid
    }

    // ──────────────────────────────────────────────────────────────────────────
    // 查询 API
    // ──────────────────────────────────────────────────────────────────────────

    /// 射线拾取：返回第一个与射线相交的面索引及 t 值。
    ///
    /// `origin`：射线起点；`dir`：射线方向（无需归一化）。
    pub fn ray_cast(&self, origin: DVec3, dir: DVec3) -> Option<(usize, f64)> {
        if self.nodes.is_empty() {
            return None;
        }
        let inv_dir = DVec3::new(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);
        let mut best: Option<(usize, f64)> = None;
        self.ray_cast_node(0, origin, inv_dir, &mut best);
        best
    }

    fn ray_cast_node(
        &self,
        node_idx: usize,
        origin: DVec3,
        inv_dir: DVec3,
        best: &mut Option<(usize, f64)>,
    ) {
        let node = &self.nodes[node_idx];
        let t_aabb = node.aabb().ray_intersect(origin, inv_dir);
        let t_hit = match t_aabb {
            None => return,
            Some(t) => t,
        };

        // 如果 AABB 交点已经比当前最佳更远，剪枝
        if let Some((_, best_t)) = best {
            if t_hit > *best_t {
                return;
            }
        }

        match node {
            BvhNode::Leaf { start, end, .. } => {
                for i in *start..*end {
                    let fi = self.face_indices[i];
                    // 简单用面 AABB 作为粗检（精确射线-面相交留给调用方）
                    if let Some(t) = self.face_aabbs[fi].ray_intersect(origin, inv_dir) {
                        let update = best.map_or(true, |(_, bt)| t < bt);
                        if update {
                            *best = Some((fi, t));
                        }
                    }
                }
            }
            BvhNode::Internal { left, right, .. } => {
                self.ray_cast_node(*left, origin, inv_dir, best);
                self.ray_cast_node(*right, origin, inv_dir, best);
            }
        }
    }

    /// 返回所有与给定 AABB 相交的面索引。
    ///
    /// 用于间隙/重叠检测的候选面对筛选。
    pub fn query_aabb(&self, query: &Aabb) -> Vec<usize> {
        let mut result = Vec::new();
        if !self.nodes.is_empty() {
            self.query_aabb_node(0, query, &mut result);
        }
        result
    }

    fn query_aabb_node(&self, node_idx: usize, query: &Aabb, result: &mut Vec<usize>) {
        let node = &self.nodes[node_idx];
        if !node.aabb().intersects(query) {
            return;
        }
        match node {
            BvhNode::Leaf { start, end, .. } => {
                for i in *start..*end {
                    let fi = self.face_indices[i];
                    if self.face_aabbs[fi].intersects(query) {
                        result.push(fi);
                    }
                }
            }
            BvhNode::Internal { left, right, .. } => {
                self.query_aabb_node(*left, query, result);
                self.query_aabb_node(*right, query, result);
            }
        }
    }

    /// 返回距离给定点最近的 k 个面索引（近似，按 AABB 距离排序）。
    ///
    /// `max_dist`：搜索半径（超出则不返回）。
    pub fn nearest_faces(&self, point: DVec3, max_dist: f64, max_k: usize) -> Vec<(usize, f64)> {
        let mut candidates: Vec<(usize, f64)> = Vec::new();
        if self.nodes.is_empty() {
            return candidates;
        }
        let max_dist_sq = max_dist * max_dist;
        self.nearest_faces_node(0, point, max_dist_sq, &mut candidates);
        candidates.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        candidates.truncate(max_k);
        candidates
    }

    fn nearest_faces_node(
        &self,
        node_idx: usize,
        point: DVec3,
        max_dist_sq: f64,
        result: &mut Vec<(usize, f64)>,
    ) {
        let node = &self.nodes[node_idx];
        let d_sq = node.aabb().point_dist_sq(point);
        if d_sq > max_dist_sq {
            return;
        }
        match node {
            BvhNode::Leaf { start, end, .. } => {
                for i in *start..*end {
                    let fi = self.face_indices[i];
                    let face_d_sq = self.face_aabbs[fi].point_dist_sq(point);
                    if face_d_sq <= max_dist_sq {
                        result.push((fi, face_d_sq.sqrt()));
                    }
                }
            }
            BvhNode::Internal { left, right, .. } => {
                self.nearest_faces_node(*left, point, max_dist_sq, result);
                self.nearest_faces_node(*right, point, max_dist_sq, result);
            }
        }
    }

    /// 返回与另一个 BVH 可能相交的面对候选列表。
    ///
    /// 用于布尔运算前的面对筛选（替代 O(n²) 暴力遍历）。
    pub fn candidate_pairs(bvh_a: &Bvh, bvh_b: &Bvh) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        if bvh_a.nodes.is_empty() || bvh_b.nodes.is_empty() {
            return pairs;
        }
        Self::candidate_pairs_node(bvh_a, 0, bvh_b, 0, &mut pairs);
        pairs
    }

    fn candidate_pairs_node(
        bvh_a: &Bvh,
        node_a: usize,
        bvh_b: &Bvh,
        node_b: usize,
        pairs: &mut Vec<(usize, usize)>,
    ) {
        let na = &bvh_a.nodes[node_a];
        let nb = &bvh_b.nodes[node_b];

        if !na.aabb().intersects(nb.aabb()) {
            return;
        }

        match (na, nb) {
            (BvhNode::Leaf { start: sa, end: ea, .. }, BvhNode::Leaf { start: sb, end: eb, .. }) => {
                for ia in *sa..*ea {
                    for ib in *sb..*eb {
                        let fa = bvh_a.face_indices[ia];
                        let fb = bvh_b.face_indices[ib];
                        if bvh_a.face_aabbs[fa].intersects(&bvh_b.face_aabbs[fb]) {
                            pairs.push((fa, fb));
                        }
                    }
                }
            }
            (BvhNode::Internal { left: la, right: ra, .. }, _) => {
                Self::candidate_pairs_node(bvh_a, *la, bvh_b, node_b, pairs);
                Self::candidate_pairs_node(bvh_a, *ra, bvh_b, node_b, pairs);
            }
            (_, BvhNode::Internal { left: lb, right: rb, .. }) => {
                Self::candidate_pairs_node(bvh_a, node_a, bvh_b, *lb, pairs);
                Self::candidate_pairs_node(bvh_a, node_a, bvh_b, *rb, pairs);
            }
        }
    }

    /// 返回 BVH 的统计信息（用于调试和性能分析）。
    pub fn stats(&self) -> BvhStats {
        let mut stats = BvhStats::default();
        if !self.nodes.is_empty() {
            self.stats_node(0, 0, &mut stats);
        }
        stats
    }

    fn stats_node(&self, node_idx: usize, depth: usize, stats: &mut BvhStats) {
        stats.node_count += 1;
        stats.max_depth = stats.max_depth.max(depth);
        match &self.nodes[node_idx] {
            BvhNode::Leaf { start, end, .. } => {
                stats.leaf_count += 1;
                stats.total_leaf_faces += end - start;
                stats.max_leaf_faces = stats.max_leaf_faces.max(end - start);
            }
            BvhNode::Internal { left, right, .. } => {
                self.stats_node(*left, depth + 1, stats);
                self.stats_node(*right, depth + 1, stats);
            }
        }
    }
}

/// BVH 统计信息。
#[derive(Debug, Default)]
pub struct BvhStats {
    pub node_count: usize,
    pub leaf_count: usize,
    pub max_depth: usize,
    pub total_leaf_faces: usize,
    pub max_leaf_faces: usize,
}

impl BvhStats {
    pub fn avg_leaf_faces(&self) -> f64 {
        if self.leaf_count == 0 {
            0.0
        } else {
            self.total_leaf_faces as f64 / self.leaf_count as f64
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::{BRep, PrimitiveSolid};

    #[test]
    fn bvh_build_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let bvh = Bvh::build(&brep);
        let stats = bvh.stats();
        // 长方体有 6 个面
        assert_eq!(stats.total_leaf_faces, 6);
        assert!(stats.node_count > 0);
    }

    #[test]
    fn bvh_query_aabb_full() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let bvh = Bvh::build(&brep);

        // 查询整个模型范围，应该返回所有 6 个面
        let big_aabb = Aabb {
            min: DVec3::splat(-10.0),
            max: DVec3::splat(10.0),
        };
        let faces = bvh.query_aabb(&big_aabb);
        assert_eq!(faces.len(), 6);
    }

    #[test]
    fn bvh_query_aabb_empty() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let bvh = Bvh::build(&brep);

        // 查询远离模型的位置，应该返回空
        let far_aabb = Aabb {
            min: DVec3::splat(100.0),
            max: DVec3::splat(200.0),
        };
        let faces = bvh.query_aabb(&far_aabb);
        assert!(faces.is_empty());
    }

    #[test]
    fn bvh_candidate_pairs() {
        let box_a = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let box_b = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let bvh_a = Bvh::build(&box_a);
        let bvh_b = Bvh::build(&box_b);

        // 两个完全重叠的长方体：所有面对都应该是候选
        let pairs = Bvh::candidate_pairs(&bvh_a, &bvh_b);
        assert!(!pairs.is_empty());
    }

    #[test]
    fn bvh_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let bvh = Bvh::build(&brep);
        let stats = bvh.stats();
        assert!(stats.total_leaf_faces > 0);
        assert!(stats.node_count > 0);
    }

    #[test]
    fn aabb_surface_area() {
        let aabb = Aabb {
            min: DVec3::ZERO,
            max: DVec3::ONE,
        };
        // 单位立方体表面积 = 6
        assert!((aabb.surface_area() - 6.0).abs() < 1e-10);
    }

    #[test]
    fn aabb_ray_intersect() {
        let aabb = Aabb {
            min: DVec3::ZERO,
            max: DVec3::ONE,
        };
        let inv_dir = DVec3::new(1.0, 1.0, 1.0); // dir = (1,1,1)
        // 从 (-1,-1,-1) 射向 (1,1,1) 方向，应该相交
        let origin = DVec3::splat(-1.0);
        assert!(aabb.ray_intersect(origin, inv_dir).is_some());

        // 从 (-1,-1,-1) 射向 (-1,-1,-1) 方向（背向），不应该相交
        let inv_dir_back = DVec3::new(-1.0, -1.0, -1.0);
        assert!(aabb.ray_intersect(origin, inv_dir_back).is_none());
    }
}
