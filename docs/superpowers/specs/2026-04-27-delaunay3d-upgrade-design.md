# Phase 1: Delaunay3D 核心升级设计

## 目标

将 `Delaunay3D` 从当前的"CentroidStar 种子 + 启发式四面体分割"升级为真正的 **Bowyer-Watson 3D 增量插入 + Shewchuk 细化**，使其能对非凸/复杂闭合曲面产生高质量的约束 Delaunay 四面体网格。

## 现状

当前 `crates/algo/src/delaunay_3d.rs`（~1930 行）的实际工作方式：

1. `mesh_3d()` → `CentroidStarMesher3D` 生成种子网格
2. `refine_bad_tetrahedra()` → 找到最差四面体，在内部插入点（centroid/barycentric/circumcenter blend），1 个 tet 拆成 4 个
3. `optimize_local_face_flips()` → 2-3 / 3-2 / 4-4 局部翻转（**已有生产级实现，直接复用**）

**缺失的能力：**
- 没有空球性质维护，不是真正的 Delaunay
- 没有 Bowyer-Watson 空洞搜索和逐点连接
- 没有约束边界恢复（边/面）
- 没有外心插入细化
- 只通过分割单个坏 tet 来改善质量，效率低、质量控制弱

## 新架构

```
Phase 1: Bowyer-Watson 3D 增量插入
  1a. 从边界点构建超四面体
  1b. 逐点插入：定位 → 空洞搜索 → 删除空洞 → 连接新点
  1c. 移除超四面体及相邻单元
  1d. 输出：点集 Delaunay 四面体化

Phase 2: 约束边界恢复
  2a. 边恢复（flip 优先，失败则 Steiner 点）
  2b. 面恢复（边交换优先，失败则 Steiner 点）

Phase 3: Shewchuk 细化 + Sliver 消除
  3a. 优先队列：每次处理 ρ 最大的四面体
  3b. 插入外心 → 局部 Delaunay 恢复
  3c. 域外外心拒绝 + 面标记
  3d. 复用已有 optimize_local_face_flips 消除 sliver
```

## 文件变更

| 文件 | 操作 | 职责 |
|---|---|---|
| `crates/algo/src/delaunay_3d.rs` | **重写** | 公开 API（`Delaunay3D` struct、`Mesher3D` impl），管线编排 |
| `crates/algo/src/delaunay_core.rs` | **新建** | 内部数据结构（`TetMesh`）、Bowyer-Watson 插入、空洞搜索、邻居管理 |
| `crates/algo/src/boundary_recovery.rs` | **新建** | 边恢复、面恢复、Steiner 点插入 |
| `crates/algo/src/geometry.rs` | **新建** | 从 `delaunay_3d.rs` 提取几何原语（circumsphere、in_sphere、二分角、体积等），消除重复定义 |

`delaunay_3d.rs` 中将保留：
- `Delaunay3D` 公共 struct（API 不变）
- `Mesher3D` trait impl（API 不变）
- 现有细化辅助函数的过渡期兼容包装

## 内部数据结构

```rust
// delaunay_core.rs

struct Tet {
    nodes: [u32; 4],       // 顶点在 TetMesh.nodes 中的索引
    neighbors: [u32; 4],   // 四个面的邻居 tet 索引（面 i 的对顶点是 nodes[i]）
                           // u32::MAX (= 0xFFFF_FFFF) 表示无边界面
}

struct TetMesh {
    nodes: Vec<[f64; 3]>,          // 所有节点坐标
    tets: Vec<Tet>,                 // 所有四面体
    node_to_surface_id: HashMap<u32, u64>,  // 内部节点映射到外部 node_id
}

impl TetMesh {
    fn bounding_box(&self) -> ([f64; 3], [f64; 3]);
    fn find_containing_tet(&self, p: [f64; 3], seed: u32) -> Option<u32>;
    fn collect_cavity(&self, start_tet: u32, p: [f64; 3]) -> Vec<u32>;
    fn cavity_boundary_faces(&self, cavity: &[u32]) -> Vec<[u32; 3]>;
    fn insert_point(&mut self, p: [f64; 3], surface_id: u64);
    fn remove_tets(&mut self, indices: &[u32]);
    fn update_neighbors(&mut self);
}
```

## Phase 1 关键操作

### 1a. 点定位（walk）

从随机种子出发。对于当前 tet 的 4 个面：若 `orient3d(face[0], face[1], face[2], point) < -eps`（point 在面的外侧），穿过该面进入邻居。终止于 point 在所有面内侧的 tet。退化时选第一个正定向面。

### 1b. 空洞搜索（BFS）

从包含点 p 的 tet 出发 BFS：对所有当前 tet 的 face，`in_sphere_test(opposite_node, face[0], face[1], face[2], p) > 0` 则将该 face 的邻居入队。收集到的所有 tet 即为空洞。

### 1c. 连接新点

空洞的外边界面 = 在空洞中出现恰好一次的面。每个外界面 `[a, b, c]` + 新点 p → 新 tet。邻居在更新阶段通过面匹配重建。

### 1d. 超四面体

算包围盒后，以对角线 d 膨胀：在 `±10d` 处放置顶点，确保超四面体包含所有输入点。删除阶段：标记所有包含任一超四面体顶点的 tet。

### 精确性

- 使用 f64 + epsilon 1e-12
- 退化四点（共面/共线）由 Shewchuk 的 `orient3d` 自适应精确谓词处理
- 体积 ≈ 0 的 tet 在细化阶段被消除

## Phase 2 关键操作

### 边恢复

1. 检测：输入 boundary faces 中每条边，检查是否对应 TetMesh 中某 tet 的一条边
2. 对每条缺失边：追踪线段穿过的面序列
3. 对这些面尝试 2-3 / 3-2 flip 消除相交
4. flip 失败 → 在交点上插入 Steiner 点
5. 最大 flip 尝试次数：8

### 面恢复

1. 检测：每个 boundary face 是否在 TetMesh 中存在
2. 对缺失面：收集穿过该面的边
3. 对每条贯穿边：尝试边交换（3 tet → 2 tet）消除
4. 若全部贯穿边消除后面仍未出现 → 在面内部插入 Steiner 点

## Phase 3 关键操作

### 劣质 tet 优先队列

```
std::collections::BinaryHeap<BadTet>
// BadTet = (tet_index, ρ_score)
// 大顶堆 → 始终处理最差 tet
```

每次插入一个外心后，只把受影响的 tet（与新点相邻的那些）重新入队。

### 外心插入

`circumsphere(a,b,c,d) → (center, radius)` 已有生产级实现，直接复用。插入本身复用 Phase 1 的 `insert_point`，保证局部 Delaunay 恢复。

### 域外拒绝

对每个外心：
1. 检查它是否在输入 boundary surface 定义的闭合域内（射线投射法 + 面半空间）
2. 若在域外，标记对应面为 "已拒绝"，不再尝试
3. 若所有相邻面的外心都被拒绝，对该 tet 使用 centroid 替代外心

### 终止条件

```
while queue not empty:
    if 全部 tet 的 ρ < ρ_max: break
    if 迭代次数 > max_iter: break
```

ρ_max 默认 2.0（Shewchuk 保证可终止），std::env::var 调试开关保留。

### 复用

`optimize_local_face_flips`（~1300 行）——已有的 2-3/3-2/4-4 翻转实现——在细化结束后额外运行 2-4 次消除残余 sliver。该函数通过面/边映射驱动，不依赖 Delaunay 性质，因此与新的 Delaunay 管线完全兼容。

## 几何提取

新建 `geometry.rs`，从 `delaunay_3d.rs` 和 `frontal_3d.rs` 中提取以下重复定义：

| 函数 | 来源 | 用途 |
|---|---|---|
| `circumsphere` | delaunay_3d.rs | Bowyer-Watson, 外心 |
| `in_sphere_test` | delaunay_3d.rs | 空洞搜索 |
| `radius_edge_ratio` / `radius_edge_ratio_points` | delaunay_3d.rs | 细化 |
| `dihedral` / `min_dihedral_points` | delaunay_3d.rs, frontal_3d.rs（重复） | 质量评估 |
| `tetra_volume` | delaunay_3d.rs, frontal_3d.rs（重复） | 退化检测 |
| `solve_3x3` | delaunay_3d.rs, frontal_3d.rs（重复） | circumsphere |
| `point_in_tetrahedron` | delaunay_3d.rs, frontal_3d.rs（重复） | 定位 |
| `select_refinement_point` + 级联函数 | delaunay_3d.rs | 可复用 |

删除原有重复定义，统一引用 `geometry.rs`。

## 测试策略

| 层级 | 内容 | 文件 |
|---|---|---|
| 单元测试 | 每个几何原语：regular tet, degenerate tet, 边界条件 | `geometry.rs` |
| 单元测试 | TetMesh 操作：插入、删除、邻居更新、空洞搜索 | `delaunay_core.rs` |
| 集成测试 | Bowyer-Watson 3D：随机点集 → 验证空球性质 | `delaunay_core.rs` |
| 集成测试 | 边界恢复：立方体各面的边和三角形 | `boundary_recovery.rs` |
| 回归测试 | 保留现有全部 25 个测试（API 兼容） | `delaunay_3d.rs` |
| 质量回归 | 立方体/拉伸盒/细长盒 → 品质度量（已有） | `mesher3d_quality_regression.rs` |
| 新质量测试 | **提高阈值**：min_dihedral > 5°, p95_radius_edge < 5.0, sliver_frac < 0.15 | `mesher3d_quality_regression.rs` |

## 风险与对策

| 风险 | 对策 |
|---|---|
| 边界恢复在某些几何体上不收敛 | 限制 flip 尝试次数 + Steiner 点插入，退化到已有 fallback |
| 外心插入导致四面体数爆炸 | Shewchuk 终止理论保证 ρ_max ≥ 2.0 可终止；添加 max_insertions 硬上限 |
| 浮点精度导致空球判断不一致 | 引入符号扰动（symbolic perturbation）处理共面/共球退化；epsilon 可配置 |
| 现有 API 用户受影响 | `Delaunay3D` struct 和 `Mesher3D` impl 的公开签名不变；新增字段通过 `Default` 兼容 |
| 性能退化（新管线可能更慢） | 步行点定位 O(n½) vs 原 O(1) centroid 选择；对典型输入（10³~10⁵ 点）可接受 |

## 依赖与串行化

本 Phase 不依赖其他算法修改。它提升 `Delaunay3D` 的质量，使得依赖它的 `Frontal3D`、`Hxt3D`、`MmgRemesh` 自然获益（它们都调用 `Delaunay3D`）。

在后续 Phase 中：
- **Frontal3D**：用真正的推进前沿 + Bowyer-Watson 局部插入替代当前"只调 Delaunay3D + 后处理"
- **Hxt3D**：在 Delaunay3D 基础上增加 Hilbert 排序 + Rayon 并行（独立变换）
- **MmgRemesh**：在 Delaunay3D 的网格上叠加度量场驱动的局部操作
