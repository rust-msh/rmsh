# rmsh ROADMAP — 开源网格产品能力补全计划

> v1 · 2026-10-12 · 状态：charter 生效
> 配套调研：pro 仓 `docs/rmsh-meshing-landscape-20261012.md`（同类产品全景与缺口依据）。

## 0. 产品定位与依赖宪章

**定位**：rmsh 是一个**相对独立的开源 mesh 组件产品**——格式层（meshio 形态）
+ 网格模型 + 可嵌入剖分算法组件（MeshGems 形态）。v1 **不做** GUI、求解器、
正面 hex-dominant 战争。

**依赖宪章（硬约束）**：rmsh 只允许依赖
1. `std`；
2. rmsh 自家 crate（`rmsh-model/geo/io/algo/...`）；
3. **开源 rcad**（单体内核侧 crate：`rcad-kernel` / `rcad-algo` / `rcad-brep` /
   `rcad-modeling` / `rcad-step` / `rcad-mesh`）。

禁止依赖 pro 侧专有 crate（rcad-io / rcad-document / rcad-cam / rcad-pdm …）。
**核心四库（model/geo/io/algo）保持零 rcad 依赖**——rcad 集成一律走依赖倒置：
查询/几何源 trait 定义在 rmsh-algo，适配实现放 rcad 侧（依赖方向永远
rcad → rmsh）。分层裁定（布尔不进 rmsh、曲率/特征查询挂 trait、边界层与
尺寸场框架自长）见调研文档 §C。

## 1. 商业基准线（"达到商业产品"的可验收定义）

对标组件形态的商业底线（MeshGems 组件能力 + Gmsh 开源基线）：

| 域 | 验收能力 |
|---|---|
| I/O | MSH v2/v4 **ASCII + 二进制**；Gmsh 全元素类型号含 2–3 阶高阶单元；Physical groups / 周期链接 / 分区语义完整往返；VTK/VTU 互操作 |
| 2D | 曲率自适应三角化、特征保持（尖角/短边捕捉）、三角→四边重组、2D 边界层 |
| 3D | 水密表面→四面体（鲁棒、多材料界面吻合）、3D 边界层（prism 增长+碰撞收缩）、体积优化器（swap/collapse/sliver） |
| 自适应 | 尺寸场框架（曲率/距离/数学源、合成、限梯度、背景网格）+ a posteriori 循环接口 |
| 质量 | 优化后硬门：min dihedral ≥ 5°、p95 radius-edge ≤ 2.0、sliver 占比 ≤ 0.1%（按场景族钉入 corpus） |
| 工程 | rayon 并行、Python 绑定、CLI、文档、golden corpus + 公开基准跑分 |

## 2. 阶段计划

### P0 — 地基与宪章（进行中）
- [x] msh 读写归属裁定：pro 侧门面退役（pro 仓 `620205eaa`）
- [x] renderer/viewer 退役（`27d0f9f`，Rule 4 删除；py 移出构建留作 P6 种子）
- [x] algo 内嵌测试 STEP 债清偿（`27d0f9f`，改 msh 夹具 + 解析器节点数守卫真缺陷修复）
- [x] 全测试面绿：`cargo test --workspace` 127/0/1（`27d0f9f`）
- 验收：零编译债、零跳测

### P1 — I/O 对齐 meshio 基线
- [x] MSH v2/v4 **二进制**读写（`7d89357`：v2 分块契约 + v4.1 头 ASCII/载荷二进制 + 端序魔数校验，LE-only 大端点名拒）
- [x] 高阶单元（2–3 阶全类型号）读写（`a5dbb6c`：70 id 全表实证钉死——本地 Gmsh 4.14 块头逐族探针，任务书 10 处错误被实证推翻；Unknown 往返 + 计数守卫收紧）
- [x] Physical groups / 周期语义（`a5dbb6c`：$Periodic 四路径全接线，布局按真 Gmsh 夹具钉死；分区（$Partitioned）延后单列）
- 验收：**Gmsh 互操作往返——已达成**（本地 Gmsh 4.14 闭环：rmsh 读真夹具逐字段精确；Gmsh 读 16 个 rmsh 产物零错、计数精确；百万节点规模压测待跑）

### P2 — 几何接口与自适应地基（依赖倒置首卡）
- [x] `GeometrySource` trait（`ce2fcac`：closest_hit + principal_curvatures v1 面；特征尺寸随特征捕捉卡）
- [ ] rcad 侧适配（rmsh-adaptor crate 或单体内核侧）实现 trait + **BRep 细分桥收编**
- [x] 尺寸场框架（`b023d45`：SizeField trait + 均匀/解析/曲率/距离四源 + Min/Max 合成 + 背景格点限梯度——26 邻域真欧氏 GS 传播、无变化早退即 Bellman 不动点；带孔板验收数字全钉：孔缘 factor·r 精确、过渡带 0.8 fp 精确、梯度超调实测最差 0.82 < 1.05 钉容差）
- [x] 曲率场 + 距离场两个几何派生源（消费 GeometrySource）
- 验收：**达成**——带孔板曲率自适应场景字段级 corpus 钉（P3 把场接入 Frontal-Delaunay 即网格级自适应）；裁剪分类器单列独立卡

### P3 — 2D 商业线（✅ 关闭，2026-10-13 质量尾卡）
- [x] Frontal-Delaunay 曲率自适应接尺寸场（P3a 卡：`with_size_field` 场驱动前缘 + corpus 钉）
- [x] 特征保持（`e2d03f8`，见上行）
- [x] 三角→四边重组（Blossom 式）+ quad paving 收口（quad_recombine 11 测 + quad_paving_2d 8 测，3 todo 已清）
- [x] 2D 边界层（boundary_layer_2d 12 测）
- ✅ 验收达成（2026-10-13 质量尾卡）：三角化根因修复（精确谓词 + Lawson 增量构造 + CDT 域分类加固）翻 concave 精确面积测试并顺带关闭边界锚定平单元尾——两条例行 Gmsh 4.14 CLI 实测以 **×1.15 相对线**收口（uniform p95 0.7058 vs Gmsh 0.7032、graded 0.7582 vs 0.6840，条 0.8086/0.7866）；corpus 孔环场离散化处方被探针驳回（7-gon 环 p95 0.758→0.818 反而恶化，5-gon 保留，证据入档）；新回归守卫：谓词模块 5 单测（精确共圆/大偏移）

### P4 — 3D 商业线
- [x] 水密预检与修复提示（`d0a1bbd`：geo/watertight.rs 942 行——六类缺陷全结构化[孔洞环/非流形/绕向/重复/退化/Möller 自交+邻接豁免] + 每类 hint；已知坏 L 形 0.4 钉共面重叠汤）
- [x] 多材料/分区界面吻合 v1（`c1ca476`：多材料分类 + 腔体减除——射线奇偶级联[Grazing 降级 + winding fallback]、闭合/定向/AABB 拒绝全点名、体积簿记面散度 1e-9）
- [x] **多材料 v2 一致性剖分（2026-10-13）**：凸腔体界面真吻合——kernel tet 按腔面平面递归闭式切割（1+3 角 tet + 棱柱 3-tet / 2+2 双楔各 3-tet，随机 200 平面守恒钉），单网格全材料打标签、界面三角形两侧 bit 一致（同边同插值一次计算 + 量化焊点 + 坐标哈希符号扰动一致切割）、体积恒等式转精确（1e-9，跨骑归零）；**实测覆盖率门**（<98% 任意腔 → 诚实回退 v1 跨骑路径，刀刃退化带[腔面掠过 kernel 顶点/边]入档）；两项 v1 限制收紧：面穿刺 enclosure（tri-tri 边贯穿双向 + 节点在面拒绝，U 形非凸外壳夹具）+ 真面不相交 disjointness（任意接触拒绝，取代 touching-boxes 廉查）；一般腔体（非凸/交切）= TetGen 约束边/面翻转恢复路线（H. Si, ACM TOMS 2015，概念级 clean-room），长期卡
  测试：嵌套立方体（双材料精确体积 0.578125/0.421875 + 界面 99% 带钉）、fully-nested 双材料、退化带回退钉、穿刺/接触拒绝、区域体积恒等 1e-9；workspace 350/0/3
- [x] 3D 边界层（`4ac10f8`：prism 增长 + 射线净空碰撞帽 v1——CCW 规范化/面积加权法向角点对角线精确钉/棱柱正体积 3-四面体分解钉/总体积守恒；Gmsh fan/cross 角点延后 P4 邻接）
- [x] 体积优化器补完（`06db6e2`：mesh_optimize 302 → mod 1540 + tests 707/25 测——度量族[最小二面角/radius-edge 纠正 √6/4/sliver/缩放雅可比/长宽比] + 塌缩[链接条件+模拟合并]/2-3 翻转[穿刺有效+元组改进或滑移放宽]/3-2 双星形有效[总体积不变量抓到 union 增长]/黄金分割重定位[Gmsh 式]；box 最小二面角 5.6→20.6、滑叶 0.033→0；p95≤2.0 未达（6.1，角对角纺锤残余——**已裁定，见下行**）+ tet_mesh 潜伏 3-2 无穿刺校验缺陷记卡）
- [x] **纺锤决策收口（P3 质量尾卡，2026-10-13）：按族换门基线**——探针实证 18/74 个 r/e>3 四面体全部锚定在输入表面契约冻结的角对角面长边上（最差 10.17 跨 (0,0,0)-(1,1,0) 面对角线），对 quad 表面有界盒体是结构性残余，union 保持局部操作集（塌缩/2-3/3-2/重定位）原理上不可达；box 族基线钉 6.107±0.05（其余宪章门全过：min dihedral 20.6≥5°、滑叶 0）；多翻转（4-4/5-6 复合）与退火记为长期路径（mesh_optimize 模块文档 "corner-diagonal spindle family" 节，与 HXT 全量并行插入同款诚实入档模式），v1 不排期
- [x] HXT 并行化（rayon）（`819718b`：图着色[6 色 box 实测钉]节点独立 + rayon 并行评估/串行升序确定性提交[任意线程数 bit-identical 钉] + 并行质量度量[to_bits 等价串行] + num_threads 生效 + 线程观测仪器化[强制 4 线程 size==4 钉]；box p95 并行反超串行 13%；全量并行插入[Hilbert/CAS]入档长期路径）
- 验收：fTetWild benchmark 子集成功率 ≥ 95%；质量门全过线；边界层案例零碰撞

### P5 — 高阶与重网格（✅ 关闭，2026-10-13）
- [x] 高阶弯曲单元（模型层 20 个 2-3 阶语义变体 + 规范双向 Gmsh id 映射
  `to_gmsh_type_id`/`from_gmsh_type_id`；io 读侧语义变体往返四 flavor、写侧
  节点数守卫泛化到全部有表类型；algo `high_order.rs` 弯曲单元机构——17 个正则
  族的参考布局全部真 Gmsh 4.14 探针实证（transfinite/extrude 单元素 + 逐角点
  仿射反推，产物 D:/tmp/rmsh-p5），等参映射 = 截断幂重心格点基（单纯形）×
  1D Lagrange 乘积（张量族）× serendip 单项式集数值构造（quad8/12、hex20、
  prism15、tri9——tri9 是不完备单纯形，直接乘积公式不构成基，单项式集 =
  Gmsh 自己的 gmshGenerateMonomialsTriangle serendip 集）；`.valid` = 采样
  有符号雅可比度量（线=弦向投影、面=角点法向投影面积、体=行列式且相对
  角点标架符号——真网格存在整体反向元）> 0；几何误差钉：90° 弧 Line3
  1.078e-2 ≤ 1.5e-2、Line4 2.727e-3 ≤ 4e-3、直线弦 2.929e-1（27 倍差），
  四分之一圆柱 Tri6 同钉；金字塔三型 io 完整、弯曲映射诚实缓议（Bergot 基，
  探针节点序已入模块文档）。**收尾：最后 3 个忽略项清零（2026-10-13）**
  ——两例 Gmsh CLI 验收钉（uniform/graded ×1.15）转默认流程实跑
  [Gmsh 4.14 缺席时诚实失败，不再 #[ignore]；实测复现
  0.7058/0.7032、0.7582/0.6840]，traits.rs 的 rust,ignore 文档示例
  转正为真 doctest；workspace **375/0/0**）附带修复：v4.1 写侧缺 `$Entities`（元素块引用
  未声明实体，Gmsh 元素读侧拒读而节点读侧会自建——真实缺陷，交叉验证探针
  抓到），ASCII/二进制双路补齐 + Gmsh-exact 布局钉 + 四 flavor 零错互操作）
- [x] MMG 式各向异性重网格补完（`0d56204`：Metric3::intersect 同时约化载波基逐方向 max[Gmsh mostaniso 形，对角平均回退被证据驳回] + collapse_edge_3d[度量空间门 + P4a 序有效性门 + 点名拒绝] + coarsen_pass 接线 + 12/12 测）
- 验收：高阶 corpus 场景 + 重网格 metric 收敛

### P6 — 产品化（2026-10-13 排期裁定更新）
- [x] CLI（convert / refine / quality 三动词起步）——早已建成
  （crates/cli，rmsh bin 26 测全绿），勾账补正 2026-10-13
- [x] **Python 绑定**（✅ 2026-10-13，gmsh 迁移路径重建）：用户裁定
  保留旧 py 的 **gmsh API 兼容层定位**（`import gmsh` → `import rmsh`
  即迁移），在当前架构上重建后端——旧 5059 行单文件（gmsh 仿制 +
  旧 BRep/viewer 债）按 Rule 4 整体退役，新后端 = 模块化
  （state/occ/mesh_api/model_api/io_api）+ 真 rcad 内核：
  occ 基元（box/sphere/cylinder/cone/torus/rectangle/disk）→
  rcad-modeling primapi；cut/fuse/intersect → OCCT 对齐的 rcad-algo
  BOP；translate/rotate/dilate/mirror → BRep::apply_transform；
  synchronize → 内核 IncrementalMesh 三角化（catch_unwind 容器）+
  rmsh-bridge 走读（welded Triangle3 + element→entity 记账）；
  generate(3) → CentroidStarMesher3D；getNodes/getElements → gmsh
  数据布局（类型号 = P5 语义变体经 to_gmsh_type_id 原生映射）；
  STEP 读写 → rcad-step；物理组/选项表/logger。诚实省略：gui/plugin
  （宪章禁域）、fragment/addVolume/extrude/revolve/fillet/chamfer/
  setOrder>1/optimize/refine/recombine（NotImplementedError 带原因）。
  e2e：init→建模→布尔（解析体积 1.86596 钉）→变换→同步→体积网格
  →getNodes/getElements→msh 往返→STEP 导出重导入（发现并钉住
  rcad-step 读侧球面裁剪缺陷，见交接缺陷账）→2D 三角→选项/logger。
  workspace 375/0/0（py 0 测 + e2e 脚本，运行方式见 e2e.py 头注）
- [ ] **$Partitioned 分区语义**（用户裁定需要）：对齐 Gmsh MSH 4.1
  的 $PartitionedEntities（numPartitions / ghost / parent /
  partitions 四路读写，真 Gmsh METIS 分区夹具探针实证），
  四 flavor 完整往返。验收 = 分区网格互操作闭环
- [ ] **VTK/VTU 互操作**（归属裁定：放 rmsh）：格式层是 meshio
  形态的本职；宪章依赖方向恒 rcad → rmsh，VTK 读写只需 rmsh-model
  （核心四库零 rcad 依赖不破）；放 rcad-io 则 rmsh 独立产品永远缺
  这块且方向反了。真 VTK cell 类型映射（含 P5 高阶单元——外仓
  vtk.rs 的 Gmsh-id 兜底臂届时退役）。验收 = ParaView 可读 +
  高阶单元类型正确
- [ ] **a posteriori 自适应循环接口**（用户裁定推进）：误差估计 ←
  网格/解反馈 → 尺寸场更新的迭代闭环 trait 面（消费 P2b 尺寸场框架）
- [ ] **golden corpus + 公开基准跑分**（用户裁定推进）：
  按场景族的语料钉 + 可复现的公开跑分
- [ ] 文档站 + 发布准备（crates.io / GitHub Release 流水线）——
  **用户缓议**（2026-10-13：以后合适的时候做；保留排队不排期）
- 验收：绑定 e2e、分区与 VTK/VTU 互操作、基准公开、文档完整、
  发布流水线绿（缓议项除外）

## 3. 里程碑依赖

P0 → P1 → P2 →（P3、P4 可并行，均依赖 P2 的尺寸场）→ P5 → P6。
P6 的 CLI/文档可与 P3/P4 并行滚动。

## 4. 风险与明示不做

- **hex-dominant**：不做（Cubit 级护城河；重组路线的 quad 已覆盖最近需求）。
- **GUI**：不做（组件形态）。
- **解算器**：不做。
- 工程纪律：C 盘常年满——rmsh 构建一律 `CARGO_TARGET_DIR=D:/tmp/rmsh-target`。
- 开源 rcad 内核在快速演进（布尔/位置合成等 wave 波动）：rmsh 消费面只钉
  trait 与格式语义，不钉内核内部行为。
