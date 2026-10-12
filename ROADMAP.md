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

### P3 — 2D 商业线
- [ ] Frontal-Delaunay 曲率自适应接尺寸场
- [x] 特征保持（`e2d03f8`，见上行）
- [ ] 三角→四边重组（Blossom 式）+ quad paving 收口（清 3 todo）
- [ ] 2D 边界层
- ⏳ 验收大幅收窄（P3 质量卡 `9c503ff`，Gmsh 2D 管线移植）：均匀 p95 6.37→**1.05**（bar 0.81）、结构病治愈（三角形数/边均值/主体分布对齐 Gmsh）；剩 ~5% 边界锚定平单元尾（需精确谓词，`[p3a-gap]` 钉新基线 1.15/1.35）；梯度盘 1.25（旧 0.98 的入档回归——corpus 孔环应场离散化，随尾卡处理）

### P4 — 3D 商业线
- [ ] 水密预检与修复提示（孔洞/自交分类——geo/topo 扩展）
- [ ] 多材料/分区界面吻合
- [ ] 3D 边界层（prism 增长 + 碰撞收缩）
- [ ] 体积优化器补完（清 5 todo + 质量 gate 联动）
- [ ] HXT 并行化（rayon）
- 验收：fTetWild benchmark 子集成功率 ≥ 95%；质量门全过线；边界层案例零碰撞

### P5 — 高阶与重网格
- [ ] 高阶弯曲单元（valid GE 映射，几何误差 ≤ 容差）
- [ ] MMG 式各向异性重网格补完（清 2 todo）
- 验收：高阶 corpus 场景 + 重网格 metric 收敛

### P6 — 产品化
- [ ] Python 绑定（py crate 移植重写）
- [ ] CLI（convert / refine / quality 三动词起步）
- [ ] 文档站 + 发布准备（crates.io / GitHub Release 流水线）
- 验收：绑定 e2e、文档完整、发布流水线绿

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
