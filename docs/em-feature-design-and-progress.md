# EMStudio 功能设计与开发进度

## 1. 项目概述

EMStudio 是一个类似 Ansys Electronics Desktop 的专业电磁仿真工具，功能对标 **HFSS（全波求解）+ Q3D Extractor（准静态寄生参数提取）**。后端依赖 [Rem](https://github.com/javagg/rem2.git) 完成电磁仿真计算。

### 1.1 技术栈

| 层次 | 技术选型 | 说明 |
|------|---------|------|
| 语言 | Rust (Edition 2024) | 高性能、内存安全 |
| GUI 框架 | egui 0.33 + eframe 0.33 | 即时模式 GUI，跨平台 |
| 3D 渲染 | wgpu 27 | Vulkan/Metal/DX12/WebGPU 统一抽象 |
| 几何内核 | rcad (git submodule) | B-Rep 内核：体素创建、布尔运算、扫掠、STEP 导入导出 |
| 数学库 | glam 0.29 | 向量/矩阵/四元数 |
| 序列化 | serde + serde_json | JSON 工程文件 |
| 2D 绘图 | egui_plot 0.33 | S 参数/RLCG 曲线 |
| 布局系统 | egui_dock 0.18 + egui_tiles 0.14 | 可停靠面板 + 分裂视图 |
| 运行平台 | Native + WASM (trunk) | 桌面应用 + 浏览器 |
| 浏览器存储 | OPFS (Origin Private File System) | Web 基本版 Local-First 模式存储层 |

### 1.2 架构总览

```
┌──────────────────────────────────────────────────────────────────┐
│                     emstudio-main (入口)                         │
│               Native (eframe) / WASM (trunk+WebGPU)              │
├──────────────────────────────────────────────────────────────────┤
│                     emstudio-app (应用层)                         │
│    EmStudioApp · Ribbon Bar · Tab 页 · Dock 面板 · Edition 门控  │
├──────────────┬───────────────────┬───────────────────────────────┤
│ emstudio-    │  emstudio-render  │  emstudio-infra               │
│ components   │  (3D 渲染引擎)    │  (后端抽象)                    │
│ (UI 组件)    │  wgpu Pipeline    │  Standalone / Cloud / OPFS    │
│ Ribbon/Dock  │  Camera/Colormap  │  Solver 调度 · Web Worker     │
├──────────────┴───────────────────┴───────────────────────────────┤
│                    emstudio-domain (领域模型)                      │
│      Project · Design · Material · Geometry · SolutionType       │
├──────────────────────────┬───────────────────────────────────────┤
│  emstudio-solver         │  emstudio-touchstone                  │
│  (求解器 trait + 调度)    │  (Touchstone .snp 解析/写入)          │
│  Rem 集成入口             │  v1.0 & v2.0 全格式支持               │
│  Native + WASM Worker    │                                       │
└──────────────────────────┴───────────────────────────────────────┘
```

### 1.3 Workspace 结构

```
emstudio/
├── Cargo.toml                  # Workspace 根
├── readme.md                   # 项目简介
├── vendor/
│   └── rcad/                   # rcad 几何内核（git submodule）
│       └── libs/
│           ├── rcad-kernel/    # B-Rep 拓扑 + 解析几何
│           ├── rcad-modeling/  # 体素创建 + 扫掠
│           ├── rcad-algorithms/# 布尔运算 + 倒角
│           ├── rcad-render/    # wgpu 细分 + 拾取
│           ├── rcad-step/      # STEP 导入导出
│           └── rcad-scene/     # 场景交互
├── crates/
│   ├── main/                   # 应用入口（Native + WASM）
│   ├── app/                    # 主 UI 应用
│   ├── components/             # 自定义 UI 组件（Ribbon、Dock）
│   ├── domain/                 # 核心领域模型（20+ 子模块）
│   ├── infra/                  # 后端抽象（Standalone / Cloud / WasmBackend）
│   ├── render/                 # GPU 3D 渲染引擎
│   ├── solver/                 # 求解器抽象 + Rem 集成
│   ├── touchstone/             # Touchstone S 参数文件解析
│   └── worker/                 # Web Worker（OPFS + 求解器调度）
├── examples/
│   └── field-vis/              # 3D 场可视化示例
└── docs/
    ├── em-project-file-design.md         # 工程文件格式设计
    ├── em-result-file-formats.md         # 仿真结果文件格式设计
    ├── em-result-visualization-design.md # 结果可视化方案设计
    └── em-feature-design-and-progress.md # 本文档
```

### 1.4 部署方式

EMStudio 支持两种部署方式，功能和行为有所差异：

| 维度 | 桌面版 (Native) | Web 版 (WASM) |
|------|-----------------|---------------|
| 运行环境 | eframe / Vulkan / Metal / DX12 | trunk + WebGPU |
| 工程管理 | 新建 / 打开 / 保存 / 另存为 | 工程预创建、默认打开，仅支持**保存**（不可新建、不可另存） |
| 文件系统 | 本地磁盘 | 服务端存储（专业版/企业版）或 OPFS（基本版 Local-First） |
| 后端通信 | 本地进程调用 | REST API（专业版/企业版）或 浏览器内 Web Worker（基本版 Local-First） |

### 1.5 功能版本（Edition）

EMStudio 分为三个功能版本，桌面版与 Web 版均适用：

| 版本 | 目标用户 | 求解能力 | 部署方式 |
|------|---------|---------|---------|
| **基本版 (Basic)** | 学生、教育、轻量评估 | 基础 HFSS 驱动模态 + Q3D 电容 | 桌面 / Web（含 Local-First） |
| **专业版 (Professional)** | 工程师、中小团队 | 全部 HFSS + Q3D 求解类型 | 桌面 / Web |
| **企业版 (Enterprise)** | 大型团队、HPC 集群 | 全部求解 + 分布式并行 + Optimetrics | 桌面 / Web |

#### 1.5.1 Web 基本版 Local-First 模式

Web 基本版提供完全在浏览器内运行的 **Local-First** 模式，无需服务端：

```
┌─────────────────────────────────────────────────────┐
│                  浏览器主线程                          │
│          egui (WASM) + WebGPU 渲染                    │
├─────────────────────────────────────────────────────┤
│              Web Worker（后端线程）                     │
│   ┌─────────────┐  ┌──────────────────────────────┐ │
│   │ Rem Solver   │  │ OPFS (Origin Private FS)     │ │
│   │ (WASM 编译)  │  │  ├── project.emsp            │ │
│   │              │  │  ├── results/                 │ │
│   │              │  │  └── materials.emsm           │ │
│   └─────────────┘  └──────────────────────────────┘ │
└─────────────────────────────────────────────────────┘
```

| 特性 | 说明 |
|------|------|
| 存储层 | OPFS (Origin Private File System)，运行在 Web Worker 内，不依赖用户可见文件系统 |
| 求解器 | Rem 编译为 WASM，在 Web Worker 中执行，不阻塞 UI 线程 |
| 离线能力 | Service Worker 缓存后可完全离线运行 |
| 数据边界 | 所有数据留在浏览器本地，不上传服务端 |
| 限制 | 受浏览器内存/线程限制，适合中小规模模型 |

### 1.6 工程管理行为差异

| 操作 | 桌面版 | Web 版 |
|------|-------|--------|
| 新建工程 | ✅ 支持 | ❌ 不支持（工程由平台预创建） |
| 打开工程 | ✅ 文件对话框选择 | ✅ 自动打开预分配工程 |
| 保存 | ✅ 支持 | ✅ 支持（写回服务端 / OPFS） |
| 另存为 | ✅ 支持 | ❌ 不支持 |
| 自动保存 | ✅ .emsp.auto | ✅ 定时写入 OPFS / 服务端 |
| 工程锁 | ✅ .emsp.lock（文件锁） | ✅ 服务端锁 / OPFS 独占 |

---

## 2. 求解器能力覆盖

EMStudio 通过后端 Rem 求解器同时支持两大类电磁仿真：

### 2.1 HFSS 全波求解（FEM）

| 求解类型 | 枚举值 | 说明 | 典型应用 |
|---------|--------|------|---------|
| 驱动模态 | `DrivenModal` | 以模态分解计算 S 参数 | 天线、滤波器、连接器 |
| 驱动端子 | `DrivenTerminal` | 以端子电压/电流计算 S 参数 | 多导体传输线 |
| 本征模 | `Eigenmode` | 计算谐振腔本征频率和 Q 值 | 谐振腔、介质谐振器 |
| 瞬态 | `Transient` | 时域有限元分析 | 宽带脉冲响应 |
| SBR+ | `SBRPlus` | 射线追踪 + 物理光学 | 大型散射体 RCS |

**核心输出**：S/Y/Z 参数矩阵、远场方向图、近场分布、天线参数（增益/效率/波束宽度）、E/H/J 场分布。

### 2.2 Q3D 准静态寄生参数提取（MoM + FMM）

| 求解类型 | 枚举值 | 说明 | 典型应用 |
|---------|--------|------|---------|
| DC 电阻 + 电感 | `Q3D_DCRL` | 直流电阻和低频电感 | PCB 走线 DC IR-Drop |
| AC 电阻 + 电感 | `Q3D_ACRL` | 含趋肤/邻近效应 | 高速互连、封装 |
| 电容 | `Q3D_C` | 拉普拉斯静电求解 | 信号间耦合电容 |
| 电容 + 电导 | `Q3D_CG` | 含介质损耗 | 有损介质中的电容 |

**核心输出**：RLCG 矩阵（频率依赖）、DC 矩阵、等效电路 SPICE 网表、电流密度/电荷分布场数据。

---

## 3. 功能模块设计

### 3.1 工程管理

| 功能 | 说明 | 详细设计 |
|------|------|---------|
| 工程文件 (.emsp) | JSON 格式，Project → Design 多层结构 | [em-project-file-design.md §1-2](em-project-file-design.md) |
| 多设计支持 | 一个工程包含多个独立设计 | [em-project-file-design.md §2](em-project-file-design.md) |
| 参数化变量 | 工程级/设计级变量，表达式求值 | [em-project-file-design.md §3](em-project-file-design.md) |
| 数据集 | 频率/温度依赖查找表 | [em-project-file-design.md §3](em-project-file-design.md) |
| 自动保存/恢复 | .emsp.auto 备份 + 崩溃恢复 | [em-project-file-design.md §1.2](em-project-file-design.md) |
| 工程锁 | .emsp.lock 防止并发写入 | [em-project-file-design.md §1.2](em-project-file-design.md) |
| 版本迁移 | 向前兼容的 JSON 格式升级策略 | [em-project-file-design.md §10](em-project-file-design.md) |

### 3.2 几何建模

| 功能 | 说明 | 详细设计 |
|------|------|---------|
| 历史记录式建模 | 操作历史 + 对象快照双层结构 | [em-project-file-design.md §3.3](em-project-file-design.md) |
| 基本体素 | Box、Cylinder、Sphere、Cone、Torus | [em-project-file-design.md §3.3.2](em-project-file-design.md) |
| 布尔运算 | Unite、Subtract、Intersect | [em-project-file-design.md §3.3.2](em-project-file-design.md) |
| 几何变换 | 平移、旋转、缩放、镜像 | [em-project-file-design.md §3.3.2](em-project-file-design.md) |
| 扫掠操作 | SweepAlongPath、SweepAroundAxis | [em-project-file-design.md §3.3.2](em-project-file-design.md) |
| 坐标系 | 全局/局部坐标系定义 | [em-project-file-design.md §3.2](em-project-file-design.md) |
| 命名选择 | 面/边/顶点/对象的命名集合 | [em-project-file-design.md §3.2](em-project-file-design.md) |
| 参数化尺寸 | 几何尺寸可引用变量表达式 | [em-project-file-design.md §3.1](em-project-file-design.md) |

### 3.3 材料系统

| 功能 | 说明 | 详细设计 |
|------|------|---------|
| 材料属性 | ε_r, μ_r, σ, tan_δ, 密度 | [em-project-file-design.md §4.3](em-project-file-design.md) |
| 属性值来源 | 常量 / 表达式 / 数据集 三种模式 | [em-project-file-design.md §4.3](em-project-file-design.md) |
| 材料分类 | Conductor / Dielectric / Magnetic / Composite / Gas | [em-project-file-design.md §4.3](em-project-file-design.md) |
| 材料库 | .emsm 可复用材料库文件 | [em-project-file-design.md §1.2](em-project-file-design.md) |
| 频率依赖 | 通过 Dataset 引用实现频率依赖材料属性 | [em-project-file-design.md §3.1](em-project-file-design.md) |

### 3.4 边界条件

| 边界类型 | 适用 | 说明 |
|---------|------|------|
| PerfectE (PEC) | HFSS / Q3D | 理想电导体 |
| PerfectH (PMC) | HFSS / Q3D | 理想磁导体 |
| Radiation | HFSS | 辐射吸收边界 |
| PML | HFSS | 完美匹配层 |
| Impedance | HFSS | 阻抗边界 |
| FiniteConductivity | HFSS / Q3D | 有限电导率 |
| Symmetry | HFSS / Q3D | 对称面 |
| MasterSlave | HFSS | 周期边界 |
| ThinConductor | Q3D | 薄导体（PCB 走线） |
| InfiniteGroundPlane | Q3D | 无限大地平面 |
| OpenBoundary | Q3D | 开放边界 |

> 详细设计：[em-project-file-design.md §4.5](em-project-file-design.md)

### 3.5 激励/端口

| 激励类型 | 适用 | 说明 |
|---------|------|------|
| WavePort | HFSS | 波端口（波导/同轴截面） |
| LumpedPort | HFSS | 集总端口（PCB 走线） |
| FloquetPort | HFSS | Floquet 端口（周期结构） |
| IncidentWave | HFSS | 入射波（散射/RCS） |
| VoltageDrop | HFSS | 电压差激励 |
| Source | Q3D | 源端子（电流注入点） |
| Sink | Q3D | 汇端子（电流回路点） |

> 详细设计：[em-project-file-design.md §4.6](em-project-file-design.md)

### 3.6 Q3D 网络定义

| 功能 | 说明 |
|------|------|
| Net 定义 | 将导体对象分组为命名电气网络 |
| 源-汇端子 | 每个 Net 上分配 Source/Sink 端子对 |
| 接地网络 | 标记参考地网络 (is_ground_reference) |
| RLCG 矩阵索引 | 所有矩阵按 Net/Terminal 索引 |

> 详细设计：[em-project-file-design.md §4.6.1](em-project-file-design.md)

### 3.7 分析设置

**HFSS 分析设置**：

| 参数 | 说明 |
|------|------|
| 求解频率 | 自适应加密的参考频率 |
| 最大轮次 / Max Delta S | 自适应收敛控制 |
| 频率扫描 | Discrete / Interpolating / Fast |
| 网格初始化 | Lambda 目标、默认网格密度 |

**Q3D 分析设置**：

| 参数 | 说明 |
|------|------|
| 求解类型 | Q3D_DCRL / Q3D_ACRL / Q3D_C / Q3D_CG |
| 自适应频率 | MoM 自适应加密参考频率 |
| Max Delta Energy | 能量收敛阈值 |
| DC 设置 | DC 电阻/电感提取开关 |
| 频率扫描 | 对数/线性离散扫描 |

> 详细设计：[em-project-file-design.md §4.8](em-project-file-design.md)

### 3.8 网格控制

| 类型 | 说明 |
|------|------|
| LengthBased | 最大单元尺寸约束 |
| SkinDepth | 趋肤深度层数控制 |
| CurvatureBased | 曲面法线偏差控制 |
| ModelResolution | 全局最小特征尺寸 |

> 详细设计：[em-project-file-design.md §4.7](em-project-file-design.md)

### 3.9 后处理报告

**报告类别**：

| 类别 | 适用 | 可用数据 |
|------|------|---------|
| SParameter | HFSS | S/Y/Z 参数、VSWR、群时延 |
| FarField | HFSS | 增益、方向性、轴比、效率 |
| NearField | HFSS | 近场 E/H 分布 |
| Fields | HFSS / Q3D | 场量表达式 |
| Eigenmode | HFSS | 谐振频率、Q 值 |
| Emission | HFSS | 辐射功率、EMI |
| RLCGMatrix | Q3D | RLCG 矩阵元素 vs 频率 |
| Q3DFields | Q3D | 电流密度、电场、电荷分布 |

**图表类型**：Rectangular | Polar | Smith | DataTable | Polar3D | MatrixTable

> 详细设计：[em-project-file-design.md §4.14](em-project-file-design.md)

### 3.10 输出变量

从仿真结果派生的命名数学表达式，可用于 Reports 和 Optimetrics。

| 适用 | 示例 |
|------|------|
| HFSS | `dB(S(Port1,Port1))`、`max(GainTotal)`、`bandwidth(S(1,1), -10)` |
| Q3D | `R(Signal1:T1, Signal1:T1)`、`L(DiffP:TP_src, DiffN:TN_src)`、`C(...)` |

> 详细设计：[em-project-file-design.md §4.11](em-project-file-design.md)

### 3.11 场叠加显示

| 场量 | 适用 | 说明 |
|------|------|------|
| E（电场） | HFSS / Q3D | 电场强度分布 |
| H（磁场） | HFSS | 磁场强度分布 |
| Jvol（体电流） | HFSS / Q3D | 体电流密度 |
| Jsurf（面电流） | HFSS / Q3D | 面电流密度 |
| Poynting | HFSS | 坡印廷矢量（功率流） |
| SAR | HFSS | 比吸收率 |
| ChargeDistribution | Q3D | 电荷面密度 |
| OhmicLoss | Q3D | 欧姆损耗密度 |

**绘图类型**：Surface | CutPlane | Volume | Line
**显示模式**：Shaded（云图）| Arrow（矢量箭头）| 相位动画

> 详细设计：[em-project-file-design.md §4.12](em-project-file-design.md)

### 3.12 Optimetrics（参数化分析与优化）

| 类型 | 说明 |
|------|------|
| ParametricSweep | 参数扫描（LinearStep / DiscreteList / LogScale） |
| Optimization | 优化（QuasiNewton / PatternSearch / GeneticAlgorithm / SNLP） |
| Sensitivity | 灵敏度分析 |
| Statistical | 统计分析（蒙特卡罗） |
| Tuning | 实时交互调参 |

> 详细设计：[em-project-file-design.md §4.13](em-project-file-design.md)

### 3.13 辐射设置（HFSS 专用）

| 类型 | 说明 |
|------|------|
| InfiniteSphere | 远场球面采样（theta-phi 方向图） |
| InfinitePlane | 远场平面采样 |
| Line | 沿直线采样近场 |
| Rectangle | 矩形平面采样近场 |
| Sphere | 球面采样近场（近远场变换） |

> 详细设计：[em-project-file-design.md §4.10](em-project-file-design.md)

---

## 4. 仿真结果体系

### 4.1 结果文件格式

| 文件 | 格式 | 适用 | 详细设计 |
|------|------|------|---------|
| validation_report.json | JSON | HFSS / Q3D | [em-result-file-formats.md §2.1](em-result-file-formats.md) |
| convergence.json | JSON | HFSS / Q3D | [em-result-file-formats.md §2.2](em-result-file-formats.md) |
| mesh_stats.json | JSON | HFSS / Q3D | [em-result-file-formats.md §2.3](em-result-file-formats.md) |
| profile.json | JSON | HFSS / Q3D | [em-result-file-formats.md §2.4](em-result-file-formats.md) |
| s_parameters.json | JSON | HFSS | [em-result-file-formats.md §2.5](em-result-file-formats.md) |
| s_parameters.snp | Touchstone | HFSS | [em-result-file-formats.md §2.6](em-result-file-formats.md) |
| far_field_*.json | JSON | HFSS | [em-result-file-formats.md §2.7](em-result-file-formats.md) |
| near_field_*.json | JSON | HFSS | [em-result-file-formats.md §2.8](em-result-file-formats.md) |
| rlcg_matrix.json | JSON | Q3D | [em-result-file-formats.md §2.9](em-result-file-formats.md) |
| equivalent_circuit.sp | SPICE | Q3D | [em-result-file-formats.md §2.10](em-result-file-formats.md) |
| summary.json | JSON | HFSS / Q3D | [em-result-file-formats.md §2.11](em-result-file-formats.md) |
| solver.log | Text | HFSS / Q3D | [em-result-file-formats.md §2.12](em-result-file-formats.md) |
| *.msh | Gmsh MSH 4.1 | HFSS / Q3D | [em-result-file-formats.md §3.2](em-result-file-formats.md) |
| *.emsfld | Binary | HFSS / Q3D | [em-result-file-formats.md §3.3](em-result-file-formats.md) |

### 4.2 结果目录结构

```
Project.emsp.results/
├── design-001/                          # HFSS 全波设计
│   ├── validation_report.json
│   ├── solver.log
│   ├── Setup1/
│   │   ├── convergence.json
│   │   ├── mesh_stats.json, profile.json, s_parameters.json
│   │   ├── mesh/                        # 各轮网格 (.msh)
│   │   ├── solutions/                   # 场解数据 (.bin)
│   │   ├── fields/                      # 导出场数据 (.emsfld)
│   │   ├── far_field/, near_field/      # 远/近场结果
│   │   └── Sweep1/                      # 频率扫描结果
│   ├── Setup1__Optimetrics/             # 参数化/优化结果
│   └── exports/                         # 用户导出 (CSV/PNG)
│
└── design-002/                          # Q3D 准静态设计
    ├── validation_report.json
    ├── solver.log
    ├── Q3D_Setup1/
    │   ├── convergence.json             # delta_energy + rlcg_snapshot
    │   ├── mesh_stats.json, profile.json
    │   ├── rlcg_matrix.json             # RLCG 矩阵
    │   ├── mesh/                        # MoM 面网格
    │   ├── fields/                      # j_field/e_field/charge
    │   └── Sweep1/                      # 频率扫描 RLCG
    │       ├── rlcg_matrix.json
    │       └── s_parameters_from_rlcg.snp
    └── exports/
        ├── equivalent_circuit.sp        # SPICE 网表
        └── report_*.csv
```

> 详细设计：[em-project-file-design.md §5](em-project-file-design.md)、[em-result-file-formats.md §4](em-result-file-formats.md)

---

## 5. 可视化系统

### 5.1 架构

```
┌─────────────────────────────────────────────────────────┐
│                    UI Layer (egui)                        │
│  ReportPanel · FieldPanel · FarFieldPanel · Properties   │
├─────────────────────────────────────────────────────────┤
│              Visualization Mapping Layer                  │
│  QuantityExpr · ColorMapper · GeometryMapper             │
├─────────────────────────────────────────────────────────┤
│                 Data Access Layer                         │
│  JsonLoader · MshLoader · FldLoader · SnpLoader          │
├─────────────────────────────────────────────────────────┤
│               GPU Render Layer (wgpu)                     │
│  MeshPipeline · FieldPipeline · FarFieldPipeline         │
│  Camera · LightSystem · PickingSystem                    │
└─────────────────────────────────────────────────────────┘
```

### 5.2 可视化类型

| 类别 | 类型 | 渲染方式 | 详细设计 |
|------|------|---------|---------|
| 2D 报告 | S 参数矩形图 | egui_plot | [em-result-visualization-design.md §3.5](em-result-visualization-design.md) |
| 2D 报告 | Smith 圆图 | egui_plot | [em-result-visualization-design.md §3.6](em-result-visualization-design.md) |
| 2D 报告 | 极坐标方向图 | egui_plot | [em-result-visualization-design.md §3.7](em-result-visualization-design.md) |
| 2D 报告 | 收敛曲线 | egui_plot | [em-result-visualization-design.md §3.8](em-result-visualization-design.md) |
| 2D 报告 | RLCG vs 频率 (Q3D) | egui_plot | [em-result-visualization-design.md §3.9](em-result-visualization-design.md) |
| 数据表 | RLCG 矩阵表格 (Q3D) | egui Table | [em-result-visualization-design.md §3.10](em-result-visualization-design.md) |
| 2D 报告 | Q3D 收敛曲线 | egui_plot | [em-result-visualization-design.md §3.11](em-result-visualization-design.md) |
| 3D 场图 | 表面云图 | wgpu | [em-result-visualization-design.md §4.2](em-result-visualization-design.md) |
| 3D 场图 | 矢量箭头 | wgpu instanced | [em-result-visualization-design.md §4.2](em-result-visualization-design.md) |
| 3D 场图 | 切面可视化 | wgpu | [em-result-visualization-design.md §4.6](em-result-visualization-design.md) |
| 3D 场图 | 等值面 | wgpu (Marching Tet) | [em-result-visualization-design.md §4.2](em-result-visualization-design.md) |
| 3D 场图 | 相位动画 | wgpu 实时 | [em-result-visualization-design.md §4.7](em-result-visualization-design.md) |
| 3D 远场 | 3D 辐射方向图 | wgpu | [em-result-visualization-design.md §5](em-result-visualization-design.md) |
| 3D 网格 | 网格质量可视化 | wgpu | [em-result-visualization-design.md §6](em-result-visualization-design.md) |
| 交互 | GPU Picking + 探针 | wgpu | [em-result-visualization-design.md §7](em-result-visualization-design.md) |
| 导出 | CSV / PNG / VTK / Touchstone | 文件 I/O | [em-result-visualization-design.md §10](em-result-visualization-design.md) |

---

## 6. 代码模块现状

### 6.1 各 Crate 实现状态

| Crate | 代码量 | 成熟度 | 说明 |
|-------|-------|--------|------|
| **emstudio-domain** | ~2,900 LOC | ⭐⭐⭐⭐ | 完整领域模型（20 个子模块）、表达式引擎、验证、依赖图、文件 I/O |
| **emstudio-render** | ~2,490 LOC | ⭐⭐⭐⭐⭐ | 完整的 wgpu 3D 渲染引擎，5 种可视化模式均已实现 |
| **emstudio-touchstone** | ~1,131 LOC | ⭐⭐⭐⭐⭐ | 生产就绪，Touchstone v1.0/v2.0 完整解析和写入 |
| **emstudio-app** | ~545 LOC | ⭐⭐⭐ | 主 UI 框架搭建完成，Ribbon + Tab + Dock 布局，WASM 后端轮询 |
| **emstudio-worker** | ~260 LOC | ⭐⭐⭐ | Web Worker OPFS 后端 + 求解器调度（Local-First 模式） |
| **emstudio-infra** | ~400 LOC | ⭐⭐⭐ | Standalone + Cloud + WasmBackend（OPFS），Backend trait 含 poll/async solve |
| **emstudio-main** | ~87 LOC | ⭐⭐⭐ | Native 入口完整，WASM 入口使用 LocalFirst 模式 |
| **emstudio-components** | ~81 LOC | ⭐⭐⭐ | Ribbon Bar 和 Dock 面板组件 |
| **emstudio-solver** | ~21 LOC | ⭐ | 仅 trait 定义 + PlaceholderSolver，Rem 未集成 |

**总计**：~7,900+ LOC Rust 代码

### 6.2 Render 引擎详细状态

渲染引擎是目前最成熟的模块，占代码库 59%：

| 子模块 | 代码量 | 状态 | 功能 |
|--------|-------|------|------|
| FieldPipeline | ~631 LOC | ✅ 完成 | wgpu 离屏渲染、WGSL 着色器、色标纹理 |
| FieldSceneState | ~576 LOC | ✅ 完成 | 场景管理、5 种可视化模式、相机控制 |
| FieldMesh | ~323 LOC | ✅ 完成 | 网格生成（UV 球体）、场数据存储 |
| Colormap | ~131 LOC | ✅ 完成 | 4 种色标（Rainbow/Viridis/CoolWarm/Grayscale）|
| OrbitCamera | ~113 LOC | ✅ 完成 | 球坐标轨道相机、旋转/缩放/平移 |
| ArrowPipeline | ~LOC | ✅ 完成 | Instanced 矢量箭头渲染 |
| SliceExtraction | ~LOC | ✅ 完成 | 切面数据提取与可视化 |
| FarFieldGen | ~LOC | ✅ 完成 | 远场方向图表面生成 |
| PhaseAnimation | ~LOC | ✅ 完成 | 相位扫描动画 |

> **注意**：目前渲染引擎仅在**合成测试数据**上运行，尚未接入真实仿真结果。

### 6.3 Touchstone 引擎详细状态

| 子模块 | 状态 | 功能 |
|--------|------|------|
| Parser | ✅ 完成 | v1.0 & v2.0 格式解析，含行号错误报告 |
| Writer | ✅ 完成 | 标准格式写入 |
| Types | ✅ 完成 | 全参数类型（S/Y/Z/H/G）、全数据格式（RI/MA/dB） |
| 复数运算 | ✅ 完成 | 阻抗归一化、格式转换 |

### 6.4 Domain 模型详细状态

领域模型已对齐设计文档，20 个子模块覆盖完整工程文件结构：

| 子模块 | 状态 | 功能 |
|--------|------|------|
| solution_type | ✅ 完成 | SolutionType 枚举（HFSS 5 种 + Q3D 4 种求解类型） |
| variable | ✅ 完成 | Variable、PropertyValue（Constant/Expression/Dataset）、DatasetDefinition |
| material | ✅ 完成 | MaterialDef、MaterialCategory、MaterialProperties（6 种 EM 属性） |
| geometry | ✅ 完成 | Geometry 历史记录 + 快照、GeometryOperation（30 种命令）、GeoObject |
| boundary | ✅ 完成 | Boundary（11 种边界类型）、Assignment |
| excitation | ✅ 完成 | Excitation（7 种激励类型，含 HFSS 端口 + Q3D 源/汇） |
| net | ✅ 完成 | Net、Terminal（Q3D 网络定义） |
| mesh | ✅ 完成 | MeshOperation（4 种网格控制） |
| analysis | ✅ 完成 | AnalysisSetup、FrequencySweep（3 种扫描）、DcSettings |
| radiation | ✅ 完成 | RadiationSetup、FarFieldSetup、NearFieldSetup（Line/Rectangle/Sphere） |
| output_variable | ✅ 完成 | OutputVariable（结果派生表达式） |
| field_overlay | ✅ 完成 | FieldOverlay（8 种场量、4 种绘图类型） |
| optimetrics | ✅ 完成 | OptimetricsSetup（5 种 tagged enum：参数扫描/优化/灵敏度/统计/调参） |
| report | ✅ 完成 | Report（8 种报告类别、6 种图表类型）、Trace、Axis、Marker |
| solution_index | ✅ 完成 | SolutionIndex、SetupSolutionStatus、SolveStatus、RlcgSummary |
| design | ✅ 完成 | Design 聚合、DesignSettings、Definitions、CoordinateSystem、NamedSelection |
| project | ✅ 完成 | EmProject 顶层、ProjectMetadata |
| expression | ✅ 完成 | 递归下降解析器（四则运算/函数/变量/单位）、变量求值、循环依赖检测 |
| validation | ✅ 完成 | 材料/对象引用检查、命名唯一性、HFSS/Q3D 特定验证 |
| dependency | ✅ 完成 | DependencyGraph、定义→引用分析、stale 标记 |
| file_io | ✅ 完成 | JSON load/save、auto-save/recovery、ProjectLock、结果目录管理 |

### 6.5 Web Worker / OPFS 状态

| 子模块 | 状态 | 功能 |
|--------|------|------|
| worker_protocol | ✅ 完成 | WorkerCommand/WorkerResponse 枚举（MessagePack 编码） |
| opfs | ✅ 完成 | OPFS API 封装（目录创建/文件读写删除/目录遍历） |
| worker entry | ✅ 完成 | onmessage 命令分发 + 项目索引管理 |
| wasm_backend | ✅ 完成 | WasmBackend 实现 Backend trait（channel 桥接 Worker） |

---

## 7. 开发里程碑

### Milestone 0：基础框架 ✅ 已完成

> **目标**：搭建可运行的桌面应用框架

| 任务 | 模块 | 状态 |
|------|------|------|
| Cargo Workspace 搭建 | 根 | ✅ 完成 |
| egui + eframe 桌面应用入口 | `main` | ✅ 完成 |
| WASM 入口代码 | `main` | ✅ 完成（未测试） |
| Ribbon Bar 工具栏 | `components` | ✅ 完成 |
| Dock 面板 + Tab 页布局 | `components` / `app` | ✅ 完成 |
| 基础领域模型 (Project/Design/Material) | `domain` | ✅ 完成 |
| Backend trait 抽象 | `infra` | ✅ 完成 |
| Standalone 后端 | `infra` | ✅ 完成 |
| Solver trait 抽象 | `solver` | ✅ 完成 |
| 内存中工程创建/重置 | `app` | ✅ 完成 |
| 操作日志 | `app` | ✅ 完成 |

### Milestone 1：3D 渲染引擎 ✅ 已完成

> **目标**：建立完整的 wgpu 3D 可视化管线

| 任务 | 模块 | 状态 |
|------|------|------|
| wgpu 设备初始化 + 离屏渲染 | `render` | ✅ 完成 |
| WGSL 着色器编写 | `render` | ✅ 完成 |
| 轨道相机（旋转/缩放/平移） | `render` | ✅ 完成 |
| 色标系统（4 种 colormap） | `render` | ✅ 完成 |
| 表面云图渲染 | `render` | ✅ 完成 |
| 矢量箭头 Instanced 渲染 | `render` | ✅ 完成 |
| 切面可视化 | `render` | ✅ 完成 |
| 远场方向图表面 | `render` | ✅ 完成 |
| 相位动画 | `render` | ✅ 完成 |
| 渲染结果 blit 到 egui | `render` | ✅ 完成 |
| field-vis 示例程序 | `examples` | ✅ 完成 |

### Milestone 2：S 参数文件支持 ✅ 已完成

> **目标**：完整的 Touchstone 文件读写能力

| 任务 | 模块 | 状态 |
|------|------|------|
| Touchstone v1.0 解析 | `touchstone` | ✅ 完成 |
| Touchstone v2.0 解析 | `touchstone` | ✅ 完成 |
| 全参数类型支持 (S/Y/Z/H/G) | `touchstone` | ✅ 完成 |
| 全数据格式支持 (RI/MA/dB) | `touchstone` | ✅ 完成 |
| 文件写入 | `touchstone` | ✅ 完成 |
| 错误恢复 + 行号定位 | `touchstone` | ✅ 完成 |

### Milestone 3：工程文件 I/O ✅ 已完成

> **目标**：实现 .emsp 工程文件的完整读写

| 任务 | 模块 | 状态 |
|------|------|------|
| 完整领域模型实现（对齐设计文档） | `domain` | ✅ 完成 |
| EmProject 序列化/反序列化 | `domain` | ✅ 完成 |
| 工程完整性验证 | `domain` | ✅ 完成 |
| 定义-引用依赖图 | `domain` | ✅ 完成 |
| 变量表达式求值引擎 | `domain` | ✅ 完成 |
| .emsp.lock 并发锁机制 | `domain` | ✅ 完成 |
| .emsp.auto 自动保存/恢复 | `domain` | ✅ 完成 |
| 结果目录创建与管理 | `domain` | ✅ 完成 |
| 结果索引 (SolutionIndex) 跟踪 | `domain` | ✅ 完成 |
| 结果过期 (is_stale) 机制 | `domain` | ✅ 完成 |
| 版本迁移策略实现 | `domain` | ✅ 完成 |

### Milestone 4：几何建模 🔄 部分完成

> **目标**：基于 rcad B-Rep 内核支持参数化 3D 几何建模
>
> **依赖**：rcad（git submodule，已引入 vendor/rcad），提供 B-Rep 拓扑、体素创建、布尔运算、扫掠、STEP 导入导出

| 任务 | 模块 | rcad 对应 | 状态 |
|------|------|-----------|------|
| 操作历史记录引擎 | `domain` | — | ✅ 完成 |
| 基本体素创建 (Box/Cylinder/Sphere/Cone/Torus) | `domain` | `rcad-modeling` (box/cylinder/sphere/cone/torus_brep) | ✅ 完成 |
| 布尔运算 (Unite/Subtract/Intersect) | `domain` | `rcad-algorithms` (boolean_op) | ✅ 完成 |
| 几何变换 (Move/Rotate/Scale/Mirror) | `domain` | `rcad-kernel` (BRep::apply_transform, DAffine3) | ✅ 完成 |
| 扫掠操作 (Extrude/Revolve/SweepPipe) | `domain` | `rcad-modeling` (extrude/revolve/sweep_pipe) | 🔲 未开始 |
| 参数化尺寸（变量引用） | `domain` | — | 🔲 未开始 |
| 3D 几何渲染（实体 + 线框） | `render` | `rcad-render` (Tessellator/WgpuRenderer) | 🔲 未开始 |
| 几何拾取与选择 | `render` | `rcad-render` (SelectionState, pick_face/pick_edge) | 🔲 未开始 |
| 坐标系可视化 | `render` | — | 🔲 未开始 |
| CAD 导入 (STEP/STL) | `domain` | `rcad-step` (STEP import/export with colors) | 🔲 未开始 |
| 几何建模 UI 面板 | `components` | — | 🔲 未开始 |

### Milestone 5：求解器集成 🔲 未开始

> **目标**：集成 Rem 求解器，实现完整仿真流程

| 任务 | 模块 | 状态 |
|------|------|------|
| Rem 库 Rust 绑定 | `solver` | 🔲 未开始 |
| HFSS FEM 求解调度 | `solver` | 🔲 未开始 |
| Q3D MoM 求解调度 | `solver` | 🔲 未开始 |
| 自适应网格加密循环 | `solver` | 🔲 未开始 |
| 频率扫描执行 | `solver` | 🔲 未开始 |
| 仿真进度回调 | `solver` → `app` | 🔲 未开始 |
| 收敛历史实时写入 | `solver` | 🔲 未开始 |
| 场数据导出 (.emsfld) | `solver` | 🔲 未开始 |
| S 参数提取与写入 | `solver` | 🔲 未开始 |
| RLCG 矩阵提取与写入 | `solver` | 🔲 未开始 |
| 远场/近场积分计算 | `solver` | 🔲 未开始 |
| 模型验证 (Validation) | `solver` | 🔲 未开始 |
| 求解器日志 | `solver` | 🔲 未开始 |

### Milestone 6：2D 报告系统 🔲 未开始

> **目标**：实现 HFSS/Q3D 结果的 2D 图表可视化

| 任务 | 模块 | 状态 |
|------|------|------|
| ResultDataStore 数据管理 | `domain` | 🔲 未开始 |
| JSON/Touchstone 结果加载 | `domain` | 🔲 未开始 |
| Trace/QuantityExpression 系统 | `domain` | 🔲 未开始 |
| S 参数矩形图 (egui_plot) | `components` | 🔲 未开始 |
| Smith 圆图 | `components` | 🔲 未开始 |
| 极坐标方向图 | `components` | 🔲 未开始 |
| 收敛曲线（HFSS: Delta S / Q3D: Delta Energy） | `components` | 🔲 未开始 |
| RLCG 矩阵 vs 频率曲线 (Q3D) | `components` | 🔲 未开始 |
| RLCG 矩阵数据表格 + 热力图 (Q3D) | `components` | 🔲 未开始 |
| Marker / Cursor 交互 | `components` | 🔲 未开始 |
| 数据表格 | `components` | 🔲 未开始 |
| 3D 矩形图（参数扫描） | `components` | 🔲 未开始 |
| ReportPanel 集成到 Dock | `app` | 🔲 未开始 |
| CSV 导出 | `components` | 🔲 未开始 |

### Milestone 7：3D 场数据管线 🔲 未开始

> **目标**：连通真实仿真结果到 3D 渲染管线

| 任务 | 模块 | 状态 |
|------|------|------|
| Gmsh MSH 4.1 加载器 | `domain` | 🔲 未开始 |
| .emsfld 二进制内存映射加载 | `domain` | 🔲 未开始 |
| 网格表面提取（四面体 → 三角形） | `render` | 🔲 未开始 |
| 场数据 → 顶点颜色映射 | `render` | 🔲 未开始 |
| 网格线框/实体渲染（真实数据） | `render` | 🔲 未开始 |
| 表面云图渲染（真实场数据） | `render` | 🔲 未开始 |
| 矢量箭头渲染（真实场数据） | `render` | 🔲 未开始 |
| 切面可视化（真实场数据） | `render` | 🔲 未开始 |
| 等值面提取（Marching Tetrahedra） | `render` | 🔲 未开始 |
| 3D 远场方向图（真实数据） | `render` | 🔲 未开始 |
| 网格质量可视化 | `render` | 🔲 未开始 |
| GPU Picking + 场值探针 | `render` | 🔲 未开始 |
| PNG 截图导出 | `render` | 🔲 未开始 |
| VTK 导出（ParaView 兼容） | `domain` | 🔲 未开始 |

### Milestone 8：Q3D 专项功能 🔲 未开始

> **目标**：完整的 Q3D 准静态分析工作流

| 任务 | 模块 | 状态 |
|------|------|------|
| RlcgMatrixData 加载/序列化 | `domain` | 🔲 未开始 |
| 等效电路 SPICE 导出 | `domain` | 🔲 未开始 |
| RLCG → S 参数转换 | `domain` | 🔲 未开始 |
| Q3D 电流密度场叠加 | `render` | 🔲 未开始 |
| Q3D 电荷分布场叠加 | `render` | 🔲 未开始 |
| Q3D 欧姆损耗场叠加 | `render` | 🔲 未开始 |
| Q3D Net/Terminal 编辑 UI | `components` | 🔲 未开始 |
| Q3D 验证检查（Net/Terminal） | `solver` | 🔲 未开始 |

### Milestone 9：Optimetrics 🔲 未开始

> **目标**：参数扫描与优化

| 任务 | 模块 | 状态 |
|------|------|------|
| 参数扫描引擎 | `solver` | 🔲 未开始 |
| 优化算法（QN/Pattern/GA/SNLP） | `solver` | 🔲 未开始 |
| 灵敏度分析 | `solver` | 🔲 未开始 |
| 统计分析（蒙特卡罗） | `solver` | 🔲 未开始 |
| 交互式调参 (Tuning) | `app` | 🔲 未开始 |
| 优化结果汇总 (summary.json) | `solver` | 🔲 未开始 |
| 优化收敛曲线 | `components` | 🔲 未开始 |
| 参数化对比图 | `components` | 🔲 未开始 |

### Milestone 10：平台与部署 🔄 部分完成

> **目标**：WASM 部署、Cloud 模式、版本分级与 Local-First 支持

| 任务 | 模块 | 状态 |
|------|------|------|
| WASM 构建验证 (trunk) | `main` | 🔲 未开始 |
| WebGPU 适配测试 | `render` | 🔲 未开始 |
| Cloud 后端 REST API | `infra` | 🔲 未开始 |
| 远程求解任务提交 | `infra` | 🔲 未开始 |
| 结果文件下载/流式加载 | `infra` | 🔲 未开始 |
| 多用户工程锁协调 | `infra` | 🔲 未开始 |
| Edition 功能门控 (Basic/Pro/Enterprise) | `app` / `domain` | 🔲 未开始 |
| Web 版工程管理适配（禁用新建/另存） | `app` | 🔲 未开始 |
| OPFS 存储层（Web Worker 内） | `worker` | ✅ 完成 |
| Worker 通信协议 + WasmBackend | `domain` / `infra` | ✅ 完成 |
| Rem WASM 编译 + Worker 集成 | `solver` | 🔲 未开始 |
| Local-First 离线缓存（Service Worker） | `main` | 🔲 未开始 |
| Web 版预创建工程加载流程 | `infra` / `app` | 🔲 未开始 |

---

## 8. 开发进度总览

```
Milestone 0: 基础框架         [████████████████████] 100%  ✅
Milestone 1: 3D 渲染引擎      [████████████████████] 100%  ✅
Milestone 2: Touchstone 支持   [████████████████████] 100%  ✅
Milestone 3: 工程文件 I/O      [████████████████████] 100%  ✅
Milestone 4: 几何建模          [████████            ]  36%  🔄
Milestone 5: 求解器集成        [                    ]   0%  🔲
Milestone 6: 2D 报告系统       [                    ]   0%  🔲
Milestone 7: 3D 场数据管线     [                    ]   0%  🔲
Milestone 8: Q3D 专项功能      [                    ]   0%  🔲
Milestone 9: Optimetrics       [                    ]   0%  🔲
Milestone 10: 平台与部署        [████                ]  20%  🔄
─────────────────────────────────────────────────────────────
整体进度                       [████████            ]  42%
```

### 建议开发优先级

```
                         ┌───────────────────┐
                         │  M3: 工程文件 I/O  │  ← 最高优先：所有功能的基础
                         └────────┬──────────┘
                    ┌─────────────┼─────────────┐
                    ▼             ▼              ▼
            ┌──────────┐  ┌─────────────┐ ┌──────────┐
            │ M4: 建模  │  │ M5: 求解器  │ │ M6: 报告 │
            └─────┬────┘  └──────┬──────┘ └─────┬────┘
                  │              │               │
                  └──────────────┼───────────────┘
                                 ▼
                    ┌────────────────────────┐
                    │ M7: 3D 场数据管线       │  ← 连通渲染引擎与真实数据
                    └────────────┬───────────┘
                    ┌────────────┼───────────┐
                    ▼            ▼            ▼
              ┌──────────┐ ┌──────────┐ ┌──────────┐
              │ M8: Q3D  │ │M9: Optim │ │M10: 部署 │
              └──────────┘ └──────────┘ └──────────┘
```

**关键路径**：M3（工程 I/O）→ M5（求解器）→ M7（场数据管线）

这三个里程碑打通后，EMStudio 将具备完整的「建模 → 仿真 → 可视化」工作流。

---

## 9. 设计文档索引

| 文档 | 内容 | 规模 |
|------|------|------|
| [em-project-file-design.md](em-project-file-design.md) | 工程文件格式规范：层次结构、定义-引用架构、完整 JSON Schema、Rust 类型映射、文件操作 API | ~3,900 行 |
| [em-result-file-formats.md](em-result-file-formats.md) | 仿真结果文件格式：JSON 结果文件（验证/收敛/S 参数/RLCG 矩阵/等效电路）、二进制格式（MSH/emsfld）、Rust 类型定义 | ~1,760 行 |
| [em-result-visualization-design.md](em-result-visualization-design.md) | 可视化系统：2D 报告（S 参数/Smith/极坐标/RLCG）、3D 场图（云图/矢量/切面/动画）、交互系统、导出功能、实现路线 | ~1,810 行 |
| [em-feature-design-and-progress.md](em-feature-design-and-progress.md) | 本文档：功能全景、代码现状、开发里程碑、进度追踪 | — |
