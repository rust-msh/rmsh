# EMStudio 仿真结果可视化方案设计

## 1. 概述

本文档定义 EMStudio 仿真结果的可视化系统设计。参考 Ansys HFSS/AEDT 的后处理体系（Reports、Field Overlays、Far Field Plots），结合 EMStudio 自身技术栈（Rust + egui + egui_plot + wgpu）进行适配设计。

### 1.1 设计目标

- **所见即所得**：结果与 3D 几何模型叠加显示，支持交互式旋转/缩放/切片
- **实时响应**：利用 wgpu GPU 加速渲染大规模场数据和网格，保证交互帧率 ≥ 30 fps
- **按需加载**：基于 `.emsfld` 二进制文件的随机访问能力，仅加载当前频点/分量的数据
- **分层架构**：数据加载、可视化映射、GPU 渲染三层分离，易于扩展新的可视化类型
- **跨平台**：Native（wgpu 硬件加速）和 WASM（WebGPU / 软件回退）统一 API

### 1.2 可视化类型总览

| 类别 | 可视化类型 | 数据来源 | 渲染方式 |
|------|-----------|----------|---------|
| 2D 报告图表 | S 参数曲线、Smith 圆图 | `s_parameters.json` / `.snp` | egui_plot 2D |
| 2D 报告图表 | 收敛曲线、网格增长曲线 | `convergence.json` | egui_plot 2D |
| 2D 报告图表 | 远场方向图（极坐标/直角坐标） | `far_field_*.json` | egui_plot 2D |
| 2D 报告图表 | 优化收敛曲线 | `summary.json` | egui_plot 2D |
| 3D 场图 | E/H/J 场幅度云图 | `.emsfld` + `.msh` | wgpu 3D |
| 3D 场图 | 矢量场箭头图 | `.emsfld` + `.msh` | wgpu 3D |
| 3D 场图 | 表面电流分布 | `.emsfld` + `.msh` | wgpu 3D |
| 3D 场图 | 动画（相位扫描） | `.emsfld` + `.msh` | wgpu 3D |
| 3D 远场 | 3D 辐射方向图 | `far_field_*.json` | wgpu 3D |
| 3D 网格 | 网格质量可视化 | `.msh` | wgpu 3D |
| 数据表格 | S 参数数据表 | `s_parameters.json` | egui Table |
| 数据表格 | 天线参数摘要 | `far_field_*.json` | egui Table |
| 2D 报告图表 | RLCG 参数 vs 频率（Q3D） | `rlcg_matrix.json` | egui_plot 2D |
| 数据表格 | RLCG 矩阵表格（Q3D） | `rlcg_matrix.json` | egui Table |
| 3D 场图 | Q3D 电流密度分布 | `.emsfld` + `.msh` | wgpu 3D |
| 3D 场图 | Q3D 电场/电荷分布 | `.emsfld` + `.msh` | wgpu 3D |

### 1.3 参考：HFSS 后处理体系

HFSS 的结果可视化分为三大子系统：

| HFSS 子系统 | 功能 | EMStudio 对应 |
|------------|------|-------------|
| **Reports** | 2D 曲线图（矩形、极坐标、Smith 圆图、数据表） | ReportPanel（egui_plot） |
| **Field Overlays** | 3D 场叠加显示（云图、矢量、等值面、切面） | FieldOverlayRenderer（wgpu） |
| **Far Field Plots** | 3D 辐射方向图、2D 极坐标方向图 | FarFieldRenderer（wgpu + egui_plot） |
| **Mesh Display** | 网格可视化、质量着色 | MeshRenderer（wgpu） |

### 1.4 参考：Q3D Extractor 后处理体系

Q3D Extractor 的结果可视化与 HFSS 有相似的框架但侧重不同：参考 [Q3D Reports](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/ReportsandPostProc/CreatingReportsQ3D.htm) 和 [Q3D Field Overlays](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Q3DExtractor/Content/ReportsandPostProc/PlottingFieldOverlaysinQ3D.htm)。

| Q3D 子系统 | 功能 | EMStudio 对应 |
|-----------|------|-------------|
| **Matrix Reports** | RLCG 矩阵元素 vs 频率曲线、矩阵数据表格 | ReportPanel（egui_plot + egui Table） |
| **Field Overlays** | 电流密度 (J)、电场 (E)、电荷分布 (ρ) 的 3D 叠加 | FieldOverlayRenderer（wgpu） |
| **Equivalent Circuit** | 等效电路模型导出、SPICE 网表生成 | 导出功能（非可视化） |
| **S-Parameter Conversion** | RLCG → S 参数转换与 Touchstone 导出 | ReportPanel（egui_plot） |

> **Q3D 与 HFSS 后处理差异**：Q3D 没有远场辐射方向图（准静态求解无远场概念），也没有 Smith 圆图（不直接输出 S 参数，但可通过 RLCG → S 参数转换后使用）。Q3D 的核心结果是 **RLCG 矩阵**，其可视化以矩阵元素随频率变化的曲线和矩阵数据表格为主。Q3D 的场叠加支持电流密度和电场分布，用于分析电流拥挤、趋肤效应和耦合路径。

---

## 2. 系统架构

### 2.1 分层架构

```
┌─────────────────────────────────────────────────────────────────┐
│                        UI Layer (egui)                          │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌───────────────────┐  │
│  │ReportPanel│ │FieldPanel│ │FarField  │ │PropertyEditor     │  │
│  │(egui_plot)│ │(3D view) │ │Panel     │ │(plot options,     │  │
│  │          │ │          │ │          │ │ colormap, range)  │  │
│  └────┬─────┘ └────┬─────┘ └────┬─────┘ └────────┬──────────┘  │
│       │            │            │                 │              │
├───────┼────────────┼────────────┼─────────────────┼──────────────┤
│       │      Visualization Mapping Layer          │              │
│  ┌────▼─────────────────────────▼─────────────────▼──────────┐  │
│  │                   VisDataPipeline                          │  │
│  │  ┌─────────────┐ ┌──────────────┐ ┌────────────────────┐  │  │
│  │  │QuantityExpr │ │ ColorMapper  │ │ GeometryMapper     │  │  │
│  │  │(dB, phase,  │ │ (colormap,   │ │ (surface extract,  │  │  │
│  │  │ mag, real)  │ │  range, log) │ │  iso-surface,      │  │  │
│  │  └─────────────┘ └──────────────┘ │  slice plane)      │  │  │
│  │                                    └────────────────────┘  │  │
│  └────────────────────────┬──────────────────────────────────┘  │
│                           │                                      │
├───────────────────────────┼──────────────────────────────────────┤
│                    Data Access Layer                              │
│  ┌────────────────────────▼──────────────────────────────────┐  │
│  │                   ResultDataStore                          │  │
│  │  ┌───────────┐ ┌──────────┐ ┌──────────┐ ┌────────────┐  │  │
│  │  │JsonLoader  │ │MshLoader │ │FldLoader │ │SnpLoader   │  │  │
│  │  │(serde)    │ │(gmsh)    │ │(mmap)    │ │(touchstone)│  │  │
│  │  └───────────┘ └──────────┘ └──────────┘ └────────────┘  │  │
│  └───────────────────────────────────────────────────────────┘  │
│                                                                  │
├──────────────────────────────────────────────────────────────────┤
│                    GPU Render Layer (wgpu)                        │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │                   SceneRenderer                           │   │
│  │  ┌────────────┐ ┌─────────────┐ ┌──────────────────┐    │   │
│  │  │MeshPipeline│ │FieldPipeline│ │FarFieldPipeline  │    │   │
│  │  │(wireframe, │ │(colormap    │ │(3D pattern       │    │   │
│  │  │ solid,     │ │ texture,    │ │ surface mesh)    │    │   │
│  │  │ quality)   │ │ arrows)     │ │                  │    │   │
│  │  └────────────┘ └─────────────┘ └──────────────────┘    │   │
│  │  ┌────────────┐ ┌─────────────┐ ┌──────────────────┐    │   │
│  │  │Camera      │ │LightSystem  │ │PickingSystem     │    │   │
│  │  │(orbit,pan, │ │(ambient +   │ │(GPU color-id     │    │   │
│  │  │ zoom)      │ │ directional)│ │ picking)         │    │   │
│  │  └────────────┘ └─────────────┘ └──────────────────┘    │   │
│  └──────────────────────────────────────────────────────────┘   │
└──────────────────────────────────────────────────────────────────┘
```

### 2.2 模块职责

| 模块 | crate | 职责 |
|------|-------|------|
| `ResultDataStore` | `emstudio-domain` | 结果文件索引、加载、缓存；提供统一的数据访问 trait |
| `VisDataPipeline` | `emstudio-domain` | 将原始数值数据映射为可视化数据（颜色、位置、箭头） |
| `ReportPanel` | `emstudio-components` | 基于 egui_plot 的 2D 报告图表（曲线、极坐标、Smith 圆图） |
| `FieldOverlayPanel` | `emstudio-components` | 3D 场叠加的 UI 控制面板（量选择、色标、范围、切面） |
| `SceneRenderer` | `emstudio-render` | wgpu GPU 渲染管线管理，网格/场/远场的 3D 渲染 |
| `EmStudioApp` | `emstudio-app` | 将上述模块集成到 DockArea Tab 系统中 |

---

## 3. 2D 报告系统（Reports）

### 3.1 设计思想

参考 HFSS Reports 模块：用户创建一个 Report，选择报告类型（矩形图/极坐标/Smith 圆图/数据表），然后向报告中添加 Trace（数据曲线），每条 Trace 对应一个量表达式（如 `dB(S(1,1))`）。

EMStudio 报告系统核心概念：

```
Report
 ├── report_type: ReportType       # 图表类型
 ├── traces: Vec<Trace>            # 数据曲线列表
 ├── x_axis: AxisConfig            # X 轴配置
 ├── y_axis: Vec<AxisConfig>       # Y 轴配置（支持双 Y 轴）
 └── display: DisplayConfig        # 显示选项
```

### 3.2 报告类型

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ReportType {
    /// 矩形坐标图（最常用：S参数 vs 频率）
    Rectangular,
    /// 极坐标图（远场方向图）
    Polar,
    /// Smith 圆图（阻抗匹配分析）
    SmithChart,
    /// 数据表格
    DataTable,
    /// 3D 矩形图（参数扫描结果：X=频率，Y=参数值，Z=S参数）
    Rectangular3D,
}
```

### 3.3 Trace 与量表达式

参考 HFSS 的量表达式体系，定义 Trace 数据源：

```rust
/// 一条数据曲线
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Trace {
    pub name: String,
    /// 量表达式，如 "dB(S(1,1))", "ang_deg(S(2,1))", "re(Z(1,1))"
    pub expression: QuantityExpression,
    /// 数据来源
    pub data_source: TraceDataSource,
    /// 显示样式
    pub style: TraceStyle,
}

/// 量表达式
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuantityExpression {
    /// 基础量：S/Y/Z/H/G 参数、远场增益、近场分量等
    pub base_quantity: BaseQuantity,
    /// 变换函数链：dB(), mag(), ang_deg(), re(), im(), vswr() 等
    pub transforms: Vec<Transform>,
}

/// 基础量枚举
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BaseQuantity {
    /// 网络参数：S(i,j), Y(i,j), Z(i,j)
    NetworkParameter {
        param_type: NetworkParamType,  // S, Y, Z, H, G
        port_i: u32,
        port_j: u32,
    },
    /// 远场量：GainTotal, GainTheta, GainPhi, DirectivityTotal, AxialRatio
    FarFieldQuantity {
        quantity: FarFieldQuantityType,
        cut_plane: CutPlane,  // phi=0, phi=90, theta=90, etc.
    },
    /// 近场量：E_mag, H_mag, Poynting_mag
    NearFieldQuantity {
        quantity: NearFieldQuantityType,
    },
    /// 天线参数：PeakGain, Efficiency, Beamwidth
    AntennaParameter {
        parameter: AntennaParamType,
    },
    /// Q3D RLCG 矩阵元素：R(i,j), L(i,j), C(i,j), G(i,j)
    RlcgParameter {
        param_type: RlcgParamType,  // R, L, C, G
        terminal_i: u32,
        terminal_j: u32,
    },
    /// Q3D 收敛数据：DeltaEnergy, NumTriangles, RLCG 极值
    Q3dConvergenceQuantity {
        quantity: Q3dConvergenceQuantityType,
    },
    /// 收敛数据：MaxDeltaS, NumTetrahedra
    ConvergenceQuantity {
        quantity: ConvergenceQuantityType,
    },
    /// 输出变量（用户自定义表达式）
    OutputVariable {
        name: String,
    },
}

/// 变换函数
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Transform {
    /// 取模 |z|
    Magnitude,
    /// 分贝 20*log10(|z|)
    DB,
    /// 相位（度）
    AngleDeg,
    /// 相位（弧度）
    AngleRad,
    /// 实部
    Real,
    /// 虚部
    Imaginary,
    /// VSWR = (1+|Γ|)/(1-|Γ|)
    VSWR,
    /// 群延迟 -d(phase)/d(freq)
    GroupDelay,
    /// 归一化
    Normalize,
    /// 对数
    Log10,
}
```

### 3.4 Trace 数据源

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraceDataSource {
    /// 分析设置名称
    pub setup: String,
    /// 频率扫描名称（None 则取自适应频率结果）
    pub sweep: Option<String>,
    /// 参数化变量值（用于参数扫描对比）
    pub variation: Option<HashMap<String, String>>,
}
```

### 3.5 S 参数矩形图

最常用的报告类型。X 轴为频率，Y 轴为 S 参数幅度/相位。

**交互功能**（参考 HFSS）：

| 功能 | 说明 | 实现方式 |
|------|------|---------|
| Marker | 在曲线上放置标记点，显示精确数值 | egui_plot 的 Points + Text |
| Delta Marker | 两个 Marker 之间的差值 | 计算两个 Marker 坐标差 |
| Cursor | 跟随鼠标的十字光标，显示最近曲线值 | egui_plot pointer coordinate |
| Zoom | 框选放大、滚轮缩放 | egui_plot 内置 |
| 多 Y 轴 | 幅度（左 Y 轴）+ 相位（右 Y 轴） | egui_plot 双坐标系叠加 |
| Trace 叠加 | 多条 Trace 同图对比 | 多条 Line 同一 Plot |
| 参数化叠加 | 不同参数值的 S 曲线对比 | 按 variation 分色 |
| 导出 | 导出为 CSV / PNG / Clipboard | CSV 写文件，PNG 截图 |

**数据流**：

```
s_parameters.json / .snp
       │
       ▼
 SnpLoader / JsonLoader
       │
       ▼
 QuantityExpression::evaluate()
   (raw complex → dB/phase/mag/...)
       │
       ▼
 Vec<PlotPoint> { freq, value }
       │
       ▼
 egui_plot::Line
       │
       ▼
 ReportPanel::ui()
```

**Rust 接口**：

```rust
pub struct ReportPanel {
    pub report: Report,
    /// 已加载的 Trace 数据缓存
    trace_cache: HashMap<String, Vec<[f64; 2]>>,
    /// Marker 列表
    markers: Vec<Marker>,
    /// 是否需要重新加载
    dirty: bool,
}

impl ReportPanel {
    /// 创建新的报告面板
    pub fn new(report: Report) -> Self { todo!() }

    /// 从 ResultDataStore 加载 Trace 数据
    pub fn load_traces(&mut self, store: &ResultDataStore) { todo!() }

    /// 渲染到 egui UI
    pub fn ui(&mut self, ui: &mut egui::Ui) { todo!() }

    /// 导出为 CSV
    pub fn export_csv(&self, path: &Path) -> Result<(), ExportError> { todo!() }
}
```

### 3.6 Smith 圆图

显示 S 参数在 Smith 圆图上的轨迹，用于阻抗匹配分析。

**渲染要素**：

| 要素 | 说明 | 绘制方式 |
|------|------|---------|
| 等电阻圆 | r = 0, 0.5, 1, 2, 5, ∞ | egui_plot 参数化曲线 |
| 等电抗弧 | x = ±0.5, ±1, ±2, ±5 | egui_plot 参数化曲线 |
| S 参数轨迹 | 频率扫描的 S₁₁ 轨迹 | egui_plot Line（实部 vs 虚部） |
| Marker | 标记特定频率点 | egui_plot Points + 标注文字 |
| 归一化阻抗 | 鼠标悬停显示 Z/Y 值 | tooltip |

**数据变换**：

```rust
/// S参数 → Smith 圆图坐标
fn s_to_smith(s_real: f64, s_imag: f64) -> (f64, f64) {
    // Smith 圆图直接使用 S 参数的实部和虚部作为坐标
    // |S| ≤ 1 时在单位圆内
    (s_real, s_imag)
}

/// S参数 → 归一化阻抗
fn s_to_z_normalized(s_real: f64, s_imag: f64, z0: f64) -> (f64, f64) {
    // z = z0 * (1 + S) / (1 - S)
    let denom_r = (1.0 - s_real) * (1.0 - s_real) + s_imag * s_imag;
    let z_real = z0 * ((1.0 - s_real * s_real - s_imag * s_imag) / denom_r);
    let z_imag = z0 * (2.0 * s_imag / denom_r);
    (z_real, z_imag)
}
```

### 3.7 极坐标方向图

显示远场增益的极坐标切面图。

**典型切面**：

| 切面 | 含义 | 典型用途 |
|------|------|---------|
| φ = 0° (E-Plane) | XZ 平面方向图 | E 面波束宽度 |
| φ = 90° (H-Plane) | YZ 平面方向图 | H 面波束宽度 |
| θ = 90° | XY 平面（水平面）方向图 | 全向性评估 |

**渲染要素**：

| 要素 | 说明 |
|------|------|
| 同心圆 | dB 刻度圆（0, -3, -10, -20, -30 dB） |
| 角度射线 | 每 30° 或 15° 一条 |
| 增益曲线 | 极坐标下的 Gain(θ) 或 Gain(φ) |
| 3dB 波束宽度 | 高亮标注 |
| HPBW 标记 | 半功率波束宽度数值 |

**数据流**：

```
far_field_*.json
      │
      ▼
 FarFieldData.derived_quantities.GainTotal
      │
      ▼
 extract_cut(phi=0°) → Vec<(angle_deg, gain_dbi)>
      │
      ▼
 polar_to_cartesian → Vec<[f64; 2]>
      │
      ▼
 egui_plot::Polygon / Line（极坐标变换为直角坐标绘制）
```

### 3.8 收敛曲线

显示自适应求解过程中的收敛趋势。

**图表内容**：

| 子图 | X 轴 | Y 轴 | 说明 |
|------|------|------|------|
| Delta S 收敛 | Pass number | Max Delta S | 带目标线（target_max_delta_s） |
| 网格增长 | Pass number | Num tetrahedra | 柱状图或折线图 |
| 求解时间 | Pass number | Solve time (sec) | 分阶段堆叠柱状图 |
| 内存占用 | Pass number | Peak memory (MB) | 折线图 |

### 3.9 Q3D RLCG 矩阵曲线图

显示 Q3D 准静态提取结果中 RLCG 矩阵元素随频率变化的趋势。参考 [Q3D Matrix Reports](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/ReportsandPostProc/CreatingReportsQ3D.htm)。

**图表内容**：

| 子图 | X 轴 | Y 轴 | 说明 |
|------|------|------|------|
| R vs Freq | Frequency | Resistance (Ω) | 自阻/互阻随频率变化（趋肤效应） |
| L vs Freq | Frequency | Inductance (nH) | 自感/互感随频率变化 |
| C Matrix | Frequency | Capacitance (pF) | 通常与频率无关，显示为水平线 |
| G vs Freq | Frequency | Conductance (mS) | 介质损耗随频率变化 |

**Trace 表达式示例**（Q3D 量表达式）：

| 表达式 | 含义 |
|--------|------|
| `R(Signal1:T1, Signal1:T1)` | Signal1 自电阻 |
| `R(Signal1:T1, Signal2:T3)` | Signal1 与 Signal2 互电阻 |
| `L(Signal1:T1, Signal1:T1)` | Signal1 自感 |
| `L(Signal1:T1, Signal2:T3)` | Signal1 与 Signal2 互感 |
| `C(Signal1:T1, Signal1:T1)` | Signal1 自电容 |
| `C(Signal1:T1, Signal2:T3)` | Signal1 与 Signal2 互电容（负值） |
| `G(Signal1:T1, Signal1:T1)` | Signal1 自电导 |

**数据流**：

```
rlcg_matrix.json
       │
       ▼
 RlcgMatrixData.matrices.R.data_per_frequency
       │
       ▼
 extract R[i][j] at each frequency
       │
       ▼
 Vec<PlotPoint> { freq, value }
       │
       ▼
 egui_plot::Line（多条 Trace：自阻 + 互阻同图对比）
       │
       ▼
 ReportPanel::ui()
```

**Rust 接口**：

```rust
/// Q3D RLCG 矩阵报告面板
pub struct RlcgReportPanel {
    pub report: Report,
    pub traces: Vec<RlcgTrace>,
    pub matrix_type_filter: RlcgParamType,  // R, L, C, G
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RlcgParamType {
    R,  // 电阻
    L,  // 电感
    C,  // 电容
    G,  // 电导
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Q3dConvergenceQuantityType {
    /// 能量变化百分比
    DeltaEnergy,
    /// MoM 面网格三角形数量
    NumTriangles,
    /// FEM 四面体数量（DC 提取时）
    NumTetrahedra,
    /// RLCG 矩阵极值变化（R_max, L_max, C_max, G_max）
    RlcgMaxR,
    RlcgMaxL,
    RlcgMaxC,
    RlcgMaxG,
}

/// 从 RLCG 矩阵数据中提取指定元素的频率曲线
pub fn extract_rlcg_trace(
    data: &RlcgMatrixData,
    param_type: RlcgParamType,
    terminal_i: usize,
    terminal_j: usize,
) -> Vec<PlotPoint> {
    let series = match param_type {
        RlcgParamType::R => &data.matrices.r,
        RlcgParamType::L => &data.matrices.l,
        RlcgParamType::C => &data.matrices.c,
        RlcgParamType::G => &data.matrices.g,
    };
    if let Some(s) = series {
        s.data_per_frequency.iter().map(|fm| {
            PlotPoint {
                x: fm.frequency,
                y: fm.matrix[terminal_i][terminal_j],
            }
        }).collect()
    } else {
        vec![]
    }
}
```

**交互功能**：

| 功能 | 说明 | 实现方式 |
|------|------|---------|
| 矩阵类型切换 | R/L/C/G 标签页或下拉菜单 | egui ComboBox |
| 端子对选择 | 选择要显示的端子对 (i,j) | 复选框矩阵或列表 |
| 自阻/互阻分离 | 对角线（自阻）与非对角线（互阻）分图显示 | 过滤 i==j / i!=j |
| DC 值标注 | 在 f=0 处标注 DC 电阻/电感值 | Marker at x=0 |
| Log/Linear 切换 | 频率轴对数/线性切换 | X 轴 scale 控制 |
| 耦合系数 | 显示 k = M / sqrt(L_ii * L_jj) | 计算派生量 |
| 导出 CSV | 矩阵数据表格导出 | CSV writer |

### 3.10 Q3D RLCG 矩阵数据表格

以矩阵表格形式显示指定频率下的 RLCG 数值。参考 [Viewing Matrix Data in Q3D](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/ViewingMatrixDatainQ3D.htm)。

**表格布局**：

```
┌──────────────────────────────────────────────────────────┐
│  Matrix Type: [R ▼]    Frequency: [1.0 GHz ▼]           │
├──────────────┬──────────┬──────────┬──────────┬──────────┤
│              │ Sig1:T1  │ Sig1:T2  │ Sig2:T3  │ Sig2:T4  │
├──────────────┼──────────┼──────────┼──────────┼──────────┤
│ Sig1:T1      │  0.285   │  0.012   │  0.005   │  0.002   │
│ Sig1:T2      │  0.012   │  0.285   │  0.002   │  0.005   │
│ Sig2:T3      │  0.005   │  0.002   │  0.268   │  0.010   │
│ Sig2:T4      │  0.002   │  0.005   │  0.010   │  0.268   │
└──────────────┴──────────┴──────────┴──────────┴──────────┘
```

**表格交互**：

| 功能 | 说明 |
|------|------|
| 频率选择 | 下拉菜单或滑块选择当前频率点 |
| 矩阵类型切换 | R / L / C / G 四个标签页 |
| 单元格着色 | 矩阵元素按值大小热力图着色（对角线自身颜色深，耦合项颜色浅） |
| 对称性高亮 | 高亮显示主对角线，区分自阻/互阻 |
| DC 对比 | 可选显示 DC 值与 AC 值的对比（△R = R_ac - R_dc） |
| 排序 | 按行/列排序，或按耦合强度排序 |
| 导出 | 导出为 CSV（完整矩阵 + 所有频率） |

**Rust 接口**：

```rust
pub struct RlcgMatrixTablePanel {
    /// 当前选择的矩阵类型
    pub param_type: RlcgParamType,
    /// 当前选择的频率索引
    pub frequency_index: usize,
    /// RLCG 数据源
    pub data: Arc<RlcgMatrixData>,
    /// 热力图色标
    pub heatmap_colormap: ColormapType,
    /// 显示选项
    pub options: MatrixTableOptions,
}

#[derive(Debug, Clone)]
pub struct MatrixTableOptions {
    /// 是否启用热力图着色
    pub enable_heatmap: bool,
    /// 是否显示 DC 值对比
    pub show_dc_comparison: bool,
    /// 数值精度（小数位数）
    pub decimal_places: u32,
    /// 是否显示单位
    pub show_unit: bool,
}
```

### 3.11 Q3D 收敛曲线

Q3D 准静态求解的收敛曲线与 HFSS 类似，但使用不同的收敛指标。

**图表内容**：

| 子图 | X 轴 | Y 轴 | 说明 |
|------|------|------|------|
| Delta Energy 收敛 | Pass number | Delta Energy (%) | 带目标线（target_max_delta_energy） |
| 网格增长 | Pass number | Num triangles | MoM 面网格增长（柱状图） |
| RLCG 值稳定性 | Pass number | R_max / L_max / C_max | 矩阵极值随自适应轮次的变化 |
| 求解时间 | Pass number | Solve time (sec) | 分阶段堆叠 |

> **与 HFSS 收敛曲线的差异**：HFSS 使用 `Max Delta S`（S 参数变化）作为收敛指标，Q3D 使用 `Delta Energy`（能量变化百分比）。Q3D 还额外显示 RLCG 矩阵极值的稳定趋势，这是 Q3D 特有的收敛可视化——当矩阵值不再随网格加密而显著变化时，即可认为结果已收敛。

---

## 4. 3D 场可视化（Field Overlays）

### 4.1 设计思想

参考 HFSS Field Overlays：用户在 3D 模型上叠加显示电磁场分布。EMStudio 将场数据从 `.emsfld` 按需加载，映射为 GPU 可渲染的顶点颜色或纹理，通过 wgpu 管线叠加到几何模型上。

### 4.2 场叠加类型

```rust
/// 场叠加定义（对应工程文件中的 field_overlays[]）
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldOverlay {
    pub name: String,
    pub overlay_type: FieldOverlayType,
    pub quantity: FieldQuantityConfig,
    pub visual: FieldVisualConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FieldOverlayType {
    /// 表面云图：在几何表面上用颜色映射场幅度
    SurfacePlot {
        /// 显示在哪些面上（None = 所有外表面）
        surfaces: Option<Vec<String>>,
    },
    /// 体积云图切面：用一个平面截取体积场数据
    SlicePlot {
        /// 切面定义
        plane: SlicePlane,
    },
    /// 矢量箭头图：在采样点显示场矢量方向和幅度
    VectorPlot {
        /// 箭头密度（采样间隔）
        spacing_mm: f64,
        /// 箭头缩放因子
        arrow_scale: f64,
    },
    /// 等值面：场幅度等于特定值的 3D 等值面
    IsosurfacePlot {
        /// 等值面的值
        iso_value: f64,
    },
    /// 动画：相位从 0° 到 360° 扫描，显示场的时域变化
    AnimatedPlot {
        /// 相位步进（度）
        phase_step_deg: f64,
        /// 每秒帧数
        fps: f64,
    },
}

/// 切面定义
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SlicePlane {
    /// 平行于 XY 平面，Z = value
    XY { z_mm: f64 },
    /// 平行于 XZ 平面，Y = value
    XZ { y_mm: f64 },
    /// 平行于 YZ 平面，X = value
    YZ { x_mm: f64 },
    /// 任意平面：法向量 + 平面上一点
    Arbitrary {
        normal: [f64; 3],
        point: [f64; 3],
    },
}
```

### 4.3 场量配置

```rust
/// 场量配置
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldQuantityConfig {
    /// 场类型
    pub field_type: FieldType,
    /// 分量选择
    pub component: FieldComponent,
    /// 频率选择
    pub frequency: String,
    /// 相位（度），用于复数场 → 实数场映射
    pub phase_deg: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FieldType {
    Electric,       // E-field (V/m) — HFSS & Q3D
    Magnetic,       // H-field (A/m) — HFSS only
    Current,        // J-field (A/m²) — 表面电流 (HFSS & Q3D)
    Poynting,       // S-field (W/m²) — 功率流密度 (HFSS only)
    VolumeCurrent,  // Jvol-field (A/m²) — 体电流密度 (Q3D: DC/AC 电流)
    ChargeDistribution, // ρ-field (C/m²) — 电荷面密度 (Q3D: 电容提取)
    OhmicLoss,      // P_loss (W/m³) — 欧姆损耗密度 (Q3D: 电阻提取)
}

/// Q3D 场量说明：
/// - `VolumeCurrent`：Q3D AC/DC 电阻和电感提取时的体电流密度分布，
///   用于识别电流拥挤、趋肤效应和邻近效应区域。
/// - `ChargeDistribution`：Q3D 电容提取时的电荷面密度分布，
///   用于分析导体间的电场耦合路径。
/// - `OhmicLoss`：Q3D 的欧姆损耗密度，用于识别热点区域。
/// - `Electric`：Q3D 电容/电导提取时的电场分布。
/// - `Current`：Q3D 面电流密度，用于薄导体的电流可视化。

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FieldComponent {
    /// 矢量幅度 |F| = sqrt(|Fx|² + |Fy|² + |Fz|²)
    Magnitude,
    /// X 分量
    X,
    /// Y 分量
    Y,
    /// Z 分量
    Z,
    /// 矢量本身（用于矢量箭头图）
    Vector,
}
```

### 4.4 色标系统（Colormap）

```rust
/// 色标可视化配置
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldVisualConfig {
    /// 色标类型
    pub colormap: ColormapType,
    /// 值域范围
    pub range: ValueRange,
    /// 色标刻度
    pub scale: ScaleType,
    /// 透明度（0.0 ~ 1.0）
    pub opacity: f64,
    /// 是否显示色标条
    pub show_legend: bool,
    /// 色标条位置
    pub legend_position: LegendPosition,
    /// 超出范围的处理
    pub out_of_range: OutOfRangePolicy,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ColormapType {
    /// 彩虹色标（HFSS 默认）：蓝 → 青 → 绿 → 黄 → 红
    Rainbow,
    /// 热力色标：黑 → 红 → 黄 → 白
    Hot,
    /// 冷暖色标：蓝 → 白 → 红（适合正负值）
    CoolWarm,
    /// 灰度
    Grayscale,
    /// Viridis（感知均匀，色盲友好）
    Viridis,
    /// Plasma
    Plasma,
    /// 自定义色标
    Custom { stops: Vec<ColorStop> },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColorStop {
    pub position: f64,  // 0.0 ~ 1.0
    pub color: [u8; 4], // RGBA
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ValueRange {
    /// 自动：使用数据的 min/max
    Auto,
    /// 手动指定范围
    Manual { min: f64, max: f64 },
    /// 对称：[-max_abs, +max_abs]（适合矢量分量）
    Symmetric,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ScaleType {
    /// 线性映射
    Linear,
    /// 对数映射（dB 刻度）
    Logarithmic,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LegendPosition {
    Left,
    Right,
    Top,
    Bottom,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum OutOfRangePolicy {
    /// 钳位到最近的端点颜色
    Clamp,
    /// 显示为透明
    Transparent,
    /// 显示为指定颜色
    FixedColor([u8; 4]),
}
```

### 4.5 GPU 渲染管线

#### 4.5.1 表面云图管线（Surface Colormap Pipeline）

场幅度数据映射为顶点颜色，通过 wgpu 渲染管线显示在网格表面上。

```
┌──────────────────────────────────────────────────────────┐
│                Surface Colormap Pipeline                   │
│                                                           │
│  CPU Side:                                                │
│  ┌──────────┐    ┌────────────┐    ┌─────────────────┐   │
│  │ .emsfld  │───▶│QuantityExpr│───▶│ per-vertex      │   │
│  │ (mmap)   │    │ evaluate() │    │ scalar value    │   │
│  └──────────┘    └────────────┘    └───────┬─────────┘   │
│                                            │              │
│  ┌──────────┐                     ┌────────▼──────────┐  │
│  │ .msh     │────────────────────▶│ Vertex Buffer     │  │
│  │ (nodes)  │                     │ [pos.xyz, value]  │  │
│  └──────────┘                     └────────┬──────────┘  │
│                                            │              │
│  GPU Side:                                 │              │
│  ┌─────────────────────────────────────────▼───────┐     │
│  │ Vertex Shader                                    │     │
│  │  - transform position by MVP matrix              │     │
│  │  - pass scalar value to fragment                 │     │
│  └─────────────────────────────────┬───────────────┘     │
│                                    │                      │
│  ┌─────────────────────────────────▼───────────────┐     │
│  │ Fragment Shader                                  │     │
│  │  - sample 1D colormap texture by scalar value    │     │
│  │  - apply lighting (optional)                     │     │
│  │  - output RGBA                                   │     │
│  └─────────────────────────────────────────────────┘     │
└──────────────────────────────────────────────────────────┘
```

**Vertex Buffer 布局**：

```rust
#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
pub struct FieldVertex {
    pub position: [f32; 3],    // 节点坐标
    pub normal: [f32; 3],      // 法向量（用于光照）
    pub field_value: f32,      // 归一化后的场值 [0, 1]
}
```

**WGSL Shader**：

```wgsl
// ---- Vertex Shader ----
struct Uniforms {
    mvp: mat4x4<f32>,
    model: mat4x4<f32>,
    light_dir: vec3<f32>,
    ambient: f32,
    value_min: f32,
    value_max: f32,
    opacity: f32,
    _pad: f32,
};

@group(0) @binding(0) var<uniform> u: Uniforms;
@group(0) @binding(1) var colormap_texture: texture_1d<f32>;
@group(0) @binding(2) var colormap_sampler: sampler;

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) field_value: f32,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) world_normal: vec3<f32>,
    @location(1) normalized_value: f32,
};

@vertex
fn vs_main(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    out.clip_position = u.mvp * vec4<f32>(in.position, 1.0);
    out.world_normal = (u.model * vec4<f32>(in.normal, 0.0)).xyz;
    // 归一化场值到 [0, 1]
    out.normalized_value = clamp(
        (in.field_value - u.value_min) / (u.value_max - u.value_min),
        0.0, 1.0
    );
    return out;
}

// ---- Fragment Shader ----
@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    // 从 1D 色标纹理采样
    let color = textureSample(colormap_texture, colormap_sampler, in.normalized_value);

    // 简单 Lambertian 光照
    let n = normalize(in.world_normal);
    let diffuse = max(dot(n, normalize(u.light_dir)), 0.0);
    let lighting = u.ambient + (1.0 - u.ambient) * diffuse;

    return vec4<f32>(color.rgb * lighting, color.a * u.opacity);
}
```

#### 4.5.2 矢量箭头管线（Vector Arrow Pipeline）

将场矢量显示为 3D 箭头。

**实现方案**：使用 instanced rendering，每个箭头是一个实例。

```rust
/// 箭头实例数据
#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
pub struct ArrowInstance {
    pub position: [f32; 3],    // 箭头起点
    pub direction: [f32; 3],   // 箭头方向（归一化）
    pub magnitude: f32,        // 幅度（决定箭头长度和颜色）
    pub _pad: f32,
}
```

**箭头几何体**：预生成一个单位箭头 Mesh（圆柱体 + 锥体），通过 instance transform 缩放和旋转。

#### 4.5.3 等值面管线（Isosurface Pipeline）

从体积场数据中提取等值面。

**算法**：Marching Tetrahedra（四面体网格天然适配）

```rust
pub struct IsosurfaceExtractor;

impl IsosurfaceExtractor {
    /// 从四面体网格 + 节点场值中提取等值面
    /// 返回三角形网格
    pub fn extract(
        mesh: &MshMesh,
        field_values: &[f64],  // 每节点一个标量值
        iso_value: f64,
    ) -> TriangleMesh {
        // Marching Tetrahedra 算法：
        // 对每个四面体，根据 4 个顶点的场值与 iso_value 的大小关系
        // 确定等值面与四面体边的交点，生成 0~2 个三角形
        todo!()
    }
}
```

### 4.6 切面可视化（Slice Plot）

在任意平面上显示体积场数据的截面分布。

**数据流**：

```
.emsfld (volume field data, per-node)
       │
       ▼
 SlicePlane 定义截面
       │
       ▼
 mesh-plane intersection → 截面多边形网格
       │
       ▼
 插值每个截面顶点的场值 (barycentric interpolation)
       │
       ▼
 colormap → 截面三角形顶点颜色
       │
       ▼
 Surface Colormap Pipeline 渲染
```

```rust
pub struct SliceMeshGenerator;

impl SliceMeshGenerator {
    /// 用平面截取四面体网格，生成截面三角形网格
    /// 每个截面顶点插值出场值
    pub fn generate(
        mesh: &MshMesh,
        field_values: &[f64],
        plane: &SlicePlane,
    ) -> SliceMesh {
        // 遍历每个四面体，计算与平面的交线
        // 交线上的点通过重心坐标插值场值
        todo!()
    }
}

pub struct SliceMesh {
    pub vertices: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    pub field_values: Vec<f32>,
    pub indices: Vec<u32>,
}
```

### 4.7 相位动画

对于复数场数据，支持相位扫描动画，展示场的时域变化：

```
E(t) = Re[ E_complex * e^{j * phase} ]
```

**实现**：

```rust
pub struct PhaseAnimator {
    /// 当前相位（度）
    pub current_phase_deg: f64,
    /// 相位步进
    pub phase_step_deg: f64,
    /// 是否播放中
    pub playing: bool,
    /// 帧率控制
    pub fps: f64,
    /// 上一帧时间戳
    last_frame: std::time::Instant,
}

impl PhaseAnimator {
    /// 计算当前相位下的实数场值
    /// E_real = Re(E) * cos(phase) - Im(E) * sin(phase)
    pub fn evaluate(
        &self,
        field_real: &[f64],
        field_imag: &[f64],
    ) -> Vec<f64> {
        let phase_rad = self.current_phase_deg.to_radians();
        let cos_p = phase_rad.cos();
        let sin_p = phase_rad.sin();
        field_real.iter().zip(field_imag.iter())
            .map(|(re, im)| re * cos_p - im * sin_p)
            .collect()
    }

    /// 推进一帧
    pub fn tick(&mut self) {
        if self.playing {
            let now = std::time::Instant::now();
            if now.duration_since(self.last_frame).as_secs_f64() >= 1.0 / self.fps {
                self.current_phase_deg =
                    (self.current_phase_deg + self.phase_step_deg) % 360.0;
                self.last_frame = now;
            }
        }
    }
}
```

---

## 5. 3D 远场方向图

### 5.1 3D 辐射方向图表面

将远场增益数据渲染为 3D 球面变形表面：每个方向 (θ, φ) 上的半径 = f(Gain)，颜色也映射增益值。

**球面 → 变形表面**：

```rust
/// 将远场增益数据生成 3D 方向图 Mesh
pub fn generate_3d_pattern_mesh(
    far_field: &FarFieldData,
    quantity: &str,  // "GainTotal", "GainTheta", etc.
    scale: PatternScale,
) -> PatternMesh {
    let theta_pts = far_field.theta.num_points;
    let phi_pts = far_field.phi.num_points;
    let gain_data = &far_field.derived_quantities[quantity].data;

    let mut vertices = Vec::new();
    let mut indices = Vec::new();

    for phi_idx in 0..phi_pts {
        for theta_idx in 0..theta_pts {
            let theta_rad = (far_field.theta.start_deg
                + theta_idx as f64 * far_field.theta.step_deg)
                .to_radians();
            let phi_rad = (far_field.phi.start_deg
                + phi_idx as f64 * far_field.phi.step_deg)
                .to_radians();

            // 增益值 → 半径
            let gain_dbi = gain_data[phi_idx][theta_idx];
            let radius = scale.gain_to_radius(gain_dbi);

            // 球坐标 → 笛卡尔坐标
            let x = radius * theta_rad.sin() * phi_rad.cos();
            let y = radius * theta_rad.sin() * phi_rad.sin();
            let z = radius * theta_rad.cos();

            vertices.push(PatternVertex {
                position: [x as f32, y as f32, z as f32],
                gain_normalized: scale.normalize(gain_dbi) as f32,
            });
        }
    }

    // 生成三角形索引（相邻 theta/phi 点连接）
    for phi_idx in 0..(phi_pts - 1) {
        for theta_idx in 0..(theta_pts - 1) {
            let i00 = (phi_idx * theta_pts + theta_idx) as u32;
            let i01 = i00 + 1;
            let i10 = ((phi_idx + 1) * theta_pts + theta_idx) as u32;
            let i11 = i10 + 1;
            indices.extend_from_slice(&[i00, i01, i10, i01, i11, i10]);
        }
    }

    PatternMesh { vertices, indices }
}

pub enum PatternScale {
    /// 线性：radius = max(gain_linear, floor)
    Linear { floor_dbi: f64 },
    /// 对数：radius = max(gain_dbi - floor, 0) / (peak - floor)
    Logarithmic { floor_dbi: f64 },
}

impl PatternScale {
    pub fn gain_to_radius(&self, gain_dbi: f64) -> f64 {
        match self {
            Self::Linear { floor_dbi } => {
                let gain_linear = 10f64.powf(gain_dbi / 10.0);
                let floor_linear = 10f64.powf(floor_dbi / 10.0);
                (gain_linear - floor_linear).max(0.0)
            }
            Self::Logarithmic { floor_dbi } => {
                (gain_dbi - floor_dbi).max(0.0) / (0.0 - floor_dbi).abs()
            }
        }
    }

    pub fn normalize(&self, gain_dbi: f64) -> f64 {
        // 归一化到 [0, 1] 用于色标映射
        self.gain_to_radius(gain_dbi)
    }
}
```

### 5.2 3D 方向图交互

| 功能 | 说明 |
|------|------|
| 旋转 | 轨迹球旋转，观察不同角度 |
| 缩放 | 滚轮缩放 |
| 切面切换 | 快捷键切换 E/H 面高亮 |
| 增益数值 | 鼠标悬停显示 (θ, φ, Gain) |
| 坐标轴 | 显示 X/Y/Z 参考轴 |
| 频率动画 | 滑块切换频率，观察方向图随频率的变化 |
| 叠加模型 | 半透明显示天线几何模型作为参考 |

---

## 6. 3D 网格可视化

### 6.1 网格显示模式

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MeshDisplayMode {
    /// 线框模式：仅显示四面体棱边
    Wireframe {
        line_color: [u8; 4],
        line_width: f32,
    },
    /// 实体模式：显示外表面三角形
    Solid {
        face_color: [u8; 4],
        show_edges: bool,
        edge_color: [u8; 4],
    },
    /// 质量着色：按四面体质量指标着色
    QualityColor {
        metric: MeshQualityMetric,
        colormap: ColormapType,
        range: ValueRange,
    },
    /// 材料区域着色：按材料分色
    MaterialColor,
    /// 透明 + 线框：外表面半透明 + 内部线框
    XRay {
        surface_opacity: f32,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MeshQualityMetric {
    /// 纵横比 (≥1, 越小越好)
    AspectRatio,
    /// 体积
    Volume,
    /// 最小二面角
    MinDihedralAngle,
}
```

### 6.2 网格表面提取

四面体体网格的外表面（用于实体/云图显示）通过边界面提取：

```rust
/// 从四面体网格中提取外表面三角形
/// 原理：一个三角面如果只属于一个四面体，则为外表面
pub fn extract_surface_triangles(mesh: &MshMesh) -> Vec<SurfaceTriangle> {
    // 1. 遍历所有四面体，列出其 4 个三角面
    // 2. 用 HashMap 统计每个三角面出现次数（需排序顶点 tag 以去重）
    // 3. 出现一次的面即为外表面
    // 4. 计算每个三角形的法向量
    todo!()
}
```

---

## 7. 交互系统

### 7.1 相机控制

```rust
pub struct OrbitCamera {
    /// 观察目标点
    pub target: [f64; 3],
    /// 距目标的距离
    pub distance: f64,
    /// 方位角（绕 Y 轴旋转，度）
    pub azimuth_deg: f64,
    /// 仰角（度，-90 ~ 90）
    pub elevation_deg: f64,
    /// 垂直视场角（度）
    pub fov_deg: f64,
    /// 近裁剪面
    pub near: f64,
    /// 远裁剪面
    pub far: f64,
}

impl OrbitCamera {
    /// 计算 View 矩阵
    pub fn view_matrix(&self) -> [[f64; 4]; 4] { todo!() }

    /// 计算 Projection 矩阵
    pub fn projection_matrix(&self, aspect_ratio: f64) -> [[f64; 4]; 4] { todo!() }

    /// 鼠标拖拽 → 旋转
    pub fn orbit(&mut self, delta_x: f64, delta_y: f64) {
        self.azimuth_deg += delta_x * 0.5;
        self.elevation_deg = (self.elevation_deg + delta_y * 0.5).clamp(-89.0, 89.0);
    }

    /// 鼠标中键拖拽 → 平移
    pub fn pan(&mut self, delta_x: f64, delta_y: f64) { todo!() }

    /// 滚轮 → 缩放
    pub fn zoom(&mut self, delta: f64) {
        self.distance *= (1.0 - delta * 0.1).max(0.01);
    }

    /// 适配到包围盒
    pub fn fit_to_bounds(&mut self, aabb_min: [f64; 3], aabb_max: [f64; 3]) {
        // 计算包围盒中心作为 target
        // 计算合适的 distance 使得包围盒完全可见
        todo!()
    }

    /// 预设视角
    pub fn set_view_preset(&mut self, preset: ViewPreset) {
        match preset {
            ViewPreset::Front  => { self.azimuth_deg = 0.0;   self.elevation_deg = 0.0; }
            ViewPreset::Back   => { self.azimuth_deg = 180.0; self.elevation_deg = 0.0; }
            ViewPreset::Left   => { self.azimuth_deg = -90.0; self.elevation_deg = 0.0; }
            ViewPreset::Right  => { self.azimuth_deg = 90.0;  self.elevation_deg = 0.0; }
            ViewPreset::Top    => { self.azimuth_deg = 0.0;   self.elevation_deg = 89.0; }
            ViewPreset::Bottom => { self.azimuth_deg = 0.0;   self.elevation_deg = -89.0; }
            ViewPreset::Iso    => { self.azimuth_deg = 45.0;  self.elevation_deg = 35.0; }
        }
    }
}

pub enum ViewPreset {
    Front, Back, Left, Right, Top, Bottom, Iso,
}
```

### 7.2 GPU Picking（拾取）

支持鼠标点击选择几何对象或查询场值。

```rust
/// 基于颜色 ID 的 GPU Picking
pub struct PickingSystem {
    /// 离屏渲染目标（每个像素存储 object_id）
    pick_texture: wgpu::Texture,
    /// 回读缓冲区
    readback_buffer: wgpu::Buffer,
}

impl PickingSystem {
    /// 渲染一帧 picking pass（每个三角形/对象用唯一颜色编码）
    pub fn render_pick_pass(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        scene: &SceneData,
        camera: &OrbitCamera,
    ) { todo!() }

    /// 查询鼠标位置的 object_id
    pub fn query(&self, x: u32, y: u32) -> Option<PickResult> { todo!() }
}

pub enum PickResult {
    GeometryObject { object_id: u64 },
    MeshElement { element_tag: u64 },
    FieldValue { node_tag: u64, value: f64 },
}
```

### 7.3 Probe（场值探针）

参考 HFSS 的 Field Calculator 探针功能：在 3D 视图中点击任意位置，查询该点的场值。

```rust
pub struct FieldProbe {
    /// 探针位置（世界坐标）
    pub position_mm: [f64; 3],
    /// 探测结果
    pub result: Option<ProbeResult>,
}

pub struct ProbeResult {
    /// 所在四面体 tag
    pub element_tag: u64,
    /// 重心坐标（用于插值验证）
    pub barycentric: [f64; 4],
    /// 场值
    pub field: ProbeFieldValue,
}

pub struct ProbeFieldValue {
    pub e_field: [Complex64; 3],   // (Ex, Ey, Ez)
    pub e_magnitude: f64,          // |E|
    pub h_field: [Complex64; 3],   // (Hx, Hy, Hz)
    pub h_magnitude: f64,          // |H|
}
```

---

## 8. 数据加载与缓存

### 8.1 ResultDataStore

统一管理所有结果数据的加载和缓存。

```rust
pub struct ResultDataStore {
    /// 结果目录路径
    result_dir: PathBuf,
    /// JSON 文件缓存
    json_cache: HashMap<String, serde_json::Value>,
    /// 场数据句柄（mmap，按需加载单频点）
    field_handles: HashMap<String, FieldFileHandle>,
    /// 网格缓存
    mesh_cache: Option<MshMesh>,
    /// Touchstone 缓存
    touchstone_cache: HashMap<String, TouchstoneData>,
}

impl ResultDataStore {
    /// 打开结果目录
    pub fn open(result_dir: &Path) -> Result<Self, DataError> { todo!() }

    /// 加载 S 参数（优先从 Touchstone，回退到 JSON）
    pub fn load_s_parameters(
        &mut self,
        setup: &str,
        sweep: &str,
    ) -> Result<&SParameterData, DataError> { todo!() }

    /// 加载远场数据
    pub fn load_far_field(
        &mut self,
        setup: &str,
        frequency: &str,
    ) -> Result<&FarFieldData, DataError> { todo!() }

    /// 加载收敛历史
    pub fn load_convergence(
        &mut self,
        setup: &str,
    ) -> Result<&ConvergenceHistory, DataError> { todo!() }

    /// 加载场数据（指定频点，按需 mmap）
    pub fn load_field(
        &mut self,
        field_file: &str,
        freq_index: usize,
    ) -> Result<FieldSlice, DataError> { todo!() }

    /// 加载网格
    pub fn load_mesh(&mut self) -> Result<&MshMesh, DataError> { todo!() }
}

/// 单频点场数据切片（零拷贝引用 mmap 区域）
pub struct FieldSlice<'a> {
    pub frequency_hz: f64,
    pub num_nodes: usize,
    pub num_components: usize,
    /// 实际数据：[node_0_re_x, node_0_im_x, node_0_re_y, ...] 或 f32 版本
    pub data_f64: Option<&'a [f64]>,
    pub data_f32: Option<&'a [f32]>,
}
```

### 8.2 场数据内存映射

`.emsfld` 文件通过 `mmap` 映射到内存，实现零拷贝按需加载：

```rust
pub struct FieldFileHandle {
    /// 内存映射
    mmap: memmap2::Mmap,
    /// 文件头（解析后的副本）
    header: EmsFldHeader,
    /// 频率表
    frequencies: Vec<f64>,
    /// 索引表
    block_index: Vec<FieldBlockInfo>,
}

impl FieldFileHandle {
    pub fn open(path: &Path) -> Result<Self, DataError> {
        let file = std::fs::File::open(path)?;
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        // 解析 header、频率表、索引表
        todo!()
    }

    /// 获取指定频点的场数据切片（零拷贝）
    pub fn slice(&self, freq_index: usize) -> Result<FieldSlice<'_>, DataError> {
        let block_info = &self.block_index[freq_index];
        let offset = block_info.offset as usize;
        let size = block_info.size_bytes as usize;
        let raw_bytes = &self.mmap[offset..offset + size];
        // 根据 data_type 转换为 &[f64] 或 &[f32]
        todo!()
    }
}
```

### 8.3 缓存策略

| 数据类型 | 缓存策略 | 说明 |
|---------|---------|------|
| JSON 元数据 | 全量缓存 | KB 级别，直接 deserialize 并持有 |
| Touchstone | 全量缓存 | 通常 < 10 MB |
| 网格 (.msh) | 全量缓存 | 仅最终网格，通常 < 100 MB |
| 场数据 (.emsfld) | mmap + LRU | 文件 mmap，频点级 GPU Buffer LRU 缓存 |
| GPU Buffer | LRU | 最多缓存 N 个频点的 Vertex Buffer，超出时回收最久未用 |
| 色标纹理 | 持久 | 1D 纹理，几乎不占 VRAM |

```rust
pub struct GpuBufferCache {
    /// 最大缓存的频点数
    max_entries: usize,
    /// LRU 缓存：key = (field_file, freq_index)
    cache: lru::LruCache<(String, usize), wgpu::Buffer>,
}
```

---

## 9. UI 面板集成

### 9.1 Tab 系统扩展

在现有的 `CenterTab` 枚举上扩展可视化 Tab：

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
enum CenterTab {
    /// 3D 建模视图
    Modeling,
    /// 2D 报告图表
    Report(ReportId),
    /// 3D 场叠加视图
    FieldOverlay(OverlayId),
    /// 3D 远场方向图
    FarField3D(FarFieldViewId),
    /// 求解器日志
    Log,
}
```

### 9.2 右侧属性面板

当选中一个可视化 Tab 时，右侧面板显示对应的属性编辑器：

| 当前 Tab | 右侧面板内容 |
|---------|------------|
| Report(S-Param) | Trace 列表、添加/删除 Trace、量表达式编辑、坐标轴范围、Marker 管理 |
| Report(Polar) | 切面选择、增益类型、dB 范围 |
| Report(Smith) | 参考阻抗、归一化选项 |
| FieldOverlay | 场类型/分量选择、频率滑块、色标类型/范围、切面位置滑块、透明度 |
| FarField3D | 增益类型、缩放模式、底噪、频率选择 |
| Modeling | 对象属性、材料、变换 |

### 9.3 属性面板示例：场叠加

```rust
pub fn field_overlay_properties(
    ui: &mut egui::Ui,
    overlay: &mut FieldOverlay,
    available_frequencies: &[String],
) {
    ui.heading("Field Overlay Properties");
    ui.separator();

    // 场类型选择
    ui.horizontal(|ui| {
        ui.label("Field:");
        egui::ComboBox::from_id_salt("field_type")
            .selected_text(format!("{:?}", overlay.quantity.field_type))
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut overlay.quantity.field_type,
                    FieldType::Electric, "E-Field",
                );
                ui.selectable_value(
                    &mut overlay.quantity.field_type,
                    FieldType::Magnetic, "H-Field",
                );
                ui.selectable_value(
                    &mut overlay.quantity.field_type,
                    FieldType::Current, "J-Surface",
                );
            });
    });

    // 分量选择
    ui.horizontal(|ui| {
        ui.label("Component:");
        egui::ComboBox::from_id_salt("component")
            .selected_text(format!("{:?}", overlay.quantity.component))
            .show_ui(ui, |ui| {
                for comp in &[
                    FieldComponent::Magnitude,
                    FieldComponent::X,
                    FieldComponent::Y,
                    FieldComponent::Z,
                ] {
                    ui.selectable_value(
                        &mut overlay.quantity.component,
                        comp.clone(),
                        format!("{comp:?}"),
                    );
                }
            });
    });

    // 频率选择
    ui.horizontal(|ui| {
        ui.label("Frequency:");
        egui::ComboBox::from_id_salt("frequency")
            .selected_text(&overlay.quantity.frequency)
            .show_ui(ui, |ui| {
                for freq in available_frequencies {
                    ui.selectable_value(
                        &mut overlay.quantity.frequency,
                        freq.clone(),
                        freq,
                    );
                }
            });
    });

    // 相位滑块
    ui.add(
        egui::Slider::new(&mut overlay.quantity.phase_deg, 0.0..=360.0)
            .text("Phase (°)")
    );

    ui.separator();
    ui.heading("Color Map");

    // 色标类型
    ui.horizontal(|ui| {
        ui.label("Colormap:");
        egui::ComboBox::from_id_salt("colormap")
            .selected_text(format!("{:?}", overlay.visual.colormap))
            .show_ui(ui, |ui| {
                for cmap in &[
                    ColormapType::Rainbow,
                    ColormapType::Hot,
                    ColormapType::CoolWarm,
                    ColormapType::Viridis,
                    ColormapType::Grayscale,
                ] {
                    ui.selectable_value(
                        &mut overlay.visual.colormap,
                        cmap.clone(),
                        format!("{cmap:?}"),
                    );
                }
            });
    });

    // 值域范围
    match &mut overlay.visual.range {
        ValueRange::Auto => {
            if ui.button("Switch to Manual Range").clicked() {
                overlay.visual.range = ValueRange::Manual {
                    min: 0.0,
                    max: 1.0,
                };
            }
        }
        ValueRange::Manual { min, max } => {
            ui.horizontal(|ui| {
                ui.label("Min:");
                ui.add(egui::DragValue::new(min).speed(0.1));
                ui.label("Max:");
                ui.add(egui::DragValue::new(max).speed(0.1));
            });
            if ui.button("Switch to Auto Range").clicked() {
                overlay.visual.range = ValueRange::Auto;
            }
        }
        _ => {}
    }

    // 透明度
    ui.add(
        egui::Slider::new(&mut overlay.visual.opacity, 0.0..=1.0)
            .text("Opacity")
    );
}
```

---

## 10. 导出功能

### 10.1 导出类型

| 导出目标 | 格式 | 数据来源 | 说明 |
|---------|------|---------|------|
| 报告数据 | CSV | ReportPanel traces | 频率 + 各 Trace 值，逗号分隔 |
| 报告截图 | PNG | ReportPanel | egui 截图或自渲染 |
| 场图截图 | PNG | SceneRenderer | wgpu offscreen render → readback → PNG |
| 3D 场数据 | VTK | .msh + .emsfld | ParaView 兼容的 VTK Unstructured Grid |
| S 参数 | Touchstone | SParameterData | 调用 emstudio-touchstone writer |
| 方向图数据 | CSV | FarFieldData | θ, φ, GainTotal, GainTheta, GainPhi |
| RLCG 矩阵 | CSV | RlcgMatrixData | Q3D：全频率完整 RLCG 矩阵 |
| 等效电路 | SPICE | RlcgMatrixData | Q3D：SPICE 网表导出 |
| RLCG → S 参数 | Touchstone | RlcgMatrixData | Q3D：RLCG 矩阵转换为 S 参数 |

### 10.2 CSV 导出

```rust
/// S 参数报告 CSV 导出
pub fn export_s_param_csv(
    s_data: &SParameterData,
    traces: &[Trace],
    path: &Path,
) -> Result<(), ExportError> {
    let mut writer = csv::Writer::from_path(path)?;

    // 写表头
    let mut header = vec!["Frequency (GHz)".to_string()];
    for trace in traces {
        header.push(trace.name.clone());
    }
    writer.write_record(&header)?;

    // 写数据行
    for (i, &freq) in s_data.frequencies.iter().enumerate() {
        let mut row = vec![freq.to_string()];
        for trace in traces {
            let value = trace.expression.evaluate_at_index(s_data, i);
            row.push(format!("{:.6}", value));
        }
        writer.write_record(&row)?;
    }

    writer.flush()?;
    Ok(())
}
```

### 10.3 PNG 截图

```rust
/// 3D 视图截图
pub fn capture_screenshot(
    renderer: &SceneRenderer,
    width: u32,
    height: u32,
) -> Result<Vec<u8>, RenderError> {
    // 1. 创建离屏纹理
    // 2. 渲染一帧到离屏纹理
    // 3. 将纹理数据回读到 CPU buffer
    // 4. 编码为 PNG
    todo!()
}
```

---

## 11. 实现路线

### Phase 1：2D 报告基础

| 任务 | 模块 | 依赖 |
|------|------|------|
| ResultDataStore + JSON/Touchstone 加载 | `emstudio-domain` | 现有 touchstone crate |
| S 参数矩形图（egui_plot） | `emstudio-components` | ResultDataStore |
| Trace/QuantityExpression 系统 | `emstudio-domain` | — |
| 收敛曲线 | `emstudio-components` | ResultDataStore |
| ReportPanel 集成到 DockArea | `emstudio-app` | ReportPanel |

### Phase 2：2D 报告进阶

| 任务 | 模块 | 依赖 |
|------|------|------|
| Smith 圆图 | `emstudio-components` | Phase 1 |
| 极坐标方向图 | `emstudio-components` | Phase 1 |
| Marker / Cursor 交互 | `emstudio-components` | Phase 1 |
| CSV 导出 | `emstudio-components` | Phase 1 |
| 数据表格 | `emstudio-components` | Phase 1 |

### Phase 3：3D 网格与场渲染

| 任务 | 模块 | 依赖 |
|------|------|------|
| MSH 加载器 | `emstudio-domain` | — |
| 网格表面提取 | `emstudio-render` | MSH 加载器 |
| 网格线框/实体渲染 | `emstudio-render` | 表面提取 |
| OrbitCamera + 交互 | `emstudio-render` | — |
| 色标系统（1D Texture） | `emstudio-render` | — |
| 表面云图管线 | `emstudio-render` | 色标系统 |
| .emsfld mmap 加载 | `emstudio-domain` | — |
| 场叠加表面渲染 | `emstudio-render` | 云图管线 + mmap |

### Phase 4：3D 高级可视化

| 任务 | 模块 | 依赖 |
|------|------|------|
| 矢量箭头渲染（instanced） | `emstudio-render` | Phase 3 |
| 切面可视化 | `emstudio-render` | Phase 3 |
| 等值面提取（Marching Tet） | `emstudio-render` | Phase 3 |
| 3D 远场方向图 | `emstudio-render` | Phase 3 |
| 相位动画 | `emstudio-render` | Phase 3 |
| GPU Picking | `emstudio-render` | Phase 3 |
| 场值探针 | `emstudio-render` | Picking + Phase 3 |
| PNG 截图导出 | `emstudio-render` | Phase 3 |

### Phase 5：Q3D 准静态可视化

| 任务 | 模块 | 依赖 |
|------|------|------|
| RlcgMatrixData JSON 加载 | `emstudio-domain` | — |
| RLCG 矩阵元素 vs 频率曲线图 | `emstudio-components` | Phase 1 + RlcgMatrixData |
| RLCG 矩阵数据表格（热力图着色） | `emstudio-components` | Phase 1 + RlcgMatrixData |
| Q3D 收敛曲线（Delta Energy） | `emstudio-components` | Phase 1 |
| Q3D 电流密度场叠加（VolumeCurrent） | `emstudio-render` | Phase 3 |
| Q3D 电荷分布场叠加（ChargeDistribution） | `emstudio-render` | Phase 3 |
| Q3D 欧姆损耗场叠加（OhmicLoss） | `emstudio-render` | Phase 3 |
| RLCG → S 参数转换 + Touchstone 导出 | `emstudio-domain` | touchstone crate |
| 等效电路 SPICE 导出 | `emstudio-domain` | RlcgMatrixData |
| RLCG 矩阵 CSV 导出 | `emstudio-components` | RlcgMatrixData |

---

## 12. 参考资料

- [Ansys HFSS Creating Reports](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/ReportsandPostProc/CreatingReports.htm) — HFSS 报告系统
- [Ansys HFSS Field Overlays](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/ReportsandPostProc/FieldOverlays.htm) — HFSS 场叠加
- [Ansys HFSS Far Field Plots](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/ReportsandPostProc/FarFieldPlots.htm) — HFSS 远场显示
- [Ansys HFSS Smith Chart](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/DesignerHelp/Content/ReportsandPostProc/SmithChart.htm) — Smith 圆图
- [egui_plot Documentation](https://docs.rs/egui_plot/latest/egui_plot/) — egui 2D 绑图库
- [wgpu Documentation](https://docs.rs/wgpu/latest/wgpu/) — Rust WebGPU 实现
- [WGSL Specification](https://www.w3.org/TR/WGSL/) — WebGPU 着色语言
- [Gmsh MSH 4.1 Format](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format) — 网格文件格式
- [Marching Tetrahedra](https://en.wikipedia.org/wiki/Marching_tetrahedra) — 等值面提取算法
- [Smith Chart Theory](https://en.wikipedia.org/wiki/Smith_chart) — Smith 圆图原理
- [Q3D Creating Reports](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/ReportsandPostProc/CreatingReportsQ3D.htm) — Q3D 报告系统
- [Q3D Viewing Matrix Data](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/ViewingMatrixDatainQ3D.htm) — Q3D 矩阵数据查看
- [Q3D Plotting Field Overlays](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Q3DExtractor/Content/ReportsandPostProc/PlottingFieldOverlaysinQ3D.htm) — Q3D 场叠加显示
- [Q3D Exporting Equivalent Circuit](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/ExportingQ3DExtractorEquivalentCircuitData.htm) — Q3D 等效电路导出
- [Q3D Exporting S-Parameters](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Q3DExtractor/Content/Q3D/ExportingSParameterData.htm) — Q3D S 参数导出
- [PyAEDT Q3D Post-Processing](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.q3d.Q3d.html) — PyAEDT Q3D 后处理 API
