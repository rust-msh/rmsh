# EMStudio 电磁仿真工程文件方案设计

## 1. 概述

本文档为 EMStudio 设计一套完整的电磁仿真工程文件方案。方案参考 Ansys HFSS/AEDT 的工程文件体系（[Ansys Electronics Desktop Files](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/Variables/ANSYSElectronicsDesktopFiles.htm)、[PyAEDT Project Configuration](https://aedt.docs.pyansys.com/version/stable/User_guide/pyaedt_file_data/project.html)），以及 Ansys Q3D Extractor 的寄生参数提取体系（[Q3D Extractor Help](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/GettingStarted/Q3DExtractorGettingStartedGuides.htm)、[PyAEDT Q3D API](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.q3d.Q3d.html)），结合 EMStudio 自身技术栈（Rust + egui + wgpu + serde）进行适配设计。

> **EMStudio 求解器覆盖范围**：EMStudio 同时支持 **全波求解**（对应 HFSS 的驱动模态/端子/本征模求解）和 **准静态寄生参数提取**（对应 Q3D Extractor 的 RLCG 矩阵提取）。两类求解共享同一工程文件格式、几何建模和材料系统，但在边界条件、激励定义、分析设置和结果类型上有所区别。

### 1.1 设计目标

- **可序列化**：基于 JSON，利用 serde 自然映射到 Rust 类型
- **自描述**：文件包含完整的版本与元数据信息，可独立解析
- **层次清晰**：参考 HFSS 的 Project → Design → Model/Boundary/Setup 层次结构
- **可扩展**：预留扩展点，便于后续增加新的求解器类型、边界条件、后处理功能
- **跨平台**：同一格式在 Native 和 WASM 环境下通用

### 1.2 文件扩展名

| 文件/目录 | 扩展名 | 说明 |
|-----------|--------|------|
| 工程文件 | `.emsp` | EMStudio Project，主工程文件（JSON 格式） |
| 结果目录 | `.emsp.results/` | 仿真结果数据目录 |
| 锁文件 | `.emsp.lock` | 工程打开锁，防止并发写入 |
| 自动保存 | `.emsp.auto` | 自动恢复用备份 |
| 材料库 | `.emsm` | EMStudio Materials，可复用的材料库文件 |

> **参考**：HFSS 使用 `.aedt` 主文件 + `.aedtresults/` 结果目录 + `.aedt.lock` 锁文件的组织方式。

---

## 2. 工程层次结构

参照 [HFSS 工程树](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/PDFs/An%20Introduction%20to%20HFSS.pdf) 和 [Q3D Extractor 工程树](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/GettingStarted/Q3DExtractorGettingStartedGuides.htm) 的层级组织：

```
EMStudio Project (.emsp)
 ├── metadata                    # 文件元数据
 ├── variables                   # 工程级参数化变量（$前缀，跨设计共享）
 ├── datasets                    # 工程级数据集（频率曲线、查找表等）
 └── designs[]                   # 设计列表（一个工程可含多个设计）
      ├── general                # 设计基本信息与求解类型
      ├── design_settings        # 设计全局设置（端口归一化、验证选项等）
      ├── local_variables        # 设计级变量（仅本设计可见）
      ├── definitions            # 集中定义层（定义-引用架构）
      │    ├── materials[]       #   材料定义
      │    ├── coordinate_systems[]  #   坐标系定义
      │    └── named_selections[]    #   命名选择（面/边/顶点集）
      ├── geometry               # 几何模型（历史记录式）
      │    ├── operations[]      #   建模操作历史（引用 definitions 中的材料/坐标系）
      │    └── objects[]         #   由 operations 回放生成的最终对象快照
      ├── boundaries[]           # 边界条件（引用 geometry.objects / named_selections）
      ├── excitations[]          # 激励/端口（HFSS: 波端口/集总端口; Q3D: 端子/源-汇）
      ├── nets[]                 # 网络定义（Q3D 专用：导体分组与源汇分配）
      ├── mesh_operations[]      # 网格控制（引用 geometry.objects / named_selections）
      ├── analysis_setups[]      # 分析设置（可引用变量作为频率值等）
      │    └── frequency_sweeps[]
      ├── radiation              # 辐射设置（HFSS 专用：远场球面、近场采样定义）
      │    ├── far_field_setups[]
      │    └── near_field_setups[]
      ├── output_variables[]     # 输出变量（从仿真结果派生的表达式）
      ├── field_overlays[]       # 场叠加显示定义（E/H/J 可视化方案）
      ├── optimetrics[]          # 参数化扫描、优化、灵敏度分析
      └── reports[]              # 后处理报告（图表、数据表、矩阵表）
```

> **HFSS vs Q3D 专用节点**：`radiation`（远场/近场）仅适用于 HFSS 全波设计；`nets`（网络/源汇分配）仅适用于 Q3D 准静态设计。其余节点（geometry、materials、boundaries、excitations、analysis_setups 等）两者共享，但各自有不同的子类型。

---

## 3. 定义-引用架构（Definition-Reference Architecture）

参考 [HFSS 变量系统](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/Variables/DesignVariables.htm) 和 [材料表达式定义](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/HFSS3DLayout/Content/Materials/DefiningMaterialPropertiesasExpressions.htm)：HFSS 采用**"集中定义，名称引用"**的模式——变量、材料、坐标系、数据集等在专门的节点下统一定义，然后被几何操作、边界条件、激励、分析设置等通过**名称字符串**引用。修改定义时，所有引用处自动跟随更新。

EMStudio 采用相同的架构，形成清晰的**依赖关系图**。

### 3.1 引用关系全景图

```
┌─────────────────────────────────────────────────────────────────────┐
│                        EMStudio Project                             │
│                                                                     │
│  ┌─────────────────┐    ┌──────────────────┐                       │
│  │  Project Vars   │    │    Datasets      │                       │
│  │  ($freq, $sub_h)│    │  (ds_perm_vs_f)  │                       │
│  └────┬───┬────────┘    └────────┬─────────┘                       │
│       │   │                      │                                  │
│  ─ ─ ─│─ ─│─ ─ ─ ─ ─ ─ ─ ─ ─ ─ │─ ─ ─ ─ ─ ─ ─ Design ─ ─ ─ ─   │
│       │   │                      │                                  │
│  ┌────▼───▼────────┐   ┌────────▼─────────┐  ┌─────────────────┐  │
│  │  Local Variables │   │   Materials      │  │ CoordinateSys   │  │
│  │  (patch_l, ...)  │   │ (copper, FR4...) │  │ (Global, Local) │  │
│  └──┬──┬──┬─────────┘   └──┬───┬───┬───────┘  └──┬──────────────┘  │
│     │  │  │                │   │   │              │                  │
│     │  │  │   ┌────────────┘   │   │              │                  │
│     │  │  │   │                │   │              │                  │
│  ┌──▼──▼──▼───▼────┐   ┌──────▼───▼──────────────▼──┐              │
│  │   Geometry      │   │  Named Selections          │              │
│  │  (operations)   │◄──│  (faces, edges, groups)     │              │
│  └──┬──────────────┘   └──┬──────────────────────────┘              │
│     │                     │                                          │
│     │  ┌──────────────────┤                                          │
│     │  │                  │                                          │
│  ┌──▼──▼───┐  ┌──────────▼──┐  ┌────────────────┐                  │
│  │Boundaries│  │ Excitations │  │ Mesh Operations│                  │
│  └──────────┘  └─────────────┘  └────────────────┘                  │
│         │              │                │                            │
│         └──────────────┼────────────────┘                            │
│                        ▼                                             │
│               ┌─────────────────┐                                   │
│               │ Analysis Setups │ ◄── 引用 $freq 等变量             │
│               └────────┬────────┘                                   │
│                        ▼                                             │
│               ┌─────────────────┐                                   │
│               │   Optimetrics   │ ◄── 驱动变量遍历                   │
│               └─────────────────┘                                   │
└─────────────────────────────────────────────────────────────────────┘
```

### 3.2 变量系统

参考 [HFSS 工程变量与设计变量](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/Variables/ProjectVariables.htm)：

| 变量类型 | 作用域 | 命名规则 | HFSS 对应 | 用途 |
|---------|--------|---------|----------|------|
| 工程变量 | 全工程（所有设计共享） | `$` 前缀 | Project Variable | 材料属性、跨设计共享尺寸 |
| 设计变量 | 单个设计内 | 无前缀 | Design Variable | 几何尺寸、本设计的参数 |
| 内置变量 | 求解器提供 | 保留名 | Intrinsic Variable | `Freq`、`Phase`、`Time` 等 |

#### 变量定义 JSON

```json
{
  "variables": {
    "$freq": {
      "value": "2.4GHz",
      "description": "工作频率",
      "unit_type": "Frequency"
    },
    "$eps_sub": {
      "value": "4.4",
      "description": "基板介电常数",
      "unit_type": "None"
    },
    "$sub_h": {
      "value": "1.6mm",
      "description": "基板厚度",
      "unit_type": "Length"
    }
  }
}
```

#### 设计级局部变量

```json
{
  "local_variables": {
    "patch_l": {
      "value": "28.5mm",
      "expression": null,
      "description": "贴片长度"
    },
    "patch_w": {
      "value": "37.0mm",
      "expression": null,
      "description": "贴片宽度"
    },
    "patch_x": {
      "value": null,
      "expression": "(60mm - patch_l) / 2",
      "description": "贴片X位置（居中表达式）"
    }
  }
}
```

> **表达式（Expression）**：变量的值可以是常量，也可以是引用其他变量的数学表达式。引擎在回放/求解时递归解析表达式，形成依赖链。

### 3.3 数据集（Datasets）

参考 [HFSS Dataset 表达式](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS3DLayout/Content/Variables/UsingDatasetExpressions.htm)：数据集是命名的**查找表**，可在材料属性、边界条件等处通过 `$ds_name(Freq)` 表达式引用。

```json
{
  "datasets": {
    "$ds_perm_vs_freq": {
      "description": "FR4 介电常数随频率变化",
      "independent_variable": "Freq",
      "independent_unit": "GHz",
      "dependent_unit": "None",
      "data": [
        { "x": 1.0,  "y": 4.50 },
        { "x": 2.0,  "y": 4.45 },
        { "x": 5.0,  "y": 4.35 },
        { "x": 10.0, "y": 4.20 },
        { "x": 20.0, "y": 4.05 }
      ],
      "interpolation": "PiecewiseLinear"
    },
    "$ds_cond_vs_temp": {
      "description": "铜电导率随温度变化",
      "independent_variable": "Temp",
      "independent_unit": "cel",
      "dependent_unit": "S/m",
      "data": [
        { "x": 20.0,  "y": 58000000.0 },
        { "x": 50.0,  "y": 54000000.0 },
        { "x": 100.0, "y": 48000000.0 }
      ],
      "interpolation": "PiecewiseLinear"
    }
  }
}
```

### 3.4 命名选择（Named Selections）

HFSS 中命名选择用于将一组面、边或顶点绑定到一个逻辑名称上，供边界、激励、网格操作引用。这样即使参数化变化导致拓扑变化，引用关系也能通过命名选择保持稳定。

```json
{
  "named_selections": [
    {
      "name": "GND_Bottom",
      "type": "Face",
      "selection": [
        { "object": "GND_Plane", "face": "ZMin" }
      ],
      "description": "地平面底面"
    },
    {
      "name": "RadiationFaces",
      "type": "Face",
      "selection": [
        { "object": "AirBox", "face": "XMin" },
        { "object": "AirBox", "face": "XMax" },
        { "object": "AirBox", "face": "YMin" },
        { "object": "AirBox", "face": "YMax" },
        { "object": "AirBox", "face": "ZMax" }
      ],
      "description": "辐射边界面集合"
    },
    {
      "name": "FeedPort",
      "type": "Face",
      "selection": [
        { "object": "FeedLine", "face": "YMin" }
      ],
      "description": "馈电端口面"
    },
    {
      "name": "AllConductors",
      "type": "Object",
      "selection": ["GND_Plane", "Patch", "FeedLine"],
      "description": "所有导体对象"
    }
  ]
}
```

### 3.5 材料定义中的引用

材料属性可以是**常量**、**变量引用**或**数据集引用**：

```json
{
  "name": "FR4_parametric",
  "category": "Dielectric",
  "properties": {
    "permittivity": {
      "type": "expression",
      "expression": "$eps_sub"
    },
    "permeability": { "type": "constant", "value": 1.0 },
    "conductivity": { "type": "constant", "value": 0.0 },
    "dielectric_loss_tangent": {
      "type": "dataset",
      "dataset": "$ds_losstangent_vs_freq",
      "independent_variable": "Freq"
    },
    "mass_density": { "type": "constant", "value": 1900.0 }
  }
}
```

> **三种属性值类型**：
> - `constant`：固定数值
> - `expression`：引用变量或数学表达式（如 `"$eps_sub * 1.02"`）
> - `dataset`：引用数据集，随独立变量（如 `Freq`）插值

### 3.6 引用在各子系统中的体现

| 被引用的定义 | 引用方 | 引用方式 | 示例 |
|-------------|--------|---------|------|
| 工程变量 `$freq` | Analysis Setup | `solution_frequency` 字段值 | `"solution_frequency": "$freq"` |
| 工程变量 `$eps_sub` | Material 属性 | expression 类型 | `"permittivity": {"type":"expression","expression":"$eps_sub"}` |
| 设计变量 `patch_l` | Geometry Operation | parameters 中的尺寸值 | `"size": ["patch_l", "patch_w", 0.035]` |
| 数据集 `$ds_perm_vs_freq` | Material 属性 | dataset 类型 | `"permittivity": {"type":"dataset","dataset":"$ds_perm_vs_freq"}` |
| 材料名 `copper` | Geometry Operation | attributes.material | `"material": "copper"` |
| 坐标系 `Local_CS1` | Geometry Operation | attributes.coordinate_system | `"coordinate_system": "Local_CS1"` |
| 命名选择 `GND_Bottom` | Boundary | assignment.targets | `"targets": ["@GND_Bottom"]` |
| 命名选择 `FeedPort` | Excitation | assignment.targets | `"targets": ["@FeedPort"]` |
| 命名选择 `AllConductors` | Mesh Operation | assignment.targets | `"targets": ["@AllConductors"]` |
| 几何对象 `AirBox` | Boundary | assignment.targets | `"targets": ["AirBox"]` |

> **`@` 前缀约定**：引用命名选择时使用 `@` 前缀（如 `@GND_Bottom`），与直接引用几何对象名（如 `AirBox`）区分。

### 3.7 依赖关系验证

修改或删除某个定义时，引擎需检查引用关系：

```rust
/// 定义-引用的依赖关系图
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DependencyGraph {
    /// 节点：所有可被引用的定义
    pub definitions: HashMap<DefinitionId, DefinitionKind>,
    /// 边：引用方 → 被引用的定义
    pub references: Vec<Reference>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DefinitionKind {
    ProjectVariable,
    DesignVariable,
    Dataset,
    Material,
    CoordinateSystem,
    NamedSelection,
    GeometryObject,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Reference {
    pub from: ReferenceSource,   // 引用方（如 "Boundary:Radiation1"）
    pub to: DefinitionId,         // 被引用的定义（如 "Material:copper"）
    pub field: String,            // 引用字段（如 "assignment.targets"）
}

impl DependencyGraph {
    /// 查找某定义的所有引用者（用于删除前检查）
    pub fn find_dependents(&self, def_id: &DefinitionId) -> Vec<&Reference> {
        self.references.iter().filter(|r| &r.to == def_id).collect()
    }

    /// 验证所有引用是否有效（定义存在且类型匹配）
    pub fn validate(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        for reference in &self.references {
            if !self.definitions.contains_key(&reference.to) {
                errors.push(ValidationError::DanglingReference {
                    from: reference.from.clone(),
                    missing: reference.to.clone(),
                    field: reference.field.clone(),
                });
            }
        }
        errors
    }
}
```

**验证规则：**

| 场景 | 行为 |
|------|------|
| 删除材料 `copper` | 检查是否有几何对象引用 → 若有则警告用户 |
| 删除变量 `$freq` | 检查材料属性 / Setup 是否引用 → 阻止删除或级联清理 |
| 删除几何对象 `AirBox` | 检查边界 / 激励 / 网格操作是否引用 → 级联删除或警告 |
| 删除命名选择 | 同上，检查边界 / 激励 / 网格操作引用 |
| 修改变量值 | 标记依赖链上所有结果为"过期"，需重新回放/求解 |
| 重命名对象 | 自动更新所有引用中的名称字符串 |

### 3.8 Rust 类型定义

```rust
/// 属性值 — 支持常量、表达式、数据集三种来源
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum PropertyValue {
    #[serde(rename = "constant")]
    Constant { value: f64 },
    #[serde(rename = "expression")]
    Expression { expression: String },
    #[serde(rename = "dataset")]
    Dataset {
        dataset: String,
        independent_variable: String,
    },
}

/// 变量定义
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Variable {
    pub value: Option<String>,
    pub expression: Option<String>,
    pub description: String,
    pub unit_type: Option<String>,
}

/// 数据集定义
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetDefinition {
    pub description: String,
    pub independent_variable: String,
    pub independent_unit: String,
    pub dependent_unit: String,
    pub data: Vec<DataPoint>,
    pub interpolation: InterpolationType,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DataPoint {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum InterpolationType {
    PiecewiseLinear,
    CubicSpline,
    Debye,
    DjordjevicSarkar,
}

/// 命名选择
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NamedSelection {
    pub name: String,
    pub selection_type: SelectionType,
    pub selection: Vec<serde_json::Value>,
    pub description: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SelectionType {
    Face,
    Edge,
    Vertex,
    Object,
}

/// 设计中的集中定义层
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Definitions {
    pub materials: Vec<Material>,
    pub coordinate_systems: Vec<CoordinateSystem>,
    pub named_selections: Vec<NamedSelection>,
}
```

---

## 4. 完整 JSON Schema 定义

### 4.1 顶层结构

```json
{
  "metadata": {
    "version": "1.0.0",
    "application": "EMStudio",
    "created_at": "2026-04-03T10:00:00Z",
    "modified_at": "2026-04-03T12:30:00Z",
    "author": "user@example.com",
    "description": "微带贴片天线仿真工程"
  },
  "variables": {
    "$freq": { "value": "2.4GHz", "description": "工作频率", "unit_type": "Frequency" },
    "$eps_sub": { "value": "4.4", "description": "基板介电常数", "unit_type": "None" },
    "$sub_h": { "value": "1.6mm", "description": "基板厚度", "unit_type": "Length" }
  },
  "datasets": {
    "$ds_losstangent_vs_freq": {
      "description": "FR4 损耗正切随频率变化",
      "independent_variable": "Freq",
      "independent_unit": "GHz",
      "dependent_unit": "None",
      "data": [
        { "x": 1.0, "y": 0.018 },
        { "x": 5.0, "y": 0.020 },
        { "x": 10.0, "y": 0.025 }
      ],
      "interpolation": "PiecewiseLinear"
    }
  },
  "designs": [ "..." ]
}
```

### 4.2 Design（设计）

每个 Design 对应一个独立的电磁仿真问题，参考 HFSS 中一个 Project 可包含多个 Design 的概念。

```json
{
  "id": "design-001",
  "name": "Patch Antenna",
  "solution_type": "DrivenModal",
  "units": "mm",
  "design_settings": { "..." },
  "local_variables": {
    "patch_l": { "value": "28.5mm", "description": "贴片长度" },
    "patch_w": { "value": "37.0mm", "description": "贴片宽度" },
    "patch_x": { "expression": "(60mm - patch_l) / 2", "description": "贴片X位置（居中）" }
  },
  "definitions": {
    "materials": [ "..." ],
    "coordinate_systems": [ "..." ],
    "named_selections": [ "..." ]
  },
  "geometry": { "..." },
  "boundaries": [ "..." ],
  "excitations": [ "..." ],
  "mesh_operations": [ "..." ],
  "analysis_setups": [ "..." ],
  "radiation": { "..." },
  "output_variables": [ "..." ],
  "field_overlays": [ "..." ],
  "optimetrics": [ "..." ],
  "reports": [ "..." ]
}
```

> **注意**：`definitions` 中的材料、坐标系、命名选择等是**定义层**，被 geometry/boundaries/excitations 等通过名称引用。详见 §3 定义-引用架构。

**`solution_type`** 枚举值（参考 HFSS 求解类型和 [Q3D 求解类型](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/GettingStarted/Q3DExtractorGettingStartedGuides.htm)）：

| 值 | 类别 | 说明 |
|----|------|------|
| `DrivenModal` | HFSS 全波 | 驱动模态求解（S 参数，模态阻抗） |
| `DrivenTerminal` | HFSS 全波 | 驱动端子求解（S 参数，端子电压/电流） |
| `Eigenmode` | HFSS 全波 | 本征模求解（谐振频率，Q 值） |
| `Transient` | HFSS 全波 | 瞬态求解（时域） |
| `SBRPlus` | HFSS 全波 | 射线追踪求解（电大尺寸问题） |
| `Q3D_DCRL` | Q3D 准静态 | 直流电阻和低频电感提取 |
| `Q3D_ACRL` | Q3D 准静态 | 交流电阻和电感提取（含趋肤效应/邻近效应） |
| `Q3D_C` | Q3D 准静态 | 电容矩阵提取（拉普拉斯方程静电求解） |
| `Q3D_CG` | Q3D 准静态 | 电容-电导矩阵提取（含介质损耗） |

> **Q3D 求解器特性**：Q3D 的 3D 寄生提取采用 **矩量法 (MoM)** 配合 **快速多极子加速 (FMM)**，是基于面的求解器而非体积网格。DC 电阻提取和 2D Extractor 模式使用 FEM 四面体网格。这与 HFSS 的全 FEM 体积求解形成互补。

### 4.3 Geometry（几何模型 — 历史记录式建模）

参考 [HFSS History Tree](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Subsystems/An%20Introduction%20to%20HFSS/Content/HistoryTree.htm)：HFSS 的 3D Modeler 采用**历史记录式建模**（History-based Modeling），模型树中每个对象下记录了完整的操作链——创建基本体、布尔运算、倒角、扫掠等。参数修改时，引擎会**按顺序回放操作历史**来重建几何体，这是参数化设计和 Optimetrics 的基础。

EMStudio 的几何模块同样采用 **操作历史（Operation History）** 作为核心数据结构，而非仅存储最终几何体。

#### 4.3.1 设计思想

```
用户操作 → 追加 Operation → 回放全部 Operations → 生成 Objects 快照 → 渲染
                                    ↑ 参数变更时从这里重新回放
```

| 概念 | HFSS 对应 | EMStudio 对应 | 说明 |
|------|----------|--------------|------|
| 操作历史 | History Tree | `geometry.operations[]` | 有序操作列表，记录每一步建模动作 |
| 最终对象 | Solid/Sheet/Line | `geometry.objects[]` | 回放后生成的当前几何体快照 |
| 操作属性编辑 | Properties 面板 | 修改 operation 参数后重新回放 | 参数化驱动 |
| 清除历史 | Purge History | 丢弃 operations，仅保留 objects | 不可逆，用于精简文件 |

#### 4.3.2 Geometry JSON 结构

```json
{
  "coordinate_systems": [
    {
      "name": "Global",
      "type": "Cartesian",
      "origin": [0.0, 0.0, 0.0],
      "x_axis": [1.0, 0.0, 0.0],
      "y_axis": [0.0, 1.0, 0.0]
    }
  ],
  "operations": [
    {
      "step": 1,
      "command": "CreateBox",
      "result_object": "Substrate",
      "parameters": {
        "position": [0.0, 0.0, 0.0],
        "size": ["$sub_w", "$sub_w", "$sub_h"]
      },
      "attributes": {
        "material": "FR4_epoxy",
        "solve_inside": true,
        "color": [128, 200, 128],
        "transparency": 0.4,
        "coordinate_system": "Global",
        "group": "Antenna"
      }
    },
    {
      "step": 2,
      "command": "CreateBox",
      "result_object": "GND_Plane",
      "parameters": {
        "position": [0.0, 0.0, 0.0],
        "size": ["$sub_w", "$sub_w", 0.035]
      },
      "attributes": {
        "material": "copper",
        "solve_inside": false,
        "color": [255, 180, 50],
        "transparency": 0.0,
        "group": "Antenna"
      }
    },
    {
      "step": 3,
      "command": "CreateBox",
      "result_object": "Patch",
      "parameters": {
        "position": ["$patch_x", "$patch_y", "$sub_h"],
        "size": ["$patch_l", "$patch_w", 0.035]
      },
      "attributes": {
        "material": "copper",
        "solve_inside": false,
        "color": [255, 180, 50],
        "transparency": 0.0,
        "group": "Antenna"
      }
    },
    {
      "step": 4,
      "command": "CreateCylinder",
      "result_object": "Via1",
      "parameters": {
        "center": [30.0, 30.0, 0.0],
        "radius": 0.5,
        "height": "$sub_h",
        "axis": "Z"
      },
      "attributes": {
        "material": "copper",
        "solve_inside": false,
        "color": [255, 180, 50],
        "group": "Feed"
      }
    },
    {
      "step": 5,
      "command": "Subtract",
      "result_object": "Substrate",
      "parameters": {
        "blank": "Substrate",
        "tool": ["Via1"],
        "keep_tool": false
      }
    },
    {
      "step": 6,
      "command": "CreatePolyline",
      "result_object": "FeedProfile",
      "parameters": {
        "points": [[28.5, 0, 1.635], [28.5, 11.5, 1.635]],
        "closed": false
      },
      "attributes": {
        "group": "Feed"
      }
    },
    {
      "step": 7,
      "command": "SweepAlongPath",
      "result_object": "FeedLine",
      "parameters": {
        "profile": "FeedProfile",
        "path_points": [[0, 0, 0], [3.0, 0, 0]],
        "draft_angle": "0deg",
        "twist_angle": "0deg"
      },
      "attributes": {
        "material": "copper",
        "solve_inside": false,
        "group": "Feed"
      }
    },
    {
      "step": 8,
      "command": "DuplicateAlongLine",
      "result_object": "Via1_Array",
      "parameters": {
        "source": "Via1",
        "direction": [10.0, 0.0, 0.0],
        "count": 5
      }
    },
    {
      "step": 9,
      "command": "SetMaterial",
      "parameters": {
        "target": "Substrate",
        "material": "Rogers_RO4003C"
      }
    },
    {
      "step": 10,
      "command": "CreateBox",
      "result_object": "AirBox",
      "parameters": {
        "position": [-30.0, -30.0, -30.0],
        "size": [120.0, 120.0, 62.0]
      },
      "attributes": {
        "material": "vacuum",
        "solve_inside": true,
        "color": [200, 200, 255],
        "transparency": 0.95,
        "group": "Environment"
      }
    }
  ],
  "objects": [
    {
      "id": 1,
      "name": "Substrate",
      "derived_from_step": 5,
      "material": "Rogers_RO4003C",
      "solve_inside": true,
      "color": [128, 200, 128],
      "transparency": 0.4,
      "group": "Antenna",
      "bounding_box": { "min": [0,0,0], "max": [60,60,1.6] }
    },
    {
      "id": 2,
      "name": "GND_Plane",
      "derived_from_step": 2,
      "material": "copper",
      "solve_inside": false,
      "color": [255, 180, 50],
      "transparency": 0.0,
      "group": "Antenna",
      "bounding_box": { "min": [0,0,0], "max": [60,60,0.035] }
    },
    {
      "id": 3,
      "name": "Patch",
      "derived_from_step": 3,
      "material": "copper",
      "solve_inside": false,
      "color": [255, 180, 50],
      "transparency": 0.0,
      "group": "Antenna",
      "bounding_box": { "min": [15.75,11.5,1.635], "max": [44.25,48.5,1.67] }
    }
  ]
}
```

> **`objects`** 是回放 `operations` 后的**生成快照**，用于快速加载和渲染。当 `operations` 存在时，`objects` 可由引擎重建；当历史被清除（Purge History）时，`objects` 成为唯一的几何数据源。

#### 4.3.3 操作命令（Operation Commands）总表

参考 [HFSS History Tree 操作](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/GettingStarted/WorkingwiththeHistoryTree.htm) 和 [HFSS 布尔运算](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/HFSS3DLayout/Content/Modeler/BooleanOperationsonObjects.htm)：

**基本体创建命令：**

| command | parameters | 说明 |
|---------|-----------|------|
| `CreateBox` | `position`, `size` | 长方体 |
| `CreateCylinder` | `center`, `radius`, `height`, `axis` | 圆柱体 |
| `CreateSphere` | `center`, `radius` | 球体 |
| `CreateCone` | `center`, `radius_bottom`, `radius_top`, `height` | 圆锥/截锥 |
| `CreateTorus` | `center`, `major_radius`, `minor_radius` | 环面 |
| `CreatePolyline` | `points`, `closed` | 折线/多边形线 |
| `CreateRectangle` | `center`, `width`, `height`, `axis` | 矩形面 |
| `CreateCircle` | `center`, `radius`, `axis` | 圆形面 |

**布尔运算命令：**

| command | parameters | 说明 |
|---------|-----------|------|
| `Unite` | `objects[]` | 合并多个对象 |
| `Subtract` | `blank`, `tool[]`, `keep_tool` | 从 blank 中减去 tool |
| `Intersect` | `objects[]` | 保留重叠部分 |

**变换/修改命令：**

| command | parameters | 说明 |
|---------|-----------|------|
| `Move` | `target`, `vector` | 平移 |
| `Rotate` | `target`, `axis`, `angle` | 旋转 |
| `Mirror` | `target`, `plane` | 镜像 |
| `Scale` | `target`, `factor` | 缩放 |
| `DuplicateAlongLine` | `source`, `direction`, `count` | 沿线阵列复制 |
| `DuplicateAroundAxis` | `source`, `axis`, `angle`, `count` | 环形阵列复制 |

**扫掠命令（参考 [HFSS Sweep](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/Modeler/SweepingAroundanAxis.htm)）：**

| command | parameters | 说明 |
|---------|-----------|------|
| `SweepAlongVector` | `profile`, `vector`, `draft_angle` | 沿向量拉伸 |
| `SweepAlongPath` | `profile`, `path_points`, `twist_angle` | 沿路径扫掠 |
| `SweepAroundAxis` | `profile`, `axis`, `angle` | 绕轴旋转扫掠 |

**属性修改命令：**

| command | parameters | 说明 |
|---------|-----------|------|
| `SetMaterial` | `target`, `material` | 修改材料 |
| `SetColor` | `target`, `color`, `transparency` | 修改显示属性 |
| `Rename` | `target`, `new_name` | 重命名对象 |
| `SetGroup` | `target`, `group` | 设置分组 |
| `SetSolveInside` | `target`, `solve_inside` | 设置是否内部求解 |

**高级命令：**

| command | parameters | 说明 |
|---------|-----------|------|
| `Fillet` | `target`, `edges[]`, `radius` | 倒圆角 |
| `Chamfer` | `target`, `edges[]`, `distance` | 倒斜角 |
| `Section` | `target`, `plane` | 截面（3D→2D） |
| `Import` | `file_path`, `format`, `result_object` | 导入外部几何（STEP/STL/OBJ） |

#### 4.3.4 操作的 Rust 类型定义

```rust
/// 几何建模操作 — 历史记录中的一步
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometryOperation {
    pub step: u32,
    pub command: OperationCommand,
    /// 操作产出的对象名（创建/布尔运算结果），None 表示修改类操作
    pub result_object: Option<String>,
    /// 命令参数，不同 command 有不同 schema
    pub parameters: serde_json::Value,
    /// 可选的对象属性（材料、颜色等），仅创建类命令使用
    pub attributes: Option<ObjectAttributes>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OperationCommand {
    // 基本体创建
    CreateBox,
    CreateCylinder,
    CreateSphere,
    CreateCone,
    CreateTorus,
    CreatePolyline,
    CreateRectangle,
    CreateCircle,
    // 布尔运算
    Unite,
    Subtract,
    Intersect,
    // 变换
    Move,
    Rotate,
    Mirror,
    Scale,
    DuplicateAlongLine,
    DuplicateAroundAxis,
    // 扫掠
    SweepAlongVector,
    SweepAlongPath,
    SweepAroundAxis,
    // 属性修改
    SetMaterial,
    SetColor,
    Rename,
    SetGroup,
    SetSolveInside,
    // 高级
    Fillet,
    Chamfer,
    Section,
    Import,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObjectAttributes {
    pub material: Option<String>,
    pub solve_inside: Option<bool>,
    pub color: Option<[u8; 3]>,
    pub transparency: Option<f32>,
    pub coordinate_system: Option<String>,
    pub group: Option<String>,
}

/// 回放 operations 后生成的对象快照
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometryObject {
    pub id: u64,
    pub name: String,
    /// 产生此对象的最后一步操作的 step 编号
    pub derived_from_step: u32,
    pub material: String,
    pub solve_inside: bool,
    pub color: [u8; 3],
    pub transparency: f32,
    pub group: Option<String>,
    pub bounding_box: Option<BoundingBox>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundingBox {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

/// 完整几何数据
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Geometry {
    pub coordinate_systems: Vec<CoordinateSystem>,
    /// 建模操作历史（有序，按 step 排列）
    pub operations: Vec<GeometryOperation>,
    /// 操作回放后的对象快照（可由 operations 重建）
    pub objects: Vec<GeometryObject>,
}
```

#### 4.3.5 操作回放引擎

```rust
impl Geometry {
    /// 回放所有 operations，重建 objects 快照
    pub fn rebuild(&mut self) -> Result<(), GeometryError> {
        let mut state = GeometryState::new();
        for op in &self.operations {
            state.apply(op)?;
        }
        self.objects = state.into_objects();
        Ok(())
    }

    /// 修改某一步的参数后，从该步开始增量回放
    pub fn modify_and_rebuild(
        &mut self,
        step: u32,
        new_params: serde_json::Value,
    ) -> Result<(), GeometryError> {
        if let Some(op) = self.operations.iter_mut().find(|o| o.step == step) {
            op.parameters = new_params;
        }
        self.rebuild()
    }

    /// 在指定位置插入新操作并回放
    pub fn insert_operation(
        &mut self,
        after_step: u32,
        operation: GeometryOperation,
    ) -> Result<(), GeometryError> {
        // 重编号后续步骤
        for op in self.operations.iter_mut() {
            if op.step > after_step {
                op.step += 1;
            }
        }
        self.operations.push(operation);
        self.operations.sort_by_key(|o| o.step);
        self.rebuild()
    }

    /// 删除指定步骤的操作并回放（参考 HFSS Delete Command）
    pub fn delete_operation(&mut self, step: u32) -> Result<(), GeometryError> {
        self.operations.retain(|o| o.step != step);
        self.rebuild()
    }

    /// 清除历史，仅保留 objects 快照（不可逆，参考 HFSS Purge History）
    pub fn purge_history(&mut self) {
        self.operations.clear();
    }
}
```

#### 4.3.6 UI 模型树展示

在 EMStudio 的左侧 Dock Panel（模型树）中，几何对象按照 HFSS 模型树的方式展示，每个对象可展开查看其操作历史：

```
▼ Geometry
  ▼ Substrate         [Rogers_RO4003C]
    ├── CreateBox      (step 1)  — 点击可编辑参数
    ├── Subtract       (step 5)  — Via1 布尔减
    └── SetMaterial    (step 9)  — FR4 → Rogers
  ▼ GND_Plane          [copper]
    └── CreateBox      (step 2)
  ▼ Patch              [copper]
    └── CreateBox      (step 3)
  ▼ Via1               [copper]
    └── CreateCylinder (step 4)
  ▼ FeedLine           [copper]
    ├── CreatePolyline (step 6)
    └── SweepAlongPath (step 7)
  ▼ AirBox             [vacuum]
    └── CreateBox      (step 10)
```

> 选中某个操作步骤时，右侧 Properties 面板显示该步的可编辑参数；修改后触发从该步开始的增量回放。

### 4.4 Materials（材料）

参考 HFSS 材料属性和 [PyAEDT 材料 API](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.modules.material_lib.Materials.add_material.html)：

```json
[
  {
    "name": "FR4_epoxy",
    "category": "Dielectric",
    "properties": {
      "permittivity": 4.4,
      "permeability": 1.0,
      "conductivity": 0.0,
      "dielectric_loss_tangent": 0.02,
      "magnetic_loss_tangent": 0.0,
      "mass_density": 1900.0
    },
    "appearance": {
      "color": [128, 200, 128],
      "transparency": 0.3
    }
  },
  {
    "name": "copper",
    "category": "Conductor",
    "properties": {
      "permittivity": 1.0,
      "permeability": 0.999991,
      "conductivity": 58000000.0,
      "dielectric_loss_tangent": 0.0,
      "magnetic_loss_tangent": 0.0,
      "mass_density": 8933.0
    },
    "appearance": {
      "color": [255, 180, 50],
      "transparency": 0.0
    }
  },
  {
    "name": "vacuum",
    "category": "Dielectric",
    "properties": {
      "permittivity": 1.0,
      "permeability": 1.0,
      "conductivity": 0.0,
      "dielectric_loss_tangent": 0.0,
      "magnetic_loss_tangent": 0.0,
      "mass_density": 0.0
    },
    "appearance": {
      "color": [200, 200, 255],
      "transparency": 0.95
    }
  }
]
```

**材料类别 `category`**：`Conductor` | `Dielectric` | `Magnetic` | `Composite` | `Gas`

**频率依赖材料**（预留扩展）：

```json
{
  "name": "Lossy_Substrate",
  "category": "Dielectric",
  "properties": {
    "permittivity": {
      "type": "frequency_dependent",
      "dataset": [
        { "frequency": "1GHz", "value": 4.5 },
        { "frequency": "5GHz", "value": 4.3 },
        { "frequency": "10GHz", "value": 4.1 }
      ]
    }
  }
}
```

### 4.5 Boundaries（边界条件）

参考 [HFSS 边界条件类型](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/HFSS/AssigningBoundariesinHFSSandHFSSIE.htm)：

```json
[
  {
    "name": "Radiation1",
    "type": "Radiation",
    "assignment": {
      "target_type": "Object",
      "targets": ["AirBox"]
    }
  },
  {
    "name": "GND",
    "type": "PerfectE",
    "assignment": {
      "target_type": "Face",
      "targets": [{"object": "Substrate", "face": "ZMin"}]
    }
  },
  {
    "name": "FiniteCond1",
    "type": "FiniteConductivity",
    "assignment": {
      "target_type": "Face",
      "targets": [{"object": "Patch", "face": "ZMax"}]
    },
    "properties": {
      "conductivity": 58000000.0,
      "roughness": "0um"
    }
  },
  {
    "name": "Sym1",
    "type": "Symmetry",
    "assignment": {
      "target_type": "Face",
      "targets": [{"object": "AirBox", "face": "YMin"}]
    },
    "properties": {
      "symmetry_type": "PerfectE"
    }
  }
]
```

**边界类型 `type`** 枚举：

| 类型 | 适用 | 说明 | 关键属性 |
|------|------|------|---------|
| `PerfectE` | HFSS / Q3D | 理想电导体（PEC） | — |
| `PerfectH` | HFSS / Q3D | 理想磁导体（PMC） | — |
| `Radiation` | HFSS | 辐射边界（吸收边界） | — |
| `PML` | HFSS | 完美匹配层 | `num_layers`, `min_frequency` |
| `Impedance` | HFSS | 阻抗边界 | `resistance`, `reactance` |
| `FiniteConductivity` | HFSS / Q3D | 有限电导率 | `conductivity`, `roughness` |
| `Symmetry` | HFSS / Q3D | 对称面 | `symmetry_type`: PerfectE / PerfectH |
| `MasterSlave` | HFSS | 主从周期边界 | `master`, `slave`, `phase_delay` |
| `ThinConductor` | Q3D | 薄导体边界（PCB 走线、薄涂层） | `conductivity`, `thickness` |
| `InfiniteGroundPlane` | Q3D | 无限大地平面 | — |
| `OpenBoundary` | Q3D | 开放边界（场延伸至无穷远） | — |

> **Q3D 边界条件说明**：`ThinConductor` 是 Q3D 特有的高效边界条件，用于建模 PCB 走线、薄膜涂层等不需要体积网格的薄导体表面。`InfiniteGroundPlane` 用于建模半无限大参考地平面。`OpenBoundary` 允许场向无穷远自然衰减。参考 [Q3D Thin Conductor Boundaries](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/AssigningThinConductorBoundaries.htm)。

### 4.6 Excitations（激励/端口）

参考 [HFSS 端口类型](https://innovationspace.ansys.com/courses/wp-content/uploads/sites/5/2021/07/HFSS_3DLGS_2019R3_EN_LE04_Ports-1.pdf) 和 [PyAEDT wave_port](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.hfss.Hfss.wave_port.html)：

```json
[
  {
    "name": "Port1",
    "type": "WavePort",
    "assignment": {
      "target_type": "Face",
      "targets": [{"object": "FeedLine", "face": "XMin"}]
    },
    "properties": {
      "num_modes": 1,
      "impedance": "50ohm",
      "deembed_distance": "0mm",
      "renormalize": true,
      "renorm_impedance": "50ohm"
    }
  },
  {
    "name": "Port2",
    "type": "LumpedPort",
    "assignment": {
      "target_type": "Face",
      "targets": [{"object": "FeedGap", "face": "ZMin"}]
    },
    "properties": {
      "impedance": "50ohm",
      "integration_line": {
        "start": [30.0, 30.0, 0.0],
        "end": [30.0, 30.0, 1.6]
      }
    }
  },
  {
    "name": "Incident1",
    "type": "IncidentWave",
    "properties": {
      "wave_type": "PlaneWave",
      "polarization": "Linear",
      "theta_inc": "0deg",
      "phi_inc": "0deg",
      "e_theta": "1V/m",
      "e_phi": "0V/m"
    }
  }
]
```

**激励类型 `type`** 枚举：

| 类型 | 适用 | 说明 | 适用场景 |
|------|------|------|---------|
| `WavePort` | HFSS | 波端口 | 波导、同轴线截面 |
| `LumpedPort` | HFSS | 集总端口 | PCB 走线、芯片封装 |
| `FloquetPort` | HFSS | Floquet 端口 | 周期结构（相控阵单元） |
| `IncidentWave` | HFSS | 入射波 | 散射/RCS 分析 |
| `VoltageDrop` | HFSS | 电压差激励 | 简单驱动 |
| `Source` | Q3D | 源端子（电流/电压注入点） | 导体的信号输入端 |
| `Sink` | Q3D | 汇端子（电流回路点） | 导体的信号回流端/接地端 |

> **Q3D 激励模型**：与 HFSS 的波端口（支持传播模式分解）不同，Q3D 使用 **集总端子模型**（源-汇对），在导体的边或面上分配源（Source）和汇（Sink），定义电流注入和回流路径。每个源-汇对构成一个提取通道，求解器计算所有通道之间的 RLCG 耦合矩阵。参考 [Q3D Terminal Assignments](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/AssigningTerminals.htm)。

### 4.6.1 Q3D 网络定义（Nets）

Q3D 的核心概念是 **网络（Net）**——将导体对象分组为命名的电气网络，并在其上分配源/汇端子。一个 Net 可包含多个几何对象（如走线、过孔、焊盘），构成一条完整的信号路径。

```json
{
  "nets": [
    {
      "name": "Signal1",
      "objects": ["Trace1", "Via1", "Pad1"],
      "terminals": [
        {
          "name": "T1",
          "type": "Source",
          "assignment": {
            "target_type": "Face",
            "targets": [{"object": "Trace1", "face": "XMin"}]
          }
        },
        {
          "name": "T2",
          "type": "Sink",
          "assignment": {
            "target_type": "Face",
            "targets": [{"object": "Pad1", "face": "ZMax"}]
          }
        }
      ]
    },
    {
      "name": "Ground",
      "objects": ["GND_Plane"],
      "is_ground_reference": true,
      "terminals": []
    },
    {
      "name": "Signal2",
      "objects": ["Trace2"],
      "terminals": [
        {
          "name": "T3",
          "type": "Source",
          "assignment": {
            "target_type": "Face",
            "targets": [{"object": "Trace2", "face": "XMin"}]
          }
        },
        {
          "name": "T4",
          "type": "Sink",
          "assignment": {
            "target_type": "Face",
            "targets": [{"object": "Trace2", "face": "XMax"}]
          }
        }
      ]
    }
  ]
}
```

**网络属性**：

| 字段 | 类型 | 说明 |
|------|------|------|
| `name` | String | 网络名称（如 "Signal1"、"VDD"、"GND"） |
| `objects` | Vec\<String\> | 属于该网络的几何对象列表 |
| `is_ground_reference` | bool | 是否为接地参考网络（默认 false） |
| `terminals` | Vec\<Terminal\> | 该网络上的源/汇端子列表 |

> **参考**：在 AEDT 中，HFSS 没有"Net"概念（端口直接分配到面），而 Q3D 以 Net 为中心组织提取——所有 RLCG 矩阵都是按 Net/Terminal 索引的。

### 4.7 Mesh Operations（网格操作）

参考 [HFSS 网格设置](https://www.ansys.com/training-center/course-catalog/electronics/ansys-hfss-3d-components-boundary-conditions-ports-and-mesh)：

```json
[
  {
    "name": "PatchMesh",
    "type": "LengthBased",
    "assignment": {
      "target_type": "Object",
      "targets": ["Patch"]
    },
    "properties": {
      "max_element_length": "2mm",
      "restrict_max_length": true
    }
  },
  {
    "name": "SkinMesh",
    "type": "SkinDepth",
    "assignment": {
      "target_type": "Object",
      "targets": ["copper_box"]
    },
    "properties": {
      "skin_depth": "0.02mm",
      "num_layers": 3
    }
  },
  {
    "name": "CurveMesh",
    "type": "CurvatureBased",
    "assignment": {
      "target_type": "Object",
      "targets": ["cylinder_1"]
    },
    "properties": {
      "normal_deviation": "15deg"
    }
  }
]
```

**网格操作类型**：`LengthBased` | `SkinDepth` | `CurvatureBased` | `ModelResolution`

### 4.8 Analysis Setups（分析设置）

参考 [PyAEDT create_setup](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.hfss.Hfss.create_setup.html) 和 [PyAEDT Q3D Setup Templates](https://aedt.docs.pyansys.com/version/stable/API/SetupTemplatesQ3D.html)：

**HFSS 分析设置**：

```json
[
  {
    "name": "Setup1",
    "enabled": true,
    "solution_frequency": "2.4GHz",
    "max_passes": 15,
    "max_delta_s": 0.02,
    "min_converged_passes": 2,
    "order_basis": "Mixed",
    "solver_type": "Direct",
    "frequency_sweeps": [
      {
        "name": "Sweep1",
        "type": "Interpolating",
        "start": "1GHz",
        "stop": "4GHz",
        "step": "0.01GHz",
        "save_fields": true,
        "save_rad_fields": true
      }
    ],
    "mesh_refinement": {
      "initial_mesh_settings": {
        "lambda_target": 0.3333,
        "use_default_lambda": true
      }
    }
  }
]
```

**Q3D 分析设置**：

```json
[
  {
    "name": "Q3D_Setup1",
    "enabled": true,
    "solution_type": "Q3D_ACRL",
    "adaptive_frequency": "1GHz",
    "max_passes": 10,
    "max_delta_energy": 0.02,
    "min_converged_passes": 2,
    "percent_refinement": 30,
    "solver_type": "Direct",
    "frequency_sweeps": [
      {
        "name": "Sweep1",
        "type": "Discrete",
        "start": "10MHz",
        "stop": "5GHz",
        "count": 50,
        "scale": "Logarithmic"
      }
    ],
    "dc_settings": {
      "compute_dc_resistance": true,
      "compute_dc_inductance": true
    }
  }
]
```

**Q3D 分析设置属性**：

| 字段 | 类型 | 说明 |
|------|------|------|
| `solution_type` | String | Q3D 求解类型：`Q3D_DCRL` / `Q3D_ACRL` / `Q3D_C` / `Q3D_CG` |
| `adaptive_frequency` | String | 自适应加密的参考频率 |
| `max_delta_energy` | f64 | 能量收敛阈值（替代 HFSS 的 max_delta_s） |
| `percent_refinement` | u32 | 每轮自适应加密的网格加密百分比 |
| `dc_settings` | Object | DC 提取控制（仅 DCRL 和 ACRL） |

**频率扫描 `type`**（HFSS 和 Q3D 共用）：

| 类型 | 说明 |
|------|------|
| `Discrete` | 逐频点求解，精度最高，速度最慢 |
| `Interpolating` | 自适应插值，精度与速度均衡（推荐） |
| `Fast` | 基于展开的快速扫描，速度最快 |

### 4.9 Design Settings（设计全局设置）

参考 [HFSS Design Settings](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/Modeler/DesignSettingsforHFSS.htm)：Design Settings 控制整个设计的求解行为全局开关，独立于各个 Analysis Setup。

```json
{
  "design_settings": {
    "port_impedance_normalization": {
      "enabled": true,
      "reference_impedance": "50ohm"
    },
    "deembedding": {
      "enabled": false,
      "default_distance": "0mm"
    },
    "s_matrix_type": "Modal",
    "environment_temperature": "22cel",
    "include_gravity": false,
    "model_validation": {
      "validate_before_solve": true,
      "check_intersections": true,
      "check_duplicate_boundaries": true,
      "check_port_on_boundary": true
    },
    "solver_options": {
      "use_shell_elements": false,
      "curved_elements_order": "Mixed",
      "allow_solver_fallback": true
    }
  }
}
```

| 设置项 | 说明 |
|--------|------|
| `port_impedance_normalization` | 端口阻抗归一化（默认 50Ω） |
| `deembedding` | 全局 de-embedding 默认值 |
| `s_matrix_type` | S 矩阵类型：`Modal` / `Terminal` |
| `environment_temperature` | 环境温度（影响温度依赖材料） |
| `model_validation` | 求解前自动检查项开关 |
| `solver_options` | 求解器高级选项 |

### 4.10 Radiation（辐射设置）

参考 [HFSS Post Processing](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/ReportsandPostProc/PostProcessingandGeneratingReports.htm)：Radiation 节点定义远场和近场的**采样设置**，是生成方向图、增益、轴比等辐射参数报告的前提。这与 Boundaries 中的辐射边界不同——辐射边界定义电磁场如何被吸收，Radiation 节点定义**在哪些方向/位置采样场数据**。

```json
{
  "radiation": {
    "far_field_setups": [
      {
        "name": "InfiniteSphere1",
        "type": "InfiniteSphere",
        "coordinate_system": "Global",
        "theta": {
          "start": "0deg",
          "stop": "180deg",
          "step": "1deg"
        },
        "phi": {
          "start": "0deg",
          "stop": "360deg",
          "step": "1deg"
        },
        "use_custom_radiation_surface": false,
        "radiation_surface": null
      },
      {
        "name": "HemiSphere_Upper",
        "type": "InfiniteSphere",
        "coordinate_system": "Global",
        "theta": {
          "start": "0deg",
          "stop": "90deg",
          "step": "2deg"
        },
        "phi": {
          "start": "0deg",
          "stop": "360deg",
          "step": "2deg"
        },
        "use_custom_radiation_surface": false,
        "radiation_surface": null
      }
    ],
    "near_field_setups": [
      {
        "name": "NearFieldLine1",
        "type": "Line",
        "start_point": [0, 0, 10],
        "end_point": [100, 0, 10],
        "num_points": 201
      },
      {
        "name": "NearFieldRect1",
        "type": "Rectangle",
        "center": [0, 0, 50],
        "width": 200,
        "height": 200,
        "axis": "Z",
        "num_points_u": 101,
        "num_points_v": 101
      }
    ],
    "antenna_parameters": {
      "reference_impedance": "50ohm",
      "calculate_antenna_params": true
    }
  }
}
```

**远场设置类型**：

| 类型 | 说明 |
|------|------|
| `InfiniteSphere` | 标准远场球面（theta-phi 采样），最常用 |
| `InfinitePlane` | 远场平面采样 |

**近场设置类型**：

| 类型 | 说明 |
|------|------|
| `Line` | 沿直线采样近场（近场分布扫描） |
| `Rectangle` | 在矩形平面上采样近场 |
| `Sphere` | 在球面上采样近场（用于近远场变换） |

### 4.11 Output Variables（输出变量）

参考 [HFSS Output Variables](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS3DLayout/Content/ReportsandPostProc/SpecifyingOutputVariables.htm)：Output Variables 是从仿真结果派生的**命名数学表达式**。与设计变量（`variables`）不同，输出变量的值不是用户输入的，而是**由求解结果计算**得到的。它们被 Reports 和 Optimetrics 引用，作为观察量或优化目标。

```json
{
  "output_variables": [
    {
      "name": "S11_at_center",
      "expression": "dB(S(Port1,Port1))",
      "description": "中心频率回波损耗"
    },
    {
      "name": "BW_10dB",
      "expression": "bandwidth(S(Port1,Port1), -10)",
      "description": "10dB 回波损耗带宽"
    },
    {
      "name": "PeakGain",
      "expression": "max(GainTotal)",
      "description": "峰值增益"
    },
    {
      "name": "Impedance_Real",
      "expression": "re(Z(Port1,Port1))",
      "description": "输入阻抗实部"
    },
    {
      "name": "Efficiency",
      "expression": "RadiatedPower / IncidentPower",
      "description": "天线辐射效率"
    }
  ]
}
```

> **与设计变量的区别**：设计变量（`$freq`、`patch_l`）是**输入参数**，驱动几何/材料/设置；输出变量是**结果观察量**，从 S 参数、场数据等求解结果中提取，可用于 Report 绘图或 Optimetrics 优化目标。

**Q3D 输出变量示例**：

```json
{
  "output_variables": [
    {
      "name": "SelfR_Signal1",
      "expression": "R(Signal1:T1, Signal1:T1)",
      "description": "Signal1 自电阻"
    },
    {
      "name": "SelfL_Signal1",
      "expression": "L(Signal1:T1, Signal1:T1)",
      "description": "Signal1 自感"
    },
    {
      "name": "MutualL_12",
      "expression": "L(Signal1:T1, Signal2:T3)",
      "description": "Signal1 与 Signal2 互感"
    },
    {
      "name": "CouplingCoeff_12",
      "expression": "L(Signal1:T1, Signal2:T3) / sqrt(L(Signal1:T1, Signal1:T1) * L(Signal2:T3, Signal2:T3))",
      "description": "Signal1 与 Signal2 电感耦合系数"
    },
    {
      "name": "MutualC_12",
      "expression": "abs(C(Signal1:T1, Signal2:T3))",
      "description": "Signal1 与 Signal2 互电容"
    },
    {
      "name": "SelfC_Signal1",
      "expression": "C(Signal1:T1, Signal1:T1)",
      "description": "Signal1 自电容（到地总电容）"
    },
    {
      "name": "R_DC_Signal1",
      "expression": "R_DC(Signal1:T1, Signal1:T1)",
      "description": "Signal1 DC 电阻"
    },
    {
      "name": "CharImpedance_Signal1",
      "expression": "sqrt(L(Signal1:T1, Signal1:T1) / C(Signal1:T1, Signal1:T1))",
      "description": "Signal1 特征阻抗估算"
    }
  ]
}
```

> **Q3D 量表达式体系**：Q3D 输出变量使用 `R(net:term, net:term)`、`L()`、`C()`、`G()` 作为基础量函数，以端子名称为索引。还支持 `R_DC()` 和 `L_DC()` 引用 DC 提取值。这些量表达式可以在 Reports 和 Optimetrics 中被引用，例如优化目标"最小化信号间耦合系数"或参数扫描"走线间距对互电容的影响"。

### 4.12 Field Overlays（场叠加显示）

参考 [HFSS Field Overlays](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/HFSS3DLayout/Content/3DLayout/PlottingFieldOverlays.htm)：Field Overlays 定义在 3D 模型上叠加显示电磁场的方案。这些定义保存在工程文件中，打开工程后可直接复现场分布可视化，无需重新配置。

```json
{
  "field_overlays": [
    {
      "name": "E_Field_Top",
      "quantity": "E",
      "component": "Mag",
      "plot_type": "Surface",
      "assignment": {
        "type": "Face",
        "targets": [{ "object": "Substrate", "face": "ZMax" }]
      },
      "solution": "Setup1 : LastAdaptive",
      "frequency": "$freq",
      "phase": "0deg",
      "scale": {
        "type": "Log",
        "min": 1.0,
        "max": 1000.0,
        "unit": "V/m"
      },
      "display": {
        "plot_style": "Shaded",
        "show_arrows": false,
        "num_colors": 256,
        "opacity": 0.8
      }
    },
    {
      "name": "J_Surf_Patch",
      "quantity": "Jsurf",
      "component": "Mag",
      "plot_type": "Surface",
      "assignment": {
        "type": "Object",
        "targets": ["Patch"]
      },
      "solution": "Setup1 : LastAdaptive",
      "frequency": "$freq",
      "phase": "0deg",
      "scale": {
        "type": "Linear",
        "min": null,
        "max": null,
        "unit": "A/m"
      },
      "display": {
        "plot_style": "Arrow",
        "show_arrows": true,
        "arrow_spacing": 5,
        "num_colors": 256,
        "opacity": 1.0
      }
    },
    {
      "name": "E_Field_CrossSection",
      "quantity": "E",
      "component": "Mag",
      "plot_type": "CutPlane",
      "assignment": {
        "type": "Plane",
        "normal": "Y",
        "position": 30.0
      },
      "solution": "Setup1 : Sweep1",
      "frequency": "2.4GHz",
      "phase": "0deg",
      "scale": {
        "type": "Log",
        "min": null,
        "max": null,
        "unit": "V/m"
      },
      "display": {
        "plot_style": "Shaded",
        "show_arrows": false,
        "num_colors": 256,
        "opacity": 0.9
      }
    }
  ]
}
```

**场量 `quantity`** 枚举：

| 量 | 说明 |
|----|------|
| `E` | 电场强度 |
| `H` | 磁场强度 |
| `Jvol` | 体电流密度 |
| `Jsurf` | 面电流密度 |
| `SAR` | 比吸收率 |
| `Poynting` | 坡印廷矢量 |
| `ChargeDistribution` | 电荷面密度（Q3D 电容提取） |
| `OhmicLoss` | 欧姆损耗密度（Q3D 电阻提取） |

**分量 `component`** 枚举：`Mag` | `MagX` | `MagY` | `MagZ` | `Real` | `Imag` | `Vector`

**绘图类型 `plot_type`**：`Surface` | `CutPlane` | `Volume` | `Line`

**Q3D 场叠加示例**：

```json
{
  "field_overlays": [
    {
      "name": "J_Current_Signal1",
      "quantity": "Jvol",
      "component": "Mag",
      "plot_type": "Surface",
      "assignment": {
        "type": "Object",
        "targets": ["Trace1", "Via1", "Pad1"]
      },
      "solution": "Q3D_Setup1 : LastAdaptive",
      "frequency": "1GHz",
      "phase": "0deg",
      "scale": {
        "type": "Log",
        "min": 100.0,
        "max": 1e7,
        "unit": "A/m2"
      },
      "display": {
        "plot_style": "Shaded",
        "show_arrows": true,
        "arrow_spacing": 3,
        "num_colors": 256,
        "opacity": 1.0
      }
    },
    {
      "name": "E_Field_Dielectric",
      "quantity": "E",
      "component": "Mag",
      "plot_type": "CutPlane",
      "assignment": {
        "type": "Plane",
        "normal": "Z",
        "position": 0.8
      },
      "solution": "Q3D_Setup1 : LastAdaptive",
      "frequency": "1GHz",
      "phase": "0deg",
      "scale": {
        "type": "Log",
        "min": null,
        "max": null,
        "unit": "V/m"
      },
      "display": {
        "plot_style": "Shaded",
        "show_arrows": false,
        "num_colors": 256,
        "opacity": 0.85
      }
    },
    {
      "name": "Charge_Distribution",
      "quantity": "ChargeDistribution",
      "component": "Mag",
      "plot_type": "Surface",
      "assignment": {
        "type": "Object",
        "targets": ["Trace1", "Trace2", "GND_Plane"]
      },
      "solution": "Q3D_Setup1 : LastAdaptive",
      "frequency": "1GHz",
      "phase": "0deg",
      "scale": {
        "type": "Linear",
        "min": null,
        "max": null,
        "unit": "C/m2"
      },
      "display": {
        "plot_style": "Shaded",
        "show_arrows": false,
        "num_colors": 256,
        "opacity": 1.0
      }
    }
  ]
}
```

> **Q3D 场叠加典型用途**：
> - **电流密度 (Jvol)**：识别电流拥挤区域（如弯折处、过孔连接处），分析趋肤效应和邻近效应在高频下对电流分布的影响。
> - **电场 (E)**：分析导体间电场耦合路径，识别电容耦合的主要贡献区域。
> - **电荷分布 (ChargeDistribution)**：直观显示导体表面电荷分布，帮助理解电容矩阵结果的物理来源。
> - **欧姆损耗 (OhmicLoss)**：识别热点区域，评估导体损耗对整体电阻的贡献。
>
> 参考 [Plotting Field Overlays in Q3D](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Q3DExtractor/Content/ReportsandPostProc/PlottingFieldOverlaysinQ3D.htm)。

### 4.13 Optimetrics（参数化分析与优化）

参考 [HFSS Parametric Overview](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/HFSS3DLayout/Content/Optimetrics/ParametricOverview.htm)：Optimetrics 驱动设计变量在指定范围内变化，自动执行多次仿真。支持参数扫描、优化、灵敏度分析等多种分析类型。

```json
{
  "optimetrics": [
    {
      "name": "LengthSweep",
      "type": "ParametricSweep",
      "enabled": true,
      "setup": "Setup1",
      "sweep_definitions": [
        {
          "variable": "patch_l",
          "type": "LinearStep",
          "start": "25mm",
          "stop": "32mm",
          "step": "0.5mm"
        }
      ],
      "constraints": [],
      "goals": []
    },
    {
      "name": "FreqEpsSweep",
      "type": "ParametricSweep",
      "enabled": true,
      "setup": "Setup1",
      "sweep_definitions": [
        {
          "variable": "$freq",
          "type": "DiscreteList",
          "values": ["2.0GHz", "2.4GHz", "2.8GHz"]
        },
        {
          "variable": "$eps_sub",
          "type": "LinearStep",
          "start": "3.5",
          "stop": "5.0",
          "step": "0.5"
        }
      ],
      "constraints": [],
      "goals": []
    },
    {
      "name": "MatchOptimize",
      "type": "Optimization",
      "enabled": true,
      "setup": "Setup1",
      "algorithm": "QuasiNewton",
      "max_iterations": 50,
      "variables": [
        {
          "variable": "patch_l",
          "min": "20mm",
          "max": "35mm",
          "starting": "28.5mm"
        },
        {
          "variable": "patch_w",
          "min": "25mm",
          "max": "45mm",
          "starting": "37mm"
        }
      ],
      "goals": [
        {
          "name": "MinS11",
          "expression": "S11_at_center",
          "condition": "Minimize",
          "weight": 1.0
        },
        {
          "name": "S11_below_10dB",
          "expression": "dB(S(Port1,Port1))",
          "condition": "LessThan",
          "target": -10.0,
          "frequency_range": { "start": "2.3GHz", "stop": "2.5GHz" },
          "weight": 2.0
        }
      ],
      "constraints": [
        {
          "expression": "patch_l",
          "condition": "LessThan",
          "target_expression": "patch_w"
        }
      ]
    },
    {
      "name": "SensitivityCheck",
      "type": "Sensitivity",
      "enabled": true,
      "setup": "Setup1",
      "variables": [
        {
          "variable": "patch_l",
          "variation": "5%",
          "distribution": "Uniform"
        },
        {
          "variable": "$eps_sub",
          "variation": "10%",
          "distribution": "Gaussian"
        }
      ],
      "output": "S11_at_center",
      "num_samples": 100
    }
  ]
}
```

**Optimetrics 类型**：

| type | 说明 | 关键参数 |
|------|------|---------|
| `ParametricSweep` | 参数扫描 | `sweep_definitions`（变量遍历规则） |
| `Optimization` | 自动优化 | `algorithm`, `goals`, `constraints` |
| `Sensitivity` | 灵敏度分析 | `variation`, `distribution`, `num_samples` |
| `Statistical` | 统计分析 | `distribution`, `num_trials` |
| `Tuning` | 交互式调参 | `variables`（实时调整） |

**扫描定义 `type`**：`LinearStep` | `LinearCount` | `LogScale` | `DiscreteList`

**优化算法 `algorithm`**：`QuasiNewton` | `PatternSearch` | `GeneticAlgorithm` | `SNLP` | `Bayesian`

### 4.14 Reports（后处理报告）

参考 [HFSS Post Processing](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/ReportsandPostProc/PostProcessingandGeneratingReports.htm)：Reports 定义如何从仿真结果中提取和展示数据。

```json
{
  "reports": [
    {
      "name": "S11_Return_Loss",
      "category": "SParameter",
      "chart_type": "Rectangular",
      "solution": "Setup1 : Sweep1",
      "domain": {
        "type": "Frequency",
        "primary_sweep": "Freq"
      },
      "traces": [
        {
          "name": "dB(S(Port1,Port1))",
          "expression": "dB(S(Port1,Port1))",
          "style": { "color": [0, 0, 255], "line_width": 2, "line_style": "Solid" }
        }
      ],
      "x_axis": {
        "label": "Frequency",
        "unit": "GHz",
        "auto_range": true
      },
      "y_axis": {
        "label": "S11 (dB)",
        "unit": "dB",
        "min": -30.0,
        "max": 0.0
      },
      "markers": [
        {
          "name": "Center",
          "trace": "dB(S(Port1,Port1))",
          "x_value": "2.4GHz"
        }
      ],
      "limit_lines": [
        {
          "name": "-10dB_Threshold",
          "y_value": -10.0,
          "style": { "color": [255, 0, 0], "line_style": "Dashed" }
        }
      ]
    },
    {
      "name": "Smith_Chart",
      "category": "SParameter",
      "chart_type": "Smith",
      "solution": "Setup1 : Sweep1",
      "domain": {
        "type": "Frequency",
        "primary_sweep": "Freq"
      },
      "traces": [
        {
          "name": "S(Port1,Port1)",
          "expression": "S(Port1,Port1)"
        }
      ],
      "markers": [
        { "name": "f0", "trace": "S(Port1,Port1)", "x_value": "2.4GHz" }
      ]
    },
    {
      "name": "Radiation_Pattern",
      "category": "FarField",
      "chart_type": "Polar",
      "solution": "Setup1 : LastAdaptive",
      "domain": {
        "type": "Angle",
        "primary_sweep": "Theta",
        "fixed_values": { "Freq": "$freq", "Phi": "0deg" }
      },
      "far_field_setup": "InfiniteSphere1",
      "traces": [
        {
          "name": "GainTotal_E",
          "expression": "GainTotal",
          "fixed_values": { "Phi": "0deg" }
        },
        {
          "name": "GainTotal_H",
          "expression": "GainTotal",
          "fixed_values": { "Phi": "90deg" }
        }
      ],
      "y_axis": {
        "label": "Gain (dBi)",
        "unit": "dBi",
        "min": -20.0,
        "max": 10.0
      }
    },
    {
      "name": "Parametric_Comparison",
      "category": "SParameter",
      "chart_type": "Rectangular",
      "solution": "Setup1 : Sweep1",
      "domain": {
        "type": "Frequency",
        "primary_sweep": "Freq"
      },
      "traces": [
        {
          "name": "S11_nominal",
          "expression": "dB(S(Port1,Port1))",
          "parametric_values": { "patch_l": "28.5mm" }
        },
        {
          "name": "S11_short",
          "expression": "dB(S(Port1,Port1))",
          "parametric_values": { "patch_l": "26mm" }
        },
        {
          "name": "S11_long",
          "expression": "dB(S(Port1,Port1))",
          "parametric_values": { "patch_l": "31mm" }
        }
      ]
    }
  ]
}
```

**报告类别 `category`**：

| 类别 | 可用数据 |
|------|---------|
| `SParameter` | S/Y/Z 参数、VSWR、群时延 |
| `FarField` | 增益、方向性、轴比、极化、效率 |
| `NearField` | 近场 E/H 分布 |
| `Fields` | 场量（从 Field Calculator 表达式） |
| `Eigenmode` | 谐振频率、Q 值 |
| `Emission` | 辐射功率、EMI |
| `RLCGMatrix` | Q3D：RLCG 矩阵元素 vs 频率 |
| `Q3DFields` | Q3D：电流密度、电场、电荷分布 |

**图表类型 `chart_type`**：

| 类型 | 说明 |
|------|------|
| `Rectangular` | 直角坐标（X-Y 图） |
| `Polar` | 极坐标（方向图） |
| `Smith` | 史密斯圆图（阻抗匹配） |
| `DataTable` | 数据表格 |
| `Polar3D` | 3D 极坐标（立体方向图） |
| `MatrixTable` | Q3D 矩阵表格（RLCG 矩阵 + 热力图着色） |

**Q3D 报告示例**：

```json
{
  "reports": [
    {
      "name": "R_vs_Frequency",
      "category": "RLCGMatrix",
      "chart_type": "Rectangular",
      "solution": "Q3D_Setup1 : Sweep1",
      "domain": {
        "type": "Frequency",
        "primary_sweep": "Freq"
      },
      "traces": [
        {
          "name": "R_self_Signal1",
          "expression": "R(Signal1:T1, Signal1:T1)",
          "style": { "color": [0, 0, 255], "line_width": 2, "line_style": "Solid" }
        },
        {
          "name": "R_self_Signal2",
          "expression": "R(Signal2:T3, Signal2:T3)",
          "style": { "color": [255, 0, 0], "line_width": 2, "line_style": "Solid" }
        },
        {
          "name": "R_mutual_12",
          "expression": "R(Signal1:T1, Signal2:T3)",
          "style": { "color": [0, 128, 0], "line_width": 1, "line_style": "Dashed" }
        }
      ],
      "x_axis": {
        "label": "Frequency",
        "unit": "GHz",
        "scale": "Logarithmic",
        "auto_range": true
      },
      "y_axis": {
        "label": "Resistance (Ω)",
        "unit": "ohm",
        "auto_range": true
      }
    },
    {
      "name": "L_Matrix_Table",
      "category": "RLCGMatrix",
      "chart_type": "MatrixTable",
      "solution": "Q3D_Setup1 : Sweep1",
      "domain": {
        "type": "Frequency",
        "primary_sweep": "Freq",
        "fixed_values": { "Freq": "1GHz" }
      },
      "matrix_type": "L",
      "display_options": {
        "heatmap_enabled": true,
        "show_unit": true,
        "decimal_places": 4
      }
    },
    {
      "name": "Coupling_Coefficient",
      "category": "RLCGMatrix",
      "chart_type": "Rectangular",
      "solution": "Q3D_Setup1 : Sweep1",
      "domain": {
        "type": "Frequency",
        "primary_sweep": "Freq"
      },
      "traces": [
        {
          "name": "k_12",
          "expression": "L(Signal1:T1, Signal2:T3) / sqrt(L(Signal1:T1, Signal1:T1) * L(Signal2:T3, Signal2:T3))"
        }
      ],
      "y_axis": {
        "label": "Coupling Coefficient",
        "unit": "",
        "min": 0.0,
        "max": 1.0
      }
    }
  ]
}
```

> **Q3D 报告说明**：Q3D 的 `RLCGMatrix` 类别报告支持 `Rectangular`（曲线图）和 `MatrixTable`（矩阵数据表格）两种图表类型。曲线图适合观察 RLCG 参数随频率的变化趋势（如趋肤效应导致的电阻增长），矩阵表格适合在特定频率下查看完整的耦合矩阵。参考 [Q3D Creating Reports](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/ReportsandPostProc/CreatingReportsQ3D.htm)。

---

## 5. 仿真生成文件与结果跟踪

参考 [Ansys Electronics Desktop Files](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/Variables/ANSYSElectronicsDesktopFiles.htm)：HFSS 的仿真结果**不存储在 `.aedt` 工程文件中**，而是放在独立的 `.aedtresults/` 目录下。工程文件只记录"指向哪里找结果"的索引信息（`.asol` 求解变量数据库）。EMStudio 采用相同策略。

### 5.1 文件生成时间线

| 阶段 | 生成的文件 | 位置 | 格式 | 用途 |
|------|-----------|------|------|------|
| **工程创建** | `*.emsp` | 工程目录 | JSON 文本 | 工程定义（几何/材料/设置等） |
| **工程打开** | `*.emsp.lock` | 工程目录 | JSON 文本 | 并发写入保护 |
| **定时/手动** | `*.emsp.auto` | 工程目录 | JSON 文本 | 自动恢复备份 |
| **验证阶段** | `validation_report.json` | results 目录 | JSON 文本 | 模型检查结果（错误/警告） |
| **初始网格** | `initial_mesh.msh` | results/Setup/mesh/ | Gmsh MSH 4.1 | 初始四面体网格 |
| **每轮自适应** | `pass_N_mesh.msh` | results/Setup/mesh/ | Gmsh MSH 4.1 | 第 N 轮细化后的网格 |
| **每轮自适应** | `pass_N_solution.bin` | results/Setup/ | 二进制 | 第 N 轮的场解（FEM 系数） |
| **每轮自适应** | `convergence.json` | results/Setup/ | JSON 文本 | 收敛历史（实时追加） |
| **自适应完成** | `mesh_stats.json` | results/Setup/ | JSON 文本 | 最终网格统计 |
| **自适应完成** | `s_parameters.json` | results/Setup/ | JSON 文本 | 自适应频率处的 S 参数 |
| **频率扫描** | `sweep_s_parameters.json` | results/Setup/Sweep/ | JSON 文本 | 全频段 S 参数 |
| **频率扫描** | `s_parameters.snp` | results/Setup/Sweep/ | Touchstone | 标准 Touchstone 格式导出 |
| **场保存** | `e_field_*.bin`, `h_field_*.bin` | results/Setup/fields/ | 二进制 | 每个频点的电/磁场数据 |
| **远场计算** | `far_field_*.json` | results/Setup/far_field/ | JSON 文本 | 方向图、增益、轴比等 |
| **近场计算** | `near_field_*.json` | results/Setup/near_field/ | JSON 文本 | 近场采样数据 |
| **求解器日志** | `solver.log` | results/design/ | 纯文本 | 求解器消息、警告、错误 |
| **性能剖面** | `profile.json` | results/Setup/ | JSON 文本 | 各阶段耗时、内存峰值 |
| **Optimetrics** | `variation_*/` | results/Setup/ | 子目录 | 每组参数变量值的独立结果 |
| **报告导出** | `report_*.csv` | results/exports/ | CSV | 手动或自动导出的报告数据 |
| **场图导出** | `field_plot_*.png` | results/exports/ | PNG | 场叠加截图 |

**Q3D 准静态求解额外生成的文件**：

| 阶段 | 生成的文件 | 位置 | 格式 | 用途 |
|------|-----------|------|------|------|
| **自适应完成** | `rlcg_matrix.json` | results/Setup/ | JSON 文本 | 自适应频率处的 RLCG 矩阵 |
| **频率扫描** | `rlcg_matrix.json` | results/Setup/Sweep/ | JSON 文本 | 全频段 RLCG 矩阵 |
| **场保存** | `j_field_*.bin` | results/Setup/fields/ | 二进制 | Q3D 电流密度场数据 |
| **场保存** | `charge_*.bin` | results/Setup/fields/ | 二进制 | Q3D 电荷分布数据 |
| **等效电路导出** | `equivalent_circuit.sp` | results/exports/ | SPICE | SPICE 等效电路网表 |
| **等效电路导出** | `circuit_export_config.json` | results/exports/ | JSON 文本 | 导出配置元数据 |
| **S 参数转换** | `s_parameters_from_rlcg.snp` | results/Setup/Sweep/ | Touchstone | RLCG → S 参数转换结果 |

### 5.2 结果目录完整结构

```
MyAntenna.emsp                          # 工程文件（JSON）
MyAntenna.emsp.lock                     # 打开锁
MyAntenna.emsp.auto                     # 自动保存

MyAntenna.emsp.results/                 # 结果根目录
 ├── solve_log.txt                      # 全局求解日志
 │
 ├── design-001/                        # 设计级目录
 │    ├── validation_report.json        # 模型验证报告
 │    ├── solver.log                    # 设计级求解日志
 │    │
 │    ├── Setup1/                       # Analysis Setup 级目录
 │    │    ├── convergence.json         # 自适应收敛历史
 │    │    ├── mesh_stats.json          # 最终网格统计
 │    │    ├── profile.json             # 性能剖面（耗时/内存）
 │    │    ├── s_parameters.json        # 自适应频率处 S 参数
 │    │    │
 │    │    ├── mesh/                    # 网格数据
 │    │    │    ├── initial_mesh.msh    # 初始网格（Gmsh MSH 4.1）
 │    │    │    ├── pass_1_mesh.msh     # 第 1 轮自适应网格
 │    │    │    ├── pass_2_mesh.msh     # 第 2 轮
 │    │    │    └── final_mesh.msh      # 最终收敛网格
 │    │    │
 │    │    ├── solutions/              # 场解数据（FEM 系数矩阵）
 │    │    │    ├── pass_1_solution.bin
 │    │    │    ├── pass_2_solution.bin
 │    │    │    └── final_solution.bin
 │    │    │
 │    │    ├── fields/                 # 导出/缓存的场数据
 │    │    │    ├── e_field_2.4GHz.bin # 电场（频点×分量）
 │    │    │    └── h_field_2.4GHz.bin # 磁场
 │    │    │
 │    │    ├── far_field/              # 远场数据
 │    │    │    ├── InfiniteSphere1_2.4GHz.json  # 方向图
 │    │    │    └── antenna_params_2.4GHz.json   # 天线参数
 │    │    │
 │    │    ├── near_field/             # 近场数据
 │    │    │    └── NearFieldLine1_2.4GHz.json
 │    │    │
 │    │    └── Sweep1/                 # 频率扫描结果
 │    │         ├── s_parameters.json  # 全频段 S 参数（JSON）
 │    │         ├── s_parameters.s1p   # Touchstone 格式
 │    │         ├── fields/            # 扫频场数据（如 save_fields=true）
 │    │         └── far_field/         # 扫频远场
 │    │
 │    ├── Setup1__Optimetrics/         # Optimetrics 结果
 │    │    ├── LengthSweep/            # 参数扫描
 │    │    │    ├── variation_001/      # patch_l=25.0mm
 │    │    │    │    ├── convergence.json
 │    │    │    │    ├── s_parameters.json
 │    │    │    │    └── ...
 │    │    │    ├── variation_002/      # patch_l=25.5mm
 │    │    │    └── summary.json       # 所有变量组合的汇总
 │    │    └── MatchOptimize/          # 优化
 │    │         ├── iteration_001/
 │    │         ├── iteration_002/
 │    │         ├── convergence.json   # 优化收敛历史
 │    │         └── best_result.json   # 最优解
 │    │
 │    ├── snapshots/                   # 3D 视图截图缓存
 │    │    └── viewport_001.png
 │    │
 │    └── exports/                     # 用户手动导出
 │         ├── report_S11.csv
 │         └── field_E_top.png
 │
 └── design-002/                       # 多设计各自独立
      └── ...
```

**Q3D 准静态设计的结果目录结构**（与 HFSS 共享框架，但内容不同）：

```
MyPCB_Parasitics.emsp.results/
 └── design-002/                        # Q3D 准静态设计
      ├── validation_report.json        # 模型验证（含 Net/Terminal 检查）
      ├── solver.log                    # 求解日志
      │
      ├── Q3D_Setup1/                   # Q3D 分析设置
      │    ├── convergence.json         # 自适应收敛（delta_energy + rlcg_snapshot）
      │    ├── mesh_stats.json          # MoM 面网格统计
      │    ├── profile.json             # 性能剖面
      │    ├── rlcg_matrix.json         # 自适应频率处的 RLCG 矩阵
      │    │
      │    ├── mesh/                    # 网格数据
      │    │    ├── initial_mesh.msh    # 初始面网格
      │    │    ├── pass_1_mesh.msh
      │    │    └── final_mesh.msh      # 最终收敛面网格
      │    │
      │    ├── fields/                  # Q3D 场数据
      │    │    ├── j_field_1GHz.bin    # 电流密度场
      │    │    ├── e_field_1GHz.bin    # 电场
      │    │    └── charge_1GHz.bin     # 电荷分布
      │    │
      │    └── Sweep1/                  # 频率扫描
      │         ├── rlcg_matrix.json    # 全频段 RLCG 矩阵
      │         ├── s_parameters_from_rlcg.snp  # RLCG→S 参数转换
      │         └── fields/             # 扫频场数据（按频点）
      │
      ├── Q3D_Setup1__Optimetrics/      # Q3D Optimetrics
      │    └── SpacingSweep/            # 走线间距扫描
      │         ├── variation_001/       # spacing=0.1mm
      │         │    ├── convergence.json
      │         │    └── rlcg_matrix.json
      │         ├── variation_002/       # spacing=0.2mm
      │         └── summary.json
      │
      └── exports/                      # 导出文件
           ├── equivalent_circuit.sp     # SPICE 等效电路
           ├── circuit_export_config.json
           ├── report_R_vs_freq.csv
           └── report_L_matrix.csv
```

### 5.3 工程文件中的结果索引

结果数据存在 `.emsp.results/` 目录中，但工程文件 `.emsp` 需要记录**哪些结果已存在**以及**结果是否过期**，这对应 HFSS 中 `.asol` 文件的角色。

在 Design JSON 中增加 `solution_index` 字段：

```json
{
  "solution_index": {
    "last_solve_time": "2026-04-04T14:30:00Z",
    "setups": {
      "Setup1": {
        "status": "Converged",
        "solved_at": "2026-04-04T14:28:00Z",
        "converged_pass": 3,
        "num_tetrahedra": 10250,
        "is_stale": false,
        "solved_variations": {
          "nominal": {
            "variables": {},
            "result_path": "design-001/Setup1/"
          }
        },
        "sweeps": {
          "Sweep1": {
            "status": "Completed",
            "num_frequency_points": 301,
            "result_path": "design-001/Setup1/Sweep1/"
          }
        }
      }
    },
    "optimetrics": {
      "LengthSweep": {
        "status": "Completed",
        "total_variations": 15,
        "completed_variations": 15,
        "result_path": "design-001/Setup1__Optimetrics/LengthSweep/"
      },
      "MatchOptimize": {
        "status": "Converged",
        "total_iterations": 23,
        "best_cost": -22.5,
        "result_path": "design-001/Setup1__Optimetrics/MatchOptimize/"
      }
    },
    "stale_reason": null
  }
}
```

**Q3D 准静态设计的结果索引**：

```json
{
  "solution_index": {
    "last_solve_time": "2026-04-04T15:30:00Z",
    "setups": {
      "Q3D_Setup1": {
        "status": "Converged",
        "solved_at": "2026-04-04T15:28:00Z",
        "solution_type": "Q3D_ACRL",
        "converged_pass": 3,
        "num_triangles": 6100,
        "num_tetrahedra": 0,
        "final_delta_energy": 0.015,
        "is_stale": false,
        "rlcg_summary": {
          "num_nets": 3,
          "num_terminals": 4,
          "r_max_ohm": 0.285,
          "l_max_nH": 2.93,
          "c_max_pF": 0.445,
          "g_max_mS": 0.0013
        },
        "solved_variations": {
          "nominal": {
            "variables": {},
            "result_path": "design-002/Q3D_Setup1/"
          }
        },
        "sweeps": {
          "Sweep1": {
            "status": "Completed",
            "num_frequency_points": 50,
            "result_path": "design-002/Q3D_Setup1/Sweep1/"
          }
        }
      }
    },
    "optimetrics": {
      "SpacingSweep": {
        "status": "Completed",
        "total_variations": 10,
        "completed_variations": 10,
        "result_path": "design-002/Q3D_Setup1__Optimetrics/SpacingSweep/"
      }
    },
    "exports": {
      "equivalent_circuit": {
        "exported_at": "2026-04-04T15:35:00Z",
        "model_type": "BroadbandLumped",
        "file_path": "design-002/exports/equivalent_circuit.sp"
      }
    },
    "stale_reason": null
  }
}
```

> **Q3D 结果索引额外字段**：
> - `solution_type`：记录 Q3D 求解类型（DCRL/ACRL/C/CG），影响哪些 RLCG 矩阵可用。
> - `num_triangles`：MoM 面网格三角形数（Q3D 主要指标，替代 HFSS 的 `num_tetrahedra`）。
> - `final_delta_energy`：能量收敛阈值（替代 HFSS 的 `max_delta_s`）。
> - `rlcg_summary`：RLCG 矩阵摘要信息，供 UI 快速显示而无需加载完整矩阵文件。
> - `exports`：记录已导出的等效电路模型信息。

**Q3D 结果过期的额外触发条件**：

| 用户操作 | Q3D 专有影响 | `is_stale` 行为 |
|---------|-------------|-----------------|
| 修改 Net 分配 | RLCG 矩阵索引变化 | 所有结果标记为 stale |
| 修改 Source/Sink 端子 | 提取通道变化 | 所有结果标记为 stale |
| 添加/移除 Net | 矩阵维度变化 | 所有结果标记为 stale |
| 修改 ThinConductor 边界 | 影响 MoM 面网格 | 标记 stale |
| 修改 Ground Net 引用 | 参考电位变化 | 所有结果标记为 stale |

| 用户操作 | 影响 | `is_stale` 行为 |
|---------|------|-----------------|
| 修改几何操作参数 | 网格和场解失效 | 该 Setup 下所有结果标记为 stale |
| 修改材料属性 | 场解失效，网格可能仍有效 | 标记 stale，保留网格 |
| 修改边界/激励 | 场解失效 | 标记 stale |
| 修改 Setup 参数 | 当前 Setup 结果失效 | 仅该 Setup 标记 stale |
| 修改变量值 | 依赖链上的结果失效 | 根据依赖图标记 |
| 修改 Report 定义 | 结果数据仍有效 | 不标记 stale，仅重新渲染报告 |
| 添加新对象（不影响已有对象） | 已有结果仍有效 | 不标记 stale |

### 5.4 收敛数据格式

```json
{
  "setup": "Setup1",
  "passes": [
    {
      "pass_number": 1,
      "num_tetrahedra": 5420,
      "max_delta_s": 0.15,
      "elapsed_time_sec": 12.3,
      "peak_memory_mb": 256
    },
    {
      "pass_number": 2,
      "num_tetrahedra": 8103,
      "max_delta_s": 0.045,
      "elapsed_time_sec": 18.7,
      "peak_memory_mb": 412
    },
    {
      "pass_number": 3,
      "num_tetrahedra": 10250,
      "max_delta_s": 0.012,
      "elapsed_time_sec": 25.1,
      "peak_memory_mb": 580
    }
  ],
  "converged": true,
  "converged_at_pass": 3,
  "total_elapsed_sec": 56.1,
  "total_peak_memory_mb": 580
}
```

### 5.5 性能剖面格式

参考 [HFSS Solution Profile](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Mechanical/Content/ReportsandPostProc/ViewingaSolutionProfile.htm)：

```json
{
  "setup": "Setup1",
  "phases": [
    {
      "name": "Validation",
      "start_time": "2026-04-04T14:25:00Z",
      "elapsed_sec": 0.5,
      "status": "OK"
    },
    {
      "name": "InitialMesh",
      "start_time": "2026-04-04T14:25:01Z",
      "elapsed_sec": 2.3,
      "num_tetrahedra": 5420,
      "peak_memory_mb": 128
    },
    {
      "name": "AdaptivePass",
      "pass_number": 1,
      "start_time": "2026-04-04T14:25:03Z",
      "sub_phases": [
        { "name": "MatrixAssembly", "elapsed_sec": 3.2, "peak_memory_mb": 210 },
        { "name": "MatrixSolve", "elapsed_sec": 8.1, "peak_memory_mb": 256 },
        { "name": "FieldRecovery", "elapsed_sec": 0.8, "peak_memory_mb": 200 },
        { "name": "ErrorEstimation", "elapsed_sec": 0.2, "peak_memory_mb": 180 }
      ],
      "elapsed_sec": 12.3,
      "peak_memory_mb": 256
    },
    {
      "name": "FrequencySweep",
      "sweep_name": "Sweep1",
      "start_time": "2026-04-04T14:26:30Z",
      "elapsed_sec": 45.0,
      "num_frequency_points": 301,
      "peak_memory_mb": 620
    },
    {
      "name": "FarFieldComputation",
      "start_time": "2026-04-04T14:27:15Z",
      "elapsed_sec": 5.2,
      "peak_memory_mb": 300
    }
  ],
  "total_elapsed_sec": 135.0,
  "total_peak_memory_mb": 620,
  "cpu_cores_used": 8
}
```

### 5.6 S 参数数据格式

```json
{
  "ports": ["Port1", "Port2"],
  "frequencies_ghz": [1.0, 1.5, 2.0, 2.4, 3.0, 4.0],
  "data": {
    "S11": {
      "magnitude_db": [-2.1, -5.3, -12.8, -25.6, -8.4, -3.2],
      "phase_deg": [170.2, 145.6, 98.3, 2.1, -85.4, -155.3]
    },
    "S21": {
      "magnitude_db": [-35.2, -28.1, -18.5, -12.3, -20.1, -30.5],
      "phase_deg": [-45.3, -92.1, -135.8, -178.2, 120.5, 60.8]
    }
  },
  "touchstone_export": "s_parameters.s2p"
}
```

### 5.7 远场数据格式

```json
{
  "setup": "InfiniteSphere1",
  "frequency": "2.4GHz",
  "coordinate_system": "Global",
  "theta_deg": [0, 1, 2, "...180"],
  "phi_deg": [0, 1, 2, "...360"],
  "gain_total_dbi": [[-3.2, -3.1, "..."], "..."],
  "gain_theta_dbi": [[-5.1, -5.0, "..."], "..."],
  "gain_phi_dbi": [[-8.3, -8.2, "..."], "..."],
  "axial_ratio_db": [[40.0, 39.5, "..."], "..."],
  "antenna_parameters": {
    "peak_gain_dbi": 7.2,
    "peak_directivity_dbi": 7.5,
    "radiation_efficiency": 0.93,
    "beamwidth_e_plane_deg": 78.0,
    "beamwidth_h_plane_deg": 85.0,
    "front_to_back_ratio_db": 15.3,
    "radiated_power_w": 0.0093,
    "accepted_power_w": 0.01,
    "incident_power_w": 0.01
  }
}
```

### 5.8 Optimetrics 汇总格式

```json
{
  "optimetrics_name": "LengthSweep",
  "type": "ParametricSweep",
  "variables_swept": ["patch_l"],
  "total_variations": 15,
  "variations": [
    {
      "index": 1,
      "variables": { "patch_l": "25.0mm" },
      "status": "Converged",
      "output_values": {
        "S11_at_center": -8.2,
        "PeakGain": 6.1,
        "BW_10dB": 0.0
      },
      "result_path": "variation_001/"
    },
    {
      "index": 2,
      "variables": { "patch_l": "25.5mm" },
      "status": "Converged",
      "output_values": {
        "S11_at_center": -12.5,
        "PeakGain": 6.8,
        "BW_10dB": 42000000.0
      },
      "result_path": "variation_002/"
    }
  ]
}
```

### 5.9 Rust 类型定义

```rust
// ========================
// 结果索引（存储在 .emsp 工程文件中）
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolutionIndex {
    pub last_solve_time: Option<String>,
    pub setups: HashMap<String, SetupSolutionStatus>,
    pub optimetrics: HashMap<String, OptimetricsSolutionStatus>,
    pub stale_reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SetupSolutionStatus {
    pub status: SolveStatus,
    pub solved_at: Option<String>,
    pub converged_pass: Option<u32>,
    pub num_tetrahedra: Option<u64>,
    pub is_stale: bool,
    pub solved_variations: HashMap<String, VariationResult>,
    pub sweeps: HashMap<String, SweepResult>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SolveStatus {
    NotSolved,
    InProgress,
    Converged,
    NotConverged,
    Failed,
    Completed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariationResult {
    pub variables: HashMap<String, String>,
    pub result_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SweepResult {
    pub status: SolveStatus,
    pub num_frequency_points: Option<u32>,
    pub result_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimetricsSolutionStatus {
    pub status: SolveStatus,
    pub total_variations: u32,
    pub completed_variations: u32,
    pub result_path: String,
}

// ========================
// 结果数据（存储在 .emsp.results/ 目录中）
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolveProfile {
    pub setup: String,
    pub phases: Vec<ProfilePhase>,
    pub total_elapsed_sec: f64,
    pub total_peak_memory_mb: u64,
    pub cpu_cores_used: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProfilePhase {
    pub name: String,
    pub start_time: String,
    pub elapsed_sec: f64,
    pub peak_memory_mb: Option<u64>,
    pub sub_phases: Option<Vec<ProfilePhase>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FarFieldData {
    pub setup: String,
    pub frequency: String,
    pub coordinate_system: String,
    pub theta_deg: Vec<f64>,
    pub phi_deg: Vec<f64>,
    pub gain_total_dbi: Vec<Vec<f64>>,
    pub antenna_parameters: AntennaResultParameters,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AntennaResultParameters {
    pub peak_gain_dbi: f64,
    pub peak_directivity_dbi: f64,
    pub radiation_efficiency: f64,
    pub beamwidth_e_plane_deg: f64,
    pub beamwidth_h_plane_deg: f64,
    pub front_to_back_ratio_db: f64,
    pub radiated_power_w: f64,
    pub accepted_power_w: f64,
    pub incident_power_w: f64,
}
```

---

## 6. Rust 类型映射

以下展示核心 Rust struct 定义，均派生 `Serialize` / `Deserialize` 以便与 JSON 互转。

```rust
use serde::{Deserialize, Serialize};

// ========================
// 顶层工程
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EmProject {
    pub metadata: ProjectMetadata,
    pub variables: HashMap<String, Variable>,
    pub designs: Vec<Design>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectMetadata {
    pub version: String,
    pub application: String,
    pub created_at: String,
    pub modified_at: String,
    pub author: String,
    pub description: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Variable {
    pub value: String,
    pub description: String,
}

// ========================
// 设计
// ========================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SolutionType {
    // HFSS 全波求解
    DrivenModal,
    DrivenTerminal,
    Eigenmode,
    Transient,
    SBRPlus,
    // Q3D 准静态寄生提取
    Q3D_DCRL,    // DC 电阻 + 低频电感
    Q3D_ACRL,    // AC 电阻 + 电感（含趋肤/邻近效应）
    Q3D_C,       // 电容矩阵（静电求解）
    Q3D_CG,      // 电容 + 电导矩阵（含介质损耗）
}

impl SolutionType {
    /// 是否为 Q3D 准静态求解类型
    pub fn is_q3d(&self) -> bool {
        matches!(self, Self::Q3D_DCRL | Self::Q3D_ACRL | Self::Q3D_C | Self::Q3D_CG)
    }

    /// 是否为 HFSS 全波求解类型
    pub fn is_hfss(&self) -> bool {
        !self.is_q3d()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Design {
    pub id: String,
    pub name: String,
    pub solution_type: SolutionType,
    pub units: String,
    pub design_settings: DesignSettings,
    pub geometry: Geometry,
    pub materials: Vec<Material>,
    pub boundaries: Vec<Boundary>,
    pub excitations: Vec<Excitation>,
    pub nets: Vec<Net>,                  // Q3D: 网络定义与源汇分配
    pub mesh_operations: Vec<MeshOperation>,
    pub analysis_setups: Vec<AnalysisSetup>,
    pub radiation: RadiationSetup,       // HFSS: 远场/近场设置
    pub output_variables: Vec<OutputVariable>,
    pub field_overlays: Vec<FieldOverlay>,
    pub optimetrics: Vec<OptimetricsSetup>,
    pub reports: Vec<Report>,
}

// ========================
// 几何（历史记录式建模）
// 完整定义见 §3.3.4，此处为摘要
// ========================

// Geometry, GeometryOperation, OperationCommand, ObjectAttributes,
// GeometryObject, BoundingBox 等类型定义见 §3.3.4

// ========================
// 材料
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Material {
    pub name: String,
    pub category: MaterialCategory,
    pub properties: MaterialProperties,
    pub appearance: MaterialAppearance,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MaterialCategory {
    Conductor,
    Dielectric,
    Magnetic,
    Composite,
    Gas,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaterialProperties {
    pub permittivity: f64,
    pub permeability: f64,
    pub conductivity: f64,
    pub dielectric_loss_tangent: f64,
    pub magnetic_loss_tangent: f64,
    pub mass_density: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaterialAppearance {
    pub color: [u8; 3],
    pub transparency: f32,
}

// ========================
// 边界条件
// ========================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoundaryType {
    PerfectE,
    PerfectH,
    Radiation,
    PML,
    Impedance,
    FiniteConductivity,
    Symmetry,
    MasterSlave,
    // Q3D 专用边界
    ThinConductor,        // 薄导体边界（PCB 走线、薄涂层）
    InfiniteGroundPlane,  // 无限大地平面
    OpenBoundary,         // 开放边界（场延伸至无穷远）
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Boundary {
    pub name: String,
    #[serde(rename = "type")]
    pub boundary_type: BoundaryType,
    pub assignment: Assignment,
    pub properties: Option<serde_json::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Assignment {
    pub target_type: AssignmentTarget,
    pub targets: Vec<serde_json::Value>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AssignmentTarget {
    Object,
    Face,
    Edge,
}

// ========================
// 激励
// ========================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExcitationType {
    WavePort,
    LumpedPort,
    FloquetPort,
    IncidentWave,
    VoltageDrop,
    // Q3D 专用激励
    Source,   // 源端子（电流/电压注入点）
    Sink,     // 汇端子（电流回路点）
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Excitation {
    pub name: String,
    #[serde(rename = "type")]
    pub excitation_type: ExcitationType,
    pub assignment: Option<Assignment>,
    pub properties: serde_json::Value,
}

// ========================
// 网格操作
// ========================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MeshOperationType {
    LengthBased,
    SkinDepth,
    CurvatureBased,
    ModelResolution,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshOperation {
    pub name: String,
    #[serde(rename = "type")]
    pub mesh_type: MeshOperationType,
    pub assignment: Assignment,
    pub properties: serde_json::Value,
}

// ========================
// 分析设置
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisSetup {
    pub name: String,
    pub enabled: bool,
    pub solution_frequency: String,
    pub max_passes: u32,
    pub max_delta_s: f64,
    pub min_converged_passes: u32,
    pub order_basis: String,
    pub solver_type: String,
    pub frequency_sweeps: Vec<FrequencySweep>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SweepType {
    Discrete,
    Interpolating,
    Fast,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrequencySweep {
    pub name: String,
    pub sweep_type: SweepType,
    pub start: String,
    pub stop: String,
    pub step: String,
    pub save_fields: bool,
    pub save_rad_fields: bool,
}

// ========================
// 仿真状态与结果索引
// 结果数据类型的完整定义见 §5.9
// ========================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SimulationStatus {
    Idle,
    Validating,
    Meshing,
    Solving,
    PostProcessing,
    Finished,
    Failed,
}

// SolutionIndex, SetupSolutionStatus, SolveStatus, VariationResult,
// SweepResult, OptimetricsSolutionStatus, SolveProfile, ProfilePhase,
// ConvergenceData, FarFieldData, AntennaResultParameters
// 等类型定义见 §5.9

// ========================
// Design Settings
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DesignSettings {
    pub port_impedance_normalization: PortNormalization,
    pub deembedding: DeembeddingSettings,
    pub s_matrix_type: String,
    pub environment_temperature: String,
    pub model_validation: ValidationSettings,
    pub solver_options: SolverOptions,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PortNormalization {
    pub enabled: bool,
    pub reference_impedance: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeembeddingSettings {
    pub enabled: bool,
    pub default_distance: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationSettings {
    pub validate_before_solve: bool,
    pub check_intersections: bool,
    pub check_duplicate_boundaries: bool,
    pub check_port_on_boundary: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolverOptions {
    pub use_shell_elements: bool,
    pub curved_elements_order: String,
    pub allow_solver_fallback: bool,
}

// ========================
// Radiation（辐射设置）
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RadiationSetup {
    pub far_field_setups: Vec<FarFieldSetup>,
    pub near_field_setups: Vec<NearFieldSetup>,
    pub antenna_parameters: AntennaParameters,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FarFieldSetup {
    pub name: String,
    pub setup_type: String,        // "InfiniteSphere" | "InfinitePlane"
    pub coordinate_system: String,
    pub theta: AngleRange,
    pub phi: AngleRange,
    pub use_custom_radiation_surface: bool,
    pub radiation_surface: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AngleRange {
    pub start: String,
    pub stop: String,
    pub step: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum NearFieldSetup {
    Line {
        name: String,
        start_point: [f64; 3],
        end_point: [f64; 3],
        num_points: u32,
    },
    Rectangle {
        name: String,
        center: [f64; 3],
        width: f64,
        height: f64,
        axis: String,
        num_points_u: u32,
        num_points_v: u32,
    },
    Sphere {
        name: String,
        center: [f64; 3],
        radius: f64,
        num_theta: u32,
        num_phi: u32,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AntennaParameters {
    pub reference_impedance: String,
    pub calculate_antenna_params: bool,
}

// ========================
// Output Variables（输出变量）
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputVariable {
    pub name: String,
    pub expression: String,
    pub description: String,
}

// ========================
// Field Overlays（场叠加显示）
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldOverlay {
    pub name: String,
    pub quantity: FieldQuantity,
    pub component: String,        // "Mag" | "MagX" | "MagY" | "MagZ" | "Real" | "Imag" | "Vector"
    pub plot_type: FieldPlotType,
    pub assignment: serde_json::Value,
    pub solution: String,
    pub frequency: String,
    pub phase: String,
    pub scale: FieldScale,
    pub display: FieldDisplay,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FieldQuantity {
    E, H, Jvol, Jsurf, SAR, Poynting,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FieldPlotType {
    Surface, CutPlane, Volume, Line,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldScale {
    pub scale_type: String,       // "Linear" | "Log"
    pub min: Option<f64>,
    pub max: Option<f64>,
    pub unit: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldDisplay {
    pub plot_style: String,       // "Shaded" | "Arrow" | "Contour"
    pub show_arrows: bool,
    pub num_colors: u32,
    pub opacity: f32,
}

// ========================
// Optimetrics（参数化分析）
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum OptimetricsSetup {
    ParametricSweep {
        name: String,
        enabled: bool,
        setup: String,
        sweep_definitions: Vec<SweepDefinition>,
        constraints: Vec<serde_json::Value>,
        goals: Vec<serde_json::Value>,
    },
    Optimization {
        name: String,
        enabled: bool,
        setup: String,
        algorithm: String,
        max_iterations: u32,
        variables: Vec<OptimizationVariable>,
        goals: Vec<OptimizationGoal>,
        constraints: Vec<serde_json::Value>,
    },
    Sensitivity {
        name: String,
        enabled: bool,
        setup: String,
        variables: Vec<SensitivityVariable>,
        output: String,
        num_samples: u32,
    },
    Statistical {
        name: String,
        enabled: bool,
        setup: String,
        variables: Vec<SensitivityVariable>,
        num_trials: u32,
    },
    Tuning {
        name: String,
        enabled: bool,
        setup: String,
        variables: Vec<String>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SweepDefinition {
    pub variable: String,
    pub sweep_type: String,       // "LinearStep" | "LinearCount" | "LogScale" | "DiscreteList"
    // LinearStep: start, stop, step
    // DiscreteList: values
    pub params: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationVariable {
    pub variable: String,
    pub min: String,
    pub max: String,
    pub starting: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationGoal {
    pub name: String,
    pub expression: String,
    pub condition: String,        // "Minimize" | "Maximize" | "LessThan" | "GreaterThan" | "EqualTo"
    pub target: Option<f64>,
    pub frequency_range: Option<FrequencyRange>,
    pub weight: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrequencyRange {
    pub start: String,
    pub stop: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SensitivityVariable {
    pub variable: String,
    pub variation: String,
    pub distribution: String,     // "Uniform" | "Gaussian" | "LogNormal"
}

// ========================
// Reports（后处理报告）
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Report {
    pub name: String,
    pub category: ReportCategory,
    pub chart_type: ChartType,
    pub solution: String,
    pub domain: ReportDomain,
    pub traces: Vec<ReportTrace>,
    pub x_axis: Option<AxisConfig>,
    pub y_axis: Option<AxisConfig>,
    pub markers: Vec<ReportMarker>,
    pub limit_lines: Vec<LimitLine>,
    pub far_field_setup: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReportCategory {
    SParameter, FarField, NearField, Fields, Eigenmode, Emission,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ChartType {
    Rectangular, Polar, Smith, DataTable, Polar3D,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportDomain {
    pub domain_type: String,
    pub primary_sweep: String,
    pub fixed_values: Option<HashMap<String, String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportTrace {
    pub name: String,
    pub expression: String,
    pub style: Option<TraceStyle>,
    pub parametric_values: Option<HashMap<String, String>>,
    pub fixed_values: Option<HashMap<String, String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraceStyle {
    pub color: [u8; 3],
    pub line_width: u32,
    pub line_style: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AxisConfig {
    pub label: String,
    pub unit: String,
    pub min: Option<f64>,
    pub max: Option<f64>,
    pub auto_range: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportMarker {
    pub name: String,
    pub trace: String,
    pub x_value: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LimitLine {
    pub name: String,
    pub y_value: f64,
    pub style: Option<TraceStyle>,
}
```

---

## 7. 完整工程文件示例

以下是一个微带贴片天线仿真工程的完整 `.emsp` 文件示例：

```json
{
  "metadata": {
    "version": "1.0.0",
    "application": "EMStudio",
    "created_at": "2026-04-03T10:00:00Z",
    "modified_at": "2026-04-03T12:30:00Z",
    "author": "engineer@example.com",
    "description": "2.4GHz 微带贴片天线仿真"
  },
  "variables": {
    "$freq": { "value": "2.4GHz", "description": "中心频率", "unit_type": "Frequency" },
    "$eps_sub": { "value": "4.4", "description": "基板介电常数", "unit_type": "None" },
    "$sub_h": { "value": "1.6mm", "description": "FR4 基板厚度", "unit_type": "Length" }
  },
  "datasets": {
    "$ds_losstangent_vs_freq": {
      "description": "FR4 损耗正切随频率变化",
      "independent_variable": "Freq", "independent_unit": "GHz", "dependent_unit": "None",
      "data": [{"x":1.0,"y":0.018},{"x":5.0,"y":0.020},{"x":10.0,"y":0.025}],
      "interpolation": "PiecewiseLinear"
    }
  },
  "designs": [
    {
      "id": "design-001",
      "name": "Patch Antenna 2.4GHz",
      "solution_type": "DrivenModal",
      "units": "mm",
      "local_variables": {
        "patch_l": { "value": "28.5mm", "description": "贴片长度" },
        "patch_w": { "value": "37.0mm", "description": "贴片宽度" },
        "patch_x": { "expression": "(60mm - patch_l) / 2", "description": "贴片X居中" }
      },
      "definitions": {
        "materials": [
          {
            "name": "vacuum", "category": "Dielectric",
            "properties": {
              "permittivity": {"type":"constant","value":1.0}, "permeability": {"type":"constant","value":1.0},
              "conductivity": {"type":"constant","value":0.0}, "dielectric_loss_tangent": {"type":"constant","value":0.0}
            },
            "appearance": { "color": [200,200,255], "transparency": 0.95 }
          },
          {
            "name": "copper", "category": "Conductor",
            "properties": {
              "permittivity": {"type":"constant","value":1.0}, "permeability": {"type":"constant","value":0.999991},
              "conductivity": {"type":"constant","value":58000000.0}, "dielectric_loss_tangent": {"type":"constant","value":0.0}
            },
            "appearance": { "color": [255,180,50], "transparency": 0.0 }
          },
          {
            "name": "FR4_epoxy", "category": "Dielectric",
            "properties": {
              "permittivity": {"type":"expression","expression":"$eps_sub"},
              "permeability": {"type":"constant","value":1.0},
              "conductivity": {"type":"constant","value":0.0},
              "dielectric_loss_tangent": {"type":"dataset","dataset":"$ds_losstangent_vs_freq","independent_variable":"Freq"}
            },
            "appearance": { "color": [128,200,128], "transparency": 0.3 }
          }
        ],
        "coordinate_systems": [
          { "name": "Global", "type": "Cartesian", "origin": [0,0,0], "x_axis": [1,0,0], "y_axis": [0,1,0] }
        ],
        "named_selections": [
          { "name": "GND_Bottom", "type": "Face", "selection": [{"object":"GND_Plane","face":"ZMin"}], "description": "地平面底面" },
          { "name": "FeedPort", "type": "Face", "selection": [{"object":"FeedLine","face":"YMin"}], "description": "馈电端口面" },
          { "name": "RadiationFaces", "type": "Object", "selection": ["AirBox"], "description": "辐射边界" }
        ]
      },
      "geometry": {
        "coordinate_systems": [
          {
            "name": "Global",
            "type": "Cartesian",
            "origin": [0.0, 0.0, 0.0],
            "x_axis": [1.0, 0.0, 0.0],
            "y_axis": [0.0, 1.0, 0.0]
          }
        ],
        "operations": [
          {
            "step": 1, "command": "CreateBox", "result_object": "GND_Plane",
            "parameters": { "position": [0, 0, 0], "size": [60, 60, 0.035] },
            "attributes": { "material": "copper", "solve_inside": false, "color": [255,180,50], "transparency": 0.0, "group": "Antenna" }
          },
          {
            "step": 2, "command": "CreateBox", "result_object": "Substrate",
            "parameters": { "position": [0, 0, 0.035], "size": [60, 60, "$sub_h"] },
            "attributes": { "material": "FR4_epoxy", "solve_inside": true, "color": [128,200,128], "transparency": 0.4, "group": "Antenna" }
          },
          {
            "step": 3, "command": "CreateBox", "result_object": "Patch",
            "parameters": { "position": [15.75, 11.5, 1.635], "size": ["$patch_l", "$patch_w", 0.035] },
            "attributes": { "material": "copper", "solve_inside": false, "color": [255,180,50], "transparency": 0.0, "group": "Antenna" }
          },
          {
            "step": 4, "command": "CreateBox", "result_object": "FeedLine",
            "parameters": { "position": [28.5, 0, 1.635], "size": [3.0, 11.5, 0.035] },
            "attributes": { "material": "copper", "solve_inside": false, "color": [255,180,50], "transparency": 0.0, "group": "Feed" }
          },
          {
            "step": 5, "command": "CreateBox", "result_object": "AirBox",
            "parameters": { "position": [-30, -30, -30], "size": [120, 120, 62] },
            "attributes": { "material": "vacuum", "solve_inside": true, "color": [200,200,255], "transparency": 0.95, "group": "Environment" }
          }
        ],
        "objects": [
          { "id": 1, "name": "GND_Plane", "derived_from_step": 1, "material": "copper", "solve_inside": false, "color": [255,180,50], "transparency": 0.0, "group": "Antenna" },
          { "id": 2, "name": "Substrate", "derived_from_step": 2, "material": "FR4_epoxy", "solve_inside": true, "color": [128,200,128], "transparency": 0.4, "group": "Antenna" },
          { "id": 3, "name": "Patch", "derived_from_step": 3, "material": "copper", "solve_inside": false, "color": [255,180,50], "transparency": 0.0, "group": "Antenna" },
          { "id": 4, "name": "FeedLine", "derived_from_step": 4, "material": "copper", "solve_inside": false, "color": [255,180,50], "transparency": 0.0, "group": "Feed" },
          { "id": 5, "name": "AirBox", "derived_from_step": 5, "material": "vacuum", "solve_inside": true, "color": [200,200,255], "transparency": 0.95, "group": "Environment" }
        ]
      },
      "boundaries": [
        {
          "name": "Radiation1",
          "type": "Radiation",
          "assignment": { "target_type": "Object", "targets": ["@RadiationFaces"] }
        },
        {
          "name": "PEC_GND",
          "type": "PerfectE",
          "assignment": { "target_type": "Face", "targets": ["@GND_Bottom"] }
        }
      ],
      "excitations": [
        {
          "name": "Port1",
          "type": "LumpedPort",
          "assignment": { "target_type": "Face", "targets": ["@FeedPort"] },
          "properties": {
            "impedance": "50ohm",
            "integration_line": {
              "start": [29.5, 0.0, 0.035],
              "end": [29.5, 0.0, 1.635]
            }
          }
        }
      ],
      "mesh_operations": [
        {
          "name": "PatchRefine",
          "type": "LengthBased",
          "assignment": { "target_type": "Object", "targets": ["Patch"] },
          "properties": { "max_element_length": "2mm" }
        }
      ],
      "analysis_setups": [
        {
          "name": "Setup1",
          "enabled": true,
          "solution_frequency": "$freq",
          "max_passes": 15,
          "max_delta_s": 0.02,
          "min_converged_passes": 2,
          "order_basis": "Mixed",
          "solver_type": "Direct",
          "frequency_sweeps": [
            {
              "name": "Sweep1",
              "type": "Interpolating",
              "start": "1GHz",
              "stop": "4GHz",
              "step": "0.01GHz",
              "save_fields": true,
              "save_rad_fields": true
            }
          ]
        }
      ],
      "reports": [
        {
          "name": "S11_Return_Loss",
          "report_type": "SParameter",
          "expressions": ["dB(S(Port1,Port1))"],
          "domain": "Sweep1"
        },
        {
          "name": "Radiation_Pattern",
          "report_type": "FarField",
          "expressions": ["GainTotal"],
          "domain": "Setup1 : LastAdaptive"
        }
      ]
    }
  ]
}
```

### 7.2 Q3D 准静态寄生参数提取示例

以下是一个 PCB 差分走线寄生参数提取的 Q3D 设计片段（仅展示与 HFSS 差异的关键部分）：

```json
{
  "metadata": {
    "version": "1.0.0",
    "application": "EMStudio",
    "created_at": "2026-04-04T09:00:00Z",
    "modified_at": "2026-04-04T10:30:00Z",
    "author": "si_engineer@example.com",
    "description": "PCB 差分走线寄生参数提取"
  },
  "variables": {
    "$trace_w": { "value": "0.1mm", "description": "走线宽度", "unit_type": "Length" },
    "$trace_s": { "value": "0.15mm", "description": "差分对间距", "unit_type": "Length" },
    "$sub_h": { "value": "0.1mm", "description": "介质层厚度", "unit_type": "Length" },
    "$trace_l": { "value": "10mm", "description": "走线长度", "unit_type": "Length" }
  },
  "designs": [
    {
      "id": "design-q3d-001",
      "name": "Differential_Pair_Q3D",
      "solution_type": "Q3D_ACRL",
      "units": "mm",
      "geometry": {
        "operations": [
          {
            "id": "op1",
            "command": "CreateBox",
            "params": {
              "name": "Substrate",
              "position": [0, 0, 0],
              "size": ["$trace_l", "2mm", "$sub_h"],
              "material": "FR4_epoxy"
            }
          },
          {
            "id": "op2",
            "command": "CreateBox",
            "params": {
              "name": "Trace_P",
              "position": [0, "1 - $trace_s/2 - $trace_w", "$sub_h"],
              "size": ["$trace_l", "$trace_w", "0.035mm"],
              "material": "copper"
            }
          },
          {
            "id": "op3",
            "command": "CreateBox",
            "params": {
              "name": "Trace_N",
              "position": [0, "1 + $trace_s/2", "$sub_h"],
              "size": ["$trace_l", "$trace_w", "0.035mm"],
              "material": "copper"
            }
          },
          {
            "id": "op4",
            "command": "CreateBox",
            "params": {
              "name": "GND_Plane",
              "position": [0, 0, "-0.035mm"],
              "size": ["$trace_l", "2mm", "0.035mm"],
              "material": "copper"
            }
          }
        ]
      },
      "boundaries": [
        {
          "name": "OpenBC",
          "type": "OpenBoundary",
          "assignment": { "target_type": "Auto", "targets": [] }
        }
      ],
      "nets": [
        {
          "name": "DiffP",
          "objects": ["Trace_P"],
          "terminals": [
            {
              "name": "TP_src",
              "type": "Source",
              "assignment": {
                "target_type": "Face",
                "targets": [{"object": "Trace_P", "face": "XMin"}]
              }
            },
            {
              "name": "TP_sink",
              "type": "Sink",
              "assignment": {
                "target_type": "Face",
                "targets": [{"object": "Trace_P", "face": "XMax"}]
              }
            }
          ]
        },
        {
          "name": "DiffN",
          "objects": ["Trace_N"],
          "terminals": [
            {
              "name": "TN_src",
              "type": "Source",
              "assignment": {
                "target_type": "Face",
                "targets": [{"object": "Trace_N", "face": "XMin"}]
              }
            },
            {
              "name": "TN_sink",
              "type": "Sink",
              "assignment": {
                "target_type": "Face",
                "targets": [{"object": "Trace_N", "face": "XMax"}]
              }
            }
          ]
        },
        {
          "name": "GND",
          "objects": ["GND_Plane"],
          "is_ground_reference": true,
          "terminals": []
        }
      ],
      "analysis_setups": [
        {
          "name": "Q3D_Setup1",
          "enabled": true,
          "solution_type": "Q3D_ACRL",
          "adaptive_frequency": "5GHz",
          "max_passes": 10,
          "max_delta_energy": 0.02,
          "min_converged_passes": 2,
          "percent_refinement": 30,
          "solver_type": "Direct",
          "frequency_sweeps": [
            {
              "name": "WideSweep",
              "type": "Discrete",
              "start": "100MHz",
              "stop": "10GHz",
              "count": 40,
              "scale": "Logarithmic"
            }
          ],
          "dc_settings": {
            "compute_dc_resistance": true,
            "compute_dc_inductance": true
          }
        }
      ],
      "output_variables": [
        {
          "name": "Z_diff",
          "expression": "sqrt((L(DiffP:TP_src, DiffP:TP_src) + L(DiffN:TN_src, DiffN:TN_src) - 2*L(DiffP:TP_src, DiffN:TN_src)) / (C(DiffP:TP_src, DiffP:TP_src) + C(DiffN:TN_src, DiffN:TN_src) - 2*abs(C(DiffP:TP_src, DiffN:TN_src))))",
          "description": "差分特征阻抗估算"
        },
        {
          "name": "Coupling_k",
          "expression": "L(DiffP:TP_src, DiffN:TN_src) / sqrt(L(DiffP:TP_src, DiffP:TP_src) * L(DiffN:TN_src, DiffN:TN_src))",
          "description": "电感耦合系数"
        }
      ],
      "field_overlays": [
        {
          "name": "J_Current_DP",
          "quantity": "Jvol",
          "component": "Mag",
          "plot_type": "Surface",
          "assignment": {
            "type": "Object",
            "targets": ["Trace_P", "Trace_N"]
          },
          "solution": "Q3D_Setup1 : LastAdaptive",
          "frequency": "5GHz",
          "phase": "0deg",
          "scale": {
            "type": "Log",
            "min": null,
            "max": null,
            "unit": "A/m2"
          }
        }
      ],
      "reports": [
        {
          "name": "RLCG_vs_Freq",
          "category": "RLCGMatrix",
          "chart_type": "Rectangular",
          "solution": "Q3D_Setup1 : WideSweep",
          "domain": { "type": "Frequency", "primary_sweep": "Freq" },
          "traces": [
            { "name": "R_DiffP", "expression": "R(DiffP:TP_src, DiffP:TP_src)" },
            { "name": "L_DiffP", "expression": "L(DiffP:TP_src, DiffP:TP_src)" },
            { "name": "L_mutual", "expression": "L(DiffP:TP_src, DiffN:TN_src)" }
          ]
        },
        {
          "name": "L_Matrix_5GHz",
          "category": "RLCGMatrix",
          "chart_type": "MatrixTable",
          "solution": "Q3D_Setup1 : WideSweep",
          "domain": {
            "type": "Frequency",
            "primary_sweep": "Freq",
            "fixed_values": { "Freq": "5GHz" }
          },
          "matrix_type": "L"
        }
      ],
      "optimetrics": [
        {
          "name": "SpacingSweep",
          "type": "ParametricSweep",
          "enabled": true,
          "setup": "Q3D_Setup1",
          "sweep_definitions": [
            {
              "variable": "$trace_s",
              "type": "LinearStep",
              "start": "0.05mm",
              "stop": "0.5mm",
              "step": "0.05mm"
            }
          ],
          "goals": []
        }
      ]
    }
  ]
}
```

---

## 8. 与现有代码的映射关系

当前 `crates/domain/src/lib.rs` 中的类型与新方案的对应关系：

| 现有类型 | 新方案对应 | 改动说明 |
|---------|-----------|---------|
| `Project` | `EmProject` + `Design` | 拆分为工程(多设计容器)和设计两层 |
| `EmModel` | `Definitions` + `Geometry` | 模型拆分为**定义层**（材料/坐标系/命名选择）+ 几何层 |
| `GeometryObject` | `GeometryOperation` + `GeometryObject` | 从单一对象升级为**操作历史 + 对象快照**双层结构（参考 HFSS History Tree） |
| `Material` | `Material`（增强版） | 属性值支持 constant/expression/dataset 三种来源 |
| `SimulationStatus` | `SimulationStatus`（扩展版） | 增加 Validating/Meshing/PostProcessing 阶段 |
| `SolveResult` | `SolutionIndex` + 结果文件体系 | 工程文件中存索引（状态/路径/过期标记），结果数据存 `.emsp.results/` |
| _(新增)_ | `Variable` + `PropertyValue` | 参数化变量系统，支持表达式引用链 |
| _(新增)_ | `DatasetDefinition` | 命名查找表（频率/温度依赖曲线等） |
| _(新增)_ | `NamedSelection` | 面/边/对象命名选择，稳定引用标识 |
| _(新增)_ | `DependencyGraph` | 定义-引用依赖关系图，支持级联验证 |
| _(新增)_ | `OperationCommand` | 30+ 种建模操作命令枚举 |
| _(新增)_ | `Boundary` | 参考 HFSS 边界条件体系 |
| _(新增)_ | `Excitation` | 参考 HFSS 端口/激励体系 |
| _(新增)_ | `MeshOperation` | 参考 HFSS 网格控制 |
| _(新增)_ | `AnalysisSetup` + `FrequencySweep` | 参考 HFSS 分析设置与扫频 |
| _(新增)_ | `DesignSettings` | 设计全局设置（端口归一化、验证、求解器选项） |
| _(新增)_ | `RadiationSetup` + `FarFieldSetup` + `NearFieldSetup` | 远场球面/近场采样定义，方向图前提 |
| _(新增)_ | `OutputVariable` | 从仿真结果派生的命名表达式 |
| _(新增)_ | `FieldOverlay` + `FieldQuantity` | 场可视化方案定义（E/H/J/SAR 叠加显示） |
| _(新增)_ | `OptimetricsSetup`（5 种变体） | 参数扫描/优化/灵敏度/统计/调参 |
| _(新增)_ | `Report`（增强版）+ `ReportTrace` + `ChartType` | 完整报告定义（图表类型、坐标轴、标记线、多 trace） |

---

## 9. 文件操作 API 设计

### 9.1 序列化/反序列化

```rust
impl EmProject {
    /// 从 .emsp 文件加载工程
    pub fn load(path: &Path) -> Result<Self, ProjectError> {
        let content = std::fs::read_to_string(path)?;
        let project: EmProject = serde_json::from_str(&content)?;
        project.validate()?;
        Ok(project)
    }

    /// 保存工程到 .emsp 文件
    pub fn save(&self, path: &Path) -> Result<(), ProjectError> {
        let content = serde_json::to_string_pretty(self)?;
        std::fs::write(path, content)?;
        Ok(())
    }

    /// 验证工程完整性
    pub fn validate(&self) -> Result<(), ValidationError> {
        for design in &self.designs {
            design.validate()?;
        }
        Ok(())
    }
}

impl Design {
    /// 验证设计完整性（材料引用、端口面合法性等）
    pub fn validate(&self) -> Result<(), ValidationError> {
        // 检查几何对象引用的材料是否存在于材料库中
        // 检查边界条件和激励引用的对象/面是否存在
        // 检查分析设置的参数合理性
        todo!()
    }
}
```

### 9.2 锁文件管理

```rust
pub struct ProjectLock {
    path: PathBuf,
}

impl ProjectLock {
    pub fn acquire(project_path: &Path) -> Result<Self, LockError> {
        let lock_path = project_path.with_extension("emsp.lock");
        if lock_path.exists() {
            return Err(LockError::AlreadyLocked);
        }
        let info = LockInfo {
            pid: std::process::id(),
            hostname: hostname::get()?.to_string_lossy().to_string(),
            locked_at: chrono::Utc::now().to_rfc3339(),
        };
        std::fs::write(&lock_path, serde_json::to_string(&info)?)?;
        Ok(Self { path: lock_path })
    }
}

impl Drop for ProjectLock {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}
```

### 9.3 自动保存

```rust
impl EmProject {
    pub fn auto_save(&self, project_path: &Path) -> Result<(), ProjectError> {
        let auto_path = project_path.with_extension("emsp.auto");
        self.save(&auto_path)
    }

    pub fn recover(project_path: &Path) -> Option<Self> {
        let auto_path = project_path.with_extension("emsp.auto");
        if auto_path.exists() {
            Self::load(&auto_path).ok()
        } else {
            None
        }
    }
}
```

---

## 10. 版本迁移策略

工程文件通过 `metadata.version` 字段标识格式版本，采用语义化版本号：

| 版本 | 变更类型 | 处理方式 |
|------|---------|---------|
| `1.0.x` → `1.0.y` | Patch: 修复/补充字段默认值 | 自动兼容，无需迁移 |
| `1.x` → `1.y` | Minor: 新增可选字段 | 向后兼容，旧文件可直接打开 |
| `N.x` → `M.x` | Major: 结构性变更 | 需要迁移函数，提示用户确认 |

```rust
pub fn migrate(project: serde_json::Value) -> Result<EmProject, MigrationError> {
    let version = project["metadata"]["version"].as_str().unwrap_or("0.0.0");
    match version {
        v if v.starts_with("1.") => serde_json::from_value(project).map_err(Into::into),
        "0.1.0" => migrate_v0_1_to_v1_0(project),
        _ => Err(MigrationError::UnsupportedVersion(version.to_string())),
    }
}
```

---

## 11. 参考资料

### HFSS 相关

- [Ansys Electronics Desktop Files](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/Variables/ANSYSElectronicsDesktopFiles.htm) - AEDT 文件体系
- [An Introduction to HFSS](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/PDFs/An%20Introduction%20to%20HFSS.pdf) - HFSS 入门与工程结构
- [Modeling Practice in HFSS](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Subsystems/An%20Introduction%20to%20HFSS/Content/ModelingPracticeinHFSS.htm) - HFSS 建模实践
- [Design Settings for HFSS](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/Modeler/DesignSettingsforHFSS.htm) - HFSS 设计设置
- [Assigning Boundaries in HFSS](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/HFSS/AssigningBoundariesinHFSSandHFSSIE.htm) - HFSS 边界条件
- [PyAEDT Hfss API](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.hfss.Hfss.html) - PyAEDT HFSS 类参考
- [PyAEDT Project Configuration](https://aedt.docs.pyansys.com/version/stable/User_guide/pyaedt_file_data/project.html) - PyAEDT 配置文件格式
- [PyAEDT Materials](https://aedt.docs.pyansys.com/version/stable/User_guide/pyaedt_file_data/materials.html) - PyAEDT 材料文件格式
- [Working with Design Variables](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/Variables/DesignVariables.htm) - HFSS 设计变量与表达式
- [Working with Project Variables](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/Variables/ProjectVariables.htm) - HFSS 工程变量（$前缀）
- [Defining Material Properties as Expressions](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/HFSS3DLayout/Content/Materials/DefiningMaterialPropertiesasExpressions.htm) - 材料属性表达式定义
- [Defining Frequency-Dependent Material Properties](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/Materials/DefiningFrequencyDependentMaterialProperties.htm) - 频率依赖材料
- [Using Dataset Expressions](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS3DLayout/Content/Variables/UsingDatasetExpressions.htm) - 数据集表达式
- [Spatially Dependent Boundaries](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/HFSS/SpatiallyDependentBoundaries.htm) - 空间依赖边界条件
- [HFSS History Tree](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Subsystems/An%20Introduction%20to%20HFSS/Content/HistoryTree.htm) - 历史记录树
- [Cleaning Up Model History](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/Modeler/CleaningUpHistory.htm) - 模型历史清理

### Q3D Extractor 相关

- [Q3D Extractor Getting Started](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/GettingStarted/Q3DExtractorGettingStartedGuides.htm) - Q3D 入门指南
- [Q3D Extractor Help (PDF)](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/PDFs/Q3D%20Extractor.pdf) - Q3D 完整帮助文档
- [Q3D Scripting Guide](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/PDFs/Q3D%20ExtractorScriptingGuide.pdf) - Q3D 脚本编程指南
- [Assigning Thin Conductor Boundaries](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/AssigningThinConductorBoundaries.htm) - Q3D 薄导体边界
- [Viewing Matrix Data in Q3D](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/ViewingMatrixDatainQ3D.htm) - Q3D 矩阵数据查看
- [Exporting Equivalent Circuit Data](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/ExportingQ3DExtractorEquivalentCircuitData.htm) - Q3D 等效电路导出
- [Frequency Sweeps in Q3D](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/FrequencySweepsinQ3DExtractor.htm) - Q3D 频率扫描
- [Plotting Field Overlays in Q3D](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Q3DExtractor/Content/ReportsandPostProc/PlottingFieldOverlaysinQ3D.htm) - Q3D 场叠加显示
- [Exporting S-Parameters from Q3D](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Q3DExtractor/Content/Q3D/ExportingSParameterData.htm) - Q3D S 参数导出
- [PyAEDT Q3D Class Reference](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.q3d.Q3d.html) - PyAEDT Q3D API
- [PyAEDT Q3D Setup Templates](https://aedt.docs.pyansys.com/version/stable/API/SetupTemplatesQ3D.html) - Q3D 分析设置模板
- [Ansys Q3D Extractor Product Page](https://www.ansys.com/products/electronics/ansys-q3d-extractor) - Q3D 产品主页
