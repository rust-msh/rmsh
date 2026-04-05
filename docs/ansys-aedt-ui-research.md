# Ansys Electronics Desktop (AEDT) UI 设计调研

## 1. 整体界面布局

AEDT 采用标准工程 CAD 应用布局，从上到下依次为：

```
┌──────────────────────────────────────────────────────────┐
│  Title Bar                                               │
├──────────────────────────────────────────────────────────┤
│  Menu Bar (File, Edit, View, Project, Draw, Modeler ...) │
├──────────────────────────────────────────────────────────┤
│  Quick Access Toolbar (QAT)                              │
├──────────────────────────────────────────────────────────┤
│  Ribbon (Desktop | View | Simulation | Automation | ...) │
├────────────┬─────────────────────────────┬───────────────┤
│  Project   │                             │               │
│  Manager   │     Central Workspace       │  (Optional    │
│  (Tree)    │     (3D Modeler /           │   Side Panel) │
│            │      Layout /               │               │
│  ──────    │      Schematic Editor)      │               │
│  Properties│                             │               │
│  Window    │                             │               │
├────────────┴─────────────────────────────┴───────────────┤
│  Message Manager / Progress Window (Bottom Dock)         │
├──────────────────────────────────────────────────────────┤
│  Status Bar (坐标, 单位, 消息/进度开关)                    │
└──────────────────────────────────────────────────────────┘
```

---

## 2. Ribbon 工具栏

### 2.1 Ribbon 整体尺寸

```
┌─────────────────────────────────────────────────────────────────────┐
│ Tab Strip (25-30px)                                                 │
│  ┌─────────┐ ┌──────┐ ┌────────────┐ ┌────────────┐                │
│  │ Desktop │ │ View │ │ Simulation │ │ Automation │ │ Draw │ ...    │
│  └─────────┘ └──────┘ └────────────┘ └────────────┘                │
├─────────────────────────────────────────────────────────────────────┤
│ Command Area (~78px)                                                │
│ ┌───────────────────────────────┬─┬──────────────────────────────┐  │
│ │          Group A              │ │          Group B             │  │
│ │ SubGrp1       │ SubGrp2      │ │ SubGrp1        │ SubGrp2    │  │
│ │ [LG][LG][LG]  │ [sm]        │ │ [LG][LG]       │ [sm][sm]   │  │
│ │                │ [sm]        │ │                 │ [sm][sm]   │  │
│ │                │ [sm]        │ │                 │ [sm]       │  │
│ │   ← subgroup separator →    │ │   ← subgroup separator →    │  │
│ ├───────────────────────────────┤ ├──────────────────────────────┤  │
│ │         Group A               │ │         Group B              │  │
│ └───────────────────────────────┴─┴──────────────────────────────┘  │
│  Group Label Row (~16px)       ↑ group separator (1px + 6px pad)   │
└─────────────────────────────────────────────────────────────────────┘
总高度: ~115-130px (Tab Strip 28px + Command Area 78px + Group Label 16px)
```

### 2.2 固定 Tab（始终可见）

| Tab | 功能 |
|-----|------|
| **Desktop** | 项目级操作：新建、打开、保存、导入/导出、通用选项 |
| **View** | 窗口布局配置；内容随当前编辑器变化（3D Modeler / Report / Layout） |
| **Simulation** | 验证、分析设置、求解、任务管理和监控 |
| **Automation** | 脚本、ACT 扩展、PyAEDT、自定义工具包 |

### 2.3 上下文相关 Tab（随设计类型动态出现）

当插入设计后，根据求解器类型出现额外 Tab：

- **HFSS / Maxwell / Icepak / Q3D** 设计 → 增加 **Draw**, **Model**, **Results**
- **HFSS 3D Layout** 设计 → 增加 **Layout**, **Results**
- **Circuit / Twin Builder** 设计 → 增加 **Schematic**, **Results**

这些 Tab 仅在相关设计编辑器活跃时才显示。上下文 Tab 在 Tab Strip 上方有一条彩色强调条（accent bar），颜色与 Tab 组类别关联，用于视觉区分常驻 Tab 和上下文 Tab。

### 2.4 Tab Strip 交互行为

| 行为 | 视觉表现 |
|------|----------|
| **Normal Tab** | 文字颜色 `#444`~`#5A5A5A`，背景透明 |
| **Hover Tab** | 背景微亮（浅灰色填充），文字颜色不变 |
| **Active Tab** | 文字变为 Ansys 蓝 `#0055A5`~`#0070C0`；Tab 底边与 Ribbon body 无边框分隔（视觉上 Tab "连通"到 Ribbon 内容区）；Tab 背景与 Ribbon body 同色 |
| **上下文 Tab 组** | Tab 上方有 ~3px 彩色条（如蓝色/绿色/橙色），标识所属功能组 |

---

### 2.5 按钮类型与视觉规格

#### 2.5.1 Large Button（大按钮）

```
┌─────────────┐
│             │
│   [32x32]   │   ← 32x32px 图标居中
│    icon     │
│             │
│  Label Ln1  │   ← 文字居中，最多 2 行
│  Label Ln2  │   ← 每行约 6-8 个英文字符
└─────────────┘
整体尺寸: 约 44-60px 宽 × 66px 高（不含 group label）
```

- **用途**: 高频/核心操作，如 New, Open, Save, Analyze, Box, Cylinder 等
- **文字**: 位于图标下方，居中对齐，9-10pt 字号，最多 2 行自动换行
- **间距**: 图标上方 padding ~4px，图标与文字间 ~2px，文字下方 ~2px

#### 2.5.2 Small Button（小按钮）

```
┌──────────────────────┐
│ [16x16] Button Label │   ← 16x16px 图标 + 文字在右侧
└──────────────────────┘
整体尺寸: 高 ~22px（3 个小按钮堆叠 = 66px ≈ 1 个大按钮高度）
```

- **用途**: 次要操作或空间紧凑时，如 Group, Ungroup, Assign Color 等
- **文字**: 位于图标右侧，左对齐，单行
- **排列**: 通常 3 个一列竖直堆叠，与 1 个大按钮等高

#### 2.5.3 Split Button（分裂按钮）

```
Large Split:                    Small Split:
┌─────────────┐                ┌──────────────────┬──┐
│   [32x32]   │                │ [16x16] Label    │▼ │
│    icon     │                └──────────────────┴──┘
├─────────────┤ ← 视觉分隔线
│ Label    ▼  │ ← 下半部分点击打开下拉菜单
└─────────────┘
```

- **上半区/主区**: 点击执行默认操作
- **下半区/箭头区**: 点击打开下拉菜单，箭头区宽度 ~12-14px
- **视觉分隔**: 上下两区之间有 1px 细线或在 hover 时才显示分隔
- **用途**: 有默认操作 + 可选变体，如 Paste (粘贴 / 选择性粘贴)

#### 2.5.4 Dropdown Button（纯下拉按钮）

```
Large Dropdown:                Small Dropdown:
┌─────────────┐               ┌──────────────────┬──┐
│   [32x32]   │               │ [16x16] Label    │▼ │
│    icon     │               └──────────────────┴──┘
│ Label    ▼  │ ← 整个按钮都是下拉触发器
└─────────────┘
```

- **与 Split Button 区别**: 没有分隔线，整个按钮点击都打开菜单，无默认直接操作
- **箭头**: ▼ 三角形紧跟标签文字右侧

#### 2.5.5 Toggle Button（切换按钮）

- **外观**: 与 Large/Small Button 相同的尺寸和图标规格
- **区别**: 有 checked/unchecked 两种状态
- **Checked 状态视觉**: 持续显示高亮背景填充 + 边框（类似 hover 但更深的蓝色），即使鼠标离开也保持
- **用途**: 如 Grid 显示开关、Snap 开关、选择模式切换（Object/Face/Edge/Vertex）

#### 2.5.6 Gallery（图库控件）

```
┌─────────────────────────┐
│ [item][item][item][▲]   │   ← Ribbon 内嵌显示 1 行
│ [item][item][item][▼]   │   ← 上/下箭头滚动
│                   [▼▼]  │   ← 展开按钮打开完整 gallery popup
└─────────────────────────┘

展开后的 Gallery Popup:
┌─────────────────────────────┐
│  Category Header            │
│ ┌─────┐┌─────┐┌─────┐      │
│ │ opt ││ opt ││ opt │ ...  │  ← 网格排列的可视化选项
│ └─────┘└─────┘└─────┘      │
│  Category Header 2          │
│ ┌─────┐┌─────┐             │
│ │ opt ││ opt │              │
│ └─────┘└─────┘             │
├─────────────────────────────┤
│  More Options...            │  ← 底部链接打开完整对话框
└─────────────────────────────┘
```

- **用途**: Automation Tab 的 Toolkit gallery、材料快速选择等

#### 2.5.7 内嵌控件

| 控件 | 说明 | 示例 |
|------|------|------|
| **Combo Box** | 内嵌下拉选择框，约 100-150px 宽，22px 高 | 绘图平面选择器 (XY/YZ/XZ) |
| **Spinner / 数值框** | 带上下箭头的数值输入，约 60-80px 宽 | 坐标入口、网格间距 |
| **Label** | 静态文字标签，常配合 Combo Box | "Plane:" "Units:" |
| **Checkbox** | 内嵌复选框 + 文字 | View Tab 中的显示选项 |

---

### 2.6 按钮交互状态（适用于所有按钮类型）

| 状态 | 视觉变化 | 近似样式 |
|------|----------|----------|
| **Normal** | 平面/透明背景，无边框 | `background: transparent; border: none` |
| **Hover** | 浅蓝色填充 + 1px 边框出现 | `background: #DCEEFB; border: 1px solid #B8D6F0` |
| **Pressed** | 更深的蓝色填充 + 边框 | `background: #B8D6F0; border: 1px solid #7EB4EA` |
| **Disabled** | 图标变灰度 + 降低透明度，文字变浅灰 | `opacity: 0.4; filter: grayscale(1); color: #A0A0A0` |
| **Checked/Active** | 持续的强调色填充 + 边框（比 hover 更饱和） | `background: #C4DFF0; border: 1px solid #7EB4EA` |
| **Keyboard Focus** | 虚线框或细线聚焦框 | `outline: 1px dotted #333` |

**状态转换说明**:
- Normal → Hover: 鼠标进入按钮区域，~0ms 延迟即时响应
- Hover → Pressed: 鼠标按下
- Pressed → Normal/Checked: 鼠标释放后，普通按钮回到 Normal，Toggle 按钮切换 Checked 状态
- Disabled 状态时不响应任何鼠标交互，无 hover 效果
- Split Button 的上半区和下半区有独立的 hover 状态

---

### 2.7 Group 与 SubGroup 布局模式

#### 层级结构

Ribbon 的按钮组织采用三级层级：

```
Tab
 └── Group (底部显示 group label，group 之间有较宽间距 + 分隔线)
      └── SubGroup (同一 group 内的按钮子组，子组之间有竖直分隔线)
           └── Item (具体的按钮/控件)
```

**Group vs SubGroup 的区别**：
- **Group**: 底部有居中文字标签（如 "Project"、"Primitives"），group 之间有间距 + 分隔线，是用户可感知的功能大类
- **SubGroup**: 同一 Group 内部的按钮分组，仅用竖线隔开，无独立标签。用于将同一功能类别内的按钮按操作性质进一步分组

**示例 — Desktop Tab 的 "Project" Group**：

```
┌───────────────────────────────────────┐
│ [New▲] [Open▲] [Save▲] │ [Save As▼] │   ← SubGroup 1 │ SubGroup 2
│                         │ [Close▼]   │      竖线分隔子组
├───────────────────────────────────────┤
│              Project                  │   ← Group label
└───────────────────────────────────────┘
```

**示例 — Draw Tab 的 "Primitives" Group**：

```
┌──────────────────────────────────────────────────────┐
│ [Box▲] [Cyl▲] [Sph▲] [Poly▲] │ [Cone▼]  [Torus▼]  │
│                                │ [Rect▼]  [Ellip▼]  │
│                                │ [Circ▼]  [Poly▼]   │
│                                │ [Arc▼]   [Spline▼] │
├──────────────────────────────────────────────────────┤
│                   Primitives                         │
└──────────────────────────────────────────────────────┘
```

#### SubGroup 内部排列 Pattern

| Pattern | 布局 | 适用场景 |
|---------|------|----------|
| **A: 单个大按钮** | `[LG]` | 独立核心操作 |
| **B: 多个大按钮并排** | `[LG][LG][LG]` | 同级高频操作（如 New/Open/Save） |
| **C: 大按钮 + 小按钮栈** | `[LG] [sm]` | 1 个主操作 + 多个辅助 |
|  |  `     [sm]` | |
|  |  `     [sm]` | |
| **D: 三行小按钮** | `[sm][sm][sm]` | 多个同级次要操作 |
|  | `[sm][sm][sm]` | |
|  | `[sm][sm]    ` | |
| **E: 大按钮 + 内嵌控件** | `[LG] [combo ▼]` | 主操作 + 选择器 |
|  |  `     [combo ▼]` | （如 Draw + Plane 选择） |
|  |  `     [sm]     ` | |

#### SubGroup 分隔线

```
... │ ...
    │       ← 1px 竖直线，颜色 #A0A0A0
... │ ...       上下各内缩 4px，左右各 3px padding
```

- SubGroup 之间的分隔线在 Group 的**按钮区域内**绘制，不延伸到底部的 Group label 区域

#### Group 分隔线

Group 之间有更明显的视觉分隔：左右各 6px padding + 1px 竖线贯穿整个 Group 高度（含 label 区域）

#### Group 标签

- **位置**: Group 底部居中
- **字号**: 8-9pt
- **颜色**: `#444444`~`#666666`
- **背景**: 与 Ribbon body 相同（无额外背景色）
- **可点击性**: 部分 Group 标签右侧有小箭头（Dialog Launcher ↗），点击打开完整设置对话框

```
┌────────────────────────────┐
│     [按钮区域内容]          │
├────────────────────────────┤
│    Group Name          [↗] │  ← Dialog Launcher 箭头（可选）
└────────────────────────────┘
```

---

### 2.8 下拉菜单规格

#### 菜单外观

```
┌─────────────────────────────────┐  ← 1px border #A0A0A0
│  Menu Item                      │     border-radius: 2-4px
│──────────────────────────────── │     box-shadow: 0 2px 8px rgba(0,0,0,0.15)
│  [16x16] Item Label    Ctrl+X  │     背景: #FFFFFF
│  [16x16] Item Label    Ctrl+C  │
│  [16x16] Item Label    Ctrl+V  │
│─────────────────────────────── │  ← 分隔线: 1px #D0D0D0, 左右各留 icon 列宽度
│  Section Header                 │  ← 不可点击，粗体或灰色文字
│  [16x16] Item Label             │
│  [16x16] Item Label         ▶  │  ← 子菜单箭头
└─────────────────────────────────┘
```

#### 菜单项详细规格

| 元素 | 规格 |
|------|------|
| **Menu Item 高度** | ~32px（标准），~24px（紧凑模式） |
| **图标区域** | 左侧 ~36px 列宽，16x16 图标居中 |
| **文字** | 图标右侧，左对齐，9-10pt |
| **快捷键文字** | 右对齐，灰色 `#888888` |
| **子菜单箭头** | 右侧 ▶ 符号，hover 后 200-400ms 延迟打开子菜单 |
| **分隔线** | 1px 水平线 `#D0D0D0`，上下各 2-4px padding |
| **Section Header** | 不可点击，粗体或颜色 `#666666`，无图标 |

#### 菜单项交互状态

| 状态 | 视觉 |
|------|------|
| **Normal** | 白色背景 |
| **Hover** | 浅蓝高亮 `#DCEEFB`，可能有 1px 边框 |
| **Disabled** | 文字和图标变灰 `#A0A0A0`，不可点击 |
| **Checked** | 图标位置显示 ✓ 勾选标记 |

---

### 2.9 Tooltip / Screentip（增强提示）

AEDT 使用 Microsoft Office 风格的增强 Screentip：

```
┌──────────────────────────────┐
│  Command Name        (Ctrl+N)│  ← 粗体标题 + 快捷键
├──────────────────────────────┤
│  Detailed description of     │  ← 普通字体描述文字
│  what this command does.     │
│  Can span multiple lines.    │
│                              │
│  [可选: 命令预览图片]          │  ← 部分命令有预览图
└──────────────────────────────┘
```

| 参数 | 值 |
|------|-----|
| **首次显示延迟** | ~500ms |
| **显示持续时间** | ~5000ms |
| **再次显示延迟** | ~100ms（鼠标移到相邻按钮时快速显示） |
| **位置** | 按钮下方，不遮挡按钮本身 |
| **背景色** | 浅黄 `#FFFFEF` 或白色 `#FFFFFF` |
| **边框** | 1px `#A0A0A0` |
| **阴影** | 轻微投影 |

---

### 2.10 自适应折叠（Responsive Collapse）

窗口变窄时，Ribbon Group 按优先级依次经历 4 个折叠阶段：

```
Stage 1 - Large (完全展开):
┌────────────────────────────────────┐
│ [LG-icon]  [LG-icon]  [sm] [sm]   │
│  Label      Label     [sm] [sm]   │
│                       [sm] [sm]   │
│         Group Name                 │
└────────────────────────────────────┘

Stage 2 - Medium (大按钮缩为小按钮):
┌──────────────────────────────┐
│ [sm] Label  [sm] Label       │
│ [sm] Label  [sm] Label       │
│ [sm] Label  [sm] Label       │
│       Group Name             │
└──────────────────────────────┘

Stage 3 - Small (隐藏文字标签，仅图标):
┌──────────────────┐
│ [ic][ic][ic][ic]  │
│ [ic][ic][ic][ic]  │
│ [ic][ic]          │
│   Group Name      │
└──────────────────┘

Stage 4 - Popup (整个 Group 折叠为 1 个按钮):
┌──────┐
│[icon]│
│  ▼   │
│Group │
│ Name │
└──────┘  ← 点击展开为浮动面板，显示完整 Group 内容
```

**折叠优先级规则**:
- 每个 Group 有 `ScalingPolicy` 定义折叠顺序
- 宽度最不重要的 Group 先折叠（通常含较少/次要按钮的 Group）
- 同一 Group 内从右到左折叠
- 最重要的 Group（如 Draw Tab 的 Primitives）最后折叠

---

### 2.11 Quick Access Toolbar (QAT)

```
位于 Ribbon 上方（可选下方）:
┌──────────────────────────────────────────────┐
│ [💾][↩][↪][▶] .... [▼]                      │  ← 16x16 icon-only buttons
│                       └── Customize dropdown │      hit target: 22x22px
└──────────────────────────────────────────────┘
```

| 元素 | 规格 |
|------|------|
| **按钮图标** | 16x16px，仅图标无文字 |
| **按钮点击区域** | 22x22px |
| **分隔线** | 1px 竖直线，上下各 2px padding |
| **自定义下拉箭头** | QAT 最右端，▼ 箭头，点击显示常用命令列表 + "Show Below/Above the Ribbon" 切换 |
| **右键菜单** | 任意 Ribbon 按钮右键 → "Add to Quick Access Toolbar" |

---

### 2.12 各 Tab 按钮详细清单

> 表格中 `│` 表示同一 Group 内的 SubGroup 分界

#### Desktop Tab

| Group | SubGroup | 按钮 | 类型 | 大小 |
|-------|----------|------|------|------|
| **Project** | 1 | New | Large Button | LG |
|  | 1 | Open | Large Button | LG |
|  | 1 | Save | Large Button | LG |
|  | │ | — | *SubGroup 分隔线* | — |
|  | 2 | Save As | Small Button | SM |
|  | 2 | Close | Small Button | SM |
| **Import/Export** | 1 | Import | Large Split (下拉选择格式) | LG |
|  | 1 | Export | Large Split | LG |

#### Draw Tab (HFSS/Maxwell 3D)

| Group | SubGroup | 按钮 | 类型 | 大小 |
|-------|----------|------|------|------|
| **Primitives** | 1 | Box | Large Split | LG |
|  | 1 | Cylinder | Large Split | LG |
|  | 1 | Sphere | Large Button | LG |
|  | 1 | Polyline | Large Button | LG |
|  | │ | — | *SubGroup 分隔线* | — |
|  | 2 | Cone | Small Button | SM |
|  | 2 | Torus | Small Button | SM |
|  | 2 | Rectangle | Small Button | SM |
|  | 2 | Ellipse | Small Button | SM |
|  | 2 | Circle | Small Button | SM |
|  | 2 | Polygon | Small Button | SM |
|  | 2 | Arc | Small Button | SM |
|  | 2 | Spline | Small Button | SM |
| **Plane** | 1 | [Plane Combo Box] | Combo Box (XY/YZ/XZ) | — |
| **Units** | 1 | Units | Large Dropdown (mm/cm/m/mil/in) | LG |
| **Material** | 1 | Assign Material | Large Button | LG |

#### Model Tab

| Group | SubGroup | 按钮 | 类型 | 大小 |
|-------|----------|------|------|------|
| **Boolean** | 1 | Unite | Large Button | LG |
|  | 1 | Subtract | Large Button | LG |
|  | 1 | Intersect | Large Button | LG |
|  | │ | — | *SubGroup 分隔线* | — |
|  | 2 | Split | Small Button | SM |
| **Object** | 1 | Group | Small Button | SM |
|  | 1 | Ungroup | Small Button | SM |
|  | 1 | Color | Small Button | SM |
|  | 1 | Transparency | Small Button | SM |

#### Simulation Tab

| Group | SubGroup | 按钮 | 类型 | 大小 |
|-------|----------|------|------|------|
| **Validate** | 1 | Validate | Large Button | LG |
| **Analysis** | 1 | Add Setup | Large Button | LG |
|  | 1 | Add Sweep | Small Button | SM |
| **Solve** | 1 | Analyze All | Large Button | LG |
|  | 1 | Solve | Large Button | LG |
|  | │ | — | *SubGroup 分隔线* | — |
|  | 2 | Abort | Small Button (disabled) | SM |

#### Results Tab

| Group | SubGroup | 按钮 | 类型 | 大小 |
|-------|----------|------|------|------|
| **Report** | 1 | Create Report | Large Split (类型选择) | LG |
|  | 1 | Solution Data | Large Button | LG |
| **Field Overlays** | 1 | Plot Fields | Large Split (E/H/SAR) | LG |
|  | 1 | Animate | Small Button | SM |

#### View Tab (3D Modeler 模式)

| Group | SubGroup | 按钮 | 类型 | 大小 |
|-------|----------|------|------|------|
| **Visibility** | 1 | Grid | Toggle (Small) | SM |
|  | 1 | Ruler | Toggle (Small) | SM |
|  | 1 | Axes | Toggle (Small) | SM |
| **Render** | 1 | Shaded | Toggle (Small) | SM |
|  | 1 | Wireframe | Toggle (Small) | SM |
| **Zoom** | 1 | Fit All | Large Button | LG |
|  | 1 | Zoom In | Small Button | SM |
|  | 1 | Zoom Out | Small Button | SM |

### 2.13 传统菜单栏（与 Ribbon 并存）

Ribbon 上方保留传统菜单栏：**File, Edit, View, Project, Draw, Modeler, [求解器菜单], Tools, Window, Help**

求解器菜单名称随设计类型变化（如 "HFSS" 或 "Maxwell3D"）。

---

## 3. Dock 系统

### 3.1 核心 Dock 面板

| 面板 | 功能 | 默认位置 |
|------|------|----------|
| **Project Manager** | 项目/设计树形视图 | 左侧 |
| **Properties Window** | 选中对象的上下文属性编辑 | 左侧（与 Project Manager 标签页组） |
| **Message Manager** | 警告、错误、信息消息 | 底部 |
| **Progress Window** | 仿真进度条 | 底部 |
| **3D Modeler / Layout** | 中央几何或布局编辑器 | 中央主区域 |

### 3.2 Dock 操作行为

| 操作 | 说明 |
|------|------|
| **拖放停靠** | 拖拽面板标题栏时，出现停靠引导菱形（上/下/左/右/中心） |
| **分割视图** | 放置到其他窗口的上/下/左/右目标，创建并排或堆叠布局 |
| **标签组合** | 放置到另一窗口的中心目标，合并为标签组（底部显示 Tab 标签）。中央 3D Modeler 不可与其他面板合并 |
| **浮动窗口** | 拖拽释放到任何停靠目标之外，可自由浮动（支持第二显示器） |
| **双击标题栏** | 切换 docked/floating 状态 |
| **标签重排** | 标签组内的 Tab 可通过拖拽左右重排 |

### 3.3 Auto-Hide（自动隐藏 / Pin 机制）

每个可停靠面板标题栏都有一个 **Pin 图标**：

- **Pinned**（竖直图钉）：面板始终可见
- **Unpinned**（水平图钉，Auto-Hide 开启）：面板折叠为边缘的窄标签；鼠标悬停或点击标签临时显示
- 右键标题栏 → "Auto Hide" 可切换
- **标签组共同折叠**：Auto-Hide 一个标签组时，组内所有标签一起折叠

### 3.4 布局持久化

- 面板布局保存在 Windows 注册表中：`HKCU\Software\Ansoft\ElectronicsDesktop\<version>\Desktop`
- 可通过 "Window > Restore Default Layout"（2023 R2+）恢复默认布局
- 也可删除注册表中的 WindowPositions/DockingLayout 键来重置

---

## 4. 配色方案与视觉风格

### 4.1 主题

从 2024 R1 开始提供三种配色方案（Dark/Light 在 2024 R2 为 Windows-only beta）：

| 主题 | 说明 |
|------|------|
| **Classic** | 传统外观，浅灰与白色，默认主题 |
| **Dark** | 深色主题，减轻视觉疲劳（beta） |
| **Light** | 现代浅色主题替代方案（beta） |

切换路径：Desktop Tab → General Options，或 Tools → Options → General Options → Desktop Configuration → User Interface（需重启）

### 4.2 Classic 主题近似配色

| 元素 | 近似颜色 |
|------|----------|
| 窗口边框 / Ribbon 背景 | 浅灰 `#F0F0F0` ~ `#E8E8E8` |
| 工具栏 / 菜单栏 | 稍深灰 `#D6D6D6` ~ `#DCDCDC` |
| 选中图标高亮 | 中灰 `#B0B0B0` ~ `#BCBCBC` |
| 活跃 Tab 文字 | Ansys 蓝 `#0055A5` ~ `#0070C0` |
| 非活跃 Tab 文字 | 深灰 `#444444` ~ `#5A5A5A` |
| 面板/画布背景 | 白色 `#FFFFFF` |
| 边框 / 分隔线 | 银灰 `#C0C0C0` ~ `#A0A0A0` |

整体风格接近传统 Windows Office 风格的银灰色 Ribbon 外观。

---

## 5. 关键界面元素

### 5.1 Project Manager 树结构

以 HFSS 设计为例的典型层级结构：

```
Project (.aedt)
  └── HFSSDesign
       ├── Model
       │    ├── 3D Objects
       │    ├── Coordinate Systems
       │    ├── 3D Components
       │    └── Materials
       ├── Boundaries (PEC, PMC, Radiation, Impedance ...)
       ├── Excitations (Wave Ports, Lumped Ports, Incident Wave ...)
       ├── Mesh Operations (Length-based, Skin Depth, Curvilinear)
       ├── Analysis
       │    └── Setup1
       │         ├── Adaptive Mesh Settings
       │         └── Sweep1 (Discrete / Interpolating / Fast)
       ├── Optimetrics (Parametric, Optimization, Sensitivity, Tuning)
       ├── Field Overlays (E-field, H-field, SAR)
       └── Results (S-params, Reports, Far-field)
```

右键任意树节点提供上下文操作菜单。

### 5.2 Properties Window

| Tab | 内容 |
|-----|------|
| **Attributes** | 对象名、材料、颜色、透明度、Solve Inside、Model/Non-Model 切换 |
| **Command/Definition** | 几何参数（位置、尺寸、轴）、坐标系、操作历史 |
| **Boundary/Excitation** | 边界/激励类型、值、方向、频率依赖设置 |

部分属性编辑通过双击弹出模态对话框（如设计属性、项目变量、某些 Setup 参数）。

### 5.3 3D Modeler 视口

| 元素 | 说明 |
|------|------|
| **Orientation Gadget（方向魔方）** | 视口角落的交互式立方体；点击面/边/角切换标准视图（Top, Front, Iso 等） |
| **坐标轴** | 原点处显示三轴（可通过 View → Coordinate System 切换） |
| **网格** | 可通过 View → Grid 开关；间距在 Modeler 设置中调整 |
| **视图控制** | 滚轮缩放；中键拖拽旋转；Ctrl+中键平移；View → Fit All 重置 |
| **选择模式** | Object / Face / Edge / Vertex 四种选择模式，从 Modeler 工具栏或菜单切换 |
| **渲染模式** | 线框、着色、透明等渲染风格选项 |

### 5.4 Status Bar（状态栏）

位于窗口最底部，包含：

- **Show/Hide Messages** 按钮 —— 切换 Message Manager
- **Show/Hide Progress** 按钮 —— 切换 Progress 窗口
- **X / Y / Z 坐标框** —— 绘图操作时显示
- **坐标模式控制** —— 绝对/相对坐标输入切换
- **模型单位显示** —— 当前单位系统
- 最顶部活跃仿真的进度指示

### 5.5 Message Manager

消息分三类显示：
- **Error**（阻止仿真的严重问题）
- **Warning**（可能影响结果的潜在问题）
- **Info**（常规状态/进度信息）

---

## 6. 设计要点总结

| 特征 | AEDT 实现方式 |
|------|---------------|
| **Ribbon + 菜单共存** | 保留传统菜单栏的同时提供 Ribbon，两套入口并行 |
| **上下文敏感 Ribbon** | Tab 随当前设计类型/编辑器动态增减 |
| **Group + SubGroup 层级** | Group 底部有 label，内含多个 SubGroup 用竖线分隔；实现按钮的两级分组 |
| **Ribbon 自适应** | 窗口缩小时图标压缩、Group 折叠为下拉按钮 |
| **灵活 Dock** | 支持拖放停靠、分割视图、标签组、浮动、Auto-Hide |
| **Pin/Unpin 机制** | 面板可 Pin（固定可见）或 Unpin（Auto-Hide 到边缘标签） |
| **Project Tree 核心地位** | 左侧树状结构是主要的项目导航和操作入口 |
| **属性面板联动** | 选中对象后 Properties Window 实时更新 |
| **多主题支持** | Classic（默认灰白）+ Dark/Light 主题（2024 R1+） |
| **布局持久化** | 用户窗口布局保存在注册表中，支持恢复默认 |

---

## 参考来源

- [Working with Ribbons - Ansys Help](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/Maxwell/Content/GettingStarted/WorkingWithRibbons.htm)
- [Electronics Desktop Windows Overview](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/Mechanical/Content/GettingStarted/ProjectWindowsOverview.htm)
- [Choosing a Color Scheme - Ansys Help](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Maxwell/Content/ChoosingColorScheme.htm)
- [Intro to AEDT User Interface - Ansys Innovation Courses](https://innovationspace.ansys.com/courses/courses/etm-using-ansys-hfss-and-icepak/lessons/intro-to-aedt-user-interface-lesson-1-7/)
- [Customizing the Automation Tab - Ansys Help](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Icepak/Content/GettingStarted/CustomizingtheAutomationTab.htm)
