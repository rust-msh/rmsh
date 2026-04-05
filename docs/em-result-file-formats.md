# EMStudio 仿真结果文件格式设计

## 1. 概述

本文档定义 EMStudio 仿真过程中产生的各类结果文件的详细格式。所有结果文件存放在 `.emsp.results/` 目录下，与工程文件 `.emsp` 分离存储。

### 1.1 设计原则

- **文本 vs 二进制分治**：元数据/索引/小体量结果用 JSON 文本（人类可读、可版本控制）；大体量数值数据用自定义二进制（紧凑、高效、支持内存映射）
- **自描述**：每个文件都包含版本号和足够的元信息，可独立解析
- **流式写入**：仿真过程中生成的文件（收敛历史、求解器日志）支持追加写入，中途崩溃不会丢失已有数据
- **按需加载**：二进制场数据支持分频点/分分量随机读取，不需要加载整个文件

### 1.2 文件格式总览

| 文件 | 格式 | 扩展名 | 体量 | 读写模式 |
|------|------|--------|------|---------|
| 验证报告 | JSON | `.json` | KB | 一次写入 |
| 收敛历史 | JSON | `.json` | KB | 追加写入 |
| 网格统计 | JSON | `.json` | KB | 一次写入 |
| 性能剖面 | JSON | `.json` | KB | 追加写入 |
| S 参数 | JSON | `.json` | KB~MB | 一次写入 |
| S 参数 | Touchstone | `.snp` | KB~MB | 一次写入 |
| 远场数据 | JSON | `.json` | MB | 一次写入 |
| 近场数据 | JSON | `.json` | MB | 一次写入 |
| 天线参数 | JSON | `.json` | KB | 一次写入 |
| Optimetrics 汇总 | JSON | `.json` | KB~MB | 追加写入 |
| 网格数据 | Gmsh MSH 4.1 | `.msh` | MB~GB | 一次写入，随机读 |
| 场解数据 | 二进制 | `.emsfld` | MB~GB | 一次写入，随机读 |
| 求解器日志 | 纯文本 | `.log` | KB~MB | 追加写入 |
| RLCG 矩阵数据 | JSON | `.json` | KB~MB | 一次写入 |
| 等效电路模型 | SPICE | `.sp` | KB | 一次写入 |
| 报告导出 | CSV | `.csv` | KB | 一次写入 |
| 场图导出 | PNG | `.png` | KB~MB | 一次写入 |

---

## 2. JSON 格式文件

### 2.1 验证报告（validation_report.json）

在 Analysis 求解前自动运行模型验证，检查几何合法性、边界/激励完整性、材料引用等问题。

```json
{
  "format_version": "1.0",
  "file_type": "ValidationReport",
  "design_id": "design-001",
  "validated_at": "2026-04-04T14:25:00Z",
  "overall_status": "Warning",
  "checks": [
    {
      "category": "Geometry",
      "check": "IntersectionCheck",
      "status": "Pass",
      "message": "No intersecting objects found"
    },
    {
      "category": "Geometry",
      "check": "SmallFeatureCheck",
      "status": "Warning",
      "message": "Object 'Via1' has feature size 0.02mm < lambda/100 at solution frequency",
      "details": {
        "object": "Via1",
        "feature_size_mm": 0.02,
        "lambda_fraction": 0.00016
      }
    },
    {
      "category": "Boundaries",
      "check": "RadiationBoundarySize",
      "status": "Pass",
      "message": "Radiation boundary is >= lambda/4 from nearest object"
    },
    {
      "category": "Boundaries",
      "check": "UnassignedFaces",
      "status": "Pass",
      "message": "All outer faces have boundary assignments"
    },
    {
      "category": "Excitations",
      "check": "PortOnBoundary",
      "status": "Pass",
      "message": "All wave ports are on outer boundary faces"
    },
    {
      "category": "Excitations",
      "check": "PortModeCount",
      "status": "Pass",
      "message": "Requested modes can be solved for all ports"
    },
    {
      "category": "Materials",
      "check": "UnresolvedReference",
      "status": "Pass",
      "message": "All material references resolved"
    },
    {
      "category": "Variables",
      "check": "ExpressionEvaluation",
      "status": "Pass",
      "message": "All variable expressions evaluate successfully"
    },
    {
      "category": "Mesh",
      "check": "MeshabilityEstimate",
      "status": "Pass",
      "message": "Estimated initial mesh: ~5000 tetrahedra",
      "details": {
        "estimated_elements": 5000,
        "estimated_memory_mb": 128
      }
    }
  ],
  "errors": 0,
  "warnings": 1,
  "info": 0
}
```

**检查类别（category）**：

| 类别 | 检查项 | 说明 |
|------|--------|------|
| `Geometry` | `IntersectionCheck` | 对象间非法相交 |
| `Geometry` | `SmallFeatureCheck` | 小于 lambda/100 的特征尺寸 |
| `Geometry` | `OpenSheet` | 未封闭的面片 |
| `Boundaries` | `RadiationBoundarySize` | 辐射边界距模型距离是否足够 |
| `Boundaries` | `UnassignedFaces` | 外表面无边界分配 |
| `Boundaries` | `DuplicateAssignment` | 同一面分配了多个冲突边界 |
| `Excitations` | `PortOnBoundary` | 波端口是否在外边界面上 |
| `Excitations` | `PortModeCount` | 端口模式数是否合理 |
| `Excitations` | `IntegrationLine` | 集总端口积分线方向检查 |
| `Materials` | `UnresolvedReference` | 引用了不存在的材料 |
| `Materials` | `NegativeProperty` | 材料属性值非物理（如负介电常数） |
| `Variables` | `ExpressionEvaluation` | 变量表达式是否可求值 |
| `Variables` | `CircularDependency` | 变量循环引用 |
| `Mesh` | `MeshabilityEstimate` | 网格生成可行性预估 |
| `Nets` | `NetConductorAssignment` | Q3D：所有导体对象是否已分配到网络 |
| `Nets` | `SourceSinkPairing` | Q3D：每个信号网络是否至少有一对源-汇端子 |
| `Nets` | `GroundNetDefined` | Q3D：是否定义了接地参考网络 |
| `Nets` | `FloatingConductor` | Q3D：检测未分配任何网络的导体对象 |
| `Boundaries` | `OpenBoundaryDistance` | Q3D：开放边界与导体的距离是否足够 |
| `Boundaries` | `ThinConductorThickness` | Q3D：薄导体厚度是否远小于导体宽度 |

**状态（status）枚举**：`Pass` | `Warning` | `Error` | `Info`

---

### 2.2 收敛历史（convergence.json）

每轮自适应完成后追加一条记录，支持仿真过程中实时读取以更新 UI 上的收敛曲线。

```json
{
  "format_version": "1.0",
  "file_type": "ConvergenceHistory",
  "design_id": "design-001",
  "setup": "Setup1",
  "solution_frequency": "2.4GHz",
  "target_max_delta_s": 0.02,
  "target_min_converged_passes": 2,
  "max_passes_limit": 15,
  "passes": [
    {
      "pass_number": 1,
      "timestamp": "2026-04-04T14:25:15Z",
      "mesh": {
        "num_tetrahedra": 5420,
        "num_nodes": 1210,
        "min_edge_length_mm": 0.05,
        "max_edge_length_mm": 8.2,
        "mean_edge_length_mm": 2.1
      },
      "solution": {
        "max_delta_s": null,
        "delta_energy": null,
        "matrix_size": 16260,
        "num_rhs": 1
      },
      "performance": {
        "mesh_time_sec": 1.2,
        "solve_time_sec": 8.5,
        "error_estimation_time_sec": 0.3,
        "total_time_sec": 10.0,
        "peak_memory_mb": 256
      }
    },
    {
      "pass_number": 2,
      "timestamp": "2026-04-04T14:25:25Z",
      "mesh": {
        "num_tetrahedra": 8103,
        "num_nodes": 1852,
        "min_edge_length_mm": 0.03,
        "max_edge_length_mm": 8.2,
        "mean_edge_length_mm": 1.7
      },
      "solution": {
        "max_delta_s": 0.045,
        "delta_energy": 0.032,
        "matrix_size": 24309,
        "num_rhs": 1
      },
      "performance": {
        "mesh_time_sec": 1.8,
        "solve_time_sec": 14.2,
        "error_estimation_time_sec": 0.5,
        "total_time_sec": 16.5,
        "peak_memory_mb": 412
      }
    },
    {
      "pass_number": 3,
      "timestamp": "2026-04-04T14:25:42Z",
      "mesh": {
        "num_tetrahedra": 10250,
        "num_nodes": 2380,
        "min_edge_length_mm": 0.02,
        "max_edge_length_mm": 8.2,
        "mean_edge_length_mm": 1.4
      },
      "solution": {
        "max_delta_s": 0.012,
        "delta_energy": 0.008,
        "matrix_size": 30750,
        "num_rhs": 1
      },
      "performance": {
        "mesh_time_sec": 2.1,
        "solve_time_sec": 20.3,
        "error_estimation_time_sec": 0.7,
        "total_time_sec": 23.1,
        "peak_memory_mb": 580
      }
    }
  ],
  "result": {
    "converged": true,
    "converged_at_pass": 3,
    "consecutive_converged_passes": 2,
    "total_elapsed_sec": 49.6,
    "final_max_delta_s": 0.012
  }
}
```

> **Q3D 收敛历史差异**：Q3D 准静态求解的收敛指标为 `max_delta_energy`（能量变化百分比）而非 HFSS 的 `max_delta_s`（S 参数变化）。此外，Q3D 收敛历史中还记录每轮提取的 RLCG 矩阵摘要，以便观察矩阵值随网格加密的稳定趋势。Q3D 的 `solution` 字段使用 `delta_energy` 替代 `max_delta_s` 作为主要收敛判据。

**Q3D 收敛历史示例**：

```json
{
  "format_version": "1.0",
  "file_type": "ConvergenceHistory",
  "design_id": "design-002",
  "setup": "Q3D_Setup1",
  "solution_type": "Q3D_ACRL",
  "adaptive_frequency": "1GHz",
  "target_max_delta_energy": 0.02,
  "target_min_converged_passes": 2,
  "max_passes_limit": 10,
  "passes": [
    {
      "pass_number": 1,
      "timestamp": "2026-04-04T15:10:00Z",
      "mesh": {
        "num_triangles": 3200,
        "num_tetrahedra": 0,
        "num_nodes": 1680,
        "min_edge_length_mm": 0.01,
        "max_edge_length_mm": 5.0,
        "mean_edge_length_mm": 0.8
      },
      "solution": {
        "max_delta_s": null,
        "delta_energy": null,
        "matrix_size": 9600,
        "num_rhs": 3
      },
      "rlcg_snapshot": {
        "R_max_ohm": 0.125,
        "L_max_nH": 2.85,
        "C_max_pF": 0.42,
        "G_max_mS": 0.0012
      },
      "performance": {
        "mesh_time_sec": 0.8,
        "solve_time_sec": 5.2,
        "error_estimation_time_sec": 0.2,
        "total_time_sec": 6.2,
        "peak_memory_mb": 180
      }
    },
    {
      "pass_number": 2,
      "timestamp": "2026-04-04T15:10:06Z",
      "mesh": {
        "num_triangles": 4800,
        "num_tetrahedra": 0,
        "num_nodes": 2520,
        "min_edge_length_mm": 0.008,
        "max_edge_length_mm": 5.0,
        "mean_edge_length_mm": 0.65
      },
      "solution": {
        "max_delta_s": null,
        "delta_energy": 0.035,
        "matrix_size": 14400,
        "num_rhs": 3
      },
      "rlcg_snapshot": {
        "R_max_ohm": 0.128,
        "L_max_nH": 2.91,
        "C_max_pF": 0.44,
        "G_max_mS": 0.0013
      },
      "performance": {
        "mesh_time_sec": 1.1,
        "solve_time_sec": 8.8,
        "error_estimation_time_sec": 0.3,
        "total_time_sec": 10.2,
        "peak_memory_mb": 290
      }
    },
    {
      "pass_number": 3,
      "timestamp": "2026-04-04T15:10:16Z",
      "mesh": {
        "num_triangles": 6100,
        "num_tetrahedra": 0,
        "num_nodes": 3200,
        "min_edge_length_mm": 0.005,
        "max_edge_length_mm": 5.0,
        "mean_edge_length_mm": 0.52
      },
      "solution": {
        "max_delta_s": null,
        "delta_energy": 0.015,
        "matrix_size": 18300,
        "num_rhs": 3
      },
      "rlcg_snapshot": {
        "R_max_ohm": 0.129,
        "L_max_nH": 2.93,
        "C_max_pF": 0.445,
        "G_max_mS": 0.0013
      },
      "performance": {
        "mesh_time_sec": 1.5,
        "solve_time_sec": 14.2,
        "error_estimation_time_sec": 0.4,
        "total_time_sec": 16.1,
        "peak_memory_mb": 420
      }
    }
  ],
  "result": {
    "converged": true,
    "converged_at_pass": 3,
    "consecutive_converged_passes": 2,
    "total_elapsed_sec": 32.5,
    "final_max_delta_energy": 0.015
  }
}
```

> **Q3D 网格特点**：Q3D 的 MoM（矩量法）求解器主要使用**三角形面网格**（`num_triangles`）而非体积四面体（`num_tetrahedra`），因此面网格单元数是 Q3D 收敛的关键指标。DC 电阻提取使用 FEM 时才会产生四面体。`rlcg_snapshot` 记录每轮的 RLCG 矩阵极值，用于直观判断矩阵值是否趋于稳定。

---

### 2.3 网格统计（mesh_stats.json）

最终收敛网格的详细统计信息，用于质量评估。

```json
{
  "format_version": "1.0",
  "file_type": "MeshStatistics",
  "design_id": "design-001",
  "setup": "Setup1",
  "generated_at": "2026-04-04T14:25:42Z",
  "converged_pass": 3,
  "global": {
    "num_tetrahedra": 10250,
    "num_nodes": 2380,
    "num_boundary_triangles": 4120,
    "num_edges": 15200,
    "element_order": "Mixed"
  },
  "edge_length": {
    "unit": "mm",
    "min": 0.02,
    "max": 8.2,
    "mean": 1.4,
    "rms": 1.8,
    "histogram": {
      "bins_mm": [0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0],
      "counts": [120, 1850, 3400, 3100, 1600, 180]
    }
  },
  "element_volume": {
    "unit": "mm3",
    "min": 0.00001,
    "max": 12.5,
    "mean": 0.82,
    "histogram": {
      "bins_mm3": [0, 0.01, 0.1, 1.0, 10.0, 100.0],
      "counts": [45, 890, 4200, 4800, 315]
    }
  },
  "element_quality": {
    "metric": "AspectRatio",
    "min": 1.0,
    "max": 28.5,
    "mean": 2.3,
    "elements_above_10": 42,
    "elements_above_20": 3
  },
  "per_object": [
    {
      "object": "Substrate",
      "num_tetrahedra": 4200,
      "mean_edge_length_mm": 1.2,
      "material": "FR4_epoxy"
    },
    {
      "object": "Patch",
      "num_tetrahedra": 1800,
      "mean_edge_length_mm": 0.8,
      "material": "copper"
    },
    {
      "object": "AirBox",
      "num_tetrahedra": 3500,
      "mean_edge_length_mm": 2.5,
      "material": "vacuum"
    }
  ]
}
```

---

### 2.4 性能剖面（profile.json）

记录仿真各阶段的精确耗时和资源消耗，用于性能分析和 HPC 规模评估。

```json
{
  "format_version": "1.0",
  "file_type": "SolveProfile",
  "design_id": "design-001",
  "setup": "Setup1",
  "started_at": "2026-04-04T14:25:00Z",
  "finished_at": "2026-04-04T14:28:30Z",
  "environment": {
    "cpu_model": "Intel Xeon W-2295",
    "cpu_cores_available": 18,
    "cpu_cores_used": 8,
    "ram_total_gb": 64,
    "os": "linux-x86_64",
    "solver_version": "emstudio-solver 0.1.0"
  },
  "phases": [
    {
      "name": "Validation",
      "elapsed_sec": 0.5,
      "peak_memory_mb": 64,
      "status": "OK"
    },
    {
      "name": "InitialMesh",
      "elapsed_sec": 2.3,
      "peak_memory_mb": 128,
      "details": {
        "num_tetrahedra": 5420
      }
    },
    {
      "name": "AdaptiveSolve",
      "elapsed_sec": 49.6,
      "peak_memory_mb": 580,
      "details": {
        "num_passes": 3,
        "converged": true,
        "per_pass": [
          {
            "pass": 1,
            "assembly_sec": 2.5,
            "factorization_sec": 5.0,
            "solve_sec": 0.8,
            "field_recovery_sec": 0.5,
            "error_estimation_sec": 0.3,
            "mesh_refinement_sec": 0.9,
            "total_sec": 10.0,
            "peak_memory_mb": 256,
            "matrix_nonzeros": 482000
          },
          {
            "pass": 2,
            "assembly_sec": 3.8,
            "factorization_sec": 9.2,
            "solve_sec": 1.2,
            "field_recovery_sec": 0.8,
            "error_estimation_sec": 0.5,
            "mesh_refinement_sec": 1.0,
            "total_sec": 16.5,
            "peak_memory_mb": 412,
            "matrix_nonzeros": 735000
          },
          {
            "pass": 3,
            "assembly_sec": 5.0,
            "factorization_sec": 12.5,
            "solve_sec": 1.8,
            "field_recovery_sec": 1.1,
            "error_estimation_sec": 0.7,
            "mesh_refinement_sec": 2.0,
            "total_sec": 23.1,
            "peak_memory_mb": 580,
            "matrix_nonzeros": 920000
          }
        ]
      }
    },
    {
      "name": "FrequencySweep",
      "sweep_name": "Sweep1",
      "elapsed_sec": 85.0,
      "peak_memory_mb": 620,
      "details": {
        "type": "Interpolating",
        "num_frequency_points": 301,
        "num_actual_solves": 12,
        "interpolation_error": 0.001
      }
    },
    {
      "name": "FarFieldComputation",
      "elapsed_sec": 5.2,
      "peak_memory_mb": 300,
      "details": {
        "setup": "InfiniteSphere1",
        "theta_points": 181,
        "phi_points": 361
      }
    }
  ],
  "totals": {
    "elapsed_sec": 142.6,
    "peak_memory_mb": 620,
    "cpu_time_sec": 980.5
  }
}
```

---

### 2.5 S 参数数据（s_parameters.json）

存储网络参数矩阵数据。自适应频率处的结果和扫频结果使用相同格式。

```json
{
  "format_version": "1.0",
  "file_type": "SParameterData",
  "design_id": "design-001",
  "setup": "Setup1",
  "sweep": "Sweep1",
  "solution_type": "DrivenModal",
  "reference_impedance_ohm": 50.0,
  "num_ports": 1,
  "port_names": ["Port1"],
  "num_frequencies": 301,
  "frequency_unit": "GHz",
  "data_format": "RealImaginary",
  "frequencies": [1.0, 1.01, 1.02, "...（301个频点）"],
  "parameters": {
    "S11": {
      "real": [-0.85, -0.84, -0.83, "..."],
      "imag": [0.12, 0.13, 0.14, "..."]
    }
  },
  "derived": {
    "S11_magnitude_db": [-1.4, -1.5, -1.6, "..."],
    "S11_phase_deg": [172.0, 171.2, 170.3, "..."],
    "S11_vswr": [13.0, 12.3, 11.5, "..."],
    "Z11_real_ohm": [22.5, 23.1, 23.8, "..."],
    "Z11_imag_ohm": [-35.2, -33.8, -32.1, "..."],
    "group_delay_ns": [0.12, 0.13, 0.14, "..."]
  }
}
```

> **多端口示例**（2 端口）：`parameters` 中包含 `S11`、`S12`、`S21`、`S22` 四组数据。

**`data_format`** 枚举：
- `RealImaginary`：实部+虚部（精度无损，推荐内部存储）
- `MagnitudeAngle`：幅度+相位角（度）
- `dBAngle`：dB 幅度+相位角（度）

---

### 2.6 Touchstone 导出（s_parameters.snp）

遵循 [Touchstone 2.0 规范](https://ibis.org/touchstone_ver2.0/touchstone_ver2_0.pdf)，用于与第三方工具互操作。

**单端口示例（.s1p）**：
```
[Version] 2.0
! EMStudio S-Parameter Export
! Project: MyAntenna.emsp
! Design: Patch Antenna 2.4GHz
! Setup: Setup1 / Sweep1
! Generated: 2026-04-04T14:30:00Z
[Number of Ports] 1
[Number of Frequencies] 301
[Reference]
50.0
# GHz S RI R 50
[Network Data]
1.000000000  -0.850000  0.120000
1.010000000  -0.840000  0.130000
1.020000000  -0.830000  0.140000
! ... 298 more lines ...
4.000000000  -0.780000  0.180000
[End]
```

**双端口示例（.s2p）**：
```
[Version] 2.0
[Number of Ports] 2
[Two-Port Data Order] 21_12
[Reference]
50.0 50.0
# GHz S RI R 50
[Network Data]
! freq       S11_re     S11_im     S21_re     S21_im     S12_re     S12_im     S22_re     S22_im
1.000000000  -0.850000  0.120000   0.005000   0.001000   0.005000   0.001000   -0.820000  0.110000
[End]
```

---

### 2.7 远场数据（far_field_*.json）

每个远场设置 × 每个频点生成一个文件。

```json
{
  "format_version": "1.0",
  "file_type": "FarFieldData",
  "design_id": "design-001",
  "setup": "Setup1",
  "far_field_setup": "InfiniteSphere1",
  "frequency": "2.4GHz",
  "coordinate_system": "Global",
  "reference_impedance_ohm": 50.0,
  "theta": {
    "start_deg": 0,
    "stop_deg": 180,
    "step_deg": 1,
    "num_points": 181
  },
  "phi": {
    "start_deg": 0,
    "stop_deg": 360,
    "step_deg": 1,
    "num_points": 361
  },
  "fields": {
    "E_theta": {
      "description": "远场 Theta 分量（复数）",
      "unit": "V",
      "data_real": [["2D array: [phi][theta]"]],
      "data_imag": [["2D array: [phi][theta]"]]
    },
    "E_phi": {
      "description": "远场 Phi 分量（复数）",
      "unit": "V",
      "data_real": [["2D array: [phi][theta]"]],
      "data_imag": [["2D array: [phi][theta]"]]
    }
  },
  "derived_quantities": {
    "GainTotal": {
      "description": "总增益",
      "unit": "dBi",
      "data": [["2D array: [phi][theta]"]]
    },
    "GainTheta": {
      "unit": "dBi",
      "data": [["2D array"]]
    },
    "GainPhi": {
      "unit": "dBi",
      "data": [["2D array"]]
    },
    "DirectivityTotal": {
      "unit": "dBi",
      "data": [["2D array"]]
    },
    "AxialRatio": {
      "description": "轴比（线极化趋近无穷大）",
      "unit": "dB",
      "data": [["2D array"]]
    },
    "Polarization": {
      "description": "极化类型标记",
      "values": "LHCP|RHCP|Linear",
      "data": [["2D string array"]]
    }
  },
  "antenna_parameters": {
    "peak_gain_dbi": 7.2,
    "peak_gain_theta_deg": 0,
    "peak_gain_phi_deg": 0,
    "peak_directivity_dbi": 7.5,
    "radiation_efficiency": 0.93,
    "total_efficiency": 0.91,
    "beamwidth_e_plane_deg": 78.0,
    "beamwidth_h_plane_deg": 85.0,
    "cross_pol_level_db": -25.3,
    "front_to_back_ratio_db": 15.3,
    "side_lobe_level_db": -12.5,
    "radiated_power_w": 0.0093,
    "accepted_power_w": 0.01,
    "incident_power_w": 0.01
  }
}
```

---

### 2.8 近场数据（near_field_*.json）

每个近场设置 × 每个频点生成一个文件。

```json
{
  "format_version": "1.0",
  "file_type": "NearFieldData",
  "design_id": "design-001",
  "setup": "Setup1",
  "near_field_setup": "NearFieldLine1",
  "frequency": "2.4GHz",
  "setup_type": "Line",
  "geometry": {
    "start_point_mm": [0, 0, 10],
    "end_point_mm": [100, 0, 10],
    "num_points": 201
  },
  "sampling_points_mm": [
    [0, 0, 10], [0.5, 0, 10], [1.0, 0, 10], "...（201 个点）"
  ],
  "fields": {
    "E": {
      "unit": "V/m",
      "x_real": [12.5, 12.3, "..."],
      "x_imag": [-3.2, -3.1, "..."],
      "y_real": [0.1, 0.1, "..."],
      "y_imag": [0.0, 0.0, "..."],
      "z_real": [5.8, 5.7, "..."],
      "z_imag": [-1.2, -1.1, "..."]
    },
    "H": {
      "unit": "A/m",
      "x_real": [0.03, 0.03, "..."],
      "x_imag": [-0.01, -0.01, "..."],
      "y_real": [0.08, 0.08, "..."],
      "y_imag": [-0.02, -0.02, "..."],
      "z_real": [0.0, 0.0, "..."],
      "z_imag": [0.0, 0.0, "..."]
    }
  },
  "derived": {
    "E_magnitude": {
      "unit": "V/m",
      "data": [13.8, 13.6, "..."]
    },
    "H_magnitude": {
      "unit": "A/m",
      "data": [0.088, 0.087, "..."]
    },
    "Poynting_magnitude": {
      "unit": "W/m2",
      "data": [0.61, 0.59, "..."]
    }
  }
}
```

**矩形面近场**（`setup_type: "Rectangle"`）的 `geometry` 和 `sampling_points` 为 2D 网格，`fields` 中的数据数组为扁平化的行优先 2D 数组。

---

### 2.9 RLCG 矩阵数据（rlcg_matrix.json）— Q3D 专用

Q3D 准静态求解的核心输出是 RLCG 矩阵——描述导体网络之间的寄生电阻（R）、电感（L）、电容（C）和电导（G）耦合关系。参考 [Viewing Matrix Data in Q3D](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/ViewingMatrixDatainQ3D.htm)。

**单频点 RLCG 矩阵数据**：

```json
{
  "format_version": "1.0",
  "file_type": "RLCGMatrixData",
  "design_id": "design-002",
  "setup": "Q3D_Setup1",
  "solution_type": "Q3D_ACRL",
  "num_nets": 3,
  "net_names": ["Signal1", "Signal2", "Ground"],
  "num_terminals": 4,
  "terminal_names": ["Signal1:T1", "Signal1:T2", "Signal2:T3", "Signal2:T4"],
  "ground_net": "Ground",
  "num_frequencies": 50,
  "frequency_unit": "GHz",
  "frequencies": [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0],
  "matrices": {
    "R": {
      "description": "电阻矩阵（AC + DC）",
      "unit": "ohm",
      "data_per_frequency": [
        {
          "frequency": 0.01,
          "matrix": [
            [0.125, 0.003, 0.001, 0.0005],
            [0.003, 0.125, 0.0005, 0.001],
            [0.001, 0.0005, 0.118, 0.002],
            [0.0005, 0.001, 0.002, 0.118]
          ]
        },
        {
          "frequency": 1.0,
          "matrix": [
            [0.285, 0.012, 0.005, 0.002],
            [0.012, 0.285, 0.002, 0.005],
            [0.005, 0.002, 0.268, 0.010],
            [0.002, 0.005, 0.010, 0.268]
          ]
        }
      ]
    },
    "L": {
      "description": "电感矩阵",
      "unit": "nH",
      "data_per_frequency": [
        {
          "frequency": 0.01,
          "matrix": [
            [2.85, 0.42, 0.15, 0.08],
            [0.42, 2.85, 0.08, 0.15],
            [0.15, 0.08, 2.62, 0.38],
            [0.08, 0.15, 0.38, 2.62]
          ]
        }
      ]
    },
    "C": {
      "description": "电容矩阵",
      "unit": "pF",
      "data_per_frequency": [
        {
          "frequency": 0.01,
          "matrix": [
            [0.445, -0.085, -0.032, -0.012],
            [-0.085, 0.445, -0.012, -0.032],
            [-0.032, -0.012, 0.420, -0.078],
            [-0.012, -0.032, -0.078, 0.420]
          ]
        }
      ]
    },
    "G": {
      "description": "电导矩阵（介质损耗）",
      "unit": "mS",
      "data_per_frequency": [
        {
          "frequency": 1.0,
          "matrix": [
            [0.0013, -0.0002, -0.0001, 0.0],
            [-0.0002, 0.0013, 0.0, -0.0001],
            [-0.0001, 0.0, 0.0012, -0.0002],
            [0.0, -0.0001, -0.0002, 0.0012]
          ]
        }
      ]
    }
  },
  "dc_data": {
    "R_dc": {
      "description": "DC 电阻矩阵",
      "unit": "ohm",
      "matrix": [
        [0.120, 0.002, 0.001, 0.0003],
        [0.002, 0.120, 0.0003, 0.001],
        [0.001, 0.0003, 0.115, 0.002],
        [0.0003, 0.001, 0.002, 0.115]
      ]
    },
    "L_dc": {
      "description": "DC 电感矩阵（低频极限）",
      "unit": "nH",
      "matrix": [
        [3.12, 0.48, 0.18, 0.09],
        [0.48, 3.12, 0.09, 0.18],
        [0.18, 0.09, 2.88, 0.42],
        [0.09, 0.18, 0.42, 2.88]
      ]
    }
  }
}
```

**矩阵数据说明**：

| 矩阵 | 求解类型 | 物理含义 | 对角线 | 非对角线 |
|------|---------|---------|--------|---------|
| R | Q3D_DCRL / Q3D_ACRL | 导体自阻/互阻 | 自电阻（Ω） | 互电阻（Ω） |
| L | Q3D_DCRL / Q3D_ACRL | 自感/互感 | 自感（nH） | 互感（nH） |
| C | Q3D_C / Q3D_CG | 电容矩阵 | 自电容（pF，正） | 互电容（pF，负） |
| G | Q3D_CG | 电导矩阵（介质损耗） | 自电导（mS） | 互电导（mS） |

> **矩阵对称性**：RLCG 矩阵均为对称矩阵（`M[i][j] == M[j][i]`），但 JSON 中存储完整矩阵以方便直接索引。C 矩阵的对角线为正值（表示该端子到所有其他端子和地的总电容），非对角线为负值（表示端子间的耦合电容）。

> **频率依赖性**：R 和 L 矩阵通常随频率变化（趋肤效应和邻近效应），因此以频率为索引存储。DC 值单独存储在 `dc_data` 中（频率 = 0 的极限值）。C 矩阵通常与频率无关（静电解），但 CG 模式下的 G 矩阵与频率成正比（介质损耗）。

---

### 2.10 等效电路模型（equivalent_circuit.sp）— Q3D 专用

参考 [Exporting Equivalent Circuit Data](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/ExportingQ3DExtractorEquivalentCircuitData.htm)：Q3D 可将提取的 RLCG 矩阵导出为 SPICE 等效电路网表，用于 SI/PI 仿真工具（如 HSPICE、Spectre、LTspice）。

**宽带集总模型（Broadband SPICE）**：

```spice
* EMStudio Q3D Equivalent Circuit Export
* Project: PCB_Interconnect.emsp
* Design: Signal_Trace_Design
* Setup: Q3D_Setup1 (Q3D_ACRL)
* Generated: 2026-04-04T15:30:00Z
* Frequency range: 10MHz - 5GHz
* Model type: Broadband Lumped

.SUBCKT Q3D_Model Signal1_src Signal1_sink Signal2_src Signal2_sink GND

* Self impedance: Signal1
R_self_1  Signal1_src  n1_1  0.285
L_self_1  n1_1         Signal1_sink  2.85n

* Self impedance: Signal2
R_self_2  Signal2_src  n2_1  0.268
L_self_2  n2_1         Signal2_sink  2.62n

* Mutual inductance: Signal1 <-> Signal2
K_12  L_self_1  L_self_2  0.147

* Capacitance: Signal1 to GND
C_1g  Signal1_src  GND  0.328p

* Capacitance: Signal2 to GND
C_2g  Signal2_src  GND  0.310p

* Mutual capacitance: Signal1 <-> Signal2
C_12  Signal1_src  Signal2_src  0.085p

.ENDS Q3D_Model
```

**导出配置**（JSON 元数据）：

```json
{
  "format_version": "1.0",
  "file_type": "EquivalentCircuitExport",
  "design_id": "design-002",
  "setup": "Q3D_Setup1",
  "export_format": "SPICE",
  "model_type": "BroadbandLumped",
  "frequency_range": {
    "start": "10MHz",
    "stop": "5GHz"
  },
  "nets_included": ["Signal1", "Signal2"],
  "ground_net": "Ground",
  "output_file": "equivalent_circuit.sp",
  "options": {
    "include_dc_resistance": true,
    "include_mutual_inductance": true,
    "include_mutual_capacitance": true,
    "include_dielectric_loss": false,
    "coupling_threshold": 0.01
  }
}
```

**导出模型类型**：

| 类型 | 说明 | 适用场景 |
|------|------|---------|
| `BroadbandLumped` | 宽带集总 RLC 模型 | 一般 SI 分析，简单互连 |
| `FrequencyDependentLumped` | 频率依赖集总模型（多级 RL 梯形网络） | 宽频 SI，需要反映趋肤效应 |
| `TLineModel` | 传输线 W-Element 模型 | 长走线、高速串行链路 |
| `SParameterBlock` | 基于 S 参数的行为模型（Touchstone） | 全频段精确匹配，信号完整性仿真 |

---

### 2.11 Optimetrics 汇总（summary.json）

参数扫描或优化完成后，汇总所有变量组合的关键输出。

**参数扫描汇总**：

```json
{
  "format_version": "1.0",
  "file_type": "OptimetricsSummary",
  "design_id": "design-001",
  "optimetrics_name": "LengthSweep",
  "type": "ParametricSweep",
  "setup": "Setup1",
  "started_at": "2026-04-04T15:00:00Z",
  "finished_at": "2026-04-04T16:30:00Z",
  "swept_variables": ["patch_l"],
  "output_variables": ["S11_at_center", "PeakGain", "BW_10dB"],
  "total_variations": 15,
  "completed_variations": 15,
  "failed_variations": 0,
  "variations": [
    {
      "index": 1,
      "variables": { "patch_l": "25.0mm" },
      "status": "Converged",
      "num_passes": 4,
      "outputs": {
        "S11_at_center": { "value": -8.2, "unit": "dB" },
        "PeakGain": { "value": 6.1, "unit": "dBi" },
        "BW_10dB": { "value": 0, "unit": "MHz" }
      },
      "result_path": "variation_001/"
    },
    {
      "index": 2,
      "variables": { "patch_l": "25.5mm" },
      "status": "Converged",
      "num_passes": 3,
      "outputs": {
        "S11_at_center": { "value": -12.5, "unit": "dB" },
        "PeakGain": { "value": 6.8, "unit": "dBi" },
        "BW_10dB": { "value": 42, "unit": "MHz" }
      },
      "result_path": "variation_002/"
    }
  ]
}
```

**优化汇总**：

```json
{
  "format_version": "1.0",
  "file_type": "OptimetricsSummary",
  "design_id": "design-001",
  "optimetrics_name": "MatchOptimize",
  "type": "Optimization",
  "setup": "Setup1",
  "algorithm": "QuasiNewton",
  "max_iterations": 50,
  "started_at": "2026-04-04T17:00:00Z",
  "finished_at": "2026-04-04T18:15:00Z",
  "optimized_variables": ["patch_l", "patch_w"],
  "goals": [
    { "name": "MinS11", "expression": "S11_at_center", "condition": "Minimize" },
    { "name": "S11_below_10dB", "expression": "dB(S(Port1,Port1))", "condition": "LessThan", "target": -10.0 }
  ],
  "total_iterations": 23,
  "converged": true,
  "convergence_history": [
    { "iteration": 1, "cost": -5.2, "variables": { "patch_l": "28.5mm", "patch_w": "37.0mm" } },
    { "iteration": 2, "cost": -10.8, "variables": { "patch_l": "29.1mm", "patch_w": "36.2mm" } },
    { "iteration": 23, "cost": -22.5, "variables": { "patch_l": "29.8mm", "patch_w": "35.5mm" } }
  ],
  "best_result": {
    "iteration": 23,
    "variables": { "patch_l": "29.8mm", "patch_w": "35.5mm" },
    "cost": -22.5,
    "outputs": {
      "S11_at_center": { "value": -22.5, "unit": "dB" },
      "PeakGain": { "value": 7.1, "unit": "dBi" },
      "BW_10dB": { "value": 95, "unit": "MHz" }
    },
    "result_path": "iteration_023/"
  }
}
```

---

### 2.12 求解器日志（solver.log）

纯文本格式，带时间戳的逐行日志，支持追加写入。

```
[2026-04-04T14:25:00.000Z] [INFO]  EMStudio Solver v0.1.0 started
[2026-04-04T14:25:00.001Z] [INFO]  Design: Patch Antenna 2.4GHz (design-001)
[2026-04-04T14:25:00.002Z] [INFO]  Setup: Setup1, Solution frequency: 2.4GHz
[2026-04-04T14:25:00.003Z] [INFO]  Solution type: DrivenModal, Max passes: 15, Max delta S: 0.02
[2026-04-04T14:25:00.100Z] [INFO]  === Validation Phase ===
[2026-04-04T14:25:00.500Z] [INFO]  Validation passed (0 errors, 1 warning)
[2026-04-04T14:25:00.501Z] [WARN]  Object 'Via1' has feature size 0.02mm < lambda/100
[2026-04-04T14:25:00.600Z] [INFO]  === Initial Mesh Generation ===
[2026-04-04T14:25:02.900Z] [INFO]  Initial mesh: 5420 tetrahedra, 1210 nodes
[2026-04-04T14:25:02.901Z] [INFO]  === Adaptive Pass 1 ===
[2026-04-04T14:25:03.100Z] [INFO]  Matrix assembly: 16260 DOFs, 482000 nonzeros
[2026-04-04T14:25:08.100Z] [INFO]  Direct solve completed (5.0 sec)
[2026-04-04T14:25:09.000Z] [INFO]  Field recovery completed
[2026-04-04T14:25:09.300Z] [INFO]  Error estimation: max delta E = 15.2%
[2026-04-04T14:25:10.200Z] [INFO]  Mesh refined: 5420 -> 8103 tetrahedra
[2026-04-04T14:25:10.201Z] [INFO]  === Adaptive Pass 2 ===
[2026-04-04T14:25:10.500Z] [INFO]  Matrix assembly: 24309 DOFs, 735000 nonzeros
[2026-04-04T14:25:19.700Z] [INFO]  Direct solve completed (9.2 sec)
[2026-04-04T14:25:20.500Z] [INFO]  Max delta S = 0.045 (target: 0.02)
[2026-04-04T14:25:21.500Z] [INFO]  Mesh refined: 8103 -> 10250 tetrahedra
[2026-04-04T14:25:21.501Z] [INFO]  === Adaptive Pass 3 ===
[2026-04-04T14:25:22.000Z] [INFO]  Matrix assembly: 30750 DOFs, 920000 nonzeros
[2026-04-04T14:25:34.500Z] [INFO]  Direct solve completed (12.5 sec)
[2026-04-04T14:25:35.600Z] [INFO]  Max delta S = 0.012 (target: 0.02) -- CONVERGED
[2026-04-04T14:25:35.601Z] [INFO]  Converged after 3 passes (2 consecutive passes below threshold)
[2026-04-04T14:25:35.602Z] [INFO]  === Frequency Sweep: Sweep1 ===
[2026-04-04T14:25:35.603Z] [INFO]  Interpolating sweep: 1.0 - 4.0 GHz, 301 points
[2026-04-04T14:27:00.600Z] [INFO]  Sweep completed: 12 actual solves, interpolation error < 0.001
[2026-04-04T14:27:00.700Z] [INFO]  Touchstone export: s_parameters.s1p
[2026-04-04T14:27:00.800Z] [INFO]  === Far Field Computation ===
[2026-04-04T14:27:06.000Z] [INFO]  InfiniteSphere1 @ 2.4GHz: peak gain = 7.2 dBi
[2026-04-04T14:27:06.001Z] [INFO]  === Solve Complete ===
[2026-04-04T14:27:06.002Z] [INFO]  Total time: 126.0 sec, Peak memory: 620 MB
```

**日志级别**：`DEBUG` | `INFO` | `WARN` | `ERROR` | `FATAL`

---

## 3. 二进制格式文件

### 3.1 设计思想

二进制文件用于存储大体量数值数据（网格、场解），设计目标：

- **紧凑**：比 JSON 节省 5~10 倍空间
- **内存映射友好**：支持 `mmap` 直接映射到结构体数组
- **随机访问**：通过文件头索引表，可按频点/分量直接跳转读取
- **跨平台**：固定使用小端（Little-Endian）字节序、IEEE 754 浮点

### 3.2 网格数据（*.msh — Gmsh MSH 4.1 格式）

采用 [Gmsh MSH 4.1](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format) 开放标准格式存储四面体有限元网格，而非自定义二进制格式。

**选择 Gmsh MSH 4.1 的理由**：

| 优势 | 说明 |
|------|------|
| **开放标准** | 由 [Gmsh](https://gmsh.info/) 定义，文档公开、格式稳定 |
| **生态成熟** | ParaView、meshio、MFEM、FEniCS、deal.II 等主流工具直接支持读写 |
| **ASCII + Binary** | 同一格式支持可读文本模式（调试）和紧凑二进制模式（生产） |
| **丰富的元素类型** | 原生支持线性/二阶四面体、三角形、六面体等，覆盖 FEM 需求 |
| **Physical Groups** | 可映射 EMStudio 的材料区域、边界标记、命名选择 |
| **免实现成本** | Rust 生态已有 `gmsh` crate 和 `nom-mesh` 等解析库 |

**EMStudio 的使用约定**：

- 默认使用 **Binary 模式**（`file-type=1`）以减小文件体积
- 调试时可切换为 **ASCII 模式**（`file-type=0`）
- 使用 `$PhysicalNames` 存储材料区域和边界标记的语义名称
- 使用 `$Entities` 关联 EMStudio 几何对象

#### 3.2.1 Physical Groups 映射约定

EMStudio 使用 Gmsh 的 Physical Groups 机制来标记网格中的材料区域和边界条件：

| Physical Group 维度 | 用途 | 命名规则 | 示例 |
|-------------------|------|---------|------|
| 3 (Volume) | 材料区域 | `mat:<material_name>` | `mat:FR4_epoxy`, `mat:copper`, `mat:vacuum` |
| 2 (Surface) | 边界条件标记 | `bc:<boundary_name>` | `bc:Radiation1`, `bc:PEC_GND`, `bc:Port1` |
| 2 (Surface) | 命名选择（面） | `ns:<selection_name>` | `ns:GND_Bottom`, `ns:FeedPort` |
| 1 (Curve) | 命名选择（边） | `ns:<selection_name>` | `ns:PatchEdge` |

#### 3.2.2 文件结构示例

**ASCII 模式**（用于说明，生产环境使用 Binary）：

```
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
5
3 1 "mat:vacuum"
3 2 "mat:FR4_epoxy"
3 3 "mat:copper"
2 4 "bc:Radiation1"
2 5 "bc:PEC_GND"
$EndPhysicalNames
$Entities
0 0 0 3
! 三个 Volume entities
3 1 0.0 0.0 0.0 120.0 120.0 62.0 1 1 0
3 2 0.0 0.0 0.035 60.0 60.0 1.635 1 2 0
3 3 0.0 0.0 0.0 60.0 60.0 0.035 1 3 0
$EndEntities
$Nodes
3 2380 1 2380
3 1 0 850
1
2
3
...
0.0 0.0 0.0
1.5 0.0 0.0
3.0 0.0 0.0
...
3 2 0 980
851
852
...
0.0 0.0 0.035
0.5 0.0 0.035
...
3 3 0 550
1831
1832
...
$EndNodes
$Elements
3 10250 1 10250
3 1 4 3500
1 1 2 3 4
2 2 3 4 5
...
3 2 4 4200
3501 851 852 853 854
...
3 3 4 2550
7701 1831 1832 1833 1834
...
$EndElements
```

**要点说明**：

| 部分 | 说明 |
|------|------|
| `$MeshFormat` | `4.1` 版本号，`0` 表示 ASCII（生产用 `1` 表示 binary），`8` 表示 `sizeof(size_t)` |
| `$PhysicalNames` | 将 physical group tag 映射到语义名称（`mat:*` / `bc:*` / `ns:*`） |
| `$Entities` | 几何实体定义，每个 volume entity 关联一个 physical group（材料） |
| `$Nodes` | 按 entity 分块存储节点，先列 tag 再列坐标 |
| `$Elements` | 按 entity 分块存储四面体（type=`4`），每行 `elementTag node1 node2 node3 node4` |

#### 3.2.3 Binary 模式细节

Binary 模式下，`$MeshFormat` 行本身仍为 ASCII，但其后紧跟一个 4 字节整数 `1`（用于字节序检测）。`$Nodes` 和 `$Elements` 的块头为 ASCII，数据部分为二进制：

- 节点坐标：`f64 × 3`（24 字节/节点）
- 节点 tag：`size_t`（8 字节/节点）
- 四面体连接：`size_t × 5`（40 字节/元素，含 element tag）

**体量对比**（10000 个四面体、~2300 个节点）：

| 模式 | 文件大小 | 说明 |
|------|---------|------|
| ASCII | ~1.2 MB | 人类可读，适合调试 |
| Binary | ~520 KB | 紧凑高效，适合生产 |
| JSON 等效 | ~2.5 MB | 对比参考 |

#### 3.2.4 可选的 $NodeData 段

Gmsh MSH 4.1 还支持 `$NodeData` 和 `$ElementData` 段，可直接在网格文件中附带场数据。EMStudio **不使用**此机制存储场解（场解单独存放在 `.emsfld` 文件中），但可利用它存储**网格质量指标**等辅助数据：

```
$NodeData
1
"element_quality"
1
0.0
3
0
1
10250
1 2.1
2 1.8
3 3.5
...
$EndNodeData
```

#### 3.2.5 Rust 读写接口

```rust
/// 网格文件读写（基于 Gmsh MSH 4.1）
pub struct MshMesh {
    pub version: String,           // "4.1"
    pub binary: bool,              // true=binary, false=ASCII
    pub physical_names: Vec<PhysicalName>,
    pub nodes: Vec<MshNode>,
    pub elements: Vec<MshElement>,
}

pub struct PhysicalName {
    pub dimension: u32,            // 1=curve, 2=surface, 3=volume
    pub tag: u32,
    pub name: String,              // "mat:copper", "bc:PEC_GND", etc.
}

pub struct MshNode {
    pub tag: u64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

pub struct MshElement {
    pub tag: u64,
    pub element_type: u32,         // 4=tet4, 11=tet10, 2=tri3, etc.
    pub node_tags: Vec<u64>,
}

impl MshMesh {
    /// 从 .msh 文件加载（自动检测 ASCII/Binary）
    pub fn load(path: &Path) -> Result<Self, MeshError> { todo!() }

    /// 保存为 .msh 文件
    pub fn save(&self, path: &Path, binary: bool) -> Result<(), MeshError> { todo!() }

    /// 按 Physical Group 名称前缀筛选元素
    pub fn elements_by_material(&self, material: &str) -> Vec<&MshElement> { todo!() }

    /// 获取边界面元素
    pub fn boundary_elements(&self, boundary: &str) -> Vec<&MshElement> { todo!() }
}
```

---

### 3.3 场解数据（*.emsfld）

存储频域 FEM 场解（复数向量场），支持多频点、多分量。

**文件布局**：

```
┌──────────────────────────────────────────────────────┐
│ Header (固定 128 字节)                                │
├──────────────────────────────────────────────────────┤
│ Frequency Table (num_frequencies × f64)               │
├──────────────────────────────────────────────────────┤
│ Field Data Index (num_frequencies × FieldBlockInfo)   │
├──────────────────────────────────────────────────────┤
│ Field Block 0 (频点 0 的场数据)                       │
├──────────────────────────────────────────────────────┤
│ Field Block 1 (频点 1 的场数据)                       │
├──────────────────────────────────────────────────────┤
│ ...                                                   │
└──────────────────────────────────────────────────────┘
```

**Header 结构（128 字节）**：

```rust
#[repr(C, packed)]
pub struct EmsFldHeader {
    pub magic: [u8; 8],           // b"EMSFLD\0\0"
    pub version: u32,             // 1
    pub byte_order: u32,          // 0x01020304
    pub field_type: u32,          // 0=E-field, 1=H-field, 2=J-field, 3=Combined
    pub data_type: u32,           // 0=complex f64, 1=complex f32
    pub num_nodes: u64,           // 场采样点数（通常等于网格节点数）
    pub num_components: u32,      // 3 (向量场 x,y,z) 或 1 (标量场)
    pub num_frequencies: u32,     // 频点数
    pub frequency_unit: u32,      // 0=Hz, 1=kHz, 2=MHz, 3=GHz
    pub freq_table_offset: u64,   // 频率表偏移量
    pub index_offset: u64,        // 索引表偏移量
    pub data_offset: u64,         // 第一个 Field Block 偏移量
    pub mesh_file: [u8; 32],      // 关联的网格文件名（如 "final_mesh.msh"）
    pub _reserved: [u8; 12],
}
```

**Field Block Info（每频点索引，16 字节）**：

```rust
#[repr(C, packed)]
pub struct FieldBlockInfo {
    pub offset: u64,              // 该频点数据块在文件中的偏移量
    pub size_bytes: u64,          // 该频点数据块大小
}
```

**Field Block 数据**（单频点）：

```
对于 complex f64 向量场 (E_x, E_y, E_z)：
  每个节点: 6 × f64 = 48 字节 (re_x, im_x, re_y, im_y, re_z, im_z)
  总大小: num_nodes × 48 字节
```

| 数据类型 | 每节点字节 | 10K 节点/1 频点 | 10K 节点/301 频点 |
|---------|-----------|----------------|------------------|
| complex f64 向量场 | 48 | 480 KB | 141 MB |
| complex f32 向量场 | 24 | 240 KB | 70 MB |

**随机访问模式**：
1. 读取 Header → 获取 `index_offset`
2. 读取 FieldBlockInfo[freq_idx] → 获取目标频点的 `offset` 和 `size`
3. Seek 到 `offset`，读取该频点数据

---

## 4. 文件间关联关系

```
convergence.json ──────────────────► mesh/pass_N_mesh.msh
  (pass_number)                       (每轮对应一个网格)

mesh_stats.json ───────────────────► mesh/final_mesh.msh
  (统计信息来源)                       (最终收敛网格)

s_parameters.json ─────────────────► solutions/final_solution.emsfld
  (从场解中提取)                       (FEM 系数)

far_field_*.json ──────────────────► solutions/final_solution.emsfld
  (从场解+积分计算)                    + mesh/final_mesh.msh

fields/e_field_*.emsfld ───────────► mesh/final_mesh.msh
  (场数据与网格配对)                    (共享节点索引)

profile.json ──────────────────────► convergence.json
  (剖面包含收敛过程耗时)                (收敛轮次对应)

summary.json (Optimetrics) ────────► variation_*/s_parameters.json
  (汇总索引指向各变量组合)              (每组变量的结果)
```

**Q3D 准静态求解文件关联**：

```
convergence.json (Q3D) ───────────► mesh/pass_N_mesh.msh
  (delta_energy + rlcg_snapshot)      (MoM 面网格)

rlcg_matrix.json ─────────────────► convergence.json
  (RLCG 矩阵随频率变化)                (收敛历史中的 rlcg_snapshot)

rlcg_matrix.json ─────────────────► equivalent_circuit.sp
  (矩阵数据)                          (导出为 SPICE 网表)

fields/j_field_*.emsfld ──────────► mesh/final_mesh.msh
  (Q3D 电流密度场)                     (MoM 面网格节点索引)

fields/e_field_*.emsfld ──────────► mesh/final_mesh.msh
  (Q3D 电场数据)                       (用于电容/电导提取)

rlcg_matrix.json ─────────────────► s_parameters.snp (可选)
  (RLCG → S 参数转换)                  (Touchstone 导出)
```

---

## 5. Rust 类型定义

```rust
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// ========================
// 通用文件头
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileHeader {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
}

// ========================
// 验证报告
// ========================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CheckStatus {
    Pass,
    Warning,
    Error,
    Info,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationCheck {
    pub category: String,
    pub check: String,
    pub status: CheckStatus,
    pub message: String,
    pub details: Option<serde_json::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationReport {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
    pub validated_at: String,
    pub overall_status: CheckStatus,
    pub checks: Vec<ValidationCheck>,
    pub errors: u32,
    pub warnings: u32,
    pub info: u32,
}

// ========================
// 收敛历史
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PassMeshInfo {
    pub num_tetrahedra: u64,
    pub num_nodes: u64,
    pub min_edge_length_mm: f64,
    pub max_edge_length_mm: f64,
    pub mean_edge_length_mm: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PassSolutionInfo {
    pub max_delta_s: Option<f64>,
    pub delta_energy: Option<f64>,
    pub matrix_size: u64,
    pub num_rhs: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PassPerformance {
    pub mesh_time_sec: f64,
    pub solve_time_sec: f64,
    pub error_estimation_time_sec: f64,
    pub total_time_sec: f64,
    pub peak_memory_mb: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdaptivePass {
    pub pass_number: u32,
    pub timestamp: String,
    pub mesh: PassMeshInfo,
    pub solution: PassSolutionInfo,
    pub performance: PassPerformance,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceResult {
    pub converged: bool,
    pub converged_at_pass: Option<u32>,
    pub consecutive_converged_passes: u32,
    pub total_elapsed_sec: f64,
    pub final_max_delta_s: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceHistory {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
    pub setup: String,
    pub solution_frequency: String,
    pub target_max_delta_s: f64,
    pub max_passes_limit: u32,
    pub passes: Vec<AdaptivePass>,
    pub result: Option<ConvergenceResult>,
}

// ========================
// S 参数
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SParameterData {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
    pub setup: String,
    pub sweep: Option<String>,
    pub reference_impedance_ohm: f64,
    pub num_ports: u32,
    pub port_names: Vec<String>,
    pub num_frequencies: u32,
    pub frequency_unit: String,
    pub data_format: String,
    pub frequencies: Vec<f64>,
    pub parameters: HashMap<String, ComplexArray>,
    pub derived: Option<HashMap<String, Vec<f64>>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComplexArray {
    pub real: Vec<f64>,
    pub imag: Vec<f64>,
}

// ========================
// Q3D RLCG 矩阵
// ========================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RlcgMatrixData {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
    pub setup: String,
    pub solution_type: String,
    pub num_nets: u32,
    pub net_names: Vec<String>,
    pub num_terminals: u32,
    pub terminal_names: Vec<String>,
    pub ground_net: String,
    pub num_frequencies: u32,
    pub frequency_unit: String,
    pub frequencies: Vec<f64>,
    pub matrices: RlcgMatrices,
    pub dc_data: Option<DcMatrixData>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RlcgMatrices {
    pub r: Option<MatrixSeries>,
    pub l: Option<MatrixSeries>,
    pub c: Option<MatrixSeries>,
    pub g: Option<MatrixSeries>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MatrixSeries {
    pub description: String,
    pub unit: String,
    pub data_per_frequency: Vec<FrequencyMatrix>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrequencyMatrix {
    pub frequency: f64,
    pub matrix: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DcMatrixData {
    pub r_dc: Option<DcMatrix>,
    pub l_dc: Option<DcMatrix>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DcMatrix {
    pub description: String,
    pub unit: String,
    pub matrix: Vec<Vec<f64>>,
}

/// Q3D 收敛历史中每轮的 RLCG 矩阵摘要
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RlcgSnapshot {
    pub r_max_ohm: Option<f64>,
    pub l_max_nh: Option<f64>,
    pub c_max_pf: Option<f64>,
    pub g_max_ms: Option<f64>,
}

/// 等效电路导出配置
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EquivalentCircuitExport {
    pub format_version: String,
    pub file_type: String,
    pub design_id: String,
    pub setup: String,
    pub export_format: String,
    pub model_type: EquivalentCircuitModelType,
    pub frequency_range: FrequencyRange,
    pub nets_included: Vec<String>,
    pub ground_net: String,
    pub output_file: String,
    pub options: CircuitExportOptions,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum EquivalentCircuitModelType {
    BroadbandLumped,
    FrequencyDependentLumped,
    TLineModel,
    SParameterBlock,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CircuitExportOptions {
    pub include_dc_resistance: bool,
    pub include_mutual_inductance: bool,
    pub include_mutual_capacitance: bool,
    pub include_dielectric_loss: bool,
    pub coupling_threshold: f64,
}

// ========================
// 二进制文件头（内存映射用）
// ========================

// ========================
// 网格文件（Gmsh MSH 4.1）
// 完整定义见 §3.2.5
// ========================

// MshMesh, PhysicalName, MshNode, MshElement 类型定义见 §3.2.5

#[repr(C, packed)]
pub struct EmsFldHeader {
    pub magic: [u8; 8],
    pub version: u32,
    pub byte_order: u32,
    pub field_type: u32,
    pub data_type: u32,
    pub num_nodes: u64,
    pub num_components: u32,
    pub num_frequencies: u32,
    pub frequency_unit: u32,
    pub freq_table_offset: u64,
    pub index_offset: u64,
    pub data_offset: u64,
    pub mesh_file: [u8; 32],
    pub _reserved: [u8; 12],
}

#[repr(C, packed)]
pub struct FieldBlockInfo {
    pub offset: u64,
    pub size_bytes: u64,
}
```

---

## 6. 参考资料

- [Gmsh MSH 4.1 File Format](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format) — 网格文件标准格式（本项目采用）
- [Gmsh 官方文档](https://gmsh.info/doc/texinfo/) — Gmsh 完整参考手册
- [Touchstone 2.0 Specification](https://ibis.org/touchstone_ver2.0/touchstone_ver2_0.pdf) — S 参数文件标准
- [Touchstone 2.1](https://www.ibis.org/touchstone_ver2.1/) — 最新版 Touchstone 规范
- [MFEM Mesh Format v1.0](https://mfem.org/mesh-format-v1.0/) — FEM 网格文件格式参考
- [VTK File Formats](https://docs.vtk.org/en/latest/vtk_file_formats/vtk_legacy_file_format.html) — 可视化工具包文件格式
- [HDF5 Data Model](https://docs.hdfgroup.org/documentation/hdf5/latest/_h5_d_m__u_g.html) — 科学数据存储格式
- [Amelet-HDF](https://www.axessim.fr/docs/amelethdf/1.5.4/) — 电磁仿真数据交换 HDF5 模式
- [HFSS Exporting Field Plots](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/ReportsandPostProc/ExportingFieldPlots.htm) — HFSS 场导出
- [HFSS Mesh Statistics](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/ReportsandPostProc/ViewingMeshStatistics.htm) — HFSS 网格统计
- [Q3D Viewing Matrix Data](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/ViewingMatrixDatainQ3D.htm) — Q3D 矩阵数据查看
- [Q3D Exporting Equivalent Circuit Data](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/ExportingQ3DExtractorEquivalentCircuitData.htm) — Q3D 等效电路导出
- [Q3D Exporting S-Parameters](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Q3DExtractor/Content/Q3D/ExportingSParameterData.htm) — Q3D S 参数导出
- [Q3D Frequency Sweeps](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/FrequencySweepsinQ3DExtractor.htm) — Q3D 频率扫描
- [Q3D Plotting Field Overlays](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Q3DExtractor/Content/ReportsandPostProc/PlottingFieldOverlaysinQ3D.htm) — Q3D 场叠加显示
- [PyAEDT Q3D Class Reference](https://aedt.docs.pyansys.com/version/stable/API/_autosummary/ansys.aedt.core.q3d.Q3d.html) — PyAEDT Q3D API
- [HFSS Solution Profile](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Mechanical/Content/ReportsandPostProc/ViewingaSolutionProfile.htm) — HFSS 求解剖面
