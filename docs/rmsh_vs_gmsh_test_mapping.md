# RMSH 测试用例与 GMSH 对照表

## 概述

rmsh (Rust Mesh) 是一个 Rust 网格生成与优化框架，其算法体系、文件格式、Python API 和拓扑分类均围绕 GMSH 生态设计。本文档以 **GMSH 功能维度** 为线索，逐项映射对应的 rmsh 测试用例。

| GMSH 功能域 | rmsh 模块 | 对照关系 |
|---|---|---|
| 2D 网格算法 | `crates/algo/src/` | 算法号 1:1 映射 |
| 3D 网格算法 | `crates/algo/src/` | 算法号 1:1 映射 |
| 网格优化 | `crates/algo/src/` | 算法移植自 GMSH 源码 |
| MSH 文件格式 | `crates/io/src/msh.rs` | v2.2/v4.1 全读写 |
| STEP I/O | `crates/io/src/step.rs` | gmsh_strict 兼容模式 |
| Python API | `crates/py/src/lib.rs` | API 签名 1:1 模仿 |
| 单元类型 ID | `crates/model/src/element.rs` | `from_gmsh_type_id()` |
| 拓扑分类 | `crates/geo/src/classify.rs` | 移植自 `classifyFaces` |
| 渲染风格 | `crates/renderer/src/scene.rs` | GMSH 默认配色 |

---

## 一、2D 网格算法 (GMSH Mesh.Algorithm)

| GMSH 算法号 | GMSH 名称 | rmsh 模块 | rmsh 测试函数 | 测试内容 |
|---|---|---|---|---|
| 1 | MeshAdapt | `mesh_adapt_2d.rs` | `mesh_adapt_handles_square_with_hole` | 带孔方形域的局部边分裂/折叠/交换 |
| 5 | Delaunay | `delaunay_2d.rs` | `delaunay_2d_meshes_square` | 方形域 Bowyer-Watson 三角剖分 |
| 5 | Delaunay | `delaunay_2d.rs` | `delaunay_2d_meshes_with_hole` | 带孔域边界约束三角剖分 |
| 5 | Delaunay | `delaunay_2d.rs` | `delaunay_2d_name_is_stable` | 模块标识 |
| 6 | Frontal-Delaunay | `frontal_delaunay_2d.rs` | `frontal_delaunay_handles_l_shape` | L 形域推进波前 |
| 6 | Frontal-Delaunay | `frontal_delaunay_2d.rs` | `frontal_delaunay_handles_rectangle` | 矩形域推进波前 |
| 6 | Frontal-Delaunay | `frontal_delaunay_2d.rs` | `frontal_delaunay_handles_hole_domain` | 带孔方域推进波前 |
| 6 | Frontal-Delaunay | `frontal_delaunay_2d.rs` | `bowyer_watson_insertion_adds_local_triangles` | 逐点插入验证 |
| 6 | Frontal-Delaunay | `frontal_delaunay_2d.rs` | `segment_intersection_detects_crossing` | 线段相交检测 |
| 6 | Frontal-Delaunay | `frontal_delaunay_2d.rs` | `frontal_quality_stays_close_to_planar_fallback` | 网格质量退化兜底 |
| 7 | BAMG | `bamg_2d.rs` | `bamg_metric_affects_density` | 各向异性度量控制密度 (自适应循环) |
| 7 | BAMG | `bamg_2d.rs` | `bamg_adaptive_loop_converges` | 方形域自适应 split/swap/smooth 收敛 |
| 7 | BAMG | `bamg_2d.rs` | `bamg_anisotropic_stretches` | 各向异性度量 (1.0, 0.1, 0°) 拉伸验证 |
| 7 | BAMG | `bamg_2d.rs` | `bamg_metric_swap_preserves_orientation` | 度量驱动的边翻转保向性 |
| 7 | BAMG | `bamg_2d.rs` | `metric2_intersect_isotropic` / `_anisotropic` | 度量张量相交插值 |
| 7 | BAMG | `bamg_2d.rs` | `edge_metric_length_isotropic` | 度量空间边长 (2 点 Gauss-Legendre) |
| 7 | BAMG | `bamg_2d.rs` | `metric_midpoint_near_endpoints` | 度量空间中点计算 |
| 8 | Frontal-Quads | `frontal_quads_2d.rs` | `frontal_quads_rectangle` | 矩形域三角形→四边形重组 |
| 8 | Frontal-Quads | `frontal_quads_2d.rs` | `frontal_quads_square` | 方形域四边形网格生成 |
| 9/11 | Quad Paving | `quad_paving_2d.rs` | `quad_paving_rectangle_produces_quads` | 矩形域四边形生成 |

**GMSH 参考源码**: `meshGFaceDelaunayInsertion.cpp`, `Mesh/meshGFaceQuadqs.cpp`

---

## 二、3D 网格算法 (GMSH Mesh.Algorithm3D)

| GMSH 算法号 | GMSH 名称 | rmsh 模块 | rmsh 测试函数 | 测试内容 |
|---|---|---|---|---|
| 1 | Delaunay | `delaunay_3d.rs` | `delaunay3d_name_is_stable` | 模块标识 |
| 1 | Delaunay | `delaunay_3d.rs` | `delaunay3d_mesh_flow_runs` | 立方体完整网格流水线 |
| 1 | Delaunay | `delaunay_3d.rs` | `delaunay3d_respects_mesh_size_density` | 网格尺寸密度响应 |
| 1 | Delaunay | `delaunay_3d.rs` | `delaunay3d_rejects_bad_mesh_params` | 非法参数拒绝 |
| 1 | Delaunay | `delaunay_3d.rs` | `delaunay3d_rejects_invalid_algo_params` | 非法算法参数拒绝 |
| 1 | Delaunay | `delaunay_3d.rs` | `validate_params_accepts_default_configuration` | 默认配置验证 |
| 1 | Delaunay | `delaunay_3d.rs` | `circumsphere_and_in_sphere_work_for_regular_tet` | 外接球/球内测试 |
| 1 | Delaunay | `delaunay_3d.rs` | `radius_edge_ratio_is_finite_*` / `*_known_value` | 半径边比计算 |
| 1 | Delaunay | `delaunay_3d.rs` | `local_flip_pass_activates_on_edge_fan` | 边扇局部翻转 |
| 1 | Delaunay | `delaunay_3d.rs` | `refinement_produces_more_elements_than_seed` | Delaunay 细化产出更多单元 |
| 1 | Delaunay | `delaunay_3d.rs` | 3 个 `solve_3x3_*` | 3x3 线性系统求解器 |
| 1 | Delaunay | `delaunay_3d.rs` | `tetra_centroid_is_average_of_four_nodes` | 四面体质心 |
| 4 | Frontal 3D | `frontal_3d.rs` | `frontal_3d_generates_mesh` | 立方体推进波前 3D |
| 10 | HXT | `hxt_3d.rs` | `hilbert_index_progresses_along_diagonal` | Hilbert 排序索引 |
| 10 | HXT | `hxt_3d.rs` | `hilbert_index_distinguishes_adjacent_points` | Hilbert 区分邻近点 |
| 10 | HXT | `hxt_3d.rs` | `hilbert_index_is_deterministic` | Hilbert 确定性 |
| 10 | HXT | `hxt_3d.rs` | `hilbert_index_out_of_range_clamps` | Hilbert 边界截断 |
| 10 | HXT | `hxt_3d.rs` | `grid_coloring_eight_colors` | 3D 8色网格着色 |
| 10 | HXT | `hxt_3d.rs` | `adjacent_cells_have_different_colors` | 相邻网格颜色不同 |
| 10 | HXT | `hxt_3d.rs` | `tet_ownership_exclusive_access` | CAS 所有权互斥 |
| 10 | HXT | `hxt_3d.rs` | `split_containing_tet_works` | 四面体 1→4 分裂 |
| 10 | HXT | `hxt_3d.rs` | `split_containing_tet_outside_point` | 外部点无分裂 |
| 10 | HXT | `hxt_3d.rs` | `hxt_3d_generates_mesh` | 立方体 HXT 并行网格 |
| 10 | HXT | `hxt_3d.rs` | `hxt_3d_single_threaded_works` | 单线程模式 |
| 7 | MMG3D | `mmg_remesh.rs` | `classify_edges_buckets_metric_lengths` | 边长分类 |
| 7 | MMG3D | `mmg_remesh.rs` | `mmg_remesh_generates_volume_mesh` | 立方体 MMG3D 网格 |

**GMSH 参考源码**: `Mesh/meshGRegion.cpp`

---

## 三、基础 3D 网格生成（Centroid Star / 兜底方案）

非独立 GMSH 算法，作为高级算法的 Seed 或退化兜底。

| rmsh 模块 | rmsh 测试函数 | 测试内容 |
|---|---|---|
| `tetrahedralize3d.rs` | `tetrahedralize_cube_surface` | 立方体闭面 → 12 四面体 |
| `tetrahedralize3d.rs` | `centroid_star_mesher_trait_path_works` | Trait 接口路径 |
| `tetrahedralize3d.rs` | `centroid_star_mesher_rejects_invalid_params` | 非法参数拒绝 |
| `tetrahedralize3d.rs` | `tetrahedralize_rejects_empty_input` | 空输入拒绝 |
| `tetrahedralize3d.rs` | `tetrahedralize_reports_missing_nodes_*` | 缺失节点报告 |
| `tetrahedralize3d.rs` | `tetrahedralize_rejects_degenerate_generated_tets` | 退化四面体拒绝 |
| `tetrahedralize3d.rs` | `collect_boundary_polygons_*` (2 个) | 边界多边形收集 |
| `tetrahedralize3d.rs` | `centroid_of_nodes_*` (2 个) | 节点质心计算 |
| `tetrahedralize3d.rs` | `tetra_signed_volume6_*` (2 个) | 四面体有符号体积 |
| `tetrahedralize3d.rs` | `output_volume_sum_equals_cube_volume` | 体积和 = 1.0 |
| `tetrahedralize3d.rs` | `every_node_used_in_at_least_one_tet` | 节点无孤立 |
| `tetrahedralize3d.rs` | `no_degenerate_tets_in_cube_output` | 无退化四面体 |
| `tetrahedralize3d.rs` | `element_count_is_correct_for_cube` | 单元计数验证 |

---

## 四、Bowyer-Watson 2D Delaunay（基础三角剖分）

| rmsh 模块 | rmsh 测试函数 | 对应 GMSH 概念 |
|---|---|---|
| `triangulate2d.rs` | `unit_square_delaunay` | GMSH 基础 Delaunay 三角剖分 |
| `triangulate2d.rs` | `mesh_unit_square` | `gmsh.model.mesh.generate(2)` |
| `triangulate2d.rs` | `mesh_l_shape` | 非凸域三角剖分 |
| `triangulate2d.rs` | `mesh_polygon_produces_planar_2d_mesh` | 平面 2D 网格约束 |
| `triangulate2d.rs` | `mesh_rejects_bad_inputs` | 非法输入拒绝 |
| `triangulate2d.rs` | `point_in_polygon` | GMSH `getBoundary` 点包含测试 |
| `triangulate2d.rs` | `fewer_than_three_points_returns_empty` | 退化输入 |
| `triangulate2d.rs` | `exactly_three_points_gives_one_triangle` | 极小输入 |
| `triangulate2d.rs` | `all_output_indices_are_valid` | 索引有效性 |
| `triangulate2d.rs` | `no_degenerate_triangles_in_output` | 无零面积三角形 |
| `triangulate2d.rs` | `no_duplicate_triangles_in_output` | 无重复三角形 |
| `triangulate2d.rs` | `all_input_points_appear_in_output_*` | 所有输入点都被使用 |
| `triangulate2d.rs` | `finer_mesh_size_produces_more_elements` | 网格密度响应（对应 `Mesh.MeshSizeFactor`） |
| `triangulate2d.rs` | `mesh_polygon_all_centroids_inside_polygon` | 单元包含性验证 |

---

## 五、网格优化

### Laplacian 平滑

| rmsh 模块 | rmsh 测试函数 | 对应 GMSH |
|---|---|---|
| `laplacian_smooth.rs` | `smooth_does_not_change_node_count` | `Mesh.Smoothing` = 1 |
| `laplacian_smooth.rs` | `smooth_does_not_change_element_count` | `Mesh.Smoothing` = 1 |
| `laplacian_smooth.rs` | `smooth_does_not_move_boundary_nodes_by_default` | 边界节点不变 |
| `laplacian_smooth.rs` | `smooth_returns_ok_on_empty_mesh` | 空网格稳定 |
| `laplacian_smooth.rs` | `uniform_smooth_on_regular_grid_converges` | 规则网格收敛 |
| `laplacian_smooth.rs` | `uniform_smooth_with_omega_zero_point_five_keeps_nodes_inside_domain` | 松弛因子 ω=0.5 |
| `laplacian_smooth.rs` | `cotangent_variant_smooths_triangle_mesh` | Cotangent 加权变体（已实现） |
| `laplacian_smooth.rs` | `taubin_variant_smooths_without_error` | Taubin λ/μ 变体（已实现） |
| `laplacian_smooth.rs` | `build_node_adjacency_two_triangles` | 节点邻接表构建 |
| `laplacian_smooth.rs` | `build_edge_triangle_map` | 边→三角形映射表 |
| `laplacian_smooth.rs` | `collect_boundary_nodes_two_triangles` | 边界节点收集 |

**GMSH 参考源码**: `Mesh/qualityMeasures.cpp`

### TetMesh 翻转

| rmsh 模块 | rmsh 测试函数 | 对应 GMSH |
|---|---|---|
| `tet_mesh.rs` | `round_trip_conversion` | TetMesh ↔ Mesh 转换 |
| `tet_mesh.rs` | `neighbor_table_two_tets` / `three_tets` | 邻居表构建 |
| `tet_mesh.rs` | `apply_2to3_produces_three_tets` | 2→3 翻转 |
| `tet_mesh.rs` | `apply_3to2_produces_two_tets` | 3→2 翻转 |
| `tet_mesh.rs` | `bistellar_flip_4_to_4_works` | 4→4 翻转 (4 四面体共享边) |
| `tet_mesh.rs` | `bistellar_flip_4_to_4_rejects_wrong_count` | 4→4 翻转参数校验 |
| `tet_mesh.rs` | `tetmesh_flip_activates_on_edge_fan` | 边扇翻转激活 |
| `tet_mesh.rs` | `quality_parity_with_mesh_version` | 聚合质量指标 |
| `tet_mesh.rs` | `quality_parity_with_mesh_version` | min_dihedral / sliver_fraction / max_radius_edge |

### 网格质量优化器（已实现）

| rmsh 模块 | rmsh 测试函数 | 对应 GMSH |
|---|---|---|
| `mesh_optimize.rs` | `test_triangle_min_angle_equilateral` | 三角形最小角指标 |
| `mesh_optimize.rs` | `test_triangle_scaled_jacobian_equilateral` | 三角形 Scaled Jacobian |
| `mesh_optimize.rs` | `test_triangle_aspect_ratio_unit` / `_degenerate` | 三角形纵横比 |
| `mesh_optimize.rs` | `test_tet_min_dihedral_regular` | 四面体最小二面角 |
| `mesh_optimize.rs` | `test_radius_edge_ratio_regular` | 四面体半径边比 |
| `mesh_optimize.rs` | `test_tet_scaled_jacobian_regular` | 四面体 Scaled Jacobian |
| `mesh_optimize.rs` | `test_tet_aspect_ratio_regular` | 四面体纵横比 |
| `mesh_optimize.rs` | `test_should_swap_2d_improves_quality` | 2D 边翻转 |
| `mesh_optimize.rs` | `test_node_insertion_split_triangle` | 节点插入（形心分裂） |
| `mesh_optimize.rs` | `test_edge_collapse_short_edge` | 短边折叠 |
| `mesh_optimize.rs` | `test_optimizer_on_empty_mesh` | 空网格稳定 |
| `mesh_optimize.rs` | `test_optimizer_triangle_quality` | 三角形质量优化 |
| `mesh_optimize.rs` | `test_quality_improvement` | 质量提升验证 |

**GMSH 参考源码**: `Mesh/qualityMeasures.cpp`, `Mesh/meshGRegionDelaunayInsertion.cpp`

### 质量回归测试

| rmsh 测试文件 | rmsh 测试函数 | 对应 GMSH 场景 |
|---|---|---|
| `tests/mesher3d_quality_regression.rs` | `mesher3d_quality_baseline_cube` | 立方体 Delaunay/Frontal/HXT 质量基线 |
| `tests/mesher3d_quality_regression.rs` | `mesher3d_quality_baseline_stretched_box` | 拉伸盒体 Frontal ≥ Delaunay 质量 |
| `tests/mesher3d_quality_regression.rs` | `mesher3d_quality_slender_box_edge_pressure` | 细长盒体边缘压力测试 |

---

## 六、MSH 文件 I/O（GMSH 原生格式）

| rmsh 模块 | rmsh 测试函数 | GMSH 格式 |
|---|---|---|
| `crates/io/src/msh.rs` | `test_parse_simple_msh_v4` | MSH v4.1 ASCII 解析 |
| `crates/io/src/msh.rs` | `test_parse_simple_msh_v2` | MSH v2.2 ASCII 解析 |
| `crates/io/src/msh.rs` | `test_parse_msh_v4_sets_entity_as_physical_tag` | MSH v4 entity tag → physical |
| `crates/io/src/msh.rs` | `test_load_msh_from_bytes_binary_v2` | MSH v2.2 二进制解析 |
| `crates/io/src/msh.rs` | `test_load_msh_from_bytes_binary_v2_first_block_type_10` | 二进制 v2 高阶三角形 (type 10) |
| `crates/io/src/msh.rs` | `test_load_msh_from_bytes_binary_v2_type_10_then_type_12` | 二进制 v2 多 element block |
| `crates/io/src/msh.rs` | `test_load_msh_from_bytes_binary_v2_preserves_unknown_high_order_elements` | 二进制 v2 未知高阶类型 (type 118) |
| `crates/io/src/msh.rs` | `test_write_roundtrip_msh_v2` | MSH v2 写出+回读 |
| `crates/io/src/msh.rs` | `test_write_roundtrip_msh_v4` | MSH v4 写出+回读 |
| `crates/io/src/msh.rs` | `test_save_msh_v4_to_path_and_load_msh_from_path_roundtrip` | MSH v4 文件 I/O 回环 |
| `crates/io/src/msh.rs` | `test_load_msh_from_path_reports_io_error_for_missing_file` | 缺失文件报错 |

**对应 GMSH**: `gmsh.model.mesh.write("file.msh")` / `gmsh.open("file.msh")`

---

## 七、STEP I/O 与 GMSH 严格模式

| rmsh 测试函数 | 对应 GMSH 功能 | 验证内容 |
|---|---|---|
| `step.rs:parse_simple_tetra_faceted_brep` | GMSH 导入 STEP | 四面体 B-Rep 解析 |
| `step.rs:roundtrip_write_then_parse` | GMSH 导出 STEP | 写出+回解析 |
| `step.rs:write_empty_mesh_fails` | — | 空网格错误处理 |
| `step.rs:parse_generated_step_test_file` | GMSH 读取 STEP 文件 | 真实 STEP 文件解析 |
| `step.rs:save_and_reload_step_file` | GMSH 文件 I/O | 临时文件回环 |
| `step.rs:write_brep_default_protocol_is_ap214` | GMSH `Mesh.StepProtocol=AP214` | 默认 AP214 协议 |
| `step.rs:write_brep_ap242_with_color_emits_style_chain` | GMSH STEP 颜色导出 | AP242 颜色实体 |
| `step.rs:strict_selection_keeps_face_edges_and_curved_standalone_edges_only` | `STEP.GmshStrict=1` | 严格模式边选择 |
| `step.rs:strict_normalize_fills_missing_line_curve_on_face_edge` | `STEP.GmshStrict=1` | 严格模式补缺失曲线 |
| `step.rs:strict_normalize_drops_degenerate_no_curve_wire_edges` | `STEP.GmshStrict=1` | 严格模式去退化边 |

**Python 端严格模式测试** (对应 GMSH `gmsh.model.occ.importShapes` / `gmsh.write`):

| rmsh 测试函数 | 对比 GMSH 脚本 | GMSH 实体验证 |
|---|---|---|
| `py::lib.rs:strict_cylinder_emits_cylindrical_surface_and_seam_curve` | `export_each_3d_step_gmsh.py` | 3 ADVANCED_FACE, 1 CYLINDRICAL_SURFACE, ≥1 SEAM_CURVE |
| `py::lib.rs:strict_frustum_cone_emits_conical_side_and_three_faces` | `export_each_3d_step_gmsh.py` | 3 ADVANCED_FACE, 1 CONICAL_SURFACE, 3 EDGE_CURVE, 无三角化 |
| `py::lib.rs:strict_standalone_line_emits_wireframe_curve_set` | `export_each_3d_step_gmsh.py` | GEOMETRIC_CURVE_SET |
| `py::lib.rs:strict_standalone_spline_emits_bspline_curve` | `export_each_3d_step_gmsh.py` | B_SPLINE_CURVE_WITH_KNOTS |
| `py::lib.rs:step_protocol_uses_default_option_value` | GMSH `Mesh.StepProtocol` | 默认 AP214 |
| `py::lib.rs:option_default_number_accepts_geometry_lines_alias` | GMSH `Geometry.Lines` | 选项别名兼容 |
| `py::lib.rs:estimate_mesh_characteristic_size_uses_bbox_diagonal` | GMSH `Mesh.MeshSizeFactor` | 包围盒对角线/20 |

**Python 端额外对比用例** (脚本级, 非单元测试):

| 脚本路径 | 对比目标 | 用途 |
|---|---|---|
| `crates/py/examples/gmsh_step_baseline.py` | `gmsh` Python 包 | 生成 GMSH STEP 基线 |
| `crates/py/examples/compare_rmsh_gmsh_strict_step_entities.py` | rmsh vs gmsh | STEP 实体计数对比 |
| `crates/py/examples/export_each_3d_step_gmsh.py` | GMSH 导出 | 逐个 3D 基元 GMSH 基线 |
| `crates/py/examples/export_each_3d_step_rmsh_gmsh_strict.py` | rmsh vs gmsh | 逐个 3D 基元严格模式对照 |
| `crates/py/examples/rmsh_to_gmsh_viewer_pipeline.py` | rmsh → GMSH rewrite | 流水线兼容性 |
| `crates/py/examples/compare_boolean_api_contract_rmsh_gmsh.py` | rmsh vs gmsh | 布尔 API 签名一致性 |
| `crates/py/examples/compare_boolean_ops_rmsh_gmsh_step_entities.py` | rmsh vs gmsh | 布尔 STEP 实体计数 |
| `crates/py/examples/compare_cut_min_rmsh_gmsh_step_entities.py` | rmsh vs gmsh | Cut-min 操作 STEP 对比 |
| `crates/py/examples/boolean_ops_to_step_gmsh.py` | GMSH 布尔 + STEP | GMSH 布尔基线 |
| `crates/py/examples/boolean_ops_to_step.py` | rmsh 布尔 + STEP | rmsh 布尔对照 |
| `crates/py/examples/meshing_algorithms.py` | GMSH 算法全覆盖 | 2D/3D 算法 1/4/5/6/7/8/9/10 |
| `crates/py/examples/tutorials.py` | GMSH 官方教程 t1-t17 | API 签名逐项对应 |

---

## 八、单元类型映射 (GMSH Element Type ID)

| rmsh 测试函数 | 测试内容 | GMSH ID 范围 |
|---|---|---|
| `element.rs:unknown_high_order_types_have_correct_dimension` | 高阶类型维度 (ID 8-14) | ID 8-14 |
| `element.rs:unknown_volume_families_expose_canonical_faces_and_edges` | 体单元面/边规范 | tet/hex/prism/pyramid 族 |
| `element.rs:all_element_types_have_correct_dimension` | 全类型维度穷举 | ID 1-7, 15 + 高阶 |
| `element.rs:element_dimension_aligns_with_node_count` | 节点数→维度一致性 | 全部 |

**GMSH Element Type 映射表** (摘自 `element.rs:from_gmsh_type_id`):

| GMSH ID | GMSH 名称 | rmsh ElementType | 维度 |
|---|---|---|---|
| 15 | POINT | Point1 | 0 |
| 1 | LINE | Line2 | 1 |
| 2 | TRIANGLE | Triangle3 | 2 |
| 3 | QUAD | Quad4 | 2 |
| 4 | TETRAHEDRON | Tetrahedron4 | 3 |
| 5 | HEXAHEDRON | Hexahedron8 | 3 |
| 6 | PRISM | Prism6 | 3 |
| 7 | PYRAMID | Pyramid5 | 3 |
| 8-14, 16+ | 高阶类型 | Unknown(id) | 按 family 推断 |

---

## 九、拓扑分类 (GMSH classifyFaces)

| rmsh 模块 | rmsh 测试函数 | 对应 GMSH 源码 |
|---|---|---|
| `geo/src/classify.rs` | `volumes_connected_by_full_face_form_one_topovolume` | `classifyFaces` 全面邻接合并 |
| `geo/src/classify.rs` | `volumes_touching_at_single_node_stay_separate` | `classifyFaces` 点接触分离 |
| `geo/src/classify.rs` | `pure_line_and_point_mesh_generates_edges_and_vertices` | `classifyFaces` 1D/0D 分类 |
| `model/src/topology.rs` | `geometric_entities_have_correct_dimensions` | GMSH GVertex/GEdge/GFace/GRegion |
| `model/src/topology.rs` | `geometric_model_contains_all_dimensions` | GMSH GModel 多维度包含 |
| `model/src/topology.rs` | `dimension_containment_relationships` | GMSH 拓扑实体结构 |

**GMSH 参考**: `classifyFaces` 二面角分类算法

---

## 十、端到端流水线 (STEP → 网格 → MSH → GMSH 回环)

| rmsh 测试函数 | 流水线步骤 | 验证 GMSH 兼容性 |
|---|---|---|
| `tetrahedralize3d.rs:step_mesh_can_be_saved_as_gmsh_v2_and_v4_after_3d_meshing` | STEP → CentroidStar3D → MSH v2/v4 → 回读 | 节点/单元计数一致 |
| `viewer/src/app.rs:viewer_step_to_3d_meshing_gmsh_roundtrip` | STEP → 拓扑分类 → 3D网格 → MSH v2/v4 → 回读 | 节点/单元计数一致 |
| `viewer/src/app.rs:viewer_step_to_3d_meshing_via_centroid_star_trait_roundtrip` | STEP → CentroidStar3D → MSH v4 → 回读 | 节点/单元计数一致 |
| `viewer/src/app.rs:viewer_step_to_3d_meshing_via_delaunay_trait_roundtrip` | STEP → Delaunay3D → MSH v2 → 回读 | 节点/单元计数一致 |
| `viewer/src/app.rs:viewer_step_to_3d_meshing_via_delaunay_respects_size` | STEP → Delaunay3D 两档尺寸 | 粗/细网格区分 |

---

## 十一、选项系统 (GMSH Option Keys)

rmsh 实现了 GMSH 风格的选项键系统，当前对齐状态记录于 [GMSH_OPTIONS_REFERENCE.md](../GMSH_OPTIONS_REFERENCE.md)。

| GMSH 选项键 | rmsh 状态 | 测试覆盖 |
|---|---|---|
| `Mesh.Algorithm` | **Critical** | `meshing_algorithms.py` 全覆盖 1/5/6/7/8/9 |
| `Mesh.Algorithm3D` | **Critical** | `meshing_algorithms.py` 全覆盖 1/4/10 |
| `Mesh.MeshSizeFactor` | **Critical** | 各算法 element_size 参数测试 |
| `Mesh.Smoothing` | **High** | `laplacian_smooth.rs` 14 个测试 |
| `Geometry.Tolerance` | **Critical** | 几何容差 |
| `STEP.GmshStrict` | **rmsh 扩展** | `step.rs` 严格模式测试 |
| `STEP.Protocol` | **rmsh 扩展** | AP214 / AP242 协议选择 |

---

## 十二、渲染与可视化

| rmsh 组件 | GMSH 对应 | 细节 |
|---|---|---|
| `renderer/src/scene.rs` - 默认颜色 | GMSH 浅钢蓝表面 + 橙色节点 | 配色 1:1 复制 |
| `renderer/src/scene.rs` - 渐变背景 | GMSH 深色主题 | 渐变背景风格 |
| `viewer/src/app.rs` - 比例尺 | GMSH 标尺绘制 | 算法移植 |

---

## 汇总统计

| 类别 | rmsh 测试函数数 | 对应 GMSH 概念数 |
|---|---|---|
| 2D 网格算法 | 21 | 6 (algo 1/5/6/7/8/9) |
| 3D 网格算法 | 31 | 5 (algo 1/4/7/10 + bistellar flips) |
| Centroid Star (兜底) | 18 | — |
| Bowyer-Watson 2D | 14 | 1 (Delaunay) |
| Laplacian 平滑 | 15 | 1 (Smooth) |
| TetMesh 翻转 | 11 | 1 (Bistellar flips: 2→3, 3→2, 4→4) |
| 网格质量优化器 | 14 | 1 (Quality度量 + 拓扑优化) |
| 质量回归测试 | 3 | 1 (多算法质量基线) |
| MSH 文件 I/O | 11 | 2 (v2.2 / v4.1) |
| STEP I/O + 严格模式 | 12 | 1 (gmsh_strict) |
| 单元类型映射 | 4 | 8 (GMSH type ID 族) |
| 拓扑分类 | 6 | 1 (classifyFaces) |
| 端到端流水线 | 5 | 1 (STEP→Mesh→MSH) |
| Python API + 对比脚本 | 19 脚本 + 8 单元测试 | 1 (gmsh Python API) |
| **合计** | **~156 测试函数 + ~19 脚本** | **~21 GMSH 功能域** |
