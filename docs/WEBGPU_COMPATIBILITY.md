# WebGPU 兼容性检查清单

**目标**：确保 EMStudio 渲染引擎在 Native、WASM (WebGPU)、浏览器上无缝运行

**版本**：wgpu 0.27+ WebGPU 适配

---

## 1. 依赖版本确认

### ✅ 已检查项目

| 依赖 | 版本要求 | Status | 备注 |
|------|---------|--------|------|
| wgpu | ^27 | ✅ | WebGPU 正式支持 |
| egui_wgpu | 0.33 | ✅ | 完全支持 WASM |
| eframe | 0.33 | ✅ | WASM runner 包含 |
| egui | 0.33 | ✅ | 跨平台兼容 |

---

## 2. 着色器 (WGSL) 兼容性

### 检查清单

- [x] 所有着色器文件使用 WGSL（不是 GLSL/SPIR-V）
  - `render/src/shaders/*.wgsl`
  
- [x] 着色器中没有平台特定的扩展
  - NO: `#version 450`
  - NO: glsl 风格的 GLSL 指令
  - YES: `@vertex`, `@fragment` 属性语法

- [x] 数据布局遵守 WGSL 规范
  - `@align(16)` 对齐限制
  - `@size()` 显式大小指定

- [x] 工作组大小 ≤ 256（WebGPU 限制）
  ```wgsl
  // ✅ OK:
  @compute @workgroup_size(32, 32)
  
  // ❌ NO:
  @compute @workgroup_size(512)  // Too large for WebGPU
  ```

### 相关文件

- `render/src/field_pipeline.rs` → 着色器加载的地方
  ```rust
  let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
      label: Some("field_shader"),
      source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(include_str!("shaders/field.wgsl"))),
  });
  ```

---

## 3. 纹理与缓冲区格式

### WebGPU 支持格式清单

#### ✅ 推荐使用的贴图格式

| 格式 | 用途 | 状态 |
|------|------|------|
| `Rgba8UnormSrgb` | LDR 颜色输出 | ✅ 完全支持 ✅ ✅ |
| `Rgba8Unorm` | 线性 RGBA | ✅ 完全支持 |
| `R32Float` | 浮点场数据 | ✅ 完全支持 |
| `Rg32Float` | 矢量字段 | ✅ 完全支持 |
| `Rgba32Float` | 高精度数据 | ⚠️ 受限（仅某些应用） |
| `Depth32Float` | 深度缓冲区 | ✅ 完全支持 |
| `Depth24Plus` | 深度缓冲区 | ❌ 仅 Native |

#### 检查项目

- [x] `render/src/field_pipeline.rs` 中的所有纹理格式在上表中
- [x] 如果使用 `Depth24Plus`，需要提供 WebGPU 替代：
  ```rust
  #[cfg(target_arch = "wasm32")]
  const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;
  
  #[cfg(not(target_arch = "wasm32"))]
  const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth24Plus;
  ```

- [x] 离屏帧缓冲区格式
  - 当前：何处定义？(`field_pipeline.rs` line ~242)
  - 应该是：`Rgba8UnormSrgb`（WebGPU 标准）

---

## 4. 管线配置兼容性

### WebGPU 限制（vs Native）

#### ✅ 已支持

| 特性 | 状态 | 说明 |
|------|------|------|
| PrimitiveTopology::TriangleList | ✅ | 标准三角形渲染 |
| 无索引顶点（隐式索引） | ✅ | `draw(0..vertex_count)` |
| 索引缓冲区 | ✅ | `draw_indexed()` |
| 实例化渲染 | ✅ | `draw_indexed_indirect()` |
| 深度测试/写入 | ✅ | Depth24Plus/Depth32Float |

#### ⚠️ 需要检查

| 特性 | 状态 | WebGPU 限制 |
|------|------|-----------|
| 双面渲染 | ⚠️ | 需要分别编译前/后 |
| 多采样 MSAA | ⚠️ | 仅 4/8/16 (不支持 2) |
| 保守光栅化 | ❌ | 不支持 |
| 动态缓冲区偏移数 | ⚠️ | max 8（通常够） |

#### 检查位置

- `render/src/field_pipeline.rs` → `RenderPipelineDescriptor` 配置
- `render/src/arrow_pipeline.rs` → arrow 实例化配置

查找代码：
```rust
// Line ~162
let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
    label: Some("field_pipeline"),
    layout: Some(&pipeline_layout),
    push_constant_ranges: &[],  // ← Check this
    vertex: wgpu::VertexState {
        module: &shader,
        entry_point: "vs_main",
        buffers: &vertex_buffers,
    },
    primitive: wgpu::PrimitiveState {
        topology: wgpu::PrimitiveTopology::TriangleList,
        strip_index_format: None,
        front_face: wgpu::FrontFace::Ccw,
        cull_mode: Some(wgpu::Face::Back),  // ← OK for WebGPU
        polygon_mode: wgpu::PolygonMode::Fill,  // ← OK
        unclipped_depth: false,  // ← OK
    },
    depth_stencil: Some(wgpu::DepthStencilState {
        format: wgpu::TextureFormat::Depth32Float,  // ← Must check format
        depth_write_enabled: true,
        depth_compare: wgpu::CompareFunction::Less,
        stencil: wgpu::StencilState::default(),
        bias: wgpu::DepthBiasState::default(),
    }),
    // ...
});
```

---

## 5. 设备初始化与功能检查

### 当前代码路径

**需要检查**：如何在 eframe + egui-wgpu 中初始化 wgpu？

eframe 在 WASM 版本中使用 `WebRunner`，它处理所有 wgpu 初始化。关键是验证：

- [x] `WebRunner::new().start()` 正确传递 WebGPU 支持
- [x] 没有硬编码的后端偏好 (Vulkan/Metal/DX12)

**main/src/lib.rs**：
```rust
WebRunner::new()
    .start(
        canvas,
        eframe::WebOptions::default(),  // ← Uses WebGPU by default on wasm32
        Box::new(|cc| Ok(Box::new(App::new(...)))),
    )
    .await
```

**可选优化**（如果需要显式控制）：

```rust
#[cfg(target_arch = "wasm32")]
let mut web_options = eframe::WebOptions::default();
web_options.prefer_dark_mode = Some(true);
// Add WebGPU instance configuration if available
```

---

## 6. 内存与性能限制

### WebGPU 限制 (vs Native)

| 限制 | WebGPU | Native | 影响 |
|------|--------|--------|------|
| 缓冲区最大大小 | 2^53 字节 | 理论无限 | 低（通常不是瓶颈） |
| 纹理最大分辨率 | 8192×8192 (典型) | 16384+ | 中（字段可视化可能受限） |
| 绑定组大小 | 640 | 无限制 | 低 |
| 存储绑定数 | 4-8 | 无限制 | 中（计算着色器） |
| 计算工作组大小 | 256 线程 | 1024+ | 高（并行场处理） |

### 缓解措施

- [x] 如果字段可视化纹理 > 4096×4096，需要分块/mipmap（当前未触发：默认离屏尺寸随视口）
- [x] 计算着色器（如果有）workgroup 大小 ≤ 256（当前渲染链路未使用 compute shader）
- [x] 将大型计算预处理到 CPU (Native) 或推迟到后端 (WASM)（已采用 Worker + 后端策略）

---

## 7. 浏览器兼容性

### 浏览器 WebGPU 支持状态（2026.04）

| 浏览器 | 版本 | WebGPU 状态 | 测试状态 |
|--------|------|-----------|---------|
| Chrome | 123+ | ✅ 稳定 | 🔲 待测 |
| Edge | 123+ | ✅ 稳定 | 🔲 待测 |
| Firefox | 121+ | ⚠️ 实验性 | 🔲 待测 |
| Safari | 18+ | ⚠️ 预览 | 🔲 待测 |
| Opera | 109+ | ✅ 稳定 | 🔲 待测 |

### 测试环境建议

```bash
# 启用 WebGPU 标志
# Chrome: chrome://flags/#enable-webgpu (already default in 123+)
# Firefox: about:config dom.webgpu.enabled = true
# Safari: Develop menu → Experimental Features → WebGPU
```

---

## 8. 构建配置

### Cargo.toml WebGPU Features

**rootCargo.toml**：
```toml
[workspace.dependencies]
wgpu = "27"  # Implicit default features include WebGPU
```

**每个 crate**（如 render）：
```toml
[dependencies]
wgpu.workspace = true  # Already supports WebGPU via workspace
```

**检查**：
```bash
cargo tree -i wgpu | grep -i webgpu
```

---

## 9. 调试 WASM 渲染问题

### Chrome DevTools

1. **Performance 标签**
   - 录制一个帧
   - 查看是否有 WebGPU 调用
   - 标色检查 GPU 时间

2. **Console**
   ```javascript
   // 检查 WebGPU 支持
   console.log(navigator.gpu);
   
   // 列举 Adapter
   const adapter = await navigator.gpu.requestAdapter();
   const info = await adapter.requestAdapterInfo();
   console.log(info);
   ```

3. **Coverage 标签**
   - 查看 WASM 模块是否被加载
   - 检查着色器是否被编译

### 常见问题排查

| 症状 | 可能原因 | 解决方案 |
|------|---------|--------|
| 黑屏 | 渲染管线未初始化 | 检查是否有 JS/WASM 错误 (Console) |
| 纹理不显示 | 格式不兼容 | 使用 `Rgba8UnormSrgb` 替代 `Rgba8Unorm` |
| 着色器编译失败 | WGSL 语法错误 | 检查着色器文件，run locally with `cargo check` |
| 低帧率 | 性能瓶颈 | 减少顶点数、降低纹理分辨率、启用 LOD |
| 工作线程崩溃 | Worker WASM 加载失败 | 检查 `window.WORKER_PATH` 和 fetch 错误 |

---

## 10. 验证清单（构建前）

- [ ] 所有 `.wgsl` 文件验证（无 GLSL 或 SPIR-V）
- [ ] 纹理格式在 `src/wasm_compat.rs` 中检查 (`Depth24Plus` → 替代方案)
- [ ] 渲染管线配置无 WebGPU 不兼容特性
- [ ] Device 初始化通过 eframe WebRunner（默认）
- [ ] WASM panic hooks 已配置 (`console_error_panic_hook`)
- [ ] main.rs 正确处理 canvas 和 Service Worker

### 运行验证

```bash
# 1. 检查 WASM 构建（不需要运行）
cd crates/main
cargo check --target wasm32-unknown-unknown

# 2. 检查着色器
find crates/render/src/shaders -name "*.wgsl" -exec wgsl-analyzer {} \;

# 3. 构建完整的 WASM
./scripts/build-wasm.sh

# 4. 启动本地服务器并测试
python3 -m http.server 8000 --directory crates/main/dist
# 访问 http://localhost:8000
```

---

## 11. 性能优化建议（如适用）

### 对于 WASM

1. **减少 CPU 往复**
   - 批量上传缓冲区数据
   - 使用 indirect 绘制以减少 JS/WASM 调用

2. **异步纹理加载**
   ```rust
   // ✅ 好：异步加载字段数据→GPU 纹理
   async fn load_field_data_async(&self, url: &str) {
       let bytes = fetch_bytes(url).await;
       self.upload_to_gpu(&bytes);
   }
   
   // ❌ 差：同步请求阻塞 UI
   let bytes = std::fs::read(path)?;
   ```

3. **纹理压缩**
   - 如果支持，使用 BC/ASTC 压缩
   - 注意：某些格式 WebGPU 可能不支持

---

## 检查流程总结

```
1. [Code Review]
   ├─ 检查着色器格式 (WGSL) ✅
   ├─ 检查纹理格式兼容性
   ├─ 检查管线配置
   └─ 检查 device 初始化路径

2. [构建验证]
   ├─ cargo check --target wasm32-unknown-unknown
   └─ ./scripts/build-wasm.sh

3. [Runtime 检查]
   ├─ Chrome DevTools → Console (错误)
   ├─ Chrome DevTools → Performance (性能)
   └─ 测试多个浏览器和操作系统

4. [部署]
   ├─ Service Worker 缓存 WASM 模块
   ├─ CORS 头配置正确
   └─ 内存使用监控
```

---

## 12. 本次实测结论（2026-04-06）

- `emstudio-render` 在 `wasm32-unknown-unknown` 目标编译通过。
- 渲染关键格式已满足 WebGPU 兼容要求：
   - color: `Rgba8UnormSrgb`
   - depth: `Depth32Float`
   - topology: `TriangleList` / `LineList`
   - 无 `Depth24Plus` 依赖。
- WASM 入口使用 `WebRunner::new().start(..., WebOptions::default(), ...)`，未硬编码 Native 后端。
- 当前阻塞项不在 WebGPU 渲染层，而在 `emstudio-main` 的依赖树（`getrandom` wasm `js` feature）。

