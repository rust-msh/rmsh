# Milestone 10 实现计划：平台与部署

**目标**：完成 WASM 部署、Cloud 模式、版本分级与 Local-First 支持

**当前进度**：20% (部分功能已完成)

---

## 1. 当前状态总结

### ✅ 已完成部分

| 功能 | 模块 | 状态 |
|------|------|------|
| Edition 功能门控 (Basic/Pro/Enterprise) | `app` / `domain` | ✅ |
| Web 版工程管理适配（禁用新建/另存） | `app` | ✅ |
| OPFS 存储层设计 | `worker` | ✅ |
| Worker 通信协议 (MessagePack) | `domain` | ✅ |
| WasmBackend trait 实现 | `infra` | ✅ |
| 主线程 WASM 入口 (`start_with_canvas_id`) | `main` | ✅ |
| Worker 基础框架 (onmessage 分发) | `worker` | ✅ |

### 🔲 需要完成的任务

#### 第一阶段：WASM 构建与基础设施 (优先级：⭐⭐⭐⭐⭐)

| 任务 | 模块 | 说明 | 预期工作量 |
|------|------|------|----------|
| **1.1** WASM 构建脚本 | `main` | 创建 `build-wasm.sh` 调用 trunk，处理多目标编译 | 2h |
| **1.2** WebGPU 双栈适配 | `render` | 确保 wgpu 在 Native/WASM 两个目标上工作，处理 feature 差异 | 3h |
| **1.3** rmsh/Rem WASM target 配置 | `solver` | 添加 `[target.'cfg(target_arch = "wasm32")']` 特殊处理，可能需要 C++ binding 适配 | 4h |
| **1.4** Worker WASM 单独编译 | `worker` | 配置 trunk 生成独立的 `worker.js` | 2h |
| **1.5** 构建验证测试 | CI | 创建 GitHub Actions 工作流验证 WASM 编译 | 2h |

**小计 M1**：~13 小时

#### 第二阶段：Local-First 离线能力 (优先级：⭐⭐⭐⭐)

| 任务 | 模块 | 说明 | 预期工作量 |
|------|------|------|----------|
| **2.1** Service Worker 脚本生成 | `main` | 生成 `dist/sw.js`，缓存 WASM 模块 + 静态资源 | 3h |
| **2.2** Service Worker 激活与清理 | `main` | 版本控制、缓存清理、降级处理 | 2h |
| **2.3** OPFS 文件系统初始化 | `worker` | 保证 OPFS 目录结构完整、权限正确 | 2h |
| **2.4** 网络离线检测 | `app` | 主线程检测网络状态，UI 显示离线提示 | 2h |
| **2.5** 离线下 Worker 通信保证 | `infra` | 处理 Worker 缺失、IDB 存储降级方案 | 2h |

**小计 M2**：~11 小时

#### 第三阶段：Web 工程加载流程 (优先级：⭐⭐⭐⭐)

| 任务 | 模块 | 说明 | 预期工作量 |
|------|------|------|----------|
| **3.1** 预创建工程加载 API | `infra` | `load_predefined_project(project_id)` 从 Worker 加载预创建的 `.emsp` | 2h |
| **3.2** 工程索引管理 | `worker` | 维护 `__projects_index.json`，快速列举工程列表 | 2h |
| **3.3** 加载进度回调 | `app` | UI 展示"加载中..."，进度条更新 | 1h |
| **3.4** 工程预热 | `app` | 加载后自动解析、验证、初始化 UI 状态 | 2h |
| **3.5** 工程切换（多 Tab 加载） | `app` | 支持同时打开多个预创建工程 | 2h |

**小计 M3**：~9 小时

#### 第四阶段：Cloud 模式设计 (优先级：⭐⭐⭐)

| 任务 | 模块 | 说明 | 预期工作量 |
|------|------|------|----------|
| **4.1** Cloud 后端 API 规范 | `docs` | 设计 RESTful API (POST /projects, GET /projects/:id, PATCH 等) | 3h |
| **4.2** 认证与授权框架 | `infra` | JWT token / API key 支持，CORS 设置 | 3h |
| **4.3** 工程上传/下载流式处理 | `infra` | 大文件分片传输、断点续传、流式解析 | 5h |
| **4.4** 多用户并发锁机制 | `infra` | 工程锁、操作冲突解决、merge 策略 | 4h |
| **4.5** 云端求解任务队列 | `infra` | 任务提交、进度查询、结果回调 | 4h |
| **4.6** 云端 Webhook + 长连接 | `infra` | 结果通知推送、实时仿真进度 | 3h |

**小计 M4**：~22 小时（可分次迭代）

---

## 2. implementation 细节

### Phase 1.1：WASM 构建脚本

**需要创建**：`scripts/build-wasm.sh`

```bash
#!/bin/bash
set -e

# Build worker WASM
cd crates/worker
wasm-pack build --target web --out-dir ../../dist/wasm/worker

# Build main WASM via trunk
cd ../main
trunk build --release

# Output: dist/index.html + dist/emstudio_main.js + dist/worker/...
```

**修改内容**：
1. [main/Trunk.toml](../crates/main/Trunk.toml) 确保 `dist` 指向正确的构建输出目录
2. [solver/Cargo.toml](../crates/solver/Cargo.toml) 添加 WASM-only feature flag

### Phase 1.2：WebGPU 适配检查清单

#### 检查项目：

- [ ] `wgpu` v27 支持 WebGPU 后端（检查 features）
- [ ] `render` 模块中 wgpu device 初始化代码是否处理 WebGPU 特有限制
- [ ] eframe + egui-wgpu 是否都支持 WASM 运行时
- [ ] 着色器 (WGSL) 代码是否兼容 WebGPU 规范
- [ ] 纹理、缓冲区格式是否 WebGPU 兼容
- [ ] 渲染管线是否处理 WebGPU 的限制 (如 workgroup size)

#### 关键代码区域：

- `render/src/lib.rs`：Device 初始化
- `render/src/shaders/*.wgsl`：着色器代码
- `main/src/lib.rs`：eframe 启动配置

### Phase 1.3：rmsh/Rem WASM 目标编译

**关键问题**：
- rmsh 和 Rem 可能包含 C/C++ 代码（通过 FFI 绑定）
- 需要检查是否支持 WASM 编译或需要提供纯 Rust 替代方案

**检查步骤**：
1. `cargo build --target wasm32-unknown-unknown --package emstudio-solver`
   - 记录所有编译错误
   - 识别 WASM 不兼容的依赖

2. 对于 C/C++ 绑定生成的代码：
   - 使用 `#[cfg(not(target_arch = "wasm32"))]` 隔离
   - 或提供 WASM-only shim 实现

**创建文件**：[solver/src/wasm_compat.rs](../crates/solver/src/wasm_compat.rs)

```rust
#[cfg(target_arch = "wasm32")]
pub mod wasm_only_solver {
    // WASM-compatible solver stub
    // 将真实求解延迟到后端（Cloud/Native）
}

#[cfg(not(target_arch = "wasm32"))]
pub use native_solver::*;
```

### Phase 1.4：Worker WASM 独立编译

修改 [main/Trunk.toml](../crates/main/Trunk.toml)：

```toml
[[hooks]]
stage = "post_build"
command = "sh"
command_args = ["-c", "cd ../worker && wasm-pack build --target web --out-dir ../../dist/worker"]
```

结果：`dist/worker_bg.wasm` 和 `worker_bg.js`，主线程加载：

```javascript
// in HTML or JavaScript
const worker = new Worker('/worker_bg.js');
```

### Phase 2.1：Service Worker 脚本模板

**创建**：`main/src/service_worker.js`（注意：这是 JavaScript，不编译）

```javascript
const CACHE_VERSION = 'emstudio-v1';
const ASSETS_TO_CACHE = [
  '/',
  '/index.html',
  '/emstudio_main.js',
  '/emstudio_main_bg.wasm',
  '/worker_bg.js',
  '/worker_bg.wasm',
];

self.addEventListener('install', (e) => {
  e.waitUntil(
    caches.open(CACHE_VERSION).then((cache) => {
      return cache.addAll(ASSETS_TO_CACHE);
    })
  );
  self.skipWaiting();
});

self.addEventListener('fetch', (e) => {
  e.respondWith(
    caches.match(e.request).then((resp) => {
      return resp || fetch(e.request);
    })
  );
});
```

修改 [main/index.html](../crates/main/index.html)：

```html
<script>
  if ('serviceWorker' in navigator) {
    navigator.serviceWorker.register('/sw.js');
  }
</script>
```

Trunk 构建后复制 `src/service_worker.js` → `dist/sw.js`

### Phase 2.3：OPFS 初始化确保

修改 [worker/src/opfs.rs](../crates/worker/src/opfs.rs)：

```rust
pub async fn init_filesystem() -> Result<(), String> {
    let StorageManager = storage_manager().map_err(|e| format!("{:?}", e))?;
    let root = ensure_root_dir().await?;
    ensure_projects_dir().await?;
    ensure_results_dir().await?;
    // ...
    Ok(())
}
```

然后在 [worker/src/lib.rs](../crates/worker/src/lib.rs) 的 entry point 调用：

```rust
#[wasm_bindgen(start)]
pub async fn main() {
    // Initialize OPFS
    if let Err(e) = opfs::init_filesystem().await {
        web_sys::console::error_1(&e.into());
    }
    // Then setup message handler...
}
```

### Phase 3.1：预创建工程加载 API

修改 [infra/src/lib.rs](../crates/infra/src/lib.rs)，添加新方法到 `Backend` trait：

```rust
pub trait Backend: Send + Sync {
    // ... existing methods ...
    
    async fn load_predefined_project(
        &self,
        project_id: &str,
    ) -> Result<Project, BackendError>;
}
```

在 `WasmBackend` 中实现：

```rust
impl Backend for WasmBackend {
    async fn load_predefined_project(&self, project_id: &str) -> Result<Project, BackendError> {
        let cmd = WorkerCommand::LoadProject { id: project_id.to_string() };
        self.send_worker_command(cmd).await
    }
}
```

### Phase 4.1：Cloud REST API 规范

**创建**：`docs/cloud-api-spec.md`

关键端点：

```
POST   /api/v1/projects                 # 创建工程
GET    /api/v1/projects                 # 列举工程
GET    /api/v1/projects/:project_id     # 获取工程
PATCH  /api/v1/projects/:project_id     # 更新工程
DELETE /api/v1/projects/:project_id     # 删除工程
POST   /api/v1/projects/:project_id/solve  # 提交仿真任务
GET    /api/v1/solve-tasks/:task_id     # 查询任务进度
GET    /api/v1/projects/:project_id/results/:result_id  # 获取结果
```

---

## 3. 依赖关系与可并行化任务

```
┌─────────────────────────────────────────────┐
│   Phase 1: WASM 基础 (序列化执行)            │
│  ├─ 1.1 构建脚本 (2h)                       │
│  ├─ 1.2 WebGPU 适配 (3h)  ┐                │
│  ├─ 1.3 rmsh/Rem WASM (4h) ├─ 可并行      │
│  └─ 1.4 Worker 编译 (2h)   ┘                │
└──────────────┬──────────────────────────────┘
               │ (depends on Phase 1 completion)
┌──────────────▼───────────────────────────────┐
│  Phase 2 & 3: 离线 + 工程加载 (可并行)       │
│  ├─ Phase 2: Local-First (11h)              │
│  └─ Phase 3: Web 工程流程 (9h)              │
└──────────────┬──────────────────────────────┘
               │ (after Phase 1 & 2/3 complete)
┌──────────────▼──────────────────────────────┐
│  Phase 4: Cloud API (按需迭代, 可与前面并行) │
│  ├─ 4.1-4.3 基础 API (11h)                  │
│  └─ 4.4-4.6 高级功能 (11h)                  │
└──────────────────────────────────────────────┘
```

---

## 4. 测试策略

### 单元测试

- [ ] `domain::worker_protocol` 消息编解码
- [ ] `infra::WasmBackend` 命令分发
- [ ] `domain::expression` 在 WASM 上运行

### 集成测试

- [ ] WASM 构建是否生成有效的二进制
- [ ] 主线程与 Worker 通信无死锁
- [ ] 离线模式下工程保存/加载不丢失
- [ ] 预创建工程加载后能正常渲染

### 浏览器兼容性测试

- [ ] Chrome 123+ (WebGPU stable)
- [ ] Firefox 121+ (experimental WebGPU)
- [ ] Edge 123+ (WebGPU)
- [ ] Safari 18+ (WebGPU preview)

---

## 5. 时间表与里程碑

| Phase | 预期完成时间 | 关键产物 |
|-------|------------|--------|
| M1 (WASM 基础) | 2-3 天 | 可用的 WASM 构建、WebGPU 适配完成 |
| M2 + M3 (离线 + 工程加载) | 3-4 天 | Service Worker、OPFS 初始化、工程预加载 |
| M4 (Cloud API) | 1-2 周 | API 规范、原型实现 |
| **总计** | **3-4 周** | **M10 完成，可发布 Web 基本版** |

---

## 6. 风险与缓解

| 风险 | 概率 | 影响 | 缓解方案 |
|------|------|------|--------|
| rmsh/Rem WASM 编译失败 | 高 | 阻塞主线程 | 提前一周调研 C++ binding，准备 shim 方案 |
| WebGPU 浏览器支持不足 | 中 | 用户基数小 | 同时支持 Canvas 2D 降级渲染 |
| OPFS 权限问题 | 低 | 工程保存失败 | 测试多浏览器，准备 IDB 备用方案 |
| Worker 与主线程通信延迟 | 低 | UI 卡顿 | 实现超时和异步进度回调 |

---

## 7. 代码审查清单

每个 Phase 完成后：

- [ ] 所有代码有清晰的 WASM/Native 边界标记
- [ ] 错误处理完整（Web API 可能返回 Promise 拒绝）
- [ ] 内存泄漏检查（特别是 JS-Rust 交界处）
- [ ] 性能分析（WebAssembly profiler：Chrome DevTools）
- [ ] 文档更新（README 新增 WASM 构建说明）

