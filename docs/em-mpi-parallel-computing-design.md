# EMStudio MPI 并行计算方案调研与设计

## 1. Ansys HFSS / Q3D 并行计算体系调研

### 1.1 HFSS 并行层次总览

HFSS 提供了**多层次、可嵌套**的并行计算架构，从线程级到集群级覆盖完整：

```
┌─────────────────────────────────────────────────────────────────────┐
│                        HFSS 并行层次                                 │
│                                                                     │
│  Level 0: 矩阵多处理 (MP)                                           │
│    └─ 单机多线程，共享内存，加速矩阵 LU 分解                          │
│                                                                     │
│  Level 1: 频谱分解 (SDM — Spectral Decomposition Method)            │
│    └─ 将各频率点分配到不同核/节点，独立求解                            │
│                                                                     │
│  Level 2: 域分解 (DDM — Domain Decomposition Method)                │
│    └─ 将单频点的 FEM 网格空间分区，MPI 分布式求解                      │
│                                                                     │
│  Level 3: 分布式求解 (DSO — Distributed Solve Option)               │
│    └─ 参数扫描/优化/阵列的各独立变量组合分配到不同机器                  │
│                                                                     │
│  Level 4: 多层嵌套 (Multi-Tier HPC)                                 │
│    └─ SDM（外层频率并行）+ DDM（内层网格并行）组合                     │
└─────────────────────────────────────────────────────────────────────┘
```

### 1.2 各并行方法详解

#### 1.2.1 矩阵多处理 MP（共享内存）

- **并行粒度**：单个矩阵求解内的线程并行
- **内存模型**：共享内存（单机多核）
- **适用场景**：单机加速，所有求解器默认启用
- **技术细节**：FEM 系统矩阵 Ax=b 的 LU 分解或 GMRES 迭代均可使用多线程加速
- **限制**：受单机内存限制，无法超越单节点 RAM 上限

#### 1.2.2 频谱分解 SDM（频率级并行）

- **并行粒度**：各频率点独立求解
- **内存模型**：各任务独立内存空间
- **适用场景**：频率扫描中有大量频率点，且单频点问题可装入单节点内存
- **扩展性**：近线性扩展（各频率点完全独立，无通信开销）
- **额外能力**：多频率自适应网格加密（Multi-Frequency Adaptive Meshing）也可并行化

#### 1.2.3 域分解 DDM（网格级并行，MPI）

这是 HFSS 最核心的 MPI 并行方法：

- **并行粒度**：将单频点的 FEM 网格划分为 N-1 个不重叠空间子域
- **内存模型**：分布式内存（MPI 进程间通信）
- **最低要求**：至少 3 个 MPI 任务（1 个 Master + 2 个子域求解器）
- **Master 节点角色**：
  - 负责域组装、网格加密决策、解管理
  - **不参与子域求解**（不分配子域给 Master）
- **子域求解**：
  - 每个子域在本地用**直接求解器（LU 分解）**独立求解
  - 子域边界通过全局**迭代求解器**耦合收敛
- **并行阶段**：网格生成、矩阵组装、求解、场恢复全部并行化
- **推荐配置**：每个 MPI 任务分配 2~3 个 CPU 核（混合 MPI+线程）
- **网络要求**：推荐 InfiniBand 高速互连；千兆以太网可用但较慢

**DDM 求解流程**：

```
┌──────────┐    划分网格      ┌───────────┐
│          │  ──────────────► │ 子域 0    │ ─── MPI Rank 1
│  Master  │                  ├───────────┤
│ (Rank 0) │  ──────────────► │ 子域 1    │ ─── MPI Rank 2
│          │                  ├───────────┤
│ 不求解    │  ──────────────► │ 子域 2    │ ─── MPI Rank 3
│ 只协调    │                  ├───────────┤
│          │  ──────────────► │  ...      │ ─── MPI Rank N
└─────┬────┘                  └─────┬─────┘
      │                             │
      │    迭代交换边界数据（MPI）      │
      ◄─────────────────────────────►
      │                             │
      │    全局迭代收敛后汇总解        │
      ◄─────────────────────────────┘
```

**DDM 技术还支撑了以下衍生功能**：

| 技术 | 说明 |
|------|------|
| FE-BI 混合求解 | 有限元 + 边界积分混合 |
| IE 积分方程区域 | 积分方程子域 |
| 3D Component Array | 阵列域分解 |
| Mesh Fusion | 独立网格融合 |
| MLFMM | 多级快速多极子加速 |

#### 1.2.4 分布式求解 DSO（任务级并行）

- **并行粒度**：完全独立的仿真任务
- **适用目标**：
  - 参数扫描的各变量组合
  - 离散频率扫描的各频率点
  - 阵列各单元贡献
  - Optimetrics 各设计点
- **特点**：「尴尬并行」——零通信，近线性扩展
- **关键约束**：**DSO 与 DDM 互斥**——启用 DDM 后，DSO 自动禁用

#### 1.2.5 多层嵌套 HPC

现代 HFSS 支持两层嵌套：

```
128 核集群示例：
┌──────────────────────────────────────────┐
│            外层：SDM（4 个频率组）          │
│  ┌─────────┬─────────┬─────────┬────────┐ │
│  │ f1      │ f2      │ f3      │ f4     │ │
│  │ 32 核   │ 32 核   │ 32 核   │ 32 核  │ │
│  │ DDM     │ DDM     │ DDM     │ DDM    │ │
│  │ (内层)  │ (内层)  │ (内层)  │ (内层) │ │
│  └─────────┴─────────┴─────────┴────────┘ │
└──────────────────────────────────────────┘
```

通过 HFSS 的 "Two-Level Distribution" 设置配置。

### 1.3 HFSS 矩阵求解器

| 求解器 | 方法 | 内存特点 | 分布式支持 | 适用场景 |
|--------|------|---------|-----------|---------|
| 直接求解器 | LU 分解 | 高内存（存储完整分解矩阵） | 支持分布式内存分解 | 中小规模，高鲁棒性 |
| 迭代求解器 | GMRES (Krylov) | 低内存（矩阵-向量乘法） | 支持分布式 | 电大尺寸问题，内存瓶颈 |
| DDM 耦合求解器 | 迭代器协调子域 | 子域各自独立 + 边界迭代 | 原生 MPI | 超大规模模型 |

### 1.4 Q3D Extractor 并行计算

Q3D 使用**三种不同的求解器**，并行策略各异：

#### Q3D 求解器架构

| 提取类型 | 求解方法 | 网格类型 | 并行方法 |
|---------|---------|---------|---------|
| CG（电容+电导） | BEM/MoM + FMM 加速 | 表面三角形网格 | MPI 分布式矩阵组装+求解 |
| DC RL（直流电阻+电感） | FEM 体积求解 | 四面体体网格 | 共享内存多线程 |
| AC RL（交流电阻+电感） | 表面积分/BEM | 表面网格 | 共享内存 + MPI |

#### Q3D CG 分布式求解（MPI）

这是 Q3D 最核心的 MPI 场景：

- CG 求解器产生**稠密矩阵**（Green 函数导致 N^2 交互），计算和内存密集
- MPI 将矩阵组装和求解分布到多个 rank
- **推荐配置**：每个 MPI 任务约 4 个 CPU 核
- **约束**：CG Matrix Solver 分布与 Solver Nets/Sources Distribution Type 互斥
- **适用场景**：问题规模超出单节点 RAM 时

```
Q3D CG 分布式求解流程：

  导体表面网格 (N 个三角形)
        │
        ▼
  Green 函数矩阵 G[i][j] (N×N 稠密矩阵)
        │
        ▼
  ┌──── MPI 分布 ────┐
  │   Rank 0: 行 0~k  │
  │   Rank 1: 行 k~2k │
  │   ...              │
  │   Rank P: 行 ..~N  │
  └────────────────────┘
        │
        ▼
  并行矩阵求解 → RLCG 矩阵
```

#### Q3D 支持的 HPC 方法

| 方法 | Q3D 支持 |
|------|---------|
| 多线程 (MP) | ✅ |
| 频谱分解 (SDM) | ✅（多频率并行） |
| 域分解 (DDM) | ✅ |
| MPI 分布式 | ✅（CG 矩阵求解） |
| DSO | ✅（参数扫描/Optimetrics） |

### 1.5 MPI 实现选型

| 平台 | 默认 MPI | 备选 MPI | 说明 |
|------|---------|---------|------|
| Windows | Intel MPI 2021.x | Microsoft MPI 10.1 | Intel MPI 随安装包捆绑 |
| Linux | Intel MPI 2021.x | Open MPI 4.0.5 | Intel MPI 为主要支持 |

- Intel MPI 单机运行**零配置**（捆绑安装）
- 多机运行需配置凭据和防火墙规则
- MS-MPI 需额外安装 `msmpisetup.exe` 并启动 `MsMpiLaunchSvc` 服务

### 1.6 HPC 作业调度集成

Ansys Electronics Desktop 原生支持以下调度器：

| 调度器 | 协议变体 |
|--------|---------|
| SLURM | SLURM / SLURMSSH / SLURMSMB |
| PBS Professional | PBS / PBSSSH / PBSSMB |
| IBM Spectrum LSF | LSF / LSFSSH / LSFSMB |
| Sun Grid Engine (SGE/UGE) | UGE / SGE |
| Windows HPC | Windows HPC |
| Ansys Cloud | Ansys Cloud |

通过 **RSM（Remote Solve Manager）** 中间件实现跨平台作业提交（如 Windows → Linux SLURM）。

### 1.7 HFSS HPC 配置模式

**自动模式**：用户仅指定可用硬件（机器列表 + 总核数），HFSS 自动决定任务分配策略。

**手动模式**：用户显式配置：
- `num_tasks`：MPI 进程数
- `cores_per_task`：每个 MPI 进程内的共享内存线程数
- 各主机的任务/核分配

**关键原则**：
- 不要超分配（线程数 ≤ 物理核数）
- 内存受限时：减少每节点 task 数（每个 task 分到更多 RAM）
- DDM 无 FE-BI/IE 时：每 task 约 2~3 核最优
- 每个 HFSS license 自带 4 个免费 HPC 增量；超出需 HPC Pack 授权

### 1.8 分布式求解的基础设施要求

| 要求 | 说明 |
|------|------|
| 软件一致性 | 所有计算节点安装相同版本 HFSS + MPI |
| 硬件同构性 | 推荐相同 OS、架构、相似硬件配置 |
| 共享存储 | NFS/SMB 共享文件系统，或各节点相同路径的本地安装 |
| 项目文件可达性 | 所有节点可访问工程文件和临时目录 |
| 网络 | 防火墙允许 MPI 通信端口；推荐 InfiniBand |
| 授权 | 足够的 HPC license token |

---

## 2. EMStudio 并行计算方案设计

基于对 Ansys HFSS/Q3D 并行体系的调研，结合 EMStudio 的 Rust 技术栈和 Rem 求解器后端，设计以下并行计算方案。

### 2.1 设计原则

| 原则 | 说明 |
|------|------|
| **渐进式实现** | 从单机多线程开始，逐步扩展到 MPI 分布式 |
| **Rem 原生能力优先** | Rem 求解器自身的并行能力是基础，EMStudio 负责调度和编排 |
| **Rust 生态适配** | 优先使用 Rust 原生并行（rayon/tokio），MPI 层通过 FFI 绑定 |
| **WASM 全能力运行** | 通过 [jsmpi](https://github.com/jsmpi/jsmpi) 在浏览器中用 Web Worker 模拟 MPI，实现 WASM 环境下的本地并行求解 |
| **配置驱动** | 并行策略通过 HPC 配置文件控制，不硬编码 |
| **单机先行** | 优先保证单机多线程/多进程体验，多机 MPI 作为高级功能 |

### 2.2 并行层次设计

对标 HFSS 的多层并行体系，EMStudio 规划以下层次：

```
┌─────────────────────────────────────────────────────────────────┐
│                     EMStudio 并行层次                             │
│                                                                 │
│  L0: 线程并行 (Thread Parallelism)                               │
│    ├─ Rem 求解器内部的多线程矩阵运算                               │
│    ├─ rayon 并行迭代器（网格生成、后处理）                          │
│    ├─ wgpu GPU 并行渲染（可视化层，已实现）                         │
│    └─ [WASM] Web Worker 线程（通过 jsmpi 模拟 MPI rank）          │
│                                                                 │
│  L1: 频率并行 (Spectral Decomposition)                           │
│    ├─ [Native] 各频率点独立 Rem 求解进程                           │
│    └─ [WASM]  各频率点独立 Web Worker（jsmpi rank）               │
│                                                                 │
│  L2: 域分解 (Domain Decomposition)                               │
│    ├─ [Native] FEM 网格空间分区，MPI 分布式求解（Rem + rsmpi）     │
│    └─ [WASM]  FEM 网格空间分区，Web Worker 分布式求解（Rem + jsmpi）│
│                                                                 │
│  L3: 任务并行 (Task Distribution)                                │
│    ├─ [Native] 参数扫描 / Optimetrics 各变量组合独立求解进程        │
│    └─ [WASM]  各变量组合独立 Web Worker                           │
│                                                                 │
│  L4: 混合并行 (Hybrid)                                           │
│    └─ L1+L2 或 L1+L3 嵌套组合                                    │
└─────────────────────────────────────────────────────────────────┘
```

> **WASM 并行模型**：在浏览器环境中，Rem 编译为 WASM 模块，每个 Web Worker 加载一个 Rem WASM 实例作为一个 MPI rank。[jsmpi](https://github.com/jsmpi/jsmpi) 提供与 rsmpi 兼容的 API（`send`/`receive`/`broadcast`/`barrier`），底层通过 `SharedArrayBuffer` + `Atomics` 实现 Worker 间同步通信。Rem 的 MPI 代码在 Native 和 WASM 之间**零修改**切换——Native 编译链接 rsmpi（真实 MPI），WASM 编译链接 jsmpi（Web Worker 模拟）。

### 2.3 Rust MPI 集成方案

#### 方案对比

| 方案 | 库 | 优点 | 缺点 |
|------|-----|------|------|
| **A: rsmpi** | [rsmpi](https://github.com/rsmpi/rsmpi) | Rust 原生 MPI 绑定，类型安全，支持 MPI-3.1 | 需系统安装 MPI（OpenMPI/MPICH），WASM 不可用 |
| **B: jsmpi** | [jsmpi](https://github.com/jsmpi/jsmpi) | 浏览器 WASM 环境 MPI，与 rsmpi API 兼容 | 仅限浏览器，性能受 Web Worker 限制，集合操作尚不完整 |
| **C: Rem FFI** | Rem 自身的 MPI 层 | 直接复用 Rem 的 MPI 并行能力 | EMStudio 只做调度，不直接操作 MPI |
| **D: 多进程** | std::process + IPC | 无 MPI 依赖，简单可靠 | 无法做细粒度的网格域分解，WASM 不可用 |
| **E: 混合** | rsmpi/jsmpi + Rem FFI + 多进程 | 各取所长，Native/WASM 统一 API | 复杂度高 |

**推荐方案：E（混合）**，分阶段实施：

```
Phase 1: 多进程 + jsmpi 基础（L1 频率并行 + L3 任务并行）
  ├─ [Native] std::process 启动多个 Rem 求解进程，文件系统 IPC
  └─ [WASM]  jsmpi Coordinator 启动多个 Web Worker，
             每个 Worker 加载 Rem WASM 实例作为独立 rank
             Worker 间通过 SharedArrayBuffer + Atomics 通信

Phase 2: Rem MPI 域分解（L2 域分解）
  ├─ [Native] Rem 链接 rsmpi → 真实 MPI，mpirun 启动
  └─ [WASM]  Rem 链接 jsmpi → Web Worker MPI 模拟
             同一份 Rem 源码，编译目标不同，MPI 层自动切换

Phase 3: 高级 HPC（L4 混合并行 + 调度器）
  ├─ [Native] rsmpi 实现 EMStudio 自身的 MPI 协调层
  │           编排 SDM + DDM 多层嵌套，SLURM/PBS 集成
  └─ [WASM]  jsmpi 协调多层 Web Worker 嵌套
             Cloud 模式可将大规模作业提交到服务端
```

### 2.4 jsmpi 技术架构

[jsmpi](https://github.com/jsmpi/jsmpi) 是专为浏览器 WASM 设计的 Rust MPI 兼容层，让编译到 WASM 的 rsmpi 风格程序可以通过 Web Worker 并行执行。

#### 核心架构

```
┌──────────────────────────────────────────────────────────────┐
│                      浏览器主线程                              │
│  ┌──────────────────────────────────────────────────────┐    │
│  │              jsmpi Coordinator                        │    │
│  │  - 创建/管理 Web Workers                              │    │
│  │  - 分配 MPI rank                                     │    │
│  │  - 协调 Worker 生命周期                                │    │
│  └─────────┬─────────┬──────────┬──────────┬────────────┘    │
│            │         │          │          │                  │
│  ┌─────────▼───┐ ┌───▼───────┐ ┌▼────────┐ ┌▼─────────────┐ │
│  │ Web Worker 0│ │Web Worker 1│ │Worker 2 │ │  Worker N    │ │
│  │ (Rank 0)   │ │(Rank 1)   │ │(Rank 2) │ │ (Rank N)     │ │
│  │            │ │           │ │         │ │              │ │
│  │ ┌────────┐ │ │ ┌────────┐│ │┌───────┐│ │ ┌──────────┐ │ │
│  │ │Rem WASM│ │ │ │Rem WASM││ ││Rem    ││ │ │Rem WASM  │ │ │
│  │ │Instance│ │ │ │Instance││ ││WASM   ││ │ │Instance  │ │ │
│  │ └────────┘ │ │ └────────┘│ │└───────┘│ │ └──────────┘ │ │
│  └─────┬──────┘ └─────┬─────┘ └────┬────┘ └──────┬──────┘ │
│        │              │            │              │         │
│  ══════╪══════════════╪════════════╪══════════════╪═══════  │
│        │     SharedArrayBuffer + Atomics          │         │
│        │       (阻塞式同步 MPI 通信)                │         │
│  ══════╪══════════════╪════════════╪══════════════╪═══════  │
└────────┴──────────────┴────────────┴──────────────┴─────────┘
```

#### 关键技术点

| 技术 | 说明 |
|------|------|
| **SharedArrayBuffer** | Worker 间共享内存区域，用于零拷贝数据传递 |
| **Atomics.wait / notify** | 实现阻塞式 MPI 操作（send/receive 语义要求阻塞等待） |
| **API 兼容性** | 保持 `use mpi::traits::*` 导入模式，`world.rank()`、`send()`、`receive()` 调用路径与 rsmpi 一致 |
| **条件编译切换** | Rem 使用 `#[cfg(target_arch = "wasm32")]` 在 rsmpi 和 jsmpi 之间自动切换 |

#### 已支持的 MPI 操作

| 操作 | 状态 | 说明 |
|------|------|------|
| `initialize()` / Universe | ✅ | MPI 初始化，创建 Communicator |
| `rank()` / `size()` | ✅ | 查询当前 rank 和总 rank 数 |
| `send()` / `receive()` | ✅ | 点对点阻塞通信 |
| `broadcast_into()` | ✅ | 广播 |
| `barrier()` | ✅ | 全局同步屏障 |
| `gather` / `scatter` | 🔲 | 待实现（Rem DDM 需要） |
| `reduce` / `all_reduce` | 🔲 | 待实现（Rem DDM 需要） |
| 非阻塞操作 | 🔲 | 待实现 |

> **Rem 在 jsmpi 上的 DDM 可行性**：DDM 的核心 MPI 模式是「子域求解 → 边界交换（send/receive）→ 全局收敛检查（all_reduce）→ 迭代」。jsmpi 已支持 send/receive/barrier，待 all_reduce 实现后即可支持 DDM 迭代收敛。gather/scatter 用于初始网格分区分发和最终结果汇集，也需要实现。

#### Rem 条件编译示例

```rust
// rem/src/mpi_backend.rs
// 同一套代码，Native 用真实 MPI，WASM 用 jsmpi Web Worker 模拟

#[cfg(not(target_arch = "wasm32"))]
use rsmpi as mpi;

#[cfg(target_arch = "wasm32")]
use jsmpi as mpi;

use mpi::traits::*;

pub fn solve_subdomain(world: &mpi::topology::SimpleCommunicator) {
    let rank = world.rank();
    let size = world.size();

    if rank == 0 {
        // Master: 分发网格子域
        for dest in 1..size {
            let subdomain = partition_mesh(dest as usize);
            world.process_at_rank(dest).send(&subdomain);
        }
    } else {
        // Worker: 接收子域并求解
        let (subdomain, _) = world.any_process().receive::<SubdomainMesh>();
        let local_solution = solve_local_fem(&subdomain);
        world.process_at_rank(0).send(&local_solution);
    }
}
```

#### WASM 部署要求

| 要求 | 说明 |
|------|------|
| 浏览器支持 | 需要 `SharedArrayBuffer` 和 `Atomics` API（Chrome 68+, Firefox 79+, Safari 15.2+） |
| HTTP 头 | 服务器需设置 `Cross-Origin-Opener-Policy: same-origin` 和 `Cross-Origin-Embedder-Policy: require-corp` |
| Worker 数量 | 受浏览器 Web Worker 数量限制（通常 ≤ `navigator.hardwareConcurrency`） |
| 内存限制 | 每个 Worker 的 WASM 内存受浏览器限制（通常 2~4GB） |
| 性能预期 | Web Worker 通信延迟高于真实 MPI（~μs vs ~ns），适合粗粒度并行 |

### 2.5 HPC 配置数据结构

```rust
/// HPC 并行计算配置
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HpcConfig {
    /// 配置模式
    pub mode: HpcMode,
    /// MPI 实现选型
    pub mpi_vendor: MpiVendor,
    /// 计算资源
    pub resources: ComputeResources,
    /// 各并行层级开关与配置
    pub parallelism: ParallelismConfig,
    /// 作业调度器配置（可选，多机场景）
    pub scheduler: Option<SchedulerConfig>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum HpcMode {
    /// 自动模式：EMStudio 根据硬件自动决定并行策略
    Auto,
    /// 手动模式：用户显式配置
    Manual,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MpiVendor {
    /// 不使用 MPI（仅多线程/多进程）
    None,
    /// OpenMPI (Native Linux/macOS)
    OpenMPI,
    /// MPICH / Intel MPI (Native Linux/Windows)
    MPICH,
    /// jsmpi — Web Worker 模拟 MPI (WASM 专用)
    Jsmpi,
    /// 系统默认：Native 自动检测 OpenMPI/MPICH，WASM 自动使用 jsmpi
    SystemDefault,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComputeResources {
    /// 本地机器可用核数（0 = 自动检测）
    pub local_cores: u32,
    /// 每个 MPI 任务分配的核数
    pub cores_per_task: u32,
    /// 每个任务的最大内存 (MB)，0 = 不限制
    pub max_memory_per_task_mb: u64,
    /// 远程机器列表（多机场景）
    pub remote_hosts: Vec<RemoteHost>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RemoteHost {
    pub hostname: String,
    pub cores: u32,
    pub max_memory_mb: u64,
    /// SSH 端口（用于远程启动 MPI 进程）
    pub ssh_port: u16,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParallelismConfig {
    /// L0: 线程并行（Rem 内部 + rayon）
    pub thread_parallelism: ThreadConfig,
    /// L1: 频率并行（SDM）
    pub spectral_decomposition: Option<SdmConfig>,
    /// L2: 域分解（DDM）
    pub domain_decomposition: Option<DdmConfig>,
    /// L3: 任务并行（DSO）
    pub task_distribution: Option<DsoConfig>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThreadConfig {
    /// 线程数（0 = 自动，使用所有可用核）
    pub num_threads: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SdmConfig {
    pub enabled: bool,
    /// 同时并行的最大频率点数（0 = 自动）
    pub max_concurrent_frequencies: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DdmConfig {
    pub enabled: bool,
    /// 子域数量（0 = 自动，= num_tasks - 1）
    pub num_subdomains: u32,
    /// 子域求解器类型
    pub subdomain_solver: SubdomainSolverType,
    /// 全局迭代求解器收敛阈值
    pub coupling_convergence: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum SubdomainSolverType {
    /// 直接求解器（LU 分解），高鲁棒性，高内存
    Direct,
    /// 迭代求解器（GMRES），低内存
    Iterative,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DsoConfig {
    pub enabled: bool,
    /// 最大并发任务数（0 = 自动）
    pub max_concurrent_tasks: u32,
}

/// 作业调度器配置
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SchedulerConfig {
    pub scheduler_type: SchedulerType,
    /// 队列/分区名称
    pub queue: String,
    /// 作业最大运行时间
    pub walltime: String,
    /// 额外调度器参数
    pub extra_args: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SchedulerType {
    /// 本地直接运行（无调度器）
    Local,
    /// SLURM
    Slurm,
    /// PBS Professional
    Pbs,
    /// IBM Spectrum LSF
    Lsf,
    /// Sun Grid Engine
    Sge,
}
```

### 2.6 工程文件集成

在 Analysis Setup 中增加 HPC 配置（对应 [em-project-file-design.md §4.8](em-project-file-design.md) 的扩展）：

```json
{
  "analysis_setups": [
    {
      "name": "Setup1",
      "solution_frequency": "2.4GHz",
      "max_passes": 15,
      "max_delta_s": 0.02,
      "hpc": {
        "mode": "Auto",
        "mpi_vendor": "SystemDefault",
        "resources": {
          "local_cores": 0,
          "cores_per_task": 2,
          "max_memory_per_task_mb": 0,
          "remote_hosts": []
        },
        "parallelism": {
          "thread_parallelism": { "num_threads": 0 },
          "spectral_decomposition": { "enabled": true, "max_concurrent_frequencies": 0 },
          "domain_decomposition": null,
          "task_distribution": null
        },
        "scheduler": null
      }
    }
  ]
}
```

**Q3D Setup HPC 示例**（CG 分布式求解）：

```json
{
  "name": "Q3D_Setup1",
  "solution_type": "Q3D_ACRL",
  "hpc": {
    "mode": "Manual",
    "mpi_vendor": "OpenMPI",
    "resources": {
      "local_cores": 32,
      "cores_per_task": 4,
      "max_memory_per_task_mb": 16384,
      "remote_hosts": [
        { "hostname": "node02", "cores": 32, "max_memory_mb": 65536, "ssh_port": 22 },
        { "hostname": "node03", "cores": 32, "max_memory_mb": 65536, "ssh_port": 22 }
      ]
    },
    "parallelism": {
      "thread_parallelism": { "num_threads": 4 },
      "spectral_decomposition": { "enabled": true, "max_concurrent_frequencies": 4 },
      "domain_decomposition": {
        "enabled": true,
        "num_subdomains": 0,
        "subdomain_solver": "Direct",
        "coupling_convergence": 1e-4
      },
      "task_distribution": null
    },
    "scheduler": {
      "scheduler_type": "Slurm",
      "queue": "compute",
      "walltime": "4:00:00",
      "extra_args": ["--mem-per-cpu=4G"]
    }
  }
}
```

### 2.7 求解调度引擎设计

```rust
/// 求解调度器：根据 HPC 配置编排并行求解
pub struct SolveOrchestrator {
    hpc_config: HpcConfig,
    solver_backend: Arc<dyn SolverBackend>,
}

impl SolveOrchestrator {
    /// 执行完整分析（自适应 + 扫频）
    pub async fn run_analysis(
        &self,
        setup: &AnalysisSetup,
        progress: &dyn ProgressCallback,
    ) -> Result<AnalysisResult, SolveError> {
        // 1. 模型验证
        self.validate(setup)?;

        // 2. 自适应网格加密循环（单频点或多频点并行）
        let adaptive_result = self.run_adaptive_loop(setup, progress).await?;

        // 3. 频率扫描（SDM 频率并行 或 DSO 任务并行）
        if let Some(sweeps) = &setup.frequency_sweeps {
            for sweep in sweeps {
                self.run_frequency_sweep(setup, sweep, &adaptive_result, progress).await?;
            }
        }

        Ok(adaptive_result)
    }

    /// 频率扫描：根据配置选择并行策略
    async fn run_frequency_sweep(
        &self,
        setup: &AnalysisSetup,
        sweep: &FrequencySweep,
        mesh: &AdaptiveResult,
        progress: &dyn ProgressCallback,
    ) -> Result<SweepResult, SolveError> {
        let freq_points = sweep.generate_frequency_points();
        let parallelism = &self.hpc_config.parallelism;

        if let Some(sdm) = &parallelism.spectral_decomposition {
            if sdm.enabled {
                // 频率并行：将频率点分组，各组并行求解
                return self.solve_frequencies_parallel(freq_points, mesh, sdm, progress).await;
            }
        }

        // 串行求解
        self.solve_frequencies_sequential(freq_points, mesh, progress).await
    }

    /// SDM 频率并行求解
    async fn solve_frequencies_parallel(
        &self,
        freq_points: Vec<f64>,
        mesh: &AdaptiveResult,
        sdm: &SdmConfig,
        progress: &dyn ProgressCallback,
    ) -> Result<SweepResult, SolveError> {
        let max_concurrent = if sdm.max_concurrent_frequencies == 0 {
            self.hpc_config.resources.available_tasks()
        } else {
            sdm.max_concurrent_frequencies as usize
        };

        // 使用信号量控制并发度
        let semaphore = Arc::new(tokio::sync::Semaphore::new(max_concurrent));
        let mut handles = Vec::new();

        for freq in freq_points {
            let permit = semaphore.clone().acquire_owned().await?;
            let solver = self.solver_backend.clone();
            let mesh = mesh.clone();

            let handle = tokio::spawn(async move {
                let result = solver.solve_at_frequency(freq, &mesh).await;
                drop(permit);
                result
            });
            handles.push(handle);
        }

        // 收集所有结果
        let results = futures::future::try_join_all(handles).await?;
        Ok(SweepResult::from_frequency_results(results))
    }
}
```

### 2.8 多进程求解调度（Phase 1 — Native）

Phase 1 Native 端采用多进程方案，通过 `std::process::Command` 启动 Rem 求解进程：

```rust
/// 多进程求解器：通过子进程启动独立 Rem 实例
pub struct MultiProcessSolver {
    rem_binary_path: PathBuf,
    scratch_dir: PathBuf,
    max_processes: usize,
}

impl MultiProcessSolver {
    /// 并行求解多个频率点
    pub async fn solve_frequencies(
        &self,
        frequencies: &[f64],
        mesh_path: &Path,
        setup: &AnalysisSetup,
    ) -> Result<Vec<FrequencyResult>, SolveError> {
        let semaphore = Arc::new(tokio::sync::Semaphore::new(self.max_processes));
        let mut tasks = Vec::new();

        for (i, &freq) in frequencies.iter().enumerate() {
            let permit = semaphore.clone().acquire_owned().await.unwrap();
            let rem_bin = self.rem_binary_path.clone();
            let scratch = self.scratch_dir.join(format!("freq_{i:04}"));
            let mesh = mesh_path.to_path_buf();
            let setup_json = serde_json::to_string(setup)?;

            let task = tokio::spawn(async move {
                std::fs::create_dir_all(&scratch)?;

                // 写入求解配置
                let config_path = scratch.join("solve_config.json");
                std::fs::write(&config_path, &setup_json)?;

                // 启动 Rem 子进程
                let output = tokio::process::Command::new(&rem_bin)
                    .args(["solve", "--mesh", mesh.to_str().unwrap()])
                    .args(["--frequency", &freq.to_string()])
                    .args(["--config", config_path.to_str().unwrap()])
                    .args(["--output-dir", scratch.to_str().unwrap()])
                    .output()
                    .await?;

                drop(permit);

                if !output.status.success() {
                    return Err(SolveError::RemFailed(
                        String::from_utf8_lossy(&output.stderr).to_string()
                    ));
                }

                // 读取结果
                let result_path = scratch.join("result.json");
                let result: FrequencyResult = serde_json::from_str(
                    &std::fs::read_to_string(result_path)?
                )?;
                Ok(result)
            });

            tasks.push(task);
        }

        let results: Vec<Result<FrequencyResult, SolveError>> =
            futures::future::join_all(tasks)
                .await
                .into_iter()
                .map(|r| r.unwrap())
                .collect();

        results.into_iter().collect()
    }
}
```

### 2.9 Web Worker 求解调度（Phase 1 — WASM）

Phase 1 WASM 端采用 jsmpi Web Worker 方案，每个 Worker 加载一个 Rem WASM 实例：

```rust
/// WASM Web Worker 求解器：通过 jsmpi 在浏览器中并行求解
#[cfg(target_arch = "wasm32")]
pub struct WebWorkerSolver {
    /// Rem WASM 模块 URL（预编译的 .wasm 文件）
    rem_wasm_url: String,
    /// 最大 Worker 数量
    max_workers: usize,
}

#[cfg(target_arch = "wasm32")]
impl WebWorkerSolver {
    pub fn new(rem_wasm_url: &str) -> Self {
        // 自动检测可用核数
        let max_workers = web_sys::window()
            .and_then(|w| w.navigator().hardware_concurrency().try_into().ok())
            .unwrap_or(4);

        Self {
            rem_wasm_url: rem_wasm_url.to_string(),
            max_workers,
        }
    }

    /// 并行求解多个频率点（SDM 模式）
    pub async fn solve_frequencies(
        &self,
        frequencies: &[f64],
        mesh_data: &[u8],       // 序列化的网格数据
        setup_json: &str,
    ) -> Result<Vec<FrequencyResult>, SolveError> {
        // 通过 jsmpi Coordinator 创建 Worker 池
        let num_workers = frequencies.len().min(self.max_workers);

        // jsmpi 初始化：创建 SharedArrayBuffer 通信通道
        let coordinator = jsmpi::Coordinator::new(num_workers)?;

        // 将频率点分组分配给各 Worker
        let chunks = distribute_evenly(frequencies, num_workers);

        let mut results = Vec::new();
        for (rank, freq_chunk) in chunks.iter().enumerate() {
            // 每个 Worker 加载 Rem WASM 模块并求解分配到的频率点
            let worker_result = coordinator.spawn_task(rank, |world| {
                // 这段代码在 Web Worker 中执行
                // world 是 jsmpi 提供的 MPI Communicator
                let mut freq_results = Vec::new();
                for &freq in freq_chunk {
                    let result = rem_wasm::solve_at_frequency(
                        world, mesh_data, freq, setup_json
                    );
                    freq_results.push(result);
                }
                freq_results
            }).await?;

            results.extend(worker_result);
        }

        Ok(results)
    }

    /// DDM 域分解求解（Phase 2 — 需要 jsmpi 完成 gather/scatter/all_reduce）
    pub async fn solve_ddm(
        &self,
        mesh_data: &[u8],
        frequency: f64,
        num_subdomains: usize,
        setup_json: &str,
    ) -> Result<SolveResult, SolveError> {
        // 需要 num_subdomains + 1 个 Worker（1 Master + N 子域）
        let num_workers = num_subdomains + 1;
        let coordinator = jsmpi::Coordinator::new(num_workers)?;

        // 所有 Worker 运行相同的 Rem DDM 代码
        // jsmpi 提供 MPI 语义，Rem 内部用 send/receive/all_reduce 协调
        let result = coordinator.run_all(|world| {
            // 这段 Rem 代码与 Native MPI 版本完全相同
            // #[cfg] 条件编译在此处已选择 jsmpi 作为 mpi 后端
            rem_wasm::solve_ddm(world, mesh_data, frequency, setup_json)
        }).await?;

        Ok(result)
    }
}
```

**Native vs WASM 统一调度接口**：

```rust
/// 平台无关的求解接口
pub trait ParallelSolver {
    async fn solve_frequencies(
        &self,
        frequencies: &[f64],
        mesh: &MeshData,
        setup: &AnalysisSetup,
    ) -> Result<Vec<FrequencyResult>, SolveError>;
}

#[cfg(not(target_arch = "wasm32"))]
impl ParallelSolver for MultiProcessSolver { /* ... */ }

#[cfg(target_arch = "wasm32")]
impl ParallelSolver for WebWorkerSolver { /* ... */ }

/// 平台自动选择
pub fn create_solver(config: &HpcConfig) -> Box<dyn ParallelSolver> {
    #[cfg(not(target_arch = "wasm32"))]
    {
        Box::new(MultiProcessSolver::from_config(config))
    }
    #[cfg(target_arch = "wasm32")]
    {
        Box::new(WebWorkerSolver::new(&config.rem_wasm_url()))
    }
}
```

### 2.10 MPI 域分解接口（Phase 2 — Native rsmpi / WASM jsmpi）

```rust
/// MPI 域分解求解器：通过 Rem 的分布式求解接口
/// Native 用 mpirun 启动，WASM 用 jsmpi Coordinator 启动
pub struct MpiDdmSolver {
    config: DdmConfig,
    num_tasks: u32,
    cores_per_task: u32,
    hosts: Vec<String>,
}

impl MpiDdmSolver {
    /// 启动 MPI 域分解求解
    pub async fn solve(
        &self,
        mesh_path: &Path,
        frequency: f64,
        setup: &AnalysisSetup,
        output_dir: &Path,
    ) -> Result<SolveResult, SolveError> {
        let hostfile = self.write_hostfile()?;
        let total_tasks = self.num_tasks;

        // 通过 mpirun/mpiexec 启动分布式 Rem
        let output = tokio::process::Command::new("mpirun")
            .args(["-np", &total_tasks.to_string()])
            .args(["--hostfile", hostfile.to_str().unwrap()])
            .args(["--map-by", &format!("slot:PE={}", self.cores_per_task)])
            .arg("rem-solver")
            .args(["--mode", "ddm"])
            .args(["--mesh", mesh_path.to_str().unwrap()])
            .args(["--frequency", &frequency.to_string()])
            .args(["--num-subdomains", &(total_tasks - 1).to_string()])
            .args(["--convergence", &self.config.coupling_convergence.to_string()])
            .args(["--output-dir", output_dir.to_str().unwrap()])
            .output()
            .await?;

        if !output.status.success() {
            return Err(SolveError::MpiLaunchFailed(
                String::from_utf8_lossy(&output.stderr).to_string()
            ));
        }

        let result_path = output_dir.join("ddm_result.json");
        let result: SolveResult = serde_json::from_str(
            &std::fs::read_to_string(result_path)?
        )?;

        Ok(result)
    }

    fn write_hostfile(&self) -> Result<PathBuf, SolveError> {
        let path = std::env::temp_dir().join("emstudio_hostfile");
        let content: String = self.hosts.iter()
            .map(|h| format!("{} slots={}", h, self.cores_per_task))
            .collect::<Vec<_>>()
            .join("\n");
        std::fs::write(&path, content)?;
        Ok(path)
    }
}
```

### 2.11 Q3D CG 分布式求解设计

```rust
/// Q3D CG 分布式矩阵求解器
pub struct Q3dDistributedCgSolver {
    config: DdmConfig,
    num_tasks: u32,
    cores_per_task: u32,  // 推荐 4
}

impl Q3dDistributedCgSolver {
    /// 分布式 CG 矩阵组装与求解
    pub async fn solve_cg(
        &self,
        surface_mesh_path: &Path,
        nets: &[Net],
        frequency: f64,
        output_dir: &Path,
    ) -> Result<RlcgMatrix, SolveError> {
        // Q3D CG 分布式流程：
        // 1. Master 进行表面网格分区（按三角形行分块）
        // 2. 各 rank 并行计算 Green 函数矩阵的局部块
        // 3. 并行求解线性系统 → 电荷分布
        // 4. 从电荷分布积分得到电容矩阵
        // 5. Master 汇总 C/G 矩阵

        let total_tasks = self.num_tasks;

        let output = tokio::process::Command::new("mpirun")
            .args(["-np", &total_tasks.to_string()])
            .args(["--map-by", &format!("slot:PE={}", self.cores_per_task)])
            .arg("rem-solver")
            .args(["--mode", "q3d-cg-distributed"])
            .args(["--mesh", surface_mesh_path.to_str().unwrap()])
            .args(["--frequency", &frequency.to_string()])
            .args(["--output-dir", output_dir.to_str().unwrap()])
            .output()
            .await?;

        if !output.status.success() {
            return Err(SolveError::MpiLaunchFailed(
                String::from_utf8_lossy(&output.stderr).to_string()
            ));
        }

        let matrix_path = output_dir.join("rlcg_matrix.json");
        let matrix: RlcgMatrix = serde_json::from_str(
            &std::fs::read_to_string(matrix_path)?
        )?;
        Ok(matrix)
    }
}
```

---

## 3. HFSS vs EMStudio 并行方案对比

| 特性 | Ansys HFSS/Q3D | EMStudio (Native) | EMStudio (WASM) |
|------|---------------|-------------------|-----------------|
| 线程并行 | 内置（矩阵库） | Rem 内部 + rayon | Web Worker (jsmpi) |
| 频率并行 (SDM) | 内置 | Phase 1: 多进程 | Phase 1: Web Worker |
| 域分解 (DDM) | 内置 MPI | Phase 2: Rem + rsmpi | Phase 2: Rem + jsmpi |
| 任务并行 (DSO) | 内置 | Phase 1: 多进程 | Phase 1: Web Worker |
| 多层嵌套 | 支持两层 | Phase 3: rsmpi 协调 | Phase 3: jsmpi 协调 |
| MPI 实现 | Intel MPI / MS-MPI | OpenMPI / MPICH | jsmpi (Web Worker) |
| MPI 通信基础 | TCP / InfiniBand | TCP / InfiniBand | SharedArrayBuffer + Atomics |
| 自动模式 | 支持 | Phase 1 起支持 | Phase 1 起支持 |
| 手动模式 | 支持 | Phase 1 起支持 | Phase 1 起支持 |
| 调度器集成 | SLURM/PBS/LSF/SGE | Phase 3: SLURM 优先 | N/A（浏览器本地） |
| Q3D CG 分布式 | 支持 | Phase 2 | Phase 2 |
| Rem 代码共享 | — | ✅ 相同源码 | ✅ 相同源码（条件编译） |

---

## 4. 实现路线

### Phase 1：多进程 + jsmpi 基础并行

**目标**：通过多进程（Native）和 Web Worker（WASM）实现频率并行和任务并行。

| 任务 | 模块 | 平台 | 优先级 |
|------|------|------|--------|
| HpcConfig 数据结构 | `domain` | 通用 | P0 |
| HPC 配置 JSON 序列化/反序列化 | `domain` | 通用 | P0 |
| HPC 配置 UI 面板 | `components` | 通用 | P1 |
| MultiProcessSolver 多进程调度 | `solver` | Native | P0 |
| SolveOrchestrator 求解编排 | `solver` | 通用 | P0 |
| 频率并行 — SDM 多进程版 | `solver` | Native | P0 |
| 频率并行 — SDM Web Worker 版 (jsmpi) | `solver` | WASM | P0 |
| jsmpi Coordinator 集成到 EMStudio | `infra` | WASM | P0 |
| Rem WASM 模块加载 + Worker 实例化 | `infra` | WASM | P0 |
| 参数扫描并行 — DSO 多进程版 | `solver` | Native | P1 |
| 参数扫描并行 — DSO Web Worker 版 | `solver` | WASM | P1 |
| 自动模式（核数检测：Native=sysinfo / WASM=hardwareConcurrency） | `solver` | 通用 | P1 |
| 进度回调 + 实时收敛显示 | `solver` → `app` | 通用 | P1 |

### Phase 2：Rem MPI 域分解（Native rsmpi + WASM jsmpi）

**目标**：实现单问题域分解，突破单节点内存限制。Native 通过真实 MPI，WASM 通过 jsmpi Web Worker。

| 任务 | 模块 | 平台 | 优先级 |
|------|------|------|--------|
| Rem rsmpi 接口适配验证 | `solver` | Native | P0 |
| Rem jsmpi 接口适配验证 | `solver` | WASM | P0 |
| jsmpi gather/scatter/all_reduce 实现推进 | 上游 jsmpi | WASM | P0（前置依赖） |
| MpiDdmSolver 域分解求解器 | `solver` | 通用 | P0 |
| Q3dDistributedCgSolver CG 分布式 | `solver` | 通用 | P0 |
| mpirun/mpiexec 启动管理 | `solver` | Native | P0 |
| jsmpi Worker 池管理 + rank 分配 | `infra` | WASM | P0 |
| Hostfile 自动生成 | `solver` | Native | P1 |
| DDM 子域数自动优化 | `solver` | 通用 | P1 |
| 分布式内存直接求解器 | `solver` | Native | P2 |
| FE-BI 混合求解支持 | `solver` | Native | P2 |

### Phase 3：高级 HPC（嵌套并行 + 调度器 + Cloud）

**目标**：支持多层嵌套并行、集群调度器集成、WASM Cloud 模式。

| 任务 | 模块 | 平台 | 优先级 |
|------|------|------|--------|
| rsmpi 协调层（EMStudio 自身 MPI 编排） | `solver` | Native | P0 |
| jsmpi 协调层（多层 Worker 嵌套编排） | `solver` | WASM | P0 |
| SDM + DDM 两层嵌套编排 | `solver` | 通用 | P0 |
| SLURM 作业提交 | `infra` | Native | P1 |
| PBS 作业提交 | `infra` | Native | P2 |
| LSF 作业提交 | `infra` | Native | P2 |
| 远程节点状态监控 | `infra` | Native | P1 |
| 集群资源自动发现 | `infra` | Native | P2 |
| WASM Cloud 模式：大规模作业提交到服务端 | `infra` | WASM | P1 |
| HPC 性能剖面与调优建议 | `solver` | 通用 | P2 |

---

## 5. 参考资料

### Ansys HFSS HPC

- [HFSS High Performance Computing](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/HPC/HighPerformanceComputing.htm)
- [Domain Decomposition Method](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/HFSS/Content/HFSS/DomainDecompositionMethod.htm)
- [Distributed Analysis in HFSS](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/HPC/DistributedAnalysis.htm)
- [DSO Behavior in HFSS](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/HPC/DSOBehaviorinHFSSandHFSSIE.htm)
- [DDM HPC Configuration Guidelines](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/HPC/DomainDecompositionSolverGuidelinesforHPCConfiguration.htm)
- [Two-Level Distribution Guidelines](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/HFSS/Content/HPC/2LevelDistributionsGuidelines.htm)
- [Distributed Memory Solutions with HFSS](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/HFSS/DistributedMemorySolutionswithHFSS.htm)
- [Enable Domain Decomposition](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/HFSS/EnableDomainDecomposition.htm)
- [Direct Matrix Solver](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/HFSS/DirectMatrixSolver.htm)
- [Iterative Matrix Solver Technical Details](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/HFSS3DLayout/Content/HFSS/IterativeMatrixSolverTechnicalDetails.htm)
- [MLFMM Usage Guidelines](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/HFSS/MLFMMUsageGuidelines.htm)
- [Installation Requirements for Distributed Memory Solutions](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/HFSS/InstallationRequirementsforDistributedMemorySolutionswithHFSS.htm)
- [Setting HPC and Analysis Options](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/Variables/SettingHPCandAnalysisOptions.htm)
- [Ansys Blog: Optimize HFSS with HPC](https://www.ansys.com/zh-tw/blog/how-to-optimize-speed-scalability-ansys-hfss-hpc)
- [Ansys Blog: HFSS HPC on Ansys Cloud](https://www.ansys.com/zh-cn/blog/ansys-hfss-siwave-hpc-capabilities-ansys-cloud)
- [Webinar: Foundation of DDM in HFSS](https://www.ansys.com/resource-center/webinar/the-foundation-of-domain-decomposition-technologies-in-ansys-hfss)

### Ansys Q3D HPC

- [Q3D High Performance Computing](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/HPC/HighPerformanceComputing.htm)
- [Q3D Distributed Memory CG Solutions](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/DistributedMemoryCGSolutions.htm)
- [Q3D Field Simulation Methods](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Q3DExtractor/Content/Q3D/FieldSimulationMethods.htm)
- [Q3D Setting HPC Options](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/Q3DExtractor/Content/Variables/SettingHPCandAnalysisOptions.htm)
- [Q3D Windows to Linux Job Submission](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Q3DExtractor/Content/HPC/WindowstoLinuxJobSubmission.htm)

### MPI 与调度器

- [Ansys 2025 R1 MPI Support Matrix](https://www.ansys.com/content/dam/release/2025-r1/platform-support/ansys-2025-r1-message-passing-interface-support-for-parallel-computing.pdf)
- [Ansys 2024 R2 MPI Support Matrix](https://www.ansys.com/content/dam/it-solutions/platform-support/2024-r2/ansys-2024-r2-message-passing-interface-support-for-parallel-computing.pdf)
- [AEDT SLURM Integration](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v251/en/Subsystems/HFSS/Content/HPC/IntegrationwithSLURMLinuxScheduler.htm)
- [Ansys RSM User's Guide](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/pdf/Ansys_Remote_Solve_Manager_Users_Guide.pdf)
- [Ansys EM Install Guide (Windows)](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/PDFs/AnsysEMInstallGuide-Windows.pdf)
- [HPC Administrator's Guide](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/PDFs/HPC_Admin.pdf)
- [AEDT SelectScheduler API](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Maxwell/Subsystems/Maxwell%20Scripting/Content/SelectScheduler.htm)

### Rust 生态

- [rsmpi — Rust MPI bindings](https://github.com/rsmpi/rsmpi) — MPI-3.1 Rust 绑定（Native）
- [jsmpi — Browser MPI via Web Workers](https://github.com/jsmpi/jsmpi) — WASM 环境 MPI 兼容层，rsmpi API 兼容，SharedArrayBuffer + Atomics 实现
- [rayon](https://docs.rs/rayon/) — Rust 数据并行库
- [tokio](https://docs.rs/tokio/) — Rust 异步运行时

### Web 平台

- [SharedArrayBuffer — MDN](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/SharedArrayBuffer) — Web Worker 共享内存
- [Atomics — MDN](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Atomics) — 原子操作（wait/notify 实现阻塞语义）
- [Web Workers API — MDN](https://developer.mozilla.org/en-US/docs/Web/API/Web_Workers_API) — 浏览器多线程
- [Cross-Origin Isolation — web.dev](https://web.dev/articles/cross-origin-isolation-guide) — SharedArrayBuffer 所需的 COOP/COEP HTTP 头配置
