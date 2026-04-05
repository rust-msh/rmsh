# EmStudio Architecture Exploration

## Overview
EmStudio is a Rust-based electromagnetic simulation and visualization application built with egui, wgpu, and a modular architecture. This document explores the render, domain, and infrastructure layers to understand available capabilities for UI layout integration.

## Technology Stack

### Core Dependencies (Cargo.toml workspace)
- **UI Framework**: `egui` 0.33 + `eframe` 0.33 (with wgpu backend)
- **GPU Rendering**: `wgpu` 27 with `egui-wgpu` 0.33
- **Math**: `glam` 0.29 (vectors, matrices)
- **Data Serialization**: `serde` 1.0 + `rmp-serde` (MessagePack format)
- **UI Layout**: `egui_dock` 0.18 (docking system)
- **Tables/Grids**: `egui_table` 0.7, `egui_tiles` 0.14
- **Plotting**: `egui_plot` 0.33
- **File Dialogs**: `rfd` 0.15

---

## 1. Domain Layer (`crates/domain/src/lib.rs`)

### Data Structures

#### SimulationStatus Enum
```rust
pub enum SimulationStatus {
    Idle,
    Solving,
    Finished,
    Failed,
}
```

#### Material
```rust
pub struct Material {
    pub name: String,
    pub relative_permittivity: f32,
    pub conductivity: f32,
}
```
- Default: "Vacuum" material (ε_r=1.0, σ=0.0)

#### GeometryObject
```rust
pub struct GeometryObject {
    pub id: u64,
    pub name: String,
    pub mesh_hint: String,
}
```
- Represents geometry entities in the model
- `mesh_hint`: e.g., "auto", "fine", "coarse"

#### EmModel
```rust
pub struct EmModel {
    pub name: String,
    pub objects: Vec<GeometryObject>,
    pub materials: Vec<Material>,
}
```
- Default: "Untitled Model" with 0 objects and 1 Vacuum material

#### SolveResult
```rust
pub struct SolveResult {
    pub field_preview: String,
    pub converged: bool,
}
```
- Contains solution convergence status
- `field_preview`: Describes the field data (currently a string placeholder)

#### Project
```rust
pub struct Project {
    pub id: String,
    pub title: String,
    pub model: EmModel,
    pub status: SimulationStatus,
    pub last_result: Option<SolveResult>,
}
```
- Default project: "New Project" with Idle status

### Key Capabilities
✅ Full serialization support (JSON and MessagePack)
✅ Domain models are clean and ready for UI binding
✅ No GPU/rendering concerns in this layer (proper separation)

---

## 2. Infrastructure Layer (`crates/infra/src/lib.rs`)

### Run Modes
```rust
pub enum RunMode {
    Standalone,    // In-memory project storage
    Cloud,         // Remote API endpoint
}
```

### Backend Trait
```rust
pub trait Backend {
    fn save_project(&mut self, project: Project) -> Result<(), BackendError>;
    fn load_project(&self, id: &str) -> Result<Project, BackendError>;
    fn solve(&self, project: &Project) -> Result<SolveResult, BackendError>;
    fn mode(&self) -> RunMode;
}
```

### Implementations
1. **StandaloneBackend**: In-memory HashMap-based project storage
2. **CloudBackend**: Placeholder for remote solving via HTTP endpoint

### File I/O
- Format: MessagePack (.emsp files)
- Functions:
  - `save_project_to_file(project: &Project, path: &Path)`
  - `load_project_from_file(path: &Path)`

### Error Handling
```rust
pub enum BackendError {
    ProjectNotFound(String),
    IoError(String),
    SerializeError(String),
    DeserializeError(String),
}
```

### Key Capabilities
✅ Async file I/O abstraction
✅ Swappable backends (standalone vs cloud)
✅ Solver abstraction (currently PlaceholderSolver)
✅ Proper error propagation

---

## 3. Render Layer (`crates/render/src/`)

### 3.1 SceneViewport (Main Viewport Container)
**Location**: `lib.rs`

```rust
pub struct SceneViewport {
    pub title: String,
    pub config: WgpuRenderConfig,
    runtime_status: RuntimeStatus,
    render_state: Option<OffscreenRenderer>,
    frame_counter: u64,
}
```

**Method**: `pub fn ui(&mut self, ui: &mut Ui, project: &Project)`
- Main integration point for rendering in egui UI
- Allocates space with `ui.allocate_exact_size()`
- Uses egui painter for rendering status overlays
- Supports drag interaction

**Key Capabilities**:
- ✅ Offscreen WGPU rendering
- ✅ GPU resource initialization (async on native, deferred on WASM)
- ✅ Frame counting and status reporting
- ✅ Drag input handling
- ✅ Repaint requests for animation

### 3.2 FieldSceneState (Advanced 3D Visualization)
**Location**: `scene.rs`

```rust
pub struct FieldSceneState {
    pub camera: OrbitCamera,
    pub colormap: ColormapType,
    pub opacity: f32,
    pub show_wireframe: bool,
    pub field_range: [f32; 2],
    pub vis_mode: VisMode,
    pub animator: PhaseAnimator,
    pub slice_z: f32,
    // ... internal GPU resources
}
```

**Public Methods**:
1. `init_gpu(&mut self, render_state: &egui_wgpu::RenderState, mesh: &FieldMesh)`
   - Initializes all GPU pipelines and resources
   - Called once from `App::new()`

2. `show_viewport(&mut self, ui: &mut egui::Ui)`
   - Renders the 3D field visualization
   - Handles mouse interaction (rotate, pan, zoom)
   - Implements mode switching logic
   - Supports phase animation playback

3. `show_controls(&mut self, ui: &mut egui::Ui)`
   - UI control panel for visualization parameters
   - Mode selector (Surface, Arrows, Slice, FarField, Animation)
   - Colormap picker
   - Opacity slider
   - Wireframe toggle
   - View preset buttons (Front, Back, Left, Right, Top, Bottom, Iso)
   - Phase animation controls

4. `show_colorbar(&self, ui: &mut egui::Ui)`
   - Vertical colorbar legend with field range annotations

**Visualization Modes** (VisMode enum):
```rust
pub enum VisMode {
    Surface,      // Colormap on UV sphere
    Arrows,       // Vector arrows on cube surface
    Slice,        // Slice plane through volume
    FarField,     // 3D radiation pattern
    Animation,    // Phase-animated complex field
}
```

**Key Capabilities**:
- ✅ 5 distinct visualization modes
- ✅ Orbit camera (azimuth, elevation, distance, zoom)
- ✅ 4 colormaps (Rainbow, Viridis, Cool-Warm, Grayscale)
- ✅ Complex field animation (phase playback at 0-360°)
- ✅ Wireframe overlay rendering
- ✅ Opacity control (0-1 range)
- ✅ Pre-generated meshes (sphere, cube, far-field, slice)

### 3.3 OrbitCamera
**Location**: `camera.rs`

```rust
pub struct OrbitCamera {
    pub target: Vec3,
    pub distance: f32,      // 0.1 to 100.0
    pub azimuth: f32,       // radians
    pub elevation: f32,     // radians, clamped ±π/2
    pub fov_y: f32,         // default: π/4
    pub near: f32,          // default: 0.01
    pub far: f32,           // default: 100.0
}
```

**Methods**:
- `rotate(dx: f32, dy: f32)` - Mouse drag rotation (sensitivity 0.005)
- `zoom(delta: f32)` - Mouse wheel zoom (exponential, 0.1x-10x)
- `pan(dx: f32, dy: f32)` - Screen-space panning
- `eye_position() -> Vec3` - Compute camera position
- `view_matrix() -> Mat4` - Compute view matrix (RH)
- `projection_matrix(aspect: f32) -> Mat4` - Compute projection (RH, perspective)
- `view_projection(aspect: f32) -> Mat4` - Combined matrix
- `set_preset(preset: ViewPreset)` - Apply predefined view angles

**View Presets**:
```rust
pub enum ViewPreset {
    Front,   // azimuth: 0°, elevation: 0°
    Back,    // azimuth: 180°, elevation: 0°
    Left,    // azimuth: -90°, elevation: 0°
    Right,   // azimuth: 90°, elevation: 0°
    Top,     // azimuth: 0°, elevation: 90°
    Bottom,  // azimuth: 0°, elevation: -90°
    Iso,     // azimuth: 34.3°, elevation: 23.0°
}
```

**Key Capabilities**:
- ✅ Full 3D orbit camera control
- ✅ Right-handed coordinate system (Y-up)
- ✅ Parameterized FOV, near/far planes
- ✅ 7 preset view angles
- ✅ Smooth zooming with exponential scaling

### 3.4 Mesh Data Structures
**Location**: `mesh_data.rs`

#### FieldVertex (GPU Vertex Format)
```rust
#[repr(C)]
pub struct FieldVertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
    pub field_value: f32,
}
```
- Corresponds to WGSL vertex shader inputs
- Properly aligned for GPU buffer uploads

#### ArrowInstance (Instanced Rendering)
```rust
#[repr(C)]
pub struct ArrowInstance {
    pub position: [f32; 3],
    pub direction: [f32; 3],
    pub magnitude: f32,
    pub _pad: f32,
}
```
- Per-instance data for arrow vector visualization
- Allows efficient rendering of 1000s of arrows

#### FieldMesh
```rust
pub struct FieldMesh {
    pub vertices: Vec<FieldVertex>,
    pub indices: Vec<u32>,
    pub wire_indices: Vec<u32>,
    pub field_range: [f32; 2],
    pub field_imag: Option<Vec<f32>>,        // For phase animation
    pub vector_field: Option<Vec<[f32; 3]>>, // For arrow visualization
}
```

**Static Methods**:
1. `uv_sphere(n_lat: u32, n_lon: u32, radius: f32) -> Self`
   - Generates UV sphere with synthetic spherical harmonic field
   - Supports complex field data (real + imaginary components)
   - Used for Surface and Animation modes

2. `cube(subdivisions: u32, half_size: f32) -> Self`
   - Generates subdivided cube with rotating vortex vector field
   - 6 faces, each subdivided into quads
   - Used for Arrows mode

3. `generate_arrows(&self, every_n: u32) -> Vec<ArrowInstance>`
   - Subsamples vector field data into arrow instances
   - Normalizes magnitudes for visual scaling

**Key Capabilities**:
- ✅ Multiple mesh topologies (sphere, cube, far-field, slice planes)
- ✅ Real and imaginary field component storage
- ✅ Vector field data for arrow visualization
- ✅ Field range computation
- ✅ Wireframe index generation

### 3.5 Mesh Generation Modules

#### Animation (`animation.rs`)
```rust
pub struct PhaseAnimator {
    pub phase_deg: f32,           // 0-360°
    pub playing: bool,
    pub speed_deg_per_sec: f32,   // 10-720° typical
}
```

**Methods**:
- `tick(dt: f32)` - Advance phase based on delta time
- `apply(field_real: &[f32], field_imag: &[f32]) -> Vec<f32>`
  - Computes: E(t) = Re(E)·cos(φ) - Im(E)·sin(φ)
- `envelope_range(field_real: &[f32], field_imag: &[f32]) -> [f32; 2]`
  - Conservative field range over all phases

**Key Capabilities**:
- ✅ Real-time phase sweeping at configurable speeds
- ✅ Proper Euler formula application for complex fields
- ✅ Field range scaling for all phases

#### Colormap (`colormap.rs`)
```rust
pub enum ColormapType {
    Rainbow,    // HSV sweep (blue → red)
    Viridis,    // Perceptually uniform (purple → yellow)
    CoolWarm,   // Diverging (blue → gray → red)
    Grayscale,  // Linear intensity
}
```

**Methods**:
- `generate_lut(n: usize) -> Vec<[u8; 4]>`
  - Returns RGBA8 lookup table with n entries
  - Typical: n=256 for GPU texture

**Key Capabilities**:
- ✅ 4 professional colormaps
- ✅ Piecewise-linear interpolation for Viridis and Cool-Warm
- ✅ HSV color space for Rainbow

#### Far-Field (`far_field.rs`)
- Generates 3D radiation pattern meshes
- Default: 60×120 latitude/longitude resolution
- Patch antenna gain pattern synthesis

#### Slice (`slice.rs`)
```rust
pub enum SliceAxis {
    X, Y, Z
}
```
- Generates 2D slice planes through volume
- Synthetic volume field sampling at configurable Z position

### 3.6 FieldPipeline (GPU Rendering Backend)
**Location**: `field_pipeline.rs`

```rust
pub struct FieldPipeline {
    // Scene rendering
    scene_pipeline: wgpu::RenderPipeline,
    wire_pipeline: wgpu::RenderPipeline,
    // Vertex/index buffers
    vertex_buf: wgpu::Buffer,
    index_buf: wgpu::Buffer,
    wire_index_buf: wgpu::Buffer,
    // Uniforms
    uniform_buf: wgpu::Buffer,
    scene_bind_group: wgpu::BindGroup,
    // Colormap texture
    colormap_texture: wgpu::Texture,
    colormap_sampler: wgpu::Sampler,
    // Offscreen framebuffer (with depth)
    color_texture: wgpu::Texture,
    depth_view: wgpu::TextureView,
    fb_size: [u32; 2],
    // Blit to screen
    blit_pipeline: wgpu::RenderPipeline,
    // ...
}
```

**Uniform Buffer**:
```rust
#[repr(C)]
pub struct FieldUniforms {
    pub mvp: [f32; 16],         // Model-view-projection matrix
    pub eye_pos: [f32; 3],      // Camera position (for lighting)
    pub _pad0: f32,
    pub light_dir: [f32; 3],    // Directional light direction
    pub _pad1: f32,
    pub field_min: f32,         // Field range min (for colormap)
    pub field_max: f32,         // Field range max
    pub opacity: f32,           // Blend opacity (0-1)
    pub _pad2: f32,
}
```

**Methods**:
- `new(device, queue, target_format, mesh, colormap)` - Constructor
- `resize_if_needed(device, size)` - Handle viewport resize
- `update_uniforms(queue, uniforms)` - Update per-frame data
- `render_scene(encoder, show_wireframe, arrow_pipeline)` - Render to offscreen buffer
- `blit(render_pass)` - Copy result to egui's render pass
- `swap_mesh(device, mesh)` - Change displayed mesh (mode switching)
- `update_colormap(device, queue, colormap)` - Update colormap texture
- `update_vertices(queue, verts)` - Update vertex data (phase animation)

**Key Capabilities**:
- ✅ Offscreen rendering with depth (prevents Z-order issues with egui)
- ✅ Colormap-based field visualization
- ✅ Dynamic mesh switching
- ✅ Wireframe overlay rendering
- ✅ Dual rendering pipelines (solid + wireframe)
- ✅ Real-time vertex updates (animation)

### 3.7 Arrow Pipeline
**Location**: `arrow_pipeline.rs`

Implements instanced rendering of vector field arrows:
- Supports up to 4096 arrow instances
- Per-instance position, direction, magnitude
- Automatic scaling based on field range

---

## 4. App Layer (`crates/app/src/lib.rs`)

### Project Architecture

#### App Structure
```rust
pub struct App {
    project: Project,                           // Domain model
    backend: Box<dyn Backend>,                  // Infra layer
    viewport: SceneViewport,                    // Simple render viewport
    dock_state: DockState<CenterTab>,           // Layout management
    ribbon_state: RibbonState,                  // Toolbar state
    ribbon_tabs: Vec<RibbonTab>,                // Toolbar configuration
    current_file: Option<PathBuf>,
    unsaved_changes: bool,
    status_text: String,
    log_text: String,
    file_dialog_rx/tx: mpsc::Receiver/Sender,   // Async dialog results
}
```

#### Layout Tabs
```rust
enum CenterTab {
    Modeling,   // Contains SceneViewport
    Result,     // Shows solve result preview
    Log,        // Application log viewer
}
```

**Layout Structure** (using egui_dock):
```
┌─────────────────────────────────┐
│    Ribbon (tabs + buttons)       │
├──────────────────────┬──────────┤
│  Modeling (Viewport) │ Result   │
│  (main 3D view)      │ (65% W)  │
├──────────────────────┤          │
│     Log              │          │
│     (50% H)          │          │
└──────────────────────┴──────────┘
```

#### Ribbon System
From `crates/components/src/ribbon.rs`:

```rust
pub enum RibbonAction {
    // Desktop
    NewProject, OpenProject, SaveProject, SaveAs, CloseProject,
    ImportStep, ImportSat, ExportStep, ExportSat,
    // View
    ToggleGrid, ToggleRuler, ToggleCoordSystem,
    RenderShaded, RenderWireframe, FitAll, ZoomIn, ZoomOut,
    // Simulation
    Validate, AddSetup, AddSweep, Solve, SolveAll, Abort,
    // Draw
    DrawBox, DrawCylinder, DrawSphere, DrawCone, DrawTorus,
    DrawRectangle, DrawEllipse, DrawCircle, DrawPolygon,
    DrawPolyline, DrawArc, DrawSpline, SetPlane, SetUnits,
    AssignMaterial,
    // Model
    BoolUnite, BoolSubtract, BoolIntersect, BoolSplit,
    GroupObjects, UngroupObjects, AssignColor, SetTransparency,
    // Results
    CreateReport, SolutionData, PlotFields,
    PlotEField, PlotHField, PlotSAR, Animate,
}
```

**Ribbon State**:
```rust
pub struct RibbonState {
    pub active_tab: usize,
    pub toggles: HashMap<String, bool>,  // Grid, Ruler, CoordSystem
    open_popup: Option<Id>,
}
```

#### File Operations
- **Save**: MessagePack format (.emsp extension)
- **Load**: Full project deserialization
- **Dialogs**: Async rfd dialogs (native + WASM compatible)
- **Polling**: `poll_file_dialogs()` checks completion each frame

### Key Capabilities
✅ Docking layout with egui_dock
✅ Ribbon-style toolbar (AEDT-inspired)
✅ Async file I/O with native OS dialogs
✅ Project state management
✅ Undo/redo ready (not yet implemented)

---

## 5. Integration Points & Capabilities Summary

### What's Available for UI Layout Integration

#### Rendering
1. **Simple Viewport** (SceneViewport)
   - Quick offscreen WGPU rendering
   - Basic status display
   - Minimal configuration

2. **Advanced 3D Viewport** (FieldSceneState)
   - Full orbit camera with 7 presets
   - 5 visualization modes
   - Real-time phase animation
   - Vector field arrow rendering
   - Slice plane visualization
   - 4 professional colormaps
   - Wireframe overlay
   - Opacity control
   - Colorbar legend

#### Data & State
1. **Domain Models**
   - Project, EmModel, GeometryObject, Material
   - SolveResult with field data
   - Fully serializable (JSON + MessagePack)

2. **Backend Abstraction**
   - Standalone in-memory projects
   - Cloud/remote solving support
   - File I/O (.emsp format)
   - Pluggable solver implementations

#### UI Framework
1. **Layout System**
   - egui_dock for docking panels
   - egui_table, egui_tiles for grids/layouts
   - egui_plot for plotting
   - Ribbon-style toolbar (AEDT-inspired)

2. **Input Handling**
   - Mouse drag/scroll interaction
   - Keyboard shortcuts (extensible)
   - Async file dialogs

### Data Flow

```
UI Layer (App)
    ↓
Domain Models (Project, EmModel, etc.)
    ↓
Infrastructure (Backend, File I/O)
    ↓
Solver (PlaceholderSolver → Real Solver)
    ↑
Results (SolveResult, field data)
    ↓
Render Layer (FieldSceneState, FieldPipeline)
    ↓
GPU (wgpu, egui-wgpu integration)
    ↓
Screen Output
```

---

## 6. Technical Constraints & Considerations

### Cross-Platform (WASM + Native)
- File dialogs: rfd (async, works both platforms)
- GPU: wgpu (universal, native offscreen + WASM canvas)
- Math: glam (no platform dependencies)
- Feature gates in Cargo.toml for platform-specific code

### GPU Rendering
- Offscreen framebuffer approach (prevents z-order issues with egui)
- Blit to egui render pass in paint callback
- wgsl shader language (WGPU standard)
- Support for complex field animations via vertex updates

### Performance Considerations
- Pre-generated meshes (sphere, cube, far-field, slice)
- Instance rendering for arrows (up to 4096)
- Colormap as GPU texture (256 or 512 LUT)
- Frame counting and status tracking

### UI Thread Constraints
- Synchronous UI rendering (egui paradigm)
- File dialogs spawn async tasks with channels
- GPU work in callback prepare() phase

---

## 7. Ready-to-Use Components

✅ **Complete 3D Visualization Stack**
- Orbit camera with presets
- 5 visualization modes (Surface, Arrows, Slice, FarField, Animation)
- Real-time phase animation
- Colormap system
- Wireframe rendering
- Vector field visualization

✅ **Project Management**
- File I/O (save/load)
- Project structure (model, materials, results)
- Undo/redo infrastructure (state-based)

✅ **UI Framework**
- Docking layout system
- Ribbon toolbar
- File dialogs
- Status/log display

✅ **Math & Physics**
- Camera matrices (view, projection, combined)
- Complex field phase computation
- Field range calculation
- Vector field synthesis

---

## 8. Next Steps for Integration

1. **Extend Domain Model** - Add simulation setup parameters
2. **Implement Actual Solver** - Replace PlaceholderSolver
3. **Field Data Pipeline** - Connect solver output to visualization
4. **UI Components** - Build model editor, parameter panels
5. **Scene Management** - Improve geometry object handling

All necessary rendering and visualization infrastructure is already in place and waiting for higher-level UI logic to drive it.
