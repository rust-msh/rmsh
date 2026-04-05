# EMStudio Render Crate - Comprehensive Exploration Report

**Date**: 2026-04-05  
**Project**: EMStudio  
**Focus**: emstudio-render crate structure, data flow, and mesh format specifications

---

## 1. CRATE FILE STRUCTURE

### Directory Layout
```
crates/render/src/
├── lib.rs                  # Main crate exports
├── scene.rs               # High-level visualization state and UI
├── field_pipeline.rs      # GPU rendering pipeline (wgpu) for field colormaps
├── arrow_pipeline.rs      # GPU pipeline for instanced arrow rendering
├── mesh_data.rs           # Core mesh data structures and synthetic generators
├── animation.rs           # Phase animation for complex field data
├── slice.rs               # Slice plane generation for volume visualization
├── far_field.rs           # Far-field radiation pattern mesh generation
├── camera.rs              # Orbit camera controller with view presets
├── colormap.rs            # Colormap lookup tables (Rainbow, Viridis, etc.)
├── field_shader.wgsl      # WGSL shader for field visualization and blitting
└── [No separate tests/examples in this crate]
```

**Total Files**: 11 source files  
**Lines of Code**: ~2,500 lines (mostly GPU pipeline setup and UI logic)

---

## 2. KEY STRUCTS AND DATA STRUCTURES

### 2.1 Mesh Data (`mesh_data.rs`)

#### `FieldVertex` - GPU Vertex Format
```rust
#[repr(C)]
pub struct FieldVertex {
    pub position: [f32; 3],      // 3D position (x, y, z)
    pub normal: [f32; 3],        // Surface normal for lighting
    pub field_value: f32,        // Scalar field value (colored by colormap)
}
// Total: 28 bytes per vertex
// Buffer layout: Float32x3 @ offset 0, Float32x3 @ offset 12, Float32 @ offset 24
```

#### `ArrowInstance` - GPU Instance Data (for Instanced Rendering)
```rust
#[repr(C)]
pub struct ArrowInstance {
    pub position: [f32; 3],       // Arrow base position
    pub direction: [f32; 3],      // Normalized direction vector
    pub magnitude: f32,           // Magnitude (for arrow length scaling)
    pub _pad: f32,                // Padding for alignment
}
// Total: 32 bytes per instance
// Step mode: Instance (one per arrow, not per vertex)
```

#### `FieldMesh` - Complete Mesh for Field Visualization
```rust
pub struct FieldMesh {
    pub vertices: Vec<FieldVertex>,           // All mesh vertices
    pub indices: Vec<u32>,                    // Triangle indices (filled mesh)
    pub wire_indices: Vec<u32>,               // Line indices (wireframe only)
    pub field_range: [f32; 2],                // [min, max] for colormap scaling
    pub field_imag: Option<Vec<f32>>,         // Imaginary part (for phase animation)
    pub vector_field: Option<Vec<[f32; 3]>>, // 3D vector field (for arrows)
}
```

### 2.2 Rendering Pipeline (`field_pipeline.rs`)

#### `FieldUniforms` - Per-Frame GPU Constants
```rust
#[repr(C)]
pub struct FieldUniforms {
    pub mvp: [f32; 16],           // Model-View-Projection matrix
    pub eye_pos: [f32; 3],        // Camera eye position (for lighting)
    pub _pad0: f32,
    pub light_dir: [f32; 3],      // Light direction (eye-relative)
    pub _pad1: f32,
    pub field_min: f32,           // Min value for colormap normalization
    pub field_max: f32,           // Max value for colormap normalization
    pub opacity: f32,             // Alpha blending opacity (0-1)
    pub _pad2: f32,
}
// Total: 80 bytes (must be 16-byte aligned)
```

#### `FieldPipeline` - GPU Rendering State Machine
**Primary Responsibility**: Manages all wgpu resources for rendering field-colored meshes

**Key Components**:
1. **Scene Pipelines**:
   - `scene_pipeline`: Renders filled triangles with field-value colormap
   - `wire_pipeline`: Renders wireframe edges on top of scene
   
2. **Buffers**:
   - `vertex_buf`: Stores all FieldVertex data (GPU VRAM)
   - `index_buf`: Triangle indices (3 per triangle)
   - `wire_index_buf`: Edge indices (2 per edge)
   - `uniform_buf`: Per-frame FieldUniforms (updated each frame)

3. **Colormaps**:
   - `colormap_texture`: 1D texture (256x1) with RGBA8 LUT
   - `colormap_sampler`: Linear filtering sampler
   - Updated dynamically when user switches colormaps

4. **Offscreen Framebuffer** (deferred rendering):
   - `color_texture`: Rgba8UnormSrgb (with depth)
   - `depth_view`: Depth32Float (for depth testing)
   - Resized to match viewport size (rounded up to 16-byte alignment)
   
5. **Blit Pipeline** (screen output):
   - `blit_pipeline`: Copies offscreen result to egui render pass
   - Uses full-screen triangle quad (no geometry needed)

### 2.3 Scene State (`scene.rs`)

#### `FieldSceneState` - High-Level Visualization State
```rust
pub struct FieldSceneState {
    pub camera: OrbitCamera,                  // Interactive 3D camera
    pub colormap: ColormapType,               // Selected colormap
    pub opacity: f32,                         // Transparency (0-1)
    pub show_wireframe: bool,                 // Toggle wireframe overlay
    pub field_range: [f32; 2],                // Min/max for colormap
    pub vis_mode: VisMode,                    // Current visualization mode
    
    // Animation
    pub animator: PhaseAnimator,              // Complex field phase animation
    last_frame_time: Option<Instant>,
    
    // Slice control
    pub slice_z: f32,                         // Z position for slice plane (-1 to +1)
    
    // Pre-generated meshes (for mode switching)
    sphere_mesh: Option<FieldMesh>,           // Spherical harmonic field demo
    cube_mesh: Option<FieldMesh>,             // Vector field (arrows) demo
    far_field_mesh: Option<FieldMesh>,        // Radiation pattern demo
    
    // Internal state
    render_state: Option<egui_wgpu::RenderState>,
    colormap_dirty: bool,
    mode_dirty: bool,
    slice_dirty: bool,
}
```

#### `VisMode` - Visualization Modes
```rust
pub enum VisMode {
    Surface,     // Colormap on UV sphere (synthetic spherical harmonic)
    Arrows,      // Vector arrows on cube surface (synthetic vortex field)
    Slice,       // Slice plane through volume (synthetic wave-like field)
    FarField,    // 3D radiation pattern (synthetic patch antenna gain)
    Animation,   // Phase-swept complex field (uses field_imag)
}
```

### 2.4 Camera (`camera.rs`)

#### `OrbitCamera` - Interactive 3D Camera
```rust
pub struct OrbitCamera {
    pub target: Vec3,                         // Look-at target
    pub distance: f32,                        // Distance from target
    pub azimuth: f32,                         // Rotation angle around Y (radians)
    pub elevation: f32,                       // Vertical angle (radians)
    pub fov_y: f32,                           // Field of view (Y axis)
    pub near: f32,                            // Near clipping plane
    pub far: f32,                             // Far clipping plane
}
```

**Methods**:
- `rotate(dx, dy)`: Update azimuth/elevation from mouse movement
- `zoom(delta)`: Adjust distance logarithmically
- `pan(dx, dy)`: Move target in screen-space (preserves distance)
- `view_matrix()`: Compute glam::Mat4 for camera transform
- `projection_matrix(aspect)`: Perspective matrix (RHS)
- `view_projection(aspect)`: Combined MVP
- `set_preset(preset)`: Jump to standard views (Front, Back, Left, Right, Top, Iso)

### 2.5 Animation (`animation.rs`)

#### `PhaseAnimator` - Complex Field Time-Domain Synthesis
```rust
pub struct PhaseAnimator {
    pub phase_deg: f32,              // Current phase [0°, 360°)
    pub playing: bool,               // Is animation active?
    pub speed_deg_per_sec: f32,      // Animation speed
}

// Computes: E(t) = Re(E) * cos(φ) - Im(E) * sin(φ)
pub fn apply(&self, field_real: &[f32], field_imag: &[f32]) -> Vec<f32>
```

---

## 3. SYNTHETIC TEST DATA GENERATION

### How Test Data is Created

The crate does **not** load real EMStudio solver results (`.msh`, `.emsfld`). Instead, it has **built-in synthetic data generators** for demonstration:

#### 3.1 UV Sphere with Spherical Harmonics (`mesh_data.rs`)
```rust
pub fn FieldMesh::uv_sphere(n_lat: u32, n_lon: u32, radius: f32) -> Self
```

**Generated Data**:
- **Vertices**: `(n_lat+1) × (n_lon+1)` vertices on sphere
- **Field Value**: `sin(3φ) × cos(2θ)` (spherical harmonic pattern)
- **Imaginary Part**: `cos(2φ) × sin(3θ)` (for phase animation)
- **Use Case**: VisMode::Surface and VisMode::Animation modes

**Default Creation** (in `scene.rs:init_gpu`):
```rust
let sphere_mesh = FieldMesh::uv_sphere(32, 64, 1.0);  // Implied default
```

#### 3.2 Subdivided Cube with Vortex Vector Field (`mesh_data.rs`)
```rust
pub fn FieldMesh::cube(subdivisions: u32, half_size: f32) -> Self
```

**Generated Data**:
- **Mesh**: 6 faces of cube, each subdivided into quads
- **Vector Field**: `synthetic_vector_field(pos) -> (f32, f32, f32)`
  ```rust
  let vx = -z * 0.8 + y * 0.3;
  let vy = (x² + z²).sqrt().sin() * 0.5;
  let vz = x * 0.8 - y * 0.3;  // Rotating vortex pattern
  ```
- **Field Value**: Magnitude of vector field
- **Use Case**: VisMode::Arrows mode

**Default Creation** (in `scene.rs:init_gpu`):
```rust
let cube_mesh = FieldMesh::cube(10, 1.0);  // 11×11 per face
```

#### 3.3 Slice Plane with Synthetic Volume Field (`slice.rs`)
```rust
pub fn generate_slice_mesh(
    axis: SliceAxis,
    value: f32,
    extent: f32,
    resolution: u32,
    field_fn: &dyn Fn(f32, f32, f32) -> f32,
) -> FieldMesh
```

**Generated Data**:
- **Grid**: `resolution × resolution` quad grid on XY/XZ/YZ plane
- **Field Function**: `synthetic_volume_field(x, y, z) -> f32`
  ```rust
  let r = sqrt(x² + y² + z²);
  sin(4r) / r  // Radial wave pattern
  ```
- **Use Case**: VisMode::Slice mode

#### 3.4 Far-Field Radiation Pattern (`far_field.rs`)
```rust
pub fn generate_pattern_mesh(
    n_theta: u32,
    n_phi: u32,
    gain_fn: &dyn Fn(f32, f32) -> f32,
) -> FieldMesh
```

**Generated Data**:
- **Parametric Surface**: Sphere whose radius is modulated by gain
- **Default Gain Function**: `patch_gain(theta, phi)` (patch antenna pattern)
- **Field Value**: Gain in dBi
- **Use Case**: VisMode::FarField mode

**Default Creation** (in `scene.rs:init_gpu`):
```rust
let far_field_mesh = far_field::generate_pattern_mesh(60, 120, &far_field::patch_gain);
```

#### 3.5 Arrow Mesh (`mesh_data.rs`)
```rust
pub fn generate_arrow_base_mesh() -> (Vec<[f32; 3]>, Vec<u32>)
```

**Generated Data**:
- **Shaft**: Cylinder with `segments=8` rings (from vertex buffer)
- **Head**: Cone tip pointing in +Y direction
- **Dimensions**:
  - Shaft radius: 0.02
  - Shaft length: 0.7
  - Head radius: 0.06
  - Head length: 0.3

**Usage**: Created once and instanced in ArrowPipeline

---

## 4. DATA FLOW: INPUT TO GPU

### Pipeline: Data → GPU → Screen

```
┌─────────────────────────────────────────────────────────────┐
│ 1. Data Creation Layer                                      │
├─────────────────────────────────────────────────────────────┤
│ FieldMesh::uv_sphere()          → Vec<FieldVertex>          │
│ FieldMesh::cube()               → Vec<FieldVertex>          │
│ generate_slice_mesh()           → Vec<FieldVertex>          │
│ generate_pattern_mesh()         → Vec<FieldVertex>          │
│ cube.generate_arrows(every_n)   → Vec<ArrowInstance>        │
└──────────────────┬──────────────────────────────────────────┘
                   │
┌──────────────────▼──────────────────────────────────────────┐
│ 2. GPU Upload Layer (FieldPipeline::new)                   │
├─────────────────────────────────────────────────────────────┤
│ device.create_buffer_init(&vertices)    → vertex_buf        │
│ device.create_buffer_init(&indices)     → index_buf         │
│ device.create_buffer_init(&wire_indices) → wire_index_buf   │
│ colormap.generate_lut(256)              → colormap_texture  │
└──────────────────┬──────────────────────────────────────────┘
                   │
┌──────────────────▼──────────────────────────────────────────┐
│ 3. Per-Frame Rendering (FieldSceneState::show_viewport)   │
├─────────────────────────────────────────────────────────────┤
│ camera.view_projection(aspect)          → FieldUniforms.mvp │
│ queue.write_buffer(&uniform_buf, ...)   → GPU constants     │
│ phase_animator.apply(real, imag)        → new field_values  │
│ pipeline.update_vertices(queue, verts)  → GPU vertex update │
└──────────────────┬──────────────────────────────────────────┘
                   │
┌──────────────────▼──────────────────────────────────────────┐
│ 4. GPU Rendering (FieldPipeline::render_scene)            │
├─────────────────────────────────────────────────────────────┤
│ render_pass.set_pipeline(&scene_pipeline)                   │
│ render_pass.set_vertex_buffer(0, vertex_buf)                │
│ render_pass.set_index_buffer(index_buf)                     │
│ render_pass.draw_indexed(0..num_indices)                    │
│   → Executes field_shader.wgsl::vs_main + fs_main          │
│                                                              │
│ [Optional] Wireframe overlay:                               │
│ render_pass.set_pipeline(&wire_pipeline)                    │
│ render_pass.draw_indexed(0..num_wire_indices)               │
│                                                              │
│ [Optional] Arrow instancing:                                │
│ arrow_pipeline.draw(&mut rpass, &scene_bind_group)          │
└──────────────────┬──────────────────────────────────────────┘
                   │
┌──────────────────▼──────────────────────────────────────────┐
│ 5. Framebuffer Blit (FieldPipeline::blit)                 │
├─────────────────────────────────────────────────────────────┤
│ Offscreen color_texture (Rgba8UnormSrgb)                    │
│   ↓ Sampled by blit_pipeline                                │
│ egui render pass (target_format)                            │
│   ↓ Displayed on screen                                     │
└─────────────────────────────────────────────────────────────┘
```

### Key Data Sizes

For a typical demo with 10K vertices:

| Item | Size | Notes |
|------|------|-------|
| `FieldVertex` (1 vertex) | 28 bytes | 3 float32 + 3 float32 + 1 float32 |
| Vertex buffer (10K) | 280 KB | Fits in GPU cache |
| Index buffer (30K triangles) | 120 KB | 3 indices per triangle |
| Wireframe buffer (30K edges) | 120 KB | 2 indices per edge |
| Colormap texture (256x1) | 1 KB | RGBA8 |
| Uniform buffer | 80 bytes | Updated per frame |

---

## 5. FieldSceneState AND FieldPipeline EXPECTATIONS

### What FieldPipeline Expects as Input

#### Constructor Input
```rust
pub fn new(
    device: &wgpu::Device,              // GPU device for resource creation
    queue: &wgpu::Queue,                // GPU command queue
    target_format: wgpu::TextureFormat, // Output texture format (e.g., Rgba8UnormSrgb)
    mesh: &FieldMesh,                   // Mesh with vertices, indices, field_range
    colormap: ColormapType,             // Color mapping function
) -> Self
```

#### Required FieldMesh Structure
- **vertices**: Vec<FieldVertex> with:
  - Valid 3D positions
  - Surface normals (for lighting)
  - Scalar field values (will be colormap-indexed)
- **indices**: Vec<u32> (triangle list, 3 per triangle)
- **wire_indices**: Vec<u32> (line list, 2 per line)
- **field_range**: [f32; 2] = [min_value, max_value] for colormap normalization
- **field_imag**: Optional Vec<f32> (for phase animation mode)
- **vector_field**: Optional Vec<[f32; 3]> (for arrow rendering)

#### Per-Frame Inputs (update_uniforms, update_colormap, etc.)
```rust
pub fn update_uniforms(&self, queue: &wgpu::Queue, uniforms: &FieldUniforms)
pub fn update_colormap(&mut self, device: &wgpu::Device, queue: &wgpu::Queue, 
                       colormap: ColormapType)
pub fn resize_if_needed(&mut self, device: &wgpu::Device, size: [u32; 2])
pub fn update_vertices(&self, queue: &wgpu::Queue, vertices: &[FieldVertex])
pub fn swap_mesh(&mut self, device: &wgpu::Device, mesh: &FieldMesh)
```

### What FieldSceneState Expects as Input

#### Constructor/Initialization
```rust
pub fn init_gpu(&mut self, render_state: &egui_wgpu::RenderState, mesh: &FieldMesh)
```

**Expects**:
- Valid wgpu RenderState (from egui integration)
- Initial FieldMesh to display
- GPU device/queue available

#### UI Control Loop (show_viewport, show_controls)
```rust
pub fn show_viewport(&mut self, ui: &mut egui::Ui)  // Renders 3D view
pub fn show_controls(&mut self, ui: &mut egui::Ui)  // Shows control panel
pub fn show_colorbar(&self, ui: &mut egui::Ui)      // Legend
```

**Expects**:
- egui::Ui context for drawing
- User input (mouse, keyboard) handled by egui
- Sufficient viewport space (allocates full available width/height)

---

## 6. MESH FORMAT SPECIFICATIONS

### 6.1 .msh File Format (Gmsh MSH 4.1)

**Source**: docs/em-result-file-formats.md § 3.2

**Purpose**: Store computational mesh (nodes, elements, boundaries)

**File Structure**:
```
$MeshFormat
  4.1 0 8            # Version 4.1, binary format, 8-byte reals
$EndMeshFormat

$Entities
  numPoints numCurves numSurfaces numVolumes
  // ... geometric entity definitions
$EndEntities

$Nodes
  numEntityBlocks numNodes minNodeTag maxNodeTag
  entityDim entityTag numNodesInEntity nodeTag x y z  [...]
  // ... per-entity node blocks
$EndNodes

$Elements
  numEntityBlocks numElements minElementTag maxElementTag
  entityDim entityTag elementType numElementsInEntity elementTag nodeTag₁ ... nodeTagₙ [...]
  // ... per-entity element blocks
$EndElements

[Optional] $NodeData, $ElementData (quality metrics)
```

**Element Types**:
- **2** = Line (2 nodes)
- **1** = Triangle (3 nodes)
- **2** = Triangle (3 nodes)
- **4** = Tetrahedron (4 nodes) ← Primary for 3D FEM
- **9** = Pyramid (5 nodes)
- **10** = Prism/Wedge (6 nodes)

**Data Types**:
- f64 for coordinates (8 bytes each)
- u32/u64 for tags/indices (depends on version)

**Typical Sizes**:
- Small mesh: 10K nodes → ~1 MB
- Medium mesh: 100K nodes → ~10 MB
- Large mesh: 1M+ nodes → ~100+ MB

**Random Access Capability**: Yes (via entity block indexing)

---

### 6.2 .emsfld File Format (EMStudio Field Binary)

**Source**: docs/em-result-file-formats.md § 3.3

**Purpose**: Store frequency-domain FEM field solutions (E/H/J fields at multiple frequencies)

**File Structure**:
```
┌─────────────────────────────────────────┐
│ Header (128 bytes, fixed)               │
├─────────────────────────────────────────┤
│ Frequency Table                         │
│ (num_frequencies × f64)                 │
├─────────────────────────────────────────┤
│ Field Block Index                       │
│ (num_frequencies × FieldBlockInfo)      │
├─────────────────────────────────────────┤
│ Field Block 0 (Frequency point 0)       │
├─────────────────────────────────────────┤
│ Field Block 1 (Frequency point 1)       │
├─────────────────────────────────────────┤
│ ...                                     │
└─────────────────────────────────────────┘
```

#### Header Structure (128 bytes)
```rust
#[repr(C, packed)]
pub struct EmsFldHeader {
    pub magic: [u8; 8],           // b"EMSFLD\0\0"
    pub version: u32,             // 1 (current)
    pub byte_order: u32,          // 0x01020304 (little-endian marker)
    pub field_type: u32,          // 0=E, 1=H, 2=J, 3=Combined
    pub data_type: u32,           // 0=complex f64, 1=complex f32
    pub num_nodes: u64,           // Number of field sample points (≈ mesh nodes)
    pub num_components: u32,      // 3 (vector), 1 (scalar)
    pub num_frequencies: u32,     // Frequency points
    pub frequency_unit: u32,      // 0=Hz, 1=kHz, 2=MHz, 3=GHz
    pub freq_table_offset: u64,   // Offset to frequency table
    pub index_offset: u64,        // Offset to field block index
    pub data_offset: u64,         // Offset to first field block
    pub mesh_file: [u8; 32],      // Associated .msh filename
    pub _reserved: [u8; 12],      // Future use
}
```

#### Field Block Layout (Single Frequency)
**For complex f64 vector field (E_x, E_y, E_z)**:
```
Per node:
  re_x (f64) | im_x (f64) | re_y (f64) | im_y (f64) | re_z (f64) | im_z (f64)
  = 48 bytes per node

Total block size: num_nodes × 48 bytes

Example (10K nodes, 1 frequency):
  10,000 × 48 = 480 KB

Example (10K nodes, 301 frequencies):
  10,000 × 48 × 301 = 144.5 MB
```

#### Random Access Pattern
```rust
1. Read header (first 128 bytes)
2. Seek to index_offset + freq_idx × 16 bytes
3. Read FieldBlockInfo {offset, size_bytes}
4. Seek to offset
5. Read size_bytes of field data
```

---

## 7. SCENE COMPOSITION

### How Scenes Are Built

#### Flow Diagram
```
FieldSceneState::init_gpu()
  ├── Creates 4 pre-rendered meshes:
  │   ├── sphere_mesh (FieldMesh::uv_sphere)
  │   ├── cube_mesh (FieldMesh::cube)
  │   ├── far_field_mesh (far_field::generate_pattern_mesh)
  │   └── [slice_mesh generated on-demand in VisMode::Slice]
  │
  ├── Creates FieldPipeline (GPU pipelines, buffers, colormaps)
  │
  ├── Creates ArrowPipeline (for VisMode::Arrows)
  │
  └── Stores in egui_wgpu::CallbackResources (thread-safe)

FieldSceneState::show_viewport()
  ├── Handles camera interaction (rotate, pan, zoom)
  ├── Detects mode/colormap/slice changes
  ├── Updates animation state (if playing)
  ├── Builds FieldUniforms from camera
  │
  ├── Creates FieldSceneCallback {uniforms, show_wireframe, show_arrows, ...}
  │
  └── Submits callback to egui_wgpu
        ├── prepare() → Upload uniforms, render to offscreen FBO
        └── paint()  → Blit offscreen texture to egui render pass
```

#### Mode Switching Logic (`apply_mode_switch()`)
```rust
match self.vis_mode {
    VisMode::Surface | VisMode::Animation => Use sphere_mesh
    VisMode::Arrows  => Use cube_mesh + generate_arrows(every_n=3)
    VisMode::Slice   => Generate new slice_mesh on-the-fly
    VisMode::FarField => Use far_field_mesh
}

Then:
  pipeline.swap_mesh(&device, &mesh)  // Replace GPU buffers
```

#### Animation Integration (PhaseAnimator)
```rust
if vis_mode == VisMode::Animation:
    animator.tick(dt)  // Increment phase_deg
    apply() → [field_real, field_imag] → [time_domain_values]
    pipeline.update_vertices(queue, modified_vertices)
```

---

## 8. FILE LOADING CODE (Current State)

### Status: NOT YET IMPLEMENTED

The crate currently **does not load** `.msh` or `.emsfld` files. It only has synthetic generators.

### Integration Points (Where Loading Would Hook)

1. **Data Acquisition Layer** (hypothetical):
   ```rust
   pub trait FieldDataSource {
       fn load_mesh(&self, freq_idx: usize) -> Result<FieldMesh, Error>;
       fn load_frequencies(&self) -> Result<Vec<f64>, Error>;
   }
   ```

2. **In FieldSceneState::init_gpu**:
   ```rust
   let mesh = if let Some(data_source) = &self.data_source {
       data_source.load_mesh(0)?
   } else {
       FieldMesh::uv_sphere(32, 64, 1.0)  // Fallback to synthetic
   };
   ```

3. **File Loaders Needed**:
   - `MshLoader`: Parse Gmsh MSH 4.1 → Convert to FieldVertex/indices
   - `FldLoader`: Memory-map or read .emsfld, seek to frequency point
   - `QuantityMapper`: Transform field data (e.g., dB scaling, phase extraction)

---

## 9. DOCUMENTATION REFERENCES

### Key Documentation Files

| File | Focus |
|------|-------|
| `docs/em-result-file-formats.md` | Complete `.msh` and `.emsfld` specs with Rust struct definitions |
| `docs/em-result-visualization-design.md` | High-level visualization architecture (Reports, Field Overlays, Far Field) |

### Key Design Principles from Docs

1. **Separation of Concerns**:
   - Data loading (IO, parsing)
   - Visualization mapping (field → color, arrows, etc.)
   - GPU rendering (wgpu pipelines)

2. **Per-Frequency Random Access**: Both `.msh` and `.emsfld` support selective loading (no need to load entire file)

3. **Complex Field Handling**: Store real + imaginary parts separately for phase animation and S-parameter plots

4. **Extensibility**: New visualization types (iso-surfaces, slice planes, etc.) added without changing core pipeline

---

## 10. SUMMARY TABLE: KEY COMPONENTS

| Component | File | Purpose | Input | Output |
|-----------|------|---------|-------|--------|
| **FieldVertex** | mesh_data.rs | GPU vertex format | (position, normal, field_value) | 28-byte packed struct |
| **FieldMesh** | mesh_data.rs | Complete mesh + metadata | vertices, indices, field_range | GPU-ready format |
| **FieldPipeline** | field_pipeline.rs | GPU rendering state | FieldMesh, colormap | Offscreen framebuffer |
| **FieldUniforms** | field_pipeline.rs | Per-frame GPU constants | Camera MVP, field min/max | 80-byte buffer |
| **FieldSceneState** | scene.rs | High-level UI state | User input, mode selection | Rendered 3D view + controls |
| **OrbitCamera** | camera.rs | Interactive camera | Mouse/scroll input | View/projection matrices |
| **PhaseAnimator** | animation.rs | Complex field time-domain | Phase angle, frequency data | Animated field values |
| **ArrowPipeline** | arrow_pipeline.rs | Instanced arrow rendering | Arrow positions/directions | Vector field overlay |
| **Colormaps** | colormap.rs | Color mapping LUTs | Field value [0,1] | RGBA8 color |

---

## 11. DATA FLOW DIAGRAM (Text)

```
User Action (UI)
  ↓
FieldSceneState::show_viewport()
  ├── Mouse/Keyboard → Camera update
  ├── Colormap combo → colormap_dirty = true
  ├── Mode selection → mode_dirty = true
  │
  ├── [if colormap_dirty] pipeline.update_colormap()
  ├── [if mode_dirty] apply_mode_switch() → pipeline.swap_mesh()
  ├── [if animation] animator.tick() → apply_phase_animation()
  │
  ├── Compute FieldUniforms from camera
  ├── pipeline.update_uniforms()
  │
  ├── Create FieldSceneCallback
  └── ui.painter().add(egui_wgpu::Callback)
        ↓
egui_wgpu CallbackTrait::prepare()
  ├── pipeline.resize_if_needed()
  ├── pipeline.render_scene(encoder, show_wireframe, arrow_pipeline)
  │   └── Render to offscreen framebuffer
  └── Return empty command buffers (already recorded in encoder)
        ↓
egui_wgpu CallbackTrait::paint()
  ├── pipeline.blit(render_pass)
  │   └── Full-screen triangle sampling color_texture
  └── Output appears on egui panel
```

---

## 12. INTEGRATION WITH OTHER CRATES

### Dependencies (from Cargo.toml)
```toml
[dependencies]
bytemuck = "..."           # Pod/Zeroable derives for GPU data
egui = "..."               # UI framework
egui-wgpu = "..."          # egui + wgpu integration
emstudio-domain = "..."    # Shared types (enums, structs)
glam = "..."               # Math (Vec3, Mat4)
wgpu = "..."               # GPU API abstraction
rcad-kernel = "..."        # CAD geometry (not used in render)
rcad-render = "..."        # CAD mesh rendering (not used in render)
```

### Crates That Import emstudio-render
- **emstudio-app**: Uses FieldSceneState for main visualization
- **emstudio-components**: May provide UI panels for controls

---

## CONCLUSION

The emstudio-render crate is a **complete, self-contained GPU rendering system** with:

✅ **Synthetic data generation** for demos and testing  
✅ **GPU pipelines** for field visualization (colormap, wireframe, arrows, slices, far-field)  
✅ **Interactive camera** with orbit/pan/zoom  
✅ **Animation support** for phase-swept complex fields  
✅ **Per-frame update** capability for interactive control  

❌ **File I/O not yet implemented** (`.msh` and `.emsfld` loaders needed)  
❌ **No loading layer** for real simulation results  

The architecture is **modular and extensible**, with clear separation between:
- **Data structures** (FieldMesh, FieldVertex, FieldUniforms)
- **GPU pipelines** (scene, wireframe, arrows, blit)
- **UI state** (FieldSceneState, camera, animation)

Ready for integration with result file loaders in a separate crate.
