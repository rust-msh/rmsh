# EMStudio Render Crate - Quick Reference Guide

## File Paths & Key Structs

### Core Data Structures

**Location**: `crates/render/src/mesh_data.rs`

```
FieldVertex (28 bytes)
├── position: [f32; 3]        (0-11)
├── normal: [f32; 3]          (12-23)
└── field_value: f32          (24-27)

FieldMesh
├── vertices: Vec<FieldVertex>
├── indices: Vec<u32>         (triangles: 3 per tri)
├── wire_indices: Vec<u32>    (edges: 2 per edge)
├── field_range: [f32; 2]
├── field_imag: Option<Vec<f32>>
└── vector_field: Option<Vec<[f32; 3]>>
```

### GPU Pipeline

**Location**: `crates/render/src/field_pipeline.rs`

```
FieldPipeline
├── scene_pipeline (filled triangles)
├── wire_pipeline (edges only)
├── arrow_pipeline (instanced arrows)
├── vertex_buf, index_buf, wire_index_buf (GPU buffers)
├── uniform_buf (FieldUniforms per frame)
├── colormap_texture (1D LUT, 256x1 RGBA8)
├── Offscreen FB: color_texture (Rgba8UnormSrgb) + depth_view
└── blit_pipeline (to egui)

FieldUniforms (80 bytes, per frame)
├── mvp: [f32; 16]
├── eye_pos: [f32; 3]
├── light_dir: [f32; 3]
├── field_min/max: f32
└── opacity: f32
```

### Scene State

**Location**: `crates/render/src/scene.rs`

```
FieldSceneState
├── camera: OrbitCamera
├── colormap: ColormapType (enum)
├── vis_mode: VisMode (enum)
├── animator: PhaseAnimator
├── sphere_mesh, cube_mesh, far_field_mesh
├── render_state: Option<egui_wgpu::RenderState>
└── [dirty flags: colormap_dirty, mode_dirty, slice_dirty]

VisMode enum
├── Surface (sphere with harmonic field)
├── Arrows (cube with vortex vectors)
├── Slice (plane through volume)
├── FarField (radiation pattern)
└── Animation (phase-swept complex field)
```

---

## Synthetic Data Generators

### All Located in: `crates/render/src/mesh_data.rs` (except noted)

| Function | Input | Output | Use Case |
|----------|-------|--------|----------|
| `FieldMesh::uv_sphere(n_lat, n_lon, r)` | Lat/lon subdivisions | Sphere mesh with sin(3φ)cos(2θ) field | Surface, Animation modes |
| `FieldMesh::cube(subdiv, half_size)` | Grid size | Cube with vortex vector field | Arrows mode |
| `generate_slice_mesh(axis, pos, extent, res, fn)` | Axis, position, resolution | Planar grid mesh | Slice mode |
| `generate_pattern_mesh(n_θ, n_φ, gain_fn)` | Theta/phi res | Radius-modulated sphere | FarField mode |
| `generate_arrow_base_mesh()` | — | Shaft + cone arrow | ArrowPipeline base |

**Locations**:
- slice.rs: `generate_slice_mesh()`, `synthetic_volume_field()`
- far_field.rs: `generate_pattern_mesh()`, `patch_gain()`
- mesh_data.rs: Everything else

---

## Data Flow (Simplified)

```
User Input
   ↓ (mouse, UI)
FieldSceneState::show_viewport()
   ├─→ Camera interaction
   ├─→ Check dirty flags
   ├─→ Build FieldUniforms
   └─→ Callback to GPU
        ↓
FieldSceneCallback::prepare()
   ├─→ Update GPU buffers
   └─→ Render to offscreen FB
        ↓
FieldSceneCallback::paint()
   ├─→ Blit offscreen → egui
   └─→ Display on panel
```

---

## File Format Specs

### .msh Format (Gmsh MSH 4.1)

**Use**: Store computational mesh (nodes, elements, boundaries)

**Key Sections**:
- `$Entities`: Geometric topology (0D points → 3D volumes)
- `$Nodes`: Node coordinates + tags
- `$Elements`: Element connectivity (triangle, tetrahedron, etc.)

**Element Types**:
- Type 1: Line (2 nodes)
- Type 2: Triangle (3 nodes)
- Type 4: Tetrahedron (4 nodes) ← Primary for 3D FEM

**Data**:
- f64 for coordinates (8 bytes each)
- u32/u64 for tags

**Size Estimate**: 10K nodes ≈ 1 MB

**Random Access**: Yes (entity block indexing)

### .emsfld Format (EMStudio Field Binary)

**Use**: Store frequency-domain FEM field solutions (E/H/J fields)

**Layout**:
```
Header (128 bytes)
├── Magic: b"EMSFLD\0\0"
├── num_nodes, num_frequencies, field_type
└── Offsets: freq_table, index, data

Frequency Table (f64 × num_freq)
Frequency Block Index (16 bytes × num_freq)
Field Blocks (frequency point data)
```

**Field Block Data** (per frequency):
- Complex f64 vector field: `48 bytes/node` (re_x, im_x, re_y, im_y, re_z, im_z)
- Complex f32 vector field: `24 bytes/node`

**Size Estimate**: 10K nodes, 301 frequencies ≈ 144 MB

**Random Access Pattern**:
```
1. Seek to (index_offset + freq_idx × 16)
2. Read FieldBlockInfo {offset, size}
3. Seek to offset
4. Read size_bytes
```

---

## How to Add Real Data Loaders

The render crate currently uses **synthetic generators only**. To load real `.msh`/`.emsfld`:

### Step 1: Create Data Loaders (new crate)
```rust
// emstudio-data-loaders
pub trait FieldDataSource {
    fn load_mesh(&self, freq_idx: usize) -> Result<FieldMesh>;
    fn frequencies(&self) -> Vec<f64>;
}

pub struct MshLoader { /* parse .msh */ }
pub struct FldLoader { /* memory-map .emsfld */ }
```

### Step 2: Integrate with FieldSceneState
```rust
pub struct FieldSceneState {
    pub data_source: Option<Arc<dyn FieldDataSource>>,
    // ... rest of fields
}

pub fn init_gpu(&mut self, render_state: ..., data_source: Option<Arc<...>>) {
    let mesh = if let Some(ds) = &self.data_source {
        ds.load_mesh(0)?
    } else {
        FieldMesh::uv_sphere(32, 64, 1.0)  // Fallback
    };
    // ...
}
```

### Step 3: Add Frequency Selector UI
```rust
pub fn show_controls(&mut self, ui: &mut egui::Ui) {
    // ... existing controls ...
    if let Some(ds) = &self.data_source {
        let freqs = ds.frequencies();
        // ComboBox to select frequency
        // On change: load_mesh(selected_freq)
    }
}
```

---

## GPU Resource Sizes (Typical 10K-vertex Scene)

| Resource | Size | Notes |
|----------|------|-------|
| Vertex buffer | 280 KB | 10K × 28 bytes |
| Index buffer (triangles) | 120 KB | 30K indices |
| Wire buffer (edges) | 120 KB | 30K indices |
| Colormap texture | 1 KB | 256×1 RGBA8 |
| Uniform buffer | 80 bytes | Per-frame constants |
| Offscreen FB (1920×1080) | ~33 MB | Color (Rgba8) + Depth (f32) |
| **Total VRAM** | ~33.5 MB | Dominated by framebuffer |

---

## Camera Controls

**Location**: `crates/render/src/camera.rs`

```
OrbitCamera
├── target: Vec3              (look-at point)
├── distance: f32             (radius from target)
├── azimuth: f32              (rotation around Y, radians)
├── elevation: f32            (vertical angle, radians)
└── fov_y: f32, near, far

Methods:
├── rotate(dx, dy)            (mouse drag → azimuth/elevation)
├── zoom(delta)               (scroll → distance, log scale)
├── pan(dx, dy)               (middle click → move target)
├── view_projection(aspect)   (→ Mat4 for GPU)
└── set_preset(preset)        (Front, Back, Left, Right, Top, Iso)
```

---

## Animation for Complex Fields

**Location**: `crates/render/src/animation.rs`

```
PhaseAnimator
├── phase_deg: f32            (current phase 0°-360°)
├── playing: bool
├── speed_deg_per_sec: f32

Apply formula:
  E(t) = Re(E) × cos(φ) - Im(E) × sin(φ)

Per-frame:
  1. tick(dt) → update phase_deg
  2. apply(real_values, imag_values) → time-domain values
  3. pipeline.update_vertices() → GPU update
```

---

## Colormaps

**Location**: `crates/render/src/colormap.rs`

```
ColormapType enum
├── Rainbow    (HSV 240°→0°)
├── Viridis    (perceptually uniform)
├── CoolWarm   (diverging: blue → gray → red)
└── Grayscale  (0.0→1.0 → black→white)

Method:
  colormap.generate_lut(256) → Vec<[u8; 4]> (RGBA)
  Then: texture = device.create_texture_init(lut_data)
```

---

## Integration Checklist

- [x] GPU pipeline setup (field_pipeline.rs)
- [x] Scene state management (scene.rs)
- [x] Interactive camera (camera.rs)
- [x] Synthetic data for demos (mesh_data.rs)
- [x] Phase animation (animation.rs)
- [x] Colormaps (colormap.rs)
- [x] Slice plane generation (slice.rs)
- [x] Far-field pattern (far_field.rs)
- [x] Arrow instancing (arrow_pipeline.rs)
- [ ] `.msh` file loader
- [ ] `.emsfld` file loader
- [ ] Quantity expression evaluator (dB, phase, etc.)
- [ ] Data caching/memory management

---

## Common Tasks

### Load a Custom Mesh
```rust
let mesh = FieldMesh {
    vertices: vec![...],  // Must be FieldVertex
    indices: vec![...],   // Triangle indices
    wire_indices: vec![...],  // Edge indices
    field_range: [min, max],
    field_imag: None,
    vector_field: None,
};
state.sphere_mesh = Some(mesh);
state.mode_dirty = true;
```

### Change Colormap
```rust
state.colormap = ColormapType::Viridis;
state.colormap_dirty = true;
```

### Play Animation
```rust
state.vis_mode = VisMode::Animation;
state.animator.playing = true;
state.animator.speed_deg_per_sec = 90.0;
```

### Toggle Wireframe
```rust
state.show_wireframe = !state.show_wireframe;
```

---

## Key Exports from lib.rs

```rust
pub use camera::{OrbitCamera, ViewPreset};
pub use colormap::ColormapType;
pub use mesh_data::{FieldMesh, FieldVertex};
pub use scene::{FieldSceneState, VisMode};
```

