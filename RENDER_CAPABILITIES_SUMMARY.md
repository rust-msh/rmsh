# EmStudio Render & Domain Layer - Capabilities at a Glance

## 🎨 Rendering Capabilities

### FieldSceneState - Complete 3D Visualization System

```
┌─────────────────────────────────────────────────────────────────────┐
│                    FieldSceneState                                   │
│  (Advanced 3D field visualization with GPU acceleration)             │
├─────────────────────────────────────────────────────────────────────┤
│                                                                       │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │  Visualization Modes (VisMode)                              │   │
│  ├──────────────────────────────────────────────────────────────┤   │
│  │  1. Surface    → Colormap on UV sphere (complex fields)     │   │
│  │  2. Arrows     → Vector field arrows on cube surface        │   │
│  │  3. Slice      → 2D slice plane through 3D volume           │   │
│  │  4. FarField   → 3D radiation pattern visualization         │   │
│  │  5. Animation  → Phase-swept animation of complex field     │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                                                                       │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │  Orbit Camera (7 presets)                                   │   │
│  ├──────────────────────────────────────────────────────────────┤   │
│  │  • Front, Back, Left, Right (orthographic views)            │   │
│  │  • Top, Bottom (plan/section views)                         │   │
│  │  • Iso (isometric 34.3°/23.0° angles)                       │   │
│  │  • Smooth zoom, pan, rotate (mouse interaction)             │   │
│  │  • Configurable FOV, near/far planes                        │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                                                                       │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │  Colormaps (4 professional schemes)                         │   │
│  ├──────────────────────────────────────────────────────────────┤   │
│  │  • Rainbow    (HSV sweep: blue→red)                         │   │
│  │  • Viridis    (Perceptually uniform: purple→yellow)         │   │
│  │  • Cool-Warm  (Diverging: blue→gray→red)                    │   │
│  │  • Grayscale  (Linear intensity)                            │   │
│  │  (GPU-accelerated lookup table textures)                    │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                                                                       │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │  Rendering Features                                         │   │
│  ├──────────────────────────────────────────────────────────────┤   │
│  │  ✓ Real-time phase animation (0-360°, configurable speed)  │   │
│  │  ✓ Wireframe overlay (toggleable)                           │   │
│  │  ✓ Opacity control (0.0 - 1.0)                             │   │
│  │  ✓ Dynamic mesh switching (sphere → cube → far-field)       │   │
│  │  ✓ Real-time vertex updates (complex field time-domain)    │   │
│  │  ✓ Offscreen framebuffer with depth (no Z-order issues)    │   │
│  │  ✓ Colorbar legend with range annotations                   │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                                                                       │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │  Pre-generated Meshes                                       │   │
│  ├──────────────────────────────────────────────────────────────┤   │
│  │  • UV Sphere   (configurable latitude/longitude resolution) │   │
│  │  • Cube        (subdivided faces with vector field)         │   │
│  │  • Far-Field   (60×120 patch antenna pattern)              │   │
│  │  • Slice Plane (X/Y/Z axis, configurable position)         │   │
│  │  (With synthetic field data for demo)                       │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                                                                       │
└─────────────────────────────────────────────────────────────────────┘
```

### GPU-Accelerated Rendering Pipeline

```
┌─────────────────────────────────────────────────────────────────┐
│  FieldPipeline (wgpu-based)                                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                   │
│  Vertex Buffers          Uniform Buffer        Colormap Texture  │
│  ├─ vertices             ├─ MVP matrix         ├─ RGBA8 LUT      │
│  ├─ indices              ├─ camera position    ├─ 256 colors     │
│  └─ wireframe indices    ├─ light direction    └─ samplers       │
│                          └─ field range                          │
│                                                                   │
│  Offscreen Render Pass (with depth buffer)                       │
│  ├─ Scene Pipeline (solid fill + colormap)                      │
│  └─ Wire Pipeline (overlay wireframe)                           │
│                                                                   │
│  Blit to egui render pass (no depth)                            │
│  └─ Integrates seamlessly with egui UI                          │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘
```

## 📊 Domain Models

```
┌─────────────────────────────────────────────────────────────────┐
│  Project (Top-level container)                                  │
├─────────────────────────────────────────────────────────────────┤
│  ├─ id: String                                                   │
│  ├─ title: String                                               │
│  ├─ status: SimulationStatus (Idle|Solving|Finished|Failed)     │
│  │                                                               │
│  ├─ model: EmModel                                              │
│  │  ├─ name: String                                             │
│  │  ├─ objects: Vec<GeometryObject>                             │
│  │  │  ├─ id: u64                                               │
│  │  │  ├─ name: String                                          │
│  │  │  └─ mesh_hint: String (auto|fine|coarse)                 │
│  │  │                                                            │
│  │  └─ materials: Vec<Material>                                 │
│  │     ├─ name: String                                          │
│  │     ├─ relative_permittivity: f32 (ε_r)                      │
│  │     └─ conductivity: f32 (σ)                                 │
│  │                                                               │
│  └─ last_result: Option<SolveResult>                            │
│     ├─ field_preview: String                                    │
│     └─ converged: bool                                          │
│                                                                   │
│  Serialization: JSON + MessagePack (.emsp)                       │
└─────────────────────────────────────────────────────────────────┘
```

## 🔄 Backend & File I/O

```
┌──────────────────────────────┐
│  Backend Trait               │
├──────────────────────────────┤
│  save_project()              │
│  load_project()              │
│  solve()                     │
│  mode()                      │
└──────────────────────────────┘
         ▲
         │ implements
         │
    ┌────┴─────────────┬──────────────┐
    │                  │              │
┌───────────┐   ┌──────────────┐   ┌────────────┐
│Standalone │   │    Cloud     │   │Custom Impl │
│Backend    │   │   Backend    │   │            │
├───────────┤   ├──────────────┤   └────────────┘
│HashMap    │   │HTTP Endpoint │
│in-memory  │   │(placeholder) │
└───────────┘   └──────────────┘

File I/O:
  • save_project_to_file(project, path) → .emsp
  • load_project_from_file(path) → Project
  • Format: MessagePack (binary, compact)
  • Cross-platform (native + WASM)
```

## 🎯 Camera & Math

### OrbitCamera Features
```
Eye Position Calculation:
  x = distance × cos(elevation) × sin(azimuth)
  y = distance × sin(elevation)
  z = distance × cos(elevation) × cos(azimuth)

Matrices (Right-handed, Y-up):
  • View Matrix       (via Mat4::look_at_rh)
  • Projection Matrix (via Mat4::perspective_rh)
  • Combined MVP      (projection × view)

Interactions:
  • Rotate: dx, dy → azimuth, elevation changes
  • Zoom:   delta → exponential distance scaling
  • Pan:    dx, dy → target translation in screen space
```

### Complex Field Animation
```
Phase Animator:
  phase_deg:       0° - 360° (configurable)
  speed_deg/sec:   10 - 720° (configurable)
  
Time-domain field computation:
  E(t) = Re(E) × cos(φ) - Im(E) × sin(φ)
  
Field envelope range (conservative):
  min = -√(Re² + Im²)
  max = +√(Re² + Im²)
```

## 🛠️ Integration Points

### From App Layer
```rust
// Initialize visualization
let mut scene = FieldSceneState::new();
scene.init_gpu(&render_state, &sphere_mesh);

// Each frame
scene.show_viewport(&mut ui);      // Render 3D view
scene.show_controls(&mut ui);      // Control panel
scene.show_colorbar(&ui);          // Legend

// Change visualization
scene.vis_mode = VisMode::Arrows;  // Set mode
scene.colormap = ColormapType::Viridis; // Change colormap
scene.opacity = 0.8;               // Change opacity
scene.camera.set_preset(ViewPreset::Iso); // Preset view
```

### From Domain Layer
```rust
// Domain models ready for binding
let project = Project::default();
println!("Model has {} objects", project.model.objects.len());
println!("Model has {} materials", project.model.materials.len());

// Serialize/deserialize
let json = serde_json::to_string(&project)?;
let saved = serde_json::from_str(&json)?;
```

### From Infra Layer
```rust
// File operations
save_project_to_file(&project, path)?;
let loaded = load_project_from_file(path)?;

// Solving
let result = backend.solve(&project)?;
println!("Converged: {}", result.converged);
```

## 📈 Mesh Data Types

```
FieldVertex (GPU format):
  position: [f32; 3]
  normal:   [f32; 3]
  field_value: f32

ArrowInstance (instanced rendering):
  position:  [f32; 3]
  direction: [f32; 3]
  magnitude: f32

FieldMesh (container):
  vertices: Vec<FieldVertex>
  indices: Vec<u32>
  wire_indices: Vec<u32>
  field_range: [f32; 2]
  field_imag: Option<Vec<f32>>        (complex field imaginary part)
  vector_field: Option<Vec<[f32; 3]>> (3D vector field for arrows)
```

## 🔧 Configuration & Debugging

```
WgpuRenderConfig:
  use_webgpu: bool        (auto-detect based on target)
  msaa_samples: u32       (1, 2, 4, 8, 16)

RuntimeStatus:
  PendingInit             (GPU resources not yet created)
  Ready                   (Rendering active)
  Unsupported(msg)        (Platform doesn't support wgpu)
  Failed(msg)             (GPU initialization failed)

Frame Counter:
  Tracks successful render frames
  Used for performance monitoring
```

## 💡 Ready-to-Use Features

| Feature | Status | Details |
|---------|--------|---------|
| **Orbit Camera** | ✅ Complete | 7 presets, full 3D control |
| **5 Vis Modes** | ✅ Complete | Surface, Arrows, Slice, FarField, Animation |
| **4 Colormaps** | ✅ Complete | Rainbow, Viridis, Cool-Warm, Grayscale |
| **Phase Animation** | ✅ Complete | Real-time playback, configurable speed |
| **Wireframe Rendering** | ✅ Complete | Toggle overlay |
| **Vector Field Arrows** | ✅ Complete | Up to 4096 instances, dynamic scaling |
| **Colorbar Legend** | ✅ Complete | With range annotations |
| **File I/O** | ✅ Complete | MessagePack format, async dialogs |
| **Project Management** | ✅ Complete | Standalone + Cloud backends |
| **Docking UI Layout** | ✅ Complete | egui_dock integration |
| **Ribbon Toolbar** | ✅ Complete | AEDT-inspired, 40+ actions |
| **Real-time Vertex Updates** | ✅ Complete | For animation support |

