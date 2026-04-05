# EmStudio Render & Domain Layer Exploration - Index

This folder contains comprehensive documentation of the EmStudio project's render, domain, and infrastructure layers. All exploration is read-only and focused on understanding the existing architecture.

## Generated Documents

### 1. **QUICK_REFERENCE.md** ⭐ START HERE
- **Best for**: Quick overview and navigation
- **Contains**: File structure, component summary, key capabilities, integration patterns
- **Length**: ~2 pages
- **Key sections**: 
  - File structure map
  - Core components at a glance
  - Common integration patterns
  - Quick checklist for getting started

### 2. **ARCHITECTURE_EXPLORATION.md** 🔬 DETAILED ANALYSIS
- **Best for**: Deep understanding of each layer
- **Contains**: Complete breakdown of all 5 layers (Domain, Infra, Render, App, Components)
- **Length**: ~8 pages
- **Key sections**:
  - Tech stack dependencies
  - Domain models with field descriptions
  - Backend abstraction and file I/O
  - Render layer (Scene, Camera, Mesh, Pipeline)
  - GPU rendering architecture
  - Integration points and data flow
  - Technical constraints and cross-platform notes

### 3. **RENDER_CAPABILITIES_SUMMARY.md** 🎨 VISUAL REFERENCE
- **Best for**: Visualizing capabilities and relationships
- **Contains**: ASCII diagrams, tables, and organized feature lists
- **Length**: ~3 pages
- **Key sections**:
  - FieldSceneState capabilities diagram
  - GPU rendering pipeline flow
  - Domain models hierarchy
  - Backend architecture
  - Camera mathematics
  - Feature readiness table

### 4. **PROJECT_EXPLORATION_REPORT.md** 📋 ORIGINAL ANALYSIS
- **Best for**: Historical context and original findings
- **Contains**: Initial project structure analysis
- **Note**: Some newer findings in ARCHITECTURE_EXPLORATION.md

### 5. **FILE_STRUCTURE.md** 📁 REFERENCE
- **Best for**: File mapping and organization
- **Contains**: Directory tree with file descriptions

---

## Quick Start Path

### If you have 5 minutes:
1. Read: **QUICK_REFERENCE.md** (File Structure + Core Components sections)
2. Focus: Figure out which files contain what

### If you have 15 minutes:
1. Read: **QUICK_REFERENCE.md** (all sections)
2. Scan: **RENDER_CAPABILITIES_SUMMARY.md** (diagrams and tables)
3. Focus: Understand what rendering capabilities exist

### If you have 30+ minutes:
1. Read: **QUICK_REFERENCE.md** (overview)
2. Read: **ARCHITECTURE_EXPLORATION.md** (detailed breakdown)
3. Reference: **RENDER_CAPABILITIES_SUMMARY.md** (visual validation)
4. Deep dive: Review actual source files as needed

---

## Layer Overview

### Domain Layer (`crates/domain/src/lib.rs`)
- Pure data models: Project, EmModel, GeometryObject, Material, SolveResult
- No rendering or UI logic
- Full serialization support (JSON + MessagePack)
- **Status**: Production-ready, minimal

### Infrastructure Layer (`crates/infra/src/lib.rs`)
- Backend trait for pluggable implementations
- File I/O (.emsp format using MessagePack)
- Solver abstraction (currently PlaceholderSolver)
- Standalone and Cloud backends
- **Status**: Production-ready framework, solver needs implementation

### Render Layer (`crates/render/src/`)
- **FieldSceneState** - Main 3D visualization system
  - 5 visualization modes (Surface, Arrows, Slice, FarField, Animation)
  - Orbit camera with 7 presets
  - 4 professional colormaps
  - Real-time phase animation
- **OrbitCamera** - Full 3D camera controller
- **FieldPipeline** - GPU rendering backend (wgpu)
- **FieldMesh** - Mesh data with field values
- **PhaseAnimator** - Complex field time-domain animation
- **Status**: Complete and feature-rich

### App Layer (`crates/app/src/lib.rs`)
- Main App struct integrating all layers
- Docking layout with 3 tabs (Modeling, Result, Log)
- Ribbon toolbar state management
- File dialog handling
- Project state tracking
- **Status**: Functional, ready for UI extensions

### Components Layer (`crates/components/src/`)
- **Ribbon** - 40+ toolbar actions (AEDT-inspired)
- **Dock** - Layout management
- **Status**: Framework-ready

---

## Key Technologies

| Category | Technology | Version |
|----------|-----------|---------|
| **UI** | egui + eframe | 0.33 |
| **GPU** | wgpu | 27 |
| **Math** | glam | 0.29 |
| **Serialization** | serde + rmp-serde | 1.0 |
| **Layout** | egui_dock | 0.18 |
| **File Dialogs** | rfd | 0.15 |
| **Platform** | Native (macOS/Linux/Windows) + WASM |

---

## Core Classes & Methods

### FieldSceneState (Main Visualization)
```rust
pub struct FieldSceneState {
    pub camera: OrbitCamera,
    pub colormap: ColormapType,
    pub vis_mode: VisMode,
    pub animator: PhaseAnimator,
    // ... internal GPU resources
}

// Key methods:
init_gpu(&mut self, render_state, mesh)     // Initialize GPU
show_viewport(&mut self, ui)                 // Render 3D view
show_controls(&mut self, ui)                 // Control panel
show_colorbar(&self, ui)                     // Legend
```

### OrbitCamera
```rust
pub struct OrbitCamera {
    pub target: Vec3,
    pub distance: f32,
    pub azimuth: f32,
    pub elevation: f32,
    // ... other camera params
}

// Key methods:
rotate(dx, dy)                               // Mouse rotation
zoom(delta)                                  // Mouse wheel
pan(dx, dy)                                  // Screen-space pan
view_projection(aspect) -> Mat4              // Combined matrix
set_preset(ViewPreset)                       // 7 presets
```

### FieldMesh
```rust
pub struct FieldMesh {
    pub vertices: Vec<FieldVertex>,
    pub indices: Vec<u32>,
    pub field_range: [f32; 2],
    pub field_imag: Option<Vec<f32>>,
    pub vector_field: Option<Vec<[f32; 3]>>,
}

// Key static methods:
uv_sphere(n_lat, n_lon, radius)              // Sphere mesh
cube(subdivisions, size)                     // Cube mesh
```

### Backend Trait
```rust
pub trait Backend {
    fn save_project(&mut self, project) -> Result<(), BackendError>;
    fn load_project(&self, id) -> Result<Project, BackendError>;
    fn solve(&self, project) -> Result<SolveResult, BackendError>;
    fn mode(&self) -> RunMode;
}
```

---

## Visualization Modes

| Mode | Purpose | Geometry |
|------|---------|----------|
| **Surface** | Colormap on field values | UV Sphere |
| **Arrows** | Vector field visualization | Subdivided Cube |
| **Slice** | 2D slice through volume | Configurable plane |
| **FarField** | 3D radiation pattern | Spherical mesh |
| **Animation** | Phase-swept complex field | UV Sphere (time-varying) |

---

## File Organization Reference

```
Key Render Files:
  scene.rs            ← FieldSceneState (START HERE)
  camera.rs           ← OrbitCamera
  mesh_data.rs        ← FieldVertex, FieldMesh, ArrowInstance
  field_pipeline.rs   ← GPU rendering backend
  animation.rs        ← PhaseAnimator
  colormap.rs         ← 4 colormaps
  arrow_pipeline.rs   ← Arrow rendering
  
Key Domain Files:
  domain/lib.rs       ← Project, EmModel, Material, etc.
  
Key Infra Files:
  infra/lib.rs        ← Backend trait, file I/O

Key App Files:
  app/lib.rs          ← App struct, integration
  components/ribbon.rs ← Ribbon toolbar (40+ actions)
```

---

## Integration Capabilities

### Rendering
- ✅ 5 visualization modes with real-time switching
- ✅ Orbit camera with 7 preset views
- ✅ 4 professional colormaps
- ✅ Complex field phase animation (0-360°)
- ✅ Real-time vertex updates
- ✅ Wireframe overlay
- ✅ Vector field arrows (up to 4096)
- ✅ Colorbar legend

### Data Management
- ✅ JSON + MessagePack serialization
- ✅ File I/O (.emsp format)
- ✅ Backend abstraction (Standalone, Cloud, Custom)
- ✅ Project state tracking
- ✅ Async file dialogs

### UI/UX
- ✅ Docking layout system
- ✅ Ribbon toolbar (40+ actions)
- ✅ Status/log display
- ✅ Mouse interaction (drag, scroll, presets)

---

## Data Flow

```
App (UI) 
  → Project (domain model)
  → Backend (file I/O, solver)
  → SolveResult (field data)
  → FieldMesh (GPU vertices)
  → FieldPipeline (GPU rendering)
  → Screen output
```

---

## Next Steps for Implementation

1. **Extend Domain Model** - Add simulation parameters (frequency, boundary conditions, etc.)
2. **Implement Real Solver** - Replace PlaceholderSolver with actual electromagnetic solver
3. **Field Data Pipeline** - Connect solver output to visualization
4. **Model Editor** - Add UI for creating/editing geometry and materials
5. **Parameter Panels** - Build simulation setup UI
6. **Results Analysis** - Enhance result visualization and export

---

## Document Maintenance Notes

- **Last Updated**: April 4, 2026
- **Scope**: Read-only exploration of render, domain, and infra layers
- **Analysis Depth**: Comprehensive coverage of all major components
- **Source**: Direct analysis of source code in `crates/render/src/`, `crates/domain/src/`, and `crates/infra/src/`

All findings are based on static code analysis. See actual source files for implementation details.

