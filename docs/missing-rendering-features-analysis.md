# Missing Rendering Features: rcad-render vs emstudio-render Analysis

## Overview

After thorough analysis of design docs and codebase, there's a clear **architectural split**:
- **rcad-render** (~2,000 LOC): CAD geometry visualization (tessellation, picking, camera, selection state, display modes, grid, axes, lighting, screenshot)
- **emstudio-render** (~250 LOC): Bridge layer connecting domain geometry engine to rcad-render via egui_wgpu callbacks

> **Updated 2026-04**: Many previously missing features have been implemented. This document reflects current status.

---

## 1. Core CAD Viewport Features (Standard in AEDT)

### 1.1 Viewport Display Modes ✅ IMPLEMENTED

rcad-render now supports 4 display modes:
- **SolidWithEdges** (default) — triangles + wireframe
- **Solid** — triangles only
- **Wireframe** — edges only
- **Transparent** — semi-transparent surfaces + wireframe

Still missing:
- **Hidden line** removal with edge emphasis (HLR exists as SVG export but not realtime)
- **Material-colored display** (per-object color API exists, needs multi-body rendering)

### 1.2 Grid, Axes, and Reference System ✅ IMPLEMENTED

rcad-render now provides:
- **Coordinate axes** (XYZ triad, red/green/blue arrows with cone heads)
- **Background grid** (XZ plane, major lines at 1.0, minor at 0.2, toggleable)

Still missing:
- **Axis labels** (X/Y/Z text — needs egui overlay or texture atlas)
- **Bounding box display**
- **Measurement ruler/scale**

### 1.3 View Presets and Navigation ⚠️ PARTIAL

OrbitCamera in rcad-render has basic orbit/pan/zoom, but missing:
- **View presets** (Front, Back, Left, Right, Top, Bottom, Isometric) - mentioned in visualization design
- **Fit to bounds** / Auto-zoom functionality
- **View history** (navigation stack)
- **Keyboard shortcuts** for view controls
- **Camera animation** (smooth transitions between views)
- **Viewport split/layout** management

**Status in code**:
- `rcad-render/Camera` has basic rotate/pan/zoom
- `emstudio-render/camera.rs` has `ViewPreset` enum (line 11: `OrbitCamera, ViewPreset`) but incomplete
- Missing integration with UI ribbon

### 1.4 Object Visibility and Layers ❌ MISSING

No layer/visibility system. Missing:
- **Per-object visibility toggle** (show/hide individual objects)
- **Group visibility** (toggle groups like "Antenna", "Environment" defined in geometry)
- **Layer system** with show/hide/select
- **Selective transparency** (hide some objects while keeping others semi-transparent)
- **Bounding box view** (display only bounding boxes, not full geometry)
- **Object outline/highlight** on selection

**Impact**: Cannot manage complex models with many objects. AEDT has this as fundamental feature.

**Design doc references**:
- `em-project-file-design.md` defines object groups (lines 591, 707, 712, 725)
- Objects have color/transparency properties but no visibility control

### 1.5 Selection and Highlighting ⚠️ PARTIAL

rcad-render has:
- SelectionState with face/edge selection
- Picking functions (pick_face, pick_edge)

Missing:
- **Vertex selection** mode
- **Object-level selection** (whole solid, not just faces)
- **Selection by box/polygon** (drag-select multiple objects)
- **Named selection visualization** (highlight sets like "@GND_Bottom")
- **Selection highlighting colors/styles** (glow, outline, etc.)
- **Selection feedback** (status bar showing "Face 5 selected", etc.)

**Status**: Foundation exists but needs UI integration and expanded selection modes.

### 1.6 Lighting and Shading ❌ MISSING

No proper lighting system. Missing:
- **Ambient + Directional light** (standard Phong/Blinn-Phong)
- **Specular highlights** on materials
- **Smooth shading** (Gouraud/Phong interpolation)
- **Flat shading** option
- **Edge highlighting/silhouette** for wireframe clarity
- **Depth cueing/fog** for depth perception
- **Light control UI** (direction, intensity, ambient level)

**Impact**: Geometry looks flat and hard to understand spatially. Basic shading is standard in CAD.

**Current state**: render pipelines exist but likely use basic lighting or none at all.

### 1.7 Picking Feedback and Tooltips ⚠️ PARTIAL

rcad-render has ray-triangle picking and edge picking. Missing:
- **Hover feedback** (highlight without clicking)
- **Tooltips** showing object name on hover
- **Status bar** with picked element info
- **Picked element color/style** differentiation
- **Intersection preview** (before clicking)
- **Pick radius/tolerance** visualization

**Status**: Core picking exists, needs UI integration.

---

## 2. Simulation Results Visualization Pipeline

These features are outlined in design docs but mostly **未开始 (not started)** in the milestone list.

### 2.1 Result Data Loading Infrastructure ❌ MISSING

Design doc defines:
- `ResultDataStore` trait (em-result-visualization-design.md §8.1)
- Loaders for JSON, Touchstone, MSH (Gmsh), .emsfld (binary)
- mmap-based zero-copy field data access
- LRU caching strategy

**Status**: 
- emstudio-domain and emstudio-touchstone exist but ResultDataStore integration is NOT implemented
- MSH loader: 🔲 Not started (Milestone 7)
- .emsfld loader: 🔲 Not started (Milestone 7)
- **These are blocking all 3D field visualization**

**Code impact**: None of the field rendering pipelines can function without this.

### 2.2 GPU-Resident Field Data Management ❌ MISSING

Missing infrastructure:
- **Vertex Buffer management** for field data
- **Field value normalization** (map to [0,1] for colormap)
- **Per-vertex colormaps** based on field magnitude
- **Component extraction** (X, Y, Z, magnitude from complex field)
- **Barycentric interpolation** for field values within elements

**Related to Milestone 7**: "3D 场数据管线" (3D Field Data Pipeline)

### 2.3 2D Reporting System ❌ MISSING

Design doc defines comprehensive 2D report types, but **Milestone 6** (2D Report System) is 0% complete:

Missing components:
1. **S-Parameter Rectangular Plots**
   - Trace system (multiple curves, legend)
   - Marker/Delta Marker interaction
   - Dual Y-axis (magnitude + phase)
   - Parameter sweep overlay
   
2. **Smith Chart**
   - Isoresistance/isoreactance circles
   - S-parameter trajectory
   - Normalized impedance calculation
   
3. **Polar Radiation Plots**
   - Cut-plane selection (E-plane, H-plane)
   - Concentric dB circles
   - 3dB beamwidth markers
   
4. **Convergence Curves**
   - Pass-by-pass convergence tracking
   - Mesh growth visualization
   - Delta-S or Delta-Energy metrics
   
5. **Q3D-Specific Reports**
   - RLCG matrix vs frequency curves
   - RLCG matrix data tables with heatmaps
   - Coupling coefficient calculation
   - DC/AC comparison

**Status**: Zero implementation (Milestone 6 shows 0% progress)
**Blocking**: All 2D visualization in Results tab

---

## 3. 3D Field Overlay Visualization

Design doc extensively covers these but **Milestone 7** (3D Field Data Pipeline) is 0% complete.

### 3.1 Surface Colormap Rendering ❌ MISSING

Design: em-result-visualization-design.md §4.5.1 with full WGSL shader code

Missing:
- Field Pipeline integration with wgpu
- Vertex shader with MVP transform + field value pass-through
- Fragment shader with 1D colormap texture sampling
- Lighting (Lambertian diffuse)
- Value range normalization (auto, manual, symmetric modes)

**Status**: emstudio-render has FieldPipeline (631 LOC) but likely incomplete

### 3.2 Vector Arrow Rendering ❌ MISSING

Design: em-result-visualization-design.md §4.5.2

Missing:
- Instanced rendering system for arrows
- Arrow geometry primitive (cylinder + cone)
- Arrow scale based on field magnitude
- Color mapping for arrow magnitude
- Dense sampling visualization

**Status**: emstudio-render has ArrowPipeline (145 LOC)

### 3.3 Isosurface Extraction ❌ MISSING

Design: em-result-visualization-design.md §4.2 & 4.5.3

Missing:
- **Marching Tetrahedra algorithm** implementation
- Isosurface value selection
- Triangle mesh generation from iso-value
- GPU rendering of extracted surface

**Complexity**: Medium-high. Not trivial algorithm.

### 3.4 Slice Plane Visualization ❌ MISSING

Design: em-result-visualization-design.md §4.6

Missing:
- Plane intersection with tetrahedral mesh
- Slice mesh generation
- Barycentric interpolation of field values on slice
- Slice plane parameter UI (position, orientation)

**Status**: emstudio-render has Slice module (101 LOC) but likely stub

### 3.5 Phase Animation ❌ MISSING

Design: em-result-visualization-design.md §4.7

Missing:
- Phase stepping logic (0° to 360°)
- Real/imaginary component extraction from complex field
- Phase-dependent field evaluation: E_real = Re(E)*cos(φ) - Im(E)*sin(φ)
- Frame rate control
- Play/pause controls

**Status**: emstudio-render has animation.rs (52 LOC) - likely incomplete

### 3.6 3D Far-Field Pattern Rendering ❌ MISSING

Design: em-result-visualization-design.md §5

Missing:
- Gain-to-radius mapping (linear, logarithmic)
- Sphere deformation based on gain values
- 3D pattern surface mesh generation
- Gain-value colormap
- Pattern scale selection UI

**Status**: emstudio-render has FarFieldGen (referenced but lines unknown)

### 3.7 GPU Picking for Field Data ❌ MISSING

Design: em-result-visualization-design.md §7.2

Missing:
- Off-screen picking pass (color-coded object IDs)
- Pixel readback buffer
- Picking result query and decoding

**Status**: Design only, no implementation

### 3.8 Field Probe System ❌ MISSING

Design: em-result-visualization-design.md §7.3

Missing:
- Point-in-tetrahedron testing
- Barycentric coordinate calculation
- Field value interpolation at arbitrary 3D points
- Probe UI (click to place, result display)
- Hover value display

---

## 4. Geometry Rendering and Visualization

From **Milestone 4**: Geometry modeling is ~36% complete.

### 4.1 Solid Body Rendering ❌ MISSING

Missing from geometry rendering pipeline:
- **Proper tessellation** of BREP solids to triangles
- **Normal calculation** for smooth shading
- **Face grouping** by object (not just flat vertex list)
- **Back-face culling** for performance
- **Smooth vs flat shading** toggle

**Current state**:
- rcad-render has Tessellator.tessellate() which flattens geometry
- No per-face normal or face tracking
- No separate render groups per solid

### 4.2 STEP/STL CAD Import Visualization ❌ MISSING

Design mentions STEP import but:
- rcad-step exists (import/export)
- No integration with visualization pipeline
- No color/transparency from STEP colors
- No object tree visualization

### 4.3 Object Tree and Hierarchy Display ❌ MISSING

Missing:
- **Outliner/Model Tree** panel showing object hierarchy
- **Expand/collapse** object groups
- **Visibility toggle** in tree
- **Selection sync** (click in tree = select in 3D view)
- **Rename/property edit** in tree

### 4.4 Geometric Measurement ❌ MISSING

Missing measurement tools:
- **Distance measurement** between points/edges/faces
- **Angle measurement** between edges/faces
- **Area/volume calculation** display
- **Coordinate display** for selected elements
- **Bounding box dimensions** displayed

---

## 5. Material and Appearance Properties

### 5.1 Material Visualization ❌ MISSING

Design doc defines:
- Object colors and transparency (em-project-file-design.md §3.3.3)
- Material properties (epsilon_r, mu_r, sigma, tan_delta)

Missing implementation:
- **Per-object color rendering** (currently all objects same color)
- **Transparency/opacity** per object
- **Material-based coloring** (show material type visually)
- **Dielectric appearance** (translucent vs opaque rendering)
- **Conductor appearance** (metallic/shiny rendering)

**Current state**: Color/transparency stored in GeoObject but not used in rendering

### 5.2 Texture/Appearance Mapping ❌ MISSING

Missing:
- **Checkerboard pattern** for transparent objects (clarity)
- **Metallic texture** for conductors
- **Dielectric sheen** for insulators

---

## 6. UI Integration and Control Panels

### 6.1 Viewport Properties Panel ❌ MISSING

Right-side panel missing:
- **View mode selector** (wireframe, solid, transparent, etc.)
- **Light direction control** (direction, intensity)
- **Clipping plane controls**
- **Background color picker**
- **Grid visibility and spacing** sliders
- **Axis/origin display toggles**

### 6.2 Field Overlay Properties Panel ❌ MISSING (Partial)

Design doc §9.3 shows UI mockup, but implementation missing:
- **Field type selector** (E, H, J, Poynting, etc.)
- **Component selector** (Magnitude, X, Y, Z, Vector)
- **Frequency slider** with available frequencies
- **Phase slider** (0-360°)
- **Colormap picker** with preview
- **Value range controls** (auto, manual, symmetric)
- **Opacity slider**
- **Slice plane controls** (position, orientation)
- **Isosurface value input**
- **Animation controls** (play, pause, speed)

### 6.3 Material/Property Editor ❌ MISSING

Missing UI for:
- **Object appearance properties** (color, transparency, group)
- **Material assignment** picker
- **Coordinate system assignment** for local CS
- **Named selection** assignment UI

---

## 7. Export and Screenshot Functionality

### 7.1 3D Screenshot/PNG Export ❌ MISSING

Design doc §10.3 defines interface but missing:
- Offscreen render target creation
- Screenshot capture at arbitrary resolution
- PNG encoding and save

**Status**: Mentioned in milestones but not implemented

### 7.2 Data Export Formats ⚠️ PARTIAL

Implemented:
- S-Parameter Touchstone (emstudio-touchstone crate)

Missing:
- CSV reports with headers
- VTK export for ParaView compatibility
- SPICE netlist (Q3D equivalent circuit)

---

## 8. Performance and Optimization Features

### 8.1 Level-of-Detail (LOD) Rendering ❌ MISSING

Missing:
- Coarse/fine mesh switching based on zoom
- View distance culling
- Frustum culling
- Occlusion culling

### 8.2 Progressive Loading ❌ MISSING

For large field data:
- Streaming field data from disk
- Progressive mesh loading
- Adaptive refinement during interaction

---

## 9. Accessibility and Interactive Feedback

### 9.1 Measurement and Dimension Display ❌ MISSING

Missing:
- On-screen dimension labels
- Parametric dimension display (from design)
- Unit display and conversion

### 9.2 Query and Information Display ❌ MISSING

Missing:
- **Element info on hover** (object name, material, properties)
- **Coordinate display** in viewport
- **Field value tooltip** at cursor position
- **Status bar** with selection/mode info

---

## Summary: Feature Matrix

| Feature Category | Subcategory | Status | Complexity | Blocker? |
|---|---|---|---|---|
| **CAD Viewport Basics** | Display modes | ❌ Missing | Low | No |
| | Grid/Axes/Origin | ❌ Missing | Low-Med | **Yes** |
| | View presets | ⚠️ Partial | Low | No |
| | Visibility/Layers | ❌ Missing | Low-Med | **Yes** |
| | Selection & Highlight | ⚠️ Partial | Low-Med | No |
| | Lighting & Shading | ❌ Missing | Medium | **Yes** |
| | Object Properties UI | ❌ Missing | Low | No |
| **Geometry Rendering** | Solid body rendering | ❌ Missing | Low-Med | **Yes** |
| | Face/edge display | ⚠️ Partial | Low | No |
| | Material colors | ❌ Missing | Low | No |
| | CAD import viz | ❌ Missing | Low-Med | No |
| **2D Reports** | S-Parameter plots | ❌ Missing | Medium | **Yes** |
| | Smith charts | ❌ Missing | Medium | **Yes** |
| | Polar plots | ❌ Missing | Medium | **Yes** |
| | Convergence curves | ❌ Missing | Low | **Yes** |
| | RLCG reports | ❌ Missing | Medium | **Yes** |
| **3D Field Viz** | Result data loading | ❌ Missing | Medium | **Critical** |
| | Surface colormaps | ❌ Missing | Medium | **Critical** |
| | Vector arrows | ❌ Missing | Medium | **Critical** |
| | Isosurfaces | ❌ Missing | High | No |
| | Slice planes | ❌ Missing | Medium | **Critical** |
| | Phase animation | ❌ Missing | Medium | No |
| | Far-field patterns | ❌ Missing | Medium | No |
| | Field probes | ❌ Missing | Medium | No |
| | GPU picking | ❌ Missing | Medium | No |
| **Exports** | Screenshots | ❌ Missing | Low | No |
| | CSV reports | ❌ Missing | Low | No |
| | VTK export | ❌ Missing | Medium | No |
| | SPICE export | ❌ Missing | Low | No |

---

## Critical Path Dependencies

To achieve functional "rendering" (both CAD and EM results):

```
MUST-HAVE (blocking everything):
1. Grid/axes/origin display ← Required for CAD basic usability
2. Solid geometry rendering ← Required for any CAD work
3. View presets + camera controls ← Navigation essential
4. Result data loaders (.msh, .emsfld) ← Without this, no field viz possible
5. 2D report infrastructure ← Core post-processing feature
6. Surface colormap rendering ← First EM result type users see

THEN-HAVE (dependent on above):
7. Selection highlighting UI
8. Material colors
9. Visibility/layers
10. Lighting system
11. Advanced field viz (iso-surfaces, arrows, etc.)
```

---

## Recommendations

1. **Priority 1 (Critical)**: Grid/axes, solid rendering, result data loaders
   - These enable basic CAD + result viewing workflows
   
2. **Priority 2 (High)**: 2D reports, surface colormaps
   - These enable analysis workflows (looking at results)
   
3. **Priority 3 (Medium)**: Selection UI, visibility, lighting
   - These improve usability but don't block core workflows
   
4. **Priority 4 (Nice-to-have)**: Advanced field viz, exports, measurements
   - These are enhancements but less critical for MVP

