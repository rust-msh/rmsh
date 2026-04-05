# EMStudio Render Crate - Exploration Index

**Date**: 2026-04-05  
**Project Root**: `/Users/alex/works/emstudio/`  
**Crate Root**: `/Users/alex/works/emstudio/crates/render/`

---

## Generated Documentation Files

Three comprehensive exploration documents have been created to help you understand the render crate:

### 1. **RENDER_CRATE_EXPLORATION.md** (770 lines)
**Most Detailed Reference**

A complete, section-by-section breakdown covering:
- Crate file structure and organization
- All key structs and data structures with field explanations
- Synthetic test data generation and patterns
- Data flow from input to GPU rendering
- FieldPipeline and FieldSceneState expectations
- Complete mesh format specifications (.msh and .emsfld)
- Scene composition and rendering architecture
- File loading integration points
- Summary tables and documentation references

**Use When**: You need the authoritative detailed reference; understanding specific data structure layouts; researching file format specifications.

**Location**: `/Users/alex/works/emstudio/RENDER_CRATE_EXPLORATION.md`

---

### 2. **RENDER_QUICK_REFERENCE.md** (380 lines)
**Developer Quick Lookup**

A structured quick-lookup guide including:
- File paths and key structs with hierarchical trees
- Synthetic data generators (table format)
- Simplified data flow diagrams
- File format specs in condensed form
- How to add real data loaders (step-by-step)
- GPU resource sizes and requirements
- Camera controls and animation API
- Colormap options
- Integration checklist
- Common tasks with code snippets

**Use When**: You need a quick lookup during development; integrating new features; remembering API signatures.

**Location**: `/Users/alex/works/emstudio/RENDER_QUICK_REFERENCE.md`

---

### 3. **RENDER_DATA_FLOW_VISUAL.md** (450 lines)
**Architectural Visual Guide**

Comprehensive ASCII diagrams and flow charts:
- Complete end-to-end data flow diagram (user input to screen)
- Data structure hierarchy tree (state ownership)
- GPU pipeline architecture (pipelines, bind groups, buffers)
- GPU buffer usage breakdown
- Detailed render pass sequence (frame-by-frame)
- Synthetic data generation flow
- Mode switching logic flowchart
- File format to struct mapping

**Use When**: Understanding overall architecture; visualizing data flow; explaining to others; debugging rendering issues.

**Location**: `/Users/alex/works/emstudio/RENDER_DATA_FLOW_VISUAL.md`

---

## Source Files Referenced

### Core Render Crate Files
```
crates/render/src/
├── lib.rs                    → Public exports
├── scene.rs                  → FieldSceneState, UI logic (577 lines)
├── field_pipeline.rs         → GPU rendering pipeline (632 lines)
├── arrow_pipeline.rs         → Instanced arrow rendering (~150 lines)
├── mesh_data.rs              → FieldMesh, FieldVertex, generators (324 lines)
├── animation.rs              → PhaseAnimator (53 lines)
├── slice.rs                  → Slice plane generation (110+ lines)
├── far_field.rs              → Far-field pattern mesh (~100 lines)
├── camera.rs                 → OrbitCamera controller (~100 lines)
├── colormap.rs               → Colormap LUTs (~150 lines)
└── field_shader.wgsl         → WGSL shaders (not detailed in docs)
```

### Documentation Files
```
docs/
├── em-result-file-formats.md             → .msh and .emsfld format specs
├── em-result-visualization-design.md     → Visualization architecture
```

---

## Key Findings Summary

### Architecture: 5-Layer Stack
```
1. User Interaction Layer
   └─ Mouse/keyboard → egui UI

2. State Management Layer
   └─ FieldSceneState → camera, colormap, mode, animation

3. GPU Command Recording Layer
   └─ FieldSceneCallback → uniforms, render commands

4. GPU Rendering Layer
   └─ FieldPipeline → 5 render pipelines + offscreen FB

5. Integration Layer
   └─ egui_wgpu → screen output
```

### Data Flow: 6-Step Process
```
1. User Input (mouse, UI) → Camera/state update
2. State Dirty Check → Update GPU uniforms/mesh if needed
3. Build FieldUniforms → MVP, field range, camera position
4. GPU Render Pass → Offscreen framebuffer with depth
5. Blit Pass → Full-screen triangle to egui's render target
6. egui Composite → Final screen display
```

### Synthetic Data: 5 Generators
| Generator | Vertices | Use Case |
|-----------|----------|----------|
| UV Sphere | ~2K | Surface + Animation modes |
| Subdivided Cube | ~4K | Arrows mode (vortex field) |
| Slice Plane | ~1.7K | Slice mode (volume field) |
| Far-Field Pattern | ~7.4K | FarField mode (antenna gain) |
| Arrow Base Mesh | ~100 | ArrowPipeline instancing base |

### GPU Resources: Typical 10K-Vertex Scene
- Vertex buffer: 280 KB
- Index buffers: 240 KB
- Colormap: 1 KB
- Uniform buffer: 80 bytes
- Offscreen FB (1920×1080): ~33 MB
- **Total VRAM**: ~33.5 MB (dominated by framebuffer)

### File Formats: Two Binary Types
1. **.msh** (Gmsh MSH 4.1): Mesh nodes, elements, geometry
   - Random access via entity blocks
   - Supports triangle, tetrahedron elements
   - Typical: 10K nodes ≈ 1 MB

2. **.emsfld** (EMStudio Binary): Frequency-domain field solutions
   - Complex f64 or f32 vector fields
   - Random access per frequency point
   - Typical: 10K nodes × 301 frequencies ≈ 144 MB

---

## What's NOT in the Render Crate

❌ **File Loaders**: No `.msh` or `.emsfld` parsers/readers
  - Integration point: Implement in separate `emstudio-data-loaders` crate

❌ **Quantity Expressions**: No dB/phase/magnitude computation
  - Would be in data mapping layer

❌ **Result Data Caching**: No persistent storage of loaded data
  - Would need LRU cache for large files

❌ **Mesh Quality Metrics**: No visualization of element quality
  - Would extend to MeshRenderer pipeline

---

## How to Use These Documents

### Scenario 1: "I want to understand the overall architecture"
**→ Read**: RENDER_DATA_FLOW_VISUAL.md § 1-2 (diagrams + hierarchy)
**→ Then**: RENDER_CRATE_EXPLORATION.md § 2-3 (data structures)

### Scenario 2: "How do I add a new visualization mode?"
**→ Read**: RENDER_QUICK_REFERENCE.md § Mode Switching
**→ Reference**: RENDER_DATA_FLOW_VISUAL.md § 7 (switching logic)
**→ Code template**: Check scene.rs::apply_mode_switch()

### Scenario 3: "I need to load real .msh/.emsfld files"
**→ Read**: RENDER_CRATE_EXPLORATION.md § 6-8 (formats + integration)
**→ Plan**: RENDER_QUICK_REFERENCE.md § How to Add Real Data Loaders
**→ Template code**: In § 8 of QUICK_REFERENCE.md

### Scenario 4: "What do I need to do to integrate with FieldSceneState?"
**→ Read**: RENDER_CRATE_EXPLORATION.md § 5 (expectations)
**→ Reference**: RENDER_DATA_FLOW_VISUAL.md § 8 (struct mapping)
**→ Code**: FieldSceneState::init_gpu() in scene.rs

### Scenario 5: "Debug: meshes not rendering, what's the data flow?"
**→ Follow**: RENDER_DATA_FLOW_VISUAL.md § 1 (complete flow)
**→ Check**: GPU resources in § 4 (buffers, textures)
**→ Verify**: Render pass sequence in § 5 (command recording)

---

## Cross-Document Index

### FieldVertex
- **EXPLORATION**: § 2.1 (definition + buffer layout)
- **QUICK_REF**: Shown in struct tree
- **VISUAL**: § 2 (in hierarchy), § 4 (in buffers)

### FieldMesh
- **EXPLORATION**: § 2.1 + § 3 (generators)
- **QUICK_REF**: § Synthetic Data Generators (table)
- **VISUAL**: § 2 (hierarchy), § 6 (generation flow)

### FieldPipeline
- **EXPLORATION**: § 2.2 (complete breakdown)
- **QUICK_REF**: § GPU Pipeline struct tree
- **VISUAL**: § 3 (GPU pipeline architecture), § 4 (buffers)

### FieldSceneState
- **EXPLORATION**: § 2.3 (fields + methods)
- **QUICK_REF**: § Scene State struct tree
- **VISUAL**: § 2 (hierarchy), § 1 (in data flow)

### Data Flow
- **EXPLORATION**: § 4 (detailed pipeline)
- **QUICK_REF**: § Data Flow (simplified)
- **VISUAL**: § 1 (complete ASCII diagram)

### File Formats
- **EXPLORATION**: § 6 (complete specs with code examples)
- **QUICK_REF**: § File Format Specs (condensed)
- **VISUAL**: § 8 (struct mapping diagrams)

### Synthetic Data
- **EXPLORATION**: § 3 (all generators detailed)
- **QUICK_REF**: § Synthetic Data Generators (table)
- **VISUAL**: § 6 (generation flow)

---

## Document Statistics

| Document | Lines | Sections | Tables | Diagrams | Code |
|----------|-------|----------|--------|----------|------|
| EXPLORATION | 770 | 12 | 8 | 1 ASCII | Yes |
| QUICK_REF | 380 | 15 | 12 | 5 simple | Many snippets |
| VISUAL | 450 | 8 | 5 | 15 ASCII | Few |

**Total**: 1,600 lines of documentation covering all aspects of the render crate

---

## Recommendations for Next Steps

1. **Read VISUAL first** (450 lines, mostly diagrams)
   - Understand overall architecture in 30 minutes

2. **Skim QUICK_REF** (380 lines, organized reference)
   - Familiarize yourself with APIs and quick lookups

3. **Deep dive EXPLORATION** (770 lines, detailed reference)
   - Reference during implementation

4. **Keep QUICK_REF open** during coding
   - Use for quick struct/function lookups

---

## Questions Answered by These Docs

✅ "What files are in the render crate?"  
✅ "How is data structured for GPU rendering?"  
✅ "What does FieldPipeline expect as input?"  
✅ "How are synthetic test meshes generated?"  
✅ "What's the complete rendering pipeline flow?"  
✅ "What are .msh and .emsfld file formats?"  
✅ "Where would I hook in file loaders?"  
✅ "How do I add a new visualization mode?"  
✅ "What's the memory footprint?"  
✅ "How does animation work?"  
✅ "What GPU pipelines are there?"  
✅ "How are scenes composed?"  

---

## Document Maintenance

These documents are READ-ONLY descriptions of the codebase as it exists on 2026-04-05.

**If the codebase changes**, update these documents in the order:
1. EXPLORATION (source of truth)
2. QUICK_REF (condensed version)
3. VISUAL (derived diagrams)

---

## Contact/Notes

- All paths are absolute, relative to `/Users/alex/works/emstudio/`
- Source files examined: 11 Rust files + 2 documentation files
- Focus: Data structures, rendering pipeline, file format specifications
- Not included: Shader details (field_shader.wgsl), low-level wgpu API
- Best viewed: EXPLORATION in text editor with good wrapping, VISUAL in monospace font

