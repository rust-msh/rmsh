# EMStudio Render Crate - Visual Data Flow & Architecture

## 1. Complete Data Flow Diagram

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                         USER INTERACTION LAYER                               │
│  ┌─────────────────┬─────────────────┬──────────────────┬─────────────────┐  │
│  │  Mouse Input    │  UI Controls    │  Keyboard Input  │  Time/Animation │  │
│  │  (Drag, Scroll) │  (Dropdowns)    │  (None yet)      │  (Frame delta)  │  │
│  └─────────────────┴─────────────────┴──────────────────┴─────────────────┘  │
│                                      │                                        │
└──────────────────────────────────────┼────────────────────────────────────────┘
                                       │
                                       ▼
                    ┌──────────────────────────────────┐
                    │   FieldSceneState::show_viewport │
                    └──────────────────────────────────┘
                     │                                 │
        ┌────────────┼────────────┬───────────────┬────┴────────┐
        │            │            │               │             │
        ▼            ▼            ▼               ▼             ▼
   [Camera]   [Dirty Flags]  [Animation]  [Colormap]   [Uniform Build]
   rotate()   check & clear   tick() &    update()    view_projection()
   pan()                       apply()     create LUT  eye_position()
   zoom()                                            light_direction()
   
        │            │            │               │             │
        └────────────┼────────────┴───────────────┴─────────────┘
                     │
                     ▼
        ┌────────────────────────────────┐
        │  Build FieldUniforms (80 bytes)│
        │  - MVP matrix                  │
        │  - Eye position                │
        │  - Light direction             │
        │  - Field min/max               │
        │  - Opacity                     │
        └────────────────────────────────┘
                     │
                     ▼
        ┌────────────────────────────────────┐
        │  Create FieldSceneCallback         │
        │  {uniforms, show_wireframe, ...}   │
        └────────────────────────────────────┘
                     │
                     ▼
        ┌────────────────────────────────────┐
        │  ui.painter().add(Callback)        │
        │  ↓ egui_wgpu integration           │
        └────────────────────────────────────┘
                     │
        ┌────────────┴────────────┐
        │                         │
        ▼                         ▼
┌────────────────────┐   ┌──────────────────┐
│ prepare()          │   │ paint()          │
├────────────────────┤   ├──────────────────┤
│ 1. Resize if      │   │ 1. Get pipeline  │
│    needed         │   │ 2. Call blit()   │
│ 2. Update uniforms│   │ 3. Output to     │
│ 3. Record render  │   │    egui FB       │
│    commands       │   └──────────────────┘
│ 4. Return empty   │        │
│    (pre-recorded) │        ▼
└────────────────────┘   ┌──────────────────┐
        │                │ egui_wgpu renders│
        │                │ to screen        │
        │                └──────────────────┘
        │
        ▼
┌────────────────────────────────────────────┐
│  FieldPipeline::render_scene()             │
├────────────────────────────────────────────┤
│ 1. begin_render_pass to offscreen FB       │
│    color: Rgba8UnormSrgb                   │
│    depth: Depth32Float                     │
│ 2. Clear: color=(0.12, 0.12, 0.15, 1.0)   │
│    depth: 1.0                              │
│                                             │
│ 3. Set scene_pipeline                      │
│    Set bind groups (uniforms, colormap)    │
│    Set vertex buffer (FieldVertex)         │
│    Set index buffer (triangle indices)     │
│    draw_indexed(0..num_indices)            │
│       ↓ field_shader.wgsl::vs_main         │
│       ↓ field_shader.wgsl::fs_main         │
│                                             │
│ 4. [Optional] Wireframe overlay            │
│    Set wire_pipeline                       │
│    Set index buffer (edge indices)         │
│    draw_indexed(0..num_wire_indices)       │
│       ↓ field_shader.wgsl::fs_wire         │
│                                             │
│ 5. [Optional] Arrow instancing             │
│    arrow_pipeline.draw(...)                │
│       ↓ field_shader.wgsl::vs_arrow        │
│       ↓ field_shader.wgsl::fs_main         │
│                                             │
│ 6. End render pass                         │
│    (Output: color_texture with depth)      │
└────────────────────────────────────────────┘
        │
        ▼
┌────────────────────────────────────────────┐
│  FieldPipeline::blit()                     │
├────────────────────────────────────────────┤
│ Render pass to egui's target               │
│ Set blit_pipeline                          │
│ Set bind group (color_texture, sampler)    │
│ draw(0..3, 0..1)  ← Full-screen triangle   │
│    ↓ field_shader.wgsl::vs_blit            │
│    ↓ field_shader.wgsl::fs_blit            │
│    Samples color_texture linearly          │
│    Outputs to egui framebuffer             │
└────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────┐
│  egui Composites on Screen       │
│  (Panel background + 3D view)    │
└──────────────────────────────────┘
```

---

## 2. Data Structure Hierarchy

```
FieldSceneState (UI state + rendering control)
│
├─ camera: OrbitCamera
│   ├─ target: Vec3
│   ├─ distance: f32
│   ├─ azimuth: f32
│   ├─ elevation: f32
│   ├─ fov_y: f32
│   ├─ near/far: f32
│   └─ [Methods: rotate, zoom, pan, view_projection]
│
├─ colormap: ColormapType (enum)
│   ├─ Rainbow
│   ├─ Viridis
│   ├─ CoolWarm
│   └─ Grayscale
│       └─ generate_lut(256) → Vec<[u8; 4]>
│
├─ animator: PhaseAnimator
│   ├─ phase_deg: f32 (0-360)
│   ├─ playing: bool
│   ├─ speed_deg_per_sec: f32
│   └─ apply(real, imag) → Vec<f32>
│
├─ sphere_mesh: Option<FieldMesh>
├─ cube_mesh: Option<FieldMesh>
├─ far_field_mesh: Option<FieldMesh>
│
└─ [GPU Resources in egui_wgpu::CallbackResources]
   └─ FieldSceneResources
      ├─ pipeline: FieldPipeline
      │   ├─ scene_pipeline: wgpu::RenderPipeline
      │   ├─ wire_pipeline: wgpu::RenderPipeline
      │   ├─ vertex_buf: wgpu::Buffer (FieldVertex)
      │   ├─ index_buf: wgpu::Buffer (u32 triangles)
      │   ├─ wire_index_buf: wgpu::Buffer (u32 edges)
      │   ├─ uniform_buf: wgpu::Buffer (FieldUniforms)
      │   ├─ colormap_texture: wgpu::Texture (256x1 RGBA8)
      │   ├─ color_view: wgpu::TextureView (offscreen FB)
      │   ├─ depth_view: wgpu::TextureView (offscreen depth)
      │   ├─ blit_pipeline: wgpu::RenderPipeline
      │   └─ blit_bind_group: wgpu::BindGroup
      │
      └─ arrow_pipeline: Option<ArrowPipeline>
         ├─ pipeline: wgpu::RenderPipeline
         ├─ base_vertex_buf: wgpu::Buffer (arrow mesh)
         ├─ base_index_buf: wgpu::Buffer
         └─ instance_buf: wgpu::Buffer (ArrowInstance)


FieldMesh (Render-ready mesh)
│
├─ vertices: Vec<FieldVertex>
│  └─ FieldVertex (28 bytes each)
│     ├─ position: [f32; 3]
│     ├─ normal: [f32; 3]
│     └─ field_value: f32
│
├─ indices: Vec<u32> (triangle indices: 3 per tri)
├─ wire_indices: Vec<u32> (edge indices: 2 per edge)
├─ field_range: [f32; 2] (min, max for colormap)
├─ field_imag: Option<Vec<f32>> (for animation)
└─ vector_field: Option<Vec<[f32; 3]>> (for arrows)


FieldUniforms (80 bytes per frame)
│
├─ mvp: [f32; 16] (Model-View-Projection matrix)
├─ eye_pos: [f32; 3] (camera eye position)
├─ _pad0: f32 (alignment)
├─ light_dir: [f32; 3] (light direction)
├─ _pad1: f32 (alignment)
├─ field_min: f32 (colormap min)
├─ field_max: f32 (colormap max)
├─ opacity: f32 (alpha blending)
└─ _pad2: f32 (alignment to 80 bytes = 5 × 16)
```

---

## 3. GPU Pipeline Architecture

```
┌──────────────────────────────────────────────────────────────┐
│                    GPU PIPELINE SETUP                         │
├──────────────────────────────────────────────────────────────┤
│                                                               │
│  Shader Module (WGSL)                                        │
│  ├─ vs_main (scene vertex shader)                            │
│  ├─ fs_main (scene fragment shader with colormap lookup)     │
│  ├─ fs_wire (wireframe fragment shader - white edges)        │
│  ├─ vs_arrow (arrow instancing vertex shader)                │
│  ├─ vs_blit (full-screen triangle VS)                        │
│  └─ fs_blit (blit fragment shader - texture sampling)        │
│                                                               │
│  Bind Group Layouts                                          │
│  ├─ SceneBGL (scene_bind_group_layout)                       │
│  │  ├─ Binding 0: uniform_buf (FieldUniforms)               │
│  │  ├─ Binding 1: colormap_texture (1D LUT)                 │
│  │  └─ Binding 2: colormap_sampler (Linear)                 │
│  │                                                            │
│  └─ BlitBGL (blit_bind_group_layout)                         │
│     ├─ Binding 0: color_texture (offscreen FB)              │
│     └─ Binding 1: blit_sampler (Linear)                     │
│                                                               │
│  Pipelines (RenderPipeline)                                  │
│  ├─ scene_pipeline                                           │
│  │  ├─ Topology: TriangleList                                │
│  │  ├─ Vertex input: FieldVertex::buffer_layout()           │
│  │  ├─ Fragment output: Rgba8UnormSrgb + Alpha blending     │
│  │  ├─ Depth: Depth32Float (write enabled, Less compare)    │
│  │  └─ Outputs to: color_view (offscreen FB)                │
│  │                                                            │
│  ├─ wire_pipeline                                            │
│  │  ├─ Topology: LineList                                    │
│  │  ├─ Depth: Depth32Float (write disabled, LessEqual)      │
│  │  ├─ Depth bias: -2 (constant), -1.0 (slope)              │
│  │  └─ Fragment output: White color                         │
│  │                                                            │
│  ├─ arrow_pipeline                                           │
│  │  ├─ Topology: TriangleList                                │
│  │  ├─ Vertex input: FieldVertex + ArrowInstance            │
│  │  ├─ Instance step mode: Instance                          │
│  │  └─ Outputs to: color_view (offscreen FB, composited)   │
│  │                                                            │
│  └─ blit_pipeline                                            │
│     ├─ Topology: TriangleList                                │
│     ├─ No vertex input (generated full-screen triangle)      │
│     ├─ Fragment output: target_format (egui's format)        │
│     └─ No depth                                              │
│                                                               │
└──────────────────────────────────────────────────────────────┘
```

---

## 4. GPU Buffer Usage

```
┌──────────────────────────────────────────────────────────┐
│              GPU BUFFERS & TEXTURES                       │
├──────────────────────────────────────────────────────────┤
│                                                           │
│ VERTEX BUFFERS                                           │
│ ├─ vertex_buf (dynamic, COPY_DST)                        │
│ │  └─ Contains: Vec<FieldVertex>                         │
│ │     Size: num_vertices × 28 bytes                      │
│ │     Ex: 10K vertices = 280 KB                          │
│ │     Updated: queue.write_buffer() per frame (if anim)  │
│ │                                                         │
│ └─ [ArrowPipeline] base_vertex_buf (static)             │
│    └─ Arrow geometry: shaft + cone                       │
│       Size: ~100-200 bytes                               │
│                                                           │
│ INDEX BUFFERS                                            │
│ ├─ index_buf (static)                                    │
│ │  └─ Triangle indices (u32, 3 per triangle)             │
│ │     Ex: 30K triangles = 120 KB                         │
│ │                                                         │
│ ├─ wire_index_buf (static)                               │
│ │  └─ Edge indices (u32, 2 per edge)                     │
│ │     Ex: 30K edges = 120 KB                             │
│ │                                                         │
│ └─ [ArrowPipeline] base_index_buf (static)              │
│    └─ Arrow geometry indices: ~300 bytes                 │
│                                                           │
│ UNIFORM BUFFERS                                          │
│ └─ uniform_buf (dynamic, COPY_DST)                       │
│    └─ FieldUniforms: 80 bytes                            │
│       Updated: queue.write_buffer() every frame          │
│                                                           │
│ TEXTURES                                                 │
│ ├─ colormap_texture (1D, 256×1, Rgba8UnormSrgb)         │
│ │  └─ Size: 1 KB (4 bytes × 256)                         │
│ │     Updated: When colormap changes (rare)              │
│ │                                                         │
│ ├─ color_texture (Rgba8UnormSrgb) — OFFSCREEN BUFFER    │
│ │  └─ Size: viewport_width × viewport_height            │
│ │     Ex: 1920×1080 = 8.3 MB                             │
│ │     Resized to 16-byte boundary (1936×1088 ≈ 8.4 MB)  │
│ │     Reused: every frame                                │
│ │                                                         │
│ ├─ depth_texture (Depth32Float) — OFFSCREEN BUFFER      │
│ │  └─ Size: 4 bytes × viewport_width × viewport_height   │
│ │     Ex: 1920×1080 = 8.3 MB                             │
│ │     Attachments: color_view, depth_view                │
│ │                                                         │
│ └─ [ArrowPipeline] instance_buf (dynamic, COPY_DST)     │
│    └─ ArrowInstance data: 32 bytes each                  │
│       Ex: 5K arrows = 160 KB                             │
│                                                           │
│ SAMPLERS                                                 │
│ ├─ colormap_sampler (Linear filtering)                   │
│ └─ blit_sampler (Linear filtering)                       │
│                                                           │
└──────────────────────────────────────────────────────────┘
```

---

## 5. Rendering Pass Sequence

```
┌─────────────────────────────────────────────────────────────┐
│                   RENDER PASS SEQUENCE                       │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│ Each Frame:                                                 │
│                                                              │
│ 1. PREPARE PHASE (FieldSceneCallback::prepare)             │
│    ├─ Resize offscreen FB if viewport changed              │
│    ├─ Update uniform_buf with current FieldUniforms        │
│    ├─ Create command encoder                               │
│    └─ Record render pass to offscreen FB                   │
│                                                              │
│ 2. SCENE RENDER PASS (to offscreen color_texture)          │
│    ├─ Clear: color=(0.12, 0.12, 0.15, 1.0), depth=1.0     │
│    ├─ Attachment: color_view (Rgba8UnormSrgb)              │
│    ├─ Attachment: depth_view (Depth32Float)                │
│    │                                                         │
│    ├─ FILLED MESH DRAW CALL                                │
│    │  ├─ set_pipeline(&scene_pipeline)                     │
│    │  ├─ set_bind_group(0, &scene_bind_group)              │
│    │  ├─ set_vertex_buffer(0, vertex_buf)                  │
│    │  ├─ set_index_buffer(index_buf, Uint32)               │
│    │  └─ draw_indexed(0..num_indices)                      │
│    │     └─ Executes: num_indices / 3 triangles            │
│    │        (Vertex: vs_main, Fragment: fs_main)          │
│    │                                                         │
│    ├─ [IF show_wireframe] WIREFRAME DRAW CALL              │
│    │  ├─ set_pipeline(&wire_pipeline)                      │
│    │  ├─ set_vertex_buffer(0, vertex_buf)                  │
│    │  ├─ set_index_buffer(wire_index_buf, Uint32)          │
│    │  └─ draw_indexed(0..num_wire_indices)                 │
│    │     └─ Executes: num_wire_indices / 2 lines           │
│    │        (Fragment: fs_wire = white)                    │
│    │                                                         │
│    └─ [IF show_arrows] ARROW DRAW CALL (Instanced)        │
│       ├─ arrow_pipeline.draw(&mut rpass, ...)              │
│       ├─ set_pipeline(&arrow_pipeline.pipeline)            │
│       ├─ set_vertex_buffer(0, base_vertex_buf)             │
│       ├─ set_vertex_buffer(1, instance_buf)                │
│       ├─ set_index_buffer(base_index_buf, Uint32)          │
│       └─ draw_indexed_instanced(                           │
│            0..num_base_indices,                             │
│            0..num_instances                                 │
│          )                                                  │
│          └─ Executes: num_instances × num_base_indices    │
│             (Vertex: vs_arrow, Fragment: fs_main)         │
│                                                              │
│    [End render pass]                                        │
│    (Output: color_texture, depth_texture)                  │
│                                                              │
│ 3. PAINT PHASE (FieldSceneCallback::paint)                │
│    └─ Receive render_pass from egui_wgpu                  │
│                                                              │
│ 4. BLIT RENDER PASS (to egui's render pass)               │
│    ├─ set_pipeline(&blit_pipeline)                         │
│    ├─ set_bind_group(0, &blit_bind_group)                  │
│    ├─ draw(0..3, 0..1)  ← Full-screen triangle            │
│    │  └─ Executes: 1 instance, 3 vertices                 │
│    │     (Vertex: vs_blit generates pos, Fragment: fs_blit)│
│    │                                                         │
│    └─ [End render pass]                                    │
│       (Output: egui framebuffer)                           │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

---

## 6. Synthetic Data Generation Flow

```
┌────────────────────────────────────────────────────────────┐
│           SYNTHETIC DATA GENERATION (Demo Data)             │
├────────────────────────────────────────────────────────────┤
│                                                             │
│ SPHERE (VisMode::Surface, VisMode::Animation)             │
│ ├─ FieldMesh::uv_sphere(32, 64, 1.0)                       │
│ ├─ Grid: 33 × 65 = 2145 vertices                           │
│ ├─ Vertices:                                               │
│ │  for lat in 0..33:                                       │
│ │    θ = π × lat / 32                                      │
│ │    for lon in 0..65:                                     │
│ │      φ = 2π × lon / 64                                   │
│ │      x = sin(θ) cos(φ), y = cos(θ), z = sin(θ) sin(φ)   │
│ │      field = sin(3φ) cos(2θ)                            │
│ │      imag = cos(2φ) sin(3θ)                             │
│ ├─ Indices: ~4K triangles                                  │
│ └─ Wire: ~6K edges                                         │
│                                                             │
│ CUBE (VisMode::Arrows)                                     │
│ ├─ FieldMesh::cube(10, 1.0)                                │
│ ├─ 6 faces × 11×11 grid = 726 vertices per face            │
│ ├─ Total: ~4356 vertices                                   │
│ ├─ Vector field: synthetic_vortex(x,y,z)                   │
│ │  vx = -0.8z + 0.3y                                       │
│ │  vy = sqrt(x²+z²).sin() × 0.5                            │
│ │  vz = 0.8x - 0.3y                                        │
│ ├─ Field value: magnitude of vector                        │
│ └─ Arrows: subsampled every 3rd vertex                     │
│                                                             │
│ SLICE (VisMode::Slice)                                     │
│ ├─ generate_slice_mesh(SliceAxis::Z, z_val, 1.0, 40)      │
│ ├─ Grid: 41 × 41 = 1681 vertices                           │
│ ├─ Plane position: z = user slider [-1, +1]               │
│ ├─ Field fn: synthetic_volume_field(x,y,z)                 │
│ │  r = sqrt(x²+y²+z²)                                      │
│ │  field = sin(4r) / r                                     │
│ ├─ Indices: ~1600 triangles                                │
│ └─ Generated on demand per frame                           │
│                                                             │
│ FAR-FIELD (VisMode::FarField)                             │
│ ├─ generate_pattern_mesh(60, 120, patch_gain)              │
│ ├─ Grid: 61 × 121 = 7381 vertices                          │
│ ├─ Parametric surface (θ, φ):                              │
│ │  r = (gain_dbi + 30) / 30, clamped [0.05, 2.0]           │
│ │  x = r sin(θ) cos(φ)                                     │
│ │  y = r cos(θ)                                            │
│ │  z = r sin(θ) sin(φ)                                     │
│ ├─ gain_fn: patch_gain(θ, φ) — synthetic antenna           │
│ ├─ Indices: ~7K triangles                                  │
│ └─ Pre-generated once                                      │
│                                                             │
│ ARROWS (Vector visualization on cube mesh)                 │
│ ├─ cube_mesh.generate_arrows(every_n=3)                    │
│ ├─ Filter vertices: take every 3rd                         │
│ ├─ ArrowInstance {pos, dir, mag}:                          │
│ │  pos = vertex position                                   │
│ │  dir = normalize(vector_field[i])                        │
│ │  mag = ||vector_field[i]|| / max_mag                     │
│ ├─ Output: ~500-1000 arrow instances                       │
│ └─ Rendered via arrow_pipeline (instancing)                │
│                                                             │
└────────────────────────────────────────────────────────────┘
```

---

## 7. Mode Switching Logic

```
┌──────────────────────────────────────────────────┐
│      VISUALIZATION MODE SWITCHING                │
├──────────────────────────────────────────────────┤
│                                                  │
│ User selects new vis_mode                       │
│ ↓                                                │
│ FieldSceneState.vis_mode = new_mode             │
│ FieldSceneState.mode_dirty = true               │
│                                                  │
│ Next frame in show_viewport():                  │
│ ↓                                                │
│ if mode_dirty:                                  │
│   mode_dirty = false                            │
│   apply_mode_switch()                           │
│     ↓                                            │
│     match vis_mode {                            │
│       Surface | Animation:                      │
│         mesh = sphere_mesh (pre-generated)      │
│         [Animation: compute envelope_range()]   │
│         ↓                                        │
│       Arrows:                                   │
│         mesh = cube_mesh (pre-generated)        │
│         generate_arrows(3)                      │
│         ↓ ArrowPipeline::upload_instances()     │
│         ↓                                        │
│       Slice:                                    │
│         Generate new slice_mesh on-the-fly      │
│         (using current slice_z value)           │
│         ↓                                        │
│       FarField:                                 │
│         mesh = far_field_mesh (pre-generated)   │
│         ↓                                        │
│     }                                           │
│     ↓                                            │
│     pipeline.swap_mesh(device, &mesh)           │
│       ├─ device.create_buffer_init(&verts)     │
│       ├─ device.create_buffer_init(&indices)   │
│       └─ Update: num_indices, num_wire_indices  │
│                                                  │
└──────────────────────────────────────────────────┘
```

---

## 8. File Format Mapping to Structs

```
┌──────────────────────────────────────────────────────┐
│      .msh (Gmsh MSH 4.1) → FieldMesh                │
├──────────────────────────────────────────────────────┤
│                                                      │
│ .msh File Sections          FieldMesh Fields        │
│ ┌──────────────────────┐    ┌─────────────────────┐│
│ │ $Nodes               │───→│ vertices:           ││
│ │ {coords, tags}       │    │   FieldVertex[]     ││
│ │                      │    │   pos, normal, fv   ││
│ └──────────────────────┘    └─────────────────────┘│
│                                                      │
│ ┌──────────────────────┐    ┌─────────────────────┐│
│ │ $Elements            │───→│ indices:            ││
│ │ {connectivity}       │    │   u32[]             ││
│ │ (type 2=triangle,    │    │   (3 per triangle)  ││
│ │  type 4=tetrahedron) │    │                     ││
│ └──────────────────────┘    └─────────────────────┘│
│                                                      │
│ Derived Data:                                       │
│ • Normals: recompute from element connectivity     │
│ • Wireframe: edges from element faces              │
│ • Field range: [min(field_values), max(...)]       │
│                                                      │
└──────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────┐
│      .emsfld → FieldMesh + metadata                 │
├──────────────────────────────────────────────────────┤
│                                                      │
│ .emsfld File Sections      Usage                    │
│ ┌──────────────────────┐                            │
│ │ Header (128 bytes)   │                            │
│ │ {num_nodes,          │──→ Allocate vertex vec    │
│ │  num_frequencies,    │    Load reference mesh    │
│ │  field_type,         │    Determine data type    │
│ │  mesh_file}          │                            │
│ └──────────────────────┘                            │
│                                                      │
│ ┌──────────────────────┐                            │
│ │ Frequency Table      │──→ Populate frequencies   │
│ │ (f64[num_freq])      │    UI frequency selector  │
│ └──────────────────────┘                            │
│                                                      │
│ ┌──────────────────────┐                            │
│ │ Field Data Blocks    │──→ field_value per vertex │
│ │ (indexed by freq)    │    field_imag (optional)  │
│ │ {complex f64 vecs}   │    vector_field (optional)│
│ └──────────────────────┘                            │
│                                                      │
│ Result:                                             │
│ FieldMesh {                                         │
│   vertices: [...],     ← coords + field_value      │
│   indices: [...],      ← from reference .msh       │
│   field_range: [min, max],                         │
│   field_imag: [...],   ← for animation             │
│   vector_field: [...], ← for arrows                │
│ }                                                   │
│                                                      │
└──────────────────────────────────────────────────────┘
```

---

## Summary: You Now Understand

✅ **File Organization**: 11 source files, modular design  
✅ **Data Structures**: FieldVertex, FieldMesh, FieldUniforms, FieldSceneState  
✅ **GPU Pipelines**: 5 render pipelines (scene, wireframe, arrows, blit) + framebuffer management  
✅ **Synthetic Data**: 5 generator functions for demo visualizations  
✅ **Animation**: Phase-swept complex field via real+imaginary parts  
✅ **File Formats**: .msh (Gmsh) and .emsfld (EMStudio binary) specifications  
✅ **Data Flow**: User input → UI state → GPU upload → Render → Blit → Screen  
✅ **Integration Points**: Where file loaders and data sources would plug in  

**Next Steps**: Implement `.msh` parser and `.emsfld` loader in a separate crate, then integrate with FieldSceneState.

