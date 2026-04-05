# M7 & M8 Quick Reference Guide

## Data Structures at a Glance

### Visualization Configuration
```
FieldOverlay (5 types)
├─ SurfacePlot (cloud map)
├─ SlicePlot (plane intersection)
├─ VectorPlot (arrow glyphs)
├─ IsosurfacePlot (contour surface)
└─ AnimatedPlot (phase sweep)

FieldQuantityConfig
├─ field_type: FieldType (E, H, J, VolumeCurrent, ChargeDistribution, OhmicLoss)
├─ component: FieldComponent (Magnitude, X, Y, Z, Vector)
├─ frequency: String
└─ phase_deg: f64

FieldVisualConfig
├─ colormap: ColormapType
├─ range: ValueRange (Auto, Manual, Symmetric)
├─ scale: ScaleType (Linear, Logarithmic)
├─ opacity: f64
└─ legend_position: LegendPosition
```

### Camera & Interaction
```
OrbitCamera (spherical coordinates)
├─ target: [f64; 3]
├─ distance: f64
├─ azimuth_deg: f64 (rotation around Y)
├─ elevation_deg: f64 (-90 to 90)
└─ fov_deg, near, far: f64

PickingSystem
├─ render_pick_pass() → off-screen color IDs
└─ query(x, y) → PickResult

FieldProbe
├─ position_mm: [f64; 3]
├─ element_tag: u64
├─ barycentric: [f64; 4] (interpolation weights)
└─ field: ProbeFieldValue (E, H vectors)
```

### GPU Vertex Layouts
```
FieldVertex (per-node)
├─ position: [f32; 3]
├─ normal: [f32; 3]
└─ field_value: f32 (normalized [0, 1])

ArrowInstance (instanced rendering)
├─ position: [f32; 3]
├─ direction: [f32; 3]
└─ magnitude: f32
```

---

## File Format Quick Lookup

| File | Format | Q3D Only? | Sections | Key Fields |
|------|--------|-----------|----------|-----------|
| `rlcg_matrix.json` | JSON | YES | R, L, C, G matrices per frequency | num_terminals, terminal_names, matrices |
| `equivalent_circuit.sp` | SPICE netlist | YES | .SUBCKT with R/L/K/C elements | Self/mutual impedance, capacitance |
| `*.emsfld` | Binary | NO | Header + Index + Field blocks | Header: magic, version, field_type, data_type |
| `final_mesh.msh` | Gmsh 4.1 | NO | $PhysicalNames, $Entities, $Nodes, $Elements | Version 4.1, binary mode recommended |

### EmsFldHeader (128 bytes)
```rust
magic: "EMSFLD\0\0"          // [u8; 8]
version: u32 = 1
byte_order: u32 = 0x01020304 (little-endian indicator)
field_type: u32             // 0=E, 1=H, 2=J, 3=Combined
data_type: u32              // 0=complex f64, 1=complex f32
num_nodes: u64
num_components: u32         // 3 (vector) or 1 (scalar)
num_frequencies: u32
frequency_unit: u32         // 0=Hz, 1=kHz, 2=MHz, 3=GHz
freq_table_offset: u64
index_offset: u64
data_offset: u64
mesh_file: [u8; 32]         // associated .msh filename
_reserved: [u8; 12]
```

### FieldBlockInfo (16 bytes × num_frequencies)
```rust
offset: u64          // file offset of this frequency block
size_bytes: u64      // block size in bytes
```

---

## M7 IMPLEMENTATION TASKS

### 1. Isosurface Extraction
**Function**: `IsosurfaceExtractor::extract(mesh, field_values, iso_value) → TriangleMesh`
**Algorithm**: Marching Tetrahedra
- Input: tetrahedral mesh + per-node scalar field
- Output: triangle mesh (vertices, normals, indices)
- For each tet: classify 4 vertices, find edge crossings, generate 0-2 triangles

### 2. Slice Mesh Generation
**Function**: `SliceMeshGenerator::generate(mesh, field_values, plane) → SliceMesh`
**Algorithm**: Tet-plane intersection
- For each tet: find face-plane intersections
- Triangulate resulting polygon
- Interpolate field values using barycentric coordinates
- Compute normals

### 3. GPU Pipeline Setup
- **Vertex Shader**: Transform position (MVP), pass field_value to fragment
- **Fragment Shader**: Sample 1D colormap texture, apply Lambertian lighting
- **Uniforms**: mvp, model, light_dir, ambient, value_min, value_max, opacity
- **Bindings**: [0] Uniforms, [1] colormap texture 1D, [2] sampler

### 4. Orbit Camera
- **view_matrix()**: Spherical to Cartesian, compute view matrix
- **projection_matrix()**: Perspective projection with aspect ratio
- **orbit(dx, dy)**: Update azimuth/elevation
- **zoom(delta)**: Update distance
- **fit_to_bounds()**: Auto-fit camera to AABB

### 5. GPU Picking
- Off-screen render pass with object IDs as colors
- Readback pixel at mouse position
- Map color → object_id

### 6. Field Probe
- Find containing tetrahedron (point-in-tet test)
- Compute barycentric coordinates [w0, w1, w2, w3]
- Interpolate field: `E = w0*E0 + w1*E1 + w2*E2 + w3*E3`
- Compute magnitude: `|E| = sqrt(|Ex|² + |Ey|² + |Ez|²)`

### 7. PNG Screenshot
1. Create off-screen texture (RGBA, dimensions)
2. Render to texture
3. Readback texture data (GPU → CPU buffer)
4. Encode PNG using `image` crate

---

## M8 IMPLEMENTATION TASKS

### 1. RLCG Matrix Parser
- Load `rlcg_matrix.json`
- Structure: matrices { R, L, C, G } with data_per_frequency
- Matrices are symmetric
- Access: `matrices["R"][freq_idx].matrix[i][j]`

### 2. SPICE Exporter
- Generate `.sp` file from RLCG data
- Model types: BroadbandLumped, FrequencyDependentLumped, TLineModel, SParameterBlock
- Elements: R (resistance), L (inductance), K (mutual coupling), C (capacitance)
- Structure: `.SUBCKT Q3D_Model <terminals> ... .ENDS`

### 3. .emsfld Header Parser
```rust
// Read 128 bytes from file start
let header: EmsFldHeader = deserialize_from_bytes(&buf[0..128]);
// Validate magic == b"EMSFLD\0\0"
// Validate version == 1
// Check byte_order for endianness
```

### 4. FieldFileHandle::open()
1. Open file, create mmap
2. Read header (128 bytes)
3. Read frequency table at freq_table_offset
4. Read index table (num_frequencies × FieldBlockInfo) at index_offset
5. Store in FieldFileHandle struct

### 5. FieldFileHandle::slice()
```rust
fn slice(&self, freq_idx: usize) -> Result<FieldSlice<'_>, Error> {
    let block_info = &self.block_index[freq_idx];
    let offset = block_info.offset as usize;
    let size = block_info.size_bytes as usize;
    let raw_bytes = &self.mmap[offset..offset + size];
    
    // Cast to &[f64] or &[f32]
    // Return FieldSlice with zero-copy reference to mmap region
}
```

### 6. VTK Exporter
- Input: mesh (.msh) + field_values (.emsfld)
- Output: VTK Unstructured Grid file
- Cell type: 10 (tetrahedron)
- Point data: field value arrays
- Format: legacy ASCII or XML (.vtu)

### 7. Phase Animator
```rust
fn evaluate(&self, re: &[f64], im: &[f64]) -> Vec<f64> {
    let phase_rad = current_phase_deg.to_radians();
    re.iter().zip(im.iter())
        .map(|(r, i)| r * cos(phase_rad) - i * sin(phase_rad))
        .collect()
}
```

---

## Storage Requirements

### Field Data (.emsfld)
- **10K nodes, 1 frequency, complex f64**: 480 KB
- **10K nodes, 301 frequencies, complex f64**: 141 MB
- **10K nodes, 301 frequencies, complex f32**: 70 MB

### Mesh Storage (10K tets)
- **ASCII MSH**: ~1.2 MB
- **Binary MSH**: ~520 KB

### RLCG Matrix (4 terminals, 50 frequencies)
- **JSON**: ~50-100 KB

---

## Key Formulas

### Barycentric Interpolation
```
weights [w0, w1, w2, w3] sum to 1.0
value_at_point = w0*v0 + w1*v1 + w2*v2 + w3*v3
```

### Phase Animation (Complex → Real)
```
E_real = Re(E) * cos(phase) - Im(E) * sin(phase)
where phase in degrees → radians
```

### Colormap Normalization
```
normalized = clamp((value - min) / (max - min), 0.0, 1.0)
color = sample_colormap_texture(normalized)
```

### Orbit Camera Position
```
x = distance * sin(azimuth_rad) * cos(elevation_rad)
y = distance * sin(elevation_rad)
z = distance * cos(azimuth_rad) * cos(elevation_rad)
position = target + (x, y, z)
```

---

## Error Handling

### File Loading
- Check magic number == `b"EMSFLD\0\0"`
- Check byte_order for endianness compatibility
- Validate index bounds (freq_idx < num_frequencies)
- Handle mmap failures gracefully

### Geometry Operations
- Point-in-tet test with barycentric coordinate bounds checking
- Plane-tet intersection edge case handling
- Division by zero in barycentric coordinate computation

### GPU Operations
- Check shader compilation errors
- Verify texture format support
- Handle off-screen render target creation failures

---

## Testing Priorities

1. **IsosurfaceExtractor** - Unit test with simple tet mesh, verify triangle count
2. **SliceMeshGenerator** - Test aligned planes (XY/XZ/YZ) and arbitrary planes
3. **OrbitCamera** - Verify view/projection matrices, test zoom bounds
4. **FieldFileHandle** - Test random access to different frequency indices
5. **FieldProbe** - Verify barycentric interpolation accuracy
6. **PNG export** - Compare screenshots with reference images

