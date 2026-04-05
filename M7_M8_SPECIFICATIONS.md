# EMStudio M7 & M8 Implementation Specifications
## Extracted from Design Documents

**Extraction Date**: 2026-04-05  
**Source Documents**:
- `/Users/alex/works/emstudio/docs/em-result-visualization-design.md`
- `/Users/alex/works/emstudio/docs/em-result-file-formats.md`

---

## I. 3D FIELD VISUALIZATION (§4)

### 4.1 Field Overlay Types

**Enum: FieldOverlayType**

```rust
pub enum FieldOverlayType {
    /// Surface cloud plot: color-mapped field magnitude on geometry surfaces
    SurfacePlot {
        surfaces: Option<Vec<String>>,  // None = all outer surfaces
    },
    /// Slice plot: volume field data cross-section through arbitrary plane
    SlicePlot {
        plane: SlicePlane,
    },
    /// Vector arrow plot: field vector arrows at sampled points
    VectorPlot {
        spacing_mm: f64,      // arrow density (sampling interval)
        arrow_scale: f64,     // arrow length scaling factor
    },
    /// Isosurface: 3D contour surface where field magnitude equals iso_value
    IsosurfacePlot {
        iso_value: f64,
    },
    /// Animation: phase sweep from 0° to 360° showing time-domain evolution
    AnimatedPlot {
        phase_step_deg: f64,
        fps: f64,
    },
}
```

**Enum: SlicePlane**

```rust
pub enum SlicePlane {
    XY { z_mm: f64 },                          // Z = constant
    XZ { y_mm: f64 },                          // Y = constant
    YZ { x_mm: f64 },                          // X = constant
    Arbitrary {
        normal: [f64; 3],                      // plane normal vector
        point: [f64; 3],                       // point on plane
    },
}
```

**Enum: FieldType**

```rust
pub enum FieldType {
    Electric,              // E-field (V/m) — HFSS & Q3D
    Magnetic,              // H-field (A/m) — HFSS only
    Current,               // J-field (A/m²) — surface current
    Poynting,              // S-field (W/m²) — power density
    VolumeCurrent,         // Jvol-field (A/m²) — volume current (Q3D)
    ChargeDistribution,    // ρ-field (C/m²) — charge density (Q3D)
    OhmicLoss,             // P_loss (W/m³) — ohmic loss density (Q3D)
}

pub enum FieldComponent {
    Magnitude,             // |F| = sqrt(|Fx|² + |Fy|² + |Fz|²)
    X, Y, Z,               // individual components
    Vector,                // full vector (for arrow plots)
}
```

---

### 4.2 ISOSURFACE / MARCHING TETRAHEDRA

**Algorithm: Marching Tetrahedra**

```rust
pub struct IsosurfaceExtractor;

impl IsosurfaceExtractor {
    /// Extract isosurface from tetrahedral mesh + scalar field values
    /// Returns triangle mesh representation of iso_value contour
    pub fn extract(
        mesh: &MshMesh,
        field_values: &[f64],  // scalar value per node
        iso_value: f64,
    ) -> TriangleMesh {
        // For each tetrahedron:
        //   1. Evaluate field values at 4 vertices
        //   2. Determine intersections between field surface and edges
        //   3. Generate 0-2 triangles per tet based on edge crossing patterns
        //   4. Accumulate into triangle mesh
        todo!()
    }
}

pub struct TriangleMesh {
    pub vertices: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    pub indices: Vec<u32>,
}
```

---

### 4.3 CUT PLANE VISUALIZATION

**Algorithm: Slice Mesh Generation**

```rust
pub struct SliceMeshGenerator;

impl SliceMeshGenerator {
    /// Generate triangle mesh from plane-tet intersections
    /// Each slice vertex interpolates field value via barycentric coords
    pub fn generate(
        mesh: &MshMesh,
        field_values: &[f64],
        plane: &SlicePlane,
    ) -> SliceMesh {
        // For each tet:
        //   1. Classify vertices (above/below/on plane)
        //   2. Compute edge-plane intersection points
        //   3. Generate polygon from intersections
        //   4. Triangulate polygon
        //   5. Interpolate field values at intersection points
        todo!()
    }
}

pub struct SliceMesh {
    pub vertices: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    pub field_values: Vec<f32>,
    pub indices: Vec<u32>,
}
```

---

## II. INTERACTIVE SYSTEMS (§7)

### 7.1 ORBIT CAMERA

```rust
pub struct OrbitCamera {
    pub target: [f64; 3],              // view target point
    pub distance: f64,                 // distance from target
    pub azimuth_deg: f64,              // rotation around Y-axis
    pub elevation_deg: f64,            // tilt (-90 to 90)
    pub fov_deg: f64,                  // vertical field of view
    pub near: f64,                     // near clipping plane
    pub far: f64,                      // far clipping plane
}

impl OrbitCamera {
    pub fn view_matrix(&self) -> [[f64; 4]; 4] { todo!() }
    pub fn projection_matrix(&self, aspect_ratio: f64) -> [[f64; 4]; 4] { todo!() }
    
    pub fn orbit(&mut self, delta_x: f64, delta_y: f64) {
        self.azimuth_deg += delta_x * 0.5;
        self.elevation_deg = (self.elevation_deg + delta_y * 0.5)
            .clamp(-89.0, 89.0);
    }
    
    pub fn pan(&mut self, delta_x: f64, delta_y: f64) { todo!() }
    
    pub fn zoom(&mut self, delta: f64) {
        self.distance *= (1.0 - delta * 0.1).max(0.01);
    }
}
```

---

### 7.2 GPU PICKING + PROBE

**Field Probe - Spatial Interpolation**

```rust
pub struct FieldProbe {
    pub position_mm: [f64; 3],     // world coordinates
    pub result: Option<ProbeResult>,
}

pub struct ProbeResult {
    pub element_tag: u64,          // containing tet
    pub barycentric: [f64; 4],     // interpolation weights
    pub field: ProbeFieldValue,
}

pub struct ProbeFieldValue {
    pub e_field: [Complex64; 3],   // Ex, Ey, Ez (complex)
    pub e_magnitude: f64,          // |E|
    pub h_field: [Complex64; 3],   // Hx, Hy, Hz (complex)
    pub h_magnitude: f64,          // |H|
}
```

---

## III. EXPORT FUNCTIONALITY (§10)

### 10.1 PNG Export

```rust
pub fn capture_screenshot(
    renderer: &SceneRenderer,
    width: u32,
    height: u32,
) -> Result<Vec<u8>, RenderError> {
    // 1. Create off-screen texture
    // 2. Render to texture
    // 3. Readback to CPU
    // 4. Encode as PNG
    todo!()
}
```

### 10.2 VTK Export

```
mesh (.msh) → VTK Unstructured Grid (tetrahedral cells)
field_values (.emsfld) → VTK Point Data Arrays
```

---

## IV. FILE FORMAT SPECIFICATIONS

### 4.1 RLCG MATRIX FORMAT (§2.9)

**File**: `rlcg_matrix.json` (Q3D-specific)

**Key Structure**:

```json
{
  "format_version": "1.0",
  "file_type": "RLCGMatrixData",
  "num_terminals": 4,
  "terminal_names": ["Signal1:T1", "Signal1:T2", "Signal2:T3", "Signal2:T4"],
  "num_frequencies": 50,
  "frequencies": [0.01, 0.02, ..., 5.0],
  "matrices": {
    "R": {
      "unit": "ohm",
      "data_per_frequency": [
        {
          "frequency": 0.01,
          "matrix": [[...], [...], ...]
        },
        ...
      ]
    },
    "L": { "unit": "nH", "data_per_frequency": [...] },
    "C": { "unit": "pF", "data_per_frequency": [...] },
    "G": { "unit": "mS", "data_per_frequency": [...] }
  },
  "dc_data": {
    "R_dc": { "unit": "ohm", "matrix": [[...], ...] },
    "L_dc": { "unit": "nH", "matrix": [[...], ...] }
  }
}
```

**Matrix Properties**:
- All matrices are symmetric: `M[i][j] == M[j][i]`
- R, L: frequency-dependent (skin effect)
- C: frequency-independent
- C diagonal: positive (total capacitance)
- C off-diagonal: negative (coupling)

---

### 4.2 EQUIVALENT CIRCUIT SPICE FORMAT (§2.10)

**File**: `equivalent_circuit.sp` (Q3D export)

```spice
.SUBCKT Q3D_Model Signal1_src Signal1_sink Signal2_src Signal2_sink GND

* Self impedance
R_self_1  Signal1_src  n1_1  0.285
L_self_1  n1_1         Signal1_sink  2.85n

* Mutual inductance
K_12  L_self_1  L_self_2  0.147

* Capacitance
C_1g   Signal1_src  GND  0.328p
C_12   Signal1_src  Signal2_src  0.085p

.ENDS Q3D_Model
```

**Model Types**:
- `BroadbandLumped`: Simple RLC (general SI)
- `FrequencyDependentLumped`: Multi-stage RL ladder (skin effect)
- `TLineModel`: Transmission line W-element
- `SParameterBlock`: S-parameter behavioral (Touchstone)

---

### 4.3 FIELD DATA BINARY FORMAT (§3.3)

**File**: `.emsfld` (frequency-domain FEM field)

**Layout**:

```
Header (128 bytes)
├─ magic: "EMSFLD\0\0"
├─ version: u32
├─ byte_order: u32 (0x01020304 for endian detection)
├─ field_type: u32 (0=E, 1=H, 2=J, 3=Combined)
├─ data_type: u32 (0=complex f64, 1=complex f32)
├─ num_nodes: u64
├─ num_components: u32 (3 or 1)
├─ num_frequencies: u32
├─ frequency_unit: u32 (0=Hz, 1=kHz, 2=MHz, 3=GHz)
├─ freq_table_offset: u64
├─ index_offset: u64
├─ data_offset: u64
└─ mesh_file: [u8; 32]

Frequency Table (num_frequencies × f64)
├─ freq[0]: f64
├─ freq[1]: f64
└─ ...

Field Data Index (num_frequencies × FieldBlockInfo)
├─ [0]: { offset: u64, size_bytes: u64 }
├─ [1]: { offset: u64, size_bytes: u64 }
└─ ...

Field Blocks (frequency-indexed)
├─ Block 0: num_nodes × (6 × f64 or 6 × f32)
├─ Block 1: num_nodes × (6 × f64 or 6 × f32)
└─ ...
```

**Per-Node Format (complex f64 vector)**:
```
[re_x: f64, im_x: f64, re_y: f64, im_y: f64, re_z: f64, im_z: f64]
= 48 bytes/node
```

**Random Access**:
```
1. Read header from offset 0
2. Seek to index_offset + freq_idx × 16
3. Read FieldBlockInfo → offset, size_bytes
4. Seek to offset, read size_bytes
5. Cast to &[f64] or &[f32] (zero-copy via mmap)
```

**Storage (10K nodes)**:
- 1 frequency complex f64: 480 KB
- 301 frequencies complex f64: 141 MB
- 1 frequency complex f32: 240 KB
- 301 frequencies complex f32: 70 MB

---

## V. GPU RENDERING PIPELINE DETAILS

**Vertex Buffer Layout**:

```rust
#[repr(C)]
pub struct FieldVertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
    pub field_value: f32,    // normalized [0, 1]
}
```

**Vector Arrow Instance**:

```rust
#[repr(C)]
pub struct ArrowInstance {
    pub position: [f32; 3],  // base
    pub direction: [f32; 3], // normalized
    pub magnitude: f32,      // length factor
    pub _pad: f32,
}
```

**WGSL Uniforms**:

```wgsl
struct Uniforms {
    mvp: mat4x4<f32>,
    model: mat4x4<f32>,
    light_dir: vec3<f32>,
    ambient: f32,
    value_min: f32,
    value_max: f32,
    opacity: f32,
    _pad: f32,
};

// Bindings:
// [0] Uniforms buffer
// [1] 1D colormap texture
// [2] Colormap sampler
```

---

## VI. IMPLEMENTATION CHECKLIST

### M7 Tasks

- [ ] IsosurfaceExtractor::extract() - Marching Tetrahedra
- [ ] SliceMeshGenerator::generate() - Tet-plane intersections
- [ ] FieldVertex layout - GPU vertex buffer
- [ ] WGSL Shaders - Colormap + lighting
- [ ] ArrowInstance rendering - Vector visualization
- [ ] OrbitCamera matrices - View/projection
- [ ] PickingSystem - GPU color-picking
- [ ] FieldProbe - Barycentric interpolation
- [ ] PNG export - Screenshot capture

### M8 Tasks

- [ ] RLCG JSON parser - Load matrix data
- [ ] SPICE exporter - Generate netlist
- [ ] EmsFldHeader parser - Binary format
- [ ] FieldFileHandle::open() - mmap + indexing
- [ ] FieldFileHandle::slice() - Zero-copy access
- [ ] VTK exporter - Mesh + fields
- [ ] Phase animator - Complex field evolution

---

## KEY ALGORITHMS & FORMULAS

**Barycentric Interpolation (Tetrahedron)**:
```
point P = w0*v0 + w1*v1 + w2*v2 + w3*v3
where w0 + w1 + w2 + w3 = 1.0

field(P) = w0*f0 + w1*f1 + w2*f2 + w3*f3
```

**Phase Animation**:
```
E(t) = Re[ E_complex * e^(j*phase_rad) ]
     = Re_part * cos(phase) - Im_part * sin(phase)
```

**Colormap Normalization**:
```
normalized = clamp((value - min) / (max - min), 0.0, 1.0)
```

**Orbit Camera Position**:
```
position = target + distance * [
  sin(azimuth_rad) * cos(elevation_rad),
  sin(elevation_rad),
  cos(azimuth_rad) * cos(elevation_rad)
]
```

