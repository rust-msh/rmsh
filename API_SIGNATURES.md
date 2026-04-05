# M7 & M8 API Signatures & Type Definitions

## M7: Visualization & Interaction APIs

### Isosurface Extraction

```rust
pub struct TriangleMesh {
    pub vertices: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    pub indices: Vec<u32>,
}

pub struct IsosurfaceExtractor;

impl IsosurfaceExtractor {
    /// Extract isosurface from tet mesh + scalar field
    /// 
    /// # Arguments
    /// * `mesh` - Tetrahedral FEM mesh
    /// * `field_values` - Per-node scalar field values
    /// * `iso_value` - Target isosurface value
    ///
    /// # Returns
    /// Triangle mesh representation of the isosurface
    pub fn extract(
        mesh: &MshMesh,
        field_values: &[f64],
        iso_value: f64,
    ) -> Result<TriangleMesh, GeometryError>;
    
    /// Compute vertex normals from triangle connectivity
    pub fn compute_vertex_normals(
        vertices: &[[f32; 3]],
        indices: &[u32],
    ) -> Vec<[f32; 3]>;
}
```

### Slice Mesh Generation

```rust
pub struct SliceMesh {
    pub vertices: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    pub field_values: Vec<f32>,  // interpolated at vertices
    pub indices: Vec<u32>,
}

pub struct SliceMeshGenerator;

impl SliceMeshGenerator {
    /// Generate triangle mesh from plane-tet intersections
    ///
    /// # Arguments
    /// * `mesh` - Tetrahedral FEM mesh
    /// * `field_values` - Per-node scalar field values
    /// * `plane` - Slicing plane definition
    ///
    /// # Returns
    /// Triangle mesh of the cross-section with interpolated field values
    pub fn generate(
        mesh: &MshMesh,
        field_values: &[f64],
        plane: &SlicePlane,
    ) -> Result<SliceMesh, GeometryError>;
    
    /// Classify vertex position relative to plane
    fn point_plane_side(point: [f64; 3], plane: &SlicePlane) -> f64;
    
    /// Compute edge-plane intersection with interpolation
    fn edge_plane_intersection(
        v1: [f64; 3], f1: f64,
        v2: [f64; 3], f2: f64,
        plane: &SlicePlane,
    ) -> Option<([f32; 3], f32)>;
}
```

### GPU Rendering Pipeline

```rust
#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
pub struct FieldVertex {
    pub position: [f32; 3],      // vertex position
    pub normal: [f32; 3],        // surface normal for lighting
    pub field_value: f32,        // normalized field [0, 1]
}

pub struct FieldPipeline {
    layout: wgpu::PipelineLayout,
    pipeline: wgpu::RenderPipeline,
    colormap_texture: wgpu::Texture,
    colormap_sampler: wgpu::Sampler,
    uniforms_buffer: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
}

impl FieldPipeline {
    /// Create field visualization pipeline
    pub fn new(
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        format: wgpu::TextureFormat,
        colormap: &ColormapType,
    ) -> Result<Self, RenderError>;
    
    /// Update uniforms (MVP, lighting, value range)
    pub fn update_uniforms(
        &mut self,
        queue: &wgpu::Queue,
        mvp: [[f32; 4]; 4],
        model: [[f32; 4]; 4],
        light_dir: [f32; 3],
        ambient: f32,
        value_min: f32,
        value_max: f32,
        opacity: f32,
    );
    
    /// Create vertex buffer from mesh
    pub fn create_vertex_buffer(
        device: &wgpu::Device,
        vertices: &[FieldVertex],
    ) -> wgpu::Buffer;
    
    /// Create index buffer
    pub fn create_index_buffer(
        device: &wgpu::Device,
        indices: &[u32],
    ) -> wgpu::Buffer;
    
    /// Render call (to be invoked from render pass)
    pub fn render<'a>(
        &'a self,
        render_pass: &mut wgpu::RenderPass<'a>,
        vertex_buffer: &'a wgpu::Buffer,
        index_buffer: &'a wgpu::Buffer,
        index_count: u32,
    );
}
```

**WGSL Shaders (embedded in code)**:

```wgsl
// vertex.wgsl
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

@group(0) @binding(0) var<uniform> u: Uniforms;
@group(0) @binding(1) var colormap_texture: texture_1d<f32>;
@group(0) @binding(2) var colormap_sampler: sampler;

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) field_value: f32,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) world_normal: vec3<f32>,
    @location(1) normalized_value: f32,
};

@vertex
fn vs_main(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    out.clip_position = u.mvp * vec4<f32>(in.position, 1.0);
    out.world_normal = (u.model * vec4<f32>(in.normal, 0.0)).xyz;
    out.normalized_value = clamp(
        (in.field_value - u.value_min) / (u.value_max - u.value_min),
        0.0, 1.0
    );
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    let color = textureSample(colormap_texture, colormap_sampler, in.normalized_value);
    let n = normalize(in.world_normal);
    let diffuse = max(dot(n, normalize(u.light_dir)), 0.0);
    let lighting = u.ambient + (1.0 - u.ambient) * diffuse;
    return vec4<f32>(color.rgb * lighting, color.a * u.opacity);
}
```

### Orbit Camera

```rust
#[derive(Clone)]
pub struct OrbitCamera {
    pub target: [f64; 3],              // view center
    pub distance: f64,                 // radius
    pub azimuth_deg: f64,              // horizontal rotation (Y-axis)
    pub elevation_deg: f64,            // vertical tilt (-90 to 90)
    pub fov_deg: f64,                  // vertical field of view
    pub near: f64,                     // near clipping plane
    pub far: f64,                      // far clipping plane
}

impl OrbitCamera {
    /// Create camera at default position
    pub fn new() -> Self {
        Self {
            target: [0.0, 0.0, 0.0],
            distance: 100.0,
            azimuth_deg: 45.0,
            elevation_deg: 35.0,
            fov_deg: 45.0,
            near: 0.1,
            far: 10000.0,
        }
    }
    
    /// Compute world position from spherical coordinates
    pub fn position(&self) -> [f64; 3] {
        let az_rad = self.azimuth_deg.to_radians();
        let el_rad = self.elevation_deg.to_radians();
        let cos_el = el_rad.cos();
        [
            self.target[0] + self.distance * az_rad.sin() * cos_el,
            self.target[1] + self.distance * el_rad.sin(),
            self.target[2] + self.distance * az_rad.cos() * cos_el,
        ]
    }
    
    /// Compute view matrix (position, target, up=[0,1,0])
    pub fn view_matrix(&self) -> [[f64; 4]; 4];
    
    /// Compute projection matrix
    pub fn projection_matrix(&self, aspect_ratio: f64) -> [[f64; 4]; 4];
    
    /// Mouse drag rotation
    pub fn orbit(&mut self, delta_x: f64, delta_y: f64) {
        self.azimuth_deg += delta_x * 0.5;
        self.elevation_deg = (self.elevation_deg + delta_y * 0.5).clamp(-89.0, 89.0);
    }
    
    /// Middle-mouse pan
    pub fn pan(&mut self, delta_x: f64, delta_y: f64);
    
    /// Scroll zoom
    pub fn zoom(&mut self, delta: f64) {
        self.distance *= (1.0 - delta * 0.1).max(0.01);
    }
    
    /// Auto-fit to bounding box
    pub fn fit_to_bounds(&mut self, aabb_min: [f64; 3], aabb_max: [f64; 3]) {
        let center = [
            (aabb_min[0] + aabb_max[0]) / 2.0,
            (aabb_min[1] + aabb_max[1]) / 2.0,
            (aabb_min[2] + aabb_max[2]) / 2.0,
        ];
        let extent = [
            (aabb_max[0] - aabb_min[0]) / 2.0,
            (aabb_max[1] - aabb_min[1]) / 2.0,
            (aabb_max[2] - aabb_min[2]) / 2.0,
        ];
        let radius = (extent[0] * extent[0] + extent[1] * extent[1] + extent[2] * extent[2]).sqrt();
        self.target = center;
        self.distance = radius / (self.fov_deg.to_radians() / 2.0).tan();
    }
    
    /// Apply preset viewpoint
    pub fn set_view_preset(&mut self, preset: ViewPreset) {
        match preset {
            ViewPreset::Front  => { self.azimuth_deg = 0.0;   self.elevation_deg = 0.0; }
            ViewPreset::Back   => { self.azimuth_deg = 180.0; self.elevation_deg = 0.0; }
            ViewPreset::Left   => { self.azimuth_deg = -90.0; self.elevation_deg = 0.0; }
            ViewPreset::Right  => { self.azimuth_deg = 90.0;  self.elevation_deg = 0.0; }
            ViewPreset::Top    => { self.azimuth_deg = 0.0;   self.elevation_deg = 89.0; }
            ViewPreset::Bottom => { self.azimuth_deg = 0.0;   self.elevation_deg = -89.0; }
            ViewPreset::Iso    => { self.azimuth_deg = 45.0;  self.elevation_deg = 35.0; }
        }
    }
}

pub enum ViewPreset {
    Front, Back, Left, Right, Top, Bottom, Iso,
}
```

### GPU Picking System

```rust
pub struct PickingSystem {
    pick_texture: wgpu::Texture,
    pick_view: wgpu::TextureView,
    readback_buffer: wgpu::Buffer,
    width: u32,
    height: u32,
}

pub enum PickResult {
    GeometryObject { object_id: u64 },
    MeshElement { element_tag: u64 },
    FieldValue { node_tag: u64, value: f64 },
}

impl PickingSystem {
    /// Create picking system with target dimensions
    pub fn new(
        device: &wgpu::Device,
        width: u32,
        height: u32,
    ) -> Result<Self, RenderError>;
    
    /// Render picking pass (object IDs as colors)
    pub fn render_pick_pass(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        device: &wgpu::Device,
        scene: &SceneData,
        camera: &OrbitCamera,
    ) -> Result<(), RenderError>;
    
    /// Query object at pixel location
    pub fn query(
        &self,
        x: u32,
        y: u32,
    ) -> Result<Option<PickResult>, RenderError>;
}
```

### Field Probe

```rust
use num_complex::Complex64;

pub struct FieldProbe {
    pub position_mm: [f64; 3],
    pub result: Option<ProbeResult>,
}

pub struct ProbeResult {
    pub element_tag: u64,
    pub barycentric: [f64; 4],
    pub field: ProbeFieldValue,
}

pub struct ProbeFieldValue {
    pub e_field: [Complex64; 3],
    pub e_magnitude: f64,
    pub h_field: [Complex64; 3],
    pub h_magnitude: f64,
}

impl FieldProbe {
    /// Query field at world position
    pub fn query(
        position_mm: [f64; 3],
        mesh: &MshMesh,
        e_field_values: &[[Complex64; 3]],  // per-node
        h_field_values: &[[Complex64; 3]],  // per-node
    ) -> Result<Self, ProbeError>;
    
    /// Find containing tetrahedron
    fn find_containing_tet(
        point: [f64; 3],
        mesh: &MshMesh,
    ) -> Option<(u64, [f64; 4])>;  // (tet_tag, barycentric)
    
    /// Compute barycentric coordinates
    fn barycentric_coordinates(
        point: [f64; 3],
        v0: [f64; 3],
        v1: [f64; 3],
        v2: [f64; 3],
        v3: [f64; 3],
    ) -> [f64; 4];
}
```

### Screenshot Export

```rust
pub fn capture_screenshot(
    renderer: &SceneRenderer,
    camera: &OrbitCamera,
    width: u32,
    height: u32,
) -> Result<Vec<u8>, ExportError>;

pub fn capture_screenshot_to_file(
    renderer: &SceneRenderer,
    camera: &OrbitCamera,
    width: u32,
    height: u32,
    path: &Path,
) -> Result<(), ExportError>;
```

---

## M8: Data Access & Export APIs

### Field File Format

```rust
#[repr(C, packed)]
pub struct EmsFldHeader {
    pub magic: [u8; 8],              // b"EMSFLD\0\0"
    pub version: u32,                // 1
    pub byte_order: u32,             // 0x01020304
    pub field_type: u32,             // 0=E, 1=H, 2=J, 3=Combined
    pub data_type: u32,              // 0=complex f64, 1=complex f32
    pub num_nodes: u64,
    pub num_components: u32,         // 3 (vector) or 1 (scalar)
    pub num_frequencies: u32,
    pub frequency_unit: u32,         // 0=Hz, 1=kHz, 2=MHz, 3=GHz
    pub freq_table_offset: u64,
    pub index_offset: u64,
    pub data_offset: u64,
    pub mesh_file: [u8; 32],
    pub _reserved: [u8; 12],
}

#[repr(C, packed)]
pub struct FieldBlockInfo {
    pub offset: u64,
    pub size_bytes: u64,
}

pub struct FieldSlice<'a> {
    pub frequency_hz: f64,
    pub num_nodes: usize,
    pub num_components: usize,
    pub data_f64: Option<&'a [f64]>,
    pub data_f32: Option<&'a [f32]>,
}
```

### Field File Access

```rust
pub struct FieldFileHandle {
    mmap: memmap2::Mmap,
    header: EmsFldHeader,
    frequencies: Vec<f64>,
    block_index: Vec<FieldBlockInfo>,
}

impl FieldFileHandle {
    /// Open and memory-map a field file
    pub fn open(path: &Path) -> Result<Self, DataError> {
        let file = std::fs::File::open(path)?;
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        
        // Parse header
        let header_bytes = &mmap[0..128];
        let header = Self::parse_header(header_bytes)?;
        
        // Validate magic and byte order
        if header.magic != b"EMSFLD\0\0" {
            return Err(DataError::InvalidMagic);
        }
        if header.byte_order != 0x01020304 {
            return Err(DataError::ByteOrderMismatch);
        }
        
        // Parse frequency table
        let freq_offset = header.freq_table_offset as usize;
        let freq_count = header.num_frequencies as usize;
        let frequencies = Self::parse_frequencies(&mmap, freq_offset, freq_count)?;
        
        // Parse block index
        let index_offset = header.index_offset as usize;
        let block_index = Self::parse_block_index(&mmap, index_offset, freq_count)?;
        
        Ok(Self {
            mmap,
            header,
            frequencies,
            block_index,
        })
    }
    
    /// Get zero-copy field slice for frequency index
    pub fn slice(&self, freq_idx: usize) -> Result<FieldSlice<'_>, DataError> {
        if freq_idx >= self.block_index.len() {
            return Err(DataError::FrequencyIndexOutOfBounds);
        }
        
        let block_info = &self.block_index[freq_idx];
        let offset = block_info.offset as usize;
        let size = block_info.size_bytes as usize;
        let raw_bytes = &self.mmap[offset..offset + size];
        
        let data_f64 = if self.header.data_type == 0 {
            Some(Self::cast_to_f64_slice(raw_bytes))
        } else {
            None
        };
        
        Ok(FieldSlice {
            frequency_hz: self.frequencies[freq_idx],
            num_nodes: self.header.num_nodes as usize,
            num_components: self.header.num_components as usize,
            data_f64,
            data_f32: None,
        })
    }
    
    pub fn frequency_list(&self) -> &[f64] {
        &self.frequencies
    }
}
```

### RLCG Matrix Handling

```rust
use serde_json::Value;

pub struct RlcgMatrixData {
    pub num_terminals: usize,
    pub terminal_names: Vec<String>,
    pub frequencies: Vec<f64>,
    pub r_matrix: Vec<Vec<Vec<f64>>>,  // [freq_idx][i][j]
    pub l_matrix: Vec<Vec<Vec<f64>>>,
    pub c_matrix: Vec<Vec<Vec<f64>>>,
    pub g_matrix: Vec<Vec<Vec<f64>>>,
    pub dc_r: Vec<Vec<f64>>,
    pub dc_l: Vec<Vec<f64>>,
}

impl RlcgMatrixData {
    /// Parse rlcg_matrix.json
    pub fn from_json(data: &Value) -> Result<Self, ParseError> {
        // Extract header fields
        let num_terminals = data["num_terminals"].as_u64().ok_or(ParseError::Missing)? as usize;
        let terminal_names = data["terminal_names"]
            .as_array()
            .ok_or(ParseError::Missing)?
            .iter()
            .filter_map(|v| v.as_str().map(|s| s.to_string()))
            .collect();
        
        let frequencies = data["frequencies"]
            .as_array()
            .ok_or(ParseError::Missing)?
            .iter()
            .filter_map(|v| v.as_f64())
            .collect();
        
        // Extract matrix data
        let r_matrix = Self::extract_matrix_per_frequency(
            &data["matrices"]["R"]
        )?;
        let l_matrix = Self::extract_matrix_per_frequency(
            &data["matrices"]["L"]
        )?;
        let c_matrix = Self::extract_matrix_per_frequency(
            &data["matrices"]["C"]
        )?;
        let g_matrix = Self::extract_matrix_per_frequency(
            &data["matrices"]["G"]
        )?;
        
        // Extract DC data
        let dc_r = Self::extract_single_matrix(&data["dc_data"]["R_dc"])?;
        let dc_l = Self::extract_single_matrix(&data["dc_data"]["L_dc"])?;
        
        Ok(Self {
            num_terminals,
            terminal_names,
            frequencies,
            r_matrix,
            l_matrix,
            c_matrix,
            g_matrix,
            dc_r,
            dc_l,
        })
    }
    
    pub fn get_r_matrix(&self, freq_idx: usize) -> Option<&Vec<Vec<f64>>> {
        self.r_matrix.get(freq_idx)
    }
    
    pub fn get_impedance_at_frequency(
        &self,
        freq_hz: f64,
        port_i: usize,
        port_j: usize,
    ) -> Complex64 {
        // Interpolate R and L to frequency
        // Z = R + j*ω*L where ω = 2π*freq_hz
        todo!()
    }
}
```

### SPICE Exporter

```rust
pub struct SpiceExportConfig {
    pub model_type: SpiceModelType,
    pub include_dc_resistance: bool,
    pub include_mutual_inductance: bool,
    pub include_mutual_capacitance: bool,
    pub include_dielectric_loss: bool,
    pub coupling_threshold: f64,
}

pub enum SpiceModelType {
    BroadbandLumped,
    FrequencyDependentLumped,
    TLineModel,
    SParameterBlock,
}

pub fn export_spice(
    rlcg: &RlcgMatrixData,
    config: &SpiceExportConfig,
    output_path: &Path,
) -> Result<(), ExportError> {
    // Open output file
    let mut file = std::fs::File::create(output_path)?;
    
    // Write header
    writeln!(file, "* EMStudio Q3D Equivalent Circuit Export")?;
    writeln!(file, "* Generated: {}", chrono::Utc::now())?;
    
    // Write subcircuit
    write!(file, ".SUBCKT Q3D_Model ")?;
    for (i, name) in rlcg.terminal_names.iter().enumerate() {
        write!(file, "{} ", name)?;
    }
    writeln!(file, "GND\n")?;
    
    // Write self impedance (R + L for each terminal pair)
    for i in 0..rlcg.num_terminals {
        let r_val = rlcg.dc_r[i][i];
        let l_val = rlcg.dc_l[i][i];
        let term_src = &rlcg.terminal_names[i];
        
        writeln!(file, "* Self impedance: {}", term_src)?;
        writeln!(file, "R_self_{}  {}_src  n{}_1  {}", i, term_src, i, r_val)?;
        writeln!(file, "L_self_{}  n{}_1   {}_sink  {}n", i, i, term_src, l_val)?;
    }
    
    // Write mutual inductances
    for i in 0..rlcg.num_terminals {
        for j in (i+1)..rlcg.num_terminals {
            let coupling = rlcg.dc_l[i][j] / (rlcg.dc_l[i][i] * rlcg.dc_l[j][j]).sqrt();
            if coupling.abs() > config.coupling_threshold {
                writeln!(file, "* Mutual inductance: {} <-> {}", i, j)?;
                writeln!(file, "K_{}{}  L_self_{}  L_self_{}  {}", i, j, i, j, coupling)?;
            }
        }
    }
    
    // Write capacitances
    for i in 0..rlcg.num_terminals {
        let c_to_gnd = rlcg.c_matrix[0][i][i];
        writeln!(file, "* Capacitance: {} to GND", i)?;
        writeln!(file, "C_{}_g  {}_src  GND  {}p", i, &rlcg.terminal_names[i], c_to_gnd)?;
    }
    
    // Write mutual capacitances
    for i in 0..rlcg.num_terminals {
        for j in (i+1)..rlcg.num_terminals {
            let c_mutual = rlcg.c_matrix[0][i][j].abs();
            if c_mutual > config.coupling_threshold {
                writeln!(file, "* Mutual capacitance: {} <-> {}", i, j)?;
                writeln!(file, "C_{}{}  {}_src  {}_src  {}p", i, j,
                    &rlcg.terminal_names[i], &rlcg.terminal_names[j], c_mutual)?;
            }
        }
    }
    
    writeln!(file, "\n.ENDS Q3D_Model")?;
    Ok(())
}
```

### VTK Export

```rust
pub fn export_vtk(
    mesh: &MshMesh,
    field_values: &[f64],
    output_path: &Path,
) -> Result<(), ExportError>;

pub fn export_vtk_xml(
    mesh: &MshMesh,
    field_values: &[f64],
    output_path: &Path,
) -> Result<(), ExportError> {
    // Generate VTK XML unstructured grid (.vtu)
    // Cell type: 10 (tetrahedron)
    // Point data: field_values
    todo!()
}
```

### Phase Animation

```rust
pub struct PhaseAnimator {
    pub current_phase_deg: f64,
    pub phase_step_deg: f64,
    pub playing: bool,
    pub fps: f64,
    last_frame: std::time::Instant,
}

impl PhaseAnimator {
    pub fn new(phase_step_deg: f64, fps: f64) -> Self {
        Self {
            current_phase_deg: 0.0,
            phase_step_deg,
            playing: false,
            fps,
            last_frame: std::time::Instant::now(),
        }
    }
    
    /// Evaluate real field from complex components at current phase
    pub fn evaluate(
        &self,
        field_real: &[f64],
        field_imag: &[f64],
    ) -> Vec<f64> {
        let phase_rad = self.current_phase_deg.to_radians();
        let cos_p = phase_rad.cos();
        let sin_p = phase_rad.sin();
        
        field_real.iter()
            .zip(field_imag.iter())
            .map(|(re, im)| re * cos_p - im * sin_p)
            .collect()
    }
    
    /// Advance one frame (if playing and time elapsed)
    pub fn tick(&mut self) {
        if self.playing {
            let now = std::time::Instant::now();
            let elapsed = now.duration_since(self.last_frame).as_secs_f64();
            if elapsed >= 1.0 / self.fps {
                self.current_phase_deg =
                    (self.current_phase_deg + self.phase_step_deg) % 360.0;
                self.last_frame = now;
            }
        }
    }
}
```

