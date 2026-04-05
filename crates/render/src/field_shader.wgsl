// ============================================================
// Scene pass: field-colored mesh with Blinn-Phong lighting
// ============================================================

struct Uniforms {
    mvp:       mat4x4<f32>,
    eye_pos:   vec3<f32>,
    _pad0:     f32,
    light_dir: vec3<f32>,
    _pad1:     f32,
    field_min: f32,
    field_max: f32,
    opacity:   f32,
    _pad2:     f32,
};

@group(0) @binding(0) var<uniform> u: Uniforms;
@group(0) @binding(1) var colormap_tex: texture_2d<f32>;
@group(0) @binding(2) var colormap_samp: sampler;

struct VsOut {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) world_normal: vec3<f32>,
    @location(1) world_pos:    vec3<f32>,
    @location(2) field_t:      f32,
};

@vertex
fn vs_main(
    @location(0) position:    vec3<f32>,
    @location(1) normal:      vec3<f32>,
    @location(2) field_value: f32,
) -> VsOut {
    var out: VsOut;
    out.clip_pos     = u.mvp * vec4<f32>(position, 1.0);
    out.world_normal = normal;
    out.world_pos    = position;
    let range = u.field_max - u.field_min;
    if range > 0.0 {
        out.field_t = clamp((field_value - u.field_min) / range, 0.0, 1.0);
    } else {
        out.field_t = 0.5;
    }
    return out;
}

@fragment
fn fs_main(in: VsOut) -> @location(0) vec4<f32> {
    let n = normalize(in.world_normal);
    let l = normalize(u.light_dir);
    let v = normalize(u.eye_pos - in.world_pos);
    let h = normalize(l + v);

    let ambient  = 0.15;
    let diffuse  = max(dot(n, l), 0.0) * 0.7;
    let specular = pow(max(dot(n, h), 0.0), 32.0) * 0.3;

    // Sample colormap (256x1 2D texture, sample at (t, 0.5))
    let color = textureSampleLevel(colormap_tex, colormap_samp, vec2<f32>(in.field_t, 0.5), 0.0).rgb;
    let lit = color * (ambient + diffuse) + vec3<f32>(specular);

    return vec4<f32>(lit, u.opacity);
}

// ============================================================
// Wireframe pass: solid dark lines
// ============================================================

@fragment
fn fs_wire(in: VsOut) -> @location(0) vec4<f32> {
    return vec4<f32>(0.1, 0.1, 0.1, 0.6);
}

// ============================================================
// Arrow instanced pass: 3D arrows colored by magnitude
// ============================================================
// Instance attributes:
//   @location(3) inst_position: vec3<f32>
//   @location(4) inst_direction: vec3<f32>  (unit vector)
//   @location(5) inst_magnitude: f32        (0..1 normalized)

@vertex
fn vs_arrow(
    @location(0) vert_pos: vec3<f32>,       // base arrow mesh vertex (oriented along +Y)
    @location(1) vert_normal: vec3<f32>,
    @location(2) _vert_field: f32,
    @location(3) inst_position: vec3<f32>,
    @location(4) inst_direction: vec3<f32>,
    @location(5) inst_magnitude: f32,
) -> VsOut {
    // Build rotation matrix from +Y to inst_direction
    let up = vec3<f32>(0.0, 1.0, 0.0);
    let dir = normalize(inst_direction);
    let cos_a = dot(up, dir);

    var world_pos: vec3<f32>;
    var world_normal: vec3<f32>;

    let arrow_scale = 0.3 * (0.3 + 0.7 * inst_magnitude);

    if cos_a > 0.9999 {
        // Nearly aligned with +Y, no rotation needed
        world_pos = vert_pos * arrow_scale + inst_position;
        world_normal = vert_normal;
    } else if cos_a < -0.9999 {
        // Opposite direction: flip
        world_pos = vec3<f32>(vert_pos.x, -vert_pos.y, -vert_pos.z) * arrow_scale + inst_position;
        world_normal = vec3<f32>(vert_normal.x, -vert_normal.y, -vert_normal.z);
    } else {
        let axis = normalize(cross(up, dir));
        let sin_a = sqrt(1.0 - cos_a * cos_a);
        // Rodrigues rotation
        let scaled = vert_pos * arrow_scale;
        world_pos = scaled * cos_a + cross(axis, scaled) * sin_a + axis * dot(axis, scaled) * (1.0 - cos_a) + inst_position;
        world_normal = vert_normal * cos_a + cross(axis, vert_normal) * sin_a + axis * dot(axis, vert_normal) * (1.0 - cos_a);
    }

    var out: VsOut;
    out.clip_pos = u.mvp * vec4<f32>(world_pos, 1.0);
    out.world_normal = normalize(world_normal);
    out.world_pos = world_pos;
    out.field_t = inst_magnitude;
    return out;
}

// ============================================================
// Blit pass: fullscreen quad to copy offscreen texture to screen
// ============================================================

struct BlitVsOut {
    @builtin(position) pos: vec4<f32>,
    @location(0) uv: vec2<f32>,
};

@vertex
fn vs_blit(@builtin(vertex_index) vi: u32) -> BlitVsOut {
    // Full-screen triangle (3 vertices cover the whole screen)
    var positions = array<vec2<f32>, 3>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 3.0, -1.0),
        vec2<f32>(-1.0,  3.0),
    );
    var uvs = array<vec2<f32>, 3>(
        vec2<f32>(0.0, 1.0),
        vec2<f32>(2.0, 1.0),
        vec2<f32>(0.0, -1.0),
    );
    var out: BlitVsOut;
    out.pos = vec4<f32>(positions[vi], 0.0, 1.0);
    out.uv  = uvs[vi];
    return out;
}

@group(0) @binding(0) var blit_tex:  texture_2d<f32>;
@group(0) @binding(1) var blit_samp: sampler;

@fragment
fn fs_blit(in: BlitVsOut) -> @location(0) vec4<f32> {
    return textureSampleLevel(blit_tex, blit_samp, in.uv, 0.0);
}
