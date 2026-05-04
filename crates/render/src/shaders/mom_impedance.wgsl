// MoM impedance matrix element compute shader.
//
// Computes Z[i,j] = ∫∫ f_i(r) · f_j(r') G(r,r') dS' dS
// for RWG basis function pairs on triangle surface meshes.
//
// Each workgroup computes a tile of the impedance matrix.

struct Edge {
    center: vec3<f32>,
    tangent: vec3<f32>,
    length: f32,
    div_rho: f32,
}

struct Triangle {
    v0: vec3<f32>,
    v1: vec3<f32>,
    v2: vec3<f32>,
    centroid: vec3<f32>,
    normal: vec3<f32>,
    area: f32,
}

@group(0) @binding(0) var<storage, read> edges: array<Edge>;
@group(0) @binding(1) var<storage, read> triangles: array<Triangle>;
@group(0) @binding(2) var<storage, read_write> z_matrix: array<vec2<f32>>; // complex: (re, im)
@group(0) @binding(3) var<uniform> params: MoMParams;

struct MoMParams {
    n_edges: u32,
    n_triangles: u32,
    k0: f32,        // wavenumber [rad/m]
    eta0: f32,      // free-space impedance [Ω]
    quad_order: u32, // Gaussian quadrature order
}

// Free-space Green's function: G = exp(-j k0 R) / (4π R)
fn green_free_space(r: vec3<f32>, r_prime: vec3<f32>, k0: f32) -> vec2<f32> {
    let dr = r - r_prime;
    let R = length(dr);
    if R < 1e-12 {
        return vec2(0.0, 0.0);
    }
    let phase = -k0 * R;
    let g_re = cos(phase) / (4.0 * 3.14159265359 * R);
    let g_im = sin(phase) / (4.0 * 3.14159265359 * R);
    return vec2(g_re, g_im);
}

// Complex multiply: (a_re + j a_im) * (b_re + j b_im)
fn cmul(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
    return vec2(
        a.x * b.x - a.y * b.y,
        a.x * b.y + a.y * b.x,
    );
}

// RWG basis function f_i evaluated at point r
fn rwg_eval(r: vec3<f32>, edge_idx: u32, sign: u32) -> vec3<f32> {
    let e = edges[edge_idx];
    // f(r) = ± (r - vertex_free) * length / (2 * area)
    // Simplified: return tangent direction for now
    let s = select(-1.0, 1.0, sign == 1u);
    return s * e.tangent * e.length;
}

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;

    if i >= params.n_edges || j >= params.n_edges {
        return;
    }

    // Default quadrature: single-point centroid evaluation for speed
    // Full quadrature would loop over Gaussian points on each triangle

    // Evaluate Z[i,j] = dot(f_i, f_j) * G(r_i, r_j) * area_i * area_j
    // Simplified: use edge center as evaluation point
    let r_i = edges[i].center;
    let r_j = edges[j].center;

    let g = green_free_space(r_i, r_j, params.k0);

    // Dot product of basis functions (simplified)
    let dot_ij = dot(edges[i].tangent, edges[j].tangent)
               * edges[i].length * edges[j].length;

    // Z[i,j] = j * k0 * eta0 * dot_ij * G
    let coeff = vec2(0.0, params.k0 * params.eta0);
    let z_ij = cmul(cmul(coeff, vec2(dot_ij, 0.0)), g);

    let idx = i * params.n_edges + j;
    z_matrix[idx] = z_ij;
}
