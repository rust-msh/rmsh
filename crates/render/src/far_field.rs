use crate::mesh_data::{FieldMesh, FieldVertex};
use emstudio_domain::result_store::FarFieldData;
use std::collections::HashSet;

/// Generate a 3D far-field radiation pattern mesh.
///
/// The pattern is a sphere whose radius at each (theta, phi) is modulated
/// by a gain function. The field_value is set to the gain for colormap display.
pub fn generate_pattern_mesh(
    n_theta: u32,
    n_phi: u32,
    gain_fn: &dyn Fn(f32, f32) -> f32,
) -> FieldMesh {
    let mut vertices = Vec::new();
    let mut field_min = f32::MAX;
    let mut field_max = f32::MIN;

    for it in 0..=n_theta {
        let theta = std::f32::consts::PI * it as f32 / n_theta as f32;
        let sin_t = theta.sin();
        let cos_t = theta.cos();

        for ip in 0..=n_phi {
            let phi = 2.0 * std::f32::consts::PI * ip as f32 / n_phi as f32;
            let sin_p = phi.sin();
            let cos_p = phi.cos();

            let gain_dbi = gain_fn(theta, phi);
            // Convert gain dBi to radius: scale so peak ≈ 1.0, floor at -30 dB
            let radius = ((gain_dbi + 30.0) / 30.0).clamp(0.05, 2.0);

            field_min = field_min.min(gain_dbi);
            field_max = field_max.max(gain_dbi);

            let x = radius * sin_t * cos_p;
            let y = radius * cos_t;
            let z = radius * sin_t * sin_p;

            // Normal approximation: use radial direction
            let n_len = (x * x + y * y + z * z).sqrt().max(1e-6);
            vertices.push(FieldVertex {
                position: [x, y, z],
                normal: [x / n_len, y / n_len, z / n_len],
                field_value: gain_dbi,
            });
        }
    }

    // Triangle indices
    let row_len = n_phi + 1;
    let mut indices = Vec::new();
    for it in 0..n_theta {
        for ip in 0..n_phi {
            let i00 = it * row_len + ip;
            let i01 = i00 + 1;
            let i10 = (it + 1) * row_len + ip;
            let i11 = i10 + 1;
            indices.extend_from_slice(&[i00, i10, i01]);
            indices.extend_from_slice(&[i01, i10, i11]);
        }
    }

    // Wireframe (sparse — every 5th edge to avoid clutter)
    let mut edge_set: HashSet<(u32, u32)> = HashSet::new();
    let mut wire_indices = Vec::new();
    for (ti, tri) in indices.chunks(3).enumerate() {
        if ti % 5 != 0 {
            continue;
        }
        for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            if edge_set.insert(key) {
                wire_indices.push(a);
                wire_indices.push(b);
            }
        }
    }

    FieldMesh {
        vertices,
        indices,
        wire_indices,
        field_range: [field_min, field_max],
        field_imag: None,
        vector_field: None,
    }
}

/// Synthetic dipole-like radiation pattern.
/// Main lobe along +Y (theta=0), with side lobes.
pub fn dipole_gain(theta: f32, _phi: f32) -> f32 {
    let cos_t = theta.cos();
    // Dipole: gain ∝ sin²(theta), but let's make it more interesting
    let main = cos_t.powi(4);
    let gain_linear = main.max(0.001);
    10.0 * gain_linear.log10() // Convert to dBi (approximate)
}

/// Synthetic patch antenna pattern with directional main beam.
pub fn patch_gain(theta: f32, phi: f32) -> f32 {
    let cos_t = theta.cos();
    let cos_p = phi.cos();
    // Patch-like: main beam at theta=0, with cos taper and some phi variation
    let main = cos_t.powi(3).max(0.0);
    let phi_factor = 1.0 + 0.3 * (2.0 * cos_p).powi(2);
    let gain_linear = (main * phi_factor).max(0.001);
    10.0 * gain_linear.log10()
}

/// Generate a 3D far-field radiation pattern mesh from real measured/simulated data.
///
/// Reads a derived quantity (e.g., "GainTotal") from the FarFieldData and uses it
/// as the gain function for the mesh generation.
pub fn generate_pattern_mesh_from_data(
    data: &FarFieldData,
    quantity: &str,
) -> Option<FieldMesh> {
    let dq = data.derived_quantities.get(quantity)?;
    let n_theta = data.theta.num_points;
    let n_phi = data.phi.num_points;

    if dq.data.len() != n_theta * n_phi || n_theta == 0 || n_phi == 0 {
        return None;
    }

    let mut vertices = Vec::new();
    let mut field_min = f32::MAX;
    let mut field_max = f32::MIN;

    for it in 0..n_theta {
        let theta = (data.theta.start_deg + it as f64 * data.theta.step_deg).to_radians() as f32;
        let sin_t = theta.sin();
        let cos_t = theta.cos();

        for ip in 0..n_phi {
            let phi = (data.phi.start_deg + ip as f64 * data.phi.step_deg).to_radians() as f32;
            let sin_p = phi.sin();
            let cos_p = phi.cos();

            let idx = it * n_phi + ip;
            let gain_dbi = dq.data[idx] as f32;
            let radius = ((gain_dbi + 30.0) / 30.0).clamp(0.05, 2.0);

            field_min = field_min.min(gain_dbi);
            field_max = field_max.max(gain_dbi);

            let x = radius * sin_t * cos_p;
            let y = radius * cos_t;
            let z = radius * sin_t * sin_p;

            let n_len = (x * x + y * y + z * z).sqrt().max(1e-6);
            vertices.push(FieldVertex {
                position: [x, y, z],
                normal: [x / n_len, y / n_len, z / n_len],
                field_value: gain_dbi,
            });
        }
    }

    // Triangle indices
    let row_len = n_phi as u32;
    let mut indices = Vec::new();
    for it in 0..n_theta as u32 - 1 {
        for ip in 0..n_phi as u32 - 1 {
            let i00 = it * row_len + ip;
            let i01 = i00 + 1;
            let i10 = (it + 1) * row_len + ip;
            let i11 = i10 + 1;
            indices.extend_from_slice(&[i00, i10, i01]);
            indices.extend_from_slice(&[i01, i10, i11]);
        }
    }

    // Sparse wireframe
    let mut edge_set: HashSet<(u32, u32)> = HashSet::new();
    let mut wire_indices = Vec::new();
    for (ti, tri) in indices.chunks(3).enumerate() {
        if ti % 5 != 0 {
            continue;
        }
        for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            let key = if a < b { (a, b) } else { (b, a) };
            if edge_set.insert(key) {
                wire_indices.push(a);
                wire_indices.push(b);
            }
        }
    }

    if field_min > field_max {
        field_min = 0.0;
        field_max = 1.0;
    }

    Some(FieldMesh {
        vertices,
        indices,
        wire_indices,
        field_range: [field_min, field_max],
        field_imag: None,
        vector_field: None,
    })
}
