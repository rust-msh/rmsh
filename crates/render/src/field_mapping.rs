// ---------------------------------------------------------------------------
// Field Mapping — Map emsfld field data to FieldMesh vertex colors
// ---------------------------------------------------------------------------
//
// Bridges the field data loader (emsfld_loader) with the render pipeline
// (FieldMesh). Maps complex vector/scalar field values to per-vertex
// scalar field_value for colormap rendering.

use std::collections::HashMap;

use emstudio_domain::emsfld_loader::FieldBlock;

use crate::mesh_data::FieldMesh;

/// Which component/quantity to extract from the field data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldComponent {
    /// sqrt(|Fx|^2 + |Fy|^2 + |Fz|^2) for vector, |F| for scalar
    Magnitude,
    /// Real part of the selected component (or scalar)
    RealPart,
    /// Imaginary part
    ImagPart,
    /// Phase angle in degrees
    Phase,
    /// X component magnitude (vector fields only)
    ComponentX,
    /// Y component magnitude
    ComponentY,
    /// Z component magnitude
    ComponentZ,
    /// Poynting vector magnitude |E × H*| / 2
    Poynting,
    /// Specific Absorption Rate [W/kg]
    SAR,
    /// Ohmic loss density |J|²/σ [W/m³]
    OhmicLoss,
}

/// Compute a derived field quantity at a single node.
///
/// `e` = (Ex_re, Ex_im, Ey_re, Ey_im, Ez_re, Ez_im)
/// `h` = (Hx_re, Hx_im, Hy_re, Hy_im, Hz_re, Hz_im)
pub fn compute_derived_field(
    component: FieldComponent,
    e: Option<(f64, f64, f64, f64, f64, f64)>,
    h: Option<(f64, f64, f64, f64, f64, f64)>,
    sigma: f64,
    _rho: f64,
) -> f64 {
    match component {
        FieldComponent::Poynting => {
            let e = match e {
                Some(e) => e,
                None => return 0.0,
            };
            let h = match h {
                Some(h) => h,
                None => return 0.0,
            };
            // S = ½ Re(E × H*)
            let sx = e.2 * h.4 + e.3 * h.5 - (e.4 * h.2 + e.5 * h.3);
            let sy = e.4 * h.0 + e.5 * h.1 - (e.0 * h.4 + e.1 * h.5);
            let sz = e.0 * h.2 + e.1 * h.3 - (e.2 * h.0 + e.3 * h.1);
            0.5 * (sx.powi(2) + sy.powi(2) + sz.powi(2)).sqrt()
        }
        FieldComponent::SAR => {
            let e = match e {
                Some(e) => e,
                None => return 0.0,
            };
            if sigma <= 0.0 {
                return 0.0;
            }
            let e2 = e.0.powi(2) + e.1.powi(2) + e.2.powi(2) + e.3.powi(2) + e.4.powi(2) + e.5.powi(2);
            sigma * e2 / (2.0 * _rho.max(1.0))
        }
        FieldComponent::OhmicLoss => {
            let e = match e {
                Some(e) => e,
                None => return 0.0,
            };
            let e2 = e.0.powi(2) + e.1.powi(2) + e.2.powi(2) + e.3.powi(2) + e.4.powi(2) + e.5.powi(2);
            sigma * e2
        }
        _ => 0.0,
    }
}

/// Map emsfld field data onto a FieldMesh's vertex field_value.
///
/// `node_to_vertex` maps mesh node tags (from .msh) to FieldMesh vertex indices.
/// Field data is node-indexed (node 0 = first node in .msh), vertex indices may
/// differ due to surface extraction.
///
/// Returns the field range [min, max].
pub fn map_field_to_mesh(
    mesh: &mut FieldMesh,
    field: &FieldBlock,
    component: FieldComponent,
    node_to_vertex: &HashMap<u64, u32>,
) -> [f32; 2] {
    let mut field_min = f32::MAX;
    let mut field_max = f32::MIN;

    // Also fill field_imag for phase animation support
    let mut imag_values = vec![0.0f32; mesh.vertices.len()];
    let mut has_imag = false;

    for (&node_tag, &vertex_idx) in node_to_vertex {
        let vi = vertex_idx as usize;
        if vi >= mesh.vertices.len() {
            continue;
        }

        // node_tag is 1-based in MSH, field data is 0-based
        let node_idx = (node_tag - 1) as usize;
        if node_idx >= field.num_nodes {
            continue;
        }

        let val = extract_value(field, node_idx, component);
        mesh.vertices[vi].field_value = val as f32;
        field_min = field_min.min(val as f32);
        field_max = field_max.max(val as f32);

        // Extract imaginary for phase animation
        if field.num_components >= 3 {
            let base = node_idx * field.num_components;
            if let (Some(vx), Some(vy), Some(vz)) =
                (field.data.get(base), field.data.get(base + 1), field.data.get(base + 2))
            {
                let imag_mag = (vx.imag * vx.imag + vy.imag * vy.imag + vz.imag * vz.imag).sqrt();
                imag_values[vi] = imag_mag as f32;
                if imag_mag.abs() > 1e-20 {
                    has_imag = true;
                }
            }
        } else if field.num_components == 1 {
            if let Some(cv) = field.data.get(node_idx) {
                imag_values[vi] = cv.imag as f32;
                if cv.imag.abs() > 1e-20 {
                    has_imag = true;
                }
            }
        }
    }

    if field_min > field_max {
        field_min = 0.0;
        field_max = 1.0;
    }

    mesh.field_range = [field_min, field_max];
    if has_imag {
        mesh.field_imag = Some(imag_values);
    }

    [field_min, field_max]
}

/// Map field data to a FieldMesh using direct node index mapping
/// (when the mesh vertices are in the same order as the field data nodes).
pub fn map_field_direct(
    mesh: &mut FieldMesh,
    field: &FieldBlock,
    component: FieldComponent,
) -> [f32; 2] {
    let mut field_min = f32::MAX;
    let mut field_max = f32::MIN;
    let mut imag_values = vec![0.0f32; mesh.vertices.len()];
    let mut has_imag = false;

    for (vi, vertex) in mesh.vertices.iter_mut().enumerate() {
        if vi >= field.num_nodes {
            break;
        }

        let val = extract_value(field, vi, component);
        vertex.field_value = val as f32;
        field_min = field_min.min(val as f32);
        field_max = field_max.max(val as f32);

        // Imaginary for phase animation
        if field.num_components == 1 {
            if let Some(cv) = field.data.get(vi) {
                imag_values[vi] = cv.imag as f32;
                if cv.imag.abs() > 1e-20 {
                    has_imag = true;
                }
            }
        }
    }

    if field_min > field_max {
        field_min = 0.0;
        field_max = 1.0;
    }

    mesh.field_range = [field_min, field_max];
    if has_imag {
        mesh.field_imag = Some(imag_values);
    }

    [field_min, field_max]
}

/// Also extract 3D vector field data for arrow visualization.
pub fn map_vector_field(
    mesh: &mut FieldMesh,
    field: &FieldBlock,
    node_to_vertex: &HashMap<u64, u32>,
) {
    if field.num_components < 3 {
        return;
    }

    let mut vectors = vec![[0.0f32; 3]; mesh.vertices.len()];

    for (&node_tag, &vertex_idx) in node_to_vertex {
        let vi = vertex_idx as usize;
        if vi >= mesh.vertices.len() {
            continue;
        }
        let node_idx = (node_tag - 1) as usize;
        if node_idx >= field.num_nodes {
            continue;
        }

        let base = node_idx * field.num_components;
        if let (Some(vx), Some(vy), Some(vz)) = (
            field.data.get(base),
            field.data.get(base + 1),
            field.data.get(base + 2),
        ) {
            vectors[vi] = [
                vx.magnitude() as f32,
                vy.magnitude() as f32,
                vz.magnitude() as f32,
            ];
        }
    }

    mesh.vector_field = Some(vectors);
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn extract_value(field: &FieldBlock, node_idx: usize, component: FieldComponent) -> f64 {
    match component {
        FieldComponent::Magnitude => {
            if field.num_components >= 3 {
                field.vector_magnitude(node_idx).unwrap_or(0.0)
            } else {
                field
                    .data
                    .get(node_idx)
                    .map(|cv| cv.magnitude())
                    .unwrap_or(0.0)
            }
        }
        FieldComponent::RealPart => {
            if field.num_components >= 3 {
                // Real part of vector magnitude
                let base = node_idx * field.num_components;
                let rx = field.data.get(base).map(|c| c.real).unwrap_or(0.0);
                let ry = field.data.get(base + 1).map(|c| c.real).unwrap_or(0.0);
                let rz = field.data.get(base + 2).map(|c| c.real).unwrap_or(0.0);
                (rx * rx + ry * ry + rz * rz).sqrt()
            } else {
                field.data.get(node_idx).map(|c| c.real).unwrap_or(0.0)
            }
        }
        FieldComponent::ImagPart => {
            if field.num_components >= 3 {
                let base = node_idx * field.num_components;
                let ix = field.data.get(base).map(|c| c.imag).unwrap_or(0.0);
                let iy = field.data.get(base + 1).map(|c| c.imag).unwrap_or(0.0);
                let iz = field.data.get(base + 2).map(|c| c.imag).unwrap_or(0.0);
                (ix * ix + iy * iy + iz * iz).sqrt()
            } else {
                field.data.get(node_idx).map(|c| c.imag).unwrap_or(0.0)
            }
        }
        FieldComponent::Phase => field
            .data
            .get(node_idx * field.num_components)
            .map(|c| c.phase_deg())
            .unwrap_or(0.0),
        FieldComponent::ComponentX => {
            let base = node_idx * field.num_components;
            field.data.get(base).map(|c| c.magnitude()).unwrap_or(0.0)
        }
        FieldComponent::ComponentY => {
            let base = node_idx * field.num_components;
            field
                .data
                .get(base + 1)
                .map(|c| c.magnitude())
                .unwrap_or(0.0)
        }
        FieldComponent::ComponentZ => {
            let base = node_idx * field.num_components;
            field
                .data
                .get(base + 2)
                .map(|c| c.magnitude())
                .unwrap_or(0.0)
        }
        // Derived quantities (Poynting, SAR, OhmicLoss) require both E and H
        // field data — use `compute_derived_field()` with both field blocks.
        FieldComponent::Poynting | FieldComponent::SAR | FieldComponent::OhmicLoss => {
            // Single-field extraction: return magnitude as fallback
            if field.num_components >= 3 {
                field.vector_magnitude(node_idx).unwrap_or(0.0)
            } else {
                field.data.get(node_idx).map(|cv| cv.magnitude()).unwrap_or(0.0)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_data::FieldVertex;
    use emstudio_domain::emsfld_loader::ComplexValue;

    fn make_test_field_mesh() -> FieldMesh {
        FieldMesh {
            vertices: vec![
                FieldVertex {
                    position: [0.0, 0.0, 0.0],
                    normal: [0.0, 0.0, 1.0],
                    field_value: 0.0,
                },
                FieldVertex {
                    position: [1.0, 0.0, 0.0],
                    normal: [0.0, 0.0, 1.0],
                    field_value: 0.0,
                },
                FieldVertex {
                    position: [0.0, 1.0, 0.0],
                    normal: [0.0, 0.0, 1.0],
                    field_value: 0.0,
                },
            ],
            indices: vec![0, 1, 2],
            wire_indices: vec![0, 1, 1, 2, 2, 0],
            field_range: [0.0, 0.0],
            field_imag: None,
            vector_field: None,
        }
    }

    fn make_test_field_block() -> FieldBlock {
        // 3-component vector field for 3 nodes
        FieldBlock {
            frequency: 2.4,
            data: vec![
                // Node 0: (3+4i, 0, 0)
                ComplexValue { real: 3.0, imag: 4.0 },
                ComplexValue { real: 0.0, imag: 0.0 },
                ComplexValue { real: 0.0, imag: 0.0 },
                // Node 1: (0, 1+0i, 0)
                ComplexValue { real: 0.0, imag: 0.0 },
                ComplexValue { real: 1.0, imag: 0.0 },
                ComplexValue { real: 0.0, imag: 0.0 },
                // Node 2: (0, 0, 2+0i)
                ComplexValue { real: 0.0, imag: 0.0 },
                ComplexValue { real: 0.0, imag: 0.0 },
                ComplexValue { real: 2.0, imag: 0.0 },
            ],
            num_nodes: 3,
            num_components: 3,
        }
    }

    #[test]
    fn map_magnitude_direct() {
        let mut mesh = make_test_field_mesh();
        let field = make_test_field_block();

        let range = map_field_direct(&mut mesh, &field, FieldComponent::Magnitude);

        // Node 0: magnitude = sqrt(5^2 + 0 + 0) = 5.0
        assert!((mesh.vertices[0].field_value - 5.0).abs() < 1e-4);
        // Node 1: magnitude = sqrt(0 + 1^2 + 0) = 1.0
        assert!((mesh.vertices[1].field_value - 1.0).abs() < 1e-4);
        // Node 2: magnitude = sqrt(0 + 0 + 2^2) = 2.0
        assert!((mesh.vertices[2].field_value - 2.0).abs() < 1e-4);

        assert!((range[0] - 1.0).abs() < 1e-4); // min
        assert!((range[1] - 5.0).abs() < 1e-4); // max
    }

    #[test]
    fn map_component_x() {
        let mut mesh = make_test_field_mesh();
        let field = make_test_field_block();

        map_field_direct(&mut mesh, &field, FieldComponent::ComponentX);

        // Node 0: |Fx| = |3+4i| = 5.0
        assert!((mesh.vertices[0].field_value - 5.0).abs() < 1e-4);
        // Node 1: |Fx| = 0
        assert!((mesh.vertices[1].field_value).abs() < 1e-4);
    }

    #[test]
    fn map_with_node_mapping() {
        let mut mesh = make_test_field_mesh();
        let field = make_test_field_block();

        // Map: MSH node tag 1 → vertex 0, tag 2 → vertex 1, tag 3 → vertex 2
        let mut node_map = HashMap::new();
        node_map.insert(1u64, 0u32);
        node_map.insert(2u64, 1u32);
        node_map.insert(3u64, 2u32);

        let range = map_field_to_mesh(&mut mesh, &field, FieldComponent::Magnitude, &node_map);

        assert!((mesh.vertices[0].field_value - 5.0).abs() < 1e-4);
        assert!(range[1] > range[0]);
    }
}
