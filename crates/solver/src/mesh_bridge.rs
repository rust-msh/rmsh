//! BRep → rmsh surface mesh → tetrahedral volume mesh → .msh file.
//!
//! This module bridges the emstudio rcad B-Rep geometry kernel with the rmsh
//! mesh generation toolkit. The conversion goes through raw coordinate arrays
//! to avoid type incompatibility between emstudio's rcad copy and rmsh's copy.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use rcad_kernel::BRep;
use rcad_render::Tessellator;
use rmsh_algo::{Delaunay3D, tetrahedralize_closed_surface};
use rmsh_algo::traits::{MeshParams as RmshMeshParams, Mesher3D};
use rmsh_model::{Element, ElementType, Mesh as RmshMesh, Node};

use emstudio_domain::Design;
use emstudio_domain::mesh::MeshOperationType;
use emstudio_domain::geometry::GeoObject;

use crate::error::SolverError;

// ---------------------------------------------------------------------------
// Intermediate types
// ---------------------------------------------------------------------------

/// Surface geometry extracted from a B-Rep object as raw coordinate arrays.
pub struct BRepSurfaceData {
    pub object_name: String,
    pub vertices: Vec<[f64; 3]>,
    pub triangles: Vec<[u64; 3]>,
    pub material: String,
}

/// Meshing parameters derived from a Design's analysis setup and mesh operations.
pub struct MeshConfig {
    /// Global maximum element size (derived from lambda/10 rule).
    pub max_element_size: f64,
    /// Per-object element size overrides from MeshOperation entries.
    pub object_overrides: HashMap<String, f64>,
    /// Design length unit → meter scale factor.
    pub unit_scale: f64,
}

/// Statistics about a generated mesh.
#[derive(Debug, Clone)]
pub struct MeshStats {
    pub num_nodes: usize,
    pub num_tetrahedra: usize,
    pub num_triangles: usize,
}

// ---------------------------------------------------------------------------
// BRep → surface extraction
// ---------------------------------------------------------------------------

/// Extract triangulated surface meshes from B-Rep objects using rcad-render's
/// Tessellator. The output is raw coordinate arrays, avoiding cross-crate
/// type incompatibility.
pub fn extract_brep_surfaces(
    breps: &HashMap<String, BRep>,
    objects: &[GeoObject],
) -> Vec<BRepSurfaceData> {
    let material_map: HashMap<&str, &str> = objects
        .iter()
        .map(|o| (o.name.as_str(), o.material.as_str()))
        .collect();

    breps
        .iter()
        .map(|(name, brep)| {
            let tessellated = Tessellator::tessellate(brep);

            let vertices: Vec<[f64; 3]> = tessellated
                .nodes
                .iter()
                .map(|v| [v[0] as f64, v[1] as f64, v[2] as f64])
                .collect();

            let triangles: Vec<[u64; 3]> = tessellated
                .indices
                .chunks(3)
                .map(|tri| [tri[0] as u64, tri[1] as u64, tri[2] as u64])
                .collect();

            let material = material_map
                .get(name.as_str())
                .unwrap_or(&"vacuum")
                .to_string();

            BRepSurfaceData {
                object_name: name.clone(),
                vertices,
                triangles,
                material,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Design → mesh parameters
// ---------------------------------------------------------------------------

/// Derive meshing parameters from a Design's analysis setup and mesh operations.
pub fn mesh_config_from_design(design: &Design) -> MeshConfig {
    let unit_scale = unit_to_meters(&design.units);

    // Find the first enabled analysis setup's solution frequency.
    let solution_freq_hz = design
        .analysis_setups
        .iter()
        .filter(|s| s.enabled)
        .find_map(|s| parse_frequency(&s.solution_frequency).ok())
        .unwrap_or(1.0e9); // default 1 GHz

    // lambda/10 rule: max element size = (c / freq) / 10, converted to design units
    let c0 = 3.0e8; // speed of light [m/s]
    let lambda_m = c0 / solution_freq_hz;
    let max_element_size = (lambda_m / 10.0) / unit_scale;

    // Per-object overrides from LengthBased mesh operations.
    let mut object_overrides = HashMap::new();
    for op in &design.mesh_operations {
        if op.mesh_type == MeshOperationType::LengthBased {
            if let Some(max_len) = op.properties.get("max_length").and_then(|v| v.as_f64()) {
                for target in &op.assignment.targets {
                    object_overrides.insert(target.clone(), max_len);
                }
            }
        }
    }

    MeshConfig {
        max_element_size,
        object_overrides,
        unit_scale,
    }
}

// ---------------------------------------------------------------------------
// Surface → volume mesh generation
// ---------------------------------------------------------------------------

/// Generate a tetrahedral volume mesh from extracted B-Rep surfaces, assign
/// physical groups, and save as Gmsh MSH v4.
///
/// Returns the path to the written .msh file and mesh statistics.
pub fn generate_mesh(
    surfaces: &[BRepSurfaceData],
    config: &MeshConfig,
    output_dir: &Path,
) -> Result<(PathBuf, MeshStats), SolverError> {
    if surfaces.is_empty() {
        return Err(SolverError::MeshGeneration(
            "No geometry objects to mesh".into(),
        ));
    }

    // Build a combined rmsh surface mesh from all objects.
    let mut combined = RmshMesh::new();
    let mut next_node_id: u64 = 1;
    let mut next_elem_id: u64 = 1;
    // Physical group tags: 1-based, one per object (for material assignment).
    let mut _physical_tag: i32 = 1;

    for surface in surfaces {
        let node_offset = next_node_id;

        // Add nodes with offset IDs.
        for (i, vtx) in surface.vertices.iter().enumerate() {
            combined.add_node(Node::new(
                node_offset + i as u64,
                vtx[0],
                vtx[1],
                vtx[2],
            ));
        }

        // Add triangle elements referencing offset node IDs.
        for tri in &surface.triangles {
            let mut elem = Element::new(
                next_elem_id,
                ElementType::Triangle3,
                vec![
                    tri[0] + node_offset,
                    tri[1] + node_offset,
                    tri[2] + node_offset,
                ],
            );
            elem.physical_tag = Some(_physical_tag);
            combined.add_element(elem);
            next_elem_id += 1;
        }

        next_node_id += surface.vertices.len() as u64;
        _physical_tag += 1;
    }

    // Determine element size: use global max, or smallest override if available.
    let element_size = if config.object_overrides.is_empty() {
        config.max_element_size
    } else {
        config
            .object_overrides
            .values()
            .copied()
            .fold(config.max_element_size, f64::min)
    };

    // Run 3D tetrahedral meshing.
    let volume_mesh = {
        // Try Delaunay3D first; fall back to centroid-star on failure.
        let delaunay = Delaunay3D::new();
        let params = RmshMeshParams::with_size(element_size);
        match delaunay.mesh_3d(&combined, &params) {
            Ok(m) => m,
            Err(_) => {
                // Fallback: simpler centroid-star algorithm.
                tetrahedralize_closed_surface(&combined)
                    .map_err(|e| SolverError::MeshGeneration(e.to_string()))?
            }
        }
    };

    let stats = MeshStats {
        num_nodes: volume_mesh.node_count(),
        num_tetrahedra: volume_mesh.elements_by_dimension(3).len(),
        num_triangles: volume_mesh.elements_by_dimension(2).len(),
    };

    // Save as .msh v4.
    let mesh_path = output_dir.join("mesh.msh");
    rmsh_io::save_msh_v4_to_path(&mesh_path, &volume_mesh)
        .map_err(|e| SolverError::MeshGeneration(format!("Failed to write .msh: {e}")))?;

    Ok((mesh_path, stats))
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Parse a frequency string like "2.4GHz", "900MHz", "100kHz", "50Hz".
pub fn parse_frequency(s: &str) -> Result<f64, SolverError> {
    let s = s.trim();

    // Try plain number first.
    if let Ok(v) = s.parse::<f64>() {
        return Ok(v);
    }

    let multipliers: &[(&str, f64)] = &[
        ("THz", 1e12),
        ("GHz", 1e9),
        ("MHz", 1e6),
        ("kHz", 1e3),
        ("Hz", 1.0),
    ];

    for (suffix, mult) in multipliers {
        if let Some(num_part) = s.strip_suffix(suffix) {
            return num_part
                .trim()
                .parse::<f64>()
                .map(|v| v * mult)
                .map_err(|_| {
                    SolverError::ConfigGeneration(format!("Invalid frequency: '{s}'"))
                });
        }
    }

    Err(SolverError::ConfigGeneration(format!(
        "Cannot parse frequency: '{s}'"
    )))
}

/// Map design unit string to meters.
pub fn unit_to_meters(unit: &str) -> f64 {
    match unit.to_lowercase().as_str() {
        "m" | "meter" | "meters" => 1.0,
        "cm" | "centimeter" => 1e-2,
        "mm" | "millimeter" => 1e-3,
        "um" | "micrometer" | "micron" => 1e-6,
        "nm" | "nanometer" => 1e-9,
        "in" | "inch" => 0.0254,
        "mil" => 2.54e-5,
        "ft" | "foot" => 0.3048,
        _ => 1e-3, // default mm
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_frequency_ghz() {
        let f = parse_frequency("2.4GHz").unwrap();
        assert!((f - 2.4e9).abs() < 1.0);
    }

    #[test]
    fn parse_frequency_mhz() {
        let f = parse_frequency("900MHz").unwrap();
        assert!((f - 9.0e8).abs() < 1.0);
    }

    #[test]
    fn parse_frequency_plain_number() {
        let f = parse_frequency("1e9").unwrap();
        assert!((f - 1e9).abs() < 1.0);
    }

    #[test]
    fn parse_frequency_invalid() {
        assert!(parse_frequency("abc").is_err());
    }

    #[test]
    fn unit_conversions() {
        assert!((unit_to_meters("mm") - 1e-3).abs() < 1e-15);
        assert!((unit_to_meters("cm") - 1e-2).abs() < 1e-15);
        assert!((unit_to_meters("m") - 1.0).abs() < 1e-15);
        assert!((unit_to_meters("in") - 0.0254).abs() < 1e-10);
        assert!((unit_to_meters("mil") - 2.54e-5).abs() < 1e-15);
    }
}
