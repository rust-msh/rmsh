// ---------------------------------------------------------------------------
// VTK Export — ParaView-compatible VTK Legacy format writer
// ---------------------------------------------------------------------------

use std::io::Write;
use std::path::Path;

use crate::emsfld_loader::FieldBlock;
use crate::msh_loader::{self, MshMesh};

/// Export mesh and field data as VTK Legacy ASCII format (.vtk).
pub fn export_vtk_ascii(
    path: &Path,
    mesh: &MshMesh,
    field: Option<&FieldBlock>,
) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;

    // Collect tet elements only
    let tets: Vec<_> = mesh
        .elements
        .iter()
        .filter(|e| e.element_type == msh_loader::element_types::TET4)
        .collect();

    // Build node index mapping: msh tag → 0-based contiguous index
    let mut tag_to_idx: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();
    let mut used_nodes: Vec<(u64, [f64; 3])> = Vec::new();

    for tet in &tets {
        for &tag in &tet.node_tags[..4] {
            if !tag_to_idx.contains_key(&tag) {
                let idx = used_nodes.len();
                let pos = mesh.node_position(tag).unwrap_or([0.0; 3]);
                tag_to_idx.insert(tag, idx);
                used_nodes.push((tag, pos));
            }
        }
    }

    let num_points = used_nodes.len();
    let num_cells = tets.len();

    // Header
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "EMStudio Field Export")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET UNSTRUCTURED_GRID")?;

    // Points
    writeln!(file, "POINTS {} float", num_points)?;
    for &(_, pos) in &used_nodes {
        writeln!(file, "{} {} {}", pos[0], pos[1], pos[2])?;
    }

    // Cells: each tet has 4 nodes, plus 1 for the count = 5 ints per cell
    let total_ints = num_cells * 5;
    writeln!(file, "CELLS {} {}", num_cells, total_ints)?;
    for tet in &tets {
        let i0 = tag_to_idx[&tet.node_tags[0]];
        let i1 = tag_to_idx[&tet.node_tags[1]];
        let i2 = tag_to_idx[&tet.node_tags[2]];
        let i3 = tag_to_idx[&tet.node_tags[3]];
        writeln!(file, "4 {} {} {} {}", i0, i1, i2, i3)?;
    }

    // Cell types: VTK_TETRA = 10
    writeln!(file, "CELL_TYPES {}", num_cells)?;
    for _ in 0..num_cells {
        writeln!(file, "10")?;
    }

    // Point data (field values)
    if let Some(field) = field {
        writeln!(file, "POINT_DATA {}", num_points)?;

        if field.num_components >= 3 {
            // Write vector field magnitude as scalar
            writeln!(file, "SCALARS field_magnitude float 1")?;
            writeln!(file, "LOOKUP_TABLE default")?;
            for &(tag, _) in &used_nodes {
                let node_idx = (tag - 1) as usize;
                let mag = if node_idx < field.num_nodes {
                    field.vector_magnitude(node_idx).unwrap_or(0.0)
                } else {
                    0.0
                };
                writeln!(file, "{}", mag)?;
            }

            // Write vector field
            writeln!(file, "VECTORS field_vector float")?;
            for &(tag, _) in &used_nodes {
                let node_idx = (tag - 1) as usize;
                if node_idx < field.num_nodes {
                    let base = node_idx * field.num_components;
                    let vx = field.data.get(base).map(|c| c.magnitude()).unwrap_or(0.0);
                    let vy = field.data.get(base + 1).map(|c| c.magnitude()).unwrap_or(0.0);
                    let vz = field.data.get(base + 2).map(|c| c.magnitude()).unwrap_or(0.0);
                    writeln!(file, "{} {} {}", vx, vy, vz)?;
                } else {
                    writeln!(file, "0 0 0")?;
                }
            }
        } else {
            // Scalar field
            writeln!(file, "SCALARS field_value float 1")?;
            writeln!(file, "LOOKUP_TABLE default")?;
            for &(tag, _) in &used_nodes {
                let node_idx = (tag - 1) as usize;
                let val = if node_idx < field.num_nodes {
                    field.data.get(node_idx).map(|c| c.magnitude()).unwrap_or(0.0)
                } else {
                    0.0
                };
                writeln!(file, "{}", val)?;
            }
        }
    }

    Ok(())
}

/// Export mesh and field data as VTK Legacy binary format (.vtk).
pub fn export_vtk_binary(
    path: &Path,
    mesh: &MshMesh,
    field: Option<&FieldBlock>,
) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;

    let tets: Vec<_> = mesh
        .elements
        .iter()
        .filter(|e| e.element_type == msh_loader::element_types::TET4)
        .collect();

    let mut tag_to_idx: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();
    let mut used_nodes: Vec<(u64, [f64; 3])> = Vec::new();

    for tet in &tets {
        for &tag in &tet.node_tags[..4] {
            if !tag_to_idx.contains_key(&tag) {
                let idx = used_nodes.len();
                let pos = mesh.node_position(tag).unwrap_or([0.0; 3]);
                tag_to_idx.insert(tag, idx);
                used_nodes.push((tag, pos));
            }
        }
    }

    let num_points = used_nodes.len();
    let num_cells = tets.len();

    // Header (always ASCII)
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "EMStudio Field Export")?;
    writeln!(file, "BINARY")?;
    writeln!(file, "DATASET UNSTRUCTURED_GRID")?;

    // Points (big-endian float32)
    writeln!(file, "POINTS {} float", num_points)?;
    for &(_, pos) in &used_nodes {
        file.write_all(&(pos[0] as f32).to_be_bytes())?;
        file.write_all(&(pos[1] as f32).to_be_bytes())?;
        file.write_all(&(pos[2] as f32).to_be_bytes())?;
    }
    writeln!(file)?;

    // Cells
    let total_ints = num_cells * 5;
    writeln!(file, "CELLS {} {}", num_cells, total_ints)?;
    for tet in &tets {
        file.write_all(&4i32.to_be_bytes())?;
        file.write_all(&(tag_to_idx[&tet.node_tags[0]] as i32).to_be_bytes())?;
        file.write_all(&(tag_to_idx[&tet.node_tags[1]] as i32).to_be_bytes())?;
        file.write_all(&(tag_to_idx[&tet.node_tags[2]] as i32).to_be_bytes())?;
        file.write_all(&(tag_to_idx[&tet.node_tags[3]] as i32).to_be_bytes())?;
    }
    writeln!(file)?;

    // Cell types
    writeln!(file, "CELL_TYPES {}", num_cells)?;
    for _ in 0..num_cells {
        file.write_all(&10i32.to_be_bytes())?;
    }
    writeln!(file)?;

    // Point data
    if let Some(field) = field {
        writeln!(file, "POINT_DATA {}", num_points)?;
        writeln!(file, "SCALARS field_magnitude float 1")?;
        writeln!(file, "LOOKUP_TABLE default")?;
        for &(tag, _) in &used_nodes {
            let node_idx = (tag - 1) as usize;
            let mag = if node_idx < field.num_nodes {
                if field.num_components >= 3 {
                    field.vector_magnitude(node_idx).unwrap_or(0.0) as f32
                } else {
                    field.data.get(node_idx).map(|c| c.magnitude() as f32).unwrap_or(0.0)
                }
            } else {
                0.0
            };
            file.write_all(&mag.to_be_bytes())?;
        }
        writeln!(file)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::emsfld_loader::ComplexValue;
    use std::io::{BufReader, Cursor};

    fn make_test_mesh() -> MshMesh {
        let msh_data = "\
$MeshFormat\n\
4.1 0 8\n\
$EndMeshFormat\n\
$Nodes\n\
1 4 1 4\n\
3 1 0 4\n\
1\n\
2\n\
3\n\
4\n\
0.0 0.0 0.0\n\
1.0 0.0 0.0\n\
0.0 1.0 0.0\n\
0.0 0.0 1.0\n\
$EndNodes\n\
$Elements\n\
1 1 1 1\n\
3 1 4 1\n\
1 1 2 3 4\n\
$EndElements\n";
        let cursor = Cursor::new(msh_data.as_bytes());
        let mut reader = BufReader::new(cursor);
        MshMesh::read_from(&mut reader).unwrap()
    }

    #[test]
    fn vtk_ascii_export() {
        let mesh = make_test_mesh();
        let field = FieldBlock {
            frequency: 1.0,
            data: vec![
                ComplexValue { real: 1.0, imag: 0.0 },
                ComplexValue { real: 2.0, imag: 0.0 },
                ComplexValue { real: 3.0, imag: 0.0 },
                ComplexValue { real: 4.0, imag: 0.0 },
            ],
            num_nodes: 4,
            num_components: 1,
        };

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.vtk");
        export_vtk_ascii(&path, &mesh, Some(&field)).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("DATASET UNSTRUCTURED_GRID"));
        assert!(content.contains("POINTS 4 float"));
        assert!(content.contains("CELLS 1 5"));
        assert!(content.contains("SCALARS field_value float 1"));
    }
}
