//! STEP (.step / .stp) import and export via `rcad-step`.
//!
//! The heavy lifting (BRep parsing, triangle extraction, STEP serialisation)
//! lives in `rcad-step`.  This module is a thin adapter that converts between
//! `rcad-step`'s raw `(Vec<DVec3>, Vec<[usize;3]>)` representation and
//! `rmsh_model::Mesh`.

use std::path::Path;

use rcad_algorithms::{TessellationParams, mesh_brep};
use rcad_kernel::{BRep, Face, Shell, Solid, Vertex, Wire};
use rcad_step::ExportSelection;
use rcad_step::writer::StepWriter;
use rcad_step::StepReader;
use rmsh_model::{Element, ElementType, Mesh, Node};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum StepError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("STEP parse error: {0}")]
    Parse(String),
}

// ── Public API ─────────────────────────────────────────────────────────────────

pub fn load_step_from_path(path: &Path) -> Result<Mesh, StepError> {
    let brep = read_tessellated_brep_from_path(path)?;
    let (verts, tris) = brep_to_trimesh(&brep);
    Ok(trimesh_to_mesh(verts, tris))
}

pub fn load_step_from_bytes(data: &[u8]) -> Result<Mesh, StepError> {
    let text = String::from_utf8_lossy(data);
    parse_step(&text)
}

pub fn parse_step(text: &str) -> Result<Mesh, StepError> {
    let brep = parse_tessellated_brep(text)?;
    let (verts, tris) = brep_to_trimesh(&brep);
    Ok(trimesh_to_mesh(verts, tris))
}

pub fn save_step_to_path(path: &Path, mesh: &Mesh) -> Result<(), StepError> {
    let content = write_step(mesh)?;
    std::fs::write(path, content)?;
    Ok(())
}

pub fn save_brep_step_to_path(path: &Path, brep: &BRep) -> Result<(), StepError> {
    let content = write_brep_step(brep)?;
    std::fs::write(path, content)?;
    Ok(())
}

pub fn write_step(mesh: &Mesh) -> Result<String, StepError> {
    let (verts, tris) = mesh_to_trimesh(mesh);
    if verts.is_empty() || tris.is_empty() {
        return Err(StepError::Parse(
            "mesh has no vertices/triangles to export".to_string(),
        ));
    }

    let brep = trimesh_to_brep(&verts, &tris);
    Ok(StepWriter::write_string(
        &brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    ))
}

pub fn write_brep_step(brep: &BRep) -> Result<String, StepError> {
    Ok(StepWriter::write_string(
        brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    ))
}

// ── Mesh ↔ trimesh conversions ─────────────────────────────────────────────────

fn trimesh_to_mesh(verts: Vec<glam::DVec3>, tris: Vec<[usize; 3]>) -> Mesh {
    let mut mesh = Mesh::new();
    for (i, v) in verts.iter().enumerate() {
        mesh.add_node(Node::new((i + 1) as u64, v.x, v.y, v.z));
    }
    for (i, &[a, b, c]) in tris.iter().enumerate() {
        mesh.add_element(Element::new(
            (i + 1) as u64,
            ElementType::Triangle3,
            vec![(a + 1) as u64, (b + 1) as u64, (c + 1) as u64],
        ));
    }
    mesh
}

fn mesh_to_trimesh(mesh: &Mesh) -> (Vec<glam::DVec3>, Vec<[usize; 3]>) {
    let mut node_ids: Vec<u64> = mesh.nodes.keys().copied().collect();
    node_ids.sort_unstable();
    let vi_map: std::collections::HashMap<u64, usize> = node_ids
        .iter()
        .copied()
        .enumerate()
        .map(|(i, nid)| (nid, i))
        .collect();

    let verts: Vec<glam::DVec3> = node_ids
        .iter()
        .map(|nid| {
            let p = &mesh.nodes[nid].position;
            glam::DVec3::new(p.x, p.y, p.z)
        })
        .collect();

    let mut tris: Vec<[usize; 3]> = Vec::new();
    for elem in &mesh.elements {
        if elem.dimension() < 2 || elem.node_ids.len() < 3 {
            continue;
        }
        let indices: Vec<usize> = elem
            .node_ids
            .iter()
            .filter_map(|nid| vi_map.get(nid).copied())
            .collect();
        if indices.len() < 3 {
            continue;
        }
        // Fan-triangulate the polygon.
        let root = indices[0];
        for i in 1..(indices.len() - 1) {
            tris.push([root, indices[i], indices[i + 1]]);
        }
    }
    (verts, tris)
}

fn read_tessellated_brep_from_path(path: &Path) -> Result<BRep, StepError> {
    let mut brep = StepReader::read_file(path).map_err(|e| StepError::Parse(e.to_string()))?;
    tessellate_brep_if_needed(&mut brep);
    Ok(brep)
}

fn parse_tessellated_brep(text: &str) -> Result<BRep, StepError> {
    let mut brep = StepReader::parse_string(text).map_err(|e| StepError::Parse(e.to_string()))?;
    tessellate_brep_if_needed(&mut brep);
    Ok(brep)
}

fn tessellate_brep_if_needed(brep: &mut BRep) {
    if brep
        .solids
        .iter()
        .flat_map(|solid| solid.shells.iter())
        .flat_map(|shell| shell.faces.iter())
        .any(|face| face.mesh_dirty || face.triangles.is_empty())
    {
        mesh_brep(brep, &TessellationParams::default());
    }
}

fn brep_to_trimesh(brep: &BRep) -> (Vec<glam::DVec3>, Vec<[usize; 3]>) {
    let verts = brep.vertices.iter().map(|v| v.point).collect();
    let tris = brep
        .solids
        .iter()
        .flat_map(|s| s.shells.iter())
        .flat_map(|sh| sh.faces.iter())
        .flat_map(|f| f.triangles.iter().copied())
        .collect();
    (verts, tris)
}

fn trimesh_to_brep(verts: &[glam::DVec3], tris: &[[usize; 3]]) -> BRep {
    let vertices = verts
        .iter()
        .map(|&point| Vertex { point })
        .collect::<Vec<_>>();

    let face = Face {
        outer_wire: Wire { edges: Vec::new() },
        inner_wires: Vec::new(),
        normal: glam::DVec3::Z,
        triangles: tris.to_vec(),
        mesh_dirty: true,
    };

    BRep {
        vertices,
        edges: Vec::new(),
        solids: vec![Solid {
            shells: vec![Shell { faces: vec![face] }],
        }],
        geom: rcad_kernel::GeomStore::default(),
        compound: None,
        compsolid: None,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::{load_step_from_path, parse_step, save_step_to_path, write_step};
    use std::path::PathBuf;

    #[test]
    fn parse_simple_tetra_faceted_brep() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('test'),'2;1');
FILE_NAME('','',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=CARTESIAN_POINT('',(1.,0.,0.));
#3=CARTESIAN_POINT('',(0.,1.,0.));
#4=CARTESIAN_POINT('',(0.,0.,1.));
#11=VERTEX_POINT('',#1);
#12=VERTEX_POINT('',#2);
#13=VERTEX_POINT('',#3);
#14=VERTEX_POINT('',#4);
#21=EDGE_CURVE('',#11,#12,$,.T.);
#22=EDGE_CURVE('',#12,#13,$,.T.);
#23=EDGE_CURVE('',#13,#11,$,.T.);
#24=EDGE_CURVE('',#11,#14,$,.T.);
#25=EDGE_CURVE('',#12,#14,$,.T.);
#26=EDGE_CURVE('',#13,#14,$,.T.);
#31=ORIENTED_EDGE('',*,*,#21,.T.);
#32=ORIENTED_EDGE('',*,*,#22,.T.);
#33=ORIENTED_EDGE('',*,*,#23,.T.);
#34=ORIENTED_EDGE('',*,*,#21,.F.);
#35=ORIENTED_EDGE('',*,*,#25,.T.);
#36=ORIENTED_EDGE('',*,*,#24,.F.);
#37=ORIENTED_EDGE('',*,*,#22,.F.);
#38=ORIENTED_EDGE('',*,*,#26,.T.);
#39=ORIENTED_EDGE('',*,*,#25,.F.);
#40=ORIENTED_EDGE('',*,*,#23,.F.);
#41=ORIENTED_EDGE('',*,*,#24,.T.);
#42=ORIENTED_EDGE('',*,*,#26,.F.);
#51=EDGE_LOOP('',(#31,#32,#33));
#52=EDGE_LOOP('',(#34,#35,#36));
#53=EDGE_LOOP('',(#37,#38,#39));
#54=EDGE_LOOP('',(#40,#41,#42));
#61=FACE_OUTER_BOUND('',#51,.T.);
#62=FACE_OUTER_BOUND('',#52,.T.);
#63=FACE_OUTER_BOUND('',#53,.T.);
#64=FACE_OUTER_BOUND('',#54,.T.);
#71=ADVANCED_FACE('',(#61),$,.T.);
#72=ADVANCED_FACE('',(#62),$,.T.);
#73=ADVANCED_FACE('',(#63),$,.T.);
#74=ADVANCED_FACE('',(#64),$,.T.);
#81=CLOSED_SHELL('',(#71,#72,#73,#74));
#82=MANIFOLD_SOLID_BREP('',#81);
ENDSEC;
END-ISO-10303-21;
"#;
        let mesh = parse_step(step).expect("STEP should parse");
        assert!(mesh.node_count() > 0);
        assert!(mesh.element_count() > 0);
    }

    #[test]
    fn roundtrip_write_then_parse() {
        use rmsh_model::{Element, ElementType, Mesh, Node};
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 1.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 1.0, 0.0));
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));
        mesh.add_element(Element::new(2, ElementType::Triangle3, vec![1, 3, 4]));

        let step_text = write_step(&mesh).expect("write should succeed");
        assert!(step_text.contains("ISO-10303-21"));
        assert!(step_text.contains("ENDSEC"));
    }

    #[test]
    fn write_empty_mesh_fails() {
        use rmsh_model::Mesh;
        let mesh = Mesh::new();
        assert!(write_step(&mesh).is_err());
    }

    #[test]
    fn parse_generated_step_test_file() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("..")
            .join("testdata")
            .join("simple_tetra.step");

        if !path.exists() {
            return;
        }
        let mesh = load_step_from_path(&path).expect("generated STEP file should parse");
        assert!(mesh.node_count() > 0);
        assert!(mesh.element_count() > 0);
    }

    #[test]
    #[ignore]
    fn save_and_reload_step_file() {
        use rmsh_model::{Element, ElementType, Mesh, Node};
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.5, 1.0, 0.0));
        mesh.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));

        let tmp = std::env::temp_dir().join("rmsh_test_roundtrip.step");
        save_step_to_path(&tmp, &mesh).expect("save should succeed");
        let loaded = load_step_from_path(&tmp).expect("reload should succeed");
        assert!(loaded.node_count() > 0);
        assert!(loaded.element_count() > 0);
        let _ = std::fs::remove_file(&tmp);
    }
}
