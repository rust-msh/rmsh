//! STEP (.step / .stp) import and export via `rcad-step`.
//!
//! The heavy lifting (BRep parsing, triangle extraction, STEP serialisation)
//! lives in `rcad-step`.  This module is a thin adapter that converts between
//! `rcad-step`'s raw `(Vec<DVec3>, Vec<[usize;3]>)` representation and
//! `rmsh_model::Mesh`.

use std::path::Path;

use rcad_algorithms::{TessellationParams, mesh_brep};
use rcad_kernel::geom::{Curve3, Line3};
use rcad_kernel::appearance::{Color, StepColor};
use rcad_kernel::{BRep, Face, Shell, Solid, Vertex, Wire};
use rcad_step::ExportSelection;
use rcad_step::writer::{StepHeader, StepProtocol, StepWriteOptions, StepWriter};
use rcad_step::StepReader;
use rmsh_model::{Element, ElementType, Mesh, Node};
use thiserror::Error;

fn strict_export_selection(brep: &BRep) -> Option<(Vec<usize>, Vec<usize>)> {
    let mut selected_faces = Vec::new();
    let mut face_edge_used = vec![false; brep.edges.len()];

    let mut face_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                selected_faces.push(face_idx);
                face_idx += 1;

                for we in &face.outer_wire.edges {
                    if we.idx < face_edge_used.len() {
                        face_edge_used[we.idx] = true;
                    }
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        if we.idx < face_edge_used.len() {
                            face_edge_used[we.idx] = true;
                        }
                    }
                }
            }
        }
    }

    let mut selected_edges = Vec::new();
    let edge_curve = &brep.geom.edge_curve;
    for (edge_idx, used_by_face) in face_edge_used.iter().copied().enumerate() {
        if used_by_face {
            continue;
        }
        let has_curve3 = edge_curve
            .get(edge_idx)
            .and_then(|v| *v)
            .is_some();
        if has_curve3 {
            selected_edges.push(edge_idx);
        }
    }

    if selected_faces.is_empty() && selected_edges.is_empty() {
        None
    } else {
        Some((selected_faces, selected_edges))
    }
}

fn normalize_for_strict_step_export(brep: &BRep) -> BRep {
    let mut out = brep.clone();

    if out.geom.edge_curve.len() < out.edges.len() {
        out.geom.edge_curve.resize(out.edges.len(), None);
    }
    if out.geom.edge_curve_range.len() < out.edges.len() {
        out.geom.edge_curve_range.resize(out.edges.len(), None);
    }

    let mut referenced = vec![false; out.edges.len()];
    for solid in &out.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    if we.idx < referenced.len() {
                        referenced[we.idx] = true;
                    }
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        if we.idx < referenced.len() {
                            referenced[we.idx] = true;
                        }
                    }
                }
            }
        }
    }

    for (edge_idx, edge) in out.edges.iter().enumerate() {
        if !referenced[edge_idx] || out.geom.edge_curve[edge_idx].is_some() {
            continue;
        }
        let Some(ps) = out.vertices.get(edge.start).map(|v| v.point) else {
            continue;
        };
        let Some(pe) = out.vertices.get(edge.end).map(|v| v.point) else {
            continue;
        };
        let d = pe - ps;
        let len = d.length();
        if !len.is_finite() || len <= 1e-12 {
            continue;
        }
        let dir = d / len;
        let cid = out.geom.curves.len();
        out.geom.curves.push(Curve3::Line(Line3 {
            origin: ps,
            direction: dir,
        }));
        out.geom.edge_curve[edge_idx] = Some(cid);
        out.geom.edge_curve_range[edge_idx] = Some([0.0, len]);
    }

    let edge_curve = out.geom.edge_curve.clone();
    for solid in &mut out.solids {
        for shell in &mut solid.shells {
            for face in &mut shell.faces {
                sanitize_wire_for_strict_export(&mut face.outer_wire, &out.edges, &edge_curve);
                for wire in &mut face.inner_wires {
                    sanitize_wire_for_strict_export(wire, &out.edges, &edge_curve);
                }
            }
        }
    }

    out
}

fn sanitize_wire_for_strict_export(
    wire: &mut Wire,
    edges: &[rcad_kernel::Edge],
    edge_curve: &[Option<usize>],
) {
    let mut filtered = Vec::with_capacity(wire.edges.len());
    for we in &wire.edges {
        if we.idx >= edges.len() {
            continue;
        }
        let has_curve = edge_curve.get(we.idx).and_then(|v| *v).is_some();
        if !has_curve {
            let e = &edges[we.idx];
            if e.start == e.end {
                continue;
            }
        }
        if filtered
            .last()
            .is_some_and(|prev: &rcad_kernel::topology::WireEdge| prev.idx == we.idx && prev.forward == we.forward)
        {
            continue;
        }
        filtered.push(*we);
    }
    wire.edges = filtered;
}

#[derive(Error, Debug)]
pub enum StepError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("STEP parse error: {0}")]
    Parse(String),
}

#[derive(Debug, Clone)]
pub struct BrepStepWriteOptions {
    pub protocol: StepProtocol,
    pub solid_color: Option<Color>,
    pub header: Option<StepHeader>,
    pub gmsh_strict: bool,
}

impl Default for BrepStepWriteOptions {
    fn default() -> Self {
        Self {
            protocol: StepProtocol::Ap242,
            solid_color: None,
            header: None,
            gmsh_strict: false,
        }
    }
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

pub fn save_brep_step_to_path_with_options(
    path: &Path,
    brep: &BRep,
    options: &BrepStepWriteOptions,
) -> Result<(), StepError> {
    let content = write_brep_step_with_options(brep, options)?;
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
    write_brep_step_with_options(
        brep,
        &BrepStepWriteOptions {
            protocol: StepProtocol::Ap214,
            ..Default::default()
        },
    )
}

pub fn write_brep_step_with_options(
    brep: &BRep,
    options: &BrepStepWriteOptions,
) -> Result<String, StepError> {
    let colors = options.solid_color.map(|c| StepColor {
        solid_color: Some(c),
        face_colors: Vec::new(),
    });
    let step_options = StepWriteOptions {
        protocol: options.protocol,
        colors,
        properties: Vec::new(),
        ap242_metadata: None,
        header: options.header.clone().unwrap_or_default(),
        export_standalone_wire_overlay: true,
    };

    Ok(StepWriter::write_string_with_options(
        brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
        &step_options,
    ))
}

// ── Mesh �?trimesh conversions ─────────────────────────────────────────────────

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
        sample_point: None,
        surface_idx: None,
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
    use super::{
        load_step_from_path, normalize_for_strict_step_export, parse_step, save_step_to_path,
        strict_export_selection, write_brep_step,
        write_brep_step_with_options, write_step, BrepStepWriteOptions,
    };
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

    #[test]
    fn write_brep_default_protocol_is_ap214() {
        let verts = vec![
            glam::DVec3::new(0.0, 0.0, 0.0),
            glam::DVec3::new(1.0, 0.0, 0.0),
            glam::DVec3::new(0.0, 1.0, 0.0),
        ];
        let tris = vec![[0, 1, 2]];
        let brep = super::trimesh_to_brep(&verts, &tris);

        let step = write_brep_step(&brep).expect("default STEP export should succeed");
        assert!(step.contains("FILE_SCHEMA"));
        assert!(step.to_ascii_lowercase().contains("214"));
    }

    #[test]
    fn write_brep_ap242_with_color_emits_style_chain() {
        let verts = vec![
            glam::DVec3::new(0.0, 0.0, 0.0),
            glam::DVec3::new(1.0, 0.0, 0.0),
            glam::DVec3::new(0.0, 1.0, 0.0),
        ];
        let tris = vec![[0, 1, 2]];
        let brep = super::trimesh_to_brep(&verts, &tris);
        let options = BrepStepWriteOptions {
            protocol: rcad_step::StepProtocol::Ap242,
            solid_color: Some(rcad_kernel::appearance::Color::from_rgb8(30, 144, 255)),
            header: None,
            gmsh_strict: false,
        };

        let step = write_brep_step_with_options(&brep, &options)
            .expect("AP242 + color STEP export should succeed");

        assert!(step.contains("AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF"));
        assert!(step.contains("COLOUR_RGB"));
        assert!(step.contains("PRESENTATION_STYLE_ASSIGNMENT"));
        assert!(step.contains("STYLED_ITEM"));
    }

    #[test]
    fn strict_selection_keeps_face_edges_and_curved_standalone_edges_only() {
        use rcad_kernel::geom::{Curve3, Line3};
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
        use rcad_kernel::BRep;

        let mut brep = BRep::new();
        brep.vertices = vec![
            Vertex { point: glam::DVec3::new(0.0, 0.0, 0.0) },
            Vertex { point: glam::DVec3::new(1.0, 0.0, 0.0) },
            Vertex { point: glam::DVec3::new(0.0, 1.0, 0.0) },
            Vertex { point: glam::DVec3::new(2.0, 0.0, 0.0) },
        ];
        brep.edges = vec![
            Edge { start: 0, end: 1 }, // used by face
            Edge { start: 1, end: 2 }, // orphan, no curve
            Edge { start: 2, end: 3 }, // orphan, has curve
        ];
        brep.solids = vec![Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![WireEdge::fwd(0)],
                    },
                    inner_wires: Vec::new(),
                    normal: glam::DVec3::Z,
                    triangles: vec![[0, 1, 2]],
                    mesh_dirty: false,
                    sample_point: None,
                    surface_idx: None,
                    sample_point: None,
                    surface_idx: None,
                }],
            }],
        }];
        brep.geom.edge_curve = vec![
            Some(0),
            None,
            Some(1),
        ];
        brep.geom.curves = vec![
            Curve3::Line(Line3 {
                origin: glam::DVec3::new(0.0, 0.0, 0.0),
                direction: glam::DVec3::X,
            }),
            Curve3::Line(Line3 {
                origin: glam::DVec3::new(0.0, 1.0, 0.0),
                direction: glam::DVec3::X,
            }),
        ];

        let (faces, edges) = strict_export_selection(&brep).expect("selection should exist");
        assert_eq!(faces, vec![0]);
        assert_eq!(edges, vec![2]);
    }

    #[test]
    fn strict_normalize_fills_missing_line_curve_on_face_edge() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
        use rcad_kernel::BRep;

        let mut brep = BRep::new();
        brep.vertices = vec![
            Vertex { point: glam::DVec3::new(0.0, 0.0, 0.0) },
            Vertex { point: glam::DVec3::new(1.0, 0.0, 0.0) },
            Vertex { point: glam::DVec3::new(0.0, 1.0, 0.0) },
        ];
        brep.edges = vec![
            Edge { start: 0, end: 1 },
            Edge { start: 1, end: 2 },
        ];
        brep.solids = vec![Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
                    },
                    inner_wires: Vec::new(),
                    normal: glam::DVec3::Z,
                    triangles: vec![[0, 1, 2]],
                    mesh_dirty: false,
                    sample_point: None,
                    surface_idx: None,
                    sample_point: None,
                    surface_idx: None,
                }],
            }],
        }];
        brep.geom.edge_curve = vec![None, None];
        brep.geom.edge_curve_range = vec![None, None];

        let normalized = normalize_for_strict_step_export(&brep);
        assert!(normalized.geom.edge_curve[0].is_some());
        assert!(normalized.geom.edge_curve[1].is_some());
        assert!(normalized.geom.edge_curve_range[0].is_some());
        assert!(normalized.geom.edge_curve_range[1].is_some());
    }

    #[test]
    fn strict_normalize_drops_degenerate_no_curve_wire_edges() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
        use rcad_kernel::BRep;

        let mut brep = BRep::new();
        brep.vertices = vec![
            Vertex { point: glam::DVec3::new(0.0, 0.0, 0.0) },
            Vertex { point: glam::DVec3::new(1.0, 0.0, 0.0) },
        ];
        brep.edges = vec![
            Edge { start: 0, end: 0 },
            Edge { start: 0, end: 1 },
        ];
        brep.solids = vec![Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(1)],
                    },
                    inner_wires: Vec::new(),
                    normal: glam::DVec3::Z,
                    triangles: vec![[0, 1, 1]],
                    mesh_dirty: false,
                    sample_point: None,
                    surface_idx: None,
                    sample_point: None,
                    surface_idx: None,
                }],
            }],
        }];
        brep.geom.edge_curve = vec![None, None];
        brep.geom.edge_curve_range = vec![None, None];

        let normalized = normalize_for_strict_step_export(&brep);
        let outer = &normalized.solids[0].shells[0].faces[0].outer_wire.edges;
        assert!(!outer.iter().any(|we| we.idx == 0));
        assert_eq!(outer.iter().filter(|we| we.idx == 1).count(), 1);
    }
}
