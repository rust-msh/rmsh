//! OBJ mesh exporter for BRep solids.
//!
//! Writes the triangulated faces of a `BRep` as a Wavefront OBJ file.
//! Each face must already be triangulated (`Face.triangles` non-empty); faces
//! without triangles are skipped.
//!
//! Analogous to the mesh-export path in OCCT `RWMesh_FaceMeshComp`.

use std::io::{self, Write};
use std::path::Path;

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::topology::{Face, Shell, Solid, Vertex, Wire};

/// Errors that can occur when reading/parsing OBJ files.
#[derive(Debug, Clone)]
pub enum ObjError {
    Io(String),
    InvalidFormat(String),
    EmptyResult(String),
}

impl std::fmt::Display for ObjError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(msg) => write!(f, "I/O error: {msg}"),
            Self::InvalidFormat(msg) => write!(f, "invalid OBJ format: {msg}"),
            Self::EmptyResult(msg) => write!(f, "OBJ parse produced empty result: {msg}"),
        }
    }
}

impl std::error::Error for ObjError {}

/// Wavefront OBJ reader (mesh-only).
pub struct ObjReader;

impl ObjReader {
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<BRep, ObjError> {
        let content = std::fs::read_to_string(path).map_err(|e| ObjError::Io(e.to_string()))?;
        Self::parse_string(&content)
    }

    pub fn parse_string(content: &str) -> Result<BRep, ObjError> {
        let mut positions: Vec<DVec3> = Vec::new();
        let mut triangles: Vec<[usize; 3]> = Vec::new();

        for (line_idx, raw) in content.lines().enumerate() {
            let line_no = line_idx + 1;
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let mut parts = line.split_whitespace();
            let Some(kind) = parts.next() else {
                continue;
            };

            match kind {
                "v" => {
                    let x = parse_f64(parts.next(), line_no, "x")?;
                    let y = parse_f64(parts.next(), line_no, "y")?;
                    let z = parse_f64(parts.next(), line_no, "z")?;
                    positions.push(DVec3::new(x, y, z));
                }
                "f" => {
                    let refs: Vec<String> = parts.map(ToString::to_string).collect();
                    if refs.len() < 3 {
                        return Err(ObjError::InvalidFormat(format!(
                            "line {line_no}: face must have at least 3 vertices"
                        )));
                    }

                    let mut polygon = Vec::with_capacity(refs.len());
                    for rf in refs {
                        polygon.push(parse_face_index(&rf, positions.len(), line_no)?);
                    }

                    for i in 1..polygon.len() - 1 {
                        triangles.push([polygon[0], polygon[i], polygon[i + 1]]);
                    }
                }
                _ => {
                    // Ignore unsupported records (vn, vt, g, o, ...).
                }
            }
        }

        if positions.is_empty() {
            return Err(ObjError::EmptyResult("no vertices found".into()));
        }
        if triangles.is_empty() {
            return Err(ObjError::EmptyResult("no faces found".into()));
        }

        let mut brep = BRep::new();
        brep.vertices = positions.into_iter().map(|point| Vertex { point }).collect();

        let mut faces = Vec::with_capacity(triangles.len());
        for tri in triangles {
            let normal = compute_triangle_normal(
                brep.vertices[tri[0]].point,
                brep.vertices[tri[1]].point,
                brep.vertices[tri[2]].point,
            );
            faces.push(Face {
                outer_wire: Wire { edges: vec![] },
                inner_wires: vec![],
                normal,
                triangles: vec![tri],
                mesh_dirty: false,
            });
        }

        brep.solids.push(Solid {
            shells: vec![Shell { faces }],
        });
        Ok(brep)
    }
}

/// Wavefront OBJ writer helpers.
pub struct ObjWriter;

impl ObjWriter {
    pub fn write_string(brep: &BRep) -> String {
        let mut out = Vec::new();
        let _ = write_obj(brep, &mut out);
        String::from_utf8_lossy(&out).into_owned()
    }

    pub fn write_file<P: AsRef<Path>>(brep: &BRep, path: P) -> Result<usize, io::Error> {
        let mut file = std::fs::File::create(path)?;
        write_obj(brep, &mut file)
    }
}

fn parse_f64(raw: Option<&str>, line_no: usize, field: &str) -> Result<f64, ObjError> {
    let Some(text) = raw else {
        return Err(ObjError::InvalidFormat(format!(
            "line {line_no}: missing vertex {field}"
        )));
    };
    text.parse::<f64>().map_err(|_| {
        ObjError::InvalidFormat(format!("line {line_no}: invalid float '{text}' for {field}"))
    })
}

fn parse_face_index(raw: &str, vertex_count: usize, line_no: usize) -> Result<usize, ObjError> {
    let Some(head) = raw.split('/').next() else {
        return Err(ObjError::InvalidFormat(format!(
            "line {line_no}: invalid face token '{raw}'"
        )));
    };

    let index = head.parse::<isize>().map_err(|_| {
        ObjError::InvalidFormat(format!("line {line_no}: invalid face index '{head}'"))
    })?;

    if index == 0 {
        return Err(ObjError::InvalidFormat(format!(
            "line {line_no}: OBJ index 0 is invalid"
        )));
    }

    let resolved = if index > 0 {
        index - 1
    } else {
        vertex_count as isize + index
    };

    if resolved < 0 || resolved as usize >= vertex_count {
        return Err(ObjError::InvalidFormat(format!(
            "line {line_no}: face index '{head}' out of range"
        )));
    }

    Ok(resolved as usize)
}

fn compute_triangle_normal(a: DVec3, b: DVec3, c: DVec3) -> DVec3 {
    let n = (b - a).cross(c - a);
    if n.length_squared() < 1.0e-24 {
        DVec3::Z
    } else {
        n.normalize()
    }
}

/// Write `brep` as Wavefront OBJ text to `writer`.
///
/// Vertex indices in the OBJ file are 1-based.  All triangles from all faces
/// of all solids are emitted.  Faces without pre-triangulated data are skipped.
///
/// Returns the number of triangles written.
pub fn write_obj(brep: &BRep, writer: &mut impl Write) -> io::Result<usize> {
    // Emit all unique vertices first.
    for v in &brep.vertices {
        writeln!(writer, "v {:.9} {:.9} {:.9}", v.point.x, v.point.y, v.point.z)?;
    }

    let mut total_tris = 0usize;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for &[i, j, k] in &face.triangles {
                    // OBJ is 1-based
                    writeln!(writer, "f {} {} {}", i + 1, j + 1, k + 1)?;
                    total_tris += 1;
                }
            }
        }
    }

    Ok(total_tris)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal triangulated BRep: one solid with a single face that has
    /// two triangles (a square split diagonally).
    fn make_triangulated_brep() -> BRep {
        let mut brep = BRep::new();
        // 4 vertices of a 1×1 square in XY plane
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 3

        let face = Face {
            outer_wire: Wire { edges: vec![] },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![[0, 1, 2], [0, 2, 3]],
            mesh_dirty: false,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });
        brep
    }

    #[test]
    fn write_obj_produces_correct_output() {
        let brep = make_triangulated_brep();
        let mut buf = Vec::new();
        let n = write_obj(&brep, &mut buf).expect("write_obj should succeed");
        let text = String::from_utf8(buf).unwrap();

        let f_lines: Vec<&str> = text.lines().filter(|l| l.starts_with('f')).collect();
        let v_lines: Vec<&str> = text.lines().filter(|l| l.starts_with('v')).collect();

        assert_eq!(n, 2, "should return 2 triangles");
        assert_eq!(f_lines.len(), 2, "should have 2 'f' lines");
        assert_eq!(v_lines.len(), 4, "should have 4 'v' lines");

        // OBJ indices are 1-based
        assert!(text.contains("f 1 2 3"), "first triangle should be f 1 2 3");
        assert!(text.contains("f 1 3 4"), "second triangle should be f 1 3 4");
    }

    #[test]
    fn write_obj_empty_brep() {
        let brep = BRep::new();
        let mut buf = Vec::new();
        let n = write_obj(&brep, &mut buf).expect("write_obj should handle empty BRep");
        assert_eq!(n, 0);
    }

    #[test]
    fn write_obj_face_without_triangles_is_skipped() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        let face = Face {
            outer_wire: Wire { edges: vec![] },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![], // no triangles
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });
        let mut buf = Vec::new();
        let n = write_obj(&brep, &mut buf).unwrap();
        assert_eq!(n, 0, "face without triangles should produce 0 triangles");
    }

    #[test]
    fn parse_obj_triangle() {
        let obj = "\
v 0 0 0
v 1 0 0
v 0 1 0
f 1 2 3
";
        let brep = ObjReader::parse_string(obj).expect("obj parse should succeed");
        assert_eq!(brep.vertices.len(), 3);
        assert_eq!(brep.solids.len(), 1);
        assert_eq!(brep.solids[0].shells[0].faces.len(), 1);
        assert_eq!(brep.solids[0].shells[0].faces[0].triangles, vec![[0, 1, 2]]);
    }

    #[test]
    fn parse_obj_negative_indices() {
        let obj = "\
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
f -4 -3 -2 -1
";
        let brep = ObjReader::parse_string(obj).expect("obj parse with negative indices");
        let faces = &brep.solids[0].shells[0].faces;
        assert_eq!(faces.len(), 2, "quad should fan-triangulate to 2 faces");
    }

    #[test]
    fn obj_round_trip() {
        let src = make_triangulated_brep();
        let text = ObjWriter::write_string(&src);
        let dst = ObjReader::parse_string(&text).expect("round-trip parse should succeed");
        let tri_count: usize = dst.solids[0].shells[0]
            .faces
            .iter()
            .map(|f| f.triangles.len())
            .sum();
        assert_eq!(tri_count, 2);
        assert_eq!(dst.vertices.len(), 4);
    }
}
