//! Minimal IGES mesh exchange (Type 106 copious-data polyline records).
//!
//! This module provides a pragmatic IGES bridge:
//! - Export: each triangle is written as a closed 3D polyline entity.
//! - Import: type-106 polyline entities are parsed and fan-triangulated.
//!
//! Scope: mesh-level interoperability, not full analytic B-Rep IGES support.

use std::io;
use std::io::Write;
use std::path::Path;

use glam::DVec3;
use rcad_kernel::topology::{Face, Shell, Solid, Vertex, Wire};
use rcad_kernel::BRep;

#[derive(Debug, Clone)]
pub enum IgesError {
    Io(String),
    InvalidFormat(String),
    EmptyResult(String),
}

impl std::fmt::Display for IgesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(msg) => write!(f, "I/O error: {msg}"),
            Self::InvalidFormat(msg) => write!(f, "invalid IGES format: {msg}"),
            Self::EmptyResult(msg) => write!(f, "IGES parse produced empty result: {msg}"),
        }
    }
}

impl std::error::Error for IgesError {}

pub struct IgesWriter;

impl IgesWriter {
    pub fn write_string(brep: &BRep) -> String {
        let mut out = Vec::new();
        let _ = Self::write_to(brep, &mut out);
        String::from_utf8_lossy(&out).into_owned()
    }

    pub fn write_file<P: AsRef<Path>>(brep: &BRep, path: P) -> Result<usize, io::Error> {
        let mut file = std::fs::File::create(path)?;
        Self::write_to(brep, &mut file)
    }

    pub fn write_to(brep: &BRep, writer: &mut impl Write) -> Result<usize, io::Error> {
        let mut entities = Vec::new();
        let mut tri_count = 0usize;

        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for &[i, j, k] in &face.triangles {
                        if i >= brep.vertices.len() || j >= brep.vertices.len() || k >= brep.vertices.len() {
                            continue;
                        }
                        let a = brep.vertices[i].point;
                        let b = brep.vertices[j].point;
                        let c = brep.vertices[k].point;
                        entities.push(format!(
                            "106,2,4,{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9};",
                            a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, a.x, a.y, a.z
                        ));
                        tri_count += 1;
                    }
                }
            }
        }

        let mut lines: Vec<String> = Vec::new();
        let mut s_count = 0usize;
        let mut g_count = 0usize;

        s_count += 1;
        lines.push(section_line(
            "RCAD IGES mesh export (Type 106 polyline triangles)",
            'S',
            s_count,
        ));

        g_count += 1;
        lines.push(section_line("1H,,1H;", 'G', g_count));

        let mut p_lines: Vec<(String, usize)> = Vec::new();
        let mut p_meta: Vec<(usize, usize)> = Vec::with_capacity(entities.len());
        let mut p_seq = 1usize;
        for (entity_idx, payload) in entities.iter().enumerate() {
            let start = p_seq;
            let chunks = chunk_ascii(payload, 64);
            for ch in chunks {
                // Column 1..64 parameter text, 65..72 directory pointer.
                let data = format!("{:<64}{:>8}", ch, entity_idx * 2 + 1);
                p_lines.push((data, p_seq));
                p_seq += 1;
            }
            p_meta.push((start, p_seq - start));
        }

        let mut d_count = 0usize;
        for (idx, (p_start, p_len)) in p_meta.iter().copied().enumerate() {
            let d_seq1 = d_count + 1;
            let d_seq2 = d_count + 2;
            let line1 = d_line([106, p_start as i32, 0, 1, 0, 0, 0, 0, d_seq2 as i32]);
            let line2 = d_line([106, 0, 0, p_len as i32, 0, 0, 0, 0, 0]);
            lines.push(section_line(&line1, 'D', d_seq1));
            lines.push(section_line(&line2, 'D', d_seq2));
            d_count += 2;
            let _ = idx;
        }

        for (data, seq) in p_lines {
            lines.push(section_line(&data, 'P', seq));
        }
        let p_count = p_seq.saturating_sub(1);

        let term = format!("S{:>7}G{:>7}D{:>7}P{:>7}", s_count, g_count, d_count, p_count);
        lines.push(section_line(&term, 'T', 1));

        for line in lines {
            writer.write_all(line.as_bytes())?;
            writer.write_all(b"\n")?;
        }

        Ok(tri_count)
    }
}

pub struct IgesReader;

impl IgesReader {
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<BRep, IgesError> {
        let content = std::fs::read_to_string(path).map_err(|e| IgesError::Io(e.to_string()))?;
        Self::parse_string(&content)
    }

    pub fn parse_string(content: &str) -> Result<BRep, IgesError> {
        let mut parameter_blob = String::new();

        for line in content.lines() {
            if is_section_line(line, 'P') {
                let head = first_n_chars(line, 64);
                let text = head.trim_end();
                parameter_blob.push_str(text);
            }
        }

        if parameter_blob.trim().is_empty() {
            parameter_blob = content.to_string();
        }

        let mut brep = BRep::new();
        let mut faces = Vec::new();

        for rec in parameter_blob.split(';') {
            let rec = rec.trim();
            if rec.is_empty() {
                continue;
            }
            let vals: Vec<&str> = rec
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .collect();
            if vals.len() < 3 {
                continue;
            }

            let Ok(entity_type) = vals[0].parse::<i32>() else {
                continue;
            };
            if entity_type != 106 {
                continue;
            }

            let n_points = vals[2].parse::<usize>().map_err(|_| {
                IgesError::InvalidFormat(format!("invalid type-106 point count in record: {rec}"))
            })?;
            if vals.len() < 3 + n_points * 3 {
                return Err(IgesError::InvalidFormat(format!(
                    "type-106 record has insufficient coordinates: {rec}"
                )));
            }

            let mut points = Vec::with_capacity(n_points);
            for i in 0..n_points {
                let x = parse_iges_float(vals[3 + i * 3])?;
                let y = parse_iges_float(vals[4 + i * 3])?;
                let z = parse_iges_float(vals[5 + i * 3])?;
                points.push(DVec3::new(x, y, z));
            }

            if points.len() >= 2 && points[0].distance(points[points.len() - 1]) <= 1e-9 {
                points.pop();
            }
            if points.len() < 3 {
                continue;
            }

            let base = brep.vertices.len();
            for p in &points {
                brep.vertices.push(Vertex { point: *p });
            }

            let mut tris = Vec::with_capacity(points.len().saturating_sub(2));
            for i in 1..points.len() - 1 {
                tris.push([base, base + i, base + i + 1]);
            }
            if tris.is_empty() {
                continue;
            }

            let normal = compute_triangle_normal(points[0], points[1], points[2]);
            faces.push(Face {
                outer_wire: Wire { edges: vec![] },
                inner_wires: vec![],
                normal,
                triangles: tris,
                mesh_dirty: false,
            });
        }

        if faces.is_empty() {
            return Err(IgesError::EmptyResult(
                "no type-106 polyline entities were found".into(),
            ));
        }

        brep.solids.push(Solid {
            shells: vec![Shell { faces }],
        });
        Ok(brep)
    }
}

fn d_line(fields: [i32; 9]) -> String {
    let mut out = String::with_capacity(72);
    for value in fields {
        out.push_str(&format!("{:>8}", value));
    }
    out
}

fn section_line(data: &str, section: char, seq: usize) -> String {
    format!("{:<72}{}{:>7}", truncate_ascii(data, 72), section, seq)
}

fn truncate_ascii(text: &str, width: usize) -> String {
    if text.len() <= width {
        return text.to_string();
    }
    text.chars().take(width).collect()
}

fn chunk_ascii(text: &str, chunk: usize) -> Vec<String> {
    if text.is_empty() {
        return vec![String::new()];
    }
    let mut out = Vec::new();
    let mut start = 0usize;
    while start < text.len() {
        let end = (start + chunk).min(text.len());
        out.push(text[start..end].to_string());
        start = end;
    }
    out
}

fn is_section_line(line: &str, section: char) -> bool {
    line.chars().nth(72) == Some(section)
}

fn first_n_chars(line: &str, n: usize) -> String {
    line.chars().take(n).collect()
}

fn parse_iges_float(raw: &str) -> Result<f64, IgesError> {
    let normalized = raw.replace('D', "E").replace('d', "e");
    normalized
        .parse::<f64>()
        .map_err(|_| IgesError::InvalidFormat(format!("invalid float '{raw}'")))
}

fn compute_triangle_normal(a: DVec3, b: DVec3, c: DVec3) -> DVec3 {
    let n = (b - a).cross(c - a);
    if n.length_squared() < 1.0e-24 {
        DVec3::Z
    } else {
        n.normalize()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_tri_brep() -> BRep {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.solids.push(Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire { edges: vec![] },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![[0, 1, 2]],
                    mesh_dirty: false,
                }],
            }],
        });
        brep
    }

    #[test]
    fn iges_round_trip_mesh() {
        let src = make_tri_brep();
        let text = IgesWriter::write_string(&src);
        assert!(text.contains("106,2,4"));

        let parsed = IgesReader::parse_string(&text).expect("IGES parse should succeed");
        assert_eq!(parsed.vertices.len(), 3);
        assert_eq!(parsed.solids.len(), 1);
        assert_eq!(parsed.solids[0].shells[0].faces.len(), 1);
        assert_eq!(parsed.solids[0].shells[0].faces[0].triangles.len(), 1);
    }

    #[test]
    fn iges_parse_invalid_106_returns_error() {
        let bad = "106,2,4,0,0,0,1,0,0;";
        let err = IgesReader::parse_string(bad).expect_err("expected parse error");
        match err {
            IgesError::InvalidFormat(_) => {}
            _ => panic!("unexpected error variant"),
        }
    }
}
