use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, Cursor, Write};
use std::path::Path;

use rmsh_model::{Element, ElementType, Mesh, Node};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum MshError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parse error at line {line}: {message}")]
    Parse { line: usize, message: String },
    #[error("Unsupported MSH format version: {0}")]
    UnsupportedVersion(String),
    #[error("Unsupported element type for MSH write: {0:?}")]
    UnsupportedElementType(ElementType),
    #[error("Element references missing node ID: {0}")]
    MissingNode(u64),
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum MshVersion {
    V2,
    V4,
}

pub fn load_msh_from_path(path: &Path) -> Result<Mesh, MshError> {
    let bytes = std::fs::read(path)?;
    load_msh_from_bytes(&bytes)
}

pub fn load_msh_from_bytes(data: &[u8]) -> Result<Mesh, MshError> {
    // Gmsh v2 binary files keep ASCII section markers but encode $Nodes/$Elements payload
    // in little-endian binary. Detect by checking the MeshFormat line near file start.
    let header = &data[..data.len().min(160)];
    let is_binary_v2 = contains_bytes(header, b"$MeshFormat")
        && (contains_bytes(header, b"2.2 1 ") || contains_bytes(header, b"2.0 1 "));
    let is_binary_v4 = contains_bytes(header, b"$MeshFormat") && contains_bytes(header, b"4.1 1 ");

    if is_binary_v2 {
        parse_msh_v2_binary(data)
    } else if is_binary_v4 {
        parse_msh_v4_binary(data)
    } else {
        parse_msh(Cursor::new(data))
    }
}

fn contains_bytes(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() {
        return true;
    }
    haystack.windows(needle.len()).any(|w| w == needle)
}

pub fn save_msh_v2_to_path(path: &Path, mesh: &Mesh) -> Result<(), MshError> {
    let mut file = File::create(path)?;
    write_msh_v2(&mut file, mesh)
}

pub fn save_msh_v4_to_path(path: &Path, mesh: &Mesh) -> Result<(), MshError> {
    let mut file = File::create(path)?;
    write_msh_v4(&mut file, mesh)
}

pub fn write_msh_v2<W: Write>(writer: &mut W, mesh: &Mesh) -> Result<(), MshError> {
    validate_mesh(mesh)?;

    let nodes = sorted_nodes(mesh);
    let elements = sorted_elements(mesh);

    writeln!(writer, "$MeshFormat")?;
    writeln!(writer, "2.2 0 8")?;
    writeln!(writer, "$EndMeshFormat")?;
    write_physical_names(writer, mesh)?;

    writeln!(writer, "$Nodes")?;
    writeln!(writer, "{}", nodes.len())?;
    for node in nodes {
        writeln!(
            writer,
            "{} {} {} {}",
            node.id, node.position.x, node.position.y, node.position.z
        )?;
    }
    writeln!(writer, "$EndNodes")?;

    writeln!(writer, "$Elements")?;
    writeln!(writer, "{}", elements.len())?;
    for element in elements {
        let etype = gmsh_type_id(element.etype)?;
        match element.physical_tag {
            Some(physical_tag) => {
                write!(writer, "{} {} 2 {} 0", element.id, etype, physical_tag)?;
            }
            None => {
                write!(writer, "{} {} 0", element.id, etype)?;
            }
        }
        for node_id in &element.node_ids {
            write!(writer, " {}", node_id)?;
        }
        writeln!(writer)?;
    }
    writeln!(writer, "$EndElements")?;

    // Gmsh writes $Periodic after $Elements; empty links emit nothing.
    crate::msh_periodic::write_periodic_v2(writer, &mesh.periodic)?;

    Ok(())
}

pub fn write_msh_v4<W: Write>(writer: &mut W, mesh: &Mesh) -> Result<(), MshError> {
    validate_mesh(mesh)?;

    let nodes = sorted_nodes(mesh);
    let elements = sorted_elements(mesh);
    let min_node_tag = nodes.first().map(|node| node.id).unwrap_or(0);
    let max_node_tag = nodes.last().map(|node| node.id).unwrap_or(0);
    let min_element_tag = elements.first().map(|element| element.id).unwrap_or(0);
    let max_element_tag = elements.last().map(|element| element.id).unwrap_or(0);
    let entity_dim = elements
        .iter()
        .map(|element| element.dimension() as i32)
        .max()
        .unwrap_or(0);

    writeln!(writer, "$MeshFormat")?;
    writeln!(writer, "4.1 0 8")?;
    writeln!(writer, "$EndMeshFormat")?;
    write_physical_names(writer, mesh)?;
    write_entities_v4(writer, &declared_entities(mesh, &elements))?;
    if let Some(info) = &mesh.partitioned {
        crate::msh_partitioned::write_partitioned_v4(writer, info)?;
    }

    writeln!(writer, "$Nodes")?;
    if nodes.is_empty() {
        writeln!(writer, "0 0 0 0")?;
    } else {
        writeln!(
            writer,
            "1 {} {} {}",
            nodes.len(),
            min_node_tag,
            max_node_tag
        )?;
        writeln!(writer, "{} 1 0 {}", entity_dim, nodes.len())?;
        for node in &nodes {
            writeln!(writer, "{}", node.id)?;
        }
        for node in &nodes {
            writeln!(
                writer,
                "{} {} {}",
                node.position.x, node.position.y, node.position.z
            )?;
        }
    }
    writeln!(writer, "$EndNodes")?;

    let mut blocks: BTreeMap<(u8, i32, i32), Vec<&rmsh_model::Element>> = BTreeMap::new();
    for element in &elements {
        let gmsh_type = gmsh_type_id(element.etype)?;
        let entity_tag = element.physical_tag.unwrap_or(1);
        blocks
            .entry((element.dimension(), entity_tag, gmsh_type))
            .or_default()
            .push(*element);
    }

    writeln!(writer, "$Elements")?;
    if elements.is_empty() {
        writeln!(writer, "0 0 0 0")?;
    } else {
        writeln!(
            writer,
            "{} {} {} {}",
            blocks.len(),
            elements.len(),
            min_element_tag,
            max_element_tag
        )?;
        for ((dimension, entity_tag, gmsh_type), block_elements) in blocks {
            writeln!(
                writer,
                "{} {} {} {}",
                dimension,
                entity_tag,
                gmsh_type,
                block_elements.len()
            )?;
            for element in block_elements {
                write!(writer, "{}", element.id)?;
                for node_id in &element.node_ids {
                    write!(writer, " {}", node_id)?;
                }
                writeln!(writer)?;
            }
        }
    }
    writeln!(writer, "$EndElements")?;
    if let Some(info) = &mesh.partitioned {
        crate::msh_partitioned::write_ghost_elements_v4(writer, &info.ghost_elements)?;
    }

    // Gmsh writes $Periodic after $Elements; empty links emit nothing.
    crate::msh_periodic::write_periodic_v4(writer, &mesh.periodic)?;

    Ok(())
}

/// Parse a Gmsh MSH file (v2.2 or v4.1 ASCII) from a reader.
pub fn parse_msh<R: BufRead>(reader: R) -> Result<Mesh, MshError> {
    let mut mesh = Mesh::new();
    let mut lines = reader.lines();
    let mut line_num: usize = 0;
    let mut version = MshVersion::V4;

    let next_line =
        |lines: &mut std::io::Lines<R>, line_num: &mut usize| -> Result<String, MshError> {
            *line_num += 1;
            lines
                .next()
                .ok_or_else(|| MshError::Parse {
                    line: *line_num,
                    message: "Unexpected end of file".into(),
                })?
                .map_err(MshError::Io)
        };

    while let Some(line_result) = lines.next() {
        line_num += 1;
        let line = line_result.map_err(MshError::Io)?;
        let trimmed = line.trim();

        match trimmed {
            "$MeshFormat" => {
                let format_line = next_line(&mut lines, &mut line_num)?;
                let parts: Vec<&str> = format_line.trim().split_whitespace().collect();
                if parts.is_empty() {
                    return Err(MshError::Parse {
                        line: line_num,
                        message: "Empty format line".into(),
                    });
                }
                let ver_str = parts[0];
                if ver_str.starts_with("2.") {
                    version = MshVersion::V2;
                } else if ver_str.starts_with("4.") {
                    version = MshVersion::V4;
                } else {
                    return Err(MshError::UnsupportedVersion(ver_str.into()));
                }

                let end = next_line(&mut lines, &mut line_num)?;
                if end.trim() != "$EndMeshFormat" {
                    return Err(MshError::Parse {
                        line: line_num,
                        message: "Expected $EndMeshFormat".into(),
                    });
                }
            }
            "$PhysicalNames" => {
                let count_line = next_line(&mut lines, &mut line_num)?;
                let count: usize = count_line.trim().parse().map_err(|_| MshError::Parse {
                    line: line_num,
                    message: "Invalid physical names count".into(),
                })?;
                for _ in 0..count {
                    let pn_line = next_line(&mut lines, &mut line_num)?;
                    let parts: Vec<&str> = pn_line.trim().splitn(3, ' ').collect();
                    if parts.len() >= 3 {
                        let dim: i32 = parts[0].parse().unwrap_or(0);
                        let tag: i32 = parts[1].parse().unwrap_or(0);
                        let name = parts[2].trim_matches('"').to_string();
                        mesh.physical_names.insert((dim, tag), name);
                    }
                }
                let end = next_line(&mut lines, &mut line_num)?;
                if end.trim() != "$EndPhysicalNames" {
                    return Err(MshError::Parse {
                        line: line_num,
                        message: "Expected $EndPhysicalNames".into(),
                    });
                }
            }
            "$Nodes" => match version {
                MshVersion::V2 => parse_nodes_v2(&mut lines, &mut line_num, &mut mesh)?,
                MshVersion::V4 => parse_nodes_v4(&mut lines, &mut line_num, &mut mesh)?,
            },
            "$Elements" => match version {
                MshVersion::V2 => parse_elements_v2(&mut lines, &mut line_num, &mut mesh)?,
                MshVersion::V4 => parse_elements_v4(&mut lines, &mut line_num, &mut mesh)?,
            },
            "$Periodic" => {
                let mut src = crate::msh_periodic::BufLineSource {
                    lines: &mut lines,
                    line_num: &mut line_num,
                };
                mesh.periodic = match version {
                    MshVersion::V2 => crate::msh_periodic::parse_periodic_v2(&mut src)?,
                    MshVersion::V4 => crate::msh_periodic::parse_periodic_v4(&mut src)?,
                };
            }
            "$PartitionedEntities" if version == MshVersion::V4 => {
                let mut section = Vec::new();
                loop {
                    let l = next_line(&mut lines, &mut line_num)?;
                    if l.trim() == "$EndPartitionedEntities" {
                        break;
                    }
                    section.push(l);
                }
                let refs: Vec<&str> = section.iter().map(|s| s.as_str()).collect();
                let mut info = crate::msh_partitioned::parse_partitioned_v4_ascii(
                    &refs,
                    line_num,
                )?;
                info.ghost_elements = mesh
                    .partitioned
                    .take()
                    .map(|p| p.ghost_elements)
                    .unwrap_or_default();
                mesh.partitioned = Some(info);
            }
            "$GhostElements" if version == MshVersion::V4 => {
                let mut section = Vec::new();
                loop {
                    let l = next_line(&mut lines, &mut line_num)?;
                    if l.trim() == "$EndGhostElements" {
                        break;
                    }
                    section.push(l);
                }
                let refs: Vec<&str> = section.iter().map(|s| s.as_str()).collect();
                let cells = crate::msh_partitioned::parse_ghost_elements_v4_ascii(
                    &refs,
                    line_num,
                )?;
                match mesh.partitioned.as_mut() {
                    Some(info) => info.ghost_elements = cells,
                    None => {
                        return Err(MshError::Parse {
                            line: line_num,
                            message: "$GhostElements without $PartitionedEntities".into(),
                        })
                    }
                }
            }
            _ => {
                if trimmed.starts_with('$') && !trimmed.starts_with("$End") {
                    let end_tag = format!("$End{}", &trimmed[1..]);
                    loop {
                        let skip_line = next_line(&mut lines, &mut line_num)?;
                        if skip_line.trim() == end_tag {
                            break;
                        }
                    }
                }
            }
        }
    }

    log::info!(
        "Parsed MSH: {} nodes, {} elements",
        mesh.node_count(),
        mesh.element_count()
    );

    Ok(mesh)
}

fn parse_nodes_v2<R: BufRead>(
    lines: &mut std::io::Lines<R>,
    line_num: &mut usize,
    mesh: &mut Mesh,
) -> Result<(), MshError> {
    let header = next_line_raw(lines, line_num)?;
    let num_nodes: usize = header.trim().parse().map_err(|_| MshError::Parse {
        line: *line_num,
        message: "Invalid node count".into(),
    })?;

    for _ in 0..num_nodes {
        let node_line = next_line_raw(lines, line_num)?;
        let parts: Vec<&str> = node_line.trim().split_whitespace().collect();
        if parts.len() < 4 {
            return Err(MshError::Parse {
                line: *line_num,
                message: "Invalid node line, expected: tag x y z".into(),
            });
        }
        let tag: u64 = parts[0].parse().map_err(|_| MshError::Parse {
            line: *line_num,
            message: "Invalid node tag".into(),
        })?;
        let x: f64 = parts[1].parse().map_err(|_| MshError::Parse {
            line: *line_num,
            message: "Invalid node x coordinate".into(),
        })?;
        let y: f64 = parts[2].parse().map_err(|_| MshError::Parse {
            line: *line_num,
            message: "Invalid node y coordinate".into(),
        })?;
        let z: f64 = parts[3].parse().map_err(|_| MshError::Parse {
            line: *line_num,
            message: "Invalid node z coordinate".into(),
        })?;
        mesh.add_node(Node::new(tag, x, y, z));
    }

    let end = next_line_raw(lines, line_num)?;
    if end.trim() != "$EndNodes" {
        return Err(MshError::Parse {
            line: *line_num,
            message: "Expected $EndNodes".into(),
        });
    }

    Ok(())
}

fn parse_elements_v2<R: BufRead>(
    lines: &mut std::io::Lines<R>,
    line_num: &mut usize,
    mesh: &mut Mesh,
) -> Result<(), MshError> {
    let header = next_line_raw(lines, line_num)?;
    let num_elements: usize = header.trim().parse().map_err(|_| MshError::Parse {
        line: *line_num,
        message: "Invalid element count".into(),
    })?;

    for _ in 0..num_elements {
        let elem_line = next_line_raw(lines, line_num)?;
        let parts: Vec<&str> = elem_line.trim().split_whitespace().collect();
        if parts.len() < 3 {
            return Err(MshError::Parse {
                line: *line_num,
                message: "Invalid element line".into(),
            });
        }
        let elem_tag: u64 = parts[0].parse().map_err(|_| MshError::Parse {
            line: *line_num,
            message: "Invalid element tag".into(),
        })?;
        let element_type_id: i32 = parts[1].parse().map_err(|_| MshError::Parse {
            line: *line_num,
            message: "Invalid element type".into(),
        })?;
        let num_tags: usize = parts[2].parse().map_err(|_| MshError::Parse {
            line: *line_num,
            message: "Invalid number of tags".into(),
        })?;

        let physical_tag = if num_tags > 0 {
            Some(parts[3].parse::<i32>().map_err(|_| MshError::Parse {
                line: *line_num,
                message: "Invalid physical tag".into(),
            })?)
        } else {
            None
        };

        let node_start = 3 + num_tags;
        if parts.len() < node_start {
            return Err(MshError::Parse {
                line: *line_num,
                message: "Element line too short for tags".into(),
            });
        }
        let node_ids: Vec<u64> = parts[node_start..]
            .iter()
            .map(|s| {
                s.parse::<u64>().map_err(|_| MshError::Parse {
                    line: *line_num,
                    message: "Invalid node id in element".into(),
                })
            })
            .collect::<Result<_, _>>()?;

        // The declared element type fixes the node count: a mismatch would
        // silently build a structurally corrupt element (e.g. a 2-node
        // Triangle3) that only detonates downstream in the meshing stack.
        let etype = ElementType::from_gmsh_type_id(element_type_id);
        let expected = etype.node_count();
        if expected != 0 && node_ids.len() != expected {
            return Err(MshError::Parse {
                line: *line_num,
                message: format!(
                    "Element {} declares type {} ({} node(s)) but lists {} node id(s)",
                    elem_tag,
                    element_type_id,
                    expected,
                    node_ids.len()
                ),
            });
        }

        let mut elem = Element::new(elem_tag, etype, node_ids);
        elem.physical_tag = physical_tag;
        mesh.add_element(elem);
    }

    let end = next_line_raw(lines, line_num)?;
    if end.trim() != "$EndElements" {
        return Err(MshError::Parse {
            line: *line_num,
            message: "Expected $EndElements".into(),
        });
    }

    Ok(())
}

fn parse_nodes_v4<R: BufRead>(
    lines: &mut std::io::Lines<R>,
    line_num: &mut usize,
    mesh: &mut Mesh,
) -> Result<(), MshError> {
    let header = next_line_raw(lines, line_num)?;
    let parts: Vec<&str> = header.trim().split_whitespace().collect();
    if parts.len() < 4 {
        return Err(MshError::Parse {
            line: *line_num,
            message: "Invalid nodes header".into(),
        });
    }
    let num_entity_blocks: usize = parts[0].parse().unwrap_or(0);

    for _ in 0..num_entity_blocks {
        let block_header = next_line_raw(lines, line_num)?;
        let bp: Vec<&str> = block_header.trim().split_whitespace().collect();
        if bp.len() < 4 {
            return Err(MshError::Parse {
                line: *line_num,
                message: "Invalid node block header".into(),
            });
        }
        let num_in_block: usize = bp[3].parse().unwrap_or(0);

        let mut tags = Vec::with_capacity(num_in_block);
        for _ in 0..num_in_block {
            let tag_line = next_line_raw(lines, line_num)?;
            let tag: u64 = tag_line.trim().parse().map_err(|_| MshError::Parse {
                line: *line_num,
                message: "Invalid node tag".into(),
            })?;
            tags.push(tag);
        }

        for tag in tags {
            let coord_line = next_line_raw(lines, line_num)?;
            let coords: Vec<f64> = coord_line
                .trim()
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if coords.len() >= 3 {
                mesh.add_node(Node::new(tag, coords[0], coords[1], coords[2]));
            }
        }
    }

    let end = next_line_raw(lines, line_num)?;
    if end.trim() != "$EndNodes" {
        return Err(MshError::Parse {
            line: *line_num,
            message: "Expected $EndNodes".into(),
        });
    }

    Ok(())
}

fn parse_elements_v4<R: BufRead>(
    lines: &mut std::io::Lines<R>,
    line_num: &mut usize,
    mesh: &mut Mesh,
) -> Result<(), MshError> {
    let header = next_line_raw(lines, line_num)?;
    let parts: Vec<&str> = header.trim().split_whitespace().collect();
    if parts.len() < 4 {
        return Err(MshError::Parse {
            line: *line_num,
            message: "Invalid elements header".into(),
        });
    }
    let num_entity_blocks: usize = parts[0].parse().unwrap_or(0);

    for _ in 0..num_entity_blocks {
        let block_header = next_line_raw(lines, line_num)?;
        let bp: Vec<&str> = block_header.trim().split_whitespace().collect();
        if bp.len() < 4 {
            return Err(MshError::Parse {
                line: *line_num,
                message: "Invalid element block header".into(),
            });
        }
        let entity_tag: i32 = bp[1].parse().map_err(|_| MshError::Parse {
            line: *line_num,
            message: "Invalid entity tag".into(),
        })?;
        let element_type_id: i32 = bp[2].parse().unwrap_or(0);
        let num_in_block: usize = bp[3].parse().unwrap_or(0);
        let etype = ElementType::from_gmsh_type_id(element_type_id);
        let expected_nodes = etype.node_count();

        for _ in 0..num_in_block {
            let elem_line = next_line_raw(lines, line_num)?;
            let values: Vec<u64> = elem_line
                .trim()
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if values.is_empty() {
                continue;
            }
            let elem_tag = values[0];
            let node_ids: Vec<u64> = values[1..].to_vec();

            if expected_nodes > 0 && node_ids.len() != expected_nodes {
                log::warn!(
                    "Element {} (type {:?}): expected {} nodes, got {}",
                    elem_tag,
                    etype,
                    expected_nodes,
                    node_ids.len()
                );
            }

            let mut elem = Element::new(elem_tag, etype, node_ids);
            // In MSH 4.x blocks are grouped by entity; rem uses this as physical tag.
            elem.physical_tag = Some(entity_tag);
            mesh.add_element(elem);
        }
    }

    let end = next_line_raw(lines, line_num)?;
    if end.trim() != "$EndElements" {
        return Err(MshError::Parse {
            line: *line_num,
            message: "Expected $EndElements".into(),
        });
    }

    Ok(())
}

fn next_line_raw<R: BufRead>(
    lines: &mut std::io::Lines<R>,
    line_num: &mut usize,
) -> Result<String, MshError> {
    *line_num += 1;
    lines
        .next()
        .ok_or_else(|| MshError::Parse {
            line: *line_num,
            message: "Unexpected end of file".into(),
        })?
        .map_err(MshError::Io)
}

fn validate_mesh(mesh: &Mesh) -> Result<(), MshError> {
    for element in &mesh.elements {
        gmsh_type_id(element.etype)?;
        // An element whose Gmsh type prescribes a node count must carry
        // exactly that many node ids — semantic high-order variants and
        // Unknown(id) forms alike: a mismatch would write a file no reader
        // can parse back.
        let expected = element.etype.node_count();
        if expected > 0 && element.node_ids.len() != expected {
            return Err(MshError::Parse {
                line: 0,
                message: format!(
                    "element {} declares type {:?} ({} node(s)) but carries {} node id(s)",
                    element.id,
                    element.etype,
                    expected,
                    element.node_ids.len()
                ),
            });
        }
        for node_id in &element.node_ids {
            if !mesh.nodes.contains_key(node_id) {
                return Err(MshError::MissingNode(*node_id));
            }
        }
    }
    Ok(())
}

fn write_physical_names<W: Write>(writer: &mut W, mesh: &Mesh) -> Result<(), MshError> {
    if mesh.physical_names.is_empty() {
        return Ok(());
    }

    let mut physical_names: Vec<_> = mesh.physical_names.iter().collect();
    physical_names.sort_by_key(|((dim, tag), _)| (*dim, *tag));

    writeln!(writer, "$PhysicalNames")?;
    writeln!(writer, "{}", physical_names.len())?;
    for ((dim, tag), name) in physical_names {
        writeln!(writer, "{} {} \"{}\"", dim, tag, name)?;
    }
    writeln!(writer, "$EndPhysicalNames")?;

    Ok(())
}

fn sorted_nodes(mesh: &Mesh) -> Vec<&rmsh_model::Node> {
    let mut nodes: Vec<_> = mesh.nodes.values().collect();
    nodes.sort_by_key(|node| node.id);
    nodes
}

/// The entities a v4.1 file must declare: every (dimension, entity tag)
/// pair the element blocks reference, with the bounding box of the tag's
/// nodes. Gmsh's element reader refuses blocks whose entity was never
/// declared in `$Entities` (its node reader auto-creates, its element
/// reader does not — pinned in GModelIO_MSH4.cpp), so the writer must
/// emit the section. Each entity carries itself as its single physical
/// tag, matching the writer's tag-as-physical convention.
fn declared_entities(
    mesh: &Mesh,
    elements: &[&rmsh_model::Element],
) -> BTreeMap<(u8, i32), ([f64; 3], [f64; 3])> {
    let mut out: BTreeMap<(u8, i32), ([f64; 3], [f64; 3])> = BTreeMap::new();
    for element in elements {
        let tag = element.physical_tag.unwrap_or(1);
        let entry = out
            .entry((element.dimension(), tag))
            .or_insert(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]));
        for node_id in &element.node_ids {
            if let Some(node) = mesh.nodes.get(node_id) {
                let p = [
                    node.position.x,
                    node.position.y,
                    node.position.z,
                ];
                for d in 0..3 {
                    entry.0[d] = entry.0[d].min(p[d]);
                    entry.1[d] = entry.1[d].max(p[d]);
                }
            }
        }
    }
    out
}

/// Write the v4.1 `$Entities` section (ASCII): the four class counts, then
/// per class one record per entity. Points carry three coordinates; every
/// other entity a six-value bounding box; each record closes with the
/// single self-referencing physical tag and (dim > 0) an empty bounding
/// list. Layout pinned against Gmsh 4.14's readMSH4EntityInfo.
fn write_entities_v4(
    writer: &mut impl Write,
    entities: &BTreeMap<(u8, i32), ([f64; 3], [f64; 3])>,
) -> Result<(), MshError> {
    writeln!(writer, "$Entities")?;
    let count = |d: u8| entities.keys().filter(|key| key.0 == d).count();
    writeln!(
        writer,
        "{} {} {} {}",
        count(0),
        count(1),
        count(2),
        count(3)
    )?;
    for ((dim, tag), (min, max)) in entities {
        if *dim == 0 {
            writeln!(writer, "{} {} {} {} 1 {}", tag, min[0], min[1], min[2], tag)?;
        } else {
            writeln!(
                writer,
                "{} {} {} {} {} {} {} 1 {} 0",
                tag, min[0], min[1], min[2], max[0], max[1], max[2], tag
            )?;
        }
    }
    writeln!(writer, "$EndEntities")?;
    Ok(())
}

/// Binary form of [`write_entities_v4`]: u64 class counts, then per entity
/// an i32 tag, 3 (points) or 6 (dim > 0) f64 box values, a u64 physical
/// count with the i32 self tag, and for dim > 0 a u64 zero bounding count.
fn write_entities_v4_binary(
    writer: &mut impl Write,
    entities: &BTreeMap<(u8, i32), ([f64; 3], [f64; 3])>,
) -> Result<(), MshError> {
    writeln!(writer, "$Entities")?;
    let count = |d: u8| entities.keys().filter(|key| key.0 == d).count();
    le_u64(writer, count(0) as u64)?;
    le_u64(writer, count(1) as u64)?;
    le_u64(writer, count(2) as u64)?;
    le_u64(writer, count(3) as u64)?;
    for ((dim, tag), (min, max)) in entities {
        le_i32(writer, i32::from(*tag))?;
        if *dim == 0 {
            le_f64(writer, min[0])?;
            le_f64(writer, min[1])?;
            le_f64(writer, min[2])?;
        } else {
            for value in min.iter().chain(max.iter()) {
                le_f64(writer, *value)?;
            }
        }
        le_u64(writer, 1)?;
        le_i32(writer, *tag)?;
        if *dim > 0 {
            le_u64(writer, 0)?;
        }
    }
    writeln!(writer, "$EndEntities")?;
    Ok(())
}

fn sorted_elements(mesh: &Mesh) -> Vec<&rmsh_model::Element> {
    let mut elements: Vec<_> = mesh.elements.iter().collect();
    elements.sort_by_key(|element| element.id);
    elements
}

fn gmsh_type_id(element_type: ElementType) -> Result<i32, MshError> {
    // The model layer owns the canonical type-id table (semantic variants
    // and the Unknown fidelity form share one mapping); an id with no
    // verified node-count entry has no valid Gmsh element mapping.
    element_type
        .to_gmsh_type_id()
        .ok_or(MshError::UnsupportedElementType(element_type))
}

fn skip_to_marker(pos: &mut usize, bytes: &[u8], marker: &[u8]) {
    while *pos + marker.len() <= bytes.len() {
        if bytes[*pos..].starts_with(marker) {
            *pos += marker.len();
            if *pos < bytes.len() && bytes[*pos] == b'\n' {
                *pos += 1;
            }
            return;
        }
        *pos += 1;
    }
}

/// Node count for a v2 chunked element type: the shared verified Gmsh
/// table first, then the legacy local fallback for the exotic ids the
/// verification sweep could not cover (94-117, 120-131 — absent from both
/// the official manual's abridged list and real Gmsh 4.14 output, kept so
/// such files still parse at their historical counts).
fn chunk_node_count(element_type_id: i32) -> Option<usize> {
    rmsh_model::gmsh_node_count(element_type_id)
        .or_else(|| v2_nodes_per_element(element_type_id))
}

/// Legacy fallback table, trimmed to the ids NOT covered by
/// [`rmsh_model::gmsh_node_count`] (everything verified there was deleted
/// from here so the two tables cannot diverge).
fn v2_nodes_per_element(element_type_id: i32) -> Option<usize> {
    match element_type_id {
        94 => Some(216),
        95 => Some(343),
        96 => Some(512),
        97 => Some(729),
        98 => Some(1000),
        99 => Some(32),
        100 => Some(44),
        101 => Some(56),
        102 => Some(68),
        103 => Some(80),
        104 => Some(92),
        105 => Some(104),
        106 => Some(126),
        107 => Some(196),
        108 => Some(288),
        109 => Some(405),
        110 => Some(550),
        111 => Some(24),
        112 => Some(33),
        113 => Some(42),
        114 => Some(51),
        115 => Some(60),
        116 => Some(69),
        117 => Some(78),
        120 => Some(91),
        121 => Some(140),
        122 => Some(204),
        123 => Some(285),
        124 => Some(385),
        125 => Some(21),
        126 => Some(29),
        127 => Some(37),
        128 => Some(45),
        129 => Some(53),
        130 => Some(61),
        131 => Some(69),
        _ => None,
    }
}
fn parse_msh_v2_binary(bytes: &[u8]) -> Result<Mesh, MshError> {
    let mut pos = 0usize;
    let mut mesh = Mesh::new();

    let read_line = |pos: &mut usize| -> Option<&str> {
        let start = *pos;
        let slice = &bytes[start..];
        if let Some(nl) = slice.iter().position(|&b| b == b'\n') {
            let line = std::str::from_utf8(&slice[..nl]).ok()?.trim();
            *pos = start + nl + 1;
            Some(line)
        } else if !slice.is_empty() {
            let line = std::str::from_utf8(slice).ok()?.trim();
            *pos = bytes.len();
            Some(line)
        } else {
            None
        }
    };

    let read_i32_le = |pos: &mut usize| -> Result<i32, MshError> {
        if *pos + 4 > bytes.len() {
            return Err(MshError::Parse {
                line: 0,
                message: "binary v2: unexpected end reading int32".into(),
            });
        }
        let v = i32::from_le_bytes([bytes[*pos], bytes[*pos + 1], bytes[*pos + 2], bytes[*pos + 3]]);
        *pos += 4;
        Ok(v)
    };

    let read_f64_le = |pos: &mut usize| -> Result<f64, MshError> {
        if *pos + 8 > bytes.len() {
            return Err(MshError::Parse {
                line: 0,
                message: "binary v2: unexpected end reading float64".into(),
            });
        }
        let v = f64::from_le_bytes([
            bytes[*pos],
            bytes[*pos + 1],
            bytes[*pos + 2],
            bytes[*pos + 3],
            bytes[*pos + 4],
            bytes[*pos + 5],
            bytes[*pos + 6],
            bytes[*pos + 7],
        ]);
        *pos += 8;
        Ok(v)
    };

    loop {
        let Some(line) = read_line(&mut pos) else {
            break;
        };
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        match line {
            "$MeshFormat" => {
                let format_line = read_line(&mut pos).unwrap_or_default();
                // Binary files carry the endianness magic right after the
                // format line: validate it (little endian only) instead of
                // letting the line scan swallow the bytes silently.
                let parts: Vec<&str> = format_line.split_whitespace().collect();
                let is_binary = parts.len() >= 2 && parts[1] == "1";
                if is_binary {
                    read_binary_magic(bytes, &mut pos, "v2")?;
                }
                skip_to_marker(&mut pos, bytes, b"$EndMeshFormat");
            }
            "$PhysicalNames" => {
                if let Some(count_line) = read_line(&mut pos) {
                    let n_names: usize = count_line.trim().parse().unwrap_or(0);
                    for _ in 0..n_names {
                        if let Some(name_line) = read_line(&mut pos) {
                            let parts: Vec<&str> = name_line.splitn(3, ' ').collect();
                            if parts.len() >= 3 {
                                let dim = parts[0].parse::<i32>().unwrap_or(0);
                                let tag = parts[1].parse::<i32>().unwrap_or(0);
                                let name = parts[2].trim().trim_matches('"').to_string();
                                mesh.physical_names.insert((dim, tag), name);
                            }
                        }
                    }
                }
                skip_to_marker(&mut pos, bytes, b"$EndPhysicalNames");
            }
            "$Periodic" => {
                // Gmsh keeps the MSH2 $Periodic section in ASCII even inside
                // binary files (only $Nodes/$Elements payloads go binary),
                // so the ASCII-form parser runs against the byte cursor.
                let mut src = crate::msh_periodic::ByteLineSource { bytes, pos: &mut pos };
                mesh.periodic = crate::msh_periodic::parse_periodic_v2(&mut src)?;
            }
            "$Nodes" => {
                let Some(count_line) = read_line(&mut pos) else {
                    break;
                };
                let n_nodes: usize = count_line.trim().parse().map_err(|_| MshError::Parse {
                    line: 0,
                    message: "binary v2: invalid node count".into(),
                })?;
                for _ in 0..n_nodes {
                    let id = read_i32_le(&mut pos)? as u64;
                    let x = read_f64_le(&mut pos)?;
                    let y = read_f64_le(&mut pos)?;
                    let z = read_f64_le(&mut pos)?;
                    mesh.add_node(Node::new(id, x, y, z));
                }

                if pos < bytes.len() && bytes[pos] == b'\n' {
                    pos += 1;
                }
                if pos + 9 <= bytes.len() && &bytes[pos..pos + 9] == b"$EndNodes" {
                    pos += 9;
                    if pos < bytes.len() && bytes[pos] == b'\n' {
                        pos += 1;
                    }
                }
            }
            "$Elements" => {
                let Some(total_line) = read_line(&mut pos) else {
                    break;
                };
                let total_elements: usize = total_line.trim().parse().map_err(|_| MshError::Parse {
                    line: 0,
                    message: "binary v2: invalid total element count".into(),
                })?;

                let mut parsed_elements = 0usize;
                while parsed_elements < total_elements {
                    if pos + b"$EndElements".len() <= bytes.len()
                        && bytes[pos..].starts_with(b"$EndElements")
                    {
                        break;
                    }

                    let element_type_id = read_i32_le(&mut pos)?;
                    let n_elems = read_i32_le(&mut pos)? as usize;
                    let n_tags = read_i32_le(&mut pos)? as usize;
                    let Some(n_nodes_per) = chunk_node_count(element_type_id) else {
                        return Err(MshError::Parse {
                            line: 0,
                            message: format!("binary v2: unsupported element type {}", element_type_id),
                        });
                    };

                    let etype = ElementType::from_gmsh_type_id(element_type_id);
                    for _ in 0..n_elems {
                        let elem_id = read_i32_le(&mut pos)? as u64;
                        let mut physical_tag: Option<i32> = None;
                        for ti in 0..n_tags {
                            let tag = read_i32_le(&mut pos)?;
                            if ti == 0 {
                                physical_tag = Some(tag);
                            }
                        }

                        let mut node_ids = Vec::with_capacity(n_nodes_per);
                        for _ in 0..n_nodes_per {
                            node_ids.push(read_i32_le(&mut pos)? as u64);
                        }

                        let mut elem = Element::new(elem_id, etype, node_ids);
                        elem.physical_tag = physical_tag;
                        mesh.add_element(elem);
                    }

                    parsed_elements += n_elems;
                }

                if parsed_elements != total_elements {
                    return Err(MshError::Parse {
                        line: 0,
                        message: format!(
                            "binary v2: expected {} elements, parsed {}",
                            total_elements, parsed_elements
                        ),
                    });
                }

                // Consume trailing newline and $EndElements marker.
                if pos < bytes.len() && bytes[pos] == b'\n' {
                    pos += 1;
                }
                if pos + b"$EndElements".len() <= bytes.len()
                    && bytes[pos..].starts_with(b"$EndElements")
                {
                    pos += b"$EndElements".len();
                    if pos < bytes.len() && bytes[pos] == b'\n' {
                        pos += 1;
                    }
                } else {
                    skip_to_marker(&mut pos, bytes, b"$EndElements");
                }
            }
            _ if line.starts_with('$') => {
                let end_marker = format!("$End{}", &line[1..]);
                skip_to_marker(&mut pos, bytes, end_marker.as_bytes());
            }
            _ => {}
        }
    }

    if mesh.nodes.is_empty() {
        // An empty mesh is a legitimate (if useless) file: the format line
        // said binary and the sections parsed — return the empty mesh.
        return Ok(mesh);
    }
    if mesh.elements.is_empty() {
        return Err(MshError::Parse {
            line: 0,
            message: "binary v2: no elements found".into(),
        });
    }

    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn sample_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 0.0, 1.0));

        let mut tri = Element::new(1, ElementType::Triangle3, vec![1, 2, 3]);
        tri.physical_tag = Some(11);
        mesh.add_element(tri);
        mesh.add_element(Element::new(2, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));
        mesh.physical_names.insert((2, 11), "surface".to_string());
        mesh
    }

    fn assert_mesh_core_eq(actual: &Mesh, expected: &Mesh) {
        assert_eq!(actual.node_count(), expected.node_count());
        assert_eq!(actual.element_count(), expected.element_count());
        assert_eq!(actual.physical_names, expected.physical_names);

        let mut actual_nodes: Vec<_> = actual.nodes.iter().collect();
        actual_nodes.sort_by_key(|(id, _)| **id);
        let mut expected_nodes: Vec<_> = expected.nodes.iter().collect();
        expected_nodes.sort_by_key(|(id, _)| **id);
        for ((actual_id, actual_node), (expected_id, expected_node)) in
            actual_nodes.into_iter().zip(expected_nodes)
        {
            assert_eq!(actual_id, expected_id);
            assert_eq!(actual_node.position, expected_node.position);
        }

        let mut actual_elements: Vec<_> = actual.elements.iter().collect();
        actual_elements.sort_by_key(|element| element.id);
        let mut expected_elements: Vec<_> = expected.elements.iter().collect();
        expected_elements.sort_by_key(|element| element.id);
        for (actual_element, expected_element) in actual_elements.into_iter().zip(expected_elements)
        {
            assert_eq!(actual_element.id, expected_element.id);
            assert_eq!(actual_element.etype, expected_element.etype);
            assert_eq!(actual_element.node_ids, expected_element.node_ids);
        }
    }

    fn temp_msh_path(stem: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time should be after unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!("{}_{}_{}.msh", env!("CARGO_PKG_NAME"), stem, unique))
    }

    #[test]
    fn test_parse_simple_msh_v4() {
        let msh_data = r#"$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
$EndNodes
$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
$EndElements
"#;
        let mesh = parse_msh(Cursor::new(msh_data.as_bytes())).unwrap();
        assert_eq!(mesh.node_count(), 4);
        assert_eq!(mesh.element_count(), 1);
        assert_eq!(
            mesh.elements[0].etype,
            rmsh_model::ElementType::Tetrahedron4
        );
    }

    #[test]
    fn test_parse_simple_msh_v2() {
        let msh_data = r#"$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 0.0 1.0 0.0
4 0.0 0.0 1.0
$EndNodes
$Elements
2
1 2 2 0 1 1 2 3
2 4 2 0 1 1 2 3 4
$EndElements
"#;
        let mesh = parse_msh(Cursor::new(msh_data.as_bytes())).unwrap();
        assert_eq!(mesh.node_count(), 4);
        assert_eq!(mesh.element_count(), 2);
        assert_eq!(mesh.elements[0].etype, rmsh_model::ElementType::Triangle3);
        assert_eq!(
            mesh.elements[1].etype,
            rmsh_model::ElementType::Tetrahedron4
        );
        assert_eq!(mesh.elements[0].node_ids, vec![1, 2, 3]);
        assert_eq!(mesh.elements[1].node_ids, vec![1, 2, 3, 4]);
        assert_eq!(mesh.elements[0].physical_tag, Some(0));
        assert_eq!(mesh.elements[1].physical_tag, Some(0));
    }

    #[test]
    fn test_parse_msh_v4_sets_entity_as_physical_tag() {
        let msh_data = r#"$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 1 3
2 1 0 3
1
2
3
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
$EndNodes
$Elements
1 1 1 1
2 42 2 1
1 1 2 3
$EndElements
"#;
        let mesh = parse_msh(Cursor::new(msh_data.as_bytes())).unwrap();
        assert_eq!(mesh.elements.len(), 1);
        assert_eq!(mesh.elements[0].physical_tag, Some(42));
    }

    #[test]
    fn test_load_msh_from_bytes_binary_v2() {
        let mut data = Vec::<u8>::new();
        data.extend_from_slice(b"$MeshFormat\n2.2 1 8\n");
        data.extend_from_slice(&1_i32.to_le_bytes());
        data.extend_from_slice(&1.0f64.to_le_bytes());
        data.extend_from_slice(b"\n$EndMeshFormat\n");

        data.extend_from_slice(b"$Nodes\n2\n");
        data.extend_from_slice(&1_i32.to_le_bytes());
        data.extend_from_slice(&0.0f64.to_le_bytes());
        data.extend_from_slice(&0.0f64.to_le_bytes());
        data.extend_from_slice(&0.0f64.to_le_bytes());
        data.extend_from_slice(&2_i32.to_le_bytes());
        data.extend_from_slice(&1.0f64.to_le_bytes());
        data.extend_from_slice(&0.0f64.to_le_bytes());
        data.extend_from_slice(&0.0f64.to_le_bytes());
        data.extend_from_slice(b"\n$EndNodes\n");

        data.extend_from_slice(b"$Elements\n1\n");
        data.extend_from_slice(&1_i32.to_le_bytes()); // type line2
        data.extend_from_slice(&1_i32.to_le_bytes()); // n elems
        data.extend_from_slice(&2_i32.to_le_bytes()); // n tags
        data.extend_from_slice(&1_i32.to_le_bytes()); // elem id
        data.extend_from_slice(&77_i32.to_le_bytes()); // physical tag
        data.extend_from_slice(&0_i32.to_le_bytes()); // geometrical tag
        data.extend_from_slice(&1_i32.to_le_bytes()); // n1
        data.extend_from_slice(&2_i32.to_le_bytes()); // n2
        data.extend_from_slice(b"$EndElements\n");

        let mesh = load_msh_from_bytes(&data).expect("binary v2 should parse");
        assert_eq!(mesh.node_count(), 2);
        assert_eq!(mesh.element_count(), 1);
        assert_eq!(mesh.elements[0].physical_tag, Some(77));
        assert_eq!(mesh.elements[0].node_ids, vec![1, 2]);
    }

    #[test]
    fn test_load_msh_from_bytes_binary_v2_first_block_type_10() {
        let mut data = Vec::<u8>::new();
        data.extend_from_slice(b"$MeshFormat\n2.2 1 8\n");
        data.extend_from_slice(&1_i32.to_le_bytes());
        data.extend_from_slice(&1.0f64.to_le_bytes());
        data.extend_from_slice(b"\n$EndMeshFormat\n");

        data.extend_from_slice(b"$Nodes\n9\n");
        for (id, x, y, z) in [
            (1_i32, 0.0_f64, 0.0_f64, 0.0_f64),
            (2, 1.0_f64, 0.0_f64, 0.0_f64),
            (3, 0.0_f64, 1.0_f64, 0.0_f64),
            (4, 0.5_f64, 0.0_f64, 0.0_f64),
            (5, 0.5_f64, 0.5_f64, 0.0_f64),
            (6, 0.0_f64, 0.5_f64, 0.0_f64),
            (7, 0.25_f64, 0.0_f64, 0.0_f64),
            (8, 0.75_f64, 0.0_f64, 0.0_f64),
            (9, 0.5_f64, 0.25_f64, 0.0_f64),
        ] {
            data.extend_from_slice(&id.to_le_bytes());
            data.extend_from_slice(&x.to_le_bytes());
            data.extend_from_slice(&y.to_le_bytes());
            data.extend_from_slice(&z.to_le_bytes());
        }
        data.extend_from_slice(b"\n$EndNodes\n");

        data.extend_from_slice(b"$Elements\n1\n");
        data.extend_from_slice(&10_i32.to_le_bytes()); // type 10, low byte is 0x0A
        data.extend_from_slice(&1_i32.to_le_bytes()); // n elems
        data.extend_from_slice(&2_i32.to_le_bytes()); // n tags
        data.extend_from_slice(&1_i32.to_le_bytes()); // elem id
        data.extend_from_slice(&88_i32.to_le_bytes()); // physical tag
        data.extend_from_slice(&0_i32.to_le_bytes()); // geometrical tag
        for node_id in [1_i32, 2, 3, 4, 5, 6, 7, 8, 9] {
            data.extend_from_slice(&node_id.to_le_bytes());
        }
        data.extend_from_slice(b"$EndElements\n");

        let mesh = load_msh_from_bytes(&data).expect("binary v2 type-10 first block should parse");
        assert_eq!(mesh.node_count(), 9);
        assert_eq!(mesh.element_count(), 1);
        assert_eq!(mesh.elements[0].physical_tag, Some(88));
        assert_eq!(mesh.elements[0].node_ids, vec![1, 2, 3, 4, 5, 6, 7, 8, 9]);
    }

    #[test]
    fn test_load_msh_from_bytes_binary_v2_type_10_then_type_12() {
        let mut data = Vec::<u8>::new();
        data.extend_from_slice(b"$MeshFormat\n2.2 1 8\n");
        data.extend_from_slice(&1_i32.to_le_bytes());
        data.extend_from_slice(&1.0f64.to_le_bytes());
        data.extend_from_slice(b"\n$EndMeshFormat\n");

        data.extend_from_slice(b"$Nodes\n27\n");
        for id in 1_i32..=27_i32 {
            let coord = id as f64;
            data.extend_from_slice(&id.to_le_bytes());
            data.extend_from_slice(&coord.to_le_bytes());
            data.extend_from_slice(&0.0_f64.to_le_bytes());
            data.extend_from_slice(&0.0_f64.to_le_bytes());
        }
        data.extend_from_slice(b"\n$EndNodes\n");

        data.extend_from_slice(b"$Elements\n2\n");

        data.extend_from_slice(&10_i32.to_le_bytes());
        data.extend_from_slice(&1_i32.to_le_bytes());
        data.extend_from_slice(&2_i32.to_le_bytes());
        data.extend_from_slice(&1_i32.to_le_bytes());
        data.extend_from_slice(&88_i32.to_le_bytes());
        data.extend_from_slice(&0_i32.to_le_bytes());
        for node_id in [1_i32, 2, 3, 4, 5, 6, 7, 8, 9] {
            data.extend_from_slice(&node_id.to_le_bytes());
        }

        data.extend_from_slice(&12_i32.to_le_bytes());
        data.extend_from_slice(&1_i32.to_le_bytes());
        data.extend_from_slice(&2_i32.to_le_bytes());
        data.extend_from_slice(&2_i32.to_le_bytes());
        data.extend_from_slice(&99_i32.to_le_bytes());
        data.extend_from_slice(&0_i32.to_le_bytes());
        for node_id in 1_i32..=27_i32 {
            data.extend_from_slice(&node_id.to_le_bytes());
        }
        data.extend_from_slice(b"$EndElements\n");

        let mesh = load_msh_from_bytes(&data).expect("binary v2 type-10 and type-12 blocks should parse");
        assert_eq!(mesh.node_count(), 27);
        assert_eq!(mesh.element_count(), 2);
        assert_eq!(mesh.elements[0].physical_tag, Some(88));
        assert_eq!(mesh.elements[1].physical_tag, Some(99));
        assert_eq!(mesh.elements[1].node_ids.len(), 27);
    }

    #[test]
    fn test_load_msh_from_bytes_binary_v2_preserves_unknown_high_order_elements() {
        // Type 119 (Pyramid55, order 4) is outside the P5 semantic variant
        // set: it must keep the Unknown(id) fidelity form.
        let mut data = Vec::<u8>::new();
        data.extend_from_slice(b"$MeshFormat\n2.2 1 8\n");
        data.extend_from_slice(&1_i32.to_le_bytes());
        data.extend_from_slice(&1.0f64.to_le_bytes());
        data.extend_from_slice(b"\n$EndMeshFormat\n");

        data.extend_from_slice(b"$Nodes\n55\n");
        for id in 1_i32..=55_i32 {
            let coord = id as f64;
            data.extend_from_slice(&id.to_le_bytes());
            data.extend_from_slice(&coord.to_le_bytes());
            data.extend_from_slice(&0.0_f64.to_le_bytes());
            data.extend_from_slice(&0.0_f64.to_le_bytes());
        }
        data.extend_from_slice(b"\n$EndNodes\n");

        data.extend_from_slice(b"$Elements\n1\n");
        data.extend_from_slice(&119_i32.to_le_bytes());
        data.extend_from_slice(&1_i32.to_le_bytes());
        data.extend_from_slice(&2_i32.to_le_bytes());
        data.extend_from_slice(&7_i32.to_le_bytes());
        data.extend_from_slice(&123_i32.to_le_bytes());
        data.extend_from_slice(&0_i32.to_le_bytes());
        for node_id in 1_i32..=55_i32 {
            data.extend_from_slice(&node_id.to_le_bytes());
        }
        data.extend_from_slice(b"$EndElements\n");

        let mesh = load_msh_from_bytes(&data)
            .expect("binary v2 unknown high-order element should be preserved");
        assert_eq!(mesh.node_count(), 55);
        assert_eq!(mesh.element_count(), 1);
        assert_eq!(mesh.elements[0].etype, rmsh_model::ElementType::Unknown(119));
        assert_eq!(mesh.elements[0].physical_tag, Some(123));
        assert_eq!(mesh.elements[0].node_ids.len(), 55);
    }

    #[test]
    fn test_write_roundtrip_msh_v2() {
        let mesh = sample_mesh();
        let mut output = Vec::new();
        write_msh_v2(&mut output, &mesh).unwrap();

        let parsed = parse_msh(Cursor::new(output)).unwrap();
        assert_mesh_core_eq(&parsed, &mesh);
    }

    #[test]
    fn test_write_roundtrip_msh_v4() {
        let mesh = sample_mesh();
        let mut output = Vec::new();
        write_msh_v4(&mut output, &mesh).unwrap();

        let parsed = parse_msh(Cursor::new(output)).unwrap();
        assert_mesh_core_eq(&parsed, &mesh);
    }

    #[test]
    fn test_save_msh_v4_to_path_and_load_msh_from_path_roundtrip() {
        let mesh = sample_mesh();
        let path = temp_msh_path("roundtrip_v4");

        save_msh_v4_to_path(&path, &mesh).expect("save_msh_v4_to_path should succeed");
        let parsed = load_msh_from_path(&path).expect("load_msh_from_path should succeed");

        assert_mesh_core_eq(&parsed, &mesh);

        std::fs::remove_file(&path).expect("temporary msh file should be removable");
    }

    #[test]
    fn test_load_msh_from_path_reports_io_error_for_missing_file() {
        let path = temp_msh_path("missing");
        let err = load_msh_from_path(&path).expect_err("missing file should return io error");

        match err {
            MshError::Io(io_err) => {
                assert_eq!(io_err.kind(), std::io::ErrorKind::NotFound);
            }
            other => panic!("expected MshError::Io, got {other:?}"),
        }
    }
}
// ── Binary MSH (v2.2 chunked + v4.1) ────────────────────────────────────────
//
// The v2 binary layout is the CHUNKED convention the reader above (and its
// four tests) already pin: the `$Elements` payload groups per element TYPE
// as `i32 type; i32 num_elements; i32 num_tags;` followed by
// `num_elements x (i32 id; num_tags x i32 tags; nodes x i32)`. The v4.1
// binary layout (verified against real Gmsh 4.14 output, which silently
// ingests nothing from a deviant file) keeps ONLY the section marker lines
// ASCII: section headers, entity block headers and payloads are all
// little-endian binary (`size_t` = u64, int = i32, coordinates = f64).
// Both flavors open with the endianness magic right after the format line:
// readers accept the documented 12-byte form (`i32 1; f64 1.0`) and the
// compact 4-byte form (`i32 1`) that real Gmsh 4.14 writes; writers emit
// the documented 12-byte form. Big-endian files are refused with a
// located, naming error (little endian only).

pub(crate) fn le_i32<W: Write>(writer: &mut W, v: i32) -> Result<(), MshError> {
    writer.write_all(&v.to_le_bytes()).map_err(MshError::Io)
}

pub(crate) fn le_u64<W: Write>(writer: &mut W, v: u64) -> Result<(), MshError> {
    writer.write_all(&v.to_le_bytes()).map_err(MshError::Io)
}

pub(crate) fn le_f64<W: Write>(writer: &mut W, v: f64) -> Result<(), MshError> {
    writer.write_all(&v.to_le_bytes()).map_err(MshError::Io)
}

pub(crate) fn bin_read_i32(bytes: &[u8], pos: &mut usize, what: &str) -> Result<i32, MshError> {
    if *pos + 4 > bytes.len() {
        return Err(MshError::Parse {
            line: 0,
            message: format!("binary: unexpected end reading int32 ({what})"),
        });
    }
    let v = i32::from_le_bytes([bytes[*pos], bytes[*pos + 1], bytes[*pos + 2], bytes[*pos + 3]]);
    *pos += 4;
    Ok(v)
}

pub(crate) fn bin_read_u64(bytes: &[u8], pos: &mut usize, what: &str) -> Result<u64, MshError> {
    if *pos + 8 > bytes.len() {
        return Err(MshError::Parse {
            line: 0,
            message: format!("binary: unexpected end reading size_t ({what})"),
        });
    }
    let mut buf = [0u8; 8];
    buf.copy_from_slice(&bytes[*pos..*pos + 8]);
    *pos += 8;
    Ok(u64::from_le_bytes(buf))
}

pub(crate) fn bin_read_f64(bytes: &[u8], pos: &mut usize, what: &str) -> Result<f64, MshError> {
    if *pos + 8 > bytes.len() {
        return Err(MshError::Parse {
            line: 0,
            message: format!("binary: unexpected end reading float64 ({what})"),
        });
    }
    let mut buf = [0u8; 8];
    buf.copy_from_slice(&bytes[*pos..*pos + 8]);
    *pos += 8;
    Ok(f64::from_le_bytes(buf))
}

fn bin_read_line<'a>(bytes: &'a [u8], pos: &mut usize) -> Result<&'a str, MshError> {
    let slice = &bytes[*pos..];
    let Some(nl) = slice.iter().position(|&b| b == b'\n') else {
        return Err(MshError::Parse {
            line: 0,
            message: "binary: unexpected end of file expecting a line".into(),
        });
    };
    let line = std::str::from_utf8(&slice[..nl]).map_err(|_| MshError::Parse {
        line: 0,
        message: "binary: section line is not UTF-8".into(),
    })?;
    *pos += nl + 1;
    Ok(line.trim())
}

fn bin_expect_marker(bytes: &[u8], pos: &mut usize, marker: &str) -> Result<(), MshError> {
    if *pos < bytes.len() && bytes[*pos] == b'\n' {
        *pos += 1;
    }
    let line = bin_read_line(bytes, pos)?;
    if line != marker {
        return Err(MshError::Parse {
            line: 0,
            message: format!("binary: expected {marker}, found {line}"),
        });
    }
    Ok(())
}

/// The binary endianness magic: `i32 1` followed by `f64 1.0` (little
/// endian), written right after every binary `format_line`. (Readers also
/// accept the compact `i32 1`-only form real Gmsh 4.14 emits; see
/// [`read_binary_magic`].)
fn write_binary_magic<W: Write>(writer: &mut W) -> Result<(), MshError> {
    le_i32(writer, 1)?;
    le_f64(writer, 1.0)
}

/// Validates the endianness magic at the cursor and refuses big-endian /
/// corrupt files with a located, naming error. Two little-endian forms are
/// accepted:
/// - the 12-byte form `i32 1; f64 1.0` (the documented magic; what this
///   crate writes),
/// - the compact 4-byte form `i32 1` that real Gmsh 4.14 writes (verified
///   against its v2.2 and v4.1 binary output; that form is followed by the
///   `$EndMeshFormat` newline, whose first byte never reads as f64 1.0, so
///   the two forms cannot be confused).
fn read_binary_magic(bytes: &[u8], pos: &mut usize, flavor: &str) -> Result<(), MshError> {
    if *pos + 4 > bytes.len() {
        return Err(MshError::Parse {
            line: 0,
            message: format!("binary {flavor}: truncated endianness magic"),
        });
    }
    let flag = i32::from_le_bytes([bytes[*pos], bytes[*pos + 1], bytes[*pos + 2], bytes[*pos + 3]]);
    if flag != 1 {
        let be_flag =
            i32::from_be_bytes([bytes[*pos], bytes[*pos + 1], bytes[*pos + 2], bytes[*pos + 3]]);
        return Err(MshError::Parse {
            line: 0,
            message: format!(
                "binary {flavor}: bad endianness magic (i32 {flag}, f64 flag reads as {be_flag}) - big-endian files are not supported; expected little-endian i32 1, f64 1.0"
            ),
        });
    }
    if *pos + 12 <= bytes.len() {
        let one = f64::from_le_bytes(bytes[*pos + 4..*pos + 12].try_into().unwrap());
        if (one - 1.0).abs() < f64::EPSILON {
            *pos += 12;
            return Ok(());
        }
    }
    *pos += 4;
    Ok(())
}

fn v2_int32_id(id: u64) -> Result<i32, MshError> {
    i32::try_from(id).map_err(|_| MshError::Parse {
        line: 0,
        message: format!("binary v2: node/element id {id} exceeds the 32-bit v2 range"),
    })
}

/// Write a Gmsh MSH v2.2 BINARY file (chunked element blocks, little
/// endian). Tag semantics mirror [`write_msh_v2`]: elements with a
/// physical tag carry 2 tags `[physical, 0]`, elements without carry 0
/// tags; a chunk groups one (type, has-physical) class so the per-chunk
/// `num_tags` is uniform, matching what [`parse_msh_v2_binary`] reads.
pub fn write_msh_v2_binary<W: Write>(writer: &mut W, mesh: &Mesh) -> Result<(), MshError> {
    validate_mesh(mesh)?;

    let nodes = sorted_nodes(mesh);
    let elements = sorted_elements(mesh);

    writeln!(writer, "$MeshFormat")?;
    writeln!(writer, "2.2 1 8")?;
    write_binary_magic(writer)?;
    writeln!(writer)?;
    writeln!(writer, "$EndMeshFormat")?;
    write_physical_names(writer, mesh)?;

    writeln!(writer, "$Nodes")?;
    // The v2 binary layout keeps the section count lines in ASCII (the
    // payloads are binary); a binary count would make the reader's line
    // scan swallow the whole payload as one non-UTF-8 line.
    writeln!(writer, "{}", nodes.len())?;
    for node in &nodes {
        le_i32(writer, v2_int32_id(node.id)?)?;
        le_f64(writer, node.position.x)?;
        le_f64(writer, node.position.y)?;
        le_f64(writer, node.position.z)?;
    }
    writeln!(writer)?;
    writeln!(writer, "$EndNodes")?;

    let mut chunks: BTreeMap<(i32, bool), Vec<&rmsh_model::Element>> = BTreeMap::new();
    for element in &elements {
        let gmsh_type = gmsh_type_id(element.etype)?;
        chunks
            .entry((gmsh_type, element.physical_tag.is_some()))
            .or_default()
            .push(*element);
    }

    writeln!(writer, "$Elements")?;
    writeln!(writer, "{}", elements.len())?;
    for ((gmsh_type, has_physical), chunk) in chunks {
        le_i32(writer, gmsh_type)?;
        le_i32(writer, chunk.len() as i32)?;
        le_i32(writer, if has_physical { 2 } else { 0 })?;
        for element in &chunk {
            le_i32(writer, v2_int32_id(element.id)?)?;
            if has_physical {
                le_i32(writer, element.physical_tag.unwrap_or(0))?;
                le_i32(writer, 0)?;
            }
            for node_id in &element.node_ids {
                le_i32(writer, v2_int32_id(*node_id)?)?;
            }
        }
    }
    writeln!(writer)?;
    writeln!(writer, "$EndElements")?;

    // The $Periodic section stays ASCII in v2 binary files (the Gmsh
    // convention); empty links emit nothing.
    crate::msh_periodic::write_periodic_v2(writer, &mesh.periodic)?;

    Ok(())
}

/// Write a Gmsh MSH v4.1 BINARY file. Pinned against real Gmsh 4.14 (which
/// silently ingests nothing from a file whose section headers are ASCII):
/// only the section MARKER lines are ASCII; the section headers, entity
/// block headers and payloads are little-endian binary (`size_t` = u64,
/// int = i32, coordinates = f64). The block structure mirrors
/// [`write_msh_v4`] exactly.
pub fn write_msh_v4_binary<W: Write>(writer: &mut W, mesh: &Mesh) -> Result<(), MshError> {
    validate_mesh(mesh)?;

    let nodes = sorted_nodes(mesh);
    let elements = sorted_elements(mesh);
    let min_node_tag = nodes.first().map(|node| node.id).unwrap_or(0);
    let max_node_tag = nodes.last().map(|node| node.id).unwrap_or(0);
    let min_element_tag = elements.first().map(|element| element.id).unwrap_or(0);
    let max_element_tag = elements.last().map(|element| element.id).unwrap_or(0);
    let entity_dim = elements
        .iter()
        .map(|element| element.dimension() as i32)
        .max()
        .unwrap_or(0);

    writeln!(writer, "$MeshFormat")?;
    writeln!(writer, "4.1 1 8")?;
    write_binary_magic(writer)?;
    writeln!(writer)?;
    writeln!(writer, "$EndMeshFormat")?;
    write_physical_names(writer, mesh)?;
    write_entities_v4_binary(writer, &declared_entities(mesh, &elements))?;
    if let Some(info) = &mesh.partitioned {
        crate::msh_partitioned::write_partitioned_v4_binary(writer, info)?;
    }

    writeln!(writer, "$Nodes")?;
    if nodes.is_empty() {
        le_u64(writer, 0)?;
        le_u64(writer, 0)?;
        le_u64(writer, 0)?;
        le_u64(writer, 0)?;
    } else {
        le_u64(writer, 1)?;
        le_u64(writer, nodes.len() as u64)?;
        le_u64(writer, min_node_tag)?;
        le_u64(writer, max_node_tag)?;
        le_i32(writer, entity_dim)?;
        le_i32(writer, 1)?;
        le_i32(writer, 0)?;
        le_u64(writer, nodes.len() as u64)?;
        for node in &nodes {
            le_u64(writer, node.id)?;
        }
        for node in &nodes {
            le_f64(writer, node.position.x)?;
            le_f64(writer, node.position.y)?;
            le_f64(writer, node.position.z)?;
        }
    }
    writeln!(writer, "$EndNodes")?;

    let mut blocks: BTreeMap<(u8, i32, i32), Vec<&rmsh_model::Element>> = BTreeMap::new();
    for element in &elements {
        let gmsh_type = gmsh_type_id(element.etype)?;
        let entity_tag = element.physical_tag.unwrap_or(1);
        blocks
            .entry((element.dimension(), entity_tag, gmsh_type))
            .or_default()
            .push(*element);
    }

    writeln!(writer, "$Elements")?;
    if elements.is_empty() {
        le_u64(writer, 0)?;
        le_u64(writer, 0)?;
        le_u64(writer, 0)?;
        le_u64(writer, 0)?;
    } else {
        le_u64(writer, blocks.len() as u64)?;
        le_u64(writer, elements.len() as u64)?;
        le_u64(writer, min_element_tag)?;
        le_u64(writer, max_element_tag)?;
        for ((dimension, entity_tag, gmsh_type), block_elements) in blocks {
            le_i32(writer, i32::from(dimension))?;
            le_i32(writer, entity_tag)?;
            le_i32(writer, gmsh_type)?;
            le_u64(writer, block_elements.len() as u64)?;
            for element in block_elements {
                le_u64(writer, element.id)?;
                for node_id in &element.node_ids {
                    le_u64(writer, *node_id)?;
                }
            }
        }
    }
    writeln!(writer, "$EndElements")?;
    if let Some(info) = &mesh.partitioned {
        crate::msh_partitioned::write_ghost_elements_v4_binary(writer, &info.ghost_elements)?;
    }

    // Gmsh writes $Periodic after $Elements; empty links emit nothing.
    crate::msh_periodic::write_periodic_v4_binary(writer, &mesh.periodic)?;

    Ok(())
}

pub fn save_msh_v2_binary_to_path(path: &Path, mesh: &Mesh) -> Result<(), MshError> {
    let mut file = File::create(path)?;
    write_msh_v2_binary(&mut file, mesh)
}

pub fn save_msh_v4_binary_to_path(path: &Path, mesh: &Mesh) -> Result<(), MshError> {
    let mut file = File::create(path)?;
    write_msh_v4_binary(&mut file, mesh)
}

/// Parse a Gmsh MSH v4.1 BINARY file. Pinned against real Gmsh 4.14
/// output: only the section MARKER lines (`$Nodes`, ...) are ASCII — the
/// section headers, entity block headers and payloads are all
/// little-endian binary (`size_t` = u64, int = i32, coordinates = f64).
/// Node blocks list ALL tags first, then ALL coordinates. The block's
/// entity tag rides as the physical tag — the same "rem uses this as
/// physical tag" convention as the ASCII v4 reader.
fn parse_msh_v4_binary(bytes: &[u8]) -> Result<Mesh, MshError> {
    let mut mesh = Mesh::new();
    let mut pos = 0usize;

    let line = bin_read_line(bytes, &mut pos)?;
    if line != "$MeshFormat" {
        return Err(MshError::Parse {
            line: 0,
            message: format!("binary v4: missing $MeshFormat, found {line}"),
        });
    }
    let format_line = bin_read_line(bytes, &mut pos)?;
    let format_parts: Vec<&str> = format_line.split_whitespace().collect();
    if format_parts.len() < 3 || !format_parts[0].starts_with("4.") || format_parts[1] != "1" {
        return Err(MshError::Parse {
            line: 0,
            message: format!("binary v4: not a v4.1 binary format line: {format_line}"),
        });
    }
    read_binary_magic(bytes, &mut pos, "v4")?;
    bin_expect_marker(bytes, &mut pos, "$EndMeshFormat")?;

    loop {
        if pos >= bytes.len() {
            break;
        }
        let line = bin_read_line(bytes, &mut pos)?;
        if line.is_empty() {
            continue;
        }
        match line {
            "$PhysicalNames" => {
                let count_line = bin_read_line(bytes, &mut pos)?;
                let n_names: usize = count_line.trim().parse().unwrap_or(0);
                for _ in 0..n_names {
                    let name_line = bin_read_line(bytes, &mut pos)?;
                    let parts: Vec<&str> = name_line.splitn(3, ' ').collect();
                    if parts.len() >= 3 {
                        let dim = parts[0].parse::<i32>().unwrap_or(0);
                        let tag = parts[1].parse::<i32>().unwrap_or(0);
                        let name = parts[2].trim().trim_matches('"').to_string();
                        mesh.physical_names.insert((dim, tag), name);
                    }
                }
                bin_expect_marker(bytes, &mut pos, "$EndPhysicalNames")?;
            }
            "$Nodes" => {
                // The section header is four binary u64s in binary mode
                // (numBlocks, numNodes, minTag, maxTag) — NOT an ASCII
                // line. Real Gmsh writes it that way and refuses nothing:
                // an ASCII header makes Gmsh silently ingest an empty
                // model.
                let num_blocks = bin_read_u64(bytes, &mut pos, "nodes block count")? as usize;
                let _num_nodes = bin_read_u64(bytes, &mut pos, "nodes count")?;
                let _min_tag = bin_read_u64(bytes, &mut pos, "nodes min tag")?;
                let _max_tag = bin_read_u64(bytes, &mut pos, "nodes max tag")?;
                for _ in 0..num_blocks {
                    // Entity block header: i32 entityDim, i32 entityTag,
                    // i32 parametric, u64 numNodesInBlock — binary.
                    let _entity_dim = bin_read_i32(bytes, &mut pos, "node block dim")?;
                    let _entity_tag = bin_read_i32(bytes, &mut pos, "node block entity tag")?;
                    let parametric: i32 = bin_read_i32(bytes, &mut pos, "node block parametric")?;
                    if parametric != 0 {
                        return Err(MshError::Parse {
                            line: 0,
                            message: "binary v4: parametric node blocks are not supported".into(),
                        });
                    }
                    let count = bin_read_u64(bytes, &mut pos, "node block count")? as usize;
                    // The v4.1 layout lists ALL node tags first, then ALL
                    // coordinates — two runs, not per-node interleaved.
                    let mut ids = Vec::with_capacity(count);
                    for _ in 0..count {
                        ids.push(bin_read_u64(bytes, &mut pos, "node tag")?);
                    }
                    for id in ids {
                        let x = bin_read_f64(bytes, &mut pos, "node x")?;
                        let y = bin_read_f64(bytes, &mut pos, "node y")?;
                        let z = bin_read_f64(bytes, &mut pos, "node z")?;
                        mesh.add_node(Node::new(id, x, y, z));
                    }
                }
                bin_expect_marker(bytes, &mut pos, "$EndNodes")?;
            }
            "$Elements" => {
                let num_blocks = bin_read_u64(bytes, &mut pos, "element block count")? as usize;
                let _num_elements = bin_read_u64(bytes, &mut pos, "element count")?;
                let _min_tag = bin_read_u64(bytes, &mut pos, "element min tag")?;
                let _max_tag = bin_read_u64(bytes, &mut pos, "element max tag")?;
                for _ in 0..num_blocks {
                    // Entity block header: i32 entityDim, i32 entityTag,
                    // i32 elementType, u64 numElementsInBlock — binary.
                    let _entity_dim = bin_read_i32(bytes, &mut pos, "element block dim")?;
                    let entity_tag: i32 = bin_read_i32(bytes, &mut pos, "element block entity tag")?;
                    let element_type_id: i32 = bin_read_i32(bytes, &mut pos, "element type")?;
                    let count = bin_read_u64(bytes, &mut pos, "element block count")? as usize;
                    let etype = ElementType::from_gmsh_type_id(element_type_id);
                    let expected = etype.node_count();
                    if expected == 0 {
                        return Err(MshError::Parse {
                            line: 0,
                            message: format!("binary v4: unsupported element type {element_type_id}"),
                        });
                    }
                    for _ in 0..count {
                        let id = bin_read_u64(bytes, &mut pos, "element tag")?;
                        let mut node_ids = Vec::with_capacity(expected);
                        for _ in 0..expected {
                            node_ids.push(bin_read_u64(bytes, &mut pos, "element node")?);
                        }
                        let mut elem = Element::new(id, etype, node_ids);
                        elem.physical_tag = Some(entity_tag);
                        mesh.add_element(elem);
                    }
                }
                bin_expect_marker(bytes, &mut pos, "$EndElements")?;
            }
            "$Periodic" => {
                mesh.periodic = crate::msh_periodic::parse_periodic_v4_binary(bytes, &mut pos)?;
                bin_expect_marker(bytes, &mut pos, "$EndPeriodic")?;
            }
            "$Entities" => {
                // Real Gmsh always writes $Entities in v4.1; its binary
                // payload must be consumed structurally (a line scan would
                // hit raw binary). The content is discarded.
                crate::msh_periodic::skip_entities_v4_binary(bytes, &mut pos)?;
                bin_expect_marker(bytes, &mut pos, "$EndEntities")?;
            }
            "$PartitionedEntities" => {
                // The partition decomposition is RETAINED (unlike the
                // discarded $Entities payload): the retained records are
                // what the v4 writers re-emit for the round trip.
                let mut info = crate::msh_partitioned::parse_partitioned_v4_binary(
                    bytes,
                    &mut pos,
                )?;
                info.ghost_elements = mesh
                    .partitioned
                    .take()
                    .map(|p| p.ghost_elements)
                    .unwrap_or_default();
                mesh.partitioned = Some(info);
                bin_expect_marker(bytes, &mut pos, "$EndPartitionedEntities")?;
            }
            "$GhostElements" => {
                let cells =
                    crate::msh_partitioned::parse_ghost_elements_v4_binary(bytes, &mut pos)?;
                match mesh.partitioned.as_mut() {
                    Some(info) => info.ghost_elements = cells,
                    None => {
                        return Err(MshError::Parse {
                            line: 0,
                            message: "$GhostElements without $PartitionedEntities".into(),
                        })
                    }
                }
                bin_expect_marker(bytes, &mut pos, "$EndGhostElements")?;
            }
            _ if line.starts_with('$') && !line.starts_with("$End") => {
                let end_tag = format!("$End{}", &line[1..]);
                loop {
                    let skip_line = bin_read_line(bytes, &mut pos)?;
                    if skip_line == end_tag {
                        break;
                    }
                }
            }
            _ => {}
        }
    }

    Ok(mesh)
}

#[cfg(test)]
mod binary_tests {
    use super::*;
    use rmsh_model::{ElementType, Node};

    /// A mixed-dimension reference mesh: one line, one triangle, one tet,
    /// each with a distinct physical tag — enough to exercise block keying
    /// in every writer/reader pair.
    fn reference_mesh() -> Mesh {
        let mut mesh = Mesh::new();
        for (id, xyz) in [
            (1, [0.0, 0.0, 0.0]),
            (2, [1.0, 0.0, 0.0]),
            (3, [0.0, 1.0, 0.0]),
            (4, [0.0, 0.0, 1.0]),
        ] {
            mesh.add_node(Node::new(id, xyz[0], xyz[1], xyz[2]));
        }
        let mut line = Element::new(1, ElementType::Line2, vec![1, 2]);
        line.physical_tag = Some(7);
        mesh.add_element(line);
        let mut tri = Element::new(2, ElementType::Triangle3, vec![1, 2, 3]);
        tri.physical_tag = Some(8);
        mesh.add_element(tri);
        let mut tet = Element::new(3, ElementType::Tetrahedron4, vec![1, 2, 3, 4]);
        tet.physical_tag = Some(9);
        mesh.add_element(tet);
        mesh
    }

    fn assert_meshes_agree(a: &Mesh, b: &Mesh) {
        assert_eq!(a.node_count(), b.node_count());
        assert_eq!(a.element_count(), b.element_count());
        let mut anodes: Vec<_> = a.nodes.values().collect();
        anodes.sort_by_key(|n| n.id);
        let mut bnodes: Vec<_> = b.nodes.values().collect();
        bnodes.sort_by_key(|n| n.id);
        for (x, y) in anodes.iter().zip(bnodes.iter()) {
            assert_eq!(x.id, y.id);
            assert_eq!(x.position.x.to_bits(), y.position.x.to_bits());
            assert_eq!(x.position.y.to_bits(), y.position.y.to_bits());
            assert_eq!(x.position.z.to_bits(), y.position.z.to_bits());
        }
        let mut aelems: Vec<_> = a.elements.iter().collect();
        aelems.sort_by_key(|e| e.id);
        let mut belems: Vec<_> = b.elements.iter().collect();
        belems.sort_by_key(|e| e.id);
        for (x, y) in aelems.iter().zip(belems.iter()) {
            assert_eq!(x.id, y.id);
            assert_eq!(x.etype, y.etype);
            assert_eq!(x.node_ids, y.node_ids);
            assert_eq!(x.physical_tag, y.physical_tag);
        }
    }

    #[test]
    fn v2_binary_round_trips_field_exact() {
        let mesh = reference_mesh();
        let mut bytes = Vec::new();
        write_msh_v2_binary(&mut bytes, &mesh).expect("v2 binary write");
        assert!(bytes.starts_with(b"$MeshFormat\n2.2 1 8\n"));
        let loaded = parse_msh_v2_binary(&bytes).expect("v2 binary read");
        assert_meshes_agree(&mesh, &loaded);
    }

    #[test]
    fn v4_binary_round_trips_field_exact() {
        let mesh = reference_mesh();
        let mut bytes = Vec::new();
        write_msh_v4_binary(&mut bytes, &mesh).expect("v4 binary write");
        assert!(bytes.starts_with(b"$MeshFormat\n4.1 1 8\n"));
        let loaded = parse_msh_v4_binary(&bytes).expect("v4 binary read");
        assert_meshes_agree(&mesh, &loaded);
    }

    #[test]
    fn binary_detection_routes_both_flavors_and_ascii_still_works() {
        let mesh = reference_mesh();

        let mut v2 = Vec::new();
        write_msh_v2_binary(&mut v2, &mesh).expect("write v2 binary");
        let loaded = load_msh_from_bytes(&v2).expect("route v2 binary");
        assert_meshes_agree(&mesh, &loaded);

        let mut v4 = Vec::new();
        write_msh_v4_binary(&mut v4, &mesh).expect("write v4 binary");
        let loaded = load_msh_from_bytes(&v4).expect("route v4 binary");
        assert_meshes_agree(&mesh, &loaded);

        let mut ascii = Vec::new();
        write_msh_v4(&mut ascii, &mesh).expect("write v4 ascii");
        let loaded = load_msh_from_bytes(&ascii).expect("route v4 ascii");
        assert_meshes_agree(&mesh, &loaded);
    }

    #[test]
    fn ascii_and_binary_flavors_parse_to_the_same_mesh() {
        let mesh = reference_mesh();
        let mut ascii = Vec::new();
        write_msh_v4(&mut ascii, &mesh).expect("ascii write");
        let mut binary = Vec::new();
        write_msh_v4_binary(&mut binary, &mesh).expect("binary write");
        let from_ascii = load_msh_from_bytes(&ascii).expect("read ascii");
        let from_binary = load_msh_from_bytes(&binary).expect("read binary");
        assert_meshes_agree(&from_ascii, &from_binary);
    }

    #[test]
    fn binary_magic_corruption_is_refused_with_a_naming_error() {
        let mesh = reference_mesh();
        let mut bytes = Vec::new();
        write_msh_v2_binary(&mut bytes, &mesh).expect("write");
        // Corrupt the magic int (the 4 bytes right after the format line).
        let magic_at = b"$MeshFormat\n2.2 1 8\n".len();
        bytes[magic_at] = 0xFF;
        let err = parse_msh_v2_binary(&bytes).expect_err("corrupt magic must fail");
        let message = err.to_string();
        assert!(
            message.contains("endianness magic"),
            "the error must name the magic: {message}"
        );
    }

    #[test]
    fn truncated_binary_payload_is_refused_not_panicked() {
        let mesh = reference_mesh();
        let mut bytes = Vec::new();
        write_msh_v4_binary(&mut bytes, &mesh).expect("write");
        bytes.truncate(bytes.len() - 5);
        let err = parse_msh_v4_binary(&bytes).expect_err("truncated payload must fail");
        assert!(
            err.to_string().contains("unexpected end"),
            "the error must name the truncation: {err}"
        );
    }

    #[test]
    fn empty_mesh_binary_round_trips_both_flavors() {
        let mesh = Mesh::new();

        let mut v2 = Vec::new();
        write_msh_v2_binary(&mut v2, &mesh).expect("empty v2 write");
        let loaded = parse_msh_v2_binary(&v2).expect("empty v2 read");
        assert_eq!(loaded.node_count(), 0);
        assert_eq!(loaded.element_count(), 0);

        let mut v4 = Vec::new();
        write_msh_v4_binary(&mut v4, &mesh).expect("empty v4 write");
        let loaded = parse_msh_v4_binary(&v4).expect("empty v4 read");
        assert_eq!(loaded.node_count(), 0);
        assert_eq!(loaded.element_count(), 0);
    }
}
