// ---------------------------------------------------------------------------
// Gmsh MSH 4.1 Loader — Parse ASCII and Binary .msh mesh files
// ---------------------------------------------------------------------------
//
// Supports the Gmsh MSH 4.1 format:
//   - $MeshFormat (version, file-type, data-size)
//   - $PhysicalNames (dimension, tag, name)
//   - $Entities (volumes with physical group tags)
//   - $Nodes (per-entity blocks with coordinates)
//   - $Elements (per-entity blocks with connectivity)
//
// EMStudio conventions for PhysicalNames:
//   - "mat:<material>"  → volume material assignment
//   - "bc:<boundary>"   → boundary condition surface
//   - "ns:<selection>"  → named selection (face/edge)

use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use thiserror::Error;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug, Error)]
pub enum MeshError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("Parse error at line {line}: {message}")]
    Parse { line: usize, message: String },
    #[error("Unexpected section: expected {expected}, found {found}")]
    UnexpectedSection { expected: String, found: String },
    #[error("Unsupported MSH version: {0}")]
    UnsupportedVersion(String),
    #[error("Invalid binary data: {0}")]
    InvalidBinary(String),
}

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

/// A parsed Gmsh MSH 4.1 mesh.
#[derive(Debug, Clone)]
pub struct MshMesh {
    pub version: String,
    pub binary: bool,
    pub physical_names: Vec<PhysicalName>,
    pub entities: Vec<MshEntity>,
    pub nodes: Vec<MshNode>,
    pub elements: Vec<MshElement>,
    /// Map from node tag to index in `nodes` vec (for fast lookup).
    node_index: HashMap<u64, usize>,
}

#[derive(Debug, Clone)]
pub struct PhysicalName {
    pub dimension: u32,
    pub tag: u32,
    pub name: String,
}

#[derive(Debug, Clone)]
pub struct MshEntity {
    pub dimension: u32,
    pub tag: i32,
    pub physical_tags: Vec<i32>,
}

#[derive(Debug, Clone, Copy)]
pub struct MshNode {
    pub tag: u64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone)]
pub struct MshElement {
    pub tag: u64,
    /// Gmsh element type: 1=line2, 2=tri3, 4=tet4, 11=tet10, etc.
    pub element_type: u32,
    /// Entity dimension this element belongs to.
    pub entity_dim: u32,
    /// Entity tag this element belongs to.
    pub entity_tag: i32,
    pub node_tags: Vec<u64>,
}

/// Common Gmsh element type codes.
pub mod element_types {
    pub const LINE2: u32 = 1;
    pub const TRI3: u32 = 2;
    pub const QUAD4: u32 = 3;
    pub const TET4: u32 = 4;
    pub const HEX8: u32 = 5;
    pub const TET10: u32 = 11;
    pub const TRI6: u32 = 9;

    /// Number of nodes for a given element type.
    pub fn num_nodes(element_type: u32) -> Option<usize> {
        match element_type {
            1 => Some(2),   // line2
            2 => Some(3),   // tri3
            3 => Some(4),   // quad4
            4 => Some(4),   // tet4
            5 => Some(8),   // hex8
            9 => Some(6),   // tri6
            11 => Some(10), // tet10
            15 => Some(1),  // point
            _ => None,
        }
    }
}

impl MshMesh {
    /// Load a .msh file (auto-detects ASCII/Binary).
    pub fn load(path: &Path) -> Result<Self, MeshError> {
        let file = std::fs::File::open(path)?;
        let mut reader = BufReader::new(file);
        Self::read_from(&mut reader)
    }

    /// Parse from a reader.
    pub fn read_from<R: Read + Seek>(reader: &mut BufReader<R>) -> Result<Self, MeshError> {
        let mut mesh = MshMesh {
            version: String::new(),
            binary: false,
            physical_names: Vec::new(),
            entities: Vec::new(),
            nodes: Vec::new(),
            elements: Vec::new(),
            node_index: HashMap::new(),
        };

        let mut line_num = 0usize;

        // Read $MeshFormat first
        let line = read_line(reader, &mut line_num)?;
        if line.trim() != "$MeshFormat" {
            return Err(MeshError::UnexpectedSection {
                expected: "$MeshFormat".into(),
                found: line.trim().into(),
            });
        }
        let fmt_line = read_line(reader, &mut line_num)?;
        let parts: Vec<&str> = fmt_line.trim().split_whitespace().collect();
        if parts.len() < 3 {
            return Err(MeshError::Parse {
                line: line_num,
                message: "MeshFormat line needs 3 fields".into(),
            });
        }
        mesh.version = parts[0].to_string();
        if !mesh.version.starts_with("4.") {
            return Err(MeshError::UnsupportedVersion(mesh.version.clone()));
        }
        let file_type: u32 = parts[1].parse().map_err(|_| MeshError::Parse {
            line: line_num,
            message: "invalid file-type".into(),
        })?;
        mesh.binary = file_type == 1;
        let _data_size: u32 = parts[2].parse().map_err(|_| MeshError::Parse {
            line: line_num,
            message: "invalid data-size".into(),
        })?;

        if mesh.binary {
            // In binary mode, after the format line there's a newline, then 4 bytes endian marker
            // We need to read a newline then the binary int
            let mut endian_buf = [0u8; 4];
            reader.read_exact(&mut endian_buf)?;
            let endian_val = u32::from_le_bytes(endian_buf);
            if endian_val != 1 {
                return Err(MeshError::InvalidBinary(
                    "endian marker mismatch (expected little-endian 0x00000001)".into(),
                ));
            }
            // Read the trailing newline after the endian marker
            let mut nl = [0u8; 1];
            reader.read_exact(&mut nl)?;
        }

        let end_fmt = read_line(reader, &mut line_num)?;
        if end_fmt.trim() != "$EndMeshFormat" {
            return Err(MeshError::UnexpectedSection {
                expected: "$EndMeshFormat".into(),
                found: end_fmt.trim().into(),
            });
        }

        // Read remaining sections in any order
        loop {
            let line = match read_line(reader, &mut line_num) {
                Ok(l) => l,
                Err(MeshError::Io(ref e)) if e.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e),
            };
            let section = line.trim().to_string();
            if section.is_empty() {
                continue;
            }

            match section.as_str() {
                "$PhysicalNames" => {
                    mesh.physical_names =
                        read_physical_names(reader, &mut line_num)?;
                }
                "$Entities" => {
                    mesh.entities = read_entities(reader, &mut line_num)?;
                }
                "$Nodes" => {
                    if mesh.binary {
                        mesh.nodes = read_nodes_binary(reader, &mut line_num)?;
                    } else {
                        mesh.nodes = read_nodes_ascii(reader, &mut line_num)?;
                    }
                }
                "$Elements" => {
                    if mesh.binary {
                        mesh.elements = read_elements_binary(reader, &mut line_num)?;
                    } else {
                        mesh.elements = read_elements_ascii(reader, &mut line_num)?;
                    }
                }
                _ => {
                    // Skip unknown sections
                    let end_tag = format!("$End{}", &section[1..]);
                    loop {
                        let skip = read_line(reader, &mut line_num)?;
                        if skip.trim() == end_tag {
                            break;
                        }
                    }
                }
            }
        }

        // Build node index
        for (i, node) in mesh.nodes.iter().enumerate() {
            mesh.node_index.insert(node.tag, i);
        }

        Ok(mesh)
    }

    /// Get node by tag.
    pub fn node_by_tag(&self, tag: u64) -> Option<&MshNode> {
        self.node_index.get(&tag).map(|&i| &self.nodes[i])
    }

    /// Get node position as [f64; 3].
    pub fn node_position(&self, tag: u64) -> Option<[f64; 3]> {
        self.node_by_tag(tag).map(|n| [n.x, n.y, n.z])
    }

    /// Filter elements by physical name prefix ("mat:", "bc:", "ns:").
    pub fn elements_by_physical(&self, prefix: &str) -> Vec<&MshElement> {
        // Build entity_tag → physical_tag mapping
        let mut entity_to_phys: HashMap<(u32, i32), Vec<i32>> = HashMap::new();
        for ent in &self.entities {
            entity_to_phys.insert((ent.dimension, ent.tag), ent.physical_tags.clone());
        }

        // Collect physical tags that match prefix
        let matching_tags: Vec<u32> = self
            .physical_names
            .iter()
            .filter(|pn| pn.name.starts_with(prefix))
            .map(|pn| pn.tag)
            .collect();

        self.elements
            .iter()
            .filter(|el| {
                if let Some(phys_tags) = entity_to_phys.get(&(el.entity_dim, el.entity_tag)) {
                    phys_tags
                        .iter()
                        .any(|&pt| matching_tags.contains(&(pt as u32)))
                } else {
                    false
                }
            })
            .collect()
    }

    /// Get all tetrahedral elements (tet4 + tet10).
    pub fn tetrahedra(&self) -> Vec<&MshElement> {
        self.elements
            .iter()
            .filter(|e| e.element_type == element_types::TET4 || e.element_type == element_types::TET10)
            .collect()
    }

    /// Get all triangle elements (tri3 + tri6).
    pub fn triangles(&self) -> Vec<&MshElement> {
        self.elements
            .iter()
            .filter(|e| e.element_type == element_types::TRI3 || e.element_type == element_types::TRI6)
            .collect()
    }

    /// Total number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Total number of elements.
    pub fn num_elements(&self) -> usize {
        self.elements.len()
    }
}

// ---------------------------------------------------------------------------
// ASCII section readers
// ---------------------------------------------------------------------------

fn read_line<R: Read>(reader: &mut BufReader<R>, line_num: &mut usize) -> Result<String, MeshError> {
    let mut line = String::new();
    let bytes = reader.read_line(&mut line)?;
    if bytes == 0 {
        return Err(MeshError::Io(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "unexpected end of file",
        )));
    }
    *line_num += 1;
    Ok(line)
}

fn read_physical_names<R: Read>(
    reader: &mut BufReader<R>,
    line_num: &mut usize,
) -> Result<Vec<PhysicalName>, MeshError> {
    let count_line = read_line(reader, line_num)?;
    let count: usize = count_line.trim().parse().map_err(|_| MeshError::Parse {
        line: *line_num,
        message: "invalid physical names count".into(),
    })?;

    let mut names = Vec::with_capacity(count);
    for _ in 0..count {
        let line = read_line(reader, line_num)?;
        let parts: Vec<&str> = line.trim().splitn(3, ' ').collect();
        if parts.len() < 3 {
            return Err(MeshError::Parse {
                line: *line_num,
                message: "physical name needs 3 fields".into(),
            });
        }
        let dimension: u32 = parts[0].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid dimension".into(),
        })?;
        let tag: u32 = parts[1].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid tag".into(),
        })?;
        let name = parts[2].trim_matches('"').to_string();
        names.push(PhysicalName {
            dimension,
            tag,
            name,
        });
    }

    let end = read_line(reader, line_num)?;
    if end.trim() != "$EndPhysicalNames" {
        return Err(MeshError::UnexpectedSection {
            expected: "$EndPhysicalNames".into(),
            found: end.trim().into(),
        });
    }
    Ok(names)
}

fn read_entities<R: Read>(
    reader: &mut BufReader<R>,
    line_num: &mut usize,
) -> Result<Vec<MshEntity>, MeshError> {
    let header = read_line(reader, line_num)?;
    let counts: Vec<usize> = header
        .trim()
        .split_whitespace()
        .map(|s| s.parse::<usize>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid entity counts".into(),
        })?;

    // counts = [numPoints, numCurves, numSurfaces, numVolumes]
    let (num_points, num_curves, num_surfaces, num_volumes) = match counts.len() {
        4 => (counts[0], counts[1], counts[2], counts[3]),
        _ => {
            return Err(MeshError::Parse {
                line: *line_num,
                message: "entity header needs 4 counts".into(),
            })
        }
    };

    let mut entities = Vec::new();

    // Skip point entities (dim=0)
    for _ in 0..num_points {
        let _line = read_line(reader, line_num)?;
    }

    // Read curve entities (dim=1)
    for _ in 0..num_curves {
        let line = read_line(reader, line_num)?;
        if let Some(ent) = parse_entity_line(&line, 1) {
            entities.push(ent);
        }
    }

    // Read surface entities (dim=2)
    for _ in 0..num_surfaces {
        let line = read_line(reader, line_num)?;
        if let Some(ent) = parse_entity_line(&line, 2) {
            entities.push(ent);
        }
    }

    // Read volume entities (dim=3)
    for _ in 0..num_volumes {
        let line = read_line(reader, line_num)?;
        if let Some(ent) = parse_entity_line(&line, 3) {
            entities.push(ent);
        }
    }

    let end = read_line(reader, line_num)?;
    if end.trim() != "$EndEntities" {
        return Err(MeshError::UnexpectedSection {
            expected: "$EndEntities".into(),
            found: end.trim().into(),
        });
    }
    Ok(entities)
}

fn parse_entity_line(line: &str, dimension: u32) -> Option<MshEntity> {
    let parts: Vec<&str> = line.trim().split_whitespace().collect();
    if parts.is_empty() {
        return None;
    }
    let tag: i32 = parts[0].parse().ok()?;

    // Entity line format:
    //   point (dim=0): tag x y z numPhysicalTags [physicalTags...]
    //   curve/surf/vol (dim>=1): tag minX minY minZ maxX maxY maxZ numPhysicalTags [physicalTags...] numBounding [boundingTags...]
    // Index of numPhysicalTags field:
    let phys_idx = if dimension == 0 { 4 } else { 7 }; // 1(tag) + 3(xyz) or 1(tag) + 6(bbox)
    if parts.len() <= phys_idx {
        return Some(MshEntity {
            dimension,
            tag,
            physical_tags: Vec::new(),
        });
    }

    let num_phys: usize = parts[phys_idx].parse().ok()?;
    let mut physical_tags = Vec::with_capacity(num_phys);
    for i in 0..num_phys {
        if let Some(pt) = parts.get(phys_idx + 1 + i) {
            if let Ok(v) = pt.parse::<i32>() {
                physical_tags.push(v);
            }
        }
    }

    Some(MshEntity {
        dimension,
        tag,
        physical_tags,
    })
}

fn read_nodes_ascii<R: Read>(
    reader: &mut BufReader<R>,
    line_num: &mut usize,
) -> Result<Vec<MshNode>, MeshError> {
    let header = read_line(reader, line_num)?;
    let parts: Vec<&str> = header.trim().split_whitespace().collect();
    // numEntityBlocks numNodes minNodeTag maxNodeTag
    if parts.len() < 4 {
        return Err(MeshError::Parse {
            line: *line_num,
            message: "nodes header needs 4 fields".into(),
        });
    }
    let num_blocks: usize = parts[0].parse().map_err(|_| MeshError::Parse {
        line: *line_num,
        message: "invalid numEntityBlocks".into(),
    })?;
    let total_nodes: usize = parts[1].parse().map_err(|_| MeshError::Parse {
        line: *line_num,
        message: "invalid numNodes".into(),
    })?;

    let mut nodes = Vec::with_capacity(total_nodes);

    for _ in 0..num_blocks {
        let block_header = read_line(reader, line_num)?;
        let bp: Vec<&str> = block_header.trim().split_whitespace().collect();
        // entityDim entityTag parametric numNodesInBlock
        if bp.len() < 4 {
            return Err(MeshError::Parse {
                line: *line_num,
                message: "node block header needs 4 fields".into(),
            });
        }
        let num_in_block: usize = bp[3].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid numNodesInBlock".into(),
        })?;

        // First, read all node tags
        let mut tags = Vec::with_capacity(num_in_block);
        for _ in 0..num_in_block {
            let tag_line = read_line(reader, line_num)?;
            let tag: u64 = tag_line.trim().parse().map_err(|_| MeshError::Parse {
                line: *line_num,
                message: "invalid node tag".into(),
            })?;
            tags.push(tag);
        }

        // Then, read all coordinates
        for tag in tags {
            let coord_line = read_line(reader, line_num)?;
            let cp: Vec<&str> = coord_line.trim().split_whitespace().collect();
            if cp.len() < 3 {
                return Err(MeshError::Parse {
                    line: *line_num,
                    message: "node needs 3 coordinates".into(),
                });
            }
            let x: f64 = cp[0].parse().map_err(|_| MeshError::Parse {
                line: *line_num,
                message: "invalid x".into(),
            })?;
            let y: f64 = cp[1].parse().map_err(|_| MeshError::Parse {
                line: *line_num,
                message: "invalid y".into(),
            })?;
            let z: f64 = cp[2].parse().map_err(|_| MeshError::Parse {
                line: *line_num,
                message: "invalid z".into(),
            })?;
            nodes.push(MshNode { tag, x, y, z });
        }
    }

    let end = read_line(reader, line_num)?;
    if end.trim() != "$EndNodes" {
        return Err(MeshError::UnexpectedSection {
            expected: "$EndNodes".into(),
            found: end.trim().into(),
        });
    }
    Ok(nodes)
}

fn read_elements_ascii<R: Read>(
    reader: &mut BufReader<R>,
    line_num: &mut usize,
) -> Result<Vec<MshElement>, MeshError> {
    let header = read_line(reader, line_num)?;
    let parts: Vec<&str> = header.trim().split_whitespace().collect();
    if parts.len() < 4 {
        return Err(MeshError::Parse {
            line: *line_num,
            message: "elements header needs 4 fields".into(),
        });
    }
    let num_blocks: usize = parts[0].parse().map_err(|_| MeshError::Parse {
        line: *line_num,
        message: "invalid numEntityBlocks".into(),
    })?;
    let total_elements: usize = parts[1].parse().map_err(|_| MeshError::Parse {
        line: *line_num,
        message: "invalid numElements".into(),
    })?;

    let mut elements = Vec::with_capacity(total_elements);

    for _ in 0..num_blocks {
        let block_header = read_line(reader, line_num)?;
        let bp: Vec<&str> = block_header.trim().split_whitespace().collect();
        // entityDim entityTag elementType numElementsInBlock
        if bp.len() < 4 {
            return Err(MeshError::Parse {
                line: *line_num,
                message: "element block header needs 4 fields".into(),
            });
        }
        let entity_dim: u32 = bp[0].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid entityDim".into(),
        })?;
        let entity_tag: i32 = bp[1].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid entityTag".into(),
        })?;
        let element_type: u32 = bp[2].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid elementType".into(),
        })?;
        let num_in_block: usize = bp[3].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid numElementsInBlock".into(),
        })?;

        let nodes_per_elem = element_types::num_nodes(element_type).unwrap_or(0);

        for _ in 0..num_in_block {
            let el_line = read_line(reader, line_num)?;
            let ep: Vec<&str> = el_line.trim().split_whitespace().collect();
            if ep.is_empty() {
                continue;
            }
            let tag: u64 = ep[0].parse().map_err(|_| MeshError::Parse {
                line: *line_num,
                message: "invalid element tag".into(),
            })?;

            let mut node_tags = Vec::with_capacity(nodes_per_elem);
            for i in 1..=nodes_per_elem {
                if let Some(s) = ep.get(i) {
                    let nt: u64 = s.parse().map_err(|_| MeshError::Parse {
                        line: *line_num,
                        message: "invalid node tag in element".into(),
                    })?;
                    node_tags.push(nt);
                }
            }

            elements.push(MshElement {
                tag,
                element_type,
                entity_dim,
                entity_tag,
                node_tags,
            });
        }
    }

    let end = read_line(reader, line_num)?;
    if end.trim() != "$EndElements" {
        return Err(MeshError::UnexpectedSection {
            expected: "$EndElements".into(),
            found: end.trim().into(),
        });
    }
    Ok(elements)
}

// ---------------------------------------------------------------------------
// Binary section readers
// ---------------------------------------------------------------------------

fn read_nodes_binary<R: Read + Seek>(
    reader: &mut BufReader<R>,
    line_num: &mut usize,
) -> Result<Vec<MshNode>, MeshError> {
    let header = read_line(reader, line_num)?;
    let parts: Vec<&str> = header.trim().split_whitespace().collect();
    if parts.len() < 4 {
        return Err(MeshError::Parse {
            line: *line_num,
            message: "nodes header needs 4 fields".into(),
        });
    }
    let num_blocks: usize = parts[0].parse().map_err(|_| MeshError::Parse {
        line: *line_num,
        message: "invalid numEntityBlocks".into(),
    })?;
    let total_nodes: usize = parts[1].parse().map_err(|_| MeshError::Parse {
        line: *line_num,
        message: "invalid numNodes".into(),
    })?;

    let mut nodes = Vec::with_capacity(total_nodes);

    for _ in 0..num_blocks {
        // Block header is ASCII: entityDim entityTag parametric numNodesInBlock\n
        let block_header = read_line(reader, line_num)?;
        let bp: Vec<&str> = block_header.trim().split_whitespace().collect();
        if bp.len() < 4 {
            return Err(MeshError::Parse {
                line: *line_num,
                message: "node block header needs 4 fields".into(),
            });
        }
        let num_in_block: usize = bp[3].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid numNodesInBlock".into(),
        })?;

        // Binary: read node tags (size_t = u64 each)
        let mut tags = Vec::with_capacity(num_in_block);
        for _ in 0..num_in_block {
            let mut buf = [0u8; 8];
            reader.read_exact(&mut buf)?;
            tags.push(u64::from_le_bytes(buf));
        }

        // Binary: read coordinates (3 × f64 each = 24 bytes per node)
        for tag in tags {
            let mut buf = [0u8; 24];
            reader.read_exact(&mut buf)?;
            let x = f64::from_le_bytes(buf[0..8].try_into().unwrap());
            let y = f64::from_le_bytes(buf[8..16].try_into().unwrap());
            let z = f64::from_le_bytes(buf[16..24].try_into().unwrap());
            nodes.push(MshNode { tag, x, y, z });
        }
    }

    // Read trailing newline + $EndNodes
    let end = read_line(reader, line_num)?;
    let end_trimmed = end.trim();
    if end_trimmed != "$EndNodes" {
        // The newline after binary data might be consumed; try next line
        if !end_trimmed.is_empty() && end_trimmed != "$EndNodes" {
            let end2 = read_line(reader, line_num)?;
            if end2.trim() != "$EndNodes" {
                return Err(MeshError::UnexpectedSection {
                    expected: "$EndNodes".into(),
                    found: end2.trim().into(),
                });
            }
        } else {
            let end2 = read_line(reader, line_num)?;
            if end2.trim() != "$EndNodes" {
                return Err(MeshError::UnexpectedSection {
                    expected: "$EndNodes".into(),
                    found: end2.trim().into(),
                });
            }
        }
    }

    Ok(nodes)
}

fn read_elements_binary<R: Read + Seek>(
    reader: &mut BufReader<R>,
    line_num: &mut usize,
) -> Result<Vec<MshElement>, MeshError> {
    let header = read_line(reader, line_num)?;
    let parts: Vec<&str> = header.trim().split_whitespace().collect();
    if parts.len() < 4 {
        return Err(MeshError::Parse {
            line: *line_num,
            message: "elements header needs 4 fields".into(),
        });
    }
    let num_blocks: usize = parts[0].parse().map_err(|_| MeshError::Parse {
        line: *line_num,
        message: "invalid numEntityBlocks".into(),
    })?;
    let total_elements: usize = parts[1].parse().map_err(|_| MeshError::Parse {
        line: *line_num,
        message: "invalid numElements".into(),
    })?;

    let mut elements = Vec::with_capacity(total_elements);

    for _ in 0..num_blocks {
        // Block header is ASCII
        let block_header = read_line(reader, line_num)?;
        let bp: Vec<&str> = block_header.trim().split_whitespace().collect();
        if bp.len() < 4 {
            return Err(MeshError::Parse {
                line: *line_num,
                message: "element block header needs 4 fields".into(),
            });
        }
        let entity_dim: u32 = bp[0].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid entityDim".into(),
        })?;
        let entity_tag: i32 = bp[1].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid entityTag".into(),
        })?;
        let element_type: u32 = bp[2].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid elementType".into(),
        })?;
        let num_in_block: usize = bp[3].parse().map_err(|_| MeshError::Parse {
            line: *line_num,
            message: "invalid numElementsInBlock".into(),
        })?;

        let nodes_per_elem = element_types::num_nodes(element_type).unwrap_or(0);
        // Binary: each element is (elementTag + node_tags) as size_t (u64)
        let vals_per_elem = 1 + nodes_per_elem; // tag + node_tags

        for _ in 0..num_in_block {
            let mut buf = vec![0u8; vals_per_elem * 8];
            reader.read_exact(&mut buf)?;

            let tag = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let mut node_tags = Vec::with_capacity(nodes_per_elem);
            for i in 0..nodes_per_elem {
                let offset = (1 + i) * 8;
                let nt = u64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap());
                node_tags.push(nt);
            }

            elements.push(MshElement {
                tag,
                element_type,
                entity_dim,
                entity_tag,
                node_tags,
            });
        }
    }

    // Read trailing newline + $EndElements
    let end = read_line(reader, line_num)?;
    let end_trimmed = end.trim();
    if end_trimmed != "$EndElements" {
        if !end_trimmed.is_empty() {
            return Err(MeshError::UnexpectedSection {
                expected: "$EndElements".into(),
                found: end_trimmed.into(),
            });
        }
        let end2 = read_line(reader, line_num)?;
        if end2.trim() != "$EndElements" {
            return Err(MeshError::UnexpectedSection {
                expected: "$EndElements".into(),
                found: end2.trim().into(),
            });
        }
    }

    Ok(elements)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn sample_ascii_msh() -> &'static str {
        "$MeshFormat\n\
         4.1 0 8\n\
         $EndMeshFormat\n\
         $PhysicalNames\n\
         3\n\
         3 1 \"mat:vacuum\"\n\
         3 2 \"mat:copper\"\n\
         2 3 \"bc:PEC_GND\"\n\
         $EndPhysicalNames\n\
         $Entities\n\
         0 0 0 2\n\
         1 0.0 0.0 0.0 10.0 10.0 10.0 1 1 0\n\
         2 2.0 2.0 2.0 8.0 8.0 8.0 1 2 0\n\
         $EndEntities\n\
         $Nodes\n\
         2 8 1 8\n\
         3 1 0 4\n\
         1\n\
         2\n\
         3\n\
         4\n\
         0.0 0.0 0.0\n\
         10.0 0.0 0.0\n\
         5.0 10.0 0.0\n\
         5.0 5.0 10.0\n\
         3 2 0 4\n\
         5\n\
         6\n\
         7\n\
         8\n\
         2.0 2.0 2.0\n\
         8.0 2.0 2.0\n\
         5.0 8.0 2.0\n\
         5.0 5.0 8.0\n\
         $EndNodes\n\
         $Elements\n\
         2 2 1 2\n\
         3 1 4 1\n\
         1 1 2 3 4\n\
         3 2 4 1\n\
         2 5 6 7 8\n\
         $EndElements\n"
    }

    #[test]
    fn parse_ascii_mesh() {
        let data = sample_ascii_msh();
        let cursor = Cursor::new(data.as_bytes());
        let mut reader = BufReader::new(cursor);
        let mesh = MshMesh::read_from(&mut reader).unwrap();

        assert_eq!(mesh.version, "4.1");
        assert!(!mesh.binary);
        assert_eq!(mesh.physical_names.len(), 3);
        assert_eq!(mesh.physical_names[0].name, "mat:vacuum");
        assert_eq!(mesh.physical_names[0].dimension, 3);
        assert_eq!(mesh.physical_names[1].name, "mat:copper");
        assert_eq!(mesh.physical_names[2].name, "bc:PEC_GND");
        assert_eq!(mesh.nodes.len(), 8);
        assert_eq!(mesh.elements.len(), 2);

        // Check node coordinates
        let n1 = mesh.node_by_tag(1).unwrap();
        assert!((n1.x - 0.0).abs() < 1e-10);
        assert!((n1.y - 0.0).abs() < 1e-10);

        let n4 = mesh.node_by_tag(4).unwrap();
        assert!((n4.z - 10.0).abs() < 1e-10);

        // Check element connectivity
        assert_eq!(mesh.elements[0].element_type, element_types::TET4);
        assert_eq!(mesh.elements[0].node_tags, vec![1, 2, 3, 4]);
        assert_eq!(mesh.elements[1].node_tags, vec![5, 6, 7, 8]);
    }

    #[test]
    fn tetrahedra_filter() {
        let data = sample_ascii_msh();
        let cursor = Cursor::new(data.as_bytes());
        let mut reader = BufReader::new(cursor);
        let mesh = MshMesh::read_from(&mut reader).unwrap();

        let tets = mesh.tetrahedra();
        assert_eq!(tets.len(), 2);
    }

    #[test]
    fn physical_name_filter() {
        let data = sample_ascii_msh();
        let cursor = Cursor::new(data.as_bytes());
        let mut reader = BufReader::new(cursor);
        let mesh = MshMesh::read_from(&mut reader).unwrap();

        let mat_elements = mesh.elements_by_physical("mat:");
        assert_eq!(mat_elements.len(), 2); // both tets belong to mat: volumes

        let bc_elements = mesh.elements_by_physical("bc:");
        assert_eq!(bc_elements.len(), 0); // no surface elements defined
    }

    #[test]
    fn node_index_lookup() {
        let data = sample_ascii_msh();
        let cursor = Cursor::new(data.as_bytes());
        let mut reader = BufReader::new(cursor);
        let mesh = MshMesh::read_from(&mut reader).unwrap();

        assert!(mesh.node_by_tag(1).is_some());
        assert!(mesh.node_by_tag(8).is_some());
        assert!(mesh.node_by_tag(99).is_none());

        let pos = mesh.node_position(2).unwrap();
        assert!((pos[0] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn entity_physical_tags() {
        let data = sample_ascii_msh();
        let cursor = Cursor::new(data.as_bytes());
        let mut reader = BufReader::new(cursor);
        let mesh = MshMesh::read_from(&mut reader).unwrap();

        assert_eq!(mesh.entities.len(), 2);
        assert_eq!(mesh.entities[0].dimension, 3);
        assert_eq!(mesh.entities[0].tag, 1);
        assert_eq!(mesh.entities[0].physical_tags, vec![1]);
        assert_eq!(mesh.entities[1].tag, 2);
        assert_eq!(mesh.entities[1].physical_tags, vec![2]);
    }

    #[test]
    fn element_types_num_nodes() {
        assert_eq!(element_types::num_nodes(element_types::TET4), Some(4));
        assert_eq!(element_types::num_nodes(element_types::TET10), Some(10));
        assert_eq!(element_types::num_nodes(element_types::TRI3), Some(3));
        assert_eq!(element_types::num_nodes(element_types::LINE2), Some(2));
        assert_eq!(element_types::num_nodes(99), None);
    }
}
