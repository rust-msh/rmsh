//! VTK/VTU unstructured-grid I/O (the rmsh format layer's VTK flavor).
//!
//! Two flavors, both ASCII in v1 (the binary legacy flavor and the
//! appended / base64 VTU payloads are documented follow-ups):
//! - legacy VTK (.vtk): `ASCII` `UNSTRUCTURED_GRID`;
//! - VTU (.vtu): inline `format="ascii"` DataArrays.
//!
//! The heart of the module is the cell-type mapping between the Gmsh
//! element space (rmsh's native [`ElementType`], Gmsh node ordering) and
//! the VTK cell types:
//! - the VTK type ids and the five node-order permutations are pinned to
//!   the ecosystem reference (meshio `_vtk_common` / `gmsh/common` — the
//!   same tables the ParaView ecosystem consumes), with the VTK wedge
//!   normal convention handled by the wedge6 flip `[0, 2, 1, 3, 5, 4]`;
//!   only five types permute at all (tetra10, hex20, hex27, wedge15,
//!   pyramid13) — for everything else Gmsh and VTK order identically;
//! - Gmsh elements with no VTK counterpart (third-order triangles/quads,
//!   Tet20, Hex64, Pyramid14, Triangle9, Quad12/16) are REFUSED by name
//!   on write and read — never silently downgraded;
//! - the read-side permutation is computed as the inverse of the
//!   write-side permutation (never transcribed by hand).
//!
//! Physical groups ride as an integer cell scalar named `physical` (one
//! tag per cell, 0 = untagged; the VTK flavors carry no physical-name
//! table in v1).

use std::io::Write;

use rmsh_model::{Element, ElementType, Mesh};

use crate::MshError;

// ─── The Gmsh <-> VTK cell mapping ────────────────────────────────────────────

/// The VTK cell type id for a Gmsh-ordered element type, plus the
/// node-order permutation: `vtk_nodes[j] = gmsh_nodes[perm[j]]` (empty =
/// identity). `None` for elements without a VTK counterpart (the error
/// path names them).
fn vtk_cell(etype: ElementType) -> Option<(i32, &'static [usize])> {
    const IDENTITY: &[usize] = &[];
    match etype {
        ElementType::Point1 => Some((1, IDENTITY)),
        ElementType::Line2 => Some((3, IDENTITY)),
        ElementType::Line3 => Some((21, IDENTITY)),
        ElementType::Line4 => Some((35, IDENTITY)),
        ElementType::Triangle3 => Some((5, IDENTITY)),
        ElementType::Triangle6 => Some((22, IDENTITY)),
        ElementType::Quad4 => Some((9, IDENTITY)),
        ElementType::Quad8 => Some((23, IDENTITY)),
        ElementType::Quad9 => Some((28, IDENTITY)),
        ElementType::Tetrahedron4 => Some((10, IDENTITY)),
        ElementType::Tetrahedron10 => Some((24, &[0, 1, 2, 3, 4, 5, 6, 7, 9, 8])),
        ElementType::Hexahedron8 => Some((12, IDENTITY)),
        ElementType::Hexahedron20 => Some((
            25,
            &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15,
            ],
        )),
        ElementType::Hexahedron27 => Some((
            29,
            &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15,
                22, 23, 21, 24, 20, 25, 26,
            ],
        )),
        ElementType::Prism6 => Some((13, &[0, 2, 1, 3, 5, 4])),
        ElementType::Prism15 => Some((
            26,
            &[0, 1, 2, 3, 4, 5, 6, 9, 7, 12, 14, 13, 8, 10, 11],
        )),
        ElementType::Pyramid5 => Some((14, IDENTITY)),
        ElementType::Pyramid13 => Some((27, &[0, 1, 2, 3, 4, 5, 8, 10, 6, 7, 9, 11, 12])),
        _ => None,
    }
}

/// Invert a permutation (`inv[value] = index`).
fn invert(perm: &[usize]) -> Vec<usize> {
    let mut inv = vec![0; perm.len()];
    for (i, &v) in perm.iter().enumerate() {
        inv[v] = i;
    }
    inv
}

/// Apply a node permutation: `out[j] = nodes[perm[j]]`, or `nodes[j]`
/// when the permutation is the identity (empty).
fn permute_nodes(perm: &[usize], nodes: &[u64]) -> Vec<u64> {
    if perm.is_empty() {
        nodes.to_vec()
    } else {
        perm.iter().map(|&j| nodes[j]).collect()
    }
}

fn vtk_type_name(vtk_type: i32) -> &'static str {
    match vtk_type {
        1 => "VTK_VERTEX",
        3 => "VTK_LINE",
        5 => "VTK_TRIANGLE",
        9 => "VTK_QUAD",
        10 => "VTK_TETRA",
        12 => "VTK_HEXAHEDRON",
        13 => "VTK_WEDGE",
        14 => "VTK_PYRAMID",
        21 => "VTK_QUADRATIC_EDGE",
        22 => "VTK_QUADRATIC_TRIANGLE",
        23 => "VTK_QUADRATIC_QUAD",
        24 => "VTK_QUADRATIC_TETRA",
        25 => "VTK_QUADRATIC_HEXAHEDRON",
        26 => "VTK_QUADRATIC_WEDGE",
        27 => "VTK_QUADRATIC_PYRAMID",
        28 => "VTK_BIQUADRATIC_QUAD",
        29 => "VTK_TRIQUADRATIC_HEXAHEDRON",
        35 => "VTK_CUBIC_LINE",
        _ => "VTK_UNKNOWN",
    }
}

// ─── Legacy VTK (.vtk, ASCII) ─────────────────────────────────────────────────

/// Write the legacy ASCII VTK unstructured grid. Physical groups ride as
/// an integer cell scalar named `physical` (0 = untagged). VTK point ids
/// are the 0-based ranks of the sorted node ids.
pub fn write_vtk_legacy<W: Write>(writer: &mut W, mesh: &Mesh) -> Result<(), MshError> {
    for element in &mesh.elements {
        if vtk_cell(element.etype).is_none() {
            return Err(MshError::UnsupportedElementType(element.etype));
        }
    }
    writeln!(writer, "# vtk DataFile Version 4.2")?;
    writeln!(writer, "rmsh unstructured grid")?;
    writeln!(writer, "ASCII")?;
    writeln!(writer, "DATASET UNSTRUCTURED_GRID")?;

    let mut nodes: Vec<&rmsh_model::Node> = mesh.nodes.values().collect();
    nodes.sort_by_key(|n| n.id);
    writeln!(writer, "POINTS {} double", nodes.len())?;
    for node in &nodes {
        writeln!(
            writer,
            "{} {} {}",
            node.position.x, node.position.y, node.position.z
        )?;
    }

    let total_values: usize = mesh.elements.iter().map(|e| e.node_ids.len() + 1).sum();
    writeln!(writer, "CELLS {} {}", mesh.elements.len(), total_values)?;
    for element in &mesh.elements {
        let (_, perm) = vtk_cell(element.etype).unwrap_or((0, &[]));
        write!(writer, "{}", element.node_ids.len())?;
        for node in permute_nodes(perm, &element.node_ids) {
            let rank = nodes
                .binary_search_by_key(&node, |n| n.id)
                .map_err(|_| MshError::MissingNode(node))?;
            write!(writer, " {rank}")?;
        }
        writeln!(writer)?;
    }

    writeln!(writer, "CELL_TYPES {}", mesh.elements.len())?;
    for element in &mesh.elements {
        let (vtk_type, _) = vtk_cell(element.etype).unwrap_or((0, &[]));
        writeln!(writer, "{vtk_type}")?;
    }

    if mesh.elements.iter().any(|e| e.physical_tag.is_some()) {
        writeln!(writer, "CELL_DATA {}", mesh.elements.len())?;
        writeln!(writer, "SCALARS physical int 1")?;
        writeln!(writer, "LOOKUP_TABLE default")?;
        for element in &mesh.elements {
            writeln!(writer, "{}", element.physical_tag.unwrap_or(0))?;
        }
    }
    Ok(())
}

/// Read a legacy ASCII VTK unstructured grid. Only `ASCII` files with a
/// `DATASET UNSTRUCTURED_GRID` are accepted (the binary legacy flavor is a
/// documented follow-up); unknown cell types are refused by name. The
/// `physical` cell scalar (when present) restores the physical tags.
pub fn read_vtk_legacy(bytes: &[u8]) -> Result<Mesh, MshError> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| MshError::Parse { line: 0, message: "legacy VTK is not UTF-8".into() })?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| MshError::Parse { line: 1, message: "empty file".into() })?;
    if !header.starts_with("# vtk DataFile Version") {
        return Err(MshError::Parse {
            line: 1,
            message: "missing '# vtk DataFile Version' header".into(),
        });
    }
    let _title = lines.next().ok_or_else(|| MshError::Parse {
        line: 2,
        message: "missing title line".into(),
    })?;
    let flavor = lines.next().ok_or_else(|| MshError::Parse {
        line: 3,
        message: "missing ASCII/BINARY line".into(),
    })?;
    if flavor.trim() != "ASCII" {
        return Err(MshError::Parse {
            line: 3,
            message: "only ASCII legacy VTK is supported in v1".into(),
        });
    }
    let dataset = lines
        .find(|l| l.starts_with("DATASET"))
        .ok_or_else(|| MshError::Parse { line: 0, message: "missing DATASET".into() })?;
    if !dataset.starts_with("DATASET UNSTRUCTURED_GRID") {
        return Err(MshError::Parse {
            line: 0,
            message: "only DATASET UNSTRUCTURED_GRID is supported".into(),
        });
    }

    // one word-vector walk over the remaining lines
    let words: Vec<&str> = lines.flat_map(|l| l.split_whitespace()).collect();
    let mut p = 0usize;
    let err = |message: String| MshError::Parse { line: 0, message };
    let mut word = || -> Result<&'static str, MshError> {
        // transmute-free: rebind lifetime through the slice is unnecessary;
        // clone into String instead
        unreachable!()
    };
    let _ = &mut word;
    let _ = err;

    let mut mesh = Mesh::new();
    let mut point_ids: Vec<u64> = Vec::new();
    let mut cell_nodes: Vec<Vec<u64>> = Vec::new();
    let mut physical: Vec<i32> = Vec::new();

    // sequential scan with an index (no lifetimes games)
    let at = |i: usize| -> Option<&str> { words.get(i).copied() };
    while p < words.len() {
        match at(p) {
            Some("POINTS") => {
                p += 1;
                let n: usize =
                    at(p).ok_or_else(|| err("POINTS: missing count".into()))?.parse().map_err(|_| err("POINTS: bad count".into()))?;
                p += 1;
                let _flavor = at(p);
                p += 1;
                point_ids.clear();
                for i in 0..n {
                    let x: f64 = at(p).ok_or_else(|| err("POINTS: missing coordinates".into()))?.parse().map_err(|_| err("POINTS: bad coordinate".into()))?;
                    let y: f64 = at(p + 1).ok_or_else(|| err("POINTS: missing coordinates".into()))?.parse().map_err(|_| err("POINTS: bad coordinate".into()))?;
                    let z: f64 = at(p + 2).ok_or_else(|| err("POINTS: missing coordinates".into()))?.parse().map_err(|_| err("POINTS: bad coordinate".into()))?;
                    p += 3;
                    let id = (i + 1) as u64;
                    mesh.add_node(rmsh_model::Node::new(id, x, y, z));
                    point_ids.push(id);
                }
            }
            Some("CELLS") => {
                p += 1;
                let n: usize = at(p).ok_or_else(|| err("CELLS: missing count".into()))?.parse().map_err(|_| err("CELLS: bad count".into()))?;
                p += 1;
                p += 1; // size field
                cell_nodes.clear();
                for _ in 0..n {
                    let count: usize = at(p).ok_or_else(|| err("CELLS: missing count".into()))?.parse().map_err(|_| err("CELLS: bad count".into()))?;
                    p += 1;
                    let mut ids = Vec::with_capacity(count);
                    for _ in 0..count {
                        let rank: usize = at(p).ok_or_else(|| err("CELLS: missing node".into()))?.parse().map_err(|_| err("CELLS: bad node".into()))?;
                        p += 1;
                        ids.push(*point_ids.get(rank).ok_or_else(|| {
                            err(format!("CELLS: point rank {rank} out of range"))
                        })?);
                    }
                    cell_nodes.push(ids);
                }
            }
            Some("CELL_TYPES") => {
                p += 1;
                let n: usize = at(p).ok_or_else(|| err("CELL_TYPES: missing count".into()))?.parse().map_err(|_| err("CELL_TYPES: bad count".into()))?;
                p += 1;
                let mut next_element_id: u64 = 1;
                if cell_nodes.len() != n {
                    return Err(err(format!(
                        "CELL_TYPES ({n}) and CELLS ({}) counts disagree",
                        cell_nodes.len()
                    )));
                }
                for ids in &cell_nodes {
                    let vtk_type: i32 = at(p).ok_or_else(|| err("CELL_TYPES: missing type".into()))?.parse().map_err(|_| err("CELL_TYPES: bad type".into()))?;
                    p += 1;
                    let (etype, perm) = cell_from_vtk(vtk_type).ok_or_else(|| {
                        err(format!(
                            "unsupported VTK cell type {vtk_type} ({})",
                            vtk_type_name(vtk_type)
                        ))
                    })?;
                    // perm is already the inverse of the write permutation
                    let node_ids = permute_nodes(&perm, ids);
                    let element = Element::new(next_element_id, etype, node_ids);
                    next_element_id += 1;
                    mesh.add_element(element);
                }
            }
            Some("SCALARS") => {
                // "SCALARS <name> <type> <count>" then "LOOKUP_TABLE default"
                // then values — retained only for the physical tag
                let name = at(p + 1).unwrap_or("").to_string();
                p += 3; // SCALARS name type
                p += 1; // count
                if at(p) == Some("LOOKUP_TABLE") {
                    p += 2; // LOOKUP_TABLE default
                }
                if name == "physical" {
                    physical.clear();
                    while p < words.len() {
                        let Ok(v) = at(p).unwrap().parse::<i32>() else { break };
                        p += 1;
                        physical.push(v);
                    }
                }
            }
            Some(_) => p += 1,
            None => break,
        }
    }

    for (element, tag) in mesh.elements.iter_mut().zip(physical.iter()) {
        if *tag != 0 {
            element.physical_tag = Some(*tag);
        }
    }
    Ok(mesh)
}

/// The VTK cell type + inverse permutation for a VTK cell type id (the
/// read direction).
fn cell_from_vtk(vtk_type: i32) -> Option<(ElementType, Vec<usize>)> {
    vtk_cell(vtk_cell_source(vtk_type)?).map(|(_, perm)| {
        let inv = invert(&perm);
        (vtk_cell_source(vtk_type).unwrap(), inv)
    })
}

/// The rmsh element type a VTK cell type maps to (None = unsupported).
fn vtk_cell_source(vtk_type: i32) -> Option<ElementType> {
    match vtk_type {
        1 => Some(ElementType::Point1),
        3 => Some(ElementType::Line2),
        21 => Some(ElementType::Line3),
        35 => Some(ElementType::Line4),
        5 => Some(ElementType::Triangle3),
        22 => Some(ElementType::Triangle6),
        9 => Some(ElementType::Quad4),
        23 => Some(ElementType::Quad8),
        28 => Some(ElementType::Quad9),
        10 => Some(ElementType::Tetrahedron4),
        24 => Some(ElementType::Tetrahedron10),
        12 => Some(ElementType::Hexahedron8),
        25 => Some(ElementType::Hexahedron20),
        29 => Some(ElementType::Hexahedron27),
        13 => Some(ElementType::Prism6),
        26 => Some(ElementType::Prism15),
        14 => Some(ElementType::Pyramid5),
        27 => Some(ElementType::Pyramid13),
        _ => None,
    }
}

// ─── VTU (.vtu, inline ASCII) ─────────────────────────────────────────────────

fn xml_attr(tag: &str, name: &str) -> Option<String> {
    let needle = format!("{name}=\"");
    let start = tag.find(&needle)? + needle.len();
    let rel_end = tag[start..].find('"')?;
    Some(tag[start..start + rel_end].to_string())
}

/// Write the inline-ASCII VTU flavor (ParaView-readable; the appended /
/// base64 binary payload forms are documented follow-ups).
pub fn write_vtu<W: Write>(writer: &mut W, mesh: &Mesh) -> Result<(), MshError> {
    for element in &mesh.elements {
        if vtk_cell(element.etype).is_none() {
            return Err(MshError::UnsupportedElementType(element.etype));
        }
    }
    writeln!(
        writer,
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
    )?;
    writeln!(writer, "  <UnstructuredGrid>")?;
    let mut nodes: Vec<&rmsh_model::Node> = mesh.nodes.values().collect();
    nodes.sort_by_key(|n| n.id);
    writeln!(
        writer,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        nodes.len(),
        mesh.elements.len()
    )?;
    writeln!(writer, "      <Points>")?;
    writeln!(writer, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")?;
    for node in &nodes {
        writeln!(
            writer,
            "          {} {} {}",
            node.position.x, node.position.y, node.position.z
        )?;
    }
    writeln!(writer, "        </DataArray>")?;
    writeln!(writer, "      </Points>")?;
    writeln!(writer, "      <Cells>")?;

    let mut connectivity = String::new();
    let mut offsets = String::new();
    let mut types = String::new();
    let mut offset: usize = 0;
    for element in &mesh.elements {
        let (vtk_type, perm) = vtk_cell(element.etype).unwrap_or((0, &[]));
        let vtk_nodes = permute_nodes(perm, &element.node_ids);
        for node in &vtk_nodes {
            let rank = nodes
                .binary_search_by_key(node, |n| n.id)
                .map_err(|_| MshError::MissingNode(*node))?;
            connectivity.push_str(&format!("{rank} "));
        }
        offset += vtk_nodes.len();
        offsets.push_str(&format!("{offset} "));
        types.push_str(&format!("{vtk_type} "));
    }
    writeln!(
        writer,
        "        <DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">{}</DataArray>",
        connectivity.trim_end()
    )?;
    writeln!(
        writer,
        "        <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">{}</DataArray>",
        offsets.trim_end()
    )?;
    writeln!(
        writer,
        "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">{}</DataArray>",
        types.trim_end()
    )?;
    writeln!(writer, "      </Cells>")?;

    if mesh.elements.iter().any(|e| e.physical_tag.is_some()) {
        writeln!(writer, "      <CellData Scalars=\"physical\">")?;
        writeln!(writer, "        <DataArray type=\"Int32\" Name=\"physical\" format=\"ascii\">")?;
        let tags: Vec<String> = mesh
            .elements
            .iter()
            .map(|e| e.physical_tag.unwrap_or(0).to_string())
            .collect();
        writeln!(writer, "          {}", tags.join(" "))?;
        writeln!(writer, "        </DataArray>")?;
        writeln!(writer, "      </CellData>")?;
    }
    writeln!(writer, "    </Piece>")?;
    writeln!(writer, "  </UnstructuredGrid>")?;
    writeln!(writer, "</VTKFile>")?;
    Ok(())
}

/// The values of the first `<DataArray ... Name="name" ...>` whose body
/// parses as `u64` (empty on absence/mismatch).
fn data_array_u64(text: &str, name: &str) -> Option<Vec<u64>> {
    let marker = format!("Name=\"{name}\"");
    let at = text.find(&marker)?;
    let open_end = text[at..].find('>')? + at + 1;
    let close = text[open_end..].find("</DataArray>")? + open_end;
    Some(
        text[open_end..close]
            .split_whitespace()
            .map(|t| t.parse::<u64>().ok())
            .collect::<Option<Vec<_>>>()?,
    )
}

/// Read an inline-ASCII VTU file (the structure the rmsh writer emits;
/// binary-appended payloads are a documented follow-up). The `physical`
/// cell scalar (when present) restores the physical tags.
pub fn read_vtu(bytes: &[u8]) -> Result<Mesh, MshError> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| MshError::Parse { line: 0, message: "VTU is not UTF-8".into() })?;
    let err = |message: String| MshError::Parse { line: 0, message };
    let piece_tag = text
        .lines()
        .find(|l| l.contains("<Piece"))
        .ok_or_else(|| err("no <Piece> element in the VTU file".into()))?;
    let num_points: usize = xml_attr(piece_tag, "NumberOfPoints")
        .ok_or_else(|| err("<Piece> without NumberOfPoints".into()))?
        .parse()
        .map_err(|_| err("bad NumberOfPoints".into()))?;
    let num_cells: usize = xml_attr(piece_tag, "NumberOfCells")
        .ok_or_else(|| err("<Piece> without NumberOfCells".into()))?
        .parse()
        .map_err(|_| err("bad NumberOfCells".into()))?;

    // points: the first DataArray inside <Points>
    let points_at = text.find("<Points>").ok_or_else(|| err("missing <Points>".into()))?;
    let points_close_rel = text[points_at..]
        .find("</Points>")
        .ok_or_else(|| err("unterminated <Points>".into()))?;
    let points_block = &text[points_at..points_at + points_close_rel];
    let coord_open = points_block
        .find('>')
        .ok_or_else(|| err("malformed Points DataArray".into()))?;
    let coord_close = points_block[coord_open..]
        .find("</DataArray>")
        .ok_or_else(|| err("unterminated Points DataArray".into()))?
        + coord_open;
    let coord_values: Vec<f64> = points_block[coord_open..coord_close]
        .split_whitespace()
        .filter_map(|t| t.parse::<f64>().ok())
        .collect();
    if coord_values.len() < num_points * 3 {
        return Err(err(format!(
            "Points carry {} values for {num_points} points",
            coord_values.len()
        )));
    }

    let mut mesh = Mesh::new();
    for i in 0..num_points {
        mesh.add_node(rmsh_model::Node::new(
            (i + 1) as u64,
            coord_values[i * 3],
            coord_values[i * 3 + 1],
            coord_values[i * 3 + 2],
        ));
    }

    let connectivity = data_array_u64(text, "connectivity")
        .ok_or_else(|| err("missing connectivity DataArray".into()))?;
    let offsets = data_array_u64(text, "offsets")
        .ok_or_else(|| err("missing offsets DataArray".into()))?;
    let types_raw = data_array_u64(text, "types")
        .ok_or_else(|| err("missing types DataArray".into()))?;
    if offsets.len() != num_cells || types_raw.len() != num_cells {
        return Err(err(
            "offsets/types counts disagree with NumberOfCells".into(),
        ));
    }

    // connectivity ranks are 0-based point ranks in file order — the same
    // order the points were added in, so rank r -> node id r + 1
    let mut next_element_id: u64 = 1;
    let mut start = 0usize;
    for (cell, &vtk_type_u) in types_raw.iter().enumerate() {
        let vtk_type = vtk_type_u as i32;
        let end = offsets[cell] as usize;
        let vtk_nodes = connectivity
            .get(start..end)
            .ok_or_else(|| err("connectivity out of bounds".into()))?;
        start = end;
        let etype = vtk_cell_source(vtk_type)
            .ok_or_else(|| err(format!("unsupported VTK cell type {vtk_type} ({})", vtk_type_name(vtk_type))))?;
        let inv = invert(vtk_cell(etype).map(|(_, p)| p).unwrap_or(&[]));
        let permuted = permute_nodes(&inv, vtk_nodes);
        let node_ids: Vec<u64> = permuted.iter().map(|r| r + 1).collect();
        mesh.add_element(Element::new(next_element_id, etype, node_ids));
        next_element_id += 1;
    }

    if let Some(spos) = text.find("Name=\"physical\"") {
        let open_end = text[spos..].find('>').map(|x| x + spos + 1).unwrap_or(spos);
        let close = text[open_end..]
            .find("</DataArray>")
            .map(|x| x + open_end)
            .unwrap_or(text.len());
        let tags: Vec<i32> = text[open_end..close]
            .split_whitespace()
            .map_while(|t| t.parse::<i32>().ok())
            .collect();
        for (element, tag) in mesh.elements.iter_mut().zip(tags.iter()) {
            if *tag != 0 {
                element.physical_tag = Some(*tag);
            }
        }
    }
    Ok(mesh)
}
