//! `$Periodic` section read/write for all four MSH flavors.
//!
//! Pinned layouts (verified against real Gmsh 4.14 output; the fixtures were
//! generated with `Periodic Surface {2} = {1} Translate {1,0,0}` and read
//! back byte-for-byte):
//!
//! v2.2 — ASCII form, used for BOTH v2 flavors because Gmsh keeps this
//! section in ASCII even inside binary files (only `$Nodes`/`$Elements`
//! payloads go binary):
//!
//! ```text
//! $Periodic
//! <numLinks>
//! <entityDim> <entityTagSlave> <entityTagMaster>   (one line, 3 ints)
//! Affine <v1> ... <vn>                             (optional line; literal
//!                                                   "Affine" keyword, values
//!                                                   verbatim on one line)
//! <numCorrespondingNodes>                          (one line)
//! <nodeSlave> <nodeMaster>                         (one pair per line)
//! ...                                              (per link)
//! $EndPeriodic
//! ```
//!
//! The slave tag is written FIRST (Gmsh: `dimension entity-tag
//! master-entity-tag`), and each correspondence pair is
//! `(slave node, master node)`. The affine count is NOT stored in v2; the
//! optional `Affine` keyword is what makes the line recognizable.
//!
//! v4.1 ASCII — same shape, but the affine count leads the values on one
//! line (a bare `0` when the link has no transform):
//!
//! ```text
//! $Periodic
//! <numLinks>
//! <entityDim> <entityTagSlave> <entityTagMaster>
//! <numAffine> <v1> ... <vn>
//! <numCorrespondingNodes>
//! <nodeSlave> <nodeMaster>
//! ...
//! $EndPeriodic
//! ```
//!
//! v4.1 binary — the `$Periodic` header stays ASCII and the payloads are
//! little-endian binary: `u64 numLinks`; per link `i32 entityDim, i32
//! entityTagSlave, i32 entityTagMaster, u64 numAffine, f64 values, u64
//! numCorrespondingNodes, u64 node pairs`; then a newline and
//! `$EndPeriodic` on their own line.

use std::io::{BufRead, Write};

use rmsh_model::PeriodicLink;

use crate::msh::{bin_read_f64, bin_read_i32, bin_read_u64, le_f64, le_i32, le_u64, MshError};

/// Line-oriented source for the ASCII periodic parsers. The v2 binary
/// reader feeds it from the byte cursor; the ASCII readers feed it from
/// their line iterators (preserving line numbers for located errors).
pub(crate) trait LineSource {
    fn next_line(&mut self) -> Result<String, MshError>;

    /// Best-effort location for error reports (0 where no line number
    /// exists, i.e. the byte-cursor binary source).
    fn line_hint(&self) -> usize;

    fn parse_error(&self, message: String) -> MshError {
        MshError::Parse {
            line: self.line_hint(),
            message,
        }
    }
}

/// Adapter over the ASCII readers' `Lines` iterator + line counter.
pub(crate) struct BufLineSource<'a, R: BufRead> {
    pub(crate) lines: &'a mut std::io::Lines<R>,
    pub(crate) line_num: &'a mut usize,
}

impl<R: BufRead> LineSource for BufLineSource<'_, R> {
    fn next_line(&mut self) -> Result<String, MshError> {
        *self.line_num += 1;
        self.lines
            .next()
            .ok_or_else(|| MshError::Parse {
                line: *self.line_num,
                message: "Unexpected end of file".into(),
            })?
            .map_err(MshError::Io)
    }

    fn line_hint(&self) -> usize {
        *self.line_num
    }
}

/// Adapter over the v2 binary reader's byte cursor.
pub(crate) struct ByteLineSource<'a> {
    pub(crate) bytes: &'a [u8],
    pub(crate) pos: &'a mut usize,
}

impl LineSource for ByteLineSource<'_> {
    fn next_line(&mut self) -> Result<String, MshError> {
        let start = *self.pos;
        let slice = &self.bytes[start..];
        let Some(nl) = slice.iter().position(|&b| b == b'\n') else {
            return Err(MshError::Parse {
                line: 0,
                message: "binary v2: unexpected end of file expecting a line".into(),
            });
        };
        let line = std::str::from_utf8(&slice[..nl]).map_err(|_| MshError::Parse {
            line: 0,
            message: "binary v2: section line is not UTF-8".into(),
        })?;
        *self.pos = start + nl + 1;
        Ok(line.trim().to_string())
    }

    fn line_hint(&self) -> usize {
        0
    }
}

/// Parse `<count>` `(slave, master)` node pairs from one-or-more lines.
/// Gmsh writes one pair per line but the layout is token-based, so pairs
/// packed several per line parse identically.
fn read_pairs(src: &mut impl LineSource, num_pairs: usize) -> Result<Vec<(u64, u64)>, MshError> {
    let mut tokens: Vec<u64> = Vec::with_capacity(num_pairs * 2);
    while tokens.len() < num_pairs * 2 {
        let line = src.next_line()?;
        for token in line.split_whitespace() {
            let value: u64 = token.parse().map_err(|_| {
                src.parse_error(format!("Invalid node tag in periodic correspondence: {token}"))
            })?;
            tokens.push(value);
        }
    }
    if tokens.len() != num_pairs * 2 {
        return Err(src.parse_error(format!(
            "Periodic correspondence declares {num_pairs} pair(s) but carries {} node tag(s)",
            tokens.len()
        )));
    }
    Ok(tokens
        .chunks(2)
        .map(|pair| (pair[0], pair[1]))
        .collect())
}

/// Parse the three-int entity header `<entityDim> <slave> <master>`.
fn parse_entity_line(
    src: &mut impl LineSource,
    line: &str,
) -> Result<(u8, i32, i32), MshError> {
    let parts: Vec<&str> = line.trim().split_whitespace().collect();
    if parts.len() != 3 {
        return Err(src.parse_error(format!(
            "Invalid periodic entity line (expected 'entityDim entityTag entityTagMaster'): {line}"
        )));
    }
    let entity_dim: u8 = parts[0].parse().map_err(|_| {
        src.parse_error(format!("Invalid periodic entity dimension: {}", parts[0]))
    })?;
    let slave: i32 = parts[1]
        .parse()
        .map_err(|_| src.parse_error(format!("Invalid periodic entity tag: {}", parts[1])))?;
    let master: i32 = parts[2].parse().map_err(|_| {
        src.parse_error(format!("Invalid periodic master entity tag: {}", parts[2]))
    })?;
    Ok((entity_dim, slave, master))
}

/// Parse the v2.2 `$Periodic` payload (after the `$Periodic` marker line,
/// up to and excluding `$EndPeriodic`).
pub(crate) fn parse_periodic_v2(src: &mut impl LineSource) -> Result<Vec<PeriodicLink>, MshError> {
    let count_line = src.next_line()?;
    let num_links: usize = count_line.trim().parse().map_err(|_| {
        src.parse_error(format!("Invalid periodic link count: {count_line}"))
    })?;

    let mut links = Vec::with_capacity(num_links);
    for _ in 0..num_links {
        let header = src.next_line()?;
        let (entity_dim, slave, master) = parse_entity_line(src, &header)?;

        // The affine line is optional and introduces itself with the
        // literal `Affine` keyword; anything else is already the
        // correspondence count.
        let mut affine = Vec::new();
        let mut next = src.next_line()?;
        if next.trim_start().starts_with("Affine") {
            let values = &next.trim_start()["Affine".len()..];
            for token in values.split_whitespace() {
                let value: f64 = token.parse().map_err(|_| {
                    src.parse_error(format!("Invalid periodic affine value: {token}"))
                })?;
                affine.push(value);
            }
            next = src.next_line()?;
        }

        let num_pairs: usize = next.trim().parse().map_err(|_| {
            src.parse_error(format!(
                "Invalid periodic correspondence count: {next}"
            ))
        })?;
        let correspondence = read_pairs(src, num_pairs)?;

        links.push(PeriodicLink {
            entity_dim,
            master,
            slave,
            affine,
            correspondence,
        });
    }
    Ok(links)
}

/// Parse the v4.1 ASCII `$Periodic` payload (after the `$Periodic` marker
/// line, up to and excluding `$EndPeriodic`).
pub(crate) fn parse_periodic_v4(src: &mut impl LineSource) -> Result<Vec<PeriodicLink>, MshError> {
    let count_line = src.next_line()?;
    let num_links: usize = count_line.trim().parse().map_err(|_| {
        src.parse_error(format!("Invalid periodic link count: {count_line}"))
    })?;

    let mut links = Vec::with_capacity(num_links);
    for _ in 0..num_links {
        let header = src.next_line()?;
        let (entity_dim, slave, master) = parse_entity_line(src, &header)?;

        // The affine count leads the values on one line; tolerate the
        // values spilling onto continuation lines.
        let mut line = src.next_line()?;
        let mut tokens = line.split_whitespace();
        let num_affine: usize = tokens
            .next()
            .unwrap_or_default()
            .parse()
            .map_err(|_| src.parse_error(format!("Invalid periodic affine count: {line}")))?;
        let mut affine: Vec<f64> = Vec::with_capacity(num_affine);
        for token in tokens.by_ref() {
            let value: f64 = token.parse().map_err(|_| {
                src.parse_error(format!("Invalid periodic affine value: {token}"))
            })?;
            affine.push(value);
        }
        while affine.len() < num_affine {
            line = src.next_line()?;
            for token in line.split_whitespace() {
                let value: f64 = token.parse().map_err(|_| {
                    src.parse_error(format!("Invalid periodic affine value: {token}"))
                })?;
                affine.push(value);
            }
        }
        if affine.len() != num_affine {
            return Err(src.parse_error(format!(
                "Periodic link declares {num_affine} affine value(s) but carries {}",
                affine.len()
            )));
        }

        let count_line = src.next_line()?;
        let num_pairs: usize = count_line.trim().parse().map_err(|_| {
            src.parse_error(format!(
                "Invalid periodic correspondence count: {count_line}"
            ))
        })?;
        let correspondence = read_pairs(src, num_pairs)?;

        links.push(PeriodicLink {
            entity_dim,
            master,
            slave,
            affine,
            correspondence,
        });
    }
    Ok(links)
}

/// Parse the v4.1 binary `$Periodic` payload (after the ASCII `$Periodic`
/// marker line, up to the newline + `$EndPeriodic` line, which the caller
/// consumes).
pub(crate) fn parse_periodic_v4_binary(
    bytes: &[u8],
    pos: &mut usize,
) -> Result<Vec<PeriodicLink>, MshError> {
    let num_links = bin_read_u64(bytes, pos, "periodic link count")? as usize;
    let mut links = Vec::with_capacity(num_links);
    for _ in 0..num_links {
        let entity_dim = bin_read_i32(bytes, pos, "periodic entity dimension")?;
        let entity_dim = u8::try_from(entity_dim).map_err(|_| MshError::Parse {
            line: 0,
            message: format!("binary v4: invalid periodic entity dimension {entity_dim}"),
        })?;
        let slave = bin_read_i32(bytes, pos, "periodic entity tag")?;
        let master = bin_read_i32(bytes, pos, "periodic master entity tag")?;
        let num_affine = bin_read_u64(bytes, pos, "periodic affine count")? as usize;
        let mut affine = Vec::with_capacity(num_affine);
        for _ in 0..num_affine {
            affine.push(bin_read_f64(bytes, pos, "periodic affine value")?);
        }
        let num_pairs = bin_read_u64(bytes, pos, "periodic node count")? as usize;
        let mut correspondence = Vec::with_capacity(num_pairs);
        for _ in 0..num_pairs {
            let slave_node = bin_read_u64(bytes, pos, "periodic slave node")?;
            let master_node = bin_read_u64(bytes, pos, "periodic master node")?;
            correspondence.push((slave_node, master_node));
        }
        links.push(PeriodicLink {
            entity_dim,
            master,
            slave,
            affine,
            correspondence,
        });
    }
    Ok(links)
}

/// Structurally skip the v4.1 binary `$Entities` payload (after the ASCII
/// `$Entities` marker line, up to the newline + `$EndEntities` line, which
/// the caller consumes). Real Gmsh always writes this section in v4.1, so
/// a reader that cannot consume its binary payload cannot read real files.
/// The layout (all little endian):
/// `u64 numPoints; u64 numCurves; u64 numSurfaces; u64 numVolumes;` then
/// per point: `i32 tag; f64 x; f64 y; f64 z; u64 numPhysicalTags;
/// i32 physicalTag x n;` and per curve/surface/volume: `i32 tag;
/// f64 minX; f64 minY; f64 minZ; f64 maxX; f64 maxY; f64 maxZ;
/// u64 numPhysicalTags; i32 physicalTag x n; u64 numBounding;
/// i32 boundingTag x n;`. The content is discarded: rmsh keys entity
/// information off the element blocks, not off `$Entities`.
pub(crate) fn skip_entities_v4_binary(bytes: &[u8], pos: &mut usize) -> Result<(), MshError> {
    let num_points = bin_read_u64(bytes, pos, "entity point count")? as usize;
    let num_curves = bin_read_u64(bytes, pos, "entity curve count")? as usize;
    let num_surfaces = bin_read_u64(bytes, pos, "entity surface count")? as usize;
    let num_volumes = bin_read_u64(bytes, pos, "entity volume count")? as usize;

    // Points carry one coordinate triple; the higher-dimensional entities
    // share one record layout (a 6-double bounding box).
    for _ in 0..num_points {
        let _tag = bin_read_i32(bytes, pos, "entity tag")?;
        let _x = bin_read_f64(bytes, pos, "entity coordinate")?;
        let _y = bin_read_f64(bytes, pos, "entity coordinate")?;
        let _z = bin_read_f64(bytes, pos, "entity coordinate")?;
        skip_physical_tags(bytes, pos)?;
    }
    for _ in 0..(num_curves + num_surfaces + num_volumes) {
        let _tag = bin_read_i32(bytes, pos, "entity tag")?;
        let _min_x = bin_read_f64(bytes, pos, "entity bounding box")?;
        let _min_y = bin_read_f64(bytes, pos, "entity bounding box")?;
        let _min_z = bin_read_f64(bytes, pos, "entity bounding box")?;
        let _max_x = bin_read_f64(bytes, pos, "entity bounding box")?;
        let _max_y = bin_read_f64(bytes, pos, "entity bounding box")?;
        let _max_z = bin_read_f64(bytes, pos, "entity bounding box")?;
        skip_physical_tags(bytes, pos)?;
        let num_bounding = bin_read_u64(bytes, pos, "entity bounding count")? as usize;
        for _ in 0..num_bounding {
            let _tag = bin_read_i32(bytes, pos, "entity bounding tag")?;
        }
    }
    Ok(())
}

/// Consume a `u64 numPhysicalTags; i32 physicalTag x n;` run.
fn skip_physical_tags(bytes: &[u8], pos: &mut usize) -> Result<(), MshError> {
    let num_physical = bin_read_u64(bytes, pos, "entity physical tag count")? as usize;
    for _ in 0..num_physical {
        let _tag = bin_read_i32(bytes, pos, "entity physical tag")?;
    }
    Ok(())
}

/// Write the v2.2 `$Periodic` section (ASCII form). A no-op when `links`
/// is empty, so periodic-free files stay byte-identical to before. Also
/// used by the v2 BINARY writer: Gmsh keeps this section ASCII there.
pub(crate) fn write_periodic_v2<W: Write>(
    writer: &mut W,
    links: &[PeriodicLink],
) -> Result<(), MshError> {
    if links.is_empty() {
        return Ok(());
    }
    writeln!(writer, "$Periodic")?;
    writeln!(writer, "{}", links.len())?;
    for link in links {
        writeln!(writer, "{} {} {}", link.entity_dim, link.slave, link.master)?;
        if !link.affine.is_empty() {
            write!(writer, "Affine")?;
            for value in &link.affine {
                write!(writer, " {value}")?;
            }
            writeln!(writer)?;
        }
        writeln!(writer, "{}", link.correspondence.len())?;
        for (slave_node, master_node) in &link.correspondence {
            writeln!(writer, "{slave_node} {master_node}")?;
        }
    }
    writeln!(writer, "$EndPeriodic")?;
    Ok(())
}

/// Write the v4.1 ASCII `$Periodic` section. A no-op when `links` is empty.
pub(crate) fn write_periodic_v4<W: Write>(
    writer: &mut W,
    links: &[PeriodicLink],
) -> Result<(), MshError> {
    if links.is_empty() {
        return Ok(());
    }
    writeln!(writer, "$Periodic")?;
    writeln!(writer, "{}", links.len())?;
    for link in links {
        writeln!(writer, "{} {} {}", link.entity_dim, link.slave, link.master)?;
        write!(writer, "{}", link.affine.len())?;
        for value in &link.affine {
            write!(writer, " {value}")?;
        }
        writeln!(writer)?;
        writeln!(writer, "{}", link.correspondence.len())?;
        for (slave_node, master_node) in &link.correspondence {
            writeln!(writer, "{slave_node} {master_node}")?;
        }
    }
    writeln!(writer, "$EndPeriodic")?;
    Ok(())
}

/// Write the v4.1 binary `$Periodic` section: ASCII header line, binary
/// payloads, then the newline + `$EndPeriodic` line. A no-op when `links`
/// is empty.
pub(crate) fn write_periodic_v4_binary<W: Write>(
    writer: &mut W,
    links: &[PeriodicLink],
) -> Result<(), MshError> {
    if links.is_empty() {
        return Ok(());
    }
    writeln!(writer, "$Periodic")?;
    le_u64(writer, links.len() as u64)?;
    for link in links {
        le_i32(writer, i32::from(link.entity_dim))?;
        le_i32(writer, link.slave)?;
        le_i32(writer, link.master)?;
        le_u64(writer, link.affine.len() as u64)?;
        for value in &link.affine {
            le_f64(writer, *value)?;
        }
        le_u64(writer, link.correspondence.len() as u64)?;
        for (slave_node, master_node) in &link.correspondence {
            le_u64(writer, *slave_node)?;
            le_u64(writer, *master_node)?;
        }
    }
    writeln!(writer)?;
    writeln!(writer, "$EndPeriodic")?;
    Ok(())
}
