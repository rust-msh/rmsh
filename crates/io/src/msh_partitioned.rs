//! `$PartitionedEntities` + `$GhostElements` parse and write (MSH 4.1).
//!
//! Layout authority: Gmsh 4.14 `GModelIO_MSH4.cpp` (the partitioned
//! read/write arms route through `readMSH4Entities(partition=true)` /
//! `writeMSH4Entities(partition=true)` and `readMSH4GhostElements`), pinned
//! against real METIS-partitioned output (a 3-partition cube,
//! `testdata/partitioned_cube_v4a.msh` ASCII and
//! `testdata/partitioned_cube_v4b.msh` binary). The retained
//! [`PartitionInfo`](rmsh_model::PartitionInfo) records round trip
//! field-exact.
//!
//! Record grammar (ASCII shown; binary mirrors with LE u64/i32/f64):
//!
//! ```text
//! $PartitionedEntities
//! numPartitions
//! ghostSize
//! ghostSize × (entityTag partition)          // ghost entities
//! numVertices numEdges numFaces numRegions
//! // per point:  4 + numParts + 3(bbox) + 1 + numPhys
//! // per dim>0:  4 + numParts + 6(bbox) + 1 + numPhys + 1 + numBound
//! tag parentDim parentTag numParts parts... bbox... numPhys phys...
//! [dim 1: numBoundVertices signedVertexTags...]
//! [dim 2: numBoundEdges    signedEdgeTags...]
//! [dim 3: numBoundSurfaces signedSurfaceTags...]
//! $EndPartitionedEntities
//!
//! $GhostElements
//! numGhostCells
//! numGhostCells × (elementTag partition numGhostParts ghostParts...)
//! $EndGhostElements
//! ```
//!
use std::io::Write;

use rmsh_model::{GhostElement, GhostEntity, PartitionInfo, PartitionedEntity};

use crate::MshError;

// ─── ASCII ────────────────────────────────────────────────────────────────────

/// Parse the `$PartitionedEntities` section from the raw lines between the
/// markers (token stream — records are whitespace-aligned, not
/// line-aligned).
pub(crate) fn parse_partitioned_v4_ascii(
    lines: &[&str],
    section_line: usize,
) -> Result<PartitionInfo, MshError> {
    let err = |message: String| MshError::Parse {
        line: section_line,
        message,
    };
    let mut tokens = lines.iter().flat_map(|l| l.split_whitespace());
    let mut next = || {
        tokens
            .next()
            .ok_or_else(|| err("unexpected end of $PartitionedEntities".into()))
    };

    let num_partitions: u64 = next()?.parse().map_err(|_| err("invalid numPartitions".into()))?;
    let ghost_size: u64 = next()?.parse().map_err(|_| err("invalid ghostSize".into()))?;
    let mut ghosts = Vec::new();
    for _ in 0..ghost_size {
        let tag: i32 = next()?.parse().map_err(|_| err("invalid ghost tag".into()))?;
        let partition: i32 =
            next()?.parse().map_err(|_| err("invalid ghost partition".into()))?;
        // ghost entities in the header are volumes (the partition's volume
        // placeholder); the dim is not stored in the section — retained as
        // 3, the only form Gmsh 4.14 writes after -part_ghosts.
        ghosts.push(GhostEntity {
            dim: 3,
            tag,
            partition,
        });
    }
    let mut counts = [0u64; 4];
    for c in &mut counts {
        *c = next()?.parse().map_err(|_| err("invalid entity count".into()))?;
    }

    let mut entities = Vec::new();
    for (dim, n) in counts.iter().enumerate() {
        let dim = dim as i32;
        for _ in 0..*n {
            let tag: i32 = next()?.parse().map_err(|_| err("invalid entity tag".into()))?;
            let parent_dim: i32 = next()?.parse().map_err(|_| err("invalid parent dim".into()))?;
            let parent_tag: i32 = next()?.parse().map_err(|_| err("invalid parent tag".into()))?;
            let num_parts: u64 = next()?.parse().map_err(|_| err("invalid partition count".into()))?;
            let mut partitions = Vec::new();
            for _ in 0..num_parts {
                partitions.push(next()?.parse().map_err(|_| err("invalid partition id".into()))?);
            }
            let nb = if dim == 0 { 3 } else { 6 };
            let mut bbox_vals = [0.0f64; 6];
            for v in bbox_vals.iter_mut().take(nb) {
                *v = next()?.parse().map_err(|_| err("invalid bounding box".into()))?;
            }
            let bbox = (
                [bbox_vals[0], bbox_vals[1], bbox_vals[2]],
                [bbox_vals[3], bbox_vals[4], bbox_vals[5]],
            );
            let num_phys: u64 = next()?.parse().map_err(|_| err("invalid physical count".into()))?;
            let mut physical_tags = Vec::new();
            for _ in 0..num_phys {
                physical_tags.push(next()?.parse().map_err(|_| err("invalid physical tag".into()))?);
            }
            let mut bounding = Vec::new();
            if dim > 0 {
                let num_bound: u64 =
                    next()?.parse().map_err(|_| err("invalid bounding count".into()))?;
                for _ in 0..num_bound {
                    let signed: i32 =
                        next()?.parse().map_err(|_| err("invalid bounding tag".into()))?;
                    let (bound_tag, sign) = if signed < 0 { (-signed, -1) } else { (signed, 1) };
                    bounding.push((bound_tag, sign));
                }
            }
            entities.push(PartitionedEntity {
                dim,
                tag,
                parent_dim,
                parent_tag,
                partitions,
                bbox: Some(bbox),
                physical_tags,
                bounding,
            });
        }
    }
    Ok(PartitionInfo {
        num_partitions,
        entities,
        ghosts,
        ghost_elements: Vec::new(),
    })
}

/// Parse the `$GhostElements` section (the 4.1 flat-cell form; the 4.2
/// block form is not produced by the writers rmsh interops with).
pub(crate) fn parse_ghost_elements_v4_ascii(
    lines: &[&str],
    section_line: usize,
) -> Result<Vec<GhostElement>, MshError> {
    let err = |message: String| MshError::Parse {
        line: section_line,
        message,
    };
    let mut tokens = lines.iter().flat_map(|l| l.split_whitespace());
    let mut next = || {
        tokens
            .next()
            .ok_or_else(|| err("unexpected end of $GhostElements".into()))
    };
    let num_cells: u64 = next()?.parse().map_err(|_| err("invalid ghost cell count".into()))?;
    let mut cells = Vec::new();
    for _ in 0..num_cells {
        let element_tag: u64 =
            next()?.parse().map_err(|_| err("invalid ghost element tag".into()))?;
        let partition: i32 = next()?.parse().map_err(|_| err("invalid ghost partition".into()))?;
        let num_ghosts: u64 =
            next()?.parse().map_err(|_| err("invalid ghost partition count".into()))?;
        let mut ghost_partitions = Vec::new();
        for _ in 0..num_ghosts {
            ghost_partitions
                .push(next()?.parse().map_err(|_| err("invalid ghost partition".into()))?);
        }
        cells.push(GhostElement {
            element_tag,
            partition,
            ghost_partitions,
        });
    }
    Ok(cells)
}

pub(crate) fn write_partitioned_v4(
    writer: &mut impl Write,
    info: &PartitionInfo,
) -> Result<(), MshError> {
    writeln!(writer, "$PartitionedEntities")?;
    writeln!(writer, "{}", info.num_partitions)?;
    writeln!(writer, "{}", info.ghosts.len())?;
    for ghost in &info.ghosts {
        writeln!(writer, "{} {}", ghost.tag, ghost.partition)?;
    }
    let count = |d: i32| info.entities.iter().filter(|e| e.dim == d).count();
    writeln!(writer, "{} {} {} {}", count(0), count(1), count(2), count(3))?;
    for e in &info.entities {
        let mut line = String::new();
        line.push_str(&format!("{} {} {} ", e.tag, e.parent_dim, e.parent_tag));
        line.push_str(&format!("{} ", e.partitions.len()));
        for p in &e.partitions {
            line.push_str(&format!("{p} "));
        }
        let (min, max) = e.bbox.unwrap_or(([0.0; 3], [0.0; 3]));
        if e.dim == 0 {
            line.push_str(&format!("{} {} {} ", min[0], min[1], min[2]));
        } else {
            line.push_str(&format!(
                "{} {} {} {} {} {} ",
                min[0], min[1], min[2], max[0], max[1], max[2]
            ));
        }
        line.push_str(&format!("{} ", e.physical_tags.len()));
        for t in &e.physical_tags {
            line.push_str(&format!("{t} "));
        }
        if e.dim > 0 {
            line.push_str(&format!("{} ", e.bounding.len()));
            for (bound_tag, sign) in &e.bounding {
                line.push_str(&format!("{} ", bound_tag * sign));
            }
        }
        writeln!(writer, "{}", line.trim_end())?;
    }
    writeln!(writer, "$EndPartitionedEntities")?;
    Ok(())
}

pub(crate) fn write_ghost_elements_v4(
    writer: &mut impl Write,
    cells: &[GhostElement],
) -> Result<(), MshError> {
    if cells.is_empty() {
        return Ok(());
    }
    writeln!(writer, "$GhostElements")?;
    writeln!(writer, "{}", cells.len())?;
    for cell in cells {
        let mut line = format!("{} {} {}", cell.element_tag, cell.partition, cell.ghost_partitions.len());
        for gp in &cell.ghost_partitions {
            line.push_str(&format!(" {gp}"));
        }
        writeln!(writer, "{line}")?;
    }
    writeln!(writer, "$EndGhostElements")?;
    Ok(())
}

// ─── Binary ───────────────────────────────────────────────────────────────────

/// Parse `$PartitionedEntities` from the binary bytes at `pos` (the marker
/// line has been consumed). Mirrors `readMSH4Entities(partition=true)`:
/// LE u64 counts, i32 tags/dims, f64 boxes (3 values for points, 6 for
/// dim > 0), bounding as u64 count + signed i32 values.
pub(crate) fn parse_partitioned_v4_binary(
    bytes: &[u8],
    pos: &mut usize,
) -> Result<PartitionInfo, MshError> {
    let err = |message: String| MshError::Parse {
        line: 0,
        message,
    };
    let u64_at = |pos: &mut usize| -> Result<u64, MshError> {
        if *pos + 8 > bytes.len() {
            return Err(err("truncated $PartitionedEntities".into()));
        }
        let v = u64::from_le_bytes(bytes[*pos..*pos + 8].try_into().unwrap());
        *pos += 8;
        Ok(v)
    };
    let i32_at = |pos: &mut usize| -> Result<i32, MshError> {
        if *pos + 4 > bytes.len() {
            return Err(err("truncated $PartitionedEntities".into()));
        }
        let v = i32::from_le_bytes(bytes[*pos..*pos + 4].try_into().unwrap());
        *pos += 4;
        Ok(v)
    };
    let f64_at = |pos: &mut usize| -> Result<f64, MshError> {
        if *pos + 8 > bytes.len() {
            return Err(err("truncated $PartitionedEntities".into()));
        }
        let v = f64::from_le_bytes(bytes[*pos..*pos + 8].try_into().unwrap());
        *pos += 8;
        Ok(v)
    };

    let num_partitions = u64_at(pos)?;
    let ghost_size = u64_at(pos)?;
    let mut ghosts = Vec::new();
    for _ in 0..ghost_size {
        let tag = i32_at(pos)?;
        let partition = i32_at(pos)?;
        ghosts.push(GhostEntity {
            dim: 3,
            tag,
            partition,
        });
    }
    let mut counts = [0u64; 4];
    for c in &mut counts {
        *c = u64_at(pos)?;
    }
    let mut entities = Vec::new();
    for (dim, n) in counts.iter().enumerate() {
        let dim = dim as i32;
        for _ in 0..*n {
            let tag = i32_at(pos)?;
            let parent_dim = i32_at(pos)?;
            let parent_tag = i32_at(pos)?;
            let num_parts = u64_at(pos)?;
            let mut partitions = Vec::new();
            for _ in 0..num_parts {
                partitions.push(i32_at(pos)?);
            }
            let nb = if dim == 0 { 3 } else { 6 };
            let mut bbox_vals = [0.0f64; 6];
            for v in bbox_vals.iter_mut().take(nb) {
                *v = f64_at(pos)?;
            }
            let bbox = (
                [bbox_vals[0], bbox_vals[1], bbox_vals[2]],
                [bbox_vals[3], bbox_vals[4], bbox_vals[5]],
            );
            let num_phys = u64_at(pos)?;
            let mut physical_tags = Vec::new();
            for _ in 0..num_phys {
                physical_tags.push(i32_at(pos)?);
            }
            let mut bounding = Vec::new();
            if dim > 0 {
                // Empirically pinned against the real binary fixture: the
                // 4.14 binary writer emits numBound SIGNED values (the
                // orientation folded into the sign), the same single-value
                // form as ASCII.
                let num_bound = u64_at(pos)?;
                for _ in 0..num_bound {
                    let signed = i32_at(pos)?;
                    let (bound_tag, sign) = if signed < 0 {
                        (-signed, -1)
                    } else {
                        (signed, 1)
                    };
                    bounding.push((bound_tag, sign));
                }
            }
            entities.push(PartitionedEntity {
                dim,
                tag,
                parent_dim,
                parent_tag,
                partitions,
                bbox: Some(bbox),
                physical_tags,
                bounding,
            });
        }
    }
    Ok(PartitionInfo {
        num_partitions,
        entities,
        ghosts,
        ghost_elements: Vec::new(),
    })
}

pub(crate) fn write_partitioned_v4_binary(
    writer: &mut impl Write,
    info: &PartitionInfo,
) -> Result<(), MshError> {
    let put_u64 = |w: &mut dyn Write, v: u64| w.write_all(&v.to_le_bytes()).map_err(MshError::Io);
    let put_i32 = |w: &mut dyn Write, v: i32| w.write_all(&v.to_le_bytes()).map_err(MshError::Io);
    let put_f64 = |w: &mut dyn Write, v: f64| w.write_all(&v.to_le_bytes()).map_err(MshError::Io);

    let mut w = writer;
    writeln!(w, "$PartitionedEntities")?;
    put_u64(&mut w, info.num_partitions)?;
    put_u64(&mut w, info.ghosts.len() as u64)?;
    for ghost in &info.ghosts {
        put_i32(&mut w, ghost.tag)?;
        put_i32(&mut w, ghost.partition)?;
    }
    let count = |d: i32| info.entities.iter().filter(|e| e.dim == d).count();
    put_u64(&mut w, count(0) as u64)?;
    put_u64(&mut w, count(1) as u64)?;
    put_u64(&mut w, count(2) as u64)?;
    put_u64(&mut w, count(3) as u64)?;
    for e in &info.entities {
        put_i32(&mut w, e.tag)?;
        put_i32(&mut w, e.parent_dim)?;
        put_i32(&mut w, e.parent_tag)?;
        put_u64(&mut w, e.partitions.len() as u64)?;
        for p in &e.partitions {
            put_i32(&mut w, *p)?;
        }
        let nb = if e.dim == 0 { 3 } else { 6 };
        let (min, max) = e.bbox.unwrap_or(([0.0; 3], [0.0; 3]));
        for v in [min[0], min[1], min[2], max[0], max[1], max[2]].into_iter().take(nb) {
            put_f64(&mut w, v)?;
        }
        put_u64(&mut w, e.physical_tags.len() as u64)?;
        for t in &e.physical_tags {
            put_i32(&mut w, *t)?;
        }
        if e.dim > 0 {
            put_u64(&mut w, e.bounding.len() as u64)?;
            for (bound_tag, sign) in &e.bounding {
                put_i32(&mut w, bound_tag * sign)?;
            }
        }
    }
    writeln!(w, "$EndPartitionedEntities")?;
    Ok(())
}

/// Parse `$GhostElements` (4.1 flat-cell form).
pub(crate) fn parse_ghost_elements_v4_binary(
    bytes: &[u8],
    pos: &mut usize,
) -> Result<Vec<GhostElement>, MshError> {
    let err = |message: String| MshError::Parse {
        line: 0,
        message,
    };
    let u64_at = |pos: &mut usize| -> Result<u64, MshError> {
        if *pos + 8 > bytes.len() {
            return Err(err("truncated $GhostElements".into()));
        }
        let v = u64::from_le_bytes(bytes[*pos..*pos + 8].try_into().unwrap());
        *pos += 8;
        Ok(v)
    };
    let i32_at = |pos: &mut usize| -> Result<i32, MshError> {
        if *pos + 4 > bytes.len() {
            return Err(err("truncated $GhostElements".into()));
        }
        let v = i32::from_le_bytes(bytes[*pos..*pos + 4].try_into().unwrap());
        *pos += 4;
        Ok(v)
    };
    let num_cells = u64_at(pos)?;
    let mut cells = Vec::new();
    for _ in 0..num_cells {
        let element_tag = u64_at(pos)?;
        let partition = i32_at(pos)?;
        let num_ghosts = u64_at(pos)?;
        let mut ghost_partitions = Vec::new();
        for _ in 0..num_ghosts {
            ghost_partitions.push(i32_at(pos)?);
        }
        cells.push(GhostElement {
            element_tag,
            partition,
            ghost_partitions,
        });
    }
    Ok(cells)
}

pub(crate) fn write_ghost_elements_v4_binary(
    writer: &mut impl Write,
    cells: &[GhostElement],
) -> Result<(), MshError> {
    if cells.is_empty() {
        return Ok(());
    }
    let put_u64 = |w: &mut dyn Write, v: u64| w.write_all(&v.to_le_bytes()).map_err(MshError::Io);
    let put_i32 = |w: &mut dyn Write, v: i32| w.write_all(&v.to_le_bytes()).map_err(MshError::Io);
    let mut w = writer;
    writeln!(w, "$GhostElements")?;
    put_u64(&mut w, cells.len() as u64)?;
    for cell in cells {
        put_u64(&mut w, cell.element_tag)?;
        put_i32(&mut w, cell.partition)?;
        put_u64(&mut w, cell.ghost_partitions.len() as u64)?;
        for gp in &cell.ghost_partitions {
            put_i32(&mut w, *gp)?;
        }
    }
    writeln!(w, "$EndGhostElements")?;
    Ok(())
}
