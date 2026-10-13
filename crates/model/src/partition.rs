//! Partitioned-mesh semantics (the MSH 4.1 `$PartitionedEntities` and
//! `$GhostElements` sections).
//!
//! The structures retain the partition decomposition verbatim — the
//! entity records carry their parent links, partition lists, bounding
//! boxes, physical tags and the signed boundary-entity references
//! exactly as Gmsh writes them (layout pinned against Gmsh 4.14
//! `GModelIO_MSH4.cpp` and real METIS-partitioned output; see the io
//! tests). The retained records are what the v4 writers re-emit, so a
//! partitioned mesh round trips field-exact.
//!
//! MSH 2.2 has no partition section (Gmsh splits partitioned v2 meshes
//! into per-partition files instead); the v2 writer ignores the
//! partition data and the v2 reader never populates it — documented
//! fidelity line, matching Gmsh's own single-file v2 behavior.

use serde::{Deserialize, Serialize};

/// One partition entity: a piece of a parent entity owned by one
/// partition (or shared by several, for boundary pieces).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PartitionedEntity {
    /// Topological dimension (0 point, 1 curve, 2 surface, 3 volume).
    pub dim: i32,
    /// The partition entity's own tag (per-dimension space).
    pub tag: i32,
    /// Dimension of the parent (unpartitioned) entity.
    pub parent_dim: i32,
    /// Tag of the parent entity (0 when the record carries no parent).
    pub parent_tag: i32,
    /// The partitions sharing this piece.
    pub partitions: Vec<i32>,
    /// The entity's bounding box as written `(min, max)`; `None` when the
    /// record carried the degenerate zero box (empty entity).
    pub bbox: Option<([f64; 3], [f64; 3])>,
    /// Physical tags attached to this partition entity.
    pub physical_tags: Vec<i32>,
    /// The entity's boundary pieces as `(tag, orientation)` pairs
    /// (orientation ±1; the ASCII flavor folds the sign into the tag and
    /// the binary flavor carries the pair — the retained form unifies
    /// both).
    pub bounding: Vec<(i32, i32)>,
}

/// A ghost entity: the placeholder on partition `partition` that holds
/// elements owned by other partitions (written with `-part_ghosts`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GhostEntity {
    pub dim: i32,
    pub tag: i32,
    /// The partition this ghost entity belongs to.
    pub partition: i32,
}

/// One `$GhostElements` cell: an element owned by `partition` that also
/// appears (as a ghost) on the partitions in `ghost_partitions`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GhostElement {
    pub element_tag: u64,
    pub partition: i32,
    pub ghost_partitions: Vec<i32>,
}

/// The retained `$PartitionedEntities` + `$GhostElements` content of a
/// partitioned mesh.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct PartitionInfo {
    /// Number of partitions the mesh was split into.
    pub num_partitions: u64,
    /// All partition entities, in file order (points, then curves,
    /// surfaces, volumes).
    pub entities: Vec<PartitionedEntity>,
    /// The ghost entities from the partitioned-entities header.
    pub ghosts: Vec<GhostEntity>,
    /// The `$GhostElements` cells (empty when the file was written
    /// without ghost cells).
    pub ghost_elements: Vec<GhostElement>,
}
