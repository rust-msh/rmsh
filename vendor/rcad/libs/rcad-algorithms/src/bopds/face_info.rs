use std::collections::HashSet;

/// Per-face intersection bookkeeping (OCCT: BOPDS_FaceInfo).
#[derive(Debug, Clone, Default)]
pub struct FaceInfo {
    /// Indices of PaveBlocks that lie ON this face (from E-F intersection).
    pub pave_blocks_on: HashSet<usize>,
    /// Indices of intersection curves that lie IN this face (from F-F intersection).
    pub curves_in: HashSet<usize>,
    /// Vertex indices that lie ON this face.
    pub vertices_on: HashSet<usize>,
    /// Vertex indices that lie IN this face (from F-F intersection).
    pub vertices_in: HashSet<usize>,
}
