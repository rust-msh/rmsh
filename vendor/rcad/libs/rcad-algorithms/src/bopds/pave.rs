/// A parametric point along an edge's curve (OCCT: BOPDS_Pave).
#[derive(Debug, Clone, Copy)]
pub struct Pave {
    /// Index of the vertex at this parametric point (in DS.vertices).
    pub vertex_idx: usize,
    /// Parametric value on the edge's curve.
    pub param: f64,
}

/// A segment of an edge between two paves (OCCT: BOPDS_PaveBlock).
/// When an edge is split by intersections, it becomes multiple PaveBlocks.
#[derive(Debug, Clone)]
pub struct PaveBlock {
    /// Original edge index in DS.edges.
    pub original_edge: usize,
    pub pave1: Pave,
    pub pave2: Pave,
    /// New edge index assigned during result building.
    pub new_edge: Option<usize>,
}

impl PaveBlock {
    pub fn new(original_edge: usize, pave1: Pave, pave2: Pave) -> Self {
        Self {
            original_edge,
            pave1,
            pave2,
            new_edge: None,
        }
    }
}
