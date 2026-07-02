use serde::{Deserialize, Serialize};

/// Supported finite element types.
///
/// Both first-order (P1) and second-order (P2) types are represented as named
/// variants.  Gmsh type IDs not listed here fall through to [`Unknown`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ElementType {
    // ── 0-D: Point ──────────────────────────────────────────────────────────
    /// 1-node point
    Point1,

    // ── 1-D: Line / Edge ────────────────────────────────────────────────────
    /// 2-node line (P1)
    Line2,
    /// 3-node line with edge midpoint (P2)
    Line3,

    // ── 2-D: Surface elements ───────────────────────────────────────────────
    /// 3-node triangle (P1)
    Triangle3,
    /// 6-node triangle with edge midpoints (P2)
    Triangle6,
    /// 4-node quadrilateral (P1)
    Quad4,
    /// 9-node quadrilateral with edge midpoints + centre (P2)
    Quad9,

    // ── 3-D: Volume elements ────────────────────────────────────────────────
    /// 4-node tetrahedron (P1)
    Tetrahedron4,
    /// 10-node tetrahedron with edge midpoints (P2)
    Tetrahedron10,
    /// 8-node hexahedron (P1)
    Hexahedron8,
    /// 27-node hexahedron with edge midpoints + face centres + interior (P2)
    Hexahedron27,
    /// 6-node prism / wedge (P1)
    Prism6,
    /// 18-node prism with edge midpoints + face centres (P2)
    Prism18,
    /// 5-node pyramid (P1)
    Pyramid5,
    /// 14-node pyramid with edge midpoints + face centre (P2)
    Pyramid14,

    /// Unknown / unsupported type (holds raw Gmsh type ID).
    Unknown(i32),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ElementFamily {
    Point,
    Line,
    Triangle,
    Quadrilateral,
    Tetrahedron,
    Hexahedron,
    Prism,
    Pyramid,
}

fn family_from_gmsh_type_id(id: i32) -> Option<ElementFamily> {
    match id {
        // Point
        15 => Some(ElementFamily::Point),
        // Line (1st+high order)
        1 | 8 | 26 | 27 | 28 => Some(ElementFamily::Line),
        // Triangle (1st+high order)
        2 | 9 | 20 | 21 | 22 | 23 | 24 | 25 => Some(ElementFamily::Triangle),
        // Quadrilateral (1st+high order)
        3 | 10 | 16 | 36 | 37 | 38 | 47 | 48 | 49 | 50 | 51 => Some(ElementFamily::Quadrilateral),
        // Tetrahedron (1st+high order)
        4 | 11 | 29 | 30 | 31 => Some(ElementFamily::Tetrahedron),
        // Hexahedron (1st+high order)
        5 | 12 | 17 | 92 | 93 => Some(ElementFamily::Hexahedron),
        // Prism (1st+high order)
        6 | 13 | 18 | 90 | 91 => Some(ElementFamily::Prism),
        // Pyramid (1st+high order)
        7 | 14 | 19 | 118 | 119 => Some(ElementFamily::Pyramid),
        _ => None,
    }
}

fn gmsh_dimension_from_type_id(id: i32) -> u8 {
    match family_from_gmsh_type_id(id) {
        Some(ElementFamily::Point) => 0,
        Some(ElementFamily::Line) => 1,
        Some(ElementFamily::Triangle) | Some(ElementFamily::Quadrilateral) => 2,
        Some(ElementFamily::Tetrahedron)
        | Some(ElementFamily::Hexahedron)
        | Some(ElementFamily::Prism)
        | Some(ElementFamily::Pyramid) => 3,
        None => 0,
    }
}

const LINE2_EDGES: &[[usize; 2]] = &[[0, 1]];
const TRI3_EDGES: &[[usize; 2]] = &[[0, 1], [1, 2], [2, 0]];
const QUAD4_EDGES: &[[usize; 2]] = &[[0, 1], [1, 2], [2, 3], [3, 0]];
const TET4_EDGES: &[[usize; 2]] = &[[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]];
const HEX8_EDGES: &[[usize; 2]] = &[
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 4],
    [0, 4],
    [1, 5],
    [2, 6],
    [3, 7],
];
const PRISM6_EDGES: &[[usize; 2]] = &[
    [0, 1],
    [1, 2],
    [2, 0],
    [3, 4],
    [4, 5],
    [5, 3],
    [0, 3],
    [1, 4],
    [2, 5],
];
const PYRAMID5_EDGES: &[[usize; 2]] = &[
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0],
    [0, 4],
    [1, 4],
    [2, 4],
    [3, 4],
];

const TET4_FACES: &[&[usize]] = &[&[0, 1, 2], &[0, 1, 3], &[1, 2, 3], &[0, 2, 3]];
const HEX8_FACES: &[&[usize]] = &[
    &[0, 1, 2, 3],
    &[4, 5, 6, 7],
    &[0, 1, 5, 4],
    &[2, 3, 7, 6],
    &[0, 3, 7, 4],
    &[1, 2, 6, 5],
];
const PRISM6_FACES: &[&[usize]] = &[
    &[0, 1, 2],
    &[3, 4, 5],
    &[0, 1, 4, 3],
    &[1, 2, 5, 4],
    &[0, 2, 5, 3],
];
const PYRAMID5_FACES: &[&[usize]] = &[
    &[0, 1, 2, 3],
    &[0, 1, 4],
    &[1, 2, 4],
    &[2, 3, 4],
    &[0, 3, 4],
];

// ─── Reference-element node positions (for P1→P2 promotion) ────────────────

/// Line3: corner node 0, corner node 1, midpoint
const REF_LINE3: &[(f64, f64, f64)] = &[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 0.0, 0.0)];

/// Triangle6: corners (0,1,2) then edge midpoints (3,4,5) opposite vertices (2,0,1).
const REF_TRI6: &[(f64, f64, f64)] = &[
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
    (0.5, 0.0, 0.0), (0.5, 0.5, 0.0), (0.0, 0.5, 0.0),
];

/// Quad9: corners (0-3) then edge midpoints (4-7) then centre (8).
const REF_QUAD9: &[(f64, f64, f64)] = &[
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.5, 0.0, 0.0), (1.0, 0.5, 0.0), (0.5, 1.0, 0.0), (0.0, 0.5, 0.0),
    (0.5, 0.5, 0.0),
];

/// Tetrahedron10: corners (0-3) then edge-midpoints (4-9) in order matching
/// TET4_EDGES (0-1, 1-2, 2-0, 0-3, 1-3, 2-3).
const REF_TET10: &[(f64, f64, f64)] = &[
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
    (0.5, 0.0, 0.0), (0.5, 0.5, 0.0), (0.0, 0.5, 0.0),
    (0.0, 0.0, 0.5), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5),
];

/// Hexahedron27: corners (0-7), edge midpoints (8-19) in HEX8_EDGES order,
/// face centres (20-25), interior (26).
const REF_HEX27: &[(f64, f64, f64)] = &[
    // corners 0-7
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    // edge midpoints 8-19 (HEX8_EDGES order)
    (0.5, 0.0, 0.0), (1.0, 0.5, 0.0), (0.5, 1.0, 0.0), (0.0, 0.5, 0.0),
    (0.5, 0.0, 1.0), (1.0, 0.5, 1.0), (0.5, 1.0, 1.0), (0.0, 0.5, 1.0),
    (0.0, 0.0, 0.5), (1.0, 0.0, 0.5), (1.0, 1.0, 0.5), (0.0, 1.0, 0.5),
    // face centres 20-25 (HEX8_FACES order)
    (0.5, 0.5, 0.0), (0.5, 0.5, 1.0),
    (0.5, 0.0, 0.5), (0.5, 1.0, 0.5),
    (0.0, 0.5, 0.5), (1.0, 0.5, 0.5),
    // interior 26
    (0.5, 0.5, 0.5),
];

impl ElementType {
    /// Convert from Gmsh element type ID (MSH v4 format).
    pub fn from_gmsh_type_id(id: i32) -> Self {
        match id {
            15 => ElementType::Point1,
            1 => ElementType::Line2,
            8 => ElementType::Line3,
            2 => ElementType::Triangle3,
            9 => ElementType::Triangle6,
            3 => ElementType::Quad4,
            10 => ElementType::Quad9,
            4 => ElementType::Tetrahedron4,
            11 => ElementType::Tetrahedron10,
            5 => ElementType::Hexahedron8,
            12 => ElementType::Hexahedron27,
            6 => ElementType::Prism6,
            13 => ElementType::Prism18,
            7 => ElementType::Pyramid5,
            14 => ElementType::Pyramid14,
            _ => ElementType::Unknown(id),
        }
    }

    /// Convert to Gmsh element type ID (MSH v4 format).
    pub fn to_gmsh_type_id(&self) -> i32 {
        match self {
            ElementType::Point1 => 15,
            ElementType::Line2 => 1,
            ElementType::Line3 => 8,
            ElementType::Triangle3 => 2,
            ElementType::Triangle6 => 9,
            ElementType::Quad4 => 3,
            ElementType::Quad9 => 10,
            ElementType::Tetrahedron4 => 4,
            ElementType::Tetrahedron10 => 11,
            ElementType::Hexahedron8 => 5,
            ElementType::Hexahedron27 => 12,
            ElementType::Prism6 => 6,
            ElementType::Prism18 => 13,
            ElementType::Pyramid5 => 7,
            ElementType::Pyramid14 => 14,
            ElementType::Unknown(id) => *id,
        }
    }

    /// Number of nodes for this element type.
    pub fn node_count(&self) -> usize {
        match self {
            ElementType::Point1 => 1,
            ElementType::Line2 => 2,
            ElementType::Line3 => 3,
            ElementType::Triangle3 => 3,
            ElementType::Triangle6 => 6,
            ElementType::Quad4 => 4,
            ElementType::Quad9 => 9,
            ElementType::Tetrahedron4 => 4,
            ElementType::Tetrahedron10 => 10,
            ElementType::Hexahedron8 => 8,
            ElementType::Hexahedron27 => 27,
            ElementType::Prism6 => 6,
            ElementType::Prism18 => 18,
            ElementType::Pyramid5 => 5,
            ElementType::Pyramid14 => 14,
            ElementType::Unknown(_) => 0,
        }
    }

    /// Topological dimension of this element (0=point, 1=edge, 2=face, 3=volume).
    pub fn dimension(&self) -> u8 {
        match self {
            ElementType::Point1 => 0,
            ElementType::Line2 | ElementType::Line3 => 1,
            ElementType::Triangle3 | ElementType::Triangle6 | ElementType::Quad4 | ElementType::Quad9 => 2,
            ElementType::Tetrahedron4
            | ElementType::Tetrahedron10
            | ElementType::Hexahedron8
            | ElementType::Hexahedron27
            | ElementType::Prism6
            | ElementType::Prism18
            | ElementType::Pyramid5
            | ElementType::Pyramid14 => 3,
            ElementType::Unknown(id) => gmsh_dimension_from_type_id(*id),
        }
    }

    /// Return the faces of a volume element as arrays of local node indices.
    /// Each face is a slice of node indices (3 for triangular faces, 4 for quad faces).
    /// Same topology for P1 and P2 variants — returns corner-node patterns.
    pub fn faces(&self) -> &[&[usize]] {
        match self {
            ElementType::Tetrahedron4 | ElementType::Tetrahedron10 => TET4_FACES,
            ElementType::Hexahedron8 | ElementType::Hexahedron27 => HEX8_FACES,
            ElementType::Prism6 | ElementType::Prism18 => PRISM6_FACES,
            ElementType::Pyramid5 | ElementType::Pyramid14 => PYRAMID5_FACES,
            ElementType::Unknown(id) => match family_from_gmsh_type_id(*id) {
                Some(ElementFamily::Tetrahedron) => TET4_FACES,
                Some(ElementFamily::Hexahedron) => HEX8_FACES,
                Some(ElementFamily::Prism) => PRISM6_FACES,
                Some(ElementFamily::Pyramid) => PYRAMID5_FACES,
                _ => &[],
            },
            _ => &[],
        }
    }

    /// Return the edges of an element as pairs of local node indices.
    /// Same topology for P1 and P2 variants — returns corner-node pairs.
    pub fn edges(&self) -> &[[usize; 2]] {
        match self {
            ElementType::Line2 | ElementType::Line3 => LINE2_EDGES,
            ElementType::Triangle3 | ElementType::Triangle6 => TRI3_EDGES,
            ElementType::Quad4 | ElementType::Quad9 => QUAD4_EDGES,
            ElementType::Tetrahedron4 | ElementType::Tetrahedron10 => TET4_EDGES,
            ElementType::Hexahedron8 | ElementType::Hexahedron27 => HEX8_EDGES,
            ElementType::Prism6 | ElementType::Prism18 => PRISM6_EDGES,
            ElementType::Pyramid5 | ElementType::Pyramid14 => PYRAMID5_EDGES,
            ElementType::Unknown(id) => match family_from_gmsh_type_id(*id) {
                Some(ElementFamily::Line) => LINE2_EDGES,
                Some(ElementFamily::Triangle) => TRI3_EDGES,
                Some(ElementFamily::Quadrilateral) => QUAD4_EDGES,
                Some(ElementFamily::Tetrahedron) => TET4_EDGES,
                Some(ElementFamily::Hexahedron) => HEX8_EDGES,
                Some(ElementFamily::Prism) => PRISM6_EDGES,
                Some(ElementFamily::Pyramid) => PYRAMID5_EDGES,
                _ => &[],
            },
            _ => &[],
        }
    }

    /// Reference-element node positions suitable for computing edge midpoints.
    ///
    /// Returns `Some(slice)` for known types with the parametric position of
    /// each node on the reference element.  Used when promoting P1→P2 to
    /// determine which nodes are edge midpoints vs corner nodes.
    pub fn reference_node_positions(&self) -> Option<&'static [(f64, f64, f64)]> {
        match self {
            ElementType::Line3 => Some(&REF_LINE3),
            ElementType::Triangle6 => Some(&REF_TRI6),
            ElementType::Quad9 => Some(&REF_QUAD9),
            ElementType::Tetrahedron10 => Some(&REF_TET10),
            ElementType::Hexahedron27 => Some(&REF_HEX27),
            _ => None,
        }
    }
}

/// A finite element consisting of a type and connectivity (node IDs).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Element {
    pub id: u64,
    pub etype: ElementType,
    /// Physical group tag
    pub physical_tag: Option<i32>,
    /// Node IDs forming this element (global IDs referencing `Node::id`).
    pub node_ids: Vec<u64>,
}

impl Element {
    pub fn new(id: u64, etype: ElementType, node_ids: Vec<u64>) -> Self {
        Self {
            id,
            etype,
            physical_tag: None,
            node_ids,
        }
    }

    pub fn dimension(&self) -> u8 {
        self.etype.dimension()
    }
}

#[cfg(test)]
mod tests {
    use super::ElementType;

    #[test]
    fn unknown_high_order_types_have_correct_dimension() {
        assert_eq!(ElementType::Unknown(8).dimension(), 1); // line3
        assert_eq!(ElementType::Unknown(9).dimension(), 2); // tri6
        assert_eq!(ElementType::Unknown(10).dimension(), 2); // quad9
        assert_eq!(ElementType::Unknown(11).dimension(), 3); // tet10
        assert_eq!(ElementType::Unknown(12).dimension(), 3); // hex27
        assert_eq!(ElementType::Unknown(13).dimension(), 3); // prism18
        assert_eq!(ElementType::Unknown(14).dimension(), 3); // pyramid14
    }

    #[test]
    fn unknown_volume_families_expose_canonical_faces_and_edges() {
        assert_eq!(ElementType::Unknown(11).faces().len(), 4);
        assert_eq!(ElementType::Unknown(11).edges().len(), 6);

        assert_eq!(ElementType::Unknown(12).faces().len(), 6);
        assert_eq!(ElementType::Unknown(12).edges().len(), 12);

        assert_eq!(ElementType::Unknown(13).faces().len(), 5);
        assert_eq!(ElementType::Unknown(13).edges().len(), 9);

        assert_eq!(ElementType::Unknown(14).faces().len(), 5);
        assert_eq!(ElementType::Unknown(14).edges().len(), 8);
    }

    #[test]
    fn all_element_types_have_correct_dimension() {
        // 0-Dimensional (Point)
        assert_eq!(ElementType::Point1.dimension(), 0);

        // 1-Dimensional (Curve/Edge)
        assert_eq!(ElementType::Line2.dimension(), 1);
        assert_eq!(ElementType::Line3.dimension(), 1);

        // 2-Dimensional (Surface/Face)
        assert_eq!(ElementType::Triangle3.dimension(), 2);
        assert_eq!(ElementType::Triangle6.dimension(), 2);
        assert_eq!(ElementType::Quad4.dimension(), 2);
        assert_eq!(ElementType::Quad9.dimension(), 2);

        // 3-Dimensional (Volume/Region)
        assert_eq!(ElementType::Tetrahedron4.dimension(), 3);
        assert_eq!(ElementType::Tetrahedron10.dimension(), 3);
        assert_eq!(ElementType::Hexahedron8.dimension(), 3);
        assert_eq!(ElementType::Hexahedron27.dimension(), 3);
        assert_eq!(ElementType::Prism6.dimension(), 3);
        assert_eq!(ElementType::Prism18.dimension(), 3);
        assert_eq!(ElementType::Pyramid5.dimension(), 3);
        assert_eq!(ElementType::Pyramid14.dimension(), 3);

        // Unknown types should infer dimension from Gmsh type ID
        // Gmsh type IDs: 15=point, 1–8,26–28=line, 2–3,9–10,16,20–25,36–51=face, 4–7,11–14,17–19,29–31,90–93,118–119=volume
        assert_eq!(ElementType::Unknown(15).dimension(), 0); // point
        assert_eq!(ElementType::Unknown(1).dimension(), 1); // line2
        assert_eq!(ElementType::Unknown(2).dimension(), 2); // tri3
        assert_eq!(ElementType::Unknown(3).dimension(), 2); // quad4
        assert_eq!(ElementType::Unknown(4).dimension(), 3); // tet4
        assert_eq!(ElementType::Unknown(5).dimension(), 3); // hex8
    }

    #[test]
    fn element_dimension_aligns_with_node_count() {
        // Point: 1 node
        assert_eq!(ElementType::Point1.node_count(), 1);
        assert_eq!(ElementType::Point1.dimension(), 0);

        // Line: 2 nodes (P1) and 3 nodes (P2)
        assert_eq!(ElementType::Line2.node_count(), 2);
        assert_eq!(ElementType::Line2.dimension(), 1);
        assert_eq!(ElementType::Line3.node_count(), 3);
        assert_eq!(ElementType::Line3.dimension(), 1);

        // Triangle: 3 (P1) and 6 (P2)
        assert_eq!(ElementType::Triangle3.node_count(), 3);
        assert_eq!(ElementType::Triangle3.dimension(), 2);
        assert_eq!(ElementType::Triangle6.node_count(), 6);
        assert_eq!(ElementType::Triangle6.dimension(), 2);

        // Quad: 4 (P1) and 9 (P2)
        assert_eq!(ElementType::Quad4.node_count(), 4);
        assert_eq!(ElementType::Quad4.dimension(), 2);
        assert_eq!(ElementType::Quad9.node_count(), 9);
        assert_eq!(ElementType::Quad9.dimension(), 2);

        // Tet: 4 (P1) and 10 (P2)
        assert_eq!(ElementType::Tetrahedron4.node_count(), 4);
        assert_eq!(ElementType::Tetrahedron4.dimension(), 3);
        assert_eq!(ElementType::Tetrahedron10.node_count(), 10);
        assert_eq!(ElementType::Tetrahedron10.dimension(), 3);

        // Hex: 8 (P1) and 27 (P2)
        assert_eq!(ElementType::Hexahedron8.node_count(), 8);
        assert_eq!(ElementType::Hexahedron8.dimension(), 3);
        assert_eq!(ElementType::Hexahedron27.node_count(), 27);
        assert_eq!(ElementType::Hexahedron27.dimension(), 3);

        // Prism: 6 (P1) and 18 (P2)
        assert_eq!(ElementType::Prism6.node_count(), 6);
        assert_eq!(ElementType::Prism6.dimension(), 3);
        assert_eq!(ElementType::Prism18.node_count(), 18);
        assert_eq!(ElementType::Prism18.dimension(), 3);

        // Pyramid: 5 (P1) and 14 (P2)
        assert_eq!(ElementType::Pyramid5.node_count(), 5);
        assert_eq!(ElementType::Pyramid5.dimension(), 3);
        assert_eq!(ElementType::Pyramid14.node_count(), 14);
        assert_eq!(ElementType::Pyramid14.dimension(), 3);
    }
}
