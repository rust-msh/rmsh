use serde::{Deserialize, Serialize};

/// Supported finite element types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ElementType {
    /// 2-node line
    Line2,
    /// 3-node triangle
    Triangle3,
    /// 4-node quadrilateral
    Quad4,
    /// 4-node tetrahedron
    Tetrahedron4,
    /// 8-node hexahedron
    Hexahedron8,
    /// 6-node prism (wedge)
    Prism6,
    /// 5-node pyramid
    Pyramid5,
    /// 1-node point
    Point1,
    // ── Semantic high-order variants (Gmsh MSH type ids in parentheses;
    // orders 2-3, the charter's high-order business line). Read back from
    // MSH by [`ElementType::from_gmsh_type_id`], written by
    // [`ElementType::to_gmsh_type_id`]. ──
    /// 3-node quadratic line (type 8)
    Line3,
    /// 4-node cubic line (type 26)
    Line4,
    /// 6-node quadratic triangle (type 9)
    Triangle6,
    /// 9-node cubic triangle, Gmsh incomplete variant (type 20)
    Triangle9,
    /// 10-node cubic triangle, complete (type 21)
    Triangle10,
    /// 8-node quadratic quadrilateral, Gmsh incomplete/serendipity variant (type 16)
    Quad8,
    /// 9-node quadratic quadrilateral, complete (type 10)
    Quad9,
    /// 12-node cubic quadrilateral, Gmsh incomplete/serendipity variant (type 39)
    Quad12,
    /// 16-node cubic quadrilateral, complete (type 36)
    Quad16,
    /// 10-node quadratic tetrahedron (type 11)
    Tetrahedron10,
    /// 20-node cubic tetrahedron (type 29)
    Tetrahedron20,
    /// 20-node quadratic hexahedron, Gmsh incomplete/serendipity variant (type 17)
    Hexahedron20,
    /// 27-node quadratic hexahedron, complete (type 12)
    Hexahedron27,
    /// 64-node cubic hexahedron, complete (type 92)
    Hexahedron64,
    /// 15-node quadratic prism, Gmsh incomplete/serendipity variant (type 18)
    Prism15,
    /// 18-node quadratic prism, complete (type 13)
    Prism18,
    /// 40-node cubic prism, complete (type 90)
    Prism40,
    /// 13-node quadratic pyramid, Gmsh incomplete variant (type 19)
    Pyramid13,
    /// 14-node quadratic pyramid, complete (type 14)
    Pyramid14,
    /// 30-node cubic pyramid (type 118)
    Pyramid30,
    /// Unknown / unsupported type
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
        // Line (orders 1-10; 62-66 are the order 6-10 edges)
        1 | 8 | 26 | 27 | 28 | 62 | 63 | 64 | 65 | 66 => Some(ElementFamily::Line),
        // Triangle (complete + Gmsh incomplete variants, orders 1-10)
        2 | 9 | 20..=25 | 32 | 33 | 42..=46 | 52..=56 => Some(ElementFamily::Triangle),
        // Quadrilateral (complete + Gmsh incomplete variants, orders 1-10)
        3 | 10 | 16 | 36..=41 | 47..=51 | 57..=61 => Some(ElementFamily::Quadrilateral),
        // Tetrahedron (orders 1-5)
        4 | 11 | 29 | 30 | 31 => Some(ElementFamily::Tetrahedron),
        // Hexahedron (8/27 complete, 20 incomplete, 64/125 for orders 3-4)
        5 | 12 | 17 | 92 | 93 => Some(ElementFamily::Hexahedron),
        // Prism (6/18 complete, 15 incomplete, 40/75 for orders 3-4)
        6 | 13 | 18 | 90 | 91 => Some(ElementFamily::Prism),
        // Pyramid (5/14 complete, 13 incomplete, 30/55 for orders 3-4)
        7 | 14 | 19 | 118 | 119 => Some(ElementFamily::Pyramid),
        _ => None,
    }
}

/// Gmsh MSH element type number -> node count, for every type whose count
/// could be verified. Two evidence sources:
/// - types 1-31, 92 and 93 appear explicitly in the official Gmsh manual's
///   MSH file format table;
/// - every other entry was read back empirically from meshes written by the
///   real Gmsh 4.14 (high-order probes of every family; counts extracted
///   from the v4.1 block headers, `dim tag type count` plus one line per
///   element).
/// Types with a variable node count (polygons, polyhedra) are intentionally
/// absent: a wrong count would corrupt every read that indexes node lists
/// by it.
pub fn gmsh_node_count(id: i32) -> Option<usize> {
    let count = match id {
        15 => 1,
        // 2- to 11-node lines (orders 1-10).
        1 => 2,
        8 => 3,
        26 => 4,
        27 => 5,
        28 => 6,
        62 => 7,
        63 => 8,
        64 => 9,
        65 => 10,
        66 => 11,
        // Triangles: complete and Gmsh incomplete variants, orders 1-10.
        2 => 3,
        9 => 6,
        20 => 9,
        21 => 10,
        22 => 12,
        23 => 15,
        24 => 15,
        25 => 21,
        32 => 22,
        33 => 28,
        42 => 28,
        43 => 36,
        44 => 45,
        45 => 55,
        46 => 66,
        52 => 18,
        53 => 21,
        54 => 24,
        55 => 27,
        56 => 30,
        // Quadrangles: complete (orders 1-10) and Gmsh incomplete variants.
        3 => 4,
        10 => 9,
        16 => 8,
        36 => 16,
        37 => 25,
        38 => 36,
        39 => 12,
        40 => 16,
        41 => 20,
        47 => 49,
        48 => 64,
        49 => 81,
        50 => 100,
        51 => 121,
        57 => 24,
        58 => 28,
        59 => 32,
        60 => 36,
        61 => 40,
        // Tetrahedra, orders 1-5.
        4 => 4,
        11 => 10,
        29 => 20,
        30 => 35,
        31 => 56,
        // Hexahedra: 8/27 complete, 20 incomplete, 64/125 for orders 3-4.
        5 => 8,
        12 => 27,
        17 => 20,
        92 => 64,
        93 => 125,
        // Prisms: 6/18 complete, 15 incomplete, 40/75 for orders 3-4.
        6 => 6,
        13 => 18,
        18 => 15,
        90 => 40,
        91 => 75,
        // Pyramids: 5/14 complete, 13 incomplete, 30/55 for orders 3-4.
        7 => 5,
        14 => 14,
        19 => 13,
        118 => 30,
        119 => 55,
        _ => return None,
    };
    Some(count)
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

impl ElementType {
    /// Convert from Gmsh element type ID (MSH v4 format). The order 2-3
    /// types load as their semantic variants; every other listed id keeps
    /// the [`ElementType::Unknown`] fidelity form.
    pub fn from_gmsh_type_id(id: i32) -> Self {
        match id {
            15 => ElementType::Point1,
            1 => ElementType::Line2,
            2 => ElementType::Triangle3,
            3 => ElementType::Quad4,
            4 => ElementType::Tetrahedron4,
            5 => ElementType::Hexahedron8,
            6 => ElementType::Prism6,
            7 => ElementType::Pyramid5,
            // Semantic high-order set (order 2-3).
            8 => ElementType::Line3,
            26 => ElementType::Line4,
            9 => ElementType::Triangle6,
            20 => ElementType::Triangle9,
            21 => ElementType::Triangle10,
            16 => ElementType::Quad8,
            10 => ElementType::Quad9,
            39 => ElementType::Quad12,
            36 => ElementType::Quad16,
            11 => ElementType::Tetrahedron10,
            29 => ElementType::Tetrahedron20,
            17 => ElementType::Hexahedron20,
            12 => ElementType::Hexahedron27,
            92 => ElementType::Hexahedron64,
            18 => ElementType::Prism15,
            13 => ElementType::Prism18,
            90 => ElementType::Prism40,
            19 => ElementType::Pyramid13,
            14 => ElementType::Pyramid14,
            118 => ElementType::Pyramid30,
            _ => ElementType::Unknown(id),
        }
    }

    /// Gmsh MSH type id of this element type; `None` only for an
    /// [`ElementType::Unknown`] whose id is absent from the verified
    /// node-count table (such an id has no valid Gmsh element mapping).
    pub fn to_gmsh_type_id(&self) -> Option<i32> {
        match self {
            ElementType::Point1 => Some(15),
            ElementType::Line2 => Some(1),
            ElementType::Triangle3 => Some(2),
            ElementType::Quad4 => Some(3),
            ElementType::Tetrahedron4 => Some(4),
            ElementType::Hexahedron8 => Some(5),
            ElementType::Prism6 => Some(6),
            ElementType::Pyramid5 => Some(7),
            ElementType::Line3 => Some(8),
            ElementType::Line4 => Some(26),
            ElementType::Triangle6 => Some(9),
            ElementType::Triangle9 => Some(20),
            ElementType::Triangle10 => Some(21),
            ElementType::Quad8 => Some(16),
            ElementType::Quad9 => Some(10),
            ElementType::Quad12 => Some(39),
            ElementType::Quad16 => Some(36),
            ElementType::Tetrahedron10 => Some(11),
            ElementType::Tetrahedron20 => Some(29),
            ElementType::Hexahedron20 => Some(17),
            ElementType::Hexahedron27 => Some(12),
            ElementType::Hexahedron64 => Some(92),
            ElementType::Prism15 => Some(18),
            ElementType::Prism18 => Some(13),
            ElementType::Prism40 => Some(90),
            ElementType::Pyramid13 => Some(19),
            ElementType::Pyramid14 => Some(14),
            ElementType::Pyramid30 => Some(118),
            ElementType::Unknown(id) => {
                if gmsh_node_count(*id).is_some() {
                    Some(*id)
                } else {
                    None
                }
            }
        }
    }

    /// Number of nodes for this element type. Every type with a verified
    /// Gmsh table entry (semantic variants included) resolves through the
    /// table; genuinely unknown ids report 0 (callers treat 0 as "not
    /// validated").
    pub fn node_count(&self) -> usize {
        match self {
            ElementType::Unknown(id) => gmsh_node_count(*id).unwrap_or(0),
            other => other
                .to_gmsh_type_id()
                .and_then(gmsh_node_count)
                .unwrap_or(0),
        }
    }

    /// Topological dimension of this element (0=point, 1=edge, 2=face, 3=volume).
    pub fn dimension(&self) -> u8 {
        let id = match self {
            ElementType::Unknown(id) => *id,
            other => match other.to_gmsh_type_id() {
                Some(id) => id,
                None => return 0,
            },
        };
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

    /// Return the faces of a volume element as arrays of local node indices.
    /// Each face is a slice of node indices (3 for triangular faces, 4 for quad faces).
    pub fn faces(&self) -> &[&[usize]] {
        let id = match self {
            ElementType::Unknown(id) => *id,
            other => match other.to_gmsh_type_id() {
                Some(id) => id,
                None => return &[],
            },
        };
        match family_from_gmsh_type_id(id) {
            Some(ElementFamily::Tetrahedron) => TET4_FACES,
            Some(ElementFamily::Hexahedron) => HEX8_FACES,
            Some(ElementFamily::Prism) => PRISM6_FACES,
            Some(ElementFamily::Pyramid) => PYRAMID5_FACES,
            _ => &[],
        }
    }

    /// Return the edges of an element as pairs of local node indices.
    pub fn edges(&self) -> &[[usize; 2]] {
        let id = match self {
            ElementType::Unknown(id) => *id,
            other => match other.to_gmsh_type_id() {
                Some(id) => id,
                None => return &[],
            },
        };
        match family_from_gmsh_type_id(id) {
            Some(ElementFamily::Line) => LINE2_EDGES,
            Some(ElementFamily::Triangle) => TRI3_EDGES,
            Some(ElementFamily::Quadrilateral) => QUAD4_EDGES,
            Some(ElementFamily::Tetrahedron) => TET4_EDGES,
            Some(ElementFamily::Hexahedron) => HEX8_EDGES,
            Some(ElementFamily::Prism) => PRISM6_EDGES,
            Some(ElementFamily::Pyramid) => PYRAMID5_EDGES,
            _ => &[],
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
    use super::{Element, ElementType};

    #[test]
    fn gmsh_node_count_table_matches_verified_sources() {
        // Doc-verified (official manual MSH table).
        assert_eq!(super::gmsh_node_count(15), Some(1));
        assert_eq!(super::gmsh_node_count(1), Some(2));
        assert_eq!(super::gmsh_node_count(11), Some(10)); // Tet10
        assert_eq!(super::gmsh_node_count(12), Some(27)); // Hex27
        assert_eq!(super::gmsh_node_count(92), Some(64)); // Hex64, order 3
        assert_eq!(super::gmsh_node_count(93), Some(125)); // Hex125, order 4
        // Empirically verified against real Gmsh 4.14 output (the manual
        // abridges its list after 31; these came from block headers of
        // meshes Gmsh itself wrote).
        assert_eq!(super::gmsh_node_count(8), Some(3)); // Line3
        assert_eq!(super::gmsh_node_count(16), Some(8)); // Quad8
        assert_eq!(super::gmsh_node_count(36), Some(16)); // order-3 complete quad
        assert_eq!(super::gmsh_node_count(39), Some(12)); // order-3 incomplete quad
        assert_eq!(super::gmsh_node_count(51), Some(121)); // order-10 complete quad
        assert_eq!(super::gmsh_node_count(66), Some(11)); // order-10 edge
        assert_eq!(super::gmsh_node_count(90), Some(40)); // order-3 prism
        assert_eq!(super::gmsh_node_count(91), Some(75)); // order-4 prism
        assert_eq!(super::gmsh_node_count(118), Some(30)); // order-3 pyramid
        assert_eq!(super::gmsh_node_count(119), Some(55)); // order-4 pyramid
        // Variable-count and unlisted types stay unlisted.
        assert_eq!(super::gmsh_node_count(34), None);
        assert_eq!(super::gmsh_node_count(35), None);
        assert_eq!(super::gmsh_node_count(999), None);
    }

    #[test]
    fn every_listed_type_count_has_a_family_and_positive_count() {
        let listed = [
            15, 1, 8, 26, 27, 28, 62, 63, 64, 65, 66, 2, 9, 20, 21, 22, 23, 24, 25, 32, 33, 42,
            43, 44, 45, 46, 52, 53, 54, 55, 56, 3, 10, 16, 36, 37, 38, 39, 40, 41, 47, 48, 49,
            50, 51, 57, 58, 59, 60, 61, 4, 11, 29, 30, 31, 5, 12, 17, 92, 93, 6, 13, 18, 90, 91,
            7, 14, 19, 118, 119,
        ];
        assert_eq!(listed.len(), 70, "the table pins 70 verified type ids");
        for id in listed {
            let count = super::gmsh_node_count(id)
                .unwrap_or_else(|| panic!("type {id} must have a count"));
            assert!(count >= 1, "type {id} count must be positive");
            assert!(
                super::family_from_gmsh_type_id(id).is_some(),
                "type {id} has a count but no family"
            );
        }
    }

    #[test]
    fn serde_json_cycle_preserves_unknown_high_order_element() {
        let etype = ElementType::Unknown(11);
        let json = serde_json::to_string(&etype).expect("serialize element type");
        let back: ElementType = serde_json::from_str(&json).expect("deserialize element type");
        assert_eq!(back, etype, "the Unknown(id) form must survive a JSON cycle");

        let element = Element::new(7, ElementType::Unknown(16), vec![1, 2, 3, 4, 5, 6, 7, 8]);
        let json = serde_json::to_string(&element).expect("serialize element");
        let back: Element = serde_json::from_str(&json).expect("deserialize element");
        assert_eq!(back.id, 7);
        assert_eq!(back.etype, ElementType::Unknown(16));
        assert_eq!(back.physical_tag, None);
        assert_eq!(back.node_ids, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    }

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

        // 2-Dimensional (Surface/Face)
        assert_eq!(ElementType::Triangle3.dimension(), 2);
        assert_eq!(ElementType::Quad4.dimension(), 2);

        // 3-Dimensional (Volume/Region)
        assert_eq!(ElementType::Tetrahedron4.dimension(), 3);
        assert_eq!(ElementType::Hexahedron8.dimension(), 3);
        assert_eq!(ElementType::Prism6.dimension(), 3);
        assert_eq!(ElementType::Pyramid5.dimension(), 3);

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

        // Line: 2 nodes
        assert_eq!(ElementType::Line2.node_count(), 2);
        assert_eq!(ElementType::Line2.dimension(), 1);

        // Triangle: 3 nodes
        assert_eq!(ElementType::Triangle3.node_count(), 3);
        assert_eq!(ElementType::Triangle3.dimension(), 2);

        // Quad: 4 nodes
        assert_eq!(ElementType::Quad4.node_count(), 4);
        assert_eq!(ElementType::Quad4.dimension(), 2);

        // Tet: 4 nodes
        assert_eq!(ElementType::Tetrahedron4.node_count(), 4);
        assert_eq!(ElementType::Tetrahedron4.dimension(), 3);

        // Hex: 8 nodes
        assert_eq!(ElementType::Hexahedron8.node_count(), 8);
        assert_eq!(ElementType::Hexahedron8.dimension(), 3);

        // Prism: 6 nodes
        assert_eq!(ElementType::Prism6.node_count(), 6);
        assert_eq!(ElementType::Prism6.dimension(), 3);

        // Pyramid: 5 nodes
        assert_eq!(ElementType::Pyramid5.node_count(), 5);
        assert_eq!(ElementType::Pyramid5.dimension(), 3);
    }

    // ── P5: semantic high-order variants (orders 2-3, Gmsh type ids pinned
    // by the verified node-count table and the real-Gmsh probe) ──

    /// The 20 semantic high-order variants with their Gmsh MSH type ids.
    const HIGH_ORDER_SEMANTIC: &[(ElementType, i32)] = &[
        (ElementType::Line3, 8),
        (ElementType::Line4, 26),
        (ElementType::Triangle6, 9),
        (ElementType::Triangle9, 20),
        (ElementType::Triangle10, 21),
        (ElementType::Quad8, 16),
        (ElementType::Quad9, 10),
        (ElementType::Quad12, 39),
        (ElementType::Quad16, 36),
        (ElementType::Tetrahedron10, 11),
        (ElementType::Tetrahedron20, 29),
        (ElementType::Hexahedron20, 17),
        (ElementType::Hexahedron27, 12),
        (ElementType::Hexahedron64, 92),
        (ElementType::Prism15, 18),
        (ElementType::Prism18, 13),
        (ElementType::Prism40, 90),
        (ElementType::Pyramid13, 19),
        (ElementType::Pyramid14, 14),
        (ElementType::Pyramid30, 118),
    ];

    #[test]
    fn high_order_semantic_ids_round_trip_bidirectionally() {
        for (etype, id) in HIGH_ORDER_SEMANTIC {
            let back = ElementType::from_gmsh_type_id(*id);
            assert_eq!(
                back, *etype,
                "type {id} must read back as its semantic variant"
            );
            let id_back = back.to_gmsh_type_id();
            assert_eq!(
                id_back,
                Some(*id),
                "semantic variant {back:?} must map back to type {id}"
            );
            assert_ne!(*etype, ElementType::Unknown(*id));
        }
    }

    #[test]
    fn high_order_semantic_node_counts_match_the_gmsh_table() {
        for (etype, id) in HIGH_ORDER_SEMANTIC {
            let expected = super::gmsh_node_count(*id)
                .unwrap_or_else(|| panic!("type {id} must be in the verified table"));
            assert_eq!(
                etype.node_count(),
                expected,
                "{etype:?} node count must match the verified Gmsh table"
            );
        }
    }

    #[test]
    fn high_order_semantic_dimensions() {
        assert_eq!(ElementType::Line3.dimension(), 1);
        assert_eq!(ElementType::Line4.dimension(), 1);
        assert_eq!(ElementType::Triangle6.dimension(), 2);
        assert_eq!(ElementType::Triangle9.dimension(), 2);
        assert_eq!(ElementType::Triangle10.dimension(), 2);
        assert_eq!(ElementType::Quad8.dimension(), 2);
        assert_eq!(ElementType::Quad9.dimension(), 2);
        assert_eq!(ElementType::Quad12.dimension(), 2);
        assert_eq!(ElementType::Quad16.dimension(), 2);
        assert_eq!(ElementType::Tetrahedron10.dimension(), 3);
        assert_eq!(ElementType::Tetrahedron20.dimension(), 3);
        assert_eq!(ElementType::Hexahedron20.dimension(), 3);
        assert_eq!(ElementType::Hexahedron27.dimension(), 3);
        assert_eq!(ElementType::Hexahedron64.dimension(), 3);
        assert_eq!(ElementType::Prism15.dimension(), 3);
        assert_eq!(ElementType::Prism18.dimension(), 3);
        assert_eq!(ElementType::Prism40.dimension(), 3);
        assert_eq!(ElementType::Pyramid13.dimension(), 3);
        assert_eq!(ElementType::Pyramid14.dimension(), 3);
        assert_eq!(ElementType::Pyramid30.dimension(), 3);
    }

    #[test]
    fn high_order_semantic_volume_families_expose_canonical_faces_and_edges() {
        // Corner topology only (the established family-table convention):
        // face/edge tables describe the linear skeleton used for adjacency.
        assert_eq!(ElementType::Tetrahedron10.faces().len(), 4);
        assert_eq!(ElementType::Tetrahedron10.edges().len(), 6);
        assert_eq!(ElementType::Tetrahedron20.faces().len(), 4);
        assert_eq!(ElementType::Tetrahedron20.edges().len(), 6);
        assert_eq!(ElementType::Hexahedron20.faces().len(), 6);
        assert_eq!(ElementType::Hexahedron20.edges().len(), 12);
        assert_eq!(ElementType::Hexahedron27.faces().len(), 6);
        assert_eq!(ElementType::Hexahedron27.edges().len(), 12);
        assert_eq!(ElementType::Hexahedron64.faces().len(), 6);
        assert_eq!(ElementType::Hexahedron64.edges().len(), 12);
        assert_eq!(ElementType::Prism15.faces().len(), 5);
        assert_eq!(ElementType::Prism15.edges().len(), 9);
        assert_eq!(ElementType::Prism18.faces().len(), 5);
        assert_eq!(ElementType::Prism18.edges().len(), 9);
        assert_eq!(ElementType::Prism40.faces().len(), 5);
        assert_eq!(ElementType::Prism40.edges().len(), 9);
        assert_eq!(ElementType::Pyramid13.faces().len(), 5);
        assert_eq!(ElementType::Pyramid13.edges().len(), 8);
        assert_eq!(ElementType::Pyramid14.faces().len(), 5);
        assert_eq!(ElementType::Pyramid14.edges().len(), 8);
        assert_eq!(ElementType::Pyramid30.faces().len(), 5);
        assert_eq!(ElementType::Pyramid30.edges().len(), 8);
        // Surface families carry edge skeletons.
        assert_eq!(ElementType::Triangle6.edges().len(), 3);
        assert_eq!(ElementType::Triangle9.edges().len(), 3);
        assert_eq!(ElementType::Triangle10.edges().len(), 3);
        assert_eq!(ElementType::Quad8.edges().len(), 4);
        assert_eq!(ElementType::Quad9.edges().len(), 4);
        assert_eq!(ElementType::Quad12.edges().len(), 4);
        assert_eq!(ElementType::Quad16.edges().len(), 4);
        // Lines are their own skeleton.
        assert_eq!(ElementType::Line3.edges().len(), 1);
        assert_eq!(ElementType::Line4.edges().len(), 1);
    }

    #[test]
    fn orders_outside_the_semantic_set_stay_unknown() {
        // Order 4+ types keep the Unknown(id) fidelity form.
        assert_eq!(ElementType::from_gmsh_type_id(22), ElementType::Unknown(22)); // Tri12, order 4
        assert_eq!(ElementType::from_gmsh_type_id(30), ElementType::Unknown(30)); // Tet35, order 4
        assert_eq!(ElementType::from_gmsh_type_id(31), ElementType::Unknown(31)); // Tet56, order 5
        assert_eq!(ElementType::from_gmsh_type_id(93), ElementType::Unknown(93)); // Hex125, order 4
        assert_eq!(ElementType::from_gmsh_type_id(91), ElementType::Unknown(91)); // Prism75, order 4
        assert_eq!(ElementType::from_gmsh_type_id(119), ElementType::Unknown(119)); // Pyramid55, order 4
        assert_eq!(ElementType::from_gmsh_type_id(27), ElementType::Unknown(27)); // Line5, order 4
    }

    #[test]
    fn high_order_semantic_survives_serde_json_cycle() {
        let etype = ElementType::Tetrahedron10;
        let json = serde_json::to_string(&etype).expect("serialize");
        let back: ElementType = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(back, ElementType::Tetrahedron10);
    }
}
