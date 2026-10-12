// P1 tranche 2: high-order element fidelity + $Periodic semantics.
//
// Contracts pinned here (all layouts verified against real Gmsh 4.14 output
// generated on this machine, see the layout comments in msh.rs / msh_periodic.rs):
// - A mesh carrying high-order elements (Gmsh types 11 = Tet10, 16 = Quad8,
//   8 = Line3) round trips field-exact through all four MSH flavors.
// - Node-count guards refuse known-type/wrong-count elements (v2 ASCII guard,
//   v4 binary structural check, writer-side validate_mesh).
// - $Periodic links round trip field-exact through all four flavors; a
//   periodic-free mesh stays byte-identical to its pre-feature output.

use rmsh_io::{
    load_msh_from_bytes, write_msh_v2, write_msh_v2_binary, write_msh_v4, write_msh_v4_binary,
    MshError,
};
use rmsh_model::{Element, ElementType, Mesh, Node, PeriodicLink};

/// A mesh with one Tet10 (type 11), one Quad8 (type 16) and one Line3
/// (type 8). Node pools are disjoint per element so the writer-side
/// node-existence check passes for every flavor.
fn high_order_mesh() -> Mesh {
    let mut mesh = Mesh::new();
    // Tet10 nodes 1..=10
    for id in 1..=10u64 {
        let c = id as f64;
        mesh.add_node(Node::new(id, c, 0.0, 0.0));
    }
    // Quad8 nodes 11..=18
    for id in 11..=18u64 {
        let c = id as f64;
        mesh.add_node(Node::new(id, 0.0, c, 0.0));
    }
    // Line3 nodes 19..=21
    for id in 19..=21u64 {
        let c = id as f64;
        mesh.add_node(Node::new(id, 0.0, 0.0, c));
    }

    let mut tet10 = Element::new(
        1,
        ElementType::Unknown(11),
        vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    );
    tet10.physical_tag = Some(21);
    mesh.add_element(tet10);

    let mut quad8 = Element::new(2, ElementType::Unknown(16), vec![11, 12, 13, 14, 15, 16, 17, 18]);
    quad8.physical_tag = Some(22);
    mesh.add_element(quad8);

    let mut line3 = Element::new(3, ElementType::Unknown(8), vec![19, 20, 21]);
    line3.physical_tag = Some(23);
    mesh.add_element(line3);

    mesh
}

fn assert_high_order_mesh_agrees(loaded: &Mesh) {
    assert_eq!(loaded.node_count(), 21, "all 21 nodes must survive");
    assert_eq!(loaded.element_count(), 3, "all 3 elements must survive");

    let mut elements: Vec<_> = loaded.elements.iter().collect();
    elements.sort_by_key(|e| e.id);

    let tet10 = elements[0];
    assert_eq!(tet10.etype, ElementType::Unknown(11), "Tet10 must stay type 11");
    assert_eq!(tet10.node_ids, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    assert_eq!(tet10.physical_tag, Some(21));

    let quad8 = elements[1];
    assert_eq!(quad8.etype, ElementType::Unknown(16), "Quad8 must stay type 16");
    assert_eq!(quad8.node_ids, vec![11, 12, 13, 14, 15, 16, 17, 18]);
    assert_eq!(quad8.physical_tag, Some(22));

    let line3 = elements[2];
    assert_eq!(line3.etype, ElementType::Unknown(8), "Line3 must stay type 8");
    assert_eq!(line3.node_ids, vec![19, 20, 21]);
    assert_eq!(line3.physical_tag, Some(23));
}

#[test]
fn high_order_mesh_round_trips_all_four_flavors() {
    let mesh = high_order_mesh();

    let mut v2 = Vec::new();
    write_msh_v2(&mut v2, &mesh).expect("v2 ascii write of high-order mesh");
    let loaded = load_msh_from_bytes(&v2).expect("v2 ascii read");
    assert_high_order_mesh_agrees(&loaded);

    let mut v2b = Vec::new();
    write_msh_v2_binary(&mut v2b, &mesh).expect("v2 binary write of high-order mesh");
    let loaded = load_msh_from_bytes(&v2b).expect("v2 binary read");
    assert_high_order_mesh_agrees(&loaded);

    let mut v4 = Vec::new();
    write_msh_v4(&mut v4, &mesh).expect("v4 ascii write of high-order mesh");
    let loaded = load_msh_from_bytes(&v4).expect("v4 ascii read");
    assert_high_order_mesh_agrees(&loaded);

    let mut v4b = Vec::new();
    write_msh_v4_binary(&mut v4b, &mesh).expect("v4 binary write of high-order mesh");
    let loaded = load_msh_from_bytes(&v4b).expect("v4 binary read");
    assert_high_order_mesh_agrees(&loaded);
}

#[test]
fn v2_ascii_refuses_wrong_node_count_for_known_high_order_type() {
    // Type 11 (Tet10) declared with only 9 node ids: the ASCII node-count
    // guard must refuse it instead of building a structurally corrupt
    // element (high-order types are no longer exempt from the guard).
    let msh_data = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n9\n1 0 0 0\n2 1 0 0\n3 0 1 0\n\
                    4 0 0 1\n5 1 1 0\n6 1 0 1\n7 0 1 1\n8 0.5 0 0\n9 0 0.5 0\n$EndNodes\n\
                    $Elements\n1\n1 11 0 1 2 3 4 5 6 7 8 9\n$EndElements\n";
    let err = parse_expect_err(msh_data.as_bytes());
    let message = err.to_string();
    assert!(
        message.contains("10 node(s)") && message.contains("9"),
        "the error must name the expected and actual node counts: {message}"
    );
}

#[test]
fn v2_binary_refuses_type_outside_the_node_count_table() {    // The chunked v2 binary layout has no per-file node count: the table IS
    // the contract. A type with no entry anywhere must be refused by name.
    let mut data = Vec::<u8>::new();
    data.extend_from_slice(b"$MeshFormat\n2.2 1 8\n");
    data.extend_from_slice(&1_i32.to_le_bytes());
    data.extend_from_slice(&1.0f64.to_le_bytes());
    data.extend_from_slice(b"\n$EndMeshFormat\n");
    data.extend_from_slice(b"$Nodes\n1\n");
    data.extend_from_slice(&1_i32.to_le_bytes());
    data.extend_from_slice(&0.0f64.to_le_bytes());
    data.extend_from_slice(&0.0f64.to_le_bytes());
    data.extend_from_slice(&0.0f64.to_le_bytes());
    data.extend_from_slice(b"\n$EndNodes\n");
    data.extend_from_slice(b"$Elements\n1\n");
    data.extend_from_slice(&999_i32.to_le_bytes()); // no such Gmsh type
    data.extend_from_slice(&1_i32.to_le_bytes());
    data.extend_from_slice(&0_i32.to_le_bytes());
    data.extend_from_slice(&1_i32.to_le_bytes());
    data.extend_from_slice(&1_i32.to_le_bytes());
    data.extend_from_slice(b"$EndElements\n");

    let err = load_msh_from_bytes(&data).expect_err("type 999 must be refused");
    assert!(
        err.to_string().contains("999"),
        "the error must name the offending type: {err}"
    );
}

#[test]
fn v2_binary_loads_the_compact_gmsh_magic() {
    // Real Gmsh 4.14 writes ONLY the i32 1 after the format line (no f64
    // half); the reader must accept that compact form byte-for-byte.
    let mut data = Vec::<u8>::new();
    data.extend_from_slice(b"$MeshFormat\n2.2 1 8\n");
    data.extend_from_slice(&1_i32.to_le_bytes());
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
    data.extend_from_slice(&1_i32.to_le_bytes()); // line2
    data.extend_from_slice(&1_i32.to_le_bytes());
    data.extend_from_slice(&0_i32.to_le_bytes()); // no tags
    data.extend_from_slice(&1_i32.to_le_bytes());
    data.extend_from_slice(&1_i32.to_le_bytes());
    data.extend_from_slice(&2_i32.to_le_bytes());
    data.extend_from_slice(b"$EndElements\n");

    let mesh = load_msh_from_bytes(&data).expect("compact-magic v2 binary should parse");
    assert_eq!(mesh.node_count(), 2);
    assert_eq!(mesh.element_count(), 1);
    assert_eq!(mesh.elements[0].node_ids, vec![1, 2]);
}

#[test]
fn v4_binary_loads_the_compact_gmsh_magic() {
    // Same compact form for the v4.1 binary flavor: strip the f64 half of
    // a written file's magic and the reader must still load it.
    let mut mesh = Mesh::new();
    for (id, xyz) in [(1, [0.0, 0.0, 0.0]), (2, [1.0, 0.0, 0.0]), (3, [0.0, 1.0, 0.0])] {
        mesh.add_node(Node::new(id, xyz[0], xyz[1], xyz[2]));
    }
    let mut tri = Element::new(1, ElementType::Triangle3, vec![1, 2, 3]);
    tri.physical_tag = Some(7);
    mesh.add_element(tri);

    let mut bytes = Vec::new();
    write_msh_v4_binary(&mut bytes, &mesh).expect("write");
    let magic = b"$MeshFormat\n4.1 1 8\n".len();
    let mut compact = Vec::new();
    compact.extend_from_slice(&bytes[..magic + 4]);
    compact.extend_from_slice(&bytes[magic + 12..]);
    let loaded = load_msh_from_bytes(&compact).expect("compact-magic v4 binary should parse");
    assert_eq!(loaded.node_count(), 3);
    assert_eq!(loaded.element_count(), 1);
    assert_eq!(loaded.elements[0].etype, ElementType::Triangle3);
    assert_eq!(loaded.elements[0].physical_tag, Some(7));
    assert_eq!(loaded.elements[0].node_ids, vec![1, 2, 3]);
}

#[test]
fn validate_mesh_refuses_mismatched_unknown_element_at_write() {
    let mut mesh = Mesh::new();
    for id in 1..=9u64 {
        let c = id as f64;
        mesh.add_node(Node::new(id, c, 0.0, 0.0));
    }
    // Type 11 carries 10 nodes; this element only has 9.
    mesh.add_element(Element::new(
        1,
        ElementType::Unknown(11),
        vec![1, 2, 3, 4, 5, 6, 7, 8, 9],
    ));

    let mut out = Vec::new();
    let err = write_msh_v2(&mut out, &mesh).expect_err("count mismatch must be refused");
    let message = err.to_string();
    assert!(
        message.contains("10 node(s)") && message.contains("9"),
        "the error must name the expected and actual node counts: {message}"
    );
}

/// Two periodic links: link A carries the 16-value affine transform only,
/// link B carries correspondence pairs only (no affine). Same shape as the
/// real Gmsh 4.14 fixtures generated for this tranche.
fn periodic_mesh() -> Mesh {
    let mut mesh = Mesh::new();
    for id in 1..=6u64 {
        let c = id as f64;
        mesh.add_node(Node::new(id, c, 0.0, 0.0));
    }
    let mut tri = Element::new(1, ElementType::Triangle3, vec![1, 2, 3]);
    tri.physical_tag = Some(5);
    mesh.add_element(tri);
    mesh.periodic.push(PeriodicLink {
        entity_dim: 2,
        master: 1,
        slave: 2,
        affine: vec![
            1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ],
        correspondence: Vec::new(),
    });
    mesh.periodic.push(PeriodicLink {
        entity_dim: 1,
        master: 1,
        slave: 5,
        affine: Vec::new(),
        correspondence: vec![(6, 2), (13, 9), (14, 10)],
    });
    mesh
}

fn assert_periodic_links_agree(loaded: &Mesh) {
    assert_eq!(loaded.periodic.len(), 2, "both links must survive");
    let (a, b) = (&loaded.periodic[0], &loaded.periodic[1]);
    assert_eq!(a.entity_dim, 2);
    assert_eq!(a.slave, 2);
    assert_eq!(a.master, 1);
    assert_eq!(a.affine.len(), 16);
    assert_eq!(a.affine[0], 1.0);
    assert_eq!(a.affine[3], 1.0);
    assert_eq!(a.affine[5], 1.0);
    assert_eq!(a.affine[10], 1.0);
    assert_eq!(a.affine[15], 1.0);
    assert!(a.correspondence.is_empty());
    assert_eq!(b.entity_dim, 1);
    assert_eq!(b.slave, 5);
    assert_eq!(b.master, 1);
    assert!(b.affine.is_empty());
    assert_eq!(b.correspondence, vec![(6, 2), (13, 9), (14, 10)]);
}

#[test]
fn periodic_links_round_trip_all_four_flavors() {
    let mesh = periodic_mesh();

    let mut v2 = Vec::new();
    write_msh_v2(&mut v2, &mesh).expect("v2 ascii write with periodic");
    let loaded = load_msh_from_bytes(&v2).expect("v2 ascii read");
    assert_periodic_links_agree(&loaded);

    let mut v2b = Vec::new();
    write_msh_v2_binary(&mut v2b, &mesh).expect("v2 binary write with periodic");
    let loaded = load_msh_from_bytes(&v2b).expect("v2 binary read");
    assert_periodic_links_agree(&loaded);

    let mut v4 = Vec::new();
    write_msh_v4(&mut v4, &mesh).expect("v4 ascii write with periodic");
    let loaded = load_msh_from_bytes(&v4).expect("v4 ascii read");
    assert_periodic_links_agree(&loaded);

    let mut v4b = Vec::new();
    write_msh_v4_binary(&mut v4b, &mesh).expect("v4 binary write with periodic");
    let loaded = load_msh_from_bytes(&v4b).expect("v4 binary read");
    assert_periodic_links_agree(&loaded);
}

#[test]
fn v2_binary_periodic_section_is_ascii_exactly_like_gmsh() {
    // Gmsh keeps the MSH2 $Periodic section in ASCII even inside binary
    // files (only $Nodes/$Elements payloads go binary). Pin those bytes.
    let mesh = periodic_mesh();
    let mut bytes = Vec::new();
    write_msh_v2_binary(&mut bytes, &mesh).expect("write");
    let at = find_subslice(&bytes, b"$Periodic").expect("v2 binary output must carry $Periodic");
    let end = find_subslice(&bytes, b"$EndPeriodic").expect("missing $EndPeriodic")
        + b"$EndPeriodic".len()
        + 1; // the newline after the marker
    let section = &bytes[at..end];
    let expected = b"$Periodic\n2\n2 2 1\nAffine 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n0\n\
                     1 5 1\n3\n6 2\n13 9\n14 10\n$EndPeriodic\n";
    assert_eq!(
        std::str::from_utf8(section).expect("v2 binary $Periodic stays ASCII"),
        std::str::from_utf8(expected).expect("literal is utf8"),
        "v2 binary $Periodic must be byte-ASCII"
    );
}

#[test]
fn v4_binary_periodic_section_has_binary_payload() {
    // The MSH4 binary $Periodic section keeps its header lines ASCII and
    // writes payloads little-endian binary: u64 numLinks; per link i32 dim,
    // i32 slave, i32 master, u64 numAffine, f64 values, u64 numPairs,
    // u64 node pairs.
    let mut mesh = Mesh::new();
    for id in 1..=3u64 {
        mesh.add_node(Node::new(id, id as f64, 0.0, 0.0));
    }
    mesh.periodic.push(PeriodicLink {
        entity_dim: 2,
        master: 1,
        slave: 2,
        affine: vec![1.0, 0.0, 0.0, 1.0],
        correspondence: vec![(3, 1)],
    });

    let mut bytes = Vec::new();
    write_msh_v4_binary(&mut bytes, &mesh).expect("write");
    let at = find_subslice(&bytes, b"$Periodic")
        .expect("v4 binary output must contain a $Periodic section");

    let mut expected = Vec::<u8>::new();
    expected.extend_from_slice(b"$Periodic\n");
    expected.extend_from_slice(&1_u64.to_le_bytes());
    expected.extend_from_slice(&2_i32.to_le_bytes());
    expected.extend_from_slice(&2_i32.to_le_bytes());
    expected.extend_from_slice(&1_i32.to_le_bytes());
    expected.extend_from_slice(&4_u64.to_le_bytes()); // numAffine
    expected.extend_from_slice(&1.0f64.to_le_bytes());
    expected.extend_from_slice(&0.0f64.to_le_bytes());
    expected.extend_from_slice(&0.0f64.to_le_bytes());
    expected.extend_from_slice(&1.0f64.to_le_bytes());
    expected.extend_from_slice(&1_u64.to_le_bytes()); // numCorrespondingNodes
    expected.extend_from_slice(&3_u64.to_le_bytes()); // slave node
    expected.extend_from_slice(&1_u64.to_le_bytes()); // master node
    expected.extend_from_slice(b"\n$EndPeriodic\n");

    assert_eq!(
        &bytes[at..at + expected.len()],
        &expected[..],
        "v4 binary $Periodic payload must be byte-exact little-endian"
    );
}

/// The full $Periodic section of a real Gmsh 4.14 MSH v2.2 file (generated
/// from two translated unit squares, Mesh.SavePeriodic behavior), verbatim.
const GMSH_V2_PERIODIC: &str = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Periodic\n9\n\
0 5 1\nAffine 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n1\n5 1\n\
0 6 2\nAffine 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n1\n6 2\n\
0 7 3\nAffine 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n1\n7 3\n\
0 8 4\nAffine 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n1\n8 4\n\
1 5 1\nAffine 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n3\n5 1\n6 2\n13 9\n\
1 6 2\nAffine 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n3\n6 2\n7 3\n14 10\n\
1 7 3\nAffine 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n3\n7 3\n8 4\n15 11\n\
1 8 4\nAffine 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n3\n5 1\n8 4\n16 12\n\
2 2 1\nAffine 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n4\n5 1\n6 2\n7 3\n8 4\n\
$EndPeriodic\n";

/// The full $Periodic section of a real Gmsh 4.14 MSH v4.1 file, verbatim.
const GMSH_V4_PERIODIC: &str = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Periodic\n9\n\
0 5 1\n16 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n1\n5 1\n\
0 6 2\n16 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n1\n6 2\n\
0 7 3\n16 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n1\n7 3\n\
0 8 4\n16 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n1\n8 4\n\
1 5 1\n16 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n3\n5 1\n6 2\n13 9\n\
1 6 2\n16 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n3\n6 2\n7 3\n14 10\n\
1 7 3\n16 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n3\n7 3\n8 4\n15 11\n\
1 8 4\n16 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n3\n5 1\n8 4\n16 12\n\
2 2 1\n16 1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1\n4\n5 1\n6 2\n7 3\n8 4\n\
$EndPeriodic\n";

fn assert_real_gmsh_links(loaded: &Mesh) {
    assert_eq!(loaded.periodic.len(), 9, "Gmsh wrote 9 periodic links");
    // Slave-first ordering, exactly as Gmsh writes it.
    let first = &loaded.periodic[0];
    assert_eq!(first.entity_dim, 0);
    assert_eq!(first.slave, 5);
    assert_eq!(first.master, 1);
    assert_eq!(first.affine.len(), 16);
    assert_eq!(first.correspondence, vec![(5, 1)]);
    // The surface-level link closes the mesh-wide correspondence.
    let last = &loaded.periodic[8];
    assert_eq!(last.entity_dim, 2);
    assert_eq!(last.slave, 2);
    assert_eq!(last.master, 1);
    assert_eq!(last.correspondence, vec![(5, 1), (6, 2), (7, 3), (8, 4)]);
}

#[test]
fn parses_real_gmsh_v2_periodic_section_verbatim() {
    let loaded = load_msh_from_bytes(GMSH_V2_PERIODIC.as_bytes()).expect("parse gmsh v2 fixture");
    assert_real_gmsh_links(&loaded);
}

#[test]
fn parses_real_gmsh_v4_periodic_section_verbatim() {
    let loaded = load_msh_from_bytes(GMSH_V4_PERIODIC.as_bytes()).expect("parse gmsh v4 fixture");
    assert_real_gmsh_links(&loaded);
}

fn sample_periodic_free_mesh() -> Mesh {
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

#[test]
fn periodic_free_v2_output_is_byte_stable() {
    // A mesh without periodic links must produce exactly the same bytes as
    // before the feature existed (no $Periodic section is emitted).
    let mut output = Vec::new();
    write_msh_v2(&mut output, &sample_periodic_free_mesh()).expect("write");
    let expected = "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n1\n2 11 \"surface\"\n\
                    $EndPhysicalNames\n$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n$EndNodes\n\
                    $Elements\n2\n1 2 2 11 0 1 2 3\n2 4 0 1 2 3 4\n$EndElements\n";
    assert_eq!(String::from_utf8(output).expect("utf8"), expected);
}

fn parse_expect_err(data: &[u8]) -> MshError {
    match load_msh_from_bytes(data) {
        Err(err) => err,
        Ok(_) => panic!("expected a parse error"),
    }
}

fn find_subslice(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack.windows(needle.len()).position(|w| w == needle)
}
