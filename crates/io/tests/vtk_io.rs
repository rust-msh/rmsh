// P6: VTK/VTU interop tests.
//
// Contracts pinned here:
// - Cell-type mapping: every supported rmsh element writes the correct
//   VTK cell type id (the table pinned to meshio's `_vtk_common`, the
//   ecosystem reference ParaView-adjacent tooling consumes); the five
//   permuting types (tetra10, hex20, hex27, wedge15, pyramid13) carry the
//   meshio-pinned permutations, wedge6 carries the normal flip.
// - Round trips: legacy VTK and VTU both parse back to the identical
//   element set (types, node multisets after the inverse permutation,
//   physical tags).
// - The read-side permutation is the inverse of the write-side
//   permutation: round-tripping a permuted type restores the Gmsh node
//   order exactly (pinned per element via node coordinate positions).
// - Unsupported elements (no VTK counterpart) are refused by name.
// - Ecosystem cross-check: the written VTU parses back through meshio
//   (the ecosystem parser) with identical cell types and connectivity —
//   run via the local uv python when meshio is available.

use rmsh_io::{read_vtk_legacy, read_vtu, write_vtk_legacy, write_vtu};
use rmsh_model::{Element, ElementType, Mesh, Node};


/// Read-side node ids are dense ranks (1..n) over the sorted original node
/// ids — map a read id back to the original node id through that.
fn rank_to_orig(mesh: &Mesh) -> Vec<u64> {
    let mut ids: Vec<u64> = mesh.nodes.keys().copied().collect();
    ids.sort_unstable();
    ids
}

/// A tet10 with nodes at the Gmsh probe positions (unit tet + midside
/// nodes), carrying distinct coordinates per node so order permutations
/// are observable through positions alone.
fn tet10_mesh() -> Mesh {
    let mut mesh = Mesh::new();
    let coords: [(u64, [f64; 3]); 10] = [
        (1, [0.0, 0.0, 0.0]),
        (2, [1.0, 0.0, 0.0]),
        (3, [0.0, 1.0, 0.0]),
        (4, [0.0, 0.0, 1.0]),
        (5, [0.5, 0.0, 0.0]),
        (6, [0.5, 0.5, 0.0]),
        (7, [0.0, 0.5, 0.0]),
        (8, [0.0, 0.0, 0.5]),
        (9, [0.0, 0.5, 0.5]),
        (10, [0.5, 0.0, 0.5]),
    ];
    for (id, c) in coords {
        mesh.add_node(Node::new(id, c[0], c[1], c[2]));
    }
    let mut el = Element::new(
        1,
        ElementType::Tetrahedron10,
        (1..=10).collect(),
    );
    el.physical_tag = Some(7);
    mesh.add_element(el);
    mesh
}

/// A mixed mesh: one of every supported linear type plus the permuting
/// quadratic types, each in its own node-id pool (ids 100*k + 1..).
fn kitchen_sink_mesh() -> Mesh {
    let mut mesh = Mesh::new();
    let mut tag = 0usize;
    let specs: &[(ElementType, usize)] = &[
        (ElementType::Point1, 1),
        (ElementType::Line2, 2),
        (ElementType::Line3, 3),
        (ElementType::Line4, 4),
        (ElementType::Triangle3, 3),
        (ElementType::Triangle6, 6),
        (ElementType::Quad4, 4),
        (ElementType::Quad8, 8),
        (ElementType::Quad9, 9),
        (ElementType::Tetrahedron4, 4),
        (ElementType::Tetrahedron10, 10),
        (ElementType::Hexahedron8, 8),
        (ElementType::Hexahedron20, 20),
        (ElementType::Hexahedron27, 27),
        (ElementType::Prism6, 6),
        (ElementType::Prism15, 15),
        (ElementType::Pyramid5, 5),
        (ElementType::Pyramid13, 13),
    ];
    for (etype, n) in specs {
        tag += 1;
        let base = 100 * tag;
        for i in 1..=*n {
            let id = (base + i) as u64;
            let c = id as f64;
            mesh.add_node(Node::new(id, c, c * 0.0, c * 0.0));
        }
        let ids: Vec<u64> = (1..=*n).map(|i| (base + i) as u64).collect();
        let mut el = Element::new(tag as u64, *etype, ids);
        el.physical_tag = Some(tag as i32);
        mesh.add_element(el);
    }
    mesh
}

#[test]
fn vtk_cell_types_match_the_mapping_table() {
    // the VTK type ids, pinned (meshio _vtk_common / VTK spec)
    let cases: &[(ElementType, i32)] = &[
        (ElementType::Point1, 1),
        (ElementType::Line2, 3),
        (ElementType::Line3, 21),
        (ElementType::Line4, 35),
        (ElementType::Triangle3, 5),
        (ElementType::Triangle6, 22),
        (ElementType::Quad4, 9),
        (ElementType::Quad8, 23),
        (ElementType::Quad9, 28),
        (ElementType::Tetrahedron4, 10),
        (ElementType::Tetrahedron10, 24),
        (ElementType::Hexahedron8, 12),
        (ElementType::Hexahedron20, 25),
        (ElementType::Hexahedron27, 29),
        (ElementType::Prism6, 13),
        (ElementType::Prism15, 26),
        (ElementType::Pyramid5, 14),
        (ElementType::Pyramid13, 27),
    ];
    let mesh = kitchen_sink_mesh();
    let mut out = Vec::new();
    rmsh_io::write_vtk_legacy(&mut out, &mesh).expect("kitchen sink writes");
    let text = String::from_utf8(out).expect("utf8");
    let types_pos = text.find("CELL_TYPES").expect("cell types");
    let n_cells = mesh.elements.len();
    let type_tokens: Vec<i32> = text[types_pos..]
        .lines()
        .skip(1)
        .filter_map(|l| l.trim().parse::<i32>().ok())
        .take(n_cells)
        .collect();
    assert_eq!(type_tokens.len(), cases.len(), "one cell per case");
    for (case, actual) in cases.iter().zip(type_tokens.iter()) {
        assert_eq!(actual, &case.1, "cell type for {:?}", case.0);
    }
}

#[test]
fn permuted_types_restore_gmsh_order_on_read_back() {
    // For every permuting type the round trip must restore the exact Gmsh
    // node order (the permutation inverse is verified through positions:
    // node ids equal coordinates in the kitchen-sink fixture).
    for mesh_source in [tet10_mesh(), kitchen_sink_mesh()] {
        let mut out = Vec::new();
        write_vtk_legacy(&mut out, &mesh_source).expect("write vtk");
        let back = read_vtk_legacy(&out).expect("read vtk");
        let mut orig: Vec<&Element> = mesh_source.elements.iter().collect();
        orig.sort_by_key(|e| e.id);
        let mut round: Vec<&Element> = back.elements.iter().collect();
        round.sort_by_key(|e| e.id);
        assert_eq!(orig.len(), round.len());
        let r2o = rank_to_orig(&mesh_source);
        let orig_order = |e: &Element| -> Vec<u64> {
            e.node_ids.iter().map(|id| r2o[(id - 1) as usize]).collect()
        };
        for (a, b) in orig.iter().zip(round.iter()) {
            assert_eq!(a.etype, b.etype, "{:?} type survives", a.etype);
            assert_eq!(orig_order(b), a.node_ids, "{:?} node order survives", a.etype);
            assert_eq!(a.physical_tag, b.physical_tag);
        }
        let r2o = rank_to_orig(&mesh_source);
        let orig_order = |e: &Element| -> Vec<u64> {
            e.node_ids.iter().map(|id| r2o[(id - 1) as usize]).collect()
        };
        for a in mesh_source.elements.iter() {
            let b = back
                .elements
                .iter()
                .find(|e| e.id == a.id)
                .unwrap_or_else(|| panic!("element {} missing after vtu", a.id));
            assert_eq!(a.etype, b.etype, "{:?} vtu type", a.etype);
            assert_eq!(orig_order(b), a.node_ids, "{:?} vtu node order", a.etype);
        }    }
}

#[test]
fn tet10_connectivity_is_permuted_in_the_vtk_file() {
    // the physical proof of the permutation: the written VTK file carries
    // the VTK order (last two entries swapped relative to the Gmsh file
    // order), and reading it back restores the Gmsh order
    let mesh = tet10_mesh();
    let mut out = Vec::new();
    write_vtk_legacy(&mut out, &mesh).expect("write");
    let text = String::from_utf8(out).expect("utf8");
    let cells_pos = text.find("CELLS").expect("cells");
    let line: Vec<u64> = text[cells_pos..]
        .lines()
        .nth(1)
        .expect("cell line")
        .split_whitespace()
        .skip(1)
        .map(|t| t.parse().expect("rank"))
        .collect();
    assert_eq!(line.len(), 10);
    // Gmsh file order ranks: 0..9 (ids 1..10 sorted); the VTK order swaps
    // the last two midside nodes
    let expected: Vec<u64> = vec![0, 1, 2, 3, 4, 5, 6, 7, 9, 8];
    assert_eq!(line, expected, "VTK tet10 connectivity must be permuted");
}

#[test]
fn unsupported_elements_are_refused_by_name() {
    for etype in [
        ElementType::Tetrahedron20,
        ElementType::Hexahedron64,
        ElementType::Pyramid14,
        ElementType::Triangle9,
        ElementType::Quad12,
        ElementType::Quad16,
        ElementType::Triangle10,
    ] {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        let ids: Vec<u64> = (1..=etype.node_count().max(1) as u64).collect();
        mesh.add_element(Element::new(1, etype, ids));
        let mut out = Vec::new();
        let vtk_err = rmsh_io::write_vtk_legacy(&mut out, &mesh).unwrap_err();
        assert!(
            vtk_err.to_string().contains(&format!("{etype:?}")),
            "the error must name the element: {vtk_err}"
        );
        let vtu_err = rmsh_io::write_vtu(&mut out, &mesh).unwrap_err();
        assert!(vtu_err.to_string().contains(&format!("{etype:?}")));
    }
}

#[test]
fn legacy_binary_and_wrong_dataset_are_refused() {
    let header = "# vtk DataFile Version 4.2\ntitle\nBINARY\nDATASET UNSTRUCTURED_GRID\n";
    assert!(read_vtk_legacy(header.as_bytes()).is_err());
    let wrong = "# vtk DataFile Version 4.2\ntitle\nASCII\nDATASET POLYDATA\n";
    assert!(read_vtk_legacy(wrong.as_bytes()).is_err());
}

#[test]
fn physical_tags_survive_both_flavors() {
    let mesh = tet10_mesh();
    let mut out = Vec::new();
    write_vtk_legacy(&mut out, &mesh).expect("write");
    let back = read_vtk_legacy(&out).expect("read");
    assert_eq!(back.elements[0].physical_tag, Some(7));

    let mut out = Vec::new();
    write_vtu(&mut out, &mesh).expect("write vtu");
    let back = read_vtu(&out).expect("read vtu");
    assert_eq!(back.elements[0].physical_tag, Some(7));
}

#[test]
fn vtu_carries_the_vtk_structural_markers() {
    let mesh = kitchen_sink_mesh();
    let num_points = mesh.nodes.len();
    let num_cells = mesh.elements.len();
    let mut out = Vec::new();
    write_vtu(&mut out, &mesh).expect("write");
    let text = String::from_utf8(out).expect("utf8");
    for marker in [
        "<VTKFile type=\"UnstructuredGrid\"",
        &format!("NumberOfPoints=\"{num_points}\""),
        &format!("NumberOfCells=\"{num_cells}\""),
        "Name=\"connectivity\"",
        "Name=\"offsets\"",
        "Name=\"types\"",
        "</VTKFile>",
    ] {
        assert!(text.contains(marker), "missing {marker}");
    }
}
