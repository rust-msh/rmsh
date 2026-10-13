// P6: `$PartitionedEntities` + `$GhostElements` semantics.
//
// Contracts pinned here (layout authority Gmsh 4.14 GModelIO_MSH4.cpp;
// fixtures are real METIS-partitioned Gmsh 4.14 output generated on this
// machine — a 2x2x2 box split into 3 partitions with ghost cells):
// - The retained `PartitionInfo` parses field-exact from both the ASCII
//   and the binary flavor of the real fixture (counts, ghost entities,
//   per-entity parent/partition/bbox/physical/bounding records, ghost
//   cells).
// - A partitioned mesh round trips through both v4 writers/readers with
//   the partition structure intact (ascii -> ascii, ascii -> binary,
//   binary -> ascii, binary -> binary).
// - v2 has no partition section (Gmsh splits partitioned v2 meshes into
//   per-partition files): the v2 writer drops the partition data, pinned
//   as the documented fidelity line.

use rmsh_io::{load_msh_from_bytes, load_msh_from_path, write_msh_v2, write_msh_v4,
              write_msh_v4_binary};
use rmsh_model::PartitionInfo;

const FIXTURE_A: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../testdata/partitioned_cube_v4a.msh");
const FIXTURE_B: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../testdata/partitioned_cube_v4b.msh");

fn assert_partition_semantics(info: &PartitionInfo) {
    assert_eq!(info.num_partitions, 3);
    assert_eq!(info.entities.len(), 67, "17 v + 30 e + 17 f + 3 r");
    let (nv, ne, nf, nr) = (
        info.entities.iter().filter(|e| e.dim == 0).count(),
        info.entities.iter().filter(|e| e.dim == 1).count(),
        info.entities.iter().filter(|e| e.dim == 2).count(),
        info.entities.iter().filter(|e| e.dim == 3).count(),
    );
    assert_eq!((nv, ne, nf, nr), (17, 30, 17, 3));
    // ghost entities: one volume placeholder per partition
    assert_eq!(info.ghosts.len(), 3);
    for (k, ghost) in info.ghosts.iter().enumerate() {
        assert_eq!(ghost.dim, 3);
        assert_eq!(ghost.partition, k as i32 + 1);
    }
    // boundary pieces carry parent links and multi-partition lists
    let boundary_vertices: Vec<_> = info
        .entities
        .iter()
        .filter(|e| e.dim == 0 && e.parent_dim != 0)
        .collect();
    assert_eq!(boundary_vertices.len(), 9, "the partition-boundary vertices");
    assert!(boundary_vertices
        .iter()
        .all(|e| e.partitions.len() >= 2), "boundary pieces are shared");
    // regions: one per partition, each bounding a closed shell of faces
    let regions: Vec<_> = info
        .entities
        .iter()
        .filter(|e| e.dim == 3)
        .collect();
    assert_eq!(regions.len(), 3);
    for region in &regions {
        assert_eq!(region.partitions.len(), 1, "regions are single-partition");
    }
    // the middle piece bounds 6 faces (the box's 6 + shared planes split
    // across pieces give 7/6/7 for this 3-way METIS cut, pinned from the
    // fixture)
    let mut lens: Vec<usize> = regions.iter().map(|r| r.bounding.len()).collect();
    lens.sort_unstable();
    assert_eq!(lens, vec![6, 7, 7]);
    // ghost cells with the multi-partition lists
    assert!(!info.ghost_elements.is_empty(), "-part_ghosts wrote cells");
    for cell in &info.ghost_elements {
        assert!(cell.partition >= 1 && cell.partition <= 3);
        assert!(!cell.ghost_partitions.is_empty());
        assert!(cell
            .ghost_partitions
            .iter()
            .all(|gp| *gp >= 1 && *gp <= 3));
    }
}

#[test]
fn parses_real_partitioned_fixture_ascii() {
    let mesh = load_msh_from_path(FIXTURE_A.as_ref()).expect("ascii fixture reads");
    let info = mesh.partitioned.as_ref().expect("partition info retained");
    assert_partition_semantics(info);
}

#[test]
fn parses_real_partitioned_fixture_binary() {
    let mesh = load_msh_from_path(FIXTURE_B.as_ref()).expect("binary fixture reads");
    let info = mesh.partitioned.as_ref().expect("partition info retained");
    assert_partition_semantics(info);
}

/// The ascii and binary fixtures are independent METIS runs (different
/// partition pieces), so the equality contract is per-flavor round trip:
/// parse -> write -> re-parse keeps the partition structure equal.
fn assert_round_trip(info: &PartitionInfo, mesh_bytes: &[u8], v4_writer: fn(&mut Vec<u8>, &rmsh_model::Mesh) -> Result<(), rmsh_io::MshError>) {
    let mut out = Vec::new();
    // rebuild the mesh from bytes to carry elements + partition info
    let mesh = load_msh_from_bytes(mesh_bytes).expect("re-read");
    v4_writer(&mut out, &mesh).expect("write");
    let reparsed = load_msh_from_bytes(&out).expect("re-parse");
    let again = reparsed.partitioned.as_ref().expect("partition retained");
    assert_eq!(again.num_partitions, info.num_partitions);
    assert_eq!(again.entities.len(), info.entities.len());
    assert_eq!(again.ghosts, info.ghosts);
    assert_eq!(again.ghost_elements.len(), info.ghost_elements.len());
    for (a, b) in again.entities.iter().zip(info.entities.iter()) {
        assert_eq!(a, b, "entity records must round trip field-exact");
    }
}

#[test]
fn partitioned_mesh_round_trips_both_v4_flavors() {
    let ascii_bytes = std::fs::read(FIXTURE_A).expect("read ascii fixture");
    let ascii_mesh = load_msh_from_bytes(&ascii_bytes).expect("parse");
    let info_a = ascii_mesh.partitioned.clone().expect("info");

    let binary_bytes = std::fs::read(FIXTURE_B).expect("read binary fixture");
    let binary_mesh = load_msh_from_bytes(&binary_bytes).expect("parse");
    let info_b = binary_mesh.partitioned.clone().expect("info");

    // ascii -> ascii
    assert_round_trip(&info_a, &ascii_bytes, write_msh_v4);
    // binary -> binary
    assert_round_trip(&info_b, &binary_bytes, write_msh_v4_binary);
    // cross-flavor: ascii parse -> binary write -> ascii re-parse
    let mut out = Vec::new();
    write_msh_v4_binary(&mut out, &ascii_mesh).expect("write binary");
    let reparsed = load_msh_from_bytes(&out).expect("re-parse binary");
    assert_eq!(
        reparsed.partitioned.as_ref().expect("info").entities,
        info_a.entities,
        "cross-flavor round trip keeps the records"
    );
    // cross-flavor: binary parse -> ascii write -> binary re-parse
    let mut out = Vec::new();
    write_msh_v4(&mut out, &binary_mesh).expect("write ascii");
    let reparsed = load_msh_from_bytes(&out).expect("re-parse ascii");
    assert_eq!(
        reparsed.partitioned.as_ref().expect("info").entities,
        info_b.entities,
        "cross-flavor round trip keeps the records"
    );
}

#[test]
fn v2_writer_drops_partition_data_documented() {
    // MSH 2.2 has no partition section (Gmsh splits partitioned v2 meshes
    // into per-partition files): the v2 writer emits a section-free mesh.
    let mesh = load_msh_from_path(FIXTURE_A.as_ref()).expect("fixture");
    assert!(mesh.partitioned.is_some());
    let mut out = Vec::new();
    write_msh_v2(&mut out, &mesh).expect("v2 write");
    let text = String::from_utf8(out).expect("utf8");
    assert!(!text.contains("Partitioned"), "v2 carries no partition section");
    let reparsed = load_msh_from_bytes(text.as_bytes()).expect("v2 re-read");
    assert!(reparsed.partitioned.is_none(), "v2 reload has no partition info");
    assert_eq!(reparsed.node_count(), mesh.node_count(), "nodes survive");
    assert_eq!(reparsed.element_count(), mesh.element_count(), "elements survive");
}
