//! Integration tests for the rcad-scene assembly/instance tree.

use std::sync::Arc;

use glam::{DAffine3, DVec3};
use rcad_kernel::BRep;
use rcad_modeling::{make_box_brep, make_sphere_brep};
use rcad_scene::assembly::{Assembly, AssemblyNode, NodeContent, assembly_from_parts};

// ─── helpers ──────────────────────────────────────────────────────────────────

fn face_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum()
}

fn make_box() -> BRep {
    make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap()
}

// ─── tests ────────────────────────────────────────────────────────────────────

/// Two-level tree: parent at (1,0,0), child at (0,2,0) in parent space.
/// World position of child leaf should be (1,2,0).
#[test]
fn flatten_world_transforms() {
    let brep = Arc::new(make_box());

    let child = AssemblyNode::new_leaf_with_transform(
        0,
        "child",
        Arc::clone(&brep),
        DAffine3::from_translation(DVec3::new(0.0, 2.0, 0.0)),
    );
    let parent_xform = DAffine3::from_translation(DVec3::new(1.0, 0.0, 0.0));
    let parent = AssemblyNode {
        id: 0,
        name: "parent".to_string(),
        transform: parent_xform,
        content: NodeContent::Assembly(vec![child]),
        metadata: Default::default(),
    };

    let mut asm = Assembly::new("test");
    asm.add_node(parent);

    let flat = asm.flatten();
    assert_eq!(flat.len(), 1, "one leaf");

    let (_b, world_xform) = &flat[0];
    let world_origin = world_xform.translation;
    assert!(
        (world_origin - DVec3::new(1.0, 2.0, 0.0)).length() < 1e-9,
        "expected world origin (1,2,0), got {:?}",
        world_origin
    );
}

/// Two boxes (6 faces each) → `to_brep()` should give 12 faces total.
#[test]
fn to_brep_face_count() {
    let box_brep = Arc::new(make_box());
    let mut asm = Assembly::new("two_boxes");
    asm.add_part("a", Arc::clone(&box_brep));
    asm.add_part_at(
        "b",
        Arc::clone(&box_brep),
        DVec3::new(3.0, 0.0, 0.0),
    );

    let merged = asm.to_brep();
    assert_eq!(
        face_count(&merged),
        12,
        "two boxes should have 12 faces total"
    );
}

/// The same Arc<BRep> used in two nodes: flatten should give 2 entries,
/// and Arc::ptr_eq confirms they share the same underlying data.
#[test]
fn instance_shared_arc() {
    let shared = Arc::new(make_box());
    let mut asm = Assembly::new("shared");
    asm.add_part("inst1", Arc::clone(&shared));
    asm.add_part_with_transform(
        "inst2",
        Arc::clone(&shared),
        DAffine3::from_translation(DVec3::new(5.0, 0.0, 0.0)),
    );

    let flat = asm.flatten();
    assert_eq!(flat.len(), 2);
    // Both entries point to the same underlying allocation.
    assert!(Arc::ptr_eq(&flat[0].0, &flat[1].0));
}

/// Serialize an Assembly to JSON and deserialize it back; verify structure integrity.
#[test]
fn serde_roundtrip() {
    let box_brep = Arc::new(make_box());
    let mut asm = Assembly::new("serde_test");
    asm.add_part("part_a", Arc::clone(&box_brep));
    asm.add_part_at("part_b", Arc::clone(&box_brep), DVec3::new(2.0, 0.0, 0.0));

    let json = serde_json::to_string(&asm).expect("serialize");
    let restored: Assembly = serde_json::from_str(&json).expect("deserialize");

    assert_eq!(restored.name, "serde_test");
    assert_eq!(restored.roots.len(), 2);
    assert_eq!(restored.roots[0].name, "part_a");
    assert_eq!(restored.roots[1].name, "part_b");

    // Verify instance count is preserved.
    assert_eq!(restored.instance_count(), 2);
}

/// `assembly_from_parts` convenience constructor: build from flat list.
#[test]
fn assembly_from_parts_helper() {
    let box_brep = make_box();
    let sphere_brep =
        make_sphere_brep(DVec3::ZERO, 1.0).unwrap();

    let parts = vec![
        ("box".to_string(), box_brep, DAffine3::IDENTITY),
        (
            "sphere".to_string(),
            sphere_brep,
            DAffine3::from_translation(DVec3::new(3.0, 0.0, 0.0)),
        ),
    ];

    let asm = assembly_from_parts("mixed", parts);
    assert_eq!(asm.instance_count(), 2);

    let flat = asm.flatten();
    // Sphere instance should be at x=3.
    let sphere_trans = flat[1].1.translation;
    assert!(
        (sphere_trans.x - 3.0).abs() < 1e-9,
        "sphere x should be 3, got {}",
        sphere_trans.x
    );
}
