use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

use crate::BRep;

/// A stable reference to a topological entity in a B-Rep.
///
/// Face indexing follows RCAD's flattened face order (solid/shell/face traversal).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum TopoEntityRef {
    Vertex(usize),
    Edge(usize),
    Face(usize),
    Solid(usize),
}

/// Baseline persistent naming table for topology entities.
///
/// This is a lightweight hook layer analogous to OCCT OCAF naming tables:
/// it provides stable user-level names and bidirectional resolution between
/// names and topology references.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PersistentNamingHooks {
    name_to_ref: BTreeMap<String, TopoEntityRef>,
    ref_to_name: BTreeMap<TopoEntityRef, String>,
}

impl PersistentNamingHooks {
    pub fn new() -> Self {
        Self::default()
    }

    /// Build a baseline naming table with deterministic default labels.
    ///
    /// Generated labels are:
    /// - vertices: `v0`, `v1`, ...
    /// - edges: `e0`, `e1`, ...
    /// - faces: `f0`, `f1`, ... (flattened index)
    /// - solids: `s0`, `s1`, ...
    pub fn with_default_labels_for_brep(brep: &BRep) -> Self {
        let mut out = Self::new();
        for i in 0..brep.vertices.len() {
            out.bind_unchecked(format!("v{i}"), TopoEntityRef::Vertex(i));
        }
        for i in 0..brep.edges.len() {
            out.bind_unchecked(format!("e{i}"), TopoEntityRef::Edge(i));
        }
        for i in 0..flat_face_count(brep) {
            out.bind_unchecked(format!("f{i}"), TopoEntityRef::Face(i));
        }
        for i in 0..brep.solids.len() {
            out.bind_unchecked(format!("s{i}"), TopoEntityRef::Solid(i));
        }
        out
    }

    /// Bind a user-visible `name` to an entity reference.
    ///
    /// If the name or entity is already bound, the old binding is replaced.
    pub fn bind(&mut self, name: impl Into<String>, target: TopoEntityRef) {
        self.bind_unchecked(name.into(), target);
    }

    /// Bind with topology bounds check against `brep`.
    pub fn bind_for_brep(
        &mut self,
        brep: &BRep,
        name: impl Into<String>,
        target: TopoEntityRef,
    ) -> Result<(), String> {
        if !is_valid_ref_for_brep(brep, target) {
            return Err(format!("invalid topology reference for BRep: {target:?}"));
        }
        self.bind_unchecked(name.into(), target);
        Ok(())
    }

    pub fn resolve(&self, name: &str) -> Option<TopoEntityRef> {
        self.name_to_ref.get(name).copied()
    }

    pub fn name_of(&self, target: TopoEntityRef) -> Option<&str> {
        self.ref_to_name.get(&target).map(String::as_str)
    }

    pub fn unbind_name(&mut self, name: &str) -> Option<TopoEntityRef> {
        let target = self.name_to_ref.remove(name)?;
        self.ref_to_name.remove(&target);
        Some(target)
    }

    pub fn unbind_ref(&mut self, target: TopoEntityRef) -> Option<String> {
        let name = self.ref_to_name.remove(&target)?;
        self.name_to_ref.remove(&name);
        Some(name)
    }

    pub fn rename(&mut self, old_name: &str, new_name: impl Into<String>) -> Result<(), String> {
        let Some(target) = self.resolve(old_name) else {
            return Err(format!("name '{old_name}' not found"));
        };
        let new_name = new_name.into();
        if new_name == old_name {
            return Ok(());
        }
        if let Some(existing) = self.resolve(&new_name) {
            return Err(format!(
                "name '{new_name}' is already bound to {existing:?}"
            ));
        }
        self.unbind_name(old_name);
        self.bind_unchecked(new_name, target);
        Ok(())
    }

    /// Returns all invalid bindings for the given `brep`.
    pub fn validate_against_brep(&self, brep: &BRep) -> Vec<String> {
        let mut issues = Vec::new();
        for (name, target) in &self.name_to_ref {
            if !is_valid_ref_for_brep(brep, *target) {
                issues.push(format!("name '{name}' points to out-of-range entity {target:?}"));
            }
        }
        issues
    }

    /// Remove bindings that no longer point to valid topology entities.
    pub fn retain_valid_for_brep(&mut self, brep: &BRep) {
        let invalid_names: Vec<String> = self
            .name_to_ref
            .iter()
            .filter_map(|(name, target)| {
                if is_valid_ref_for_brep(brep, *target) {
                    None
                } else {
                    Some(name.clone())
                }
            })
            .collect();
        for name in invalid_names {
            self.unbind_name(&name);
        }
    }

    pub fn len(&self) -> usize {
        self.name_to_ref.len()
    }

    pub fn is_empty(&self) -> bool {
        self.name_to_ref.is_empty()
    }

    fn bind_unchecked(&mut self, name: String, target: TopoEntityRef) {
        if let Some(old_target) = self.name_to_ref.remove(&name) {
            self.ref_to_name.remove(&old_target);
        }
        if let Some(old_name) = self.ref_to_name.remove(&target) {
            self.name_to_ref.remove(&old_name);
        }
        self.name_to_ref.insert(name.clone(), target);
        self.ref_to_name.insert(target, name);
    }

    /// Propagate names from a pre-operation naming table through an index-level
    /// mapping to produce a post-operation naming table.
    ///
    /// `face_map[old_face_idx]` → `Some(new_face_idx)` if the face survived the
    /// operation and was remapped.  `None` means the face was consumed/deleted.
    /// New entities (generated by the operation) can be named via the
    /// `new_face_names` slice: `(new_face_idx, name)` pairs.
    ///
    /// The same logic applies independently to vertices and edges via their
    /// respective maps.  This method updates `self` in place and returns the set
    /// of names that were dropped (because the entity was removed).
    ///
    /// Analogous to OCCT `TNaming_NamedShape` propagation after a BRep rebuild.
    pub fn propagate_through_remap(
        &mut self,
        face_map: &[Option<usize>],
        edge_map: &[Option<usize>],
        vertex_map: &[Option<usize>],
        new_face_names: &[(usize, String)],
        new_edge_names: &[(usize, String)],
        new_vertex_names: &[(usize, String)],
    ) -> Vec<String> {
        let mut dropped = Vec::new();

        // Collect existing bindings so we can update them.
        let snapshot: Vec<(String, TopoEntityRef)> = self
            .name_to_ref
            .iter()
            .map(|(n, r)| (n.clone(), *r))
            .collect();

        // Clear tables; we rebuild them below.
        self.name_to_ref.clear();
        self.ref_to_name.clear();

        for (name, old_ref) in snapshot {
            let new_ref = match old_ref {
                TopoEntityRef::Face(i) => {
                    if face_map.is_empty() {
                        Some(old_ref)
                    } else {
                        face_map.get(i).and_then(|r| *r).map(TopoEntityRef::Face)
                    }
                }
                TopoEntityRef::Edge(i) => {
                    if edge_map.is_empty() {
                        Some(old_ref)
                    } else {
                        edge_map.get(i).and_then(|r| *r).map(TopoEntityRef::Edge)
                    }
                }
                TopoEntityRef::Vertex(i) => {
                    if vertex_map.is_empty() {
                        Some(old_ref)
                    } else {
                        vertex_map.get(i).and_then(|r| *r).map(TopoEntityRef::Vertex)
                    }
                }
                TopoEntityRef::Solid(_) => Some(old_ref), // solids not remapped
            };
            match new_ref {
                Some(r) => self.bind_unchecked(name, r),
                None => dropped.push(name),
            }
        }

        // Register names for new entities.
        for (idx, name) in new_face_names {
            self.bind_unchecked(name.clone(), TopoEntityRef::Face(*idx));
        }
        for (idx, name) in new_edge_names {
            self.bind_unchecked(name.clone(), TopoEntityRef::Edge(*idx));
        }
        for (idx, name) in new_vertex_names {
            self.bind_unchecked(name.clone(), TopoEntityRef::Vertex(*idx));
        }

        dropped
    }

    /// Convenience: propagate names from an operation that only remaps faces
    /// (e.g. fillet, chamfer).
    ///
    /// `face_map[old_idx]` → `Some(new_idx)` or `None` (removed).
    /// `new_face_names` names any newly generated faces.
    pub fn propagate_face_remap(
        &mut self,
        face_map: &[Option<usize>],
        new_face_names: &[(usize, String)],
    ) -> Vec<String> {
        self.propagate_through_remap(face_map, &[], &[], new_face_names, &[], &[])
    }

    /// Build a simple identity remap for `n` entities (nothing moved or removed).
    pub fn identity_map(n: usize) -> Vec<Option<usize>> {
        (0..n).map(Some).collect()
    }

    /// Iterate over all (name, entity_ref) bindings.
    pub fn iter(&self) -> impl Iterator<Item = (&str, TopoEntityRef)> {
        self.name_to_ref.iter().map(|(n, r)| (n.as_str(), *r))
    }
}

fn flat_face_count(brep: &BRep) -> usize {
    brep
        .solids
        .iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum()
}

fn is_valid_ref_for_brep(brep: &BRep, target: TopoEntityRef) -> bool {
    match target {
        TopoEntityRef::Vertex(i) => i < brep.vertices.len(),
        TopoEntityRef::Edge(i) => i < brep.edges.len(),
        TopoEntityRef::Face(i) => i < flat_face_count(brep),
        TopoEntityRef::Solid(i) => i < brep.solids.len(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BRep, PrimitiveSolid};

    fn unit_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        })
    }

    #[test]
    fn default_labels_cover_basic_topology() {
        let brep = unit_box();
        let hooks = PersistentNamingHooks::with_default_labels_for_brep(&brep);
        assert_eq!(hooks.resolve("v0"), Some(TopoEntityRef::Vertex(0)));
        assert_eq!(hooks.resolve("e0"), Some(TopoEntityRef::Edge(0)));
        assert_eq!(hooks.resolve("f0"), Some(TopoEntityRef::Face(0)));
        assert_eq!(hooks.resolve("s0"), Some(TopoEntityRef::Solid(0)));
    }

    #[test]
    fn bind_and_rename_roundtrip() {
        let brep = unit_box();
        let mut hooks = PersistentNamingHooks::new();
        hooks
            .bind_for_brep(&brep, "mount_hole", TopoEntityRef::Edge(1))
            .expect("bind should succeed");
        assert_eq!(hooks.resolve("mount_hole"), Some(TopoEntityRef::Edge(1)));

        hooks
            .rename("mount_hole", "outer_profile")
            .expect("rename should succeed");
        assert_eq!(hooks.resolve("mount_hole"), None);
        assert_eq!(hooks.resolve("outer_profile"), Some(TopoEntityRef::Edge(1)));
    }

    #[test]
    fn validate_and_retain_invalid_bindings() {
        let brep = unit_box();
        let mut hooks = PersistentNamingHooks::new();
        hooks.bind("bad_edge", TopoEntityRef::Edge(9999));
        hooks.bind("good_vertex", TopoEntityRef::Vertex(0));

        let issues = hooks.validate_against_brep(&brep);
        assert_eq!(issues.len(), 1);
        assert!(issues[0].contains("bad_edge"));

        hooks.retain_valid_for_brep(&brep);
        assert_eq!(hooks.resolve("bad_edge"), None);
        assert_eq!(hooks.resolve("good_vertex"), Some(TopoEntityRef::Vertex(0)));
    }

    // ── propagate_through_remap ───────────────────────────────────────────────

    #[test]
    fn propagate_face_remap_shifts_kept_faces() {
        // Simulate a fillet on a 6-face box: the 2 affected faces get new
        // indices (0→0, 1→1) and 3 new faces are appended at positions 6,7,8.
        let brep = unit_box();
        let mut hooks = PersistentNamingHooks::with_default_labels_for_brep(&brep);

        // Before: f0..f5 are bound.
        assert_eq!(hooks.resolve("f0"), Some(TopoEntityRef::Face(0)));
        assert_eq!(hooks.resolve("f5"), Some(TopoEntityRef::Face(5)));

        // After fillet on edge 0: faces 0 and 1 are trimmed (preserved but remapped
        // to same positions), all others remain. 3 new faces are appended.
        let face_map: Vec<Option<usize>> = vec![Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)];
        let new_faces: Vec<(usize, String)> = vec![
            (6, "fillet_face".to_string()),
            (7, "closing_0".to_string()),
            (8, "closing_1".to_string()),
        ];

        let dropped = hooks.propagate_face_remap(&face_map, &new_faces);
        assert!(dropped.is_empty(), "no faces should be dropped with identity map");

        // Original names still resolve to their (unchanged) positions.
        assert_eq!(hooks.resolve("f0"), Some(TopoEntityRef::Face(0)));
        assert_eq!(hooks.resolve("f5"), Some(TopoEntityRef::Face(5)));
        // New names are bound.
        assert_eq!(hooks.resolve("fillet_face"), Some(TopoEntityRef::Face(6)));
        assert_eq!(hooks.resolve("closing_0"), Some(TopoEntityRef::Face(7)));
    }

    #[test]
    fn propagate_face_remap_drops_removed_faces() {
        let brep = unit_box();
        let mut hooks = PersistentNamingHooks::new();
        hooks.bind("top", TopoEntityRef::Face(5));
        hooks.bind("bottom", TopoEntityRef::Face(0));

        // Simulate removing face 5 (top removed) and keeping face 0 at new index 0.
        let face_map: Vec<Option<usize>> = vec![Some(0), Some(1), Some(2), Some(3), Some(4), None];
        let dropped = hooks.propagate_face_remap(&face_map, &[]);

        assert!(dropped.contains(&"top".to_string()), "top should be in dropped list");
        assert_eq!(hooks.resolve("top"), None);
        assert_eq!(hooks.resolve("bottom"), Some(TopoEntityRef::Face(0)));
    }

    #[test]
    fn propagate_remaps_edge_and_vertex_indices() {
        let brep = unit_box();
        let mut hooks = PersistentNamingHooks::new();
        hooks.bind("edge_a", TopoEntityRef::Edge(3));
        hooks.bind("vert_b", TopoEntityRef::Vertex(5));

        // Simulate a rebuild where edge 3 → 2 and vertex 5 → 4.
        let edge_map: Vec<Option<usize>> = vec![Some(0), Some(1), Some(2), Some(2), Some(4), Some(5),
                                                 Some(6), Some(7), Some(8), Some(9), Some(10), Some(11)];
        let vert_map: Vec<Option<usize>> = vec![Some(0), Some(1), Some(2), Some(3), Some(4), Some(4),
                                                 Some(6), Some(7)];

        let dropped = hooks.propagate_through_remap(&[], &edge_map, &vert_map, &[], &[], &[]);
        assert!(dropped.is_empty());
        assert_eq!(hooks.resolve("edge_a"), Some(TopoEntityRef::Edge(2)));
        assert_eq!(hooks.resolve("vert_b"), Some(TopoEntityRef::Vertex(4)));
    }

    #[test]
    fn identity_map_preserves_all() {
        let n = 6;
        let map = PersistentNamingHooks::identity_map(n);
        assert_eq!(map.len(), n);
        for (i, m) in map.iter().enumerate() {
            assert_eq!(*m, Some(i));
        }
    }

    #[test]
    fn iter_visits_all_bindings() {
        let brep = unit_box();
        let hooks = PersistentNamingHooks::with_default_labels_for_brep(&brep);
        let count = hooks.iter().count();
        // v0..v7 (8) + e0..e11 (12) + f0..f5 (6) + s0 (1) = 27
        assert_eq!(count, 27);
    }
}
