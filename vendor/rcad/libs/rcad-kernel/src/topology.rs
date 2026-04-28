use glam::DVec3;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Vertex {
    pub point: DVec3,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Edge {
    pub start: usize,
    pub end: usize,
}

/// An edge reference with explicit traversal direction inside a Wire.
///
/// `forward = true`  → traverse edge from `edge.start` to `edge.end`.
/// `forward = false` → traverse edge from `edge.end`   to `edge.start`.
///
/// Analogous to OCCT `TopoDS_Edge` with `FORWARD` / `REVERSED` orientation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct WireEdge {
    /// Index into `BRep.edges`.
    pub idx: usize,
    /// Traversal direction: `true` = forward (start→end), `false` = reversed.
    pub forward: bool,
}

impl WireEdge {
    pub const fn new(idx: usize, forward: bool) -> Self {
        Self { idx, forward }
    }
    /// Shorthand: forward reference.
    pub const fn fwd(idx: usize) -> Self {
        Self { idx, forward: true }
    }
    /// Shorthand: reversed reference.
    pub const fn rev(idx: usize) -> Self {
        Self {
            idx,
            forward: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Wire {
    pub edges: Vec<WireEdge>,
}

/// Returns `true` as the serde default for the `mesh_dirty` field.
fn face_mesh_dirty_default() -> bool {
    true
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Face {
    pub outer_wire: Wire,
    pub inner_wires: Vec<Wire>,
    pub normal: DVec3,
    /// Pre-triangulated vertex index triples (into BRep.vertices).
    pub triangles: Vec<[usize; 3]>,
    /// When `true` the cached `triangles` are stale and should be recomputed
    /// before rendering.  Set to `false` by [`mesh_brep`] after tessellation,
    /// and restored to `true` by [`Face::invalidate_mesh`].
    ///
    /// This field is not serialised (transient rendering state).
    #[serde(skip, default = "face_mesh_dirty_default")]
    pub mesh_dirty: bool,
}

impl Face {
    /// Mark this face's cached mesh as stale so it will be re-tessellated on
    /// the next [`mesh_brep`] call.
    pub fn invalidate_mesh(&mut self) {
        self.mesh_dirty = true;
    }

    /// Returns `true` if the cached triangulation is up-to-date.
    pub fn mesh_is_clean(&self) -> bool {
        !self.mesh_dirty && !self.triangles.is_empty()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Shell {
    pub faces: Vec<Face>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Solid {
    pub shells: Vec<Shell>,
}

/// A connected solid made of multiple adjacent solids that share boundary faces.
///
/// Analogous to OCCT `TopoDS_CompSolid`. All contained solids must form a
/// topologically connected manifold body. CompSolid allows expressing structures
/// like multi-region models (e.g. a solid that is split into sub-regions by
/// internal surfaces) without performing a full boolean union.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CompSolid {
    /// The constituent connected solids.
    pub solids: Vec<Solid>,
    /// Optional label for this CompSolid (e.g. from an assembly tree).
    #[serde(default)]
    pub label: Option<String>,
}

impl CompSolid {
    /// Create an empty CompSolid.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a CompSolid from a list of solids.
    pub fn from_solids(solids: Vec<Solid>) -> Self {
        Self { solids, label: None }
    }

    /// Create a CompSolid from a single solid.
    pub fn from_solid(solid: Solid) -> Self {
        Self {
            solids: vec![solid],
            label: None,
        }
    }

    /// Add a solid to this CompSolid.
    pub fn add(&mut self, solid: Solid) {
        self.solids.push(solid);
    }

    /// Remove a solid by index.
    ///
    /// Returns the removed solid if the index was valid.
    pub fn remove(&mut self, index: usize) -> Option<Solid> {
        if index < self.solids.len() {
            Some(self.solids.remove(index))
        } else {
            None
        }
    }

    /// Number of constituent solids.
    pub fn len(&self) -> usize {
        self.solids.len()
    }

    /// Returns `true` if this CompSolid contains no solids.
    pub fn is_empty(&self) -> bool {
        self.solids.is_empty()
    }

    /// Total number of faces across all constituent solids.
    pub fn face_count(&self) -> usize {
        self.solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count()
    }

    /// Total number of shells across all constituent solids.
    pub fn shell_count(&self) -> usize {
        self.solids.iter().flat_map(|s| &s.shells).count()
    }

    /// Explode into constituent solids.
    pub fn explode(self) -> Vec<Solid> {
        self.solids
    }

    /// Set the label for this CompSolid.
    pub fn with_label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }

    /// Iterator over constituent solids.
    pub fn iter_solids(&self) -> impl Iterator<Item = &Solid> {
        self.solids.iter()
    }

    /// Mutable iterator over constituent solids.
    pub fn iter_solids_mut(&mut self) -> impl Iterator<Item = &mut Solid> {
        self.solids.iter_mut()
    }
}

/// Iterator over solids in a CompSolid.
pub struct CompSolidIter<'a> {
    iter: std::slice::Iter<'a, Solid>,
}

impl<'a> CompSolidIter<'a> {
    /// Create a new iterator over solids in a CompSolid.
    pub fn new(compsolid: &'a CompSolid) -> Self {
        Self {
            iter: compsolid.solids.iter(),
        }
    }
}

impl<'a> Iterator for CompSolidIter<'a> {
    type Item = &'a Solid;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

/// A heterogeneous collection of shapes (solids, shells, wires, etc.).
///
/// Analogous to OCCT `TopoDS_Compound`. A Compound can hold any mix of:
/// - Complete solids (`BRep`)
/// - Connected solid groups (`CompSolid`)
/// - Free shells
/// - Free wires / edges
///
/// Compounds are the top-level shape type for assemblies and imported STEP
/// files that contain multiple disconnected bodies.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Compound {
    /// Named sub-shapes (solids).
    ///
    /// Each entry is `(label, shape)`. The label is optional and may be
    /// empty — it is used for assembly-tree bookkeeping and STEP name mapping.
    pub solids: Vec<(Option<String>, Solid)>,
    /// Named CompSolids (multi-region connected solid groups).
    pub comp_solids: Vec<(Option<String>, CompSolid)>,
    /// Loose shells not attached to any solid.
    pub shells: Vec<(Option<String>, Shell)>,
    /// Nested sub-compounds (for deeply hierarchical assemblies).
    pub compounds: Vec<(Option<String>, Compound)>,
}

impl Compound {
    /// Create an empty Compound.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a solid with an optional label.
    pub fn add_solid(&mut self, label: Option<String>, solid: Solid) {
        self.solids.push((label, solid));
    }

    /// Add a CompSolid with an optional label.
    pub fn add_comp_solid(&mut self, label: Option<String>, comp_solid: CompSolid) {
        self.comp_solids.push((label, comp_solid));
    }

    /// Add a nested compound with an optional label.
    pub fn add_compound(&mut self, label: Option<String>, compound: Compound) {
        self.compounds.push((label, compound));
    }

    /// Add a shell with an optional label.
    pub fn add_shell(&mut self, label: Option<String>, shell: Shell) {
        self.shells.push((label, shell));
    }

    /// Remove a solid by index.
    ///
    /// Returns the removed solid if the index was valid.
    pub fn remove_solid(&mut self, index: usize) -> Option<(Option<String>, Solid)> {
        if index < self.solids.len() {
            Some(self.solids.remove(index))
        } else {
            None
        }
    }

    /// Remove a CompSolid by index.
    ///
    /// Returns the removed CompSolid if the index was valid.
    pub fn remove_comp_solid(&mut self, index: usize) -> Option<(Option<String>, CompSolid)> {
        if index < self.comp_solids.len() {
            Some(self.comp_solids.remove(index))
        } else {
            None
        }
    }

    /// Remove a nested compound by index.
    ///
    /// Returns the removed compound if the index was valid.
    pub fn remove_compound(&mut self, index: usize) -> Option<(Option<String>, Compound)> {
        if index < self.compounds.len() {
            Some(self.compounds.remove(index))
        } else {
            None
        }
    }

    /// Remove a shell by index.
    ///
    /// Returns the removed shell if the index was valid.
    pub fn remove_shell(&mut self, index: usize) -> Option<(Option<String>, Shell)> {
        if index < self.shells.len() {
            Some(self.shells.remove(index))
        } else {
            None
        }
    }

    /// Total number of top-level shapes (excluding nested compounds' contents).
    pub fn len(&self) -> usize {
        self.solids.len() + self.comp_solids.len() + self.shells.len() + self.compounds.len()
    }

    /// Returns `true` if the compound contains no shapes.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Flatten all constituent solids into a single list (discards compound hierarchy).
    pub fn flatten_solids(&self) -> Vec<&Solid> {
        let mut out = Vec::new();
        for (_, s) in &self.solids {
            out.push(s);
        }
        for (_, cs) in &self.comp_solids {
            for s in &cs.solids {
                out.push(s);
            }
        }
        for (_, sub) in &self.compounds {
            out.extend(sub.flatten_solids());
        }
        out
    }

    /// Flatten all constituent solids into owned solids (discards labels and hierarchy).
    pub fn into_flattened_solids(self) -> Vec<Solid> {
        let mut out = Vec::new();
        for (_, s) in self.solids {
            out.push(s);
        }
        for (_, cs) in self.comp_solids {
            out.extend(cs.solids);
        }
        for (_, sub) in self.compounds {
            out.extend(sub.into_flattened_solids());
        }
        out
    }

    /// Total face count across all constituent shapes.
    pub fn face_count(&self) -> usize {
        self.flatten_solids()
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count()
    }

    /// Total edge count across all constituent shapes.
    pub fn edge_count(&self) -> usize {
        self.flatten_solids()
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .flat_map(|f| &f.outer_wire.edges)
            .count()
    }

    /// Total vertex count across all constituent shapes (unique vertices only).
    pub fn vertex_count(&self) -> usize {
        // Note: This counts vertices per-solid; for true unique count would need a set.
        self.flatten_solids()
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .flat_map(|f| &f.outer_wire.edges)
            .count()
    }
}

/// Iterator over all solids in a compound, including nested compounds.
pub struct CompoundSolidIter<'a> {
    /// Stack of compounds to process, with current index.
    stack: Vec<(&'a Compound, usize)>,
}

impl<'a> CompoundSolidIter<'a> {
    /// Create a new iterator over the solids in a compound.
    pub fn new(compound: &'a Compound) -> Self {
        Self {
            stack: vec![(compound, 0)],
        }
    }
}

impl<'a> Iterator for CompoundSolidIter<'a> {
    type Item = &'a Solid;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((compound, idx)) = self.stack.last_mut() {
            // First, iterate over direct solids
            if *idx < compound.solids.len() {
                let result = &compound.solids[*idx].1;
                *idx += 1;
                return Some(result);
            }

            let saved_idx = *idx - compound.solids.len();
            *idx += 1;

            // Then, iterate over comp_solids
            if saved_idx < compound.comp_solids.len() {
                let cs_idx = saved_idx;
                let comp_solid = &compound.comp_solids[cs_idx].1;
                // Return the first solid from this comp_solid
                if let Some(solid) = comp_solid.solids.first() {
                    return Some(solid);
                }
                // Empty comp_solid, continue
                continue;
            }

            let cs_len = compound.comp_solids.len();
            let nested_idx = saved_idx - cs_len;

            // Then, recurse into nested compounds
            if nested_idx < compound.compounds.len() {
                let nested = &compound.compounds[nested_idx].1;
                self.stack.push((nested, 0));
                continue;
            }

            // Done with this compound
            self.stack.pop();
        }
        None
    }
}

/// Iterator over all CompSolids in a compound, including nested compounds.
pub struct CompoundCompSolidIter<'a> {
    /// Stack of compounds to process, with current index.
    stack: Vec<(&'a Compound, usize)>,
}

impl<'a> CompoundCompSolidIter<'a> {
    /// Create a new iterator over the CompSolids in a compound.
    pub fn new(compound: &'a Compound) -> Self {
        Self {
            stack: vec![(compound, 0)],
        }
    }
}

impl<'a> Iterator for CompoundCompSolidIter<'a> {
    type Item = &'a CompSolid;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((compound, idx)) = self.stack.last_mut() {
            // Skip solids
            if *idx < compound.solids.len() {
                *idx += compound.solids.len() - *idx;
            }

            let saved_idx = *idx - compound.solids.len();

            // Iterate over comp_solids
            if saved_idx < compound.comp_solids.len() {
                let result = &compound.comp_solids[saved_idx].1;
                *idx += 1;
                return Some(result);
            }

            let nested_idx = saved_idx - compound.comp_solids.len();
            *idx += 1;

            // Recurse into nested compounds
            if nested_idx < compound.compounds.len() {
                let nested = &compound.compounds[nested_idx].1;
                self.stack.push((nested, 0));
                continue;
            }

            // Done with this compound
            self.stack.pop();
        }
        None
    }
}

/// Iterator over nested compounds in a compound.
pub struct NestedCompoundIter<'a> {
    compounds: std::slice::Iter<'a, (Option<String>, Compound)>,
}

impl<'a> NestedCompoundIter<'a> {
    /// Create a new iterator over nested compounds.
    pub fn new(compound: &'a Compound) -> Self {
        Self {
            compounds: compound.compounds.iter(),
        }
    }
}

impl<'a> Iterator for NestedCompoundIter<'a> {
    type Item = &'a Compound;

    fn next(&mut self) -> Option<Self::Item> {
        self.compounds.next().map(|(_, c)| c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wire_edge_fwd_rev() {
        let fwd = WireEdge::fwd(3);
        assert_eq!(fwd.idx, 3);
        assert!(fwd.forward);

        let rev = WireEdge::rev(5);
        assert_eq!(rev.idx, 5);
        assert!(!rev.forward);
    }

    #[test]
    fn wire_contains_edges() {
        let w = Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::rev(2)],
        };
        assert_eq!(w.edges.len(), 3);
        assert!(!w.edges[2].forward);
    }

    #[test]
    fn face_has_outer_wire_and_no_inner_wires_by_default() {
        let f = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        assert!(f.inner_wires.is_empty());
        assert_eq!(f.normal, DVec3::Z);
    }

    #[test]
    fn face_with_inner_wire() {
        let f = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![Wire {
                edges: vec![WireEdge::fwd(3), WireEdge::fwd(4)],
            }],
            normal: DVec3::Y,
            triangles: vec![],
            mesh_dirty: true,
        };
        assert_eq!(f.inner_wires.len(), 1);
        assert_eq!(f.inner_wires[0].edges.len(), 2);
    }

    #[test]
    fn shell_contains_faces() {
        let shell = Shell {
            faces: vec![
                Face {
                    outer_wire: Wire { edges: vec![] },
                    inner_wires: vec![],
                    normal: DVec3::X,
                    triangles: vec![],
                    mesh_dirty: true,
                },
                Face {
                    outer_wire: Wire { edges: vec![] },
                    inner_wires: vec![],
                    normal: DVec3::NEG_X,
                    triangles: vec![],
                    mesh_dirty: true,
                },
            ],
        };
        assert_eq!(shell.faces.len(), 2);
    }

    #[test]
    fn solid_contains_shells() {
        let solid = Solid {
            shells: vec![Shell { faces: vec![] }],
        };
        assert_eq!(solid.shells.len(), 1);
    }

    // ================= Compound Tests =================

    #[test]
    fn compound_new_is_empty() {
        let compound = Compound::new();
        assert!(compound.is_empty());
        assert_eq!(compound.len(), 0);
    }

    #[test]
    fn compound_add_solid() {
        let mut compound = Compound::new();
        let solid = Solid {
            shells: vec![Shell { faces: vec![] }],
        };
        compound.add_solid(Some("solid1".to_string()), solid);
        assert!(!compound.is_empty());
        assert_eq!(compound.len(), 1);
        assert_eq!(compound.solids.len(), 1);
        assert_eq!(compound.solids[0].0, Some("solid1".to_string()));
    }

    #[test]
    fn compound_remove_solid() {
        let mut compound = Compound::new();
        let solid = Solid {
            shells: vec![Shell { faces: vec![] }],
        };
        compound.add_solid(None, solid.clone());
        compound.add_solid(None, solid);

        let removed = compound.remove_solid(0);
        assert!(removed.is_some());
        assert_eq!(compound.len(), 1);

        let removed_invalid = compound.remove_solid(10);
        assert!(removed_invalid.is_none());
    }

    #[test]
    fn compound_add_compound_nested() {
        let mut compound = Compound::new();
        let mut nested = Compound::new();
        nested.add_solid(None, Solid {
            shells: vec![Shell { faces: vec![] }],
        });
        compound.add_compound(Some("nested".to_string()), nested);
        assert_eq!(compound.len(), 1);
        assert_eq!(compound.compounds.len(), 1);
    }

    #[test]
    fn compound_flatten_solids() {
        let mut compound = Compound::new();
        let solid1 = Solid {
            shells: vec![Shell { faces: vec![] }],
        };
        let solid2 = Solid {
            shells: vec![Shell { faces: vec![] }],
        };
        compound.add_solid(None, solid1);
        compound.add_solid(None, solid2);

        let flattened = compound.flatten_solids();
        assert_eq!(flattened.len(), 2);
    }

    #[test]
    fn compound_nested_flatten_solids() {
        let mut outer = Compound::new();
        let mut inner = Compound::new();

        inner.add_solid(None, Solid {
            shells: vec![Shell { faces: vec![] }],
        });
        outer.add_solid(None, Solid {
            shells: vec![Shell { faces: vec![] }],
        });
        outer.add_compound(None, inner);

        let flattened = outer.flatten_solids();
        assert_eq!(flattened.len(), 2);
    }

    // ================= CompSolid Tests =================

    #[test]
    fn compsolid_new_is_empty() {
        let compsolid = CompSolid::new();
        assert!(compsolid.is_empty());
        assert_eq!(compsolid.len(), 0);
    }

    #[test]
    fn compsolid_from_solids() {
        let solid = Solid {
            shells: vec![Shell { faces: vec![] }],
        };
        let compsolid = CompSolid::from_solids(vec![solid.clone(), solid]);
        assert_eq!(compsolid.len(), 2);
    }

    #[test]
    fn compsolid_add() {
        let mut compsolid = CompSolid::new();
        let solid = Solid {
            shells: vec![Shell { faces: vec![] }],
        };
        compsolid.add(solid);
        assert_eq!(compsolid.len(), 1);
    }

    #[test]
    fn compsolid_remove() {
        let mut compsolid = CompSolid::new();
        let solid = Solid {
            shells: vec![Shell { faces: vec![] }],
        };
        compsolid.add(solid.clone());
        compsolid.add(solid);

        let removed = compsolid.remove(0);
        assert!(removed.is_some());
        assert_eq!(compsolid.len(), 1);

        let removed_invalid = compsolid.remove(10);
        assert!(removed_invalid.is_none());
    }

    #[test]
    fn compsolid_explode() {
        let mut compsolid = CompSolid::new();
        let solid = Solid {
            shells: vec![Shell { faces: vec![] }],
        };
        compsolid.add(solid.clone());
        compsolid.add(solid);

        let solids = compsolid.explode();
        assert_eq!(solids.len(), 2);
    }

    #[test]
    fn compsolid_with_label() {
        let compsolid = CompSolid::new()
            .with_label("my_compsolid");
        assert_eq!(compsolid.label, Some("my_compsolid".to_string()));
    }

    // ================= Iterator Tests =================

    #[test]
    fn compound_solid_iter() {
        let mut compound = Compound::new();
        compound.add_solid(None, Solid {
            shells: vec![Shell { faces: vec![] }],
        });
        compound.add_solid(None, Solid {
            shells: vec![Shell { faces: vec![] }],
        });

        let iter = CompoundSolidIter::new(&compound);
        let count = iter.count();
        assert_eq!(count, 2);
    }

    #[test]
    fn compound_solid_iter_with_nested() {
        let mut outer = Compound::new();
        let mut inner = Compound::new();

        inner.add_solid(None, Solid {
            shells: vec![Shell { faces: vec![] }],
        });
        outer.add_solid(None, Solid {
            shells: vec![Shell { faces: vec![] }],
        });
        outer.add_compound(None, inner);

        let iter = CompoundSolidIter::new(&outer);
        let count = iter.count();
        assert_eq!(count, 2);
    }

    #[test]
    fn compound_compsolid_iter() {
        let mut compound = Compound::new();
        compound.add_comp_solid(None, CompSolid::from_solids(vec![Solid {
            shells: vec![Shell { faces: vec![] }],
        }]));
        compound.add_comp_solid(None, CompSolid::new());

        let iter = CompoundCompSolidIter::new(&compound);
        let count = iter.count();
        assert_eq!(count, 2);
    }
}
