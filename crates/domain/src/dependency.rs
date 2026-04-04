//! Definition-reference dependency graph.
//!
//! Tracks which definitions (materials, variables, coordinate systems, etc.)
//! are referenced by which consumers (boundaries, excitations, geometry, etc.).

use std::collections::HashSet;

use crate::design::Design;

/// A unique identifier for a definition.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum DefinitionId {
    ProjectVariable(String),
    DesignVariable(String),
    Dataset(String),
    Material(String),
    CoordinateSystem(String),
    NamedSelection(String),
    GeometryObject(String),
}

impl std::fmt::Display for DefinitionId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ProjectVariable(n) => write!(f, "Variable:{n}"),
            Self::DesignVariable(n) => write!(f, "Variable:{n}"),
            Self::Dataset(n) => write!(f, "Dataset:{n}"),
            Self::Material(n) => write!(f, "Material:{n}"),
            Self::CoordinateSystem(n) => write!(f, "CoordSys:{n}"),
            Self::NamedSelection(n) => write!(f, "NamedSelection:{n}"),
            Self::GeometryObject(n) => write!(f, "Object:{n}"),
        }
    }
}

/// A reference from a consumer to a definition.
#[derive(Debug, Clone)]
pub struct Reference {
    /// Source description, e.g. `"Boundary:Radiation1"`.
    pub from: String,
    /// The definition being referenced.
    pub to: DefinitionId,
    /// Which field contains the reference, e.g. `"material"`.
    pub field: String,
}

/// Dependency graph of all definition→reference relationships in a design.
#[derive(Debug, Default)]
pub struct DependencyGraph {
    /// All known definitions.
    pub definitions: HashSet<DefinitionId>,
    /// All references from consumers to definitions.
    pub references: Vec<Reference>,
}

impl DependencyGraph {
    pub fn new() -> Self {
        Self::default()
    }

    /// Build the dependency graph from a design.
    pub fn from_design(design: &Design) -> Self {
        let mut graph = Self::new();

        // Register definitions
        for mat in &design.definitions.materials {
            graph.definitions.insert(DefinitionId::Material(mat.name.clone()));
        }
        for cs in &design.definitions.coordinate_systems {
            graph.definitions.insert(DefinitionId::CoordinateSystem(cs.name.clone()));
        }
        for ns in &design.definitions.named_selections {
            graph.definitions.insert(DefinitionId::NamedSelection(ns.name.clone()));
        }
        for (name, _) in &design.local_variables {
            graph.definitions.insert(DefinitionId::DesignVariable(name.clone()));
        }
        for obj in &design.geometry.objects {
            graph.definitions.insert(DefinitionId::GeometryObject(obj.name.clone()));
        }

        // Scan references from geometry objects → materials
        for obj in &design.geometry.objects {
            graph.references.push(Reference {
                from: format!("GeometryObject:{}", obj.name),
                to: DefinitionId::Material(obj.material.clone()),
                field: "material".to_string(),
            });
        }

        // Scan references from boundaries → geometry objects
        for boundary in &design.boundaries {
            for target in &boundary.assignment.targets {
                graph.references.push(Reference {
                    from: format!("Boundary:{}", boundary.name),
                    to: DefinitionId::GeometryObject(target.clone()),
                    field: "assignment.targets".to_string(),
                });
            }
        }

        // Scan references from excitations → geometry objects
        for excitation in &design.excitations {
            if let Some(assignment) = &excitation.assignment {
                for target in &assignment.targets {
                    graph.references.push(Reference {
                        from: format!("Excitation:{}", excitation.name),
                        to: DefinitionId::GeometryObject(target.clone()),
                        field: "assignment.targets".to_string(),
                    });
                }
            }
        }

        // Scan references from nets → geometry objects (Q3D)
        for net in &design.nets {
            for obj_name in &net.objects {
                graph.references.push(Reference {
                    from: format!("Net:{}", net.name),
                    to: DefinitionId::GeometryObject(obj_name.clone()),
                    field: "objects".to_string(),
                });
            }
        }

        graph
    }

    /// Find all references that depend on a given definition.
    pub fn find_dependents(&self, def_id: &DefinitionId) -> Vec<&Reference> {
        self.references
            .iter()
            .filter(|r| &r.to == def_id)
            .collect()
    }

    /// Check if a definition can be safely deleted (no dependents).
    pub fn can_delete(&self, def_id: &DefinitionId) -> Result<(), Vec<String>> {
        let dependents = self.find_dependents(def_id);
        if dependents.is_empty() {
            Ok(())
        } else {
            Err(dependents.iter().map(|r| r.from.clone()).collect())
        }
    }

    /// Mark analysis setups as stale when a definition changes.
    /// Returns the names of setups that should be marked stale.
    pub fn affected_setups(&self, changed_def: &DefinitionId, design: &Design) -> Vec<String> {
        // For now, if any geometry/material/boundary changes, all setups are stale
        // (a more granular approach would trace the dependency chain)
        let dependents = self.find_dependents(changed_def);
        if dependents.is_empty() {
            return Vec::new();
        }
        design
            .analysis_setups
            .iter()
            .map(|s| s.name.clone())
            .collect()
    }
}

/// Mark setups as stale in the solution index when a definition changes.
pub fn mark_stale(design: &mut Design, changed_def: &DefinitionId) {
    let graph = DependencyGraph::from_design(design);
    let affected = graph.affected_setups(changed_def, design);

    for setup_name in &affected {
        if let Some(status) = design.solution_index.setups.get_mut(setup_name) {
            status.is_stale = true;
        }
    }

    if !affected.is_empty() {
        design.solution_index.stale_reason =
            Some(format!("Modified: {changed_def}"));
    }
}
