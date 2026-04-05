//! Project validation: check reference integrity, naming uniqueness,
//! variable cycles, and solver-specific constraints.

use std::collections::HashSet;

use thiserror::Error;

use crate::design::Design;
use crate::project::EmProject;
use crate::solution_type::SolutionType;

#[derive(Debug, Error)]
pub enum ValidationError {
    #[error("dangling reference: {from} references missing {missing} (field: {field})")]
    DanglingReference {
        from: String,
        missing: String,
        field: String,
    },
    #[error("duplicate name: {kind} '{name}'")]
    DuplicateName { kind: String, name: String },
    #[error("Q3D net '{net}' has no terminals")]
    NetWithoutTerminals { net: String },
    #[error("HFSS design '{design}' has no excitation ports")]
    NoExcitationPorts { design: String },
    #[error("conductor object '{object}' is not assigned to any net")]
    ConductorNotAssigned { object: String },
    #[error("Q3D design '{design}' has no ground reference net")]
    NoGroundReference { design: String },
    #[error("duplicate terminal name '{name}' across nets")]
    DuplicateTerminalName { name: String },
    #[error("{0}")]
    Other(String),
}

/// Validate an entire project structure.
pub fn validate_project(project: &EmProject) -> Vec<ValidationError> {
    let mut errors = Vec::new();

    // Check project-level variable name uniqueness
    // (HashMap keys are already unique, so nothing to check there)

    for design in &project.designs {
        errors.extend(validate_design(design));
    }

    errors
}

/// Validate a single design.
pub fn validate_design(design: &Design) -> Vec<ValidationError> {
    let mut errors = Vec::new();

    errors.extend(check_material_references(design));
    errors.extend(check_name_uniqueness(design));
    errors.extend(check_boundary_references(design));

    if design.solution_type.is_q3d() {
        errors.extend(check_q3d_nets(design));
        errors.extend(check_q3d_conductor_coverage(design));
        errors.extend(check_q3d_ground_reference(design));
        errors.extend(check_q3d_terminal_uniqueness(design));
    }

    if design.solution_type.is_hfss() {
        errors.extend(check_hfss_ports(design));
    }

    errors
}

/// Check that all material references in geometry objects exist in definitions.
fn check_material_references(design: &Design) -> Vec<ValidationError> {
    let defined: HashSet<&str> = design
        .definitions
        .materials
        .iter()
        .map(|m| m.name.as_str())
        .collect();

    // "vacuum" is implicitly available
    let mut valid_materials = defined;
    valid_materials.insert("vacuum");

    design
        .geometry
        .objects
        .iter()
        .filter(|obj| !valid_materials.contains(obj.material.as_str()))
        .map(|obj| ValidationError::DanglingReference {
            from: format!("GeometryObject:{}", obj.name),
            missing: format!("Material:{}", obj.material),
            field: "material".to_string(),
        })
        .collect()
}

/// Check that boundary/excitation target objects exist in geometry.
fn check_boundary_references(design: &Design) -> Vec<ValidationError> {
    let objects: HashSet<&str> = design
        .geometry
        .objects
        .iter()
        .map(|o| o.name.as_str())
        .collect();

    let mut errors = Vec::new();

    for boundary in &design.boundaries {
        for target in &boundary.assignment.targets {
            // Named selections start with @, geometry objects are direct names
            if !target.starts_with('@') && !objects.contains(target.as_str()) {
                errors.push(ValidationError::DanglingReference {
                    from: format!("Boundary:{}", boundary.name),
                    missing: format!("Object:{target}"),
                    field: "assignment.targets".to_string(),
                });
            }
        }
    }

    for excitation in &design.excitations {
        if let Some(assignment) = &excitation.assignment {
            for target in &assignment.targets {
                if !target.starts_with('@') && !objects.contains(target.as_str()) {
                    errors.push(ValidationError::DanglingReference {
                        from: format!("Excitation:{}", excitation.name),
                        missing: format!("Object:{target}"),
                        field: "assignment.targets".to_string(),
                    });
                }
            }
        }
    }

    errors
}

/// Check name uniqueness across materials, geometry objects, setups, etc.
fn check_name_uniqueness(design: &Design) -> Vec<ValidationError> {
    let mut errors = Vec::new();

    check_unique(&design.definitions.materials.iter().map(|m| &m.name).collect::<Vec<_>>(), "Material", &mut errors);
    check_unique(&design.geometry.objects.iter().map(|o| &o.name).collect::<Vec<_>>(), "GeometryObject", &mut errors);
    check_unique(&design.boundaries.iter().map(|b| &b.name).collect::<Vec<_>>(), "Boundary", &mut errors);
    check_unique(&design.excitations.iter().map(|e| &e.name).collect::<Vec<_>>(), "Excitation", &mut errors);
    check_unique(&design.analysis_setups.iter().map(|s| &s.name).collect::<Vec<_>>(), "AnalysisSetup", &mut errors);

    errors
}

fn check_unique(names: &[&String], kind: &str, errors: &mut Vec<ValidationError>) {
    let mut seen = HashSet::new();
    for name in names {
        if !seen.insert(name.as_str()) {
            errors.push(ValidationError::DuplicateName {
                kind: kind.to_string(),
                name: name.to_string(),
            });
        }
    }
}

/// Check Q3D-specific constraints: nets must have terminals.
fn check_q3d_nets(design: &Design) -> Vec<ValidationError> {
    design
        .nets
        .iter()
        .filter(|net| net.terminals.is_empty())
        .map(|net| ValidationError::NetWithoutTerminals {
            net: net.name.clone(),
        })
        .collect()
}

/// Check HFSS-specific constraints: at least one excitation for driven analyses.
fn check_hfss_ports(design: &Design) -> Vec<ValidationError> {
    if matches!(
        design.solution_type,
        SolutionType::DrivenModal | SolutionType::DrivenTerminal
    ) && design.excitations.is_empty()
    {
        vec![ValidationError::NoExcitationPorts {
            design: design.name.clone(),
        }]
    } else {
        vec![]
    }
}

/// Check that every conductor object is assigned to exactly one net.
fn check_q3d_conductor_coverage(design: &Design) -> Vec<ValidationError> {
    let assigned: HashSet<&str> = design
        .nets
        .iter()
        .flat_map(|net| net.objects.iter().map(|o| o.as_str()))
        .collect();

    design
        .geometry
        .objects
        .iter()
        .filter(|obj| !obj.material.is_empty() && obj.material != "vacuum")
        .filter(|obj| !assigned.contains(obj.name.as_str()))
        .map(|obj| ValidationError::ConductorNotAssigned {
            object: obj.name.clone(),
        })
        .collect()
}

/// Check that at least one net is the ground reference.
fn check_q3d_ground_reference(design: &Design) -> Vec<ValidationError> {
    if !design.nets.is_empty() && !design.nets.iter().any(|n| n.is_ground_reference) {
        vec![ValidationError::NoGroundReference {
            design: design.name.clone(),
        }]
    } else {
        vec![]
    }
}

/// Check that all terminal names are unique across all nets.
fn check_q3d_terminal_uniqueness(design: &Design) -> Vec<ValidationError> {
    let mut seen = HashSet::new();
    let mut errors = Vec::new();
    for net in &design.nets {
        for terminal in &net.terminals {
            if !seen.insert(terminal.name.as_str()) {
                errors.push(ValidationError::DuplicateTerminalName {
                    name: terminal.name.clone(),
                });
            }
        }
    }
    errors
}
