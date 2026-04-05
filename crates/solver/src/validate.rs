//! Pre-solve validation: checks that a Design is ready for simulation.

use emstudio_domain::Design;
use emstudio_domain::solution_type::SolutionType;
use emstudio_domain::validation::validate_design;

use crate::error::SolverError;

/// Validate that a Design is ready for solver execution.
///
/// This combines the domain-level validation with solver-specific checks.
pub fn validate_for_solve(design: &Design) -> Result<(), SolverError> {
    // 1. Run domain-level validation (material refs, name uniqueness, etc.).
    let errors = validate_design(design);
    if !errors.is_empty() {
        let messages: Vec<String> = errors.iter().map(|e: &emstudio_domain::validation::ValidationError| e.to_string()).collect();
        return Err(SolverError::Validation(format!(
            "Design validation failed:\n  - {}",
            messages.join("\n  - ")
        )));
    }

    // 2. Check geometry is non-empty.
    if design.geometry.objects.is_empty() {
        return Err(SolverError::Validation(
            "Design has no geometry objects".into(),
        ));
    }

    // 3. Check at least one enabled analysis setup exists.
    let has_enabled_setup = design.analysis_setups.iter().any(|s| s.enabled);
    if !has_enabled_setup {
        return Err(SolverError::Validation(
            "No enabled analysis setup found".into(),
        ));
    }

    // 4. HFSS-specific: driven analyses need at least one excitation.
    if matches!(
        design.solution_type,
        SolutionType::DrivenModal | SolutionType::DrivenTerminal
    ) && design.excitations.is_empty()
    {
        return Err(SolverError::Validation(
            "HFSS driven analysis requires at least one excitation port".into(),
        ));
    }

    // 5. Q3D-specific: needs nets with terminals.
    if design.solution_type.is_q3d() && design.nets.is_empty() {
        return Err(SolverError::Validation(
            "Q3D analysis requires at least one net definition".into(),
        ));
    }

    Ok(())
}
