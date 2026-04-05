use serde::{Deserialize, Serialize};

use crate::solution_type::SolutionType;

/// Product edition controlling feature availability.
///
/// - **Basic**: Limited HFSS (Driven Modal only) + Q3D (Capacitance only), single setup.
/// - **Professional**: Full HFSS + Q3D capabilities, unlimited setups.
/// - **Enterprise**: Professional + Optimetrics + Distributed Solving.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Edition {
    Basic,
    #[default]
    Professional,
    Enterprise,
}

impl Edition {
    /// Whether this edition permits the given solution type.
    pub fn allows_solution_type(&self, st: SolutionType) -> bool {
        match self {
            Edition::Basic => matches!(st, SolutionType::DrivenModal | SolutionType::Q3D_C),
            Edition::Professional | Edition::Enterprise => true,
        }
    }

    /// Whether "Solve All" (analyze all setups) is available.
    pub fn allows_solve_all(&self) -> bool {
        !matches!(self, Edition::Basic)
    }

    /// Whether the Optimetrics module (parametric sweeps, optimization, etc.) is available.
    pub fn allows_optimetrics(&self) -> bool {
        matches!(self, Edition::Enterprise)
    }

    /// Whether distributed (multi-machine) solving is available.
    pub fn allows_distributed_solve(&self) -> bool {
        matches!(self, Edition::Enterprise)
    }

    /// Maximum number of analysis setups allowed. `None` means unlimited.
    pub fn max_setups(&self) -> Option<usize> {
        match self {
            Edition::Basic => Some(1),
            _ => None,
        }
    }

    /// Human-readable name for display in UI.
    pub fn display_name(&self) -> &'static str {
        match self {
            Edition::Basic => "Basic",
            Edition::Professional => "Professional",
            Edition::Enterprise => "Enterprise",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_allows_driven_modal_and_q3d_c() {
        let e = Edition::Basic;
        assert!(e.allows_solution_type(SolutionType::DrivenModal));
        assert!(e.allows_solution_type(SolutionType::Q3D_C));
    }

    #[test]
    fn basic_blocks_other_solution_types() {
        let e = Edition::Basic;
        assert!(!e.allows_solution_type(SolutionType::DrivenTerminal));
        assert!(!e.allows_solution_type(SolutionType::Eigenmode));
        assert!(!e.allows_solution_type(SolutionType::Transient));
        assert!(!e.allows_solution_type(SolutionType::SBRPlus));
        assert!(!e.allows_solution_type(SolutionType::Q3D_DCRL));
        assert!(!e.allows_solution_type(SolutionType::Q3D_ACRL));
        assert!(!e.allows_solution_type(SolutionType::Q3D_CG));
    }

    #[test]
    fn professional_allows_all_solution_types() {
        let e = Edition::Professional;
        assert!(e.allows_solution_type(SolutionType::DrivenTerminal));
        assert!(e.allows_solution_type(SolutionType::Eigenmode));
        assert!(e.allows_solution_type(SolutionType::Q3D_DCRL));
        assert!(e.allows_solution_type(SolutionType::Q3D_CG));
    }

    #[test]
    fn basic_no_solve_all() {
        assert!(!Edition::Basic.allows_solve_all());
        assert!(Edition::Professional.allows_solve_all());
        assert!(Edition::Enterprise.allows_solve_all());
    }

    #[test]
    fn optimetrics_enterprise_only() {
        assert!(!Edition::Basic.allows_optimetrics());
        assert!(!Edition::Professional.allows_optimetrics());
        assert!(Edition::Enterprise.allows_optimetrics());
    }

    #[test]
    fn basic_max_one_setup() {
        assert_eq!(Edition::Basic.max_setups(), Some(1));
        assert_eq!(Edition::Professional.max_setups(), None);
        assert_eq!(Edition::Enterprise.max_setups(), None);
    }
}
