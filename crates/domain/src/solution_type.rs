use serde::{Deserialize, Serialize};

/// Electromagnetic simulation solution type.
///
/// EMStudio supports two solver families:
/// - **HFSS**: Full-wave FEM (Finite Element Method)
/// - **Q3D**: Quasi-static parasitic extraction (MoM + FMM)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum SolutionType {
    // HFSS Full-wave
    /// Driven modal analysis — S-parameters via mode decomposition.
    DrivenModal,
    /// Driven terminal analysis — S-parameters via terminal V/I.
    DrivenTerminal,
    /// Eigenmode analysis — resonant frequencies and Q factors.
    Eigenmode,
    /// Transient analysis — time-domain FEM.
    Transient,
    /// SBR+ — ray tracing + physical optics for large scatterers.
    SBRPlus,

    // Q3D Quasi-static
    /// DC resistance + low-frequency inductance.
    Q3D_DCRL,
    /// AC resistance + inductance (with skin/proximity effects).
    Q3D_ACRL,
    /// Capacitance matrix (electrostatics).
    Q3D_C,
    /// Capacitance + conductance matrix (with dielectric loss).
    Q3D_CG,
}

impl SolutionType {
    /// Returns true if this is an HFSS (full-wave) solution type.
    pub fn is_hfss(&self) -> bool {
        matches!(
            self,
            Self::DrivenModal
                | Self::DrivenTerminal
                | Self::Eigenmode
                | Self::Transient
                | Self::SBRPlus
        )
    }

    /// Returns true if this is a Q3D (quasi-static) solution type.
    pub fn is_q3d(&self) -> bool {
        matches!(self, Self::Q3D_DCRL | Self::Q3D_ACRL | Self::Q3D_C | Self::Q3D_CG)
    }
}
