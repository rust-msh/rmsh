use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnsysDesignKind {
    Hfss,
    Q3d,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnsysSolutionType {
    // HFSS
    DrivenModal,
    DrivenTerminal,
    Eigenmode,
    Transient,
    SbrPlus,

    // Q3D
    Q3dDcrl,
    Q3dAcrl,
    Q3dC,
    Q3dCg,
    Unknown,
}

impl AnsysSolutionType {
    pub fn as_ansys_label(self) -> &'static str {
        match self {
            Self::DrivenModal => "DrivenModal",
            Self::DrivenTerminal => "DrivenTerminal",
            Self::Eigenmode => "Eigenmode",
            Self::Transient => "Transient",
            Self::SbrPlus => "SBR+",
            Self::Q3dDcrl => "Q3D DC RL",
            Self::Q3dAcrl => "Q3D AC RL",
            Self::Q3dC => "Q3D Capacitance",
            Self::Q3dCg => "Q3D Capacitance + Conductance",
            Self::Unknown => "Unknown",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnsysDesign {
    pub name: String,
    pub kind: AnsysDesignKind,
    pub solution_type: AnsysSolutionType,
    #[serde(default)]
    pub variables: BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnsysProject {
    pub name: String,
    #[serde(default)]
    pub designs: Vec<AnsysDesign>,
}
