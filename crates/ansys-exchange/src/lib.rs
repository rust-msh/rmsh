pub mod error;
pub mod export;
pub mod import;
pub mod model;

pub use error::AnsysExchangeError;
pub use export::{export_manifest_json_string, export_pyaedt_script};
pub use import::{import_aedt_str, import_manifest_json_str};
#[cfg(not(target_arch = "wasm32"))]
pub use export::export_pyaedt_script_file;
#[cfg(not(target_arch = "wasm32"))]
pub use import::import_aedt_file;
pub use model::{AnsysDesign, AnsysDesignKind, AnsysProject, AnsysSolutionType};

#[cfg(test)]
mod integration_tests {
    use super::*;

    #[test]
    fn import_then_export_script() {
        let input = r#"
oProject = oDesktop.NewProject("RFSystem")
oProject.InsertDesign("HFSS", "Antenna", "DrivenModal", "")
"#;

        let project = import_aedt_str(input).unwrap();
        let script = export_pyaedt_script(&project);

        assert!(script.contains("RFSystem"));
        assert!(script.contains("Antenna"));
        assert!(script.contains("DrivenModal"));
    }
}
