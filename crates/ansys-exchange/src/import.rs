use crate::error::AnsysExchangeError;
use crate::model::{AnsysDesign, AnsysDesignKind, AnsysProject, AnsysSolutionType};

#[cfg(not(target_arch = "wasm32"))]
use std::path::Path;

pub fn import_aedt_str(content: &str) -> Result<AnsysProject, AnsysExchangeError> {
    let trimmed = content.trim();
    if trimmed.is_empty() {
        return Err(AnsysExchangeError::Parse("empty input".to_string()));
    }

    let project_name = parse_project_name(content).unwrap_or_else(|| "ImportedProject".to_string());
    let mut designs = parse_insert_designs(content);

    if designs.is_empty() {
        designs = parse_keyword_designs(content);
    }

    if designs.is_empty() {
        return Err(AnsysExchangeError::Unsupported(
            "no HFSS/Q3D design markers found in input".to_string(),
        ));
    }

    Ok(AnsysProject {
        name: project_name,
        designs,
    })
}

#[cfg(not(target_arch = "wasm32"))]
pub fn import_aedt_file(path: &Path) -> Result<AnsysProject, AnsysExchangeError> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| AnsysExchangeError::Io(e.to_string()))?;
    import_aedt_str(&content)
}

pub fn import_manifest_json_str(content: &str) -> Result<AnsysProject, AnsysExchangeError> {
    serde_json::from_str(content).map_err(|e| AnsysExchangeError::Json(e.to_string()))
}

fn parse_project_name(content: &str) -> Option<String> {
    for line in content.lines() {
        if let Some(v) = take_after_equals(line, "ProjectName") {
            if !v.is_empty() {
                return Some(v);
            }
        }
        if line.contains("NewProject") {
            if let Some(v) = first_quoted(line) {
                return Some(v);
            }
        }
    }
    None
}

fn parse_insert_designs(content: &str) -> Vec<AnsysDesign> {
    let mut out = Vec::new();
    for line in content.lines() {
        if !line.contains("InsertDesign") {
            continue;
        }
        let quoted = all_quoted(line);
        if quoted.len() < 2 {
            continue;
        }
        let tool = quoted[0].to_ascii_lowercase();
        let name = quoted[1].clone();
        let solution = quoted
            .get(2)
            .map(|s| parse_solution_type(s))
            .unwrap_or(AnsysSolutionType::Unknown);

        if tool.contains("hfss") {
            out.push(AnsysDesign {
                name,
                kind: AnsysDesignKind::Hfss,
                solution_type: if solution == AnsysSolutionType::Unknown {
                    AnsysSolutionType::DrivenModal
                } else {
                    solution
                },
                variables: Default::default(),
            });
        } else if tool.contains("q3d") {
            out.push(AnsysDesign {
                name,
                kind: AnsysDesignKind::Q3d,
                solution_type: if solution == AnsysSolutionType::Unknown {
                    AnsysSolutionType::Q3dC
                } else {
                    solution
                },
                variables: Default::default(),
            });
        }
    }
    out
}

fn parse_keyword_designs(content: &str) -> Vec<AnsysDesign> {
    let mut out = Vec::new();
    let mut hfss_count = 0usize;
    let mut q3d_count = 0usize;

    for line in content.lines() {
        let lower = line.to_ascii_lowercase();
        if lower.contains("hfss") && !lower.contains("q3d") {
            hfss_count += 1;
            out.push(AnsysDesign {
                name: format!("HFSSDesign{}", hfss_count),
                kind: AnsysDesignKind::Hfss,
                solution_type: parse_solution_type(&lower),
                variables: Default::default(),
            });
        }
        if lower.contains("q3d") {
            q3d_count += 1;
            out.push(AnsysDesign {
                name: format!("Q3DDesign{}", q3d_count),
                kind: AnsysDesignKind::Q3d,
                solution_type: parse_solution_type(&lower),
                variables: Default::default(),
            });
        }
    }

    dedup_by_name_kind(out)
}

fn dedup_by_name_kind(items: Vec<AnsysDesign>) -> Vec<AnsysDesign> {
    let mut out: Vec<AnsysDesign> = Vec::new();
    for d in items {
        if out.iter().any(|x| x.name == d.name && x.kind == d.kind) {
            continue;
        }
        out.push(d);
    }
    out
}

fn parse_solution_type(s: &str) -> AnsysSolutionType {
    let lower = s.to_ascii_lowercase();

    if lower.contains("drivenmodal") || lower.contains("driven modal") {
        return AnsysSolutionType::DrivenModal;
    }
    if lower.contains("driventerminal") || lower.contains("driven terminal") {
        return AnsysSolutionType::DrivenTerminal;
    }
    if lower.contains("eigenmode") {
        return AnsysSolutionType::Eigenmode;
    }
    if lower.contains("transient") {
        return AnsysSolutionType::Transient;
    }
    if lower.contains("sbr") {
        return AnsysSolutionType::SbrPlus;
    }
    if lower.contains("ac rl") || lower.contains("acrl") {
        return AnsysSolutionType::Q3dAcrl;
    }
    if lower.contains("dc rl") || lower.contains("dcrl") {
        return AnsysSolutionType::Q3dDcrl;
    }
    if lower.contains("conductance") || lower.contains("cg") {
        return AnsysSolutionType::Q3dCg;
    }
    if lower.contains("capacitance") || lower.contains(" q3d c") {
        return AnsysSolutionType::Q3dC;
    }

    AnsysSolutionType::Unknown
}

fn take_after_equals(line: &str, key: &str) -> Option<String> {
    let key_pos = line.find(key)?;
    let rem = &line[key_pos + key.len()..];
    let eq_pos = rem.find('=')?;
    let val = rem[eq_pos + 1..].trim().trim_matches(';').trim();
    if val.is_empty() {
        return None;
    }
    let unquoted = val.trim_matches('"').trim_matches('\'').to_string();
    Some(unquoted)
}

fn first_quoted(line: &str) -> Option<String> {
    all_quoted(line).into_iter().next()
}

fn all_quoted(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut chars = line.char_indices();
    while let Some((i, c)) = chars.next() {
        if c != '"' {
            continue;
        }
        if let Some((j, _)) = chars.find(|(_, ch)| *ch == '"') {
            out.push(line[i + 1..j].to_string());
        } else {
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn import_from_insert_design_lines() {
        let input = r#"
oProject = oDesktop.NewProject("BoardChannel")
oProject.InsertDesign("HFSS", "Antenna", "DrivenModal", "")
oProject.InsertDesign("Q3D Extractor", "Busbar", "Q3D AC RL", "")
"#;

        let project = import_aedt_str(input).unwrap();
        assert_eq!(project.name, "BoardChannel");
        assert_eq!(project.designs.len(), 2);
        assert_eq!(project.designs[0].kind, AnsysDesignKind::Hfss);
        assert_eq!(project.designs[0].solution_type, AnsysSolutionType::DrivenModal);
        assert_eq!(project.designs[1].kind, AnsysDesignKind::Q3d);
        assert_eq!(project.designs[1].solution_type, AnsysSolutionType::Q3dAcrl);
    }

    #[test]
    fn import_falls_back_to_keyword_scan() {
        let input = "HFSS DrivenTerminal setup\nQ3D capacitance matrix\n";
        let project = import_aedt_str(input).unwrap();
        assert_eq!(project.designs.len(), 2);
        assert_eq!(project.designs[0].name, "HFSSDesign1");
        assert_eq!(project.designs[0].solution_type, AnsysSolutionType::DrivenTerminal);
        assert_eq!(project.designs[1].name, "Q3DDesign1");
        assert_eq!(project.designs[1].solution_type, AnsysSolutionType::Q3dC);
    }

    #[test]
    fn import_rejects_empty_input() {
        let err = import_aedt_str("   ").unwrap_err();
        assert!(matches!(err, AnsysExchangeError::Parse(_)));
    }
}
