//! Map emstudio Design to rem PalaceConfig JSON.
//!
//! We generate a JSON file on disk that rem's `load_config()` reads.
//! This provides insulation against rem API changes and makes the config
//! independently inspectable and testable.

use std::path::{Path, PathBuf};

use serde_json::{json, Value};

use emstudio_domain::Design;
use emstudio_domain::solution_type::SolutionType;
use emstudio_domain::boundary::BoundaryType;
use emstudio_domain::excitation::ExcitationType;
use emstudio_domain::variable::PropertyValue;

use crate::error::SolverError;
use crate::mesh_bridge::parse_frequency;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Generate a rem-compatible Palace JSON config from an EMStudio Design and
/// write it to `output_dir/palace_config.json`.
///
/// Returns the path to the written config file.
pub fn write_palace_config(
    design: &Design,
    mesh_path: &Path,
    output_dir: &Path,
) -> Result<PathBuf, SolverError> {
    let config_json = build_palace_json(design, mesh_path, output_dir)?;

    let config_path = output_dir.join("palace_config.json");
    let content = serde_json::to_string_pretty(&config_json)
        .map_err(|e| SolverError::ConfigGeneration(e.to_string()))?;
    std::fs::write(&config_path, content)
        .map_err(|e| SolverError::io(&config_path, e))?;

    Ok(config_path)
}

// ---------------------------------------------------------------------------
// JSON construction
// ---------------------------------------------------------------------------

fn build_palace_json(
    design: &Design,
    mesh_path: &Path,
    output_dir: &Path,
) -> Result<Value, SolverError> {
    let problem_type = map_solution_type(&design.solution_type);
    let unit_scale = crate::mesh_bridge::unit_to_meters(&design.units);

    let mut config = json!({
        "Problem": {
            "Type": problem_type,
            "Verbose": 1,
            "Output": output_dir.to_string_lossy(),
        },
        "Model": {
            "Mesh": mesh_path.to_string_lossy(),
            "L0": unit_scale,
        },
    });

    // ── Domains (materials) ──────────────────────────────────────────────
    let materials = build_materials(design);
    if !materials.is_empty() {
        config["Domains"] = json!({ "Materials": materials });
    }

    // ── Boundaries ───────────────────────────────────────────────────────
    let boundaries = build_boundaries(design);
    config["Boundaries"] = boundaries;

    // ── Solver ───────────────────────────────────────────────────────────
    let solver = build_solver(design)?;
    config["Solver"] = solver;

    Ok(config)
}

/// Map EMStudio SolutionType to rem ProblemType string.
fn map_solution_type(st: &SolutionType) -> &'static str {
    match st {
        SolutionType::DrivenModal | SolutionType::DrivenTerminal => "Driven",
        SolutionType::Eigenmode => "Eigenmode",
        SolutionType::Transient => "Transient",
        SolutionType::SBRPlus => "SBR",
        SolutionType::Q3D_C | SolutionType::Q3D_CG => "Electrostatic",
        SolutionType::Q3D_DCRL | SolutionType::Q3D_ACRL => "Magnetostatic",
    }
}

// ---------------------------------------------------------------------------
// Materials
// ---------------------------------------------------------------------------

fn build_materials(design: &Design) -> Vec<Value> {
    design
        .definitions
        .materials
        .iter()
        .enumerate()
        .map(|(i, mat)| {
            let eps = extract_constant(&mat.properties.permittivity, 1.0);
            let mu = extract_constant(&mat.properties.permeability, 1.0);
            let sigma = extract_constant(&mat.properties.conductivity, 0.0);
            let loss_tan = extract_constant(&mat.properties.dielectric_loss_tangent, 0.0);

            json!({
                "Attributes": [i as u32 + 1],
                "Permittivity": eps,
                "Permeability": mu,
                "Conductivity": sigma,
                "LossTan": loss_tan,
            })
        })
        .collect()
}

fn extract_constant(pv: &PropertyValue, default: f64) -> f64 {
    match pv {
        PropertyValue::Constant { value } => *value,
        _ => default,
    }
}

// ---------------------------------------------------------------------------
// Boundaries
// ---------------------------------------------------------------------------

fn build_boundaries(design: &Design) -> Value {
    let mut pec_attrs: Vec<u32> = Vec::new();
    let mut pmc_attrs: Vec<u32> = Vec::new();
    let mut absorbing_attrs: Vec<u32> = Vec::new();
    let mut ground_attrs: Vec<u32> = Vec::new();
    let mut lumped_ports: Vec<Value> = Vec::new();
    let mut wave_ports: Vec<Value> = Vec::new();
    let mut terminals: Vec<Value> = Vec::new();

    // Assign sequential boundary surface tags.
    let mut next_bc_tag: u32 = 100;

    for bc in &design.boundaries {
        let tag = next_bc_tag;
        next_bc_tag += 1;

        match bc.boundary_type {
            BoundaryType::PerfectE => pec_attrs.push(tag),
            BoundaryType::PerfectH => pmc_attrs.push(tag),
            BoundaryType::Radiation | BoundaryType::PML => absorbing_attrs.push(tag),
            BoundaryType::Symmetry => pec_attrs.push(tag), // treated as PEC symmetry
            BoundaryType::InfiniteGroundPlane => ground_attrs.push(tag),
            _ => {} // FiniteConductivity, Impedance, MasterSlave, etc. — extend later
        }
    }

    // Excitations → ports.
    let mut port_index: u32 = 1;
    for exc in &design.excitations {
        let tag = next_bc_tag;
        next_bc_tag += 1;

        match exc.excitation_type {
            ExcitationType::LumpedPort => {
                let r = exc.properties.get("impedance")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(50.0);
                lumped_ports.push(json!({
                    "Index": port_index,
                    "Attributes": [tag],
                    "R": r,
                    "Excitation": true,
                }));
                port_index += 1;
            }
            ExcitationType::WavePort => {
                wave_ports.push(json!({
                    "Index": port_index,
                    "Attributes": [tag],
                    "Excitation": true,
                    "Mode": 1,
                }));
                port_index += 1;
            }
            ExcitationType::Source => {
                terminals.push(json!({
                    "Index": port_index,
                    "Attributes": [tag],
                }));
                port_index += 1;
            }
            _ => {}
        }
    }

    let mut boundaries = json!({});
    if !pec_attrs.is_empty() {
        boundaries["PEC"] = json!({ "Attributes": pec_attrs });
    }
    if !pmc_attrs.is_empty() {
        boundaries["PMC"] = json!({ "Attributes": pmc_attrs });
    }
    if !absorbing_attrs.is_empty() {
        boundaries["Absorbing"] = json!({ "Attributes": absorbing_attrs, "Order": 1 });
    }
    if !ground_attrs.is_empty() {
        boundaries["Ground"] = json!({ "Attributes": ground_attrs });
    }
    if !lumped_ports.is_empty() {
        boundaries["LumpedPort"] = Value::Array(lumped_ports);
    }
    if !wave_ports.is_empty() {
        boundaries["WavePort"] = Value::Array(wave_ports);
    }
    if !terminals.is_empty() {
        boundaries["Terminal"] = Value::Array(terminals);
    }

    boundaries
}

// ---------------------------------------------------------------------------
// Solver section
// ---------------------------------------------------------------------------

fn build_solver(design: &Design) -> Result<Value, SolverError> {
    let mut solver = json!({ "Order": 1 });

    // Find the first enabled setup.
    let setup = design
        .analysis_setups
        .iter()
        .find(|s| s.enabled)
        .ok_or_else(|| SolverError::ConfigGeneration("No enabled analysis setup".into()))?;

    // Linear solver defaults.
    solver["Linear"] = json!({
        "Type": "GMRES",
        "Tol": 1e-6,
        "MaxIter": 200,
    });

    match &design.solution_type {
        SolutionType::DrivenModal | SolutionType::DrivenTerminal => {
            // Build Driven solver section from frequency sweeps.
            if let Some(sweep) = setup.frequency_sweeps.first() {
                let min_freq = parse_frequency(&sweep.start)?;
                let max_freq = parse_frequency(&sweep.stop)?;
                let freq_step = sweep
                    .step
                    .as_deref()
                    .map(parse_frequency)
                    .transpose()?
                    .unwrap_or((max_freq - min_freq) / 10.0);

                solver["Driven"] = json!({
                    "MinFreq": min_freq,
                    "MaxFreq": max_freq,
                    "FreqStep": freq_step,
                    "SaveStep": 1,
                });
            } else {
                // Single-frequency solve at solution frequency.
                let freq = parse_frequency(&setup.solution_frequency)?;
                solver["Driven"] = json!({
                    "MinFreq": freq,
                    "MaxFreq": freq,
                    "FreqStep": 0.0,
                    "SaveStep": 1,
                });
            }
        }
        SolutionType::Eigenmode => {
            solver["Eigenmode"] = json!({
                "N": 10,
                "Tol": 1e-6,
                "MaxIter": 200,
                "Target": parse_frequency(&setup.solution_frequency).unwrap_or(0.0),
            });
        }
        SolutionType::Transient => {
            solver["Transient"] = json!({
                "Type": "GeneralizedAlpha",
                "MaxTime": 1e-9,
                "TimeStep": 1e-12,
                "SaveStep": 10,
            });
        }
        SolutionType::SBRPlus => {
            let freq = parse_frequency(&setup.solution_frequency).unwrap_or(1e9);
            solver["SBR"] = json!({
                "FreqMin": freq,
                "FreqMax": freq,
                "FreqStep": 0.0,
                "RayDensity": 1e4,
                "MaxBounces": 5,
                "ThetaInc": 0.0,
                "PhiInc": 0.0,
                "Polarization": "theta",
            });
        }
        SolutionType::Q3D_C | SolutionType::Q3D_CG => {
            solver["Electrostatic"] = json!({ "Save": 1 });
        }
        SolutionType::Q3D_DCRL | SolutionType::Q3D_ACRL => {
            solver["Magnetostatic"] = json!({ "Save": 1 });
        }
    }

    Ok(solver)
}

#[cfg(test)]
mod tests {
    use super::*;
    use emstudio_domain::SolutionType;

    #[test]
    fn solution_type_mapping() {
        assert_eq!(map_solution_type(&SolutionType::DrivenModal), "Driven");
        assert_eq!(map_solution_type(&SolutionType::Q3D_C), "Electrostatic");
        assert_eq!(map_solution_type(&SolutionType::SBRPlus), "SBR");
        assert_eq!(map_solution_type(&SolutionType::Eigenmode), "Eigenmode");
    }
}
