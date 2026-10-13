//! The `rmsh` command-line interface library (P6: convert / refine / quality).
//!
//! Adjudications pinned here:
//!
//! - **Hand-rolled argument parsing, no clap.** The dependency charter
//!   (rmsh ROADMAP.md section 0, hard constraint) whitelists only `std`,
//!   rmsh-owned crates and open-source rcad kernel crates. clap is none of
//!   those, and a three-verb CLI is small enough that hand-rolled subcommand
//!   matching keeps the charter intact.
//! - **Output flavor by flags, never inferred from the file name.** `.msh`
//!   carries neither the MSH version nor the encoding, so inference would be
//!   guesswork. `--v2/--v4` (default `--v4`) and `--ascii/--binary`
//!   (default `--ascii`) select the writer; the defaults match what current
//!   Gmsh writes (MSH 4.1 ASCII).
//! - **Pinned stdout contract.** `refine` prints before/after quality stats
//!   tables; `quality` prints the verdict, per-defect-class counts, hints
//!   and a stats table (see [`stats_table`] / [`quality_report`]). The
//!   integration tests (`tests/cli_verbs.rs`) pin the bytes.
//! - **Exit codes:** 0 success (quality: watertight), 1 quality verdict
//!   NOT watertight (script-friendly), 2 usage error, 3 runtime error
//!   (I/O, parse, empty mesh, optimizer failure).

use std::io::Write;
use std::path::Path;

use rmsh_algo::mesh_optimize::{MeshQualityOptimizer, TetMeshQualityStats, mesh_quality_stats};
use rmsh_algo::traits::{MeshParams, OptimizeParams};
use rmsh_geo::watertight::{WatertightReport, audit_watertight};
use rmsh_io::{
    load_msh_from_path, save_msh_v2_binary_to_path, save_msh_v2_to_path,
    save_msh_v4_binary_to_path, save_msh_v4_to_path,
};
use rmsh_model::Mesh;

// ── exit codes ───────────────────────────────────────────────────────────────

/// Success (for `quality`: the mesh is watertight).
pub const EXIT_OK: i32 = 0;
/// `quality`: the mesh is NOT watertight.
pub const EXIT_NOT_WATERTIGHT: i32 = 1;
/// Usage error: unknown command/flag, missing arguments, bad values.
pub const EXIT_USAGE: i32 = 2;
/// Runtime error: I/O, parse, empty mesh, optimizer failure.
pub const EXIT_RUNTIME: i32 = 3;

// ── parsed command line ──────────────────────────────────────────────────────

/// MSH file version to write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MshVersion {
    /// MSH 2.2.
    V2,
    /// MSH 4.1 (default).
    V4,
}

/// Payload encoding to write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Encoding {
    /// ASCII payload (default).
    Ascii,
    /// Binary (little-endian) payload.
    Binary,
}

/// The output MSH flavor: version + encoding, default v4 ASCII.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MshFlavor {
    /// MSH file version.
    pub version: MshVersion,
    /// Payload encoding.
    pub encoding: Encoding,
}

impl Default for MshFlavor {
    fn default() -> Self {
        MshFlavor {
            version: MshVersion::V4,
            encoding: Encoding::Ascii,
        }
    }
}

impl MshFlavor {
    /// Human-readable flavor label used on the convert summary line.
    pub fn label(&self) -> &'static str {
        match (self.version, self.encoding) {
            (MshVersion::V2, Encoding::Ascii) => "msh v2.2 ascii",
            (MshVersion::V2, Encoding::Binary) => "msh v2.2 binary",
            (MshVersion::V4, Encoding::Ascii) => "msh v4.1 ascii",
            (MshVersion::V4, Encoding::Binary) => "msh v4.1 binary",
        }
    }
}

/// Size/round overrides for `rmsh refine`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct RefineParams {
    /// `--passes N`: optimization rounds (default 3; 0 skips).
    pub passes: Option<u32>,
    /// `--min S`: collapse threshold / minimum size (default auto).
    pub min: Option<f64>,
    /// `--max S`: maximum element size (default auto).
    pub max: Option<f64>,
    /// `--size S`: target element size; fills unset `--min` (S/10) and
    /// `--max` (2S), mirroring `MeshParams::with_size`.
    pub size: Option<f64>,
}

/// A parsed verb invocation.
#[derive(Debug, Clone, PartialEq)]
pub enum Command {
    /// `rmsh convert <in> <out> [flavor flags]`.
    Convert {
        /// Input mesh path (any MSH v2/v4, ASCII/binary).
        input: String,
        /// Output mesh path.
        output: String,
        /// Output flavor (flags; default v4 ASCII).
        flavor: MshFlavor,
    },
    /// `rmsh refine <in> <out> [size flags] [flavor flags]`.
    Refine {
        /// Input mesh path (must carry Tetrahedron4 elements).
        input: String,
        /// Output mesh path.
        output: String,
        /// Size/round overrides.
        params: RefineParams,
        /// Output flavor (flags; default v4 ASCII).
        flavor: MshFlavor,
    },
    /// `rmsh quality <in>`.
    Quality {
        /// Input mesh path.
        input: String,
    },
}

/// What `run` should do after a successful parse.
#[derive(Debug, Clone, PartialEq)]
pub enum Parsed {
    /// Execute a verb.
    Run(Command),
    /// Print the top-level help (stdout, exit 0).
    HelpTop,
    /// Print a verb's help (stdout, exit 0).
    HelpVerb(&'static str),
    /// Print the version line (stdout, exit 0).
    Version,
}

/// A usage error (always exit code 2). The message is printed after the
/// `error: ` prefix on stderr.
#[derive(Debug, Clone, PartialEq)]
pub struct UsageError(pub String);

impl std::fmt::Display for UsageError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

// ── argument parsing ─────────────────────────────────────────────────────────

/// Parse the process argv (including argv[0], the program name).
pub fn parse_args(argv: &[String]) -> Result<Parsed, UsageError> {
    let Some(first) = argv.first() else {
        return Err(UsageError("no command given (try 'rmsh --help')".into()));
    };
    match first.as_str() {
        "-h" | "--help" | "help" => Ok(Parsed::HelpTop),
        "-V" | "--version" => Ok(Parsed::Version),
        "convert" => parse_convert(&argv[1..]),
        "refine" => parse_refine(&argv[1..]),
        "quality" => parse_quality(&argv[1..]),
        other => Err(UsageError(format!(
            "unknown command '{other}' (try 'rmsh --help')"
        ))),
    }
}

/// Flags accepted by every verb: `-h/--help`; convert/refine also take the
/// output-flavor flags; refine additionally the numeric flags.
fn parse_verb_flags(
    verb: &'static str,
    rest: &[String],
    required: &str,
    arity: usize,
) -> Result<(Vec<String>, MshFlavor, RefineParams), UsageError> {
    // Normalize `--opt=value` into two tokens so a single loop handles both
    // spellings for the value flags.
    let mut tokens: Vec<String> = Vec::with_capacity(rest.len());
    for tok in rest {
        if let Some((name, value)) = tok.split_once('=') {
            if name.starts_with("--") {
                tokens.push(name.to_string());
                tokens.push(value.to_string());
                continue;
            }
        }
        tokens.push(tok.clone());
    }
    let rest: &[String] = &tokens;

    let mut positional = Vec::new();
    let mut flavor = MshFlavor::default();
    let mut seen_v2: Option<&str> = None;
    let mut seen_v4: Option<&str> = None;
    let mut seen_ascii: Option<&str> = None;
    let mut seen_binary: Option<&str> = None;
    let mut params = RefineParams::default();
    let takes_flavor = verb == "convert" || verb == "refine";
    let takes_numeric = verb == "refine";

    let mut i = 0usize;
    while i < rest.len() {
        let tok = rest[i].as_str();
        match tok {
            "--v2" | "--v4" | "--ascii" | "--binary" => {
                if !takes_flavor {
                    return Err(unrecognized(verb, tok));
                }
                let slot = match tok {
                    "--v2" => &mut seen_v2,
                    "--v4" => &mut seen_v4,
                    "--ascii" => &mut seen_ascii,
                    _ => &mut seen_binary,
                };
                if let Some(previous) = slot.replace(tok) {
                    return Err(UsageError(format!(
                        "conflicting options '{previous}' and '{tok}'"
                    )));
                }
            }
            "--passes" | "--min" | "--max" | "--size" => {
                if !takes_numeric {
                    return Err(unrecognized(verb, tok));
                }
                let Some(value) = rest.get(i + 1) else {
                    return Err(UsageError(format!(
                        "option '{tok}' requires a value (try 'rmsh {verb} --help')"
                    )));
                };
                let value = value.trim();
                match tok {
                    "--passes" => {
                        let parsed = value
                            .parse::<u32>()
                            .map_err(|_| bad_value(tok, value, "a non-negative integer"))?;
                        params.passes = Some(parsed);
                    }
                    _ => {
                        let parsed = value
                            .parse::<f64>()
                            .ok()
                            .filter(|v| v.is_finite() && *v > 0.0)
                            .ok_or_else(|| bad_value(tok, value, "a positive number"))?;
                        match tok {
                            "--min" => params.min = Some(parsed),
                            "--max" => params.max = Some(parsed),
                            _ => params.size = Some(parsed),
                        }
                    }
                }
                i += 1;
            }
            _ if tok.starts_with("--") => {
                return Err(unrecognized(verb, tok));
            }
            _ if tok.starts_with('-') && tok.len() > 1 => {
                return Err(unrecognized(verb, tok));
            }
            _ => positional.push(tok.to_string()),
        }
        i += 1;
    }

    if let (Some(a), Some(b)) = (seen_v2, seen_v4) {
        return Err(UsageError(format!(
            "conflicting options '{a}' and '{b}'"
        )));
    }
    if let (Some(a), Some(b)) = (seen_ascii, seen_binary) {
        return Err(UsageError(format!(
            "conflicting options '{a}' and '{b}'"
        )));
    }
    flavor.version = if seen_v2.is_some() {
        MshVersion::V2
    } else {
        MshVersion::V4
    };
    flavor.encoding = if seen_binary.is_some() {
        Encoding::Binary
    } else {
        Encoding::Ascii
    };

    if positional.len() != arity {
        return Err(UsageError(format!(
            "{verb} requires {required} (try 'rmsh {verb} --help')"
        )));
    }
    Ok((positional, flavor, params))
}

fn unrecognized(verb: &str, opt: &str) -> UsageError {
    UsageError(format!(
        "unrecognized option '{opt}' for '{verb}' (try 'rmsh {verb} --help')"
    ))
}

fn bad_value(opt: &str, value: &str, expected: &str) -> UsageError {
    UsageError(format!(
        "invalid value '{value}' for '{opt}': expected {expected}"
    ))
}

fn parse_convert(rest: &[String]) -> Result<Parsed, UsageError> {
    if wants_help(rest) {
        return Ok(Parsed::HelpVerb("convert"));
    }
    let (positional, flavor, _) = parse_verb_flags("convert", rest, "<in> <out>", 2)?;
    Ok(Parsed::Run(Command::Convert {
        input: positional[0].clone(),
        output: positional[1].clone(),
        flavor,
    }))
}

fn parse_refine(rest: &[String]) -> Result<Parsed, UsageError> {
    if wants_help(rest) {
        return Ok(Parsed::HelpVerb("refine"));
    }
    let (positional, flavor, params) = parse_verb_flags("refine", rest, "<in> <out>", 2)?;
    Ok(Parsed::Run(Command::Refine {
        input: positional[0].clone(),
        output: positional[1].clone(),
        params,
        flavor,
    }))
}

fn parse_quality(rest: &[String]) -> Result<Parsed, UsageError> {
    if wants_help(rest) {
        return Ok(Parsed::HelpVerb("quality"));
    }
    let (positional, _, _) = parse_verb_flags("quality", rest, "<in>", 1)?;
    Ok(Parsed::Run(Command::Quality {
        input: positional[0].clone(),
    }))
}

/// A verb help request anywhere in the verb's arguments.
fn wants_help(rest: &[String]) -> bool {
    rest.iter().any(|t| t == "-h" || t == "--help")
}

// ── help texts ───────────────────────────────────────────────────────────────

/// The top-level help (first line pinned by tests).
pub const MAIN_HELP: &str = "\
rmsh - mesh toolchain CLI (convert / refine / quality)

Usage:
  rmsh convert <in> <out> [--v2|--v4] [--ascii|--binary]
  rmsh refine <in> <out> [--passes N] [--min S] [--max S] [--size S] [--v2|--v4] [--ascii|--binary]
  rmsh quality <in>

Commands:
  convert  Convert a mesh between MSH formats (v2/v4, ASCII/binary)
  refine   Optimize tetrahedral mesh quality and print before/after stats
  quality  Audit watertightness and print mesh quality stats

Options:
  -h, --help     Print help
  -V, --version  Print version

Run 'rmsh <command> --help' for command-specific help.
";

const CONVERT_HELP: &str = "\
rmsh convert - convert a mesh between MSH formats

Usage: rmsh convert <in> <out> [--v2|--v4] [--ascii|--binary]

Reads MSH v2/v4, ASCII or binary (auto-detected from the payload).
Writes the flavor selected by flags; the default output is MSH v4.1
ASCII. Neither the version nor the encoding is inferred from the
output file name: '.msh' carries neither.

Options:
  --v2       Write MSH v2.2
  --v4       Write MSH v4.1 (default)
  --ascii    Write an ASCII payload (default)
  --binary   Write a binary payload
  -h, --help Print this help

Exit codes:
  0  success
  2  usage error
  3  read/write/parse error
";

const REFINE_HELP: &str = "\
rmsh refine - optimize tetrahedral mesh quality

Usage: rmsh refine <in> <out> [--passes N] [--min S] [--max S] [--size S] [--v2|--v4] [--ascii|--binary]

Runs the rmsh-algo MeshQualityOptimizer (edge collapse, 2-3/3-2 flips,
quality-gated relocation) over the Tetrahedron4 elements of the input
and writes the optimized mesh in the selected flavor. Prints a
before/after quality stats table to stdout (the pinned CLI contract):

  Quality stats (before)
    min_dihedral_deg : <f64, 6 decimals>
    p95_radius_edge  : <f64, 6 decimals>
    sliver_fraction  : <f64, 6 decimals>
    tet_count        : <integer>

Options:
  --passes N  Optimization rounds (default 3; 0 skips optimization)
  --min S     Collapse threshold / minimum size (default: auto, 0.5x mean edge)
  --max S     Maximum element size carried on the mesh params (default: auto)
  --size S    Target element size; fills unset --min (S/10) and --max (2S)
  --v2/--v4, --ascii/--binary  Output flavor (default: MSH v4.1 ASCII)
  -h, --help  Print this help

Exit codes:
  0  success
  2  usage error
  3  read/write/parse/optimizer error
";

const QUALITY_HELP: &str = "\
rmsh quality - audit watertightness and mesh quality

Usage: rmsh quality <in>

Reads MSH v2/v4 (ASCII or binary, auto-detected), runs the rmsh-geo
watertight audit (holes, non-manifold edges, orientation, duplicates,
degenerates, self-intersections) and the tet quality stats. Prints the
verdict, per-defect-class counts, one actionable hint per defect and
the quality stats table.

Exit codes:
  0  mesh is watertight
  1  mesh is NOT watertight (script-friendly signal)
  2  usage error
  3  read/parse error
";

fn verb_help(verb: &str) -> &'static str {
    match verb {
        "convert" => CONVERT_HELP,
        "refine" => REFINE_HELP,
        "quality" => QUALITY_HELP,
        _ => MAIN_HELP,
    }
}

// ── pinned table rendering ───────────────────────────────────────────────────

const STATS_LABEL_WIDTH: usize = 16;
const DEFECT_LABEL_WIDTH: usize = 22;

fn push_row(out: &mut String, label: &str, value: &str, width: usize) {
    out.push_str("  ");
    out.push_str(&format!("{:<width$}", label, width = width));
    out.push_str(" : ");
    out.push_str(value);
    out.push('\n');
}

/// Render the pinned quality stats table. Byte format:
/// title line, then one row per metric as
/// `  {label:<16} : {value}` with the three ratios in `{:.6}` form and
/// `tet_count` as a plain integer (`inf` prints as `inf`).
pub fn stats_table(title: &str, s: &TetMeshQualityStats) -> String {
    let mut out = String::new();
    out.push_str(title);
    out.push('\n');
    push_row(
        &mut out,
        "min_dihedral_deg",
        &format!("{:.6}", s.min_dihedral_deg),
        STATS_LABEL_WIDTH,
    );
    push_row(
        &mut out,
        "p95_radius_edge",
        &format!("{:.6}", s.p95_radius_edge),
        STATS_LABEL_WIDTH,
    );
    push_row(
        &mut out,
        "sliver_fraction",
        &format!("{:.6}", s.sliver_fraction),
        STATS_LABEL_WIDTH,
    );
    push_row(
        &mut out,
        "tet_count",
        &s.tet_count.to_string(),
        STATS_LABEL_WIDTH,
    );
    out
}

/// Render the full `quality` report: verdict line, per-defect-class counts
/// (`  {label:<22} : {count}`), one hint per line when defects exist, then
/// the pinned stats table.
pub fn quality_report(report: &WatertightReport, stats: &TetMeshQualityStats) -> String {
    let mut out = String::new();
    out.push_str(if report.is_watertight {
        "watertight: yes\n"
    } else {
        "watertight: no\n"
    });
    if !report.is_watertight {
        out.push_str("defects\n");
        push_row(
            &mut out,
            "boundary_loops",
            &report.boundary_loops.len().to_string(),
            DEFECT_LABEL_WIDTH,
        );
        push_row(
            &mut out,
            "non_manifold_edges",
            &report.non_manifold_edges.len().to_string(),
            DEFECT_LABEL_WIDTH,
        );
        push_row(
            &mut out,
            "orientation_mismatches",
            &report.orientation_mismatches.len().to_string(),
            DEFECT_LABEL_WIDTH,
        );
        push_row(
            &mut out,
            "duplicate_triangles",
            &report.duplicate_triangles.len().to_string(),
            DEFECT_LABEL_WIDTH,
        );
        push_row(
            &mut out,
            "degenerate_triangles",
            &report.degenerate_triangles.len().to_string(),
            DEFECT_LABEL_WIDTH,
        );
        push_row(
            &mut out,
            "self_intersections",
            &report.self_intersections.len().to_string(),
            DEFECT_LABEL_WIDTH,
        );
        if !report.hints.is_empty() {
            out.push_str("hints\n");
            for hint in &report.hints {
                out.push_str("  ");
                out.push_str(hint);
                out.push('\n');
            }
        }
    }
    out.push_str(&stats_table("Quality stats", stats));
    out
}

// ── execution ────────────────────────────────────────────────────────────────

/// Run with explicit argv and writers; returns the process exit code.
/// `argv` includes the program name at index 0 (as `std::env::args` does).
pub fn run(argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> i32 {
    let args: &[String] = if argv.is_empty() { &[] } else { &argv[1..] };
    match parse_args(args) {
        Err(e) => {
            let _ = writeln!(err, "error: {e}");
            EXIT_USAGE
        }
        Ok(Parsed::Version) => {
            let _ = writeln!(out, "rmsh {}", env!("CARGO_PKG_VERSION"));
            EXIT_OK
        }
        Ok(Parsed::HelpTop) => {
            let _ = write!(out, "{MAIN_HELP}");
            EXIT_OK
        }
        Ok(Parsed::HelpVerb(verb)) => {
            let _ = write!(out, "{}", verb_help(verb));
            EXIT_OK
        }
        Ok(Parsed::Run(command)) => execute(&command, out, err),
    }
}

/// Run against the process environment (argv + stdout/stderr), flushing
/// before returning the exit code.
pub fn run_from_env() -> i32 {
    let argv: Vec<String> = std::env::args().collect();
    let stdout = std::io::stdout();
    let stderr = std::io::stderr();
    let mut out = stdout.lock();
    let mut err = stderr.lock();
    let code = run(&argv, &mut out, &mut err);
    let _ = out.flush();
    let _ = err.flush();
    code
}

/// Load a mesh for any verb: shared error shape and the empty-mesh guard
/// (the MSH parser accepts section-less files as empty meshes, which the
/// CLI must reject as garbage input).
fn load_mesh(path: &str, err: &mut dyn Write) -> Option<Mesh> {
    match load_msh_from_path(Path::new(path)) {
        Ok(mesh) => {
            if mesh.node_count() == 0 && mesh.element_count() == 0 {
                let _ = writeln!(err, "error: no mesh data found in '{path}'");
                None
            } else {
                Some(mesh)
            }
        }
        Err(e) => {
            let _ = writeln!(err, "error: cannot read '{path}': {e}");
            None
        }
    }
}

fn save_mesh(path: &str, flavor: &MshFlavor, mesh: &Mesh) -> Result<(), rmsh_io::MshError> {
    let p = Path::new(path);
    match (flavor.version, flavor.encoding) {
        (MshVersion::V2, Encoding::Ascii) => save_msh_v2_to_path(p, mesh),
        (MshVersion::V2, Encoding::Binary) => save_msh_v2_binary_to_path(p, mesh),
        (MshVersion::V4, Encoding::Ascii) => save_msh_v4_to_path(p, mesh),
        (MshVersion::V4, Encoding::Binary) => save_msh_v4_binary_to_path(p, mesh),
    }
}

fn execute(command: &Command, out: &mut dyn Write, err: &mut dyn Write) -> i32 {
    match command {
        Command::Convert { input, output, flavor } => run_convert(input, output, *flavor, out, err),
        Command::Refine { input, output, params, flavor } => {
            run_refine(input, output, params, *flavor, out, err)
        }
        Command::Quality { input } => run_quality(input, out, err),
    }
}

fn run_convert(
    input: &str,
    output: &str,
    flavor: MshFlavor,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    let Some(mesh) = load_mesh(input, err) else {
        return EXIT_RUNTIME;
    };
    if let Err(e) = save_mesh(output, &flavor, &mesh) {
        let _ = writeln!(err, "error: cannot write '{output}': {e}");
        return EXIT_RUNTIME;
    }
    let _ = writeln!(
        out,
        "wrote {output}: {} nodes, {} elements ({})",
        mesh.node_count(),
        mesh.element_count(),
        flavor.label()
    );
    EXIT_OK
}

fn run_refine(
    input: &str,
    output: &str,
    params: &RefineParams,
    flavor: MshFlavor,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    let Some(mesh) = load_mesh(input, err) else {
        return EXIT_RUNTIME;
    };
    let before = mesh_quality_stats(&mesh);

    // MeshParams assembly: explicit flags win; --size fills unset --min
    // (S/10) and --max (2S) per MeshParams::with_size semantics; everything
    // unset stays 0 so the optimizer applies its documented auto defaults
    // (collapse threshold 0.5x mean edge length).
    let size = params.size.unwrap_or(0.0);
    let mesh_params = MeshParams {
        element_size: size,
        min_size: params.min.or(params.size.map(|s| s / 10.0)).unwrap_or(0.0),
        max_size: params.max.or(params.size.map(|s| s * 2.0)).unwrap_or(0.0),
        optimize_passes: params.passes.unwrap_or(3),
    };

    let mut mesh = mesh;
    let optimizer = MeshQualityOptimizer::new();
    if let Err(e) = optimizer.optimize_with_params(&mut mesh, &OptimizeParams::default(), &mesh_params)
    {
        let _ = writeln!(err, "error: optimization failed: {e}");
        return EXIT_RUNTIME;
    }
    let after = mesh_quality_stats(&mesh);

    let _ = write!(
        out,
        "{}{}",
        stats_table("Quality stats (before)", &before),
        stats_table("Quality stats (after)", &after)
    );

    if let Err(e) = save_mesh(output, &flavor, &mesh) {
        let _ = writeln!(err, "error: cannot write '{output}': {e}");
        return EXIT_RUNTIME;
    }
    EXIT_OK
}

fn run_quality(input: &str, out: &mut dyn Write, err: &mut dyn Write) -> i32 {
    let Some(mesh) = load_mesh(input, err) else {
        return EXIT_RUNTIME;
    };
    let report = audit_watertight(&mesh);
    let stats = mesh_quality_stats(&mesh);
    let _ = write!(out, "{}", quality_report(&report, &stats));
    if report.is_watertight {
        EXIT_OK
    } else {
        EXIT_NOT_WATERTIGHT
    }
}

// ── unit tests (rendering + parsing; the process contract lives in
//    tests/cli_verbs.rs) ──────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(args: &[&str]) -> Vec<String> {
        // argv INCLUDING the program name at index 0 (as run() expects).
        std::iter::once("rmsh".to_string())
            .chain(args.iter().map(|s| s.to_string()))
            .collect()
    }

    fn args(args: &[&str]) -> Vec<String> {
        // Arguments AFTER the program name (as parse_args expects).
        args.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn stats_table_is_byte_pinned() {
        let stats = TetMeshQualityStats {
            min_dihedral_deg: 20.605215,
            p95_radius_edge: 1.2345678,
            sliver_fraction: 0.0,
            tet_count: 96,
        };
        let expected = "\
Quality stats
  min_dihedral_deg : 20.605215
  p95_radius_edge  : 1.234568
  sliver_fraction  : 0.000000
  tet_count        : 96
";
        assert_eq!(stats_table("Quality stats", &stats), expected);
    }

    #[test]
    fn empty_tet_stats_render_inf_p95() {
        let stats = TetMeshQualityStats {
            min_dihedral_deg: 0.0,
            p95_radius_edge: f64::INFINITY,
            sliver_fraction: 1.0,
            tet_count: 0,
        };
        let table = stats_table("Quality stats", &stats);
        assert!(table.contains("  p95_radius_edge  : inf\n"), "{table:?}");
        assert!(table.contains("  sliver_fraction  : 1.000000\n"));
        assert!(table.contains("  tet_count        : 0\n"));
    }

    #[test]
    fn quality_report_watertight_shape() {
        let report = WatertightReport {
            is_watertight: true,
            ..WatertightReport::default()
        };
        let stats = TetMeshQualityStats {
            min_dihedral_deg: 0.0,
            p95_radius_edge: f64::INFINITY,
            sliver_fraction: 1.0,
            tet_count: 0,
        };
        let text = quality_report(&report, &stats);
        assert!(text.starts_with("watertight: yes\nQuality stats\n"));
        assert!(!text.contains("defects"), "no defects section when clean");
    }

    #[test]
    fn usage_error_paths() {
        let mut out = Vec::new();
        let mut err = Vec::new();
        assert_eq!(run(&argv(&[]), &mut out, &mut err), EXIT_USAGE);
        assert!(out.is_empty());
        let msg = String::from_utf8(err).unwrap();
        assert!(msg.contains("no command given"), "{msg:?}");

        let mut err = Vec::new();
        assert_eq!(run(&argv(&["wat"]), &mut out, &mut err), EXIT_USAGE);
        assert!(String::from_utf8(err).unwrap().contains("unknown command 'wat'"));
    }

    #[test]
    fn version_and_help_in_memory() {
        let mut out = Vec::new();
        let mut err = Vec::new();
        assert_eq!(run(&argv(&["--version"]), &mut out, &mut err), EXIT_OK);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            format!("rmsh {}\n", env!("CARGO_PKG_VERSION"))
        );
        let mut out = Vec::new();
        assert_eq!(run(&argv(&["--help"]), &mut out, &mut err), EXIT_OK);
        assert!(
            String::from_utf8(out)
                .unwrap()
                .starts_with("rmsh - mesh toolchain CLI (convert / refine / quality)\n")
        );
    }

    #[test]
    fn flavor_flags_parse_and_conflict() {
        let parsed = parse_args(&args(&[
            "convert",
            "a.msh",
            "b.msh",
            "--v2",
            "--binary",
        ]))
        .expect("parses");
        match parsed {
            Parsed::Run(Command::Convert { flavor, .. }) => {
                assert_eq!(flavor, MshFlavor { version: MshVersion::V2, encoding: Encoding::Binary });
                assert_eq!(flavor.label(), "msh v2.2 binary");
            }
            other => panic!("unexpected parse: {other:?}"),
        }
        assert!(parse_args(&args(&["convert", "a", "b", "--v2", "--v4"])).is_err());
        assert!(parse_args(&args(&["quality", "a", "--v2"])).is_err());
    }

    #[test]
    fn refine_numeric_flags_parse() {
        let parsed = parse_args(&args(&[
            "refine",
            "in.msh",
            "out.msh",
            "--passes",
            "5",
            "--min=0.25",
            "--size",
            "1.5",
        ]))
        .expect("parses");
        match parsed {
            Parsed::Run(Command::Refine { params, .. }) => {
                assert_eq!(params.passes, Some(5));
                assert_eq!(params.min, Some(0.25));
                assert_eq!(params.size, Some(1.5));
                assert_eq!(params.max, None);
            }
            other => panic!("unexpected parse: {other:?}"),
        }
        // Bad values are usage errors.
        assert!(parse_args(&args(&["refine", "a", "b", "--passes", "many"])).is_err());
        assert!(parse_args(&args(&["refine", "a", "b", "--min", "0"])).is_err());
        assert!(parse_args(&args(&["refine", "a", "b", "--min", "-1"])).is_err());
    }
}
