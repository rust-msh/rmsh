//! Integration tests for the `rmsh` CLI (P6: convert / refine / quality).
//!
//! These tests are the PINNED product-surface contract: they drive the built
//! binary through `std::process::Command` and assert byte-exact stdout table
//! lines, stderr error messages and exit codes. Fixtures are built
//! programmatically with rmsh-model + rmsh-io writers into a per-test temp
//! directory (no binary fixtures in the repo).
//!
//! Exit-code contract:
//!   0  success (quality: watertight)
//!   1  quality: mesh is NOT watertight
//!   2  usage error (unknown command/flag, missing args, bad values)
//!   3  runtime error (I/O, parse, empty mesh, optimizer failure)

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use rmsh_algo::mesh_optimize::mesh_quality_stats;
use rmsh_io::load_msh_from_path;
use rmsh_model::{Element, ElementType, Mesh, Node};

// ── helpers ──────────────────────────────────────────────────────────────────

/// Run the built `rmsh` binary with the given arguments.
fn run_cli(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rmsh"))
        .args(args)
        .output()
        .expect("the rmsh binary must be built and runnable")
}

fn stdout_of(output: &Output) -> String {
    String::from_utf8(output.stdout.clone()).expect("stdout is UTF-8")
}

fn stderr_of(output: &Output) -> String {
    String::from_utf8(output.stderr.clone()).expect("stderr is UTF-8")
}

/// Per-test scratch directory under the system temp dir.
struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn new(name: &str) -> Self {
        let path = std::env::temp_dir().join(format!(
            "rmsh_cli_test_{}_{}_{}",
            std::process::id(),
            name,
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("clock is after epoch")
                .as_nanos()
        ));
        fs::create_dir_all(&path).expect("create temp dir");
        TempDir { path }
    }

    fn file(&self, name: &str) -> PathBuf {
        self.path.join(name)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

/// Write `text` to `path` as UTF-8.
fn write_text(path: &Path, text: &str) {
    let mut f = fs::File::create(path).expect("create fixture file");
    f.write_all(text.as_bytes()).expect("write fixture file");
}

/// A closed, consistently outward-oriented unit cube: 8 nodes, 12 triangles.
/// Every edge is traversed exactly once in each direction, so the watertight
/// audit reports zero defects.
fn watertight_cube() -> Mesh {
    let mut mesh = Mesh::new();
    for (id, (x, y, z)) in [
        (1u64, (0.0, 0.0, 0.0)),
        (2, (1.0, 0.0, 0.0)),
        (3, (1.0, 1.0, 0.0)),
        (4, (0.0, 1.0, 0.0)),
        (5, (0.0, 0.0, 1.0)),
        (6, (1.0, 0.0, 1.0)),
        (7, (1.0, 1.0, 1.0)),
        (8, (0.0, 1.0, 1.0)),
    ] {
        mesh.add_node(Node::new(id, x, y, z));
    }
    for (id, tri) in [
        (1u64, [1u64, 3, 2]), // bottom (z=0), outward -z
        (2, [1, 4, 3]),
        (3, [5, 6, 7]), // top (z=1), outward +z
        (4, [5, 7, 8]),
        (5, [1, 2, 6]), // front (y=0), outward -y
        (6, [1, 6, 5]),
        (7, [2, 3, 7]), // right (x=1), outward +x
        (8, [2, 7, 6]),
        (9, [3, 4, 8]), // back (y=1), outward +y
        (10, [3, 8, 7]),
        (11, [4, 1, 5]), // left (x=0), outward -x
        (12, [4, 5, 8]),
    ] {
        mesh.add_element(Element::new(id, ElementType::Triangle3, tri.to_vec()));
    }
    mesh
}

/// The watertight cube with the top-face triangle element 4 removed: its
/// 3 edges become a closed boundary loop (a hole the audit must report).
fn punched_cube() -> Mesh {
    let mut mesh = watertight_cube();
    mesh.elements.retain(|e| e.id != 4);
    mesh
}

/// Two tets sharing face (1,2,3): a flat cap (node 4 at height 0.1) and a
/// tall tet below (node 5). The shared face admits a quality-improving 2-3
/// swap (the same construction the mesh_optimize unit tests prove out), and
/// no edge is shorter than 0.6, so `--min 0.3` disables collapses entirely.
fn sliver_pair() -> Mesh {
    let mut mesh = Mesh::new();
    for (id, (x, y, z)) in [
        (1u64, (0.0, 0.0, 0.0)),
        (2, (1.0, 0.0, 0.0)),
        (3, (0.5, 0.9, 0.0)),
        (4, (0.5, 0.3, 0.1)),
        (5, (0.5, 0.3, -1.5)),
    ] {
        mesh.add_node(Node::new(id, x, y, z));
    }
    mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));
    mesh.add_element(Element::new(2, ElementType::Tetrahedron4, vec![1, 2, 3, 5]));
    mesh
}

fn write_msh_v4_ascii(path: &Path, mesh: &Mesh) {
    rmsh_io::save_msh_v4_to_path(path, mesh).expect("write v4 ascii fixture");
}

// ── --version / --help ───────────────────────────────────────────────────────

#[test]
fn version_prints_the_workspace_package_version() {
    for flag in ["--version", "-V"] {
        let out = run_cli(&[flag]);
        assert_eq!(out.status.code(), Some(0), "-V/--version exits 0");
        let expected = format!("rmsh {}\n", env!("CARGO_PKG_VERSION"));
        assert_eq!(stdout_of(&out), expected, "{flag} prints the version line");
    }
}

#[test]
fn main_help_first_line_is_pinned() {
    for flag in ["--help", "-h"] {
        let out = run_cli(&[flag]);
        assert_eq!(out.status.code(), Some(0), "--help exits 0");
        let stdout = stdout_of(&out);
        let first = stdout.lines().next().expect("help has a first line");
        assert_eq!(
            first, "rmsh - mesh toolchain CLI (convert / refine / quality)",
            "main help first line is a pinned contract"
        );
        assert!(stdout.contains("Usage:"), "main help shows the usage block");
        for verb in ["convert", "refine", "quality"] {
            assert!(stdout.contains(verb), "main help lists {verb}");
        }
    }
}

#[test]
fn verb_help_first_lines_are_pinned() {
    let expected_first = [
        ("convert", "rmsh convert - convert a mesh between MSH formats"),
        ("refine", "rmsh refine - optimize tetrahedral mesh quality"),
        ("quality", "rmsh quality - audit watertightness and mesh quality"),
    ];
    for (verb, first_line) in expected_first {
        let out = run_cli(&[verb, "--help"]);
        assert_eq!(out.status.code(), Some(0), "{verb} --help exits 0");
        let stdout = stdout_of(&out);
        let first = stdout.lines().next().expect("verb help has a first line");
        assert_eq!(first, first_line, "{verb} help first line is pinned");
    }
}

// ── usage errors (exit 2) ────────────────────────────────────────────────────

#[test]
fn no_args_is_a_usage_error_on_stderr() {
    let out = run_cli(&[]);
    assert_eq!(out.status.code(), Some(2), "no args exits 2");
    let stderr = stderr_of(&out);
    assert!(stderr.contains("error:"), "stderr explains the error");
    assert!(stderr.contains("--help"), "stderr points at --help");
    assert!(stdout_of(&out).is_empty(), "no stdout chatter on usage error");
}

#[test]
fn unknown_command_is_a_usage_error() {
    let out = run_cli(&["frobnicate"]);
    assert_eq!(out.status.code(), Some(2));
    let stderr = stderr_of(&out);
    assert!(stderr.contains("unknown command 'frobnicate'"));
    assert!(stderr.contains("--help"));
}

#[test]
fn unrecognized_flag_is_a_usage_error() {
    let dir = TempDir::new("refine_bad_flag");
    let input = dir.file("in.msh");
    write_msh_v4_ascii(&input, &sliver_pair());
    let out = run_cli(&[
        "refine",
        input.to_str().unwrap(),
        dir.file("out.msh").to_str().unwrap(),
        "--wat",
    ]);
    assert_eq!(out.status.code(), Some(2), "unrecognized option exits 2");
    assert!(stderr_of(&out).contains("unrecognized option '--wat'"));
}

#[test]
fn conflicting_flavor_flags_are_a_usage_error() {
    let dir = TempDir::new("convert_conflict");
    let input = dir.file("in.msh");
    write_msh_v4_ascii(&input, &watertight_cube());
    let out = run_cli(&[
        "convert",
        input.to_str().unwrap(),
        dir.file("out.msh").to_str().unwrap(),
        "--v2",
        "--v4",
    ]);
    assert_eq!(out.status.code(), Some(2), "--v2 --v4 conflict exits 2");
    assert!(stderr_of(&out).contains("conflicting"));
}

#[test]
fn non_numeric_flag_value_is_a_usage_error() {
    let dir = TempDir::new("refine_bad_number");
    let input = dir.file("in.msh");
    write_msh_v4_ascii(&input, &sliver_pair());
    let out = run_cli(&[
        "refine",
        input.to_str().unwrap(),
        dir.file("out.msh").to_str().unwrap(),
        "--passes",
        "many",
    ]);
    assert_eq!(out.status.code(), Some(2), "bad --passes exits 2");
    let stderr = stderr_of(&out);
    assert!(stderr.contains("--passes"), "stderr names the offending flag");
    assert!(stderr.contains("many"), "stderr echoes the bad value");
}

// ── convert ──────────────────────────────────────────────────────────────────

#[test]
fn convert_v4_ascii_to_v2_binary_roundtrips_the_mesh() {
    let dir = TempDir::new("convert_v4a_v2b");
    let input = dir.file("cube_v4a.msh");
    write_msh_v4_ascii(&input, &watertight_cube());
    let output = dir.file("cube_v2b.msh");

    let out = run_cli(&[
        "convert",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "--v2",
        "--binary",
    ]);
    assert_eq!(out.status.code(), Some(0), "convert succeeds");
    assert!(
        stdout_of(&out).contains("wrote"),
        "convert prints a summary line: {:?}",
        stdout_of(&out)
    );
    assert!(output.exists(), "the output file is written");

    // Load back through rmsh-io and compare the mesh content.
    let original = load_msh_from_path(&input).expect("load fixture");
    let roundtripped = load_msh_from_path(&output).expect("load converted output");
    assert_eq!(roundtripped.node_count(), original.node_count());
    assert_eq!(roundtripped.element_count(), original.element_count());

    // Every node position survives the ASCII -> binary conversion exactly.
    for (id, node) in &original.nodes {
        let back = roundtripped
            .nodes
            .get(id)
            .unwrap_or_else(|| panic!("node {id} missing after conversion"));
        assert_eq!(back.position, node.position, "node {id} position exact");
    }
    // Every element's id, type and connectivity survive (sorted by id).
    let mut a: Vec<_> = original
        .elements
        .iter()
        .map(|e| (e.id, e.etype, e.node_ids.clone()))
        .collect();
    let mut b: Vec<_> = roundtripped
        .elements
        .iter()
        .map(|e| (e.id, e.etype, e.node_ids.clone()))
        .collect();
    a.sort_by_key(|e| e.0);
    b.sort_by_key(|e| e.0);
    assert_eq!(a, b, "elements survive conversion exactly");
}

#[test]
fn convert_defaults_to_v4_ascii_output() {
    let dir = TempDir::new("convert_default");
    let input = dir.file("in.msh");
    write_msh_v4_ascii(&input, &watertight_cube());
    let output = dir.file("out.msh");
    let out = run_cli(&["convert", input.to_str().unwrap(), output.to_str().unwrap()]);
    assert_eq!(out.status.code(), Some(0));
    // The default flavor is named on the summary line and loadable as v4.
    assert!(
        stdout_of(&out).contains("v4.1 ascii"),
        "summary names the flavor: {:?}",
        stdout_of(&out)
    );
    let data = fs::read(&output).expect("read default output");
    let text = String::from_utf8(data).expect("default output is ASCII text");
    assert!(text.contains("4.1 0 8"), "default payload is MSH v4.1 ASCII");
}

#[test]
fn convert_garbage_input_fails_cleanly() {
    let dir = TempDir::new("convert_garbage");
    let input = dir.file("garbage.msh");
    write_text(&input, "this is not a mesh\n");
    let out = run_cli(&[
        "convert",
        input.to_str().unwrap(),
        dir.file("out.msh").to_str().unwrap(),
    ]);
    assert_eq!(out.status.code(), Some(3), "garbage input exits 3");
    let stderr = stderr_of(&out);
    assert!(stderr.contains("error:"), "stderr carries the error: {stderr:?}");
    assert!(
        stderr.contains("garbage.msh"),
        "stderr names the input file: {stderr:?}"
    );
    assert!(stdout_of(&out).is_empty(), "no stdout on failure");
}

#[test]
fn convert_missing_input_fails_cleanly() {
    let dir = TempDir::new("convert_missing");
    let missing = dir.file("nope.msh");
    let out = run_cli(&[
        "convert",
        missing.to_str().unwrap(),
        dir.file("out.msh").to_str().unwrap(),
    ]);
    assert_eq!(out.status.code(), Some(3), "missing input exits 3");
    let stderr = stderr_of(&out);
    assert!(stderr.contains("error:"));
    assert!(stderr.contains("nope.msh"));
}

#[test]
fn convert_missing_arguments_is_a_usage_error() {
    let out = run_cli(&["convert", "only-one-arg"]);
    assert_eq!(out.status.code(), Some(2));
    assert!(stderr_of(&out).contains("convert requires"));
}

// ── quality ──────────────────────────────────────────────────────────────────

#[test]
fn quality_watertight_cube_exits_0_with_pinned_table() {
    let dir = TempDir::new("quality_ok");
    let input = dir.file("cube.msh");
    write_msh_v4_ascii(&input, &watertight_cube());
    let out = run_cli(&["quality", input.to_str().unwrap()]);
    assert_eq!(out.status.code(), Some(0), "watertight cube exits 0");
    let expected = "\
watertight: yes
Quality stats
  min_dihedral_deg : 0.000000
  p95_radius_edge  : inf
  sliver_fraction  : 1.000000
  tet_count        : 0
";
    assert_eq!(stdout_of(&out), expected, "quality stdout is byte-pinned");
}

#[test]
fn quality_punched_cube_exits_1_with_hole_hint() {
    let dir = TempDir::new("quality_hole");
    let input = dir.file("punched.msh");
    write_msh_v4_ascii(&input, &punched_cube());
    let out = run_cli(&["quality", input.to_str().unwrap()]);
    assert_eq!(out.status.code(), Some(1), "non-watertight exits 1");
    let expected = "\
watertight: no
defects
  boundary_loops         : 1
  non_manifold_edges     : 0
  orientation_mismatches : 0
  duplicate_triangles    : 0
  degenerate_triangles   : 0
  self_intersections     : 0
hints
  hole: cap the boundary loop of 3 edges [(5, 7), (7, 8), (8, 5)] (total length 3.414214)
Quality stats
  min_dihedral_deg : 0.000000
  p95_radius_edge  : inf
  sliver_fraction  : 1.000000
  tet_count        : 0
";
    assert_eq!(stdout_of(&out), expected, "quality stdout is byte-pinned");
}

#[test]
fn quality_missing_input_exits_3() {
    let dir = TempDir::new("quality_missing");
    let missing = dir.file("gone.msh");
    let out = run_cli(&["quality", missing.to_str().unwrap()]);
    assert_eq!(out.status.code(), Some(3));
    assert!(stderr_of(&out).contains("error:"));
}

// ── refine ───────────────────────────────────────────────────────────────────

#[test]
fn refine_improves_min_dihedral_and_prints_before_after_tables() {
    let dir = TempDir::new("refine_pair");
    let input = dir.file("sliver.msh");
    write_msh_v4_ascii(&input, &sliver_pair());
    let output = dir.file("refined.msh");

    let out = run_cli(&[
        "refine",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "--min",
        "0.3",
    ]);
    assert_eq!(out.status.code(), Some(0), "refine succeeds");
    let stdout = stdout_of(&out);
    assert!(
        stdout.contains("Quality stats (before)"),
        "before table present: {stdout:?}"
    );
    assert!(
        stdout.contains("Quality stats (after)"),
        "after table present: {stdout:?}"
    );
    // The pinned row labels appear in both tables.
    assert_eq!(stdout.matches("  min_dihedral_deg : ").count(), 2);
    assert_eq!(stdout.matches("  p95_radius_edge  : ").count(), 2);
    assert_eq!(stdout.matches("  sliver_fraction  : ").count(), 2);
    assert_eq!(stdout.matches("  tet_count        : ").count(), 2);

    // Measure: the saved output loads and its min dihedral beats the input's.
    let before_mesh = load_msh_from_path(&input).expect("load input");
    let after_mesh = load_msh_from_path(&output).expect("load refined output");
    let before = mesh_quality_stats(&before_mesh);
    let after = mesh_quality_stats(&after_mesh);
    assert!(
        after.min_dihedral_deg > before.min_dihedral_deg,
        "refine must improve the min dihedral: before={} after={}",
        before.min_dihedral_deg,
        after.min_dihedral_deg
    );
    // The optimizer's quality ledger: tets 2 -> 3 via the 2-3 swap.
    assert_eq!(before.tet_count, 2);
    assert_eq!(after.tet_count, 3);
}

#[test]
fn refine_output_flavor_flags_apply() {
    let dir = TempDir::new("refine_flavor");
    let input = dir.file("sliver.msh");
    write_msh_v4_ascii(&input, &sliver_pair());
    let output = dir.file("refined_v2b.msh");
    let out = run_cli(&[
        "refine",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "--v2",
        "--binary",
    ]);
    assert_eq!(out.status.code(), Some(0));
    let mesh = load_msh_from_path(&output).expect("v2 binary refined output loads");
    assert_eq!(mesh.element_count(), 3, "optimized mesh written");
}

#[test]
fn refine_zero_passes_is_a_noop() {
    let dir = TempDir::new("refine_zero_passes");
    let input = dir.file("sliver.msh");
    write_msh_v4_ascii(&input, &sliver_pair());
    let output = dir.file("same.msh");
    let out = run_cli(&[
        "refine",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "--passes",
        "0",
    ]);
    assert_eq!(out.status.code(), Some(0), "--passes 0 is accepted");
    let before_mesh = load_msh_from_path(&input).expect("load input");
    let after_mesh = load_msh_from_path(&output).expect("load output");
    let before = mesh_quality_stats(&before_mesh);
    let after = mesh_quality_stats(&after_mesh);
    assert_eq!(
        before.min_dihedral_deg, after.min_dihedral_deg,
        "zero passes must not touch the mesh"
    );
    assert_eq!(before.tet_count, after.tet_count);
}
