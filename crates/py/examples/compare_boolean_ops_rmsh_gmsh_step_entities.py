from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RMSH_DIR = ROOT / "out_boolean_ops_step"
GMSH_DIR = ROOT / "out_boolean_ops_step_gmsh"

PAIRS = [
    ("fuse", RMSH_DIR / "bool_fuse.step", GMSH_DIR / "gmsh_bool_fuse.step"),
    ("cut", RMSH_DIR / "bool_cut.step", GMSH_DIR / "gmsh_bool_cut.step"),
    ("intersect", RMSH_DIR / "bool_intersect.step", GMSH_DIR / "gmsh_bool_intersect.step"),
    ("fragment", RMSH_DIR / "bool_fragment.step", GMSH_DIR / "gmsh_bool_fragment.step"),
]

TOKENS = [
    "ADVANCED_FACE",
    "CLOSED_SHELL",
    "MANIFOLD_SOLID_BREP",
    "CONICAL_SURFACE",
    "CYLINDRICAL_SURFACE",
    "SPHERICAL_SURFACE",
    "TOROIDAL_SURFACE",
    "PLANE",
    "EDGE_CURVE",
    "ORIENTED_EDGE",
    "SEAM_CURVE",
    "SURFACE_CURVE",
    "CIRCLE",
    "LINE",
    "VERTEX_POINT",
]

ENTITY_RE = re.compile(r"^\s*#\d+\s*=\s*([A-Z0-9_]+)\s*\(")


def run_script(script_name: str) -> None:
    script = ROOT / script_name
    subprocess.run([sys.executable, str(script)], check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare rmsh/gmsh STEP entity counts for boolean examples. "
            "By default, mismatches fail with exit code 1."
        )
    )
    parser.add_argument(
        "--allow-differences",
        action="store_true",
        help="Do not fail process when entity mismatches are found.",
    )
    parser.add_argument(
        "--ops",
        nargs="+",
        choices=[name for name, _, _ in PAIRS],
        help="Compare only selected operations (default: all).",
    )
    return parser.parse_args()


def count_entities_exact(step_path: Path) -> dict[str, int]:
    counts = {token: 0 for token in TOKENS}
    for line in step_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = ENTITY_RE.match(line)
        if not m:
            continue
        token = m.group(1)
        if token in counts:
            counts[token] += 1
    return counts


def count_total_entities(step_path: Path) -> int:
    total = 0
    for line in step_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if ENTITY_RE.match(line):
            total += 1
    return total


def print_table(rmsh_counts: dict[str, int], gmsh_counts: dict[str, int]) -> None:
    print("entity                          rmsh   gmsh   diff(rmsh-gmsh)")
    print("--------------------------------------------------------------")
    for token in TOKENS:
        r = rmsh_counts.get(token, 0)
        g = gmsh_counts.get(token, 0)
        print(f"{token:30} {r:5d} {g:6d} {r - g:16d}")


def main() -> None:
    args = parse_args()

    RMSH_DIR.mkdir(parents=True, exist_ok=True)
    GMSH_DIR.mkdir(parents=True, exist_ok=True)

    # Regenerate baseline files so comparison reflects current backend behavior.
    run_script("boolean_ops_to_step.py")
    run_script("boolean_ops_to_step_gmsh.py")

    selected_ops = set(args.ops) if args.ops else {name for name, _, _ in PAIRS}
    selected_pairs = [pair for pair in PAIRS if pair[0] in selected_ops]
    if not selected_pairs:
        raise SystemExit("No operations selected")

    any_mismatch = False

    for op, rmsh_step, gmsh_step in selected_pairs:
        print(f"\n=== {op.upper()} ===")
        print(f"rmsh: {rmsh_step}")
        print(f"gmsh: {gmsh_step}")

        rmsh_counts = count_entities_exact(rmsh_step)
        gmsh_counts = count_entities_exact(gmsh_step)
        print_table(rmsh_counts, gmsh_counts)

        rmsh_total = count_total_entities(rmsh_step)
        gmsh_total = count_total_entities(gmsh_step)
        print(f"TOTAL_ENTITIES                 {rmsh_total:5d} {gmsh_total:6d} {rmsh_total - gmsh_total:16d}")

        mismatches = [t for t in TOKENS if rmsh_counts[t] != gmsh_counts[t]]
        if mismatches or rmsh_total != gmsh_total:
            any_mismatch = True
            print("Mismatches:")
            if mismatches:
                for t in mismatches:
                    print(f"- {t}")
            if rmsh_total != gmsh_total:
                print("- TOTAL_ENTITIES")
        else:
            print("No mismatches in tracked entity counts.")

    print("\n" + "=" * 62)
    if any_mismatch:
        print("Comparison finished with differences.")
        if not args.allow_differences:
            raise SystemExit(1)
    else:
        print("Comparison finished: tracked entity counts all match.")


if __name__ == "__main__":
    main()
