from __future__ import annotations

import re
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "out_point_curve_surface_solid_step"
RMSH_STEP = OUT_DIR / "rmsh_baseline_strict.step"
GMSH_STEP = OUT_DIR / "gmsh_baseline.step"

TOKENS = [
    "ADVANCED_FACE",
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
    "ELLIPSE",
    "B_SPLINE_CURVE_WITH_KNOTS",
    "LINE",
    "VERTEX_POINT",
    "GEOMETRIC_CURVE_SET",
]


def run_script(script_name: str) -> None:
    script = ROOT / script_name
    cmd = [sys.executable, str(script)]
    subprocess.run(cmd, check=True)


def count_tokens(step_path: Path) -> dict[str, int]:
    text = step_path.read_text(encoding="utf-8", errors="ignore")
    return {token: text.count(token) for token in TOKENS}


def count_entities_exact(step_path: Path) -> dict[str, int]:
    # Parse STEP entities as lines like: #123=ENTITY_TYPE(...);
    # This avoids substring collisions, e.g. PLANE vs PLANE_ANGLE_*.
    entity_re = re.compile(r"^\s*#\d+\s*=\s*([A-Z0-9_]+)\s*\(")
    counts = {token: 0 for token in TOKENS}
    for line in step_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = entity_re.match(line)
        if not m:
            continue
        t = m.group(1)
        if t in counts:
            counts[t] += 1
    return counts


def print_table(rmsh_counts: dict[str, int], gmsh_counts: dict[str, int]) -> None:
    print("token                           rmsh   gmsh   diff(rmsh-gmsh)")
    print("-------------------------------------------------------------")
    for token in TOKENS:
        r = rmsh_counts.get(token, 0)
        g = gmsh_counts.get(token, 0)
        d = r - g
        print(f"{token:30} {r:5d} {g:6d} {d:16d}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    run_script("export_point_curve_surface_solid_step_rmsh_gmsh_strict.py")
    run_script("gmsh_step_baseline.py")

    rmsh_legacy = count_tokens(RMSH_STEP)
    gmsh_legacy = count_tokens(GMSH_STEP)
    rmsh_exact = count_entities_exact(RMSH_STEP)
    gmsh_exact = count_entities_exact(GMSH_STEP)

    print(f"rmsh: {RMSH_STEP}")
    print(f"gmsh: {GMSH_STEP}")
    print("\n[Exact STEP entity counts]")
    print_table(rmsh_exact, gmsh_exact)

    print("\n[Legacy substring counts]")
    print_table(rmsh_legacy, gmsh_legacy)

    print(
        "\nNote: mismatch verdict uses exact STEP entity counts; "
        "substring counts are informative only."
    )

    mismatches = [t for t in TOKENS if rmsh_exact[t] != gmsh_exact[t]]
    if mismatches:
        print("\nMismatched entities (exact):")
        for t in mismatches:
            print(f"- {t}")
    else:
        print("\nAll tracked STEP entity counts match.")


if __name__ == "__main__":
    main()
