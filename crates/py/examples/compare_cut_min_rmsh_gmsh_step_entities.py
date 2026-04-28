from __future__ import annotations

import re
from pathlib import Path

import gmsh
import rmsh

ROOT = Path(__file__).resolve().parent
OUT_RMSH = ROOT / "out_cut_min_step" / "rmsh_cut_min.step"
OUT_GMSH = ROOT / "out_cut_min_step" / "gmsh_cut_min.step"

TOKENS = [
    "ADVANCED_FACE",
    "CLOSED_SHELL",
    "MANIFOLD_SOLID_BREP",
    "SPHERICAL_SURFACE",
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


def export_rmsh_cut_min() -> dict:
    OUT_RMSH.parent.mkdir(parents=True, exist_ok=True)

    rmsh.initialize()
    try:
        rmsh.clear()
        rmsh.model.add("rmsh_cut_min")

        box = rmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        sphere = rmsh.model.occ.addSphere(0.5, 0.5, 0.5, 0.45)
        out_dim_tags, _ = rmsh.model.occ.cut([(3, box)], [(3, sphere)])
        if not out_dim_tags:
            raise RuntimeError("rmsh cut returned no output")

        out_tag = out_dim_tags[0][1]
        geom_stats = rmsh._rmsh.model_occ_debug_shape_geom(out_tag)
        print(f"[rmsh cut-min] pre-sync geom tag={out_tag}: {geom_stats}")

        rmsh.model.occ.synchronize()
        rmsh.write(str(OUT_RMSH))
    finally:
        rmsh.finalize()

    print(f"[rmsh cut-min] wrote {OUT_RMSH}")
    return geom_stats


def export_gmsh_cut_min() -> None:
    OUT_GMSH.parent.mkdir(parents=True, exist_ok=True)

    gmsh.initialize()
    try:
        gmsh.clear()
        gmsh.model.add("gmsh_cut_min")

        box = gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        sphere = gmsh.model.occ.addSphere(0.5, 0.5, 0.5, 0.45)
        gmsh.model.occ.cut([(3, box)], [(3, sphere)])

        gmsh.model.occ.synchronize()
        gmsh.write(str(OUT_GMSH))
    finally:
        gmsh.finalize()

    print(f"[gmsh cut-min] wrote {OUT_GMSH}")


def main() -> None:
    _ = export_rmsh_cut_min()
    export_gmsh_cut_min()

    rmsh_counts = count_entities_exact(OUT_RMSH)
    gmsh_counts = count_entities_exact(OUT_GMSH)

    print("\n=== CUT_MIN ===")
    print(f"rmsh: {OUT_RMSH}")
    print(f"gmsh: {OUT_GMSH}")
    print_table(rmsh_counts, gmsh_counts)

    rmsh_total = count_total_entities(OUT_RMSH)
    gmsh_total = count_total_entities(OUT_GMSH)
    print(f"TOTAL_ENTITIES                 {rmsh_total:5d} {gmsh_total:6d} {rmsh_total - gmsh_total:16d}")


if __name__ == "__main__":
    main()
