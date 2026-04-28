from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import gmsh
import rmsh

ROOT = Path(__file__).resolve().parent
OUT_RMSH = ROOT / "out_boolean_ops_step" / "bool_intersect.step"
OUT_GMSH = ROOT / "out_boolean_ops_step_gmsh" / "gmsh_bool_intersect.step"

ENTITY_RE = re.compile(r"^\s*#(\d+)\s*=\s*([A-Z0-9_]+)\s*\((.*)\);\s*$")
ID_RE = re.compile(r"#(\d+)")


def parse(path: Path):
    ents = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = ENTITY_RE.match(line)
        if m:
            ents[int(m.group(1))] = (m.group(2), m.group(3))
    return ents


def step_entity_counts(ents) -> Counter:
    return Counter(t for t, _ in ents.values())


def resolve_edge_basis(ents, curve_id: int) -> tuple[int | None, str]:
    ctyp, carg = ents.get(curve_id, ("?", ""))
    if ctyp in ("SURFACE_CURVE", "SEAM_CURVE"):
        ids = [int(x) for x in ID_RE.findall(carg)]
        if ids:
            basis_id = ids[0]
            basis_typ = ents.get(basis_id, ("?", ""))[0]
            return basis_id, basis_typ
        return None, "?"
    return curve_id, ctyp


def analyze_edge_curve_sources(ents) -> dict:
    edge_curve_wrapper_types = Counter()
    edge_curve_basis_types = Counter()
    basis_circle_ids = set()
    edge_total = 0

    for _eid, (typ, args) in ents.items():
        if typ != "EDGE_CURVE":
            continue
        ids = [int(x) for x in ID_RE.findall(args)]
        if len(ids) < 3:
            continue
        edge_total += 1
        curve_id = ids[2]
        wrapper_typ = ents.get(curve_id, ("?", ""))[0]
        edge_curve_wrapper_types[wrapper_typ] += 1

        basis_id, basis_typ = resolve_edge_basis(ents, curve_id)
        edge_curve_basis_types[basis_typ] += 1
        if basis_typ == "CIRCLE" and basis_id is not None:
            basis_circle_ids.add(basis_id)

    counts = step_entity_counts(ents)
    total_circles = counts.get("CIRCLE", 0)
    referenced_basis_circles = len(basis_circle_ids)

    return {
        "edge_total": edge_total,
        "edge_curve_wrapper_types": dict(edge_curve_wrapper_types),
        "edge_curve_basis_types": dict(edge_curve_basis_types),
        "total_circles": total_circles,
        "referenced_basis_circles": referenced_basis_circles,
        "non_basis_circles": max(total_circles - referenced_basis_circles, 0),
        "counts": counts,
    }


def export_rmsh_intersect(path: Path) -> dict:
    rmsh.initialize()
    try:
        rmsh.clear()
        rmsh.model.add("rmsh_intersect_probe")
        a = rmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.6, 1.2, 1.0)
        b = rmsh.model.occ.addSphere(1.0, 0.6, 0.5, 0.7)
        out_dim_tags, _ = rmsh.model.occ.intersect([(3, a)], [(3, b)])
        if not out_dim_tags:
            raise RuntimeError("rmsh intersect returned no output")
        out_tag = out_dim_tags[0][1]
        pre_sync = rmsh._rmsh.model_occ_debug_shape_geom(out_tag)
        rmsh.model.occ.synchronize()
        path.parent.mkdir(parents=True, exist_ok=True)
        rmsh.write(str(path))
        return pre_sync
    finally:
        rmsh.finalize()


def export_gmsh_intersect(path: Path) -> None:
    gmsh.initialize()
    try:
        gmsh.clear()
        gmsh.model.add("gmsh_intersect_probe")
        a = gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.6, 1.2, 1.0)
        b = gmsh.model.occ.addSphere(1.0, 0.6, 0.5, 0.7)
        gmsh.model.occ.intersect([(3, a)], [(3, b)])
        gmsh.model.occ.synchronize()
        path.parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()


def dump(label: str, step: Path, pre_sync: dict | None = None) -> None:
    ents = parse(step)
    info = analyze_edge_curve_sources(ents)

    print(f"== {label} ==")
    if pre_sync is not None:
        print("pre-sync curve3_kinds:", pre_sync.get("curve3_kinds", {}))
        print("pre-sync surface3_kinds:", pre_sync.get("surface3_kinds", {}))
        print("pre-sync edge_curve_some:", pre_sync.get("edge_curve_some"))
        print("pre-sync edges:", pre_sync.get("edges"))
    print("step EDGE_CURVE count:", info["edge_total"])
    print("step wrapper types for EDGE_CURVE basis ref:", info["edge_curve_wrapper_types"])
    print("step resolved basis types used by EDGE_CURVE:", info["edge_curve_basis_types"])
    print("step total CIRCLE entities:", info["total_circles"])
    print("step CIRCLE used as EDGE_CURVE basis:", info["referenced_basis_circles"])
    print("step non-basis CIRCLE entities:", info["non_basis_circles"])
    print("step total SURFACE_CURVE:", info["counts"].get("SURFACE_CURVE", 0))
    print("step total SEAM_CURVE:", info["counts"].get("SEAM_CURVE", 0))
    print()


def main() -> None:
    rmsh_pre_sync = export_rmsh_intersect(OUT_RMSH)
    export_gmsh_intersect(OUT_GMSH)
    dump("rmsh_intersect", OUT_RMSH, rmsh_pre_sync)
    dump("gmsh_intersect", OUT_GMSH, None)


if __name__ == "__main__":
    main()
