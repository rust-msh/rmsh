from pathlib import Path

import gmsh


def _build_overlap_pair() -> tuple[int, int]:
    a = gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.6, 1.2, 1.0)
    b = gmsh.model.occ.addSphere(1.0, 0.6, 0.5, 0.7)
    return a, b


def _export_case(out_dir: Path, name: str, op: str) -> None:
    out_file = out_dir / f"{name}.step"

    gmsh.initialize()
    try:
        gmsh.clear()
        gmsh.model.add(name)

        a, b = _build_overlap_pair()
        obj = [(3, a)]
        tool = [(3, b)]

        if op == "fuse":
            gmsh.model.occ.fuse(obj, tool)
        elif op == "cut":
            gmsh.model.occ.cut(obj, tool)
        elif op == "intersect":
            gmsh.model.occ.intersect(obj, tool)
        elif op == "fragment":
            gmsh.clear()
            gmsh.model.add(name)
            p = gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.3, 1.0, 1.0)
            q = gmsh.model.occ.addBox(0.7, 0.2, 0.0, 1.3, 1.0, 1.0)
            gmsh.model.occ.fragment([(3, p)], [(3, q)])
        else:
            raise ValueError(f"Unknown operation: {op}")

        gmsh.model.occ.synchronize()
        gmsh.write(str(out_file))
    finally:
        gmsh.finalize()

    print(f"[{op}] wrote {out_file}")


def export_boolean_ops_to_step_gmsh(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    _export_case(out_dir, "gmsh_bool_fuse", "fuse")
    _export_case(out_dir, "gmsh_bool_cut", "cut")
    _export_case(out_dir, "gmsh_bool_intersect", "intersect")
    _export_case(out_dir, "gmsh_bool_fragment", "fragment")


if __name__ == "__main__":
    output_dir = Path(__file__).resolve().parent / "out_boolean_ops_step_gmsh"
    export_boolean_ops_to_step_gmsh(output_dir)
