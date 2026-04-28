from pathlib import Path

import rmsh


def _build_overlap_pair() -> tuple[int, int]:
    """Create two overlapping solids for boolean operations."""
    a = rmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.6, 1.2, 1.0)
    b = rmsh.model.occ.addSphere(1.0, 0.6, 0.5, 0.7)
    return a, b


def _export_case(out_dir: Path, name: str, op: str) -> None:
    out_file = out_dir / f"{name}.step"

    rmsh.initialize()
    try:
        rmsh.clear()
        rmsh.model.add(name)

        a, b = _build_overlap_pair()
        obj = [(3, a)]
        tool = [(3, b)]

        if op == "fuse":
            out_dim_tags, _ = rmsh.model.occ.fuse(obj, tool)
            if out_dim_tags:
                out_tag = out_dim_tags[0][1]
                stats = rmsh._rmsh.model_occ_debug_shape_geom(out_tag)
                print(f"[fuse] pre-sync geom tag={out_tag}: {stats}")
        elif op == "cut":
            out_dim_tags, _ = rmsh.model.occ.cut(obj, tool)
            if out_dim_tags:
                out_tag = out_dim_tags[0][1]
                stats = rmsh._rmsh.model_occ_debug_shape_geom(out_tag)
                print(f"[cut] pre-sync geom tag={out_tag}: {stats}")
        elif op == "intersect":
            out_dim_tags, _ = rmsh.model.occ.intersect(obj, tool)
            if out_dim_tags:
                out_tag = out_dim_tags[0][1]
                stats = rmsh._rmsh.model_occ_debug_shape_geom(out_tag)
                print(f"[intersect] pre-sync geom tag={out_tag}: {stats}")
        elif op == "fragment":
            # Fragment-like example: split overlap region between two boxes.
            rmsh.clear()
            rmsh.model.add(name)
            p = rmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.3, 1.0, 1.0)
            q = rmsh.model.occ.addBox(0.7, 0.2, 0.0, 1.3, 1.0, 1.0)
            rmsh.model.occ.fragment([(3, p)], [(3, q)])
        else:
            raise ValueError(f"Unknown operation: {op}")

        rmsh.model.occ.synchronize()
        rmsh.write(str(out_file))
    finally:
        rmsh.finalize()

    print(f"[{op}] wrote {out_file}")


def export_boolean_ops_to_step(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    _export_case(out_dir, "bool_fuse", "fuse")
    _export_case(out_dir, "bool_cut", "cut")
    _export_case(out_dir, "bool_intersect", "intersect")
    _export_case(out_dir, "bool_fragment", "fragment")


if __name__ == "__main__":
    output_dir = Path(__file__).resolve().parent / "out_boolean_ops_step"
    export_boolean_ops_to_step(output_dir)
