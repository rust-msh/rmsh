from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Any

import gmsh
import rmsh


@dataclass
class OpCase:
    name: str
    op: str
    kwargs: dict[str, Any]
    required_match: bool = True


def _build_pair_gmsh() -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    a = gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.6, 1.2, 1.0)
    b = gmsh.model.occ.addSphere(1.0, 0.6, 0.5, 0.7)
    return [(3, a)], [(3, b)]


def _build_pair_rmsh() -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    a = rmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.6, 1.2, 1.0)
    b = rmsh.model.occ.addSphere(1.0, 0.6, 0.5, 0.7)
    return [(3, a)], [(3, b)]


def _run_gmsh(case: OpCase) -> tuple[list[tuple[int, int]], list[list[tuple[int, int]]]]:
    gmsh.initialize()
    try:
        gmsh.clear()
        gmsh.model.add(f"gmsh_{case.name}")
        obj, tool = _build_pair_gmsh()
        fn = getattr(gmsh.model.occ, case.op)
        out_dim_tags, out_dim_tags_map = fn(obj, tool, **case.kwargs)
        return out_dim_tags, out_dim_tags_map
    finally:
        gmsh.finalize()


def _run_rmsh(case: OpCase) -> tuple[list[tuple[int, int]], list[list[tuple[int, int]]]]:
    rmsh.initialize()
    try:
        rmsh.clear()
        rmsh.model.add(f"rmsh_{case.name}")
        obj, tool = _build_pair_rmsh()
        fn = getattr(rmsh.model.occ, case.op)
        out_dim_tags, out_dim_tags_map = fn(obj, tool, **case.kwargs)
        return out_dim_tags, out_dim_tags_map
    finally:
        rmsh.finalize()


def _shape(out: list[tuple[int, int]], mapping: list[list[tuple[int, int]]]) -> tuple[int, list[int]]:
    return len(out), [len(x) for x in mapping]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare rmsh and gmsh boolean API output contracts. "
            "By default, fragment mismatches are treated as known differences."
        )
    )
    parser.add_argument(
        "--core-only",
        action="store_true",
        help="Check only fuse/cut/intersect (skip fragment cases).",
    )
    parser.add_argument(
        "--strict-fragment",
        action="store_true",
        help="Treat fragment mismatches as failures.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    cases = [
        OpCase("fuse_default", "fuse", {}),
        OpCase("fuse_keep_inputs", "fuse", {"removeObject": False, "removeTool": False}),
        OpCase("cut_default", "cut", {}),
        OpCase("cut_keep_tool", "cut", {"removeTool": False}),
        OpCase("intersect_default", "intersect", {}),
        OpCase("intersect_tag_keep", "intersect", {"tag": 99, "removeObject": False, "removeTool": False}),
        OpCase("fragment_default", "fragment", {}, required_match=args.strict_fragment),
        OpCase(
            "fragment_keep_inputs",
            "fragment",
            {"removeObject": False, "removeTool": False},
            required_match=args.strict_fragment,
        ),
    ]

    if args.core_only:
        cases = [case for case in cases if case.op != "fragment"]

    all_required_ok = True

    for case in cases:
        g_out, g_map = _run_gmsh(case)
        r_out, r_map = _run_rmsh(case)

        g_shape = _shape(g_out, g_map)
        r_shape = _shape(r_out, r_map)

        same = g_shape == r_shape
        if case.required_match and not same:
            all_required_ok = False

        print(f"\n=== {case.name} ({case.op}) ===")
        print(f"kwargs: {case.kwargs}")
        print(f"gmsh: out={g_shape[0]} map_sizes={g_shape[1]}")
        print(f"rmsh: out={r_shape[0]} map_sizes={r_shape[1]}")
        if same:
            status = "MATCH"
        elif case.required_match:
            status = "DIFF"
        else:
            status = "KNOWN_DIFF"
        print("status:", status)

        if not same:
            print(f"gmsh out tags: {g_out}")
            print(f"rmsh out tags: {r_out}")
            print(f"gmsh map: {g_map}")
            print(f"rmsh map: {r_map}")

    print("\n" + "=" * 64)
    if all_required_ok:
        print("Boolean API contract comparison: all required checks passed.")
    else:
        print("Boolean API contract comparison: required differences detected.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
