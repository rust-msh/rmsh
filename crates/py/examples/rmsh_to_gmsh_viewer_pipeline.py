from __future__ import annotations

import argparse
from pathlib import Path

import gmsh
import rmsh


def build_demo_geometry(include_torus: bool) -> None:
    # 0D points
    p0 = rmsh.model.occ.addPoint(0.0, 0.0, 0.0)
    p1 = rmsh.model.occ.addPoint(1.0, 0.0, 0.0)
    p2 = rmsh.model.occ.addPoint(1.2, 0.8, 0.0)
    p3 = rmsh.model.occ.addPoint(0.2, 1.0, 0.0)

    # 1D curves
    rmsh.model.occ.addLine(p0, p1)
    rmsh.model.occ.addSpline([p0, p2, p3])

    c_center = rmsh.model.occ.addPoint(2.5, 0.0, 0.0)
    c_start = rmsh.model.occ.addPoint(3.0, 0.0, 0.0)
    c_end = rmsh.model.occ.addPoint(2.5, 0.5, 0.0)
    rmsh.model.occ.addCircleArc(c_start, c_center, c_end)

    # 2D faces
    rmsh.model.occ.addRectangle(-1.5, 0.2, 0.0, 0.8, 0.6)
    rmsh.model.occ.addDisk(-0.2, 1.0, 0.0, 0.4, 0.2)

    # 3D solids
    rmsh.model.occ.addBox(0.0, -1.2, 0.0, 0.8, 0.6, 0.5)
    rmsh.model.occ.addCylinder(1.4, -1.2, 0.0, 0.0, 0.0, 0.8, 0.25)
    rmsh.model.occ.addSphere(2.5, -0.9, 0.3, 0.3)
    # gmsh-compatible signatures:
    #   addCone(x, y, z, dx, dy, dz, r1, r2)
    #   addTorus(x, y, z, r1, r2)
    rmsh.model.occ.addCone(3.4, -1.3, 0.0, 0.0, 0.0, 0.9, 0.35, 0.15)

    if include_torus:
        rmsh.model.occ.addTorus(4.4, -1.0, 0.3, 0.7, 0.15)


def export_rmsh_step(path: Path, strict: bool, include_torus: bool) -> None:
    rmsh.initialize()
    try:
        rmsh.clear()
        rmsh.model.add("rmsh_pipeline_demo")

        if strict:
            rmsh.option.setNumber("STEP.GmshStrict", 1)
            rmsh.option.setString("STEP.Protocol", "AP214")

        build_demo_geometry(include_torus=include_torus)
        rmsh.model.occ.synchronize()
        rmsh.write(str(path))
    finally:
        rmsh.finalize()


def rewrite_with_gmsh(src_step: Path, dst_step: Path, dst_msh: Path | None) -> None:
    gmsh.initialize()
    try:
        gmsh.model.add("gmsh_rewrite")
        gmsh.model.occ.importShapes(str(src_step))
        gmsh.model.occ.synchronize()

        # Re-emit STEP through OCC to maximize viewer compatibility.
        gmsh.write(str(dst_step))

        if dst_msh is not None:
            gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 50)
            gmsh.model.mesh.generate(2)
            gmsh.write(str(dst_msh))
    finally:
        gmsh.finalize()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export rmsh STEP, then rewrite through gmsh for viewer-friendly output."
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "out_all_dimensions_step",
        help="Output directory",
    )
    parser.add_argument(
        "--name",
        default="all_dimensions_demo_pipeline",
        help="Base file name (without extension)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Enable rmsh strict STEP mode (off by default)",
    )
    parser.add_argument(
        "--include-torus",
        action="store_true",
        help="Include torus solid in source model (off by default)",
    )
    parser.add_argument(
        "--write-msh",
        action="store_true",
        help="Also write a .msh file from the rewritten model",
    )
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    rmsh_step = args.out_dir / f"{args.name}_rmsh.step"
    gmsh_step = args.out_dir / f"{args.name}_gmsh_rewrite.step"
    gmsh_msh = args.out_dir / f"{args.name}.msh" if args.write_msh else None

    export_rmsh_step(rmsh_step, strict=args.strict, include_torus=args.include_torus)
    rewrite_with_gmsh(rmsh_step, gmsh_step, gmsh_msh)

    print(f"rmsh_step: {rmsh_step.resolve()}")
    print(f"gmsh_rewrite_step: {gmsh_step.resolve()}")
    if gmsh_msh is not None:
        print(f"gmsh_mesh: {gmsh_msh.resolve()}")


if __name__ == "__main__":
    main()
