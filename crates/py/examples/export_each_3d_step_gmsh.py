from pathlib import Path

import gmsh


def export_one(name: str, builder, out_dir: Path) -> Path:
    out_path = out_dir / f"{name}.step"
    gmsh.initialize()
    try:
        gmsh.model.add(name)
        builder()
        gmsh.model.occ.synchronize()
        gmsh.write(str(out_path))
    finally:
        gmsh.finalize()
    return out_path


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "out_each_3d_step_gmsh"
    out_dir.mkdir(parents=True, exist_ok=True)

    exports = [
        (
            "box",
            lambda: gmsh.model.occ.addBox(0.0, 0.0, 0.0, 2.0, 1.2, 1.0),
        ),
        (
            "sphere",
            lambda: gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0),
        ),
        (
            "cylinder",
            lambda: gmsh.model.occ.addCylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.5),
        ),
        (
            "cone",
            lambda: gmsh.model.occ.addCone(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.8, 0.2),
        ),
        (
            "torus",
            lambda: gmsh.model.occ.addTorus(0.0, 0.0, 0.0, 1.0, 0.25),
        ),
    ]

    for name, builder in exports:
        out = export_one(name, builder, out_dir)
        print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
