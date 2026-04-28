from pathlib import Path

import rmsh


def export_one(name: str, build_fn, out_dir: Path) -> Path:
    out_path = out_dir / f"{name}.step"
    rmsh.initialize()
    try:
        rmsh.clear()
        build_fn()
        rmsh.write(str(out_path))
    finally:
        rmsh.finalize()
    return out_path


def main() -> None:
    out_dir = Path(__file__).resolve().parent / "out_each_3d_step"
    out_dir.mkdir(parents=True, exist_ok=True)

    # gmsh-compatible signatures used below:
    #   addCone(x, y, z, dx, dy, dz, r1, r2)
    #   addTorus(x, y, z, r1, r2)
    exports = [
        (
            "box",
            lambda: rmsh.model.occ.addBox(0.0, 0.0, 0.0, 2.0, 1.2, 1.0),
        ),
        (
            "sphere",
            lambda: rmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0),
        ),
        (
            "cylinder",
            lambda: rmsh.model.occ.addCylinder(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.5),
        ),
        (
            "cone",
            lambda: rmsh.model.occ.addCone(0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.8, 0.2),
        ),
        (
            "torus",
            lambda: rmsh.model.occ.addTorus(0.0, 0.0, 0.0, 1.0, 0.25),
        ),
    ]

    for name, builder in exports:
        out = export_one(name, builder, out_dir)
        print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
