from pathlib import Path

import rmsh


def export_surface_and_solids_one_step(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "surface_and_solids.step"

    rmsh.initialize()
    try:
        rmsh.clear()

        # Surface: analytic rectangle face.
        rmsh.model.occ.addRectangle(0.0, 5.0, 0.0, 2.0, 1.0)

        # Solids: box + sphere + cone + torus.
        rmsh.model.occ.addBox(4.0, 5.0, 0.0, 1.0, 2.0, 3.0)
        rmsh.model.occ.addSphere(8.0, 5.8, 1.0, 0.9)
        # gmsh-compatible signatures:
        #   addCone(x, y, z, dx, dy, dz, r1, r2)
        #   addTorus(x, y, z, r1, r2)
        rmsh.model.occ.addCone(10.5, 5.0, 0.0, 0.0, 0.0, 2.5, 0.9, 0.2)
        rmsh.model.occ.addTorus(14.0, 6.0, 1.2, 1.0, 0.3)

        # One-shot export: all shapes in current CAD scene are written as one STEP.
        rmsh.write(str(out_file))

    finally:
        rmsh.finalize()

    print(f"Exported: {out_file}")


if __name__ == "__main__":
    output_dir = Path(__file__).resolve().parent / "out_point_curve_surface_solid_step"
    export_surface_and_solids_one_step(output_dir)
