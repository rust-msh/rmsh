from pathlib import Path

import rmsh


def export_point_curve_surface_solid_one_step(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "point_curve_surface_solid.step"

    rmsh.initialize()
    try:
        rmsh.clear()

        # Point: represented as a tiny sphere centered at the point.
        rmsh.model.occ.addSphere(1.0, 2.0, 3.0, 1e-3)

        # Curve: represented as a very thin cylinder along the segment.
        rmsh.model.occ.addCylinder(4.0, 0.0, 0.0, 2.0, 1.0, 0.5, 1e-3)

        # Surface: analytic rectangle face.
        rmsh.model.occ.addRectangle(0.0, 5.0, 0.0, 2.0, 1.0)

        # Solids: box + sphere + cone + torus.
        rmsh.model.occ.addBox(4.0, 5.0, 0.0, 1.0, 2.0, 3.0)
        rmsh.model.occ.addSphere(8.0, 5.8, 1.0, 0.9)
        rmsh.model.occ.addCone(10.5, 5.0, 0.0, 0.0, 0.0, 2.5, 0.9, 0.2)
        rmsh.model.occ.addTorus(14.0, 6.0, 1.2, 0.0, 0.0, 1.0, 1.0, 0.3)

        # One-shot export: all shapes in current CAD scene are written as one STEP.
        rmsh.write(str(out_file))

    finally:
        rmsh.finalize()

    print(f"Exported: {out_file}")


if __name__ == "__main__":
    output_dir = Path(__file__).resolve().parent / "out_point_curve_surface_solid_step"
    export_point_curve_surface_solid_one_step(output_dir)
