from pathlib import Path

import rmsh


def export_rmsh_step(out_path: Path) -> None:
    rmsh.initialize()
    try:
        rmsh.clear()
        rmsh.model.add("rmsh_baseline")

        # Strict gmsh-compatible STEP writer profile.
        rmsh.option.setNumber("STEP.GmshStrict", 1)
        rmsh.option.setString("STEP.Protocol", "AP214")

        # Solids
        rmsh.model.occ.addBox(-1.2, -0.6, 0.0, 0.5, 0.2, 0.2)
        rmsh.model.occ.addCylinder(-0.3, -0.5, 0.0, 0.0, 0.0, 0.4, 0.1)
        rmsh.model.occ.addSphere(0.1, -0.4, 0.0, 0.223606798)
        # gmsh-compatible signatures:
        #   addTorus(x, y, z, r1, r2)
        #   addCone(x, y, z, dx, dy, dz, r1, r2)
        rmsh.model.occ.addTorus(0.8, -0.5, 0.0, 0.25, 0.05)
        rmsh.model.occ.addCone(1.3, -0.5, 0.0, 0.0, 0.0, -0.3, 0.1, 0.2)

        # Surfaces
        rmsh.model.occ.addRectangle(-1.2, 0.2, 0.0, 0.4, 0.3)
        rmsh.model.occ.addDisk(-0.5, 0.5, 0.0, 0.2, 0.2)
        rmsh.model.occ.addDisk(0.1, 0.3, 0.0, 0.2, 0.1)

        # Standalone curves
        p1 = rmsh.model.occ.addPoint(0.9, 0.5, 0.0)
        p2 = rmsh.model.occ.addPoint(0.3, 0.6, 0.0)
        rmsh.model.occ.addLine(p1, p2)

        p3 = rmsh.model.occ.addPoint(1.1, 0.2, 0.0)
        p4 = rmsh.model.occ.addPoint(1.8, -0.2, 0.0)
        p5 = rmsh.model.occ.addPoint(1.6, 0.6, 0.0)
        rmsh.model.occ.addSpline([p3, p4, p5])

        p6 = rmsh.model.occ.addPoint(1.2, 0.6, 0.0)
        p7 = rmsh.model.occ.addPoint(1.7, 0.6, 0.0)
        p8 = rmsh.model.occ.addPoint(1.2, 1.1, 0.0)
        rmsh.model.occ.addCircleArc(p7, p6, p8)

        rmsh.model.occ.synchronize()
        rmsh.write(str(out_path))
        print(f"Wrote: {out_path.resolve()}")
    finally:
        rmsh.finalize()


if __name__ == "__main__":
    out = (
        Path(__file__).resolve().parent
        / "out_point_curve_surface_solid_step"
        / "rmsh_baseline_strict.step"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    export_rmsh_step(out)
