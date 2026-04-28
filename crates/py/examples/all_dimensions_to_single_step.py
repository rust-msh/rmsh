from pathlib import Path

import rmsh


def export_all_dimensions_to_step(out_path: Path) -> None:
    rmsh.initialize()
    try:
        rmsh.clear()
        rmsh.model.add("all_dimensions_demo")

        # Optional: strict gmsh-compatible STEP profile.
        rmsh.option.setNumber("STEP.GmshStrict", 1)
        rmsh.option.setString("STEP.Protocol", "AP214")

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
        rmsh.model.occ.addTorus(4.4, -1.0, 0.3, 0.7, 0.15)

        rmsh.model.occ.synchronize()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        rmsh.write(str(out_path))
        print(f"Wrote STEP: {out_path.resolve()}")
    finally:
        rmsh.finalize()


if __name__ == "__main__":
    output = Path(__file__).resolve().parent / "out_all_dimensions_step" / "all_dimensions_demo.step"
    export_all_dimensions_to_step(output)
