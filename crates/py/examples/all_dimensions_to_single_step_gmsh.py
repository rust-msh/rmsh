from pathlib import Path

import gmsh


def export_all_dimensions_to_step_gmsh(out_path: Path) -> None:
    gmsh.initialize()
    try:
        gmsh.model.add("all_dimensions_demo_gmsh")

        # 0D points
        p0 = gmsh.model.occ.addPoint(0.0, 0.0, 0.0)
        p1 = gmsh.model.occ.addPoint(1.0, 0.0, 0.0)
        p2 = gmsh.model.occ.addPoint(1.2, 0.8, 0.0)
        p3 = gmsh.model.occ.addPoint(0.2, 1.0, 0.0)

        # 1D curves
        gmsh.model.occ.addLine(p0, p1)
        gmsh.model.occ.addSpline([p0, p2, p3])

        c_center = gmsh.model.occ.addPoint(2.5, 0.0, 0.0)
        c_start = gmsh.model.occ.addPoint(3.0, 0.0, 0.0)
        c_end = gmsh.model.occ.addPoint(2.5, 0.5, 0.0)
        gmsh.model.occ.addCircleArc(c_start, c_center, c_end)

        # 2D faces
        gmsh.model.occ.addRectangle(-1.5, 0.2, 0.0, 0.8, 0.6)
        gmsh.model.occ.addDisk(-0.2, 1.0, 0.0, 0.4, 0.2)

        # 3D solids
        gmsh.model.occ.addBox(0.0, -1.2, 0.0, 0.8, 0.6, 0.5)
        gmsh.model.occ.addCylinder(1.4, -1.2, 0.0, 0.0, 0.0, 0.8, 0.25)
        gmsh.model.occ.addSphere(2.5, -0.9, 0.3, 0.3)
        gmsh.model.occ.addCone(3.4, -1.3, 0.0, 0.0, 0.0, 0.9, 0.35, 0.15)
        gmsh.model.occ.addTorus(4.4, -1.0, 0.3, 0.7, 0.15)

        gmsh.model.occ.synchronize()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(out_path))
        print(f"Wrote STEP: {out_path.resolve()}")
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    output = Path(__file__).resolve().parent / "out_all_dimensions_step_gmsh" / "all_dimensions_demo_gmsh.step"
    export_all_dimensions_to_step_gmsh(output)
