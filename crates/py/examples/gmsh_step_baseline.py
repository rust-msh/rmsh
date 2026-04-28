from pathlib import Path

import gmsh


def export_gmsh_step(out_path: Path) -> None:
    gmsh.initialize()
    try:
        gmsh.model.add("gmsh_baseline")

        # Solids
        gmsh.model.occ.addBox(-1.2, -0.6, 0.0, 0.5, 0.2, 0.2)
        gmsh.model.occ.addCylinder(-0.3, -0.5, 0.0, 0.0, 0.0, 0.4, 0.1)
        gmsh.model.occ.addSphere(0.1, -0.4, 0.0, 0.223606798)
        gmsh.model.occ.addTorus(0.8, -0.5, 0.0, 0.25, 0.05)
        gmsh.model.occ.addCone(1.3, -0.5, 0.0, 0.0, 0.0, -0.3, 0.1, 0.2)

        # Surfaces
        gmsh.model.occ.addRectangle(-1.2, 0.2, 0.0, 0.4, 0.3)
        gmsh.model.occ.addDisk(-0.5, 0.5, 0.0, 0.2, 0.2)
        gmsh.model.occ.addDisk(0.1, 0.3, 0.0, 0.2, 0.1)

        # Standalone curves
        p1 = gmsh.model.occ.addPoint(0.9, 0.5, 0.0)
        p2 = gmsh.model.occ.addPoint(0.3, 0.6, 0.0)
        gmsh.model.occ.addLine(p1, p2)

        p3 = gmsh.model.occ.addPoint(1.1, 0.2, 0.0)
        p4 = gmsh.model.occ.addPoint(1.8, -0.2, 0.0)
        p5 = gmsh.model.occ.addPoint(1.6, 0.6, 0.0)
        gmsh.model.occ.addSpline([p3, p4, p5])

        p6 = gmsh.model.occ.addPoint(1.2, 0.6, 0.0)
        p7 = gmsh.model.occ.addPoint(1.7, 0.6, 0.0)
        p8 = gmsh.model.occ.addPoint(1.2, 1.1, 0.0)
        gmsh.model.occ.addCircleArc(p7, p6, p8)

        gmsh.model.occ.synchronize()
        gmsh.write(str(out_path))
        print(f"Wrote: {out_path.resolve()}")
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    out = Path(__file__).resolve().parent / "out_point_curve_surface_solid_step" / "gmsh_baseline.step"
    out.parent.mkdir(parents=True, exist_ok=True)
    export_gmsh_step(out)
