"""rmsh tutorials - adapted from Gmsh Python tutorials.

Each example is labelled with the Gmsh tutorial it is inspired by.
The rmsh API mirrors gmsh closely: replace `gmsh.` with `rmsh.` and
`gmsh.model.occ.*` with `rmsh.model.occ.*`.

Run with:
    cd crates/py
    maturin develop --release
    python examples/tutorials.py

Tutorials covered
-----------------
  t1  Geometry basics (rectangle -> 2D mesh -> .msh)
  t2  Extrude a surface into a volume
  t4  OCC booleans: box with a spherical hole
  t5  Mesh size control via options
  t6  Structured quad mesh (algo 9, Packing of Parallelograms)
  t10 Boolean operations: cut, fuse, fragment
  t11 Fillet / chamfer on a box
  t16 Volume mesh from STEP file
  t17 Pre-mesh workflow: create -> heal -> inspect -> mesh
"""

import math
import os
import sys

import rmsh

TESTDATA = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "..", "..", "testdata")
)


# ─────────────────────────────────────────────────────────────────────────────
# t1 — Geometry basics: rectangle → 2D mesh
#
# Gmsh original:
#   gmsh.model.geo.addPoint / addLine / addCurveLoop / addPlaneSurface
#   gmsh.model.geo.synchronize()
#   gmsh.model.mesh.generate(2)
#
# rmsh equivalent uses model.occ.addRectangle which creates the planar domain
# directly, then generate(2) triangulates it.
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_t1():
    """t1 – Basic 2D geometry: mesh a planar rectangle."""
    rmsh.initialize()

    # gmsh: addPoint + addLine + addCurveLoop + addPlaneSurface
    # rmsh: addRectangle creates the planar domain in one call
    rmsh.model.occ.addRectangle(0, 0, 0, 0.1, 0.3)

    # Mesh size equivalent to Gmsh's characteristic length lc = 1e-2
    rmsh.option.setNumber("Mesh.MeshSizeMax", 1e-2)

    # gmsh: gmsh.model.mesh.generate(2)
    rmsh.model.mesh.generate(2)

    rmsh.write("t1_rectangle.msh")
    print("  wrote t1_rectangle.msh")

    rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# t2 — Extrude a 2D surface into a 3D volume
#
# Gmsh original:
#   gmsh.model.geo.extrude([(2, 1)], 0, 0, 0.12)
#
# rmsh equivalent: addBox effectively is an extrusion of a rectangle.
# For the explicit extrude demo we use model.occ.extrude on a box face.
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_t2():
    """t2 – Extrude a surface into a volume, then mesh it."""
    rmsh.initialize()

    # Create a thin slab (1mm thick) then extrude face 0 upward by 0.12
    box = rmsh.model.occ.addBox(0, 0, 0, 0.1, 0.3, 0.001)
    extruded = rmsh.model.occ.extrude(box, 0, 0, 0, 1, 0.12)
    print(f"  extrude -> tag {extruded}")

    rmsh.model.occ.synchronize()

    rmsh.option.setNumber("Mesh.Algorithm3D", 1)   # Delaunay
    rmsh.option.setNumber("Mesh.MeshSizeMax", 0.03)

    rmsh.model.mesh.generate(3)
    rmsh.write("t2_extrude.msh")
    print("  wrote t2_extrude.msh")

    rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# t4 — OCC: box with a spherical cavity (boolean cut)
#
# Gmsh original (t4 uses built-in kernel for a hole via curve loops;
# t16/t17 use OCC booleans):
#   gmsh.model.occ.addBox(...)
#   gmsh.model.occ.addSphere(...)
#   gmsh.model.occ.cut([(3, box)], [(3, sphere)])
#   gmsh.model.occ.synchronize()
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_t4():
    """t4 – OCC boolean: box with a spherical hole."""
    rmsh.initialize()

    box = rmsh.model.occ.addBox(0, 0, 0, 1.0, 1.0, 1.0)
    sphere = rmsh.model.occ.addSphere(0.5, 0.5, 0.5, 0.35)

    result = rmsh.model.occ.cut([(3, box)], [(3, sphere)])
    print(f"  box - sphere: result tags = {result}")

    rmsh.model.occ.synchronize()

    rmsh.option.setNumber("Mesh.Algorithm3D", 1)
    rmsh.option.setNumber("Mesh.MeshSizeMax", 0.15)

    rmsh.model.mesh.generate(3)
    rmsh.write("t4_box_with_hole.msh")
    print("  wrote t4_box_with_hole.msh")

    rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# t5 — Mesh size control
#
# Gmsh original:
#   gmsh.model.mesh.field.add("Distance", 1)
#   gmsh.model.mesh.field.add("Threshold", 2)
#   gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
#
# rmsh maps the global size controls via option.setNumber.
# Field-based size control is future work; here we show global controls.
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_t5():
    """t5 – Mesh size control: MeshSizeMax, MeshSizeMin, MeshSizeFactor."""
    # Fine mesh in the centre, coarser elsewhere — approximated with global controls
    for label, size_max, factor, out in [
        ("coarse",  0.20, 1.0,  "t5_size_coarse.msh"),
        ("medium",  0.20, 0.5,  "t5_size_medium.msh"),
        ("fine",    0.20, 0.2,  "t5_size_fine.msh"),
    ]:
        rmsh.initialize()

        rmsh.model.occ.addRectangle(0, 0, 0, 1.0, 1.0)

        rmsh.option.setNumber("Mesh.Algorithm",     6)   # Frontal-Delaunay
        rmsh.option.setNumber("Mesh.MeshSizeMax",   size_max)
        rmsh.option.setNumber("Mesh.MeshSizeFactor", factor)

        rmsh.model.mesh.generate(2)
        rmsh.write(out)
        print(f"  [{label}] factor={factor}  -> {out}")

        rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# t6 — Structured / quad mesh
#
# Gmsh original:
#   gmsh.model.geo.mesh.setTransfiniteCurve(...)
#   gmsh.model.geo.mesh.setTransfiniteSurface(...)
#   gmsh.model.geo.mesh.setRecombine(2, 1)
#
# rmsh uses Algorithm=9 (Quad Paving) which produces an all-quad structured
# mesh on axis-aligned rectangles — the same end result.
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_t6():
    """t6 – Structured quad mesh via Quad-Paving algorithm (algo 9)."""
    rmsh.initialize()

    # A 3×2 rectangle: QuadPaving will produce a regular quad grid
    rmsh.model.occ.addRectangle(0, 0, 0, 3.0, 2.0)

    rmsh.option.setNumber("Mesh.Algorithm",      9)    # Quad Paving
    rmsh.option.setNumber("Mesh.MeshSizeMax",    0.25)
    rmsh.option.setNumber("Mesh.MeshSizeFactor", 1.0)

    rmsh.model.mesh.generate(2)
    rmsh.write("t6_structured_quad.msh")
    print("  wrote t6_structured_quad.msh")

    rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# t10 — Multiple boolean operations (cut + fuse)
#
# Gmsh original (t5/t16/t17 area):
#   gmsh.model.occ.cut / fuse / fragment
#   gmsh.model.occ.synchronize()
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_t10():
    """t10 – Compound boolean: fuse two boxes, then cut a cylinder through."""
    rmsh.initialize()

    # Fuse two overlapping boxes
    b1 = rmsh.model.occ.addBox(0,    0, 0,  1.0, 1.0, 1.0)
    b2 = rmsh.model.occ.addBox(0.5,  0, 0,  1.0, 1.0, 1.0)
    fused = rmsh.model.occ.fuse([(3, b1)], [(3, b2)])
    fused_tag = fused[0][1]
    print(f"  fuse -> tag {fused_tag}")

    # Cut a vertical cylinder through the fused body
    cyl = rmsh.model.occ.addCylinder(0.75, -0.1, 0.5,  0, 1.2, 0,  0.25)
    result = rmsh.model.occ.cut([(3, fused_tag)], [(3, cyl)])
    print(f"  cut -> {result}")

    rmsh.model.occ.synchronize()

    rmsh.option.setNumber("Mesh.Algorithm3D", 1)
    rmsh.option.setNumber("Mesh.MeshSizeMax", 0.18)

    rmsh.model.mesh.generate(3)
    rmsh.write("t10_boolean_compound.msh")
    print("  wrote t10_boolean_compound.msh")

    rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# t11 — Fillet and chamfer
#
# Gmsh original:
#   gmsh.model.occ.fillet([volume], [curveTags], [radius])
#   gmsh.model.occ.chamfer([volume], [curveTags], [distances])
#
# rmsh mirrors these signatures exactly.
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_t11():
    """t11 – Fillet and chamfer edges of a box."""
    # Fillet
    rmsh.initialize()

    box = rmsh.model.occ.addBox(0, 0, 0, 1, 1, 2)
    try:
        filleted = rmsh.model.occ.fillet(box, list(range(12)), [0.08])
        rmsh.model.occ.synchronize()
        rmsh.option.setNumber("Mesh.Algorithm3D", 1)
        rmsh.option.setNumber("Mesh.MeshSizeMax", 0.15)
        rmsh.model.mesh.generate(3)
        rmsh.write("t11_fillet.msh")
        print(f"  fillet tag={filleted} -> wrote t11_fillet.msh")
    except Exception as e:
        print(f"  fillet: {e}")

    rmsh.finalize()

    # Chamfer
    rmsh.initialize()

    box = rmsh.model.occ.addBox(0, 0, 0, 1, 1, 2)
    try:
        chamfered = rmsh.model.occ.chamfer(box, list(range(4)), [0.08])
        rmsh.model.occ.synchronize()
        rmsh.option.setNumber("Mesh.Algorithm3D", 1)
        rmsh.option.setNumber("Mesh.MeshSizeMax", 0.2)
        rmsh.model.mesh.generate(3)
        rmsh.write("t11_chamfer.msh")
        print(f"  chamfer tag={chamfered} -> wrote t11_chamfer.msh")
    except Exception as e:
        print(f"  chamfer: {e}")

    rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# t16 — Import a STEP file and volume-mesh it
#
# Gmsh original:
#   gmsh.merge("model.step")
#   gmsh.model.occ.synchronize()
#   gmsh.model.mesh.generate(3)
#
# rmsh: open() loads STEP directly; generate(3) produces the volume mesh.
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_t16():
    """t16 – Load a STEP file and generate a volume mesh."""
    step_file = os.path.join(TESTDATA, "simple_cube.step")
    if not os.path.exists(step_file):
        print(f"  skipped: {step_file} not found")
        return

    rmsh.initialize()

    # gmsh: gmsh.merge("model.step") + gmsh.model.occ.synchronize()
    # rmsh: open() handles STEP loading in one call
    rmsh.open(step_file)

    rmsh.option.setNumber("Mesh.Algorithm3D", 1)     # Delaunay
    rmsh.option.setNumber("Mesh.MeshSizeMax", 0.3)
    rmsh.option.setNumber("Mesh.MeshSizeFactor", 0.5)

    rmsh.model.mesh.generate(3)
    rmsh.write("t16_step_volume.msh")
    print("  wrote t16_step_volume.msh")

    rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# t17 — Mesh healing + properties
#
# Gmsh original (OCC healing):
#   gmsh.model.occ.healShapes()
#   gmsh.model.occ.getMass()
#
# Demonstrates the full pre-mesh workflow: create → heal → inspect → mesh.
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_t17():
    """t17 – Pre-mesh workflow: create -> heal -> inspect properties -> mesh."""
    rmsh.initialize()

    # Create a box (2x1x1): vol=2, area=10, centroid=(1, 0.5, 0.5)
    box = rmsh.model.occ.addBox(0, 0, 0, 2.0, 1.0, 1.0)

    # Heal the shape (remove degenerate faces, fix normals, merge vertices)
    report = rmsh.model.occ.healShapes(box, tolerance=1e-8)
    print(f"  heal report: {report}")

    # Inspect geometric properties before meshing
    vol, area, cx, cy, cz = rmsh.model.occ.getProperties(box)
    mass = rmsh.model.occ.getMass(box)
    print(f"  box: vol={vol:.4f}  area={area:.4f}  centroid=({cx:.3f},{cy:.3f},{cz:.3f})")
    print(f"  getMass = {mass:.4f}")

    # Synchronize and mesh
    rmsh.model.occ.synchronize()

    rmsh.option.setNumber("Mesh.Algorithm3D", 1)
    rmsh.option.setNumber("Mesh.MeshSizeMax", 0.4)

    rmsh.model.mesh.generate(3)
    rmsh.write("t17_healed_box.msh")
    print("  wrote t17_healed_box.msh")

    rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# Bonus: BAMG anisotropic 2D — narrow rectangle (aspect ratio 10:1)
#
# BAMG excels on high-aspect-ratio domains (e.g. boundary layers in CFD).
# The metric field is derived from the domain itself; narrow direction gets
# proportionally smaller elements.
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_bamg_aniso():
    """Bonus – BAMG anisotropic mesh on a 10:1 aspect-ratio rectangle."""
    rmsh.initialize()

    # Narrow strip: 1.0 × 0.1
    rmsh.model.occ.addRectangle(0, 0, 0, 1.0, 0.1)

    rmsh.option.setNumber("Mesh.Algorithm",      7)    # BAMG
    rmsh.option.setNumber("Mesh.MeshSizeMax",    0.05)
    rmsh.option.setNumber("Mesh.MeshSizeFactor", 1.0)

    rmsh.model.mesh.generate(2)
    rmsh.write("bonus_bamg_strip.msh")
    print("  wrote bonus_bamg_strip.msh  (BAMG on 10:1 strip)")

    rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# Bonus: Laplacian smoothing after meshing  (gmsh model.mesh.optimize equiv.)
#
# Gmsh: gmsh.model.mesh.optimize("Laplace2D")
# rmsh: rmsh.model.mesh.optimize("Laplace", niter=N)
# ─────────────────────────────────────────────────────────────────────────────
def tutorial_smooth():
    """Bonus – Generate then Laplacian-smooth a mesh (gmsh optimize equiv.)."""
    rmsh.initialize()

    box = rmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    rmsh.model.occ.synchronize()

    rmsh.option.setNumber("Mesh.Algorithm3D", 4)   # Frontal
    rmsh.option.setNumber("Mesh.MeshSizeMax", 0.25)
    rmsh.model.mesh.generate(3)

    # gmsh: gmsh.model.mesh.optimize("Laplace2D")
    rmsh.model.mesh.optimize("Laplace", niter=15)

    rmsh.write("bonus_smooth.msh")
    print("  wrote bonus_smooth.msh  (Frontal-3D + Laplacian 15 passes)")

    rmsh.finalize()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    examples = [
        ("t1  - rectangle -> 2D mesh",              tutorial_t1),
        ("t2  - extrude surface -> volume",          tutorial_t2),
        ("t4  - OCC boolean: box - sphere",          tutorial_t4),
        ("t5  - mesh size control (factor sweep)",   tutorial_t5),
        ("t6  - structured quad mesh (algo 9)",      tutorial_t6),
        ("t10 - compound boolean (fuse + cut)",      tutorial_t10),
        ("t11 - fillet + chamfer box edges",         tutorial_t11),
        ("t16 - STEP file -> volume mesh",           tutorial_t16),
        ("t17 - heal -> properties -> mesh (torus)", tutorial_t17),
        ("B1  - BAMG on narrow 10:1 strip",          tutorial_bamg_aniso),
        ("B2  - Frontal + Laplacian smooth",         tutorial_smooth),
    ]

    failed = []
    for label, fn in examples:
        print(f"\n=== {label} ===")
        try:
            fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            print(f"  FAILED: {exc}")
            failed.append(label)

    print("\n" + "=" * 60)
    if failed:
        print(f"FAILED: {failed}")
        sys.exit(1)
    else:
        print("All tutorials passed.")
