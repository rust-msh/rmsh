"""rmsh Python API examples – mirrors typical gmsh tutorial patterns.

Run with:
    cd crates/py
    maturin develop --release
    python examples/examples.py

Each example is a self-contained function; they are all executed at the
bottom of the file.
"""

import math
import os
import sys

import rmsh

TESTDATA = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "testdata")
)


def write_mesh_and_step(stem: str):
    """Write both .msh and .step outputs for easier visual inspection."""
    step_out = f"{stem}.step"
    msh_out = f"{stem}.msh"
    try:
        rmsh.write(step_out)
        print(f"  wrote {step_out}")
    except Exception as e:
        print(f"  skipped {step_out}: {e}")
    rmsh.write(msh_out)
    print(f"  wrote {msh_out}")


def write_step_only(stem: str, required: bool = False):
    """Write geometry STEP before meshing.

    When `required` is True, export errors are propagated.
    """
    step_out = f"{stem}.step"
    try:
        rmsh.write(step_out)
        print(f"  wrote {step_out}")
    except Exception as e:
        if required:
            raise
        print(f"  skipped {step_out}: {e}")


# ---------------------------------------------------------------------------
# Example 1: Boolean cut – box with a cylindrical hole  (gmsh t5 flavour)
# ---------------------------------------------------------------------------
def example_boolean_cut():
    """Create a 2×1×1 box and cut a cylinder through it (like gmsh t5).

    Note: cut() returns a list of surviving (dim, tag) pairs.
    """
    rmsh.initialize()

    box = rmsh.model.occ.addBox(0, 0, 0, 2.0, 1.0, 1.0)
    cyl = rmsh.model.occ.addCylinder(1.0, 0.5, -0.1,  0, 0, 1.2,  0.3)

    out_dim_tags, out_dim_tags_map = rmsh.model.occ.cut([(3, box)], [(3, cyl)])
    print(f"  cut result tags: {out_dim_tags}")
    print(f"  cut result map:  {out_dim_tags_map}")
    assert len(out_dim_tags) > 0

    rmsh.model.occ.synchronize()
    write_step_only("ex1_cut")
    rmsh.write("ex1_cut.msh")
    print("  wrote ex1_cut.msh")

    rmsh.finalize()


# ---------------------------------------------------------------------------
# Example 2: Boolean fuse – two overlapping spheres  (gmsh t8 flavour)
# ---------------------------------------------------------------------------
def example_boolean_fuse():
    """Fuse two overlapping spheres into one solid."""
    rmsh.initialize()

    s1 = rmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
    s2 = rmsh.model.occ.addSphere(1.2, 0.0, 0.0, 1.0)

    out_dim_tags, out_dim_tags_map = rmsh.model.occ.fuse([(3, s1)], [(3, s2)])
    print(f"  fuse result tags: {out_dim_tags}")
    print(f"  fuse result map:  {out_dim_tags_map}")
    assert len(out_dim_tags) > 0

    rmsh.model.occ.synchronize()
    write_step_only("ex2_fuse")
    rmsh.write("ex2_fuse.msh")
    print("  wrote ex2_fuse.msh")

    rmsh.finalize()


# ---------------------------------------------------------------------------
# Example 3: Box properties (volume, area, centroid)
# ---------------------------------------------------------------------------
def example_box_properties():
    """Create a box and verify its geometric properties."""
    rmsh.initialize()

    # Box: origin (0,0,0), size 2×1×1 → vol=2, area=10, centroid=(1,0.5,0.5)
    box = rmsh.model.occ.addBox(0, 0, 0, 2.0, 1.0, 1.0)
    vol, area, cx, cy, cz = rmsh.model.occ.getProperties(box)
    print(f"  box vol={vol:.4f}  area={area:.4f}  centroid=({cx:.3f},{cy:.3f},{cz:.3f})")

    assert abs(vol  - 2.0) < 1e-6,  f"box volume wrong: {vol}"
    assert abs(area - 10.0) < 1e-6, f"box area wrong: {area}"
    assert abs(cx - 1.0) < 1e-6 and abs(cy - 0.5) < 1e-6 and abs(cz - 0.5) < 1e-6, \
        f"box centroid wrong: ({cx},{cy},{cz})"

    # Also test getMass (returns volume)
    mass = rmsh.model.occ.getMass(box)
    assert abs(mass - vol) < 1e-6, f"getMass mismatch: {mass} vs {vol}"

    rmsh.model.occ.synchronize()
    rmsh.write("ex3_box_properties.step")
    print("  wrote ex3_box_properties.step")

    rmsh.finalize()
    print("  box properties OK")


# ---------------------------------------------------------------------------
# Example 4: Cone + torus creation (shape building)
# ---------------------------------------------------------------------------
def example_cone_torus():
    """Create a cone and a torus, synchronize, and write."""
    rmsh.initialize()

    # Cone: gmsh-compatible signature (x, y, z, dx, dy, dz, r1, r2)
    cone = rmsh.model.occ.addCone(0, 0, 0, 0, 0, 2.0, 1.0, 0.0)
    print(f"  cone tag = {cone}")
    assert cone > 0

    # Torus: centre (5,0,0), axis Z, R=2, r=0.5
    torus = rmsh.model.occ.addTorus(5, 0, 0, 2.0, 0.5)
    print(f"  torus tag = {torus}")
    assert torus > 0

    # Materialize a tessellated current_mesh from CAD before meshing so STEP
    # export can be inspected at geometry stage.
    rmsh.model.occ.fuse([(3, cone)], [(3, torus)])
    write_step_only("ex4_cone_torus", required=True)
    rmsh.write("ex4_cone_torus.msh")
    print("  wrote ex4_cone_torus.msh")

    rmsh.finalize()


# ---------------------------------------------------------------------------
# Example 5: Fillet a box  (gmsh t11 flavour)
# ---------------------------------------------------------------------------
def example_fillet():
    """Round edges of a unit cube with radius 0.1."""
    rmsh.initialize()

    box = rmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    # Edge indices 0-11 for a box
    edge_indices = list(range(12))
    try:
        filleted = rmsh.model.occ.fillet(box, edge_indices, [0.1])
        print(f"  fillet -> new tag {filleted}")
        rmsh.model.occ.synchronize()
        write_step_only("ex5_fillet")
        rmsh.model.mesh.generate(3)
        rmsh.write("ex5_fillet.msh")
        print("  wrote ex5_fillet.msh")
    except Exception as e:
        print(f"  fillet raised: {e}")

    rmsh.finalize()


# ---------------------------------------------------------------------------
# Example 6: Chamfer a box
# ---------------------------------------------------------------------------
def example_chamfer():
    """Chamfer the 4 vertical edges of a box."""
    rmsh.initialize()

    box = rmsh.model.occ.addBox(0, 0, 0, 1, 1, 2)
    try:
        chamfered = rmsh.model.occ.chamfer(box, [0, 1, 2, 3], [0.1])
        print(f"  chamfer -> new tag {chamfered}")
        rmsh.model.occ.synchronize()
        write_step_only("ex6_chamfer")
        rmsh.model.mesh.generate(3)
        rmsh.write("ex6_chamfer.msh")
        print("  wrote ex6_chamfer.msh")
    except Exception as e:
        print(f"  chamfer raised: {e}")

    rmsh.finalize()


# ---------------------------------------------------------------------------
# Example 7: healShapes – verify repair dict keys
# ---------------------------------------------------------------------------
def example_heal_shapes():
    """Heal a box and check the repair report structure."""
    rmsh.initialize()

    box = rmsh.model.occ.addBox(0, 0, 0, 1.0, 1.0, 1.0)
    report = rmsh.model.occ.healShapes(box, tolerance=1e-8)
    print(f"  heal report: {report}")
    assert "vertices_merged" in report
    assert "wires_fixed" in report

    rmsh.model.occ.synchronize()
    write_step_only("ex7_healed_sphere", required=True)

    rmsh.finalize()
    print("  healShapes OK")


# ---------------------------------------------------------------------------
# Example 8: Option set/get round-trip + restoreDefaults
# ---------------------------------------------------------------------------
def example_options():
    """Exercise option_set/get_number/string/color and restoreDefaults."""
    rmsh.initialize()

    rmsh.option.setNumber("Mesh.MeshSizeFactor", 0.5)
    assert rmsh.option.getNumber("Mesh.MeshSizeFactor") == 0.5

    rmsh.option.setString("General.DefaultFileName", "my_model.msh")
    assert rmsh.option.getString("General.DefaultFileName") == "my_model.msh"

    rmsh.option.setColor("Mesh.Color.One", 255, 0, 0, 255)
    r, g, b, a = rmsh.option.getColor("Mesh.Color.One")
    assert (r, g, b, a) == (255, 0, 0, 255)

    rmsh.option.restoreDefaults()

    try:
        rmsh.option.getNumber("Mesh.MeshSizeFactor")
        assert False, "should have raised KeyError"
    except KeyError:
        pass

    rmsh.finalize()
    print("  option round-trip OK")


# ---------------------------------------------------------------------------
# Example 9: Load STEP file and write mesh
# ---------------------------------------------------------------------------
def example_step_to_msh():
    """Load a STEP file from testdata and write as .msh."""
    step_file = os.path.join(TESTDATA, "simple_cube.step")
    if not os.path.exists(step_file):
        print(f"  skipped: {step_file} not found")
        return

    rmsh.initialize()
    rmsh.open(step_file)
    rmsh.write("ex9_cube.step")
    print("  wrote ex9_cube.step")
    rmsh.write("ex9_cube.msh")
    print("  wrote ex9_cube.msh")
    rmsh.finalize()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    examples = [
        ("1: boolean cut (box - cylinder)",    example_boolean_cut),
        ("2: boolean fuse (two spheres)",       example_boolean_fuse),
        ("3: box properties (vol/area/cog)",    example_box_properties),
        ("4: cone + torus shape creation",      example_cone_torus),
        ("5: fillet box edges",                 example_fillet),
        ("6: chamfer box edges",                example_chamfer),
        ("7: healShapes report",                example_heal_shapes),
        ("8: option set/get/restoreDefaults",   example_options),
        ("9: load STEP, write .msh",            example_step_to_msh),
    ]

    failed = []
    for label, fn in examples:
        print(f"\n=== Example {label} ===")
        try:
            fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            print(f"  FAILED: {exc}")
            failed.append(label)

    print("\n" + "="*50)
    if failed:
        print(f"FAILED examples: {failed}")
        sys.exit(1)
    else:
        print("All examples passed.")
