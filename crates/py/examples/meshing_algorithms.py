"""rmsh meshing algorithm examples.

Demonstrates how to select different 2D and 3D meshing algorithms using
`option.setNumber("Mesh.Algorithm", ...)` and `option.setNumber("Mesh.Algorithm3D", ...)`.
Mirrors the Gmsh option API so scripts can be ported with minimal changes.

Run with:
    cd crates/py
    maturin develop --release
    python examples/meshing_algorithms.py

Algorithm codes (Gmsh-compatible):
    2D  -- Mesh.Algorithm
        1  MeshAdapt          (local edge-split/collapse/swap)
        5  Delaunay           (standard Delaunay triangulation)
        6  Frontal-Delaunay   (advancing front + Delaunay, default)
        7  BAMG               (anisotropic, metric-field driven)
        8  Frontal-Quads      (advancing front, quad-dominant)
        9  Packing of Parallelograms / Quad-paving (all-quad on rectangles)

    3D  -- Mesh.Algorithm3D
        1  Delaunay           (Bowyer-Watson + Delaunay refinement, default)
        4  Frontal            (advancing-front tetrahedralization)
        10 HXT                (high-performance parallel Delaunay)
"""

import math
import os
import sys

import rmsh

TESTDATA = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "testdata")
)

STEP_CUBE = os.path.join(TESTDATA, "simple_cube.step")


def _count_elements(path: str) -> dict:
    """Return element/node statistics for a written .msh file by re-opening it."""
    rmsh.initialize()
    rmsh.open(path)
    # Count via mesh stats from a second open (mesh is in current_mesh after open)
    rmsh.finalize()
    return {}   # placeholder – actual element counts available via model.mesh.getElements (stub)


def _write_msh_and_step(stem: str):
    """Write both .msh and .step outputs (used by legacy call sites)."""
    step_out = f"{stem}.step"
    msh_out = f"{stem}.msh"
    rmsh.write(step_out)
    rmsh.write(msh_out)
    return msh_out, step_out


def _write_step_before_meshing(stem: str):
    """Write STEP at geometry stage before meshing (best-effort)."""
    step_out = f"{stem}.step"
    try:
        rmsh.write(step_out)
        print(f"  wrote {step_out}")
    except Exception as e:
        print(f"  skipped {step_out}: {e}")


def _run_2d(algo_id: int, algo_name: str, h: float = 0.15) -> str:
    """
    Build a 2×1 planar rectangle and mesh it with the given 2D algorithm.
    Returns the output file path.
    """
    rmsh.initialize()

    # addRectangle creates a planar surface mesh directly (Gmsh-compatible)
    rmsh.model.occ.addRectangle(0, 0, 0, 2.0, 1.0)

    # Set algorithm and size options, then generate the 2D surface mesh
    rmsh.option.setNumber("Mesh.Algorithm", algo_id)
    rmsh.option.setNumber("Mesh.MeshSizeMax", h)
    rmsh.option.setNumber("Mesh.MeshSizeFactor", 1.0)

    stem = f"algo2d_{algo_id}_{algo_name.replace(' ', '_').lower()}"
    _write_step_before_meshing(stem)

    rmsh.model.mesh.generate(2)

    out = f"{stem}.msh"
    rmsh.write(out)
    rmsh.finalize()
    return out


def _run_3d(algo_id: int, algo_name: str, h: float = 0.4) -> str:
    """
    Build a unit-cube solid and mesh its volume with the given 3D algorithm.
    Returns the output file path.
    """
    rmsh.initialize()

    box = rmsh.model.occ.addBox(0, 0, 0, 1.0, 1.0, 1.0)
    rmsh.model.occ.synchronize()

    rmsh.option.setNumber("Mesh.Algorithm3D", algo_id)
    rmsh.option.setNumber("Mesh.MeshSizeMax", h)
    rmsh.option.setNumber("Mesh.MeshSizeFactor", 1.0)

    stem = f"algo3d_{algo_id}_{algo_name.replace(' ', '_').lower()}"
    _write_step_before_meshing(stem)

    rmsh.model.mesh.generate(3)

    out = f"{stem}.msh"
    rmsh.write(out)
    rmsh.finalize()
    return out


# ---------------------------------------------------------------------------
# 2D algorithm examples
# ---------------------------------------------------------------------------

def example_2d_frontal_delaunay():
    """2D: Frontal-Delaunay (algo 6, Gmsh default)."""
    out = _run_2d(6, "Frontal-Delaunay")
    print(f"  wrote {out}")


def example_2d_delaunay():
    """2D: Pure Delaunay triangulation (algo 5)."""
    out = _run_2d(5, "Delaunay")
    print(f"  wrote {out}")


def example_2d_bamg():
    """2D: BAMG anisotropic meshing (algo 7).

    For a simple domain without an explicit metric field, BAMG falls back
    to an isotropic triangulation with element sizes determined by the
    diagonal metric (equal stretch in X and Y).
    """
    out = _run_2d(7, "BAMG")
    print(f"  wrote {out}")


def example_2d_quad_paving():
    """2D: Quad paving / packing of parallelograms (algo 9).

    On an axis-aligned rectangle, QuadPaving produces a fully structured
    all-quad mesh.  On other domains it falls back to triangles.
    """
    rmsh.initialize()

    # Use a clean rectangle so the structured quad mesher fires
    rmsh.model.occ.addRectangle(0, 0, 0, 3.0, 2.0)

    rmsh.option.setNumber("Mesh.Algorithm", 9)
    rmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)

    _write_step_before_meshing("algo2d_9_quad_paving")

    rmsh.model.mesh.generate(2)
    out = "algo2d_9_quad_paving.msh"
    rmsh.write(out)
    rmsh.finalize()
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# 3D algorithm examples
# ---------------------------------------------------------------------------

def example_3d_delaunay():
    """3D: Delaunay tetrahedralization (algo 1, Gmsh default for 3D)."""
    out = _run_3d(1, "Delaunay")
    print(f"  wrote {out}")


def example_3d_frontal():
    """3D: Frontal-Delaunay 3D (algo 4).

    Delegates to Delaunay internally with slightly tightened quality targets
    (max radius-edge ratio 2.2, min dihedral angle > 0°).
    """
    out = _run_3d(4, "Frontal")
    print(f"  wrote {out}")


def example_3d_hxt():
    """3D: HXT parallel Delaunay (algo 10).

    HXT uses Hilbert-curve-ordered parallel insertion for high throughput.
    With refinement enabled it delegates to Delaunay3D; without refinement
    it uses the centroid-star decomposition.
    """
    out = _run_3d(10, "HXT")
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# Mesh quality: size control + Laplacian smoothing
# ---------------------------------------------------------------------------

def example_size_control():
    """Demonstrate MeshSizeFactor and MeshSizeMin/Max for density control."""
    rmsh.initialize()
    box = rmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    rmsh.model.occ.synchronize()

    # Coarse mesh first
    rmsh.option.setNumber("Mesh.Algorithm3D", 1)
    rmsh.option.setNumber("Mesh.MeshSizeMax", 0.6)
    rmsh.option.setNumber("Mesh.MeshSizeFactor", 1.0)
    _write_step_before_meshing("size_coarse")
    rmsh.model.mesh.generate(3)
    rmsh.write("size_coarse.msh")
    print("  wrote size_coarse.msh")
    rmsh.finalize()

    # Fine mesh
    rmsh.initialize()
    box = rmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    rmsh.model.occ.synchronize()

    rmsh.option.setNumber("Mesh.Algorithm3D", 1)
    rmsh.option.setNumber("Mesh.MeshSizeMax", 0.6)
    rmsh.option.setNumber("Mesh.MeshSizeFactor", 0.3)   # 3× finer
    _write_step_before_meshing("size_fine")
    rmsh.model.mesh.generate(3)
    rmsh.write("size_fine.msh")
    print("  wrote size_fine.msh")
    rmsh.finalize()


def example_laplacian_smoothing():
    """Apply Laplacian smoothing after mesh generation.

    Corresponds to: gmsh.model.mesh.optimize("Laplace2D") or ("Laplace").
    """
    rmsh.initialize()

    box = rmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
    rmsh.model.occ.synchronize()

    rmsh.option.setNumber("Mesh.Algorithm3D", 1)
    rmsh.option.setNumber("Mesh.MeshSizeMax", 0.35)
    _write_step_before_meshing("smooth_before")
    rmsh.model.mesh.generate(3)
    rmsh.write("smooth_before.msh")
    print("  wrote smooth_before.msh")

    # Smooth with 20 Laplacian iterations
    rmsh.model.mesh.optimize("Laplace", niter=20)
    _write_step_before_meshing("smooth_after")
    rmsh.write("smooth_after.msh")
    print("  wrote smooth_after.msh  (Laplacian-smoothed, 20 passes)")

    rmsh.finalize()


# ---------------------------------------------------------------------------
# Comparison: all 2D algorithms on the same domain
# ---------------------------------------------------------------------------

def example_compare_2d_algorithms():
    """Run all available 2D algorithms on the same L-shaped domain and compare."""
    algorithms = [
        (5, "Delaunay"),
        (6, "Frontal-Delaunay"),
        (7, "BAMG"),
        (9, "Quad-Paving"),
    ]

    for algo_id, name in algorithms:
        try:
            out = _run_2d(algo_id, name, h=0.2)
            print(f"  [{algo_id:2d}] {name:22s} -> {out}")
        except Exception as e:
            print(f"  [{algo_id:2d}] {name:22s} -> ERROR: {e}")


# ---------------------------------------------------------------------------
# Comparison: all 3D algorithms on the same box domain
# ---------------------------------------------------------------------------

def example_compare_3d_algorithms():
    """Run all available 3D algorithms on the same box domain and compare."""
    algorithms = [
        (1,  "Delaunay"),
        (4,  "Frontal"),
        (10, "HXT"),
    ]

    for algo_id, name in algorithms:
        try:
            out = _run_3d(algo_id, name, h=0.4)
            print(f"  [{algo_id:2d}] {name:22s} -> {out}")
        except Exception as e:
            print(f"  [{algo_id:2d}] {name:22s} -> ERROR: {e}")


# ---------------------------------------------------------------------------
# STEP file meshing with algorithm selection
# ---------------------------------------------------------------------------

def example_step_with_algorithm():
    """Load a STEP file and mesh with explicitly chosen algorithms."""
    if not os.path.exists(STEP_CUBE):
        print(f"  skipped: {STEP_CUBE} not found")
        return

    rmsh.initialize()
    rmsh.open(STEP_CUBE)

    # Use Frontal 3D with a fine size factor
    rmsh.option.setNumber("Mesh.Algorithm3D", 4)
    rmsh.option.setNumber("Mesh.MeshSizeFactor", 0.5)

    _write_step_before_meshing("step_frontal3d")

    rmsh.model.mesh.generate(3)
    rmsh.write("step_frontal3d.msh")
    print("  wrote step_frontal3d.msh")

    rmsh.finalize()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    examples = [
        # 2D algorithms
        ("2D Frontal-Delaunay (algo 6)",         example_2d_frontal_delaunay),
        ("2D Delaunay (algo 5)",                  example_2d_delaunay),
        ("2D BAMG anisotropic (algo 7)",          example_2d_bamg),
        ("2D Quad paving (algo 9)",               example_2d_quad_paving),
        # 3D algorithms
        ("3D Delaunay (algo 1)",                  example_3d_delaunay),
        ("3D Frontal (algo 4)",                   example_3d_frontal),
        ("3D HXT parallel (algo 10)",             example_3d_hxt),
        # Size control & smoothing
        ("Size control (coarse vs fine)",         example_size_control),
        ("Laplacian smoothing (20 passes)",       example_laplacian_smoothing),
        # Comparison sweeps
        ("Compare all 2D algorithms",             example_compare_2d_algorithms),
        ("Compare all 3D algorithms",             example_compare_3d_algorithms),
        # STEP
        ("STEP + Frontal-3D + size factor 0.5",  example_step_with_algorithm),
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
        print("All meshing algorithm examples passed.")
