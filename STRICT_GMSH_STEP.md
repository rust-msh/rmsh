# Strict Gmsh STEP Notes

This note records the current strict STEP export invariants used for gmsh/OCCT compatibility.

## Scope

- Applies to the strict mode path enabled by:
  - `STEP.GmshStrict = 1`
  - optional `STEP.Protocol = AP214`
- Focuses on analytic primitives exported from `rmsh.model.occ.add*`.

## Primitive Invariants

### Cylinder

Expected STEP markers in strict mode:

- `ADVANCED_FACE`: 3
- `CYLINDRICAL_SURFACE`: 1
- `SEAM_CURVE`: at least 1

Intent:

- Keep analytic side + two planar caps.
- Keep explicit seam representation for stable viewer behavior.

### Cone / Frustum (r1 != r2)

Expected STEP markers in strict mode:

- `ADVANCED_FACE`: 3
- `CONICAL_SURFACE`: 1
- `EDGE_CURVE`: 3
- no `TRIANGULATED` or `TESSELLATED` fallback entities

Intent:

- Preserve analytic conical side.
- Preserve two cap faces with stable topology.

## Regression Tests

Rust tests are located in:

- `crates/py/src/lib.rs`

Key tests:

- `strict_cylinder_emits_cylindrical_surface_and_seam_curve`
- `strict_frustum_cone_emits_conical_side_and_three_faces`

## 2D Meshing Fallback (Python Backend)

This behavior is related to algorithm examples and `model.mesh.generate(2)`
in the Python backend (`crates/py/src/lib.rs`).

Current rule:

- First try boundary extraction from the temporary surface mesh triangles.
- If mesh-based extraction fails, fallback to extracting a polygon from
  the first available CAD face `outer_wire` (OCC topology).

Intent:

- Keep 2D algorithm examples (`Mesh.Algorithm` 5/6/7/9) runnable even when
  a temporary tessellation does not provide a usable boundary edge set.

Known limitation:

- Fallback currently uses the first available outer wire and does not yet
  represent multiple loops/holes.

## Quick Validation Command

Run the strict primitive exporter and inspect generated files:

```powershell
$env:PYTHONPATH = "c:\Users\lilu\works\rmsh\crates\py\python"
& "c:\Users\lilu\works\rmsh\crates\py\.venv\Scripts\python.exe" "c:\Users\lilu\works\rmsh\crates\py\examples\export_each_3d_step_rmsh_gmsh_strict.py"
```

Generated output folder:

- `crates/py/examples/out_each_3d_step_rmsh_gmsh_strict`
