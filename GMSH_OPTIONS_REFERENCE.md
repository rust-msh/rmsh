# Gmsh Option Keys Reference

This document lists all known Gmsh option keys organized by namespace. rmsh will progressively support these via `option_set_number`, `option_set_string`, `option_set_color`, and their `get_*` counterparts.

For each key: **type**, **default**, **description**, and **rmsh priority**.

---

## General.*

| Key | Type | Default | Description | rmsh Priority |
|-----|------|---------|-------------|---------------|
| `General.Verbosity` | number | 5 | Verbosity level (0=silent, 5=normal, 10=debug) | High |
| `General.Terminal` | number | 0 | Log to terminal (1=yes) | Low |
| `General.NoPopup` | number | 0 | Disable popup dialogs | Low |
| `General.AbortOnError` | number | 0 | Abort on error (0=no, 1=warn, 2=abort) | Medium |
| `General.ExpertMode` | number | 0 | Enable expert mode (suppress certain warnings) | Low |
| `General.NumThreads` | number | 0 | Number of OpenMP threads (0=system default) | High |
| `General.GraphicsWidth` | number | 800 | Window width in pixels | Low |
| `General.GraphicsHeight` | number | 600 | Window height in pixels | Low |
| `General.GraphicsFontSize` | number | 12 | Font size for graphics | Low |
| `General.GraphicsFontSizeTitle` | number | 14 | Font size for titles | Low |
| `General.TrackballHyperbolicSheet` | number | 1 | Use hyperbolic sheet for trackball rotation | Low |
| `General.RotationCenterX` | number | 0 | X coordinate of rotation center | Low |
| `General.RotationCenterY` | number | 0 | Y coordinate of rotation center | Low |
| `General.RotationCenterZ` | number | 0 | Z coordinate of rotation center | Low |
| `General.ScaleX` | number | 1 | X scale factor for display | Low |
| `General.ScaleY` | number | 1 | Y scale factor for display | Low |
| `General.ScaleZ` | number | 1 | Z scale factor for display | Low |
| `General.TranslationX` | number | 0 | X translation for display | Low |
| `General.TranslationY` | number | 0 | Y translation for display | Low |
| `General.FileName` | string | "" | Current model file name (read-only) | Low |
| `General.DefaultFileName` | string | "untitled.geo" | Default file name | Low |

---

## Geometry.*

| Key | Type | Default | Description | rmsh Priority |
|-----|------|---------|-------------|---------------|
| `Geometry.Tolerance` | number | 1e-8 | Geometric tolerance for model operations | **Critical** |
| `Geometry.ToleranceBoolean` | number | 0 | Tolerance for boolean operations (0=auto) | High |
| `Geometry.OCCTargetMinimumDistance` | number | 0 | Minimum distance target for OCC healing | High |
| `Geometry.OCCBooleanPreserveNumbering` | number | 1 | Try to preserve entity tags after booleans | High |
| `Geometry.OCCImportLabels` | number | 1 | Import labels/names from OCC/STEP files | High |
| `Geometry.OCCSewFaces` | number | 1 | Sew OCC faces when creating geometry | Medium |
| `Geometry.OCCFixDegenerated` | number | 0 | Fix degenerated edges/faces in OCC shapes | High |
| `Geometry.OCCFixSmallEdges` | number | 0 | Fix small edges in OCC shapes | High |
| `Geometry.OCCFixSmallFaces` | number | 0 | Fix small faces in OCC shapes | High |
| `Geometry.OCCScaling` | number | 1 | Scale factor when reading OCC/STEP files | High |
| `Geometry.AutoCoherence` | number | 1 | Auto merge geometrically identical entities | Medium |
| `Geometry.MatchGeomAndMesh` | number | 0 | Match geometry and mesh | Low |
| `Geometry.Points` | number | 1 | Display geometry points | Low |
| `Geometry.Lines` | number | 1 | Display geometry lines/curves | Low |
| `Geometry.Surfaces` | number | 0 | Display geometry surfaces | Low |
| `Geometry.Volumes` | number | 0 | Display geometry volumes | Low |
| `Geometry.PointNumbers` | number | 0 | Display point entity tags | Low |
| `Geometry.LineNumbers` | number | 0 | Display curve entity tags | Low |
| `Geometry.SurfaceNumbers` | number | 0 | Display surface entity tags | Low |
| `Geometry.VolumeNumbers` | number | 0 | Display volume entity tags | Low |

---

## Mesh.*

| Key | Type | Default | Description | rmsh Priority |
|-----|------|---------|-------------|---------------|
| `Mesh.Algorithm` | number | 6 | 2D meshing algorithm: 1=MeshAdapt, 2=Auto, 5=Delaunay, 6=Frontal-Delaunay, 7=BAMG, 8=Frontal-Delaunay-Quads, 9=Packing-Parallelograms | **Critical** |
| `Mesh.Algorithm3D` | number | 1 | 3D meshing algorithm: 1=Delaunay, 3=Initial-mesh-only, 4=Frontal, 7=MMG3D, 9=R-tree, 10=HXT | **Critical** |
| `Mesh.MeshSizeFactor` | number | 1 | Global mesh size scale factor | **Critical** |
| `Mesh.MeshSizeMin` | number | 0 | Minimum mesh element size (0=no limit) | **Critical** |
| `Mesh.MeshSizeMax` | number | 1e22 | Maximum mesh element size | **Critical** |
| `Mesh.MeshSizeFromPoints` | number | 1 | Compute mesh size from point sizes | High |
| `Mesh.MeshSizeFromCurvature` | number | 0 | Compute size based on local curvature | High |
| `Mesh.MeshSizeExtendFromBoundary` | number | 1 | Extend mesh sizes from boundary | High |
| `Mesh.CharacteristicLengthFactor` | number | 1 | *Deprecated*, use MeshSizeFactor | Low |
| `Mesh.CharacteristicLengthMin` | number | 0 | *Deprecated*, use MeshSizeMin | Low |
| `Mesh.CharacteristicLengthMax` | number | 1e22 | *Deprecated*, use MeshSizeMax | Low |
| `Mesh.Smoothing` | number | 1 | Number of Laplacian smoothing steps | High |
| `Mesh.SmoothNormals` | number | 0 | Smooth mesh normals for display | Low |
| `Mesh.ElementOrder` | number | 1 | Polynomial order of mesh elements (1=linear, 2=quadratic, ...) | **Critical** |
| `Mesh.HighOrderOptimize` | number | 0 | Optimize high-order meshes (0=none, 1=optimization, 2=elastic+optimization) | High |
| `Mesh.SecondOrderIncomplete` | number | 0 | Use incomplete second-order elements (8-node hex, etc.) | Medium |
| `Mesh.SecondOrderLinear` | number | 0 | Move second-order points using linear interpolation | Medium |
| `Mesh.Optimize` | number | 1 | Optimize mesh to improve quality | High |
| `Mesh.OptimizeNetgen` | number | 0 | Optimize mesh using Netgen library | Medium |
| `Mesh.OptimizeThreshold` | number | 0.3 | Optimize elements with quality below this threshold | High |
| `Mesh.RecombineAll` | number | 0 | Recombine triangles into quads (all surfaces) | High |
| `Mesh.RecombinationAlgorithm` | number | 1 | Recombination algorithm: 0=simple, 1=blossom, 2=simple-full, 3=blossom-full | Medium |
| `Mesh.Recombine3DAll` | number | 0 | Recombine 3D mesh into hexahedra | Medium |
| `Mesh.SubdivisionAlgorithm` | number | 0 | Mesh subdivision: 0=none, 1=all-quads, 2=all-hexas, 3=barycentric | Medium |
| `Mesh.Format` | number | 10 | Output format: 1=.msh1, 2=.msh2, 4=.med, 10=auto, 16=.vtk, 19=.vrml, 21=.mail, 26=.bdf, 27=.cgns, 28=.med, 30=.mesh, 31=.bdf, 32=.diff, 33=.inp, 39=.su2 | **Critical** |
| `Mesh.MshFileVersion` | number | 4.1 | Version of the MSH file format (2.2 or 4.1) | **Critical** |
| `Mesh.MshFilePartitioned` | number | 0 | Write mesh partitions in separate files | Low |
| `Mesh.Binary` | number | 0 | Write mesh in binary format | High |
| `Mesh.SaveAll` | number | 0 | Save all elements (including duplicates) | Medium |
| `Mesh.SaveParametric` | number | 0 | Save parametric coordinates in mesh | Low |
| `Mesh.SaveGroupsOfNodes` | number | 0 | Save groups of nodes as physical entities | Medium |
| `Mesh.Partitioner` | number | 1 | Mesh partitioner: 1=Metis, 2=SimplePartition | Medium |
| `Mesh.NbPartitions` | number | 1 | Number of mesh partitions | Medium |
| `Mesh.ColorCarousel` | number | 1 | Color mesh by: 0=solid, 1=by element type, 2=by elementary entity, 3=by partition | Low |
| `Mesh.AngleSmoothNormals` | number | 30 | Threshold angle (deg) for computing smooth normals | Low |
| `Mesh.AngleToleranceFacetOverlap` | number | 0.1 | Tolerance angle for mesh facet overlap | Low |
| `Mesh.NbNodes` | number | — | Number of mesh nodes (read-only) | Low |
| `Mesh.NbTriangles` | number | — | Number of triangle elements (read-only) | Low |
| `Mesh.NbTetrahedra` | number | — | Number of tetrahedral elements (read-only) | Low |
| `Mesh.NbQuadrangles` | number | — | Number of quad elements (read-only) | Low |
| `Mesh.NbHexahedra` | number | — | Number of hex elements (read-only) | Low |
| `Mesh.NbPrisms` | number | — | Number of prism elements (read-only) | Low |
| `Mesh.NbPyramids` | number | — | Number of pyramid elements (read-only) | Low |
| `Mesh.Points` | number | 0 | Display mesh nodes | Low |
| `Mesh.Lines` | number | 0 | Display mesh line elements | Low |
| `Mesh.SurfaceEdges` | number | 1 | Display surface mesh edges | Low |
| `Mesh.SurfaceFaces` | number | 0 | Display surface mesh faces (filled) | Low |
| `Mesh.VolumeEdges` | number | 1 | Display volume mesh edges | Low |
| `Mesh.VolumeFaces` | number | 0 | Display volume mesh faces (filled) | Low |
| `Mesh.PointNumbers` | number | 0 | Display node tags | Low |
| `Mesh.LineNumbers` | number | 0 | Display line element tags | Low |
| `Mesh.SurfaceNumbers` | number | 0 | Display surface element tags | Low |
| `Mesh.VolumeNumbers` | number | 0 | Display volume element tags | Low |

---

## Solver.*

| Key | Type | Default | Description | rmsh Priority |
|-----|------|---------|-------------|---------------|
| `Solver.AutoMesh` | number | 1 | Automatically (re-)mesh before solving | Low |
| `Solver.AutoMeshThreshold` | number | 0.3 | Quality threshold triggering auto-remesh | Low |
| `Solver.AlwaysSave` | number | 0 | Always save mesh before launching solver | Low |
| `Solver.AutoCheck` | number | 1 | Auto-check solver input | Low |
| `Solver.Executable0`–`Solver.Executable9` | string | "" | Path to solver executable | Low |
| `Solver.Name0`–`Solver.Name9` | string | "" | Solver name | Low |
| `Solver.RemoteLogin` | string | "" | Login for remote solver | Low |

---

## PostProcessing.*

| Key | Type | Default | Description | rmsh Priority |
|-----|------|---------|-------------|---------------|
| `PostProcessing.Format` | number | 1 | Default post-processing format | Low |
| `PostProcessing.Binary` | number | 0 | Write post-processing in binary format | Low |
| `PostProcessing.AnimationDelay` | number | 0.1 | Delay between animation frames (s) | Low |
| `PostProcessing.AnimationCycle` | number | 0 | Animation cycle mode | Low |
| `PostProcessing.CombineRemoveOriginal` | number | 1 | Remove original views after combining | Low |
| `PostProcessing.Smoothing` | number | 0 | Apply smoothing to post-processing data | Low |
| `PostProcessing.HorizontalScales` | number | 1 | Show horizontal color scales | Low |

---

## View.*

These apply to individual views (post-processing data sets). Access as `View[i].Key` or `View.Key` for defaults.

| Key | Type | Default | Description | rmsh Priority |
|-----|------|---------|-------------|---------------|
| `View.Type` | number | 1 | View type: 1=3D, 2=2D space, 3=2D time | Low |
| `View.Visible` | number | 1 | Show/hide view | Low |
| `View.Format` | string | "%.3g" | Number display format | Low |
| `View.Name` | string | "" | View name | Low |
| `View.Axes` | number | 0 | Axes display: 0=none, 1=simple axes, 2=box, 3=full grid, 4=open grid, 5=ruler | Low |
| `View.IntervalsType` | number | 2 | Contour interval type: 1=iso, 2=continuous, 3=discrete, 4=numeric | Low |
| `View.NbIso` | number | 10 | Number of iso-contours/intervals | Low |
| `View.RangeType` | number | 0 | Range type: 0=default, 1=custom, 2=per step | Low |
| `View.CustomMin` | number | 0 | Custom min for color scale | Low |
| `View.CustomMax` | number | 1 | Custom max for color scale | Low |
| `View.ColormapNumber` | number | 2 | Colormap: 0=grey, 1=hot, 2=cool, 3=rainbow, ... | Low |
| `View.ColormapInvert` | number | 0 | Invert colormap | Low |
| `View.ShowScale` | number | 1 | Show color scale legend | Low |
| `View.ScaleType` | number | 1 | Scale type: 1=linear, 2=logarithmic, 3=double log | Low |
| `View.VectorType` | number | 1 | Vector display: 1=segment, 2=arrow, 3=pyramid, 4=3D arrow, 5=displacement, 6=comet | Low |
| `View.ArrowSizeMax` | number | 60 | Maximum arrow size | Low |
| `View.ArrowSizeMin` | number | 0 | Minimum arrow size | Low |
| `View.Explode` | number | 1 | Explode factor for display | Low |
| `View.Light` | number | 1 | Enable lighting for view | Low |
| `View.SmoothNormals` | number | 0 | Smooth normals in view | Low |
| `View.OffsetX/Y/Z` | number | 0 | Translation offset for view | Low |
| `View.RaiseX/Y/Z` | number | 0 | Extrusion raise for 2D view | Low |
| `View.TimeStep` | number | 0 | Current time step | Low |
| `View.DrawStrings` | number | 1 | Draw string annotations | Low |
| `View.DrawPoints` | number | 1 | Draw point data | Low |
| `View.DrawLines` | number | 1 | Draw line data | Low |
| `View.DrawTriangles` | number | 1 | Draw triangle data | Low |
| `View.DrawTetrahedra` | number | 1 | Draw tetrahedral data | Low |
| `View.DrawScalars` | number | 1 | Draw scalar field data | Low |
| `View.DrawVectors` | number | 1 | Draw vector field data | Low |
| `View.DrawTensors` | number | 1 | Draw tensor field data | Low |

---

## rmsh Implementation Priority

### Phase 1 — Critical (mesh workflow)
- `Mesh.Algorithm`, `Mesh.Algorithm3D`
- `Mesh.MeshSizeFactor`, `Mesh.MeshSizeMin`, `Mesh.MeshSizeMax`
- `Mesh.ElementOrder`
- `Mesh.Format`, `Mesh.MshFileVersion`, `Mesh.Binary`
- `Geometry.Tolerance`
- `General.Verbosity`, `General.NumThreads`

### Phase 2 — High (quality and CAD)
- `Mesh.Optimize`, `Mesh.OptimizeThreshold`, `Mesh.OptimizeNetgen`
- `Mesh.RecombineAll`, `Mesh.RecombinationAlgorithm`
- `Mesh.MeshSizeFromCurvature`, `Mesh.MeshSizeExtendFromBoundary`
- `Mesh.Smoothing`
- `Geometry.OCCTargetMinimumDistance`, `Geometry.OCCBooleanPreserveNumbering`
- `Geometry.OCCImportLabels`, `Geometry.OCCFixDegenerated/SmallEdges/SmallFaces`
- `Geometry.OCCScaling`

### Phase 3 — Medium (advanced meshing)
- `Mesh.HighOrderOptimize`
- `Mesh.SecondOrderIncomplete`, `Mesh.SecondOrderLinear`
- `Mesh.SubdivisionAlgorithm`
- `Mesh.Partitioner`, `Mesh.NbPartitions`
- `Geometry.AutoCoherence`

### Phase 4 — Low (display / UI options)
- All `General.Graphics*`, `Mesh.*Numbers`, `Geometry.*Numbers` display options
- All `View.*` options (post-processing display)
- All `Solver.*` options

---

## Notes

- Option types: **number** (double), **string**, **color** (RGB int)
- `option_restore_defaults` resets all number/string/color options to defaults
- `Mesh.NbNodes`, `Mesh.NbTriangles`, etc. are read-only stats populated after meshing
- Deprecated `CharacteristicLength*` keys map to corresponding `MeshSize*` keys
- `View[i].Key` syntax accesses per-view options by index (future: rmsh may use named views)

---

## Current Alignment Snapshot (2026-04-19)

This section tracks the observed alignment between gmsh options and rmsh options
for keys that rmsh currently consumes or exposes in common workflows.

Legend:

- `Aligned`: behavior matches gmsh sufficiently for normal scripts.
- `Compatible`: behavior differs, but still predictable and acceptable.
- `Divergent`: behavior mismatch that can break gmsh-oriented scripts.

### Option API Surface

| Area | gmsh | rmsh | Status |
|------|------|------|--------|
| `setNumber/getNumber` | yes | yes | Aligned |
| `setString/getString` | yes | yes | Aligned |
| `setColor/getColor` | yes | yes | Aligned |
| `restoreDefaults` | yes | yes | Aligned |

### Default Read Semantics (Known Keys)

| Key | gmsh default read | rmsh default read | Status |
|-----|-------------------|-------------------|--------|
| `Mesh.MeshSizeMax` | `1e22` | `1e22` | Aligned |
| `Mesh.MeshSizeMin` | `0` | `0` | Aligned |
| `Mesh.Algorithm` | `6` | `6` | Aligned |
| `Mesh.Algorithm3D` | `1` | `1` | Aligned |
| `Geometry.Points` | `1` | `1` | Aligned |
| `General.Background` (color) | `(255,255,255,255)` after restore | `(255,255,255,255)` after restore | Aligned |

### rmsh STEP Extension Keys

These keys are rmsh extensions and not native gmsh options.

| Key | gmsh | rmsh | Status |
|-----|------|------|--------|
| `STEP.GmshStrict` | unknown option | default `1` + overridable | Compatible |
| `STEP.Protocol` | unknown option | default `AP214` + overridable | Compatible |
| `STEP.WriterStyle` | unknown option | default `gmsh-strict` + overridable | Compatible |

### Lifecycle Semantics

| Scenario | gmsh | rmsh | Status |
|----------|------|------|--------|
| Option API before `initialize` | requires initialized context | requires initialized context | Aligned |
| `restoreDefaults` then `get*` known key | returns built-in default | returns built-in default | Aligned |

### Option Consumption in rmsh Runtime

The following keys are actively consumed in code paths (not just stored):

- STEP export profile: `STEP.Protocol`, `Geometry.OCCStepProtocol`, `Geometry.OCCSTEPProtocol`
- STEP strict toggle: `STEP.GmshStrict`, `Geometry.OCCStepGmshStrict`
- STEP color: `Geometry.OCCSolidColor`, `STEP.SolidColor`
- Mesh size and aliases: `Mesh.MeshSizeMin/Max/Factor`, `Mesh.CharacteristicLengthMin/Max/Factor`
- Mesh algorithms: `Mesh.Algorithm`, `Mesh.Algorithm3D`
- Viewer toggles: `Geometry.Points`, `Mesh.Points`, `Geometry.Curves`, `Geometry.Lines`, `Geometry.Surfaces`, `Geometry.Volumes`

When `Mesh.MeshSizeMax` is not explicitly set, rmsh now derives a conservative
default target size from the current geometry scale (surface mesh AABB diagonal)
instead of using a fixed constant.

Keys not consumed by runtime logic are currently stored and retrievable, but
have no behavioral effect.

### Remaining Gaps

- The global gmsh option universe is much larger than rmsh's actively consumed
	subset; many keys are currently no-op in rmsh.
- Read-only statistics keys (`Mesh.Nb*`) are not yet modeled as live option
	values in the option API.
