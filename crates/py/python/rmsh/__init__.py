"""The rmsh Python API - a gmsh-compatible procedural surface.

Migration from gmsh: replace `import gmsh` with `import rmsh` (the
namespace tree - `rmsh.model.occ.addBox`, `rmsh.model.mesh.generate`,
`rmsh.option.setNumber`, ... - mirrors the gmsh layout). The backend is
the real rcad CAD kernel (primitives, the OCCT-aligned boolean pipeline,
STEP I/O) plus the rmsh mesh model / I-O / meshing crates; element type
numbers in `getElements` are the gmsh MSH type ids.

Deliberate omissions versus the gmsh API (the rmsh charter excludes a
GUI/viewer): `gui.*` and `plugin.*` do not exist here, and gmsh's
built-in `model.geo.*` kernel is not part of v1 - build geometry through
`model.occ.*`. Functions that have no backend yet raise
NotImplementedError with an honest reason instead of silently doing
nothing.

Run with:
    maturin develop --release      (or copy the built cdylib next to
                                    python/rmsh/_rmsh)
    python examples/cad_example.py
"""

from __future__ import annotations

from . import _rmsh

__version__ = getattr(_rmsh, "__version__", "0.0.0")


def _invoke(symbol: str, *args, **kwargs):
    return getattr(_rmsh, symbol)(*args, **kwargs)


# ── root ─────────────────────────────────────────────────────────────


def initialize(argv0="", readConfigFiles=True, runInBackground=False):
    return _invoke("initialize_impl", argv0, readConfigFiles, runInBackground)


def finalize():
    return _invoke("finalize_impl")


def clear():
    return _invoke("clear_impl")


def open(fileName):
    return _invoke("open_impl", fileName)


def merge(fileName):
    return _invoke("merge_impl", fileName)


def write(fileName):
    return _invoke("write_impl", fileName)


# ── model ────────────────────────────────────────────────────────────

class _Namespace:
    pass


class Model:
    def add(self, name):
        return _invoke("model_add", name)

    def remove(self):
        return _invoke("model_remove")

    def getCurrent(self):
        return _invoke("model_get_current")

    def setCurrent(self, name):
        return _invoke("model_set_current", name)

    def getDimension(self):
        return _invoke("model_get_dimension")

    def getEntities(self, dim=-1):
        return _invoke("model_get_entities", dim)

    def getEntitiesInBoundingBox(self, xmin, ymin, zmin, xmax, ymax, zmax, dim=-1):
        return _invoke(
            "model_get_entities_in_bounding_box",
            xmin, ymin, zmin, xmax, ymax, zmax, dim,
        )

    def getBoundingBox(self, dim=-1, tag=-1):
        return _invoke("model_get_bounding_box", dim, tag)

    def getEntityName(self, dim, tag):
        return _invoke("model_get_entity_name", dim, tag)

    def setEntityName(self, dim, tag, name):
        return _invoke("model_set_entity_name", dim, tag, name)

    def addPhysicalGroup(self, dim, tags, tag=-1, name=""):
        return _invoke("model_add_physical_group", dim, tags, tag, name)

    def getPhysicalGroups(self, dim=-1):
        return _invoke("model_get_physical_groups", dim)

    def setPhysicalName(self, dim, tag, name):
        return _invoke("model_set_physical_name", dim, tag, name)

    def getPhysicalName(self, dim, tag):
        return _invoke("model_get_physical_name", dim, tag)


class _Occ:
    def addBox(self, x, y, z, dx, dy, dz, tag=-1):
        return _invoke("model_occ_add_box", x, y, z, dx, dy, dz, tag)

    def addSphere(self, x, y, z, r, tag=-1):
        return _invoke("model_occ_add_sphere", x, y, z, r, tag)

    def addCylinder(self, x, y, z, dx, dy, dz, r, tag=-1):
        return _invoke("model_occ_add_cylinder", x, y, z, dx, dy, dz, r, tag)

    def addCone(self, x, y, z, dx, dy, dz, r1, r2, tag=-1, angle=2 * 3.141592653589793):
        if angle is None or abs(angle - 2 * 3.141592653589793) < 1e-12:
            return _invoke("model_occ_add_cone", x, y, z, dx, dy, dz, r1, r2, tag)
        raise NotImplementedError("occ.addCone with a partial angle is not wired in v1")

    def addTorus(self, *args, **kwargs):
        # legacy rmsh signature: (x, y, z, dx, dy, dz, r1, r2, tag=-1)
        return _invoke("model_occ_add_torus", *args, **kwargs)

    def addRectangle(self, x, y, z, dx, dy, tag=-1):
        return _invoke("model_occ_add_rectangle", x, y, z, dx, dy, tag)

    def addDisk(self, xc, yc, zc, rx, ry, tag=-1):
        return _invoke("model_occ_add_disk", xc, yc, zc, rx, ry, tag)

    def addPoint(self, x, y, z, meshSize=0.0, tag=-1):
        return _invoke("model_occ_add_point", x, y, z, meshSize, tag)

    def addLine(self, startTag, endTag, tag=-1):
        return _invoke("model_occ_add_line", startTag, endTag, tag)

    def addCircleArc(self, startTag, centerTag, endTag, tag=-1):
        return _invoke("model_occ_add_circle_arc", startTag, centerTag, endTag, tag)

    def addSpline(self, pointTags, tag=-1):
        return _invoke("model_occ_add_spline", pointTags, tag)

    def addCurveLoop(self, curveTags, tag=-1):
        return _invoke("model_occ_add_curve_loop", curveTags, tag)

    def addPlaneSurface(self, wireTags, tag=-1):
        return _invoke("model_occ_add_plane_surface", wireTags, tag)

    def addSurfaceLoop(self, surfaceTags, tag=-1):
        return _invoke("model_occ_add_surface_loop", surfaceTags, tag)

    def addVolume(self, shellTags, tag=-1):
        return _invoke("model_occ_add_volume", shellTags, tag)

    def cut(self, objectDimTags, toolDimTags, tag=-1, removeObject=True, removeTool=True):
        return _invoke(
            "model_occ_cut",
            objectDimTags, toolDimTags, tag, removeObject, removeTool,
        )

    def fuse(self, objectDimTags, toolDimTags, tag=-1, removeObject=True, removeTool=True):
        return _invoke(
            "model_occ_fuse",
            objectDimTags, toolDimTags, tag, removeObject, removeTool,
        )

    def intersect(self, objectDimTags, toolDimTags, tag=-1, removeObject=True, removeTool=True):
        return _invoke(
            "model_occ_intersect",
            objectDimTags, toolDimTags, tag, removeObject, removeTool,
        )

    def fragment(self, objectDimTags, toolDimTags, tag=-1, removeObject=True, removeTool=True):
        return _invoke(
            "model_occ_fragment",
            objectDimTags, toolDimTags, tag, removeObject, removeTool,
        )

    def copy(self, dimTags):
        return _invoke("model_occ_copy", dimTags)

    def remove(self, dimTags, recursive=False):
        return _invoke("model_occ_remove", dimTags, recursive)

    def translate(self, dimTags, dx, dy, dz):
        return _invoke("model_occ_translate", dimTags, dx, dy, dz)

    def rotate(self, dimTags, x, y, z, ax, ay, az, angle):
        return _invoke("model_occ_rotate", dimTags, x, y, z, ax, ay, az, angle)

    def dilate(self, dimTags, x, y, z, a, b, c):
        return _invoke("model_occ_dilate", dimTags, x, y, z, a, b, c)

    def mirror(self, dimTags, a, b, c, d):
        return _invoke("model_occ_mirror", dimTags, a, b, c, d)

    def importShapes(self, fileName, highestDimOnly=True, format=""):
        return _invoke("model_occ_import_shapes", fileName, highestDimOnly, format)

    def synchronize(self):
        return _invoke("model_occ_synchronize")

    def getMass(self, tag):
        return _invoke("model_occ_get_mass", tag)

    def getProperties(self, tag):
        return _invoke("model_occ_get_properties", tag)

    def extrude(self, tag, face_idx, dx, dy, dz, distance):
        return _invoke("model_occ_extrude", tag, face_idx, dx, dy, dz, distance)

    def revolve(self, *args, **kwargs):
        return _invoke("model_occ_revolve", *args, **kwargs)

    def fillet(self, *args, **kwargs):
        return _invoke("model_occ_fillet", *args, **kwargs)

    def chamfer(self, *args, **kwargs):
        return _invoke("model_occ_chamfer", *args, **kwargs)

    def healShapes(self, tag, tolerance=1e-8):
        return _invoke("model_occ_heal_shapes", tag, tolerance)


class _Mesh:
    def setSize(self, dimTags, size):
        return _invoke("model_mesh_set_size", dimTags, size)

    def generate(self, dim):
        return _invoke("model_mesh_generate", dim)

    def setOrder(self, order):
        return _invoke("model_mesh_set_order", order)

    def getNodes(self, dim=-1, tag=-1, includeBoundary=False, returnParametricCoord=False):
        return _invoke(
            "model_mesh_get_nodes",
            dim, tag, includeBoundary, returnParametricCoord,
        )

    def getElements(self, dim=-1, tag=-1):
        return _invoke("model_mesh_get_elements", dim, tag)

    def clear(self, dimTags=(), dimTags_n=0):
        return _invoke("model_mesh_clear", list(dimTags), dimTags_n)

    def optimize(self, method="", force=False, niter=1, dimTags=(), dimTags_n=0):
        return _invoke("model_mesh_optimize", method, force, niter, list(dimTags), dimTags_n)

    def refine(self):
        return _invoke("model_mesh_refine")

    def recombine(self, dim=2, tag=-1, angle=45.0):
        return _invoke("model_mesh_recombine", dim, tag, angle)


class _Option:
    def setNumber(self, name, value):
        return _invoke("option_set_number", name, value)

    def getNumber(self, name):
        return _invoke("option_get_number", name)

    def setString(self, name, value):
        return _invoke("option_set_string", name, value)

    def getString(self, name):
        return _invoke("option_get_string", name)

    def setColor(self, name, r, g, b, a=255):
        return _invoke("option_set_color", name, r, g, b, a)

    def getColor(self, name):
        return _invoke("option_get_color", name)

    def restoreDefaults(self):
        return _invoke("option_restore_defaults")


class _Logger:
    def start(self):
        return _invoke("logger_start")

    def stop(self):
        return _invoke("logger_stop")

    def get(self):
        return _invoke("logger_get")


model = Model()
model.occ = _Occ()
model.mesh = _Mesh()
option = _Option()
logger = _Logger()

__all__ = [
    "initialize",
    "finalize",
    "clear",
    "open",
    "merge",
    "write",
    "model",
    "option",
    "logger",
]
