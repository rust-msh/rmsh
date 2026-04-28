"""Auto-generated rmsh Python API.

This file is generated from api/rmshc.h by crates/rmsh-py/build.rs.
Do not edit manually.
"""

from __future__ import annotations

from rmsh import _rmsh

__version__ = getattr(_rmsh, "__version__", "0.0.0")

def _invoke(symbol: str, *args, **kwargs):
    return getattr(_rmsh, symbol)(*args, **kwargs)

class _Namespace:
    def __init__(self, name: str):
        self._name = name

def clear(*args, **kwargs):
    return _invoke("clear", *args, **kwargs)

def finalize(*args, **kwargs):
    return _invoke("finalize", *args, **kwargs)

def initialize(*args, **kwargs):
    return _invoke("initialize", *args, **kwargs)

def merge(*args, **kwargs):
    return _invoke("merge", *args, **kwargs)

def open(*args, **kwargs):
    return _invoke("open", *args, **kwargs)

def write(*args, **kwargs):
    return _invoke("write", *args, **kwargs)

class _NsGui(_Namespace):
    pass

class _NsLogger(_Namespace):
    pass

class _NsModel(_Namespace):
    pass

class _NsModelGeo(_Namespace):
    pass

class _NsModelMesh(_Namespace):
    pass

class _NsModelOcc(_Namespace):
    pass

class _NsOption(_Namespace):
    pass

class _NsPlugin(_Namespace):
    pass

gui = _NsGui("gui")
logger = _NsLogger("logger")
model = _NsModel("model")
model_geo = _NsModelGeo("model.geo")
model_mesh = _NsModelMesh("model.mesh")
model_occ = _NsModelOcc("model.occ")
option = _NsOption("option")
plugin = _NsPlugin("plugin")

model.geo = model_geo
model.mesh = model_mesh
model.occ = model_occ

def _gui_initialize(*args, **kwargs):
    """rmshGuiInitialize(ierr)"""
    return _invoke("gui_initialize", *args, **kwargs)

gui.initialize = _gui_initialize

def _gui_run(*args, **kwargs):
    """rmshGuiRun(ierr)"""
    return _invoke("gui_run", *args, **kwargs)

gui.run = _gui_run

def _gui_wait(*args, **kwargs):
    """rmshGuiWait(time, ierr)"""
    return _invoke("gui_wait", *args, **kwargs)

gui.wait = _gui_wait

def _logger_get(*args, **kwargs):
    """rmshLoggerGet(log, log_n, ierr)"""
    return _invoke("logger_get", *args, **kwargs)

logger.get = _logger_get

def _logger_start(*args, **kwargs):
    """rmshLoggerStart(ierr)"""
    return _invoke("logger_start", *args, **kwargs)

logger.start = _logger_start

def _logger_stop(*args, **kwargs):
    """rmshLoggerStop(ierr)"""
    return _invoke("logger_stop", *args, **kwargs)

logger.stop = _logger_stop

def _model_add(*args, **kwargs):
    """rmshModelAdd(name, ierr)"""
    return _invoke("model_add", *args, **kwargs)

model.add = _model_add

def _model_addPhysicalGroup(*args, **kwargs):
    """rmshModelAddPhysicalGroup(dim, tags, tags_n, tag, name, ierr)"""
    return _invoke("model_add_physical_group", *args, **kwargs)

model.addPhysicalGroup = _model_addPhysicalGroup

def _model_getBoundingBox(*args, **kwargs):
    """rmshModelGetBoundingBox(dim, tag, xmin, ymin, zmin, xmax, ymax, zmax, ierr)"""
    return _invoke("model_get_bounding_box", *args, **kwargs)

model.getBoundingBox = _model_getBoundingBox

def _model_getCurrent(*args, **kwargs):
    """rmshModelGetCurrent(name, ierr)"""
    return _invoke("model_get_current", *args, **kwargs)

model.getCurrent = _model_getCurrent

def _model_getDimension(*args, **kwargs):
    """rmshModelGetDimension(dim, ierr)"""
    return _invoke("model_get_dimension", *args, **kwargs)

model.getDimension = _model_getDimension

def _model_getEntities(*args, **kwargs):
    """rmshModelGetEntities(dimTags, dimTags_n, dim, ierr)"""
    return _invoke("model_get_entities", *args, **kwargs)

model.getEntities = _model_getEntities

def _model_getEntitiesInBoundingBox(*args, **kwargs):
    """getEntitiesInBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax, dim=-1) -> dimTags"""
    return _invoke("model_get_entities_in_bounding_box", *args, **kwargs)

model.getEntitiesInBoundingBox = _model_getEntitiesInBoundingBox

def _model_getEntityName(*args, **kwargs):
    """rmshModelGetEntityName(dim, tag, name, ierr)"""
    return _invoke("model_get_entity_name", *args, **kwargs)

model.getEntityName = _model_getEntityName

def _model_getPhysicalGroups(*args, **kwargs):
    """rmshModelGetPhysicalGroups(dimTags, dimTags_n, dim, ierr)"""
    return _invoke("model_get_physical_groups", *args, **kwargs)

model.getPhysicalGroups = _model_getPhysicalGroups

def _model_getPhysicalName(*args, **kwargs):
    """rmshModelGetPhysicalName(dim, tag, name, ierr)"""
    return _invoke("model_get_physical_name", *args, **kwargs)

model.getPhysicalName = _model_getPhysicalName

def _model_remove(*args, **kwargs):
    """rmshModelRemove(ierr)"""
    return _invoke("model_remove", *args, **kwargs)

model.remove = _model_remove

def _model_setCurrent(*args, **kwargs):
    """rmshModelSetCurrent(name, ierr)"""
    return _invoke("model_set_current", *args, **kwargs)

model.setCurrent = _model_setCurrent

def _model_setEntityName(*args, **kwargs):
    """rmshModelSetEntityName(dim, tag, name, ierr)"""
    return _invoke("model_set_entity_name", *args, **kwargs)

model.setEntityName = _model_setEntityName

def _model_setPhysicalName(*args, **kwargs):
    """rmshModelSetPhysicalName(dim, tag, name, ierr)"""
    return _invoke("model_set_physical_name", *args, **kwargs)

model.setPhysicalName = _model_setPhysicalName

def _model_geo_addCurveLoop(*args, **kwargs):
    """rmshModelGeoAddCurveLoop(curveTags, curveTags_n, tag, ierr)"""
    return _invoke("model_geo_add_curve_loop", *args, **kwargs)

model_geo.addCurveLoop = _model_geo_addCurveLoop

def _model_geo_addLine(*args, **kwargs):
    """rmshModelGeoAddLine(startTag, endTag, tag, ierr)"""
    return _invoke("model_geo_add_line", *args, **kwargs)

model_geo.addLine = _model_geo_addLine

def _model_geo_addPlaneSurface(*args, **kwargs):
    """rmshModelGeoAddPlaneSurface(wireTags, wireTags_n, tag, ierr)"""
    return _invoke("model_geo_add_plane_surface", *args, **kwargs)

model_geo.addPlaneSurface = _model_geo_addPlaneSurface

def _model_geo_addPoint(*args, **kwargs):
    """rmshModelGeoAddPoint(x, y, z, meshSize, tag, ierr)"""
    return _invoke("model_geo_add_point", *args, **kwargs)

model_geo.addPoint = _model_geo_addPoint

def _model_geo_synchronize(*args, **kwargs):
    """rmshModelGeoSynchronize(ierr)"""
    return _invoke("model_geo_synchronize", *args, **kwargs)

model_geo.synchronize = _model_geo_synchronize

def _model_mesh_clear(*args, **kwargs):
    """rmshModelMeshClear(dimTags, dimTags_n, ierr)"""
    return _invoke("model_mesh_clear", *args, **kwargs)

model_mesh.clear = _model_mesh_clear

def _model_mesh_generate(*args, **kwargs):
    """rmshModelMeshGenerate(dim, ierr)"""
    return _invoke("model_mesh_generate", *args, **kwargs)

model_mesh.generate = _model_mesh_generate

def _model_mesh_getElements(*args, **kwargs):
    """rmshModelMeshGetElements(elementTypes_n, elementTags_n, nodeTags_n, dim, tag, ierr)"""
    return _invoke("model_mesh_get_elements", *args, **kwargs)

model_mesh.getElements = _model_mesh_getElements

def _model_mesh_getNodes(*args, **kwargs):
    """rmshModelMeshGetNodes(nodeTags_n, coord_n, parametricCoord_n, dim, tag, includeBoundary, returnParametricCoord, ierr)"""
    return _invoke("model_mesh_get_nodes", *args, **kwargs)

model_mesh.getNodes = _model_mesh_getNodes

def _model_mesh_optimize(*args, **kwargs):
    """rmshModelMeshOptimize(method, force, niter, dimTags, dimTags_n, ierr)"""
    return _invoke("model_mesh_optimize", *args, **kwargs)

model_mesh.optimize = _model_mesh_optimize

def _model_mesh_recombine(*args, **kwargs):
    """rmshModelMeshRecombine(dim, tag, angle, ierr)"""
    return _invoke("model_mesh_recombine", *args, **kwargs)

model_mesh.recombine = _model_mesh_recombine

def _model_mesh_refine(*args, **kwargs):
    """rmshModelMeshRefine(ierr)"""
    return _invoke("model_mesh_refine", *args, **kwargs)

model_mesh.refine = _model_mesh_refine

def _model_mesh_setOrder(*args, **kwargs):
    """rmshModelMeshSetOrder(order, ierr)"""
    return _invoke("model_mesh_set_order", *args, **kwargs)

model_mesh.setOrder = _model_mesh_setOrder

def _model_mesh_setSize(*args, **kwargs):
    """rmshModelMeshSetSize(dimTags, dimTags_n, size, ierr)"""
    return _invoke("model_mesh_set_size", *args, **kwargs)

model_mesh.setSize = _model_mesh_setSize

def _model_occ_addBox(*args, **kwargs):
    """rmshModelOccAddBox(x, y, z, dx, dy, dz, tag, ierr)"""
    return _invoke("model_occ_add_box", *args, **kwargs)

model_occ.addBox = _model_occ_addBox

def _model_occ_addPoint(*args, **kwargs):
    """addPoint(x, y, z, meshSize=0, tag=-1) -> tag"""
    return _invoke("model_occ_add_point", *args, **kwargs)

model_occ.addPoint = _model_occ_addPoint

def _model_occ_addLine(*args, **kwargs):
    """addLine(startTag, endTag, tag=-1) -> tag"""
    return _invoke("model_occ_add_line", *args, **kwargs)

model_occ.addLine = _model_occ_addLine

def _model_occ_addCircleArc(*args, **kwargs):
    """addCircleArc(startTag, centerTag, endTag, tag=-1) -> tag"""
    return _invoke("model_occ_add_circle_arc", *args, **kwargs)

model_occ.addCircleArc = _model_occ_addCircleArc

def _model_occ_addSpline(*args, **kwargs):
    """addSpline(pointTags, tag=-1) -> tag"""
    return _invoke("model_occ_add_spline", *args, **kwargs)

model_occ.addSpline = _model_occ_addSpline

def _model_occ_addCurveLoop(*args, **kwargs):
    """addCurveLoop(curveTags, tag=-1) -> tag"""
    return _invoke("model_occ_add_curve_loop", *args, **kwargs)

model_occ.addCurveLoop = _model_occ_addCurveLoop

def _model_occ_addPlaneSurface(*args, **kwargs):
    """addPlaneSurface(wireTags, tag=-1) -> tag"""
    return _invoke("model_occ_add_plane_surface", *args, **kwargs)

model_occ.addPlaneSurface = _model_occ_addPlaneSurface

def _model_occ_addSurfaceLoop(*args, **kwargs):
    """addSurfaceLoop(surfaceTags, tag=-1) -> tag"""
    return _invoke("model_occ_add_surface_loop", *args, **kwargs)

model_occ.addSurfaceLoop = _model_occ_addSurfaceLoop

def _model_occ_addVolume(*args, **kwargs):
    """addVolume(shellTags, tag=-1) -> tag"""
    return _invoke("model_occ_add_volume", *args, **kwargs)

model_occ.addVolume = _model_occ_addVolume

def _model_occ_addDisk(*args, **kwargs):
    """addDisk(xc, yc, zc, rx, ry, tag=-1) -> tag"""
    return _invoke("model_occ_add_disk", *args, **kwargs)

model_occ.addDisk = _model_occ_addDisk

def _model_occ_addRectangle(*args, **kwargs):
    """addRectangle(x, y, z, dx, dy, tag=-1) -> tag
    Create a planar rectangle surface in the XY plane (z offset by z).
    Matches gmsh.model.occ.addRectangle.
    """
    return _invoke("model_occ_add_rectangle", *args, **kwargs)

model_occ.addRectangle = _model_occ_addRectangle

def _model_occ_addCylinder(*args, **kwargs):
    """rmshModelOccAddCylinder(x, y, z, dx, dy, dz, r, tag, ierr)"""
    return _invoke("model_occ_add_cylinder", *args, **kwargs)

model_occ.addCylinder = _model_occ_addCylinder

def _model_occ_addSphere(*args, **kwargs):
    """rmshModelOccAddSphere(x, y, z, r, tag, ierr)"""
    return _invoke("model_occ_add_sphere", *args, **kwargs)

model_occ.addSphere = _model_occ_addSphere

def _model_occ_cut(*args, **kwargs):
    """cut(objectDimTags, toolDimTags, tag=-1, removeObject=True, removeTool=True) -> (outDimTags, outDimTagsMap)"""
    return _invoke("model_occ_cut", *args, **kwargs)

model_occ.cut = _model_occ_cut

def _model_occ_fragment(*args, **kwargs):
    """fragment(objectDimTags, toolDimTags, tag=-1, removeObject=True, removeTool=True) -> (outDimTags, outDimTagsMap)"""
    return _invoke("model_occ_fragment", *args, **kwargs)

model_occ.fragment = _model_occ_fragment

def _model_occ_intersect(*args, **kwargs):
    """intersect(objectDimTags, toolDimTags, tag=-1, removeObject=True, removeTool=True) -> (outDimTags, outDimTagsMap)"""
    return _invoke("model_occ_intersect", *args, **kwargs)

model_occ.intersect = _model_occ_intersect

def _model_occ_fuse(*args, **kwargs):
    """fuse(objectDimTags, toolDimTags, tag=-1, removeObject=True, removeTool=True) -> (outDimTags, outDimTagsMap)"""
    return _invoke("model_occ_fuse", *args, **kwargs)

model_occ.fuse = _model_occ_fuse

def _model_occ_copy(*args, **kwargs):
    """copy(dimTags) -> outDimTags"""
    return _invoke("model_occ_copy", *args, **kwargs)

model_occ.copy = _model_occ_copy

def _model_occ_remove(*args, **kwargs):
    """remove(dimTags, recursive=False)"""
    return _invoke("model_occ_remove", *args, **kwargs)

model_occ.remove = _model_occ_remove

def _model_occ_translate(*args, **kwargs):
    """translate(dimTags, dx, dy, dz)"""
    return _invoke("model_occ_translate", *args, **kwargs)

model_occ.translate = _model_occ_translate

def _model_occ_rotate(*args, **kwargs):
    """rotate(dimTags, x, y, z, ax, ay, az, angle)"""
    return _invoke("model_occ_rotate", *args, **kwargs)

model_occ.rotate = _model_occ_rotate

def _model_occ_dilate(*args, **kwargs):
    """dilate(dimTags, x, y, z, a, b, c)"""
    return _invoke("model_occ_dilate", *args, **kwargs)

model_occ.dilate = _model_occ_dilate

def _model_occ_mirror(*args, **kwargs):
    """mirror(dimTags, a, b, c, d)"""
    return _invoke("model_occ_mirror", *args, **kwargs)

model_occ.mirror = _model_occ_mirror

def _model_occ_importShapes(*args, **kwargs):
    """importShapes(fileName, highestDimOnly=True, format='') -> outDimTags"""
    return _invoke("model_occ_import_shapes", *args, **kwargs)

model_occ.importShapes = _model_occ_importShapes

def _model_occ_synchronize(*args, **kwargs):
    """rmshModelOccSynchronize(ierr)"""
    return _invoke("model_occ_synchronize", *args, **kwargs)

model_occ.synchronize = _model_occ_synchronize

def _model_occ_addCone(*args, **kwargs):
    """addCone(x, y, z, dx, dy, dz, r1, r2, tag=-1, angle=2*pi) -> tag"""
    return _invoke("model_occ_add_cone", *args, **kwargs)

model_occ.addCone = _model_occ_addCone

def _model_occ_addTorus(*args, **kwargs):
    """addTorus(x, y, z, r1, r2, tag=-1, angle=2*pi, zAxis=[]) -> tag
    Also supports legacy rmsh signature:
    addTorus(x, y, z, dx, dy, dz, r1, r2, tag=-1)."""
    return _invoke("model_occ_add_torus", *args, **kwargs)

model_occ.addTorus = _model_occ_addTorus

def _model_occ_fillet(*args, **kwargs):
    """fillet(tag, curveTags, radii) -> new_tag"""
    return _invoke("model_occ_fillet", *args, **kwargs)

model_occ.fillet = _model_occ_fillet

def _model_occ_chamfer(*args, **kwargs):
    """chamfer(tag, curveTags, distances) -> new_tag"""
    return _invoke("model_occ_chamfer", *args, **kwargs)

model_occ.chamfer = _model_occ_chamfer

def _model_occ_healShapes(*args, **kwargs):
    """healShapes(tag, tolerance=1e-8) -> report_dict"""
    return _invoke("model_occ_heal_shapes", *args, **kwargs)

model_occ.healShapes = _model_occ_healShapes

def _model_occ_getMass(*args, **kwargs):
    """getMass(tag) -> volume"""
    return _invoke("model_occ_get_mass", *args, **kwargs)

model_occ.getMass = _model_occ_getMass

def _model_occ_getProperties(*args, **kwargs):
    """getProperties(tag) -> (volume, area, cx, cy, cz)"""
    return _invoke("model_occ_get_properties", *args, **kwargs)

model_occ.getProperties = _model_occ_getProperties

def _model_occ_extrude(*args, **kwargs):
    """extrude(tag, face_idx, dx, dy, dz, distance) -> new_tag"""
    return _invoke("model_occ_extrude", *args, **kwargs)

model_occ.extrude = _model_occ_extrude

def _model_occ_revolve(*args, **kwargs):
    """revolve(tag, face_idx, ax, ay, az, dx, dy, dz, angle) -> new_tag"""
    return _invoke("model_occ_revolve", *args, **kwargs)

model_occ.revolve = _model_occ_revolve

def _option_getColor(*args, **kwargs):
    """rmshOptionGetColor(name, r, g, b, a, ierr)"""
    return _invoke("option_get_color", *args, **kwargs)

option.getColor = _option_getColor

def _option_getNumber(*args, **kwargs):
    """rmshOptionGetNumber(name, value, ierr)"""
    return _invoke("option_get_number", *args, **kwargs)

option.getNumber = _option_getNumber

def _option_getString(*args, **kwargs):
    """rmshOptionGetString(name, value, ierr)"""
    return _invoke("option_get_string", *args, **kwargs)

option.getString = _option_getString

def _option_setColor(*args, **kwargs):
    """rmshOptionSetColor(name, r, g, b, a, ierr)"""
    return _invoke("option_set_color", *args, **kwargs)

option.setColor = _option_setColor

def _option_setNumber(*args, **kwargs):
    """rmshOptionSetNumber(name, value, ierr)"""
    return _invoke("option_set_number", *args, **kwargs)

option.setNumber = _option_setNumber

def _option_setString(*args, **kwargs):
    """rmshOptionSetString(name, value, ierr)"""
    return _invoke("option_set_string", *args, **kwargs)

option.setString = _option_setString

def _option_restoreDefaults(*args, **kwargs):
    """restoreDefaults() - clear all option overrides"""
    return _invoke("option_restore_defaults", *args, **kwargs)

option.restoreDefaults = _option_restoreDefaults

def _plugin_run(*args, **kwargs):
    """rmshPluginRun(name, ierr)"""
    return _invoke("plugin_run", *args, **kwargs)

plugin.run = _plugin_run

def _plugin_setNumber(*args, **kwargs):
    """rmshPluginSetNumber(name, option, value, ierr)"""
    return _invoke("plugin_set_number", *args, **kwargs)

plugin.setNumber = _plugin_setNumber

def _plugin_setString(*args, **kwargs):
    """rmshPluginSetString(name, option, value, ierr)"""
    return _invoke("plugin_set_string", *args, **kwargs)

plugin.setString = _plugin_setString

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
    "plugin",
    "gui",
]

# Public namespace aliases
model = model
option = option
logger = logger
plugin = plugin
gui = gui
