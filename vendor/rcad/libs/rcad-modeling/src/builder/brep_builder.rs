//! Free-form BRep construction (OCCT `BRepBuilderAPI` equivalent).
//!
//! These functions provide a low-level API to incrementally build a BRep by
//! appending curves, edges, wires, faces, and solids. The caller is responsible
//! for topological consistency.

use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Surface3, SurfaceEval};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

use crate::builder::BuildError;

/// Adds a new vertex to the BRep and returns its index.
pub fn make_vertex(brep: &mut BRep, point: glam::DVec3) -> usize {
    let idx = brep.vertices.len();
    brep.vertices.push(Vertex { point });
    idx
}

/// Adds a new edge (with associated curve and parameter range) to the BRep.
///
/// Returns the index of the new edge.
pub fn make_edge(
    brep: &mut BRep,
    curve: Curve3,
    t1: f64,
    t2: f64,
    v0: usize,
    v1: usize,
) -> Result<usize, BuildError> {
    if v0 >= brep.vertices.len() {
        return Err(BuildError::InvalidIndex(v0));
    }
    if v1 >= brep.vertices.len() {
        return Err(BuildError::InvalidIndex(v1));
    }

    let edge_idx = brep.edges.len();
    brep.edges.push(Edge { start: v0, end: v1 });

    // Register curve in GeomStore
    let curve_idx = brep.geom.curves.len();
    brep.geom.curves.push(curve);

    // Align parallel GeomStore vecs to edge count
    while brep.geom.edge_curve.len() <= edge_idx {
        brep.geom.edge_curve.push(None);
    }
    while brep.geom.edge_curve_range.len() <= edge_idx {
        brep.geom.edge_curve_range.push(None);
    }
    while brep.geom.edge_degenerated.len() <= edge_idx {
        brep.geom.edge_degenerated.push(false);
    }

    brep.geom.edge_curve[edge_idx] = Some(curve_idx);
    brep.geom.edge_curve_range[edge_idx] = Some([t1, t2]);

    Ok(edge_idx)
}

/// Constructs a `Wire` from a list of `WireEdge`s without modifying the BRep.
pub fn make_wire(edges: Vec<WireEdge>) -> Wire {
    Wire { edges }
}

/// Adds a new face to the BRep (within the first solid's first shell, which is
/// created on demand) and returns the face index.
///
/// The face normal is derived from `surface.normal_at(0.0, 0.0)`.
/// Triangles are left empty; call a triangulation pass afterwards if needed.
pub fn make_face(
    brep: &mut BRep,
    surface: Surface3,
    outer: Wire,
    inner_wires: Vec<Wire>,
) -> Result<usize, BuildError> {
    if outer.edges.is_empty() {
        return Err(BuildError::DegenerateGeometry("outer wire has no edges"));
    }

    let normal = surface.normal_at(0.0, 0.0);

    // Register surface in GeomStore
    let surf_idx = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surface);

    // Create solid/shell structure on demand
    if brep.solids.is_empty() {
        brep.solids.push(Solid {
            shells: vec![Shell { faces: Vec::new() }],
        });
    }
    if brep.solids[0].shells.is_empty() {
        brep.solids[0].shells.push(Shell { faces: Vec::new() });
    }

    let face_idx = brep.solids[0].shells[0].faces.len();

    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: outer,
        inner_wires,
        normal,
        triangles: Vec::new(),
        mesh_dirty: true,
    });

    // Align face_surface vec
    while brep.geom.face_surface.len() <= face_idx {
        brep.geom.face_surface.push(None);
    }
    brep.geom.face_surface[face_idx] = Some(surf_idx);

    Ok(face_idx)
}

/// Appends a new solid (composed of the given shells) to the BRep and returns
/// its index.
pub fn make_solid(brep: &mut BRep, shells: Vec<Shell>) -> usize {
    let idx = brep.solids.len();
    brep.solids.push(Solid { shells });
    idx
}
