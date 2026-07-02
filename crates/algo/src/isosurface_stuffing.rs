//! Simplified isosurface stuffing — BCC-lattice initial tetrahedral mesh.
//!
//! Builds a Body-Centered Cubic (BCC) lattice covering the surface bounding
//! box, classifies each lattice point as inside/outside via ray casting, then
//! generates tetrahedra for fully-interior cells using a 12-tet template.
//!
//! This replaces `CentroidStarMesher3D` for domains where a single centroid
//! produces degenerate tets (non-star-shaped domains).
//!
//! Reference: Labelle & Shewchuk, "Isosurface Stuffing", SIGGRAPH 2007.

use rmsh_model::{Element, ElementType, Mesh, Node};
use thiserror::Error;

use crate::traits::{MeshAlgoError, MeshParams, Mesher3D};

#[derive(Error, Debug)]
pub enum BccError {
    #[error("BCC meshing failed: {0}")]
    Generation(String),
}

impl From<BccError> for MeshAlgoError {
    fn from(e: BccError) -> Self {
        MeshAlgoError::Generation(e.to_string())
    }
}

/// BCC lattice mesher — produces a tetrahedral mesh from a closed surface.
///
/// Resolution is determined by the BCC cell size, computed from
/// `element_size` × a density factor.
#[derive(Debug, Clone)]
pub struct BccMesher {
    /// BCC cell edge length (cube size).  Smaller = finer mesh.
    pub cell_size: f64,
}

impl Default for BccMesher {
    fn default() -> Self {
        Self { cell_size: 0.3 }
    }
}

impl Mesher3D for BccMesher {
    fn name(&self) -> &'static str {
        "BCC Isosurface Stuffing"
    }

    fn mesh_3d(&self, surface: &Mesh, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        if !params.element_size.is_finite() || params.element_size <= 0.0 {
            return Err(MeshAlgoError::InvalidInput(
                "element_size must be positive".into(),
            ));
        }
        let cell = self.cell_size.max(params.element_size * 0.5);
        bcc_tetrahedralize(surface, cell).map_err(Into::into)
    }
}

// ─── BCC tetrahedralization ────────────────────────────────────────────────

/// Maximum number of BCC cells per axis (memory cap).
const MAX_CELLS: usize = 48;

/// Generate a tetrahedral mesh by stuffing a BCC lattice inside `surface`.
fn bcc_tetrahedralize(surface: &Mesh, cell_size: f64) -> Result<Mesh, BccError> {
    // 1. Compute AABB
    let mut bmin = [f64::MAX; 3];
    let mut bmax = [f64::MIN; 3];
    let mut has_node = false;
    for n in surface.nodes.values() {
        has_node = true;
        bmin[0] = bmin[0].min(n.position.x);
        bmin[1] = bmin[1].min(n.position.y);
        bmin[2] = bmin[2].min(n.position.z);
        bmax[0] = bmax[0].max(n.position.x);
        bmax[1] = bmax[1].max(n.position.y);
        bmax[2] = bmax[2].max(n.position.z);
    }
    if !has_node {
        return Err(BccError::Generation("empty surface mesh".into()));
    }
    // Expand AABB by one cell on each side
    for i in 0..3 {
        let pad = cell_size;
        bmin[i] -= pad;
        bmax[i] += pad;
    }

    let nx = ((bmax[0] - bmin[0]) / cell_size).ceil() as usize;
    let ny = ((bmax[1] - bmin[1]) / cell_size).ceil() as usize;
    let nz = ((bmax[2] - bmin[2]) / cell_size).ceil() as usize;
    if nx > MAX_CELLS || ny > MAX_CELLS || nz > MAX_CELLS {
        return Err(BccError::Generation(format!(
            "grid too large: {nx}×{ny}×{nz} (max {MAX_CELLS} per axis); increase cell_size"
        )));
    }
    let total_cells = nx * ny * nz;
    if total_cells == 0 {
        return Err(BccError::Generation("zero-volume bounding box".into()));
    }

    // 2. Classify BCC lattice points
    // Corners: (ix, iy, iz) for ix=0..nx, etc.
    let n_corners = (nx + 1) * (ny + 1) * (nz + 1);
    // Centers: (cx, cy, cz) for cx=0..nx, etc. (one center per cell)
    let n_centers = total_cells;
    let total_pts = n_corners + n_centers;
    let mut pt_inside: Vec<bool> = Vec::with_capacity(total_pts);

    // Allocate lazily: compute and store
    let mut inside_count = 0usize;

    // Corners first
    for iz in 0..=nz {
        for iy in 0..=ny {
            for ix in 0..=nx {
                let p = [
                    bmin[0] + ix as f64 * cell_size,
                    bmin[1] + iy as f64 * cell_size,
                    bmin[2] + iz as f64 * cell_size,
                ];
                let inside = point_inside_surface(p, surface);
                pt_inside.push(inside);
                if inside { inside_count += 1; }
            }
        }
    }
    // Centers
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let p = [
                    bmin[0] + (ix as f64 + 0.5) * cell_size,
                    bmin[1] + (iy as f64 + 0.5) * cell_size,
                    bmin[2] + (iz as f64 + 0.5) * cell_size,
                ];
                let inside = point_inside_surface(p, surface);
                pt_inside.push(inside);
                if inside { inside_count += 1; }
            }
        }
    }

    if inside_count == 0 {
        return Err(BccError::Generation("no lattice points inside surface".into()));
    }

    // Helpers: corner index and center index for (ix, iy, iz)
    let cidx = |ix: usize, iy: usize, iz: usize| -> usize {
        iz * (ny + 1) * (nx + 1) + iy * (nx + 1) + ix
    };
    let ctr_idx = |ix: usize, iy: usize, iz: usize| -> usize {
        n_corners + iz * ny * nx + iy * nx + ix
    };

    // 3. Generate tets for fully-interior cells
    let mut out = Mesh::new();
    let mut node_id: u64 = 1;
    let mut elem_id: u64 = 1;
    // Map BCC point index → mesh node id
    let mut pt_to_node: Vec<u64> = vec![0; total_pts];
    let mut added_count = 0usize;

    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let c0 = cidx(ix, iy, iz);
                let c1 = cidx(ix + 1, iy, iz);
                let c2 = cidx(ix + 1, iy + 1, iz);
                let c3 = cidx(ix, iy + 1, iz);
                let c4 = cidx(ix, iy, iz + 1);
                let c5 = cidx(ix + 1, iy, iz + 1);
                let c6 = cidx(ix + 1, iy + 1, iz + 1);
                let c7 = cidx(ix, iy + 1, iz + 1);
                let center = ctr_idx(ix, iy, iz);

                // All 9 points (8 corners + center) must be inside
                if !pt_inside[c0] || !pt_inside[c1] || !pt_inside[c2] || !pt_inside[c3]
                    || !pt_inside[c4] || !pt_inside[c5] || !pt_inside[c6] || !pt_inside[c7]
                    || !pt_inside[center]
                {
                    continue;
                }

                // Add nodes only when a cell uses them (lazy)
                let add_node = |out: &mut Mesh, pt: [f64; 3], id: u64| {
                    if !out.nodes.contains_key(&id) {
                        out.add_node(Node::new(id, pt[0], pt[1], pt[2]));
                    }
                };

                // Generate node IDs lazily via closure
                let mut next_local_id = node_id;
                let mut get_id = |pi: usize| -> u64 {
                    if pt_to_node[pi] == 0 {
                        pt_to_node[pi] = next_local_id;
                        next_local_id += 1;
                    }
                    pt_to_node[pi]
                };

                let n0 = get_id(c0);
                let n1 = get_id(c1);
                let n2 = get_id(c2);
                let n3 = get_id(c3);
                let n4 = get_id(c4);
                let n5 = get_id(c5);
                let n6 = get_id(c6);
                let n7 = get_id(c7);
                let n8 = get_id(center);
                node_id = next_local_id;

                let ids = [n0, n1, n2, n3, n4, n5, n6, n7, n8];

                // Compute coordinates
                let x0 = bmin[0] + ix as f64 * cell_size;
                let y0 = bmin[1] + iy as f64 * cell_size;
                let z0 = bmin[2] + iz as f64 * cell_size;
                let coords: [[f64; 3]; 9] = [
                    [x0, y0, z0],
                    [x0 + cell_size, y0, z0],
                    [x0 + cell_size, y0 + cell_size, z0],
                    [x0, y0 + cell_size, z0],
                    [x0, y0, z0 + cell_size],
                    [x0 + cell_size, y0, z0 + cell_size],
                    [x0 + cell_size, y0 + cell_size, z0 + cell_size],
                    [x0, y0 + cell_size, z0 + cell_size],
                    [x0 + cell_size * 0.5, y0 + cell_size * 0.5, z0 + cell_size * 0.5],
                ];
                for i in 0..9 {
                    add_node(&mut out, coords[i], ids[i]);
                }

                // 12-tet template for a BCC cell.
                // Each tet = (face_triangle, center).
                // Consistent diagonal for face splits to avoid gaps.
                // Face order: -X, +X, -Y, +Y, -Z, +Z
                // Each face split into two triangles.
                #[rustfmt::skip]
                let tets: [[usize; 4]; 12] = [
                    // -X face (0,3,7,4) → triangles (0,3,7) & (0,7,4)
                    [0, 3, 7, 8], [0, 7, 4, 8],
                    // +X face (1,2,6,5) → triangles (1,2,6) & (1,6,5)
                    [1, 2, 6, 8], [1, 6, 5, 8],
                    // -Y face (0,1,5,4) → triangles (0,1,5) & (0,5,4)
                    [0, 1, 5, 8], [0, 5, 4, 8],
                    // +Y face (3,2,6,7) → triangles (3,2,6) & (3,6,7)
                    [3, 2, 6, 8], [3, 6, 7, 8],
                    // -Z face (0,1,2,3) → triangles (0,1,2) & (0,2,3)
                    [0, 1, 2, 8], [0, 2, 3, 8],
                    // +Z face (4,5,6,7) → triangles (4,5,6) & (4,6,7)
                    [4, 5, 6, 8], [4, 6, 7, 8],
                ];

                for tet in &tets {
                    let n0 = ids[tet[0]];
                    let n1 = ids[tet[1]];
                    let n2 = ids[tet[2]];
                    let n3 = ids[tet[3]];
                    out.add_element(Element::new(
                        elem_id, ElementType::Tetrahedron4, vec![n0, n1, n2, n3],
                    ));
                    elem_id += 1;
                    added_count += 1;
                }
            }
        }
    }

    if added_count == 0 {
        Err(BccError::Generation("no interior cells found".into()))
    } else {
        Ok(out)
    }
}

// ─── Ray-cast point-in-polyhedron ──────────────────────────────────────────

/// Returns `true` if `p` is inside the closed surface mesh.
fn point_inside_surface(p: [f64; 3], surface: &Mesh) -> bool {
    let mut intersections = 0usize;
    for elt in &surface.elements {
        if elt.dimension() != 2 || elt.node_ids.len() < 3 { continue; }
        let n = &elt.node_ids;

        let Some(a) = surface.nodes.get(&n[0]).map(|no| [no.position.x, no.position.y, no.position.z]) else { continue };
        let Some(b) = surface.nodes.get(&n[1]).map(|no| [no.position.x, no.position.y, no.position.z]) else { continue };
        let Some(c) = surface.nodes.get(&n[2]).map(|no| [no.position.x, no.position.y, no.position.z]) else { continue };

        // Möller–Trumbore: ray = o + t*dir, dir = [1,0,0]
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        // pvec = dir × e2 = [0, -e2[2], e2[1]] where dir=[1,0,0]
        let pvec = [0.0, -e2[2], e2[1]];
        // det = e1 · (dir × e2) = e1[1]*(-e2[2]) + e1[2]*e2[1]
        let det = e1[1] * (-e2[2]) + e1[2] * e2[1];
        if det.abs() < 1e-20 { continue; }
        let inv_det = 1.0 / det;

        let tvec = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
        let u = (tvec[0] * pvec[0] + tvec[1] * pvec[1] + tvec[2] * pvec[2]) * inv_det;
        if u < 0.0 || u > 1.0 { continue; }

        // qvec = tvec × e1 (but we only need the dot with dir=[1,0,0])
        // qvec × dir = qvec × [1,0,0] → really we want (tvec × e1) · dir
        // = tvec[1]*e1[2] - tvec[2]*e1[1]
        // Wait, for barycentric coords:
        // v = (dir · qvec) * inv_det where qvec = tvec × e1
        // But `dir = [1,0,0]`, so `dir · (tvec × e1)` = (tvec × e1)[0] = tvec[1]*e1[2] - tvec[2]*e1[1]
        let v = (tvec[1] * e1[2] - tvec[2] * e1[1]) * inv_det;
        if v < 0.0 || u + v > 1.0 { continue; }

        // t = (e2 · qvec) * inv_det where qvec = tvec × e1
        // = (e2 · (tvec × e1)) * inv_det
        // = (e2[0]*(tvec[1]*e1[2] - tvec[2]*e1[1]) + e2[1]*(tvec[2]*e1[0] - tvec[0]*e1[2]) + e2[2]*(tvec[0]*e1[1]-tvec[1]*e1[0])) * inv_det
        let t = (e2[0] * (tvec[1] * e1[2] - tvec[2] * e1[1])
               + e2[1] * (tvec[2] * e1[0] - tvec[0] * e1[2])
               + e2[2] * (tvec[0] * e1[1] - tvec[1] * e1[0])) * inv_det;
        if t > 1e-12 { intersections += 1; }
    }
    intersections % 2 == 1
}

// ─── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::MeshParams;

    fn cube_surface() -> Mesh {
        let mut mesh = Mesh::new();
        for (id, x) in [
            (1, [0.0, 0.0, 0.0]), (2, [1.0, 0.0, 0.0]), (3, [1.0, 1.0, 0.0]), (4, [0.0, 1.0, 0.0]),
            (5, [0.0, 0.0, 1.0]), (6, [1.0, 0.0, 1.0]), (7, [1.0, 1.0, 1.0]), (8, [0.0, 1.0, 1.0]),
        ] { mesh.add_node(Node::new(id, x[0], x[1], x[2])); }
        mesh.add_element(Element::new(1, ElementType::Quad4, vec![1, 2, 3, 4]));
        mesh.add_element(Element::new(2, ElementType::Quad4, vec![5, 6, 7, 8]));
        mesh.add_element(Element::new(3, ElementType::Quad4, vec![1, 2, 6, 5]));
        mesh.add_element(Element::new(4, ElementType::Quad4, vec![2, 3, 7, 6]));
        mesh.add_element(Element::new(5, ElementType::Quad4, vec![3, 4, 8, 7]));
        mesh.add_element(Element::new(6, ElementType::Quad4, vec![4, 1, 5, 8]));
        mesh
    }

    #[test]
    fn bcc_meshes_cube() {
        let surface = cube_surface();
        let mesh = BccMesher { cell_size: 0.4 }
            .mesh_3d(&surface, &MeshParams::with_size(0.4))
            .expect("BCC should mesh cube");
        assert!(mesh.elements_by_dimension(3).len() > 0, "should produce tets");
        eprintln!("cube: {} tets, {} nodes", mesh.elements.len(), mesh.nodes.len());
    }

    #[test]
    fn bcc_inside_test() {
        let surface = cube_surface();
        assert!(point_inside_surface([0.5, 0.5, 0.5], &surface));
        assert!(!point_inside_surface([-0.5, 0.5, 0.5], &surface));
        assert!(!point_inside_surface([2.0, 0.5, 0.5], &surface));
    }
}
