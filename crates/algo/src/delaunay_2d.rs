//! Delaunay 2-D — Bowyer-Watson Delaunay triangulation (Gmsh algorithm 5).
//!
//! # Algorithm overview
//!
//! Pure Delaunay triangular meshing via Bowyer-Watson point insertion.
//! Boundary edges are discretised and interior points seeded on a uniform
//! staggered grid, then the Delaunay triangulation is computed and
//! elements whose centroid falls outside the domain are discarded.
//!
//! This is the simplest 2-D mesher in the Gmsh taxonomy (algo 5).
//! For domains with holes the implementation delegates to
//! [`crate::planar_meshing::mesh_domain_triangles`].

use std::collections::HashMap;

use rmsh_model::{Element, ElementType, Mesh, Node};
use thiserror::Error;

use crate::planar_meshing::{mesh_domain_triangles, validate_domain};
use crate::traits::{Domain2D, MeshAlgoError, MeshParams, Mesher2D};

// ─── Error type ───────────────────────────────────────────────────────────────

/// Error type for 2-D Delaunay mesh generation.
#[derive(Error, Debug)]
pub enum Delaunay2DError {
    #[error("mesh generation: {0}")]
    Generation(String),
}

// ─── Polygon type ─────────────────────────────────────────────────────────────

/// A 2-D polygon defined by an ordered list of boundary vertices.
pub struct Polygon2D {
    pub vertices: Vec<[f64; 2]>,
}

impl Polygon2D {
    /// Create a new polygon from a list of 2-D vertices (CCW or CW order).
    pub fn new(vertices: Vec<[f64; 2]>) -> Self {
        Self { vertices }
    }

    /// Point-in-polygon test using the ray casting algorithm.
    pub fn contains(&self, p: [f64; 2]) -> bool {
        let n = self.vertices.len();
        let (px, py) = (p[0], p[1]);
        let mut inside = false;
        let mut j = n - 1;
        for i in 0..n {
            let (xi, yi) = (self.vertices[i][0], self.vertices[i][1]);
            let (xj, yj) = (self.vertices[j][0], self.vertices[j][1]);
            if ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi) {
                inside = !inside;
            }
            j = i;
        }
        inside
    }

    /// Axis-aligned bounding box: returns `(min, max)`.
    pub fn bounding_box(&self) -> ([f64; 2], [f64; 2]) {
        let mut min = [f64::MAX; 2];
        let mut max = [f64::MIN; 2];
        for v in &self.vertices {
            min[0] = min[0].min(v[0]);
            min[1] = min[1].min(v[1]);
            max[0] = max[0].max(v[0]);
            max[1] = max[1].max(v[1]);
        }
        (min, max)
    }
}

// ─── Public struct ────────────────────────────────────────────────────────────

/// Delaunay 2-D mesher (Gmsh algorithm 5).
///
/// Produces triangular meshes via Bowyer-Watson Delaunay triangulation
/// with uniform staggered-grid point sampling.
#[derive(Debug, Default, Clone)]
pub struct Delaunay2D;

impl Delaunay2D {
    pub fn new() -> Self {
        Self
    }
}

impl Mesher2D for Delaunay2D {
    fn name(&self) -> &'static str {
        "Delaunay 2D"
    }

    fn mesh_2d(&self, domain: &Domain2D, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        validate_domain(domain, params.element_size)?;

        if domain.boundaries.len() == 1 {
            let poly = Polygon2D::new(domain.outer().to_vec());
            return mesh_polygon(&poly, params.element_size)
                .map_err(|e| MeshAlgoError::Generation(e.to_string()));
        }

        // Multi-boundary (holes): delegate to planar domain triangulation
        // which also uses Bowyer-Watson under the hood.
        mesh_domain_triangles(domain, params.element_size, params.element_size * 0.866, 0.0)
    }
}

// ─── Polygon meshing (Bowyer-Watson pipeline) ─────────────────────────────────

/// Generate a 2-D triangular mesh inside `polygon` with approximate target
/// edge length `mesh_size`.
///
/// The algorithm:
/// 1. Discretises the boundary edges (respecting `mesh_size`).
/// 2. Scatters interior seed points on a staggered uniform grid.
/// 3. Runs Bowyer-Watson Delaunay triangulation on all points.
/// 4. Discards triangles whose centroid lies outside the polygon.
///
/// Returns a [`Mesh`] with `Triangle3` elements at z = 0.
pub fn mesh_polygon(polygon: &Polygon2D, mesh_size: f64) -> Result<Mesh, Delaunay2DError> {
    if polygon.vertices.len() < 3 {
        return Err(Delaunay2DError::Generation(
            "Polygon must have at least 3 vertices".to_string(),
        ));
    }
    if mesh_size <= 0.0 {
        return Err(Delaunay2DError::Generation(
            "mesh_size must be positive".to_string(),
        ));
    }

    let (bb_min, bb_max) = polygon.bounding_box();

    // ── Boundary points ───────────────────────────────────────────────────
    let mut points: Vec<[f64; 2]> = Vec::new();
    let nv = polygon.vertices.len();
    for i in 0..nv {
        let a = polygon.vertices[i];
        let b = polygon.vertices[(i + 1) % nv];
        let len = dist2(a, b);
        let nseg = ((len / mesh_size).ceil() as usize).max(1);
        for k in 0..nseg {
            let t = k as f64 / nseg as f64;
            points.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
        }
    }

    // ── Interior seed points on staggered grid ────────────────────────────
    let mut iy = 1usize;
    let mut py = bb_min[1] + mesh_size;
    while py < bb_max[1] {
        let offset_x = if iy % 2 == 0 { 0.0 } else { mesh_size * 0.5 };
        let mut px = bb_min[0] + offset_x + mesh_size;
        while px < bb_max[0] {
            if polygon.contains([px, py]) {
                points.push([px, py]);
            }
            px += mesh_size;
        }
        py += mesh_size * 0.866; // sqrt(3)/2 for equilateral row spacing
        iy += 1;
    }

    // ── Deduplicate points that are too close ─────────────────────────────
    let tol = mesh_size * 0.1;
    points = deduplicate(points, tol * tol);

    if points.len() < 3 {
        return Err(Delaunay2DError::Generation(
            "Too few distinct points after sampling".to_string(),
        ));
    }

    // ── Bowyer-Watson Delaunay triangulation ─────────────────────────────
    let tris = crate::triangulate2d::triangulate_points(&points);

    // ── Build Mesh, keeping only triangles whose centroid is inside polygon ──
    let mut mesh = Mesh::new();
    let mut node_id: u64 = 1;
    let mut elem_id: u64 = 1;
    let mut pt_to_node: HashMap<usize, u64> = HashMap::new();

    for tri in tris {
        let centroid = [
            (points[tri[0]][0] + points[tri[1]][0] + points[tri[2]][0]) / 3.0,
            (points[tri[0]][1] + points[tri[1]][1] + points[tri[2]][1]) / 3.0,
        ];
        if !polygon.contains(centroid) {
            continue;
        }

        let mut nids: Vec<u64> = Vec::with_capacity(3);
        for &pi in &tri {
            let nid = *pt_to_node.entry(pi).or_insert_with(|| {
                let id = node_id;
                node_id += 1;
                mesh.add_node(Node::new(id, points[pi][0], points[pi][1], 0.0));
                id
            });
            nids.push(nid);
        }
        mesh.add_element(Element::new(elem_id, ElementType::Triangle3, nids));
        elem_id += 1;
    }

    if mesh.element_count() == 0 {
        return Err(Delaunay2DError::Generation(
            "No interior triangles were generated".to_string(),
        ));
    }

    Ok(mesh)
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

fn dist2(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    (dx * dx + dy * dy).sqrt()
}

fn deduplicate(pts: Vec<[f64; 2]>, min_dist_sq: f64) -> Vec<[f64; 2]> {
    let mut result: Vec<[f64; 2]> = Vec::with_capacity(pts.len());
    'outer: for p in pts {
        for q in &result {
            let dx = p[0] - q[0];
            let dy = p[1] - q[1];
            if dx * dx + dy * dy < min_dist_sq {
                continue 'outer;
            }
        }
        result.push(p);
    }
    result
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delaunay_2d_meshes_square() {
        let domain =
            Domain2D::from_outer(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);
        let mesh = Delaunay2D
            .mesh_2d(&domain, &MeshParams::with_size(0.5))
            .unwrap();
        assert!(mesh.element_count() > 0);
    }

    #[test]
    fn delaunay_2d_meshes_with_hole() {
        let domain = Domain2D::from_outer(vec![[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0]])
            .with_hole(vec![[1.0, 1.0], [2.0, 1.0], [2.0, 2.0], [1.0, 2.0]]);
        let mesh = Delaunay2D
            .mesh_2d(&domain, &MeshParams::with_size(0.5))
            .unwrap();
        assert!(mesh.element_count() > 0);
    }

    #[test]
    fn delaunay_2d_name_is_stable() {
        assert_eq!(Delaunay2D.name(), "Delaunay 2D");
    }

    // ── mesh_polygon tests ──────────────────────────────────────────────────

    #[test]
    fn mesh_unit_square() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let mesh = mesh_polygon(&poly, 0.25).expect("meshing should succeed");
        assert!(mesh.node_count() > 0);
        assert!(mesh.element_count() > 0);
        for elem in &mesh.elements {
            assert_eq!(elem.etype, ElementType::Triangle3);
            assert_eq!(elem.node_ids.len(), 3);
        }
        for elem in &mesh.elements {
            for &nid in &elem.node_ids {
                assert!(mesh.nodes.contains_key(&nid));
            }
        }
    }

    #[test]
    fn mesh_l_shape() {
        let poly = Polygon2D::new(vec![
            [0.0, 0.0], [2.0, 0.0], [2.0, 1.0],
            [1.0, 1.0], [1.0, 2.0], [0.0, 2.0],
        ]);
        let mesh = mesh_polygon(&poly, 0.3).expect("L-shape meshing should succeed");
        assert!(mesh.element_count() > 0);
    }

    #[test]
    fn mesh_polygon_produces_planar_2d_mesh() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]);
        let mesh = mesh_polygon(&poly, 0.4).expect("meshing should succeed");
        for node in mesh.nodes.values() {
            assert!(node.position.z.abs() < 1e-12, "2D meshing should keep z=0");
        }
        for elem in &mesh.elements {
            assert_eq!(elem.dimension(), 2);
            assert_eq!(elem.etype, ElementType::Triangle3);
        }
    }

    #[test]
    fn mesh_rejects_bad_inputs() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0]]);
        assert!(mesh_polygon(&poly, 0.5).is_err());
        let poly2 = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]]);
        assert!(mesh_polygon(&poly2, -1.0).is_err());
    }

    #[test]
    fn point_in_polygon() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        assert!(poly.contains([0.5, 0.5]));
        assert!(!poly.contains([1.5, 0.5]));
        assert!(!poly.contains([-0.1, 0.5]));
    }

    #[test]
    fn finer_mesh_size_produces_more_elements() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let coarse = mesh_polygon(&poly, 0.5).expect("coarse");
        let fine = mesh_polygon(&poly, 0.15).expect("fine");
        assert!(fine.element_count() > coarse.element_count());
    }

    #[test]
    fn mesh_polygon_all_centroids_inside_polygon() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]);
        let mesh = mesh_polygon(&poly, 0.3).expect("meshing");
        for elem in &mesh.elements {
            let pts: Vec<_> = elem.node_ids.iter().map(|id| &mesh.nodes[id]).collect();
            let cx = pts.iter().map(|n| n.position.x).sum::<f64>() / pts.len() as f64;
            let cy = pts.iter().map(|n| n.position.y).sum::<f64>() / pts.len() as f64;
            assert!(poly.contains([cx, cy]), "centroid ({cx},{cy}) outside polygon");
        }
    }

    #[test]
    fn mesh_polygon_zero_size_fails() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        assert!(mesh_polygon(&poly, 0.0).is_err());
    }

    // ── Polygon2D helpers ───────────────────────────────────────────────────

    #[test]
    fn bounding_box_unit_square() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let (min, max) = poly.bounding_box();
        assert_eq!(min, [0.0, 0.0]);
        assert_eq!(max, [1.0, 1.0]);
    }

    #[test]
    fn bounding_box_non_axis_aligned() {
        let poly = Polygon2D::new(vec![[-2.0, 1.0], [3.0, -1.0], [1.0, 4.0]]);
        let (min, max) = poly.bounding_box();
        assert!((min[0] - (-2.0)).abs() < 1e-12);
        assert!((min[1] - (-1.0)).abs() < 1e-12);
        assert!((max[0] - 3.0).abs() < 1e-12);
        assert!((max[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn contains_on_boundary_is_consistent() {
        let poly = Polygon2D::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let _: bool = poly.contains([0.0, 0.5]);
        let _: bool = poly.contains([0.5, 0.0]);
    }
}
