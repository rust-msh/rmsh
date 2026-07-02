use std::collections::HashMap;

use rmsh_model::{Element, ElementType, Mesh, Node};

use crate::traits::{Domain2D, MeshAlgoError};
use crate::delaunay_2d::Polygon2D;
use crate::triangulate2d::triangulate_points;

pub(crate) fn validate_domain(domain: &Domain2D, mesh_size: f64) -> Result<(), MeshAlgoError> {
    if !mesh_size.is_finite() || mesh_size <= 0.0 {
        return Err(MeshAlgoError::InvalidInput(
            "mesh size must be a positive finite value".to_string(),
        ));
    }
    if domain.boundaries.is_empty() || domain.boundaries[0].len() < 3 {
        return Err(MeshAlgoError::InvalidInput(
            "2D domain must have an outer boundary with at least 3 vertices".to_string(),
        ));
    }
    if domain.boundaries.iter().any(|b| b.len() < 3) {
        return Err(MeshAlgoError::InvalidInput(
            "all hole boundaries must have at least 3 vertices".to_string(),
        ));
    }
    Ok(())
}

pub(crate) fn point_in_domain(domain: &Domain2D, point: [f64; 2]) -> bool {
    let outer = Polygon2D::new(domain.boundaries[0].clone());
    if !outer.contains(point) {
        return false;
    }
    !domain
        .boundaries
        .iter()
        .skip(1)
        .any(|hole| Polygon2D::new(hole.clone()).contains(point))
}

pub(crate) fn sample_domain_points(
    domain: &Domain2D,
    spacing_x: f64,
    spacing_y: f64,
    row_offset_factor: f64,
) -> Result<Vec<[f64; 2]>, MeshAlgoError> {
    validate_domain(domain, spacing_x.min(spacing_y))?;
    let outer = Polygon2D::new(domain.boundaries[0].clone());
    let (bb_min, bb_max) = outer.bounding_box();

    let mut points: Vec<[f64; 2]> = Vec::new();

    for boundary in &domain.boundaries {
        for i in 0..boundary.len() {
            let a = boundary[i];
            let b = boundary[(i + 1) % boundary.len()];
            let len = distance(a, b);
            let nseg = ((len / spacing_x.min(spacing_y)).ceil() as usize).max(1);
            for k in 0..nseg {
                let t = k as f64 / nseg as f64;
                points.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
            }
        }
    }

    let mut row = 0usize;
    let mut py = bb_min[1] + spacing_y;
    while py < bb_max[1] {
        let offset_x = if row % 2 == 0 {
            0.0
        } else {
            spacing_x * row_offset_factor
        };
        let mut px = bb_min[0] + spacing_x + offset_x;
        while px < bb_max[0] {
            let p = [px, py];
            if point_in_domain(domain, p) {
                points.push(p);
            }
            px += spacing_x;
        }
        py += spacing_y;
        row += 1;
    }

    let tol = spacing_x.min(spacing_y) * 0.1;
    Ok(deduplicate(points, tol * tol))
}

pub(crate) fn mesh_domain_triangles(
    domain: &Domain2D,
    spacing_x: f64,
    spacing_y: f64,
    row_offset_factor: f64,
) -> Result<Mesh, MeshAlgoError> {
    let points = sample_domain_points(domain, spacing_x, spacing_y, row_offset_factor)?;
    if points.len() < 3 {
        return Err(MeshAlgoError::Generation(
            "too few distinct points after domain sampling".to_string(),
        ));
    }

    let tris = triangulate_points(&points);
    let mut mesh = Mesh::new();
    let mut node_id = 1u64;
    let mut elem_id = 1u64;
    let mut point_to_node: HashMap<usize, u64> = HashMap::new();

    for tri in tris {
        let centroid = [
            (points[tri[0]][0] + points[tri[1]][0] + points[tri[2]][0]) / 3.0,
            (points[tri[0]][1] + points[tri[1]][1] + points[tri[2]][1]) / 3.0,
        ];
        if !point_in_domain(domain, centroid) {
            continue;
        }

        let mut nids = Vec::with_capacity(3);
        for &pi in &tri {
            let nid = *point_to_node.entry(pi).or_insert_with(|| {
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
        return Err(MeshAlgoError::Generation(
            "no interior triangles were generated".to_string(),
        ));
    }

    Ok(mesh)
}

pub(crate) fn is_axis_aligned_rectangle(poly: &[[f64; 2]]) -> Option<([f64; 2], [f64; 2])> {
    if poly.len() != 4 {
        return None;
    }
    let min_x = poly.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
    let max_x = poly.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
    let min_y = poly.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
    let max_y = poly.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max);
    let expected = [
        [min_x, min_y],
        [max_x, min_y],
        [max_x, max_y],
        [min_x, max_y],
    ];
    if poly.iter().all(|p| expected.contains(p)) {
        Some(([min_x, min_y], [max_x, max_y]))
    } else {
        None
    }
}

pub(crate) fn structured_quad_mesh_rectangle(min: [f64; 2], max: [f64; 2], mesh_size: f64) -> Mesh {
    let mut mesh = Mesh::new();
    let width = (max[0] - min[0]).max(mesh_size);
    let height = (max[1] - min[1]).max(mesh_size);
    let nx = ((width / mesh_size).ceil() as usize).max(1);
    let ny = ((height / mesh_size).ceil() as usize).max(1);
    let dx = width / nx as f64;
    let dy = height / ny as f64;

    let mut next_node_id = 1u64;
    let mut ids = vec![vec![0u64; nx + 1]; ny + 1];
    for iy in 0..=ny {
        for ix in 0..=nx {
            let x = min[0] + ix as f64 * dx;
            let y = min[1] + iy as f64 * dy;
            ids[iy][ix] = next_node_id;
            mesh.add_node(Node::new(next_node_id, x, y, 0.0));
            next_node_id += 1;
        }
    }

    let mut next_elem_id = 1u64;
    for iy in 0..ny {
        for ix in 0..nx {
            mesh.add_element(Element::new(
                next_elem_id,
                ElementType::Quad4,
                vec![
                    ids[iy][ix],
                    ids[iy][ix + 1],
                    ids[iy + 1][ix + 1],
                    ids[iy + 1][ix],
                ],
            ));
            next_elem_id += 1;
        }
    }

    mesh
}

fn distance(a: [f64; 2], b: [f64; 2]) -> f64 {
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
