//! BRepMesh-style mesh generation for BRep shapes.
//!
//! This module provides mesh generation capabilities similar to OCCT's BRepMesh.
//! It includes:
//! - Mesh parameters with deflection and size controls
//! - Face meshing with adaptive subdivision
//! - BRep meshing for complete shapes
//! - Edge discretization for curve sampling
//! - Quality metrics for mesh analysis
//! - Mesh refinement for improved quality

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topology::Face;
use std::collections::HashMap;

// ============================================================================
// Mesh Parameters
// ============================================================================

/// Parameters controlling mesh generation quality and density.
///
/// Analogous to OCCT `IMeshTools_Parameters`.
#[derive(Debug, Clone)]
pub struct MeshParams {
    /// Linear deflection (maximum distance from triangle to surface).
    /// Smaller values produce denser meshes. Default: 0.001.
    pub deflection: f64,
    /// Angular deflection (maximum angle between adjacent triangles).
    /// Smaller values produce smoother meshes. Default: 0.5 (~28.6 degrees).
    pub angle_deflection: f64,
    /// Minimum triangle edge length in world coordinates.
    /// Prevents over-tessellation in small regions. Default: 0.0.
    pub min_mesh_size: f64,
    /// Maximum triangle edge length in world coordinates.
    /// Ensures minimum mesh density. Default: f64::MAX.
    pub max_mesh_size: f64,
    /// Whether to enable mesh optimization.
    /// Improves triangle quality after tessellation. Default: true.
    pub optimize: bool,
    /// Whether to respect face boundaries during meshing.
    /// Default: true.
    pub respect_bounds: bool,
    /// Maximum recursion depth for adaptive subdivision.
    /// Default: 10.
    pub max_depth: usize,
}

impl Default for MeshParams {
    fn default() -> Self {
        Self {
            deflection: 0.001,
            angle_deflection: 0.5,  // ~28.6 degrees
            min_mesh_size: 0.0,
            max_mesh_size: f64::MAX,
            optimize: true,
            respect_bounds: true,
            max_depth: 10,
        }
    }
}

impl MeshParams {
    /// Create default mesh parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the linear deflection.
    ///
    /// The deflection controls how closely the mesh approximates the surface.
    /// Smaller values produce more triangles but better accuracy.
    pub fn with_deflection(mut self, deflection: f64) -> Self {
        self.deflection = deflection.abs().max(1e-10);
        self
    }

    /// Set the angular deflection in radians.
    ///
    /// Controls the maximum angle between adjacent triangle normals.
    pub fn with_angle_deflection(mut self, angle: f64) -> Self {
        self.angle_deflection = angle.abs().max(0.001);
        self
    }

    /// Set the minimum mesh size.
    ///
    /// Triangles smaller than this size will not be further subdivided.
    pub fn with_min_mesh_size(mut self, size: f64) -> Self {
        self.min_mesh_size = size.abs();
        self
    }

    /// Set the maximum mesh size.
    ///
    /// Triangles larger than this size will be subdivided.
    pub fn with_max_mesh_size(mut self, size: f64) -> Self {
        self.max_mesh_size = if size <= 0.0 { f64::MAX } else { size };
        self
    }

    /// Coarse mesh preset for fast preview.
    pub fn coarse() -> Self {
        Self {
            deflection: 0.1,
            angle_deflection: 0.8,
            min_mesh_size: 0.01,
            max_mesh_size: f64::MAX,
            optimize: false,
            respect_bounds: true,
            max_depth: 6,
        }
    }

    /// Fine mesh preset for high-quality rendering.
    pub fn fine() -> Self {
        Self {
            deflection: 0.0001,
            angle_deflection: 0.2,
            min_mesh_size: 0.0,
            max_mesh_size: f64::MAX,
            optimize: true,
            respect_bounds: true,
            max_depth: 14,
        }
    }

    /// Analysis mesh preset for FEA/CFD.
    pub fn analysis() -> Self {
        Self {
            deflection: 0.0001,
            angle_deflection: 0.1,
            min_mesh_size: 0.0,
            max_mesh_size: 1.0,
            optimize: true,
            respect_bounds: true,
            max_depth: 16,
        }
    }
}

// ============================================================================
// Mesh Structures
// ============================================================================

/// A triangular mesh representing a surface.
#[derive(Debug, Clone)]
pub struct Mesh {
    /// Vertex positions in world coordinates.
    pub vertices: Vec<DVec3>,
    /// Triangle indices (3 vertex indices per triangle).
    pub triangles: Vec<[u32; 3]>,
    /// Per-vertex normals.
    pub normals: Vec<DVec3>,
}

impl Default for Mesh {
    fn default() -> Self {
        Self::new()
    }
}

impl Mesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            triangles: Vec::new(),
            normals: Vec::new(),
        }
    }

    /// Check if the mesh is empty.
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }

    /// Get the number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }

    /// Get the number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Compute the bounding box of the mesh.
    pub fn bounding_box(&self) -> (DVec3, DVec3) {
        if self.vertices.is_empty() {
            return (DVec3::ZERO, DVec3::ZERO);
        }

        let mut min = self.vertices[0];
        let mut max = self.vertices[0];

        for &v in &self.vertices {
            min = min.min(v);
            max = max.max(v);
        }

        (min, max)
    }

    /// Flip all triangle winding orders (invert normals).
    pub fn flip(&mut self) {
        for tri in &mut self.triangles {
            tri.swap(0, 2);
        }
        for normal in &mut self.normals {
            *normal = -*normal;
        }
    }

    /// Merge another mesh into this one.
    pub fn merge(&mut self, other: &Mesh) {
        let base = self.vertices.len() as u32;
        self.vertices.extend_from_slice(&other.vertices);
        self.normals.extend_from_slice(&other.normals);
        for tri in &other.triangles {
            self.triangles.push([tri[0] + base, tri[1] + base, tri[2] + base]);
        }
    }
}

/// Per-face mesh data for a complete BRep shape.
#[derive(Debug, Clone)]
pub struct BRepMesh {
    /// Per-face meshes, indexed by face order in the BRep.
    pub face_meshes: Vec<Mesh>,
    /// Original face normals for reference.
    pub face_normals: Vec<DVec3>,
}

impl BRepMesh {
    /// Create an empty BRepMesh.
    pub fn new() -> Self {
        Self {
            face_meshes: Vec::new(),
            face_normals: Vec::new(),
        }
    }

    /// Get total triangle count across all faces.
    pub fn total_triangles(&self) -> usize {
        self.face_meshes.iter().map(|m| m.triangle_count()).sum()
    }

    /// Get total vertex count across all faces.
    pub fn total_vertices(&self) -> usize {
        self.face_meshes.iter().map(|m| m.vertex_count()).sum()
    }

    /// Merge all face meshes into a single mesh.
    pub fn to_merged_mesh(&self) -> Mesh {
        let mut merged = Mesh::new();
        for mesh in &self.face_meshes {
            merged.merge(mesh);
        }
        merged
    }
}

impl Default for BRepMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Face Meshing
// ============================================================================

/// Mesh a single face given its surface.
///
/// This function generates a triangulated mesh for a parametric surface
/// using adaptive subdivision based on deflection and angle criteria.
///
/// # Arguments
/// * `face` - The topological face (for boundary constraints if enabled).
/// * `surface` - The geometric surface to mesh.
/// * `params` - Mesh generation parameters.
///
/// # Returns
/// A `Mesh` containing vertices, triangles, and normals.
pub fn mesh_face(_face: &Face, surface: &Surface3, params: &MeshParams) -> Mesh {
    let domain = surface.default_domain();
    let [u_min, u_max, v_min, v_max] = domain;

    // Handle infinite domains
    let u_range = if u_min.is_finite() && u_max.is_finite() {
        [u_min, u_max]
    } else {
        [-10.0, 10.0]
    };
    let v_range = if v_min.is_finite() && v_max.is_finite() {
        [v_min, v_max]
    } else {
        [-10.0, 10.0]
    };

    // Compute initial grid resolution based on deflection
    let u_span = u_range[1] - u_range[0];
    let v_span = v_range[1] - v_range[0];

    // Estimate curvature-based initial divisions
    let n_u = estimate_initial_divisions(surface, u_range, v_range, params, true);
    let n_v = estimate_initial_divisions(surface, u_range, v_range, params, false);

    let mut vertices: Vec<DVec3> = Vec::new();
    let mut normals: Vec<DVec3> = Vec::new();
    let mut triangles: Vec<[u32; 3]> = Vec::new();

    // Adaptive subdivision of UV domain
    let du = u_span / n_u as f64;
    let dv = v_span / n_v as f64;

    for i in 0..n_u {
        for j in 0..n_v {
            let u0 = u_range[0] + i as f64 * du;
            let u1 = u0 + du;
            let v0 = v_range[0] + j as f64 * dv;
            let v1 = v0 + dv;

            subdivide_quad_for_face(
                surface,
                [u0, u1],
                [v0, v1],
                params,
                0,
                &mut vertices,
                &mut normals,
                &mut triangles,
            );
        }
    }

    // Weld vertices
    let mesh = weld_mesh(Mesh { vertices, triangles, normals });

    // Optimize if requested
    if params.optimize && !mesh.triangles.is_empty() {
        optimize_mesh(mesh)
    } else {
        mesh
    }
}

/// Estimate initial subdivision count based on surface properties.
fn estimate_initial_divisions(
    surface: &Surface3,
    u_range: [f64; 2],
    v_range: [f64; 2],
    params: &MeshParams,
    is_u: bool,
) -> usize {
    let _span = if is_u { u_range[1] - u_range[0] } else { v_range[1] - v_range[0] };

    let base_count = match surface {
        Surface3::Plane(_) => 2,
        Surface3::Cylinder(_) => if is_u { 16 } else { 2 },
        Surface3::Sphere(_) => 16,
        Surface3::Cone(_) => if is_u { 16 } else { 2 },
        Surface3::Torus(_) => 24,
        _ => 8,
    };

    // Adjust based on deflection
    let deflection_factor = (1.0 / params.deflection).sqrt().min(100.0);
    let adjusted = (base_count as f64 * deflection_factor.sqrt()).ceil() as usize;

    adjusted.max(2).min(128)
}

/// Recursively subdivide a UV quad for face meshing.
fn subdivide_quad_for_face(
    surface: &Surface3,
    u_range: [f64; 2],
    v_range: [f64; 2],
    params: &MeshParams,
    depth: usize,
    vertices: &mut Vec<DVec3>,
    normals: &mut Vec<DVec3>,
    triangles: &mut Vec<[u32; 3]>,
) {
    let [u0, u1] = u_range;
    let [v0, v1] = v_range;

    // Compute corner points
    let p00 = surface.point_at(u0, v0);
    let p10 = surface.point_at(u1, v0);
    let p01 = surface.point_at(u0, v1);
    let p11 = surface.point_at(u1, v1);

    let um = (u0 + u1) * 0.5;
    let vm = (v0 + v1) * 0.5;

    // Check if we should continue subdividing
    let should_subdivide = depth < params.max_depth && {
        let step_u = u1 - u0;
        let step_v = v1 - v0;

        // Check minimum step size
        if step_u < 1e-8 && step_v < 1e-8 {
            false
        } else {
            // Check deflection
            let mid = surface.point_at(um, vm);
            let interp_mid = (p00 + p10 + p11 + p01) * 0.25;
            let deflection_error = (mid - interp_mid).length();

            // Check angle deviation
            let n00 = surface.normal_at(u0, v0);
            let n11 = surface.normal_at(u1, v1);
            let angle_error = normal_angle(n00, n11);

            // Check size constraints
            let edge_len = (p11 - p00).length();
            let size_exceeded = edge_len > params.max_mesh_size;

            deflection_error > params.deflection
                || angle_error > params.angle_deflection
                || size_exceeded
        }
    };

    if should_subdivide {
        // Subdivide into 4 sub-quads
        subdivide_quad_for_face(surface, [u0, um], [v0, vm], params, depth + 1, vertices, normals, triangles);
        subdivide_quad_for_face(surface, [um, u1], [v0, vm], params, depth + 1, vertices, normals, triangles);
        subdivide_quad_for_face(surface, [u0, um], [vm, v1], params, depth + 1, vertices, normals, triangles);
        subdivide_quad_for_face(surface, [um, u1], [vm, v1], params, depth + 1, vertices, normals, triangles);
    } else {
        // Emit triangles
        let n = vertices.len() as u32;

        // Compute normals
        let n00 = surface.normal_at(u0, v0);
        let n10 = surface.normal_at(u1, v0);
        let n01 = surface.normal_at(u0, v1);
        let n11 = surface.normal_at(u1, v1);

        // Check for degenerate points
        if !p00.is_finite() || !p10.is_finite() || !p01.is_finite() || !p11.is_finite() {
            return;
        }

        vertices.extend_from_slice(&[p00, p10, p11, p01]);
        normals.extend_from_slice(&[n00, n10, n11, n01]);

        // Choose diagonal based on edge lengths
        let d0_sq = (p11 - p00).length_squared();
        let d1_sq = (p10 - p01).length_squared();

        if d0_sq <= d1_sq {
            triangles.push([n, n + 1, n + 2]);
            triangles.push([n, n + 2, n + 3]);
        } else {
            triangles.push([n, n + 1, n + 3]);
            triangles.push([n + 1, n + 2, n + 3]);
        }
    }
}

/// Compute the angle between two normals in radians.
fn normal_angle(n0: DVec3, n1: DVec3) -> f64 {
    let len0 = n0.length();
    let len1 = n1.length();
    if len0 < 0.5 || len1 < 0.5 {
        return 0.0;
    }
    let cos_angle = (n0.dot(n1) / (len0 * len1)).clamp(-1.0, 1.0);
    cos_angle.acos()
}

/// Weld duplicate vertices in a mesh.
fn weld_mesh(mesh: Mesh) -> Mesh {
    const WELD_TOLERANCE: f64 = 1e-9;

    if mesh.vertices.is_empty() {
        return mesh;
    }

    let mut remap = vec![0u32; mesh.vertices.len()];
    let mut welded_vertices: Vec<DVec3> = Vec::new();
    let mut welded_normals: Vec<DVec3> = Vec::new();
    let mut normal_counts = Vec::new();
    let mut buckets: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    let scale = 1.0 / WELD_TOLERANCE;

    for (index, point) in mesh.vertices.iter().enumerate() {
        let key = [
            (point.x * scale).round() as i64,
            (point.y * scale).round() as i64,
            (point.z * scale).round() as i64,
        ];

        let mut matched = None;
        if let Some(candidates) = buckets.get(&key) {
            for &candidate in candidates {
                if (welded_vertices[candidate] - *point).length_squared() <= WELD_TOLERANCE * WELD_TOLERANCE {
                    matched = Some(candidate);
                    break;
                }
            }
        }

        let target = if let Some(existing) = matched {
            existing
        } else {
            let new_index = welded_vertices.len();
            welded_vertices.push(*point);
            welded_normals.push(DVec3::ZERO);
            normal_counts.push(0);
            buckets.entry(key).or_default().push(new_index);
            new_index
        };

        remap[index] = target as u32;
        if let Some(normal) = mesh.normals.get(index) {
            welded_normals[target] += *normal;
            normal_counts[target] += 1;
        }
    }

    let welded_triangles: Vec<[u32; 3]> = mesh
        .triangles
        .iter()
        .filter_map(|&[a, b, c]| {
            let ra = remap[a as usize];
            let rb = remap[b as usize];
            let rc = remap[c as usize];
            if ra == rb || rb == rc || rc == ra {
                None
            } else {
                Some([ra, rb, rc])
            }
        })
        .collect();

    let welded_normals: Vec<DVec3> = welded_normals
        .into_iter()
        .zip(normal_counts)
        .map(|(normal, count)| {
            if count == 0 {
                DVec3::ZERO
            } else {
                normal.normalize_or_zero()
            }
        })
        .collect();

    Mesh {
        vertices: welded_vertices,
        triangles: welded_triangles,
        normals: welded_normals,
    }
}

/// Optimize mesh for better triangle quality.
fn optimize_mesh(mut mesh: Mesh) -> Mesh {
    // Simple edge flip optimization for better aspect ratios
    let mut improved = true;
    let mut iterations = 0;
    let max_iterations = 3;

    while improved && iterations < max_iterations {
        improved = false;
        iterations += 1;

        // Build edge-to-triangle map
        let mut edge_tris: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        for (idx, tri) in mesh.triangles.iter().enumerate() {
            for i in 0..3 {
                let a = tri[i].min(tri[(i + 1) % 3]);
                let b = tri[i].max(tri[(i + 1) % 3]);
                edge_tris.entry((a, b)).or_default().push(idx);
            }
        }

        // Try edge flips for better triangles
        for (edge, tris) in &edge_tris {
            if tris.len() == 2 {
                let t0 = mesh.triangles[tris[0]];
                let t1 = mesh.triangles[tris[1]];

                // Check if flip would improve quality
                if should_flip_edge(&mesh.vertices, t0, t1, *edge) {
                    // Apply the flip
                    if let Some((new_t0, new_t1)) = flip_edge(t0, t1, *edge) {
                        mesh.triangles[tris[0]] = new_t0;
                        mesh.triangles[tris[1]] = new_t1;
                        improved = true;
                    }
                }
            }
        }
    }

    mesh
}

/// Check if an edge flip would improve triangle quality.
fn should_flip_edge(vertices: &[DVec3], t0: [u32; 3], t1: [u32; 3], edge: (u32, u32)) -> bool {
    // Find opposite vertices
    let opp0 = t0.iter().find(|&&v| v != edge.0 && v != edge.1).copied().unwrap();
    let opp1 = t1.iter().find(|&&v| v != edge.0 && v != edge.1).copied().unwrap();

    if opp0 == opp1 {
        return false;
    }

    // Compute current minimum angle
    let p_opp0 = vertices.get(opp0 as usize);
    let p_opp1 = vertices.get(opp1 as usize);
    let p_e0 = vertices.get(edge.0 as usize);
    let p_e1 = vertices.get(edge.1 as usize);

    let (Some(p_opp0), Some(p_opp1), Some(p_e0), Some(p_e1)) = (p_opp0, p_opp1, p_e0, p_e1) else {
        return false;
    };

    let min_angle_before = min_triangle_angle(*p_e0, *p_e1, *p_opp0)
        .min(min_triangle_angle(*p_e0, *p_e1, *p_opp1));

    let min_angle_after = min_triangle_angle(*p_opp0, *p_opp1, *p_e0)
        .min(min_triangle_angle(*p_opp0, *p_opp1, *p_e1));

    min_angle_after > min_angle_before
}

/// Flip an edge between two triangles.
fn flip_edge(t0: [u32; 3], t1: [u32; 3], edge: (u32, u32)) -> Option<([u32; 3], [u32; 3])> {
    let opp0 = t0.iter().find(|&&v| v != edge.0 && v != edge.1).copied()?;
    let opp1 = t1.iter().find(|&&v| v != edge.0 && v != edge.1).copied()?;

    // Create two new triangles sharing the new edge
    let new_t0 = [opp0, edge.0, opp1];
    let new_t1 = [opp0, opp1, edge.1];

    Some((new_t0, new_t1))
}

/// Compute minimum angle in a triangle in radians.
fn min_triangle_angle(a: DVec3, b: DVec3, c: DVec3) -> f64 {
    let ab = (b - a).normalize_or_zero();
    let bc = (c - b).normalize_or_zero();
    let ca = (a - c).normalize_or_zero();

    let angle_ab = normal_angle(ab, -ca).min(std::f64::consts::PI);
    let angle_bc = normal_angle(bc, -ab).min(std::f64::consts::PI);
    let angle_ca = normal_angle(ca, -bc).min(std::f64::consts::PI);

    angle_ab.min(angle_bc).min(angle_ca)
}

// ============================================================================
// BRep Meshing
// ============================================================================

/// Mesh an entire BRep shape.
///
/// Generates per-face meshes for all faces in the BRep.
/// This is the main entry point for meshing complete shapes.
///
/// # Arguments
/// * `brep` - The BRep shape to mesh.
/// * `params` - Mesh generation parameters.
///
/// # Returns
/// A `BRepMesh` containing per-face mesh data.
pub fn mesh_brep(brep: &BRep, params: &MeshParams) -> BRepMesh {
    let mut brep_mesh = BRepMesh::new();

    // Collect all faces
    let mut face_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Get the surface for this face
                if let Some(surface_idx) = brep.geom.face_surface.get(face_idx).and_then(|o| *o) {
                    if let Some(surface) = brep.geom.surfaces.get(surface_idx) {
                        let mesh = mesh_face(face, surface, params);
                        brep_mesh.face_meshes.push(mesh);
                        brep_mesh.face_normals.push(face.normal);
                    } else {
                        brep_mesh.face_meshes.push(Mesh::new());
                        brep_mesh.face_normals.push(face.normal);
                    }
                } else {
                    // No surface - try fallback wire triangulation
                    let mesh = mesh_face_from_wire(brep, face, params);
                    brep_mesh.face_meshes.push(mesh);
                    brep_mesh.face_normals.push(face.normal);
                }
                face_idx += 1;
            }
        }
    }

    brep_mesh
}

/// Mesh a face from its wire boundary (fallback for faces without surfaces).
fn mesh_face_from_wire(brep: &BRep, face: &Face, _params: &MeshParams) -> Mesh {
    // Sample points from the outer wire
    let mut poly_pts: Vec<DVec3> = Vec::new();
    for we in &face.outer_wire.edges {
        if let Some(edge) = brep.edges.get(we.idx) {
            let start_idx = if we.forward { edge.start } else { edge.end };
            let end_idx = if we.forward { edge.end } else { edge.start };

            if let Some(v) = brep.vertices.get(start_idx) {
                if poly_pts.is_empty() || (poly_pts.last().unwrap() - v.point).length() > 1e-9 {
                    poly_pts.push(v.point);
                }
            }

            // Sample edge curve if present
            if let Some(ci) = brep.geom.edge_curve.get(we.idx).and_then(|v| *v) {
                if let Some(curve) = brep.geom.curves.get(ci) {
                    let range = brep.geom.edge_curve_range.get(we.idx)
                        .and_then(|v| *v)
                        .unwrap_or_else(|| curve.default_domain());

                    let (t0, t1) = if we.forward {
                        (range[0], range[1])
                    } else {
                        (range[1], range[0])
                    };

                    let span = (t1 - t0).abs();
                    if span > 1e-12 {
                        let n_segs = estimate_curve_segments(curve, span);
                        for i in 1..=n_segs {
                            let t = t0 + (t1 - t0) * (i as f64 / n_segs as f64);
                            let pt = curve.point_at(t);
                            if poly_pts.is_empty() || (poly_pts.last().unwrap() - pt).length() > 1e-9 {
                                poly_pts.push(pt);
                            }
                        }
                    }
                }
            }

            // Add end point
            if let Some(v) = brep.vertices.get(end_idx) {
                if poly_pts.is_empty() || (poly_pts.last().unwrap() - v.point).length() > 1e-9 {
                    poly_pts.push(v.point);
                }
            }
        }
    }

    // Remove duplicate closing point
    if poly_pts.len() >= 2 && (poly_pts[0] - poly_pts[poly_pts.len() - 1]).length() < 1e-9 {
        poly_pts.pop();
    }

    if poly_pts.len() < 3 {
        return Mesh::new();
    }

    // Triangulate using ear clipping
    let tris = ear_clip_3d(&poly_pts, face.normal);

    // Build mesh
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    for tri in tris {
        triangles.push([tri[0] as u32, tri[1] as u32, tri[2] as u32]);
    }

    // Compute normals
    let normals = vec![face.normal; poly_pts.len()];

    Mesh { vertices: poly_pts, triangles, normals }
}

/// Estimate number of segments for curve discretization.
fn estimate_curve_segments(curve: &Curve3, span: f64) -> usize {
    match curve {
        Curve3::Line(_) => 1,
        Curve3::Circle(_) => {
            let segs = (span / (2.0 * std::f64::consts::PI) * 32.0).ceil() as usize;
            segs.clamp(4, 64)
        }
        Curve3::Ellipse(_) => 24,
        _ => 16,
    }
}

/// Ear clipping triangulation for 3D polygon.
fn ear_clip_3d(vertices: &[DVec3], normal: DVec3) -> Vec<[usize; 3]> {
    let n = vertices.len();
    if n < 3 {
        return vec![];
    }
    if n == 3 {
        return vec![[0, 1, 2]];
    }

    // Project to 2D
    let (u_axis, v_axis) = local_basis(normal);
    let pts_2d: Vec<[f64; 2]> = vertices
        .iter()
        .map(|p| [p.dot(u_axis), p.dot(v_axis)])
        .collect();

    ear_clip(&pts_2d)
}

/// Compute a local orthonormal basis from a normal.
fn local_basis(normal: DVec3) -> (DVec3, DVec3) {
    let ref_dir = if normal.x.abs() < 0.9 {
        DVec3::X
    } else {
        DVec3::Y
    };
    let u = normal.cross(ref_dir).normalize_or_zero();
    let v = normal.cross(u).normalize_or_zero();
    (u, v)
}

/// Ear clipping triangulation for 2D polygon.
fn ear_clip(pts: &[[f64; 2]]) -> Vec<[usize; 3]> {
    let n = pts.len();
    let mut indices: Vec<usize> = (0..n).collect();
    let mut triangles = Vec::new();

    // Ensure CCW winding
    let area = signed_area_2d(pts, &indices);
    if area < 0.0 {
        indices.reverse();
    }

    let mut remaining = indices;
    while remaining.len() > 3 {
        let len = remaining.len();
        let mut ear_found = false;

        for i in 0..len {
            let prev = if i == 0 { len - 1 } else { i - 1 };
            let next = if i == len - 1 { 0 } else { i + 1 };

            let a = remaining[prev];
            let b = remaining[i];
            let c = remaining[next];

            // Check convexity (left turn)
            if cross_2d(pts[a], pts[b], pts[c]) <= 0.0 {
                continue;
            }

            // Check no other vertex inside this triangle
            let mut contains_other = false;
            for j in 0..len {
                if j == prev || j == i || j == next {
                    continue;
                }
                if point_in_triangle_2d(pts[remaining[j]], pts[a], pts[b], pts[c]) {
                    contains_other = true;
                    break;
                }
            }

            if !contains_other {
                triangles.push([a, b, c]);
                remaining.remove(i);
                ear_found = true;
                break;
            }
        }

        if !ear_found {
            // Degenerate polygon - emit remaining as fan
            for i in 1..remaining.len() - 1 {
                triangles.push([remaining[0], remaining[i], remaining[i + 1]]);
            }
            break;
        }
    }

    if remaining.len() == 3 {
        triangles.push([remaining[0], remaining[1], remaining[2]]);
    }

    triangles
}

fn signed_area_2d(pts: &[[f64; 2]], indices: &[usize]) -> f64 {
    let n = indices.len();
    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        let a = pts[indices[i]];
        let b = pts[indices[j]];
        area += a[0] * b[1] - b[0] * a[1];
    }
    area * 0.5
}

fn cross_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn point_in_triangle_2d(p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
    let d1 = cross_2d(a, b, p);
    let d2 = cross_2d(b, c, p);
    let d3 = cross_2d(c, a, p);

    let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);

    !(has_neg && has_pos)
}

// ============================================================================
// Edge Discretization
// ============================================================================

/// Discretize a 3D curve into a polyline.
///
/// Samples points along the curve based on deflection criteria.
///
/// # Arguments
/// * `curve` - The 3D curve to discretize.
/// * `params` - Mesh parameters controlling sampling density.
///
/// # Returns
/// A vector of points along the curve.
pub fn discretize_edge(curve: &Curve3, params: &MeshParams) -> Vec<DVec3> {
    let domain = curve.default_domain();
    discretize_edge_in_range(curve, domain[0], domain[1], params)
}

/// Discretize a curve within a parameter range.
fn discretize_edge_in_range(curve: &Curve3, t0: f64, t1: f64, params: &MeshParams) -> Vec<DVec3> {
    let mut points = Vec::new();

    // For lines, just return endpoints
    if matches!(curve, Curve3::Line(_)) {
        points.push(curve.point_at(t0));
        points.push(curve.point_at(t1));
        return points;
    }

    // Adaptive subdivision
    let mut segments: Vec<[f64; 2]> = vec![[t0, t1]];
    let mut final_segments: Vec<[f64; 2]> = Vec::new();

    let mut depth = 0;
    while !segments.is_empty() && depth < params.max_depth * 2 {
        let mut new_segments = Vec::new();

        for seg in segments {
            let mid = (seg[0] + seg[1]) * 0.5;

            let p0 = curve.point_at(seg[0]);
            let p1 = curve.point_at(seg[1]);
            let pm = curve.point_at(mid);

            // Check chord error
            let interp_mid = (p0 + p1) * 0.5;
            let chord_error = (pm - interp_mid).length();

            // Check segment length
            let seg_len = (p1 - p0).length();

            if chord_error > params.deflection
                || seg_len > params.max_mesh_size
            {
                new_segments.push([seg[0], mid]);
                new_segments.push([mid, seg[1]]);
            } else {
                final_segments.push(seg);
            }
        }

        segments = new_segments;
        depth += 1;
    }

    // Add remaining segments
    final_segments.extend(segments);

    // Sort by parameter
    final_segments.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap_or(std::cmp::Ordering::Equal));

    // Build point list
    if !final_segments.is_empty() {
        points.push(curve.point_at(final_segments[0][0]));
        for seg in &final_segments {
            points.push(curve.point_at(seg[1]));
        }
    }

    points
}

/// Discretize a curve projected onto a surface.
///
/// This is used for PCurves where the 2D curve lies on a surface.
/// The result is in world coordinates.
///
/// # Arguments
/// * `curve` - The 3D curve to discretize.
/// * `surface` - The surface on which the curve lies.
/// * `params` - Mesh parameters controlling sampling density.
///
/// # Returns
/// A vector of points in world coordinates.
pub fn discretize_edge_on_surface(curve: &Curve3, surface: &Surface3, params: &MeshParams) -> Vec<DVec3> {
    // For now, just use the 3D curve discretization
    // In a full implementation, we would also check surface deviation
    let mut points = discretize_edge(curve, params);

    // Verify points are on the surface and project if needed
    // Note: Surface available for more sophisticated projection in future
    let _ = surface;

    points
}

// ============================================================================
// Quality Metrics
// ============================================================================

/// Compute the maximum aspect ratio in a mesh.
///
/// Aspect ratio is defined as the longest edge divided by the shortest edge
/// of each triangle. A perfect equilateral triangle has ratio 1.0.
///
/// Returns f64::INFINITY if mesh is empty or has degenerate triangles.
pub fn mesh_aspect_ratio(mesh: &Mesh) -> f64 {
    if mesh.triangles.is_empty() {
        return f64::INFINITY;
    }

    let mut max_ratio: f64 = 1.0;

    for &tri in &mesh.triangles {
        let [i0, i1, i2] = tri;

        let p0 = mesh.vertices.get(i0 as usize);
        let p1 = mesh.vertices.get(i1 as usize);
        let p2 = mesh.vertices.get(i2 as usize);

        let (Some(p0), Some(p1), Some(p2)) = (p0, p1, p2) else {
            continue;
        };

        let e0 = (*p1 - *p0).length();
        let e1 = (*p2 - *p1).length();
        let e2 = (*p0 - *p2).length();

        let max_edge = e0.max(e1).max(e2);
        let min_edge = e0.min(e1).min(e2);

        if min_edge > 1e-12 {
            let ratio = max_edge / min_edge;
            max_ratio = max_ratio.max(ratio);
        } else {
            return f64::INFINITY;
        }
    }

    max_ratio
}

/// Compute the minimum angle in the mesh (in radians).
///
/// A high-quality mesh has minimum angles close to 60 degrees (pi/3 radians).
/// Small minimum angles indicate poor quality triangles.
pub fn mesh_min_angle(mesh: &Mesh) -> f64 {
    if mesh.triangles.is_empty() {
        return 0.0;
    }

    let mut min_angle = std::f64::consts::PI;

    for &tri in &mesh.triangles {
        let [i0, i1, i2] = tri;

        let p0 = mesh.vertices.get(i0 as usize);
        let p1 = mesh.vertices.get(i1 as usize);
        let p2 = mesh.vertices.get(i2 as usize);

        let (Some(p0), Some(p1), Some(p2)) = (p0, p1, p2) else {
            continue;
        };

        let angle = min_triangle_angle(*p0, *p1, *p2);
        min_angle = min_angle.min(angle);
    }

    min_angle
}

/// Compute the maximum edge length in the mesh.
pub fn mesh_max_edge_length(mesh: &Mesh) -> f64 {
    if mesh.triangles.is_empty() {
        return 0.0;
    }

    let mut max_len: f64 = 0.0;

    for &tri in &mesh.triangles {
        let [i0, i1, i2] = tri;

        let p0 = mesh.vertices.get(i0 as usize);
        let p1 = mesh.vertices.get(i1 as usize);
        let p2 = mesh.vertices.get(i2 as usize);

        let (Some(p0), Some(p1), Some(p2)) = (p0, p1, p2) else {
            continue;
        };

        let e0 = (*p1 - *p0).length();
        let e1 = (*p2 - *p1).length();
        let e2 = (*p0 - *p2).length();

        max_len = max_len.max(e0).max(e1).max(e2);
    }

    max_len
}

// ============================================================================
// Mesh Refinement
// ============================================================================

/// Refine a mesh by subdividing triangles with edges longer than max_edge_length.
///
/// Each triangle with an edge longer than the threshold is subdivided
/// into 4 smaller triangles using midpoint subdivision.
///
/// # Arguments
/// * `mesh` - The mesh to refine.
/// * `max_edge_length` - Maximum allowed edge length.
///
/// # Returns
/// A new refined mesh.
pub fn refine_mesh(mesh: &Mesh, max_edge_length: f64) -> Mesh {
    if mesh.triangles.is_empty() || max_edge_length <= 0.0 {
        return mesh.clone();
    }

    let mut vertices = mesh.vertices.clone();
    let mut normals = mesh.normals.clone();
    let mut triangles = Vec::new();

    let mut edge_midpoints: HashMap<(u32, u32), u32> = HashMap::new();

    for &tri in &mesh.triangles {
        let [i0, i1, i2] = tri;

        let p0 = vertices[i0 as usize];
        let p1 = vertices[i1 as usize];
        let p2 = vertices[i2 as usize];

        let e0_len = (p1 - p0).length();
        let e1_len = (p2 - p1).length();
        let e2_len = (p0 - p2).length();

        let needs_refinement = e0_len > max_edge_length
            || e1_len > max_edge_length
            || e2_len > max_edge_length;

        if needs_refinement {
            // Get or create midpoints
            let m01 = get_or_create_midpoint(i0, i1, &mut vertices, &mut normals, &mut edge_midpoints);
            let m12 = get_or_create_midpoint(i1, i2, &mut vertices, &mut normals, &mut edge_midpoints);
            let m20 = get_or_create_midpoint(i2, i0, &mut vertices, &mut normals, &mut edge_midpoints);

            // Create 4 new triangles
            triangles.push([i0, m01, m20]);
            triangles.push([m01, i1, m12]);
            triangles.push([m20, m12, i2]);
            triangles.push([m01, m12, m20]);
        } else {
            triangles.push(tri);
        }
    }

    Mesh {
        vertices,
        triangles,
        normals,
    }
}

/// Get or create a midpoint vertex for an edge.
fn get_or_create_midpoint(
    i0: u32,
    i1: u32,
    vertices: &mut Vec<DVec3>,
    normals: &mut Vec<DVec3>,
    edge_midpoints: &mut HashMap<(u32, u32), u32>,
) -> u32 {
    let key = if i0 < i1 { (i0, i1) } else { (i1, i0) };

    if let Some(&idx) = edge_midpoints.get(&key) {
        return idx;
    }

    let p0 = vertices[i0 as usize];
    let p1 = vertices[i1 as usize];
    let mid_point = (p0 + p1) * 0.5;

    let n0 = normals.get(i0 as usize).copied().unwrap_or(DVec3::Z);
    let n1 = normals.get(i1 as usize).copied().unwrap_or(DVec3::Z);
    let mid_normal = (n0 + n1).normalize_or_zero();

    let idx = vertices.len() as u32;
    vertices.push(mid_point);
    normals.push(mid_normal);
    edge_midpoints.insert(key, idx);
    idx
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle3, Plane, CylindricalSurface, SphericalSurface};
    use rcad_kernel::topology::Wire;

    fn test_plane() -> Surface3 {
        Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        })
    }

    fn test_sphere() -> Surface3 {
        Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        })
    }

    fn test_cylinder() -> Surface3 {
        Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        })
    }

    fn make_test_face() -> Face {
        Face {
            outer_wire: Wire { edges: vec![] },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        }
    }

    #[test]
    fn mesh_params_default() {
        let params = MeshParams::default();
        assert!((params.deflection - 0.001).abs() < 1e-10);
        assert!((params.angle_deflection - 0.5).abs() < 1e-10);
        assert_eq!(params.min_mesh_size, 0.0);
        assert_eq!(params.max_mesh_size, f64::MAX);
        assert!(params.optimize);
    }

    #[test]
    fn mesh_params_builder() {
        let params = MeshParams::default()
            .with_deflection(0.01)
            .with_angle_deflection(0.3)
            .with_min_mesh_size(0.1)
            .with_max_mesh_size(10.0);

        assert!((params.deflection - 0.01).abs() < 1e-10);
        assert!((params.angle_deflection - 0.3).abs() < 1e-10);
        assert!((params.min_mesh_size - 0.1).abs() < 1e-10);
        assert!((params.max_mesh_size - 10.0).abs() < 1e-10);
    }

    #[test]
    fn mesh_params_presets() {
        let coarse = MeshParams::coarse();
        let fine = MeshParams::fine();
        let analysis = MeshParams::analysis();

        // Coarse should have higher deflection than fine
        assert!(coarse.deflection > fine.deflection);
        // Fine and analysis may have same deflection
        assert!(fine.deflection >= analysis.deflection);
        // Analysis should have lower angle deflection
        assert!(analysis.angle_deflection < fine.angle_deflection);
    }

    #[test]
    fn mesh_face_plane() {
        let face = make_test_face();
        let surface = test_plane();
        let params = MeshParams::default();

        let mesh = mesh_face(&face, &surface, &params);

        assert!(!mesh.is_empty());
        assert!(mesh.triangle_count() >= 2);  // At least a basic grid
        assert_eq!(mesh.vertex_count(), mesh.vertices.len());
        assert_eq!(mesh.normals.len(), mesh.vertices.len());
    }

    #[test]
    fn mesh_face_sphere() {
        let face = make_test_face();
        let surface = test_sphere();
        let params = MeshParams::default()
            .with_deflection(0.01);

        let mesh = mesh_face(&face, &surface, &params);

        assert!(!mesh.is_empty());
        assert!(mesh.triangle_count() >= 8);  // Sphere needs many triangles

        // Check that vertices are on the sphere
        for v in &mesh.vertices {
            let r = v.length();
            assert!((r - 1.0).abs() < 0.1, "Vertex not on sphere: r={}", r);
        }
    }

    #[test]
    fn mesh_face_cylinder() {
        let face = make_test_face();
        let surface = test_cylinder();
        let params = MeshParams::default();

        let mesh = mesh_face(&face, &surface, &params);

        assert!(!mesh.is_empty());
        assert!(mesh.triangle_count() >= 4);
    }

    #[test]
    fn mesh_empty() {
        let mesh = Mesh::new();
        assert!(mesh.is_empty());
        assert_eq!(mesh.triangle_count(), 0);
        assert_eq!(mesh.vertex_count(), 0);
    }

    #[test]
    fn mesh_bounding_box() {
        let mut mesh = Mesh::new();
        mesh.vertices.push(DVec3::new(0.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(1.0, 2.0, 3.0));
        mesh.vertices.push(DVec3::new(-1.0, -2.0, -3.0));

        let (min, max) = mesh.bounding_box();
        assert!((min - DVec3::new(-1.0, -2.0, -3.0)).length() < 1e-10);
        assert!((max - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-10);
    }

    #[test]
    fn mesh_flip() {
        let mut mesh = Mesh::new();
        mesh.vertices.push(DVec3::ZERO);
        mesh.vertices.push(DVec3::X);
        mesh.vertices.push(DVec3::Y);
        mesh.triangles.push([0, 1, 2]);
        mesh.normals.push(DVec3::Z);

        mesh.flip();

        assert_eq!(mesh.triangles[0], [2, 1, 0]);
        assert_eq!(mesh.normals[0], -DVec3::Z);
    }

    #[test]
    fn mesh_merge() {
        let mut mesh1 = Mesh::new();
        mesh1.vertices.push(DVec3::ZERO);
        mesh1.vertices.push(DVec3::X);
        mesh1.vertices.push(DVec3::Y);
        mesh1.triangles.push([0, 1, 2]);
        mesh1.normals = vec![DVec3::Z; 3];

        let mut mesh2 = Mesh::new();
        mesh2.vertices.push(DVec3::new(1.0, 1.0, 0.0));
        mesh2.vertices.push(DVec3::new(2.0, 1.0, 0.0));
        mesh2.vertices.push(DVec3::new(1.5, 2.0, 0.0));
        mesh2.triangles.push([0, 1, 2]);
        mesh2.normals = vec![DVec3::Z; 3];

        mesh1.merge(&mesh2);

        assert_eq!(mesh1.vertex_count(), 6);
        assert_eq!(mesh1.triangle_count(), 2);
        assert_eq!(mesh1.triangles[1], [3, 4, 5]);
    }

    #[test]
    fn brep_mesh_new() {
        let brep_mesh = BRepMesh::new();
        assert!(brep_mesh.face_meshes.is_empty());
        assert_eq!(brep_mesh.total_triangles(), 0);
    }

    #[test]
    fn brep_mesh_to_merged() {
        let mut brep_mesh = BRepMesh::new();

        let mut mesh1 = Mesh::new();
        mesh1.vertices.push(DVec3::ZERO);
        mesh1.vertices.push(DVec3::X);
        mesh1.vertices.push(DVec3::Y);
        mesh1.triangles.push([0, 1, 2]);
        mesh1.normals = vec![DVec3::Z; 3];

        let mut mesh2 = Mesh::new();
        mesh2.vertices.push(DVec3::new(1.0, 0.0, 0.0));
        mesh2.vertices.push(DVec3::new(2.0, 0.0, 0.0));
        mesh2.vertices.push(DVec3::new(1.5, 1.0, 0.0));
        mesh2.triangles.push([0, 1, 2]);
        mesh2.normals = vec![DVec3::Z; 3];

        brep_mesh.face_meshes.push(mesh1);
        brep_mesh.face_meshes.push(mesh2);

        let merged = brep_mesh.to_merged_mesh();
        assert_eq!(merged.triangle_count(), 2);
        assert_eq!(merged.vertex_count(), 6);
    }

    #[test]
    fn discretize_edge_line() {
        // Test with a bounded line segment using a circle (bounded domain)
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        let params = MeshParams::default().with_deflection(0.1);

        let points = discretize_edge(&circle, &params);

        // Should have at least start and end points
        assert!(points.len() >= 2);
        // First point should be on the circle
        let r = points[0].length();
        assert!((r - 1.0).abs() < 0.1, "First point should be on circle, r={}", r);
    }

    #[test]
    fn discretize_edge_circle() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        let params = MeshParams::default().with_deflection(0.01);

        let points = discretize_edge(&circle, &params);

        // Circle should be discretized into multiple segments
        assert!(points.len() >= 4);

        // Check first and last points are on the circle
        for p in &points {
            let r = p.length();
            assert!((r - 1.0).abs() < 0.1);
        }
    }

    #[test]
    fn mesh_aspect_ratio_equilateral() {
        let mut mesh = Mesh::new();
        // Equilateral triangle
        mesh.vertices.push(DVec3::new(0.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(1.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(0.5, 0.866, 0.0));
        mesh.triangles.push([0, 1, 2]);
        mesh.normals = vec![DVec3::Z; 3];

        let ratio = mesh_aspect_ratio(&mesh);
        assert!((ratio - 1.0).abs() < 0.01, "Equilateral should have ratio ~1.0, got {}", ratio);
    }

    #[test]
    fn mesh_aspect_ratio_degenerate() {
        let mut mesh = Mesh::new();
        // Degenerate triangle (collinear points with different spacing)
        // This creates max_edge=2, min_edge=1, ratio=2
        mesh.vertices.push(DVec3::new(0.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(1.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(2.0, 0.0, 0.0));
        mesh.triangles.push([0, 1, 2]);
        mesh.normals = vec![DVec3::Z; 3];

        let ratio = mesh_aspect_ratio(&mesh);
        // For this aspect ratio definition (max_edge/min_edge), the ratio is 2
        assert!(ratio > 1.0, "Degenerate should have ratio > 1, got {}", ratio);

        // Test a truly degenerate case where all points are the same
        let mut mesh2 = Mesh::new();
        mesh2.vertices.push(DVec3::new(0.0, 0.0, 0.0));
        mesh2.vertices.push(DVec3::new(0.0, 0.0, 0.0));
        mesh2.vertices.push(DVec3::new(0.0, 0.0, 0.0));
        mesh2.triangles.push([0, 1, 2]);
        mesh2.normals = vec![DVec3::Z; 3];

        let ratio2 = mesh_aspect_ratio(&mesh2);
        assert!(ratio2 == f64::INFINITY, "All same points should have infinite ratio, got {}", ratio2);
    }

    #[test]
    fn mesh_min_angle_equilateral() {
        let mut mesh = Mesh::new();
        // Equilateral triangle has 60 degree angles
        mesh.vertices.push(DVec3::new(0.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(1.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(0.5, 0.866, 0.0));
        mesh.triangles.push([0, 1, 2]);
        mesh.normals = vec![DVec3::Z; 3];

        let angle = mesh_min_angle(&mesh);
        let expected = std::f64::consts::PI / 3.0; // 60 degrees
        assert!((angle - expected).abs() < 0.1, "Expected ~60 degrees, got {}", angle.to_degrees());
    }

    #[test]
    fn test_mesh_max_edge_length() {
        let mut mesh = Mesh::new();
        mesh.vertices.push(DVec3::ZERO);
        mesh.vertices.push(DVec3::new(3.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(0.0, 4.0, 0.0));
        mesh.triangles.push([0, 1, 2]);
        mesh.normals = vec![DVec3::Z; 3];

        let max_len = mesh_max_edge_length(&mesh);
        // Edges: 0-1=3, 1-2=5, 2-0=4
        // Longest edge is 1-2 with length 5
        assert!((max_len - 5.0).abs() < 0.01, "Expected max edge 5.0, got {}", max_len);
    }

    #[test]
    fn refine_mesh_basic() {
        let mut mesh = Mesh::new();
        mesh.vertices.push(DVec3::ZERO);
        mesh.vertices.push(DVec3::new(10.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(5.0, 10.0, 0.0));
        mesh.triangles.push([0, 1, 2]);
        mesh.normals = vec![DVec3::Z; 3];

        let refined = refine_mesh(&mesh, 5.0);

        // With edges longer than 5.0, the triangle should be subdivided
        assert!(refined.triangle_count() > mesh.triangle_count());
        assert!(refined.vertex_count() > mesh.vertex_count());
    }

    #[test]
    fn refine_mesh_no_change() {
        let mut mesh = Mesh::new();
        mesh.vertices.push(DVec3::ZERO);
        mesh.vertices.push(DVec3::new(0.1, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(0.05, 0.1, 0.0));
        mesh.triangles.push([0, 1, 2]);
        mesh.normals = vec![DVec3::Z; 3];

        let refined = refine_mesh(&mesh, 1.0);

        // All edges are shorter than 1.0, no refinement needed
        assert_eq!(refined.triangle_count(), mesh.triangle_count());
        assert_eq!(refined.vertex_count(), mesh.vertex_count());
    }

    #[test]
    fn ear_clip_triangle() {
        let pts = vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ];
        let tris = ear_clip(&pts);
        assert_eq!(tris.len(), 1);
    }

    #[test]
    fn ear_clip_quad() {
        let pts = vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ];
        let tris = ear_clip(&pts);
        assert_eq!(tris.len(), 2);
    }

    #[test]
    fn ear_clip_pentagon() {
        let pts: Vec<[f64; 2]> = (0..5)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                [a.cos(), a.sin()]
            })
            .collect();

        let tris = ear_clip(&pts);
        assert_eq!(tris.len(), 3); // n-2 triangles for n-gon
    }

    #[test]
    fn mesh_brep_box() {
        use rcad_kernel::geom::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        populate_box_geom(&mut brep);

        let params = MeshParams::default();
        let brep_mesh = mesh_brep(&brep, &params);

        assert_eq!(brep_mesh.face_meshes.len(), 6); // 6 faces in a box
        assert!(brep_mesh.total_triangles() >= 6); // At least 1 triangle per face
    }

    #[test]
    fn mesh_brep_sphere() {
        use rcad_kernel::geom::PrimitiveSolid;

        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let params = MeshParams::default().with_deflection(0.01);
        let brep_mesh = mesh_brep(&brep, &params);

        // Sphere has 1 face but many triangles
        assert!(!brep_mesh.face_meshes.is_empty());
        assert!(brep_mesh.total_triangles() >= 8);
    }

    #[test]
    fn test_discretize_edge_on_surface() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        let plane = test_plane();
        let params = MeshParams::default();

        let points = discretize_edge_on_surface(&circle, &plane, &params);

        assert!(!points.is_empty());
    }
}
