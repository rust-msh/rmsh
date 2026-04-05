// ---------------------------------------------------------------------------
// Mesh Quality — Tet mesh quality metrics for visualization
// ---------------------------------------------------------------------------

use std::collections::HashMap;

use emstudio_domain::msh_loader::{self, MshMesh};

use crate::mesh_data::{FieldMesh, FieldVertex};
use crate::surface_extraction;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QualityMetric {
    /// Ratio of circumradius to inradius (1.0 = ideal regular tet).
    AspectRatio,
    /// Minimum dihedral angle in degrees (ideal = ~70.5°).
    MinDihedralAngle,
    /// Scaled Jacobian (1.0 = ideal, 0.0 = degenerate).
    ConditionNumber,
}

impl QualityMetric {
    pub const ALL: &[QualityMetric] = &[
        QualityMetric::AspectRatio,
        QualityMetric::MinDihedralAngle,
        QualityMetric::ConditionNumber,
    ];

    pub fn name(&self) -> &'static str {
        match self {
            Self::AspectRatio => "Aspect Ratio",
            Self::MinDihedralAngle => "Min Dihedral Angle",
            Self::ConditionNumber => "Condition Number",
        }
    }
}

/// Compute mesh quality and produce a FieldMesh where field_value = quality metric.
pub fn compute_mesh_quality(mesh: &MshMesh, metric: QualityMetric) -> FieldMesh {
    // First, extract surface for rendering
    let mut field_mesh = surface_extraction::extract_surface(mesh);

    // Compute quality per tet, then average onto surface nodes
    let mut node_quality: HashMap<u64, (f64, u32)> = HashMap::new(); // (sum, count)

    for elem in &mesh.elements {
        if elem.element_type != msh_loader::element_types::TET4 {
            continue;
        }
        if elem.node_tags.len() < 4 {
            continue;
        }

        let p: Vec<[f64; 3]> = elem
            .node_tags
            .iter()
            .take(4)
            .map(|&tag| mesh.node_position(tag).unwrap_or([0.0; 3]))
            .collect();

        let quality = match metric {
            QualityMetric::AspectRatio => tet_aspect_ratio(&p[0], &p[1], &p[2], &p[3]),
            QualityMetric::MinDihedralAngle => tet_min_dihedral_angle(&p[0], &p[1], &p[2], &p[3]),
            QualityMetric::ConditionNumber => tet_condition_number(&p[0], &p[1], &p[2], &p[3]),
        };

        for &tag in &elem.node_tags[..4] {
            let entry = node_quality.entry(tag).or_insert((0.0, 0));
            entry.0 += quality;
            entry.1 += 1;
        }
    }

    // Map quality values to vertices using position matching
    let node_map = surface_extraction::build_node_to_vertex_map(mesh, &field_mesh);
    let mut field_min = f32::MAX;
    let mut field_max = f32::MIN;

    for (&tag, &vi) in &node_map {
        if let Some(&(sum, count)) = node_quality.get(&tag) {
            let avg = (sum / count as f64) as f32;
            field_mesh.vertices[vi as usize].field_value = avg;
            field_min = field_min.min(avg);
            field_max = field_max.max(avg);
        }
    }

    if field_min > field_max {
        field_min = 0.0;
        field_max = 1.0;
    }
    field_mesh.field_range = [field_min, field_max];
    field_mesh
}

fn edge_length(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Aspect ratio: ratio of longest edge to shortest altitude.
/// Returns a value where 1.0 = ideal regular tet.
fn tet_aspect_ratio(p0: &[f64; 3], p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3]) -> f64 {
    let edges = [
        edge_length(p0, p1),
        edge_length(p0, p2),
        edge_length(p0, p3),
        edge_length(p1, p2),
        edge_length(p1, p3),
        edge_length(p2, p3),
    ];

    let max_edge = edges.iter().cloned().fold(0.0f64, f64::max);
    let min_edge = edges.iter().cloned().fold(f64::MAX, f64::min);

    if min_edge < 1e-15 {
        return 0.0;
    }

    // Normalize: ideal regular tet has ratio 1.0
    min_edge / max_edge
}

/// Minimum dihedral angle of a tetrahedron (in degrees).
fn tet_min_dihedral_angle(p0: &[f64; 3], p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3]) -> f64 {
    // 4 faces of the tet, each defined by 3 vertices
    let faces: [(&[f64; 3], &[f64; 3], &[f64; 3]); 4] = [
        (p0, p2, p1),
        (p0, p1, p3),
        (p1, p2, p3),
        (p0, p3, p2),
    ];

    // Compute outward face normals
    let normals: Vec<[f64; 3]> = faces
        .iter()
        .map(|&(a, b, c)| {
            let e1 = sub(b, a);
            let e2 = sub(c, a);
            let n = cross(e1, e2);
            let len = vec_len(n);
            if len > 1e-15 {
                [n[0] / len, n[1] / len, n[2] / len]
            } else {
                [0.0, 0.0, 1.0]
            }
        })
        .collect();

    // 6 dihedral angles: one for each edge (shared by 2 faces)
    let face_pairs: [(usize, usize); 6] = [
        (0, 1), // edge p0-p1
        (0, 2), // edge p1-p2
        (0, 3), // edge p0-p2
        (1, 2), // edge p1-p3
        (1, 3), // edge p0-p3
        (2, 3), // edge p2-p3
    ];

    let mut min_angle = 180.0f64;
    for &(fi, fj) in &face_pairs {
        let cos_angle = dot(normals[fi], normals[fj]);
        // Dihedral angle = π - angle between outward normals
        let angle_deg = (std::f64::consts::PI - cos_angle.clamp(-1.0, 1.0).acos()).to_degrees();
        min_angle = min_angle.min(angle_deg);
    }

    min_angle
}

/// Condition number (scaled Jacobian) of a tetrahedron.
/// Returns 1.0 for ideal regular tet, 0.0 for degenerate.
fn tet_condition_number(p0: &[f64; 3], p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3]) -> f64 {
    let e0 = sub(p1, p0);
    let e1 = sub(p2, p0);
    let e2 = sub(p3, p0);

    // Volume = |det(e0, e1, e2)| / 6
    let det = dot(e0, cross(e1, e2));
    let volume = det.abs() / 6.0;

    if volume < 1e-20 {
        return 0.0;
    }

    // All 6 edge lengths
    let edges = [
        edge_length(p0, p1),
        edge_length(p0, p2),
        edge_length(p0, p3),
        edge_length(p1, p2),
        edge_length(p1, p3),
        edge_length(p2, p3),
    ];
    let sum_sq: f64 = edges.iter().map(|e| e * e).sum();

    // Normalized: for regular tet, 6*sqrt(2)*V = edge^3 * sqrt(2)/sqrt(3)
    // Use the quality metric: Q = 12 * (3*V)^(2/3) / sum_edge_sq
    // For a regular tet with edge a: V = a^3/(6*sqrt(2)), sum_sq = 6*a^2
    // Q = 12*(3*a^3/(6*sqrt(2)))^(2/3) / (6*a^2) = 12*(a^3/(2*sqrt(2)))^(2/3)/(6*a^2)
    // = 12 * a^2 / (2*sqrt(2))^(2/3) / (6*a^2) = 2 / (2*sqrt(2))^(2/3) ≈ 1.0
    let quality = 12.0 * (3.0 * volume).powf(2.0 / 3.0) / sum_sq;
    quality.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regular_tet_quality() {
        // Regular tetrahedron
        let p0 = [1.0, 1.0, 1.0];
        let p1 = [1.0, -1.0, -1.0];
        let p2 = [-1.0, 1.0, -1.0];
        let p3 = [-1.0, -1.0, 1.0];

        let ar = tet_aspect_ratio(&p0, &p1, &p2, &p3);
        assert!((ar - 1.0).abs() < 0.01, "aspect ratio should be ~1.0 for regular tet, got {}", ar);

        let angle = tet_min_dihedral_angle(&p0, &p1, &p2, &p3);
        assert!((angle - 70.5).abs() < 1.0, "min dihedral angle should be ~70.5° for regular tet, got {}", angle);

        let cn = tet_condition_number(&p0, &p1, &p2, &p3);
        assert!(cn > 0.9, "condition number should be ~1.0 for regular tet, got {}", cn);
    }
}
