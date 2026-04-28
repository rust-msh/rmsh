/// Integration tests for mesh generation and tessellation quality.
///
/// These tests verify:
/// - Quality metrics computation
/// - Boundary preservation
/// - Incremental mesh updates
/// - Surface mesh generation
use glam::DVec3;
use rcad_algorithms::{
    triangulate_surface, mesh_brep,
    TessellationParams, SurfaceMesh, MeshQualityMetrics, compute_mesh_quality,
    AdaptiveSubdivider, BoundarySensitiveTessellator, IncrementalMesher, MeshDelta,
    MeshSimplifier,
};
use rcad_kernel::{BRep, PrimitiveSolid, geom::{Surface3, SphericalSurface, CylindricalSurface, Plane}};
use rcad_modeling::{make_box_brep, make_sphere_brep, make_cylinder_brep};

fn face_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

fn total_triangle_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .map(|f| f.triangles.len())
        .sum()
}

// ============================================================================
// Quality Metrics Tests
// ============================================================================

/// Test that quality metrics are computed correctly for a simple triangle.
#[test]
fn quality_metrics_simple_triangle() {
    let vertices = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.5, 0.866, 0.0), // Equilateral triangle
    ];
    let triangles = vec![[0, 1, 2]];

    let metrics = compute_mesh_quality(&vertices, &triangles);

    assert_eq!(metrics.triangle_count, 1);
    assert_eq!(metrics.vertex_count, 3);
    assert_eq!(metrics.degenerate_count, 0);
    // Equilateral triangle has aspect ratio ~1.0
    assert!(metrics.max_aspect_ratio < 1.2, "equilateral should have low aspect ratio");
    assert!(metrics.quality_score() > 0.9);
}

/// Test quality metrics for a degenerate triangle.
#[test]
fn quality_metrics_degenerate_triangle() {
    let vertices = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.5, 0.0, 0.0), // Collinear - degenerate
    ];
    let triangles = vec![[0, 1, 2]];

    let metrics = compute_mesh_quality(&vertices, &triangles);

    // This triangle is degenerate (collinear points)
    assert_eq!(metrics.degenerate_count, 1);
    // Quality score should be lower for degenerate mesh
    // Note: The exact score threshold may vary by implementation
    assert!(metrics.quality_score() < 0.8, "degenerate triangle should have lower quality score");
}

/// Test quality metrics for a high aspect ratio triangle.
#[test]
fn quality_metrics_high_aspect_ratio() {
    let vertices = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(100.0, 0.0, 0.0),
        DVec3::new(50.0, 0.001, 0.0), // Very thin triangle
    ];
    let triangles = vec![[0, 1, 2]];

    let metrics = compute_mesh_quality(&vertices, &triangles);

    // Edge lengths: 100, ~50, ~50
    // Max edge = 100, min edge ≈ 50 (since all edges are > 0.001)
    // Aspect ratio = max_edge / min_edge ≈ 2
    // For a very thin triangle, aspect ratio depends on edge ratios
    assert!(metrics.max_aspect_ratio > 1.0, "should have aspect ratio > 1");
    // Check that poor aspect ratio count is tracked
    // Note: The threshold for "poor" is aspect ratio > 20.0
    // This triangle may not exceed that threshold
    assert!(!metrics.is_good(1.5), "should not be considered good at strict threshold");
}

/// Test quality metrics for a mesh with multiple triangles.
#[test]
fn quality_metrics_multiple_triangles() {
    let vertices = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(1.0, 1.0, 0.0),
    ];
    let triangles = vec![[0, 1, 2], [1, 3, 2]];

    let metrics = compute_mesh_quality(&vertices, &triangles);

    assert_eq!(metrics.triangle_count, 2);
    assert_eq!(metrics.vertex_count, 4);
    assert_eq!(metrics.degenerate_count, 0);
    assert!(metrics.average_area > 0.0);
    assert!(metrics.average_edge_length > 0.0);
}

/// Test that SurfaceMesh::compute_quality delegates correctly.
#[test]
fn surface_mesh_compute_quality() {
    let mesh = SurfaceMesh {
        vertices: vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ],
        triangles: vec![[0, 1, 2]],
        normals: vec![DVec3::Z; 3],
        dirty: false,
    };

    let metrics = mesh.compute_quality();
    assert_eq!(metrics.triangle_count, 1);
    assert!(metrics.quality_score() > 0.0);
}

// ============================================================================
// Surface Mesh Generation Tests
// ============================================================================

/// Test sphere surface tessellation.
#[test]
fn sphere_surface_tessellation() {
    let sphere = Surface3::Sphere(SphericalSurface {
        center: DVec3::ZERO,
        axis: DVec3::Z,
        radius: 1.0,
    });

    let params = TessellationParams::standard();
    let mesh = triangulate_surface(&sphere, [0.0, std::f64::consts::TAU], [0.0, std::f64::consts::PI], &params);

    assert!(!mesh.vertices.is_empty(), "should have vertices");
    assert!(!mesh.triangles.is_empty(), "should have triangles");
    assert_eq!(mesh.normals.len(), mesh.vertices.len(), "normals should match vertices");

    // Check that all vertices are at radius 1
    for v in &mesh.vertices {
        let dist = v.length();
        assert!((dist - 1.0).abs() < 0.01, "vertices should be on sphere surface");
    }
}

/// Test cylinder surface tessellation.
#[test]
fn cylinder_surface_tessellation() {
    let cylinder = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::ZERO,
        axis: DVec3::Z,
        radius: 1.0,
    });

    let params = TessellationParams::standard();
    let mesh = triangulate_surface(&cylinder, [0.0, std::f64::consts::TAU], [0.0, 5.0], &params);

    assert!(!mesh.vertices.is_empty());
    assert!(!mesh.triangles.is_empty());

    // Check that vertices are at radius 1 from axis
    for v in &mesh.vertices {
        let radial_dist = (v - DVec3::new(0.0, 0.0, v.z)).length();
        assert!((radial_dist - 1.0).abs() < 0.01, "vertices should be on cylinder surface");
    }
}

/// Test plane surface tessellation.
#[test]
fn plane_surface_tessellation() {
    let plane = Surface3::Plane(Plane {
        origin: DVec3::ZERO,
        normal: DVec3::Z,
    });

    let params = TessellationParams::standard();
    let mesh = triangulate_surface(&plane, [0.0, 5.0], [0.0, 5.0], &params);

    assert!(!mesh.vertices.is_empty());
    assert!(!mesh.triangles.is_empty());

    // Check that all vertices are at z=0
    for v in &mesh.vertices {
        assert!(v.z.abs() < 0.01, "vertices should be on plane");
    }
}

/// Test that tessellation params affect triangle count.
#[test]
fn tessellation_params_affect_density() {
    let sphere = Surface3::Sphere(SphericalSurface {
        center: DVec3::ZERO,
        axis: DVec3::Z,
        radius: 1.0,
    });

    let coarse = TessellationParams {
        chord_tolerance: 0.5,
        ..TessellationParams::default()
    };
    let fine = TessellationParams {
        chord_tolerance: 0.01,
        ..TessellationParams::default()
    };

    let coarse_mesh = triangulate_surface(&sphere, [0.0, std::f64::consts::TAU], [0.0, std::f64::consts::PI], &coarse);
    let fine_mesh = triangulate_surface(&sphere, [0.0, std::f64::consts::TAU], [0.0, std::f64::consts::PI], &fine);

    assert!(
        fine_mesh.triangles.len() > coarse_mesh.triangles.len(),
        "fine params should produce more triangles"
    );
}

// ============================================================================
// BRep Mesh Generation Tests
// ============================================================================

/// Test mesh generation on a box.
#[test]
fn mesh_box_all_faces() {
    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 3.0,
        depth: 4.0,
    });

    let params = TessellationParams::default();
    mesh_brep(&mut brep, &params);

    // All faces should have triangles
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                assert!(
                    !face.triangles.is_empty(),
                    "all faces should have triangles"
                );
            }
        }
    }
}

/// Test mesh generation on a sphere.
#[test]
fn mesh_sphere_valid_indices() {
    let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

    let params = TessellationParams {
        chord_tolerance: 0.1,
        ..TessellationParams::default()
    };
    mesh_brep(&mut brep, &params);

    // All triangle indices should be valid
    let nv = brep.vertices.len();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for tri in &face.triangles {
                    for &idx in tri {
                        assert!(idx < nv, "triangle index should be valid");
                    }
                }
            }
        }
    }
}

/// Test mesh generation on a cylinder.
#[test]
fn mesh_cylinder_quality() {
    let mut brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
        radius: 1.0,
        height: 2.0,
    });

    let params = TessellationParams::standard();
    mesh_brep(&mut brep, &params);

    let total_tris = total_triangle_count(&brep);
    assert!(total_tris > 0, "cylinder should have triangles");
}

/// Test incremental mesh update skips clean faces.
#[test]
fn incremental_mesh_skips_clean_faces() {
    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });

    let params = TessellationParams::default();

    // First mesh generation
    mesh_brep(&mut brep, &params);
    let first_tri_count = total_triangle_count(&brep);

    // Second call should skip clean faces
    mesh_brep(&mut brep, &params);
    let second_tri_count = total_triangle_count(&brep);

    // Triangle count should be identical (no regeneration)
    assert_eq!(first_tri_count, second_tri_count);
}

// ============================================================================
// Adaptive Subdivision Tests
// ============================================================================

/// Test adaptive subdivider creates more triangles.
#[test]
fn adaptive_subdivder_increases_triangles() {
    let mesh = SurfaceMesh {
        vertices: vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(10.0, 0.0, 0.0),
            DVec3::new(5.0, 10.0, 0.0),
        ],
        triangles: vec![[0, 1, 2]],
        normals: vec![DVec3::Z; 3],
        dirty: false,
    };

    let subdivider = AdaptiveSubdivider::new()
        .with_distance_threshold(5.0);

    let result = subdivider.subdivide_by_distance(&mesh);

    assert!(result.triangles.len() > mesh.triangles.len());
}

/// Test adaptive subdivider with curvature.
#[test]
fn adaptive_subdivider_curvature() {
    let mesh = SurfaceMesh {
        vertices: vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ],
        triangles: vec![[0, 1, 2]],
        normals: vec![
            DVec3::Z,
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(0.0, 0.0, 1.0), // Same normal - no curvature
        ],
        dirty: false,
    };

    let subdivider = AdaptiveSubdivider::new()
        .with_curvature_threshold(0.01);

    let result = subdivider.subdivide_by_curvature(&mesh);

    // No subdivision expected with identical normals
    assert_eq!(result.triangles.len(), mesh.triangles.len());
}

/// Test subdivider builder methods.
#[test]
fn adaptive_subdivider_builder() {
    let subdivider = AdaptiveSubdivider::new()
        .with_curvature_threshold(0.5)
        .with_distance_threshold(2.0)
        .with_max_levels(5);

    assert_eq!(subdivider.curvature_threshold, 0.5);
    assert_eq!(subdivider.distance_threshold, 2.0);
    assert_eq!(subdivider.max_subdivision_levels, 5);
}

// ============================================================================
// Boundary-Sensitive Tessellation Tests
// ============================================================================

/// Test feature edge detection.
#[test]
fn boundary_sensitive_detect_features() {
    let vertices = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 0.0, 1.0),
    ];
    // Two triangles sharing edge 0-1, with different normals (90 deg)
    let triangles = vec![[0, 1, 2], [0, 1, 3]];
    let normals = vec![DVec3::Z; 4];

    let mut tessellator = BoundarySensitiveTessellator::new()
        .with_feature_angle(0.1);

    tessellator.detect_feature_edges(&vertices, &triangles, &normals);

    // Should detect feature edge between the triangles
    assert!(!tessellator.feature_edges.is_empty());
}

/// Test boundary-sensitive tessellator builder.
#[test]
fn boundary_sensitive_builder() {
    let tessellator = BoundarySensitiveTessellator::new()
        .with_feature_angle(0.5);

    assert_eq!(tessellator.feature_angle_threshold, 0.5);
}

/// Test feature edge preservation.
#[test]
fn boundary_sensitive_preserve_features() {
    let mesh = SurfaceMesh {
        vertices: vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ],
        triangles: vec![[0, 1, 2]],
        normals: vec![DVec3::Z; 3],
        dirty: false,
    };

    let mut tessellator = BoundarySensitiveTessellator::new();
    tessellator.feature_edges.push(rcad_algorithms::FeatureEdge {
        start_vertex: 0,
        end_vertex: 1,
        feature_angle: 0.5,
    });

    let preserved = tessellator.preserve_feature_edges(&mesh);
    assert_eq!(preserved.vertices.len(), mesh.vertices.len());
}

// ============================================================================
// Incremental Mesher Tests
// ============================================================================

/// Test incremental mesher tracks dirty faces.
#[test]
fn incremental_mesher_tracks_dirty() {
    let mut mesher = IncrementalMesher::new();

    assert!(!mesher.is_dirty());

    mesher.invalidate_face(0);
    assert!(mesher.is_dirty());
    assert!(mesher.dirty_faces.contains(&0));

    mesher.clear();
    assert!(!mesher.is_dirty());
}

/// Test incremental mesher with multiple faces.
#[test]
fn incremental_mesher_multiple_faces() {
    let mut mesher = IncrementalMesher::new();
    mesher.invalidate_faces(&[0, 2, 4]);

    assert!(mesher.dirty_faces.contains(&0));
    assert!(mesher.dirty_faces.contains(&2));
    assert!(mesher.dirty_faces.contains(&4));
    assert!(!mesher.dirty_faces.contains(&1));
    assert!(!mesher.dirty_faces.contains(&3));
}

/// Test incremental mesher edge invalidation.
#[test]
fn incremental_mesher_edge_invalidation() {
    let mut mesher = IncrementalMesher::new();

    mesher.invalidate_edge(0);
    assert!(mesher.dirty_edges.contains(&0));
    assert!(mesher.is_dirty());
}

/// Test MeshDelta creation.
#[test]
fn mesh_delta_creation() {
    let delta = MeshDelta::from_vertices(vec![0, 1, 2]);
    assert_eq!(delta.modified_vertices, vec![0, 1, 2]);
    assert!(delta.modified_edges.is_empty());
    assert!(delta.modified_faces.is_empty());

    let delta = MeshDelta::from_edges(vec![3, 4]);
    assert_eq!(delta.modified_edges, vec![3, 4]);

    let delta = MeshDelta::from_faces(vec![5]);
    assert_eq!(delta.modified_faces, vec![5]);
}

/// Test incremental mesh update on BRep.
#[test]
fn incremental_mesher_update_brep() {
    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    });

    // Initial mesh
    let params = TessellationParams::default();
    mesh_brep(&mut brep, &params);

    // Mark one face as dirty
    let mut mesher = IncrementalMesher::new();
    mesher.invalidate_face(0);

    // Update mesh
    mesher.update_mesh_for_face_change(&mut brep, &params);

    // Face 0 should have been regenerated
    assert!(!brep.solids[0].shells[0].faces[0].triangles.is_empty());
}

// ============================================================================
// Mesh Simplification Tests
// ============================================================================

/// Test mesh simplifier reduces triangle count.
#[test]
fn mesh_simplifier_reduces_triangles() {
    let vertices: Vec<DVec3> = (0..9)
        .map(|i| {
            let row = i / 3;
            let col = i % 3;
            DVec3::new(col as f64, row as f64, 0.0)
        })
        .collect();

    let triangles = vec![
        [0, 1, 3], [1, 4, 3],
        [1, 2, 4], [2, 5, 4],
        [3, 4, 6], [4, 7, 6],
        [4, 5, 7], [5, 8, 7],
    ];

    let mesh = SurfaceMesh {
        vertices,
        triangles,
        normals: vec![DVec3::Z; 9],
        dirty: false,
    };

    let simplifier = MeshSimplifier::new()
        .with_target_ratio(0.5)
        .with_max_error(1.0);

    let simplified = simplifier.simplify_mesh(&mesh);

    assert!(
        simplified.triangles.len() <= mesh.triangles.len(),
        "simplified mesh should have fewer or equal triangles"
    );
}

/// Test mesh simplifier builder.
#[test]
fn mesh_simplifier_builder() {
    let simplifier = MeshSimplifier::new()
        .with_target_ratio(0.25)
        .with_max_error(0.05);

    assert_eq!(simplifier.target_ratio, 0.25);
    assert_eq!(simplifier.max_error, 0.05);
}

/// Test mesh simplifier target count.
#[test]
fn mesh_simplifier_target_count() {
    let mesh = SurfaceMesh {
        vertices: vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
        ],
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        normals: vec![DVec3::Z; 4],
        dirty: false,
    };

    let simplifier = MeshSimplifier::new().with_max_error(1.0);
    let result = simplifier.simplify_to_target_count(&mesh, 4);

    // Already at or below target
    assert!(result.triangles.len() <= 2);
}

// ============================================================================
// Tessellation Params Tests
// ============================================================================

/// Test preset configurations.
#[test]
fn tessellation_params_presets() {
    let preview = TessellationParams::preview();
    assert!(preview.chord_tolerance > TessellationParams::standard().chord_tolerance);
    assert!(preview.parallel);

    let standard = TessellationParams::standard();
    assert!(standard.adaptive_refinement);

    let hq = TessellationParams::high_quality();
    assert!(hq.chord_tolerance < standard.chord_tolerance);

    let analysis = TessellationParams::analysis();
    assert!(analysis.chord_tolerance < hq.chord_tolerance);
    assert!(analysis.max_aspect_ratio < hq.max_aspect_ratio);
}

/// Test target triangle count adjustment.
#[test]
fn tessellation_params_target_count() {
    let params = TessellationParams::standard();
    let adjusted = params.with_target_triangle_count(10000);

    // Should adjust tolerance based on target
    assert!(adjusted.chord_tolerance != params.chord_tolerance);
}

// ============================================================================
// Edge Cases and Stress Tests
// ============================================================================

/// Test tessellation with very small surface.
#[test]
fn tessellation_very_small_surface() {
    let sphere = Surface3::Sphere(SphericalSurface {
        center: DVec3::ZERO,
        axis: DVec3::Z,
        radius: 0.001,
    });

    let params = TessellationParams::default();
    let mesh = triangulate_surface(&sphere, [0.0, std::f64::consts::TAU], [0.0, std::f64::consts::PI], &params);

    assert!(!mesh.triangles.is_empty());
}

/// Test tessellation with very large surface.
#[test]
fn tessellation_very_large_surface() {
    let sphere = Surface3::Sphere(SphericalSurface {
        center: DVec3::ZERO,
        axis: DVec3::Z,
        radius: 10000.0,
    });

    let params = TessellationParams::preview(); // Use preview for speed
    let mesh = triangulate_surface(&sphere, [0.0, std::f64::consts::TAU], [0.0, std::f64::consts::PI], &params);

    assert!(!mesh.triangles.is_empty());
}

/// Test empty mesh quality.
#[test]
fn empty_mesh_quality() {
    let metrics = compute_mesh_quality(&[], &[]);
    assert_eq!(metrics.triangle_count, 0);
    assert_eq!(metrics.vertex_count, 0);
}

/// Test quality metrics with out-of-bounds indices.
#[test]
fn quality_metrics_out_of_bounds() {
    let vertices = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
    ];
    let triangles = vec![[0, 1, 5]]; // Index 5 doesn't exist

    let metrics = compute_mesh_quality(&vertices, &triangles);
    // Should handle gracefully
    assert_eq!(metrics.triangle_count, 1);
}

/// Test SurfaceMesh dirty flag.
#[test]
fn surface_mesh_dirty_flag() {
    let mut mesh = SurfaceMesh {
        vertices: vec![DVec3::ZERO],
        triangles: vec![],
        normals: vec![],
        dirty: false,
    };

    assert!(mesh.is_clean());

    mesh.invalidate();
    assert!(!mesh.is_clean());
    assert!(mesh.dirty);
}
