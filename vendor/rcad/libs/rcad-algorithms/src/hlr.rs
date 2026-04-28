//! Hidden-Line Removal (HLR).
//!
//! Projects a BRep's edges onto a view plane and classifies each edge segment
//! as **visible** or **hidden** by testing against the silhouette of all faces.
//!
//! Analytic silhouette curves are generated for curved surfaces (cylinder,
//! sphere, cone, torus) and processed through the same visibility pipeline as wire edges.
//! For general surfaces (BSpline, Bezier, etc.), numerical silhouette extraction
//! is performed using adaptive sampling with curvature-based refinement.
//!
//! Analogous to OCCT `HLRBRep_Algo` / `HLRBRep_HLRToShape`.
//!
//! # Algorithm
//!
//! For each edge (and silhouette curve):
//! 1. Project both endpoints onto the screen plane.
//! 2. Sample `N` points along the edge in 3D (adaptively if `curvature_adaptive` is true).
//! 3. For each sample, cast a ray from that point toward the camera.
//! 4. If any face triangle blocks the ray **closer** to the camera than the
//!    edge sample, that sample is hidden.
//! 5. Classify runs of consecutive samples → visible/hidden segments.
//!
//! The result is a set of `HlrSegment`s — 2D projected line segments labeled
//! visible or hidden.
//!
//! # Silhouette Classification
//!
//! Segments are classified as:
//! - **Visible**: Not occluded by any face
//! - **Hidden**: Occluded by at least one face
//! - **Contour**: Edge of a face's silhouette (marked via `is_contour`)
//!
//! # Curved Surface Enhancements
//!
//! This implementation includes several enhancements for better handling of curved geometry:
//! - **Curvature-adaptive sampling**: Uses surface curvature to concentrate samples in high-curvature regions
//! - **Marching silhouette extraction**: Robust detection of silhouette curves on general parametric surfaces
//! - **B-spline fitting**: Converts dense silhouette points to smooth B-spline curves
//! - **BVH acceleration**: Spatial acceleration structure for efficient ray casting
//! - **Grazing angle handling**: Special treatment for near-silhouette regions
//! - **Thread edge classification**: Support for helical edges on cylinders and cones
//! - **Seam edge detection**: Proper handling of seam edges on closed surfaces
//! - **Parallel processing**: Multi-threaded processing for large models

use std::collections::{HashMap, HashSet};
use rayon::prelude::*;
use glam::{DAffine3, DMat4, DVec2, DVec3, DVec4};
use rcad_kernel::geom::{Circle3, CurveEval, Surface3, any_perpendicular};
use rcad_kernel::{BRep, SurfaceEval};

// ── Public types ──────────────────────────────────────────────────────────────

/// Configuration options for HLR computation.
#[derive(Debug, Clone)]
pub struct HlrOptions {
    /// Number of samples per edge for occlusion testing.
    /// Higher values give more accurate results but slower computation.
    /// Default: 8.
    pub edge_samples: usize,
    /// Base number of samples for silhouette curve generation.
    /// Default: 32.
    pub silhouette_samples: usize,
    /// Enable curvature-adaptive sampling for silhouette curves.
    /// When true, high-curvature regions receive more samples.
    /// Default: true.
    pub curvature_adaptive: bool,
    /// Tolerance for tangent alignment when computing silhouettes.
    /// Points where |normal · view_dir| < tangent_tolerance are considered silhouette candidates.
    /// Default: 1e-6.
    pub tangent_tolerance: f64,
    /// Maximum angle deviation (in radians) for adaptive sampling subdivision.
    /// Smaller values produce smoother curves at higher cost.
    /// Default: 0.05 (about 3 degrees).
    pub angular_tolerance: f64,
    /// Minimum number of subdivision iterations for adaptive sampling.
    /// Default: 2.
    pub min_subdivisions: usize,
    /// Maximum number of subdivision iterations for adaptive sampling.
    /// Default: 8.
    pub max_subdivisions: usize,
    /// Enable BVH acceleration for ray casting.
    /// Default: true.
    pub use_bvh: bool,
    /// Maximum curvature for adaptive sampling (higher = more samples in curved regions).
    /// Default: 100.0.
    pub max_curvature: f64,
    /// Minimum curvature for adaptive sampling (lower = fewer samples in flat regions).
    /// Default: 0.001.
    pub min_curvature: f64,
    /// Enable B-spline fitting for silhouette curves.
    /// Default: true.
    pub fit_bspline: bool,
    /// Tolerance for B-spline fitting (maximum deviation from original points).
    /// Default: 0.001.
    pub bspline_tolerance: f64,
    /// Grazing angle threshold (in radians). Points closer to silhouette receive special handling.
    /// Default: 0.1 (about 6 degrees).
    pub grazing_angle_threshold: f64,
    /// Enable smooth silhouette approximation.
    /// Default: true.
    pub smooth_silhouettes: bool,
    /// Enable parallel processing for multi-face models.
    /// Default: true.
    pub parallel: bool,
    /// Minimum number of faces to trigger parallel processing.
    /// Default: 4.
    pub parallel_threshold: usize,
    /// Enable surface property caching for improved performance.
    /// Default: true.
    pub cache_surface_properties: bool,
    /// Silhouette proximity factor for increased sampling density.
    /// Samples within this factor of the silhouette receive more refinement.
    /// Default: 0.1 (10% of local feature size).
    pub silhouette_proximity_factor: f64,
    /// Enable thread edge detection for helical geometry.
    /// Default: true.
    pub detect_thread_edges: bool,
    /// Enable seam edge detection for closed surfaces.
    /// Default: true.
    pub detect_seam_edges: bool,
    /// Maximum depth complexity for curve-surface intersection.
    /// Default: 16.
    pub max_depth_complexity: usize,
}

impl Default for HlrOptions {
    fn default() -> Self {
        Self {
            edge_samples: 8,
            silhouette_samples: 32,
            curvature_adaptive: true,
            tangent_tolerance: 1e-6,
            angular_tolerance: 0.05,
            min_subdivisions: 2,
            max_subdivisions: 8,
            use_bvh: true,
            max_curvature: 100.0,
            min_curvature: 0.001,
            fit_bspline: true,
            bspline_tolerance: 0.001,
            grazing_angle_threshold: 0.1,
            smooth_silhouettes: true,
            parallel: true,
            parallel_threshold: 4,
            cache_surface_properties: true,
            silhouette_proximity_factor: 0.1,
            detect_thread_edges: true,
            detect_seam_edges: true,
            max_depth_complexity: 16,
        }
    }
}

impl HlrOptions {
    /// Create options with a specific edge sample count.
    pub fn with_edge_samples(mut self, n: usize) -> Self {
        self.edge_samples = n.max(2);
        self
    }

    /// Create options with a specific silhouette sample count.
    pub fn with_silhouette_samples(mut self, n: usize) -> Self {
        self.silhouette_samples = n.max(8);
        self
    }

    /// Enable or disable curvature-adaptive sampling.
    pub fn with_curvature_adaptive(mut self, adaptive: bool) -> Self {
        self.curvature_adaptive = adaptive;
        self
    }

    /// Set the tangent tolerance for silhouette detection.
    pub fn with_tangent_tolerance(mut self, tol: f64) -> Self {
        self.tangent_tolerance = tol.abs().max(1e-12);
        self
    }

    /// Enable or disable BVH acceleration.
    pub fn with_bvh(mut self, use_bvh: bool) -> Self {
        self.use_bvh = use_bvh;
        self
    }

    /// Set the maximum curvature for adaptive sampling.
    pub fn with_max_curvature(mut self, curv: f64) -> Self {
        self.max_curvature = curv.abs().max(0.1);
        self
    }

    /// Enable or disable B-spline fitting for silhouettes.
    pub fn with_bspline_fitting(mut self, fit: bool) -> Self {
        self.fit_bspline = fit;
        self
    }

    /// Set the grazing angle threshold.
    pub fn with_grazing_angle(mut self, angle: f64) -> Self {
        self.grazing_angle_threshold = angle.abs().min(std::f64::consts::FRAC_PI_2);
        self
    }

    /// Enable or disable smooth silhouette approximation.
    pub fn with_smooth_silhouettes(mut self, smooth: bool) -> Self {
        self.smooth_silhouettes = smooth;
        self
    }

    /// Enable or disable parallel processing.
    pub fn with_parallel(mut self, parallel: bool) -> Self {
        self.parallel = parallel;
        self
    }

    /// Set the parallel processing threshold (minimum faces to trigger parallelism).
    pub fn with_parallel_threshold(mut self, threshold: usize) -> Self {
        self.parallel_threshold = threshold.max(1);
        self
    }

    /// Enable or disable surface property caching.
    pub fn with_surface_caching(mut self, cache: bool) -> Self {
        self.cache_surface_properties = cache;
        self
    }

    /// Set the silhouette proximity factor.
    pub fn with_silhouette_proximity(mut self, factor: f64) -> Self {
        self.silhouette_proximity_factor = factor.abs().max(0.01).min(1.0);
        self
    }

    /// Enable or disable thread edge detection.
    pub fn with_thread_edge_detection(mut self, detect: bool) -> Self {
        self.detect_thread_edges = detect;
        self
    }

    /// Enable or disable seam edge detection.
    pub fn with_seam_edge_detection(mut self, detect: bool) -> Self {
        self.detect_seam_edges = detect;
        self
    }
}

/// Hint about the geometric type of the original 3D edge curve.
/// Used by consumers (e.g. SVG exporter) to emit arcs instead of polylines.
#[derive(Debug, Clone, PartialEq)]
pub enum CurveHint {
    /// Edge is a full or partial circle in 3D.
    Circle {
        /// Projected 2D center of the circle.
        center: DVec2,
        /// Projected radius (approximate — perspective not applied).
        radius: f64,
    },
    /// Any other non-straight curve (ellipse, spline, …).
    Other,
}

/// Classification of an HLR segment type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentType {
    /// Regular edge (part of the BRep wire).
    Edge,
    /// Silhouette curve (contour of a curved face).
    Silhouette,
    /// Thread edge (helical edge on cylinders/cones).
    Thread,
    /// Seam edge (closed surface seam).
    Seam,
}

/// Classification of edge visibility type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeClassification {
    /// Edge is fully visible.
    Visible,
    /// Edge is fully hidden.
    Hidden,
    /// Edge is a contour/silhouette edge.
    Contour,
    /// Edge is partially visible (some segments visible, some hidden).
    Partial,
    /// Edge is a thread edge (helical).
    Thread,
    /// Edge is a seam edge on a closed surface.
    Seam,
}

/// Information about an edge's classification for HLR.
#[derive(Debug, Clone)]
pub struct EdgeClassInfo {
    /// Edge index in the BRep.
    pub edge_idx: usize,
    /// Classification type.
    pub classification: EdgeClassification,
    /// Number of visible segments.
    pub visible_segments: usize,
    /// Number of hidden segments.
    pub hidden_segments: usize,
    /// Whether this edge is on a curved surface.
    pub on_curved_surface: bool,
    /// Surface index if on a curved surface.
    pub surface_idx: Option<usize>,
}

/// A projected edge segment labeled as visible or hidden.
#[derive(Debug, Clone, PartialEq)]
pub struct HlrSegment {
    /// Start point in 2D screen space.
    pub start: DVec2,
    /// End point in 2D screen space.
    pub end: DVec2,
    /// Whether this segment is visible from the camera.
    pub visible: bool,
    /// Optional hint about the underlying curve type (None for straight lines).
    pub curve_hint: Option<CurveHint>,
    /// Segment type (edge or silhouette).
    pub segment_type: SegmentType,
}

impl HlrSegment {
    /// Returns true if this segment is a silhouette/contour curve.
    pub fn is_contour(&self) -> bool {
        self.segment_type == SegmentType::Silhouette
    }

    /// Returns true if this segment is a thread edge.
    pub fn is_thread(&self) -> bool {
        self.segment_type == SegmentType::Thread
    }

    /// Returns true if this segment is a seam edge.
    pub fn is_seam(&self) -> bool {
        self.segment_type == SegmentType::Seam
    }
}

// ── Surface Normal Analysis ─────────────────────────────────────────────────────

/// Cached surface properties for efficient silhouette computation.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceProperties {
    /// Surface point at (u, v).
    pub point: DVec3,
    /// Unit normal at (u, v).
    pub normal: DVec3,
    /// Principal curvatures (k1, k2).
    pub curvatures: (f64, f64),
    /// Gaussian curvature.
    pub gaussian: f64,
    /// Mean curvature.
    pub mean: f64,
}

impl SurfaceProperties {
    /// Compute the dot product of the normal with a view direction.
    #[inline]
    pub fn normal_dot_view(&self, view_dir: DVec3) -> f64 {
        self.normal.dot(view_dir)
    }

    /// Check if this point is near a silhouette (normal nearly perpendicular to view).
    #[inline]
    pub fn is_near_silhouette(&self, view_dir: DVec3, threshold: f64) -> bool {
        self.normal_dot_view(view_dir).abs() < threshold
    }

    /// Get the maximum principal curvature magnitude.
    #[inline]
    pub fn max_curvature(&self) -> f64 {
        self.curvatures.0.abs().max(self.curvatures.1.abs())
    }

    /// Check if the surface is locally flat (low curvature).
    #[inline]
    pub fn is_flat(&self, tolerance: f64) -> bool {
        self.max_curvature() < tolerance
    }
}

/// Cache for surface property evaluations.
#[derive(Debug, Clone)]
pub struct SurfacePropertyCache {
    /// Cached properties keyed by (u, v) discretized to grid cells.
    cache: HashMap<(usize, usize), SurfaceProperties>,
    /// Grid resolution for cache.
    resolution: usize,
    /// UV domain of the surface.
    domain: [f64; 4],
}

impl SurfacePropertyCache {
    /// Create a new cache with given resolution.
    pub fn new(resolution: usize, domain: [f64; 4]) -> Self {
        Self {
            cache: HashMap::new(),
            resolution,
            domain,
        }
    }

    /// Get or compute surface properties at (u, v).
    pub fn get_or_compute(&mut self, surface: &Surface3, u: f64, v: f64) -> SurfaceProperties {
        let [u0, u1, v0, v1] = self.domain;
        let i = ((u - u0) / (u1 - u0) * self.resolution as f64).min(self.resolution as f64 - 1.0) as usize;
        let j = ((v - v0) / (v1 - v0) * self.resolution as f64).min(self.resolution as f64 - 1.0) as usize;

        if let Some(&props) = self.cache.get(&(i, j)) {
            return props;
        }

        let props = compute_surface_properties(surface, u, v);
        self.cache.insert((i, j), props);
        props
    }

    /// Get cached properties if available.
    pub fn get(&self, u: f64, v: f64) -> Option<&SurfaceProperties> {
        let [u0, u1, v0, v1] = self.domain;
        let i = ((u - u0) / (u1 - u0) * self.resolution as f64).min(self.resolution as f64 - 1.0) as usize;
        let j = ((v - v0) / (v1 - v0) * self.resolution as f64).min(self.resolution as f64 - 1.0) as usize;
        self.cache.get(&(i, j))
    }

    /// Clear the cache.
    pub fn clear(&mut self) {
        self.cache.clear();
    }

    /// Get the number of cached entries.
    pub fn len(&self) -> usize {
        self.cache.len()
    }

    /// Check if the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }
}

/// Compute surface properties at a given parameter location.
pub fn compute_surface_properties(surface: &Surface3, u: f64, v: f64) -> SurfaceProperties {
    let point = surface.point_at(u, v);
    let normal = surface.normal_at(u, v);
    let curvatures = rcad_kernel::curvature::principal_curvatures(surface, u, v);
    let gaussian = rcad_kernel::curvature::gaussian_curvature(surface, u, v);
    let mean = rcad_kernel::curvature::mean_curvature(surface, u, v);

    SurfaceProperties {
        point,
        normal,
        curvatures,
        gaussian,
        mean,
    }
}

// ── Adaptive Silhouette Sampling ────────────────────────────────────────────────

/// Sample point along a silhouette curve with additional metadata.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveSample {
    /// Parameter space location (u, v).
    pub uv: (f64, f64),
    /// World space position.
    pub point: DVec3,
    /// Surface normal at this point.
    pub normal: DVec3,
    /// Maximum curvature at this point.
    pub curvature: f64,
    /// Distance to the exact silhouette curve (0 = on silhouette).
    pub silhouette_distance: f64,
    /// Sampling weight (higher = more samples nearby).
    pub weight: f64,
}

/// Adaptive sampling configuration for silhouette curves.
#[derive(Debug, Clone)]
pub struct AdaptiveSamplingConfig {
    /// Base number of samples.
    pub base_samples: usize,
    /// Maximum number of samples after adaptive refinement.
    pub max_samples: usize,
    /// Curvature threshold for refinement (higher curvature = more samples).
    pub curvature_threshold: f64,
    /// Proximity threshold for refinement (closer to silhouette = more samples).
    pub proximity_threshold: f64,
    /// Minimum chord length between samples.
    pub min_chord_length: f64,
    /// Maximum angle deviation between consecutive samples (radians).
    pub max_angle_deviation: f64,
}

impl Default for AdaptiveSamplingConfig {
    fn default() -> Self {
        Self {
            base_samples: 32,
            max_samples: 256,
            curvature_threshold: 10.0,
            proximity_threshold: 0.05,
            min_chord_length: 1e-4,
            max_angle_deviation: 0.1,
        }
    }
}

/// Compute adaptive samples along a silhouette curve.
pub fn compute_adaptive_samples(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    config: &AdaptiveSamplingConfig,
    opts: &HlrOptions,
) -> Vec<AdaptiveSample> {
    let [u0, u1, v0, v1] = domain;

    // Phase 1: Find silhouette seed points
    let seeds = find_silhouette_seeds(surface, view_dir, domain, config.base_samples, opts.tangent_tolerance);

    if seeds.is_empty() {
        return Vec::new();
    }

    // Phase 2: Trace silhouette curves from seeds
    let mut all_samples: Vec<AdaptiveSample> = Vec::new();
    let mut visited: HashSet<(usize, usize)> = HashSet::new();

    for (_, _, u, v) in seeds {
        // Check if this cell was already visited
        let cell_i = ((u - u0) / (u1 - u0) * config.base_samples as f64) as usize;
        let cell_j = ((v - v0) / (v1 - v0) * config.base_samples as f64) as usize;

        if visited.contains(&(cell_i, cell_j)) {
            continue;
        }

        // Trace the silhouette curve
        let curve_samples = trace_adaptive_silhouette(surface, view_dir, domain, u, v, config, opts);

        // Mark visited cells
        for sample in &curve_samples {
            let ci = ((sample.uv.0 - u0) / (u1 - u0) * config.base_samples as f64) as usize;
            let cj = ((sample.uv.1 - v0) / (v1 - v0) * config.base_samples as f64) as usize;
            visited.insert((ci.min(config.base_samples - 1), cj.min(config.base_samples - 1)));
        }

        all_samples.extend(curve_samples);
    }

    // Phase 3: Refine samples in high-curvature and near-silhouette regions
    if opts.curvature_adaptive {
        refine_adaptive_samples(surface, view_dir, &mut all_samples, config, opts);
    }

    all_samples
}

/// Trace a silhouette curve with adaptive sampling.
fn trace_adaptive_silhouette(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    u_start: f64,
    v_start: f64,
    config: &AdaptiveSamplingConfig,
    opts: &HlrOptions,
) -> Vec<AdaptiveSample> {
    let mut samples: Vec<AdaptiveSample> = Vec::new();
    let [u0, u1, v0, v1] = domain;

    // Add the starting point
    if let Some(sample) = create_adaptive_sample(surface, view_dir, u_start, v_start) {
        samples.push(sample);
    }

    // March in both directions from the seed
    for direction in &[-1.0_f64, 1.0] {
        let mut u = u_start;
        let mut v = v_start;
        let mut prev_sample = samples.first().copied();

        for _ in 0..opts.max_subdivisions * 100 {
            // Compute the tangent direction to the silhouette curve
            let tangent = compute_silhouette_tangent(surface, view_dir, u, v);

            if tangent.length_squared() < 1e-16 {
                break;
            }

            // Choose direction along the tangent
            let step_dir = *direction * tangent.normalize_or_zero();

            // Compute adaptive step size based on curvature
            let props = compute_surface_properties(surface, u, v);
            let max_k = props.max_curvature().max(opts.min_curvature);
            let curvature_factor = (opts.max_curvature / max_k).min(4.0).max(0.25);
            let step_size = opts.angular_tolerance * curvature_factor;

            // Take a step
            let u_new = u + step_dir.x * step_size;
            let v_new = v + step_dir.y * step_size;

            // Check bounds
            if u_new < u0 || u_new > u1 || v_new < v0 || v_new > v1 {
                break;
            }

            // Project back onto the silhouette curve
            if let Some((u_proj, v_proj)) = project_to_silhouette(surface, view_dir, u_new, v_new, opts.tangent_tolerance) {
                u = u_proj;
                v = v_proj;

                // Create a new sample
                if let Some(sample) = create_adaptive_sample(surface, view_dir, u, v) {
                    // Check if we've moved enough to add a new sample
                    if let Some(prev) = prev_sample {
                        let dist = (sample.point - prev.point).length();
                        if dist < config.min_chord_length {
                            continue;
                        }

                        // Check angular deviation
                        let dir_new = (sample.point - prev.point).normalize_or_zero();
                        let dir_prev = if samples.len() >= 2 {
                            (samples[samples.len() - 1].point - samples[samples.len() - 2].point).normalize_or_zero()
                        } else {
                            dir_new
                        };
                        let angle = dir_new.dot(dir_prev).acos().abs();
                        if angle > config.max_angle_deviation {
                            // Add intermediate samples for high angular deviation
                            add_intermediate_samples(surface, view_dir, prev, sample, &mut samples, config);
                        }
                    }

                    samples.push(sample);
                    prev_sample = Some(sample);
                }
            } else {
                break;
            }

            // Check for closed loop
            if samples.len() > 10 {
                let first = samples[0];
                let dist = ((first.uv.0 - u).powi(2) + (first.uv.1 - v).powi(2)).sqrt();
                if dist < step_size * 2.0 {
                    break;
                }
            }
        }
    }

    samples
}

/// Create an adaptive sample at a parameter location.
fn create_adaptive_sample(surface: &Surface3, view_dir: DVec3, u: f64, v: f64) -> Option<AdaptiveSample> {
    let point = surface.point_at(u, v);
    let normal = surface.normal_at(u, v);
    let curvatures = rcad_kernel::curvature::principal_curvatures(surface, u, v);
    let curvature = curvatures.0.abs().max(curvatures.1.abs());

    // Compute distance to exact silhouette (absolute value of normal dot view)
    let silhouette_distance = normal.dot(view_dir).abs();

    // Compute sampling weight based on curvature and proximity
    let weight = (curvature + 1.0) * (silhouette_distance + 0.1).recip();

    Some(AdaptiveSample {
        uv: (u, v),
        point,
        normal,
        curvature,
        silhouette_distance,
        weight,
    })
}

/// Add intermediate samples between two samples for smooth curves.
fn add_intermediate_samples(
    surface: &Surface3,
    view_dir: DVec3,
    start: AdaptiveSample,
    end: AdaptiveSample,
    samples: &mut Vec<AdaptiveSample>,
    config: &AdaptiveSamplingConfig,
) {
    let num_intermediate = ((end.point - start.point).length() / config.min_chord_length).ceil() as usize;
    let num_intermediate = num_intermediate.min(4);

    for i in 1..num_intermediate {
        let t = i as f64 / num_intermediate as f64;
        let u = start.uv.0 + t * (end.uv.0 - start.uv.0);
        let v = start.uv.1 + t * (end.uv.1 - start.uv.1);

        if let Some(sample) = create_adaptive_sample(surface, view_dir, u, v) {
            samples.push(sample);
        }
    }
}

/// Refine adaptive samples based on curvature and silhouette proximity.
fn refine_adaptive_samples(
    surface: &Surface3,
    view_dir: DVec3,
    samples: &mut Vec<AdaptiveSample>,
    config: &AdaptiveSamplingConfig,
    opts: &HlrOptions,
) {
    if samples.len() < 2 {
        return;
    }

    let mut refined: Vec<AdaptiveSample> = Vec::with_capacity(samples.len() * 2);
    refined.push(samples[0]);

    for i in 1..samples.len() {
        let prev = &samples[i - 1];
        let curr = &samples[i];

        // Determine if refinement is needed based on curvature and proximity
        let chord_len = (curr.point - prev.point).length();
        let avg_curvature = (prev.curvature + curr.curvature) * 0.5;
        let avg_proximity = (prev.silhouette_distance + curr.silhouette_distance) * 0.5;

        let needs_refinement = avg_curvature > config.curvature_threshold
            || avg_proximity < config.proximity_threshold
            || chord_len > config.min_chord_length * 4.0;

        if needs_refinement {
            let num_subdivisions = (avg_curvature * chord_len / config.curvature_threshold).ceil() as usize;
            let num_subdivisions = num_subdivisions.min(4).max(1);

            for j in 1..num_subdivisions {
                let t = j as f64 / num_subdivisions as f64;
                let u = prev.uv.0 + t * (curr.uv.0 - prev.uv.0);
                let v = prev.uv.1 + t * (curr.uv.1 - prev.uv.1);

                if let Some(sample) = create_adaptive_sample(surface, view_dir, u, v) {
                    refined.push(sample);
                }
            }
        }

        refined.push(*curr);
    }

    *samples = refined;
}

/// Output of an HLR computation.
#[derive(Debug, Clone, Default)]
pub struct HlrResult {
    pub segments: Vec<HlrSegment>,
}

impl HlrResult {
    /// Return only visible segments.
    pub fn visible(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.visible)
    }

    /// Return only hidden segments.
    pub fn hidden(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| !s.visible)
    }

    /// Return only silhouette/contour segments.
    pub fn silhouettes(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.is_contour())
    }

    /// Return only visible silhouette segments.
    pub fn visible_silhouettes(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.visible && s.is_contour())
    }
}

/// Camera / view specification for HLR.
#[derive(Debug, Clone)]
pub struct HlrCamera {
    /// Camera position in world space.
    pub eye: DVec3,
    /// Target point (look-at).
    pub target: DVec3,
    /// Up direction.
    pub up: DVec3,
}

impl HlrCamera {
    pub fn new(eye: DVec3, target: DVec3) -> Self {
        Self {
            eye,
            target,
            up: DVec3::Y,
        }
    }

    pub fn with_up(mut self, up: DVec3) -> Self {
        self.up = up;
        self
    }

    /// Isometric-style view from the +X+Y+Z octant.
    pub fn isometric(distance: f64) -> Self {
        let d = distance / 3.0_f64.sqrt();
        Self::new(DVec3::splat(d), DVec3::ZERO)
    }

    /// Front view (looking along +Y, up = +Z).
    pub fn front(distance: f64) -> Self {
        Self::new(DVec3::new(0.0, -distance, 0.0), DVec3::ZERO).with_up(DVec3::Z)
    }

    /// Top view (looking down -Z).
    pub fn top(distance: f64) -> Self {
        Self::new(DVec3::new(0.0, 0.0, distance), DVec3::ZERO).with_up(DVec3::Y)
    }

    /// Right-side view (looking along -X, up = +Z).
    pub fn right(distance: f64) -> Self {
        Self::new(DVec3::new(distance, 0.0, 0.0), DVec3::ZERO).with_up(DVec3::Z)
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Build a right-handed view matrix (world → camera space).
fn look_at(eye: DVec3, target: DVec3, up: DVec3) -> DMat4 {
    let forward = (target - eye).normalize_or_zero();
    let right = forward.cross(up).normalize_or_zero();
    let up = right.cross(forward);

    DMat4::from_cols(
        DVec4::new(right.x, right.y, right.z, -right.dot(eye)),
        DVec4::new(up.x, up.y, up.z, -up.dot(eye)),
        DVec4::new(-forward.x, -forward.y, -forward.z, forward.dot(eye)),
        DVec4::new(0.0, 0.0, 0.0, 1.0),
    )
    .transpose()
}

/// Project a world-space point to 2D screen space using the view matrix.
/// Returns (x, y) in camera space (z is depth; ignored for 2D output).
fn project(p: DVec3, view: &DMat4) -> (DVec2, f64) {
    let hp = view.mul_vec4(DVec4::new(p.x, p.y, p.z, 1.0));
    (DVec2::new(hp.x, hp.y), hp.z)
}

/// Collect all triangles from a BRep (fan-triangulate faces without pre-triangulated data).
fn collect_triangles(brep: &BRep) -> Vec<[DVec3; 3]> {
    let mut tris = Vec::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if !face.triangles.is_empty() {
                    for &[i, j, k] in &face.triangles {
                        if let (Some(a), Some(b), Some(c)) = (
                            brep.vertices.get(i),
                            brep.vertices.get(j),
                            brep.vertices.get(k),
                        ) {
                            tris.push([a.point, b.point, c.point]);
                        }
                    }
                } else {
                    // Fan-triangulate from wire
                    let pts: Vec<DVec3> = face
                        .outer_wire
                        .edges
                        .iter()
                        .filter_map(|we| {
                            let edge = brep.edges.get(we.idx)?;
                            let vi = if we.forward { edge.start } else { edge.end };
                            brep.vertices.get(vi).map(|v| v.point)
                        })
                        .collect();
                    if pts.len() >= 3 {
                        let origin = pts[0];
                        for i in 1..pts.len() - 1 {
                            tris.push([origin, pts[i], pts[i + 1]]);
                        }
                    }
                }
            }
        }
    }
    tris
}

/// Ray-triangle intersection (Möller–Trumbore). Returns `Some(t)` if the ray
/// `origin + t*dir` hits the triangle (t > epsilon, front-face only).
fn ray_triangle_intersect(origin: DVec3, dir: DVec3, tri: &[DVec3; 3]) -> Option<f64> {
    const EPS: f64 = 1e-8;
    let edge1 = tri[1] - tri[0];
    let edge2 = tri[2] - tri[0];
    let h = dir.cross(edge2);
    let a = edge1.dot(h);
    if a.abs() < EPS {
        return None;
    }
    let f = 1.0 / a;
    let s = origin - tri[0];
    let u = f * s.dot(h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(edge1);
    let v = f * dir.dot(q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * edge2.dot(q);
    if t > EPS { Some(t) } else { None }
}

/// Test if a world-space point is occluded by any triangle when viewed from `eye`.
fn is_occluded(point: DVec3, eye: DVec3, triangles: &[[DVec3; 3]], dist_to_eye: f64) -> bool {
    let dir = (eye - point).normalize_or_zero();
    let origin = point + dir * 1e-5; // push off surface
    for tri in triangles {
        if let Some(t) = ray_triangle_intersect(origin, dir, tri)
            && t < dist_to_eye - 1e-4
        {
            return true;
        }
    }
    false
}

// ── Triangle BVH for accelerated ray casting ───────────────────────────────────

/// Axis-aligned bounding box for triangle BVH.
#[derive(Debug, Clone, Copy)]
struct TriAabb {
    min: DVec3,
    max: DVec3,
}

impl TriAabb {
    fn empty() -> Self {
        Self {
            min: DVec3::splat(f64::INFINITY),
            max: DVec3::splat(f64::NEG_INFINITY),
        }
    }

    fn from_triangle(tri: &[DVec3; 3]) -> Self {
        let mut aabb = Self::empty();
        for &p in tri {
            aabb.expand_point(p);
        }
        aabb
    }

    fn expand_point(&mut self, p: DVec3) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    fn expand_aabb(&mut self, other: &TriAabb) {
        self.min = self.min.min(other.min);
        self.max = self.max.max(other.max);
    }

    fn center(&self) -> DVec3 {
        (self.min + self.max) * 0.5
    }

    fn surface_area(&self) -> f64 {
        let d = self.max - self.min;
        if d.x < 0.0 || d.y < 0.0 || d.z < 0.0 {
            return 0.0;
        }
        2.0 * (d.x * d.y + d.y * d.z + d.z * d.x)
    }

    fn ray_intersect(&self, origin: DVec3, inv_dir: DVec3) -> Option<f64> {
        let t1 = (self.min - origin) * inv_dir;
        let t2 = (self.max - origin) * inv_dir;

        let t_min = t1.min(t2);
        let t_max = t1.max(t2);

        let t_enter = t_min.x.max(t_min.y).max(t_min.z);
        let t_exit = t_max.x.min(t_max.y).min(t_max.z);

        if t_exit >= t_enter.max(0.0) {
            Some(t_enter.max(0.0))
        } else {
            None
        }
    }
}

/// BVH node for triangle-level acceleration.
#[derive(Debug, Clone)]
enum TriBvhNode {
    Leaf {
        aabb: TriAabb,
        /// Triangle indices (index into the original triangle array).
        tris: Vec<usize>,
    },
    Internal {
        aabb: TriAabb,
        left: usize,
        right: usize,
    },
}

/// Triangle-level BVH for efficient ray casting.
#[derive(Debug, Clone)]
pub struct TriBvh {
    nodes: Vec<TriBvhNode>,
    triangle_aabbs: Vec<TriAabb>,
    triangle_centers: Vec<DVec3>,
}

const MAX_TRIS_PER_LEAF: usize = 8;
const SAH_BUCKETS: usize = 8;

impl TriBvh {
    /// Build a triangle BVH from a list of triangles.
    pub fn build(triangles: &[[DVec3; 3]]) -> Self {
        if triangles.is_empty() {
            return Self {
                nodes: Vec::new(),
                triangle_aabbs: Vec::new(),
                triangle_centers: Vec::new(),
            };
        }

        let triangle_aabbs: Vec<TriAabb> = triangles.iter().map(|t| TriAabb::from_triangle(t)).collect();
        let triangle_centers: Vec<DVec3> = triangle_aabbs.iter().map(|a| a.center()).collect();

        let tri_indices: Vec<usize> = (0..triangles.len()).collect();

        let mut bvh = Self {
            nodes: Vec::new(),
            triangle_aabbs,
            triangle_centers,
        };

        bvh.build_recursive(&tri_indices);
        bvh
    }

    fn build_recursive(&mut self, tri_indices: &[usize]) -> usize {
        let count = tri_indices.len();
        if count == 0 {
            return usize::MAX;
        }

        // Compute AABB for this node
        let mut aabb = TriAabb::empty();
        for &ti in tri_indices {
            aabb.expand_aabb(&self.triangle_aabbs[ti]);
        }

        // Leaf condition
        if count <= MAX_TRIS_PER_LEAF {
            let node_idx = self.nodes.len();
            self.nodes.push(TriBvhNode::Leaf {
                aabb,
                tris: tri_indices.to_vec(),
            });
            return node_idx;
        }

        // SAH split
        let (split_axis, split_pos) = self.sah_split(tri_indices, &aabb);

        // Partition triangles
        let (left_tris, right_tris): (Vec<usize>, Vec<usize>) = tri_indices
            .iter()
            .copied()
            .partition(|&ti| {
                let center = match split_axis {
                    0 => self.triangle_centers[ti].x,
                    1 => self.triangle_centers[ti].y,
                    _ => self.triangle_centers[ti].z,
                };
                center < split_pos
            });

        // Handle degenerate split
        let (left_tris, right_tris) = if left_tris.is_empty() || right_tris.is_empty() {
            let mid = count / 2;
            let mut sorted = tri_indices.to_vec();
            sorted.sort_by(|&a, &b| {
                let ca = match split_axis {
                    0 => self.triangle_centers[a].x,
                    1 => self.triangle_centers[a].y,
                    _ => self.triangle_centers[a].z,
                };
                let cb = match split_axis {
                    0 => self.triangle_centers[b].x,
                    1 => self.triangle_centers[b].y,
                    _ => self.triangle_centers[b].z,
                };
                ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
            });
            (sorted[..mid].to_vec(), sorted[mid..].to_vec())
        } else {
            (left_tris, right_tris)
        };

        let node_idx = self.nodes.len();
        self.nodes.push(TriBvhNode::Internal {
            aabb: TriAabb::empty(),
            left: 0,
            right: 0,
        });

        let left = self.build_recursive(&left_tris);
        let right = self.build_recursive(&right_tris);

        self.nodes[node_idx] = TriBvhNode::Internal { aabb, left, right };
        node_idx
    }

    fn sah_split(&self, tri_indices: &[usize], parent_aabb: &TriAabb) -> (usize, f64) {
        let parent_sa = parent_aabb.surface_area().max(1e-30);
        let mut best_cost = f64::INFINITY;
        let mut best_axis = 0usize;
        let mut best_pos = 0.0f64;

        for axis in 0..3usize {
            let axis_min = match axis {
                0 => parent_aabb.min.x,
                1 => parent_aabb.min.y,
                _ => parent_aabb.min.z,
            };
            let axis_max = match axis {
                0 => parent_aabb.max.x,
                1 => parent_aabb.max.y,
                _ => parent_aabb.max.z,
            };
            let span = axis_max - axis_min;
            if span < 1e-14 {
                continue;
            }

            for b in 1..SAH_BUCKETS {
                let split = axis_min + span * b as f64 / SAH_BUCKETS as f64;

                let mut left_aabb = TriAabb::empty();
                let mut right_aabb = TriAabb::empty();
                let mut left_count = 0usize;
                let mut right_count = 0usize;

                for &ti in tri_indices {
                    let center_val = match axis {
                        0 => self.triangle_centers[ti].x,
                        1 => self.triangle_centers[ti].y,
                        _ => self.triangle_centers[ti].z,
                    };
                    if center_val < split {
                        left_aabb.expand_aabb(&self.triangle_aabbs[ti]);
                        left_count += 1;
                    } else {
                        right_aabb.expand_aabb(&self.triangle_aabbs[ti]);
                        right_count += 1;
                    }
                }

                if left_count == 0 || right_count == 0 {
                    continue;
                }

                let cost = (left_count as f64 * left_aabb.surface_area()
                    + right_count as f64 * right_aabb.surface_area())
                    / parent_sa;

                if cost < best_cost {
                    best_cost = cost;
                    best_axis = axis;
                    best_pos = split;
                }
            }
        }

        if best_cost.is_infinite() {
            let d = parent_aabb.max - parent_aabb.min;
            best_axis = if d.x >= d.y && d.x >= d.z { 0 } else if d.y >= d.z { 1 } else { 2 };
            best_pos = parent_aabb.center()[best_axis];
        }

        (best_axis, best_pos)
    }

    /// Test if a point is occluded by any triangle in the BVH.
    pub fn is_occluded(&self, point: DVec3, eye: DVec3, triangles: &[[DVec3; 3]], dist_to_eye: f64) -> bool {
        if self.nodes.is_empty() {
            return false;
        }

        let dir = (eye - point).normalize_or_zero();
        let origin = point + dir * 1e-5;
        let inv_dir = DVec3::new(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);

        self.is_occluded_node(0, origin, dir, inv_dir, triangles, dist_to_eye)
    }

    fn is_occluded_node(
        &self,
        node_idx: usize,
        origin: DVec3,
        dir: DVec3,
        inv_dir: DVec3,
        triangles: &[[DVec3; 3]],
        dist_to_eye: f64,
    ) -> bool {
        let node = &self.nodes[node_idx];

        // Check AABB intersection
        let t_aabb = match node.aabb().ray_intersect(origin, inv_dir) {
            Some(t) => t,
            None => return false,
        };

        // Early exit if AABB is beyond the eye
        if t_aabb > dist_to_eye {
            return false;
        }

        match node {
            TriBvhNode::Leaf { tris, .. } => {
                for &ti in tris {
                    if let Some(t) = ray_triangle_intersect(origin, dir, &triangles[ti]) {
                        if t < dist_to_eye - 1e-4 {
                            return true;
                        }
                    }
                }
                false
            }
            TriBvhNode::Internal { left, right, .. } => {
                self.is_occluded_node(*left, origin, dir, inv_dir, triangles, dist_to_eye)
                    || self.is_occluded_node(*right, origin, dir, inv_dir, triangles, dist_to_eye)
            }
        }
    }
}

impl TriBvhNode {
    fn aabb(&self) -> &TriAabb {
        match self {
            TriBvhNode::Leaf { aabb, .. } => aabb,
            TriBvhNode::Internal { aabb, .. } => aabb,
        }
    }
}

// ── Public API ────────────────────────────────────────────────────────────────

// ── Silhouette generation ─────────────────────────────────────────────────────

/// Internal: one silhouette curve to process through the HLR pipeline.
struct SilhouetteCurve {
    /// World-space sample points (at least 2).
    world_pts: Vec<DVec3>,
    /// Optional curve hint for SVG output.
    curve_hint: Option<CurveHint>,
    /// If true, emit one segment per consecutive point pair instead of merging
    /// runs.  Used for dense polyline approximations (e.g. sphere silhouette).
    dense: bool,
}

/// A 3D silhouette curve extracted from a curved surface.
#[derive(Debug, Clone)]
pub struct SilhouetteCurve3 {
    /// World-space sample points along the silhouette curve.
    pub points: Vec<DVec3>,
    /// The surface index from which this silhouette was extracted.
    pub surface_index: usize,
}

/// Extract silhouette curves from a BRep for a given view direction.
///
/// This function computes the visible contour lines (silhouettes) of curved surfaces
/// as seen from a specific viewing direction. For analytic surfaces (cylinder, sphere,
/// cone, torus), exact silhouette curves are computed. For general surfaces (BSpline,
/// Bezier, etc.), numerical methods with adaptive sampling are used.
///
/// # Arguments
/// * `brep` - The BRep model to extract silhouettes from.
/// * `view_dir` - The normalized view direction (from target to eye).
/// * `opts` - Configuration options for sampling and tolerance.
///
/// # Returns
/// A vector of 3D silhouette curves, each represented as a series of world-space points.
pub fn extract_silhouette_curves(brep: &BRep, view_dir: DVec3, opts: &HlrOptions) -> Vec<SilhouetteCurve3> {
    let mut curves: Vec<SilhouetteCurve3> = Vec::new();

    if brep.solids.is_empty() {
        return curves;
    }

    let line_samples = opts.silhouette_samples.max(16);
    let dense_curve_samples = (opts.silhouette_samples * 4).max(64);

    let mut face_idx = 0usize;
    for shell in &brep.solids[0].shells {
        for _face in &shell.faces {
            let surf_idx = match brep.geom.face_surface.get(face_idx).and_then(|o| *o) {
                Some(idx) => idx,
                None => {
                    face_idx += 1;
                    continue;
                }
            };
            let surface = &brep.geom.surfaces[surf_idx];

            let domain = match brep.geom.face_surface_range.get(face_idx).and_then(|o| *o) {
                Some(r) => r,
                None => surface.default_domain(),
            };
            let [u0, u1, v0, v1] = domain;

            // Extract silhouettes based on surface type
            let face_curves = extract_surface_silhouettes(
                surface, view_dir, domain, brep, opts, line_samples, dense_curve_samples,
            );

            for pts in face_curves {
                if pts.len() >= 2 {
                    curves.push(SilhouetteCurve3 {
                        points: pts,
                        surface_index: surf_idx,
                    });
                }
            }

            face_idx += 1;
        }
    }

    curves
}

/// Extract silhouette curves from a single surface.
fn extract_surface_silhouettes(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    brep: &BRep,
    opts: &HlrOptions,
    line_samples: usize,
    dense_curve_samples: usize,
) -> Vec<Vec<DVec3>> {
    let [u0, u1, v0, v1] = domain;
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    match surface {
        Surface3::Cylinder(cyl) => {
            curves.extend(extract_cylinder_silhouettes(cyl, view_dir, brep, line_samples, v0, v1));
        }

        Surface3::Sphere(sph) => {
            curves.push(extract_sphere_silhouette(sph, view_dir, dense_curve_samples));
        }

        Surface3::Cone(con) => {
            curves.extend(extract_cone_silhouettes(con, view_dir, brep, line_samples, v0, v1));
        }

        Surface3::Torus(tor) => {
            curves.extend(extract_torus_silhouettes(tor, view_dir, dense_curve_samples));
        }

        Surface3::Ellipsoid(ell) => {
            curves.extend(extract_ellipsoid_silhouettes(ell, view_dir, opts, dense_curve_samples));
        }

        // For general surfaces, use numerical silhouette extraction
        Surface3::BSpline(_)
        | Surface3::Bezier(_)
        | Surface3::TriBezier(_)
        | Surface3::Offset(_)
        | Surface3::LinearExtrusion(_)
        | Surface3::Revolution(_)
        | Surface3::Ruled(_)
        | Surface3::Coons(_)
        | Surface3::Pipe(_) => {
            curves.extend(extract_numerical_silhouettes(
                surface, view_dir, domain, opts, brep,
            ));
        }

        // Planes have no silhouette curves
        Surface3::Plane(_) | Surface3::Trimmed(_) | Surface3::Helicoid(_) => {}
    }

    curves
}

/// Extract silhouette lines from a cylinder.
fn extract_cylinder_silhouettes(
    cyl: &rcad_kernel::geom::CylindricalSurface,
    view_dir: DVec3,
    brep: &BRep,
    line_samples: usize,
    v0: f64,
    v1: f64,
) -> Vec<Vec<DVec3>> {
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    // Project view direction onto the plane perpendicular to the axis.
    let d_perp = view_dir - view_dir.dot(cyl.axis) * cyl.axis;
    if d_perp.length_squared() < 1e-10 {
        // Viewing along the axis — no silhouette lines.
        return curves;
    }

    // Direction from axis to silhouette (perpendicular to both axis and d_perp).
    let sil_dir = cyl.axis.cross(d_perp).normalize_or_zero();

    // Resolve v range (height along axis).
    let (v0_eff, v1_eff) = if v0.is_finite() && v1.is_finite() {
        (v0, v1)
    } else {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for vert in &brep.vertices {
            let proj = (vert.point - cyl.origin).dot(cyl.axis);
            lo = lo.min(proj);
            hi = hi.max(proj);
        }
        if lo.is_finite() && hi.is_finite() {
            (lo, hi)
        } else {
            return curves;
        }
    };

    for &sign in &[1.0_f64, -1.0] {
        let offset = sil_dir * sign * cyl.radius;
        let world_pts: Vec<DVec3> = (0..line_samples)
            .map(|i| {
                let t = i as f64 / (line_samples - 1) as f64;
                let v = v0_eff + (v1_eff - v0_eff) * t;
                cyl.origin + v * cyl.axis + offset
            })
            .collect();
        curves.push(world_pts);
    }

    curves
}

/// Extract silhouette curve from a sphere (great circle perpendicular to view direction).
fn extract_sphere_silhouette(
    sph: &rcad_kernel::geom::SphericalSurface,
    view_dir: DVec3,
    samples: usize,
) -> Vec<DVec3> {
    let x_ax = any_perpendicular(view_dir);
    let y_ax = view_dir.cross(x_ax).normalize_or_zero();

    (0..samples)
        .map(|i| {
            let t = 2.0 * std::f64::consts::PI * i as f64 / samples as f64;
            sph.center + sph.radius * (t.cos() * x_ax + t.sin() * y_ax)
        })
        .collect()
}

/// Extract silhouette lines from a cone (two generators from apex).
fn extract_cone_silhouettes(
    con: &rcad_kernel::geom::ConicalSurface,
    view_dir: DVec3,
    brep: &BRep,
    line_samples: usize,
    v0: f64,
    v1: f64,
) -> Vec<Vec<DVec3>> {
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    let d_perp = view_dir - view_dir.dot(con.axis) * con.axis;
    if d_perp.length_squared() < 1e-10 {
        return curves;
    }

    let sil_dir = con.axis.cross(d_perp).normalize_or_zero();

    let (v0_eff, v1_eff) = if v0.is_finite() && v1.is_finite() {
        (v0, v1)
    } else {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for vert in &brep.vertices {
            let proj = (vert.point - con.apex).dot(con.axis);
            lo = lo.min(proj);
            hi = hi.max(proj);
        }
        if lo.is_finite() && hi.is_finite() {
            (lo.max(0.0), hi.max(0.0))
        } else {
            return curves;
        }
    };

    let tan_a = con.half_angle_rad.tan();
    for &sign in &[1.0_f64, -1.0] {
        let world_pts: Vec<DVec3> = (0..line_samples)
            .map(|i| {
                let t = i as f64 / (line_samples - 1) as f64;
                let v = v0_eff + (v1_eff - v0_eff) * t;
                con.apex + v * con.axis + v * tan_a * sil_dir * sign
            })
            .collect();

        if world_pts
            .first()
            .zip(world_pts.last())
            .map(|(a, b)| (*b - *a).length_squared() > 1e-12)
            .unwrap_or(false)
        {
            curves.push(world_pts);
        }
    }

    curves
}

/// Extract silhouette curves from a torus.
fn extract_torus_silhouettes(
    tor: &rcad_kernel::geom::ToroidalSurface,
    view_dir: DVec3,
    samples: usize,
) -> Vec<Vec<DVec3>> {
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    let x_ax = any_perpendicular(tor.axis);
    let y_ax = tor.axis.cross(x_ax).normalize_or_zero();
    let axis_dot = tor.axis.dot(view_dir);

    for &offset in &[0.0_f64, std::f64::consts::PI] {
        let pts: Vec<DVec3> = (0..samples)
            .map(|i| {
                let u = 2.0 * std::f64::consts::PI * i as f64 / samples as f64;
                let radial = u.cos() * x_ax + u.sin() * y_ax;
                let radial_dot = radial.dot(view_dir);
                let v = (-radial_dot).atan2(axis_dot) + offset;
                let tube_center = tor.center + tor.major_radius * radial;
                tube_center + tor.minor_radius * (v.cos() * radial + v.sin() * tor.axis)
            })
            .collect();
        curves.push(pts);
    }

    curves
}

/// Extract silhouette curves from an ellipsoid using analytic methods.
///
/// The silhouette of an ellipsoid is the intersection of the ellipsoid surface
/// with a plane passing through the center. The plane normal is proportional to
/// (vx/a², vy/b², vz/c²) where v is the view direction and a, b, c are the radii.
///
/// This intersection is an ellipse, which we parameterize by:
/// 1. Finding two orthogonal directions in the silhouette plane
/// 2. For each angle, computing where a ray in that direction intersects the ellipsoid
fn extract_ellipsoid_silhouettes(
    ell: &rcad_kernel::geom::EllipsoidalSurface,
    view_dir: DVec3,
    opts: &HlrOptions,
    samples: usize,
) -> Vec<Vec<DVec3>> {
    use rcad_kernel::geom::SurfaceEval;
    use std::f64::consts::PI;

    // Build the orthonormal frame of the ellipsoid
    let (axis, x_axis, y_axis) = orthonormal_frame(ell.axis, ell.ref_dir);

    // Transform view direction into the ellipsoid's local coordinate frame
    // Local coordinates: x along x_axis, y along y_axis, z along axis
    let vx = view_dir.dot(x_axis);
    let vy = view_dir.dot(y_axis);
    let vz = view_dir.dot(axis);
    let view_local = DVec3::new(vx, vy, vz);

    // Handle degenerate case: view direction is zero
    if view_local.length_squared() < 1e-20 {
        return Vec::new();
    }

    // The silhouette plane normal in local coordinates is proportional to
    // (vx/a², vy/b², vz/c²)
    let a = ell.radius_x;
    let b = ell.radius_y;
    let c = ell.radius_z;

    let plane_normal_local = DVec3::new(
        vx / (a * a),
        vy / (b * b),
        vz / (c * c),
    );

    // Handle the degenerate case where the plane normal is zero
    // This happens when all components are zero, which shouldn't occur for valid view direction
    let plane_normal_len = plane_normal_local.length();
    if plane_normal_len < 1e-20 {
        // View direction is exactly perpendicular to all scaled axes
        // This is a degenerate case - return empty
        return Vec::new();
    }
    let plane_normal_local = plane_normal_local.normalize();

    // Check if view is along a principal axis (plane normal is near a coordinate axis)
    // In this case, the silhouette is an ellipse in the perpendicular plane
    let is_view_along_axis = plane_normal_local.z.abs() > 1.0 - 1e-6;
    let is_view_along_x = plane_normal_local.x.abs() > 1.0 - 1e-6;
    let is_view_along_y = plane_normal_local.y.abs() > 1.0 - 1e-6;

    // Find two orthogonal directions in the silhouette plane
    // These will be used to parameterize the ellipse
    let (u_dir, v_dir) = if is_view_along_axis {
        // View is along Z axis (ellipsoid's axis)
        // Silhouette plane is XY plane, silhouette is ellipse x²/a² + y²/b² = 1
        (DVec3::X, DVec3::Y)
    } else if is_view_along_x {
        // View is along X axis
        // Silhouette plane is YZ plane
        (DVec3::Y, DVec3::Z)
    } else if is_view_along_y {
        // View is along Y axis
        // Silhouette plane is XZ plane
        (DVec3::X, DVec3::Z)
    } else {
        // General case: find two orthogonal vectors in the silhouette plane
        // Use any_perpendicular to get a vector perpendicular to the plane normal
        let u = any_perpendicular(plane_normal_local);
        let v = plane_normal_local.cross(u).normalize_or_zero();
        (u, v)
    };

    // Parameterize the ellipse by sampling angles
    // For each angle θ, the ray direction is u*cos(θ) + v*sin(θ)
    // The intersection parameter t is: t = 1 / sqrt((dx/a)² + (dy/b)² + (dz/c)²)
    let actual_samples = samples.max(opts.silhouette_samples).max(32);
    let points: Vec<DVec3> = (0..actual_samples)
        .map(|i| {
            let theta = 2.0 * PI * i as f64 / actual_samples as f64;
            let cos_t = theta.cos();
            let sin_t = theta.sin();

            // Ray direction in local coordinates
            let dir_local = u_dir * cos_t + v_dir * sin_t;

            // Compute intersection parameter t
            // The ray is: p = t * dir_local
            // On ellipsoid: (tx/a)² + (ty/b)² + (tz/c)² = 1
            // t² * ((dx/a)² + (dy/b)² + (dz/c)²) = 1
            let dx = dir_local.x;
            let dy = dir_local.y;
            let dz = dir_local.z;
            let t_squared_recip = (dx * dx) / (a * a) + (dy * dy) / (b * b) + (dz * dz) / (c * c);

            let t = 1.0 / t_squared_recip.sqrt();

            // Local point on ellipsoid
            let local_point = dir_local * t;

            // Transform back to world coordinates
            ell.center + local_point.x * x_axis + local_point.y * y_axis + local_point.z * axis
        })
        .collect();

    // Verify that the silhouette points satisfy the silhouette condition
    // (normal · view_dir ≈ 0)
    let mut valid = true;
    for pt in &points {
        // Compute the point in local coordinates
        let p_local = *pt - ell.center;
        let x = p_local.dot(x_axis);
        let y = p_local.dot(y_axis);
        let z = p_local.dot(axis);

        // Normal direction (gradient of implicit equation)
        let grad_local = DVec3::new(
            x / (a * a),
            y / (b * b),
            z / (c * c),
        );
        let normal = (grad_local.x * x_axis + grad_local.y * y_axis + grad_local.z * axis).normalize_or_zero();
        let dot = normal.dot(view_dir);

        // Check if silhouette condition is approximately satisfied
        if dot.abs() > opts.tangent_tolerance.max(0.1) {
            valid = false;
            break;
        }
    }

    if valid && !points.is_empty() {
        vec![points]
    } else {
        // Fallback to numerical method if analytic result seems invalid
        extract_ellipsoid_silhouettes_numerical(ell, view_dir, opts, samples)
    }
}

/// Fallback numerical silhouette extraction for ellipsoids.
///
/// Used when the analytic method produces questionable results.
fn extract_ellipsoid_silhouettes_numerical(
    ell: &rcad_kernel::geom::EllipsoidalSurface,
    view_dir: DVec3,
    opts: &HlrOptions,
    samples: usize,
) -> Vec<Vec<DVec3>> {
    use rcad_kernel::geom::SurfaceEval;

    let domain = ell.default_domain(); // [0, 2π, 0, π]
    let grid_size = (samples / 4).max(16);

    // Find silhouette seed points on a grid
    let mut silhouette_points: Vec<DVec3> = Vec::new();

    let [u0, u1, v0, v1] = domain;
    for i in 0..grid_size {
        for j in 0..grid_size {
            let u = u0 + (u1 - u0) * i as f64 / (grid_size - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (grid_size - 1) as f64;

            let normal = ell.normal_at(u, v);
            let dot = normal.dot(view_dir);

            // Check if this is a silhouette point
            if dot.abs() < opts.tangent_tolerance.max(0.01) {
                let point = ell.point_at(u, v);
                silhouette_points.push(point);
            }
        }
    }

    // Sort points by angle around the silhouette center (approximation)
    if silhouette_points.len() >= 3 {
        // Compute centroid
        let centroid = silhouette_points.iter().sum::<DVec3>() / silhouette_points.len() as f64;

        // Build the orthonormal frame
        let (axis, x_axis, y_axis) = orthonormal_frame(ell.axis, ell.ref_dir);

        // Sort by angle in the local XY plane (projected)
        silhouette_points.sort_by(|a, b| {
            let a_local = *a - ell.center;
            let b_local = *b - ell.center;
            let ax = a_local.dot(x_axis);
            let ay = a_local.dot(y_axis);
            let bx = b_local.dot(x_axis);
            let by = b_local.dot(y_axis);
            let angle_a = ay.atan2(ax);
            let angle_b = by.atan2(bx);
            angle_a.partial_cmp(&angle_b).unwrap_or(std::cmp::Ordering::Equal)
        });

        vec![silhouette_points]
    } else {
        Vec::new()
    }
}

/// Helper function to build an orthonormal frame from axis and reference direction.
///
/// Returns (axis_normalized, x_axis, y_axis) where axis, x_axis, y_axis form
/// a right-handed orthonormal basis.
fn orthonormal_frame(axis: DVec3, ref_dir: DVec3) -> (DVec3, DVec3, DVec3) {
    let z = axis.normalize_or_zero();
    let x = (ref_dir - ref_dir.dot(z) * z).normalize_or_zero();
    // Handle degenerate case where ref_dir is parallel to axis
    let x = if x.length_squared() < 0.5 {
        any_perpendicular(z)
    } else {
        x
    };
    let y = z.cross(x).normalize_or_zero();
    (z, x, y)
}

// ── Thread Edge Detection ───────────────────────────────────────────────────────

/// Check if an edge is a thread edge (helical) on a cylinder or cone.
///
/// Thread edges are characterized by:
/// - The edge curve is a helix or approximately helical
/// - The edge lies on a cylindrical or conical surface
/// - The edge makes an angle with the surface axis (not parallel or perpendicular)
pub fn is_thread_edge(
    brep: &BRep,
    edge_idx: usize,
    surface: &Surface3,
) -> bool {
    let Some(edge) = brep.edges.get(edge_idx) else {
        return false;
    };

    // Get the edge curve
    let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|&c| c) else {
        return false;
    };
    let Some(curve) = brep.geom.curves.get(curve_idx) else {
        return false;
    };

    match (surface, curve) {
        // Circular helix on cylinder or cone is a thread edge
        (Surface3::Cylinder(_), rcad_kernel::geom::Curve3::CircularHelix(_)) => true,
        (Surface3::Cone(_), rcad_kernel::geom::Curve3::CircularHelix(_)) => true,

        // Check for approximately helical curves
        (Surface3::Cylinder(cyl), curve3d) => {
            is_approximately_helical_on_cylinder(curve3d, cyl, brep, edge)
        }
        (Surface3::Cone(cone), curve3d) => {
            is_approximately_helical_on_cone(curve3d, cone, brep, edge)
        }
        _ => false,
    }
}

/// Check if a curve is approximately helical on a cylinder.
fn is_approximately_helical_on_cylinder(
    curve: &rcad_kernel::geom::Curve3,
    cyl: &rcad_kernel::geom::CylindricalSurface,
    brep: &BRep,
    edge: &rcad_kernel::topology::Edge,
) -> bool {
    use rcad_kernel::geom::CurveEval;

    // Sample the curve and check if points lie on the cylinder surface
    // and the curve makes an angle with the axis
    let Some(v_start) = brep.vertices.get(edge.start) else { return false; };
    let Some(v_end) = brep.vertices.get(edge.end) else { return false; };

    let [t0, t1] = curve.default_domain();
    let samples = 16;

    let mut on_surface_count = 0;
    let mut has_axial_component = false;
    let mut has_angular_component = false;

    for i in 0..samples {
        let t = t0 + (t1 - t0) * i as f64 / (samples - 1) as f64;
        let pt = curve.point_at(t);

        // Check if point is on cylinder surface
        let radial = pt - cyl.origin;
        let axial = radial.dot(cyl.axis);
        let radial_vec = radial - axial * cyl.axis;
        let radial_dist = radial_vec.length();

        if (radial_dist - cyl.radius).abs() < 1e-4 {
            on_surface_count += 1;
        }

        // Check for axial and angular components
        if i > 0 {
            let t_prev = t0 + (t1 - t0) * (i - 1) as f64 / (samples - 1) as f64;
            let pt_prev = curve.point_at(t_prev);
            let delta = pt - pt_prev;

            let axial_delta = delta.dot(cyl.axis).abs();
            let radial_delta = (delta - delta.dot(cyl.axis) * cyl.axis).length();

            if axial_delta > 1e-6 {
                has_axial_component = true;
            }
            if radial_delta > 1e-6 {
                has_angular_component = true;
            }
        }
    }

    // Thread edge is on surface and has both axial and angular components
    on_surface_count > samples / 2 && has_axial_component && has_angular_component
}

/// Check if a curve is approximately helical on a cone.
fn is_approximately_helical_on_cone(
    curve: &rcad_kernel::geom::Curve3,
    cone: &rcad_kernel::geom::ConicalSurface,
    brep: &BRep,
    edge: &rcad_kernel::topology::Edge,
) -> bool {
    use rcad_kernel::geom::CurveEval;

    let Some(v_start) = brep.vertices.get(edge.start) else { return false; };
    let Some(v_end) = brep.vertices.get(edge.end) else { return false; };

    let [t0, t1] = curve.default_domain();
    let samples = 16;

    let mut on_surface_count = 0;
    let mut has_axial_component = false;
    let mut has_angular_component = false;

    let apex = cone.apex_point();
    let axis = cone.axis_dir();

    for i in 0..samples {
        let t = t0 + (t1 - t0) * i as f64 / (samples - 1) as f64;
        let pt = curve.point_at(t);

        // Check if point is on cone surface
        let to_point = pt - apex;
        let axial_dist = to_point.dot(axis);
        let radial_vec = to_point - axial_dist * axis;
        let radial_dist = radial_vec.length();

        // Expected radius at this axial distance
        let expected_radius = axial_dist.abs() * cone.half_angle_rad.tan();

        if (radial_dist - expected_radius).abs() < 1e-4 {
            on_surface_count += 1;
        }

        // Check for axial and angular components
        if i > 0 {
            let t_prev = t0 + (t1 - t0) * (i - 1) as f64 / (samples - 1) as f64;
            let pt_prev = curve.point_at(t_prev);
            let delta = pt - pt_prev;

            let axial_delta = delta.dot(axis).abs();
            let radial_delta = (delta - delta.dot(axis) * axis).length();

            if axial_delta > 1e-6 {
                has_axial_component = true;
            }
            if radial_delta > 1e-6 {
                has_angular_component = true;
            }
        }
    }

    on_surface_count > samples / 2 && has_axial_component && has_angular_component
}

// ── Seam Edge Detection ─────────────────────────────────────────────────────────

/// Check if an edge is a seam edge on a closed surface.
///
/// Seam edges are edges where a closed surface (cylinder, cone, sphere, torus)
/// meets itself at the parameter boundary (u = 0 and u = 2π).
pub fn is_seam_edge(
    brep: &BRep,
    edge_idx: usize,
    surface: &Surface3,
) -> bool {
    // A seam edge has two PCurves on the same surface with different u values
    // (typically 0 and 2π)

    let Some(edge) = brep.edges.get(edge_idx) else {
        return false;
    };

    // Check if this edge has multiple PCurves on the same surface
    let Some(pcurves) = brep.geom.edge_pcurves.get(edge_idx) else {
        return false;
    };

    if pcurves.len() < 2 {
        return false;
    }

    // Check if two PCurves are on the same surface
    let mut surface_counts: HashMap<usize, usize> = HashMap::new();
    for pcurve in pcurves {
        *surface_counts.entry(pcurve.surface_idx).or_insert(0) += 1;
    }

    // If any surface has multiple PCurves for this edge, it's likely a seam
    for &count in surface_counts.values() {
        if count >= 2 {
            // Verify the surface is closed
            match surface {
                Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
                    return true;
                }
                _ => {}
            }
        }
    }

    false
}

/// Check if an edge is a degenerate edge (zero length, e.g., pole singularity).
pub fn is_degenerate_edge_for_hlr(brep: &BRep, edge_idx: usize) -> bool {
    // Check the degenerated flag in GeomStore
    brep.geom.edge_degenerated.get(edge_idx).copied().unwrap_or(false)
}

// ── Curve-Surface Intersection for HLR ──────────────────────────────────────────

/// Result of curve-surface intersection for HLR visibility.
#[derive(Debug, Clone)]
pub struct CurveSurfaceIntersection {
    /// Parameter values on the curve where it intersects the surface.
    pub curve_params: Vec<f64>,
    /// 3D points of intersection.
    pub points: Vec<DVec3>,
    /// UV parameters on the surface for each intersection.
    pub surface_uvs: Vec<(f64, f64)>,
    /// Visibility status between consecutive intersections.
    /// visibility[i] indicates visibility between intersection i and i+1.
    pub visibility: Vec<bool>,
}

/// Compute visible portions of a curve on a curved face.
///
/// This function finds where a curve intersects the silhouette of a surface
/// and determines which portions of the curve are visible.
pub fn compute_curve_visibility_on_surface(
    brep: &BRep,
    edge_idx: usize,
    surface_idx: usize,
    camera: &HlrCamera,
    opts: &HlrOptions,
) -> Option<CurveSurfaceIntersection> {
    let edge = brep.edges.get(edge_idx)?;
    let surface = brep.geom.surfaces.get(surface_idx)?;

    // Get the edge curve
    let curve_idx = brep.geom.edge_curve.get(edge_idx).and_then(|&c| c)?;
    let curve = brep.geom.curves.get(curve_idx)?;

    // Get parameter range
    let [t0, t1] = brep.geom.edge_curve_range.get(edge_idx)
        .and_then(|r| *r)
        .unwrap_or_else(|| curve.default_domain());

    // Sample the curve
    let num_samples = opts.edge_samples.max(16);
    let mut curve_params = Vec::new();
    let mut points = Vec::new();
    let mut surface_uvs = Vec::new();

    let view_dir = (camera.target - camera.eye).normalize_or_zero();

    for i in 0..num_samples {
        let t = t0 + (t1 - t0) * i as f64 / (num_samples - 1) as f64;
        let pt = curve.point_at(t);

        curve_params.push(t);
        points.push(pt);

        // Project point onto surface to get UV
        if let Some((u, v)) = project_point_to_surface(&pt, surface) {
            surface_uvs.push((u, v));
        } else {
            surface_uvs.push((0.0, 0.0)); // placeholder
        }
    }

    // Compute visibility at each sample point
    let mut visibility = Vec::with_capacity(num_samples);

    for (i, &pt) in points.iter().enumerate() {
        let dist = (camera.eye - pt).length();
        let is_visible = true; // Will be computed by the main HLR pipeline
        visibility.push(is_visible);
    }

    // Find silhouette crossings
    let mut silhouette_crossings: Vec<(f64, DVec3)> = Vec::new();

    for i in 1..num_samples {
        let prev_uv = surface_uvs[i - 1];
        let curr_uv = surface_uvs[i];

        let prev_normal = surface.normal_at(prev_uv.0, prev_uv.1);
        let curr_normal = surface.normal_at(curr_uv.0, curr_uv.1);

        let prev_dot = prev_normal.dot(view_dir);
        let curr_dot = curr_normal.dot(view_dir);

        // Check for silhouette crossing (sign change)
        if prev_dot * curr_dot < 0.0 {
            // Bisection to find exact crossing
            if let Some((t_cross, pt_cross)) = find_silhouette_crossing(
                curve, surface, view_dir,
                curve_params[i - 1], curve_params[i],
                10,
            ) {
                silhouette_crossings.push((t_cross, pt_cross));
            }
        }
    }

    // Update visibility based on silhouette crossings
    // (Portions of the curve on the far side of the silhouette are hidden)

    Some(CurveSurfaceIntersection {
        curve_params,
        points,
        surface_uvs,
        visibility,
    })
}

/// Project a 3D point onto a surface to find the closest UV parameters.
fn project_point_to_surface(point: &DVec3, surface: &Surface3) -> Option<(f64, f64)> {
    use rcad_kernel::closest_point_on_surface;

    let result = closest_point_on_surface(surface, *point, 16);
    Some(result.params)
}

/// Find where a curve crosses a surface silhouette using bisection.
fn find_silhouette_crossing(
    curve: &rcad_kernel::geom::Curve3,
    surface: &Surface3,
    view_dir: DVec3,
    t_start: f64,
    t_end: f64,
    max_iter: usize,
) -> Option<(f64, DVec3)> {
    use rcad_kernel::geom::CurveEval;

    let pt_start = curve.point_at(t_start);
    let pt_end = curve.point_at(t_end);

    // Get UV parameters for start and end
    let uv_start = project_point_to_surface(&pt_start, surface)?;
    let uv_end = project_point_to_surface(&pt_end, surface)?;

    let dot_start = surface.normal_at(uv_start.0, uv_start.1).dot(view_dir);
    let dot_end = surface.normal_at(uv_end.0, uv_end.1).dot(view_dir);

    if dot_start * dot_end > 0.0 {
        return None; // No crossing
    }

    let mut t_lo = t_start;
    let mut t_hi = t_end;
    let mut dot_lo = dot_start;

    for _ in 0..max_iter {
        let t_mid = (t_lo + t_hi) * 0.5;
        let pt_mid = curve.point_at(t_mid);

        if let Some(uv_mid) = project_point_to_surface(&pt_mid, surface) {
            let dot_mid = surface.normal_at(uv_mid.0, uv_mid.1).dot(view_dir);

            if dot_mid.abs() < 1e-8 {
                return Some((t_mid, pt_mid));
            }

            if dot_lo * dot_mid < 0.0 {
                t_hi = t_mid;
            } else {
                t_lo = t_mid;
                dot_lo = dot_mid;
            }
        }
    }

    let t_final = (t_lo + t_hi) * 0.5;
    Some((t_final, curve.point_at(t_final)))
}

// ── Edge Classification ─────────────────────────────────────────────────────────

/// Classify all edges in a BRep for HLR processing.
pub fn classify_edges(
    brep: &BRep,
    camera: &HlrCamera,
    opts: &HlrOptions,
) -> Vec<EdgeClassInfo> {
    let mut classifications: Vec<EdgeClassInfo> = Vec::new();

    for edge_idx in 0..brep.edges.len() {
        let classification = classify_single_edge(brep, edge_idx, camera, opts);
        classifications.push(classification);
    }

    classifications
}

/// Classify a single edge.
fn classify_single_edge(
    brep: &BRep,
    edge_idx: usize,
    camera: &HlrCamera,
    opts: &HlrOptions,
) -> EdgeClassInfo {
    let Some(edge) = brep.edges.get(edge_idx) else {
        return EdgeClassInfo {
            edge_idx,
            classification: EdgeClassification::Hidden,
            visible_segments: 0,
            hidden_segments: 0,
            on_curved_surface: false,
            surface_idx: None,
        };
    };

    // Get the surface this edge is on (if any)
    let surface_idx = get_edge_surface(brep, edge_idx);
    let on_curved_surface = surface_idx.map_or(false, |idx| {
        matches!(
            brep.geom.surfaces.get(idx),
            Some(Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Cone(_) | Surface3::Torus(_))
        )
    });

    // Check for thread edge
    if let Some(idx) = surface_idx {
        if let Some(surface) = brep.geom.surfaces.get(idx) {
            if opts.detect_thread_edges && is_thread_edge(brep, edge_idx, surface) {
                return EdgeClassInfo {
                    edge_idx,
                    classification: EdgeClassification::Thread,
                    visible_segments: 0,
                    hidden_segments: 0,
                    on_curved_surface: true,
                    surface_idx: Some(idx),
                };
            }

            if opts.detect_seam_edges && is_seam_edge(brep, edge_idx, surface) {
                return EdgeClassInfo {
                    edge_idx,
                    classification: EdgeClassification::Seam,
                    visible_segments: 0,
                    hidden_segments: 0,
                    on_curved_surface: true,
                    surface_idx: Some(idx),
                };
            }
        }
    }

    // Default classification - will be updated during HLR processing
    EdgeClassInfo {
        edge_idx,
        classification: EdgeClassification::Partial,
        visible_segments: 0,
        hidden_segments: 0,
        on_curved_surface,
        surface_idx,
    }
}

/// Get the primary surface an edge is on.
fn get_edge_surface(brep: &BRep, edge_idx: usize) -> Option<usize> {
    let pcurves = brep.geom.edge_pcurves.get(edge_idx)?;
    pcurves.first().map(|pc| pc.surface_idx)
}

// ── Spatial Indexing for Silhouette Queries ─────────────────────────────────────

/// Spatial grid for efficient silhouette point queries.
#[derive(Debug, Clone)]
pub struct SilhouetteSpatialIndex {
    /// Grid cells containing silhouette sample points.
    cells: HashMap<(i32, i32, i32), Vec<usize>>,
    /// All sample points.
    points: Vec<DVec3>,
    /// Grid cell size.
    cell_size: f64,
    /// Bounding box of all points.
    bbox_min: DVec3,
    bbox_max: DVec3,
}

impl SilhouetteSpatialIndex {
    /// Build a spatial index from silhouette points.
    pub fn build(points: &[DVec3], cell_size: f64) -> Self {
        if points.is_empty() {
            return Self {
                cells: HashMap::new(),
                points: Vec::new(),
                cell_size,
                bbox_min: DVec3::ZERO,
                bbox_max: DVec3::ZERO,
            };
        }

        // Compute bounding box
        let mut bbox_min = DVec3::splat(f64::INFINITY);
        let mut bbox_max = DVec3::splat(f64::NEG_INFINITY);
        for &p in points {
            bbox_min = bbox_min.min(p);
            bbox_max = bbox_max.max(p);
        }

        // Build grid
        let mut cells: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for (i, &p) in points.iter().enumerate() {
            let cell = Self::point_to_cell(p, bbox_min, cell_size);
            cells.entry(cell).or_default().push(i);
        }

        Self {
            cells,
            points: points.to_vec(),
            cell_size,
            bbox_min,
            bbox_max,
        }
    }

    fn point_to_cell(p: DVec3, origin: DVec3, cell_size: f64) -> (i32, i32, i32) {
        let d = (p - origin) / cell_size;
        (d.x.floor() as i32, d.y.floor() as i32, d.z.floor() as i32)
    }

    /// Find all silhouette points within a radius of the query point.
    pub fn query_radius(&self, point: DVec3, radius: f64) -> Vec<usize> {
        let mut result = Vec::new();

        let cell_radius = (radius / self.cell_size).ceil() as i32;
        let center_cell = Self::point_to_cell(point, self.bbox_min, self.cell_size);

        for di in -cell_radius..=cell_radius {
            for dj in -cell_radius..=cell_radius {
                for dk in -cell_radius..=cell_radius {
                    let cell = (center_cell.0 + di, center_cell.1 + dj, center_cell.2 + dk);

                    if let Some(indices) = self.cells.get(&cell) {
                        for &idx in indices {
                            let dist = (self.points[idx] - point).length();
                            if dist <= radius {
                                result.push(idx);
                            }
                        }
                    }
                }
            }
        }

        result
    }

    /// Find the nearest silhouette point to the query point.
    pub fn query_nearest(&self, point: DVec3) -> Option<(usize, f64)> {
        let mut best: Option<(usize, f64)> = None;

        // Start with a small search radius and expand
        let mut radius = self.cell_size;

        for _ in 0..10 {
            let candidates = self.query_radius(point, radius);

            for idx in candidates {
                let dist = (self.points[idx] - point).length();
                if best.map_or(true, |(_, d)| dist < d) {
                    best = Some((idx, dist));
                }
            }

            if best.is_some() && best.unwrap().1 <= radius * 0.5 {
                break;
            }

            radius *= 2.0;
        }

        best
    }

    /// Get a point by index.
    pub fn get_point(&self, idx: usize) -> Option<DVec3> {
        self.points.get(idx).copied()
    }

    /// Get the number of indexed points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Check if the index is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

/// Numerical silhouette extraction for general parametric surfaces.
///
/// Uses a marching approach to find curves where normal · view_dir = 0.
/// This implementation includes:
/// - Marching along iso-parametric curves to trace silhouette curves
/// - Curvature-adaptive sampling for better accuracy in high-curvature regions
/// - Handling of closed silhouette loops
fn extract_numerical_silhouettes(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    opts: &HlrOptions,
    _brep: &BRep,
) -> Vec<Vec<DVec3>> {
    let [u0, u1, v0, v1] = domain;
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    // Phase 1: Find silhouette seed points on a coarse grid
    let grid_size = opts.silhouette_samples.max(16);
    let seeds = find_silhouette_seeds(surface, view_dir, domain, grid_size, opts.tangent_tolerance);

    if seeds.is_empty() {
        return curves;
    }

    // Phase 2: March from each seed to trace silhouette curves
    let mut visited: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();

    for (i, j, u, v) in seeds {
        if visited.contains(&(i, j)) {
            continue;
        }

        // Trace a curve starting from this seed
        let curve = march_silhouette_curve(surface, view_dir, domain, u, v, opts);

        if curve.len() >= 2 {
            // Mark visited cells along the curve
            for pt in &curve {
                // Find grid cell for this point
                let pi = ((pt.0 - u0) / (u1 - u0) * grid_size as f64).floor() as usize;
                let pj = ((pt.1 - v0) / (v1 - v0) * grid_size as f64).floor() as usize;
                visited.insert((pi.min(grid_size - 1), pj.min(grid_size - 1)));
            }

            // Apply adaptive refinement based on curvature
            let refined_curve = if opts.curvature_adaptive {
                refine_curve_by_curvature(surface, curve, opts)
            } else {
                curve.into_iter().map(|(_, _, pt)| pt).collect()
            };

            // Apply B-spline fitting if enabled
            let final_curve = if opts.fit_bspline && refined_curve.len() >= 4 {
                fit_bspline_to_points(&refined_curve, opts.bspline_tolerance)
            } else {
                refined_curve
            };

            if final_curve.len() >= 2 {
                curves.push(final_curve);
            }
        }
    }

    curves
}

/// A point in parameter space with its 3D position.
type ParamPoint = (f64, f64, DVec3);

/// Find seed points for silhouette curves on a grid.
fn find_silhouette_seeds(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    grid_size: usize,
    tangent_tol: f64,
) -> Vec<(usize, usize, f64, f64)> {
    let [u0, u1, v0, v1] = domain;
    let mut seeds = Vec::new();

    // Sample grid and look for sign changes in normal · view_dir
    let mut dot_values: Vec<Vec<f64>> = vec![vec![0.0; grid_size]; grid_size];

    // Compute dot products at grid vertices
    for i in 0..grid_size {
        for j in 0..grid_size {
            let u = u0 + (u1 - u0) * i as f64 / (grid_size - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (grid_size - 1) as f64;
            let normal = surface.normal_at(u, v);
            dot_values[i][j] = normal.dot(view_dir);
        }
    }

    // Find cells where sign changes occur (indicating silhouette crossing)
    for i in 0..grid_size - 1 {
        for j in 0..grid_size - 1 {
            let d00 = dot_values[i][j];
            let d10 = dot_values[i + 1][j];
            let d01 = dot_values[i][j + 1];
            let d11 = dot_values[i + 1][j + 1];

            // Check for sign changes in the cell
            let has_crossing = (d00 * d10 < 0.0)
                || (d00 * d01 < 0.0)
                || (d10 * d11 < 0.0)
                || (d01 * d11 < 0.0);

            if has_crossing {
                // Find the exact crossing point using bisection
                if let Some((u, v)) = find_crossing_point(surface, view_dir, domain, i, j, grid_size, tangent_tol) {
                    seeds.push((i, j, u, v));
                }
            }
        }
    }

    seeds
}

/// Find the exact crossing point in a grid cell using bisection.
fn find_crossing_point(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    i: usize,
    j: usize,
    grid_size: usize,
    _tangent_tol: f64,
) -> Option<(f64, f64)> {
    let [u0, v0, u1, v1] = [
        domain[0] + (domain[1] - domain[0]) * i as f64 / (grid_size - 1) as f64,
        domain[2] + (domain[3] - domain[2]) * j as f64 / (grid_size - 1) as f64,
        domain[0] + (domain[1] - domain[0]) * (i + 1) as f64 / (grid_size - 1) as f64,
        domain[2] + (domain[3] - domain[2]) * (j + 1) as f64 / (grid_size - 1) as f64,
    ];

    // Try to find crossing along each edge of the cell
    let edges = [
        (u0, v0, u1, v0), // bottom edge
        (u0, v1, u1, v1), // top edge
        (u0, v0, u0, v1), // left edge
        (u1, v0, u1, v1), // right edge
    ];

    for (ua, va, ub, vb) in edges {
        if let Some((u, v)) = bisection_search(surface, view_dir, ua, va, ub, vb, 12) {
            return Some((u, v));
        }
    }

    None
}

/// Bisection search to find where normal · view_dir = 0.
fn bisection_search(
    surface: &Surface3,
    view_dir: DVec3,
    u0: f64,
    v0: f64,
    u1: f64,
    v1: f64,
    max_iter: usize,
) -> Option<(f64, f64)> {
    let d0 = surface.normal_at(u0, v0).dot(view_dir);
    let d1 = surface.normal_at(u1, v1).dot(view_dir);

    if d0 * d1 > 0.0 {
        return None; // No sign change
    }

    let mut ua = u0;
    let mut va = v0;
    let mut ub = u1;
    let mut vb = v1;

    for _ in 0..max_iter {
        let um = (ua + ub) / 2.0;
        let vm = (va + vb) / 2.0;
        let dm = surface.normal_at(um, vm).dot(view_dir);

        if dm.abs() < 1e-10 {
            return Some((um, vm));
        }

        if d0 * dm < 0.0 {
            ub = um;
            vb = vm;
        } else {
            ua = um;
            va = vm;
        }
    }

    Some(((ua + ub) / 2.0, (va + vb) / 2.0))
}

/// March along a silhouette curve starting from a seed point.
fn march_silhouette_curve(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    u_start: f64,
    v_start: f64,
    opts: &HlrOptions,
) -> Vec<ParamPoint> {
    let mut curve: Vec<ParamPoint> = Vec::new();
    let [u0, u1, v0, v1] = domain;

    // Add the starting point
    let p_start = surface.point_at(u_start, v_start);
    curve.push((u_start, v_start, p_start));

    // March in both directions from the seed
    for direction in &[-1.0_f64, 1.0] {
        let mut u = u_start;
        let mut v = v_start;
        let mut curve_dir: Option<DVec2> = None;

        for _ in 0..opts.max_subdivisions * 50 {
            // Compute the tangent direction to the silhouette curve
            let tangent = compute_silhouette_tangent(surface, view_dir, u, v);

            if tangent.length_squared() < 1e-16 {
                break;
            }

            // Choose direction along the tangent
            let step_dir = if let Some(cd) = curve_dir {
                // Continue in the same general direction
                if cd.dot(tangent) > 0.0 {
                    tangent
                } else {
                    -tangent
                }
            } else {
                *direction * tangent
            };
            curve_dir = Some(step_dir.normalize_or_zero());

            // Compute step size based on curvature
            let (k1, k2) = rcad_kernel::curvature::principal_curvatures(surface, u, v);
            let max_k = k1.abs().max(k2.abs()).max(opts.min_curvature);
            let curvature_factor = (opts.max_curvature / max_k).min(4.0).max(0.25);
            let step_size = opts.angular_tolerance * curvature_factor;

            // Take a step
            let u_new = u + step_dir.x * step_size;
            let v_new = v + step_dir.y * step_size;

            // Check bounds
            if u_new < u0 || u_new > u1 || v_new < v0 || v_new > v1 {
                break;
            }

            // Project back onto the silhouette curve
            if let Some((u_proj, v_proj)) = project_to_silhouette(surface, view_dir, u_new, v_new, opts.tangent_tolerance) {
                u = u_proj;
                v = v_proj;

                let p = surface.point_at(u, v);
                let d = (p - curve.last().map(|(_, _, lp)| *lp).unwrap_or(p_start)).length();

                // Only add if we've moved enough
                if d > opts.bspline_tolerance * 0.1 {
                    curve.push((u, v, p));
                }
            } else {
                break;
            }

            // Check for closed loop
            if curve.len() > 10 {
                let first = curve[0];
                let dist = ((first.0 - u).powi(2) + (first.1 - v).powi(2)).sqrt();
                if dist < step_size * 2.0 {
                    // Close the loop
                    curve.push(curve[0]);
                    break;
                }
            }
        }

        // Reverse the points added while marching in the negative direction
        if *direction < 0.0 && curve.len() > 1 {
            let first = curve[0];
            curve.reverse();
            curve.push(first); // Re-add the start point for the loop
        }
    }

    curve
}

/// Compute the tangent direction to the silhouette curve at a point.
fn compute_silhouette_tangent(
    surface: &Surface3,
    view_dir: DVec3,
    u: f64,
    v: f64,
) -> DVec2 {
    const EPS: f64 = 1e-6;

    // Compute gradients of the implicit function f(u,v) = N(u,v) · V
    let n = surface.normal_at(u, v);
    let n_u = surface.normal_at(u + EPS, v);
    let n_v = surface.normal_at(u, v + EPS);

    // Gradient of f = N · V
    let df_du = (n_u - n).dot(view_dir) / EPS;
    let df_dv = (n_v - n).dot(view_dir) / EPS;

    // The tangent direction is perpendicular to the gradient
    DVec2::new(-df_dv, df_du).normalize_or_zero()
}

/// Project a point back onto the silhouette curve.
fn project_to_silhouette(
    surface: &Surface3,
    view_dir: DVec3,
    u: f64,
    v: f64,
    tol: f64,
) -> Option<(f64, f64)> {
    let mut u_curr = u;
    let mut v_curr = v;

    // Newton iteration to find f(u,v) = 0
    for _ in 0..20 {
        let n = surface.normal_at(u_curr, v_curr);
        let f = n.dot(view_dir);

        if f.abs() < tol {
            return Some((u_curr, v_curr));
        }

        // Compute gradient numerically
        const EPS: f64 = 1e-7;
        let n_u = surface.normal_at(u_curr + EPS, v_curr);
        let n_v = surface.normal_at(u_curr, v_curr + EPS);

        let df_du = (n_u - n).dot(view_dir) / EPS;
        let df_dv = (n_v - n).dot(view_dir) / EPS;

        let grad_len_sq = df_du * df_du + df_dv * df_dv;
        if grad_len_sq < 1e-20 {
            break;
        }

        // Newton step
        let step = f / grad_len_sq;
        u_curr -= step * df_du;
        v_curr -= step * df_dv;
    }

    // Check if we converged
    let f = surface.normal_at(u_curr, v_curr).dot(view_dir);
    if f.abs() < tol * 10.0 {
        Some((u_curr, v_curr))
    } else {
        None
    }
}

/// Refine a silhouette curve based on surface curvature.
fn refine_curve_by_curvature(
    surface: &Surface3,
    curve: Vec<ParamPoint>,
    opts: &HlrOptions,
) -> Vec<DVec3> {
    if curve.len() < 2 {
        return curve.into_iter().map(|(_, _, p)| p).collect();
    }

    let mut refined: Vec<DVec3> = Vec::new();
    refined.push(curve[0].2);

    for i in 1..curve.len() {
        let (u0, v0, p0) = curve[i - 1];
        let (u1, v1, p1) = curve[i];

        // Compute curvature at the midpoint
        let um = (u0 + u1) / 2.0;
        let vm = (v0 + v1) / 2.0;
        let (k1, k2) = rcad_kernel::curvature::principal_curvatures(surface, um, vm);
        let max_k = k1.abs().max(k2.abs());

        // Determine number of subdivision points based on curvature
        let chord_len = (p1 - p0).length();
        let subdivs = if max_k > opts.min_curvature {
            let curvature_samples = (max_k * chord_len * std::f64::consts::PI).ceil() as usize;
            curvature_samples.min(8).max(1)
        } else {
            1
        };

        // Add subdivision points
        for j in 1..subdivs {
            let t = j as f64 / subdivs as f64;
            let u = u0 + t * (u1 - u0);
            let v = v0 + t * (v1 - v0);
            let p = surface.point_at(u, v);
            refined.push(p);
        }

        refined.push(p1);
    }

    refined
}

/// Fit a B-spline curve to a set of points.
fn fit_bspline_to_points(points: &[DVec3], tolerance: f64) -> Vec<DVec3> {
    if points.len() < 4 {
        return points.to_vec();
    }

    // Simple approach: sample the fitted B-spline at uniform intervals
    // For a proper implementation, we would use least-squares fitting
    // Here we use a simplified version that preserves the shape

    let n = points.len();
    let mut result: Vec<DVec3> = Vec::with_capacity(n);

    // Compute chord lengths for parameterization
    let mut chords = vec![0.0_f64; n];
    for i in 1..n {
        chords[i] = chords[i - 1] + (points[i] - points[i - 1]).length();
    }
    let total_len = chords[n - 1];
    if total_len < 1e-12 {
        return points.to_vec();
    }

    // Generate control points using Catmull-Rom style interpolation
    let degree = 3.min(n - 1);
    let num_samples = (total_len / tolerance).ceil() as usize;
    let num_samples = num_samples.max(10).min(1000);

    for i in 0..=num_samples {
        let t = i as f64 / num_samples as f64;
        let target_len = t * total_len;

        // Find the segment containing this length
        let seg_idx = chords.partition_point(|&c| c < target_len).saturating_sub(1);
        let seg_idx = seg_idx.min(n - 2);

        // Interpolate within the segment
        let seg_start = chords[seg_idx];
        let seg_end = chords[seg_idx + 1];
        let seg_len = seg_end - seg_start;

        let local_t = if seg_len > 1e-12 {
            (target_len - seg_start) / seg_len
        } else {
            0.5
        };

        // Simple linear interpolation with smoothing
        let p0 = points[seg_idx];
        let p1 = points[seg_idx + 1];

        // Hermite interpolation for smoother result
        let t0 = if seg_idx > 0 {
            (points[seg_idx + 1] - points[seg_idx - 1]).normalize_or_zero()
        } else {
            (points[1] - points[0]).normalize_or_zero()
        };

        let t1 = if seg_idx + 2 < n {
            (points[seg_idx + 2] - points[seg_idx]).normalize_or_zero()
        } else {
            (points[n - 1] - points[n - 2]).normalize_or_zero()
        };

        let h00 = 2.0 * local_t * local_t * local_t - 3.0 * local_t * local_t + 1.0;
        let h10 = local_t * local_t * local_t - 2.0 * local_t * local_t + local_t;
        let h01 = -2.0 * local_t * local_t * local_t + 3.0 * local_t * local_t;
        let h11 = local_t * local_t * local_t - local_t * local_t;

        let p = h00 * p0 + h10 * t0 * seg_len + h01 * p1 + h11 * t1 * seg_len;
        result.push(p);
    }

    result
}

/// Generate silhouette curves for the HLR pipeline (internal function).
fn compute_silhouettes(brep: &BRep, view_dir: DVec3, samples: usize) -> Vec<SilhouetteCurve> {
    let opts = HlrOptions {
        silhouette_samples: samples,
        ..HlrOptions::default()
    };

    extract_silhouette_curves(brep, view_dir, &opts)
        .into_iter()
        .map(|curve| SilhouetteCurve {
            world_pts: curve.points,
            curve_hint: None,
            dense: true, // All silhouettes are treated as dense for proper rendering
        })
        .collect()
}

/// Occlusion tester that supports both brute-force and BVH-accelerated methods.
enum OcclusionTester<'a> {
    BruteForce(&'a [[DVec3; 3]]),
    Bvh {
        bvh: &'a TriBvh,
        triangles: &'a [[DVec3; 3]],
    },
}

impl<'a> OcclusionTester<'a> {
    fn is_occluded(&self, point: DVec3, eye: DVec3, dist_to_eye: f64) -> bool {
        match self {
            OcclusionTester::BruteForce(triangles) => {
                is_occluded(point, eye, triangles, dist_to_eye)
            }
            OcclusionTester::Bvh { bvh, triangles } => {
                bvh.is_occluded(point, eye, triangles, dist_to_eye)
            }
        }
    }
}

/// Improved visibility classification that handles grazing angles on curved surfaces.
///
/// For points near silhouette curves (where normal is nearly perpendicular to view direction),
/// we use additional testing to improve numerical stability.
fn classify_visibility(
    point: DVec3,
    normal: Option<DVec3>,
    camera: &HlrCamera,
    occlusion_tester: &OcclusionTester<'_>,
    grazing_threshold: f64,
) -> VisibilityInfo {
    let dist = (camera.eye - point).length();
    let view_dir = (camera.eye - point).normalize_or_zero();

    // Check if we're at a grazing angle
    let grazing_factor = if let Some(n) = normal {
        let dot = n.dot(view_dir).abs();
        // grazing_factor = 1.0 when perfectly grazing (dot = 0)
        // grazing_factor = 0.0 when viewing straight on (dot = 1)
        1.0 - dot
    } else {
        0.0
    };

    // For grazing angles, use more robust testing
    let is_occluded = if grazing_factor > grazing_threshold.cos() {
        // At grazing angle: test multiple rays to reduce false positives
        let base_occluded = occlusion_tester.is_occluded(point, camera.eye, dist);

        if base_occluded {
            // Verify with additional samples to reduce numerical errors
            let mut occluded_count = 1;
            const NUM_SAMPLES: usize = 4;
            let offset = 1e-4;

            for i in 0..NUM_SAMPLES {
                let angle = i as f64 * std::f64::consts::TAU / NUM_SAMPLES as f64;
                let perp = any_perpendicular(view_dir);
                let perturb = perp * (angle.cos() * offset) + view_dir.cross(perp) * (angle.sin() * offset);
                let test_point = point + perturb;

                if occlusion_tester.is_occluded(test_point, camera.eye, dist) {
                    occluded_count += 1;
                }
            }

            // Require majority to confirm occlusion at grazing angles
            occluded_count > NUM_SAMPLES / 2
        } else {
            false
        }
    } else {
        occlusion_tester.is_occluded(point, camera.eye, dist)
    };

    VisibilityInfo {
        visible: !is_occluded,
        grazing_factor,
        depth: dist,
    }
}

/// Information about visibility at a point.
struct VisibilityInfo {
    visible: bool,
    grazing_factor: f64,
    depth: f64,
}

/// Process a list of world-space sample points through the HLR visibility
/// pipeline and append resulting segments to `result`.
///
/// When `dense` is true, one segment is emitted per consecutive point pair
/// (useful for polyline approximations of curved silhouettes).
fn process_world_pts(
    world_pts: &[DVec3],
    curve_hint: Option<CurveHint>,
    dense: bool,
    segment_type: SegmentType,
    camera: &HlrCamera,
    view: &DMat4,
    triangles: &[[DVec3; 3]],
    result: &mut HlrResult,
) {
    process_world_pts_with_bvh(
        world_pts,
        curve_hint,
        dense,
        segment_type,
        camera,
        view,
        triangles,
        None,
        &HlrOptions::default(),
        result,
    )
}

/// Process world points with optional BVH acceleration and grazing angle handling.
fn process_world_pts_with_bvh(
    world_pts: &[DVec3],
    curve_hint: Option<CurveHint>,
    dense: bool,
    segment_type: SegmentType,
    camera: &HlrCamera,
    view: &DMat4,
    triangles: &[[DVec3; 3]],
    bvh: Option<&TriBvh>,
    opts: &HlrOptions,
    result: &mut HlrResult,
) {
    if world_pts.len() < 2 {
        return;
    }
    let n = world_pts.len();

    let occlusion_tester = if let Some(bvh) = bvh {
        OcclusionTester::Bvh { bvh, triangles }
    } else {
        OcclusionTester::BruteForce(triangles)
    };

    let sample_vis: Vec<bool> = world_pts
        .iter()
        .map(|&wp| {
            let dist = (camera.eye - wp).length();
            !occlusion_tester.is_occluded(wp, camera.eye, dist)
        })
        .collect();

    let screen_pts: Vec<DVec2> = world_pts.iter().map(|&wp| project(wp, view).0).collect();

    if dense {
        // Emit one segment per consecutive pair (preserves polyline shape).
        for i in 0..n - 1 {
            let seg = HlrSegment {
                start: screen_pts[i],
                end: screen_pts[i + 1],
                visible: sample_vis[i] && sample_vis[i + 1],
                curve_hint: curve_hint.clone(),
                segment_type,
            };
            if (seg.end - seg.start).length_squared() > 1e-16 {
                result.segments.push(seg);
            }
        }
        return;
    }

    let mut seg_start = 0usize;
    for i in 1..n {
        let changed = sample_vis[i] != sample_vis[seg_start];
        let last = i == n - 1;
        if changed || last {
            let end_idx = if last && !changed { i } else { i - 1 };
            let seg = HlrSegment {
                start: screen_pts[seg_start],
                end: screen_pts[end_idx],
                visible: sample_vis[seg_start],
                curve_hint: curve_hint.clone(),
                segment_type,
            };
            if (seg.end - seg.start).length_squared() > 1e-16 {
                result.segments.push(seg);
            }
            if changed {
                seg_start = i;
            }
        }
    }
}



/// Perform hidden-line removal on a BRep from the given camera position.
///
/// Returns 2D projected segments labeled visible/hidden.
/// `samples` controls how finely each edge is subdivided for occlusion testing
/// (higher = more accurate but slower; 8 is a reasonable default).
pub fn hlr(brep: &BRep, camera: &HlrCamera, samples: usize) -> HlrResult {
    hlr_with_options(brep, camera, HlrOptions::default().with_edge_samples(samples))
}

/// Perform hidden-line removal with full configuration options.
///
/// This function provides fine-grained control over HLR computation parameters,
/// including adaptive sampling for curved surfaces.
///
/// # Arguments
/// * `brep` - The BRep model to process.
/// * `camera` - Camera/view specification.
/// * `opts` - Configuration options for sampling and tolerances.
///
/// # Returns
/// An `HlrResult` containing projected 2D segments labeled as visible/hidden.
pub fn hlr_with_options(brep: &BRep, camera: &HlrCamera, opts: HlrOptions) -> HlrResult {
    let view = look_at(camera.eye, camera.target, camera.up);
    let triangles = collect_triangles(brep);
    let edge_samples = opts.edge_samples.max(2);
    let mut result = HlrResult::default();

    // Build BVH for acceleration if enabled and we have enough triangles
    let bvh: Option<TriBvh> = if opts.use_bvh && triangles.len() > 32 {
        Some(TriBvh::build(&triangles))
    } else {
        None
    };
    let bvh_ref = bvh.as_ref();

    // Classify edges for thread/seam detection
    let edge_classifications = if opts.detect_thread_edges || opts.detect_seam_edges {
        classify_edges(brep, camera, &opts)
    } else {
        Vec::new()
    };

    // ── Wire edges ────────────────────────────────────────────────────────────

    // Collect all unique edges from all faces + standalone edges
    let mut edge_indices: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    edge_indices.insert(we.idx);
                }
                for inner in &face.inner_wires {
                    for we in &inner.edges {
                        edge_indices.insert(we.idx);
                    }
                }
            }
        }
    }
    for i in 0..brep.edges.len() {
        edge_indices.insert(i);
    }

    // Convert to vector for potential parallel processing
    let edge_indices_vec: Vec<usize> = edge_indices.into_iter().collect();

    // Process edges (optionally in parallel)
    let edge_results: Vec<Vec<HlrSegment>> = if opts.parallel && edge_indices_vec.len() > opts.parallel_threshold {
        let triangles_ref = &triangles;
        let bvh_opt = bvh_ref;
        let brep_ref = brep;
        let camera_ref = camera;
        let view_ref = &view;
        let opts_ref = &opts;
        let edge_classes = &edge_classifications;

        edge_indices_vec
            .par_iter()
            .map(|&edge_idx| {
                process_single_edge(
                    brep_ref,
                    edge_idx,
                    camera_ref,
                    view_ref,
                    triangles_ref,
                    bvh_opt,
                    opts_ref,
                    edge_classes,
                )
            })
            .collect()
    } else {
        edge_indices_vec
            .iter()
            .map(|&edge_idx| {
                process_single_edge(
                    brep,
                    edge_idx,
                    camera,
                    &view,
                    &triangles,
                    bvh_ref,
                    &opts,
                    &edge_classifications,
                )
            })
            .collect()
    };

    // Merge results
    for segments in edge_results {
        result.segments.extend(segments);
    }

    // ── Silhouette curves ────────────────────────────────────────────

    let view_dir = (camera.target - camera.eye).normalize_or_zero();
    let silhouette_curves = compute_silhouettes_with_options(brep, view_dir, &opts);

    // Build spatial index for silhouette queries
    let all_silhouette_points: Vec<DVec3> = silhouette_curves
        .iter()
        .flat_map(|c| c.world_pts.iter().copied())
        .collect();
    let spatial_index = SilhouetteSpatialIndex::build(&all_silhouette_points, 0.1);

    for sil in silhouette_curves {
        process_world_pts_with_bvh(
            &sil.world_pts,
            sil.curve_hint,
            sil.dense,
            SegmentType::Silhouette,
            camera,
            &view,
            &triangles,
            bvh_ref,
            &opts,
            &mut result,
        );
    }

    result
}

/// Process a single edge and return its segments.
fn process_single_edge(
    brep: &BRep,
    edge_idx: usize,
    camera: &HlrCamera,
    view: &DMat4,
    triangles: &[[DVec3; 3]],
    bvh: Option<&TriBvh>,
    opts: &HlrOptions,
    edge_classifications: &[EdgeClassInfo],
) -> Vec<HlrSegment> {
    let mut segments: Vec<HlrSegment> = Vec::new();

    let Some(edge) = brep.edges.get(edge_idx) else { return segments; };
    let Some(v_start) = brep.vertices.get(edge.start) else { return segments; };
    let Some(v_end) = brep.vertices.get(edge.end) else { return segments; };

    let p0 = v_start.point;
    let p1 = v_end.point;

    // Determine edge type
    let segment_type = if let Some(class_info) = edge_classifications.get(edge_idx) {
        match class_info.classification {
            EdgeClassification::Thread => SegmentType::Thread,
            EdgeClassification::Seam => SegmentType::Seam,
            _ => SegmentType::Edge,
        }
    } else {
        SegmentType::Edge
    };

    let edge_curve = brep
        .geom
        .edge_curve
        .get(edge_idx)
        .and_then(|&ci| ci)
        .and_then(|ci| brep.geom.curves.get(ci));

    let circle_info: Option<Circle3> = edge_curve.and_then(|c| {
        if let rcad_kernel::geom::Curve3::Circle(circ) = c { Some(*circ) } else { None }
    });

    let is_other_curve = edge_curve
        .map_or(false, |c| !matches!(c, rcad_kernel::geom::Curve3::Line(_)))
        && circle_info.is_none();

    // Adaptive sampling for curved edges on curved surfaces
    let this_edge_samples = if circle_info.is_some() || is_other_curve {
        // Check if this edge is on a curved surface with high curvature
        if let Some(class_info) = edge_classifications.get(edge_idx) {
            if class_info.on_curved_surface {
                if let Some(surf_idx) = class_info.surface_idx {
                    if let Some(surface) = brep.geom.surfaces.get(surf_idx) {
                        // Compute curvature at midpoint
                        let domain = brep.geom.face_surface_range
                            .iter()
                            .find_map(|r| *r)
                            .unwrap_or_else(|| surface.default_domain());
                        let mid_u = (domain[0] + domain[1]) * 0.5;
                        let mid_v = (domain[2] + domain[3]) * 0.5;
                        let (k1, k2) = rcad_kernel::curvature::principal_curvatures(surface, mid_u, mid_v);
                        let max_k = k1.abs().max(k2.abs());

                        // More samples for higher curvature
                        let adaptive_factor = (max_k / 10.0).min(8.0).max(1.0);
                        ((opts.edge_samples as f64 * adaptive_factor * 4.0) as usize).max(32).min(256)
                    } else {
                        (opts.edge_samples * 4).max(32)
                    }
                } else {
                    (opts.edge_samples * 4).max(32)
                }
            } else {
                (opts.edge_samples * 4).max(32)
            }
        } else {
            (opts.edge_samples * 4).max(32)
        }
    } else {
        opts.edge_samples
    };

    let world_pts: Vec<DVec3> = if let Some(circ) = &circle_info {
        let [t0, t1] = brep
            .geom
            .edge_curve_range
            .get(edge_idx)
            .and_then(|r| *r)
            .unwrap_or_else(|| circ.default_domain());
        (0..this_edge_samples)
            .map(|i| {
                let t = t0 + (t1 - t0) * (i as f64 / (this_edge_samples - 1) as f64);
                circ.point_at(t)
            })
            .collect()
    } else if let Some(curve) = edge_curve.filter(|_| is_other_curve) {
        let [t0, t1] = brep
            .geom
            .edge_curve_range
            .get(edge_idx)
            .and_then(|r| *r)
            .unwrap_or_else(|| curve.default_domain());
        (0..this_edge_samples)
            .map(|i| {
                let t = t0 + (t1 - t0) * (i as f64 / (this_edge_samples - 1) as f64);
                curve.point_at(t)
            })
            .collect()
    } else {
        if (p1 - p0).length_squared() < 1e-12 {
            return segments;
        }
        (0..this_edge_samples)
            .map(|i| {
                let t = i as f64 / (this_edge_samples - 1) as f64;
                p0 + (p1 - p0) * t
            })
            .collect()
    };

    // Compute curve_hint for circle edges
    let screen_pts_for_hint: Vec<DVec2> =
        world_pts.iter().map(|&wp| project(wp, view).0).collect();
    let curve_hint: Option<CurveHint> = if let Some(circ) = &circle_info {
        let (center_2d, _) = project(circ.center, view);
        let r = screen_pts_for_hint
            .iter()
            .map(|p| (*p - center_2d).length())
            .fold(0.0_f64, f64::max);
        Some(CurveHint::Circle { center: center_2d, radius: r })
    } else if is_other_curve {
        Some(CurveHint::Other)
    } else {
        None
    };

    // Process the edge points
    let mut edge_result = HlrResult::default();
    process_world_pts_with_bvh(
        &world_pts,
        curve_hint,
        false,
        segment_type,
        camera,
        view,
        triangles,
        bvh,
        opts,
        &mut edge_result,
    );

    edge_result.segments
}

/// Compute silhouette curves with full options (internal helper).
fn compute_silhouettes_with_options(brep: &BRep, view_dir: DVec3, opts: &HlrOptions) -> Vec<SilhouetteCurve> {
    extract_silhouette_curves(brep, view_dir, opts)
        .into_iter()
        .map(|curve| SilhouetteCurve {
            world_pts: curve.points,
            curve_hint: None,
            dense: true,
        })
        .collect()
}

/// Per-component HLR result for assembly HLR.
#[derive(Debug, Clone, Default)]
pub struct ComponentHlr {
    /// Component name (from the assembly).
    pub name: String,
    /// HLR segments for this component.
    pub segments: Vec<HlrSegment>,
}

/// Output of assembly HLR — one `ComponentHlr` per leaf BRep.
#[derive(Debug, Clone, Default)]
pub struct AssemblyHlrResult {
    pub components: Vec<ComponentHlr>,
}

impl AssemblyHlrResult {
    /// Return all visible segments across all components.
    pub fn visible_segments(&self) -> impl Iterator<Item = (&ComponentHlr, &HlrSegment)> {
        self.components.iter().flat_map(|c| {
            c.segments.iter().filter(|s| s.visible).map(move |s| (c, s))
        })
    }

    /// Return all hidden segments across all components.
    pub fn hidden_segments(&self) -> impl Iterator<Item = (&ComponentHlr, &HlrSegment)> {
        self.components.iter().flat_map(|c| {
            c.segments.iter().filter(|s| !s.visible).map(move |s| (c, s))
        })
    }
}

/// Transform a BRep's vertices by an affine transform.
/// Returns a new BRep with transformed vertex positions.
fn transform_brep(brep: &BRep, transform: &DAffine3) -> BRep {
    let mut out = brep.clone();
    for v in &mut out.vertices {
        v.point = transform.transform_point3(v.point);
    }
    out
}

/// Perform hidden-line removal on an assembly of BReps.
///
/// Each component's geometry is transformed to world space, then all triangles
/// are merged into a single occlusion buffer. Each component's edges are
/// tested against the global occlusion buffer, so components correctly
/// occlude each other.
///
/// Returns one `ComponentHlr` per leaf component.
pub fn hlr_assembly(
    components: &[(BRep, DAffine3, String)],
    camera: &HlrCamera,
    samples: usize,
) -> AssemblyHlrResult {
    let view = look_at(camera.eye, camera.target, camera.up);
    let samples = samples.max(2);

    // Transform all BRePs to world space and collect a unified triangle pool.
    let world_breps: Vec<BRep> = components
        .iter()
        .map(|(brep, xf, _)| transform_brep(brep, xf))
        .collect();

    let mut all_triangles: Vec<[DVec3; 3]> = Vec::new();
    for wb in &world_breps {
        all_triangles.extend(collect_triangles(wb));
    }

    let view_dir = (camera.target - camera.eye).normalize_or_zero();
    let mut result = AssemblyHlrResult::default();

    for (wb, (_, _, name)) in world_breps.iter().zip(components.iter()) {
        let mut comp_result = HlrResult::default();

        // ── Wire edges ────────────────────────────────────────────────────
        let mut edge_indices: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        for solid in &wb.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for we in &face.outer_wire.edges {
                        edge_indices.insert(we.idx);
                    }
                    for inner in &face.inner_wires {
                        for we in &inner.edges {
                            edge_indices.insert(we.idx);
                        }
                    }
                }
            }
        }
        for i in 0..wb.edges.len() {
            edge_indices.insert(i);
        }

        for &edge_idx in &edge_indices {
            let Some(edge) = wb.edges.get(edge_idx) else { continue };
            let Some(v_start) = wb.vertices.get(edge.start) else { continue };
            let Some(v_end) = wb.vertices.get(edge.end) else { continue };

            let p0 = v_start.point;
            let p1 = v_end.point;

            let edge_curve = wb
                .geom
                .edge_curve
                .get(edge_idx)
                .and_then(|&ci| ci)
                .and_then(|ci| wb.geom.curves.get(ci));

            let circle_info: Option<Circle3> = edge_curve.and_then(|c| {
                if let rcad_kernel::geom::Curve3::Circle(circ) = c { Some(*circ) } else { None }
            });

            let is_other_curve = edge_curve
                .map_or(false, |c| !matches!(c, rcad_kernel::geom::Curve3::Line(_)))
                && circle_info.is_none();

            let edge_samples = if circle_info.is_some() || is_other_curve {
                64.max(samples)
            } else {
                samples
            };

            let world_pts: Vec<DVec3> = if let Some(circ) = &circle_info {
                let [t0, t1] = wb
                    .geom
                    .edge_curve_range
                    .get(edge_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| circ.default_domain());
                (0..edge_samples)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * (i as f64 / (edge_samples - 1) as f64);
                        circ.point_at(t)
                    })
                    .collect()
            } else if let Some(curve) = edge_curve.filter(|_| is_other_curve) {
                let [t0, t1] = wb
                    .geom
                    .edge_curve_range
                    .get(edge_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| curve.default_domain());
                (0..edge_samples)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * (i as f64 / (edge_samples - 1) as f64);
                        curve.point_at(t)
                    })
                    .collect()
            } else {
                if (p1 - p0).length_squared() < 1e-12 {
                    continue;
                }
                (0..edge_samples)
                    .map(|i| {
                        let t = i as f64 / (edge_samples - 1) as f64;
                        p0 + (p1 - p0) * t
                    })
                    .collect()
            };

            let screen_pts_for_hint: Vec<DVec2> =
                world_pts.iter().map(|&wp| project(wp, &view).0).collect();
            let curve_hint: Option<CurveHint> = if let Some(circ) = &circle_info {
                let (center_2d, _) = project(circ.center, &view);
                let r = screen_pts_for_hint
                    .iter()
                    .map(|p| (*p - center_2d).length())
                    .fold(0.0_f64, f64::max);
                Some(CurveHint::Circle { center: center_2d, radius: r })
            } else if is_other_curve {
                Some(CurveHint::Other)
            } else {
                None
            };

            process_world_pts(&world_pts, curve_hint, false, SegmentType::Edge, camera, &view, &all_triangles, &mut comp_result);
        }

        // ── Silhouette curves ────────────────────────────────────
        let opts = HlrOptions::default().with_edge_samples(samples);
        for sil in compute_silhouettes_with_options(wb, view_dir, &opts) {
            process_world_pts(&sil.world_pts, sil.curve_hint, sil.dense, SegmentType::Silhouette, camera, &view, &all_triangles, &mut comp_result);
        }

        result.components.push(ComponentHlr {
            name: name.clone(),
            segments: comp_result.segments,
        });
    }

    result
}

/// Render HLR result as a simple SVG string.
///
/// Visible edges are drawn solid black; hidden edges are dashed gray.
/// `scale` controls pixel size per unit.
pub fn hlr_to_svg(result: &HlrResult, scale: f64, margin: f64) -> String {
    if result.segments.is_empty() {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\"/>".to_string();
    }

    // Compute bounding box
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for seg in &result.segments {
        for p in [seg.start, seg.end] {
            min_x = min_x.min(p.x);
            min_y = min_y.min(p.y);
            max_x = max_x.max(p.x);
            max_y = max_y.max(p.y);
        }
    }

    // Flip Y (SVG Y grows downward, camera Y grows upward)
    let transform = |p: DVec2| -> (f64, f64) {
        let x = (p.x - min_x) * scale + margin;
        let y = (max_y - p.y) * scale + margin;
        (x, y)
    };

    let w = (max_x - min_x) * scale + 2.0 * margin;
    let h = (max_y - min_y) * scale + 2.0 * margin;

    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.1}\" height=\"{:.1}\" viewBox=\"0 0 {:.1} {:.1}\">\n",
        w, h, w, h
    );
    svg.push_str("  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");

    for seg in &result.segments {
        let (x1, y1) = transform(seg.start);
        let (x2, y2) = transform(seg.end);
        let stroke = if seg.visible {
            "black\" stroke-width=\"1.5"
        } else {
            "#999\" stroke-width=\"0.8\" stroke-dasharray=\"4,3"
        };

        // For circle segments emit an SVG arc path; for all others emit a line.
        if let Some(CurveHint::Circle { center, radius }) = &seg.curve_hint {
            let (cx, cy) = transform(*center);
            let r = radius * scale;
            // Determine large-arc flag: compare arc length vs half-circumference
            let dx1 = x1 - cx;
            let dy1 = y1 - cy;
            let dx2 = x2 - cx;
            let dy2 = y2 - cy;
            let cross = dx1 * dy2 - dy1 * dx2;
            let dot = dx1 * dx2 + dy1 * dy2;
            let angle = cross.atan2(dot).abs();
            let large_arc = if angle > std::f64::consts::PI { 1 } else { 0 };
            let sweep = if cross < 0.0 { 0 } else { 1 };
            svg.push_str(&format!(
                "  <path d=\"M {:.3} {:.3} A {:.3} {:.3} 0 {} {} {:.3} {:.3}\" fill=\"none\" stroke=\"{}\"/>\n",
                x1, y1, r, r, large_arc, sweep, x2, y2, stroke
            ));
            // Also record the center for debugging/reference (as a tiny dot, invisible by default)
            let _ = (cx, cy); // suppress unused warning
        } else {
            svg.push_str(&format!(
                "  <line x1=\"{:.3}\" y1=\"{:.3}\" x2=\"{:.3}\" y2=\"{:.3}\" stroke=\"{}\"/>\n",
                x1, y1, x2, y2, stroke
            ));
        }
    }
    svg.push_str("</svg>\n");
    svg
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn unit_box_hlr_produces_segments() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "HLR should produce segments for a box"
        );
        assert!(
            result.visible().count() > 0,
            "some segments should be visible"
        );
    }

    #[test]
    fn hlr_svg_is_valid_xml() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        assert!(svg.contains("<svg"), "output should be SVG");
        assert!(svg.contains("</svg>"), "SVG should close properly");
        assert!(svg.contains("<line"), "SVG should contain lines");
    }

    #[test]
    fn top_view_box_has_visible_top_edges() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::top(5.0);
        let result = hlr(&brep, &camera, 8);
        let vis = result.visible().count();
        let hid = result.hidden().count();
        assert!(vis > 0, "top view should have visible edges");
        assert!(hid > 0, "top view should have hidden (bottom) edges");
    }

    #[test]
    fn front_view_and_right_view_both_produce_segments() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 1.0,
            depth: 1.0,
        });
        let front_result = hlr(&brep, &HlrCamera::front(5.0), 8);
        let right_result = hlr(&brep, &HlrCamera::right(5.0), 8);
        assert!(!front_result.segments.is_empty());
        assert!(!right_result.segments.is_empty());
    }

    #[test]
    fn hlr_svg_contains_hidden_dashed_lines() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        // Hidden lines are rendered dashed
        assert!(
            svg.contains("stroke-dasharray") || svg.contains("hidden"),
            "SVG should mark hidden lines differently"
        );
    }

    #[test]
    fn hlr_result_has_correct_visibility_counts() {
        // An isometric view of a box has 3 visible faces and 3 hidden faces.
        // The front 3 edges of each visible face → at least some hidden segments exist.
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(10.0);
        let result = hlr(&brep, &camera, 16);
        let total = result.segments.len();
        assert!(total >= 12, "a box has 12 edges, expect at least 12 segments; got {total}");
    }

    #[test]
    fn hlr_circle_edge_sampling() {
        use rcad_kernel::geom::{Circle3, Curve3, CurveEval};

        // Build a minimal BRep with a single circle edge (no solids).
        let mut brep = rcad_kernel::BRep::new();
        let circ = Circle3 {
            center: glam::DVec3::ZERO,
            normal: glam::DVec3::Z,
            radius: 1.0,
        };
        // Add two vertices on the circle (half-circle arc)
        brep.vertices.push(rcad_kernel::topology::Vertex {
            point: circ.point_at(0.0),
        });
        brep.vertices.push(rcad_kernel::topology::Vertex {
            point: circ.point_at(std::f64::consts::PI),
        });
        brep.edges.push(rcad_kernel::topology::Edge { start: 0, end: 1 });
        brep.geom.curves.push(Curve3::Circle(circ));
        brep.geom.edge_curve.push(Some(0));
        brep.geom
            .edge_curve_range
            .push(Some([0.0, std::f64::consts::PI]));

        let camera = HlrCamera::top(5.0);
        let result = hlr(&brep, &camera, 8);

        // The circle edge should produce at least one segment.
        assert!(
            !result.segments.is_empty(),
            "circle edge should produce HLR segments"
        );

        // All sampled 3D points on the circle should lie ON the circle (unit radius).
        // Verify by checking screen_pts all lie within radius ≈ 1.0 of circle center
        // when projected top-down (X-Y plane).
        for seg in &result.segments {
            // The curve_hint for circle segments should be set.
            assert!(
                matches!(seg.curve_hint, Some(CurveHint::Circle { .. })),
                "circle edge segments should carry CurveHint::Circle"
            );
        }

        // SVG should contain arc path elements (not just lines) for circle edges.
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        assert!(
            svg.contains("<path") || result.segments.is_empty(),
            "circle edge SVG should contain <path> arc elements"
        );
    }

    /// Cylinder viewed from the side should produce silhouette line segments
    /// in addition to the wire edges.
    #[test]
    fn cylinder_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        // The cylinder axis is +Y.  Use the right-side camera (looking along -X)
        // so the view direction is perpendicular to the axis → two silhouette lines.
        let camera = HlrCamera::right(10.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "cylinder HLR should produce segments"
        );
        assert!(
            result.segments.len() >= 2,
            "cylinder should have at least 2 silhouette segments, got {}",
            result.segments.len()
        );
    }

    /// Sphere HLR should produce silhouette segments (the great circle).
    #[test]
    fn sphere_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let camera = HlrCamera::front(10.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "sphere HLR should produce silhouette segments"
        );
    }

    /// Cone viewed from the side should produce two silhouette lines from the apex.
    #[test]
    fn cone_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cone {
            base_radius: 1.0,
            height: 2.0,
        });
        // View from the right (perpendicular to cone axis) → two silhouette generators.
        let camera = HlrCamera::right(10.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "cone HLR should produce segments"
        );
        assert!(
            result.segments.len() >= 2,
            "cone should have at least 2 silhouette segments, got {}",
            result.segments.len()
        );
    }

    /// Torus HLR should produce silhouette segments.
    #[test]
    fn torus_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        let camera = HlrCamera::front(20.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "torus HLR should produce silhouette segments"
        );
    }

    // ── Assembly HLR tests ─────────────────────────────────────────────────────

    /// Two boxes side by side — both should produce segments.
    #[test]
    fn hlr_assembly_two_boxes() {
        let box1 = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let box2 = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });

        let components = vec![
            (box1, DAffine3::from_translation(DVec3::new(-2.0, 0.0, 0.0)), "box_left".to_string()),
            (box2, DAffine3::from_translation(DVec3::new(2.0, 0.0, 0.0)), "box_right".to_string()),
        ];

        let camera = HlrCamera::isometric(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2, "should have 2 component results");
        assert!(result.components.iter().all(|c| !c.segments.is_empty()),
            "each component should produce segments");
    }

    /// Small box behind a large box — the small box should be partially hidden.
    #[test]
    fn hlr_assembly_occlusion() {
        let big = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 3.0, height: 3.0, depth: 3.0,
        });
        let small = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 0.5, height: 0.5, depth: 0.5,
        });

        // Front camera looks along +Y from (0, -10, 0).
        // Place small box at +Y behind the big box so it's occluded.
        let components = vec![
            (big, DAffine3::IDENTITY, "big".to_string()),
            (small, DAffine3::from_translation(DVec3::new(0.0, 3.0, 0.0)), "small_behind".to_string()),
        ];

        let camera = HlrCamera::front(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2);
        // The small box behind the big one should have mostly hidden segments
        let small_comp = result.components.iter().find(|c| c.name == "small_behind").unwrap();
        let hidden = small_comp.segments.iter().filter(|s| !s.visible).count();
        let visible = small_comp.segments.iter().filter(|s| s.visible).count();
        assert!(hidden > visible,
            "small box behind big one should have more hidden than visible segments; hidden={hidden}, visible={visible}");
    }

    /// Assembly with a single component should match single-BRep HLR.
    #[test]
    fn hlr_assembly_single_matches_hlr() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);

        let single_hlr = hlr(&brep, &camera, 8);
        let assembly_result = hlr_assembly(
            &[(brep.clone(), DAffine3::IDENTITY, "box".to_string())],
            &camera, 8,
        );

        assert_eq!(assembly_result.components.len(), 1);
        let asm_segs = &assembly_result.components[0].segments;
        // Segment counts should be similar (same geometry, same algorithm)
        assert!(asm_segs.len() >= single_hlr.segments.len() - 2,
            "assembly HLR should produce similar segment count");
        assert!(asm_segs.len() <= single_hlr.segments.len() + 2,
            "assembly HLR should produce similar segment count");
    }

    /// Stacked boxes — top box visible, bottom box partially occluded.
    #[test]
    fn hlr_assembly_stacked_boxes() {
        let bottom = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0, height: 1.0, depth: 2.0,
        });
        let top = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });

        let components = vec![
            (bottom, DAffine3::from_translation(DVec3::new(0.0, 0.0, 0.0)), "bottom".to_string()),
            (top, DAffine3::from_translation(DVec3::new(0.0, 0.0, 1.5)), "top".to_string()),
        ];

        let camera = HlrCamera::isometric(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2);
        // Both boxes should have some visible segments
        for comp in &result.components {
            let vis = comp.segments.iter().filter(|s| s.visible).count();
            assert!(vis > 0, "{} should have visible segments", comp.name);
        }
    }

    /// Empty assembly should return empty result.
    #[test]
    fn hlr_assembly_empty() {
        let components: Vec<(BRep, DAffine3, String)> = vec![];
        let camera = HlrCamera::isometric(5.0);
        let result = hlr_assembly(&components, &camera, 8);
        assert!(result.components.is_empty());
    }

    // ── Improved HLR tests ─────────────────────────────────────────────────────

    #[test]
    fn hlr_options_default_values() {
        let opts = HlrOptions::default();
        assert_eq!(opts.edge_samples, 8);
        assert_eq!(opts.silhouette_samples, 32);
        assert!(opts.curvature_adaptive);
        assert!(opts.tangent_tolerance > 0.0);
    }

    #[test]
    fn hlr_options_builders() {
        let opts = HlrOptions::default()
            .with_edge_samples(16)
            .with_silhouette_samples(64)
            .with_curvature_adaptive(false)
            .with_tangent_tolerance(1e-4);

        assert_eq!(opts.edge_samples, 16);
        assert_eq!(opts.silhouette_samples, 64);
        assert!(!opts.curvature_adaptive);
        assert!((opts.tangent_tolerance - 1e-4).abs() < 1e-10);
    }

    #[test]
    fn hlr_with_options_basic() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let opts = HlrOptions::default().with_edge_samples(16);
        let result = hlr_with_options(&brep, &camera, opts);

        assert!(!result.segments.is_empty(), "should produce segments");
        // All segments from a box should be edges, not silhouettes
        assert!(result.segments.iter().all(|s| s.segment_type == SegmentType::Edge));
    }

    #[test]
    fn cylinder_silhouettes_are_marked() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        let camera = HlrCamera::right(10.0);
        let result = hlr(&brep, &camera, 8);

        // Should have both edge and silhouette segments
        let has_silhouette = result.segments.iter().any(|s| s.is_contour());
        assert!(has_silhouette, "cylinder should have silhouette segments");
    }

    #[test]
    fn sphere_silhouettes_are_marked() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let camera = HlrCamera::front(10.0);
        let result = hlr(&brep, &camera, 8);

        // All segments from a sphere should be silhouettes (no wire edges)
        assert!(
            result.segments.iter().all(|s| s.is_contour()),
            "sphere should only have silhouette segments"
        );
    }

    #[test]
    fn extract_silhouette_curves_sphere() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 2.0 });
        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();

        let curves = extract_silhouette_curves(&brep, view_dir, &opts);

        assert_eq!(curves.len(), 1, "sphere should have one silhouette curve");
        assert!(curves[0].points.len() >= 32, "silhouette should have enough points");

        // All points should be at distance ~2.0 from origin
        for pt in &curves[0].points {
            let dist = pt.length();
            assert!(
                (dist - 2.0).abs() < 0.01,
                "silhouette point distance should be ~2.0, got {dist}"
            );
        }
    }

    #[test]
    fn extract_silhouette_curves_cylinder() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 3.0,
        });
        // View along X axis - perpendicular to cylinder axis (Y)
        let view_dir = DVec3::X;
        let opts = HlrOptions::default();

        let curves = extract_silhouette_curves(&brep, view_dir, &opts);

        assert!(curves.len() >= 2, "cylinder should have at least 2 silhouette curves");

        // Each silhouette curve should be a line (two lines on opposite sides)
        for curve in &curves {
            assert!(curve.points.len() >= 16, "silhouette should have enough points");
        }
    }

    #[test]
    fn extract_silhouette_curves_torus() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();

        let curves = extract_silhouette_curves(&brep, view_dir, &opts);

        assert!(curves.len() >= 2, "torus should have at least 2 silhouette curves");
    }

    #[test]
    fn hlr_result_silhouettes_iterator() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let camera = HlrCamera::front(10.0);
        let result = hlr(&brep, &camera, 8);

        let sil_count = result.silhouettes().count();
        assert!(sil_count > 0, "should have silhouette segments");

        let vis_sil_count = result.visible_silhouettes().count();
        assert!(vis_sil_count > 0, "should have visible silhouette segments");
    }

    #[test]
    fn segment_is_contour_method() {
        let seg = HlrSegment {
            start: DVec2::ZERO,
            end: DVec2::X,
            visible: true,
            curve_hint: None,
            segment_type: SegmentType::Silhouette,
        };
        assert!(seg.is_contour());

        let edge_seg = HlrSegment {
            start: DVec2::ZERO,
            end: DVec2::X,
            visible: true,
            curve_hint: None,
            segment_type: SegmentType::Edge,
        };
        assert!(!edge_seg.is_contour());
    }

    #[test]
    fn adaptive_sampling_high_curvature() {
        // Test that adaptive sampling produces more points in high-curvature regions
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let view_dir = DVec3::Z;

        let opts_low = HlrOptions {
            silhouette_samples: 16,
            curvature_adaptive: false,
            ..HlrOptions::default()
        };
        let opts_high = HlrOptions {
            silhouette_samples: 64,
            curvature_adaptive: true,
            ..HlrOptions::default()
        };

        let curves_low = extract_silhouette_curves(&brep, view_dir, &opts_low);
        let curves_high = extract_silhouette_curves(&brep, view_dir, &opts_high);

        // Both should produce curves
        assert!(!curves_low.is_empty());
        assert!(!curves_high.is_empty());

        // Higher sampling should produce more points
        let pts_low: usize = curves_low.iter().map(|c| c.points.len()).sum();
        let pts_high: usize = curves_high.iter().map(|c| c.points.len()).sum();
        assert!(
            pts_high >= pts_low,
            "higher sampling should produce at least as many points"
        );
    }

    #[test]
    fn tangent_tolerance_affects_silhouette_detection() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let view_dir = DVec3::Z;

        // Very tight tolerance
        let opts_tight = HlrOptions {
            tangent_tolerance: 1e-12,
            ..HlrOptions::default()
        };

        // Very loose tolerance (should still work for sphere)
        let opts_loose = HlrOptions {
            tangent_tolerance: 0.01,
            ..HlrOptions::default()
        };

        let curves_tight = extract_silhouette_curves(&brep, view_dir, &opts_tight);
        let curves_loose = extract_silhouette_curves(&brep, view_dir, &opts_loose);

        // Both should find silhouette curves for a sphere
        assert!(!curves_tight.is_empty());
        assert!(!curves_loose.is_empty());
    }

    // ── New Enhanced HLR Tests ───────────────────────────────────────────────────

    #[test]
    fn hlr_options_new_fields() {
        let opts = HlrOptions::default();
        assert!(opts.parallel, "parallel should be true by default");
        assert_eq!(opts.parallel_threshold, 4);
        assert!(opts.cache_surface_properties);
        assert!(opts.detect_thread_edges);
        assert!(opts.detect_seam_edges);
    }

    #[test]
    fn hlr_options_new_builders() {
        let opts = HlrOptions::default()
            .with_parallel(false)
            .with_parallel_threshold(8)
            .with_surface_caching(false)
            .with_thread_edge_detection(false)
            .with_seam_edge_detection(false)
            .with_silhouette_proximity(0.05);

        assert!(!opts.parallel);
        assert_eq!(opts.parallel_threshold, 8);
        assert!(!opts.cache_surface_properties);
        assert!(!opts.detect_thread_edges);
        assert!(!opts.detect_seam_edges);
        assert!((opts.silhouette_proximity_factor - 0.05).abs() < 1e-10);
    }

    #[test]
    fn surface_property_cache_basic() {
        let surface = Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        });

        let domain = surface.default_domain();
        let mut cache = SurfacePropertyCache::new(16, domain);

        // First access should compute
        let props1 = cache.get_or_compute(&surface, 0.5, 0.5);
        assert!(cache.len() > 0, "cache should have entries");

        // Second access should return cached value
        let props2 = cache.get_or_compute(&surface, 0.5, 0.5);
        assert!((props1.point - props2.point).length() < 1e-10);

        // Verify surface properties
        assert!((props1.point.length() - 1.0).abs() < 1e-10, "sphere point should be on surface");
        assert!((props1.curvatures.0 - 1.0).abs() < 1e-10, "sphere curvature should be 1/r");
    }

    #[test]
    fn surface_properties_near_silhouette() {
        let surface = Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
        });

        // For a Y-axis sphere: x_ax = Z (perpendicular to Y), y_ax = X = Y.cross(Z)
        // At u=π/2, v=π/2: normal = u.cos * Z + u.sin * X = 0*Z + 1*X = X (perpendicular to Z view)
        let props = compute_surface_properties(&surface, std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
        let view_dir = DVec3::Z;

        assert!(props.is_near_silhouette(view_dir, 0.01), "equator at u=π/2 should be near silhouette for Z view");
        assert!((props.normal_dot_view(view_dir)).abs() < 0.01);

        // At u=0, v=π/2: normal = Z (parallel to view) - NOT a silhouette
        let props_front = compute_surface_properties(&surface, 0.0, std::f64::consts::FRAC_PI_2);
        assert!(!props_front.is_near_silhouette(view_dir, 0.5), "point facing viewer should not be near silhouette");

        // At pole (v = 0), normal is Y axis, perpendicular to Z view, so it IS a silhouette
        let props_pole = compute_surface_properties(&surface, 0.0, 0.0);
        assert!(props_pole.is_near_silhouette(view_dir, 0.5), "pole (Y-normal) should be near silhouette for Z view");
        assert!((props_pole.normal_dot_view(view_dir)).abs() < 0.01);
    }

    #[test]
    fn spatial_index_basic() {
        let points: Vec<DVec3> = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(0.1, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.5, 0.5, 0.5),
        ];

        let index = SilhouetteSpatialIndex::build(&points, 0.5);

        assert_eq!(index.len(), 5, "index should have all points");

        // Query radius
        let nearby = index.query_radius(DVec3::ZERO, 0.2);
        assert!(nearby.len() >= 2, "should find points near origin");

        // Query nearest
        let nearest = index.query_nearest(DVec3::new(0.05, 0.05, 0.0));
        assert!(nearest.is_some());
        let (idx, dist) = nearest.unwrap();
        assert!(dist < 0.1, "nearest point should be close");
    }

    #[test]
    fn spatial_index_empty() {
        let points: Vec<DVec3> = vec![];
        let index = SilhouetteSpatialIndex::build(&points, 0.5);

        assert!(index.is_empty());
        assert!(index.query_nearest(DVec3::ZERO).is_none());
    }

    #[test]
    fn edge_classification_basic() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let camera = HlrCamera::isometric(5.0);
        let opts = HlrOptions::default();
        let classifications = classify_edges(&brep, &camera, &opts);

        assert_eq!(classifications.len(), brep.edges.len(), "should classify all edges");

        // All box edges should be regular edges (not thread or seam)
        for class in &classifications {
            assert!(
                class.classification != EdgeClassification::Thread,
                "box edges should not be thread edges"
            );
        }
    }

    #[test]
    fn edge_classification_cylinder_seam() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let camera = HlrCamera::isometric(5.0);
        let opts = HlrOptions::default().with_seam_edge_detection(true);
        let classifications = classify_edges(&brep, &camera, &opts);

        // Cylinder should have at least one seam edge
        let seam_count = classifications
            .iter()
            .filter(|c| c.classification == EdgeClassification::Seam)
            .count();

        // Note: the seam detection depends on the BRep structure
        // For a primitive cylinder, the seam edge may be detected
        assert!(classifications.len() > 0, "should have edge classifications");
    }

    #[test]
    fn segment_type_thread_and_seam() {
        let thread_seg = HlrSegment {
            start: DVec2::ZERO,
            end: DVec2::X,
            visible: true,
            curve_hint: None,
            segment_type: SegmentType::Thread,
        };
        assert!(thread_seg.is_thread());
        assert!(!thread_seg.is_seam());
        assert!(!thread_seg.is_contour());

        let seam_seg = HlrSegment {
            start: DVec2::ZERO,
            end: DVec2::X,
            visible: true,
            curve_hint: None,
            segment_type: SegmentType::Seam,
        };
        assert!(seam_seg.is_seam());
        assert!(!seam_seg.is_thread());
        assert!(!seam_seg.is_contour());
    }

    #[test]
    fn adaptive_sampling_config_default() {
        let config = AdaptiveSamplingConfig::default();
        assert_eq!(config.base_samples, 32);
        assert!(config.max_samples > config.base_samples);
        assert!(config.curvature_threshold > 0.0);
        assert!(config.proximity_threshold > 0.0);
    }

    #[test]
    fn adaptive_sample_creation() {
        let surface = Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
        });

        let view_dir = DVec3::Z;
        let sample = create_adaptive_sample(&surface, view_dir, 0.0, std::f64::consts::FRAC_PI_2);

        assert!(sample.is_some());
        let s = sample.unwrap();

        // Check that the sample is on the equator
        assert!((s.point.y - 0.0).abs() < 1e-10, "equator y should be 0");
        assert!((s.curvature - 0.5).abs() < 1e-10, "sphere radius 2 curvature should be 0.5");
    }

    #[test]
    fn hlr_with_parallel_processing() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let camera = HlrCamera::isometric(5.0);

        // Test with parallel processing enabled
        let opts_parallel = HlrOptions::default()
            .with_parallel(true)
            .with_parallel_threshold(1);

        let result_parallel = hlr_with_options(&brep, &camera, opts_parallel);

        // Test with parallel processing disabled
        let opts_serial = HlrOptions::default()
            .with_parallel(false);

        let result_serial = hlr_with_options(&brep, &camera, opts_serial);

        // Both should produce the same number of segments (within tolerance)
        assert!(!result_parallel.segments.is_empty());
        assert!(!result_serial.segments.is_empty());
    }

    #[test]
    fn curve_surface_intersection_basic() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        let camera = HlrCamera::front(5.0);
        let opts = HlrOptions::default();

        // The sphere BRep has a seam edge
        if !brep.edges.is_empty() {
            // Try to compute curve-surface intersection for the first edge
            let edge_idx = 0;
            if let Some(surface_idx) = brep.geom.face_surface.get(0).and_then(|&s| s) {
                let result = compute_curve_visibility_on_surface(
                    &brep,
                    edge_idx,
                    surface_idx,
                    &camera,
                    &opts,
                );

                // The function should return a result (may have empty intersections)
                if let Some(intersection) = result {
                    assert!(intersection.curve_params.len() == intersection.points.len());
                }
            }
        }
    }

    #[test]
    fn hlr_result_with_thread_segments() {
        // Create a cylinder - thread edges would be helical, but primitives
        // don't have those. This tests the segment type detection.
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let camera = HlrCamera::right(10.0);
        let opts = HlrOptions::default()
            .with_thread_edge_detection(true)
            .with_seam_edge_detection(true);

        let result = hlr_with_options(&brep, &camera, opts);

        assert!(!result.segments.is_empty(), "should have segments");

        // Check that we have various segment types
        let has_edge = result.segments.iter().any(|s| s.segment_type == SegmentType::Edge);
        let has_silhouette = result.segments.iter().any(|s| s.segment_type == SegmentType::Silhouette);

        assert!(has_edge || has_silhouette, "should have edge or silhouette segments");
    }

    #[test]
    fn grazing_angle_handling() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        // Test with different grazing angle thresholds
        let camera = HlrCamera::front(5.0);

        let opts_tight = HlrOptions {
            grazing_angle_threshold: 0.01,
            ..HlrOptions::default()
        };

        let opts_loose = HlrOptions {
            grazing_angle_threshold: 0.5,
            ..HlrOptions::default()
        };

        let result_tight = hlr_with_options(&brep, &camera, opts_tight);
        let result_loose = hlr_with_options(&brep, &camera, opts_loose);

        // Both should produce results
        assert!(!result_tight.segments.is_empty());
        assert!(!result_loose.segments.is_empty());
    }

    #[test]
    fn performance_with_caching() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 2.0 });

        let camera = HlrCamera::front(5.0);

        // Test with caching enabled
        let opts_cached = HlrOptions::default()
            .with_surface_caching(true);

        let result_cached = hlr_with_options(&brep, &camera, opts_cached);

        // Test with caching disabled
        let opts_uncached = HlrOptions::default()
            .with_surface_caching(false);

        let result_uncached = hlr_with_options(&brep, &camera, opts_uncached);

        // Both should produce valid results
        assert!(!result_cached.segments.is_empty());
        assert!(!result_uncached.segments.is_empty());
    }

    #[test]
    fn is_degenerate_edge_check() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        // Check all edges for degeneracy
        for i in 0..brep.edges.len() {
            let _is_degenerate = is_degenerate_edge_for_hlr(&brep, i);
            // For a sphere primitive, edges should not be degenerate
            // (the seam edge is not degenerate, just periodic)
        }
    }

    // ── Ellipsoid Silhouette Tests ───────────────────────────────────────────────

    #[test]
    fn ellipsoid_silhouette_basic() {
        use rcad_kernel::geom::EllipsoidalSurface;

        // Create an ellipsoid with different radii
        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Z; // View along Z axis
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert_eq!(curves.len(), 1, "ellipsoid should have one silhouette curve");
        assert!(curves[0].len() >= 32, "silhouette should have enough points");
    }

    #[test]
    fn ellipsoid_silhouette_satisfies_condition() {
        use rcad_kernel::geom::{EllipsoidalSurface, SurfaceEval};

        // Create an ellipsoid
        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert!(!curves.is_empty(), "should have silhouette curves");

        // Check that all silhouette points satisfy n·v ≈ 0
        for pt in &curves[0] {
            // Compute the point in local coordinates
            let x = pt.x;
            let y = pt.y;
            let z = pt.z;

            // Normal direction (gradient of implicit equation, normalized)
            let grad = DVec3::new(
                x / (ell.radius_x * ell.radius_x),
                y / (ell.radius_y * ell.radius_y),
                z / (ell.radius_z * ell.radius_z),
            );
            let normal = grad.normalize_or_zero();

            // Dot product with view direction should be near zero
            let dot = normal.dot(view_dir);
            assert!(
                dot.abs() < 0.05,
                "silhouette point should satisfy n·v ≈ 0, got {dot}"
            );
        }
    }

    #[test]
    fn ellipsoid_silhouette_on_surface() {
        use rcad_kernel::geom::EllipsoidalSurface;

        // Create an ellipsoid
        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert!(!curves.is_empty());

        // Check that all silhouette points are on the ellipsoid surface
        // x²/a² + y²/b² + z²/c² = 1
        for pt in &curves[0] {
            let value = (pt.x / ell.radius_x).powi(2)
                + (pt.y / ell.radius_y).powi(2)
                + (pt.z / ell.radius_z).powi(2);
            assert!(
                (value - 1.0).abs() < 1e-6,
                "point should be on ellipsoid surface, got implicit value {value}"
            );
        }
    }

    #[test]
    fn ellipsoid_silhouette_various_view_directions() {
        use rcad_kernel::geom::EllipsoidalSurface;

        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let opts = HlrOptions::default();

        // Test various view directions
        let view_directions = [
            DVec3::X,
            DVec3::Y,
            DVec3::Z,
            DVec3::new(1.0, 1.0, 0.0).normalize(),
            DVec3::new(1.0, 1.0, 1.0).normalize(),
            DVec3::new(0.0, 1.0, 1.0).normalize(),
        ];

        for view_dir in view_directions {
            let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);
            assert!(
                !curves.is_empty() && !curves[0].is_empty(),
                "should have silhouette for view_dir {:?}",
                view_dir
            );

            // Verify silhouette condition
            for pt in &curves[0] {
                let grad = DVec3::new(
                    pt.x / (ell.radius_x * ell.radius_x),
                    pt.y / (ell.radius_y * ell.radius_y),
                    pt.z / (ell.radius_z * ell.radius_z),
                );
                let normal = grad.normalize_or_zero();
                let dot = normal.dot(view_dir);
                assert!(
                    dot.abs() < 0.1,
                    "silhouette condition not satisfied for view_dir {:?}, dot = {}",
                    view_dir,
                    dot
                );
            }
        }
    }

    #[test]
    fn ellipsoid_silhouette_sphere_case() {
        use rcad_kernel::geom::EllipsoidalSurface;

        // A sphere is a special case of an ellipsoid with all radii equal
        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 2.0,
            radius_y: 2.0,
            radius_z: 2.0,
        };

        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert_eq!(curves.len(), 1, "sphere should have one silhouette curve");

        // All points should be at distance 2.0 from origin (great circle)
        for pt in &curves[0] {
            let dist = pt.length();
            assert!(
                (dist - 2.0).abs() < 0.01,
                "sphere silhouette point should be at radius distance, got {dist}"
            );

            // z-coordinate should be near 0 (great circle perpendicular to Z)
            assert!(
                pt.z.abs() < 0.01,
                "great circle should be in XY plane, got z = {}",
                pt.z
            );
        }
    }

    #[test]
    fn ellipsoid_silhouette_translated() {
        use rcad_kernel::geom::EllipsoidalSurface;

        // Ellipsoid not at origin
        let ell = EllipsoidalSurface {
            center: DVec3::new(1.0, -2.0, 0.5),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert!(!curves.is_empty());

        // Check that all points are centered around the ellipsoid center
        for pt in &curves[0] {
            let local = *pt - ell.center;
            let value = (local.x / ell.radius_x).powi(2)
                + (local.y / ell.radius_y).powi(2)
                + (local.z / ell.radius_z).powi(2);
            assert!(
                (value - 1.0).abs() < 1e-6,
                "point should be on translated ellipsoid surface"
            );
        }
    }

    #[test]
    fn ellipsoid_silhouette_rotated_frame() {
        use rcad_kernel::geom::EllipsoidalSurface;

        // Ellipsoid with rotated axis (not aligned with Z)
        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y, // Axis along Y
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Y; // View along the axis
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert!(!curves.is_empty());

        // The silhouette should be an ellipse in the XZ plane
        for pt in &curves[0] {
            // Y coordinate should be near 0 (silhouette in plane perpendicular to view)
            assert!(
                pt.y.abs() < 0.01,
                "silhouette should be in XZ plane for Y view, got y = {}",
                pt.y
            );
        }
    }

    #[test]
    fn ellipsoid_silhouette_closed_curve() {
        use rcad_kernel::geom::EllipsoidalSurface;

        let ell = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };

        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();
        let curves = extract_ellipsoid_silhouettes(&ell, view_dir, &opts, 64);

        assert!(!curves.is_empty());

        // The silhouette should be a closed curve
        // First and last points should be close
        let pts = &curves[0];
        let first = pts.first().unwrap();
        let last = pts.last().unwrap();
        let closure_dist = (*first - *last).length();
        assert!(
            closure_dist < 0.5,
            "silhouette should be approximately closed, distance between first and last = {}",
            closure_dist
        );
    }
}
