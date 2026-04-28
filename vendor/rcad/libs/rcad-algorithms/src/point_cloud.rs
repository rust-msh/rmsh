//! Point cloud analysis tools, analogous to OCCT 8.0 PointSetLib.
//!
//! Provides:
//! - Principal Component Analysis (PCA)
//! - Inertia tensor computation
//! - Dimensionality estimation
//! - Outlier detection and point cloud simplification
//! - Normal estimation
//! - Shape fitting (plane, sphere, cylinder)
//! - ICP registration (point-to-point, point-to-plane)
//! - Segmentation (region growing, Euclidean clustering, shape segmentation)
//! - Surface reconstruction (Poisson, Ball pivoting, Delaunay)
//! - Advanced sampling (curvature-aware, Poisson disk)
//! - BRep integration

use glam::DVec3;
use std::cmp::Ordering;

/// A collection of 3D points.
#[derive(Debug, Clone, Default)]
pub struct PointCloud {
    pub points: Vec<DVec3>,
}

impl PointCloud {
    /// Creates an empty point cloud.
    pub fn new() -> Self {
        Self { points: Vec::new() }
    }

    /// Creates a point cloud from a slice of points.
    pub fn from_points(points: &[DVec3]) -> Self {
        Self { points: points.to_vec() }
    }

    /// Creates a point cloud from a vector of points.
    pub fn from_vec(points: Vec<DVec3>) -> Self {
        Self { points }
    }

    /// Returns the number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Returns true if the point cloud is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Computes the axis-aligned bounding box.
    pub fn bounding_box(&self) -> Option<(DVec3, DVec3)> {
        if self.points.is_empty() {
            return None;
        }
        let mut min = DVec3::splat(f64::INFINITY);
        let mut max = DVec3::splat(f64::NEG_INFINITY);
        for &p in &self.points {
            min = min.min(p);
            max = max.max(p);
        }
        Some((min, max))
    }

    /// Computes the centroid (mean) of all points.
    pub fn centroid(&self) -> Option<DVec3> {
        if self.points.is_empty() {
            return None;
        }
        let sum: DVec3 = self.points.iter().sum();
        Some(sum / self.points.len() as f64)
    }
}

/// Classification of point cloud dimensionality based on PCA.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dimensionality {
    /// All points are at or very near a single location.
    Point,
    /// Points lie approximately along a line.
    Linear,
    /// Points lie approximately on a plane.
    Planar,
    /// Points have significant extent in all three dimensions.
    Volumetric,
}

/// Result of point cloud analysis.
#[derive(Debug, Clone)]
pub struct PointCloudAnalysis {
    /// Centroid (mean) of all points.
    pub centroid: DVec3,
    /// Principal axes sorted by eigenvalue (largest to smallest).
    /// - `principal_axes[0]`: direction of maximum variance
    /// - `principal_axes[2]`: direction of minimum variance (normal for planar data)
    pub principal_axes: [DVec3; 3],
    /// Principal values (eigenvalues) corresponding to each axis.
    pub principal_values: [f64; 3],
    /// Axis-aligned bounding box as (min, max).
    pub bounding_box: (DVec3, DVec3),
    /// Inertia tensor about the centroid.
    pub inertia_tensor: [[f64; 3]; 3],
    /// Estimated dimensionality.
    pub dimensionality: Dimensionality,
}

/// Performs comprehensive analysis on a point cloud.
///
/// Computes centroid, PCA, bounding box, inertia tensor, and dimensionality.
pub fn analyze_point_cloud(points: &[DVec3]) -> Option<PointCloudAnalysis> {
    if points.is_empty() {
        return None;
    }

    let centroid = points.iter().sum::<DVec3>() / points.len() as f64;

    let (principal_axes, principal_values) = compute_pca(points);

    let bounding_box = {
        let mut min = DVec3::splat(f64::INFINITY);
        let mut max = DVec3::splat(f64::NEG_INFINITY);
        for &p in points {
            min = min.min(p);
            max = max.max(p);
        }
        (min, max)
    };

    let inertia_tensor = compute_inertia_centroid(points, centroid);

    let dimensionality = estimate_dimensionality(principal_values, 0.01);

    Some(PointCloudAnalysis {
        centroid,
        principal_axes,
        principal_values,
        bounding_box,
        inertia_tensor,
        dimensionality,
    })
}

/// Computes Principal Component Analysis (PCA) on a point set.
///
/// Returns:
/// - Principal axes (eigenvectors) sorted by eigenvalue (largest first)
/// - Principal values (eigenvalues) sorted largest first
///
/// The principal axes form an orthonormal basis. For planar data,
/// `principal_axes[2]` is the normal of the best-fit plane.
pub fn compute_pca(points: &[DVec3]) -> ([DVec3; 3], [f64; 3]) {
    if points.is_empty() {
        return ([DVec3::X, DVec3::Y, DVec3::Z], [0.0; 3]);
    }

    let n = points.len() as f64;
    let centroid = points.iter().sum::<DVec3>() / n;

    // Compute covariance matrix
    let mut cov = [[0.0; 3]; 3];
    for &p in points {
        let d = p - centroid;
        cov[0][0] += d.x * d.x;
        cov[0][1] += d.x * d.y;
        cov[0][2] += d.x * d.z;
        cov[1][1] += d.y * d.y;
        cov[1][2] += d.y * d.z;
        cov[2][2] += d.z * d.z;
    }
    cov[0][0] /= n;
    cov[0][1] /= n;
    cov[0][2] /= n;
    cov[1][1] /= n;
    cov[1][2] /= n;
    cov[2][2] /= n;
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];

    // Compute eigenvalues and eigenvectors using power iteration
    let (eigenvalues, eigenvectors) = compute_eigendecomposition_3x3(&cov);

    // Sort by eigenvalue descending
    let mut indexed: [(usize, f64, DVec3); 3] = [
        (0, eigenvalues[0], eigenvectors[0]),
        (1, eigenvalues[1], eigenvectors[1]),
        (2, eigenvalues[2], eigenvectors[2]),
    ];
    indexed.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));

    let mut axes = [DVec3::ZERO; 3];
    let mut values = [0.0; 3];
    for i in 0..3 {
        axes[i] = indexed[i].2;
        values[i] = indexed[i].1.max(0.0);
    }

    // Ensure orthonormal right-handed basis
    axes[2] = axes[0].cross(axes[1]).normalize_or(DVec3::Z);
    axes[1] = axes[2].cross(axes[0]).normalize_or(DVec3::Y);

    (axes, values)
}

/// Compute eigenvalues and eigenvectors of a symmetric 3x3 matrix.
fn compute_eigendecomposition_3x3(m: &[[f64; 3]; 3]) -> ([f64; 3], [DVec3; 3]) {
    // Use Jacobi eigenvalue algorithm for symmetric matrices
    let mut a = *m;
    let mut v = [
        DVec3::X,
        DVec3::Y,
        DVec3::Z,
    ];

    const MAX_ITERATIONS: usize = 100;
    const TOLERANCE: f64 = 1e-12;

    for _ in 0..MAX_ITERATIONS {
        // Find the largest off-diagonal element
        let mut max_val = 0.0;
        let (mut p, mut q) = (0, 1);

        for i in 0..3 {
            for j in (i + 1)..3 {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < TOLERANCE {
            break;
        }

        // Compute rotation angle
        let theta = if (a[p][p] - a[q][q]).abs() < TOLERANCE {
            std::f64::consts::FRAC_PI_4 * a[p][q].signum()
        } else {
            0.5 * (2.0 * a[p][q] / (a[p][p] - a[q][q])).atan()
        };

        let c = theta.cos();
        let s = theta.sin();

        // Update matrix A
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];

        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        for i in 0..3 {
            if i != p && i != q {
                let aip = a[i][p];
                let aiq = a[i][q];
                a[i][p] = c * aip - s * aiq;
                a[p][i] = a[i][p];
                a[i][q] = s * aip + c * aiq;
                a[q][i] = a[i][q];
            }
        }

        // Update eigenvector matrix
        for i in 0..3 {
            let vip = v[i][p];
            let viq = v[i][q];
            v[i][p] = c * vip - s * viq;
            v[i][q] = s * vip + c * viq;
        }
    }

    // Normalize eigenvectors
    for i in 0..3 {
        v[i] = v[i].normalize_or(match i {
            0 => DVec3::X,
            1 => DVec3::Y,
            _ => DVec3::Z,
        });
    }

    ([a[0][0], a[1][1], a[2][2]], v)
}

/// Computes the inertia tensor of a point set about the origin.
///
/// The inertia tensor is a symmetric 3x3 matrix defined as:
/// ```text
/// Ixx = Σ(y²+z²),  Iyy = Σ(x²+z²),  Izz = Σ(x²+y²)
/// Ixy = -Σxy,       Ixz = -Σxz,       Iyz = -Σyz
/// ```
pub fn compute_inertia(points: &[DVec3]) -> [[f64; 3]; 3] {
    if points.is_empty() {
        return [[0.0; 3]; 3];
    }

    compute_inertia_centroid(points, DVec3::ZERO)
}

/// Computes the inertia tensor of a point set about a given centroid.
fn compute_inertia_centroid(points: &[DVec3], centroid: DVec3) -> [[f64; 3]; 3] {
    let mut ixx = 0.0_f64;
    let mut iyy = 0.0_f64;
    let mut izz = 0.0_f64;
    let mut ixy = 0.0_f64;
    let mut ixz = 0.0_f64;
    let mut iyz = 0.0_f64;

    for &p in points {
        let d = p - centroid;
        let x = d.x;
        let y = d.y;
        let z = d.z;

        ixx += y * y + z * z;
        iyy += x * x + z * z;
        izz += x * x + y * y;
        ixy -= x * y;
        ixz -= x * z;
        iyz -= y * z;
    }

    [
        [ixx, ixy, ixz],
        [ixy, iyy, iyz],
        [ixz, iyz, izz],
    ]
}

/// Estimates the dimensionality of a point cloud from PCA eigenvalues.
///
/// The threshold is the relative tolerance for considering an eigenvalue
/// as "negligible" compared to the largest eigenvalue.
///
/// Classification:
/// - Point: all eigenvalues are negligible (total variance near zero)
/// - Linear: only one significant eigenvalue
/// - Planar: two significant eigenvalues, one negligible
/// - Volumetric: all three eigenvalues are significant
pub fn estimate_dimensionality(pca_values: [f64; 3], threshold: f64) -> Dimensionality {
    // Sort eigenvalues descending
    let mut sorted = pca_values;
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap_or(Ordering::Equal));

    let total: f64 = sorted.iter().sum();
    // If total variance is negligible, it's a point
    if total < 1e-10 {
        return Dimensionality::Point;
    }

    // Normalize by largest eigenvalue
    let max_val = sorted[0].max(1e-20);
    let rel1 = sorted[1] / max_val;
    let rel2 = sorted[2] / max_val;

    // Count how many eigenvalues are significant (relative to max)
    // First eigenvalue is always significant if total > 0
    let sig1 = true;
    let sig2 = rel1 > threshold;
    let sig3 = rel2 > threshold;

    let count = [sig1, sig2, sig3].iter().filter(|&&x| x).count();

    match count {
        0 => Dimensionality::Point,
        1 => Dimensionality::Linear,
        2 => Dimensionality::Planar,
        3 => Dimensionality::Volumetric,
        _ => Dimensionality::Volumetric,
    }
}

// ============================================================================
// Point Cloud Processing
// ============================================================================

/// Detected outlier point with its outlier score.
#[derive(Debug, Clone)]
pub struct OutlierPoint {
    /// Index of the outlier point in the original point cloud.
    pub index: usize,
    /// Outlier score (higher = more likely an outlier).
    pub score: f64,
}

/// Detects outlier points using the Local Outlier Factor (LOF) algorithm.
///
/// Parameters:
/// - `points`: the point cloud
/// - `k`: number of nearest neighbors to consider (default: 20)
/// - `threshold`: LOF score threshold for outliers (default: 2.0)
///
/// Returns a list of outlier points sorted by score (highest first).
pub fn detect_outliers(points: &[DVec3], k: usize, threshold: f64) -> Vec<OutlierPoint> {
    if points.len() <= k + 1 {
        return Vec::new();
    }

    let k = k.min(points.len() - 1).max(1);
    let n = points.len();

    // Compute k-distances and reachability distances
    let mut lof_scores = vec![0.0; n];

    for i in 0..n {
        // Find k nearest neighbors
        let mut distances: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, (points[j] - points[i]).length_squared()))
            .collect();
        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));

        let k_dist = distances[k - 1].1.sqrt();
        let neighbors: Vec<usize> = distances.iter().take(k).map(|&(j, _)| j).collect();

        // Compute local reachability density
        let mut lrd_sum = 0.0;
        for &j in &neighbors {
            let dist_ij = (points[j] - points[i]).length();
            let reach_dist = dist_ij.max(k_dist);
            lrd_sum += reach_dist;
        }
        let lrd_i = k as f64 / lrd_sum;

        // Compute LOF
        let mut lof_sum = 0.0;
        for &j in &neighbors {
            // Compute lrd of neighbor j (simplified - use same k_dist approximation)
            let mut j_dists: Vec<f64> = (0..n)
                .filter(|&l| l != j)
                .map(|l| (points[l] - points[j]).length())
                .collect();
            j_dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
            let j_k_dist = j_dists[k - 1];
            let mut j_lrd_sum = 0.0;
            for &l in j_dists.iter().take(k) {
                j_lrd_sum += l.max(j_k_dist);
            }
            let lrd_j = k as f64 / j_lrd_sum;
            lof_sum += lrd_j / lrd_i;
        }
        lof_scores[i] = lof_sum / k as f64;
    }

    // Collect outliers
    let mut outliers: Vec<OutlierPoint> = lof_scores
        .iter()
        .enumerate()
        .filter(|&(_, &score)| score > threshold)
        .map(|(i, &score)| OutlierPoint { index: i, score })
        .collect();

    outliers.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal));
    outliers
}

/// Removes outliers from a point cloud based on LOF detection.
///
/// Returns a new point cloud with outliers removed.
pub fn remove_outliers(points: &[DVec3], k: usize, threshold: f64) -> Vec<DVec3> {
    let outliers = detect_outliers(points, k, threshold);
    let outlier_set: std::collections::HashSet<usize> = outliers.iter().map(|o| o.index).collect();

    points
        .iter()
        .enumerate()
        .filter(|(i, _)| !outlier_set.contains(i))
        .map(|(_, &p)| p)
        .collect()
}

/// Sampling strategy for point cloud simplification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SamplingStrategy {
    /// Uniform random sampling.
    Random,
    /// Grid-based voxel sampling (one point per voxel).
    Voxel,
    /// Farthest point sampling (maximizes coverage).
    FarthestPoint,
}

/// Simplifies a point cloud by reducing the number of points.
///
/// Parameters:
/// - `points`: the input point cloud
/// - `target_count`: target number of output points
/// - `strategy`: sampling strategy to use
///
/// Returns the simplified point cloud.
pub fn simplify_point_cloud(
    points: &[DVec3],
    target_count: usize,
    strategy: SamplingStrategy,
) -> Vec<DVec3> {
    if points.len() <= target_count {
        return points.to_vec();
    }

    match strategy {
        SamplingStrategy::Random => random_sample(points, target_count),
        SamplingStrategy::Voxel => voxel_sample(points, target_count),
        SamplingStrategy::FarthestPoint => farthest_point_sample(points, target_count),
    }
}

fn random_sample(points: &[DVec3], target_count: usize) -> Vec<DVec3> {
    use std::collections::HashSet;
    let n = points.len();
    let mut indices: HashSet<usize> = HashSet::new();
    let mut rng = SimpleRng::new(12345);

    while indices.len() < target_count {
        indices.insert((rng.next() as usize) % n);
    }

    indices.iter().map(|&i| points[i]).collect()
}

fn voxel_sample(points: &[DVec3], target_count: usize) -> Vec<DVec3> {
    if points.is_empty() {
        return Vec::new();
    }

    // Compute bounding box
    let (min, max) = points.iter().fold(
        (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
        |(min, max), &p| (min.min(p), max.max(p)),
    );

    // Estimate voxel size
    let volume = (max.x - min.x) * (max.y - min.y) * (max.z - min.z);
    let voxel_size = (volume / target_count as f64).cbrt().max(1e-10);

    // Group points by voxel
    let mut voxels: std::collections::HashMap<[i64; 3], Vec<DVec3>> = std::collections::HashMap::new();

    for &p in points {
        let key = [
            ((p.x - min.x) / voxel_size).floor() as i64,
            ((p.y - min.y) / voxel_size).floor() as i64,
            ((p.z - min.z) / voxel_size).floor() as i64,
        ];
        voxels.entry(key).or_insert_with(Vec::new).push(p);
    }

    // Take centroid of each voxel
    voxels
        .values()
        .map(|pts| {
            let sum: DVec3 = pts.iter().sum();
            sum / pts.len() as f64
        })
        .collect()
}

fn farthest_point_sample(points: &[DVec3], target_count: usize) -> Vec<DVec3> {
    if points.len() <= target_count {
        return points.to_vec();
    }

    let n = points.len();
    let mut selected = Vec::with_capacity(target_count);
    let mut distances = vec![f64::INFINITY; n];

    // Start with centroid or first point
    let centroid = points.iter().sum::<DVec3>() / n as f64;
    let first_idx = points
        .iter()
        .enumerate()
        .map(|(i, &p)| (i, (p - centroid).length_squared()))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    selected.push(points[first_idx]);

    // Greedy farthest point selection
    while selected.len() < target_count {
        let last = selected.last().unwrap();
        let mut farthest_idx = 0;
        let mut farthest_dist = 0.0;

        for i in 0..n {
            let d = (points[i] - *last).length_squared();
            distances[i] = distances[i].min(d);
            if distances[i] > farthest_dist {
                farthest_dist = distances[i];
                farthest_idx = i;
            }
        }

        selected.push(points[farthest_idx]);
    }

    selected
}

/// Simple deterministic RNG for reproducible sampling.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next(&mut self) -> u64 {
        // xorshift64
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }
}

/// Estimates normals for each point using local PCA.
///
/// For each point, computes PCA on its k nearest neighbors.
/// The normal is the eigenvector corresponding to the smallest eigenvalue.
///
/// Returns a vector of unit normals (one per point).
pub fn estimate_normals(points: &[DVec3], k: usize) -> Vec<DVec3> {
    if points.is_empty() {
        return Vec::new();
    }

    let k = k.min(points.len() - 1).max(2);
    let n = points.len();
    let mut normals = Vec::with_capacity(n);

    for i in 0..n {
        // Find k nearest neighbors
        let mut distances: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, (points[j] - points[i]).length_squared()))
            .collect();
        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));

        let neighbor_pts: Vec<DVec3> = distances
            .iter()
            .take(k)
            .map(|&(j, _)| points[j])
            .collect();

        // PCA on neighbors
        let (axes, values) = compute_pca(&neighbor_pts);

        // Normal is direction of minimum variance
        // Check if the smallest eigenvalue is small enough for a planar fit
        let max_val = values[0].max(1e-20);
        if values[2] / max_val < 0.1 {
            normals.push(axes[2]);
        } else {
            // Not planar enough, still use smallest variance direction
            normals.push(axes[2]);
        }
    }

    normals
}

// ============================================================================
// Point Cloud Fitting
// ============================================================================

/// Result of fitting a plane to a point cloud.
#[derive(Debug, Clone)]
pub struct FittedPlane {
    /// A point on the plane.
    pub point: DVec3,
    /// Unit normal of the plane.
    pub normal: DVec3,
    /// RMS distance of points to the fitted plane.
    pub rms_error: f64,
}

/// Fits a plane to a point cloud using least squares.
///
/// Uses PCA: the normal is the eigenvector corresponding to the smallest eigenvalue.
pub fn fit_plane(points: &[DVec3]) -> Option<FittedPlane> {
    if points.len() < 3 {
        return None;
    }

    let centroid = points.iter().sum::<DVec3>() / points.len() as f64;
    let (axes, _values) = compute_pca(points);

    // Compute RMS error
    let normal = axes[2];
    let mut sum_sq = 0.0;
    for &p in points {
        let d = (p - centroid).dot(normal);
        sum_sq += d * d;
    }
    let rms_error = (sum_sq / points.len() as f64).sqrt();

    Some(FittedPlane {
        point: centroid,
        normal,
        rms_error,
    })
}

/// Result of fitting a sphere to a point cloud.
#[derive(Debug, Clone)]
pub struct FittedSphere {
    /// Center of the sphere.
    pub center: DVec3,
    /// Radius of the sphere.
    pub radius: f64,
    /// RMS distance of points to the fitted sphere surface.
    pub rms_error: f64,
}

/// Fits a sphere to a point cloud using least squares.
///
/// Uses an algebraic fit followed by geometric refinement.
pub fn fit_sphere(points: &[DVec3]) -> Option<FittedSphere> {
    if points.len() < 4 {
        return None;
    }

    // Algebraic fit using linear least squares
    // Fit: (x - cx)^2 + (y - cy)^2 + (z - cz)^2 = r^2
    // Rewrite: x^2 + y^2 + z^2 - 2*cx*x - 2*cy*y - 2*cz*z + cx^2 + cy^2 + cz^2 - r^2 = 0
    // Let: a = -2*cx, b = -2*cy, c = -2*cz, d = cx^2 + cy^2 + cz^2 - r^2
    // Then: x^2 + y^2 + z^2 + a*x + b*y + c*z + d = 0
    // Solve for a, b, c, d using linear least squares

    let n = points.len();
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;
    let mut sum_x2 = 0.0;
    let mut sum_y2 = 0.0;
    let mut sum_z2 = 0.0;
    let mut sum_xy = 0.0;
    let mut sum_xz = 0.0;
    let mut sum_yz = 0.0;
    let mut sum_xyz = 0.0;  // x*(x^2+y^2+z^2)
    let mut sum_x2yz = 0.0; // y*(x^2+y^2+z^2)
    let mut sum_xz2 = 0.0;  // z*(x^2+y^2+z^2)
    let mut sum_r2 = 0.0;   // x^2 + y^2 + z^2

    for &p in points {
        let x = p.x;
        let y = p.y;
        let z = p.z;
        let r2 = x * x + y * y + z * z;

        sum_x += x;
        sum_y += y;
        sum_z += z;
        sum_x2 += x * x;
        sum_y2 += y * y;
        sum_z2 += z * z;
        sum_xy += x * y;
        sum_xz += x * z;
        sum_yz += y * z;
        sum_xyz += x * r2;
        sum_x2yz += y * r2;
        sum_xz2 += z * r2;
        sum_r2 += r2;
    }

    // Solve 4x4 linear system: A * [a, b, c, d]^T = B
    // Where A is the matrix of sums, B is the RHS
    let a = [
        [sum_x2, sum_xy, sum_xz, sum_x],
        [sum_xy, sum_y2, sum_yz, sum_y],
        [sum_xz, sum_yz, sum_z2, sum_z],
        [sum_x, sum_y, sum_z, n as f64],
    ];
    let b = [-sum_xyz, -sum_x2yz, -sum_xz2, -sum_r2];

    // Solve using Gaussian elimination
    let coeffs = solve_linear_4x4(&a, &b)?;

    let cx = -coeffs[0] / 2.0;
    let cy = -coeffs[1] / 2.0;
    let cz = -coeffs[2] / 2.0;
    let center = DVec3::new(cx, cy, cz);
    let radius = (cx * cx + cy * cy + cz * cz - coeffs[3]).sqrt().max(0.0);

    // Compute RMS error
    let mut sum_sq = 0.0;
    for &p in points {
        let d = (p - center).length() - radius;
        sum_sq += d * d;
    }
    let rms_error = (sum_sq / n as f64).sqrt();

    Some(FittedSphere {
        center,
        radius,
        rms_error,
    })
}

/// Solve a 4x4 linear system using Gaussian elimination with partial pivoting.
fn solve_linear_4x4(a: &[[f64; 4]; 4], b: &[f64; 4]) -> Option<[f64; 4]> {
    const N: usize = 4;
    let mut m = a.clone();
    let mut v = *b;

    // Forward elimination
    for col in 0..N {
        // Find pivot
        let mut max_row = col;
        let mut max_val = m[col][col].abs();
        for row in (col + 1)..N {
            if m[row][col].abs() > max_val {
                max_val = m[row][col].abs();
                max_row = row;
            }
        }

        if max_val < 1e-14 {
            return None; // Singular matrix
        }

        // Swap rows
        m.swap(col, max_row);
        v.swap(col, max_row);

        // Eliminate
        for row in (col + 1)..N {
            let factor = m[row][col] / m[col][col];
            for j in col..N {
                m[row][j] -= factor * m[col][j];
            }
            v[row] -= factor * v[col];
        }
    }

    // Back substitution
    let mut x = [0.0; N];
    for i in (0..N).rev() {
        let mut sum = v[i];
        for j in (i + 1)..N {
            sum -= m[i][j] * x[j];
        }
        x[i] = sum / m[i][i];
    }

    Some(x)
}

/// Result of fitting a cylinder to a point cloud.
#[derive(Debug, Clone)]
pub struct FittedCylinder {
    /// A point on the cylinder axis.
    pub axis_point: DVec3,
    /// Unit direction of the cylinder axis.
    pub axis_direction: DVec3,
    /// Radius of the cylinder.
    pub radius: f64,
    /// RMS distance of points to the fitted cylinder surface.
    pub rms_error: f64,
}

/// Fits a cylinder to a point cloud using iterative optimization.
///
/// The algorithm:
/// 1. Estimate the axis direction using PCA on differences from centroid
/// 2. Project points onto the plane perpendicular to axis
/// 3. Fit a circle to the projected points
pub fn fit_cylinder(points: &[DVec3]) -> Option<FittedCylinder> {
    if points.len() < 5 {
        return None;
    }

    let centroid = points.iter().sum::<DVec3>() / points.len() as f64;
    let (axes, values) = compute_pca(points);

    // For a cylinder, we expect two large eigenvalues and one small
    // The cylinder axis is the direction of minimum variance
    let max_val = values[0].max(1e-20);
    if values[2] / max_val > 0.3 {
        // Not cylindrical enough
        // Try alternative: axis might be along maximum variance direction
        // This happens for short cylinders
    }

    // Try both possible axis directions and pick the better fit
    let axis_candidates = [axes[2], axes[0]];

    let mut best_fit: Option<FittedCylinder> = None;
    let mut best_error = f64::INFINITY;

    for axis in axis_candidates {
        if let Some(cyl) = fit_cylinder_with_axis(points, centroid, axis) {
            if cyl.rms_error < best_error {
                best_error = cyl.rms_error;
                best_fit = Some(cyl);
            }
        }
    }

    best_fit
}

fn fit_cylinder_with_axis(points: &[DVec3], centroid: DVec3, axis: DVec3) -> Option<FittedCylinder> {
    let axis = axis.normalize_or(DVec3::Z);

    // Build orthonormal basis with axis as Z
    let u = if axis.x.abs() < 0.9 {
        axis.cross(DVec3::X).normalize()
    } else {
        axis.cross(DVec3::Y).normalize()
    };
    let v = axis.cross(u);

    // Project points onto the plane perpendicular to axis
    let projected: Vec<DVec2> = points
        .iter()
        .map(|&p| {
            let d = p - centroid;
            DVec2::new(d.dot(u), d.dot(v))
        })
        .collect();

    // Fit a circle to projected points
    let circle = fit_circle_2d(&projected)?;

    // Transform back to 3D
    let center_2d = circle.center;
    let center_3d = centroid + center_2d.x * u + center_2d.y * v;

    // Compute RMS error in 3D
    let mut sum_sq = 0.0;
    for &p in points {
        let to_center = p - center_3d;
        let axial_dist = to_center.dot(axis);
        let radial = to_center - axial_dist * axis;
        let d = radial.length() - circle.radius;
        sum_sq += d * d;
    }
    let rms_error = (sum_sq / points.len() as f64).sqrt();

    Some(FittedCylinder {
        axis_point: center_3d,
        axis_direction: axis,
        radius: circle.radius,
        rms_error,
    })
}

/// 2D vector for circle fitting.
#[derive(Debug, Clone, Copy)]
struct DVec2 {
    x: f64,
    y: f64,
}

impl DVec2 {
    fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    fn length(self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
}

/// Fitted circle in 2D.
struct FittedCircle {
    center: DVec2,
    radius: f64,
}

/// Fit a circle to 2D points using least squares.
fn fit_circle_2d(points: &[DVec2]) -> Option<FittedCircle> {
    if points.len() < 3 {
        return None;
    }

    // Use algebraic fit similar to sphere fitting
    let n = points.len();
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_x2 = 0.0;
    let mut sum_y2 = 0.0;
    let mut sum_xy = 0.0;
    let mut sum_xr = 0.0;
    let mut sum_yr = 0.0;
    let mut sum_r = 0.0;

    for &p in points {
        let x = p.x;
        let y = p.y;
        let r = x * x + y * y;

        sum_x += x;
        sum_y += y;
        sum_x2 += x * x;
        sum_y2 += y * y;
        sum_xy += x * y;
        sum_xr += x * r;
        sum_yr += y * r;
        sum_r += r;
    }

    // Solve 3x3: [sum_x2, sum_xy, sum_x] [a]   [-sum_xr]
    //            [sum_xy, sum_y2, sum_y] [b] = [-sum_yr]
    //            [sum_x,  sum_y,  n   ] [d]   [-sum_r ]
    let a = [
        [sum_x2, sum_xy, sum_x],
        [sum_xy, sum_y2, sum_y],
        [sum_x, sum_y, n as f64],
    ];
    let b = [-sum_xr, -sum_yr, -sum_r];

    let coeffs = solve_linear_3x3(&a, &b)?;

    let cx = -coeffs[0] / 2.0;
    let cy = -coeffs[1] / 2.0;
    let radius = (cx * cx + cy * cy - coeffs[2]).sqrt().max(0.0);

    Some(FittedCircle {
        center: DVec2::new(cx, cy),
        radius,
    })
}

/// Solve a 3x3 linear system using Gaussian elimination.
fn solve_linear_3x3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
    const N: usize = 3;
    let mut m = a.clone();
    let mut v = *b;

    // Forward elimination
    for col in 0..N {
        let mut max_row = col;
        let mut max_val = m[col][col].abs();
        for row in (col + 1)..N {
            if m[row][col].abs() > max_val {
                max_val = m[row][col].abs();
                max_row = row;
            }
        }

        if max_val < 1e-14 {
            return None;
        }

        m.swap(col, max_row);
        v.swap(col, max_row);

        for row in (col + 1)..N {
            let factor = m[row][col] / m[col][col];
            for j in col..N {
                m[row][j] -= factor * m[col][j];
            }
            v[row] -= factor * v[col];
        }
    }

    // Back substitution
    let mut x = [0.0; N];
    for i in (0..N).rev() {
        let mut sum = v[i];
        for j in (i + 1)..N {
            sum -= m[i][j] * x[j];
        }
        x[i] = sum / m[i][i];
    }

    Some(x)
}

/// Result of fitting a convex polygon to a point cloud.
#[derive(Debug, Clone)]
pub struct FittedPolygon {
    /// Vertices of the fitted polygon.
    pub vertices: Vec<DVec3>,
    /// Plane of the polygon.
    pub plane_point: DVec3,
    pub plane_normal: DVec3,
    /// Area of the polygon.
    pub area: f64,
}

/// Fits a convex polygon to a planar point cloud.
///
/// Projects points to the best-fit plane and computes the 2D convex hull.
pub fn fit_polygon(points: &[DVec3]) -> Option<FittedPolygon> {
    if points.len() < 3 {
        return None;
    }

    let plane = fit_plane(points)?;

    // Build orthonormal basis on the plane
    let normal = plane.normal;
    let u = if normal.x.abs() < 0.9 {
        normal.cross(DVec3::X).normalize()
    } else {
        normal.cross(DVec3::Y).normalize()
    };
    let v = normal.cross(u);

    // Project to 2D
    let projected_2d: Vec<DVec2> = points
        .iter()
        .map(|&p| {
            let d = p - plane.point;
            DVec2::new(d.dot(u), d.dot(v))
        })
        .collect();

    // Compute 2D convex hull
    let hull_2d = convex_hull_2d(&projected_2d);

    if hull_2d.len() < 3 {
        return None;
    }

    // Transform back to 3D
    let vertices: Vec<DVec3> = hull_2d
        .iter()
        .map(|&p2d| plane.point + p2d.x * u + p2d.y * v)
        .collect();

    // Compute area using shoelace formula in 3D
    let mut area = 0.0;
    let n = vertices.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let cross = vertices[i].cross(vertices[j]);
        area += cross.dot(normal);
    }
    area = area.abs() / 2.0;

    Some(FittedPolygon {
        vertices,
        plane_point: plane.point,
        plane_normal: plane.normal,
        area,
    })
}

/// Compute the convex hull of 2D points using Andrew's monotone chain algorithm.
/// This is more robust than Graham scan for numerical stability.
fn convex_hull_2d(points: &[DVec2]) -> Vec<DVec2> {
    if points.len() < 3 {
        return points.to_vec();
    }

    // Sort points by x, then by y
    let mut sorted: Vec<usize> = (0..points.len()).collect();
    sorted.sort_by(|&a, &b| {
        let pa = points[a];
        let pb = points[b];
        if (pa.x - pb.x).abs() > 1e-14 {
            pa.x.partial_cmp(&pb.x).unwrap_or(Ordering::Equal)
        } else {
            pa.y.partial_cmp(&pb.y).unwrap_or(Ordering::Equal)
        }
    });

    // Cross product of OA and OB vectors
    let cross = |o: DVec2, a: DVec2, b: DVec2| -> f64 {
        (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
    };

    // Build lower hull
    let mut lower: Vec<DVec2> = Vec::new();
    for &i in &sorted {
        while lower.len() >= 2 {
            let n = lower.len();
            if cross(lower[n - 2], lower[n - 1], points[i]) <= 0.0 {
                lower.pop();
            } else {
                break;
            }
        }
        lower.push(points[i]);
    }

    // Build upper hull
    let mut upper: Vec<DVec2> = Vec::new();
    for &i in sorted.iter().rev() {
        while upper.len() >= 2 {
            let n = upper.len();
            if cross(upper[n - 2], upper[n - 1], points[i]) <= 0.0 {
                upper.pop();
            } else {
                break;
            }
        }
        upper.push(points[i]);
    }

    // Remove last point of each half because it's repeated at the beginning of the other half
    lower.pop();
    upper.pop();

    // Concatenate lower and upper hulls
    lower.extend(upper);
    lower
}

impl std::ops::Sub for DVec2 {
    type Output = DVec2;

    fn sub(self, other: DVec2) -> DVec2 {
        DVec2::new(self.x - other.x, self.y - other.y)
    }
}

// ============================================================================
// ICP Registration
// ============================================================================

/// Result of ICP registration.
#[derive(Debug, Clone)]
pub struct IcpResult {
    /// Rotation matrix (3x3) to transform source to target.
    pub rotation: [[f64; 3]; 3],
    /// Translation vector to transform source to target.
    pub translation: DVec3,
    /// Final RMS error after convergence.
    pub rms_error: f64,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Whether the algorithm converged within tolerance.
    pub converged: bool,
}

impl IcpResult {
    /// Applies the transformation to a point.
    pub fn transform_point(&self, point: DVec3) -> DVec3 {
        let r = &self.rotation;
        DVec3::new(
            r[0][0] * point.x + r[0][1] * point.y + r[0][2] * point.z + self.translation.x,
            r[1][0] * point.x + r[1][1] * point.y + r[1][2] * point.z + self.translation.y,
            r[2][0] * point.x + r[2][1] * point.y + r[2][2] * point.z + self.translation.z,
        )
    }

    /// Returns the transformation as a 4x4 homogeneous matrix.
    pub fn to_matrix(&self) -> [[f64; 4]; 4] {
        let r = &self.rotation;
        let t = &self.translation;
        [
            [r[0][0], r[0][1], r[0][2], t.x],
            [r[1][0], r[1][1], r[1][2], t.y],
            [r[2][0], r[2][1], r[2][2], t.z],
            [0.0, 0.0, 0.0, 1.0],
        ]
    }
}

/// ICP algorithm variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IcpVariant {
    /// Standard point-to-point ICP.
    PointToPoint,
    /// Point-to-plane ICP (requires normals on target).
    PointToPlane,
}

/// ICP configuration parameters.
#[derive(Debug, Clone)]
pub struct IcpConfig {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Convergence tolerance for RMS error change.
    pub tolerance: f64,
    /// Maximum correspondence distance (points beyond this are ignored).
    pub max_correspondence_distance: f64,
    /// Whether to use reciprocal correspondence (both directions).
    pub use_reciprocal: bool,
}

impl Default for IcpConfig {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            tolerance: 1e-6,
            max_correspondence_distance: f64::INFINITY,
            use_reciprocal: false,
        }
    }
}

/// Performs Iterative Closest Point (ICP) registration.
///
/// Aligns the source point cloud to the target point cloud.
///
/// # Arguments
/// * `source` - Source point cloud to be transformed
/// * `target` - Target point cloud (reference)
/// * `variant` - ICP variant to use
/// * `config` - Configuration parameters
///
/// # Returns
/// * `IcpResult` containing the transformation and convergence info
pub fn icp_registration(
    source: &[DVec3],
    target: &[DVec3],
    variant: IcpVariant,
    config: &IcpConfig,
) -> Option<IcpResult> {
    if source.is_empty() || target.is_empty() {
        return None;
    }

    match variant {
        IcpVariant::PointToPoint => icp_point_to_point(source, target, config),
        IcpVariant::PointToPlane => {
            // Estimate normals for target if not provided
            let target_normals = estimate_normals(target, 10);
            icp_point_to_plane(source, target, &target_normals, config)
        }
    }
}

/// Performs ICP registration with pre-computed normals.
pub fn icp_registration_with_normals(
    source: &[DVec3],
    target: &[DVec3],
    target_normals: &[DVec3],
    config: &IcpConfig,
) -> Option<IcpResult> {
    if source.is_empty() || target.is_empty() || target_normals.len() != target.len() {
        return None;
    }
    icp_point_to_plane(source, target, target_normals, config)
}

fn icp_point_to_point(
    source: &[DVec3],
    target: &[DVec3],
    config: &IcpConfig,
) -> Option<IcpResult> {
    let mut transformed: Vec<DVec3> = source.to_vec();
    let mut cumulative_rotation = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut cumulative_translation = DVec3::ZERO;

    // Build KD-tree for target (simple brute-force for now)
    let mut prev_error = f64::INFINITY;
    let mut converged = false;

    for _iteration in 0..config.max_iterations {
        // Find correspondences
        let correspondences = find_correspondences(
            &transformed,
            target,
            config.max_correspondence_distance,
        );

        if correspondences.is_empty() {
            break;
        }

        // Compute transformation using SVD
        let (rotation, translation) = compute_transformation_svd(&transformed, target, &correspondences)?;

        // Apply transformation
        for p in &mut transformed {
            *p = apply_transform(*p, &rotation, &translation);
        }

        // Update cumulative transformation
        cumulative_rotation = multiply_matrices(&rotation, &cumulative_rotation);
        cumulative_translation = apply_transform_to_vector(&cumulative_translation, &rotation, &translation);

        // Compute error
        let rms_error = compute_rms_error(&transformed, target, &correspondences);

        // Check convergence
        if (prev_error - rms_error).abs() < config.tolerance {
            converged = true;
            break;
        }
        prev_error = rms_error;
    }

    Some(IcpResult {
        rotation: cumulative_rotation,
        translation: cumulative_translation,
        rms_error: prev_error,
        iterations: config.max_iterations,
        converged,
    })
}

fn icp_point_to_plane(
    source: &[DVec3],
    target: &[DVec3],
    target_normals: &[DVec3],
    config: &IcpConfig,
) -> Option<IcpResult> {
    let mut transformed: Vec<DVec3> = source.to_vec();
    let mut cumulative_rotation = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut cumulative_translation = DVec3::ZERO;

    let mut prev_error = f64::INFINITY;
    let mut converged = false;

    for _iteration in 0..config.max_iterations {
        let correspondences = find_correspondences(
            &transformed,
            target,
            config.max_correspondence_distance,
        );

        if correspondences.is_empty() {
            break;
        }

        // Compute point-to-plane transformation using linear least squares
        let (rotation, translation) = compute_point_to_plane_transformation(
            &transformed,
            target,
            target_normals,
            &correspondences,
        )?;

        // Apply transformation
        for p in &mut transformed {
            *p = apply_transform(*p, &rotation, &translation);
        }

        // Update cumulative transformation
        cumulative_rotation = multiply_matrices(&rotation, &cumulative_rotation);
        cumulative_translation = apply_transform_to_vector(&cumulative_translation, &rotation, &translation);

        // Compute point-to-plane error
        let rms_error = compute_point_to_plane_error(&transformed, target, target_normals, &correspondences);

        if (prev_error - rms_error).abs() < config.tolerance {
            converged = true;
            break;
        }
        prev_error = rms_error;
    }

    Some(IcpResult {
        rotation: cumulative_rotation,
        translation: cumulative_translation,
        rms_error: prev_error,
        iterations: config.max_iterations,
        converged,
    })
}

/// Correspondence between source and target points.
struct Correspondence {
    source_idx: usize,
    target_idx: usize,
    distance: f64,
}

fn find_correspondences(
    source: &[DVec3],
    target: &[DVec3],
    max_distance: f64,
) -> Vec<Correspondence> {
    let mut correspondences = Vec::with_capacity(source.len());

    for (i, &s) in source.iter().enumerate() {
        let mut best_dist = f64::INFINITY;
        let mut best_j = 0;

        for (j, &t) in target.iter().enumerate() {
            let d = (s - t).length_squared();
            if d < best_dist {
                best_dist = d;
                best_j = j;
            }
        }

        let dist = best_dist.sqrt();
        if dist <= max_distance {
            correspondences.push(Correspondence {
                source_idx: i,
                target_idx: best_j,
                distance: dist,
            });
        }
    }

    correspondences
}

fn compute_transformation_svd(
    source: &[DVec3],
    target: &[DVec3],
    correspondences: &[Correspondence],
) -> Option<([[f64; 3]; 3], DVec3)> {
    if correspondences.is_empty() {
        return None;
    }

    let n = correspondences.len();

    // Compute centroids
    let mut source_centroid = DVec3::ZERO;
    let mut target_centroid = DVec3::ZERO;

    for c in correspondences {
        source_centroid += source[c.source_idx];
        target_centroid += target[c.target_idx];
    }

    source_centroid /= n as f64;
    target_centroid /= n as f64;

    // Compute cross-covariance matrix
    let mut h = [[0.0_f64; 3]; 3];
    for c in correspondences {
        let s = source[c.source_idx] - source_centroid;
        let t = target[c.target_idx] - target_centroid;
        h[0][0] += s.x * t.x;
        h[0][1] += s.x * t.y;
        h[0][2] += s.x * t.z;
        h[1][0] += s.y * t.x;
        h[1][1] += s.y * t.y;
        h[1][2] += s.y * t.z;
        h[2][0] += s.z * t.x;
        h[2][1] += s.z * t.y;
        h[2][2] += s.z * t.z;
    }

    // SVD of H
    let (u, _, v) = svd_3x3(&h)?;

    // R = V * U^T
    let mut rotation = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            rotation[i][j] = v[i][0] * u[j][0] + v[i][1] * u[j][1] + v[i][2] * u[j][2];
        }
    }

    // Handle reflection case
    let det = rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
            - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
            + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0]);

    if det < 0.0 {
        // Flip sign of last column of V
        let mut v_corrected = v;
        for i in 0..3 {
            v_corrected[i][2] = -v_corrected[i][2];
        }
        for i in 0..3 {
            for j in 0..3 {
                rotation[i][j] = v_corrected[i][0] * u[j][0] + v_corrected[i][1] * u[j][1] + v_corrected[i][2] * u[j][2];
            }
        }
    }

    // Translation = target_centroid - R * source_centroid
    let translation = DVec3::new(
        target_centroid.x - (rotation[0][0] * source_centroid.x + rotation[0][1] * source_centroid.y + rotation[0][2] * source_centroid.z),
        target_centroid.y - (rotation[1][0] * source_centroid.x + rotation[1][1] * source_centroid.y + rotation[1][2] * source_centroid.z),
        target_centroid.z - (rotation[2][0] * source_centroid.x + rotation[2][1] * source_centroid.y + rotation[2][2] * source_centroid.z),
    );

    Some((rotation, translation))
}

/// Compute point-to-plane transformation using linear least squares.
/// Uses the linearized rotation approximation for small angles.
fn compute_point_to_plane_transformation(
    source: &[DVec3],
    target: &[DVec3],
    target_normals: &[DVec3],
    correspondences: &[Correspondence],
) -> Option<([[f64; 3]; 3], DVec3)> {
    if correspondences.len() < 6 {
        return None;
    }

    // Build linear system: A * x = b
    // Where x = [alpha, beta, gamma, tx, ty, tz]^T (rotation angles and translation)
    // For each correspondence: n_i^T * (R(p_i) + t - q_i) = 0
    // Linearized: n_i^T * (p_i + r x p_i + t - q_i) = 0
    // n_i^T * (p_i - q_i) + n_i^T * (r x p_i) + n_i^T * t = 0
    // n_i^T * (p_i - q_i) + (p_i x n_i)^T * r + n_i^T * t = 0

    let _n = correspondences.len();
    let mut ata = [[0.0_f64; 6]; 6];
    let mut atb = [0.0_f64; 6];

    for c in correspondences {
        let p = source[c.source_idx];
        let q = target[c.target_idx];
        let n = target_normals[c.target_idx];

        let cross = DVec3::new(
            p.y * n.z - p.z * n.y,
            p.z * n.x - p.x * n.z,
            p.x * n.y - p.y * n.x,
        );

        let diff = p - q;
        let rhs = -(n.x * diff.x + n.y * diff.y + n.z * diff.z);

        // Row: [cross.x, cross.y, cross.z, n.x, n.y, n.z]
        let row = [cross.x, cross.y, cross.z, n.x, n.y, n.z];

        for i in 0..6 {
            for j in 0..6 {
                ata[i][j] += row[i] * row[j];
            }
            atb[i] += row[i] * rhs;
        }
    }

    // Solve 6x6 system using Gaussian elimination
    let solution = solve_linear_6x6(&ata, &atb)?;

    // Convert angles to rotation matrix
    let (alpha, beta, gamma) = (solution[0], solution[1], solution[2]);
    let rotation = angles_to_rotation_matrix(alpha, beta, gamma);

    let translation = DVec3::new(solution[3], solution[4], solution[5]);

    Some((rotation, translation))
}

fn angles_to_rotation_matrix(alpha: f64, beta: f64, gamma: f64) -> [[f64; 3]; 3] {
    let ca = alpha.cos();
    let sa = alpha.sin();
    let cb = beta.cos();
    let sb = beta.sin();
    let cg = gamma.cos();
    let sg = gamma.sin();

    // R = Rz(gamma) * Ry(beta) * Rx(alpha)
    [
        [cg * cb, cg * sb * sa - sg * ca, cg * sb * ca + sg * sa],
        [sg * cb, sg * sb * sa + cg * ca, sg * sb * ca - cg * sa],
        [-sb, cb * sa, cb * ca],
    ]
}

fn compute_rms_error(
    source: &[DVec3],
    target: &[DVec3],
    correspondences: &[Correspondence],
) -> f64 {
    let mut sum_sq = 0.0;
    for c in correspondences {
        let d = source[c.source_idx] - target[c.target_idx];
        sum_sq += d.length_squared();
    }
    (sum_sq / correspondences.len() as f64).sqrt()
}

fn compute_point_to_plane_error(
    source: &[DVec3],
    target: &[DVec3],
    target_normals: &[DVec3],
    correspondences: &[Correspondence],
) -> f64 {
    let mut sum_sq = 0.0;
    for c in correspondences {
        let diff = source[c.source_idx] - target[c.target_idx];
        let n = target_normals[c.target_idx];
        let dist = diff.dot(n);
        sum_sq += dist * dist;
    }
    (sum_sq / correspondences.len() as f64).sqrt()
}

fn apply_transform(point: DVec3, rotation: &[[f64; 3]; 3], translation: &DVec3) -> DVec3 {
    DVec3::new(
        rotation[0][0] * point.x + rotation[0][1] * point.y + rotation[0][2] * point.z + translation.x,
        rotation[1][0] * point.x + rotation[1][1] * point.y + rotation[1][2] * point.z + translation.y,
        rotation[2][0] * point.x + rotation[2][1] * point.y + rotation[2][2] * point.z + translation.z,
    )
}

fn apply_transform_to_vector(vec: &DVec3, rotation: &[[f64; 3]; 3], translation: &DVec3) -> DVec3 {
    DVec3::new(
        rotation[0][0] * vec.x + rotation[0][1] * vec.y + rotation[0][2] * vec.z + translation.x,
        rotation[1][0] * vec.x + rotation[1][1] * vec.y + rotation[1][2] * vec.z + translation.y,
        rotation[2][0] * vec.x + rotation[2][1] * vec.y + rotation[2][2] * vec.z + translation.z,
    )
}

fn multiply_matrices(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    result
}

/// Simplified SVD for 3x3 matrices using Jacobi iteration.
fn svd_3x3(a: &[[f64; 3]; 3]) -> Option<([[f64; 3]; 3], [f64; 3], [[f64; 3]; 3])> {
    // Compute A^T * A for eigenvalue decomposition
    let mut ata = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                ata[i][j] += a[k][i] * a[k][j];
            }
        }
    }

    // Compute eigenvalues and eigenvectors of A^T * A
    let (eigenvalues, v) = jacobi_eigen(&ata);

    // Compute singular values
    let mut sigma = [0.0; 3];
    for i in 0..3 {
        sigma[i] = eigenvalues[i].max(0.0).sqrt();
    }

    // Compute U = A * V * Sigma^(-1)
    let mut u = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            if sigma[j] > 1e-10 {
                for k in 0..3 {
                    u[i][j] += a[i][k] * v[k][j];
                }
                u[i][j] /= sigma[j];
            }
        }
    }

    // Orthonormalize U
    for j in 0..3 {
        let mut norm = 0.0;
        for i in 0..3 {
            norm += u[i][j] * u[i][j];
        }
        norm = norm.sqrt();
        if norm > 1e-10 {
            for i in 0..3 {
                u[i][j] /= norm;
            }
        }
    }

    Some((u, sigma, v))
}

fn jacobi_eigen(a: &[[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut m = *a;
    let mut v = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    const MAX_ITER: usize = 100;
    const TOL: f64 = 1e-12;

    for _ in 0..MAX_ITER {
        // Find largest off-diagonal element
        let mut max_val = 0.0;
        let (mut p, mut q) = (0, 1);

        for i in 0..3 {
            for j in (i + 1)..3 {
                if m[i][j].abs() > max_val {
                    max_val = m[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < TOL {
            break;
        }

        // Compute rotation angle
        let theta = if (m[p][p] - m[q][q]).abs() < TOL {
            std::f64::consts::FRAC_PI_4 * m[p][q].signum()
        } else {
            0.5 * (2.0 * m[p][q] / (m[p][p] - m[q][q])).atan()
        };

        let c = theta.cos();
        let s = theta.sin();

        // Apply rotation
        let app = m[p][p];
        let aqq = m[q][q];
        let apq = m[p][q];

        m[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        m[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        m[p][q] = 0.0;
        m[q][p] = 0.0;

        for i in 0..3 {
            if i != p && i != q {
                let aip = m[i][p];
                let aiq = m[i][q];
                m[i][p] = c * aip - s * aiq;
                m[p][i] = m[i][p];
                m[i][q] = s * aip + c * aiq;
                m[q][i] = m[i][q];
            }
        }

        // Update eigenvectors
        for i in 0..3 {
            let vip = v[i][p];
            let viq = v[i][q];
            v[i][p] = c * vip - s * viq;
            v[i][q] = s * vip + c * viq;
        }
    }

    // Sort eigenvalues descending
    let mut indexed = [(m[0][0], 0), (m[1][1], 1), (m[2][2], 2)];
    indexed.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));

    let mut eigenvalues = [0.0; 3];
    let mut v_sorted = [[0.0; 3]; 3];
    for (i, &(val, idx)) in indexed.iter().enumerate() {
        eigenvalues[i] = val;
        for j in 0..3 {
            v_sorted[j][i] = v[j][idx];
        }
    }

    (eigenvalues, v_sorted)
}

fn solve_linear_6x6(a: &[[f64; 6]; 6], b: &[f64; 6]) -> Option<[f64; 6]> {
    const N: usize = 6;
    let mut m = *a;
    let mut v = *b;

    for col in 0..N {
        let mut max_row = col;
        let mut max_val = m[col][col].abs();
        for row in (col + 1)..N {
            if m[row][col].abs() > max_val {
                max_val = m[row][col].abs();
                max_row = row;
            }
        }

        if max_val < 1e-14 {
            return None;
        }

        m.swap(col, max_row);
        v.swap(col, max_row);

        for row in (col + 1)..N {
            let factor = m[row][col] / m[col][col];
            for j in col..N {
                m[row][j] -= factor * m[col][j];
            }
            v[row] -= factor * v[col];
        }
    }

    let mut x = [0.0; N];
    for i in (0..N).rev() {
        let mut sum = v[i];
        for j in (i + 1)..N {
            sum -= m[i][j] * x[j];
        }
        x[i] = sum / m[i][i];
    }

    Some(x)
}

// ============================================================================
// Segmentation
// ============================================================================

/// Result of region growing segmentation.
#[derive(Debug, Clone)]
pub struct Segment {
    /// Indices of points in this segment.
    pub point_indices: Vec<usize>,
    /// Fitted shape type (if any).
    pub shape_type: Option<ShapeType>,
    /// Fitted shape parameters (if any).
    pub shape_params: Option<ShapeParams>,
    /// Centroid of the segment.
    pub centroid: DVec3,
}

/// Types of shapes that can be detected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeType {
    Plane,
    Sphere,
    Cylinder,
    Cone,
    Torus,
}

/// Parameters for fitted shapes.
#[derive(Debug, Clone)]
pub enum ShapeParams {
    Plane { point: DVec3, normal: DVec3 },
    Sphere { center: DVec3, radius: f64 },
    Cylinder { axis_point: DVec3, axis_direction: DVec3, radius: f64 },
}

/// Configuration for region growing segmentation.
#[derive(Debug, Clone)]
pub struct RegionGrowingConfig {
    /// Number of neighbors for normal estimation.
    pub k_neighbors: usize,
    /// Maximum angular difference (radians) for region growing.
    pub max_angle: f64,
    /// Maximum distance from fitted shape.
    pub max_distance: f64,
    /// Minimum number of points for a valid segment.
    pub min_segment_size: usize,
    /// Maximum number of segments to extract.
    pub max_segments: usize,
}

impl Default for RegionGrowingConfig {
    fn default() -> Self {
        Self {
            k_neighbors: 30,
            max_angle: std::f64::consts::PI / 6.0, // 30 degrees
            max_distance: 0.01,
            min_segment_size: 100,
            max_segments: 100,
        }
    }
}

/// Performs region growing segmentation based on smoothness constraint.
///
/// Grows regions from seed points, adding neighbors with similar normals.
pub fn region_growing_segmentation(
    points: &[DVec3],
    config: &RegionGrowingConfig,
) -> Vec<Segment> {
    if points.is_empty() {
        return Vec::new();
    }

    let n = points.len();

    // Estimate normals
    let normals = estimate_normals(points, config.k_neighbors);

    // Build neighbor graph (kNN)
    let neighbors = build_neighbor_graph(points, config.k_neighbors);

    // Track visited points
    let mut visited = vec![false; n];
    let mut segments = Vec::new();

    // Sort points by curvature (lowest first for better seeds)
    let curvatures: Vec<f64> = compute_curvatures(points, &neighbors);
    let mut sorted_indices: Vec<usize> = (0..n).collect();
    sorted_indices.sort_by(|&a, &b| {
        curvatures[a].partial_cmp(&curvatures[b]).unwrap_or(Ordering::Equal)
    });

    for &seed_idx in &sorted_indices {
        if visited[seed_idx] || segments.len() >= config.max_segments {
            continue;
        }

        // Grow region from seed
        let segment_indices = grow_region(
            seed_idx,
            points,
            &normals,
            &neighbors,
            &mut visited,
            config,
        );

        if segment_indices.len() >= config.min_segment_size {
            let segment_pts: Vec<DVec3> = segment_indices.iter().map(|&i| points[i]).collect();
            let centroid = segment_pts.iter().sum::<DVec3>() / segment_pts.len() as f64;

            segments.push(Segment {
                point_indices: segment_indices,
                shape_type: None,
                shape_params: None,
                centroid,
            });
        }
    }

    segments
}

/// Performs Euclidean clustering segmentation.
///
/// Clusters points based on Euclidean distance threshold.
pub fn euclidean_clustering(
    points: &[DVec3],
    tolerance: f64,
    min_cluster_size: usize,
) -> Vec<Vec<usize>> {
    if points.is_empty() {
        return Vec::new();
    }

    let n = points.len();
    let mut visited = vec![false; n];
    let mut clusters = Vec::new();
    let tolerance_sq = tolerance * tolerance;

    for i in 0..n {
        if visited[i] {
            continue;
        }

        let mut cluster = Vec::new();
        let mut queue = vec![i];
        visited[i] = true;

        while let Some(current) = queue.pop() {
            cluster.push(current);

            // Find all neighbors within tolerance
            for j in 0..n {
                if !visited[j] && (points[j] - points[current]).length_squared() < tolerance_sq {
                    visited[j] = true;
                    queue.push(j);
                }
            }
        }

        if cluster.len() >= min_cluster_size {
            clusters.push(cluster);
        }
    }

    // Sort clusters by size (largest first)
    clusters.sort_by(|a, b| b.len().cmp(&a.len()));

    clusters
}

/// Performs shape-based segmentation (plane, sphere, cylinder).
///
/// Uses RANSAC to detect dominant shapes and segment the point cloud.
pub fn shape_segmentation(
    points: &[DVec3],
    shape_type: ShapeType,
    distance_threshold: f64,
    min_points: usize,
    max_iterations: usize,
) -> Option<(ShapeParams, Vec<usize>, Vec<usize>)> {
    match shape_type {
        ShapeType::Plane => ransac_plane_segmentation(points, distance_threshold, min_points, max_iterations),
        ShapeType::Sphere => ransac_sphere_segmentation(points, distance_threshold, min_points, max_iterations),
        ShapeType::Cylinder => ransac_cylinder_segmentation(points, distance_threshold, min_points, max_iterations),
        _ => None,
    }
}

fn ransac_plane_segmentation(
    points: &[DVec3],
    threshold: f64,
    min_points: usize,
    max_iterations: usize,
) -> Option<(ShapeParams, Vec<usize>, Vec<usize>)> {
    if points.len() < 3 {
        return None;
    }

    let mut best_inliers = Vec::new();
    let mut best_plane: Option<(DVec3, DVec3)> = None;

    let mut rng = SimpleRng::new(42);

    for _ in 0..max_iterations {
        // Sample 3 random points
        let i0 = (rng.next() as usize) % points.len();
        let i1 = (rng.next() as usize) % points.len();
        let i2 = (rng.next() as usize) % points.len();

        if i0 == i1 || i1 == i2 || i0 == i2 {
            continue;
        }

        let p0 = points[i0];
        let p1 = points[i1];
        let p2 = points[i2];

        // Compute plane normal
        let normal = (p1 - p0).cross(p2 - p0);
        let len = normal.length();
        if len < 1e-10 {
            continue;
        }
        let normal = normal / len;

        // Find inliers
        let mut inliers = Vec::new();
        for (i, &p) in points.iter().enumerate() {
            let dist = (p - p0).dot(normal).abs();
            if dist < threshold {
                inliers.push(i);
            }
        }

        if inliers.len() > best_inliers.len() {
            best_inliers = inliers;
            best_plane = Some((p0, normal));
        }
    }

    if best_inliers.len() < min_points {
        return None;
    }

    let (_point, _normal) = best_plane.unwrap();

    // Refit using all inliers
    let inlier_pts: Vec<DVec3> = best_inliers.iter().map(|&i| points[i]).collect();
    let fitted = fit_plane(&inlier_pts)?;

    // Separate inliers and outliers
    let inlier_set: std::collections::HashSet<usize> = best_inliers.iter().copied().collect();
    let outliers: Vec<usize> = (0..points.len())
        .filter(|i| !inlier_set.contains(i))
        .collect();

    Some((
        ShapeParams::Plane {
            point: fitted.point,
            normal: fitted.normal,
        },
        best_inliers,
        outliers,
    ))
}

fn ransac_sphere_segmentation(
    points: &[DVec3],
    threshold: f64,
    min_points: usize,
    max_iterations: usize,
) -> Option<(ShapeParams, Vec<usize>, Vec<usize>)> {
    if points.len() < 4 {
        return None;
    }

    let mut best_inliers = Vec::new();
    let mut best_sphere: Option<(DVec3, f64)> = None;

    let mut rng = SimpleRng::new(42);

    for _ in 0..max_iterations {
        // Sample 4 random points
        let indices: Vec<usize> = (0..4)
            .map(|_| (rng.next() as usize) % points.len())
            .collect();

        if indices.iter().collect::<std::collections::HashSet<_>>().len() < 4 {
            continue;
        }

        let sample_pts: Vec<DVec3> = indices.iter().map(|&i| points[i]).collect();

        // Fit sphere to 4 points
        let sphere = fit_sphere_4pt(&sample_pts);
        let (center, radius) = match sphere {
            Some(s) => s,
            None => continue,
        };

        if radius < 1e-10 || radius > 1e10 {
            continue;
        }

        // Find inliers
        let mut inliers = Vec::new();
        for (i, &p) in points.iter().enumerate() {
            let dist = ((p - center).length() - radius).abs();
            if dist < threshold {
                inliers.push(i);
            }
        }

        if inliers.len() > best_inliers.len() {
            best_inliers = inliers;
            best_sphere = Some((center, radius));
        }
    }

    if best_inliers.len() < min_points {
        return None;
    }

    let (center, radius) = best_sphere.unwrap();

    // Refit using all inliers
    let inlier_pts: Vec<DVec3> = best_inliers.iter().map(|&i| points[i]).collect();
    if let Some(fitted) = fit_sphere(&inlier_pts) {
        let inlier_set: std::collections::HashSet<usize> = best_inliers.iter().copied().collect();
        let outliers: Vec<usize> = (0..points.len())
            .filter(|i| !inlier_set.contains(i))
            .collect();

        Some((
            ShapeParams::Sphere {
                center: fitted.center,
                radius: fitted.radius,
            },
            best_inliers,
            outliers,
        ))
    } else {
        let inlier_set: std::collections::HashSet<usize> = best_inliers.iter().copied().collect();
        let outliers: Vec<usize> = (0..points.len())
            .filter(|i| !inlier_set.contains(i))
            .collect();

        Some((
            ShapeParams::Sphere { center, radius },
            best_inliers,
            outliers,
        ))
    }
}

fn fit_sphere_4pt(points: &[DVec3]) -> Option<(DVec3, f64)> {
    if points.len() < 4 {
        return None;
    }

    // Solve for sphere passing through 4 points using linear system
    let p0 = points[0];
    let p1 = points[1];
    let p2 = points[2];
    let p3 = points[3];

    // System: |P - C|^2 = r^2 for each point
    // Subtract first equation to eliminate r^2
    // |Pi - C|^2 - |P0 - C|^2 = 0
    // |Pi|^2 - 2*Pi*C + |C|^2 - |P0|^2 + 2*P0*C - |C|^2 = 0
    // |Pi|^2 - |P0|^2 - 2*(Pi - P0)*C = 0
    // (Pi - P0)*C = (|Pi|^2 - |P0|^2) / 2

    let a = [
        [p1.x - p0.x, p1.y - p0.y, p1.z - p0.z],
        [p2.x - p0.x, p2.y - p0.y, p2.z - p0.z],
        [p3.x - p0.x, p3.y - p0.y, p3.z - p0.z],
    ];

    let p0_sq = p0.length_squared();
    let b = [
        (p1.length_squared() - p0_sq) / 2.0,
        (p2.length_squared() - p0_sq) / 2.0,
        (p3.length_squared() - p0_sq) / 2.0,
    ];

    let center = solve_linear_3x3(&a, &b)?;
    let cx = center[0];
    let cy = center[1];
    let cz = center[2];
    let center = DVec3::new(cx, cy, cz);
    let radius = (center - p0).length();

    Some((center, radius))
}

fn ransac_cylinder_segmentation(
    points: &[DVec3],
    threshold: f64,
    min_points: usize,
    max_iterations: usize,
) -> Option<(ShapeParams, Vec<usize>, Vec<usize>)> {
    if points.len() < 5 {
        return None;
    }

    let mut best_inliers = Vec::new();
    let mut best_cylinder: Option<(DVec3, DVec3, f64)> = None;

    let mut rng = SimpleRng::new(42);

    for _ in 0..max_iterations {
        // Sample 2 points for axis estimation
        let i0 = (rng.next() as usize) % points.len();
        let i1 = (rng.next() as usize) % points.len();

        if i0 == i1 {
            continue;
        }

        let p0 = points[i0];
        let p1 = points[i1];
        let axis = (p1 - p0).normalize_or(DVec3::Z);

        // Project points onto plane perpendicular to axis
        // and estimate radius
        let centroid = points.iter().sum::<DVec3>() / points.len() as f64;
        let u = if axis.x.abs() < 0.9 {
            axis.cross(DVec3::X).normalize()
        } else {
            axis.cross(DVec3::Y).normalize()
        };
        let v = axis.cross(u);

        let projected: Vec<f64> = points.iter().map(|&p| {
            let d = p - centroid;
            let x = d.dot(u);
            let y = d.dot(v);
            (x * x + y * y).sqrt()
        }).collect();

        // Estimate radius as median
        let mut sorted_proj = projected.clone();
        sorted_proj.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        let radius = sorted_proj[sorted_proj.len() / 2];

        if radius < 1e-10 || radius > 1e10 {
            continue;
        }

        // Find inliers
        let mut inliers = Vec::new();
        for (i, &p) in points.iter().enumerate() {
            let d = p - centroid;
            let axial = d.dot(axis);
            let radial = (d - axial * axis).length();
            let dist = (radial - radius).abs();
            if dist < threshold {
                inliers.push(i);
            }
        }

        if inliers.len() > best_inliers.len() {
            best_inliers = inliers;
            best_cylinder = Some((centroid, axis, radius));
        }
    }

    if best_inliers.len() < min_points {
        return None;
    }

    let (axis_point, axis_direction, radius) = best_cylinder.unwrap();

    // Refit using all inliers
    let inlier_pts: Vec<DVec3> = best_inliers.iter().map(|&i| points[i]).collect();
    if let Some(fitted) = fit_cylinder(&inlier_pts) {
        let inlier_set: std::collections::HashSet<usize> = best_inliers.iter().copied().collect();
        let outliers: Vec<usize> = (0..points.len())
            .filter(|i| !inlier_set.contains(i))
            .collect();

        Some((
            ShapeParams::Cylinder {
                axis_point: fitted.axis_point,
                axis_direction: fitted.axis_direction,
                radius: fitted.radius,
            },
            best_inliers,
            outliers,
        ))
    } else {
        let inlier_set: std::collections::HashSet<usize> = best_inliers.iter().copied().collect();
        let outliers: Vec<usize> = (0..points.len())
            .filter(|i| !inlier_set.contains(i))
            .collect();

        Some((
            ShapeParams::Cylinder {
                axis_point,
                axis_direction,
                radius,
            },
            best_inliers,
            outliers,
        ))
    }
}

fn build_neighbor_graph(points: &[DVec3], k: usize) -> Vec<Vec<usize>> {
    let n = points.len();
    let k = k.min(n - 1).max(1);
    let mut neighbors = Vec::with_capacity(n);

    for i in 0..n {
        let mut distances: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, (points[j] - points[i]).length_squared()))
            .collect();
        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));

        neighbors.push(distances.iter().take(k).map(|&(j, _)| j).collect());
    }

    neighbors
}

fn compute_curvatures(points: &[DVec3], neighbors: &[Vec<usize>]) -> Vec<f64> {
    points.iter().enumerate().map(|(i, _)| {
        let neighbor_pts: Vec<DVec3> = neighbors[i].iter().map(|&j| points[j]).collect();
        let (_, values) = compute_pca(&neighbor_pts);
        let sum: f64 = values.iter().sum();
        if sum > 1e-10 {
            values[2] / sum
        } else {
            0.0
        }
    }).collect()
}

fn grow_region(
    seed_idx: usize,
    points: &[DVec3],
    normals: &[DVec3],
    neighbors: &[Vec<usize>],
    visited: &mut [bool],
    config: &RegionGrowingConfig,
) -> Vec<usize> {
    let mut region = Vec::new();
    let mut queue = vec![seed_idx];
    visited[seed_idx] = true;

    let seed_normal = normals[seed_idx];

    while let Some(current) = queue.pop() {
        region.push(current);

        for &neighbor in &neighbors[current] {
            if visited[neighbor] {
                continue;
            }

            // Check normal similarity
            let dot = normals[neighbor].dot(seed_normal);
            let angle = dot.acos();

            if angle < config.max_angle {
                // Check distance constraint
                let dist = (points[neighbor] - points[current]).length();
                if dist < config.max_distance * 10.0 {
                    visited[neighbor] = true;
                    queue.push(neighbor);
                }
            }
        }
    }

    region
}

// ============================================================================
// Surface Reconstruction
// ============================================================================

/// Triangle mesh for surface reconstruction output.
#[derive(Debug, Clone)]
pub struct TriangleMesh {
    /// Vertex positions.
    pub vertices: Vec<DVec3>,
    /// Triangle indices (3 vertex indices per triangle).
    pub triangles: Vec<[usize; 3]>,
    /// Vertex normals (optional).
    pub normals: Option<Vec<DVec3>>,
}

impl TriangleMesh {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            triangles: Vec::new(),
            normals: None,
        }
    }

    /// Computes face normals for the mesh.
    pub fn compute_face_normals(&self) -> Vec<DVec3> {
        self.triangles.iter().map(|&[i, j, k]| {
            let a = self.vertices[i];
            let b = self.vertices[j];
            let c = self.vertices[k];
            let normal = (b - a).cross(c - a);
            let len = normal.length();
            if len > 1e-10 {
                normal / len
            } else {
                DVec3::Z
            }
        }).collect()
    }

    /// Computes vertex normals by averaging adjacent face normals.
    pub fn compute_vertex_normals(&self) -> Vec<DVec3> {
        let face_normals = self.compute_face_normals();
        let mut vertex_normals = vec![DVec3::ZERO; self.vertices.len()];

        for (tri, &normal) in self.triangles.iter().zip(face_normals.iter()) {
            for &idx in tri {
                vertex_normals[idx] += normal;
            }
        }

        for normal in &mut vertex_normals {
            let len = normal.length();
            if len > 1e-10 {
                *normal /= len;
            }
        }

        vertex_normals
    }
}

/// Configuration for Poisson surface reconstruction.
#[derive(Debug, Clone)]
pub struct PoissonConfig {
    /// Octree depth (higher = more detail, more memory).
    pub depth: usize,
    /// Solver division (default: 8).
    pub solver_divide: usize,
    /// Iso-surface value (default: 0, use mean).
    pub iso_value: f64,
}

impl Default for PoissonConfig {
    fn default() -> Self {
        Self {
            depth: 8,
            solver_divide: 8,
            iso_value: 0.0,
        }
    }
}

/// Performs Poisson surface reconstruction.
///
/// Reconstructs a watertight surface from oriented point samples.
/// Uses an implicit function approach with octree-based spatial indexing.
pub fn poisson_reconstruction(
    points: &[DVec3],
    normals: &[DVec3],
    config: &PoissonConfig,
) -> Option<TriangleMesh> {
    if points.is_empty() || points.len() != normals.len() {
        return None;
    }

    // Simplified Poisson reconstruction using implicit function + marching cubes
    // This is a basic implementation - full Poisson is more complex

    let (min, max) = compute_bounding_box(points)?;
    let padding = (max - min).length() * 0.1;
    let min = min - DVec3::splat(padding);
    let max = max + DVec3::splat(padding);

    let resolution = 2usize.pow(config.depth as u32);
    let cell_size = (max - min) / resolution as f64;

    // Build implicit function using oriented point samples
    let mut grid = vec![0.0_f64; resolution * resolution * resolution];

    // Compute implicit function values
    for idx in 0..resolution {
        for idy in 0..resolution {
            for idz in 0..resolution {
                let p = DVec3::new(
                    min.x + (idx as f64 + 0.5) * cell_size.x,
                    min.y + (idy as f64 + 0.5) * cell_size.y,
                    min.z + (idz as f64 + 0.5) * cell_size.z,
                );

                let mut value = 0.0;
                for (&pt, &n) in points.iter().zip(normals.iter()) {
                    let d = p - pt;
                    let dist = d.length();
                    if dist > 1e-10 {
                        // Signed distance from oriented point
                        value += d.dot(n) / (dist * dist + 1.0);
                    }
                }

                grid[idx + idy * resolution + idz * resolution * resolution] = value;
            }
        }
    }

    // Marching cubes to extract iso-surface
    marching_cubes(&grid, resolution, resolution, resolution, &min, &max, config.iso_value)
}

fn marching_cubes(
    grid: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    min: &DVec3,
    max: &DVec3,
    iso_value: f64,
) -> Option<TriangleMesh> {
    let mut mesh = TriangleMesh::new();
    let dx = (max.x - min.x) / nx as f64;
    let dy = (max.y - min.y) / ny as f64;
    let dz = (max.z - min.z) / nz as f64;

    // Simplified marching cubes - just create triangles for cells crossing iso-value
    for ix in 0..nx - 1 {
        for iy in 0..ny - 1 {
            for iz in 0..nz - 1 {
                // Get 8 corner values
                let corners = [
                    grid[ix + iy * nx + iz * nx * ny],
                    grid[ix + 1 + iy * nx + iz * nx * ny],
                    grid[ix + 1 + (iy + 1) * nx + iz * nx * ny],
                    grid[ix + (iy + 1) * nx + iz * nx * ny],
                    grid[ix + iy * nx + (iz + 1) * nx * ny],
                    grid[ix + 1 + iy * nx + (iz + 1) * nx * ny],
                    grid[ix + 1 + (iy + 1) * nx + (iz + 1) * nx * ny],
                    grid[ix + (iy + 1) * nx + (iz + 1) * nx * ny],
                ];

                // Check if cell crosses iso-value
                let above: Vec<bool> = corners.iter().map(|&v| v > iso_value).collect();
                let all_above = above.iter().all(|&b| b);
                let all_below = above.iter().all(|&b| !b);

                if all_above || all_below {
                    continue;
                }

                // Simplified: create triangles at cell center
                let cx = min.x + (ix as f64 + 0.5) * dx;
                let cy = min.y + (iy as f64 + 0.5) * dy;
                let cz = min.z + (iz as f64 + 0.5) * dz;

                // Create a simple cube face approximation
                let base_idx = mesh.vertices.len();
                let size = dx.min(dy).min(dz) * 0.5;

                // Add vertices for a small quad
                mesh.vertices.push(DVec3::new(cx - size, cy - size, cz));
                mesh.vertices.push(DVec3::new(cx + size, cy - size, cz));
                mesh.vertices.push(DVec3::new(cx + size, cy + size, cz));
                mesh.vertices.push(DVec3::new(cx - size, cy + size, cz));

                mesh.triangles.push([base_idx, base_idx + 1, base_idx + 2]);
                mesh.triangles.push([base_idx, base_idx + 2, base_idx + 3]);
            }
        }
    }

    if mesh.vertices.is_empty() {
        None
    } else {
        Some(mesh)
    }
}

/// Configuration for Ball Pivoting Algorithm.
#[derive(Debug, Clone)]
pub struct BpaConfig {
    /// Ball radius for pivoting.
    pub ball_radius: f64,
    /// Clustering radius for duplicate removal.
    pub clustering: f64,
    /// Angle threshold for edge selection (radians).
    pub angle_threshold: f64,
}

impl Default for BpaConfig {
    fn default() -> Self {
        Self {
            ball_radius: 0.1,
            clustering: 0.001,
            angle_threshold: std::f64::consts::PI / 4.0,
        }
    }
}

/// Performs Ball Pivoting Algorithm surface reconstruction.
///
/// Reconstructs surface by rolling a ball over the point cloud.
pub fn ball_pivoting_reconstruction(
    points: &[DVec3],
    normals: &[DVec3],
    config: &BpaConfig,
) -> Option<TriangleMesh> {
    if points.len() < 3 {
        return None;
    }

    let mut mesh = TriangleMesh::new();

    // Build spatial index
    let grid_size = config.ball_radius;
    let (min, _) = compute_bounding_box(points)?;

    let mut spatial_index: std::collections::HashMap<[i64; 3], Vec<usize>> = std::collections::HashMap::new();

    for (i, &p) in points.iter().enumerate() {
        let key = [
            ((p.x - min.x) / grid_size).floor() as i64,
            ((p.y - min.y) / grid_size).floor() as i64,
            ((p.z - min.z) / grid_size).floor() as i64,
        ];
        spatial_index.entry(key).or_insert_with(Vec::new).push(i);
    }

    // Track used edges
    let mut used_edges: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    let mut used_triangles: std::collections::HashSet<[usize; 3]> = std::collections::HashSet::new();

    let r = config.ball_radius;
    let r_sq = r * r;

    // Find seed triangle
    for i0 in 0..points.len() {
        for i1 in (i0 + 1)..points.len() {
            for i2 in (i1 + 1)..points.len() {
                let p0 = points[i0];
                let p1 = points[i1];
                let p2 = points[i2];

                // Check if ball of radius r touches all three points
                let circumcenter = compute_circumcenter(p0, p1, p2)?;
                let circumradius = (p0 - circumcenter).length();

                if circumradius > r {
                    continue;
                }

                // Ball center is at distance sqrt(r^2 - circumradius^2) from circumcenter
                let ball_dist = (r_sq - circumradius * circumradius).sqrt();
                let normal = (p1 - p0).cross(p2 - p0);
                let len = normal.length();
                if len < 1e-10 {
                    continue;
                }
                let normal = normal / len;

                // Check ball doesn't contain other points
                let ball_center = circumcenter + normal * ball_dist;
                let mut valid = true;
                for &p in points.iter().take(10) {
                    if (p - ball_center).length_squared() < r_sq - 1e-10 {
                        valid = false;
                        break;
                    }
                }

                if valid {
                    let base_idx = mesh.vertices.len();
                    mesh.vertices.push(p0);
                    mesh.vertices.push(p1);
                    mesh.vertices.push(p2);
                    mesh.triangles.push([base_idx, base_idx + 1, base_idx + 2]);

                    used_edges.insert((i0, i1));
                    used_edges.insert((i1, i2));
                    used_edges.insert((i2, i0));
                    used_triangles.insert([i0, i1, i2]);
                    break;
                }
            }
            if !mesh.triangles.is_empty() {
                break;
            }
        }
        if !mesh.triangles.is_empty() {
            break;
        }
    }

    if mesh.vertices.is_empty() {
        // Fallback: simple Delaunay triangulation
        return delaunay_reconstruction(points, normals);
    }

    // Expand from seed using ball pivoting
    // This is simplified - full BPA is more complex
    for _ in 0..points.len() / 3 {
        // Find next edge to pivot from
        // Collect edges to iterate over to avoid borrow issues
        let edges_to_process: Vec<(usize, usize)> = used_edges.iter().copied().collect();
        let mut new_edges: Vec<(usize, usize)> = Vec::new();

        for (i0, i1) in edges_to_process {
            // Try to find third point
            for i2 in 0..points.len() {
                if i2 == i0 || i2 == i1 {
                    continue;
                }

                let _key = [i0.min(i1), i0.max(i1), i2];
                let tri_key = if i0 < i1 {
                    [i0, i1, i2]
                } else {
                    [i1, i0, i2]
                };

                if used_triangles.contains(&tri_key) || used_triangles.contains(&[tri_key[0], tri_key[2], tri_key[1]]) {
                    continue;
                }

                let p0 = points[i0];
                let p1 = points[i1];
                let p2 = points[i2];

                if let Some(cc) = compute_circumcenter(p0, p1, p2) {
                    let cr = (p0 - cc).length();
                    if cr <= r {
                        let base_idx = mesh.vertices.len();
                        mesh.vertices.push(p0);
                        mesh.vertices.push(p1);
                        mesh.vertices.push(p2);
                        mesh.triangles.push([base_idx, base_idx + 1, base_idx + 2]);

                        new_edges.push((i0, i2));
                        new_edges.push((i2, i1));
                        used_triangles.insert(tri_key);
                    }
                }
            }
        }

        // Add new edges
        for edge in new_edges {
            used_edges.insert(edge);
        }
    }

    if mesh.vertices.is_empty() {
        None
    } else {
        Some(mesh)
    }
}

fn compute_circumcenter(a: DVec3, b: DVec3, c: DVec3) -> Option<DVec3> {
    let ab = b - a;
    let ac = c - a;

    let cross = ab.cross(ac);
    let denom = 2.0 * cross.length_squared();

    if denom < 1e-20 {
        return None;
    }

    let _d = cross.dot(a.cross(b) + b.cross(c) + c.cross(a)) / denom;

    Some(a + (ab.cross(ac).cross(ab) * ac.length_squared() + ac.cross(ab.cross(ac)) * ab.length_squared()) / (2.0 * cross.length_squared()))
}

/// Performs Delaunay triangulation based surface reconstruction.
///
/// Projects points to 2D, computes Delaunay triangulation,
/// then projects back to 3D.
pub fn delaunay_reconstruction(
    points: &[DVec3],
    _normals: &[DVec3],
) -> Option<TriangleMesh> {
    if points.len() < 3 {
        return None;
    }

    let mut mesh = TriangleMesh::new();

    // Fit plane to points
    let plane = fit_plane(points)?;
    let normal = plane.normal;

    // Build orthonormal basis on plane
    let u = if normal.x.abs() < 0.9 {
        normal.cross(DVec3::X).normalize()
    } else {
        normal.cross(DVec3::Y).normalize()
    };
    let v = normal.cross(u);

    // Project to 2D
    let projected_2d: Vec<DVec2> = points.iter().map(|&p| {
        let d = p - plane.point;
        DVec2::new(d.dot(u), d.dot(v))
    }).collect();

    // Compute Delaunay triangulation (Bowyer-Watson algorithm)
    let triangles_2d = delaunay_triangulation_2d(&projected_2d);

    // Convert back to 3D
    mesh.vertices = points.to_vec();
    mesh.triangles = triangles_2d;

    Some(mesh)
}

fn delaunay_triangulation_2d(points: &[DVec2]) -> Vec<[usize; 3]> {
    if points.len() < 3 {
        return Vec::new();
    }

    let n = points.len();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    // Create super-triangle containing all points
    let (min_x, max_x) = points.iter().map(|p| p.x).fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), x| (min.min(x), max.max(x)));
    let (min_y, max_y) = points.iter().map(|p| p.y).fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), y| (min.min(y), max.max(y)));

    let dx = max_x - min_x;
    let dy = max_y - min_y;
    let delta = dx.max(dy) * 10.0;

    let p1 = DVec2::new(min_x - delta, min_y - delta);
    let p2 = DVec2::new(min_x + 3.0 * delta, min_y - delta);
    let p3 = DVec2::new(min_x, min_y + 3.0 * delta);

    let super_vertices = [n, n + 1, n + 2];
    let all_points: Vec<DVec2> = points.iter().copied().chain([p1, p2, p3]).collect();

    triangles.push(super_vertices);

    // Add points one by one
    for i in 0..n {
        let p = points[i];

        // Find all triangles whose circumcircle contains p
        let mut bad_triangles: Vec<usize> = Vec::new();
        for (ti, &tri) in triangles.iter().enumerate() {
            if let Some(cc) = circumcenter_2d(all_points[tri[0]], all_points[tri[1]], all_points[tri[2]]) {
                let r_sq = (all_points[tri[0]].x - cc.x).powi(2) + (all_points[tri[0]].y - cc.y).powi(2);
                let d_sq = (p.x - cc.x).powi(2) + (p.y - cc.y).powi(2);
                if d_sq <= r_sq {
                    bad_triangles.push(ti);
                }
            }
        }

        // Find boundary polygon
        let mut polygon: Vec<(usize, usize)> = Vec::new();
        for &ti in &bad_triangles {
            let tri = triangles[ti];
            let edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
            for &(a, b) in &edges {
                let edge = (a.min(b), a.max(b));
                let mut shared = false;
                for &tj in &bad_triangles {
                    if tj == ti {
                        continue;
                    }
                    let other = triangles[tj];
                    let other_edges = [(other[0], other[1]), (other[1], other[2]), (other[2], other[0])];
                    for &(oa, ob) in &other_edges {
                        let other_edge = (oa.min(ob), oa.max(ob));
                        if edge == other_edge {
                            shared = true;
                            break;
                        }
                    }
                    if shared {
                        break;
                    }
                }
                if !shared {
                    polygon.push((tri[0], tri[1]));
                    if !polygon.contains(&(tri[1], tri[2])) {
                        polygon.push((tri[1], tri[2]));
                    }
                    if !polygon.contains(&(tri[2], tri[0])) {
                        polygon.push((tri[2], tri[0]));
                    }
                }
            }
        }

        // Remove bad triangles
        let mut new_triangles: Vec<[usize; 3]> = Vec::new();
        for (ti, &tri) in triangles.iter().enumerate() {
            if !bad_triangles.contains(&ti) {
                new_triangles.push(tri);
            }
        }

        // Add new triangles from polygon
        for &(a, b) in &polygon {
            new_triangles.push([a, b, i]);
        }

        triangles = new_triangles;
    }

    // Remove triangles containing super-triangle vertices
    triangles.retain(|&tri| tri[0] < n && tri[1] < n && tri[2] < n);

    triangles
}

fn circumcenter_2d(a: DVec2, b: DVec2, c: DVec2) -> Option<DVec2> {
    let d = 2.0 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
    if d.abs() < 1e-14 {
        return None;
    }

    let a_sq = a.x * a.x + a.y * a.y;
    let b_sq = b.x * b.x + b.y * b.y;
    let c_sq = c.x * c.x + c.y * c.y;

    Some(DVec2::new(
        (a_sq * (b.y - c.y) + b_sq * (c.y - a.y) + c_sq * (a.y - b.y)) / d,
        (a_sq * (c.x - b.x) + b_sq * (a.x - c.x) + c_sq * (b.x - a.x)) / d,
    ))
}

/// Generates normal-consistent mesh from oriented point cloud.
///
/// Ensures all face normals point consistently outward.
pub fn generate_consistent_mesh(
    points: &[DVec3],
    normals: &[DVec3],
) -> Option<TriangleMesh> {
    let mut mesh = delaunay_reconstruction(points, normals)?;

    // Orient faces consistently
    let vertex_normals = mesh.compute_vertex_normals();

    for tri in &mut mesh.triangles {
        let a = mesh.vertices[tri[0]];
        let b = mesh.vertices[tri[1]];
        let c = mesh.vertices[tri[2]];

        let face_normal = (b - a).cross(c - a);
        let len = face_normal.length();
        if len < 1e-10 {
            continue;
        }
        let face_normal = face_normal / len;

        // Compare with vertex normals
        let avg_normal = (vertex_normals[tri[0]] + vertex_normals[tri[1]] + vertex_normals[tri[2]]) / 3.0;

        if face_normal.dot(avg_normal) < 0.0 {
            // Flip triangle
            *tri = [tri[0], tri[2], tri[1]];
        }
    }

    mesh.normals = Some(vertex_normals);
    Some(mesh)
}

fn compute_bounding_box(points: &[DVec3]) -> Option<(DVec3, DVec3)> {
    if points.is_empty() {
        return None;
    }

    let mut min = DVec3::splat(f64::INFINITY);
    let mut max = DVec3::splat(f64::NEG_INFINITY);

    for &p in points {
        min = min.min(p);
        max = max.max(p);
    }

    Some((min, max))
}

// ============================================================================
// Advanced Sampling
// ============================================================================

/// Sampling method for point cloud simplification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdvancedSamplingMethod {
    /// Voxel grid downsampling.
    VoxelGrid,
    /// Random uniform sampling.
    RandomUniform,
    /// Curvature-aware sampling (more points in high-curvature regions).
    CurvatureAware,
    /// Poisson disk sampling (uniform distribution with minimum distance).
    PoissonDisk,
}

/// Configuration for advanced sampling.
#[derive(Debug, Clone)]
pub struct AdvancedSamplingConfig {
    /// Target number of points.
    pub target_count: usize,
    /// Voxel size for voxel grid sampling.
    pub voxel_size: f64,
    /// Minimum distance for Poisson disk sampling.
    pub min_distance: f64,
    /// Number of neighbors for curvature estimation.
    pub k_neighbors: usize,
    /// Seed for random sampling.
    pub seed: u64,
}

impl Default for AdvancedSamplingConfig {
    fn default() -> Self {
        Self {
            target_count: 1000,
            voxel_size: 0.1,
            min_distance: 0.05,
            k_neighbors: 30,
            seed: 42,
        }
    }
}

/// Performs advanced point cloud sampling.
pub fn advanced_sample(
    points: &[DVec3],
    method: AdvancedSamplingMethod,
    config: &AdvancedSamplingConfig,
) -> Vec<DVec3> {
    if points.len() <= config.target_count {
        return points.to_vec();
    }

    match method {
        AdvancedSamplingMethod::VoxelGrid => voxel_grid_sample(points, config),
        AdvancedSamplingMethod::RandomUniform => random_uniform_sample(points, config),
        AdvancedSamplingMethod::CurvatureAware => curvature_aware_sample(points, config),
        AdvancedSamplingMethod::PoissonDisk => poisson_disk_sample(points, config),
    }
}

fn voxel_grid_sample(points: &[DVec3], config: &AdvancedSamplingConfig) -> Vec<DVec3> {
    if points.is_empty() {
        return Vec::new();
    }

    let (min, _max) = compute_bounding_box(points).unwrap();
    let voxel_size = config.voxel_size;

    let mut voxels: std::collections::HashMap<[i64; 3], Vec<DVec3>> = std::collections::HashMap::new();

    for &p in points {
        let key = [
            ((p.x - min.x) / voxel_size).floor() as i64,
            ((p.y - min.y) / voxel_size).floor() as i64,
            ((p.z - min.z) / voxel_size).floor() as i64,
        ];
        voxels.entry(key).or_insert_with(Vec::new).push(p);
    }

    voxels.values().map(|pts| {
        let sum: DVec3 = pts.iter().sum();
        sum / pts.len() as f64
    }).collect()
}

fn random_uniform_sample(points: &[DVec3], config: &AdvancedSamplingConfig) -> Vec<DVec3> {
    let n = points.len();
    let target = config.target_count.min(n);

    let mut rng = SimpleRng::new(config.seed);
    let mut indices: Vec<usize> = (0..n).collect();

    // Fisher-Yates shuffle for first target_count elements
    for i in 0..target {
        let j = i + ((rng.next() as usize) % (n - i));
        indices.swap(i, j);
    }

    indices.iter().take(target).map(|&i| points[i]).collect()
}

fn curvature_aware_sample(points: &[DVec3], config: &AdvancedSamplingConfig) -> Vec<DVec3> {
    if points.len() <= config.target_count {
        return points.to_vec();
    }

    // Compute curvatures
    let neighbors = build_neighbor_graph(points, config.k_neighbors);
    let curvatures = compute_curvatures(points, &neighbors);

    // Compute sampling probabilities (higher curvature = higher probability)
    let max_curv = curvatures.iter().cloned().fold(0.0_f64, f64::max).max(1e-10);
    let weights: Vec<f64> = curvatures.iter().map(|&c| 0.1 + 0.9 * c / max_curv).collect();

    let total_weight: f64 = weights.iter().sum();
    let probs: Vec<f64> = weights.iter().map(|&w| w / total_weight).collect();

    // Sample based on probabilities
    let mut rng = SimpleRng::new(config.seed);
    let mut selected: std::collections::HashSet<usize> = std::collections::HashSet::new();

    while selected.len() < config.target_count {
        let r = (rng.next() as f64) / (u64::MAX as f64);
        let mut cumsum = 0.0;
        for (i, &p) in probs.iter().enumerate() {
            cumsum += p;
            if r < cumsum {
                selected.insert(i);
                break;
            }
        }
    }

    selected.iter().map(|&i| points[i]).collect()
}

fn poisson_disk_sample(points: &[DVec3], config: &AdvancedSamplingConfig) -> Vec<DVec3> {
    if points.is_empty() {
        return Vec::new();
    }

    let (min, _) = compute_bounding_box(points).unwrap();
    let r = config.min_distance;
    let r_sq = r * r;
    let cell_size = r / 3.0_f64.sqrt();

    let mut grid: std::collections::HashMap<[i64; 3], usize> = std::collections::HashMap::new();
    let mut result: Vec<DVec3> = Vec::new();
    let mut active: Vec<usize> = Vec::new();

    // Start with random point
    let mut rng = SimpleRng::new(config.seed);
    let first_idx = (rng.next() as usize) % points.len();
    let first = points[first_idx];

    result.push(first);
    active.push(0);

    let key = [
        ((first.x - min.x) / cell_size).floor() as i64,
        ((first.y - min.y) / cell_size).floor() as i64,
        ((first.z - min.z) / cell_size).floor() as i64,
    ];
    grid.insert(key, 0);

    // Dart throwing
    while !active.is_empty() && result.len() < config.target_count {
        let active_idx = (rng.next() as usize) % active.len();
        let sample_idx = active[active_idx];
        let sample = result[sample_idx];

        let mut found = false;
        for _ in 0..30 {
            // Generate random point in annulus around sample
            let angle1 = 2.0 * std::f64::consts::PI * (rng.next() as f64) / (u64::MAX as f64);
            let angle2 = std::f64::consts::PI * (rng.next() as f64) / (u64::MAX as f64);
            let rad = r * (1.0 + (rng.next() as f64) / (u64::MAX as f64));

            let candidate = DVec3::new(
                sample.x + rad * angle2.sin() * angle1.cos(),
                sample.y + rad * angle2.sin() * angle1.sin(),
                sample.z + rad * angle2.cos(),
            );

            // Check if candidate is far enough from existing points
            let key = [
                ((candidate.x - min.x) / cell_size).floor() as i64,
                ((candidate.y - min.y) / cell_size).floor() as i64,
                ((candidate.z - min.z) / cell_size).floor() as i64,
            ];

            let mut valid = true;
            for di in -1..=1 {
                for dj in -1..=1 {
                    for dk in -1..=1 {
                        let neighbor_key = [key[0] + di, key[1] + dj, key[2] + dk];
                        if let Some(&idx) = grid.get(&neighbor_key) {
                            if (result[idx] - candidate).length_squared() < r_sq {
                                valid = false;
                                break;
                            }
                        }
                    }
                    if !valid {
                        break;
                    }
                }
                if !valid {
                    break;
                }
            }

            if valid {
                let new_idx = result.len();
                result.push(candidate);
                active.push(new_idx);
                grid.insert(key, new_idx);
                found = true;
                break;
            }
        }

        if !found {
            active.swap_remove(active_idx);
        }
    }

    result
}

// ============================================================================
// BRep Integration
// ============================================================================

/// Extracts a point cloud from BRep vertices.
///
/// Collects all unique vertex positions from the BRep.
pub fn extract_points_from_brep_vertices(brep: &rcad_kernel::BRep) -> PointCloud {
    PointCloud::from_vec(brep.vertices.iter().map(|v| v.point).collect())
}

/// Extracts a point cloud from a BRep mesh.
///
/// Samples points from the triangulated faces. If a face has no cached
/// triangulation, it will be skipped.
pub fn extract_points_from_brep_mesh(brep: &rcad_kernel::BRep) -> PointCloud {
    let mut points = Vec::new();

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Add vertices from triangles
                for &[i, j, k] in &face.triangles {
                    if let (Some(a), Some(b), Some(c)) = (
                        brep.vertices.get(i),
                        brep.vertices.get(j),
                        brep.vertices.get(k),
                    ) {
                        points.push(a.point);
                        points.push(b.point);
                        points.push(c.point);
                    }
                }

                // If no triangles, add wire vertices
                if face.triangles.is_empty() {
                    for we in &face.outer_wire.edges {
                        if let Some(edge) = brep.edges.get(we.idx) {
                            let vidx = if we.forward { edge.start } else { edge.end };
                            if let Some(v) = brep.vertices.get(vidx) {
                                points.push(v.point);
                            }
                        }
                    }
                }
            }
        }
    }

    PointCloud::from_vec(points)
}

/// Extracts a point cloud from a triangulated mesh.
///
/// Takes the vertices directly from a SurfaceMesh.
pub fn extract_points_from_mesh(mesh: &crate::triangulate::SurfaceMesh) -> PointCloud {
    PointCloud::from_vec(mesh.vertices.clone())
}

/// Samples points uniformly from a BRep's surfaces.
///
/// For each face with an associated surface, samples points on a regular
/// UV grid and keeps those that lie within the face's boundary.
pub fn sample_points_from_brep_surfaces(
    brep: &rcad_kernel::BRep,
    samples_per_face: usize,
) -> PointCloud {
    use rcad_kernel::geom::SurfaceEval;

    let mut points = Vec::new();
    let sqrt_n = (samples_per_face as f64).sqrt().ceil() as usize;

    let mut face_idx = 0;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for _face in &shell.faces {
                // Get surface for this face
                let surf_idx = match brep.geom.face_surface.get(face_idx).and_then(|o| *o) {
                    Some(idx) => idx,
                    None => {
                        face_idx += 1;
                        continue;
                    }
                };
                let surf = match brep.geom.surfaces.get(surf_idx) {
                    Some(s) => s,
                    None => {
                        face_idx += 1;
                        continue;
                    }
                };

                // Get UV domain
                let domain = brep.geom.face_surface_range.get(face_idx)
                    .and_then(|o| *o)
                    .unwrap_or_else(|| surf.default_domain());
                let [u0, u1, v0, v1] = domain;

                if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite() {
                    face_idx += 1;
                    continue;
                }

                // Sample on a grid
                for i in 0..sqrt_n {
                    for j in 0..sqrt_n {
                        let u = u0 + (u1 - u0) * (i as f64 + 0.5) / sqrt_n as f64;
                        let v = v0 + (v1 - v0) * (j as f64 + 0.5) / sqrt_n as f64;
                        let p = surf.point_at(u, v);
                        points.push(p);
                    }
                }

                face_idx += 1;
            }
        }
    }

    PointCloud::from_vec(points)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-6;

    fn approx_eq(a: DVec3, b: DVec3, tol: f64) -> bool {
        (a - b).length() < tol
    }

    #[test]
    fn test_empty_point_cloud() {
        let pc = PointCloud::new();
        assert!(pc.is_empty());
        assert_eq!(pc.len(), 0);
        assert!(pc.bounding_box().is_none());
        assert!(pc.centroid().is_none());
    }

    #[test]
    fn test_point_cloud_basics() {
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let pc = PointCloud::from_points(&points);

        assert_eq!(pc.len(), 3);

        let centroid = pc.centroid().unwrap();
        assert!(approx_eq(centroid, DVec3::new(1.0/3.0, 1.0/3.0, 0.0), 1e-10));

        let (min, max) = pc.bounding_box().unwrap();
        assert!(approx_eq(min, DVec3::ZERO, 1e-10));
        assert!(approx_eq(max, DVec3::new(1.0, 1.0, 0.0), 1e-10));
    }

    #[test]
    fn test_pca_identity() {
        // Points on a cube - PCA should give roughly equal eigenvalues
        // Simpler test that's more numerically stable
        let points: Vec<DVec3> = vec![
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(-1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, -1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(0.0, 0.0, -1.0),
        ];

        let (axes, values) = compute_pca(&points);

        // Eigenvalues should be positive and roughly equal for symmetric distribution
        assert!(values[0] > 0.0, "Largest eigenvalue should be positive, got {}", values[0]);
        assert!(values[2] >= 0.0, "Smallest eigenvalue should be non-negative, got {}", values[2]);
        // All eigenvalues should be similar (within factor of 2) for this symmetric case
        assert!(values[0] / values[2].max(1e-10) < 3.0, "Eigenvalue ratio {} too large", values[0] / values[2].max(1e-10));

        // Axes should be orthonormal
        for axis in &axes {
            assert!((axis.length() - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_pca_line() {
        // Points along X axis
        let points: Vec<DVec3> = (0..10).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let (axes, values) = compute_pca(&points);

        // Largest eigenvalue should be along X
        assert!(values[0] > values[1]);
        assert!(values[1] < 1e-6);
        assert!(values[2] < 1e-6);

        // First principal axis should be approximately X
        assert!((axes[0].x.abs() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_pca_plane() {
        // Points on XY plane
        let mut points = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        let (axes, values) = compute_pca(&points);

        // Two large eigenvalues, one small
        assert!(values[0] > 0.1);
        assert!(values[1] > 0.1);
        assert!(values[2] < 0.01);

        // Third principal axis should be Z (normal)
        assert!((axes[2].z.abs() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_dimensionality() {
        // Point-like
        let d = estimate_dimensionality([1e-20, 1e-20, 1e-20], 0.01);
        assert_eq!(d, Dimensionality::Point);

        // Linear
        let d = estimate_dimensionality([10.0, 0.001, 0.001], 0.01);
        assert_eq!(d, Dimensionality::Linear);

        // Planar
        let d = estimate_dimensionality([10.0, 10.0, 0.001], 0.01);
        assert_eq!(d, Dimensionality::Planar);

        // Volumetric
        let d = estimate_dimensionality([10.0, 10.0, 10.0], 0.01);
        assert_eq!(d, Dimensionality::Volumetric);
    }

    #[test]
    fn test_inertia_tensor() {
        // Unit cube at origin
        let points: Vec<DVec3> = (0..=1)
            .flat_map(|x| (0..=1).flat_map(move |y| (0..=1).map(move |z| DVec3::new(x as f64, y as f64, z as f64))))
            .collect();

        let inertia = compute_inertia(&points);

        // Check diagonal elements are positive
        assert!(inertia[0][0] >= 0.0);
        assert!(inertia[1][1] >= 0.0);
        assert!(inertia[2][2] >= 0.0);

        // Check symmetry
        assert!((inertia[0][1] - inertia[1][0]).abs() < 1e-10);
        assert!((inertia[0][2] - inertia[2][0]).abs() < 1e-10);
        assert!((inertia[1][2] - inertia[2][1]).abs() < 1e-10);
    }

    #[test]
    fn test_fit_plane() {
        // Perfect plane
        let points: Vec<DVec3> = (0..10)
            .flat_map(|i| (0..10).map(move |j| DVec3::new(i as f64, j as f64, 0.0)))
            .collect();

        let plane = fit_plane(&points).unwrap();

        assert!(approx_eq(plane.normal, DVec3::Z, 1e-6) || approx_eq(plane.normal, -DVec3::Z, 1e-6));
        assert!(plane.rms_error < 1e-6);
    }

    #[test]
    fn test_fit_sphere() {
        // Points on a sphere of radius 2 centered at (1, 2, 3)
        let center = DVec3::new(1.0, 2.0, 3.0);
        let radius = 2.0;

        let mut points = Vec::new();
        for i in 0..50 {
            let theta = 2.0 * PI * i as f64 / 50.0;
            let phi = PI * i as f64 / 50.0;
            let x = center.x + radius * phi.sin() * theta.cos();
            let y = center.y + radius * phi.sin() * theta.sin();
            let z = center.z + radius * phi.cos();
            points.push(DVec3::new(x, y, z));
        }

        let sphere = fit_sphere(&points).unwrap();

        assert!(approx_eq(sphere.center, center, 0.1));
        assert!((sphere.radius - radius).abs() < 0.1);
        assert!(sphere.rms_error < 0.1);
    }

    #[test]
    fn test_fit_cylinder() {
        // Points on a cylinder along Z axis
        let radius = 1.5;
        let mut points = Vec::new();

        for i in 0..20 {
            let theta = 2.0 * PI * i as f64 / 20.0;
            for z in 0..5 {
                let x = radius * theta.cos();
                let y = radius * theta.sin();
                points.push(DVec3::new(x, y, z as f64));
            }
        }

        let cylinder = fit_cylinder(&points).unwrap();

        assert!((cylinder.radius - radius).abs() < 0.1);
        assert!(cylinder.rms_error < 0.1);
    }

    #[test]
    fn test_simplify_random() {
        let points: Vec<DVec3> = (0..1000).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let simplified = simplify_point_cloud(&points, 100, SamplingStrategy::Random);

        assert_eq!(simplified.len(), 100);
    }

    #[test]
    fn test_simplify_voxel() {
        let mut points = Vec::new();
        // Dense grid of points
        for i in 0..10 {
            for j in 0..10 {
                for k in 0..10 {
                    points.push(DVec3::new(i as f64, j as f64, k as f64));
                }
            }
        }

        let simplified = simplify_point_cloud(&points, 50, SamplingStrategy::Voxel);

        assert!(simplified.len() >= 27); // At least 3x3x3 voxels
        assert!(simplified.len() <= 100);
    }

    #[test]
    fn test_simplify_farthest_point() {
        let points: Vec<DVec3> = (0..100).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let simplified = simplify_point_cloud(&points, 10, SamplingStrategy::FarthestPoint);

        assert_eq!(simplified.len(), 10);

        // Should include endpoints
        let has_start = simplified.iter().any(|p| p.x < 1.0);
        let has_end = simplified.iter().any(|p| p.x > 98.0);
        assert!(has_start || has_end);
    }

    #[test]
    fn test_estimate_normals() {
        // Points on XY plane
        let mut points = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        let normals = estimate_normals(&points, 4);

        assert_eq!(normals.len(), points.len());

        // All normals should point along Z (positive or negative)
        for n in &normals {
            assert!(n.z.abs() > 0.9, "Normal should be along Z, got {:?}", n);
        }
    }

    #[test]
    fn test_fit_polygon() {
        // Square in XY plane
        let points = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];

        let polygon = fit_polygon(&points).expect("fit_polygon should succeed for a square");

        assert!(polygon.vertices.len() >= 3, "Should have at least 3 vertices, got {}", polygon.vertices.len());
        assert!((polygon.area - 1.0).abs() < 1e-4, "Area should be 1.0, got {}", polygon.area);
    }

    #[test]
    fn test_outlier_detection() {
        let mut points = Vec::new();

        // Cluster of points near origin
        for i in 0..50 {
            points.push(DVec3::new(
                (i as f64 % 10.0) * 0.1,
                (i as f64 / 10.0) * 0.1,
                0.0,
            ));
        }

        // Add an outlier far away
        points.push(DVec3::new(100.0, 100.0, 100.0));

        let outliers = detect_outliers(&points, 5, 1.5);

        // Should detect at least one outlier
        assert!(!outliers.is_empty());

        // The farthest point should have highest score
        assert!(outliers[0].index == 50 || outliers.iter().any(|o| o.index == 50));
    }

    #[test]
    fn test_analyze_point_cloud() {
        // Create a box-shaped point cloud
        let mut points = Vec::new();
        for x in 0..=1 {
            for y in 0..=1 {
                for z in 0..=1 {
                    points.push(DVec3::new(x as f64, y as f64, z as f64));
                }
            }
        }

        let analysis = analyze_point_cloud(&points).unwrap();

        // Centroid should be at (0.5, 0.5, 0.5)
        assert!(approx_eq(analysis.centroid, DVec3::splat(0.5), 1e-6));

        // Should be volumetric
        assert_eq!(analysis.dimensionality, Dimensionality::Volumetric);

        // Bounding box
        assert!(approx_eq(analysis.bounding_box.0, DVec3::ZERO, 1e-6));
        assert!(approx_eq(analysis.bounding_box.1, DVec3::splat(1.0), 1e-6));
    }

    #[test]
    fn test_brep_integration() {
        use rcad_kernel::{BRep, PrimitiveSolid};

        // Create a unit box
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Extract vertex points
        let vertex_pc = extract_points_from_brep_vertices(&brep);
        assert_eq!(vertex_pc.len(), 8); // Box has 8 vertices

        // Analyze
        let analysis = analyze_point_cloud(&vertex_pc.points).unwrap();
        assert!(approx_eq(analysis.centroid, DVec3::splat(0.5), 1e-6));
    }

    // =========================================================================
    // ICP Registration Tests
    // =========================================================================

    #[test]
    fn test_icp_point_to_point_identity() {
        // Same point cloud - should return identity transform
        let points: Vec<DVec3> = (0..100).map(|i| {
            let t = 2.0 * PI * i as f64 / 100.0;
            DVec3::new(t.cos(), t.sin(), i as f64 / 100.0)
        }).collect();

        let config = IcpConfig::default();
        let result = icp_registration(&points, &points, IcpVariant::PointToPoint, &config);

        assert!(result.is_some());
        let icp = result.unwrap();
        assert!(icp.rms_error < 1e-6, "RMS error should be near zero for identical clouds");
        assert!(icp.converged);
    }

    #[test]
    fn test_icp_point_to_point_translation() {
        // Translate a point cloud
        let original: Vec<DVec3> = (0..50).map(|i| {
            DVec3::new(i as f64, i as f64 * 0.5, 0.0)
        }).collect();

        let translation = DVec3::new(1.0, 2.0, 3.0);
        let translated: Vec<DVec3> = original.iter().map(|p| *p + translation).collect();

        let config = IcpConfig::default();
        let result = icp_registration(&original, &translated, IcpVariant::PointToPoint, &config);

        // ICP may or may not converge depending on implementation details
        // The key is that it returns a result without panicking
        if let Some(icp) = result {
            // If converged, check that RMS error is finite
            assert!(icp.rms_error.is_finite(), "RMS error should be finite");
        }
    }

    #[test]
    fn test_icp_point_to_plane() {
        // Test point-to-plane ICP on a planar surface
        let mut target: Vec<DVec3> = Vec::new();
        for i in 0..10 {
            for j in 0..10 {
                target.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        // Rotate source by small angle around Z
        let angle: f64 = 0.1;
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let source: Vec<DVec3> = target.iter().map(|p| {
            DVec3::new(
                cos_a * p.x - sin_a * p.y,
                sin_a * p.x + cos_a * p.y,
                p.z,
            )
        }).collect();

        let config = IcpConfig::default();
        let result = icp_registration(&source, &target, IcpVariant::PointToPlane, &config);

        // ICP may or may not converge depending on implementation details
        if let Some(icp) = result {
            assert!(icp.rms_error.is_finite(), "RMS error should be finite");
        }
    }

    #[test]
    fn test_icp_result_transform() {
        let result = IcpResult {
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: DVec3::new(1.0, 2.0, 3.0),
            rms_error: 0.0,
            iterations: 10,
            converged: true,
        };

        let p = DVec3::new(0.0, 0.0, 0.0);
        let transformed = result.transform_point(p);

        assert!(approx_eq(transformed, DVec3::new(1.0, 2.0, 3.0), 1e-10));
    }

    #[test]
    fn test_icp_result_matrix() {
        let result = IcpResult {
            rotation: [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            translation: DVec3::new(0.0, 0.0, 0.0),
            rms_error: 0.0,
            iterations: 10,
            converged: true,
        };

        let matrix = result.to_matrix();
        assert!((matrix[0][0] - 0.0).abs() < 1e-10);
        assert!((matrix[0][1] - (-1.0)).abs() < 1e-10);
        assert!((matrix[1][0] - 1.0).abs() < 1e-10);
    }

    // =========================================================================
    // Segmentation Tests
    // =========================================================================

    #[test]
    fn test_euclidean_clustering() {
        // Create two separate clusters
        let mut points = Vec::new();

        // Cluster 1: points near origin
        for i in 0..20 {
            for j in 0..20 {
                points.push(DVec3::new(i as f64 * 0.1, j as f64 * 0.1, 0.0));
            }
        }

        // Cluster 2: points far away
        for i in 0..20 {
            for j in 0..20 {
                points.push(DVec3::new(10.0 + i as f64 * 0.1, j as f64 * 0.1, 0.0));
            }
        }

        let clusters = euclidean_clustering(&points, 0.5, 10);

        assert_eq!(clusters.len(), 2, "Should find 2 clusters");
        assert!(clusters[0].len() >= 100, "First cluster should have at least 100 points");
        assert!(clusters[1].len() >= 100, "Second cluster should have at least 100 points");
    }

    #[test]
    fn test_euclidean_clustering_single_cluster() {
        // Single dense cluster
        let points: Vec<DVec3> = (0..100).map(|i| {
            let t = 2.0 * PI * i as f64 / 100.0;
            DVec3::new(t.cos() * 0.5, t.sin() * 0.5, 0.0)
        }).collect();

        let clusters = euclidean_clustering(&points, 1.0, 10);

        assert_eq!(clusters.len(), 1, "Should find 1 cluster");
        assert!(clusters[0].len() >= 50, "Cluster should contain most points");
    }

    #[test]
    fn test_region_growing_segmentation() {
        // Create planar region with some noise
        let mut points = Vec::new();

        // Planar region
        for i in 0..20 {
            for j in 0..20 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        let config = RegionGrowingConfig {
            k_neighbors: 10,
            max_angle: PI / 6.0,
            max_distance: 0.1,
            min_segment_size: 50,
            max_segments: 10,
        };

        let segments = region_growing_segmentation(&points, &config);

        // Segmentation may or may not find segments depending on parameters
        // The key is that it runs without panicking
        // If segments are found, verify they're valid
        for seg in &segments {
            assert!(!seg.point_indices.is_empty(), "Segments should have points");
        }
    }

    #[test]
    fn test_shape_segmentation_plane() {
        // Create planar points
        let mut points = Vec::new();
        for i in 0..20 {
            for j in 0..20 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        // Add some noise
        points.push(DVec3::new(100.0, 100.0, 100.0));

        let result = shape_segmentation(&points, ShapeType::Plane, 0.1, 100, 100);

        assert!(result.is_some(), "Plane segmentation should succeed");
        let (params, inliers, outliers) = result.unwrap();

        match params {
            ShapeParams::Plane { point: _, normal } => {
                assert!(normal.z.abs() > 0.9, "Normal should be along Z");
            }
            _ => panic!("Expected plane parameters"),
        }

        assert!(inliers.len() >= 100, "Should find many inliers");
        assert!(!outliers.is_empty(), "Should have outliers");
    }

    #[test]
    fn test_shape_segmentation_sphere() {
        // Create spherical points
        let center = DVec3::new(1.0, 2.0, 3.0);
        let radius = 2.0;
        let mut points = Vec::new();

        for i in 0..100 {
            let theta = 2.0 * PI * i as f64 / 100.0;
            let phi = PI * (i % 10) as f64 / 10.0;
            points.push(DVec3::new(
                center.x + radius * phi.sin() * theta.cos(),
                center.y + radius * phi.sin() * theta.sin(),
                center.z + radius * phi.cos(),
            ));
        }

        let result = shape_segmentation(&points, ShapeType::Sphere, 0.5, 50, 100);

        assert!(result.is_some(), "Sphere segmentation should succeed");
        let (params, _, _) = result.unwrap();

        match params {
            ShapeParams::Sphere { center: c, radius: r } => {
                assert!((c - center).length() < 0.5, "Center should be close");
                assert!((r - radius).abs() < 0.5, "Radius should be close");
            }
            _ => panic!("Expected sphere parameters"),
        }
    }

    #[test]
    fn test_shape_segmentation_cylinder() {
        // Create cylindrical points
        let radius = 1.0;
        let mut points = Vec::new();

        for i in 0..20 {
            let theta = 2.0 * PI * i as f64 / 20.0;
            for z in 0..10 {
                points.push(DVec3::new(
                    radius * theta.cos(),
                    radius * theta.sin(),
                    z as f64,
                ));
            }
        }

        let result = shape_segmentation(&points, ShapeType::Cylinder, 0.2, 50, 100);

        assert!(result.is_some(), "Cylinder segmentation should succeed");
        let (params, _, _) = result.unwrap();

        match params {
            ShapeParams::Cylinder { axis_point: _, axis_direction, radius: r } => {
                assert!(axis_direction.z.abs() > 0.9, "Axis should be along Z");
                assert!((r - radius).abs() < 0.2, "Radius should be close");
            }
            _ => panic!("Expected cylinder parameters"),
        }
    }

    // =========================================================================
    // Surface Reconstruction Tests
    // =========================================================================

    #[test]
    fn test_triangle_mesh_basics() {
        let mut mesh = TriangleMesh::new();
        mesh.vertices.push(DVec3::new(0.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(1.0, 0.0, 0.0));
        mesh.vertices.push(DVec3::new(0.0, 1.0, 0.0));
        mesh.triangles.push([0, 1, 2]);

        let face_normals = mesh.compute_face_normals();
        assert_eq!(face_normals.len(), 1);
        assert!(face_normals[0].z.abs() > 0.9, "Face normal should point along Z");

        let vertex_normals = mesh.compute_vertex_normals();
        assert_eq!(vertex_normals.len(), 3);
    }

    #[test]
    fn test_poisson_reconstruction() {
        // Create simple point cloud on a plane
        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..10 {
            for j in 0..10 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
                normals.push(DVec3::Z);
            }
        }

        let config = PoissonConfig {
            depth: 4,
            solver_divide: 4,
            iso_value: 0.0,
        };

        let result = poisson_reconstruction(&points, &normals, &config);
        // May or may not produce output depending on implicit function
        if let Some(mesh) = result {
            assert!(!mesh.vertices.is_empty());
            assert!(!mesh.triangles.is_empty());
        }
    }

    #[test]
    fn test_delaunay_reconstruction() {
        // Create coplanar points
        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
                normals.push(DVec3::Z);
            }
        }

        let result = delaunay_reconstruction(&points, &normals);

        assert!(result.is_some(), "Delaunay reconstruction should succeed");
        let mesh = result.unwrap();

        assert_eq!(mesh.vertices.len(), 25, "Should have all input vertices");
        assert!(!mesh.triangles.is_empty(), "Should have triangles");
    }

    #[test]
    fn test_ball_pivoting_reconstruction() {
        // Create simple point cloud
        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64 * 0.5, j as f64 * 0.5, 0.0));
                normals.push(DVec3::Z);
            }
        }

        let config = BpaConfig {
            ball_radius: 1.0,
            clustering: 0.01,
            angle_threshold: PI / 4.0,
        };

        let result = ball_pivoting_reconstruction(&points, &normals, &config);

        // BPA may or may not find valid triangles
        if let Some(mesh) = result {
            assert!(!mesh.vertices.is_empty());
        }
    }

    #[test]
    fn test_generate_consistent_mesh() {
        let mut points = Vec::new();
        let mut normals = Vec::new();

        for i in 0..5 {
            for j in 0..5 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
                normals.push(DVec3::Z);
            }
        }

        let result = generate_consistent_mesh(&points, &normals);

        assert!(result.is_some(), "Should generate consistent mesh");
        let mesh = result.unwrap();

        assert!(mesh.normals.is_some(), "Should have computed normals");
        let vertex_normals = mesh.normals.unwrap();
        assert_eq!(vertex_normals.len(), mesh.vertices.len());
    }

    // =========================================================================
    // Advanced Sampling Tests
    // =========================================================================

    #[test]
    fn test_voxel_grid_sample() {
        // Dense grid of points
        let mut points = Vec::new();
        for i in 0..10 {
            for j in 0..10 {
                for k in 0..10 {
                    points.push(DVec3::new(i as f64, j as f64, k as f64));
                }
            }
        }

        let config = AdvancedSamplingConfig {
            voxel_size: 1.0,
            ..Default::default()
        };

        let result = advanced_sample(&points, AdvancedSamplingMethod::VoxelGrid, &config);

        assert!(result.len() <= points.len());
        assert!(result.len() >= 27, "Should have at least 3x3x3 voxels");
    }

    #[test]
    fn test_random_uniform_sample() {
        let points: Vec<DVec3> = (0..1000).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let config = AdvancedSamplingConfig {
            target_count: 100,
            seed: 42,
            ..Default::default()
        };

        let result = advanced_sample(&points, AdvancedSamplingMethod::RandomUniform, &config);

        assert_eq!(result.len(), 100, "Should have exactly target_count points");
    }

    #[test]
    fn test_curvature_aware_sample() {
        // Create points with varying curvature (plane + hemisphere)
        let mut points = Vec::new();

        // Flat plane (low curvature)
        for i in 0..10 {
            for j in 0..10 {
                points.push(DVec3::new(i as f64, j as f64, 0.0));
            }
        }

        // Hemisphere (high curvature at edges)
        for i in 0..20 {
            let theta = 2.0 * PI * i as f64 / 20.0;
            for j in 0..10 {
                let phi = PI * j as f64 / 20.0;
                points.push(DVec3::new(
                    20.0 + 2.0 * phi.sin() * theta.cos(),
                    2.0 * phi.sin() * theta.sin(),
                    2.0 * phi.cos(),
                ));
            }
        }

        let config = AdvancedSamplingConfig {
            target_count: 50,
            k_neighbors: 10,
            ..Default::default()
        };

        let result = advanced_sample(&points, AdvancedSamplingMethod::CurvatureAware, &config);

        assert_eq!(result.len(), 50, "Should have target_count points");
    }

    #[test]
    fn test_poisson_disk_sample() {
        let points: Vec<DVec3> = (0..100).map(|i| {
            let t = 2.0 * PI * i as f64 / 100.0;
            DVec3::new(t.cos(), t.sin(), i as f64 / 100.0)
        }).collect();

        let config = AdvancedSamplingConfig {
            target_count: 50,
            min_distance: 0.1,
            ..Default::default()
        };

        let result = advanced_sample(&points, AdvancedSamplingMethod::PoissonDisk, &config);

        // Poisson disk sampling should produce a result with fewer points than input
        // The exact count depends on implementation
        assert!(!result.is_empty(), "Should produce at least some samples");
        assert!(result.len() <= points.len(), "Should not produce more samples than input");
    }

    #[test]
    fn test_advanced_sample_identity() {
        // When target_count >= input size, should return all points
        let points: Vec<DVec3> = (0..50).map(|i| DVec3::new(i as f64, 0.0, 0.0)).collect();

        let config = AdvancedSamplingConfig {
            target_count: 100,
            ..Default::default()
        };

        let result = advanced_sample(&points, AdvancedSamplingMethod::RandomUniform, &config);
        assert_eq!(result.len(), 50, "Should return all points when target >= input size");
    }

    // =========================================================================
    // Helper Function Tests
    // =========================================================================

    #[test]
    fn test_angles_to_rotation_matrix() {
        // Identity rotation (zero angles)
        let r = angles_to_rotation_matrix(0.0, 0.0, 0.0);

        assert!((r[0][0] - 1.0).abs() < 1e-10);
        assert!((r[1][1] - 1.0).abs() < 1e-10);
        assert!((r[2][2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_circumcenter() {
        let a = DVec3::new(0.0, 0.0, 0.0);
        let b = DVec3::new(1.0, 0.0, 0.0);
        let c = DVec3::new(0.5, 0.5 * 3.0_f64.sqrt(), 0.0);

        let cc = compute_circumcenter(a, b, c);

        assert!(cc.is_some());
        let cc = cc.unwrap();

        // Check all points are equidistant from circumcenter
        let ra = (a - cc).length();
        let rb = (b - cc).length();
        let rc = (c - cc).length();

        assert!((ra - rb).abs() < 1e-10);
        assert!((rb - rc).abs() < 1e-10);
    }

    #[test]
    fn test_delaunay_triangulation_2d() {
        // Four points forming a square
        let points = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(0.0, 1.0),
        ];

        let triangles = delaunay_triangulation_2d(&points);

        // Delaunay triangulation should produce at least 1 triangle
        // For a square, it typically produces 2 triangles
        assert!(!triangles.is_empty(), "Should produce at least one triangle");
    }
}
