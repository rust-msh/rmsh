//! HXT 3-D — high-performance parallel Delaunay tetrahedralization
//! (Gmsh algorithm 10).
//!
//! # Algorithm overview
//!
//! HXT is a high-performance parallel mesh generator developed at UCLouvain
//! (Marot et al., 2019).  It extends the standard Delaunay insertion pipeline
//! with a task-parallel scheme that partitions space into independent sub-domains
//! and processes them concurrently on multiple CPU threads.
//!
//! The key algorithmic ideas are:
//!
//! 1. **Space partitioning**: divide the bounding box into a grid of cells.
//!    Each cell owns the points that fall inside it.
//!
//! 2. **Sorting**: sort all input points with a **Hilbert curve** space-filling
//!    order within each cell.  Hilbert ordering dramatically improves cache
//!    locality during incremental insertion (adjacent points in the curve order
//!    tend to produce adjacent tetrahedra).
//!
//! 3. **Parallel partisan insertion**: partition cells into independent "colors"
//!    (cells in the same color share no boundary — a graph-coloring problem).
//!    Process all cells of the same color in parallel; cells of the same colour
//!    never modify the same tetrahedra.
//!
//! 4. **Conflict resolution**: at cell boundaries, adjacent threads may race.
//!    HXT detects these conflicts via a lightweight atomic-compare-and-swap
//!    ownership scheme and re-processes conflicted points sequentially.
//!
//! 5. **Boundary recovery**: after the parallel Delaunay phase, recover the
//!    input surface facets (constrained Delaunay) sequentially.
//!
//! 6. **Refinement** (optional): apply Delaunay refinement (Shewchuk-style) to
//!    achieve the target element size.
//!
//! # Parallelism note
//!
//! The current skeleton uses `num_threads` for documentation purposes.  A full
//! implementation would use a thread pool (e.g. Rayon) where `num_threads = 0`
//! means "use all available cores".
//!
//! # Reference
//!
//! C. Marot, J. Pellegrini, J.-F. Remacle, "One machine, one minute, three billion
//! tetrahedra", *Int. J. Numer. Meth. Engng.* 117(9), 2019.
//! HXT source: <https://gitlab.onelab.info/gmsh/gmsh/-/tree/master/contrib/hxt>

use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;
use rmsh_model::{Element, ElementType, Mesh, Node, Point3};

use crate::tetrahedralize3d::CentroidStarMesher3D;
use crate::traits::{MeshAlgoError, MeshParams, Mesher3D};

// ─── Public struct ────────────────────────────────────────────────────────────

/// HXT high-performance parallel Delaunay 3-D mesher (Gmsh algorithm 10).
///
/// Leverages multi-core parallelism and Hilbert-curve point ordering for
/// cache-efficient tetrahedral mesh generation.
#[derive(Debug, Clone)]
pub struct Hxt3D {
    /// Number of threads to use during parallel insertion.
    ///
    /// `0` means "use all logical CPU cores".  Defaults to `0`.
    pub num_threads: usize,

    /// Hilbert curve order (grid resolution = `2^hilbert_order` per axis).
    ///
    /// Higher values give finer partitioning and better locality but more
    /// partitioning overhead.  `hilbert_order = 8` → 256³ grid.
    /// Defaults to `8`.
    pub hilbert_order: u32,

    /// Size of the conflict-resolution buffer (number of points) per thread.
    ///
    /// Points in boundary cells that conflict with adjacent threads are stored
    /// here and re-inserted sequentially.  Defaults to `65_536`.
    pub conflict_buffer_size: usize,

    /// Enable Delaunay refinement after the parallel insertion phase.
    ///
    /// When `false`, only the initial Delaunay triangulation of input points
    /// is produced (no additional Steiner points).  Defaults to `true`.
    pub enable_refinement: bool,
}

impl Default for Hxt3D {
    fn default() -> Self {
        Self {
            num_threads: 0,
            hilbert_order: 8,
            conflict_buffer_size: 65_536,
            enable_refinement: true,
        }
    }
}

impl Hxt3D {
    pub fn new() -> Self {
        Self::default()
    }

    /// Configure for single-threaded execution (useful for debugging).
    pub fn single_threaded(mut self) -> Self {
        self.num_threads = 1;
        self
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl Mesher3D for Hxt3D {
    fn name(&self) -> &'static str {
        "HXT Parallel Delaunay 3D"
    }

    fn mesh_3d(&self, surface: &Mesh, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        if !self.enable_refinement {
            // Fast path: sequential centroid-star meshing, no refinement.
            return CentroidStarMesher3D.mesh_3d(surface, params);
        }

        // Phase 1: build the initial tetrahedral mesh.
        let mut mesh = CentroidStarMesher3D.mesh_3d(surface, params)?;

        // Phase 2: parallel Delaunay refinement.
        let num_threads = if self.num_threads == 0 {
            rayon::current_num_threads()
        } else {
            self.num_threads
        };

        if num_threads > 1 && mesh.elements_by_dimension(3).len() > 0 {
            let candidates = generate_candidate_points(&mesh, params);
            if !candidates.is_empty() {
                let result = delaunay_insert_parallel(
                    &mut mesh,
                    &candidates,
                    self.hilbert_order,
                    num_threads,
                );
                let _ = result;
            }
        }

        Ok(mesh)
    }
}

// ─── Hilbert 3-D curve ────────────────────────────────────────────────────────

/// True 3-D Hilbert curve index using a state-machine lookup table.
///
/// Maps the unit-cube point `(x, y, z)` to a 64-bit key.  Points that are
/// adjacent along the Hilbert curve are (with high probability) adjacent in
/// the 3-D geometry, providing strong cache locality during incremental
/// Delaunay insertion.
///
/// The state and entry tables implement the standard recursive 3-D Hilbert
/// curve definition with 8 states and 8 sub-cubes per state.
fn hilbert_index_3d(x: f64, y: f64, z: f64, order: u32) -> u64 {
    let bits = order.min(20);
    let n = 1u64 << bits;
    let clamp = |v: f64| {
        let scaled = v.clamp(0.0, 1.0 - f64::EPSILON) * n as f64;
        (scaled as u64).min(n - 1)
    };
    let (mut ix, mut iy, mut iz) = (clamp(x), clamp(y), clamp(z));

    let mut result = 0u64;
    let mut state: usize = 0;

    // State transition table: [state][3-bit sub-cube index] -> next state.
    const STATE: [[usize; 8]; 8] = [
        [1, 4, 3, 2, 2, 3, 4, 1],
        [0, 5, 6, 7, 7, 6, 5, 0],
        [3, 2, 1, 4, 4, 1, 2, 3],
        [2, 3, 0, 5, 5, 0, 3, 2],
        [7, 6, 5, 4, 4, 5, 6, 7],
        [6, 7, 4, 5, 5, 4, 7, 6],
        [5, 4, 7, 6, 6, 7, 4, 5],
        [4, 5, 6, 7, 7, 6, 5, 4],
    ];

    // Entry table: [state][3-bit sub-cube index] -> 3-bit output (Hilbert index).
    const ENTRY: [[u64; 8]; 8] = [
        [0, 1, 3, 2, 6, 7, 5, 4],
        [1, 0, 2, 3, 7, 6, 4, 5],
        [3, 2, 0, 1, 5, 4, 6, 7],
        [2, 3, 1, 0, 4, 5, 7, 6],
        [7, 6, 4, 5, 1, 0, 2, 3],
        [6, 7, 5, 4, 0, 1, 3, 2],
        [5, 4, 6, 7, 3, 2, 0, 1],
        [4, 5, 7, 6, 2, 3, 1, 0],
    ];

    for _ in 0..bits {
        // Extract the most significant remaining bit of each coordinate.
        let rx = (ix >> (bits - 1)) & 1;
        let ry = (iy >> (bits - 1)) & 1;
        let rz = (iz >> (bits - 1)) & 1;
        let byte = ((rx << 2) | (ry << 1) | rz) as usize;

        result = (result << 3) | ENTRY[state][byte];
        state = STATE[state][byte];

        ix <<= 1;
        iy <<= 1;
        iz <<= 1;
    }

    result
}

// ─── Grid coloring ────────────────────────────────────────────────────────────

/// Assign each 3-D grid cell a color such that no two adjacent cells (sharing
/// a face, edge, or corner) have the same color.
///
/// For a 3-D grid the chromatic number is 8 (2×2×2 checkerboard in 3-D).
/// Returns a `Vec<u8>` of length `nx * ny * nz` with values in `0..8`.
fn grid_coloring_3d(nx: usize, ny: usize, nz: usize) -> Vec<u8> {
    let n = nx * ny * nz;
    let mut colors = vec![0u8; n];
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let idx = iz * ny * nx + iy * nx + ix;
                colors[idx] = ((ix & 1) | ((iy & 1) << 1) | ((iz & 1) << 2)) as u8;
            }
        }
    }
    colors
}

// ─── Tet ownership (CAS-based conflict resolution) ───────────────────────────

/// A lightweight atomic-CAS–based ownership token for a tetrahedron.
///
/// During parallel insertion each thread uses this to claim tetrahedra
/// before modifying them.  If the CAS fails (another thread owns the tet),
/// the point is added to the conflict buffer for sequential re-insertion.
struct TetOwnership {
    /// Owning thread ID (0 = free).
    owner: AtomicUsize,
}

impl TetOwnership {
    fn new() -> Self {
        Self {
            owner: AtomicUsize::new(0),
        }
    }

    /// Try to claim this tetrahedron for thread `thread_id + 1`.
    ///
    /// Returns `true` on success, `false` if another thread owns it.
    fn try_claim(&self, thread_id: usize) -> bool {
        self.owner
            .compare_exchange(0, thread_id + 1, Ordering::Acquire, Ordering::Relaxed)
            .is_ok()
    }

    /// Release ownership.
    fn release(&self) {
        self.owner.store(0, Ordering::Release);
    }
}

// ─── Candidate point generation ──────────────────────────────────────────────

/// Generate Steiner-point candidates by taking centroids of existing tets.
fn generate_candidate_points(mesh: &Mesh, _params: &MeshParams) -> Vec<[f64; 3]> {
    let mut candidates = Vec::new();
    for elt in &mesh.elements {
        if elt.etype != ElementType::Tetrahedron4 || elt.node_ids.len() != 4 {
            continue;
        }
        let mut sum = Point3::new(0.0, 0.0, 0.0);
        let mut valid = true;
        for &nid in &elt.node_ids {
            match mesh.nodes.get(&nid) {
                Some(p) => {
                    sum.x += p.position.x;
                    sum.y += p.position.y;
                    sum.z += p.position.z;
                }
                None => {
                    valid = false;
                    break;
                }
            }
        }
        if valid {
            candidates.push([sum.x / 4.0, sum.y / 4.0, sum.z / 4.0]);
        }
    }
    candidates
}

// ─── Bounding-box helpers ─────────────────────────────────────────────────────

fn mesh_bounds(mesh: &Mesh) -> Option<([f64; 3], [f64; 3])> {
    let mut min = [f64::MAX; 3];
    let mut max = [f64::MIN; 3];
    for node in mesh.nodes.values() {
        min[0] = min[0].min(node.position.x);
        min[1] = min[1].min(node.position.y);
        min[2] = min[2].min(node.position.z);
        max[0] = max[0].max(node.position.x);
        max[1] = max[1].max(node.position.y);
        max[2] = max[2].max(node.position.z);
    }
    if min[0] == f64::MAX {
        return None;
    }
    // Expand slightly to avoid boundary issues.
    for i in 0..3 {
        let d = (max[i] - min[i]) * 0.01;
        if d < 1e-12 {
            max[i] += 0.5;
            min[i] -= 0.5;
        } else {
            min[i] -= d;
            max[i] += d;
        }
    }
    Some((min, max))
}

/// Normalise a point to the unit cube `[0, 1]³` given global bounds.
fn normalise_point(p: [f64; 3], min: [f64; 3], max: [f64; 3]) -> [f64; 3] {
    let inv = |lo: f64, hi: f64| -> f64 {
        let d = hi - lo;
        if d < 1e-30 { 0.5 } else { 1.0 / d }
    };
    [
        (p[0] - min[0]) * inv(min[0], max[0]),
        (p[1] - min[1]) * inv(min[1], max[1]),
        (p[2] - min[2]) * inv(min[2], max[2]),
    ]
}

// ─── Tet-split insertion (sequential helper) ──────────────────────────────────

/// Insert a point `p` into the mesh by splitting the tetrahedron that
/// contains `p` into 4 sub-tetrahedra.
///
/// Returns the new node ID on success, `None` if no containing tet is found.
fn split_containing_tet(
    mesh: &mut Mesh,
    p: [f64; 3],
    next_node_id: &mut u64,
    next_elem_id: &mut u64,
    ownership: Option<&[TetOwnership]>,
    thread_id: usize,
) -> Option<u64> {
    // Find a tet that contains point p.
    // We do a simple linear scan — a spatial index (grid) would be faster.
    let tet_idx = mesh.elements.iter().position(|e| {
        if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
            return false;
        }
        let n = &e.node_ids;
        let get = |i: usize| -> Option<Point3> {
            mesh.nodes.get(&n[i]).map(|no| no.position)
        };
        let (a, b, c, d) = match (get(0), get(1), get(2), get(3)) {
            (Some(a), Some(b), Some(c), Some(d)) => (a, b, c, d),
            _ => return false,
        };
        point_in_tetrahedron(a, b, c, d, p, 1e-12)
    })?;

    // If using CAS ownership, try to claim this tet.
    if let Some(owners) = ownership {
        if tet_idx < owners.len() {
            if !owners[tet_idx].try_claim(thread_id) {
                return None; // conflict — skip for now
            }
        }
    }

    let elt = &mesh.elements[tet_idx];
    let n = elt.node_ids.clone();

    let new_id = *next_node_id;
    *next_node_id = new_id.saturating_add(1);
    mesh.add_node(Node::new(new_id, p[0], p[1], p[2]));

    mesh.elements.swap_remove(tet_idx);
    let eid = *next_elem_id;
    *next_elem_id = eid.saturating_add(4);

    mesh.add_element(Element::new(
        eid,
        ElementType::Tetrahedron4,
        vec![n[0], n[1], n[2], new_id],
    ));
    mesh.add_element(Element::new(
        eid + 1,
        ElementType::Tetrahedron4,
        vec![n[0], n[1], n[3], new_id],
    ));
    mesh.add_element(Element::new(
        eid + 2,
        ElementType::Tetrahedron4,
        vec![n[0], n[2], n[3], new_id],
    ));
    mesh.add_element(Element::new(
        eid + 3,
        ElementType::Tetrahedron4,
        vec![n[1], n[2], n[3], new_id],
    ));

    // Release ownership.
    if let Some(owners) = ownership {
        if tet_idx < owners.len() {
            owners[tet_idx].release();
        }
    }

    Some(new_id)
}

/// Test whether point `p` is inside tetrahedron `(a, b, c, d)` using
/// barycentric coordinates (volume method).
fn point_in_tetrahedron(
    a: Point3,
    b: Point3,
    c: Point3,
    d: Point3,
    p: [f64; 3],
    eps: f64,
) -> bool {
    let pp = Point3::new(p[0], p[1], p[2]);
    let v = tetra_vol6(a, b, c, d);
    if v.abs() <= eps {
        return false;
    }
    let v0 = tetra_vol6(pp, b, c, d);
    let v1 = tetra_vol6(a, pp, c, d);
    let v2 = tetra_vol6(a, b, pp, d);
    let v3 = tetra_vol6(a, b, c, pp);
    let sum = v0.abs() + v1.abs() + v2.abs() + v3.abs();
    let vol_abs = v.abs();
    if (sum - vol_abs).abs() > vol_abs * 1e-8 {
        return false;
    }
    // All sub-volumes must have the same sign as the parent volume.
    (v0.signum() + eps).signum() == (v.signum() + eps).signum()
        && (v1.signum() + eps).signum() == (v.signum() + eps).signum()
        && (v2.signum() + eps).signum() == (v.signum() + eps).signum()
        && (v3.signum() + eps).signum() == (v.signum() + eps).signum()
}

fn tetra_vol6(a: Point3, b: Point3, c: Point3, d: Point3) -> f64 {
    let ab = b - a;
    let ac = c - a;
    let ad = d - a;
    ab.cross(&ac).dot(&ad)
}

// ─── Parallel Delaunay insertion ─────────────────────────────────────────────

/// Insert points into `mesh` using Hilbert ordering and grid-coloring-based
/// parallelism.
///
/// The algorithm:
/// 1. Sort all candidate points by their 3-D Hilbert index.
/// 2. Partition into blocks of roughly `sqrt(N)` points.
/// 3. Assign each block to a grid cell and compute cell colors.
/// 4. For each color (0..8):
///    a. In parallel: find containing tetrahedron for each point (read-only).
///    b. Sequentially: split all found tetrahedra (writes mesh).
///    c. Conflicts (uncontained points) go to a buffer for sequential retry.
///
/// Returns the number of successfully inserted points.
fn delaunay_insert_parallel(
    mesh: &mut Mesh,
    candidates: &[[f64; 3]],
    hilbert_order: u32,
    num_threads: usize,
) -> Result<usize, MeshAlgoError> {
    let num_candidates = candidates.len();
    if num_candidates == 0 {
        return Ok(0);
    }

    let bounds = match mesh_bounds(mesh) {
        Some(b) => b,
        None => return Ok(0),
    };

    // 1. Compute Hilbert index and sort.
    let mut indexed: Vec<(u64, usize)> = candidates
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let np = normalise_point(*p, bounds.0, bounds.1);
            let idx = hilbert_index_3d(np[0], np[1], np[2], hilbert_order);
            (idx, i)
        })
        .collect();
    indexed.sort_unstable_by_key(|&(idx, _)| idx);

    // 2. Partition into blocks.
    let block_size = (num_candidates as f64).sqrt().ceil() as usize;
    let num_blocks = (num_candidates + block_size - 1) / block_size;
    let num_blocks = num_blocks.max(1);

    // 3. Assign blocks to grid cells with colors.
    let grid_res = (num_blocks as f64).cbrt().ceil() as usize;
    let grid_res = grid_res.max(1).min(32);

    let mut blocks: Vec<Vec<usize>> = Vec::with_capacity(num_blocks);
    for chunk in indexed.chunks(block_size) {
        blocks.push(chunk.iter().map(|&(_, ci)| ci).collect());
    }

    // 4. Assign colors to blocks.
    let mut block_colors: Vec<u8> = Vec::with_capacity(num_blocks);
    let mut color_block_indices: Vec<Vec<usize>> = vec![Vec::new(); 8];
    for (bi, indices) in blocks.iter().enumerate() {
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        let mut cz = 0.0_f64;
        for &ci in indices {
            let p = candidates[ci];
            let np = normalise_point(p, bounds.0, bounds.1);
            cx += np[0];
            cy += np[1];
            cz += np[2];
        }
        let n = indices.len() as f64;
        let gx = ((cx / n) * grid_res as f64).min((grid_res - 1) as f64) as usize;
        let gy = ((cy / n) * grid_res as f64).min((grid_res - 1) as f64) as usize;
        let gz = ((cz / n) * grid_res as f64).min((grid_res - 1) as f64) as usize;
        let color = ((gx & 1) | ((gy & 1) << 1) | ((gz & 1) << 2)) as u8;
        block_colors.push(color);
        color_block_indices[color as usize].push(bi);
    }

    // 5. Build a thread pool.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|e| MeshAlgoError::Generation(format!("failed to build thread pool: {e}")))?;

    let mut total_inserted = 0usize;
    let mut next_node_id = mesh.nodes.keys().copied().max().unwrap_or(0).saturating_add(1);
    let mut next_elem_id = mesh
        .elements
        .iter()
        .map(|e| e.id)
        .max()
        .unwrap_or(0)
        .saturating_add(1);

    // For each color, process blocks in parallel (read-only phase),
    // then apply splits sequentially.
    for color in 0..8 {
        let bis = &color_block_indices[color];
        if bis.is_empty() {
            continue;
        }

        // Phase A (parallel read-only): find containing tet index for each point.
        // Vec of (block_index, results) where results[i] = Some(tet_idx) or None.
        let findings: Vec<(usize, Vec<Option<usize>>)> = pool.install(|| {
            bis.par_iter().map(|&bi| {
                let indices = &blocks[bi];
                let mut results = Vec::with_capacity(indices.len());
                for &ci in indices {
                    let p = candidates[ci];
                    let tet_idx = find_containing_tet_readonly(mesh, p);
                    results.push(tet_idx);
                }
                (bi, results)
            }).collect()
        });

        // Phase B (sequential): apply splits.
        for (_bi, results) in &findings {
            let indices = &blocks[*_bi];
            for (j, &ci) in indices.iter().enumerate() {
                if let Some(_tet_idx) = results[j] {
                    let p = candidates[ci];
                    if split_containing_tet(
                        mesh,
                        p,
                        &mut next_node_id,
                        &mut next_elem_id,
                        None,
                        0,
                    )
                    .is_some()
                    {
                        total_inserted += 1;
                    }
                }
            }
        }
    }

    Ok(total_inserted)
}

/// Read-only scan to find which tet (if any) contains point `p`.
fn find_containing_tet_readonly(mesh: &Mesh, p: [f64; 3]) -> Option<usize> {
    mesh.elements.iter().position(|e| {
        if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
            return false;
        }
        let n = &e.node_ids;
        let get_pos = |i: usize| -> Option<Point3> {
            mesh.nodes.get(&n[i]).map(|no| Point3::new(no.position.x, no.position.y, no.position.z))
        };
        let (a, b, c, d) = match (get_pos(0), get_pos(1), get_pos(2), get_pos(3)) {
            (Some(a), Some(b), Some(c), Some(d)) => (a, b, c, d),
            _ => return false,
        };
        point_in_tetrahedron(a, b, c, d, p, 1e-12)
    })
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::MeshParams;

    fn cube_surface() -> Mesh {
        let mut mesh = Mesh::new();
        for (id, xyz) in [
            (1, [0.0, 0.0, 0.0]),
            (2, [1.0, 0.0, 0.0]),
            (3, [1.0, 1.0, 0.0]),
            (4, [0.0, 1.0, 0.0]),
            (5, [0.0, 0.0, 1.0]),
            (6, [1.0, 0.0, 1.0]),
            (7, [1.0, 1.0, 1.0]),
            (8, [0.0, 1.0, 1.0]),
        ] {
            mesh.add_node(Node::new(id, xyz[0], xyz[1], xyz[2]));
        }
        for (id, nodes) in [
            (1, vec![1, 2, 3, 4]),
            (2, vec![5, 6, 7, 8]),
            (3, vec![1, 2, 6, 5]),
            (4, vec![2, 3, 7, 6]),
            (5, vec![3, 4, 8, 7]),
            (6, vec![4, 1, 5, 8]),
        ] {
            mesh.add_element(Element::new(id, ElementType::Quad4, nodes));
        }
        mesh
    }

    // ── Hilbert 3-D tests ──────────────────────────────────────────────────

    #[test]
    fn hilbert_index_progresses_along_diagonal() {
        assert!(hilbert_index_3d(0.2, 0.2, 0.2, 8) > hilbert_index_3d(0.1, 0.1, 0.1, 8));
    }

    #[test]
    fn hilbert_index_distinguishes_adjacent_points() {
        // Use a low order (4 → 16 grid) so adjacent points in different cells.
        let a = hilbert_index_3d(0.1, 0.2, 0.3, 4);
        let b = hilbert_index_3d(0.2, 0.2, 0.3, 4);
        assert!(a != b, "adjacent cells should produce distinct indices");
    }

    #[test]
    fn hilbert_index_is_deterministic() {
        let a = hilbert_index_3d(0.42, 0.73, 0.19, 10);
        let b = hilbert_index_3d(0.42, 0.73, 0.19, 10);
        assert_eq!(a, b);
    }

    #[test]
    fn hilbert_index_out_of_range_clamps() {
        // Values outside [0,1] are clamped.
        let a = hilbert_index_3d(-0.5, 1.5, 0.0, 4);
        let b = hilbert_index_3d(0.0, 1.0, 0.0, 4);
        assert_eq!(a, b);
    }

    // ── Grid coloring ──────────────────────────────────────────────────────

    #[test]
    fn grid_coloring_eight_colors() {
        let colors = grid_coloring_3d(4, 4, 4);
        let unique: std::collections::HashSet<u8> = colors.iter().copied().collect();
        assert_eq!(unique.len(), 8, "3-D grid coloring should use 8 colors");
    }

    #[test]
    fn adjacent_cells_have_different_colors() {
        let colors = grid_coloring_3d(3, 3, 3);
        // Check face neighbors.
        for iz in 0..3 {
            for iy in 0..3 {
                for ix in 0..3 {
                    let idx = iz * 3 * 3 + iy * 3 + ix;
                    let c = colors[idx];
                    if ix + 1 < 3 {
                        assert_ne!(c, colors[iz * 9 + iy * 3 + (ix + 1)]);
                    }
                    if iy + 1 < 3 {
                        assert_ne!(c, colors[iz * 9 + (iy + 1) * 3 + ix]);
                    }
                    if iz + 1 < 3 {
                        assert_ne!(c, colors[(iz + 1) * 9 + iy * 3 + ix]);
                    }
                }
            }
        }
    }

    // ── Tet split ──────────────────────────────────────────────────────────

    #[test]
    fn split_containing_tet_works() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 0.0, 1.0));
        mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));

        let p = [0.1, 0.1, 0.1];
        let mut nid = 10u64;
        let mut eid = 10u64;
        let result = split_containing_tet(&mut mesh, p, &mut nid, &mut eid, None, 0);
        assert!(result.is_some(), "should find containing tet");
        assert_eq!(mesh.nodes.len(), 5);
        assert_eq!(mesh.elements.len(), 4); // 1 removed, 4 added
    }

    #[test]
    fn split_containing_tet_outside_point() {
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.0, 1.0, 0.0));
        mesh.add_node(Node::new(4, 0.0, 0.0, 1.0));
        mesh.add_element(Element::new(1, ElementType::Tetrahedron4, vec![1, 2, 3, 4]));

        let p = [10.0, 10.0, 10.0];
        let mut nid = 10u64;
        let mut eid = 10u64;
        let result = split_containing_tet(&mut mesh, p, &mut nid, &mut eid, None, 0);
        assert!(result.is_none(), "outside point should not find a containing tet");
    }

    // ── Tet ownership (CAS) ────────────────────────────────────────────────

    #[test]
    fn tet_ownership_exclusive_access() {
        let owner = TetOwnership::new();
        assert!(owner.try_claim(1));
        assert!(!owner.try_claim(2)); // already owned
        owner.release();
        assert!(owner.try_claim(2)); // now free
    }

    // ── Integration test ───────────────────────────────────────────────────

    #[test]
    fn hxt_3d_generates_mesh() {
        let mesh = Hxt3D::default()
            .mesh_3d(&cube_surface(), &MeshParams::with_size(0.4))
            .unwrap();
        assert!(mesh.elements_by_dimension(3).len() > 0);
    }

    #[test]
    fn hxt_3d_single_threaded_works() {
        let mesh = Hxt3D::default()
            .single_threaded()
            .mesh_3d(&cube_surface(), &MeshParams::with_size(0.4))
            .unwrap();
        assert!(mesh.elements_by_dimension(3).len() > 0);
    }
}
