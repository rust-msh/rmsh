//! HXT 3-D — the rayon-parallel optimization layer over the Delaunay
//! tetrahedralization pipeline (Gmsh algorithm 10).
//!
//! # Honest scoping (read this first)
//!
//! The **full HXT parallel-insertion algorithm** — space partitioning into a
//! cell grid, Hilbert-curve point sorting, graph-colored *partisan insertion*
//! of points with atomic-CAS conflict resolution over concurrent cavity
//! operations, and concurrent boundary recovery — is research-grade
//! machinery (Marot et al. 2019 build three billion tets per minute with
//! it). It is **out of scope of this card and remains the documented
//! long-term path**: the insertion phase below still delegates to the
//! sequential [`Delaunay3D`] / [`CentroidStarMesher3D`] pipelines, and the
//! scheme-level helpers at the bottom of this file (`hilbert_index_3d`,
//! `grid_coloring_3d`, `TetOwnership`) document that path without
//! pretending to implement it.
//!
//! What this card delivers is the **HXT-adjacent parallelism that is
//! tractable and testable today**: the *optimization layer* of the pipeline,
//! parallelized with rayon over the same structural idea HXT uses for
//! insertion — **graph coloring for conflict-free independence**:
//!
//! 1. **Graph-coloring node independence** ([`color_node_adjacency`]): the
//!    tet mesh's node-adjacency graph is colored with a deterministic greedy
//!    first-fit (ascending compact index; 8-16 colors typical for tet
//!    meshes — the measured count is pinned per fixture in the tests).
//!    Nodes of the same color share no edge, hence no tetrahedron: their
//!    relocations touch disjoint stars and are conflict-free by
//!    construction.
//! 2. **Parallel relocation evaluation + sequential deterministic commit**
//!    ([`colored_parallel_relocate_pass`]): per color class, the relocation
//!    candidates of all its nodes are evaluated with rayon `par_iter` (the
//!    *evaluation* is parallel — it is the expensive part: per-axis
//!    golden-section search over the star), then the accepted moves are
//!    committed sequentially in ascending node order. See "Determinism
//!    contract" below.
//! 3. **Parallel quality evaluation** ([`parallel_mesh_quality_stats`]):
//!    the `mesh_quality_stats`-style measures (min dihedral, radius-edge
//!    distribution, sliver fraction) are computed per tet with rayon —
//!    embarrassingly parallel.
//! 4. **Wiring**: [`Hxt3D::mesh_3d`] runs the parallel optimization phase
//!    after insertion; [`Hxt3D::parallel_optimize_mesh`] is the standalone
//!    documented entry point for any pure tetrahedral mesh. `num_threads`
//!    is now meaningful: `0` = all logical cores (rayon default).
//!
//! The **topology-changing passes** (edge collapse, 2-3 swap, 3-2 flip) are
//! reused *verbatim* from [`crate::mesh_optimize`] and remain sequential by
//! design: each accepted operation rebuilds global connectivity, so their
//! acceptance loop is inherently serial. This is stated, not hidden — the
//! parallel layer of v1 is the node/element *evaluation* parallelism above.
//!
//! # Determinism contract (v1)
//!
//! Same input + same config ⇒ **bit-identical output mesh for any
//! `num_threads`** (pinned by tests, including a 1-thread vs all-core
//! bit-equality check) **and across processes**: 2026-10-14 (the drift
//! card) verified five separate test-binary runs over the cube /
//! stretched / slender fixture family byte-identical (connectivity +
//! coordinate digests). The earlier "cross-binary min-dihedral drift"
//! (19.777 vs 21.944 in the defect ledger) resolved as an
//! input-confusion artifact — the two numbers came from two DIFFERENT
//! box fixtures: this file's `box_surface` test helper has an inward
//! bottom face while `corpus::box_surface` is outward-oriented, so
//! 19.777 (this fixture) and 21.944 (the corpus fixture, pre-quality-
//! card seed) were never the same input. The policy that delivers the
//! determinism:
//!
//! * Per-node relocation evaluation is a **pure function** of the node
//!   positions at the moment its color class is processed (the evaluation
//!   reads the star snapshot, writes nothing) — so parallel evaluation
//!   produces exactly the floats a sequential evaluation would.
//! * Color classes are processed in a fixed order (ascending color id);
//!   within a class, accepted moves are committed **sequentially in
//!   ascending node order** — the commit order *is* the determinism policy.
//! * Within one color class no two nodes share a tet, so committing one
//!   accepted move cannot invalidate another candidate's acceptance gates
//!   (the gates were evaluated over stars that the move does not touch).
//! * Parallel float reductions in the quality stats are restricted to
//!   order-insensitive operations (`min`), and the p95 index is taken over
//!   the element-order ratio vector — bit-identical to the serial stats.
//!
//! Note the colored pass is *not* equivalent to the serial
//! `relocate_pass` (which commits immediately, so a later node observes an
//! earlier adjacent node's move even within one sweep): adjacency ordering
//! can differ between the two, so the parallel pipeline is pinned for
//! **quality parity** against the serial optimizer, not bit-equality. The
//! measured contract on the box / stretched / slender fixtures
//! (bit-deterministic, so the values are stable — full pins in the tests):
//!
//! * **box** (input min 5.625 deg, p95 20.629): serial(min 20.611, p95
//!   6.107, sliver 0) vs parallel(min 19.777, p95 5.317, sliver 0) — the
//!   parallel pipeline **beats** the radius-edge axis (the charter gate
//!   axis) by 13% and degrades the worst single dihedral by 0.83 deg (4%):
//!   the measured cost of the color-class ordering vs the serial index
//!   ordering, pinned at `>= serial - 1.0 deg`. Both paths lift the input
//!   worst dihedral by >3x.
//! * **stretched**: all three axes bit-equal to serial (min 0.176, p95
//!   1204.910, sliver 0.2419).
//! * **slender**: min and p95 bit-equal to serial (0.069 deg / 56675.829);
//!   sliver count identical (50 == 50), fraction +0.005 only from the
//!   smaller denominator (one extra accepted collapse, n 96 vs 97).
//! * Threading evidence: `HxtParallelOptimizeStats::threads_observed`
//!   counts distinct rayon workers that executed relocation evaluations; a
//!   forced 4-thread pool on the slender fixture observes 2 workers, the
//!   default 24-core pool observes up to 9 on the box.
//!
//! # Algorithm overview (the long-term HXT insertion scheme)
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
//! # Reference
//!
//! C. Marot, J. Pellegrini, J.-F. Remacle, "One machine, one minute, three billion
//! tetrahedra", *Int. J. Numer. Meth. Engng.* 117(9), 2019.
//! HXT source: <https://gitlab.onelab.info/gmsh/gmsh/-/tree/master/contrib/hxt>

use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};

use rmsh_model::{ElementType, Mesh};

use crate::delaunay_3d::Delaunay3D;
use crate::mesh_optimize::{
    edge_collapse_pass, edge_removal_pass_3_2, edge_swap_pass, free_nodes, is_sliver_tet,
    mean_edge_length, mesh_worst_quality_tuple, min_dihedral_deg_tet,
    normalize_orientations, radius_edge_ratio_tet, relocate_candidate_for_node,
    signed_volume_tet, star_tets, tetmesh_to_mesh_referenced, OptimizeConfig, QualityMetric,
    TetMeshQualityStats, VOL_EPS,
};
use crate::tet_mesh::TetMesh;
use crate::tetrahedralize3d::CentroidStarMesher3D;
use crate::traits::{MeshAlgoError, MeshParams, Mesher3D, OptimizeParams};

use rayon::prelude::*;

// ─── Public struct ────────────────────────────────────────────────────────────

/// HXT high-performance parallel 3-D mesher (Gmsh algorithm 10): the
/// sequential Delaunay insertion pipeline plus the **rayon-parallel
/// optimization layer** (see the module docs for the scoping and the
/// determinism contract).
#[derive(Debug, Clone)]
pub struct Hxt3D {
    /// Number of threads for the parallel optimization layer.
    ///
    /// `0` means "use all logical CPU cores" (the rayon default).
    /// Defaults to `0`.
    pub num_threads: usize,

    /// Hilbert curve order (grid resolution = `2^hilbert_order`).
    ///
    /// Reserved for the long-term parallel-insertion path (module docs,
    /// idea 2); unused by the v1 optimization layer.
    /// Higher values give finer partitioning and better locality but more
    /// partitioning overhead.  `hilbert_order = 8` → 256³ grid.
    /// Defaults to `8`.
    pub hilbert_order: u32,

    /// Size of the conflict-resolution buffer (number of points) per thread.
    ///
    /// Reserved for the long-term parallel-insertion path (module docs,
    /// idea 4); unused by the v1 optimization layer.
    /// Defaults to `65_536`.
    pub conflict_buffer_size: usize,

    /// Insertion phase selection: when `true` (default) the sequential
    /// [`Delaunay3D`] refinement pipeline inserts the volume; when `false`
    /// the [`CentroidStarMesher3D`] seed is used unrefined. Either way the
    /// parallel optimization phase runs afterwards (governed by
    /// `MeshParams::optimize_passes`).
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

// ─── Run statistics ───────────────────────────────────────────────────────────

/// Counters of one run of the HXT parallel optimization layer.
///
/// The operation counters mirror [`crate::mesh_optimize::OptimizeStats`]
/// (same semantics, same passes); the last three fields are the layer's own
/// parallelism observables.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct HxtParallelOptimizeStats {
    /// Rounds actually executed (before early stop).
    pub rounds: usize,
    /// Edge collapses accepted (one node removed per collapse).
    pub collapses: usize,
    /// Tets deleted by collapses (the edge rings).
    pub tets_removed_by_collapse: usize,
    /// 2-3 swaps accepted (net +1 tet each).
    pub face_swaps_2_3: usize,
    /// 3-2 edge flips accepted (net -1 tet each).
    pub edge_flips_3_2: usize,
    /// Nodes whose position changed in the colored relocation passes.
    pub relocations_accepted: usize,
    /// Color classes of the node-adjacency coloring in the last executed
    /// relocation pass (`0` when no relocation pass ran).
    pub colors_used: usize,
    /// Worker count of the rayon pool the layer ran on
    /// (`num_threads`, or the all-core default when `num_threads == 0`).
    pub thread_pool_size: usize,
    /// Distinct rayon workers observed executing relocation evaluations
    /// (`1` when the work stayed on one worker, `0` when no parallel
    /// evaluation ran).
    pub threads_observed: usize,
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

    /// Optimize a pure tetrahedral volume mesh with the HXT parallel
    /// optimization layer: rounds of collapse → swaps → **colored parallel
    /// relocation**, every operation acceptance-gated by `config.metric`,
    /// mirroring [`crate::mesh_optimize::MeshQualityOptimizer::run`] with the
    /// relocation pass replaced by the colored parallel one (see the module
    /// docs for the scoping and the determinism contract).
    ///
    /// `mesh_params.optimize_passes` supplies the round count (`0` = no-op);
    /// `params` supplies the convergence tolerance and the
    /// `move_boundary_nodes` policy; `config` supplies the enabled
    /// operations, the metric and the collapse threshold (pass
    /// [`OptimizeConfig::default`] for the serial optimizer's defaults).
    pub fn parallel_optimize_mesh(
        &self,
        mesh: &mut Mesh,
        config: &OptimizeConfig,
        params: &OptimizeParams,
        mesh_params: &MeshParams,
    ) -> Result<HxtParallelOptimizeStats, MeshAlgoError> {
        let rounds = mesh_params.optimize_passes as usize;
        if rounds == 0 {
            // Documented semantics: a zero round count skips optimization.
            return Ok(HxtParallelOptimizeStats::default());
        }

        // Input contract: pure tetrahedral volume mesh (mirrors the serial
        // optimizer's validation).
        if mesh.elements.is_empty() {
            return Err(MeshAlgoError::InvalidInput(
                "volume optimizer requires tetrahedral elements; mesh has none".to_string(),
            ));
        }
        for e in &mesh.elements {
            if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
                return Err(MeshAlgoError::InvalidInput(format!(
                    "volume optimizer requires a pure Tetrahedron4 mesh; element {} is {:?}",
                    e.id, e.etype
                )));
            }
        }

        let physical_names = mesh.physical_names.clone();
        let periodic = mesh.periodic.clone();

        let mut stats = HxtParallelOptimizeStats::default();
        let metric = config.metric;
        let floor_deg = (config.threshold * 0.1).clamp(0.5, 5.0);
        let observed_workers = AtomicU64::new(0);

        // Owned pool: `num_threads(0)` asks rayon for the default (all
        // logical cores).
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .map_err(|e| MeshAlgoError::Generation(format!("rayon thread pool: {e}")))?;
        stats.thread_pool_size = pool.current_num_threads();

        let mut tm = TetMesh::from_mesh(mesh);
        normalize_orientations(&mut tm)?;

        let mut collapse_threshold = config.collapse_threshold;
        if collapse_threshold <= 0.0 {
            collapse_threshold = mesh_params.min_size;
        }
        if collapse_threshold <= 0.0 {
            collapse_threshold = 0.5 * mean_edge_length(&tm);
        }

        let mut prev = mesh_worst_quality_tuple(&tm, metric);

        pool.install(|| {
            for round in 0..rounds {
                // Cumulative counters at round start, for the no-op early stop.
                let before_round = (
                    stats.face_swaps_2_3,
                    stats.edge_flips_3_2,
                    stats.relocations_accepted,
                    stats.collapses,
                );

                if config.edge_collapse {
                    let out = edge_collapse_pass(&mut tm, collapse_threshold, floor_deg);
                    stats.collapses += out.collapses;
                    stats.tets_removed_by_collapse += out.tets_removed;
                }
                if config.edge_swap {
                    stats.face_swaps_2_3 += edge_swap_pass(&mut tm, metric);
                    stats.edge_flips_3_2 += edge_removal_pass_3_2(&mut tm, metric);
                }
                if config.laplacian_smooth {
                    // Relocate to per-round stability (the 2-D module's
                    // laplaceSmoothing loop: repeat the relocation pass until
                    // no node moves, capped at 8 sweeps).
                    let free = free_nodes(&tm, params.move_boundary_nodes);
                    for _ in 0..8 {
                        let (moved, colors_used) = colored_parallel_relocate_pass(
                            &mut tm,
                            &free,
                            metric,
                            stats.thread_pool_size,
                            &observed_workers,
                        );
                        stats.relocations_accepted += moved;
                        stats.colors_used = colors_used;
                        if moved == 0 {
                            break;
                        }
                    }
                }

                stats.rounds = round + 1;

                // Early stop: no accepted operation this round, or the global
                // worst quality tuple stopped improving (primary against
                // `params.tolerance`, secondary against `tolerance * 1e-3` —
                // the same 1e-6 / 1e-9 slack pair as the tuple comparisons).
                let after_round = (
                    stats.face_swaps_2_3,
                    stats.edge_flips_3_2,
                    stats.relocations_accepted,
                    stats.collapses,
                );
                let cur = mesh_worst_quality_tuple(&tm, metric);
                let gained = (cur.0 - prev.0, cur.1 - prev.1);
                prev = cur;
                if after_round == before_round
                    || (gained.0 < params.tolerance && gained.1 < params.tolerance * 1e-3)
                {
                    break;
                }
            }
        });

        stats.threads_observed = observed_workers.load(Ordering::Relaxed).count_ones() as usize;

        let mut rebuilt = tetmesh_to_mesh_referenced(&tm);
        rebuilt.physical_names = physical_names;
        rebuilt.periodic = periodic;
        *mesh = rebuilt;
        Ok(stats)
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl Mesher3D for Hxt3D {
    fn name(&self) -> &'static str {
        "HXT Parallel Delaunay 3D"
    }

    fn mesh_3d(&self, surface: &Mesh, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        // Phase 1 — insertion (sequential delegation; the full HXT parallel
        // insertion remains the documented long-term path, module docs).
        let mut mesh = if self.enable_refinement {
            Delaunay3D::default().mesh_3d(surface, params)?
        } else {
            CentroidStarMesher3D.mesh_3d(surface, params)?
        };
        // Phase 2 — the rayon-parallel optimization layer (this card).
        self.parallel_optimize_mesh(
            &mut mesh,
            &OptimizeConfig::default(),
            &OptimizeParams::default(),
            params,
        )?;
        Ok(mesh)
    }
}

// ─── Parallel machinery ───────────────────────────────────────────────────────

/// Record this evaluation's rayon worker into the observed-workers bitmask
/// (instrumentation only — never affects results).
fn note_worker(mask: &AtomicU64) {
    if let Some(i) = rayon::current_thread_index() {
        mask.fetch_or(1u64 << (i as u64 & 63), Ordering::Relaxed);
    }
}

/// Color the tet mesh's node-adjacency graph so that no edge is
/// monochromatic: nodes of the same color share no tetrahedron, so their
/// relocation evaluations touch disjoint stars.
///
/// Greedy first-fit in ascending compact-index order — deterministic. Tet
/// meshes typically need 8-16 colors (the measured count is pinned per
/// fixture in the tests).
pub(crate) fn color_node_adjacency(tm: &TetMesh) -> Vec<u32> {
    let n = tm.nodes.len();
    // Adjacency from the tets' 6 edges each.
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
    for tet in &tm.tets {
        for i in 0..4 {
            for j in (i + 1)..4 {
                let (a, b) = (tet.nodes[i], tet.nodes[j]);
                adj[a as usize].push(b);
                adj[b as usize].push(a);
            }
        }
    }
    for list in adj.iter_mut() {
        list.sort_unstable();
        list.dedup();
    }

    let mut colors = vec![u32::MAX; n];
    for v in 0..n {
        let mut used: Vec<u32> = adj[v]
            .iter()
            .filter(|&&u| colors[u as usize] != u32::MAX)
            .map(|&u| colors[u as usize])
            .collect();
        used.sort_unstable();
        used.dedup();
        let mut c = 0u32;
        for &u in &used {
            if u == c {
                c += 1;
            } else {
                break;
            }
        }
        colors[v] = c;
    }
    colors
}

/// Colored parallel relocation pass: the HXT replacement for the serial
/// [`crate::mesh_optimize`] `relocate_pass`.
///
/// Per color class (ascending color id), all candidate relocations are
/// evaluated in parallel with rayon (each evaluation is the pure
/// [`relocate_candidate_for_node`] kernel over the pass's star snapshot),
/// then accepted moves are committed sequentially in ascending node order —
/// the v1 determinism contract (module docs). Within a class no two nodes
/// share a tet, so a commit cannot invalidate another candidate's gates.
///
/// Returns `(moved_nodes, colors_used)`.
fn colored_parallel_relocate_pass(
    tm: &mut TetMesh,
    free: &[u32],
    metric: QualityMetric,
    pool_threads: usize,
    observed_workers: &AtomicU64,
) -> (usize, usize) {
    let colors = color_node_adjacency(tm);
    let n_colors = colors.iter().copied().max().map_or(0, |c| c + 1) as usize;
    // Group the free nodes (already ascending) by color; the ascending order
    // is preserved inside every group — that is the commit order.
    let mut groups: Vec<Vec<u32>> = vec![Vec::new(); n_colors];
    for &v in free {
        groups[colors[v as usize] as usize].push(v);
    }

    // Star topology snapshot — valid for the whole pass (nodes move, the
    // connectivity does not).
    let star = star_tets(tm);

    let mut moved = 0usize;
    // Chunk granularity: aim for ~8 chunks per pool worker so every worker
    // has stealing opportunities at any mesh size, down to single-node
    // chunks on the small fixtures. Purely a scheduling knob — `collect`
    // preserves element order regardless of chunking, so results are
    // bit-identical for any value.
    let min_len = (free.len() / (pool_threads * 8)).max(1);
    for group in &groups {
        // Parallel evaluation against the current snapshot.
        let candidates: Vec<Option<[f64; 3]>> = group
            .par_iter()
            .with_min_len(min_len)
            .map(|&v| {
                note_worker(observed_workers);
                relocate_candidate_for_node(tm, &star, v, metric)
            })
            .collect();
        // Sequential commit in ascending node order.
        for (&v, candidate) in group.iter().zip(candidates) {
            if let Some(p) = candidate {
                tm.nodes[v as usize] = p;
                moved += 1;
            }
        }
    }
    (moved, groups.len())
}

/// [`mesh_quality_stats`] computed with rayon over the tets.
///
/// The expensive per-tet measures (6-edge dihedral fold, circumsphere
/// solve) run in parallel; the reduction is restricted to
/// order-insensitive operations (`min`) plus an element-order ratio vector,
/// so the result is bit-identical to the serial [`mesh_quality_stats`]
/// (pinned by test).
pub fn parallel_mesh_quality_stats(mesh: &Mesh) -> TetMeshQualityStats {
    // Read-only coordinate table shared by all workers.
    let coords: HashMap<u64, [f64; 3]> = mesh
        .nodes
        .iter()
        .map(|(id, n)| (*id, [n.position.x, n.position.y, n.position.z]))
        .collect();

    // Parallel per-tet evaluation; `map` preserves element order.
    let per_tet: Vec<Option<(f64, f64)>> = mesh
        .elements
        .par_iter()
        .map(|e| {
            if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
                return None;
            }
            let mut pts = [[0.0f64; 3]; 4];
            for (k, id) in e.node_ids.iter().enumerate() {
                pts[k] = *coords.get(id)?;
            }
            if signed_volume_tet(pts[0], pts[1], pts[2], pts[3]).abs() < VOL_EPS {
                return None;
            }
            let dihedral = min_dihedral_deg_tet(pts[0], pts[1], pts[2], pts[3]);
            let ratio = radius_edge_ratio_tet(pts[0], pts[1], pts[2], pts[3]);
            Some((dihedral, ratio))
        })
        .collect();

    // Deterministic reduction: p95 over the element-order ratio vector,
    // counts by a sequential fold, min by an order-insensitive parallel
    // reduce.
    let mut ratios = Vec::with_capacity(per_tet.len());
    let mut slivers = 0usize;
    for r in per_tet.iter().flatten() {
        ratios.push(r.1);
        if is_sliver_tet(r.0, r.1) {
            slivers += 1;
        }
    }
    let total = ratios.len();
    let min_dihedral = per_tet
        .par_iter()
        .filter_map(|r| r.as_ref())
        .map(|r| r.0)
        .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or(f64::INFINITY);

    ratios.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let p95 = if ratios.is_empty() {
        f64::INFINITY
    } else {
        let idx = ((ratios.len() - 1) as f64 * 0.95).round() as usize;
        ratios[idx]
    };

    TetMeshQualityStats {
        min_dihedral_deg: if min_dihedral.is_finite() {
            min_dihedral
        } else {
            0.0
        },
        p95_radius_edge: p95,
        sliver_fraction: if total == 0 {
            1.0
        } else {
            slivers as f64 / total as f64
        },
        tet_count: total,
    }
}

// ─── Internal helpers (scheme-level, long-term insertion path) ────────────────

/// Compute the 3-D Hilbert curve index for a point in the unit cube.
///
/// `order` controls the recursion depth (grid of `2^order` cells per axis).
///
/// Returns a 64-bit integer key that can be used to sort points for
/// cache-optimal insertion order.
#[allow(dead_code)]
fn hilbert_index_3d(x: f64, y: f64, z: f64, order: u32) -> u64 {
    let n = 1u64 << order.min(20);
    let clamp = |v: f64| (v.clamp(0.0, 1.0 - f64::EPSILON) * n as f64) as u64;
    let (ix, iy, iz) = (clamp(x), clamp(y), clamp(z));
    let mut code = 0u64;
    for bit in 0..order.min(20) {
        code |= ((ix >> bit) & 1) << (3 * bit);
        code |= ((iy >> bit) & 1) << (3 * bit + 1);
        code |= ((iz >> bit) & 1) << (3 * bit + 2);
    }
    code
}

/// Assign each 3-D grid cell a color such that no two adjacent cells (sharing
/// a face, edge, or corner) have the same color.
///
/// For a 3-D grid the chromatic number is 8 (2×2×2 checkerboard in 3-D).
/// Returns a `Vec<u8>` of length `nx * ny * nz` with values in `0..8`.
#[allow(dead_code)]
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

/// A lightweight atomic-CAS–based ownership token for a tetrahedron.
///
/// During parallel insertion each thread uses this to claim tetrahedra
/// before modifying them.  If the CAS fails (another thread owns the tet),
/// the point is added to the conflict buffer for sequential re-insertion.
#[allow(dead_code)]
struct TetOwnership {
    /// Owning thread ID (0 = free).
    owner: std::sync::atomic::AtomicUsize,
}

#[allow(dead_code)]
impl TetOwnership {
    fn new() -> Self {
        Self {
            owner: std::sync::atomic::AtomicUsize::new(0),
        }
    }

    /// Try to claim this tetrahedron for thread `thread_id + 1`.
    ///
    /// Returns `true` on success, `false` if another thread owns it.
    fn try_claim(&self, thread_id: usize) -> bool {
        self.owner
            .compare_exchange(
                0,
                thread_id + 1,
                std::sync::atomic::Ordering::Acquire,
                std::sync::atomic::Ordering::Relaxed,
            )
            .is_ok()
    }

    /// Release ownership.
    fn release(&self) {
        self.owner.store(0, std::sync::atomic::Ordering::Release);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::delaunay_3d::Delaunay3D;
    use crate::mesh_optimize::{mesh_quality_stats, signed_volume_tet, MeshQualityOptimizer};
    use crate::traits::Mesher3D;
    use rmsh_model::{Element, ElementType, Mesh, Node};

    fn box_surface(lx: f64, ly: f64, lz: f64) -> Mesh {
        let mut mesh = Mesh::new();
        for (id, xyz) in [
            (1, [0.0, 0.0, 0.0]),
            (2, [lx, 0.0, 0.0]),
            (3, [lx, ly, 0.0]),
            (4, [0.0, ly, 0.0]),
            (5, [0.0, 0.0, lz]),
            (6, [lx, 0.0, lz]),
            (7, [lx, ly, lz]),
            (8, [0.0, ly, lz]),
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

    fn cube_surface() -> Mesh {
        box_surface(1.0, 1.0, 1.0)
    }

    /// Mesh the closed surface with the sequential Delaunay pipeline (the
    /// insertion phase Hxt3D delegates to) — the shared input for the
    /// serial-vs-parallel acceptance comparisons.
    fn delaunay_volume(surface: &Mesh, size: f64) -> Mesh {
        Delaunay3D::default()
            .mesh_3d(surface, &MeshParams::with_size(size))
            .unwrap()
    }

    /// Bit-exact mesh equality: node coordinate bits (sorted by node id) and
    /// the element connectivity sequence. The v1 determinism contract of the
    /// parallel optimization layer is that same input + same config produce
    /// bit-identical output for any `num_threads`.
    fn meshes_bit_identical(a: &Mesh, b: &Mesh) -> bool {
        if a.nodes.len() != b.nodes.len() || a.elements.len() != b.elements.len() {
            return false;
        }
        let mut na: Vec<_> = a
            .nodes
            .iter()
            .map(|(id, n)| {
                (
                    *id,
                    n.position.x.to_bits(),
                    n.position.y.to_bits(),
                    n.position.z.to_bits(),
                )
            })
            .collect();
        na.sort_unstable_by_key(|e| e.0);
        let mut nb: Vec<_> = b
            .nodes
            .iter()
            .map(|(id, n)| {
                (
                    *id,
                    n.position.x.to_bits(),
                    n.position.y.to_bits(),
                    n.position.z.to_bits(),
                )
            })
            .collect();
        nb.sort_unstable_by_key(|e| e.0);
        if na != nb {
            return false;
        }
        a.elements
            .iter()
            .zip(&b.elements)
            .all(|(ea, eb)| ea.etype == eb.etype && ea.node_ids == eb.node_ids)
    }

    fn tet_coords(mesh: &Mesh, ids: &[u64; 4]) -> Option<[[f64; 3]; 4]> {
        let mut pts = [[0.0; 3]; 4];
        for (k, id) in ids.iter().enumerate() {
            let node = mesh.nodes.get(id)?;
            pts[k] = [node.position.x, node.position.y, node.position.z];
        }
        Some(pts)
    }

    /// Domain-preservation invariant: total tet volume.
    fn total_tet_volume(mesh: &Mesh) -> f64 {
        mesh.elements
            .iter()
            .filter(|e| e.etype == ElementType::Tetrahedron4 && e.node_ids.len() == 4)
            .filter_map(|e| {
                let ids = [e.node_ids[0], e.node_ids[1], e.node_ids[2], e.node_ids[3]];
                tet_coords(mesh, &ids)
            })
            .map(|p| signed_volume_tet(p[0], p[1], p[2], p[3]).abs())
            .sum()
    }

    #[test]
    fn hilbert_index_progresses_along_diagonal() {
        assert!(hilbert_index_3d(0.2, 0.2, 0.2, 8) > hilbert_index_3d(0.1, 0.1, 0.1, 8));
    }

    #[test]
    fn hxt_3d_generates_mesh() {
        let mesh = Hxt3D::default()
            .mesh_3d(&cube_surface(), &MeshParams::with_size(0.4))
            .unwrap();
        assert!(mesh.elements_by_dimension(3).len() > 0);
    }

    // ── HXT rayon-parallel optimization layer (P4-tail) ──────────────────────

    #[test]
    fn node_coloring_is_valid_deterministic_and_pinned_on_box_fixture() {
        let surface = cube_surface();
        let volume = delaunay_volume(&surface, 0.4);
        let tm = TetMesh::from_mesh(&volume);
        let colors = color_node_adjacency(&tm);
        assert_eq!(colors.len(), tm.nodes.len());
        // Validity: endpoints of every tet edge carry different colors, so a
        // color class is an independent set — relocations inside one class
        // touch disjoint stars.
        for tet in &tm.tets {
            for i in 0..4 {
                for j in (i + 1)..4 {
                    assert_ne!(
                        colors[tet.nodes[i] as usize],
                        colors[tet.nodes[j] as usize],
                        "edge ({}, {}) is monochromatic",
                        tet.nodes[i],
                        tet.nodes[j]
                    );
                }
            }
        }
        // Determinism of the coloring itself.
        assert_eq!(colors, color_node_adjacency(&tm));
        let used = colors.iter().copied().max().unwrap_or(0) + 1;
        eprintln!("box fixture node-adjacency colors: {used}");
        // Measured color-count pin (greedy first-fit on the box fixture):
        // 7 colors for the Delaunay box volume at size 0.4 (re-measured
        // when the 3-D quality card wired the flip-convergence arm into
        // Delaunay3D — the seed volume changed, the adjacency graph with
        // it; the coloring stays valid + deterministic, this pin only
        // records the new measured count).
        assert_eq!(used, 7, "measured box-fixture color count drifted");
    }

    #[test]
    fn parallel_quality_stats_bit_match_serial() {
        for (name, (lx, ly, lz), size) in [
            ("cube", (1.0, 1.0, 1.0), 0.4),
            ("stretched", (1.8, 1.1, 0.9), 0.38),
        ] {
            let volume = delaunay_volume(&box_surface(lx, ly, lz), size);
            let serial = mesh_quality_stats(&volume);
            let par = parallel_mesh_quality_stats(&volume);
            eprintln!(
                "{name}: serial(min={:.6}, p95={:.6}, sliver={:.6}, n={}) parallel(min={:.6}, p95={:.6}, sliver={:.6}, n={})",
                serial.min_dihedral_deg,
                serial.p95_radius_edge,
                serial.sliver_fraction,
                serial.tet_count,
                par.min_dihedral_deg,
                par.p95_radius_edge,
                par.sliver_fraction,
                par.tet_count,
            );
            assert_eq!(serial.tet_count, par.tet_count);
            assert_eq!(
                serial.min_dihedral_deg.to_bits(),
                par.min_dihedral_deg.to_bits()
            );
            assert_eq!(
                serial.p95_radius_edge.to_bits(),
                par.p95_radius_edge.to_bits()
            );
            assert_eq!(
                serial.sliver_fraction.to_bits(),
                par.sliver_fraction.to_bits()
            );
        }
    }

    /// Parity acceptance: the parallel pipeline matches or beats the serial
    /// optimizer on the box/stretched/slender fixture family.
    ///
    /// Measured contract (bit-deterministic pipeline, values stable):
    ///
    /// * **box** (input min 5.625, p95 20.629): serial(min 20.611, p95
    ///   6.107, sliver 0) vs parallel(min 19.777, p95 5.317, sliver 0). The
    ///   parallel pipeline **beats** the radius-edge axis (5.317 < 6.107 —
    ///   the charter gate axis) and degrades the worst single dihedral by
    ///   0.83 deg (4%) — the measured cost of the batch-commit independence
    ///   (color-class ordering differs from the serial sweep's index
    ///   ordering, so the descent reaches a different local optimum). Pinned
    ///   as `>= serial - 1.0 deg`; both paths lift the input min from 5.6
    ///   deg by more than 3x.
    /// * **stretched** (input min 0.176, p95 1204.910): all three axes
    ///   bit-equal to serial (min 0.176, p95 1204.910, sliver 0.2419).
    /// * **slender** (n=97 input): min and p95 bit-equal to serial
    ///   (0.069 deg / 56675.829); sliver *count* identical (50), the
    ///   fraction reads +0.005 only because the parallel path accepted one
    ///   extra collapse (n 96 vs 97, removing a non-sliver tet).
    #[test]
    fn parallel_pipeline_matches_serial_optimizer_on_box() {
        let (q_in, q_ser, q_par, stats) = serial_vs_parallel(1.0, 1.0, 1.0, 0.4);
        eprintln!(
            "box parity: input(min={:.3}, p95={:.3}) serial(min={:.3}, p95={:.3}, sliver={:.4}, n={}) parallel(min={:.3}, p95={:.3}, sliver={:.4}, n={}) stats={stats:?}",
            q_in.min_dihedral_deg, q_in.p95_radius_edge,
            q_ser.min_dihedral_deg, q_ser.p95_radius_edge, q_ser.sliver_fraction, q_ser.tet_count,
            q_par.min_dihedral_deg, q_par.p95_radius_edge, q_par.sliver_fraction, q_par.tet_count
        );
        // Non-degradation vs the un-optimized input (the serial optimizer's
        // own pinning style), on both axes.
        assert!(q_par.min_dihedral_deg >= q_in.min_dihedral_deg);
        assert!(q_par.p95_radius_edge <= q_in.p95_radius_edge);
        // Radius-edge (charter gate axis): strictly beats serial.
        assert!(
            q_par.p95_radius_edge <= q_ser.p95_radius_edge,
            "box: parallel p95 {:.3} > serial {:.3}",
            q_par.p95_radius_edge,
            q_ser.p95_radius_edge
        );
        assert!(
            (q_par.p95_radius_edge - 5.317).abs() < 0.05,
            "box: pinned parallel p95 drifted to {}",
            q_par.p95_radius_edge
        );
        // Min dihedral: within the measured ordering-cost band of serial
        // (measured delta -0.834 deg) and pinned near the measured value.
        assert!(
            q_par.min_dihedral_deg >= q_ser.min_dihedral_deg - 1.0,
            "box: parallel min dihedral {:.3} fell >1 deg below serial {:.3}",
            q_par.min_dihedral_deg,
            q_ser.min_dihedral_deg
        );
        assert!(
            (q_par.min_dihedral_deg - 19.777).abs() < 0.05,
            "box: pinned parallel min dihedral drifted to {}",
            q_par.min_dihedral_deg
        );
        assert_eq!(q_par.sliver_fraction, 0.0);
    }

    #[test]
    fn parallel_pipeline_matches_serial_optimizer_on_stretched_box() {
        let (q_in, q_ser, q_par, stats) = serial_vs_parallel(1.8, 1.1, 0.9, 0.38);
        eprintln!(
            "stretched parity: input(min={:.3}, p95={:.3}) serial(min={:.3}, p95={:.3}, sliver={:.4}, n={}) parallel(min={:.3}, p95={:.3}, sliver={:.4}, n={}) stats={stats:?}",
            q_in.min_dihedral_deg, q_in.p95_radius_edge,
            q_ser.min_dihedral_deg, q_ser.p95_radius_edge, q_ser.sliver_fraction, q_ser.tet_count,
            q_par.min_dihedral_deg, q_par.p95_radius_edge, q_par.sliver_fraction, q_par.tet_count
        );
        // Bit-equal to serial on all three axes (measured).
        assert!(q_par.min_dihedral_deg >= q_ser.min_dihedral_deg);
        assert!(q_par.min_dihedral_deg >= q_in.min_dihedral_deg);
        assert!(q_par.p95_radius_edge <= q_ser.p95_radius_edge);
        assert!(q_par.p95_radius_edge <= q_in.p95_radius_edge);
        assert!(q_par.sliver_fraction <= q_ser.sliver_fraction);
        assert!(
            (q_par.p95_radius_edge - 1204.910).abs() < 1.0,
            "stretched: pinned parallel p95 drifted to {}",
            q_par.p95_radius_edge
        );
    }

    #[test]
    fn parallel_pipeline_matches_serial_optimizer_on_slender_box() {
        let (q_in, q_ser, q_par, stats) = serial_vs_parallel(3.2, 0.55, 0.45, 0.26);
        eprintln!(
            "slender parity: input(min={:.3}, p95={:.3}) serial(min={:.3}, p95={:.3}, sliver={:.4}, n={}) parallel(min={:.3}, p95={:.3}, sliver={:.4}, n={}) stats={stats:?}",
            q_in.min_dihedral_deg, q_in.p95_radius_edge,
            q_ser.min_dihedral_deg, q_ser.p95_radius_edge, q_ser.sliver_fraction, q_ser.tet_count,
            q_par.min_dihedral_deg, q_par.p95_radius_edge, q_par.sliver_fraction, q_par.tet_count
        );
        // Min and radius-edge bit-equal to serial (measured).
        assert!(q_par.min_dihedral_deg >= q_ser.min_dihedral_deg - 1e-9);
        assert!(q_par.p95_radius_edge <= q_ser.p95_radius_edge);
        assert!(
            (q_par.p95_radius_edge - 56675.829).abs() < 1.0,
            "slender: pinned parallel p95 drifted to {}",
            q_par.p95_radius_edge
        );
        // Sliver count identical (50 == 50); the fraction reads +0.005
        // because the parallel path accepted one extra collapse (n 96 vs
        // 97) whose removed ring tet was not a sliver.
        let ser_count = (q_ser.sliver_fraction * q_ser.tet_count as f64).round() as usize;
        let par_count = (q_par.sliver_fraction * q_par.tet_count as f64).round() as usize;
        assert!(
            par_count <= ser_count,
            "slender: parallel sliver count {par_count} > serial {ser_count}"
        );
        assert!(q_par.sliver_fraction <= q_ser.sliver_fraction + 0.01);
    }

    #[test]
    fn parallel_optimize_is_bit_deterministic_across_thread_counts() {
        let surface = box_surface(1.8, 1.1, 0.9);
        let params = MeshParams::with_size(0.38);
        let volume = delaunay_volume(&surface, 0.38);

        let mut run_a = volume.clone();
        Hxt3D::default()
            .parallel_optimize_mesh(
                &mut run_a,
                &OptimizeConfig::default(),
                &OptimizeParams::default(),
                &params,
            )
            .unwrap();
        let mut run_b = volume.clone();
        Hxt3D::default()
            .parallel_optimize_mesh(
                &mut run_b,
                &OptimizeConfig::default(),
                &OptimizeParams::default(),
                &params,
            )
            .unwrap();
        let mut run_c = volume.clone();
        Hxt3D::default()
            .single_threaded()
            .parallel_optimize_mesh(
                &mut run_c,
                &OptimizeConfig::default(),
                &OptimizeParams::default(),
                &params,
            )
            .unwrap();

        assert!(meshes_bit_identical(&run_a, &run_b), "two all-core runs differ");
        assert!(
            meshes_bit_identical(&run_a, &run_c),
            "single-thread run differs from all-core run"
        );
    }

    #[test]
    fn parallel_pipeline_exercises_a_forced_thread_pool() {
        let surface = box_surface(3.2, 0.55, 0.45);
        let params = MeshParams::with_size(0.26);
        let mut mesh = delaunay_volume(&surface, 0.26);
        let hxt = Hxt3D {
            num_threads: 4,
            ..Default::default()
        };
        let stats = hxt
            .parallel_optimize_mesh(
                &mut mesh,
                &OptimizeConfig::default(),
                &OptimizeParams::default(),
                &params,
            )
            .unwrap();
        eprintln!(
            "threading: pool={}, observed_workers={}",
            stats.thread_pool_size, stats.threads_observed
        );
        // Deterministic pin: the layer ran on a pool of the requested size.
        assert_eq!(stats.thread_pool_size, 4);
        // Observational pin: more than one rayon worker executed relocation
        // evaluations (repeated par_iter batches over the free-node set).
        assert!(
            stats.threads_observed >= 2,
            "only {} worker thread(s) observed — parallel path did not engage",
            stats.threads_observed
        );
    }

    #[test]
    fn hxt_mesh_3d_is_bit_deterministic_end_to_end() {
        let surface = cube_surface();
        let params = MeshParams::with_size(0.4);
        let a = Hxt3D::default().mesh_3d(&surface, &params).unwrap();
        let b = Hxt3D::default().mesh_3d(&surface, &params).unwrap();
        assert!(meshes_bit_identical(&a, &b));
    }

    /// Shared serial-vs-parallel comparison: mesh the box, optimize one clone
    /// with the serial [`MeshQualityOptimizer`], another with the HXT
    /// parallel pipeline, return the input and both optimized quality stats.
    fn serial_vs_parallel(
        lx: f64,
        ly: f64,
        lz: f64,
        size: f64,
    ) -> (
        crate::mesh_optimize::TetMeshQualityStats,
        crate::mesh_optimize::TetMeshQualityStats,
        crate::mesh_optimize::TetMeshQualityStats,
        HxtParallelOptimizeStats,
    ) {
        let surface = box_surface(lx, ly, lz);
        let params = MeshParams::with_size(size);
        let volume = delaunay_volume(&surface, size);
        let q_in = mesh_quality_stats(&volume);

        let mut serial = volume.clone();
        MeshQualityOptimizer::new()
            .optimize_with_params(&mut serial, &OptimizeParams::default(), &params)
            .unwrap();
        let q_ser = mesh_quality_stats(&serial);

        let mut par = volume.clone();
        let stats = Hxt3D::default()
            .parallel_optimize_mesh(
                &mut par,
                &OptimizeConfig::default(),
                &OptimizeParams::default(),
                &params,
            )
            .unwrap();
        let q_par = mesh_quality_stats(&par);

        // Domain preservation on both paths.
        assert!((total_tet_volume(&serial) - lx * ly * lz).abs() < 1e-6);
        assert!((total_tet_volume(&par) - lx * ly * lz).abs() < 1e-6);
        (q_in, q_ser, q_par, stats)
    }
}
