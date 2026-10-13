//! The a posteriori adaptive loop (P6): error estimate → size-field
//! update → remesh, closed over the existing pieces — the P2b size-field
//! framework consumes the produced field and Frontal-Delaunay-2D
//! (`with_size_field`) re-meshes from it.
//!
//! # The surface
//!
//! - [`ErrorEstimator`]: per-element error indicators η_K for a mesh and
//!   a nodal solution u_h (the caller supplies the solution — an analytic
//!   reference, an interpolated computed field, anything nodal).
//! - [`GradientRecovery`]: the concrete Zienkiewicz-Zhu-style estimator —
//!   per-element constant gradients, area-weighted nodal recovery, and
//!   η_K = h_K · sqrt(mean over the element's nodes of |∇rec − ∇_K|²).
//!   A linear field recovers exactly (indicator exactly zero).
//! - [`AdaptiveSizeField`]: the indicator → target-size map
//!   (h_new = h_old · tol/η, clamped to [h_min, h_max], nodally averaged
//!   over incident elements), answering `SizeField::size_at` through
//!   nearest-node interpolation on the OLD mesh — the standard
//!   piecewise-constant transfer to the next mesh.
//! - [`AdaptiveLoop2D`]: the driver — mesh, estimate, update, re-mesh,
//!   until the max indicator drops under the tolerance or max_steps is
//!   reached. Deterministic: every step is a pure function of the
//!   previous mesh.
//!
//! # Fidelity lines (v1)
//!
//! - 2-D only (the Frontal-Delaunay-2D route); the 3-D loop is a
//!   follow-up once a 3-D size-field mesher consumes the same trait.
//! - The solution is provided nodally per iteration (a closure); a
//!   persistent discrete solution transported across meshes is the
//!   caller's business, not the loop's.
//!
//! # Empirical loop contracts (probed, not assumed)
//!
//! - **The ZZ indicator is not a monotone convergence gauge at a sharp
//!   front.** While h enters the front's feature scale, the band-mean
//!   indicator GROWS first (the estimator sees more of the front as the
//!   elements resolve it) and only then falls. The global max carries
//!   boundary-anchored noise on top (the P3 lesson: recovery at boundary
//!   nodes averages fewer samples). The honest loop contract is therefore
//!   GEOMETRIC CONCENTRATION — the feature band's triangle count grows
//!   and the band's mean triangle size drops under half the far field's —
//!   not indicator monotonicity; the loop terminates on
//!   `max_indicator <= error_tolerance` OR `max_steps`, and tests pin the
//!   geometric contract.
//! - **The solution order is part of the estimator's contract.** `u_h`
//!   values are consumed in dense sorted-node-id order (the order
//!   [`AdaptiveLoop2D::run`] samples in, and the order direct callers
//!   must build). The id -> index map is built once per call with the
//!   keys sorted first — a raw HashMap-order enumerate scrambles the
//!   solution across nodes, every indicator becomes noise, the derived
//!   size field collapses to `h_min` domain-wide, and the next mesh
//!   explodes (32k seeds on a domain that needs 70). The linear-field
//!   exactness test doubles as the regression pin for the ordering.
//! - **The loop's `MeshParams` must carry the clamp range** `[h_min,
//!   h_max]` as `min_size`/`max_size`: the mesher clamps every size-field
//!   answer into `[min_size, max_size]`, so `element_size: h_max` alone
//!   would silently floor the refinement at `h_max / 10` (the mesher's
//!   uniform default floor) instead of the loop's `h_min`.
//! - **Performance line (v1):** a loop step is a full re-mesh, and the
//!   mesher's per-iteration structures (the `worst_active` scan, the
//!   frozen size-snapshot clones, the relocation passes) are quadratic
//!   overall in debug builds — invisible at the P3a mesh sizes, but a
//!   loop step on a few-thousand-triangle mesh costs tens of seconds.
//!   Fine for the interface card; a mesher perf card if the loop ever
//!   lands on a hot path.

use std::collections::HashMap;

use glam::DVec3;

use rmsh_model::{Element, ElementType, Mesh};

use crate::size_field::SizeField;
use crate::frontal_delaunay_2d::FrontalDelaunay2D;
use crate::traits::{MeshParams, Mesher2D};

// ─── Errors ───────────────────────────────────────────────────────────────────

/// Every way the adaptive loop can refuse or fail.
#[derive(Debug, thiserror::Error)]
pub enum AdaptiveError {
    #[error("solution length {len} != node count {expected}")]
    SolutionMismatch { len: usize, expected: usize },
    #[error("error_tolerance must be finite and positive, got {value}")]
    InvalidTolerance { value: f64 },
    #[error("invalid clamp range: h_min={h_min}, h_max={h_max}")]
    InvalidClamp { h_min: f64, h_max: f64 },
    #[error("meshing failed at step {step}: {detail}")]
    MeshingFailed {
        step: usize,
        detail: String,
    },
}

// ─── The estimator trait ──────────────────────────────────────────────────────

/// Per-element error indicators for a nodal solution on a mesh. The
/// indicators are relative magnitudes — the loop only consumes their
/// ordering and their ratio against the tolerance.
pub trait ErrorEstimator {
    fn indicators(&self, mesh: &Mesh, u_h: &[f64]) -> Result<Vec<f64>, AdaptiveError>;
}

// ─── Gradient recovery (Zienkiewicz-Zhu form) ────────────────────────────────

/// The recovered-gradient estimator: per-element constant gradients from
/// the nodal P1 field, area-weighted nodal recovery, and
/// η_K = h_K · sqrt(mean over the element's nodes of |∇rec_i − ∇_K|²).
/// A field that is linear on an element recovers exactly (η_K = 0).
#[derive(Debug, Clone, Copy, Default)]
pub struct GradientRecovery;

impl ErrorEstimator for GradientRecovery {
    fn indicators(&self, mesh: &Mesh, u_h: &[f64]) -> Result<Vec<f64>, AdaptiveError> {
        if u_h.len() != mesh.nodes.len() {
            return Err(AdaptiveError::SolutionMismatch {
                len: u_h.len(),
                expected: mesh.nodes.len(),
            });
        }
        let index = build_index(mesh);
        // element gradients (P1: constant per triangle)
        struct ElemGrad {
            nodes: [u64; 3],
            grad: [f64; 2],
            area: f64,
            h: f64,
        }
        let mut elems: Vec<ElemGrad> = Vec::new();
        let mut node_sizes: HashMap<u64, f64> = HashMap::new();
        for element in &mesh.elements {
            if element.etype != ElementType::Triangle3 {
                continue; // the estimator consumes the triangle skeleton
            }
            let ids = [
                element.node_ids[0],
                element.node_ids[1],
                element.node_ids[2],
            ];
            let p = [
                mesh.nodes.get(&ids[0]).expect("node").position,
                mesh.nodes.get(&ids[1]).expect("node").position,
                mesh.nodes.get(&ids[2]).expect("node").position,
            ];
            let ux = p[1].x - p[0].x;
            let uy = p[1].y - p[0].y;
            let vx = p[2].x - p[0].x;
            let vy = p[2].y - p[0].y;
            let det = ux * vy - uy * vx;
            if det.abs() <= f64::EPSILON {
                continue; // degenerate triangle: no gradient, no indicator
            }
            let area = det.abs() * 0.5;
            let du = [
                u_h[index[&ids[1]]] - u_h[index[&ids[0]]],
                u_h[index[&ids[2]]] - u_h[index[&ids[0]]],
            ];
            // [gx gy]ᵀ = (1/det) [vy -uy; -vx ux] [du0 du1]ᵀ
            let grad = [
                (vy * du[0] - uy * du[1]) / det,
                (-vx * du[0] + ux * du[1]) / det,
            ];
            let h = area.max(1e-300).sqrt() * 2.0f64.sqrt(); // characteristic edge scale
            elems.push(ElemGrad {
                nodes: ids,
                grad,
                area,
                h,
            });
            // element "size" for the recovery area weights and the size map
            let h_elem = 2.0 * area.sqrt();
            for id in ids {
                *node_sizes.entry(id).or_insert(0.0) += h_elem * area;
            }
        }

        // nodal recovered gradients: area-weighted average
        let mut rec: HashMap<u64, [f64; 2]> = HashMap::new();
        let mut weight: HashMap<u64, f64> = HashMap::new();
        for e in &elems {
            for &id in &e.nodes {
                let entry = rec.entry(id).or_insert([0.0, 0.0]);
                entry[0] += e.grad[0] * e.area;
                entry[1] += e.grad[1] * e.area;
                *weight.entry(id).or_insert(0.0) += e.area;
            }
        }
        for (id, g) in rec.iter_mut() {
            let w = weight[id];
            if w > 0.0 {
                g[0] /= w;
                g[1] /= w;
            }
        }

        // indicators
        let mut indicators = Vec::new();
        for e in &elems {
            let mut acc = 0.0;
            for &id in &e.nodes {
                if let Some(g) = rec.get(&id) {
                    let d0 = g[0] - e.grad[0];
                    let d1 = g[1] - e.grad[1];
                    acc += d0 * d0 + d1 * d1;
                }
            }
            let eta = e.h * (acc / 3.0).sqrt();
            indicators.push(eta);
        }
        Ok(indicators)
    }
}

/// The node-index for a nodal-solution lookup: the solution array is
/// indexed by the DENSE node order (sorted ids), the same convention the
/// loop driver uses when it samples `u_h` (and the one this module
/// documents for direct callers). The map is built once per call (a
/// per-lookup scan would be quadratic); the keys are SORTED first — a
/// raw `keys().enumerate()` walks the HashMap's internal order, which is
/// neither the sorted order nor stable across runs, and a mismatch
/// scrambles the solution across nodes (every indicator becomes noise,
/// the derived size field collapses to h_min domain-wide, and the next
/// mesh explodes).
fn build_index(mesh: &Mesh) -> HashMap<u64, usize> {
    let mut ids: Vec<u64> = mesh.nodes.keys().copied().collect();
    ids.sort_unstable();
    ids.into_iter()
        .enumerate()
        .map(|(i, k)| (k, i))
        .collect()
}

// ─── The indicator → size map ─────────────────────────────────────────────────

/// The target-size field built from indicators: per-node
/// h_new = h_old · tol/η (clamped), nodally averaged over incident
/// elements, answered through nearest-node interpolation.
#[derive(Debug)]
pub struct AdaptiveSizeField {
    /// (position, size) per retained node of the old mesh.
    nodes: Vec<(DVec3, f64)>,
    h_min: f64,
    h_max: f64,
}

impl AdaptiveSizeField {
    /// Build the size map: `h_new = h_old · (tol/η_K)^(1/p)` per element
    /// (p = 1, the P1 degree), nodally averaged over incident elements,
    /// clamped to `[h_min, h_max]`.
    pub fn build(
        mesh: &Mesh,
        indicators: &[f64],
        tolerance: f64,
        h_min: f64,
        h_max: f64,
    ) -> Result<Self, AdaptiveError> {
        if !(tolerance.is_finite() && tolerance > 0.0) {
            return Err(AdaptiveError::InvalidTolerance { value: tolerance });
        }
        if !(h_min.is_finite() && h_max.is_finite() && h_min > 0.0 && h_min <= h_max) {
            return Err(AdaptiveError::InvalidClamp { h_min, h_max });
        }
        // per-element target sizes
        let mut elem_target: HashMap<usize, f64> = HashMap::new();
        let mut k = 0usize;
        for element in &mesh.elements {
            if element.etype != ElementType::Triangle3 {
                continue;
            }
            if let Some(&eta) = indicators.get(k) {
                let ratio = (tolerance / eta.max(1e-300)).powf(1.0);
                elem_target.insert(k, ratio);
            }
            k += 1;
        }
        // nodal value = MIN over the incident elements' targets: a node
        // touched by any element that wants fine sizes stays fine — the
        // refinement-preserving choice (an average lets one large far-field
        // neighbor wash the signal out on coarse meshes)
        let mut tri_index = 0usize;
        let mut node_min: HashMap<u64, f64> = HashMap::new();
        for element in &mesh.elements {
            if element.etype != ElementType::Triangle3 {
                continue;
            }
            let h_old = triangle_size(mesh, element);
            let target = elem_target.get(&tri_index).copied().unwrap_or(1.0);
            let h_new = (h_old * target).clamp(h_min, h_max);
            for &id in &element.node_ids {
                let entry = node_min.entry(id).or_insert(h_max);
                *entry = (*entry).min(h_new);
            }
            tri_index += 1;
        }
        // deterministic retention order (sorted by node id — a HashMap walk
        // would make the nearest-node ties non-reproducible)
        let mut node_ids: Vec<u64> = mesh.nodes.keys().copied().collect();
        node_ids.sort_unstable();
        let mut nodes = Vec::new();
        for id in &node_ids {
            if let Some(h) = node_min.get(id) {
                let node = &mesh.nodes[id];
                nodes.push((
                    DVec3::new(node.position.x, node.position.y, node.position.z),
                    h.clamp(h_min, h_max),
                ));
            }
        }
        Ok(AdaptiveSizeField {
            nodes,
            h_min,
            h_max,
        })
    }
}

impl SizeField for AdaptiveSizeField {
    fn size_at(&self, p: DVec3) -> f64 {
        // nearest-node interpolation on the old mesh (documented v1
        // transfer); ties break toward the first node in retention order
        let mut best: Option<(f64, f64)> = None;
        for (pos, h) in &self.nodes {
            let dx = pos.x - p.x;
            let dy = pos.y - p.y;
            let dz = pos.z - p.z;
            let d2 = dx * dx + dy * dy + dz * dz;
            if best.map(|(bd, _)| d2 < bd).unwrap_or(true) {
                best = Some((d2, *h));
            }
        }
        best.map(|(_, h)| h.clamp(self.h_min, self.h_max))
            .unwrap_or(self.h_max)
    }
}

fn triangle_size(mesh: &Mesh, element: &Element) -> f64 {
    let p = [
        mesh.nodes.get(&element.node_ids[0]).map(|n| n.position),
        mesh.nodes.get(&element.node_ids[1]).map(|n| n.position),
        mesh.nodes.get(&element.node_ids[2]).map(|n| n.position),
    ];
    if let [Some(a), Some(b), Some(c)] = p {
        let ux = b.x - a.x;
        let uy = b.y - a.y;
        let vx = c.x - a.x;
        let vy = c.y - a.y;
        (0.5 * (ux * vy - uy * vx).abs()).sqrt() * 2.0f64.sqrt()
    } else {
        1.0
    }
}

// ─── The driver loop ──────────────────────────────────────────────────────────

/// Loop configuration. Defaults: `max_steps = 5`, `error_tolerance = 0.05`,
/// `h_min = 0.01`, `h_max = 1.0`.
#[derive(Debug, Clone)]
pub struct AdaptiveConfig {
    pub max_steps: usize,
    pub error_tolerance: f64,
    pub h_min: f64,
    pub h_max: f64,
}

impl Default for AdaptiveConfig {
    fn default() -> Self {
        AdaptiveConfig {
            max_steps: 5,
            error_tolerance: 0.05,
            h_min: 0.01,
            h_max: 1.0,
        }
    }
}

/// One loop step's bookkeeping.
#[derive(Debug, Clone)]
pub struct AdaptiveStep {
    pub step: usize,
    pub nodes: usize,
    pub elements: usize,
    pub max_indicator: f64,
    /// The per-element indicators of this step's mesh (element order =
    /// mesh element order), retained so consumers can pin any statistic
    /// (the global max carries boundary-anchored noise — the P3 lesson).
    pub indicators: Vec<f64>,
}

/// The adaptive loop report: per-step bookkeeping plus the per-step
/// meshes (renumbered per step) and the final mesh.
#[derive(Debug, Clone)]
pub struct AdaptiveReport {
    pub steps: Vec<AdaptiveStep>,
    pub meshes: Vec<Mesh>,
    pub mesh: Mesh,
    pub converged: bool,
}

/// The 2-D driver: mesh → estimate → update → re-mesh.
pub struct AdaptiveLoop2D<'a> {
    domain: &'a crate::traits::Domain2D,
    config: AdaptiveConfig,
}

impl<'a> AdaptiveLoop2D<'a> {
    pub fn new(domain: &'a crate::traits::Domain2D, config: &AdaptiveConfig) -> Self {
        AdaptiveLoop2D {
            domain,
            config: config.clone(),
        }
    }

    /// Run the loop. `u_h` samples the solution at a node position —
    /// supply the analytic reference, an interpolant of the computed
    /// field, whatever nodal solution drives the adaptivity.
    pub fn run(
        &self,
        u_h: impl Fn([f64; 3]) -> f64,
    ) -> Result<AdaptiveReport, AdaptiveError> {
        if !(self.config.error_tolerance.is_finite() && self.config.error_tolerance > 0.0) {
            return Err(AdaptiveError::InvalidTolerance {
                value: self.config.error_tolerance,
            });
        }
        let estimator = GradientRecovery;
        let mut steps = Vec::new();
        let mut meshes = Vec::new();
        let mut size_field: Option<AdaptiveSizeField> = None;
        let mut current: Option<Mesh> = None;

        for step in 0..=self.config.max_steps {
            // the mesher clamps every size-field answer into
            // [min_size, max_size] — the params MUST carry the loop's
            // clamp range or the h_min refinement never reaches the mesher
            let params = MeshParams {
                element_size: self.config.h_max,
                min_size: self.config.h_min,
                max_size: self.config.h_max,
                optimize_passes: 3,
            };
            let meshed = match &size_field {
                Some(field) => FrontalDelaunay2D::new()
                    .with_size_field(field)
                    .mesh_2d(self.domain, &params)
                    .map_err(|e| AdaptiveError::MeshingFailed { step, detail: e.to_string() })?,
                None => FrontalDelaunay2D::new()
                    .mesh_2d(self.domain, &params)
                    .map_err(|e| AdaptiveError::MeshingFailed { step, detail: e.to_string() })?,
            };
            current = Some(meshed);
            meshes.push(current.clone().expect("mesh"));
            // nodal solution in the dense sorted-node order
            let mesh_now = current.as_ref().expect("mesh");
            let mut sorted: Vec<&rmsh_model::Node> = mesh_now.nodes.values().collect();
            sorted.sort_by_key(|n| n.id);
            let u: Vec<f64> = sorted
                .iter()
                .map(|n| u_h([n.position.x, n.position.y, n.position.z]))
                .collect();

            let indicators = estimator.indicators(mesh_now, &u)?;
            let max_indicator = indicators.iter().copied().fold(0.0, f64::max);
            steps.push(AdaptiveStep {
                step,
                nodes: mesh_now.node_count(),
                elements: mesh_now.elements.len(),
                max_indicator,
                indicators: indicators.clone(),
            });

            let converged = max_indicator <= self.config.error_tolerance;
            if converged || step == self.config.max_steps {
                break;
            }

            let field = AdaptiveSizeField::build(
                current.as_ref().expect("mesh"),
                &indicators,
                self.config.error_tolerance,
                self.config.h_min,
                self.config.h_max,
            )?;
            size_field = Some(field);
        }

        let converged = steps
            .last()
            .map(|s| s.max_indicator <= self.config.error_tolerance)
            .unwrap_or(false);
        Ok(AdaptiveReport {
            steps,
            meshes,
            mesh: current.unwrap_or_default(),
            converged,
        })
    }
}
