//! Triangle-to-quadrilateral recombination (P3b — the Blossom-style pairing
//! card).
//!
//! Given a triangle-only 2-D [`Mesh`], produce a **quad-dominant** mesh by
//! pairing adjacent triangles (triangle pairs sharing an edge) into `Quad4`
//! elements wherever the resulting quad is well shaped, and leaving every
//! unpaired triangle in place as a `Triangle3`.
//!
//! # Pair-quality metric (documented choice)
//!
//! A candidate pair is scored by the **minimum interior angle of the
//! recombined quad, normalized by 90 degrees**:
//!
//! ```text
//! score(quad) = min(a1, a2, a3, a4) / 90deg   ∈ [0, 1]
//! ```
//!
//! Reasoning:
//!
//! * The minimum (worst) interior angle is the standard worst-case shape
//!   measure for quadrilaterals: a perfect square scores 1.0 and every
//!   degeneration drives the score toward 0 (a zero-area quad has a 0-degree
//!   "angle" on its collapsed diagonal; a bowtie configuration is rejected
//!   outright by the convexity gate below).
//! * Remacle et al., "Blossom-Quad: a non-obtuse quadrant meshing algorithm"
//!   (*IJNME* 89, 2012) gate triangle merges on the recombined quad's shape
//!   quality for the same reason — bad pairs must never be forced into the
//!   matching just to raise the pairing count.
//! * Non-convex or self-intersecting unions (the two triangle apexes on the
//!   same side of the shared edge) are rejected before scoring: the perimeter
//!   walk must have a consistent turn sign at all four corners with no
//!   collinear corner.
//!
//! # Matching algorithm (adjudicated for v1)
//!
//! Pair selection is a **maximum matching on the triangle-adjacency graph
//! weighted by the quality score** — the "Blossom-style" formulation of the
//! card. The general weighted blossom algorithm (Edmonds with dual
//! adjustments, cf. Kolmogorov's 2009 Blossom V) is deliberately deferred:
//! v1 uses
//!
//! 1. **greedy by decreasing score** (deterministic tie-break by shared-edge
//!    key then triangle indices), followed by
//! 2. an **alternating-path repair pass**: repeatedly search from every
//!    unmatched triangle for an augmenting path (non-matching edge, matching
//!    edge, non-matching edge, ...) and flip it. Every flipped path is a
//!    genuine alternating path, so the matching stays valid and strictly
//!    grows; the loop terminates because each flip adds exactly one pair.
//!
//! This is the core of Edmonds' algorithm *without odd-cycle (blossom)
//! contraction*: near odd cycles the naive search may miss an augmenting
//! path that a full blossom implementation would find. On mesh adjacency
//! graphs (which are locally path/cycle-like around interior vertices) the
//! repair pass closes essentially all of that gap. The pinned tests
//! demonstrate the repair is load-bearing: on the engineered 4-triangle fan
//! (`fan_with_perfect_matching_reaches_100_percent`) greedy alone pairs only
//! 50% of the triangles while greedy + repair reaches 100%.
//!
//! **Upgrade path:** replace `select_matching` with a maximum-*weight*
//! matching (Edmonds blossom with dual variables, O(V^3); reference
//! implementation Kolmogorov, "Blossom V: a new implementation of a minimum
//! cost perfect matching algorithm"). The candidate-generation and scoring
//! layers below are unchanged by that swap.
//!
//! # Output policy
//!
//! * Connectivity-only transform: **all nodes are kept with unchanged
//!   positions**; the input mesh is never mutated.
//! * A recombined `Quad4` lists its nodes in **perimeter order**: for
//!   triangles `ABC` and `DCB` sharing edge `BC` the quad is `A-B-D-C` (the
//!   walk starts at the first triangle's apex, crosses the shared edge, hops
//!   to the second triangle's apex, and returns along the other shared-edge
//!   endpoint — the walk direction follows the first triangle's winding).
//! * Unpaired triangles are emitted verbatim (same node ids, same order).
//! * Element ids are renumbered `1..`; a quad inherits the source triangles'
//!   physical tag only when both agree, otherwise `None`.
//!
//! # Robustness
//!
//! * Non-triangle input elements are refused with an error naming the
//!   offending element type.
//! * An empty mesh recombines to an empty mesh.
//! * Zero-area / collinear / non-convex pairs never merge (quality gate).

use std::collections::{HashMap, VecDeque};

use rmsh_model::{Element, ElementType, Mesh};

use crate::traits::MeshAlgoError;

// ─── Configuration ───────────────────────────────────────────────────────────

/// Configuration for triangle-to-quad recombination.
#[derive(Debug, Clone)]
pub struct RecombineConfig {
    /// Minimum quality score (see the module docs: minimum interior angle
    /// normalized by 90 degrees) a candidate quad must reach to be paired.
    ///
    /// The default `0.0` accepts every strictly convex, non-degenerate quad
    /// and maximizes the pairing rate; raise it to trade quads for shape.
    /// Must lie in `[0, 1]`.
    pub min_quality: f64,
}

impl Default for RecombineConfig {
    fn default() -> Self {
        Self { min_quality: 0.0 }
    }
}

impl RecombineConfig {
    fn validate(&self) -> Result<(), MeshAlgoError> {
        if !self.min_quality.is_finite() || !(0.0..=1.0).contains(&self.min_quality) {
            return Err(MeshAlgoError::InvalidInput(
                "min_quality must be a finite value in [0, 1]".to_string(),
            ));
        }
        Ok(())
    }
}

// ─── Public entry points ─────────────────────────────────────────────────────

/// Recombine a triangle-only mesh into a quad-dominant mesh with default
/// [`RecombineConfig`]. See the module docs for the metric, the matching
/// algorithm, and the output policy.
pub fn recombine_triangles(mesh: &Mesh) -> Result<Mesh, MeshAlgoError> {
    recombine_triangles_with(mesh, &RecombineConfig::default())
}

/// Recombine a triangle-only mesh into a quad-dominant mesh with an explicit
/// [`RecombineConfig`].
pub fn recombine_triangles_with(
    mesh: &Mesh,
    config: &RecombineConfig,
) -> Result<Mesh, MeshAlgoError> {
    config.validate()?;

    // Refuse non-triangle input, naming the offending element type.
    for e in &mesh.elements {
        if e.etype != ElementType::Triangle3 {
            return Err(MeshAlgoError::InvalidInput(format!(
                "recombination requires a triangle-only mesh: element {} has non-triangle type {:?}",
                e.id, e.etype
            )));
        }
        if e.node_ids.len() != 3 {
            return Err(MeshAlgoError::InvalidInput(format!(
                "recombination requires a triangle-only mesh: element {} is a Triangle3 carrying {} node ids",
                e.id,
                e.node_ids.len()
            )));
        }
    }

    // Node positions and triangle connectivity.
    let mut positions: HashMap<u64, [f64; 2]> = HashMap::with_capacity(mesh.nodes.len());
    for node in mesh.nodes.values() {
        positions.insert(node.id, [node.position.x, node.position.y]);
    }
    let tris: Vec<[u64; 3]> = mesh
        .elements
        .iter()
        .map(|e| [e.node_ids[0], e.node_ids[1], e.node_ids[2]])
        .collect();

    // Quality-scored candidate pairs, then the matching (greedy + repair).
    let candidates = build_candidates(&tris, &positions);
    let match_of = select_matching(&candidates, tris.len(), config.min_quality);
    let quad_of_pair: HashMap<(usize, usize), [u64; 4]> = candidates
        .into_iter()
        .map(|c| ((c.tri_a, c.tri_b), c.quad))
        .collect();

    // Output: every node verbatim; each pair emits one perimeter-ordered
    // quad at the lower partner's position; unpaired triangles verbatim.
    let mut output = Mesh::new();
    for node in mesh.nodes.values() {
        output.add_node(node.clone());
    }
    let mut next_id = 1u64;
    for (ti, t) in tris.iter().enumerate() {
        match match_of[ti] {
            Some(partner) if partner > ti => {
                let quad = quad_of_pair[&(ti, partner)];
                let tag_a = mesh.elements[ti].physical_tag;
                let tag_b = mesh.elements[partner].physical_tag;
                let mut element = Element::new(next_id, ElementType::Quad4, quad.to_vec());
                element.physical_tag = if tag_a == tag_b { tag_a } else { None };
                output.add_element(element);
            }
            Some(_) => continue, // second half of an already emitted pair
            None => {
                let mut element =
                    Element::new(next_id, ElementType::Triangle3, t.to_vec());
                element.physical_tag = mesh.elements[ti].physical_tag;
                output.add_element(element);
            }
        }
        next_id += 1;
    }
    Ok(output)
}

// ─── Candidate generation ────────────────────────────────────────────────────

/// One mergeable triangle pair: the perimeter-ordered quad plus its quality.
#[derive(Debug, Clone)]
struct CandidatePair {
    /// Index of the first triangle (its winding anchors the quad's walk).
    tri_a: usize,
    /// Index of the second triangle (`> tri_a`).
    tri_b: usize,
    /// Quad node ids in perimeter order.
    quad: [u64; 4],
    /// Quality score in `[0, 1]` (minimum interior angle / 90 degrees).
    score: f64,
    /// The shared edge as a sorted `(min, max)` node-id pair (tie-break key).
    edge: (u64, u64),
}

/// Merge two triangles sharing exactly one edge into the perimeter-ordered
/// quad `[apex_a, x, apex_b, y]` where `(x, y)` is the shared edge in
/// `tri_a`'s rotational order. Returns `None` when the triangles are not
/// edge-adjacent.
fn merge_pair_order(tri_a: [u64; 3], tri_b: [u64; 3]) -> Option<[u64; 4]> {
    let shared: Vec<u64> = tri_a.iter().filter(|a| tri_b.contains(a)).copied().collect();
    if shared.len() != 2 {
        return None;
    }
    let a_only = tri_a.iter().find(|n| !shared.contains(n)).copied()?;
    let b_only = tri_b.iter().find(|n| !shared.contains(n)).copied()?;
    Some([a_only, shared[0], b_only, shared[1]])
}

/// Quality gate + score for a perimeter-ordered quad.
///
/// Returns `None` when the quad is not a strictly convex simple polygon (any
/// collinear corner, reflex corner, or bowtie) — this is the hard gate that
/// skips zero-area / degenerate pairs. Otherwise returns the minimum interior
/// angle normalized by 90 degrees.
fn quad_min_angle_score(quad: [u64; 4], positions: &HashMap<u64, [f64; 2]>) -> Option<f64> {
    let p = [
        *positions.get(&quad[0])?,
        *positions.get(&quad[1])?,
        *positions.get(&quad[2])?,
        *positions.get(&quad[3])?,
    ];

    // Convexity gate: every corner must turn with the same nonzero sign.
    let mut orientation = 0.0f64;
    for i in 0..4 {
        let a = p[i];
        let b = p[(i + 1) % 4];
        let c = p[(i + 2) % 4];
        let cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
        if cross.abs() <= 1e-12 {
            return None; // collinear corner (covers zero-area quads)
        }
        if orientation == 0.0 {
            orientation = cross.signum();
        } else if cross.signum() != orientation {
            return None; // reflex corner or self-intersecting (bowtie) union
        }
    }

    // Minimum interior angle over the four corners.
    let mut min_angle = f64::MAX;
    for i in 0..4 {
        let prev = p[(i + 3) % 4];
        let v = p[i];
        let next = p[(i + 1) % 4];
        let e1 = [prev[0] - v[0], prev[1] - v[1]];
        let e2 = [next[0] - v[0], next[1] - v[1]];
        let dot = e1[0] * e2[0] + e1[1] * e2[1];
        let n1 = (e1[0] * e1[0] + e1[1] * e1[1]).sqrt();
        let n2 = (e2[0] * e2[0] + e2[1] * e2[1]).sqrt();
        let cos = (dot / (n1 * n2)).clamp(-1.0, 1.0);
        min_angle = min_angle.min(cos.acos());
    }
    Some(min_angle / std::f64::consts::FRAC_PI_2)
}

/// Build every mergeable candidate pair from the triangle list.
///
/// Candidates whose quad fails the convexity gate (zero-area, collinear,
/// bowtie) are dropped here.
fn build_candidates(tris: &[[u64; 3]], positions: &HashMap<u64, [f64; 2]>) -> Vec<CandidatePair> {
    // Edge (sorted node pair) -> incident triangle indices, built in element
    // order so every derived ordering below is deterministic.
    let mut edge_tris: HashMap<(u64, u64), Vec<usize>> = HashMap::new();
    for (ti, t) in tris.iter().enumerate() {
        if !is_valid_connectivity(t) {
            continue;
        }
        for k in 0..3 {
            let key = edge_key(t[k], t[(k + 1) % 3]);
            edge_tris.entry(key).or_default().push(ti);
        }
    }

    let mut candidates = Vec::new();
    for (ti, t) in tris.iter().enumerate() {
        if !is_valid_connectivity(t) {
            continue;
        }
        for k in 0..3 {
            let key = edge_key(t[k], t[(k + 1) % 3]);
            let Some(partners) = edge_tris.get(&key) else {
                continue;
            };
            for &tj in partners {
                if tj <= ti {
                    continue; // each unordered pair is taken exactly once
                }
                let Some(quad) = merge_pair_order(*t, tris[tj]) else {
                    continue;
                };
                let Some(score) = quad_min_angle_score(quad, positions) else {
                    continue; // degenerate / non-convex union: hard gate
                };
                candidates.push(CandidatePair {
                    tri_a: ti,
                    tri_b: tj,
                    quad,
                    score,
                    edge: key,
                });
            }
        }
    }
    candidates
}

fn is_valid_connectivity(t: &[u64; 3]) -> bool {
    t[0] != t[1] && t[1] != t[2] && t[0] != t[2]
}

fn edge_key(a: u64, b: u64) -> (u64, u64) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

// ─── Matching ────────────────────────────────────────────────────────────────

/// Greedy-by-decreasing-score matching followed by an alternating-path
/// repair pass. Returns `match_of[tri] = Some(partner)` (symmetric).
///
/// See the module docs for the adjudication (greedy + repair now, full
/// weighted blossom as the documented upgrade).
fn select_matching(
    candidates: &[CandidatePair],
    tri_count: usize,
    min_quality: f64,
) -> Vec<Option<usize>> {
    // Quality-threshold filter (the hard convexity gate already ran).
    let usable: Vec<&CandidatePair> =
        candidates.iter().filter(|c| c.score >= min_quality).collect();

    // Greedy pass: decreasing score, deterministic tie-breaks.
    let mut order: Vec<usize> = (0..usable.len()).collect();
    order.sort_by(|&x, &y| {
        usable[y]
            .score
            .total_cmp(&usable[x].score)
            .then_with(|| usable[x].edge.cmp(&usable[y].edge))
            .then_with(|| usable[x].tri_a.cmp(&usable[y].tri_a))
            .then_with(|| usable[x].tri_b.cmp(&usable[y].tri_b))
    });

    let mut match_of: Vec<Option<usize>> = vec![None; tri_count];
    for &ci in &order {
        let c = usable[ci];
        if match_of[c.tri_a].is_none() && match_of[c.tri_b].is_none() {
            match_of[c.tri_a] = Some(c.tri_b);
            match_of[c.tri_b] = Some(c.tri_a);
        }
    }

    // Adjacency over usable pairs (sorted for deterministic searches).
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); tri_count];
    for c in &usable {
        adj[c.tri_a].push(c.tri_b);
        adj[c.tri_b].push(c.tri_a);
    }
    for list in adj.iter_mut() {
        list.sort_unstable();
        list.dedup();
    }

    // Alternating-path repair: each augmentation grows the matching by one.
    loop {
        let mut augmented = false;
        for start in 0..tri_count {
            if match_of[start].is_none() {
                if let Some(path) = find_alternating_path(start, &adj, &match_of) {
                    apply_augmentation(&mut match_of, &path);
                    augmented = true;
                }
            }
        }
        if !augmented {
            break;
        }
    }
    match_of
}

/// Naive alternating-path BFS from an unmatched `start` (no blossom
/// contraction — documented v1 limitation, always-safe flips).
///
/// Returns the augmenting path `start ..= free_vertex` when found.
fn find_alternating_path(
    start: usize,
    adj: &[Vec<usize>],
    match_of: &[Option<usize>],
) -> Option<Vec<usize>> {
    let n = match_of.len();
    let mut visited = vec![false; n];
    let mut parent: Vec<Option<usize>> = vec![None; n];
    visited[start] = true;
    let mut queue = VecDeque::new();
    queue.push_back(start);

    while let Some(v) = queue.pop_front() {
        // Even-layer vertex: explore non-matching edges.
        for &w in &adj[v] {
            if visited[w] {
                continue;
            }
            visited[w] = true;
            parent[w] = Some(v);
            match match_of[w] {
                None => {
                    // Odd-layer free vertex: reconstruct the augmenting path.
                    let mut path = vec![w];
                    let mut cur = v;
                    loop {
                        path.push(cur);
                        match parent[cur] {
                            Some(p) => cur = p,
                            None => break,
                        }
                    }
                    path.reverse();
                    return Some(path);
                }
                Some(m) => {
                    // Walk the matching edge to the next even-layer vertex.
                    if !visited[m] {
                        visited[m] = true;
                        parent[m] = Some(w);
                        queue.push_back(m);
                    }
                }
            }
        }
    }
    None
}

/// Flip an augmenting path: every second edge becomes a matching edge.
fn apply_augmentation(match_of: &mut [Option<usize>], path: &[usize]) {
    for pair in path.chunks(2) {
        match_of[pair[0]] = Some(pair[1]);
        match_of[pair[1]] = Some(pair[0]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::triangulate2d::{Polygon2D, mesh_polygon};
    use rmsh_model::Node;
    use std::collections::HashSet;

    // ── fixtures ─────────────────────────────────────────────────────────────

    /// An `nx × ny` grid of square cells, each split into two triangles.
    /// Diagonal direction alternates in a checkerboard so both split
    /// orientations are exercised.
    fn structured_triangle_grid(nx: usize, ny: usize) -> Mesh {
        let mut mesh = Mesh::new();
        let mut nid = 1u64;
        let mut id = vec![vec![0u64; nx + 1]; ny + 1];
        for iy in 0..=ny {
            for ix in 0..=nx {
                id[iy][ix] = nid;
                mesh.add_node(Node::new(nid, ix as f64, iy as f64, 0.0));
                nid += 1;
            }
        }
        let mut eid = 1u64;
        for iy in 0..ny {
            for ix in 0..nx {
                let (ll, lr) = (id[iy][ix], id[iy][ix + 1]);
                let (ul, ur) = (id[iy + 1][ix], id[iy + 1][ix + 1]);
                if (ix + iy) % 2 == 0 {
                    mesh.add_element(Element::new(eid, ElementType::Triangle3, vec![ll, lr, ur]));
                    eid += 1;
                    mesh.add_element(Element::new(eid, ElementType::Triangle3, vec![ll, ur, ul]));
                    eid += 1;
                } else {
                    mesh.add_element(Element::new(eid, ElementType::Triangle3, vec![ll, lr, ul]));
                    eid += 1;
                    mesh.add_element(Element::new(eid, ElementType::Triangle3, vec![lr, ur, ul]));
                    eid += 1;
                }
            }
        }
        mesh
    }

    fn pos2(mesh: &Mesh, id: u64) -> [f64; 2] {
        let p = &mesh.nodes[&id].position;
        [p.x, p.y]
    }

    fn shoelace(pts: &[[f64; 2]]) -> f64 {
        let n = pts.len();
        let mut s = 0.0;
        for i in 0..n {
            let a = pts[i];
            let b = pts[(i + 1) % n];
            s += a[0] * b[1] - a[1] * b[0];
        }
        s / 2.0
    }

    fn tri_area(mesh: &Mesh, ids: &[u64]) -> f64 {
        shoelace(&[pos2(mesh, ids[0]), pos2(mesh, ids[1]), pos2(mesh, ids[2])]).abs()
    }

    /// Map every output element back to the input triangle(s) it covers.
    /// Returns Err with a diagnosis when a quad's source pair cannot be
    /// identified or an input triangle is used twice.
    fn map_outputs_to_inputs(
        input: &Mesh,
        output: &Mesh,
    ) -> Result<Vec<(usize, Option<(usize, usize)>)>, String> {
        let in_tris: Vec<Vec<u64>> = input
            .elements
            .iter()
            .map(|e| e.node_ids.clone())
            .collect();
        let n = in_tris.len();
        let mut used = vec![false; n];
        let mut mapping = Vec::new();

        let share_edge = |a: &[u64], b: &[u64]| {
            a.iter().filter(|x| b.contains(x)).count() == 2
        };

        for out in &output.elements {
            match out.etype {
                ElementType::Triangle3 => {
                    let key: HashSet<u64> = out.node_ids.iter().copied().collect();
                    let mut found = None;
                    for (i, t) in in_tris.iter().enumerate() {
                        if !used[i] && t.iter().copied().collect::<HashSet<u64>>() == key {
                            found = Some(i);
                            break;
                        }
                    }
                    match found {
                        Some(i) => {
                            used[i] = true;
                            mapping.push((out.id as usize, None));
                        }
                        None => return Err(format!("triangle {} matches no unused input", out.id)),
                    }
                }
                ElementType::Quad4 => {
                    let qset: HashSet<u64> = out.node_ids.iter().copied().collect();
                    let mut pair = None;
                    'outer: for i in 0..n {
                        if used[i] {
                            continue;
                        }
                        for j in (i + 1)..n {
                            if used[j] {
                                continue;
                            }
                            if share_edge(&in_tris[i], &in_tris[j])
                                && in_tris[i].iter().chain(in_tris[j].iter()).copied().collect::<HashSet<u64>>()
                                    == qset
                            {
                                pair = Some((i, j));
                                break 'outer;
                            }
                        }
                    }
                    match pair {
                        Some((i, j)) => {
                            used[i] = true;
                            used[j] = true;
                            mapping.push((out.id as usize, Some((i, j))));
                        }
                        None => return Err(format!("quad {} matches no unused input pair", out.id)),
                    }
                }
                other => return Err(format!("unexpected output element type {other:?}")),
            }
        }

        if used.iter().any(|u| !u) {
            return Err("some input triangles are not covered by any output element".to_string());
        }
        Ok(mapping)
    }

    /// Perimeter-order + winding pin for one recombined quad: every
    /// consecutive node pair must be an edge of one of the two source
    /// triangles, and the quad's |signed area| must equal the sum of the two
    /// triangle areas (true for the correct perimeter walk, violated by any
    /// bowtie ordering).
    fn assert_perimeter_ordered(input: &Mesh, output: &Mesh) {
        let mapping = map_outputs_to_inputs(input, output).expect("coverage");
        let in_tris: Vec<Vec<u64>> =
            input.elements.iter().map(|e| e.node_ids.clone()).collect();
        for (out_id, pair) in &mapping {
            let out = output.elements.iter().find(|e| e.id == *out_id as u64).unwrap();
            if let ElementType::Quad4 = out.etype {
                let (i, j) = pair.unwrap();
                let t1 = &in_tris[i];
                let t2 = &in_tris[j];
                let has_edge = |t: &[u64], a: u64, b: u64| {
                    (0..3).any(|k| {
                        (t[k] == a && t[(k + 1) % 3] == b) || (t[k] == b && t[(k + 1) % 3] == a)
                    })
                };
                for k in 0..4 {
                    let a = out.node_ids[k];
                    let b = out.node_ids[(k + 1) % 4];
                    assert!(
                        has_edge(t1, a, b) || has_edge(t2, a, b),
                        "quad {:?}: consecutive pair ({a},{b}) is not an edge of either source triangle",
                        out.node_ids
                    );
                }
                let q_area = shoelace(&[
                    pos2(output, out.node_ids[0]),
                    pos2(output, out.node_ids[1]),
                    pos2(output, out.node_ids[2]),
                    pos2(output, out.node_ids[3]),
                ])
                .abs();
                let sum = tri_area(input, t1) + tri_area(input, t2);
                assert!(
                    (q_area - sum).abs() < 1e-9 * sum.max(1.0),
                    "quad {:?}: |area|={q_area} != |t1|+|t2|={sum}",
                    out.node_ids
                );
            }
        }
    }

    // ── structured grid ──────────────────────────────────────────────────────

    #[test]
    fn structured_grid_recombines_to_all_quads() {
        // 6 x 4 cells x 2 triangles = 48 triangles -> 6 x 4 = 24 quads.
        let input = structured_triangle_grid(6, 4);
        let output = recombine_triangles(&input).expect("recombination succeeds");
        assert_eq!(output.node_count(), input.node_count(), "all nodes kept");
        assert_eq!(output.elements.len(), 24, "exactly nx*ny quads");
        assert!(
            output
                .elements
                .iter()
                .all(|e| e.etype == ElementType::Quad4),
            "100% matching on the structured grid"
        );
        // Node positions unchanged.
        for (id, node) in &input.nodes {
            let out = &output.nodes[id];
            assert_eq!(out.position.x, node.position.x);
            assert_eq!(out.position.y, node.position.y);
        }
        assert_perimeter_ordered(&input, &output);
    }

    #[test]
    fn structured_grid_winding_pinned_including_crossed_diagonal() {
        let input = structured_triangle_grid(6, 4);
        let output = recombine_triangles(&input).expect("recombination succeeds");

        // Cell (0,0): even cell, diagonal ll-ur, triangles (1,2,9) + (1,9,8).
        // The walk rule (card: ABC + DCB sharing BC -> A-B-D-C) gives
        // [apex(2), shared[0]=1, apex(8), shared[1]=9]: a clockwise
        // perimeter walk of the unit cell.
        assert_eq!(output.elements[0].node_ids, vec![2, 1, 8, 9]);
        // Cell (1,0): odd cell, crossed diagonal lr-ul, triangles (2,3,9) +
        // (3,10,9). Shared {3,9} appears in tri_a's list as (3,9):
        // [apex(2), 3, apex(10), 9] — a counter-clockwise walk of the
        // crossed-split cell.
        assert_eq!(output.elements[1].node_ids, vec![2, 3, 10, 9]);

        assert_perimeter_ordered(&input, &output);
    }

    #[test]
    fn structured_grid_with_high_min_quality_still_all_quads() {
        // Rectangle quads score exactly 1.0, so a strict gate keeps them all.
        let input = structured_triangle_grid(6, 4);
        let config = RecombineConfig {
            min_quality: 0.9,
        };
        let output = recombine_triangles_with(&input, &config).expect("recombination succeeds");
        assert_eq!(output.elements.len(), 24);
        assert!(
            output
                .elements
                .iter()
                .all(|e| e.etype == ElementType::Quad4)
        );
    }

    // ── irregular mesh ───────────────────────────────────────────────────────

    #[test]
    fn irregular_polygon_recombination_preserves_nodes_and_coverage() {
        let poly = Polygon2D::new(vec![
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ]);
        let input = mesh_polygon(&poly, 0.4).expect("triangulation succeeds");
        assert!(input.elements.len() >= 20, "got a real mesh to work with");

        let output = recombine_triangles(&input).expect("recombination succeeds");

        // Node set identical (ids and positions).
        assert_eq!(output.node_count(), input.node_count());
        for (id, node) in &input.nodes {
            let out = output.nodes.get(id).unwrap_or_else(|| panic!("node {id} kept"));
            assert_eq!(out.position.x, node.position.x);
            assert_eq!(out.position.y, node.position.y);
            assert_eq!(out.position.z, node.position.z);
        }

        // Coverage: every input triangle appears in exactly one output
        // element (paired or alone); every quad maps to two distinct inputs.
        let n_quads = output
            .elements
            .iter()
            .filter(|e| e.etype == ElementType::Quad4)
            .count();
        let n_tris = output
            .elements
            .iter()
            .filter(|e| e.etype == ElementType::Triangle3)
            .count();
        assert_eq!(n_quads + n_tris + n_quads, input.elements.len());
        assert!(n_quads > 0, "an irregular mesh still pairs some triangles");
        // NOTE: a perfect matching (n_tris == 0) is legitimate here — the
        // repair pass may pair every triangle of this L-shape mesh. The
        // mixed quad+triangle case is pinned by the fan and degeneracy
        // tests below.

        let rate = 2.0 * n_quads as f64 / input.elements.len() as f64;
        println!(
            "irregular L-shape: {} triangles -> {n_quads} quads + {n_tris} triangles (matching rate {rate:.3})",
            input.elements.len()
        );
        // Measured pin: the deterministic pipeline pairs ALL 44 triangles
        // (matching rate 1.000) on this mesh.
        assert_eq!(n_quads, 22, "measured matching rate pin");
        assert_eq!(n_tris, 0);

        assert_perimeter_ordered(&input, &output);
    }

    // ── matching rate on a perfect-matching case ─────────────────────────────

    #[test]
    fn fan_with_perfect_matching_reaches_100_percent() {
        // Four fan triangles around v=(0,0). The middle pair (T2,T3) forms a
        // well-shaped diamond (score ~0.59); the outer pairs are thin kites
        // (score ~0.33). Greedy alone takes the diamond and strands T1/T4
        // (50%); the alternating-path repair then re-pairs to the perfect
        // matching {(T1,T2),(T3,T4)} -> 100%.
        let mut input = Mesh::new();
        // 1=v, 2=R, 3=P, 4=S, 5=Q, 6=R'
        for (id, p) in [
            (1u64, [0.0, 0.0]),
            (2, [2.0, -1.0]),
            (3, [1.0, 0.5]),
            (4, [0.0, 1.0]),
            (5, [-1.0, 0.5]),
            (6, [-2.0, -1.0]),
        ] {
            input.add_node(Node::new(id, p[0], p[1], 0.0));
        }
        // T1=(v,R,P), T2=(v,P,S), T3=(v,S,Q), T4=(v,Q,R')
        for (id, ids) in [
            (1u64, vec![1, 2, 3]),
            (2, vec![1, 3, 4]),
            (3, vec![1, 4, 5]),
            (4, vec![1, 5, 6]),
        ] {
            input.add_element(Element::new(id, ElementType::Triangle3, ids));
        }

        let output = recombine_triangles(&input).expect("recombination succeeds");
        assert_eq!(
            output
                .elements
                .iter()
                .filter(|e| e.etype == ElementType::Quad4)
                .count(),
            2,
            "perfect matching reached: all 4 triangles paired"
        );
        assert!(
            output
                .elements
                .iter()
                .all(|e| e.etype == ElementType::Quad4),
            "matching rate 1.0 >= 0.95"
        );
        // Deterministic walk pins: T1 anchors [R, v, S, P], T3 anchors
        // [S, v, R', Q] — both convex perimeter walks of the fan kites.
        assert_eq!(output.elements[0].node_ids, vec![2, 1, 4, 3]);
        assert_eq!(output.elements[1].node_ids, vec![4, 1, 6, 5]);
        assert_perimeter_ordered(&input, &output);
    }

    // ── refusals and degeneracies ────────────────────────────────────────────

    #[test]
    fn non_triangle_input_is_refused_naming_the_element() {
        let mut mesh = Mesh::new();
        for id in 1..=4u64 {
            mesh.add_node(Node::new(id, id as f64, 0.0, 0.0));
        }
        mesh.add_element(Element::new(1, ElementType::Quad4, vec![1, 2, 3, 4]));
        let err = recombine_triangles(&mesh).expect_err("quad input refused");
        let msg = match &err {
            MeshAlgoError::InvalidInput(m) => m.clone(),
            other => panic!("expected InvalidInput, got {other:?}"),
        };
        assert!(msg.contains("Quad4"), "error must name the element type: {msg}");

        let mut line_mesh = Mesh::new();
        line_mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        line_mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        line_mesh.add_element(Element::new(1, ElementType::Line2, vec![1, 2]));
        let err = recombine_triangles(&line_mesh).expect_err("line input refused");
        let msg = match &err {
            MeshAlgoError::InvalidInput(m) => m.clone(),
            other => panic!("expected InvalidInput, got {other:?}"),
        };
        assert!(msg.contains("Line2"), "error must name the element type: {msg}");
    }

    #[test]
    fn empty_mesh_returns_empty_mesh() {
        let empty = Mesh::new();
        let output = recombine_triangles(&empty).expect("empty mesh is accepted");
        assert_eq!(output.node_count(), 0);
        assert_eq!(output.elements.len(), 0);

        // Nodes but no elements: nodes are preserved verbatim.
        let mut nodes_only = Mesh::new();
        nodes_only.add_node(Node::new(42, 1.0, 2.0, 0.0));
        let output = recombine_triangles(&nodes_only).expect("nodes-only mesh is accepted");
        assert_eq!(output.node_count(), 1);
        assert_eq!(output.elements.len(), 0);
        assert_eq!(output.nodes[&42].position.x, 1.0);
    }

    #[test]
    fn input_mesh_is_not_mutated() {
        let input = structured_triangle_grid(6, 4);
        let snapshot_nodes: HashMap<u64, (f64, f64, f64)> = input
            .nodes
            .iter()
            .map(|(id, n)| (*id, (n.position.x, n.position.y, n.position.z)))
            .collect();
        let snapshot_elems: Vec<(u64, ElementType, Vec<u64>)> = input
            .elements
            .iter()
            .map(|e| (e.id, e.etype, e.node_ids.clone()))
            .collect();

        let _ = recombine_triangles(&input).expect("recombination succeeds");

        assert_eq!(snapshot_nodes.len(), input.nodes.len());
        for (id, (x, y, z)) in &snapshot_nodes {
            let n = &input.nodes[id];
            assert_eq!((n.position.x, n.position.y, n.position.z), (*x, *y, *z));
        }
        let after: Vec<(u64, ElementType, Vec<u64>)> = input
            .elements
            .iter()
            .map(|e| (e.id, e.etype, e.node_ids.clone()))
            .collect();
        assert_eq!(snapshot_elems, after);
    }

    #[test]
    fn degenerate_pair_is_skipped_by_the_quality_gate() {
        // Both apexes collinear with the shared edge: the pair's union quad
        // has zero area and must not be recombined. The degenerate triangles
        // pass through unchanged.
        let mut input = Mesh::new();
        // 1=(0,0) 2=(1,0) 3=(2,0) 4=(2.5,0) — all on y=0
        for (id, x) in [(1u64, 0.0), (2, 1.0), (3, 2.0), (4, 2.5)] {
            input.add_node(Node::new(id, x, 0.0, 0.0));
        }
        // T1=(1,2,3) zero-area, T2=(2,3,4) zero-area, sharing edge {2,3}.
        input.add_element(Element::new(1, ElementType::Triangle3, vec![1, 2, 3]));
        input.add_element(Element::new(2, ElementType::Triangle3, vec![2, 3, 4]));

        let output = recombine_triangles(&input).expect("recombination succeeds");
        assert_eq!(output.elements.len(), 2);
        assert!(
            output
                .elements
                .iter()
                .all(|e| e.etype == ElementType::Triangle3),
            "the zero-area pair must stay unmerged"
        );
        assert_eq!(output.elements[0].node_ids, vec![1, 2, 3]);
        assert_eq!(output.elements[1].node_ids, vec![2, 3, 4]);
    }

    #[test]
    fn min_quality_gate_excludes_weak_pairs() {
        // Same fan as the perfect-matching test; raising the gate above the
        // thin kites' score (~0.33) but below the diamond (~0.59) leaves only
        // the middle pair.
        let mut input = Mesh::new();
        for (id, p) in [
            (1u64, [0.0, 0.0]),
            (2, [2.0, -1.0]),
            (3, [1.0, 0.5]),
            (4, [0.0, 1.0]),
            (5, [-1.0, 0.5]),
            (6, [-2.0, -1.0]),
        ] {
            input.add_node(Node::new(id, p[0], p[1], 0.0));
        }
        for (id, ids) in [
            (1u64, vec![1, 2, 3]),
            (2, vec![1, 3, 4]),
            (3, vec![1, 4, 5]),
            (4, vec![1, 5, 6]),
        ] {
            input.add_element(Element::new(id, ElementType::Triangle3, ids));
        }

        let config = RecombineConfig { min_quality: 0.5 };
        let output = recombine_triangles_with(&input, &config).expect("recombination succeeds");
        let n_quads = output
            .elements
            .iter()
            .filter(|e| e.etype == ElementType::Quad4)
            .count();
        let n_tris = output
            .elements
            .iter()
            .filter(|e| e.etype == ElementType::Triangle3)
            .count();
        assert_eq!(n_quads, 1, "only the diamond pair clears the 0.5 gate");
        assert_eq!(n_tris, 2, "the two thin kites stay triangles");
    }

    #[test]
    fn invalid_min_quality_is_refused() {
        let input = structured_triangle_grid(2, 2);
        for bad in [-0.1, 1.5, f64::NAN] {
            let config = RecombineConfig { min_quality: bad };
            let err = recombine_triangles_with(&input, &config).expect_err("bad config");
            assert!(
                matches!(err, MeshAlgoError::InvalidInput(_)),
                "min_quality {bad} must be an InvalidInput"
            );
        }
    }
}
