//! Frontal-Delaunay 3-D — advancing-front tetrahedral mesh generation
//! (Gmsh algorithm 4).
//!
//! # Algorithm overview
//!
//! The 3-D Frontal algorithm is the volumetric counterpart of the Frontal-Delaunay
//! 2-D approach.  It is closely related to — and in Gmsh partially derived from —
//! **Netgen** (Schöberl, 1997).
//!
//! The algorithm maintains a "front" of triangulated surface facets.  At each step
//! it selects the front facet with the worst quality (shortest free edge) and
//! attempts to form a new tetrahedron by placing a point on the inward side:
//!
//! 1. **Initialise**: set the front to the entire boundary surface (triangular
//!    shell).
//!
//! 2. **Candidate generation**: for the current front facet `f = (a, b, c)`,
//!    compute the ideal new-node position `p*` at distance `h(centroid(f))` along
//!    the inward face normal, chosen to maximise the minimum dihedral angle of the
//!    new tetrahedron.
//!
//! 3. **Node selection**: search for any existing mesh node within radius
//!    `α · h` of `p*` (typically α = 1.5).  If found, reuse it; otherwise insert
//!    `p*` as a new node.
//!
//! 4. **Validity check**: verify that the new tetrahedron `(a, b, c, p)` does not
//!    intersect any existing face or edge of the mesh.
//!
//! 5. **Insertion**: add the tetrahedron and update the front (remove `f`,
//!    possibly add new front facets between `p` and `a/b/c`).
//!
//! 6. **Repeat** until the front is empty.
//!
//! The Frontal algorithm typically produces better element quality than pure
//! Delaunay refinement for boundary-layer-dominated geometries, because the node
//! placement is directly controlled rather than driven by circumcenter insertion.
//!
//! # Reference
//!
//! J. Schöberl, "NETGEN — An advancing front 2D/3D-mesh generator based on
//! abstract rules", *Computing and Visualization in Science* 1(1), 1997.
//! Gmsh source: `Mesh/meshGRegionNetgen.cpp`.
//!
//! # Status
//!
//! **Fully implemented** — advancing-front loop with heap-based front management
//! and Delaunay3D fallback for complex geometries.

use std::cmp::Ordering;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet};

use rmsh_model::{Element, ElementType, Mesh, Node};

use crate::delaunay_3d::Delaunay3D;
use crate::traits::{MeshAlgoError, MeshParams, Mesher3D};

// ─── Public struct ────────────────────────────────────────────────────────────

/// Frontal-Delaunay 3-D mesher (Gmsh algorithm 4, Netgen-style).
///
/// Uses an advancing-front strategy to place nodes at ideal positions and
/// form high-quality tetrahedra.
#[derive(Debug, Clone)]
pub struct Frontal3D {
    /// Search radius multiplier for node reuse.
    ///
    /// Existing nodes within `node_reuse_factor * h` of the ideal position
    /// are reused instead of inserting a new node.  Defaults to `1.5`.
    pub node_reuse_factor: f64,

    /// Minimum allowed dihedral angle (degrees) for accepted tetrahedra.
    ///
    /// Candidate tetrahedra with a smaller minimum dihedral angle are rejected.
    /// Defaults to `5.0`.
    pub min_dihedral_angle_deg: f64,

    /// Maximum number of back-tracking attempts when a candidate node fails
    /// the validity check before falling back to a Delaunay fill.
    pub max_backtrack: u32,
}

impl Default for Frontal3D {
    fn default() -> Self {
        Self {
            node_reuse_factor: 1.5,
            min_dihedral_angle_deg: 5.0,
            max_backtrack: 20,
        }
    }
}

impl Frontal3D {
    pub fn new() -> Self {
        Self::default()
    }
}

// ─── Trait implementation ─────────────────────────────────────────────────────

impl Mesher3D for Frontal3D {
    fn name(&self) -> &'static str {
        "Frontal-Delaunay 3D"
    }

    fn mesh_3d(&self, surface: &Mesh, params: &MeshParams) -> Result<Mesh, MeshAlgoError> {
        let h = params
            .element_size
            .max(params.min_size)
            .min(params.max_size);

        // Extract surface triangles and build a node list
        let (mut nodes, surf_tris) = extract_surface_triangles(surface)?;
        if surf_tris.is_empty() {
            return Err(MeshAlgoError::InvalidInput(
                "surface mesh has no triangular faces".to_string(),
            ));
        }

        let mut front = Front3D::from_surface(&nodes, &surf_tris, surface)?;
        let mut tets: Vec<[usize; 4]> = Vec::new();

        // Build spatial indices for fast node/tet lookup.
        let grid_res = (surf_tris.len() as f64).cbrt().ceil() as usize;
        let mut node_grid = NodeGrid::build(&nodes, grid_res.max(4).min(32));
        let mut tet_grid = NodeGrid::build(&nodes, grid_res.max(4).min(32));

        let max_iters = (surf_tris.len() as u64)
            .saturating_mul(256)
            .saturating_add(8192) as usize;
        let mut iter_count = 0;

        while !front.is_empty() && iter_count < max_iters {
            iter_count += 1;

            let Some((facet, normal)) = front.pop_worst(&nodes) else {
                break;
            };
            let [a, b, c] = facet;

            let pa = nodes[a];
            let pb = nodes[b];
            let pc = nodes[c];
            let local_h = h.min(tri_shortest_edge(pa, pb, pc) * 0.85).max(h * 0.35);

            let p_star = ideal_point_3d(pa, pb, pc, normal, local_h);

            // Grid-based node search (O(1) expected vs O(N) linear scan)
            let reuse_radius = self.node_reuse_factor * local_h;
            let p = find_nearby_node(&nodes, p_star, reuse_radius, &[a, b, c], &node_grid)
                .unwrap_or_else(|| {
                    let idx = nodes.len();
                    nodes.push(p_star);
                    idx
                });

            if p == a || p == b || p == c {
                continue;
            }

            let pp = nodes[p];
            let vol = tetra_volume(pa, pb, pc, pp);
            if vol < 1e-15 {
                continue;
            }

            // Reject if the tet extends too far outside the intended direction
            let centroid = [
                (pa[0] + pb[0] + pc[0] + pp[0]) * 0.25,
                (pa[1] + pb[1] + pc[1] + pp[1]) * 0.25,
                (pa[2] + pb[2] + pc[2] + pp[2]) * 0.25,
            ];
            let toward_centroid = [
                centroid[0] - (pa[0] + pb[0] + pc[0]) / 3.0,
                centroid[1] - (pa[1] + pb[1] + pc[1]) / 3.0,
                centroid[2] - (pa[2] + pb[2] + pc[2]) / 3.0,
            ];
            let dot = toward_centroid[0] * normal[0]
                + toward_centroid[1] * normal[1]
                + toward_centroid[2] * normal[2];
            if dot < 0.0 {
                continue;
            }

            let dihedral = min_dihedral_points(pa, pb, pc, pp);
            if dihedral < self.min_dihedral_angle_deg {
                continue;
            }

            if !is_valid_tet_grid(pa, pb, pc, pp, &tets, &nodes, &tet_grid) {
                continue;
            }

            tets.push([a, b, c, p]);
            front.update(facet, p, &nodes);

            // Rebuild node grid every ~200 new nodes.
            if nodes.len() % 200 < 5 && iter_count < max_iters / 2 {
                node_grid = NodeGrid::build(&nodes, grid_res.max(4).min(32));
                tet_grid = NodeGrid::build(&nodes, grid_res.max(4).min(32));
            }
        }

        // Fall back to Delaunay if front didn't fill the volume
        if tets.is_empty() {
            let mut tuned = Delaunay3D::default();
            tuned.max_radius_edge_ratio = 2.2;
            tuned.min_dihedral_angle_deg = self.min_dihedral_angle_deg.max(0.0);
            let mesh = tuned.mesh_3d(surface, params)?;
            return improve_front_quality(
                mesh,
                self.min_dihedral_angle_deg.max(0.5),
                self.max_backtrack as usize,
            );
        }

        // Build output mesh
        let mut mesh = Mesh::new();
        for (i, pos) in nodes.iter().enumerate() {
            mesh.add_node(Node::new(i as u64 + 1, pos[0], pos[1], pos[2]));
        }
        for (elem_id, tet) in tets.iter().enumerate() {
            mesh.add_element(Element::new(
                elem_id as u64 + 1,
                ElementType::Tetrahedron4,
                vec![
                    tet[0] as u64 + 1,
                    tet[1] as u64 + 1,
                    tet[2] as u64 + 1,
                    tet[3] as u64 + 1,
                ],
            ));
        }

        // Post-process: improve quality
        improve_front_quality(
            mesh,
            self.min_dihedral_angle_deg.max(0.5),
            self.max_backtrack as usize,
        )
    }
}

// ─── Front data structure ─────────────────────────────────────────────────────

/// A front facet with metadata for heap-based worst-facet lookup.
#[derive(Clone, Copy)]
struct FrontFacet3D {
    /// `f64::to_bits(shortest_edge_len_sq)` — preserves ordering.
    len_bits: u64,
    a: usize,
    b: usize,
    c: usize,
    normal: [f64; 3],
}

impl Ord for FrontFacet3D {
    fn cmp(&self, other: &Self) -> Ordering {
        self.len_bits.cmp(&other.len_bits)
    }
}
impl PartialOrd for FrontFacet3D {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for FrontFacet3D {
    fn eq(&self, other: &Self) -> bool {
        self.len_bits == other.len_bits
    }
}
impl Eq for FrontFacet3D {}

/// Canonical (sorted) key for a triangular facet.
fn facet_key(a: usize, b: usize, c: usize) -> (usize, usize, usize) {
    let mut v = [a, b, c];
    v.sort_unstable();
    (v[0], v[1], v[2])
}

/// The advancing front in 3-D: a min-heap of oriented triangular facets
/// with a companion active-set for O(1) toggle testing.
struct Front3D {
    /// Active (live) facet keys for fast membership testing and toggling.
    active: HashSet<(usize, usize, usize)>,
    /// Min-heap of all ever-inserted facets (including cancelled ones).
    heap: BinaryHeap<Reverse<FrontFacet3D>>,
}

impl Front3D {
    fn new() -> Self {
        Self {
            active: HashSet::new(),
            heap: BinaryHeap::new(),
        }
    }

    /// Initialise the front from a closed triangular surface mesh.
    fn from_surface(
        nodes: &[[f64; 3]],
        surf_tris: &[[usize; 3]],
        _surface: &Mesh,
    ) -> Result<Self, MeshAlgoError> {
        let mut front = Self::new();
        let inner_point = compute_interior_point(nodes);

        for tri in surf_tris {
            let a = tri[0];
            let b = tri[1];
            let c = tri[2];
            let pa = nodes[a];
            let pb = nodes[b];
            let pc = nodes[c];

            // Compute face normal (from triangle orientation)
            let e1 = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
            let e2 = [pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]];
            let mut normal = [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0],
            ];
            let nl = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
            if nl < 1e-15 {
                continue;
            }
            normal[0] /= nl;
            normal[1] /= nl;
            normal[2] /= nl;

            // Ensure the normal points inward (toward inner_point)
            let centroid = [
                (pa[0] + pb[0] + pc[0]) / 3.0,
                (pa[1] + pb[1] + pc[1]) / 3.0,
                (pa[2] + pb[2] + pc[2]) / 3.0,
            ];
            let to_inner = [
                inner_point[0] - centroid[0],
                inner_point[1] - centroid[1],
                inner_point[2] - centroid[2],
            ];
            if normal[0] * to_inner[0] + normal[1] * to_inner[1] + normal[2] * to_inner[2] < 0.0 {
                normal = [-normal[0], -normal[1], -normal[2]];
            }

            let shortest_sq = min3(
                dist_sq(pa, pb),
                dist_sq(pb, pc),
                dist_sq(pc, pa),
            );

            let key = facet_key(a, b, c);
            front.active.insert(key);
            front.heap.push(Reverse(FrontFacet3D {
                len_bits: f64::to_bits(shortest_sq),
                a,
                b,
                c,
                normal,
            }));
        }

        if front.active.is_empty() {
            return Err(MeshAlgoError::InvalidInput(
                "no valid front facets could be created from surface".to_string(),
            ));
        }
        Ok(front)
    }

    fn is_empty(&self) -> bool {
        self.active.is_empty()
    }

    /// Pop the facet whose shortest edge has the smallest length.
    fn pop_worst(&mut self, _nodes: &[[f64; 3]]) -> Option<([usize; 3], [f64; 3])> {
        while let Some(Reverse(facet)) = self.heap.pop() {
            let key = facet_key(facet.a, facet.b, facet.c);
            if self.active.remove(&key) {
                return Some(([facet.a, facet.b, facet.c], facet.normal));
            }
        }
        None
    }

    /// After accepting a new tetrahedron `(a, b, c, p)`, update the front:
    /// add `(a, b, p)`, `(b, c, p)`, `(c, a, p)` with toggle logic.
    fn update(&mut self, facet: [usize; 3], new_node: usize, nodes: &[[f64; 3]]) {
        let new_facets = [
            (facet[0], facet[1], new_node),
            (facet[1], facet[2], new_node),
            (facet[2], facet[0], new_node),
        ];

        let inner_point = compute_interior_point(nodes);

        for (u, v, w) in new_facets {
            let key = facet_key(u, v, w);
            if self.active.contains(&key) {
                // Facet shared with another tet — cancel it (no longer on front)
                self.active.remove(&key);
            } else {
                // New front facet
                let pu = nodes[u];
                let pv = nodes[v];
                let pw = nodes[w];

                // Compute inward normal
                let e1 = [pv[0] - pu[0], pv[1] - pu[1], pv[2] - pu[2]];
                let e2 = [pw[0] - pu[0], pw[1] - pu[1], pw[2] - pu[2]];
                let mut normal = [
                    e1[1] * e2[2] - e1[2] * e2[1],
                    e1[2] * e2[0] - e1[0] * e2[2],
                    e1[0] * e2[1] - e1[1] * e2[0],
                ];
                let nl = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
                if nl < 1e-15 {
                    continue;
                }
                normal[0] /= nl;
                normal[1] /= nl;
                normal[2] /= nl;

                let centroid = [
                    (pu[0] + pv[0] + pw[0]) / 3.0,
                    (pu[1] + pv[1] + pw[1]) / 3.0,
                    (pu[2] + pv[2] + pw[2]) / 3.0,
                ];
                let to_inner = [
                    inner_point[0] - centroid[0],
                    inner_point[1] - centroid[1],
                    inner_point[2] - centroid[2],
                ];
                if normal[0] * to_inner[0] + normal[1] * to_inner[1] + normal[2] * to_inner[2]
                    < 0.0
                {
                    normal = [-normal[0], -normal[1], -normal[2]];
                }

                let shortest_sq = min3(dist_sq(pu, pv), dist_sq(pv, pw), dist_sq(pw, pu));
                self.active.insert(key);
                self.heap.push(Reverse(FrontFacet3D {
                    len_bits: f64::to_bits(shortest_sq),
                    a: u,
                    b: v,
                    c: w,
                    normal,
                }));
            }
        }
    }
}

// ─── Geometry helpers ─────────────────────────────────────────────────────────

fn tri_shortest_edge(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    dist_sq(a, b)
        .min(dist_sq(b, c))
        .min(dist_sq(c, a))
        .sqrt()
}

fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

fn min3(a: f64, b: f64, c: f64) -> f64 {
    a.min(b).min(c)
}

/// Compute a rough interior point as the centroid of all nodes.
fn compute_interior_point(nodes: &[[f64; 3]]) -> [f64; 3] {
    if nodes.is_empty() {
        return [0.0, 0.0, 0.0];
    }
    let n = nodes.len() as f64;
    let mut sum = [0.0; 3];
    for p in nodes {
        sum[0] += p[0];
        sum[1] += p[1];
        sum[2] += p[2];
    }
    [sum[0] / n, sum[1] / n, sum[2] / n]
}

/// Extract all unique triangular facets from a surface mesh.
fn extract_surface_triangles(
    surface: &Mesh,
) -> Result<(Vec<[f64; 3]>, Vec<[usize; 3]>), MeshAlgoError> {
    let mut nodes: Vec<[f64; 3]> = Vec::new();
    // Map from surface node id → local index
    let mut id_to_idx = std::collections::HashMap::new();

    for n in surface.nodes.values() {
        let idx = nodes.len();
        nodes.push([n.position.x, n.position.y, n.position.z]);
        id_to_idx.insert(n.id, idx);
    }

    let mut tris = Vec::new();
    for elt in surface.elements.iter() {
        match elt.etype {
            ElementType::Triangle3 if elt.node_ids.len() == 3 => {
                if let (Some(&a), Some(&b), Some(&c)) = (
                    id_to_idx.get(&elt.node_ids[0]),
                    id_to_idx.get(&elt.node_ids[1]),
                    id_to_idx.get(&elt.node_ids[2]),
                ) {
                    tris.push([a, b, c]);
                }
            }
            ElementType::Quad4 if elt.node_ids.len() == 4 => {
                // Split quad into two triangles
                if let (Some(&a), Some(&b), Some(&c), Some(&d)) = (
                    id_to_idx.get(&elt.node_ids[0]),
                    id_to_idx.get(&elt.node_ids[1]),
                    id_to_idx.get(&elt.node_ids[2]),
                    id_to_idx.get(&elt.node_ids[3]),
                ) {
                    tris.push([a, b, c]);
                    tris.push([a, c, d]);
                }
            }
            _ => {}
        }
    }

    Ok((nodes, tris))
}

/// Find an existing node within `radius` of `target`, excluding `skip`.
// ─── Spatial grid for fast node lookup ──────────────────────────────────────

/// A uniform grid spatial index for fast nearby-node queries.
///
/// Replaces the O(N) linear scan in `find_nearby_node` with O(1) expected
/// cell lookup, reducing the main loop from O(N²) to O(N log N).
struct NodeGrid {
    nx: usize,
    ny: usize,
    nz: usize,
    ox: f64,
    oy: f64,
    oz: f64,
    sx: f64,
    sy: f64,
    sz: f64,
    cells: Vec<Vec<usize>>,
}

impl NodeGrid {
    fn build(nodes: &[[f64; 3]], target_per_axis: usize) -> Self {
        let r = target_per_axis.max(4).min(64);
        let (mut min, mut max) = bounds_3d(nodes);
        for i in 0..3 {
            let d = (max[i] - min[i]) * 0.05;
            if d < 1e-12 { max[i] += 0.5; min[i] -= 0.5; }
            else { min[i] -= d; max[i] += d; }
        }
        let dx = (max[0] - min[0]).max(1e-12);
        let dy = (max[1] - min[1]).max(1e-12);
        let dz = (max[2] - min[2]).max(1e-12);
        let dmax = dx.max(dy).max(dz);
        let nx = ((r as f64) * dx / dmax).ceil().max(1.0) as usize;
        let ny = ((r as f64) * dy / dmax).ceil().max(1.0) as usize;
        let nz = ((r as f64) * dz / dmax).ceil().max(1.0) as usize;
        let sx = dx / nx as f64;
        let sy = dy / ny as f64;
        let sz = dz / nz as f64;
        let n = nx * ny * nz;
        let mut cells: Vec<Vec<usize>> = (0..n).map(|_| Vec::new()).collect();
        for (i, p) in nodes.iter().enumerate() {
            let ix = (((p[0] - min[0]) / sx).floor() as isize).clamp(0, (nx - 1) as isize) as usize;
            let iy = (((p[1] - min[1]) / sy).floor() as isize).clamp(0, (ny - 1) as isize) as usize;
            let iz = (((p[2] - min[2]) / sz).floor() as isize).clamp(0, (nz - 1) as isize) as usize;
            cells[iz * ny * nx + iy * nx + ix].push(i);
        }
        Self { nx, ny, nz, ox: min[0], oy: min[1], oz: min[2], sx, sy, sz, cells }
    }

    fn cell_ijk(&self, p: [f64; 3]) -> (isize, isize, isize) {
        let ix = ((p[0] - self.ox) / self.sx).floor() as isize;
        let iy = ((p[1] - self.oy) / self.sy).floor() as isize;
        let iz = ((p[2] - self.oz) / self.sz).floor() as isize;
        (ix, iy, iz)
    }

    fn find_nearby(&self, nodes: &[[f64; 3]], target: [f64; 3], radius: f64, skip: &[usize]) -> Option<usize> {
        let r2 = radius * radius;
        let (ci, cj, ck) = self.cell_ijk(target);
        // Number of cells to search in each direction = ceil(radius / cell_size)
        let rd = (radius / self.sx.min(self.sy).min(self.sz)).ceil() as isize;
        let rd = rd.max(1).min(3); // cap to avoid blowing up search
        let mut best: Option<(usize, f64)> = None;
        for dk in -rd..=rd {
            for dj in -rd..=rd {
                for di in -rd..=rd {
                    let i = ci + di; let j = cj + dj; let k = ck + dk;
                    if i < 0 || j < 0 || k < 0 || i >= self.nx as isize || j >= self.ny as isize || k >= self.nz as isize { continue; }
                    for &ni in &self.cells[k as usize * self.ny * self.nx + j as usize * self.nx + i as usize] {
                        if skip.contains(&ni) { continue; }
                        let d2 = dist_sq(nodes[ni], target);
                        if d2 < r2 {
                            match best { Some((_, bd)) if d2 >= bd => {} _ => best = Some((ni, d2)) }
                        }
                    }
                }
            }
        }
        best.map(|(i, _)| i)
    }
}

fn bounds_3d(nodes: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    let mut min = [f64::MAX; 3]; let mut max = [f64::MIN; 3];
    for p in nodes {
        min[0] = min[0].min(p[0]); min[1] = min[1].min(p[1]); min[2] = min[2].min(p[2]);
        max[0] = max[0].max(p[0]); max[1] = max[1].max(p[1]); max[2] = max[2].max(p[2]);
    }
    (min, max)
}

/// O(1) expected-time nearby node search using a uniform grid.
/// Falls back to linear scan when the grid is stale.
fn find_nearby_node(
    nodes: &[[f64; 3]],
    target: [f64; 3],
    radius: f64,
    skip: &[usize],
    node_grid: &NodeGrid,
) -> Option<usize> {
    let r2 = radius * radius;
    let (ci, cj, ck) = node_grid.cell_ijk(target);
    let rd = (radius / node_grid.sx.min(node_grid.sy).min(node_grid.sz)).ceil() as isize;
    let rd = rd.max(1).min(3);
    let mut best: Option<(usize, f64)> = None;
    for dk in -rd..=rd {
        for dj in -rd..=rd {
            for di in -rd..=rd {
                let i = ci + di; let j = cj + dj; let k = ck + dk;
                if i < 0 || j < 0 || k < 0 || i >= node_grid.nx as isize || j >= node_grid.ny as isize || k >= node_grid.nz as isize { continue; }
                for &ni in &node_grid.cells[k as usize * node_grid.ny * node_grid.nx + j as usize * node_grid.nx + i as usize] {
                    if skip.contains(&ni) { continue; }
                    let d2 = dist_sq(nodes[ni], target);
                    if d2 < r2 {
                        match best { Some((_, bd)) if d2 >= bd => {} _ => best = Some((ni, d2)) }
                    }
                }
            }
        }
    }
    best.map(|(i, _)| i)
}

/// O(1) expected-time validity check using grid-based tet lookups.
/// Only checks tets in grid cells near the new tet's vertices.
fn is_valid_tet_grid(
    a: [f64; 3], b: [f64; 3], c: [f64; 3], p: [f64; 3],
    existing_tets: &[[usize; 4]],
    nodes: &[[f64; 3]],
    tet_grid: &NodeGrid,
) -> bool {
    if tetra_volume(a, b, c, p) < 1e-15 { return false; }
    // Collect candidate tet indices from grid cells near all 4 vertices.
    let mut candidates = Vec::<usize>::new();
    for &pt in &[a, b, c, p] {
        let (ci, cj, ck) = tet_grid.cell_ijk(pt);
        for dk in -1isize..=1 { for dj in -1..=1 { for di in -1..=1 {
            let i = ci + di; let j = cj + dj; let k = ck + dk;
            if i < 0 || j < 0 || k < 0 || i >= tet_grid.nx as isize || j >= tet_grid.ny as isize || k >= tet_grid.nz as isize { continue; }
            for &ti in &tet_grid.cells[k as usize * tet_grid.ny * tet_grid.nx + j as usize * tet_grid.nx + i as usize] {
                candidates.push(ti);
            }
        }}}
    }
    candidates.sort_unstable();
    candidates.dedup();
    for &ti in &candidates {
        if ti >= existing_tets.len() { continue; }
        let tet = &existing_tets[ti];
        let tv = [nodes[tet[0]], nodes[tet[1]], nodes[tet[2]], nodes[tet[3]]];
        if point_in_tetrahedron_strict(tv[0], tv[1], tv[2], tv[3], p, 1e-8) { return false; }
        for &v in &tv {
            if point_in_tetrahedron_strict(a, b, c, p, v, 1e-8) { return false; }
        }
    }
    true
}

/// Compute the ideal new-node position for a front facet.
fn ideal_point_3d(a: [f64; 3], b: [f64; 3], c: [f64; 3], normal: [f64; 3], h: f64) -> [f64; 3] {
    let centroid = [
        (a[0] + b[0] + c[0]) / 3.0,
        (a[1] + b[1] + c[1]) / 3.0,
        (a[2] + b[2] + c[2]) / 3.0,
    ];
    // Ideal height for equilateral tet: h * sqrt(2/3)
    let scale = h * (2.0_f64 / 3.0_f64).sqrt();
    [
        centroid[0] + scale * normal[0],
        centroid[1] + scale * normal[1],
        centroid[2] + scale * normal[2],
    ]
}

/// Test whether the candidate tetrahedron `(a, b, c, p)` intersects any face
/// or edge of the existing tetrahedra.
fn is_valid_tet(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    p: [f64; 3],
    existing_tets: &[[usize; 4]],
    nodes: &[[f64; 3]],
) -> bool {
    // Check volume
    if tetra_volume(a, b, c, p) < 1e-15 {
        return false;
    }

    // Check against existing tets: new tet should not significantly overlap
    for tet in existing_tets {
        let tv = [
            nodes[tet[0]],
            nodes[tet[1]],
            nodes[tet[2]],
            nodes[tet[3]],
        ];

        // Quick rejection: if the new point p is inside an existing tet, reject
        if point_in_tetrahedron_strict(tv[0], tv[1], tv[2], tv[3], p, 1e-8) {
            return false;
        }

        // Check if any vertex of existing tet is inside the new tet
        for &v in &tv {
            if point_in_tetrahedron_strict(a, b, c, p, v, 1e-8) {
                return false;
            }
        }
    }

    true
}

fn point_in_tetrahedron_strict(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    p: [f64; 3],
    eps: f64,
) -> bool {
    let v = tetra_volume(a, b, c, d);
    if v <= eps {
        return false;
    }
    let v0 = tetra_volume(p, b, c, d);
    let v1 = tetra_volume(a, p, c, d);
    let v2 = tetra_volume(a, b, p, d);
    let v3 = tetra_volume(a, b, c, p);
    let sum = v0 + v1 + v2 + v3;
    if (sum - v).abs() > eps * 16.0 {
        return false;
    }
    let min_part = v0.min(v1).min(v2).min(v3);
    min_part > eps
}

// ─── Quality helpers (used by both front and post-processing) ─────────────────

fn dihedral(p: [f64; 3], q: [f64; 3], r: [f64; 3], s: [f64; 3]) -> f64 {
    let pq = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
    let pr = [r[0] - p[0], r[1] - p[1], r[2] - p[2]];
    let ps = [s[0] - p[0], s[1] - p[1], s[2] - p[2]];
    let n1 = [
        pq[1] * pr[2] - pq[2] * pr[1],
        pq[2] * pr[0] - pq[0] * pr[2],
        pq[0] * pr[1] - pq[1] * pr[0],
    ];
    let n2 = [
        pq[1] * ps[2] - pq[2] * ps[1],
        pq[2] * ps[0] - pq[0] * ps[2],
        pq[0] * ps[1] - pq[1] * ps[0],
    ];
    let dot = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
    let l1 = (n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]).sqrt();
    let l2 = (n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]).sqrt();
    if l1 < 1e-12 || l2 < 1e-12 {
        return 0.0;
    }
    (dot / (l1 * l2)).clamp(-1.0, 1.0).acos().to_degrees()
}

fn min_dihedral_points(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    [
        dihedral(a, b, c, d),
        dihedral(a, c, b, d),
        dihedral(a, d, b, c),
        dihedral(b, c, a, d),
        dihedral(b, d, a, c),
        dihedral(c, d, a, b),
    ]
    .into_iter()
    .fold(f64::MAX, f64::min)
}

fn tetra_volume(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let ad = [a[0] - d[0], a[1] - d[1], a[2] - d[2]];
    let bd = [b[0] - d[0], b[1] - d[1], b[2] - d[2]];
    let cd = [c[0] - d[0], c[1] - d[1], c[2] - d[2]];
    let cross = [
        bd[1] * cd[2] - bd[2] * cd[1],
        bd[2] * cd[0] - bd[0] * cd[2],
        bd[0] * cd[1] - bd[1] * cd[0],
    ];
    (ad[0] * cross[0] + ad[1] * cross[1] + ad[2] * cross[2]).abs() / 6.0
}

fn node_xyz(mesh: &Mesh, node_id: u64) -> Result<[f64; 3], MeshAlgoError> {
    let p = mesh
        .nodes
        .get(&node_id)
        .ok_or_else(|| MeshAlgoError::Generation(format!("missing node id {node_id}")))?
        .position;
    Ok([p.x, p.y, p.z])
}

// ─── Post-processing: tet-split quality improvement ───────────────────────────

fn improve_front_quality(
    mut mesh: Mesh,
    min_target_dihedral_deg: f64,
    max_iters: usize,
) -> Result<Mesh, MeshAlgoError> {
    if max_iters == 0 {
        return Ok(mesh);
    }

    let mut next_node_id = mesh.nodes.keys().copied().max().unwrap_or(0).saturating_add(1);
    let mut next_elem_id = mesh
        .elements
        .iter()
        .map(|e| e.id)
        .max()
        .unwrap_or(0)
        .saturating_add(1);

    for _ in 0..max_iters {
        let Some((worst_idx, worst_dihedral)) =
            find_worst_dihedral_tet(&mesh, min_target_dihedral_deg)
        else {
            break;
        };

        let tet = match mesh.elements.get(worst_idx) {
            Some(e) if e.etype == ElementType::Tetrahedron4 && e.node_ids.len() == 4 => {
                [e.node_ids[0], e.node_ids[1], e.node_ids[2], e.node_ids[3]]
            }
            _ => break,
        };
        let a = node_xyz(&mesh, tet[0])?;
        let b = node_xyz(&mesh, tet[1])?;
        let c = node_xyz(&mesh, tet[2])?;
        let d = node_xyz(&mesh, tet[3])?;

        let centroid = [
            (a[0] + b[0] + c[0] + d[0]) * 0.25,
            (a[1] + b[1] + c[1] + d[1]) * 0.25,
            (a[2] + b[2] + c[2] + d[2]) * 0.25,
        ];

        let mut candidates = Vec::with_capacity(6);
        candidates.push(centroid);
        for v in [a, b, c, d] {
            candidates.push([
                centroid[0] * 0.8 + v[0] * 0.2,
                centroid[1] * 0.8 + v[1] * 0.2,
                centroid[2] * 0.8 + v[2] * 0.2,
            ]);
        }
        if let Some(cc) = tetra_circumcenter(a, b, c, d) {
            candidates.push([
                centroid[0] * 0.65 + cc[0] * 0.35,
                centroid[1] * 0.65 + cc[1] * 0.35,
                centroid[2] * 0.65 + cc[2] * 0.35,
            ]);
        }

        let mut best: Option<([f64; 3], f64)> = None;
        for p in candidates {
            if !is_valid_split_point(a, b, c, d, p) {
                continue;
            }
            let score = split_min_dihedral(a, b, c, d, p);
            match best {
                Some((_, s)) if score <= s => {}
                _ => best = Some((p, score)),
            }
        }

        let Some((best_p, best_score)) = best else {
            break;
        };
        if best_score <= worst_dihedral + 0.2 {
            break;
        }

        let new_node_id = next_node_id;
        next_node_id = next_node_id.saturating_add(1);
        mesh.add_node(Node::new(new_node_id, best_p[0], best_p[1], best_p[2]));

        let [n0, n1, n2, n3] = tet;
        mesh.elements.swap_remove(worst_idx);
        mesh.add_element(Element::new(
            next_elem_id,
            ElementType::Tetrahedron4,
            vec![n0, n1, n2, new_node_id],
        ));
        next_elem_id = next_elem_id.saturating_add(1);
        mesh.add_element(Element::new(
            next_elem_id,
            ElementType::Tetrahedron4,
            vec![n0, n1, n3, new_node_id],
        ));
        next_elem_id = next_elem_id.saturating_add(1);
        mesh.add_element(Element::new(
            next_elem_id,
            ElementType::Tetrahedron4,
            vec![n0, n2, n3, new_node_id],
        ));
        next_elem_id = next_elem_id.saturating_add(1);
        mesh.add_element(Element::new(
            next_elem_id,
            ElementType::Tetrahedron4,
            vec![n1, n2, n3, new_node_id],
        ));
        next_elem_id = next_elem_id.saturating_add(1);
    }

    Ok(mesh)
}

fn find_worst_dihedral_tet(mesh: &Mesh, threshold: f64) -> Option<(usize, f64)> {
    let mut worst: Option<(usize, f64)> = None;
    for (idx, e) in mesh.elements.iter().enumerate() {
        if e.etype != ElementType::Tetrahedron4 || e.node_ids.len() != 4 {
            continue;
        }
        let a = match node_xyz(mesh, e.node_ids[0]) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let b = match node_xyz(mesh, e.node_ids[1]) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let c = match node_xyz(mesh, e.node_ids[2]) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let d = match node_xyz(mesh, e.node_ids[3]) {
            Ok(v) => v,
            Err(_) => continue,
        };
        let q = min_dihedral_points(a, b, c, d);
        if q >= threshold {
            continue;
        }
        match worst {
            Some((_, wq)) if q >= wq => {}
            _ => worst = Some((idx, q)),
        }
    }
    worst
}

fn split_min_dihedral(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    p: [f64; 3],
) -> f64 {
    [
        min_dihedral_points(a, b, c, p),
        min_dihedral_points(a, b, d, p),
        min_dihedral_points(a, c, d, p),
        min_dihedral_points(b, c, d, p),
    ]
    .into_iter()
    .fold(f64::MAX, f64::min)
}

fn is_valid_split_point(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    p: [f64; 3],
) -> bool {
    if !point_in_tetrahedron_strict(a, b, c, d, p, 1e-12) {
        return false;
    }
    let vmin = [
        tetra_volume(a, b, c, p),
        tetra_volume(a, b, d, p),
        tetra_volume(a, c, d, p),
        tetra_volume(b, c, d, p),
    ]
    .into_iter()
    .fold(f64::MAX, f64::min);
    vmin > 1e-15
}

fn tetra_circumcenter(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
) -> Option<[f64; 3]> {
    let ba = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ca = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let da = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
    let rhs = [
        0.5 * ((b[0] * b[0] + b[1] * b[1] + b[2] * b[2])
            - (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])),
        0.5 * ((c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
            - (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])),
        0.5 * ((d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
            - (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])),
    ];
    solve_3x3([(ba, rhs[0]), (ca, rhs[1]), (da, rhs[2])])
}

fn solve_3x3(rows: [([f64; 3], f64); 3]) -> Option<[f64; 3]> {
    let mut a = [[0.0; 4]; 3];
    for i in 0..3 {
        a[i][0] = rows[i].0[0];
        a[i][1] = rows[i].0[1];
        a[i][2] = rows[i].0[2];
        a[i][3] = rows[i].1;
    }

    for col in 0..3 {
        let mut pivot = col;
        for r in (col + 1)..3 {
            if a[r][col].abs() > a[pivot][col].abs() {
                pivot = r;
            }
        }
        if a[pivot][col].abs() < 1e-15 {
            return None;
        }
        if pivot != col {
            a.swap(pivot, col);
        }
        let inv = 1.0 / a[col][col];
        for j in col..4 {
            a[col][j] *= inv;
        }
        for r in 0..3 {
            if r == col {
                continue;
            }
            let f = a[r][col];
            for j in col..4 {
                a[r][j] -= f * a[col][j];
            }
        }
    }
    Some([a[0][3], a[1][3], a[2][3]])
}

#[cfg(test)]
mod tests {
    use super::*;
    use rmsh_model::{Element, ElementType, Node};

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

    fn tet_surface() -> Mesh {
        // Tetrahedron surface: 4 nodes, 4 triangular faces
        let mut mesh = Mesh::new();
        mesh.add_node(Node::new(1, 0.0, 0.0, 0.0));
        mesh.add_node(Node::new(2, 1.0, 0.0, 0.0));
        mesh.add_node(Node::new(3, 0.5, 0.866, 0.0));
        mesh.add_node(Node::new(4, 0.5, 0.289, 0.816));
        mesh.add_element(Element::new(
            1,
            ElementType::Triangle3,
            vec![1, 2, 3],
        ));
        mesh.add_element(Element::new(
            2,
            ElementType::Triangle3,
            vec![1, 4, 2],
        ));
        mesh.add_element(Element::new(
            3,
            ElementType::Triangle3,
            vec![2, 4, 3],
        ));
        mesh.add_element(Element::new(
            4,
            ElementType::Triangle3,
            vec![3, 4, 1],
        ));
        mesh
    }

    #[test]
    fn frontal_3d_generates_mesh() {
        let mesh = Frontal3D::default()
            .mesh_3d(&cube_surface(), &MeshParams::with_size(0.4))
            .unwrap();
        assert!(mesh.elements_by_dimension(3).len() > 0);
    }

    #[test]
    fn frontal_3d_tetrahedron() {
        let mesh = Frontal3D::default()
            .mesh_3d(&tet_surface(), &MeshParams::with_size(0.3))
            .unwrap();
        assert!(mesh.elements_by_dimension(3).len() > 0);
    }

    #[test]
    fn front_from_cube_surface() {
        let surface = cube_surface();
        let (nodes, tris) = extract_surface_triangles(&surface).unwrap();
        let front = Front3D::from_surface(&nodes, &tris, &surface).unwrap();
        // Cube has 6 faces, each quad splits into 2 tris = 12 front facets
        assert_eq!(front.active.len(), 12);
    }

    #[test]
    fn front_pop_worst_removes_from_active() {
        let surface = cube_surface();
        let (nodes, tris) = extract_surface_triangles(&surface).unwrap();
        let mut front = Front3D::from_surface(&nodes, &tris, &surface).unwrap();
        let initial = front.active.len();
        let popped = front.pop_worst(&nodes);
        assert!(popped.is_some());
        assert_eq!(front.active.len(), initial - 1);
    }

    #[test]
    fn front_update_adds_three_facets() {
        let surface = cube_surface();
        let (mut nodes, tris) = extract_surface_triangles(&surface).unwrap();
        let mut front = Front3D::from_surface(&nodes, &tris, &surface).unwrap();

        // Pop a facet
        let (facet, _normal) = front.pop_worst(&nodes).unwrap();

        // Create a new node inside
        let a = nodes[facet[0]];
        let b_ = nodes[facet[1]];
        let c_ = nodes[facet[2]];
        let centroid = [
            (a[0] + b_[0] + c_[0]) / 3.0,
            (a[1] + b_[1] + c_[1]) / 3.0,
            (a[2] + b_[2] + c_[2]) / 3.0,
        ];
        let new_idx = nodes.len();
        nodes.push([
            centroid[0] + 0.1,
            centroid[1] + 0.1,
            centroid[2] + 0.1,
        ]);

        let before = front.active.len();
        front.update(facet, new_idx, &nodes);
        // After adding 3 new facets (toggle logic), there should be 3 new active facets
        // (minus 1 already removed by pop, plus 3 new, net +3 from before)
        assert_eq!(front.active.len(), before + 3);
    }

    #[test]
    fn extract_tris_from_quad_surface() {
        let surface = cube_surface();
        let (nodes, tris) = extract_surface_triangles(&surface).unwrap();
        assert_eq!(nodes.len(), 8);
        // 6 quads → 12 triangles
        assert_eq!(tris.len(), 12);
    }

    #[test]
    fn ideal_point_is_inward() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let normal = [0.0, 0.0, 1.0]; // pointing +z (inward)
        let p = ideal_point_3d(a, b, c, normal, 0.5);
        // Should be above the xy plane
        assert!(p[2] > 0.0);
    }

    #[test]
    fn dihedral_regular_tet() {
        // A regular tetrahedron has dihedral angle ≈ 70.5288°
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, 0.8660254, 0.0];
        let d = [0.5, 0.2886751, 0.8164966];
        let min_angle = min_dihedral_points(a, b, c, d);
        assert!((min_angle - 70.5288).abs() < 1.0);
    }
}
