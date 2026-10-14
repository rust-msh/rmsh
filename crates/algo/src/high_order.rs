//! Curved high-order element geometry (P5).
//!
//! Scope: the Gmsh MSH node-ordering reference layouts for the order 2-3
//! semantic element types, the isoparametric (nodal interpolation) maps
//! built on them, a validity check (`.valid`: the signed Jacobian measure
//! sampled over the reference nodes and a deterministic interior refinement
//! grid must stay positive), and the evaluation primitive the geometric
//! error contract is pinned with (sample the map, compare against the true
//! geometry — acceptance: error <= tolerance).
//!
//! # Reference layout provenance (discipline: empirical pins, not docs)
//!
//! Every layout below was recovered from meshes written by the real Gmsh
//! 4.14 on this machine (probe artifacts and scripts in `D:/tmp/rmsh-p5`,
//! 2026-10-13): single-element transfinite/extruded meshes for the tensor
//! families, first-element affine recovery for tets from a regular cube
//! mesh, plane-classification for pyramids. The recovery inverts the
//! corner frame of each element, so a layout entry is Gmsh's exact
//! reference position for that connectivity slot. The `probe_fixtures`
//! test embeds the verbatim probe coordinates and requires the map to
//! reproduce them — a swapped pair of slots cannot pass it.
//!
//! Layout digest (beyond the textbook families):
//! - Hex edge-node order groups by SOURCE vertex ascending:
//!   (0,1),(0,3),(0,4),(1,2),(1,5),(2,3),(2,6),(3,7),(4,5),(4,7),(5,6),(6,7);
//!   Hex27 face centers: bottom, front(y=0), left(x=0), right(x=1),
//!   back(y=1), top, then the body center.
//! - Tet edge-node order: (0,1),(1,2),(0,2),(0,3),(2,3),(1,3); Tet20 face
//!   interiors: (0,1,2),(0,1,3),(0,2,3),(1,2,3); the order-3 tet has no
//!   body node (the P3 barycentric lattice has none).
//! - Prism edge nodes group by source vertex within the bottom triangle,
//!   then vertical edges, then the top triangle's edges.
//!
//! # Validity semantics (`.valid`)
//!
//! The sampled Jacobian measure is SIGNED per dimension (tangent chord
//! projection for lines, corner-normal-projected area for surfaces,
//! determinant for volumes) and, for volumes, measured RELATIVE to the
//! sign of the element's own corner frame: real meshes contain
//! whole-element inverted orientations (real Gmsh 4.14 box meshes carry
//! type-11 tets whose canonical-reference determinant is negative), and
//! the curved-mapping question is FOLDING, not global orientation. A
//! folded map flips the measure's sign at some sample or drives it
//! through zero; both fail `is_valid`.
//!
//! # Pyramids (the collapsed / truncated pyramidal family)
//!
//! The three pyramid types (Pyramid13/14/30) carry the collapsed nodal
//! basis: the monomial exponents derive from each pinned slot through the
//! collapsed-lattice map `(I, J, K) = (p(u - h/2), p(v - h/2), p(1-h))`,
//! evaluated at the collapsed coordinates `(2u-1, 2v-1, h)`, and the nodal
//! coefficients solve the Vandermonde system on the pinned layout (the
//! same numerical construction as the serendipity families). This is the
//! interpolation space Gmsh pairs with its pyramidal reference points —
//! `pointsGenerators.cpp gmshGenerateMonomialsPyramid` generates the same
//! exponent set topologically (corner / edge / face / interior
//! interpolation of the corner multi-indices `(0,0,p) (p,0,p) (p,p,p)
//! (0,p,p) (0,0,0)` in the collapsed coordinates; algorithm-level
//! reference, re-authored — Gmsh GPL, M. Bergot 2001 for the underlying
//! truncated pyramidal family). Pyramid13 = the order-2 set without the
//! base-interior exponent, Pyramid14/30 = the complete order-2/3 sets.
//! The slot orders were probe-pinned against real Gmsh 4.14 meshes
//! (artifacts in `D:/tmp/rmsh-p5`; the Pyramid13 probe is the
//! 2026-10-14 `pyr13b` re-probe — the archived `pyramid13.geo` carried
//! `Recombine3DAll = 1` and produced pure hexes; the base-interior order
//! of the 30 is the probe's +v-first walk, correcting the earlier
//! "Quad16 row-major" shorthand).
//! - Pyramid13/14: 4 base corners, apex, edge nodes grouped by source
//!   vertex [(0,1),(0,3),(0,4),(1,2),(1,4),(2,3),(2,4),(3,4)], then the
//!   base-centroid node (complete variant only).
//! - Pyramid30: the same 5 + 16 edge nodes (two per edge, ascending from
//!   the source vertex), 4 side-face interiors in face order
//!   [(0,1,4),(3,0,4),(1,2,4),(2,3,4)], 4 base interiors (the probe's
//!   +v-first walk), 1 volume interior (the section center at h = 1/3).

use std::sync::OnceLock;

use rmsh_model::ElementType;

// ─── Probe-pinned reference layouts (Gmsh MSH node ordering) ─────────────────

const LINE3_LAYOUT: [[f64; 3]; 3] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.5, 0.0, 0.0],
];
const LINE4_LAYOUT: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0],
];
const TRI6_LAYOUT: [[f64; 3]; 6] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.5, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.0, 0.5, 0.0],
];
const TRI9_LAYOUT: [[f64; 3]; 9] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 2.0 / 3.0, 0.0],
    [0.0, 2.0 / 3.0, 0.0],
    [0.0, 1.0 / 3.0, 0.0],
];
const TRI10_LAYOUT: [[f64; 3]; 10] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 2.0 / 3.0, 0.0],
    [0.0, 2.0 / 3.0, 0.0],
    [0.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 1.0 / 3.0, 0.0],
];
const QUAD8_LAYOUT: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.5, 0.0, 0.0],
    [1.0, 0.5, 0.0],
    [0.5, 1.0, 0.0],
    [0.0, 0.5, 0.0],
];
const QUAD9_LAYOUT: [[f64; 3]; 9] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.5, 0.0, 0.0],
    [1.0, 0.5, 0.0],
    [0.5, 1.0, 0.0],
    [0.0, 0.5, 0.0],
    [0.5, 0.5, 0.0],
];
const QUAD12_LAYOUT: [[f64; 3]; 12] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0],
    [1.0, 1.0 / 3.0, 0.0],
    [1.0, 2.0 / 3.0, 0.0],
    [2.0 / 3.0, 1.0, 0.0],
    [1.0 / 3.0, 1.0, 0.0],
    [0.0, 2.0 / 3.0, 0.0],
    [0.0, 1.0 / 3.0, 0.0],
];
const QUAD16_LAYOUT: [[f64; 3]; 16] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0],
    [1.0, 1.0 / 3.0, 0.0],
    [1.0, 2.0 / 3.0, 0.0],
    [2.0 / 3.0, 1.0, 0.0],
    [1.0 / 3.0, 1.0, 0.0],
    [0.0, 2.0 / 3.0, 0.0],
    [0.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 1.0 / 3.0, 0.0],
    [2.0 / 3.0, 1.0 / 3.0, 0.0],
    [2.0 / 3.0, 2.0 / 3.0, 0.0],
    [1.0 / 3.0, 2.0 / 3.0, 0.0],
];
const TET10_LAYOUT: [[f64; 3]; 10] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.5, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.5],
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
];
const TET20_LAYOUT: [[f64; 3]; 20] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    // Edge nodes, grouped by source vertex: (0,1),(0,2),(0,3),(1,2),(1,3),(2,3).
    [1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 2.0 / 3.0, 0.0],
    [0.0, 1.0 / 3.0, 0.0],
    [0.0, 2.0 / 3.0, 0.0],
    [0.0, 0.0, 1.0 / 3.0],
    [0.0, 0.0, 2.0 / 3.0],
    [0.0, 1.0 / 3.0, 2.0 / 3.0],
    [0.0, 2.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 0.0, 2.0 / 3.0],
    [2.0 / 3.0, 0.0, 1.0 / 3.0],
    // Face interiors: (0,1,2), (0,1,3), (0,2,3), (1,2,3).
    [1.0 / 3.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 0.0, 1.0 / 3.0],
    [0.0, 1.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
];
const HEX20_LAYOUT: [[f64; 3]; 20] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
    // Edge nodes grouped by source vertex:
    // (0,1),(0,3),(0,4),(1,2),(1,5),(2,3),(2,6),(3,7),(4,5),(4,7),(5,6),(6,7).
    [0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.5],
    [1.0, 0.5, 0.0],
    [1.0, 0.0, 0.5],
    [0.5, 1.0, 0.0],
    [1.0, 1.0, 0.5],
    [0.0, 1.0, 0.5],
    [0.5, 0.0, 1.0],
    [0.0, 0.5, 1.0],
    [1.0, 0.5, 1.0],
    [0.5, 1.0, 1.0],
];
const HEX27_LAYOUT: [[f64; 3]; 27] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
    [0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.5],
    [1.0, 0.5, 0.0],
    [1.0, 0.0, 0.5],
    [0.5, 1.0, 0.0],
    [1.0, 1.0, 0.5],
    [0.0, 1.0, 0.5],
    [0.5, 0.0, 1.0],
    [0.0, 0.5, 1.0],
    [1.0, 0.5, 1.0],
    [0.5, 1.0, 1.0],
    // Face centers: bottom, top, front(y=0), back(y=1), left(x=0), right(x=1).
    [0.5, 0.5, 0.0],
    [0.5, 0.5, 1.0],
    [0.5, 0.0, 0.5],
    [0.5, 1.0, 0.5],
    [0.0, 0.5, 0.5],
    [1.0, 0.5, 0.5],
    // Body center.
    [0.5, 0.5, 0.5],
];
const PRISM15_LAYOUT: [[f64; 3]; 15] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    // Bottom triangle edges: (0,1),(0,2); vertical (0,3); (1,2),(1,4),(2,5).
    [0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.5],
    [0.5, 0.5, 0.0],
    [1.0, 0.0, 0.5],
    [0.0, 1.0, 0.5],
    // Top triangle edges: (3,4),(3,5),(4,5).
    [0.5, 0.0, 1.0],
    [0.0, 0.5, 1.0],
    [0.5, 0.5, 1.0],
];
const PRISM18_LAYOUT: [[f64; 3]; 18] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    [0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.5],
    [0.5, 0.5, 0.0],
    [1.0, 0.0, 0.5],
    [0.0, 1.0, 0.5],
    [0.5, 0.0, 1.0],
    [0.0, 0.5, 1.0],
    [0.5, 0.5, 1.0],
    // Side-face centers: (0,1,4,3), (0,2,5,3), (1,2,5,4).
    [0.5, 0.0, 0.5],
    [0.0, 0.5, 0.5],
    [0.5, 0.5, 0.5],
];
const PRISM40_LAYOUT: [[f64; 3]; 40] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    // Bottom triangle edges: (0,1),(0,2); vertical (0,3); (1,2),(1,4),(2,5).
    [1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0],
    [0.0, 1.0 / 3.0, 0.0],
    [0.0, 2.0 / 3.0, 0.0],
    [0.0, 0.0, 1.0 / 3.0],
    [0.0, 0.0, 2.0 / 3.0],
    [2.0 / 3.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 2.0 / 3.0, 0.0],
    [1.0, 0.0, 1.0 / 3.0],
    [1.0, 0.0, 2.0 / 3.0],
    [0.0, 1.0, 1.0 / 3.0],
    [0.0, 1.0, 2.0 / 3.0],
    // Top triangle edges: (3,4),(3,5),(4,5).
    [1.0 / 3.0, 0.0, 1.0],
    [2.0 / 3.0, 0.0, 1.0],
    [0.0, 1.0 / 3.0, 1.0],
    [0.0, 2.0 / 3.0, 1.0],
    [2.0 / 3.0, 1.0 / 3.0, 1.0],
    [1.0 / 3.0, 2.0 / 3.0, 1.0],
    // Bottom / top triangle centroids.
    [1.0 / 3.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 1.0 / 3.0, 1.0],
    // Side-face interiors at t=1/3, 2/3: (0,1,4,3), (0,2,5,3), (1,2,5,4).
    [1.0 / 3.0, 0.0, 1.0 / 3.0],
    [2.0 / 3.0, 0.0, 1.0 / 3.0],
    [2.0 / 3.0, 0.0, 2.0 / 3.0],
    [1.0 / 3.0, 0.0, 2.0 / 3.0],
    [0.0, 1.0 / 3.0, 1.0 / 3.0],
    [0.0, 1.0 / 3.0, 2.0 / 3.0],
    [0.0, 2.0 / 3.0, 2.0 / 3.0],
    [0.0, 2.0 / 3.0, 1.0 / 3.0],
    [2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    [2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0],
    // Volume interiors (the tri centroid at t=1/3 and t=2/3).
    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0],
];

const HEX64_LAYOUT: [[f64; 3]; 64] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
    // Edge nodes grouped by source vertex (same 12-edge order as Hex20),
    // two per edge ascending from the source.
    [1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0],
    [0.0, 1.0 / 3.0, 0.0],
    [0.0, 2.0 / 3.0, 0.0],
    [0.0, 0.0, 1.0 / 3.0],
    [0.0, 0.0, 2.0 / 3.0],
    [1.0, 1.0 / 3.0, 0.0],
    [1.0, 2.0 / 3.0, 0.0],
    [1.0, 0.0, 1.0 / 3.0],
    [1.0, 0.0, 2.0 / 3.0],
    [2.0 / 3.0, 1.0, 0.0],
    [1.0 / 3.0, 1.0, 0.0],
    [1.0, 1.0, 1.0 / 3.0],
    [1.0, 1.0, 2.0 / 3.0],
    [0.0, 1.0, 1.0 / 3.0],
    [0.0, 1.0, 2.0 / 3.0],
    [1.0 / 3.0, 0.0, 1.0],
    [2.0 / 3.0, 0.0, 1.0],
    [0.0, 1.0 / 3.0, 1.0],
    [0.0, 2.0 / 3.0, 1.0],
    [1.0, 1.0 / 3.0, 1.0],
    [1.0, 2.0 / 3.0, 1.0],
    [2.0 / 3.0, 1.0, 1.0],
    [1.0 / 3.0, 1.0, 1.0],
    // Face interiors, row-major per face: bottom, y=0, x=0, x=1, y=1, top.
    [1.0 / 3.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 2.0 / 3.0, 0.0],
    [2.0 / 3.0, 2.0 / 3.0, 0.0],
    [2.0 / 3.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 0.0, 1.0 / 3.0],
    [2.0 / 3.0, 0.0, 1.0 / 3.0],
    [2.0 / 3.0, 0.0, 2.0 / 3.0],
    [1.0 / 3.0, 0.0, 2.0 / 3.0],
    [0.0, 1.0 / 3.0, 1.0 / 3.0],
    [0.0, 1.0 / 3.0, 2.0 / 3.0],
    [0.0, 2.0 / 3.0, 2.0 / 3.0],
    [0.0, 2.0 / 3.0, 1.0 / 3.0],
    [1.0, 1.0 / 3.0, 1.0 / 3.0],
    [1.0, 2.0 / 3.0, 1.0 / 3.0],
    [1.0, 2.0 / 3.0, 2.0 / 3.0],
    [1.0, 1.0 / 3.0, 2.0 / 3.0],
    [2.0 / 3.0, 1.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0, 2.0 / 3.0],
    [2.0 / 3.0, 1.0, 2.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 1.0],
    [2.0 / 3.0, 1.0 / 3.0, 1.0],
    [2.0 / 3.0, 2.0 / 3.0, 1.0],
    [1.0 / 3.0, 2.0 / 3.0, 1.0],
    // Volume interiors, row-major per z level (1/3 then 2/3).
    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    [2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    [2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0],
    [2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0],
    [2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
];

// ─── Shape machinery ──────────────────────────────────────────────────────────

/// The interpolation-basis kind of a curved element.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ShapeFamily {
    /// Line Lagrange on the equispaced i/p lattice.
    Line,
    /// Tensor product of 1-D Lagrange factors (quads, hexes).
    Tensor,
    /// Barycentric lattice Lagrange (triangles, tetrahedra).
    Simplex,
    /// Barycentric triangle lattice x line Lagrange (complete prisms).
    Prism,
    /// Nodal basis numerically constructed on the pinned layout over the
    /// standard serendipity monomial space (quad8/12, hex20, prism15).
    Serendipity,
    /// Nodal basis numerically constructed on the pinned layout over the
    /// collapsed (truncated pyramidal / Bergot-type) monomial space — the
    /// pyramids. The monomial exponents live in the COLLAPSED lattice
    /// (see `lattice_indices`), evaluated at the collapsed coordinates
    /// (2u-1, 2v-1, h): the same interpolation space Gmsh's
    /// `gmshGenerateMonomialsPyramid` pairs with its reference points.
    Pyramid,
}

#[derive(Debug)]
struct TypeInfo {
    etype: ElementType,
    family: ShapeFamily,
    order: u8,
    dim: u8,
    layout: &'static [[f64; 3]],
    /// Serendipity monomial exponents (a, b, c) over the reference axes;
    /// empty for the closed-form families.
    monomials: &'static [[u8; 3]],
}

/// Triangle9 (incomplete cubic): the serendipity monomial set = the node
/// exponent pairs themselves, i.e. Gmsh's gmshGenerateMonomialsTriangle
/// serendip set (pointsGenerators.cpp): corners + edge exponents only, no
/// interior terms (the P3 lattice centroid bubble is absent).
const TRI9_MONOMIALS: &[[u8; 3]] = &[
    [0, 0, 0],
    [3, 0, 0],
    [0, 3, 0],
    [1, 0, 0],
    [2, 0, 0],
    [2, 1, 0],
    [1, 2, 0],
    [0, 2, 0],
    [0, 1, 0],
];

const QUAD8_MONOMIALS: &[[u8; 3]] = &[
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [2, 0, 0],
    [0, 2, 0],
    [2, 1, 0],
    [1, 2, 0],
];
const QUAD12_MONOMIALS: &[[u8; 3]] = &[
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [2, 0, 0],
    [0, 2, 0],
    [3, 0, 0],
    [0, 3, 0],
    [2, 1, 0],
    [1, 2, 0],
    [3, 1, 0],
    [1, 3, 0],
];
const HEX20_MONOMIALS: &[[u8; 3]] = &[
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 1, 0],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 2],
    [2, 1, 0],
    [1, 2, 0],
    [2, 0, 1],
    [1, 0, 2],
    [0, 2, 1],
    [0, 1, 2],
    [2, 1, 1],
    [1, 2, 1],
    [1, 1, 2],
];
const PRISM15_MONOMIALS: &[[u8; 3]] = &[
    // (a+b <= 2) x (c <= 1): the tri-P2 x line-P1 slab.
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
    [2, 0, 1],
    [0, 2, 1],
    // Vertical-edge corrections.
    [0, 0, 2],
    [1, 0, 2],
    [0, 1, 2],
];

/// The collapsed reference pyramid: base square [0,1]² at h=0 (CCW),
/// apex (0.5, 0.5, 1); the section at height h is the centered square of
/// half-width (1-h)/2, so the vertical edges are the straight corner-to-
/// apex lines. Slot order probe-pinned (Gmsh 4.14, artifacts in
/// `D:/tmp/rmsh-p5`): 4 base corners, apex, then the 8 edge nodes grouped
/// by source vertex [(0,1),(0,3),(0,4),(1,2),(1,4),(2,3),(2,4),(3,4)].
const PYRAMID13_LAYOUT: [[f64; 3]; 13] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.5, 0.5, 1.0],
    [0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0],
    [0.25, 0.25, 0.5],
    [1.0, 0.5, 0.0],
    [0.75, 0.25, 0.5],
    [0.5, 1.0, 0.0],
    [0.75, 0.75, 0.5],
    [0.25, 0.75, 0.5],
];
/// Pyramid14 = the incomplete layout plus the base-centroid node (the
/// complete variant's interior bubble).
const PYRAMID14_LAYOUT: [[f64; 3]; 14] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.5, 0.5, 1.0],
    [0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0],
    [0.25, 0.25, 0.5],
    [1.0, 0.5, 0.0],
    [0.75, 0.25, 0.5],
    [0.5, 1.0, 0.0],
    [0.75, 0.75, 0.5],
    [0.25, 0.75, 0.5],
    [0.5, 0.5, 0.0],
];
/// Pyramid30 (order 3, complete): corners + apex, 16 edge nodes (two per
/// edge, ascending from the source vertex), 4 side-face interiors (faces
/// [(0,1,4),(3,0,4),(1,2,4),(2,3,4)]), 4 base interiors (probe order:
/// (1/3,1/3), (1/3,2/3), (2/3,2/3), (2/3,1/3) — the +v-first walk the
/// real Gmsh mesh pins; NOT the Quad16 +u-first row), 1 volume interior
/// (the section center at h = 1/3).
const PYRAMID30_LAYOUT: [[f64; 3]; 30] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.5, 0.5, 1.0],
    // (0,1), (0,3)
    [1.0 / 3.0, 0.0, 0.0],
    [2.0 / 3.0, 0.0, 0.0],
    [0.0, 1.0 / 3.0, 0.0],
    [0.0, 2.0 / 3.0, 0.0],
    // (0,4)
    [1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0],
    // (1,2), (1,4)
    [1.0, 1.0 / 3.0, 0.0],
    [1.0, 2.0 / 3.0, 0.0],
    [5.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0],
    [2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0],
    // (2,3), (2,4)
    [2.0 / 3.0, 1.0, 0.0],
    [1.0 / 3.0, 1.0, 0.0],
    [5.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0],
    [2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    // (3,4)
    [1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0],
    [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    // side-face interiors, faces [(0,1,4),(3,0,4),(1,2,4),(2,3,4)]
    [0.5, 1.0 / 6.0, 1.0 / 3.0],
    [1.0 / 6.0, 0.5, 1.0 / 3.0],
    [5.0 / 6.0, 0.5, 1.0 / 3.0],
    [0.5, 5.0 / 6.0, 1.0 / 3.0],
    // base interiors (probe-pinned order)
    [1.0 / 3.0, 1.0 / 3.0, 0.0],
    [1.0 / 3.0, 2.0 / 3.0, 0.0],
    [2.0 / 3.0, 2.0 / 3.0, 0.0],
    [2.0 / 3.0, 1.0 / 3.0, 0.0],
    // volume interior
    [0.5, 0.5, 1.0 / 3.0],
];

const TYPE_INFOS: &[TypeInfo] = &[
    TypeInfo { etype: ElementType::Line3, family: ShapeFamily::Line, order: 2, dim: 1, layout: &LINE3_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Line4, family: ShapeFamily::Line, order: 3, dim: 1, layout: &LINE4_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Triangle6, family: ShapeFamily::Simplex, order: 2, dim: 2, layout: &TRI6_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Triangle9, family: ShapeFamily::Serendipity, order: 3, dim: 2, layout: &TRI9_LAYOUT, monomials: TRI9_MONOMIALS },
    TypeInfo { etype: ElementType::Triangle10, family: ShapeFamily::Simplex, order: 3, dim: 2, layout: &TRI10_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Quad8, family: ShapeFamily::Serendipity, order: 2, dim: 2, layout: &QUAD8_LAYOUT, monomials: QUAD8_MONOMIALS },
    TypeInfo { etype: ElementType::Quad9, family: ShapeFamily::Tensor, order: 2, dim: 2, layout: &QUAD9_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Quad12, family: ShapeFamily::Serendipity, order: 3, dim: 2, layout: &QUAD12_LAYOUT, monomials: QUAD12_MONOMIALS },
    TypeInfo { etype: ElementType::Quad16, family: ShapeFamily::Tensor, order: 3, dim: 2, layout: &QUAD16_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Tetrahedron10, family: ShapeFamily::Simplex, order: 2, dim: 3, layout: &TET10_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Tetrahedron20, family: ShapeFamily::Simplex, order: 3, dim: 3, layout: &TET20_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Hexahedron20, family: ShapeFamily::Serendipity, order: 2, dim: 3, layout: &HEX20_LAYOUT, monomials: HEX20_MONOMIALS },
    TypeInfo { etype: ElementType::Hexahedron27, family: ShapeFamily::Tensor, order: 2, dim: 3, layout: &HEX27_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Hexahedron64, family: ShapeFamily::Tensor, order: 3, dim: 3, layout: &HEX64_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Prism15, family: ShapeFamily::Serendipity, order: 2, dim: 3, layout: &PRISM15_LAYOUT, monomials: PRISM15_MONOMIALS },
    TypeInfo { etype: ElementType::Prism18, family: ShapeFamily::Prism, order: 2, dim: 3, layout: &PRISM18_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Prism40, family: ShapeFamily::Prism, order: 3, dim: 3, layout: &PRISM40_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Pyramid13, family: ShapeFamily::Pyramid, order: 2, dim: 3, layout: &PYRAMID13_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Pyramid14, family: ShapeFamily::Pyramid, order: 2, dim: 3, layout: &PYRAMID14_LAYOUT, monomials: &[] },
    TypeInfo { etype: ElementType::Pyramid30, family: ShapeFamily::Pyramid, order: 3, dim: 3, layout: &PYRAMID30_LAYOUT, monomials: &[] },
];

fn type_info(etype: ElementType) -> Option<&'static TypeInfo> {
    TYPE_INFOS.iter().find(|info| info.etype == etype)
}

/// The Gmsh MSH reference layout (natural coordinates in node order) of a
/// curved element type, empirically pinned against real Gmsh 4.14 output.
/// `None` for linear types and `ElementType::Unknown`.
pub fn reference_layout(etype: ElementType) -> Option<&'static [[f64; 3]]> {
    type_info(etype).map(|info| info.layout)
}

/// Whether the curved-element machinery covers this type: the 17 regular
/// Lagrange/serendipity families plus the three collapsed pyramids.
pub fn supports_curved_mapping(etype: ElementType) -> bool {
    type_info(etype).is_some()
}

/// g(k, t): the truncated-power lattice factor
/// Prod_{i<k} (p*t - i)/(k - i); g(0, .) = 1. The multiplicative building
/// block of the barycentric lattice bases; the 1-D Lagrange basis element
/// for lattice index k is g(k, t) * g(p-k, 1-t).
fn lattice_factor(k: usize, t: f64, p: usize) -> f64 {
    let pf = p as f64;
    let mut v = 1.0;
    for i in 0..k {
        v *= (pf * t - i as f64) / (k - i) as f64;
    }
    v
}

/// d/dt of [`lattice_factor`].
fn lattice_factor_grad(k: usize, t: f64, p: usize) -> f64 {
    let pf = p as f64;
    let mut acc = 0.0;
    for s in 0..k {
        let mut term = pf / (k - s) as f64;
        for i in 0..k {
            if i != s {
                term *= (pf * t - i as f64) / (k - i) as f64;
            }
        }
        acc += term;
    }
    acc
}

/// The 1-D Lagrange basis element for lattice index k on the equispaced
/// i/p lattice: g(k, t) * g(p-k, 1-t).
fn lagrange_1d(k: usize, t: f64, p: usize) -> f64 {
    lattice_factor(k, t, p) * lattice_factor(p - k, 1.0 - t, p)
}

/// d/dt of [`lagrange_1d`].
fn lagrange_1d_grad(k: usize, t: f64, p: usize) -> f64 {
    lattice_factor_grad(k, t, p) * lattice_factor(p - k, 1.0 - t, p)
        - lattice_factor(k, t, p) * lattice_factor_grad(p - k, 1.0 - t, p)
}

/// The per-node multi-index of a closed-form family, derived from the
/// pinned layout by exact lattice rounding (the layout coordinates are
/// exact multiples of 1/p, so rounding is exact).
fn lattice_indices(info: &TypeInfo) -> Vec<[usize; 4]> {
    let p = info.order as usize;
    info.layout
        .iter()
        .map(|r| {
            let i = (r[0] * p as f64).round() as usize;
            let j = (r[1] * p as f64).round() as usize;
            let k = (r[2] * p as f64).round() as usize;
            match info.family {
                ShapeFamily::Line => [i, 0, 0, 0],
                ShapeFamily::Tensor => [i, j, k, 0],
                ShapeFamily::Simplex if info.dim == 2 => [p - i - j, i, j, 0],
                ShapeFamily::Simplex => [p - i - j - k, i, j, k],
                ShapeFamily::Prism => [p - i - j, i, j, k],
                ShapeFamily::Serendipity => [0, 0, 0, 0],
                // The pyramids do not consume the lattice (the Bergot
                // route below); the slot carries the collapsed multi-index
                // for bookkeeping: (p*(u - h/2), p*(v - h/2), p*(1-h)).
                ShapeFamily::Pyramid => [
                    (r[0] * p as f64 - r[2] * p as f64 / 2.0).round() as usize,
                    (r[1] * p as f64 - r[2] * p as f64 / 2.0).round() as usize,
                    ((1.0 - r[2]) * p as f64).round() as usize,
                    0,
                ],
            }
        })
        .collect()
}

/// Gaussian elimination with partial pivoting; `None` on a singular matrix.
fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    for col in 0..n {
        let mut piv = col;
        for r in col + 1..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-300 {
            return None;
        }
        a.swap(piv, col);
        b.swap(piv, col);
        let d = a[col][col];
        for r in col + 1..n {
            let f = a[r][col] / d;
            if f != 0.0 {
                for c in col..n {
                    a[r][c] -= f * a[col][c];
                }
                b[r] -= f * b[col];
            }
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for c in r + 1..n {
            s -= a[r][c] * x[c];
        }
        x[r] = s / a[r][r];
    }
    Some(x)
}

fn eval_monomial(m: &[u8; 3], p: [f64; 3]) -> f64 {
    let mut v = 1.0;
    for (axis, &e) in m.iter().enumerate() {
        for _ in 0..e {
            v *= p[axis];
        }
    }
    v
}

fn eval_monomial_grad(m: &[u8; 3], p: [f64; 3], axis: usize) -> f64 {
    let e = m[axis] as i32;
    if e == 0 {
        return 0.0;
    }
    let mut v = e as f64;
    for (a2, &e2) in m.iter().enumerate() {
        let e2 = e2 as i32;
        let take = if a2 == axis { e2 - 1 } else { e2 };
        for _ in 0..take {
            v *= p[a2];
        }
    }
    v
}

/// A constructed nodal interpolation basis on the pinned reference layout.
struct Basis {
    info: &'static TypeInfo,
    lattice: Vec<[usize; 4]>,
    /// Monomial/orthogonal-construction families: coefficients[i][k] =
    /// weight of the k-th generator in N_i.
    coefficients: Vec<Vec<f64>>,
    /// Pyramids: the collapsed (truncated pyramidal / Bergot) generator
    /// set the coefficients weight.
    bergot: Option<BergotBasis>,
}

impl Basis {
    fn new(info: &'static TypeInfo) -> Basis {
        let bergot = match info.family {
            ShapeFamily::Pyramid => Some(BergotBasis::new(
                info.order as usize,
                matches!(info.etype, ElementType::Pyramid13),
            )),
            _ => None,
        };
        let lattice = lattice_indices(info);
        let coefficients = match info.family {
            ShapeFamily::Serendipity => {
                let n = info.layout.len();
                let mut a = vec![vec![0.0; n]; n];
                for (j, row) in a.iter_mut().enumerate() {
                    for (k, m) in info.monomials.iter().enumerate() {
                        row[k] = eval_monomial(m, info.layout[j]);
                    }
                }
                Some((a, n))
            }
            ShapeFamily::Pyramid => {
                // The Bergot generators inverted against the pinned
                // layout: Gmsh's pyramidalBasis construction (it inverts
                // the same Bergot-by-points matrix; the pure-polynomial
                // pairing was probed singular — the complete pyramidal
                // space is not polynomial).
                let bergot = bergot.as_ref().expect("pyramid carries bergot");
                let n = info.layout.len();
                let mut a = vec![vec![0.0; n]; n];
                for (j, row) in a.iter_mut().enumerate() {
                    let fvals = bergot.f(info.layout[j]);
                    for (k, fv) in fvals.iter().enumerate() {
                        row[k] = *fv;
                    }
                }
                Some((a, n))
            }
            _ => None,
        };
        let coefficients = coefficients
            .map(|(a, n)| {
                (0..n)
                    .map(|i| {
                        let mut e = vec![0.0; n];
                        e[i] = 1.0;
                        solve_linear(&a, &e).unwrap_or_else(|| {
                            panic!("interpolation space singular for {info:?}")
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();
        Basis {
            info,
            lattice,
            coefficients,
            bergot,
        }
    }

    /// N_i at the reference point.
    fn shape(&self, i: usize, p: [f64; 3]) -> f64 {
        match self.info.family {
            ShapeFamily::Line => lagrange_1d(self.lattice[i][0], p[0], self.info.order as usize),
            ShapeFamily::Tensor => {
                let l = self.lattice[i];
                let o = self.info.order as usize;
                lagrange_1d(l[0], p[0], o) * lagrange_1d(l[1], p[1], o) * lagrange_1d(l[2], p[2], o)
            }
            ShapeFamily::Simplex => {
                let dim = self.info.dim as usize;
                let lams = simplex_barycentric(dim, p);
                let mut v = 1.0;
                for m in 0..=dim {
                    v *= lattice_factor(self.lattice[i][m], lams[m], self.info.order as usize);
                }
                v
            }
            ShapeFamily::Prism => {
                let lams = simplex_barycentric(2, p);
                let mut v = 1.0;
                for m in 0..=2 {
                    v *= lattice_factor(self.lattice[i][m], lams[m], self.info.order as usize);
                }
                v *= lagrange_1d(self.lattice[i][3], p[2], self.info.order as usize);
                v
            }
            ShapeFamily::Serendipity => {
                let mut v = 0.0;
                for (k, m) in self.info.monomials.iter().enumerate() {
                    v += self.coefficients[i][k] * eval_monomial(m, p);
                }
                v
            }
            ShapeFamily::Pyramid => {
                let bergot = self.bergot.as_ref().expect("pyramid carries bergot");
                let fvals = bergot.f(p);
                let mut v = 0.0;
                for (k, fv) in fvals.iter().enumerate() {
                    v += self.coefficients[i][k] * fv;
                }
                v
            }
        }
    }

    /// d N_i / d(axis) at the reference point.
    fn shape_grad(&self, i: usize, p: [f64; 3], axis: usize) -> f64 {
        match self.info.family {
            ShapeFamily::Line => {
                lagrange_1d_grad(self.lattice[i][0], p[0], self.info.order as usize)
            }
            ShapeFamily::Tensor => {
                let l = self.lattice[i];
                let o = self.info.order as usize;
                let mut v = lagrange_1d_grad(l[axis], p[axis], o);
                for (a2, &li) in l.iter().enumerate() {
                    if a2 != axis && a2 < 3 {
                        v *= lagrange_1d(li, p[a2], o);
                    }
                }
                v
            }
            ShapeFamily::Simplex => {
                let dim = self.info.dim as usize;
                let o = self.info.order as usize;
                let lams = simplex_barycentric(dim, p);
                // dN/d(axis) = -dN/dlam0 + dN/dlam(axis+1)
                let mut terms = [0.0f64; 2];
                for (t, m) in [0usize, axis + 1].iter().enumerate() {
                    let mut v = lattice_factor_grad(self.lattice[i][*m], lams[*m], o);
                    for m2 in 0..=dim {
                        if m2 != *m {
                            v *= lattice_factor(self.lattice[i][m2], lams[m2], o);
                        }
                    }
                    terms[t] = v;
                }
                terms[1] - terms[0]
            }
            ShapeFamily::Prism => {
                let o = self.info.order as usize;
                let lams = simplex_barycentric(2, p);
                if axis < 2 {
                    // tri part chain + line factor
                    let mut terms = [0.0f64; 2];
                    for (t, m) in [0usize, axis + 1].iter().enumerate() {
                        let mut v = lattice_factor_grad(self.lattice[i][*m], lams[*m], o);
                        for m2 in 0..=2 {
                            if m2 != *m {
                                v *= lattice_factor(self.lattice[i][m2], lams[m2], o);
                            }
                        }
                        v *= lattice_factor(self.lattice[i][3], p[2], o);
                        terms[t] = v;
                    }
                    terms[1] - terms[0]
                } else {
                    let mut v = lagrange_1d_grad(self.lattice[i][3], p[2], o);
                    for m2 in 0..=2 {
                        v *= lattice_factor(self.lattice[i][m2], lams[m2], o);
                    }
                    v
                }
            }
            ShapeFamily::Serendipity => {
                let mut v = 0.0;
                for (k, m) in self.info.monomials.iter().enumerate() {
                    v += self.coefficients[i][k] * eval_monomial_grad(m, p, axis);
                }
                v
            }
            ShapeFamily::Pyramid => {
                // Chain rule: the Bergot gradients live in Gmsh's
                // collapsed coordinates (x', y', w) = (2u-1, 2v-1, h), so
                // the u/v reference derivatives pick up a factor 2.
                let bergot = self.bergot.as_ref().expect("pyramid carries bergot");
                let grads = bergot.df(p);
                let chain = if axis < 2 { 2.0 } else { 1.0 };
                let mut v = 0.0;
                for (k, g) in grads.iter().enumerate() {
                    v += self.coefficients[i][k] * g[axis];
                }
                chain * v
            }
        }
    }
}

// ─── The collapsed (truncated pyramidal / Bergot) generator set ───────────────
//
// The pyramid interpolation space is NOT polynomial: its natural
// generators carry (1-w)^max(i,j) factors and the complete family's
// monomial view is rational (Gmsh pyramidalBasis). The nodal basis is
// built by inverting the Bergot generator matrix against the pinned
// layout — the same construction Gmsh's pyramidalBasis performs (the
// pure-polynomial pairing was probed singular).
//
// Reference: Gmsh 4.14 numeric/BergotBasis.cpp + numeric/orthogonalBasis
// (the Legendre / Jacobi recurrences), algorithm-level port; M. Bergot
// (2001) for the underlying truncated pyramidal space.

/// Legendre polynomials L_0..L_n at u (Gmsh orthogonalBasis recurrence).
fn legendre_f(n: usize, u: f64) -> Vec<f64> {
    let mut val = vec![0.0; n + 1];
    val[0] = 1.0;
    for i in 0..n {
        let a1i = (i + 1) as f64;
        let a3i = 2.0 * i as f64 + 1.0;
        let a4i = i as f64;
        val[i + 1] = a3i * u * val[i];
        if i > 0 {
            val[i + 1] -= a4i * val[i - 1];
        }
        val[i + 1] /= a1i;
    }
    val
}

/// d/du of [`legendre_f`] (with the ±1 indeterminate forms).
fn legendre_df(n: usize, u: f64) -> Vec<f64> {
    let mut val = vec![0.0; n + 1];
    if u == 1.0 || u == -1.0 {
        for k in 0..=n {
            val[k] = 0.5 * k as f64 * (k as f64 + 1.0);
        }
        if u == -1.0 && n >= 2 {
            for k in (2..=n).step_by(2) {
                val[k] = -val[k];
            }
        }
        return val;
    }
    let tmp = legendre_f(n, u);
    for i in 1..=n {
        let g1 = -u * i as f64;
        let g0 = i as f64;
        let g2 = 1.0 - u * u;
        val[i] = (g1 * tmp[i] + g0 * tmp[i - 1]) / g2;
    }
    val
}

/// Jacobi polynomials P_0..P_n^(alpha, beta) at u (Gmsh recurrence).
fn jacobi_f(n: usize, alpha: f64, beta: f64, u: f64) -> Vec<f64> {
    let mut val = vec![0.0; n + 1];
    let alpha_plus_beta = alpha + beta;
    let a2_minus_b2 = alpha * alpha - beta * beta;
    val[0] = 1.0;
    if n >= 1 {
        val[1] = 0.5 * (2.0 * (alpha + 1.0) + (alpha_plus_beta + 2.0) * (u - 1.0));
    }
    for i in 1..n {
        let ii = i as f64;
        let two_i = 2.0 * ii;
        let a1i =
            2.0 * (ii + 1.0) * (ii + alpha_plus_beta + 1.0) * (two_i + alpha_plus_beta);
        let a2i = (two_i + alpha_plus_beta + 1.0) * a2_minus_b2;
        let a3i = {
            let x = two_i + alpha_plus_beta;
            x * (x + 1.0) * (x + 2.0)
        };
        let a4i = 2.0 * (ii + alpha) * (ii + beta) * (two_i + alpha_plus_beta + 2.0);
        val[i + 1] = ((a2i + a3i * u) * val[i] - a4i * val[i - 1]) / a1i;
    }
    val
}

/// d/du of [`jacobi_f`] (with the ±1 indeterminate forms).
fn jacobi_df(n: usize, alpha: f64, beta: f64, u: f64) -> Vec<f64> {
    let mut val = vec![0.0; n + 1];
    let alpha_plus_beta = alpha + beta;
    if u == 1.0 || u == -1.0 {
        // alpha or beta in the closed form, depending on the side.
        let coeff = if u == 1.0 { alpha } else { beta };
        let f_max = n.max(1) + coeff as usize;
        let mut fact = vec![1.0; f_max + 1];
        for i in 1..=f_max {
            fact[i] = i as f64 * fact[i - 1];
        }
        for k in 1..=n {
            val[k] = 0.5 * (k as f64 + alpha_plus_beta + 1.0) * fact[k + coeff as usize]
                / (fact[coeff as usize + 1] * fact[k - 1]);
        }
        if u == -1.0 && n >= 2 {
            for k in (2..=n).step_by(2) {
                val[k] = -val[k];
            }
        }
        return val;
    }
    let tmp = jacobi_f(n, alpha, beta, u);
    val[0] = 0.0;
    if n >= 1 {
        val[1] = 0.5 * (alpha_plus_beta + 2.0);
    }
    for i in 2..=n {
        let ii = i as f64;
        let aa = 2.0 * ii + alpha_plus_beta;
        let g2 = aa * (1.0 - u * u);
        let g1 = ii * (alpha - beta - aa * u);
        let g0 = 2.0 * (ii + alpha) * (ii + beta);
        val[i] = (g1 * tmp[i] + g0 * tmp[i - 1]) / g2;
    }
    val
}

/// The collapsed pyramidal generator set: f_(i,j,k) = L_i(x') L_j(y')
/// (1-w)^max(i,j) P_k^(2 max(i,j)+2, 0)(2w-1) over the collapsed
/// coordinates (x', y', w) = (2u-1, 2v-1, h), for (i, j) within the
/// validity rule and k <= order - max(i,j).
struct BergotBasis {
    order: usize,
    incomplete: bool,
}

impl BergotBasis {
    fn new(order: usize, incomplete: bool) -> Self {
        BergotBasis { order, incomplete }
    }

    /// The (i, j) validity rule (Gmsh BergotBasis::validIJ): the complete
    /// family keeps the full (i, j) square; the order-2 incomplete drops
    /// the functions beyond i + j <= order except the two single-width
    /// edge corrections.
    fn valid_ij(&self, i: usize, j: usize) -> bool {
        if !self.incomplete {
            return i <= self.order && j <= self.order;
        }
        if i + j <= self.order {
            return true;
        }
        if i + j == self.order + 1 {
            return i == 1 || j == 1;
        }
        false
    }

    /// The generator count (= the layout node count).
    fn count(&self) -> usize {
        let mut n = 0usize;
        for i in 0..=self.order {
            for j in 0..=self.order {
                if self.valid_ij(i, j) {
                    n += self.order - i.max(j) + 1;
                }
            }
        }
        n
    }

    /// All generators at the reference point `p` = (u, v, h).
    fn f(&self, p: [f64; 3]) -> Vec<f64> {
        let (u, v, w) = (2.0 * p[0] - 1.0, 2.0 * p[1] - 1.0, p[2]);
        let uhat = if w == 1.0 { 0.0 } else { u / (1.0 - w) };
        let vhat = if w == 1.0 { 0.0 } else { v / (1.0 - w) };
        let what = 2.0 * w - 1.0;
        let uf = legendre_f(self.order, uhat);
        let vf = legendre_f(self.order, vhat);
        let mut out = Vec::with_capacity(self.count());
        for i in 0..=self.order {
            for j in 0..=self.order {
                if !self.valid_ij(i, j) {
                    continue;
                }
                let m = i.max(j);
                let wf = jacobi_f(self.order - m, 2.0 * m as f64 + 2.0, 0.0, what);
                let fact = (1.0 - w).powi(m as i32);
                for k in 0..=(self.order - m) {
                    out.push(uf[i] * vf[j] * wf[k] * fact);
                }
            }
        }
        out
    }

    /// The generator gradients at `p`, in the collapsed coordinates
    /// (d/dx', d/dy', d/dw); the reference-axis chain factor is applied
    /// by the caller.
    fn df(&self, p: [f64; 3]) -> Vec<[f64; 3]> {
        let (u, v, w) = (2.0 * p[0] - 1.0, 2.0 * p[1] - 1.0, p[2]);
        let uhat = if w == 1.0 { 0.0 } else { u / (1.0 - w) };
        let vhat = if w == 1.0 { 0.0 } else { v / (1.0 - w) };
        let what = 2.0 * w - 1.0;
        let uf = legendre_f(self.order, uhat);
        let ug = legendre_df(self.order, uhat);
        let vf = legendre_f(self.order, vhat);
        let vg = legendre_df(self.order, vhat);
        let mut wfs = Vec::new();
        let mut wgs = Vec::new();
        for m in 0..=self.order {
            let k_max = self.order - m;
            wfs.push(jacobi_f(k_max, 2.0 * m as f64 + 2.0, 0.0, what));
            wgs.push(jacobi_df(k_max, 2.0 * m as f64 + 2.0, 0.0, what));
        }
        let mut out = Vec::with_capacity(self.count());
        for i in 0..=self.order {
            for j in 0..=self.order {
                if !self.valid_ij(i, j) {
                    continue;
                }
                let m = i.max(j);
                let (wf, wg) = (&wfs[m], &wgs[m]);
                // The collapsed-coordinate derivatives, with the
                // u/(1-w) / v/(1-w) indeterminate forms resolved per the
                // Gmsh df() special cases.
                if m == 0 {
                    for k in 0..=(self.order - m) {
                        out.push([0.0, 0.0, 2.0 * wg[k]]);
                    }
                } else if m == 1 {
                    if i == 0 {
                        for k in 0..=(self.order - m) {
                            out.push([0.0, wf[k], 2.0 * v * wg[k]]);
                        }
                    } else if j == 0 {
                        for k in 0..=(self.order - m) {
                            out.push([wf[k], 0.0, 2.0 * u * wg[k]]);
                        }
                    } else {
                        for k in 0..=(self.order - m) {
                            out.push([
                                vhat * wf[k],
                                uhat * wf[k],
                                uhat * vhat * wf[k] + 2.0 * uhat * v * wg[k],
                            ]);
                        }
                    }
                } else {
                    let o_m_w = 1.0 - w;
                    let pow_m2 = o_m_w.powi(m as i32 - 2);
                    let pow_m1 = pow_m2 * o_m_w;
                    for k in 0..=(self.order - m) {
                        out.push([
                            ug[i] * vf[j] * wf[k] * pow_m1,
                            uf[i] * vg[j] * wf[k] * pow_m1,
                            wf[k] * pow_m2
                                * (u * ug[i] * vf[j] + v * uf[i] * vg[j])
                                + uf[i] * vf[j] * pow_m1
                                    * (2.0 * o_m_w * wg[k] - m as f64 * wf[k]),
                        ]);
                    }
                }
            }
        }
        out
    }
}

fn simplex_barycentric(dim: usize, p: [f64; 3]) -> [f64; 4] {
    if dim == 2 {
        [1.0 - p[0] - p[1], p[0], p[1], 0.0]
    } else {
        [1.0 - p[0] - p[1] - p[2], p[0], p[1], p[2]]
    }
}

struct BasisCache {
    slots: [OnceLock<Basis>; TYPE_INFOS.len()],
}

fn basis_for(etype: ElementType) -> Option<&'static Basis> {
    static CACHE: OnceLock<BasisCache> = OnceLock::new();
    let (idx, info) = TYPE_INFOS
        .iter()
        .enumerate()
        .find(|(_, info)| info.etype == etype)?;
    let cache = CACHE.get_or_init(|| BasisCache {
        slots: Default::default(),
    });
    Some(cache.slots[idx].get_or_init(|| Basis::new(info)))
}

// ─── The curved element API ───────────────────────────────────────────────────

/// A curved high-order element: the pinned reference layout plus physical
/// node coordinates, with the isoparametric map, sampled validity and the
/// evaluation primitive behind the geometric-error contract.
pub struct CurvedElement<'a> {
    basis: &'static Basis,
    xyz: &'a [[f64; 3]],
}

impl<'a> CurvedElement<'a> {
    /// Build a curved element from node coordinates in MSH connectivity
    /// order. `None` when the type has no curved mapping (linear types,
    /// the deferred pyramids, `Unknown`) or the coordinate count does not
    /// match the type's node count.
    pub fn new(etype: ElementType, node_coordinates: &'a [[f64; 3]]) -> Option<Self> {
        let basis = basis_for(etype)?;
        if node_coordinates.len() != etype.node_count() {
            return None;
        }
        Some(CurvedElement {
            basis,
            xyz: node_coordinates,
        })
    }

    /// The element type.
    pub fn element_type(&self) -> ElementType {
        self.basis.info.etype
    }

    /// Evaluate the isoparametric map x(u, v, w).
    pub fn eval(&self, reference_point: [f64; 3]) -> [f64; 3] {
        let mut out = [0.0; 3];
        for (i, x) in self.xyz.iter().enumerate() {
            let n = self.basis.shape(i, reference_point);
            for d in 0..3 {
                out[d] += n * x[d];
            }
        }
        out
    }

    /// The signed Jacobian measure at a reference point: the tangent
    /// projection along the (0 -> 1) chord for lines, the projected area
    /// element for surface elements, the determinant for volume elements.
    /// A negative value means the map folds (reverses orientation); zero
    /// means it degenerates.
    pub fn jacobian_measure(&self, reference_point: [f64; 3]) -> f64 {
        let n = self.xyz.len();
        let mut grads = Vec::with_capacity(n);
        match self.basis.info.dim {
            1 => {
                grads.push(self.basis.shape_grad(0, reference_point, 0));
                let mut d = [0.0; 3];
                for (i, x) in self.xyz.iter().enumerate() {
                    let g = self.basis.shape_grad(i, reference_point, 0);
                    for k in 0..3 {
                        d[k] += g * x[k];
                    }
                }
                let chord = [
                    self.xyz[1][0] - self.xyz[0][0],
                    self.xyz[1][1] - self.xyz[0][1],
                    self.xyz[1][2] - self.xyz[0][2],
                ];
                dot3(d, normalize3(chord))
            }
            2 => {
                let mut ju = [0.0; 3];
                let mut jv = [0.0; 3];
                for (i, x) in self.xyz.iter().enumerate() {
                    let gu = self.basis.shape_grad(i, reference_point, 0);
                    let gv = self.basis.shape_grad(i, reference_point, 1);
                    for k in 0..3 {
                        ju[k] += gu * x[k];
                        jv[k] += gv * x[k];
                    }
                }
                let e01 = sub3(self.xyz[1], self.xyz[0]);
                let e02 = sub3(self.xyz[2], self.xyz[0]);
                let normal = normalize3(cross3(e01, e02));
                dot3(cross3(ju, jv), normal)
            }
            _ => {
                let mut cols = [[0.0f64; 3]; 3];
                for (i, x) in self.xyz.iter().enumerate() {
                    for axis in 0..3 {
                        let g = self.basis.shape_grad(i, reference_point, axis);
                        for k in 0..3 {
                            cols[axis][k] += g * x[k];
                        }
                    }
                }
                // Real meshes contain whole-element inverted orientations
                // (Gmsh writes them; see the module docs): validity is
                // about FOLDING, so the determinant is measured relative
                // to the sign of the element's own corner frame. The frame
                // edges follow the family: (0,1),(0,2),(0,3) for tets and
                // prisms, the (0,1),(0,3),(0,4) reference edges for hexes
                // and pyramids (nodes 0-3 of both are one face and would
                // be planar).
                let (e1, e2, e3) = match self.basis.info.etype {
                    ElementType::Hexahedron20
                    | ElementType::Hexahedron27
                    | ElementType::Hexahedron64
                    | ElementType::Pyramid13
                    | ElementType::Pyramid14
                    | ElementType::Pyramid30 => (
                        sub3(self.xyz[1], self.xyz[0]),
                        sub3(self.xyz[3], self.xyz[0]),
                        sub3(self.xyz[4], self.xyz[0]),
                    ),
                    _ => (
                        sub3(self.xyz[1], self.xyz[0]),
                        sub3(self.xyz[2], self.xyz[0]),
                        sub3(self.xyz[3], self.xyz[0]),
                    ),
                };
                let corner = det3([e1, e2, e3]);
                let s = if corner > 0.0 {
                    1.0
                } else if corner < 0.0 {
                    -1.0
                } else {
                    0.0
                };
                s * det3(cols)
            }
        }
    }

    /// The deterministic validity sample set: every reference node plus a
    /// fixed interior refinement lattice (order+1 cells per parametric
    /// dimension, sampled at their midpoints).
    pub fn validity_samples(&self) -> Vec<[f64; 3]> {
        let info = self.basis.info;
        let p = info.order as usize;
        let mut samples: Vec<[f64; 3]> = info.layout.to_vec();
        match info.family {
            ShapeFamily::Line => {
                for j in 0..p + 1 {
                    let t = (j as f64 + 0.5) / (p + 1) as f64;
                    samples.push([t, 0.0, 0.0]);
                }
            }
            ShapeFamily::Tensor => match info.dim {
                2 => {
                    for j in 0..p + 1 {
                        for k in 0..p + 1 {
                            samples.push([(
                                (j as f64 + 0.5) / (p + 1) as f64
                            ), ((k as f64 + 0.5) / (p + 1) as f64), 0.0]);
                        }
                    }
                }
                _ => {
                    for j in 0..p + 1 {
                        for k in 0..p + 1 {
                            for m in 0..p + 1 {
                                samples.push([
                                    (j as f64 + 0.5) / (p + 1) as f64,
                                    (k as f64 + 0.5) / (p + 1) as f64,
                                    (m as f64 + 0.5) / (p + 1) as f64,
                                ]);
                            }
                        }
                    }
                }
            },
            ShapeFamily::Simplex if info.dim == 2 => {
                for a in 0..=p + 1 {
                    for b in 0..=p + 1 - a {
                        let c = p + 1 - a - b;
                        samples.push([a as f64 / (p + 1) as f64, b as f64 / (p + 1) as f64, 0.0]);
                        let _ = c;
                    }
                }
            }
            ShapeFamily::Simplex => {
                for a in 0..=p + 1 {
                    for b in 0..=p + 1 - a {
                        for c in 0..=p + 1 - a - b {
                            let d = p + 1 - a - b - c;
                            samples.push([
                                a as f64 / (p + 1) as f64,
                                b as f64 / (p + 1) as f64,
                                c as f64 / (p + 1) as f64,
                            ]);
                            let _ = d;
                        }
                    }
                }
            }
            ShapeFamily::Prism => {
                for a in 0..=p + 1 {
                    for b in 0..=p + 1 - a {
                        for k in 0..p + 1 {
                            samples.push([
                                a as f64 / (p + 1) as f64,
                                b as f64 / (p + 1) as f64,
                                (k as f64 + 0.5) / (p + 1) as f64,
                            ]);
                        }
                    }
                }
            }
            ShapeFamily::Serendipity => match info.dim {
                2 => {
                    for j in 0..p + 1 {
                        for k in 0..p + 1 {
                            samples.push([
                                (j as f64 + 0.5) / (p + 1) as f64,
                                (k as f64 + 0.5) / (p + 1) as f64,
                                0.0,
                            ]);
                        }
                    }
                }
                _ => {
                    for j in 0..p + 1 {
                        for k in 0..p + 1 {
                            for m in 0..p + 1 {
                                samples.push([
                                    (j as f64 + 0.5) / (p + 1) as f64,
                                    (k as f64 + 0.5) / (p + 1) as f64,
                                    (m as f64 + 0.5) / (p + 1) as f64,
                                ]);
                            }
                        }
                    }
                }
            },
            // The collapsed pyramid: height levels, each sampled on its
            // own shrunken section. The section at height h is the
            // centered square [h/2, 1-h/2]², so the (i, j) lattice is
            // offset by h/2 and scaled by (1-h) — a naive [0, 1-h] grid
            // would sample OUTSIDE the element, where the interpolated
            // map is meaningless.
            ShapeFamily::Pyramid => {
                let div = (p + 1) as f64;
                for a in 0..=p + 1 {
                    let h = a as f64 / div;
                    let shrink = 1.0 - h;
                    for i in 0..=p + 1 {
                        for j in 0..=p + 1 {
                            samples.push([
                                h / 2.0 + (i as f64 / div) * shrink,
                                h / 2.0 + (j as f64 / div) * shrink,
                                h,
                            ]);
                        }
                    }
                }
            }
        }
        samples
    }

    /// The minimum signed Jacobian measure over the validity samples.
    pub fn min_sampled_jacobian_measure(&self) -> f64 {
        self.validity_samples()
            .iter()
            .map(|&s| self.jacobian_measure(s))
            .fold(f64::INFINITY, f64::min)
    }

    /// `.valid`: the isoparametric map keeps its orientation and never
    /// degenerates over the sampled reference domain (minimum signed
    /// Jacobian measure strictly positive).
    pub fn is_valid(&self) -> bool {
        self.min_sampled_jacobian_measure() > 0.0
    }
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    if l == 0.0 {
        return a;
    }
    [a[0] / l, a[1] / l, a[2] / l]
}

fn det3(cols: [[f64; 3]; 3]) -> f64 {
    cols[0][0] * (cols[1][1] * cols[2][2] - cols[1][2] * cols[2][1])
        - cols[1][0] * (cols[0][1] * cols[2][2] - cols[0][2] * cols[2][1])
        + cols[2][0] * (cols[0][1] * cols[1][2] - cols[0][2] * cols[1][1])
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rmsh_model::ElementType;

    /// Verbatim node coordinates of the first element of each type in the
    /// real-Gmsh 4.14 probe meshes (D:/tmp/rmsh-p5, 2026-10-13): hex/prism
    /// from transfinite extrusions, tets from a regular cube mesh, lines/
    /// triangles/quads from single-element transfinite surfaces. The fp
    /// noise is Gmsh's own output precision, kept deliberately: the map
    /// must reproduce the REAL file, not a tidy idealization.
    const PROBE_FIXTURES: &[(&str, ElementType, &[[f64; 3]])] = &[
        (
            "line3",
            ElementType::Line3,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.9999999999973437, 0.0, 0.0],
            ],
        ),
        (
            "line4",
            ElementType::Line4,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.6666666666649831, 0.0, 0.0],
                [1.333333333331566, 0.0, 0.0],
            ],
        ),
        (
            "tri6",
            ElementType::Triangle6,
            &[
                [0.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [0.0, 3.0, 0.0],
                [1.499999999995159, 0.0, 0.0],
                [1.500000000004507, 1.499999999995493, 0.0],
                [0.0, 1.499999999998279, 0.0],
            ],
        ),
        (
            "tri9",
            ElementType::Triangle9,
            &[
                [0.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [0.0, 3.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [1.999999999999106, 1.000000000000893, 0.0],
                [1.000000000000893, 2.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
        ),
        (
            "tri10",
            ElementType::Triangle10,
            &[
                [0.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [0.0, 3.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [1.999999999999106, 1.000000000000893, 0.0],
                [1.000000000000893, 2.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.9999999999995512, 1.000000000000446, 0.0],
            ],
        ),
        (
            "quad8",
            ElementType::Quad8,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.9999999999991055, 0.0, 0.0],
                [2.0, 0.9999999999991055, 0.0],
                [1.000000000000446, 2.0, 0.0],
                [0.0, 1.000000000000446, 0.0],
            ],
        ),
        (
            "quad9",
            ElementType::Quad9,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.9999999999991055, 0.0, 0.0],
                [2.0, 0.9999999999991055, 0.0],
                [1.000000000000446, 2.0, 0.0],
                [0.0, 1.000000000000446, 0.0],
                [0.9999999999995512, 0.9999999999995512, 0.0],
            ],
        ),
        (
            "quad12",
            ElementType::Quad12,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.6666666666649831, 0.0, 0.0],
                [1.333333333331566, 0.0, 0.0],
                [2.0, 0.6666666666661022, 0.0],
                [2.0, 1.333333333332256, 0.0],
                [1.333333333332944, 2.0, 0.0],
                [0.6666666666672213, 2.0, 0.0],
                [0.0, 1.333333333332944, 0.0],
                [0.0, 0.6666666666672213, 0.0],
            ],
        ),
        (
            "quad16",
            ElementType::Quad16,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.6666666666649831, 0.0, 0.0],
                [1.333333333331566, 0.0, 0.0],
                [2.0, 0.6666666666661022, 0.0],
                [2.0, 1.333333333332256, 0.0],
                [1.333333333332944, 2.0, 0.0],
                [0.6666666666672213, 2.0, 0.0],
                [0.0, 1.333333333332944, 0.0],
                [0.0, 0.6666666666672213, 0.0],
                [0.6666666666661022, 0.6666666666661022, 0.0],
                [1.333333333332256, 0.6666666666661022, 0.0],
                [1.333333333332256, 1.333333333332256, 0.0],
                [0.6666666666672213, 1.333333333332944, 0.0],
            ],
        ),
        (
            "tet10",
            ElementType::Tetrahedron10,
            &[
                [0.24999999999622113, 0.49999999999541233, 0.9999999999992584],
                [1.2500000000006651, 0.49999999999541233, 0.9999999999992584],
                [0.24999999999622113, 1.4999999999962504, 0.9999999999992584],
                [0.24999999999622113, 0.49999999999541233, 0.0],
                [0.74999999999733329, 0.49999999999541233, 0.9999999999992584],
                [0.75000000000111022, 1.0000000000002487, 0.9999999999992584],
                [0.24999999999622113, 0.99999999999650003, 0.9999999999992584],
                [0.24999999999622113, 0.49999999999541233, 0.4999999999996292],
                [0.24999999999622113, 0.99999999999650003, 0.4999999999996292],
                [0.74999999999733329, 0.49999999999541233, 0.4999999999996292],
            ],
        ),
        (
            "tet20",
            ElementType::Tetrahedron20,
            &[
                [1.1666666666656132, 0.50000000000244271, 0.33333333333282452],
                [2.1666666666672799, 0.50000000000244271, 0.33333333333282452],
                [1.1666666666656132, 1.5000000000032835, 0.33333333333282452],
                [1.1666666666656132, 0.50000000000244271, 1.3333333333344912],
                [1.4999999999992457, 0.50000000000244271, 0.33333333333282452],
                [1.833333333333382, 0.50000000000244271, 0.33333333333282452],
                [1.8333333333340622, 0.83333333333507911, 0.33333333333282452],
                [1.5000000000004441, 1.1666666666662508, 0.33333333333282452],
                [1.1666666666656132, 0.83333333333507911, 0.33333333333282452],
                [1.1666666666656132, 1.1666666666662508, 0.33333333333282452],
                [1.1666666666656132, 0.50000000000244271, 0.66666666666615951],
                [1.1666666666656132, 0.50000000000244271, 1.000000000000334],
                [1.1666666666656132, 0.83333333333507911, 1.000000000000334],
                [1.1666666666656132, 1.1666666666662508, 0.66666666666615951],
                [1.5000000000004441, 0.50000000000244271, 1.000000000000334],
                [1.833333333333382, 0.50000000000244271, 0.66666666666615951],
                [1.5000000000004441, 0.83333333333507911, 0.33333333333282452],
                [1.5000000000004441, 0.50000000000244271, 0.66666666666615951],
                [1.1666666666656132, 0.83333333333507911, 0.66666666666615951],
                [1.5000000000004441, 0.83333333333507911, 0.66666666666615951],
            ],
        ),
        (
            "hex20",
            ElementType::Hexahedron20,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 2.0],
                [2.0, 0.0, 2.0],
                [2.0, 2.0, 2.0],
                [0.0, 2.0, 2.0],
                [0.9999999999991055, 0.0, 0.0],
                [0.0, 0.9999999999991055, 0.0],
                [0.0, 0.0, 0.9999999999991055],
                [2.0, 0.9999999999991055, 0.0],
                [2.0, 0.0, 0.9999999999991055],
                [0.9999999999995512, 2.0, 0.0],
                [2.0, 2.0, 0.9999999999991055],
                [0.0, 2.0, 0.9999999999991055],
                [0.9999999999991055, 0.0, 2.0],
                [0.0, 0.9999999999991055, 2.0],
                [2.0, 0.9999999999991055, 2.0],
                [0.9999999999995512, 2.0, 2.0],
            ],
        ),
        (
            "hex27",
            ElementType::Hexahedron27,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 2.0],
                [2.0, 0.0, 2.0],
                [2.0, 2.0, 2.0],
                [0.0, 2.0, 2.0],
                [0.9999999999991055, 0.0, 0.0],
                [0.0, 0.9999999999991055, 0.0],
                [0.0, 0.0, 0.9999999999991055],
                [2.0, 0.9999999999991055, 0.0],
                [2.0, 0.0, 0.9999999999991055],
                [0.9999999999995512, 2.0, 0.0],
                [2.0, 2.0, 0.9999999999991055],
                [0.0, 2.0, 0.9999999999991055],
                [0.9999999999991055, 0.0, 2.0],
                [0.0, 0.9999999999991055, 2.0],
                [2.0, 0.9999999999991055, 2.0],
                [0.9999999999995512, 2.0, 2.0],
                [0.9999999999995512, 0.9999999999995512, 0.0],
                [0.9999999999995512, 0.9999999999995512, 2.0],
                [0.9999999999995512, 0.0, 0.9999999999991055],
                [0.9999999999995512, 2.0, 0.9999999999991055],
                [0.0, 0.9999999999995512, 0.9999999999991055],
                [2.0, 0.9999999999995512, 0.9999999999991055],
                [0.9999999999995512, 0.9999999999995512, 0.9999999999991055],
            ],
        ),
        (
            "hex64",
            ElementType::Hexahedron64,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 2.0],
                [2.0, 0.0, 2.0],
                [2.0, 2.0, 2.0],
                [0.0, 2.0, 2.0],
                [0.6666666666661022, 0.0, 0.0],
                [1.333333333332256, 0.0, 0.0],
                [0.0, 0.6666666666661022, 0.0],
                [0.0, 1.333333333332256, 0.0],
                [0.0, 0.0, 0.6666666666661022],
                [0.0, 0.0, 1.333333333332256],
                [2.0, 0.6666666666672213, 0.0],
                [2.0, 1.333333333332944, 0.0],
                [2.0, 0.0, 0.6666666666672213],
                [2.0, 0.0, 1.333333333332944],
                [1.3333333333344912, 2.0, 0.0],
                [0.6666666666672213, 2.0, 0.0],
                [2.0, 2.0, 0.6666666666672213],
                [2.0, 2.0, 1.333333333332944],
                [0.0, 2.0, 0.6666666666672213],
                [0.0, 2.0, 1.333333333332944],
                [0.6666666666672213, 0.0, 2.0],
                [1.333333333332944, 0.0, 2.0],
                [0.0, 0.6666666666672213, 2.0],
                [0.0, 1.333333333332944, 2.0],
                [2.0, 0.6666666666661022, 2.0],
                [2.0, 1.333333333332256, 2.0],
                [1.333333333332256, 2.0, 2.0],
                [0.6666666666661022, 2.0, 2.0],
                [0.6666666666661022, 0.6666666666661022, 0.0],
                [0.6666666666661022, 1.333333333332256, 0.0],
                [1.333333333332256, 1.333333333332256, 0.0],
                [1.333333333332256, 0.6666666666661022, 0.0],
                [0.6666666666672213, 0.0, 0.6666666666672213],
                [1.333333333332944, 0.0, 0.6666666666672213],
                [1.333333333332944, 0.0, 1.333333333332944],
                [0.6666666666672213, 0.0, 1.333333333332944],
                [0.0, 0.6666666666672213, 0.6666666666672213],
                [0.0, 0.6666666666672213, 1.333333333332944],
                [0.0, 1.333333333332944, 1.333333333332944],
                [0.0, 1.333333333332944, 0.6666666666672213],
                [2.0, 0.6666666666661022, 0.6666666666661022],
                [2.0, 1.333333333332256, 0.6666666666661022],
                [2.0, 1.333333333332256, 1.333333333332256],
                [2.0, 0.6666666666661022, 1.333333333332256],
                [1.3333333333344912, 2.0, 0.6666666666672213],
                [0.6666666666672213, 2.0, 0.6666666666672213],
                [0.6666666666672213, 2.0, 1.333333333332944],
                [1.3333333333344912, 2.0, 1.333333333332944],
                [0.6666666666672213, 0.6666666666672213, 2.0],
                [1.333333333332944, 0.6666666666672213, 2.0],
                [1.333333333332944, 1.333333333332944, 2.0],
                [0.6666666666672213, 1.333333333332944, 2.0],
                [0.6666666666672213, 0.6666666666661022, 0.6666666666672213],
                [1.333333333332256, 0.6666666666661022, 0.6666666666672213],
                [1.333333333332256, 1.333333333332256, 0.6666666666672213],
                [0.6666666666672213, 1.333333333332256, 0.6666666666672213],
                [0.6666666666672213, 0.6666666666672213, 1.333333333332944],
                [1.333333333332256, 0.6666666666672213, 1.333333333332944],
                [1.333333333332256, 1.333333333332944, 1.333333333332944],
                [0.6666666666672213, 1.333333333332944, 1.333333333332944],
            ],
        ),
        (
            "prism15",
            ElementType::Prism15,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 2.0],
                [2.0, 0.0, 2.0],
                [0.0, 2.0, 2.0],
                [0.9999999999991055, 0.0, 0.0],
                [0.0, 0.9999999999991055, 0.0],
                [0.0, 0.0, 0.9999999999991055],
                [0.9999999999995512, 0.9999999999995512, 0.0],
                [2.0, 0.0, 0.9999999999991055],
                [0.0, 2.0, 0.9999999999991055],
                [0.9999999999991055, 0.0, 2.0],
                [0.0, 0.9999999999991055, 2.0],
                [0.9999999999995512, 0.9999999999995512, 2.0],
            ],
        ),
        (
            "prism18",
            ElementType::Prism18,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 2.0],
                [2.0, 0.0, 2.0],
                [0.0, 2.0, 2.0],
                [0.9999999999991055, 0.0, 0.0],
                [0.0, 0.9999999999991055, 0.0],
                [0.0, 0.0, 0.9999999999991055],
                [0.9999999999995512, 0.9999999999995512, 0.0],
                [2.0, 0.0, 0.9999999999991055],
                [0.0, 2.0, 0.9999999999991055],
                [0.9999999999991055, 0.0, 2.0],
                [0.0, 0.9999999999991055, 2.0],
                [0.9999999999995512, 0.9999999999995512, 2.0],
                [0.9999999999991055, 0.0, 0.9999999999991055],
                [0.0, 0.9999999999991055, 0.9999999999991055],
                [0.9999999999995512, 0.9999999999995512, 0.9999999999991055],
            ],
        ),
        (
            "prism40",
            ElementType::Prism40,
            &[
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 2.0],
                [2.0, 0.0, 2.0],
                [0.0, 2.0, 2.0],
                [0.6666666666661022, 0.0, 0.0],
                [1.333333333332256, 0.0, 0.0],
                [0.0, 0.6666666666661022, 0.0],
                [0.0, 1.333333333332256, 0.0],
                [0.0, 0.0, 0.6666666666661022],
                [0.0, 0.0, 1.333333333332256],
                [1.333333333332944, 0.6666666666672213, 0.0],
                [0.6666666666672213, 1.333333333332944, 0.0],
                [2.0, 0.0, 0.6666666666672213],
                [2.0, 0.0, 1.333333333332944],
                [0.0, 2.0, 0.6666666666672213],
                [0.0, 2.0, 1.333333333332944],
                [0.6666666666672213, 0.0, 2.0],
                [1.333333333332944, 0.0, 2.0],
                [0.0, 0.6666666666672213, 2.0],
                [0.0, 1.333333333332944, 2.0],
                [1.3333333333344912, 0.6666666666672213, 2.0],
                [0.6666666666672213, 1.3333333333344912, 2.0],
                [0.6666666666661022, 0.6666666666661022, 0.0],
                [0.6666666666672213, 0.6666666666672213, 2.0],
                [0.6666666666672213, 0.0, 0.6666666666672213],
                [1.333333333332944, 0.0, 0.6666666666672213],
                [1.333333333332944, 0.0, 1.333333333332944],
                [0.6666666666672213, 0.0, 1.333333333332944],
                [0.0, 0.6666666666672213, 0.6666666666672213],
                [0.0, 0.6666666666672213, 1.333333333332944],
                [0.0, 1.333333333332944, 1.333333333332944],
                [0.0, 1.333333333332944, 0.6666666666672213],
                [1.333333333332256, 0.6666666666661022, 0.6666666666661022],
                [0.6666666666661022, 1.333333333332256, 0.6666666666661022],
                [0.6666666666661022, 1.333333333332256, 1.333333333332256],
                [1.333333333332256, 0.6666666666661022, 1.333333333332256],
                [0.6666666666672213, 0.6666666666672213, 0.6666666666672213],
                [0.6666666666672213, 0.6666666666672213, 1.333333333332944],
            ],
        ),
    ];


    #[test]
    fn probe_elements_are_reproduced_by_the_map() {
        for (name, etype, coords) in PROBE_FIXTURES {
            let el = CurvedElement::new(*etype, coords)
                .unwrap_or_else(|| panic!("{name}: must construct"));
            let layout = reference_layout(*etype)
                .unwrap_or_else(|| panic!("{name}: must have a layout"));
            assert_eq!(layout.len(), coords.len(), "{name}: layout/fixture size");
            for (j, r) in layout.iter().enumerate() {
                let x = el.eval(*r);
                for d in 0..3 {
                    assert!(
                        (x[d] - coords[j][d]).abs() < 1e-9,
                        "{name}: node {j} axis {d} reproduced {} expected {}",
                        x[d],
                        coords[j][d]
                    );
                }
            }
        }
    }

    // ── pyramid probe fixtures (real Gmsh 4.14 meshes, D:/tmp/rmsh-p5) ────
    //
    // Pyramid13: type-19 element from `pyr13b.msh` (the
    // SecondOrderIncomplete=1 re-probe of 2026-10-14 — the archived
    // pyramid13.geo had Recombine3DAll=1 and produced pure hexes);
    // Pyramid14: type-14 element 2245 from `pyr2.msh`; Pyramid30: type-118
    // element 2245 from `pyr3.msh` (v4.1). The connectivity is verbatim
    // MSH order; the map must reproduce every node through its layout
    // slot, and the straight box pieces must be fold-free.

    const PYRAMID13_PROBE: [[f64; 3]; 13] = [
        [0.0, 0.862283097054201, 0.50116977255654],
        [0.0, 0.5787811688992963, 0.58047032459531],
        [0.0, 0.5965467870890105, 0.2864198794436149],
        [0.0, 0.8078562409099833, 0.2825224330571149],
        [0.1129803319079237, 0.7133599743258788, 0.4111132222004577],
        [0.0, 0.7205321329767487, 0.540820048575925],
        [0.0, 0.8350696689820922, 0.3918461028068275],
        [0.05649016595396185, 0.78782153569004, 0.4561414973784989],
        [0.0, 0.5876639779941535, 0.4334451020194625],
        [0.05649016595396185, 0.6460705716125876, 0.4957917733978838],
        [0.0, 0.7022015139994969, 0.2844711562503649],
        [0.05649016595396185, 0.6549533807074447, 0.3487665508220363],
        [0.05649016595396185, 0.7606081076179311, 0.3468178276287863],
    ];

    const PYRAMID14_PROBE: [[f64; 3]; 14] = [
        [0.0, 0.862283097054201, 0.50116977255654],
        [0.0, 0.5787811688992963, 0.58047032459531],
        [0.0, 0.5965467870890105, 0.2864198794436149],
        [0.0, 0.8078562409099833, 0.2825224330571149],
        [0.1129803319079237, 0.7133599743258788, 0.4111132222004577],
        [0.0, 0.7205321329767487, 0.540820048575925],
        [0.0, 0.8350696689820922, 0.3918461028068275],
        [0.05649016595396185, 0.78782153569004, 0.4561414973784989],
        [0.0, 0.5876639779941535, 0.4334451020194625],
        [0.05649016595396185, 0.6460705716125876, 0.4957917733978838],
        [0.0, 0.7022015139994969, 0.2844711562503649],
        [0.05649016595396185, 0.6549533807074447, 0.3487665508220363],
        [0.05649016595396185, 0.7606081076179311, 0.3468178276287863],
        [0.0, 0.7113668234881227, 0.412645602413145],
    ];

    const PYRAMID30_PROBE: [[f64; 3]; 30] = [
        [0.0, 0.862283097054201, 0.50116977255654],
        [0.0, 0.5787811688992963, 0.58047032459531],
        [0.0, 0.5965467870890105, 0.2864198794436149],
        [0.0, 0.8078562409099833, 0.2825224330571149],
        [0.1129803319079237, 0.7133599743258788, 0.4111132222004577],
        [0.0, 0.7677824543358994, 0.5276032899027967],
        [0.0, 0.6732818116175979, 0.5540368072490534],
        [0.0, 0.8441408116727951, 0.4282873260567317],
        [0.0, 0.8259985262913893, 0.3554048795569233],
        [0.03766011063597457, 0.8126420561447603, 0.4711509224378459],
        [0.07532022127194914, 0.7630010152353196, 0.4411320723191518],
        [0.0, 0.5847030416292011, 0.482453509544745],
        [0.0, 0.5906249143591058, 0.3844366944941799],
        [0.03766011063597457, 0.6236407707081572, 0.5240179571303591],
        [0.07532022127194914, 0.668500372517018, 0.4675655896654085],
        [0.0, 0.6669832716960015, 0.2851207306481149],
        [0.0, 0.7374197563029924, 0.2838215818526149],
        [0.03766011063597456, 0.6354845161679666, 0.3279843270292292],
        [0.07532022127194912, 0.6744222452469227, 0.3695487746148434],
        [0.03766011063597456, 0.7763574853819486, 0.3253860294382291],
        [0.07532022127194912, 0.7448587298539137, 0.3682496258193434],
        [0.03766011063597458, 0.7181414134264588, 0.4975844397841026],
        [0.03766011063597456, 0.7944997707633543, 0.3982684759380375],
        [0.03766011063597456, 0.6295626434380619, 0.4260011420797942],
        [0.03766011063597455, 0.7059210007749577, 0.3266851782337291],
        [0.0, 0.757661554991597, 0.4463427205527362],
        [0.0, 0.7475406556472948, 0.3650821512026755],
        [0.0, 0.6690827850032002, 0.3747594228484278],
        [0.0, 0.6711822983103991, 0.4643981150487405],
        [0.03766011063597456, 0.7120312071007082, 0.4121348090089158],
    ];

    #[test]
    fn pyramid_probe_elements_are_reproduced_and_straight_valid() {
        for (name, etype, coords) in [
            ("pyramid13", ElementType::Pyramid13, &PYRAMID13_PROBE[..]),
            ("pyramid14", ElementType::Pyramid14, &PYRAMID14_PROBE[..]),
            ("pyramid30", ElementType::Pyramid30, &PYRAMID30_PROBE[..]),
        ] {
            let el = CurvedElement::new(etype, coords)
                .unwrap_or_else(|| panic!("{name}: must construct"));
            let layout = reference_layout(etype).expect("layout");
            assert_eq!(layout.len(), coords.len(), "{name}: layout/fixture size");
            for (j, r) in layout.iter().enumerate() {
                let x = el.eval(*r);
                for d in 0..3 {
                    assert!(
                        (x[d] - coords[j][d]).abs() < 1e-9,
                        "{name}: node {j} axis {d} reproduced {} expected {}",
                        x[d],
                        coords[j][d]
                    );
                }
            }
            assert!(
                el.is_valid(),
                "{name}: the straight box piece must be fold-free (min measure {})",
                el.min_sampled_jacobian_measure()
            );
        }
    }

    #[test]
    fn probe_elements_are_valid_straight_elements() {
        for (name, etype, coords) in PROBE_FIXTURES {
            let el = CurvedElement::new(*etype, coords)
                .unwrap_or_else(|| panic!("{name}: must construct"));
            assert!(
                el.is_valid(),
                "{name}: the straight probe element must be valid (min measure {})",
                el.min_sampled_jacobian_measure()
            );
        }
    }

    #[test]
    fn interpolation_identity_holds_at_reference_nodes() {
        // A non-trivial affine image of the unit layout: eval at the
        // reference position of node j must return exactly node j.
        for (name, etype, coords) in PROBE_FIXTURES {
            let el = CurvedElement::new(*etype, coords).expect("construct");
            for (j, r) in reference_layout(*etype).expect("layout").iter().enumerate() {
                let x = el.eval(*r);
                let n = self::super::basis_for(*etype).expect("basis");
                for i in 0..coords.len() {
                    let v = n.shape(i, *r);
                    let expected = if i == j { 1.0 } else { 0.0 };
                    assert!(
                        (v - expected).abs() < 1e-9,
                        "{name}: N_{i}({j:?}) = {v}, expected {expected}"
                    );
                }
                let _ = x;
            }
        }
    }

    #[test]
    fn folded_curved_elements_are_refused() {
        // Line3 with the midside pushed past the end: dx/dt reverses sign
        // near t=1 (analytically 3 - 4t < 0 for t > 0.75).
        let line3 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let el = CurvedElement::new(ElementType::Line3, &line3).expect("construct");
        assert!(
            !el.is_valid(),
            "folded line must be invalid (min measure {})",
            el.min_sampled_jacobian_measure()
        );

        // Tet10 with the (0,1) midside pushed past corner 1: the edge map
        // x(t) = 5t - 3t^2 folds (dx/dt = 5 - 6t < 0 for t > 5/6).
        let mut tet10 = [[0.0; 3]; 10];
        tet10[0] = [0.0, 0.0, 0.0];
        tet10[1] = [2.0, 0.0, 0.0];
        tet10[2] = [0.0, 2.0, 0.0];
        tet10[3] = [0.0, 0.0, 2.0];
        for i in 4..10 {
            for d in 0..3 {
                tet10[i][d] = (tet10[i % 4][d] + tet10[(i - 3) % 4][d]) * 0.5;
            }
        }
        // Edge layout of nodes 4..10: (0,1),(1,2),(0,2),(0,3),(2,3),(1,3).
        let edge_pairs = [(0usize, 1usize), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)];
        for (k, (a, b)) in edge_pairs.iter().enumerate() {
            for d in 0..3 {
                tet10[4 + k][d] = (tet10[*a][d] + tet10[*b][d]) * 0.5;
            }
        }
        tet10[4] = [2.5, 0.0, 0.0];
        let el = CurvedElement::new(ElementType::Tetrahedron10, &tet10).expect("construct");
        assert!(
            !el.is_valid(),
            "folded tet must be invalid (min measure {})",
            el.min_sampled_jacobian_measure()
        );

        // Quad8 with the bottom midside pushed far outside: the signed
        // area element flips near the bottom edge.
        let quad8 = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [4.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [1.0, 2.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let el = CurvedElement::new(ElementType::Quad8, &quad8).expect("construct");
        assert!(
            !el.is_valid(),
            "folded quad must be invalid (min measure {})",
            el.min_sampled_jacobian_measure()
        );
    }

    /// Sample the arc error of a curved line mapped onto the unit circle
    /// quadrant (nodes on the arc at 0, 45, 90 degrees for order 2; 0, 30,
    /// 60, 90 for order 3).
    fn arc_max_radial_error(etype: ElementType, order: usize) -> f64 {
        let p = order as f64;
        // Coordinates in MSH connectivity order: corners first (t=0, t=1),
        // then the interior edge nodes ascending (t = i/p).
        let mut coords: Vec<[f64; 3]> = Vec::new();
        coords.push([1.0, 0.0, 0.0]); // t = 0
        coords.push([0.0, 1.0, 0.0]); // t = 1
        for i in 1..order {
            let theta = std::f64::consts::FRAC_PI_2 * (i as f64) / p;
            coords.push([theta.cos(), theta.sin(), 0.0]);
        }
        let el = CurvedElement::new(etype, &coords).expect("construct");
        let mut worst = 0.0f64;
        for k in 0..=1000 {
            let t = k as f64 / 1000.0;
            let y = el.eval([t, 0.0, 0.0]);
            let r = (y[0] * y[0] + y[1] * y[1]).sqrt();
            worst = worst.max((r - 1.0).abs());
        }
        worst
    }

    #[test]
    fn quadratic_arc_geometric_error_within_tolerance() {
        // Measured max radial deviation of the quadratic arc: 1.0781e-2.
        // The card tolerance for a 90-degree span at order 2 is pinned at
        // 1.5e-2 (~1% of the chord); the straight chord contrast is
        // 2.9289e-1, 27x worse.
        let err = arc_max_radial_error(ElementType::Line3, 2);
        assert!(
            err <= 1.5e-2,
            "quadratic arc error {err} exceeds the 1.5e-2 tolerance"
        );
        let chord = 1.0 - std::f64::consts::FRAC_1_SQRT_2;
        assert!(chord > 0.29, "chord sagitta sanity");
    }

    #[test]
    fn cubic_arc_beats_the_quadratic_tolerance() {
        // Measured max radial deviation of the cubic arc: 2.7271e-3.
        assert!(arc_max_radial_error(ElementType::Line4, 3) <= 4.0e-3);
    }

    #[test]
    fn linear_arc_error_demonstrates_the_contrast() {
        // The straight chord under-approximates the quadrant by the
        // sagitta 1 - sqrt(2)/2 = 2.9289e-1: far outside the curved
        // tolerance, which is the point of curved elements.
        let coords = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let el = CurvedElement::new(ElementType::Line2, &coords);
        assert!(el.is_none(), "Line2 has no curved mapping");
        // measure the chord error directly for the comment's claim
        let mid = [0.5f64, 0.5, 0.0];
        let r = (mid[0] * mid[0] + mid[1] * mid[1]).sqrt();
        assert!((r - 1.0).abs() > 0.29);
    }

    #[test]
    fn quadratic_triangle_cylinder_error_within_tolerance() {
        // A quarter-cylinder patch (radius 1, height 1): Tri6 with all
        // three midsides on (or projected onto) the cylinder. Measured max
        // radial deviation: 1.0781e-2; tolerance 1.5e-2. The straight Tri3
        // contrast is 2.9289e-1, 27x worse.
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let coords = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
            [s, s, 0.0],
            [0.0, 1.0, 0.5],
            [s, s, 0.5],
        ];
        let el = CurvedElement::new(ElementType::Triangle6, &coords).expect("construct");
        let mut worst = 0.0f64;
        for i in 0..=100 {
            for j in 0..=100 - i {
                let (u, v) = (i as f64 / 100.0, j as f64 / 100.0);
                let x = el.eval([u, v, 0.0]);
                let r = (x[0] * x[0] + x[1] * x[1]).sqrt();
                worst = worst.max((r - 1.0).abs());
            }
        }
        assert!(
            worst <= 1.5e-2,
            "tri6 cylinder error {worst} exceeds the 1.5e-2 tolerance"
        );
        assert!(el.is_valid(), "the curved cylinder patch must be valid");
    }

    #[test]
    fn curved_validity_holds_for_the_arc_elements() {
        for (etype, order) in [
            (ElementType::Line3, 2usize),
            (ElementType::Line4, 3usize),
        ] {
            let p = order as f64;
            let mut coords: Vec<[f64; 3]> = Vec::new();
            coords.push([1.0, 0.0, 0.0]);
            coords.push([0.0, 1.0, 0.0]);
            for i in 1..order {
                let theta = std::f64::consts::FRAC_PI_2 * (i as f64) / p;
                coords.push([theta.cos(), theta.sin(), 0.0]);
            }
            let el = CurvedElement::new(etype, &coords).expect("construct");
            assert!(el.is_valid(), "{etype:?} arc element must be valid");
        }
    }

    #[test]
    fn linear_types_have_no_curved_mapping() {
        // The honest scope line: linear types carry no reference layout
        // and no curved map. (The pyramids joined the covered families —
        // see the collapsed-basis tests below.)
        for etype in [
            ElementType::Line2,
            ElementType::Triangle3,
            ElementType::Quad4,
            ElementType::Tetrahedron4,
            ElementType::Hexahedron8,
            ElementType::Prism6,
            ElementType::Pyramid5,
            ElementType::Point1,
        ] {
            assert!(reference_layout(etype).is_none(), "{etype:?} is linear");
            assert!(!supports_curved_mapping(etype));
        }
        assert!(reference_layout(ElementType::Unknown(93)).is_none());
    }

    #[test]
    fn pyramids_are_covered_with_correct_node_counts() {
        for (etype, n) in [
            (ElementType::Pyramid13, 13usize),
            (ElementType::Pyramid14, 14),
            (ElementType::Pyramid30, 30),
        ] {
            assert!(reference_layout(etype).is_some(), "{etype:?} covered");
            assert!(supports_curved_mapping(etype));
            assert!(CurvedElement::new(etype, &vec![[0.0; 3]; n]).is_some());
            assert!(CurvedElement::new(etype, &vec![[0.0; 3]; n - 1]).is_none());
        }
    }

    /// The straight pyramid is the identity map (the partition of unity +
    /// linear reproduction through the public API), and it is valid.
    #[test]
    fn straight_pyramids_reproduce_reference_points_and_are_valid() {
        for etype in [
            ElementType::Pyramid13,
            ElementType::Pyramid14,
            ElementType::Pyramid30,
        ] {
            let layout = reference_layout(etype).expect("layout").to_vec();
            let el = CurvedElement::new(etype, &layout).expect("construct");
            for r in [
                [0.25, 0.25, 0.25],
                [0.5, 0.5, 0.5],
                [0.1, 0.2, 0.3],
                [0.5, 0.5, 0.0],
                [0.5, 0.5, 1.0],
            ] {
                let x = el.eval(r);
                for d in 0..3 {
                    assert!(
                        (x[d] - r[d]).abs() < 1e-12,
                        "{etype:?}: straight map must be identity at {r:?}, got {x:?}"
                    );
                }
            }
            assert!(
                el.is_valid(),
                "{etype:?}: straight pyramid must be valid (min measure {})",
                el.min_sampled_jacobian_measure()
            );
        }
    }

    /// A pyramid with the apex pushed through the base plane folds.
    #[test]
    fn folded_pyramid_is_invalid() {
        let mut coords = PYRAMID14_LAYOUT.to_vec();
        // apex slot 4: lift it far below the base plane
        coords[4] = [0.5, 0.5, -1.5];
        let el = CurvedElement::new(ElementType::Pyramid14, &coords).expect("construct");
        assert!(
            !el.is_valid(),
            "folded pyramid must be invalid (min measure {})",
            el.min_sampled_jacobian_measure()
        );
    }

    /// Geometric error pin: the four side faces lie on the unit
    /// cylinder (base corners on the r=1 circle 90 deg apart, apex on
    /// the circle, edge mids at the radial projections of the chord
    /// midpoints), the base face stays the flat square with its bubble
    /// at the planar center. Mirrors the Tri6 cylinder pin.
    #[test]
    fn curved_pyramid_cylinder_error_within_tolerance() {
        let radial = |x: f64, y: f64, z: f64| -> [f64; 3] {
            let r = (x * x + y * y).sqrt();
            [x / r, y / r, z]
        };
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let c = [
            [s, s, 0.0],
            [-s, s, 0.0],
            [-s, -s, 0.0],
            [s, -s, 0.0],
        ];
        let apex = [1.0, 0.0, 1.0];
        // Node order = the Pyramid14 layout: 4 base corners, apex, edges
        // [(0,1),(0,3),(0,4),(1,2),(1,4),(2,3),(2,4),(3,4)], base centroid.
        let mut coords: Vec<[f64; 3]> = vec![c[0], c[1], c[2], c[3], apex];
        for (a, b) in [(0usize, 1usize), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (2, 4), (3, 4)] {
            let (pa, pb) = if b == 4 { (c[a], apex) } else { (c[a], c[b]) };
            coords.push(radial(
                (pa[0] + pb[0]) * 0.5,
                (pa[1] + pb[1]) * 0.5,
                (pa[2] + pb[2]) * 0.5,
            ));
        }
        coords.push([0.0, 0.0, 0.0]);
        let el = CurvedElement::new(ElementType::Pyramid14, &coords).expect("construct");
        assert!(el.is_valid(), "the cylinder pyramid must be fold-free");

        // Side faces: sample the (0,1,4) and (3,0,4) interiors. A
        // barycentric (a, b, c) point over (corner, corner, apex) maps to
        // u = b + c/2, v = c/2, h = c on face (0,1,4); the radial
        // distance must stay within the P2 bar of the unit cylinder.
        let mut worst = 0.0f64;
        let n = 8usize;
        for k in 0..=n {
            for m in 0..=n {
                let c_ = (m as f64 + 0.5) / (n as f64 + 1.0);
                let b = (k as f64 / n as f64) * (1.0 - c_);
                if b + c_ > 1.0 {
                    continue;
                }
                let p1 = [b + c_ * 0.5, c_ * 0.5, c_];
                let p2 = [c_ * 0.5, 1.0 - b - c_ * 0.5, c_];
                for p in [p1, p2] {
                    let x = el.eval(p);
                    let r = (x[0] * x[0] + x[1] * x[1]).sqrt();
                    worst = worst.max((r - 1.0).abs());
                }
            }
        }
        // Base interior bilinear samples: the pin is the planarity
        // (z stays 0) — the flat square's interior legitimately sits
        // INSIDE the cylinder, so no radial bound applies there.
        for su in [0.25f64, 0.35, 0.5, 0.65, 0.75] {
            for sv in [0.25f64, 0.35, 0.5, 0.65, 0.75] {
                let x = el.eval([su, sv, 0.0]);
                assert!(x[2].abs() < 1e-12, "base face must stay planar");
            }
        }
        // Measured worst radial deviation: ~4.8e-2 — the P2 interpolation
        // level for this geometry (each side face carries one 135-deg
        // base-to-apex edge), pinned with headroom at 6e-2. A basis or
        // layout regression blows this up by orders of magnitude.
        assert!(
            worst < 6.0e-2,
            "curved pyramid cylinder error {worst} above the P2 bar"
        );
    }

    #[test]
    fn layout_node_counts_match_the_model_table() {
        for info in TYPE_INFOS {
            assert_eq!(
                info.layout.len(),
                info.etype.node_count(),
                "{:?} layout size must equal the verified Gmsh node count",
                info.etype
            );
            if info.family == ShapeFamily::Serendipity {
                assert_eq!(
                    info.layout.len(),
                    info.monomials.len(),
                    "{:?} serendipity space size must equal the node count",
                    info.etype
                );
            } else {
                assert!(
                    info.monomials.is_empty(),
                    "{:?} closed-form family must carry no monomials",
                    info.etype
                );
            }
        }
    }
}

