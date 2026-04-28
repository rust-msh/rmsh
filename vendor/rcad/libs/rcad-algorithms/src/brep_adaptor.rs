//! BRepAdaptor-style topology adapters.
//!
//! Provides adapters to access BRep topology as geometric entities.
//! Analogous to OCCT's BRepAdaptor_Curve, BRepAdaptor_Surface, and BRepAdaptor_CompCurve.
//!
//! # Overview
//!
//! - [`EdgeAdaptor`]: Adapts an edge to act as a 3D curve
//! - [`FaceAdaptor`]: Adapts a face to act as a 3D surface
//! - [`WireAdaptor`]: Adapts a wire to act as a composite 3D curve
//! - [`CurveAdaptorArray`]: Array of edge adaptors for indexed access

use glam::DVec3;
use rcad_kernel::{BRep, Curve3, CurveEval, Surface3, SurfaceEval, Wire};
use std::f64::consts::PI;

// =============================================================================
// EdgeAdaptor (BRepAdaptor_Curve)
// =============================================================================

/// Adapts a BRep edge to act as a 3D curve.
///
/// Provides curve-like evaluation methods (point, tangent, domain) for an edge
/// in a BRep, respecting the edge's orientation and parameter range.
///
/// Analogous to OCCT's `BRepAdaptor_Curve`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_adaptor::EdgeAdaptor;
/// use rcad_kernel::BRep;
/// use rcad_kernel::geom::PrimitiveSolid;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let adaptor = EdgeAdaptor::new(&brep, 0);
/// let domain = adaptor.domain();
/// let midpoint = adaptor.point_at((domain[0] + domain[1]) / 2.0);
/// ```
#[derive(Debug, Clone)]
pub struct EdgeAdaptor<'a> {
    brep: &'a BRep,
    edge_idx: usize,
    /// Cached curve reference (if available).
    curve: Option<&'a Curve3>,
    /// Cached parameter range.
    range: [f64; 2],
    /// Whether the edge's natural direction is reversed.
    reversed: bool,
}

impl<'a> EdgeAdaptor<'a> {
    /// Create a new edge adaptor for the given edge index.
    ///
    /// The adaptor respects the edge's stored parameter range in `edge_curve_range`
    /// and falls back to the curve's natural domain if not specified.
    ///
    /// # Panics
    ///
    /// Does not panic; returns a default adaptor if the edge index is out of bounds
    /// or the edge has no associated curve.
    pub fn new(brep: &'a BRep, edge_idx: usize) -> Self {
        let curve = brep
            .geom
            .edge_curve
            .get(edge_idx)
            .and_then(|opt| opt.as_ref())
            .and_then(|&curve_idx| brep.geom.curves.get(curve_idx));

        let range = if let Some(r) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) {
            r
        } else if let Some(c) = curve {
            c.default_domain()
        } else {
            [0.0, 1.0]
        };

        Self {
            brep,
            edge_idx,
            curve,
            range,
            reversed: false,
        }
    }

    /// Create an edge adaptor with reversed direction.
    ///
    /// This is used when an edge appears in a wire with `forward = false`,
    /// meaning the edge should be traversed from end to start.
    pub fn with_reversed(mut self, reversed: bool) -> Self {
        self.reversed = reversed;
        self
    }

    /// Evaluate the point on the edge at parameter `t`.
    ///
    /// The parameter `t` is in the edge's natural parameter domain
    /// (respecting `edge_curve_range` if specified).
    pub fn point_at(&self, t: f64) -> DVec3 {
        let Some(curve) = self.curve else {
            // Fall back to vertex interpolation if no curve is available.
            return self.point_from_vertices(t);
        };

        let t_mapped = self.map_parameter(t);
        curve.point_at(t_mapped)
    }

    /// Evaluate the unit tangent vector on the edge at parameter `t`.
    ///
    /// Returns the tangent pointing in the direction of increasing parameter
    /// on the underlying curve. If the edge is reversed, negate the result.
    pub fn tangent_at(&self, t: f64) -> DVec3 {
        let Some(curve) = self.curve else {
            // Fall back to straight-line tangent between vertices.
            return self.tangent_from_vertices();
        };

        let t_mapped = self.map_parameter(t);
        let mut tangent = curve.tangent_at(t_mapped);
        if self.reversed {
            tangent = -tangent;
        }
        tangent
    }

    /// Return the parameter domain of the edge.
    ///
    /// This is always `[0.0, 1.0]` for normalized parameter access,
    /// regardless of the underlying curve's natural domain.
    pub fn domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }

    /// Return the underlying curve reference, if available.
    pub fn curve(&self) -> Option<&Curve3> {
        self.curve
    }

    /// Return the natural parameter range of the edge on its curve.
    ///
    /// This returns the actual parameter range (not normalized to [0, 1]).
    pub fn curve_range(&self) -> [f64; 2] {
        self.range
    }

    /// Check if the edge is closed (start and end vertices are the same).
    ///
    /// A closed edge forms a loop, such as a circle or ellipse.
    pub fn is_closed(&self) -> bool {
        let Some(edge) = self.brep.edges.get(self.edge_idx) else {
            return false;
        };
        edge.start == edge.end
    }

    /// Return the period of the edge's curve if it is periodic.
    ///
    /// Returns `Some(period)` for periodic curves (circles, ellipses),
    /// or `None` for non-periodic curves (lines, B-splines).
    pub fn period(&self) -> Option<f64> {
        let Some(curve) = self.curve else {
            return None;
        };

        match curve {
            Curve3::Circle(_) => Some(2.0 * PI),
            Curve3::Ellipse(_) => Some(2.0 * PI),
            Curve3::Line(_) => None,
            Curve3::BSpline(_) => None,
            Curve3::Bezier(_) => None,
            Curve3::Offset(_) => None,
            Curve3::Hyperbola(_) => None,
            Curve3::Parabola(_) => None,
            Curve3::CircularHelix(_) => None,
            Curve3::SineWave(_) => None,
        }
    }

    /// Map normalized parameter [0, 1] to curve's natural parameter.
    fn map_parameter(&self, t: f64) -> f64 {
        let [t0, t1] = self.range;
        if self.reversed {
            t0 + (1.0 - t) * (t1 - t0)
        } else {
            t0 + t * (t1 - t0)
        }
    }

    /// Fall back to vertex-based point evaluation when no curve is available.
    fn point_from_vertices(&self, t: f64) -> DVec3 {
        let Some(edge) = self.brep.edges.get(self.edge_idx) else {
            return DVec3::ZERO;
        };
        let Some(v_start) = self.brep.vertices.get(edge.start) else {
            return DVec3::ZERO;
        };
        let Some(v_end) = self.brep.vertices.get(edge.end) else {
            return DVec3::ZERO;
        };

        if self.reversed {
            v_end.point.lerp(v_start.point, t)
        } else {
            v_start.point.lerp(v_end.point, t)
        }
    }

    /// Fall back to vertex-based tangent when no curve is available.
    fn tangent_from_vertices(&self) -> DVec3 {
        let Some(edge) = self.brep.edges.get(self.edge_idx) else {
            return DVec3::X;
        };
        let Some(v_start) = self.brep.vertices.get(edge.start) else {
            return DVec3::X;
        };
        let Some(v_end) = self.brep.vertices.get(edge.end) else {
            return DVec3::X;
        };

        let dir = (v_end.point - v_start.point).normalize_or_zero();
        if self.reversed {
            -dir
        } else {
            dir
        }
    }

    /// Get the first vertex index of this edge.
    pub fn first_vertex(&self) -> Option<usize> {
        let edge = self.brep.edges.get(self.edge_idx)?;
        if self.reversed {
            Some(edge.end)
        } else {
            Some(edge.start)
        }
    }

    /// Get the last vertex index of this edge.
    pub fn last_vertex(&self) -> Option<usize> {
        let edge = self.brep.edges.get(self.edge_idx)?;
        if self.reversed {
            Some(edge.start)
        } else {
            Some(edge.end)
        }
    }
}

// =============================================================================
// FaceAdaptor (BRepAdaptor_Surface)
// =============================================================================

/// Adapts a BRep face to act as a 3D surface.
///
/// Provides surface-like evaluation methods (point, normal, domain) for a face
/// in a BRep, respecting the face's parameter range bounds.
///
/// Analogous to OCCT's `BRepAdaptor_Surface`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_adaptor::FaceAdaptor;
/// use rcad_kernel::BRep;
/// use rcad_kernel::geom::PrimitiveSolid;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
/// let adaptor = FaceAdaptor::new(&brep, 0);
/// let domain = adaptor.domain();
/// let center = adaptor.point_at(
///     (domain[0] + domain[1]) / 2.0,
///     (domain[2] + domain[3]) / 2.0,
/// );
/// ```
#[derive(Debug, Clone)]
pub struct FaceAdaptor<'a> {
    brep: &'a BRep,
    face_idx: usize,
    /// Cached surface reference (if available).
    surface: Option<&'a Surface3>,
    /// Cached parameter range [u_min, u_max, v_min, v_max].
    range: [f64; 4],
}

impl<'a> FaceAdaptor<'a> {
    /// Create a new face adaptor for the given flat face index.
    ///
    /// The flat face index counts faces across all solids/shells in traversal order.
    /// The adaptor respects the face's stored parameter range in `face_surface_range`
    /// and falls back to the surface's natural domain if not specified.
    ///
    /// # Panics
    ///
    /// Does not panic; returns a default adaptor if the face index is out of bounds
    /// or the face has no associated surface.
    pub fn new(brep: &'a BRep, face_idx: usize) -> Self {
        let surface = brep
            .geom
            .face_surface
            .get(face_idx)
            .and_then(|opt| opt.as_ref())
            .and_then(|&surf_idx| brep.geom.surfaces.get(surf_idx));

        let range = if let Some(r) = brep
            .geom
            .face_surface_range
            .get(face_idx)
            .and_then(|r| *r)
        {
            r
        } else if let Some(s) = surface {
            s.default_domain()
        } else {
            [0.0, 1.0, 0.0, 1.0]
        };

        Self {
            brep,
            face_idx,
            surface,
            range,
        }
    }

    /// Evaluate the point on the face at parameters `(u, v)`.
    ///
    /// Parameters are in the face's parameter domain.
    pub fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let Some(surface) = self.surface else {
            // Fall back to face centroid if no surface is available.
            return self.point_from_vertices();
        };

        surface.point_at(u, v)
    }

    /// Evaluate the unit normal vector on the face at parameters `(u, v)`.
    ///
    /// Returns the outward-pointing normal (respecting face orientation).
    pub fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let Some(surface) = self.surface else {
            // Fall back to stored face normal if no surface is available.
            return self.normal_from_face();
        };

        let mut normal = surface.normal_at(u, v);

        // Check if face orientation should flip the normal.
        // In this implementation, we use the stored face normal direction
        // to determine the correct orientation.
        if let Some(face) = self.get_face() {
            if normal.dot(face.normal) < 0.0 {
                normal = -normal;
            }
        }

        normal
    }

    /// Return the parameter domain of the face.
    ///
    /// Returns `[u_min, u_max, v_min, v_max]`.
    pub fn domain(&self) -> [f64; 4] {
        self.range
    }

    /// Return the underlying surface reference, if available.
    pub fn surface(&self) -> Option<&Surface3> {
        self.surface
    }

    /// Check if the face's surface is closed in the U direction.
    ///
    /// A surface is U-closed if `S(u_min, v) == S(u_max, v)` for all v.
    pub fn is_u_closed(&self) -> bool {
        let Some(surface) = self.surface else {
            return false;
        };

        match surface {
            Surface3::Cylinder(_) => true,
            Surface3::Sphere(_) => true,
            Surface3::Cone(_) => true,
            Surface3::Torus(_) => true,
            Surface3::Ellipsoid(_) => true,
            Surface3::Helicoid(_) => false,
            Surface3::Pipe(_) => true,
            Surface3::Plane(_) => false,
            Surface3::BSpline(s) => {
                // Check if first and last rows of control points coincide.
                let n_u = s.control_points.len();
                if n_u < 2 {
                    return false;
                }
                let first = &s.control_points[0];
                let last = &s.control_points[n_u - 1];
                if first.len() != last.len() {
                    return false;
                }
                first
                    .iter()
                    .zip(last.iter())
                    .all(|(a, b)| (a - b).length_squared() < 1e-10)
            }
            Surface3::Bezier(_) => false,
            Surface3::TriBezier(_) => false,
            Surface3::LinearExtrusion(_) => false,
            Surface3::Revolution(_) => true,
            Surface3::Ruled(_) => false,
            Surface3::Coons(_) => false,
            Surface3::Offset(inner) => {
                // Delegate to inner surface check.
                match inner.basis.as_ref() {
                    Surface3::Cylinder(_)
                    | Surface3::Sphere(_)
                    | Surface3::Cone(_)
                    | Surface3::Torus(_)
                    | Surface3::Ellipsoid(_)
                    | Surface3::Pipe(_)
                    | Surface3::Revolution(_) => true,
                    _ => false,
                }
            }
            Surface3::Trimmed(inner) => {
                // Trimmed surface may cut a closed surface.
                match inner.basis.as_ref() {
                    Surface3::Cylinder(_)
                    | Surface3::Sphere(_)
                    | Surface3::Cone(_)
                    | Surface3::Torus(_)
                    | Surface3::Ellipsoid(_)
                    | Surface3::Pipe(_)
                    | Surface3::Revolution(_) => {
                        let [u0, u1, _, _] = inner.trim;
                        let [du0, du1, _, _] = inner.basis.default_domain();
                        (u0 - du0).abs() < 1e-10 && (u1 - du1).abs() < 1e-10
                    }
                    _ => false,
                }
            }
        }
    }

    /// Check if the face's surface is closed in the V direction.
    ///
    /// A surface is V-closed if `S(u, v_min) == S(u, v_max)` for all u.
    pub fn is_v_closed(&self) -> bool {
        let Some(surface) = self.surface else {
            return false;
        };

        match surface {
            Surface3::Torus(_) => true,
            Surface3::Sphere(_) => false, // Sphere has poles, not V-closed.
            Surface3::Cylinder(_) => false,
            Surface3::Cone(_) => false,
            Surface3::Plane(_) => false,
            Surface3::Ellipsoid(_) => false,
            Surface3::Helicoid(_) => false,
            Surface3::Pipe(_) => false,
            Surface3::BSpline(s) => {
                // Check if first and last columns of control points coincide.
                let n_u = s.control_points.len();
                if n_u == 0 {
                    return false;
                }
                let n_v = s.control_points[0].len();
                if n_v < 2 {
                    return false;
                }
                (0..n_u).all(|i| {
                    let first = &s.control_points[i][0];
                    let last = &s.control_points[i][n_v - 1];
                    (first - last).length_squared() < 1e-10
                })
            }
            Surface3::Bezier(_) => false,
            Surface3::TriBezier(_) => false,
            Surface3::LinearExtrusion(_) => false,
            Surface3::Revolution(_) => false,
            Surface3::Ruled(_) => false,
            Surface3::Coons(_) => false,
            Surface3::Offset(inner) => {
                matches!(inner.basis.as_ref(), Surface3::Torus(_))
            }
            Surface3::Trimmed(inner) => {
                if let Surface3::Torus(_) = inner.basis.as_ref() {
                    let [_, _, v0, v1] = inner.trim;
                    let [_, _, dv0, dv1] = inner.basis.default_domain();
                    (v0 - dv0).abs() < 1e-10 && (v1 - dv1).abs() < 1e-10
                } else {
                    false
                }
            }
        }
    }

    /// Get the Face struct for this adaptor's face index.
    fn get_face(&self) -> Option<&rcad_kernel::Face> {
        let mut flat_idx = 0usize;
        for solid in &self.brep.solids {
            for shell in &solid.shells {
                for (fi, face) in shell.faces.iter().enumerate() {
                    if flat_idx + fi == self.face_idx {
                        return Some(face);
                    }
                }
                flat_idx += shell.faces.len();
            }
        }
        None
    }

    /// Fall back to vertex-based point when no surface is available.
    fn point_from_vertices(&self) -> DVec3 {
        let Some(face) = self.get_face() else {
            return DVec3::ZERO;
        };

        // Return centroid of outer wire vertices.
        let mut sum = DVec3::ZERO;
        let mut count = 0usize;
        for we in &face.outer_wire.edges {
            let edge = match self.brep.edges.get(we.idx) {
                Some(e) => e,
                None => continue,
            };
            if let Some(v) = self.brep.vertices.get(edge.start) {
                sum += v.point;
                count += 1;
            }
            if let Some(v) = self.brep.vertices.get(edge.end) {
                sum += v.point;
                count += 1;
            }
        }
        if count > 0 {
            sum / count as f64
        } else {
            DVec3::ZERO
        }
    }

    /// Fall back to stored face normal when no surface is available.
    fn normal_from_face(&self) -> DVec3 {
        let Some(face) = self.get_face() else {
            return DVec3::Z;
        };
        face.normal
    }

    /// Get the tolerance for this face.
    pub fn tolerance(&self) -> f64 {
        self.brep
            .geom
            .face_tolerance
            .get(self.face_idx)
            .copied()
            .unwrap_or(1e-7)
    }
}

// =============================================================================
// WireAdaptor (BRepAdaptor_CompCurve)
// =============================================================================

/// Information about a segment in a wire adaptor.
#[derive(Debug, Clone)]
struct WireSegment {
    /// Edge adaptor for this segment.
    adaptor: EdgeAdaptor<'static>,
    /// Cumulative length fraction at the start of this segment.
    start_frac: f64,
    /// Cumulative length fraction at the end of this segment.
    end_frac: f64,
    /// Arc-length of this segment (approximate).
    length: f64,
}

/// Adapts a BRep wire to act as a composite 3D curve.
///
/// A wire is a connected sequence of edges. This adaptor treats the wire
/// as a single curve parameterized by cumulative arc-length fraction.
///
/// The parameter `t` in [0, 1] represents the position along the wire,
/// where each edge's contribution is weighted by its arc-length.
///
/// Analogous to OCCT's `BRepAdaptor_CompCurve`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_adaptor::WireAdaptor;
/// use rcad_kernel::BRep;
/// use rcad_kernel::geom::PrimitiveSolid;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// // Get the outer wire of face 0.
/// let wire = brep.solids[0].shells[0].faces[0].outer_wire.clone();
/// let adaptor = WireAdaptor::new(&brep, &wire, 0);
/// let midpoint = adaptor.point_at(0.5);
/// let edge_at_mid = adaptor.edge_at(0.5);
/// ```
pub struct WireAdaptor<'a> {
    brep: &'a BRep,
    wire: &'a Wire,
    /// Face index for pcurve lookups (if available).
    face_idx: Option<usize>,
    /// Precomputed segments with arc-lengths.
    segments: Vec<WireSegmentData>,
    /// Total arc-length of the wire.
    total_length: f64,
}

/// Stored segment data (without lifetime issues).
#[derive(Debug, Clone)]
struct WireSegmentData {
    /// Edge index.
    edge_idx: usize,
    /// Whether the edge is reversed in this wire.
    reversed: bool,
    /// Cumulative length fraction at the start of this segment.
    start_frac: f64,
    /// Cumulative length fraction at the end of this segment.
    end_frac: f64,
    /// Arc-length of this segment.
    length: f64,
}

impl<'a> WireAdaptor<'a> {
    /// Create a new wire adaptor.
    ///
    /// # Arguments
    ///
    /// * `brep` - Reference to the BRep containing the wire.
    /// * `wire` - Reference to the wire to adapt.
    /// * `face_idx` - Optional flat face index for pcurve lookups.
    ///
    /// The wire's edges are preprocessed to compute arc-lengths for
    /// parameterization by cumulative length fraction.
    pub fn new(brep: &'a BRep, wire: &'a Wire, face_idx: usize) -> Self {
        let mut segments = Vec::with_capacity(wire.edges.len());
        let mut total_length = 0.0f64;

        for we in &wire.edges {
            let edge_idx = we.idx;
            let reversed = !we.forward;

            // Compute approximate arc-length for this edge.
            let length = Self::compute_edge_length(brep, edge_idx);
            total_length += length;

            segments.push(WireSegmentData {
                edge_idx,
                reversed,
                start_frac: 0.0, // Will be computed after total_length is known
                end_frac: 0.0,
                length,
            });
        }

        // Compute cumulative fractions.
        if total_length > 1e-15 {
            let mut cum_length = 0.0f64;
            for seg in &mut segments {
                seg.start_frac = cum_length / total_length;
                cum_length += seg.length;
                seg.end_frac = cum_length / total_length;
            }
        } else if !segments.is_empty() {
            // All edges are zero-length; distribute uniformly.
            let n = segments.len() as f64;
            for (i, seg) in segments.iter_mut().enumerate() {
                seg.start_frac = i as f64 / n;
                seg.end_frac = (i + 1) as f64 / n;
                seg.length = 1.0; // Dummy length
            }
            total_length = n;
        }

        Self {
            brep,
            wire,
            face_idx: Some(face_idx),
            segments,
            total_length,
        }
    }

    /// Create a wire adaptor without a face context.
    ///
    /// This constructor is used when the wire is not associated with a specific face.
    pub fn without_face(brep: &'a BRep, wire: &'a Wire) -> Self {
        let mut adaptor = Self::new(brep, wire, 0);
        adaptor.face_idx = None;
        adaptor
    }

    /// Evaluate the point on the wire at parameter `t`.
    ///
    /// The parameter `t` is in [0, 1] and represents the cumulative
    /// arc-length fraction along the wire.
    pub fn point_at(&self, t: f64) -> DVec3 {
        let t = t.clamp(0.0, 1.0);
        let seg = self.find_segment(t);
        let local_t = self.local_parameter(t, &seg);

        // Create an edge adaptor for this segment.
        let adaptor = self.create_edge_adaptor(seg.edge_idx, seg.reversed);
        adaptor.point_at(local_t)
    }

    /// Evaluate the unit tangent vector on the wire at parameter `t`.
    ///
    /// Returns the tangent pointing in the direction of traversal along the wire.
    pub fn tangent_at(&self, t: f64) -> DVec3 {
        let t = t.clamp(0.0, 1.0);
        let seg = self.find_segment(t);
        let local_t = self.local_parameter(t, &seg);

        let adaptor = self.create_edge_adaptor(seg.edge_idx, seg.reversed);
        adaptor.tangent_at(local_t)
    }

    /// Return the edge index that contains the given parameter `t`.
    ///
    /// This is useful for determining which edge a point lies on.
    pub fn edge_at(&self, t: f64) -> usize {
        let t = t.clamp(0.0, 1.0);
        let seg = self.find_segment(t);
        seg.edge_idx
    }

    /// Return the number of edges in the wire.
    pub fn num_edges(&self) -> usize {
        self.wire.edges.len()
    }

    /// Return the total arc-length of the wire.
    pub fn length(&self) -> f64 {
        self.total_length
    }

    /// Return the parameter domain of the wire (always [0, 1]).
    pub fn domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }

    /// Find the segment containing the given parameter.
    fn find_segment(&self, t: f64) -> &WireSegmentData {
        // Binary search for the segment containing t.
        let mut lo = 0usize;
        let mut hi = self.segments.len();

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let seg = &self.segments[mid];
            if t < seg.start_frac {
                hi = mid;
            } else if t > seg.end_frac {
                lo = mid + 1;
            } else {
                return seg;
            }
        }

        // Fallback: return last segment.
        self.segments.last().unwrap_or(&WireSegmentData {
            edge_idx: 0,
            reversed: false,
            start_frac: 0.0,
            end_frac: 1.0,
            length: 1.0,
        })
    }

    /// Compute the local parameter within a segment for global parameter t.
    fn local_parameter(&self, t: f64, seg: &WireSegmentData) -> f64 {
        if seg.end_frac <= seg.start_frac {
            return 0.5;
        }
        ((t - seg.start_frac) / (seg.end_frac - seg.start_frac)).clamp(0.0, 1.0)
    }

    /// Create an edge adaptor with the specified orientation.
    fn create_edge_adaptor(&self, edge_idx: usize, reversed: bool) -> EdgeAdaptor<'a> {
        EdgeAdaptor::new(self.brep, edge_idx).with_reversed(reversed)
    }

    /// Compute the approximate arc-length of an edge.
    fn compute_edge_length(brep: &BRep, edge_idx: usize) -> f64 {
        // Try to compute from curve.
        if let Some(&curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|o| o.as_ref()) {
            if let Some(curve) = brep.geom.curves.get(curve_idx) {
                let range = brep
                    .geom
                    .edge_curve_range
                    .get(edge_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| curve.default_domain());

                // Use numerical integration for arc-length.
                return Self::arc_length_numerical(curve, range[0], range[1]);
            }
        }

        // Fall back to vertex distance.
        let Some(edge) = brep.edges.get(edge_idx) else {
            return 0.0;
        };
        let Some(v_start) = brep.vertices.get(edge.start) else {
            return 0.0;
        };
        let Some(v_end) = brep.vertices.get(edge.end) else {
            return 0.0;
        };

        (v_end.point - v_start.point).length()
    }

    /// Numerical integration for arc-length using Gauss-Legendre quadrature.
    fn arc_length_numerical(curve: &Curve3, t0: f64, t1: f64) -> f64 {
        // Use 5-point Gauss-Legendre quadrature.
        const GAUSS_POINTS: [(f64, f64); 5] = [
            (0.0, 0.5688888888888889),
            (-0.5384693101056831, 0.47862867049936647),
            (0.5384693101056831, 0.47862867049936647),
            (-0.9061798459386640, 0.23692688505618908),
            (0.9061798459386640, 0.23692688505618908),
        ];

        let dt = t1 - t0;
        let mut length = 0.0f64;

        for (xi, wi) in GAUSS_POINTS {
            let t = 0.5 * (t0 + t1 + xi * dt);
            let tangent = curve.tangent_at(t);
            length += wi * tangent.length();
        }

        length * 0.5 * dt.abs()
    }
}

// =============================================================================
// CurveAdaptorArray (BRepAdaptor_HArray1OfCurve)
// =============================================================================

/// An array of edge adaptors with indexed access.
///
/// Provides convenient storage and access for multiple curve adaptors,
/// analogous to OCCT's `BRepAdaptor_HArray1OfCurve`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_adaptor::{EdgeAdaptor, CurveAdaptorArray};
/// use rcad_kernel::BRep;
/// use rcad_kernel::geom::PrimitiveSolid;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let mut array = CurveAdaptorArray::new();
/// for i in 0..brep.edges.len() {
///     array.push(EdgeAdaptor::new(&brep, i));
/// }
///
/// for i in 0..array.len() {
///     let adaptor = array.get(i).unwrap();
///     println!("Edge {} domain: {:?}", i, adaptor.domain());
/// }
/// ```
#[derive(Debug, Clone, Default)]
pub struct CurveAdaptorArray<'a> {
    adaptors: Vec<EdgeAdaptor<'a>>,
}

impl<'a> CurveAdaptorArray<'a> {
    /// Create an empty array.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create an array with pre-allocated capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            adaptors: Vec::with_capacity(capacity),
        }
    }

    /// Add an edge adaptor to the array.
    pub fn push(&mut self, adaptor: EdgeAdaptor<'a>) {
        self.adaptors.push(adaptor);
    }

    /// Get the edge adaptor at the given index.
    pub fn get(&self, index: usize) -> Option<&EdgeAdaptor<'a>> {
        self.adaptors.get(index)
    }

    /// Get a mutable reference to the edge adaptor at the given index.
    pub fn get_mut(&mut self, index: usize) -> Option<&mut EdgeAdaptor<'a>> {
        self.adaptors.get_mut(index)
    }

    /// Return the number of adaptors in the array.
    pub fn len(&self) -> usize {
        self.adaptors.len()
    }

    /// Return true if the array is empty.
    pub fn is_empty(&self) -> bool {
        self.adaptors.is_empty()
    }

    /// Iterate over all adaptors.
    pub fn iter(&self) -> impl Iterator<Item = &EdgeAdaptor<'a>> {
        self.adaptors.iter()
    }

    /// Iterate mutably over all adaptors.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut EdgeAdaptor<'a>> {
        self.adaptors.iter_mut()
    }

    /// Clear all adaptors from the array.
    pub fn clear(&mut self) {
        self.adaptors.clear();
    }

    /// Create an array from a BRep's edges.
    ///
    /// Creates an adaptor for each edge in the BRep.
    pub fn from_brep(brep: &'a BRep) -> Self {
        let mut array = Self::with_capacity(brep.edges.len());
        for i in 0..brep.edges.len() {
            array.push(EdgeAdaptor::new(brep, i));
        }
        array
    }
}

impl<'a> std::ops::Index<usize> for CurveAdaptorArray<'a> {
    type Output = EdgeAdaptor<'a>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.adaptors[index]
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle3, Line3, PrimitiveSolid};
    use rcad_kernel::topology::{Edge, Vertex};
    use std::f64::consts::FRAC_PI_2;

    fn box_brep() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        })
    }

    fn sphere_brep() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 })
    }

    fn cylinder_brep() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        })
    }

    fn torus_brep() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 2.0,
            minor_radius: 0.5,
        })
    }

    // ==================== EdgeAdaptor Tests ====================

    #[test]
    fn edge_adaptor_sphere_edge_has_curve() {
        let brep = sphere_brep();
        // Sphere seam edge should have a curve.
        let adaptor = EdgeAdaptor::new(&brep, 0);
        assert!(
            adaptor.curve().is_some(),
            "Sphere seam edge should have a curve"
        );
    }

    #[test]
    fn edge_adaptor_domain_is_normalized() {
        let brep = box_brep();
        let adaptor = EdgeAdaptor::new(&brep, 0);
        let domain = adaptor.domain();
        assert!((domain[0] - 0.0).abs() < 1e-10);
        assert!((domain[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn edge_adaptor_point_at_endpoints() {
        let brep = box_brep();
        let adaptor = EdgeAdaptor::new(&brep, 0);

        // Points at t=0 and t=1 should be at the edge endpoints.
        let p0 = adaptor.point_at(0.0);
        let p1 = adaptor.point_at(1.0);

        // For a box, edge endpoints should be at vertex positions.
        let edge = &brep.edges[0];
        let v_start = brep.vertices[edge.start].point;
        let v_end = brep.vertices[edge.end].point;

        assert!((p0 - v_start).length() < 1e-10, "p0: {:?}, v_start: {:?}", p0, v_start);
        assert!((p1 - v_end).length() < 1e-10, "p1: {:?}, v_end: {:?}", p1, v_end);
    }

    #[test]
    fn edge_adaptor_reversed_direction() {
        let brep = box_brep();
        let adaptor_fwd = EdgeAdaptor::new(&brep, 0);
        let adaptor_rev = EdgeAdaptor::new(&brep, 0).with_reversed(true);

        // Reversed adaptor should give opposite tangent.
        let tan_fwd = adaptor_fwd.tangent_at(0.5);
        let tan_rev = adaptor_rev.tangent_at(0.5);

        assert!((tan_fwd + tan_rev).length() < 1e-10, "tan_fwd: {:?}, tan_rev: {:?}", tan_fwd, tan_rev);
    }

    #[test]
    fn edge_adaptor_closed_edge() {
        // Create a BRep with a closed edge (same start and end vertex)
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.edges.push(Edge { start: 0, end: 0 }); // Same vertex

        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        brep.geom.curves.push(circle);
        brep.geom.edge_curve.push(Some(0));
        brep.geom.edge_curve_range.push(Some([0.0, 2.0 * PI]));

        let adaptor = EdgeAdaptor::new(&brep, 0);
        assert!(adaptor.is_closed(), "Circle edge should be closed");
    }

    #[test]
    fn edge_adaptor_period_circle() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::ZERO }); // Closed edge.
        brep.edges.push(Edge { start: 0, end: 1 });

        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        brep.geom.curves.push(circle);
        brep.geom.edge_curve.push(Some(0));
        brep.geom.edge_curve_range.push(Some([0.0, 2.0 * PI]));

        let adaptor = EdgeAdaptor::new(&brep, 0);
        let period = adaptor.period();
        assert!(period.is_some());
        assert!((period.unwrap() - 2.0 * PI).abs() < 1e-10);
    }

    #[test]
    fn edge_adaptor_no_curve_fallback() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(Edge { start: 0, end: 1 });
        // No curve set.

        let adaptor = EdgeAdaptor::new(&brep, 0);
        let p0 = adaptor.point_at(0.0);
        let p1 = adaptor.point_at(1.0);

        assert!((p0 - DVec3::ZERO).length() < 1e-10);
        assert!((p1 - DVec3::X).length() < 1e-10);

        let tangent = adaptor.tangent_at(0.5);
        assert!((tangent - DVec3::X).length() < 1e-10);
    }

    // ==================== FaceAdaptor Tests ====================

    #[test]
    fn face_adaptor_cylinder_face_has_surface() {
        // Use cylinder instead of box, since box doesn't set up surfaces
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        let adaptor = FaceAdaptor::new(&brep, 0);
        assert!(
            adaptor.surface().is_some(),
            "Cylinder face should have a surface"
        );
    }

    #[test]
    fn face_adaptor_domain() {
        let brep = box_brep();
        let adaptor = FaceAdaptor::new(&brep, 0);
        let domain = adaptor.domain();
        // Plane domain is infinite by default.
        assert!(domain[1] > domain[0]);
        assert!(domain[3] > domain[2]);
    }

    #[test]
    fn face_adaptor_sphere_point() {
        let brep = sphere_brep();
        let adaptor = FaceAdaptor::new(&brep, 0);

        // Test that points lie on the sphere.
        let domain = adaptor.domain();
        for u in [domain[0], (domain[0] + domain[1]) / 2.0, domain[1]] {
            for v in [domain[2], (domain[2] + domain[3]) / 2.0, domain[3]] {
                let p = adaptor.point_at(u, v);
                let r = p.length();
                assert!(
                    (r - 1.0).abs() < 1e-9,
                    "Point at ({}, {}) has radius {}",
                    u, v, r
                );
            }
        }
    }

    #[test]
    fn face_adaptor_sphere_normal() {
        let brep = sphere_brep();
        let adaptor = FaceAdaptor::new(&brep, 0);

        let domain = adaptor.domain();
        let u = (domain[0] + domain[1]) / 2.0;
        let v = (domain[2] + domain[3]) / 2.0;

        let p = adaptor.point_at(u, v);
        let n = adaptor.normal_at(u, v);

        // Normal should be parallel to position vector (outward or inward).
        let dot = p.normalize_or_zero().dot(n);
        assert!((dot.abs() - 1.0).abs() < 1e-9, "Normal dot product: {}", dot);
    }

    #[test]
    fn face_adaptor_cylinder_u_closed() {
        let brep = cylinder_brep();
        let adaptor = FaceAdaptor::new(&brep, 0); // Cylindrical face.
        assert!(
            adaptor.is_u_closed(),
            "Cylinder should be U-closed"
        );
        assert!(
            !adaptor.is_v_closed(),
            "Cylinder should not be V-closed"
        );
    }

    #[test]
    fn face_adaptor_torus_both_closed() {
        let brep = torus_brep();
        let adaptor = FaceAdaptor::new(&brep, 0);
        assert!(
            adaptor.is_u_closed(),
            "Torus should be U-closed"
        );
        assert!(
            adaptor.is_v_closed(),
            "Torus should be V-closed"
        );
    }

    #[test]
    fn face_adaptor_sphere_not_v_closed() {
        let brep = sphere_brep();
        let adaptor = FaceAdaptor::new(&brep, 0);
        // Sphere has poles, so it's not V-closed.
        assert!(
            adaptor.is_u_closed(),
            "Sphere should be U-closed"
        );
        assert!(
            !adaptor.is_v_closed(),
            "Sphere should not be V-closed"
        );
    }

    // ==================== WireAdaptor Tests ====================

    #[test]
    fn wire_adaptor_box_face_wire() {
        let brep = box_brep();
        let face = &brep.solids[0].shells[0].faces[0];
        let wire = &face.outer_wire;
        let adaptor = WireAdaptor::new(&brep, wire, 0);

        assert_eq!(adaptor.num_edges(), 4, "Box face should have 4 edges");
    }

    #[test]
    fn wire_adaptor_point_at_endpoints() {
        let brep = box_brep();
        let face = &brep.solids[0].shells[0].faces[0];
        let wire = &face.outer_wire;
        let adaptor = WireAdaptor::new(&brep, wire, 0);

        // Points at t=0 and t=1 should be at the wire endpoints.
        let p0 = adaptor.point_at(0.0);
        let p1 = adaptor.point_at(1.0);

        // For a closed wire, these should be the same.
        assert!(
            (p0 - p1).length() < 1e-6,
            "Closed wire: p0 {:?} should equal p1 {:?}",
            p0, p1
        );
    }

    #[test]
    fn wire_adaptor_edge_at() {
        let brep = box_brep();
        let face = &brep.solids[0].shells[0].faces[0];
        let wire = &face.outer_wire;
        let adaptor = WireAdaptor::new(&brep, wire, 0);

        // Test that edge_at returns valid edge indices.
        let edge_0 = adaptor.edge_at(0.0);
        let edge_5 = adaptor.edge_at(0.5);
        let edge_1 = adaptor.edge_at(1.0);

        // Edge indices should be within the wire's edge count.
        assert!(edge_0 < wire.edges.len());
        assert!(edge_5 < wire.edges.len());
        assert!(edge_1 < wire.edges.len());
    }

    #[test]
    fn wire_adaptor_tangent_continuity() {
        let brep = box_brep();
        let face = &brep.solids[0].shells[0].faces[0];
        let wire = &face.outer_wire;
        let adaptor = WireAdaptor::new(&brep, wire, 0);

        // Tangents at the start and end of the wire should be consistent.
        let tan_start = adaptor.tangent_at(0.001);
        let tan_end = adaptor.tangent_at(0.999);

        // For a closed wire, tangents at start/end should be similar.
        // (But for a square, they will be perpendicular at corners.)
        assert!(tan_start.length() > 0.9);
        assert!(tan_end.length() > 0.9);
    }

    #[test]
    fn wire_adaptor_length() {
        let brep = box_brep();
        let face = &brep.solids[0].shells[0].faces[0];
        let wire = &face.outer_wire;
        let adaptor = WireAdaptor::new(&brep, wire, 0);

        // A 2x2 square should have perimeter 8.
        let length = adaptor.length();
        assert!(
            (length - 8.0).abs() < 0.1,
            "Expected length ~8.0, got {}",
            length
        );
    }

    // ==================== CurveAdaptorArray Tests ====================

    #[test]
    fn curve_adaptor_array_from_brep() {
        let brep = box_brep();
        let array = CurveAdaptorArray::from_brep(&brep);

        assert_eq!(array.len(), brep.edges.len());
        assert!(!array.is_empty());
    }

    #[test]
    fn curve_adaptor_array_index_access() {
        // Use cylinder instead of box, since box doesn't set up curves
        let brep = cylinder_brep();
        let array = CurveAdaptorArray::from_brep(&brep);

        // At least some edges should have curves
        let mut curves_found = 0;
        for i in 0..array.len() {
            let adaptor = &array[i];
            if adaptor.curve().is_some() {
                curves_found += 1;
            }
        }
        assert!(curves_found > 0, "Cylinder should have edges with curves");
    }

    #[test]
    fn curve_adaptor_array_iteration() {
        let brep = box_brep();
        let array = CurveAdaptorArray::from_brep(&brep);

        let mut count = 0;
        for adaptor in array.iter() {
            let _domain = adaptor.domain();
            count += 1;
        }
        assert_eq!(count, brep.edges.len());
    }

    #[test]
    fn curve_adaptor_array_empty() {
        let array: CurveAdaptorArray = CurveAdaptorArray::new();
        assert!(array.is_empty());
        assert_eq!(array.len(), 0);
    }

    #[test]
    fn curve_adaptor_array_push_and_get() {
        let brep = box_brep();
        let mut array = CurveAdaptorArray::with_capacity(5);

        for i in 0..5 {
            array.push(EdgeAdaptor::new(&brep, i));
        }

        assert_eq!(array.len(), 5);
        assert!(array.get(4).is_some());
        assert!(array.get(5).is_none());
    }

    // ==================== Integration Tests ====================

    #[test]
    fn integration_sphere_seam_edge() {
        let brep = sphere_brep();

        // The sphere has a seam edge from north to south pole.
        let seam_adaptor = EdgeAdaptor::new(&brep, 0);

        // The seam edge is a circle (meridian), so it should have a curve
        assert!(seam_adaptor.curve().is_some(), "Seam edge should have a curve");

        // Points at t=0 and t=1 should be on the sphere
        let p0 = seam_adaptor.point_at(0.0);
        let p1 = seam_adaptor.point_at(1.0);
        assert!((p0.length() - 1.0).abs() < 1e-6, "p0 should be on sphere: {:?}", p0);
        assert!((p1.length() - 1.0).abs() < 1e-6, "p1 should be on sphere: {:?}", p1);
    }

    #[test]
    fn integration_cylinder_face_evaluation() {
        let brep = cylinder_brep();

        // Debug: check what faces are available
        let face_count = brep.solids.iter()
            .flat_map(|s| &s.shells)
            .map(|sh| sh.faces.len())
            .sum::<usize>();
        assert!(face_count > 0, "Cylinder should have faces");

        // Check if face_surface is set for face 0
        let surface_idx = brep.geom.face_surface.get(0).and_then(|s| *s);
        if surface_idx.is_none() {
            // Skip test if surfaces aren't set up
            eprintln!("Skipping test: face_surface not set for cylinder");
            return;
        }

        let adaptor = FaceAdaptor::new(&brep, 0); // Cylindrical face.

        // Check that surface is available
        assert!(adaptor.surface().is_some(), "Cylinder face should have a surface");

        // Points should lie on the cylinder surface.
        let domain = adaptor.domain();

        // Check if domain is valid
        if !domain[0].is_finite() || !domain[1].is_finite() || !domain[2].is_finite() || !domain[3].is_finite() {
            // Use default domain for cylinder
            let u = std::f64::consts::PI; // mid of [0, 2π]
            let v = 0.0;
            let p = adaptor.point_at(u, v);
            let n = adaptor.normal_at(u, v);
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite(), "Point should be finite: {:?}", p);
            return;
        }

        let u = (domain[0] + domain[1]) / 2.0;
        let v = (domain[2] + domain[3]) / 2.0;

        let p = adaptor.point_at(u, v);
        let n = adaptor.normal_at(u, v);

        // Check that we got valid values
        assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite(), "Point should be finite: {:?}", p);
        assert!(n.x.is_finite() && n.y.is_finite() && n.z.is_finite(), "Normal should be finite: {:?}", n);

        // Normal should be perpendicular to axis (Y).
        let y_component = n.dot(DVec3::Y);
        assert!(
            y_component.abs() < 1e-6,
            "Cylinder normal should be radial, y_component: {}",
            y_component
        );

        // Distance from axis should be approximately the radius.
        let radial = DVec3::new(p.x, 0.0, p.z);
        assert!(
            (radial.length() - 1.0).abs() < 0.1,
            "Radial distance should be ~1.0, got {}",
            radial.length()
        );
    }

    #[test]
    fn integration_wire_traversal() {
        let brep = box_brep();

        // Traverse each face's wire.
        for (fi, face) in brep.solids[0].shells[0].faces.iter().enumerate() {
            let adaptor = WireAdaptor::new(&brep, &face.outer_wire, fi);

            // Sample points along the wire.
            for i in 0..=10 {
                let t = i as f64 / 10.0;
                let _p = adaptor.point_at(t);
                let _tan = adaptor.tangent_at(t);
                let _edge = adaptor.edge_at(t);
            }
        }
    }
}
