//! Exact geometric predicates via Shewchuk adaptive arithmetic.
//!
//! Thin wrappers around the [`robust`] crate that match rmsh conventions:
//!
//! | Function | Convention | Positive return |
//! |---|---|---|
//! | [`orient3d`] | `orient3d(a,b,c,d)` → `6 × signed_volume(a,b,c,d)` | d is below (a,b,c) plane, abc ccw from above |
//! | [`in_sphere`] | d is inside circumsphere of (a,b,c,d) | inside, 0 = cocircular, negative = outside |
//! | [`orient2d`] | `orient2d(a,b,c)` → `2 × signed_area(a,b,c)` | c is left of directed line (a,b) |
//!
//! These are drop-in replacements for float-only versions, guaranteed to
//! return the correct sign for any _well-separated_ input (the full exact
//! expansion is computed on demand when the approximate result is within
//! the error bound).

use robust::{Coord, Coord3D};

// ─── 3-D predicates ──────────────────────────────────────────────────────────

/// Robust 3-D orientation test (equivalent to `6 × tetrahedron signed volume`).
///
/// Returns a positive value if the point `d` lies **below** the oriented plane
/// defined by `a, b, c` (i.e. the tetrahedron `(a,b,c,d)` has positive signed
/// volume).  Returns `0` when the four points are coplanar.
pub fn orient3d(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    robust::orient3d(
        Coord3D { x: a[0], y: a[1], z: a[2] },
        Coord3D { x: b[0], y: b[1], z: b[2] },
        Coord3D { x: c[0], y: c[1], z: c[2] },
        Coord3D { x: d[0], y: d[1], z: d[2] },
    )
}

/// Robust 3-D in-sphere test.
///
/// Returns a positive value if point `e` lies **inside** the circumsphere of
/// tetrahedron `(a,b,c,d)`, negative if outside, zero if exactly on the sphere.
///
/// `a,b,c,d` must be oriented so the tetrahedron has **positive** signed volume
/// (`orient3d(a,b,c,d) > 0`).  This wrapper handles arbitrary orientation by
/// swapping two points when necessary and negating the result accordingly.
pub fn in_sphere(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    e: [f64; 3],
) -> f64 {
    // Normalise orientation: if (a,b,c,d) has negative orientation,
    // swap b ↔ c (which flips the sign of both orient3d and insphere).
    // The robust::insphere call then receives positively-oriented input.
    if orient3d(a, b, c, d) >= 0.0 {
        robust::insphere(
            Coord3D { x: a[0], y: a[1], z: a[2] },
            Coord3D { x: b[0], y: b[1], z: b[2] },
            Coord3D { x: c[0], y: c[1], z: c[2] },
            Coord3D { x: d[0], y: d[1], z: d[2] },
            Coord3D { x: e[0], y: e[1], z: e[2] },
        )
    } else {
        robust::insphere(
            Coord3D { x: a[0], y: a[1], z: a[2] },
            Coord3D { x: c[0], y: c[1], z: c[2] },
            Coord3D { x: b[0], y: b[1], z: b[2] },
            Coord3D { x: d[0], y: d[1], z: d[2] },
            Coord3D { x: e[0], y: e[1], z: e[2] },
        )
    }
}

// ─── 2-D predicates ──────────────────────────────────────────────────────────

/// Robust 2-D orientation test (equivalent to `2 × triangle signed area`).
///
/// Returns a positive value if `(a,b,c)` is counter-clockwise (c lies to the
/// left of directed line ab).  `0` if collinear.
pub fn orient2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    robust::orient2d(
        Coord { x: a[0], y: a[1] },
        Coord { x: b[0], y: b[1] },
        Coord { x: c[0], y: c[1] },
    )
}

/// Robust 2-D in-circle test.
///
/// Returns a positive value if point `d` lies **inside** the circumcircle of
/// triangle `(a,b,c)`, negative if outside. `a,b,c` must be counter-clockwise
/// (`orient2d(a,b,c) > 0`).
///
/// The wrapper normalises arbitrary orientation by swapping `b ↔ c`.
pub fn in_circle(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    d: [f64; 2],
) -> f64 {
    if orient2d(a, b, c) >= 0.0 {
        robust::incircle(
            Coord { x: a[0], y: a[1] },
            Coord { x: b[0], y: b[1] },
            Coord { x: c[0], y: c[1] },
            Coord { x: d[0], y: d[1] },
        )
    } else {
        robust::incircle(
            Coord { x: a[0], y: a[1] },
            Coord { x: c[0], y: c[1] },
            Coord { x: b[0], y: b[1] },
            Coord { x: d[0], y: d[1] },
        )
    }
}

// ─── Helper: triangle area ───────────────────────────────────────────────────

/// Area of triangle `(a,b,c)` via robust orientation (more accurate than cross
/// product for nearly-degenerate triangles).
pub fn triangle_area_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    orient2d(a, b, c).abs() * 0.5
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orient3d_regular_tet_is_positive() {
        // a,b,c in xy-plane (z=0), ccw from +z.
        // robust::orient3d: positive if d is BELOW plane (negative-z side).
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d_below = [0.0, 0.0, -1.0];
        assert!(orient3d(a, b, c, d_below) > 0.0);
        // Point above plane → negative
        let d_above = [0.0, 0.0, 1.0];
        assert!(orient3d(a, b, c, d_above) < 0.0);
        // Swapped orientation flips sign
        assert!(orient3d(a, c, b, d_below) < 0.0);
        assert!(orient3d(a, c, b, d_above) > 0.0);
    }

    #[test]
    fn orient3d_coplanar_is_zero() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.5, 0.5, 0.0]; // coplanar (z=0)
        let res = orient3d(a, b, c, d);
        assert_eq!(res, 0.0, "coplanar should be exactly 0, got {res}");
    }

    #[test]
    fn orient3d_nearly_coplanar_is_not_misclassified() {
        // Large-scale base, tiny Z displacement — robust predicate should not
        // produce zero (which pure f64 would).
        let a = [0.0, 0.0, 0.0];
        let b = [1e10, 0.0, 0.0];
        let c = [0.0, 1e10, 0.0];
        let d = [0.0, 0.0, -1e-10]; // slightly below the plane
        let res = orient3d(a, b, c, d);
        assert_ne!(res, 0.0, "nearly-coplanar must not be exactly 0");
        // Above → positive (below plane)
        assert!(res > 0.0);
    }

    #[test]
    fn in_sphere_regular_tet_center_is_inside() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];
        // Circumcenter: (0.5, 0.5, 0.5)
        let e = [0.5, 0.5, 0.5];
        assert!(in_sphere(a, b, c, d, e) > 0.0);
        // Far away: outside
        let e_far = [10.0, 10.0, 10.0];
        assert!(in_sphere(a, b, c, d, e_far) < 0.0);
    }

    #[test]
    fn in_sphere_degenerate_orientation_still_works() {
        // Tet with negative orientation — wrapper must normalise.
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0]; // swapped relative to regular
        let c = [1.0, 0.0, 0.0];
        let d = [0.0, 0.0, 1.0];
        let e = [0.5, 0.5, 0.5];
        // Should still work regardless of orientation
        let res = in_sphere(a, b, c, d, e);
        assert!(res > 0.0 || res == 0.0, "should be inside, got {res}");
    }

    #[test]
    fn in_sphere_cosphere_is_zero() {
        // Five cocircular points in 2D (z=0 for all) → cospherical degenerate
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = [-1.0, 0.0, 0.0];
        let d = [0.0, -1.0, 0.0];
        let e = [0.0, 0.0, 0.0]; // on the sphere (circle in 2D through a,b,c,d)
        // in_sphere should return 0 (exactly) or very close
        let res = in_sphere(a, b, c, d, e);
        assert!(res.abs() < 1e-15, "cospheric should be ~0, got {res}");
    }

    #[test]
    fn orient2d_ccw_is_positive() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        assert!(orient2d(a, b, c) > 0.0);
        assert!(orient2d(a, c, b) < 0.0);
    }

    #[test]
    fn orient2d_collinear_is_zero() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [2.0, 0.0];
        assert_eq!(orient2d(a, b, c), 0.0);
    }

    #[test]
    fn triangle_area_2d_is_half_of_absolute_orient2d() {
        let a = [0.0, 0.0];
        let b = [2.0, 0.0];
        let c = [1.0, 1.0];
        let area = triangle_area_2d(a, b, c);
        assert!((area - 1.0).abs() < 1e-15, "area should be 1.0, got {area}");
    }

    #[test]
    fn orient3d_agrees_with_signed_volume_for_random_point_set() {
        // Deterministic point set: orient3d/6 should match scalar-triple-product vol.
        let vol = |a:[f64;3],b:[f64;3],c:[f64;3],d:[f64;3]| -> f64 {
            let ad = [a[0]-d[0], a[1]-d[1], a[2]-d[2]];
            let bd = [b[0]-d[0], b[1]-d[1], b[2]-d[2]];
            let cd = [c[0]-d[0], c[1]-d[1], c[2]-d[2]];
            let cross = [
                bd[1]*cd[2] - bd[2]*cd[1],
                bd[2]*cd[0] - bd[0]*cd[2],
                bd[0]*cd[1] - bd[1]*cd[0],
            ];
            (ad[0]*cross[0] + ad[1]*cross[1] + ad[2]*cross[2]) / 6.0
        };
        let cases = [
            // regular unit tet
            ([0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]),
            // large scale
            ([0.,0.,0.],[1e8,0.,0.],[0.,1e8,0.],[0.,0.,1e8]),
            // near-coplanar
            ([0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.5,0.5,1e-12]),
            // negative orientation
            ([0.,0.,0.],[0.,1.,0.],[1.,0.,0.],[0.,0.,1.]),
            // very slender
            ([0.,0.,0.],[10.,0.,0.],[0.,0.1,0.],[0.,0.,0.01]),
        ];
        for &(a,b,c,d) in &cases {
            let o3 = orient3d(a, b, c, d);
            let six_vol = vol(a, b, c, d) * 6.0;
            let diff = (o3 - six_vol).abs();
            let scale = o3.abs().max(six_vol.abs()).max(1.0);
            assert!(
                diff / scale < 1e-14 || (o3 == 0.0 && six_vol == 0.0),
                "orient3d={o3} != 6*vol={six_vol}, diff/scale={}",
                diff / scale
            );
        }
    }
}
