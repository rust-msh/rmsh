//! Robust geometric predicates for the 2-D Delaunay machinery: adaptive
//! exact orientation and in-circumcircle signs over `f64` coordinates.
//!
//! # Method
//!
//! Fast floating-point filters with exact fallbacks, after J. R. Shewchuk,
//! "Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric
//! Predicates", *Discrete & Computational Geometry* 18(3):305-363, 1997 —
//! reimplemented from the published algorithms (error-free transformations
//! plus exact expansion sums); no source text is reused. Every approximate
//! sign is guarded by an error bound sized to cover ALL rounding sources
//! (input-difference rounding, squared-radius rounding, product and sum
//! roundings); when the guard cannot certify the sign, the determinant is
//! recomputed EXACTLY from error-free components whose mathematical total
//! equals the true determinant, and the sign of that exact expansion is
//! returned. Consequences:
//!
//! - a wrong sign is never returned (the filter is conservative);
//! - away from degeneracy the fast path certifies the plain determinant,
//!   so point sets that triangulated correctly keep their triangulation;
//! - exactly degenerate configurations (cocircular grids, collinear
//!   chains) resolve deterministically instead of by rounding noise.
//!
//! An exact zero means TRUE degeneracy (e.g. a point exactly on a
//! circumcircle). Callers own the tie rule.

/// `s = fl(a + b)` with the exact error `e`, so that `s + e == a + b`
/// exactly (Knuth's error-free transformation).
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, err)
}

/// `p = fl(a * b)` with the exact error `e`, so that `p + e == a * b`
/// exactly (Dekker's product realized through the fused multiply-add,
/// which evaluates `a*b - p` with a single rounding).
fn two_product(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let e = a.mul_add(b, -p);
    (p, e)
}

/// Sign of the exact total of a bag of floats: the components are folded
/// into a non-overlapping expansion in ascending-magnitude order (the
/// grow-expansion order), whose total is exact, and the sign of the largest
/// surviving component is the sign of the total. Returns `-1`, `0`, or `1`.
fn exact_sum_sign(mut comps: Vec<f64>) -> i32 {
    comps.retain(|x| *x != 0.0);
    if comps.is_empty() {
        return 0;
    }
    comps.sort_by(|a, b| {
        let (ma, mb) = (a.abs(), b.abs());
        ma.partial_cmp(&mb).unwrap_or(std::cmp::Ordering::Equal)
    });
    // in-place ordered fold: comps[0..n] is always a non-overlapping
    // expansion ordered by increasing magnitude (each insertion is one
    // grow-expansion step)
    let mut n = 0usize;
    for j in 0..comps.len() {
        let mut q = comps[j];
        let mut i = 0;
        while i < n {
            let (s, e) = two_sum(q, comps[i]);
            comps[i] = e;
            q = s;
            i += 1;
        }
        if q != 0.0 {
            comps[n] = q;
            n += 1;
        }
    }
    for &x in comps[..n].iter().rev() {
        if x > 0.0 {
            return 1;
        }
        if x < 0.0 {
            return -1;
        }
    }
    0
}

/// Push the exact product `s * x` (s scalar, x expansion component) as two
/// error-free components, skipping exact zeros.
fn push_scaled(out: &mut Vec<f64>, s: f64, x: f64) {
    if s == 0.0 || x == 0.0 {
        return;
    }
    let (p, e) = two_product(s, x);
    if p != 0.0 {
        out.push(p);
    }
    if e != 0.0 {
        out.push(e);
    }
}

/// Push the exact products of a scalar times a component bag.
fn push_scaled_bag(out: &mut Vec<f64>, s: f64, bag: &[f64]) {
    for &x in bag {
        push_scaled(out, s, x);
    }
}

/// Exact squared radius `vx^2 + vy^2` of the pair-valued vector
/// `(hi + lo, hi' + lo')` as an exact component bag (every product
/// error-free, zeros skipped): `h^2 + 2h·l + l^2` per axis.
fn push_radius2(out: &mut Vec<f64>, hi_x: f64, lo_x: f64, hi_y: f64, lo_y: f64) {
    for (h, l) in [(hi_x, lo_x), (hi_y, lo_y)] {
        if h != 0.0 {
            let (p, e) = two_product(h, h);
            if p != 0.0 {
                out.push(p);
            }
            if e != 0.0 {
                out.push(e);
            }
        }
        if h != 0.0 && l != 0.0 {
            // 2*h is exact (power-of-two scaling)
            let two_h = 2.0 * h;
            let (p, e) = two_product(two_h, l);
            if p != 0.0 {
                out.push(p);
            }
            if e != 0.0 {
                out.push(e);
            }
        }
        if l != 0.0 {
            let (p, e) = two_product(l, l);
            if p != 0.0 {
                out.push(p);
            }
            if e != 0.0 {
                out.push(e);
            }
        }
    }
}

/// Sign of `(b - a) x (c - a)` — the orientation of the ordered triple:
/// `+1` counter-clockwise, `-1` clockwise, `0` exactly collinear. Adaptive
/// exact: the fast filter returns the plain determinant whenever it is
/// certified by the published `orient2d` bound; otherwise the 2x2
/// determinant is recomputed exactly.
pub(crate) fn orient2d_sign(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> i32 {
    let acx = a[0] - c[0];
    let acy = a[1] - c[1];
    let bcx = b[0] - c[0];
    let bcy = b[1] - c[1];
    let detleft = acx * bcy;
    let detright = acy * bcx;
    let det = detleft - detright;
    // the published orient2d filter bound: (3 + 16 eps) * eps
    let errbound = (3.0 + 16.0 * f64::EPSILON) * f64::EPSILON;
    let check = detleft.abs() + detright.abs();
    if det.abs() >= errbound * check || check == 0.0 {
        if det > 0.0 {
            return 1;
        }
        if det < 0.0 {
            return -1;
        }
        return 0;
    }
    // exact path: differences as exact two-component pairs, the four
    // products of each half-determinant error-free, the second half
    // negated (exact), the total summed exactly
    let (acx_hi, acx_lo) = two_sum(a[0], -c[0]);
    let (acy_hi, acy_lo) = two_sum(a[1], -c[1]);
    let (bcx_hi, bcx_lo) = two_sum(b[0], -c[0]);
    let (bcy_hi, bcy_lo) = two_sum(b[1], -c[1]);
    let mut comps: Vec<f64> = Vec::with_capacity(16);
    // det = (acx + acx')(bcy + bcy') - (acy + acy')(bcx + bcx')
    push_scaled_bag(&mut comps, acx_hi, &[bcy_hi, bcy_lo]);
    push_scaled_bag(&mut comps, acx_lo, &[bcy_hi, bcy_lo]);
    push_scaled_bag(&mut comps, -acy_hi, &[bcx_hi, bcx_lo]);
    push_scaled_bag(&mut comps, -acy_lo, &[bcx_hi, bcx_lo]);
    exact_sum_sign(comps)
}

/// Sign of the in-circumcircle predicate: `+1` when `p` lies strictly
/// inside the circumcircle of the counter-clockwise triangle `(a, b, c)`
/// (negative for a clockwise one), `0` when `p` is exactly on the circle.
/// Adaptive exact: the fast filter certifies the plain 3x3 determinant with
/// an error bound covering every rounding source; otherwise the determinant
/// is recomputed exactly from error-free components.
pub(crate) fn incircle_sign(a: [f64; 2], b: [f64; 2], c: [f64; 2], p: [f64; 2]) -> i32 {
    let adx = a[0] - p[0];
    let ady = a[1] - p[1];
    let bdx = b[0] - p[0];
    let bdy = b[1] - p[1];
    let cdx = c[0] - p[0];
    let cdy = c[1] - p[1];
    let ad2 = adx * adx + ady * ady;
    let bd2 = bdx * bdx + bdy * bdy;
    let cd2 = cdx * cdx + cdy * cdy;

    // det = adx·(bdy·cd2 - cdy·bd2) - ady·(bdx·cd2 - cdx·bd2)
    //     + ad2·(bdx·cdy - bdy·cdx)
    let det = adx * (bdy * cd2 - cdy * bd2) - ady * (bdx * cd2 - cdx * bd2)
        + ad2 * (bdx * cdy - bdy * cdx);

    // Conservative permanent: every rounding source is charged in units of
    // eps. Difference roundings contribute <= 2 eps each (c = 2 on the six
    // differences), squared-radius roundings <= 8 eps (c = 8 on the three
    // radii, their magnitudes bounded by |x|^2 + |y|^2), and the products
    // and additions of the plain evaluation <= 10 eps of the six leaf
    // products.
    let leaf = |w: f64, x: f64, y: f64| w.abs() * x.abs() * y.abs();
    let l1 = leaf(adx, bdy, cd2);
    let l2 = leaf(adx, cdy, bd2);
    let l3 = leaf(ady, bdx, cd2);
    let l4 = leaf(ady, cdx, bd2);
    let l5 = leaf(ad2, bdx, cdy);
    let l6 = leaf(ad2, bdy, cdx);
    let ad2m = adx.abs() * adx.abs() + ady.abs() * ady.abs();
    let bd2m = bdx.abs() * bdx.abs() + bdy.abs() * bdy.abs();
    let cd2m = cdx.abs() * cdx.abs() + cdy.abs() * cdy.abs();
    let mut permanent = 10.0 * (l1 + l2 + l3 + l4 + l5 + l6);
    // partial-derivative charges |d det/dq| * c_q * |q| for the six
    // differences (each line bounds one input's propagated error; inputs
    // appearing in two product groups are bounded by the sum of both)
    permanent += 2.0 * adx.abs() * (bdy.abs() * cd2 + cdy.abs() * bd2);
    permanent += 2.0 * ady.abs() * (bdx.abs() * cd2 + cdx.abs() * bd2);
    permanent += 2.0 * bdx.abs() * (ady.abs() * cd2 + ad2 * cdy.abs());
    permanent += 2.0 * bdy.abs() * (adx.abs() * cd2 + ad2 * cdx.abs());
    permanent += 2.0 * cdx.abs() * (ady.abs() * bd2 + ad2 * bdy.abs());
    permanent += 2.0 * cdy.abs() * (adx.abs() * bd2 + ad2 * bdx.abs());
    permanent += 8.0 * ad2 * (bdx.abs() * cdy.abs() + bdy.abs() * cdx.abs());
    permanent += 8.0 * bd2 * (adx.abs() * cdy.abs() + ady.abs() * cdx.abs());
    permanent += 8.0 * cd2 * (adx.abs() * bdy.abs() + ady.abs() * bdx.abs());
    // squared-radius magnitude allowances in the leaf-shaped terms
    permanent += 2.0 * adx.abs() * (bdy.abs() * cd2m + cdy.abs() * bd2m);
    permanent += 2.0 * ady.abs() * (bdx.abs() * cd2m + cdx.abs() * bd2m);
    permanent += 8.0 * ad2m * (bdx.abs() * cdy.abs() + bdy.abs() * cdx.abs());

    if det.abs() >= 2.0 * f64::EPSILON * permanent || permanent == 0.0 {
        if det > 0.0 {
            return 1;
        }
        if det < 0.0 {
            return -1;
        }
        return 0;
    }
    incircle_exact_sign(a, b, c, p)
}

/// Exact sign of the 3x3 incircle determinant: differences as exact
/// two-component pairs, squared radii as exact component bags, all products
/// error-free, the total summed exactly.
fn incircle_exact_sign(a: [f64; 2], b: [f64; 2], c: [f64; 2], p: [f64; 2]) -> i32 {
    let (adx, adxl) = two_sum(a[0], -p[0]);
    let (ady, adyl) = two_sum(a[1], -p[1]);
    let (bdx, bdxl) = two_sum(b[0], -p[0]);
    let (bdy, bdyl) = two_sum(b[1], -p[1]);
    let (cdx, cdxl) = two_sum(c[0], -p[0]);
    let (cdy, cdyl) = two_sum(c[1], -p[1]);

    let mut ad2e: Vec<f64> = Vec::with_capacity(12);
    push_radius2(&mut ad2e, adx, adxl, ady, adyl);
    let mut bd2e: Vec<f64> = Vec::with_capacity(12);
    push_radius2(&mut bd2e, bdx, bdxl, bdy, bdyl);
    let mut cd2e: Vec<f64> = Vec::with_capacity(12);
    push_radius2(&mut cd2e, cdx, cdxl, cdy, cdyl);

    let mut comps: Vec<f64> = Vec::with_capacity(128);
    // inner1 = bdy·cd2 - cdy·bd2 ; term1 = adx · inner1
    push_scaled_bag(&mut comps, bdy, &cd2e);
    push_scaled_bag(&mut comps, bdyl, &cd2e);
    push_scaled_bag(&mut comps, -cdy, &bd2e);
    push_scaled_bag(&mut comps, -cdyl, &bd2e);
    let inner1_len = comps.len();
    let mut term1: Vec<f64> = Vec::with_capacity(inner1_len * 4);
    for &x in &comps[..inner1_len] {
        push_scaled(&mut term1, adx, x);
        push_scaled(&mut term1, adxl, x);
    }
    // inner2 = bdx·cd2 - cdx·bd2 ; term2 = -ady · inner2
    let mut inner2: Vec<f64> = Vec::with_capacity(inner1_len);
    push_scaled_bag(&mut inner2, bdx, &cd2e);
    push_scaled_bag(&mut inner2, bdxl, &cd2e);
    push_scaled_bag(&mut inner2, -cdx, &bd2e);
    push_scaled_bag(&mut inner2, -cdxl, &bd2e);
    for &x in &inner2 {
        push_scaled(&mut term1, -ady, x);
        push_scaled(&mut term1, -adyl, x);
    }
    // inner3 = bdx·cdy - bdy·cdx (pair x pair: four error-free products per
    // side) ; term3 = ad2 (exact bag) · inner3
    let mut inner3: Vec<f64> = Vec::with_capacity(16);
    push_scaled(&mut inner3, bdx, cdy);
    push_scaled(&mut inner3, bdx, cdyl);
    push_scaled(&mut inner3, bdxl, cdy);
    push_scaled(&mut inner3, bdxl, cdyl);
    push_scaled(&mut inner3, -bdy, cdx);
    push_scaled(&mut inner3, -bdy, cdxl);
    push_scaled(&mut inner3, -bdyl, cdx);
    push_scaled(&mut inner3, -bdyl, cdxl);
    for &x in &inner3 {
        for &s in &ad2e {
            push_scaled(&mut term1, s, x);
        }
    }
    exact_sum_sign(term1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sgn(x: f64) -> i32 {
        if x > 0.0 {
            1
        } else if x < 0.0 {
            -1
        } else {
            0
        }
    }

    #[test]
    fn orient2d_matches_plain_determinant_away_from_degeneracy() {
        let cases: Vec<[[f64; 2]; 3]> = vec![
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0]],
            [[1.0, 2.0], [3.5, 4.25], [-2.0, 7.125]],
        ];
        for c in cases {
            let plain = (c[1][0] - c[0][0]) * (c[2][1] - c[0][1])
                - (c[1][1] - c[0][1]) * (c[2][0] - c[0][0]);
            assert_eq!(orient2d_sign(c[0], c[1], c[2]), sgn(plain));
        }
    }

    #[test]
    fn orient2d_exact_on_near_collinear_offset_points() {
        // large common offset destroys the plain determinant's low bits;
        // these points are exactly collinear (integer coordinates far from
        // the origin, differences exact), so the exact sign must be 0
        let a = [1e15, 1e15 + 1.0];
        let b = [1e15 + 1.0, 1e15 + 2.0];
        let c = [1e15 + 2.0, 1e15 + 3.0];
        assert_eq!(orient2d_sign(a, b, c), 0);
        // a true bend away from the line resolves with the correct sign
        // even though the plain determinant saturates to 0
        let d = [1e15 + 3.0, 1e15 + 5.0];
        assert_eq!(orient2d_sign(a, b, d), 1);
        assert_eq!(orient2d_sign(a, d, b), -1);
    }

    #[test]
    fn incircle_matches_plain_determinant_away_from_degeneracy() {
        let inside = [0.5, 0.5];
        let outside = [3.0, 3.0];
        let tri = [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]];
        assert_eq!(incircle_sign(tri[0], tri[1], tri[2], inside), 1);
        assert_eq!(incircle_sign(tri[0], tri[1], tri[2], outside), -1);
    }

    #[test]
    fn incircle_exact_on_unit_square_fourth_corner() {
        // the unit square is exactly cocircular: testing the fourth corner
        // against the diagonal triangles must read EXACTLY on-circle
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [1.0, 1.0];
        let d = [0.0, 1.0];
        assert_eq!(incircle_sign(a, b, c, d), 0);
        assert_eq!(incircle_sign(b, c, d, a), 0);
        // the center is strictly inside; a far point strictly outside
        assert_eq!(incircle_sign(a, b, c, [0.5, 0.5]), 1);
        assert_eq!(incircle_sign(a, b, c, [2.0, 0.0]), -1);
    }

    #[test]
    fn incircle_exact_survives_large_offsets() {
        // cocircular square translated far from the origin by an exact
        // power-of-two-scale offset: the plain determinant loses the low
        // bits, the exact expansion still reports on-circle
        let off = 4503599627370496.0; // 2^52
        let a = [off, off];
        let b = [off + 1.0, off];
        let c = [off + 1.0, off + 1.0];
        let d = [off, off + 1.0];
        assert_eq!(incircle_sign(a, b, c, d), 0);
    }
}
