//! The collapsed (truncated pyramidal / Bergot) generator set for the
//! pyramid families: Legendre x Legendre x Jacobi products with the
//! (1-w)^max(i,j) collapse factor, evaluated over the collapsed
//! coordinates (x', y', w) = (2u-1, 2v-1, h). The nodal basis inverts
//! this generator matrix against the pinned layout — Gmsh's
//! pyramidalBasis construction (the pure-polynomial pairing was probed
//! singular: the complete pyramidal space is not polynomial).
//!
//! Reference: Gmsh 4.14 numeric/BergotBasis.cpp +
//! numeric/orthogonalBasis (the Legendre / Jacobi recurrences),
//! algorithm-level port; M. Bergot (2001) for the underlying space.

pub(super) fn legendre_f(n: usize, u: f64) -> Vec<f64> {
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
pub(super) fn legendre_df(n: usize, u: f64) -> Vec<f64> {
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
pub(super) fn jacobi_f(n: usize, alpha: f64, beta: f64, u: f64) -> Vec<f64> {
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
pub(super) fn jacobi_df(n: usize, alpha: f64, beta: f64, u: f64) -> Vec<f64> {
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
pub(super) struct BergotBasis {
    order: usize,
    incomplete: bool,
}

impl BergotBasis {
    pub(super) fn new(order: usize, incomplete: bool) -> Self {
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
    pub(super) fn f(&self, p: [f64; 3]) -> Vec<f64> {
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
    pub(super) fn df(&self, p: [f64; 3]) -> Vec<[f64; 3]> {
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

