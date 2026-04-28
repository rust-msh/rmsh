//! OCCT math-style mathematical utilities.
//!
//! Provides mathematical algorithms and solvers including:
//! - Root finding (Newton-Raphson, bisection, secant)
//! - Multi-dimensional Newton methods
//! - Polynomial solvers (linear through quartic)
//! - Eigenvalue/matrix utilities
//! - Numerical integration
//! - Optimization (golden section search)

use glam::{DMat2, DMat3, DVec2, DVec3};
use std::f64::consts::FRAC_1_SQRT_2;

// =============================================================================
// Root Finding
// =============================================================================

/// Newton-Raphson method for finding roots of a function.
///
/// # Arguments
/// * `f` - The function to find root of
/// * `df` - The derivative of the function
/// * `x0` - Initial guess
/// * `tol` - Tolerance for convergence
/// * `max_iter` - Maximum number of iterations
///
/// # Returns
/// The root if found within tolerance and iteration limit
pub fn newton_raphson(
    f: fn(f64) -> f64,
    df: fn(f64) -> f64,
    x0: f64,
    tol: f64,
    max_iter: usize,
) -> Option<f64> {
    let mut x = x0;

    for _ in 0..max_iter {
        let fx = f(x);
        if fx.abs() < tol {
            return Some(x);
        }

        let dfx = df(x);
        if dfx.abs() < 1e-15 {
            return None; // Derivative too small
        }

        let x_new = x - fx / dfx;

        // Check for convergence
        if (x_new - x).abs() < tol {
            return Some(x_new);
        }

        x = x_new;
    }

    // Check final value
    if f(x).abs() < tol {
        Some(x)
    } else {
        None
    }
}

/// Bisection method for finding roots of a function.
///
/// # Arguments
/// * `f` - The function to find root of
/// * `a` - Lower bound of interval
/// * `b` - Upper bound of interval
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The root if found within interval and tolerance
pub fn bisection(f: fn(f64) -> f64, a: f64, b: f64, tol: f64) -> Option<f64> {
    let mut lo = a;
    let mut hi = b;

    let f_lo = f(lo);
    let f_hi = f(hi);

    // Check if interval brackets a root
    if f_lo * f_hi > 0.0 {
        return None;
    }

    // Check if bounds are already roots
    if f_lo.abs() < tol {
        return Some(lo);
    }
    if f_hi.abs() < tol {
        return Some(hi);
    }

    let max_iter = ((hi - lo) / tol).ceil() as usize;

    for _ in 0..max_iter {
        let mid = (lo + hi) / 2.0;
        let f_mid = f(mid);

        if f_mid.abs() < tol || (hi - lo) / 2.0 < tol {
            return Some(mid);
        }

        if f_lo * f_mid < 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    Some((lo + hi) / 2.0)
}

/// Secant method for finding roots of a function.
///
/// # Arguments
/// * `f` - The function to find root of
/// * `x0` - First initial guess
/// * `x1` - Second initial guess
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The root if found within tolerance
pub fn secant(f: fn(f64) -> f64, x0: f64, x1: f64, tol: f64) -> Option<f64> {
    let mut x_prev = x0;
    let mut x_curr = x1;

    let max_iter = 100;

    for _ in 0..max_iter {
        let f_prev = f(x_prev);
        let f_curr = f(x_curr);

        if f_curr.abs() < tol {
            return Some(x_curr);
        }

        let denom = f_curr - f_prev;
        if denom.abs() < 1e-15 {
            return None;
        }

        let x_new = x_curr - f_curr * (x_curr - x_prev) / denom;

        if (x_new - x_curr).abs() < tol {
            return Some(x_new);
        }

        x_prev = x_curr;
        x_curr = x_new;
    }

    if f(x_curr).abs() < tol {
        Some(x_curr)
    } else {
        None
    }
}

// =============================================================================
// Multi-dimensional Newton Methods
// =============================================================================

/// Newton-Raphson method for 2D systems of equations.
///
/// Solves the system F(x) = 0 where F: R^2 -> R^2
///
/// # Arguments
/// * `f` - The function vector F(x)
/// * `jacobian` - The Jacobian matrix of F
/// * `x0` - Initial guess
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The root vector if found
pub fn newton_2d(
    f: fn(DVec2) -> DVec2,
    jacobian: fn(DVec2) -> DMat2,
    x0: DVec2,
    tol: f64,
) -> Option<DVec2> {
    let mut x = x0;
    let max_iter = 50;

    for _ in 0..max_iter {
        let fx = f(x);

        if fx.length() < tol {
            return Some(x);
        }

        let j = jacobian(x);
        let det = j.determinant();

        if det.abs() < 1e-15 {
            return None; // Singular Jacobian
        }

        let j_inv = j.inverse();
        let delta = j_inv * fx;

        let x_new = x - delta;

        if delta.length() < tol {
            return Some(x_new);
        }

        x = x_new;
    }

    if f(x).length() < tol {
        Some(x)
    } else {
        None
    }
}

/// Newton-Raphson method for 3D systems of equations.
///
/// Solves the system F(x) = 0 where F: R^3 -> R^3
///
/// # Arguments
/// * `f` - The function vector F(x)
/// * `jacobian` - The Jacobian matrix of F
/// * `x0` - Initial guess
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The root vector if found
pub fn newton_3d(
    f: fn(DVec3) -> DVec3,
    jacobian: fn(DVec3) -> DMat3,
    x0: DVec3,
    tol: f64,
) -> Option<DVec3> {
    let mut x = x0;
    let max_iter = 50;

    for _ in 0..max_iter {
        let fx = f(x);

        if fx.length() < tol {
            return Some(x);
        }

        let j = jacobian(x);
        let det = j.determinant();

        if det.abs() < 1e-15 {
            return None; // Singular Jacobian
        }

        if let Some(j_inv) = inverse_3x3(j) {
            let delta = j_inv * fx;
            let x_new = x - delta;

            if delta.length() < tol {
                return Some(x_new);
            }

            x = x_new;
        } else {
            return None;
        }
    }

    if f(x).length() < tol {
        Some(x)
    } else {
        None
    }
}

// =============================================================================
// Polynomial Solvers
// =============================================================================

/// Solve linear equation ax + b = 0
pub fn solve_linear(a: f64, b: f64) -> Option<f64> {
    if a.abs() < 1e-15 {
        if b.abs() < 1e-15 {
            Some(0.0) // Infinite solutions, return 0
        } else {
            None // No solution
        }
    } else {
        Some(-b / a)
    }
}

/// Solve quadratic equation ax^2 + bx + c = 0
///
/// Returns real roots in ascending order
pub fn solve_quadratic(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() < 1e-15 {
        // Linear case
        return solve_linear(b, c).into_iter().collect();
    }

    let disc = b * b - 4.0 * a * c;

    if disc < 0.0 {
        return Vec::new(); // No real roots
    }

    if disc.abs() < 1e-15 {
        // Single root (double root)
        return vec![-b / (2.0 * a)];
    }

    let sqrt_disc = disc.sqrt();
    let q = if b >= 0.0 {
        -0.5 * (b + sqrt_disc)
    } else {
        -0.5 * (b - sqrt_disc)
    };

    let mut roots = vec![q / a, c / q];
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots
}

/// Solve cubic equation ax^3 + bx^2 + cx + d = 0
///
/// Uses Cardano's formula and returns real roots in ascending order
pub fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    if a.abs() < 1e-15 {
        return solve_quadratic(b, c, d);
    }

    // Normalize to x^3 + px^2 + qx + r = 0
    let p = b / a;
    let q = c / a;
    let r = d / a;

    // Substitute x = t - p/3 to get t^3 + at + b = 0
    let a_coef = q - p * p / 3.0;
    let b_coef = 2.0 * p * p * p / 27.0 - p * q / 3.0 + r;

    let disc = b_coef * b_coef / 4.0 + a_coef * a_coef * a_coef / 27.0;

    let offset = p / 3.0;

    if disc.abs() < 1e-15 {
        // One or two roots (discriminant near zero)
        if a_coef.abs() < 1e-15 {
            // Triple root
            return vec![-offset];
        }
        let t = 3.0 * b_coef / a_coef;
        let mut roots = vec![-offset + t, -offset - t / 2.0];
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        roots
    } else if disc > 0.0 {
        // One real root
        let sqrt_disc = disc.sqrt();
        let u = cube_root(-b_coef / 2.0 + sqrt_disc);
        let v = cube_root(-b_coef / 2.0 - sqrt_disc);
        vec![u + v - offset]
    } else {
        // Three real roots (trigonometric solution)
        let m = 2.0 * (-a_coef / 3.0).sqrt();
        let theta = (-b_coef / 2.0) / ((-a_coef * a_coef * a_coef / 27.0).sqrt());

        let theta = theta.clamp(-1.0, 1.0);
        let theta = theta.acos() / 3.0;

        let mut roots = vec![
            m * theta.cos() - offset,
            m * (theta + std::f64::consts::TAU / 3.0).cos() - offset,
            m * (theta + 2.0 * std::f64::consts::TAU / 3.0).cos() - offset,
        ];
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        roots
    }
}

/// Real cube root
fn cube_root(x: f64) -> f64 {
    if x >= 0.0 {
        x.cbrt()
    } else {
        -(-x).cbrt()
    }
}

/// Solve quartic equation ax^4 + bx^3 + cx^2 + dx + e = 0
///
/// Uses Ferrari's method and returns real roots in ascending order
pub fn solve_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Vec<f64> {
    if a.abs() < 1e-15 {
        return solve_cubic(b, c, d, e);
    }

    // Normalize to x^4 + px^3 + qx^2 + rx + s = 0
    let p = b / a;
    let q = c / a;
    let r = d / a;
    let s = e / a;

    // Substitute x = y - p/4 to get depressed quartic y^4 + py^2 + qy + r = 0
    let a1 = q - 3.0 * p * p / 8.0;
    let b1 = r + p * p * p / 8.0 - p * q / 2.0;
    let c1 = s - 3.0 * p * p * p * p / 256.0 + p * p * q / 16.0 - p * r / 4.0;

    // Handle b1 = 0 case (quartic with only even powers)
    if b1.abs() < 1e-12 {
        let disc = a1 * a1 - 4.0 * c1;
        if disc < -1e-10 {
            return Vec::new();
        }
        if disc.abs() < 1e-10 {
            let y = (-a1 / 2.0).sqrt();
            return vec![y - p / 4.0, -y - p / 4.0];
        }
        let sqrt_disc = disc.sqrt();
        let y1_sq = (-a1 + sqrt_disc) / 2.0;
        let y2_sq = (-a1 - sqrt_disc) / 2.0;

        let mut roots = Vec::new();
        if y1_sq >= -1e-10 {
            let y1 = y1_sq.max(0.0).sqrt();
            roots.push(y1 - p / 4.0);
            roots.push(-y1 - p / 4.0);
        }
        if y2_sq >= -1e-10 {
            let y2 = y2_sq.max(0.0).sqrt();
            roots.push(y2 - p / 4.0);
            roots.push(-y2 - p / 4.0);
        }
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        return roots;
    }

    // Solve resolvent cubic: t^3 + 2*a1*t^2 + (a1^2 - 4*c1)*t - b1^2 = 0
    let resolvent_roots = solve_cubic(1.0, 2.0 * a1, a1 * a1 - 4.0 * c1, -b1 * b1);

    if resolvent_roots.is_empty() {
        return Vec::new();
    }

    // Find a positive root of the resolvent
    let t = resolvent_roots
        .iter()
        .find(|&&t| t > 1e-10)
        .copied()
        .unwrap_or(resolvent_roots[0]);

    let sqrt_t = t.max(0.0).sqrt();

    let mut roots = Vec::new();

    if sqrt_t > 1e-10 {
        let inner1 = -(a1 + t + b1 / sqrt_t);
        let inner2 = -(a1 + t - b1 / sqrt_t);

        if inner1 >= -1e-10 {
            let s1 = inner1.max(0.0).sqrt();
            roots.push((sqrt_t + s1) / 2.0 - p / 4.0);
            roots.push((sqrt_t - s1) / 2.0 - p / 4.0);
        }
        if inner2 >= -1e-10 {
            let s2 = inner2.max(0.0).sqrt();
            roots.push((-sqrt_t + s2) / 2.0 - p / 4.0);
            roots.push((-sqrt_t - s2) / 2.0 - p / 4.0);
        }
    } else {
        // t is nearly zero, use alternative formula
        let inner = -(a1 + t);
        if inner >= -1e-10 {
            let s = inner.max(0.0).sqrt();
            roots.push(s / 2.0 - p / 4.0);
            roots.push(-s / 2.0 - p / 4.0);
        }
    }

    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots.dedup_by(|a, b| (*a - *b).abs() < 1e-10);
    roots
}

// =============================================================================
// Eigenvalue/Matrix Utilities
// =============================================================================

/// Compute eigenvalues of a 2x2 matrix.
///
/// Returns eigenvalues in descending order (largest first)
pub fn eigenvalues_2x2(m: DMat2) -> (f64, f64) {
    // For matrix [[a, b], [c, d]]:
    // eigenvalues satisfy: lambda^2 - (a+d)*lambda + (ad-bc) = 0
    let trace = m.x_axis.x + m.y_axis.y; // trace = sum of diagonal
    let det = m.determinant();

    let disc = trace * trace - 4.0 * det;

    if disc < 0.0 {
        // Complex eigenvalues - return real parts
        let real_part = trace / 2.0;
        (real_part, real_part)
    } else {
        let sqrt_disc = disc.sqrt();
        let e1 = (trace + sqrt_disc) / 2.0;
        let e2 = (trace - sqrt_disc) / 2.0;
        if e1 >= e2 {
            (e1, e2)
        } else {
            (e2, e1)
        }
    }
}

/// Compute eigenvalues of a 3x3 matrix using characteristic polynomial.
///
/// Returns eigenvalues sorted in descending order (largest first)
pub fn eigenvalues_3x3(m: DMat3) -> (f64, f64, f64) {
    // Characteristic polynomial: det(A - lambda*I) = 0
    // For 3x3: -lambda^3 + tr(A)*lambda^2 - S*lambda + det(A) = 0
    // where S = sum of principal minors

    let trace = m.x_axis.x + m.y_axis.y + m.z_axis.z; // trace = sum of diagonal
    let det = m.determinant();

    // Sum of principal 2x2 minors:
    // M11 = (a22*a33 - a23*a32), M22 = (a11*a33 - a13*a31), M33 = (a11*a22 - a12*a21)
    let s = (m.y_axis.y * m.z_axis.z - m.y_axis.z * m.z_axis.y)  // M11
          + (m.x_axis.x * m.z_axis.z - m.x_axis.z * m.z_axis.x)  // M22
          + (m.x_axis.x * m.y_axis.y - m.x_axis.y * m.y_axis.x); // M33

    // Coefficients of characteristic polynomial: lambda^3 - trace*lambda^2 + s*lambda - det = 0
    let roots = solve_cubic(1.0, -trace, s, -det);

    match roots.len() {
        0 => (0.0, 0.0, 0.0),
        1 => (roots[0], roots[0], roots[0]),
        2 => {
            if roots[0] >= roots[1] {
                (roots[0], roots[1], roots[1])
            } else {
                (roots[1], roots[0], roots[0])
            }
        }
        3 => {
            let mut r = roots;
            r.sort_by(|a, b| b.partial_cmp(a).unwrap()); // Descending
            (r[0], r[1], r[2])
        }
        _ => (0.0, 0.0, 0.0),
    }
}

/// Compute the inverse of a 3x3 matrix.
///
/// Returns None if the matrix is singular
pub fn inverse_3x3(m: DMat3) -> Option<DMat3> {
    let det = determinant_3x3(m);
    if det.abs() < 1e-15 {
        return None;
    }

    // Cofactor matrix (transpose of adjugate)
    let cofactor00 = m.y_axis.y * m.z_axis.z - m.y_axis.z * m.z_axis.y;
    let cofactor01 = -(m.y_axis.x * m.z_axis.z - m.y_axis.z * m.z_axis.x);
    let cofactor02 = m.y_axis.x * m.z_axis.y - m.y_axis.y * m.z_axis.x;

    let cofactor10 = -(m.x_axis.y * m.z_axis.z - m.x_axis.z * m.z_axis.y);
    let cofactor11 = m.x_axis.x * m.z_axis.z - m.x_axis.z * m.z_axis.x;
    let cofactor12 = -(m.x_axis.x * m.z_axis.y - m.x_axis.y * m.z_axis.x);

    let cofactor20 = m.x_axis.y * m.y_axis.z - m.x_axis.z * m.y_axis.y;
    let cofactor21 = -(m.x_axis.x * m.y_axis.z - m.x_axis.z * m.y_axis.x);
    let cofactor22 = m.x_axis.x * m.y_axis.y - m.x_axis.y * m.y_axis.x;

    let inv_det = 1.0 / det;

    Some(DMat3::from_cols(
        DVec3::new(cofactor00 * inv_det, cofactor10 * inv_det, cofactor20 * inv_det),
        DVec3::new(cofactor01 * inv_det, cofactor11 * inv_det, cofactor21 * inv_det),
        DVec3::new(cofactor02 * inv_det, cofactor12 * inv_det, cofactor22 * inv_det),
    ))
}

/// Compute the determinant of a 3x3 matrix.
pub fn determinant_3x3(m: DMat3) -> f64 {
    m.determinant()
}

// =============================================================================
// Numerical Integration
// =============================================================================

/// Simpson's rule for numerical integration.
///
/// Integrates f from a to b using n subintervals (n must be even)
pub fn simpson_integrate(f: fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let n = if n % 2 != 0 { n + 1 } else { n }; // Ensure even
    let h = (b - a) / n as f64;

    let mut sum = f(a) + f(b);

    for i in 1..n {
        let x = a + i as f64 * h;
        let coef = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += coef * f(x);
    }

    sum * h / 3.0
}

/// Gaussian quadrature nodes and weights for various orders
fn gaussian_nodes_weights(n: usize) -> Vec<(f64, f64)> {
    match n {
        1 => vec![(0.0, 2.0)],
        2 => vec![
            (-FRAC_1_SQRT_2, 1.0),
            (FRAC_1_SQRT_2, 1.0),
        ],
        3 => vec![
            (0.0, 8.0 / 9.0),
            (-0.7745966692414834, 5.0 / 9.0),
            (0.7745966692414834, 5.0 / 9.0),
        ],
        4 => vec![
            (-0.8611363115940526, 0.3478548451374538),
            (-0.3399810435848563, 0.6521451548625461),
            (0.3399810435848563, 0.6521451548625461),
            (0.8611363115940526, 0.3478548451374538),
        ],
        5 => vec![
            (0.0, 0.5688888888888889),
            (-0.5384693101056831, 0.4786286704993665),
            (0.5384693101056831, 0.4786286704993665),
            (-0.9061798459386640, 0.2369268850561891),
            (0.9061798459386640, 0.2369268850561891),
        ],
        6 => vec![
            (-0.9324695142031521, 0.1713244923791704),
            (-0.6612093864662645, 0.3607615730481386),
            (-0.2386191860831969, 0.4679139345726910),
            (0.2386191860831969, 0.4679139345726910),
            (0.6612093864662645, 0.3607615730481386),
            (0.9324695142031521, 0.1713244923791704),
        ],
        _ => gaussian_nodes_weights(6), // Default to 6-point rule
    }
}

/// Gaussian quadrature for numerical integration.
///
/// Integrates f from a to b using n-point Gaussian quadrature
pub fn gaussian_quadrature(f: fn(f64) -> f64, a: f64, b: f64, n_points: usize) -> f64 {
    let nodes_weights = gaussian_nodes_weights(n_points);

    // Transform from [-1, 1] to [a, b]
    let scale = (b - a) / 2.0;
    let shift = (a + b) / 2.0;

    let mut sum = 0.0;
    for (node, weight) in nodes_weights {
        let x = shift + scale * node;
        sum += weight * f(x);
    }

    sum * scale
}

// =============================================================================
// Optimization
// =============================================================================

/// Golden ratio constant
const PHI: f64 = 1.618033988749895;
const RESPHI: f64 = 0.3819660112501051; // 1/phi^2

/// Golden section search for finding minimum.
///
/// # Arguments
/// * `f` - Function to minimize
/// * `a` - Lower bound
/// * `b` - Upper bound
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The x value that minimizes f in [a, b]
pub fn golden_section_min<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, tol: f64) -> f64 {
    let mut lo = a;
    let mut hi = b;

    let mut c = lo + RESPHI * (hi - lo);
    let mut d = hi - RESPHI * (hi - lo);

    let mut fc = f(c);
    let mut fd = f(d);

    while (hi - lo).abs() > tol {
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = lo + RESPHI * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = hi - RESPHI * (hi - lo);
            fd = f(d);
        }
    }

    (lo + hi) / 2.0
}

/// Golden section search for finding maximum.
///
/// # Arguments
/// * `f` - Function to maximize
/// * `a` - Lower bound
/// * `b` - Upper bound
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The x value that maximizes f in [a, b]
pub fn golden_section_max<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, tol: f64) -> f64 {
    golden_section_min(|x| -f(x), a, b, tol)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // --- Root Finding Tests ---

    #[test]
    fn test_newton_raphson_simple() {
        // Solve x^2 - 4 = 0, root is 2
        let f = |x: f64| x * x - 4.0;
        let df = |x: f64| 2.0 * x;

        let root = newton_raphson(f, df, 3.0, 1e-10, 100);
        assert!(root.is_some());
        assert!((root.unwrap() - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_newton_raphson_cubic() {
        // Solve x^3 - x - 2 = 0
        let f = |x: f64| x * x * x - x - 2.0;
        let df = |x: f64| 3.0 * x * x - 1.0;

        let root = newton_raphson(f, df, 2.0, 1e-10, 100);
        assert!(root.is_some());
        assert!((root.unwrap() - 1.521379706804567).abs() < 1e-6);
    }

    #[test]
    fn test_bisection() {
        // Solve x^2 - 4 = 0 in [0, 3]
        let f = |x: f64| x * x - 4.0;

        let root = bisection(f, 0.0, 3.0, 1e-10);
        assert!(root.is_some());
        assert!((root.unwrap() - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_bisection_no_root() {
        // Try to find root of x^2 + 1 = 0 (no real roots)
        let f = |x: f64| x * x + 1.0;

        let root = bisection(f, -2.0, 2.0, 1e-10);
        assert!(root.is_none());
    }

    #[test]
    fn test_secant() {
        // Solve x^2 - 4 = 0
        let f = |x: f64| x * x - 4.0;

        let root = secant(f, 1.0, 3.0, 1e-10);
        assert!(root.is_some());
        assert!((root.unwrap() - 2.0).abs() < 1e-8);
    }

    // --- Multi-dimensional Newton Tests ---

    #[test]
    fn test_newton_2d() {
        // Solve system:
        // x^2 + y^2 = 4
        // x - y = 0
        // Solution: x = y = sqrt(2)
        let f = |v: DVec2| DVec2::new(v.x * v.x + v.y * v.y - 4.0, v.x - v.y);
        let jacobian = |v: DVec2| DMat2::from_cols(DVec2::new(2.0 * v.x, 1.0), DVec2::new(2.0 * v.y, -1.0));

        let root = newton_2d(f, jacobian, DVec2::new(1.5, 1.5), 1e-10);
        assert!(root.is_some());
        let r = root.unwrap();
        assert!((r.x - 2_f64.sqrt()).abs() < 1e-6);
        assert!((r.y - 2_f64.sqrt()).abs() < 1e-6);
    }

    #[test]
    fn test_newton_3d() {
        // Solve system:
        // x + y + z = 6
        // x - y = 0
        // y - z = 0
        // Solution: x = y = z = 2
        let f = |v: DVec3| DVec3::new(v.x + v.y + v.z - 6.0, v.x - v.y, v.y - v.z);
        let jacobian = |_v: DVec3| {
            DMat3::from_cols(
                DVec3::new(1.0, 1.0, 0.0),
                DVec3::new(1.0, -1.0, 1.0),
                DVec3::new(1.0, 0.0, -1.0),
            )
        };

        let root = newton_3d(f, jacobian, DVec3::new(1.0, 1.0, 1.0), 1e-10);
        assert!(root.is_some());
        let r = root.unwrap();
        assert!((r.x - 2.0).abs() < 1e-6);
        assert!((r.y - 2.0).abs() < 1e-6);
        assert!((r.z - 2.0).abs() < 1e-6);
    }

    // --- Polynomial Solver Tests ---

    #[test]
    fn test_solve_linear() {
        let root = solve_linear(2.0, -4.0);
        assert!(root.is_some());
        assert!((root.unwrap() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_solve_quadratic_two_roots() {
        // x^2 - 5x + 6 = 0 has roots 2 and 3
        let roots = solve_quadratic(1.0, -5.0, 6.0);
        assert_eq!(roots.len(), 2);
        assert!((roots[0] - 2.0).abs() < 1e-10);
        assert!((roots[1] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_solve_quadratic_one_root() {
        // x^2 - 4x + 4 = 0 has double root 2
        let roots = solve_quadratic(1.0, -4.0, 4.0);
        assert_eq!(roots.len(), 1);
        assert!((roots[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_solve_quadratic_no_roots() {
        // x^2 + 1 = 0 has no real roots
        let roots = solve_quadratic(1.0, 0.0, 1.0);
        assert!(roots.is_empty());
    }

    #[test]
    fn test_solve_cubic_one_real_root() {
        // x^3 - x - 2 = 0 has one real root
        let roots = solve_cubic(1.0, 0.0, -1.0, -2.0);
        assert_eq!(roots.len(), 1);
        assert!((roots[0] - 1.521379706804567).abs() < 1e-6);
    }

    #[test]
    fn test_solve_cubic_three_real_roots() {
        // x^3 - 6x^2 + 11x - 6 = 0 has roots 1, 2, 3
        let roots = solve_cubic(1.0, -6.0, 11.0, -6.0);
        assert_eq!(roots.len(), 3);
        assert!((roots[0] - 1.0).abs() < 1e-6);
        assert!((roots[1] - 2.0).abs() < 1e-6);
        assert!((roots[2] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_solve_quartic_two_roots() {
        // x^4 - 5x^2 + 4 = 0 has roots -2, -1, 1, 2
        let roots = solve_quartic(1.0, 0.0, -5.0, 0.0, 4.0);
        assert!(roots.len() >= 2);
        // Check that 1 and 2 are among the roots
        assert!(roots.iter().any(|&r| (r - 1.0).abs() < 1e-6));
        assert!(roots.iter().any(|&r| (r - 2.0).abs() < 1e-6));
    }

    #[test]
    fn test_solve_quartic_four_roots() {
        // (x-1)(x-2)(x-3)(x-4) = x^4 - 10x^3 + 35x^2 - 50x + 24
        let roots = solve_quartic(1.0, -10.0, 35.0, -50.0, 24.0);
        assert_eq!(roots.len(), 4);
        assert!((roots[0] - 1.0).abs() < 1e-4);
        assert!((roots[1] - 2.0).abs() < 1e-4);
        assert!((roots[2] - 3.0).abs() < 1e-4);
        assert!((roots[3] - 4.0).abs() < 1e-4);
    }

    // --- Eigenvalue/Matrix Tests ---

    #[test]
    fn test_eigenvalues_2x2() {
        let m = DMat2::from_cols(DVec2::new(2.0, 0.0), DVec2::new(0.0, 3.0));
        let (e1, e2) = eigenvalues_2x2(m);
        assert!((e1 - 3.0).abs() < 1e-10);
        assert!((e2 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_eigenvalues_2x2_rotation() {
        // Rotation by 45 degrees has complex eigenvalues
        let angle = PI / 4.0;
        let m = DMat2::from_cols(
            DVec2::new(angle.cos(), angle.sin()),
            DVec2::new(-angle.sin(), angle.cos()),
        );
        let (e1, e2) = eigenvalues_2x2(m);
        // Real parts should be cos(45deg) ~ 0.707
        assert!((e1 - 0.7071067811865476).abs() < 1e-10);
        assert!((e2 - 0.7071067811865476).abs() < 1e-10);
    }

    #[test]
    fn test_eigenvalues_3x3_diagonal() {
        let m = DMat3::from_cols(
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 2.0, 0.0),
            DVec3::new(0.0, 0.0, 3.0),
        );
        let (e1, e2, e3) = eigenvalues_3x3(m);
        assert!((e1 - 3.0).abs() < 1e-6);
        assert!((e2 - 2.0).abs() < 1e-6);
        assert!((e3 - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_inverse_3x3() {
        let m = DMat3::from_cols(
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(0.0, 1.0, 4.0),
            DVec3::new(5.0, 6.0, 0.0),
        );

        let inv = inverse_3x3(m);
        assert!(inv.is_some());

        let inv = inv.unwrap();
        let product = m * inv;
        let identity = DMat3::IDENTITY;

        for i in 0..3 {
            for j in 0..3 {
                assert!((product.col(i)[j] - identity.col(i)[j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_inverse_3x3_singular() {
        let m = DMat3::from_cols(
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(2.0, 4.0, 6.0), // Linearly dependent
            DVec3::new(5.0, 6.0, 0.0),
        );

        let inv = inverse_3x3(m);
        assert!(inv.is_none());
    }

    #[test]
    fn test_determinant_3x3() {
        let m = DMat3::from_cols(
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(4.0, 5.0, 6.0),
            DVec3::new(7.0, 8.0, 9.0),
        );

        // det of this matrix is 0 (columns are linearly dependent)
        let det = determinant_3x3(m);
        assert!(det.abs() < 1e-10);
    }

    // --- Integration Tests ---

    #[test]
    fn test_simpson_integrate() {
        // Integrate x^2 from 0 to 1, should be 1/3
        let f = |x: f64| x * x;
        let result = simpson_integrate(f, 0.0, 1.0, 100);
        assert!((result - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_simpson_integrate_sin() {
        // Integrate sin(x) from 0 to pi, should be 2
        // Simpson's rule with n=100 has error O(h^4) ≈ 1e-8
        let result = simpson_integrate(f64::sin, 0.0, PI, 100);
        assert!((result - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_gaussian_quadrature() {
        // Integrate x^2 from -1 to 1, should be 2/3
        let f = |x: f64| x * x;
        let result = gaussian_quadrature(f, -1.0, 1.0, 3);
        assert!((result - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_gaussian_quadrature_sin() {
        // Integrate sin(x) from 0 to pi, should be 2
        // 5-point Gaussian quadrature has error ~1e-7 for this case
        let result = gaussian_quadrature(f64::sin, 0.0, PI, 5);
        assert!((result - 2.0).abs() < 1e-6);
    }

    // --- Optimization Tests ---

    #[test]
    fn test_golden_section_min() {
        // Find minimum of (x-2)^2 in [0, 4]
        let f = |x: f64| (x - 2.0) * (x - 2.0);

        let min_x = golden_section_min(f, 0.0, 4.0, 1e-10);
        assert!((min_x - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_golden_section_min_cubic() {
        // Find minimum of x^3 - 3x + 2 in [-2, 2]
        // Minimum is at x = 1
        let f = |x: f64| x * x * x - 3.0 * x + 2.0;

        let min_x = golden_section_min(f, -2.0, 2.0, 1e-10);
        assert!((min_x - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_golden_section_max() {
        // Find maximum of -x^2 + 4 in [-2, 2]
        // Maximum is at x = 0
        let f = |x: f64| -x * x + 4.0;

        let max_x = golden_section_max(f, -2.0, 2.0, 1e-10);
        assert!(max_x.abs() < 1e-6);
    }
}
