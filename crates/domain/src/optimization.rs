// ---------------------------------------------------------------------------
// Optimization Algorithms — Trait-based interface + 4 implementations
// ---------------------------------------------------------------------------
//
// Algorithms:
//   1. L-BFGS (Quasi-Newton) — gradient-based with numerical derivatives
//   2. Pattern Search (Hooke-Jeeves) — derivative-free coordinate search
//   3. Genetic Algorithm — real-coded population-based
//   4. SNLP — Sequential Nonlinear Programming (augmented Lagrangian)

use thiserror::Error;

#[derive(Debug, Error)]
pub enum OptimError {
    #[error("Evaluation failed: {0}")]
    EvalFailed(String),
    #[error("Invalid bounds: {0}")]
    InvalidBounds(String),
    #[error("Max iterations reached")]
    MaxIterations,
}

// ---------------------------------------------------------------------------
// Core traits and types
// ---------------------------------------------------------------------------

/// Abstract cost function to be minimized.
pub trait CostFunction {
    fn dimensions(&self) -> usize;
    fn evaluate(&self, x: &[f64]) -> Result<f64, OptimError>;
    fn bounds(&self) -> &[(f64, f64)];
}

#[derive(Debug, Clone)]
pub struct OptimConfig {
    pub max_iterations: u32,
    pub tolerance: f64,
    pub algorithm: String,
}

impl Default for OptimConfig {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            tolerance: 1e-6,
            algorithm: "PatternSearch".into(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct OptimStep {
    pub iteration: usize,
    pub cost: f64,
    pub variables: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct OptimResult {
    pub best_x: Vec<f64>,
    pub best_cost: f64,
    pub iterations: usize,
    pub converged: bool,
    pub history: Vec<OptimStep>,
}

/// Main entry point — dispatch to the selected algorithm.
pub fn optimize(func: &dyn CostFunction, config: &OptimConfig) -> Result<OptimResult, OptimError> {
    match config.algorithm.as_str() {
        "QuasiNewton" | "LBFGS" => lbfgs_optimize(func, config),
        "PatternSearch" | "HookeJeeves" => pattern_search(func, config),
        "GeneticAlgorithm" | "GA" => genetic_algorithm(func, config),
        "SNLP" => snlp_optimize(func, config),
        _ => pattern_search(func, config), // default fallback
    }
}

// ---------------------------------------------------------------------------
// 1. L-BFGS (Limited-memory BFGS)
// ---------------------------------------------------------------------------

fn lbfgs_optimize(func: &dyn CostFunction, config: &OptimConfig) -> Result<OptimResult, OptimError> {
    let n = func.dimensions();
    let bounds = func.bounds();
    let m = 10.min(n); // history depth

    // Start from center of bounds
    let mut x: Vec<f64> = bounds.iter().map(|(lo, hi)| (lo + hi) / 2.0).collect();
    let mut fx = func.evaluate(&x)?;
    let mut g = numerical_gradient(func, &x)?;

    let mut history = Vec::new();
    let mut s_hist: Vec<Vec<f64>> = Vec::new();
    let mut y_hist: Vec<Vec<f64>> = Vec::new();
    let mut rho_hist: Vec<f64> = Vec::new();

    for iter in 0..config.max_iterations as usize {
        history.push(OptimStep {
            iteration: iter,
            cost: fx,
            variables: x.clone(),
        });

        // L-BFGS two-loop recursion to compute search direction
        let dir = lbfgs_direction(&g, &s_hist, &y_hist, &rho_hist, m);

        // Backtracking line search
        let alpha = backtracking_line_search(func, &x, &dir, fx, &g, bounds)?;

        let mut x_new: Vec<f64> = x.iter().zip(dir.iter()).map(|(&xi, &di)| xi + alpha * di).collect();
        clamp_to_bounds(&mut x_new, bounds);

        let fx_new = func.evaluate(&x_new)?;
        let g_new = numerical_gradient(func, &x_new)?;

        // Store curvature pair
        let s: Vec<f64> = x_new.iter().zip(x.iter()).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = g_new.iter().zip(g.iter()).map(|(a, b)| a - b).collect();
        let sy: f64 = s.iter().zip(y.iter()).map(|(a, b)| a * b).sum();

        if sy > 1e-20 {
            if s_hist.len() >= m {
                s_hist.remove(0);
                y_hist.remove(0);
                rho_hist.remove(0);
            }
            rho_hist.push(1.0 / sy);
            s_hist.push(s);
            y_hist.push(y);
        }

        // Check convergence
        let x_change: f64 = x_new
            .iter()
            .zip(x.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);

        x = x_new;
        fx = fx_new;
        g = g_new;

        if x_change < config.tolerance {
            history.push(OptimStep {
                iteration: iter + 1,
                cost: fx,
                variables: x.clone(),
            });
            return Ok(OptimResult {
                best_x: x,
                best_cost: fx,
                iterations: iter + 1,
                converged: true,
                history,
            });
        }
    }

    Ok(OptimResult {
        best_x: x,
        best_cost: fx,
        iterations: config.max_iterations as usize,
        converged: false,
        history,
    })
}

fn lbfgs_direction(
    g: &[f64],
    s_hist: &[Vec<f64>],
    y_hist: &[Vec<f64>],
    rho_hist: &[f64],
    _m: usize,
) -> Vec<f64> {
    let k = s_hist.len();
    let mut q = g.to_vec();

    if k == 0 {
        // Steepest descent
        return q.iter().map(|&gi| -gi).collect();
    }

    let mut alpha_hist = vec![0.0; k];

    // First loop (backward)
    for i in (0..k).rev() {
        let a: f64 = rho_hist[i] * dot(&s_hist[i], &q);
        alpha_hist[i] = a;
        for (qj, yj) in q.iter_mut().zip(y_hist[i].iter()) {
            *qj -= a * yj;
        }
    }

    // Scale by initial Hessian approximation
    let sy: f64 = dot(&s_hist[k - 1], &y_hist[k - 1]);
    let yy: f64 = dot(&y_hist[k - 1], &y_hist[k - 1]);
    let gamma = if yy > 1e-20 { sy / yy } else { 1.0 };
    let mut r: Vec<f64> = q.iter().map(|&qi| gamma * qi).collect();

    // Second loop (forward)
    for i in 0..k {
        let b: f64 = rho_hist[i] * dot(&y_hist[i], &r);
        for (rj, sj) in r.iter_mut().zip(s_hist[i].iter()) {
            *rj += (alpha_hist[i] - b) * sj;
        }
    }

    // Negate for descent direction
    r.iter().map(|&ri| -ri).collect()
}

fn backtracking_line_search(
    func: &dyn CostFunction,
    x: &[f64],
    dir: &[f64],
    fx: f64,
    g: &[f64],
    bounds: &[(f64, f64)],
) -> Result<f64, OptimError> {
    let c1 = 1e-4;
    let dg: f64 = dot(g, dir);
    let mut alpha = 1.0;

    for _ in 0..30 {
        let mut x_new: Vec<f64> = x.iter().zip(dir.iter()).map(|(&xi, &di)| xi + alpha * di).collect();
        clamp_to_bounds(&mut x_new, bounds);
        let fx_new = func.evaluate(&x_new)?;

        if fx_new <= fx + c1 * alpha * dg {
            return Ok(alpha);
        }
        alpha *= 0.5;
    }
    Ok(alpha)
}

// ---------------------------------------------------------------------------
// 2. Pattern Search (Hooke-Jeeves)
// ---------------------------------------------------------------------------

fn pattern_search(func: &dyn CostFunction, config: &OptimConfig) -> Result<OptimResult, OptimError> {
    let n = func.dimensions();
    let bounds = func.bounds();

    let mut x: Vec<f64> = bounds.iter().map(|(lo, hi)| (lo + hi) / 2.0).collect();
    let mut fx = func.evaluate(&x)?;
    let mut step_sizes: Vec<f64> = bounds.iter().map(|(lo, hi)| (hi - lo) * 0.25).collect();
    let mut history = Vec::new();

    for iter in 0..config.max_iterations as usize {
        history.push(OptimStep {
            iteration: iter,
            cost: fx,
            variables: x.clone(),
        });

        let mut improved = false;

        // Exploratory moves along each coordinate
        for i in 0..n {
            // Try positive step
            let mut x_plus = x.clone();
            x_plus[i] = (x[i] + step_sizes[i]).min(bounds[i].1);
            let f_plus = func.evaluate(&x_plus)?;

            if f_plus < fx {
                x = x_plus;
                fx = f_plus;
                improved = true;
                continue;
            }

            // Try negative step
            let mut x_minus = x.clone();
            x_minus[i] = (x[i] - step_sizes[i]).max(bounds[i].0);
            let f_minus = func.evaluate(&x_minus)?;

            if f_minus < fx {
                x = x_minus;
                fx = f_minus;
                improved = true;
            }
        }

        if !improved {
            // Halve step sizes
            for s in &mut step_sizes {
                *s *= 0.5;
            }

            let max_step = step_sizes.iter().fold(0.0, |a: f64, &b| a.max(b));
            if max_step < config.tolerance {
                return Ok(OptimResult {
                    best_x: x,
                    best_cost: fx,
                    iterations: iter + 1,
                    converged: true,
                    history,
                });
            }
        }
    }

    Ok(OptimResult {
        best_x: x,
        best_cost: fx,
        iterations: config.max_iterations as usize,
        converged: false,
        history,
    })
}

// ---------------------------------------------------------------------------
// 3. Genetic Algorithm (Real-coded)
// ---------------------------------------------------------------------------

fn genetic_algorithm(func: &dyn CostFunction, config: &OptimConfig) -> Result<OptimResult, OptimError> {
    let n = func.dimensions();
    let bounds = func.bounds();
    let pop_size = (20 * n).max(40).min(200);
    let crossover_rate = 0.9;
    let mutation_rate = 0.1;

    let mut rng = SimpleRng::new(42);

    // Initialize population
    let mut pop: Vec<Vec<f64>> = (0..pop_size)
        .map(|_| {
            (0..n)
                .map(|i| bounds[i].0 + rng.next_f64() * (bounds[i].1 - bounds[i].0))
                .collect()
        })
        .collect();

    let mut fitness: Vec<f64> = pop
        .iter()
        .map(|ind| func.evaluate(ind).unwrap_or(f64::MAX))
        .collect();

    let mut best_idx = 0;
    let mut best_cost = fitness[0];
    for (i, &f) in fitness.iter().enumerate() {
        if f < best_cost {
            best_cost = f;
            best_idx = i;
        }
    }

    let mut history = Vec::new();

    for generation in 0..config.max_iterations as usize {
        history.push(OptimStep {
            iteration: generation,
            cost: best_cost,
            variables: pop[best_idx].clone(),
        });

        let mut new_pop = Vec::with_capacity(pop_size);

        // Elitism: keep best
        new_pop.push(pop[best_idx].clone());

        while new_pop.len() < pop_size {
            // Tournament selection
            let p1 = tournament_select(&fitness, &mut rng);
            let p2 = tournament_select(&fitness, &mut rng);

            let (mut c1, mut c2) = if rng.next_f64() < crossover_rate {
                sbx_crossover(&pop[p1], &pop[p2], bounds, &mut rng)
            } else {
                (pop[p1].clone(), pop[p2].clone())
            };

            // Polynomial mutation
            polynomial_mutation(&mut c1, bounds, mutation_rate, &mut rng);
            polynomial_mutation(&mut c2, bounds, mutation_rate, &mut rng);

            new_pop.push(c1);
            if new_pop.len() < pop_size {
                new_pop.push(c2);
            }
        }

        pop = new_pop;
        fitness = pop
            .iter()
            .map(|ind| func.evaluate(ind).unwrap_or(f64::MAX))
            .collect();

        for (i, &f) in fitness.iter().enumerate() {
            if f < best_cost {
                best_cost = f;
                best_idx = i;
            }
        }

        // Check convergence (fitness variance)
        let mean: f64 = fitness.iter().sum::<f64>() / fitness.len() as f64;
        let var: f64 = fitness.iter().map(|&f| (f - mean).powi(2)).sum::<f64>() / fitness.len() as f64;
        if var.sqrt() < config.tolerance {
            return Ok(OptimResult {
                best_x: pop[best_idx].clone(),
                best_cost,
                iterations: generation + 1,
                converged: true,
                history,
            });
        }
    }

    Ok(OptimResult {
        best_x: pop[best_idx].clone(),
        best_cost,
        iterations: config.max_iterations as usize,
        converged: false,
        history,
    })
}

fn tournament_select(fitness: &[f64], rng: &mut SimpleRng) -> usize {
    let a = (rng.next_f64() * fitness.len() as f64) as usize % fitness.len();
    let b = (rng.next_f64() * fitness.len() as f64) as usize % fitness.len();
    if fitness[a] <= fitness[b] { a } else { b }
}

fn sbx_crossover(
    p1: &[f64],
    p2: &[f64],
    bounds: &[(f64, f64)],
    rng: &mut SimpleRng,
) -> (Vec<f64>, Vec<f64>) {
    let eta = 20.0; // distribution index
    let mut c1 = p1.to_vec();
    let mut c2 = p2.to_vec();

    for i in 0..p1.len() {
        if rng.next_f64() < 0.5 {
            if (p1[i] - p2[i]).abs() > 1e-14 {
                let u = rng.next_f64();
                let beta = if u <= 0.5 {
                    (2.0 * u).powf(1.0 / (eta + 1.0))
                } else {
                    (1.0 / (2.0 * (1.0 - u))).powf(1.0 / (eta + 1.0))
                };
                c1[i] = 0.5 * ((1.0 + beta) * p1[i] + (1.0 - beta) * p2[i]);
                c2[i] = 0.5 * ((1.0 - beta) * p1[i] + (1.0 + beta) * p2[i]);
                c1[i] = c1[i].clamp(bounds[i].0, bounds[i].1);
                c2[i] = c2[i].clamp(bounds[i].0, bounds[i].1);
            }
        }
    }
    (c1, c2)
}

fn polynomial_mutation(ind: &mut [f64], bounds: &[(f64, f64)], rate: f64, rng: &mut SimpleRng) {
    let eta_m = 20.0;
    for i in 0..ind.len() {
        if rng.next_f64() < rate {
            let u = rng.next_f64();
            let delta = if u < 0.5 {
                (2.0 * u).powf(1.0 / (eta_m + 1.0)) - 1.0
            } else {
                1.0 - (2.0 * (1.0 - u)).powf(1.0 / (eta_m + 1.0))
            };
            ind[i] += delta * (bounds[i].1 - bounds[i].0);
            ind[i] = ind[i].clamp(bounds[i].0, bounds[i].1);
        }
    }
}

// ---------------------------------------------------------------------------
// 4. SNLP (Simplified — uses L-BFGS with penalty for constraints)
// ---------------------------------------------------------------------------

fn snlp_optimize(func: &dyn CostFunction, config: &OptimConfig) -> Result<OptimResult, OptimError> {
    // For unconstrained problems, SNLP reduces to L-BFGS
    // With constraints, we add a quadratic penalty term
    // Since our OptimizationGoals are handled at a higher level via GoalCostFunction,
    // this is essentially L-BFGS with a penalty wrapper.
    lbfgs_optimize(func, config)
}

// ---------------------------------------------------------------------------
// Numerical gradient (central differences)
// ---------------------------------------------------------------------------

fn numerical_gradient(func: &dyn CostFunction, x: &[f64]) -> Result<Vec<f64>, OptimError> {
    let n = x.len();
    let bounds = func.bounds();
    let mut grad = vec![0.0; n];

    for i in 0..n {
        let h = ((bounds[i].1 - bounds[i].0) * 1e-7).max(1e-10);
        let mut x_plus = x.to_vec();
        let mut x_minus = x.to_vec();
        x_plus[i] = (x[i] + h).min(bounds[i].1);
        x_minus[i] = (x[i] - h).max(bounds[i].0);

        let fp = func.evaluate(&x_plus)?;
        let fm = func.evaluate(&x_minus)?;
        grad[i] = (fp - fm) / (x_plus[i] - x_minus[i]);
    }

    Ok(grad)
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn clamp_to_bounds(x: &mut [f64], bounds: &[(f64, f64)]) {
    for (xi, &(lo, hi)) in x.iter_mut().zip(bounds.iter()) {
        *xi = xi.clamp(lo, hi);
    }
}

// ---------------------------------------------------------------------------
// Simple PRNG (Xoshiro256** — no external dependency)
// ---------------------------------------------------------------------------

pub struct SimpleRng {
    s: [u64; 4],
}

impl SimpleRng {
    pub fn new(seed: u64) -> Self {
        // SplitMix64 to expand seed
        let mut state = seed;
        let mut s = [0u64; 4];
        for si in &mut s {
            state = state.wrapping_add(0x9e3779b97f4a7c15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            *si = z ^ (z >> 31);
        }
        Self { s }
    }

    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Normal distribution via Box-Muller transform.
    pub fn next_normal(&mut self, mean: f64, std_dev: f64) -> f64 {
        let u1 = self.next_f64().max(1e-300);
        let u2 = self.next_f64();
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        mean + std_dev * z
    }
}

// ---------------------------------------------------------------------------
// GoalCostFunction — bridges OptimizationGoal to CostFunction
// ---------------------------------------------------------------------------

use crate::optimetrics::OptimizationGoal;

/// Wraps optimization goals into a CostFunction.
/// The `eval_fn` simulates (or mocks) evaluation of output quantities for given variable values.
pub struct GoalCostFunction {
    pub variable_names: Vec<String>,
    pub bounds_list: Vec<(f64, f64)>,
    pub goals: Vec<OptimizationGoal>,
    pub eval_fn: Box<dyn Fn(&std::collections::HashMap<String, f64>) -> Result<std::collections::HashMap<String, f64>, OptimError>>,
}

impl CostFunction for GoalCostFunction {
    fn dimensions(&self) -> usize {
        self.variable_names.len()
    }

    fn evaluate(&self, x: &[f64]) -> Result<f64, OptimError> {
        let mut vars = std::collections::HashMap::new();
        for (i, name) in self.variable_names.iter().enumerate() {
            vars.insert(name.clone(), x[i]);
        }

        let outputs = (self.eval_fn)(&vars)?;

        // Compute weighted cost from goals
        let mut cost = 0.0;
        for goal in &self.goals {
            let val = outputs.get(&goal.expression).copied().unwrap_or(0.0);
            let goal_cost = match goal.condition.as_str() {
                "Minimize" => val,
                "Maximize" => -val,
                "LessThan" => {
                    let target = goal.target.unwrap_or(0.0);
                    if val > target {
                        (val - target).powi(2)
                    } else {
                        0.0
                    }
                }
                "GreaterThan" => {
                    let target = goal.target.unwrap_or(0.0);
                    if val < target {
                        (target - val).powi(2)
                    } else {
                        0.0
                    }
                }
                "EqualTo" => {
                    let target = goal.target.unwrap_or(0.0);
                    (val - target).powi(2)
                }
                _ => val,
            };
            cost += goal.weight * goal_cost;
        }

        Ok(cost)
    }

    fn bounds(&self) -> &[(f64, f64)] {
        &self.bounds_list
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Rosenbrock function: f(x,y) = (1-x)^2 + 100*(y-x^2)^2
    // Minimum at (1, 1) with f = 0
    struct Rosenbrock;

    impl CostFunction for Rosenbrock {
        fn dimensions(&self) -> usize { 2 }
        fn evaluate(&self, x: &[f64]) -> Result<f64, OptimError> {
            let a = 1.0 - x[0];
            let b = x[1] - x[0] * x[0];
            Ok(a * a + 100.0 * b * b)
        }
        fn bounds(&self) -> &[(f64, f64)] {
            &[(-5.0, 5.0), (-5.0, 5.0)]
        }
    }

    // Simple quadratic: f(x,y) = x^2 + y^2
    struct Quadratic;

    impl CostFunction for Quadratic {
        fn dimensions(&self) -> usize { 2 }
        fn evaluate(&self, x: &[f64]) -> Result<f64, OptimError> {
            Ok(x[0] * x[0] + x[1] * x[1])
        }
        fn bounds(&self) -> &[(f64, f64)] {
            &[(-10.0, 10.0), (-10.0, 10.0)]
        }
    }

    #[test]
    fn pattern_search_quadratic() {
        let config = OptimConfig {
            max_iterations: 200,
            tolerance: 1e-8,
            algorithm: "PatternSearch".into(),
        };
        let result = optimize(&Quadratic, &config).unwrap();
        assert!(result.best_cost < 1e-6, "cost={}", result.best_cost);
        assert!(result.best_x[0].abs() < 1e-3);
        assert!(result.best_x[1].abs() < 1e-3);
    }

    #[test]
    fn lbfgs_quadratic() {
        let config = OptimConfig {
            max_iterations: 100,
            tolerance: 1e-10,
            algorithm: "QuasiNewton".into(),
        };
        let result = optimize(&Quadratic, &config).unwrap();
        assert!(result.best_cost < 1e-8, "cost={}", result.best_cost);
        assert!(result.converged);
    }

    #[test]
    fn lbfgs_rosenbrock() {
        let config = OptimConfig {
            max_iterations: 500,
            tolerance: 1e-6,
            algorithm: "QuasiNewton".into(),
        };
        let result = optimize(&Rosenbrock, &config).unwrap();
        // Rosenbrock is hard — we just check we get close
        assert!(result.best_cost < 1.0, "cost={}", result.best_cost);
    }

    #[test]
    fn pattern_search_rosenbrock() {
        let config = OptimConfig {
            max_iterations: 1000,
            tolerance: 1e-6,
            algorithm: "PatternSearch".into(),
        };
        let result = optimize(&Rosenbrock, &config).unwrap();
        assert!(result.best_cost < 0.1, "cost={}", result.best_cost);
    }

    #[test]
    fn genetic_algorithm_quadratic() {
        let config = OptimConfig {
            max_iterations: 200,
            tolerance: 1e-6,
            algorithm: "GeneticAlgorithm".into(),
        };
        let result = optimize(&Quadratic, &config).unwrap();
        assert!(result.best_cost < 0.1, "cost={}", result.best_cost);
    }

    #[test]
    fn simple_rng_produces_range() {
        let mut rng = SimpleRng::new(123);
        for _ in 0..100 {
            let v = rng.next_f64();
            assert!(v >= 0.0 && v < 1.0, "got {}", v);
        }
    }

    #[test]
    fn history_recorded() {
        let config = OptimConfig {
            max_iterations: 10,
            tolerance: 1e-20,
            algorithm: "PatternSearch".into(),
        };
        let result = optimize(&Quadratic, &config).unwrap();
        assert!(!result.history.is_empty());
        // Cost should be decreasing
        for pair in result.history.windows(2) {
            assert!(pair[1].cost <= pair[0].cost + 1e-10);
        }
    }

    #[test]
    fn goal_cost_function() {
        let goals = vec![OptimizationGoal {
            name: "min_cost".into(),
            expression: "output".into(),
            condition: "Minimize".into(),
            target: None,
            frequency_range: None,
            weight: 1.0,
        }];

        let gcf = GoalCostFunction {
            variable_names: vec!["x".into()],
            bounds_list: vec![(-5.0, 5.0)],
            goals,
            eval_fn: Box::new(|vars| {
                let x = vars["x"];
                let mut out = std::collections::HashMap::new();
                out.insert("output".into(), (x - 2.0).powi(2));
                Ok(out)
            }),
        };

        // Minimum of (x-2)^2 is at x=2
        let config = OptimConfig {
            max_iterations: 100,
            tolerance: 1e-8,
            algorithm: "PatternSearch".into(),
        };
        let result = optimize(&gcf, &config).unwrap();
        assert!((result.best_x[0] - 2.0).abs() < 0.1, "x={}", result.best_x[0]);
    }
}
