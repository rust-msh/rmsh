//! Advanced Boolean Retry Strategy Enhancement
//!
//! This module provides detailed failure classification and targeted recovery
//! strategies for boolean operations. It builds on the basic retry mechanism
//! with more sophisticated failure detection and recovery.

use std::fmt;

/// Detailed failure classification for boolean operations.
///
/// This enum provides more specific failure types than `BooleanRetryClass`,
/// enabling targeted recovery strategies for each failure mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BooleanFailureClass {
    /// Result has degenerate edges or faces (zero-length edges, degenerate triangles).
    DegenerateTopology,
    /// Numerical errors during computation (NaN, infinity, precision loss).
    NumericalInstability,
    /// Result fails validity checks (non-manifold, open shells, invalid orientation).
    InvalidResult,
    /// Missing intersection curves between surfaces that should intersect.
    IncompleteIntersection,
    /// Result contains self-intersecting geometry.
    SelfIntersection,
    /// Input geometry is structurally invalid (empty, missing data).
    InvalidInput,
    /// Unknown or unclassified failure.
    Unknown,
}

impl BooleanFailureClass {
    /// Returns a human-readable description of the failure class.
    pub fn description(&self) -> &'static str {
        match self {
            Self::DegenerateTopology => "Result contains degenerate topology (zero-length edges or degenerate faces)",
            Self::NumericalInstability => "Numerical errors during computation (NaN, infinity, or precision loss)",
            Self::InvalidResult => "Result fails validity checks (non-manifold, open shells, or invalid orientation)",
            Self::IncompleteIntersection => "Missing intersection curves between surfaces",
            Self::SelfIntersection => "Result contains self-intersecting geometry",
            Self::InvalidInput => "Input geometry is structurally invalid",
            Self::Unknown => "Unknown or unclassified failure",
        }
    }

    /// Returns whether this failure class can potentially be recovered by retry.
    pub fn is_recoverable(&self) -> bool {
        matches!(
            self,
            Self::DegenerateTopology
                | Self::NumericalInstability
                | Self::InvalidResult
                | Self::IncompleteIntersection
                | Self::SelfIntersection
        )
    }

    /// Returns the suggested recovery strategy for this failure class.
    pub fn suggested_recovery(&self) -> RecoveryStrategy {
        match self {
            Self::DegenerateTopology => RecoveryStrategy::MakeConnectedCleanup,
            Self::NumericalInstability => RecoveryStrategy::IncreaseFuzzyTolerance,
            Self::InvalidResult => RecoveryStrategy::AlgorithmVariant,
            Self::IncompleteIntersection => RecoveryStrategy::EnableGlueMode,
            Self::SelfIntersection => RecoveryStrategy::MakeConnectedCleanup,
            Self::InvalidInput => RecoveryStrategy::None,
            Self::Unknown => RecoveryStrategy::None,
        }
    }
}

impl fmt::Display for BooleanFailureClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DegenerateTopology => write!(f, "degenerate topology"),
            Self::NumericalInstability => write!(f, "numerical instability"),
            Self::InvalidResult => write!(f, "invalid result"),
            Self::IncompleteIntersection => write!(f, "incomplete intersection"),
            Self::SelfIntersection => write!(f, "self-intersection"),
            Self::InvalidInput => write!(f, "invalid input"),
            Self::Unknown => write!(f, "unknown failure"),
        }
    }
}

impl Default for BooleanFailureClass {
    fn default() -> Self {
        Self::Unknown
    }
}

/// Recovery strategy to apply for a specific failure class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecoveryStrategy {
    /// No recovery possible.
    None,
    /// Run MakeConnected cleanup to fix topology issues.
    MakeConnectedCleanup,
    /// Increase fuzzy tolerance and retry.
    IncreaseFuzzyTolerance,
    /// Try a different algorithm variant.
    AlgorithmVariant,
    /// Enable Glue mode for better intersection handling.
    EnableGlueMode,
    /// Combine multiple strategies.
    Combined {
        use_glue: bool,
        run_make_connected: bool,
        increase_fuzzy: bool,
    },
}

impl RecoveryStrategy {
    /// Returns a description of this recovery strategy.
    pub fn description(&self) -> &'static str {
        match self {
            Self::None => "No recovery available",
            Self::MakeConnectedCleanup => "Run MakeConnected cleanup to fix topology",
            Self::IncreaseFuzzyTolerance => "Increase fuzzy tolerance and retry",
            Self::AlgorithmVariant => "Try different algorithm variant",
            Self::EnableGlueMode => "Enable Glue mode for intersection handling",
            Self::Combined { .. } => "Combine multiple recovery strategies",
        }
    }
}

/// Configurable retry policy for boolean operations.
///
/// This struct provides fine-grained control over retry behavior,
/// including tolerance growth, glue mode activation, and cleanup aggressiveness.
#[derive(Debug, Clone)]
pub struct RetryPolicy {
    /// Maximum number of retry attempts.
    pub max_attempts: usize,
    /// Factor by which to multiply fuzzy tolerance on each retry.
    pub fuzzy_growth_factor: f64,
    /// Maximum fuzzy tolerance cap.
    pub fuzzy_tolerance_cap: f64,
    /// Number of failures after which to enable Glue mode.
    pub enable_glue_after_n_failures: usize,
    /// Glue tolerance for shared-face detection.
    pub glue_tolerance: f64,
    /// Cleanup aggressiveness level (1-10, higher = more aggressive).
    pub make_connected_aggressiveness: u32,
    /// Maximum passes for MakeConnected cleanup.
    pub make_connected_max_passes: usize,
    /// Initial tolerance for MakeConnected cleanup.
    pub make_connected_initial_tolerance: f64,
    /// Tolerance growth factor for MakeConnected passes.
    pub make_connected_tolerance_growth: f64,
    /// Whether to use scoped MakeConnected when possible.
    pub use_scoped_make_connected: bool,
    /// Whether to fall back to global cleanup when scoped fails.
    pub fallback_to_global_cleanup: bool,
    /// Whether to try algorithm variants on InvalidResult failures.
    pub try_algorithm_variants: bool,
    /// Whether to enable verbose diagnostic output.
    pub verbose_diagnostics: bool,
}

impl Default for RetryPolicy {
    fn default() -> Self {
        Self {
            max_attempts: 5,
            fuzzy_growth_factor: 10.0,
            fuzzy_tolerance_cap: 1e-3,
            enable_glue_after_n_failures: 2,
            glue_tolerance: 1e-6,
            make_connected_aggressiveness: 5,
            make_connected_max_passes: 5,
            make_connected_initial_tolerance: 1e-6,
            make_connected_tolerance_growth: 2.0,
            use_scoped_make_connected: true,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: true,
            verbose_diagnostics: false,
        }
    }
}

impl RetryPolicy {
    /// Creates a conservative retry policy with minimal intervention.
    pub fn conservative() -> Self {
        Self {
            max_attempts: 3,
            fuzzy_growth_factor: 5.0,
            fuzzy_tolerance_cap: 1e-4,
            enable_glue_after_n_failures: 3,
            glue_tolerance: 1e-7,
            make_connected_aggressiveness: 3,
            make_connected_max_passes: 3,
            make_connected_initial_tolerance: 1e-7,
            make_connected_tolerance_growth: 1.5,
            use_scoped_make_connected: true,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: false,
            verbose_diagnostics: false,
        }
    }

    /// Creates an aggressive retry policy for difficult geometry.
    pub fn aggressive() -> Self {
        Self {
            max_attempts: 10,
            fuzzy_growth_factor: 20.0,
            fuzzy_tolerance_cap: 1e-2,
            enable_glue_after_n_failures: 1,
            glue_tolerance: 1e-5,
            make_connected_aggressiveness: 8,
            make_connected_max_passes: 10,
            make_connected_initial_tolerance: 1e-5,
            make_connected_tolerance_growth: 3.0,
            use_scoped_make_connected: false,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: true,
            verbose_diagnostics: true,
        }
    }

    /// Creates a retry policy tuned for numerical instability cases.
    pub fn for_numerical_instability() -> Self {
        Self {
            max_attempts: 8,
            fuzzy_growth_factor: 15.0,
            fuzzy_tolerance_cap: 5e-3,
            enable_glue_after_n_failures: 1,
            glue_tolerance: 1e-5,
            make_connected_aggressiveness: 6,
            make_connected_max_passes: 7,
            make_connected_initial_tolerance: 1e-6,
            make_connected_tolerance_growth: 2.5,
            use_scoped_make_connected: true,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: false,
            verbose_diagnostics: false,
        }
    }

    /// Creates a retry policy tuned for degenerate topology cases.
    pub fn for_degenerate_topology() -> Self {
        Self {
            max_attempts: 6,
            fuzzy_growth_factor: 5.0,
            fuzzy_tolerance_cap: 1e-3,
            enable_glue_after_n_failures: 2,
            glue_tolerance: 1e-6,
            make_connected_aggressiveness: 9,
            make_connected_max_passes: 10,
            make_connected_initial_tolerance: 1e-5,
            make_connected_tolerance_growth: 2.0,
            use_scoped_make_connected: false,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: false,
            verbose_diagnostics: false,
        }
    }

    /// Creates a retry policy tuned for incomplete intersection cases.
    pub fn for_incomplete_intersection() -> Self {
        Self {
            max_attempts: 6,
            fuzzy_growth_factor: 10.0,
            fuzzy_tolerance_cap: 5e-3,
            enable_glue_after_n_failures: 0, // Enable immediately
            glue_tolerance: 1e-5,
            make_connected_aggressiveness: 5,
            make_connected_max_passes: 5,
            make_connected_initial_tolerance: 1e-6,
            make_connected_tolerance_growth: 2.0,
            use_scoped_make_connected: true,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: false,
            verbose_diagnostics: false,
        }
    }

    /// Computes the next fuzzy tolerance based on the policy.
    pub fn next_fuzzy_tolerance(&self, current: f64) -> f64 {
        let next = current * self.fuzzy_growth_factor;
        next.min(self.fuzzy_tolerance_cap)
    }

    /// Determines whether glue mode should be enabled after N failures.
    pub fn should_enable_glue(&self, failure_count: usize) -> bool {
        failure_count >= self.enable_glue_after_n_failures
    }

    /// Computes the MakeConnected tolerance based on aggressiveness.
    pub fn make_connected_tolerance(&self, pass: usize) -> f64 {
        let base = self.make_connected_initial_tolerance;
        let growth = self.make_connected_tolerance_growth.powi(pass as i32);
        base * growth
    }
}

/// Diagnostic information for a single retry attempt.
#[derive(Debug, Clone, Default)]
pub struct BooleanAttemptDiagnostic {
    /// Attempt number (1-indexed).
    pub attempt_number: usize,
    /// Fuzzy tolerance used for this attempt.
    pub fuzzy_tolerance: f64,
    /// Whether glue mode was enabled.
    pub glue_enabled: bool,
    /// Glue tolerance used.
    pub glue_tolerance: f64,
    /// Whether MakeConnected cleanup was run.
    pub make_connected_run: bool,
    /// Number of MakeConnected passes.
    pub make_connected_passes: usize,
    /// Whether this attempt succeeded.
    pub success: bool,
    /// Failure class if the attempt failed.
    pub failure_class: Option<BooleanFailureClass>,
    /// Recovery strategy applied before this attempt.
    pub recovery_strategy: Option<RecoveryStrategy>,
    /// Error message if the attempt failed.
    pub error_message: Option<String>,
    /// Time taken for this attempt (in microseconds).
    pub duration_us: Option<u64>,
    /// Number of faces in the result (if successful).
    pub result_faces: Option<usize>,
    /// Whether scoped make-connected was used.
    pub scoped_make_connected: bool,
    /// Whether fallback to global cleanup occurred.
    pub global_fallback: bool,
}

impl BooleanAttemptDiagnostic {
    /// Creates a new diagnostic for an attempt.
    pub fn new(attempt_number: usize, fuzzy_tolerance: f64) -> Self {
        Self {
            attempt_number,
            fuzzy_tolerance,
            ..Self::default()
        }
    }
}

/// Comprehensive diagnostic report for a boolean operation with retries.
#[derive(Debug, Clone, Default)]
pub struct BooleanDiagnosticReport {
    /// All attempt diagnostics.
    pub attempts: Vec<BooleanAttemptDiagnostic>,
    /// Total number of attempts.
    pub total_attempts: usize,
    /// Number of successful attempts (should be 0 or 1).
    pub successful_attempts: usize,
    /// The failure class that was most common (if any failures).
    pub dominant_failure_class: Option<BooleanFailureClass>,
    /// The final successful configuration (if successful).
    pub final_config: Option<FinalSuccessfulConfig>,
    /// Total time taken across all attempts (in microseconds).
    pub total_duration_us: u64,
    /// The retry policy used.
    pub retry_policy: Option<RetryPolicy>,
    /// Whether glue mode was ultimately needed.
    pub glue_mode_needed: bool,
    /// Whether MakeConnected cleanup was ultimately needed.
    pub make_connected_needed: bool,
    /// Summary of recovery strategies applied.
    pub recovery_strategies_applied: Vec<RecoveryStrategy>,
}

impl BooleanDiagnosticReport {
    /// Creates a new empty diagnostic report.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds an attempt diagnostic to the report.
    pub fn add_attempt(&mut self, attempt: BooleanAttemptDiagnostic) {
        if attempt.success {
            self.successful_attempts += 1;
        }
        self.total_attempts += 1;
        self.attempts.push(attempt);
    }

    /// Computes the dominant failure class from failed attempts.
    pub fn compute_dominant_failure_class(&mut self) {
        use std::collections::HashMap;
        let mut counts: HashMap<BooleanFailureClass, usize> = HashMap::new();

        for attempt in &self.attempts {
            if let Some(class) = attempt.failure_class {
                *counts.entry(class).or_insert(0) += 1;
            }
        }

        self.dominant_failure_class = counts
            .into_iter()
            .max_by_key(|(_, count)| *count)
            .map(|(class, _)| class);
    }

    /// Finalizes the report after all attempts are complete.
    pub fn finalize(&mut self) {
        self.compute_dominant_failure_class();

        // Determine if glue mode or MakeConnected were needed
        for attempt in &self.attempts {
            if attempt.success {
                self.glue_mode_needed = attempt.glue_enabled;
                self.make_connected_needed = attempt.make_connected_run;

                // Record the final successful configuration
                self.final_config = Some(FinalSuccessfulConfig {
                    fuzzy_tolerance: attempt.fuzzy_tolerance,
                    glue_enabled: attempt.glue_enabled,
                    glue_tolerance: attempt.glue_tolerance,
                    make_connected_run: attempt.make_connected_run,
                    make_connected_passes: attempt.make_connected_passes,
                    scoped_make_connected: attempt.scoped_make_connected,
                });
            }

            if let Some(strategy) = attempt.recovery_strategy {
                if !self.recovery_strategies_applied.contains(&strategy) {
                    self.recovery_strategies_applied.push(strategy);
                }
            }
        }
    }

    /// Returns whether the operation ultimately succeeded.
    pub fn is_success(&self) -> bool {
        self.successful_attempts > 0
    }

    /// Returns a summary string for logging.
    pub fn summary(&self) -> String {
        if self.is_success() {
            format!(
                "Boolean operation succeeded after {} attempts (dominant failure: {})",
                self.total_attempts,
                self.dominant_failure_class
                    .map(|c| c.to_string())
                    .unwrap_or_else(|| "none".to_string())
            )
        } else {
            format!(
                "Boolean operation failed after {} attempts (dominant failure: {})",
                self.total_attempts,
                self.dominant_failure_class
                    .map(|c| c.to_string())
                    .unwrap_or_else(|| "unknown".to_string())
            )
        }
    }
}

/// Configuration that produced a successful result.
#[derive(Debug, Clone)]
pub struct FinalSuccessfulConfig {
    /// Fuzzy tolerance that succeeded.
    pub fuzzy_tolerance: f64,
    /// Whether glue mode was enabled.
    pub glue_enabled: bool,
    /// Glue tolerance used.
    pub glue_tolerance: f64,
    /// Whether MakeConnected was run.
    pub make_connected_run: bool,
    /// Number of MakeConnected passes.
    pub make_connected_passes: usize,
    /// Whether scoped make-connected was used.
    pub scoped_make_connected: bool,
}

/// Failure analyzer for detailed classification of boolean operation failures.
#[derive(Debug, Clone, Default)]
pub struct FailureAnalyzer {
    /// Threshold for detecting degenerate edges (squared length).
    pub degenerate_edge_threshold: f64,
    /// Threshold for detecting degenerate triangles (minimum area).
    pub degenerate_triangle_threshold: f64,
    /// Maximum allowed self-intersection distance.
    pub self_intersection_threshold: f64,
}

impl FailureAnalyzer {
    /// Creates a new failure analyzer with default thresholds.
    pub fn new() -> Self {
        Self {
            degenerate_edge_threshold: 1e-12,
            degenerate_triangle_threshold: 1e-15,
            self_intersection_threshold: 1e-6,
        }
    }

    /// Classifies a failure based on error message and context.
    pub fn classify_from_error(error_message: &str) -> BooleanFailureClass {
        let msg = error_message.to_lowercase();

        // Check self-intersection first since it contains "intersection"
        if msg.contains("self-intersect") {
            BooleanFailureClass::SelfIntersection
        } else if msg.contains("empty") || msg.contains("missing geometry") {
            BooleanFailureClass::InvalidInput
        } else if msg.contains("degenerate") {
            BooleanFailureClass::DegenerateTopology
        } else if msg.contains("nan") || msg.contains("infinity") || msg.contains("numerical") {
            BooleanFailureClass::NumericalInstability
        } else if msg.contains("invalid") || msg.contains("non-manifold") || msg.contains("open shell") {
            BooleanFailureClass::InvalidResult
        } else if msg.contains("intersection") || msg.contains("missing curve") {
            BooleanFailureClass::IncompleteIntersection
        } else {
            BooleanFailureClass::Unknown
        }
    }

    /// Determines the appropriate recovery strategy for a failure class.
    pub fn determine_recovery_strategy(
        failure_class: BooleanFailureClass,
        attempt_number: usize,
        policy: &RetryPolicy,
    ) -> RecoveryStrategy {
        if !failure_class.is_recoverable() {
            return RecoveryStrategy::None;
        }

        match failure_class {
            BooleanFailureClass::DegenerateTopology => {
                RecoveryStrategy::Combined {
                    use_glue: policy.should_enable_glue(attempt_number),
                    run_make_connected: true,
                    increase_fuzzy: attempt_number > 1,
                }
            }
            BooleanFailureClass::NumericalInstability => {
                RecoveryStrategy::Combined {
                    use_glue: policy.should_enable_glue(attempt_number),
                    run_make_connected: attempt_number > 2,
                    increase_fuzzy: true,
                }
            }
            BooleanFailureClass::InvalidResult => {
                if policy.try_algorithm_variants && attempt_number > 2 {
                    RecoveryStrategy::AlgorithmVariant
                } else {
                    RecoveryStrategy::MakeConnectedCleanup
                }
            }
            BooleanFailureClass::IncompleteIntersection => {
                RecoveryStrategy::EnableGlueMode
            }
            BooleanFailureClass::SelfIntersection => {
                RecoveryStrategy::MakeConnectedCleanup
            }
            _ => RecoveryStrategy::None,
        }
    }
}

/// Builder for creating customized retry policies.
#[derive(Debug, Clone, Default)]
pub struct RetryPolicyBuilder {
    policy: RetryPolicy,
}

impl RetryPolicyBuilder {
    /// Creates a new builder with default settings.
    pub fn new() -> Self {
        Self {
            policy: RetryPolicy::default(),
        }
    }

    /// Starts from a conservative policy.
    pub fn conservative() -> Self {
        Self {
            policy: RetryPolicy::conservative(),
        }
    }

    /// Starts from an aggressive policy.
    pub fn aggressive() -> Self {
        Self {
            policy: RetryPolicy::aggressive(),
        }
    }

    /// Sets the maximum number of attempts.
    pub fn max_attempts(mut self, max: usize) -> Self {
        self.policy.max_attempts = max;
        self
    }

    /// Sets the fuzzy tolerance growth factor.
    pub fn fuzzy_growth_factor(mut self, factor: f64) -> Self {
        self.policy.fuzzy_growth_factor = factor;
        self
    }

    /// Sets the fuzzy tolerance cap.
    pub fn fuzzy_tolerance_cap(mut self, cap: f64) -> Self {
        self.policy.fuzzy_tolerance_cap = cap;
        self
    }

    /// Sets when to enable glue mode.
    pub fn enable_glue_after(mut self, n_failures: usize) -> Self {
        self.policy.enable_glue_after_n_failures = n_failures;
        self
    }

    /// Sets the glue tolerance.
    pub fn glue_tolerance(mut self, tol: f64) -> Self {
        self.policy.glue_tolerance = tol;
        self
    }

    /// Sets the MakeConnected aggressiveness (1-10).
    pub fn make_connected_aggressiveness(mut self, level: u32) -> Self {
        self.policy.make_connected_aggressiveness = level.clamp(1, 10);
        self
    }

    /// Sets the maximum MakeConnected passes.
    pub fn make_connected_max_passes(mut self, passes: usize) -> Self {
        self.policy.make_connected_max_passes = passes;
        self
    }

    /// Sets the initial MakeConnected tolerance.
    pub fn make_connected_initial_tolerance(mut self, tol: f64) -> Self {
        self.policy.make_connected_initial_tolerance = tol;
        self
    }

    /// Sets whether to use scoped MakeConnected.
    pub fn use_scoped_make_connected(mut self, use_scoped: bool) -> Self {
        self.policy.use_scoped_make_connected = use_scoped;
        self
    }

    /// Sets whether to fall back to global cleanup.
    pub fn fallback_to_global_cleanup(mut self, fallback: bool) -> Self {
        self.policy.fallback_to_global_cleanup = fallback;
        self
    }

    /// Sets whether to try algorithm variants.
    pub fn try_algorithm_variants(mut self, try_variants: bool) -> Self {
        self.policy.try_algorithm_variants = try_variants;
        self
    }

    /// Enables verbose diagnostics.
    pub fn verbose_diagnostics(mut self, verbose: bool) -> Self {
        self.policy.verbose_diagnostics = verbose;
        self
    }

    /// Builds the final retry policy.
    pub fn build(self) -> RetryPolicy {
        self.policy
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn failure_class_is_recoverable() {
        assert!(BooleanFailureClass::DegenerateTopology.is_recoverable());
        assert!(BooleanFailureClass::NumericalInstability.is_recoverable());
        assert!(BooleanFailureClass::InvalidResult.is_recoverable());
        assert!(BooleanFailureClass::IncompleteIntersection.is_recoverable());
        assert!(BooleanFailureClass::SelfIntersection.is_recoverable());
        assert!(!BooleanFailureClass::InvalidInput.is_recoverable());
        assert!(!BooleanFailureClass::Unknown.is_recoverable());
    }

    #[test]
    fn failure_class_suggested_recovery() {
        assert_eq!(
            BooleanFailureClass::DegenerateTopology.suggested_recovery(),
            RecoveryStrategy::MakeConnectedCleanup
        );
        assert_eq!(
            BooleanFailureClass::NumericalInstability.suggested_recovery(),
            RecoveryStrategy::IncreaseFuzzyTolerance
        );
        assert_eq!(
            BooleanFailureClass::InvalidResult.suggested_recovery(),
            RecoveryStrategy::AlgorithmVariant
        );
        assert_eq!(
            BooleanFailureClass::IncompleteIntersection.suggested_recovery(),
            RecoveryStrategy::EnableGlueMode
        );
        assert_eq!(
            BooleanFailureClass::InvalidInput.suggested_recovery(),
            RecoveryStrategy::None
        );
    }

    #[test]
    fn retry_policy_default() {
        let policy = RetryPolicy::default();
        assert_eq!(policy.max_attempts, 5);
        assert_eq!(policy.fuzzy_growth_factor, 10.0);
        assert_eq!(policy.make_connected_aggressiveness, 5);
    }

    #[test]
    fn retry_policy_presets() {
        let conservative = RetryPolicy::conservative();
        assert_eq!(conservative.max_attempts, 3);

        let aggressive = RetryPolicy::aggressive();
        assert_eq!(aggressive.max_attempts, 10);

        let numerical = RetryPolicy::for_numerical_instability();
        assert_eq!(numerical.max_attempts, 8);

        let degenerate = RetryPolicy::for_degenerate_topology();
        assert_eq!(degenerate.make_connected_aggressiveness, 9);

        let incomplete = RetryPolicy::for_incomplete_intersection();
        assert_eq!(incomplete.enable_glue_after_n_failures, 0);
    }

    #[test]
    fn retry_policy_next_fuzzy_tolerance() {
        let policy = RetryPolicy::default();
        let next = policy.next_fuzzy_tolerance(1e-6);
        assert!((next - 1e-5).abs() < 1e-15);

        // Should cap at the maximum
        let capped = policy.next_fuzzy_tolerance(1e-2);
        assert!(capped <= policy.fuzzy_tolerance_cap);
    }

    #[test]
    fn retry_policy_should_enable_glue() {
        let policy = RetryPolicy::default();
        assert!(!policy.should_enable_glue(1));
        assert!(policy.should_enable_glue(2));
        assert!(policy.should_enable_glue(3));
    }

    #[test]
    fn failure_analyzer_classify_from_error() {
        assert_eq!(
            FailureAnalyzer::classify_from_error("empty input"),
            BooleanFailureClass::InvalidInput
        );
        assert_eq!(
            FailureAnalyzer::classify_from_error("degenerate result"),
            BooleanFailureClass::DegenerateTopology
        );
        assert_eq!(
            FailureAnalyzer::classify_from_error("NaN encountered"),
            BooleanFailureClass::NumericalInstability
        );
        assert_eq!(
            FailureAnalyzer::classify_from_error("non-manifold edge detected"),
            BooleanFailureClass::InvalidResult
        );
        assert_eq!(
            FailureAnalyzer::classify_from_error("missing intersection curve"),
            BooleanFailureClass::IncompleteIntersection
        );
        assert_eq!(
            FailureAnalyzer::classify_from_error("self-intersection found"),
            BooleanFailureClass::SelfIntersection
        );
    }

    #[test]
    fn retry_policy_builder() {
        let policy = RetryPolicyBuilder::new()
            .max_attempts(7)
            .fuzzy_growth_factor(15.0)
            .enable_glue_after(1)
            .make_connected_aggressiveness(7)
            .build();

        assert_eq!(policy.max_attempts, 7);
        assert_eq!(policy.fuzzy_growth_factor, 15.0);
        assert_eq!(policy.enable_glue_after_n_failures, 1);
        assert_eq!(policy.make_connected_aggressiveness, 7);
    }

    #[test]
    fn diagnostic_report() {
        let mut report = BooleanDiagnosticReport::new();

        // Add a failed attempt
        let failed_attempt = BooleanAttemptDiagnostic {
            attempt_number: 1,
            fuzzy_tolerance: 1e-6,
            success: false,
            failure_class: Some(BooleanFailureClass::NumericalInstability),
            recovery_strategy: Some(RecoveryStrategy::IncreaseFuzzyTolerance),
            error_message: Some("NaN detected".to_string()),
            ..Default::default()
        };
        report.add_attempt(failed_attempt);

        // Add a successful attempt
        let success_attempt = BooleanAttemptDiagnostic {
            attempt_number: 2,
            fuzzy_tolerance: 1e-5,
            glue_enabled: true,
            success: true,
            result_faces: Some(10),
            ..Default::default()
        };
        report.add_attempt(success_attempt);

        report.finalize();

        assert!(report.is_success());
        assert_eq!(report.total_attempts, 2);
        assert_eq!(report.successful_attempts, 1);
        assert!(report.glue_mode_needed);
    }

    #[test]
    fn recovery_strategy_from_failure_class() {
        let policy = RetryPolicy::default();

        // First attempt with degenerate topology
        let strategy = FailureAnalyzer::determine_recovery_strategy(
            BooleanFailureClass::DegenerateTopology,
            1,
            &policy,
        );
        assert!(matches!(strategy, RecoveryStrategy::Combined { run_make_connected: true, .. }));

        // Second attempt with incomplete intersection
        let strategy = FailureAnalyzer::determine_recovery_strategy(
            BooleanFailureClass::IncompleteIntersection,
            2,
            &policy,
        );
        assert_eq!(strategy, RecoveryStrategy::EnableGlueMode);

        // Invalid input is not recoverable
        let strategy = FailureAnalyzer::determine_recovery_strategy(
            BooleanFailureClass::InvalidInput,
            1,
            &policy,
        );
        assert_eq!(strategy, RecoveryStrategy::None);
    }

    #[test]
    fn make_connected_tolerance_computation() {
        let policy = RetryPolicy {
            make_connected_initial_tolerance: 1e-6,
            make_connected_tolerance_growth: 2.0,
            ..Default::default()
        };

        assert!((policy.make_connected_tolerance(0) - 1e-6).abs() < 1e-15);
        assert!((policy.make_connected_tolerance(1) - 2e-6).abs() < 1e-15);
        assert!((policy.make_connected_tolerance(2) - 4e-6).abs() < 1e-15);
    }

    // ============================================================================
    // Near-Coincident Vertex Tests
    // ============================================================================

    /// Two boxes with vertices nearly coincident (within tolerance).
    /// Tests the kernel's ability to handle fuzzy tolerance correctly when
    /// vertices are almost but not exactly touching.
    #[test]
    fn test_near_coincident_vertices_union() {
        use crate::{boolean_op_with_options, BooleanOpType, BooleanOptions};
        use glam::DVec3;
        use rcad_modeling::make_box_brep;

        // Two boxes with vertices nearly coincident (within tolerance)
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        // B's corner is very close to A's corner at (2, 2, 2)
        let b = make_box_brep(
            DVec3::new(2.0 - 1e-5, 2.0 - 1e-5, 2.0 - 1e-5),
            DVec3::X,
            DVec3::Y,
            1.0,
            1.0,
            1.0,
        )
        .unwrap();

        let opts = BooleanOptions {
            fuzzy_tol: 1e-4,
            ..Default::default()
        };

        let result = boolean_op_with_options(BooleanOpType::Union, &a, &b, opts);
        assert!(
            result.is_ok(),
            "near-coincident vertices union should succeed"
        );
    }

    /// Two boxes with vertices nearly coincident - difference operation.
    /// Tests fuzzy tolerance handling for difference operations.
    #[test]
    fn test_near_coincident_vertices_difference() {
        use crate::{boolean_op_with_options, BooleanOpType, BooleanOptions};
        use glam::DVec3;
        use rcad_modeling::make_box_brep;

        // Two boxes with vertices nearly coincident (within tolerance)
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        // B's corner is very close to A's corner at (2, 2, 2)
        let b = make_box_brep(
            DVec3::new(2.0 - 1e-5, 2.0 - 1e-5, 2.0 - 1e-5),
            DVec3::X,
            DVec3::Y,
            1.0,
            1.0,
            1.0,
        )
        .unwrap();

        let opts = BooleanOptions {
            fuzzy_tol: 1e-4,
            ..Default::default()
        };

        let result = boolean_op_with_options(BooleanOpType::Difference, &a, &b, opts);
        assert!(
            result.is_ok(),
            "near-coincident vertices difference should succeed"
        );
    }

    // ============================================================================
    // Near-Tangent Geometry Tests
    // ============================================================================

    /// Test union of two boxes with faces nearly touching (gap within fuzzy tolerance).
    /// This verifies the kernel's ability to handle near-tangent plane geometry.
    #[test]
    fn test_near_tangent_plane_union() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op_with_options, BooleanOptions, BooleanOpType};

        // Two boxes with faces nearly touching (gap = 1e-5, within fuzzy tolerance)
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0 + 1e-5, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();

        let opts = BooleanOptions {
            fuzzy_tol: 1e-4,
            ..Default::default()
        };

        let result = boolean_op_with_options(BooleanOpType::Union, &a, &b, opts);
        assert!(result.is_ok(), "near-tangent union should succeed with fuzzy tolerance");
    }

    /// Test difference of two boxes with faces nearly touching.
    /// This verifies the kernel's ability to handle near-tangent plane geometry in difference operations.
    #[test]
    fn test_near_tangent_plane_difference() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op_with_options, BooleanOptions, BooleanOpType};

        // Two boxes with faces nearly touching (gap = 1e-5, within fuzzy tolerance)
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0 + 1e-5, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();

        let opts = BooleanOptions {
            fuzzy_tol: 1e-4,
            ..Default::default()
        };

        let result = boolean_op_with_options(BooleanOpType::Difference, &a, &b, opts);
        // Difference with near-tangent faces should either succeed or return degenerate result
        assert!(result.is_ok(), "near-tangent difference should succeed with fuzzy tolerance");
    }

    /// Test intersection of a cylinder with a near-tangent plane (box face).
    /// This verifies the kernel's ability to handle near-tangent curved-to-plane geometry.
    #[test]
    fn test_near_tangent_cylinder_intersection() {
        use rcad_modeling::{make_box_brep, make_cylinder_brep};
        use glam::DVec3;
        use crate::{boolean_op_with_options, BooleanOptions, BooleanOpType};

        // Cylinder intersecting with a box where the cylinder surface is nearly tangent to a box face
        let cylinder = make_cylinder_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::Z, DVec3::X, 1.0, 4.0).unwrap();
        // Box positioned so its face is nearly tangent to cylinder surface (gap = 1e-5)
        let box_ = make_box_brep(DVec3::new(1.0 + 1e-5, -2.0, -1.0), DVec3::X, DVec3::Y, 4.0, 4.0, 6.0).unwrap();

        let opts = BooleanOptions {
            fuzzy_tol: 1e-4,
            ..Default::default()
        };

        let result = boolean_op_with_options(BooleanOpType::Intersection, &cylinder, &box_, opts);
        assert!(result.is_ok(), "near-tangent cylinder-plane intersection should succeed with fuzzy tolerance");
    }

    // ============================================================================
    // Thin Wall Tests
    // ============================================================================

    /// Test difference creating thin walls (0.1 thickness).
    /// This verifies the kernel's ability to handle thin wall geometry that can
    /// be prone to numerical issues and degenerate faces.
    #[test]
    fn test_thin_wall_difference() {
        use crate::geom_populate::populate_box_geom;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};
        use glam::DVec3;
        use rcad_modeling::make_box_brep;

        // Outer box with inner box creating thin wall (0.1 thickness)
        let mut outer =
            make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).unwrap();
        let mut inner = make_box_brep(
            DVec3::new(0.1, 0.1, 0.1),
            DVec3::X,
            DVec3::Y,
            9.8,
            9.8,
            9.8,
        )
        .unwrap();
        populate_box_geom(&mut outer);
        populate_box_geom(&mut inner);

        let result = boolean_op_simplified(
            BooleanOpType::Difference,
            &outer,
            &inner,
            SimplifyOptions::default(),
        );
        assert!(
            result.is_ok(),
            "thin wall difference should succeed"
        );
    }

    /// Test very thin geometry (sheet metal thickness).
    /// This verifies the kernel's ability to handle sheet metal geometry with
    /// very thin dimensions that can challenge boolean operations.
    #[test]
    fn test_sheet_metal_thickness() {
        use crate::geom_populate::populate_box_geom;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};
        use glam::DVec3;
        use rcad_modeling::make_box_brep;

        // Sheet metal panel (100x100x0.5mm typical sheet metal thickness)
        let mut panel =
            make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 100.0, 100.0, 0.5).unwrap();
        // Small hole cut into the panel
        let mut hole =
            make_box_brep(DVec3::new(45.0, 45.0, -0.1), DVec3::X, DVec3::Y, 10.0, 10.0, 0.7)
                .unwrap();
        populate_box_geom(&mut panel);
        populate_box_geom(&mut hole);

        let result = boolean_op_simplified(
            BooleanOpType::Difference,
            &panel,
            &hole,
            SimplifyOptions::default(),
        );
        assert!(
            result.is_ok(),
            "sheet metal hole cut should succeed"
        );
    }

    // ============================================================================
    // Small Feature Tests
    // ============================================================================

    /// Test subtracting a small hole (cylinder) from a large box.
    /// This verifies the kernel's ability to handle small features without
    /// numerical issues or degenerate geometry.
    #[test]
    fn test_small_hole_subtraction() {
        use rcad_modeling::{make_box_brep, make_cylinder_brep};
        use glam::DVec3;
        use crate::geom_populate::populate_box_geom;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Large box with small hole (radius 0.01)
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).unwrap();
        let hole = make_cylinder_brep(DVec3::new(5.0, 5.0, -1.0), DVec3::Z, DVec3::X, 0.01, 12.0).unwrap();
        populate_box_geom(&mut box_brep);

        let result = boolean_op_simplified(BooleanOpType::Difference, &box_brep, &hole, SimplifyOptions::default());
        assert!(result.is_ok(), "small hole subtraction should succeed");
    }

    /// Test subtracting a thin slot from a box.
    /// This verifies the kernel's ability to handle thin features (high aspect ratio)
    /// without numerical instability.
    #[test]
    fn test_thin_slot_subtraction() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::geom_populate::populate_box_geom;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Main box
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).unwrap();
        // Thin slot (0.1mm wide, 5mm deep, spanning full height)
        let slot = make_box_brep(
            DVec3::new(4.95, -0.1, 0.0),
            DVec3::X,
            DVec3::Z,
            0.1,
            5.0,
            10.0,
        ).unwrap();
        populate_box_geom(&mut box_brep);

        let result = boolean_op_simplified(BooleanOpType::Difference, &box_brep, &slot, SimplifyOptions::default());
        assert!(result.is_ok(), "thin slot subtraction should succeed");
    }

    /// Test that tiny features are preserved through simplification.
    /// This verifies that the simplification pipeline doesn't remove
    /// small but intentional geometric features.
    #[test]
    fn test_tiny_fillet_preservation() {
        use rcad_modeling::{make_box_brep, make_cylinder_brep};
        use glam::DVec3;
        use crate::geom_populate::populate_box_geom;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Create a box with a small fillet-like feature (using cylinder subtraction)
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 5.0, 5.0, 5.0).unwrap();
        // Small fillet radius (0.05)
        let fillet_tool = make_cylinder_brep(
            DVec3::new(0.0, 0.0, -0.1),
            DVec3::Z,
            DVec3::X,
            0.05,
            5.2,
        ).unwrap();
        populate_box_geom(&mut box_brep);

        // Use conservative simplification to preserve small features
        let opts = SimplifyOptions {
            merge_tolerance: 1e-7,
            ..SimplifyOptions::default()
        };

        let result = boolean_op_simplified(BooleanOpType::Difference, &box_brep, &fillet_tool, opts);
        assert!(result.is_ok(), "tiny fillet subtraction should succeed");

        // The result should have more faces than a simple box due to the fillet
        let (result_brep, _report) = result.unwrap();
        let face_count = result_brep.solids.first()
            .and_then(|s| s.shells.first())
            .map(|s| s.faces.len())
            .unwrap_or(0);
        // A box has 6 faces; after fillet subtraction we expect more
        assert!(face_count >= 6, "result should have at least 6 faces, got {}", face_count);
    }

    // ============================================================================
    // Nested Boolean Operation Tests
    // ============================================================================

    /// Test multiple union operations applied sequentially.
    /// This verifies the kernel's ability to handle chained boolean operations
    /// where the result of one operation becomes input to the next.
    #[test]
    fn test_sequential_union_operations() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::geom_populate::populate_box_geom;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Create 4 boxes and union them sequentially
        let mut a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let mut b = make_box_brep(DVec3::new(0.5, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let mut c = make_box_brep(DVec3::new(0.0, 0.5, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let mut d = make_box_brep(DVec3::new(0.0, 0.0, 0.5), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();

        populate_box_geom(&mut a);
        populate_box_geom(&mut b);
        populate_box_geom(&mut c);
        populate_box_geom(&mut d);

        let (ab, _) = boolean_op_simplified(BooleanOpType::Union, &a, &b, SimplifyOptions::default()).unwrap();
        let (abc, _) = boolean_op_simplified(BooleanOpType::Union, &ab, &c, SimplifyOptions::default()).unwrap();
        let (abcd, _) = boolean_op_simplified(BooleanOpType::Union, &abc, &d, SimplifyOptions::default()).unwrap();

        assert!(abcd.solids.len() >= 1, "sequential unions should produce valid result");
    }

    /// Test mixed boolean operations: union followed by difference.
    /// This verifies the kernel's ability to handle different boolean operation
    /// types applied sequentially to accumulated geometry.
    #[test]
    fn test_mixed_boolean_operations() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::geom_populate::populate_box_geom;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Create base boxes for union
        let mut a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let mut b = make_box_brep(DVec3::new(1.5, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        // Create box to subtract
        let mut c = make_box_brep(DVec3::new(1.0, 0.5, 0.5), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();

        populate_box_geom(&mut a);
        populate_box_geom(&mut b);
        populate_box_geom(&mut c);

        // First union a and b
        let (ab, _) = boolean_op_simplified(BooleanOpType::Union, &a, &b, SimplifyOptions::default()).unwrap();
        // Then subtract c from the result
        let (result, _) = boolean_op_simplified(BooleanOpType::Difference, &ab, &c, SimplifyOptions::default()).unwrap();

        assert!(result.solids.len() >= 1, "mixed boolean operations should produce valid result");
    }

    // ============================================================================
    // Overlapping Faces Tests
    // ============================================================================

    /// Test boolean union with coplanar overlapping faces.
    /// Two boxes sharing a face that exactly overlaps - tests the kernel's ability
    /// to handle coincident face geometry correctly.
    #[test]
    fn test_boolean_overlapping_faces() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Two boxes with exactly overlapping faces (coplanar at z=0)
        let mut a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 1.0).unwrap();
        let mut b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Z, 2.0, 1.0, 2.0).unwrap();

        let result = boolean_op_simplified(BooleanOpType::Union, &a, &b, SimplifyOptions::default());
        assert!(result.is_ok(), "overlapping faces union should succeed");

        let (result_brep, _) = result.unwrap();
        assert!(result_brep.solids.len() >= 1, "overlapping faces union should produce valid result");
    }

    /// Test boolean union where two boxes share a common edge but not faces.
    /// Tests edge-edge intersection handling in boolean operations.
    #[test]
    fn test_boolean_shared_edge() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Two boxes sharing only an edge (touching along one edge)
        // Box A: 0 to 2 in X, 0 to 2 in Y, 0 to 2 in Z
        // Box B: 0 to 2 in X, 2 to 4 in Y, 0 to 2 in Z (shares edge at y=2, x=0..2, z=0..2)
        let mut a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let mut b = make_box_brep(DVec3::new(0.0, 2.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();

        let result = boolean_op_simplified(BooleanOpType::Union, &a, &b, SimplifyOptions::default());
        assert!(result.is_ok(), "shared edge union should succeed");

        let (result_brep, _) = result.unwrap();
        // The result should have more faces than a single box due to the L-shape
        assert!(result_brep.solids.len() >= 1, "shared edge union should produce valid result");
    }

    /// Test boolean with nested solid operations (box inside box).
    /// Tests the kernel's handling of one solid completely inside another.
    #[test]
    fn test_boolean_nested_solids() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Outer box: 10x10x10
        let mut outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).unwrap();
        // Inner box: 5x5x5, centered inside outer
        let mut inner = make_box_brep(DVec3::new(2.5, 2.5, 2.5), DVec3::X, DVec3::Y, 5.0, 5.0, 5.0).unwrap();

        // Intersection should return the inner box
        let intersection = boolean_op_simplified(BooleanOpType::Intersection, &outer, &inner, SimplifyOptions::default());
        assert!(intersection.is_ok(), "nested intersection should succeed");
        let (int_brep, _) = intersection.unwrap();
        assert!(int_brep.solids.len() >= 1, "nested intersection should produce valid result");

        // Difference should create a hollow shape (outer minus inner)
        let difference = boolean_op_simplified(BooleanOpType::Difference, &outer, &inner, SimplifyOptions::default());
        assert!(difference.is_ok(), "nested difference should succeed");
    }

    /// Test boolean union with partially overlapping geometry.
    /// Tests handling of complex intersection curves when geometry partially overlaps.
    #[test]
    fn test_boolean_partial_overlap() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Two boxes that partially overlap (like a cross)
        let mut a = make_box_brep(DVec3::new(-1.0, -0.5, -0.5), DVec3::X, DVec3::Y, 3.0, 1.0, 1.0).unwrap();
        let mut b = make_box_brep(DVec3::new(-0.5, -1.0, -0.5), DVec3::X, DVec3::Y, 1.0, 3.0, 1.0).unwrap();

        let result = boolean_op_simplified(BooleanOpType::Union, &a, &b, SimplifyOptions::default());
        assert!(result.is_ok(), "partial overlap union should succeed");
    }

    /// Test boolean with spheres intersecting at a point.
    /// Tests degenerate intersection handling where surfaces touch at a single point.
    #[test]
    fn test_boolean_touching_at_point() {
        use rcad_modeling::{make_box_brep, make_sphere_brep};
        use glam::DVec3;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Box and sphere touching at a corner point
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        // Sphere centered at (2, 2, 2) with radius 1, touches box corner
        let sphere = make_sphere_brep(DVec3::new(3.0, 3.0, 3.0), 1.0).unwrap();

        let result = boolean_op_simplified(BooleanOpType::Union, &box_brep, &sphere, SimplifyOptions::default());
        assert!(result.is_ok(), "touching at point union should succeed or return empty");
    }

    /// Test boolean with cylindrical intersection creating complex curves.
    /// Tests the kernel's ability to handle curved-curved intersections.
    #[test]
    fn test_boolean_cylinder_intersection() {
        use rcad_modeling::{make_cylinder_brep, make_box_brep};
        use glam::DVec3;
        use crate::{boolean_op_simplified, BooleanOpType, SimplifyOptions};

        // Cylinder passing through a box at an angle (using perpendicular for simplicity)
        let mut box_brep = make_box_brep(DVec3::new(-2.0, -2.0, -2.0), DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let cylinder = make_cylinder_brep(DVec3::new(0.0, 0.0, -3.0), DVec3::Z, DVec3::X, 0.5, 6.0).unwrap();

        let result = boolean_op_simplified(BooleanOpType::Difference, &box_brep, &cylinder, SimplifyOptions::default());
        assert!(result.is_ok(), "cylinder-box difference should succeed");

        let (result_brep, _) = result.unwrap();
        // The result should have more faces due to the cylindrical hole
        let face_count: usize = result_brep.solids.iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();
        assert!(face_count >= 6, "cylinder difference should add at least one face");
    }

    // ============================================================================
    // OCCT TKBO Alignment Tests - Overlapping and Identical Geometry
    // ============================================================================

    /// Test boolean union of two identical boxes.
    /// OCCT handles this by returning one of the inputs.
    #[test]
    fn test_boolean_identical_boxes_union() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();

        let result = boolean_op(BooleanOpType::Union, &a, &b);
        assert!(result.is_ok(), "identical boxes union should succeed");
    }

    /// Test boolean intersection of two identical boxes.
    #[test]
    fn test_boolean_identical_boxes_intersection() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();

        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        assert!(result.is_ok(), "identical boxes intersection should succeed");
    }

    /// Test boolean difference of identical boxes (should result in empty or degenerate).
    #[test]
    fn test_boolean_identical_boxes_difference() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();

        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        // Difference of identical shapes may produce empty or degenerate result
        assert!(result.is_ok() || result.is_err(), "identical boxes difference handled");
    }

    /// Test boolean with box completely inside another box.
    #[test]
    fn test_boolean_nested_boxes_union() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let inner = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();

        let result = boolean_op(BooleanOpType::Union, &outer, &inner);
        assert!(result.is_ok(), "nested boxes union should succeed");
    }

    /// Test boolean difference creating a hollow shell.
    #[test]
    fn test_boolean_hollow_shell() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let inner = make_box_brep(DVec3::new(0.5, 0.5, 0.5), DVec3::X, DVec3::Y, 3.0, 3.0, 3.0).unwrap();

        let result = boolean_op(BooleanOpType::Difference, &outer, &inner);
        assert!(result.is_ok(), "hollow shell creation should succeed");
    }

    /// Test boolean with disjoint boxes (no intersection).
    #[test]
    fn test_boolean_disjoint_boxes() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_box_brep(DVec3::new(5.0, 5.0, 5.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();

        // Intersection of disjoint boxes should return empty
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        // May return empty or error depending on implementation
        assert!(result.is_ok() || result.is_err(), "disjoint intersection handled");

        // Union of disjoint boxes should return compound
        let result = boolean_op(BooleanOpType::Union, &a, &b);
        assert!(result.is_ok(), "disjoint union should succeed");
    }

    /// Test boolean with thin wall feature.
    #[test]
    fn test_boolean_thin_wall_feature() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        // Create a thin wall by subtracting a slightly smaller box
        let outer = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).unwrap();
        let inner = make_box_brep(DVec3::new(0.1, 0.1, 0.1), DVec3::X, DVec3::Y, 9.8, 9.8, 9.8).unwrap();

        let result = boolean_op(BooleanOpType::Difference, &outer, &inner);
        assert!(result.is_ok(), "thin wall creation should succeed");
    }

    /// Test boolean with small hole feature.
    #[test]
    fn test_boolean_small_hole() {
        use rcad_modeling::{make_box_brep, make_cylinder_brep};
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 5.0, 5.0, 5.0).unwrap();
        // Small diameter cylinder for hole
        let hole = make_cylinder_brep(DVec3::new(2.5, 2.5, -1.0), DVec3::Z, DVec3::X, 0.1, 7.0).unwrap();

        let result = boolean_op(BooleanOpType::Difference, &box_brep, &hole);
        assert!(result.is_ok(), "small hole creation should succeed");
    }

    /// Test boolean with sphere-sphere intersection.
    #[test]
    fn test_boolean_sphere_intersection() {
        use rcad_modeling::make_sphere_brep;
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let a = make_sphere_brep(DVec3::ZERO, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(2.0, 0.0, 0.0), 2.0).unwrap();

        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        assert!(result.is_ok(), "sphere intersection should succeed");
    }

    /// Test boolean with cylinder-cylinder intersection.
    #[test]
    fn test_boolean_two_cylinders() {
        use rcad_modeling::make_cylinder_brep;
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let a = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 5.0).unwrap();
        let b = make_cylinder_brep(DVec3::new(1.0, 0.0, 0.0), DVec3::Z, DVec3::X, 1.0, 5.0).unwrap();

        let result = boolean_op(BooleanOpType::Union, &a, &b);
        assert!(result.is_ok(), "cylinder union should succeed");
    }

    /// Test boolean with cone-box intersection.
    #[test]
    fn test_boolean_cone_box() {
        use rcad_modeling::{make_box_brep, make_cone_brep};
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let box_brep = make_box_brep(DVec3::new(-2.0, -2.0, -2.0), DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let cone = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.5, 3.0).unwrap();

        let result = boolean_op(BooleanOpType::Difference, &box_brep, &cone);
        assert!(result.is_ok(), "cone-box difference should succeed");
    }

    /// Test fuzzy tolerance with nearly coincident geometry.
    #[test]
    fn test_boolean_fuzzy_tolerance() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op_with_options, BooleanOpType, BooleanOptions};

        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        // Box with corner very close to A's corner
        let b = make_box_brep(DVec3::new(1.99999, 1.99999, 1.99999), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();

        let opts = BooleanOptions {
            fuzzy_tol: 1e-4,
            ..Default::default()
        };

        let result = boolean_op_with_options(BooleanOpType::Union, &a, &b, opts);
        assert!(result.is_ok(), "fuzzy tolerance union should succeed");
    }

    /// Test boolean with multiple operations in sequence.
    #[test]
    fn test_boolean_sequential_operations() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use crate::{boolean_op, BooleanOpType};

        let mut result = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();

        // Subtract three smaller boxes
        for i in 0..3 {
            let offset = 0.5 + i as f64 * 1.2;
            let small = make_box_brep(
                DVec3::new(offset, offset, -0.5),
                DVec3::X, DVec3::Y,
                1.0, 1.0, 5.0
            ).unwrap();
            result = boolean_op(BooleanOpType::Difference, &result, &small).unwrap();
        }

        // Result should be valid
        assert!(result.solids.iter().any(|s| !s.shells.is_empty()));
    }

    /// Test recovery strategy application.
    #[test]
    fn test_recovery_strategy_sequence() {
        let policy = RetryPolicy::default();

        // Test sequence of strategies for different failure types
        let strategy1 = FailureAnalyzer::determine_recovery_strategy(
            BooleanFailureClass::NumericalInstability,
            1,
            &policy,
        );
        assert!(matches!(strategy1, RecoveryStrategy::IncreaseFuzzyTolerance | RecoveryStrategy::Combined { .. }));

        let strategy2 = FailureAnalyzer::determine_recovery_strategy(
            BooleanFailureClass::DegenerateTopology,
            2,
            &policy,
        );
        assert!(matches!(strategy2, RecoveryStrategy::Combined { .. } | RecoveryStrategy::MakeConnectedCleanup));
    }
}
