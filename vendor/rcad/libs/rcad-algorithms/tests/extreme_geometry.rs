/// Tests for extreme geometry detection and handling in boolean operations.
use glam::DVec3;
use rcad_algorithms::{
    AspectRatioAdaptiveTolerance, DegenerateGeometryHandler, DegenerateType,
    HighAspectRatioEdge, NearDegenerateGeometry, NearTangentHandler, NearTangentSeverity,
    SizeDifferenceHandler, SizeDifferenceAnalysis,
    ExtremeGeometryAnalysis, ExtremeGeometryAnalysisOptions,
    analyze_extreme_geometry, analyze_size_difference,
    detect_high_aspect_ratio_edges, detect_near_degenerate_geometry,
    detect_near_tangent_configurations,
    ExtremeGeometryRetryConfig, ExtremeGeometryRetryPolicy,
    ASPECT_RATIO_THRESHOLD, ASPECT_RATIO_VERY_HIGH, SIZE_RATIO_THRESHOLD,
};
use rcad_kernel::BRep;
use rcad_kernel::geom::{Surface3, Plane, SphericalSurface};
use rcad_algorithms::tolerance::{TOLERANCE_ABS, AdaptiveTolerance, ToleranceLevel};

// ── Near-Tangent Geometry Tests ────────────────────────────────────────────────

#[test]
fn near_tangent_severity_critical() {
    let handler = NearTangentHandler::default();
    // 0 degrees = exactly tangent
    assert_eq!(handler.classify_severity(0.0), NearTangentSeverity::Critical);
}

#[test]
fn near_tangent_severity_near_tangent() {
    let handler = NearTangentHandler::default();
    // 0.001 degrees = very close to tangent
    let angle = 0.001_f64.to_radians();
    assert_eq!(handler.classify_severity(angle), NearTangentSeverity::NearTangent);
}

#[test]
fn near_tangent_severity_marginal() {
    let handler = NearTangentHandler::default();
    // 0.05 degrees = marginal
    let angle = 0.05_f64.to_radians();
    assert_eq!(handler.classify_severity(angle), NearTangentSeverity::Marginal);
}

#[test]
fn near_tangent_severity_not_tangent() {
    let handler = NearTangentHandler::default();
    // 1 degree = clearly not tangent
    let angle = 1.0_f64.to_radians();
    assert_eq!(handler.classify_severity(angle), NearTangentSeverity::NotTangent);
}

#[test]
fn near_tangent_handler_default() {
    let handler = NearTangentHandler::default();
    assert!(handler.base_tolerance > 0.0);
    assert!(handler.angular_threshold > 0.0);
    assert!(handler.fuzzy_multiplier > 1.0);
    assert!(handler.max_fuzzy > handler.base_tolerance);
}

#[test]
fn near_tangent_handler_from_adaptive() {
    let adaptive = AdaptiveTolerance::from_scale(100.0);
    let handler = NearTangentHandler::from_adaptive(adaptive);
    assert!((handler.base_tolerance - adaptive.coincidence()).abs() < 1e-15);
}

#[test]
fn near_tangent_fuzzy_adjustment_critical() {
    let handler = NearTangentHandler::default();
    // Critical cases should get large fuzzy adjustment
    let adjustment = handler.compute_fuzzy_adjustment(0.0);
    assert!(adjustment > handler.base_tolerance * 100.0);
}

#[test]
fn near_tangent_fuzzy_adjustment_marginal() {
    let handler = NearTangentHandler::default();
    // Marginal cases should get smaller fuzzy adjustment
    let angle = 0.05_f64.to_radians();
    let adjustment = handler.compute_fuzzy_adjustment(angle);
    assert!(adjustment >= handler.base_tolerance);
    // Marginal adjustment should be less than critical (1000x)
    assert!(adjustment < handler.base_tolerance * 1000.0);
}

#[test]
fn near_tangent_adjust_tolerance_empty_configs() {
    let handler = NearTangentHandler::default();
    let base_fuzzy = TOLERANCE_ABS;
    let adjusted = handler.adjust_tolerance_for_tangency(base_fuzzy, &[]);
    assert!((adjusted - base_fuzzy).abs() < 1e-15);
}

#[test]
fn near_tangent_adjust_tolerance_with_configs() {
    let handler = NearTangentHandler::default();
    let base_fuzzy = TOLERANCE_ABS;

    // Create a mock config with high fuzzy adjustment
    let config = rcad_algorithms::NearTangentConfig {
        point: DVec3::ZERO,
        normal_a: DVec3::Z,
        normal_b: DVec3::Z,
        angle: 0.0,
        severity: NearTangentSeverity::Critical,
        suggested_fuzzy_adjustment: TOLERANCE_ABS * 100.0,
    };

    let adjusted = handler.adjust_tolerance_for_tangency(base_fuzzy, &[config]);
    assert!(adjusted >= TOLERANCE_ABS * 100.0);
}

// ── High Aspect Ratio Tests ────────────────────────────────────────────────────

#[test]
fn aspect_ratio_thresholds() {
    assert!(ASPECT_RATIO_THRESHOLD > 0.0);
    assert!(ASPECT_RATIO_VERY_HIGH > ASPECT_RATIO_THRESHOLD);
}

#[test]
fn aspect_ratio_adaptive_tolerance_default() {
    let aat = AspectRatioAdaptiveTolerance::default();
    assert!(aat.base_tolerance > 0.0);
    assert!(aat.aspect_ratio_threshold > 0.0);
    assert!(aat.max_multiplier > 1.0);
}

#[test]
fn aspect_ratio_tolerance_multiplier_normal() {
    let aat = AspectRatioAdaptiveTolerance::default();
    // Below threshold should return 1.0
    let mult = aat.compute_tolerance_multiplier(10.0);
    assert!((mult - 1.0).abs() < 1e-10);
}

#[test]
fn aspect_ratio_tolerance_multiplier_high() {
    let aat = AspectRatioAdaptiveTolerance::default();
    // Above threshold should return > 1.0
    let mult = aat.compute_tolerance_multiplier(ASPECT_RATIO_THRESHOLD * 10.0);
    assert!(mult > 1.0);
    assert!(mult < aat.max_multiplier);
}

#[test]
fn aspect_ratio_effective_tolerance() {
    let aat = AspectRatioAdaptiveTolerance::default();
    let tol = aat.effective_tolerance(ASPECT_RATIO_THRESHOLD * 100.0);
    assert!(tol > aat.base_tolerance);
}

// ── Near-Degenerate Geometry Tests ──────────────────────────────────────────────

#[test]
fn degenerate_geometry_handler_default() {
    let handler = DegenerateGeometryHandler::default();
    assert!(handler.zero_tolerance > 0.0);
    assert!(handler.collinear_tolerance > 0.0);
    assert!(handler.min_area > 0.0);
    assert!(handler.min_edge_length > 0.0);
}

#[test]
fn degenerate_geometry_handler_from_adaptive() {
    let adaptive = AdaptiveTolerance::from_scale(100.0);
    let handler = DegenerateGeometryHandler::from_adaptive(adaptive);
    assert!((handler.zero_tolerance - adaptive.coincidence()).abs() < 1e-15);
}

#[test]
fn degenerate_type_variants() {
    // Just ensure all variants exist and can be compared
    assert!(DegenerateType::NearZeroLengthEdge != DegenerateType::NearZeroCurvature);
    assert!(DegenerateType::NearZeroAreaFace != DegenerateType::NearCollinearBoundary);
}

// ── Size Difference Tests ────────────────────────────────────────────────────────

#[test]
fn size_ratio_threshold() {
    assert!(SIZE_RATIO_THRESHOLD > 0.0);
}

#[test]
fn size_difference_handler_default() {
    let handler = SizeDifferenceHandler::default();
    assert!(handler.base_tolerance > 0.0);
    assert!(handler.size_ratio_threshold > 0.0);
    assert!(handler.max_multiplier > 1.0);
}

#[test]
fn size_difference_handler_from_adaptive() {
    let adaptive = AdaptiveTolerance::from_scale(100.0);
    let handler = SizeDifferenceHandler::from_adaptive(adaptive);
    assert!((handler.base_tolerance - adaptive.coincidence()).abs() < 1e-15);
}

#[test]
fn size_difference_compute_characteristic_size_empty() {
    let handler = SizeDifferenceHandler::default();
    let brep = BRep::default();
    let size = handler.compute_characteristic_size(&brep);
    assert!((size - 1.0).abs() < 1e-10);
}

// ── Extreme Geometry Retry Policy Tests ──────────────────────────────────────────

#[test]
fn extreme_geometry_retry_policy_variants() {
    assert!(ExtremeGeometryRetryPolicy::None != ExtremeGeometryRetryPolicy::PreAnalyze);
    assert!(ExtremeGeometryRetryPolicy::AdaptiveTolerance != ExtremeGeometryRetryPolicy::GeometryAware);
}

#[test]
fn extreme_geometry_retry_config_default() {
    let config = ExtremeGeometryRetryConfig::default();
    assert_eq!(config.policy, ExtremeGeometryRetryPolicy::AdaptiveTolerance);
    assert!(config.check_near_tangent);
    assert!(config.check_aspect_ratio);
    assert!(config.check_degenerate);
    assert!(config.check_size_difference);
    assert!(config.max_fuzzy_multiplier > 1.0);
    assert!(config.extra_retry_steps > 0);
}

#[test]
fn extreme_geometry_retry_config_none() {
    let config = ExtremeGeometryRetryConfig::none();
    assert_eq!(config.policy, ExtremeGeometryRetryPolicy::None);
    assert!(!config.check_near_tangent);
    assert!(!config.check_aspect_ratio);
    assert!(!config.check_degenerate);
    assert!(!config.check_size_difference);
}

#[test]
fn extreme_geometry_retry_config_geometry_aware() {
    let config = ExtremeGeometryRetryConfig::geometry_aware();
    assert_eq!(config.policy, ExtremeGeometryRetryPolicy::GeometryAware);
}

#[test]
fn extreme_geometry_build_retry_ladder_none() {
    let config = ExtremeGeometryRetryConfig::none();
    let base_ladder = vec![TOLERANCE_ABS * 10.0, TOLERANCE_ABS * 100.0];
    let analysis = ExtremeGeometryAnalysis::default();

    let ladder = config.build_retry_ladder(&base_ladder, &analysis);
    assert_eq!(ladder.len(), base_ladder.len());
}

#[test]
fn extreme_geometry_build_retry_ladder_with_near_tangent() {
    let config = ExtremeGeometryRetryConfig::default();
    let base_ladder = vec![TOLERANCE_ABS * 10.0];

    let mut analysis = ExtremeGeometryAnalysis::default();
    analysis.near_tangent_configs.push(rcad_algorithms::NearTangentConfig {
        point: DVec3::ZERO,
        normal_a: DVec3::Z,
        normal_b: DVec3::Z,
        angle: 0.0,
        severity: NearTangentSeverity::Critical,
        suggested_fuzzy_adjustment: TOLERANCE_ABS * 500.0,
    });

    let ladder = config.build_retry_ladder(&base_ladder, &analysis);
    // Should include the suggested fuzzy adjustment
    let expected = TOLERANCE_ABS * 500.0;
    assert!(ladder.len() > base_ladder.len() || ladder.iter().any(|&v| (v - expected).abs() < TOLERANCE_ABS));
}

#[test]
fn extreme_geometry_build_retry_ladder_with_size_difference() {
    let config = ExtremeGeometryRetryConfig::default();
    let base_ladder = vec![TOLERANCE_ABS * 10.0];

    let mut analysis = ExtremeGeometryAnalysis::default();
    analysis.size_difference = Some(SizeDifferenceAnalysis {
        size_a: 1000.0,
        size_b: 1.0,
        size_ratio: 1000.0,
        is_extreme: true,
        suggested_tolerance_multiplier: 100.0,
        use_relative_tolerances: true,
    });

    let ladder = config.build_retry_ladder(&base_ladder, &analysis);
    // Should include tolerance for size difference
    assert!(ladder.len() >= base_ladder.len());
}

// ── Comprehensive Analysis Tests ────────────────────────────────────────────────

#[test]
fn extreme_geometry_analysis_default() {
    let analysis = ExtremeGeometryAnalysis::default();
    assert!(analysis.near_tangent_configs.is_empty());
    assert!(analysis.high_aspect_ratio_edges.is_empty());
    assert!(analysis.degenerate_geometry.is_empty());
    assert!(analysis.size_difference.is_none());
    assert!(!analysis.has_extreme_geometry);
    assert!(analysis.issues_summary.is_empty());
}

#[test]
fn extreme_geometry_analysis_options_default() {
    let options = ExtremeGeometryAnalysisOptions::default();
    assert!(options.tolerance > 0.0);
    assert!(options.check_near_tangent);
    assert!(options.check_aspect_ratio);
    assert!(options.check_degenerate);
    assert!(options.check_size_difference);
}

#[test]
fn analyze_extreme_geometry_empty_brep() {
    let brep = BRep::default();
    let options = ExtremeGeometryAnalysisOptions::default();

    let analysis = analyze_extreme_geometry(&brep, None::<&BRep>, &options);
    // Empty BRep should not have extreme geometry issues
    assert!(!analysis.has_extreme_geometry);
}

// ── Integration with Adaptive Tolerance Tests ────────────────────────────────────

#[test]
fn adaptive_tolerance_integration() {
    let tol = AdaptiveTolerance::from_scale(1000.0);

    let near_tangent_handler = NearTangentHandler::from_adaptive(tol);
    let aspect_handler = AspectRatioAdaptiveTolerance::from_adaptive(tol);
    let degenerate_handler = DegenerateGeometryHandler::from_adaptive(tol);
    let size_handler = SizeDifferenceHandler::from_adaptive(tol);

    // All handlers should use the adaptive tolerance
    assert!((near_tangent_handler.base_tolerance - tol.coincidence()).abs() < 1e-15);
    assert!((aspect_handler.base_tolerance - tol.coincidence()).abs() < 1e-15);
    assert!((degenerate_handler.zero_tolerance - tol.coincidence()).abs() < 1e-15);
    assert!((size_handler.base_tolerance - tol.coincidence()).abs() < 1e-15);
}

#[test]
fn tolerance_level_consistency() {
    let tol = AdaptiveTolerance::from_scale(10.0);

    // Verify tolerance levels are consistent
    assert!(tol.tolerance(ToleranceLevel::Strict) < tol.tolerance(ToleranceLevel::Normal));
    assert!(tol.tolerance(ToleranceLevel::Normal) < tol.tolerance(ToleranceLevel::Relaxed));
    assert!(tol.tolerance(ToleranceLevel::Relaxed) < tol.tolerance(ToleranceLevel::Coarse));
}
