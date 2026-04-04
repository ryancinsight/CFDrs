//! Analysis module tests — pressure, resistance, flow, and performance metrics.
//!
//! Tests verify the correctness of accumulation invariants, series/parallel
//! resistance laws, flow rate aggregation, and blood shear safety limits.

use approx::assert_relative_eq;
use cfd_1d::solver::analysis::{
    BloodShearLimits, FlowAnalysis, PressureAnalysis, ResistanceAnalysis,
};
use std::collections::HashMap;

// ============================================================================
// PressureAnalysis: max/min tracking invariants
// ============================================================================

/// Pressure extrema tracking must be exact.
#[test]
fn test_pressure_max_min_tracking() {
    let mut analysis = PressureAnalysis::<f64>::new();
    analysis.add_pressure("n1".into(), 300.0);
    analysis.add_pressure("n2".into(), 100.0);
    analysis.add_pressure("n3".into(), 200.0);

    assert_relative_eq!(analysis.max_pressure, 300.0, epsilon = 1e-12);
    assert_relative_eq!(analysis.min_pressure, 100.0, epsilon = 1e-12);
}

/// Pressure range must equal max − min.
#[test]
fn test_pressure_range_equals_max_minus_min() {
    let mut analysis = PressureAnalysis::<f64>::new();
    analysis.add_pressure("a".into(), 500.0);
    analysis.add_pressure("b".into(), 200.0);

    assert_relative_eq!(analysis.pressure_range(), 300.0, epsilon = 1e-12);
}

/// Empty pressure analysis: range must return 0.
#[test]
fn test_pressure_range_empty_returns_zero() {
    let analysis = PressureAnalysis::<f64>::new();
    assert_eq!(analysis.pressure_range(), 0.0);
}

/// Average pressure must equal arithmetic mean of inputs.
#[test]
fn test_pressure_average_correct() {
    let mut analysis = PressureAnalysis::<f64>::new();
    analysis.add_pressure("a".into(), 100.0);
    analysis.add_pressure("b".into(), 200.0);
    analysis.add_pressure("c".into(), 300.0);

    assert_relative_eq!(analysis.average_pressure(), 200.0, epsilon = 1e-10);
}

// ============================================================================
// ResistanceAnalysis: series/parallel duality theorem
// ============================================================================

/// `series_resistance` must equal sum and be ≥ max individual resistance.
#[test]
fn test_resistance_series_greater_than_max() {
    let mut analysis = ResistanceAnalysis::<f64>::new();
    analysis.add_resistance("r1".into(), 100.0);
    analysis.add_resistance("r2".into(), 200.0);
    analysis.add_resistance("r3".into(), 300.0);

    let series = analysis.series_resistance();
    let max_r = 300.0_f64;

    assert_relative_eq!(series, 600.0, epsilon = 1e-10);
    assert!(series >= max_r, "R_series must be ≥ max(R_i)");
}

/// `parallel_resistance` must equal 1/Σ(1/Rᵢ) and be ≤ min individual resistance.
#[test]
fn test_resistance_parallel_less_than_min() {
    let mut analysis = ResistanceAnalysis::<f64>::new();
    analysis.add_resistance("r1".into(), 100.0);
    analysis.add_resistance("r2".into(), 200.0);
    analysis.add_resistance("r3".into(), 300.0);

    let parallel = analysis.parallel_resistance();
    let expected = 1.0 / (1.0 / 100.0 + 1.0 / 200.0 + 1.0 / 300.0);
    let min_r = 100.0_f64;

    assert_relative_eq!(parallel, expected, max_relative = 1e-10);
    assert!(parallel <= min_r, "R_parallel must be ≤ min(R_i)");
}

/// Series/Parallel duality bound: R_parallel ≤ min ≤ max ≤ R_series.
#[test]
fn test_resistance_duality_bound() {
    let mut analysis = ResistanceAnalysis::<f64>::new();
    for (i, &r) in [50.0, 150.0, 400.0, 800.0].iter().enumerate() {
        analysis.add_resistance(format!("r{i}"), r);
    }

    let r_par = analysis.parallel_resistance();
    let (_, r_min) = analysis.min_resistance().expect("test invariant");
    let (_, r_max) = analysis.max_resistance().expect("test invariant");
    let r_ser = analysis.series_resistance();

    assert!(
        r_par <= r_min,
        "R_parallel ({}) must be ≤ min(R_i) ({})",
        r_par,
        r_min
    );
    assert!(r_min <= r_max, "min(R_i) must be ≤ max(R_i)");
    assert!(
        r_max <= r_ser,
        "max(R_i) ({}) must be ≤ R_series ({})",
        r_max,
        r_ser
    );
}

/// `add_resistance` must accumulate `total_resistance` correctly.
#[test]
fn test_resistance_total_accumulates_on_add() {
    let mut analysis = ResistanceAnalysis::<f64>::new();
    analysis.add_resistance("a".into(), 10.0);
    analysis.add_resistance("b".into(), 20.0);
    assert_relative_eq!(analysis.total_resistance, 30.0, epsilon = 1e-12);
}

// ============================================================================
// FlowAnalysis: aggregation and blood shear limits
// ============================================================================

/// Average flow rate must equal arithmetic mean.
#[test]
fn test_flow_average_correct() {
    let mut analysis = FlowAnalysis::<f64>::new();
    analysis.add_component_flow("c1".into(), 1e-9_f64);
    analysis.add_component_flow("c2".into(), 3e-9_f64);

    assert_relative_eq!(analysis.average_flow_rate(), 2e-9_f64, epsilon = 1e-21);
}

/// `flag_fda_shear_limit_violations` must flag components exceeding stress limit.
#[test]
fn test_blood_shear_violation_stress_limit() {
    let mut analysis = FlowAnalysis::<f64>::new();
    analysis.add_wall_shear_stress("edge_1".into(), 180.0_f64);
    analysis.add_wall_shear_rate("edge_1".into(), 20_000.0);

    let limits = BloodShearLimits {
        max_wall_shear_stress_pa: 150.0,
        max_wall_shear_rate_per_s: None,
        max_giersiepen_hi: None,
        max_taskin_hi: None,
    };

    let violations = analysis.flag_fda_shear_limit_violations(&limits);
    assert_eq!(violations.len(), 1);
    assert_eq!(violations[0].component_id, "edge_1");
    assert!(violations[0].stress_exceedance_ratio > 1.0);
}

/// `flag_fda_shear_limit_violations` must flag components exceeding the optional rate limit.
#[test]
fn test_blood_shear_violation_rate_limit() {
    let mut analysis = FlowAnalysis::<f64>::new();
    analysis.add_wall_shear_stress("edge_2".into(), 80.0_f64);
    analysis.add_wall_shear_rate("edge_2".into(), 60_000.0);

    let limits = BloodShearLimits {
        max_wall_shear_stress_pa: 150.0,
        max_wall_shear_rate_per_s: Some(40_000.0),
        max_giersiepen_hi: None,
        max_taskin_hi: None,
    };

    let violations = analysis.flag_fda_shear_limit_violations(&limits);
    assert_eq!(violations.len(), 1, "Rate-limit violation must be flagged");
}

/// Components within both stress and rate limits must produce zero violations.
#[test]
fn test_blood_shear_no_violation_within_limits() {
    let mut analysis = FlowAnalysis::<f64>::new();
    analysis.add_wall_shear_stress("safe_edge".into(), 50.0_f64);
    analysis.add_wall_shear_rate("safe_edge".into(), 10_000.0);

    let limits = BloodShearLimits {
        max_wall_shear_stress_pa: 150.0,
        max_wall_shear_rate_per_s: Some(40_000.0),
        max_giersiepen_hi: None,
        max_taskin_hi: None,
    };

    let violations = analysis.flag_fda_shear_limit_violations(&limits);
    assert!(
        violations.is_empty(),
        "Safe conditions must produce zero violations"
    );
}

/// `flag_hemolysis_limit_violations` must account for exposure duration.
#[test]
fn test_blood_damage_violation_respects_residence_time() {
    let mut analysis = FlowAnalysis::<f64>::new();
    analysis.add_wall_shear_stress("edge_damage".into(), 180.0_f64);

    let limits = BloodShearLimits::<f64>::fda_conservative_whole_blood()
        .with_hemolysis_limits(Some(1e-3), Some(5e-3));
    let residence_times = HashMap::from([("edge_damage".to_string(), 0.5_f64)]);

    let violations = analysis.flag_hemolysis_limit_violations(&limits, &residence_times);
    assert_eq!(violations.len(), 1);
    assert_eq!(violations[0].component_id, "edge_damage");
    assert!(violations[0]
        .giersiepen_exceedance_ratio
        .is_some_and(|ratio| ratio > 1.0));
    assert!(violations[0]
        .taskin_exceedance_ratio
        .is_some_and(|ratio| ratio > 1.0));
}
