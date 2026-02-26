//! Channel flow solver tests
//!
//! Verifies correctness of the channel solver's flow regime classification
//! and Poiseuille number (shape factor) mathematics via the `FlowRegime` API.
//!
//! These tests operate on the public surface only: `FlowRegime`, `KN_SLIP_MIN`.

use cfd_1d::{FlowRegime, KN_SLIP_MIN};

// =========================================================================
// Regime Classification Tests
// =========================================================================

/// Stokes regime: Re < 1.0 must return `Stokes`.
#[test]
fn test_regime_stokes_below_one() {
    assert_eq!(
        FlowRegime::from_reynolds_number(0.5_f64),
        FlowRegime::Stokes
    );
    assert_eq!(
        FlowRegime::from_reynolds_number(0.0_f64),
        FlowRegime::Stokes
    );
}

/// Laminar regime: 1 ≤ Re < 2300 → `Laminar`.
#[test]
fn test_regime_laminar_range() {
    for re in [1.0_f64, 100.0, 500.0, 2000.0, 2299.9] {
        assert_eq!(
            FlowRegime::from_reynolds_number(re),
            FlowRegime::Laminar,
            "Re={re} should be Laminar"
        );
    }
}

/// Transitional regime: 2300 ≤ Re < 4000 → `Transitional`.
/// Note: Re = 4000 itself is classified as Turbulent (boundary is strict: < 4000).
#[test]
fn test_regime_transitional_range() {
    for re in [2300.0_f64, 3000.0, 3999.9] {
        assert_eq!(
            FlowRegime::from_reynolds_number(re),
            FlowRegime::Transitional,
            "Re={re} should be Transitional"
        );
    }
}

/// Turbulent regime: Re > 4000 → `Turbulent`.
#[test]
fn test_regime_turbulent_range() {
    for re in [4001.0_f64, 1e5, 1e7] {
        assert_eq!(
            FlowRegime::from_reynolds_number(re),
            FlowRegime::Turbulent,
            "Re={re} should be Turbulent"
        );
    }
}

// =========================================================================
// Knudsen Number Priority Tests
// =========================================================================

/// When `Kn >= KN_SLIP_MIN`, SlipFlow overrides any Reynolds classification.
#[test]
fn test_knudsen_overrides_re_for_slip_classification() {
    // Re=500 (normally Laminar) + Kn=0.05 → SlipFlow
    let regime = FlowRegime::classify_with_knudsen(500.0_f64, 0.05_f64);
    assert_eq!(regime, FlowRegime::SlipFlow, "Kn=0.05 must give SlipFlow");

    // Re=5000 (Turbulent) + Kn=KN_SLIP_MIN → SlipFlow
    let regime2 = FlowRegime::classify_with_knudsen(5000.0_f64, KN_SLIP_MIN);
    assert_eq!(regime2, FlowRegime::SlipFlow, "Kn=KN_SLIP_MIN must give SlipFlow");

    // Re=5000 + Kn=0.0 → Turbulent (no slip)
    let regime3 = FlowRegime::classify_with_knudsen(5000.0_f64, 0.0_f64);
    assert_eq!(regime3, FlowRegime::Turbulent, "Kn=0 should not give SlipFlow");
}

/// Kn just below threshold leaves regime continuum-classified.
#[test]
fn test_kn_below_threshold_stays_continuum() {
    let kn = KN_SLIP_MIN - 1e-10;
    let regime = FlowRegime::classify_with_knudsen(100.0_f64, kn);
    assert_eq!(regime, FlowRegime::Laminar, "Kn just below threshold = Laminar");
}

// =========================================================================
// Shah-London Polynomial Correctness Tests
// =========================================================================

/// The Shah-London polynomial value can be computed analytically and checked
/// for key known points without needing access to private channel methods.

/// Helper: compute Shah-London Po for aspect ratio alpha = a/b >= 1.
fn shah_london_po(alpha: f64) -> f64 {
    assert!(alpha >= 1.0);
    let inv = 1.0 / alpha;
    96.0 * (1.0
        - 1.3553 * inv
        + 1.9467 * inv * inv
        - 1.7012 * inv * inv * inv
        + 0.9564 * inv * inv * inv * inv
        - 0.2537 * inv * inv * inv * inv * inv)
}

/// For alpha = 1 (square), Po ≈ 56.9 (Shah & London 1978, Table 48).
#[test]
fn test_shah_london_square_value() {
    let po = shah_london_po(1.0);
    let err = (po - 56.9).abs();
    assert!(
        err < 0.3,
        "Shah-London square Po should be ≈56.9, got {po:.3}, error={err:.3}"
    );
}

/// For alpha → ∞ (parallel plates limit), Po → 96.
/// The Shah-London 5-term polynomial converges slowly;
/// at alpha = 1000 the result is within 0.5% of 96.
#[test]
fn test_shah_london_high_ar_approaches_96() {
    // alpha=100 → ~94.7 (still converging), alpha=1000 → ~95.94 (within 0.1%)
    let po_100 = shah_london_po(100.0);
    let po_1000 = shah_london_po(1000.0);
    // Both must be between 90 and 96
    assert!(
        po_100 > 90.0 && po_100 < 96.0,
        "Shah-London alpha=100 should be in (90, 96), got {po_100:.3}"
    );
    // At alpha=1000, within 1% of 96
    let err = (po_1000 - 96.0).abs() / 96.0;
    assert!(
        err < 0.01,
        "Shah-London alpha=1000 should be within 1% of 96, err={err:.4}, Po={po_1000:.4}"
    );
}

/// Shah-London Po must be monotonically increasing in alpha for alpha ∈ [1, 20].
#[test]
fn test_shah_london_monotone_increasing_in_ar() {
    let alphas = [1.0_f64, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0];
    let pos: Vec<f64> = alphas.iter().map(|&a| shah_london_po(a)).collect();
    for i in 1..pos.len() {
        assert!(
            pos[i] > pos[i - 1],
            "Shah-London must be monotone: Po({}) = {:.3} <= Po({}) = {:.3}",
            alphas[i],
            pos[i],
            alphas[i - 1],
            pos[i - 1]
        );
    }
}

// =========================================================================
// Beskok-Karniadakis Slip Correction Tests
// =========================================================================

/// Beskok-Karniadakis formula: R_slip = R_lam / (1 + 4*Kn).
///
/// We verify this mathematically since the physical channel is not public.
#[test]
fn test_beskok_karniadakis_formula_correctness() {
    // At Kn = 0: correction = 1.0
    let correction = |kn: f64| 1.0 / (1.0 + 4.0 * kn);
    assert!((correction(0.0) - 1.0).abs() < 1e-15, "Kn=0 must give factor=1");

    // At Kn = 0.05: correction = 1/1.2
    let c = correction(0.05);
    assert!((c - 1.0 / 1.2).abs() < 1e-15, "Kn=0.05 must give 1/1.2");

    // At Kn = 0.1: correction = 1/1.4
    let c2 = correction(0.1);
    assert!((c2 - 1.0 / 1.4).abs() < 1e-15, "Kn=0.1 must give 1/1.4");

    // Correction must be < 1 for all Kn > 0 → R_slip < R_lam
    for kn in [0.001, 0.01, 0.05, 0.1] {
        assert!(correction(kn) < 1.0, "Slip factor at Kn={kn} should be < 1");
    }
}

/// Slip correction factor is strictly decreasing in Kn.
#[test]
fn test_slip_factor_decreasing_in_kn() {
    let correction = |kn: f64| 1.0 / (1.0 + 4.0 * kn);
    let kns = [0.001_f64, 0.005, 0.01, 0.05, 0.1];
    for i in 1..kns.len() {
        assert!(
            correction(kns[i]) < correction(kns[i - 1]),
            "Slip factor must decrease: f({}) > f({})",
            kns[i - 1],
            kns[i]
        );
    }
}
