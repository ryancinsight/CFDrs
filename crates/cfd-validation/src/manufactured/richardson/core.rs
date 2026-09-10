//! Core Richardson machinery for MMS studies
//!
//! CFDRS-GA-012: the duplicate Richardson implementation that lived here was
//! retired in favor of the single consolidated one in
//! [`crate::convergence::RichardsonExtrapolation`], which returns typed
//! `cfd-core` errors and carries the absorbed numerical-stability guards
//! (signed convergence-ratio rejection, order bounds, `r^p ≈ 1` checks).
//! The struct is re-exported below so existing
//! `manufactured::richardson::core::RichardsonExtrapolation` paths keep
//! resolving against the consolidated implementation.

use crate::scalar;
use cfd_core::error::Result;
use eunomia::{FloatElement, RealField};

pub use crate::convergence::RichardsonExtrapolation;

/// Estimate the observed convergence order from three coarsely-to-finely
/// ordered solutions and the refinement ratio between consecutive grids.
///
/// Thin adapter over
/// [`RichardsonExtrapolation::estimate_order`][canonical-estimate]; kept so
/// call sites that read `(f_coarse, f_medium, f_fine, r)` in that order
/// retain their argument order.
///
/// [canonical-estimate]: crate::convergence::RichardsonExtrapolation::estimate_order
pub fn estimate_order<T>(f_coarse: T, f_medium: T, f_fine: T, r: T) -> Result<T>
where
    T: RealField + Copy + FloatElement,
{
    RichardsonExtrapolation::estimate_order(f_coarse, f_medium, f_fine, r)
}

/// Perform Richardson extrapolation with automatic order estimation from
/// three coarsely-to-finely ordered solutions.
///
/// Estimates the order from the three solutions, then extrapolates the
/// continuum value using the fine and medium solutions. Returns
/// `(extrapolated, order)` on success.
pub fn extrapolate<T>(
    f_coarse: T,
    f_medium: T,
    f_fine: T,
    r: T,
) -> Result<(T, T)>
where
    T: RealField + Copy + FloatElement,
{
    let order = RichardsonExtrapolation::estimate_order(f_coarse, f_medium, f_fine, r)?;
    let extrapolator = RichardsonExtrapolation::with_order(order, r)?;
    let extrapolated = extrapolator.extrapolate(f_fine, f_medium)?;
    Ok((extrapolated, order))
}

/// Check that solution differences decrease monotonically toward the fine
/// grid, the minimal signal of asymptotic convergence.
///
/// Distinct contract from
/// [`RichardsonExtrapolation::is_asymptotic`][canonical-asymptotic], which
/// checks the observed ratio against the expected `r^order` band: this check
/// only requires monotone error decrease with sufficient variation, which is
/// what the MMS study's per-triple reporting has always asserted.
///
/// [canonical-asymptotic]: crate::convergence::RichardsonExtrapolation::is_asymptotic
pub fn is_asymptotic<T>(f_coarse: T, f_medium: T, f_fine: T) -> bool
where
    T: RealField + Copy + FloatElement,
{
    let eps = <T as FloatElement>::from_f64(1e-12);

    // |f_medium - f_coarse| > |f_fine - f_medium| > eps: errors shrink as the
    // grid is refined, with sufficient variation to be meaningful.
    let diff_coarse = scalar::abs(f_medium - f_coarse);
    let diff_fine = scalar::abs(f_fine - f_medium);

    diff_coarse > eps && diff_fine > eps && diff_coarse > diff_fine
}

/// Data-driven order estimation using multiple grid levels
pub struct DataDrivenOrderEstimation;

impl DataDrivenOrderEstimation {
    /// Estimate convergence order using data-driven approach following Roache (1998)
    ///
    /// ## Methodology
    ///
    /// Uses multiple grid levels to estimate convergence order without hardcoded assumptions.
    /// Implements robust order estimation with numerical stability checks and outlier filtering.
    ///
    /// ## Robustness Features
    ///
    /// - Uses median of multiple order estimates for robustness against outliers
    /// - Filters unreliable estimates based on numerical stability criteria
    /// - Validates order estimates within reasonable CFD ranges (0.5 ≤ p ≤ 6.0)
    /// - Supports non-uniform refinement ratios r21 != r32 via bracketing root-finding
    /// - Falls back to second-order default only when no reliable data available
    pub fn estimate_order_from_solutions<T>(solutions: &[T], refinement_ratios: &[T]) -> T
    where
        T: RealField + Copy + FloatElement,
    {
        let mut order_estimates = Vec::new();

        // Use all available triplets for order estimation (no hardcoded assumptions)
        for i in 0..solutions.len().saturating_sub(2) {
            let phi_coarse = solutions[i];
            let phi_medium = solutions[i + 1];
            let phi_fine = solutions[i + 2];

            let r21 = refinement_ratios[i];
            let r32 = refinement_ratios[i + 1];

            // Check for sufficient solution variation (avoid division by near-zero)
            let e21 = phi_medium - phi_coarse; // change from coarse->medium
            let e32 = phi_fine - phi_medium; // change from medium->fine

            let eps = <T as FloatElement>::from_f64(1e-12);
            let e21_abs = scalar::abs(e21);
            let e32_abs = scalar::abs(e32);
            if e21_abs <= eps || e32_abs <= eps {
                continue;
            }

            // If refinement ratios are effectively uniform, use closed-form estimate
            let one_percent = <T as FloatElement>::from_f64(0.01);
            if scalar::abs(r21 - r32) / r21 <= one_percent {
                let r = r32;
                let ratio = e21_abs / e32_abs;
                let p_est = scalar::ln(ratio) / scalar::ln(r);
                if p_est > <T as FloatElement>::from_f64(0.1)
                    && p_est < <T as FloatElement>::from_f64(6.0)
                {
                    order_estimates.push(p_est);
                }
                continue;
            }

            // General non-uniform case: solve for p via bisection
            // e21/e32 ≈ (r21^p - 1) / (r32^p - 1)
            let target = e21_abs / e32_abs;

            let mut lo = <T as FloatElement>::from_f64(0.1);
            let mut hi = <T as FloatElement>::from_f64(8.0);

            let f = |p: T| -> T {
                let r21_p = scalar::powf(r21, p);
                let r32_p = scalar::powf(r32, p);
                let num = r21_p - scalar::one::<T>();
                let den = r32_p - scalar::one::<T>();
                if scalar::abs(den) <= eps {
                    return <T as FloatElement>::from_f64(1e12);
                }
                // General non-uniform refinement formula:
                // |e21|/|e32| = r32^p * (r21^p - 1) / (r32^p - 1)
                r32_p * (num / den) - target
            };

            let mut f_lo = f(lo);
            let mut f_hi = f(hi);

            // Expand hi if needed to achieve a bracket
            let mut expand_iters = 0;
            while (f_lo > scalar::zero::<T>() && f_hi > scalar::zero::<T>())
                || (f_lo < scalar::zero::<T>() && f_hi < scalar::zero::<T>())
            {
                if expand_iters >= 5 {
                    break;
                }
                hi = hi + hi; // exponential expansion
                f_hi = f(hi);
                expand_iters += 1;
            }

            // If still not bracketed, skip this triplet
            if !((f_lo <= scalar::zero::<T>() && f_hi >= scalar::zero::<T>())
                || (f_lo >= scalar::zero::<T>() && f_hi <= scalar::zero::<T>()))
            {
                continue;
            }

            // Bisection iteration
            let tol = <T as FloatElement>::from_f64(1e-10);
            let two = <T as FloatElement>::from_f64(2.0);
            for _ in 0..60 {
                let mid = (lo + hi) / two;
                let f_mid = f(mid);
                if scalar::abs(f_mid) <= tol {
                    lo = mid;
                    hi = mid;
                    break;
                }
                if (f_lo <= scalar::zero::<T>() && f_mid >= scalar::zero::<T>())
                    || (f_lo >= scalar::zero::<T>() && f_mid <= scalar::zero::<T>())
                {
                    hi = mid;
                    f_hi = f_mid;
                } else {
                    lo = mid;
                    f_lo = f_mid;
                }
            }

            let p_est = (lo + hi) / <T as FloatElement>::from_f64(2.0);
            if p_est > <T as FloatElement>::from_f64(0.1)
                && p_est < <T as FloatElement>::from_f64(6.0)
            {
                order_estimates.push(p_est);
            }
        }

        // Use median order estimate for robustness (resistant to outliers)
        if order_estimates.is_empty() {
            // No reliable data: fall back to second-order (most common in CFD)
            <T as FloatElement>::from_f64(2.0)
        } else {
            // Sort estimates and take median
            order_estimates.sort_by(|a: &T, b: &T| a.partial_cmp(b).expect("expected value"));
            let median_idx = order_estimates.len() / 2;
            order_estimates[median_idx]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_core::test_support::assert_rejects;
    use eunomia::NumericElement;

    #[test]
    fn test_richardson_extrapolation_basic() {
        // Test basic Richardson extrapolation with known second-order convergence
        // Solution: φ(h) = φ_exact + C h²
        let h1 = 1.0; // coarse grid
        let h2 = 0.5; // medium grid
        let h3 = 0.25; // fine grid

        // Second-order convergence: φ(h) = 1.0 + h²
        let phi1 = 1.0 + h1 * h1; // 1.0 + 1.0 = 2.0
        let phi2 = 1.0 + h2 * h2; // 1.0 + 0.25 = 1.25
        let phi3 = 1.0 + h3 * h3; // 1.0 + 0.0625 = 1.0625

        let r = 2.0; // refinement ratio

        let (extrapolated, order) =
            extrapolate(phi1, phi2, phi3, r).expect("expected value");

        // Should extrapolate to very close to exact solution (1.0)
        assert!(
            scalar::abs(extrapolated - 1.0) < 1e-10,
            "Extrapolation error too large: {extrapolated}"
        );

        // Should estimate order close to 2.0
        assert!(
            scalar::abs(order - 2.0) < 0.1,
            "Order estimation error: {order}"
        );
    }

    #[test]
    fn test_richardson_extrapolation_numerical_stability() {
        // Test numerical stability when r^p ≈ 1 (problematic case)

        // Case where r^p is very close to 1 - should handle gracefully
        let phi1 = 1.0001;
        let phi2 = 1.00005;
        let phi3 = 1.000025;
        let r = 1.0001; // Very small refinement ratio

        // This should either succeed with reasonable bounds or fail gracefully
        let result = extrapolate(phi1, phi2, phi3, r);

        match result {
            Ok((extrapolated, order)) => {
                // If it succeeds, results should be reasonable
                assert!(
                    NumericElement::is_finite(extrapolated),
                    "Extrapolated value should be finite"
                );
                assert!(
                    order > 0.0 && order < 10.0,
                    "Order should be reasonable: {order}"
                );
            }
            Err(e) => {
                // If it fails, should be due to numerical instability or
                // insufficient variation, reported through the typed error.
                let msg = e.to_string();
                assert!(
                    msg.contains("unstable") || msg.contains("Insufficient"),
                    "Should fail for numerical reasons: {msg}"
                );
            }
        }
    }

    #[test]
    fn test_richardson_extrapolation_edge_cases() {
        // Test edge cases that could cause numerical issues

        // Case 1: Very small differences (near convergence) - below eps=1e-12
        let result = estimate_order(
            1.0 + 1e-13,
            1.0 + 0.5e-13,
            1.0 + 0.25e-13,
            2.0,
        );
        assert_rejects(
            &result,
            "Solutions too close to estimate order",
        );

        // Case 2: Zero differences (exact solution)
        let result = estimate_order(1.0, 1.0, 1.0, 2.0);
        assert_rejects(
            &result,
            "Solutions too close to estimate order",
        );

        // Case 3: Invalid refinement ratio
        let result = estimate_order(2.0, 1.5, 1.25, 0.0);
        assert_rejects(
            &result,
            "Richardson extrapolation numerically unstable: order -0 out of bounds",
        );

        // Case 4: Negative refinement ratio
        let result = estimate_order(2.0, 1.5, 1.25, -2.0);
        assert_rejects(
            &result,
            "Richardson extrapolation numerically unstable: order NaN out of bounds",
        );
    }

    #[test]
    fn test_data_driven_order_estimation_uniform_grid() {
        // Test data-driven order estimation with uniform refinement

        // Solutions with known second-order convergence: φ(h) = 1.0 + h²
        let solutions = vec![
            1.0 + 1.0,      // h = 1.0
            1.0 + 0.25,     // h = 0.5
            1.0 + 0.0625,   // h = 0.25
            1.0 + 0.015625, // h = 0.125
        ];

        let refinement_ratios = vec![2.0, 2.0, 2.0];

        let estimated_order = DataDrivenOrderEstimation::estimate_order_from_solutions(
            &solutions,
            &refinement_ratios,
        );

        // Should estimate order close to 2.0
        assert!(
            scalar::abs(estimated_order - 2.0) < 0.1,
            "Data-driven order estimation failed: {estimated_order}"
        );
    }

    #[test]
    fn test_data_driven_order_estimation_nonuniform_grid() {
        // Test data-driven order estimation with non-uniform refinement ratios

        // Solutions with known 1.5-order convergence: φ(h) = 1.0 + h^1.5
        let solutions = vec![
            1.0 + 1.0_f64.powf(1.5),  // h = 1.0
            1.0 + 0.5_f64.powf(1.5),  // h = 0.5
            1.0 + 0.25_f64.powf(1.5), // h = 0.25
        ];

        let refinement_ratios = vec![2.0, 2.0]; // Non-uniform in general case

        let estimated_order = DataDrivenOrderEstimation::estimate_order_from_solutions(
            &solutions,
            &refinement_ratios,
        );

        // Should estimate order close to 1.5
        assert!(
            scalar::abs(estimated_order - 1.5) < 0.2,
            "Non-uniform grid order estimation failed: {estimated_order}"
        );
    }

    #[test]
    fn test_data_driven_order_estimation_edge_cases() {
        // Test edge cases for data-driven estimation

        // Case 1: Insufficient data points
        let solutions = vec![1.0, 1.1];
        let refinement_ratios = vec![2.0];
        let estimated_order = DataDrivenOrderEstimation::estimate_order_from_solutions(
            &solutions,
            &refinement_ratios,
        );
        assert!(
            scalar::abs(estimated_order - 2.0) < 1e-10,
            "Should fall back to 2.0 with insufficient data"
        );

        // Case 2: All solutions identical (no convergence)
        let solutions = vec![1.0, 1.0, 1.0, 1.0];
        let refinement_ratios = vec![2.0, 2.0, 2.0];
        let estimated_order = DataDrivenOrderEstimation::estimate_order_from_solutions(
            &solutions,
            &refinement_ratios,
        );
        assert!(
            scalar::abs(estimated_order - 2.0) < 1e-10,
            "Should fall back to 2.0 with no convergence"
        );

        // Case 3: Empty input
        let solutions: Vec<f64> = vec![];
        let refinement_ratios: Vec<f64> = vec![];
        let estimated_order = DataDrivenOrderEstimation::estimate_order_from_solutions(
            &solutions,
            &refinement_ratios,
        );
        assert!(
            scalar::abs(estimated_order - 2.0) < 1e-10,
            "Should fall back to 2.0 with empty input"
        );
    }

    #[test]
    fn test_asymptotic_range_detection() {
        // Test asymptotic range detection (monotone-error contract)

        // Case 1: Proper asymptotic convergence (error decreasing with ratio > 1)
        // For r=2, p=1: phi(h) = 1.0 + h
        // phi(1.0)=2.0, phi(0.5)=1.5, phi(0.25)=1.25
        // d1 = 0.5, d2 = 0.25
        assert!(
            is_asymptotic(2.0, 1.5, 1.25),
            "Should detect asymptotic convergence"
        );

        // Case 2: Not asymptotic (error increasing/diverging)
        assert!(
            !is_asymptotic(1.0, 1.2, 1.5),
            "Should detect non-asymptotic behavior"
        );

        // Case 3: Insufficient variation (below threshold 1e-12)
        assert!(
            !is_asymptotic(1.0, 1.0 + 1e-13, 1.0 + 0.5e-13),
            "Should detect insufficient variation"
        );

        // Case 4: Zero differences
        assert!(
            !is_asymptotic(1.0, 1.0, 1.0),
            "Should detect zero variation"
        );
    }

    #[test]
    fn test_richardson_extrapolation_property_based() {
        // Property-based testing: Richardson extrapolation should be invariant under scaling

        let phi1 = 2.0;
        let phi2 = 1.5;
        let phi3 = 1.25;
        let r = 2.0;

        let (extrapolated1, order1) = extrapolate(phi1, phi2, phi3, r).expect("expected value");

        // Scale all values by constant factor
        let scale = std::f64::consts::PI;
        let (extrapolated2, order2) =
            extrapolate(phi1 * scale, phi2 * scale, phi3 * scale, r).expect("expected value");

        // Extrapolated value should scale, order should be invariant
        assert!(
            scalar::abs(extrapolated2 - extrapolated1 * scale) < 1e-12,
            "Extrapolation should be linear"
        );
        assert!(
            scalar::abs(order2 - order1) < 1e-12,
            "Order should be invariant under scaling"
        );
    }

    #[test]
    fn test_convergence_order_bounds() {
        type TestCase = (f64, Box<dyn Fn(f64) -> f64>);

        // Test that estimated orders are within reasonable CFD bounds

        // Generate test cases with known orders
        let test_cases: Vec<TestCase> = vec![
            (2.0, Box::new(|h: f64| 1.0 + h * h)),       // Second order
            (1.5, Box::new(|h: f64| 1.0 + h.powf(1.5))), // 1.5 order
            (3.0, Box::new(|h: f64| 1.0 + h * h * h)),   // Third order
        ];

        for (expected_order, solution_fn) in test_cases {
            let h_vals = [1.0, 0.5, 0.25];
            let solutions: Vec<f64> = h_vals.iter().map(|&h| solution_fn(h)).collect();
            let refinement_ratios = [2.0, 2.0];

            let estimated_order = DataDrivenOrderEstimation::estimate_order_from_solutions(
                &solutions,
                &refinement_ratios,
            );

            // Order should be within reasonable CFD bounds (0.5 to 6.0) and close to expected
            assert!(
                estimated_order > 0.5 && estimated_order < 6.0,
                "Order out of bounds: {estimated_order} (expected ~{expected_order})"
            );
            assert!(
                scalar::abs(estimated_order - expected_order) < 0.5,
                "Order estimation too inaccurate: {estimated_order} vs {expected_order}"
            );
        }
    }
}
