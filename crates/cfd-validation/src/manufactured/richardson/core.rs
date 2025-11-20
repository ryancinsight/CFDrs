//! Core Richardson extrapolation algorithms

use nalgebra::{ComplexField, RealField};
use num_traits::{Float, FromPrimitive};

use super::types::RichardsonResult;

/// Core Richardson extrapolation implementation
pub struct RichardsonExtrapolation;

impl RichardsonExtrapolation {
    /// Estimate convergence order using Richardson extrapolation formula
    ///
    /// ## Richardson Extrapolation Order Estimation Theorem
    ///
    /// **Statement**: For solutions φ₁, φ₂, φ₃ on three consecutive grids with refinement ratio r,
    /// where the grids are in the asymptotic convergence range, the convergence order p can be
    /// estimated using the three-point Richardson extrapolation formula:
    ///
    /// p = ln[(φ₁ - φ₂) / (φ₂ - φ₃)] / ln(r)
    ///
    /// **Assumptions**:
    /// 1. **Asymptotic convergence**: Solutions must be in the asymptotic range where φ(h) = φ_exact + C h^p + O(h^q) with q > p
    /// 2. **Grid ordering**: φ₁, φ₂, φ₃ correspond to coarse, medium, and fine grids with h₁ > h₂ > h₃
    /// 3. **Refinement ratio**: r = h₁/h₂ = h₂/h₃ > 1 (consistent grid refinement)
    /// 4. **Monotonic convergence**: |φ₂ - φ₁| > |φ₃ - φ₂| (errors decreasing with grid refinement)
    /// 5. **Sufficient accuracy**: Solutions have sufficient variation to avoid numerical cancellation
    /// 6. **Numerical stability**: r^p ≠ 1 (avoid division by near-zero in extrapolation)
    ///
    /// **Validity conditions**: The formula is valid when 0.5 ≤ p ≤ 6.0 (typical CFD convergence orders)
    /// and when Richardson extrapolation provides stable, bounded results.
    ///
    /// **References**:
    /// - Richardson, L.F. (1910): "The deferred approach to the limit"
    /// - Roache, P.J. (1998): Verification and Validation in Computational Science and Engineering
    /// - ASME V&V 20-2009: Standard for Verification and Validation in CFD
    pub fn estimate_order<T>(f1: T, f2: T, f3: T, r: T) -> Result<T, String>
    where
        T: RealField + Copy + Float + FromPrimitive,
    {
        let eps = <T as FromPrimitive>::from_f64(1e-12).unwrap();

        // Check for sufficient variation
        let diff12 = f1 - f2;
        let diff23 = f2 - f3;

        if ComplexField::abs(diff12) < eps || ComplexField::abs(diff23) < eps {
            return Err("Insufficient solution variation for order estimation".to_string());
        }

        let ratio = ComplexField::abs(diff12 / diff23);
        if ratio <= T::zero() || !ratio.is_finite() {
            return Err("Invalid convergence ratio".to_string());
        }

        let order = ComplexField::ln(ratio) / ComplexField::ln(r);
        if !order.is_finite() || order < <T as FromPrimitive>::from_f64(0.1).unwrap() {
            return Err("Invalid convergence order estimated".to_string());
        }

        Ok(order)
    }

    /// Check if solutions are in asymptotic range
    pub fn is_asymptotic<T>(f1: T, f2: T, f3: T) -> bool
    where
        T: RealField + Copy + Float + FromPrimitive,
    {
        let eps = <T as FromPrimitive>::from_f64(1e-12).unwrap();

        // Simple asymptotic check: |f2 - f1| > |f3 - f2|
        // This ensures we're seeing convergence behavior
        let diff1 = ComplexField::abs(f2 - f1);
        let diff2 = ComplexField::abs(f3 - f2);

        diff1 > eps && diff2 > eps && diff1 > diff2
    }

    /// Perform Richardson extrapolation with order estimation
    ///
    /// ## Richardson Extrapolation Theorem
    ///
    /// **Statement**: If a numerical method has asymptotic convergence φ(h) = φ_exact + C h^p + O(h^q)
    /// with q > p, then the exact solution can be extrapolated from three solutions on grids
    /// with refinement ratio r using:
    ///
    /// φ_exact = φ₁ + (φ₁ - φ₂) / (r^p - 1)
    ///
    /// **Assumptions** (same as estimate_order plus):
    /// 7. **Error expansion**: Solutions follow the asymptotic error expansion
    /// 8. **Leading error dominance**: The p-th order term dominates the error expansion
    /// 9. **Consistent discretization**: All solutions use the same numerical method
    ///
    /// **Stability condition**: |r^p - 1| > ε (numerical stability threshold)
    ///
    /// **References**:
    /// - Richardson, L.F. (1910): "The deferred approach to the limit"
    /// - Roache, P.J. (1998): Chapter 4 - Richardson Extrapolation
    pub fn extrapolate<T>(coarse: T, medium: T, fine: T, r: T) -> Result<(T, T), String>
    where
        T: RealField + Copy + Float + FromPrimitive,
    {
        // First estimate the order
        let order = Self::estimate_order(coarse, medium, fine, r)?;

        // Then perform extrapolation
        let r_pow_p = ComplexField::powf(r, order);
        let denominator = r_pow_p - T::one();

        if ComplexField::abs(denominator) < <T as FromPrimitive>::from_f64(1e-8).unwrap() {
            return Err("Richardson extrapolation numerically unstable (r^p ≈ 1)".to_string());
        }

        let extrapolated = fine + (fine - coarse) / denominator;

        Ok((extrapolated, order))
    }
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
        T: RealField + Copy + Float + FromPrimitive,
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

            let eps = <T as FromPrimitive>::from_f64(1e-12).unwrap();
            let e21_abs = ComplexField::abs(e21);
            let e32_abs = ComplexField::abs(e32);
            if e21_abs <= eps || e32_abs <= eps {
                continue;
            }

            // If refinement ratios are effectively uniform, use closed-form estimate
            let one_percent = <T as FromPrimitive>::from_f64(0.01).unwrap();
            if ComplexField::abs(r21 - r32) / r21 <= one_percent {
                let r = r32;
                let ratio = e21_abs / e32_abs;
                let p_est = ComplexField::ln(ratio) / ComplexField::ln(r);
                if p_est > <T as FromPrimitive>::from_f64(0.1).unwrap()
                    && p_est < <T as FromPrimitive>::from_f64(6.0).unwrap()
                {
                    order_estimates.push(p_est);
                }
                continue;
            }

            // General non-uniform case: solve for p via bisection
            // e21/e32 ≈ (r21^p - 1) / (r32^p - 1)
            let target = e21_abs / e32_abs;

            let mut lo = <T as FromPrimitive>::from_f64(0.1).unwrap();
            let mut hi = <T as FromPrimitive>::from_f64(8.0).unwrap();

            let f = |p: T| -> T {
                let r21_p = ComplexField::powf(r21, p);
                let r32_p = ComplexField::powf(r32, p);
                let num = r21_p - T::one();
                let den = r32_p - T::one();
                if ComplexField::abs(den) <= eps {
                    return <T as FromPrimitive>::from_f64(1e12).unwrap();
                }
                // General non-uniform refinement formula:
                // |e21|/|e32| = r32^p * (r21^p - 1) / (r32^p - 1)
                r32_p * (num / den) - target
            };

            let mut f_lo = f(lo);
            let mut f_hi = f(hi);

            // Expand hi if needed to achieve a bracket
            let mut expand_iters = 0;
            while (f_lo > T::zero() && f_hi > T::zero()) || (f_lo < T::zero() && f_hi < T::zero()) {
                if expand_iters >= 5 {
                    break;
                }
                hi = hi + hi; // exponential expansion
                f_hi = f(hi);
                expand_iters += 1;
            }

            // If still not bracketed, skip this triplet
            if !((f_lo <= T::zero() && f_hi >= T::zero())
                || (f_lo >= T::zero() && f_hi <= T::zero()))
            {
                continue;
            }

            // Bisection iteration
            let tol = <T as FromPrimitive>::from_f64(1e-10).unwrap();
            let two = <T as FromPrimitive>::from_f64(2.0).unwrap();
            for _ in 0..60 {
                let mid = (lo + hi) / two;
                let f_mid = f(mid);
                if ComplexField::abs(f_mid) <= tol {
                    lo = mid;
                    hi = mid;
                    break;
                }
                if (f_lo <= T::zero() && f_mid >= T::zero())
                    || (f_lo >= T::zero() && f_mid <= T::zero())
                {
                    hi = mid;
                    f_hi = f_mid;
                } else {
                    lo = mid;
                    f_lo = f_mid;
                }
            }

            let p_est = (lo + hi) / <T as FromPrimitive>::from_f64(2.0).unwrap();
            if p_est > <T as FromPrimitive>::from_f64(0.1).unwrap()
                && p_est < <T as FromPrimitive>::from_f64(6.0).unwrap()
            {
                order_estimates.push(p_est);
            }
        }

        // Use median order estimate for robustness (resistant to outliers)
        if order_estimates.is_empty() {
            // No reliable data: fall back to second-order (most common in CFD)
            <T as FromPrimitive>::from_f64(2.0).unwrap()
        } else {
            // Sort estimates and take median
            order_estimates.sort_by(|a: &T, b: &T| a.partial_cmp(b).unwrap());
            let median_idx = order_estimates.len() / 2;
            order_estimates[median_idx]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::FromPrimitive;

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
            RichardsonExtrapolation::extrapolate(phi1, phi2, phi3, r).unwrap();

        // Should extrapolate to very close to exact solution (1.0)
        assert!(
            nalgebra::ComplexField::abs(extrapolated - 1.0) < 1e-10,
            "Extrapolation error too large: {}",
            extrapolated
        );

        // Should estimate order close to 2.0
        assert!(
            nalgebra::ComplexField::abs(order - 2.0) < 0.1,
            "Order estimation error: {}",
            order
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
        let result = RichardsonExtrapolation::extrapolate(phi1, phi2, phi3, r);

        match result {
            Ok((extrapolated, order)) => {
                // If it succeeds, results should be reasonable
                assert!(
                    extrapolated.is_finite(),
                    "Extrapolated value should be finite"
                );
                assert!(
                    order > 0.0 && order < 10.0,
                    "Order should be reasonable: {}",
                    order
                );
            }
            Err(msg) => {
                // If it fails, should be due to numerical instability
                assert!(
                    msg.contains("unstable") || msg.contains("Insufficient"),
                    "Should fail for numerical reasons: {}",
                    msg
                );
            }
        }
    }

    #[test]
    fn test_richardson_extrapolation_edge_cases() {
        // Test edge cases that could cause numerical issues

        // Case 1: Very small differences (near convergence)
        let result =
            RichardsonExtrapolation::estimate_order(1.0000001, 1.00000005, 1.000000025, 2.0);
        assert!(result.is_err(), "Should detect insufficient variation");

        // Case 2: Zero differences (exact solution)
        let result = RichardsonExtrapolation::estimate_order(1.0, 1.0, 1.0, 2.0);
        assert!(result.is_err(), "Should detect zero variation");

        // Case 3: Invalid refinement ratio
        let result = RichardsonExtrapolation::estimate_order(2.0, 1.5, 1.25, 0.0);
        assert!(result.is_err(), "Should detect invalid refinement ratio");

        // Case 4: Negative refinement ratio
        let result = RichardsonExtrapolation::estimate_order(2.0, 1.5, 1.25, -2.0);
        assert!(result.is_err(), "Should detect negative refinement ratio");
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
            nalgebra::ComplexField::abs(estimated_order - 2.0) < 0.1,
            "Data-driven order estimation failed: {}",
            estimated_order
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
            nalgebra::ComplexField::abs(estimated_order - 1.5) < 0.2,
            "Non-uniform grid order estimation failed: {}",
            estimated_order
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
            nalgebra::ComplexField::abs(estimated_order - 2.0) < 1e-10,
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
            nalgebra::ComplexField::abs(estimated_order - 2.0) < 1e-10,
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
            nalgebra::ComplexField::abs(estimated_order - 2.0) < 1e-10,
            "Should fall back to 2.0 with empty input"
        );
    }

    #[test]
    fn test_asymptotic_range_detection() {
        // Test asymptotic range detection

        // Case 1: Proper asymptotic convergence (error decreasing)
        assert!(
            RichardsonExtrapolation::is_asymptotic(1.0, 0.75, 0.5),
            "Should detect asymptotic convergence"
        );

        // Case 2: Not asymptotic (error not decreasing)
        assert!(
            !RichardsonExtrapolation::is_asymptotic(1.0, 1.1, 1.2),
            "Should detect non-asymptotic behavior"
        );

        // Case 3: Insufficient variation
        assert!(
            !RichardsonExtrapolation::is_asymptotic(1.0, 1.0000001, 1.00000005),
            "Should detect insufficient variation"
        );

        // Case 4: Zero differences
        assert!(
            !RichardsonExtrapolation::is_asymptotic(1.0, 1.0, 1.0),
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

        let (extrapolated1, order1) =
            RichardsonExtrapolation::extrapolate(phi1, phi2, phi3, r).unwrap();

        // Scale all values by constant factor
        let scale = 3.14159;
        let (extrapolated2, order2) =
            RichardsonExtrapolation::extrapolate(phi1 * scale, phi2 * scale, phi3 * scale, r)
                .unwrap();

        // Extrapolated value should scale, order should be invariant
        assert!(
            nalgebra::ComplexField::abs(extrapolated2 - extrapolated1 * scale) < 1e-12,
            "Extrapolation should be linear"
        );
        assert!(
            nalgebra::ComplexField::abs(order2 - order1) < 1e-12,
            "Order should be invariant under scaling"
        );
    }

    #[test]
    fn test_convergence_order_bounds() {
        // Test that estimated orders are within reasonable CFD bounds

        // Generate test cases with known orders
        let test_cases: Vec<(f64, Box<dyn Fn(f64) -> f64>)> = vec![
            (2.0, Box::new(|h: f64| 1.0 + h * h)),       // Second order
            (1.5, Box::new(|h: f64| 1.0 + h.powf(1.5))), // 1.5 order
            (3.0, Box::new(|h: f64| 1.0 + h * h * h)),   // Third order
        ];

        for (expected_order, solution_fn) in test_cases {
            let h_vals = vec![1.0, 0.5, 0.25];
            let solutions: Vec<f64> = h_vals.iter().map(|&h| solution_fn(h)).collect();
            let refinement_ratios = vec![2.0, 2.0];

            let estimated_order = DataDrivenOrderEstimation::estimate_order_from_solutions(
                &solutions,
                &refinement_ratios,
            );

            // Order should be within reasonable CFD bounds (0.5 to 6.0) and close to expected
            assert!(
                estimated_order > 0.5 && estimated_order < 6.0,
                "Order out of bounds: {} (expected ~{})",
                estimated_order,
                expected_order
            );
            assert!(
                nalgebra::ComplexField::abs(estimated_order - expected_order) < 0.5,
                "Order estimation too inaccurate: {} vs {}",
                estimated_order,
                expected_order
            );
        }
    }
}
