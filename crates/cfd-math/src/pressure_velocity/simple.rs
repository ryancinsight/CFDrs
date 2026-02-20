//! SIMPLE algorithm configuration and result types.
//!
//! ## Algorithm Reference
//! Patankar, S.V. & Spalding, D.B. (1972). A calculation procedure for heat, mass and
//! momentum transfer in three-dimensional parabolic flows.
//! *International Journal of Heat and Mass Transfer*, 15(10), 1787–1806.
//!
//! ## Mathematical Invariant
//! Under-relaxation factors α satisfy 0 < α ≤ 1. The product α_u · α_p should be
//! tuned such that the SIMPLE iteration converges monotonically:
//! α_p ≈ 1 - α_u is the classical recommendation (Patankar 1980, §6.7).

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Configuration for the SIMPLE pressure-velocity coupling algorithm.
///
/// Controls convergence criteria and under-relaxation for the iterative
/// SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) procedure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SIMPLEConfig<T: RealField + Copy> {
    /// Maximum number of outer SIMPLE iterations
    pub max_iterations: usize,
    /// Convergence tolerance on velocity residuals (L∞ norm)
    pub tolerance: T,
    /// Under-relaxation factor for u- and v-momentum equations (0 < α_u ≤ 1)
    pub alpha_u: T,
    /// Under-relaxation factor for pressure correction (0 < α_p ≤ 1)
    pub alpha_p: T,
    /// Under-relaxation factor for viscosity (0 < α_μ ≤ 1)
    pub alpha_mu: T,
    /// Non-Newtonian viscosity update every N outer iterations (1 = every iter)
    pub viscosity_update_interval: usize,
}

impl<T: RealField + Copy + FromPrimitive> Default for SIMPLEConfig<T> {
    /// Standard SIMPLE relaxation factors following Patankar (1980) §6.7
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            tolerance: T::from_f64(1e-6).expect("1e-6 fits in T"),
            alpha_u: T::from_f64(0.7).expect("0.7 fits in T"),
            alpha_p: T::from_f64(0.3).expect("0.3 fits in T"),
            alpha_mu: T::one(),
            viscosity_update_interval: 1,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> SIMPLEConfig<T> {
    /// Construct with explicit parameters.
    ///
    /// # Panics
    /// Panics if α_u or α_p are not in (0, 1].
    #[must_use]
    pub fn new(
        max_iterations: usize,
        tolerance: T,
        alpha_u: T,
        alpha_p: T,
        alpha_mu: T,
        viscosity_update_interval: usize,
    ) -> Self {
        assert!(
            alpha_u > T::zero() && alpha_u <= T::one(),
            "alpha_u must be in (0, 1]"
        );
        assert!(
            alpha_p > T::zero() && alpha_p <= T::one(),
            "alpha_p must be in (0, 1]"
        );
        assert!(
            alpha_mu > T::zero() && alpha_mu <= T::one(),
            "alpha_mu must be in (0, 1]"
        );
        Self { max_iterations, tolerance, alpha_u, alpha_p, alpha_mu, viscosity_update_interval }
    }
}

/// Result returned by a completed SIMPLE solve.
#[derive(Debug, Clone)]
pub struct SolveResult<T> {
    /// Number of outer SIMPLE iterations performed
    pub iterations: usize,
    /// Final velocity residual (L∞ norm of momentum imbalance)
    pub residual: T,
    /// Whether the solve converged within `SIMPLEConfig::tolerance`
    pub converged: bool,
}

impl<T> SolveResult<T> {
    /// Construct a solve result.
    #[must_use]
    pub fn new(iterations: usize, residual: T, converged: bool) -> Self {
        Self { iterations, residual, converged }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_relaxation_satisfies_patankar_recommendation() {
        let cfg = SIMPLEConfig::<f64>::default();
        // Patankar: α_u + α_p ≈ 1
        let sum = cfg.alpha_u + cfg.alpha_p;
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn solve_result_reports_correctly() {
        let r = SolveResult::new(42, 1e-8, true);
        assert!(r.converged);
        assert_eq!(r.iterations, 42);
    }
}
