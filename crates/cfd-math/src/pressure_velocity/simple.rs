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

use eunomia::{FloatElement, RealField};
use serde::{Deserialize, Serialize};

fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

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
    /// Number of pressure correction loops per outer iteration (1 = SIMPLE, >1 = PISO)
    pub n_correctors: usize,
}

impl<T: RealField + FloatElement + Copy> Default for SIMPLEConfig<T> {
    /// Standard SIMPLE relaxation factors following Patankar (1980) §6.7
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            tolerance: from_f64(1e-6),
            alpha_u: from_f64(0.7),
            alpha_p: from_f64(0.3),
            alpha_mu: T::ONE,
            viscosity_update_interval: 1,
            n_correctors: 1, // Default to traditional SIMPLE
        }
    }
}

impl<T: RealField + Copy> SIMPLEConfig<T> {
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
        n_correctors: usize,
    ) -> Self {
        assert!(
            alpha_u > T::ZERO && alpha_u <= T::ONE,
            "alpha_u must be in (0, 1]"
        );
        assert!(
            alpha_p > T::ZERO && alpha_p <= T::ONE,
            "alpha_p must be in (0, 1]"
        );
        assert!(
            alpha_mu > T::ZERO && alpha_mu <= T::ONE,
            "alpha_mu must be in (0, 1]"
        );
        Self {
            max_iterations,
            tolerance,
            alpha_u,
            alpha_p,
            alpha_mu,
            viscosity_update_interval,
            n_correctors,
        }
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
        Self {
            iterations,
            residual,
            converged,
        }
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

    #[test]
    fn explicit_config_preserves_relaxation_values() {
        let cfg = SIMPLEConfig::new(25, 1e-9, 0.6, 0.4, 0.8, 2, 3);

        assert_eq!(cfg.max_iterations, 25);
        assert_eq!(cfg.tolerance, 1e-9);
        assert_eq!(cfg.alpha_u, 0.6);
        assert_eq!(cfg.alpha_p, 0.4);
        assert_eq!(cfg.alpha_mu, 0.8);
        assert_eq!(cfg.viscosity_update_interval, 2);
        assert_eq!(cfg.n_correctors, 3);
    }
}
