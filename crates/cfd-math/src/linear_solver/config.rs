//! Convergence configuration for CFD iterative linear solves.
//!
//! Athena owns the convergence policy that the Krylov recurrences enforce
//! (Atlas ADR 0033). This type is the CFD-domain-facing form of that policy:
//! the tolerances and iteration budget a CFD solver carries in its own
//! configuration, converted once at the solve boundary by
//! [`krylov::convergence_policy`](super::krylov::convergence_policy).
//!
//! Keeping the domain-side form here rather than reaching for Athena's
//! validated `ConvergencePolicy` directly serves two properties Athena's type
//! deliberately does not: it is constructible without a fallible step, so a
//! CFD configuration struct can hold one as a plain field, and it defaults,
//! so partially specified configurations read as struct literals.

use eunomia::{FloatElement, NumericElement, RealField};

/// Absolute residual tolerance a configuration defaults to.
///
/// Chosen as the coarsest tolerance at which a double-precision CFD residual
/// is still dominated by discretization rather than round-off; every solver in
/// the workspace overrides it with a problem-derived value.
const DEFAULT_TOLERANCE: f64 = 1e-6;

/// Iteration budget a configuration defaults to.
const DEFAULT_MAX_ITERATIONS: usize = 1000;

/// Convergence configuration shared by the CFD Krylov entry points.
///
/// A solve stops when the residual norm falls to
///
/// ```text
/// threshold = max(tolerance, relative_tolerance · ‖b‖₂)
/// ```
///
/// or when `max_iterations` is exhausted, whichever comes first.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IterativeSolverConfig<T> {
    /// Maximum number of iterations before the solve reports exhaustion.
    pub max_iterations: usize,
    /// Absolute residual norm tolerance for convergence.
    pub tolerance: T,
    /// Relative residual tolerance, measured against `‖b‖₂`.
    ///
    /// Zero (the default) reduces the threshold to the purely absolute test.
    /// A non-zero value makes convergence invariant under a uniform rescaling
    /// of `A·x = b`, which the absolute test alone is not: scaling the system
    /// by `α` scales every residual by `α`, so a fixed absolute tolerance
    /// becomes unreachable for `α ≫ 1` and is satisfied by `x = 0` for
    /// `α ≪ 1`.
    pub relative_tolerance: T,
}

impl<T: RealField + Copy> IterativeSolverConfig<T> {
    /// Create a configuration with the given absolute convergence tolerance.
    ///
    /// The relative tolerance starts at zero and the iteration budget at the
    /// same default [`Self::default`] uses.
    #[must_use]
    pub fn new(tolerance: T) -> Self {
        Self {
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance,
            relative_tolerance: <T as NumericElement>::ZERO,
        }
    }

    /// Convergence threshold for a right-hand side of norm `rhs_norm`.
    ///
    /// A solve has converged when `‖b − A·x‖₂ ≤ threshold(‖b‖₂)`. This is the
    /// same rule Athena applies internally; it is exposed so a caller that
    /// validates a returned iterate against its own residual uses one
    /// definition rather than a restated one.
    #[must_use]
    pub fn threshold(&self, rhs_norm: T) -> T {
        let relative = self.relative_tolerance * rhs_norm;
        if relative > self.tolerance {
            relative
        } else {
            self.tolerance
        }
    }

    /// Override the iteration budget.
    #[must_use]
    pub const fn with_max_iterations(mut self, max_iterations: usize) -> Self {
        self.max_iterations = max_iterations;
        self
    }

    /// Set the relative residual tolerance (see
    /// [`Self::relative_tolerance`]).
    #[must_use]
    pub const fn with_relative_tolerance(mut self, relative_tolerance: T) -> Self {
        self.relative_tolerance = relative_tolerance;
        self
    }
}

impl<T: RealField + Copy + FloatElement> Default for IterativeSolverConfig<T> {
    fn default() -> Self {
        Self {
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance: <T as FloatElement>::from_f64(DEFAULT_TOLERANCE),
            relative_tolerance: <T as NumericElement>::ZERO,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::IterativeSolverConfig;

    #[test]
    fn the_threshold_is_the_larger_of_the_absolute_and_relative_tests() {
        let config = IterativeSolverConfig::new(1e-8_f64).with_relative_tolerance(1e-10);
        // ‖b‖₂ = 1: the relative test contributes 1e-10, below the absolute.
        assert_eq!(config.threshold(1.0), 1e-8);
        // ‖b‖₂ = 1e6: the relative test contributes 1e-4 and now dominates.
        assert_eq!(config.threshold(1e6), 1e-4);
    }

    #[test]
    fn a_zero_relative_tolerance_reduces_the_threshold_to_the_absolute_test() {
        let config = IterativeSolverConfig::new(1e-8_f64);
        assert_eq!(config.threshold(1e12), 1e-8);
    }

    #[test]
    fn the_builders_leave_the_untouched_fields_at_their_defaults() {
        let config = IterativeSolverConfig::<f64>::default().with_max_iterations(7);
        assert_eq!(config.max_iterations, 7);
        assert_eq!(config.tolerance, 1e-6);
        assert_eq!(config.relative_tolerance, 0.0);
    }
}
