//! Convergence criteria and status definitions
//!
//! Implements convergence assessment following CFD best practices.

use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Convergence status for iterative solvers
#[derive(Debug, Clone, PartialEq)]
pub enum ConvergenceStatus<T: RealField + Copy> {
    /// Converged to specified tolerance
    Converged {
        /// Final error or residual
        final_error: T,
        /// Number of iterations taken
        iterations: usize,
        /// Convergence criterion met
        criterion: ConvergenceCriterion,
    },
    /// Not yet converged
    NotConverged {
        /// Current error or residual
        current_error: T,
        /// Current iteration count
        iterations: usize,
    },
    /// Diverging solution
    Diverging {
        /// Error growth rate
        growth_rate: T,
        /// Iteration where divergence detected
        iterations: usize,
    },
    /// Stalled convergence
    Stalled {
        /// Error at stall point
        stall_error: T,
        /// Iterations since stall
        stall_iterations: usize,
    },
}

impl<T: RealField + Copy> ConvergenceStatus<T> {
    /// Check if converged
    pub fn is_converged(&self) -> bool {
        matches!(self, Self::Converged { .. })
    }

    /// Get current error
    pub fn error(&self) -> T {
        match self {
            Self::Converged { final_error, .. } => *final_error,
            Self::NotConverged { current_error, .. } => *current_error,
            Self::Diverging { .. } => <T as SafeFromF64>::from_f64_or_zero(f64::INFINITY),
            Self::Stalled { stall_error, .. } => *stall_error,
        }
    }
}

/// Convergence criterion type
#[derive(Debug, Clone, PartialEq)]
pub enum ConvergenceCriterion {
    /// Absolute tolerance on residual or error
    Absolute,
    /// Relative tolerance on change
    Relative,
    /// Combined absolute and relative
    Combined,
    /// Machine precision reached
    MachinePrecision,
    /// Maximum iterations reached
    MaxIterations,
}

/// Grid Convergence Index (GCI) calculator
///
/// Implements Roache's GCI method for uncertainty quantification
#[derive(Debug, Clone)]
pub struct GridConvergenceIndex<T: RealField + Copy> {
    /// Safety factor (1.25 for 3+ grids, 3.0 for 2 grids)
    pub safety_factor: T,
    /// Observed order of accuracy
    pub order: T,
    /// Grid refinement ratio
    pub refinement_ratio: T,
}

impl<T: RealField + Copy + FromPrimitive> GridConvergenceIndex<T> {
    /// Create GCI calculator with recommended safety factor
    pub fn new(num_grids: usize, order: T, refinement_ratio: T) -> Self {
        let safety_factor = if num_grids >= 3 {
            T::from_f64_or_one(1.25) // Recommended for systematic studies
        } else {
            T::from_f64_or_one(3.0) // Conservative for limited grids
        };

        Self {
            safety_factor,
            order,
            refinement_ratio,
        }
    }

    /// Compute GCI for fine grid solution
    pub fn compute_fine(&self, f_fine: T, f_coarse: T) -> T {
        let epsilon = (f_coarse - f_fine).abs() / f_fine.abs();
        let r_p = self.refinement_ratio.powf(self.order);

        self.safety_factor * epsilon / (r_p - T::one())
    }

    /// Compute GCI for coarse grid solution
    pub fn compute_coarse(&self, f_fine: T, f_coarse: T) -> T {
        let r_p = self.refinement_ratio.powf(self.order);
        r_p * self.compute_fine(f_fine, f_coarse)
    }

    /// Check if solutions are in asymptotic range
    ///
    /// Returns true if `GCI_coarse` / (r^p * `GCI_fine`) ≈ 1
    pub fn is_asymptotic(&self, gci_fine: T, gci_coarse: T) -> bool {
        let r_p = self.refinement_ratio.powf(self.order);
        let ratio = gci_coarse / (r_p * gci_fine);

        // Should be within 3% of unity for asymptotic range
        (ratio - T::one()).abs() < <T as SafeFromF64>::from_f64_or_zero(0.03)
    }

    /// Compute uncertainty band for solution
    pub fn uncertainty_band(&self, f_fine: T, gci: T) -> (T, T) {
        (f_fine - gci, f_fine + gci)
    }
}

/// Convergence monitor for tracking solver progress
#[derive(Debug, Clone)]
pub struct ConvergenceMonitor<T: RealField + Copy> {
    /// History of errors/residuals
    pub history: Vec<T>,
    /// Absolute tolerance
    pub abs_tolerance: T,
    /// Relative tolerance
    pub rel_tolerance: T,
    /// Maximum allowed iterations
    pub max_iterations: usize,
    /// Stall detection window
    pub stall_window: usize,
}

impl<T: RealField + Copy + FromPrimitive + std::iter::Sum> ConvergenceMonitor<T> {
    /// Create a new convergence monitor
    pub fn new(abs_tol: T, rel_tol: T, max_iter: usize) -> Self {
        Self {
            history: Vec::with_capacity(max_iter),
            abs_tolerance: abs_tol,
            rel_tolerance: rel_tol,
            max_iterations: max_iter,
            stall_window: 10,
        }
    }

    /// Update monitor with new error value
    pub fn update(&mut self, error: T) {
        self.history.push(error);
    }

    /// Check convergence status
    pub fn check_status(&self) -> ConvergenceStatus<T> {
        if self.history.is_empty() {
            return ConvergenceStatus::NotConverged {
                current_error: <T as SafeFromF64>::from_f64_or_zero(f64::INFINITY),
                iterations: 0,
            };
        }

        let Some(&current_error) = self.history.last() else {
            return ConvergenceStatus::NotConverged {
                current_error: <T as SafeFromF64>::from_f64_or_zero(f64::INFINITY),
                iterations: 0,
            };
        };
        let iterations = self.history.len();

        // Check for stalled convergence FIRST using coefficient of variation
        // This must be checked before relative convergence to distinguish stalling from convergence
        if iterations >= self.stall_window {
            let window_start = iterations - self.stall_window;
            let window_errors = &self.history[window_start..];
            let mean_error = window_errors.iter().copied().sum::<T>()
                / T::from_usize(self.stall_window).unwrap();

            // Avoid division by zero
            if mean_error > T::zero() {
                let variance = window_errors
                    .iter()
                    .map(|e| (*e - mean_error).powi(2))
                    .sum::<T>()
                    / T::from_usize(self.stall_window).unwrap();

                let std_dev = variance.sqrt();
                // Use coefficient of variation (CV) for scale-invariant stall detection
                let cv = std_dev / mean_error;

                // Stalled if CV is very small (< 1% of relative tolerance)
                // AND error is still above absolute tolerance (otherwise it's converged)
                if cv < self.rel_tolerance * <T as SafeFromF64>::from_f64_or_zero(0.01)
                    && current_error > self.abs_tolerance
                {
                    return ConvergenceStatus::Stalled {
                        stall_error: current_error,
                        stall_iterations: self.stall_window,
                    };
                }
            }
        }

        // Check relative convergence
        if iterations > 1 {
            let prev_error = self.history[iterations - 2];

            // Avoid division by zero
            if prev_error > T::zero() {
                let rel_change = (current_error - prev_error).abs() / prev_error;

                if rel_change < self.rel_tolerance {
                    return ConvergenceStatus::Converged {
                        final_error: current_error,
                        iterations,
                        criterion: ConvergenceCriterion::Relative,
                    };
                }

                // Check for divergence (error growing by more than 10%)
                if current_error > prev_error * <T as SafeFromF64>::from_f64_or_one(1.1) {
                    let growth_rate = current_error / prev_error;
                    return ConvergenceStatus::Diverging {
                        growth_rate,
                        iterations,
                    };
                }
            }
        }

        // Check absolute convergence only if error is small
        if current_error < self.abs_tolerance {
            return ConvergenceStatus::Converged {
                final_error: current_error,
                iterations,
                criterion: ConvergenceCriterion::Absolute,
            };
        }

        // Check max iterations
        if iterations >= self.max_iterations {
            return ConvergenceStatus::Converged {
                final_error: current_error,
                iterations,
                criterion: ConvergenceCriterion::MaxIterations,
            };
        }

        ConvergenceStatus::NotConverged {
            current_error,
            iterations,
        }
    }

    /// Get convergence rate over recent iterations
    pub fn convergence_rate(&self, window: usize) -> Option<T> {
        if self.history.len() < window + 1 {
            return None;
        }

        let start = self.history.len() - window - 1;
        let previous_error = self.history[start];
        let current_error = *self.history.last().unwrap();

        if previous_error <= T::zero() || current_error <= T::zero() {
            return None;
        }

        Some((current_error / previous_error).powf(T::one() / T::from_usize(window).unwrap()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gci_calculation() {
        let gci = GridConvergenceIndex::<f64>::new(3, 2.0, 2.0);

        let f_fine = 1.234;
        let f_coarse = 1.256;

        let gci_fine = gci.compute_fine(f_fine, f_coarse);
        assert!(gci_fine > 0.0);
        assert!(gci_fine < 0.01); // Should be small for good convergence
    }

    #[test]
    fn test_convergence_monitor() {
        let mut monitor = ConvergenceMonitor::<f64>::new(1e-5, 1e-3, 100);

        // Simulate convergence - final error should be below absolute tolerance
        let errors = vec![1.0, 0.1, 0.01, 0.001, 0.0001, 0.000001];
        for e in errors {
            monitor.update(e);
        }

        let status = monitor.check_status();
        assert!(status.is_converged());
    }

    #[test]
    fn test_divergence_detection() {
        let mut monitor = ConvergenceMonitor::<f64>::new(1e-6, 1e-3, 100);

        // Simulate divergence
        monitor.update(0.1);
        monitor.update(0.5);
        monitor.update(2.0);

        let status = monitor.check_status();
        assert!(matches!(status, ConvergenceStatus::Diverging { .. }));
    }
}
