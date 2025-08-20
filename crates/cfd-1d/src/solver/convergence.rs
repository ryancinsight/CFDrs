//! Convergence checking for network solvers

use nalgebra::{RealField, DVector};
use cfd_core::Result;

/// Convergence checker for iterative solutions
pub struct ConvergenceChecker<T: RealField> {
    tolerance: T,
    max_iterations: usize,
}

impl<T: RealField> ConvergenceChecker<T> {
    /// Create a new convergence checker
    pub fn new(tolerance: T) -> Self {
        Self {
            tolerance,
            max_iterations: 1000,
        }
    }

    /// Update tolerance
    pub fn update_tolerance(&mut self, tolerance: T) {
        self.tolerance = tolerance;
    }

    /// Check if solution has converged
    pub fn check(&self, solution: &DVector<T>) -> Result<()> {
        // For now, we assume the linear solver has already checked convergence
        // This is a placeholder for more sophisticated convergence checks
        if solution.iter().any(|x| !x.is_finite()) {
            return Err(cfd_core::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Diverged
            ));
        }
        Ok(())
    }

    /// Check residual convergence
    pub fn check_residual(&self, residual: T) -> bool {
        residual < self.tolerance
    }
}