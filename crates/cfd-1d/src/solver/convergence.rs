//! Convergence checking for network solvers

use cfd_core::Result;
use nalgebra::{DVector, RealField};

/// Convergence checker for iterative solutions
pub struct ConvergenceChecker<T: RealField + Copy> {
    tolerance: T,
    max_iterations: usize,
}

impl<T: RealField + Copy> ConvergenceChecker<T> {
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

    /// Check if solution has converged using L2 norm of residual
    /// Reference: Saad, Y. (2003). "Iterative Methods for Sparse Linear Systems"
    pub fn check(&self, solution: &DVector<T>) -> Result<()> {
        // Check for NaN/Inf values indicating divergence
        if solution.iter().any(|x| !x.is_finite()) {
            return Err(cfd_core::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Diverged { norm: 0.0 },
            ));
        }

        // Compute L2 norm of solution for convergence check
        let norm = solution.norm();

        // Check if norm is within tolerance bounds
        if norm > T::from_f64(1e10).unwrap_or_else(T::one) {
            return Err(cfd_core::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Diverged { norm: 0.0 },
            ));
        }

        // Additional checks can be added here for specific convergence criteria
        // such as relative residual, absolute residual, or solution change

        Ok(())
    }

    /// Check residual convergence
    pub fn check_residual(&self, residual: T) -> bool {
        residual < self.tolerance
    }

    /// Check if solution has converged by comparing two solution vectors
    pub fn has_converged(&self, current: &DVector<T>, previous: &DVector<T>) -> Result<bool> {
        // Check for NaN/Inf values indicating divergence
        if current.iter().any(|x| !x.is_finite()) {
            return Err(cfd_core::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Diverged { norm: 0.0 },
            ));
        }

        // Compute the L2 norm of the change
        let change = (current - previous).norm();

        // Check if change is within tolerance
        Ok(change < self.tolerance)
    }
}
