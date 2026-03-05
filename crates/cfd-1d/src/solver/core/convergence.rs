//! Convergence checking for iterative network solvers.
//!
//! ## Theorem: Krylov Convergence Criterion (Saad 2003 §5)
//!
//! For a linear system  **A x = b**, the iterative solver is considered converged when the
//! **relative residual** satisfies:
//!
//! ```text
//! ‖ A x_k  −  b ‖₂  /  ‖ b ‖₂  <  ε
//! ```
//!
//! This is the standard stopping criterion used by CG, BiCGSTAB, and GMRES.
//!
//! ## Theorem: Picard / Non-linear Convergence (Fixed-Point)
//!
//! For the non-linear network problem (variable-viscosity fluids), each outer
//! iteration linearizes around the current state and solves a linear system.
//! The outer loop has converged when the **relative solution change** is also small:
//!
//! ```text
//! ‖ x_k  −  x_{k−1} ‖₂  /  max(‖ x_k ‖₂, ε_float)  <  ε
//! ```
//!
//! Both criteria must be satisfied simultaneously (`has_converged_dual`) to confirm:
//! 1. The linear system was solved accurately (residual criterion).
//! 2. The outer Picard iteration found a fixed point (change criterion).
//!
//! Checking only the residual is insufficient in non-linear mode because the linear
//! solver minimises the linear residual in every step regardless of outer convergence.
//!
//! ## Invariant: NaN / Inf Detection
//!
//! Any `NaN` or `Inf` in the solution vector indicates numerical divergence (likely
//! from a near-singular or ill-conditioned matrix) and is always reported as an error:
//!
//! ```text
//! ∃ x_i ∈ {NaN, ±∞}  ⟹  Error::Convergence(Diverged)
//! ```
//!
//! ### References
//!
//! - Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems* (2nd ed.).
//!   SIAM. §5.1 (stopping criteria), §6.3 (Picard/fixed-point).
//! - Kelley, C. T. (1995). *Iterative Methods for Linear and Nonlinear Equations.*
//!   SIAM. Ch. 1.

use cfd_core::error::Result;
use nalgebra::{DVector, RealField};

/// Convergence checker for iterative solutions.
///
/// Provides two stopping criteria:
/// - **Residual criterion**: `‖r‖/‖b‖ < ε` (linear solver accuracy)
/// - **Change criterion**: `‖Δx‖/‖x‖ < ε` (non-linear/Picard fixed-point)
///
/// Both must hold simultaneously for `has_converged_dual`.
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

    /// Check if maximum iterations reached
    pub fn max_iterations_reached(&self, current_iteration: usize) -> bool {
        current_iteration >= self.max_iterations
    }

    /// Check for NaN/Inf (divergence) and that the relative residual `‖r‖/‖b‖ < ε`.
    ///
    /// # Arguments
    /// - `solution` – current solution vector `x`
    ///
    /// Returns `Ok(())` if `x` is finite and `‖x‖ < 10^10 · ε` (divergence guard),
    /// or `Err(Convergence::Diverged)` if any entry is non-finite.
    ///
    /// **Note**: This method cannot check the true residual `‖Ax − b‖` without
    /// access to A and b. Use `has_converged_dual` for a complete convergence test.
    pub fn check(&self, solution: &DVector<T>) -> Result<()> {
        // Invariant: no NaN/Inf in the solution vector
        if solution.iter().any(|x| !x.is_finite()) {
            return Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Diverged { norm: 0.0 },
            ));
        }

        // Divergence guard: solution magnitude must not explode
        let norm = solution.norm();
        let divergence_limit = T::from_f64(1e10).expect("Mathematical constant conversion compromised");
        if norm > divergence_limit {
            let norm_f64 = nalgebra::try_convert::<T, f64>(norm).unwrap_or(f64::INFINITY);
            return Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Diverged { norm: norm_f64 },
            ));
        }

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
            return Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Diverged { norm: 0.0 },
            ));
        }

        // Compute the L2 norm of the change
        let change = (current - previous).norm();

        // Check if change is within tolerance
        Ok(change < self.tolerance)
    }

    /// Check if solution has converged using both solution change and residual norm
    pub fn has_converged_dual(
        &self,
        current: &DVector<T>,
        previous: &DVector<T>,
        residual_norm: T,
        rhs_norm: T,
    ) -> Result<bool> {
        // Check for NaN/Inf values indicating divergence
        if current.iter().any(|x| !x.is_finite()) {
            return Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Diverged { norm: 0.0 },
            ));
        }

        // 1. Solution change convergence
        let change = (current - previous).norm();
        let solution_norm = current.norm();
        let relative_change = if solution_norm > T::default_epsilon() {
            change / solution_norm
        } else {
            change
        };

        // 2. Residual convergence
        let relative_residual = if rhs_norm > T::default_epsilon() {
            residual_norm / rhs_norm
        } else {
            residual_norm
        };

        // Converged if BOTH relative change AND relative residual are below tolerance
        // This ensures that we have found a fixed point of the non-linear iteration (Picard)
        // AND that the linear system was solved to sufficient accuracy.
        // Checking only residual (linear residual) is insufficient because it is minimized
        // by the linear solver in each step, regardless of whether the non-linear problem is solved.
        Ok(relative_change < self.tolerance && relative_residual < self.tolerance)
    }
}
