//! Core traits for linear solvers and preconditioners

use super::config::IterativeSolverConfig;
use cfd_core::error::Result;
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;

/// Configuration trait for solvers
pub trait Configurable<T: RealField + Copy> {
    /// Configuration type
    type Config;

    /// Get solver configuration
    fn config(&self) -> &Self::Config;
}

/// Object-safe trait for linear solvers used with trait objects
pub trait LinearSolver<T: RealField + Copy>: Send + Sync {
    /// Solve Ax = b, returning the solution vector
    fn solve_system(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>>;
}

/// Trait for iterative linear solvers
/// Operates on pre-allocated vectors to avoid repeated allocations
pub trait IterativeLinearSolver<T: RealField + Copy>:
    Send + Sync + Configurable<T, Config = IterativeSolverConfig<T>>
{
    /// Solve Ax = b in-place
    ///
    /// # Arguments
    /// * `a` - System matrix
    /// * `b` - Right-hand side vector
    /// * `x` - Solution vector (also serves as initial guess)
    /// * `preconditioner` - Optional preconditioner
    /// 
    /// # Errors
    /// Returns error if solver fails to converge or encounters numerical issues
    fn solve<P: Preconditioner<T>>(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x: &mut DVector<T>,
        preconditioner: Option<&P>,
    ) -> Result<()>;
}

/// Preconditioner trait for linear system solvers
///
/// This API enforces explicit memory management and avoids hidden allocations
/// by requiring the user to provide the output vector.
pub trait Preconditioner<T: RealField + Copy>: Send + Sync {
    /// Apply the preconditioner: M(z) = r
    ///
    /// Solves the preconditioning system and stores the result in `z`.
    /// This approach makes memory management explicit and avoids hidden allocations.
    ///
    /// # Errors
    /// Returns error if preconditioning fails or matrices are incompatible
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()>;
}

/// Convergence monitoring and theoretical bounds for linear solvers
#[derive(Debug, Clone)]
pub struct ConvergenceMonitor<T: RealField + Copy> {
    /// Initial residual norm
    pub initial_residual: T,
    /// Current iteration
    pub iteration: usize,
    /// Residual norms at each iteration
    pub residual_history: Vec<T>,
    /// Theoretical convergence bound (if available)
    pub theoretical_bound: Option<T>,
    /// Estimated condition number (if available)
    pub condition_number_estimate: Option<f64>,
}

impl<T: RealField + Copy> ConvergenceMonitor<T> {
    /// Create new convergence monitor
    pub fn new(initial_residual: T) -> Self {
        Self {
            initial_residual,
            iteration: 0,
            residual_history: vec![initial_residual],
            theoretical_bound: None,
            condition_number_estimate: None,
        }
    }

    /// Record residual at current iteration
    pub fn record_residual(&mut self, residual: T) {
        self.iteration += 1;
        self.residual_history.push(residual);
    }

    /// Set theoretical convergence bound
    pub fn set_theoretical_bound(&mut self, bound: T) {
        self.theoretical_bound = Some(bound);
    }

    /// Set condition number estimate
    pub fn set_condition_number_estimate(&mut self, kappa: f64) {
        self.condition_number_estimate = Some(kappa);
    }

    /// Get convergence factor (reduction per iteration)
    pub fn convergence_factor(&self) -> Option<T> {
        if self.residual_history.len() >= 2 {
            let r0 = self.initial_residual;
            let r_final = *self.residual_history.last()?;
            // For simplicity, compute geometric mean of reduction factors
            let ratio = r_final / r0;
            Some(ratio.powf(T::from_f64(1.0 / self.iteration as f64)?))
        } else {
            None
        }
    }

    /// Get theoretical CG convergence bound: O(√κ)
    pub fn cg_theoretical_bound(&self, kappa: f64) -> T {
        let sqrt_kappa = kappa.sqrt();
        let bound = 2.0 * ((sqrt_kappa - 1.0) / (sqrt_kappa + 1.0));
        T::from_f64(bound).unwrap_or_else(|| T::one())
    }

    /// Get theoretical GMRES convergence bound (approximate)
    pub fn gmres_theoretical_bound(&self, _kappa: f64) -> Option<T> {
        // GMRES bound is more complex, depends on field of values
        // For now, return a conservative estimate
        Some(T::from_f64(0.8).unwrap_or_else(|| T::one())) // Conservative bound for non-symmetric systems
    }

    /// Check if convergence is within theoretical expectations
    pub fn validate_convergence(&self) -> cfd_core::error::Result<()> {
        if let (Some(factor), Some(theoretical)) = (self.convergence_factor(), self.theoretical_bound) {
            if factor > theoretical * T::from_f64(1.5).unwrap_or_else(|| T::one() + T::one()) {
                return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                    "Convergence factor exceeds theoretical bound"
                )));
            }
        }
        Ok(())
    }
}