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

/// Trait for iterative linear solvers
/// Operates on pre-allocated vectors to avoid repeated allocations
pub trait IterativeLinearSolver<T: RealField + Copy>: Send + Sync + Configurable<T, Config = IterativeSolverConfig<T>> {
    /// Solve Ax = b in-place
    /// 
    /// # Arguments
    /// * `a` - System matrix
    /// * `b` - Right-hand side vector
    /// * `x` - Solution vector (also serves as initial guess)
    /// * `preconditioner` - Optional preconditioner
    fn solve<P: Preconditioner<T>>(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x: &mut DVector<T>,
        preconditioner: Option<&P>,
    ) -> Result<()>;
}

/// Helper function for convergence checking
pub fn has_converged<T: RealField>(residual_norm: T, tolerance: T) -> bool {
    residual_norm < tolerance
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
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()>;
}
