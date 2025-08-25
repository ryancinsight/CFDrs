//! Core traits for linear solvers and preconditioners

use cfd_core::error::Result;
use cfd_core::solver::LinearSolverConfig;
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;

/// Trait for linear solvers
pub trait LinearSolver<T: RealField + Copy>: Send + Sync {
    /// Solve Ax = b
    fn solve(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>>;

    /// Get solver configuration
    fn config(&self) -> &LinearSolverConfig<T>;

    /// Check if residual satisfies convergence criteria
    fn is_converged(&self, residual_norm: T) -> bool {
        residual_norm < self.config().tolerance
    }
}

/// Simplified, efficient preconditioner trait
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
