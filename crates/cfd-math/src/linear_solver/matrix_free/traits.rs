//! Traits for matrix-free linear solvers.
//!
//! This module defines the interface for solvers that work with
//! `LinearOperator` instances instead of explicit sparse matrices.

use super::operator::LinearOperator;
use crate::error::Result;
use nalgebra::RealField;

/// Trait for matrix-free linear solvers.
///
/// This trait defines the interface for iterative solvers that work
/// with linear operators instead of explicit sparse matrices.
pub trait MatrixFreeSolver<T: RealField + Copy> {
    /// Solve Ax = b using a matrix-free operator.
    ///
    /// This is a simplified interface for demonstration. Full preconditioner
    /// support will be added in future iterations.
    ///
    /// # Arguments
    ///
    /// * `operator` - Linear operator representing A
    /// * `b` - Right-hand side vector
    /// * `x` - Initial guess (input) / solution (output)
    ///
    /// # Errors
    ///
    /// Returns error if convergence fails or operator application fails.
    fn solve(&self, operator: &dyn LinearOperator<T>, b: &[T], x: &mut [T]) -> Result<()>;
}
