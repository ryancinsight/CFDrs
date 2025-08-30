//! Linear solver implementations for FDM.
//!
//! Provides iterative solvers for the linear systems arising from FDM discretization.

use cfd_core::error::{Error, Result};
use cfd_math::SparseMatrix;
use nalgebra::{DVector, RealField};
use num_traits::{FromPrimitive, Zero};

use super::config::FdmConfig;

/// Shared Gauss-Seidel linear solver implementation
///
/// Solves the linear system Ax = b using Gauss-Seidel iteration with relaxation.
/// Returns an error if convergence is not achieved within `max_iterations`.
pub fn solve_gauss_seidel<T: RealField + Copy + FromPrimitive + Copy>(
    matrix: &SparseMatrix<T>,
    rhs: &DVector<T>,
    config: &FdmConfig<T>,
    solver_name: &str,
) -> Result<DVector<T>> {
    let n = rhs.len();
    let mut solution: DVector<T> = DVector::from_element(n, T::zero());

    for iteration in 0..config.max_iterations() {
        let mut max_residual = T::zero();

        for (row_idx, row) in matrix.row_iter().enumerate() {
            let mut sum = T::zero();
            let mut diagonal = T::one();

            // Sum contributions from other variables and find diagonal
            for (col_idx, value) in row.col_indices().iter().zip(row.values()) {
                if row_idx == *col_idx {
                    diagonal = *value;
                } else {
                    sum += *value * solution[*col_idx];
                }
            }

            // Check for zero diagonal (singular matrix)
            if diagonal.abs() < T::from_f64(1e-14).unwrap_or_else(|| T::zero()) {
                return Err(Error::InvalidConfiguration(format!(
                    "{solver_name}: Singular matrix detected (zero diagonal)"
                )));
            }

            // Update solution
            let current_value = (rhs[row_idx] - sum) / diagonal;
            let residual = (current_value - solution[row_idx]).abs();

            if residual > max_residual {
                max_residual = residual;
            }

            // Apply relaxation
            solution[row_idx] = solution[row_idx]
                + config.relaxation_factor() * (current_value - solution[row_idx]);
        }

        if config.verbose() && iteration % crate::constants::solver::LOG_INTERVAL == 0 {
            println!("{solver_name} iteration {iteration}: residual = {max_residual:?}");
        }

        if max_residual < config.tolerance() {
            if config.verbose() {
                println!("{} converged in {} iterations", solver_name, iteration + 1);
            }
            return Ok(solution);
        }
    }

    // Convergence failure
    Err(Error::InvalidConfiguration(format!(
        "{}: Failed to converge after {} iterations",
        solver_name,
        config.max_iterations()
    )))
}
