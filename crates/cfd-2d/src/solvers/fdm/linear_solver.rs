//! Linear solver implementations for FDM.
//!
//! Provides iterative solvers for the linear systems arising from FDM discretization.

use cfd_core::error::{Error, Result};
use cfd_math::SparseMatrix;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

use super::config::FdmConfig;

/// Shared Gauss-Seidel linear solver implementation
///
/// Solves the linear system Ax = b using Gauss-Seidel iteration with relaxation.
///
/// True Gauss-Seidel method: Uses most recently updated solution values (not all old values).
/// For each row i, the update is: x_i^(k+1) = (b_i - Σ_{j<i} a_ij*x_j^(k+1) - Σ_{j>i} a_ij*x_j^(k)) / a_ii
///
/// This provides better convergence than Jacobi (which uses all old values).
///
/// Returns an error if convergence is not achieved within `max_iterations`.
pub fn solve_gauss_seidel<T: RealField + Copy + FromPrimitive>(
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

            // Sum contributions from other variables using UPDATED values (Gauss-Seidel)
            // For j < i: use x_j^(k+1) (already updated in this iteration)
            // For j > i: use x_j^(k) (not yet updated in this iteration)
            for (col_idx, value) in row.col_indices().iter().zip(row.values()) {
                if row_idx == *col_idx {
                    diagonal = *value;
                } else {
                    // Use the current solution value (updated for j < row_idx, old for j > row_idx)
                    sum += *value * solution[*col_idx];
                }
            }

            // Check for zero diagonal (singular matrix)
            if diagonal.abs() < T::from_f64(1e-14).unwrap_or_else(|| T::zero()) {
                return Err(Error::InvalidConfiguration(format!(
                    "{solver_name}: Singular matrix detected (zero diagonal)"
                )));
            }

            // Compute new value: x_i = (b_i - sum) / a_ii
            let new_value = (rhs[row_idx] - sum) / diagonal;

            // Calculate residual before applying relaxation
            let residual = (new_value - solution[row_idx]).abs();
            if residual > max_residual {
                max_residual = residual;
            }

            // Apply successive over-relaxation (SOR): x_i = x_i + ω(x_i_new - x_i)
            solution[row_idx] =
                solution[row_idx] + config.relaxation_factor() * (new_value - solution[row_idx]);
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
