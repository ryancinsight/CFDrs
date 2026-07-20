//! Direct sparse solver for linear systems.
//!
//! Uses [`leto_ops::SparseLuSolver`] — the atlas-native sparse direct solver
//! backed by dense partial-pivoting LU — for systems up to `max_size`. The
//! dense LU path in `leto-ops` serves as both the primary sparse solver and
//! the fallback, eliminating the external `rsparse` dependency.
//!
//! # Theorem — LU Factorisation Uniqueness
//!
//! For a nonsingular matrix $A$, the LU factorisation $PA = LU$ exists and is
//! unique up to diagonal scaling when $A$ is strongly regular. The
//! `leto-ops` implementation uses partial pivoting (largest-magnitude pivot
//! selection) to maintain numerical stability.

use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto_ops::{CsrMatrix, RealScalar as LetoRealScalar, SparseLuSolver as LetoCsrLuSolver};

use super::dense_bridge::solve_leto_csr_with_leto_dense_array;

/// Direct sparse solver configuration.
#[derive(Debug, Clone)]
pub struct DirectSparseSolver {
    /// Maximum system size (number of rows) to attempt with direct solve.
    pub max_size: usize,
    /// Ordering strategy (reserved for API compatibility; the atlas-native solver
    /// uses partial-column pivoting regardless of this value).
    pub ordering: i8,
    /// Pivot tolerance for LU factorization.
    pub pivot_tolerance: f64,
}

impl Default for DirectSparseSolver {
    fn default() -> Self {
        Self {
            max_size: 2048,
            ordering: 1,
            pivot_tolerance: 1e-12,
        }
    }
}

impl DirectSparseSolver {
    /// Returns true if the solver should attempt a direct solve for this size.
    pub fn can_handle_size(&self, size: usize) -> bool {
        size <= self.max_size
    }

    /// Solve Ax = b using the atlas-native sparse LU factorization.
    pub fn solve<T>(&self, matrix: &CsrMatrix<T>, rhs: &Array1<T>) -> Result<Array1<T>>
    where
        T: RealField + Copy + FloatElement + LetoRealScalar,
    {
        let nrows = matrix.nrows();
        let ncols = matrix.ncols();
        let rhs_len = rhs.shape()[0];

        if nrows != ncols {
            return Err(Error::InvalidConfiguration(
                "Direct solver requires a square matrix".to_string(),
            ));
        }
        if rhs_len != nrows {
            return Err(Error::InvalidConfiguration(
                "RHS dimension does not match matrix size".to_string(),
            ));
        }

        let atlas_solver = LetoCsrLuSolver {
            max_size: self.max_size,
            pivot_tolerance: self.pivot_tolerance,
        };

        let rhs_slice: Vec<T> = rhs.iter().copied().collect();

        match atlas_solver.solve(matrix, &rhs_slice) {
            Ok(solution) => {
                let mut x = Array1::from_elem([nrows], <T as NumericElement>::ZERO);
                for (i, value) in solution.into_iter().enumerate() {
                    x[i] = value;
                }
                Self::ensure_finite_solution(&x)?;
                Ok(x)
            }
            Err(e) => {
                let sparse_error = e.to_string();
                self.retry_dense_or_error(matrix, rhs, nrows, sparse_error)
            }
        }
    }

    fn retry_dense_or_error<T>(
        &self,
        matrix: &CsrMatrix<T>,
        rhs: &Array1<T>,
        nrows: usize,
        sparse_error: String,
    ) -> Result<Array1<T>>
    where
        T: RealField + Copy + FloatElement + LetoRealScalar,
    {
        let dense_threshold = 1024usize;
        if nrows <= dense_threshold {
            tracing::warn!(
                size = nrows,
                ordering = self.ordering,
                pivot_tolerance = self.pivot_tolerance,
                "Direct sparse LU failed; retrying with dense LU"
            );
            let dense_solution =
                solve_leto_csr_with_leto_dense_array(matrix, rhs).map_err(|dense_error| {
                    Error::Solver(format!(
                        "Direct sparse LU failed: {sparse_error}; dense LU fallback failed: \
                         {dense_error}"
                    ))
                })?;
            Self::ensure_finite_solution(&dense_solution)?;
            return Ok(dense_solution);
        }

        Err(Error::Solver(format!(
            "Direct sparse LU failed: {sparse_error}"
        )))
    }

    fn ensure_finite_solution<T>(solution: &Array1<T>) -> Result<()>
    where
        T: RealField + Copy + FloatElement,
    {
        for (index, value) in solution.iter().enumerate() {
            if !NumericElement::to_f64(*value).is_finite() {
                return Err(Error::ConversionError(format!(
                    "Direct LU solution entry {index} is not finite"
                )));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::SparseMatrixBuilder;
    use leto::Array1;
    use std::panic::{catch_unwind, AssertUnwindSafe};

    #[test]
    fn direct_solver_solves_small_system() {
        let mut builder = SparseMatrixBuilder::new(2, 2);
        builder.add_entry(0, 0, 3.0_f64).unwrap();
        builder.add_entry(0, 1, 1.0_f64).unwrap();
        builder.add_entry(1, 0, 1.0_f64).unwrap();
        builder.add_entry(1, 1, 2.0_f64).unwrap();

        let mut assembly_rhs = Array1::from_shape_vec([2], vec![9.0_f64, 8.0_f64]).unwrap();
        let matrix = builder.build_with_rhs(&mut assembly_rhs).unwrap();

        let solver = DirectSparseSolver::default();
        let rhs = Array1::from_shape_vec([2], vec![9.0_f64, 8.0_f64]).unwrap();
        let x = solver.solve(&matrix, &rhs).unwrap();

        assert!((x[0] - 2.0_f64).abs() < 1e-8_f64);
        assert!((x[1] - 3.0_f64).abs() < 1e-8_f64);
    }

    #[test]
    fn direct_solver_singular_system_does_not_panic() {
        let builder = SparseMatrixBuilder::new(2, 2);
        let mut assembly_rhs = Array1::from_shape_vec([2], vec![0.0_f64, 0.0_f64]).unwrap();
        let matrix = builder.build_with_rhs(&mut assembly_rhs).unwrap();

        let solver = DirectSparseSolver::default();
        let rhs = Array1::from_shape_vec([2], vec![0.0_f64, 0.0_f64]).unwrap();
        let result = catch_unwind(AssertUnwindSafe(|| solver.solve(&matrix, &rhs)));

        assert!(
            result.is_ok(),
            "direct solver must not panic on singular input"
        );
        assert!(result.unwrap().is_err());
    }

    #[test]
    fn direct_solver_falls_back_to_dense_for_small_system_exceeding_atlas_limit() {
        // With max_size=1, a 2×2 system exceeds the atlas-LU limit but
        // falls through to the dense fallback (n=2 ≤ dense_threshold=1024).
        let solver = DirectSparseSolver {
            max_size: 1,
            ordering: 1,
            pivot_tolerance: 1e-12,
        };
        let mut builder = SparseMatrixBuilder::new(2, 2);
        builder.add_entry(0, 0, 2.0_f64).unwrap();
        builder.add_entry(1, 1, 4.0_f64).unwrap();
        let mut assembly_rhs = Array1::from_shape_vec([2], vec![6.0_f64, 8.0_f64]).unwrap();
        let matrix = builder.build_with_rhs(&mut assembly_rhs).unwrap();
        let rhs = Array1::from_shape_vec([2], vec![6.0_f64, 8.0_f64]).unwrap();
        // Should succeed via dense fallback.
        let x = solver.solve(&matrix, &rhs).expect("dense fallback should handle small n");
        assert!((x[0] - 3.0_f64).abs() < 1e-10);
        assert!((x[1] - 2.0_f64).abs() < 1e-10);
    }

    #[test]
    fn dense_lu_fallback_uses_leto_provider() {
        let mut builder = SparseMatrixBuilder::new(2, 2);
        builder.add_entry(0, 0, 4.0_f64).unwrap();
        builder.add_entry(0, 1, 1.0_f64).unwrap();
        builder.add_entry(1, 0, 2.0_f64).unwrap();
        builder.add_entry(1, 1, 3.0_f64).unwrap();

        let mut assembly_rhs = Array1::from_shape_vec([2], vec![1.0_f64, 7.0_f64]).unwrap();
        let matrix = builder.build_with_rhs(&mut assembly_rhs).unwrap();

        let rhs = Array1::from_shape_vec([2], vec![1.0_f64, 7.0_f64]).unwrap();
        let x = solve_leto_csr_with_leto_dense_array(&matrix, &rhs).unwrap();

        assert!((x[0] + 0.4_f64).abs() < 1e-12_f64);
        assert!((x[1] - 2.6_f64).abs() < 1e-12_f64);
    }
}
