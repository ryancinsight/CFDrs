//! Incomplete Cholesky factorization preconditioner

use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;

// Tolerance for symmetry checking
const SYMMETRY_TOLERANCE: f64 = 1e-10;

/// Incomplete Cholesky factorization preconditioner (IC(0))
///
/// For symmetric positive definite matrices, computes L such that L*L^T ≈ A
/// where L has the same sparsity pattern as the lower triangular part of A.
pub struct IncompleteCholesky<T: RealField + Copy> {
    /// Lower triangular factor stored in CSR format
    l_factor: CsrMatrix<T>,
}

impl<T: RealField + Copy> IncompleteCholesky<T> {
    /// Construct IC(0) preconditioner from symmetric positive definite matrix
    pub fn new(a: &CsrMatrix<T>) -> Result<Self> {
        // Validate matrix is square
        if a.nrows() != a.ncols() {
            return Err(Error::InvalidInput(format!(
                "Matrix must be square, got {}x{}",
                a.nrows(),
                a.ncols()
            )));
        }

        // Check symmetry
        Self::check_symmetry(a)?;

        // Perform IC(0) factorization
        let l_factor = Self::factorize(a)?;

        Ok(Self { l_factor })
    }

    /// Check matrix symmetry
    fn check_symmetry(a: &CsrMatrix<T>) -> Result<()> {
        let tol = T::from_f64(SYMMETRY_TOLERANCE).unwrap_or_else(T::zero);
        let n = a.nrows();
        let mut max_asymmetry = T::zero();
        let mut max_i = 0;
        let mut max_j = 0;

        // Helper to get value from sparse entry
        let get_value = |i: usize, j: usize| -> T {
            match a.get_entry(i, j) {
                Some(entry) => entry.into_value(),
                None => T::zero(),
            }
        };

        for i in 0..n {
            for j in i + 1..n {
                let a_ij = get_value(i, j);
                let a_ji = get_value(j, i);
                let diff = (a_ij - a_ji).abs();

                if diff > max_asymmetry {
                    max_asymmetry = diff;
                    max_i = i;
                    max_j = j;
                }
            }
        }

        if max_asymmetry > tol {
            return Err(Error::InvalidInput(format!(
                "Matrix is not symmetric: max asymmetry {:.2e} at ({}, {})",
                max_asymmetry.to_subset().unwrap_or(0.0),
                max_i,
                max_j
            )));
        }

        Ok(())
    }

    /// Perform IC(0) factorization
    fn factorize(a: &CsrMatrix<T>) -> Result<CsrMatrix<T>> {
        let n = a.nrows();
        let mut l_vals = Vec::new();
        let mut l_indices = Vec::new();
        let mut l_offsets = vec![0];

        // Storage buffer for rows of L
        let mut l_rows: Vec<Vec<(usize, T)>> = vec![Vec::new(); n];

        // Helper to get value from sparse entry
        let get_value = |i: usize, j: usize| -> T {
            match a.get_entry(i, j) {
                Some(entry) => entry.into_value(),
                None => T::zero(),
            }
        };

        // IC(0) factorization algorithm
        for i in 0..n {
            let mut l_ii = get_value(i, i);

            // Subtract contributions from previous columns
            for &(k, l_ik) in &l_rows[i] {
                if k < i {
                    l_ii -= l_ik * l_ik;
                }
            }

            // Check for positive definiteness
            if l_ii <= T::zero() {
                return Err(Error::Numerical(NumericalErrorKind::NotPositiveDefinite));
            }

            l_ii = l_ii.sqrt();
            l_rows[i].push((i, l_ii));

            // Compute entries in column i below diagonal
            let row_start = a.row_offsets()[i];
            let row_end = a.row_offsets()[i + 1];

            for idx in row_start..row_end {
                let j = a.col_indices()[idx];
                if j > i {
                    // Only process lower triangular part
                    let mut l_ji = get_value(j, i);

                    // Subtract contributions from previous columns
                    for &(k, l_ik) in &l_rows[i] {
                        if k < i {
                            // Find l_jk in row j
                            if let Some(&(_, l_jk)) = l_rows[j].iter().find(|(col, _)| *col == k) {
                                l_ji -= l_ik * l_jk;
                            }
                        }
                    }

                    l_ji /= l_ii;
                    l_rows[j].push((i, l_ji));
                }
            }
        }

        // Convert to CSR format
        for row in &l_rows {
            let mut sorted_row = row.clone();
            sorted_row.sort_by_key(|(col, _)| *col);

            for &(col, val) in &sorted_row {
                l_indices.push(col);
                l_vals.push(val);
            }
            l_offsets.push(l_indices.len());
        }

        CsrMatrix::try_from_csr_data(n, n, l_offsets, l_indices, l_vals).map_err(|_| {
            Error::InvalidInput("Failed to create CSR matrix from IC(0) factorization".to_string())
        })
    }

    /// Forward substitution: solve L*y = b
    fn forward_substitution(&self, b: &DVector<T>, y: &mut DVector<T>) {
        let n = self.l_factor.nrows();

        for i in 0..n {
            let mut sum = b[i];

            let row_start = self.l_factor.row_offsets()[i];
            let row_end = self.l_factor.row_offsets()[i + 1];

            for idx in row_start..row_end {
                let j = self.l_factor.col_indices()[idx];
                if j < i {
                    sum -= self.l_factor.values()[idx] * y[j];
                } else if j == i {
                    sum /= self.l_factor.values()[idx];
                    break;
                }
            }

            y[i] = sum;
        }
    }

    /// Backward substitution: solve L^T*x = y
    fn backward_substitution(&self, y: &DVector<T>, x: &mut DVector<T>) {
        let n = self.l_factor.nrows();

        for i in (0..n).rev() {
            let mut sum = y[i];

            // Process contributions from later rows
            for j in i + 1..n {
                let row_start = self.l_factor.row_offsets()[j];
                let row_end = self.l_factor.row_offsets()[j + 1];

                for idx in row_start..row_end {
                    if self.l_factor.col_indices()[idx] == i {
                        sum -= self.l_factor.values()[idx] * x[j];
                        break;
                    }
                }
            }

            // Divide by diagonal element
            let row_start = self.l_factor.row_offsets()[i];
            let row_end = self.l_factor.row_offsets()[i + 1];

            for idx in row_start..row_end {
                if self.l_factor.col_indices()[idx] == i {
                    x[i] = sum / self.l_factor.values()[idx];
                    break;
                }
            }
        }
    }
}

impl<T: RealField + Copy> Preconditioner<T> for IncompleteCholesky<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        let n = r.len();

        // Allocate intermediate result for two-stage solve
        // NOTE: Cannot use pre-allocated workspace with immutable &self API
        let mut y = DVector::zeros(n);

        // Solve L*L^T*z = r via forward and backward substitution
        self.forward_substitution(r, &mut y);
        self.backward_substitution(&y, z);

        Ok(())
    }
}
