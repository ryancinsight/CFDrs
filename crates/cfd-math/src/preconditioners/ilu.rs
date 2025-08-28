//! Incomplete LU factorization preconditioner

use crate::linear_solver::{LinearSolverError, Preconditioner, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;

/// Incomplete LU factorization preconditioner (ILU(0))
///
/// Computes L and U factors such that L*U ≈ A with the same sparsity pattern as A.
pub struct IncompleteLU<T: RealField + Copy> {
    /// Combined LU factors (L has unit diagonal)
    lu_factor: CsrMatrix<T>,
    /// Workspace for forward/backward substitution
    workspace: Vec<T>,
}

impl<T: RealField + Copy> IncompleteLU<T> {
    /// Construct ILU(0) preconditioner
    pub fn new(a: &CsrMatrix<T>) -> Result<Self> {
        if a.nrows() != a.ncols() {
            return Err(LinearSolverError::NonSquareMatrix {
                rows: a.nrows(),
                cols: a.ncols(),
            });
        }

        let n = a.nrows();
        let lu_factor = Self::factorize(a)?;

        Ok(Self {
            lu_factor,
            workspace: vec![T::zero(); n],
        })
    }

    /// Perform ILU(0) factorization
    fn factorize(a: &CsrMatrix<T>) -> Result<CsrMatrix<T>> {
        // Create a mutable copy of the matrix
        let mut lu_vals = a.values().to_vec();
        let lu_indices = a.col_indices().to_vec();
        let lu_offsets = a.row_offsets().to_vec();

        let n = a.nrows();

        // ILU(0) factorization
        for k in 0..n - 1 {
            let diag_idx = Self::find_diagonal_index(&lu_offsets, &lu_indices, k)?;
            let a_kk = lu_vals[diag_idx];

            if a_kk.abs() <= T::default_epsilon() {
                return Err(LinearSolverError::SingularMatrix { pivot_row: k });
            }

            // Process rows below k
            for i in k + 1..n {
                let row_start = lu_offsets[i];
                let row_end = lu_offsets[i + 1];

                // Find a_ik
                let mut a_ik_idx = None;
                for idx in row_start..row_end {
                    if lu_indices[idx] == k {
                        a_ik_idx = Some(idx);
                        break;
                    }
                }

                if let Some(idx) = a_ik_idx {
                    // Compute multiplier
                    let mult = lu_vals[idx] / a_kk;
                    lu_vals[idx] = mult;

                    // Update row i
                    for j_idx in row_start..row_end {
                        let j = lu_indices[j_idx];
                        if j > k {
                            // Find a_kj in row k
                            let k_row_start = lu_offsets[k];
                            let k_row_end = lu_offsets[k + 1];

                            for k_idx in k_row_start..k_row_end {
                                if lu_indices[k_idx] == j {
                                    lu_vals[j_idx] -= mult * lu_vals[k_idx];
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        CsrMatrix::try_from_csr_data(n, n, lu_offsets, lu_indices, lu_vals)
            .map_err(|_| LinearSolverError::InvalidMatrix)
    }

    /// Find index of diagonal element in row
    fn find_diagonal_index(offsets: &[usize], indices: &[usize], row: usize) -> Result<usize> {
        let row_start = offsets[row];
        let row_end = offsets[row + 1];

        for idx in row_start..row_end {
            if indices[idx] == row {
                return Ok(idx);
            }
        }

        Err(LinearSolverError::InvalidMatrix)
    }

    /// Forward substitution with L (unit diagonal)
    fn forward_substitution(&self, b: &DVector<T>, y: &mut DVector<T>) {
        let n = self.lu_factor.nrows();

        for i in 0..n {
            let mut sum = b[i];

            let row_start = self.lu_factor.row_offsets()[i];
            let row_end = self.lu_factor.row_offsets()[i + 1];

            for idx in row_start..row_end {
                let j = self.lu_factor.col_indices()[idx];
                if j < i {
                    sum -= self.lu_factor.values()[idx] * y[j];
                }
            }

            y[i] = sum;
        }
    }

    /// Backward substitution with U
    fn backward_substitution(&self, y: &DVector<T>, x: &mut DVector<T>) {
        let n = self.lu_factor.nrows();

        for i in (0..n).rev() {
            let mut sum = y[i];

            let row_start = self.lu_factor.row_offsets()[i];
            let row_end = self.lu_factor.row_offsets()[i + 1];

            let mut diag_val = T::one();

            for idx in row_start..row_end {
                let j = self.lu_factor.col_indices()[idx];
                if j > i {
                    sum -= self.lu_factor.values()[idx] * x[j];
                } else if j == i {
                    diag_val = self.lu_factor.values()[idx];
                }
            }

            x[i] = sum / diag_val;
        }
    }
}

impl<T: RealField + Copy> Preconditioner<T> for IncompleteLU<T> {
    fn apply(&mut self, b: &DVector<T>) -> DVector<T> {
        let n = b.len();
        let mut y = DVector::zeros(n);
        let mut x = DVector::zeros(n);

        self.forward_substitution(b, &mut y);
        self.backward_substitution(&y, &mut x);

        x
    }

    fn name(&self) -> &str {
        "Incomplete LU (ILU(0))"
    }
}