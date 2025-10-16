//! Triangular system solvers for ILU preconditioner

use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;

/// Forward substitution with L (unit diagonal)
///
/// Solves L*y = b where L is lower triangular with unit diagonal
pub fn forward_substitution<T: RealField + Copy>(
    lu_factor: &CsrMatrix<T>,
    b: &DVector<T>,
    y: &mut DVector<T>,
) {
    let n = lu_factor.nrows();

    for i in 0..n {
        let mut sum = b[i];

        let row_start = lu_factor.row_offsets()[i];
        let row_end = lu_factor.row_offsets()[i + 1];

        for idx in row_start..row_end {
            let j = lu_factor.col_indices()[idx];
            if j < i {
                sum -= lu_factor.values()[idx] * y[j];
            }
        }

        y[i] = sum;
    }
}

/// Backward substitution with U
///
/// Solves U*x = y where U is upper triangular
pub fn backward_substitution<T: RealField + Copy>(
    lu_factor: &CsrMatrix<T>,
    y: &DVector<T>,
    x: &mut DVector<T>,
) {
    let n = lu_factor.nrows();

    for i in (0..n).rev() {
        let mut sum = y[i];

        let row_start = lu_factor.row_offsets()[i];
        let row_end = lu_factor.row_offsets()[i + 1];

        let mut diag_val = T::one();

        for idx in row_start..row_end {
            let j = lu_factor.col_indices()[idx];
            if j > i {
                sum -= lu_factor.values()[idx] * x[j];
            } else if j == i {
                diag_val = lu_factor.values()[idx];
            }
        }

        x[i] = sum / diag_val;
    }
}
