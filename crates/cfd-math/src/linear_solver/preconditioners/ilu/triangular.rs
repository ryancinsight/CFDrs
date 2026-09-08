//! Triangular system solvers for ILU preconditioner

use eunomia::{NumericElement, RealField};
use leto::Array1;
use leto_ops::{CsrMatrix, Scalar as LetoScalar};

/// Forward substitution with L (unit diagonal)
///
/// Solves L*y = b where L is lower triangular with unit diagonal
pub fn forward_substitution<T: RealField + Copy + LetoScalar>(
    lu_factor: &CsrMatrix<T>,
    b: &Array1<T>,
    y: &mut Array1<T>,
) {
    let n = lu_factor.nrows();

    for i in 0..n {
        let mut sum = b[i];

        let row = lu_factor.row(i);
        for (&j, &value) in row.col_indices().iter().zip(row.values()) {
            if j < i {
                sum -= value * y[j];
            }
        }

        y[i] = sum;
    }
}

/// Backward substitution with U
///
/// Solves U*x = y where U is upper triangular
pub fn backward_substitution<T: RealField + Copy + NumericElement + LetoScalar>(
    lu_factor: &CsrMatrix<T>,
    y: &Array1<T>,
    x: &mut Array1<T>,
) {
    let n = lu_factor.nrows();

    for i in (0..n).rev() {
        let mut sum = y[i];

        let mut diag_val = <T as NumericElement>::ONE;

        let row = lu_factor.row(i);
        for (&j, &value) in row.col_indices().iter().zip(row.values()) {
            if j > i {
                sum -= value * x[j];
            } else if j == i {
                diag_val = value;
            }
        }

        x[i] = sum / diag_val;
    }
}
