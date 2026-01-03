//! ILU(0) factorization - no fill beyond original sparsity
//!
//! Reference: Saad, Y. (2003). Iterative Methods for Sparse Linear Systems (2nd ed.). SIAM, §10.4.

use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::RealField;
use nalgebra_sparse::CsrMatrix;

use super::utils;

/// Perform ILU(0) factorization (original algorithm, optimized)
///
/// ILU(0) maintains the same sparsity pattern as the input matrix A.
/// The algorithm performs in-place factorization without adding any fill-in.
pub fn factorize<T: RealField + Copy>(a: &CsrMatrix<T>) -> Result<CsrMatrix<T>> {
    // Create a mutable copy of the matrix
    let mut lu_vals = a.values().to_vec();
    let lu_indices = a.col_indices().to_vec();
    let lu_offsets = a.row_offsets().to_vec();

    let n = a.nrows();

    // ILU(0) factorization
    for k in 0..n - 1 {
        let diag_idx = utils::find_diagonal_index(&lu_offsets, &lu_indices, k)?;
        let a_kk = lu_vals[diag_idx];

        if a_kk.abs() <= T::default_epsilon() {
            return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
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
                                let k_val = lu_vals[k_idx];
                                lu_vals[j_idx] -= mult * k_val;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    CsrMatrix::try_from_csr_data(n, n, lu_offsets, lu_indices, lu_vals).map_err(|_| {
        Error::InvalidInput("Failed to create CSR matrix from ILU(0) factorization".to_string())
    })
}
