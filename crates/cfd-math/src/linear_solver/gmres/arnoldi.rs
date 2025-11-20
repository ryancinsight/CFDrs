//! Arnoldi iteration for GMRES
//!
//! Modified Gram-Schmidt orthogonalization to build Krylov subspace basis

use crate::sparse::spmv;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::{DMatrix, DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;

/// Arnoldi iteration with Modified Gram-Schmidt orthogonalization
///
/// Builds an orthonormal basis for the Krylov subspace K_k(A, r0) using
/// Modified Gram-Schmidt for numerical stability.
///
/// # Arguments
///
/// * `a` - System matrix
/// * `v` - Orthonormal basis vectors (n × (m+1))
/// * `h` - Upper Hessenberg matrix ((m+1) × m)
/// * `k` - Current iteration index
/// * `work` - Workspace vector for matrix-vector product
///
/// # Returns
///
/// Norm of the newly orthogonalized vector, or error if breakdown occurs
///
/// # Algorithm
///
/// 1. Compute w = A * v_k
/// 2. For j = 0..k: orthogonalize w against v_j using MGS
/// 3. Normalize w to obtain v_{k+1}
/// 4. Store coefficients in Hessenberg matrix H
pub fn arnoldi_iteration<T: RealField + Copy + FromPrimitive>(
    a: &CsrMatrix<T>,
    v: &mut DMatrix<T>,
    h: &mut DMatrix<T>,
    k: usize,
    work: &mut DVector<T>,
) -> Result<T> {
    let n = a.nrows();

    // Extract k-th basis vector
    let v_k = v.column(k).clone_owned();

    // Compute w = A * v_k via efficient sparse matrix-vector product
    spmv(a, &v_k, work);

    // Modified Gram-Schmidt orthogonalization
    for j in 0..=k {
        let v_j = v.column(j);
        let h_jk = work.dot(&v_j);
        h[(j, k)] = h_jk;

        // w = w - h_jk * v_j
        for i in 0..n {
            work[i] -= h_jk * v_j[i];
        }
    }

    // Compute norm and normalize
    let norm = work.norm();
    if norm < T::from_f64(1e-14).unwrap_or(T::zero()) {
        // Breakdown: residual is in Krylov subspace (happy breakdown or singular matrix)
        return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
    }

    h[(k + 1, k)] = norm;

    // Store normalized vector as (k+1)-th basis
    let inv_norm = T::one() / norm;
    for i in 0..n {
        v[(i, k + 1)] = work[i] * inv_norm;
    }

    Ok(norm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_arnoldi_identity_matrix() {
        // Test Arnoldi on identity matrix
        let n = 3;
        let a = CsrMatrix::identity(n);

        let mut v = DMatrix::zeros(n, 3);
        let mut h = DMatrix::zeros(3, 2);
        let mut work = DVector::zeros(n);

        // Initialize first basis vector
        v[(0, 0)] = 1.0;

        // First Arnoldi step: A*e_1 = e_1
        let result = arnoldi_iteration(&a, &mut v, &mut h, 0, &mut work);

        // For identity: A*v_0 = v_0, so MGS gives h[0,0]=1, residual=0 (happy breakdown)
        assert!(result.is_err()); // Expect happy breakdown (singular matrix error)
        assert_relative_eq!(h[(0, 0)], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_spmv() {
        // Test sparse matrix-vector product
        let n = 3;
        let a = CsrMatrix::identity(n);
        let x = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut y = DVector::zeros(n);

        spmv(&a, &x, &mut y);

        assert_relative_eq!(y[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(y[1], 2.0, epsilon = 1e-10);
        assert_relative_eq!(y[2], 3.0, epsilon = 1e-10);
    }
}
