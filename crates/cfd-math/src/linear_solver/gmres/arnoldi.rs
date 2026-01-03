//! Arnoldi iteration for GMRES
//!
//! Modified Gram-Schmidt orthogonalization to build Krylov subspace basis

use crate::linear_solver::traits::{LinearOperator, Preconditioner};
use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::{DMatrix, DVector, RealField};
use num_traits::FromPrimitive;

/// Arnoldi iteration with Modified Gram-Schmidt orthogonalization
///
/// Builds an orthonormal basis for the Krylov subspace K_k(A, r0) (or K_k(M^{-1}A, r0)) using
/// Modified Gram-Schmidt for numerical stability.
///
/// # Arguments
///
/// * `a` - System operator
/// * `v` - Orthonormal basis vectors (n × (m+1))
/// * `h` - Upper Hessenberg matrix ((m+1) × m)
/// * `k` - Current iteration index
/// * `work` - Workspace vector for matrix-vector product (stores A*v_k initially, then orthogonalized)
/// * `preconditioner` - Optional preconditioner for Left Preconditioning (M^{-1})
/// * `precond_work` - Optional workspace for preconditioning result (required if preconditioner is Some)
///
/// # Returns
///
/// Norm of the newly orthogonalized vector, or error if breakdown occurs
pub fn arnoldi_iteration<T, Op, P>(
    a: &Op,
    v: &mut DMatrix<T>,
    h: &mut DMatrix<T>,
    k: usize,
    work: &mut DVector<T>,
    preconditioner: Option<&P>,
    precond_work: Option<&mut DVector<T>>,
) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
    Op: LinearOperator<T> + ?Sized,
    P: Preconditioner<T> + ?Sized,
{
    let n = work.len();

    // Extract k-th basis vector
    let v_k = v.column(k).clone_owned();

    // Compute w = A * v_k
    a.apply(&v_k, work)?;

    // Apply Left Preconditioning if provided: w <- M^{-1} * w
    if let Some(precond) = preconditioner {
        if let Some(p_work) = precond_work {
            precond.apply_to(work, p_work)?;
            work.copy_from(p_work);
        } else {
            return Err(Error::InvalidConfiguration(
                "Preconditioner workspace required but not provided".to_string(),
            ));
        }
    }

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
    h[(k + 1, k)] = norm;

    if norm > T::from_f64(1e-14).unwrap_or(T::zero()) {
        // Store normalized vector as (k+1)-th basis
        let inv_norm = T::one() / norm;
        for i in 0..n {
            v[(i, k + 1)] = work[i] * inv_norm;
        }
    }

    Ok(norm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::spmv;
    use approx::assert_relative_eq;
    use nalgebra_sparse::CsrMatrix;

    #[test]
    fn test_arnoldi_identity_matrix() {
        use crate::linear_solver::traits::Preconditioner;
        // Mock preconditioner to satisfy type inference if needed, or just use None with explicit type
        struct NoOpPrecond;
        impl Preconditioner<f64> for NoOpPrecond {
            fn apply_to(&self, _r: &DVector<f64>, _z: &mut DVector<f64>) -> Result<()> { Ok(()) }
        }
        let no_precond: Option<&NoOpPrecond> = None;

        // Test Arnoldi on identity matrix
        let n = 3;
        let a = CsrMatrix::identity(n);

        let mut v = DMatrix::zeros(n, 3);
        let mut h = DMatrix::zeros(3, 2);
        let mut work = DVector::zeros(n);

        // Initialize first basis vector
        v[(0, 0)] = 1.0;

        // First Arnoldi step: A*e_1 = e_1
        let result = arnoldi_iteration(&a, &mut v, &mut h, 0, &mut work, no_precond, None);

        // For identity: A*v_0 = v_0, so MGS gives h[0,0]=1, residual=0 (happy breakdown)
        assert!(result.is_ok()); 
        let norm = result.unwrap();
        assert!(norm < 1e-10);
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
