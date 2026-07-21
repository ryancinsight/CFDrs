//! Arnoldi iteration for GMRES
//!
//! Modified Gram-Schmidt orthogonalization to build Krylov subspace basis

use crate::linear_solver::traits::{LinearOperator, Preconditioner};
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::{Array1, Array2};

type KrylovBasis<T> = Array2<T>;
type Hessenberg<T> = Array2<T>;

#[inline]
fn basis_entry<T: Copy>(basis: &KrylovBasis<T>, row: usize, col: usize) -> T {
    basis[[row, col]]
}

#[inline]
fn set_basis_entry<T>(basis: &mut KrylovBasis<T>, row: usize, col: usize, value: T) {
    basis[[row, col]] = value;
}

#[inline]
fn set_hessenberg_entry<T>(hessenberg: &mut Hessenberg<T>, row: usize, col: usize, value: T) {
    hessenberg[[row, col]] = value;
}

fn fill_basis_column<T: Copy>(basis: &KrylovBasis<T>, col: usize, column: &mut Array1<T>) {
    for row in 0..column.shape()[0] {
        column[row] = basis_entry(basis, row, col);
    }
}

fn dot_work_column<T: RealField + Copy>(work: &Array1<T>, basis: &KrylovBasis<T>, col: usize) -> T {
    (0..work.shape()[0])
        .map(|row| work[row] * basis_entry(basis, row, col))
        .fold(<T as NumericElement>::ZERO, |acc, value| acc + value)
}

fn vector_norm<T: NumericElement>(work: &Array1<T>) -> T {
    let mut sum = T::ZERO;
    for row in 0..work.shape()[0] {
        sum += work[row] * work[row];
    }
    sum.sqrt()
}

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
/// * `basis_work` - Workspace vector for the extracted k-th basis column
/// * `work` - Workspace vector for matrix-vector product (stores A*v_k initially, then orthogonalized)
/// * `preconditioner` - Optional preconditioner for Left Preconditioning (M^{-1})
/// * `precond_work` - Optional workspace for preconditioning result (required if preconditioner is Some)
///
/// # Returns
///
/// Norm of the newly orthogonalized vector, or error if breakdown occurs
pub fn arnoldi_iteration<T, Op, P>(
    a: &Op,
    v: &mut KrylovBasis<T>,
    h: &mut Hessenberg<T>,
    k: usize,
    basis_work: &mut Array1<T>,
    work: &mut Array1<T>,
    preconditioner: Option<&P>,
    precond_work: Option<&mut Array1<T>>,
) -> Result<T>
where
    T: RealField + Copy + FloatElement,
    Op: LinearOperator<T> + ?Sized,
    P: Preconditioner<T> + ?Sized,
{
    let n = work.shape()[0];

    // Extract k-th basis vector
    fill_basis_column(v, k, basis_work);

    // Compute w = A * v_k
    a.apply(basis_work, work)?;

    // Apply Left Preconditioning if provided: w <- M^{-1} * w
    if let Some(precond) = preconditioner {
        if let Some(p_work) = precond_work {
            precond.apply_to(work, p_work)?;
            for i in 0..n {
                work[i] = p_work[i];
            }
        } else {
            return Err(Error::InvalidConfiguration(
                "Preconditioner workspace required but not provided".to_string(),
            ));
        }
    }

    // Modified Gram-Schmidt orthogonalization
    for j in 0..=k {
        let h_jk = dot_work_column(work, v, j);
        set_hessenberg_entry(h, j, k, h_jk);

        // w = w - h_jk * v_j
        for i in 0..n {
            work[i] -= h_jk * basis_entry(v, i, j);
        }
    }

    // Compute norm and normalize
    let norm = vector_norm(work);
    set_hessenberg_entry(h, k + 1, k, norm);

    if norm > <T as RealField>::EPSILON {
        // Store normalized vector as (k+1)-th basis
        let inv_norm = <T as NumericElement>::ONE / norm;
        for i in 0..n {
            set_basis_entry(v, i, k + 1, work[i] * inv_norm);
        }
    }

    Ok(norm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::spmv;
    use eunomia::assert_relative_eq;
    use leto::Array1;
    use leto_ops::CsrMatrix;

    #[test]
    fn test_arnoldi_identity_matrix() {
        use crate::linear_solver::traits::Preconditioner;
        // Mock preconditioner to satisfy type inference if needed, or just use None with explicit type
        struct NoOpPrecond;
        impl Preconditioner<f64> for NoOpPrecond {
            fn apply_to(&self, r: &Array1<f64>, z: &mut Array1<f64>) -> Result<()> {
                for i in 0..z.shape()[0] {
                    z[i] = r[i];
                }
                Ok(())
            }
        }
        let no_precond: Option<&NoOpPrecond> = None;

        // Test Arnoldi on identity matrix
        let n = 3;
        let a = CsrMatrix::from_parts(vec![1.0; n], (0..n).collect(), (0..=n).collect(), n, n)
            .expect("valid identity CSR");

        let mut v = KrylovBasis::zeros([n, 3]);
        let mut h = Hessenberg::zeros([3, 2]);
        let mut basis_work = Array1::zeros([n]);
        let mut work = Array1::zeros([n]);

        // Initialize first basis vector
        v[[0, 0]] = 1.0;

        // First Arnoldi step: A*e_1 = e_1
        let result = arnoldi_iteration(
            &a,
            &mut v,
            &mut h,
            0,
            &mut basis_work,
            &mut work,
            no_precond,
            None,
        );

        // For identity: A*v_0 = v_0, so MGS gives h[0,0]=1, residual=0 (happy breakdown)
        assert!(result.is_ok());
        let norm = result.unwrap();
        assert!(norm < 1e-10);
        assert_relative_eq!(h[[0, 0]], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_spmv() {
        // Test sparse matrix-vector product
        let n = 3;
        let a = CsrMatrix::from_parts(vec![1.0; n], (0..n).collect(), (0..=n).collect(), n, n)
            .expect("valid identity CSR");
        let x = Array1::from_shape_vec([n], vec![1.0, 2.0, 3.0]).unwrap();
        let mut y = Array1::zeros([n]);

        spmv(&a, &x, &mut y);

        assert_relative_eq!(y[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(y[1], 2.0, epsilon = 1e-10);
        assert_relative_eq!(y[2], 3.0, epsilon = 1e-10);
    }
}
