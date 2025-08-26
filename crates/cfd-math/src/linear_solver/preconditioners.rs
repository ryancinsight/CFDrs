//! Preconditioner implementations for iterative solvers

use super::traits::Preconditioner;
use crate::sparse::SparseMatrixExt;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;

/// Identity preconditioner (no preconditioning)
#[derive(Default)]
pub struct IdentityPreconditioner;

impl<T: RealField + Copy> Preconditioner<T> for IdentityPreconditioner {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        z.copy_from(r);
        Ok(())
    }
}

/// Jacobi (diagonal) preconditioner with optimized memory usage
pub struct JacobiPreconditioner<T: RealField + Copy> {
    inv_diagonal: DVector<T>,
}

impl<T: RealField + From<f64> + FromPrimitive + Copy> JacobiPreconditioner<T> {
    /// Create Jacobi preconditioner from matrix diagonal
    pub fn new(a: &CsrMatrix<T>) -> Result<Self> {
        let n = a.nrows();
        if n != a.ncols() {
            return Err(Error::InvalidConfiguration(
                "Jacobi preconditioner requires square matrix".to_string(),
            ));
        }

        // Extract diagonal directly
        let diag = a.diagonal();
        let mut inv_diagonal = DVector::zeros(n);

        for (i, val) in diag.iter().enumerate() {
            if val.abs() < T::from_f64(1e-14).unwrap_or_else(|| T::zero()) {
                return Err(Error::Numerical(NumericalErrorKind::DivisionByZero));
            }
            inv_diagonal[i] = T::one() / *val;
        }

        Ok(Self { inv_diagonal })
    }
}

impl<T: RealField + Copy> Preconditioner<T> for JacobiPreconditioner<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        z.copy_from(r);
        z.component_mul_assign(&self.inv_diagonal);
        Ok(())
    }
}

/// SOR (Successive Over-Relaxation) preconditioner with validation
pub struct SORPreconditioner<T: RealField + Copy> {
    matrix: CsrMatrix<T>,
    pub(crate) omega: T,
}

impl<T: RealField + From<f64> + FromPrimitive + Copy> SORPreconditioner<T> {
    /// Create SOR preconditioner with specified relaxation parameter
    pub fn new(a: &CsrMatrix<T>, omega: T) -> Result<Self> {
        let n = a.nrows();
        if n != a.ncols() {
            return Err(Error::InvalidConfiguration(
                "SOR preconditioner requires square matrix".to_string(),
            ));
        }

        // Validate omega range for stability
        let zero = T::zero();
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        if omega <= zero || omega >= two {
            return Err(Error::InvalidConfiguration(
                "SOR omega parameter must be in range (0, 2) for stability".to_string(),
            ));
        }

        Ok(Self {
            matrix: a.clone(),
            omega,
        })
    }

    /// Get the relaxation parameter omega
    pub fn omega(&self) -> T {
        self.omega
    }

    /// Create SOR preconditioner with omega optimized for 1D Poisson problems
    ///
    /// ## Warning
    /// This function computes an optimal omega value **specifically** for matrices
    /// arising from 1D Poisson equations discretized with second-order finite
    /// differences on uniform grids. Using this for any other matrix type will
    /// likely result in suboptimal performance or convergence issues.
    ///
    /// ## Validation
    /// This function validates that the matrix has the expected tridiagonal
    /// structure before computing the optimal omega.
    pub fn with_omega_for_1d_poisson(a: &CsrMatrix<T>) -> Result<Self> {
        // Validate matrix structure for 1D Poisson
        Self::validate_1d_poisson_structure(a)?;

        let n = a.nrows() as f64;
        let omega_opt = 2.0 / (1.0 + (std::f64::consts::PI / n).sin());
        let omega = T::from_f64(omega_opt).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::ConversionFailed {
                from_type: "f64",
                to_type: std::any::type_name::<T>(),
                value: omega_opt.to_string(),
            })
        })?;

        Self::new(a, omega)
    }

    /// Validate that matrix has 1D Poisson structure (tridiagonal with specific pattern)
    fn validate_1d_poisson_structure(a: &CsrMatrix<T>) -> Result<()> {
        let n = a.nrows();

        // Check basic structure: each row should have at most 3 non-zeros
        for i in 0..n {
            let row = a.row(i);
            if row.nnz() > 3 {
                return Err(Error::InvalidConfiguration(format!(
                    "Row {} has {} non-zeros; 1D Poisson should have at most 3",
                    i,
                    row.nnz()
                )));
            }

            // Check that non-zeros are in expected positions (diagonal and adjacent)
            for &j in row.col_indices() {
                if (j as i32 - i as i32).abs() > 1 {
                    return Err(Error::InvalidConfiguration(format!(
                        "Non-zero at ({}, {}) violates tridiagonal structure",
                        i, j
                    )));
                }
            }
        }

        Ok(())
    }
}

impl<T: RealField + Copy> Preconditioner<T> for SORPreconditioner<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        let n = r.len();
        z.fill(T::zero());

        // Forward SOR sweep with in-place operations (standard SOR: (D/ω - L) z = r)
        for i in 0..n {
            let mut sum = T::zero();
            let row = self.matrix.row(i);
            let mut diag = T::one();

            // Process row entries
            for (j, val) in row.col_indices().iter().zip(row.values()) {
                if *j < i {
                    // Accumulate strictly lower part L z
                    sum = sum + *val * z[*j];
                } else if *j == i {
                    diag = *val;
                }
            }

            // Standard SOR preconditioner solves (D/ω - L) z = r
            // => (diag/ω) z_i = r_i + (L z)_i
            z[i] = (r[i] + sum) * self.omega / diag;
        }

        Ok(())
    }
}
