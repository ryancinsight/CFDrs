//! Basic preconditioners (Identity, Jacobi, SOR)
//!
//! These preconditioners are simple to implement and provide a baseline
//! for more advanced techniques.

use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto_ops::{
    CsrMatrix, IdentityPreconditioner as LetoIdentityPreconditioner,
    JacobiPreconditioner as LetoJacobiPreconditioner, Preconditioner as LetoPreconditioner,
    SORPreconditioner as LetoSORPreconditioner, Scalar as LetoScalar,
};

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[inline]
fn vector_len<T>(vector: &Array1<T>) -> usize {
    vector.shape()[0]
}

fn validate_vector_len<T>(name: &str, vector: &Array1<T>, expected: usize) -> Result<()> {
    let actual = vector_len(vector);
    if actual != expected {
        return Err(Error::InvalidConfiguration(format!(
            "{name} length mismatch: expected {expected}, got {actual}"
        )));
    }
    Ok(())
}

/// Identity preconditioner (no preconditioning)
#[derive(Default)]
pub struct IdentityPreconditioner;

impl<T: RealField + Copy> Preconditioner<T> for IdentityPreconditioner {
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> Result<()> {
        validate_vector_len("identity preconditioner output", z, vector_len(r))?;
        LetoPreconditioner::apply_to(&LetoIdentityPreconditioner, r, z).map_err(Into::into)
    }
}

/// Jacobi (diagonal) preconditioner with memory management
pub struct JacobiPreconditioner<T: RealField + Copy> {
    inner: LetoJacobiPreconditioner<T>,
}

impl<T: RealField + FloatElement + Copy + LetoScalar> JacobiPreconditioner<T> {
    /// Create Jacobi preconditioner from matrix diagonal
    ///
    /// # Zero Diagonal Handling
    ///
    /// When a diagonal entry is zero or near-zero (< 1e-14), it is replaced with 1.0
    /// to avoid division by zero. This effectively makes the preconditioner act as
    /// identity for those rows, which is appropriate for:
    /// - DOFs with no element contributions (e.g., unused nodes)
    /// - Pressure DOFs in mixed formulations before stabilization
    pub fn new(a: &CsrMatrix<T>) -> Result<Self> {
        let n = a.nrows();
        if n != a.ncols() {
            return Err(Error::InvalidConfiguration(
                "Jacobi preconditioner requires square matrix".to_string(),
            ));
        }

        Ok(Self {
            inner: LetoJacobiPreconditioner::from_matrix_identity_on_zero(a),
        })
    }
}

impl<T: RealField + FloatElement + Copy + LetoScalar> Preconditioner<T> for JacobiPreconditioner<T> {
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> Result<()> {
        let n = self.inner.nrows();
        validate_vector_len("Jacobi residual", r, n)?;
        validate_vector_len("Jacobi output", z, n)?;
        LetoPreconditioner::apply_to(&self.inner, r, z).map_err(Into::into)
    }
}

/// SOR (Successive Over-Relaxation) preconditioner with validation.
pub struct SORPreconditioner<T: RealField + Copy + LetoScalar> {
    inner: LetoSORPreconditioner<T>,
}

impl<T: RealField + FloatElement + Copy + LetoScalar> SORPreconditioner<T> {
    /// Create SOR preconditioner with specified relaxation parameter
    pub fn new(a: &CsrMatrix<T>, omega: T) -> Result<Self> {
        let n = a.nrows();
        if n != a.ncols() {
            return Err(Error::InvalidConfiguration(
                "SOR preconditioner requires square matrix".to_string(),
            ));
        }

        // Validate omega range for stability
        let zero = <T as NumericElement>::ZERO;
        let two = from_f64(2.0);
        if omega <= zero || omega >= two {
            return Err(Error::InvalidConfiguration(
                "SOR omega parameter must be in range (0, 2) for stability".to_string(),
            ));
        }

        Ok(Self {
            inner: LetoSORPreconditioner::new(a.clone(), omega)?,
        })
    }

    /// Get the relaxation parameter omega
    pub fn omega(&self) -> T {
        self.inner.omega()
    }

    /// Create SOR preconditioner with omega tuned for 1D Poisson problems
    pub fn with_omega_for_1d_poisson(a: &CsrMatrix<T>) -> Result<Self> {
        // Validate matrix structure for 1D Poisson
        Self::validate_1d_poisson_structure(a)?;

        let n = a.nrows() as f64;
        let omega_opt = 2.0 / (1.0 + (std::f64::consts::PI / n).sin());
        let omega = from_f64(omega_opt);

        Self::new(a, omega)
    }

    /// Validate that matrix has 1D Poisson structure (tridiagonal with specific pattern)
    fn validate_1d_poisson_structure(a: &CsrMatrix<T>) -> Result<()> {
        let n = a.nrows();

        // Check structure: each row should have at most 3 non-zeros
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
                        "Non-zero at ({i}, {j}) violates tridiagonal structure"
                    )));
                }
            }
        }

        Ok(())
    }
}

impl<T: RealField + Copy + NumericElement + LetoScalar> Preconditioner<T> for SORPreconditioner<T> {
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> Result<()> {
        let n = self.inner.nrows();
        validate_vector_len("SOR residual", r, n)?;
        validate_vector_len("SOR output", z, n)?;
        LetoPreconditioner::apply_to(&self.inner, r, z).map_err(Into::into)
    }
}
