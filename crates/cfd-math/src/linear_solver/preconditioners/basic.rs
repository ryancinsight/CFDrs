//! Basic preconditioners (Identity, Jacobi, SOR)
//!
//! These preconditioners are simple to implement and provide a baseline
//! for more advanced techniques.

use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto_ops::{CsrMatrix, Scalar as LetoScalar};

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[inline]
fn diagonal_epsilon<T: FloatElement>() -> T {
    from_f64(1e-14)
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
        for idx in 0..vector_len(r) {
            z[idx] = r[idx];
        }
        Ok(())
    }
}

/// Jacobi (diagonal) preconditioner with memory management
pub struct JacobiPreconditioner<T: RealField + Copy> {
    inv_diagonal: Array1<T>,
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

        // Extract diagonal directly
        let diag = a.diagonal();
        let mut inv_diagonal = Array1::zeros([n]);

        let zero_tol = diagonal_epsilon();
        for i in 0..diag.len() {
            let val = diag[i];
            if NumericElement::abs(val) < zero_tol {
                // Replace zero diagonal with 1.0 (identity for this row)
                // This handles DOFs with no element contributions
                inv_diagonal[i] = <T as NumericElement>::ONE;
            } else {
                inv_diagonal[i] = <T as NumericElement>::ONE / val;
            }
        }

        Ok(Self { inv_diagonal })
    }
}

impl<T: RealField + Copy> Preconditioner<T> for JacobiPreconditioner<T> {
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> Result<()> {
        validate_vector_len("Jacobi residual", r, vector_len(&self.inv_diagonal))?;
        validate_vector_len("Jacobi output", z, vector_len(&self.inv_diagonal))?;
        for idx in 0..vector_len(r) {
            z[idx] = r[idx] * self.inv_diagonal[idx];
        }
        Ok(())
    }
}

/// SOR (Successive Over-Relaxation) preconditioner with validation.
pub struct SORPreconditioner<T: RealField + Copy + LetoScalar> {
    matrix: CsrMatrix<T>,
    pub(crate) omega: T,
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
            matrix: a.clone(),
            omega,
        })
    }

    /// Get the relaxation parameter omega
    pub fn omega(&self) -> T {
        self.omega
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
        let n = self.matrix.nrows();
        validate_vector_len("SOR residual", r, n)?;
        validate_vector_len("SOR output", z, n)?;
        for idx in 0..n {
            z[idx] = <T as NumericElement>::ZERO;
        }

        // Forward SOR sweep with in-place operations (standard SOR: (D/ω - L) z = r)
        for i in 0..n {
            let mut sum = <T as NumericElement>::ZERO;
            let row = self.matrix.row(i);
            let mut diag = <T as NumericElement>::ONE;

            // Process row entries
            for (j, val) in row.col_indices().iter().zip(row.values()) {
                if *j < i {
                    // Accumulate strictly lower part L z
                    sum += *val * z[*j];
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
