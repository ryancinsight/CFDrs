//! Core types for Incomplete LU factorization

use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, Result};
use eunomia::{NumericElement, RealField};
use leto::Array1;
use leto_ops::{CsrMatrix, Scalar as LetoScalar};

use super::ilu0;
use super::iluk;
use super::triangular;

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

/// Incomplete LU factorization preconditioner (ILU(k))
///
/// Computes L and U factors such that L*U ≈ A with controlled fill-in.
/// The parameter k controls the level of fill: ILU(0) has the same sparsity as A,
/// while ILU(k) allows fill up to level k.
///
/// Reference: Saad, Y. (2003). Iterative Methods for Sparse Linear Systems (2nd ed.). SIAM, §10.4.
pub struct IncompleteLU<T: RealField + Copy + LetoScalar> {
    /// Combined LU factors (L has unit diagonal)
    pub(super) lu_factor: CsrMatrix<T>,
    /// Fill level k (0 means no fill beyond original sparsity)
    fill_level: usize,
}

impl<T: RealField + Copy + LetoScalar> IncompleteLU<T> {
    /// Construct ILU(0) preconditioner (no fill beyond original sparsity)
    pub fn new(a: &CsrMatrix<T>) -> Result<Self> {
        Self::with_fill_level(a, 0)
    }

    /// Construct ILU(k) preconditioner with specified fill level
    ///
    /// # Arguments
    ///
    /// * `a` - Input matrix (must be square)
    /// * `k` - Fill level (0 = no fill, 1 = level-1 fill, etc.)
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use cfd_math::preconditioners::IncompleteLU;
    /// use leto_ops::CsrMatrix;
    ///
    /// // Create ILU(2) preconditioner
    /// let ilu = IncompleteLU::with_fill_level(&matrix, 2)?;
    /// ```
    pub fn with_fill_level(a: &CsrMatrix<T>, k: usize) -> Result<Self> {
        use cfd_core::error::Error;

        if a.nrows() != a.ncols() {
            return Err(Error::InvalidInput(format!(
                "Matrix must be square, got {}x{}",
                a.nrows(),
                a.ncols()
            )));
        }

        let lu_factor = if k == 0 {
            ilu0::factorize(a)?
        } else {
            iluk::factorize(a, k)?
        };

        Ok(Self {
            lu_factor,
            fill_level: k,
        })
    }

    /// Get the fill level k
    pub fn fill_level(&self) -> usize {
        self.fill_level
    }
}

impl<T: RealField + Copy + NumericElement + LetoScalar> Preconditioner<T> for IncompleteLU<T> {
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> Result<()> {
        let n = self.lu_factor.nrows();
        validate_vector_len("ILU residual", r, n)?;
        validate_vector_len("ILU output", z, n)?;

        // Use workspace for intermediate result
        let mut y = Array1::zeros([n]);
        let mut solution = Array1::zeros([n]);

        // Solve LU*z = r via forward and backward substitution
        triangular::forward_substitution(&self.lu_factor, r, &mut y);
        triangular::backward_substitution(&self.lu_factor, &y, &mut solution);

        for idx in 0..n {
            z[idx] = solution[idx];
        }

        Ok(())
    }
}
