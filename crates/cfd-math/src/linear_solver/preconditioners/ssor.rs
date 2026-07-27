//! Symmetric Successive Over-Relaxation preconditioner.
//!
//! This is a compatibility wrapper around `leto_ops::SSORPreconditioner`
//! so `cfd-math` keeps a stable API while SSOT lives in `leto-ops`.

use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, RealField};
use leto::Array1;
use leto_ops::{
    CsrMatrix, Preconditioner as LetoPreconditioner, SSORPreconditioner, Scalar as LetoScalar,
};

/// Symmetric Successive Over-Relaxation (SSOR) preconditioner.
pub struct SSOR<T: RealField + Copy + LetoScalar> {
    inner: SSORPreconditioner<T>,
}

#[inline]
fn vector_len<T>(vector: &Array1<T>) -> usize {
    vector.shape()[0]
}

impl<T: RealField + Copy + FloatElement + LetoScalar> SSOR<T> {
    /// Create SSOR preconditioner with default relaxation parameter.
    pub fn new(matrix: CsrMatrix<T>) -> Result<Self> {
        Ok(Self {
            inner: SSORPreconditioner::new(matrix)?,
        })
    }

    /// Create SSOR preconditioner with specified relaxation parameter.
    pub fn with_omega(matrix: CsrMatrix<T>, omega: T) -> Result<Self> {
        Ok(Self {
            inner: SSORPreconditioner::with_omega(matrix, omega)?,
        })
    }
}

impl<T: RealField + Copy + FloatElement + LetoScalar> Preconditioner<T> for SSOR<T> {
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> Result<()> {
        let n = self.inner.nrows();
        let r_len = vector_len(r);
        if r_len != n {
            return Err(Error::InvalidConfiguration(format!(
                "SSOR residual length mismatch: expected {n}, got {r_len}"
            )));
        }

        let z_len = vector_len(z);
        if z_len != n {
            return Err(Error::InvalidConfiguration(format!(
                "SSOR output length mismatch: expected {n}, got {z_len}"
            )));
        }

        LetoPreconditioner::apply_to(&self.inner, r, z).map_err(Into::into)
    }
}
