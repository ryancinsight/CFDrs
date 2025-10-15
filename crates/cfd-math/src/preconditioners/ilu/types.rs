//! Core types for Incomplete LU factorization

use crate::linear_solver::Preconditioner;
use cfd_core::error::Result;
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;

use super::ilu0;
use super::iluk;
use super::triangular;

/// Incomplete LU factorization preconditioner (ILU(k))
///
/// Computes L and U factors such that L*U ≈ A with controlled fill-in.
/// The parameter k controls the level of fill: ILU(0) has the same sparsity as A,
/// while ILU(k) allows fill up to level k.
///
/// Reference: Saad, Y. (2003). Iterative Methods for Sparse Linear Systems (2nd ed.). SIAM, §10.4.
pub struct IncompleteLU<T: RealField + Copy> {
    /// Combined LU factors (L has unit diagonal)
    pub(super) lu_factor: CsrMatrix<T>,
    /// Fill level k (0 means no fill beyond original sparsity)
    fill_level: usize,
}

impl<T: RealField + Copy> IncompleteLU<T> {
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
    /// use nalgebra_sparse::CsrMatrix;
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

impl<T: RealField + Copy> Preconditioner<T> for IncompleteLU<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        let n = r.len();

        // Use workspace for intermediate result
        let mut y = DVector::zeros(n);

        // Solve LU*z = r via forward and backward substitution
        triangular::forward_substitution(&self.lu_factor, r, &mut y);
        triangular::backward_substitution(&self.lu_factor, &y, z);

        Ok(())
    }
}
