//! Symmetric Successive Over-Relaxation preconditioner

use crate::linear_solver::Preconditioner;
use cfd_core::error::Result;
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;

// Default relaxation parameter
const DEFAULT_OMEGA: f64 = 1.0;

/// Symmetric Successive Over-Relaxation (SSOR) preconditioner
///
/// Applies forward and backward SOR sweeps for preconditioning.
pub struct SSOR<T: RealField + Copy> {
    /// System matrix
    matrix: CsrMatrix<T>,
    /// Relaxation parameter (0 < ω < 2)
    omega: T,
}

impl<T: RealField + Copy + FromPrimitive> SSOR<T> {
    /// Create SSOR preconditioner with default relaxation parameter
    pub fn new(matrix: CsrMatrix<T>) -> Result<Self> {
        Self::with_omega(matrix, T::from_f64(DEFAULT_OMEGA).unwrap_or_else(T::one))
    }

    /// Create SSOR preconditioner with specified relaxation parameter
    pub fn with_omega(matrix: CsrMatrix<T>, omega: T) -> Result<Self> {
        Ok(Self {
            matrix,
            omega,
        })
    }

    /// Forward SOR sweep
    fn forward_sweep(&self, b: &DVector<T>, x: &mut DVector<T>) {
        let n = self.matrix.nrows();

        for i in 0..n {
            let mut sum = b[i];
            let mut diag = T::one();

            let row_start = self.matrix.row_offsets()[i];
            let row_end = self.matrix.row_offsets()[i + 1];

            for idx in row_start..row_end {
                let j = self.matrix.col_indices()[idx];
                let val = self.matrix.values()[idx];

                if j < i {
                    sum -= val * x[j];
                } else if j == i {
                    diag = val;
                } else {
                    sum -= val * x[j];
                }
            }

            x[i] = (T::one() - self.omega) * x[i] + self.omega * sum / diag;
        }
    }

    /// Backward SOR sweep
    fn backward_sweep(&self, b: &DVector<T>, x: &mut DVector<T>) {
        let n = self.matrix.nrows();

        for i in (0..n).rev() {
            let mut sum = b[i];
            let mut diag = T::one();

            let row_start = self.matrix.row_offsets()[i];
            let row_end = self.matrix.row_offsets()[i + 1];

            for idx in row_start..row_end {
                let j = self.matrix.col_indices()[idx];
                let val = self.matrix.values()[idx];

                if j < i {
                    sum -= val * x[j];
                } else if j == i {
                    diag = val;
                } else {
                    sum -= val * x[j];
                }
            }

            x[i] = (T::one() - self.omega) * x[i] + self.omega * sum / diag;
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> Preconditioner<T> for SSOR<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        // Initialize z to zero for the first sweep
        z.fill(T::zero());

        // Forward sweep: solve (D + ωL)z = ωr
        self.forward_sweep(r, z);
        
        // Backward sweep: solve (D + ωU)z = ωr using result from forward sweep
        self.backward_sweep(r, z);

        Ok(())
    }
}