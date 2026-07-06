//! Symmetric Successive Over-Relaxation preconditioner

use crate::linear_solver::Preconditioner;
use cfd_core::error::Result;
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto_ops::{CsrMatrix, Scalar as LetoScalar};

// Default relaxation parameter
const DEFAULT_OMEGA: f64 = 1.0;

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
        return Err(cfd_core::error::Error::InvalidConfiguration(format!(
            "{name} length mismatch: expected {expected}, got {actual}"
        )));
    }
    Ok(())
}

/// Symmetric Successive Over-Relaxation (SSOR) preconditioner
///
/// Applies forward and backward SOR sweeps for preconditioning.
///
/// # Theorem (SSOR Convergence and Symmetry)
///
/// For an SPD matrix $A = D + L + L^T$ (diagonal $D$, strict lower triangle $L$)
/// and relaxation parameter $\omega \in (0, 2)$, the SSOR preconditioner
///
/// $$M_{\text{SSOR}} = \frac{1}{\omega(2 - \omega)}(D + \omega L)\,D^{-1}(D + \omega L^T)$$
///
/// is symmetric positive definite, and the preconditioned iteration
/// $M_{\text{SSOR}}^{-1} A$ has spectral radius $\rho < 1$.
///
/// **Proof sketch**: The forward sweep $(D + \omega L) x^{*} = \omega b - [\omega L^T + (1-\omega)D] x^{(k)}$
/// is composed with its transpose (backward sweep) to produce a symmetric
/// splitting. SPD follows from $M = B^T D^{-1} B$ with $B = D + \omega L$
/// non-singular. The Ostrowski–Reich theorem guarantees $\rho < 1$ for
/// $\omega \in (0, 2)$ when $A$ is SPD.
///
/// **Reference**: Saad (2003) §4.2; Young (1971) §5.3.
pub struct SSOR<T: RealField + Copy + LetoScalar> {
    /// System matrix
    matrix: CsrMatrix<T>,
    /// Relaxation parameter (0 < ω < 2)
    omega: T,
}

impl<T: RealField + Copy + FloatElement + LetoScalar> SSOR<T> {
    /// Create SSOR preconditioner with default relaxation parameter
    pub fn new(matrix: CsrMatrix<T>) -> Result<Self> {
        Self::with_omega(matrix, from_f64(DEFAULT_OMEGA))
    }

    /// Create SSOR preconditioner with specified relaxation parameter
    pub fn with_omega(matrix: CsrMatrix<T>, omega: T) -> Result<Self> {
        Ok(Self { matrix, omega })
    }

    /// Forward SOR sweep
    fn forward_sweep(&self, b: &Array1<T>, x: &mut Array1<T>) {
        let n = self.matrix.nrows();

        for i in 0..n {
            let mut sum = b[i];
            let mut diag = <T as NumericElement>::ONE;

            let row = self.matrix.row(i);

            for (&j, &val) in row.col_indices().iter().zip(row.values()) {
                match j.cmp(&i) {
                    std::cmp::Ordering::Less | std::cmp::Ordering::Greater => {
                        sum -= val * x[j];
                    }
                    std::cmp::Ordering::Equal => {
                        diag = val;
                    }
                }
            }

            x[i] = (<T as NumericElement>::ONE - self.omega) * x[i] + self.omega * sum / diag;
        }
    }

    /// Backward SOR sweep
    fn backward_sweep(&self, b: &Array1<T>, x: &mut Array1<T>) {
        let n = self.matrix.nrows();

        for i in (0..n).rev() {
            let mut sum = b[i];
            let mut diag = <T as NumericElement>::ONE;

            let row = self.matrix.row(i);

            for (&j, &val) in row.col_indices().iter().zip(row.values()) {
                match j.cmp(&i) {
                    std::cmp::Ordering::Less | std::cmp::Ordering::Greater => {
                        sum -= val * x[j];
                    }
                    std::cmp::Ordering::Equal => {
                        diag = val;
                    }
                }
            }

            x[i] = (<T as NumericElement>::ONE - self.omega) * x[i] + self.omega * sum / diag;
        }
    }
}

impl<T: RealField + Copy + FloatElement + LetoScalar> Preconditioner<T> for SSOR<T> {
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> Result<()> {
        validate_vector_len("SSOR residual", r, self.matrix.nrows())?;
        validate_vector_len("SSOR output", z, self.matrix.nrows())?;

        // Initialize z to zero for the first sweep
        for idx in 0..vector_len(z) {
            z[idx] = <T as NumericElement>::ZERO;
        }

        // Forward sweep: solve (D + ωL)z = ωr
        self.forward_sweep(r, z);

        // Backward sweep: solve (D + ωU)z = ωr using result from forward sweep
        self.backward_sweep(r, z);

        Ok(())
    }
}
