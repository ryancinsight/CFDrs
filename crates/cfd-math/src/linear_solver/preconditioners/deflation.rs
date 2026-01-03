//! Deflation preconditioners for eigenvalue problems
//!
//! Deflation techniques improve convergence for eigenvalue problems by removing
//! known eigenvectors from the spectrum, forcing the iterative solver to focus
//! on the remaining eigenvalues.

use nalgebra::{DVector, RealField};
use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, Result};

/// Deflation preconditioner for eigenvalue problems
pub struct DeflationPreconditioner<T: RealField + Copy> {
    /// Base preconditioner (e.g., ILU, Jacobi)
    base_preconditioner: Box<dyn Preconditioner<T>>,
    /// Known eigenvectors to deflate
    eigenvectors: Vec<DVector<T>>,
    /// Corresponding eigenvalues
    eigenvalues: Vec<T>,
}

impl<T: RealField + Copy> DeflationPreconditioner<T> {
    /// Create deflation preconditioner
    ///
    /// # Arguments
    ///
    /// * `base_preconditioner` - Underlying preconditioner
    /// * `eigenvectors` - Known eigenvectors to deflate
    /// * `eigenvalues` - Corresponding eigenvalues
    ///
    /// # Returns
    ///
    /// Deflation preconditioner that removes known eigenmodes
    pub fn new(
        base_preconditioner: Box<dyn Preconditioner<T>>,
        eigenvectors: Vec<DVector<T>>,
        eigenvalues: Vec<T>,
    ) -> Result<Self> {
        if eigenvectors.len() != eigenvalues.len() {
            return Err(Error::InvalidConfiguration(
                "Number of eigenvectors must match number of eigenvalues".to_string(),
            ));
        }

        // Check eigenvector dimensions
        if let Some(first_vec) = eigenvectors.first() {
            let n = first_vec.len();
            for (i, vec) in eigenvectors.iter().enumerate() {
                if vec.len() != n {
                    return Err(Error::InvalidConfiguration(format!(
                        "Eigenvector {} has wrong dimension: expected {}, got {}",
                        i,
                        n,
                        vec.len()
                    )));
                }
            }
        }

        Ok(Self {
            base_preconditioner,
            eigenvectors,
            eigenvalues,
        })
    }

    /// Add a new eigenpair to deflate
    pub fn add_eigenpair(&mut self, eigenvector: DVector<T>, eigenvalue: T) -> Result<()> {
        // Check dimension against first eigenvector if available
        if let Some(first_vec) = self.eigenvectors.first() {
            if eigenvector.len() != first_vec.len() {
                return Err(Error::InvalidConfiguration(format!(
                    "Eigenvector dimension mismatch: expected {}, got {}",
                    first_vec.len(),
                    eigenvector.len()
                )));
            }
        }

        self.eigenvectors.push(eigenvector);
        self.eigenvalues.push(eigenvalue);
        Ok(())
    }
}

impl<T: RealField + Copy> Preconditioner<T> for DeflationPreconditioner<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        // First apply base preconditioner
        self.base_preconditioner.apply_to(r, z)?;

        // Apply deflation correction
        // z_deflated = z - sum_i (φ_i^T * z) * φ_i / λ_i
        for (eigenvec, &eigenval) in self.eigenvectors.iter().zip(&self.eigenvalues) {
            let coeff = eigenvec.dot(z) / eigenval;
            // Simple deflation: subtract projection onto eigenvector
            for i in 0..z.len() {
                z[i] -= coeff * eigenvec[i];
            }
        }

        Ok(())
    }
}
