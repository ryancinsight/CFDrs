//! Deflation preconditioners for eigenvalue problems
//!
//! Deflation techniques improve convergence for eigenvalue problems by removing
//! known eigenvectors from the spectrum, forcing the iterative solver to focus
//! on the remaining eigenvalues.

use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, Result};
use eunomia::{NumericElement, RealField};
use leto::Array1;

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

fn dot_arrays<T: RealField + Copy + NumericElement>(lhs: &Array1<T>, rhs: &Array1<T>) -> T {
    let mut sum = <T as NumericElement>::ZERO;
    for idx in 0..vector_len(lhs) {
        sum += lhs[idx] * rhs[idx];
    }
    sum
}

/// Deflation preconditioner for eigenvalue problems
pub struct DeflationPreconditioner<T: RealField + Copy> {
    /// Base preconditioner (e.g., ILU, Jacobi)
    base_preconditioner: Box<dyn Preconditioner<T>>,
    /// Known eigenvectors to deflate
    eigenvectors: Vec<Array1<T>>,
    /// Corresponding eigenvalues
    eigenvalues: Vec<T>,
}

impl<T: RealField + Copy + NumericElement> DeflationPreconditioner<T> {
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
        eigenvectors: Vec<Array1<T>>,
        eigenvalues: Vec<T>,
    ) -> Result<Self> {
        if eigenvectors.len() != eigenvalues.len() {
            return Err(Error::InvalidConfiguration(
                "Number of eigenvectors must match number of eigenvalues".to_string(),
            ));
        }

        // Check eigenvector dimensions
        if let Some(first_vec) = eigenvectors.first() {
            let n = vector_len(first_vec);
            for (i, vec) in eigenvectors.iter().enumerate() {
                if vector_len(vec) != n {
                    return Err(Error::InvalidConfiguration(format!(
                        "Eigenvector {} has wrong dimension: expected {}, got {}",
                        i,
                        n,
                        vector_len(vec)
                    )));
                }
            }
        }

        for (idx, eigenvalue) in eigenvalues.iter().enumerate() {
            if *eigenvalue == <T as NumericElement>::ZERO {
                return Err(Error::InvalidConfiguration(format!(
                    "Eigenvalue {idx} is zero; deflation requires nonzero eigenvalues"
                )));
            }
        }

        Ok(Self {
            base_preconditioner,
            eigenvectors,
            eigenvalues,
        })
    }

    /// Add a new eigenpair to deflate
    pub fn add_eigenpair(&mut self, eigenvector: Array1<T>, eigenvalue: T) -> Result<()> {
        // Check dimension against first eigenvector if available
        if let Some(first_vec) = self.eigenvectors.first() {
            if vector_len(&eigenvector) != vector_len(first_vec) {
                return Err(Error::InvalidConfiguration(format!(
                    "Eigenvector dimension mismatch: expected {}, got {}",
                    vector_len(first_vec),
                    vector_len(&eigenvector)
                )));
            }
        }

        if eigenvalue == <T as NumericElement>::ZERO {
            return Err(Error::InvalidConfiguration(
                "Eigenvalue is zero; deflation requires nonzero eigenvalues".to_string(),
            ));
        }

        self.eigenvectors.push(eigenvector);
        self.eigenvalues.push(eigenvalue);
        Ok(())
    }
}

impl<T: RealField + Copy + NumericElement> Preconditioner<T> for DeflationPreconditioner<T> {
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> Result<()> {
        validate_vector_len("deflation output", z, vector_len(r))?;
        for (idx, eigenvector) in self.eigenvectors.iter().enumerate() {
            validate_vector_len(
                &format!("deflation eigenvector {idx}"),
                eigenvector,
                vector_len(r),
            )?;
        }

        // First apply base preconditioner
        let mut base_solution = Array1::from_elem(r.shape(), <T as NumericElement>::ZERO);
        self.base_preconditioner.apply_to(r, &mut base_solution)?;
        for idx in 0..vector_len(r) {
            z[idx] = base_solution[idx];
        }

        // Apply deflation correction
        // z_deflated = z - sum_i (φ_i^T * z) * φ_i / λ_i
        for (eigenvec, &eigenval) in self.eigenvectors.iter().zip(&self.eigenvalues) {
            let coeff = dot_arrays(eigenvec, z) / eigenval;
            // Simple deflation: subtract projection onto eigenvector
            for i in 0..vector_len(z) {
                z[i] -= coeff * eigenvec[i];
            }
        }

        Ok(())
    }
}
