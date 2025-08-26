//! Error metrics for CFD validation.
//!
//! This module provides various error metrics to quantify the difference
//! between numerical solutions and analytical or reference solutions.

mod analysis;
mod normalized;
mod norms;
mod statistics;

#[cfg(test)]
mod tests;

// Re-export main components
pub use analysis::ErrorAnalysis;
pub use normalized::{
    MeanAbsoluteError, NormalizationMethod, NormalizedRMSE, RelativeError, RootMeanSquareError,
};
pub use norms::{L1Norm, L2Norm, LInfNorm};
pub use statistics::ErrorStatistics;

use cfd_core::error::Result;
use nalgebra::{RealField, Vector3};

/// Trait for error metrics
pub trait ErrorMetric<T: RealField + Copy> {
    /// Compute the error between numerical and reference solutions
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T>;

    /// Compute the error for vector fields
    fn compute_vector_error(
        &self,
        numerical: &[Vector3<T>],
        reference: &[Vector3<T>],
    ) -> Result<T> {
        use cfd_core::error::Error;

        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Numerical and reference solutions must have the same length".to_string(),
            ));
        }

        // Extract magnitude for each vector
        let num_magnitudes: Vec<T> = numerical.iter().map(nalgebra::Matrix::norm).collect();
        let ref_magnitudes: Vec<T> = reference.iter().map(nalgebra::Matrix::norm).collect();

        self.compute_error(&num_magnitudes, &ref_magnitudes)
    }

    /// Get the name of the error metric
    fn name(&self) -> &str;
}
