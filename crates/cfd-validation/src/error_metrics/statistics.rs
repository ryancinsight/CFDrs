//! Statistical error metrics and comprehensive error analysis

use super::{ErrorMetric, norms::{L1Norm, L2Norm, LInfNorm}, normalized::{MeanAbsoluteError, RootMeanSquareError, RelativeError}};
use cfd_core::error::{Error, Result};
use nalgebra::{RealField, Vector3};
use num_traits::cast::FromPrimitive;

/// Error statistics for comprehensive analysis
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ErrorStatistics<T: RealField + Copy> {
    /// L1 norm error
    pub l1_norm: T,
    /// L2 norm error
    pub l2_norm: T,
    /// L∞ norm error
    pub linf_norm: T,
    /// Mean absolute error
    pub mae: T,
    /// Root mean square error
    pub rmse: T,
    /// Relative L2 error
    pub relative_l2: T,
    /// Number of data points
    pub num_points: usize,
}

impl<T: RealField + Copy + FromPrimitive + Copy> ErrorStatistics<T> {
    /// Compute comprehensive error statistics
    pub fn compute(numerical: &[T], reference: &[T]) -> Result<Self> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string()
            ));
        }

        let l1 = L1Norm;
        let l2 = L2Norm;
        let linf = LInfNorm;
        let mae = MeanAbsoluteError;
        let rmse = RootMeanSquareError;
        let rel_l2 = RelativeError::new(L2Norm, 1e-12);

        Ok(Self {
            l1_norm: l1.compute_error(numerical, reference)?,
            l2_norm: l2.compute_error(numerical, reference)?,
            linf_norm: linf.compute_error(numerical, reference)?,
            mae: mae.compute_error(numerical, reference)?,
            rmse: rmse.compute_error(numerical, reference)?,
            relative_l2: rel_l2.compute_error(numerical, reference)?,
            num_points: numerical.len(),
        })
    }

    /// Compute error statistics for vector fields
    pub fn compute_vector(numerical: &[Vector3<T>], reference: &[Vector3<T>]) -> Result<Self> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string()
            ));
        }

        // Extract magnitudes
        let num_magnitudes: Vec<T> = numerical.iter().map(nalgebra::Matrix::norm).collect();
        let ref_magnitudes: Vec<T> = reference.iter().map(nalgebra::Matrix::norm).collect();

        Self::compute(&num_magnitudes, &ref_magnitudes)
    }
}