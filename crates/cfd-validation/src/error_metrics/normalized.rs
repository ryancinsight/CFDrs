//! Normalized and relative error metrics

use super::{
    norms::{L1Norm, L2Norm},
    ErrorMetric,
};
use crate::scalar;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, RealField};

/// Relative error metric
pub struct RelativeError<M> {
    base_metric: M,
    tolerance: f64,
}

impl<M> RelativeError<M> {
    /// Create new relative error metric
    pub fn new(base_metric: M, tolerance: f64) -> Self {
        Self {
            base_metric,
            tolerance,
        }
    }
}

impl<T, M> ErrorMetric<T> for RelativeError<M>
where
    T: RealField + FloatElement + Copy,
    M: ErrorMetric<T>,
{
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string(),
            ));
        }

        let absolute_error = self.base_metric.compute_error(numerical, reference)?;
        let reference_norm = self
            .base_metric
            .compute_error(reference, &vec![scalar::zero(); reference.len()])?;

        let tolerance_t = scalar::from_f64(self.tolerance);
        if reference_norm < tolerance_t {
            // Reference is essentially zero, return absolute error
            Ok(absolute_error)
        } else {
            Ok(absolute_error / reference_norm)
        }
    }

    fn name(&self) -> &'static str {
        "Relative Error"
    }
}

/// Root Mean Square Error (RMSE)
pub struct RootMeanSquareError;

impl<T: RealField + FloatElement + Copy> ErrorMetric<T> for RootMeanSquareError {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        // RMSE is the same as L2 norm
        let l2 = L2Norm;
        l2.compute_error(numerical, reference)
    }

    fn name(&self) -> &'static str {
        "Root Mean Square Error"
    }
}

/// Mean Absolute Error (MAE)
pub struct MeanAbsoluteError;

impl<T: RealField + FloatElement + Copy> ErrorMetric<T> for MeanAbsoluteError {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        // MAE is the same as L1 norm
        let l1 = L1Norm;
        l1.compute_error(numerical, reference)
    }

    fn name(&self) -> &'static str {
        "Mean Absolute Error"
    }
}

/// Normalized Root Mean Square Error (NRMSE)
pub struct NormalizedRMSE {
    normalization_method: NormalizationMethod,
}

/// Methods for normalizing RMSE
#[derive(Debug, Clone, Copy)]
pub enum NormalizationMethod {
    /// Normalize by the range (max - min) of reference values
    Range,
    /// Normalize by the mean of reference values
    Mean,
    /// Normalize by the standard deviation of reference values
    StandardDeviation,
}

impl NormalizedRMSE {
    /// Create new normalized RMSE metric
    #[must_use]
    pub fn new(method: NormalizationMethod) -> Self {
        Self {
            normalization_method: method,
        }
    }
}

impl<T: RealField + FloatElement + Copy> ErrorMetric<T> for NormalizedRMSE {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string(),
            ));
        }

        if numerical.is_empty() {
            return Ok(scalar::zero::<T>());
        }

        // Compute RMSE
        let rmse = RootMeanSquareError;
        let rmse_value = rmse.compute_error(numerical, reference)?;

        // Compute normalization factor
        let normalization = match self.normalization_method {
            NormalizationMethod::Range => {
                let min_val = reference
                    .iter()
                    .fold(reference[0], |acc, x| if *x < acc { *x } else { acc });
                let max_val = reference
                    .iter()
                    .fold(reference[0], |acc, x| if *x > acc { *x } else { acc });
                max_val - min_val
            }
            NormalizationMethod::Mean => {
                let sum: T = reference
                    .iter()
                    .fold(scalar::zero::<T>(), |acc, x| acc + *x);
                let n = scalar::from_usize(reference.len());
                sum / n
            }
            NormalizationMethod::StandardDeviation => {
                // Compute mean
                let sum: T = reference
                    .iter()
                    .fold(scalar::zero::<T>(), |acc, x| acc + *x);
                let n = scalar::from_usize(reference.len());
                let mean = sum / n;

                // Compute variance
                let variance: T = reference
                    .iter()
                    .map(|x| {
                        let diff = *x - mean;
                        diff * diff
                    })
                    .fold(scalar::zero::<T>(), |acc, x| acc + x)
                    / n;

                scalar::sqrt(variance)
            }
        };

        if normalization == scalar::zero::<T>() {
            Ok(rmse_value)
        } else {
            Ok(rmse_value / normalization)
        }
    }

    fn name(&self) -> &'static str {
        "Normalized RMSE"
    }
}
