//! Error metrics for CFD validation.
//!
//! This module provides various error metrics to quantify the difference
//! between numerical solutions and analytical or reference solutions.

use cfd_core::error::{Error, Result};
use nalgebra::{RealField, Vector3};
use num_traits::cast::FromPrimitive;

/// Trait for error metrics
pub trait ErrorMetric<T: RealField> {
    /// Compute the error between numerical and reference solutions
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T>;

    /// Compute the error for vector fields
    fn compute_vector_error(&self, numerical: &[Vector3<T>], reference: &[Vector3<T>]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Numerical and reference solutions must have the same length".to_string()
            ));
        }

        // Extract magnitude for each vector
        let num_magnitudes: Vec<T> = numerical.iter().map(|v| v.norm()).collect();
        let ref_magnitudes: Vec<T> = reference.iter().map(|v| v.norm()).collect();

        self.compute_error(&num_magnitudes, &ref_magnitudes)
    }

    /// Get the name of the error metric
    fn name(&self) -> &str;
}

/// L2 (Euclidean) norm error metric
pub struct L2Norm;

impl<T: RealField + FromPrimitive> ErrorMetric<T> for L2Norm {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string()
            ));
        }

        if numerical.is_empty() {
            return Ok(T::zero());
        }

        let sum_squared_diff: T = numerical
            .iter()
            .zip(reference.iter())
            .map(|(num, ref_val)| {
                let diff = num.clone() - ref_val.clone();
                diff.clone() * diff
            })
            .fold(T::zero(), |acc, x| acc + x);

        let n = T::from_usize(numerical.len()).ok_or_else(|| {
            Error::InvalidConfiguration("Failed to convert array length to target type".to_string())
        })?;
        Ok((sum_squared_diff / n).sqrt())
    }

    fn name(&self) -> &str {
        "L2 Norm"
    }
}

/// L∞ (maximum) norm error metric
pub struct LInfNorm;

impl<T: RealField> ErrorMetric<T> for LInfNorm {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string()
            ));
        }

        if numerical.is_empty() {
            return Ok(T::zero());
        }

        let max_diff = numerical
            .iter()
            .zip(reference.iter())
            .map(|(num, ref_val)| (num.clone() - ref_val.clone()).abs())
            .fold(T::zero(), |acc, x| if x > acc { x } else { acc });

        Ok(max_diff)
    }

    fn name(&self) -> &str {
        "L∞ Norm"
    }
}

/// L1 (Manhattan) norm error metric
pub struct L1Norm;

impl<T: RealField + FromPrimitive> ErrorMetric<T> for L1Norm {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string()
            ));
        }

        if numerical.is_empty() {
            return Ok(T::zero());
        }

        let sum_abs_diff: T = numerical
            .iter()
            .zip(reference.iter())
            .map(|(num, ref_val)| (num.clone() - ref_val.clone()).abs())
            .fold(T::zero(), |acc, x| acc + x);

        let n = T::from_usize(numerical.len()).ok_or_else(|| {
            Error::InvalidConfiguration("Failed to convert array length to target type".to_string())
        })?;
        Ok(sum_abs_diff / n)
    }

    fn name(&self) -> &str {
        "L1 Norm"
    }
}

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
    T: RealField + FromPrimitive,
    M: ErrorMetric<T>,
{
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string()
            ));
        }

        let absolute_error = self.base_metric.compute_error(numerical, reference)?;
        let reference_norm = self.base_metric.compute_error(reference, &vec![T::zero(); reference.len()])?;

        let tolerance_t = T::from_f64(self.tolerance).ok_or_else(|| {
            Error::InvalidConfiguration("Failed to convert tolerance to target type".to_string())
        })?;
        if reference_norm < tolerance_t {
            // Reference is essentially zero, return absolute error
            Ok(absolute_error)
        } else {
            Ok(absolute_error / reference_norm)
        }
    }

    fn name(&self) -> &str {
        "Relative Error"
    }
}

/// Root Mean Square Error (RMSE)
pub struct RootMeanSquareError;

impl<T: RealField + FromPrimitive> ErrorMetric<T> for RootMeanSquareError {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        // RMSE is the same as L2 norm
        let l2 = L2Norm;
        l2.compute_error(numerical, reference)
    }

    fn name(&self) -> &str {
        "Root Mean Square Error"
    }
}

/// Mean Absolute Error (MAE)
pub struct MeanAbsoluteError;

impl<T: RealField + FromPrimitive> ErrorMetric<T> for MeanAbsoluteError {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        // MAE is the same as L1 norm
        let l1 = L1Norm;
        l1.compute_error(numerical, reference)
    }

    fn name(&self) -> &str {
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
    pub fn new(method: NormalizationMethod) -> Self {
        Self {
            normalization_method: method,
        }
    }
}

impl<T: RealField + FromPrimitive> ErrorMetric<T> for NormalizedRMSE {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string()
            ));
        }

        if numerical.is_empty() {
            return Ok(T::zero());
        }

        // Compute RMSE
        let rmse = RootMeanSquareError;
        let rmse_value = rmse.compute_error(numerical, reference)?;

        // Compute normalization factor
        let normalization = match self.normalization_method {
            NormalizationMethod::Range => {
                let min_val = reference.iter().fold(reference[0].clone(), |acc, x| {
                    if *x < acc { x.clone() } else { acc }
                });
                let max_val = reference.iter().fold(reference[0].clone(), |acc, x| {
                    if *x > acc { x.clone() } else { acc }
                });
                max_val - min_val
            },
            NormalizationMethod::Mean => {
                let sum: T = reference.iter().fold(T::zero(), |acc, x| acc + x.clone());
                let n = T::from_usize(reference.len()).expect("CRITICAL: Add proper error handling");
                sum / n
            },
            NormalizationMethod::StandardDeviation => {
                // Compute mean
                let sum: T = reference.iter().fold(T::zero(), |acc, x| acc + x.clone());
                let n = T::from_usize(reference.len()).expect("CRITICAL: Add proper error handling");
                let mean = sum / n.clone();

                // Compute variance
                let variance: T = reference.iter()
                    .map(|x| {
                        let diff = x.clone() - mean.clone();
                        diff.clone() * diff
                    })
                    .fold(T::zero(), |acc, x| acc + x) / n;

                variance.sqrt()
            },
        };

        if normalization == T::zero() {
            Ok(rmse_value)
        } else {
            Ok(rmse_value / normalization)
        }
    }

    fn name(&self) -> &str {
        "Normalized RMSE"
    }
}

/// Error statistics for comprehensive analysis
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ErrorStatistics<T: RealField> {
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

impl<T: RealField + FromPrimitive> ErrorStatistics<T> {
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
        let num_magnitudes: Vec<T> = numerical.iter().map(|v| v.norm()).collect();
        let ref_magnitudes: Vec<T> = reference.iter().map(|v| v.norm()).collect();

        Self::compute(&num_magnitudes, &ref_magnitudes)
    }
}

/// Utility functions for error analysis
pub struct ErrorAnalysis;

impl ErrorAnalysis {
    /// Compute convergence rate from error measurements at different grid sizes
    pub fn convergence_rate<T: RealField + FromPrimitive>(
        grid_sizes: &[T],
        errors: &[T],
    ) -> Result<T> {
        if grid_sizes.len() != errors.len() || grid_sizes.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 grid sizes and corresponding errors".to_string()
            ));
        }

        // Use least squares fit to log(error) = log(C) + p*log(h)
        // where p is the convergence rate
        let n = T::from_usize(grid_sizes.len()).expect("CRITICAL: Add proper error handling");

        let log_h: Vec<T> = grid_sizes.iter().map(|h| h.clone().ln()).collect();
        let log_e: Vec<T> = errors.iter().map(|e| e.clone().ln()).collect();

        // Compute means
        let mean_log_h = log_h.iter().fold(T::zero(), |acc, x| acc + x.clone()) / n.clone();
        let mean_log_e = log_e.iter().fold(T::zero(), |acc, x| acc + x.clone()) / n.clone();

        // Compute slope (convergence rate)
        let numerator: T = log_h.iter().zip(log_e.iter())
            .map(|(h, e)| (h.clone() - mean_log_h.clone()) * (e.clone() - mean_log_e.clone()))
            .fold(T::zero(), |acc, x| acc + x);

        let denominator: T = log_h.iter()
            .map(|h| {
                let diff = h.clone() - mean_log_h.clone();
                diff.clone() * diff
            })
            .fold(T::zero(), |acc, x| acc + x);

        if denominator == T::zero() {
            return Err(Error::InvalidConfiguration(
                "Cannot compute convergence rate: grid sizes are identical".to_string()
            ));
        }

        Ok(numerator / denominator)
    }

    /// Check if error is within acceptable tolerance
    pub fn is_acceptable<T: RealField>(error: T, tolerance: T) -> bool {
        error <= tolerance
    }

    /// Compute error reduction factor between two measurements
    pub fn error_reduction_factor<T: RealField>(coarse_error: T, fine_error: T) -> T {
        if fine_error == T::zero() {
            return T::max_value().unwrap_or_else(|| T::from_f64(1e10).unwrap_or_else(|| T::from_f64(1000000.0).unwrap_or_else(|| T::one())));
        }
        coarse_error / fine_error
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    #[test]
    fn test_l2_norm() {
        let l2 = L2Norm;

        // Test identical arrays
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.0, 2.0, 3.0];
        let error = l2.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        assert_relative_eq!(error, 0.0, epsilon = 1e-10);

        // Test known error
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![0.0, 0.0, 0.0];
        let error = l2.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        let expected = (1.0_f64 + 4.0 + 9.0).sqrt() / 3.0_f64.sqrt(); // sqrt(14/3)
        assert_relative_eq!(error, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_linf_norm() {
        let linf = LInfNorm;

        // Test maximum difference
        let numerical = vec![1.0, 5.0, 3.0];
        let reference = vec![0.0, 2.0, 3.0];
        let error = linf.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        assert_relative_eq!(error, 3.0, epsilon = 1e-10); // max(|1-0|, |5-2|, |3-3|) = 3
    }

    #[test]
    fn test_l1_norm() {
        let l1 = L1Norm;

        // Test mean absolute difference
        let numerical = vec![1.0, 5.0, 3.0];
        let reference = vec![0.0, 2.0, 1.0];
        let error = l1.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        let expected = (1.0 + 3.0 + 2.0) / 3.0; // (|1-0| + |5-2| + |3-1|) / 3 = 2.0
        assert_relative_eq!(error, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_relative_error() {
        let rel_error = RelativeError::new(L2Norm, 1e-12);

        // Test relative error
        let numerical = vec![1.1, 2.1, 3.1];
        let reference = vec![1.0, 2.0, 3.0];
        let error = rel_error.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");

        let abs_error = L2Norm.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        let ref_norm = L2Norm.compute_error(&reference, &vec![0.0; 3]).expect("CRITICAL: Add proper error handling");
        let expected = abs_error / ref_norm;

        assert_relative_eq!(error, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_relative_error_zero_reference() {
        let rel_error = RelativeError::new(L2Norm, 1e-6);

        // Test with near-zero reference (should return absolute error)
        let numerical = vec![1e-8, 2e-8, 3e-8];
        let reference = vec![0.0, 0.0, 0.0];
        let error = rel_error.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");

        let abs_error = L2Norm.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        assert_relative_eq!(error, abs_error, epsilon = 1e-10);
    }

    #[test]
    fn test_rmse() {
        let rmse = RootMeanSquareError;
        let l2 = L2Norm;

        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![0.5, 1.5, 2.5];

        let rmse_error = rmse.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        let l2_error = l2.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");

        // RMSE should be the same as L2 norm
        assert_relative_eq!(rmse_error, l2_error, epsilon = 1e-10);
    }

    #[test]
    fn test_mae() {
        let mae = MeanAbsoluteError;
        let l1 = L1Norm;

        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![0.5, 1.5, 2.5];

        let mae_error = mae.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        let l1_error = l1.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");

        // MAE should be the same as L1 norm
        assert_relative_eq!(mae_error, l1_error, epsilon = 1e-10);
    }

    #[test]
    fn test_normalized_rmse_range() {
        let nrmse = NormalizedRMSE::new(NormalizationMethod::Range);

        let numerical = vec![1.0, 3.0, 5.0];
        let reference = vec![0.0, 2.0, 4.0]; // range = 4.0

        let error = nrmse.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        let rmse = RootMeanSquareError.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        let expected = rmse / 4.0; // normalized by range

        assert_relative_eq!(error, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_normalized_rmse_mean() {
        let nrmse = NormalizedRMSE::new(NormalizationMethod::Mean);

        let numerical = vec![1.0, 3.0, 5.0];
        let reference = vec![0.0, 2.0, 4.0]; // mean = 2.0

        let error = nrmse.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        let rmse = RootMeanSquareError.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        let expected = rmse / 2.0; // normalized by mean

        assert_relative_eq!(error, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_error_statistics() {
        let numerical = vec![1.1, 2.1, 3.1];
        let reference = vec![1.0, 2.0, 3.0];

        let stats = ErrorStatistics::compute(&numerical, &reference).expect("CRITICAL: Add proper error handling");

        assert_eq!(stats.num_points, 3);
        assert!(stats.l1_norm > 0.0);
        assert!(stats.l2_norm > 0.0);
        assert!(stats.linf_norm > 0.0);
        assert!(stats.mae > 0.0);
        assert!(stats.rmse > 0.0);
        assert!(stats.relative_l2 > 0.0);

        // L2 and RMSE should be the same
        assert_relative_eq!(stats.l2_norm, stats.rmse, epsilon = 1e-10);
        // L1 and MAE should be the same
        assert_relative_eq!(stats.l1_norm, stats.mae, epsilon = 1e-10);
    }

    #[test]
    fn test_error_statistics_vector() {
        let numerical = vec![
            Vector3::new(1.1, 0.0, 0.0),
            Vector3::new(2.1, 0.0, 0.0),
            Vector3::new(3.1, 0.0, 0.0),
        ];
        let reference = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
        ];

        let stats = ErrorStatistics::compute_vector(&numerical, &reference).expect("CRITICAL: Add proper error handling");

        assert_eq!(stats.num_points, 3);
        assert!(stats.l2_norm > 0.0);
    }

    #[test]
    fn test_vector_error_metric() {
        let l2 = L2Norm;

        let numerical = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 2.0, 0.0),
        ];
        let reference = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
        ];

        let error = l2.compute_vector_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");

        // Magnitudes are [1.0, 2.0], reference magnitudes are [0.0, 0.0]
        // L2 error = sqrt((1^2 + 2^2) / 2) = sqrt(2.5)
        let expected = (2.5_f64).sqrt();
        assert_relative_eq!(error, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_convergence_rate() {
        // Test theoretical second-order convergence
        let grid_sizes = vec![0.1, 0.05, 0.025];
        let errors = vec![0.01, 0.0025, 0.000625]; // errors ∝ h^2

        let rate = ErrorAnalysis::convergence_rate(&grid_sizes, &errors).expect("CRITICAL: Add proper error handling");
        assert_relative_eq!(rate, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_convergence_rate_first_order() {
        // Test theoretical first-order convergence
        let grid_sizes = vec![0.1, 0.05, 0.025];
        let errors = vec![0.1, 0.05, 0.025]; // errors ∝ h^1

        let rate = ErrorAnalysis::convergence_rate(&grid_sizes, &errors).expect("CRITICAL: Add proper error handling");
        assert_relative_eq!(rate, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_error_reduction_factor() {
        let coarse_error = 0.01;
        let fine_error = 0.0025;

        let factor = ErrorAnalysis::error_reduction_factor(coarse_error, fine_error);
        assert_relative_eq!(factor, 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_is_acceptable() {
        assert!(ErrorAnalysis::is_acceptable(0.001, 0.01));
        assert!(!ErrorAnalysis::is_acceptable(0.01, 0.001));
        assert!(ErrorAnalysis::is_acceptable(0.01, 0.01));
    }

    #[test]
    fn test_empty_arrays() {
        let l2 = L2Norm;
        let numerical: Vec<f64> = vec![];
        let reference: Vec<f64> = vec![];

        let error = l2.compute_error(&numerical, &reference).expect("CRITICAL: Add proper error handling");
        assert_eq!(error, 0.0);
    }

    #[test]
    fn test_mismatched_lengths() {
        let l2 = L2Norm;
        let numerical = vec![1.0, 2.0];
        let reference = vec![1.0, 2.0, 3.0];

        assert!(l2.compute_error(&numerical, &reference).is_err());
    }

    #[test]
    fn test_convergence_rate_insufficient_data() {
        let grid_sizes = vec![0.1];
        let errors = vec![0.01];

        assert!(ErrorAnalysis::convergence_rate(&grid_sizes, &errors).is_err());
    }

    #[test]
    fn test_convergence_rate_identical_grid_sizes() {
        let grid_sizes = vec![0.1, 0.1, 0.1];
        let errors = vec![0.01, 0.01, 0.01];

        assert!(ErrorAnalysis::convergence_rate(&grid_sizes, &errors).is_err());
    }
}
