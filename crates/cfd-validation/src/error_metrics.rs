//! Error metrics for CFD validation.
//!
//! This module provides various error metrics to quantify the difference
//! between numerical solutions and analytical or reference solutions.

use cfd_core::error::{Error, Result};
use nalgebra::{RealField, Vector3};
use num_traits::cast::FromPrimitive;

/// Trait for error metrics
pub trait ErrorMetric<T: RealField + Copy> {
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
        let num_magnitudes: Vec<T> = numerical.iter().map(nalgebra::Matrix::norm).collect();
        let ref_magnitudes: Vec<T> = reference.iter().map(nalgebra::Matrix::norm).collect();

        self.compute_error(&num_magnitudes, &ref_magnitudes)
    }

    /// Get the name of the error metric
    fn name(&self) -> &str;
}

/// L2 (Euclidean) norm error metric
pub struct L2Norm;

impl<T: RealField + Copy + FromPrimitive + Copy> ErrorMetric<T> for L2Norm {
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
                let diff = *num - *ref_val;
                diff * diff
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

impl<T: RealField + Copy> ErrorMetric<T> for LInfNorm {
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
            .map(|(num, ref_val)| (*num - *ref_val).abs())
            .fold(T::zero(), |acc, x| if x > acc { x } else { acc });

        Ok(max_diff)
    }

    fn name(&self) -> &str {
        "L∞ Norm"
    }
}

/// L1 (Manhattan) norm error metric
pub struct L1Norm;

impl<T: RealField + Copy + FromPrimitive + Copy> ErrorMetric<T> for L1Norm {
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
            .map(|(num, ref_val)| (*num - *ref_val).abs())
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
    T: RealField + FromPrimitive + Copy,
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

impl<T: RealField + Copy + FromPrimitive + Copy> ErrorMetric<T> for RootMeanSquareError {
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

impl<T: RealField + Copy + FromPrimitive + Copy> ErrorMetric<T> for MeanAbsoluteError {
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
    #[must_use] pub fn new(method: NormalizationMethod) -> Self {
        Self {
            normalization_method: method,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> ErrorMetric<T> for NormalizedRMSE {
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
                let min_val = reference.iter().fold(reference[0], |acc, x| {
                    if *x < acc { *x } else { acc }
                });
                let max_val = reference.iter().fold(reference[0], |acc, x| {
                    if *x > acc { *x } else { acc }
                });
                max_val - min_val
            },
            NormalizationMethod::Mean => {
                let sum: T = reference.iter().fold(T::zero(), |acc, x| acc + *x);
                let n = T::from_usize(reference.len()).ok_or_else(|| {
                    Error::InvalidConfiguration("Failed to convert array length to target type".to_string())
                })?;
                sum / n
            },
            NormalizationMethod::StandardDeviation => {
                // Compute mean
                let sum: T = reference.iter().fold(T::zero(), |acc, x| acc + *x);
                let n = T::from_usize(reference.len()).ok_or_else(|| {
                    Error::InvalidConfiguration("Failed to convert array length to target type".to_string())
                })?;
                let mean = sum / n;

                // Compute variance
                let variance: T = reference.iter()
                    .map(|x| {
                        let diff = *x - mean;
                        diff * diff
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

/// Utility functions for error analysis
pub struct ErrorAnalysis;

impl ErrorAnalysis {
    /// Compute convergence rate from error measurements at different grid sizes
    pub fn convergence_rate<T: RealField + Copy + FromPrimitive + Copy>(
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
        let n = T::from_usize(grid_sizes.len()).ok_or_else(|| {
            Error::InvalidConfiguration("Failed to convert grid size length to target type".to_string())
        })?;

        let log_h: Vec<T> = grid_sizes.iter().map(|h| h.ln()).collect();
        let log_e: Vec<T> = errors.iter().map(|e| e.ln()).collect();

        // Compute means
        let mean_log_h = log_h.iter().fold(T::zero(), |acc, x| acc + *x) / n;
        let mean_log_e = log_e.iter().fold(T::zero(), |acc, x| acc + *x) / n;

        // Compute slope (convergence rate)
        let numerator: T = log_h.iter().zip(log_e.iter())
            .map(|(h, e)| (*h - mean_log_h) * (*e - mean_log_e))
            .fold(T::zero(), |acc, x| acc + x);

        let denominator: T = log_h.iter()
            .map(|h| {
                let diff = *h - mean_log_h;
                diff * diff
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
    pub fn is_acceptable<T: RealField + Copy>(error: T, tolerance: T) -> bool {
        error <= tolerance
    }

    /// Compute error reduction factor between two measurements
    pub fn error_reduction_factor<T: RealField + Copy>(coarse_error: T, fine_error: T) -> T {
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
    fn test_l2_norm() -> Result<()> {
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.1, 1.9, 3.2];
        
        let l2 = L2Norm;
        let error = l2.compute_error(&numerical, &reference)?;
        
        // Expected: sqrt(((0.1)^2 + (0.1)^2 + (0.2)^2) / 3) ≈ 0.1414
        assert_relative_eq!(error, 0.1414213, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn test_l2_norm_vector() -> Result<()> {
        // Test L2 norm with scalar values - Vector3 doesn't implement RealField
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.1, 1.9, 3.2];
        
        let l2 = L2Norm;
        let error = l2.compute_error(&numerical, &reference)?;
        
        // Expected: sqrt((0.1^2 + 0.1^2 + 0.2^2)/3) ≈ 0.1414
        assert_relative_eq!(error, 0.141421, epsilon = 1e-5);
        Ok(())
    }

    #[test]
    fn test_linf_norm() -> Result<()> {
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.1, 1.9, 3.5];
        
        let linf = LInfNorm;
        let error = linf.compute_error(&numerical, &reference)?;
        
        // Expected: max(0.1, 0.1, 0.5) = 0.5
        assert_relative_eq!(error, 0.5, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_l1_norm() -> Result<()> {
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.1, 1.9, 3.2];
        
        let l1 = L1Norm;
        let error = l1.compute_error(&numerical, &reference)?;
        
        // Expected: (0.1 + 0.1 + 0.2) / 3 ≈ 0.1333
        assert_relative_eq!(error, 0.1333333, epsilon = 1e-5);
        Ok(())
    }

    #[test]
    fn test_relative_l2_error() -> Result<()> {
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.1, 2.2, 3.3];
        
        let rel_error = RelativeError::new(L2Norm, 1e-12);
        let error = rel_error.compute_error(&numerical, &reference)?;
        
        let abs_error = L2Norm.compute_error(&numerical, &reference)?;
        let ref_norm = L2Norm.compute_error(&reference, &vec![0.0; 3])?;
        let expected = abs_error / ref_norm;
        
        assert_relative_eq!(error, expected, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_relative_error_normalization() -> Result<()> {
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.1, 2.2, 3.3];
        
        let rel_error = RelativeError::new(L2Norm, 1e-12);
        let error = rel_error.compute_error(&numerical, &reference)?;
        
        let abs_error = L2Norm.compute_error(&numerical, &reference)?;
        let expected = abs_error / 3.3; // Max of reference
        
        assert_relative_eq!(error, expected, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_rmse() -> Result<()> {
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.1, 1.9, 3.2];
        
        let rmse = RootMeanSquareError;
        let rmse_error = rmse.compute_error(&numerical, &reference)?;
        let l2_error = L2Norm.compute_error(&numerical, &reference)?;
        
        // RMSE should equal L2 norm for this case
        assert_relative_eq!(rmse_error, l2_error, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_mae() -> Result<()> {
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.1, 1.9, 3.2];
        
        let mae = MeanAbsoluteError;
        let mae_error = mae.compute_error(&numerical, &reference)?;
        let l1_error = L1Norm.compute_error(&numerical, &reference)?;
        
        // MAE should equal L1 norm
        assert_relative_eq!(mae_error, l1_error, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_nrmse_range() -> Result<()> {
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.1, 1.9, 3.2];
        
        let nrmse = NormalizedRMSE::new(NormalizationMethod::Range);
        let error = nrmse.compute_error(&numerical, &reference)?;
        let rmse = RootMeanSquareError.compute_error(&numerical, &reference)?;
        let expected = rmse / (3.2 - 1.1); // Range of reference
        
        assert_relative_eq!(error, expected, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_nrmse_mean() -> Result<()> {
        let numerical = vec![1.0, 2.0, 3.0];
        let reference = vec![1.1, 1.9, 3.2];
        
        let nrmse = NormalizedRMSE::new(NormalizationMethod::Mean);
        let error = nrmse.compute_error(&numerical, &reference)?;
        let rmse = RootMeanSquareError.compute_error(&numerical, &reference)?;
        let mean = (1.1 + 1.9 + 3.2) / 3.0;
        let expected = rmse / mean;
        
        assert_relative_eq!(error, expected, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_error_statistics() -> Result<()> {
        let numerical = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let reference = vec![1.1, 1.9, 3.2, 3.8, 5.3];
        
        let stats = ErrorStatistics::compute(&numerical, &reference)?;
        
        assert!(stats.num_points == 5);
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
        Ok(())
    }

    #[test]
    fn test_error_statistics_vector() -> Result<()> {
        use nalgebra::Vector3;
        let numerical = vec![
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(3.0, 4.0, 5.0),
        ];
        let reference = vec![
            Vector3::new(1.1, 1.9, 3.2),
            Vector3::new(3.2, 3.8, 5.3),
        ];
        
        let stats = ErrorStatistics::compute_vector(&numerical, &reference)?;
        
        assert!(stats.num_points == 2);
        assert!(stats.l2_norm > 0.0);
        Ok(())
    }

    #[test]
    fn test_vector_error_metric() -> Result<()> {
        let l2 = L2Norm;

        let numerical = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 2.0, 0.0),
        ];
        let reference = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
        ];

        let error = l2.compute_vector_error(&numerical, &reference)?;

        // Magnitudes are [1.0, 2.0], reference magnitudes are [0.0, 0.0]
        // L2 error = sqrt((1^2 + 2^2) / 2) = sqrt(2.5)
        let expected = (2.5_f64).sqrt();
        assert_relative_eq!(error, expected, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_convergence_rate() -> Result<()> {
        // Test theoretical second-order convergence
        let grid_sizes = vec![0.1, 0.05, 0.025];
        let errors = vec![0.01, 0.0025, 0.000625]; // errors ∝ h^2

        let rate = ErrorAnalysis::convergence_rate(&grid_sizes, &errors)?;
        assert_relative_eq!(rate, 2.0, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_convergence_rate_first_order() -> Result<()> {
        // Test theoretical first-order convergence
        let grid_sizes = vec![0.1, 0.05, 0.025];
        let errors = vec![0.1, 0.05, 0.025]; // errors ∝ h^1

        let rate = ErrorAnalysis::convergence_rate(&grid_sizes, &errors)?;
        assert_relative_eq!(rate, 1.0, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_error_reduction_factor() -> Result<()> {
        let coarse_error = 0.01;
        let fine_error = 0.0025;

        let factor = ErrorAnalysis::error_reduction_factor(coarse_error, fine_error);
        assert_relative_eq!(factor, 4.0, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_is_acceptable() -> Result<()> {
        assert!(ErrorAnalysis::is_acceptable(0.001, 0.01));
        assert!(!ErrorAnalysis::is_acceptable(0.01, 0.001));
        assert!(ErrorAnalysis::is_acceptable(0.01, 0.01));
        Ok(())
    }

    #[test]
    fn test_empty_arrays() -> Result<()> {
        let l2 = L2Norm;
        let numerical: Vec<f64> = vec![];
        let reference: Vec<f64> = vec![];

        let error = l2.compute_error(&numerical, &reference)?;
        assert_eq!(error, 0.0);
        Ok(())
    }

    #[test]
    fn test_mismatched_lengths() -> Result<()> {
        let l2 = L2Norm;
        let numerical = vec![1.0, 2.0];
        let reference = vec![1.0, 2.0, 3.0];

        assert!(l2.compute_error(&numerical, &reference).is_err());
        Ok(())
    }

    #[test]
    fn test_convergence_rate_insufficient_data() -> Result<()> {
        let grid_sizes = vec![0.1];
        let errors = vec![0.01];

        assert!(ErrorAnalysis::convergence_rate(&grid_sizes, &errors).is_err());
        Ok(())
    }

    #[test]
    fn test_convergence_rate_identical_grid_sizes() -> Result<()> {
        let grid_sizes = vec![0.1, 0.1, 0.1];
        let errors = vec![0.01, 0.01, 0.01];

        assert!(ErrorAnalysis::convergence_rate(&grid_sizes, &errors).is_err());
        Ok(())
    }
}
