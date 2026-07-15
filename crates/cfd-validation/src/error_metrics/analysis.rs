//! Error analysis utilities and convergence rate computation

use crate::scalar;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, RealField};

/// Utility functions for error analysis
pub struct ErrorAnalysis;

impl ErrorAnalysis {
    /// Compute convergence rate from error measurements at different grid sizes
    pub fn convergence_rate<T: RealField + Copy + FloatElement>(
        grid_sizes: &[T],
        errors: &[T],
    ) -> Result<T> {
        if grid_sizes.len() != errors.len() || grid_sizes.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 grid sizes and corresponding errors".to_string(),
            ));
        }

        // Use least squares fit to log(error) = log(C) + p*log(h)
        // where p is the convergence rate
        let n = scalar::from_usize::<T>(grid_sizes.len());

        let log_h: Vec<T> = grid_sizes.iter().map(|h| scalar::ln(*h)).collect();
        let log_e: Vec<T> = errors.iter().map(|e| scalar::ln(*e)).collect();

        // Compute means
        let mean_log_h = log_h.iter().fold(scalar::zero::<T>(), |acc, x| acc + *x) / n;
        let mean_log_e = log_e.iter().fold(scalar::zero::<T>(), |acc, x| acc + *x) / n;

        // Compute slope (convergence rate)
        let numerator: T = log_h
            .iter()
            .zip(log_e.iter())
            .map(|(h, e)| (*h - mean_log_h) * (*e - mean_log_e))
            .fold(scalar::zero::<T>(), |acc, x| acc + x);

        let denominator: T = log_h
            .iter()
            .map(|h| {
                let diff = *h - mean_log_h;
                diff * diff
            })
            .fold(scalar::zero::<T>(), |acc, x| acc + x);

        if denominator == scalar::zero::<T>() {
            return Err(Error::InvalidConfiguration(
                "Cannot compute convergence rate: grid sizes are identical".to_string(),
            ));
        }

        Ok(numerator / denominator)
    }

    /// Check if error is within acceptable tolerance
    pub fn is_acceptable<T: RealField + Copy>(error: T, tolerance: T) -> bool {
        error <= tolerance
    }

    /// Compute error reduction factor between two measurements
    pub fn error_reduction_factor<T: RealField + Copy + FloatElement>(
        coarse_error: T,
        fine_error: T,
    ) -> T {
        if fine_error == scalar::zero::<T>() {
            return <T as RealField>::max_value();
        }
        coarse_error / fine_error
    }
}
