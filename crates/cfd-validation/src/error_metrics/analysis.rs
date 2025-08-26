//! Error analysis utilities and convergence rate computation

use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
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
                "Need at least 2 grid sizes and corresponding errors".to_string(),
            ));
        }
        // Use least squares fit to log(error) = log(C) + p*log(h)
        // where p is the convergence rate
        let n = T::from_usize(grid_sizes.len()).ok_or_else(|| {
            Error::InvalidConfiguration(
                "Failed to convert grid size length to target type".to_string(),
            )
        })?;
        let log_h: Vec<T> = grid_sizes.iter().map(|h| h.ln()).collect();
        let log_e: Vec<T> = errors.iter().map(|e| e.ln()).collect();
        // Compute means
        let mean_log_h = log_h.iter().fold(T::zero(), |acc, x| acc + *x) / n;
        let mean_log_e = log_e.iter().fold(T::zero(), |acc, x| acc + *x) / n;
        // Compute slope (convergence rate)
        let numerator: T = log_h
            .iter()
            .zip(log_e.iter())
            .map(|(h, e)| (*h - mean_log_h) * (*e - mean_log_e))
            .fold(T::zero(), |acc, x| acc + x);
        let denominator: T = log_h
            .map(|h| {
                let diff = *h - mean_log_h;
                diff * diff
            })
        if denominator == T::zero() {
                "Cannot compute convergence rate: grid sizes are identical".to_string(),
        Ok(numerator / denominator)
    }
    /// Check if error is within acceptable tolerance
    pub fn is_acceptable<T: RealField + Copy>(error: T, tolerance: T) -> bool {
        error <= tolerance
    /// Compute error reduction factor between two measurements
    pub fn error_reduction_factor<T: RealField + Copy>(coarse_error: T, fine_error: T) -> T {
        if fine_error == T::zero() {
            return T::max_value().unwrap_or_else(|| {
                T::from_f64(1e10)
                    .unwrap_or_else(|| T::from_f64(1000000.0).unwrap_or_else(|| T::one()))
            });
        coarse_error / fine_error
}
