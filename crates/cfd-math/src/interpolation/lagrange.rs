//! Lagrange polynomial interpolation.

use super::traits::Interpolation;
use cfd_core::error::{Error, Result};
use nalgebra::RealField;

/// Lagrange polynomial interpolation
pub struct LagrangeInterpolation<T: RealField + Copy> {
    x_data: Vec<T>,
    y_data: Vec<T>,
}

impl<T: RealField + Copy> LagrangeInterpolation<T> {
    /// Create new Lagrange interpolation
    pub fn new(x_data: Vec<T>, y_data: Vec<T>) -> Result<Self> {
        if x_data.len() != y_data.len() {
            return Err(Error::InvalidConfiguration(
                "x_data and y_data must have the same length".to_string(),
            ));
        }

        if x_data.is_empty() {
            return Err(Error::InvalidConfiguration(
                "Need at least 1 point for interpolation".to_string(),
            ));
        }

        Ok(Self { x_data, y_data })
    }

    /// Compute Lagrange basis polynomial L_i(x)
    fn lagrange_basis(&self, i: usize, x: &T) -> T {
        self.x_data
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .fold(T::one(), |acc, (_j, xj)| {
                acc * (*x - *xj) / (self.x_data[i] - *xj)
            })
    }
}

impl<T: RealField + Copy> Interpolation<T> for LagrangeInterpolation<T> {
    fn interpolate(&self, x: T) -> Result<T> {
        // Use iterator combinators with enumerate for zero-copy optimization
        Ok(self
            .y_data
            .iter()
            .enumerate()
            .map(|(i, yi)| *yi * self.lagrange_basis(i, &x))
            .reduce(|acc, term| acc + term)
            .unwrap_or_else(|| T::zero()))
    }

    /// Get the bounds of the interpolation domain
    fn bounds(&self) -> (T, T) {
        let min = self
            .x_data
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let max = self
            .x_data
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        match (min, max) {
            (Some(min_val), Some(max_val)) => (*min_val, *max_val),
            _ => (T::zero(), T::zero()), // Should never happen due to constructor validation
        }
    }
}
