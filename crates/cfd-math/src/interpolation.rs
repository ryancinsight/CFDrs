//! Interpolation methods for CFD simulations.
//!
//! This module provides various interpolation algorithms optimized for CFD applications
//! with support for both regular and irregular grids.

use cfd_core::error::{Error, Result};
use cfd_core::numeric;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use std::cmp::Ordering;
/// Trait for interpolation methods
pub trait Interpolation<T: RealField + Copy>: Send + Sync {
    /// Interpolate at a given point
    fn interpolate(&self, x: T) -> Result<T>;
    /// Interpolate at multiple points using iterator
    fn interpolate_many<I>(&self, points: I) -> impl Iterator<Item = Result<T>>
    where
        I: IntoIterator<Item = T>,
        Self: Clone,
    {
        let interpolator = self;
        points.into_iter().map(move |x| interpolator.interpolate(x))
    }
    /// Get the domain bounds
    fn bounds(&self) -> (T, T);
}
/// Linear interpolation between data points
#[derive(Clone)]
pub struct LinearInterpolation<T: RealField + Copy> {
    /// X coordinates (must be sorted)
    x_data: Vec<T>,
    /// Y values corresponding to x_data
    y_data: Vec<T>,
impl<T: RealField + Copy> LinearInterpolation<T> {
    /// Create new linear interpolation from data points
    }

    pub fn new(x_data: Vec<T>, y_data: Vec<T>) -> Result<Self> {
        if x_data.len() != y_data.len() {
            return Err(Error::InvalidConfiguration(
                "x_data and y_data must have the same length".to_string(),
            ));
        }
        if x_data.len() < 2 {
                "Need at least 2 points for interpolation".to_string(),
        // Verify x_data is sorted
        if !x_data.windows(2).all(|w| w[0] <= w[1]) {
                "x_data must be sorted in ascending order".to_string(),
        Ok(Self { x_data, y_data })
    /// Find the interval containing x using binary search
    fn find_interval(&self, x: &T) -> Option<usize> {
        match self.x_data.binary_search_by(|probe| {
            if probe < x {
                Ordering::Less
            } else if probe > x {
                Ordering::Greater
            } else {
                Ordering::Equal
            }
        }) {
            Ok(idx) => Some(idx),
            Err(idx) => {
                if idx == 0 || idx >= self.x_data.len() {
                    None
                } else {
                    Some(idx - 1)
                }
impl<T: RealField + Copy> Interpolation<T> for LinearInterpolation<T> {
    fn interpolate(&self, x: T) -> Result<T> {
        // Check bounds
        if x < self.x_data[0] || x > self.x_data[self.x_data.len() - 1] {
                "Interpolation point is outside the data range".to_string(),
        // Find interval
        let idx = self.find_interval(&x).ok_or_else(|| {
            Error::InvalidConfiguration("Failed to find interpolation interval".to_string())
        })?;
        // Handle exact match
        if self.x_data[idx] == x {
            return Ok(self.y_data[idx]);
        // Linear interpolation
        let x0 = &self.x_data[idx];
        let x1 = &self.x_data[idx + 1];
        let y0 = &self.y_data[idx];
        let y1 = &self.y_data[idx + 1];
        let t = (x - *x0) / (*x1 - *x0);
        Ok(*y0 + t * (*y1 - *y0))
    fn bounds(&self) -> (T, T) {
        (self.x_data[0], self.x_data[self.x_data.len() - 1])
/// Cubic spline interpolation
///
/// Natural cubic spline with continuous second derivatives.
/// Reference: Burden & Faires, "Numerical Analysis"
    }

pub struct CubicSplineInterpolation<T: RealField + Copy> {
    /// Y values
    #[allow(dead_code)]
    /// Spline coefficients
    coefficients: SplineCoefficients<T>,
struct SplineCoefficients<T> {
    a: Vec<T>, // y values
    b: Vec<T>, // first derivatives
    c: Vec<T>, // second derivatives / 2
    d: Vec<T>, // third derivatives / 6
impl<T: RealField + FromPrimitive + Copy> CubicSplineInterpolation<T> {
    /// Create new cubic spline interpolation
        if x_data.len() < 3 {
                "Need at least 3 points for cubic spline interpolation".to_string(),
        if !x_data.windows(2).all(|w| w[0] < w[1]) {
                "x_data must be sorted in ascending order with no duplicates".to_string(),
        let coefficients = Self::compute_coefficients(&x_data, &y_data)?;
        Ok(Self {
            x_data,
            y_data,
            coefficients,
        })
    /// Compute spline coefficients using Thomas algorithm
    fn compute_coefficients(x_data: &[T], y_data: &[T]) -> Result<SplineCoefficients<T>> {
        let n = x_data.len();
        let mut h = Vec::with_capacity(n - 1);
        let mut alpha = Vec::with_capacity(n - 1);
        // Compute intervals and divided differences
        for i in 0..n - 1 {
            h.push(x_data[i + 1] - x_data[i]);
        for i in 1..n - 1 {
            let term1 = (y_data[i + 1] - y_data[i]) / h[i];
            let term2 = (y_data[i] - y_data[i - 1]) / h[i - 1];
            alpha.push(cfd_core::numeric::from_f64(3.0)? * (term1 - term2));
        // Solve tridiagonal system for c coefficients
        let mut l = vec![T::one(); n];
        let mut mu = vec![T::zero(); n];
        let mut z = vec![T::zero(); n];
            l[i] = cfd_core::numeric::from_f64(2.0)? * (x_data[i + 1] - x_data[i - 1])
                - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i - 1] - h[i - 1] * z[i - 1]) / l[i];
        // Back substitution
        let mut c = vec![T::zero(); n];
        let mut b = vec![T::zero(); n - 1];
        let mut d = vec![T::zero(); n - 1];
        for j in (0..n - 1).rev() {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y_data[j + 1] - y_data[j]) / h[j]
                - h[j] * (c[j + 1] + cfd_core::numeric::from_f64(2.0)? * c[j])
                    / cfd_core::numeric::from_f64(3.0)?;
            d[j] = (c[j + 1] - c[j]) / (cfd_core::numeric::from_f64(3.0)? * h[j]);
        Ok(SplineCoefficients {
            a: y_data.to_vec(),
            b,
            c: c[..n - 1].to_vec(),
            d,
    /// Find the interval containing x
        self.x_data
            .windows(2)
            .position(|w| x >= &w[0] && x <= &w[1])
impl<T: RealField + FromPrimitive + Copy> Interpolation<T> for CubicSplineInterpolation<T> {
        // Evaluate cubic polynomial
        let dx = x - self.x_data[idx];
        let result = self.coefficients.a[idx]
            + dx * (self.coefficients.b[idx]
                + dx * (self.coefficients.c[idx] + dx * self.coefficients.d[idx]));
        Ok(result)
/// Lagrange polynomial interpolation
    }

}

pub struct LagrangeInterpolation<T: RealField + Copy> {
impl<T: RealField + Copy> LagrangeInterpolation<T> {
    /// Create new Lagrange interpolation
        if x_data.is_empty() {
                "Need at least 1 point for interpolation".to_string(),
    /// Compute Lagrange basis polynomial L_i(x)
    fn lagrange_basis(&self, i: usize, x: &T) -> T {
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .fold(T::one(), |acc, (_j, xj)| {
                acc * (*x - *xj) / (self.x_data[i] - *xj)
            })
impl<T: RealField + Copy> Interpolation<T> for LagrangeInterpolation<T> {
        // Use iterator combinators with enumerate for zero-copy optimization
        Ok(self
            .y_data
            .map(|(i, yi)| *yi * self.lagrange_basis(i, &x))
            .reduce(|acc, term| acc + term)
            .unwrap_or_else(|| T::zero()))
    /// Get the bounds of the interpolation domain
        let min = self
            .x_data
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let max = self
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        match (min, max) {
            (Some(min_val), Some(max_val)) => (*min_val, *max_val),
            _ => (T::zero(), T::zero()), // Should never happen due to constructor validation
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    fn test_linear_interpolation() -> Result<()> {
        let x_data = vec![0.0, 1.0, 2.0];
        let y_data = vec![0.0, 1.0, 4.0];
        let interp = LinearInterpolation::new(x_data, y_data)?;
        // Test exact points
        assert_eq!(interp.interpolate(0.0)?, 0.0);
        assert_eq!(interp.interpolate(1.0)?, 1.0);
        assert_eq!(interp.interpolate(2.0)?, 4.0);
        // Test interpolation
        assert_relative_eq!(interp.interpolate(0.5)?, 0.5, epsilon = 1e-10);
        assert_relative_eq!(interp.interpolate(1.5)?, 2.5, epsilon = 1e-10);
        Ok(())
    }

    fn test_cubic_spline_interpolation() -> Result<()> {
        let x_data = vec![0.0, 1.0, 2.0, 3.0];
        let y_data = vec![0.0, 1.0, 4.0, 9.0];
        let spline = CubicSplineInterpolation::new(x_data, y_data)?;
        assert_relative_eq!(spline.interpolate(0.0)?, 0.0, epsilon = 1e-10);
        assert_relative_eq!(spline.interpolate(2.0)?, 4.0, epsilon = 1e-10);
        // Test smoothness (approximate for quadratic)
        assert_relative_eq!(spline.interpolate(1.5)?, 2.25, epsilon = 0.1);
        assert_relative_eq!(spline.interpolate(2.5)?, 6.25, epsilon = 0.1);
    }

    fn test_lagrange_interpolation() -> Result<()> {
        let y_data = vec![1.0, 3.0, 7.0];
        let interp = LagrangeInterpolation::new(x_data, y_data)?;
        assert_relative_eq!(interp.interpolate(0.0)?, 1.0, epsilon = 1e-10);
        assert_relative_eq!(interp.interpolate(1.0)?, 3.0, epsilon = 1e-10);
        assert_relative_eq!(interp.interpolate(2.0)?, 7.0, epsilon = 1e-10);
        // Test interpolation (quadratic: y = x^2 + x + 1)
        assert_relative_eq!(interp.interpolate(0.5)?, 1.75, epsilon = 1e-10);

    }


}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
