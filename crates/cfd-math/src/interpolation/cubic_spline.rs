//! Cubic spline interpolation.

use super::traits::Interpolation;
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

// Named constants for spline calculations
const TWO: f64 = 2.0;
const THREE: f64 = 3.0;

/// Cubic spline interpolation
///
/// Natural cubic spline with continuous second derivatives.
/// Reference: Burden & Faires, "Numerical Analysis"
#[derive(Clone)]
pub struct CubicSplineInterpolation<T: RealField + Copy> {
    /// X coordinates (must be sorted)
    x_data: Vec<T>,
    /// Y values
    y_data: Vec<T>,
    /// Spline coefficients
    coefficients: SplineCoefficients<T>,
}

#[derive(Clone)]
struct SplineCoefficients<T> {
    a: Vec<T>, // y values
    b: Vec<T>, // first derivatives
    c: Vec<T>, // second derivatives / 2
    d: Vec<T>, // third derivatives / 6
}

impl<T: RealField + FromPrimitive + Copy> CubicSplineInterpolation<T> {
    /// Create new cubic spline interpolation
    pub fn new(x_data: Vec<T>, y_data: Vec<T>) -> Result<Self> {
        if x_data.len() != y_data.len() {
            return Err(Error::InvalidConfiguration(
                "x_data and y_data must have the same length".to_string(),
            ));
        }

        if x_data.len() < 3 {
            return Err(Error::InvalidConfiguration(
                "Need at least 3 points for cubic spline interpolation".to_string(),
            ));
        }

        // Verify x_data is sorted
        if !x_data.windows(2).all(|w| w[0] < w[1]) {
            return Err(Error::InvalidConfiguration(
                "x_data must be sorted in ascending order with no duplicates".to_string(),
            ));
        }

        let coefficients = Self::compute_coefficients(&x_data, &y_data)?;

        Ok(Self {
            x_data,
            y_data,
            coefficients,
        })
    }

    /// Compute spline coefficients using Thomas algorithm
    fn compute_coefficients(x_data: &[T], y_data: &[T]) -> Result<SplineCoefficients<T>> {
        let n = x_data.len();
        let mut h = Vec::with_capacity(n - 1);
        let mut alpha = Vec::with_capacity(n - 1);

        // Compute intervals and divided differences
        for i in 0..n - 1 {
            h.push(x_data[i + 1] - x_data[i]);
        }

        for i in 1..n - 1 {
            let term1 = (y_data[i + 1] - y_data[i]) / h[i];
            let term2 = (y_data[i] - y_data[i - 1]) / h[i - 1];
            alpha.push(T::from_f64(THREE).unwrap_or_else(|| T::zero()) * (term1 - term2));
        }

        // Solve tridiagonal system for c coefficients
        let mut l = vec![T::one(); n];
        let mut mu = vec![T::zero(); n];
        let mut z = vec![T::zero(); n];

        for i in 1..n - 1 {
            l[i] = T::from_f64(TWO).unwrap_or_else(|| T::zero()) * (x_data[i + 1] - x_data[i - 1])
                - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i - 1] - h[i - 1] * z[i - 1]) / l[i];
        }

        // Back substitution
        let mut c = vec![T::zero(); n];
        let mut b = vec![T::zero(); n - 1];
        let mut d = vec![T::zero(); n - 1];

        for j in (0..n - 1).rev() {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y_data[j + 1] - y_data[j]) / h[j]
                - h[j] * (c[j + 1] + T::from_f64(TWO).unwrap_or_else(|| T::zero()) * c[j])
                    / T::from_f64(THREE).unwrap_or_else(|| T::zero());
            d[j] = (c[j + 1] - c[j]) / (T::from_f64(THREE).unwrap_or_else(|| T::zero()) * h[j]);
        }

        Ok(SplineCoefficients {
            a: y_data.to_vec(),
            b,
            c: c[..n - 1].to_vec(),
            d,
        })
    }

    /// Find the interval containing x
    fn find_interval(&self, x: &T) -> Option<usize> {
        self.x_data
            .windows(2)
            .position(|w| x >= &w[0] && x <= &w[1])
    }
}

impl<T: RealField + FromPrimitive + Copy> Interpolation<T> for CubicSplineInterpolation<T> {
    fn interpolate(&self, x: T) -> Result<T> {
        // Check bounds
        if x < self.x_data[0] || x > self.x_data[self.x_data.len() - 1] {
            return Err(Error::InvalidConfiguration(
                "Interpolation point is outside the data range".to_string(),
            ));
        }

        // Find interval
        let idx = self.find_interval(&x).ok_or_else(|| {
            Error::InvalidConfiguration("Failed to find interpolation interval".to_string())
        })?;

        // Evaluate cubic polynomial
        let dx = x - self.x_data[idx];
        let result = self.coefficients.a[idx]
            + dx * (self.coefficients.b[idx]
                + dx * (self.coefficients.c[idx] + dx * self.coefficients.d[idx]));

        Ok(result)
    }

    fn bounds(&self) -> (T, T) {
        (self.x_data[0], self.x_data[self.x_data.len() - 1])
    }
}
