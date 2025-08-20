//! Interpolation methods for CFD simulations.
//!
//! This module provides various interpolation algorithms optimized for CFD applications
//! with support for both regular and irregular grids.

use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use std::cmp::Ordering;

/// Trait for interpolation methods
pub trait Interpolation<T: RealField>: Send + Sync {
    /// Interpolate at a given point
    fn interpolate(&self, x: T) -> Result<T>;

    /// Interpolate at multiple points using iterator
    fn interpolate_many<I>(&self, points: I) -> impl Iterator<Item = Result<T>>
    where
        I: IntoIterator<Item = T>,
        Self: Clone,
    {
        let interpolator = self.clone();
        points.into_iter().map(move |x| interpolator.interpolate(x))
    }

    /// Get the domain bounds
    fn bounds(&self) -> (T, T);
}

/// Linear interpolation between data points
#[derive(Clone)]
pub struct LinearInterpolation<T: RealField> {
    /// X coordinates (must be sorted)
    x_data: Vec<T>,
    /// Y values corresponding to x_data
    y_data: Vec<T>,
}

impl<T: RealField> LinearInterpolation<T> {
    /// Create new linear interpolation from data points
    pub fn new(x_data: Vec<T>, y_data: Vec<T>) -> Result<Self> {
        if x_data.len() != y_data.len() {
            return Err(Error::InvalidConfiguration(
                "x_data and y_data must have the same length".to_string(),
            ));
        }
        
        if x_data.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 points for interpolation".to_string(),
            ));
        }

        // Verify x_data is sorted
        if !x_data.windows(2).all(|w| w[0] <= w[1]) {
            return Err(Error::InvalidConfiguration(
                "x_data must be sorted in ascending order".to_string(),
            ));
        }

        Ok(Self { x_data, y_data })
    }

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
            }
        }
    }
}

impl<T: RealField> Interpolation<T> for LinearInterpolation<T> {
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

        // Handle exact match
        if self.x_data[idx] == x {
            return Ok(self.y_data[idx].clone());
        }

        // Linear interpolation
        let x0 = &self.x_data[idx];
        let x1 = &self.x_data[idx + 1];
        let y0 = &self.y_data[idx];
        let y1 = &self.y_data[idx + 1];

        let t = (x - x0.clone()) / (x1.clone() - x0.clone());
        Ok(y0.clone() + t * (y1.clone() - y0.clone()))
    }

    fn bounds(&self) -> (T, T) {
        (
            self.x_data[0].clone(),
            self.x_data[self.x_data.len() - 1].clone(),
        )
    }
}

/// Cubic spline interpolation
///
/// Natural cubic spline with continuous second derivatives.
/// Reference: Burden & Faires, "Numerical Analysis"
#[derive(Clone)]
pub struct CubicSplineInterpolation<T: RealField> {
    /// X coordinates (must be sorted)
    x_data: Vec<T>,
    /// Y values
    #[allow(dead_code)]
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

impl<T: RealField + FromPrimitive> CubicSplineInterpolation<T> {
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
            h.push(x_data[i + 1].clone() - x_data[i].clone());
        }

        for i in 1..n - 1 {
            let term1 = (y_data[i + 1].clone() - y_data[i].clone()) / h[i].clone();
            let term2 = (y_data[i].clone() - y_data[i - 1].clone()) / h[i - 1].clone();
            alpha.push(T::from_f64(3.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? * (term1 - term2));
        }

        // Solve tridiagonal system for c coefficients
        let mut l = vec![T::one(); n];
        let mut mu = vec![T::zero(); n];
        let mut z = vec![T::zero(); n];

        for i in 1..n - 1 {
            l[i] = T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? * (x_data[i + 1].clone() - x_data[i - 1].clone())
                - h[i - 1].clone() * mu[i - 1].clone();
            mu[i] = h[i].clone() / l[i].clone();
            z[i] = (alpha[i - 1].clone() - h[i - 1].clone() * z[i - 1].clone()) / l[i].clone();
        }

        // Back substitution
        let mut c = vec![T::zero(); n];
        let mut b = vec![T::zero(); n - 1];
        let mut d = vec![T::zero(); n - 1];

        for j in (0..n - 1).rev() {
            c[j] = z[j].clone() - mu[j].clone() * c[j + 1].clone();
            b[j] = (y_data[j + 1].clone() - y_data[j].clone()) / h[j].clone()
                - h[j].clone() * (c[j + 1].clone() + T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? * c[j].clone())
                    / T::from_f64(3.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
            d[j] = (c[j + 1].clone() - c[j].clone()) / (T::from_f64(3.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? * h[j].clone());
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

impl<T: RealField + FromPrimitive> Interpolation<T> for CubicSplineInterpolation<T> {
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
        let dx = x - self.x_data[idx].clone();
        let result = self.coefficients.a[idx].clone()
            + dx.clone() * (self.coefficients.b[idx].clone()
                + dx.clone() * (self.coefficients.c[idx].clone()
                    + dx.clone() * self.coefficients.d[idx].clone()));

        Ok(result)
    }

    fn bounds(&self) -> (T, T) {
        (
            self.x_data[0].clone(),
            self.x_data[self.x_data.len() - 1].clone(),
        )
    }
}

/// Lagrange polynomial interpolation
pub struct LagrangeInterpolation<T: RealField> {
    x_data: Vec<T>,
    y_data: Vec<T>,
}

impl<T: RealField> LagrangeInterpolation<T> {
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
                acc * (x.clone() - xj.clone()) / (self.x_data[i].clone() - xj.clone())
            })
    }
}

impl<T: RealField> Interpolation<T> for LagrangeInterpolation<T> {
    fn interpolate(&self, x: T) -> Result<T> {
        // Use iterator combinators with enumerate for zero-copy optimization
        Ok(self.y_data
            .iter()
            .enumerate()
            .map(|(i, yi)| yi.clone() * self.lagrange_basis(i, &x))
            .reduce(|acc, term| acc + term)
            .unwrap_or_else(|| T::zero()))
    }

    fn bounds(&self) -> (T, T) {
        // Use iterator min/max with partial_cmp for better performance
        let (min, max) = self.x_data
            .iter()
            .fold((None, None), |(min_acc, max_acc), x| {
                let new_min = min_acc.map_or(Some(x), |m| if x < m { Some(x) } else { Some(m) });
                let new_max = max_acc.map_or(Some(x), |m| if x > m { Some(x) } else { Some(m) });
                (new_min, new_max)
            });

        (min.expect("x_data is guaranteed to be non-empty by constructor").clone(), max.expect("x_data is guaranteed to be non-empty by constructor").clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_linear_interpolation() {
        let x_data = vec![0.0, 1.0, 2.0, 3.0];
        let y_data = vec![0.0, 1.0, 4.0, 9.0];
        
        let interp = LinearInterpolation::new(x_data, y_data).unwrap();
        
        // Test exact points
        assert_eq!(interp.interpolate(0.0).unwrap(), 0.0);
        assert_eq!(interp.interpolate(1.0).unwrap(), 1.0);
        assert_eq!(interp.interpolate(2.0).unwrap(), 4.0);
        
        // Test interpolation
        assert_relative_eq!(interp.interpolate(0.5).unwrap(), 0.5, epsilon = 1e-10);
        assert_relative_eq!(interp.interpolate(1.5).unwrap(), 2.5, epsilon = 1e-10);
    }

    #[test]
    fn test_cubic_spline() {
        // Test with a known function: y = x^2
        let x_data = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y_data = vec![0.0, 1.0, 4.0, 9.0, 16.0];
        
        let spline = CubicSplineInterpolation::new(x_data, y_data).unwrap();
        
        // Test at data points
        assert_relative_eq!(spline.interpolate(0.0).unwrap(), 0.0, epsilon = 1e-10);
        assert_relative_eq!(spline.interpolate(2.0).unwrap(), 4.0, epsilon = 1e-10);
        
        // Test interpolation (should be close to x^2)
        assert_relative_eq!(spline.interpolate(1.5).unwrap(), 2.25, epsilon = 0.1);
        assert_relative_eq!(spline.interpolate(2.5).unwrap(), 6.25, epsilon = 0.1);
    }

    #[test]
    fn test_lagrange_interpolation() {
        // Test with a quadratic function
        let x_data = vec![0.0, 1.0, 2.0];
        let y_data = vec![1.0, 3.0, 7.0]; // y = x^2 + x + 1
        
        let interp = LagrangeInterpolation::new(x_data, y_data).unwrap();
        
        // Test exact recovery at nodes
        assert_relative_eq!(interp.interpolate(0.0).unwrap(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(interp.interpolate(1.0).unwrap(), 3.0, epsilon = 1e-10);
        assert_relative_eq!(interp.interpolate(2.0).unwrap(), 7.0, epsilon = 1e-10);
        
        // Test at intermediate point
        assert_relative_eq!(interp.interpolate(0.5).unwrap(), 1.75, epsilon = 1e-10);
    }
}