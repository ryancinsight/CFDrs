//! Linear interpolation between data points.

use super::traits::Interpolation;
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use std::cmp::Ordering;

/// Linear interpolation between data points
#[derive(Clone)]
pub struct LinearInterpolation<T: RealField + Copy> {
    /// X coordinates (must be sorted)
    x_data: Vec<T>,
    /// Y values corresponding to x_data
    y_data: Vec<T>,
}

impl<T: RealField + Copy> LinearInterpolation<T> {
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

impl<T: RealField + Copy> Interpolation<T> for LinearInterpolation<T> {
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
            return Ok(self.y_data[idx]);
        }

        // Linear interpolation
        let x0 = &self.x_data[idx];
        let x1 = &self.x_data[idx + 1];
        let y0 = &self.y_data[idx];
        let y1 = &self.y_data[idx + 1];

        let t = (x - *x0) / (*x1 - *x0);
        Ok(*y0 + t * (*y1 - *y0))
    }

    fn bounds(&self) -> (T, T) {
        (self.x_data[0], self.x_data[self.x_data.len() - 1])
    }
}
