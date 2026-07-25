//! Cubic spline interpolation — thin wrapper around `leto-ops` SSOT implementation.

use super::traits::Interpolation;
use cfd_core::error::Result;
use eunomia::{FloatElement, RealField};
use leto_ops::application::interpolation::Interpolation1D;

/// Natural cubic spline interpolant.
///
/// Delegates to the Atlas-canonical [`leto_ops::application::interpolation::CubicSplineInterpolation`].
#[derive(Clone)]
pub struct CubicSplineInterpolation<T: RealField + FloatElement + Copy>(
    leto_ops::application::interpolation::CubicSplineInterpolation<T>,
);

impl<T: RealField + FloatElement + Copy> CubicSplineInterpolation<T> {
    /// Construct from sorted `x_data` and corresponding `y_data` (≥ 3 points required).
    ///
    /// # Errors
    /// Returns `Err` when data has fewer than 3 strictly-increasing nodes.
    pub fn new(x_data: Vec<T>, y_data: Vec<T>) -> Result<Self> {
        leto_ops::application::interpolation::CubicSplineInterpolation::new(x_data, y_data)
            .map(Self)
            .map_err(cfd_core::error::Error::from)
    }
}

impl<T: RealField + FloatElement + Copy> Interpolation<T> for CubicSplineInterpolation<T> {
    fn interpolate(&self, x: T) -> Result<T> {
        self.0.interpolate(x).map_err(cfd_core::error::Error::from)
    }

    fn bounds(&self) -> (T, T) {
        self.0.bounds()
    }
}

