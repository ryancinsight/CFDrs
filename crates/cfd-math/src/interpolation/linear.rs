//! Linear interpolation — thin wrapper around `leto-ops` SSOT implementation.
//!
//! All algorithmic logic lives in `leto_ops::application::interpolation`.
//! This module re-exports the canonical type and adapts its error to
//! `cfd_core::Error` via the `From<LetoError>` bridge.

use super::traits::Interpolation;
use cfd_core::error::Result;
use eunomia::{FloatElement, RealField};
use leto_ops::application::interpolation::Interpolation1D;

/// Linear interpolation between data points.
///
/// Delegates to the Atlas-canonical [`leto_ops::application::interpolation::LinearInterpolation`].
#[derive(Clone)]
pub struct LinearInterpolation<T: RealField + FloatElement + Copy>(
    leto_ops::application::interpolation::LinearInterpolation<T>,
);

impl<T: RealField + FloatElement + Copy> LinearInterpolation<T> {
    /// Construct from sorted `x_data` and corresponding `y_data`.
    ///
    /// # Errors
    /// Returns `Err` when data is empty, has fewer than 2 points, or nodes are not strictly increasing.
    pub fn new(x_data: Vec<T>, y_data: Vec<T>) -> Result<Self> {
        leto_ops::application::interpolation::LinearInterpolation::new(x_data, y_data)
            .map(Self)
            .map_err(cfd_core::error::Error::from)
    }
}

impl<T: RealField + FloatElement + Copy> Interpolation<T> for LinearInterpolation<T> {
    fn interpolate(&self, x: T) -> Result<T> {
        self.0.interpolate(x).map_err(cfd_core::error::Error::from)
    }

    fn bounds(&self) -> (T, T) {
        self.0.bounds()
    }
}
