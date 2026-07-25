//! Lagrange polynomial interpolation — thin wrapper around `leto-ops` SSOT implementation.

use super::traits::Interpolation;
use cfd_core::error::Result;
use eunomia::{FloatElement, RealField};
use leto_ops::application::interpolation::Interpolation1D;

/// Barycentric Lagrange polynomial interpolation.
///
/// Delegates to the Atlas-canonical [`leto_ops::application::interpolation::LagrangeInterpolation`].
#[derive(Clone)]
pub struct LagrangeInterpolation<T: RealField + FloatElement + Copy>(
    leto_ops::application::interpolation::LagrangeInterpolation<T>,
);

impl<T: RealField + FloatElement + Copy> LagrangeInterpolation<T> {
    /// Construct from sorted `x_data` and corresponding `y_data` (≥ 2 points required).
    ///
    /// # Errors
    /// Returns `Err` when nodes are not strictly increasing.
    pub fn new(x_data: Vec<T>, y_data: Vec<T>) -> Result<Self> {
        leto_ops::application::interpolation::LagrangeInterpolation::new(x_data, y_data)
            .map(Self)
            .map_err(cfd_core::error::Error::from)
    }
}

impl<T: RealField + FloatElement + Copy> Interpolation<T> for LagrangeInterpolation<T> {
    fn interpolate(&self, x: T) -> Result<T> {
        self.0.interpolate(x).map_err(cfd_core::error::Error::from)
    }

    fn bounds(&self) -> (T, T) {
        self.0.bounds()
    }
}
