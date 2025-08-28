//! Traits for interpolation methods.

use cfd_core::error::Result;
use nalgebra::RealField;

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
