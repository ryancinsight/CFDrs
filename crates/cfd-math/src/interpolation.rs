//! Interpolation methods for CFD simulations.

use nalgebra::RealField;

/// Trait for interpolation methods
pub trait Interpolation<T: RealField> {
    /// Interpolate at a given point
    fn interpolate(&self, x: T) -> T;
}

/// Linear interpolation
pub struct LinearInterpolation<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> Interpolation<T> for LinearInterpolation<T> {
    fn interpolate(&self, _x: T) -> T {
        T::zero()
    }
}

/// Cubic spline interpolation
pub struct CubicSplineInterpolation<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> Interpolation<T> for CubicSplineInterpolation<T> {
    fn interpolate(&self, _x: T) -> T {
        T::zero()
    }
}