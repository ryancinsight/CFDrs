//! Numerical integration methods.

use nalgebra::RealField;

/// Trait for quadrature methods
pub trait Quadrature<T: RealField> {
    /// Integrate a function
    fn integrate(&self, f: impl Fn(T) -> T, a: T, b: T) -> T;
}

/// Gauss quadrature
pub struct GaussQuadrature<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> Quadrature<T> for GaussQuadrature<T> {
    fn integrate(&self, _f: impl Fn(T) -> T, _a: T, _b: T) -> T {
        T::zero()
    }
}