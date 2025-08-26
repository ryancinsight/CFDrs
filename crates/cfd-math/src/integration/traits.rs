//! Core traits for numerical integration

use nalgebra::RealField;

/// Trait for quadrature methods
pub trait Quadrature<T: RealField + Copy> {
    /// Integrate a function over interval [a, b]
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T;

    /// Get the order of accuracy
    fn order(&self) -> usize;

    /// Get the number of quadrature points
    fn num_points(&self) -> usize;
}
