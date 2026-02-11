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

/// Trait for 3D quadrature methods (specifically for tetrahedra)
pub trait Quadrature3D<T: RealField + Copy> {
    /// Get quadrature points in barycentric coordinates (L1, L2, L3, L4)
    fn points(&self) -> &[[T; 4]];

    /// Get quadrature weights (sum should be 1.0)
    fn weights(&self) -> &[T];

    /// Get the order (degree) of accuracy
    fn order(&self) -> usize;

    /// Get the number of quadrature points
    fn num_points(&self) -> usize {
        self.weights().len()
    }
}
