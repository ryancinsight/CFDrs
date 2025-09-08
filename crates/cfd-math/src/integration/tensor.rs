//! Multi-dimensional integration using tensor products

use crate::integration::traits::Quadrature;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Multi-dimensional integration using tensor products
pub struct TensorProductQuadrature<Q> {
    base_rule: Q,
    dimension: usize,
}

impl<Q> TensorProductQuadrature<Q> {
    /// Create tensor product quadrature
    pub fn new(base_rule: Q, dimension: usize) -> Self {
        Self {
            base_rule,
            dimension,
        }
    }

    /// Get the dimension of the quadrature
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Integrate over 2D rectangle [ax, bx] × [ay, by]
    pub fn integrate_2d<T, F>(&self, f: F, ax: T, bx: T, ay: T, by: T) -> T
    where
        T: RealField + From<f64> + FromPrimitive + Copy + Clone,
        F: Fn(T, T) -> T,
        Q: Quadrature<T>,
    {
        // Use the base quadrature rule for each dimension
        // Integrate first over y, then over x using the base rule
        let integral_over_y = |x: T| -> T { self.base_rule.integrate(|y| f(x, y), ay, by) };

        self.base_rule.integrate(integral_over_y, ax, bx)
    }
}
