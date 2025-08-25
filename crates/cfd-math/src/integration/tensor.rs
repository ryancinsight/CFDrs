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
        // Use the base quadrature rule for tensor product integration
        // For now, use composite Simpson's rule for 2D
        let n = 10; // Number of intervals in each direction
        let hx = (bx - ax) / T::from_usize(n).unwrap_or_else(|| T::zero());
        let hy = (by - ay) / T::from_usize(n).unwrap_or_else(|| T::zero());

        let mut result = T::zero();
        let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        let nine = T::from_f64(9.0).unwrap_or_else(|| T::zero());

        for i in 0..=n {
            for j in 0..=n {
                let x = ax + T::from_usize(i).unwrap_or_else(|| T::zero()) * hx;
                let y = ay + T::from_usize(j).unwrap_or_else(|| T::zero()) * hy;

                // Correct 2D Simpson's rule weights using tensor product of 1D weights
                let weight_i = if i == 0 || i == n {
                    T::one()
                } else if i % 2 == 1 {
                    four
                } else {
                    two
                };

                let weight_j = if j == 0 || j == n {
                    T::one()
                } else if j % 2 == 1 {
                    four
                } else {
                    two
                };

                let weight = weight_i * weight_j;

                result += weight * f(x, y);
            }
        }

        result * hx * hy / nine
    }
}
