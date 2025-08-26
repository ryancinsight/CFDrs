//! Variable (adaptive) quadrature with error control

use crate::integration::traits::Quadrature;
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Variable quadrature with error control
pub struct VariableQuadrature<Q> {
    base_rule: Q,
    tolerance: f64,
    max_depth: usize,
}

impl<Q> VariableQuadrature<Q> {
    /// Create variable quadrature with given tolerance
    pub fn new(base_rule: Q, tolerance: f64, max_depth: usize) -> Self {
        Self {
            base_rule,
            tolerance,
            max_depth,
        }
    }

    /// Variable integration with recursive subdivision
    pub fn integrate_adaptive<T, F>(&self, f: F, a: T, b: T) -> Result<T>
    where
        T: RealField + From<f64> + FromPrimitive + Copy,
        F: Fn(T) -> T + Copy,
        Q: Quadrature<T>,
    {
        self.integrate_recursive(f, a, b, 0)
    }

    fn integrate_recursive<T, F>(&self, f: F, a: T, b: T, depth: usize) -> Result<T>
    where
        T: RealField + From<f64> + FromPrimitive + Copy,
        F: Fn(T) -> T + Copy,
        Q: Quadrature<T>,
    {
        if depth > self.max_depth {
            return Err(Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded {
                    max: self.max_depth,
                },
            ));
        }

        // Compute integral over whole interval
        let whole = self.base_rule.integrate(f, a, b);

        // Compute integral over two halves
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        let mid = (a + b) / two;
        let left = self.base_rule.integrate(f, a, mid);
        let right = self.base_rule.integrate(f, mid, b);
        let halves = left + right;

        // Estimate error
        let error_estimate = (halves - whole).abs();
        let tolerance_t = T::from_f64(self.tolerance).unwrap_or_else(|| T::zero());

        if error_estimate < tolerance_t {
            // Accept the more accurate estimate from halves
            Ok(halves)
        } else {
            // Recursively refine both halves
            let left_refined = self.integrate_recursive(f, a, mid, depth + 1)?;
            let right_refined = self.integrate_recursive(f, mid, b, depth + 1)?;
            Ok(left_refined + right_refined)
        }
    }
}
