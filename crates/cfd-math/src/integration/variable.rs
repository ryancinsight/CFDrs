//! Variable (adaptive) quadrature with error control

use crate::integration::traits::Quadrature;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};

fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

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
        T: RealField + FloatElement + Copy,
        F: Fn(T) -> T + Copy,
        Q: Quadrature<T>,
    {
        self.integrate_recursive(f, a, b, 0)
    }

    fn integrate_recursive<T, F>(&self, f: F, a: T, b: T, depth: usize) -> Result<T>
    where
        T: RealField + FloatElement + Copy,
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
        let two = from_f64::<T>(2.0);
        let mid = (a + b) / two;
        let left = self.base_rule.integrate(f, a, mid);
        let right = self.base_rule.integrate(f, mid, b);
        let halves = left + right;

        // Estimate error
        let error_estimate = <T as NumericElement>::abs(halves - whole);
        let tolerance_t = from_f64::<T>(self.tolerance);

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integration::quadrature::GaussQuadrature;
    use approx::assert_relative_eq;

    #[test]
    fn adaptive_quadrature_refines_until_tolerance() {
        let base = GaussQuadrature::<f64>::new(2).unwrap();
        let adaptive = VariableQuadrature::new(base, 1e-12, 20);

        let integral = adaptive.integrate_adaptive(|x| x * x, 0.0, 1.0).unwrap();

        assert_relative_eq!(integral, 1.0 / 3.0, epsilon = 1e-12);
    }
}
