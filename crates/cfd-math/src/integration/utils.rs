//! Utility functions for common integration patterns

use super::{
    CompositeQuadrature, GaussQuadrature, Quadrature, SimpsonsRule, TrapezoidalRule,
    VariableQuadrature,
};
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Utility functions for common integration patterns
pub struct IntegrationUtils;

impl IntegrationUtils {
    /// Integrate using trapezoidal rule with n intervals
    pub fn trapezoidal<T, F>(f: F, a: T, b: T, n: usize) -> T
    where
        T: RealField + From<f64> + FromPrimitive + Copy,
        F: Fn(T) -> T,
    {
        let composite = CompositeQuadrature::new(TrapezoidalRule, n);
        composite.integrate(f, a, b)
    }

    /// Integrate using Simpson's rule with n intervals (must be even)
    pub fn simpsons<T, F>(f: F, a: T, b: T, n: usize) -> Result<T>
    where
        T: RealField + From<f64> + FromPrimitive + Copy,
        F: Fn(T) -> T,
    {
        if n % 2 != 0 {
            return Err(Error::InvalidConfiguration(
                "Simpson's rule requires even number of intervals".to_string(),
            ));
        }

        let composite = CompositeQuadrature::new(SimpsonsRule, n / 2);
        Ok(composite.integrate(f, a, b))
    }

    /// Integrate using Gauss-Legendre quadrature
    pub fn gauss_legendre<T, F>(f: F, a: T, b: T, order: usize) -> Result<T>
    where
        T: RealField + From<f64> + FromPrimitive + Copy,
        F: Fn(T) -> T,
    {
        let gauss = GaussQuadrature::new(order)?;
        Ok(gauss.integrate(f, a, b))
    }

    /// Variable integration with automatic error control
    pub fn variable<T, F>(f: F, a: T, b: T, tolerance: f64) -> Result<T>
    where
        T: RealField + From<f64> + FromPrimitive + Copy,
        F: Fn(T) -> T + Copy,
    {
        let gauss = GaussQuadrature::new(3)?; // Use 3-point Gauss rule
        let variable = VariableQuadrature::new(gauss, tolerance, 20);
        variable.integrate_adaptive(f, a, b)
    }
}
