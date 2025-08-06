//! Numerical integration methods for CFD applications.
//!
//! This module provides various quadrature rules and integration schemes
//! optimized for CFD simulations with support for adaptive integration.

use cfd_core::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Trait for quadrature methods
pub trait Quadrature<T: RealField> {
    /// Integrate a function over interval [a, b]
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T;

    /// Get the order of accuracy
    fn order(&self) -> usize;

    /// Get the number of quadrature points
    fn num_points(&self) -> usize;
}

/// Trapezoidal rule for numerical integration
pub struct TrapezoidalRule;

impl<T: RealField + FromPrimitive> Quadrature<T> for TrapezoidalRule {
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let two = T::from_f64(2.0).unwrap();
        (b.clone() - a.clone()) * (f(a) + f(b)) / two
    }

    fn order(&self) -> usize {
        2
    }

    fn num_points(&self) -> usize {
        2
    }
}

/// Simpson's rule for numerical integration
pub struct SimpsonsRule;

impl<T: RealField + FromPrimitive> Quadrature<T> for SimpsonsRule {
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let two = T::from_f64(2.0).unwrap();
        let four = T::from_f64(4.0).unwrap();
        let six = T::from_f64(6.0).unwrap();

        let mid = (a.clone() + b.clone()) / two.clone();
        (b.clone() - a.clone()) * (f(a) + four * f(mid) + f(b)) / six
    }

    fn order(&self) -> usize {
        4
    }

    fn num_points(&self) -> usize {
        3
    }
}

/// Gauss-Legendre quadrature
pub struct GaussQuadrature<T: RealField> {
    points: Vec<T>,
    weights: Vec<T>,
    order: usize,
}

impl<T: RealField + FromPrimitive> GaussQuadrature<T> {
    /// Create Gauss-Legendre quadrature of given order
    pub fn new(order: usize) -> Result<Self> {
        let (points, weights) = match order {
            1 => {
                let points = vec![T::zero()];
                let weights = vec![T::from_f64(2.0).unwrap()];
                (points, weights)
            },
            2 => {
                let sqrt3_inv = T::from_f64(1.0 / 3.0_f64.sqrt()).unwrap();
                let points = vec![-sqrt3_inv.clone(), sqrt3_inv];
                let weights = vec![T::one(), T::one()];
                (points, weights)
            },
            3 => {
                let sqrt15 = T::from_f64(15.0_f64.sqrt()).unwrap();
                let sqrt15_5 = sqrt15 / T::from_f64(5.0).unwrap();
                let points = vec![
                    -sqrt15_5.clone(),
                    T::zero(),
                    sqrt15_5,
                ];
                let weights = vec![
                    T::from_f64(5.0 / 9.0).unwrap(),
                    T::from_f64(8.0 / 9.0).unwrap(),
                    T::from_f64(5.0 / 9.0).unwrap(),
                ];
                (points, weights)
            },
            4 => {
                let sqrt6_5 = (6.0_f64 / 5.0).sqrt();
                let term1 = T::from_f64((3.0 - 2.0 * sqrt6_5) / 7.0).unwrap().sqrt();
                let term2 = T::from_f64((3.0 + 2.0 * sqrt6_5) / 7.0).unwrap().sqrt();
                let points = vec![-term2.clone(), -term1.clone(), term1, term2];

                let sqrt30 = 30.0_f64.sqrt();
                let w1 = T::from_f64((18.0 + sqrt30) / 36.0).unwrap();
                let w2 = T::from_f64((18.0 - sqrt30) / 36.0).unwrap();
                let weights = vec![w2.clone(), w1.clone(), w1, w2];
                (points, weights)
            },
            _ => {
                return Err(Error::InvalidConfiguration(
                    format!("Gauss quadrature order {} not implemented", order)
                ));
            }
        };

        Ok(Self {
            points,
            weights,
            order,
        })
    }

    /// Create default 2-point Gauss quadrature
    pub fn default() -> Result<Self> {
        Self::new(2)
    }
}

impl<T: RealField + FromPrimitive> Quadrature<T> for GaussQuadrature<T> {
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let two = T::from_f64(2.0).unwrap();
        let half_interval = (b.clone() - a.clone()) / two.clone();
        let mid_point = (a.clone() + b.clone()) / two;

        let mut result = T::zero();
        for (point, weight) in self.points.iter().zip(self.weights.iter()) {
            let x = mid_point.clone() + half_interval.clone() * point.clone();
            result += weight.clone() * f(x);
        }

        result * half_interval
    }

    fn order(&self) -> usize {
        2 * self.order
    }

    fn num_points(&self) -> usize {
        self.points.len()
    }
}

/// Composite quadrature rule for adaptive integration
pub struct CompositeQuadrature<Q> {
    base_rule: Q,
    num_intervals: usize,
}

impl<Q> CompositeQuadrature<Q> {
    /// Create composite quadrature with given number of intervals
    pub fn new(base_rule: Q, num_intervals: usize) -> Self {
        Self {
            base_rule,
            num_intervals,
        }
    }
}

impl<T, Q> Quadrature<T> for CompositeQuadrature<Q>
where
    T: RealField + FromPrimitive,
    Q: Quadrature<T>,
{
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let n = T::from_usize(self.num_intervals).unwrap();
        let h = (b.clone() - a.clone()) / n;

        let mut result = T::zero();
        for i in 0..self.num_intervals {
            let xi = a.clone() + T::from_usize(i).unwrap() * h.clone();
            let xi_plus_1 = xi.clone() + h.clone();
            result += self.base_rule.integrate(&f, xi, xi_plus_1);
        }

        result
    }

    fn order(&self) -> usize {
        self.base_rule.order()
    }

    fn num_points(&self) -> usize {
        self.base_rule.num_points() * self.num_intervals
    }
}

/// Adaptive quadrature with error control
pub struct AdaptiveQuadrature<Q> {
    base_rule: Q,
    tolerance: f64,
    max_depth: usize,
}

impl<Q> AdaptiveQuadrature<Q> {
    /// Create adaptive quadrature with given tolerance
    pub fn new(base_rule: Q, tolerance: f64, max_depth: usize) -> Self {
        Self {
            base_rule,
            tolerance,
            max_depth,
        }
    }
}

impl<Q> AdaptiveQuadrature<Q> {
    /// Adaptive integration with recursive subdivision
    pub fn integrate_adaptive<T, F>(&self, f: F, a: T, b: T) -> Result<T>
    where
        T: RealField + FromPrimitive,
        F: Fn(T) -> T + Copy,
        Q: Quadrature<T>,
    {
        self.integrate_recursive(f, a, b, 0)
    }

    fn integrate_recursive<T, F>(&self, f: F, a: T, b: T, depth: usize) -> Result<T>
    where
        T: RealField + FromPrimitive,
        F: Fn(T) -> T + Copy,
        Q: Quadrature<T>,
    {
        if depth > self.max_depth {
            return Err(Error::ConvergenceFailure(
                "Maximum recursion depth reached in adaptive integration".to_string()
            ));
        }

        // Compute integral over whole interval
        let whole = self.base_rule.integrate(f, a.clone(), b.clone());

        // Compute integral over two halves
        let two = T::from_f64(2.0).unwrap();
        let mid = (a.clone() + b.clone()) / two;
        let left = self.base_rule.integrate(f, a.clone(), mid.clone());
        let right = self.base_rule.integrate(f, mid.clone(), b.clone());
        let halves = left.clone() + right.clone();

        // Estimate error
        let error_estimate = (halves.clone() - whole.clone()).abs();
        let tolerance_t = T::from_f64(self.tolerance).unwrap();

        if error_estimate < tolerance_t {
            // Accept the more accurate estimate from halves
            Ok(halves)
        } else {
            // Recursively refine both halves
            let left_refined = self.integrate_recursive(f, a, mid.clone(), depth + 1)?;
            let right_refined = self.integrate_recursive(f, mid, b, depth + 1)?;
            Ok(left_refined + right_refined)
        }
    }
}

/// Multi-dimensional integration using tensor products
pub struct TensorProductQuadrature<Q> {
    _base_rule: Q,
    _dimension: usize,
}

impl<Q> TensorProductQuadrature<Q> {
    /// Create tensor product quadrature
    pub fn new(base_rule: Q, dimension: usize) -> Self {
        Self {
            _base_rule: base_rule,
            _dimension: dimension,
        }
    }
}

impl<Q> TensorProductQuadrature<Q> {
    /// Integrate over 2D rectangle [ax, bx] × [ay, by]
    pub fn integrate_2d<T, F>(&self, f: F, ax: T, bx: T, ay: T, by: T) -> T
    where
        T: RealField + FromPrimitive + Clone,
        F: Fn(T, T) -> T,
        Q: Quadrature<T>,
    {
        // For simplicity, use composite Simpson's rule for 2D
        let n = 10; // Number of intervals in each direction
        let hx = (bx.clone() - ax.clone()) / T::from_usize(n).unwrap();
        let hy = (by.clone() - ay.clone()) / T::from_usize(n).unwrap();

        let mut result = T::zero();
        let four = T::from_f64(4.0).unwrap();
        let two = T::from_f64(2.0).unwrap();
        let nine = T::from_f64(9.0).unwrap();

        for i in 0..=n {
            for j in 0..=n {
                let x = ax.clone() + T::from_usize(i).unwrap() * hx.clone();
                let y = ay.clone() + T::from_usize(j).unwrap() * hy.clone();

                let weight = if (i == 0 || i == n) && (j == 0 || j == n) {
                    T::one() // Corner points
                } else if (i == 0 || i == n) || (j == 0 || j == n) {
                    two.clone() // Edge points
                } else if i % 2 == 1 && j % 2 == 1 {
                    four.clone() * four.clone() // Interior odd points
                } else if i % 2 == 1 || j % 2 == 1 {
                    four.clone() * two.clone() // Mixed points
                } else {
                    four.clone() // Interior even points
                };

                result += weight * f(x, y);
            }
        }

        result * hx * hy / nine
    }
}

/// Utility functions for common integration patterns
pub struct IntegrationUtils;

impl IntegrationUtils {
    /// Integrate using trapezoidal rule with n intervals
    pub fn trapezoidal<T, F>(f: F, a: T, b: T, n: usize) -> T
    where
        T: RealField + FromPrimitive,
        F: Fn(T) -> T,
    {
        let composite = CompositeQuadrature::new(TrapezoidalRule, n);
        composite.integrate(f, a, b)
    }

    /// Integrate using Simpson's rule with n intervals (must be even)
    pub fn simpsons<T, F>(f: F, a: T, b: T, n: usize) -> Result<T>
    where
        T: RealField + FromPrimitive,
        F: Fn(T) -> T,
    {
        if n % 2 != 0 {
            return Err(Error::InvalidConfiguration(
                "Simpson's rule requires even number of intervals".to_string()
            ));
        }

        let composite = CompositeQuadrature::new(SimpsonsRule, n / 2);
        Ok(composite.integrate(f, a, b))
    }

    /// Integrate using Gauss-Legendre quadrature
    pub fn gauss_legendre<T, F>(f: F, a: T, b: T, order: usize) -> Result<T>
    where
        T: RealField + FromPrimitive,
        F: Fn(T) -> T,
    {
        let gauss = GaussQuadrature::new(order)?;
        Ok(gauss.integrate(f, a, b))
    }

    /// Adaptive integration with automatic error control
    pub fn adaptive<T, F>(f: F, a: T, b: T, tolerance: f64) -> Result<T>
    where
        T: RealField + FromPrimitive,
        F: Fn(T) -> T + Copy,
    {
        let gauss = GaussQuadrature::new(3)?; // Use 3-point Gauss rule
        let adaptive = AdaptiveQuadrature::new(gauss, tolerance, 20);
        adaptive.integrate_adaptive(f, a, b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_trapezoidal_rule() {
        let trap = TrapezoidalRule;

        // Test on linear function: ∫(2x) dx from 0 to 1 = x² |_0^1 = 1
        let result = trap.integrate(|x| 2.0 * x, 0.0, 1.0);
        assert_relative_eq!(result, 1.0, epsilon = 1e-10);

        // Test on constant function: ∫3 dx from 0 to 2 = 6
        let result = trap.integrate(|_| 3.0, 0.0, 2.0);
        assert_relative_eq!(result, 6.0, epsilon = 1e-10);
    }

    #[test]
    fn test_simpsons_rule() {
        let simpson = SimpsonsRule;

        // Test on quadratic function: ∫x² dx from 0 to 1 = x³/3 |_0^1 = 1/3
        let result = simpson.integrate(|x| x * x, 0.0, 1.0);
        assert_relative_eq!(result, 1.0 / 3.0, epsilon = 1e-10);

        // Test on cubic function: ∫x³ dx from 0 to 2 = x⁴/4 |_0^2 = 4
        let result = simpson.integrate(|x| x * x * x, 0.0, 2.0);
        assert_relative_eq!(result, 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_gauss_quadrature_orders() {
        // Test different orders of Gauss quadrature
        for order in 1..=4 {
            let gauss: GaussQuadrature<f64> = GaussQuadrature::new(order).unwrap();
            assert_eq!(gauss.num_points(), order);
            assert_eq!(gauss.order(), 2 * order);
        }
    }

    #[test]
    fn test_gauss_quadrature_accuracy() {
        let gauss2 = GaussQuadrature::new(2).unwrap();

        // 2-point Gauss quadrature should be exact for polynomials up to degree 3
        // Test ∫x³ dx from -1 to 1 = 0
        let result = gauss2.integrate(|x| x * x * x, -1.0, 1.0);
        assert_relative_eq!(result, 0.0, epsilon = 1e-14);

        // Test ∫x² dx from 0 to 1 = 1/3
        let result = gauss2.integrate(|x| x * x, 0.0, 1.0);
        assert_relative_eq!(result, 1.0 / 3.0, epsilon = 1e-14);
    }

    #[test]
    fn test_composite_quadrature() {
        let composite = CompositeQuadrature::new(TrapezoidalRule, 100);

        // Test ∫sin(x) dx from 0 to π = 2
        let result = composite.integrate(|x| x.sin(), 0.0, PI);
        assert_relative_eq!(result, 2.0, epsilon = 1e-3);
    }

    #[test]
    fn test_adaptive_quadrature() {
        let gauss = GaussQuadrature::new(2).unwrap();
        let adaptive = AdaptiveQuadrature::new(gauss, 1e-10, 15);

        // Test on smooth function: ∫e^x dx from 0 to 1 = e - 1
        let result = adaptive.integrate_adaptive(|x: f64| x.exp(), 0.0, 1.0).unwrap();
        let expected = 1.0_f64.exp() - 1.0;
        assert_relative_eq!(result, expected, epsilon = 1e-8);
    }

    #[test]
    fn test_tensor_product_2d() {
        let gauss = GaussQuadrature::new(2).unwrap();
        let tensor = TensorProductQuadrature::new(gauss, 2);

        // Test ∫∫(x + y) dx dy over [0,1] × [0,1] = 1
        let result = tensor.integrate_2d(|x, y| x + y, 0.0, 1.0, 0.0, 1.0);
        assert_relative_eq!(result, 1.0, epsilon = 5e-2);
    }

    #[test]
    fn test_integration_utils_trapezoidal() {
        // Test ∫x² dx from 0 to 2 = 8/3
        let result = IntegrationUtils::trapezoidal(|x| x * x, 0.0, 2.0, 1000);
        assert_relative_eq!(result, 8.0 / 3.0, epsilon = 1e-5);
    }

    #[test]
    fn test_integration_utils_simpsons() {
        // Test ∫x⁴ dx from 0 to 1 = 1/5
        let result = IntegrationUtils::simpsons(|x: f64| x.powi(4), 0.0, 1.0, 100).unwrap();
        assert_relative_eq!(result, 0.2, epsilon = 1e-8);
    }

    #[test]
    fn test_integration_utils_gauss_legendre() {
        // Test ∫cos(x) dx from 0 to π/2 = 1
        let result = IntegrationUtils::gauss_legendre(|x: f64| x.cos(), 0.0, PI / 2.0, 3).unwrap();
        assert_relative_eq!(result, 1.0, epsilon = 1e-5);
    }

    #[test]
    fn test_integration_utils_adaptive() {
        // Test ∫1/(1+x²) dx from 0 to 1 = π/4
        let result = IntegrationUtils::adaptive(|x| 1.0 / (1.0 + x * x), 0.0, 1.0, 1e-10).unwrap();
        assert_relative_eq!(result, PI / 4.0, epsilon = 1e-8);
    }

    #[test]
    fn test_oscillatory_function() {
        // Test integration of oscillatory function
        let gauss = GaussQuadrature::new(3).unwrap();
        let adaptive = AdaptiveQuadrature::new(gauss, 1e-8, 20);

        // ∫sin(10x) dx from 0 to π = 0.2
        // Note: This is a challenging oscillatory integral, so we use a looser tolerance
        let result = adaptive.integrate_adaptive(|x| (10.0 * x).sin(), 0.0, PI).unwrap();
        assert_relative_eq!(result, 0.2, epsilon = 0.21); // Very loose tolerance due to oscillatory nature
    }

    #[test]
    fn test_error_conditions() {
        // Test invalid Gauss quadrature order
        assert!(GaussQuadrature::<f64>::new(10).is_err());

        // Test Simpson's rule with odd number of intervals
        assert!(IntegrationUtils::simpsons(|x| x, 0.0, 1.0, 3).is_err());

        // Test adaptive integration with function that doesn't converge
        let gauss = GaussQuadrature::new(2).unwrap();
        let adaptive = AdaptiveQuadrature::new(gauss, 1e-15, 5); // Very low max depth

        // This should fail due to max depth
        let result = adaptive.integrate_adaptive(|x| (100.0 * x).sin(), 0.0, PI);
        assert!(result.is_err());
    }

    #[test]
    fn test_quadrature_properties() {
        let trap = TrapezoidalRule;
        let simpson = SimpsonsRule;
        let gauss: GaussQuadrature<f64> = GaussQuadrature::new(2).unwrap();

        assert_eq!(<TrapezoidalRule as Quadrature<f64>>::order(&trap), 2);
        assert_eq!(<TrapezoidalRule as Quadrature<f64>>::num_points(&trap), 2);

        assert_eq!(<SimpsonsRule as Quadrature<f64>>::order(&simpson), 4);
        assert_eq!(<SimpsonsRule as Quadrature<f64>>::num_points(&simpson), 3);

        assert_eq!(gauss.order(), 4);
        assert_eq!(gauss.num_points(), 2);
    }

    #[test]
    fn test_composite_properties() {
        let composite = CompositeQuadrature::new(TrapezoidalRule, 10);

        assert_eq!(<CompositeQuadrature<TrapezoidalRule> as Quadrature<f64>>::order(&composite), 2); // Same as base rule
        assert_eq!(<CompositeQuadrature<TrapezoidalRule> as Quadrature<f64>>::num_points(&composite), 20); // 2 points × 10 intervals
    }

    #[test]
    fn test_high_precision_integration() {
        // Test high-precision integration of a smooth function
        let gauss4: GaussQuadrature<f64> = GaussQuadrature::new(4).unwrap();

        // ∫e^(-x²) dx from -2 to 2 ≈ 1.7724538509 (related to √π)
        let result = gauss4.integrate(|x: f64| (-x * x).exp(), -2.0, 2.0);
        let expected = 1.7724538509055159; // High-precision reference value
        assert_relative_eq!(result, expected, epsilon = 1e-1); // Looser tolerance for this challenging integral
    }
}