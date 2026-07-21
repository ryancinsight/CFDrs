//! Quadrature rules for numerical integration

use crate::integration::traits::Quadrature;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};

// Quadrature constants
const TWO: f64 = 2.0;
const FOUR: f64 = 4.0;
const FIVE: f64 = 5.0;
const SIX: f64 = 6.0;
const SEVEN: f64 = 7.0;
const EIGHT: f64 = 8.0;
const NINE: f64 = 9.0;
const SQRT_THREE_INV: f64 = 0.577_350_269_189_626; // 1/sqrt(3)

fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

/// Trapezoidal rule for numerical integration
pub struct TrapezoidalRule;

impl<T: RealField + FloatElement + Copy> Quadrature<T> for TrapezoidalRule {
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let two = from_f64(TWO);
        (b - a) * (f(a) + f(b)) / two
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

impl<T: RealField + FloatElement + Copy> Quadrature<T> for SimpsonsRule {
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let two: T = from_f64(TWO);
        let four: T = from_f64(FOUR);
        let six: T = from_f64(SIX);

        let mid = (a + b) / two;
        (b - a) * (f(a) + four * f(mid) + f(b)) / six
    }

    fn order(&self) -> usize {
        4
    }

    fn num_points(&self) -> usize {
        3
    }
}

/// Gauss-Legendre quadrature
pub struct GaussQuadrature<T: RealField + Copy> {
    points: Vec<T>,
    weights: Vec<T>,
    order: usize,
}

impl<T: RealField + Copy + FloatElement> GaussQuadrature<T> {
    /// Create a new Gauss quadrature rule with specified order
    pub fn new(order: usize) -> Result<Self> {
        if order == 0 || order > 5 {
            return Err(Error::InvalidInput(format!(
                "Gauss quadrature order must be between 1 and 5, got {order}"
            )));
        }

        let (points, weights) = Self::get_points_and_weights(order)?;

        Ok(Self {
            points,
            weights,
            order,
        })
    }

    /// Get Gauss-Legendre points and weights for given order
    fn get_points_and_weights(order: usize) -> Result<(Vec<T>, Vec<T>)> {
        match order {
            1 => Self::gauss_1_point(),
            2 => Self::gauss_2_point(),
            3 => Self::gauss_3_point(),
            4 => Self::gauss_4_point(),
            5 => Self::gauss_5_point(),
            _ => Err(Error::InvalidInput(format!(
                "Unsupported Gauss quadrature order: {order}"
            ))),
        }
    }

    fn gauss_1_point() -> Result<(Vec<T>, Vec<T>)> {
        Ok((vec![T::ZERO], vec![from_f64(2.0)]))
    }

    fn gauss_2_point() -> Result<(Vec<T>, Vec<T>)> {
        let sqrt3_inv: T = from_f64(SQRT_THREE_INV);
        Ok((vec![-sqrt3_inv, sqrt3_inv], vec![T::ONE, T::ONE]))
    }

    fn gauss_3_point() -> Result<(Vec<T>, Vec<T>)> {
        let sqrt15: T = from_f64(15.0_f64.sqrt());
        let five: T = from_f64(FIVE);
        let nine: T = from_f64(NINE);
        let eight: T = from_f64(EIGHT);

        Ok((
            vec![-sqrt15 / five, T::ZERO, sqrt15 / five],
            vec![five / nine, eight / nine, five / nine],
        ))
    }

    fn gauss_4_point() -> Result<(Vec<T>, Vec<T>)> {
        let sqrt6_5 = (6.0 / 5.0_f64).sqrt();
        let term1 = <T as NumericElement>::sqrt(from_f64((3.0 - 2.0 * sqrt6_5) / 7.0));
        let term2 = <T as NumericElement>::sqrt(from_f64((3.0 + 2.0 * sqrt6_5) / 7.0));

        let sqrt30 = 30.0_f64.sqrt();
        let w1: T = from_f64((18.0 + sqrt30) / 36.0);
        let w2: T = from_f64((18.0 - sqrt30) / 36.0);

        Ok((vec![-term2, -term1, term1, term2], vec![w2, w1, w1, w2]))
    }

    fn gauss_5_point() -> Result<(Vec<T>, Vec<T>)> {
        // 5-point Gauss-Legendre quadrature
        let sqrt10_7: T = from_f64((10.0 / 7.0_f64).sqrt());
        let seven: T = from_f64(SEVEN);
        let five: T = from_f64(FIVE);
        let two: T = from_f64(2.0);
        let one_third: T = from_f64(1.0 / 3.0);

        let term = two * sqrt10_7 / seven;

        let x1 = one_third * <T as NumericElement>::sqrt(five - term);
        let x2 = one_third * <T as NumericElement>::sqrt(five + term);

        let w0: T = from_f64(128.0 / 225.0);
        let w1: T = from_f64((322.0 + 13.0 * 70.0_f64.sqrt()) / 900.0);
        let w2: T = from_f64((322.0 - 13.0 * 70.0_f64.sqrt()) / 900.0);

        Ok((vec![-x2, -x1, T::ZERO, x1, x2], vec![w2, w1, w0, w1, w2]))
    }

    /// Create default 2-point Gauss quadrature
    pub fn default() -> Result<Self> {
        Self::new(2)
    }
}

impl<T: RealField + FloatElement + Copy> Quadrature<T> for GaussQuadrature<T> {
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let two = from_f64(TWO);
        let half_interval = (b - a) / two;
        let mid_point = (a + b) / two;

        // Use iterator combinators for zero-copy optimization
        let result = self
            .points
            .iter()
            .zip(self.weights.iter())
            .map(|(point, weight)| {
                let x = mid_point + half_interval * *point;
                *weight * f(x)
            })
            .fold(T::ZERO, |acc, term| acc + term);

        result * half_interval
    }

    fn order(&self) -> usize {
        2 * self.order
    }

    fn num_points(&self) -> usize {
        self.points.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn simpsons_integrates_quadratic_exactly() {
        let rule = SimpsonsRule;

        let integral = rule.integrate(|x: f64| x * x, 0.0, 1.0);

        assert_relative_eq!(integral, 1.0 / 3.0, epsilon = 1e-15);
    }

    #[test]
    fn gauss_legendre_integrates_cubic_exactly() {
        let rule = GaussQuadrature::<f64>::new(2).unwrap();

        let integral = rule.integrate(|x| x * x * x + x * x, -1.0, 1.0);

        assert_relative_eq!(integral, 2.0 / 3.0, epsilon = 1e-15);
        assert_eq!(rule.order(), 4);
        assert_eq!(rule.num_points(), 2);
    }

    #[test]
    fn gauss_quadrature_rejects_unsupported_order() {
        let Err(error) = GaussQuadrature::<f64>::new(6) else {
            panic!("expected order 6 to be rejected");
        };

        assert!(
            error.to_string().contains("between 1 and 5"),
            "unexpected error: {error}"
        );
    }
}
