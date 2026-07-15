//! Composite quadrature rules for integration over multiple intervals

use crate::integration::traits::Quadrature;
use eunomia::{FloatElement, NumericElement, RealField};

fn from_usize<T: FloatElement>(value: usize) -> T {
    let value_u64 = u64::try_from(value).expect("invariant: usize value fits in u64");
    <T as FloatElement>::from_f64(<u64 as NumericElement>::to_f64(value_u64))
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
    T: RealField + FloatElement + Copy,
    Q: Quadrature<T>,
{
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let n = from_usize::<T>(self.num_intervals);
        let h = (b - a) / n;

        // Use iterator range with fold for zero-copy optimization
        (0..self.num_intervals)
            .map(|i| {
                let xi = a + from_usize::<T>(i) * h;
                let xi_plus_1 = xi + h;
                self.base_rule.integrate(&f, xi, xi_plus_1)
            })
            .fold(T::ZERO, |acc, integral| acc + integral)
    }

    fn order(&self) -> usize {
        self.base_rule.order()
    }

    fn num_points(&self) -> usize {
        self.base_rule.num_points() * self.num_intervals
    }
}
