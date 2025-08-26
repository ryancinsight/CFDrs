//! Composite quadrature rules for integration over multiple intervals
use cfd_core::error::{Error, NumericalErrorKind};
use cfd_core::numeric;

use crate::integration::traits::Quadrature;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
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
impl<T, Q> Quadrature<T> for CompositeQuadrature<Q>
where
    T: RealField + From<f64> + FromPrimitive + Copy,
    Q: Quadrature<T>,
{
    fn integrate<F>(&self, f: F, a: T, b: T) -> T
    where
        F: Fn(T) -> T,
    {
        let n = cfd_core::numeric::from_usize(self.num_intervals)?;
        let h = (b - a) / n;
        // Use iterator range with fold for zero-copy optimization
        (0..self.num_intervals)
            .map(|i| {
                let xi = a + cfd_core::numeric::from_usize(i)? * h;
                let xi_plus_1 = xi + h;
                self.base_rule.integrate(&f, xi, xi_plus_1)
            })
            .fold(T::zero(), |acc, integral| acc + integral)
    fn order(&self) -> usize {
        self.base_rule.order()
    fn num_points(&self) -> usize {
        self.base_rule.num_points() * self.num_intervals
