//! Statistical operations for iterators

use nalgebra::RealField;
use num_traits::FromPrimitive;
/// Extension trait for statistical operations on iterators
pub trait StatisticsIteratorExt: Iterator {
    /// Compute mean using single-pass iterator
    fn mean<T>(self) -> Option<T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy + FromPrimitive,
    {
        let (sum, count) = self.fold((T::zero(), 0), |(sum, count), x| (sum + x, count + 1));
        if count > 0 {
            T::from_usize(count).map(|n| sum / n)
        } else {
            None
        }
    }
    /// Compute variance using Welford's online algorithm
    fn variance<T>(self) -> Option<T>
        let (count, _mean, m2) =
            self.fold((0usize, T::zero(), T::zero()), |(count, mean, m2), x| {
                let current_count = count + 1;
                let delta = x - mean;
                let current_mean =
                    mean + delta / T::from_usize(current_count).unwrap_or_else(T::zero);
                let delta2 = x - current_mean;
                let current_m2 = m2 + delta * delta2;
                (current_count, current_mean, current_m2)
            });
        if count > 1 {
            T::from_usize(count - 1).map(|n| m2 / n)
    /// Compute standard deviation
    fn std_dev<T>(self) -> Option<T>
        self.variance().map(|v| v.sqrt())
}
impl<T: Iterator> StatisticsIteratorExt for T {}
