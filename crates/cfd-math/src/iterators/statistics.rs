//! Statistical operations for iterators

use eunomia::{FloatElement, NumericElement, RealField};

fn from_usize<T: FloatElement>(value: usize) -> T {
    let value_u64 =
        u64::try_from(value).expect("invariant: usize fits into u64 on supported targets");
    <T as FloatElement>::from_f64(<u64 as NumericElement>::to_f64(value_u64))
}

/// Extension trait for statistical operations on iterators
pub trait StatisticsIteratorExt: Iterator {
    /// Compute mean using single-pass iterator
    fn mean<T>(self) -> Option<T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + FloatElement + Copy,
    {
        let (sum, count) = self.fold((T::ZERO, 0), |(sum, count), x| (sum + x, count + 1));
        if count > 0 {
            Some(sum / from_usize(count))
        } else {
            None
        }
    }

    /// Compute variance using Welford's online algorithm
    fn variance<T>(self) -> Option<T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + FloatElement + Copy,
    {
        let (count, _mean, m2) = self.fold((0usize, T::ZERO, T::ZERO), |(count, mean, m2), x| {
            let current_count = count + 1;
            let delta = x - mean;
            let current_mean = mean + delta / from_usize(current_count);
            let delta2 = x - current_mean;
            let current_m2 = m2 + delta * delta2;
            (current_count, current_mean, current_m2)
        });

        if count > 1 {
            Some(m2 / from_usize(count - 1))
        } else {
            None
        }
    }

    /// Compute standard deviation
    fn std_dev<T>(self) -> Option<T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + FloatElement + Copy,
    {
        self.variance().map(NumericElement::sqrt)
    }
}

impl<T: Iterator> StatisticsIteratorExt for T {}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn mean_returns_none_for_empty_iterator() {
        let values: Vec<f64> = Vec::new();

        assert_eq!(values.into_iter().mean::<f64>(), None);
    }

    #[test]
    fn mean_computes_arithmetic_average() {
        let values = vec![1.0, 2.0, 4.0, 5.0];

        assert_relative_eq!(
            values.into_iter().mean::<f64>().unwrap(),
            3.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn variance_uses_sample_denominator() {
        let values = vec![1.0, 2.0, 3.0, 4.0];

        assert_relative_eq!(
            values.into_iter().variance::<f64>().unwrap(),
            5.0 / 3.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn standard_deviation_is_variance_square_root() {
        let values = vec![1.0, 2.0, 3.0, 4.0];

        assert_relative_eq!(
            values.into_iter().std_dev::<f64>().unwrap(),
            (5.0_f64 / 3.0).sqrt(),
            epsilon = 1e-15
        );
    }
}
