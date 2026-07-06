//! Norm operations for iterators

use eunomia::RealField;

/// Extension trait for norm operations on iterators
pub trait NormIteratorExt: Iterator {
    /// Compute L2 norm without cloning
    fn l2_norm<T>(self) -> T
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy,
    {
        self.map(|x| x * x).fold(T::ZERO, |acc, x| acc + x).sqrt()
    }

    /// Compute L1 norm without cloning
    fn l1_norm<T>(self) -> T
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy,
    {
        self.map(|x| x.abs()).fold(T::ZERO, |acc, x| acc + x)
    }

    /// Compute L∞ norm without cloning
    fn linf_norm<T>(self) -> T
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy,
    {
        self.map(|x| x.abs())
            .fold(T::ZERO, |acc, x| if x > acc { x } else { acc })
    }
}

impl<T: Iterator> NormIteratorExt for T {}
