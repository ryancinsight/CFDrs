//! Parallel iterator operations for CFD computations

use nalgebra::RealField;
use rayon::prelude::*;

/// Extension trait for parallel operations on iterators
pub trait ParallelIteratorExt: ParallelIterator {
    /// Parallel reduction for sum
    fn parallel_sum<T>(self) -> T
    where
        Self: ParallelIterator<Item = T> + Sized,
        T: RealField + Copy + Send + Sync,
    {
        self.reduce(T::zero, |acc, x| acc + x)
    }

    /// Parallel reduction for maximum
    fn parallel_max<T>(self) -> Option<T>
    where
        Self: ParallelIterator<Item = T> + Sized,
        T: RealField + Copy + Send + Sync,
    {
        self.reduce_with(|a, b| if a > b { a } else { b })
    }

    /// Parallel reduction for minimum
    fn parallel_min<T>(self) -> Option<T>
    where
        Self: ParallelIterator<Item = T> + Sized,
        T: RealField + Copy + Send + Sync,
    {
        self.reduce_with(|a, b| if a < b { a } else { b })
    }

    /// Parallel L2 norm computation
    fn parallel_l2_norm<T>(self) -> T
    where
        Self: ParallelIterator<Item = T> + Sized,
        T: RealField + Copy + Send + Sync,
    {
        self.map(|x| x * x).reduce(T::zero, |acc, x| acc + x).sqrt()
    }
}

impl<T: ParallelIterator> ParallelIteratorExt for T {}
