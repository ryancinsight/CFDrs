//! Zero-copy iterator utilities for CFD operations
//!
//! This module provides efficient iterator combinators optimized for CFD computations,
//! strictly following zero-copy principles with proper borrowing.

use nalgebra::{RealField, DVector};
use rayon::prelude::*;
use std::iter::Iterator;
use num_traits::FromPrimitive;

/// Extension trait for mathematical operations on iterators
pub trait MathIteratorExt: Iterator
where
    Self::Item: RealField + Send + Sync + Copy,
    Self: Sized,
{
    /// Compute L2 norm without cloning
    fn l2_norm(self) -> Self::Item {
        self.map(|x| x * x)
            .fold(Self::Item::zero(), |acc, x| acc + x)
            .sqrt()
    }

    /// Compute L1 norm without cloning
    fn l1_norm(self) -> Self::Item {
        self.map(|x| x.abs())
            .fold(Self::Item::zero(), |acc, x| acc + x)
    }

    /// Compute L∞ norm without cloning
    fn linf_norm(self) -> Self::Item {
        self.map(|x| x.abs())
            .fold(Self::Item::zero(), |acc, x| if x > acc { x } else { acc })
    }

    /// Compute mean using single-pass iterator
    fn mean(self) -> Option<Self::Item> {
        let (sum, count) = self.fold((Self::Item::zero(), 0), |(sum, count), x| (sum + x, count + 1));
        if count > 0 {
            Some(sum / Self::Item::from_usize(count).unwrap_or_else(|| Self::Item::zero()))
        } else {
            None
        }
    }

    /// Compute variance using Welford's online algorithm
    fn variance(self) -> Option<Self::Item> {
        let (count, _mean, m2) = self.fold(
            (0, Self::Item::zero(), Self::Item::zero()),
            |(count, mean, m2), x| {
                let current_count = count + 1;
                let delta = x - mean;
                let current_mean = mean + delta / Self::Item::from_usize(current_count).unwrap_or_else(|| Self::Item::zero());
                let delta2 = x - current_mean;
                let current_m2 = m2 + delta * delta2;
                (current_count, current_mean, current_m2)
            }
        );
        
        if count > 1 {
            Some(m2 / Self::Item::from_usize(count - 1).unwrap_or_else(|| Self::Item::zero()))
        } else {
            None
        }
    }
}

// Blanket implementation for all iterators with appropriate items
impl<I> MathIteratorExt for I
where
    I: Iterator,
    I::Item: RealField + Send + Sync + Copy,
{}

/// Efficient vector operations using zero-copy principles
pub struct VectorOps;

impl VectorOps {
    /// Add two vectors without cloning - uses borrowing
    pub fn add<T: RealField + Copy>(a: &DVector<T>, b: &DVector<T>) -> Result<DVector<T>, &'static str> {
        if a.len() != b.len() {
            return Err("Vector dimensions must match");
        }
        
        let result: Vec<T> = a.iter()
            .zip(b.iter())
            .map(|(&x, &y)| x + y)  // Use references, no cloning!
            .collect();
        
        Ok(DVector::from_vec(result))
    }

    /// Subtract vectors without cloning
    pub fn sub<T: RealField + Copy>(a: &DVector<T>, b: &DVector<T>) -> Result<DVector<T>, &'static str> {
        if a.len() != b.len() {
            return Err("Vector dimensions must match");
        }
        
        let result: Vec<T> = a.iter()
            .zip(b.iter())
            .map(|(&x, &y)| x - y)
            .collect();
        
        Ok(DVector::from_vec(result))
    }

    /// Parallel add using rayon for large vectors
    pub fn parallel_add<T>(a: &DVector<T>, b: &DVector<T>) -> Result<DVector<T>, &'static str>
    where
        T: RealField + Copy + Send + Sync,
    {
        if a.len() != b.len() {
            return Err("Vector dimensions must match");
        }
        
        let result: Vec<T> = a.as_slice()
            .par_iter()
            .zip(b.as_slice().par_iter())
            .map(|(&x, &y)| x + y)
            .collect();
        
        Ok(DVector::from_vec(result))
    }

    /// Dot product without cloning
    pub fn dot<T: RealField + Copy>(a: &DVector<T>, b: &DVector<T>) -> Result<T, &'static str> {
        if a.len() != b.len() {
            return Err("Vector dimensions must match");
        }
        
        Ok(a.iter()
            .zip(b.iter())
            .map(|(&x, &y)| x * y)
            .fold(T::zero(), |acc, x| acc + x))
    }
}

/// Windowed operations for finite differences
pub struct WindowedOps;

impl WindowedOps {
    /// Apply a windowed operation without cloning
    pub fn apply<T, F>(data: &[T], window_size: usize, f: F) -> Vec<T>
    where
        T: RealField + Copy,
        F: Fn(&[T]) -> T,
    {
        data.windows(window_size)
            .map(|window| f(window))
            .collect()
    }

    /// Central difference without cloning
    pub fn central_diff<T: RealField + Copy>(data: &[T], spacing: T) -> Vec<T> {
        if data.len() < 3 {
            return vec![];
        }
        
        data.windows(3)
            .map(|w| (w[2] - w[0]) / (spacing + spacing))
            .collect()
    }

    /// Forward difference without cloning
    pub fn forward_diff<T: RealField + Copy>(data: &[T], spacing: T) -> Vec<T> {
        if data.len() < 2 {
            return vec![];
        }
        
        data.windows(2)
            .map(|w| (w[1] - w[0]) / spacing)
            .collect()
    }
}

/// Parallel reduction operations
pub struct ParallelOps;

impl ParallelOps {
    /// Parallel sum without cloning
    pub fn sum<T>(data: &[T]) -> T
    where
        T: RealField + Copy + Send + Sync,
    {
        data.par_iter()
            .copied()
            .reduce(|| T::zero(), |a, b| a + b)
    }

    /// Parallel norm computation
    pub fn l2_norm<T>(data: &[T]) -> T
    where
        T: RealField + Copy + Send + Sync,
    {
        data.par_iter()
            .map(|&x| x * x)
            .reduce(|| T::zero(), |a, b| a + b)
            .sqrt()
    }

    /// Parallel maximum
    pub fn max<T>(data: &[T]) -> Option<T>
    where
        T: RealField + Copy + Send + Sync,
    {
        data.par_iter()
            .copied()
            .reduce_with(|a, b| if a > b { a } else { b })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_l2_norm() {
        let data = vec![3.0, 4.0];
        let norm = data.into_iter().l2_norm();
        assert_relative_eq!(norm, 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_vector_add() {
        let a = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let b = DVector::from_vec(vec![4.0, 5.0, 6.0]);
        let c = VectorOps::add(&a, &b).unwrap();
        assert_eq!(c[0], 5.0);
        assert_eq!(c[1], 7.0);
        assert_eq!(c[2], 9.0);
    }

    #[test]
    fn test_central_diff() {
        let data = vec![1.0, 4.0, 9.0, 16.0];
        let diff = WindowedOps::central_diff(&data, 1.0);
        assert_eq!(diff.len(), 2);
        assert_relative_eq!(diff[0], 4.0, epsilon = 1e-10); // (9-1)/2
        assert_relative_eq!(diff[1], 6.0, epsilon = 1e-10); // (16-4)/2
    }
}