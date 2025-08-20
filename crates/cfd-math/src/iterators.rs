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

    /// Apply windowed operations for finite difference computations
    fn windowed_diff<F>(self, window_size: usize, f: F) -> WindowedDiff<Self, F, Self::Item>
    where
        F: Fn(&[Self::Item]) -> Self::Item,
    {
        WindowedDiff::new(self, window_size, f)
    }

    /// Compute running average without cloning
    fn running_average(self) -> impl Iterator<Item = Self::Item> {
        let mut count = 0;
        let mut sum = Self::Item::zero();

        self.map(move |x| {
            count += 1;
            sum = sum + x;
            sum / Self::Item::from_usize(count).unwrap_or_else(|| Self::Item::zero())
        })
    }

    /// Compute exponential moving average with given alpha
    fn exponential_moving_average(self, alpha: Self::Item) -> impl Iterator<Item = Self::Item> {
        let mut ema = None;

        self.map(move |x| {
            match ema {
                None => {
                    ema = Some(x);
                    x
                }
                Some(ref mut prev) => {
                    *prev = alpha * x + (Self::Item::one() - alpha) * *prev;
                    *prev
                }
            }
        })
    }
}

// Blanket implementation for all iterators with appropriate items
impl<I> MathIteratorExt for I
where
    I: Iterator,
    I::Item: RealField + Send + Sync + Copy,
{}

/// Windowed difference iterator for finite difference operations
pub struct WindowedDiff<I, F, T> {
    iter: I,
    window_size: usize,
    buffer: Vec<T>,
    buffer_start: usize,
    buffer_len: usize,
    func: F,
}

impl<I, F, T> WindowedDiff<I, F, T>
where
    I: Iterator<Item = T>,
    F: Fn(&[T]) -> T,
    T: Clone,
{
    fn new(iter: I, window_size: usize, func: F) -> Self {
        Self {
            iter,
            window_size,
            buffer: Vec::with_capacity(window_size),
            buffer_start: 0,
            buffer_len: 0,
            func,
        }
    }
}

impl<I, F, T> Iterator for WindowedDiff<I, F, T>
where
    I: Iterator<Item = T>,
    F: Fn(&[T]) -> T,
    T: Clone,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        // Fill buffer to window size
        while self.buffer_len < self.window_size {
            match self.iter.next() {
                Some(val) => {
                    if self.buffer.len() < self.window_size {
                        self.buffer.push(val);
                    } else {
                        self.buffer[self.buffer_len] = val;
                    }
                    self.buffer_len += 1;
                }
                None => {
                    if self.buffer_len == 0 {
                        return None;
                    }
                    break;
                }
            }
        }

        if self.buffer_len < self.window_size {
            return None;
        }

        // Apply function to window
        let result = (self.func)(&self.buffer[self.buffer_start..self.buffer_start + self.window_size]);

        // Slide window
        self.buffer_start = (self.buffer_start + 1) % self.buffer.capacity();
        self.buffer_len -= 1;

        Some(result)
    }
}

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

    /// Compute norm without cloning
    pub fn norm<T: RealField + Copy>(v: &DVector<T>) -> T {
        v.iter()
            .map(|&x| x * x)
            .fold(T::zero(), |acc, x| acc + x)
            .sqrt()
    }

    /// Scale vector without cloning
    pub fn scale<T: RealField + Copy>(v: &DVector<T>, scalar: T) -> DVector<T> {
        let result: Vec<T> = v.iter()
            .map(|&x| x * scalar)
            .collect();
        DVector::from_vec(result)
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

    /// Backward difference without cloning
    pub fn backward_diff<T: RealField + Copy>(data: &[T], spacing: T) -> Vec<T> {
        if data.len() < 2 {
            return vec![];
        }
        
        data.windows(2)
            .map(|w| (w[1] - w[0]) / spacing)
            .collect()
    }

    /// Second derivative using central difference
    pub fn second_derivative<T: RealField + Copy>(data: &[T], spacing: T) -> Vec<T> {
        if data.len() < 3 {
            return vec![];
        }
        
        let spacing_sq = spacing * spacing;
        data.windows(3)
            .map(|w| (w[2] - w[1] - w[1] + w[0]) / spacing_sq)
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

    /// Parallel minimum
    pub fn min<T>(data: &[T]) -> Option<T>
    where
        T: RealField + Copy + Send + Sync,
    {
        data.par_iter()
            .copied()
            .reduce_with(|a, b| if a < b { a } else { b })
    }

    /// Parallel mean computation
    pub fn mean<T>(data: &[T]) -> T
    where
        T: RealField + Copy + Send + Sync + FromPrimitive,
    {
        let sum = Self::sum(data);
        sum / T::from_usize(data.len()).unwrap_or_else(T::one)
    }
}

/// Composable iterator chains for CFD operations
pub trait CfdIteratorChain<T>: Iterator<Item = T>
where
    T: RealField,
    Self: Sized,
{
    /// Chain with finite difference operation
    fn finite_diff(self, spacing: T) -> impl Iterator<Item = T>
    where
        T: Copy,
    {
        let mut prev = None;
        self.filter_map(move |current| {
            match prev {
                None => {
                    prev = Some(current);
                    None
                }
                Some(p) => {
                    let diff = (current - p) / spacing;
                    prev = Some(current);
                    Some(diff)
                }
            }
        })
    }

    /// Apply smoothing using moving average
    fn smooth(self, window_size: usize) -> impl Iterator<Item = T>
    where
        T: Copy + FromPrimitive,
    {
        let mut buffer = Vec::with_capacity(window_size);
        let mut sum = T::zero();
        
        self.map(move |x| {
            if buffer.len() == window_size {
                sum = sum - buffer[0];
                buffer.remove(0);
            }
            buffer.push(x);
            sum = sum + x;
            sum / T::from_usize(buffer.len()).unwrap_or_else(T::one)
        })
    }
}

impl<I, T> CfdIteratorChain<T> for I
where
    I: Iterator<Item = T>,
    T: RealField,
{}

/// Field operations trait for CFD fields
pub trait CfdFieldOps<T>: Iterator<Item = T>
where
    T: RealField,
    Self: Sized,
{
    /// Compute gradient using central differences
    fn gradient(self, spacing: T) -> impl Iterator<Item = T>
    where
        T: Copy,
    {
        let mut buffer = Vec::with_capacity(3);
        
        self.filter_map(move |x| {
            buffer.push(x);
            if buffer.len() == 3 {
                let grad = (buffer[2] - buffer[0]) / (spacing + spacing);
                buffer.remove(0);
                Some(grad)
            } else {
                None
            }
        })
    }

    /// Compute divergence (for vector fields)
    fn divergence(self, spacing: T) -> T
    where
        T: Copy,
    {
        self.gradient(spacing).fold(T::zero(), |acc, x| acc + x)
    }
}

impl<I, T> CfdFieldOps<T> for I
where
    I: Iterator<Item = T>,
    T: RealField,
{}

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

    #[test]
    fn test_parallel_sum() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let sum = ParallelOps::sum(&data);
        assert_relative_eq!(sum, 15.0, epsilon = 1e-10);
    }

    #[test]
    fn test_running_average() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let avgs: Vec<f64> = data.into_iter().running_average().collect();
        assert_relative_eq!(avgs[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(avgs[1], 1.5, epsilon = 1e-10);
        assert_relative_eq!(avgs[2], 2.0, epsilon = 1e-10);
        assert_relative_eq!(avgs[3], 2.5, epsilon = 1e-10);
    }
}