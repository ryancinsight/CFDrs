//! Zero-copy iterator utilities for CFD operations
//!
//! This module provides efficient iterator combinators optimized for CFD computations,
//! strictly following zero-copy principles with proper borrowing.

use nalgebra::{RealField, DVector};
use rayon::prelude::*;
use std::iter::Iterator;
use num_traits::FromPrimitive;

/// Extension trait for mathematical operations on iterators
pub trait MathIteratorExt: Iterator {
    /// Compute L2 norm without cloning
    fn l2_norm<T>(self) -> T
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy,
    {
        self.map(|x| x * x)
            .fold(T::zero(), |acc, x| acc + x)
            .sqrt()
    }

    /// Compute L1 norm without cloning
    fn l1_norm<T>(self) -> T
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy,
    {
        self.map(|x| x.abs())
            .fold(T::zero(), |acc, x| acc + x)
    }

    /// Compute L∞ norm without cloning
    fn linf_norm<T>(self) -> T
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy,
    {
        self.map(|x| x.abs())
            .fold(T::zero(), |acc, x| if x > acc { x } else { acc })
    }

    /// Compute mean using single-pass iterator
    fn mean<T>(self) -> Option<T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy + FromPrimitive,
    {
        let (sum, count) = self.fold((T::zero(), 0), |(sum, count), x| (sum + x, count + 1));
        if count > 0 {
            Some(sum / T::from_usize(count).unwrap_or_else(T::zero))
        } else {
            None
        }
    }

    /// Compute variance using Welford's online algorithm
    fn variance<T>(self) -> Option<T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy + FromPrimitive,
    {
        let (count, _mean, m2) = self.fold(
            (0, T::zero(), T::zero()),
            |(count, mean, m2), x| {
                let current_count = count + 1;
                let delta = x - mean;
                let current_mean = mean + delta / T::from_usize(current_count).unwrap_or_else(T::zero);
                let delta2 = x - current_mean;
                let current_m2 = m2 + delta * delta2;
                (current_count, current_mean, current_m2)
            }
        );
        
        if count > 1 {
            Some(m2 / T::from_usize(count - 1).unwrap_or_else(T::zero))
        } else {
            None
        }
    }

    /// Apply windowed operations for finite difference computations
    fn windowed_diff<T, F>(self, window_size: usize, f: F) -> WindowedDiff<Self, F, T>
    where
        Self: Iterator<Item = T> + Sized,
        F: Fn(&[T]) -> T,
        T: Clone,
    {
        WindowedDiff::new(self, window_size, f)
    }

    /// Compute running average without cloning
    fn running_average<T>(self) -> RunningAverage<Self, T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy + FromPrimitive,
    {
        RunningAverage::new(self)
    }

    /// Compute exponential moving average with given alpha
    fn exponential_moving_average<T>(self, alpha: T) -> ExponentialMovingAverage<Self, T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy,
    {
        ExponentialMovingAverage::new(self, alpha)
    }
}

// Blanket implementation for all iterators
impl<I: Iterator> MathIteratorExt for I {}

/// Running average iterator
pub struct RunningAverage<I, T> {
    iter: I,
    count: usize,
    sum: T,
}

impl<I, T> RunningAverage<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy + FromPrimitive,
{
    fn new(iter: I) -> Self {
        Self {
            iter,
            count: 0,
            sum: T::zero(),
        }
    }
}

impl<I, T> Iterator for RunningAverage<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy + FromPrimitive,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|x| {
            self.count += 1;
            self.sum = self.sum + x;
            self.sum / T::from_usize(self.count).unwrap_or_else(T::zero)
        })
    }
}

/// Exponential moving average iterator
pub struct ExponentialMovingAverage<I, T> {
    iter: I,
    alpha: T,
    ema: Option<T>,
}

impl<I, T> ExponentialMovingAverage<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy,
{
    fn new(iter: I, alpha: T) -> Self {
        Self {
            iter,
            alpha,
            ema: None,
        }
    }
}

impl<I, T> Iterator for ExponentialMovingAverage<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|x| {
            match self.ema {
                None => {
                    self.ema = Some(x);
                    x
                }
                Some(ref mut prev) => {
                    *prev = self.alpha * x + (T::one() - self.alpha) * *prev;
                    *prev
                }
            }
        })
    }
}

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
pub trait CfdIteratorChain: Iterator {
    /// Chain with finite difference operation
    fn finite_diff<T>(self, spacing: T) -> FiniteDiff<Self, T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy,
    {
        FiniteDiff::new(self, spacing)
    }

    /// Apply smoothing using moving average
    fn smooth<T>(self, window_size: usize) -> Smooth<Self, T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy + FromPrimitive,
    {
        Smooth::new(self, window_size)
    }
}

impl<I: Iterator> CfdIteratorChain for I {}

/// Finite difference iterator
pub struct FiniteDiff<I, T> {
    iter: I,
    spacing: T,
    prev: Option<T>,
}

impl<I, T> FiniteDiff<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy,
{
    fn new(iter: I, spacing: T) -> Self {
        Self {
            iter,
            spacing,
            prev: None,
        }
    }
}

impl<I, T> Iterator for FiniteDiff<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.iter.next() {
                Some(current) => {
                    match self.prev {
                        None => {
                            self.prev = Some(current);
                            continue;
                        }
                        Some(p) => {
                            let diff = (current - p) / self.spacing;
                            self.prev = Some(current);
                            return Some(diff);
                        }
                    }
                }
                None => return None,
            }
        }
    }
}

/// Smoothing iterator
pub struct Smooth<I, T> {
    iter: I,
    window_size: usize,
    buffer: Vec<T>,
    sum: T,
}

impl<I, T> Smooth<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy + FromPrimitive,
{
    fn new(iter: I, window_size: usize) -> Self {
        Self {
            iter,
            window_size,
            buffer: Vec::with_capacity(window_size),
            sum: T::zero(),
        }
    }
}

impl<I, T> Iterator for Smooth<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy + FromPrimitive,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|x| {
            if self.buffer.len() == self.window_size {
                self.sum = self.sum - self.buffer[0];
                self.buffer.remove(0);
            }
            self.buffer.push(x);
            self.sum = self.sum + x;
            self.sum / T::from_usize(self.buffer.len()).unwrap_or_else(T::one)
        })
    }
}

/// Field operations trait for CFD fields
pub trait CfdFieldOps: Iterator {
    /// Compute gradient using central differences
    fn gradient<T>(self, spacing: T) -> Gradient<Self, T>
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy,
    {
        Gradient::new(self, spacing)
    }

    /// Compute divergence (for vector fields)
    fn divergence<T>(self, spacing: T) -> T
    where
        Self: Iterator<Item = T> + Sized,
        T: RealField + Copy,
    {
        self.gradient(spacing).fold(T::zero(), |acc, x| acc + x)
    }
}

impl<I: Iterator> CfdFieldOps for I {}

/// Gradient iterator
pub struct Gradient<I, T> {
    iter: I,
    spacing: T,
    buffer: Vec<T>,
}

impl<I, T> Gradient<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy,
{
    fn new(iter: I, spacing: T) -> Self {
        Self {
            iter,
            spacing,
            buffer: Vec::with_capacity(3),
        }
    }
}

impl<I, T> Iterator for Gradient<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.iter.next() {
                Some(x) => {
                    self.buffer.push(x);
                    if self.buffer.len() == 3 {
                        let grad = (self.buffer[2] - self.buffer[0]) / (self.spacing + self.spacing);
                        self.buffer.remove(0);
                        return Some(grad);
                    }
                }
                None => return None,
            }
        }
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