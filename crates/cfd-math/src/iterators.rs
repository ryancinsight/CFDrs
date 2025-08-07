//! Advanced iterator utilities for zero-copy CFD operations.
//!
//! This module provides iterator combinators and utilities optimized for CFD computations,
//! following CUPID principles (Composable, Unix Philosophy, Predictable, Idiomatic, Domain-based).

use nalgebra::{RealField, DVector};
use rayon::prelude::*;
use std::iter::{Iterator, IntoIterator};

/// Extension trait for iterator-based mathematical operations
pub trait MathIteratorExt<T>: Iterator<Item = T>
where
    T: RealField + Send + Sync,
    Self: Sized,
{
    /// Compute L2 norm using zero-copy iterator operations
    fn l2_norm(self) -> T {
        self.map(|x| x.clone() * x.clone())
            .fold(T::zero(), |acc, x| acc + x)
            .sqrt()
    }

    /// Compute L1 norm using iterator combinators
    fn l1_norm(self) -> T {
        self.map(|x| x.abs())
            .fold(T::zero(), |acc, x| acc + x)
    }

    /// Compute L∞ norm using iterator operations
    fn linf_norm(self) -> T {
        self.map(|x| x.abs())
            .fold(T::zero(), |acc, x| if x > acc { x } else { acc })
    }

    /// Compute mean using single-pass iterator
    fn mean(self) -> Option<T> {
        let (sum, count) = self.fold((T::zero(), 0), |(sum, count), x| (sum + x, count + 1));
        if count > 0 {
            Some(sum / T::from_usize(count).unwrap())
        } else {
            None
        }
    }

    /// Compute variance using Welford's online algorithm
    fn variance(self) -> Option<T> {
        let (count, _mean, m2) = self.fold(
            (0, T::zero(), T::zero()),
            |(count, mean, m2), x| {
                let new_count = count + 1;
                let delta = x.clone() - mean.clone();
                let new_mean = mean + delta.clone() / T::from_usize(new_count).unwrap();
                let delta2 = x - new_mean.clone();
                let new_m2 = m2 + delta * delta2;
                (new_count, new_mean, new_m2)
            }
        );
        
        if count > 1 {
            Some(m2 / T::from_usize(count - 1).unwrap())
        } else {
            None
        }
    }

    /// Apply windowed operations for finite difference computations
    fn windowed_diff<F>(self, window_size: usize, f: F) -> WindowedDiff<Self, F, T>
    where
        F: Fn(&[T]) -> T,
    {
        WindowedDiff::new(self, window_size, f)
    }

    // Note: Parallel reduction is available through rayon's ParallelIterator trait directly
}

impl<I, T> MathIteratorExt<T> for I
where
    I: Iterator<Item = T>,
    T: RealField + Send + Sync,
{
}

/// Windowed difference iterator for finite difference operations
pub struct WindowedDiff<I, F, T> {
    iter: I,
    window_size: usize,
    buffer: Vec<T>,
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
        while self.buffer.len() < self.window_size {
            if let Some(item) = self.iter.next() {
                self.buffer.push(item);
            } else {
                return None;
            }
        }

        if self.buffer.len() == self.window_size {
            let result = (self.func)(&self.buffer);
            // Slide window
            self.buffer.remove(0);
            if let Some(item) = self.iter.next() {
                self.buffer.push(item);
            }
            Some(result)
        } else {
            None
        }
    }
}

/// Zero-copy vector operations using iterator combinators
pub struct VectorOps;

/// Zero-copy slice operations for CFD computations
pub struct SliceOps;

impl VectorOps {
    /// Element-wise addition using iterators
    pub fn add<T: RealField>(a: &DVector<T>, b: &DVector<T>) -> Result<DVector<T>, &'static str> {
        if a.len() != b.len() {
            return Err("Vector dimensions must match");
        }

        let result: Vec<T> = a.iter()
            .zip(b.iter())
            .map(|(x, y)| x.clone() + y.clone())
            .collect();

        Ok(DVector::from_vec(result))
    }

    /// Element-wise multiplication using iterators
    pub fn hadamard<T: RealField>(a: &DVector<T>, b: &DVector<T>) -> Result<DVector<T>, &'static str> {
        if a.len() != b.len() {
            return Err("Vector dimensions must match");
        }

        let result: Vec<T> = a.iter()
            .zip(b.iter())
            .map(|(x, y)| x.clone() * y.clone())
            .collect();

        Ok(DVector::from_vec(result))
    }

    /// Dot product using iterator fold
    pub fn dot<T: RealField>(a: &DVector<T>, b: &DVector<T>) -> Result<T, &'static str> {
        if a.len() != b.len() {
            return Err("Vector dimensions must match");
        }

        let result = a.iter()
            .zip(b.iter())
            .map(|(x, y)| x.clone() * y.clone())
            .fold(T::zero(), |acc, x| acc + x);

        Ok(result)
    }

    /// Parallel vector operations for large vectors
    pub fn parallel_add<T>(a: &DVector<T>, b: &DVector<T>) -> Result<DVector<T>, &'static str>
    where
        T: RealField + Send + Sync,
    {
        if a.len() != b.len() {
            return Err("Vector dimensions must match");
        }

        let result: Vec<T> = a.as_slice()
            .par_iter()
            .zip(b.as_slice().par_iter())
            .map(|(x, y)| x.clone() + y.clone())
            .collect();

        Ok(DVector::from_vec(result))
    }

    /// In-place vector addition (zero-copy for destination)
    pub fn add_inplace<T: RealField>(dest: &mut DVector<T>, src: &DVector<T>) -> Result<(), &'static str> {
        if dest.len() != src.len() {
            return Err("Vector dimensions must match");
        }

        dest.iter_mut()
            .zip(src.iter())
            .for_each(|(d, s)| *d += s.clone());

        Ok(())
    }

    /// Zero-copy vector scaling
    pub fn scale_inplace<T: RealField>(vector: &mut DVector<T>, scalar: T) {
        vector.iter_mut().for_each(|x| *x *= scalar.clone());
    }
}

impl SliceOps {
    /// Zero-copy slice-based dot product
    pub fn dot_slice<T: RealField>(a: &[T], b: &[T]) -> Result<T, &'static str> {
        if a.len() != b.len() {
            return Err("Slice dimensions must match");
        }

        let result = a.iter()
            .zip(b.iter())
            .map(|(x, y)| x.clone() * y.clone())
            .fold(T::zero(), |acc, x| acc + x);

        Ok(result)
    }

    /// Zero-copy slice-based L2 norm
    pub fn l2_norm_slice<T: RealField>(slice: &[T]) -> T {
        slice.iter()
            .map(|x| x.clone() * x.clone())
            .fold(T::zero(), |acc, x| acc + x)
            .sqrt()
    }

    /// Zero-copy slice-based element-wise addition
    pub fn add_slice_inplace<T: RealField>(dest: &mut [T], src: &[T]) -> Result<(), &'static str> {
        if dest.len() != src.len() {
            return Err("Slice dimensions must match");
        }

        dest.iter_mut()
            .zip(src.iter())
            .for_each(|(d, s)| *d += s.clone());

        Ok(())
    }

    /// Zero-copy slice-based scaling
    pub fn scale_slice_inplace<T: RealField>(slice: &mut [T], scalar: T) {
        slice.iter_mut().for_each(|x| *x *= scalar.clone());
    }

    /// Zero-copy windowed operations on slices
    pub fn windowed_operation<T, F, R>(slice: &[T], window_size: usize, op: F) -> Vec<R>
    where
        T: Clone,
        F: Fn(&[T]) -> R,
    {
        slice.windows(window_size)
            .map(op)
            .collect()
    }

    /// Zero-copy chunked operations on slices
    pub fn chunked_operation<T, F, R>(slice: &[T], chunk_size: usize, op: F) -> Vec<R>
    where
        T: Clone,
        F: Fn(&[T]) -> R,
    {
        slice.chunks(chunk_size)
            .map(op)
            .collect()
    }
}

/// Composable iterator chains for CFD operations
pub trait CfdIteratorChain<T>: Iterator<Item = T>
where
    T: RealField,
    Self: Sized,
{
    /// Chain with gradient computation
    fn with_gradient(self, spacing: T) -> impl Iterator<Item = T> {
        self.collect::<Vec<_>>()
            .windows(2)
            .map(move |w| (w[1].clone() - w[0].clone()) / spacing.clone())
            .collect::<Vec<_>>()
            .into_iter()
    }

    /// Chain with smoothing filter
    fn with_smoothing(self, window: usize) -> impl Iterator<Item = T> {
        let data: Vec<_> = self.collect();
        data.windows(window)
            .map(|w| {
                let sum = w.iter().fold(T::zero(), |acc, x| acc + x.clone());
                sum / T::from_usize(w.len()).unwrap()
            })
            .collect::<Vec<_>>()
            .into_iter()
    }
}

impl<I, T> CfdIteratorChain<T> for I
where
    I: Iterator<Item = T>,
    T: RealField,
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DVector;

    #[test]
    fn test_l2_norm() {
        let data = vec![3.0, 4.0];
        let norm = data.into_iter().l2_norm();
        assert!((norm - 5.0_f64).abs() < 1e-10);
    }

    #[test]
    fn test_vector_ops() {
        let a = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let b = DVector::from_vec(vec![4.0, 5.0, 6.0]);
        
        let sum = VectorOps::add(&a, &b).unwrap();
        assert_eq!(sum[0], 5.0);
        assert_eq!(sum[1], 7.0);
        assert_eq!(sum[2], 9.0);
        
        let dot = VectorOps::dot(&a, &b).unwrap();
        assert_eq!(dot, 32.0); // 1*4 + 2*5 + 3*6 = 32
    }

    #[test]
    fn test_windowed_diff() {
        let data = vec![1.0, 2.0, 4.0, 7.0, 11.0];
        let diffs: Vec<_> = data.into_iter()
            .windowed_diff(2, |w| w[1] - w[0])
            .collect();
        
        assert_eq!(diffs, vec![1.0, 2.0, 3.0, 4.0]);
    }
}
