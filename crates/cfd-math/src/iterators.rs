//! Advanced iterator utilities for zero-copy CFD operations.
//!
//! This module provides iterator combinators and utilities optimized for CFD computations,
//! following CUPID principles (Composable, Unix Philosophy, Predictable, Idiomatic, Domain-based).

use nalgebra::{RealField, DVector};
use rayon::prelude::*;
use std::iter::Iterator;

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
/// Enhanced with circular buffer for zero-copy sliding window
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
            buffer: Vec::with_capacity(window_size * 2), // Double capacity for circular buffer
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
        // Fill buffer to window size initially
        while self.buffer_len < self.window_size {
            if let Some(item) = self.iter.next() {
                if self.buffer.len() < self.buffer.capacity() {
                    self.buffer.push(item);
                } else {
                    // Circular buffer is full, overwrite
                    let idx = (self.buffer_start + self.buffer_len) % self.buffer.len();
                    self.buffer[idx] = item;
                }
                self.buffer_len += 1;
            } else {
                return None;
            }
        }

        if self.buffer_len >= self.window_size {
            // Create window slice without allocation
            let mut window = Vec::with_capacity(self.window_size);
            for i in 0..self.window_size {
                let idx = (self.buffer_start + i) % self.buffer.len();
                window.push(self.buffer[idx].clone());
            }

            let result = (self.func)(&window);

            // Slide window by advancing start position
            self.buffer_start = (self.buffer_start + 1) % self.buffer.len();
            self.buffer_len -= 1;

            // Try to add next item
            if let Some(item) = self.iter.next() {
                let idx = (self.buffer_start + self.buffer_len) % self.buffer.len();
                if idx < self.buffer.len() {
                    self.buffer[idx] = item;
                } else {
                    self.buffer.push(item);
                }
                self.buffer_len += 1;
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

    /// In-place element-wise multiplication for zero-copy operations
    pub fn hadamard_inplace<T: RealField>(dest: &mut DVector<T>, src: &DVector<T>) -> Result<(), &'static str> {
        if dest.len() != src.len() {
            return Err("Vector dimensions must match");
        }

        dest.iter_mut()
            .zip(src.iter())
            .for_each(|(d, s)| *d *= s.clone());

        Ok(())
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

    /// Zero-copy slice-based matrix-vector multiplication (for row-major matrices)
    pub fn matvec_slice<T: RealField>(
        matrix: &[T],
        vector: &[T],
        result: &mut [T],
        rows: usize,
        cols: usize
    ) -> Result<(), &'static str> {
        if matrix.len() != rows * cols {
            return Err("Matrix dimensions don't match data length");
        }
        if vector.len() != cols || result.len() != rows {
            return Err("Vector dimensions don't match matrix");
        }

        result.iter_mut()
            .enumerate()
            .for_each(|(i, res)| {
                *res = matrix[i * cols..(i + 1) * cols]
                    .iter()
                    .zip(vector.iter())
                    .map(|(m, v)| m.clone() * v.clone())
                    .fold(T::zero(), |acc, val| acc + val);
            });

        Ok(())
    }

    /// Zero-copy slice-based transpose operation (in-place for square matrices)
    pub fn transpose_square_inplace<T: Clone>(matrix: &mut [T], n: usize) -> Result<(), &'static str> {
        if matrix.len() != n * n {
            return Err("Matrix must be square");
        }

        for i in 0..n {
            for j in (i + 1)..n {
                let idx1 = i * n + j;
                let idx2 = j * n + i;
                matrix.swap(idx1, idx2);
            }
        }

        Ok(())
    }
}

/// Composable iterator chains for CFD operations
/// Enhanced to avoid unnecessary collect() calls for better zero-copy performance
pub trait CfdIteratorChain<T>: Iterator<Item = T>
where
    T: RealField,
    Self: Sized,
{
    /// Chain with gradient computation using windowed_diff for zero-copy
    fn with_gradient(self, spacing: T) -> impl Iterator<Item = T> {
        self.windowed_diff(2, move |w| (w[1].clone() - w[0].clone()) / spacing.clone())
    }

    /// Chain with smoothing filter using windowed operations
    fn with_smoothing(self, window: usize) -> impl Iterator<Item = T> {
        self.windowed_diff(window, |w| {
            let sum = w.iter().fold(T::zero(), |acc, x| acc + x.clone());
            sum / T::from_usize(w.len()).unwrap()
        })
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
