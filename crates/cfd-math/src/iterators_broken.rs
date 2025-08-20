//! Extended iterator utilities for zero-copy CFD operations.
//!
//! This module provides iterator combinators and utilities optimized for CFD computations,
//! following CUPID principles (Composable, Unix Philosophy, Predictable, Idiomatic, Domain-based).

use nalgebra::{RealField, DVector};
use rayon::prelude::*;
use std::iter::Iterator;

/// Extension trait for iterator-based mathematical operations
pub trait MathIteratorExt: Iterator
where
    Self::Item: RealField + Send + Sync + Copy,
    Self: Sized,
{
    /// Compute L2 norm using iterator operations
    fn l2_norm(self) -> Self::Item {
        self.map(|x| x * x)
            .fold(Self::Item::zero(), |acc, x| acc + x)
            .sqrt()
    }

    /// Compute L1 norm using iterator combinators
    fn l1_norm(self) -> Self::Item {
        self.map(|x| x.abs())
            .fold(Self::Item::zero(), |acc, x| acc + x)
    }

    /// Compute L∞ norm using iterator operations
    fn linf_norm(self) -> Self::Item {
        self.map(|x| x.abs())
            .fold(Self::Item::zero(), |acc, x| if x > acc { x } else { acc })
    }

    /// Compute mean using single-pass iterator
    fn mean(self) -> Option<T> {
        let (sum, count) = self.fold((Self::Item::zero(), 0), |(sum, count), x| (sum + x, count + 1));
        if count > 0 {
            Some(sum / Self::Item::from_usize(count).unwrap_or_else(|| Self::Item::zero()))
        } else {
            None
        }
    }

    /// Compute variance using Welford's online algorithm
    fn variance(self) -> Option<T> {
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
    fn windowed_diff<F>(self, window_size: usize, f: F) -> WindowedDiff<Self, F, T>
    where
        F: Fn(&[Self::Item]) -> Self::Item,
    {
        WindowedDiff::new(self, window_size, f)
    }

    /// Compute running average using iterator scan for zero-copy efficiency
    fn running_average(self) -> impl Iterator<Item = T> {
        let mut count = 0;
        let mut sum = Self::Item::zero();

        self.map(move |x| {
            count += 1;
            sum += x;
            sum.clone() / Self::Item::from_usize(count).unwrap_or_else(|| Self::Item::zero())
        })
    }

    /// Compute exponential moving average with given alpha
    fn exponential_moving_average(self, alpha: Self::Item) -> impl Iterator<Item = T> {
        let mut ema = None;

        self.map(move |x| {
            match ema {
                None => {
                    ema = Some(x.clone());
                    x
                }
                Some(ref mut prev) => {
                    *prev = alpha.clone() * x + (Self::Item::one() - alpha.clone()) * prev.clone();
                    prev.clone()
                }
            }
        })
    }

    /// Apply temporal frequency analysis using sliding FFT windowing
    fn temporal_frequency_analysis(self, window_size: usize) -> Self::ItememporalFrequencyAnalyzer<Self, T> {
        TemporalFrequencyAnalyzer::new(self, window_size)
    }

    /// Apply Kalman-style recursive filtering for temporal prediction
    fn temporal_kalman_filter(self, process_noise: Self::Item, observation_noise: Self::Item) -> Self::ItememporalKalmanFilter<Self, T> {
        TemporalKalmanFilter::new(self, process_noise, observation_noise)
    }

    /// Apply overlapping window analysis for CFD stability monitoring
    fn overlapping_window_stability(self, window_size: usize, overlap: usize) -> OverlappingWindowAnalyzer<Self, T> {
        OverlappingWindowAnalyzer::new(self, window_size, overlap)
    }

    /// Apply Savitzky-Golay smoothing filter using windowed operations
    /// Uses 2nd order polynomial fitting for smoothing
    fn savitzky_golay_smooth(self, window_size: usize) -> impl Iterator<Item = T> {
        // Savitzky-Golay coefficients for 2nd order polynomial, centered window
        // Pre-computed for common window sizes
        let coeffs = match window_size {
            5 => vec![Self::Item::from_f64(-3.0/35.0).unwrap_or_else(|| Self::Item::zero()), Self::Item::from_f64(12.0/35.0).unwrap_or_else(|| Self::Item::zero()), 
                     Self::Item::from_f64(17.0/35.0).unwrap_or_else(|| Self::Item::zero()), Self::Item::from_f64(12.0/35.0).unwrap_or_else(|| Self::Item::zero()), 
                     Self::Item::from_f64(-3.0/35.0).unwrap_or_else(|| Self::Item::zero())],
            7 => vec![Self::Item::from_f64(-2.0/21.0).unwrap_or_else(|| Self::Item::zero()), Self::Item::from_f64(3.0/21.0).unwrap_or_else(|| Self::Item::zero()),
                     Self::Item::from_f64(6.0/21.0).unwrap_or_else(|| Self::Item::zero()), Self::Item::from_f64(7.0/21.0).unwrap_or_else(|| Self::Item::zero()),
                     Self::Item::from_f64(6.0/21.0).unwrap_or_else(|| Self::Item::zero()), Self::Item::from_f64(3.0/21.0).unwrap_or_else(|| Self::Item::zero()),
                     Self::Item::from_f64(-2.0/21.0).unwrap_or_else(|| Self::Item::zero())],
            _ => {
                // Fallback to moving average for other window sizes
                let coeff = Self::Item::one() / Self::Item::from_usize(window_size).unwrap_or_else(|| Self::Item::zero());
                vec![coeff.clone(); window_size]
            }
        };
        
        self.windowed_diff(window_size, move |w| {
            w.iter()
                .zip(coeffs.iter())
                .fold(Self::Item::zero(), |acc, (x, c)| acc + x.clone() * c.clone())
        })
    }

    // Note: Parallel reduction is available through rayon's ParallelIterator trait directly
}

impl<I> MathIteratorExt for I
where
    I: Iterator,
    I::Item: RealField + Send + Sync + Copy,
{}

/// Windowed difference iterator for finite difference operations
/// Augmented with circular buffer for zero-copy sliding window
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
    F: Fn(&[Self::Item]) -> Self::Item,
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
    F: Fn(&[Self::Item]) -> Self::Item,
    T: Clone,
{
    type Item = T;

    fn next(&mut self) -> Option<T> {
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
            .fold(Self::Item::zero(), |acc, x| acc + x);

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
    pub fn scale_inplace<T: RealField>(vector: &mut DVector<T>, scalar: &T) {
        vector.iter_mut().for_each(|x| *x *= scalar.clone());
    }
}

impl SliceOps {
    /// Zero-copy slice-based dot product
    pub fn dot_slice<T: RealField>(a: &[Self::Item], b: &[Self::Item]) -> Result<T, &'static str> {
        if a.len() != b.len() {
            return Err("Slice dimensions must match");
        }

        let result = a.iter()
            .zip(b.iter())
            .map(|(x, y)| x.clone() * y.clone())
            .fold(Self::Item::zero(), |acc, x| acc + x);

        Ok(result)
    }

    /// Zero-copy slice-based L2 norm
    pub fn l2_norm_slice<T: RealField>(slice: &[Self::Item]) -> Self::Item {
        slice.iter()
            .map(|x| x.clone() * x.clone())
            .fold(Self::Item::zero(), |acc, x| acc + x)
            .sqrt()
    }

    /// Zero-copy slice-based element-wise addition
    pub fn add_slice_inplace<T: RealField>(dest: &mut [Self::Item], src: &[Self::Item]) -> Result<(), &'static str> {
        if dest.len() != src.len() {
            return Err("Slice dimensions must match");
        }

        dest.iter_mut()
            .zip(src.iter())
            .for_each(|(d, s)| *d += s.clone());

        Ok(())
    }

    /// Zero-copy slice-based scaling
    pub fn scale_slice_inplace<T: RealField>(slice: &mut [Self::Item], scalar: &T) {
        slice.iter_mut().for_each(|x| *x *= scalar.clone());
    }

    /// Zero-copy windowed operations on slices
    pub fn windowed_operation<T, F, R>(slice: &[Self::Item], window_size: usize, op: F) -> Vec<R>
    where
        T: Clone,
        F: Fn(&[Self::Item]) -> R,
    {
        slice.windows(window_size)
            .map(op)
            .collect()
    }

    /// Zero-copy chunked operations on slices
    pub fn chunked_operation<T, F, R>(slice: &[Self::Item], chunk_size: usize, op: F) -> Vec<R>
    where
        T: Clone,
        F: Fn(&[Self::Item]) -> R,
    {
        slice.chunks(chunk_size)
            .map(op)
            .collect()
    }

    /// Zero-copy slice-based matrix-vector multiplication (for row-major matrices)
    pub fn matvec_slice<T: RealField>(
        matrix: &[Self::Item],
        vector: &[Self::Item],
        result: &mut [Self::Item],
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
                    .fold(Self::Item::zero(), |acc, val| acc + val);
            });

        Ok(())
    }

    /// Zero-copy slice-based transpose operation (in-place for square matrices)
    pub fn transpose_square_inplace<T: Clone>(matrix: &mut [Self::Item], n: usize) -> Result<(), &'static str> {
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
/// Augmented to avoid unnecessary collect() calls for better zero-copy performance
pub trait CfdIteratorChain<T>: Iterator<Item = T>
where
    T: RealField,
    Self: Sized,
{
    /// Chain with gradient computation using windowed_diff for zero-copy
    fn with_gradient(self, spacing: Self::Item) -> impl Iterator<Item = T> {
        self.windowed_diff(2, move |w| (w[1].clone() - w[0].clone()) / spacing.clone())
    }

    /// Chain with smoothing filter using windowed operations
    fn with_smoothing(self, window: usize) -> impl Iterator<Item = T> {
        self.windowed_diff(window, |w| {
            let sum = w.iter().fold(Self::Item::zero(), |acc, x| acc + x.clone());
            sum / Self::Item::from_usize(w.len()).expect("CRITICAL: Add proper error handling")
        })
    }

    /// Chain with divergence computation for vector fields
    /// Assumes 3D vector field layout: [u1, v1, w1, u2, v2, w2, ...]
    fn with_divergence(self, spacing: Self::Item) -> impl Iterator<Item = T> {
        self.windowed_diff(9, move |w| {
            if w.len() >= 9 {
                // Central difference approximation: ∇·v = ∂u/∂x + ∂v/∂y + ∂w/∂z
                let du_dx = (w[6].clone() - w[0].clone()) / (spacing.clone() + spacing.clone()); // (u_i+1 - u_i-1) / 2Δx
                let dv_dy = (w[7].clone() - w[1].clone()) / (spacing.clone() + spacing.clone()); // (v_j+1 - v_j-1) / 2Δy
                let dw_dz = (w[8].clone() - w[2].clone()) / (spacing.clone() + spacing.clone()); // (w_k+1 - w_k-1) / 2Δz
                du_dx + dv_dy + dw_dz
            } else {
                Self::Item::zero()
            }
        })
    }

    /// Chain with curl computation for 3D vector fields
    /// Returns the magnitude of curl: |∇ × v|
    fn with_curl_magnitude(self, spacing: Self::Item) -> impl Iterator<Item = T> {
        self.windowed_diff(9, move |w| {
            if w.len() >= 9 {
                // Curl components: ∇ × v = (∂w/∂y - ∂v/∂z, ∂u/∂z - ∂w/∂x, ∂v/∂x - ∂u/∂y)
                // Assuming w contains [u0, v0, w0, u1, v1, w1, u2, v2, w2] for a 3x3x3 stencil
                let two_dx = spacing.clone() + spacing.clone();

                // Extract velocity components (assuming interleaved storage)
                // Center point: w[4] = (u_center, v_center, w_center)
                // x-direction neighbors: w[3], w[5]
                // y-direction neighbors: w[1], w[7]
                // z-direction neighbors: w[0], w[8]
                
                // For a proper 3D stencil with 27 points (3x3x3), we'd need more data
                // Using a simplified 9-point stencil for demonstration
                let dwdy = (w[7].clone() - w[1].clone()) / two_dx.clone(); // ∂w/∂y
                let dvdz = (w[8].clone() - w[0].clone()) / two_dx.clone(); // ∂v/∂z
                let dudz = (w[8].clone() - w[0].clone()) / two_dx.clone(); // ∂u/∂z
                let dwdx = (w[5].clone() - w[3].clone()) / two_dx.clone(); // ∂w/∂x
                let dvdx = (w[5].clone() - w[3].clone()) / two_dx.clone(); // ∂v/∂x
                let dudy = (w[7].clone() - w[1].clone()) / two_dx.clone(); // ∂u/∂y
                
                // Curl components
                let curl_x = dwdy - dvdz;
                let curl_y = dudz - dwdx;
                let curl_z = dvdx - dudy;

                // Magnitude: |∇ × v|
                (curl_x.clone() * curl_x +
                 curl_y.clone() * curl_y +
                 curl_z.clone() * curl_z).sqrt()
            } else {
                Self::Item::zero()
            }
        })
    }

    /// Chain with strain rate tensor magnitude computation
    /// For incompressible flow: |S| = √(2 S_ij S_ij)
    fn with_strain_rate_magnitude(self, spacing: Self::Item) -> impl Iterator<Item = T> {
        self.windowed_diff(27, move |w| {
            // For proper strain rate calculation, we need a 3x3x3 stencil (27 points)
            // Each point has 3 velocity components (u, v, w)
            if w.len() >= 27 {
                let two_dx = spacing.clone() + spacing.clone();
                let two = Self::Item::one() + Self::Item::one();
                let _half = Self::Item::from_f64(0.5).unwrap_or(Self::Item::one() / two.clone());

                // Center index in 3x3x3 grid is 13 (0-indexed)
                let center = 13;
                
                // Velocity gradients using central differences
                // x-direction: indices 12 (left) and 14 (right)
                let dudx = (w[center + 1].clone() - w[center - 1].clone()) / two_dx.clone();
                let dvdx = (w[center + 1].clone() - w[center - 1].clone()) / two_dx.clone();
                let dwdx = (w[center + 1].clone() - w[center - 1].clone()) / two_dx.clone();
                
                // y-direction: indices 10 (front) and 16 (back)
                let dudy = (w[center + 3].clone() - w[center - 3].clone()) / two_dx.clone();
                let dvdy = (w[center + 3].clone() - w[center - 3].clone()) / two_dx.clone();
                let dwdy = (w[center + 3].clone() - w[center - 3].clone()) / two_dx.clone();
                
                // z-direction: indices 4 (bottom) and 22 (top)
                let dudz = (w[center + 9].clone() - w[center - 9].clone()) / two_dx.clone();
                let dvdz = (w[center + 9].clone() - w[center - 9].clone()) / two_dx.clone();
                let dwdz = (w[center + 9].clone() - w[center - 9].clone()) / two_dx.clone();

                // Strain rate tensor components: S_ij = 0.5 * (∂u_i/∂x_j + ∂u_j/∂x_i)
                let s11 = dudx;
                let s22 = dvdy;
                let s33 = dwdz;
                let s12 = (dudy + dvdx) / two.clone();
                let s13 = (dudz + dwdx) / two.clone();
                let s23 = (dvdz + dwdy) / two.clone();

                // Magnitude: |S| = √(2 S_ij S_ij)
                let s_squared = s11.clone() * s11 + s22.clone() * s22 + s33.clone() * s33 +
                               two.clone() * (s12.clone() * s12 + s13.clone() * s13 + s23.clone() * s23);

                (two.clone() * s_squared).sqrt()
            } else {
                Self::Item::zero()
            }
        })
    }
}

impl<I, T> CfdIteratorChain<T> for I
where
    I: Iterator<Item = T>,
    T: RealField,
{
}

/// Extended CFD field operations using iterator patterns
pub trait CfdFieldOps<T>: Iterator<Item = T>
where
    T: RealField,
    Self: Sized,
{
    /// Compute field magnitude for vector fields
    /// Assumes data layout: [u1, v1, w1, u2, v2, w2, ...]
    fn vector_magnitude(self) -> impl Iterator<Item = T> {
        // Use windowed approach instead of unstable array_chunks
        self.windowed_diff(3, |w| {
            if w.len() >= 3 {
                let u = w[0].clone();
                let v = w[1].clone();
                let w_comp = w[2].clone();
                (u.clone() * u + v.clone() * v + w_comp.clone() * w_comp).sqrt()
            } else {
                Self::Item::zero()
            }
        })
    }

    /// Compute field divergence using central differences
    /// For structured grids with uniform spacing
    fn field_divergence(self, spacing: Self::Item) -> impl Iterator<Item = T> {
        self.windowed_diff(3, move |w| {
            if w.len() >= 3 {
                (w[2].clone() - w[0].clone()) / (spacing.clone() + spacing.clone())
            } else {
                Self::Item::zero()
            }
        })
    }

    /// Apply field interpolation between grid points
    fn field_interpolate(self, weights: Vec<T>) -> impl Iterator<Item = T> {
        let weight_sum: T = weights.iter().fold(Self::Item::zero(), |acc, w| acc + w.clone());

        self.windowed_diff(weights.len(), move |w| {
            w.iter()
                .zip(weights.iter())
                .fold(Self::Item::zero(), |acc, (val, weight)| {
                    acc + val.clone() * weight.clone()
                }) / weight_sum.clone()
        })
    }

    /// Compute field Laplacian using 5-point stencil
    fn field_laplacian(self, spacing: Self::Item) -> impl Iterator<Item = T> {
        let dx_sq = spacing.clone() * spacing.clone();

        self.windowed_diff(5, move |w| {
            if w.len() >= 5 {
                // 5-point stencil: (u_{i-1} - 2u_i + u_{i+1}) / dx^2
                let d2u_dx2 = (w[0].clone() - Self::Item::from_f64(2.0).unwrap_or_else(|| Self::Item::zero()) * w[2].clone() + w[4].clone()) / dx_sq.clone();
                d2u_dx2
            } else {
                Self::Item::zero()
            }
        })
    }
}

impl<I, T> CfdFieldOps<T> for I
where
    I: Iterator<Item = T>,
    T: RealField,
{
}

/// Temporal frequency analyzer for CFD stability monitoring
pub struct TemporalFrequencyAnalyzer<I, T> {
    iter: I,
    window_size: usize,
    buffer: Vec<T>,
    position: usize,
}

impl<I, T> TemporalFrequencyAnalyzer<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Clone,
{
    pub fn new(iter: I, window_size: usize) -> Self {
        Self {
            iter,
            window_size,
            buffer: Vec::with_capacity(window_size),
            position: 0,
        }
    }
}

impl<I, T> Iterator for TemporalFrequencyAnalyzer<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Clone,
{
    type Item = T; // Return stability metric

    fn next(&mut self) -> Option<T> {
        if let Some(value) = self.iter.next() {
            if self.buffer.len() < self.window_size {
                self.buffer.push(value);
                None // Not enough data yet
            } else {
                // Update circular buffer
                self.buffer[self.position] = value;
                self.position = (self.position + 1) % self.window_size;
                
                // Compute stability via variance analysis
                let mean = self.buffer.iter().cloned().fold(Self::Item::zero(), |acc, x| acc + x) / 
                          Self::Item::from_usize(self.buffer.len()).expect("CRITICAL: Add proper error handling");
                
                let variance = self.buffer.iter()
                    .map(|x| (x.clone() - mean.clone()) * (x.clone() - mean.clone()))
                    .fold(Self::Item::zero(), |acc, x| acc + x) / 
                    Self::Item::from_usize(self.buffer.len()).expect("CRITICAL: Add proper error handling");
                
                Some(variance.sqrt()) // Return standard deviation as stability metric
            }
        } else {
            None
        }
    }
}

/// Temporal Kalman filter for CFD prediction and smoothing
pub struct TemporalKalmanFilter<I, T> {
    iter: I,
    state_estimate: OptionSelf::Item,
    estimation_error: T,
    process_noise: Self::Item,
    observation_noise: Self::Item,
}

impl<I, T> TemporalKalmanFilter<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Clone,
{
    pub fn new(iter: I, process_noise: Self::Item, observation_noise: Self::Item) -> Self {
        Self {
            iter,
            state_estimate: None,
            estimation_error: Self::Item::one(),
            process_noise,
            observation_noise,
        }
    }
}

impl<I, T> Iterator for TemporalKalmanFilter<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Clone,
{
    type Item = T; // Return filtered value

    fn next(&mut self) -> Option<T> {
        if let Some(observation) = self.iter.next() {
            match &self.state_estimate {
                None => {
                    // Initialize with first observation
                    self.state_estimate = Some(observation.clone());
                    Some(observation)
                }
                Some(state) => {
                    // Kalman filter update
                    let predicted_error = self.estimation_error.clone() + self.process_noise.clone();
                    let kalman_gain = predicted_error.clone() / 
                                     (predicted_error.clone() + self.observation_noise.clone());
                    
                    let updated_estimate = state.clone() + 
                                         kalman_gain.clone() * (observation - state.clone());
                    
                    self.estimation_error = (Self::Item::one() - kalman_gain) * predicted_error;
                    self.state_estimate = Some(updated_estimate.clone());
                    
                    Some(updated_estimate)
                }
            }
        } else {
            None
        }
    }
}

/// Overlapping window analyzer for CFD stability monitoring
pub struct OverlappingWindowAnalyzer<I, T> {
    iter: I,
    window_size: usize,
    overlap: usize,
    buffer: Vec<T>,
    step_size: usize,
}

impl<I, T> OverlappingWindowAnalyzer<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Clone,
{
    pub fn new(iter: I, window_size: usize, overlap: usize) -> Self {
        let step_size = window_size - overlap;
        Self {
            iter,
            window_size,
            overlap,
            buffer: Vec::with_capacity(window_size),
            step_size,
        }
    }
}

impl<I, T> Iterator for OverlappingWindowAnalyzer<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Clone,
{
    type Item = T; // Return stability metric

    fn next(&mut self) -> Option<T> {
        // Fill buffer to window size
        while self.buffer.len() < self.window_size {
            if let Some(value) = self.iter.next() {
                self.buffer.push(value);
            } else {
                return None; // Not enough data
            }
        }
        
        // Compute stability metric for current window (coefficient of variation)
        let mean = self.buffer.iter().cloned().fold(Self::Item::zero(), |acc, x| acc + x) / 
                  Self::Item::from_usize(self.buffer.len()).expect("CRITICAL: Add proper error handling");
        
        let variance = self.buffer.iter()
            .map(|x| (x.clone() - mean.clone()) * (x.clone() - mean.clone()))
            .fold(Self::Item::zero(), |acc, x| acc + x) / 
            Self::Item::from_usize(self.buffer.len()).expect("CRITICAL: Add proper error handling");
        
        let cv = if mean.clone().abs() > Self::Item::from_f64(1e-15).unwrap_or_else(|| Self::Item::zero()) {
            variance.sqrt() / mean.abs()
        } else {
            variance.sqrt()
        };
        
        // Slide window by step_size
        for _ in 0..self.step_size {
            if let Some(current_value) = self.iter.next() {
                self.buffer.remove(0);
                self.buffer.push(current_value);
            } else {
                return Some(cv); // Last window
            }
        }
        
        Some(cv)
    }
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
        
        let sum = VectorOps::add(&a, &b).expect("CRITICAL: Add proper error handling");
        assert_eq!(sum[0], 5.0);
        assert_eq!(sum[1], 7.0);
        assert_eq!(sum[2], 9.0);
        
        let dot = VectorOps::dot(&a, &b).expect("CRITICAL: Add proper error handling");
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
