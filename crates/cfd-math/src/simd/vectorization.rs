//! Vectorization utilities with SIMD optimizations for CFD computations.
//!
//! This module provides vectorized operations using moirai's data-parallel
//! primitives. Architecture-specific SIMD dispatch is handled internally by
//! hermes-simd (via `SimdOps`), but the generic `T: RealField` operations
//! here use moirai's `par_mut()` which already routes small arrays to the
//! sequential path automatically.
//!
//! # Soundness
//!
//! All operations use safe Rust. The prior `TypeId` + `transmute_copy` patterns
//! have been eliminated — generic operations delegate to moirai's parallel
//! iteration, and concrete f32/f64 operations delegate to `SimdOps` (hermes).

use crate::error::Result;
use crate::simd::SimdOps;
use moirai::prelude::ParallelSliceMut;
use moirai::{reduce_index_with, Adaptive};
use nalgebra::RealField;

/// Vectorized operations for CFD computations with SIMD support
pub struct VectorizedOps {
    simd_ops: SimdOps,
}

impl VectorizedOps {
    /// Create a new vectorized operations handler
    pub fn new() -> Self {
        Self {
            simd_ops: SimdOps::new(),
        }
    }

    /// Vectorized element-wise addition.
    ///
    /// For f32/f64, delegates to hermes SIMD dispatch. For other types,
    /// uses moirai's data-parallel iteration.
    pub fn add_vectorized_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        self.simd_ops.add(a, b, result)
    }

    /// Vectorized element-wise addition (f64).
    pub fn add_vectorized_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        self.simd_ops.add_f64(a, b, result)
    }

    /// Generic element-wise addition using moirai parallel iteration.
    pub fn add_generic<T: RealField + Copy + Send + Sync>(
        &self,
        a: &[T],
        b: &[T],
        result: &mut [T],
    ) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string()));
        }
        result.par_mut().enumerate(|i, r| *r = a[i] + b[i]);
        Ok(())
    }

    /// Vectorized element-wise multiplication (f32).
    pub fn mul_vectorized_f32(&self, a: &[f32], b: &[f32], result: &mut [f32]) -> Result<()> {
        self.simd_ops.mul(a, b, result)
    }

    /// Vectorized element-wise multiplication (f64).
    pub fn mul_vectorized_f64(&self, a: &[f64], b: &[f64], result: &mut [f64]) -> Result<()> {
        self.simd_ops.mul_f64(a, b, result)
    }

    /// Generic element-wise multiplication using moirai parallel iteration.
    pub fn mul_generic<T: RealField + Copy + Send + Sync>(
        &self,
        a: &[T],
        b: &[T],
        result: &mut [T],
    ) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string()));
        }
        result.par_mut().enumerate(|i, r| *r = a[i] * b[i]);
        Ok(())
    }

    /// Vectorized scalar multiplication (f32).
    pub fn scale_vectorized_f32(&self, input: &[f32], scalar: f32, result: &mut [f32]) -> Result<()> {
        self.simd_ops.scale(input, scalar, result)
    }

    /// Vectorized scalar multiplication (f64).
    pub fn scale_vectorized_f64(&self, input: &[f64], scalar: f64, result: &mut [f64]) -> Result<()> {
        self.simd_ops.scale_f64(input, scalar, result)
    }

    /// Generic scalar multiplication using moirai parallel iteration.
    pub fn scale_generic<T: RealField + Copy + Send + Sync>(
        &self,
        input: &[T],
        scalar: T,
        result: &mut [T],
    ) -> Result<()> {
        if input.len() != result.len() {
            return Err(cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string()));
        }
        result.par_mut().enumerate(|i, r| *r = scalar * input[i]);
        Ok(())
    }

    /// Compute dot product (f32) via hermes SIMD.
    pub fn dot_product_f32(&self, a: &[f32], b: &[f32]) -> Result<f32> {
        self.simd_ops.dot(a, b)
    }

    /// Compute dot product (f64) via hermes SIMD.
    pub fn dot_product_f64(&self, a: &[f64], b: &[f64]) -> Result<f64> {
        self.simd_ops.dot_f64(a, b)
    }

    /// Generic dot product using moirai parallel reduction.
    pub fn dot_product_generic<T: RealField + Copy + Send + Sync>(
        &self,
        a: &[T],
        b: &[T],
    ) -> Result<T> {
        if a.len() != b.len() {
            return Err(cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string()));
        }
        let sum =
            reduce_index_with::<Adaptive, _, _, _>(a.len(), T::zero(), |i| a[i] * b[i], |acc, v| acc + v);
        Ok(sum)
    }

    /// Compute L2 norm using moirai parallel reduction.
    pub fn l2_norm<T: RealField + Copy + Send + Sync>(
        &self,
        input: &[T],
    ) -> Result<T> {
        let dot = self.dot_product_generic(input, input)?;
        Ok(dot.sqrt())
    }

    /// Broadcasting addition: adds scalar to each element of vector.
    pub fn broadcast_add<T: RealField + Copy + Send + Sync>(
        &self,
        input: &[T],
        scalar: T,
        result: &mut [T],
    ) -> Result<()> {
        if input.len() != result.len() {
            return Err(cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string()));
        }
        result.par_mut().enumerate(|i, r| *r = input[i] + scalar);
        Ok(())
    }

    /// Broadcasting multiplication: multiplies each element by broadcasted vector (row-wise).
    pub fn broadcast_mul_vector<T: RealField + Copy + Send + Sync>(
        &self,
        matrix: &[T],
        matrix_cols: usize,
        vector: &[T],
        result: &mut [T],
    ) -> Result<()> {
        if !matrix.len().is_multiple_of(matrix_cols) {
            return Err(cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string()));
        }
        if vector.len() != matrix_cols {
            return Err(cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string()));
        }
        if result.len() != matrix.len() {
            return Err(cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string()));
        }
        result
            .par_mut()
            .enumerate(|idx, r| *r = matrix[idx] * vector[idx % matrix_cols]);
        Ok(())
    }

    /// Vectorized matrix-vector multiplication for CSR sparse matrices.
    pub fn matvec_csr<T: RealField + Copy + Send + Sync>(
        &self,
        values: &[T],
        col_indices: &[usize],
        row_ptr: &[usize],
        x: &[T],
        y: &mut [T],
    ) -> Result<()> {
        if row_ptr.len() != y.len() + 1 {
            return Err(cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string()));
        }
        y.par_mut().enumerate(|i, y_i| {
            let start = row_ptr[i];
            let end = row_ptr[i + 1];
            *y_i = (start..end)
                .map(|j| values[j] * x[col_indices[j]])
                .fold(T::zero(), |acc, val| acc + val);
        });
        Ok(())
    }

    /// Vectorized convolution for filtering operations.
    pub fn convolution<T: RealField + Copy + Send + Sync>(
        &self,
        signal: &[T],
        kernel: &[T],
        result: &mut [T],
    ) -> Result<()> {
        let signal_len = signal.len();
        let kernel_len = kernel.len();
        let result_len = signal_len + kernel_len - 1;

        if result.len() != result_len {
            return Err(cfd_core::error::Error::InvalidInput("Dimension mismatch".to_string()));
        }

        result.par_mut().enumerate(|n, output| {
            let k_min = n.saturating_sub(signal_len - 1);
            let k_max = n.min(kernel_len - 1);
            let mut sum = T::zero();
            if k_min <= k_max {
                for k in k_min..=k_max {
                    sum += signal[n - k] * kernel[k];
                }
            }
            *output = sum;
        });
        Ok(())
    }
}

impl Default for VectorizedOps {
    fn default() -> Self {
        Self::new()
    }
}

/// Stencil operations with SIMD optimization
pub struct StencilOps {
    _simd_ops: SimdOps,
}

impl StencilOps {
    /// Create a new stencil operations handler
    pub fn new() -> Self {
        Self {
            _simd_ops: SimdOps::new(),
        }
    }

    /// Apply 3-point stencil with SIMD optimization
    pub fn apply_stencil_3point<T: RealField + Copy + Send + Sync>(
        &self,
        input: &[T],
        coeffs: &[T; 3],
        result: &mut [T],
    ) -> Result<()> {
        if input.len() < 3 || result.len() != input.len() - 2 {
            return Err(cfd_core::error::Error::InvalidInput(
                "Invalid dimensions for 3-point stencil".to_string(),
            ));
        }
        result.par_mut().enumerate(|i, r| {
            let idx = i + 1;
            *r = coeffs[0] * input[idx - 1] + coeffs[1] * input[idx] + coeffs[2] * input[idx + 1];
        });
        Ok(())
    }

    /// Apply 5-point stencil with SIMD optimization
    pub fn apply_stencil_5point<T: RealField + Copy + Send + Sync>(
        &self,
        input: &[T],
        coeffs: &[T; 5],
        result: &mut [T],
    ) -> Result<()> {
        if input.len() < 5 || result.len() != input.len() - 4 {
            return Err(cfd_core::error::Error::InvalidInput(
                "Invalid dimensions for 5-point stencil".to_string(),
            ));
        }
        result.par_mut().enumerate(|i, r| {
            let idx = i + 2;
            *r = coeffs[0] * input[idx - 2]
                + coeffs[1] * input[idx - 1]
                + coeffs[2] * input[idx]
                + coeffs[3] * input[idx + 1]
                + coeffs[4] * input[idx + 2];
        });
        Ok(())
    }
}

impl Default for StencilOps {
    fn default() -> Self {
        Self::new()
    }
}
