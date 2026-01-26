//! Vectorization utilities with SIMD optimizations for CFD computations.
//!
//! This module provides vectorized operations that leverage architecture-specific
//! SIMD instructions with SWAR fallbacks for improved performance.

use crate::error::Result;
use crate::simd::{SimdOps, VectorOps as SimdVectorOps};
use nalgebra::RealField;
use num_traits::cast::ToPrimitive;
use rayon::prelude::*;

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

    /// Vectorized element-wise addition with SIMD optimization
    pub fn add_vectorized<T: RealField + Copy + Send + Sync + ToPrimitive>(
        &self,
        a: &[T],
        b: &[T],
        result: &mut [T],
    ) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        // Try to use SIMD for f32/f64
        if std::any::TypeId::of::<T>() == std::any::TypeId::of::<f32>() {
            // Safe transmute for f32
            let a_f32 = unsafe { std::slice::from_raw_parts(a.as_ptr().cast::<f32>(), a.len()) };
            let b_f32 = unsafe { std::slice::from_raw_parts(b.as_ptr().cast::<f32>(), b.len()) };
            let result_f32 = unsafe {
                std::slice::from_raw_parts_mut(result.as_mut_ptr().cast::<f32>(), result.len())
            };
            return self.simd_ops.add(a_f32, b_f32, result_f32);
        }

        if std::any::TypeId::of::<T>() == std::any::TypeId::of::<f64>() {
            let a_f64 = unsafe { std::slice::from_raw_parts(a.as_ptr().cast::<f64>(), a.len()) };
            let b_f64 = unsafe { std::slice::from_raw_parts(b.as_ptr().cast::<f64>(), b.len()) };
            let result_f64 = unsafe {
                std::slice::from_raw_parts_mut(result.as_mut_ptr().cast::<f64>(), result.len())
            };
            return self.simd_ops.add_f64(a_f64, b_f64, result_f64);
        }

        // Fallback to parallel iterator for other types
        result
            .par_iter_mut()
            .zip(a.par_iter().zip(b.par_iter()))
            .for_each(|(r, (a_val, b_val))| {
                *r = *a_val + *b_val;
            });

        Ok(())
    }

    /// Vectorized element-wise multiplication with SIMD
    pub fn mul_vectorized<T: RealField + Copy + Send + Sync + ToPrimitive>(
        &self,
        a: &[T],
        b: &[T],
        result: &mut [T],
    ) -> Result<()> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        // Try to use SIMD for f32/f64
        if std::any::TypeId::of::<T>() == std::any::TypeId::of::<f32>() {
            let a_f32 = unsafe { std::slice::from_raw_parts(a.as_ptr().cast::<f32>(), a.len()) };
            let b_f32 = unsafe { std::slice::from_raw_parts(b.as_ptr().cast::<f32>(), b.len()) };
            let result_f32 = unsafe {
                std::slice::from_raw_parts_mut(result.as_mut_ptr().cast::<f32>(), result.len())
            };
            return self.simd_ops.mul(a_f32, b_f32, result_f32);
        }

        if std::any::TypeId::of::<T>() == std::any::TypeId::of::<f64>() {
            let a_f64 = unsafe { std::slice::from_raw_parts(a.as_ptr().cast::<f64>(), a.len()) };
            let b_f64 = unsafe { std::slice::from_raw_parts(b.as_ptr().cast::<f64>(), b.len()) };
            let result_f64 = unsafe {
                std::slice::from_raw_parts_mut(result.as_mut_ptr().cast::<f64>(), result.len())
            };
            return self.simd_ops.mul_f64(a_f64, b_f64, result_f64);
        }

        // Fallback
        result
            .par_iter_mut()
            .zip(a.par_iter().zip(b.par_iter()))
            .for_each(|(r, (a_val, b_val))| {
                *r = *a_val * *b_val;
            });

        Ok(())
    }

    /// Vectorized scalar multiplication with broadcasting and SIMD
    pub fn scale_vectorized<T: RealField + Copy + Send + Sync + ToPrimitive>(
        &self,
        input: &[T],
        scalar: T,
        result: &mut [T],
    ) -> Result<()> {
        if input.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        // Try to use SIMD for f32/f64
        if std::any::TypeId::of::<T>() == std::any::TypeId::of::<f32>() {
            if let Some(scalar_f32) = scalar.to_f32() {
                let input_f32 = unsafe {
                    std::slice::from_raw_parts(input.as_ptr().cast::<f32>(), input.len())
                };
                let result_f32 = unsafe {
                    std::slice::from_raw_parts_mut(result.as_mut_ptr().cast::<f32>(), result.len())
                };
                return self.simd_ops.scale(input_f32, scalar_f32, result_f32);
            }
        }

        if std::any::TypeId::of::<T>() == std::any::TypeId::of::<f64>() {
            if let Some(scalar_f64) = scalar.to_f64() {
                let input_f64 = unsafe {
                    std::slice::from_raw_parts(input.as_ptr().cast::<f64>(), input.len())
                };
                let result_f64 = unsafe {
                    std::slice::from_raw_parts_mut(result.as_mut_ptr().cast::<f64>(), result.len())
                };
                return self.simd_ops.scale_f64(input_f64, scalar_f64, result_f64);
            }
        }

        // Fallback
        result
            .par_iter_mut()
            .zip(input.par_iter())
            .for_each(|(r, val)| {
                *r = scalar * *val;
            });

        Ok(())
    }

    /// Compute dot product with SIMD optimization
    pub fn dot_product<T: RealField + Copy + Send + Sync + ToPrimitive>(
        &self,
        a: &[T],
        b: &[T],
    ) -> Result<T> {
        if a.len() != b.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        // Try to use SIMD for f32
        if std::any::TypeId::of::<T>() == std::any::TypeId::of::<f32>() {
            let a_f32 = unsafe { std::slice::from_raw_parts(a.as_ptr().cast::<f32>(), a.len()) };
            let b_f32 = unsafe { std::slice::from_raw_parts(b.as_ptr().cast::<f32>(), b.len()) };
            let result = self.simd_ops.dot(a_f32, b_f32)?;
            // Safe because we checked the type
            return Ok(unsafe { std::mem::transmute_copy(&result) });
        }

        // Try to use SIMD for f64
        if std::any::TypeId::of::<T>() == std::any::TypeId::of::<f64>() {
            let a_f64 = unsafe { std::slice::from_raw_parts(a.as_ptr().cast::<f64>(), a.len()) };
            let b_f64 = unsafe { std::slice::from_raw_parts(b.as_ptr().cast::<f64>(), b.len()) };
            let result = self.simd_ops.dot_f64(a_f64, b_f64)?;
            return Ok(unsafe { std::mem::transmute_copy(&result) });
        }

        // Fallback for other types
        let sum = a
            .par_iter()
            .zip(b.par_iter())
            .map(|(a_val, b_val)| *a_val * *b_val)
            .reduce(|| T::zero(), |acc, val| acc + val);

        Ok(sum)
    }

    /// Compute L2 norm with SIMD optimization
    pub fn l2_norm<T: RealField + Copy + Send + Sync + ToPrimitive>(
        &self,
        input: &[T],
    ) -> Result<T> {
        let dot = self.dot_product(input, input)?;
        Ok(dot.sqrt())
    }

    /// Broadcasting addition: adds scalar to each element of vector
    pub fn broadcast_add<T: RealField + Copy + Send + Sync + ToPrimitive>(
        &self,
        input: &[T],
        scalar: T,
        result: &mut [T],
    ) -> Result<()> {
        if input.len() != result.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        result
            .par_iter_mut()
            .zip(input.par_iter())
            .for_each(|(r, val)| {
                *r = *val + scalar;
            });

        Ok(())
    }

    /// Broadcasting multiplication: multiplies each element by broadcasted vector (row-wise)
    pub fn broadcast_mul_vector<T: RealField + Copy + Send + Sync + ToPrimitive>(
        &self,
        matrix: &[T],
        matrix_cols: usize,
        vector: &[T],
        result: &mut [T],
    ) -> Result<()> {
        if !matrix.len().is_multiple_of(matrix_cols) {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }
        if vector.len() != matrix_cols {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }
        if result.len() != matrix.len() {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        result
            .par_chunks_mut(matrix_cols)
            .zip(matrix.par_chunks(matrix_cols))
            .for_each(|(result_row, matrix_row)| {
                // Interior row multiplication can use SIMD if we had a specific op for it,
                // for now we use zip/iter which LLVM can often vectorize.
                result_row
                    .iter_mut()
                    .zip(matrix_row.iter())
                    .zip(vector.iter())
                    .for_each(|((r, m), v)| {
                        *r = *m * *v;
                    });
            });

        Ok(())
    }

    /// Vectorized matrix-vector multiplication for CSR sparse matrices
    pub fn matvec_csr<T: RealField + Copy + Send + Sync + ToPrimitive>(
        &self,
        values: &[T],
        col_indices: &[usize],
        row_ptr: &[usize],
        x: &[T],
        y: &mut [T],
    ) -> Result<()> {
        if row_ptr.len() != y.len() + 1 {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        y.par_iter_mut().enumerate().for_each(|(i, y_i)| {
            let start = row_ptr[i];
            let end = row_ptr[i + 1];

            *y_i = (start..end)
                .map(|j| {
                    let col = col_indices[j];
                    values[j] * x[col]
                })
                .fold(T::zero(), |acc, val| acc + val);
        });

        Ok(())
    }

    /// Vectorized convolution for filtering operations
    pub fn convolution<T: RealField + Copy + Send + Sync + ToPrimitive>(
        &self,
        signal: &[T],
        kernel: &[T],
        result: &mut [T],
    ) -> Result<()> {
        let signal_len = signal.len();
        let kernel_len = kernel.len();
        let result_len = signal_len + kernel_len - 1;

        if result.len() != result_len {
            return Err(crate::error::MathError::DimensionMismatch.into());
        }

        result.par_iter_mut().enumerate().for_each(|(n, output)| {
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
            return Err(crate::error::MathError::InvalidInput(
                "Invalid dimensions for 3-point stencil".to_string(),
            )
            .into());
        }

        // Process interior points with SIMD where possible
        result.par_iter_mut().enumerate().for_each(|(i, r)| {
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
            return Err(crate::error::MathError::InvalidInput(
                "Invalid dimensions for 5-point stencil".to_string(),
            )
            .into());
        }

        result.par_iter_mut().enumerate().for_each(|(i, r)| {
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
