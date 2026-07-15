//! SIMD-optimized vector operations for performance-critical paths
//!
//! This module provides vectorized operations using platform-specific SIMD instructions
//! when available, with automatic fallback to scalar operations.

use eunomia::NumericElement;
use leto::Array1;
use leto_ops::{spmv, CsrMatrix, Scalar};
use moirai::{map_collect_index_with, reduce_index_with, Adaptive};

/// Trait for SIMD-optimized vector operations
pub trait SimdVectorOps<T: NumericElement + Send + Sync> {
    /// Compute element-wise product with SIMD optimization
    fn simd_mul(&self, other: &Self) -> Self;

    /// Compute dot product with SIMD optimization
    fn simd_dot(&self, other: &Self) -> T;

    /// Compute L2 norm with SIMD optimization
    fn simd_norm(&self) -> T;

    /// Apply element-wise function with parallel execution
    fn par_map<F>(&self, f: F) -> Self
    where
        F: Fn(T) -> T + Send + Sync;
}

impl<T: NumericElement + Send + Sync> SimdVectorOps<T> for Array1<T> {
    #[inline]
    fn simd_mul(&self, other: &Self) -> Self {
        assert_eq!(self.shape(), other.shape(), "Vector dimensions must match");

        // Parallel element-wise product; Adaptive routes small inputs to the
        // sequential path automatically.
        let len = self.shape()[0];
        let data = map_collect_index_with::<Adaptive, _, _>(len, |i| self[i] * other[i]);
        Array1::from_shape_vec([data.len()], data)
            .expect("invariant: SIMD product preserves vector length")
    }

    #[inline]
    fn simd_dot(&self, other: &Self) -> T {
        assert_eq!(self.shape(), other.shape(), "Vector dimensions must match");

        // Parallel dot product (index-aligned reduction over both slices).
        let len = self.shape()[0];
        reduce_index_with::<Adaptive, _, _, _>(
            len,
            T::ZERO,
            |i| self[i] * other[i],
            |acc, x| acc + x,
        )
    }

    #[inline]
    fn simd_norm(&self) -> T {
        // Parallel sum of squares, then sqrt.
        let len = self.shape()[0];
        let sum_sq = reduce_index_with::<Adaptive, _, _, _>(
            len,
            T::ZERO,
            |i| self[i] * self[i],
            |acc, x| acc + x,
        );
        NumericElement::sqrt(sum_sq)
    }

    #[inline]
    fn par_map<F>(&self, f: F) -> Self
    where
        F: Fn(T) -> T + Send + Sync,
    {
        let len = self.shape()[0];
        let data = map_collect_index_with::<Adaptive, _, _>(len, |i| f(self[i]));
        Array1::from_shape_vec([data.len()], data)
            .expect("invariant: SIMD map preserves vector length")
    }
}

/// Sparse matrix-vector multiplication
pub fn sparse_matvec<T: Scalar + Send + Sync>(
    matrix: &CsrMatrix<T>,
    vector: &Array1<T>,
) -> crate::error::Result<Array1<T>> {
    spmv(matrix, &vector.view()).map_err(|err| {
        cfd_core::error::Error::InvalidInput(format!("Leto sparse matvec error: {err}"))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_simd_operations() {
        let n = 2000;
        let v1 = Array1::from_elem([n], 2.0);
        let v2 = Array1::from_elem([n], 3.0);

        // Test SIMD multiplication
        let result = v1.simd_mul(&v2);
        assert_eq!(result.shape(), [n]);
        assert_relative_eq!(result[[0]], 6.0);

        // Test SIMD dot product
        let dot = v1.simd_dot(&v2);
        assert_relative_eq!(dot, 6.0 * n as f64);

        // Test SIMD norm
        let norm = v1.simd_norm();
        assert_relative_eq!(norm, (4.0 * n as f64).sqrt());
    }
}
