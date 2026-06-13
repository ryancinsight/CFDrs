//! SIMD-optimized vector operations for performance-critical paths
//!
//! This module provides vectorized operations using platform-specific SIMD instructions
//! when available, with automatic fallback to scalar operations.

use moirai::prelude::ParallelSlice;
use moirai::{map_collect_index_with, reduce_index_with, Adaptive};
use nalgebra::{DVector, RealField};

/// Trait for SIMD-optimized vector operations
pub trait SimdVectorOps<T: RealField + Copy + Send + Sync> {
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

impl<T: RealField + Copy + Send + Sync> SimdVectorOps<T> for DVector<T> {
    #[inline]
    fn simd_mul(&self, other: &Self) -> Self {
        assert_eq!(self.len(), other.len(), "Vector dimensions must match");

        // Parallel element-wise product; Adaptive routes small inputs to the
        // sequential path automatically.
        let a = self.data.as_slice();
        let b = other.data.as_slice();
        let data = map_collect_index_with::<Adaptive, _, _>(a.len(), |i| a[i] * b[i]);
        DVector::from_vec(data)
    }

    #[inline]
    fn simd_dot(&self, other: &Self) -> T {
        assert_eq!(self.len(), other.len(), "Vector dimensions must match");

        // Parallel dot product (index-aligned reduction over both slices).
        let a = self.data.as_slice();
        let b = other.data.as_slice();
        reduce_index_with::<Adaptive, _, _, _>(a.len(), T::zero(), |i| a[i] * b[i], |acc, x| acc + x)
    }

    #[inline]
    fn simd_norm(&self) -> T {
        // Parallel sum of squares, then sqrt.
        let sum_sq = self
            .data
            .as_slice()
            .par()
            .map_reduce(T::zero(), |&x| x * x, |acc, x| acc + x);
        sum_sq.sqrt()
    }

    #[inline]
    fn par_map<F>(&self, f: F) -> Self
    where
        F: Fn(T) -> T + Send + Sync,
    {
        let data: Vec<T> = self.data.as_slice().par().map_collect(|&x| f(x));
        DVector::from_vec(data)
    }
}

/// Sparse matrix-vector multiplication
pub fn sparse_matvec<T: RealField + Copy + Send + Sync>(
    matrix: &nalgebra_sparse::CsrMatrix<T>,
    vector: &DVector<T>,
) -> DVector<T> {
    // Use nalgebra's sparse matrix-vector multiplication
    matrix * vector
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_simd_operations() {
        let n = 2000;
        let v1 = DVector::from_element(n, 2.0);
        let v2 = DVector::from_element(n, 3.0);

        // Test SIMD multiplication
        let result = v1.simd_mul(&v2);
        assert_eq!(result.len(), n);
        assert_relative_eq!(result[0], 6.0);

        // Test SIMD dot product
        let dot = v1.simd_dot(&v2);
        assert_relative_eq!(dot, 6.0 * n as f64);

        // Test SIMD norm
        let norm = v1.simd_norm();
        assert_relative_eq!(norm, (4.0 * n as f64).sqrt());
    }
}
