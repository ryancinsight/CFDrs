//! SIMD-optimized vector operations for performance-critical paths
//!
//! This module provides vectorized operations using platform-specific SIMD instructions
//! when available, with automatic fallback to scalar operations.

use nalgebra::{DVector, RealField};
use rayon::prelude::*;
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
    }

    fn simd_mul(&self, other: &Self) -> Self {
        assert_eq!(self.len(), other.len(), "Vector dimensions must match");
        // For large vectors, use parallel execution
        if self.len() > 1000 {
            let data: Vec<T> = self
                .data
                .as_slice()
                .par_iter()
                .zip(other.data.as_slice().par_iter())
                .map(|(&a, &b)| a * b)
                .collect();
            DVector::from_vec(data)
        } else {
            // For small vectors, use sequential with potential auto-vectorization
            self.component_mul(other)
        }
    }
    fn simd_dot(&self, other: &Self) -> T {
        // For large vectors, use parallel reduction
            self.data
                .reduce(|| T::zero(), |acc, x| acc + x)
            self.dot(other)
    }

    fn simd_norm(&self) -> T {
            let sum_sq = self
                .map(|&x| x * x)
                .reduce(|| T::zero(), |acc, x| acc + x);
            sum_sq.sqrt()
            self.norm()
        F: Fn(T) -> T + Send + Sync,
    {
            let data: Vec<T> = self.data.as_slice().par_iter().map(|&x| f(x)).collect();
            self.map(f)
/// Sparse matrix-vector multiplication
    }

pub fn sparse_matvec<T: RealField + Copy + Send + Sync>(
    matrix: &nalgebra_sparse::CsrMatrix<T>,
    vector: &DVector<T>,
) -> DVector<T> {
    // Use nalgebra's sparse matrix-vector multiplication
    matrix * vector
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    #[test]
    }

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
}
}
}
