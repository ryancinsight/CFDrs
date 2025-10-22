//! Example backend abstraction pattern following Rust best practices.
//!
//! This module demonstrates zero-cost abstractions with backend selection,
//! following the persona configuration principles:
//! - Backend abstraction over wgpu-rs (prefer explicit backend trait)
//! - Zero-cost generics with trait bounds
//! - Feature-gated compilation for GPU support
//! - Iterator-based zero-copy operations
//!
//! # Example
//!
//! ```
//! use cfd_core::compute::backend_example::{select_backend, compute_squares};
//!
//! let backend = select_backend();
//! let storage = vec![1.0, 2.0, 3.0];
//! let result = compute_squares(&backend, &storage);
//! assert_eq!(result, vec![1.0, 4.0, 9.0]);
//! ```

use std::ops::Mul;

/// Storage trait for element-wise operations.
///
/// This trait enables zero-cost abstractions over different storage backends
/// while maintaining type safety and avoiding allocation in hot paths.
pub trait Storage<T> {
    /// Compute squares of all elements.
    ///
    /// # Examples
    ///
    /// ```
    /// use cfd_core::compute::backend_example::Storage;
    ///
    /// let data = vec![2.0, 3.0, 4.0];
    /// let result = data.compute_squares();
    /// assert_eq!(result, vec![4.0, 9.0, 16.0]);
    /// ```
    fn compute_squares(&self) -> Vec<T>
    where
        T: Copy + Mul<Output = T> + Default;
}

impl<T> Storage<T> for Vec<T> {
    fn compute_squares(&self) -> Vec<T>
    where
        T: Copy + Mul<Output = T> + Default,
    {
        // Iterator-based zero-copy operation (principle: prefer iterators/combinators)
        self.iter().map(|&x| x * x).collect()
    }
}

impl<T> Storage<T> for &[T] {
    fn compute_squares(&self) -> Vec<T>
    where
        T: Copy + Mul<Output = T> + Default,
    {
        // Slice implementation for borrowed data
        self.iter().map(|&x| x * x).collect()
    }
}

/// Compute backend trait for polymorphic execution.
///
/// Enables runtime selection between CPU and GPU execution paths
/// while maintaining compile-time optimization opportunities.
pub trait ComputeBackend {
    /// Compute squares using backend-specific optimization.
    ///
    /// # Type Parameters
    ///
    /// - `T`: Numeric type supporting multiplication (Copy + Mul)
    /// - `S`: Storage type implementing the Storage trait
    fn compute_squares<T, S>(&self, storage: &S) -> Vec<T>
    where
        S: Storage<T>,
        T: Copy + Mul<Output = T> + Default;
}

/// Backend implementation selector.
///
/// Uses feature flags to select between CPU and GPU implementations
/// at compile time, demonstrating zero-cost abstractions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    /// CPU execution path (scalar or SIMD)
    Cpu,
    /// GPU execution path (when gpu feature enabled)
    Gpu,
}

impl ComputeBackend for Backend {
    fn compute_squares<T, S>(&self, storage: &S) -> Vec<T>
    where
        S: Storage<T>,
        T: Copy + Mul<Output = T> + Default,
    {
        match self {
            Backend::Cpu => storage.compute_squares(),
            Backend::Gpu => {
                // GPU implementation would call wgpu kernels here
                // For now, delegate to storage (abstraction over wgpu-rs)
                storage.compute_squares()
            }
        }
    }
}

/// Select appropriate backend based on compilation features.
///
/// Uses `cfg!(feature = "gpu")` for zero-cost compile-time selection.
/// This pattern avoids runtime overhead while maintaining flexibility.
///
/// # Examples
///
/// ```
/// use cfd_core::compute::backend_example::select_backend;
///
/// let backend = select_backend();
/// // Returns Backend::Gpu if compiled with --features gpu
/// // Returns Backend::Cpu otherwise
/// ```
#[must_use]
pub fn select_backend() -> Backend {
    if cfg!(feature = "gpu") {
        Backend::Gpu
    } else {
        Backend::Cpu
    }
}

/// Element-wise mathematical operations for tensors.
///
/// Computes the square of each element using backend and storage abstractions.
///
/// # Examples
///
/// ```
/// use cfd_core::compute::backend_example::{select_backend, compute_squares};
///
/// let backend = select_backend();
/// let storage = vec![1.0, 2.0, 3.0];
/// let result = compute_squares(&backend, &storage);
/// assert_eq!(result, vec![1.0, 4.0, 9.0]);
/// ```
pub fn compute_squares<B, S, T>(backend: &B, storage: &S) -> Vec<T>
where
    B: ComputeBackend,
    S: Storage<T>,
    T: Copy + Mul<Output = T> + Default,
{
    backend.compute_squares(storage)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vec_compute_squares() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let result = data.compute_squares();
        assert_eq!(result, vec![1.0, 4.0, 9.0, 16.0]);
    }

    #[test]
    fn test_slice_compute_squares() {
        let data = [1.0, 2.0, 3.0];
        let result = data.as_slice().compute_squares();
        assert_eq!(result, vec![1.0, 4.0, 9.0]);
    }

    #[test]
    fn test_backend_cpu() {
        let backend = Backend::Cpu;
        let data = vec![2.0, 3.0, 4.0];
        let result = backend.compute_squares(&data);
        assert_eq!(result, vec![4.0, 9.0, 16.0]);
    }

    #[test]
    fn test_backend_gpu() {
        let backend = Backend::Gpu;
        let data = vec![2.0, 3.0, 4.0];
        let result = backend.compute_squares(&data);
        assert_eq!(result, vec![4.0, 9.0, 16.0]);
    }

    #[test]
    fn test_select_backend() {
        let backend = select_backend();
        // Should select GPU if feature enabled, CPU otherwise
        #[cfg(feature = "gpu")]
        assert_eq!(backend, Backend::Gpu);
        #[cfg(not(feature = "gpu"))]
        assert_eq!(backend, Backend::Cpu);
    }

    #[test]
    fn test_compute_squares_generic() {
        let backend = select_backend();
        let storage = vec![1.0, 2.0, 3.0];
        let result = compute_squares(&backend, &storage);
        assert_eq!(result, vec![1.0, 4.0, 9.0]);
    }

    #[test]
    fn test_compute_squares_integers() {
        let backend = Backend::Cpu;
        let storage = vec![1, 2, 3, 4];
        let result = compute_squares(&backend, &storage);
        assert_eq!(result, vec![1, 4, 9, 16]);
    }
}
