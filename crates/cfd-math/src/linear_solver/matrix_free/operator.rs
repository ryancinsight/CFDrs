//! Linear operator trait for matrix-free computations.
//!
//! This trait defines the interface for linear operators that can compute
//! matrix-vector products without explicitly storing the matrix. This enables
//! memory-efficient implementations of iterative solvers for large-scale problems.
//!
//! Supports both CPU and GPU execution with automatic dispatch based on
//! problem size and available hardware.

use cfd_core::error::{Error, Result};
use nalgebra::RealField;

/// Trait for linear operators that can compute matrix-vector products.
///
/// Implementors of this trait represent linear transformations `A: R^n -> R^n`
/// that can compute `y = A*x` efficiently without storing the matrix `A` explicitly.
///
/// This is fundamental for matrix-free iterative methods where the operator
/// is defined implicitly through the physics (e.g., discretized PDEs).
pub trait LinearOperator<T: RealField + Copy> {
    /// Apply the linear operator to a vector.
    ///
    /// Computes `y = A*x` where `A` is the implicit matrix represented by this operator.
    ///
    /// # Arguments
    ///
    /// * `x` - Input vector of size `self.size()`
    /// * `y` - Output vector of size `self.size()` (modified in place)
    ///
    /// # Errors
    ///
    /// Returns an error if the vector dimensions are incompatible or if
    /// the operator application fails numerically.
    fn apply(&self, x: &[T], y: &mut [T]) -> Result<()>;

    /// Return the dimension of the operator (n for A: R^n -> R^n).
    fn size(&self) -> usize;

    /// Optional: Return an estimate of the operator norm.
    ///
    /// Used for preconditioner scaling and convergence monitoring.
    /// Default implementation returns None (no estimate available).
    fn norm_estimate(&self) -> Option<T> {
        None
    }

    /// Optional: Apply the transpose operator.
    ///
    /// For non-symmetric operators, this enables transpose-free methods.
    /// Default implementation returns an error (not implemented).
    fn apply_transpose(&self, _x: &[T], _y: &mut [T]) -> Result<()> {
        Err(Error::InvalidConfiguration(
            "Transpose operator not implemented".to_string(),
        ))
    }

    /// Optional: Check if the operator is symmetric.
    ///
    /// Important for algorithm selection (e.g., CG requires symmetry).
    /// Default assumes non-symmetric.
    fn is_symmetric(&self) -> bool {
        false
    }

    /// Optional: Check if the operator is positive definite.
    ///
    /// Important for convergence guarantees in CG.
    /// Default assumes unknown.
    fn is_positive_definite(&self) -> Option<bool> {
        None
    }
}

/// Extended trait for GPU-accelerated linear operators.
///
/// This trait extends LinearOperator with GPU capabilities, allowing
/// operators to be executed on GPU hardware when available and beneficial.
pub trait GpuLinearOperator<T: RealField + Copy + bytemuck::Pod + bytemuck::Zeroable>:
    LinearOperator<T>
{
    /// Apply the operator using GPU acceleration.
    ///
    /// # Arguments
    ///
    /// * `gpu_context` - GPU context for execution
    /// * `input_buffer` - GPU buffer containing input vector
    /// * `output_buffer` - GPU buffer for output vector (modified in place)
    ///
    /// # Errors
    ///
    /// Returns an error if GPU execution fails or buffers are incompatible.
    fn apply_gpu(
        &self,
        gpu_context: &std::sync::Arc<cfd_core::compute::gpu::GpuContext>,
        input_buffer: &cfd_core::compute::gpu::GpuBuffer<T>,
        output_buffer: &mut cfd_core::compute::gpu::GpuBuffer<T>,
    ) -> Result<()>;

    /// Check if GPU acceleration is supported for this operator.
    ///
    /// Default implementation returns true if GPU context is available.
    fn supports_gpu(&self) -> bool {
        true
    }

    /// Estimate if GPU execution would be beneficial for given problem size.
    ///
    /// # Arguments
    ///
    /// * `problem_size` - Size of the linear system (n for n×n matrix)
    ///
    /// # Returns
    ///
    /// Returns true if GPU execution is recommended, false for CPU execution.
    fn should_use_gpu(&self, problem_size: usize) -> bool {
        // Default: use GPU for problems larger than 10,000 elements
        problem_size > 10_000
    }
}

/// Identity operator for testing and preconditioning.
pub struct IdentityOperator {
    size: usize,
}

impl IdentityOperator {
    /// Create a new identity operator of given size.
    pub fn new(size: usize) -> Self {
        Self { size }
    }
}

impl<T: RealField + Copy> LinearOperator<T> for IdentityOperator {
    fn apply(&self, x: &[T], y: &mut [T]) -> Result<()> {
        if x.len() != self.size || y.len() != self.size {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Vector dimensions don't match operator size".to_string(),
            ));
        }

        y.copy_from_slice(x);
        Ok(())
    }

    fn size(&self) -> usize {
        self.size
    }

    fn is_symmetric(&self) -> bool {
        true
    }

    fn is_positive_definite(&self) -> Option<bool> {
        Some(true)
    }
}

/// Scaled operator for preconditioning.
pub struct ScaledOperator<'a, T, Op> {
    operator: &'a Op,
    scale: T,
}

impl<'a, T, Op> ScaledOperator<'a, T, Op> {
    /// Create a new scaled operator.
    pub fn new(operator: &'a Op, scale: T) -> Self {
        Self { operator, scale }
    }
}

impl<T, Op> LinearOperator<T> for ScaledOperator<'_, T, Op>
where
    T: RealField + Copy,
    Op: LinearOperator<T>,
{
    fn apply(&self, x: &[T], y: &mut [T]) -> Result<()> {
        // First apply the base operator
        self.operator.apply(x, y)?;

        // Then scale the result
        for yi in y.iter_mut() {
            *yi *= self.scale;
        }

        Ok(())
    }

    fn size(&self) -> usize {
        self.operator.size()
    }

    fn norm_estimate(&self) -> Option<T> {
        self.operator
            .norm_estimate()
            .map(|norm| norm * self.scale.abs())
    }

    fn apply_transpose(&self, x: &[T], y: &mut [T]) -> Result<()> {
        // For scaled operators, transpose is scale * transpose(op)
        self.operator.apply_transpose(x, y)?;
        for yi in y.iter_mut() {
            *yi *= self.scale;
        }
        Ok(())
    }

    fn is_symmetric(&self) -> bool {
        self.operator.is_symmetric()
    }

    fn is_positive_definite(&self) -> Option<bool> {
        // Scaling by positive factor preserves positive definiteness
        if self.scale > T::zero() {
            self.operator.is_positive_definite()
        } else {
            Some(false)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_identity_operator() {
        let op: IdentityOperator = IdentityOperator::new(5);
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut y = vec![0.0; 5];

        <IdentityOperator as LinearOperator<f64>>::apply(&op, &x, &mut y).unwrap();

        for i in 0..5 {
            assert_relative_eq!(y[i], x[i], epsilon = 1e-10);
        }

        assert_eq!(<IdentityOperator as LinearOperator<f64>>::size(&op), 5);
        assert!(<IdentityOperator as LinearOperator<f64>>::is_symmetric(&op));
        assert_eq!(
            <IdentityOperator as LinearOperator<f64>>::is_positive_definite(&op),
            Some(true)
        );
    }

    #[test]
    fn test_scaled_operator() {
        let base_op = IdentityOperator::new(3);
        let op = ScaledOperator::new(&base_op, 2.0);
        let x = vec![1.0, 2.0, 3.0];
        let mut y = vec![0.0; 3];

        op.apply(&x, &mut y).unwrap();

        assert_relative_eq!(y[0], 2.0, epsilon = 1e-10);
        assert_relative_eq!(y[1], 4.0, epsilon = 1e-10);
        assert_relative_eq!(y[2], 6.0, epsilon = 1e-10);

        assert_eq!(op.size(), 3);
        assert!(op.is_symmetric());
        assert_eq!(op.is_positive_definite(), Some(true));
    }

    #[test]
    fn test_operator_dimension_mismatch() {
        let op = IdentityOperator::new(3);
        let x = vec![1.0, 2.0]; // Wrong size
        let mut y = vec![0.0; 3];

        assert!(op.apply(&x, &mut y).is_err());
    }
}
