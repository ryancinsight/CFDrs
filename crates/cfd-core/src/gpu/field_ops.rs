//! Refactored GPU field operations with proper separation of concerns
//!
//! This module delegates to specialized kernels for individual operations
//! rather than maintaining a monolithic implementation.

use crate::compute::gpu::GpuContext;
use crate::gpu::kernels::{FieldAddKernel, FieldMulKernel, Laplacian2DKernel};
use std::sync::Arc;

/// Facade for GPU-accelerated field operations
pub struct GpuFieldOps {
    add_kernel: FieldAddKernel,
    mul_kernel: FieldMulKernel,
    laplacian_kernel: Laplacian2DKernel,
}

impl GpuFieldOps {
    /// Create new GPU field operations handler
    #[must_use]
    pub fn new(context: Arc<GpuContext>) -> Self {
        Self {
            add_kernel: FieldAddKernel::new(context.clone()),
            mul_kernel: FieldMulKernel::new(context.clone()),
            laplacian_kernel: Laplacian2DKernel::new(context),
        }
    }

    /// Add two fields element-wise: result = a + b
    pub fn add_fields(&self, a: &[f32], b: &[f32], result: &mut [f32]) {
        self.add_kernel.execute(a, b, result);
    }

    /// Multiply field by scalar: result = field * scalar
    pub fn multiply_field(&self, field: &[f32], scalar: f32, result: &mut [f32]) {
        self.mul_kernel.execute(field, scalar, result);
    }

    /// Compute 2D Laplacian
    pub fn laplacian_2d(
        &self,
        field: &[f32],
        nx: usize,
        ny: usize,
        dx: f32,
        dy: f32,
        result: &mut [f32],
    ) {
        self.laplacian_kernel.execute(field, nx, ny, dx, dy, result);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_field_operations() {
        // Test only runs if GPU is available
        if let Ok(context) = GpuContext::create() {
            let ops = GpuFieldOps::new(Arc::new(context));

            // Test field addition
            let a = vec![1.0_f32; 1024];
            let b = vec![2.0_f32; 1024];
            let mut result = vec![0.0_f32; 1024];

            ops.add_fields(&a, &b, &mut result);

            for val in &result {
                assert!((val - 3.0).abs() < 1e-6, "Field addition failed");
            }

            // Test scalar multiplication
            let field = vec![2.0_f32; 1024];
            let mut result = vec![0.0_f32; 1024];

            ops.multiply_field(&field, 3.0, &mut result);

            for val in &result {
                assert!((val - 6.0).abs() < 1e-6, "Scalar multiplication failed");
            }
        }
    }
}
