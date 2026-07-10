//! GPU kernel validation tests
//! Verifies GPU implementations match CPU results within numerical tolerance

#[cfg(test)]
mod tests {
    use super::super::field_ops::GpuFieldOps;
    use crate::compute::gpu::GpuContext;
    use approx::assert_relative_eq;
    use std::sync::Arc;

    #[test]
    fn test_gpu_laplacian_2d() {
        let context = if let Ok(ctx) = GpuContext::create() {
            Arc::new(ctx)
        } else {
            eprintln!("GPU not available, skipping test");
            return;
        };

        let gpu_ops = GpuFieldOps::new(context);

        // 4x4 test field
        let nx = 4u32;
        let ny = 4u32;
        let dx = 1.0f32;
        let dy = 1.0f32;

        // Simple test field: f(x,y) = x^2 + y^2
        // Laplacian should be 4 everywhere (d²f/dx² + d²f/dy² = 2 + 2 = 4)
        let field: Vec<f32> = (0..16)
            .map(|i| {
                let x = (i % 4) as f32;
                let y = (i / 4) as f32;
                x * x + y * y
            })
            .collect();

        let mut gpu_result = vec![0.0; field.len()];
        gpu_ops.laplacian_2d(&field, nx as usize, ny as usize, dx, dy, &mut gpu_result);

        // Interior points should have Laplacian ≈ 4
        // Boundary handling may differ, so only check interior
        for j in 1..3 {
            for i in 1..3 {
                let idx = (j * nx + i) as usize;
                // Finite difference approximation
                let laplacian_x = (field[idx - 1] - 2.0 * field[idx] + field[idx + 1]) / (dx * dx);
                let laplacian_y = (field[idx - nx as usize] - 2.0 * field[idx]
                    + field[idx + nx as usize])
                    / (dy * dy);
                let expected = laplacian_x + laplacian_y;

                assert_relative_eq!(gpu_result[idx], expected, epsilon = 0.1);
            }
        }
    }
}
