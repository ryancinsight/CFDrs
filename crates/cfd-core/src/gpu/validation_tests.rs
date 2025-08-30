//! GPU kernel validation tests
//! Verifies GPU implementations match CPU results within numerical tolerance

#[cfg(test)]
mod tests {
    use super::super::field_ops::GpuFieldOps;
    use crate::compute::gpu::GpuContext;
    use approx::assert_relative_eq;
    use std::sync::Arc;

    const TOLERANCE: f32 = 1e-5;

    #[test]
    fn test_gpu_cpu_parity_add() {
        // Create GPU context
        let context = match GpuContext::create() {
            Ok(ctx) => Arc::new(ctx),
            Err(_) => {
                eprintln!("GPU not available, skipping test");
                return;
            }
        };

        let gpu_ops = GpuFieldOps::new(context);

        // Test data
        let a = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let b = vec![8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];

        // CPU reference
        let cpu_result: Vec<f32> = a.iter().zip(b.iter()).map(|(x, y)| x + y).collect();

        // GPU computation
        let mut gpu_result = vec![0.0; a.len()];
        gpu_ops.add_fields(&a, &b, &mut gpu_result);

        // Verify parity
        for (cpu, gpu) in cpu_result.iter().zip(gpu_result.iter()) {
            assert_relative_eq!(cpu, gpu, epsilon = TOLERANCE);
        }
    }

    #[test]
    fn test_gpu_laplacian_2d() {
        let context = match GpuContext::create() {
            Ok(ctx) => Arc::new(ctx),
            Err(_) => {
                eprintln!("GPU not available, skipping test");
                return;
            }
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
        gpu_ops.laplacian_2d(&field, nx, ny, dx, dy, &mut gpu_result);

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

    #[test]
    fn test_gpu_performance_characteristics() {
        // This test verifies GPU operations complete without error
        // Performance benchmarking would be done separately
        let context = match GpuContext::create() {
            Ok(ctx) => Arc::new(ctx),
            Err(_) => {
                eprintln!("GPU not available, skipping test");
                return;
            }
        };

        let gpu_ops = GpuFieldOps::new(context);

        // Large arrays to test GPU handling
        let size = 1024 * 1024; // 1M elements
        let a: Vec<f32> = (0..size).map(|i| i as f32).collect();
        let b: Vec<f32> = (0..size).map(|i| (size - i) as f32).collect();
        let mut result = vec![0.0; size];

        // Should complete without panic
        gpu_ops.add_fields(&a, &b, &mut result);

        // Spot check some values
        assert_relative_eq!(result[0], size as f32, epsilon = TOLERANCE);
        assert_relative_eq!(result[size - 1], size as f32, epsilon = TOLERANCE);
    }
}
