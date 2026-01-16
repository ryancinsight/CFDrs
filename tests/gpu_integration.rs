//! GPU integration tests
//!
//! Tests GPU compute functionality with wgpu

#[cfg(feature = "gpu")]
mod gpu_tests {
    use cfd_core::compute::gpu::GpuContext;
    use cfd_suite::compute_unified::{Backend, UnifiedCompute};

    #[test]
    fn test_gpu_context_creation() {
        // Test GPU context creation (will use software renderer if no GPU)
        let result = GpuContext::create();

        // GPU context creation can fail in headless environments, which is OK
        if let Ok(context) = result {
            println!("GPU context created successfully");
            // Just check that we got a valid device
            assert!(context.limits.max_texture_dimension_2d > 0);
        } else {
            println!("GPU not available (expected in CI/headless environments)");
            // This is acceptable - GPU might not be available
        }
    }

    #[test]
    fn test_unified_compute_backend_selection() {
        let compute = UnifiedCompute::new();
        assert!(compute.is_ok());

        if let Ok(compute) = compute {
            let backend = compute.backend();
            match backend {
                Backend::Gpu => println!("Using GPU backend"),
                Backend::Simd => println!("Using SIMD backend"),
                Backend::Swar => println!("Using SWAR backend"),
            }
        }
    }

    #[test]
    fn test_simd_vector_operations() {
        // TODO: Replace unwrap-based error handling with proper Result types and error propagation
        // DEPENDENCIES: Add comprehensive error handling framework for GPU compute initialization
        // BLOCKED BY: Limited understanding of GPU compute failure modes and recovery strategies
        // PRIORITY: High - Essential for robust GPU testing and debugging
        let compute = UnifiedCompute::new().unwrap();

        let a = vec![1.0f32; 1000];
        let b = vec![2.0f32; 1000];
        let mut result = vec![0.0f32; 1000];

        // Test addition
        compute.vector_add_f32(&a, &b, &mut result).unwrap();
        assert_eq!(result[0], 3.0);
        assert_eq!(result[999], 3.0);

        // Test multiplication
        compute.vector_mul_f32(&a, &b, &mut result).unwrap();
        assert_eq!(result[0], 2.0);
        assert_eq!(result[999], 2.0);
    }

    #[test]
    fn test_matrix_vector_multiplication() {
        // TODO: Replace unwrap-based error handling with proper Result types and error propagation
        // DEPENDENCIES: Add comprehensive error handling framework for GPU compute initialization
        // BLOCKED BY: Limited understanding of GPU compute failure modes and recovery strategies
        // PRIORITY: High - Essential for robust GPU testing and debugging
        let compute = UnifiedCompute::new().unwrap();

        // 3x3 matrix (row-major)
        let matrix = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let vector = vec![1.0, 1.0, 1.0];
        let mut result = vec![0.0f32; 3];

        compute
            .matvec_f32(&matrix, &vector, &mut result, 3, 3)
            // TODO: Replace unwrap-based error handling with proper Result types and error propagation
            // DEPENDENCIES: Add comprehensive error handling framework for GPU compute operations
            // BLOCKED BY: Limited understanding of GPU compute operation failure modes and recovery strategies
            // PRIORITY: High - Essential for robust GPU testing and debugging
            .unwrap();

        // Expected: [6.0, 15.0, 24.0]
        assert_eq!(result[0], 6.0);
        assert_eq!(result[1], 15.0);
        assert_eq!(result[2], 24.0);
    }
}

#[cfg(not(feature = "gpu"))]
mod gpu_tests {
    #[test]
    fn test_gpu_feature_disabled() {
        println!("GPU feature is disabled, skipping GPU tests");
        assert!(true);
    }
}
