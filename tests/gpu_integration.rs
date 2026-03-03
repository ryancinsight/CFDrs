//! GPU integration tests
//!
//! Tests GPU compute functionality with wgpu.
//! Detailed GPU kernel tests live in `crates/cfd-core/tests/gpu_integration_test.rs`.

#[cfg(feature = "gpu")]
mod gpu_tests {
    use cfd_core::compute::backend::ComputeCapability;
    use cfd_core::compute::gpu::GpuContext;
    use cfd_core::compute::traits::ComputeBackend;

    #[test]
    fn test_gpu_context_creation() {
        // GPU context creation can fail in headless environments, which is OK
        if let Ok(context) = GpuContext::create() {
            println!("GPU context created successfully");
            assert!(context.limits.max_texture_dimension_2d > 0);
        } else {
            println!("GPU not available (expected in CI/headless environments)");
        }
    }

    #[test]
    fn test_compute_capability_detection() {
        let caps = ComputeCapability::detect();
        assert!(
            caps.backends.contains(&ComputeBackend::Cpu),
            "CPU backend must always be available"
        );
        assert!(caps.compute_units > 0);
        assert!(caps.available_memory > 0);
    }

    #[test]
    fn test_backend_selection_by_problem_size() {
        let caps = ComputeCapability::detect();

        // Tiny problem → always CPU
        let small = caps.select_backend(100);
        assert_eq!(small, ComputeBackend::Cpu);

        // Large problem → GPU if available, else SIMD or CPU
        let large = caps.select_backend(200_000);
        if caps.backends.contains(&ComputeBackend::Gpu) {
            assert_eq!(large, ComputeBackend::Gpu);
        } else if caps.backends.contains(&ComputeBackend::Simd) {
            assert_eq!(large, ComputeBackend::Simd);
        } else {
            assert_eq!(large, ComputeBackend::Cpu);
        }
    }
}

#[cfg(not(feature = "gpu"))]
mod gpu_tests {
    #[test]
    fn test_gpu_feature_disabled() {
        // GPU feature is off — kernel-level tests live in cfd-core
        assert!(true);
    }
}
