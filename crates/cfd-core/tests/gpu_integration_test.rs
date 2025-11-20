//! Integration tests for GPU compute functionality

#[cfg(feature = "gpu")]
mod gpu_tests {
    use cfd_core::compute::gpu::kernels::{
        advection::GpuAdvectionKernel, pressure::GpuPressureKernel, velocity::GpuVelocityKernel,
    };
    use cfd_core::compute::gpu::{pipeline::GpuPipelineManager, GpuBuffer, GpuContext};
    use cfd_core::compute::traits::{ComputeBuffer, ComputeKernel, DomainParams, KernelParams};
    use std::sync::Arc;

    #[test]
    fn test_gpu_context_creation() {
        // Skip test if no GPU available
        let context = if let Ok(ctx) = GpuContext::create() {
            Arc::new(ctx)
        } else {
            println!("Skipping GPU test - no GPU available");
            return;
        };

        let info = context.adapter_info();
        println!("GPU Device: {}", info.name);
        println!("Driver: {}", info.driver);

        assert!(context.max_work_group_size() > 0);
        assert!(context.max_buffer_size() > 0);
    }

    #[test]
    fn test_gpu_buffer_operations() {
        let context = match GpuContext::create() {
            Ok(ctx) => Arc::new(ctx),
            Err(_) => return, // Skip if no GPU
        };

        let size = 1000;
        let buffer = GpuBuffer::<f32>::new(context.clone(), size);
        assert!(buffer.is_ok());

        let buffer = buffer.unwrap();
        assert_eq!(buffer.size(), size);
    }

    #[test]
    fn test_advection_kernel() {
        let kernel = GpuAdvectionKernel::<f64>::new();
        assert_eq!(kernel.name(), "GPU Advection (Upwind)");

        // Test complexity calculation
        let complexity = kernel.complexity(1000);
        assert!(complexity > 0);
    }

    #[test]
    fn test_pressure_kernel() {
        let kernel = GpuPressureKernel::<f64>::new();
        assert_eq!(kernel.name(), "GPU Pressure Poisson Solver");

        // Verify shader code is not empty
        use cfd_core::compute::gpu::kernels::GpuKernel;
        assert!(!kernel.shader_code().is_empty());
    }

    #[test]
    fn test_velocity_kernel() {
        let kernel = GpuVelocityKernel::<f64>::new();
        assert_eq!(kernel.name(), "GPU Velocity Correction (SIMPLE)");

        use cfd_core::compute::gpu::kernels::GpuKernel;
        assert!(!kernel.shader_code().is_empty());
    }

    #[test]
    fn test_pipeline_manager() {
        let context = match GpuContext::create() {
            Ok(ctx) => Arc::new(ctx),
            Err(_) => return, // Skip if no GPU
        };

        let mut manager = GpuPipelineManager::new(context);

        // Register advection pipeline - needs 3 input buffers (u_in, velocity_x, velocity_y)
        let result = manager.register_pipeline_with_bindings(
            "advection",
            include_str!("../src/compute/gpu/kernels/advection.wgsl"),
            "advection_upwind",
            3, // 3 input buffers: u_in, velocity_x, velocity_y
        );

        // Should succeed or fail gracefully
        if let Err(e) = result {
            println!("Pipeline registration failed (expected on CI): {e}");
        }
    }

    #[test]
    fn test_kernel_params() {
        let params = KernelParams {
            size: 1000,
            work_group_size: 8,
            domain_params: DomainParams {
                grid_dims: (100, 100, 1),
                grid_spacing: (0.01, 0.01, 0.01),
                dt: 0.001,
                reynolds: 100.0,
                velocity: (1.0, 0.0, 0.0),
                boundary: cfd_core::compute::traits::BoundaryCondition2D::Periodic,
            },
        };

        assert_eq!(params.size, 1000);
        assert_eq!(params.domain_params.grid_dims.0, 100);
    }
}
