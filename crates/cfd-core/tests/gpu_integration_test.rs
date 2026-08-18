//! Integration tests for GPU compute functionality

#[cfg(feature = "gpu")]
mod gpu_tests {
    use cfd_core::compute::gpu::kernels::{
        advection::{AdvectionConfig, GpuAdvectionKernel},
        pressure::{GpuPressureKernel, PressureConfig},
        velocity::{GpuVelocityKernel, VelocityConfig},
    };
    use cfd_core::compute::gpu::{GpuBuffer, GpuContext};
    use cfd_core::compute::traits::{ComputeBuffer, DomainParams, KernelParams};
    use std::sync::Arc;

    #[test]
    fn test_gpu_context_creation() {
        // Skip test if no GPU available
        let context = if let Ok(ctx) = GpuContext::create() {
            Arc::new(ctx)
        } else {
            tracing::debug!("GPU unavailable for context integration test");
            return;
        };

        assert_eq!(context.backend_name(), "wgpu");

        assert!(context.max_work_group_size() > 0);
        assert!(context.max_buffer_size() > 0);
    }

    #[test]
    fn test_gpu_buffer_operations() {
        let context = match GpuContext::create() {
            Ok(ctx) => Arc::new(ctx),
            Err(error) => {
                tracing::debug!(error = %error, "GPU unavailable for buffer integration test");
                return;
            }
        };

        let size = 1000;
        let buffer = GpuBuffer::<f32>::new(context.clone(), size)
            .expect("invariant: GPU integration buffer allocates on a valid context");
        assert_eq!(buffer.size(), size);
    }

    #[test]
    fn test_advection_kernel() {
        let context = match GpuContext::create() {
            Ok(context) => Arc::new(context),
            Err(error) => {
                tracing::debug!(error = %error, "GPU unavailable for advection integration test");
                return;
            }
        };
        let kernel = GpuAdvectionKernel::new(context).expect("advection shader must compile");
        let config = AdvectionConfig::new([3, 3, 1], [1.0; 3], 0.5)
            .expect("invariant: valid advection integration grid is accepted");
        let scalar: Vec<f32> = (0..config.element_count())
            .map(|index| index as f32)
            .collect();
        let velocity = vec![0.0; config.element_count()];
        let mut output = vec![0.0; config.element_count()];

        kernel
            .execute(&scalar, &velocity, &velocity, config, &mut output)
            .expect("invariant: advection kernel executes on valid inputs");

        assert_eq!(output, scalar);
    }

    #[test]
    fn test_pressure_kernel() {
        let context = match GpuContext::create() {
            Ok(context) => Arc::new(context),
            Err(error) => {
                tracing::debug!(error = %error, "GPU unavailable for pressure integration test");
                return;
            }
        };
        let kernel = GpuPressureKernel::new(context).expect("pressure shaders must compile");
        let config = PressureConfig::new([3, 3, 3], [1.0; 3], 1.0)
            .expect("invariant: valid pressure integration grid is accepted");
        let pressure = vec![0.0; config.element_count()];
        let source = vec![0.0; config.element_count()];
        let mut residual = vec![1.0; config.element_count()];

        kernel
            .residual(&pressure, &source, config, &mut residual)
            .expect("invariant: pressure kernel executes on valid inputs");

        assert_eq!(residual, pressure);
    }

    #[test]
    fn test_velocity_kernel() {
        let context = match GpuContext::create() {
            Ok(context) => Arc::new(context),
            Err(error) => {
                tracing::debug!(error = %error, "GPU unavailable for velocity integration test");
                return;
            }
        };
        let kernel = GpuVelocityKernel::new(context).expect("velocity shaders must compile");
        let config = VelocityConfig::new([3, 3, 3], [1.0; 3], 0.5, 2.0)
            .expect("invariant: valid velocity integration grid is accepted");
        let velocity = vec![0.0; config.element_count()];
        let mut source = vec![1.0; config.element_count()];

        kernel
            .divergence_source(&velocity, &velocity, &velocity, config, &mut source)
            .expect("invariant: velocity kernel executes on valid inputs");

        assert_eq!(source, velocity);
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
