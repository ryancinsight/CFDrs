//! Tests for compute backend

use super::*;
use crate::compute::cpu::CpuAdvectionKernel;
use crate::compute::traits::{BoundaryCondition2D, ComputeBackend, DomainParams, KernelParams};

#[test]
fn test_backend_detection() {
        let capability = ComputeCapability::detect();

        // CPU should always be available
        assert!(capability.backends.contains(&ComputeBackend::Cpu));

        // Check if detected backends are actually available
        for backend in &capability.backends {
            assert!(backend.is_available());
        }

        println!("Detected backends: {:?}", capability.backends);
        println!("Preferred backend: {:?}", capability.preferred_backend);
        println!("Compute units: {}", capability.compute_units);
    }

    #[test]
    fn test_backend_selection() {
        let capability = ComputeCapability::detect();

        // Small problem should use CPU
        let small_backend = capability.select_backend(100);
        assert_eq!(small_backend, ComputeBackend::Cpu);

        // Medium problem should use SIMD if available
        let medium_backend = capability.select_backend(10_000);
        if capability.backends.contains(&ComputeBackend::Simd) {
            assert_eq!(medium_backend, ComputeBackend::Simd);
        } else {
            assert_eq!(medium_backend, ComputeBackend::Cpu);
        }

        // Large problem should use GPU if available
        let large_backend = capability.select_backend(1_000_000);
        if capability.backends.contains(&ComputeBackend::Gpu) {
            assert_eq!(large_backend, ComputeBackend::Gpu);
        } else if capability.backends.contains(&ComputeBackend::Simd) {
            assert_eq!(large_backend, ComputeBackend::Simd);
        } else {
            assert_eq!(large_backend, ComputeBackend::Cpu);
        }
    }

    #[test]
    fn test_cpu_advection_kernel_linear_exactness() {
        let kernel = CpuAdvectionKernel::<f64>::new();

        // Grid and parameters
        let nx = 16;
        let ny = 12;
        let size = nx * ny;
        let dx = 0.1;
        let dy = 0.2;
        let dt = 0.01;
        let u = 1.5; // constant velocity x
        let v = -0.7; // constant velocity y
        let a = 2.0; // linear coefficient in x
        let b = -1.0; // linear coefficient in y

        // Construct linear field: phi(i,j) = a * (i*dx) + b * (j*dy)
        let mut input = vec![0.0f64; size];
        for j in 0..ny {
            for i in 0..nx {
                input[j * nx + i] = a * (i as f64 * dx) + b * (j as f64 * dy);
            }
        }
        let mut output = vec![0.0f64; size];

        let params = KernelParams {
            size,
            work_group_size: 8,
            domain_params: DomainParams {
                grid_dims: (nx, ny, 1),
                grid_spacing: (dx, dy, 1.0),
                dt,
                reynolds: 100.0,
                velocity: (u, v, 0.0),
                boundary: BoundaryCondition2D::DirichletZero,
            },
        };

        // Execute one step
        kernel.execute(&input, &mut output, params).unwrap();

        // Analytic update for linear field under constant velocity
        // phi_new = phi_old - dt * (u*a + v*b)
        let shift = dt * (u * a + v * b);

        // Validate exactness for interior points (away from boundaries)
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;
                let expected = input[idx] - shift;
                assert!((output[idx] - expected).abs() < 1e-12);
            }
        }

        // Check complexity calculation
        let complexity = kernel.complexity(size);
        assert_eq!(complexity, size * 10);
    }

    #[test]
    fn test_dispatcher_creation() {
        // Test automatic backend selection
        let dispatcher = ComputeDispatcher::new().unwrap();
        println!(
            "Dispatcher using backend: {:?}",
            dispatcher.current_backend()
        );

        // Test specific backend selection
        let cpu_dispatcher = ComputeDispatcher::with_backend(ComputeBackend::Cpu).unwrap();
        assert_eq!(cpu_dispatcher.current_backend(), ComputeBackend::Cpu);
    }

    #[cfg(feature = "gpu")]
    #[test]
    fn test_gpu_context_creation() {
        use crate::compute::gpu::GpuContext;

        // Try to create GPU context
        match GpuContext::create() {
            Ok(context) => {
                let info = context.adapter_info();
                println!("GPU Adapter: {:?}", info.name);
                println!("GPU Backend: {:?}", info.backend);
                println!("Max work group size: {}", context.max_work_group_size());
                println!(
                    "Max buffer size: {} MB",
                    context.max_buffer_size() / (1024 * 1024)
                );
            }
            Err(e) => {
                println!("GPU not available: {e}");
            }
        }
    }

    #[cfg(feature = "gpu")]
    #[test]
    fn test_gpu_buffer() {
        use crate::compute::gpu::{GpuBuffer, GpuContext};
        use crate::compute::traits::ComputeBuffer;
        use std::sync::Arc;

        // Skip if GPU not available
        let context = match GpuContext::create() {
            Ok(ctx) => Arc::new(ctx),
            Err(_) => return,
        };

        // Test buffer creation
        let size = 1000;
        let buffer = GpuBuffer::<f32>::new(context.clone(), size).unwrap();
        assert_eq!(buffer.size(), size);

        // Test buffer with data
        let data: Vec<f32> = (0..100).map(|i| i as f32).collect();
        let mut buffer = GpuBuffer::from_data(context, &data).unwrap();
        assert_eq!(buffer.size(), data.len());

        // Test read
        let read_data = buffer.read().unwrap();
        assert_eq!(read_data, data);

        // Test write
        let updated_data: Vec<f32> = (0..100).map(|i| (i * 2) as f32).collect();
        buffer.write(&updated_data).unwrap();
        let read_data = buffer.read().unwrap();
        assert_eq!(read_data, updated_data);
    }
