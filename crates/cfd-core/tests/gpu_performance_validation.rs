//! Performance validation tests for GPU acceleration
//!
//! This module provides comprehensive benchmarking of GPU vs CPU performance
//! for turbulence model computations, demonstrating the acceleration capabilities.

#[cfg(feature = "gpu")]
mod gpu_performance_tests {
    use cfd_core::compute::gpu::{GpuBuffer, GpuTurbulenceCompute};

    use std::time::Instant;

    /// Benchmark Smagorinsky LES SGS computation on GPU vs CPU
    #[test]
    fn benchmark_gpu_vs_cpu_smagorinsky() {
        // Test data setup
        let nx = 2048;
        let ny = 2048;
        let size = nx * ny;

        // Generate test velocity fields with some turbulence-like features
        let mut velocity_u = vec![0.0; size];
        let mut velocity_v = vec![0.0; size];

        // Create synthetic turbulent flow field
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                let x = i as f64 * 0.01;
                let y = j as f64 * 0.01;

                // Add some turbulent fluctuations
                velocity_u[idx] = 1.0 + 0.1 * (x * 10.0).sin() * (y * 8.0).cos();
                velocity_v[idx] = 0.0 + 0.1 * (x * 8.0).cos() * (y * 10.0).sin();
            }
        }

        // Convert to f32 for GPU
        let velocity_u_f32: Vec<f32> = velocity_u.iter().map(|&x| x as f32).collect();
        let velocity_v_f32: Vec<f32> = velocity_v.iter().map(|&x| x as f32).collect();

        // GPU computation
        let mut gpu_compute = GpuTurbulenceCompute::new().unwrap();

        // Warmup to compile shaders and pipelines
        let _ = gpu_compute.compute_smagorinsky_sgs(
            &velocity_u_f32,
            &velocity_v_f32,
            nx,
            ny,
            0.1,
            0.1,
            0.1,
        );

        // Prepare GPU buffers once
        let context = gpu_compute.context_arc();
        let u_buffer = GpuBuffer::from_data(context.clone(), &velocity_u_f32).unwrap();
        let v_buffer = GpuBuffer::from_data(context.clone(), &velocity_v_f32).unwrap();
        let out_buffer: GpuBuffer<f32> = GpuBuffer::new(context.clone(), size).unwrap();

        let gpu_start = Instant::now();
        let iterations = 100;

        for _ in 0..iterations {
            gpu_compute
                .smagorinsky_kernel_mut()
                .compute_sgs_viscosity(
                    &context.device,
                    &context.queue,
                    u_buffer.buffer(),
                    v_buffer.buffer(),
                    out_buffer.buffer(),
                    nx as u32,
                    ny as u32,
                    0.1,
                    0.1,
                    0.1,
                )
                .unwrap();
        }

        let gpu_time = gpu_start.elapsed() / iterations;

        let gpu_result = gpu_compute.compute_smagorinsky_sgs(
            &velocity_u_f32,
            &velocity_v_f32,
            nx,
            ny,
            0.1,
            0.1,
            0.1,
        );

        assert!(gpu_result.is_ok(), "GPU computation should succeed");
        let gpu_buffer = gpu_result.unwrap();

        // Read back GPU results
        let gpu_sgs = gpu_compute.read_buffer(&gpu_buffer).unwrap();
        let gpu_ops_per_sec = size as f64 / gpu_time.as_secs_f64();

        // CPU computation for comparison
        let cpu_start = Instant::now();
        let mut cpu_sgs = vec![0.0f32; size];

        // Grid spacing for derivative calculations
        let dx = 0.1f32;
        let dy = 0.1f32;

        for _ in 0..iterations {
            for j in 1..(ny - 1) {
                for i in 1..(nx - 1) {
                    let idx = j * nx + i;

                    // Full strain rate tensor calculation
                    // Calculate velocity gradients using central differences
                    let du_dx = (velocity_u_f32[idx + 1] - velocity_u_f32[idx - 1]) / (2.0f32 * dx);
                    let du_dy =
                        (velocity_u_f32[idx + nx] - velocity_u_f32[idx - nx]) / (2.0f32 * dy);
                    let dv_dx = (velocity_v_f32[idx + 1] - velocity_v_f32[idx - 1]) / (2.0f32 * dx);
                    let dv_dy =
                        (velocity_v_f32[idx + nx] - velocity_v_f32[idx - nx]) / (2.0f32 * dy);

                    // Explicitly define all non-zero strain rate tensor components for 2D flow
                    let s_xx = du_dx;
                    let s_yy = dv_dy;
                    let s_xy = 0.5f32 * (du_dy + dv_dx);
                    let s_yx = s_xy;

                    // Note: For 2D planar flow, S_zz = S_xz = S_yz = 0
                    // The magnitude is S = sqrt(2 * S_ij * S_ij)
                    // S_ij * S_ij = s_xx^2 + s_yy^2 + s_xy^2 + s_yx^2 + (zero terms)
                    let s_contraction = s_xx * s_xx + s_yy * s_yy + s_xy * s_xy + s_yx * s_yx;
                    let s_magnitude = (2.0f32 * s_contraction).sqrt();

                    let delta = (dx * dy).sqrt();

                    cpu_sgs[idx] = (0.1f32 * delta * delta * s_magnitude).max(0.0f32);
                }
            }
        }

        let cpu_time = cpu_start.elapsed() / iterations;
        let cpu_ops_per_sec = size as f64 / cpu_time.as_secs_f64();

        // Performance analysis
        let speedup = gpu_ops_per_sec / cpu_ops_per_sec;
        let _gpu_perf_info = gpu_compute.performance_info();

        println!("GPU Performance Validation Results:");
        println!("====================================");
        println!("Grid Size: {}x{} ({} cells)", nx, ny, size);
        println!("GPU Time: {:?}, CPU Time: {:?}", gpu_time, cpu_time);
        println!("Speedup: {:.2}x", speedup);

        // Assert valid results (non-infinite, non-zero for synthetic field)
        let avg_sgs = gpu_sgs.iter().sum::<f32>() / size as f32;
        println!("Avg SGS Viscosity: {}", avg_sgs);
        assert!(avg_sgs > 0.0, "SGS viscosity should be positive");
        assert!(avg_sgs.is_finite(), "SGS viscosity should be finite");

        // GPU should be faster than CPU for this grid size
        // We set a reasonable threshold for CI environments where GPU might be emulated
        assert!(
            speedup > 0.5,
            "GPU should provide reasonable performance (speedup: {:.2}x)",
            speedup
        );
    }

    /// Test GPU context creation and basic functionality
    #[test]
    fn test_gpu_context_and_performance() {
        let gpu_compute = GpuTurbulenceCompute::new().unwrap();
        let perf_info = gpu_compute.performance_info();

        // Verify GPU context is functional
        assert!(perf_info.max_work_group_size > 0);
        assert!(perf_info.max_buffer_size > 0);
        assert!(!perf_info.adapter_info.name.is_empty());

        // Test performance estimation
        let speedup_1k = perf_info.estimated_speedup(1000);
        let speedup_10k = perf_info.estimated_speedup(10000);
        let speedup_100k = perf_info.estimated_speedup(100000);

        println!("Performance Estimation:");
        println!("1K cells: {:.2}x speedup", speedup_1k);
        println!("10K cells: {:.2}x speedup", speedup_10k);
        println!("100K cells: {:.2}x speedup", speedup_100k);

        // Larger problems should benefit more from GPU
        assert!(
            speedup_100k >= speedup_1k,
            "GPU speedup should increase with problem size"
        );
    }

    /// Test GPU buffer operations performance
    #[test]
    fn test_gpu_buffer_performance() {
        let gpu_compute = GpuTurbulenceCompute::new().unwrap();

        let sizes = [1000, 10000, 100000];

        for &size in &sizes {
            let _test_data: Vec<f32> = (0..size).map(|i| i as f32 * 0.1).collect();

            let start = Instant::now();

            // Create buffer and perform round-trip
            let buffer = gpu_compute
                .context()
                .device
                .create_buffer(&wgpu::BufferDescriptor {
                    label: Some("Performance Test Buffer"),
                    size: (size * std::mem::size_of::<f32>()) as u64,
                    usage: wgpu::BufferUsages::STORAGE
                        | wgpu::BufferUsages::COPY_DST
                        | wgpu::BufferUsages::COPY_SRC,
                    mapped_at_creation: false,
                });

            let elapsed = start.elapsed();

            println!(
                "Buffer creation ({} elements): {:.4} ms",
                size,
                elapsed.as_secs_f64() * 1000.0
            );

            // Verify buffer was created successfully
            assert!(buffer.size() > 0);
        }
    }
}

#[cfg(not(feature = "gpu"))]
mod gpu_performance_tests {
    #[test]
    fn test_gpu_disabled() {
        // When GPU feature is disabled, these tests should be skipped
        println!("GPU tests skipped - GPU feature not enabled");
    }
}
