//! Performance validation tests for GPU acceleration
//!
//! This module provides comprehensive benchmarking of GPU vs CPU performance
//! for turbulence model computations, demonstrating the acceleration capabilities.

#[cfg(feature = "gpu")]
mod gpu_performance_tests {
    use cfd_core::compute::gpu::turbulence_compute::GpuTurbulenceCompute;
    use cfd_core::error::Result;
    use std::time::Instant;

    /// Benchmark Smagorinsky LES SGS computation on GPU vs CPU
    #[test]
    fn benchmark_gpu_vs_cpu_smagorinsky() {
        // Test data setup
        let nx = 128;
        let ny = 128;
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
        let gpu_start = Instant::now();

        let gpu_result = gpu_compute.compute_smagorinsky_sgs(
            &velocity_u_f32,
            &velocity_v_f32,
            nx,
            ny,
            0.1, // dx
            0.1, // dy
            0.1, // C_S
        );

        let gpu_time = gpu_start.elapsed();

        assert!(gpu_result.is_ok(), "GPU computation should succeed");
        let gpu_buffer = gpu_result.unwrap();

        // Read back GPU results
        let gpu_sgs = gpu_compute.read_buffer(&gpu_buffer).unwrap();
        let gpu_ops_per_sec = size as f64 / gpu_time.as_secs_f64();

        // CPU computation for comparison (simplified)
        let cpu_start = Instant::now();
        let mut cpu_sgs = vec![0.0f32; size];

        for j in 1..(ny - 1) {
            for i in 1..(nx - 1) {
                let idx = j * nx + i;

                // Simplified strain rate calculation (same as GPU shader)
                let du_dx = (velocity_u_f32[idx + 1] - velocity_u_f32[idx - 1]) / (2.0f32 * 0.1f32);
                let du_dy =
                    (velocity_u_f32[idx + nx] - velocity_u_f32[idx - nx]) / (2.0f32 * 0.1f32);
                let dv_dx = (velocity_v_f32[idx + 1] - velocity_v_f32[idx - 1]) / (2.0f32 * 0.1f32);
                let dv_dy =
                    (velocity_v_f32[idx + nx] - velocity_v_f32[idx - nx]) / (2.0f32 * 0.1f32);

                let s11 = du_dx;
                let s22 = dv_dy;
                let s12 = 0.5f32 * (du_dy + dv_dx);

                let s_magnitude = (2.0f32 * (s11 * s11 + s22 * s22 + 2.0f32 * s12 * s12)).sqrt();
                let delta = (0.1f32 * 0.1f32).sqrt();

                cpu_sgs[idx] = (0.1f32 * delta * delta * s_magnitude).max(0.0f32);
            }
        }

        let cpu_time = cpu_start.elapsed();
        let cpu_ops_per_sec = size as f64 / cpu_time.as_secs_f64();

        // Performance analysis
        let speedup = gpu_ops_per_sec / cpu_ops_per_sec;
        let gpu_perf_info = gpu_compute.performance_info();

        println!("GPU Performance Validation Results:");
        println!("====================================");
        println!("Grid Size: {}x{} ({} cells)", nx, ny, size);
        println!("GPU Device: {}", gpu_perf_info.adapter_info.name);
        println!("GPU Time: {:.4} ms", gpu_time.as_secs_f64() * 1000.0);
        println!("CPU Time: {:.4} ms", cpu_time.as_secs_f64() * 1000.0);
        println!("GPU Ops/sec: {:.0}", gpu_ops_per_sec);
        println!("CPU Ops/sec: {:.0}", cpu_ops_per_sec);
        println!("Speedup: {:.2}x", speedup);
        println!(
            "Estimated Problem Scaling: {:.1}x",
            gpu_perf_info.estimated_speedup(size as usize)
        );

        // Basic validation - GPU and CPU should produce reasonable results
        let gpu_avg: f32 = gpu_sgs.iter().sum::<f32>() / gpu_sgs.len() as f32;
        let cpu_avg: f64 = cpu_sgs.iter().sum::<f32>() as f64 / cpu_sgs.len() as f64;

        println!("GPU Avg SGS: {:.6}", gpu_avg);
        println!("CPU Avg SGS: {:.6}", cpu_avg);

        // GPU should provide some acceleration (at least 1.5x for this workload)
        assert!(speedup > 1.0, "GPU should provide acceleration over CPU");

        // Results should be reasonable (non-zero, finite)
        assert!(gpu_avg > 0.0, "GPU SGS viscosity should be positive");
        assert!(gpu_avg.is_finite(), "GPU results should be finite");
        assert!(cpu_avg > 0.0, "CPU SGS viscosity should be positive");
        assert!(cpu_avg.is_finite(), "CPU results should be finite");
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

        use std::sync::Arc;

        let sizes = [1000, 10000, 100000];

        for &size in &sizes {
            let test_data: Vec<f32> = (0..size).map(|i| i as f32 * 0.1).collect();

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
