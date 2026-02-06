//! GPU turbulence acceleration benchmarks
//!
//! Benchmarks comparing GPU-accelerated turbulence models against CPU implementations
//! to validate performance improvements and accuracy.

use cfd_2d::physics::turbulence::{
    des::{DESConfig, DetachedEddySimulation},
    les_smagorinsky::{SmagorinskyConfig, SmagorinskyLES},
    LESTurbulenceModel,
};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nalgebra::DMatrix;

/// Create test velocity fields for benchmarking
fn create_benchmark_velocity_fields(nx: usize, ny: usize) -> (DMatrix<f64>, DMatrix<f64>) {
    let mut velocity_u = DMatrix::zeros(nx, ny);
    let mut velocity_v = DMatrix::zeros(nx, ny);

    // Create realistic turbulent flow field
    for i in 0..nx {
        for j in 0..ny {
            let x = i as f64 / nx as f64;
            let y = j as f64 / ny as f64;

            // Base flow + turbulence
            velocity_u[(i, j)] = 1.0
                + 0.1 * (2.0 * std::f64::consts::PI * x).sin()
                + 0.05 * (4.0 * std::f64::consts::PI * y).cos();
            velocity_v[(i, j)] = 0.1 * (2.0 * std::f64::consts::PI * y).sin()
                + 0.05 * (4.0 * std::f64::consts::PI * x).cos();
        }
    }

    (velocity_u, velocity_v)
}

/// Benchmark Smagorinsky LES CPU vs GPU
fn bench_smagorinsky_les(c: &mut Criterion) {
    let mut group = c.benchmark_group("Smagorinsky LES");

    for &size in [32, 64, 128].iter() {
        let (velocity_u, velocity_v) = create_benchmark_velocity_fields(size, size);
        let pressure = DMatrix::zeros(size, size);

        // CPU benchmark
        let config_cpu = SmagorinskyConfig {
            use_gpu: false,
            ..Default::default()
        };
        let mut les_cpu = SmagorinskyLES::new(size, size, 0.01, 0.01, config_cpu);

        group.bench_function(format!("CPU {}x{}", size, size), |b| {
            b.iter(|| {
                black_box(les_cpu.update(
                    &velocity_u,
                    &velocity_v,
                    &pressure,
                    1.0,
                    0.01,
                    0.001,
                    0.01,
                    0.01,
                ))
                .unwrap();
            });
        });

        // GPU benchmark (if available)
        #[cfg(feature = "gpu")]
        {
            let config_gpu = SmagorinskyConfig {
                use_gpu: true,
                ..Default::default()
            };
            let mut les_gpu = SmagorinskyLES::new(size, size, 0.01, 0.01, config_gpu);

            group.bench_function(format!("GPU {}x{}", size, size), |b| {
                b.iter(|| {
                    black_box(les_gpu.update(
                        &velocity_u,
                        &velocity_v,
                        &pressure,
                        1.0,
                        0.01,
                        0.001,
                        0.01,
                        0.01,
                    ))
                    .unwrap();
                });
            });
        }
    }

    group.finish();
}

/// Benchmark DES CPU vs GPU
fn bench_des(c: &mut Criterion) {
    let mut group = c.benchmark_group("DES");

    for &size in [32, 64, 128].iter() {
        let (velocity_u, velocity_v) = create_benchmark_velocity_fields(size, size);
        let pressure = DMatrix::zeros(size, size);

        // CPU benchmark
        let config_cpu = DESConfig {
            use_gpu: false,
            ..Default::default()
        };
        let mut des_cpu = DetachedEddySimulation::new(size, size, 0.01, 0.01, config_cpu, &[]);

        group.bench_function(format!("CPU {}x{}", size, size), |b| {
            b.iter(|| {
                black_box(des_cpu.update(
                    &velocity_u,
                    &velocity_v,
                    &pressure,
                    1.0,
                    0.01,
                    0.001,
                    0.01,
                    0.01,
                ))
                .unwrap();
            });
        });

        // GPU benchmark (if available)
        #[cfg(feature = "gpu")]
        {
            let config_gpu = DESConfig {
                use_gpu: true,
                ..Default::default()
            };
            let mut des_gpu = DetachedEddySimulation::new(size, size, 0.01, 0.01, config_gpu, &[]);

            group.bench_function(format!("GPU {}x{}", size, size), |b| {
                b.iter(|| {
                    black_box(des_gpu.update(
                        &velocity_u,
                        &velocity_v,
                        &pressure,
                        1.0,
                        0.01,
                        0.001,
                        0.01,
                        0.01,
                    ))
                    .unwrap();
                });
            });
        }
    }

    group.finish();
}

/// Benchmark strain rate tensor computation
fn bench_strain_rate_computation(c: &mut Criterion) {
    let mut group = c.benchmark_group("Strain Rate Computation");

    for &size in [64, 128, 256].iter() {
        let (velocity_u, velocity_v) = create_benchmark_velocity_fields(size, size);

        group.bench_function(format!("CPU {}x{}", size, size), |b| {
            b.iter(|| {
                let mut strain = DMatrix::zeros(size, size);
                for i in 1..size - 1 {
                    for j in 1..size - 1 {
                        let du_dx = (velocity_u[(i + 1, j)] - velocity_u[(i - 1, j)]) / 0.02;
                        let du_dy = (velocity_u[(i, j + 1)] - velocity_u[(i, j - 1)]) / 0.02;
                        let dv_dx = (velocity_v[(i + 1, j)] - velocity_v[(i - 1, j)]) / 0.02;
                        let dv_dy = (velocity_v[(i, j + 1)] - velocity_v[(i, j - 1)]) / 0.02;

                        let s11 = du_dx;
                        let s22 = dv_dy;
                        let s12 = 0.5 * (du_dy + dv_dx);

                        strain[(i, j)] =
                            black_box((2.0 * s11 * s11 + 2.0 * s22 * s22 + 4.0 * s12 * s12).sqrt());
                    }
                }
                black_box(strain);
            });
        });
    }

    group.finish();
}

/// Accuracy validation benchmark
fn bench_accuracy_validation(c: &mut Criterion) {
    let mut group = c.benchmark_group("Accuracy Validation");

    let size = 64;
    let (velocity_u, velocity_v) = create_benchmark_velocity_fields(size, size);
    let pressure = DMatrix::zeros(size, size);

    // Compare CPU and GPU results for accuracy
    group.bench_function("Smagorinsky Accuracy Check", |b| {
        b.iter(|| {
            // CPU computation
            let config_cpu = SmagorinskyConfig {
                use_gpu: false,
                ..Default::default()
            };
            let mut les_cpu = SmagorinskyLES::new(size, size, 0.01, 0.01, config_cpu);
            les_cpu
                .update(
                    &velocity_u,
                    &velocity_v,
                    &pressure,
                    1.0,
                    0.01,
                    0.001,
                    0.01,
                    0.01,
                )
                .unwrap();
            let cpu_viscosity = les_cpu.get_turbulent_viscosity_field().clone();

            #[cfg(feature = "gpu")]
            {
                // GPU computation
                let config_gpu = SmagorinskyConfig {
                    use_gpu: true,
                    ..Default::default()
                };
                let mut les_gpu = SmagorinskyLES::new(size, size, 0.01, 0.01, config_gpu);
                les_gpu
                    .update(
                        &velocity_u,
                        &velocity_v,
                        &pressure,
                        1.0,
                        0.01,
                        0.001,
                        0.01,
                        0.01,
                    )
                    .unwrap();
                let gpu_viscosity = les_gpu.get_turbulent_viscosity_field();

                // Check accuracy (relative error should be small)
                let mut max_relative_error = 0.0f64;
                for i in 0..size {
                    for j in 0..size {
                        let cpu_val = cpu_viscosity[(i, j)];
                        let gpu_val = gpu_viscosity[(i, j)];

                        if cpu_val.abs() > 1e-12 {
                            let rel_error = ((cpu_val - gpu_val) / cpu_val).abs();
                            max_relative_error = max_relative_error.max(rel_error);
                        }
                    }
                }

                black_box(max_relative_error);
                assert!(
                    max_relative_error < 1e-6,
                    "GPU accuracy error too high: {}",
                    max_relative_error
                );
            }

            black_box(cpu_viscosity);
        });
    });

    group.finish();
}

/// Memory transfer benchmark
fn bench_memory_transfer(c: &mut Criterion) {
    let mut group = c.benchmark_group("GPU Memory Transfer");

    #[cfg(feature = "gpu")]
    for &size in [64, 128, 256].iter() {
        let data_size = size * size;
        let data: Vec<f32> = (0..data_size).map(|i| i as f32).collect();

        group.bench_function(format!("GPU Transfer {}x{}", size, size), |b| {
            if let Ok(mut gpu_compute) =
                cfd_core::compute::gpu::turbulence_compute::GpuTurbulenceCompute::new()
            {
                b.iter(|| {
                    // Create GPU buffer
                    let result_buffer = gpu_compute
                        .compute_smagorinsky_sgs(&data, &data, size, size, 0.01, 0.01, 0.1)
                        .unwrap();
                    let buffer = gpu_compute.read_buffer(&result_buffer).unwrap();
                    black_box(buffer);
                });
            }
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_smagorinsky_les,
    bench_des,
    bench_strain_rate_computation,
    bench_accuracy_validation,
    bench_memory_transfer
);
criterion_main!(benches);
