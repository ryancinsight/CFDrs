//! Core CFD operation benchmarks
//!
//! This module contains benchmarks for fundamental CFD operations including
//! grid management, linear solvers, turbulence modeling, and GPU acceleration.

use crate::PerformanceMetrics;
use cfd_validation::benchmarking::BenchmarkConfig;
use cfd_2d::grid::{Grid2D, StructuredGrid2D};
use cfd_math::linear_solver::{ConjugateGradient, IterativeLinearSolver,
};
use cfd_math::sparse::SparseMatrixBuilder;
use criterion::{black_box, BenchmarkId, Criterion, Throughput};
use nalgebra::DVector;
use std::time::Instant;

/// Benchmark grid operations across different problem sizes
pub fn benchmark_grid_operations(c: &mut Criterion, config: &BenchmarkConfig) {
    let mut group = c.benchmark_group("grid_operations");

    for &size in &config.problem_sizes {
        // Grid creation benchmark
        group.bench_with_input(
            BenchmarkId::new("grid_creation", size),
            &size,
            |b, &size| {
                b.iter(|| {
                    let grid =
                        StructuredGrid2D::<f64>::new(size, size, 0.0, 1.0, 0.0, 1.0).unwrap();
                    black_box(grid);
                });
            },
        );

        // Grid indexing benchmark
        let grid = StructuredGrid2D::<f64>::new(size, size, 0.0, 1.0, 0.0, 1.0).unwrap();
        let total_elements = size * size;

        group.throughput(Throughput::Elements(total_elements as u64));
        group.bench_with_input(BenchmarkId::new("grid_indexing", size), &size, |b, _| {
            b.iter(|| {
                let mut sum = 0.0;
                for i in 0..grid.nx() {
                    for j in 0..grid.ny() {
                        sum += black_box(i as f64 * j as f64);
                    }
                }
                black_box(sum);
            });
        });

        // Grid interpolation benchmark
        group.bench_with_input(
            BenchmarkId::new("grid_interpolation", size),
            &size,
            |b, _| {
                b.iter(|| {
                    let mut interpolated = 0.0;
                    for i in 0..size.saturating_sub(1) {
                        for j in 0..size.saturating_sub(1) {
                            let x = i as f64 / size as f64;
                            let y = j as f64 / size as f64;
                            // Bilinear interpolation simulation
                            interpolated += black_box(x * y + (1.0 - x) * (1.0 - y));
                        }
                    }
                    black_box(interpolated);
                });
            },
        );
    }

    group.finish();
}

/// Benchmark linear solver operations
pub fn benchmark_linear_solvers(c: &mut Criterion, config: &BenchmarkConfig) {
    let mut group = c.benchmark_group("solver_operations");

    for &size in &config.problem_sizes {
        // Sparse matrix assembly
        group.bench_with_input(
            BenchmarkId::new("sparse_matrix_assembly", size),
            &size,
            |b, &size| {
                b.iter(|| {
                    let mut builder = SparseMatrixBuilder::new(size, size);
                    // Create tridiagonal matrix
                    for i in 0..size {
                        builder.add_entry(i, i, 2.0).unwrap();
                        if i > 0 {
                            builder.add_entry(i, i - 1, -1.0).unwrap();
                        }
                        if i < size - 1 {
                            builder.add_entry(i, i + 1, -1.0).unwrap();
                        }
                    }
                    let matrix = builder.build().unwrap();
                    black_box(matrix);
                });
            },
        );

        // Matrix-vector multiplication
        let mut builder = SparseMatrixBuilder::new(size, size);
        for i in 0..size {
            builder.add_entry(i, i, 2.0).unwrap();
            if i > 0 {
                builder.add_entry(i, i - 1, -1.0).unwrap();
            }
            if i < size - 1 {
                builder.add_entry(i, i + 1, -1.0).unwrap();
            }
        }
        let matrix = builder.build().unwrap();
        let x = DVector::from_element(size, 1.0);
        let mut y = DVector::zeros(size);

        group.throughput(Throughput::Elements(size as u64));
        group.bench_with_input(BenchmarkId::new("spmv", size), &size, |b, _| {
            b.iter(|| {
                y.copy_from(&(&matrix * &x));
                black_box(&y);
            });
        });

        // Conjugate gradient solver
        let b_vec = DVector::from_element(size, 1.0);
        let mut x = DVector::from_element(size, 0.0);
        group.bench_with_input(
            BenchmarkId::new("conjugate_gradient", size),
            &size,
            |b, _| {
                use cfd_math::linear_solver::IdentityPreconditioner;
                let solver = ConjugateGradient::<f64>::default();
                b.iter(|| {
                    x.fill(0.0);
                    let solution = solver
                        .solve::<IdentityPreconditioner>(&matrix, &b_vec, &mut x, None)
                        .unwrap();
                    black_box(solution);
                });
            },
        );
    }

    group.finish();
}

/// Benchmark turbulence modeling operations
pub fn benchmark_turbulence_operations(c: &mut Criterion, config: &BenchmarkConfig) {
    let mut group = c.benchmark_group("turbulence_operations");

    for &size in &config.problem_sizes {
        // Strain rate tensor computation
        group.bench_with_input(
            BenchmarkId::new("strain_rate_computation", size),
            &size,
            |b, _| {
                b.iter(|| {
                    let mut strain_rate = vec![0.0; size * size];
                    for i in 1..size.saturating_sub(1) {
                        for j in 1..size.saturating_sub(1) {
                            let idx = i * size + j;
                            // Simulate velocity field derivatives
                            let du_dx = 0.1 * (i as f64).sin();
                            let du_dy = 0.1 * (j as f64).cos();
                            let dv_dx = 0.1 * (j as f64).sin();
                            let dv_dy = 0.1 * (i as f64).cos();

                            // Strain rate magnitude
                            strain_rate[idx] = black_box(
                                (2.0 * du_dx * du_dx
                                    + 2.0 * dv_dy * dv_dy
                                    + (du_dy + dv_dx).powi(2))
                                .sqrt(),
                            );
                        }
                    }
                    black_box(strain_rate);
                });
            },
        );

        // Turbulent viscosity computation (k-epsilon model simulation)
        group.bench_with_input(
            BenchmarkId::new("turbulent_viscosity_k_epsilon", size),
            &size,
            |b, _| {
                b.iter(|| {
                    let mut viscosity = vec![0.0; size * size];
                    for i in 0..size {
                        for j in 0..size {
                            let idx = i * size + j;
                            // Simulate k-epsilon model calculations
                            let k = 0.1 + 0.01 * (i as f64 * j as f64).sin();
                            let epsilon = 0.01 + 0.001 * (i as f64 + j as f64).cos();
                            let c_mu = 0.09;

                            viscosity[idx] = black_box(c_mu * k * k / (epsilon + 1e-12));
                        }
                    }
                    black_box(viscosity);
                });
            },
        );

        // Wall distance computation (for wall functions)
        group.bench_with_input(
            BenchmarkId::new("wall_distance_computation", size),
            &size,
            |b, _| {
                b.iter(|| {
                    let mut wall_distance = vec![0.0; size * size];
                    // Assume walls at j=0 and j=size-1
                    for i in 0..size {
                        for j in 0..size {
                            let idx = i * size + j;
                            let dist_to_bottom = j as f64 * 0.01;
                            let dist_to_top = (size - 1 - j) as f64 * 0.01;
                            wall_distance[idx] = black_box(dist_to_bottom.min(dist_to_top));
                        }
                    }
                    black_box(wall_distance);
                });
            },
        );
    }

    group.finish();
}

/// Benchmark GPU-accelerated operations
pub fn benchmark_gpu_operations(c: &mut Criterion, config: &BenchmarkConfig) {
    let mut group = c.benchmark_group("gpu_operations");

    // Only run GPU benchmarks if GPU feature is enabled
    #[cfg(feature = "gpu")]
    {
        for &size in &config.problem_sizes {
            // GPU memory allocation and transfer
            group.bench_with_input(
                BenchmarkId::new("gpu_memory_transfer", size),
                &size,
                |b, &size| {
                    b.iter(|| {
                        // Simulate GPU buffer operations
                        let data_size = size * size;
                        let host_data: Vec<f32> = (0..data_size).map(|i| i as f32).collect();

                        // Simulate GPU buffer creation and data transfer
                        black_box(host_data.len());
                        // In real implementation, this would create GPU buffers
                    });
                },
            );

            // GPU kernel execution simulation
            group.throughput(Throughput::Elements((size * size) as u64));
            group.bench_with_input(
                BenchmarkId::new("gpu_kernel_execution", size),
                &size,
                |b, _| {
                    b.iter(|| {
                        // Simulate GPU kernel for CFD operations
                        let mut result = vec![0.0f32; size * size];
                        for i in 0..size {
                            for j in 0..size {
                                let idx = i * size + j;
                                // Simulate a simple CFD kernel (e.g., Laplacian)
                                result[idx] = black_box(
                                    (i as f32).sin()
                                        + (j as f32).cos()
                                        + 0.25
                                            * ((i.saturating_sub(1) as f32).sin()
                                                + (i.saturating_add(1).min(size - 1) as f32).sin()
                                                + (j.saturating_sub(1) as f32).cos()
                                                + (j.saturating_add(1).min(size - 1) as f32).cos()),
                                );
                            }
                        }
                        black_box(result);
                    });
                },
            );
        }
    }

    #[cfg(not(feature = "gpu"))]
    {
        group.bench_function("gpu_unavailable", |b| {
            b.iter(|| {
                // Placeholder benchmark when GPU is not available
                black_box(());
            });
        });
    }

    group.finish();
}

/// Collect performance metrics for a benchmark operation
pub fn collect_performance_metrics<F, T>(
    operation_name: &str,
    problem_size: usize,
    operation: F,
) -> PerformanceMetrics
where
    F: FnOnce() -> T,
{
    let start = Instant::now();
    let _result = operation();
    let duration = start.elapsed();

    PerformanceMetrics::new(operation_name.to_string(), problem_size)
        .with_timing(duration)
        .with_throughput(problem_size as f64 / duration.as_secs_f64())
}

pub fn bench_core_operations(c: &mut Criterion) {
    let config = BenchmarkConfig::default();
    benchmark_grid_operations(c, &config);
    benchmark_linear_solvers(c, &config);
}
