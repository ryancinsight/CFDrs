//! Performance benchmarks for CFD Suite
//!
//! These benchmarks measure actual performance of core operations

use cfd_core::physics::fluid_dynamics::{FlowField, FlowOperations};
use cfd_math::linear_solver::{ConjugateGradient, IterativeLinearSolver};
use cfd_math::sparse::SparseMatrixBuilder;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use nalgebra::DVector;

/// Benchmark flow field operations at different scales
fn benchmark_flow_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("flow_operations");

    for size in [16, 32, 64].iter() {
        let flow_field = FlowField::<f64>::new(*size, *size, *size);
        let total_points = size * size * size;

        group.throughput(Throughput::Elements(total_points as u64));

        group.bench_with_input(BenchmarkId::new("divergence", size), size, |b, _| {
            b.iter(|| black_box(FlowOperations::divergence(&flow_field.velocity)))
        });

        group.bench_with_input(BenchmarkId::new("vorticity", size), size, |b, _| {
            b.iter(|| black_box(FlowOperations::vorticity(&flow_field.velocity)))
        });

        group.bench_with_input(BenchmarkId::new("kinetic_energy", size), size, |b, _| {
            b.iter(|| black_box(FlowOperations::kinetic_energy(&flow_field.velocity)))
        });
    }

    group.finish();
}

/// Benchmark linear solver performance
fn benchmark_linear_solver(c: &mut Criterion) {
    let mut group = c.benchmark_group("linear_solver");

    for size in [10, 50, 100, 200].iter() {
        // Create a tridiagonal system (common in CFD)
        let mut builder = SparseMatrixBuilder::new(*size, *size);
        for i in 0..*size {
            builder.add_entry(i, i, 2.0).unwrap();
            if i > 0 {
                builder.add_entry(i, i - 1, -1.0).unwrap();
            }
            if i < size - 1 {
                builder.add_entry(i, i + 1, -1.0).unwrap();
            }
        }
        let matrix = builder.build().unwrap();
        let b = DVector::from_element(*size, 1.0);

        group.throughput(Throughput::Elements(*size as u64));

        group.bench_with_input(
            BenchmarkId::new("conjugate_gradient", size),
            size,
            |bench, _| {
                let solver = ConjugateGradient::<f64>::default();
                let mut x = DVector::from_element(*size, 0.0);
                bench.iter(|| {
                    black_box(
                        solver
                            .solve(
                                &matrix,
                                &b,
                                &mut x,
                                None::<&cfd_math::linear_solver::IdentityPreconditioner>,
                            )
                            .unwrap(),
                    )
                })
            },
        );
    }

    group.finish();
}

/// Benchmark memory allocation patterns
fn benchmark_memory_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_patterns");

    // Benchmark different allocation patterns
    group.bench_function("vec_allocation_1mb", |b| {
        b.iter(|| {
            let v: Vec<f64> = black_box(vec![0.0; 131072]); // 1MB
            v
        })
    });

    group.bench_function("vec_resize_1mb", |b| {
        let mut v = Vec::with_capacity(131072);
        b.iter(|| {
            v.clear();
            v.resize(131072, 0.0);
            black_box(&v);
        })
    });

    group.bench_function("matrix_allocation_1000x1000", |b| {
        b.iter(|| {
            let m = nalgebra::DMatrix::<f64>::zeros(1000, 1000);
            black_box(m)
        })
    });

    group.finish();
}

/// Benchmark scaling behavior
fn benchmark_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("scaling");

    // Test how performance scales with problem size
    let sizes = vec![10, 20, 40, 80];

    for size in sizes {
        let flow_field = FlowField::<f64>::new(size, size, size);
        let total_ops = (size * size * size) as u64;

        group.throughput(Throughput::Elements(total_ops));

        group.bench_with_input(BenchmarkId::new("full_computation", size), &size, |b, _| {
            b.iter(|| {
                // Simulate a full CFD step
                let div = FlowOperations::divergence(&flow_field.velocity);
                let vort = FlowOperations::vorticity(&flow_field.velocity);
                let ke = FlowOperations::kinetic_energy(&flow_field.velocity);
                black_box((div, vort, ke))
            })
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_flow_operations,
    benchmark_linear_solver,
    benchmark_memory_patterns,
    benchmark_scaling
);

criterion_main!(benches);
