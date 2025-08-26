//! Performance benchmarks for critical CFD operations

use cfd_2d::grid::{Grid2D, StructuredGrid2D};
use cfd_math::sparse::SparseMatrixBuilder;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nalgebra::DVector;

fn benchmark_grid_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("grid_operations");

    for size in [10, 50, 100, 200].iter() {
        group.bench_with_input(BenchmarkId::new("grid_creation", size), size, |b, &size| {
            b.iter(|| {
                let grid = StructuredGrid2D::<f64>::new(size, size, 0.0, 1.0, 0.0, 1.0).unwrap();
                black_box(grid);
            });
        });

        let grid = StructuredGrid2D::<f64>::new(*size, *size, 0.0, 1.0, 0.0, 1.0).unwrap();

        group.bench_with_input(BenchmarkId::new("grid_indexing", size), size, |b, _| {
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
    }

    group.finish();
}

fn benchmark_sparse_matrix_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("sparse_matrix");

    for size in [100, 500, 1000, 2000].iter() {
        // Benchmark matrix assembly
        group.bench_with_input(BenchmarkId::new("assembly", size), size, |b, &size| {
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
        });

        // Create matrix for matvec benchmark
        let mut builder = SparseMatrixBuilder::new(*size, *size);
        for i in 0..*size {
            builder.add_entry(i, i, 2.0).unwrap();
            if i > 0 {
                builder.add_entry(i, i - 1, -1.0).unwrap();
            }
            if i < *size - 1 {
                builder.add_entry(i, i + 1, -1.0).unwrap();
            }
        }
        let matrix = builder.build().unwrap();
        let x = DVector::from_element(*size, 1.0);
        let mut y = DVector::zeros(*size);

        // Benchmark matrix-vector multiplication
        group.bench_with_input(BenchmarkId::new("matvec", size), size, |b, _| {
            b.iter(|| {
                // Replace custom matvec with standard multiplication
                y.copy_from(&(&matrix * &x));
                black_box(&y);
            });
        });
    }

    group.finish();
}

fn benchmark_fluid_calculations(c: &mut Criterion) {
    use cfd_core::Fluid;

    let mut group = c.benchmark_group("fluid_calculations");

    let fluid = Fluid::constant_viscosity("water", 1000.0, 0.001);

    group.bench_function("reynolds_number", |b| {
        b.iter(|| {
            let velocity = black_box(1.0);
            let length = black_box(0.1);
            let re = fluid.reynolds_number(velocity, length);
            black_box(re);
        });
    });

    group.bench_function("kinematic_viscosity", |b| {
        b.iter(|| {
            let nu = fluid.kinematic_viscosity();
            black_box(nu);
        });
    });

    group.finish();
}

fn benchmark_analytical_solutions(c: &mut Criterion) {
    use cfd_validation::analytical::{AnalyticalSolution, PoiseuilleFlow};

    let mut group = c.benchmark_group("analytical_solutions");

    let poiseuille = PoiseuilleFlow::create(
        2.5,    // u_max
        0.05,   // characteristic length (half-width), was 0.1 width
        -100.0, // pressure gradient
        0.001,  // viscosity
        cfd_validation::analytical::PoiseuilleGeometry::Plates,
    );

    group.bench_function("poiseuille_velocity", |b| {
        b.iter(|| {
            let x = black_box(0.5);
            let y = black_box(0.05);
            let z = black_box(0.0);
            let t = black_box(0.0);
            let v = poiseuille.velocity(x, y, z, t);
            black_box(v);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    benchmark_grid_operations,
    benchmark_sparse_matrix_operations,
    benchmark_fluid_calculations,
    benchmark_analytical_solutions
);
criterion_main!(benches);
