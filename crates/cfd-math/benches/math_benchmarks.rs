use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use cfd_math::{
    linear_solver::{LinearSolver, ConjugateGradient, BiCGSTAB},
    sparse::{SparseMatrix, SparseMatrixBuilder},
    interpolation::{Interpolation, LinearInterpolation, CubicSplineInterpolation},
    integration::{Quadrature, GaussQuadrature},
    differentiation::FiniteDifference,
    iterators::MathIteratorExt,

};
use nalgebra::DVector;

fn benchmark_linear_solvers(c: &mut Criterion) {
    let mut group = c.benchmark_group("linear_solvers");
    
    for size in [100, 500, 1000].iter() {
        let (matrix, rhs) = create_test_linear_system(*size);
        
        let cg_solver = ConjugateGradient::default();
        let bicgstab_solver = BiCGSTAB::default();
        
        group.bench_with_input(
            BenchmarkId::new("conjugate_gradient", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(cg_solver.solve(&matrix, &rhs, None).unwrap())
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("bicgstab", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(bicgstab_solver.solve(&matrix, &rhs, None).unwrap())
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_sparse_matrix_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("sparse_matrix");
    
    for size in [1000, 5000, 10000].iter() {
        let sparse_matrix = create_test_sparse_matrix(*size);
        let vector = DVector::from_fn(*size, |i, _| (i as f64).sin());
        
        group.bench_with_input(
            BenchmarkId::new("matrix_vector_multiply", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(&sparse_matrix * &vector)
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("matrix_construction", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(create_test_sparse_matrix(*size))
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_interpolation(c: &mut Criterion) {
    let mut group = c.benchmark_group("interpolation");
    
    for size in [100, 1000, 10000].iter() {
        let x_data: Vec<f64> = (0..*size).map(|i| i as f64 / *size as f64).collect();
        let y_data: Vec<f64> = x_data.iter().map(|&x| x.sin()).collect();
        let query_points: Vec<f64> = (0..(*size/2)).map(|i| (i as f64 + 0.5) / *size as f64).collect();
        
        let linear_interp = LinearInterpolation::new(x_data.clone(), y_data.clone()).unwrap();
        let cubic_interp = CubicSplineInterpolation::new(x_data, y_data).unwrap();
        
        group.bench_with_input(
            BenchmarkId::new("linear_interpolation", size),
            size,
            |b, _| {
                b.iter(|| {
                    for &x in &query_points {
                        black_box(linear_interp.interpolate(x));
                    }
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("cubic_spline_interpolation", size),
            size,
            |b, _| {
                b.iter(|| {
                    for &x in &query_points {
                        black_box(cubic_interp.interpolate(x));
                    }
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_integration(c: &mut Criterion) {
    let mut group = c.benchmark_group("integration");
    
    let test_function = |x: f64| x.sin() * x.exp();
    
    for order in [2, 4, 8, 16].iter() {
        let gauss_quad = GaussQuadrature::new(*order).unwrap();
        
        group.bench_with_input(
            BenchmarkId::new("gauss_quadrature", order),
            order,
            |b, _| {
                b.iter(|| {
                    black_box(gauss_quad.integrate(test_function, 0.0, 1.0))
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_differentiation(c: &mut Criterion) {
    let mut group = c.benchmark_group("differentiation");
    
    for size in [1000, 10000, 100000].iter() {
        let field: Vec<f64> = (0..*size).map(|i| (i as f64 * 0.01).sin()).collect();
        let dx = 0.01;
        
        let finite_diff = FiniteDifference::central(dx);
        
        group.bench_with_input(
            BenchmarkId::new("finite_difference_first_derivative", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(finite_diff.first_derivative(&field).unwrap())
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("finite_difference_second_derivative", size),
            size,
            |b, _| {
                b.iter(|| {
                    black_box(finite_diff.second_derivative(&field).unwrap())
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_vectorized_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("vectorized_operations");
    
    for size in [10000, 100000, 1000000].iter() {
        let vec1: Vec<f64> = (0..*size).map(|i| i as f64).collect();
        let vec2: Vec<f64> = (0..*size).map(|i| (i as f64).sin()).collect();
        
        group.bench_with_input(
            BenchmarkId::new("vector_add", size),
            size,
            |b, _| {
                b.iter(|| {
                    let result: Vec<f64> = vec1.iter()
                        .zip(vec2.iter())
                        .map(|(a, b)| a + b)
                        .collect();
                    black_box(result)
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("dot_product", size),
            size,
            |b, _| {
                b.iter(|| {
                    let result: f64 = vec1.iter()
                        .zip(vec2.iter())
                        .map(|(a, b)| a * b)
                        .sum();
                    black_box(result)
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("l2_norm", size),
            size,
            |b, _| {
                b.iter(|| {
                    let result = vec1.iter().copied().l2_norm();
                    black_box(result)
                })
            },
        );
    }
    
    group.finish();
}

fn benchmark_iterator_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("iterator_operations");
    
    for size in [10000, 100000, 1000000].iter() {
        let data: Vec<f64> = (0..*size).map(|i| (i as f64).sin()).collect();
        
        group.bench_with_input(
            BenchmarkId::new("windowed_operations", size),
            size,
            |b, _| {
                b.iter(|| {
                    let result: Vec<f64> = data.windows(3)
                        .map(|w| w[2] - w[0])
                        .collect();
                    black_box(result)
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("cumulative_sum", size),
            size,
            |b, _| {
                b.iter(|| {
                    let result: Vec<f64> = data.iter()
                        .scan(0.0, |acc, &x| {
                            *acc += x;
                            Some(*acc)
                        })
                        .collect();
                    black_box(result)
                })
            },
        );
    }
    
    group.finish();
}

// Helper functions
fn create_test_linear_system(size: usize) -> (SparseMatrix<f64>, DVector<f64>) {
    let mut builder = SparseMatrixBuilder::new(size, size);
    
    // Create a symmetric positive definite matrix (tridiagonal)
    for i in 0..size {
        let _ = builder.add_entry(i, i, 2.0);
        if i > 0 {
            let _ = builder.add_entry(i, i-1, -1.0);
        }
        if i < size - 1 {
            let _ = builder.add_entry(i, i+1, -1.0);
        }
    }
    
    let rhs = DVector::from_fn(size, |i, _| (i as f64).sin());
    let matrix = builder.build().expect("Failed to build sparse matrix");
    
    (matrix, rhs)
}

fn create_test_sparse_matrix(size: usize) -> SparseMatrix<f64> {
    let mut builder = SparseMatrixBuilder::new(size, size);
    
    // Create a 5-point stencil pattern
    for i in 0..size {
        let _ = builder.add_entry(i, i, 4.0);
        if i > 0 {
            let _ = builder.add_entry(i, i-1, -1.0);
        }
        if i < size - 1 {
            let _ = builder.add_entry(i, i+1, -1.0);
        }
        if i >= 10 {
            let _ = builder.add_entry(i, i-10, -1.0);
        }
        if i < size - 10 {
            let _ = builder.add_entry(i, i+10, -1.0);
        }
    }
    
    builder.build().expect("Failed to build sparse matrix")
}

criterion_group!(
    benches,
    benchmark_linear_solvers,
    benchmark_sparse_matrix_operations,
    benchmark_interpolation,
    benchmark_integration,
    benchmark_differentiation,
    benchmark_vectorized_operations,
    benchmark_iterator_operations
);
criterion_main!(benches);
