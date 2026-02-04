use cfd_math::linear_solver::{ConjugateGradient, IdentityPreconditioner, IterativeSolverConfig};
use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use nalgebra::DVector;
use nalgebra_sparse::CsrMatrix;

fn bench_cg(c: &mut Criterion) {
    let mut group = c.benchmark_group("cg_solver");

    for size in [100, 500, 1000].iter() {
        let n = *size;
        // Create 1D Laplacian matrix (SPD)
        let mut row_offsets = Vec::with_capacity(n + 1);
        let mut col_indices = Vec::new();
        let mut values = Vec::new();
        row_offsets.push(0);
        for i in 0..n {
            if i > 0 {
                col_indices.push(i - 1);
                values.push(-1.0);
            }
            col_indices.push(i);
            values.push(2.0);
            if i + 1 < n {
                col_indices.push(i + 1);
                values.push(-1.0);
            }
            row_offsets.push(col_indices.len());
        }
        let a = CsrMatrix::try_from_csr_data(n, n, row_offsets, col_indices, values).unwrap();

        let b = DVector::from_element(n, 1.0);

        // Limit iterations to ensure we measure setup overhead significantly
        let config = IterativeSolverConfig {
            max_iterations: 5,
            tolerance: 1e-10,
            ..Default::default()
        };
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |bench, &_| {
            bench.iter_batched_ref(
                || DVector::zeros(n),
                |x_mut| {
                    let _ = solver.solve_preconditioned(
                        black_box(&a),
                        black_box(&b),
                        black_box(&precond),
                        black_box(x_mut),
                    );
                },
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

criterion_group!(benches, bench_cg);
criterion_main!(benches);
