#![expect(
    missing_docs,
    reason = "criterion_group generates public benchmark entry points"
)]

use cfd_math::iterative::{ConjugateGradient, IdentityPreconditioner, IterativeSolverConfig};
use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use leto::Array1;
use leto_ops::CsrMatrix;

fn bench_cg(c: &mut Criterion) {
    let mut group = c.benchmark_group("cg_solver");

    for &size in &[100, 500, 1000] {
        let n = size;
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
        let a = CsrMatrix::from_parts(values, col_indices, row_offsets, n, n)
            .expect("invariant: CG benchmark matrix is structurally valid");

        let b = Array1::from_elem([n], 1.0);

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
                || Array1::zeros([n]),
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

/// Full-convergence regime (contrast with `bench_cg`'s capped-iteration
/// setup-overhead regime): diagonally dominant tridiagonal SPD systems
/// solved to 1e-8, absorbed from the retired `math_benchmarks` binary.
fn bench_cg_convergence(c: &mut Criterion) {
    let mut group = c.benchmark_group("cg_spd");
    for &n in &[128usize, 256usize, 512usize] {
        let mut row_offsets = Vec::with_capacity(n + 1);
        let mut col_indices = Vec::new();
        let mut values = Vec::new();
        row_offsets.push(0);
        for i in 0..n {
            if i > 0 {
                col_indices.push(i - 1);
                values.push(1.0);
            }
            col_indices.push(i);
            values.push(4.0);
            if i + 1 < n {
                col_indices.push(i + 1);
                values.push(1.0);
            }
            row_offsets.push(col_indices.len());
        }
        let a = CsrMatrix::from_parts(values, col_indices, row_offsets, n, n)
            .expect("invariant: CG convergence matrix is structurally valid");
        let b = Array1::from_elem([n], 1.0);
        let solver =
            ConjugateGradient::new(IterativeSolverConfig::new(1e-8).with_max_iterations(1000));
        let precond = IdentityPreconditioner;
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |bmk, &_| {
            bmk.iter_batched_ref(
                || Array1::zeros([n]),
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

criterion_group!(benches, bench_cg, bench_cg_convergence);
criterion_main!(benches);
