use cfd_math::linear_solver::DirectSparseSolver;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use leto::Array1;
use leto_ops::{lu_decompose, qr_decompose, CooMatrix, CsrMatrix as LetoCsrMatrix};

fn tridiagonal_spd_matrix(n: usize) -> (LetoCsrMatrix<f64>, Array1<f64>) {
    let mut coo = CooMatrix::new(n, n);
    for i in 0..n {
        coo.push(i, i, 2.0);
        if i > 0 {
            coo.push(i, i - 1, -1.0);
        }
        if i + 1 < n {
            coo.push(i, i + 1, -1.0);
        }
    }
    let rhs = Array1::from_elem([n], 1.0_f64);
    (coo.to_csr(), rhs)
}

fn solve_dense_cutoff(matrix: &LetoCsrMatrix<f64>, rhs: &Array1<f64>) -> Array1<f64> {
    let dense = matrix.to_dense();

    if let Ok(lu) = lu_decompose(&dense.view()) {
        if let Ok(x) = lu.solve(&rhs.view()) {
            return x;
        }
    }

    qr_decompose(&dense.view())
        .and_then(|qr| qr.solve_least_squares(&rhs.view()))
        .expect("dense QR fallback must solve SPD benchmark matrix")
}

fn solve_sparse_spd_direct(matrix: &LetoCsrMatrix<f64>, rhs: &Array1<f64>) -> Array1<f64> {
    DirectSparseSolver {
        max_size: 256,
        ordering: 0,
        pivot_tolerance: 1e-12,
    }
    .solve(matrix, rhs)
    .expect("sparse direct SPD solver must solve benchmark matrix")
}

fn bench_small_system_direct_solvers(c: &mut Criterion) {
    let mut group = c.benchmark_group("small_system_direct_solver_cutoff");
    for &n in &[32usize, 64, 128, 256] {
        let (matrix, rhs) = tridiagonal_spd_matrix(n);
        group.bench_with_input(BenchmarkId::new("dense_cutoff", n), &n, |b, _| {
            b.iter(|| black_box(solve_dense_cutoff(black_box(&matrix), black_box(&rhs))));
        });
        group.bench_with_input(BenchmarkId::new("sparse_spd_direct", n), &n, |b, _| {
            b.iter(|| black_box(solve_sparse_spd_direct(black_box(&matrix), black_box(&rhs))));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_small_system_direct_solvers);
criterion_main!(benches);
