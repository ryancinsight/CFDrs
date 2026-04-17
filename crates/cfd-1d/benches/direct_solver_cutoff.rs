use cfd_math::linear_solver::DirectSparseSolver;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nalgebra::{DMatrix, DVector};
use nalgebra_sparse::{coo::CooMatrix, CsrMatrix};

fn tridiagonal_spd_matrix(n: usize) -> (CsrMatrix<f64>, DVector<f64>) {
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
    let rhs = DVector::from_element(n, 1.0);
    (CsrMatrix::from(&coo), rhs)
}

fn solve_dense_cutoff(matrix: &CsrMatrix<f64>, rhs: &DVector<f64>) -> DVector<f64> {
    let mut dense = DMatrix::zeros(matrix.nrows(), matrix.ncols());
    for row_idx in 0..matrix.nrows() {
        let row = matrix.row(row_idx);
        for (col_idx, value) in row.col_indices().iter().zip(row.values()) {
            dense[(row_idx, *col_idx)] = *value;
        }
    }
    if let Some(cholesky) = dense.clone().cholesky() {
        return cholesky.solve(rhs);
    }
    if let Some(solution) = dense.clone().lu().solve(rhs) {
        return solution;
    }
    dense
        .qr()
        .solve(rhs)
        .expect("dense QR fallback must solve SPD benchmark matrix")
}

fn solve_sparse_spd_direct(matrix: &CsrMatrix<f64>, rhs: &DVector<f64>) -> DVector<f64> {
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
