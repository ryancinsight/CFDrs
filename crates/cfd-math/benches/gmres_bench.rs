use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use cfd_math::linear_solver::{
    preconditioners::IdentityPreconditioner, GMRES, IterativeSolverConfig,
};
use nalgebra::DVector;
use nalgebra_sparse::CsrMatrix;

fn bench_gmres_spd(c: &mut Criterion) {
    fn spd(n: usize) -> CsrMatrix<f64> {
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
        CsrMatrix::try_from_csr_data(n, n, row_offsets, col_indices, values).unwrap()
    }

    let mut group = c.benchmark_group("gmres_spd");
    for &n in &[128usize, 256usize] {
        let a = spd(n);
        let b = DVector::from_element(n, 1.0);
        // Use fewer iterations to make bench faster, but enough to trigger restarts if m is small
        // Here m=30 (default). n=128.
        let solver =
            GMRES::new(IterativeSolverConfig::new(1e-8).with_max_iterations(100), 30);
        let pre = IdentityPreconditioner;
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |bmk, &_| {
            bmk.iter(|| {
                let mut x = DVector::zeros(n);
                black_box(solver.solve_preconditioned(&a, &b, &pre, &mut x)).unwrap();
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_gmres_spd);
criterion_main!(benches);
