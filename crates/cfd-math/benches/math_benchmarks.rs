use cfd_math::linear_solver::matrix_free::{LaplacianOperator2D, LinearOperator};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

fn bench_laplacian_cpu(c: &mut Criterion) {
    let mut group = c.benchmark_group("laplacian_cpu");
    for &n in &[64usize, 128usize, 256usize] {
        let nx = n;
        let ny = n;
        let dx = 1.0f64 / (nx as f64 - 1.0);
        let dy = 1.0f64 / (ny as f64 - 1.0);
        let op = LaplacianOperator2D::new(nx, ny, dx, dy);
        let mut field = vec![0.0f64; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * dx;
                let y = j as f64 * dy;
                field[j * nx + i] = x * (1.0 - x) + y * (1.0 - y);
            }
        }
        let mut out = vec![0.0f64; nx * ny];
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &_| {
            b.iter(|| {
                black_box(op.apply(&field, &mut out)).unwrap();
            });
        });
    }
    group.finish();
}

fn bench_cg_small_spd(c: &mut Criterion) {
    use cfd_math::linear_solver::{
        preconditioners::IdentityPreconditioner, ConjugateGradient, IterativeSolverConfig,
    };
    use nalgebra::DVector;
    use nalgebra_sparse::CsrMatrix;

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

    let mut group = c.benchmark_group("cg_spd");
    for &n in &[128usize, 256usize, 512usize] {
        let a = spd(n);
        let b = DVector::from_element(n, 1.0);
        let solver =
            ConjugateGradient::new(IterativeSolverConfig::new(1e-8).with_max_iterations(1000));
        let pre = IdentityPreconditioner;
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |bmk, &_| {
            bmk.iter(|| {
                let _x = black_box(solver.solve_preconditioned(&a, &b, &pre, None)).unwrap();
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_laplacian_cpu, bench_cg_small_spd);
criterion_main!(benches);
