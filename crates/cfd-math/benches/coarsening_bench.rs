use cfd_math::linear_solver::preconditioners::multigrid::falgout_coarsening;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nalgebra_sparse::{CooMatrix, CsrMatrix};

fn bench_falgout_coarsening(c: &mut Criterion) {
    let mut group = c.benchmark_group("falgout_coarsening");
    // Test with increasing grid sizes
    for &n in &[50, 100] {
        let size = n * n;

        // Construct 2D Poisson matrix
        let mut coo = CooMatrix::new(size, size);
        for i in 0..n {
            for j in 0..n {
                let row = i * n + j;
                coo.push(row, row, 4.0);
                if i > 0 {
                    coo.push(row, (i - 1) * n + j, -1.0);
                }
                if i < n - 1 {
                    coo.push(row, (i + 1) * n + j, -1.0);
                }
                if j > 0 {
                    coo.push(row, i * n + j - 1, -1.0);
                }
                if j < n - 1 {
                    coo.push(row, i * n + j + 1, -1.0);
                }
            }
        }
        let matrix = CsrMatrix::from(&coo);
        let threshold = 0.25;

        group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, &_| {
            b.iter(|| {
                black_box(falgout_coarsening(&matrix, threshold).unwrap());
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_falgout_coarsening);
criterion_main!(benches);
