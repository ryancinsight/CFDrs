use cfd_math::linear_solver::preconditioners::multigrid::{falgout_coarsening, AlgebraicDistances};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use leto_ops::CsrMatrix;

fn poisson_matrix(n: usize) -> CsrMatrix<f64> {
    let size = n * n;
    let mut row_ptr = Vec::with_capacity(size + 1);
    let mut col_indices = Vec::new();
    let mut values = Vec::new();
    row_ptr.push(0);

    for i in 0..n {
        for j in 0..n {
            let row = i * n + j;
            if i > 0 {
                col_indices.push((i - 1) * n + j);
                values.push(-1.0);
            }
            if j > 0 {
                col_indices.push(i * n + j - 1);
                values.push(-1.0);
            }
            col_indices.push(row);
            values.push(4.0);
            if j < n - 1 {
                col_indices.push(i * n + j + 1);
                values.push(-1.0);
            }
            if i < n - 1 {
                col_indices.push((i + 1) * n + j);
                values.push(-1.0);
            }
            row_ptr.push(col_indices.len());
        }
    }

    CsrMatrix::from_parts(values, col_indices, row_ptr, size, size).unwrap()
}

fn bench_algebraic_distances(c: &mut Criterion) {
    let mut group = c.benchmark_group("algebraic_distances");

    // Test with grid size 20x20 (400 points)
    let n = 20;
    let size = n * n;

    let matrix = poisson_matrix(n);
    let threshold = 0.25;
    let coarsening = falgout_coarsening(&matrix, threshold).unwrap();

    group.bench_function(BenchmarkId::new("compute", size), |b| {
        b.iter(|| {
            black_box(AlgebraicDistances::compute(&coarsening, &matrix));
        })
    });

    group.finish();
}

criterion_group!(benches, bench_algebraic_distances);
criterion_main!(benches);
