use cfd_math::linear_solver::LinearOperator;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use leto::Array1;
use leto_ops::CsrMatrix;

fn bench_spmv(c: &mut Criterion) {
    let mut group = c.benchmark_group("spmv_leto_provider");

    for size in [100, 200, 500].iter() {
        let n = *size;
        // Create a simple tri-diagonal matrix
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
        let a = CsrMatrix::from_parts(values, col_indices, row_offsets, n, n).unwrap();
        let x = Array1::from_shape_vec([n], vec![1.0; n]).unwrap();
        let mut y = Array1::zeros([n]);

        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &_| {
            b.iter(|| {
                black_box(a.apply(black_box(&x), black_box(&mut y))).unwrap();
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_spmv);
criterion_main!(benches);
