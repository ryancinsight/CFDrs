use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use cfd_math::sparse::spmv;

fn bench_spmv(c: &mut Criterion) {
    let mut group = c.benchmark_group("spmv_serial");

    // We want a size that falls into the serial path.
    // Based on parallel_threshold logic, small matrices (e.g. 100) should be serial.
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
        let a = CsrMatrix::try_from_csr_data(n, n, row_offsets, col_indices, values).unwrap();
        let x = DVector::from_element(n, 1.0);
        let mut y = DVector::zeros(n);

        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &_| {
            b.iter(|| {
                spmv(black_box(&a), black_box(&x), black_box(&mut y));
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_spmv);
criterion_main!(benches);
