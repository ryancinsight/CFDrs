use criterion::{black_box, criterion_group, criterion_main, Criterion};
use cfd_math::sparse::SparseMatrixBuilder;
use rand::Rng;

fn bench_sparse_builder(c: &mut Criterion) {
    let mut group = c.benchmark_group("SparseMatrixBuilder");
    let rows = 1000;
    let cols = 1000;
    let num_entries = 50000;

    // Generate random entries with duplicates
    let mut rng = rand::thread_rng();
    let mut entries = Vec::with_capacity(num_entries);
    for _ in 0..num_entries {
        let row = rng.gen_range(0..rows);
        let col = rng.gen_range(0..cols);
        let val: f64 = rng.gen();
        entries.push((row, col, val));
    }

    group.bench_function("build_sequential", |b| {
        b.iter(|| {
            let mut builder = SparseMatrixBuilder::<f64>::new(rows, cols);

            for &(r, c, v) in &entries {
                builder.add_entry(r, c, v).unwrap();
            }

            let _ = black_box(builder.build().unwrap());
        })
    });

    group.bench_function("build_parallel", |b| {
        b.iter(|| {
            let mut builder = SparseMatrixBuilder::<f64>::new(rows, cols);

            for &(r, c, v) in &entries {
                builder.add_entry(r, c, v).unwrap();
            }

            let _ = black_box(builder.build_parallel().unwrap());
        })
    });

    group.finish();
}

criterion_group!(benches, bench_sparse_builder);
criterion_main!(benches);
