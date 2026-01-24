use cfd_math::iterators::StridedWindowIterator;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_strided_window(c: &mut Criterion) {
    let mut group = c.benchmark_group("strided_window");

    let window_size = 1000;
    let stride = 100;
    let data_len = 100_000;
    let data: Vec<f64> = (0..data_len).map(|i| i as f64).collect();

    group.bench_function("strided_window_large", |b| {
        b.iter(|| {
            let iter = StridedWindowIterator::new(data.iter().cloned(), window_size, stride);
            for window in iter {
                black_box(window);
            }
        });
    });

    group.finish();
}

criterion_group!(benches, bench_strided_window);
criterion_main!(benches);
