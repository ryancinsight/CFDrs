use cfd_math::simd::SimdOps;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_simd_add_f64(c: &mut Criterion) {
    let ops = SimdOps::new();
    let size = 100_000;
    let a = vec![1.0; size];
    let b = vec![2.0; size];
    let mut result = vec![0.0; size];

    let mut group = c.benchmark_group("simd_f64_ops");

    group.bench_function("add_f64", |bencher| {
        bencher.iter(|| {
            ops.add_f64(black_box(&a), black_box(&b), black_box(&mut result))
                .unwrap();
        })
    });

    group.bench_function("dot_f64", |bencher| {
        bencher.iter(|| {
            ops.dot_f64(black_box(&a), black_box(&b)).unwrap();
        })
    });

    group.finish();
}

criterion_group!(benches, bench_simd_add_f64);
criterion_main!(benches);
