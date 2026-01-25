use cfd_math::simd::SwarOps;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_swar_add_f64(c: &mut Criterion) {
    let swar = SwarOps::new();
    let size = 100_000;
    let a = vec![1.0; size];
    let b = vec![2.0; size];
    let mut result = vec![0.0; size];

    let mut group = c.benchmark_group("swar_f64_ops");

    group.bench_function("add_f64_arrays", |bencher| {
        bencher.iter(|| {
            swar.add_f64_arrays(black_box(&a), black_box(&b), black_box(&mut result))
                .unwrap();
        })
    });

    group.finish();
}

criterion_group!(benches, bench_swar_add_f64);
criterion_main!(benches);
