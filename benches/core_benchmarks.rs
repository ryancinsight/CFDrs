use cfd_core::domain::Domain;
use cfd_math::differentiation::finite_difference::FiniteDifference;
use cfd_math::differentiation::schemes::FiniteDifferenceScheme;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

fn benchmark_finite_difference(c: &mut Criterion) {
    let mut group = c.benchmark_group("finite_difference");

    for size in [100, 1000, 10000].iter() {
        let values: Vec<f64> = (0..*size).map(|i| (i as f64).sin()).collect();
        let fd = FiniteDifference::new(FiniteDifferenceScheme::Central, 0.01);

        group.bench_with_input(BenchmarkId::new("first_derivative", size), size, |b, _| {
            b.iter(|| fd.first_derivative(black_box(&values)));
        });

        group.bench_with_input(BenchmarkId::new("second_derivative", size), size, |b, _| {
            b.iter(|| fd.second_derivative(black_box(&values)));
        });
    }
    group.finish();
}

fn benchmark_domain_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("domain");

    let domain = Domain::<f64>::new(
        nalgebra::Point3::new(0.0, 0.0, 0.0),
        nalgebra::Point3::new(1.0, 1.0, 1.0),
    );

    group.bench_function("volume", |b| {
        b.iter(|| black_box(domain.volume()));
    });

    group.bench_function("contains_point", |b| {
        let point = nalgebra::Point3::new(0.5, 0.5, 0.5);
        b.iter(|| black_box(domain.contains_point(&point)));
    });

    group.finish();
}

#[cfg(feature = "simd")]
fn benchmark_simd_operations(c: &mut Criterion) {
    use cfd_math::simd::SimdOps;

    let mut group = c.benchmark_group("simd");
    let simd_ops = SimdOps::new();

    for size in [128, 1024, 8192].iter() {
        let a: Vec<f32> = (0..*size).map(|i| i as f32).collect();
        let b: Vec<f32> = (0..*size).map(|i| (i as f32) * 2.0).collect();
        let mut result = vec![0.0f32; *size];

        group.bench_with_input(BenchmarkId::new("add_f32", size), size, |bench, _| {
            bench.iter(|| simd_ops.add(&a, &b, &mut result).unwrap());
        });

        group.bench_with_input(BenchmarkId::new("dot_f32", size), size, |bench, _| {
            bench.iter(|| simd_ops.dot(&a, &b).unwrap());
        });
    }

    group.finish();
}

#[cfg(feature = "simd")]
criterion_group!(
    benches,
    benchmark_finite_difference,
    benchmark_domain_operations,
    benchmark_simd_operations
);

#[cfg(not(feature = "simd"))]
criterion_group!(
    benches,
    benchmark_finite_difference,
    benchmark_domain_operations
);
criterion_main!(benches);
