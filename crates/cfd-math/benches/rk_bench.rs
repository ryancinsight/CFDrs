use cfd_math::time_stepping::runge_kutta::LowStorageRK4;
use cfd_math::time_stepping::traits::TimeStepper;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nalgebra::DVector;

fn bench_rk_step(c: &mut Criterion) {
    let mut group = c.benchmark_group("rk_step");
    // Test with different vector sizes
    for &n in &[1000, 10000, 100000] {
        let rk4 = LowStorageRK4::<f64>::new();
        let u0 = DVector::from_element(n, 1.0);
        let dt = 0.01;
        let t = 0.0;

        // Simple function: f(u) = -u
        let f = |_t: f64, u: &DVector<f64>| Ok(-u.clone());

        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &_| {
            b.iter(|| {
                black_box(rk4.step(f, t, &u0, dt).unwrap())
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_rk_step);
criterion_main!(benches);
