use cfd_math::time_stepping::imex::IMEXTimeStepper;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nalgebra::{DMatrix, DVector};

fn bench_imex_step(c: &mut Criterion) {
    let mut group = c.benchmark_group("imex_step");
    // Test with different vector sizes
    for &n in &[1000, 10000] {
        let imex = IMEXTimeStepper::<f64>::ars343();
        let u0 = DVector::from_element(n, 1.0);
        let dt = 0.01;
        let t = 0.0;

        // Simple explicit function: f(u) = -u
        let f_explicit = |_t: f64, u: &DVector<f64>| Ok(-u.clone());
        // Simple implicit function: f(u) = 0
        let f_implicit = |_t: f64, _u: &DVector<f64>| Ok(DVector::zeros(n));
        let jacobian_zero = |_t: f64, _u: &DVector<f64>| Ok(DMatrix::zeros(n, n));

        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &_| {
            b.iter(|| {
                black_box(
                    imex.imex_step(f_explicit, f_implicit, jacobian_zero, t, &u0, dt)
                        .unwrap(),
                )
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_imex_step);
criterion_main!(benches);
