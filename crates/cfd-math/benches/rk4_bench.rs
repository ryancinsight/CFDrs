use criterion::{black_box, criterion_group, criterion_main, Criterion};
use cfd_math::time_stepping::runge_kutta::RungeKutta4;
use cfd_math::time_stepping::traits::TimeStepper;
use nalgebra::DVector;

fn rk4_benchmark(c: &mut Criterion) {
    let n = 10000;
    let u0 = DVector::from_element(n, 1.0);
    let dt = 0.01;
    let t = 0.0;
    let rk4 = RungeKutta4::<f64>::new();

    // Simple RHS: du/dt = -u
    // Note: This closure allocates, which contributes to the cost.
    // Ideally we'd use a non-allocating RHS if possible, but the API requires allocation.
    // However, we are testing the stepper overhead, so we want to see reduction in stepper allocations.
    let f = |_t: f64, u: &DVector<f64>| Ok(-u.clone());

    c.bench_function("rk4_step_10k", |b| {
        b.iter(|| {
            rk4.step(black_box(f), black_box(t), black_box(&u0), black_box(dt)).unwrap()
        })
    });
}

criterion_group!(benches, rk4_benchmark);
criterion_main!(benches);
