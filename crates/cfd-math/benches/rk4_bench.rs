use cfd_math::time_stepping::runge_kutta::RungeKutta4;
use cfd_math::time_stepping::traits::{TimeState, TimeStepper};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use leto::Array1;

fn state_from_value(len: usize, value: f64) -> TimeState<f64> {
    Array1::from_elem([len], value)
}

fn state_neg(state: &TimeState<f64>) -> TimeState<f64> {
    let len = state.shape()[0];
    Array1::from_shape_vec([len], (0..len).map(|i| -state[i]).collect())
        .expect("invariant: benchmark state length matches shape")
}

fn rk4_benchmark(c: &mut Criterion) {
    let n = 10000;
    let u0 = state_from_value(n, 1.0);
    let dt = 0.01;
    let t = 0.0;
    let rk4 = RungeKutta4::<f64>::new();

    // Simple RHS: du/dt = -u
    // Note: This closure allocates, which contributes to the cost.
    // Ideally we'd use a non-allocating RHS if possible, but the API requires allocation.
    // However, we are testing the stepper overhead, so we want to see reduction in stepper allocations.
    let f = |_t: f64, u: &TimeState<f64>| Ok(state_neg(u));

    c.bench_function("rk4_step_10k", |b| {
        b.iter(|| {
            rk4.step(black_box(f), black_box(t), black_box(&u0), black_box(dt))
                .unwrap()
        })
    });
}

criterion_group!(benches, rk4_benchmark);
criterion_main!(benches);
