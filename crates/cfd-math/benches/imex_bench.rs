use cfd_math::time_stepping::{IMEXTimeStepper, TimeMatrix, TimeState};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use leto::Array1;

fn state_from_value(len: usize, value: f64) -> TimeState<f64> {
    Array1::from_elem([len], value)
}

fn state_zeros(len: usize) -> TimeState<f64> {
    state_from_value(len, 0.0)
}

fn state_neg(state: &TimeState<f64>) -> TimeState<f64> {
    let len = state.shape()[0];
    Array1::from_shape_vec([len], (0..len).map(|i| -state[i]).collect())
        .expect("invariant: benchmark state length matches shape")
}

fn matrix_zeros(rows: usize, cols: usize) -> TimeMatrix<f64> {
    TimeMatrix::from_elem([rows, cols], 0.0)
}

fn bench_imex_step(c: &mut Criterion) {
    let mut group = c.benchmark_group("imex_step");
    // Test with different vector sizes
    for &n in &[1000, 10000] {
        let imex = IMEXTimeStepper::<f64>::ars343();
        let u0 = state_from_value(n, 1.0);
        let dt = 0.01;
        let t = 0.0;

        // Simple explicit function: f(u) = -u
        let f_explicit = |_t: f64, u: &TimeState<f64>| Ok(state_neg(u));
        // Simple implicit function: f(u) = 0
        let f_implicit = |_t: f64, u: &TimeState<f64>| Ok(state_zeros(u.shape()[0]));
        let jacobian_zero =
            |_t: f64, u: &TimeState<f64>| Ok(matrix_zeros(u.shape()[0], u.shape()[0]));

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
