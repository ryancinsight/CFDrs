//! Criterion benchmark for the provider-backed turbulence length-scale kernel.

#![expect(
    missing_docs,
    reason = "criterion_group generates the public benchmark-group function"
)]

use cfd_core::compute::gpu::{GpuTurbulenceCompute, TurbulenceGrid};
use criterion::{criterion_group, criterion_main, Criterion};

/// Measures the DES length-scale kernel on a fixed two-dimensional grid.
fn bench_compute_des(c: &mut Criterion) {
    // Initialize GPU compute
    // We do this inside the benchmark function but outside the loop
    let compute = match GpuTurbulenceCompute::new() {
        Ok(c) => c,
        Err(e) => {
            tracing::warn!(error = %e, "Skipping turbulence benchmark: GPU initialization failed");
            return;
        }
    };

    let nx = 128;
    let ny = 128;
    let size = nx * ny;
    let grid = TurbulenceGrid::new([nx, ny], [0.1, 0.1]).expect("valid benchmark grid");
    let mut output = vec![0.0; size];

    c.bench_function("compute_des_length_scale", |b| {
        b.iter(|| {
            compute
                .compute_des_length_scale(grid, 0.65, &mut output)
                .expect("Computation failed");
        });
    });
}

criterion_group!(benches, bench_compute_des);
criterion_main!(benches);
