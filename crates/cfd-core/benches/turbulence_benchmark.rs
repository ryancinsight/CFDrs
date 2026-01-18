use cfd_core::compute::gpu::turbulence_compute::GpuTurbulenceCompute;
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_compute_des(c: &mut Criterion) {
    // Initialize GPU compute
    // We do this inside the benchmark function but outside the loop
    let mut compute = match GpuTurbulenceCompute::new() {
        Ok(c) => c,
        Err(e) => {
            println!("Skipping benchmark: GPU initialization failed: {}", e);
            return;
        }
    };

    let nx = 128;
    let ny = 128;
    let size = nx * ny;
    let velocity_u = vec![1.0; size];
    let velocity_v = vec![0.5; size];

    c.bench_function("compute_des_length_scale", |b| {
        b.iter(|| {
            // This calls the function that creates new buffers every time
            compute
                .compute_des_length_scale(&velocity_u, &velocity_v, nx, ny, 0.1, 0.1, 0.65)
                .expect("Computation failed");
        })
    });
}

criterion_group!(benches, bench_compute_des);
criterion_main!(benches);
