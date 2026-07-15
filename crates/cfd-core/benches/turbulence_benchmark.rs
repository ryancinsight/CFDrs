use cfd_core::compute::gpu::{GpuTurbulenceCompute, TurbulenceGrid};
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_compute_des(c: &mut Criterion) {
    // Initialize GPU compute
    // We do this inside the benchmark function but outside the loop
    let compute = match GpuTurbulenceCompute::new() {
        Ok(c) => c,
        Err(e) => {
            println!("Skipping benchmark: GPU initialization failed: {}", e);
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
        })
    });
}

criterion_group!(benches, bench_compute_des);
criterion_main!(benches);
