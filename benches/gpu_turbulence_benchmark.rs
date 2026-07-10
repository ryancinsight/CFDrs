//! Hephaestus turbulence-kernel benchmarks.

use cfd_core::compute::gpu::{GpuTurbulenceCompute, TurbulenceGrid};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

fn bench_smagorinsky(c: &mut Criterion) {
    let compute = match GpuTurbulenceCompute::new() {
        Ok(compute) => compute,
        Err(error) => {
            eprintln!("GPU turbulence benchmark unavailable: {error}");
            return;
        }
    };
    let mut group = c.benchmark_group("Hephaestus Smagorinsky");
    for size in [64, 128, 256] {
        let grid = TurbulenceGrid::new([size, size], [0.01; 2]).unwrap();
        let velocity_x: Vec<f32> = (0..grid.element_count())
            .map(|index| (index % size) as f32 * 0.01)
            .collect();
        let velocity_y: Vec<f32> = (0..grid.element_count())
            .map(|index| (index / size) as f32 * 0.01)
            .collect();
        let mut output = vec![0.0; grid.element_count()];
        group.bench_with_input(BenchmarkId::from_parameter(size), &size, |b, _| {
            b.iter(|| {
                compute
                    .compute_smagorinsky_sgs(
                        black_box(&velocity_x),
                        black_box(&velocity_y),
                        grid,
                        black_box(0.1),
                        black_box(&mut output),
                    )
                    .unwrap();
            });
        });
    }
    group.finish();
}

fn bench_grid_quantities(c: &mut Criterion) {
    let compute = match GpuTurbulenceCompute::new() {
        Ok(compute) => compute,
        Err(error) => {
            eprintln!("GPU turbulence benchmark unavailable: {error}");
            return;
        }
    };
    let grid = TurbulenceGrid::new([256, 256], [0.01; 2]).unwrap();
    let mut output = vec![0.0; grid.element_count()];
    c.bench_function("Hephaestus DES grid scale", |b| {
        b.iter(|| {
            compute
                .compute_des_length_scale(grid, black_box(0.65), black_box(&mut output))
                .unwrap();
        });
    });
    c.bench_function("Hephaestus wall distance", |b| {
        b.iter(|| {
            compute
                .compute_wall_distance(grid, black_box(&mut output))
                .unwrap();
        });
    });
}

criterion_group!(benches, bench_smagorinsky, bench_grid_quantities);
criterion_main!(benches);
