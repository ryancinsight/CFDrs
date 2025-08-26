use cfd_2d::{
    grid::StructuredGrid2D,
    solvers::lbm::{LbmConfig, LbmSolver},
};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nalgebra::Vector2;
use std::collections::HashMap;

fn benchmark_lbm_solver(c: &mut Criterion) {
    let mut group = c.benchmark_group("lbm_solver");
    for size in [32, 64, 128].iter() {
        let grid = StructuredGrid2D::new(*size, *size, 0.0, 1.0, 0.0, 1.0).unwrap();
        let config = LbmConfig::default();
        let mut solver = LbmSolver::new(config, &grid);
        // Initialize with some density and velocity
        solver
            .initialize(
                |_x: f64, _y: f64| 1.0,
                |_x: f64, _y: f64| Vector2::new(0.0, 0.0),
            )
            .unwrap();
        group.bench_with_input(BenchmarkId::new("step", size), size, |b, _| {
            b.iter(|| {
                let boundaries = HashMap::new();
                black_box(solver.step(&boundaries).unwrap())
            })
        });
    }
    group.finish();
}
fn benchmark_grid_creation(c: &mut Criterion) {
    let mut group = c.benchmark_group("grid_creation");
    for size in [50, 100, 200].iter() {
        group.bench_with_input(
            BenchmarkId::new("structured_grid", size),
            size,
            |b, &size| {
                b.iter(|| black_box(StructuredGrid2D::new(size, size, 0.0, 1.0, 0.0, 1.0).unwrap()))
            },
        );
criterion_group!(benches, benchmark_lbm_solver, benchmark_grid_creation);
criterion_main!(benches);


}
}
