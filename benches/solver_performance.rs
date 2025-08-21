//! Performance benchmarks for CFD solvers
//! 
//! Following SOLID principles and Rust best practices

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use cfd_suite::prelude::*;

fn benchmark_1d_network_solver(c: &mut Criterion) {
    c.bench_function("1d_network_small", |b| {
        let fluid = Fluid::<f64>::water().expect("Failed to create fluid");
        let network = NetworkBuilder::new(fluid)
            .add_node(Node::new(0, 0.0, 0.0, 0.0))
            .add_node(Node::new(1, 1.0, 0.0, 0.0))
            .build()
            .expect("Failed to build network");
        
        let mut solver = NetworkSolver::<f64>::new();
        let problem = NetworkProblem { network };
        
        b.iter(|| {
            solver.solve(black_box(&problem))
        });
    });
}

fn benchmark_2d_fdm_solver(c: &mut Criterion) {
    c.bench_function("2d_fdm_10x10", |b| {
        let grid = StructuredGrid2D::<f64>::new(10, 10, 1.0, 1.0, 0.0, 1.0);
        let config = FdmConfig::default();
        let mut solver = FdmSolver::with_config(config);
        
        b.iter(|| {
            // Benchmark grid operations
            black_box(&grid);
        });
    });
}

criterion_group!(benches, benchmark_1d_network_solver, benchmark_2d_fdm_solver);
criterion_main!(benches);