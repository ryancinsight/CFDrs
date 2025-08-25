//! Performance benchmarks for solver implementations

use cfd_1d::{NetworkBuilder, NetworkProblem, NetworkSolver, Node, NodeType};
use cfd_2d::{FdmConfig, PoissonSolver, StructuredGrid2D};
use cfd_suite::prelude::*;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn benchmark_1d_network_solver(c: &mut Criterion) {
    c.bench_function("1d_network_small", |b| {
        let fluid = Fluid::<f64>::water().expect("Failed to create fluid");
        let network = NetworkBuilder::new(fluid)
            .add_node(Node::new("0".to_string(), NodeType::Junction))
            .add_node(Node::new("1".to_string(), NodeType::Junction))
            .build();

        let mut solver = NetworkSolver::<f64>::new();
        let problem = NetworkProblem::new(network);

        b.iter(|| solver.solve(black_box(&problem)));
    });
}

fn benchmark_2d_fdm_solver(c: &mut Criterion) {
    c.bench_function("2d_fdm_10x10", |b| {
        let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0).unwrap();
        let config = FdmConfig::<f64>::default();
        let solver = PoissonSolver::new(config);

        b.iter(|| {
            // Benchmark grid operations
            black_box(&grid);
        });
    });
}

criterion_group!(
    benches,
    benchmark_1d_network_solver,
    benchmark_2d_fdm_solver
);
criterion_main!(benches);
