//! Performance benchmarks for solver implementations

use cfd_1d::network::{Network, NetworkBuilder};
use cfd_1d::NetworkProblem;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::solvers::fdm::{FdmConfig, PoissonSolver};
use cfd_core::error::Result;
use cfd_core::physics::fluid::Fluid;
use cfd_suite::prelude::*;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

struct NetworkBenchmarkContext {
    solver: cfd_1d::solver::NetworkSolver<f64>,
    problem: NetworkProblem<f64>,
}

fn build_network_benchmark_context() -> Result<NetworkBenchmarkContext> {
    let fluid = Fluid::<f64>::water_20c()?;
    let mut builder = NetworkBuilder::new();
    let n0 = builder.add_inlet("0".to_string());
    let n1 = builder.add_outlet("1".to_string());
    builder.connect_with_pipe(n0, n1, "pipe".to_string());
    let graph = builder.build()?;
    let network = Network::new(graph, fluid);

    let solver = cfd_1d::solver::NetworkSolver::<f64>::new();
    let problem = NetworkProblem::new(network);

    Ok(NetworkBenchmarkContext { solver, problem })
}

fn build_fdm_benchmark_grid() -> Result<StructuredGrid2D<f64>> {
    StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0)
}

fn benchmark_1d_network_solver(c: &mut Criterion) {
    c.bench_function("1d_network_small", |b| {
        let context = match build_network_benchmark_context() {
            Ok(context) => context,
            Err(error) => panic!("Failed to build 1D network benchmark context: {error}"),
        };
        let solver = context.solver;
        let problem = context.problem;

        b.iter(|| solver.solve(black_box(&problem)));
    });
}

fn benchmark_2d_fdm_solver(c: &mut Criterion) {
    c.bench_function("2d_fdm_10x10", |b| {
        let grid = match build_fdm_benchmark_grid() {
            Ok(grid) => grid,
            Err(error) => panic!("Failed to build 2D FDM benchmark grid: {error}"),
        };
        let config = FdmConfig::<f64>::default();
        let _solver = PoissonSolver::new(config);

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
