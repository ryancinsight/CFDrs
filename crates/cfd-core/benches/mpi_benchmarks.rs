//! Comprehensive benchmarks for MPI parallelization components
//!
//! Run with: cargo bench --features mpi --bench mpi_benchmarks

#[cfg(feature = "mpi")]
use cfd_core::compute::mpi::*;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

#[cfg(feature = "mpi")]
fn bench_domain_decomposition_creation(c: &mut Criterion) {
    let universe = MpiUniverse::new().unwrap();
    let world = universe.world();

    c.bench_function("domain_decomposition_creation", |b| {
        b.iter(|| {
            let global_extents = GlobalExtents::new_2d(1000, 1000, (0.0, 1.0, 0.0, 1.0));
            let strategy = DecompositionStrategy::Cartesian2D;
            black_box(DomainDecomposition::new(global_extents, &world, strategy).unwrap())
        })
    });
}

#[cfg(feature = "mpi")]
fn bench_load_balancer_assessment(c: &mut Criterion) {
    let universe = MpiUniverse::new().unwrap();
    let world = universe.world();
    let global_extents = GlobalExtents::new_2d(1000, 1000, (0.0, 1.0, 0.0, 1.0));
    let strategy = DecompositionStrategy::Cartesian2D;
    let decomp = DomainDecomposition::new(global_extents, &world, strategy).unwrap();
    let load_balancer = LoadBalancer::new(&world, decomp, 1.2, 100).unwrap();

    c.bench_function("load_balance_assessment", |b| {
        b.iter(|| {
            let local_workload = 10000;
            black_box(load_balancer.assess_load_balance(local_workload).unwrap())
        })
    });
}

#[cfg(feature = "mpi")]
fn bench_ghost_cell_exchange(c: &mut Criterion) {
    let universe = MpiUniverse::new().unwrap();
    let world = universe.world();
    let global_extents = GlobalExtents::new_2d(100, 100, (0.0, 1.0, 0.0, 1.0));
    let strategy = DecompositionStrategy::Cartesian2D;
    let decomp = DomainDecomposition::new(global_extents, &world, strategy).unwrap();

    // Create test data
    let mut velocity_u = vec![vec![nalgebra::Vector2::new(1.0, 2.0); 102]; 102];
    let mut velocity_v = vec![vec![nalgebra::Vector2::new(3.0, 4.0); 102]; 102];
    let mut pressure = vec![vec![5.0; 102]; 102];

    c.bench_function("ghost_cell_exchange", |b| {
        b.iter(|| {
            black_box(
                update_ghost_cells(
                    &world,
                    &decomp,
                    &mut velocity_u,
                    &mut velocity_v,
                    &mut pressure,
                )
                .unwrap(),
            )
        })
    });
}

#[cfg(feature = "mpi")]
fn bench_performance_metrics_calculation(c: &mut Criterion) {
    use std::time::Duration;

    c.bench_function("performance_metrics_calculation", |b| {
        b.iter(|| {
            let mut metrics = PerformanceMetrics::new();
            metrics.num_processes = 4;
            metrics.total_time = Duration::from_secs(100);
            metrics.communication_time = Duration::from_secs(10);
            metrics.computation_time = Duration::from_secs(85);
            metrics.load_imbalance_ratio = 1.1;

            black_box(metrics.calculate_derived_metrics());
            metrics
        })
    });
}

#[cfg(feature = "mpi")]
fn bench_scaling_assessment(c: &mut Criterion) {
    use performance_validation::ScalingTestResult;

    let results = ScalingTestResult {
        core_counts: vec![1, 2, 4, 8, 16, 32, 64],
        metrics: vec![PerformanceMetrics::new(); 7],
        scaling_efficiency: vec![1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.82],
        communication_trend: vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0],
        load_imbalance_trend: vec![1.0, 1.01, 1.02, 1.03, 1.05, 1.08, 1.15],
        test_type: performance_validation::ScalingTestType::Strong,
        assessment: performance_validation::ScalingAssessment {
            grade: performance_validation::ScalingGrade::F,
            efficiency_at_64_cores: 0.0,
            comm_overhead_at_64_cores: 0.0,
            load_imbalance_at_64_cores: 0.0,
            notes: Vec::new(),
            recommendations: Vec::new(),
        },
    };

    c.bench_function("scaling_assessment", |b| {
        b.iter(|| {
            let mut results_copy = results.clone();
            black_box(results_copy.assess_scaling());
            results_copy
        })
    });
}

#[cfg(feature = "mpi")]
criterion_group! {
    name = mpi_benches;
    config = Criterion::default().sample_size(10);
    targets = bench_domain_decomposition_creation,
             bench_load_balancer_assessment,
             bench_ghost_cell_exchange,
             bench_performance_metrics_calculation,
             bench_scaling_assessment
}

#[cfg(feature = "mpi")]
criterion_main!(mpi_benches);

#[cfg(not(feature = "mpi"))]
fn main() {
    eprintln!("MPI benchmarks require the 'mpi' feature to be enabled.");
    eprintln!("Run with: cargo bench --features mpi --bench mpi_benchmarks");
    std::process::exit(1);
}
