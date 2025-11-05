//! Performance benchmark example for MPI scaling validation
//!
//! This example demonstrates how to run comprehensive scaling benchmarks
//! for the CFD suite's MPI parallelization.
//!
//! Run with:
//! ```bash
//! cargo run --example performance_benchmark --features mpi -- --cores 1,2,4,8
//! ```

#[cfg(feature = "mpi")]
use cfd_core::compute::mpi::*;
use std::time::Duration;

#[cfg(feature = "mpi")]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Initialize MPI
    let universe = MpiUniverse::new()?;
    let world = universe.world();

    // Parse command line arguments
    let args: Vec<String> = std::env::args().collect();
    let core_counts: Vec<usize> = if args.len() > 2 && args[1] == "--cores" {
        args[2].split(',')
            .map(|s| s.trim().parse().unwrap_or(1))
            .collect()
    } else {
        vec![1, 2, 4, 8, 16] // Default core counts
    };

    println!("MPI Performance Benchmark");
    println!("=========================");
    println!("World size: {}", world.size());
    println!("World rank: {}", world.rank());
    println!("Testing core counts: {:?}", core_counts);
    println!();

    // Create performance validator
    let validator = PerformanceValidator::<f64>::new(&world);

    // Run strong scaling benchmark
    if world.is_root() {
        println!("Running Strong Scaling Benchmark...");
        println!("===================================");
    }

    let strong_results = validator.run_strong_scaling_test(
        |comm, cores, tolerance| run_cfd_simulation(comm, cores, tolerance, ScalingTestType::Strong),
        &core_counts,
        1e-6,
    )?;

    if world.is_root() {
        println!("Strong Scaling Results:");
        print_scaling_results(&strong_results);
        println!();
    }

    // Run weak scaling benchmark
    if world.is_root() {
        println!("Running Weak Scaling Benchmark...");
        println!("=================================");
    }

    let weak_results = validator.run_weak_scaling_test(
        |comm, cores, tolerance| run_cfd_simulation(comm, cores, tolerance, ScalingTestType::Weak),
        &core_counts,
        1e-6,
    )?;

    if world.is_root() {
        println!("Weak Scaling Results:");
        print_scaling_results(&weak_results);
        println!();
    }

    // Assess production readiness
    if world.is_root() {
        println!("Production Readiness Assessment");
        println!("==============================");
        let readiness_report = validator.assess_production_readiness()?;

        println!("Overall Score: {}%", readiness_report.overall_score);
        println!("Component Scores:");
        for (component, score) in &readiness_report.component_scores {
            println!("  {}: {}%", component, score);
        }

        if !readiness_report.critical_issues.is_empty() {
            println!("Critical Issues:");
            for issue in &readiness_report.critical_issues {
                println!("  - {}", issue);
            }
        }

        println!();
        println!("Recommended Configuration:");
        println!("  Cores per node: {}", readiness_report.recommended_config.cores_per_node);
        println!("  MPI implementation: {}", readiness_report.recommended_config.mpi_implementation);
        println!("  Memory per core: {} MB", readiness_report.recommended_config.memory_per_core_mb);
        println!("  Network: {}", readiness_report.recommended_config.network_requirements);
        println!("  Recommended problem sizes: {:?}", readiness_report.recommended_config.recommended_problem_sizes);
    }

    Ok(())
}

#[cfg(feature = "mpi")]
fn run_cfd_simulation(
    communicator: &MpiCommunicator,
    cores: usize,
    tolerance: f64,
    test_type: ScalingTestType,
) -> MpiResult<PerformanceMetrics> {
    use super::performance_validation::{PerformanceTimer, SimulationData};

    let mut metrics = PerformanceMetrics::new();
    let mut comp_timer = PerformanceTimer::new();
    let mut comm_timer = PerformanceTimer::new();
    let mut io_timer = PerformanceTimer::new();

    metrics.num_processes = cores as usize;

    // Determine problem size based on test type
    let base_size = 1000; // Base problem size
    metrics.problem_size = match test_type {
        ScalingTestType::Strong => base_size, // Fixed size for strong scaling
        ScalingTestType::Weak => base_size * cores, // Scale with cores for weak scaling
    };

    // Simulate domain decomposition
    let global_extents = GlobalExtents::new_2d(
        metrics.problem_size,
        metrics.problem_size,
        (0.0, 1.0, 0.0, 1.0)
    );
    let strategy = DecompositionStrategy::Cartesian2D;
    let decomp = DomainDecomposition::new(global_extents, communicator, strategy)?;

    // Simulate CFD computation
    comp_timer.start();

    // Simulate several time steps
    let num_iterations = 10;
    for iteration in 0..num_iterations {
        // Simulate computation on local subdomain
        simulate_computation_step(&decomp)?;

        // Simulate communication (ghost cell exchange)
        comm_timer.start();
        simulate_communication(&decomp)?;
        comm_timer.stop();

        // Simulate I/O (occasional checkpoints)
        if iteration % 5 == 0 {
            io_timer.start();
            simulate_io()?;
            io_timer.stop();
        }
    }

    comp_timer.stop();

    // Fill in metrics
    metrics.total_time = comp_timer.total_time() + comm_timer.total_time() + io_timer.total_time();
    metrics.computation_time = comp_timer.total_time();
    metrics.communication_time = comm_timer.total_time();
    metrics.io_time = io_timer.total_time();

    // Simulate load imbalance (would be measured in real simulation)
    metrics.load_imbalance_ratio = 1.0 + (communicator.rank() as f64 * 0.1);

    // Calculate derived metrics
    metrics.calculate_derived_metrics();

    // Simulate memory usage
    metrics.memory_usage_mb = (metrics.problem_size as f64 / cores as f64 * 8.0) / (1024.0 * 1024.0);

    Ok(metrics)
}

#[cfg(feature = "mpi")]
fn simulate_computation_step(decomp: &DomainDecomposition) -> MpiResult<()> {
    // Simulate computational work proportional to local subdomain size
    let local_cells = decomp.local_subdomain().nx_local * decomp.local_subdomain().ny_local;
    let work_iterations = local_cells / 100; // Simulate some computational work

    for _ in 0..work_iterations {
        // Simulate floating point operations
        let _dummy = (0..100).map(|x| (x as f64).sin()).sum::<f64>();
    }

    Ok(())
}

#[cfg(feature = "mpi")]
fn simulate_communication(decomp: &DomainDecomposition) -> MpiResult<()> {
    // Simulate MPI communication overhead
    // In a real simulation, this would be actual ghost cell exchanges

    // Small delay to simulate communication latency
    if decomp.communicator().size() > 1 {
        std::thread::sleep(Duration::from_micros(100));
    }

    Ok(())
}

#[cfg(feature = "mpi")]
fn simulate_io() -> MpiResult<()> {
    // Simulate I/O operations (checkpointing, etc.)
    std::thread::sleep(Duration::from_millis(5));
    Ok(())
}

#[cfg(feature = "mpi")]
fn print_scaling_results(results: &ScalingTestResult) {
    println!("Core Count | Efficiency | Comm Overhead | Load Imbalance");
    println!("-----------|------------|---------------|----------------");

    for i in 0..results.core_counts.len() {
        let cores = results.core_counts[i];
        let efficiency = results.scaling_efficiency.get(i).copied().unwrap_or(0.0);
        let comm_overhead = results.communication_trend.get(i).copied().unwrap_or(0.0);
        let load_imbalance = results.load_imbalance_trend.get(i).copied().unwrap_or(1.0);

        println!("{:>10} | {:>9.3} | {:>12.1}% | {:>14.3}",
                cores, efficiency, comm_overhead, load_imbalance);
    }

    println!();
    println!("Assessment Grade: {:?}", results.assessment.grade);
    println!("Efficiency at 64 cores: {:.3}", results.assessment.efficiency_at_64_cores);
    println!("Comm overhead at 64 cores: {:.1}%", results.assessment.comm_overhead_at_64_cores);
    println!("Load imbalance at 64 cores: {:.3}", results.assessment.load_imbalance_at_64_cores);

    if !results.assessment.notes.is_empty() {
        println!("Notes:");
        for note in &results.assessment.notes {
            println!("  - {}", note);
        }
    }

    if !results.assessment.recommendations.is_empty() {
        println!("Recommendations:");
        for rec in &results.assessment.recommendations {
            println!("  - {}", rec);
        }
    }
}

#[cfg(not(feature = "mpi"))]
fn main() {
    eprintln!("This example requires the 'mpi' feature to be enabled.");
    eprintln!("Run with: cargo run --example performance_benchmark --features mpi");
    std::process::exit(1);
}
