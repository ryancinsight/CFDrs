//! Example demonstrating benchmark validation for CFD problems

use cfd_validation::benchmarks::{Benchmark, LidDrivenCavity, FlowOverCylinder};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== CFD Benchmark Validation Example ===\n");
    
    // Test 1: Lid-Driven Cavity
    println!("1. Running Lid-Driven Cavity Benchmark (Re=100)...");
    let mut cavity = LidDrivenCavity::<f64>::new(100.0, (32, 32), 1.0);
    
    println!("   Setting up benchmark...");
    cavity.setup(Default::default())?;
    
    println!("   Running SIMPLE solver (this may take a moment)...");
    let start = std::time::Instant::now();
    let solution = cavity.run()?;
    let elapsed = start.elapsed();
    println!("   Solver completed in {:.2} seconds", elapsed.as_secs_f64());
    
    // Check that we have a non-zero solution
    let max_velocity = solution.iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max);
    println!("   Maximum velocity in solution: {:.4}", max_velocity);
    
    if max_velocity > 0.0 {
        println!("   ✓ Solution contains non-zero velocities");
    } else {
        println!("   ✗ Warning: Solution appears to be zero");
    }
    
    println!("   Validating against Ghia et al. reference data...");
    let result = cavity.validate(&solution)?;
    if result.passed {
        println!("   ✓ Validation passed!");
    } else {
        println!("   ✗ Validation failed (may need more iterations)");
    }
    
    // Print some statistics
    for (key, stats) in &result.error_statistics {
        println!("   {}: L2 error = {:.4}", key, stats.l2_norm);
    }
    
    println!();
    
    // Test 2: Flow Over Cylinder
    println!("2. Running Flow Over Cylinder Benchmark (Re=40)...");
    let mut cylinder = FlowOverCylinder::<f64>::new(40.0, 1.0);
    
    println!("   Setting up benchmark...");
    cylinder.setup(Default::default())?;
    
    println!("   Running SIMPLE solver (this may take a moment)...");
    let start = std::time::Instant::now();
    let solution = cylinder.run()?;
    let elapsed = start.elapsed();
    println!("   Solver completed in {:.2} seconds", elapsed.as_secs_f64());
    
    // Check solution
    let max_velocity = solution.iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max);
    println!("   Maximum velocity in solution: {:.4}", max_velocity);
    
    if max_velocity > 0.0 {
        println!("   ✓ Solution contains non-zero velocities");
    } else {
        println!("   ✗ Warning: Solution appears to be zero");
    }
    
    println!("   Validating drag coefficient...");
    let result = cylinder.validate(&solution)?;
    if result.passed {
        println!("   ✓ Validation passed!");
    } else {
        println!("   ✗ Validation failed (may need more iterations)");
    }
    
    // Print drag coefficient
    for (key, value) in &result.metadata {
        if key.contains("cd") {
            println!("   {}: {}", key, value);
        }
    }
    
    println!("\n=== Benchmark Validation Complete ===");
    println!("Note: Full convergence may require more iterations than the default settings.");
    println!("The benchmarks are now using proper SIMPLE algorithm for Navier-Stokes equations.");
    
    Ok(())
}