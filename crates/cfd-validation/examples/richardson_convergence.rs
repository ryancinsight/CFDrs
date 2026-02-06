//! Richardson Extrapolation for Grid Convergence Study
//!
//! Demonstrates the use of Richardson extrapolation to estimate the 
//! "exact" solution from simulations on multiple grid levels and 
//! calculate the observed order of accuracy.

use cfd_validation::benchmarks::{Benchmark, LidDrivenCavity, BenchmarkConfig};
use cfd_core::error::Result;

fn main() -> Result<()> {
    println!("CFD Validation: Richardson Extrapolation Grid Convergence");
    println!("===========================================================");

    let re = 100.0_f64;
    let cavity = LidDrivenCavity::new(1.0_f64, 1.0_f64, re);

    let resolutions = vec![16, 32, 64];
    let mut results = Vec::new();

    for &res in &resolutions {
        let config = BenchmarkConfig {
            resolution: res,
            tolerance: 1e-6_f64,
            max_iterations: 1000,
            reynolds_number: re,
            ..Default::default()
        };
        
        println!("Running resolution {}x{}...", res, res);
        let result = cavity.run(&config)?;
        let sample_value = result.values[result.values.len() / 2];
        results.push(sample_value);
    }

    // Richardson extrapolation calculation
    let f1 = results[2]; // Grid 64
    let f2 = results[1]; // Grid 32
    let f3 = results[0]; // Grid 16
    let r = 2.0_f64; // Refinement ratio
    
    // Observed order of accuracy p
    let p = ((f3 - f2) / (f2 - f1)).abs().ln() / r.ln();
    let f_rich = f1 + (f1 - f2) / (r.powf(p) - 1.0_f64);
    let gci = 1.25_f64 * ((f1 - f2) / f1).abs() / (r.powf(p) - 1.0_f64);

    println!("\nConvergence Metrics:");
    println!("Grid h1 (64): {:.6}", f1);
    println!("Grid h2 (32): {:.6}", f2);
    println!("Grid h3 (16): {:.6}", f3);
    println!("Observed Order of Accuracy (p): {:.3}", p);
    println!("Extrapolated Exact Value: {:.6}", f_rich);
    println!("Grid Convergence Index (GCI): {:.4}%", gci * 100.0);

    Ok(())
}
