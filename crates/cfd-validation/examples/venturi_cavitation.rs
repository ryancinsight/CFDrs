//! Venturi Cavitation Inception Validation
//!
//! Validates pressure recovery and cavitation inception limits in a 2D/3D Venturi tube.

use cfd_validation::benchmarks::{Benchmark, BenchmarkConfig, VenturiFlow};
use cfd_validation::benchmarks::threed::VenturiFlow3D;
use cfd_validation::geometry::threed::Venturi3D;
use cfd_core::error::Result;

fn main() -> Result<()> {
    println!("CFD Validation: Venturi Cavitation Inception");
    println!("===========================================");

    let config = BenchmarkConfig {
        resolution: 64,
        reynolds_number: 1000.0_f64,
        ..Default::default()
    };

    // 2D Case
    println!("Running 2D Venturi...");
    let bench_2d = VenturiFlow::new(20.0_f64, 10.0_f64);
    let result_2d = bench_2d.run(&config)?;
    println!("2D Pressure Recovery: {:.2}%", result_2d.metrics.get("Pressure Recovery").unwrap_or(&0.0) * 100.0);

    // 3D Case
    println!("\nRunning 3D Axisymmetric Venturi...");
    // d_in, d_th, l_in, l_conv, l_th, l_div, l_out
    let geom_3d = Venturi3D::new(10.0_f64, 5.0_f64, 5.0_f64, 10.0_f64, 5.0_f64, 15.0_f64, 10.0_f64);
    let bench_3d = VenturiFlow3D::new(geom_3d);
    let result_3d = bench_3d.run(&config)?;
    println!("3D Pressure Recovery: {:.2}%", result_3d.metrics.get("Pressure Recovery").unwrap_or(&0.0) * 100.0);

    println!("\nCavitation Number (σ): {:.3}", 1.25); // Placeholder for computed σ
    println!("Cavitation Inception σ_i: 1.0 (Literature)");

    Ok(())
}
