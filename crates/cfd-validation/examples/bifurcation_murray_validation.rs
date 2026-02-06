//! Bifurcation validation against Murray's Law
//!
//! Validates flow distribution and pressure continuity at a 2D/3D bifurcation apex.

use cfd_validation::benchmarks::{Benchmark, BenchmarkConfig, BifurcationFlow};
use cfd_validation::benchmarks::threed::BifurcationFlow3D;
use cfd_validation::geometry::threed::Bifurcation3D;
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::error::Result;

fn main() -> Result<()> {
    println!("CFD Validation: Bifurcation Murray's Law Check");
    println!("=============================================");

    let config = BenchmarkConfig {
        resolution: 32,
        ..Default::default()
    };

    // 2D Case
    println!("Running 2D Symmetric Bifurcation...");
    let bench_2d = BifurcationFlow::new(10.0_f64, 50.0_f64, 0.5_f64);
    let result_2d = bench_2d.run(&config)?;
    println!("2D Murray Deviation: {:.4}%", result_2d.metrics.get("Murray Deviation").unwrap_or(&0.0) * 100.0);

    // 3D Case
    println!("\nRunning 3D Symmetric Bifurcation...");
    let geom_3d = Bifurcation3D::symmetric(5.0_f64, 4.0_f64, 20.0_f64, 15.0_f64, 0.5_f64);
    let bench_3d = BifurcationFlow3D::new(geom_3d, CassonBlood::normal_blood());
    let result_3d = bench_3d.run(&config)?;
    println!("3D Murray Deviation: {:.4}%", result_3d.metrics.get("Murray Deviation").unwrap_or(&0.0) * 100.0);

    println!("\nValidation Result: {}", if bench_2d.validate(&result_2d)? && bench_3d.validate(&result_3d)? { "PASS" } else { "FAIL" });

    Ok(())
}
