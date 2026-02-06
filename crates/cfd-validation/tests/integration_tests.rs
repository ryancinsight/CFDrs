//! Integration tests for the CFD Validation Suite
//!
//! Verifies that the complete validation chain (Geometry -> Solver -> Benchmark -> Result)
//! works correctly for representative 2D and 3D cases.

use cfd_validation::benchmarks::{Benchmark, BenchmarkConfig, LidDrivenCavity, BifurcationFlow3D};
use cfd_validation::geometry::threed::Bifurcation3D;
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::error::Result;
use nalgebra::RealField;

#[test]
fn test_lid_driven_cavity_integration() -> Result<()> {
    let re = 100.0_f64;
    let cavity = LidDrivenCavity::new(1.0_f64, 1.0_f64, re);
    
    let config = BenchmarkConfig {
        resolution: 32,
        tolerance: 1e-4_f64,
        max_iterations: 100,
        reynolds_number: re,
        time_step: None,
        parallel: false,
    };
    
    let result = cavity.run(&config)?;
    
    assert_eq!(result.name, "Lid-Driven Cavity");
    assert!(!result.convergence.is_empty());
    assert!(cavity.validate(&result)?);
    
    Ok(())
}

#[test]
fn test_3d_bifurcation_integration() -> Result<()> {
    // Small scale 3D bifurcation test
    let parent_radius = 5.0_f64;
    let daughter_radius = 4.0_f64;
    let l_parent = 20.0_f64;
    let l_daughter = 15.0_f64;
    let angle = 0.5_f64; // rad
    
    let geom = Bifurcation3D::symmetric(parent_radius, daughter_radius, l_parent, l_daughter, angle);
    let fluid = CassonBlood::normal_blood();
    let bench = BifurcationFlow3D::new(geom, fluid);
    
    let config = BenchmarkConfig {
        resolution: 16, // Low res for fast test
        tolerance: 1e-2_f64,
        max_iterations: 10,
        reynolds_number: 50.0_f64,
        time_step: Some(0.01_f64),
        parallel: false,
    };
    
    let result = bench.run(&config)?;
    
    assert!(result.name.contains("Bifurcation"));
    assert!(result.metrics.contains_key("Murray Deviation"));
    
    Ok(())
}
