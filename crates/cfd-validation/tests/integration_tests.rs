//! Integration tests for the CFD Validation Suite
//!
//! Verifies that the complete validation chain (Geometry -> Solver -> Benchmark -> Result)
//! works correctly for representative 2D and 3D cases.

use cfd_core::error::Result;
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_validation::benchmarks::{Benchmark, BenchmarkConfig, BifurcationFlow3D, LidDrivenCavity};
use cfd_validation::geometry::threed::Bifurcation3D;

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
    // 3D bifurcation test using SI metre dimensions that are consistent with
    // BifurcationGeometry3D (100 µm parent, 80 µm daughters, 1 mm lengths).
    // These match the BifurcationFlow3D benchmark's internal expectations and
    // the mesh builder's SI-metre coordinate system.
    let d_parent = 100e-6_f64; // 100 µm diameter
    let d_daughter = 80e-6_f64; // 80 µm diameter (Murray-optimal: 0.8^3 + 0.8^3 ≈ 1.02^3...)
    let l_parent = 1e-3_f64; // 1 mm
    let l_daughter = 1e-3_f64; // 1 mm
    let angle = 0.5_f64; // 0.5 rad branching half-angle

    let geom = Bifurcation3D::symmetric(d_parent, d_daughter, l_parent, l_daughter, angle);
    let fluid = CassonBlood::normal_blood();
    let bench = BifurcationFlow3D::new(geom, fluid);

    let config = BenchmarkConfig {
        resolution: 4, // Coarse mesh for fast CI
        tolerance: 1e-2_f64,
        max_iterations: 10,
        reynolds_number: 1.0_f64, // Re ≈ 1 (Stokes regime) for stable FEM solve
        time_step: Some(0.01_f64),
        parallel: false,
    };

    let result = bench.run(&config)?;

    assert!(result.name.contains("Bifurcation"));
    assert!(result.metrics.contains_key("Murray Deviation"));

    Ok(())
}
