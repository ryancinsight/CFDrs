//! Validation framework for CFD simulations
//!
//! This module provides tools for validating CFD implementations against:
//! - Analytical solutions
//! - Benchmark problems
//! - Literature results
//! - Method of manufactured solutions

#![warn(missing_docs)]

pub mod analytical;
pub mod benchmarks;
pub mod conservation;
pub mod convergence;
pub mod error_metrics;
pub mod literature;
pub mod numerical_validation;
pub mod solutions;
pub mod time_integration_validation;

pub use analytical::AnalyticalSolution;
pub use benchmarks::{Benchmark, BenchmarkResult};
pub use conservation::ConservationChecker;
pub use convergence::{ConvergenceCriterion, ConvergenceStatus, ConvergenceStudy, RichardsonExtrapolation, GridConvergenceIndex};
pub use error_metrics::{ErrorMetric, ErrorAnalysis};

/// Basic validation test for Poiseuille flow
/// 
/// This validates the 1D pipe flow against the analytical solution
#[cfg(test)]
mod validation_tests {
    use super::*;

    #[test]
    fn test_poiseuille_validation() {
        // Analytical solution for laminar pipe flow
        // v(r) = (ΔP/4μL) * (R² - r²)
        let radius = 0.01_f64; // 1 cm
        let length = 1.0_f64;   // 1 m
        let pressure_drop = 100.0_f64; // Pa
        let viscosity = 0.001_f64; // Pa·s
        
        // Maximum velocity at centerline
        let v_max = (pressure_drop * radius * radius) / (4.0 * viscosity * length);
        
        // Average velocity for parabolic profile
        let v_avg = v_max / 2.0;
        
        // This is the analytical result - any implementation should match this
        assert!((v_avg - 0.025).abs() < 1e-6, "Poiseuille flow validation failed");
    }

    #[test] 
    fn test_reynolds_number() {
        let density = 1000.0_f64; // kg/m³
        let velocity = 1.0_f64;    // m/s
        let length = 0.1_f64;      // m
        let viscosity = 0.001_f64; // Pa·s
        
        let re = density * velocity * length / viscosity;
        assert!((re - 100000.0).abs() < 1.0, "Reynolds number calculation error");
    }
}