//! Validation framework for CFD simulations
//!
//! This module provides tools for validating CFD implementations against:
//! - Analytical solutions
//! - Benchmark problems
//! - Literature results
//! - Method of manufactured solutions

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// CFD validation allows
#![allow(clippy::similar_names)]           // Test variables often similar (u1,u2; err1,err2)
#![allow(clippy::cast_precision_loss)]     // Acceptable in validation calculations  
#![allow(clippy::must_use_candidate)]      // Validation utilities often used in expressions

pub mod analytical;
pub mod analytical_benchmarks;
pub mod benchmarks;
pub mod conservation;
pub mod convergence;
pub mod error_metrics;
pub mod literature;
pub mod manufactured;
pub mod numerical;
pub mod solutions;
pub mod time_integration;

// The public modules are the primary API.
// Users should access types through the module hierarchy:
//   use cfd_validation::analytical::AnalyticalSolution;
//   use cfd_validation::convergence::GridConvergenceIndex;
//   use cfd_validation::conservation::ConservationChecker;
//
// This hierarchical structure provides clear organization:
// - analytical: Known exact solutions (Poiseuille, Couette, etc.)
// - benchmarks: Standard test cases (lid-driven cavity, etc.)
// - conservation: Mass/momentum/energy conservation checks
// - convergence: Grid/temporal convergence studies
// - error_metrics: L2 norm, L∞ norm, etc.
// - literature: Comparison with published results
// - manufactured: Method of manufactured solutions
// - numerical: Numerical analysis tools
// - solutions: Solution comparison utilities
// - time_integration: Temporal accuracy analysis
