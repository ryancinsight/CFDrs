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
#![allow(clippy::missing_errors_doc)]      // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)]      // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)]          // Signed to unsigned casts common in CFD indexing
#![allow(clippy::cast_possible_wrap)]      // Wrap-around acceptable for grid indices
#![allow(clippy::too_many_arguments)]      // CFD functions often need many physical parameters
#![allow(clippy::float_cmp)]               // Float comparisons necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)]        // Result types maintained for API consistency
#![allow(clippy::items_after_statements)]  // Helper functions after statements improve readability
#![allow(clippy::many_single_char_names)]       // Mathematical notation (i,j,k,x,y,z) is standard
#![allow(clippy::unreadable_literal)]      // Long literals used for precise physical constants
#![allow(clippy::redundant_closure_for_method_calls)] // Closures improve readability in numerical pipelines
#![allow(clippy::doc_markdown)]            // Math notation doesn't need backticks
#![allow(clippy::needless_pass_by_value)]  // Pass by value for Copy types is idiomatic
#![allow(clippy::return_self_not_must_use)]  // Builder patterns used internally
#![allow(clippy::ptr_arg)]                 // &Vec used for API compatibility
#![allow(clippy::should_implement_trait)]  // CFD-specific trait implementations

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
