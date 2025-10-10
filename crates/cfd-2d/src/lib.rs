//! 2D CFD simulations with domain-based organization.
//!
//! This crate provides 2D computational fluid dynamics functionality organized by domain:
//! - `solvers`: Numerical methods (FDM, FVM, LBM)
//! - `physics`: Physical models (energy, momentum, turbulence)
//! - `discretization`: Numerical schemes for discretization

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// 2D CFD simulation allows
#![allow(clippy::similar_names)]           // CFD variables (u,v,p; nx,ny; dx,dy; i,j) often similar
#![allow(clippy::cast_precision_loss)]     // Performance-critical numerical loops
#![allow(clippy::cast_possible_truncation)] // Grid indices and array sizes typically small
#![allow(clippy::unused_self)]             // Solver trait methods maintain consistent interfaces
#![allow(clippy::must_use_candidate)]      // Solver utilities and getters used in computational contexts
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
#![allow(clippy::too_many_lines)]           // Complex solver implementations require detailed methods
#![allow(clippy::needless_range_loop)]      // Explicit indexing clearer for multi-dimensional CFD arrays
#![allow(clippy::struct_field_names)]       // Field names like field_* common in computational contexts
#![allow(clippy::used_underscore_binding)]  // Underscore prefixed bindings used for intentional partial use

// Core modules
pub mod constants;
pub mod error;
pub mod fields;
pub mod grid;
pub mod problem;

// Domain-organized modules
pub mod discretization;
pub mod physics;
pub mod solvers;

// Algorithm modules
pub mod piso_algorithm;
pub mod pressure_velocity;
pub mod schemes;
pub mod stability;

// The crate's public API is its module hierarchy.
// Users should access types with clear, logical paths:
//   use cfd_2d::solvers::fvm::FvmSolver;
//   use cfd_2d::physics::turbulence::KEpsilonModel;
//   use cfd_2d::discretization::ConvectionScheme;
//   use cfd_2d::fields::SimulationFields;
//   use cfd_2d::grid::StructuredGrid2D;
//
// This hierarchical structure is self-documenting and aligns with Rust best practices.

// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface
