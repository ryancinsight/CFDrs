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
