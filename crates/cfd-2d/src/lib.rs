#![allow(dead_code)]
#![cfg_attr(not(test), allow(unused))]
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
// Core modules
pub mod constants;
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
// Re-export main types from domain modules
pub use discretization::{
    CentralDifference, ConvectionScheme, ConvectionSchemeFactory, FirstOrderUpwind, HybridScheme,
    PowerLawScheme, QuickScheme,
};
pub use physics::{
    EnergyEquationSolver, KEpsilonModel, MomentumCoefficients, MomentumComponent, MomentumSolver,
    VorticityStreamSolver, WallFunction,
pub use solvers::{
    AdvectionDiffusionSolver, FdmConfig, FluxScheme, FvmConfig, FvmSolver, LbmConfig, LbmSolver,
    PoissonSolver, D2Q9,
// Re-export core types
pub use fields::{Field2D, SimulationFields};
pub use grid::{BoundaryType, Grid2D, GridEdge, StructuredGrid2D};
pub use pressure_velocity::{PressureVelocityConfig, PressureVelocitySolver};
pub use problem::{IncompressibleFlowProblem, IncompressibleFlowSolution};
pub use schemes::{FluxLimiter, Grid2D as SchemeGrid2D, SpatialScheme, TimeScheme};
// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface
