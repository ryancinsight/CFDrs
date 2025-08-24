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
pub mod pressure_velocity;
pub mod piso_algorithm;
pub mod schemes;

// Re-export main types from domain modules
pub use discretization::{
    ConvectionScheme, ConvectionSchemeFactory,
    FirstOrderUpwind, CentralDifference,
    HybridScheme, PowerLawScheme, QuickScheme
};

pub use physics::{
    EnergyEquationSolver,
    MomentumSolver, MomentumComponent, MomentumCoefficients,
    KEpsilonModel, WallFunction,
    VorticityStreamSolver
};

pub use solvers::{
    PoissonSolver, AdvectionDiffusionSolver, FdmConfig,
    FvmSolver, FvmConfig, FluxScheme,
    LbmSolver, LbmConfig, D2Q9
};

// Re-export core types
pub use fields::{Field2D, SimulationFields};
pub use grid::{Grid2D, StructuredGrid2D, BoundaryType, GridEdge};
pub use pressure_velocity::{PressureVelocitySolver, PressureVelocityConfig};
pub use problem::{IncompressibleFlowProblem, IncompressibleFlowSolution};
pub use schemes::{SpatialScheme, FluxLimiter, TimeScheme, Grid2D as SchemeGrid2D};

// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface