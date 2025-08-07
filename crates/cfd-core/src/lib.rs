//! Core abstractions and plugin system for the CFD simulation suite.
//!
//! This crate provides the fundamental traits, types, and plugin infrastructure
//! that all other crates in the suite build upon.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod boundary;
pub mod domain;
pub mod error;
pub mod fluid;
pub mod plugin;
pub mod problem;
pub mod solver;
pub mod state;
pub mod time;

pub use boundary::BoundaryCondition;
pub use domain::Domain;
pub use error::{Error, Result};
pub use fluid::Fluid;
pub use plugin::{Plugin, PluginRegistry, SimulationPlugin};
pub use problem::Problem;
pub use solver::Solver;
pub use state::SimulationState;
pub use time::TimeIntegrator;

/// Core-specific prelude for internal use and advanced scenarios
///
/// This prelude is primarily for internal use within the CFD suite.
/// Most users should use `cfd_suite::prelude::*` instead for a unified interface.
///
/// Use this prelude only when you need core-specific functionality not exposed
/// in the main prelude, or when developing plugins/extensions.
pub mod prelude {
    pub use crate::{
        boundary::BoundaryCondition,
        domain::{Domain, Domain1D, Domain2D, Domain3D},
        error::{Error, Result},
        fluid::Fluid,
        plugin::{Plugin, PluginRegistry, SimulationPlugin, SolverFactory, PluginMetadata},
        problem::{Problem, ProblemBuilder, ProblemParameters},
        solver::{Solver, SolverConfig, IterativeSolver, DirectSolver},
        state::SimulationState,
        time::{TimeIntegrator, ForwardEuler, BackwardEuler, CrankNicolson},
    };
}