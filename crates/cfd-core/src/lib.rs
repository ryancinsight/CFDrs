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

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::{
        boundary::BoundaryCondition,
        domain::Domain,
        error::{Error, Result},
        fluid::Fluid,
        plugin::{Plugin, PluginRegistry, SimulationPlugin},
        problem::Problem,
        solver::Solver,
        state::SimulationState,
        time::TimeIntegrator,
    };
}