//! Core abstractions and plugin system for the CFD simulation suite.
//!
//! This crate provides the fundamental traits, types, and plugin infrastructure
//! that all other crates in the suite build upon.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod aggregates;
pub mod boundary;
pub mod domain;
pub mod error;
pub mod factory;
pub mod fluid;
pub mod orchestration;
pub mod plugin;
pub mod problem;
pub mod services;
pub mod solver;
pub mod state;
pub mod time;
pub mod values;

pub use aggregates::{SimulationAggregate, MeshAggregate, SimulationMetadata, PhysicalParameters};
pub use boundary::BoundaryCondition;
pub use domain::Domain;
pub use error::{Error, Result};
pub use fluid::Fluid;
pub use plugin::{Plugin, PluginRegistry, SimulationPlugin};
pub use problem::Problem;
pub use services::{FluidDynamicsService, MeshQualityService, FlowRegime, QualityLevel};
pub use solver::{Solver, SolverConfiguration, SolverConfig, LinearSolverConfig, NetworkSolverConfig};
pub use factory::{SolverFactory, SolverFactoryRegistry, Builder, ConfigurationBuilder, ResourceManager};
pub use orchestration::{SimulationOrchestrator, ExecutionContext, ExecutionContextBuilder, ExecutionMode, PerformanceMetrics};
pub use state::SimulationState;
pub use time::TimeIntegrator;
pub use values::{ReynoldsNumber, Pressure, Velocity, Temperature, DimensionlessNumber};

// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface