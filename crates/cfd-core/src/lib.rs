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
pub mod domains;
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
pub use solver::{Solver, SolverConfiguration, Configurable, Validatable, IterativeSolver, DirectSolver, SolverConfig, LinearSolverConfig, NetworkSolverConfig};
pub use factory::{AbstractSolverFactory, ConcreteSolverFactory, FactoryCapability, SolverFactoryRegistry, Builder, ConfigurationBuilder, ResourceManager};
pub use orchestration::{SimulationOrchestrator, ExecutionContext, ExecutionContextBuilder, ExecutionMode, PerformanceMetrics};
pub use state::SimulationState;
pub use time::TimeIntegrator;
pub use values::{ReynoldsNumber, Pressure, Velocity, Temperature, DimensionlessNumber};
pub use domains::{
    FlowField, VelocityField, PressureField, TurbulenceModel,
    DiscretizationScheme, TimeIntegrationScheme, LinearSystemSolver,
    MeshGeneration, MeshRefinement, MeshQuality,
    BoundaryConditionType, BoundaryConditionApplicator,
    FluidProperties, SolidProperties, InterfaceProperties
};

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