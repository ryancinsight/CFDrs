//! Core abstractions and plugin system for the CFD simulation suite.
//!
//! This crate provides the fundamental traits, types, and plugin infrastructure
//! that all other crates in the suite build upon.
//!
//! # Intended Usage
//!
//! This crate is a foundational component of the `cfd-suite`. For the best
//! experience and a unified API, it is recommended to use the `cfd-suite`
//! crate directly, which provides a curated prelude module.
//!
//! If using `cfd-core` as a standalone library, you can either:
//! - Import types via their full paths (e.g., `cfd_core::solver::Solver`)
//! - Use the local prelude: `use cfd_core::prelude::*;`

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod aggregates;
pub mod boundary;
pub mod broadcast;
pub mod cavitation;
pub mod compute;
pub mod constants;
pub mod conversion;
pub mod domain;
pub mod domains;
pub mod error;
pub mod factory;
pub mod fluid;
pub mod interpolation;
pub mod plugin;
pub mod problem;
pub mod services;
pub mod solver;
pub mod state;
pub mod time;
pub mod values;

/// Prelude module for convenient imports of commonly used types
///
/// # Example
/// ```rust
/// use cfd_core::prelude::*;
/// ```
pub mod prelude {
    //! Common imports for CFD core functionality

    // Essential types that users will directly interact with
    pub use crate::boundary::{BoundaryCondition, WallType};
    pub use crate::domain::Domain;
    pub use crate::error::{Error, Result};
    pub use crate::fluid::Fluid;
    pub use crate::problem::Problem;
    pub use crate::solver::{Solver, SolverConfig, SolverConfiguration};
    pub use crate::state::SimulationState;
    pub use crate::time::TimeIntegrator;
    pub use crate::values::{Pressure, ReynoldsNumber, Temperature, Velocity};

    // Plugin system - only expose the main trait
    pub use crate::plugin::{Plugin, SimulationPlugin};
}

// Extended API - for plugin developers and advanced users
// These are intentionally not in the prelude to avoid cluttering the namespace

/// Factory system for dynamic solver creation (advanced usage)
pub mod factories {
    pub use crate::factory::{ConcreteSolverFactory, FactoryCapability, SolverFactoryRegistry};
}

/// Plugin system internals (for plugin developers)
pub mod plugins {
    pub use crate::plugin::{
        PluginHealthStatus, PluginMetrics, PluginRegistry, SystemHealthSummary, SystemStatus,
    };
}

/// Extended solver traits (for solver implementors)
pub mod solvers {
    pub use crate::solver::configuration::{LinearSolverConfig, NetworkSolverConfig};
    pub use crate::solver::{Configurable, DirectSolver, IterativeSolver, Validatable};
}

// Extended API - Advanced types in organized namespaces
// These are intentionally not in the prelude to avoid cluttering the namespace

/// Domain-specific abstractions, traits, and types.
///
/// This module provides a curated set of traits that define the
/// core physics and numerical behaviors for CFD simulations.
pub mod domain_traits {
    pub use crate::domains::{
        BoundaryConditionApplicator, DiscretizationScheme, FlowField, FluidProperties,
        InterfaceProperties, LinearSystemSolver, MeshGeneration, MeshQuality, MeshRefinement,
        PressureField, SolidProperties, TimeIntegrationScheme, TurbulenceModel, VelocityField,
    };
}

/// Aggregate types for complex simulations
pub mod aggregates_api {
    pub use crate::aggregates::{
        PhysicalParameters, ProblemAggregate, SimulationAggregate, SimulationMetadata,
    };
}

/// Service layer abstractions
pub mod services_api {
    pub use crate::services::{FlowRegime, FluidDynamicsService, MeshQualityService, QualityLevel};
}

// Note: TypeErasedFactory and TypeErasedSolver are internal implementation details
// and should NOT be exposed in the public API. They remain accessible only through
// the factory module for those who need to implement custom factories.

// IMPORTANT: Core types should be imported via the prelude module or their full paths.
// Direct re-exports at the crate root have been removed to enforce a single, clear import path.
