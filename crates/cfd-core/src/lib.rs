#![allow(dead_code)]
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

// #![warn(missing_docs)] // TODO: Enable after documentation sprint
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod aggregates;
pub mod boundary;
pub mod cavitation;
pub mod constants;
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
    pub use crate::solver::{Solver, SolverConfiguration, SolverConfig};
    pub use crate::state::SimulationState;
    pub use crate::time::TimeIntegrator;
    pub use crate::values::{Pressure, Velocity, Temperature, ReynoldsNumber};
    
    // Plugin system - only expose the main trait
    pub use crate::plugin::{Plugin, SimulationPlugin};
}

// Extended API - for plugin developers and advanced users
// These are intentionally not in the prelude to avoid cluttering the namespace

/// Factory system for dynamic solver creation (advanced usage)
pub mod factories {
    pub use crate::factory::{
        ConcreteSolverFactory, FactoryCapability, SolverFactoryRegistry,
    };
}

/// Plugin system internals (for plugin developers)
pub mod plugins {
    pub use crate::plugin::{
        PluginRegistry, PluginHealthStatus,
        PluginMetrics, SystemStatus, SystemHealthSummary
    };
}

/// Extended solver traits (for solver implementors)
pub mod solvers {
    pub use crate::solver::{
        Configurable, Validatable, IterativeSolver, DirectSolver,
        LinearSolverConfig, NetworkSolverConfig
    };
}

// Re-export advanced types in organized namespaces
// These use different names to avoid conflicts with the actual modules

/// Domain-specific abstractions
pub use domains::{
    FlowField, VelocityField, PressureField, TurbulenceModel,
    DiscretizationScheme, TimeIntegrationScheme, LinearSystemSolver,
    MeshGeneration, MeshRefinement, MeshQuality,
    BoundaryConditionApplicator,
    FluidProperties, SolidProperties, InterfaceProperties
};

/// Aggregate types for complex simulations  
pub use aggregates::{
    SimulationAggregate, MeshAggregate, 
    SimulationMetadata, PhysicalParameters
};

/// Service layer abstractions
pub use services::{
    FluidDynamicsService, MeshQualityService, 
    FlowRegime, QualityLevel
};

// Note: TypeErasedFactory and TypeErasedSolver are internal implementation details
// and should NOT be exposed in the public API. They remain accessible only through
// the factory module for those who need to implement custom factories.
// Re-export commonly used types
pub use error::{Error, Result};
pub use fluid::Fluid;
pub use boundary::BoundaryCondition;
pub use problem::Problem;
pub use solver::{Solver, SolverConfiguration, Configurable, Validatable, NetworkSolverConfig, SolverConfig};
