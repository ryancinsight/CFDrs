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
// Field names like field_* common in kernel contexts
// CFD-specific allows for production numerical code

pub mod abstractions;
pub mod compute;
pub mod error;
pub mod geometry;
pub use crate::management::conversion;
/// Domain management, coordination, and plugin system
pub mod management;
pub mod physics;
pub mod scalar;

pub use crate::scalar::CfdScalar;

/// Prelude module for convenient imports of commonly used types
///
/// # Example
/// ```rust
/// use cfd_core::prelude::*;
/// ```
/// Common imports for CFD core functionality
pub mod prelude {

    // Essential types that users will directly interact with
    pub use crate::abstractions::problem::Problem;
    pub use crate::abstractions::state::SimulationState;
    pub use crate::compute::solver::{Solver, SolverConfig, SolverConfiguration};
    pub use crate::compute::time::TimeIntegrator;
    pub use crate::error::{Error, Result};
    pub use crate::geometry::Domain;
    pub use crate::physics::boundary::{BoundaryCondition, WallType};
    pub use crate::physics::fluid::ConstantPropertyFluid;
    pub use crate::physics::values::{Pressure, ReynoldsNumber, Temperature, Velocity};

    // Plugin system - only expose the main trait
    pub use crate::management::plugin::{Plugin, SimulationPlugin};
}

// Extended API - for plugin developers and advanced users
// These are intentionally not in the prelude to avoid cluttering the namespace

/// Factory system for dynamic solver creation (advanced usage)
pub mod factories {
    pub use crate::management::factory::{
        ConcreteSolverFactory, FactoryCapability, SolverFactoryRegistry,
    };
}

/// Plugin system internals (for plugin developers)
pub mod plugins {
    pub use crate::management::plugin::{
        PluginHealthStatus, PluginMetrics, PluginRegistry, SystemHealthSummary, SystemStatus,
    };
}

/// Extended solver traits (for solver implementors)
pub mod solvers {
    pub use crate::compute::solver::{Configurable, DirectSolver, IterativeSolver, Validatable};
}

/// Aggregate types for complex simulations
pub mod aggregates_api {
    pub use crate::management::aggregates::{
        PhysicalParameters, ProblemAggregate, SimulationAggregate, SimulationMetadata,
    };
}

/// Service layer abstractions
pub mod services_api {
    pub use crate::geometry::mesh::{MeshQualityService, QualityLevel};
    pub use crate::physics::fluid_dynamics::flow_regimes::FlowRegime;
    pub use crate::physics::fluid_dynamics::service::FluidDynamicsService;
}

// Note: TypeErasedFactory and TypeErasedSolver are internal implementation details
// and should NOT be exposed in the public API. They remain accessible only through
// the factory module for those who need to implement custom factories.

// IMPORTANT: Core types should be imported via the prelude module or their full paths.
// Direct re-exports at the crate root have been removed to enforce a single, clear import path.
