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
#![allow(clippy::too_many_lines)] // GPU kernel implementations require detailed logic
#![allow(clippy::struct_field_names)]
// Field names like field_* common in kernel contexts
// CFD-specific allows for production numerical code
#![allow(clippy::similar_names)] // Mathematical variables often have similar names (u,v,w; p,q,r)
#![allow(clippy::cast_precision_loss)] // Precision loss acceptable in CFD for performance
#![allow(clippy::cast_possible_truncation)] // GPU buffer sizes and indices are typically small
#![allow(clippy::unused_self)] // Trait methods may not use self but maintain interface consistency
#![allow(clippy::must_use_candidate)] // CFD functions often have side effects or are utilities
#![allow(clippy::missing_errors_doc)] // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)] // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)] // Signed to unsigned casts common in CFD indexing
#![allow(clippy::cast_possible_wrap)] // Wrap-around acceptable for grid indices
#![allow(clippy::too_many_arguments)] // CFD functions often need many physical parameters
#![allow(clippy::float_cmp)] // Float comparisons necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)] // Result types maintained for API consistency
#![allow(clippy::items_after_statements)] // Helper functions after statements improve readability
#![allow(clippy::many_single_char_names)] // Mathematical notation (i,j,k,x,y,z) is standard
#![allow(clippy::unreadable_literal)] // Long literals used for precise physical constants
#![allow(clippy::redundant_closure_for_method_calls)] // Closures improve readability in numerical pipelines
#![allow(clippy::doc_markdown)] // Math notation doesn't need backticks
#![allow(clippy::needless_pass_by_value)] // Pass by value for Copy types is idiomatic
#![allow(clippy::return_self_not_must_use)] // Builder patterns used internally
#![allow(clippy::ptr_arg)] // &Vec used for API compatibility
#![allow(clippy::should_implement_trait)] // CFD-specific trait implementations

pub mod abstractions;
pub mod compute;
pub mod domain;
pub mod error;
pub mod interpolation;
pub use crate::management::conversion;
/// Domain management, coordination, and plugin system
pub mod management;
pub mod physics;

/// Prelude module for convenient imports of commonly used types
///
/// # Example
/// ```rust
/// use cfd_core::prelude::*;
/// ```
/// Common imports for CFD core functionality
pub mod prelude {

    // Essential types that users will directly interact with
    pub use crate::error::{Error, Result};
    pub use crate::domain::Domain;
    pub use crate::physics::boundary::{BoundaryCondition, WallType};
    pub use crate::physics::fluid::{ConstantPropertyFluid, Fluid};
    pub use crate::abstractions::problem::Problem;
    pub use crate::compute::solver::{Solver, SolverConfig, SolverConfiguration};
    pub use crate::abstractions::state::SimulationState;
    pub use crate::compute::time::TimeIntegrator;
    pub use crate::physics::values::{Pressure, ReynoldsNumber, Temperature, Velocity};

    // Plugin system - only expose the main trait
    pub use crate::management::plugin::{Plugin, SimulationPlugin};
}

// Extended API - for plugin developers and advanced users
// These are intentionally not in the prelude to avoid cluttering the namespace

/// Factory system for dynamic solver creation (advanced usage)
pub mod factories {
    pub use crate::management::factory::{ConcreteSolverFactory, FactoryCapability, SolverFactoryRegistry};
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
    pub use crate::physics::fluid_dynamics::flow_regimes::FlowRegime;
    pub use crate::physics::fluid_dynamics::service::FluidDynamicsService;
    pub use crate::domain::mesh::{MeshQualityService, QualityLevel};
}

// Note: TypeErasedFactory and TypeErasedSolver are internal implementation details
// and should NOT be exposed in the public API. They remain accessible only through
// the factory module for those who need to implement custom factories.

// IMPORTANT: Core types should be imported via the prelude module or their full paths.
// Direct re-exports at the crate root have been removed to enforce a single, clear import path.
