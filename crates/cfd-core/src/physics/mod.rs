//! Physics module consolidating fluid properties, boundary conditions, and constants.
//!
//! This module provides a unified interface to the physical models and properties
//! used in CFD simulations, following a deep vertical module structure.

pub mod boundary;
pub mod cavitation;
pub mod constants;
pub mod fluid;
pub mod fluid_dynamics;
pub mod material;
pub mod values;

// --- Curated Top-Level API ---

/// Re-export of core boundary condition types
pub mod boundary_api {
    pub use super::boundary::{
        BoundaryCondition, BoundaryConditionManager, BoundaryConditionSet, WallType,
    };
}

/// Re-export of core fluid and material types
pub mod material_api {
    pub use super::fluid::{ConstantPropertyFluid, Fluid, FluidProperties};
    pub use super::material::{MaterialDatabase, SolidProperties};
}

/// Re-export of physical and mathematical constants
pub mod constants_api {
    pub use super::constants::mathematical;
    pub use super::constants::physics;
}

/// Re-export of physical value objects
pub mod values_api {
    pub use super::values::{Pressure, ReynoldsNumber, Temperature, Velocity};
}

/// Re-export of fluid dynamics concepts and operations
pub mod dynamics_api {
    pub use super::fluid_dynamics::{
        FlowClassifier, FlowField, FlowOperations, FlowRegime, FluidDynamicsService,
        KEpsilonConstants, KEpsilonModel, MixingLengthModel, PressureField, RANSModel,
        RhieChowInterpolation, ScalarField, SmagorinskyModel, TurbulenceModel, VelocityField,
    };
}

// Re-export most commonly used types directly for convenience
pub use boundary::{BoundaryCondition, BoundaryConditionManager, BoundaryConditionSet};
pub use fluid::{ConstantPropertyFluid, Fluid, FluidProperties};
pub use fluid_dynamics::FluidDynamicsService;
pub use values::{Pressure, ReynoldsNumber, Temperature, Velocity};
