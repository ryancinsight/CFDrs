//! Physics module consolidating fluid properties, boundary conditions, and constants.
//!
//! This module provides a unified interface to the physical models and properties
//! used in CFD simulations, following a deep vertical module structure.

pub mod boundary;
pub mod cavitation;
pub mod constants;
pub mod fluid;
pub mod fluid_dynamics;
pub mod hemolysis;
pub mod material;
pub mod microfluidics;
pub mod values;

#[allow(missing_docs)]
pub mod api {
    pub use super::boundary::{
        BoundaryCondition, BoundaryConditionManager, BoundaryConditionSet, WallType,
    };
    pub use super::constants::{mathematical, physics};
    pub use super::fluid::{ConstantPropertyFluid, Fluid, FluidProperties};
    pub use super::fluid_dynamics::{
        FlowClassifier, FlowField, FlowOperations, FlowRegime, FluidDynamicsService, PressureField,
        RANSModel, RhieChowInterpolation, ScalarField, TurbulenceModel, VelocityField,
    };
    pub use super::hemolysis::{
        BloodTrauma, BloodTraumaSeverity, HemolysisCalculator, HemolysisModel, PlateletActivation,
    };
    pub use super::material::{MaterialDatabase, SolidProperties};
    pub use super::values::{Pressure, ReynoldsNumber, Temperature, Velocity};
}
