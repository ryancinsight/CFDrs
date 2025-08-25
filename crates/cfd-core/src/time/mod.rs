//! Time module: integration schemes and time step controllers
//!
//! This module groups time integration algorithms and controllers under a
//! domain-based hierarchy to uphold SSOT/SPOT and SOC.

pub mod integrators;
pub mod controllers;

// Re-export primary interfaces for backwards-compatible API surface
pub use controllers::{AdaptiveTimeStepController, VariableTimeStep};
pub use integrators::{
    BackwardEuler, CrankNicolson, ForwardEuler, RungeKutta2, RungeKutta4, TimeIntegrator,
};

