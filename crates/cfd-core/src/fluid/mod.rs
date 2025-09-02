//! Fluid properties and models with trait-based extensibility
//!
//! This module provides a comprehensive framework for fluid modeling in CFD simulations,
//! supporting Newtonian, non-Newtonian, and temperature-dependent fluid behaviors.

use crate::error::Error;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

pub mod database;
pub mod newtonian;
pub mod non_newtonian;
pub mod properties;
pub mod temperature;
pub mod traits;
pub mod validation;

// Re-export core types
pub use database::{air_20c, water_20c};
pub use newtonian::ConstantPropertyFluid;
pub use properties::FluidProperties;
pub use traits::{
    CompressibleFluid, ConstantFluid, Fluid as FluidTrait, FluidState, NonNewtonianFluid,
};

/// Type alias for backward compatibility
pub type Fluid<T> = ConstantPropertyFluid<T>;

// Legacy FluidModel trait removed - use traits::Fluid instead
