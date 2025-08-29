//! Boundary conditions for CFD simulations
//!
//! This module provides a comprehensive set of boundary condition types
//! following established CFD literature conventions.
//!
//! # References
//! - Versteeg & Malalasekera (2007). An Introduction to Computational Fluid Dynamics
//! - Patankar (1980). Numerical Heat Transfer and Fluid Flow

mod fundamental;
mod inlet_outlet;
mod set;
mod wall;

pub use fundamental::{BoundaryCondition, FundamentalBCType};
pub use inlet_outlet::{InletCondition, OutletCondition};
pub use set::BoundaryConditionSet;
pub use wall::WallType;
