//! Boundary conditions for CFD simulations
//!
//! This module provides a comprehensive set of boundary condition types
//! following established CFD literature conventions.
//!
//! # References
//! - Versteeg & Malalasekera (2007). An Introduction to Computational Fluid Dynamics
//! - Patankar (1980). Numerical Heat Transfer and Fluid Flow

mod applicator;
mod applicators;
mod error;
mod fundamental;
mod geometry;
mod ghost_cells;
mod inlet_outlet;
mod manager;
mod set;
mod specification;
mod time_dependent;
mod wall;

#[cfg(test)]
mod edge_case_tests;

pub use applicator::BoundaryConditionApplicator;
pub use applicators::{DirichletApplicator, NeumannApplicator, RobinApplicator};
pub use error::BoundaryError;
pub use fundamental::{BoundaryCondition, FundamentalBCType};
pub use geometry::{BoundaryGeometry, BoundaryRegion};
pub use ghost_cells::GhostCellCalculator;
pub use inlet_outlet::{InletCondition, OutletCondition};
pub use manager::BoundaryConditionManager;
pub use set::BoundaryConditionSet;
pub use specification::BoundaryConditionSpec;
pub use time_dependent::{TimeDependentSpec, TimeFunctionType};
pub use wall::WallType;
