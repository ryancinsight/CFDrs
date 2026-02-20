//! Boundary condition types and blood rheology model selector.
//!
//! All types are re-exported from their canonical homes in `cfd-core` and `cfd-math`.
//! This module is a thin pass-through; no types are defined here.

// Canonical boundary condition types from cfd-core
pub use cfd_core::physics::boundary::{BoundaryCondition, FundamentalBCType as BCType};

// Canonical blood rheology dispatch from cfd-core
pub use cfd_core::physics::fluid::BloodModel;
