//! Branching-junction models and tools.
//!
//! This module organizes branching concerns by responsibility:
//! - `physics`: two-way/three-way branch conservation and constitutive relations.
//! - `solver`: branching-network solver orchestration.
//! - `validation`: analytical and conservation validation utilities.
//!
//! This is the canonical SSOT location for branching junction logic.

pub mod physics;

pub mod solver;

pub mod validation;

pub use physics::{
    ThreeWayBranchJunction, ThreeWayBranchSolution, TwoWayBranchJunction, TwoWayBranchSolution,
};
pub use solver::{BranchingNetworkConfig, BranchingNetworkSolver};
pub use validation::{BranchingValidationResult, BranchingValidator};
