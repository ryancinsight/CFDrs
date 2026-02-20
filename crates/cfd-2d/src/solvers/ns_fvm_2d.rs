//! Backward-compatibility shim for `ns_fvm_2d`.
//!
//! All types have been moved to [`super::ns_fvm`].
//! This file exists solely to avoid breaking downstream `use crate::solvers::ns_fvm_2d::*` imports.
//! **Prefer importing directly from `ns_fvm` for new code.**

pub use super::ns_fvm::{
    BCType, BloodModel, BoundaryCondition, FlowField2D, NavierStokesSolver2D,
    SIMPLEConfig, SolveResult, StaggeredGrid2D,
};
