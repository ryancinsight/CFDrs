//! Non-linear solvers for CFD systems.
//!
//! ## Generic algorithms (SSOT: `leto-ops`)
//!
//! Re-exported from `leto_ops::nonlinear` so CFDrs consumers get the same
//! canonical implementation as kwavers and helios:
//!
//! - [`AndersonAccelerator`]: Type-II MGS-QR Anderson Acceleration
//!
//! ## CFD-specific algorithms
//!
//! - [`JfnkSolver`]: Jacobian-Free Newton-Krylov (uses `cfd_core::error::Result`)

// Generic algorithms — SSOT in leto-ops
pub use leto_ops::{AndersonAccelerator, AndersonConfig, AndersonMethod};

// CFD-specific: JFNK uses cfd_core::error::Result
pub mod jfnk;
pub(crate) mod linalg;

pub use jfnk::{JfnkConfig, JfnkConvergence, JfnkSolver};