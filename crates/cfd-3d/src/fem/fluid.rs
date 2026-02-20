//! Fluid properties for FEM calculations.
//!
//! Re-exported from `cfd-core` — the single source of truth for fluid material models.
//!
//! The `stress_tensor` and `strain_rate_tensor` computations that require
//! `nalgebra::Matrix3` live in [`crate::fem::stress`] as FEM-specific extensions.

pub use cfd_core::physics::fluid::ConstantPropertyFluid as FluidProperties;
