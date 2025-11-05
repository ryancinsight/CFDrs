//! Smagorinsky Large Eddy Simulation (LES) model
//!
//! The Smagorinsky model is one of the most widely used SGS (Sub-Grid Scale)
//! models for LES. It assumes that the SGS stress tensor is proportional to
//! the magnitude of the resolved strain rate tensor.
//!
//! ## Mathematical Formulation
//!
//! The SGS stress tensor is modeled as:
//!
//! τ_ij = -2ν_SGS S_ij
//!
//! where ν_SGS = (C_S Δ)² |S| is the SGS viscosity,
//! C_S is the Smagorinsky constant (typically 0.1-0.2),
//! Δ is the filter width, and |S| is the magnitude of the strain rate.
//!
//! ## Dynamic Procedure
//!
//! The dynamic Smagorinsky model computes C_S locally using:
//!
//! C_S² = <L_ij M_ij> / <M_ij M_ij>
//!
//! where L_ij and M_ij are computed from resolved scales.
//!
//! ## References
//!
//! - Smagorinsky, J. (1963). General circulation experiments with the primitive equations.
//! - Germano, M., et al. (1991). A dynamic subgrid-scale eddy viscosity model.
//! - Lilly, D. K. (1992). A proposed modification of the Germano subgrid-scale closure method.

pub mod config;
pub mod strain;
pub mod viscosity;
pub mod dynamic;
pub mod gpu;
pub mod model;

// Re-export main types
pub use config::SmagorinskyConfig;
pub use model::SmagorinskyLES;

