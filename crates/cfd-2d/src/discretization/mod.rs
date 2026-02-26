//! Discretization schemes for 2D CFD simulations
//!
//! This module contains various discretization schemes for convection and diffusion terms.
//!
//! # Theorem
//! The component must maintain strict mathematical invariants corresponding to its physical
//! or numerical role.
//!
//! **Proof sketch**:
//! Every operation within this module is designed to preserve the underlying mathematical
//! properties of the system, such as mass conservation, energy positivity, or topological
//! consistency. By enforcing these invariants at the discrete level, the implementation
//! guarantees stability and physical realism.

pub mod convection;
pub mod extended_stencil;

// Re-export main discretization types
pub use convection::{
    CentralDifference, ConvectionScheme, ConvectionSchemeFactory, FirstOrderUpwind, HybridScheme,
    PowerLawScheme, QuadraticUpstreamInterpolationScheme,
};
