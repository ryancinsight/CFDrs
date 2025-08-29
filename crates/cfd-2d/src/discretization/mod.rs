//! Discretization schemes for 2D CFD simulations
//!
//! This module contains various discretization schemes for convection and diffusion terms.

pub mod convection;

// Re-export main discretization types
pub use convection::{
    CentralDifference, ConvectionScheme, ConvectionSchemeFactory, FirstOrderUpwind, HybridScheme,
    PowerLawScheme, QuadraticUpstreamInterpolationScheme,
};
