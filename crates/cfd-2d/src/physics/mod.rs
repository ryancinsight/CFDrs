//! Physics models for 2D CFD simulations
//!
//! This module contains physical models including energy, momentum, turbulence, and vorticity.

pub mod energy;
pub mod immersed_boundary;
pub mod momentum;
pub mod turbulence;
pub mod vorticity_stream;

// Re-export main physics types
pub use energy::EnergyEquationSolver;
pub use immersed_boundary::{ImmersedBoundaryMethod, ImmersedBoundaryConfig, BoundaryPoint};
pub use momentum::{MomentumCoefficients, MomentumComponent, MomentumSolver};
pub use turbulence::{KEpsilonModel, WallFunction};
pub use vorticity_stream::VorticityStreamSolver;
