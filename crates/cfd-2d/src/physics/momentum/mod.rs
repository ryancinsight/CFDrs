//! Momentum equation solver with clean domain separation

mod boundary;
mod coefficients;
mod discretization;
mod interpolation;
mod solver;

pub use boundary::apply_momentum_boundaries;
pub use coefficients::MomentumCoefficients;
pub use discretization::{CentralDifference, DiscretizationScheme, Upwind};
pub use interpolation::rhie_chow_interpolation;
pub use solver::{MomentumComponent, MomentumSolver};
