//! Momentum equation solver with clean domain separation

mod boundary;
mod coefficient_corrections;
mod coefficients;
mod discretization;
mod interpolation;
pub mod muscl;
mod setup;
mod solver;
pub mod tvd_limiters;

#[cfg(test)]
mod tvd_limiter_edge_cases;

#[cfg(test)]
mod algorithm_validation_tests;

pub use boundary::apply_momentum_boundaries;
pub use coefficients::{ConvectionScheme, MomentumCoefficients};
pub use discretization::{CentralDifference, DiscretizationScheme, MusclDiscretization, Upwind};
pub use interpolation::rhie_chow_interpolation;
pub use muscl::{schemes, MusclOrder, MusclReconstruction, MusclScheme};
pub use setup::BoundarySetup;
pub use solver::{MomentumComponent, MomentumSolver};
pub use tvd_limiters::{Minmod, MonotonizedCentral, Superbee, TvdLimiter, VanLeer};
