//! Conservation checking for CFD simulations.
//!
//! This module provides tools to verify that CFD simulations satisfy fundamental
//! conservation laws such as mass, momentum, and energy conservation.

mod energy;
mod history;
mod mass;
mod momentum;
mod report;
mod tolerance;
mod traits;

pub use energy::EnergyConservationChecker;
pub use history::ConservationHistory;
pub use mass::MassConservationChecker;
pub use momentum::MomentumConservationChecker;
pub use report::ConservationReport;
pub use tolerance::ConservationTolerance;
pub use traits::ConservationChecker;
