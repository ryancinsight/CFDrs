//! Lagrangian cell-tracking solver for 2D velocity fields.

/// Documentation of population boundary configurations.
pub mod population;
/// Contains velocity fields and analytical boundaries.
pub mod physics;
/// Trajectory solving engine.
pub mod tracker;

pub use population::{CellPopulation, CellRoutingSummary, CellTrajectory, OutletZone, TrackedCell};
pub use physics::{CellTrackerConfig, PoiseuilleFlow2D, PsmBifurcationParams, VelocityFieldInterpolator, AsymmetricBifurcationFlow};
pub use tracker::CellTracker;
