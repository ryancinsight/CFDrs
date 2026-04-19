//! Lagrangian cell-tracking solver for 2D velocity fields.

/// Contains velocity fields and analytical boundaries.
pub mod physics;
/// Documentation of population boundary configurations.
pub mod population;
/// Trajectory solving engine.
pub mod tracker;

pub use physics::{
    AsymmetricBifurcationFlow, CellTrackerConfig, PoiseuilleFlow2D, PsmBifurcationParams,
    VelocityFieldInterpolator,
};
pub use population::{CellPopulation, CellRoutingSummary, CellTrajectory, OutletZone, TrackedCell};
pub use tracker::CellTracker;
