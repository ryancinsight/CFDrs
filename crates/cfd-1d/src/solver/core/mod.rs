//! Modularized network solver for 1D CFD analysis
//!
//! This module provides a comprehensive solver for analyzing fluid flow in microfluidic
//! networks using sparse linear algebra and circuit analogies.

mod anderson_acceleration;
mod config;
mod convergence;
mod geometry;
mod linear_system;
mod matrix_assembly;
mod network_solver;
mod problem;
mod solver_detection;
mod state;
mod status;
/// Transient solvers for composition and droplet tracking over time.
pub mod transient;
mod vector_bridge;
/// Workspace state and allocation for the 1D solver loop.
pub mod workspace;

pub use config::SolverConfig;
pub use convergence::ConvergenceChecker;
pub use geometry::NetworkDomain;
pub use linear_system::{LinearSolverMethod, LinearSystemSolver};
pub use matrix_assembly::MatrixAssembler;
pub use network_solver::NetworkSolver;
pub use problem::NetworkProblem;
pub use state::NetworkState;
pub use status::{PrimarySolveDiagnostics, PrimarySolveError, SolveFailureReason, SolvePathStatus};
pub use workspace::SolverWorkspace;

pub use transient::composition::{
    BloodEdgeTransportConfig, CompositionState, EdgeFlowEvent, InletCompositionEvent,
    InletHematocritEvent, MixtureComposition, PressureBoundaryEvent, SimulationTimeConfig,
    TransientCompositionSimulator, BLOOD_PLASMA_FLUID_ID, BLOOD_RBC_FLUID_ID,
};
pub use transient::droplets::{
    ChannelOccupancy, DropletBoundary, DropletInjection, DropletPosition, DropletSnapshot,
    DropletSplitPolicy, DropletState, DropletTrackingState, SplitMode, TransientDropletSimulator,
};
