//! Validation and benchmarking tools for CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod analytical;
pub mod benchmarks;
pub mod error_metrics;
pub mod convergence;
pub mod conservation;

// TODO: Implement these exports
// pub use analytical::{
//     AnalyticalSolution, PoiseuilleFlow, CouetteFlow, StokesFlow, TaylorGreenVortex,
// };
// pub use benchmarks::{
//     Benchmark, LidDrivenCavity, FlowOverCylinder, BackwardFacingStep,
// };
// pub use error_metrics::{ErrorMetric, L2Norm, LInfNorm, RelativeError};
// pub use convergence::{ConvergenceAnalysis, RichardsonExtrapolation};
// pub use conservation::{ConservationChecker, MassConservation, EnergyConservation};

/// Common validation types and traits
pub mod prelude {
    // TODO: Add exports when implemented
    // pub use crate::{
    //     analytical::{AnalyticalSolution, PoiseuilleFlow},
    //     benchmarks::{Benchmark, LidDrivenCavity},
    //     error_metrics::{ErrorMetric, L2Norm},
    //     convergence::ConvergenceAnalysis,
    //     conservation::ConservationChecker,
    // };
}