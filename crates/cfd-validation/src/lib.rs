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

// Export analytical solutions
pub use analytical::{
    AnalyticalSolution, PoiseuilleFlow, CouetteFlow, StokesFlow, TaylorGreenVortex,
    AnalyticalUtils,
};

// Export error metrics
pub use error_metrics::{
    ErrorMetric, L2Norm, L1Norm, LInfNorm, RelativeError, RootMeanSquareError,
    MeanAbsoluteError, NormalizedRMSE, NormalizationMethod, ErrorStatistics, ErrorAnalysis,
};

// Export convergence analysis
pub use convergence::{
    ConvergenceStudy, ConvergenceOrder, RichardsonExtrapolation, GridConvergenceIndex,
    ConvergenceAnalysis, ConvergenceStatus, ConvergenceCriterion,
};

// TODO: Implement these exports when modules are completed
// pub use benchmarks::{
//     Benchmark, LidDrivenCavity, FlowOverCylinder, BackwardFacingStep,
// };
// pub use conservation::{ConservationChecker, MassConservation, EnergyConservation};

/// Common validation types and traits
pub mod prelude {
    pub use crate::{
        analytical::{AnalyticalSolution, PoiseuilleFlow, CouetteFlow, TaylorGreenVortex, AnalyticalUtils},
        error_metrics::{ErrorMetric, L2Norm, L1Norm, LInfNorm, ErrorStatistics, ErrorAnalysis},
        convergence::{ConvergenceAnalysis, ConvergenceStudy, RichardsonExtrapolation},
        // TODO: Add when implemented
        // benchmarks::{Benchmark, LidDrivenCavity},
        // conservation::ConservationChecker,
    };
}