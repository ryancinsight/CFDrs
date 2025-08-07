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

// Export benchmark and conservation modules
pub use benchmarks::{
    Benchmark, BenchmarkResult, BenchmarkSuite, LidDrivenCavity, FlowOverCylinder, BackwardFacingStep,
};
pub use conservation::{ConservationChecker, MassConservation, EnergyConservation, ConservationReport};

/// Validation domain-specific prelude for advanced analysis
///
/// This prelude exports validation-specific functionality not available in the main prelude.
/// Use this when performing detailed validation studies, implementing custom benchmarks,
/// or conducting comprehensive error analysis.
///
/// For basic validation functionality, prefer `cfd_suite::prelude::*`.
pub mod prelude {
    // === Advanced Error Analysis ===
    // Detailed error metrics and statistical analysis
    pub use crate::{
        error_metrics::{
            RelativeError, RootMeanSquareError, MeanAbsoluteError, NormalizedRMSE,
            NormalizationMethod, ErrorStatistics, ErrorAnalysis
        },
        convergence::{
            ConvergenceStatus, GridConvergenceIndex, RichardsonExtrapolation,
            ConvergenceOrder
        },
    };

    // === Benchmark Framework ===
    // Advanced benchmarking and validation tools
    pub use crate::{
        benchmarks::{BenchmarkSuite, BenchmarkRunner, BenchmarkConfig, ValidationReport},
        conservation::{ConservationReport, ConservationTolerance, ConservationHistory},
    };
}