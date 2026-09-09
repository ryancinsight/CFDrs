//! Time integration validation module.

pub mod integrators;
pub mod results;
pub mod stability_analysis;
pub mod validation;

pub use integrators::{
    ForwardEuler, RungeKutta2, RungeKutta4, TimeIntegratorEnum, TimeIntegratorTrait,
};
pub use results::TimeIntegrationResult;
pub use stability_analysis::{StabilityAnalysisReport, StabilityAnalysisRunner};
pub use validation::TimeIntegrationValidator;

// Edge case tests
#[cfg(test)]
mod edge_case_tests;
