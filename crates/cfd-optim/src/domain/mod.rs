pub mod candidate;
mod goal;

#[cfg(test)]
pub(crate) use candidate::fixtures;
pub use candidate::{BlueprintCandidate, OperatingPoint, PatientContext};
pub use goal::OptimizationGoal;
