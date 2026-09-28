//! Milestone 12 blueprint-native candidate-space generation.
//!
//! ## Module hierarchy
//!
//! | Module | Responsibility |
//! |--------|----------------|
//! | [`builder`] | Primitive selective blueprint candidate construction |
//! | [`dimensions`] | Milestone 12 selective sweep dimensions |
//! | `sampling` | Tyche-backed stratified Milestone 12 candidate generation |
//! | [`sweep`] | Milestone 12 selective sweep |

mod builder;
mod dimensions;
mod sampling;
mod sweep;

use crate::domain::BlueprintCandidate;
pub use sampling::{generate_milestone12_lhs_params, generate_milestone12_lhs_params_seeded};
pub use sweep::milestone12::CandidateParams;

/// Return lightweight parameter tuples (~100 bytes each) for lazy
/// materialization. This is the memory-efficient alternative to
/// [`build_milestone12_blueprint_candidate_space`].
#[must_use]
pub fn build_milestone12_candidate_params() -> Vec<CandidateParams> {
    sweep::milestone12::build_milestone12_candidate_params_only()
}

/// Build only the canonical primitive selective candidates needed by the
/// Milestone 12 report pipeline.
///
/// This avoids allocating and scanning unrelated topology families when the
/// report only consumes the selective-routing lineage:
/// `Option 1 base -> Option 2 venturi -> GA refinement`.
#[must_use]
pub(crate) fn build_milestone12_candidate_space() -> Vec<BlueprintCandidate> {
    sweep::milestone12::build_milestone12_candidate_space()
}

/// Build the Milestone 12 selective-routing candidate space directly as
/// blueprint-native candidates.
///
/// **Warning:** eagerly materializes ~500K candidates (~7.5 GB). Prefer
/// [`build_milestone12_candidate_params`] for memory-constrained pipelines.
pub fn build_milestone12_blueprint_candidate_space()
-> Result<Vec<BlueprintCandidate>, crate::error::OptimError> {
    Ok(build_milestone12_candidate_space())
}

#[cfg(test)]
mod tests;
