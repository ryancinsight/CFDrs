pub mod generator;
pub mod materializer;
pub mod parameters;

pub use generator::generate_milestone12_candidate_params;
pub use parameters::CandidateParams;

use crate::domain::BlueprintCandidate;

/// Return lightweight parameter tuples (~100 bytes each) without
/// materializing full `BlueprintCandidate` objects (~10-15 KB each).
/// Callers should materialize on-demand inside evaluation loops to avoid
/// holding 500K × 15 KB ≈ 7.5 GB in memory at once.
#[must_use]
pub fn build_milestone12_candidate_params_only() -> Vec<CandidateParams> {
    generate_milestone12_candidate_params()
}

/// Build the canonical Milestone 12 candidate space by generating parameter
/// combinations and materializing each into a full `BlueprintCandidate`.
///
/// **Warning:** This eagerly materializes ALL candidates (~500K × 15 KB ≈
/// 7.5 GB). Prefer [`build_milestone12_candidate_params_only`] with
/// on-demand materialization for memory-constrained pipelines.
#[must_use]
pub fn build_milestone12_candidate_space() -> Vec<BlueprintCandidate> {
    generate_milestone12_candidate_params()
        .into_iter()
        .map(|params| params.materialize())
        .collect()
}
